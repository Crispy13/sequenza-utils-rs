use crate::cli::Bam2SeqzArgs;
use crate::errors::{AppError, Result};
use crate::external_tools::{CommandStream, ExternalTools};
use crate::seqz_core::{SeqzParams, do_seqz, seqz_header};
use crate::writer;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::process::Child;
use std::sync::Arc;
use std::time::Duration;
use tempfile::NamedTempFile;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PipelinePlan {
    pub use_bam_input: bool,
    pub regions: Vec<String>,
    pub output_paths: Vec<String>,
    pub qlimit_ascii: i32,
}

pub fn build_plan(args: &Bam2SeqzArgs) -> Result<PipelinePlan> {
    let output_paths = if args.nproc > 1 {
        args.chr
            .iter()
            .map(|region| output_for_region(&args.out, region))
            .collect()
    } else {
        vec![args.out.clone()]
    };

    Ok(PipelinePlan {
        use_bam_input: !args.pileup,
        regions: args.chr.clone(),
        output_paths,
        qlimit_ascii: args.qlimit_ascii(),
    })
}

pub fn run(args: &Bam2SeqzArgs) -> Result<()> {
    let plan = build_plan(args)?;
    let tools = ExternalTools::from_args(args);

    if args.nproc > 1 {
        for (index, region) in args.chr.iter().enumerate() {
            let output_path = plan
                .output_paths
                .get(index)
                .ok_or_else(|| AppError::ParseError {
                    message: "missing output path for region".to_string(),
                })?;
            run_one(args, &tools, std::slice::from_ref(region), output_path)?;
        }
    } else {
        run_one(args, &tools, &args.chr, &args.out)?;
    }

    Ok(())
}

fn run_one(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    regions: &[String],
    output: &str,
) -> Result<()> {
    let mut tumor_stream = open_pileup_stream(args, tools, &args.tumor, regions)?;
    let mut alt_stream = if let Some(normal2_path) = &args.normal2 {
        Some(open_pileup_stream(args, tools, normal2_path, regions)?)
    } else {
        None
    };

    let mut gc_filter: HashSet<String> = HashSet::new();
    if !regions.is_empty() {
        gc_filter.extend(
            regions
                .iter()
                .map(|region| region_chromosome(region).to_string()),
        );
    } else if !args.pileup {
        gc_filter.extend(tools.list_bam_chromosomes(&args.tumor)?);
        if let Some(normal2_path) = &args.normal2 {
            gc_filter.extend(tools.list_bam_chromosomes(normal2_path)?);
        }
    }

    let gc_intervals = parse_gc_intervals_from_file(
        &args.gc,
        if gc_filter.is_empty() {
            None
        } else {
            Some(&gc_filter)
        },
    )?;

    let mut normal_stream = open_pileup_stream(args, tools, &args.normal, regions)?;
    let mut tumor_current = tumor_stream.next_record()?;
    let mut alt_current = if let Some(stream) = alt_stream.as_mut() {
        stream.next_record()?
    } else {
        None
    };

    let params = SeqzParams {
        depth_sum: args.depth_sum,
        qlimit: args.qlimit_ascii().clamp(0, i32::from(u8::MAX)) as u8,
        hom_t: args.hom,
        het_t: args.het,
        het_f: args.het_f,
        het_only: false,
    };

    let mut progress = PipelineProgress::new(args.progress, output, regions);

    writer::with_text_output_writer(output, tools, move |out| {
        out.write_all(seqz_header().join("\t").as_bytes())?;
        out.write_all(b"\n")?;

        let mut normal_line = String::with_capacity(96);
        let mut tumor_line = String::with_capacity(96);
        let mut alt_line = String::with_capacity(96);

        while let Some(normal) = normal_stream.next_record()? {
            progress.on_processed(&normal.chromosome, normal.position);
            advance_to_target(&mut tumor_current, &mut tumor_stream, &normal)?;
            let Some(tumor) = tumor_current.as_ref() else {
                break;
            };
            if compare_coordinates(tumor, &normal) != Ordering::Equal {
                continue;
            }

            let Some(gc) = gc_value_at(&gc_intervals, &normal.chromosome, normal.position) else {
                continue;
            };

            normal.write_data_line(&mut normal_line);
            tumor.write_data_line(&mut tumor_line);

            let seqz_fields = if let Some(stream) = alt_stream.as_mut() {
                advance_to_target(&mut alt_current, stream, &normal)?;
                if let Some(alt) = alt_current.as_ref() {
                    if compare_coordinates(alt, &normal) == Ordering::Equal {
                        alt.write_data_line(&mut alt_line);
                        do_seqz(
                            &[
                                normal_line.as_str(),
                                tumor_line.as_str(),
                                gc,
                                alt_line.as_str(),
                            ],
                            &params,
                        )
                    } else {
                        do_seqz(&[normal_line.as_str(), tumor_line.as_str(), gc], &params)
                    }
                } else {
                    do_seqz(&[normal_line.as_str(), tumor_line.as_str(), gc], &params)
                }
            } else {
                do_seqz(&[normal_line.as_str(), tumor_line.as_str(), gc], &params)
            };

            if let Some(seqz) = seqz_fields {
                progress.on_emitted();
                let _ = write!(out, "{}\t{}", normal.chromosome, normal.position);
                for field in seqz {
                    out.write_all(b"\t")?;
                    out.write_all(field.as_bytes())?;
                }
                out.write_all(b"\n")?;
            }
        }
        progress.finish();
        Ok(())
    })?;

    if output.ends_with(".gz") {
        tools.tabix_index_seqz(output)?;
    }
    Ok(())
}

#[derive(Debug)]
struct PipelineProgress {
    progress_bar: Option<ProgressBar>,
    processed_records: u64,
    emitted_records: u64,
    finished: bool,
}

impl PipelineProgress {
    const UPDATE_EVERY: u64 = 10_000;

    fn new(enabled: bool, output: &str, regions: &[String]) -> Self {
        let progress_bar = if enabled {
            let bar = ProgressBar::new_spinner();
            bar.set_draw_target(ProgressDrawTarget::stderr_with_hz(4));
            let style = ProgressStyle::with_template("{spinner:.green} {elapsed_precise} {msg}")
                .unwrap_or_else(|_| ProgressStyle::default_spinner());
            bar.set_style(style);
            bar.enable_steady_tick(Duration::from_millis(200));

            let scope = if regions.is_empty() {
                "all".to_string()
            } else {
                regions.join(",")
            };
            bar.set_message(format!("starting output={output} regions={scope}"));
            Some(bar)
        } else {
            None
        };

        Self {
            progress_bar,
            processed_records: 0,
            emitted_records: 0,
            finished: false,
        }
    }

    fn on_processed(&mut self, chromosome: &str, position: i32) {
        self.processed_records += 1;
        if self.processed_records == 1 || self.processed_records.is_multiple_of(Self::UPDATE_EVERY)
        {
            self.set_message(chromosome, position);
        }
    }

    fn on_emitted(&mut self) {
        self.emitted_records += 1;
    }

    fn finish(&mut self) {
        if let Some(bar) = &self.progress_bar {
            bar.finish_with_message(format!(
                "done processed={} emitted={}",
                self.processed_records, self.emitted_records
            ));
        }
        self.finished = true;
    }

    fn set_message(&self, chromosome: &str, position: i32) {
        if let Some(bar) = &self.progress_bar {
            bar.set_message(format!(
                "processed={} emitted={} locus={chromosome}:{position}",
                self.processed_records, self.emitted_records
            ));
        }
    }
}

impl Drop for PipelineProgress {
    fn drop(&mut self) {
        if !self.finished
            && let Some(bar) = &self.progress_bar
        {
            bar.finish_and_clear();
        }
    }
}

fn open_pileup_stream(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    input_path: &str,
    regions: &[String],
) -> Result<PileupStream> {
    if args.pileup {
        if regions.is_empty() {
            PileupStream::from_path(input_path)
        } else {
            let tempfile = tools.run_tabix_pileup_to_tempfile(input_path, regions)?;
            PileupStream::from_tempfile(tempfile)
        }
    } else {
        let fasta = args
            .fasta
            .as_deref()
            .ok_or_else(|| AppError::MissingRequired {
                field: "--fasta (required when input is BAM)".to_string(),
            })?;
        if regions.len() <= 1 {
            let stream = tools.spawn_samtools_mpileup_stream(
                input_path,
                fasta,
                regions.first().map(String::as_str),
            )?;
            PileupStream::from_command_stream(stream)
        } else {
            let tempfile = tools.run_samtools_mpileup_to_tempfile(input_path, fasta, regions)?;
            PileupStream::from_tempfile(tempfile)
        }
    }
}

struct PileupStream {
    reader: Box<dyn BufRead>,
    line_buffer: String,
    _tempfile: Option<NamedTempFile>,
    process: Option<StreamProcess>,
}

#[derive(Debug)]
struct StreamProcess {
    child: Child,
    stderr_capture: NamedTempFile,
    command: String,
    finished: bool,
}

impl StreamProcess {
    fn read_stderr(&self) -> String {
        let mut bytes = Vec::new();
        if let Ok(mut file) = self.stderr_capture.reopen() {
            let _ = file.read_to_end(&mut bytes);
        }
        String::from_utf8_lossy(&bytes).to_string()
    }
}

impl PileupStream {
    fn from_path(path: &str) -> Result<Self> {
        if path.ends_with(".gz") {
            let file = File::open(path)?;
            let decoder = GzDecoder::new(file);
            Ok(Self {
                reader: Box::new(BufReader::new(decoder)),
                line_buffer: String::new(),
                _tempfile: None,
                process: None,
            })
        } else {
            let file = File::open(path)?;
            Ok(Self {
                reader: Box::new(BufReader::new(file)),
                line_buffer: String::new(),
                _tempfile: None,
                process: None,
            })
        }
    }

    fn from_tempfile(tempfile: NamedTempFile) -> Result<Self> {
        let reader_file = tempfile.reopen()?;
        Ok(Self {
            reader: Box::new(BufReader::new(reader_file)),
            line_buffer: String::new(),
            _tempfile: Some(tempfile),
            process: None,
        })
    }

    fn from_command_stream(stream: CommandStream) -> Result<Self> {
        let process = StreamProcess {
            child: stream.child,
            stderr_capture: stream.stderr_capture,
            command: stream.command,
            finished: false,
        };
        Ok(Self {
            reader: Box::new(BufReader::new(stream.stdout)),
            line_buffer: String::new(),
            _tempfile: None,
            process: Some(process),
        })
    }

    fn ensure_process_completed(&mut self) -> Result<()> {
        let Some(process) = self.process.as_mut() else {
            return Ok(());
        };
        if process.finished {
            return Ok(());
        }

        let status = process.child.wait()?;
        process.finished = true;
        if !status.success() {
            return Err(AppError::CommandFailed {
                command: process.command.clone(),
                code: status.code(),
                stderr: process.read_stderr(),
            });
        }
        Ok(())
    }

    fn next_record(&mut self) -> Result<Option<PileupRecord>> {
        loop {
            self.line_buffer.clear();
            let read = self.reader.read_line(&mut self.line_buffer)?;
            if read == 0 {
                self.ensure_process_completed()?;
                return Ok(None);
            }
            let line = self.line_buffer.trim_end_matches(['\n', '\r']);
            if line.trim().is_empty() {
                continue;
            }
            return parse_pileup_record(line).map(Some);
        }
    }
}

impl Drop for PileupStream {
    fn drop(&mut self) {
        let Some(process) = self.process.as_mut() else {
            return;
        };
        if process.finished {
            return;
        }

        match process.child.try_wait() {
            Ok(Some(_)) => {
                process.finished = true;
            }
            Ok(None) => {
                let _ = process.child.kill();
                let _ = process.child.wait();
                process.finished = true;
            }
            Err(_) => {
                let _ = process.child.kill();
            }
        }
    }
}

#[derive(Debug, Clone)]
struct PileupRecord {
    chromosome: String,
    position: i32,
    reference: String,
    depth: i32,
    pileup: String,
    quality: String,
}

impl PileupRecord {
    fn write_data_line(&self, buffer: &mut String) {
        buffer.clear();
        buffer.push_str(&self.reference);
        buffer.push('\t');
        let _ = write!(buffer, "{}", self.depth);
        buffer.push('\t');
        buffer.push_str(&self.pileup);
        buffer.push('\t');
        buffer.push_str(&self.quality);
    }
}

fn parse_pileup_record(line: &str) -> Result<PileupRecord> {
    let parts = line.splitn(6, '\t').collect::<Vec<_>>();
    if parts.len() < 6 {
        return Err(AppError::ParseError {
            message: format!("invalid pileup line: {line}"),
        });
    }

    let position = parts[1].parse::<i32>().map_err(|_| AppError::ParseError {
        message: format!("invalid pileup position: {}", parts[1]),
    })?;
    let depth = parts[3].parse::<i32>().map_err(|_| AppError::ParseError {
        message: format!("invalid pileup depth: {}", parts[3]),
    })?;

    Ok(PileupRecord {
        chromosome: parts[0].to_string(),
        position,
        reference: parts[2].to_string(),
        depth,
        pileup: parts[4].to_string(),
        quality: parts[5].to_string(),
    })
}

fn compare_coordinates(left: &PileupRecord, right: &PileupRecord) -> Ordering {
    match compare_chromosomes(&left.chromosome, &right.chromosome) {
        Ordering::Equal => left.position.cmp(&right.position),
        other => other,
    }
}

fn compare_chromosomes(left: &str, right: &str) -> Ordering {
    let left_key = chromosome_sort_key(left);
    let right_key = chromosome_sort_key(right);
    left_key.cmp(&right_key)
}

fn chromosome_sort_key(chrom: &str) -> (u8, u32, &str) {
    let raw = chrom.trim();
    let normalized = raw.strip_prefix("chr").unwrap_or(raw);
    let upper = normalized.to_ascii_uppercase();

    if let Ok(num) = upper.parse::<u32>() {
        return (0, num, raw);
    }

    match upper.as_str() {
        "X" => (1, 23, raw),
        "Y" => (1, 24, raw),
        "M" | "MT" => (1, 25, raw),
        _ => (2, 0, raw),
    }
}

fn advance_to_target(
    current: &mut Option<PileupRecord>,
    stream: &mut PileupStream,
    target: &PileupRecord,
) -> Result<()> {
    loop {
        let Some(record) = current.as_ref() else {
            return Ok(());
        };
        if compare_coordinates(record, target) == Ordering::Less {
            *current = stream.next_record()?;
            continue;
        }
        return Ok(());
    }
}

#[derive(Debug, Clone)]
struct GcInterval {
    start: i32,
    end: i32,
    gc: Arc<str>,
}

#[cfg(test)]
fn parse_gc_intervals(lines: Vec<String>) -> HashMap<String, Vec<GcInterval>> {
    parse_gc_intervals_filtered(lines, None)
}

fn parse_gc_intervals_from_file(
    path: &str,
    chromosome_filter: Option<&HashSet<String>>,
) -> Result<HashMap<String, Vec<GcInterval>>> {
    if path.ends_with(".gz") {
        let file = File::open(path)?;
        let decoder = GzDecoder::new(file);
        parse_gc_intervals_from_reader(BufReader::new(decoder), chromosome_filter)
    } else {
        let file = File::open(path)?;
        parse_gc_intervals_from_reader(BufReader::new(file), chromosome_filter)
    }
}

fn parse_gc_intervals_from_reader<R: BufRead>(
    mut reader: R,
    chromosome_filter: Option<&HashSet<String>>,
) -> Result<HashMap<String, Vec<GcInterval>>> {
    let mut current_chromosome: Option<String> = None;
    let mut current_span: i32 = 0;
    let mut keep_current = false;
    let mut map: HashMap<String, Vec<GcInterval>> = HashMap::new();
    let mut gc_values: HashMap<String, Arc<str>> = HashMap::new();

    let mut buf = String::new();
    loop {
        buf.clear();
        let read = reader.read_line(&mut buf)?;
        if read == 0 {
            break;
        }

        let line = buf.trim_end_matches(['\n', '\r']);
        if line.starts_with("variableStep") {
            let mut chromosome = None;
            let mut span = None;
            for token in line.split_whitespace() {
                if let Some(value) = token.strip_prefix("chrom=") {
                    chromosome = Some(value.to_string());
                }
                if let Some(value) = token.strip_prefix("span=") {
                    span = value.parse::<i32>().ok();
                }
            }
            current_chromosome = chromosome;
            current_span = span.unwrap_or(0);
            keep_current = current_chromosome
                .as_ref()
                .is_some_and(|chr| chromosome_filter.is_none_or(|filter| filter.contains(chr)));
            continue;
        }

        if !keep_current {
            continue;
        }

        let Some(chromosome) = &current_chromosome else {
            continue;
        };
        let mut fields = line.split_whitespace();
        let Some(start_str) = fields.next() else {
            continue;
        };
        let Some(gc) = fields.next() else {
            continue;
        };
        let Some(start) = start_str.parse::<i32>().ok() else {
            continue;
        };
        let interned_gc = intern_gc_value(&mut gc_values, gc);

        map.entry(chromosome.clone()).or_default().push(GcInterval {
            start,
            end: start + current_span,
            gc: interned_gc,
        });
    }
    Ok(map)
}

#[cfg(test)]
fn parse_gc_intervals_filtered(
    lines: Vec<String>,
    chromosome_filter: Option<&HashSet<String>>,
) -> HashMap<String, Vec<GcInterval>> {
    let mut current_chromosome: Option<String> = None;
    let mut current_span: i32 = 0;
    let mut map: HashMap<String, Vec<GcInterval>> = HashMap::new();
    let mut keep_current = false;
    let mut gc_values: HashMap<String, Arc<str>> = HashMap::new();

    for line in lines {
        if line.starts_with("variableStep") {
            let mut chromosome = None;
            let mut span = None;
            for token in line.split_whitespace() {
                if let Some(value) = token.strip_prefix("chrom=") {
                    chromosome = Some(value.to_string());
                }
                if let Some(value) = token.strip_prefix("span=") {
                    span = value.parse::<i32>().ok();
                }
            }
            current_chromosome = chromosome;
            current_span = span.unwrap_or(0);
            keep_current = current_chromosome
                .as_ref()
                .is_some_and(|chr| chromosome_filter.is_none_or(|filter| filter.contains(chr)));
            continue;
        }

        if !keep_current {
            continue;
        }

        let Some(chromosome) = &current_chromosome else {
            continue;
        };
        let mut fields = line.split_whitespace();
        let Some(start_str) = fields.next() else {
            continue;
        };
        let Some(gc) = fields.next() else {
            continue;
        };
        let Some(start) = start_str.parse::<i32>().ok() else {
            continue;
        };
        let interned_gc = intern_gc_value(&mut gc_values, gc);

        map.entry(chromosome.clone()).or_default().push(GcInterval {
            start,
            end: start + current_span,
            gc: interned_gc,
        });
    }

    map
}

fn intern_gc_value(cache: &mut HashMap<String, Arc<str>>, value: &str) -> Arc<str> {
    if let Some(existing) = cache.get(value) {
        return Arc::clone(existing);
    }

    let key = value.to_string();
    let arc: Arc<str> = Arc::from(key.as_str());
    cache.insert(key, Arc::clone(&arc));
    arc
}

fn region_chromosome(region: &str) -> &str {
    region.split(':').next().unwrap_or(region)
}

fn gc_value_at<'a>(
    map: &'a HashMap<String, Vec<GcInterval>>,
    chromosome: &str,
    position: i32,
) -> Option<&'a str> {
    let intervals = map.get(chromosome)?;
    let index = intervals.partition_point(|interval| interval.start <= position);
    if index == 0 {
        return None;
    }
    let interval = &intervals[index - 1];
    (position < interval.end).then_some(interval.gc.as_ref())
}

fn output_for_region(output: &str, region: &str) -> String {
    let sanitized_region = region.replace([':', '-'], "_");
    if output == "-" {
        return output.to_string();
    }
    let (prefix, extension) = split_ext_like_python(output);
    format!("{prefix}_{sanitized_region}{extension}")
}

fn split_ext_like_python(path: &str) -> (String, String) {
    for compound in [".seqz.gz", ".txt.gz", ".bz2", ".gz", ".txt"] {
        if let Some(prefix) = path.strip_suffix(compound) {
            return (prefix.to_string(), compound.to_string());
        }
    }
    if let Some((prefix, suffix)) = path.rsplit_once('.') {
        return (prefix.to_string(), format!(".{suffix}"));
    }
    (path.to_string(), String::new())
}

#[cfg(test)]
mod tests {
    use crate::cli::parse_args;

    #[test]
    fn build_plan_for_parallel_outputs() {
        let args = parse_args([
            "bam2seqz_rs",
            "-n",
            "n.bam",
            "-t",
            "t.bam",
            "-gc",
            "gc.wig",
            "-F",
            "ref.fa",
            "-C",
            "7",
            "12",
            "--parallel",
            "2",
            "-o",
            "out.seqz.gz",
        ])
        .expect("expected parse success");

        let plan = super::build_plan(&args).expect("expected plan generation success");
        assert_eq!(plan.output_paths.len(), 2);
        assert!(plan.output_paths.iter().any(|path| path == "out_7.seqz.gz"));
        assert!(
            plan.output_paths
                .iter()
                .any(|path| path == "out_12.seqz.gz")
        );
    }

    #[test]
    fn build_plan_sets_ascii_quality_threshold() {
        let args = parse_args([
            "bam2seqz_rs",
            "-n",
            "n.bam",
            "-t",
            "t.bam",
            "-gc",
            "gc.wig",
            "-F",
            "ref.fa",
            "-f",
            "illumina",
        ])
        .expect("expected parse success");

        let plan = super::build_plan(&args).expect("expected plan generation success");
        assert_eq!(plan.qlimit_ascii, 84);
    }

    #[test]
    fn parses_gc_intervals_and_lookup() {
        let lines = vec![
            "variableStep chrom=7 span=50".to_string(),
            "101\t26".to_string(),
            "151\t28".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        assert_eq!(super::gc_value_at(&map, "7", 101), Some("26"));
        assert_eq!(super::gc_value_at(&map, "7", 149), Some("26"));
        assert_eq!(super::gc_value_at(&map, "7", 151), Some("28"));
    }

    #[test]
    fn parses_gc_intervals_with_trailing_whitespace() {
        let lines = vec![
            "variableStep chrom=12 span=50".to_string(),
            "101\t40\t".to_string(),
            "151 42".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        assert_eq!(super::gc_value_at(&map, "12", 101), Some("40"));
        assert_eq!(super::gc_value_at(&map, "12", 151), Some("42"));
    }
}
