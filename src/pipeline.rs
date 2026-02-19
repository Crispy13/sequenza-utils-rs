use crate::cli::{Bam2SeqzArgs, BamBackend};
use crate::errors::{AppError, Result};
use crate::external_tools::{CommandStream, ExternalTools};
use crate::seqz_core::{SeqzParams, do_seqz};
use crate::writer;
use crossbeam_channel::unbounded;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
#[cfg(feature = "htslib-prototype")]
use rust_htslib::bam::Read as _;
use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::process::Child;
use std::sync::Arc;
use std::time::Duration;
use tempfile::{Builder, NamedTempFile};
use tracing::info;

const AUTO_BIN_SIZE_BP: i32 = 5_000_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PipelinePlan {
    pub use_bam_input: bool,
    pub regions: Vec<String>,
    pub output_paths: Vec<String>,
    pub qlimit_ascii: i32,
}

#[derive(Debug)]
struct RegionChunk {
    index: usize,
    region: String,
    chunk_file: NamedTempFile,
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
    info!(
        pileup = args.pileup,
        regions = args.chr.len(),
        parallel = args.nproc,
        output = %args.out,
        "starting pipeline run"
    );
    let plan = build_plan(args)?;
    let base_tools = ExternalTools::from_args(args);
    let gc_intervals = Arc::new(load_gc_intervals_for_run(args, &base_tools)?);

    if args.nproc > 1 {
        if should_auto_bin_parallel(args) {
            let auto_regions = derive_auto_binned_regions(args, &base_tools, gc_intervals.as_ref())?;
            info!(
                bins = auto_regions.len(),
                bin_size_bp = AUTO_BIN_SIZE_BP,
                "running auto-binned parallel mode with ordered single output"
            );
            return run_parallel_single_output(args, gc_intervals, auto_regions);
        }

        if args.parallel_single_output {
            let effective_regions = if args.chr.is_empty() {
                Vec::new()
            } else {
                let filtered = filter_regions_with_gc(&args.chr, gc_intervals.as_ref());
                let skipped_regions = args.chr.len().saturating_sub(filtered.len());
                if skipped_regions > 0 {
                    info!(
                        requested = args.chr.len(),
                        retained = filtered.len(),
                        skipped = skipped_regions,
                        "filtered regions with no GC overlap before mpileup"
                    );
                }
                filtered
            };

            if !args.chr.is_empty() && effective_regions.is_empty() {
                let tools = ExternalTools::from_args(args);
                return write_header_only_output(&args.out, &tools);
            }

            return run_parallel_single_output(args, gc_intervals, effective_regions);
        }

        let work = args
            .chr
            .iter()
            .enumerate()
            .map(|(index, region)| {
                let output_path = plan
                    .output_paths
                    .get(index)
                    .ok_or_else(|| AppError::ParseError {
                        message: "missing output path for region".to_string(),
                    })?
                    .clone();
                let overlaps_gc = region_overlaps_gc(gc_intervals.as_ref(), region);
                Ok((region.clone(), output_path, overlaps_gc))
            })
            .collect::<Result<Vec<_>>>()?;

        let skipped_regions = work.iter().filter(|(_, _, overlaps_gc)| !*overlaps_gc).count();
        if skipped_regions > 0 {
            info!(
                requested = work.len(),
                retained = work.len().saturating_sub(skipped_regions),
                skipped = skipped_regions,
                "filtered regions with no GC overlap before mpileup"
            );
        }

        let pool = ThreadPoolBuilder::new()
            .num_threads(args.nproc)
            .build()
            .map_err(|err| AppError::ParseError {
                message: format!("failed to initialize rayon thread pool: {err}"),
            })?;

        if !args.pileup && args.bam_backend == BamBackend::RustHtslib {
            #[cfg(feature = "htslib-prototype")]
            {
                pool.install(|| {
                    work.par_iter().try_for_each_init(
                        || crate::htslib_mpileup::HtslibWorkerContext::from_args(args),
                        |context, (region, output_path, overlaps_gc)| {
                            let tools = ExternalTools::from_args(args);
                            if *overlaps_gc {
                                run_one_with_htslib_context(
                                    args,
                                    &tools,
                                    context,
                                    std::slice::from_ref(region),
                                    output_path,
                                    gc_intervals.as_ref(),
                                )
                            } else {
                                write_header_only_output(output_path, &tools)
                            }
                        },
                    )
                })?;
            }
            #[cfg(not(feature = "htslib-prototype"))]
            {
                return Err(AppError::InvalidValue {
                    flag: "--bam-backend".to_string(),
                    value: "rust-htslib".to_string(),
                    reason:
                        "binary built without feature \"htslib-prototype\"; rebuild with --features htslib-prototype"
                            .to_string(),
                });
            }
        } else {
            pool.install(|| {
                work.par_iter().try_for_each(|(region, output_path, overlaps_gc)| {
                    let tools = ExternalTools::from_args(args);
                    if *overlaps_gc {
                        run_one(
                            args,
                            &tools,
                            std::slice::from_ref(region),
                            output_path,
                            gc_intervals.as_ref(),
                        )
                    } else {
                        write_header_only_output(output_path, &tools)
                    }
                })
            })?;
        }
    } else {
        run_one(args, &base_tools, &args.chr, &args.out, gc_intervals.as_ref())?;
    }

    Ok(())
}

fn run_parallel_single_output(
    args: &Bam2SeqzArgs,
    gc_intervals: Arc<HashMap<String, Vec<GcInterval>>>,
    regions: Vec<String>,
) -> Result<()> {
    if regions.is_empty() {
        let tools = ExternalTools::from_args(args);
        return write_header_only_output(&args.out, &tools);
    }

    let work = regions
        .iter()
        .enumerate()
        .map(|(index, region)| (index, region.clone()))
        .collect::<Vec<_>>();

    let pool = ThreadPoolBuilder::new()
        .num_threads(args.nproc)
        .build()
        .map_err(|err| AppError::ParseError {
            message: format!("failed to initialize rayon thread pool: {err}"),
        })?;
    let mut chunks: Vec<Option<RegionChunk>> = (0..work.len()).map(|_| None).collect();

    if !args.pileup && args.bam_backend == BamBackend::RustHtslib {
        #[cfg(feature = "htslib-prototype")]
        {
            let chunk_results = pool.install(|| {
                work.par_iter()
                    .map_init(
                        || crate::htslib_mpileup::HtslibWorkerContext::from_args(args),
                        |context, (index, region)| -> Result<RegionChunk> {
                            let tools = ExternalTools::from_args(args);
                            let chunk_file = Builder::new()
                                .prefix("bam2seqz_rs_parallel_chunk_")
                                .suffix(".seqz")
                                .tempfile_in("tmp")?;
                            let chunk_path = chunk_file.path().to_string_lossy().into_owned();
                            run_one_with_htslib_context(
                                args,
                                &tools,
                                context,
                                std::slice::from_ref(region),
                                &chunk_path,
                                gc_intervals.as_ref(),
                            )?;
                            Ok(RegionChunk {
                                index: *index,
                                region: region.clone(),
                                chunk_file,
                            })
                        },
                    )
                    .collect::<Vec<_>>()
            });

            for chunk_result in chunk_results {
                let region_chunk = chunk_result?;
                let index = region_chunk.index;
                chunks[index] = Some(region_chunk);
            }
        }
        #[cfg(not(feature = "htslib-prototype"))]
        {
            return Err(AppError::InvalidValue {
                flag: "--bam-backend".to_string(),
                value: "rust-htslib".to_string(),
                reason:
                    "binary built without feature \"htslib-prototype\"; rebuild with --features htslib-prototype"
                        .to_string(),
            });
        }
    } else {
        let (sender, receiver) = unbounded::<Result<RegionChunk>>();
        pool.install(|| {
            work.par_iter().for_each_with(sender, |tx, (index, region)| {
                let result = (|| {
                    let tools = ExternalTools::from_args(args);
                    let chunk_file = Builder::new()
                        .prefix("bam2seqz_rs_parallel_chunk_")
                        .suffix(".seqz")
                        .tempfile_in("tmp")?;
                    let chunk_path = chunk_file.path().to_string_lossy().into_owned();
                    run_one(
                        args,
                        &tools,
                        std::slice::from_ref(region),
                        &chunk_path,
                        gc_intervals.as_ref(),
                    )?;
                    Ok(RegionChunk {
                        index: *index,
                        region: region.clone(),
                        chunk_file,
                    })
                })();
                let _ = tx.send(result);
            });
        });

        for _ in 0..work.len() {
            let region_chunk = receiver.recv().map_err(|err| AppError::ParseError {
                message: format!("failed to receive parallel region result: {err}"),
            })??;
            let index = region_chunk.index;
            chunks[index] = Some(region_chunk);
        }
    }

    let tools = ExternalTools::from_args(args);
    writer::with_text_output_writer(&args.out, &tools, |out| {
        writer::write_seqz_header(out)?;

        for (index, maybe_chunk) in chunks.into_iter().enumerate() {
            let chunk = maybe_chunk.ok_or_else(|| AppError::ParseError {
                message: format!("missing region chunk result for index {index}"),
            })?;
            info!(
                region = %chunk.region,
                chunk = %chunk.chunk_file.path().to_string_lossy(),
                "merging ordered region chunk"
            );

            let mut reader = BufReader::new(chunk.chunk_file.reopen()?);
            let mut line = String::new();
            let mut skip_header = true;
            loop {
                line.clear();
                let read = reader.read_line(&mut line)?;
                if read == 0 {
                    break;
                }
                if skip_header {
                    skip_header = false;
                    continue;
                }
                out.write_all(line.as_bytes())?;
            }
        }
        Ok(())
    })?;

    if args.out.ends_with(".gz") {
        info!(output = %args.out, "indexing compressed seqz with tabix");
        tools.tabix_index_seqz(&args.out)?;
    }
    info!(output = %args.out, "completed parallel single-output merge");
    Ok(())
}

fn should_auto_bin_parallel(args: &Bam2SeqzArgs) -> bool {
    args.nproc > 1 && !args.pileup && !args.has_explicit_ranged_regions()
}

fn derive_auto_binned_regions(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
) -> Result<Vec<String>> {
    let chromosomes = if args.chr.is_empty() {
        tools.list_bam_chromosomes(&args.tumor)?
    } else {
        dedup_chromosomes_preserving_order(&args.chr)
    };

    let bins = auto_bin_regions_for_chromosomes(&chromosomes, gc_intervals, AUTO_BIN_SIZE_BP);
    Ok(filter_regions_with_gc(&bins, gc_intervals))
}

fn dedup_chromosomes_preserving_order(chromosomes: &[String]) -> Vec<String> {
    let mut seen = HashSet::new();
    let mut ordered = Vec::with_capacity(chromosomes.len());
    for chromosome in chromosomes {
        if seen.insert(chromosome.clone()) {
            ordered.push(chromosome.clone());
        }
    }
    ordered
}

fn auto_bin_regions_for_chromosomes(
    chromosomes: &[String],
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
    bin_size_bp: i32,
) -> Vec<String> {
    let mut regions = Vec::new();
    for chromosome in chromosomes {
        let Some(intervals) = gc_intervals.get(chromosome) else {
            continue;
        };
        let bins = gc_overlapping_bins(intervals, bin_size_bp);
        for (start, end) in bins {
            regions.push(format!("{chromosome}:{start}-{end}"));
        }
    }
    regions
}

fn gc_overlapping_bins(intervals: &[GcInterval], bin_size_bp: i32) -> Vec<(i32, i32)> {
    if bin_size_bp <= 0 {
        return Vec::new();
    }

    let mut bin_ids = BTreeSet::new();
    for interval in intervals {
        if interval.end <= interval.start {
            continue;
        }
        let start = interval.start.max(1);
        let end = (interval.end - 1).max(start);
        let start_bin = (start - 1) / bin_size_bp;
        let end_bin = (end - 1) / bin_size_bp;
        for bin_id in start_bin..=end_bin {
            let _ = bin_ids.insert(bin_id);
        }
    }

    bin_ids
        .into_iter()
        .map(|bin_id| {
            let start = bin_id * bin_size_bp + 1;
            let end = (bin_id + 1) * bin_size_bp;
            (start, end)
        })
        .collect()
}

fn run_one(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    regions: &[String],
    output: &str,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
) -> Result<()> {
    info!(
        output = %output,
        regions = %if regions.is_empty() { "all".to_string() } else { regions.join(",") },
        pileup_mode = args.pileup,
        "starting region pipeline"
    );
    let effective_regions = if regions.is_empty() {
        Vec::new()
    } else {
        filter_regions_with_gc(regions, gc_intervals)
    };

    if !regions.is_empty() {
        let skipped_regions = regions.len().saturating_sub(effective_regions.len());
        if skipped_regions > 0 {
            info!(
                requested = regions.len(),
                retained = effective_regions.len(),
                skipped = skipped_regions,
                "filtered regions with no GC overlap before mpileup"
            );
        }
        if effective_regions.is_empty() {
            return write_header_only_output(output, tools);
        }
    }

    let region_scope: &[String] = if regions.is_empty() {
        regions
    } else {
        effective_regions.as_slice()
    };

    if !args.pileup && args.bam_backend == BamBackend::RustHtslib {
        #[cfg(feature = "htslib-prototype")]
        {
            let mut context = crate::htslib_mpileup::HtslibWorkerContext::from_args(args);
            return run_one_with_htslib_context(
                args,
                tools,
                &mut context,
                region_scope,
                output,
                gc_intervals,
            );
        }
        #[cfg(not(feature = "htslib-prototype"))]
        {
            return Err(AppError::InvalidValue {
                flag: "--bam-backend".to_string(),
                value: "rust-htslib".to_string(),
                reason:
                    "binary built without feature \"htslib-prototype\"; rebuild with --features htslib-prototype"
                        .to_string(),
            });
        }
    }

    let mut tumor_stream = open_pileup_stream(args, tools, &args.tumor, region_scope)?;
    let mut alt_stream = if let Some(normal2_path) = &args.normal2 {
        Some(open_pileup_stream(args, tools, normal2_path, region_scope)?)
    } else {
        None
    };

    let mut normal_stream = open_pileup_stream(args, tools, &args.normal, region_scope)?;
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

    let mut progress = PipelineProgress::new(args.progress, output, region_scope);

    writer::with_text_output_writer(output, tools, move |out| {
        writer::write_seqz_header(out)?;

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

            let Some(gc) = gc_value_at(gc_intervals, &normal.chromosome, normal.position) else {
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
        info!(output = %output, "indexing compressed seqz with tabix");
        tools.tabix_index_seqz(output)?;
    }
    info!(output = %output, "completed region pipeline");
    Ok(())
}

#[cfg(feature = "htslib-prototype")]
fn run_one_with_htslib_context(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    context: &mut crate::htslib_mpileup::HtslibWorkerContext,
    regions: &[String],
    output: &str,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
) -> Result<()> {
    let params = SeqzParams {
        depth_sum: args.depth_sum,
        qlimit: args.qlimit_ascii().clamp(0, i32::from(u8::MAX)) as u8,
        hom_t: args.hom,
        het_t: args.het,
        het_f: args.het_f,
        het_only: false,
    };

    let mut progress = PipelineProgress::new(args.progress, output, regions);

    let mut resolved_regions = context.resolve_regions(regions)?;
    if regions.is_empty() {
        let filtered = filter_regions_with_gc(&resolved_regions, gc_intervals);
        let skipped_regions = resolved_regions.len().saturating_sub(filtered.len());
        if skipped_regions > 0 {
            info!(
                requested = resolved_regions.len(),
                retained = filtered.len(),
                skipped = skipped_regions,
                "filtered regions with no GC overlap before htslib processing"
            );
        }
        resolved_regions = filtered;
    }

    writer::with_text_output_writer(output, tools, |out| {
        writer::write_seqz_header(out)?;

        let mut normal_line = String::with_capacity(96);
        let mut tumor_line = String::with_capacity(96);
        let mut alt_line = String::with_capacity(96);

        for region in resolved_regions {
            let mut readers = context.readers_mut()?;
            let region_spec = crate::htslib_mpileup::parse_region_spec(readers.tumor, &region)?;
            let reference_bases =
                crate::htslib_mpileup::fetch_reference_bases(readers.fasta, &region_spec)?;

            crate::htslib_mpileup::fetch_region(readers.normal, &region_spec)?;
            crate::htslib_mpileup::fetch_region(readers.tumor, &region_spec)?;
            if let Some(reader) = readers.normal2.as_deref_mut() {
                crate::htslib_mpileup::fetch_region(reader, &region_spec)?;
            }

            let mut normal_iter = readers.normal.pileup();
            let mut tumor_iter = readers.tumor.pileup();
            let mut alt_iter = readers.normal2.as_deref_mut().map(|reader| reader.pileup());

            let mut tumor_current = crate::htslib_mpileup::next_record_from_pileups(
                &mut tumor_iter,
                &region_spec,
                &reference_bases,
            )?;
            let mut alt_current = if let Some(iter) = alt_iter.as_mut() {
                crate::htslib_mpileup::next_record_from_pileups(
                    iter,
                    &region_spec,
                    &reference_bases,
                )?
            } else {
                None
            };

            while let Some(normal) = crate::htslib_mpileup::next_record_from_pileups(
                &mut normal_iter,
                &region_spec,
                &reference_bases,
            )? {
                progress.on_processed(&normal.chromosome, normal.position);

                while let Some(tumor) = tumor_current.as_ref() {
                    if compare_htslib_coordinates(tumor, &normal) == Ordering::Less {
                        tumor_current = crate::htslib_mpileup::next_record_from_pileups(
                            &mut tumor_iter,
                            &region_spec,
                            &reference_bases,
                        )?;
                        continue;
                    }
                    break;
                }

                let Some(tumor) = tumor_current.as_ref() else {
                    break;
                };
                if compare_htslib_coordinates(tumor, &normal) != Ordering::Equal {
                    continue;
                }

                let Some(gc) = gc_value_at(gc_intervals, &normal.chromosome, normal.position)
                else {
                    continue;
                };

                normal.write_data_line(&mut normal_line);
                tumor.write_data_line(&mut tumor_line);

                let seqz_fields = if let Some(iter) = alt_iter.as_mut() {
                    while let Some(alt) = alt_current.as_ref() {
                        if compare_htslib_coordinates(alt, &normal) == Ordering::Less {
                            alt_current = crate::htslib_mpileup::next_record_from_pileups(
                                iter,
                                &region_spec,
                                &reference_bases,
                            )?;
                            continue;
                        }
                        break;
                    }

                    if let Some(alt) = alt_current.as_ref() {
                        if compare_htslib_coordinates(alt, &normal) == Ordering::Equal {
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
        }

        progress.finish();
        Ok(())
    })?;

    if output.ends_with(".gz") {
        info!(output = %output, "indexing compressed seqz with tabix");
        tools.tabix_index_seqz(output)?;
    }
    info!(output = %output, "completed region pipeline");
    Ok(())
}

#[cfg(feature = "htslib-prototype")]
fn compare_htslib_coordinates(
    left: &crate::htslib_mpileup::HtslibPileupRecord,
    right: &crate::htslib_mpileup::HtslibPileupRecord,
) -> Ordering {
    match compare_chromosomes(&left.chromosome, &right.chromosome) {
        Ordering::Equal => left.position.cmp(&right.position),
        other => other,
    }
}

fn write_header_only_output(output: &str, tools: &ExternalTools) -> Result<()> {
    info!(
        output = %output,
        "no requested regions overlap GC intervals; writing header-only output"
    );
    writer::with_text_output_writer(output, tools, |out| {
        writer::write_seqz_header(out)?;
        Ok(())
    })?;
    if output.ends_with(".gz") {
        info!(output = %output, "indexing compressed seqz with tabix");
        tools.tabix_index_seqz(output)?;
    }
    info!(output = %output, "completed region pipeline");
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
        info!(input = %input_path, regions = regions.len(), "opening pileup input stream");
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
        match args.bam_backend {
            BamBackend::Samtools => {
                info!(
                    input = %input_path,
                    regions = regions.len(),
                    backend = "samtools",
                    "opening BAM mpileup stream"
                );
                PileupStream::from_bam_mpileup_stream(tools, input_path, fasta, regions)
            }
            BamBackend::RustHtslib => Err(AppError::ParseError {
                message:
                    "internal error: rust-htslib backend must use direct worker-context execution"
                        .to_string(),
            }),
        }
    }
}

struct PileupStream {
    reader: Box<dyn BufRead>,
    line_buffer: Vec<u8>,
    _tempfile: Option<NamedTempFile>,
    process: Option<StreamProcess>,
    region_source: Option<RegionStreamSource>,
}

#[derive(Debug)]
struct StreamProcess {
    child: Child,
    stderr_capture: NamedTempFile,
    command: String,
    finished: bool,
}

#[derive(Debug, Clone)]
struct RegionStreamSource {
    tools: ExternalTools,
    bam: String,
    fasta: String,
    pending_regions: VecDeque<String>,
}

impl RegionStreamSource {
    fn new(tools: ExternalTools, bam: String, fasta: String, regions: &[String]) -> Self {
        Self {
            tools,
            bam,
            fasta,
            pending_regions: regions.iter().cloned().collect(),
        }
    }

    fn spawn_next(&mut self) -> Result<Option<CommandStream>> {
        let Some(region) = self.pending_regions.pop_front() else {
            return Ok(None);
        };
        info!(region = %region, "spawning next region mpileup stream");
        self.tools
            .spawn_samtools_mpileup_stream(&self.bam, &self.fasta, Some(region.as_str()))
            .map(Some)
    }
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
                line_buffer: Vec::with_capacity(256),
                _tempfile: None,
                process: None,
                region_source: None,
            })
        } else {
            let file = File::open(path)?;
            Ok(Self {
                reader: Box::new(BufReader::new(file)),
                line_buffer: Vec::with_capacity(256),
                _tempfile: None,
                process: None,
                region_source: None,
            })
        }
    }

    fn from_tempfile(tempfile: NamedTempFile) -> Result<Self> {
        let reader_file = tempfile.reopen()?;
        Ok(Self {
            reader: Box::new(BufReader::new(reader_file)),
            line_buffer: Vec::with_capacity(256),
            _tempfile: Some(tempfile),
            process: None,
            region_source: None,
        })
    }

    fn from_bam_mpileup_stream(
        tools: &ExternalTools,
        bam: &str,
        fasta: &str,
        regions: &[String],
    ) -> Result<Self> {
        if regions.is_empty() {
            let stream = tools.spawn_samtools_mpileup_stream(bam, fasta, None)?;
            return Self::from_command_stream(stream);
        }

        let mut source =
            RegionStreamSource::new(tools.clone(), bam.to_string(), fasta.to_string(), regions);
        let first = source.spawn_next()?.ok_or_else(|| AppError::ParseError {
            message: "failed to initialize region mpileup stream".to_string(),
        })?;
        let mut stream = Self::from_command_stream(first)?;
        stream.region_source = Some(source);
        Ok(stream)
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
            line_buffer: Vec::with_capacity(256),
            _tempfile: None,
            process: Some(process),
            region_source: None,
        })
    }

    fn set_command_stream(&mut self, stream: CommandStream) {
        self.reader = Box::new(BufReader::new(stream.stdout));
        self.process = Some(StreamProcess {
            child: stream.child,
            stderr_capture: stream.stderr_capture,
            command: stream.command,
            finished: false,
        });
        self.line_buffer.clear();
    }

    fn advance_region_stream(&mut self) -> Result<bool> {
        let Some(source) = self.region_source.as_mut() else {
            return Ok(false);
        };
        let Some(next_stream) = source.spawn_next()? else {
            return Ok(false);
        };
        self.set_command_stream(next_stream);
        Ok(true)
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
        info!(command = %process.command, code = ?status.code(), "mpileup stream completed");
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
            let read = self.reader.read_until(b'\n', &mut self.line_buffer)?;
            if read == 0 {
                self.ensure_process_completed()?;
                if self.advance_region_stream()? {
                    continue;
                }
                return Ok(None);
            }
            let line = trim_line_end(&self.line_buffer);
            if line.is_empty() || line.iter().all(|byte| byte.is_ascii_whitespace()) {
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

fn parse_pileup_record(line: &[u8]) -> Result<PileupRecord> {
    let parts = line.splitn(6, |byte| *byte == b'\t').collect::<Vec<_>>();
    if parts.len() < 6 {
        return Err(AppError::ParseError {
            message: format!("invalid pileup line: {}", String::from_utf8_lossy(line)),
        });
    }

    let position = parse_i32_ascii(parts[1]).ok_or_else(|| AppError::ParseError {
        message: format!(
            "invalid pileup position: {}",
            String::from_utf8_lossy(parts[1])
        ),
    })?;
    let depth = parse_i32_ascii(parts[3]).ok_or_else(|| AppError::ParseError {
        message: format!("invalid pileup depth: {}", String::from_utf8_lossy(parts[3])),
    })?;

    Ok(PileupRecord {
        chromosome: String::from_utf8_lossy(parts[0]).into_owned(),
        position,
        reference: String::from_utf8_lossy(parts[2]).into_owned(),
        depth,
        pileup: String::from_utf8_lossy(parts[4]).into_owned(),
        quality: String::from_utf8_lossy(parts[5]).into_owned(),
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
    gc: String,
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

    let mut buf = Vec::with_capacity(128);
    loop {
        buf.clear();
        let read = reader.read_until(b'\n', &mut buf)?;
        if read == 0 {
            break;
        }

        let line = trim_line_end(&buf);
        if line.is_empty() {
            continue;
        }

        if line.starts_with(b"variableStep") {
            let (chromosome, span) = parse_variable_step_header(line);
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
        let Some((start, gc)) = parse_gc_data_line(line) else {
            continue;
        };

        map.entry(chromosome.clone()).or_default().push(GcInterval {
            start,
            end: start + current_span,
            gc,
        });
    }
    Ok(map)
}

fn load_gc_intervals_for_run(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
) -> Result<HashMap<String, Vec<GcInterval>>> {
    let mut gc_filter: HashSet<String> = HashSet::new();
    if !args.chr.is_empty() {
        gc_filter.extend(args.chr.iter().map(|region| region_chromosome(region).to_string()));
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
    let gc_interval_count = gc_intervals.values().map(Vec::len).sum::<usize>();
    info!(
        chromosomes = gc_intervals.len(),
        intervals = gc_interval_count,
        "loaded gc intervals"
    );
    Ok(gc_intervals)
}

fn trim_line_end(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 {
        let value = line[end - 1];
        if value == b'\n' || value == b'\r' {
            end -= 1;
        } else {
            break;
        }
    }
    &line[..end]
}

fn parse_variable_step_header(line: &[u8]) -> (Option<String>, Option<i32>) {
    let mut chromosome = None;
    let mut span = None;

    for token in line
        .split(|byte| byte.is_ascii_whitespace())
        .filter(|token| !token.is_empty())
    {
        if let Some(value) = token.strip_prefix(b"chrom=") {
            chromosome = Some(String::from_utf8_lossy(value).into_owned());
        }
        if let Some(value) = token.strip_prefix(b"span=") {
            span = parse_i32_ascii(value);
        }
    }

    (chromosome, span)
}

fn parse_gc_data_line(line: &[u8]) -> Option<(i32, String)> {
    let mut fields = line
        .split(|byte| byte.is_ascii_whitespace())
        .filter(|token| !token.is_empty());
    let start = parse_i32_ascii(fields.next()?)?;
    let gc = String::from_utf8_lossy(fields.next()?).into_owned();
    Some((start, gc))
}

fn parse_i32_ascii(bytes: &[u8]) -> Option<i32> {
    std::str::from_utf8(bytes)
        .ok()
        .and_then(|value| value.parse::<i32>().ok())
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

        map.entry(chromosome.clone()).or_default().push(GcInterval {
            start,
            end: start + current_span,
            gc: gc.to_string(),
        });
    }

    map
}

fn filter_regions_with_gc(
    regions: &[String],
    gc_map: &HashMap<String, Vec<GcInterval>>,
) -> Vec<String> {
    regions
        .iter()
        .filter(|region| region_overlaps_gc(gc_map, region))
        .cloned()
        .collect()
}

fn region_overlaps_gc(gc_map: &HashMap<String, Vec<GcInterval>>, region: &str) -> bool {
    let chromosome = region_chromosome(region);
    let Some(intervals) = gc_map.get(chromosome) else {
        return false;
    };
    if intervals.is_empty() {
        return false;
    }

    let Some((start, end)) = parse_region_bounds(region) else {
        return true;
    };

    let index = intervals.partition_point(|interval| interval.start <= end);
    intervals[..index]
        .iter()
        .rev()
        .take_while(|interval| interval.end > start)
        .any(|interval| interval.start <= end)
}

fn parse_region_bounds(region: &str) -> Option<(i32, i32)> {
    let (_, bounds) = region.split_once(':')?;
    let (start_raw, end_raw) = bounds.split_once('-')?;
    let start = start_raw.replace(',', "").parse::<i32>().ok()?;
    let end = end_raw.replace(',', "").parse::<i32>().ok()?;
    (end >= start).then_some((start, end))
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
    (position < interval.end).then_some(interval.gc.as_str())
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

    #[test]
    fn filters_regions_without_gc_overlap() {
        let lines = vec![
            "variableStep chrom=chr20 span=50".to_string(),
            "100\t48".to_string(),
            "200\t50".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        let regions = vec![
            "chr20:100-140".to_string(),
            "chr20:151-190".to_string(),
            "chr21:100-140".to_string(),
            "chr20:bad-range".to_string(),
        ];

        let filtered = super::filter_regions_with_gc(&regions, &map);
        assert_eq!(
            filtered,
            vec![
                "chr20:100-140".to_string(),
                "chr20:bad-range".to_string(),
            ]
        );
    }

    #[test]
    fn auto_bins_gc_intervals_for_requested_chromosomes() {
        let lines = vec![
            "variableStep chrom=chr20 span=50".to_string(),
            "100\t48".to_string(),
            "4_999_960\t49".replace('_', ""),
            "5_000_100\t50".to_string(),
            "variableStep chrom=chr21 span=50".to_string(),
            "100\t51".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        let chromosomes = vec!["chr20".to_string(), "chr21".to_string()];

        let bins = super::auto_bin_regions_for_chromosomes(&chromosomes, &map, 5_000_000);
        assert_eq!(
            bins,
            vec![
                "chr20:1-5000000".to_string(),
                "chr20:5000001-10000000".to_string(),
                "chr21:1-5000000".to_string(),
            ]
        );
    }

    #[test]
    fn auto_bin_respects_input_chromosome_order() {
        let lines = vec![
            "variableStep chrom=chr1 span=50".to_string(),
            "100\t40".to_string(),
            "variableStep chrom=chr2 span=50".to_string(),
            "100\t41".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        let chromosomes = vec!["chr2".to_string(), "chr1".to_string()];

        let bins = super::auto_bin_regions_for_chromosomes(&chromosomes, &map, 5_000_000);
        assert_eq!(
            bins,
            vec!["chr2:1-5000000".to_string(), "chr1:1-5000000".to_string()]
        );
    }

    #[test]
    fn dedup_chromosomes_preserves_first_seen_order() {
        let chromosomes = vec![
            "chr20".to_string(),
            "chr20".to_string(),
            "chr21".to_string(),
            "chr20".to_string(),
            "chr1".to_string(),
            "chr21".to_string(),
        ];

        let deduped = super::dedup_chromosomes_preserving_order(&chromosomes);
        assert_eq!(
            deduped,
            vec![
                "chr20".to_string(),
                "chr21".to_string(),
                "chr1".to_string(),
            ]
        );
    }
}
