use crate::cli::{Bam2SeqzArgs, BamBackend};
use crate::errors::{AppError, Result};
use crate::external_tools::{CommandStream, ExternalTools};
use crate::seqz_core::{SeqzInput, SeqzParams, do_seqz_typed};
use crate::writer;
use crossbeam_channel::{Receiver, RecvTimeoutError, bounded};
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
#[cfg(feature = "htslib-prototype")]
use rust_htslib::bam::Read as _;
use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write, stdout};
use std::path::Path;
use std::process::Child;
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};
use tempfile::{Builder, NamedTempFile};
use tracing::{info, warn};

const AUTO_BIN_SIZE_BP: i32 = 5_000_000;
const AUTO_BIN_MIN_SIZE_BP: i32 = 500_000;
const AUTO_BIN_TARGET_TASKS_PER_THREAD: usize = 4;
const PARALLEL_CHUNK_CHANNEL_FACTOR: usize = 2;
const PARALLEL_CHUNK_INITIAL_CAPACITY: usize = 512 * 1024;
const PARALLEL_CHUNK_RECV_TIMEOUT: Duration = Duration::from_secs(5);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PipelinePlan {
    pub use_bam_input: bool,
    pub regions: Vec<String>,
    pub output_paths: Vec<String>,
    pub qlimit_ascii: i32,
}

#[derive(Debug)]
struct RegionChunkBuffer {
    index: usize,
    region: String,
    tempfile: NamedTempFile,
    elapsed: Duration,
}

#[derive(Debug)]
enum ParallelChunkMessage {
    Ready(RegionChunkBuffer),
}

#[derive(Debug)]
struct AutoBinnedRegions {
    regions: Vec<String>,
    bin_size_bp: i32,
    covered_bases: u64,
    target_bins: usize,
}

#[derive(Debug)]
struct SamtoolsWorkerContext {
    tools: ExternalTools,
    body: Vec<u8>,
}

impl SamtoolsWorkerContext {
    fn from_args(args: &Bam2SeqzArgs) -> Self {
        Self {
            tools: ExternalTools::from_args(args),
            body: Vec::with_capacity(PARALLEL_CHUNK_INITIAL_CAPACITY),
        }
    }
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
    validate_bam_fasta_chromosome_naming(args, &base_tools)?;
    let gc_intervals = Arc::new(load_gc_intervals_for_run(args, &base_tools)?);

    if args.nproc > 1 {
        if should_auto_bin_parallel(args) {
            let auto_regions =
                derive_auto_binned_regions(args, &base_tools, gc_intervals.as_ref())?;
            info!(
                bins = auto_regions.regions.len(),
                bin_size_bp = auto_regions.bin_size_bp,
                covered_bases = auto_regions.covered_bases,
                target_bins = auto_regions.target_bins,
                "running auto-binned parallel mode with ordered single output"
            );
            return run_parallel_single_output(args, gc_intervals, auto_regions.regions);
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

        let skipped_regions = work
            .iter()
            .filter(|(_, _, overlaps_gc)| !*overlaps_gc)
            .count();
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
                            run_parallel_region_output_with_htslib_context(
                                args,
                                &tools,
                                context,
                                region,
                                output_path,
                                *overlaps_gc,
                                gc_intervals.as_ref(),
                            )
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
                work.par_iter().try_for_each_init(
                    || SamtoolsWorkerContext::from_args(args),
                    |context, (region, output_path, overlaps_gc)| {
                        run_parallel_region_output(
                            args,
                            context,
                            region,
                            output_path,
                            *overlaps_gc,
                            gc_intervals.as_ref(),
                        )
                    },
                )
            })?;
        }
    } else {
        run_one(
            args,
            &base_tools,
            &args.chr,
            &args.out,
            gc_intervals.as_ref(),
        )?;
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
    let effective_parallelism = work.len().min(args.nproc);
    info!(
        regions = work.len(),
        requested_threads = args.nproc,
        effective_parallelism,
        "starting parallel single-output run"
    );
    info!(
        compression_backend = "external_bgzip_cli",
        worker_chunk_compression = true,
        "parallel chunk compression backend"
    );
    if effective_parallelism < args.nproc {
        info!(
            regions = work.len(),
            requested_threads = args.nproc,
            hint = "CPU utilization is bounded by available region tasks; use smaller bins or more regions to saturate threads",
            "parallel task shape may limit CPU saturation"
        );
    }

    let channel_capacity = (args.nproc * PARALLEL_CHUNK_CHANNEL_FACTOR).max(4);
    let (chunk_sender, chunk_receiver) = bounded::<ParallelChunkMessage>(channel_capacity);

    let writer_tools = ExternalTools::from_args(args);
    let writer_output = args.out.clone();
    let expected_chunks = work.len();
    let writer_started = Instant::now();
    let writer_handle = thread::spawn(move || {
        write_ordered_parallel_output(
            &writer_output,
            &writer_tools,
            expected_chunks,
            chunk_receiver,
        )
    });

    let workers_started = Instant::now();
    let worker_result = if !args.pileup && args.bam_backend == BamBackend::RustHtslib {
        #[cfg(feature = "htslib-prototype")]
        {
            pool.install(|| {
                work.par_iter().try_for_each_init(
                    || {
                        (
                            crate::htslib_mpileup::HtslibWorkerContext::from_args(args),
                            ExternalTools::from_args(args),
                            chunk_sender.clone(),
                        )
                    },
                    |(context, tools, tx), (index, region)| -> Result<()> {
                        let started = Instant::now();
                        let chunk_tempfile = run_one_with_htslib_context_to_bgzip_tempfile(
                            args,
                            tools,
                            context,
                            region,
                            gc_intervals.as_ref(),
                            *index == 0,
                        )?;
                        let chunk = RegionChunkBuffer {
                            index: *index,
                            region: region.clone(),
                            tempfile: chunk_tempfile,
                            elapsed: started.elapsed(),
                        };
                        tx.send(ParallelChunkMessage::Ready(chunk)).map_err(|err| {
                            AppError::ParseError {
                                message: format!(
                                    "failed to send parallel region chunk for {}: {err}",
                                    region
                                ),
                            }
                        })
                    },
                )
            })
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
            work.par_iter().try_for_each_init(
                || (SamtoolsWorkerContext::from_args(args), chunk_sender.clone()),
                |(context, tx), (index, region)| -> Result<()> {
                    let started = Instant::now();
                    let chunk_tempfile = run_one_region_to_bgzip_tempfile_with_context(
                        args,
                        context,
                        region,
                        gc_intervals.as_ref(),
                        *index == 0,
                    )?;
                    let chunk = RegionChunkBuffer {
                        index: *index,
                        region: region.clone(),
                        tempfile: chunk_tempfile,
                        elapsed: started.elapsed(),
                    };
                    tx.send(ParallelChunkMessage::Ready(chunk)).map_err(|err| {
                        AppError::ParseError {
                            message: format!(
                                "failed to send parallel region chunk for {}: {err}",
                                region
                            ),
                        }
                    })
                },
            )
        })
    };
    let workers_elapsed = workers_started.elapsed();
    drop(chunk_sender);

    let writer_result = writer_handle.join().map_err(|_| AppError::ParseError {
        message: "parallel writer thread panicked".to_string(),
    })?;
    let writer_elapsed = writer_started.elapsed();

    if let Err(error) = worker_result {
        if let Err(writer_error) = writer_result {
            info!(
                worker_error = %error,
                writer_error = %writer_error,
                "parallel run failed in both worker and writer"
            );
        }
        return Err(error);
    }
    writer_result?;

    info!(
        output = %args.out,
        regions = work.len(),
        worker_elapsed_ms = workers_elapsed.as_millis(),
        writer_elapsed_ms = writer_elapsed.as_millis(),
        "completed parallel single-output merge"
    );
    Ok(())
}

fn write_ordered_parallel_output(
    output: &str,
    tools: &ExternalTools,
    expected_chunks: usize,
    chunk_receiver: Receiver<ParallelChunkMessage>,
) -> Result<()> {
    let emit_compressed = output.ends_with(".gz");

    if output == "-" {
        let mut out = stdout().lock();
        write_ordered_parallel_chunks(&mut out, expected_chunks, emit_compressed, &chunk_receiver)?;
        out.flush()?;
    } else {
        let mut out = BufWriter::new(File::create(output)?);
        write_ordered_parallel_chunks(&mut out, expected_chunks, emit_compressed, &chunk_receiver)?;
        out.flush()?;
    }

    if emit_compressed {
        info!(output = %output, "indexing compressed seqz with tabix");
        tools.tabix_index_seqz(output)?;
    }

    Ok(())
}

fn write_ordered_parallel_chunks(
    out: &mut dyn Write,
    expected_chunks: usize,
    emit_compressed: bool,
    chunk_receiver: &Receiver<ParallelChunkMessage>,
) -> Result<()> {
    let mut next_expected = 0usize;
    let mut pending = HashMap::new();
    let mut wait_cycles = 0u64;

    while next_expected < expected_chunks {
        let message = match chunk_receiver.recv_timeout(PARALLEL_CHUNK_RECV_TIMEOUT) {
            Ok(message) => {
                wait_cycles = 0;
                message
            }
            Err(RecvTimeoutError::Timeout) => {
                wait_cycles += 1;
                warn!(
                    next_expected,
                    expected_chunks,
                    pending_chunks = pending.len(),
                    waited_seconds = wait_cycles * PARALLEL_CHUNK_RECV_TIMEOUT.as_secs(),
                    "writer still waiting for next chunk"
                );
                continue;
            }
            Err(RecvTimeoutError::Disconnected) => {
                return Err(AppError::ParseError {
                    message: format!(
                        "parallel writer channel closed before all chunks were received (next={next_expected}, expected={expected_chunks})"
                    ),
                });
            }
        };

        let ParallelChunkMessage::Ready(chunk) = message;
        let index = chunk.index;
        let _ = pending.insert(index, chunk);

        while let Some(ordered_chunk) = pending.remove(&next_expected) {
            let chunk_bytes = ordered_chunk.tempfile.as_file().metadata()?.len();
            info!(
                index = ordered_chunk.index,
                region = %ordered_chunk.region,
                chunk_bytes,
                worker_elapsed_ms = ordered_chunk.elapsed.as_millis(),
                "merging ordered region chunk"
            );

            if emit_compressed {
                append_file_to_writer(ordered_chunk.tempfile.path(), out)?;
            } else {
                append_bgzip_as_text_to_writer(ordered_chunk.tempfile.path(), out)?;
            }
            next_expected += 1;
        }
    }

    Ok(())
}

fn append_file_to_writer(path: &std::path::Path, out: &mut dyn Write) -> Result<()> {
    let mut input = File::open(path)?;
    std::io::copy(&mut input, out)?;
    Ok(())
}

fn append_bgzip_as_text_to_writer(path: &std::path::Path, out: &mut dyn Write) -> Result<()> {
    let input = File::open(path)?;
    let mut decoder = GzDecoder::new(input);
    std::io::copy(&mut decoder, out)?;
    Ok(())
}

fn run_parallel_region_output(
    args: &Bam2SeqzArgs,
    context: &mut SamtoolsWorkerContext,
    region: &str,
    output: &str,
    overlaps_gc: bool,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
) -> Result<()> {
    info!(
        output = %output,
        regions = %region,
        pileup_mode = args.pileup,
        "starting region pipeline"
    );

    let chunk_tempfile = if overlaps_gc {
        run_one_region_to_bgzip_tempfile_with_context(args, context, region, gc_intervals, true)?
    } else {
        write_bgzip_chunk_to_tempfile(&context.tools, true, &[])?
    };

    finalize_region_chunk_output(chunk_tempfile, output, &context.tools)?;
    info!(output = %output, "completed region pipeline");
    Ok(())
}

#[cfg(feature = "htslib-prototype")]
fn run_parallel_region_output_with_htslib_context(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    context: &mut crate::htslib_mpileup::HtslibWorkerContext,
    region: &str,
    output: &str,
    overlaps_gc: bool,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
) -> Result<()> {
    info!(
        output = %output,
        regions = %region,
        pileup_mode = args.pileup,
        "starting region pipeline"
    );

    let chunk_tempfile = if overlaps_gc {
        run_one_with_htslib_context_to_bgzip_tempfile(
            args,
            tools,
            context,
            region,
            gc_intervals,
            true,
        )?
    } else {
        write_bgzip_chunk_to_tempfile(tools, true, &[])?
    };

    finalize_region_chunk_output(chunk_tempfile, output, tools)?;
    info!(output = %output, "completed region pipeline");
    Ok(())
}

fn run_one_region_to_bgzip_tempfile_with_context(
    args: &Bam2SeqzArgs,
    context: &mut SamtoolsWorkerContext,
    region: &str,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
    include_header: bool,
) -> Result<NamedTempFile> {
    context.body.clear();
    run_one_region_body_to_buffer(
        args,
        &context.tools,
        region,
        gc_intervals,
        &mut context.body,
    )?;
    write_bgzip_chunk_to_tempfile(&context.tools, include_header, &context.body)
}

#[cfg(feature = "htslib-prototype")]
fn run_one_with_htslib_context_to_bgzip_tempfile(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    context: &mut crate::htslib_mpileup::HtslibWorkerContext,
    region: &str,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
    include_header: bool,
) -> Result<NamedTempFile> {
    let mut body = Vec::with_capacity(PARALLEL_CHUNK_INITIAL_CAPACITY);
    run_one_with_htslib_context_body_to_buffer(
        args,
        tools,
        context,
        region,
        gc_intervals,
        &mut body,
    )?;
    write_bgzip_chunk_to_tempfile(tools, include_header, &body)
}

fn write_bgzip_chunk_to_tempfile(
    tools: &ExternalTools,
    include_header: bool,
    body: &[u8],
) -> Result<NamedTempFile> {
    let temp = Builder::new()
        .prefix("bam2seqz_parallel_chunk_")
        .suffix(".seqz.gz")
        .tempfile_in("tmp")?;
    let temp_path = temp.path().to_string_lossy().into_owned();

    writer::with_text_output_writer(&temp_path, tools, |out| {
        if include_header {
            writer::write_seqz_header(out)?;
        }
        out.write_all(body)?;
        Ok(())
    })?;

    Ok(temp)
}

fn finalize_region_chunk_output(
    chunk_tempfile: NamedTempFile,
    output: &str,
    tools: &ExternalTools,
) -> Result<()> {
    if output.ends_with(".gz") {
        persist_tempfile_to_path(chunk_tempfile, output)?;
        info!(output = %output, "indexing compressed seqz with tabix");
        tools.tabix_index_seqz(output)?;
        return Ok(());
    }

    if output == "-" {
        let mut out = stdout().lock();
        append_bgzip_as_text_to_writer(chunk_tempfile.path(), &mut out)?;
        out.flush()?;
    } else {
        let mut out = BufWriter::new(File::create(output)?);
        append_bgzip_as_text_to_writer(chunk_tempfile.path(), &mut out)?;
        out.flush()?;
    }

    Ok(())
}

fn persist_tempfile_to_path(chunk_tempfile: NamedTempFile, output: &str) -> Result<()> {
    match chunk_tempfile.persist(output) {
        Ok(_) => Ok(()),
        Err(error) => {
            let file = error.file;
            let mut input = File::open(file.path())?;
            let mut out = BufWriter::new(File::create(output)?);
            std::io::copy(&mut input, &mut out)?;
            out.flush()?;
            fs::remove_file(file.path())?;
            Ok(())
        }
    }
}

fn run_one_region_body_to_buffer(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    region: &str,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
    out: &mut dyn Write,
) -> Result<()> {
    if !region_overlaps_gc(gc_intervals, region) {
        return Ok(());
    }

    let region_scope = [region.to_string()];
    let gc_target_bed = create_gc_target_bed_for_regions(&region_scope, gc_intervals)?;
    let gc_target_path = gc_target_bed
        .as_ref()
        .map(|bed| bed.path().to_string_lossy().into_owned());

    let mut tumor_stream = open_pileup_stream(
        args,
        tools,
        &args.tumor,
        &region_scope,
        gc_target_path.as_deref(),
    )?;
    let mut alt_stream = if let Some(normal2_path) = &args.normal2 {
        Some(open_pileup_stream(
            args,
            tools,
            normal2_path,
            &region_scope,
            gc_target_path.as_deref(),
        )?)
    } else {
        None
    };

    let mut normal_stream = open_pileup_stream(
        args,
        tools,
        &args.normal,
        &region_scope,
        gc_target_path.as_deref(),
    )?;
    let mut normal = PileupRecord::default();
    let mut tumor_current = PileupRecord::default();
    let mut tumor_present = tumor_stream.next_record_into(&mut tumor_current)?;
    let mut alt_current = PileupRecord::default();
    let mut alt_present = if let Some(stream) = alt_stream.as_mut() {
        stream.next_record_into(&mut alt_current)?
    } else {
        false
    };

    let params = SeqzParams {
        depth_sum: args.depth_sum,
        qlimit: args.qlimit_ascii().clamp(0, i32::from(u8::MAX)) as u8,
        hom_t: args.hom,
        het_t: args.het,
        het_f: args.het_f,
        het_only: false,
    };

    let mut gc_cursor = GcLookupCursor::new(gc_intervals);

    while normal_stream.next_record_into(&mut normal)? {
        advance_to_target(
            &mut tumor_present,
            &mut tumor_current,
            &mut tumor_stream,
            &normal,
        )?;
        if !tumor_present {
            break;
        }
        if compare_coordinates(&tumor_current, &normal) != Ordering::Equal {
            continue;
        }

        let Some(gc) = gc_cursor.value_at(&normal.chromosome, normal.position) else {
            continue;
        };

        let normal2_depth_override = if let Some(stream) = alt_stream.as_mut() {
            advance_to_target(&mut alt_present, &mut alt_current, stream, &normal)?;
            if alt_present && compare_coordinates(&alt_current, &normal) == Ordering::Equal {
                Some(alt_current.depth)
            } else {
                None
            }
        } else {
            None
        };

        let seqz_fields = do_seqz_typed(
            &SeqzInput {
                reference: &normal.reference,
                normal_depth: normal.depth,
                normal_pileup: &normal.pileup,
                normal_quality: &normal.quality,
                tumor_depth: tumor_current.depth,
                tumor_pileup: &tumor_current.pileup,
                tumor_quality: &tumor_current.quality,
                gc,
                normal2_depth_override,
            },
            &params,
        );

        if let Some(seqz) = seqz_fields {
            let _ = write!(out, "{}\t{}", normal.chromosome, normal.position);
            for field in seqz {
                out.write_all(b"\t")?;
                out.write_all(field.as_bytes())?;
            }
            out.write_all(b"\n")?;
        }
    }

    Ok(())
}

#[cfg(feature = "htslib-prototype")]
fn run_one_with_htslib_context_body_to_buffer(
    args: &Bam2SeqzArgs,
    _tools: &ExternalTools,
    context: &mut crate::htslib_mpileup::HtslibWorkerContext,
    region: &str,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
    out: &mut dyn Write,
) -> Result<()> {
    let params = SeqzParams {
        depth_sum: args.depth_sum,
        qlimit: args.qlimit_ascii().clamp(0, i32::from(u8::MAX)) as u8,
        hom_t: args.hom,
        het_t: args.het,
        het_f: args.het_f,
        het_only: false,
    };

    let region_value = region.to_string();
    let mut resolved_regions = context.resolve_regions(std::slice::from_ref(&region_value))?;
    resolved_regions = filter_regions_with_gc(&resolved_regions, gc_intervals);
    if resolved_regions.is_empty() {
        return Ok(());
    }

    let mut gc_cursor = GcLookupCursor::new(gc_intervals);

    for region_spec_text in resolved_regions {
        let mut readers = context.readers_mut()?;
        let region_spec =
            crate::htslib_mpileup::parse_region_spec(readers.tumor, &region_spec_text)?;
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
            crate::htslib_mpileup::next_record_from_pileups(iter, &region_spec, &reference_bases)?
        } else {
            None
        };

        while let Some(normal) = crate::htslib_mpileup::next_record_from_pileups(
            &mut normal_iter,
            &region_spec,
            &reference_bases,
        )? {
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

            let Some(gc) = gc_cursor.value_at(&normal.chromosome, normal.position) else {
                continue;
            };

            let normal2_depth_override = if let Some(iter) = alt_iter.as_mut() {
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

                alt_current
                    .as_ref()
                    .filter(|alt| compare_htslib_coordinates(alt, &normal) == Ordering::Equal)
                    .map(|alt| alt.depth)
            } else {
                None
            };

            let seqz_fields = do_seqz_typed(
                &SeqzInput {
                    reference: &normal.reference,
                    normal_depth: normal.depth,
                    normal_pileup: &normal.pileup,
                    normal_quality: &normal.quality,
                    tumor_depth: tumor.depth,
                    tumor_pileup: &tumor.pileup,
                    tumor_quality: &tumor.quality,
                    gc,
                    normal2_depth_override,
                },
                &params,
            );

            if let Some(seqz) = seqz_fields {
                let _ = write!(out, "{}\t{}", normal.chromosome, normal.position);
                for field in seqz {
                    out.write_all(b"\t")?;
                    out.write_all(field.as_bytes())?;
                }
                out.write_all(b"\n")?;
            }
        }
    }

    Ok(())
}

fn should_auto_bin_parallel(args: &Bam2SeqzArgs) -> bool {
    args.nproc > 1 && !args.pileup && !args.has_explicit_ranged_regions()
}

fn derive_auto_binned_regions(
    args: &Bam2SeqzArgs,
    tools: &ExternalTools,
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
) -> Result<AutoBinnedRegions> {
    let chromosomes = if args.chr.is_empty() {
        tools.list_bam_chromosomes(&args.tumor)?
    } else {
        dedup_chromosomes_preserving_order(&args.chr)
    };

    let covered_bases = estimate_gc_covered_bases(&chromosomes, gc_intervals);
    let target_bins = (args.nproc * AUTO_BIN_TARGET_TASKS_PER_THREAD).max(args.nproc);
    let bin_size_bp = choose_auto_bin_size_bp(covered_bases, target_bins);

    let bins = auto_bin_regions_for_chromosomes(&chromosomes, gc_intervals, bin_size_bp);
    let regions = filter_regions_with_gc(&bins, gc_intervals);

    Ok(AutoBinnedRegions {
        regions,
        bin_size_bp,
        covered_bases,
        target_bins,
    })
}

fn estimate_gc_covered_bases(
    chromosomes: &[String],
    gc_intervals: &HashMap<String, Vec<GcInterval>>,
) -> u64 {
    chromosomes
        .iter()
        .filter_map(|chromosome| find_gc_intervals_for_chromosome(gc_intervals, chromosome))
        .flat_map(|intervals| intervals.iter())
        .filter(|interval| interval.end > interval.start)
        .map(|interval| u64::try_from(interval.end - interval.start).unwrap_or(0))
        .sum()
}

fn choose_auto_bin_size_bp(covered_bases: u64, target_bins: usize) -> i32 {
    if covered_bases == 0 || target_bins == 0 {
        return AUTO_BIN_SIZE_BP;
    }

    let target_bins_u64 = u64::try_from(target_bins).unwrap_or(1);
    let raw = covered_bases.div_ceil(target_bins_u64);
    let clamped = raw.clamp(
        u64::try_from(AUTO_BIN_MIN_SIZE_BP).unwrap_or(1),
        u64::try_from(AUTO_BIN_SIZE_BP).unwrap_or(5_000_000),
    );
    i32::try_from(clamped).unwrap_or(AUTO_BIN_SIZE_BP)
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
        let Some(intervals) = find_gc_intervals_for_chromosome(gc_intervals, chromosome) else {
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

    let gc_target_bed = create_gc_target_bed_for_regions(region_scope, gc_intervals)?;
    let gc_target_path = gc_target_bed
        .as_ref()
        .map(|bed| bed.path().to_string_lossy().into_owned());

    let mut tumor_stream = open_pileup_stream(
        args,
        tools,
        &args.tumor,
        region_scope,
        gc_target_path.as_deref(),
    )?;
    let mut alt_stream = if let Some(normal2_path) = &args.normal2 {
        Some(open_pileup_stream(
            args,
            tools,
            normal2_path,
            region_scope,
            gc_target_path.as_deref(),
        )?)
    } else {
        None
    };

    let mut normal_stream = open_pileup_stream(
        args,
        tools,
        &args.normal,
        region_scope,
        gc_target_path.as_deref(),
    )?;
    let mut normal = PileupRecord::default();
    let mut tumor_current = PileupRecord::default();
    let mut tumor_present = tumor_stream.next_record_into(&mut tumor_current)?;
    let mut alt_current = PileupRecord::default();
    let mut alt_present = if let Some(stream) = alt_stream.as_mut() {
        stream.next_record_into(&mut alt_current)?
    } else {
        false
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

        let mut gc_cursor = GcLookupCursor::new(gc_intervals);

        while normal_stream.next_record_into(&mut normal)? {
            progress.on_processed(&normal.chromosome, normal.position);
            advance_to_target(
                &mut tumor_present,
                &mut tumor_current,
                &mut tumor_stream,
                &normal,
            )?;
            if !tumor_present {
                break;
            }
            if compare_coordinates(&tumor_current, &normal) != Ordering::Equal {
                continue;
            }

            let Some(gc) = gc_cursor.value_at(&normal.chromosome, normal.position) else {
                continue;
            };

            let normal2_depth_override = if let Some(stream) = alt_stream.as_mut() {
                advance_to_target(&mut alt_present, &mut alt_current, stream, &normal)?;
                if alt_present && compare_coordinates(&alt_current, &normal) == Ordering::Equal {
                    Some(alt_current.depth)
                } else {
                    None
                }
            } else {
                None
            };

            let seqz_fields = do_seqz_typed(
                &SeqzInput {
                    reference: &normal.reference,
                    normal_depth: normal.depth,
                    normal_pileup: &normal.pileup,
                    normal_quality: &normal.quality,
                    tumor_depth: tumor_current.depth,
                    tumor_pileup: &tumor_current.pileup,
                    tumor_quality: &tumor_current.quality,
                    gc,
                    normal2_depth_override,
                },
                &params,
            );

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

        let mut gc_cursor = GcLookupCursor::new(gc_intervals);

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

                let Some(gc) = gc_cursor.value_at(&normal.chromosome, normal.position) else {
                    continue;
                };

                let normal2_depth_override = if let Some(iter) = alt_iter.as_mut() {
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

                    alt_current
                        .as_ref()
                        .filter(|alt| compare_htslib_coordinates(alt, &normal) == Ordering::Equal)
                        .map(|alt| alt.depth)
                } else {
                    None
                };

                let seqz_fields = do_seqz_typed(
                    &SeqzInput {
                        reference: &normal.reference,
                        normal_depth: normal.depth,
                        normal_pileup: &normal.pileup,
                        normal_quality: &normal.quality,
                        tumor_depth: tumor.depth,
                        tumor_pileup: &tumor.pileup,
                        tumor_quality: &tumor.quality,
                        gc,
                        normal2_depth_override,
                    },
                    &params,
                );

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
    gc_target_bed: Option<&str>,
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
                    gc_target_bed = gc_target_bed.is_some(),
                    backend = "samtools",
                    "opening BAM mpileup stream"
                );
                PileupStream::from_bam_mpileup_stream(
                    tools,
                    input_path,
                    fasta,
                    regions,
                    gc_target_bed,
                )
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

struct RegionStreamSource {
    tools: ExternalTools,
    bam: String,
    fasta: String,
    target_bed: Option<String>,
    pending_regions: VecDeque<String>,
}

struct StreamProcess {
    child: Child,
    stderr_capture: NamedTempFile,
    command: String,
    finished: bool,
}

impl RegionStreamSource {
    fn new(
        tools: ExternalTools,
        bam: String,
        fasta: String,
        regions: &[String],
        target_bed: Option<&str>,
    ) -> Self {
        Self {
            tools,
            bam,
            fasta,
            target_bed: target_bed.map(str::to_string),
            pending_regions: regions.iter().cloned().collect(),
        }
    }

    fn spawn_next(&mut self) -> Result<Option<CommandStream>> {
        let Some(region) = self.pending_regions.pop_front() else {
            return Ok(None);
        };
        info!(region = %region, "spawning next region mpileup stream");
        self.tools
            .spawn_samtools_mpileup_stream(
                &self.bam,
                &self.fasta,
                Some(region.as_str()),
                self.target_bed.as_deref(),
            )
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
        target_bed: Option<&str>,
    ) -> Result<Self> {
        if regions.is_empty() {
            let stream = tools.spawn_samtools_mpileup_stream(bam, fasta, None, target_bed)?;
            return Self::from_command_stream(stream);
        }

        let mut source = RegionStreamSource::new(
            tools.clone(),
            bam.to_string(),
            fasta.to_string(),
            regions,
            target_bed,
        );
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

    fn next_record_into(&mut self, record: &mut PileupRecord) -> Result<bool> {
        loop {
            self.line_buffer.clear();
            let read = self.reader.read_until(b'\n', &mut self.line_buffer)?;
            if read == 0 {
                self.ensure_process_completed()?;
                if self.advance_region_stream()? {
                    continue;
                }
                return Ok(false);
            }
            let line = trim_line_end(&self.line_buffer);
            if line.is_empty() || line.iter().all(|byte| byte.is_ascii_whitespace()) {
                continue;
            }
            parse_pileup_record_into(line, record)?;
            return Ok(true);
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

#[derive(Debug, Clone, Default)]
struct PileupRecord {
    chromosome: String,
    chromosome_group: u8,
    chromosome_rank: u32,
    position: i32,
    reference: String,
    depth: i32,
    pileup: String,
    quality: String,
}

fn parse_pileup_record_into(line: &[u8], record: &mut PileupRecord) -> Result<()> {
    let mut parts = line.splitn(6, |byte| *byte == b'\t');
    let chromosome_part = parts.next().ok_or_else(|| AppError::ParseError {
        message: format!("invalid pileup line: {}", String::from_utf8_lossy(line)),
    })?;
    let position_part = parts.next().ok_or_else(|| AppError::ParseError {
        message: format!("invalid pileup line: {}", String::from_utf8_lossy(line)),
    })?;
    let reference_part = parts.next().ok_or_else(|| AppError::ParseError {
        message: format!("invalid pileup line: {}", String::from_utf8_lossy(line)),
    })?;
    let depth_part = parts.next().ok_or_else(|| AppError::ParseError {
        message: format!("invalid pileup line: {}", String::from_utf8_lossy(line)),
    })?;
    let pileup_part = parts.next().ok_or_else(|| AppError::ParseError {
        message: format!("invalid pileup line: {}", String::from_utf8_lossy(line)),
    })?;
    let quality_part = parts.next().ok_or_else(|| AppError::ParseError {
        message: format!("invalid pileup line: {}", String::from_utf8_lossy(line)),
    })?;

    let position = parse_i32_ascii(position_part).ok_or_else(|| AppError::ParseError {
        message: format!(
            "invalid pileup position: {}",
            String::from_utf8_lossy(position_part)
        ),
    })?;
    let depth = parse_i32_ascii(depth_part).ok_or_else(|| AppError::ParseError {
        message: format!(
            "invalid pileup depth: {}",
            String::from_utf8_lossy(depth_part)
        ),
    })?;
    decode_field_into(chromosome_part, &mut record.chromosome);
    let (chromosome_group, chromosome_rank) = chromosome_sort_class(&record.chromosome);
    record.chromosome_group = chromosome_group;
    record.chromosome_rank = chromosome_rank;
    record.position = position;
    decode_field_into(reference_part, &mut record.reference);
    record.depth = depth;
    decode_field_into(pileup_part, &mut record.pileup);
    decode_field_into(quality_part, &mut record.quality);
    Ok(())
}

fn decode_field(field: &[u8]) -> String {
    match std::str::from_utf8(field) {
        Ok(text) => text.to_owned(),
        Err(_) => String::from_utf8_lossy(field).into_owned(),
    }
}

fn decode_field_into(field: &[u8], target: &mut String) {
    target.clear();
    match std::str::from_utf8(field) {
        Ok(text) => target.push_str(text),
        Err(_) => target.push_str(&String::from_utf8_lossy(field)),
    }
}

fn compare_coordinates(left: &PileupRecord, right: &PileupRecord) -> Ordering {
    match compare_chromosome_keys(
        left.chromosome_group,
        left.chromosome_rank,
        &left.chromosome,
        right.chromosome_group,
        right.chromosome_rank,
        &right.chromosome,
    ) {
        Ordering::Equal => left.position.cmp(&right.position),
        other => other,
    }
}

#[cfg(feature = "htslib-prototype")]
fn compare_chromosomes(left: &str, right: &str) -> Ordering {
    let (left_group, left_rank) = chromosome_sort_class(left);
    let (right_group, right_rank) = chromosome_sort_class(right);
    compare_chromosome_keys(left_group, left_rank, left, right_group, right_rank, right)
}

fn compare_chromosome_keys(
    left_group: u8,
    left_rank: u32,
    left: &str,
    right_group: u8,
    right_rank: u32,
    right: &str,
) -> Ordering {
    match left_group.cmp(&right_group) {
        Ordering::Equal => match left_rank.cmp(&right_rank) {
            Ordering::Equal => left.cmp(right),
            other => other,
        },
        other => other,
    }
}

fn chromosome_sort_class(chrom: &str) -> (u8, u32) {
    let raw = chrom.trim();
    let normalized = strip_chr_prefix(raw);

    if let Some(num) = parse_u32_ascii_fast(normalized.as_bytes()) {
        return (0, num);
    }

    if normalized.eq_ignore_ascii_case("X") {
        (1, 23)
    } else if normalized.eq_ignore_ascii_case("Y") {
        (1, 24)
    } else if normalized.eq_ignore_ascii_case("M") || normalized.eq_ignore_ascii_case("MT") {
        (1, 25)
    } else {
        (2, 0)
    }
}

fn strip_chr_prefix(raw: &str) -> &str {
    if raw.len() >= 3 {
        let bytes = raw.as_bytes();
        if bytes[0].eq_ignore_ascii_case(&b'c')
            && bytes[1].eq_ignore_ascii_case(&b'h')
            && bytes[2].eq_ignore_ascii_case(&b'r')
        {
            return &raw[3..];
        }
    }
    raw
}

fn parse_u32_ascii_fast(raw: &[u8]) -> Option<u32> {
    if raw.is_empty() {
        return None;
    }

    let mut value: u32 = 0;
    for &byte in raw {
        if !byte.is_ascii_digit() {
            return None;
        }
        value = value.checked_mul(10)?.checked_add(u32::from(byte - b'0'))?;
    }
    Some(value)
}

fn advance_to_target(
    current_present: &mut bool,
    current: &mut PileupRecord,
    stream: &mut PileupStream,
    target: &PileupRecord,
) -> Result<()> {
    while *current_present && compare_coordinates(current, target) == Ordering::Less {
        *current_present = stream.next_record_into(current)?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct GcInterval {
    start: i32,
    end: i32,
    gc: String,
}

struct GcLookupCursor<'a> {
    map: &'a HashMap<String, Vec<GcInterval>>,
    chromosome: Option<String>,
    intervals: Option<&'a [GcInterval]>,
    index: usize,
    last_position: Option<i32>,
}

impl<'a> GcLookupCursor<'a> {
    fn new(map: &'a HashMap<String, Vec<GcInterval>>) -> Self {
        Self {
            map,
            chromosome: None,
            intervals: None,
            index: 0,
            last_position: None,
        }
    }

    fn value_at(&mut self, chromosome: &str, position: i32) -> Option<&'a str> {
        if self.chromosome.as_deref() != Some(chromosome) {
            self.chromosome = Some(chromosome.to_string());
            self.intervals =
                find_gc_intervals_for_chromosome(self.map, chromosome).map(Vec::as_slice);
            self.index = 0;
            self.last_position = None;
        }

        let intervals = self.intervals?;
        if intervals.is_empty() {
            return None;
        }

        if let Some(previous_position) = self.last_position {
            if position < previous_position {
                let idx = intervals.partition_point(|interval| interval.start <= position);
                if idx == 0 {
                    self.index = 0;
                    self.last_position = Some(position);
                    return None;
                }
                self.index = idx - 1;
            } else {
                while self.index < intervals.len() && position >= intervals[self.index].end {
                    self.index += 1;
                }
            }
        } else {
            let idx = intervals.partition_point(|interval| interval.start <= position);
            if idx == 0 {
                self.last_position = Some(position);
                return None;
            }
            self.index = idx - 1;
        }

        self.last_position = Some(position);

        if self.index >= intervals.len() {
            return None;
        }

        let interval = &intervals[self.index];
        if position < interval.start {
            return None;
        }
        (position < interval.end).then_some(interval.gc.as_str())
    }
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
    let mut preferred_chromosomes: Vec<String> = if !args.chr.is_empty() {
        let requested = args
            .chr
            .iter()
            .map(|region| region_chromosome(region).to_string())
            .collect::<Vec<_>>();
        dedup_chromosomes_preserving_order(&requested)
    } else if !args.pileup {
        let mut chromosomes = tools.list_bam_chromosomes(&args.tumor)?;
        if let Some(normal2_path) = &args.normal2 {
            chromosomes.extend(tools.list_bam_chromosomes(normal2_path)?);
        }
        dedup_chromosomes_preserving_order(&chromosomes)
    } else {
        let mut pileup_paths = vec![args.normal.as_str(), args.tumor.as_str()];
        if let Some(normal2_path) = args.normal2.as_deref() {
            pileup_paths.push(normal2_path);
        }
        collect_pileup_chromosomes(&pileup_paths)?
    };

    preferred_chromosomes = dedup_chromosomes_preserving_order(&preferred_chromosomes);

    let mut gc_filter: HashSet<String> = HashSet::new();
    for chromosome in &preferred_chromosomes {
        extend_chromosome_alias_filter(&mut gc_filter, chromosome);
    }

    let mut gc_intervals = parse_gc_intervals_from_file(
        &args.gc,
        if gc_filter.is_empty() {
            None
        } else {
            Some(&gc_filter)
        },
    )?;
    if !preferred_chromosomes.is_empty() {
        gc_intervals =
            remap_gc_intervals_to_preferred_chromosomes(gc_intervals, &preferred_chromosomes);
    }
    let gc_interval_count = gc_intervals.values().map(Vec::len).sum::<usize>();
    info!(
        chromosomes = gc_intervals.len(),
        intervals = gc_interval_count,
        "loaded gc intervals"
    );
    Ok(gc_intervals)
}

fn collect_pileup_chromosomes(paths: &[&str]) -> Result<Vec<String>> {
    let mut seen: HashSet<String> = HashSet::new();
    let mut ordered = Vec::new();

    for path in paths {
        if path.is_empty() {
            continue;
        }

        if path.ends_with(".gz") {
            let file = File::open(path)?;
            let decoder = GzDecoder::new(file);
            collect_pileup_chromosomes_from_reader(
                BufReader::new(decoder),
                &mut seen,
                &mut ordered,
            )?;
        } else {
            let file = File::open(path)?;
            collect_pileup_chromosomes_from_reader(BufReader::new(file), &mut seen, &mut ordered)?;
        }
    }

    Ok(ordered)
}

fn collect_pileup_chromosomes_from_reader<R: BufRead>(
    mut reader: R,
    seen: &mut HashSet<String>,
    ordered: &mut Vec<String>,
) -> Result<()> {
    let mut buf = Vec::with_capacity(256);
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

        let chromosome = line.split(|byte| *byte == b'\t').next().unwrap_or(line);
        if chromosome.is_empty() {
            continue;
        }

        let chromosome_name = decode_field(chromosome);
        if seen.insert(chromosome_name.clone()) {
            ordered.push(chromosome_name);
        }
    }
    Ok(())
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
    if bytes.is_empty() {
        return None;
    }

    let mut index = 0;
    let mut is_negative = false;
    match bytes[0] {
        b'-' => {
            is_negative = true;
            index = 1;
        }
        b'+' => {
            index = 1;
        }
        _ => {}
    }

    if index == bytes.len() {
        return None;
    }

    let mut value: i32 = 0;
    while index < bytes.len() {
        let byte = bytes[index];
        if !byte.is_ascii_digit() {
            return None;
        }
        let digit = i32::from(byte - b'0');
        value = if is_negative {
            value.checked_mul(10)?.checked_sub(digit)?
        } else {
            value.checked_mul(10)?.checked_add(digit)?
        };
        index += 1;
    }

    Some(value)
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
    let Some(intervals) = find_gc_intervals_for_chromosome(gc_map, chromosome) else {
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

fn create_gc_target_bed_for_regions(
    regions: &[String],
    gc_map: &HashMap<String, Vec<GcInterval>>,
) -> Result<Option<NamedTempFile>> {
    let mut by_chromosome: HashMap<String, Vec<(i32, i32)>> = HashMap::new();

    if regions.is_empty() {
        for (chromosome, intervals) in gc_map {
            let targets = by_chromosome.entry(chromosome.clone()).or_default();
            for interval in intervals {
                if interval.end > interval.start {
                    targets.push((interval.start, interval.end));
                }
            }
        }
    } else {
        for region in regions {
            let chromosome = region_chromosome(region);
            let Some(intervals) = find_gc_intervals_for_chromosome(gc_map, chromosome) else {
                continue;
            };

            let targets = by_chromosome.entry(chromosome.to_string()).or_default();
            if let Some((region_start, region_end_inclusive)) = parse_region_bounds(region) {
                let region_end_exclusive = region_end_inclusive.saturating_add(1);
                let start_index =
                    intervals.partition_point(|interval| interval.end <= region_start);
                for interval in intervals[start_index..]
                    .iter()
                    .take_while(|interval| interval.start < region_end_exclusive)
                {
                    let overlap_start = interval.start.max(region_start);
                    let overlap_end = interval.end.min(region_end_exclusive);
                    if overlap_end > overlap_start {
                        targets.push((overlap_start, overlap_end));
                    }
                }
            } else {
                for interval in intervals {
                    if interval.end > interval.start {
                        targets.push((interval.start, interval.end));
                    }
                }
            }
        }
    }

    let mut total_intervals = 0usize;
    for intervals in by_chromosome.values_mut() {
        if intervals.is_empty() {
            continue;
        }
        intervals.sort_unstable_by(|left, right| {
            left.0.cmp(&right.0).then_with(|| left.1.cmp(&right.1))
        });

        let mut merged: Vec<(i32, i32)> = Vec::with_capacity(intervals.len());
        for (start, end) in intervals.iter().copied() {
            if let Some(last) = merged.last_mut()
                && start <= last.1
            {
                if end > last.1 {
                    last.1 = end;
                }
            } else {
                merged.push((start, end));
            }
        }
        total_intervals += merged.len();
        *intervals = merged;
    }

    if total_intervals == 0 {
        return Ok(None);
    }

    let mut chromosome_names = by_chromosome.keys().cloned().collect::<Vec<_>>();
    chromosome_names.sort();

    let mut bed_file = Builder::new()
        .prefix("bam2seqz_gc_targets_")
        .suffix(".bed")
        .tempfile_in("tmp")?;
    {
        let out = bed_file.as_file_mut();
        for chromosome in chromosome_names {
            if let Some(intervals) = by_chromosome.get(&chromosome) {
                for (start, end) in intervals {
                    let bed_start = start.saturating_sub(1).max(0);
                    let bed_end = end.saturating_sub(1).max(bed_start + 1);
                    writeln!(out, "{chromosome}\t{bed_start}\t{bed_end}")?;
                }
            }
        }
        out.flush()?;
    }

    info!(
        regions = regions.len(),
        chromosomes = by_chromosome.len(),
        intervals = total_intervals,
        path = %bed_file.path().to_string_lossy(),
        "prepared GC target bed for mpileup"
    );
    Ok(Some(bed_file))
}

fn region_chromosome(region: &str) -> &str {
    region.split(':').next().unwrap_or(region)
}

#[cfg(test)]
fn gc_value_at<'a>(
    map: &'a HashMap<String, Vec<GcInterval>>,
    chromosome: &str,
    position: i32,
) -> Option<&'a str> {
    let intervals = find_gc_intervals_for_chromosome(map, chromosome)?;
    let index = intervals.partition_point(|interval| interval.start <= position);
    if index == 0 {
        return None;
    }
    let interval = &intervals[index - 1];
    (position < interval.end).then_some(interval.gc.as_str())
}

fn validate_bam_fasta_chromosome_naming(args: &Bam2SeqzArgs, tools: &ExternalTools) -> Result<()> {
    if args.pileup {
        return Ok(());
    }

    let Some(fasta_path) = args.fasta.as_deref() else {
        return Ok(());
    };

    let fasta_chromosomes = read_fasta_index_chromosomes(fasta_path)?;
    if fasta_chromosomes.is_empty() {
        return Ok(());
    }

    let bam_targets = if args.chr.is_empty() {
        tools.list_bam_chromosomes(&args.tumor)?
    } else {
        let requested = args
            .chr
            .iter()
            .map(|region| region_chromosome(region).to_string())
            .collect::<Vec<_>>();
        dedup_chromosomes_preserving_order(&requested)
    };

    let mut alias_only_mismatches = Vec::new();
    for chromosome in bam_targets {
        if fasta_chromosomes
            .iter()
            .any(|fasta_chromosome| fasta_chromosome == &chromosome)
        {
            continue;
        }

        if fasta_chromosomes
            .iter()
            .any(|fasta_chromosome| chromosome_alias_match(&chromosome, fasta_chromosome))
        {
            alias_only_mismatches.push(chromosome);
        }
    }

    if alias_only_mismatches.is_empty() {
        return Ok(());
    }

    let preview = alias_only_mismatches
        .iter()
        .take(3)
        .cloned()
        .collect::<Vec<_>>()
        .join(", ");

    Err(AppError::InvalidValue {
        flag: "--fasta".to_string(),
        value: fasta_path.to_string(),
        reason: format!(
            "BAM and FASTA chromosome names are prefix-mismatched (e.g. {preview}). Use a FASTA/.fai with exactly matching contig names (chr-prefix differences are not safe to auto-rewrite for mpileup)."
        ),
    })
}

fn read_fasta_index_chromosomes(fasta_path: &str) -> Result<Vec<String>> {
    let index_path = format!("{fasta_path}.fai");
    if !Path::new(&index_path).exists() {
        return Ok(Vec::new());
    }

    let file = File::open(index_path)?;
    let reader = BufReader::new(file);
    let mut chromosomes = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let chromosome = line.split('\t').next().map(str::trim).unwrap_or_default();
        if !chromosome.is_empty() && !chromosomes.iter().any(|entry| entry == chromosome) {
            chromosomes.push(chromosome.to_string());
        }
    }
    Ok(chromosomes)
}

fn remap_gc_intervals_to_preferred_chromosomes(
    mut gc_intervals: HashMap<String, Vec<GcInterval>>,
    preferred_chromosomes: &[String],
) -> HashMap<String, Vec<GcInterval>> {
    if preferred_chromosomes.is_empty() {
        return gc_intervals;
    }

    let mut remapped: HashMap<String, Vec<GcInterval>> = HashMap::new();
    for (chromosome, intervals) in gc_intervals.drain() {
        let target = resolve_preferred_chromosome(&chromosome, preferred_chromosomes)
            .unwrap_or(chromosome.as_str())
            .to_string();
        remapped.entry(target).or_default().extend(intervals);
    }

    for intervals in remapped.values_mut() {
        intervals.sort_unstable_by(|left, right| {
            left.start
                .cmp(&right.start)
                .then_with(|| left.end.cmp(&right.end))
                .then_with(|| left.gc.cmp(&right.gc))
        });
        intervals.dedup_by(|left, right| {
            left.start == right.start && left.end == right.end && left.gc == right.gc
        });
    }

    remapped
}

fn extend_chromosome_alias_filter(filter: &mut HashSet<String>, chromosome: &str) {
    for alias in chromosome_aliases(chromosome) {
        let _ = filter.insert(alias);
    }
}

fn resolve_preferred_chromosome<'a>(
    chromosome: &str,
    preferred_chromosomes: &'a [String],
) -> Option<&'a str> {
    preferred_chromosomes
        .iter()
        .find(|candidate| candidate.as_str() == chromosome)
        .map(String::as_str)
        .or_else(|| {
            preferred_chromosomes
                .iter()
                .find(|candidate| chromosome_alias_match(candidate, chromosome))
                .map(String::as_str)
        })
}

fn find_gc_intervals_for_chromosome<'a>(
    gc_map: &'a HashMap<String, Vec<GcInterval>>,
    chromosome: &str,
) -> Option<&'a Vec<GcInterval>> {
    if let Some(intervals) = gc_map.get(chromosome) {
        return Some(intervals);
    }

    for alias in chromosome_aliases(chromosome) {
        if alias == chromosome {
            continue;
        }
        if let Some(intervals) = gc_map.get(alias.as_str()) {
            return Some(intervals);
        }
    }

    gc_map
        .iter()
        .find(|(candidate, _)| chromosome_alias_match(candidate, chromosome))
        .map(|(_, intervals)| intervals)
}

fn chromosome_alias_match(left: &str, right: &str) -> bool {
    canonical_chromosome_name(left) == canonical_chromosome_name(right)
}

fn canonical_chromosome_name(raw: &str) -> String {
    let normalized = strip_chr_prefix(raw.trim());
    if normalized.eq_ignore_ascii_case("MT") {
        "M".to_string()
    } else {
        normalized.to_ascii_uppercase()
    }
}

fn chromosome_aliases(chromosome: &str) -> Vec<String> {
    let raw = chromosome.trim();
    if raw.is_empty() {
        return Vec::new();
    }

    let mut aliases: Vec<String> = Vec::with_capacity(6);
    let mut push_unique = |candidate: String| {
        if !aliases
            .iter()
            .any(|existing| existing.eq_ignore_ascii_case(&candidate))
        {
            aliases.push(candidate);
        }
    };

    push_unique(raw.to_string());
    let stripped = strip_chr_prefix(raw);
    if !stripped.is_empty() {
        push_unique(stripped.to_string());
        push_unique(format!("chr{stripped}"));

        if stripped.eq_ignore_ascii_case("M") {
            push_unique("MT".to_string());
            push_unique("chrMT".to_string());
        } else if stripped.eq_ignore_ascii_case("MT") {
            push_unique("M".to_string());
            push_unique("chrM".to_string());
        }
    }

    aliases
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
    use std::collections::HashSet;

    #[test]
    fn build_plan_for_parallel_outputs() {
        let args = parse_args([
            "bam2seqz",
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
            "bam2seqz", "-n", "n.bam", "-t", "t.bam", "-gc", "gc.wig", "-F", "ref.fa", "-f",
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
            vec!["chr20:100-140".to_string(), "chr20:bad-range".to_string(),]
        );
    }

    #[test]
    fn gc_lookup_supports_chr_prefix_aliases() {
        let lines = vec![
            "variableStep chrom=chr20 span=50".to_string(),
            "100\t48".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);

        assert_eq!(super::gc_value_at(&map, "20", 120), Some("48"));
        assert_eq!(super::gc_value_at(&map, "chr20", 120), Some("48"));
    }

    #[test]
    fn filter_regions_supports_chr_prefix_aliases() {
        let lines = vec![
            "variableStep chrom=20 span=50".to_string(),
            "100\t48".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        let regions = vec!["chr20:100-140".to_string(), "chr21:100-140".to_string()];

        let filtered = super::filter_regions_with_gc(&regions, &map);
        assert_eq!(filtered, vec!["chr20:100-140".to_string()]);
    }

    #[test]
    fn auto_bins_support_chr_prefix_aliases() {
        let lines = vec![
            "variableStep chrom=chr20 span=50".to_string(),
            "100\t48".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        let chromosomes = vec!["20".to_string()];

        let bins = super::auto_bin_regions_for_chromosomes(&chromosomes, &map, 5_000_000);
        assert_eq!(bins, vec!["20:1-5000000".to_string()]);
    }

    #[test]
    fn remaps_gc_intervals_to_preferred_chromosome_style() {
        let lines = vec![
            "variableStep chrom=chr20 span=50".to_string(),
            "100\t48".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        let preferred = vec!["20".to_string()];

        let remapped = super::remap_gc_intervals_to_preferred_chromosomes(map, &preferred);
        assert!(remapped.contains_key("20"));
        assert!(!remapped.contains_key("chr20"));
        assert_eq!(super::gc_value_at(&remapped, "20", 120), Some("48"));
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
            vec!["chr20".to_string(), "chr21".to_string(), "chr1".to_string(),]
        );
    }

    #[test]
    fn adaptive_auto_bin_shrinks_for_small_coverage() {
        let covered = 60_000_000_u64;
        let bin_size = super::choose_auto_bin_size_bp(covered, 64);
        assert!(bin_size < super::AUTO_BIN_SIZE_BP);
        assert!(bin_size >= super::AUTO_BIN_MIN_SIZE_BP);
    }

    #[test]
    fn adaptive_auto_bin_keeps_default_for_large_coverage() {
        let covered = 3_000_000_000_u64;
        let bin_size = super::choose_auto_bin_size_bp(covered, 64);
        assert_eq!(bin_size, super::AUTO_BIN_SIZE_BP);
    }

    #[test]
    fn parse_i32_ascii_supports_signed_values_and_limits() {
        assert_eq!(super::parse_i32_ascii(b"0"), Some(0));
        assert_eq!(super::parse_i32_ascii(b"+42"), Some(42));
        assert_eq!(super::parse_i32_ascii(b"-42"), Some(-42));
        assert_eq!(
            super::parse_i32_ascii(i32::MAX.to_string().as_bytes()),
            Some(i32::MAX)
        );
        assert_eq!(
            super::parse_i32_ascii(i32::MIN.to_string().as_bytes()),
            Some(i32::MIN)
        );
    }

    #[test]
    fn parse_i32_ascii_rejects_invalid_or_overflow() {
        assert_eq!(super::parse_i32_ascii(b""), None);
        assert_eq!(super::parse_i32_ascii(b"-"), None);
        assert_eq!(super::parse_i32_ascii(b"+"), None);
        assert_eq!(super::parse_i32_ascii(b"12a"), None);
        assert_eq!(super::parse_i32_ascii(b"2147483648"), None);
        assert_eq!(super::parse_i32_ascii(b"-2147483649"), None);
    }

    #[test]
    fn collect_pileup_chromosomes_from_reader_dedups_preserving_order() {
        let data = b"chr20\t1\tA\t1\t.\tB\nchr20\t2\tA\t1\t.\tB\n21\t1\tC\t1\t.\tB\n";
        let mut seen = HashSet::new();
        let mut ordered = Vec::new();

        super::collect_pileup_chromosomes_from_reader(&data[..], &mut seen, &mut ordered)
            .expect("expected parse success");

        assert_eq!(ordered, vec!["chr20".to_string(), "21".to_string()]);
    }

    #[test]
    fn gc_lookup_cursor_matches_gc_value_at() {
        let lines = vec![
            "variableStep chrom=chr20 span=50".to_string(),
            "100\t48".to_string(),
            "200\t49".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        let mut cursor = super::GcLookupCursor::new(&map);

        for position in [90, 100, 120, 149, 150, 199, 200, 220, 249, 260] {
            assert_eq!(
                cursor.value_at("chr20", position),
                super::gc_value_at(&map, "chr20", position)
            );
        }
    }

    #[test]
    fn gc_lookup_cursor_handles_chromosome_switches_and_backtracking() {
        let lines = vec![
            "variableStep chrom=chr20 span=50".to_string(),
            "100\t48".to_string(),
            "200\t49".to_string(),
            "variableStep chrom=21 span=50".to_string(),
            "100\t42".to_string(),
        ];
        let map = super::parse_gc_intervals(lines);
        let mut cursor = super::GcLookupCursor::new(&map);

        assert_eq!(cursor.value_at("chr20", 220), Some("49"));
        assert_eq!(cursor.value_at("chr20", 120), Some("48"));
        assert_eq!(cursor.value_at("21", 120), Some("42"));
        assert_eq!(cursor.value_at("chr20", 120), Some("48"));
    }
}
