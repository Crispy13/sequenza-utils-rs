use bio::io::fasta;
use clap::Parser;
#[cfg(feature = "mimalloc-allocator")]
use mimalloc::MiMalloc;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::sync::Mutex;
use std::time::Instant;

#[cfg(feature = "mimalloc-allocator")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[derive(Debug, Parser)]
#[command(name = "htslib_thread_local_bench")]
struct Args {
    #[arg(long)]
    bam: String,
    #[arg(long)]
    fasta: String,
    #[arg(short = 'C', long = "region", num_args = 1.., required = true)]
    regions: Vec<String>,
    #[arg(long, default_value_t = 1)]
    threads: usize,
    #[arg(long, default_value_t = 1)]
    bam_threads: usize,
}

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    start0: u64,
    end0: u64,
}

#[derive(Debug, Default, Clone, Copy)]
struct RegionStats {
    loci: u64,
    aligned_bases: u64,
    quality_sum: u64,
    checksum: u64,
}

impl RegionStats {
    fn merge(self, other: Self) -> Self {
        Self {
            loci: self.loci + other.loci,
            aligned_bases: self.aligned_bases + other.aligned_bases,
            quality_sum: self.quality_sum + other.quality_sum,
            checksum: self.checksum.wrapping_add(other.checksum),
        }
    }
}

struct WorkerContext {
    bam_reader: bam::IndexedReader,
    fasta_reader: fasta::IndexedReader<std::fs::File>,
}

impl WorkerContext {
    fn new(bam_path: &str, fasta_path: &str, bam_threads: usize) -> Result<Self, String> {
        let mut bam_reader = bam::IndexedReader::from_path(bam_path)
            .map_err(|err| format!("failed to open BAM index reader: {err}"))?;
        if bam_threads > 0 {
            bam_reader
                .set_threads(bam_threads)
                .map_err(|err| format!("failed to set BAM worker threads: {err}"))?;
        }

        let fasta_reader = fasta::IndexedReader::from_file(&fasta_path)
            .map_err(|err| format!("failed to open indexed FASTA reader: {err}"))?;

        Ok(Self {
            bam_reader,
            fasta_reader,
        })
    }

    fn process_region(&mut self, region: &Region) -> Result<RegionStats, String> {
        let tid = self
            .bam_reader
            .header()
            .tid(region.chrom.as_bytes())
            .ok_or_else(|| format!("chromosome not found in BAM header: {}", region.chrom))?;

        self.bam_reader
            .fetch((tid, region.start0, region.end0))
            .map_err(|err| {
                format!(
                    "failed BAM fetch for {}:{}-{}: {err}",
                    region.chrom,
                    region.start0 + 1,
                    region.end0
                )
            })?;

        self.fasta_reader
            .fetch(&region.chrom, region.start0, region.end0)
            .map_err(|err| {
                format!(
                    "failed FASTA fetch for {}:{}-{}: {err}",
                    region.chrom,
                    region.start0 + 1,
                    region.end0
                )
            })?;

        let mut reference_bases = Vec::new();
        self.fasta_reader
            .read(&mut reference_bases)
            .map_err(|err| {
                format!(
                    "failed FASTA read for {}:{}-{}: {err}",
                    region.chrom,
                    region.start0 + 1,
                    region.end0
                )
            })?;

        let mut stats = RegionStats::default();
        for pileup_result in self.bam_reader.pileup() {
            let pileup = pileup_result.map_err(|err| format!("failed pileup iteration: {err}"))?;
            let position0 = u64::from(pileup.pos());
            if position0 < region.start0 || position0 >= region.end0 {
                continue;
            }

            let Some(reference_base) = reference_bases.get((position0 - region.start0) as usize)
            else {
                continue;
            };
            stats.loci += 1;

            for alignment in pileup.alignments() {
                if alignment.is_del() || alignment.is_refskip() {
                    continue;
                }

                let Some(query_position) = alignment.qpos() else {
                    continue;
                };

                let record = alignment.record();
                if record.is_unmapped()
                    || record.is_duplicate()
                    || record.is_secondary()
                    || record.is_supplementary()
                    || record.is_quality_check_failed()
                {
                    continue;
                }

                let sequence = record.seq().as_bytes();
                let Some(base) = sequence.get(query_position).copied() else {
                    continue;
                };
                let Some(qual) = record.qual().get(query_position).copied() else {
                    continue;
                };

                stats.aligned_bases += 1;
                stats.quality_sum += u64::from(qual);
                stats.checksum = stats.checksum.wrapping_mul(1_099_511_628_211).wrapping_add(
                    u64::from(base) + u64::from(*reference_base) + u64::from(qual) + position0,
                );
            }
        }

        Ok(stats)
    }
}

fn parse_region(region: &str) -> Result<Region, String> {
    let (chrom, bounds) = region
        .split_once(':')
        .ok_or_else(|| format!("invalid region (missing ':'): {region}"))?;
    let (start_raw, end_raw) = bounds
        .split_once('-')
        .ok_or_else(|| format!("invalid region (missing '-'): {region}"))?;

    let start1 = start_raw
        .replace(',', "")
        .parse::<u64>()
        .map_err(|_| format!("invalid start coordinate: {start_raw}"))?;
    let end1 = end_raw
        .replace(',', "")
        .parse::<u64>()
        .map_err(|_| format!("invalid end coordinate: {end_raw}"))?;

    if start1 == 0 || end1 < start1 {
        return Err(format!("invalid region bounds: {region}"));
    }

    Ok(Region {
        chrom: chrom.to_string(),
        start0: start1 - 1,
        end0: end1,
    })
}

fn main() -> Result<(), String> {
    let args = Args::parse();

    if args.threads == 0 {
        return Err("--threads must be >= 1".to_string());
    }

    let regions = args
        .regions
        .iter()
        .map(|region| parse_region(region))
        .collect::<Result<Vec<_>, _>>()?;

    let started = Instant::now();

    let per_region = Mutex::new(vec![RegionStats::default(); regions.len()]);
    let pool = ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .map_err(|err| format!("failed to create rayon thread pool: {err}"))?;

    pool.install(|| {
        regions.par_iter().enumerate().try_for_each_init(
            || WorkerContext::new(&args.bam, &args.fasta, args.bam_threads),
            |context_result, (index, region)| -> Result<(), String> {
                let context = match context_result {
                    Ok(context) => context,
                    Err(err) => return Err(err.clone()),
                };
                let stats = context.process_region(region)?;
                let mut guard = per_region
                    .lock()
                    .map_err(|_| "failed to lock region-stats mutex".to_string())?;
                guard[index] = stats;
                Ok(())
            },
        )
    })?;

    let combined = per_region
        .into_inner()
        .map_err(|_| "failed to take region stats from mutex".to_string())?
        .into_iter()
        .fold(RegionStats::default(), RegionStats::merge);

    let elapsed_seconds = started.elapsed().as_secs_f64();

    println!("impl\trust_htslib_thread_local");
    println!("threads\t{}", args.threads);
    println!("bam_threads\t{}", args.bam_threads);
    println!("regions\t{}", regions.len());
    println!("elapsed_sec\t{elapsed_seconds:.4}");
    println!("loci\t{}", combined.loci);
    println!("aligned_bases\t{}", combined.aligned_bases);
    println!("quality_sum\t{}", combined.quality_sum);
    println!("checksum\t{}", combined.checksum);

    Ok(())
}
