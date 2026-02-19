# sequenza-utils-rs / bam2seqz (v0.1)

Rust implementation of `sequenza-utils` `bam2seqz` with parity-first behavior against Python `sequenza-utils` and an experimental BAM backend switch.

## Version 0.1 scope

- Default BAM backend: external `samtools mpileup` compatibility path.
- Experimental BAM backend: `--bam-backend rust-htslib` (feature-gated by `htslib-prototype`).
- Parallel support:
  - region-parallel multi-output,
  - ordered region-parallel single-output (`--parallel-single-output`),
  - auto-binning for BAM mode when using `--parallel` without explicit ranged regions.

## Build and runtime settings used for benchmarks

- Rust profile: `debug-release`.
- Feature flag for experimental backend benchmarking: `--features htslib-prototype`.
- Environment: conda env `rust_build_env`.
- Reference and fixtures:
  - FASTA: `testdata/hg19.fa`
  - GC wig: `testdata/hg19_gc50.wig.gz`
  - Normal BAM: `testdata/NA12878.chrom20.synthetic.seed20260218.normal.chr.bam`
  - Tumor BAM: `testdata/NA12878.chrom20.synthetic.seed20260218.tumor.chr.bam`

## FAQ: what is `bam_default` test?

`bam_default` means a BAM-input run without explicit region arguments (`-C`) and without forcing parallel options. In other words, it is the default whole-input behavior:

- command shape: `bam2seqz -n <normal.bam> -t <tumor.bam> -gc <gc.wig.gz> -F <ref.fa> -o <out>`
- backend defaults to `samtools` unless `--bam-backend rust-htslib` is set
- output includes all loci emitted by the active backend and downstream seqz filters

This is the baseline scenario used to compare end-to-end behavior and cost against region-specific and parallel scenarios.

## Benchmark artifacts

- Main benchmark (Python vs Rust):
  - artifact: `tmp/parity_resource_metrics_parallel8.tsv`
  - scenarios: `bam_default`, `pileup_default`, `region_12`, `parallel_regions_p8`, `normal2`, and error-semantic checks
- Extra backend benchmark (rust-htslib vs mpileup/samtools):
  - focused p8 auto-bin script: `scripts/benchmark_backend_binning_p8.sh`
  - focused outputs: `tmp/backend_bench_p8/{summary.tsv,details.tsv,comparison.txt}`
  - full matrix script: `scripts/benchmark_bam_backend_matrix.sh`
  - full matrix outputs: `tmp/backend_bam_matrix/{summary.tsv,details.tsv,comparison.tsv,correctness.tsv}`

## Benchmark results

### 1) Main benchmark: Python vs Rust (parallel-8 benchmark batch)

Source: `tmp/parity_resource_metrics_parallel8.tsv`

| scenario | Python elapsed | Rust elapsed | Notes |
|---|---:|---:|---|
| bam_default | 2:52.66 | 0:47.15 | Default BAM-input baseline |
| pileup_default | 0:28.71 | 0:08.23 | Pileup-input mode |
| region_12 | 0:27.79 | 0:01.72 | Single explicit region |
| parallel_regions_p8 | 0:44.58 | 0:01.85 | 8 explicit regions, `--parallel 8` |
| normal2 | 3:29.10 | 1:10.51 | BAM with `--normal2` |

This is the primary benchmark section for v0.1 release messaging.

### 2) Extra benchmark: rust-htslib vs mpileup/samtools backend comparison

#### 2.1 Auto-binning backend comparison (`--parallel 8`, 3 runs)

Source: `tmp/backend_bench_p8/summary.tsv`

| backend | median elapsed (s) | p95 elapsed (s) | median CPU (%) | max RSS (kB) |
|---|---:|---:|---:|---:|
| samtools | 13.67 | 14.03 | 839 | 82,168 |
| rust-htslib | 12.62 | 12.90 | 528 | 140,660 |

Observed ratio from `tmp/backend_bench_p8/comparison.txt`:

- `speedup_ratio_samtools_over_rust_htslib=1.083` (this means rust-htslib is faster on this scenario, with higher RSS).

#### 2.2 Full BAM scenario matrix (latest run)

Source: `tmp/backend_bam_matrix/comparison.tsv`

| scenario | samtools median (s) | rust-htslib median (s) | samtools/rust ratio | samtools RSS (kB) | rust RSS (kB) | rust/sam RSS ratio |
|---|---:|---:|---:|---:|---:|---:|
| single_default | 48.31 | 50.88 | 0.949 | 80,668 | 145,696 | 1.806 |
| single_region | 1.68 | 1.64 | 1.024 | 82,296 | 83,236 | 1.011 |
| parallel_multi_output_p8 | 1.85 | 1.67 | 1.108 | 81,176 | 87,224 | 1.075 |
| parallel_single_output_p8 | 1.78 | 2.50 | 0.712 | 82,136 | 88,920 | 1.083 |
| auto_bin_implicit_p8 | 13.15 | 12.22 | 1.076 | 82,624 | 137,816 | 1.668 |
| auto_bin_contig_p8 | 14.32 | 12.02 | 1.191 | 81,644 | 137,892 | 1.689 |
| normal2_default | 68.83 | 70.24 | 0.980 | 81,276 | 145,912 | 1.795 |

## Interpretation notes

- Main release benchmark remains Rust vs Python parity scenarios (`tmp/parity_resource_metrics_parallel8.tsv`).
- Backend comparison (`samtools` vs `rust-htslib`) is an extra section for implementation evaluation.
- Runtime is mixed by scenario: rust-htslib is faster in some parallel/auto-bin scenarios and slower in others.
- Memory is consistently higher with rust-htslib in current implementation.
- Semantics checks in `tmp/backend_bam_matrix/correctness.tsv` show:
  - ordered merge consistency passes,
  - auto-binning invariants pass per backend,
  - cross-backend output hashes still differ in several scenarios.

Because of those correctness deltas, `samtools` remains the default backend in v0.1.
