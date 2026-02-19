#!/usr/bin/env bash
set -euo pipefail

cd /home/eck/workspace/bam2seqz_rs

export PATH="/home/eck/software/miniconda3/envs/rust_build_env/bin:$PATH"

NORMAL="testdata/NA12878.chrom20.synthetic.seed20260218.normal.chr.bam"
TUMOR="testdata/NA12878.chrom20.synthetic.seed20260218.tumor.chr.bam"
GC="testdata/hg19_gc50.wig.gz"
FASTA="testdata/hg19.fa"
SAMTOOLS="/home/eck/software/miniconda3/envs/rust_build_env/bin/samtools"
TABIX="/home/eck/software/miniconda3/envs/rust_build_env/bin/tabix"
RUNS=3

OUT_DIR="tmp/backend_bench_p8"
mkdir -p "$OUT_DIR"

conda run -n rust_build_env cargo build --profile debug-release --features htslib-prototype --bin bam2seqz

for i in $(seq 1 "$RUNS"); do
  /usr/bin/time -v -o "$OUT_DIR/samtools_run_${i}.time.txt" \
    target/debug-release/bam2seqz \
    -n "$NORMAL" -t "$TUMOR" -gc "$GC" -F "$FASTA" \
    --parallel 8 -S "$SAMTOOLS" -T "$TABIX" \
    -o "$OUT_DIR/samtools_run_${i}.seqz"

  /usr/bin/time -v -o "$OUT_DIR/rust_htslib_run_${i}.time.txt" \
    target/debug-release/bam2seqz \
    -n "$NORMAL" -t "$TUMOR" -gc "$GC" -F "$FASTA" \
    --parallel 8 --bam-backend rust-htslib \
    -o "$OUT_DIR/rust_htslib_run_${i}.seqz"
done

python - <<'PY'
from pathlib import Path
import re
import statistics

base = Path("tmp/backend_bench_p8")
runs = 3


def parse_elapsed_seconds(text: str) -> float:
    match = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)", text)
    if not match:
        return 0.0
    value = match.group(1).strip()
    parts = value.split(":")
    if len(parts) == 2:
        minute, sec = parts
        return int(minute) * 60 + float(sec)
    if len(parts) == 3:
        hour, minute, sec = parts
        return int(hour) * 3600 + int(minute) * 60 + float(sec)
    return float(value)


def parse_cpu(text: str) -> int:
    match = re.search(r"Percent of CPU this job got:\s*(\d+)%", text)
    return int(match.group(1)) if match else 0


def parse_rss(text: str) -> int:
    match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    return int(match.group(1)) if match else 0


def percentile(values: list[float], p: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = (len(ordered) - 1) * p
    lower = int(rank)
    upper = min(lower + 1, len(ordered) - 1)
    frac = rank - lower
    return ordered[lower] * (1.0 - frac) + ordered[upper] * frac


detail_rows = []
stats = {}
for mode in ["samtools", "rust_htslib"]:
    elapsed_values = []
    rss_values = []
    cpu_values = []
    for run_idx in range(1, runs + 1):
        timing_path = base / f"{mode}_run_{run_idx}.time.txt"
        text = timing_path.read_text()
        elapsed = parse_elapsed_seconds(text)
        cpu = parse_cpu(text)
        rss = parse_rss(text)
        detail_rows.append((mode, run_idx, elapsed, cpu, rss))
        elapsed_values.append(elapsed)
        cpu_values.append(cpu)
        rss_values.append(rss)

    stats[mode] = {
        "median_elapsed": statistics.median(elapsed_values),
        "p95_elapsed": percentile(elapsed_values, 0.95),
        "median_cpu": statistics.median(cpu_values),
        "max_rss": max(rss_values),
    }

with (base / "details.tsv").open("w") as out:
    out.write("mode\trun\telapsed_sec\tcpu_pct\tmax_rss_kb\n")
    for mode, run_idx, elapsed, cpu, rss in detail_rows:
        out.write(f"{mode}\t{run_idx}\t{elapsed:.2f}\t{cpu}\t{rss}\n")

with (base / "summary.tsv").open("w") as out:
    out.write("mode\tmedian_elapsed_sec\tp95_elapsed_sec\tmedian_cpu_pct\tmax_rss_kb\n")
    for mode in ["samtools", "rust_htslib"]:
        row = stats[mode]
        out.write(
            f"{mode}\t{row['median_elapsed']:.2f}\t{row['p95_elapsed']:.2f}\t{int(row['median_cpu'])}\t{row['max_rss']}\n"
        )

ratio = 0.0
if stats["rust_htslib"]["median_elapsed"] > 0:
    ratio = stats["samtools"]["median_elapsed"] / stats["rust_htslib"]["median_elapsed"]

with (base / "comparison.txt").open("w") as out:
    out.write(
        "speedup_ratio_samtools_over_rust_htslib="
        f"{ratio:.3f}\n"
    )

print(base / "summary.tsv")
print(base / "details.tsv")
print(base / "comparison.txt")
PY

echo "[done] benchmark artifacts under $OUT_DIR"