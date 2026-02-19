#!/usr/bin/env bash
set -euo pipefail

cd /home/eck/workspace/bam2seqz_rs

export PATH="/home/eck/software/miniconda3/envs/rust_build_env/bin:$PATH"

BAM="testdata/NA12878.chrom20.synthetic.seed20260218.tumor.chr.bam"
FASTA="testdata/hg19.fa"
SAMTOOLS="/home/eck/software/miniconda3/envs/rust_build_env/bin/samtools"

REGIONS=(
  "chr20:59900-62000"
  "chr20:62001-64000"
  "chr20:64001-66000"
  "chr20:66001-68000"
  "chr20:68001-70000"
  "chr20:70001-72000"
  "chr20:72001-74000"
  "chr20:74001-76000"
)

REGION_ARGS=()
for region in "${REGIONS[@]}"; do
  REGION_ARGS+=("-C" "$region")
done

mkdir -p tmp/htslib_bench

conda run -n rust_build_env cargo build --profile debug-release --features htslib-prototype --bin htslib_thread_local_bench

/usr/bin/time -v -o tmp/htslib_bench/htslib_t1.time.txt \
  target/debug-release/htslib_thread_local_bench \
  --bam "$BAM" --fasta "$FASTA" --threads 1 --bam-threads 1 \
  "${REGION_ARGS[@]}" \
  > tmp/htslib_bench/htslib_t1.out

/usr/bin/time -v -o tmp/htslib_bench/htslib_t8.time.txt \
  target/debug-release/htslib_thread_local_bench \
  --bam "$BAM" --fasta "$FASTA" --threads 8 --bam-threads 1 \
  "${REGION_ARGS[@]}" \
  > tmp/htslib_bench/htslib_t8.out

/usr/bin/time -v -o tmp/htslib_bench/samtools_mpileup_per_region.time.txt \
  bash -lc 'set -euo pipefail; SAMTOOLS="$1"; FASTA="$2"; BAM="$3"; shift 3; for region in "$@"; do "$SAMTOOLS" mpileup -f "$FASTA" -q 20 -Q 20 -r "$region" "$BAM" >/dev/null; done' _ \
  "$SAMTOOLS" "$FASTA" "$BAM" "${REGIONS[@]}"

/usr/bin/time -v -o tmp/htslib_bench/samtools_view_pipe_mpileup_per_region.time.txt \
  bash -lc 'set -euo pipefail; SAMTOOLS="$1"; FASTA="$2"; BAM="$3"; shift 3; for region in "$@"; do "$SAMTOOLS" view -u "$BAM" "$region" | "$SAMTOOLS" mpileup -f "$FASTA" -q 20 -Q 20 - >/dev/null; done' _ \
  "$SAMTOOLS" "$FASTA" "$BAM" "${REGIONS[@]}"

python - <<'PY'
from pathlib import Path
import re

base = Path('tmp/htslib_bench')
items = [
    ('htslib_t1', 'htslib_t1.time.txt'),
    ('htslib_t8', 'htslib_t8.time.txt'),
    ('samtools_mpileup_per_region', 'samtools_mpileup_per_region.time.txt'),
    ('samtools_view_pipe_mpileup_per_region', 'samtools_view_pipe_mpileup_per_region.time.txt'),
]


def parse_elapsed_seconds(text: str) -> float:
    match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)', text)
    if not match:
        return 0.0
    value = match.group(1).strip()
    parts = value.split(':')
    if len(parts) == 2:
        minute, sec = parts
        return int(minute) * 60 + float(sec)
    if len(parts) == 3:
        hour, minute, sec = parts
        return int(hour) * 3600 + int(minute) * 60 + float(sec)
    return float(value)


def parse_cpu(text: str) -> int:
    match = re.search(r'Percent of CPU this job got:\s*(\d+)%', text)
    return int(match.group(1)) if match else 0


def parse_rss(text: str) -> int:
    match = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
    return int(match.group(1)) if match else 0

rows = []
for name, filename in items:
    data = (base / filename).read_text()
    rows.append((name, parse_elapsed_seconds(data), parse_cpu(data), parse_rss(data)))

summary = base / 'summary.tsv'
with summary.open('w') as out:
    out.write('mode\telapsed_sec\tcpu_pct\tmax_rss_kb\n')
    for mode, elapsed, cpu, rss in rows:
        out.write(f'{mode}\t{elapsed:.2f}\t{cpu}\t{rss}\n')

print(summary)
PY

echo "[done] htslib benchmark summary: tmp/htslib_bench/summary.tsv"