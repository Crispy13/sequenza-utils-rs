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
RUNS="${RUNS:-1}"

OUT_DIR="tmp/backend_bam_matrix"
DETAIL_DIR="$OUT_DIR/details"
OUTPUT_DIR="$OUT_DIR/outputs"
mkdir -p "$DETAIL_DIR" "$OUTPUT_DIR"

REGIONS=(
    "chr20:322500-323500"
    "chr20:5792200-5793200"
    "chr20:20229300-20230300"
    "chr20:29436200-29437200"
    "chr20:29505400-29506400"
    "chr20:29628300-29629300"
    "chr20:36381800-36382800"
    "chr20:48254500-48255500"
)

REGION_ARGS=()
for region in "${REGIONS[@]}"; do
  REGION_ARGS+=("-C" "$region")
done

printf '%s\n' "${REGIONS[@]}" > "$OUT_DIR/regions.txt"

cat > "$OUT_DIR/runs.tsv" <<'EOF'
scenario	kind	backend	rep	out_base
EOF

split_ext_like_python() {
  local path="$1"
  case "$path" in
    *.seqz.gz)
      printf '%s\t%s\n' "${path%.seqz.gz}" ".seqz.gz"
      ;;
    *.txt.gz)
      printf '%s\t%s\n' "${path%.txt.gz}" ".txt.gz"
      ;;
    *.bz2)
      printf '%s\t%s\n' "${path%.bz2}" ".bz2"
      ;;
    *.gz)
      printf '%s\t%s\n' "${path%.gz}" ".gz"
      ;;
    *.txt)
      printf '%s\t%s\n' "${path%.txt}" ".txt"
      ;;
    *.*)
      printf '%s\t.%s\n' "${path%.*}" "${path##*.}"
      ;;
    *)
      printf '%s\t%s\n' "$path" ""
      ;;
  esac
}

output_for_region() {
  local output="$1"
  local region="$2"
  local sanitized="${region//:/_}"
  sanitized="${sanitized//-/_}"
  local split
  split="$(split_ext_like_python "$output")"
  local prefix="${split%%$'\t'*}"
  local extension="${split#*$'\t'}"
  printf '%s_%s%s\n' "$prefix" "$sanitized" "$extension"
}

cleanup_single_output() {
  local out_base="$1"
  rm -f "$out_base" "$out_base.tbi"
}

cleanup_multi_output() {
  local out_base="$1"
  cleanup_single_output "$out_base"
  for region in "${REGIONS[@]}"; do
    local region_path
    region_path="$(output_for_region "$out_base" "$region")"
    rm -f "$region_path" "$region_path.tbi"
  done
}

run_case() {
  local scenario="$1"
  local kind="$2"
  local backend="$3"
  local rep="$4"
  local out_base="$5"
  shift 5

  if [[ "$kind" == "multi" ]]; then
    cleanup_multi_output "$out_base"
  else
    cleanup_single_output "$out_base"
  fi

  local time_file="$DETAIL_DIR/${scenario}.${backend}.r${rep}.time.txt"
  local stdout_file="$DETAIL_DIR/${scenario}.${backend}.r${rep}.stdout.txt"
  local stderr_file="$DETAIL_DIR/${scenario}.${backend}.r${rep}.stderr.txt"

  /usr/bin/time -v -o "$time_file" \
    target/debug-release/bam2seqz_rs \
    -n "$NORMAL" -t "$TUMOR" -gc "$GC" -F "$FASTA" \
    -S "$SAMTOOLS" -T "$TABIX" --bam-backend "$backend" \
    "$@" -o "$out_base" \
    > "$stdout_file" 2> "$stderr_file"

  printf '%s\t%s\t%s\t%s\t%s\n' "$scenario" "$kind" "$backend" "$rep" "$out_base" >> "$OUT_DIR/runs.tsv"
}

echo "[build] compiling bam2seqz_rs with htslib-prototype feature"
conda run -n rust_build_env cargo build --profile debug-release --features htslib-prototype --bin bam2seqz_rs

for rep in $(seq 1 "$RUNS"); do
  for backend in samtools rust-htslib; do
    run_case "single_default" "single" "$backend" "$rep" "$OUTPUT_DIR/single_default.${backend}.r${rep}.seqz"
    run_case "single_region" "single" "$backend" "$rep" "$OUTPUT_DIR/single_region.${backend}.r${rep}.seqz" -C "chr20:322500-323500"
    run_case "parallel_multi_output_p8" "multi" "$backend" "$rep" "$OUTPUT_DIR/parallel_multi_output_p8.${backend}.r${rep}.seqz" "${REGION_ARGS[@]}" --parallel 8
    run_case "parallel_single_output_p8" "single" "$backend" "$rep" "$OUTPUT_DIR/parallel_single_output_p8.${backend}.r${rep}.seqz" "${REGION_ARGS[@]}" --parallel 8 --parallel-single-output
    run_case "auto_bin_implicit_p8" "single" "$backend" "$rep" "$OUTPUT_DIR/auto_bin_implicit_p8.${backend}.r${rep}.seqz" --parallel 8
    run_case "auto_bin_contig_p8" "single" "$backend" "$rep" "$OUTPUT_DIR/auto_bin_contig_p8.${backend}.r${rep}.seqz" -C "chr20" --parallel 8
    run_case "normal2_default" "single" "$backend" "$rep" "$OUTPUT_DIR/normal2_default.${backend}.r${rep}.seqz" --normal2 "$NORMAL"
  done
done

python - <<'PY'
from __future__ import annotations

import csv
import hashlib
import re
import statistics
from pathlib import Path

base = Path("tmp/backend_bam_matrix")
detail_dir = base / "details"
runs_tsv = base / "runs.tsv"
regions = [line.strip() for line in (base / "regions.txt").read_text().splitlines() if line.strip()]


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


def parse_user_seconds(text: str) -> float:
    match = re.search(r"User time \(seconds\):\s*([0-9.]+)", text)
    return float(match.group(1)) if match else 0.0


def parse_sys_seconds(text: str) -> float:
    match = re.search(r"System time \(seconds\):\s*([0-9.]+)", text)
    return float(match.group(1)) if match else 0.0


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


def split_ext_like_python(path: str) -> tuple[str, str]:
    for compound in (".seqz.gz", ".txt.gz", ".bz2", ".gz", ".txt"):
        if path.endswith(compound):
            return path[: -len(compound)], compound
    if "." in path:
        prefix, suffix = path.rsplit(".", 1)
        return prefix, f".{suffix}"
    return path, ""


def output_for_region(output: str, region: str) -> str:
    prefix, extension = split_ext_like_python(output)
    sanitized = region.replace(":", "_").replace("-", "_")
    return f"{prefix}_{sanitized}{extension}"


def sha256_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def compare_files(lhs: Path, rhs: Path) -> tuple[str, str]:
    if not lhs.exists() or not rhs.exists():
        return "missing", f"lhs_exists={lhs.exists()} rhs_exists={rhs.exists()}"
    lhs_hash = sha256_file(lhs)
    rhs_hash = sha256_file(rhs)
    return ("pass", lhs_hash) if lhs_hash == rhs_hash else ("fail", f"{lhs_hash} != {rhs_hash}")


rows: list[dict[str, str]] = []
with runs_tsv.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        rows.append(row)

details_rows: list[dict[str, str | int | float]] = []
for row in rows:
    scenario = row["scenario"]
    backend = row["backend"]
    rep = row["rep"]
    time_file = detail_dir / f"{scenario}.{backend}.r{rep}.time.txt"
    timing_text = time_file.read_text()
    details_rows.append(
        {
            "scenario": scenario,
            "kind": row["kind"],
            "backend": backend,
            "rep": int(rep),
            "out_base": row["out_base"],
            "elapsed_sec": parse_elapsed_seconds(timing_text),
            "cpu_pct": parse_cpu(timing_text),
            "max_rss_kb": parse_rss(timing_text),
            "user_sec": parse_user_seconds(timing_text),
            "sys_sec": parse_sys_seconds(timing_text),
        }
    )

details_path = base / "details.tsv"
with details_path.open("w") as out:
    out.write(
        "scenario\tkind\tbackend\trep\telapsed_sec\tcpu_pct\tmax_rss_kb\tuser_sec\tsys_sec\tout_base\n"
    )
    for item in details_rows:
        out.write(
            f"{item['scenario']}\t{item['kind']}\t{item['backend']}\t{item['rep']}\t"
            f"{item['elapsed_sec']:.2f}\t{item['cpu_pct']}\t{item['max_rss_kb']}\t"
            f"{item['user_sec']:.2f}\t{item['sys_sec']:.2f}\t{item['out_base']}\n"
        )

grouped: dict[tuple[str, str], list[dict[str, str | int | float]]] = {}
for item in details_rows:
    grouped.setdefault((str(item["scenario"]), str(item["backend"])), []).append(item)

summary_rows: list[dict[str, str | int | float]] = []
for (scenario, backend), items in sorted(grouped.items()):
    elapsed = [float(entry["elapsed_sec"]) for entry in items]
    cpu = [float(entry["cpu_pct"]) for entry in items]
    rss = [float(entry["max_rss_kb"]) for entry in items]
    summary_rows.append(
        {
            "scenario": scenario,
            "backend": backend,
            "runs": len(items),
            "median_elapsed_sec": statistics.median(elapsed),
            "p95_elapsed_sec": percentile(elapsed, 0.95),
            "median_cpu_pct": statistics.median(cpu),
            "max_rss_kb": max(rss) if rss else 0.0,
        }
    )

summary_path = base / "summary.tsv"
with summary_path.open("w") as out:
    out.write(
        "scenario\tbackend\truns\tmedian_elapsed_sec\tp95_elapsed_sec\tmedian_cpu_pct\tmax_rss_kb\n"
    )
    for item in summary_rows:
        out.write(
            f"{item['scenario']}\t{item['backend']}\t{item['runs']}\t"
            f"{item['median_elapsed_sec']:.2f}\t{item['p95_elapsed_sec']:.2f}\t"
            f"{int(item['median_cpu_pct'])}\t{int(item['max_rss_kb'])}\n"
        )

summary_index = {(str(item["scenario"]), str(item["backend"])): item for item in summary_rows}
comparison_path = base / "comparison.tsv"
with comparison_path.open("w") as out:
    out.write(
        "scenario\tsamtools_median_elapsed_sec\trust_htslib_median_elapsed_sec\t"
        "speedup_ratio_samtools_over_rust_htslib\tsamtools_max_rss_kb\trust_htslib_max_rss_kb\trss_ratio_rust_over_samtools\n"
    )
    for scenario in sorted({str(item["scenario"]) for item in summary_rows}):
        sam = summary_index.get((scenario, "samtools"))
        rust = summary_index.get((scenario, "rust-htslib"))
        if not sam or not rust:
            continue
        sam_elapsed = float(sam["median_elapsed_sec"])
        rust_elapsed = float(rust["median_elapsed_sec"])
        speedup = sam_elapsed / rust_elapsed if rust_elapsed > 0 else 0.0
        sam_rss = float(sam["max_rss_kb"])
        rust_rss = float(rust["max_rss_kb"])
        rss_ratio = rust_rss / sam_rss if sam_rss > 0 else 0.0
        out.write(
            f"{scenario}\t{sam_elapsed:.2f}\t{rust_elapsed:.2f}\t{speedup:.3f}\t"
            f"{int(sam_rss)}\t{int(rust_rss)}\t{rss_ratio:.3f}\n"
        )

rep1_index = {
    (row["scenario"], row["backend"]): Path(row["out_base"])
    for row in rows
    if row["rep"] == "1"
}

correctness_rows: list[tuple[str, str, str, str]] = []


def add_correctness(scenario: str, check_name: str, status: str, detail: str) -> None:
    correctness_rows.append((scenario, check_name, status, detail))


single_scenarios = [
    "single_default",
    "single_region",
    "parallel_single_output_p8",
    "auto_bin_implicit_p8",
    "auto_bin_contig_p8",
    "normal2_default",
]

for scenario in single_scenarios:
    sam = rep1_index.get((scenario, "samtools"))
    rust = rep1_index.get((scenario, "rust-htslib"))
    if not sam or not rust:
        add_correctness(scenario, "samtools_vs_rust_file_hash", "missing", "missing rep1 output")
        continue
    status, detail = compare_files(sam, rust)
    add_correctness(scenario, "samtools_vs_rust_file_hash", status, detail)

multi_scenario = "parallel_multi_output_p8"
sam_multi = rep1_index.get((multi_scenario, "samtools"))
rust_multi = rep1_index.get((multi_scenario, "rust-htslib"))
if sam_multi and rust_multi:
    for region in regions:
        sam_region = Path(output_for_region(str(sam_multi), region))
        rust_region = Path(output_for_region(str(rust_multi), region))
        status, detail = compare_files(sam_region, rust_region)
        add_correctness(multi_scenario, f"samtools_vs_rust_region_hash:{region}", status, detail)
else:
    add_correctness(multi_scenario, "samtools_vs_rust_region_hash", "missing", "missing rep1 output")

for backend in ("samtools", "rust-htslib"):
    merged = rep1_index.get(("parallel_single_output_p8", backend))
    split_base = rep1_index.get(("parallel_multi_output_p8", backend))
    if not merged or not split_base or not merged.exists():
        add_correctness(
            "parallel_single_output_p8",
            f"ordered_merge_consistency:{backend}",
            "missing",
            "missing merged or split outputs",
        )
        continue

    merged_text = merged.read_text()
    concatenated_parts: list[str] = []
    for index, region in enumerate(regions):
        region_path = Path(output_for_region(str(split_base), region))
        if not region_path.exists():
            concatenated_parts = []
            break
        region_text = region_path.read_text()
        lines = region_text.splitlines(keepends=True)
        if index == 0:
            concatenated_parts.extend(lines)
        else:
            concatenated_parts.extend(lines[1:])

    if not concatenated_parts:
        add_correctness(
            "parallel_single_output_p8",
            f"ordered_merge_consistency:{backend}",
            "missing",
            "missing one or more split region files",
        )
        continue

    expected_text = "".join(concatenated_parts)
    status = "pass" if merged_text == expected_text else "fail"
    add_correctness(
        "parallel_single_output_p8",
        f"ordered_merge_consistency:{backend}",
        status,
        f"merged_lines={len(merged_text.splitlines())} expected_lines={len(expected_text.splitlines())}",
    )

for backend in ("samtools", "rust-htslib"):
    baseline = rep1_index.get(("single_default", backend))
    implicit = rep1_index.get(("auto_bin_implicit_p8", backend))
    contig = rep1_index.get(("auto_bin_contig_p8", backend))

    if baseline and implicit:
        status, detail = compare_files(baseline, implicit)
        add_correctness("auto_bin_implicit_p8", f"matches_single_default:{backend}", status, detail)
    else:
        add_correctness(
            "auto_bin_implicit_p8",
            f"matches_single_default:{backend}",
            "missing",
            "missing baseline or implicit output",
        )

    if baseline and contig:
        status, detail = compare_files(baseline, contig)
        add_correctness("auto_bin_contig_p8", f"matches_single_default:{backend}", status, detail)
    else:
        add_correctness(
            "auto_bin_contig_p8",
            f"matches_single_default:{backend}",
            "missing",
            "missing baseline or contig output",
        )

correctness_path = base / "correctness.tsv"
with correctness_path.open("w") as out:
    out.write("scenario\tcheck\tstatus\tdetail\n")
    for scenario, check, status, detail in correctness_rows:
        out.write(f"{scenario}\t{check}\t{status}\t{detail}\n")

print(details_path)
print(summary_path)
print(comparison_path)
print(correctness_path)
PY

echo "[done] backend benchmark matrix artifacts under $OUT_DIR"