#!/usr/bin/env bash
set -euo pipefail

cd /home/eck/workspace/bam2seqz_rs

export PATH="/home/eck/software/miniconda3/envs/rust_build_env/bin:$PATH"
export PYTHONPATH="/home/eck/workspace/bam2seqz_rs/sequenza-utils:${PYTHONPATH:-}"

RUN_STANDARD_SCENARIOS="${RUN_STANDARD_SCENARIOS:-1}"
RUN_LARGE_BAM_ALLCHR="${RUN_LARGE_BAM_ALLCHR:-0}"
RUN_LARGE_BAM_ALLCHR_PY="${RUN_LARGE_BAM_ALLCHR_PY:-0}"

mkdir -p tmp/python_parity_outputs_parallel8 tmp/rust_parity_outputs_parallel8 \
  tmp/python_parity_timings_parallel8 tmp/rust_parity_timings_parallel8
rm -f tmp/python_parity_timings_parallel8/*.time.txt \
  tmp/rust_parity_timings_parallel8/*.time.txt

conda run -n rust_build_env cargo build --profile debug-release --bin bam2seqz

N="testdata/NA12878.chrom20.synthetic.seed20260218.normal.chr.bam"
T="testdata/NA12878.chrom20.synthetic.seed20260218.tumor.chr.bam"
NP="testdata/NA12878.chrom20.synthetic.seed20260218.normal.region59900_64000.pileup.gz"
TP="testdata/NA12878.chrom20.synthetic.seed20260218.tumor.region59900_64000.pileup.gz"
GC="testdata/hg19_gc50.wig.gz"
FA="testdata/hg19.fa"
R1="chr20:59900-60412"; R2="chr20:60413-60925"; R3="chr20:60926-61438"; R4="chr20:61439-61951"
R5="chr20:61952-62464"; R6="chr20:62465-62977"; R7="chr20:62978-63490"; R8="chr20:63491-64000"

PY=/home/eck/software/miniconda3/envs/rust_build_env/bin/python

N_LARGE="${N_LARGE:-tmp/fixtures/NA12878.mapped.synthetic.seed20260221.normal.bam}"
T_LARGE="${T_LARGE:-tmp/fixtures/NA12878.mapped.synthetic.seed20260221.tumor.bam}"
FA_LARGE="${FA_LARGE:-$FA}"
GC_LARGE="${GC_LARGE:-$GC}"

if [[ "$RUN_STANDARD_SCENARIOS" == "1" ]]; then
  /usr/bin/time -v -o tmp/python_parity_timings_parallel8/bam_default.time.txt "$PY" -m sequenza.commands bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" -o tmp/python_parity_outputs_parallel8/bam_default.seqz.gz
  /usr/bin/time -v -o tmp/python_parity_timings_parallel8/pileup_default.time.txt "$PY" -m sequenza.commands bam2seqz --pileup -n "$NP" -t "$TP" -gc "$GC" -o tmp/python_parity_outputs_parallel8/pileup_default.seqz.gz
  /usr/bin/time -v -o tmp/python_parity_timings_parallel8/region_12.time.txt "$PY" -m sequenza.commands bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" -C chr20:62001-64000 -o tmp/python_parity_outputs_parallel8/region_12.seqz.gz
  /usr/bin/time -v -o tmp/python_parity_timings_parallel8/parallel_regions_p8.time.txt "$PY" -m sequenza.commands bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" -C "$R1" "$R2" "$R3" "$R4" "$R5" "$R6" "$R7" "$R8" --parallel 8 -o tmp/python_parity_outputs_parallel8/parallel_p8.seqz.gz
  /usr/bin/time -v -o tmp/python_parity_timings_parallel8/normal2.time.txt "$PY" -m sequenza.commands bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" --normal2 "$N" -o tmp/python_parity_outputs_parallel8/normal2.seqz.gz
  /usr/bin/time -v -o tmp/python_parity_timings_parallel8/bad_qformat_error.time.txt bash -lc 'PY=/home/eck/software/miniconda3/envs/rust_build_env/bin/python; export PYTHONPATH="/home/eck/workspace/bam2seqz_rs/sequenza-utils:$PYTHONPATH"; "$PY" -m sequenza.commands bam2seqz -n "'"$N"'" -t "'"$T"'" -gc "'"$GC"'" -F "'"$FA"'" -f invalid -o tmp/python_parity_outputs_parallel8/bad_q.seqz.gz >/dev/null 2>&1; test $? -ne 0'
  /usr/bin/time -v -o tmp/python_parity_timings_parallel8/parallel_stdout_error_p8.time.txt bash -lc 'PY=/home/eck/software/miniconda3/envs/rust_build_env/bin/python; export PYTHONPATH="/home/eck/workspace/bam2seqz_rs/sequenza-utils:$PYTHONPATH"; "$PY" -m sequenza.commands bam2seqz -n "'"$N"'" -t "'"$T"'" -gc "'"$GC"'" -F "'"$FA"'" -C "'"$R1"'" "'"$R2"'" "'"$R3"'" "'"$R4"'" "'"$R5"'" "'"$R6"'" "'"$R7"'" "'"$R8"'" --parallel 8 -o - >/dev/null 2>&1; test $? -ne 0'
fi

if [[ "$RUN_STANDARD_SCENARIOS" == "1" ]]; then
  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/bam_default.time.txt target/debug-release/bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" -o tmp/rust_parity_outputs_parallel8/bam_default.seqz.gz
  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/pileup_default.time.txt target/debug-release/bam2seqz --pileup -n "$NP" -t "$TP" -gc "$GC" -o tmp/rust_parity_outputs_parallel8/pileup_default.seqz.gz
  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/region_12.time.txt target/debug-release/bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" -C chr20:62001-64000 -o tmp/rust_parity_outputs_parallel8/region_12.seqz.gz
  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/parallel_regions_p8.time.txt target/debug-release/bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" -C "$R1" "$R2" "$R3" "$R4" "$R5" "$R6" "$R7" "$R8" --parallel 8 -o tmp/rust_parity_outputs_parallel8/parallel_p8.seqz.gz
  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/auto_bin_single_output_p8.time.txt target/debug-release/bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" --parallel 8 -o tmp/rust_parity_outputs_parallel8/auto_bin_single_output_p8.seqz.gz
  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/normal2.time.txt target/debug-release/bam2seqz -n "$N" -t "$T" -gc "$GC" -F "$FA" --normal2 "$N" -o tmp/rust_parity_outputs_parallel8/normal2.seqz.gz
  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/bad_qformat_error.time.txt bash -lc 'target/debug-release/bam2seqz -n "'"$N"'" -t "'"$T"'" -gc "'"$GC"'" -F "'"$FA"'" -f invalid -o tmp/rust_parity_outputs_parallel8/bad_q.seqz.gz >/dev/null 2>&1; test $? -ne 0'
  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/parallel_stdout_error_p8.time.txt bash -lc 'target/debug-release/bam2seqz -n "'"$N"'" -t "'"$T"'" -gc "'"$GC"'" -F "'"$FA"'" -C "'"$R1"'" "'"$R2"'" "'"$R3"'" "'"$R4"'" "'"$R5"'" "'"$R6"'" "'"$R7"'" "'"$R8"'" --parallel 8 -o - >/dev/null 2>&1; test $? -ne 0'
fi

if [[ "$RUN_LARGE_BAM_ALLCHR" == "1" ]]; then
  if [[ ! -f "$N_LARGE" || ! -f "$T_LARGE" ]]; then
    echo "large BAM fixture files not found: N_LARGE=$N_LARGE T_LARGE=$T_LARGE" >&2
    exit 1
  fi

  mapfile -t AVAILABLE_CHROMS < <(samtools idxstats "$N_LARGE" | awk '$1 != "*" { print $1 }')

  HAS_CHR_PREFIX=1
  for chrom in chr{1..22} chrX chrY; do
    if ! printf '%s\n' "${AVAILABLE_CHROMS[@]}" | grep -Fxq "$chrom"; then
      HAS_CHR_PREFIX=0
      break
    fi
  done

  if [[ "$HAS_CHR_PREFIX" == "1" ]]; then
    LARGE_CHROMS=(chr{1..22} chrX chrY)
  else
    HAS_NO_PREFIX=1
    for chrom in {1..22} X Y; do
      if ! printf '%s\n' "${AVAILABLE_CHROMS[@]}" | grep -Fxq "$chrom"; then
        HAS_NO_PREFIX=0
        break
      fi
    done

    if [[ "$HAS_NO_PREFIX" == "1" ]]; then
      LARGE_CHROMS=({1..22} X Y)
    else
      echo "expected canonical 24 chromosomes (chr1-chr22, chrX, chrY) or (1-22, X, Y)" >&2
      exit 1
    fi
  fi

  if [[ "$RUN_LARGE_BAM_ALLCHR_PY" == "1" ]]; then
    /usr/bin/time -v -o tmp/python_parity_timings_parallel8/large_bam_chr1_22_xy_default.time.txt "$PY" -m sequenza.commands bam2seqz -n "$N_LARGE" -t "$T_LARGE" -gc "$GC_LARGE" -F "$FA_LARGE" -C "${LARGE_CHROMS[@]}" --parallel 8 -o tmp/python_parity_outputs_parallel8/large_bam_chr1_22_xy_default.seqz.gz
  fi

  /usr/bin/time -v -o tmp/rust_parity_timings_parallel8/large_bam_chr1_22_xy_default.time.txt target/debug-release/bam2seqz -n "$N_LARGE" -t "$T_LARGE" -gc "$GC_LARGE" -F "$FA_LARGE" -C "${LARGE_CHROMS[@]}" --parallel 8 -o tmp/rust_parity_outputs_parallel8/large_bam_chr1_22_xy_default.seqz.gz
fi

conda run -n rust_build_env /home/eck/software/miniconda3/envs/rust_build_env/bin/python scripts/generate_parity_metrics.py

echo "benchmark_python_rust_parallel8_done"
