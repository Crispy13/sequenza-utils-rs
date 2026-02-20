#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

usage() {
  cat <<'USAGE'
Usage: scripts/benchmark_release_snapshot.sh --version <semver> [--tag <git-tag>]

Runs release benchmark snapshots, archives artifacts in tmp/benchmarks/<version>,
and appends benchmark metadata into copilot-desk/benchmark-history.tsv (if missing).
USAGE
}

VERSION=""
TAG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --version)
      VERSION="${2:-}"
      shift 2
      ;;
    --tag)
      TAG="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$VERSION" ]]; then
  echo "--version is required" >&2
  usage
  exit 1
fi

if [[ -z "$TAG" ]]; then
  TAG="v${VERSION}"
fi

export PATH="/home/eck/software/miniconda3/envs/rust_build_env/bin:$PATH"
export PYTHONPATH="/home/eck/workspace/bam2seqz_rs/sequenza-utils:${PYTHONPATH:-}"

OUT_DIR="tmp/benchmarks/${VERSION}"
mkdir -p "$OUT_DIR"

scripts/benchmark_python_rust_parallel8.sh

if [[ "${RUN_BACKEND_BENCH:-1}" == "1" ]]; then
  scripts/benchmark_backend_binning_p8.sh
fi

copy_if_exists() {
  local src="$1"
  local dst="$2"
  if [[ -e "$src" ]]; then
    cp -f "$src" "$dst"
  fi
}

copy_if_exists "tmp/parity_resource_metrics_parallel8.tsv" "$OUT_DIR/parity_resource_metrics_parallel8.tsv"
copy_if_exists "tmp/parity_resource_metrics.tsv" "$OUT_DIR/parity_resource_metrics.tsv"

if [[ -d "tmp/python_parity_timings_parallel8" ]]; then
  rm -rf "$OUT_DIR/python_parity_timings_parallel8"
  cp -r "tmp/python_parity_timings_parallel8" "$OUT_DIR/python_parity_timings_parallel8"
fi

if [[ -d "tmp/rust_parity_timings_parallel8" ]]; then
  rm -rf "$OUT_DIR/rust_parity_timings_parallel8"
  cp -r "tmp/rust_parity_timings_parallel8" "$OUT_DIR/rust_parity_timings_parallel8"
fi

copy_if_exists "tmp/backend_bench_p8/summary.tsv" "$OUT_DIR/backend_bench_p8_summary.tsv"
copy_if_exists "tmp/backend_bench_p8/details.tsv" "$OUT_DIR/backend_bench_p8_details.tsv"
copy_if_exists "tmp/backend_bench_p8/comparison.txt" "$OUT_DIR/backend_bench_p8_comparison.txt"

if git rev-parse --verify "$TAG" >/dev/null 2>&1; then
  COMMIT_SHA="$(git rev-parse "${TAG}^{}")"
  RELEASE_DATE="$(git show -s --format=%cI "${TAG}^{}")"
else
  COMMIT_SHA="$(git rev-parse HEAD)"
  RELEASE_DATE="$(git show -s --format=%cI HEAD)"
fi

GENERATED_AT="$(date -Iseconds)"
SAMTOOLS_VERSION="$(samtools --version 2>/dev/null | head -n 1 || true)"
RUST_VERSION="$(cargo --version 2>/dev/null || true)"
PYTHON_VERSION="$(/home/eck/software/miniconda3/envs/rust_build_env/bin/python --version 2>&1 || true)"

cat > "$OUT_DIR/metadata.tsv" <<EOF
key\tvalue
version\t$VERSION
tag\t$TAG
commit\t$COMMIT_SHA
release_date\t$RELEASE_DATE
generated_at\t$GENERATED_AT
samtools\t$SAMTOOLS_VERSION
rust\t$RUST_VERSION
python\t$PYTHON_VERSION
EOF

HISTORY_FILE="copilot-desk/benchmark-history.tsv"
if [[ -f "$HISTORY_FILE" ]] && ! grep -q "^${VERSION}[[:space:]]" "$HISTORY_FILE"; then
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$VERSION" "$TAG" "$COMMIT_SHA" "$RELEASE_DATE" \
    "$OUT_DIR/parity_resource_metrics_parallel8.tsv" \
    "$OUT_DIR/backend_bench_p8_summary.tsv" \
    "release snapshot" >> "$HISTORY_FILE"
fi

echo "[done] release benchmark snapshot: $OUT_DIR"
