#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-parity}"
MEM_MB="${MAX_MEM_MB:-6144}"
BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
TEST_THREADS="${TEST_THREADS:-1}"
PROFILE="${TEST_PROFILE:-debug-release}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

mkdir -p tmp

case "${MODE}" in
    parity)
        TEST_CMD=(cargo test --profile "${PROFILE}" --test parity_golden -- --nocapture --test-threads="${TEST_THREADS}")
        LOG_FILE="tmp/limited_parity_$(date +%Y%m%d_%H%M%S).log"
        ;;
    lib)
        TEST_CMD=(cargo test --profile "${PROFILE}" --lib -- --test-threads="${TEST_THREADS}")
        LOG_FILE="tmp/limited_lib_$(date +%Y%m%d_%H%M%S).log"
        ;;
    all)
        TEST_CMD=(cargo test --profile "${PROFILE}" -q -- --test-threads="${TEST_THREADS}")
        LOG_FILE="tmp/limited_all_$(date +%Y%m%d_%H%M%S).log"
        ;;
    *)
        echo "Usage: $0 [parity|lib|all]"
        echo "Environment overrides: MAX_MEM_MB=6144 CARGO_BUILD_JOBS=1 TEST_THREADS=1 TEST_PROFILE=debug-release"
        exit 2
        ;;
esac

echo "[limited-tests] mode=${MODE} profile=${PROFILE} mem_mb=${MEM_MB} build_jobs=${BUILD_JOBS} test_threads=${TEST_THREADS}" | tee "${LOG_FILE}"

if command -v prlimit >/dev/null 2>&1; then
    set -o pipefail
    /usr/bin/time -v prlimit --as="$((MEM_MB * 1024 * 1024))" -- env CARGO_BUILD_JOBS="${BUILD_JOBS}" "${TEST_CMD[@]}" 2>&1 | tee -a "${LOG_FILE}"
else
    set -o pipefail
    ulimit -v "$((MEM_MB * 1024))"
    /usr/bin/time -v env CARGO_BUILD_JOBS="${BUILD_JOBS}" "${TEST_CMD[@]}" 2>&1 | tee -a "${LOG_FILE}"
fi

echo "[limited-tests] log=${LOG_FILE}"