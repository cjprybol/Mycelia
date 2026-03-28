#!/usr/bin/env bash
# Mycelia HPC CI driver
#
# Runs on an HPC machine (e.g., via sbatch) against a checked-out Mycelia repo.
# Phases:
#   1. Run Pkg.test() with coverage.
#   2. Summarize coverage to LCOV for Codecov.
#   3. Optionally run extended tutorials + benchmarks.
#   4. Emit hpc-ci/hpc-results.json for downstream ingestion.
#   5. Optionally upload coverage to Codecov with an hpc-extended flag.
#
# Environment knobs:
#   HPC_RESULTS_DIR   Where to write hpc-results.json and artifacts (default: <repo>/hpc-ci).
#   CODECOV_TOKEN     If set, coverage uploads to Codecov.
#   HPC_CODECOV_FLAG  Codecov flag (default: hpc-extended).
#   HPC_CI_COMMIT     Commit SHA (default: git rev-parse HEAD).
#   HPC_CI_BRANCH     Branch name (default: git rev-parse --abbrev-ref HEAD).
#   HPC_CI_CHECKOUT   If set to 1, fetch/checkout HPC_CI_COMMIT/HPC_CI_BRANCH first.
#   HPC_CLUSTER_NAME  Optional logical cluster name for metadata.
#
# CLI options:
#   --tests-only       Run tests + coverage only.
#   --benchmarks-only  Run extended tutorials/benchmarks only.
#   --with-benchmarks  Run tests + extended tutorials/benchmarks.
#   --no-codecov       Skip Codecov upload even if CODECOV_TOKEN is set.

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --tests-only       Run tests + coverage only.
  --benchmarks-only  Run extended tutorials/benchmarks only.
  --with-benchmarks  Run tests + extended tutorials/benchmarks.
  --upload-only      Upload existing coverage only (no tests/benchmarks).
  --no-codecov       Skip Codecov upload even if CODECOV_TOKEN is set.
  --help             Show this help message.

This script assumes it is located under ci/hpc/ inside the Mycelia repo.
EOF
}

# Resolve repo root relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$REPO_ROOT"

# Optional checkout to a specific commit/branch for reproducible HPC runs
HPC_CI_CHECKOUT="${HPC_CI_CHECKOUT:-0}"
HPC_CI_COMMIT="${HPC_CI_COMMIT:-}"
HPC_CI_BRANCH="${HPC_CI_BRANCH:-}"

if [[ "${HPC_CI_CHECKOUT}" -eq 1 ]]; then
  if [[ -n "$(git status --porcelain)" ]]; then
    echo "Error: working tree is not clean; refusing to checkout." >&2
    exit 1
  fi

  if [[ -n "${HPC_CI_BRANCH}" ]]; then
    git fetch origin "${HPC_CI_BRANCH}"
  else
    git fetch origin
  fi

  if [[ -n "${HPC_CI_COMMIT}" ]]; then
    git checkout "${HPC_CI_COMMIT}"
  elif [[ -n "${HPC_CI_BRANCH}" ]]; then
    git checkout "${HPC_CI_BRANCH}"
  fi
fi

# Default behaviour: run tests only and upload to Codecov if configured
RUN_TESTS=1
RUN_BENCHMARKS=0
UPLOAD_CODECOV=1
RUN_UPLOAD_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tests-only)
      RUN_BENCHMARKS=0
      shift
      ;;
    --benchmarks-only)
      RUN_TESTS=0
      RUN_BENCHMARKS=1
      shift
      ;;
    --with-benchmarks)
      RUN_BENCHMARKS=1
      shift
      ;;
    --upload-only)
      RUN_TESTS=0
      RUN_BENCHMARKS=0
      RUN_UPLOAD_ONLY=1
      shift
      ;;
    --no-codecov)
      UPLOAD_CODECOV=0
      shift
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# Basic metadata
export HPC_CI_COMMIT="${HPC_CI_COMMIT:-$(git rev-parse HEAD)}"
export HPC_CI_BRANCH="${HPC_CI_BRANCH:-$(git rev-parse --abbrev-ref HEAD)}"
export HPC_CI_REPO="${HPC_CI_REPO:-$(basename "$(git rev-parse --show-toplevel)")}"

export HPC_CLUSTER_NAME="${HPC_CLUSTER_NAME:-${SLURM_CLUSTER_NAME:-unknown}}"
export HPC_JOB_ID="${HPC_JOB_ID:-${SLURM_JOB_ID:-}}"
export HPC_NODE="${HPC_NODE:-${SLURMD_NODENAME:-$(hostname)}}"

# Julia version
export HPC_JULIA_VERSION="$(julia --version 2>/dev/null | awk '{print $3}')"

# Results + artifact locations
export HPC_RESULTS_DIR="${HPC_RESULTS_DIR:-${REPO_ROOT}/hpc-ci}"
export HPC_CODECOV_FLAG="${HPC_CODECOV_FLAG:-hpc-extended}"

COVERAGE_DIR="${HPC_RESULTS_DIR}/coverage"
LOG_DIR="${HPC_RESULTS_DIR}/logs"

mkdir -p "${HPC_RESULTS_DIR}" "${COVERAGE_DIR}" "${LOG_DIR}"

echo "=== Mycelia HPC CI ==="
echo "Repo       : ${HPC_CI_REPO}"
echo "Commit     : ${HPC_CI_COMMIT}"
echo "Branch     : ${HPC_CI_BRANCH}"
echo "Cluster    : ${HPC_CLUSTER_NAME}"
echo "Job ID     : ${HPC_JOB_ID}"
echo "Node       : ${HPC_NODE}"
echo "Julia      : ${HPC_JULIA_VERSION}"
echo "Results dir: ${HPC_RESULTS_DIR}"
echo

# Helper vars for JSON summary
export HPC_TEST_STATUS="skipped"
export HPC_TEST_DURATION=0
export HPC_BENCH_STATUS="skipped"
export HPC_BENCH_DURATION=0
export HPC_COVERAGE_SUMMARY="${COVERAGE_DIR}/summary.json"

TEST_EXIT=0
BENCH_EXIT=0

###############################################################################
# Phase 1: Pkg.test() with coverage
###############################################################################
if [[ "${RUN_TESTS}" -eq 1 ]]; then
  echo ">>> Phase 1: Running Pkg.test() with coverage (MYCELIA_RUN_ALL=true)"
  TEST_START=$(date +%s)

  set +e
  # LD_LIBRARY_PATH is left as-is; clear it before invoking if cluster needs it.
  JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-auto}" \
  MYCELIA_RUN_ALL=true \
  # Coverage needs non-precompiled code; Pkg.test inherits this flag.
  julia --project=. --compiled-modules=no \
        -e 'import Pkg; Pkg.test(coverage=true)' \
        2>&1 | tee "${LOG_DIR}/tests.log"
  TEST_EXIT=${PIPESTATUS[0]}
  set -e

  TEST_END=$(date +%s)
  export HPC_TEST_DURATION=$((TEST_END - TEST_START))

  if [[ "${TEST_EXIT}" -eq 0 ]]; then
    export HPC_TEST_STATUS="passed"
  else
    export HPC_TEST_STATUS="failed"
  fi

  echo ">>> Phase 1: Pkg.test() finished with status: ${HPC_TEST_STATUS}"

  echo ">>> Phase 1b: Converting coverage to LCOV and JSON summary"
  export HPC_COVERAGE_DIR="${COVERAGE_DIR}"

  julia --project=. <<'EOF'
import Pkg

temp_env = mktempdir()
Pkg.activate(temp_env)
Pkg.add(["Coverage", "JSON"])

using Coverage
import JSON

results_dir = get(ENV, "HPC_RESULTS_DIR", "hpc-ci")
cov_dir     = get(ENV, "HPC_COVERAGE_DIR", joinpath(results_dir, "coverage"))
mkpath(cov_dir)

# Process coverage from src/ and write LCOV
coverage = process_folder("src")
LCOV.writefile(joinpath(cov_dir, "lcov.info"), coverage)

covered, total = get_summary(coverage)
pct = total == 0 ? 0.0 : 100 * covered / total

summary = Dict(
    "total_lines"      => total,
    "covered_lines"    => covered,
    "coverage_percent" => round(pct; digits=2),
)

summary_path = get(ENV, "HPC_COVERAGE_SUMMARY", joinpath(cov_dir, "summary.json"))
mkpath(dirname(summary_path))
open(summary_path, "w") do io
    JSON.print(io, summary, 4)
end

println("Coverage summary:")
println(JSON.json(summary))
EOF

fi

###############################################################################
# Phase 2: Extended tutorials + benchmarks (optional, no coverage requirement)
###############################################################################
if [[ "${RUN_BENCHMARKS}" -eq 1 ]]; then
  echo
  echo ">>> Phase 2: Running extended tutorials + benchmarks"
  BENCH_START=$(date +%s)

  set +e
  julia --project=. run_extended_tests.jl tutorials \
        2>&1 | tee "${LOG_DIR}/tutorials.log"
  TUTORIAL_EXIT=${PIPESTATUS[0]}

  julia --project=. run_extended_tests.jl benchmarks \
        2>&1 | tee "${LOG_DIR}/benchmarks.log"
  BENCHMARK_EXIT=${PIPESTATUS[0]}

  # Optionally build HTML reports for local inspection
  julia --project=. run_extended_tests.jl reports \
        2>&1 | tee "${LOG_DIR}/reports.log"
  REPORTS_EXIT=${PIPESTATUS[0]}
  set -e

  BENCH_END=$(date +%s)
  export HPC_BENCH_DURATION=$((BENCH_END - BENCH_START))

  if [[ "${TUTORIAL_EXIT}" -eq 0 && "${BENCHMARK_EXIT}" -eq 0 ]]; then
    export HPC_BENCH_STATUS="passed"
    BENCH_EXIT=0
  else
    export HPC_BENCH_STATUS="failed"
    BENCH_EXIT=1
  fi

  echo ">>> Phase 2: Extended tests finished with status: ${HPC_BENCH_STATUS}"
fi

###############################################################################
# Phase 3: Optional Codecov upload
###############################################################################
if [[ ( "${RUN_TESTS}" -eq 1 || "${RUN_UPLOAD_ONLY}" -eq 1 ) && "${UPLOAD_CODECOV}" -eq 1 ]]; then
  if [[ -n "${CODECOV_TOKEN:-}" && -f "${COVERAGE_DIR}/lcov.info" ]]; then
    CODECOV_LOG="${LOG_DIR}/codecov.log"
    : > "${CODECOV_LOG}"
    {
      echo "Codecov upload context:"
      echo "  Coverage file: ${COVERAGE_DIR}/lcov.info"
      echo "  Coverage size: $(wc -c < "${COVERAGE_DIR}/lcov.info") bytes"
      echo "  Coverage lines: $(wc -l < "${COVERAGE_DIR}/lcov.info")"
      echo "  Commit: ${HPC_CI_COMMIT}"
      echo "  Branch: ${HPC_CI_BRANCH}"
      echo "  Flag: ${HPC_CODECOV_FLAG}"
    } >> "${CODECOV_LOG}"

    CODECOV_BRANCH_ARGS=()
    if [[ -n "${HPC_CI_BRANCH}" && "${HPC_CI_BRANCH}" != "HEAD" ]]; then
      CODECOV_BRANCH_ARGS=(-B "${HPC_CI_BRANCH}")
    fi

    if command -v codecov >/dev/null 2>&1; then
      echo
      echo ">>> Phase 3: Uploading coverage to Codecov (flag=${HPC_CODECOV_FLAG})"
      set +e
      codecov \
        -t "${CODECOV_TOKEN}" \
        -F "${HPC_CODECOV_FLAG}" \
        -f "${COVERAGE_DIR}/lcov.info" \
        -C "${HPC_CI_COMMIT}" \
        "${CODECOV_BRANCH_ARGS[@]}" \
        2>&1 | tee -a "${CODECOV_LOG}"
      CODECOV_EXIT=${PIPESTATUS[0]}
      set -e
      if [[ "${CODECOV_EXIT}" -ne 0 ]]; then
        echo "Warning: Codecov upload failed (non-fatal), exit code: ${CODECOV_EXIT}" >&2
      fi
    else
      echo "Info: 'codecov' CLI not found; attempting Coverage.jl upload."
      set +e
      julia --project=. <<'EOF' 2>&1 | tee -a "${CODECOV_LOG}"
import Pkg

temp_env = mktempdir()
Pkg.activate(temp_env)
Pkg.add("Coverage")

import Coverage

coverage = Coverage.process_folder("src")
flag = get(ENV, "HPC_CODECOV_FLAG", "")

if isempty(flag)
    Coverage.Codecov.submit(coverage)
else
    try
        Coverage.Codecov.submit(coverage; flags=flag)
    catch err
        if err isa MethodError
            Coverage.Codecov.submit(coverage)
        else
            rethrow()
        end
    end
end
EOF
      CODECOV_EXIT=${PIPESTATUS[0]}
      set -e
      if [[ "${CODECOV_EXIT}" -ne 0 ]]; then
        echo "Warning: Coverage.jl Codecov upload failed (non-fatal), exit code: ${CODECOV_EXIT}" >&2
      fi
    fi
  else
    echo "Info: CODECOV_TOKEN not set or coverage file missing; skipping Codecov upload." >&2
  fi
fi

###############################################################################
# Phase 4: Emit hpc-results.json summary
###############################################################################
echo
echo ">>> Phase 4: Writing hpc-results.json summary"

julia --project=. <<'EOF'
import JSON, Dates

# Utility to parse integers from ENV
function parse_int(name::String, default::Int)
    val = get(ENV, name, "")
    try
        return isempty(val) ? default : parse(Int, val)
    catch
        return default
    end
end

results_dir = get(ENV, "HPC_RESULTS_DIR", "hpc-ci")
cov_summary_path = get(ENV, "HPC_COVERAGE_SUMMARY", "")

coverage = nothing
if !isempty(cov_summary_path) && isfile(cov_summary_path)
    try
        coverage = JSON.parsefile(cov_summary_path)
    catch err
        coverage = nothing
        println(stderr, "Warning: failed to parse coverage summary at ", cov_summary_path, ": ", err)
    end
end

hpc = Dict(
    "cluster"       => get(ENV, "HPC_CLUSTER_NAME", "unknown"),
    "job_id"        => get(ENV, "HPC_JOB_ID", ""),
    "node"          => get(ENV, "HPC_NODE", ""),
    "julia_version" => get(ENV, "HPC_JULIA_VERSION", ""),
)

cov_obj = coverage === nothing ? nothing :
    Dict(
        "flag"             => get(ENV, "HPC_CODECOV_FLAG", "hpc-extended"),
        "coverage_percent" => coverage["coverage_percent"],
        "total_lines"      => coverage["total_lines"],
        "covered_lines"    => coverage["covered_lines"],
        # Relative to HPC_RESULTS_DIR
        "report_file"      => joinpath("coverage", "lcov.info"),
        "summary_file"     => joinpath("coverage", "summary.json"),
    )

tests = Dict(
    "status"           => get(ENV, "HPC_TEST_STATUS", "skipped"),
    "duration_seconds" => parse_int("HPC_TEST_DURATION", 0),
    "log_file"         => joinpath("logs", "tests.log"),
)

benchmarks = Dict(
    "status"           => get(ENV, "HPC_BENCH_STATUS", "skipped"),
    "duration_seconds" => parse_int("HPC_BENCH_DURATION", 0),
    "tutorials_log"    => joinpath("logs", "tutorials.log"),
    "benchmarks_log"   => joinpath("logs", "benchmarks.log"),
)

result = Dict(
    "schema_version" => 1,
    "commit"         => get(ENV, "HPC_CI_COMMIT", "unknown"),
    "branch"         => get(ENV, "HPC_CI_BRANCH", "unknown"),
    "generated_at"   => Dates.format(Dates.now(Dates.UTC),
                                     Dates.dateformat"yyyy-mm-ddTHH:MM:SSZ"),
    "hpc"            => hpc,
    "coverage"       => cov_obj,
    "tests"          => tests,
    "benchmarks"     => benchmarks,
)

mkpath(results_dir)
out_path = joinpath(results_dir, "hpc-results.json")
open(out_path, "w") do io
    JSON.print(io, result, 4)
end

println("Wrote HPC results summary to: ", out_path)
EOF

echo
echo "=== Mycelia HPC CI complete ==="

# Surface failure if any phase failed
if [[ "${TEST_EXIT}" -ne 0 || "${BENCH_EXIT}" -ne 0 ]]; then
  exit 1
fi

exit 0
