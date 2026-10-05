#!/bin/bash
# Self-test for run_ont_read_identity_lrc.sbatch: runs the wrapper against a
# stub `julia` in a scratch HOME and repo, so its result checks can be exercised
# without a cluster, conda, or a real Julia.
#
# Each case isolates ONE check: the stub writes exactly the artifacts that case
# needs, so a case can only fail for the reason it names. Not covered: the
# setup guards that exit before julia runs (julia not runnable, wrong cwd,
# mkdir failure). Exits non-zero if any case returns the wrong verdict.
#
# Usage:
#   /bin/bash benchmarking/run_ont_read_identity_lrc.selftest.sh [wrapper]
# Runs under /bin/bash 3.2 (macOS) and bash 4+/5 (Linux).
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WRAPPER="${1:-${SCRIPT_DIR}/run_ont_read_identity_lrc.sbatch}"
[[ -f "${WRAPPER}" ]] || { echo "wrapper not found: ${WRAPPER}" >&2; exit 2; }
WRAPPER="$(cd "$(dirname "${WRAPPER}")" && pwd)/$(basename "${WRAPPER}")"

SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/ont-identity-selftest.XXXXXX")"
trap 'rm -rf "${SCRATCH}"' EXIT

# --- stub julia ---------------------------------------------------------------
# --version answers; -e is the preflight (fails when STUB_PREFLIGHT=fail);
# otherwise it is the main run, which writes into --output-dir according to
# STUB_MODE and exits with STUB_EXIT.
mkdir -p "${SCRATCH}/home/.juliaup/bin" "${SCRATCH}/repo/benchmarking"
: > "${SCRATCH}/repo/benchmarking/ont_read_identity.jl"
cat > "${SCRATCH}/home/.juliaup/bin/julia" <<'STUB'
#!/bin/bash
for a in "$@"; do [[ "$a" == "--version" ]] && { echo "julia version 1.10.10"; exit 0; }; done
for a in "$@"; do
  [[ "$a" == "-e" ]] && { [[ "${STUB_PREFLIGHT:-ok}" == "fail" ]] && exit 4; exit 0; }
done
out=""; prev=""
for a in "$@"; do [[ "$prev" == "--output-dir" ]] && out="$a"; prev="$a"; done
header() { printf 'organism\tbadread_version\tsource\n'; }
row() { printf 'Lambda\t%s\t%s\n' "$1" "$2"; }
summary_good() { { header; row "Badread v0.4.1" blast; row "Badread v0.4.1" gap; } > "$out/read_identity_summary.tsv"; }
per_read() { printf 'read_id\tmapped\nr1\ttrue\n' > "$out/per_read_identity.tsv"; }
ladder() { printf 'k\tp\n31\t0.17\n' > "$out/kmer_survival_ladder.tsv"; }
case "${STUB_MODE}" in
  none) ;;
  no_ladder) per_read; summary_good ;;
  no_per_read) summary_good; ladder ;;
  no_summary) per_read; ladder ;;
  unknown) per_read; ladder; { header; row unknown blast; } > "$out/read_identity_summary.tsv" ;;
  empty_first) per_read; ladder; { header; row "" blast; row "Badread v0.4.1" gap; } > "$out/read_identity_summary.tsv" ;;
  empty_second) per_read; ladder; { header; row "Badread v0.4.1" blast; row "" gap; } > "$out/read_identity_summary.tsv" ;;
  header_only) per_read; ladder; header > "$out/read_identity_summary.tsv" ;;
  no_version_col) per_read; ladder; printf 'organism\tsource\nLambda\tblast\n' > "$out/read_identity_summary.tsv" ;;
  good) per_read; ladder; summary_good ;;
esac
exit "${STUB_EXIT:-0}"
STUB
chmod +x "${SCRATCH}/home/.juliaup/bin/julia"

FAILURES=0
run_case() {
  local name="$1" expect="$2" jobid="$3"
  shift 3
  local rc verdict
  (cd "${SCRATCH}/repo" && env HOME="${SCRATCH}/home" USER=selftest SLURM_JOB_ID="${jobid}" \
    SLURM_NODELIST=n0 SLURM_CPUS_PER_TASK=1 ONT_IDENTITY_OUT_ROOT="${SCRATCH}/out" \
    "$@" /bin/bash "${WRAPPER}" > "${SCRATCH}/case-${jobid}.log" 2>&1)
  rc=$?
  verdict=ok
  if { [[ "${expect}" == fail ]] && [[ "${rc}" -eq 0 ]]; } ||
    { [[ "${expect}" == pass ]] && [[ "${rc}" -ne 0 ]]; }; then
    verdict=WRONG
    FAILURES=$((FAILURES + 1))
  fi
  printf '%-34s expect=%-4s rc=%s  %s\n' "${name}" "${expect}" "${rc}" "${verdict}"
}

run_case "no TSVs" fail 101 STUB_MODE=none
run_case "ladder missing" fail 102 STUB_MODE=no_ladder
run_case "per_read missing" fail 103 STUB_MODE=no_per_read
run_case "summary missing" fail 104 STUB_MODE=no_summary
run_case "summary records unknown" fail 105 STUB_MODE=unknown
run_case "badread_version empty, row 1" fail 106 STUB_MODE=empty_first
run_case "badread_version empty, row 2" fail 107 STUB_MODE=empty_second
run_case "summary header only" fail 108 STUB_MODE=header_only
run_case "summary lacks version column" fail 109 STUB_MODE=no_version_col
run_case "julia non-zero, TSVs all good" fail 110 STUB_MODE=good STUB_EXIT=3
run_case "preflight fails" fail 111 STUB_MODE=good STUB_PREFLIGHT=fail
run_case "good run" pass 112 STUB_MODE=good
run_case "reused (non-empty) OUT_DIR" fail 112 STUB_MODE=good

echo
if [[ "${FAILURES}" -eq 0 ]]; then
  echo "ALL CASES CORRECT"
else
  echo "${FAILURES} WRONG"
  exit 1
fi
