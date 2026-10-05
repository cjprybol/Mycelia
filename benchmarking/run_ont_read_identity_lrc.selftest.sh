#!/bin/bash
# Self-test for run_ont_read_identity_lrc.sbatch: runs the wrapper against a
# stub `julia` in a scratch HOME and repo, so every success and failure path can
# be exercised without a cluster, conda, or a real Julia.
#
# Each case isolates ONE guard: the stub writes exactly the artifacts that case
# needs, so a case can only fail for the reason it names. Exits non-zero if any
# case returns the wrong verdict.
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
# --version prints the version on stdout and an update notice on stderr (as
# juliaup can); -e is the preflight; otherwise it is the main run, which writes
# into --output-dir according to STUB_MODE and exits with STUB_EXIT.
mkdir -p "${SCRATCH}/home/.juliaup/bin"
cat > "${SCRATCH}/home/.juliaup/bin/julia" <<'STUB'
#!/bin/bash
for a in "$@"; do
  if [[ "$a" == "--version" ]]; then
    echo "julia version 1.10.10"
    echo "notice: a newer lts version is available" >&2
    exit 0
  fi
done
for a in "$@"; do
  [[ "$a" == "-e" ]] && { [[ "${STUB_PREFLIGHT:-ok}" == "fail" ]] && exit 4; exit 0; }
done
out=""; prev=""
for a in "$@"; do [[ "$prev" == "--output-dir" ]] && out="$a"; prev="$a"; done
summary() { printf 'organism\tbadread_version\tsource\nLambda\t%s\tblast\n' "$1" > "$out/read_identity_summary.tsv"; }
per_read() { printf 'read_id\tmapped\nr1\ttrue\n' > "$out/per_read_identity.tsv"; }
ladder() { printf 'k\tp\n31\t0.17\n' > "$out/kmer_survival_ladder.tsv"; }
case "${STUB_MODE}" in
  none) ;;
  no_ladder) per_read; summary "Badread v0.4.1" ;;
  no_per_read) summary "Badread v0.4.1"; ladder ;;
  no_summary) per_read; ladder ;;
  unknown) per_read; ladder; summary "unknown" ;;
  empty_version) per_read; ladder; summary "" ;;
  good) per_read; ladder; summary "Badread v0.4.1" ;;
esac
exit "${STUB_EXIT:-0}"
STUB
chmod +x "${SCRATCH}/home/.juliaup/bin/julia"

# --- scratch repos --------------------------------------------------------------
make_repo() {
  mkdir -p "$1/benchmarking"
  : > "$1/benchmarking/ont_read_identity.jl"
  printf '[deps]\n' > "$1/Manifest.toml"
}
make_repo "${SCRATCH}/repo"
git -C "${SCRATCH}/repo" init -q
git -C "${SCRATCH}/repo" add benchmarking
git -C "${SCRATCH}/repo" -c user.name=selftest -c user.email=selftest@invalid commit -q -m init
make_repo "${SCRATCH}/norepo"   # same layout, not a git repository

# A git that resolves HEAD but fails `status` (e.g. a safe.directory refusal
# that only bites some subcommands). Isolates the tree-state fallback from the
# HEAD fallback, which "not a git repo" trips first.
REAL_GIT="$(command -v git)"
mkdir -p "${SCRATCH}/gitshim"
# shellcheck disable=SC2016  # $1/$@ must stay literal: they belong to the shim
printf '#!/bin/bash\n[[ "$1" == status ]] && exit 128\nexec "%s" "$@"\n' "${REAL_GIT}" \
  > "${SCRATCH}/gitshim/git"
chmod +x "${SCRATCH}/gitshim/git"

FAILURES=0
LAST_OUT=""
run_case() {
  local name="$1" expect="$2" repo="$3" jobid="$4"
  shift 4
  local rc verdict
  (cd "${repo}" && env HOME="${SCRATCH}/home" USER=selftest SLURM_JOB_ID="${jobid}" \
    SLURM_NODELIST=n0 SLURM_CPUS_PER_TASK=1 ONT_IDENTITY_OUT_ROOT="${SCRATCH}/out" \
    "$@" /bin/bash "${WRAPPER}" > "${SCRATCH}/case-${jobid}.log" 2>&1)
  rc=$?
  LAST_OUT="${SCRATCH}/out/${jobid}"
  verdict=ok
  if { [[ "${expect}" == fail ]] && [[ "${rc}" -eq 0 ]]; } ||
    { [[ "${expect}" == pass ]] && [[ "${rc}" -ne 0 ]]; }; then
    verdict=WRONG
    FAILURES=$((FAILURES + 1))
  fi
  printf '%-34s expect=%-4s rc=%s  %s\n' "${name}" "${expect}" "${rc}" "${verdict}"
}
check() {
  local label="$1"
  shift
  if "$@"; then
    printf '%-34s ok\n' "${label}"
  else
    printf '%-34s WRONG\n' "${label}"
    FAILURES=$((FAILURES + 1))
  fi
}

R="${SCRATCH}/repo"
run_case "no TSVs" fail "$R" 101 STUB_MODE=none
run_case "ladder missing" fail "$R" 102 STUB_MODE=no_ladder
run_case "per_read missing" fail "$R" 103 STUB_MODE=no_per_read
run_case "summary missing" fail "$R" 104 STUB_MODE=no_summary
run_case "summary records unknown" fail "$R" 105 STUB_MODE=unknown
run_case "summary badread_version empty" fail "$R" 106 STUB_MODE=empty_version
run_case "julia non-zero, TSVs all good" fail "$R" 107 STUB_MODE=good STUB_EXIT=3
run_case "not a git repo" fail "${SCRATCH}/norepo" 108 STUB_MODE=good
run_case "git status fails, HEAD resolves" fail "$R" 112 STUB_MODE=good PATH="${SCRATCH}/gitshim:${PATH}"
check "  status failure: tree state unknown" grep -q '^git_tracked_files_modified=unknown$' "${LAST_OUT}/PROVENANCE"
run_case "preflight fails" fail "$R" 109 STUB_MODE=good STUB_PREFLIGHT=fail
check "  preflight: exit_status=1 recorded" grep -q '^exit_status=1$' "${LAST_OUT}/PROVENANCE"
run_case "good run" pass "$R" 110 STUB_MODE=good
GOOD="${LAST_OUT}/PROVENANCE"
check "  good: git_head is a sha" grep -q -E '^git_head=[0-9a-f]{40}$' "${GOOD}"
check "  good: tree clean" grep -q '^git_tracked_files_modified=0$' "${GOOD}"
check "  good: manifest fingerprinted" grep -q -E '^manifest_sha256=[0-9a-f]{64}$' "${GOOD}"
check "  good: julia on one line" grep -q '^julia=julia version 1.10.10$' "${GOOD}"
check "  good: every line is key=value" awk '!/^[a-z_0-9]+=/{bad=1} END{exit bad}' "${GOOD}"
check "  good: exit_status=0 recorded" grep -q '^exit_status=0$' "${GOOD}"
run_case "reused (non-empty) OUT_DIR" fail "$R" 110 STUB_MODE=good
printf 'x\n' > "$R/benchmarking/ont_read_identity.jl"
run_case "dirty tree still runs" pass "$R" 111 STUB_MODE=good
check "  dirty: modification counted" grep -q '^git_tracked_files_modified=1$' "${LAST_OUT}/PROVENANCE"

echo
if [[ "${FAILURES}" -eq 0 ]]; then
  echo "ALL CASES CORRECT"
else
  echo "${FAILURES} WRONG"
  exit 1
fi
