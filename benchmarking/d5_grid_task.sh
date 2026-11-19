# One D5 grid array task (td-yddpr.5): the body shared by the per-site sbatch
# wrappers (run_d5_grid_nersc.sbatch, run_d5_grid_lrc.sbatch). Sourced, not
# executed; it always ends in `exit`. GNU coreutils (long flags).
#
# The wrapper sets, before sourcing:
#   D5_SITE_NAME   site label for the log
#   D5_SITE_DEPOT  default JULIA_DEPOT_PATH on that site
#   D5_BENCH_DIR   the checkout's benchmarking/ directory
#   d5_load_julia  function that puts julia 1.10 on PATH (returns non-zero if it
#                  cannot)
# Environment (both sites): D5_PHASE (A | A-cells | B | merge), D5_OUTDIR,
# D5_ARGS, D5_FULL, D5_TMPDIR. See the wrappers for USAGE.
set -uo pipefail

PHASE="${D5_PHASE:-A}"
EXTRA_ARGS="${D5_ARGS:-}"
TASK_INDEX="${SLURM_ARRAY_TASK_ID:-0}"

# Toy/pilot is array indices 0-2 (read sets), or 0-3 for A-cells (the toy's
# :iterative cells). Refuse anything wider before loading julia, mkdir, or
# mktemp, so a mis-pasted full submit runs nothing: each task still queues and
# starts briefly, then exits 1 here.
# D5_FULL must be exactly 1. SLURM_ARRAY_TASK_MAX is the submitted array's last
# index (a bare sbatch uses the #SBATCH --array default of 2). A non-integer
# bound fails closed.
case "${PHASE}" in
  A | B | merge) D5_TOY_MAX=2 ;;
  A-cells) D5_TOY_MAX=3 ;;
  *) echo "ERROR: D5_PHASE must be A, A-cells, B or merge, got ${PHASE}" >&2; exit 1 ;;
esac
if [[ "${D5_FULL:-0}" != "1" ]]; then
  ARRAY_MAX="${SLURM_ARRAY_TASK_MAX:-${TASK_INDEX}}"
  if ! [[ "${ARRAY_MAX}" =~ ^[0-9]+$ && "${TASK_INDEX}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: D5 array bounds must be non-negative integers (task=${TASK_INDEX}, array max=${ARRAY_MAX}). Set D5_FULL=1 only for an approved wider array." >&2
    exit 1
  fi
  if (( ARRAY_MAX > D5_TOY_MAX || TASK_INDEX > D5_TOY_MAX )); then
    echo "ERROR: D5 array task ${TASK_INDEX} (array max ${ARRAY_MAX}) is outside the toy/pilot 0-${D5_TOY_MAX} for D5_PHASE=${PHASE}. Set D5_FULL=1 with the approved --array (see USAGE) to run past the toy. A bare sbatch submits the toy only." >&2
    exit 1
  fi
fi

OUTDIR="${D5_OUTDIR:?set D5_OUTDIR to the grid output root (one root per grid)}"

echo "=== D5 grid (${D5_SITE_NAME}) phase=${PHASE} task=${TASK_INDEX} ==="
echo "Job: ${SLURM_JOB_ID:-none}  Array: ${SLURM_ARRAY_JOB_ID:-none}  Node: ${SLURM_NODELIST:-none}  CPUs: ${SLURM_CPUS_PER_TASK:-unset}  Start: $(date)"

d5_load_julia
JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]] || ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable after d5_load_julia (${D5_SITE_NAME})." >&2
  exit 1
fi
echo "Julia: ${JULIA} ($("${JULIA}" --version 2>&1))"

# No `set -e` here, so every step that defines a root must fail loudly: an
# empty OUTDIR would otherwise make the runner write the grid into the checkout.
mkdir -p "${OUTDIR}" || { echo "ERROR: cannot create D5_OUTDIR=${OUTDIR}" >&2; exit 1; }
OUTDIR="$(realpath "${OUTDIR}")" || { echo "ERROR: cannot resolve D5_OUTDIR" >&2; exit 1; }
[[ -n "${OUTDIR}" ]] || { echo "ERROR: D5_OUTDIR resolved to an empty path" >&2; exit 1; }
# Four Julia threads on every site, as in the toy run (array 61473595), so the
# grid is one configuration whatever --cpus-per-task a site's wrapper requests.
export JULIA_NUM_THREADS=4
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-${D5_SITE_DEPOT}}"
export LD_LIBRARY_PATH=""

cd "${D5_BENCH_DIR}" || exit 1
if [[ ! -f d5_ont_adaptive_grid_run.jl ]]; then
  echo "ERROR: submit from the benchmarking/ directory of the checkout to run." >&2
  exit 1
fi
echo "Checkout: $(git -C .. rev-parse HEAD)  dirty: $(git -C .. status --porcelain --untracked-files=no | wc -l | tr -d ' ') files"

# The :iterative corrector removes its temp directories with
# renameat2(RENAME_NOREPLACE). Lustre (Perlmutter $SCRATCH, Lawrencium
# /global/scratch) rejects that with EINVAL and the corrector stops, so its temp
# root is node-local /tmp (on Perlmutter, tmpfs that counts against --mem). The
# runner re-checks the rename before any cell runs. Override the base with
# D5_TMPDIR. mktemp creates a fresh 0700 directory and fails if the name exists,
# so a pre-planted path can never become what the EXIT trap removes. Cleanup, the
# evidence copy on failure and signal handling live in d5_temp_root.sh.
# shellcheck source=d5_temp_root.sh
source ./d5_temp_root.sh || { echo "ERROR: cannot source d5_temp_root.sh" >&2; exit 1; }
# The helper is read from the live checkout at run time, not snapshotted at
# submit: refuse to run without its functions rather than run without traps.
declare -F d5_install_temp_root_traps d5_temp_root_cleanup >/dev/null ||
  { echo "ERROR: d5_temp_root.sh does not define the temp-root functions" >&2; exit 1; }
[[ "${D5_TMPDIR:-/tmp}" == /* ]] ||
  { echo "ERROR: D5_TMPDIR must be an absolute path, got ${D5_TMPDIR}" >&2; exit 1; }
D5_TMP_MADE="$(mktemp --directory "${D5_TMPDIR:-/tmp}/d5_${SLURM_JOB_ID:-$$}.XXXXXX")" ||
  { echo "ERROR: cannot create a temp root under ${D5_TMPDIR:-/tmp}" >&2; exit 1; }
# A requeued job keeps its SLURM_JOB_ID, so the restart count keeps each
# attempt's evidence apart (an existing evidence dir is never overwritten).
d5_install_temp_root_traps "${D5_TMP_MADE}" \
  "${OUTDIR}/tmp_evidence/${SLURM_JOB_ID:-$$}.r${SLURM_RESTART_COUNT:-0}" ||
  { trap - EXIT TERM INT HUP; rmdir -- "${D5_TMP_MADE}"
    echo "ERROR: cannot install the temp-root traps" >&2; exit 1; }
D5_TMP="$(realpath "${D5_TMP_MADE}")" || { echo "ERROR: cannot resolve the temp root" >&2; exit 1; }
export TMPDIR="${D5_TMP}"

# The runner's --phase, plus the index flag that selects this task's work.
INDEX_ARGS=()
RUNNER_PHASE="${PHASE}"
case "${PHASE}" in
  A) INDEX_ARGS=(--read-set-index "${TASK_INDEX}" --skip-iterative) ;;
  A-cells) RUNNER_PHASE=A; INDEX_ARGS=(--cell-index "${TASK_INDEX}") ;;
  B) INDEX_ARGS=(--read-set-index "${TASK_INDEX}") ;;
esac

# Exit codes from the runner: 0 ok; 3 = cells recorded but some need a resubmit
# (A/B) or the merged grid is incomplete or mixed (merge); anything else = the
# runner itself failed (bad args, manifest refusal, simulator or reference
# error, or a cell task whose read set is not prepared yet).
# shellcheck disable=SC2086  # D5_ARGS is intentionally word-split
"${JULIA}" --project=.. d5_ont_adaptive_grid_run.jl --phase "${RUNNER_PHASE}" \
  --output-dir "${OUTDIR}" --tmp-dir "${D5_TMP}" \
  "${INDEX_ARGS[@]+"${INDEX_ARGS[@]}"}" ${EXTRA_ARGS}
status=$?
if [[ "${status}" -eq 0 ]]; then
  echo ">>> OK phase=${PHASE} task=${TASK_INDEX}  $(date)"
elif [[ "${status}" -eq 3 && "${PHASE}" == "merge" ]]; then
  echo ">>> BLOCKED merge: read merge_provenance.json 'blocking_reasons' (absent/retryable cells -> resubmit those indices; mixed tau or code versions -> a resubmit will not help)  $(date)"
elif [[ "${status}" -eq 3 ]]; then
  echo ">>> RETRY phase=${PHASE} task=${TASK_INDEX}: some cells need a resubmit of this index  $(date)"
else
  echo ">>> FAIL phase=${PHASE} task=${TASK_INDEX} (exit ${status}); read the .err log before resubmitting  $(date)"
fi
exit "${status}"
