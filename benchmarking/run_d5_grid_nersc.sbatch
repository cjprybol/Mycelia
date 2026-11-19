#!/bin/bash
#SBATCH --job-name=d5_grid
#SBATCH --account=m4269
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --array=0-2
#SBATCH --output=d5_grid_%A_%a.out
#SBATCH --error=d5_grid_%A_%a.err
#
# The #SBATCH --array directive above is the toy/pilot (read sets 0-2). The
# full arrays are not directives: pass `--array=...` and D5_FULL=1 together.
# See USAGE below.
#
# D5 ONT adaptive grid (td-yddpr.5) on NERSC Perlmutter (CPU). Phase A is two
# arrays, every cell checkpointed; resubmitting an index resumes it:
#   D5_PHASE=A        one task per read set: simulate -> measure identity ->
#                     23-rung survival (once) -> the read set's fast cells
#                     (:none, and k the assembler rejects). Passes --skip-iterative.
#   D5_PHASE=A-cells  one task per :iterative cell (--cell-index), on a read set
#                     its D5_PHASE=A task has prepared; submit it with
#                     --dependency=afterany on the read-set array.
# One :iterative cell took 8.9 h (Lambda, ONT 99%/100x, toy job 61473595_1), so
# a read set's seven or more :iterative cells cannot share one task. --time is
# 48 h, the shared QOS maximum; shorten it on the command line for read-set,
# merge and toy tasks.
#
# Resources mirror run_track_a_baseline_nersc.sbatch (shared QOS, 200 GB). The
# Lambda 100x :iterative cell peaked at 55.9 GB RSS (sacct MaxRSS 64,657,632 KiB,
# 61.7 GiB; sacct reports KiB); T4 is 3.5x larger and unmeasured. Re-size from the cgroup memory.peak logged to the
# .err file at exit. peak_rss_bytes alone is too low: the node-local temp root
# is tmpfs, charged to --mem, not to RSS.
# Flags on the sbatch command line override the directives above.
#
# The #SBATCH --array default is the toy/pilot (0-2). A bare
# `sbatch run_d5_grid_nersc.sbatch` therefore submits three tasks, not the
# 96-task Phase A array. Any submitted array whose max index or task id is past
# the toy is refused before module load unless D5_FULL=1: above 2 for read-set
# tasks, above 3 for cell tasks (the toy has 4 :iterative cells). That includes
# the full commands below. D5_FULL does not resize --mem or --time.
#
# PRE-STAGE ONCE on a login node (compute nodes have no outbound internet):
#   module load julia/1.10.10
#   export JULIA_DEPOT_PATH=/global/cfs/projectdirs/m4269/cjprybol/.julia_depot
#   export LD_LIBRARY_PATH=""
#   julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
#   julia --project=. -e 'import Mycelia; foreach(Mycelia.add_bioconda_env,
#       ["badread", "art", "minimap2", "quast"])'
#   # references: --list does not download; Phase A downloads into <root>/refs on
#   # first use, so run one toy read-set task interactively with
#   # --skip-iterative (without it the read set's :iterative cells run on the
#   # login node, up to 8.9 h each), adding --tmp-dir "$(mktemp --directory)"
#   # (the default <root>/tmp is on Lustre and is refused), or pre-copy refs/.
#
# USAGE (submit from benchmarking/):
#   The default array is indices 0-2; only D5_ARGS="--toy" makes them the toy
#   read sets (without it they are the first three full-grid read sets, at the
#   48 h directive unless --time is given). The launch host is chosen at launch
#   by allocation burn-down (todo hpc-burndown.py); this file is the NERSC entry
#   point. Read-set task walltime is unmeasured above Lambda 100x (simulation,
#   identity and survival of up to 17 Mb of reads); a timed-out read-set task
#   restarts its read set from scratch.
#   Toy: 3 read sets (the default array), then its 4 :iterative cells:
#     D5_ARGS="--toy" D5_OUTDIR="$SCRATCH/d5_grid_toy" \
#       sbatch --parsable --time=02:00:00 run_d5_grid_nersc.sbatch
#     D5_PHASE=A-cells D5_ARGS="--toy" D5_OUTDIR="$SCRATCH/d5_grid_toy" \
#       sbatch --array=0-3 --dependency=afterany:<read_set_job_id> \
#       run_d5_grid_nersc.sbatch
#   Full Phase A, only after the user approves the launch. Both the array
#   override and D5_FULL=1 are required; either one alone is not the full array
#   (the directive default is 0-2, and a wider array without D5_FULL=1 exits 1
#   in every task). `--list` prints both array sizes (96 read sets, 720 cells):
#     D5_FULL=1 D5_OUTDIR="$SCRATCH/d5_grid" \
#       sbatch --parsable --array=0-95%16 --time=06:00:00 run_d5_grid_nersc.sbatch
#     D5_FULL=1 D5_PHASE=A-cells D5_OUTDIR="$SCRATCH/d5_grid" \
#       sbatch --array=0-719%<cap> --dependency=afterany:<read_set_job_id> \
#       run_d5_grid_nersc.sbatch
#   A cell task whose read set is not prepared exits 1 and writes no row;
#   resubmit its index once the read-set task has finished.
#   Merge after both arrays finish (one task; index 0 is inside the toy):
#     D5_PHASE=merge D5_OUTDIR=... sbatch --array=0 --time=01:00:00 \
#       --dependency=afterany:<cell_job_id> run_d5_grid_nersc.sbatch
#   Phase B toy, only after a human sets tau / tau_hat (lock step 2):
#     D5_PHASE=B D5_ARGS="--tau <x> --tau-hat <y>" D5_OUTDIR=... \
#       sbatch --array=0-2 run_d5_grid_nersc.sbatch
#   Phase B full uses the same D5_FULL=1 and --array=0-95%16 gate as Phase A.
set -uo pipefail

PHASE="${D5_PHASE:-A}"
EXTRA_ARGS="${D5_ARGS:-}"
TASK_INDEX="${SLURM_ARRAY_TASK_ID:-0}"

# Toy/pilot is array indices 0-2 (read sets), or 0-3 for A-cells (the toy's
# :iterative cells). Refuse anything wider before module load, mkdir, or
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

echo "=== D5 grid (Perlmutter) phase=${PHASE} task=${TASK_INDEX} ==="
echo "Job: ${SLURM_JOB_ID:-none}  Array: ${SLURM_ARRAY_JOB_ID:-none}  Node: ${SLURM_NODELIST:-none}  Start: $(date)"

module load julia/1.10.10
JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]] || ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable after 'module load julia/1.10.10'." >&2
  exit 1
fi

# No `set -e` here, so every step that defines a root must fail loudly: an
# empty OUTDIR would otherwise make the runner write the grid into the checkout.
mkdir -p "${OUTDIR}" || { echo "ERROR: cannot create D5_OUTDIR=${OUTDIR}" >&2; exit 1; }
OUTDIR="$(realpath "${OUTDIR}")" || { echo "ERROR: cannot resolve D5_OUTDIR" >&2; exit 1; }
[[ -n "${OUTDIR}" ]] || { echo "ERROR: D5_OUTDIR resolved to an empty path" >&2; exit 1; }
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/cfs/projectdirs/m4269/cjprybol/.julia_depot}"
export LD_LIBRARY_PATH=""

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${SCRIPT_DIR}" || exit 1
if [[ ! -f d5_ont_adaptive_grid_run.jl ]]; then
  echo "ERROR: submit from the benchmarking/ directory of the checkout to run." >&2
  exit 1
fi
echo "Checkout: $(git -C .. rev-parse HEAD)  dirty: $(git -C .. status --porcelain --untracked-files=no | wc -l | tr -d ' ') files"

# The :iterative corrector removes its temp directories with
# renameat2(RENAME_NOREPLACE). Lustre ($SCRATCH) rejects that with EINVAL and the
# corrector stops, so its temp root is node-local: /tmp on a Perlmutter compute
# node is tmpfs (it counts against --mem) and is not a symlink. The runner
# re-checks the rename before any cell runs. Override the base with D5_TMPDIR.
# mktemp creates a fresh 0700 directory and fails if the name exists, so a
# pre-planted path can never become what the EXIT trap removes. Cleanup, the
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
