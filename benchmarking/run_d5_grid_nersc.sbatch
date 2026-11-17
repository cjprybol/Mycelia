#!/bin/bash
#SBATCH --job-name=d5_grid
#SBATCH --account=m4269
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --array=0-95%16
#SBATCH --output=d5_grid_%A_%a.out
#SBATCH --error=d5_grid_%A_%a.err
#
# D5 ONT adaptive grid (td-yddpr.5) on NERSC Perlmutter (CPU), one array task per
# read set: simulate -> measure identity -> 23-rung survival (once) -> that read
# set's cells, each checkpointed. Resubmitting the same index resumes it.
#
# Resources mirror run_track_a_baseline_nersc.sbatch (shared QOS, 200 GB covers a
# T4 100x :iterative cell); re-size before the full launch from the toy run's
# peak_rss_bytes PLUS the node-local temp root: it is tmpfs, charged to --mem,
# and not in peak_rss_bytes (its size at exit is logged to the .err file).
# Flags on the sbatch command line override the directives above.
#
# PRE-STAGE ONCE on a login node (compute nodes have no outbound internet):
#   module load julia/1.10.10
#   export JULIA_DEPOT_PATH=/global/cfs/projectdirs/m4269/cjprybol/.julia_depot
#   export LD_LIBRARY_PATH=""
#   julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
#   julia --project=. -e 'import Mycelia; foreach(Mycelia.add_bioconda_env,
#       ["badread", "art", "minimap2", "quast"])'
#   # references: --list does not download; Phase A downloads into <root>/refs on
#   # first use, so run one toy task interactively (add --tmp-dir /tmp/$USER-d5:
#   # the default <root>/tmp is on Lustre and is refused) or pre-copy refs/.
#
# USAGE (submit from benchmarking/):
#   Toy (3 read sets):
#     D5_ARGS="--toy" D5_OUTDIR="$SCRATCH/d5_grid_toy" \
#       sbatch --array=0-2 --time=04:00:00 run_d5_grid_nersc.sbatch
#   Full Phase A (96 read sets) -- only after the user approves the launch:
#     D5_OUTDIR="$SCRATCH/d5_grid" sbatch run_d5_grid_nersc.sbatch
#   Merge after an array finishes (one task):
#     D5_PHASE=merge D5_OUTDIR=... sbatch --array=0 --time=01:00:00 \
#       --dependency=afterany:<array_job_id> run_d5_grid_nersc.sbatch
#   Phase B -- only after a human sets tau / tau_hat (lock step 2):
#     D5_PHASE=B D5_ARGS="--tau <x> --tau-hat <y>" D5_OUTDIR=... \
#       sbatch run_d5_grid_nersc.sbatch
set -uo pipefail

PHASE="${D5_PHASE:-A}"
OUTDIR="${D5_OUTDIR:?set D5_OUTDIR to the grid output root (one root per grid)}"
EXTRA_ARGS="${D5_ARGS:-}"
TASK_INDEX="${SLURM_ARRAY_TASK_ID:-0}"

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
# The :iterative corrector removes its temp directories with
# renameat2(RENAME_NOREPLACE). Lustre ($SCRATCH) rejects that with EINVAL and the
# corrector stops, so its temp root is node-local: /tmp on a Perlmutter compute
# node is tmpfs (it counts against --mem) and is not a symlink. The runner
# re-checks the rename before any cell runs. Override the base with D5_TMPDIR.
# mktemp creates a fresh 0700 directory and fails if the name exists, so a
# pre-planted path can never become what the EXIT trap removes.
D5_TMP="$(mktemp --directory "${D5_TMPDIR:-/tmp}/d5_${SLURM_JOB_ID:-$$}.XXXXXX")" ||
  { echo "ERROR: cannot create a temp root under ${D5_TMPDIR:-/tmp}" >&2; exit 1; }
cleanup_temp_root() {
  local rc=$?
  echo "temp root ${D5_TMP} at exit: $(du --summarize --bytes "${D5_TMP}" 2>/dev/null | cut --fields=1) bytes" >&2
  # On failure keep a copy: the corrector leaves fail-closed evidence under
  # TMPDIR, cell error rows cite those paths, and the node is gone after the job.
  if [[ "${rc}" -ne 0 && -n "$(ls --almost-all "${D5_TMP}" 2>/dev/null)" ]]; then
    local keep="${OUTDIR}/tmp_evidence/${SLURM_JOB_ID:-$$}"
    if mkdir --parents "${keep}" && cp --archive "${D5_TMP}/." "${keep}/"; then
      echo "temp root preserved at ${keep}" >&2
    else
      echo "WARN: could not preserve ${D5_TMP}" >&2
    fi
  fi
  rm --recursive --force -- "${D5_TMP}"
}
trap cleanup_temp_root EXIT
D5_TMP="$(realpath "${D5_TMP}")" || { echo "ERROR: cannot resolve the temp root" >&2; exit 1; }
export TMPDIR="${D5_TMP}"
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

INDEX_ARGS=()
if [[ "${PHASE}" == "A" || "${PHASE}" == "B" ]]; then
  INDEX_ARGS=(--read-set-index "${TASK_INDEX}")
fi

# Exit codes from the runner: 0 ok; 3 = cells recorded but some need a resubmit
# (A/B) or the merged grid is incomplete or mixed (merge); anything else = the
# runner itself failed (bad args, manifest refusal, simulator or reference error).
# shellcheck disable=SC2086  # D5_ARGS is intentionally word-split
"${JULIA}" --project=.. d5_ont_adaptive_grid_run.jl --phase "${PHASE}" \
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
