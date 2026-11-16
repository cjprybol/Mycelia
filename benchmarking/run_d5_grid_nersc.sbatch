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
# T4 100x :iterative cell); re-size from the toy run's peak_rss_bytes before the
# full launch. Flags on the sbatch command line override the directives above.
#
# PRE-STAGE ONCE on a login node (compute nodes have no outbound internet):
#   module load julia/1.10.10
#   export JULIA_DEPOT_PATH=/global/cfs/projectdirs/m4269/cjprybol/.julia_depot
#   export LD_LIBRARY_PATH=""
#   julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
#   julia --project=. -e 'import Mycelia; foreach(Mycelia.add_bioconda_env,
#       ["badread", "art", "minimap2", "quast"])'
#   # references: --list does not download; Phase A downloads into <root>/refs on
#   # first use, so run one toy task interactively or pre-copy refs/ into <root>.
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

# shellcheck disable=SC2086  # D5_ARGS is intentionally word-split
if "${JULIA}" --project=.. d5_ont_adaptive_grid_run.jl --phase "${PHASE}" \
     --output-dir "${OUTDIR}" "${INDEX_ARGS[@]+"${INDEX_ARGS[@]}"}" ${EXTRA_ARGS}; then
  echo ">>> OK phase=${PHASE} task=${TASK_INDEX}  $(date)"
  status=0
else
  echo ">>> FAIL phase=${PHASE} task=${TASK_INDEX} (resubmit this index to resume)  $(date)"
  status=1
fi
exit "${status}"
