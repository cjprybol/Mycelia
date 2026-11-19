#!/bin/bash
#SBATCH --job-name=d5_grid
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --array=0-2
#SBATCH --output=d5_grid_%A_%a.out
#SBATCH --error=d5_grid_%A_%a.err
#
# LBNL Lawrencium wrapper for one D5 grid array task (td-yddpr.5). The task
# itself (toy gate, temp root, runner call) is d5_grid_task.sh, shared with
# run_d5_grid_nersc.sbatch; this file holds only Lawrencium's directives, its
# julia, and its depot. Phases, the toy gate, D5_FULL and the two-array Phase A
# are as described in run_d5_grid_nersc.sbatch.
#
# Why these directives:
#   lr7 shares nodes (OverSubscribe=FORCE:1; 56 cores, 257-503 GB). lr6 is
#   whole-node only with 95 GB on most nodes, below a T4 100x :iterative cell.
#   lr_normal allows 72 h, above the projected ~41 h of the slowest cell
#   (T4 90-95%/100x :iterative, extrapolated from the toy run).
#   8 CPUs for 4 Julia threads (pinned in d5_grid_task.sh): the toy's tasks
#   averaged ~9.7 busy cores each (sacct TotalCPU / Elapsed) with 4 threads.
#   --mem=200G is the un-sized envelope; lower it per stratum from the cgroup
#   memory.peak logged at exit.
# /global/scratch is Lustre, which rejects the corrector's rename, so the temp
# root is node-local /tmp (D5_TMPDIR overrides); the runner checks it first.
#
# PRE-STAGE ONCE on a login node:
#   export JULIA_DEPOT_PATH=/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot
#   export LD_LIBRARY_PATH=""
#   julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
#   julia --project=. -e 'import Mycelia; foreach(Mycelia.add_bioconda_env,
#       ["badread", "art", "minimap2", "quast"])'
#   julia is juliaup's (~/.juliaup/bin, channel lts = 1.10.10 as of 2026-11-19);
#   every row records julia_version, so a juliaup update mid-grid is visible.
#
# USAGE (submit from benchmarking/; output roots under /global/scratch/users/$USER):
#   Toy: 3 read sets (the default array), then its 4 :iterative cells:
#     D5_ARGS="--toy" D5_OUTDIR=/global/scratch/users/$USER/d5_grid_toy \
#       sbatch --parsable --time=02:00:00 run_d5_grid_lrc.sbatch
#     D5_PHASE=A-cells D5_ARGS="--toy" \
#       D5_OUTDIR=/global/scratch/users/$USER/d5_grid_toy \
#       sbatch --array=0-3 --time=24:00:00 \
#       --dependency=afterany:<read_set_job_id> run_d5_grid_lrc.sbatch
#   Full Phase A, only after the user approves the launch (D5_FULL=1 and the
#   array override together, as on Perlmutter):
#     D5_FULL=1 D5_OUTDIR=/global/scratch/users/$USER/d5_grid \
#       sbatch --parsable --array=0-95 --time=06:00:00 --mem=64G \
#       run_d5_grid_lrc.sbatch
#     D5_FULL=1 D5_PHASE=A-cells D5_OUTDIR=/global/scratch/users/$USER/d5_grid \
#       sbatch --array=0-719%<cap> --dependency=afterany:<read_set_job_id> \
#       run_d5_grid_lrc.sbatch
#   Merge after both arrays finish:
#     D5_PHASE=merge D5_OUTDIR=... sbatch --array=0 --time=01:00:00 \
#       --dependency=afterany:<cell_job_id> run_d5_grid_lrc.sbatch
set -uo pipefail

D5_SITE_NAME="Lawrencium"
D5_SITE_DEPOT="/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot"
# juliaup's julia, unless one is already on PATH.
d5_load_julia() {
  command -v julia >/dev/null 2>&1 || export PATH="${HOME}/.juliaup/bin:${PATH}"
}
# sbatch runs a spooled copy of this file, so the shared body is found from the
# submit directory (benchmarking/); `bash <this file>` finds it next to the file.
D5_BENCH_DIR="${SLURM_SUBMIT_DIR:-$(dirname -- "${BASH_SOURCE[0]}")}"
# shellcheck source=d5_grid_task.sh
source "${D5_BENCH_DIR}/d5_grid_task.sh"
# d5_grid_task.sh always exits; reaching this line means it could not be read.
echo "ERROR: cannot run ${D5_BENCH_DIR}/d5_grid_task.sh; submit from the benchmarking/ directory of the checkout." >&2
exit 1
