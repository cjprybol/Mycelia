#!/bin/bash
#SBATCH --job-name=d5_grid
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --array=0-2
#SBATCH --output=d5_grid_%A_%a.out
#SBATCH --error=d5_grid_%A_%a.err
#
# LBNL Lawrencium wrapper for one D5 grid array task (td-yddpr.5). The task
# itself (toy gate, julia version check, temp root, runner call) is
# d5_grid_task.sh, shared with run_d5_grid_nersc.sbatch; this file holds only
# Lawrencium's directives, its julia, and its depot. Phases, the toy gate,
# D5_FULL and the two-array Phase A are as described in run_d5_grid_nersc.sbatch.
#
# Why these directives (cluster facts read with sinfo / scontrol / sacctmgr on
# 2026-11-19):
#   lr7 shares nodes (OverSubscribe=FORCE:1; 56 cores, 257-503 GB). lr6 is
#   whole-node only, with 95 GB on most nodes: under the 113 GB peak RSS
#   projected (not measured) for a T4 100x :iterative cell.
#   lr_normal allows 72 h, above the ~41 h projected (not measured) for the
#   slowest cell, T4 90-95%/100x :iterative, from the toy run.
#   12 CPUs for 4 Julia threads (pinned in d5_grid_task.sh): the toy's tasks
#   averaged ~9.7 busy cores each (sacct TotalCPU / Elapsed), with GC, QUAST
#   (2 threads) and minimap2 beside the Julia threads, on Perlmutter nodes that
#   were not CPU-limited. Fewer CPUs would throttle and stretch the 41 h.
#   --mem=200G is the un-sized envelope; lower it per stratum from the cgroup
#   memory.peak logged at exit. KillWait is 30 s, as on Perlmutter.
# /global/scratch is Lustre, so the temp root is node-local /tmp (D5_TMPDIR
# overrides); the runner checks the corrector's rename there before any cell.
#
# PRE-STAGE ONCE on a login node (from the checkout root):
#   export JULIA_DEPOT_PATH=/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot
#   export LD_LIBRARY_PATH=""
#   julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
#   julia --project=. -e 'import Mycelia; foreach(Mycelia.add_bioconda_env,
#       ["badread", "art", "minimap2", "quast"])'
#   # references: copy <toy root>/refs/ into a new full root before its read-set
#   # array starts, so 96 tasks do not download the same files at once.
#   julia is juliaup's (~/.juliaup/bin, channel lts = 1.10.10 on 2026-11-19).
#   d5_grid_task.sh refuses any version but 1.10.10, so a juliaup update, or
#   another julia first on the submit shell's PATH, stops the task instead of
#   changing the grid.
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
#   array override together, as on Perlmutter). Read-set task time and memory
#   are unmeasured above Lambda 100x; keep the 200G envelope until measured:
#     D5_FULL=1 D5_OUTDIR=/global/scratch/users/$USER/d5_grid \
#       sbatch --parsable --array=0-95%16 --time=06:00:00 run_d5_grid_lrc.sbatch
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
[[ -r "${D5_BENCH_DIR}/d5_grid_task.sh" ]] || {
  echo "ERROR: no readable ${D5_BENCH_DIR}/d5_grid_task.sh; submit from the benchmarking/ directory of the checkout." >&2
  exit 1
}
# shellcheck source=d5_grid_task.sh
source "${D5_BENCH_DIR}/d5_grid_task.sh"
# d5_grid_task.sh always ends in exit; reaching this line means it stopped early.
echo "ERROR: ${D5_BENCH_DIR}/d5_grid_task.sh returned without exiting (status $?)." >&2
exit 1
