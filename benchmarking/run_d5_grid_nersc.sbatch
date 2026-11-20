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
# NERSC Perlmutter (CPU) wrapper for one D5 grid array task (td-yddpr.5). The
# task itself (toy gate, temp root, runner call) is d5_grid_task.sh, shared with
# run_d5_grid_lrc.sbatch; this file holds only Perlmutter's directives, its
# julia, and its depot.
#
# The #SBATCH --array directive above is the toy/pilot (read sets 0-2). The
# full arrays are not directives: pass `--array=...` and D5_FULL=1 together.
# See USAGE below.
#
# Phase A is two arrays, every cell checkpointed; resubmitting an index resumes it:
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
# Lambda 100x :iterative cell's sampled peak_rss_bytes was 55.9 GB; the task's
# sacct MaxRSS was 64,657,632 KiB (66.2 GB; sacct reports KiB). T4 is 3.5x
# larger and unmeasured. Re-size from the cgroup memory.peak logged to the .err
# file at exit. peak_rss_bytes alone is too low: the node-local temp root is
# tmpfs, charged to --mem, not to RSS.
# On the shared QOS --mem also sets the CPUs charged: 200G allocated 108 of 256.
# Flags on the sbatch command line override the directives above.
#
# The #SBATCH --array default is the toy/pilot (0-2). A bare
# `sbatch run_d5_grid_nersc.sbatch` therefore submits three tasks, not the
# 96-task Phase A array. Any submitted array whose max index or task id is past
# the toy is refused before julia loads unless D5_FULL=1: above 2 for read-set
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
#   # references: compute nodes have no outbound internet, and Phase A refuses
#   # to download. Stage once on a login node (per-accession lock, then rename):
#   #   julia --project=. benchmarking/d5_ont_adaptive_grid_run.jl --stage-refs \
#   #       --toy --output-dir "$SCRATCH/d5_grid_toy"
#   # or copy an existing refs/ into the output root. Array tasks only read it.
#
# USAGE (submit from benchmarking/):
#   The default array is indices 0-2; only D5_ARGS="--toy" makes them the toy
#   read sets (without it they are the first three full-grid read sets, at the
#   48 h directive unless --time is given). The launch host is chosen at launch
#   by allocation burn-down (todo hpc-burndown.py). Read-set task walltime is
#   unmeasured above Lambda 100x (simulation, identity and survival of up to
#   17 Mb of reads); a timed-out read-set task restarts its read set from scratch.
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

D5_SITE_NAME="Perlmutter"
D5_SITE_DEPOT="/global/cfs/projectdirs/m4269/cjprybol/.julia_depot"
d5_load_julia() { module load julia/1.10.10; }
# sbatch runs a spooled copy of this file, so the shared body is found from the
# submit directory (benchmarking/); `bash <this file>` finds it next to the file.
D5_BENCH_DIR="${SLURM_SUBMIT_DIR:-$(dirname -- "${BASH_SOURCE[0]}")}"
[[ -r "${D5_BENCH_DIR}/d5_grid_task.sh" ]] || {
  echo "ERROR: no readable ${D5_BENCH_DIR}/d5_grid_task.sh; submit from the benchmarking/ directory of the checkout." >&2
  exit 1
}
# shellcheck source=d5_grid_task.sh
source "${D5_BENCH_DIR}/d5_grid_task.sh"
# d5_grid_task.sh always ends in exit; reaching this line means it stopped early
# (a return, or a syntax error, which bash reports just above).
echo "ERROR: ${D5_BENCH_DIR}/d5_grid_task.sh stopped without exiting (status $?)." >&2
exit 1
