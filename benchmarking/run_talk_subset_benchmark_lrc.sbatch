#!/bin/bash
#SBATCH --job-name=mycelia_talk_bench
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=talk_bench_lrc_%j.out
#SBATCH --error=talk_bench_lrc_%j.err

# Self-contained Julia benchmark subset for the LBNL talk (td-756dg) — Lawrencium.
# Pure-Julia steps only (no conda/external-tool deps): k-mer scaling + data processing.
#
# Run from benchmarking/:  sbatch run_talk_subset_benchmark_lrc.sbatch
# Prereq: deps instantiated once on a login node (repairs a stale manifest that
# is missing HDF5, then precompiles):  benchmarking/hpc-setup.sh
set -uo pipefail
echo "=== Mycelia talk-subset benchmark (Lawrencium) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

module load julia/1.10.2-11.4

# Resolve the julia binary explicitly so the job fails loudly if the module is broken.
JULIA="$(command -v julia || true)"
echo "JULIA: ${JULIA}"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Use the group depot (where hpc-setup.sh precompiled the deps); also set by .bashrc.
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

run_one () {
  local f="$1"
  echo ""
  echo ">>> RUN ${f}  $(date)"
  if LD_LIBRARY_PATH="" "${JULIA}" --project=.. "${f}"; then
    echo ">>> OK  ${f}"
  else
    echo ">>> FAIL ${f} (continuing)"
  fi
}

run_one 02_kmer_analysis_benchmark.jl
run_one 01_data_processing_benchmark.jl

echo ""
echo "=== complete: $(date) ==="
echo "--- result JSONs ---"
# `|| true`: a non-matching glob makes `ls` exit non-zero, which `pipefail`
# would otherwise propagate as the job's exit code (false FAILED in sacct).
ls -t results/*.json 2>/dev/null | head || true

# Per-benchmark status is reported via the >>> OK/FAIL markers above; exit clean
# so SLURM's job state reflects the wrapper running to completion.
exit 0
