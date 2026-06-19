#!/bin/bash
#SBATCH --job-name=mycelia_talk_bench
#SBATCH --account=m4269
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=talk_bench_%j.out
#SBATCH --error=talk_bench_%j.err

# Self-contained Julia benchmark subset for the LBNL talk (td-756dg).
# Pure-Julia steps only (no conda/external-tool deps): k-mer scaling + data processing.
#
# NERSC/Perlmutter wrapper. Run from benchmarking/:  sbatch run_talk_subset_benchmark.sbatch
# Prereq: Mycelia deps instantiated/precompiled once on a login node:
#   cd ~/workspace/Mycelia && LD_LIBRARY_PATH="" julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
set -uo pipefail
echo "=== Mycelia talk-subset benchmark ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

# Do NOT `module purge` here: it strips the Cray PE defaults, after which NERSC's
# juliaup modulefile errors ("attempt to index a nil value") and `julia` vanishes
# from PATH. Restore the PE defaults defensively, then load julia.
source /opt/cray/pe/cpe/25.09/restore_lmod_system_defaults.sh 2>/dev/null || true
module load julia/1.10.10

# Resolve the julia binary explicitly so the job fails loudly if the module is broken.
JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]]; then
  JULIA=/global/common/software/nersc9/julia/1.10.10/bin/julia
fi
echo "JULIA: ${JULIA}  ($(${JULIA} --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
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
