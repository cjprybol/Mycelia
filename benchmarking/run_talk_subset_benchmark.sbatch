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
# Prereq: deps instantiated once on a login node (repairs a stale manifest that
# is missing HDF5, then precompiles):  benchmarking/hpc-setup.sh
#   (load the julia module first, e.g. `module load julia/1.10.10`)
set -uo pipefail
echo "=== Mycelia talk-subset benchmark ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

# Do NOT `module purge` here: it strips the Cray PE defaults, after which NERSC's
# juliaup modulefile errors ("attempt to index a nil value") and `julia` vanishes
# from PATH. Restore the PE defaults defensively, then load julia.
source /opt/cray/pe/cpe/25.09/restore_lmod_system_defaults.sh 2>/dev/null || true
module load julia/1.10.10

# Resolve the julia binary explicitly, with a fallback if the module is broken.
JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]]; then
  JULIA=/global/common/software/nersc9/julia/1.10.10/bin/julia
fi
# Validate it actually runs — a broken `module load` (NERSC's juliaup modulefile
# can error) or a relocated fallback path must fail the job LOUDLY here, not
# silently produce a green job with no results downstream.
if ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable at '${JULIA}' — module load failed or path stale." >&2
  exit 1
fi
echo "JULIA: ${JULIA}  ($("${JULIA}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Pin the depot so the batch job reuses what hpc-setup.sh precompiled on the login
# node (SLURM does not reliably inherit interactive-shell exports). `:-` keeps any
# value already exported by the environment/.bashrc.
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-${HOME}/.julia}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

fail_count=0
run_one () {
  local f="$1"
  echo ""
  echo ">>> RUN ${f}  $(date)"
  if LD_LIBRARY_PATH="" "${JULIA}" --project=.. "${f}"; then
    echo ">>> OK  ${f}"
  else
    echo ">>> FAIL ${f} (continuing)"
    fail_count=$((fail_count + 1))
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

# Exit non-zero if any benchmark failed so SLURM's job state is truthful — a
# partial run dropping a result must NOT read as COMPLETED. (The empty-glob above
# is neutralized separately by `|| true`, so it cannot flip the job red.)
if [[ "${fail_count}" -gt 0 ]]; then
  echo "=== ${fail_count} benchmark(s) FAILED ==="
  exit 1
fi
exit 0
