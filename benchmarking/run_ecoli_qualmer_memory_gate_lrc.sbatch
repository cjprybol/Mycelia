#!/bin/bash
#SBATCH --job-name=ecoli_qualmer_memgate
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=ecoli_qualmer_memgate_%j.out
#SBATCH --error=ecoli_qualmer_memgate_%j.err

# td-n8ax E. coli qualmer-path MEMORY gate on lr6 (~95GB node — the node the
# un-fixed runs OOM'd on at 74-200GB). GATE: run COMPLETES and peak RSS <= 32GB
# (aspiration 16GB) with the COMPOSED config (:ultralight_quality aggregate storage
# + qualmer_prefilter_min_count=2). Driver: benchmarking/ecoli_qualmer_memory_gate.jl.
#
# Run from the branch worktree on lrc:
#   cd ~/workspace/Mycelia/.worktrees/td-n8ax-gate/benchmarking
#   sbatch run_ecoli_qualmer_memory_gate_lrc.sbatch
# Mycelia bootstraps ART + QUAST via Bioconda into the group depot on first use.
set -uo pipefail
echo "=== E. coli qualmer memory gate (Lawrencium lr6) ==="
echo "Job ${SLURM_JOB_ID}  Node ${SLURM_NODELIST}  Start $(date)"

module load julia/1.10.2-11.4
JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]] || ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable after module load" >&2
  exit 1
fi
echo "JULIA: ${JULIA}  ($("${JULIA}" --version 2>&1 | head -1))"

# Group depot (precompiled deps + Bioconda envs install here).
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
export LD_LIBRARY_PATH=""
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

cd "${HOME}/workspace/Mycelia/.worktrees/td-n8ax-gate/benchmarking" || exit 1

echo ">>> instantiate (branch may recompile Mycelia from source) $(date)"
if ! LD_LIBRARY_PATH="" "${JULIA}" --project=.. -e 'import Pkg; Pkg.instantiate()'; then
  echo ">>> instantiate FAILED" >&2
  exit 1
fi

echo ">>> RUN ecoli_qualmer_memory_gate.jl $(date)"
LD_LIBRARY_PATH="" "${JULIA}" --project=.. ecoli_qualmer_memory_gate.jl
status=$?

echo "=== gate exit=${status} complete $(date) ==="
# exit 0 = PASS (completes, <=32GB); 1 = memory FAIL (OOM/over-ceiling); 2 = code error.
exit "${status}"
