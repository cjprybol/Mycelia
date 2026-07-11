#!/bin/bash
#SBATCH --job-name=mycelia_rgv_sweep
#SBATCH --account=m4269
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=rgv_sweep_nersc_%j.out
#SBATCH --error=rgv_sweep_nersc_%j.err

# Rhizomorph correction-core real-genome validation sweep (td-wlfq / td-4e19d) —
# NERSC/Perlmutter variant of run_rgv_sweep_lrc.sbatch. Same sweep, NERSC-site
# directives (account m4269, --constraint=cpu, --qos=regular, julia/1.10.10,
# ~/.julia depot). See the Lawrencium version for full rationale.
#
# Run from benchmarking/:  sbatch run_rgv_sweep_nersc.sbatch
# Prereq (ONCE on a login node, same checkout):
#   module load julia/1.10.10
#   benchmarking/hpc-setup.sh
set -uo pipefail
echo "=== Mycelia Rhizomorph correction-validation sweep (NERSC) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

# Do NOT `module purge`: it strips Cray PE defaults and breaks NERSC's juliaup
# modulefile. Restore PE defaults defensively, then load julia.
source /opt/cray/pe/cpe/25.09/restore_lmod_system_defaults.sh 2>/dev/null || true
module load julia/1.10.10

JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]]; then
  JULIA=/global/common/software/nersc9/julia/1.10.10/bin/julia
fi
if ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable at '${JULIA}' — module load failed or path stale." >&2
  exit 1
fi
echo "JULIA: ${JULIA}  ($("${JULIA}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-${HOME}/.julia}"
export LD_LIBRARY_PATH=""

# === Sweep configuration (real Lambda-phage validation scale) =============
export MYCELIA_RGV_SMOKE="${MYCELIA_RGV_SMOKE:-false}"
export MYCELIA_RGV_ACCESSION="${MYCELIA_RGV_ACCESSION:-NC_001416}"
export MYCELIA_RGV_ERR="${MYCELIA_RGV_ERR:-0.01,0.05,0.10}"
export MYCELIA_RGV_READLEN="${MYCELIA_RGV_READLEN:-150,5000}"
export MYCELIA_RGV_COVERAGE="${MYCELIA_RGV_COVERAGE:-30}"
export MYCELIA_RGV_K="${MYCELIA_RGV_K:-21}"
export MYCELIA_RGV_SEED="${MYCELIA_RGV_SEED:-42}"
export MYCELIA_RUN_EXTERNAL="${MYCELIA_RUN_EXTERNAL:-true}"

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

echo ""
echo ">>> RUN rhizomorph_correction_validation_sweep.jl  $(date)"
echo "    ERR=${MYCELIA_RGV_ERR}  READLEN=${MYCELIA_RGV_READLEN}  COVERAGE=${MYCELIA_RGV_COVERAGE}x  K=${MYCELIA_RGV_K}  QUAST=${MYCELIA_RUN_EXTERNAL}"

if LD_LIBRARY_PATH="" "${JULIA}" --project=.. rhizomorph_correction_validation_sweep.jl; then
  echo ">>> OK  sweep complete  $(date)"; rc=0
else
  echo ">>> FAIL sweep  $(date)" >&2; rc=1
fi

echo ""
echo "=== complete: $(date) ==="
ls -t results/rhizomorph_correction_validation_sweep_*.csv 2>/dev/null | head || true
exit "${rc}"
