#!/bin/bash
#SBATCH --job-name=mycelia_rgv_sweep
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=rgv_sweep_lrc_%j.out
#SBATCH --error=rgv_sweep_lrc_%j.err

# Rhizomorph correction-core real-genome validation sweep (td-wlfq / td-4e19d) —
# Lawrencium. Runs benchmarking/rhizomorph_correction_validation_sweep.jl on the
# REAL Lambda-phage reference (NC_001416, ~48.5 kb) across the
# (error_rate x read-regime) x {naive, iterative} matrix at validation-scale
# coverage, so the scale guard emits a VERDICT (not SMOKE-ONLY).
#
# This answers the manuscript's central preliminary question: does iterative
# correction improve final ASSEMBLY quality over naive across error regimes, and
# in which regime does correction deliver value?
#
# Run from benchmarking/:  sbatch run_rgv_sweep_lrc.sbatch
# Prereq (ONCE on a login node, in the same checkout):
#   module load julia/1.10.2-11.4
#   benchmarking/hpc-setup.sh          # resolve + instantiate + precompile
set -uo pipefail
echo "=== Mycelia Rhizomorph correction-validation sweep (Lawrencium) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

module load julia/1.10.2-11.4

# Resolve the julia binary explicitly and fail LOUDLY if the module did not put
# julia on PATH — otherwise an empty ${JULIA} would run the sweep args as a
# command, fail, and still exit green with no results.
JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]] || ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable after 'module load julia/1.10.2-11.4'." >&2
  exit 1
fi
echo "JULIA: ${JULIA}  ($("${JULIA}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Reuse the group depot where hpc-setup.sh precompiled the deps (SLURM does not
# reliably inherit interactive-shell exports). `:-` keeps any pre-exported value.
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

# === Sweep configuration (real Lambda-phage validation scale) =============
# Override any of these on the sbatch line with `--export=ALL,VAR=val,...`.
export MYCELIA_RGV_SMOKE="${MYCELIA_RGV_SMOKE:-false}"          # FULL run: download NC_001416
export MYCELIA_RGV_ACCESSION="${MYCELIA_RGV_ACCESSION:-NC_001416}"
export MYCELIA_RGV_ERR="${MYCELIA_RGV_ERR:-0.01,0.05,0.10}"     # low / mid / high error regimes
export MYCELIA_RGV_READLEN="${MYCELIA_RGV_READLEN:-150,5000}"   # short-Illumina / long-Nanopore
export MYCELIA_RGV_COVERAGE="${MYCELIA_RGV_COVERAGE:-30}"       # 48.5kb x 30x ~= 1.46 Mbase > scale floor
export MYCELIA_RGV_K="${MYCELIA_RGV_K:-21}"
export MYCELIA_RGV_SEED="${MYCELIA_RGV_SEED:-42}"
# Enable QUAST (Mycelia auto-provisions it via Conda.jl); harness degrades
# gracefully to internal metrics if QUAST/Conda is unavailable.
export MYCELIA_RUN_EXTERNAL="${MYCELIA_RUN_EXTERNAL:-true}"

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

echo ""
echo ">>> RUN rhizomorph_correction_validation_sweep.jl  $(date)"
echo "    ERR=${MYCELIA_RGV_ERR}  READLEN=${MYCELIA_RGV_READLEN}  COVERAGE=${MYCELIA_RGV_COVERAGE}x  K=${MYCELIA_RGV_K}  QUAST=${MYCELIA_RUN_EXTERNAL}"

if LD_LIBRARY_PATH="" "${JULIA}" --project=.. rhizomorph_correction_validation_sweep.jl; then
  echo ">>> OK  sweep complete  $(date)"
  rc=0
else
  echo ">>> FAIL sweep  $(date)" >&2
  rc=1
fi

echo ""
echo "=== complete: $(date) ==="
echo "--- result CSVs ---"
ls -t results/rhizomorph_correction_validation_sweep_*.csv 2>/dev/null | head || true

exit "${rc}"
