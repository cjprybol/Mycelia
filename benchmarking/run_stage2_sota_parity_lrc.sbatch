#!/bin/bash
#SBATCH --job-name=rhizo_stage2_parity
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=rhizo_stage2_parity_lrc_%j.out
#SBATCH --error=rhizo_stage2_parity_lrc_%j.err

# Stage-2 SOTA-parity benchmark (rhizomorph td-3e54 / td-dgpb / td-pw7p) — Lawrencium.
#
# Runs benchmarking/real_data_corrector_validation.jl on the REPEAT-BEARING
# E. coli K-12 fixture with the standalone SPAdes SOTA arm enabled, emitting
# NGA50 + misassembly-proxy for every arm (naive / scalable / hybrid-olc /
# sota-spades). This produces the evidence a "contiguity parity vs SOTA with
# superior per-base accuracy" comparison would rest on. NO parity claim is
# asserted by running this — the committed CSV is interim engineering evidence.
#
# External tools (ART, MEGAHIT, SPAdes, MUMmer/dnadiff) bootstrap via Bioconda on
# first use into the group depot below (writable), so the first run pays a
# one-time env-install cost; later runs reuse the envs.
#
# WARNING: this harness is NOT per-arm checkpointed — a timeout loses the whole
# run. 24h is generous for 4 arms on a 4.64 Mb genome @ 50x plus first-run conda
# installs; SPAdes is the long pole. Re-submit to restart from scratch.
#
# Prereq (once, on a login node): benchmarking/hpc-setup.sh instantiates deps so
# Julia can load Mycelia (avoids the "Arrow not installed" first-run failure).
# Run from benchmarking/:  sbatch run_stage2_sota_parity_lrc.sbatch
# Override targets/scale:  RDV_TARGETS=lambda RDV_SOTA=false sbatch run_stage2_sota_parity_lrc.sbatch
set -uo pipefail
echo "=== Rhizomorph Stage-2 SOTA-parity benchmark (Lawrencium) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

module load julia/1.10.2-11.4

# Resolve julia explicitly and fail LOUDLY if the module did not put it on PATH.
JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]] || ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable after 'module load julia/1.10.2-11.4'." >&2
  exit 1
fi
echo "JULIA: ${JULIA}  ($("${JULIA}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Group depot (where hpc-setup.sh precompiled deps; conda envs install here too).
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

# Repeat-bearing fixture + SPAdes SOTA arm on by default; overridable via env.
export MYCELIA_RDV_TARGETS="${RDV_TARGETS:-ecoli_k12}"
export MYCELIA_RDV_SOTA="${RDV_SOTA:-true}"
export MYCELIA_RDV_COVERAGE="${RDV_COVERAGE:-50}"
export MYCELIA_RDV_SEED="${RDV_SEED:-42}"

echo ">>> targets=${MYCELIA_RDV_TARGETS} sota=${MYCELIA_RDV_SOTA} " \
     "coverage=${MYCELIA_RDV_COVERAGE}x seed=${MYCELIA_RDV_SEED}  $(date)"
if LD_LIBRARY_PATH="" "${JULIA}" --project=.. real_data_corrector_validation.jl; then
  echo ">>> OK  real_data_corrector_validation.jl"
  status=0
else
  echo ">>> FAIL real_data_corrector_validation.jl"
  status=1
fi

echo ""
echo "=== complete: $(date) ==="
echo "--- result artifacts ---"
ls -td results/real_data_corrector_validation_* 2>/dev/null | head -1 || true

exit "${status}"
