#!/bin/bash
#SBATCH --job-name=track_a_baseline
#SBATCH --account=m4269
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=track_a_baseline_nersc_%j.out
#SBATCH --error=track_a_baseline_nersc_%j.err

# Track A baseline benchmark (rhizomorph-paper td-bblmi) — NERSC Perlmutter (CPU).
# Runs the 288-cell viral-tier greedy-Viterbi matrix + power-analysis CV check.
#
# IMPORTANT — Perlmutter compute nodes have NO outbound internet, but the
# benchmark downloads reference genomes (NCBI) and installs Bioconda tools
# (ART/Badread/QUAST) on first use. Both MUST be pre-staged on a LOGIN node
# before submitting this job:
#
#   module load julia/1.10.10
#   export JULIA_DEPOT_PATH=/global/cfs/projectdirs/m4269/cjprybol/.julia_depot
#   export LD_LIBRARY_PATH=""
#   benchmarking/hpc-setup.sh        # once: instantiate/precompile Julia deps
#   # stage refs + conda envs (downloads 4 genomes, installs art/badread/quast):
#   LD_LIBRARY_PATH="" julia --project=. benchmarking/track_a_baseline_benchmark.jl \
#       --organisms Lambda,T4,phi29,SARS-CoV-2 --technologies illumina,pacbio,ont \
#       --coverages 10 --seeds 42 --arms qualmer \
#       --output-dir "$SCRATCH/track_a_baseline"
#
# The runner caches refs (results/.../refs) and conda envs persist in the CFS
# depot, so the compute job below runs fully offline and resumes via per-cell
# checkpoints (re-submit after any timeout).
#
# Submit from benchmarking/:  sbatch run_track_a_baseline_nersc.sbatch
# Smoke / shard override:     TRACK_A_ARGS="--smoke" sbatch run_track_a_baseline_nersc.sbatch
set -uo pipefail
echo "=== Track A baseline benchmark (NERSC Perlmutter) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

module load julia/1.10.10

JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]] || ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable after 'module load julia/1.10.10'." >&2
  exit 1
fi
echo "JULIA: ${JULIA}  ($("${JULIA}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Match the depot used during login-node pre-staging so cached conda envs +
# precompiled deps are found (NERSC default for this project; on CFS, shared
# between login and compute).
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/cfs/projectdirs/m4269/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

# Outputs on scratch (keep $HOME clean); must match the pre-stage --output-dir
# so the downloaded refs are reused offline.
OUTDIR="${TRACK_A_OUTDIR:-$SCRATCH/track_a_baseline}"

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

# shellcheck disable=SC2086  # intentional word-splitting of optional args
echo ">>> RUN track_a_baseline_benchmark.jl --output-dir ${OUTDIR} ${TRACK_A_ARGS:-}  $(date)"
if LD_LIBRARY_PATH="" "${JULIA}" --project=.. track_a_baseline_benchmark.jl \
      --output-dir "${OUTDIR}" ${TRACK_A_ARGS:-}; then
  echo ">>> OK  track_a_baseline_benchmark.jl"
  status=0
else
  echo ">>> FAIL track_a_baseline_benchmark.jl (re-submit to resume from checkpoints)"
  status=1
fi

echo ""
echo "=== complete: $(date) ==="
echo "--- result artifacts ---"
ls -t "${OUTDIR}"/*.tsv "${OUTDIR}"/*.md 2>/dev/null | head || true

exit "${status}"
