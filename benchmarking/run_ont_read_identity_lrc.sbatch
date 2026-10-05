#!/bin/bash
#SBATCH --job-name=ont_read_identity
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=ont_read_identity_lrc_%j.out
#SBATCH --error=ont_read_identity_lrc_%j.err

# ONT read-identity measurement (td-4e19d.29) — Lawrencium.
#
# Regenerates the three TSVs under results/ont_read_identity/ with the script's
# DEFAULT parameters (Lambda NC_001416, simulated ONT reads at 30x, seed 42),
# which are the parameters of the committed run. Do not pass --coverage/--seed
# here unless the intent is a DIFFERENT measurement.
#
# Output goes to a fresh per-job directory, never to the repo's results/
# directory: the script writes fixed filenames, so pointing two runs at one
# directory would let the later silently replace the earlier. Copy the three
# TSVs into the repo deliberately once they have been checked.
#
# Submit from the repo root (so --project=. resolves):
#   sbatch --chdir="$WT" benchmarking/run_ont_read_identity_lrc.sbatch
# Override the output root:
#   ONT_IDENTITY_OUT_ROOT=/some/dir sbatch --chdir="$WT" benchmarking/run_ont_read_identity_lrc.sbatch
# Prereq: deps instantiated once on a login node (JULIA_PKG_PRECOMPILE_AUTO=0).
set -uo pipefail
echo "=== ONT read-identity measurement (Lawrencium) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"
echo "Repo: $(pwd)  HEAD: $(git rev-parse HEAD 2>/dev/null || echo unknown)"

# Pinned juliaup lts channel; see run_track_a_baseline_lrc.sbatch for why the
# julia module and a bare PATH lookup are both avoided (td-j8bh).
JULIA="${HOME}/.juliaup/bin/julia"
JULIA_ARGS=(+lts)
if [[ ! -x "${JULIA}" ]] || ! "${JULIA}" "${JULIA_ARGS[@]}" --version >/dev/null 2>&1; then
  echo "ERROR: '${JULIA} ${JULIA_ARGS[*]}' not runnable (juliaup lts channel missing/broken on this host)." >&2
  exit 1
fi
echo "JULIA: ${JULIA} ${JULIA_ARGS[*]}  ($("${JULIA}" "${JULIA_ARGS[@]}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

if [[ ! -f benchmarking/ont_read_identity.jl ]]; then
  echo "ERROR: run from the repo root (sbatch --chdir=<repo>); not found: $(pwd)/benchmarking/ont_read_identity.jl" >&2
  exit 1
fi

OUT_ROOT="${ONT_IDENTITY_OUT_ROOT:-/global/scratch/users/${USER}/ont_read_identity}"
OUT_DIR="${OUT_ROOT}/${SLURM_JOB_ID}"
mkdir -p "${OUT_DIR}" || { echo "ERROR: cannot create ${OUT_DIR}" >&2; exit 1; }
echo "Output: ${OUT_DIR}"

# The script records `badread --version` BEFORE it simulates, and simulation is
# what installs the env. On a cold depot that probe would record "unknown" in
# read_identity_summary.tsv, so the envs are provisioned first.
echo ">>> PREFLIGHT conda envs  $(date)"
if ! "${JULIA}" "${JULIA_ARGS[@]}" --project=. -e \
  'import Mycelia; Mycelia.add_bioconda_env("badread"); Mycelia.add_bioconda_env("minimap2")'; then
  echo ">>> FAIL preflight (conda env provisioning)" >&2
  exit 1
fi

echo ">>> RUN ont_read_identity.jl  $(date)"
if "${JULIA}" "${JULIA_ARGS[@]}" --project=. benchmarking/ont_read_identity.jl --output-dir "${OUT_DIR}"; then
  echo ">>> OK  ont_read_identity.jl"
  status=0
else
  echo ">>> FAIL ont_read_identity.jl"
  status=1
fi

echo ""
echo "=== complete: $(date) ==="
echo "--- result artifacts ---"
ls -l "${OUT_DIR}"/*.tsv || true

exit "${status}"
