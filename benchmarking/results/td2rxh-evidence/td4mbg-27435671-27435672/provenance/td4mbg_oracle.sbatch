#!/bin/bash
#SBATCH --job-name=td4mbg_oracle
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --output=/global/scratch/users/cjprybol/td4mbg_oracle_%j.out
#SBATCH --error=/global/scratch/users/cjprybol/td4mbg_oracle_%j.err

# PR #487 amendment 4, G0 (c): fixed-toy oracle at FIX (dae40fea) then BASE (ca9eb736),
# same shared Manifest, same node. Each arm writes its own output dir; exit codes recorded.
set -uo pipefail
R=/global/home/groups-sw/pc_mfnanofabio/cjprybol/Mycelia/.worktrees
JULIA="${HOME}/.juliaup/bin/julia"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export JULIA_DEPOT_PATH=/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot
export LD_LIBRARY_PATH=""
unset MYCELIA_PROJECT JULIA_PROJECT
export JULIA_LOAD_PATH="@:@stdlib"
tmp_root="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/td4mbg_oracle_${SLURM_JOB_ID}.XXXXXX")" || exit 1
export TMPDIR="$(cd "${tmp_root}" && pwd -P)"
trap 'rm -rf "${tmp_root}"' EXIT
echo "Job ${SLURM_JOB_ID} node ${SLURM_NODELIST} start $(date) TMPDIR=${TMPDIR}"

overall=0
for arm in fix base; do
  W="${R}/td-4mbg-${arm}"
  head="$(git -C "${W}" rev-parse HEAD)" || exit 1
  dirty="$(git -C "${W}" status --porcelain --untracked-files=no)" || exit 1
  [[ -z "${dirty}" ]] || { echo "ERROR: ${arm} worktree dirty: ${dirty}" >&2; exit 1; }
  msha="$(sha256sum "${W}/Manifest.toml" | cut -c1-12)"
  out="/global/scratch/users/cjprybol/td4mbg_oracle/${SLURM_JOB_ID}/${arm}_${head:0:8}"
  mkdir -p "${out}"
  echo ">>> ORACLE ${arm} head=${head} manifest=${msha} out=${out} $(date)"
  (cd "${W}" && "${JULIA}" +lts --startup-file=no --project=. \
      benchmarking/indel_frontier_fixed_toy_proof.jl --output-dir "${out}")
  rc=$?
  echo ">>> ORACLE ${arm} rc=${rc} $(date)"
  [[ "${rc}" -eq 0 ]] || overall=1
done
echo "=== done overall=${overall} $(date) ==="
exit "${overall}"
