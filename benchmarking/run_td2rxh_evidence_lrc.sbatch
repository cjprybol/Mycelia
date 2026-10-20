#!/bin/bash
#SBATCH --job-name=mycelia_td2rxh
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-8
#SBATCH --time=1-00:00:00
#SBATCH --output=td2rxh_lrc_%A_%a.out
#SBATCH --error=td2rxh_lrc_%A_%a.err

# td-2rxh pre-registered evidence run on Lawrencium: does the merged sparse-rung
# indel staging (#420) plus the #450 reroute fix make the :nanopore corrector
# profile beat :illumina at 1000bp+? Runs
# indel_reroute_evidence.jl --cells td2rxh (6 kb x {1000,2000} bp x 20x x
# err {0.01,0.02,0.05}). One array task = one seed over all six cells, so both
# arms of every pair run back to back on the same node, which keeps the paired
# wall-time difference (H1') free of node-to-node variation.
#
# Usage (from benchmarking/ of the checkout to measure):
#   sbatch run_td2rxh_evidence_lrc.sbatch                    # seeds 1-8
#   sbatch --array=1 --export=ALL,TD2RXH_CELLS=smoke run_td2rxh_evidence_lrc.sbatch  # smoke
# Then analyse with:
#   julia +lts --project=.. td2rxh_paired_analysis.jl \
#     --csv "$(ls results/td2rxh/<ARRAY_JOB_ID>/seed*.csv | paste -sd, -)" \
#     --output-dir results/td2rxh/<ARRAY_JOB_ID>/analysis
#
# Prereq (once on a login node, same checkout): benchmarking/hpc-setup.sh.
set -uo pipefail
echo "=== td-2rxh evidence (Lawrencium) ==="
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

# Pinned juliaup lts launcher; see run_rgv_sweep_lrc.sbatch (td-j8bh) for why a
# PATH lookup or `module load julia` is not used.
JULIA="${HOME}/.juliaup/bin/julia"
JULIA_ARGS=(+lts)
if [[ ! -x "${JULIA}" ]] || ! "${JULIA}" "${JULIA_ARGS[@]}" --version >/dev/null 2>&1; then
  echo "ERROR: '${JULIA} ${JULIA_ARGS[*]}' not runnable on this host." >&2
  exit 1
fi

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
export LD_LIBRARY_PATH=""

# Measure the checkout this was SUBMITTED from, never a hardcoded path: the
# provenance below must describe the code that actually ran.
cd "${SLURM_SUBMIT_DIR}" || exit 1
if [[ ! -f indel_reroute_evidence.jl ]]; then
  echo "ERROR: submit from a checkout's benchmarking/ directory (got ${SLURM_SUBMIT_DIR})." >&2
  exit 1
fi
commit="$(git -C .. rev-parse HEAD)"
dirty="$(git -C .. status --porcelain --untracked-files=no | wc -l | tr -d ' ')"
echo "Checkout: $(cd .. && pwd)  HEAD=${commit}  tracked-dirty-files=${dirty}"
if [[ "${dirty}" != "0" ]]; then
  echo "ERROR: tracked files are modified; results would not be attributable to ${commit}." >&2
  exit 1
fi

cells="${TD2RXH_CELLS:-td2rxh}"
seed="${SLURM_ARRAY_TASK_ID}"
outdir="results/td2rxh/${SLURM_ARRAY_JOB_ID}"
mkdir -p "${outdir}"
out="${outdir}/seed${seed}.csv"

echo ">>> RUN cells=${cells} seed=${seed} -> ${out}  $(date)"
if "${JULIA}" "${JULIA_ARGS[@]}" --project=.. indel_reroute_evidence.jl \
    --cells "${cells}" --seeds "${seed}" --out "${out}" \
    --label "td2rxh@${commit:0:8}"; then
  rc=0
else
  rc=1
fi

# A green exit with a header-only CSV means every arm threw (each failure is a
# @warn + skip in the harness). Treat that as a failure, not a result.
rows=0
[[ -f "${out}" ]] && rows=$(( $(wc -l < "${out}") - 1 ))
expected=12
[[ "${cells}" == "smoke" ]] && expected=4
echo "rows=${rows} expected=${expected}"
if [[ "${rows}" -lt "${expected}" ]]; then
  echo ">>> FAIL: ${out} has ${rows} of ${expected} arm rows (arms threw; see .err)" >&2
  rc=1
fi
echo "=== done rc=${rc} $(date) ==="
exit "${rc}"
