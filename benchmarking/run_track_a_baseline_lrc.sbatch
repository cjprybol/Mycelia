#!/bin/bash
#SBATCH --job-name=track_a_baseline
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=track_a_baseline_lrc_%j.out
#SBATCH --error=track_a_baseline_lrc_%j.err

# Track A baseline benchmark (rhizomorph-paper td-bblmi) — Lawrencium.
# Runs the 288-cell viral-tier greedy-Viterbi matrix and the power-analysis CV check.
#
# UNLIKE the talk subset, this job DOES need external tools (ART, Badread, QUAST):
# Mycelia bootstraps them via Bioconda on first use, so the depot must allow conda
# installs (the group depot below is writable). First run pays a one-time env-install
# cost; later runs reuse the envs.
#
# The driver checkpoints every cell to results/track_a_baseline/cells/<id>/cell_result.json,
# so a timeout or crash is fully resumable: just re-submit this script and completed
# cells are skipped.
#
# Run from benchmarking/:  sbatch run_track_a_baseline_lrc.sbatch
# Smoke / shard:           TRACK_A_ARGS="--smoke" sbatch run_track_a_baseline_lrc.sbatch
#                          TRACK_A_ARGS="--organisms Lambda,T4 --arms kmer" sbatch run_track_a_baseline_lrc.sbatch
# Prereq: deps instantiated once on a login node:  benchmarking/hpc-setup.sh
set -uo pipefail
echo "=== Track A baseline benchmark (Lawrencium) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

# The julia/1.10.2-11.4 module DOWNGRADES to a Julia that cannot load current
# master's Manifest (extension-trigger KeyError; td-j8bh) -- the group depot's
# 144 cells that DID succeed under it predate the drift. juliaup's lts channel
# (1.10.10 at last check, matching NERSC) is already on this host. Invoke the
# launcher DIRECTLY with `+lts` rather than relying on PATH + `command -v`:
# the latter silently falls through to whatever julia is already on the
# inherited PATH if juliaup is missing/misconfigured on a given node, and even
# when it finds the juliaup shim, the shim resolves to the host's mutable
# `juliaup default` channel -- not necessarily lts -- so a `juliaup default`
# change made for an unrelated project on the same shared account would
# silently flip every one of these jobs onto a different Julia (td-j8bh
# review). `+lts` bypasses both failure modes: it is resolved by juliaup
# itself, not PATH order, and it is pinned regardless of the default.
JULIA="${HOME}/.juliaup/bin/julia"
JULIA_ARGS=(+lts)

# Fail LOUDLY if the pinned channel is not runnable -- otherwise an empty
# ${JULIA} would run the args as a command and still exit green.
if [[ ! -x "${JULIA}" ]] || ! "${JULIA}" "${JULIA_ARGS[@]}" --version >/dev/null 2>&1; then
  echo "ERROR: '${JULIA} ${JULIA_ARGS[*]}' not runnable (juliaup lts channel missing/broken on this host)." >&2
  exit 1
fi
echo "JULIA: ${JULIA} ${JULIA_ARGS[*]}  ($("${JULIA}" "${JULIA_ARGS[@]}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Group depot (where hpc-setup.sh precompiled deps; conda envs install here too).
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

# shellcheck disable=SC2086  # intentional word-splitting of optional args
echo ">>> RUN track_a_baseline_benchmark.jl ${TRACK_A_ARGS:-}  $(date)"
if LD_LIBRARY_PATH="" "${JULIA}" "${JULIA_ARGS[@]}" --project=.. track_a_baseline_benchmark.jl ${TRACK_A_ARGS:-}; then
  echo ">>> OK  track_a_baseline_benchmark.jl"
  status=0
else
  echo ">>> FAIL track_a_baseline_benchmark.jl (re-submit to resume from checkpoints)"
  status=1
fi

echo ""
echo "=== complete: $(date) ==="
echo "--- result artifacts ---"
ls -t results/track_a_baseline/*.tsv results/track_a_baseline/*.md 2>/dev/null | head || true

exit "${status}"
