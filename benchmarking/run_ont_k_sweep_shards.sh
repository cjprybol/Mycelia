#!/usr/bin/env bash
# Parallel driver for benchmarking/ont_k_sweep.jl.
#
# The 96-cell grid is embarrassingly parallel — every cell is independent and
# writes its own checkpoint directory — but run serially it is dominated by the
# ONT/100x cells, where a single assembly takes tens of minutes. This partitions
# the grid by (k, coverage) into 16 shards of 6 cells each (2 technologies x 3
# seeds) and runs them concurrently, so wall-clock is set by the slowest single
# shard rather than by the sum of all of them.
#
# SAFETY OF THE SHARED OUTPUT TREE
# --------------------------------
# All shards share one --output-dir on purpose, so checkpoints accumulate in a
# single tree and a resumed run sees every completed cell. That is safe because:
#
#   * Shards partition the grid, so no two ever compute the same cell_id and no
#     two ever write the same cells/<id>/ directory.
#   * The reference download is guarded by `if !isfile(outfile)` in
#     download_genome_by_accession, and this script pre-warms it serially below
#     so no shard ever races on it.
#
# The ONE thing that is NOT safe under concurrency is the aggregate TSV: every
# shard rewrites ont_k_sweep_results.tsv after each of its cells, so during the
# run that file is whatever the last writer produced and is NOT a complete or
# necessarily well-formed table. It is scratch. The final aggregation pass at
# the end runs ONE process over the full grid, hits every checkpoint on the
# resume path, and regenerates both the results and summary TSVs from the
# per-cell JSONs. Only that regenerated pair is the deliverable.
#
# Usage:
#   ./benchmarking/run_ont_k_sweep_shards.sh
#   OUTPUT_DIR=/scratch/ont_k_sweep ./benchmarking/run_ont_k_sweep_shards.sh
#
# Re-running is safe and cheap: completed cells are skipped via their
# checkpoints, so this doubles as the crash-recovery path.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO_ROOT}/benchmarking/results/ont_k_sweep}"
LOG_DIR="${LOG_DIR:-${OUTPUT_DIR}/shard-logs}"

KS="${KS:-11 15 21 31}"
COVERAGES="${COVERAGES:-10 30 50 100}"
# Sharding by technology as well as (k, coverage) gives 32 shards instead of 16.
# On a large host that halves wall-clock, because the grid's long pole is the
# three ONT cells at each high coverage and splitting the technologies stops
# them from queueing behind their Illumina siblings. On a laptop, set
# TECHNOLOGIES="illumina,ont" to get the coarser 16-shard layout back.
TECHNOLOGIES="${TECHNOLOGIES:-illumina ont}"

mkdir -p "${LOG_DIR}"

echo "=== ONT k-sweep parallel driver ==="
echo "repo:       ${REPO_ROOT}"
echo "output dir: ${OUTPUT_DIR}"
echo "log dir:    ${LOG_DIR}"

# Pre-warm the reference serially. Every shard calls
# download_genome_by_accession, which is idempotent but would otherwise have 16
# processes reach the not-yet-downloaded state simultaneously on a cold tree.
echo "--- pre-warming reference (serial) ---"
julia --project="${REPO_ROOT}" "${REPO_ROOT}/benchmarking/ont_k_sweep.jl" \
    --technologies illumina --ks 31 --coverages 10 --seeds 42 \
    --output-dir "${OUTPUT_DIR}" > "${LOG_DIR}/prewarm.log" 2>&1
echo "reference ready"

echo "--- launching shards ---"
pids=""
for k in ${KS}; do
    for coverage in ${COVERAGES}; do
        for technology in ${TECHNOLOGIES}; do
            log="${LOG_DIR}/shard_${technology}_k${k}_${coverage}x.log"
            julia --project="${REPO_ROOT}" "${REPO_ROOT}/benchmarking/ont_k_sweep.jl" \
                --technologies "${technology}" --ks "${k}" --coverages "${coverage}" \
                --output-dir "${OUTPUT_DIR}" > "${log}" 2>&1 &
            pids="${pids} $!"
            echo "  launched ${technology} k=${k} coverage=${coverage}x (pid $!) -> ${log}"
        done
    done
done

echo "--- waiting for $(echo ${pids} | wc -w | tr -d ' ') shards ---"
failed=0
for pid in ${pids}; do
    if ! wait "${pid}"; then
        echo "  shard pid ${pid} exited non-zero"
        failed=$((failed + 1))
    fi
done
echo "shards complete (${failed} non-zero exits)"

# Single-process aggregation over the full grid. Every cell is already
# checkpointed, so this is a pure resume pass: it recomputes nothing and writes
# the authoritative results + summary TSVs.
echo "--- final aggregation (serial, resume-only) ---"
julia --project="${REPO_ROOT}" "${REPO_ROOT}/benchmarking/ont_k_sweep.jl" \
    --output-dir "${OUTPUT_DIR}" 2>&1 | tail -40

echo "=== done ==="
