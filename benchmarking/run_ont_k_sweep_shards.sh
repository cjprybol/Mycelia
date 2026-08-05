#!/usr/bin/env bash
# Parallel driver for benchmarking/ont_k_sweep.jl.
#
# The 96-cell grid is embarrassingly parallel — every cell is independent and
# writes its own checkpoint directory — but run serially it is dominated by the
# ONT/100x cells, where a single assembly takes tens of minutes. This partitions
# the grid by (organism, k, coverage, technology) and runs the shards
# concurrently, so wall-clock is set by the slowest single shard rather than by
# the sum of all of them. With the defaults that is 32 shards of 3 cells each
# (3 seeds); see the TECHNOLOGIES note below for the coarser 16-shard layout.
#
# SAFETY OF THE SHARED OUTPUT TREE
# --------------------------------
# All shards share one --output-dir on purpose, so checkpoints accumulate in a
# single tree and a resumed run sees every completed cell. That is safe because:
#
#   * Shards partition the grid, so no two ever compute the same cell_id and no
#     two ever write the same cells/<id>/ directory.
#   * The reference is fetched by the SERIAL PRE-WARM below before any shard
#     starts. That pre-warm is load-bearing and must not be deleted:
#     download_genome_by_accession guards on `if !isfile(outfile)`, which is
#     check-then-act and NOT concurrency-safe on its own — N cold shards can all
#     observe the file missing and race. The pre-warm also warms the conda
#     environments the shards need, for the same reason (add_bioconda_env has no
#     locking, so simultaneous `conda create` calls can collide).
#
# The aggregate TSVs are rewritten by every shard after each of its cells, so
# during the run they are mid-write and should not be read. They cannot SHRINK,
# though: write_aggregate unions the in-memory rows with every checkpoint on
# disk, so a shard covering 3 cells still emits the whole tree. The final pass
# below re-emits them from the checkpoints alone.
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
# Organisms default to Lambda, matching the sweep script's own default. T4 is
# opt-in because it is 3.5x larger and its ONT high-coverage cells dominate the
# grid's cost: ORGANISMS="Lambda T4" ./run_ont_k_sweep_shards.sh
ORGANISMS="${ORGANISMS:-Lambda}"

mkdir -p "${LOG_DIR}"

echo "=== ONT k-sweep parallel driver ==="
echo "repo:       ${REPO_ROOT}"
echo "output dir: ${OUTPUT_DIR}"
echo "log dir:    ${LOG_DIR}"

# Pre-warm the reference serially. Every shard calls
# download_genome_by_accession, which is idempotent but would otherwise have 16
# processes reach the not-yet-downloaded state simultaneously on a cold tree.
echo "--- pre-warming reference (serial) ---"
# Pre-warm BOTH technologies: the illumina cell installs `art` + `quast`, and
# the ont cell installs `badread`. Warming only illumina left every ONT shard to
# race on `conda create -n badread`, which has no locking.
julia --project="${REPO_ROOT}" "${REPO_ROOT}/benchmarking/ont_k_sweep.jl" \
    --organisms "$(echo ${ORGANISMS} | tr ' ' ',')" \
    --technologies illumina,ont --ks 31 --coverages 10 --seeds 42 \
    --output-dir "${OUTPUT_DIR}" > "${LOG_DIR}/prewarm.log" 2>&1
echo "references and conda environments ready"

echo "--- launching shards ---"
pids=""
for organism in ${ORGANISMS}; do
    for k in ${KS}; do
        for coverage in ${COVERAGES}; do
            for technology in ${TECHNOLOGIES}; do
                log="${LOG_DIR}/shard_${organism}_${technology}_k${k}_${coverage}x.log"
                julia --project="${REPO_ROOT}" "${REPO_ROOT}/benchmarking/ont_k_sweep.jl" \
                    --organisms "${organism}" --technologies "${technology}" \
                    --ks "${k}" --coverages "${coverage}" \
                    --output-dir "${OUTPUT_DIR}" > "${log}" 2>&1 &
                pids="${pids} $!"
                echo "  launched ${organism} ${technology} k=${k} coverage=${coverage}x (pid $!) -> ${log}"
            done
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

# Single-process aggregation with --aggregate-only, which reads the cells
# directory and computes NOTHING.
#
# Passing a grid here instead would be wrong, and was: the flags describe a
# RECTANGLE, while the sub-grids launched above may deliberately cover only part
# of it. A run scoped to "T4 at k in {11,15,21,31}, coverages 10/30/50" handed a
# rectangle of Lambda,T4 x 7 k-values x 4 coverages had the "aggregation" begin
# computing the 96-cell difference — silently re-expanding a scope that had been
# narrowed on purpose because those cells run hours each.
echo "--- final aggregation (serial, reads checkpoints, computes nothing) ---"
if ! julia --project="${REPO_ROOT}" "${REPO_ROOT}/benchmarking/ont_k_sweep.jl" \
    --aggregate-only --output-dir "${OUTPUT_DIR}" 2>&1 | tee "${LOG_DIR}/final-aggregate.log" | tail -60; then
    echo "final aggregation reported non-ok cells (see ${LOG_DIR}/final-aggregate.log)"
    failed=$((failed + 1))
fi

if [ "${failed}" -gt 0 ]; then
    echo "=== done WITH ${failed} FAILURE(S) — the TSVs may be incomplete ==="
    exit 1
fi
echo "=== done ==="
