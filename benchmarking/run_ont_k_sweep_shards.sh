#!/usr/bin/env bash
# Parallel driver for benchmarking/ont_k_sweep.jl.
#
# The 96-cell grid is embarrassingly parallel — every cell is independent and
# writes its own checkpoint directory — but run serially it is dominated by the
# ONT/100x cells, where a single assembly takes tens of minutes. This partitions
# the grid by (organism, k, coverage, technology) and runs the shards
# concurrently, so wall-clock is set by the slowest single shard rather than by
# the sum of all of them. With the defaults that is 32 shards of 3 cells each
# (3 seeds), of which only MAX_PARALLEL run at a time -- see the concurrency cap
# below, which sizes itself from the host's memory and cores.
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
# them from queueing behind their Illumina siblings.
#
# Shard COUNT is no longer the same thing as concurrency: the cap below decides
# how many run at once, so a bigger grid now means a longer queue rather than a
# heavier host. This comment used to advise laptop users to set TECHNOLOGIES
# manually to halve the fan-out, which is exactly the kind of advisory that does
# not fire -- on 2026-08-04 it did not.
TECHNOLOGIES="${TECHNOLOGIES:-illumina ont}"
# Organisms default to Lambda, matching the sweep script's own default. T4 is
# opt-in because it is 3.5x larger and its ONT high-coverage cells dominate the
# grid's cost: ORGANISMS="Lambda T4" ./run_ont_k_sweep_shards.sh
ORGANISMS="${ORGANISMS:-Lambda}"

# --- concurrency cap ---------------------------------------------------------
# The grid is embarrassingly parallel, but "no data dependencies" is not "any
# host can hold N copies of the working set". Each shard is a full Julia process
# holding an assembly graph; measured at 2.4-5.9 GB RSS on the 2026-08-04 run
# regardless of genome size, because the footprint is runtime + graph, not input.
#
# An uncapped `for ... &` therefore scales memory with GRID SIZE rather than host
# capacity, and the grid keeps growing: 16 shards, then 32 when technology
# sharding was added, and 64 with ORGANISMS="Lambda T4". On 2026-08-04 sixteen
# shards on an 18-core laptop drove swap to 95% and completed zero cells in ten
# hours -- slower than serial would have been.
#
# The cap is derived from usable memory at a conservative GB-per-shard, then
# bounded by core count. Override explicitly with MAX_PARALLEL=N.
#
# AVAILABLE memory, never total. Sizing from total RAM is what would make this
# cap useless on the host that needs it: the incident laptop reports 128 GB and
# 18 cores, so a total-derived cap yields 18 -- MORE than the 16 shards that
# thrashed it. That machine is never idle; at incident time it was also running
# dozens of agent processes, Docker, and an fseventsd pinned at 93% CPU by the
# checkpoint writes. Plan from what is free and leave headroom.
#
# The budget is resolved most-authoritative source first, and a SLURM allocation
# outranks node memory. Inside a sub-node job (`--mem=64G`, which is this repo's
# own sub-node recipe) /proc/meminfo reports the whole NODE, so a node-derived
# cap launches enough shards to blow the allocation and get the job OOM-killed.
# A cap written to protect a laptop must not be blind on the cluster it routes
# work to.
detect_avail_gb() {
    # 1. SLURM per-node allocation (MB).
    #    ZERO MEANS ALL NODE MEMORY, not zero memory -- `--mem=0` is what this
    #    repo's own recommended Lawrencium recipe passes (`--exclusive --mem=0`).
    #    Reading it literally gave a 0 GB budget and MAX_PARALLEL=1, serializing
    #    the whole grid on the exact command the protocol tells you to run. Fall
    #    through to node detection instead, which is correct under --exclusive.
    if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "${SLURM_MEM_PER_NODE}" -gt 0 ] 2>/dev/null; then
        printf '%d' $(( SLURM_MEM_PER_NODE / 1024 )); return
    fi
    # 2. SLURM per-CPU allocation x allocated CPUs (MB); same zero semantics.
    if [ -n "${SLURM_MEM_PER_CPU:-}" ] && [ "${SLURM_MEM_PER_CPU}" -gt 0 ] 2>/dev/null \
       && [ -n "${SLURM_CPUS_ON_NODE:-}" ]; then
        printf '%d' $(( SLURM_MEM_PER_CPU * SLURM_CPUS_ON_NODE / 1024 )); return
    fi
    # 3. cgroup v2 / v1 limit, when finite
    for f in /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory/memory.limit_in_bytes; do
        if [ -r "$f" ]; then
            v=$(cat "$f" 2>/dev/null)
            case "$v" in
                max|""|*[!0-9]*) : ;;
                *) if [ "$v" -lt 1099511627776 ]; then printf '%d' $(( v / 1073741824 )); return; fi ;;
            esac
        fi
    done
    # 4. node availability
    if [ -r /proc/meminfo ]; then
        awk '/^MemAvailable:/ {printf "%d", $2/1048576; found=1}
             END {if (!found) print 8}' /proc/meminfo
    elif command -v vm_stat >/dev/null 2>&1; then
        vm_stat | awk '
            /page size of/ { for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) ps=$i }
            /^Pages free:/         { gsub(/\./,"",$3); f=$3 }
            /^Pages inactive:/     { gsub(/\./,"",$3); n=$3 }
            /^Pages speculative:/  { gsub(/\./,"",$3); sp=$3 }
            END { if (ps=="") ps=4096; printf "%d", (f+n+sp)*ps/1073741824 }'
    else
        printf '8'
    fi
}

detect_cores() {
    # SLURM tells us how many CPUs we actually hold; nproc reports the node.
    if [ -n "${SLURM_CPUS_ON_NODE:-}" ]; then printf '%s' "$SLURM_CPUS_ON_NODE"
    elif command -v nproc >/dev/null 2>&1; then nproc
    elif command -v sysctl >/dev/null 2>&1; then sysctl -n hw.ncpu 2>/dev/null || printf '4'
    else printf '4'; fi
}

# Fraction of AVAILABLE memory this sweep may claim. The rest is headroom for
# the interactive session, the editor, and whatever else shares the host.
MEM_FRACTION_PCT="${MEM_FRACTION_PCT:-60}"
GB_PER_SHARD="${GB_PER_SHARD:-6}"
_cores_for_threads="$(detect_cores)"
if [ -z "${MAX_PARALLEL:-}" ]; then
    _budget_gb=$(( $(detect_avail_gb) * MEM_FRACTION_PCT / 100 ))
    MAX_PARALLEL=$(( _budget_gb / GB_PER_SHARD ))
    # Clamps are written as `if` rather than `[ ... ] && x=y` on purpose. Under
    # `set -e` the bare AND-list form survives at top level but aborts the whole
    # script the moment it is moved inside a function -- silently, before any
    # shard launches. detect_avail_gb and detect_cores directly above are
    # already functions, so that refactor is one step away.
    if [ "${MAX_PARALLEL}" -gt "${_cores_for_threads}" ]; then
        MAX_PARALLEL="${_cores_for_threads}"
    fi
    if [ "${MAX_PARALLEL}" -lt 1 ]; then
        MAX_PARALLEL=1
    fi
fi

# Threads PER SHARD. The protocol's sbatch wrapper exports
# JULIA_NUM_THREADS=auto, which each backgrounded shard inherits -- so capping
# PROCESSES alone still oversubscribes the CPU: MAX_PARALLEL=5 on a 56-core node
# becomes 5 x 56 = 280 threads. Divide the allocation instead. Memory and CPU are
# two dimensions of the same cap and fixing only one leaves the other fighting it.
if [ -z "${JULIA_NUM_THREADS:-}" ] || [ "${JULIA_NUM_THREADS}" = auto ]; then
    JULIA_NUM_THREADS=$(( _cores_for_threads / MAX_PARALLEL ))
    if [ "${JULIA_NUM_THREADS}" -lt 1 ]; then
        JULIA_NUM_THREADS=1
    fi
    export JULIA_NUM_THREADS
fi

# Block until fewer than MAX_PARALLEL shards are running.
#
# `jobs -pr` rather than `wait -n` because macOS ships bash 3.2, which has no
# `wait -n`. Calling `jobs` here is safe alongside the `wait "${pid}"` failure
# accounting below: bash retains a terminated child's exit status for an
# explicit `wait` by pid even after `jobs` has reported it, so throttling does
# not cost the driver its ability to detect and exit 1 on failed shards. That is
# asserted directly in run_ont_k_sweep_shards_test.jl.
THROTTLE_POLL_SECONDS="${THROTTLE_POLL_SECONDS:-2}"
throttle() {
    while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "${MAX_PARALLEL}" ]; do
        sleep "${THROTTLE_POLL_SECONDS}"
    done
}

mkdir -p "${LOG_DIR}"

echo "=== ONT k-sweep parallel driver ==="
echo "repo:       ${REPO_ROOT}"
echo "output dir: ${OUTPUT_DIR}"
echo "log dir:    ${LOG_DIR}"
echo "max parallel: ${MAX_PARALLEL} shards x ${JULIA_NUM_THREADS} threads (${GB_PER_SHARD} GB/shard, ${MEM_FRACTION_PCT}% of available RAM)"

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
                # Throttle BEFORE launching, not after: calling it after admits
                # one shard over the cap on every iteration, which is exactly the
                # margin that matters when the constrained resource is memory.
                throttle
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
