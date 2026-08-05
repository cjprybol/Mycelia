# Concurrency-cap tests for benchmarking/run_ont_k_sweep_shards.sh.
#
# Run:
#   julia benchmarking/run_ont_k_sweep_shards_test.jl
#
# No Mycelia dependency and no real assembly: the driver is shell, and the thing
# under test is how many julia processes it holds open at once. A stub `julia`
# placed earlier on PATH than the real one records how many copies of itself are
# live, so every assertion here is about OBSERVED concurrency rather than about
# the text of the script. A test that grepped for `MAX_PARALLEL` would pass
# against a driver that computes a cap and then ignores it.
#
# The uncapped driver fails the first testset: an 8-cell grid launches 8
# concurrent shards regardless of what the host can hold.

import Test

const DRIVER = joinpath(@__DIR__, "run_ont_k_sweep_shards.sh")

# Stub `julia`. Registers itself in a live/ directory, folds the resulting count
# into a peak file under a mkdir mutex, holds for a beat so shards genuinely
# overlap, then deregisters.
#
# Two invocation shapes are NOT shards and are excluded from the failure
# injection below: the serial pre-warm (identified by --seeds, which no shard
# passes) and the final --aggregate-only pass. Failing the pre-warm would abort
# the driver under `set -e` before any shard launched.
const STUB_JULIA = raw"""
#!/usr/bin/env bash
set -u
TRACK="${MYCELIA_TEST_TRACK_DIR}"

is_shard=1
for arg in "$@"; do
    case "$arg" in
        --seeds|--aggregate-only) is_shard=0 ;;
    esac
done

printf 'threads=%s shard=%s args=%s\n' "${JULIA_NUM_THREADS:-unset}" "$is_shard" "$*" \
    >> "${TRACK}/calls.log"

mkdir -p "${TRACK}/live"
mkdir "${TRACK}/live/$$" 2>/dev/null || true

n=$(ls "${TRACK}/live" | wc -l | tr -d ' ')
while ! mkdir "${TRACK}/peak.lock" 2>/dev/null; do :; done
peak=$(cat "${TRACK}/peak" 2>/dev/null || echo 0)
if [ "$n" -gt "$peak" ]; then printf '%s' "$n" > "${TRACK}/peak"; fi
rmdir "${TRACK}/peak.lock"

sleep "${MYCELIA_TEST_SHARD_SECONDS:-0.4}"
rmdir "${TRACK}/live/$$" 2>/dev/null || true

if [ "$is_shard" = 1 ]; then exit "${MYCELIA_TEST_SHARD_EXIT:-0}"; fi
exit 0
"""

"""
    stage_sandbox() -> (root, track)

Stage a throwaway repo holding a copy of the driver, a stub target script, and a
stub `julia` on a bin/ directory meant to be prepended to PATH.
"""
function stage_sandbox()
    root = mktempdir()
    mkpath(joinpath(root, "benchmarking"))
    mkpath(joinpath(root, "bin"))
    track = joinpath(root, "track")
    mkpath(track)

    cp(DRIVER, joinpath(root, "benchmarking", "run_ont_k_sweep_shards.sh"))
    write(joinpath(root, "benchmarking", "ont_k_sweep.jl"), "# stub target — never executed\n")

    stub = joinpath(root, "bin", "julia")
    write(stub, STUB_JULIA)
    chmod(stub, 0o755)

    return (root, track)
end

"""
    run_driver(root, track; env...) -> (output, exitcode)

Run the staged driver with the stub `julia` first on PATH. Returns combined
stdout/stderr and the exit code; a non-zero exit is a result to assert on, not
an error.
"""
function run_driver(root, track; env::Dict{String, String} = Dict{String, String}())
    base = Dict(
        "PATH" => joinpath(root, "bin") * ":" * ENV["PATH"],
        "MYCELIA_TEST_TRACK_DIR" => track,
        "OUTPUT_DIR" => joinpath(root, "out"),
        # Poll fast so the test is seconds rather than minutes. The production
        # default is deliberately slower.
        "THROTTLE_POLL_SECONDS" => "0.05",
        # A 2x2x2 grid: 8 shards, small enough to be quick and large enough that
        # an uncapped run is unmistakable.
        "KS" => "11 15",
        "COVERAGES" => "10 30",
        "TECHNOLOGIES" => "illumina ont",
        "ORGANISMS" => "Lambda"
    )
    merge!(base, env)
    # Inherit only what bash needs; a stray SLURM_* from the caller's shell would
    # otherwise silently change the cap the driver computes.
    for key in ("HOME", "SHELL", "TMPDIR", "LANG")
        haskey(ENV, key) && (base[key] = ENV[key])
    end

    script = joinpath(root, "benchmarking", "run_ont_k_sweep_shards.sh")
    cmd = setenv(`bash $script`, base)
    out = IOBuffer()
    proc = run(pipeline(ignorestatus(cmd), stdout = out, stderr = out))
    return (String(take!(out)), proc.exitcode)
end

function peak_concurrency(track)
    isfile(joinpath(track, "peak")) ?
    parse(Int, strip(read(joinpath(track, "peak"), String))) : 0
end

function shard_calls(track)
    isfile(joinpath(track, "calls.log")) ?
    filter(l -> occursin("shard=1", l), split(strip(read(joinpath(track, "calls.log"), String)), '\n')) :
    String[]
end

Test.@testset "run_ont_k_sweep_shards.sh concurrency cap" begin
    Test.@testset "an explicit MAX_PARALLEL bounds live shards" begin
        # THE REGRESSION GUARD. Uncapped, all 8 shards of this grid run at once
        # and peak == 8; the cap must hold it to 2 no matter how big the grid is.
        root, track = stage_sandbox()
        out, code = run_driver(root, track; env = Dict("MAX_PARALLEL" => "2"))
        Test.@test length(shard_calls(track)) == 8      # whole grid still covered
        Test.@test peak_concurrency(track) <= 2         # ... but never 8 at once
        Test.@test code == 0
        Test.@test occursin("max parallel", lowercase(out))
    end

    Test.@testset "a cap of 1 serializes the grid" begin
        root, track = stage_sandbox()
        run_driver(root, track; env = Dict("MAX_PARALLEL" => "1"))
        Test.@test length(shard_calls(track)) == 8
        Test.@test peak_concurrency(track) == 1
    end

    Test.@testset "the derived cap is bounded by core count" begin
        # With no MAX_PARALLEL the driver must derive one from the host rather
        # than fanning out to grid size. Bounded above by cores; never below 1.
        root, track = stage_sandbox()
        out, _ = run_driver(root, track)
        m = match(r"max parallel:\s*(\d+)"i, out)
        Test.@test m !== nothing
        derived = parse(Int, m.captures[1])
        Test.@test 1 <= derived <= Sys.CPU_THREADS
        Test.@test peak_concurrency(track) <= derived
    end

    Test.@testset "a SLURM allocation caps below the node" begin
        # Inside `--mem=8G` the driver must size from the ALLOCATION. Reading
        # node memory here is what gets a sub-node job OOM-killed.
        root, track = stage_sandbox()
        out,
        _ = run_driver(root,
            track;
            env = Dict(
                "SLURM_MEM_PER_NODE" => "8192",   # 8 GB
                "SLURM_CPUS_ON_NODE" => "16",
                "GB_PER_SHARD" => "4",
                "MEM_FRACTION_PCT" => "100"
            ))
        m = match(r"max parallel:\s*(\d+)"i, out)
        Test.@test m !== nothing
        Test.@test parse(Int, m.captures[1]) == 2       # 8 GB / 4 GB, not 16 cores
    end

    Test.@testset "--mem=0 means the whole node, not zero memory" begin
        # SLURM reports an exclusive allocation as 0. Read literally that yields
        # a 0 GB budget and a cap of 1 — serializing the grid on the exact
        # command the LBNL recipe tells you to run.
        root, track = stage_sandbox()
        out,
        _ = run_driver(root, track; env = Dict(
            "SLURM_MEM_PER_NODE" => "0",
            "SLURM_CPUS_ON_NODE" => "16"
        ))
        m = match(r"max parallel:\s*(\d+)"i, out)
        Test.@test m !== nothing
        Test.@test parse(Int, m.captures[1]) > 1
    end

    Test.@testset "threads per shard stay within the core allocation" begin
        # JULIA_NUM_THREADS=auto is inherited by every backgrounded shard, so
        # capping processes alone still oversubscribes the CPU: 5 shards on a
        # 56-core node is 5 x 56 threads. Processes x threads must fit.
        root, track = stage_sandbox()
        run_driver(root,
            track;
            env = Dict(
                "MAX_PARALLEL" => "4",
                "SLURM_CPUS_ON_NODE" => "16",
                "JULIA_NUM_THREADS" => "auto"
            ))
        calls = shard_calls(track)
        Test.@test !isempty(calls)
        threads = [parse(Int, match(r"threads=(\S+)", c).captures[1]) for c in calls]
        Test.@test all(t -> t >= 1, threads)
        Test.@test maximum(threads) * 4 <= 16
    end

    Test.@testset "an explicit thread count is respected" begin
        root, track = stage_sandbox()
        run_driver(root, track; env = Dict(
            "MAX_PARALLEL" => "2",
            "JULIA_NUM_THREADS" => "3"
        ))
        calls = shard_calls(track)
        Test.@test !isempty(calls)
        Test.@test all(c -> occursin("threads=3", c), calls)
    end

    Test.@testset "throttling does not break failure reporting" begin
        # `throttle` calls `jobs`, which reaps finished children. The driver's
        # existing `wait "$pid"` failure detection and its exit 1 must survive
        # that — the cap sits alongside those, it does not replace them.
        root, track = stage_sandbox()
        out,
        code = run_driver(root, track; env = Dict(
            "MAX_PARALLEL" => "2",
            "MYCELIA_TEST_SHARD_EXIT" => "1"
        ))
        Test.@test code == 1
        Test.@test occursin("FAILURE", out)
        Test.@test occursin("exited non-zero", out)
    end
end
