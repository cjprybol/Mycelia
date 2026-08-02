# Track A baseline benchmark contract tests.
#
# These pin the pure helpers in benchmarking/track_a_baseline_benchmark.jl inside the
# Pkg.test() sweep, since benchmarking/ is not otherwise run in CI — the same pattern
# as gap_calibration_fitters_test.jl for benchmarking/calibration_metrics.jl.
#
# Each testset below corresponds to a defect that reached review in PR #438. All three
# were reachable only through paths the author's test plan could not exercise: a
# parse check cannot see a runtime keyword error, and `--smoke` runs exactly one
# deliberately-successful cell, so neither the error path nor a resume-from-existing-
# tree was ever executed.
#
# Run:
#   julia --project=. test/4_assembly/track_a_baseline_benchmark_test.jl

# Wrapped in a module. runtests.jl includes every test file into one shared `Main`
# (include_all_tests), and the driver defines top-level consts — COVERAGES, SEEDS,
# ORGANISMS — that collide with benchmarking/rhizomorph_stage2_toy_benchmark.jl, which
# another test already includes there. A bare include therefore died with
# "invalid redefinition of constant Main.COVERAGES" and took the whole suite with it.
# A module gives the driver its own namespace; gap_calibration_fitters_test.jl can
# include directly only because calibration_metrics.jl declares no such consts.
module TrackABaselineBenchmarkTest

import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "track_a_baseline_benchmark.jl"))

Test.@testset "track A baseline benchmark helpers" begin

    # --- A failing cell must produce a row, not throw from inside the handler ----
    # cell_row gained a required `peak_rss_method` keyword; the catch-block call was
    # not updated, so the error handler raised UndefKeywordError and the first failing
    # cell aborted the whole matrix — the opposite of the degradation it implements.
    Test.@testset "error path builds a complete row" begin
        row = error_row("Lambda", "NC_001416", "ont", 10, 42, "kmer")
        Test.@test row.status == "error"
        Test.@test row.n_contigs == 0
        Test.@test row.peak_rss_method == "unknown"
        # Rectangular aggregate: an error row must carry every schema column, or the
        # DataFrame built from mixed rows is ragged.
        Test.@test Set(keys(row)) == Set(ROW_KEYS)
    end

    # --- Default-k cell ids must not move ---------------------------------------
    # Interpolating k unconditionally renamed the default k=31 cells, orphaning the
    # 288 Lovelace + 144 LRC checkpoints written before --k existed and silently
    # converting the wrappers' "re-submit to resume" into a full recompute.
    Test.@testset "cell_id is backward compatible at default k" begin
        Test.@test cell_id_for("Lambda", "ont", 10, 42, "kmer"; k = 31) ==
                   "Lambda__ont__10x__seed42__kmer"
        # A non-default k must still be distinguishable, or two sweeps collide on one
        # checkpoint and the second silently republishes the first's result.
        Test.@test cell_id_for("Lambda", "ont", 10, 42, "kmer"; k = 19) ==
                   "Lambda__ont__10x__seed42__kmer__k19"
        Test.@test cell_id_for("Lambda", "ont", 10, 42, "kmer"; k = 19) !=
                   cell_id_for("Lambda", "ont", 10, 42, "kmer"; k = 31)
    end

    # --- Resume tolerance is provenance-only ------------------------------------
    Test.@testset "canonical tolerates provenance keys but not measurements" begin
        complete = Dict("organism" => "Lambda", "accession" => "NC_001416",
            "technology" => "ont", "coverage" => 10, "seed" => 42,
            "decoder_arm" => "kmer", "k" => 31, "n_reads" => 100, "n_contigs" => 1,
            "NGA50" => 0.0, "misassemblies" => 0.0, "genome_fraction" => 0.0,
            "duplication_ratio" => 0.0, "largest_contig" => 0, "wall_seconds" => 1.0,
            "peak_rss_bytes" => 0, "status" => "ok")

        # A pre-peak_rss_method checkpoint still reloads, tagged so it can never be
        # mistaken for a measurement.
        row = canonical(complete)
        Test.@test row.peak_rss_method == "unknown"
        Test.@test row.rss_baseline_bytes == -1
        Test.@test Set(keys(row)) == Set(ROW_KEYS)

        # A missing MEASUREMENT must fail loudly. Zero-filling is unsafe here
        # specifically: NGA50 = 0.0 is the genuine value for every ONT cell at
        # 10x/30x, so a zero-filled corrupt checkpoint is byte-identical to a real
        # result and would enter the CV that gates the pre-registration.
        for key in ("NGA50", "genome_fraction", "status", "coverage")
            truncated = copy(complete)
            delete!(truncated, key)
            Test.@test_throws ErrorException canonical(truncated)
        end
    end

    # --- Memory measurement labels itself ---------------------------------------
    Test.@testset "timed_with_peak_rss labels its measurement" begin
        result = timed_with_peak_rss(() -> sum(rand(200_000)))
        Test.@test result.method in ("sampled", "sampled-degraded", "highwater-delta")
        Test.@test result.wall_seconds >= 0
        Test.@test result.value isa Real
        if Sys.islinux() && Threads.nthreads() > 1
            # CI runs a --threads=4 ubuntu job, so this is the production path there.
            Test.@test result.method in ("sampled", "sampled-degraded")
            # Nonzero is the whole point: every kmer cell read exactly 0 before the fix.
            Test.@test result.peak_rss_bytes > 0
        end
        # A throwing workload must propagate its OWN error and must not hang: the
        # sampler is reaped in the finally, so telemetry can never discard science.
        Test.@test_throws ErrorException timed_with_peak_rss(() -> error("boom"))
    end

    # --- A mid-run sampler failure must not be labelled a clean measurement ------
    # The `try` used to wrap the WHOLE polling loop, so the first probe exception
    # ended sampling permanently for that cell. If any poll had already succeeded,
    # `reads > 0` labelled the frozen, silently-truncated peak "sampled" — a
    # plausible number presented as a measurement, which is precisely the failure the
    # read count was introduced to make visible.
    #
    # `probe` is injected because there is no portable way to make a real /proc read
    # fail on the third poll, and `current_rss_bytes` returns `nothing` on macOS.
    Test.@testset "a truncated sampler downgrades its own label" begin
        # Three good reads (the first is consumed by the baseline), then every read
        # throws forever.
        calls = Threads.Atomic{Int}(0)
        flaky = () -> begin
            n = Threads.atomic_add!(calls, 1) + 1
            n <= 3 && return 1_000_000 * n
            error("simulated /proc read failure")
        end
        result = timed_with_peak_rss(() -> (sleep(0.3); 42);
            interval_seconds = 0.005, probe = flaky)
        # Telemetry never discards science, whichever path was taken.
        Test.@test result.value == 42

        if Threads.nthreads() > 1
            Test.@test result.method == "sampled-partial"
            # The peak is the last SUCCESSFUL observation, and it is reported as a
            # lower bound rather than as "sampled".
            Test.@test result.peak_rss_bytes == 3_000_000
            Test.@test result.rss_baseline_bytes == 1_000_000
            # Sampling CONTINUED past the first failure instead of exiting the loop:
            # 3 successes + 1 failure would be the whole story under the old code.
            Test.@test calls[] > 4
        else
            # Single-threaded, so the sampler is not used at all and the row is
            # honestly labelled with the known-broken fallback semantics.
            Test.@test result.method == "highwater-delta"
        end

        # A probe that never fails is still a clean "sampled" measurement — the
        # downgrade must be caused by the failures, not by injecting a probe.
        rising = Threads.Atomic{Int}(0)
        healthy = () -> 1_000_000 * (Threads.atomic_add!(rising, 1) + 1)
        clean = timed_with_peak_rss(() -> (sleep(0.05); 7);
            interval_seconds = 0.005, probe = healthy)
        Test.@test clean.value == 7
        Test.@test clean.method ==
                   (Threads.nthreads() > 1 ? "sampled" : "highwater-delta")

        # A probe that fails on EVERY call, including the baseline, cannot sample at
        # all and must not claim to have.
        always_bad = () -> error("simulated /proc read failure")
        dead = timed_with_peak_rss(() -> 3; interval_seconds = 0.005,
            probe = always_bad)
        Test.@test dead.value == 3
        Test.@test dead.method == "highwater-delta"
        Test.@test dead.rss_baseline_bytes == -1
    end
end

end  # module TrackABaselineBenchmarkTest
