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
        # Assert the MESSAGE, not just the type. `@test_throws ErrorException` passes on
        # any unrelated ErrorException — including one from a future refactor that
        # happens to throw for a different reason — so it would keep passing while
        # silently no longer testing missing-key handling.
        for key in ("NGA50", "genome_fraction", "status", "coverage")
            truncated = copy(complete)
            delete!(truncated, key)
            thrown = nothing
            try
                canonical(truncated)
            catch e
                thrown = e
            end
            Test.@test thrown isa ErrorException
            Test.@test occursin("missing required key $(key)", thrown.msg)
        end

        # A key that IS provenance-only must be defaultable purely by being listed,
        # which is the remediation the error message advertises.
        Test.@test Set(OPTIONAL_KEYS) == Set(keys(OPTIONAL_KEY_DEFAULTS))
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

    Test.@testset "power analysis counts only measurements" begin
        # Two ways a row carries zeros that are not measurements, and `status` sees
        # only the first. On a QUAST exception run_cell substitutes empty_metrics()
        # and then derives status from n_contigs ALONE, so a QUAST failure on a
        # NON-EMPTY assembly is written as status="ok" with NGA50 = 0. Pooling that
        # into a variance estimate is the same defect the status filter was added to
        # close, reached through a different door. Mirrors ok_cells() in
        # track_a_harvest_figures.jl — the two consumers must not disagree about
        # which rows are measurements.
        row(; nga,
            status,
            n_contigs,
            largest_contig) = (
            organism = "Lambda", technology = "ont", coverage = 30,
            decoder_arm = "kmer", k = 31, NGA50 = nga, status = status,
            n_contigs = n_contigs, largest_contig = largest_contig)
        measured = [
            row(nga = 1000.0, status = "ok", n_contigs = 12,
                largest_contig = 1000),
            row(nga = 1010.0, status = "ok", n_contigs = 12,
                largest_contig = 1010),
            row(nga = 990.0, status = "ok", n_contigs = 12,
                largest_contig = 990)]
        # status="ok", non-empty assembly, but QUAST never scored it.
        unscored = row(nga = 0.0, status = "ok", n_contigs = 1200,
            largest_contig = 0)

        mktempdir() do dir
            clean = write_power_analysis(dir, DataFrames.DataFrame(measured))
            mixed = write_power_analysis(dir,
                DataFrames.DataFrame(vcat(measured, [unscored])))
            # The unscored row must not enter the group at all...
            Test.@test only(clean.n) == 3
            Test.@test only(mixed.n) == 3
            # ...and therefore must not move the variance it would otherwise inflate.
            Test.@test only(mixed.cv_nga50) == only(clean.cv_nga50)
        end
    end
    Test.@testset "aggregate publish is atomic and shard-safe" begin
        # write_aggregate runs once per CELL, and sharding into one results tree is the
        # documented workflow (the header advertises "an HPC array job can split the
        # matrix and share one results tree"; the LRC wrapper shows a shard invocation
        # with no --output-dir). A FIXED temp path therefore let two shards open the
        # same inode and both rename it, publishing interleaved content as a complete
        # table — worse than the truncation it replaced, because a short table is
        # detectable and a corrupt one is not.
        mktempdir() do dir
            target = joinpath(dir, "results.tsv")
            temps = String[]
            for i in 1:5
                _publish_atomically(target) do tmp
                    push!(temps, tmp)
                    write(tmp, "payload-$(i)")
                end
            end
            # The load-bearing property: a distinct temp per call. With a fixed name
            # this is 1, and concurrent shards collide.
            Test.@test length(unique(temps)) == 5
            Test.@test read(target, String) == "payload-5"

            # A failed write leaves no temp behind and does not damage the target.
            threw = try
                _publish_atomically(target) do tmp
                    write(tmp, "half-written")
                    error("simulated writer failure")
                end
                nothing
            catch e
                e
            end
            Test.@test threw isa ErrorException
            Test.@test occursin("simulated writer failure", threw.msg)
            Test.@test read(target, String) == "payload-5"
            Test.@test readdir(dir) == ["results.tsv"]
        end
    end
end

end  # module TrackABaselineBenchmarkTest
