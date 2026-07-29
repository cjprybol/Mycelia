# Track A harvest figures contract tests.
#
# benchmarking/ is not otherwise exercised by Pkg.test(), so these pin the pure helpers
# in benchmarking/track_a_harvest_figures.jl inside the sweep — the same pattern as
# gap_calibration_fitters_test.jl and track_a_baseline_benchmark_test.jl.
#
# Every assertion below corresponds to a defect that reached review on PR #438. They
# share one shape: a value that did not vary with its input. A caption that always
# warns, a scope caveat that never fires, a reduction seeded with `init` that always
# returns the seed. Inspection missed all of them; only running the same code against
# two different inputs exposed them, which is what these tests do.
#
# Run:
#   julia --project=. test/4_assembly/track_a_harvest_figures_test.jl

# Module-wrapped: runtests.jl includes every test file into one shared `Main`, and this
# script declares CV_THRESHOLD / TECHNOLOGY_ORDER / REFERENCE_LENGTHS at top level.
# A bare include risks the same "invalid redefinition of constant" failure that took
# down the whole suite when track_a_baseline_benchmark_test.jl was first added.
module TrackAHarvestFiguresTest

import Test
import DataFrames

include(joinpath(@__DIR__, "..", "..", "benchmarking", "track_a_harvest_figures.jl"))

# One row of the harness schema, overridable per field.
function cell(; organism = "Lambda", technology = "illumina", coverage = 30,
        seed = 42, decoder_arm = "kmer", k = 31, NGA50 = 1000.0,
        genome_fraction = 99.0, n_contigs = 10, status = "ok")
    return (organism = organism, technology = technology, coverage = coverage,
        seed = seed, decoder_arm = decoder_arm, k = k, NGA50 = NGA50,
        genome_fraction = genome_fraction, n_contigs = n_contigs, status = status)
end

Test.@testset "track A harvest figures helpers" begin
    Test.@testset "ok_cells is missing-safe" begin
        # `df.status .== "ok"` yields Union{Missing,Bool} when any status is missing,
        # and DataFrames throws on indexing with that mask — so ONE missing status
        # aborted the figure instead of dropping one row.
        df = DataFrames.DataFrame([cell(), cell(seed = 123, status = "error")])
        allowmissing = DataFrames.allowmissing(df, :status)
        allowmissing[2, :status] = missing
        Test.@test DataFrames.nrow(ok_cells(allowmissing)) == 1
        # A frame with no status column is passed through, not emptied.
        Test.@test DataFrames.nrow(ok_cells(DataFrames.select(df, DataFrames.Not(:status)))) == 2
    end

    Test.@testset "CV table groups by k and keeps non-computable groups" begin
        # Same organism/tech/coverage/arm at two k. Pooling them turns a between-SEED
        # CV into a between-K CV, which measured ~39x inflation and flipped the
        # pre-registration verdict under an axis still labelled "n = 3 seeds".
        rows = [cell(seed = s, k = kk, NGA50 = kk == 31 ? 1000.0 : 5000.0)
                for kk in (31, 19) for s in (42, 123, 456)]
        table = coefficient_of_variation_table(DataFrames.DataFrame(rows))
        Test.@test DataFrames.nrow(table) == 2          # one group per k, not one pooled
        Test.@test Set(table.k) == Set([31, 19])
        Test.@test all(table.n_seeds .== 3)

        # An all-zero NGA50 group is a RESULT (every ONT cell at 10x/30x), so it is
        # retained and marked, never silently dropped.
        zeros = DataFrames.DataFrame([cell(seed = s, NGA50 = 0.0) for s in (42, 123, 456)])
        zt = coefficient_of_variation_table(zeros)
        Test.@test DataFrames.nrow(zt) == 1
        Test.@test !only(zt.computable)
    end

    Test.@testset "an all-failed harvest names its cause" begin
        # Previously this produced a 0x0 frame with no :decoder_arm and threw a bare
        # ArgumentError from deep inside plotting.
        failed = DataFrames.DataFrame([cell(seed = s, status = "error") for s in (42, 123)])
        thrown = nothing
        try
            coefficient_of_variation_table(failed)
        catch e
            thrown = e
        end
        Test.@test thrown isa ErrorException
        Test.@test occursin("status == \"ok\"", thrown.msg)
    end

    Test.@testset "figure 1 refuses an empty kmer arm by name" begin
        # The all-rows-dropped guard above fires only when EVERY row goes; figure 1 then
        # narrows to the kmer arm, and THAT can be empty while the table is not. It threw
        # "reducing over an empty collection ... consider supplying `init`" — advice which,
        # taken literally, is what produced the caption bug this suite also pins.
        qualmer_only = DataFrames.DataFrame([cell(seed = s, decoder_arm = "qualmer")
                                             for s in (42, 123, 456)])
        table = coefficient_of_variation_table(qualmer_only)
        Test.@test DataFrames.nrow(table) == 1          # the table itself is fine
        thrown = nothing
        try
            mktempdir(dir -> figure_cv_vs_threshold(table, dir))
        catch e
            thrown = e
        end
        Test.@test thrown isa ErrorException
        Test.@test occursin("no kmer-arm rows", thrown.msg)
    end

    Test.@testset "describe_k reports the k actually present" begin
        one_k = DataFrames.DataFrame([cell(k = 31)])
        Test.@test occursin("k = 31", describe_k(one_k))
        mixed = DataFrames.DataFrame([cell(k = 31), cell(k = 19, seed = 123)])
        # A mixed table must announce itself rather than print one k over another's data.
        Test.@test occursin("MIXED", describe_k(mixed))
        Test.@test describe_k(DataFrames.DataFrame([(a = 1,)])) == ""
    end

    Test.@testset "require_columns names what is absent" begin
        df = DataFrames.DataFrame([cell()])
        Test.@test require_columns(df, "t.tsv", [:organism, :NGA50]) === nothing
        thrown = nothing
        try
            require_columns(df, "t.tsv", [:organism, :quast_nga50])
        catch e
            thrown = e
        end
        Test.@test thrown isa ErrorException
        Test.@test occursin("quast_nga50", thrown.msg)
    end

    Test.@testset "technology_color tolerates an unknown technology" begin
        # A KeyError here aborted rendering partway through a figure.
        Test.@test technology_color("illumina") == :dodgerblue3
        Test.@test technology_color("nanopore") == :grey40
    end
end

end  # module TrackAHarvestFiguresTest
