# Unit test for the metric-definition consistency guard (bead td-9p91).
#
# The bead's requirement is specific: "any aggregation spanning more than one src
# value RAISES rather than silently proceeding. Opt-out must be explicit (e.g.
# allow_mixed_src=true) and should emit a loud warning naming the src values being
# mixed."
#
# So the central assertion is not "the function exists" — it is that a mixed-src
# aggregation ACTUALLY RAISES, that the message names the mixed values, and that
# the documented opt-out warns rather than passing quietly. A guard without a test
# proving it fires is not a guard.
#
# Dependency-free apart from DataFrames (a direct Mycelia dep that loads without
# importing Mycelia), so the real column type the harness writes is exercised.

import Test
import DataFrames
import Logging

include(joinpath(@__DIR__, "..", "..", "benchmarking", "metric_source_guard.jl"))

function _msg_throws_with(f, needles::Vector{String})
    try
        f()
    catch e
        msg = sprint(showerror, e)
        for needle in needles
            Test.@test occursin(needle, msg)
        end
        return msg
    end
    Test.@test false  # expected a throw, got none
    return ""
end

# Capture @warn output so "warns loudly" is asserted, not assumed.
function _msg_capture_logs(f)
    logger = Test.TestLogger()
    result = Logging.with_logger(logger) do
        f()
    end
    return result, logger.logs
end

# A results table shaped like the RGV sweep's, parameterized by metric definition.
function _msg_frame(sources, thresholds; nga50 = nothing)
    n = length(sources)
    return DataFrames.DataFrame(
        arm = [iseven(i) ? "naive" : "iterative" for i in 1:n],
        seed = fill(42, n),
        quast_nga50 = nga50 === nothing ? collect(1.0:1.0:n) : nga50,
        metric_source = collect(sources),
        quast_min_contig = collect(thresholds)
    )
end

Test.@testset "metric-definition guard (td-9p91)" begin
    Test.@testset "single definition passes and reports the summary" begin
        df = _msg_frame(["quast", "quast", "quast"], [4_850, 4_850, 4_850])
        summary = assert_single_metric_definition(df; context = "H1 NGA50 aggregate")
        Test.@test Dict(summary) == Dict(
            :metric_source => ["quast"], :quast_min_contig => ["4850"])
    end

    Test.@testset "MIXED metric_source RAISES and names the values" begin
        # The exact situation observed on Lawrencium job 24247925: T4 rows carry
        # internal metrics while Lambda rows carry QUAST metrics.
        df = _msg_frame(
            ["quast", "internal:quast-failed", "quast"], [4_850, 4_850, 4_850])
        msg = _msg_throws_with(
            () -> assert_single_metric_definition(df; context = "H1 NGA50 aggregate"),
            ["MixedMetricDefinitionError",
                "H1 NGA50 aggregate",
                "metric_source",
                "\"internal:quast-failed\"",   # both mixed values named
                "\"quast\"",
                "allow_mixed_src"])            # the opt-out is discoverable
        # Not a bare type check: the message must carry the actionable content.
        Test.@test occursin("ALIGNMENT-VALIDATED", msg)
    end

    Test.@testset "MIXED quast_min_contig RAISES too" begin
        # Two rows both labelled "quast" are still incomparable if QUAST filtered
        # them at different thresholds. This is why the sweep records the threshold.
        df = _msg_frame(["quast", "quast"], [4_850, 5_000])
        _msg_throws_with(
            () -> assert_single_metric_definition(df; context = "threshold check"),
            ["quast_min_contig", "\"4850\"", "\"5000\""])
    end

    Test.@testset "missing values participate rather than being skipped" begin
        df = DataFrames.DataFrame(
            metric_source = ["quast", missing],
            quast_min_contig = [4_850, 4_850])
        _msg_throws_with(
            () -> assert_single_metric_definition(df; context = "missing-source check"),
            ["\"missing\"", "\"quast\""])
    end

    Test.@testset "explicit opt-out warns loudly and names the mix" begin
        df = _msg_frame(["quast", "internal:quast-failed"], [4_850, 4_850])
        summary,
        logs = _msg_capture_logs() do
            assert_single_metric_definition(df;
                context = "deliberate override", allow_mixed_src = true)
        end
        Test.@test Dict(summary)[:metric_source] == ["internal:quast-failed", "quast"]
        warns = filter(l -> l.level == Logging.Warn, logs)
        Test.@test length(warns) == 1
        Test.@test occursin("MIXED METRIC DEFINITION", warns[1].message)
        Test.@test occursin("allow_mixed_src=true", warns[1].message)
        Test.@test occursin("NOT validation-grade", warns[1].message)
        mixed = string(Dict(warns[1].kwargs)[:mixed])
        Test.@test occursin("internal:quast-failed", mixed)
        Test.@test occursin("quast", mixed)
    end

    Test.@testset "fail-closed when no definition column is present" begin
        # A table with no provenance column has not been shown to be consistent;
        # it is a table that cannot be checked. Reporting "OK" would be a lie.
        df = DataFrames.DataFrame(arm = ["naive", "iterative"], quast_nga50 = [1.0, 2.0])
        _msg_throws_with(
            () -> assert_single_metric_definition(df; context = "schemaless table"),
            ["no metric-definition column present", "schemaless table", ":metric_source"])
    end

    Test.@testset "custom provenance column (the merge tool's quast_evidence)" begin
        ok = DataFrames.DataFrame(quast_evidence = ["quast:scored", "quast:scored"])
        Test.@test Dict(assert_single_metric_definition(ok;
            columns = (:quast_evidence,), context = "merged matrix")) ==
                   Dict(:quast_evidence => ["quast:scored"])
        mixed = DataFrames.DataFrame(
            quast_evidence = ["quast:scored", "unknown:quast-unscored"])
        _msg_throws_with(
            () -> assert_single_metric_definition(mixed;
                columns = (:quast_evidence,), context = "merged matrix"),
            ["quast_evidence", "\"quast:scored\"", "\"unknown:quast-unscored\""])
    end

    Test.@testset "grouped check: per-group, and every offender reported" begin
        # The realistic shape: a groupby by organism where ONE organism's rows are
        # internally mixed. A whole-table check would flag it too, but the grouped
        # check names WHICH group, which is what an analyst needs.
        df = DataFrames.DataFrame(
            organism = ["Lambda", "Lambda", "T4", "T4", "phi29", "phi29"],
            metric_source = ["quast", "quast",
                "quast", "internal:quast-failed",
                "quast", "internal:quast-failed"],
            quast_min_contig = fill(5_000, 6))
        msg = _msg_throws_with(
            () -> assert_single_metric_definition_per_group(df, (:organism,);
                context = "per-organism NGA50"),
            ["[group T4]", "[group phi29]", "metric_source"])
        # Lambda is clean, so it must NOT be reported as an offender.
        Test.@test !occursin("[group Lambda]", msg)

        clean = DataFrames.DataFrame(
            organism = ["Lambda", "Lambda", "T4", "T4"],
            metric_source = ["quast", "quast", "internal:quast-failed",
                "internal:quast-failed"],
            quast_min_contig = fill(5_000, 4))
        # Each GROUP is internally consistent even though the table as a whole is
        # not — correct for a groupby, since each group becomes one aggregate.
        Test.@test assert_single_metric_definition_per_group(clean, (:organism,)) == 2
        _msg_throws_with(
            () -> assert_single_metric_definition(clean; context = "whole table"),
            ["metric_source"])
    end

    Test.@testset "grouped check validates its group columns" begin
        df = _msg_frame(["quast"], [5_000])
        _msg_throws_with(
            () -> assert_single_metric_definition_per_group(df, (:organism,);
                context = "bad group col"),
            ["group column :organism not present"])
    end

    Test.@testset "guarded_aggregate runs the body only when the guard passes" begin
        ran = Ref(0)
        ok = _msg_frame(["quast", "quast"], [5_000, 5_000])
        value = guarded_aggregate(ok; context = "mean NGA50") do d
            ran[] += 1
            sum(d.quast_nga50)
        end
        Test.@test value == 3.0
        Test.@test ran[] == 1

        mixed = _msg_frame(["quast", "internal:quast-failed"], [5_000, 5_000])
        _msg_throws_with(
            () -> guarded_aggregate(d -> (ran[] += 1), mixed; context = "mean NGA50"),
            ["MixedMetricDefinitionError"])
        Test.@test ran[] == 1  # body never executed on the mixed table
    end

    Test.@testset "metric_definition_summary omits absent columns" begin
        df = DataFrames.DataFrame(metric_source = ["quast", "quast"])
        Test.@test Dict(metric_definition_summary(df)) ==
                   Dict(:metric_source => ["quast"])
    end

    Test.@testset "the sweep records both definition columns" begin
        # The guard is only useful if the harness actually writes what it checks.
        sweep = read(
            joinpath(@__DIR__, "..", "..", "benchmarking",
                "rhizomorph_correction_validation_sweep.jl"), String)
        for col in METRIC_DEFINITION_COLUMNS
            Test.@test occursin("$col = ", sweep)
        end
    end

    Test.@testset "C4/I2: a PARTIALLY absent axis fails closed too" begin
        # The doctrine is per-AXIS. An earlier implementation was per-SET: it raised
        # only when EVERY definition column was missing, so a table with a uniform
        # metric_source and no quast_min_contig — the shape of every CSV produced
        # before this policy — passed and reported a single consistent definition
        # over an axis it never examined.
        partial = DataFrames.DataFrame(metric_source = ["quast", "quast"])
        _msg_throws_with(
            () -> assert_single_metric_definition(partial; context = "pre-policy CSV"),
            ["not checkable", "quast_min_contig", "pre-policy CSV"])

        # The narrower check is still available, but it must be asked for, and it
        # says out loud which axis went unchecked.
        summary, logs = _msg_capture_logs() do
            assert_single_metric_definition(partial;
                context = "pre-policy CSV", require_all_columns = false)
        end
        Test.@test Dict(summary) == Dict(:metric_source => ["quast"])
        warns = filter(l -> l.level == Logging.Warn, logs)
        Test.@test length(warns) == 1
        Test.@test occursin("NOT checked", warns[1].message)
        Test.@test occursin("quast_min_contig", string(Dict(warns[1].kwargs)[:absent]))

        # Naming the axes explicitly records the narrower claim deliberately.
        Test.@test Dict(assert_single_metric_definition(partial;
            columns = (:metric_source,), context = "explicit narrow")) ==
                   Dict(:metric_source => ["quast"])

        # Grouped form fails closed on the same grounds.
        pg = DataFrames.DataFrame(
            organism = ["Lambda", "Lambda"], metric_source = ["quast", "quast"])
        _msg_throws_with(
            () -> assert_single_metric_definition_per_group(pg, (:organism,);
                context = "grouped partial"),
            ["not checkable", "quast_min_contig"])
    end
end
