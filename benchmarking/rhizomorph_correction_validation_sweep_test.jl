# Unit test for the Rhizomorph correction-validation SCALE GUARD.
#
# Verifies the core safety property of the sweep: a toy-scale configuration
# triggers SMOKE-ONLY (verdict NOT allowed), while a realistic large-scale
# configuration would pass (verdict allowed). This test is intentionally
# dependency-free — it includes only the pure guard helpers, so it runs in
# milliseconds without loading Mycelia, hitting the network, or assembling.
#
# Run:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/rhizomorph_correction_validation_sweep_test.jl

import Test
# DataFrames is a direct dep of the Mycelia project and loads WITHOUT importing
# Mycelia, so this test stays dependency-free (no Mycelia, no network, no QUAST)
# while still exercising the production DataFrame push! schema.
import DataFrames

include(joinpath(@__DIR__, "rhizomorph_scale_guard.jl"))
# Dependency-free QUAST parser + per-arm attribution helper (same file the sweep
# includes), so the parse/wiring logic is exercised without running QUAST.
include(joinpath(@__DIR__, "quast_report_parsing.jl"))

Test.@testset "rhizomorph correction validation scale guard" begin
    # --- Toy-scale config triggers SMOKE-ONLY --------------------------------
    # Synthetic 2 kb genome at 10x effective coverage = 20,000 sequenced bases,
    # far below the 1 Mbase floor. This is exactly the local smoke configuration.
    toy_coverage = 10.0
    toy_genome_len = 2_000
    Test.@test scale_metric_bases(toy_coverage, toy_genome_len) == 20_000.0
    Test.@test scale_verdict_allowed(toy_coverage, toy_genome_len) == false

    # A real genome at trivially low coverage is also SMOKE-ONLY: Lambda (48,502 bp)
    # at 1x = 48,502 bases < floor. Guards against "big genome" being mistaken for
    # "enough data".
    Test.@test scale_verdict_allowed(1.0, 48_502) == false

    # --- Large-scale config passes -------------------------------------------
    # Lambda phage (48,502 bp) at 30x = 1,455,060 sequenced bases, above the floor.
    lambda_len = 48_502
    Test.@test scale_metric_bases(30.0, lambda_len) > SCALE_FLOOR_BASES
    Test.@test scale_verdict_allowed(30.0, lambda_len) == true

    # --- Boundary behaviour ---------------------------------------------------
    # Exactly at the floor is allowed (>=), just below is not.
    Test.@test scale_verdict_allowed(1.0, SCALE_FLOOR_BASES) == true
    Test.@test scale_verdict_allowed(1.0, SCALE_FLOOR_BASES - 1) == false

    # --- Custom floor override respected -------------------------------------
    # The harness exposes MYCELIA_RGV_SCALE_FLOOR; the guard must honor an
    # explicit floor argument in both directions.
    Test.@test scale_verdict_allowed(10.0, 2_000; floor = 10_000) == true
    Test.@test scale_verdict_allowed(10.0, 2_000; floor = 100_000) == false
end

Test.@testset "QUAST report parsing + per-arm metric wiring" begin
    mktempdir() do dir
        # --- Fixture: a minimal QUAST report.tsv with the four metrics we wire --
        report_tsv = joinpath(dir, "report.tsv")
        open(report_tsv, "w") do io
            println(io, "Assembly\tnaive_contigs")
            println(io, "# contigs\t12")
            println(io, "Genome fraction (%)\t94.37")
            println(io, "Duplication ratio\t1.02")
            println(io, "NGA50\t8421")
            println(io, "# misassemblies\t3")
            # QUAST prints "-" for a metric it could not compute; must -> missing.
            println(io, "LGA50\t-")
        end

        # --- parse_quast_metric extracts each metric value correctly -----------
        Test.@test parse_quast_metric(report_tsv, "Genome fraction (%)") == 94.37
        Test.@test parse_quast_metric(report_tsv, "NGA50") == 8421.0
        Test.@test parse_quast_metric(report_tsv, "# misassemblies") == 3.0
        Test.@test parse_quast_metric(report_tsv, "Duplication ratio") == 1.02
        # Non-numeric "-" and absent metrics both yield missing.
        Test.@test parse_quast_metric(report_tsv, "LGA50") === missing
        Test.@test parse_quast_metric(report_tsv, "Nonexistent metric") === missing

        # --- quast_metrics_for_report: QUAST-present path populates the row cols -
        q = quast_metrics_for_report(report_tsv)
        Test.@test q.metric_source == "quast"
        Test.@test q.quast_genome_fraction == 94.37
        Test.@test q.quast_nga50 == 8421.0
        Test.@test q.quast_num_misassemblies == 3.0
        Test.@test q.quast_duplication_ratio == 1.02

        # --- Fallback: no report.tsv -> internal source, all quast_* missing ----
        missing_report = joinpath(dir, "does_not_exist", "report.tsv")
        f = quast_metrics_for_report(missing_report)
        Test.@test f.metric_source == "internal"
        Test.@test f.quast_genome_fraction === missing
        Test.@test f.quast_nga50 === missing
        Test.@test f.quast_num_misassemblies === missing
        Test.@test f.quast_duplication_ratio === missing
    end
end

Test.@testset "DataFrame schema push! guards column-name drift" begin
    # This test constructs the EXACT production DataFrame schema from the sweep
    # (run_sweep's `rows`) and push!es the QUAST helper NamedTuples into it. It
    # catches drift between the helper's field names and the DataFrame columns —
    # a mismatch that would only surface at full-run time otherwise. Kept
    # Mycelia-free: only DataFrames (a Mycelia-project dep) is needed.
    df = DataFrames.DataFrame(
        reference = String[], genome_len = Int[], error_rate = Float64[],
        regime = String[], readlen = Int[], tech = String[], target_coverage = Float64[],
        effective_coverage = Float64[], k = Int[], arm = String[], ok = Bool[],
        n_contigs = Int[], total_length = Int[], largest_contig = Int[], n50 = Int[],
        genome_fraction = Float64[], runtime_s = Float64[],
        quast_genome_fraction = Union{Missing, Float64}[],
        quast_nga50 = Union{Missing, Float64}[],
        quast_num_misassemblies = Union{Missing, Float64}[],
        quast_duplication_ratio = Union{Missing, Float64}[],
        metric_source = String[]
    )

    # Non-QUAST columns for a row (mirrors what run_arm + the loop supply).
    base = (
        reference = "synthetic_2000bp", genome_len = 2_000, error_rate = 0.01,
        regime = "short-low-error", readlen = 150, tech = "illumina",
        target_coverage = 30.0, effective_coverage = 29.5, k = 21, arm = "naive",
        ok = true, n_contigs = 5, total_length = 1_980, largest_contig = 620,
        n50 = 410, genome_fraction = 99.0, runtime_s = 1.234
    )

    mktempdir() do dir
        # QUAST-present fixture (same shape as the parsing testset above).
        report_tsv = joinpath(dir, "report.tsv")
        open(report_tsv, "w") do io
            println(io, "Assembly\tnaive_contigs")
            println(io, "Genome fraction (%)\t94.37")
            println(io, "Duplication ratio\t1.02")
            println(io, "NGA50\t8421")
            println(io, "# misassemblies\t3")
        end

        # Row 1: a QUAST-scored arm. merge(base, quast) must slot cleanly into df.
        q = quast_metrics_for_report(report_tsv)
        Test.@test_nowarn push!(df, merge(base, q))

        # Row 2: an internal-fallback arm via the shared helper.
        e = empty_quast_metrics("internal:quast-disabled")
        Test.@test_nowarn push!(df, merge(base, e))
    end

    # Both pushes succeeded with no field-name mismatch, and provenance labels
    # landed in the metric_source column exactly as the helpers set them.
    Test.@test DataFrames.nrow(df) == 2
    Test.@test df.metric_source == ["quast", "internal:quast-disabled"]
    Test.@test df.quast_nga50[1] == 8421.0
    Test.@test df.quast_genome_fraction[2] === missing
end
