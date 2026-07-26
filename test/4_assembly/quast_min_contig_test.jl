# Unit test for the QUAST --min-contig threshold policy (bead td-28o0).
#
# The load-bearing property is NOT "the new formula equals the old one" — it
# deliberately does not, for T4. It is:
#
#   (a) every reference the old expression already handled keeps its EXACT
#       threshold, so results already committed and the Lambda runs in flight stay
#       comparable, and
#   (b) T4 (and anything larger) now gets a threshold a fragmented assembly can
#       actually satisfy, so QUAST runs instead of silently degrading to internal
#       metrics.
#
# Both are asserted against the specific documented values, not against a
# re-implementation of the old expression.
#
# Dependency-free: includes only the pure policy helper. No Mycelia, no QUAST.

import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "quast_min_contig.jl"))

# The expression this policy replaces, kept ONLY so the test can state exactly
# which references move and which do not.
_qmc_legacy(glen) = max(50, glen ÷ 10)

function _qmc_throws_with(f, needle)
    try
        f()
    catch e
        msg = sprint(showerror, e)
        Test.@test occursin(needle, msg)
        return nothing
    end
    Test.@test false  # expected a throw, got none
    return nothing
end

Test.@testset "QUAST --min-contig policy (td-28o0)" begin
    lambda_len = 48_502
    t4_len = 168_903
    phi29_len = 19_282
    sars_len = 29_903
    ecoli_len = 4_641_652

    Test.@testset "preserved thresholds — comparability with committed results" begin
        # These are the values the harness has ALREADY used. A change here would
        # silently redefine the metric for genomes that were working, and would
        # break comparability with Lawrencium jobs 24247923 / 24247924 (Lambda
        # 50x / 100x) that are in flight.
        Test.@test quast_min_contig(lambda_len) == 4_850
        Test.@test quast_min_contig(phi29_len) == 1_928
        Test.@test quast_min_contig(sars_len) == 2_990

        # Same assertion stated the other way: unchanged relative to the legacy
        # expression for every reference the legacy expression could serve.
        for glen in (246, 300, 359, 1_000, 5_000, phi29_len, sars_len, lambda_len)
            Test.@test quast_min_contig(glen) == _qmc_legacy(glen)
        end
    end

    Test.@testset "T4 is fixed — threshold is now reachable" begin
        # The bug: 168_903 ÷ 10 = 16_890, but the naive arm's largest contig is
        # 11_622 bp, so QUAST exits "None of the assembly files contains correct
        # contigs... decrease --min-contig threshold".
        Test.@test _qmc_legacy(t4_len) == 16_890
        t4_naive_largest_contig = 11_622
        Test.@test _qmc_legacy(t4_len) > t4_naive_largest_contig      # the failure
        Test.@test quast_min_contig(t4_len) == 5_000
        Test.@test quast_min_contig(t4_len) < t4_naive_largest_contig # the fix

        # And it stays reachable as references grow, rather than getting worse.
        Test.@test quast_min_contig(ecoli_len) == MIN_CONTIG_CEILING_BP
        Test.@test quast_min_contig(typemax(Int) ÷ 2) == MIN_CONTIG_CEILING_BP
    end

    Test.@testset "viroid-scale down-scaling (the original intent) still works" begin
        # The whole reason the ÷10 term exists: bring the threshold BELOW QUAST's
        # 500 bp default for references shorter than QUAST's default assumes.
        for glen in (246, 300, 359)
            Test.@test quast_min_contig(glen) == MIN_CONTIG_FLOOR_BP
            Test.@test quast_min_contig(glen) < 500
        end
        Test.@test quast_min_contig(0) == MIN_CONTIG_FLOOR_BP
        # 5_000 bp is where the ÷10 term crosses QUAST's own default.
        Test.@test quast_min_contig(5_000) == 500
    end

    Test.@testset "monotone, and bounded on both sides" begin
        prev = 0
        for glen in 0:997:2_000_000
            v = quast_min_contig(glen)
            Test.@test MIN_CONTIG_FLOOR_BP <= v <= MIN_CONTIG_CEILING_BP
            Test.@test v >= prev            # non-decreasing in genome length
            prev = v
        end
    end

    Test.@testset "arm independence — the property that keeps arms comparable" begin
        # The rejected "derive from the observed contig distribution" option would
        # make the threshold a function of the assembly. This policy takes only
        # the reference, so both arms of a cell are filtered identically. Asserted
        # structurally: the function has no way to see an assembly.
        Test.@test quast_min_contig(t4_len) == quast_min_contig(t4_len)
        Test.@test length(methods(quast_min_contig)) == 1
    end

    Test.@testset "keyword overrides" begin
        Test.@test quast_min_contig(lambda_len; ceiling_bp = 500) == 500
        Test.@test quast_min_contig(lambda_len; divisor = 100) == 485
        Test.@test quast_min_contig(100; floor_bp = 1) == 10
        Test.@test quast_min_contig(1; floor_bp = 1) == 1
    end

    Test.@testset "argument validation" begin
        _qmc_throws_with(() -> quast_min_contig(-1), "genome_len must be non-negative")
        _qmc_throws_with(() -> quast_min_contig(1_000; divisor = 0), "divisor must be >= 1")
        _qmc_throws_with(() -> quast_min_contig(1_000; floor_bp = 0), "floor_bp must be >= 1")
        _qmc_throws_with(
            () -> quast_min_contig(1_000; floor_bp = 900, ceiling_bp = 100),
            "must be >= floor_bp")
    end

    Test.@testset "constants are the documented ones" begin
        Test.@test MIN_CONTIG_GENOME_DIVISOR == 10
        Test.@test MIN_CONTIG_FLOOR_BP == 50
        # Anchored to the longest read length on the sweep's read-regime axis
        # (MYCELIA_RGV_READLEN default "150,5000").
        Test.@test MIN_CONTIG_CEILING_BP == 5_000
        # The ceiling must sit at or above Lambda's threshold, or Lambda moves.
        Test.@test MIN_CONTIG_CEILING_BP >= _qmc_legacy(lambda_len)
    end

    Test.@testset "the sweep harness uses the policy, not an inline expression" begin
        sweep = read(
            joinpath(@__DIR__, "..", "..", "benchmarking",
                "rhizomorph_correction_validation_sweep.jl"), String)
        Test.@test occursin("include(joinpath(@__DIR__, \"quast_min_contig.jl\"))", sweep)
        Test.@test occursin("min_contig = quast_min_contig(glen)", sweep)
        # The inline expression that caused the T4 failure must be gone.
        Test.@test !occursin("min_contig = max(50, glen ÷ 10)", sweep)
    end
end
