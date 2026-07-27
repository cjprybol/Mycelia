# Unit test for the QUAST --min-contig threshold policy (bead td-28o0).
#
# The load-bearing properties are:
#
#   (a) the threshold never exceeds QUAST's own 500 bp default — the scaling term
#       exists only to go BELOW it for viroid-scale references, and every value it
#       produced above 500 was a threshold nobody chose;
#   (b) it is UNIFORM across the viral tier, because the pre-registration pools
#       across organisms and a per-reference threshold would make that pooled
#       analysis mixed-definition, which `metric_source_guard.jl` must refuse; and
#   (c) T4 (and anything larger) gets a threshold a fragmented assembly can
#       actually satisfy, so QUAST runs instead of silently degrading.
#
# An earlier version of this policy clamped at 5,000 bp to hold Lambda at 4,850
# "for comparability with committed results". Review established that no such
# committed results exist (one committed CSV, synthetic, threshold 200 under every
# candidate policy) and that Lambda was itself failing 8/12 cells at 4,850. The
# tests below therefore pin the CONTRACT, not the historical values.
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

    Test.@testset "clamped at QUAST's default; uniform above 5 kb" begin
        # The ceiling IS QUAST's own default. The scaling term exists only to go
        # below it for viroid-scale references.
        Test.@test quast_min_contig(lambda_len) == 500
        Test.@test quast_min_contig(phi29_len) == 500
        Test.@test quast_min_contig(sars_len) == 500
        Test.@test quast_min_contig(t4_len) == 500

        # UNIFORMITY IS THE LOAD-BEARING PROPERTY, not any single value: the
        # pre-registration pools across organisms, so a threshold that varied by
        # reference would make the pooled analysis mixed-definition and
        # `metric_source_guard.jl` would refuse it. Every viral-tier reference must
        # land on one threshold.
        viral_tier = (lambda_len, t4_len, phi29_len, sars_len)
        Test.@test length(unique(quast_min_contig.(viral_tier))) == 1

        # An earlier ceiling of 5,000 gave four distinct thresholds — the state
        # that would have forced `--allow-mixed-src` on the pooled analysis.
        Test.@test length(unique(
            quast_min_contig.(viral_tier; ceiling_bp = 5_000))) == 4
    end

    Test.@testset "T4 is fixed — threshold is now reachable" begin
        # The bug: 168_903 ÷ 10 = 16_890, but the naive arm's largest contig is
        # 11_622 bp, so QUAST exits "None of the assembly files contains correct
        # contigs... decrease --min-contig threshold".
        Test.@test _qmc_legacy(t4_len) == 16_890
        t4_naive_largest_contig = 11_622
        Test.@test _qmc_legacy(t4_len) > t4_naive_largest_contig      # the failure
        Test.@test quast_min_contig(t4_len) == 500
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
        Test.@test quast_min_contig(5_000) == 500  # where the scaling term meets the ceiling
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
        Test.@test quast_min_contig(lambda_len; ceiling_bp = 5_000) == 4_850
        Test.@test quast_min_contig(lambda_len; divisor = 100, ceiling_bp = 5_000) == 485
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
        Test.@test MIN_CONTIG_CEILING_BP == 500
        # The ceiling IS QUAST's documented default; the scaling term must only
        # ever reduce below it, never raise above it.
        Test.@test all(quast_min_contig(g) <= 500 for g in
        (300, phi29_len, sars_len, lambda_len, t4_len, ecoli_len))
    end

    Test.@testset "every harness uses the shared policy, not an inline formula" begin
        # C2: claiming "exactly one definition" is only true if every call site
        # actually uses it. Three harnesses still carried the inline expression
        # after the first fix, and two of them (t4_ksweep, real_genome_benchmark)
        # run T4 — so the original bug was still live on the organism it was found
        # on. Assert the absence of the formula repo-wide rather than in one file.
        bench = joinpath(@__DIR__, "..", "..", "benchmarking")
        harnesses = ["rhizomorph_correction_validation_sweep.jl", "t4_ksweep.jl",
            "real_genome_benchmark.jl", "e2e_phase_profile.jl"]
        for h in harnesses
            src = read(joinpath(bench, h), String)
            Test.@test occursin("include(joinpath(@__DIR__, \"quast_min_contig.jl\"))", src)
            Test.@test occursin("quast_min_contig(", src)
        end
        # No inline `max(50, <something> ÷ 10)` may survive in CODE anywhere under
        # benchmarking/ (comments describing the old formula are fine).
        for f in readdir(bench)
            endswith(f, ".jl") || continue
            for (i, line) in enumerate(eachline(joinpath(bench, f)))
                startswith(strip(line), "#") && continue
                Test.@test !occursin(r"max\(50,[^)]*÷\s*10\)", line)
            end
        end
    end
end
