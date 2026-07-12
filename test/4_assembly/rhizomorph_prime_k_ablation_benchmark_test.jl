# Test for the multi-period prime(coprime)-vs-composite(factor-aligning) k
# ablation harness (bead td-tjym).
#
# This test encodes the HONEST result of an experiment DESIGNED to detect
# factor-alignment. It asserts:
#   (1) the k-SIZE threshold structure — below it recovery collapses for every
#       fixture (incl. the non-repetitive control), above it recovery saturates;
#   (2) the ISOLATING COMPARISON — above the k-size threshold, every factor-
#       sharing composite k (p|k) recovers within a small margin of its size-
#       matched coprime prime k (deltas ~0), i.e. factor-alignment is NOT isolated
#       from k-size; and
#   (3) the HARNESS-SENSITIVITY collision control fires — proving the recovery
#       metric CAN detect an aliasing-driven deficit when one truly exists, so the
#       factor-alignment null in (2) is a real negative, not metric blindness.
#
# NB the coordinator's expected "dark diagonal where p|k" pattern does NOT appear;
# per the honesty directive the test asserts the observed null (small deltas) plus
# the sensitivity control, and the benchmark reports the full (period x k) table.

import DataFrames
import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_prime_k_ablation_benchmark.jl"))

Test.@testset "rhizomorph prime(coprime)-vs-composite(factor-aligning) k ablation" begin
    output_dir = mktempdir()
    artifacts = run_prime_k_ablation_benchmark(output_dir; write_plots = true)

    fixtures = ablation_fixtures()
    n_fixtures = length(fixtures)

    Test.@testset "structure" begin
        Test.@test artifacts.per_run_rows == n_fixtures * length(ABLATION_K)
        Test.@test isfile(artifacts.per_run_csv)
        Test.@test isfile(artifacts.factor_alignment_csv)
        Test.@test isfile(artifacts.collision_control_csv)
        Test.@test isfile(artifacts.figure_png)
        Test.@test isfile(artifacts.figure_svg)
    end

    per_run = artifacts.per_run_table
    Test.@testset "per-run values" begin
        Test.@test all(0.0 .<= per_run.mean_reference_recovery .<= 1.0)
        Test.@test all(0.0 .<= per_run.min_reference_recovery .<= 1.0)
        Test.@test Set(unique(per_run.k)) == Set(ABLATION_K)
        Test.@test Set(unique(per_run.k_class)) == Set(["prime", "composite"])
        for group in DataFrames.groupby(per_run, :dataset_id)
            Test.@test DataFrames.nrow(group) == length(ABLATION_K)
        end
    end

    Test.@testset "k-size threshold (harness is not blind; k-size dominates low k)" begin
        threshold = artifacts.size_threshold_k
        Test.@test threshold in ABLATION_K
        control = per_run[per_run.dataset_id .== CONTROL_NONREPETITIVE_ID, :]
        # Below the threshold the non-repetitive control collapses (pure k-size /
        # error-tolerance effect); at/above it, it and every tandem fixture recover.
        below = control[control.k .< threshold, :]
        Test.@test all(below.mean_reference_recovery .< ABLATION_SIZE_THRESHOLD_RECOVERY)
        Test.@test only(control[control.k .== threshold, :mean_reference_recovery]) >=
                   ABLATION_SIZE_THRESHOLD_RECOVERY
        top = per_run[per_run.k .== maximum(ABLATION_K), :]
        Test.@test all(top.mean_reference_recovery .>= 0.95)
    end

    Test.@testset "factor-alignment is NOT isolated from k-size (observed null)" begin
        fa = artifacts.factor_alignment_table
        # The experiment IS designed to detect factor-alignment: it must contain
        # the odd-period (p|k composite vs size-matched coprime prime) comparisons.
        Test.@test DataFrames.nrow(fa) > 0
        Test.@test all(isodd, fa.period)
        Test.@test all(fa.composite_k .% fa.period .== 0)
        Test.@test all(fa.composite_k .>= artifacts.size_threshold_k)
        Test.@test all(fa.prime_k .>= artifacts.size_threshold_k)
        # OBSERVED: every factor-sharing composite recovers within a small margin
        # of its size-matched coprime prime -> no isolated factor-alignment effect.
        Test.@test all(abs.(fa.recovery_delta_prime_minus_composite) .<
                       FACTOR_ALIGNMENT_MIN_DELTA)
        Test.@test artifacts.factor_alignment_isolated == false
    end

    Test.@testset "harness-sensitivity collision control fires" begin
        cc = artifacts.collision_control_table
        penalty_9 = only(cc[cc.k .== 9, :collision_penalty])
        penalty_11 = only(cc[cc.k .== 11, :collision_penalty])
        Test.@test penalty_9 >= COLLISION_MIN_GAP
        Test.@test (penalty_9 - penalty_11) >= COLLISION_MIN_GAP
        Test.@test artifacts.collision_control_fires
    end

    # Report the isolating comparison (informational).
    println("\nISOLATING COMPARISON (factor-sharing composite p|k vs size-matched coprime prime):")
    for row in eachrow(artifacts.factor_alignment_table)
        println("  p=$(row.period) composite k=$(row.composite_k) rec=$(round(row.composite_recovery; digits = 3)) " *
                "vs prime k=$(row.prime_k) rec=$(round(row.prime_recovery; digits = 3)) " *
                "delta=$(round(row.recovery_delta_prime_minus_composite; digits = 3))")
    end
    println("factor-alignment isolated: $(artifacts.factor_alignment_isolated) | " *
            "k-size threshold: k=$(artifacts.size_threshold_k) | " *
            "collision control fires: $(artifacts.collision_control_fires)")
end
