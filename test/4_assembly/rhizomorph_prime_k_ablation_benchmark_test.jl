# Test for the de-novo prime-vs-composite k ablation harness (bead td-tjym).
#
# This is a VALID-harness test: beyond structural checks it asserts the two
# properties that make a null meaningful — (1) the harness produces a strong
# POSITIVE recovery signal where recovery is expected (large-k, non-repetitive
# reads assemble to the ceiling), and (2) the POSITIVE CONTROL fires — an
# engineered composite-k (k=9) k-mer collision demonstrably degrades de-novo
# recovery relative to a size-matched prime k (k=11), and the deficit is isolated
# from k-size by a structurally-identical no-collision fixture. It does NOT
# hard-assert a primality effect on the natural fixtures — that is reported.

import DataFrames
import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_prime_k_ablation_benchmark.jl"))

Test.@testset "rhizomorph prime-vs-composite k ablation (de-novo)" begin
    output_dir = mktempdir()
    artifacts = run_prime_k_ablation_benchmark(output_dir; write_plots = true)

    fixtures = ablation_fixtures()
    n_fixtures = length(fixtures)
    expected_per_run = n_fixtures * length(ABLATION_K)

    Test.@testset "structure" begin
        Test.@test artifacts.per_run_rows == expected_per_run
        Test.@test artifacts.delta_rows == n_fixtures * length(SIZE_MATCHED_PAIRS)
        Test.@test isfile(artifacts.per_run_csv)
        Test.@test isfile(artifacts.delta_csv)
        Test.@test isfile(artifacts.positive_control_csv)
        Test.@test isfile(artifacts.figure_png)
        Test.@test isfile(artifacts.figure_svg)
    end

    per_run = artifacts.per_run_table
    Test.@testset "per-run values" begin
        Test.@test all(0.0 .<= per_run.mean_reference_recovery .<= 1.0)
        Test.@test all(0.0 .<= per_run.min_reference_recovery .<= 1.0)
        Test.@test all(0.0 .<= per_run.max_reference_recovery .<= 1.0)
        Test.@test all(per_run.mean_contig_count .>= 0)
        Test.@test Set(unique(per_run.k)) == Set(ABLATION_K)
        Test.@test Set(unique(per_run.k_class)) == Set(["prime", "composite"])
        # Every fixture is covered at every k.
        for group in DataFrames.groupby(per_run, :dataset_id)
            Test.@test DataFrames.nrow(group) == length(ABLATION_K)
        end
    end

    Test.@testset "positive recovery signal (harness is not blind)" begin
        # At the largest k, every fixture must assemble to near-complete recovery.
        # A harness that could not recover anything makes a null uninformative.
        max_k = maximum(ABLATION_K)
        top = per_run[per_run.k .== max_k, :]
        Test.@test all(top.all_seeds_assembled)
        Test.@test minimum(top.mean_reference_recovery) >= 0.95
        # The non-repetitive control recovers fully at large k.
        control = per_run[
        (per_run.dataset_id .== "control_nonrepetitive") .& (per_run.k .== max_k), :]
        Test.@test only(control.mean_reference_recovery) >= 0.95
    end

    Test.@testset "positive control fires (composite-k collision detected)" begin
        summary = artifacts.positive_control_summary
        penalty_9 = only(summary[summary.k .== 9, :collision_penalty])
        penalty_11 = only(summary[summary.k .== 11, :collision_penalty])
        # Composite k=9 collision imposes a real recovery deficit (structure held
        # identical, so this is not a k-size artifact) ...
        Test.@test penalty_9 >= POSITIVE_CONTROL_MIN_GAP
        # ... that is largely resolved at the size-matched prime k=11.
        Test.@test (penalty_9 - penalty_11) >= POSITIVE_CONTROL_MIN_GAP
        Test.@test artifacts.positive_control_fires
    end

    delta = artifacts.delta_table
    Test.@testset "delta table" begin
        Test.@test all(0.0 .<= delta.composite_recovery .<= 1.0)
        Test.@test all(0.0 .<= delta.prime_recovery .<= 1.0)
        Test.@test all(delta.recovery_delta_prime_minus_composite .==
                       delta.prime_recovery .- delta.composite_recovery)
        Test.@test Set(delta.dataset_id) == Set(f.dataset_id for f in fixtures)
    end

    # Report the observed size-matched prime-minus-composite recovery deltas
    # (informational — no primality effect is asserted on the natural fixtures).
    println("\nObserved size-matched prime-minus-composite recovery deltas:")
    for drow in eachrow(delta)
        println("  $(rpad(drow.dataset_id, 28)) c=$(drow.composite_k) p=$(drow.prime_k) " *
                "composite=$(round(drow.composite_recovery; digits = 3)) " *
                "prime=$(round(drow.prime_recovery; digits = 3)) " *
                "delta=$(round(drow.recovery_delta_prime_minus_composite; digits = 3))")
    end
    println("\nPositive-control collision penalty (isolates aliasing from k-size):")
    for srow in eachrow(artifacts.positive_control_summary)
        println("  k=$(srow.k) ($(srow.k_class)) shared=$(round(srow.shared_recovery; digits = 3)) " *
                "distinct=$(round(srow.distinct_recovery; digits = 3)) " *
                "penalty=$(round(srow.collision_penalty; digits = 3))")
    end
    println("Positive control fires: $(artifacts.positive_control_fires)")
end
