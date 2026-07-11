# Structural test for the prime-vs-composite k ablation harness (bead td-tjym).
#
# Asserts the harness runs end-to-end at SMOKE scale and produces the expected
# shape of results: one row per (fixture x k x error rate x error mode), recall in
# [0, 1], one delta row per (fixture x error mode). It does NOT hard-assert a
# specific prime-advantage magnitude — that is the observed result, reported by
# the benchmark itself.

import DataFrames
import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_prime_k_ablation_benchmark.jl"))

Test.@testset "rhizomorph prime-vs-composite k ablation" begin
    output_dir = mktempdir()
    artifacts = run_prime_k_ablation_benchmark(output_dir; write_plots = true)

    fixtures = ablation_fixtures()
    n_fixtures = length(fixtures)
    n_modes = length(ABLATION_ERROR_MODES)
    expected_per_run = n_fixtures * length(ABLATION_K) * length(ABLATION_ERROR_RATES) *
                       n_modes

    Test.@testset "structure" begin
        Test.@test artifacts.per_run_rows == expected_per_run
        Test.@test artifacts.delta_rows == n_fixtures * n_modes
        Test.@test isfile(artifacts.per_run_csv)
        Test.@test isfile(artifacts.delta_csv)
        Test.@test isfile(artifacts.figure_png)
        Test.@test isfile(artifacts.figure_svg)
    end

    per_run = artifacts.per_run_table
    Test.@testset "per-run values" begin
        Test.@test all(0.0 .<= per_run.injected_error_recall .<= 1.0)
        Test.@test all(per_run.edit_distance_reduction .<= 1.0)
        Test.@test all(per_run.injected_error_count .>= 1)
        Test.@test all(per_run.corrected_mismatch .>= 0)
        Test.@test all(per_run.over_correction_positions .>= 0)
        Test.@test all(per_run.reconstructed_length .== per_run.sequence_length)
        Test.@test Set(unique(per_run.k)) == Set(ABLATION_K)
        Test.@test Set(unique(per_run.k_class)) == Set(["prime", "composite"])
        Test.@test Set(unique(per_run.error_mode)) ==
                   Set(string.(ABLATION_ERROR_MODES))
        # Every fixture is covered at every k, error rate, and error mode.
        for group in DataFrames.groupby(per_run, :dataset_id)
            Test.@test DataFrames.nrow(group) ==
                       length(ABLATION_K) * length(ABLATION_ERROR_RATES) * n_modes
        end
    end

    delta = artifacts.delta_table
    Test.@testset "delta table" begin
        Test.@test all(0.0 .<= delta.prime_mean_recall .<= 1.0)
        Test.@test all(0.0 .<= delta.composite_mean_recall .<= 1.0)
        Test.@test all(delta.recall_delta_prime_minus_composite .==
                       delta.prime_mean_recall .- delta.composite_mean_recall)
        Test.@test Set(delta.dataset_id) == Set(f.dataset_id for f in fixtures)
    end

    # Report the observed prime-vs-composite deltas (informational, not asserted).
    println("\nObserved prime-vs-composite deltas per fixture x error mode:")
    for drow in eachrow(delta)
        println("  $(rpad(drow.dataset_id, 24)) $(rpad(drow.error_mode, 9)) " *
                "prime_recall=$(round(drow.prime_mean_recall; digits = 4)) " *
                "composite_recall=$(round(drow.composite_mean_recall; digits = 4)) " *
                "recall_delta=$(round(drow.recall_delta_prime_minus_composite; digits = 4)) " *
                "edr_delta=$(round(drow.edit_distance_reduction_delta; digits = 4))")
    end
end
