# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_accuracy_benchmark_test.jl")'
# ```

import CSV
import DataFrames
import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "viterbi_accuracy_benchmark.jl"))

Test.@testset "B8 Viterbi accuracy benchmark artifacts" begin
    output_dir = mktempdir(prefix = "viterbi_accuracy_b8_test_")
    artifacts = run_viterbi_accuracy_benchmark(output_dir; write_plots = false)
    summary = CSV.read(artifacts.summary_csv, DataFrames.DataFrame)

    Test.@test artifacts.rows == 9
    Test.@test DataFrames.nrow(summary) == 9
    Test.@test Set(summary.dataset_id) == Set([
        "pstvd_nc002030", "phix174_nc001422", "austen_pride_prejudice_excerpt"
    ])
    Test.@test Set(summary.target_error_rate) == Set([0.01, 0.05, 0.10])
    Test.@test all(summary.baseline_edit_distance .> 0)
    Test.@test all(summary.corrected_edit_distance .<= summary.baseline_edit_distance)
    Test.@test all(summary.edit_distance_reduction .>= 0)
    Test.@test all(summary.injected_error_recall .>= 0)
    Test.@test isfile(artifacts.index)
    Test.@test isfile(artifacts.provenance)
end
