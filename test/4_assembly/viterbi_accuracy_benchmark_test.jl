# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_accuracy_benchmark_test.jl")'
# ```

import CSV
import DataFrames
import MetaGraphsNext
import Mycelia
import Statistics
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

    # Control A — over-correction on un-corrupted input.
    overcorrection = CSV.read(artifacts.overcorrection_csv, DataFrames.DataFrame)
    Test.@test artifacts.overcorrection_rows == 9
    Test.@test DataFrames.nrow(overcorrection) == 9
    Test.@test all(overcorrection.injected_error_count .== 0)
    Test.@test all(overcorrection.over_correction_edit_distance .>= 0)
    Test.@test all(0.0 .<= overcorrection.over_correction_rate .<= 1.0)
    # The specificity claim the manuscript needs: the corrector does not damage
    # already-correct input (abstains). Asserted as an aggregate so a single
    # tiny-fixture edit does not flake the suite.
    Test.@test Statistics.mean(overcorrection.over_correction_rate) < 0.01

    # Control B — shuffled-weight + random-rewire nulls.
    null_control = CSV.read(artifacts.null_control_csv, DataFrames.DataFrame)
    Test.@test artifacts.null_control_rows == 9
    Test.@test DataFrames.nrow(null_control) == 9
    Test.@test all(null_control.null_seed .== 20260711)
    for col in (:real_injected_error_recall, :weight_null_injected_error_recall,
        :rewire_null_injected_error_recall)
        Test.@test all(0.0 .<= null_control[!, col] .<= 1.0)
    end
    # Robust invariant: neither null may recover injected errors BETTER than the
    # real graph (a degraded graph cannot help more than the true one). The
    # magnitude of the rewire-null collapse is a reported finding, not asserted
    # here, so a surprising non-collapse surfaces as data rather than a red test.
    Test.@test Statistics.mean(null_control.real_injected_error_recall) >=
               Statistics.mean(null_control.weight_null_injected_error_recall) - 1e-9
    Test.@test Statistics.mean(null_control.real_injected_error_recall) >=
               Statistics.mean(null_control.rewire_null_injected_error_recall) - 1e-9
end

Test.@testset "shuffled-weight null preserves topology + weight multiset" begin
    # _shuffle_weighted_graph_weights must permute weights only: same edge set,
    # same multiset of weights, so any recovery collapse is attributable purely
    # to the weight->edge reassignment, not to a changed graph.
    fixture = first(viterbi_accuracy_fixtures())
    config = Mycelia.ViterbiCorrectionConfig(error_rate = 0.05)
    weighted = Mycelia.build_correction_weighted_graph(fixture.graph; config = config)
    before_edges = Set(collect(MetaGraphsNext.edge_labels(weighted)))
    before_weights = sort([weighted[s, d].weight for (s, d) in before_edges])

    shuffled = _shuffle_weighted_graph_weights(weighted, 20260711)
    after_edges = Set(collect(MetaGraphsNext.edge_labels(shuffled)))
    after_weights = sort([shuffled[s, d].weight for (s, d) in after_edges])

    Test.@test before_edges == after_edges
    Test.@test before_weights ≈ after_weights
end

Test.@testset "random-rewire null preserves vertices/edge-count/weight multiset" begin
    # _random_rewired_weighted_graph must keep the same vertex set, the same edge
    # COUNT, and the same weight multiset — only the topology (which vertices an
    # edge joins) is randomized, so a recovery collapse is attributable purely to
    # the destroyed adjacency, not to a smaller/heavier graph.
    fixture = first(viterbi_accuracy_fixtures())
    config = Mycelia.ViterbiCorrectionConfig(error_rate = 0.05)
    weighted = Mycelia.build_correction_weighted_graph(fixture.graph; config = config)
    before_vertices = Set(collect(MetaGraphsNext.labels(weighted)))
    before_edges = collect(MetaGraphsNext.edge_labels(weighted))
    before_weights = sort([weighted[s, d].weight for (s, d) in before_edges])

    rewired = _random_rewired_weighted_graph(weighted, 20260711)
    after_vertices = Set(collect(MetaGraphsNext.labels(rewired)))
    after_edges = collect(MetaGraphsNext.edge_labels(rewired))
    after_weights = sort([rewired[s, d].weight for (s, d) in after_edges])

    Test.@test before_vertices == after_vertices
    Test.@test length(after_edges) == length(before_edges)
    Test.@test before_weights ≈ after_weights
    # Topology genuinely changed (not a no-op permutation): at least some edges
    # differ from the real graph. (With hundreds of vertices this is ~certain.)
    Test.@test Set(after_edges) != Set(before_edges)
end
