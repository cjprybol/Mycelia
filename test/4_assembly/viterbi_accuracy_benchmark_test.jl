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
    # Cross-field accounting invariant: a row has zero changed observations iff it
    # has zero over-correction edit distance. Exercises the numerator logic even
    # though every observed row is the (expected) zero case — a mis-indexed or
    # inverted count would break this equivalence. (review: over-correction path)
    Test.@test all(
        (overcorrection.changed_observations .== 0) .==
        (overcorrection.over_correction_edit_distance .== 0))
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
    # real graph (a degraded graph cannot help more than the true one).
    Test.@test Statistics.mean(null_control.real_injected_error_recall) >=
               Statistics.mean(null_control.weight_null_injected_error_recall) - 1e-9
    Test.@test Statistics.mean(null_control.real_injected_error_recall) >=
               Statistics.mean(null_control.rewire_null_injected_error_recall) - 1e-9
    # POSITIVE CONTROL: without this, the "null <= real" invariant is vacuous —
    # it also passes if the corrector itself regressed to zero recovery (real ==
    # null == 0). The real graph must actually recover on these fixtures, and it
    # must never fall into the decode-collapse fallback.
    Test.@test Statistics.mean(null_control.real_injected_error_recall) > 0.5
    Test.@test all(null_control.real_decoded)
    # Fixture-scoped regression guard on the DISCRIMINATING signal: on these
    # pinned fixtures (RefSeq NC_002030.1 / NC_001422.1 + fixed text excerpt,
    # deterministic injector, fixed seed) destroying the true adjacency abolishes
    # recovery entirely. If pointed at new fixtures, re-baseline this one line;
    # the portable "null <= real" invariant above carries over unchanged.
    Test.@test Statistics.mean(null_control.rewire_null_injected_error_recall) == 0.0
    Test.@test all(.!null_control.rewire_null_decoded)
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
    # Aliasing guard: the shuffle must NOT mutate the caller's graph (it deepcopies
    # internally), so the null can never leak into the real/Control-A arms even if
    # a future fixture ships pre-weighted. Verify `weighted` is untouched and the
    # returned graph is a distinct object. (review: aliasing-landmine)
    Test.@test shuffled !== weighted
    Test.@test sort([weighted[s, d].weight
                     for (s, d) in
                         MetaGraphsNext.edge_labels(weighted)]) ≈ before_weights
    # build_correction_weighted_graph must itself return a fresh graph for these
    # fixtures (the precondition that makes the in-place null isolated).
    Test.@test Mycelia.build_correction_weighted_graph(fixture.graph; config = config) !==
               fixture.graph
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
    # No self-loops: the rewire rejects src == dst, so no vertex points to itself.
    # (The vertex/edge-count/weight tests above would all still pass with a
    # self-loop, so assert it explicitly.)
    Test.@test all(s != d for (s, d) in after_edges)
end
