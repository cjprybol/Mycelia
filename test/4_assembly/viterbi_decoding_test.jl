# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_decoding_test.jl")'
# ```

import Test
import Mycelia
import MetaGraphsNext

Test.@testset "Viterbi Decoding" begin
    function create_ambiguous_decoding_graph()
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type = String,
            vertex_data_type = Any,
            edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
            weight_function = Mycelia.Rhizomorph.edge_data_weight,
            default_weight = 0.0
        )

        for label in ("S", "A", "B", "T", "X")
            graph[label] = nothing
        end

        graph["S", "A"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            9.0,
            Mycelia.Rhizomorph.Forward,
            Mycelia.Rhizomorph.Forward
        )
        graph["S", "B"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            4.0,
            Mycelia.Rhizomorph.Forward,
            Mycelia.Rhizomorph.Forward
        )
        graph["A", "T"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            1.0,
            Mycelia.Rhizomorph.Forward,
            Mycelia.Rhizomorph.Forward
        )
        graph["A", "X"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            9.0,
            Mycelia.Rhizomorph.Forward,
            Mycelia.Rhizomorph.Forward
        )
        graph["B", "T"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            1.0,
            Mycelia.Rhizomorph.Forward,
            Mycelia.Rhizomorph.Forward
        )

        return graph
    end

    path_labels(path) = [step.vertex_label for step in path.steps]

    Test.@testset "exact decoding beats local greedy choice" begin
        graph = create_ambiguous_decoding_graph()

        greedy = Mycelia.Rhizomorph.maximum_weight_walk_next(graph, "S", 2)
        Test.@test path_labels(greedy) == ["S", "A", "X"]

        exact = Mycelia.Rhizomorph.viterbi_decode_next(graph, "S", 2; target_vertex = "T")
        Test.@test exact isa Mycelia.Rhizomorph.ViterbiDecodingResult
        Test.@test exact.path !== nothing

        exact_path = something(exact.path)
        Test.@test path_labels(exact_path) == ["S", "B", "T"]
        Test.@test exact.score ≈ log(4.0 / 13.0)
        Test.@test exact_path.total_probability ≈ 4.0 / 13.0
        Test.@test exact.diagnostics[:algorithm] == :viterbi_exact
        Test.@test exact.diagnostics[:score_domain] == :log_probability
        Test.@test exact.diagnostics[:transition_scoring] == :normalized_edge_weight
        Test.@test exact.diagnostics[:reached_target] == true
        Test.@test exact.diagnostics[:path_length] == 3
        Test.@test exact.diagnostics[:expanded_states] > 0
    end

    Test.@testset "beam-pruned decoding exposes pruning diagnostics" begin
        graph = create_ambiguous_decoding_graph()

        exact = Mycelia.Rhizomorph.viterbi_decode_next(graph, "S", 2; target_vertex = "T")
        beam = Mycelia.Rhizomorph.beam_pruned_viterbi_decode_next(
            graph,
            "S",
            2;
            target_vertex = "T"
        )

        Test.@test beam.path !== nothing
        Test.@test path_labels(something(beam.path)) == path_labels(something(exact.path))
        Test.@test beam.score ≈ exact.score
        Test.@test beam.diagnostics[:algorithm] == :viterbi_beam_pruned
        Test.@test beam.diagnostics[:beam_width] == 64
        Test.@test beam.diagnostics[:beam_log_threshold] == Inf
        Test.@test haskey(beam.diagnostics, :retained_states)
        Test.@test haskey(beam.diagnostics, :pruned_states)

        narrow = Mycelia.Rhizomorph.beam_pruned_viterbi_decode_next(
            graph,
            "S",
            2;
            target_vertex = "T",
            beam_width = 1
        )
        Test.@test narrow.path === nothing
        Test.@test narrow.score == -Inf
        Test.@test narrow.diagnostics[:reached_target] == false
        Test.@test narrow.diagnostics[:reason] == :target_unreachable
        Test.@test narrow.diagnostics[:pruned_states] > 0
    end

    Test.@testset "input validation and unreachable targets" begin
        graph = create_ambiguous_decoding_graph()

        self_target = Mycelia.Rhizomorph.viterbi_decode_next(
            graph,
            "S",
            0;
            target_vertex = "S"
        )
        Test.@test self_target.path !== nothing
        Test.@test path_labels(something(self_target.path)) == ["S"]
        Test.@test self_target.score == 0.0
        Test.@test self_target.diagnostics[:reached_target] == true

        missing = Mycelia.Rhizomorph.viterbi_decode_next(
            graph,
            "S",
            2;
            target_vertex = "missing"
        )
        Test.@test missing.path === nothing
        Test.@test missing.score == -Inf
        Test.@test missing.diagnostics[:reason] == :target_not_found

        Test.@test_throws ArgumentError Mycelia.Rhizomorph.viterbi_decode_next(
            graph,
            "S",
            -1
        )
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.viterbi_decode_next(
            graph,
            "missing",
            2
        )
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.beam_pruned_viterbi_decode_next(
            graph,
            "S",
            2;
            beam_width = 0
        )
    end
end
