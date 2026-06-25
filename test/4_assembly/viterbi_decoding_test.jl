# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_decoding_test.jl")'
# ```

import Test
import Mycelia
import MetaGraphsNext

Test.@testset "Viterbi Decoding" begin
    function create_weighted_graph(labels)
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type = eltype(collect(labels)),
            vertex_data_type = Any,
            edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
            weight_function = Mycelia.Rhizomorph.edge_data_weight,
            default_weight = 0.0
        )

        for label in labels
            graph[label] = nothing
        end

        return graph
    end

    function add_weighted_edge!(
            graph,
            source,
            target,
            weight;
            source_strand = Mycelia.Rhizomorph.Forward,
            target_strand = Mycelia.Rhizomorph.Forward
    )
        graph[source, target] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            weight,
            source_strand,
            target_strand
        )
        return graph
    end

    function create_ambiguous_decoding_graph()
        graph = create_weighted_graph(["S", "A", "B", "T", "X"])

        add_weighted_edge!(graph, "S", "A", 9.0)
        add_weighted_edge!(graph, "S", "B", 4.0)
        add_weighted_edge!(graph, "A", "T", 1.0)
        add_weighted_edge!(graph, "A", "X", 9.0)
        add_weighted_edge!(graph, "B", "T", 1.0)

        return graph
    end

    path_labels(path) = [step.vertex_label for step in path.steps]
    path_strands(path) = [step.strand for step in path.steps]
    path_step_probabilities(path) = [step.probability for step in path.steps]


    Test.@testset "correction interface delegates to graph-agnostic decoder" begin
        graph = create_weighted_graph(["S", "T"])
        add_weighted_edge!(graph, "S", "T", 1.0)

        result = Mycelia.correct_observations(graph, [["S", "T"]])

        Test.@test result isa Mycelia.ViterbiCorrectionResult
        Test.@test result.diagnostics[:interface] == :metagraphs_next
        Test.@test result.diagnostics[:algorithm] == :rhizomorph_viterbi_decode_next
        Test.@test length(result.paths) == 1
        Test.@test result.paths[1] isa Mycelia.Rhizomorph.ViterbiDecodingResult
        Test.@test result.paths[1].path !== nothing
        Test.@test path_labels(something(result.paths[1].path)) == ["S", "T"]
    end

    Test.@testset "exact decoding beats local greedy choice" begin
        graph = create_ambiguous_decoding_graph()

        greedy = Mycelia.Rhizomorph.maximum_weight_walk_next(graph, "S", 2)
        Test.@test path_labels(greedy) == ["S", "A", "X"]

        exact = Mycelia.Rhizomorph.viterbi_decode_next(graph, "S", 2; target_vertex = "T")
        Test.@test exact isa Mycelia.Rhizomorph.ViterbiDecodingResult
        Test.@test exact.path !== nothing

        exact_path = something(exact.path)
        Test.@test path_labels(exact_path) == ["S", "B", "T"]
        Test.@test path_strands(exact_path) == fill(Mycelia.Rhizomorph.Forward, 3)
        Test.@test path_step_probabilities(exact_path) ≈ [1.0, 4.0 / 13.0, 1.0]
        Test.@test exact.score ≈ log(4.0 / 13.0)
        Test.@test exact_path.total_probability ≈ exp(exact.score)
        Test.@test exact_path.total_probability ≈ 4.0 / 13.0
        Test.@test exact.diagnostics[:algorithm] == :viterbi_exact
        Test.@test exact.diagnostics[:score_domain] == :log_probability
        Test.@test exact.diagnostics[:transition_scoring] == :normalized_edge_weight
        Test.@test exact.diagnostics[:reached_target] == true
        Test.@test exact.diagnostics[:path_length] == 3
        Test.@test exact.diagnostics[:expanded_states] == 3
        Test.@test exact.diagnostics[:generated_states] == 5
        Test.@test exact.diagnostics[:retained_states] == 2
        Test.@test exact.diagnostics[:cumulative_retained_states] == 5
        Test.@test exact.diagnostics[:max_retained_states] == 2
        Test.@test exact.diagnostics[:completed_steps] == 2
        Test.@test exact.diagnostics[:skipped_transitions] == 0
    end

    Test.@testset "no-target decoding returns best deepest frontier path" begin
        graph = create_ambiguous_decoding_graph()

        exact = Mycelia.Rhizomorph.viterbi_decode_next(graph, "S", 2)

        Test.@test exact.path !== nothing
        exact_path = something(exact.path)
        Test.@test path_labels(exact_path) == ["S", "A", "X"]
        Test.@test exact.score ≈ log((9.0 / 13.0) * (9.0 / 10.0))
        Test.@test exact_path.total_probability ≈ exp(exact.score)
        Test.@test exact.diagnostics[:target_vertex] === nothing
        Test.@test exact.diagnostics[:reached_target] === nothing
        Test.@test exact.diagnostics[:path_length] == 3
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
        Test.@test beam.diagnostics[:retained_states] == 2
        Test.@test beam.diagnostics[:pruned_states] == 0

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
        Test.@test narrow.diagnostics[:retained_states] == 1
        Test.@test narrow.diagnostics[:pruned_states] == 2

        threshold_only = Mycelia.Rhizomorph.beam_pruned_viterbi_decode_next(
            graph,
            "S",
            2;
            target_vertex = "T",
            beam_width = 64,
            beam_log_threshold = 0.0
        )
        Test.@test threshold_only.path === nothing
        Test.@test threshold_only.score == -Inf
        Test.@test threshold_only.diagnostics[:beam_width] == 64
        Test.@test threshold_only.diagnostics[:beam_log_threshold] == 0.0
        Test.@test threshold_only.diagnostics[:pruned_states] == 2
    end

    Test.@testset "finite beam prunes states yet recovers the exact optimum" begin
        # S branches to THREE successors; a width-2 beam keeps only the two
        # best-scoring partial paths (A, B) and prunes the third (C). Because the
        # pruned branch is not on the optimal path, beam must still return the
        # exact optimum — the meaningful beam==exact equivalence the trivial
        # no-pruning case (beam_width >= frontier) does not exercise.
        graph = create_weighted_graph(["S", "A", "B", "C", "T"])
        add_weighted_edge!(graph, "S", "A", 9.0)
        add_weighted_edge!(graph, "S", "B", 4.0)
        add_weighted_edge!(graph, "S", "C", 1.0)
        add_weighted_edge!(graph, "A", "T", 1.0)
        add_weighted_edge!(graph, "B", "T", 1.0)
        add_weighted_edge!(graph, "C", "T", 1.0)

        exact = Mycelia.Rhizomorph.viterbi_decode_next(graph, "S", 2; target_vertex = "T")
        beam = Mycelia.Rhizomorph.beam_pruned_viterbi_decode_next(
            graph, "S", 2; target_vertex = "T", beam_width = 2)

        Test.@test beam.path !== nothing
        Test.@test path_labels(something(beam.path)) == ["S", "A", "T"]
        Test.@test path_labels(something(beam.path)) == path_labels(something(exact.path))
        Test.@test beam.score ≈ exact.score
        # Pruning genuinely occurred (3-state frontier, width 2) yet the optimum survived.
        Test.@test beam.diagnostics[:pruned_states] >= 1
    end

    Test.@testset "input validation and target boundary cases" begin
        graph = create_ambiguous_decoding_graph()
        empty_graph = create_weighted_graph(String[])

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
        Test.@test self_target.diagnostics[:completed_steps] == 0

        self_target_with_horizon = Mycelia.Rhizomorph.viterbi_decode_next(
            graph,
            "S",
            2;
            target_vertex = "S"
        )
        Test.@test self_target_with_horizon.path !== nothing
        Test.@test path_labels(something(self_target_with_horizon.path)) == ["S"]
        Test.@test self_target_with_horizon.score == 0.0

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
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.viterbi_decode_next(
            empty_graph,
            "S",
            0
        )
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.beam_pruned_viterbi_decode_next(
            graph,
            "S",
            2;
            beam_width = 0
        )
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.beam_pruned_viterbi_decode_next(
            graph,
            "S",
            2;
            beam_log_threshold = -0.1
        )
    end

    Test.@testset "dead ends and unreachable targets are explicit" begin
        graph = create_weighted_graph(["S", "T"])

        unreachable = Mycelia.Rhizomorph.viterbi_decode_next(
            graph,
            "S",
            3;
            target_vertex = "T"
        )
        Test.@test unreachable.path === nothing
        Test.@test unreachable.score == -Inf
        Test.@test unreachable.diagnostics[:reason] == :target_unreachable
        Test.@test unreachable.diagnostics[:expanded_states] == 1
        Test.@test unreachable.diagnostics[:generated_states] == 0
        Test.@test unreachable.diagnostics[:completed_steps] == 0

        no_target = Mycelia.Rhizomorph.viterbi_decode_next(graph, "S", 3)
        Test.@test no_target.path !== nothing
        Test.@test path_labels(something(no_target.path)) == ["S"]
        Test.@test no_target.score == 0.0
        Test.@test no_target.diagnostics[:expanded_states] == 1
        Test.@test no_target.diagnostics[:completed_steps] == 0
    end

    Test.@testset "strand transitions are reconstructed" begin
        graph = create_weighted_graph(["S", "A", "T"])
        add_weighted_edge!(
            graph,
            "S",
            "A",
            2.0;
            source_strand = Mycelia.Rhizomorph.Forward,
            target_strand = Mycelia.Rhizomorph.Reverse
        )
        add_weighted_edge!(
            graph,
            "S",
            "T",
            1.0;
            source_strand = Mycelia.Rhizomorph.Forward,
            target_strand = Mycelia.Rhizomorph.Forward
        )
        add_weighted_edge!(
            graph,
            "A",
            "T",
            3.0;
            source_strand = Mycelia.Rhizomorph.Reverse,
            target_strand = Mycelia.Rhizomorph.Reverse
        )

        decoded = Mycelia.Rhizomorph.viterbi_decode_next(
            graph,
            "S",
            2;
            target_vertex = "T",
            start_strand = Mycelia.Rhizomorph.Forward
        )
        Test.@test decoded.path !== nothing
        decoded_path = something(decoded.path)
        Test.@test path_labels(decoded_path) == ["S", "A", "T"]
        Test.@test path_strands(decoded_path) == [
            Mycelia.Rhizomorph.Forward,
            Mycelia.Rhizomorph.Reverse,
            Mycelia.Rhizomorph.Reverse
        ]
        Test.@test decoded.score ≈ log(2.0 / 3.0)
        Test.@test decoded_path.total_probability ≈ 2.0 / 3.0

        reverse_start = Mycelia.Rhizomorph.viterbi_decode_next(
            graph,
            "S",
            2;
            target_vertex = "T",
            start_strand = Mycelia.Rhizomorph.Reverse
        )
        Test.@test reverse_start.path === nothing
        Test.@test reverse_start.diagnostics[:reason] == :target_unreachable
    end

    Test.@testset "weight flooring and tie-breaking are deterministic" begin
        graph = create_weighted_graph(["S", "B", "A"])
        add_weighted_edge!(graph, "S", "B", 0.0)
        add_weighted_edge!(graph, "S", "A", 1.0e-12)

        decoded = Mycelia.Rhizomorph.viterbi_decode_next(graph, "S", 1)

        Test.@test decoded.path !== nothing
        decoded_path = something(decoded.path)
        Test.@test path_labels(decoded_path) == ["S", "A"]
        Test.@test decoded.score ≈ log(0.5)
        Test.@test decoded_path.total_probability ≈ 0.5
        Test.@test decoded.diagnostics[:generated_states] == 2
        Test.@test decoded.diagnostics[:skipped_transitions] == 0
    end

    Test.@testset "path reconstruction validates predecessor layers" begin
        graph = create_ambiguous_decoding_graph()

        Test.@test_throws ArgumentError Mycelia.Rhizomorph._reconstruct_viterbi_path(
            graph,
            ("T", Mycelia.Rhizomorph.Forward),
            1,
            Dict{Tuple{String, Mycelia.Rhizomorph.StrandOrientation},
                Tuple{Tuple{String, Mycelia.Rhizomorph.StrandOrientation}, Float64}}[]
        )

        Test.@test_throws ArgumentError Mycelia.Rhizomorph._reconstruct_viterbi_path(
            graph,
            ("T", Mycelia.Rhizomorph.Forward),
            1,
            [
                Dict(
                ("X",
                Mycelia.Rhizomorph.Forward) => (("A", Mycelia.Rhizomorph.Forward), 0.5)
            ),
            ]
        )
    end
end
