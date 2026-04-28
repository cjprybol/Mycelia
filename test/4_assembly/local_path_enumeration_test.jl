import Test
import Mycelia
import MetaGraphsNext

function _local_path_enumeration_graph()
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph();
        label_type = String,
        vertex_data_type = Any,
        edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
        weight_function = Mycelia.Rhizomorph.edge_data_weight,
        default_weight = 0.0
    )

    for vertex in ["A", "B", "C", "E", "D"]
        graph[vertex] = nothing
    end

    graph["A", "B"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
        9.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)
    graph["A", "C"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
        3.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)
    graph["A", "E"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
        1.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)
    graph["B", "D"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
        1.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)
    graph["C", "D"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
        1.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)
    graph["E", "D"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
        1.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)

    return graph
end

function _linear_depth_guard_graph()
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph();
        label_type = String,
        vertex_data_type = Any,
        edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
        weight_function = Mycelia.Rhizomorph.edge_data_weight,
        default_weight = 0.0
    )

    for vertex in ["A", "B", "C", "D"]
        graph[vertex] = nothing
    end

    for (src, dst) in [("A", "B"), ("B", "C"), ("C", "D")]
        graph[src, dst] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            1.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)
    end

    return graph
end

function _string_superbubble_graph()
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph();
        label_type = String,
        vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
        edge_data_type = Mycelia.Rhizomorph.StringEdgeData,
        weight_function = Mycelia.Rhizomorph.compute_edge_weight
    )

    for vertex in ["A", "B", "C", "D"]
        graph[vertex] = Mycelia.Rhizomorph.StringVertexData(vertex)
    end

    for (src, dst) in [("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")]
        graph[src, dst] = Mycelia.Rhizomorph.StringEdgeData(1)
    end

    return graph
end

function _alternative_labels(alternative)
    return [step.vertex_label for step in alternative.path.steps]
end

Test.@testset "Local Path Enumeration" begin
    Test.@testset "ranked alternatives include score and provenance" begin
        graph = _local_path_enumeration_graph()

        result = Mycelia.Rhizomorph.enumerate_local_paths(
            graph, "A", "D"; max_paths = 5, max_depth = 3)

        Test.@test result isa Mycelia.Rhizomorph.LocalPathEnumerationResult{String}
        Test.@test length(result.alternatives) == 3
        Test.@test [alternative.rank for alternative in result.alternatives] == [1, 2, 3]
        Test.@test _alternative_labels.(result.alternatives) == [
            ["A", "B", "D"],
            ["A", "C", "D"],
            ["A", "E", "D"]
        ]
        Test.@test result.alternatives[1].score ≈ 9.0 / 13.0
        Test.@test result.alternatives[2].score ≈ 3.0 / 13.0
        Test.@test result.alternatives[3].score ≈ 1.0 / 13.0
        Test.@test all(
            alternative.score == alternative.path.total_probability
            for alternative in result.alternatives
        )
        Test.@test all(
            alternative.provenance == result.provenance for alternative in result.alternatives)
        Test.@test result.provenance.method == :bounded_best_first
        Test.@test result.provenance.entry_vertex == "A"
        Test.@test result.provenance.exit_vertex == "D"
        Test.@test !result.provenance.truncated
        Test.@test result.provenance.reason == :complete
    end

    Test.@testset "max_paths guard reports bounded behavior" begin
        graph = _local_path_enumeration_graph()

        result = Mycelia.Rhizomorph.enumerate_local_paths(
            graph, "A", "D"; max_paths = 2, max_depth = 3)

        Test.@test length(result.alternatives) == 2
        Test.@test _alternative_labels.(result.alternatives) == [
            ["A", "B", "D"],
            ["A", "C", "D"]
        ]
        Test.@test result.provenance.truncated
        Test.@test result.provenance.reason == :max_paths
    end

    Test.@testset "max_depth guard reports incomplete repeat-like regions" begin
        graph = _linear_depth_guard_graph()

        result = Mycelia.Rhizomorph.enumerate_local_paths(
            graph, "A", "D"; max_paths = 5, max_depth = 2)

        Test.@test isempty(result.alternatives)
        Test.@test result.provenance.truncated
        Test.@test result.provenance.reason == :max_depth
    end

    Test.@testset "max_expansions guard reports bounded behavior" begin
        graph = _local_path_enumeration_graph()

        result = Mycelia.Rhizomorph.enumerate_local_paths(
            graph, "A", "D"; max_paths = 5, max_depth = 3, max_expansions = 1)

        Test.@test isempty(result.alternatives)
        Test.@test result.provenance.expansions == 1
        Test.@test result.provenance.truncated
        Test.@test result.provenance.reason == :max_expansions
    end

    Test.@testset "raw Rhizomorph graphs are converted before enumeration" begin
        graph = _string_superbubble_graph()

        result = Mycelia.Rhizomorph.enumerate_local_paths(
            graph, "A", "D"; max_paths = 3, max_depth = 2)

        Test.@test length(result.alternatives) == 2
        Test.@test _alternative_labels.(result.alternatives) == [
            ["A", "B", "D"],
            ["A", "C", "D"]
        ]
        Test.@test result.alternatives[1].score ≈ 0.5
        Test.@test result.alternatives[2].score ≈ 0.5
        Test.@test !result.provenance.truncated
    end

    Test.@testset "superbubble wrappers enumerate multi-branch alternatives" begin
        graph = _local_path_enumeration_graph()
        bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
            graph; min_bubble_length = 1, max_bubble_length = 5)

        Test.@test length(bubbles) == 1

        result = Mycelia.Rhizomorph.enumerate_superbubble_paths(
            graph, first(bubbles); max_paths = 5, max_depth = 3)

        Test.@test length(result.alternatives) == 3
        Test.@test result.provenance.method == :superbubble
        Test.@test _alternative_labels.(result.alternatives) == [
            ["A", "B", "D"],
            ["A", "C", "D"],
            ["A", "E", "D"]
        ]

        detected_results = Mycelia.Rhizomorph.enumerate_superbubble_paths(
            graph; min_bubble_length = 1, max_bubble_length = 5, max_bubbles = 1)

        Test.@test length(detected_results) == 1
        Test.@test length(first(detected_results).alternatives) == 3
    end
end
