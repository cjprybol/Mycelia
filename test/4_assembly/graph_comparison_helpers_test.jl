import Test
import Mycelia
import MetaGraphsNext
import Graphs

function make_string_path_graph(labels::Vector{String})
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
        edge_data_type = Mycelia.Rhizomorph.StringEdgeData
    )

    for label in labels
        graph[label] = Mycelia.Rhizomorph.StringVertexData(label)
    end

    if length(labels) >= 2
        for (src, dst) in zip(labels[1:(end - 1)], labels[2:end])
            overlap_length = min(length(src), length(dst)) - 1
            graph[src, dst] = Mycelia.Rhizomorph.StringEdgeData(overlap_length)
        end
    end

    return graph
end

Test.@testset "Graph Comparison Helpers" begin
    Test.@testset "Topological equality ignores labels" begin
        graph_a = make_string_path_graph(["ABC", "BCD", "CDE"])
        graph_b = make_string_path_graph(["JKL", "KLM", "LMN"])

        Test.@test is_topologically_equal(graph_a, graph_b)
        Test.@test !is_semantically_equal(graph_a, graph_b)
    end

    Test.@testset "Semantic equality tolerates different decompositions" begin
        graph_a = make_string_path_graph(["ABC", "BCD", "CDE"])
        graph_b = make_string_path_graph(["ABCD", "BCDE"])

        Test.@test !is_topologically_equal(graph_a, graph_b)
        Test.@test is_semantically_equal(graph_a, graph_b)
    end

    Test.@testset "Topology changes break both helpers when assembly changes" begin
        graph_a = make_string_path_graph(["ABC", "BCD", "CDE"])
        graph_b = make_string_path_graph(["ABC", "BCF", "CFG"])

        Test.@test is_topologically_equal(graph_a, graph_b)
        Test.@test !is_semantically_equal(graph_a, graph_b)
    end

    Test.@testset "Empty graphs compare equal" begin
        empty_a = make_string_path_graph(String[])
        empty_b = make_string_path_graph(String[])

        Test.@test is_topologically_equal(empty_a, empty_b)
        Test.@test is_semantically_equal(empty_a, empty_b)
    end
end
