import Test
import Mycelia
import MetaGraphsNext
import Graphs

Test.@testset "String Graph Utilities" begin
    Test.@testset "Token graph construction" begin
        tokens = [
            ["the", "cat", "sat"],
            ["the", "dog"]
        ]
        graph = Mycelia.Rhizomorph.build_token_graph(tokens; dataset_id = "tokens")

        Test.@test MetaGraphsNext.haskey(graph, "the")
        Test.@test MetaGraphsNext.haskey(graph, "cat")
        Test.@test MetaGraphsNext.haskey(graph, "the", "cat")
        Test.@test MetaGraphsNext.haskey(graph, "cat", "sat")
        Test.@test Mycelia.Rhizomorph.count_evidence(graph["the"]) == 2
    end

    Test.@testset "Token graph errors" begin
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_token_graph(Vector{Vector{String}}())
    end

    Test.@testset "String graph from file(s)" begin
        mktempdir() do dir
            path1 = joinpath(dir, "strings1.txt")
            open(path1, "w") do io
                println(io, "ATCG")
                println(io, "TCGA")
                println(io, "")
            end

            graph = Mycelia.Rhizomorph.build_string_graph_from_file(path1; min_overlap = 3)
            Test.@test MetaGraphsNext.haskey(graph, "ATCG")
            Test.@test MetaGraphsNext.haskey(graph, "TCGA")

            stats = Mycelia.Rhizomorph.get_string_graph_statistics(graph)
            Test.@test stats[:num_vertices] == 2
            Test.@test stats[:num_edges] == 1

            sources = Mycelia.Rhizomorph.find_source_strings(graph)
            sinks = Mycelia.Rhizomorph.find_sink_strings(graph)
            Test.@test "ATCG" in sources
            Test.@test "TCGA" in sinks

            path2 = joinpath(dir, "strings2.txt")
            open(path2, "w") do io
                println(io, "TCGA")
                println(io, "CGAT")
            end

            graph_multi = Mycelia.Rhizomorph.build_string_graph_from_files([path1, path2]; min_overlap = 3)
            datasets = Mycelia.Rhizomorph.get_all_dataset_ids(graph_multi["TCGA"])
            Test.@test length(datasets) == 2
        end
    end

    Test.@testset "String graph statistics on empty graph" begin
        graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph();
            label_type = String,
            vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
            edge_data_type = Mycelia.Rhizomorph.StringEdgeData
        )

        stats = Mycelia.Rhizomorph.get_string_graph_statistics(graph)
        Test.@test stats[:num_vertices] == 0
        Test.@test stats[:num_edges] == 0
        Test.@test stats[:mean_overlap_length] === nothing
    end

    Test.@testset "String graph error cases" begin
        mktempdir() do dir
            empty_path = joinpath(dir, "empty.txt")
            open(empty_path, "w") do io
                println(io, "")
            end
            Test.@test_throws ErrorException Mycelia.Rhizomorph.build_string_graph_from_file(empty_path)
        end

        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_string_graph_from_files(String[])
    end
end
