import Test
import Mycelia

Test.@testset "N-gram Graph Utilities" begin
    Test.@testset "Build from file and stats" begin
        mktempdir() do dir
            path = joinpath(dir, "doc1.txt")
            open(path, "w") do io
                println(io, "ABCD")
                println(io, "BCDA")
            end

            graph = Mycelia.Rhizomorph.build_ngram_graph_from_file(path, 2)
            stats = Mycelia.Rhizomorph.get_ngram_statistics(graph)

            Test.@test stats[:num_vertices] == 4
            Test.@test stats[:num_edges] == 3
            Test.@test stats[:n] == 2
            Test.@test stats[:most_common_ngram] in ["BC", "CD"]
        end
    end

    Test.@testset "Multi-file comparisons" begin
        mktempdir() do dir
            path1 = joinpath(dir, "doc1.txt")
            path2 = joinpath(dir, "doc2.txt")
            open(path1, "w") do io
                println(io, "ABCD")
            end
            open(path2, "w") do io
                println(io, "BCXY")
            end

            graph = Mycelia.Rhizomorph.build_ngram_graph_from_files([path1, path2], 2)
            dataset_id1 = Mycelia.Rhizomorph.get_dataset_id_from_file(path1)
            dataset_id2 = Mycelia.Rhizomorph.get_dataset_id_from_file(path2)

            shared = Mycelia.Rhizomorph.find_shared_ngrams(graph, [
                dataset_id1, dataset_id2])
            Test.@test "BC" in shared

            unique_doc1 = Mycelia.Rhizomorph.find_unique_ngrams(graph, dataset_id1)
            Test.@test "AB" in unique_doc1
            Test.@test "CD" in unique_doc1

            common = Mycelia.Rhizomorph.find_high_coverage_ngrams(graph, 2)
            Test.@test "BC" in common
        end
    end

    Test.@testset "N-gram graph error cases" begin
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_ngram_graph_from_files(String[], 2)
    end
end
