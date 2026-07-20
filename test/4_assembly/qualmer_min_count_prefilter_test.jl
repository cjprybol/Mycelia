import Test
import FASTX
import Kmers
import Mycelia

Test.@testset "Qualmer minimum-count prefilter" begin
    records = [
        FASTX.FASTQ.Record("shared_1", "AAAC", "IIII"),
        FASTX.FASTQ.Record("shared_2", "AAAG", "IIII"),
    ]

    Test.@testset "default preserves all observations" begin
        default_graph = Mycelia.Rhizomorph.build_qualmer_graph(records, 3)
        explicit_graph = Mycelia.Rhizomorph.build_qualmer_graph(records, 3; min_count = 1)

        Test.@test Set(Mycelia.Rhizomorph.labels(default_graph)) ==
                   Set(Mycelia.Rhizomorph.labels(explicit_graph))
        Test.@test Mycelia.Rhizomorph.vertex_count(default_graph) == 3
        Test.@test Mycelia.Rhizomorph.edge_count(default_graph) == 2
    end

    Test.@testset "filters singleton vertices and incident edges before allocation" begin
        graph = Mycelia.Rhizomorph.build_qualmer_graph(records, 3; min_count = 2)
        shared = Kmers.DNAKmer{3}("AAA")

        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 1
        Test.@test Mycelia.Rhizomorph.edge_count(graph) == 0
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, shared)
        Test.@test Mycelia.Rhizomorph.get_vertex_observation_count(graph, shared) == 2
    end

    Test.@testset "applies before strand conversion" begin
        expected_vertex_counts = Dict(:doublestrand => 2, :canonical => 1)
        for mode in keys(expected_vertex_counts)
            graph = Mycelia.Rhizomorph.build_qualmer_graph(
                records, 3; mode = mode, min_count = 2)
            Test.@test Mycelia.Rhizomorph.vertex_count(graph) == expected_vertex_counts[mode]
            Test.@test Mycelia.Rhizomorph.edge_count(graph) == 0
        end
    end

    Test.@testset "rejects invalid thresholds" begin
        error = try
            Mycelia.Rhizomorph.build_qualmer_graph(records, 3; min_count = 0)
            nothing
        catch exception
            exception
        end
        Test.@test error isa ArgumentError
        Test.@test occursin("min_count must be at least 1", sprint(showerror, error))
    end
end
