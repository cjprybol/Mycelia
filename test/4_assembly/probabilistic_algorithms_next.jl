# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/probabilistic_algorithms_next.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/probabilistic_algorithms_next.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import Kmers

Test.@testset "Probabilistic Algorithms Next-Generation Tests" begin
    # Helper function to create a simple test graph
    function create_test_graph()
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type = String,
            vertex_data_type = Any,
            edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
            weight_function = Mycelia.Rhizomorph.edge_data_weight,
            default_weight = 0.0
        )

        # Add vertices: ATC -> TCG -> CGA
        graph["ATC"] = nothing
        graph["TCG"] = nothing
        graph["CGA"] = nothing

        # Add edges with different coverage (weight is calculated automatically)
        # High coverage edge (3 observations)
        graph["ATC", "TCG"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            3.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)

        # Low coverage edge (1 observation)
        graph["TCG", "CGA"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            1.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)

        return graph
    end

    inverse_edge_data_weight(edge_data) = 1.0 / edge_data.weight

    Test.@testset "WalkStep and GraphPath Construction" begin
        # Test WalkStep creation
        step = Mycelia.Rhizomorph.WalkStep("ATC", Mycelia.Rhizomorph.Forward, 0.8, 0.8)
        Test.@test step.vertex_label == "ATC"
        Test.@test step.strand == Mycelia.Rhizomorph.Forward
        Test.@test step.probability == 0.8
        Test.@test step.cumulative_probability == 0.8

        # Test GraphPath creation
        steps = [
            Mycelia.Rhizomorph.WalkStep("ATC", Mycelia.Rhizomorph.Forward, 1.0, 1.0),
            Mycelia.Rhizomorph.WalkStep("TCG", Mycelia.Rhizomorph.Forward, 0.9, 0.9),
            Mycelia.Rhizomorph.WalkStep("CGA", Mycelia.Rhizomorph.Forward, 0.7, 0.63)
        ]

        path = Mycelia.Rhizomorph.GraphPath(steps)
        Test.@test path.steps == steps
        Test.@test path.total_probability == 0.63
        graph = create_test_graph()
        sequence = Mycelia.Rhizomorph.path_to_sequence(path, graph)
        Test.@test sequence isa String
        Test.@test length(sequence) > 0
    end

    Test.@testset "Probabilistic Walk" begin
        graph = create_test_graph()

        # Test basic walk
        path = Mycelia.Rhizomorph.probabilistic_walk_next(graph, "ATC", 10; seed = 42)
        Test.@test path isa Mycelia.Rhizomorph.GraphPath
        Test.@test !isempty(path.steps)
        Test.@test first(path.steps).vertex_label == "ATC"
        Test.@test first(path.steps).strand == Mycelia.Rhizomorph.Forward
        Test.@test first(path.steps).probability == 1.0

        # Test that path respects max_steps
        short_path = Mycelia.Rhizomorph.probabilistic_walk_next(graph, "ATC", 1; seed = 42)
        Test.@test length(short_path.steps) <= 2  # Start vertex + 1 step

        # Test with non-existent start vertex
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.probabilistic_walk_next(graph, "XYZ", 5)

        # Test reproducibility with same seed
        path1 = Mycelia.Rhizomorph.probabilistic_walk_next(graph, "ATC", 5; seed = 123)
        path2 = Mycelia.Rhizomorph.probabilistic_walk_next(graph, "ATC", 5; seed = 123)
        sequence1 = Mycelia.Rhizomorph.path_to_sequence(path1, graph)
        sequence2 = Mycelia.Rhizomorph.path_to_sequence(path2, graph)
        Test.@test sequence1 == sequence2
    end

    Test.@testset "Maximum Weight Walk" begin
        graph = create_test_graph()

        # Test maximum weight walk (should prefer high-weight edges)
        path = Mycelia.Rhizomorph.maximum_weight_walk_next(graph, "ATC", 10)
        Test.@test path isa Mycelia.Rhizomorph.GraphPath
        Test.@test !isempty(path.steps)
        Test.@test first(path.steps).vertex_label == "ATC"

        # With our test graph, should follow ATC -> TCG (weight 3.0) -> CGA (weight 1.0)
        if length(path.steps) >= 2
            Test.@test path.steps[2].vertex_label == "TCG"
        end
        if length(path.steps) >= 3
            Test.@test path.steps[3].vertex_label == "CGA"
        end

        # Test custom weight function
        custom_weight_path = Mycelia.Rhizomorph.maximum_weight_walk_next(
            graph, "ATC", 10;
            weight_function = inverse_edge_data_weight
        )
        Test.@test custom_weight_path isa Mycelia.Rhizomorph.GraphPath

        # Test with non-existent start vertex
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.maximum_weight_walk_next(graph, "XYZ", 5)
    end

    Test.@testset "Shortest Probability Path" begin
        graph = create_test_graph()

        # Test finding path between existing vertices
        path = Mycelia.Rhizomorph.shortest_probability_path_next(graph, "ATC", "CGA")
        Test.@test path isa Mycelia.Rhizomorph.GraphPath
        Test.@test !isempty(path.steps)
        Test.@test first(path.steps).vertex_label == "ATC"
        Test.@test last(path.steps).vertex_label == "CGA"

        # Test path to same vertex
        same_path = Mycelia.Rhizomorph.shortest_probability_path_next(graph, "ATC", "ATC")
        Test.@test same_path isa Mycelia.Rhizomorph.GraphPath
        Test.@test length(same_path.steps) == 1
        Test.@test same_path.steps[1].vertex_label == "ATC"

        # Test path to unreachable vertex (reverse direction)
        no_path = Mycelia.Rhizomorph.shortest_probability_path_next(graph, "CGA", "ATC")
        Test.@test no_path === nothing

        # Test with non-existent vertices
        Test.@test Mycelia.Rhizomorph.shortest_probability_path_next(graph, "XYZ", "ATC") ===
                   nothing
        Test.@test Mycelia.Rhizomorph.shortest_probability_path_next(graph, "ATC", "XYZ") ===
                   nothing
    end

    Test.@testset "Integration with Real Graph" begin
        # Create a real k-mer graph from sequences
        seq1 = FASTX.FASTA.Record("test1", "ATCGATCG")
        seq2 = FASTX.FASTA.Record("test2", "TCGATCGA")
        observations = [seq1, seq2]

        kmer_type = Kmers.DNAKmer{3}
        base_graph = Mycelia.Rhizomorph.build_kmer_graph(
            observations, 3; dataset_id = "test", mode = :singlestrand)
        graph = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(base_graph)

        if !isempty(MetaGraphsNext.labels(graph))
            start_vertex = first(MetaGraphsNext.labels(graph))

            # Test probabilistic walk on real graph
            prob_path = Mycelia.Rhizomorph.probabilistic_walk_next(graph, start_vertex, 5; seed = 42)
            Test.@test prob_path isa Mycelia.Rhizomorph.GraphPath
            Test.@test !isempty(prob_path.steps)
            Test.@test prob_path.total_probability > 0
            sequence = Mycelia.Rhizomorph.path_to_sequence(prob_path, graph)
            Test.@test !isempty(sequence)

            # Test maximum weight walk
            max_path = Mycelia.Rhizomorph.maximum_weight_walk_next(graph, start_vertex, 5)
            Test.@test max_path isa Mycelia.Rhizomorph.GraphPath
            Test.@test !isempty(max_path.steps)

            # Test shortest path (if there are multiple vertices)
            labels = collect(MetaGraphsNext.labels(graph))
            if length(labels) >= 2
                source = labels[1]
                target = labels[2]
                short_path = Mycelia.Rhizomorph.shortest_probability_path_next(graph, source, target)
                # May be nothing if no path exists, which is valid
                if short_path !== nothing
                    Test.@test short_path isa Mycelia.Rhizomorph.GraphPath
                    Test.@test first(short_path.steps).vertex_label == source
                    Test.@test last(short_path.steps).vertex_label == target
                end
            end
        end
    end

    Test.@testset "Strand Awareness in Algorithms" begin
        # Create a graph with mixed strand orientations
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type = String,
            vertex_data_type = Any,
            edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
            weight_function = Mycelia.Rhizomorph.edge_data_weight,
            default_weight = 0.0
        )

        graph["ATC"] = nothing
        graph["GAT"] = nothing

        # Add edge requiring strand compatibility
        graph["ATC", "GAT"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            2.0, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse)

        # Test that algorithms respect strand constraints
        path = Mycelia.Rhizomorph.probabilistic_walk_next(graph, "ATC", 5; seed = 42)
        Test.@test path isa Mycelia.Rhizomorph.GraphPath
        Test.@test first(path.steps).strand == Mycelia.Rhizomorph.Forward

        # If path continues to second vertex, check strand consistency
        if length(path.steps) >= 2
            Test.@test path.steps[2].strand == Mycelia.Rhizomorph.Reverse
        end
    end

    Test.@testset "Edge Cases and Error Handling" begin
        # Empty graph
        empty_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type = String,
            vertex_data_type = Any,
            edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData
        )

        # Should handle empty graphs gracefully
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.probabilistic_walk_next(empty_graph, "ATC", 5)
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.maximum_weight_walk_next(empty_graph, "ATC", 5)

        # Single vertex graph
        single_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type = String,
            vertex_data_type = Any,
            edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData
        )
        single_graph["ATC"] = nothing

        # Should work with single vertex
        single_path = Mycelia.Rhizomorph.probabilistic_walk_next(single_graph, "ATC", 5; seed = 42)
        Test.@test length(single_path.steps) == 1
        Test.@test single_path.steps[1].vertex_label == "ATC"

        # Shortest path in single vertex graph
        self_path = Mycelia.Rhizomorph.shortest_probability_path_next(single_graph, "ATC", "ATC")
        Test.@test self_path isa Mycelia.Rhizomorph.GraphPath
        Test.@test length(self_path.steps) == 1
    end
end
