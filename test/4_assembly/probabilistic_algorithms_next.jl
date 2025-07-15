import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX

Test.@testset "Probabilistic Algorithms Next-Generation Tests" begin
    
    # Helper function to create a simple test graph
    function create_test_graph()
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData,
            weight_function=edge_data -> edge_data.weight,
            default_weight=0.0
        )
        
        # Add vertices: ATC -> TCG -> CGA
        graph["ATC"] = Mycelia.KmerVertexData("ATC")
        graph["TCG"] = Mycelia.KmerVertexData("TCG")  
        graph["CGA"] = Mycelia.KmerVertexData("CGA")
        
        # Add edges with different weights
        graph["ATC", "TCG"] = Mycelia.KmerEdgeData(
            [(Tuple{Int,Int,Mycelia.StrandOrientation}, Tuple{Int,Int,Mycelia.StrandOrientation})[]],
            3.0,  # High weight
            Mycelia.Forward, Mycelia.Forward
        )
        
        graph["TCG", "CGA"] = Mycelia.KmerEdgeData(
            [(Tuple{Int,Int,Mycelia.StrandOrientation}, Tuple{Int,Int,Mycelia.StrandOrientation})[]],
            1.0,  # Lower weight
            Mycelia.Forward, Mycelia.Forward
        )
        
        return graph
    end
    
    Test.@testset "WalkStep and GraphPath Construction" begin
        # Test WalkStep creation
        step = Mycelia.WalkStep("ATC", Mycelia.Forward, 0.8, 0.8)
        Test.@test step.vertex_label == "ATC"
        Test.@test step.strand == Mycelia.Forward
        Test.@test step.probability == 0.8
        Test.@test step.cumulative_probability == 0.8
        
        # Test GraphPath creation
        steps = [
            Mycelia.WalkStep("ATC", Mycelia.Forward, 1.0, 1.0),
            Mycelia.WalkStep("TCG", Mycelia.Forward, 0.9, 0.9),
            Mycelia.WalkStep("CGA", Mycelia.Forward, 0.7, 0.63)
        ]
        
        path = Mycelia.GraphPath(steps)
        Test.@test path.steps == steps
        Test.@test path.total_probability == 0.63
        Test.@test path.sequence isa String
        Test.@test length(path.sequence) > 0
    end
    
    Test.@testset "Probabilistic Walk" begin
        graph = create_test_graph()
        
        # Test basic walk
        path = Mycelia.probabilistic_walk_next(graph, "ATC", 10; seed=42)
        Test.@test path isa Mycelia.GraphPath
        Test.@test !isempty(path.steps)
        Test.@test first(path.steps).vertex_label == "ATC"
        Test.@test first(path.steps).strand == Mycelia.Forward
        Test.@test first(path.steps).probability == 1.0
        
        # Test that path respects max_steps
        short_path = Mycelia.probabilistic_walk_next(graph, "ATC", 1; seed=42)
        Test.@test length(short_path.steps) <= 2  # Start vertex + 1 step
        
        # Test with non-existent start vertex
        Test.@test_throws ArgumentError Mycelia.probabilistic_walk_next(graph, "XYZ", 5)
        
        # Test reproducibility with same seed
        path1 = Mycelia.probabilistic_walk_next(graph, "ATC", 5; seed=123)
        path2 = Mycelia.probabilistic_walk_next(graph, "ATC", 5; seed=123)
        Test.@test path1.sequence == path2.sequence
    end
    
    Test.@testset "Maximum Weight Walk" begin
        graph = create_test_graph()
        
        # Test maximum weight walk (should prefer high-weight edges)
        path = Mycelia.maximum_weight_walk_next(graph, "ATC", 10)
        Test.@test path isa Mycelia.GraphPath
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
        custom_weight_path = Mycelia.maximum_weight_walk_next(
            graph, "ATC", 10;
            weight_function = edge_data -> 1.0 / edge_data.weight  # Inverse weight
        )
        Test.@test custom_weight_path isa Mycelia.GraphPath
        
        # Test with non-existent start vertex
        Test.@test_throws ArgumentError Mycelia.maximum_weight_walk_next(graph, "XYZ", 5)
    end
    
    Test.@testset "Shortest Probability Path" begin
        graph = create_test_graph()
        
        # Test finding path between existing vertices
        path = Mycelia.shortest_probability_path_next(graph, "ATC", "CGA")
        Test.@test path isa Mycelia.GraphPath
        Test.@test !isempty(path.steps)
        Test.@test first(path.steps).vertex_label == "ATC"
        Test.@test last(path.steps).vertex_label == "CGA"
        
        # Test path to same vertex
        same_path = Mycelia.shortest_probability_path_next(graph, "ATC", "ATC")
        Test.@test same_path isa Mycelia.GraphPath
        Test.@test length(same_path.steps) == 1
        Test.@test same_path.steps[1].vertex_label == "ATC"
        
        # Test path to unreachable vertex (reverse direction)
        no_path = Mycelia.shortest_probability_path_next(graph, "CGA", "ATC")
        Test.@test no_path === nothing
        
        # Test with non-existent vertices
        Test.@test Mycelia.shortest_probability_path_next(graph, "XYZ", "ATC") === nothing
        Test.@test Mycelia.shortest_probability_path_next(graph, "ATC", "XYZ") === nothing
    end
    
    Test.@testset "Integration with Real Graph" begin
        # Create a real k-mer graph from sequences
        seq1 = FASTX.FASTA.Record("test1", "ATCGATCG")
        seq2 = FASTX.FASTA.Record("test2", "TCGATCGA")
        observations = [seq1, seq2]
        
        kmer_type = BioSequences.DNAKmer{3}
        graph = Mycelia.build_kmer_graph_next(kmer_type, observations)
        
        if !isempty(MetaGraphsNext.labels(graph))
            start_vertex = first(MetaGraphsNext.labels(graph))
            
            # Test probabilistic walk on real graph
            prob_path = Mycelia.probabilistic_walk_next(graph, start_vertex, 5; seed=42)
            Test.@test prob_path isa Mycelia.GraphPath
            Test.@test !isempty(prob_path.steps)
            Test.@test prob_path.total_probability > 0
            Test.@test !isempty(prob_path.sequence)
            
            # Test maximum weight walk
            max_path = Mycelia.maximum_weight_walk_next(graph, start_vertex, 5)
            Test.@test max_path isa Mycelia.GraphPath
            Test.@test !isempty(max_path.steps)
            
            # Test shortest path (if there are multiple vertices)
            labels = collect(MetaGraphsNext.labels(graph))
            if length(labels) >= 2
                source = labels[1]
                target = labels[2]
                short_path = Mycelia.shortest_probability_path_next(graph, source, target)
                # May be nothing if no path exists, which is valid
                if short_path !== nothing
                    Test.@test short_path isa Mycelia.GraphPath
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
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData,
            weight_function=edge_data -> edge_data.weight,
            default_weight=0.0
        )
        
        graph["ATC"] = Mycelia.KmerVertexData("ATC")
        graph["GAT"] = Mycelia.KmerVertexData("GAT")  # Reverse complement of ATC
        
        # Add edge requiring strand compatibility
        graph["ATC", "GAT"] = Mycelia.KmerEdgeData(
            [(Tuple{Int,Int,Mycelia.StrandOrientation}, Tuple{Int,Int,Mycelia.StrandOrientation})[]],
            2.0,
            Mycelia.Forward, Mycelia.Reverse  # Forward ATC -> Reverse GAT
        )
        
        # Test that algorithms respect strand constraints
        path = Mycelia.probabilistic_walk_next(graph, "ATC", 5; seed=42)
        Test.@test path isa Mycelia.GraphPath
        Test.@test first(path.steps).strand == Mycelia.Forward
        
        # If path continues to second vertex, check strand consistency
        if length(path.steps) >= 2
            Test.@test path.steps[2].strand == Mycelia.Reverse
        end
    end
    
    Test.@testset "Edge Cases and Error Handling" begin
        # Empty graph
        empty_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        
        # Should handle empty graphs gracefully
        Test.@test_throws ArgumentError Mycelia.probabilistic_walk_next(empty_graph, "ATC", 5)
        Test.@test_throws ArgumentError Mycelia.maximum_weight_walk_next(empty_graph, "ATC", 5)
        
        # Single vertex graph
        single_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        single_graph["ATC"] = Mycelia.KmerVertexData("ATC")
        
        # Should work with single vertex
        single_path = Mycelia.probabilistic_walk_next(single_graph, "ATC", 5; seed=42)
        Test.@test length(single_path.steps) == 1
        Test.@test single_path.steps[1].vertex_label == "ATC"
        
        # Shortest path in single vertex graph
        self_path = Mycelia.shortest_probability_path_next(single_graph, "ATC", "ATC")
        Test.@test self_path isa Mycelia.GraphPath
        Test.@test length(self_path.steps) == 1
    end
end
