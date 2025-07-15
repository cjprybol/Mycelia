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

Test.@testset "Advanced Graph Algorithms Next-Generation Tests" begin
    
    # Helper function to create test graph with known structure
    function create_complex_test_graph()
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData,
            weight_function=edge_data -> edge_data.weight,
            default_weight=0.0
        )
        
        # Create a graph with linear path + bubble + repeat
        # Linear: A -> B -> C
        # Bubble: C -> D -> F, C -> E -> F  
        # Continue: F -> G -> H
        vertices = ["AAA", "AAT", "ATC", "TCG", "CGA", "GAT", "ATG", "TGC"]
        
        for v in vertices
            graph[v] = Mycelia.KmerVertexData(v)
        end
        
        # Linear path
        graph["AAA", "AAT"] = Mycelia.KmerEdgeData([], 1.0, Mycelia.Forward, Mycelia.Forward)
        graph["AAT", "ATC"] = Mycelia.KmerEdgeData([], 1.0, Mycelia.Forward, Mycelia.Forward)
        
        # Bubble structure
        graph["ATC", "TCG"] = Mycelia.KmerEdgeData([], 0.8, Mycelia.Forward, Mycelia.Forward)  # Path 1
        graph["ATC", "CGA"] = Mycelia.KmerEdgeData([], 0.6, Mycelia.Forward, Mycelia.Forward)  # Path 2
        graph["TCG", "GAT"] = Mycelia.KmerEdgeData([], 0.9, Mycelia.Forward, Mycelia.Forward)
        graph["CGA", "GAT"] = Mycelia.KmerEdgeData([], 0.7, Mycelia.Forward, Mycelia.Forward)
        
        # Continue linear
        graph["GAT", "ATG"] = Mycelia.KmerEdgeData([], 1.0, Mycelia.Forward, Mycelia.Forward)
        graph["ATG", "TGC"] = Mycelia.KmerEdgeData([], 1.0, Mycelia.Forward, Mycelia.Forward)
        
        return graph
    end
    
    function create_eulerian_test_graph()
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData,
            weight_function=edge_data -> edge_data.weight,
            default_weight=0.0
        )
        
        # Create simple cycle: A -> B -> C -> A
        vertices = ["ATC", "TCG", "CGA"]
        for v in vertices
            graph[v] = Mycelia.KmerVertexData(v)
        end
        
        graph["ATC", "TCG"] = Mycelia.KmerEdgeData([], 1.0, Mycelia.Forward, Mycelia.Forward)
        graph["TCG", "CGA"] = Mycelia.KmerEdgeData([], 1.0, Mycelia.Forward, Mycelia.Forward)
        graph["CGA", "ATC"] = Mycelia.KmerEdgeData([], 1.0, Mycelia.Forward, Mycelia.Forward)
        
        return graph
    end
    
    Test.@testset "BubbleStructure Construction" begin
        bubble = Mycelia.BubbleStructure("entry", "exit", 
                                       ["path1a", "path1b"], 
                                       ["path2a", "path2b"],
                                       5, 3, 0.25)
        
        Test.@test bubble.entry_vertex == "entry"
        Test.@test bubble.exit_vertex == "exit"
        Test.@test bubble.path1 == ["path1a", "path1b"]
        Test.@test bubble.path2 == ["path2a", "path2b"]
        Test.@test bubble.path1_support == 5
        Test.@test bubble.path2_support == 3
        Test.@test bubble.complexity_score == 0.25
    end
    
    Test.@testset "RepeatRegion Construction" begin
        repeat_region = Mycelia.RepeatRegion(
            ["repeat1", "repeat2"],
            [("incoming1", "repeat1")],
            [("repeat2", "outgoing1")],
            2.5, :tandem, 0.8
        )
        
        Test.@test repeat_region.repeat_vertices == ["repeat1", "repeat2"]
        Test.@test repeat_region.copy_number_estimate == 2.5
        Test.@test repeat_region.repeat_type == :tandem
        Test.@test repeat_region.confidence == 0.8
        
        # Test invalid repeat type
        Test.@test_throws AssertionError Mycelia.RepeatRegion(
            ["repeat1"], [], [], 1.0, :invalid, 0.5
        )
        
        # Test invalid confidence
        Test.@test_throws AssertionError Mycelia.RepeatRegion(
            ["repeat1"], [], [], 1.0, :tandem, 1.5
        )
    end
    
    Test.@testset "ContigPath Construction" begin
        contig = Mycelia.ContigPath(
            ["kmer1", "kmer2", "kmer3"],
            "ATCGAT",
            [10.0, 12.0, 8.0]
        )
        
        Test.@test contig.vertices == ["kmer1", "kmer2", "kmer3"]
        Test.@test contig.sequence == "ATCGAT"
        Test.@test contig.coverage_profile == [10.0, 12.0, 8.0]
        Test.@test contig.length == 6
    end
    
    Test.@testset "Degree Calculation" begin
        graph = create_complex_test_graph()
        in_degrees, out_degrees = Mycelia.calculate_degrees(graph)
        
        # Check specific vertex degrees
        Test.@test out_degrees["AAA"] == 1  # Start vertex
        Test.@test in_degrees["TGC"] == 1   # End vertex
        Test.@test out_degrees["ATC"] == 2  # Bubble entry (branches out)
        Test.@test in_degrees["GAT"] == 2   # Bubble exit (merges in)
        
        # All vertices should have entries
        for vertex in MetaGraphsNext.labels(graph)
            Test.@test haskey(in_degrees, vertex)
            Test.@test haskey(out_degrees, vertex)
        end
    end
    
    Test.@testset "Eulerian Path Detection" begin
        # Test Eulerian cycle
        cycle_graph = create_eulerian_test_graph()
        eulerian_paths = Mycelia.find_eulerian_paths_next(cycle_graph)
        
        Test.@test !isempty(eulerian_paths)
        
        # Each path should visit all vertices
        for path in eulerian_paths
            Test.@test length(path) == 4  # 3 vertices + return to start
        end
        
        # Test non-Eulerian graph
        non_eulerian = create_complex_test_graph()
        non_eulerian_paths = Mycelia.find_eulerian_paths_next(non_eulerian)
        # This graph has unbalanced degrees, so no Eulerian path
        
        # Test empty graph
        empty_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        empty_paths = Mycelia.find_eulerian_paths_next(empty_graph)
        Test.@test isempty(empty_paths)
    end
    
    Test.@testset "Bubble Detection" begin
        graph = create_complex_test_graph()
        bubbles = Mycelia.detect_bubbles_next(graph, min_bubble_length=1, max_bubble_length=5)
        
        Test.@test !isempty(bubbles)
        
        # Should detect the bubble from ATC -> GAT
        bubble_found = false
        for bubble in bubbles
            if bubble.entry_vertex == "ATC" && bubble.exit_vertex == "GAT"
                bubble_found = true
                Test.@test length(bubble.path1) >= 1
                Test.@test length(bubble.path2) >= 1
                Test.@test bubble.path1 != bubble.path2
                Test.@test bubble.complexity_score >= 0
            end
        end
        Test.@test bubble_found
        
        # Test with restrictive parameters
        no_bubbles = Mycelia.detect_bubbles_next(graph, min_bubble_length=10, max_bubble_length=5)
        Test.@test isempty(no_bubbles)
        
        # Test empty graph
        empty_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        empty_bubbles = Mycelia.detect_bubbles_next(empty_graph)
        Test.@test isempty(empty_bubbles)
    end
    
    Test.@testset "Path Support Calculation" begin
        graph = create_complex_test_graph()
        
        # Test path support calculation
        path = ["AAA", "AAT", "ATC"]
        support = Mycelia.calculate_path_support(graph, path)
        Test.@test support >= 0
        
        # Empty path
        empty_support = Mycelia.calculate_path_support(graph, String[])
        Test.@test empty_support == 0
        
        # Non-existent vertices
        fake_path = ["XXX", "YYY"]
        fake_support = Mycelia.calculate_path_support(graph, fake_path)
        Test.@test fake_support >= 0
    end
    
    Test.@testset "Repeat Detection" begin
        graph = create_complex_test_graph()
        repeats = Mycelia.resolve_repeats_next(graph, min_repeat_length=1)
        
        # Should be able to run without errors
        Test.@test repeats isa Vector{Mycelia.RepeatRegion}
        
        # Each repeat should be valid
        for repeat in repeats
            Test.@test !isempty(repeat.repeat_vertices)
            Test.@test repeat.copy_number_estimate > 0
            Test.@test repeat.repeat_type in [:tandem, :interspersed, :palindromic]
            Test.@test 0.0 <= repeat.confidence <= 1.0
        end
        
        # Test with high threshold (should find fewer/no repeats)
        strict_repeats = Mycelia.resolve_repeats_next(graph, min_repeat_length=20)
        Test.@test length(strict_repeats) <= length(repeats)
    end
    
    Test.@testset "Contig Finding" begin
        graph = create_complex_test_graph()
        contigs = Mycelia.find_contigs_next(graph, min_contig_length=1)
        
        Test.@test !isempty(contigs)
        
        # Each contig should be valid
        for contig in contigs
            Test.@test !isempty(contig.vertices)
            Test.@test !isempty(contig.sequence)
            Test.@test length(contig.coverage_profile) == length(contig.vertices)
            Test.@test contig.length == length(contig.sequence)
            Test.@test contig.length >= 1  # Due to min_contig_length=1
        end
        
        # Contigs should be sorted by length (descending)
        for i in 1:(length(contigs)-1)
            Test.@test contigs[i].length >= contigs[i+1].length
        end
        
        # Test with high threshold
        long_contigs = Mycelia.find_contigs_next(graph, min_contig_length=100)
        Test.@test length(long_contigs) <= length(contigs)
    end
    
    Test.@testset "Graph Simplification" begin
        graph = create_complex_test_graph()
        bubbles = Mycelia.detect_bubbles_next(graph)
        
        simplified = Mycelia.simplify_graph_next(graph, bubbles)
        
        Test.@test simplified isa typeof(graph)
        
        # Simplified graph should have same or fewer vertices
        Test.@test length(MetaGraphsNext.labels(simplified)) <= length(MetaGraphsNext.labels(graph))
        
        # Should still be a valid graph
        Test.@test length(MetaGraphsNext.labels(simplified)) >= 0
        
        # Test with empty bubbles
        no_bubbles = Mycelia.BubbleStructure[]
        unchanged = Mycelia.simplify_graph_next(graph, no_bubbles)
        # Should be essentially unchanged (except isolated vertex removal)
    end
    
    Test.@testset "Neighbor Finding Functions" begin
        graph = create_complex_test_graph()
        
        # Test outgoing neighbors
        out_neighbors = Mycelia.get_out_neighbors(graph, "ATC")
        Test.@test length(out_neighbors) == 2  # Should have 2 due to bubble
        Test.@test "TCG" in out_neighbors
        Test.@test "CGA" in out_neighbors
        
        # Test incoming neighbors
        in_neighbors = Mycelia.get_in_neighbors(graph, "GAT")
        Test.@test length(in_neighbors) == 2  # Should have 2 due to bubble convergence
        Test.@test "TCG" in in_neighbors
        Test.@test "CGA" in in_neighbors
        
        # Test vertex with no neighbors
        single_neighbors = Mycelia.get_out_neighbors(graph, "TGC")
        Test.@test isempty(single_neighbors)
        
        # Test non-existent vertex
        no_neighbors = Mycelia.get_out_neighbors(graph, "XXX")
        Test.@test isempty(no_neighbors)
    end
    
    Test.@testset "Path Validation" begin
        graph = create_complex_test_graph()
        
        # Valid path
        valid_path = ["AAA", "AAT", "ATC"]
        Test.@test Mycelia.is_valid_path(graph, valid_path)
        
        # Invalid path (missing edge)
        invalid_path = ["AAA", "TCG"]  # No direct edge
        Test.@test !Mycelia.is_valid_path(graph, invalid_path)
        
        # Single vertex (always valid)
        single_path = ["AAA"]
        Test.@test Mycelia.is_valid_path(graph, single_path)
        
        # Empty path (always valid)
        empty_path = String[]
        Test.@test Mycelia.is_valid_path(graph, empty_path)
    end
    
    Test.@testset "Sequence Generation" begin
        graph = create_complex_test_graph()
        
        # Test contig sequence generation
        path = ["AAA", "AAT", "ATC"]
        sequence = Mycelia.generate_contig_sequence(graph, path)
        
        Test.@test !isempty(sequence)
        Test.@test length(sequence) >= length(path)  # Should be at least as long as k-mer count
        Test.@test startswith(sequence, "AAA")  # Should start with first k-mer
        
        # Empty path
        empty_sequence = Mycelia.generate_contig_sequence(graph, String[])
        Test.@test empty_sequence == ""
        
        # Single k-mer
        single_sequence = Mycelia.generate_contig_sequence(graph, ["AAA"])
        Test.@test single_sequence == "AAA"
    end
    
    Test.@testset "Coverage Profile Generation" begin
        graph = create_complex_test_graph()
        
        path = ["AAA", "AAT", "ATC"]
        coverage = Mycelia.generate_coverage_profile(graph, path)
        
        Test.@test length(coverage) == length(path)
        Test.@test all(c >= 0 for c in coverage)
        
        # Empty path
        empty_coverage = Mycelia.generate_coverage_profile(graph, String[])
        Test.@test isempty(empty_coverage)
    end
    
    Test.@testset "Integration Test with Real Graph" begin
        # Create a real k-mer graph and test algorithms
        seq1 = FASTX.FASTA.Record("test1", "ATCGATCGAT")
        seq2 = FASTX.FASTA.Record("test2", "TCGATCGATC")
        observations = [seq1, seq2]
        
        kmer_type = BioSequences.DNAKmer{3}
        graph = Mycelia.build_kmer_graph_next(kmer_type, observations)
        
        if !isempty(MetaGraphsNext.labels(graph))
            # Test bubble detection
            bubbles = Mycelia.detect_bubbles_next(graph)
            Test.@test bubbles isa Vector{Mycelia.BubbleStructure}
            
            # Test repeat detection
            repeats = Mycelia.resolve_repeats_next(graph)
            Test.@test repeats isa Vector{Mycelia.RepeatRegion}
            
            # Test contig finding
            contigs = Mycelia.find_contigs_next(graph, min_contig_length=1)
            Test.@test contigs isa Vector{Mycelia.ContigPath}
            Test.@test !isempty(contigs)
            
            # Test graph simplification
            simplified = Mycelia.simplify_graph_next(graph, bubbles)
            Test.@test simplified isa typeof(graph)
            
            # Test Eulerian path finding
            paths = Mycelia.find_eulerian_paths_next(graph)
            Test.@test paths isa Vector{Vector{String}}
        end
    end
    
    Test.@testset "Edge Cases and Error Handling" begin
        # Empty graph tests
        empty_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        
        # All functions should handle empty graphs gracefully
        Test.@test isempty(Mycelia.find_eulerian_paths_next(empty_graph))
        Test.@test isempty(Mycelia.detect_bubbles_next(empty_graph))
        Test.@test isempty(Mycelia.resolve_repeats_next(empty_graph))
        Test.@test isempty(Mycelia.find_contigs_next(empty_graph))
        
        # Single vertex graph
        single_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        single_graph["ATC"] = Mycelia.KmerVertexData("ATC")
        
        # Should handle single vertex graphs
        Test.@test isempty(Mycelia.detect_bubbles_next(single_graph))
        Test.@test isempty(Mycelia.resolve_repeats_next(single_graph))
        single_contigs = Mycelia.find_contigs_next(single_graph, min_contig_length=1)
        Test.@test isempty(single_contigs)  # Single vertex doesn't form a contig
    end
end
