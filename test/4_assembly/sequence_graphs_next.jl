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

Test.@testset "Next-Generation Sequence Graphs Tests" begin
    
    Test.@testset "Type-Stable Metadata Structures" begin
        # Test StrandOrientation enum
        Test.@test Mycelia.Forward isa Mycelia.StrandOrientation
        Test.@test Mycelia.Reverse isa Mycelia.StrandOrientation
        Test.@test Bool(Mycelia.Forward) == true
        Test.@test Bool(Mycelia.Reverse) == false
        
        # Test GraphMode enum
        Test.@test Mycelia.SingleStrand isa Mycelia.GraphMode
        Test.@test Mycelia.DoubleStrand isa Mycelia.GraphMode
        
        # Test KmerVertexData construction
        vertex_data = Mycelia.KmerVertexData("ATCG")
        Test.@test vertex_data.canonical_kmer == "ATCG"
        Test.@test isempty(vertex_data.coverage)
        Test.@test vertex_data.coverage isa Vector{Tuple{Int, Int, Mycelia.StrandOrientation}}
        
        # Test adding coverage data with strand information
        push!(vertex_data.coverage, (1, 1, Mycelia.Forward))
        push!(vertex_data.coverage, (2, 5, Mycelia.Reverse))
        Test.@test length(vertex_data.coverage) == 2
        Test.@test vertex_data.coverage[1] == (1, 1, Mycelia.Forward)
        Test.@test vertex_data.coverage[2] == (2, 5, Mycelia.Reverse)
        
        # Test KmerEdgeData construction
        edge_data = Mycelia.KmerEdgeData(Mycelia.Forward, Mycelia.Reverse)
        Test.@test isempty(edge_data.coverage)
        Test.@test edge_data.weight == 0.0
        Test.@test edge_data.src_strand == Mycelia.Forward
        Test.@test edge_data.dst_strand == Mycelia.Reverse
        
        # Test edge data with coverage
        coverage = [((1, 1, Mycelia.Forward), (1, 2, Mycelia.Forward)), 
                   ((2, 5, Mycelia.Reverse), (2, 6, Mycelia.Reverse))]
        edge_data_with_coverage = Mycelia.KmerEdgeData(coverage, Mycelia.Forward, Mycelia.Forward)
        Test.@test length(edge_data_with_coverage.coverage) == 2
        Test.@test edge_data_with_coverage.weight == 2.0  # weight equals coverage count
        Test.@test edge_data_with_coverage.src_strand == Mycelia.Forward
        Test.@test edge_data_with_coverage.dst_strand == Mycelia.Forward
    end
    
    Test.@testset "Next-Generation K-mer Graph Construction" begin
        # Create test sequences
        seq1 = FASTX.FASTA.Record("test1", "ATCGATCG")
        seq2 = FASTX.FASTA.Record("test2", "TCGATCGA")
        observations = [seq1, seq2]
        
        # Test with 3-mers
        kmer_type = BioSequences.DNAKmer{3}
        graph = Mycelia.build_kmer_graph_next(kmer_type, observations)
        
        Test.@test graph isa MetaGraphsNext.MetaGraph
        Test.@test !isempty(MetaGraphsNext.labels(graph))
        
        # Check that vertices exist for expected k-mers
        expected_kmers = ["ATC", "TCG", "CGA", "GAT", "ATG", "GTC"]  # some from sequences
        found_kmers = collect(MetaGraphsNext.labels(graph))
        
        # At least some expected k-mers should be present
        Test.@test !isempty(intersect(expected_kmers, found_kmers))
        
        # Test vertex metadata
        for label in found_kmers
            vertex_data = graph[label]
            Test.@test vertex_data isa Mycelia.KmerVertexData
            Test.@test vertex_data.coverage isa Vector{Tuple{Int, Int, Bool}}
        end
        
        # Test edge metadata (if edges exist)
        for edge in MetaGraphsNext.edge_labels(graph)
            edge_data = graph[edge...]
            Test.@test edge_data isa Mycelia.KmerEdgeData
            Test.@test edge_data.weight >= 0.0
        end
    end
    
    Test.@testset "Graph Modes: SingleStrand vs DoubleStrand" begin
        seq = FASTX.FASTA.Record("test", "ATCG")
        observations = [seq]
        kmer_type = BioSequences.DNAKmer{3}
        
        # Test DoubleStrand graph (default) - uses canonical k-mers
        double_strand_graph = Mycelia.build_kmer_graph_next(kmer_type, observations; graph_mode=Mycelia.DoubleStrand)
        double_strand_labels = collect(MetaGraphsNext.labels(double_strand_graph))
        
        # Test SingleStrand graph - uses k-mers as-is
        single_strand_graph = Mycelia.build_kmer_graph_next(kmer_type, observations; graph_mode=Mycelia.SingleStrand)
        single_strand_labels = collect(MetaGraphsNext.labels(single_strand_graph))
        
        # Both should be valid graphs
        Test.@test double_strand_graph isa MetaGraphsNext.MetaGraph
        Test.@test single_strand_graph isa MetaGraphsNext.MetaGraph
        
        # In DoubleStrand mode, vertices should be canonical k-mers
        # In SingleStrand mode, vertices should be observed k-mers
        Test.@test !isempty(double_strand_labels)
        Test.@test !isempty(single_strand_labels)
        
        # Test that vertices contain canonical k-mers in DoubleStrand mode
        for label in double_strand_labels
            vertex_data = double_strand_graph[label]
            Test.@test vertex_data isa Mycelia.KmerVertexData
            Test.@test vertex_data.canonical_kmer == label
        end
    end
    
    Test.@testset "Empty and Edge Cases" begin
        kmer_type = BioSequences.DNAKmer{3}
        
        # Test with empty observations
        empty_graph = Mycelia.build_kmer_graph_next(kmer_type, FASTX.FASTA.Record[])
        Test.@test empty_graph isa MetaGraphsNext.MetaGraph
        Test.@test isempty(MetaGraphsNext.labels(empty_graph))
        
        # Test with sequence shorter than k
        short_seq = FASTX.FASTA.Record("short", "AT")  # length 2 < k=3
        short_graph = Mycelia.build_kmer_graph_next(kmer_type, [short_seq])
        Test.@test short_graph isa MetaGraphsNext.MetaGraph
        # Should warn but not crash
    end
    
    Test.@testset "Strand-Aware Coverage Tracking" begin
        # Create sequences with overlapping k-mers that test strand awareness
        seq1 = FASTX.FASTA.Record("seq1", "ATCGATCG")  # Forward: ATC, TCG, CGA, GAT, ATC, TCG
        seq2 = FASTX.FASTA.Record("seq2", "CGATCG")    # Forward: CGA, GAT, ATC, TCG  
        observations = [seq1, seq2]
        
        kmer_type = BioSequences.DNAKmer{3}
        graph = Mycelia.build_kmer_graph_next(kmer_type, observations; graph_mode=Mycelia.DoubleStrand)
        
        # All vertices should represent canonical k-mers
        for label in MetaGraphsNext.labels(graph)
            vertex_data = graph[label]
            Test.@test vertex_data isa Mycelia.KmerVertexData
            Test.@test vertex_data.canonical_kmer == label
            
            # Check that coverage includes strand information
            for (obs_id, pos, strand) in vertex_data.coverage
                Test.@test obs_id in [1, 2]  # From one of our sequences
                Test.@test pos >= 1
                Test.@test strand isa Mycelia.StrandOrientation
            end
        end
        
        # Check that edges have strand-aware metadata
        for edge_label in MetaGraphsNext.edge_labels(graph)
            if !isempty(edge_label)
                edge_data = graph[edge_label...]
                Test.@test edge_data isa Mycelia.KmerEdgeData
                Test.@test edge_data.src_strand isa Mycelia.StrandOrientation
                Test.@test edge_data.dst_strand isa Mycelia.StrandOrientation
                Test.@test edge_data.weight >= 0.0
            end
        end
    end
end

Test.@testset "Legacy Compatibility Layer" begin
    Test.@testset "Legacy Graph Detection" begin
        # This would test actual MetaGraphs if we had a legacy graph
        # For now, test the detection function with mock inputs
        
        # Test with next-generation graph
        next_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        
        Test.@test !Mycelia.is_legacy_graph(next_graph)
        Test.@test Mycelia.ensure_next_graph(next_graph) === next_graph
    end
    
    Test.@testset "Ensure Next Graph Function" begin
        # Test that ensure_next_graph passes through next-generation graphs unchanged
        next_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        
        result = Mycelia.ensure_next_graph(next_graph)
        Test.@test result === next_graph  # Should be identical object
        Test.@test result isa MetaGraphsNext.MetaGraph
    end
end

Test.@testset "Integration with Existing String Graphs" begin
    # Test that both string graphs and k-mer graphs use compatible MetaGraphsNext structures
    
    # Create a string graph
    string_graph = Mycelia.string_to_ngram_graph(s="ABCABC", n=2)
    Test.@test string_graph isa MetaGraphsNext.MetaGraph
    
    # Create a k-mer graph
    seq = FASTX.FASTA.Record("test", "ATCGATCG")
    kmer_type = BioSequences.DNAKmer{3}
    kmer_graph = Mycelia.build_kmer_graph_next(kmer_type, [seq])
    Test.@test kmer_graph isa MetaGraphsNext.MetaGraph
    
    # Both should use MetaGraphsNext - this confirms API unification
    Test.@test typeof(string_graph).name.module == MetaGraphsNext
    Test.@test typeof(kmer_graph).name.module == MetaGraphsNext
end
