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
import Tempfile

Test.@testset "GFA I/O Next-Generation Tests" begin
    
    Test.@testset "GFA Writing" begin
        # Create a simple test graph
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData,
            weight_function=edge_data -> edge_data.weight,
            default_weight=0.0
        )
        
        # Add vertices
        graph["ATC"] = Mycelia.KmerVertexData("ATC")
        graph["TCG"] = Mycelia.KmerVertexData("TCG")
        
        # Add coverage to vertices
        push!(graph["ATC"].coverage, (1, 1, Mycelia.Forward))
        push!(graph["ATC"].coverage, (2, 3, Mycelia.Reverse))
        push!(graph["TCG"].coverage, (1, 2, Mycelia.Forward))
        
        # Add an edge
        graph["ATC", "TCG"] = Mycelia.KmerEdgeData(Mycelia.Forward, Mycelia.Forward)
        
        # Write to temporary file
        Tempfile.mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "test.gfa")
            result_file = Mycelia.write_gfa_next(graph, gfa_file)
            
            Test.@test result_file == gfa_file
            Test.@test isfile(gfa_file)
            
            # Read and verify content
            content = read(gfa_file, String)
            lines = split(strip(content), '\n')
            
            # Check header
            Test.@test startswith(lines[1], "H\t")
            Test.@test occursin("VN:Z:1.0", lines[1])
            
            # Check segments
            segment_lines = filter(line -> startswith(line, "S\t"), lines)
            Test.@test length(segment_lines) == 2
            
            # Check that segments contain our k-mers
            segment_content = join(segment_lines, " ")
            Test.@test occursin("ATC", segment_content)
            Test.@test occursin("TCG", segment_content)
            
            # Check coverage information
            Test.@test occursin("DP:f:2", segment_content)  # ATC has 2 coverage entries
            Test.@test occursin("DP:f:1", segment_content)  # TCG has 1 coverage entry
            
            # Check links
            link_lines = filter(line -> startswith(line, "L\t"), lines)
            Test.@test length(link_lines) == 1
            
            # Check link format and orientations
            link_line = link_lines[1]
            link_fields = split(link_line, '\t')
            Test.@test length(link_fields) == 6
            Test.@test link_fields[3] == "+"  # source orientation (Forward)
            Test.@test link_fields[5] == "+"  # destination orientation (Forward)
            Test.@test link_fields[6] == "2M"  # overlap (k-1 = 3-1 = 2)
        end
    end
    
    Test.@testset "GFA Reading" begin
        # Create a test GFA file content
        gfa_content = """
        H	VN:Z:1.0
        S	1	ATC	DP:f:2.0
        S	2	TCG	DP:f:1.0
        L	1	+	2	+	2M
        """
        
        Tempfile.mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "test_input.gfa")
            write(gfa_file, gfa_content)
            
            # Test reading in DoubleStrand mode (default)
            graph = Mycelia.read_gfa_next(gfa_file)
            
            Test.@test graph isa MetaGraphsNext.MetaGraph
            
            # Check vertices
            labels = collect(MetaGraphsNext.labels(graph))
            Test.@test "1" in labels
            Test.@test "2" in labels
            Test.@test length(labels) == 2
            
            # Check vertex data
            vertex1 = graph["1"]
            vertex2 = graph["2"]
            Test.@test vertex1 isa Mycelia.KmerVertexData
            Test.@test vertex2 isa Mycelia.KmerVertexData
            Test.@test vertex1.canonical_kmer == "ATC"
            Test.@test vertex2.canonical_kmer == "TCG"
            
            # Check edges
            Test.@test MetaGraphsNext.has_edge(graph, "1", "2")
            edge_data = graph["1", "2"]
            Test.@test edge_data isa Mycelia.KmerEdgeData
            Test.@test edge_data.src_strand == Mycelia.Forward
            Test.@test edge_data.dst_strand == Mycelia.Forward
        end
    end
    
    Test.@testset "GFA Round-trip" begin
        # Create original graph
        seq1 = FASTX.FASTA.Record("test1", "ATCG")
        observations = [seq1]
        kmer_type = BioSequences.DNAKmer{3}
        
        original_graph = Mycelia.build_kmer_graph_next(kmer_type, observations)
        
        Tempfile.mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "roundtrip.gfa")
            
            # Write graph to GFA
            Mycelia.write_gfa_next(original_graph, gfa_file)
            
            # Read it back
            restored_graph = Mycelia.read_gfa_next(gfa_file)
            
            # Compare basic properties
            original_labels = Set(MetaGraphsNext.labels(original_graph))
            restored_labels = Set(MetaGraphsNext.labels(restored_graph))
            
            # Note: Labels might be different (vertex IDs vs k-mer strings)
            # but the structure should be preserved
            Test.@test length(original_labels) == length(restored_labels)
            
            # Check that we have the same number of edges
            original_edges = length(collect(MetaGraphsNext.edge_labels(original_graph)))
            restored_edges = length(collect(MetaGraphsNext.edge_labels(restored_graph)))
            Test.@test original_edges == restored_edges
        end
    end
    
    Test.@testset "Graph Modes in GFA I/O" begin
        gfa_content = """
        H	VN:Z:1.0
        S	1	ATC
        S	2	GAT
        L	1	+	2	-	2M
        """
        
        Tempfile.mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "modes_test.gfa")
            write(gfa_file, gfa_content)
            
            # Test SingleStrand mode
            single_graph = Mycelia.read_gfa_next(gfa_file, Mycelia.SingleStrand)
            Test.@test single_graph isa MetaGraphsNext.MetaGraph
            
            vertex1_single = single_graph["1"]
            vertex2_single = single_graph["2"]
            Test.@test vertex1_single.canonical_kmer == "ATC"  # unchanged
            Test.@test vertex2_single.canonical_kmer == "GAT"  # unchanged
            
            # Test DoubleStrand mode
            double_graph = Mycelia.read_gfa_next(gfa_file, Mycelia.DoubleStrand)
            Test.@test double_graph isa MetaGraphsNext.MetaGraph
            
            vertex1_double = double_graph["1"]
            vertex2_double = double_graph["2"]
            Test.@test vertex1_double.canonical_kmer == "ATC"  # canonical
            Test.@test vertex2_double.canonical_kmer == "ATC"  # GAT -> ATC (canonical)
            
            # Check edge orientations
            if MetaGraphsNext.has_edge(single_graph, "1", "2")
                edge_single = single_graph["1", "2"]
                Test.@test edge_single.src_strand == Mycelia.Forward
                Test.@test edge_single.dst_strand == Mycelia.Reverse
            end
        end
    end
    
    Test.@testset "Error Handling" begin
        # Test with non-existent file
        Test.@test_throws SystemError Mycelia.read_gfa_next("nonexistent.gfa")
        
        # Test with empty graph
        empty_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData
        )
        
        Tempfile.mktempdir() do tmpdir
            empty_gfa = joinpath(tmpdir, "empty.gfa")
            result = Mycelia.write_gfa_next(empty_graph, empty_gfa)
            Test.@test result == empty_gfa
            Test.@test isfile(empty_gfa)
            
            # Should contain only header
            content = read(empty_gfa, String)
            lines = filter(!isempty, split(strip(content), '\n'))
            Test.@test length(lines) == 1  # Only header
            Test.@test startswith(lines[1], "H\t")
        end
    end
end
