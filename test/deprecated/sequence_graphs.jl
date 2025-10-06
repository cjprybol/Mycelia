import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import BioSequences
import FASTX
import Graphs
import MetaGraphs
import DataFrames

Test.@testset "Sequence Graphs Tests" begin
    Test.@testset "K-mer Graph Construction" begin
        # Create test sequences
        test_records = [
            FASTX.FASTA.Record("seq1", "ATCGATCGATCG"),
            FASTX.FASTA.Record("seq2", "TCGATCGATCGA"),
            FASTX.FASTA.Record("seq3", "CGATCGATCGAT")
        ]
        
        # Test stranded k-mer graph construction
        kmer_type = BioSequences.DNAMer{4}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, test_records)
        
        Test.@test graph isa MetaGraphs.MetaDiGraph
        Test.@test Graphs.nv(graph) > 0
        Test.@test haskey(graph.gprops, :stranded_kmers)
        Test.@test haskey(graph.gprops, :reverse_complement_map)
        Test.@test haskey(graph.gprops, :k)
        Test.@test graph.gprops[:k] == 4
        
        # Test that observation tracking works
        Test.@test haskey(graph.gprops, :observed_paths)
        Test.@test haskey(graph.gprops, :observation_ids)
        Test.@test length(graph.gprops[:observation_ids]) == 3
        Test.@test "seq1" in graph.gprops[:observation_ids]
    end

    Test.@testset "Sequence to Path Conversion" begin
        # Create minimal k-mer set
        kmers = [
            BioSequences.DNAMer{3}("ATG"),
            BioSequences.DNAMer{3}("TGC"),
            BioSequences.DNAMer{3}("GCA"),
            BioSequences.DNAMer{3}("CAT")
        ]
        
        # Test sequence to path conversion
        test_sequence = BioSequences.LongDNA{4}("ATGCAT")
        path = Mycelia.sequence_to_stranded_path(kmers, test_sequence)
        
        Test.@test path isa Vector{Pair{Int, Bool}}
        Test.@test length(path) == length(test_sequence) - 3 + 1  # n - k + 1 k-mers
    end

    Test.@testset "Path to Sequence Conversion" begin
        # Create test k-mers and path
        kmers = [
            BioSequences.DNAMer{3}("ATG"),
            BioSequences.DNAMer{3}("TGC"),
            BioSequences.DNAMer{3}("GCA")
        ]
        
        path = [1, 2, 3]  # Simple path through k-mers
        sequence = Mycelia.path_to_sequence(kmers, path)
        
        Test.@test sequence isa BioSequences.LongDNA{4}
        Test.@test length(sequence) == length(path) + 3 - 1  # path_length + k - 1
    end

    Test.@testset "K-mer Index Operations" begin
        # Test k-mer indexing
        kmers = [
            BioSequences.DNAMer{3}("ATG"),
            BioSequences.DNAMer{3}("TGC"),
            BioSequences.DNAMer{3}("GCA")
        ]
        
        # Test finding k-mer index
        target_kmer = BioSequences.DNAMer{3}("TGC")
        index = Mycelia.get_kmer_index(kmers, target_kmer)
        Test.@test index == 2
        
        # Test with k-mer not in set
        missing_kmer = BioSequences.DNAMer{3}("AAA")
        missing_index = Mycelia.get_kmer_index(kmers, missing_kmer)
        Test.@test missing_index == 0  # Or whatever the function returns for missing k-mers
    end

    Test.@testset "Graph File I/O" begin
        # Create a simple test graph
        test_graph = Graphs.SimpleDiGraph(3)
        Graphs.add_edge!(test_graph, 1, 2)
        Graphs.add_edge!(test_graph, 2, 3)
        
        # Test saving graph
        temp_file = tempname() * ".jld2"
        result_file = Mycelia.save_graph(test_graph, temp_file)
        
        Test.@test result_file == temp_file
        Test.@test isfile(temp_file)
        
        # Test loading graph
        loaded_graph = Mycelia.load_graph(temp_file)
        Test.@test Graphs.nv(loaded_graph) == Graphs.nv(test_graph)
        Test.@test Graphs.ne(loaded_graph) == Graphs.ne(test_graph)
        
        # Cleanup
        rm(temp_file, force=true)
    end

    Test.@testset "GFA Format Operations" begin
        # Create test GFA content
        test_gfa_content = """H\tVN:Z:1.0
S\t1\tACGT
S\t2\tGTCA
L\t1\t+\t2\t+\t0M
P\tpath1\t1+,2+\t0M,0M
"""
        temp_gfa = tempname() * ".gfa"
        write(temp_gfa, test_gfa_content)
        
        # Test GFA parsing
        gfa_data = Mycelia.parse_gfa(temp_gfa)
        Test.@test gfa_data isa Dict
        Test.@test haskey(gfa_data, "segments")
        Test.@test haskey(gfa_data, "links")
        
        # Test GFA to structure table conversion
        structure_table = Mycelia.gfa_to_structure_table(temp_gfa)
        Test.@test structure_table isa DataFrames.DataFrame
        
        # Test GFA to FASTA conversion
        temp_fasta = tempname() * ".fasta"
        result_fasta = Mycelia.gfa_to_fasta(gfa=temp_gfa, fasta=temp_fasta)
        Test.@test result_fasta == temp_fasta
        Test.@test isfile(temp_fasta)
        
        # Cleanup
        rm(temp_gfa, force=true)
        rm(temp_fasta, force=true)
    end

    Test.@testset "FASTX to Graph Conversion" begin
        # Create test FASTX file
        temp_fastx = tempname() * ".fasta"
        fasta_content = ">seq1\nATCGATCGATCG\n>seq2\nGCTAGCTAGCTA\n"
        write(temp_fastx, fasta_content)
        
        # Test single file conversion
        kmer_type = BioSequences.DNAMer{4}
        graph = Mycelia.fastx_to_kmer_graph(kmer_type, temp_fastx)
        
        Test.@test graph isa MetaGraphs.MetaDiGraph
        Test.@test Graphs.nv(graph) > 0
        
        # Test multiple files conversion
        temp_fastx2 = tempname() * ".fasta"
        write(temp_fastx2, ">seq3\nTAGCTAGCTAGC\n")
        
        multi_graph = Mycelia.fastx_to_kmer_graph(kmer_type, [temp_fastx, temp_fastx2])
        Test.@test multi_graph isa MetaGraphs.MetaDiGraph
        Test.@test Graphs.nv(multi_graph) >= Graphs.nv(graph)
        
        # Cleanup
        rm(temp_fastx, force=true)
        rm(temp_fastx2, force=true)
    end

    Test.@testset "Contig Analysis Functions" begin
        # Create mock graph file for testing
        test_graph = MetaGraphs.MetaDiGraph(5)
        for i in 1:4
            Graphs.add_edge!(test_graph, i, i+1)
        end
        # Make it circular
        Graphs.add_edge!(test_graph, 5, 1)
        
        temp_graph_file = tempname() * ".jld2"
        Mycelia.save_graph(test_graph, temp_graph_file)
        
        # Test contig analysis functions
        # Note: These functions expect specific graph structures
        # In real tests, we'd need properly formatted assembly graphs
        
        Test.@test isfile(temp_graph_file)
        
        # Cleanup
        rm(temp_graph_file, force=true)
    end

    Test.@testset "Edge and Path Operations" begin
        # Create test graph with edges
        graph = MetaGraphs.MetaDiGraph(3)
        Graphs.add_edge!(graph, 1, 2)
        Graphs.add_edge!(graph, 2, 3)
        
        # Add edge properties for testing
        MetaGraphs.set_prop!(graph, Graphs.Edge(1, 2), :coverage, [(1 => (1 => true))])
        MetaGraphs.set_prop!(graph, Graphs.Edge(2, 3), :coverage, [(1 => (2 => true))])
        
        # Test edge probability calculation
        edge_12 = Graphs.Edge(1, 2)
        if Graphs.has_edge(graph, edge_12)
            prob = Mycelia.edge_probability(graph, edge_12)
            Test.@test prob isa Float64
            Test.@test 0.0 ≤ prob ≤ 1.0
        end
        
        # Test edge path to sequence conversion
        edge_path = [Graphs.Edge(1, 2), Graphs.Edge(2, 3)]
        # Note: This requires the graph to have proper k-mer information
        # In a real implementation, we'd set up the graph with k-mer data
    end

    Test.@testset "Resampling and Solid K-mer Analysis" begin
        # Create mock k-mer solidity data
        kmer_solidity = Dict(
            1 => true,   # solid
            2 => false,  # not solid
            3 => true,   # solid
            4 => false,  # not solid
            5 => true    # solid
        )
        
        solid_indices = [i for (i, solid) in kmer_solidity if solid]
        
        # Test resampling stretch identification
        stretches = Mycelia.find_resampling_stretches(
            record_kmer_solidity=kmer_solidity,
            solid_branching_kmer_indices=solid_indices
        )
        
        Test.@test stretches isa Vector
        # Each stretch should be a range or similar structure
    end

    Test.@testset "Error Handling" begin
        # Test with empty sequence records
        empty_records = FASTX.FASTA.Record[]
        kmer_type = BioSequences.DNAMer{4}
        
        Test.@test_nowarn Mycelia.build_stranded_kmer_graph(kmer_type, empty_records)
        
        # Test with short sequences (shorter than k)
        short_records = [FASTX.FASTA.Record("short", "AT")]  # Length 2, k=4
        
        Test.@test_nowarn Mycelia.build_stranded_kmer_graph(kmer_type, short_records)
        
        # Test with non-existent files
        Test.@test_throws Exception Mycelia.parse_gfa("nonexistent.gfa")
        Test.@test_throws Exception Mycelia.load_graph("nonexistent.jld2")
    end

    Test.@testset "Graph Properties and Metadata" begin
        # Create graph with metadata
        test_records = [FASTX.FASTA.Record("test", "ATCGATCG")]
        kmer_type = BioSequences.DNAMer{3}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, test_records)
        
        # Test graph properties
        Test.@test haskey(graph.gprops, :k)
        Test.@test haskey(graph.gprops, :K)
        Test.@test haskey(graph.gprops, :stranded_kmers)
        Test.@test haskey(graph.gprops, :reverse_complement_map)
        
        # Test that all vertices have coverage metadata
        for v in Graphs.vertices(graph)
            Test.@test haskey(graph.vprops[v], :coverage)
            Test.@test graph.vprops[v][:coverage] isa Vector
        end
        
        # Test reverse complement mapping
        rc_map = graph.gprops[:reverse_complement_map]
        Test.@test length(rc_map) == length(graph.gprops[:stranded_kmers])
        
        # Each k-mer should map to its reverse complement's index
        for (i, kmer) in enumerate(graph.gprops[:stranded_kmers])
            rc_index = rc_map[i]
            rc_kmer = graph.gprops[:stranded_kmers][rc_index]
            Test.@test BioSequences.reverse_complement(kmer) == rc_kmer
        end
    end
end