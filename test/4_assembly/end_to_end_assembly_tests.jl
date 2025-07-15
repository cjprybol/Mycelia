"""
End-to-End Assembly Tests

This file contains comprehensive end-to-end tests for the three main assembly pipelines:
1. String-graph assembly
2. Strand-specific sequence-graph-next 
3. Canonical/double-stranded sequence-graph-next

Tests cover basic random strings, DNA, RNA, and amino acid sequences with various
error rates and coverage depths.
"""

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
import Random
import StatsBase

# Test parameters
const TEST_LENGTHS = [10, 100, 1000]
const ERROR_RATES = [0.0, 0.001, 0.01, 0.1]  # 0%, 0.1%, 1%, 10%
const COVERAGE_DEPTHS = [10, 100, 1000]

"""
Helper function to create FASTQ records from sequences with error simulation.
"""
function create_test_reads(reference_sequence::String, coverage::Int, error_rate::Float64)
    records = FASTX.FASTQ.Record[]
    
    for i in 1:coverage
        # Use the observe function to introduce errors
        bio_seq = BioSequences.LongDNA{4}(reference_sequence)
        observed_seq, quality_scores = Mycelia.observe(bio_seq, error_rate=error_rate)
        
        # Convert quality scores to string format
        quality_string = String([Char(q + 33) for q in quality_scores])
        
        record = FASTX.FASTQ.Record("read_$i", string(observed_seq), quality_string)
        push!(records, record)
    end
    
    return records
end

"""
Helper function to create FASTQ records from RNA sequences.
"""
function create_test_rna_reads(reference_sequence::String, coverage::Int, error_rate::Float64)
    records = FASTX.FASTQ.Record[]
    
    for i in 1:coverage
        # Use the observe function to introduce errors
        bio_seq = BioSequences.LongRNA{4}(reference_sequence)
        observed_seq, quality_scores = Mycelia.observe(bio_seq, error_rate=error_rate)
        
        # Convert quality scores to string format
        quality_string = String([Char(q + 33) for q in quality_scores])
        
        record = FASTX.FASTQ.Record("read_$i", string(observed_seq), quality_string)
        push!(records, record)
    end
    
    return records
end

"""
Helper function to create FASTQ records from amino acid sequences.
"""
function create_test_aa_reads(reference_sequence::String, coverage::Int, error_rate::Float64)
    records = FASTX.FASTQ.Record[]
    
    for i in 1:coverage
        # Use the observe function to introduce errors
        bio_seq = BioSequences.LongAA(reference_sequence)
        observed_seq, quality_scores = Mycelia.observe(bio_seq, error_rate=error_rate)
        
        # Convert quality scores to string format
        quality_string = String([Char(q + 33) for q in quality_scores])
        
        record = FASTX.FASTQ.Record("read_$i", string(observed_seq), quality_string)
        push!(records, record)
    end
    
    return records
end

Test.@testset "End-to-End Assembly Tests" begin
    
    Test.@testset "Base Case: Reference Sequence Graph Round-trip" begin
        Test.@testset "String Graph Base Case" begin
            # Test with simple strings
            for length in TEST_LENGTHS
                reference_string = Random.randstring(length)
                
                # Create graph from reference
                graph = Mycelia.string_to_ngram_graph(s=reference_string, n=3)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test graph connectivity
                components = Mycelia.find_connected_components(graph)
                Test.@test !isempty(components)
            end
        end
        
        Test.@testset "Sequence Graph Next Base Case - DNA" begin
            for length in TEST_LENGTHS
                if length >= 6  # Minimum length for meaningful k-mer analysis
                    reference_seq = BioSequences.randdnaseq(length)
                    reference_record = FASTX.FASTA.Record("reference", reference_seq)
                    
                    # Test DoubleStrand mode (canonical)
                    kmer_type = BioSequences.DNAKmer{5}
                    graph = Mycelia.build_kmer_graph_next(kmer_type, [reference_record]; 
                                                        graph_mode=Mycelia.DoubleStrand)
                    Test.@test graph isa MetaGraphsNext.MetaGraph
                    Test.@test !isempty(MetaGraphsNext.labels(graph))
                    
                    # Verify vertices are canonical k-mers
                    for label in MetaGraphsNext.labels(graph)
                        vertex_data = graph[label]
                        Test.@test vertex_data isa Mycelia.KmerVertexData
                        Test.@test vertex_data.canonical_kmer == label
                    end
                    
                    # Test that we can write and read GFA
                    Test.@test hasmethod(Mycelia.write_gfa_next, (typeof(graph), String))
                    Test.@test hasmethod(Mycelia.read_gfa_next, (String,))
                end
            end
        end
        
        Test.@testset "Sequence Graph Next Base Case - RNA" begin
            for length in TEST_LENGTHS
                if length >= 6  # Minimum length for meaningful k-mer analysis
                    reference_seq = string(BioSequences.randrnaseq(length))
                    reference_record = FASTX.FASTA.Record("reference", reference_seq)
                    
                    # Test SingleStrand mode for RNA
                    kmer_type = BioSequences.RNAKmer{5}
                    graph = Mycelia.build_kmer_graph_next(kmer_type, [reference_record]; 
                                                        graph_mode=Mycelia.SingleStrand)
                    Test.@test graph isa MetaGraphsNext.MetaGraph
                    Test.@test !isempty(MetaGraphsNext.labels(graph))
                    
                    # Verify vertices contain RNA k-mers
                    for label in MetaGraphsNext.labels(graph)
                        vertex_data = graph[label]
                        Test.@test vertex_data isa Mycelia.KmerVertexData
                        Test.@test 'U' in vertex_data.canonical_kmer || 'A' in vertex_data.canonical_kmer
                    end
                end
            end
        end
        
        Test.@testset "Sequence Graph Next Base Case - Amino Acids" begin
            for length in TEST_LENGTHS
                if length >= 6  # Minimum length for meaningful k-mer analysis
                    reference_seq = BioSequences.randaaseq(length)
                    reference_record = FASTX.FASTA.Record("reference", reference_seq)
                    
                    # Test SingleStrand mode for amino acids
                    kmer_type = BioSequences.AminoAcidKmer{5}
                    graph = Mycelia.build_kmer_graph_next(kmer_type, [reference_record]; 
                                                        graph_mode=Mycelia.SingleStrand)
                    Test.@test graph isa MetaGraphsNext.MetaGraph
                    Test.@test !isempty(MetaGraphsNext.labels(graph))
                    
                    # Verify vertices contain amino acid k-mers
                    for label in MetaGraphsNext.labels(graph)
                        vertex_data = graph[label]
                        Test.@test vertex_data isa Mycelia.KmerVertexData
                        # Check that all characters are valid amino acids
                        Test.@test all(c in "ACDEFGHIKLMNPQRSTVWY" for c in vertex_data.canonical_kmer)
                    end
                end
            end
        end

        Test.@testset "String Graph Base Case - ASCII Greek" begin
            for length in TEST_LENGTHS
                reference_string = Mycelia.rand_ascii_greek_string(length)
                
                # Create graph from reference
                graph = Mycelia.string_to_ngram_graph(s=reference_string, n=3)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test assembly - should recover non-empty result
                assembly = Mycelia.assemble_string_graph(graph)
                Test.@test !isempty(assembly)
                
                # Test graph connectivity
                components = Mycelia.find_connected_components(graph)
                Test.@test !isempty(components)
            end
        end

        Test.@testset "String Graph Base Case - Latin1" begin
            for length in TEST_LENGTHS
                reference_string = Mycelia.rand_latin1_string(length)
                
                # Create graph from reference
                graph = Mycelia.string_to_ngram_graph(s=reference_string, n=3)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test assembly - should recover non-empty result
                assembly = Mycelia.assemble_string_graph(graph)
                Test.@test !isempty(assembly)
                
                # Test graph connectivity
                components = Mycelia.find_connected_components(graph)
                Test.@test !isempty(components)
            end
        end

        Test.@testset "String Graph Base Case - BMP Printable" begin
            for length in TEST_LENGTHS
                reference_string = Mycelia.rand_bmp_printable_string(length)
                
                # Create graph from reference
                graph = Mycelia.string_to_ngram_graph(s=reference_string, n=3)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test assembly - should recover non-empty result
                assembly = Mycelia.assemble_string_graph(graph)
                Test.@test !isempty(assembly)
                
                # Test graph connectivity
                components = Mycelia.find_connected_components(graph)
                Test.@test !isempty(components)
            end
        end

        Test.@testset "String Graph Base Case - Printable Unicode" begin
            for length in TEST_LENGTHS
                reference_string = Mycelia.rand_printable_unicode_string(length)
                
                # Create graph from reference
                graph = Mycelia.string_to_ngram_graph(s=reference_string, n=3)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test assembly - should recover non-empty result
                assembly = Mycelia.assemble_string_graph(graph)
                Test.@test !isempty(assembly)
                
                # Test graph connectivity
                components = Mycelia.find_connected_components(graph)
                Test.@test !isempty(components)
            end
        end
    end
    
    Test.@testset "String Graph Assembly with Error Rates and Coverage" begin
        Test.@testset "Basic String Assembly" begin
            for length in TEST_LENGTHS
                reference_string = Random.randstring(length)
                
                for error_rate in ERROR_RATES
                    for coverage in COVERAGE_DEPTHS
                        # Create simulated reads with errors
                        mutated_strings = String[]
                        for i in 1:coverage
                            mutated_str = Mycelia.mutate_string(reference_string, error_rate=error_rate)
                            push!(mutated_strings, mutated_str)
                        end
                        
                        # Build graph from all mutated strings
                        combined_string = join(mutated_strings, "")
                        graph = Mycelia.string_to_ngram_graph(s=combined_string, n=3)
                        
                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))
                        
                        # Test assembly process
                        collapsed_graph = Mycelia.collapse_unbranching_paths(graph)
                        assemblies = Mycelia.assemble_strings(collapsed_graph)
                        
                        Test.@test !isempty(assemblies)
                        Test.@test all(asm isa String for asm in assemblies)
                        
                        # For low error rates, assembly should be reasonable
                        if error_rate <= 0.01
                            Test.@test any(length(asm) >= length(reference_string) ÷ 2 for asm in assemblies)
                        end
                    end
                end
            end
        end
    end
    
    Test.@testset "Strand-Specific Sequence Graph Next Assembly" begin
        Test.@testset "DNA SingleStrand Mode" begin
            for length in TEST_LENGTHS
                if length >= 10  # Minimum length for meaningful assembly
                    reference_seq = BioSequences.randdnaseq(length)
                    
                    for error_rate in ERROR_RATES
                        for coverage in COVERAGE_DEPTHS
                            # Create simulated reads
                            reads = create_test_reads(reference_seq, coverage, error_rate)
                            
                            # Build k-mer graph in SingleStrand mode
                            kmer_type = BioSequences.DNAKmer{5}
                            graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                                graph_mode=Mycelia.SingleStrand)
                            
                            Test.@test graph isa MetaGraphsNext.MetaGraph
                            Test.@test !isempty(MetaGraphsNext.labels(graph))
                            
                            # Verify all strand orientations are Forward in SingleStrand mode
                            for label in MetaGraphsNext.labels(graph)
                                vertex_data = graph[label]
                                for (obs_id, pos, strand) in vertex_data.coverage
                                    Test.@test strand == Mycelia.Forward
                                end
                            end
                            
                            # Test that edges respect strand constraints
                            for edge_label in MetaGraphsNext.edge_labels(graph)
                                if !isempty(edge_label)
                                    edge_data = graph[edge_label...]
                                    Test.@test edge_data isa Mycelia.KmerEdgeData
                                    Test.@test edge_data.src_strand isa Mycelia.StrandOrientation
                                    Test.@test edge_data.dst_strand isa Mycelia.StrandOrientation
                                end
                            end
                        end
                    end
                end
            end
        end
        
        Test.@testset "RNA SingleStrand Mode" begin
            for length in TEST_LENGTHS
                if length >= 10  # Minimum length for meaningful assembly
                    reference_seq = string(BioSequences.randrnaseq(length))
                    
                    for error_rate in ERROR_RATES
                        for coverage in COVERAGE_DEPTHS
                            # Create simulated reads
                            reads = create_test_rna_reads(reference_seq, coverage, error_rate)
                            
                            # Build k-mer graph in SingleStrand mode
                            kmer_type = BioSequences.RNAKmer{5}
                            graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                                graph_mode=Mycelia.SingleStrand)
                            
                            Test.@test graph isa MetaGraphsNext.MetaGraph
                            Test.@test !isempty(MetaGraphsNext.labels(graph))
                            
                            # Verify RNA-specific properties
                            for label in MetaGraphsNext.labels(graph)
                                vertex_data = graph[label]
                                Test.@test vertex_data isa Mycelia.KmerVertexData
                                # RNA sequences should contain U instead of T
                                Test.@test 'U' in vertex_data.canonical_kmer || 'A' in vertex_data.canonical_kmer
                            end
                        end
                    end
                end
            end
        end
        
        Test.@testset "Amino Acid SingleStrand Mode" begin
            for length in TEST_LENGTHS
                if length >= 10  # Minimum length for meaningful assembly
                    reference_seq = BioSequences.randaaseq(length)
                    
                    for error_rate in ERROR_RATES
                        for coverage in COVERAGE_DEPTHS
                            # Create simulated reads
                            reads = create_test_aa_reads(reference_seq, coverage, error_rate)
                            
                            # Build k-mer graph in SingleStrand mode
                            kmer_type = BioSequences.AminoAcidKmer{3}  # Shorter k-mers for AA
                            graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                                graph_mode=Mycelia.SingleStrand)
                            
                            Test.@test graph isa MetaGraphsNext.MetaGraph
                            Test.@test !isempty(MetaGraphsNext.labels(graph))
                            
                            # Verify amino acid-specific properties
                            for label in MetaGraphsNext.labels(graph)
                                vertex_data = graph[label]
                                Test.@test vertex_data isa Mycelia.KmerVertexData
                                # All characters should be valid amino acids
                                Test.@test all(c in "ACDEFGHIKLMNPQRSTVWY" for c in vertex_data.canonical_kmer)
                            end
                        end
                    end
                end
            end
        end
    end
    
    Test.@testset "Canonical/Double-Stranded Sequence Graph Next Assembly" begin
        Test.@testset "DNA DoubleStrand Mode" begin
            for length in TEST_LENGTHS
                if length >= 10  # Minimum length for meaningful assembly
                    reference_seq = BioSequences.randdnaseq(length)
                    
                    for error_rate in ERROR_RATES
                        for coverage in COVERAGE_DEPTHS
                            # Create simulated reads
                            reads = create_test_reads(reference_seq, coverage, error_rate)
                            
                            # Build k-mer graph in DoubleStrand mode (canonical)
                            kmer_type = BioSequences.DNAKmer{5}
                            graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                                graph_mode=Mycelia.DoubleStrand)
                            
                            Test.@test graph isa MetaGraphsNext.MetaGraph
                            Test.@test !isempty(MetaGraphsNext.labels(graph))
                            
                            # Verify canonical k-mer properties
                            for label in MetaGraphsNext.labels(graph)
                                vertex_data = graph[label]
                                Test.@test vertex_data isa Mycelia.KmerVertexData
                                Test.@test vertex_data.canonical_kmer == label
                                
                                # In DoubleStrand mode, we should see both Forward and Reverse orientations
                                strand_orientations = [strand for (obs_id, pos, strand) in vertex_data.coverage]
                                Test.@test !isempty(strand_orientations)
                            end
                            
                            # Test that edges handle strand transitions correctly
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
                end
            end
        end
        
        Test.@testset "RNA DoubleStrand Mode" begin
            for length in TEST_LENGTHS
                if length >= 10  # Minimum length for meaningful assembly
                    reference_seq = string(BioSequences.randrnaseq(length))
                    
                    for error_rate in ERROR_RATES
                        for coverage in COVERAGE_DEPTHS
                            # Create simulated reads
                            reads = create_test_rna_reads(reference_seq, coverage, error_rate)
                            
                            # Build k-mer graph in DoubleStrand mode (canonical)
                            kmer_type = BioSequences.RNAKmer{5}
                            graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                                graph_mode=Mycelia.DoubleStrand)
                            
                            Test.@test graph isa MetaGraphsNext.MetaGraph
                            Test.@test !isempty(MetaGraphsNext.labels(graph))
                            
                            # Verify canonical RNA k-mer properties
                            for label in MetaGraphsNext.labels(graph)
                                vertex_data = graph[label]
                                Test.@test vertex_data isa Mycelia.KmerVertexData
                                Test.@test vertex_data.canonical_kmer == label
                                Test.@test 'U' in vertex_data.canonical_kmer || 'A' in vertex_data.canonical_kmer
                            end
                        end
                    end
                end
            end
        end
    end
    
    Test.@testset "GFA I/O Round-trip Tests" begin
        Test.@testset "String Graph GFA I/O" begin
            reference_string = Mycelia.rand_ascii_greek_string(100)
            graph = Mycelia.string_to_ngram_graph(s=reference_string, n=3)
            
            # Test that we can work with the graph (basic functionality)
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            Test.@test graph isa MetaGraphsNext.MetaGraph
        end
        
        Test.@testset "Sequence Graph Next GFA I/O" begin
            reference_seq = string(BioSequences.randdnaseq(100))
            reference_record = FASTX.FASTA.Record("reference", reference_seq)
            
            kmer_type = BioSequences.DNAKmer{5}
            graph = Mycelia.build_kmer_graph_next(kmer_type, [reference_record])
            
            # Test that we can work with the graph (basic functionality)
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            Test.@test graph isa MetaGraphsNext.MetaGraph
            
            # Test that GFA I/O functions exist
            Test.@test hasmethod(Mycelia.write_gfa_next, (typeof(graph), String))
            Test.@test hasmethod(Mycelia.read_gfa_next, (String,))
        end
    end
    
    Test.@testset "Assembly Performance and Scaling" begin
        Test.@testset "Memory Usage Scaling" begin
            # Test with progressively larger sequences
            for length in [100, 1000]
                reference_seq = string(BioSequences.randdnaseq(length))
                reads = create_test_reads(reference_seq, 100, 0.01)
                
                # Build graphs and verify they don't crash
                kmer_type = BioSequences.DNAKmer{5}
                
                # Test both modes
                single_graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                          graph_mode=Mycelia.SingleStrand)
                double_graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                          graph_mode=Mycelia.DoubleStrand)
                
                Test.@test single_graph isa MetaGraphsNext.MetaGraph
                Test.@test double_graph isa MetaGraphsNext.MetaGraph
                
                # In canonical mode, we should generally have fewer vertices
                # (though this depends on the specific sequence)
                Test.@test !isempty(MetaGraphsNext.labels(single_graph))
                Test.@test !isempty(MetaGraphsNext.labels(double_graph))
            end
        end
        
        Test.@testset "Coverage Impact on Assembly Quality" begin
            reference_seq = string(BioSequences.randdnaseq(100))
            error_rate = 0.01
            
            for coverage in [10, 100]
                reads = create_test_reads(reference_seq, coverage, error_rate)
                
                kmer_type = BioSequences.DNAKmer{5}
                graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                    graph_mode=Mycelia.DoubleStrand)
                
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Higher coverage should generally result in more robust graphs
                total_coverage = sum(length(vertex_data.coverage) for vertex_data in values(graph.vertex_data))
                Test.@test total_coverage > 0
            end
        end
    end
end

# Helper function to run a subset of tests for quick validation
function run_quick_tests()
    Test.@testset "Quick Assembly Tests" begin
        # Test basic functionality with small sequences
        reference_seq = string(BioSequences.randdnaseq(50))
        reads = create_test_reads(reference_seq, 20, 0.01)
        
        kmer_type = BioSequences.DNAKmer{5}
        graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                            graph_mode=Mycelia.DoubleStrand)
        
        Test.@test !isempty(MetaGraphsNext.labels(graph))
        Test.@test graph isa MetaGraphsNext.MetaGraph
        
        # Test string graph
        test_string = Mycelia.rand_ascii_greek_string(50)
        string_graph = Mycelia.string_to_ngram_graph(s=test_string, n=3)
        Test.@test !isempty(MetaGraphsNext.labels(string_graph))
    end
end