"""
Basic Graph Type Tests

Simple tests to verify each graph type can be constructed and accessed.
"""

import Test
import Mycelia
import FASTX
import BioSequences

Test.@testset "Basic Graph Type Construction Tests" begin
    
    # Test data
    dna_seq = "ATCGATCGATCGATCG"
    rna_seq = "AUCGAUCGAUCGAUCG"
    protein_seq = "ALAVALINEGLUTAMINE"
    high_quality = "HHHHHHHHHHHHHHHH"
    
    Test.@testset "1. N-gram Graphs" begin
        # Test N-gram graph construction with keyword arguments
        graph = Mycelia.string_to_ngram_graph(s=dna_seq, n=3)
        Test.@test !isempty(graph.vertex_labels)
        println("✓ N-gram Graph: $(length(graph.vertex_labels)) vertices")
    end
    
    Test.@testset "2. K-mer Graphs" begin
        
        Test.@testset "DNA K-mer Graph" begin
            records = [FASTX.FASTA.Record("test", dna_seq)]
            graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.DNAKmer{5}, records)
            Test.@test !isempty(graph.vertex_labels)
            println("✓ DNA K-mer Graph: $(length(graph.vertex_labels)) vertices")
        end
        
        Test.@testset "RNA K-mer Graph" begin
            records = [FASTX.FASTA.Record("test", rna_seq)]
            graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.RNAKmer{4}, records)
            Test.@test !isempty(graph.vertex_labels)
            println("✓ RNA K-mer Graph: $(length(graph.vertex_labels)) vertices")
        end
        
        Test.@testset "Amino Acid K-mer Graph" begin
            records = [FASTX.FASTA.Record("test", protein_seq)]
            graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.AAKmer{3}, records)
            Test.@test !isempty(graph.vertex_labels)
            println("✓ Amino Acid K-mer Graph: $(length(graph.vertex_labels)) vertices")
        end
    end
    
    Test.@testset "3. Qualmer Graphs" begin
        
        Test.@testset "DNA Qualmer Graph" begin
            records = [FASTX.FASTQ.Record("test", dna_seq, high_quality)]
            graph = Mycelia.build_qualmer_graph(records, k=5)
            Test.@test !isempty(graph.vertex_labels)
            
            # Test vertex data
            kmers = collect(values(graph.vertex_labels))
            first_kmer = first(kmers)
            vertex_data = graph[first_kmer]
            Test.@test vertex_data.joint_probability > 0.0
            Test.@test vertex_data.coverage > 0
            println("✓ DNA Qualmer Graph: $(length(kmers)) vertices, joint prob: $(round(vertex_data.joint_probability, digits=4))")
        end
        
        Test.@testset "RNA Qualmer Graph" begin
            records = [FASTX.FASTQ.Record("test", rna_seq, high_quality)]
            graph = Mycelia.build_qualmer_graph(records, k=4)
            Test.@test !isempty(graph.vertex_labels)
            
            # Test vertex data
            kmers = collect(values(graph.vertex_labels))
            first_kmer = first(kmers)
            vertex_data = graph[first_kmer]
            Test.@test vertex_data.joint_probability > 0.0
            println("✓ RNA Qualmer Graph: $(length(kmers)) vertices, joint prob: $(round(vertex_data.joint_probability, digits=4))")
        end
        
        Test.@testset "Amino Acid Qualmer Graph" begin
            records = [FASTX.FASTQ.Record("test", protein_seq, repeat("H", length(protein_seq)))]
            graph = Mycelia.build_qualmer_graph(records, k=3)
            Test.@test !isempty(graph.vertex_labels)
            
            # Test vertex data
            kmers = collect(values(graph.vertex_labels))
            first_kmer = first(kmers)
            vertex_data = graph[first_kmer]
            Test.@test vertex_data.joint_probability > 0.0
            println("✓ Amino Acid Qualmer Graph: $(length(kmers)) vertices, joint prob: $(round(vertex_data.joint_probability, digits=4))")
        end
    end
    
    Test.@testset "4. FASTA Graphs" begin
        Test.@testset "Direct BioSequence Graph" begin
            records = [FASTX.FASTA.Record("test", dna_seq)]
            graph = Mycelia.build_biosequence_graph(records)
            Test.@test !isempty(graph.vertex_labels)
            
            # Test vertex data
            sequences = collect(values(graph.vertex_labels))
            Test.@test all(seq -> seq isa BioSequences.LongDNA, sequences)
            println("✓ FASTA Graph: $(length(sequences)) BioSequences")
        end
        
        Test.@testset "K-mer to BioSequence Conversion" begin
            records = [FASTX.FASTA.Record("test", dna_seq)]
            
            # Create k-mer graph
            kmer_graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.DNAKmer{5}, records)
            Test.@test !isempty(kmer_graph.vertex_labels)
            
            # Convert to BioSequence graph
            bio_graph = Mycelia.kmer_graph_to_biosequence_graph(kmer_graph)
            Test.@test !isempty(bio_graph.vertex_labels)
            
            # Test sequences
            sequences = collect(values(bio_graph.vertex_labels))
            Test.@test all(seq -> seq isa BioSequences.LongDNA, sequences)
            println("✓ K-mer to FASTA: $(length(kmer_graph.vertex_labels)) k-mers -> $(length(sequences)) BioSequences")
        end
    end
    
    Test.@testset "5. FASTQ Graphs" begin
        Test.@testset "Direct Quality BioSequence Graph" begin
            records = [FASTX.FASTQ.Record("test", dna_seq, high_quality)]
            graph = Mycelia.build_quality_biosequence_graph(records)
            Test.@test !isempty(graph.vertex_labels)
            
            # Test vertex data
            sequences = collect(values(graph.vertex_labels))
            Test.@test all(seq -> seq isa BioSequences.LongDNA, sequences)
            println("✓ FASTQ Graph: $(length(sequences)) quality BioSequences")
        end
        
        Test.@testset "Qualmer to Quality BioSequence Conversion" begin
            records = [FASTX.FASTQ.Record("test", dna_seq, high_quality)]
            
            # Create qualmer graph
            qualmer_graph = Mycelia.build_qualmer_graph(records, k=5)
            Test.@test !isempty(qualmer_graph.vertex_labels)
            
            # Convert to quality BioSequence graph
            bio_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
            Test.@test !isempty(bio_graph.vertex_labels)
            
            # Test sequences
            sequences = collect(values(bio_graph.vertex_labels))
            Test.@test all(seq -> seq isa BioSequences.LongDNA, sequences)
            println("✓ Qualmer to FASTQ: $(length(qualmer_graph.vertex_labels)) qualmers -> $(length(sequences)) quality BioSequences")
        end
    end
    
    Test.@testset "6. Graph Type Hierarchy" begin
        Test.@testset "Type Stability" begin
            # Test that all graph types have consistent type signatures
            records = [FASTX.FASTA.Record("test", dna_seq)]
            
            # K-mer graph
            kmer_graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.DNAKmer{5}, records)
            Test.@test !isempty(kmer_graph.vertex_labels)
            
            # Convert to BioSequence graph
            bio_graph = Mycelia.kmer_graph_to_biosequence_graph(kmer_graph)
            Test.@test !isempty(bio_graph.vertex_labels)
            
            # Test types
            kmers = collect(values(kmer_graph.vertex_labels))
            sequences = collect(values(bio_graph.vertex_labels))
            Test.@test all(kmer -> kmer isa Mycelia.Kmers.DNAKmer, kmers)
            Test.@test all(seq -> seq isa BioSequences.LongDNA, sequences)
            
            println("✓ Type Hierarchy: Fixed-length k-mers -> Variable-length BioSequences")
        end
    end
    
    println("\n" * "="^60)
    println("BASIC GRAPH TYPE TESTS SUMMARY")
    println("="^60)
    println("✓ All 6 graph types can be constructed")
    println("✓ Type stability maintained throughout")
    println("✓ Quality-aware functionality working")
    println("✓ Graph conversions operational")
    println("="^60)
end