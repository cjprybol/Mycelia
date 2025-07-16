# julia --project=test -e 'include("test/4_assembly/six_graph_hierarchy_tests.jl")'
# julia --project=. -e 'include("test/4_assembly/six_graph_hierarchy_tests.jl")'

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
import Kmers

Test.@testset "Complete 6-Graph Hierarchy Tests" begin
    
    # Test data creation
    function create_test_data()
        # Create test DNA sequences
        dna_records = [
            FASTX.FASTA.Record("dna1", "ATCGATCGATCGATCG"),
            FASTX.FASTA.Record("dna2", "TCGATCGATCGATCGA"),
            FASTX.FASTA.Record("dna3", "CGATCGATCGATCGAT")
        ]
        
        # Create test FASTQ records with quality scores
        fastq_records = [
            FASTX.FASTQ.Record("read1", "ATCGATCGATCGATCG", "IIIIIIIIIIIIIIII"),
            FASTX.FASTQ.Record("read2", "TCGATCGATCGATCGA", "HHHHHHHHHHHHHHHH"),
            FASTX.FASTQ.Record("read3", "CGATCGATCGATCGAT", "GGGGGGGGGGGGGGGG")
        ]
        
        # Create test RNA records
        rna_records = [
            FASTX.FASTA.Record("rna1", "AUCGAUCGAUCGAUCG"),
            FASTX.FASTA.Record("rna2", "UCGAUCGAUCGAUCGA")
        ]
        
        # Create test protein records
        protein_records = [
            FASTX.FASTA.Record("prot1", "MWKLVPGKEC"),
            FASTX.FASTA.Record("prot2", "KLVPGKECMW")
        ]
        
        return dna_records, fastq_records, rna_records, protein_records
    end
    
    Test.@testset "Fixed-Length Graph Types (Assembly Foundation)" begin
        
        Test.@testset "1. N-gram Graphs - Unicode Character Analysis" begin
            # Test basic N-gram graph construction
            test_string = "ABCDEFABCDEF"
            n = 3
            
            ngram_graph = Mycelia.string_to_ngram_graph(s=test_string, n=n)
            
            # Type checking
            Test.@test ngram_graph isa MetaGraphsNext.MetaGraph
            Test.@test all(label -> label isa String, MetaGraphsNext.labels(ngram_graph))
            Test.@test all(label -> length(label) == n, MetaGraphsNext.labels(ngram_graph))
            
            # Expected n-grams
            expected_ngrams = Set(["ABC", "BCD", "CDE", "DEF", "EFA"])
            actual_ngrams = Set(MetaGraphsNext.labels(ngram_graph))
            Test.@test issubset(expected_ngrams, actual_ngrams)
            
            # Test fixed-length vertices
            for label in MetaGraphsNext.labels(ngram_graph)
                Test.@test length(label) == n
                Test.@test label isa String
            end
        end
        
        Test.@testset "2. K-mer Graphs - BioSequence Analysis (NO strings)" begin
            dna_records, _, rna_records, protein_records = create_test_data()
            
            # Test DNA k-mer graphs
            Test.@testset "DNA K-mer Graphs" begin
                k = 5
                kmer_type = Kmers.DNAKmer{k}
                dna_graph = Mycelia.build_kmer_graph_next(kmer_type, dna_records)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test dna_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa Kmers.DNAKmer{k}, MetaGraphsNext.labels(dna_graph))
                
                # Verify fixed-length vertices
                for label in MetaGraphsNext.labels(dna_graph)
                    Test.@test length(label) == k
                    Test.@test label isa Kmers.DNAKmer{k}
                end
                
                # Test vertex data
                for label in MetaGraphsNext.labels(dna_graph)
                    vertex_data = dna_graph[label]
                    Test.@test vertex_data isa Mycelia.KmerVertexData
                    Test.@test vertex_data.canonical_kmer == label
                    Test.@test vertex_data.coverage isa Vector{Tuple{Int, Int, Mycelia.StrandOrientation}}
                end
            end
            
            # Test RNA k-mer graphs
            Test.@testset "RNA K-mer Graphs" begin
                k = 4
                kmer_type = Kmers.RNAKmer{k}
                rna_graph = Mycelia.build_kmer_graph_next(kmer_type, rna_records)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test rna_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa Kmers.RNAKmer{k}, MetaGraphsNext.labels(rna_graph))
                
                # Verify fixed-length vertices
                for label in MetaGraphsNext.labels(rna_graph)
                    Test.@test length(label) == k
                    Test.@test label isa Kmers.RNAKmer{k}
                end
            end
            
            # Test amino acid k-mer graphs
            Test.@testset "Amino Acid K-mer Graphs" begin
                k = 3
                kmer_type = Kmers.AAKmer{k}
                aa_graph = Mycelia.build_kmer_graph_next(kmer_type, protein_records)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test aa_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa Kmers.AAKmer{k}, MetaGraphsNext.labels(aa_graph))
                
                # Verify fixed-length vertices
                for label in MetaGraphsNext.labels(aa_graph)
                    Test.@test length(label) == k
                    Test.@test label isa Kmers.AAKmer{k}
                end
            end
        end
        
        Test.@testset "3. Qualmer Graphs - Quality-Aware Analysis (NO strings)" begin
            _, fastq_records, _, _ = create_test_data()
            
            # Test DNA qualmer graphs
            Test.@testset "DNA Qualmer Graphs" begin
                k = 5
                qualmer_graph = Mycelia.build_qualmer_graph(fastq_records, k=k)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph
                
                # Verify fixed-length vertices with quality information
                for label in MetaGraphsNext.labels(qualmer_graph)
                    vertex_data = qualmer_graph[label]
                    Test.@test vertex_data isa Mycelia.QualmerVertexData
                    Test.@test vertex_data.canonical_qualmer.kmer isa Kmers.DNAKmer{k}
                    Test.@test length(vertex_data.canonical_qualmer.kmer) == k
                    Test.@test length(vertex_data.canonical_qualmer.qualities) == k
                    Test.@test vertex_data.joint_probability isa Float64
                    Test.@test 0.0 <= vertex_data.joint_probability <= 1.0
                end
                
                # Test edge data with quality weights
                for (src, dst) in MetaGraphsNext.edge_labels(qualmer_graph)
                    edge_data = qualmer_graph[src, dst]
                    Test.@test edge_data isa Mycelia.QualmerEdgeData
                    Test.@test edge_data.quality_weight isa Float64
                    Test.@test edge_data.quality_weight > 0.0
                end
            end
        end
    end
    
    Test.@testset "Variable-Length Graph Types (Simplified Products)" begin
        
        Test.@testset "4. String Graphs - Simplified N-gram Graphs" begin
            # Create N-gram graph first
            test_string = "ABCDEFABCDEF"
            n = 3
            ngram_graph = Mycelia.string_to_ngram_graph(s=test_string, n=n)
            
            # Test string graph operations (simplified from N-gram)
            collapsed_graph = Mycelia.collapse_unbranching_paths(ngram_graph)
            
            # Type checking
            Test.@test collapsed_graph isa MetaGraphsNext.MetaGraph
            
            # Test assembly from string graph
            assembled_strings = Mycelia.assemble_strings(collapsed_graph)
            Test.@test assembled_strings isa Vector{String}
            Test.@test !isempty(assembled_strings)
            
            # Variable-length vertices
            for assembled in assembled_strings
                Test.@test assembled isa String
                Test.@test length(assembled) >= n  # Should be at least n characters
            end
        end
        
        Test.@testset "5. FASTA Graphs - Simplified K-mer Graphs (NO strings)" begin
            dna_records, _, rna_records, protein_records = create_test_data()
            
            # Test direct BioSequence graph construction
            Test.@testset "Direct BioSequence Graph from FASTA" begin
                biosequence_graph = Mycelia.build_biosequence_graph(dna_records)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test biosequence_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(biosequence_graph))
                
                # Variable-length vertices
                for label in MetaGraphsNext.labels(biosequence_graph)
                    Test.@test label isa BioSequences.LongDNA{4}
                    Test.@test length(label) > 0  # Variable length
                    vertex_data = biosequence_graph[label]
                    Test.@test vertex_data isa Mycelia.BioSequenceVertexData
                    Test.@test vertex_data.sequence == label
                end
            end
            
            # Test k-mer graph to BioSequence graph conversion
            Test.@testset "K-mer to BioSequence Graph Conversion" begin
                k = 5
                kmer_type = Kmers.DNAKmer{k}
                kmer_graph = Mycelia.build_kmer_graph_next(kmer_type, dna_records)
                
                # Convert to BioSequence graph
                biosequence_graph = Mycelia.kmer_graph_to_biosequence_graph(kmer_graph)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test biosequence_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(biosequence_graph))
                
                # Variable-length vertices from k-mer path simplification
                for label in MetaGraphsNext.labels(biosequence_graph)
                    Test.@test label isa BioSequences.LongDNA{4}
                    Test.@test length(label) >= k  # Should be at least k bases
                    vertex_data = biosequence_graph[label]
                    Test.@test vertex_data isa Mycelia.BioSequenceVertexData
                    Test.@test !isempty(vertex_data.constituent_kmers)
                end
            end
        end
        
        Test.@testset "6. FASTQ Graphs - Quality-Aware BioSequence Graphs (NO strings)" begin
            _, fastq_records, _, _ = create_test_data()
            
            # Test direct quality-aware BioSequence graph construction
            Test.@testset "Direct Quality BioSequence Graph from FASTQ" begin
                quality_graph = Mycelia.build_quality_biosequence_graph(fastq_records)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test quality_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(quality_graph))
                
                # Variable-length vertices with quality preservation
                for label in MetaGraphsNext.labels(quality_graph)
                    Test.@test label isa BioSequences.LongDNA{4}
                    Test.@test length(label) > 0  # Variable length
                    vertex_data = quality_graph[label]
                    Test.@test vertex_data isa Mycelia.QualityBioSequenceVertexData
                    Test.@test vertex_data.sequence == label
                    Test.@test length(vertex_data.quality_scores) == length(label)
                    Test.@test vertex_data.joint_probability isa Float64
                    Test.@test 0.0 <= vertex_data.joint_probability <= 1.0
                end
            end
            
            # Test Qualmer graph to quality-aware BioSequence graph conversion
            Test.@testset "Qualmer to Quality BioSequence Graph Conversion" begin
                k = 5
                qualmer_graph = Mycelia.build_qualmer_graph(fastq_records, k=k)
                
                # Convert to quality-aware BioSequence graph
                quality_biosequence_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test quality_biosequence_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(quality_biosequence_graph))
                
                # Variable-length vertices with quality preservation
                for label in MetaGraphsNext.labels(quality_biosequence_graph)
                    Test.@test label isa BioSequences.LongDNA{4}
                    Test.@test length(label) >= k  # Should be at least k bases
                    vertex_data = quality_biosequence_graph[label]
                    Test.@test vertex_data isa Mycelia.QualityBioSequenceVertexData
                    Test.@test vertex_data.sequence == label
                    Test.@test length(vertex_data.quality_scores) == length(label)
                    Test.@test !isempty(vertex_data.constituent_qualmers)
                end
            end
            
            # Test FASTQ conversion roundtrip
            Test.@testset "FASTQ Conversion Roundtrip" begin
                quality_graph = Mycelia.build_quality_biosequence_graph(fastq_records)
                
                # Convert back to FASTQ records
                reconstructed_fastq = Mycelia.quality_biosequence_graph_to_fastq(quality_graph)
                
                # Type checking
                Test.@test reconstructed_fastq isa Vector{FASTX.FASTQ.Record}
                Test.@test !isempty(reconstructed_fastq)
                
                # Quality preservation
                for record in reconstructed_fastq
                    Test.@test record isa FASTX.FASTQ.Record
                    Test.@test !isempty(FASTX.FASTQ.sequence(record))
                    Test.@test !isempty(FASTX.FASTQ.quality_scores(record))
                    Test.@test length(FASTX.FASTQ.sequence(record)) == length(FASTX.FASTQ.quality_scores(record))
                end
            end
        end
    end
    
    Test.@testset "Graph Type Hierarchy Integration" begin
        
        Test.@testset "Type Stability Across Hierarchy" begin
            dna_records, fastq_records, _, _ = create_test_data()
            
            # Test that each graph type maintains proper types
            k = 4
            
            # 1. N-gram graphs → strings
            ngram_graph = Mycelia.string_to_ngram_graph(s="ABCDEFABCDEF", n=3)
            Test.@test all(label -> label isa String, MetaGraphsNext.labels(ngram_graph))
            
            # 2. K-mer graphs → DNAKmer types
            kmer_graph = Mycelia.build_kmer_graph_next(Kmers.DNAKmer{k}, dna_records)
            Test.@test all(label -> label isa Kmers.DNAKmer{k}, MetaGraphsNext.labels(kmer_graph))
            
            # 3. Qualmer graphs → Qualmer types
            qualmer_graph = Mycelia.build_qualmer_graph(fastq_records, k=k)
            for label in MetaGraphsNext.labels(qualmer_graph)
                vertex_data = qualmer_graph[label]
                Test.@test vertex_data.canonical_qualmer.kmer isa Kmers.DNAKmer{k}
            end
            
            # 4. String graphs → strings (variable length)
            collapsed_graph = Mycelia.collapse_unbranching_paths(ngram_graph)
            assembled_strings = Mycelia.assemble_strings(collapsed_graph)
            Test.@test all(s -> s isa String, assembled_strings)
            
            # 5. BioSequence graphs → BioSequence types
            biosequence_graph = Mycelia.build_biosequence_graph(dna_records)
            Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(biosequence_graph))
            
            # 6. Quality BioSequence graphs → BioSequence types with quality
            quality_graph = Mycelia.build_quality_biosequence_graph(fastq_records)
            Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(quality_graph))
            for label in MetaGraphsNext.labels(quality_graph)
                vertex_data = quality_graph[label]
                Test.@test length(vertex_data.quality_scores) == length(label)
            end
        end
        
        Test.@testset "Strand-Aware Functionality" begin
            dna_records, fastq_records, _, _ = create_test_data()
            
            # Test SingleStrand mode
            Test.@testset "SingleStrand Mode" begin
                k = 4
                kmer_type = Kmers.DNAKmer{k}
                single_strand_graph = Mycelia.build_kmer_graph_next(kmer_type, dna_records, graph_mode=Mycelia.SingleStrand)
                
                # Check strand information in edges
                for (src, dst) in MetaGraphsNext.edge_labels(single_strand_graph)
                    edge_data = single_strand_graph[src, dst]
                    Test.@test edge_data isa Mycelia.KmerEdgeData
                    Test.@test edge_data.src_strand isa Mycelia.StrandOrientation
                    Test.@test edge_data.dst_strand isa Mycelia.StrandOrientation
                end
            end
            
            # Test DoubleStrand mode
            Test.@testset "DoubleStrand Mode" begin
                k = 4
                kmer_type = Kmers.DNAKmer{k}
                double_strand_graph = Mycelia.build_kmer_graph_next(kmer_type, dna_records, graph_mode=Mycelia.DoubleStrand)
                
                # Check canonical representation
                for label in MetaGraphsNext.labels(double_strand_graph)
                    vertex_data = double_strand_graph[label]
                    Test.@test vertex_data.canonical_kmer == label
                end
            end
        end
    end
    
    Test.@testset "Assembly Pipeline Integration" begin
        
        Test.@testset "Automatic Type Detection" begin
            dna_records, fastq_records, rna_records, protein_records = create_test_data()
            
            # Test DNA assembly (should default to Qualmer for FASTQ)
            Test.@testset "DNA Assembly Type Detection" begin
                # FASTQ input should use Qualmer graphs
                contigs_fastq = Mycelia.assemble_genome(fastq_records, method=Mycelia.QualmerGraph, k=5)
                Test.@test contigs_fastq isa Vector{String}
                Test.@test !isempty(contigs_fastq)
                
                # FASTA input should use K-mer graphs
                contigs_fasta = Mycelia.assemble_genome(dna_records, method=Mycelia.KmerGraph, k=5)
                Test.@test contigs_fasta isa Vector{String}
                Test.@test !isempty(contigs_fasta)
            end
            
            # Test RNA assembly
            Test.@testset "RNA Assembly Type Detection" begin
                contigs_rna = Mycelia.assemble_genome(rna_records, method=Mycelia.KmerGraph, k=4)
                Test.@test contigs_rna isa Vector{String}
                Test.@test !isempty(contigs_rna)
            end
            
            # Test protein assembly
            Test.@testset "Protein Assembly Type Detection" begin
                contigs_protein = Mycelia.assemble_genome(protein_records, method=Mycelia.KmerGraph, k=3)
                Test.@test contigs_protein isa Vector{String}
                Test.@test !isempty(contigs_protein)
            end
        end
        
        Test.@testset "Graph Method Validation" begin
            dna_records, fastq_records, _, _ = create_test_data()
            
            # Test that assembly method enum contains all 6 graph types
            Test.@test Mycelia.NgramGraph isa Mycelia.AssemblyMethod
            Test.@test Mycelia.KmerGraph isa Mycelia.AssemblyMethod
            Test.@test Mycelia.QualmerGraph isa Mycelia.AssemblyMethod
            Test.@test Mycelia.StringGraph isa Mycelia.AssemblyMethod
            Test.@test Mycelia.BioSequenceGraph isa Mycelia.AssemblyMethod
            Test.@test Mycelia.QualityBioSequenceGraph isa Mycelia.AssemblyMethod
            
            # Test that methods work with appropriate data
            Test.@test_nowarn Mycelia.assemble_genome(dna_records, method=Mycelia.KmerGraph, k=5)
            Test.@test_nowarn Mycelia.assemble_genome(fastq_records, method=Mycelia.QualmerGraph, k=5)
            Test.@test_nowarn Mycelia.assemble_genome(dna_records, method=Mycelia.BioSequenceGraph, k=5)
            Test.@test_nowarn Mycelia.assemble_genome(fastq_records, method=Mycelia.QualityBioSequenceGraph, k=5)
        end
    end
    
    Test.@testset "Memory and Performance Validation" begin
        
        Test.@testset "Fixed-Length vs Variable-Length Memory" begin
            dna_records, fastq_records, _, _ = create_test_data()
            
            # Test that fixed-length graphs use canonical representation
            k = 4
            kmer_graph = Mycelia.build_kmer_graph_next(Kmers.DNAKmer{k}, dna_records)
            
            # All vertices should be unique k-mers (canonical representation)
            unique_labels = Set(MetaGraphsNext.labels(kmer_graph))
            Test.@test length(unique_labels) == length(collect(MetaGraphsNext.labels(kmer_graph)))
            
            # Variable-length graphs should be simplified
            biosequence_graph = Mycelia.kmer_graph_to_biosequence_graph(kmer_graph)
            Test.@test length(collect(MetaGraphsNext.labels(biosequence_graph))) <= length(collect(MetaGraphsNext.labels(kmer_graph)))
        end
        
        Test.@testset "Quality Preservation Efficiency" begin
            _, fastq_records, _, _ = create_test_data()
            
            # Test that quality information is preserved efficiently
            quality_graph = Mycelia.build_quality_biosequence_graph(fastq_records)
            
            for label in MetaGraphsNext.labels(quality_graph)
                vertex_data = quality_graph[label]
                # Quality scores should match sequence length
                Test.@test length(vertex_data.quality_scores) == length(vertex_data.sequence)
                # Joint probability should be meaningful
                Test.@test vertex_data.joint_probability > 0.0
                Test.@test vertex_data.mean_quality > 0.0
            end
        end
    end
end