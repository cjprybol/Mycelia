# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/six_graph_hierarchy_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/six_graph_hierarchy_tests.jl", "test/4_assembly", execute=false)'
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
            
            ngram_graph = Mycelia.Rhizomorph.build_ngram_graph([test_string], n; dataset_id="test")
            
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
                dna_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, k; dataset_id="test", mode=:doublestrand)
                
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
                    Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                    Test.@test vertex_data.Kmer == label
                    evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                    Test.@test !isempty(evidence_entries)
                    Test.@test all(entry -> entry.strand in (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse), evidence_entries)
                end
            end
            
            # Test RNA k-mer graphs
            Test.@testset "RNA K-mer Graphs" begin
                k = 4
                kmer_type = Kmers.RNAKmer{k}
                rna_graph = Mycelia.Rhizomorph.build_kmer_graph(rna_records, k; dataset_id="test", mode=:doublestrand)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test rna_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa Kmers.RNAKmer{k}, MetaGraphsNext.labels(rna_graph))
                
                # Verify fixed-length vertices
                for label in MetaGraphsNext.labels(rna_graph)
                    Test.@test length(label) == k
                    Test.@test label isa Kmers.RNAKmer{k}
                    vertex_data = rna_graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                    Test.@test vertex_data.Kmer == label
                end
            end
            
            # Test amino acid k-mer graphs
            Test.@testset "Amino Acid K-mer Graphs" begin
                k = 3
                kmer_type = Kmers.AAKmer{k}
                aa_graph = Mycelia.Rhizomorph.build_kmer_graph(protein_records, k; dataset_id="test", mode=:singlestrand)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test aa_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa Kmers.AAKmer{k}, MetaGraphsNext.labels(aa_graph))
                
                # Verify fixed-length vertices
                for label in MetaGraphsNext.labels(aa_graph)
                    Test.@test length(label) == k
                    Test.@test label isa Kmers.AAKmer{k}
                    vertex_data = aa_graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                    Test.@test vertex_data.Kmer == label
                end
            end
        end
        
        Test.@testset "3. Qualmer Graphs - Quality-Aware Analysis (NO strings)" begin
            _, fastq_records, _, _ = create_test_data()
            
            # Test DNA qualmer graphs
            Test.@testset "DNA Qualmer Graphs" begin
                k = 5
                qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, k; dataset_id="test", mode=:doublestrand)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph
                
                # Verify fixed-length vertices with quality information
                for label in MetaGraphsNext.labels(qualmer_graph)
                    vertex_data = qualmer_graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.QualmerVertexData
                    Test.@test vertex_data.Kmer isa Kmers.DNAKmer{k}
                    Test.@test length(vertex_data.Kmer) == k
                    evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                    Test.@test !isempty(evidence_entries)
                    Test.@test all(entry -> length(entry.quality_scores) == k, evidence_entries)
                end
                
                # Test edge data with quality evidence
                for (src, dst) in MetaGraphsNext.edge_labels(qualmer_graph)
                    edge_data = qualmer_graph[src, dst]
                    Test.@test edge_data isa Mycelia.Rhizomorph.QualmerEdgeData
                    evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(edge_data.evidence)
                    Test.@test !isempty(evidence_entries)
                end
            end
        end
    end
    
    Test.@testset "Variable-Length Graph Types (Simplified Products)" begin
        
        Test.@testset "4. String Graphs - Simplified N-gram Graphs" begin
            # Create N-gram graph first
            test_string = "ABCDEFABCDEF"
            n = 3
            ngram_graph = Mycelia.Rhizomorph.build_ngram_graph([test_string], n; dataset_id="test")
            
            # Type checking
            Test.@test ngram_graph isa MetaGraphsNext.MetaGraph
            
            # Test assembly from string graph
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(ngram_graph)
            assembled_strings = [Mycelia.Rhizomorph.path_to_sequence(path, ngram_graph) for path in paths]
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
                biosequence_graph = Mycelia.Rhizomorph.build_fasta_graph(dna_records; dataset_id="test", min_overlap=3)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test biosequence_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(biosequence_graph))
                
                # Variable-length vertices
                for label in MetaGraphsNext.labels(biosequence_graph)
                    Test.@test label isa BioSequences.LongDNA{4}
                    Test.@test length(label) > 0  # Variable length
                    vertex_data = biosequence_graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.BioSequenceVertexData
                    Test.@test vertex_data.sequence == label
                end
            end
            
            # Test k-mer graph to BioSequence graph conversion
            Test.@testset "K-mer to BioSequence Graph Conversion" begin
                k = 5
                kmer_type = Kmers.DNAKmer{k}
                kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, k; dataset_id="test", mode=:doublestrand)
                
                # Convert to BioSequence graph
                biosequence_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(kmer_graph)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test biosequence_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(biosequence_graph))
                
                # Variable-length vertices from k-mer path simplification
                for label in MetaGraphsNext.labels(biosequence_graph)
                    Test.@test label isa BioSequences.LongDNA{4}
                    Test.@test length(label) >= k  # Should be at least k bases
                    vertex_data = biosequence_graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.BioSequenceVertexData
                    evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                    Test.@test !isempty(evidence_entries)
                end
            end
        end
        
        Test.@testset "6. FASTQ Graphs - Quality-Aware BioSequence Graphs (NO strings)" begin
            _, fastq_records, _, _ = create_test_data()
            
            # Test direct quality-aware BioSequence graph construction
            Test.@testset "Direct Quality BioSequence Graph from FASTQ" begin
                quality_graph = Mycelia.Rhizomorph.build_fastq_graph(fastq_records; dataset_id="test", min_overlap=3)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test quality_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(quality_graph))
                
                # Variable-length vertices with quality preservation
                for label in MetaGraphsNext.labels(quality_graph)
                    Test.@test label isa BioSequences.LongDNA{4}
                    Test.@test length(label) > 0  # Variable length
                    vertex_data = quality_graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
                    Test.@test vertex_data.sequence == label
                    evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                    Test.@test !isempty(evidence_entries)
                    Test.@test all(entry -> length(entry.quality_scores) == length(label), evidence_entries)
                end
            end
            
            # Test Qualmer graph to quality-aware BioSequence graph conversion
            Test.@testset "Qualmer to Quality BioSequence Graph Conversion" begin
                k = 5
                qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, k; dataset_id="test", mode=:doublestrand)
                
                # Convert to quality-aware BioSequence graph
                quality_biosequence_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_graph)
                
                # Type checking - NO STRING CONVERSIONS
                Test.@test quality_biosequence_graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(quality_biosequence_graph))
                
                # Variable-length vertices with quality preservation
                for label in MetaGraphsNext.labels(quality_biosequence_graph)
                    Test.@test label isa BioSequences.LongDNA{4}
                    Test.@test length(label) >= k  # Should be at least k bases
                    vertex_data = quality_biosequence_graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
                    Test.@test vertex_data.sequence == label
                    evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                    Test.@test !isempty(evidence_entries)
                end
            end
            
            # Test FASTQ evidence preservation
            Test.@testset "FASTQ Evidence Preservation" begin
                quality_graph = Mycelia.Rhizomorph.build_fastq_graph(fastq_records; dataset_id="test", min_overlap=3)
                
                # Quality preservation
                for label in MetaGraphsNext.labels(quality_graph)
                    vertex_data = quality_graph[label]
                    evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                    Test.@test !isempty(evidence_entries)
                    Test.@test all(entry -> length(entry.quality_scores) == length(label), evidence_entries)
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
            ngram_graph = Mycelia.Rhizomorph.build_ngram_graph(["ABCDEFABCDEF"], 3; dataset_id="test")
            Test.@test all(label -> label isa String, MetaGraphsNext.labels(ngram_graph))
            
            # 2. K-mer graphs → DNAKmer types
            kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, k; dataset_id="test", mode=:doublestrand)
            Test.@test all(label -> label isa Kmers.DNAKmer{k}, MetaGraphsNext.labels(kmer_graph))
            
            # 3. Qualmer graphs → Qualmer types
            qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, k; dataset_id="test", mode=:doublestrand)
            for label in MetaGraphsNext.labels(qualmer_graph)
                vertex_data = qualmer_graph[label]
                Test.@test vertex_data.Kmer isa Kmers.DNAKmer{k}
            end
            
            # 4. String graphs → strings (variable length)
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(ngram_graph)
            assembled_strings = [Mycelia.Rhizomorph.path_to_sequence(path, ngram_graph) for path in paths]
            Test.@test all(s -> s isa String, assembled_strings)
            
            # 5. BioSequence graphs → BioSequence types
            biosequence_graph = Mycelia.Rhizomorph.build_fasta_graph(dna_records; dataset_id="test", min_overlap=3)
            Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(biosequence_graph))
            
            # 6. Quality BioSequence graphs → BioSequence types with quality
            quality_graph = Mycelia.Rhizomorph.build_fastq_graph(fastq_records; dataset_id="test", min_overlap=3)
            Test.@test all(label -> label isa BioSequences.LongDNA{4}, MetaGraphsNext.labels(quality_graph))
            for label in MetaGraphsNext.labels(quality_graph)
                vertex_data = quality_graph[label]
                evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                Test.@test !isempty(evidence_entries)
            end
        end
        
        Test.@testset "Strand-Aware Functionality" begin
            dna_records, fastq_records, _, _ = create_test_data()
            
            # Test SingleStrand mode
            Test.@testset "SingleStrand Mode" begin
                k = 4
                kmer_type = Kmers.DNAKmer{k}
                single_strand_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, k; dataset_id="test", mode=:singlestrand)
                
                # Check strand information in edges
                for (src, dst) in MetaGraphsNext.edge_labels(single_strand_graph)
                    edge_data = single_strand_graph[src, dst]
                    Test.@test edge_data isa Mycelia.Rhizomorph.KmerEdgeData
                    evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(edge_data.evidence)
                    Test.@test !isempty(evidence_entries)
                    Test.@test all(entry -> entry.strand == Mycelia.Rhizomorph.Forward, evidence_entries)
                end
            end
            
            # Test DoubleStrand mode
            Test.@testset "DoubleStrand Mode" begin
                k = 4
                kmer_type = Kmers.DNAKmer{k}
                double_strand_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, k; dataset_id="test", mode=:doublestrand)
                
                # Check label consistency
                for label in MetaGraphsNext.labels(double_strand_graph)
                    vertex_data = double_strand_graph[label]
                    Test.@test vertex_data.Kmer == label
                end
            end
        end
    end
    
    Test.@testset "Assembly Pipeline Integration" begin
        
        Test.@testset "Automatic Type Detection" begin
            dna_records, fastq_records, rna_records, protein_records = create_test_data()
            
            # Test DNA assembly (should auto-detect method based on input format)
            Test.@testset "DNA Assembly Type Detection" begin
                # FASTQ input should auto-detect to Qualmer graphs
                result_fastq = Mycelia.Rhizomorph.assemble_genome(fastq_records; k=5)
                Test.@test result_fastq isa Mycelia.Rhizomorph.AssemblyResult
                Test.@test !isempty(result_fastq.contigs)

                # FASTA input should auto-detect to K-mer graphs
                result_fasta = Mycelia.Rhizomorph.assemble_genome(dna_records; k=5)
                Test.@test result_fasta isa Mycelia.Rhizomorph.AssemblyResult
                Test.@test !isempty(result_fasta.contigs)
            end

            # Test RNA assembly
            Test.@testset "RNA Assembly Type Detection" begin
                result_rna = Mycelia.Rhizomorph.assemble_genome(rna_records; k=4)
                Test.@test result_rna isa Mycelia.Rhizomorph.AssemblyResult
                Test.@test !isempty(result_rna.contigs)
            end

            # Test protein assembly (should auto-detect SingleStrand mode)
            Test.@testset "Protein Assembly Type Detection" begin
                result_protein = Mycelia.Rhizomorph.assemble_genome(protein_records; k=3)
                Test.@test result_protein isa Mycelia.Rhizomorph.AssemblyResult
                Test.@test !isempty(result_protein.contigs)
            end
        end
        
        Test.@testset "Graph Method Validation" begin
            dna_records, fastq_records, _, _ = create_test_data()
            
            # Test that assembly method enum contains all 6 graph types
            Test.@test Mycelia.Rhizomorph.NgramGraph isa Mycelia.Rhizomorph.AssemblyMethod
            Test.@test Mycelia.Rhizomorph.KmerGraph isa Mycelia.Rhizomorph.AssemblyMethod
            Test.@test Mycelia.Rhizomorph.QualmerGraph isa Mycelia.Rhizomorph.AssemblyMethod
            Test.@test Mycelia.Rhizomorph.StringGraph isa Mycelia.Rhizomorph.AssemblyMethod
            Test.@test Mycelia.Rhizomorph.BioSequenceGraph isa Mycelia.Rhizomorph.AssemblyMethod
            Test.@test Mycelia.Rhizomorph.QualityBioSequenceGraph isa Mycelia.Rhizomorph.AssemblyMethod
            
            # Test that auto-detection works without errors for different input types
            Test.@test_nowarn Mycelia.Rhizomorph.assemble_genome(dna_records; k=5)  # Auto-detects k-mer graph for FASTA
            Test.@test_nowarn Mycelia.Rhizomorph.assemble_genome(fastq_records; k=5)  # Auto-detects qualmer graph for FASTQ
            Test.@test_nowarn Mycelia.Rhizomorph.assemble_genome(dna_records; min_overlap=10)  # Auto-detects BioSequence graph
            Test.@test_nowarn Mycelia.Rhizomorph.assemble_genome(fastq_records; min_overlap=10)  # Auto-detects Quality BioSequence graph
        end
    end
    
    Test.@testset "Memory and Performance Validation" begin
        
        Test.@testset "Fixed-Length vs Variable-Length Memory" begin
            dna_records, fastq_records, _, _ = create_test_data()
            
            # Test that canonical fixed-length graphs use canonical representation
            k = 4
            kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, k; dataset_id="test", mode=:canonical)
            
            # All vertices should be unique k-mers (canonical representation)
            unique_labels = Set(MetaGraphsNext.labels(kmer_graph))
            Test.@test length(unique_labels) == length(collect(MetaGraphsNext.labels(kmer_graph)))
            
            # Variable-length graphs should be simplified
            biosequence_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(kmer_graph)
            Test.@test length(collect(MetaGraphsNext.labels(biosequence_graph))) <= length(collect(MetaGraphsNext.labels(kmer_graph)))
        end
        
        Test.@testset "Quality Preservation Efficiency" begin
            _, fastq_records, _, _ = create_test_data()
            
            # Test that quality information is preserved efficiently
            quality_graph = Mycelia.Rhizomorph.build_fastq_graph(fastq_records; dataset_id="test", min_overlap=3)
            
            for label in MetaGraphsNext.labels(quality_graph)
                vertex_data = quality_graph[label]
                evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                Test.@test !isempty(evidence_entries)
                Test.@test all(entry -> length(entry.quality_scores) == length(vertex_data.sequence), evidence_entries)

                dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
                if !isempty(dataset_ids)
                    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, first(dataset_ids))
                    Test.@test mean_quality !== nothing
                    Test.@test all(q -> q >= 0.0, mean_quality)
                end
            end
        end
    end
end
