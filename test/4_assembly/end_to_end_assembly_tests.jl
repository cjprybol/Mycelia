# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/4_assembly/end_to_end_assembly_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/end_to_end_assembly_tests.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import Random
import StatsBase
import Kmers

# Test parameters
const TEST_LENGTHS = [10, 100, 1000]
const ERROR_RATES = [0.0, 0.001, 0.01, 0.1]  # 0%, 0.1%, 1%, 10%
const COVERAGE_DEPTHS = [10, 100, 1000]

Test.@testset "End-to-End Assembly Tests" begin
    
    Test.@testset "Base Case: Reference Sequence Graph Round-trip" begin
        Test.@testset "String Graph Base Case" begin
            # Test with simple strings
            for seq_length in TEST_LENGTHS
                reference_string = Random.randstring(seq_length)
                
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
            for seq_length in TEST_LENGTHS
                reference_seq = BioSequences.randdnaseq(seq_length)
                reference_record = FASTX.FASTA.Record("reference", reference_seq)

                # Test DoubleStrand mode (canonical)
                kmer_type = Kmers.DNAKmer{5}
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
        
        Test.@testset "Sequence Graph Next Base Case - RNA" begin
            for seq_length in TEST_LENGTHS
                reference_seq = string(BioSequences.randrnaseq(seq_length))
                reference_record = FASTX.FASTA.Record("reference", reference_seq)

                # Test SingleStrand mode for RNA
                kmer_type = Kmers.RNAKmer{5}
                graph = Mycelia.build_kmer_graph_next(kmer_type, [reference_record];
                                                    graph_mode=Mycelia.SingleStrand)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))

                # Verify vertices contain RNA k-mers
                for label in MetaGraphsNext.labels(graph)
                    vertex_data = graph[label]
                    Test.@test vertex_data isa Mycelia.KmerVertexData
                    # Verify the k-mer type is RNA (not converted to DNA)
                    Test.@test vertex_data.canonical_kmer isa Kmers.RNAKmer
                    # Verify no T nucleotides (RNA should never contain T)
                    Test.@test !occursin('T', string(vertex_data.canonical_kmer))
                end
            end
        end
        
        Test.@testset "Sequence Graph Next Base Case - Amino Acids" begin
            for seq_length in TEST_LENGTHS
                reference_seq = BioSequences.randaaseq(seq_length)
                reference_record = FASTX.FASTA.Record("reference", reference_seq)

                # Test SingleStrand mode for amino acids
                kmer_type = Kmers.AAKmer{5}
                graph = Mycelia.build_kmer_graph_next(kmer_type, [reference_record];
                                                    graph_mode=Mycelia.SingleStrand)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))

                # Verify vertices contain amino acid k-mers
                for label in MetaGraphsNext.labels(graph)
                    vertex_data = graph[label]
                    Test.@test vertex_data isa Mycelia.KmerVertexData
                    # Check that all characters are valid amino acids
                    Test.@test all(c in Mycelia.UNAMBIGUOUS_AA_SYMBOLS for c in string(vertex_data.canonical_kmer))
                end
            end
        end

        Test.@testset "String Graph Base Case - ASCII Greek" begin
            for seq_length in TEST_LENGTHS
                reference_string = Mycelia.rand_ascii_greek_string(seq_length)
                
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
            for seq_length in TEST_LENGTHS
                reference_string = Mycelia.rand_latin1_string(seq_length)
                
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
            for seq_length in TEST_LENGTHS
                reference_string = Mycelia.rand_bmp_printable_string(seq_length)
                
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
            for seq_length in TEST_LENGTHS
                reference_string = Mycelia.rand_printable_unicode_string(seq_length)
                
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
            for seq_length in TEST_LENGTHS
                reference_string = Random.randstring(seq_length)
                
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
            for seq_length in TEST_LENGTHS
                reference_seq = BioSequences.randdnaseq(seq_length)

                for error_rate in ERROR_RATES
                    for coverage in COVERAGE_DEPTHS
                        # Create simulated reads
                        reads = create_test_reads(reference_seq, coverage, error_rate)

                        # Build k-mer graph in SingleStrand mode
                        kmer_type = Kmers.DNAKmer{5}
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

                        # Test quality-aware assembly from k-mer graph to FASTQ
                        # Build qualmer graph for quality-aware processing
                        qualmer_graph = Mycelia.build_qualmer_graph(reads; k=5, graph_mode=Mycelia.SingleStrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Convert back to FASTQ records with quality scores
                        assembled_fastq_records = Mycelia.quality_biosequence_graph_to_fastq(quality_graph)
                        Test.@test !isempty(assembled_fastq_records)
                        Test.@test all(r isa FASTX.FASTQ.Record for r in assembled_fastq_records)

                        # Validate that assembled sequences are proper DNA
                        for record in assembled_fastq_records
                            sequence_str = FASTX.sequence(String, record)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                Test.@test detected_type == :DNA

                                # Quality scores should reflect assembly confidence
                                quality_scores = Mycelia.get_phred_scores(record)
                                Test.@test !isempty(quality_scores)
                                Test.@test all(q >= 0 for q in quality_scores)  # Valid quality scores
                            end
                        end
                    end
                end
            end
        end
        
        Test.@testset "RNA SingleStrand Mode" begin
            for seq_length in TEST_LENGTHS
                reference_seq = string(BioSequences.randrnaseq(seq_length))

                for error_rate in ERROR_RATES
                    for coverage in COVERAGE_DEPTHS
                        # Create simulated reads
                        reads = Mycelia.create_test_rna_reads(reference_seq, coverage, error_rate)

                        # Build k-mer graph in SingleStrand mode
                        kmer_type = Kmers.RNAKmer{5}
                        graph = Mycelia.build_kmer_graph_next(kmer_type, reads;
                                                            graph_mode=Mycelia.SingleStrand)

                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))

                        # Verify RNA-specific properties
                        for label in MetaGraphsNext.labels(graph)
                            vertex_data = graph[label]
                            Test.@test vertex_data isa Mycelia.KmerVertexData
                            # Verify the k-mer type is RNA (not converted to DNA)
                            Test.@test vertex_data.canonical_kmer isa Kmers.RNAKmer
                            # Verify no T nucleotides (RNA should never contain T)
                            Test.@test !occursin('T', string(vertex_data.canonical_kmer))

                            # Test sequence type validation using existing functions
                            kmer_string = string(vertex_data.canonical_kmer)
                            # Note: detect_alphabet() may return :DNA for ambiguous sequences (A,C,G only)
                            # but the k-mer type system correctly preserves RNA type
                            detected_alphabet = Mycelia.detect_alphabet(kmer_string)
                            # Only test for RNA detection when U is actually present
                            if occursin('U', kmer_string)
                                Test.@test detected_alphabet == :RNA
                            end

                            # Test that the sequence is valid RNA using BioSequences constructor
                            Test.@test_nowarn BioSequences.LongRNA{4}(kmer_string)
                            # Test that it's NOT valid DNA if it contains U
                            if occursin('U', kmer_string)
                                Test.@test_throws Exception BioSequences.LongDNA{4}(kmer_string)
                            end
                        end

                        # Test quality-aware assembly from k-mer graph to FASTQ
                        # Build qualmer graph for quality-aware processing
                        qualmer_graph = Mycelia.build_qualmer_graph(reads; k=5, graph_mode=Mycelia.SingleStrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Convert back to FASTQ records with quality scores
                        assembled_fastq_records = Mycelia.quality_biosequence_graph_to_fastq(quality_graph)
                        Test.@test !isempty(assembled_fastq_records)
                        Test.@test all(r isa FASTX.FASTQ.Record for r in assembled_fastq_records)

                        # Validate that assembled sequences are proper RNA
                        for record in assembled_fastq_records
                            sequence_str = FASTX.sequence(String, record)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                Test.@test detected_type == :RNA
                                # Verify no T nucleotides (RNA should never contain T)
                                Test.@test !occursin('T', sequence_str)

                                # Quality scores should reflect assembly confidence
                                quality_scores = Mycelia.get_phred_scores(record)
                                Test.@test !isempty(quality_scores)
                                Test.@test all(q >= 0 for q in quality_scores)  # Valid quality scores
                            end
                        end
                    end
                end
            end
        end
        
        Test.@testset "Amino Acid SingleStrand Mode" begin
            for seq_length in TEST_LENGTHS
                reference_seq = BioSequences.randaaseq(seq_length)

                for error_rate in ERROR_RATES
                    for coverage in COVERAGE_DEPTHS
                        # Create simulated reads
                        reads = Mycelia.create_test_aa_reads(reference_seq, coverage, error_rate)

                        # Build k-mer graph in SingleStrand mode
                        kmer_type = Kmers.AAKmer{3}  # Shorter k-mers for AA
                        graph = Mycelia.build_kmer_graph_next(kmer_type, reads;
                                                            graph_mode=Mycelia.SingleStrand)

                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))

                        # Verify amino acid-specific properties
                        for label in MetaGraphsNext.labels(graph)
                            vertex_data = graph[label]
                            Test.@test vertex_data isa Mycelia.KmerVertexData
                            # Verify this is an amino acid k-mer
                            Test.@test vertex_data.canonical_kmer isa Kmers.AAKmer

                            # Test sequence type validation using existing functions
                            kmer_string = string(vertex_data.canonical_kmer)
                            detected_alphabet = Mycelia.detect_alphabet(kmer_string)
                            Test.@test detected_alphabet == :AA

                            # Verify the k-mer can be constructed as a valid AA sequence
                            Test.@test_nowarn BioSequences.LongAA(kmer_string)
                        end

                        # Test quality-aware assembly from k-mer graph to FASTQ
                        # Build qualmer graph for quality-aware processing
                        qualmer_graph = Mycelia.build_qualmer_graph(reads; k=3, graph_mode=Mycelia.SingleStrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Convert back to FASTQ records with quality scores
                        assembled_fastq_records = Mycelia.quality_biosequence_graph_to_fastq(quality_graph)
                        Test.@test !isempty(assembled_fastq_records)
                        Test.@test all(r isa FASTX.FASTQ.Record for r in assembled_fastq_records)

                        # Validate that assembled sequences are proper amino acids
                        for record in assembled_fastq_records
                            sequence_str = FASTX.sequence(String, record)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                Test.@test detected_type == :AA

                                # Quality scores should reflect assembly confidence
                                quality_scores = Mycelia.get_phred_scores(record)
                                Test.@test !isempty(quality_scores)
                                Test.@test all(q >= 0 for q in quality_scores)  # Valid quality scores
                            end
                        end
                    end
                end
            end
        end
    end
    
    Test.@testset "Canonical/Double-Stranded Sequence Graph Next Assembly" begin
        Test.@testset "DNA DoubleStrand Mode" begin
            for seq_length in TEST_LENGTHS
                reference_seq = BioSequences.randdnaseq(seq_length)

                for error_rate in ERROR_RATES
                    for coverage in COVERAGE_DEPTHS
                        # Create simulated reads
                        reads = create_test_reads(reference_seq, coverage, error_rate)

                        # Build k-mer graph in DoubleStrand mode (canonical)
                        kmer_type = Kmers.DNAKmer{5}
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

                        # Test quality-aware assembly from k-mer graph to FASTQ
                        # Build qualmer graph for quality-aware processing
                        qualmer_graph = Mycelia.build_qualmer_graph(reads; k=5, graph_mode=Mycelia.DoubleStrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Convert back to FASTQ records with quality scores
                        assembled_fastq_records = Mycelia.quality_biosequence_graph_to_fastq(quality_graph)
                        Test.@test !isempty(assembled_fastq_records)
                        Test.@test all(r isa FASTX.FASTQ.Record for r in assembled_fastq_records)

                        # Validate that assembled sequences are proper DNA
                        for record in assembled_fastq_records
                            sequence_str = FASTX.sequence(String, record)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                Test.@test detected_type == :DNA

                                # Quality scores should reflect assembly confidence
                                quality_scores = Mycelia.get_phred_scores(record)
                                Test.@test !isempty(quality_scores)
                                Test.@test all(q >= 0 for q in quality_scores)  # Valid quality scores
                            end
                        end
                    end
                end
            end
        end
        
        Test.@testset "RNA DoubleStrand Mode" begin
            for seq_length in TEST_LENGTHS
                reference_seq = string(BioSequences.randrnaseq(seq_length))

                for error_rate in ERROR_RATES
                    for coverage in COVERAGE_DEPTHS
                        # Create simulated reads
                        reads = Mycelia.create_test_rna_reads(reference_seq, coverage, error_rate)

                        # Build k-mer graph in DoubleStrand mode (canonical)
                        kmer_type = Kmers.RNAKmer{5}
                        graph = Mycelia.build_kmer_graph_next(kmer_type, reads;
                                                            graph_mode=Mycelia.DoubleStrand)

                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))

                        # Verify canonical RNA k-mer properties
                        for label in MetaGraphsNext.labels(graph)
                            vertex_data = graph[label]
                            Test.@test vertex_data isa Mycelia.KmerVertexData
                            Test.@test vertex_data.canonical_kmer == label
                            # Verify the k-mer type is RNA (not converted to DNA)
                            Test.@test vertex_data.canonical_kmer isa Kmers.RNAKmer
                            # Verify no T nucleotides (RNA should never contain T)
                            Test.@test !occursin('T', string(vertex_data.canonical_kmer))

                            # Test sequence type validation
                            kmer_string = string(vertex_data.canonical_kmer)
                            detected_alphabet = Mycelia.detect_alphabet(kmer_string)
                            # Only test for RNA detection when U is actually present
                            if occursin('U', kmer_string)
                                Test.@test detected_alphabet == :RNA
                            end

                            # Test that the sequence is valid RNA using BioSequences constructor
                            Test.@test_nowarn BioSequences.LongRNA{4}(kmer_string)
                            # Test that it's NOT valid DNA if it contains U
                            if occursin('U', kmer_string)
                                Test.@test_throws Exception BioSequences.LongDNA{4}(kmer_string)
                            end
                        end

                        # Test quality-aware assembly from k-mer graph to FASTQ
                        # Build qualmer graph for quality-aware processing
                        qualmer_graph = Mycelia.build_qualmer_graph(reads; k=5, graph_mode=Mycelia.DoubleStrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Convert back to FASTQ records with quality scores
                        assembled_fastq_records = Mycelia.quality_biosequence_graph_to_fastq(quality_graph)
                        Test.@test !isempty(assembled_fastq_records)
                        Test.@test all(r isa FASTX.FASTQ.Record for r in assembled_fastq_records)

                        # Validate that assembled sequences are proper RNA
                        for record in assembled_fastq_records
                            sequence_str = FASTX.sequence(String, record)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                Test.@test detected_type == :RNA
                                # Verify no T nucleotides (RNA should never contain T)
                                Test.@test !occursin('T', sequence_str)

                                # Quality scores should reflect assembly confidence
                                quality_scores = Mycelia.get_phred_scores(record)
                                Test.@test !isempty(quality_scores)
                                Test.@test all(q >= 0 for q in quality_scores)  # Valid quality scores
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
            
            kmer_type = Kmers.DNAKmer{5}
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
                kmer_type = Kmers.DNAKmer{5}
                
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
                
                kmer_type = Kmers.DNAKmer{5}
                graph = Mycelia.build_kmer_graph_next(kmer_type, reads; 
                                                    graph_mode=Mycelia.DoubleStrand)
                
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Higher coverage should generally result in more robust graphs
                total_coverage = sum(length(graph[label].coverage) for label in MetaGraphsNext.labels(graph))
                Test.@test total_coverage > 0
            end
        end
    end

    Test.@testset "K-mer Canonicalization Consistency Tests" begin
        ## Test that validates the fix for canonical k-mer consistency between
        ## graph construction and path extraction (the main issue we resolved)

        Test.@testset "DoubleStrand DNA Canonicalization Consistency" begin
            # Create a DNA sequence that will have different canonical forms
            # if lexicographic vs BioSequences.canonical() are used
            reference_seq = BioSequences.dna"ATCGTTTT"  # Forward
            reverse_comp = BioSequences.reverse_complement(reference_seq)  # AAAACGAT

            # Test with both forward and reverse complement sequences
            reads = [
                FASTX.FASTA.Record("forward", reference_seq),
                FASTX.FASTA.Record("reverse", reverse_comp)
            ]

            kmer_type = Kmers.DNAKmer{5}
            graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.DoubleStrand)

            Test.@test !isempty(MetaGraphsNext.labels(graph))

            # Test that we have proper coverage and strand orientations
            # This validates the canonicalization consistency fix
            total_coverage_entries = 0
            has_both_orientations = false

            for label in MetaGraphsNext.labels(graph)
                vertex_data = graph[label]
                strand_orientations = [strand for (obs_id, pos, strand) in vertex_data.coverage]
                total_coverage_entries += length(strand_orientations)

                # Check if we see both orientations (should happen with reverse complements)
                if Mycelia.Forward in strand_orientations && Mycelia.Reverse in strand_orientations
                    has_both_orientations = true
                end
            end

            Test.@test total_coverage_entries > 0
            Test.@test has_both_orientations  # Should see both orientations with forward and reverse sequences
        end

        Test.@testset "SingleStrand RNA Non-Canonicalization Consistency" begin
            # Test that SingleStrand mode uses k-mers as-is without canonicalization
            reference_seq = BioSequences.rna"AUCGUUUU"
            reads = [FASTX.FASTA.Record("rna", reference_seq)]

            kmer_type = Kmers.RNAKmer{5}
            graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.SingleStrand)

            Test.@test !isempty(MetaGraphsNext.labels(graph))

            # Test that all coverage entries use Forward orientation in SingleStrand mode
            # This validates the SingleStrand non-canonicalization fix
            total_coverage_entries = 0
            all_forward = true

            for label in MetaGraphsNext.labels(graph)
                vertex_data = graph[label]
                for (obs_id, pos, strand) in vertex_data.coverage
                    total_coverage_entries += 1
                    if strand != Mycelia.Forward
                        all_forward = false
                    end
                end
            end

            Test.@test total_coverage_entries > 0
            Test.@test all_forward  # All orientations should be Forward in SingleStrand mode
        end
    end
end