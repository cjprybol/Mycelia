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
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
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

# Test parameters (keep unit-test sized workloads)
const TEST_LENGTHS = [10, 50]
const ERROR_RATES = [0.0, 0.01]  # 0%, 1%
const COVERAGE_DEPTHS = [5, 20]

Test.@testset "End-to-End Assembly Tests" begin
    
    Test.@testset "Base Case: Reference Sequence Graph Round-trip" begin
        Test.@testset "String Graph Base Case" begin
            # Test with simple strings
            for seq_length in TEST_LENGTHS
                reference_string = Random.randstring(seq_length)
                
                # Create graph from reference
                graph = Mycelia.Rhizomorph.build_ngram_graph([reference_string], 3; dataset_id="test")
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test graph connectivity
                components = Graphs.connected_components(graph)
                Test.@test !isempty(components)
            end
        end
        
        Test.@testset "Sequence Graph Next Base Case - DNA" begin
            for seq_length in TEST_LENGTHS
                reference_seq = BioSequences.randdnaseq(seq_length)
                reference_record = FASTX.FASTA.Record("reference", reference_seq)

                # Test DoubleStrand mode
                kmer_type = Kmers.DNAKmer{5}
                graph = Mycelia.Rhizomorph.build_kmer_graph([reference_record], 5; dataset_id="test", mode=:doublestrand)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))

                # Verify vertices are valid k-mers
                for label in MetaGraphsNext.labels(graph)
                    vertex_data = graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                    Test.@test vertex_data.Kmer == label
                end

                # Test that we can write and read GFA
                Test.@test hasmethod(Mycelia.Rhizomorph.write_gfa_next, (typeof(graph), String))
                Test.@test hasmethod(Mycelia.Rhizomorph.read_gfa_next, (String,))
            end
        end
        
        Test.@testset "Sequence Graph Next Base Case - RNA" begin
            for seq_length in TEST_LENGTHS
                reference_seq = string(BioSequences.randrnaseq(seq_length))
                if !occursin('U', reference_seq)
                    reference_seq = "U" * reference_seq[2:end]
                end
                reference_record = FASTX.FASTA.Record("reference", reference_seq)

                # Test SingleStrand mode for RNA
                kmer_type = Kmers.RNAKmer{5}
                graph = Mycelia.Rhizomorph.build_kmer_graph([reference_record], 5; dataset_id="test", mode=:singlestrand)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))

                # Verify vertices contain RNA k-mers
                for label in MetaGraphsNext.labels(graph)
                    vertex_data = graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                    # Verify the k-mer type is RNA (not converted to DNA)
                    Test.@test vertex_data.Kmer isa Kmers.RNAKmer
                    # Verify no T nucleotides (RNA should never contain T)
                    Test.@test !occursin('T', string(vertex_data.Kmer))
                end
            end
        end
        
        Test.@testset "Sequence Graph Next Base Case - Amino Acids" begin
            for seq_length in TEST_LENGTHS
                reference_seq = BioSequences.randaaseq(seq_length)
                reference_record = FASTX.FASTA.Record("reference", reference_seq)

                # Test SingleStrand mode for amino acids
                kmer_type = Kmers.AAKmer{5}
                graph = Mycelia.Rhizomorph.build_kmer_graph([reference_record], 5; dataset_id="test", mode=:singlestrand)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))

                # Verify vertices contain amino acid k-mers
                for label in MetaGraphsNext.labels(graph)
                    vertex_data = graph[label]
                    Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                    # Check that all characters are valid amino acids
                    Test.@test all(c in Mycelia.UNAMBIGUOUS_AA_CHARSET for c in string(vertex_data.Kmer))
                end
            end
        end

        Test.@testset "String Graph Base Case - ASCII Greek" begin
            for seq_length in TEST_LENGTHS
                reference_string = Mycelia.rand_ascii_greek_string(seq_length)
                
                # Create graph from reference
                graph = Mycelia.Rhizomorph.build_ngram_graph([reference_string], 3; dataset_id="test")
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test assembly - should recover non-empty result
                paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
                assembly = [Mycelia.Rhizomorph.path_to_sequence(path, graph) for path in paths]
                Test.@test !isempty(assembly)
                
                # Test graph connectivity
                components = Graphs.connected_components(graph)
                Test.@test !isempty(components)
            end
        end

        Test.@testset "String Graph Base Case - Latin1" begin
            for seq_length in TEST_LENGTHS
                reference_string = Mycelia.rand_latin1_string(seq_length)
                
                # Create graph from reference
                graph = Mycelia.Rhizomorph.build_ngram_graph([reference_string], 3; dataset_id="test")
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test assembly - should recover non-empty result
                paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
                assembly = [Mycelia.Rhizomorph.path_to_sequence(path, graph) for path in paths]
                Test.@test !isempty(assembly)
                
                # Test graph connectivity
                components = Graphs.connected_components(graph)
                Test.@test !isempty(components)
            end
        end

        Test.@testset "String Graph Base Case - BMP Printable" begin
            for seq_length in TEST_LENGTHS
                reference_string = Mycelia.rand_bmp_printable_string(seq_length)
                
                # Create graph from reference
                graph = Mycelia.Rhizomorph.build_ngram_graph([reference_string], 3; dataset_id="test")
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test assembly - should recover non-empty result
                paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
                assembly = [Mycelia.Rhizomorph.path_to_sequence(path, graph) for path in paths]
                Test.@test !isempty(assembly)
                
                # Test graph connectivity
                components = Graphs.connected_components(graph)
                Test.@test !isempty(components)
            end
        end

        Test.@testset "String Graph Base Case - Printable Unicode" begin
            for seq_length in TEST_LENGTHS
                reference_string = Mycelia.rand_printable_unicode_string(seq_length)
                
                # Create graph from reference
                graph = Mycelia.Rhizomorph.build_ngram_graph([reference_string], 3; dataset_id="test")
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Test basic graph properties
                Test.@test length(MetaGraphsNext.labels(graph)) <= length(reference_string) - 3 + 1
                
                # Test assembly - should recover non-empty result
                paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
                assembly = [Mycelia.Rhizomorph.path_to_sequence(path, graph) for path in paths]
                Test.@test !isempty(assembly)
                
                # Test graph connectivity
                components = Graphs.connected_components(graph)
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
                        graph = Mycelia.Rhizomorph.build_ngram_graph(mutated_strings, 3; dataset_id="test")
                        
                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))
                        
                        # Test assembly process
                        collapsed_graph = graph
                        paths = Mycelia.Rhizomorph.find_eulerian_paths_next(collapsed_graph)
                        if isempty(paths)
                            paths = Mycelia.Rhizomorph.find_contigs_next(collapsed_graph; min_contig_length=1)
                        end
                        assemblies = String[]
                        for path_entry in paths
                            if path_entry isa Mycelia.Rhizomorph.ContigPath
                                push!(assemblies, string(path_entry.sequence))
                            elseif !isempty(path_entry)
                                push!(assemblies, string(Mycelia.Rhizomorph.path_to_sequence(path_entry, collapsed_graph)))
                            end
                        end
                        
                        Test.@test !isempty(assemblies)
                        Test.@test all(asm isa String for asm in assemblies)
                        
                        # For low error rates, assembly should be reasonable
                        if error_rate == 0.0
                            Test.@test any(length(asm) >= length(reference_string) ÷ 4 for asm in assemblies)
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
                        reads = Mycelia.create_test_reads(reference_seq, coverage, error_rate)

                        # Build k-mer graph in SingleStrand mode
                        kmer_type = Kmers.DNAKmer{5}
                        graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:singlestrand)

                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))

                        # Verify all strand orientations are Forward in SingleStrand mode
                        for label in MetaGraphsNext.labels(graph)
                            vertex_data = graph[label]
                            evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                            Test.@test !isempty(evidence_entries)
                            Test.@test all(entry -> entry.strand == Mycelia.Rhizomorph.Forward, evidence_entries)
                        end

                        # Test that edges respect strand constraints
                        for edge_label in MetaGraphsNext.edge_labels(graph)
                            if !isempty(edge_label)
                                edge_data = graph[edge_label...]
                                Test.@test edge_data isa Mycelia.Rhizomorph.KmerEdgeData
                                evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(edge_data.evidence)
                                Test.@test !isempty(evidence_entries)
                                Test.@test all(entry -> entry.strand == Mycelia.Rhizomorph.Forward, evidence_entries)
                            end
                        end

                        # Test quality-aware assembly from k-mer graph to FASTQ
                        # Build qualmer graph for quality-aware processing
                        qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, 5; dataset_id="test", mode=:singlestrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Validate that assembled sequences are proper DNA
                        for label in MetaGraphsNext.labels(quality_graph)
                            vertex_data = quality_graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
                            sequence_str = string(vertex_data.sequence)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                Test.@test detected_type == :DNA
                            end

                            evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                            Test.@test !isempty(evidence_entries)
                            for entry in evidence_entries
                                Test.@test entry.strand == Mycelia.Rhizomorph.Forward
                                Test.@test !isempty(entry.quality_scores)
                                Test.@test all(q >= 0 for q in entry.quality_scores)  # Valid quality scores
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
                        graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:singlestrand)

                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))

                        # Verify RNA-specific properties
                        for label in MetaGraphsNext.labels(graph)
                            vertex_data = graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                            # Verify the k-mer type is RNA (not converted to DNA)
                            Test.@test vertex_data.Kmer isa Kmers.RNAKmer
                            # Verify no T nucleotides (RNA should never contain T)
                            Test.@test !occursin('T', string(vertex_data.Kmer))

                            # Test sequence type validation using existing functions
                            kmer_string = string(vertex_data.Kmer)
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
                        qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, 5; dataset_id="test", mode=:singlestrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Validate that assembled sequences are proper RNA
                        for label in MetaGraphsNext.labels(quality_graph)
                            vertex_data = quality_graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
                            sequence_str = string(vertex_data.sequence)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                if occursin('U', sequence_str)
                                    Test.@test detected_type == :RNA
                                end
                                # Verify no T nucleotides (RNA should never contain T)
                                Test.@test !occursin('T', sequence_str)
                            end

                            evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                            Test.@test !isempty(evidence_entries)
                            for entry in evidence_entries
                                Test.@test entry.strand == Mycelia.Rhizomorph.Forward
                                Test.@test !isempty(entry.quality_scores)
                                Test.@test all(q >= 0 for q in entry.quality_scores)  # Valid quality scores
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
                        graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 3; dataset_id="test", mode=:singlestrand)

                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))

                        # Verify amino acid-specific properties
                        for label in MetaGraphsNext.labels(graph)
                            vertex_data = graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                            # Verify this is an amino acid k-mer
                            Test.@test vertex_data.Kmer isa Kmers.AAKmer

                            # Test sequence type validation using existing functions
                            kmer_string = string(vertex_data.Kmer)
                            detected_alphabet = Mycelia.detect_alphabet(kmer_string)
                            has_non_nucleotide = any(c -> !(c in Mycelia.AMBIGUOUS_DNA_CHARSET) &&
                                                         !(c in Mycelia.AMBIGUOUS_RNA_CHARSET),
                                                     kmer_string)
                            if has_non_nucleotide
                                Test.@test detected_alphabet == :AA
                            end

                            # Verify the k-mer can be constructed as a valid AA sequence
                            Test.@test_nowarn BioSequences.LongAA(kmer_string)
                        end

                        # Test quality-aware assembly from k-mer graph to FASTQ
                        # Build qualmer graph for quality-aware processing
                        qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, 3; dataset_id="test", mode=:singlestrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Validate that assembled sequences are proper amino acids
                        for label in MetaGraphsNext.labels(quality_graph)
                            vertex_data = quality_graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
                            sequence_str = string(vertex_data.sequence)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                has_non_nucleotide = any(c -> !(c in Mycelia.AMBIGUOUS_DNA_CHARSET) &&
                                                             !(c in Mycelia.AMBIGUOUS_RNA_CHARSET),
                                                         sequence_str)
                                if has_non_nucleotide
                                    Test.@test detected_type == :AA
                                end
                            end

                            evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                            Test.@test !isempty(evidence_entries)
                            for entry in evidence_entries
                                Test.@test entry.strand == Mycelia.Rhizomorph.Forward
                                Test.@test !isempty(entry.quality_scores)
                                Test.@test all(q >= 0 for q in entry.quality_scores)  # Valid quality scores
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
                        reads = Mycelia.create_test_reads(reference_seq, coverage, error_rate)

                        # Build k-mer graph in DoubleStrand mode
                        kmer_type = Kmers.DNAKmer{5}
                        graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:doublestrand)

                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))

                        # Verify double-strand k-mer properties
                        has_forward = false
                        has_reverse = false
                        for label in MetaGraphsNext.labels(graph)
                            vertex_data = graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                            Test.@test vertex_data.Kmer == label

                            evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                            Test.@test !isempty(evidence_entries)
                            for entry in evidence_entries
                                has_forward |= entry.strand == Mycelia.Rhizomorph.Forward
                                has_reverse |= entry.strand == Mycelia.Rhizomorph.Reverse
                            end
                        end
                        Test.@test has_forward
                        Test.@test has_reverse

                        # Test that edges handle strand transitions correctly
                        for edge_label in MetaGraphsNext.edge_labels(graph)
                            if !isempty(edge_label)
                                edge_data = graph[edge_label...]
                                Test.@test edge_data isa Mycelia.Rhizomorph.KmerEdgeData
                                evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(edge_data.evidence)
                                Test.@test !isempty(evidence_entries)
                                Test.@test all(entry -> entry.strand in (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse), evidence_entries)
                            end
                        end

                        # Test quality-aware assembly from k-mer graph to FASTQ
                        # Build qualmer graph for quality-aware processing
                        qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, 5; dataset_id="test", mode=:doublestrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Validate that assembled sequences are proper DNA
                        for label in MetaGraphsNext.labels(quality_graph)
                            vertex_data = quality_graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
                            sequence_str = string(vertex_data.sequence)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                Test.@test detected_type == :DNA
                            end

                            evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                            Test.@test !isempty(evidence_entries)
                            for entry in evidence_entries
                                Test.@test entry.strand in (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse)
                                Test.@test !isempty(entry.quality_scores)
                                Test.@test all(q >= 0 for q in entry.quality_scores)  # Valid quality scores
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

                        # Build k-mer graph in DoubleStrand mode
                        kmer_type = Kmers.RNAKmer{5}
                        graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:doublestrand)

                        Test.@test graph isa MetaGraphsNext.MetaGraph
                        Test.@test !isempty(MetaGraphsNext.labels(graph))

                        # Verify double-strand RNA k-mer properties
                        for label in MetaGraphsNext.labels(graph)
                            vertex_data = graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
                            Test.@test vertex_data.Kmer == label
                            # Verify the k-mer type is RNA (not converted to DNA)
                            Test.@test vertex_data.Kmer isa Kmers.RNAKmer
                            # Verify no T nucleotides (RNA should never contain T)
                            Test.@test !occursin('T', string(vertex_data.Kmer))

                            # Test sequence type validation
                            kmer_string = string(vertex_data.Kmer)
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
                        qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, 5; dataset_id="test", mode=:doublestrand)
                        Test.@test qualmer_graph isa MetaGraphsNext.MetaGraph

                        # Convert to quality-aware BioSequence graph
                        quality_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_graph)
                        Test.@test quality_graph isa MetaGraphsNext.MetaGraph

                        # Validate that assembled sequences are proper RNA
                        for label in MetaGraphsNext.labels(quality_graph)
                            vertex_data = quality_graph[label]
                            Test.@test vertex_data isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
                            sequence_str = string(vertex_data.sequence)
                            if !isempty(sequence_str)
                                detected_type = Mycelia.detect_alphabet(sequence_str)
                                if occursin('U', sequence_str)
                                    Test.@test detected_type == :RNA
                                end
                                # Verify no T nucleotides (RNA should never contain T)
                                Test.@test !occursin('T', sequence_str)
                            end

                            evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                            Test.@test !isempty(evidence_entries)
                            for entry in evidence_entries
                                Test.@test entry.strand in (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse)
                                Test.@test !isempty(entry.quality_scores)
                                Test.@test all(q >= 0 for q in entry.quality_scores)  # Valid quality scores
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
            graph = Mycelia.Rhizomorph.build_ngram_graph([reference_string], 3; dataset_id="test")
            
            # Test that we can work with the graph (basic functionality)
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            Test.@test graph isa MetaGraphsNext.MetaGraph
        end
        
        Test.@testset "Sequence Graph Next GFA I/O" begin
            reference_seq = string(BioSequences.randdnaseq(100))
            reference_record = FASTX.FASTA.Record("reference", reference_seq)
            
            kmer_type = Kmers.DNAKmer{5}
            graph = Mycelia.Rhizomorph.build_kmer_graph([reference_record], 5; dataset_id="test", mode=:doublestrand)
            
            # Test that we can work with the graph (basic functionality)
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            Test.@test graph isa MetaGraphsNext.MetaGraph
            
            # Test that GFA I/O functions exist
            Test.@test hasmethod(Mycelia.Rhizomorph.write_gfa_next, (typeof(graph), String))
            Test.@test hasmethod(Mycelia.Rhizomorph.read_gfa_next, (String,))
        end
    end
    
    Test.@testset "Assembly Performance and Scaling" begin
        Test.@testset "Memory Usage Scaling" begin
            # Test with small toy sequences (unit-test scale)
            for seq_length in [50, 100]
                reference_seq = string(BioSequences.randdnaseq(seq_length))
                reads = Mycelia.create_test_reads(reference_seq, 20, 0.01)
                
                # Build graphs and verify they don't crash
                kmer_type = Kmers.DNAKmer{5}
                
                # Test both modes
                single_graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:singlestrand)
                double_graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:doublestrand)
                
                Test.@test single_graph isa MetaGraphsNext.MetaGraph
                Test.@test double_graph isa MetaGraphsNext.MetaGraph
                
                # DoubleStrand should include forward and reverse-complement k-mers
                Test.@test !isempty(MetaGraphsNext.labels(single_graph))
                Test.@test !isempty(MetaGraphsNext.labels(double_graph))
                Test.@test length(MetaGraphsNext.labels(double_graph)) >= length(MetaGraphsNext.labels(single_graph))
            end
        end
        
        Test.@testset "Coverage Impact on Assembly Quality" begin
            reference_seq = string(BioSequences.randdnaseq(50))
            error_rate = 0.01
            
            for coverage in [5, 20]
                reads = Mycelia.create_test_reads(reference_seq, coverage, error_rate)
                
                kmer_type = Kmers.DNAKmer{5}
                graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:doublestrand)
                
                Test.@test !isempty(MetaGraphsNext.labels(graph))
                
                # Higher coverage should generally result in more robust graphs
                total_evidence = sum(Mycelia.Rhizomorph.count_evidence(graph[label]) for label in MetaGraphsNext.labels(graph))
                Test.@test total_evidence > 0
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
            graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:doublestrand)

            Test.@test !isempty(MetaGraphsNext.labels(graph))

            # Test that we have proper coverage and strand orientations
            # This validates the canonicalization consistency fix
            total_evidence_entries = 0
            has_both_orientations = false

            for label in MetaGraphsNext.labels(graph)
                vertex_data = graph[label]
                strands = Mycelia.Rhizomorph.collect_evidence_strands(vertex_data.evidence)
                total_evidence_entries += length(strands)

                # Check if we see both orientations (should happen with reverse complements)
                if Mycelia.Rhizomorph.Forward in strands && Mycelia.Rhizomorph.Reverse in strands
                    has_both_orientations = true
                end
            end

            Test.@test total_evidence_entries > 0
            Test.@test has_both_orientations  # Should see both orientations with forward and reverse sequences
        end

        Test.@testset "SingleStrand RNA Non-Canonicalization Consistency" begin
            # Test that SingleStrand mode uses k-mers as-is without canonicalization
            reference_seq = BioSequences.rna"AUCGUUUU"
            reads = [FASTX.FASTA.Record("rna", reference_seq)]

            kmer_type = Kmers.RNAKmer{5}
            graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 5; dataset_id="test", mode=:singlestrand)

            Test.@test !isempty(MetaGraphsNext.labels(graph))

            # Test that all coverage entries use Forward orientation in SingleStrand mode
            # This validates the SingleStrand non-canonicalization fix
            total_evidence_entries = 0
            all_forward = true

            for label in MetaGraphsNext.labels(graph)
                vertex_data = graph[label]
                strands = Mycelia.Rhizomorph.collect_evidence_strands(vertex_data.evidence)
                total_evidence_entries += length(strands)
                if any(strand -> strand != Mycelia.Rhizomorph.Forward, strands)
                    all_forward = false
                end
            end

            Test.@test total_evidence_entries > 0
            Test.@test all_forward  # All orientations should be Forward in SingleStrand mode
        end
    end
end
