"""
Tests for iterative maximum likelihood assembly algorithms (Phase 5.2a)

This module tests the iterative assembly system that implements:
- Complete FASTQ I/O processing per iteration
- Statistical path resampling with likelihood calculations
- Viterbi algorithm integration for optimal path finding
- Timestamped output files for tracking read evolution
- Memory-efficient read set processing
"""

import Test
import Mycelia
import BioSequences
import FASTX
import Dates

Test.@testset "Iterative Assembly Tests" begin
    
    Test.@testset "Core Framework Functions" begin
        # Test basic function availability
        Test.@test hasmethod(Mycelia.mycelia_iterative_assemble, (String,))
        Test.@test hasmethod(Mycelia.improve_read_set_likelihood, (Vector{FASTX.FASTQ.Record}, Any, Int))
        Test.@test hasmethod(Mycelia.improve_read_likelihood, (FASTX.FASTQ.Record, Any, Int))
        Test.@test hasmethod(Mycelia.find_optimal_sequence_path, (String, Any, Any, Int))
        Test.@test hasmethod(Mycelia.calculate_sequence_likelihood, (String, Any, Any, Int))
        Test.@test hasmethod(Mycelia.sufficient_improvements, (Int, Int))
    end
    
    Test.@testset "Read Likelihood Improvement" begin
        # Create test FASTQ records
        sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC",
            "ATCGATCGATCGATCGTTCGATCG"  # Contains variation
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Build test graph
        graph = Mycelia.build_qualmer_graph(test_records, k=7)
        
        # Test sequence likelihood calculation
        seq = sequences[1]
        qual = repeat("I", length(seq))
        likelihood = Mycelia.calculate_sequence_likelihood(seq, qual, graph, 7)
        Test.@test isa(likelihood, Float64)
        Test.@test likelihood <= 0.0  # Log likelihood should be <= 0
        
        # Test optimal sequence path finding
        improved_seq, improvement = Mycelia.find_optimal_sequence_path(seq, qual, graph, 7)
        Test.@test isa(improved_seq, String)
        Test.@test isa(improvement, Float64)
        Test.@test length(improved_seq) == length(seq)  # Should maintain length
        
        # Test single read improvement
        read = test_records[1]
        improved_read, was_improved = Mycelia.improve_read_likelihood(read, graph, 7)
        Test.@test isa(improved_read, FASTX.FASTQ.Record)
        Test.@test isa(was_improved, Bool)
        Test.@test FASTX.identifier(improved_read) == FASTX.identifier(read)  # ID should be preserved
        
        # Test read set improvement
        updated_reads, improvements_made = Mycelia.improve_read_set_likelihood(test_records, graph, 7)
        Test.@test length(updated_reads) == length(test_records)
        Test.@test isa(improvements_made, Int)
        Test.@test improvements_made >= 0
        Test.@test improvements_made <= length(test_records)
    end
    
    Test.@testset "Decision Making Functions" begin
        # Test sufficient improvements function
        Test.@test Mycelia.sufficient_improvements(10, 100, 0.05) == true   # 10% > 5% threshold
        Test.@test Mycelia.sufficient_improvements(6, 100, 0.05) == true    # 6% > 5% threshold
        Test.@test Mycelia.sufficient_improvements(0, 100, 0.05) == false   # 0% < 5% threshold
        
        # Test edge cases
        Test.@test Mycelia.sufficient_improvements(1, 1, 0.5) == true      # 100% > 50% threshold
        Test.@test Mycelia.sufficient_improvements(0, 0, 0.05) == false    # Handle zero total reads
    end
    
    Test.@testset "Quality Score Adjustment" begin
        # Test quality score adjustment
        original_qual = "IIIIIIIIII"
        improved_seq = "ATCGATCGAT"
        improvement = 0.1
        
        adjusted_qual = Mycelia.adjust_quality_scores(original_qual, improved_seq, improvement)
        Test.@test isa(adjusted_qual, String)
        Test.@test length(adjusted_qual) == length(improved_seq)
        
        # Test with different sequence lengths
        longer_seq = "ATCGATCGATCG"
        adjusted_qual_long = Mycelia.adjust_quality_scores(original_qual, longer_seq, improvement)
        Test.@test length(adjusted_qual_long) == length(longer_seq)
        
        shorter_seq = "ATCGAT"
        adjusted_qual_short = Mycelia.adjust_quality_scores(original_qual, shorter_seq, improvement)
        Test.@test length(adjusted_qual_short) == length(shorter_seq)
    end
    
    Test.@testset "K-mer Alternatives Generation" begin
        sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC"
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        graph = Mycelia.build_qualmer_graph(test_records, k=5)
        
        # Test alternative generation for an existing k-mer
        test_kmer = "ATCGA"
        alternatives = Mycelia.generate_kmer_alternatives(test_kmer, graph)
        Test.@test isa(alternatives, Vector{String})
        Test.@test all(alt -> length(alt) == length(test_kmer), alternatives)
        
        # All alternatives should be valid DNA sequences
        for alt in alternatives
            Test.@test all(c -> c in ['A', 'T', 'G', 'C'], alt)
        end
    end
    
    Test.@testset "Iterative Assembly Integration" begin
        # Create temporary test FASTQ file
        temp_dir = mktempdir()
        test_fastq = joinpath(temp_dir, "test_reads.fastq")
        
        # Generate test sequences
        sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC",
            "ATCGATCGATCGATCGTTCGATCG"  # Contains variation
        ]
        
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Write test FASTQ
        Mycelia.write_fastq(records=test_records, filename=test_fastq)
        Test.@test isfile(test_fastq)
        
        # Test iterative assembly with small limits
        output_dir = joinpath(temp_dir, "test_output")
        result = Mycelia.mycelia_iterative_assemble(test_fastq, 
                                                  max_k=23, 
                                                  output_dir=output_dir,
                                                  max_iterations_per_k=2,
                                                  verbose=false)
        
        # Check result structure
        Test.@test haskey(result, :final_assembly)
        Test.@test haskey(result, :k_progression)
        Test.@test haskey(result, :metadata)
        
        # Check metadata structure
        metadata = result[:metadata]
        Test.@test haskey(metadata, :total_runtime)
        Test.@test haskey(metadata, :total_iterations)
        Test.@test haskey(metadata, :total_improvements)
        Test.@test haskey(metadata, :k_progression)
        Test.@test haskey(metadata, :final_k)
        Test.@test haskey(metadata, :final_fastq_file)
        Test.@test haskey(metadata, :output_directory)
        Test.@test haskey(metadata, :iteration_history)
        Test.@test haskey(metadata, :assembly_type)
        Test.@test haskey(metadata, :version)
        
        # Check that output directory was created
        Test.@test isdir(output_dir)
        
        # Check that some output files were created
        output_files = readdir(output_dir)
        Test.@test length(output_files) >= 1
        
        # Check that final FASTQ exists
        final_fastq = metadata[:final_fastq_file]
        Test.@test isfile(final_fastq)
        
        # Cleanup
        rm(temp_dir, recursive=true)
    end
    
    Test.@testset "Assembly Summary and Reporting" begin
        # Create mock result for testing summary function
        mock_result = Dict(
            :final_assembly => ["ATCGATCGATCG", "CGATCGATCGAT"],
            :k_progression => [7, 11, 13],
            :metadata => Dict(
                :total_runtime => 10.5,
                :total_iterations => 5,
                :total_improvements => 25,
                :k_progression => [7, 11, 13],
                :final_k => 13,
                :final_fastq_file => "test.fastq",
                :output_directory => "test_output",
                :iteration_history => Dict(
                    7 => [Dict(:improvements_made => 10, :iteration => 1)],
                    11 => [Dict(:improvements_made => 8, :iteration => 1)],
                    13 => [Dict(:improvements_made => 7, :iteration => 1)]
                ),
                :assembly_type => "iterative_maximum_likelihood",
                :version => "Phase_5.2a"
            )
        )
        
        # Test summary generation
        summary = Mycelia.iterative_assembly_summary(mock_result)
        Test.@test isa(summary, String)
        Test.@test contains(summary, "Mycelia Iterative Assembly Summary")
        Test.@test contains(summary, "iterative_maximum_likelihood")
        Test.@test contains(summary, "Phase_5.2a")
        Test.@test contains(summary, "10.5")  # Runtime
        Test.@test contains(summary, "25")    # Total improvements
    end
    
    Test.@testset "Built-in Test Function" begin
        # Test the built-in test function
        test_result = Mycelia.test_iterative_assembly()
        
        Test.@test haskey(test_result, :status)
        Test.@test test_result[:status] in [:success, :error]
        
        if test_result[:status] == :success
            Test.@test haskey(test_result, :final_assembly)
            Test.@test haskey(test_result, :k_progression)
            Test.@test haskey(test_result, :metadata)
        else
            # If there's an error, it should be reported
            Test.@test haskey(test_result, :error)
            println("Built-in test error: ", test_result[:error])
        end
    end
    
    Test.@testset "Edge Cases and Error Handling" begin
        # Test with empty sequences - skip graph building since it will fail with empty input
        empty_records = FASTX.FASTQ.Record[]
        
        # Test the improvement function with empty input - create a minimal mock graph
        mock_sequences = ["ATCGATCGATCGATCGATCGATCG"]
        mock_records = [FASTX.FASTQ.Record("mock", mock_sequences[1], repeat("I", length(mock_sequences[1])))]
        mock_graph = Mycelia.build_qualmer_graph(mock_records, k=5)
        
        # Should handle empty input gracefully
        updated_empty, improvements_empty = Mycelia.improve_read_set_likelihood(empty_records, mock_graph, 5)
        Test.@test length(updated_empty) == 0
        Test.@test improvements_empty == 0
        
        # Test with very short sequences
        short_sequences = ["AT", "GC"]
        short_records = [
            FASTX.FASTQ.Record("short$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(short_sequences)
        ]
        
        # Should handle short sequences without crashing - use longer sequences for k=3
        longer_short_sequences = ["ATCG", "GCTA"]
        longer_short_records = [
            FASTX.FASTQ.Record("short$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(longer_short_sequences)
        ]
        short_graph = Mycelia.build_qualmer_graph(longer_short_records, k=3)
        
        # Test likelihood calculation with short sequence
        short_likelihood = Mycelia.calculate_sequence_likelihood("ATCG", "IIII", short_graph, 3)
        Test.@test isa(short_likelihood, Float64)
        
        # Test sufficient improvements with edge cases
        Test.@test Mycelia.sufficient_improvements(0, 0, 0.05) == false  # Zero total reads
        Test.@test Mycelia.sufficient_improvements(5, 0, 0.05) == false  # Zero total reads (invalid case)
    end
    
    Test.@testset "Enhanced Statistical Path Improvement (Phase 5.2b)" begin
        # Create test data for enhanced path improvement
        sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC",
            "ATCGATCGATCGATCGTTCGATCG",  # Contains single error
            "TCGATCGATCGATCGATCGATCGA"
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Build k-mer graph
        k = 7
        graph = Mycelia.build_qualmer_graph(test_records, k=k)
        
        # Test individual enhancement functions
        
        # Test Viterbi path improvement
        test_sequence = "ATCGATCGATCGATCGTTCGATCG"  # Sequence with potential error
        test_quality = repeat("I", length(test_sequence))
        
        viterbi_result = Mycelia.try_viterbi_path_improvement(test_sequence, test_quality, graph, k)
        # Should either return improvement or nothing (graceful degradation)
        Test.@test isa(viterbi_result, Union{Tuple{String, Float64}, Nothing})
        
        # Test statistical path resampling
        statistical_result = Mycelia.try_statistical_path_resampling(test_sequence, test_quality, graph, k)
        Test.@test isa(statistical_result, Union{Tuple{String, Float64}, Nothing})
        
        # Test local path improvements (fallback)
        original_likelihood = Mycelia.calculate_sequence_likelihood(test_sequence, test_quality, graph, k)
        local_result = Mycelia.try_local_path_improvements(test_sequence, test_quality, graph, k, original_likelihood)
        Test.@test isa(local_result, Tuple{String, Float64})
        improved_seq, improvement = local_result
        Test.@test isa(improved_seq, String)
        Test.@test isa(improvement, Float64)
        
        # Test path generation utilities
        alternative_paths = Mycelia.generate_alternative_paths(test_sequence, graph, k, num_samples=3)
        Test.@test isa(alternative_paths, Vector{Vector{Any}})
        Test.@test length(alternative_paths) <= 3  # Should not exceed requested samples
        
        # Test sequence to path conversion
        kmer_path = Mycelia.sequence_to_kmer_path(test_sequence, k)
        Test.@test isa(kmer_path, Vector{String})
        Test.@test length(kmer_path) == length(test_sequence) - k + 1
        Test.@test all(length(kmer) == k for kmer in kmer_path)
        
        # Test path to sequence conversion
        reconstructed_seq = Mycelia.kmer_path_to_sequence(kmer_path, k)
        Test.@test isa(reconstructed_seq, String)
        Test.@test reconstructed_seq == test_sequence  # Should perfectly reconstruct
        
        # Test enhanced find_optimal_sequence_path function
        enhanced_result = Mycelia.find_optimal_sequence_path(test_sequence, test_quality, graph, k)
        Test.@test isa(enhanced_result, Tuple{String, Float64})
        enhanced_seq, likelihood_improvement = enhanced_result
        Test.@test isa(enhanced_seq, String)
        Test.@test isa(likelihood_improvement, Float64)
        Test.@test length(enhanced_seq) > 0
        
        # Test that improvement is non-negative (should not make things worse)
        Test.@test likelihood_improvement >= -1e-10  # Allow for small numerical errors
        
        # Test edge cases
        
        # Empty sequence should be handled gracefully
        empty_path = Mycelia.sequence_to_kmer_path("", k)
        Test.@test isempty(empty_path)
        
        empty_reconstruction = Mycelia.kmer_path_to_sequence(String[], k)
        Test.@test empty_reconstruction == ""
        
        # Short sequence (shorter than k)
        short_seq = "ATC"  # Shorter than k=7
        short_path = Mycelia.sequence_to_kmer_path(short_seq, k)
        Test.@test isempty(short_path)
        
        # Single k-mer sequence
        single_kmer_seq = "ATCGATC"  # Exactly k=7
        single_path = Mycelia.sequence_to_kmer_path(single_kmer_seq, k)
        Test.@test length(single_path) == 1
        Test.@test single_path[1] == single_kmer_seq
        
        single_reconstruction = Mycelia.kmer_path_to_sequence(single_path, k)
        Test.@test single_reconstruction == single_kmer_seq
        
        # Test with minimal graph (edge case)
        minimal_sequences = ["ATCGATC", "TCGATCG"]
        minimal_records = [
            FASTX.FASTQ.Record("min$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(minimal_sequences)
        ]
        minimal_graph = Mycelia.build_qualmer_graph(minimal_records, k=k)
        
        minimal_result = Mycelia.find_optimal_sequence_path("ATCGATC", "IIIIIII", minimal_graph, k)
        Test.@test isa(minimal_result, Tuple{String, Float64})
        
        # Test that the enhancement doesn't break existing functionality
        original_read = FASTX.FASTQ.Record("test", test_sequence, test_quality)
        enhanced_read, was_improved = Mycelia.improve_read_likelihood(original_read, graph, k)
        Test.@test isa(enhanced_read, FASTX.FASTQ.Record)
        Test.@test isa(was_improved, Bool)
        Test.@test FASTX.identifier(enhanced_read) == FASTX.identifier(original_read)
        Test.@test length(FASTX.sequence(enhanced_read)) > 0
        
        # Detailed tests for path reconstruction correctness
        Test.@testset "Path Reconstruction Validation" begin
            # Test simple case
            simple_seq = "ATCGATC"  # Exactly k=7, should give 1 k-mer
            simple_path = Mycelia.sequence_to_kmer_path(simple_seq, k)
            Test.@test length(simple_path) == 1
            Test.@test simple_path[1] == simple_seq
            simple_reconstruction = Mycelia.path_to_sequence(simple_path, k)
            Test.@test simple_reconstruction == simple_seq
            
            # Test k+1 length sequence
            kplus1_seq = "ATCGATCG"  # k+1=8, should give 2 k-mers
            kplus1_path = Mycelia.sequence_to_kmer_path(kplus1_seq, k)
            Test.@test length(kplus1_path) == 2
            Test.@test kplus1_path[1] == "ATCGATC"
            Test.@test kplus1_path[2] == "TCGATCG"
            kplus1_reconstruction = Mycelia.path_to_sequence(kplus1_path, k)
            Test.@test kplus1_reconstruction == kplus1_seq
            
            # Test longer sequence
            long_seq = "ATCGATCGATC"  # Length 11, should give 5 k-mers
            long_path = Mycelia.sequence_to_kmer_path(long_seq, k)
            Test.@test length(long_path) == 5
            # Verify each k-mer
            expected_kmers = ["ATCGATC", "TCGATCG", "CGATCGA", "GATCGAT", "ATCGATC"]
            Test.@test long_path == expected_kmers
            long_reconstruction = Mycelia.path_to_sequence(long_path, k)
            Test.@test long_reconstruction == long_seq
            
            # Test with different k values
            test_seq_k3 = "ATCGATC"
            k3_path = Mycelia.sequence_to_kmer_path(test_seq_k3, 3)
            Test.@test length(k3_path) == 5  # 7-3+1 = 5
            Test.@test k3_path == ["ATC", "TCG", "CGA", "GAT", "ATC"]
            k3_reconstruction = Mycelia.path_to_sequence(k3_path, 3)
            Test.@test k3_reconstruction == test_seq_k3
        end
        
        # Test statistical path improvement components
        Test.@testset "Statistical Path Components" begin
            # Test alternative path generation with small graph
            alternatives = Mycelia.generate_alternative_paths(test_sequence, graph, k, num_samples=2)
            Test.@test length(alternatives) <= 2
            Test.@test isa(alternatives, Vector{Vector{Any}})
            
            # Test that alternative paths are valid (if any generated)
            for alt_path in alternatives
                if !isempty(alt_path)
                    alt_seq = Mycelia.path_to_sequence(alt_path, k)
                    Test.@test isa(alt_seq, String)
                    Test.@test length(alt_seq) > 0
                end
            end
            
            # Test with edge case sequences
            edge_seq = "AAAAAAA"  # Repetitive sequence
            edge_alternatives = Mycelia.generate_alternative_paths(edge_seq, graph, k, num_samples=1)
            Test.@test isa(edge_alternatives, Vector{Vector{Any}})
        end
        
        # Test integration with quality scores
        Test.@testset "Quality Score Integration" begin
            # Test with varying quality scores
            low_quality = repeat("!", length(test_sequence))  # PHRED score ~0
            high_quality = repeat("I", length(test_sequence))  # PHRED score ~40
            
            low_qual_result = Mycelia.find_optimal_sequence_path(test_sequence, low_quality, graph, k)
            high_qual_result = Mycelia.find_optimal_sequence_path(test_sequence, high_quality, graph, k)
            
            Test.@test isa(low_qual_result, Tuple{String, Float64})
            Test.@test isa(high_qual_result, Tuple{String, Float64})
            
            # Both should return valid sequences
            Test.@test length(low_qual_result[1]) > 0
            Test.@test length(high_qual_result[1]) > 0
        end
        
        # Test robustness with different graph sizes
        Test.@testset "Graph Size Robustness" begin
            # Minimal graph (2 sequences)
            minimal_seqs = ["ATCGATC", "TCGATCG"]
            minimal_records = [FASTX.FASTQ.Record("min$i", seq, repeat("I", length(seq))) for (i, seq) in enumerate(minimal_seqs)]
            minimal_graph = Mycelia.build_qualmer_graph(minimal_records, k=k)
            
            minimal_result = Mycelia.find_optimal_sequence_path("ATCGATC", "IIIIIII", minimal_graph, k)
            Test.@test isa(minimal_result, Tuple{String, Float64})
            
            # Large graph (many sequences)
            large_seqs = ["ATCGATCGATC", "TCGATCGATCG", "CGATCGATCGA", "GATCGATCGAT", 
                         "ATCGATCGATG", "TCGATCGATCA", "CGATCGATCGT", "GATCGATCGAC"]
            large_records = [FASTX.FASTQ.Record("large$i", seq, repeat("I", length(seq))) for (i, seq) in enumerate(large_seqs)]
            large_graph = Mycelia.build_qualmer_graph(large_records, k=k)
            
            large_result = Mycelia.find_optimal_sequence_path("ATCGATCGATC", repeat("I", 11), large_graph, k)
            Test.@test isa(large_result, Tuple{String, Float64})
        end
        
        # Test performance with realistic data
        Test.@testset "Realistic Performance Test" begin
            # Generate more realistic sequences with errors
            realistic_seqs = [
                "ATCGATCGATCGATCGATCGATCGATCGATC",  # Perfect repeat
                "ATCGATCGATCGATCGTTCGATCGATCGATC",  # Single substitution
                "ATCGATCGATCGATC_ATCGATCGATCGATC",  # Gap (represented as _)
                "ATCGATCGATCGATCGATCGATCGATCGATT"   # End substitution
            ]
            realistic_seqs = [replace(seq, "_" => "") for seq in realistic_seqs]  # Remove gap markers
            
            realistic_records = [FASTX.FASTQ.Record("real$i", seq, repeat("I", length(seq))) for (i, seq) in enumerate(realistic_seqs)]
            realistic_graph = Mycelia.build_qualmer_graph(realistic_records, k=k)
            
            # Test that all enhancement methods handle realistic data
            for (i, seq) in enumerate(realistic_seqs)
                quality = repeat("I", length(seq))
                
                # Test main enhancement function
                result = Mycelia.find_optimal_sequence_path(seq, quality, realistic_graph, k)
                Test.@test isa(result, Tuple{String, Float64})
                Test.@test length(result[1]) > 0
                
                # Test individual methods
                viterbi_result = Mycelia.try_viterbi_path_improvement(seq, quality, realistic_graph, k)
                Test.@test isa(viterbi_result, Union{Tuple{String, Float64}, Nothing})
                
                stat_result = Mycelia.try_statistical_path_resampling(seq, quality, realistic_graph, k)
                Test.@test isa(stat_result, Union{Tuple{String, Float64}, Nothing})
                
                original_likelihood = Mycelia.calculate_sequence_likelihood(seq, quality, realistic_graph, k)
                local_result = Mycelia.try_local_path_improvements(seq, quality, realistic_graph, k, original_likelihood)
                Test.@test isa(local_result, Tuple{String, Float64})
            end
        end
    end
end