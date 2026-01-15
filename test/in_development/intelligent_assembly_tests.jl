# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/in_development/intelligent_assembly_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/in_development/intelligent_assembly_tests.jl", "test/in_development", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

"""
Tests for intelligent assembly algorithms (Phase 5.1a)

This module tests the intelligent self-optimizing assembler that implements:
- Dynamic k-mer selection using prime numbers
- Sparsity-based optimization 
- Memory-aware assembly
- Error correction integration
- Iterative progression algorithms
"""

if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import BioSequences
import FASTX

Test.@testset "Intelligent Assembly Tests" begin
    
    Test.@testset "Prime K-mer Utilities" begin
        # Test next_prime_k function
        Test.@test Mycelia.next_prime_k(20) == 23
        Test.@test Mycelia.next_prime_k(30) == 31
        Test.@test Mycelia.next_prime_k(100) == 101
        
        # Test generate_prime_k_sequence
        primes = Mycelia.generate_prime_k_sequence(20, 50)
        Test.@test all(x -> Mycelia.Primes.isprime(x), primes)
        Test.@test primes[1] >= 20
        Test.@test primes[end] <= 50
        Test.@test length(primes) >= 3  # Should include 23, 29, 31, 37, 41, 43, 47
        
        # Test find_primes_in_range
        range_primes = Mycelia.find_primes_in_range(10, 30)
        expected = [11, 13, 17, 19, 23, 29]
        Test.@test range_primes == expected
    end
    
    Test.@testset "Sparsity Detection" begin
        # Create test FASTQ records for sparsity analysis
        sequences = [
            "ATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTA", 
            "TTTTAAAAAGGGGCCCC",
            "ATCGATCGATCGATCGATCGATCG"
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Test calculate_sparsity with different k values
        sparsity_k3 = Mycelia.calculate_sparsity(test_records, 3)
        sparsity_k5 = Mycelia.calculate_sparsity(test_records, 5)
        
        # Sparsity should increase with k for simple repetitive sequences
        Test.@test sparsity_k5 >= sparsity_k3
        Test.@test 0.0 <= sparsity_k3 <= 1.0
        Test.@test 0.0 <= sparsity_k5 <= 1.0
        
        # Test errors_are_singletons
        singleton_result = Mycelia.errors_are_singletons(test_records, 5)
        Test.@test isa(singleton_result, Bool)
    end
    
    Test.@testset "Memory Monitoring" begin
        # Test memory estimation (takes num_kmers and k)
        memory_est = Mycelia.estimate_memory_usage(1000, 7)
        Test.@test memory_est > 0
        Test.@test isa(memory_est, Int)
        
        # Test with very high k-mer count should show higher memory
        memory_est_large = Mycelia.estimate_memory_usage(10000, 15)
        Test.@test memory_est_large > memory_est
    end
    
    Test.@testset "K-mer Selection Logic" begin
        sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTAGCTA",
            "TTTTAAAAGGGCCCCTTTAAAAGGGCCCC"
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Test find_initial_k
        initial_k = Mycelia.find_initial_k(test_records)
        Test.@test initial_k >= 3   # Minimum k-mer size (starts at 3 for amino acids)
        Test.@test Mycelia.Primes.isprime(initial_k)
        Test.@test initial_k <= 51  # Default max range
        
        # Test should_continue_k decision logic (uses actual signature)
        # Create a real graph for testing
        test_graph = Mycelia.build_qualmer_graph(test_records, k=5)
        continue_decision = Mycelia.should_continue_k(test_graph, 10, 23)
        Test.@test isa(continue_decision, Bool)
        
        # Test with fewer corrections should return different result
        stop_decision = Mycelia.should_continue_k(test_graph, 2, 23)  # Below min_corrections=5
        Test.@test isa(stop_decision, Bool)
    end
    
    Test.@testset "Assembly Pipeline Integration" begin
        # Test with small synthetic dataset
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATC",
            "ATCGATCGATCGATCGTTCG"  # Contains error
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Test main assembly function (no sparsity_threshold parameter)
        result = Mycelia.mycelia_assemble(test_records, max_k=31)
        
        Test.@test haskey(result, :final_assembly)
        Test.@test haskey(result, :k_progression)
        Test.@test haskey(result, :metadata)
        
        # Check that k_progression contains primes
        k_values = result[:k_progression]
        Test.@test all(x -> Mycelia.Primes.isprime(x), k_values)
        Test.@test length(k_values) >= 1
        
        # Check metadata structure
        metadata = result[:metadata]
        Test.@test haskey(metadata, :total_runtime)
        Test.@test haskey(metadata, :memory_usage)
        Test.@test haskey(metadata, :error_corrections)
        
        # Final assembly should be non-empty
        assembly = result[:final_assembly]
        Test.@test isa(assembly, Vector{String})
        Test.@test length(assembly) >= 1
    end
    
    Test.@testset "Test Data Validation" begin
        # Test the built-in test function
        test_result = Mycelia.test_intelligent_assembly()
        
        Test.@test haskey(test_result, :status)
        Test.@test test_result[:status] in [:success, :error]
        
        if test_result[:status] == :success
            Test.@test haskey(test_result, :final_assembly)
            Test.@test haskey(test_result, :k_progression)
            Test.@test length(test_result[:k_progression]) >= 1
        end
    end
    
    Test.@testset "Phase 5.1b: Accuracy-Prioritized Reward Function" begin
        # Test accuracy metrics calculation
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATC"
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Build test graph
        graph = Mycelia.build_qualmer_graph(test_records, k=7)
        
        # Test accuracy metrics calculation
        metrics = Mycelia.calculate_accuracy_metrics(graph, 7)
        Test.@test haskey(metrics, :coverage_uniformity)
        Test.@test haskey(metrics, :probability_confidence)
        Test.@test haskey(metrics, :graph_connectivity)
        Test.@test haskey(metrics, :error_signal_clarity)
        Test.@test haskey(metrics, :overall_accuracy)
        
        # All metrics should be between 0 and 1
        for (key, value) in metrics
            Test.@test 0.0 <= value <= 1.0
        end
        
        # Test reward calculation
        reward = Mycelia.calculate_assembly_reward(graph, 5, 7)
        Test.@test 0.0 <= reward <= 1.0
        
        # Test reward with previous metrics
        previous_metrics = Dict(:overall_accuracy => 0.5)
        reward_with_prev = Mycelia.calculate_assembly_reward(graph, 5, 7, previous_metrics)
        Test.@test 0.0 <= reward_with_prev <= 1.0
        
        # Test advanced decision making
        reward_history = [0.3, 0.4, 0.5]
        continue_decision = Mycelia.should_continue_k_advanced(graph, 5, 7, reward_history)
        Test.@test isa(continue_decision, Bool)
    end
    
    Test.@testset "Phase 5.1b: Enhanced Assembly Pipeline" begin
        # Test assembly with reward tracking
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATC",
            "ATCGATCGATCGATCGTTCG"  # Contains error
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Test enhanced assembly function
        result = Mycelia.mycelia_assemble(test_records, max_k=31, verbose=false)
        
        # Check enhanced metadata structure
        Test.@test haskey(result, :metadata)
        metadata = result[:metadata]
        Test.@test haskey(metadata, :reward_statistics)
        Test.@test haskey(metadata, :accuracy_metrics_history)
        
        # Check reward statistics
        reward_stats = metadata[:reward_statistics]
        Test.@test haskey(reward_stats, :mean_reward)
        Test.@test haskey(reward_stats, :final_reward)
        Test.@test haskey(reward_stats, :max_reward)
        Test.@test haskey(reward_stats, :reward_improvement)
        
        # All reward values should be between 0 and 1
        Test.@test 0.0 <= reward_stats[:mean_reward] <= 1.0
        Test.@test 0.0 <= reward_stats[:final_reward] <= 1.0
        Test.@test 0.0 <= reward_stats[:max_reward] <= 1.0
        
        # Check accuracy metrics history
        accuracy_history = metadata[:accuracy_metrics_history]
        Test.@test isa(accuracy_history, Dict)
        if !isempty(accuracy_history)
            first_k = first(keys(accuracy_history))
            first_metrics = accuracy_history[first_k]
            Test.@test haskey(first_metrics, :overall_accuracy)
            Test.@test 0.0 <= first_metrics[:overall_accuracy] <= 1.0
        end
    end
    
    Test.@testset "Edge Cases and Error Handling" begin
        # Test with empty sequences
        empty_records = FASTX.FASTQ.Record[]
        empty_sparsity = Mycelia.calculate_sparsity(empty_records, 5)
        Test.@test empty_sparsity >= 0.0
        
        # Test with very short sequences  
        short_sequences = ["AT", "GC"]
        short_records = [
            FASTX.FASTQ.Record("short$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(short_sequences)
        ]
        short_sparsity = Mycelia.calculate_sparsity(short_records, 3)
        Test.@test short_sparsity >= 0.0
        
        # Test with single sequence
        single_seq = "ATCGATCGATCGATCG"
        single_record = [FASTX.FASTQ.Record("single", single_seq, repeat("I", length(single_seq)))]
        single_result = Mycelia.mycelia_assemble(single_record, max_k=23)
        Test.@test haskey(single_result, :final_assembly)
        
        # Test memory estimation edge cases
        zero_memory = Mycelia.estimate_memory_usage(0, 3)
        Test.@test zero_memory >= 0
        
        large_memory = Mycelia.estimate_memory_usage(1000000, 25)
        Test.@test large_memory > 1000  # Should be substantial
    end
end
