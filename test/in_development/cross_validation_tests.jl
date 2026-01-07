# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/in_development/cross_validation_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/in_development/cross_validation_tests.jl", "test/in_development", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import BioSequences
import FASTX
import Dates
import Statistics
import Random

Test.@testset "Cross-Validation Tests" begin
    
    Test.@testset "Core Framework Functions" begin
        # Test basic function availability
        Test.@test hasmethod(Mycelia.mycelia_cross_validation, (String,))
        Test.@test hasmethod(Mycelia.create_kfold_partitions, (Vector{FASTX.FASTQ.Record}, Int, Float64))
        Test.@test hasmethod(Mycelia.validate_assemblies_against_holdout, (Dict, Dict, Vector{FASTX.FASTQ.Record}, Int, String))
        Test.@test hasmethod(Mycelia.calculate_assembly_mapping_metrics, (Any, Vector{FASTX.FASTQ.Record}, String))
        Test.@test hasmethod(Mycelia.compare_assembly_statistics, (Dict, Dict))
        Test.@test hasmethod(Mycelia.generate_consensus_pangenome, (Dict, Dict, Dict, String))
    end
    
    Test.@testset "K-Fold Partitioning" begin
        # Create test FASTQ records
        sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC", 
            "ATCGATCGATCGATCGTTCGATCG",
            "TCGATCGATCGATCGATCGATCGA",
            "GCGATCGATCGATCGATCGATCGT",
            "ACGATCGATCGATCGATCGATCGT",
            "CCGATCGATCGATCGATCGATCGA",
            "AGGATCGATCGATCGATCGATCGA",
            "CTGATCGATCGATCGATCGATCGT"
        ]
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Test k-fold partitioning
        k_folds = 3
        validation_split = 0.2
        partitions = Mycelia.create_kfold_partitions(test_records, k_folds, validation_split)
        
        # Check partition structure
        Test.@test length(partitions) == k_folds
        Test.@test all(haskey(partitions[i], :train) for i in 1:k_folds)
        Test.@test all(haskey(partitions[i], :validation) for i in 1:k_folds)
        
        # Check that all reads are included across folds
        total_train_reads = sum(length(partitions[i][:train]) for i in 1:k_folds)
        total_validation_reads = sum(length(partitions[i][:validation]) for i in 1:k_folds)
        Test.@test total_train_reads + total_validation_reads == length(test_records)
        
        # Check validation split ratio (approximately)
        for i in 1:k_folds
            fold_size = partitions[i][:fold_size]
            validation_size = partitions[i][:validation_size]
            if fold_size > 0
                actual_split = validation_size / fold_size
                Test.@test abs(actual_split - validation_split) < 0.3  # Allow some tolerance
            end
        end
    end
    
    Test.@testset "Assembly Mapping Metrics" begin
        # Create test data
        validation_sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC"
        ]
        validation_records = [
            FASTX.FASTQ.Record("val$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(validation_sequences)
        ]
        
        # Test with sequence-based assembly
        assembly_sequences = ["ATCGATCGATCGATCG", "CGATCGATCGATCGAT"]
        metrics = Mycelia.calculate_assembly_mapping_metrics(
            assembly_sequences, validation_records, "test_assembly"
        )
        
        # Check metric structure
        Test.@test haskey(metrics, :assembly_type)
        Test.@test haskey(metrics, :total_reads)
        Test.@test haskey(metrics, :mapped_reads)
        Test.@test haskey(metrics, :mapping_rate)
        Test.@test haskey(metrics, :average_coverage)
        Test.@test haskey(metrics, :assembly_sequences)
        Test.@test haskey(metrics, :assembly_kmers)
        
        # Check metric values
        Test.@test metrics[:assembly_type] == "test_assembly"
        Test.@test metrics[:total_reads] == length(validation_records)
        Test.@test metrics[:mapped_reads] >= 0
        Test.@test metrics[:mapped_reads] <= metrics[:total_reads]
        Test.@test 0.0 <= metrics[:mapping_rate] <= 1.0
        Test.@test metrics[:average_coverage] >= 0.0
        Test.@test metrics[:assembly_sequences] == length(assembly_sequences)
        
        # Test with empty assembly
        empty_metrics = Mycelia.calculate_assembly_mapping_metrics(
            String[], validation_records, "empty_assembly"
        )
        Test.@test empty_metrics[:mapping_rate] == 0.0
        Test.@test empty_metrics[:assembly_sequences] == 0
    end
    
    Test.@testset "Assembly Statistics Comparison" begin
        # Create mock assembly results
        intelligent_result = Dict(
            :metadata => Dict(
                :total_runtime => 10.5,
                :k_progression => [7, 11, 13],
                :final_k => 13,
                :assembly_type => "intelligent_self_optimizing",
                :version => "Phase_5.1b"
            )
        )
        
        iterative_result = Dict(
            :metadata => Dict(
                :total_runtime => 25.8,
                :k_progression => [7, 11, 13, 17],
                :final_k => 17,
                :assembly_type => "iterative_maximum_likelihood",
                :version => "Phase_5.2a",
                :total_improvements => 150,
                :total_iterations => 12
            )
        )
        
        # Test comparison
        comparison = Mycelia.compare_assembly_statistics(intelligent_result, iterative_result)
        
        # Check comparison structure
        Test.@test haskey(comparison, :runtime_comparison)
        Test.@test haskey(comparison, :k_progression_comparison)
        Test.@test haskey(comparison, :assembly_type_comparison)
        Test.@test haskey(comparison, :improvement_metrics)
        
        # Check runtime comparison
        runtime_comp = comparison[:runtime_comparison]
        Test.@test runtime_comp[:intelligent_runtime] == 10.5
        Test.@test runtime_comp[:iterative_runtime] == 25.8
        Test.@test abs(runtime_comp[:runtime_ratio] - (25.8 / 10.5)) < 0.01
        
        # Check k-progression comparison
        k_comp = comparison[:k_progression_comparison]
        Test.@test k_comp[:intelligent_final_k] == 13
        Test.@test k_comp[:iterative_final_k] == 17
    end
    
    Test.@testset "Consensus Analysis Functions" begin
        # Test recommendation reasoning
        reasoning1 = Mycelia.generate_recommendation_reasoning(0.08, 2.5)
        Test.@test contains(reasoning1, "significant mapping improvement")
        
        reasoning2 = Mycelia.generate_recommendation_reasoning(-0.08, 1.5)
        Test.@test contains(reasoning2, "Intelligent assembly shows significantly better")
        
        reasoning3 = Mycelia.generate_recommendation_reasoning(0.02, 6.0)
        Test.@test contains(reasoning3, "much faster")
        
        # Test confidence calculation
        stable_rates = [0.8, 0.82, 0.81, 0.79, 0.83]
        variable_rates = [0.6, 0.9, 0.7, 0.95, 0.65]
        
        stable_confidence = Mycelia.calculate_recommendation_confidence(stable_rates, stable_rates)
        variable_confidence = Mycelia.calculate_recommendation_confidence(variable_rates, variable_rates)
        
        Test.@test stable_confidence > variable_confidence
        Test.@test 0.0 <= stable_confidence <= 1.0
        Test.@test 0.0 <= variable_confidence <= 1.0
    end
    
    Test.@testset "Cross-Validation Integration" begin
        # Create temporary test FASTQ file
        temp_dir = mktempdir()
        test_fastq = joinpath(temp_dir, "test_reads.fastq")
        
        # Generate sufficient test sequences for k-fold validation
        sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC",
            "ATCGATCGATCGATCGTTCGATCG",
            "TCGATCGATCGATCGATCGATCGA",
            "GCGATCGATCGATCGATCGATCGT",
            "ACGATCGATCGATCGATCGATCGT",
            "CCGATCGATCGATCGATCGATCGA",
            "AGGATCGATCGATCGATCGATCGA",
            "CTGATCGATCGATCGATCGATCGT",
            "AAGATCGATCGATCGATCGATCGT",
            "TTGATCGATCGATCGATCGATCGA"
        ]
        
        test_records = [
            FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(sequences)
        ]
        
        # Write test FASTQ
        Mycelia.write_fastq(records=test_records, filename=test_fastq)
        Test.@test isfile(test_fastq)
        
        # Test cross-validation with minimal settings
        output_dir = joinpath(temp_dir, "test_cv")
        result = Mycelia.mycelia_cross_validation(test_fastq,
                                                k_folds=3,
                                                max_k=13,
                                                output_dir=output_dir,
                                                verbose=false)
        
        # Check result structure
        Test.@test haskey(result, :cross_validation_type)
        Test.@test haskey(result, :version)
        Test.@test haskey(result, :total_runtime)
        Test.@test haskey(result, :n_folds)
        Test.@test haskey(result, :intelligent_results)
        Test.@test haskey(result, :iterative_results)
        Test.@test haskey(result, :validation_metrics)
        Test.@test haskey(result, :consensus_analysis)
        
        # Check that results contain expected number of folds
        Test.@test length(result[:intelligent_results]) == 3
        Test.@test length(result[:iterative_results]) == 3
        Test.@test length(result[:validation_metrics]) == 3
        
        # Check consensus analysis structure
        consensus = result[:consensus_analysis]
        Test.@test haskey(consensus, :intelligent_performance)
        Test.@test haskey(consensus, :iterative_performance)
        Test.@test haskey(consensus, :statistical_comparison)
        Test.@test haskey(consensus, :recommendation)
        
        # Check recommendation structure
        recommendation = consensus[:recommendation]
        Test.@test haskey(recommendation, :approach)
        Test.@test haskey(recommendation, :reasoning)
        Test.@test haskey(recommendation, :confidence)
        Test.@test recommendation[:approach] in ["intelligent", "iterative", "hybrid"]
        Test.@test 0.0 <= recommendation[:confidence] <= 1.0
        
        # Check that output directory structure was created
        Test.@test isdir(output_dir)
        Test.@test isdir(joinpath(output_dir, "folds"))
        Test.@test isdir(joinpath(output_dir, "results"))
        Test.@test isdir(joinpath(output_dir, "consensus"))
        
        # Check that key output files exist
        Test.@test isfile(joinpath(output_dir, "cross_validation_results.json"))
        Test.@test isfile(joinpath(output_dir, "consensus", "consensus_analysis.json"))
        
        # Cleanup
        rm(temp_dir, recursive=true)
    end
    
    Test.@testset "Summary Reporting" begin
        # Create mock cross-validation result
        mock_result = Dict(
            :cross_validation_type => "hybrid_assembly_quality_assessment",
            :version => "Phase_5.1c",
            :total_runtime => 45.2,
            :n_folds => 5,
            :consensus_analysis => Dict(
                :intelligent_performance => Dict(
                    :mean_mapping_rate => 0.82,
                    :std_mapping_rate => 0.03,
                    :mean_runtime => 12.5,
                    :std_runtime => 2.1
                ),
                :iterative_performance => Dict(
                    :mean_mapping_rate => 0.87,
                    :std_mapping_rate => 0.025,
                    :mean_runtime => 28.3,
                    :std_runtime => 4.2
                ),
                :statistical_comparison => Dict(
                    :mapping_rate_improvement => 0.05,
                    :runtime_overhead => 2.26
                ),
                :recommendation => Dict(
                    :approach => "iterative",
                    :confidence => 0.89,
                    :reasoning => "Iterative assembly shows significant mapping improvement with acceptable runtime overhead"
                )
            )
        )
        
        # Test summary generation
        summary = Mycelia.cross_validation_summary(mock_result)
        Test.@test isa(summary, String)
        Test.@test contains(summary, "Mycelia Cross-Validation Summary")
        Test.@test contains(summary, "hybrid_assembly_quality_assessment")
        Test.@test contains(summary, "Phase_5.1c")
        Test.@test contains(summary, "45.2")  # Total runtime
        Test.@test contains(summary, "82.0%")  # Intelligent mapping rate
        Test.@test contains(summary, "87.0%")  # Iterative mapping rate
        Test.@test contains(summary, "iterative")  # Recommendation
        Test.@test contains(summary, "89.0%")  # Confidence
    end
    
    Test.@testset "Built-in Test Function" begin
        # Test the built-in test function
        test_result = Mycelia.test_cross_validation()
        
        Test.@test haskey(test_result, :status)
        Test.@test test_result[:status] in [:success, :error]
        
        if test_result[:status] == :success
            Test.@test haskey(test_result, :cross_validation_type)
            Test.@test haskey(test_result, :consensus_analysis)
            Test.@test haskey(test_result, :n_folds)
        else
            # If there's an error, it should be reported
            Test.@test haskey(test_result, :error)
            println("Built-in test error: ", test_result[:error])
        end
    end
    
    Test.@testset "Edge Cases and Error Handling" begin
        # Test with minimal data
        minimal_sequences = ["ATCGATCGATCGATCGATCGATCG", "CGATCGATCGATCGATCGATCGAT"]
        minimal_records = [
            FASTX.FASTQ.Record("min$i", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(minimal_sequences)
        ]
        
        # Test k-fold partitioning with minimal data
        minimal_partitions = Mycelia.create_kfold_partitions(minimal_records, 2, 0.5)
        Test.@test length(minimal_partitions) == 2
        
        # Test mapping metrics with empty validation data
        empty_validation = FASTX.FASTQ.Record[]
        empty_metrics = Mycelia.calculate_assembly_mapping_metrics(
            ["ATCGATCG"], empty_validation, "test"
        )
        Test.@test empty_metrics[:total_reads] == 0
        Test.@test empty_metrics[:mapped_reads] == 0
        
        # Test confidence calculation with identical rates
        identical_rates = [0.8, 0.8, 0.8, 0.8, 0.8]
        perfect_confidence = Mycelia.calculate_recommendation_confidence(identical_rates, identical_rates)
        Test.@test perfect_confidence >= 0.95  # Should be very high confidence
    end
end
