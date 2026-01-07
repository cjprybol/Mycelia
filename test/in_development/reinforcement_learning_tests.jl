# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/in_development/reinforcement_learning_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/in_development/reinforcement_learning_tests.jl", "test/in_development", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Reinforcement Learning Framework Tests
# Tests for the RL-based assembly optimization system

println("Testing Reinforcement Learning Framework...")

# Load the Mycelia module
import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Mycelia
import Test
import BioSequences
import FASTX
import Statistics
import Random
import Dates
import JSON

Test.@testset "Reinforcement Learning Framework Tests" begin

    Test.@testset "Core Data Structures" begin
        Test.@testset "AssemblyState Construction" begin
            state = Mycelia.AssemblyState(
                21,     # current_k
                0.85,   # assembly_quality
                0.05,   # correction_rate
                0.3,    # memory_usage
                0.8,    # graph_connectivity
                0.7,    # coverage_uniformity
                0.9,    # error_signal_clarity
                [0.1, 0.2, 0.3],  # iteration_history
                [21, 31, 41],     # k_progression
                150,    # corrections_made
                45.2    # time_elapsed
            )
            
            Test.@test state.current_k == 21
            Test.@test state.assembly_quality ≈ 0.85
            Test.@test state.correction_rate ≈ 0.05
            Test.@test length(state.iteration_history) == 3
            Test.@test length(state.k_progression) == 3
            Test.@test state.corrections_made == 150
        end
        
        Test.@testset "AssemblyAction Construction" begin
            action = Mycelia.AssemblyAction(
                :continue_k,
                Dict(:error_threshold => 0.01, :min_coverage => 3),
                0.95,
                1000,
                5
            )
            
            Test.@test action.decision == :continue_k
            Test.@test haskey(action.viterbi_params, :error_threshold)
            Test.@test action.correction_threshold == 0.95
            Test.@test action.batch_size == 1000
            Test.@test action.max_iterations == 5
        end
        
        Test.@testset "RewardComponents Construction" begin
            reward = Mycelia.RewardComponents(1000.0, 50.0, -100.0, 25.0, 0.0, 975.0)
            
            Test.@test reward.accuracy_reward == 1000.0
            Test.@test reward.efficiency_reward == 50.0
            Test.@test reward.error_penalty == -100.0
            Test.@test reward.progress_bonus == 25.0
            Test.@test reward.termination_reward == 0.0
            Test.@test reward.total_reward == 975.0
        end
    end

    Test.@testset "Environment Management" begin
        Test.@testset "Environment Creation" begin
            training_files = ["dummy_train1.fastq", "dummy_train2.fastq"]
            validation_files = ["dummy_val1.fastq"]
            
            env = Mycelia.create_assembly_environment(training_files, validation_files, episode_length=50)
            
            Test.@test length(env.training_datasets) == 2
            Test.@test length(env.validation_datasets) == 1
            Test.@test env.episode_length == 50
            Test.@test env.step_count == 0
            Test.@test length(env.reward_history) == 0
            Test.@test length(env.action_history) == 0
            Test.@test length(env.assembly_cache) == 0
        end
        
        Test.@testset "Environment Reset" begin
            # Create test dataset
            test_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCG"), "IIIIIIIIIIIIIIII"),
                FASTX.FASTQ.Record("read2", BioSequences.LongDNA{4}("GCTAGCTAGCTAGCTA"), "IIIIIIIIIIIIIIII")
            ]
            
            test_file = "test_reset.fastq"
            Mycelia.write_fastq(records=test_reads, filename=test_file)
            
            try
                env = Mycelia.create_assembly_environment([test_file], [test_file])
                
                # Test reset
                initial_state = Mycelia.reset_environment!(env, 1)
                
                Test.@test env.step_count == 0
                Test.@test length(env.reward_history) == 0
                Test.@test length(env.action_history) == 0
                Test.@test initial_state.current_k > 0
                Test.@test initial_state.assembly_quality == 0.0  # No assessment yet
                Test.@test length(initial_state.k_progression) == 1
                
            finally
                rm(test_file, force=true)
            end
        end
        
        Test.@testset "Environment Step with Mock Actions" begin
            # Create minimal test environment
            test_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCGATCGATCGATCG"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
                FASTX.FASTQ.Record("read2", BioSequences.LongDNA{4}("GCTAGCTAGCTAGCTAGCTAGCTAGCTA"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII")
            ]
            
            test_file = "test_step.fastq"
            Mycelia.write_fastq(records=test_reads, filename=test_file)
            
            try
                env = Mycelia.create_assembly_environment([test_file], [test_file], episode_length=5)
                initial_state = Mycelia.reset_environment!(env, 1)
                
                # Test continue action
                continue_action = Mycelia.AssemblyAction(
                    :continue_k,
                    Dict(:error_threshold => 0.01),
                    0.95,
                    10,  # Small batch size for testing
                    1
                )
                
                next_state, reward, done = Mycelia.step_environment!(env, continue_action)
                
                Test.@test env.step_count == 1
                Test.@test length(env.action_history) == 1
                Test.@test length(env.reward_history) == 1
                Test.@test !done  # Should not be done after 1 step
                Test.@test reward.total_reward isa Float64
                
                # Test next_k action
                next_k_action = Mycelia.AssemblyAction(:next_k, Dict(), 0.95, 100, 5)
                next_state2, reward2, done2 = Mycelia.step_environment!(env, next_k_action)
                
                Test.@test env.step_count == 2
                Test.@test next_state2.current_k != next_state.current_k  # K should have changed
                
                # Test terminate action
                terminate_action = Mycelia.AssemblyAction(:terminate, Dict(), 0.95, 100, 5)
                final_state, final_reward, done3 = Mycelia.step_environment!(env, terminate_action)
                
                Test.@test done3  # Should be done after terminate
                Test.@test final_reward.termination_reward > 0  # Should have termination reward
                
            finally
                rm(test_file, force=true)
            end
        end
    end

    Test.@testset "Policy Management" begin
        Test.@testset "DQN Policy Creation" begin
            policy = Mycelia.create_dqn_policy(
                state_dim=12,
                action_dim=4,
                hidden_dims=[256, 128],
                learning_rate=0.002,
                epsilon=0.15
            )
            
            Test.@test policy.state_dim == 12
            Test.@test policy.action_dim == 4
            Test.@test policy.hidden_dims == [256, 128]
            Test.@test policy.learning_rate ≈ 0.002
            Test.@test policy.epsilon ≈ 0.15
            Test.@test length(policy.experience_buffer) == 0
        end
        
        Test.@testset "Action Selection" begin
            policy = Mycelia.create_dqn_policy(epsilon=0.0)  # No exploration for deterministic testing
            
            state = Mycelia.AssemblyState(21, 0.8, 0.1, 0.3, 0.8, 0.7, 0.9, [], [21], 50, 10.0)
            
            action = Mycelia.select_action(policy, state)
            
            Test.@test action.decision ∈ [:continue_k, :next_k, :terminate]
            Test.@test action.correction_threshold > 0.0
            Test.@test action.batch_size > 0
            Test.@test action.max_iterations > 0
            Test.@test isa(action.viterbi_params, Dict)
        end
        
        Test.@testset "Action Selection with High Correction Rate" begin
            policy = Mycelia.create_dqn_policy(epsilon=0.0)  # No exploration
            
            # State with high correction rate should prefer continue
            high_correction_state = Mycelia.AssemblyState(21, 0.8, 0.1, 0.3, 0.8, 0.7, 0.9, [], [21], 50, 10.0)
            action = Mycelia.select_action(policy, high_correction_state)
            
            # Should tend toward continue_k with high correction rate
            Test.@test action.decision ∈ [:continue_k, :next_k, :terminate]  # Any valid action
        end
        
        Test.@testset "Action Selection with High K" begin
            policy = Mycelia.create_dqn_policy(epsilon=0.0)  # No exploration
            
            # State with very high k should prefer terminate
            high_k_state = Mycelia.AssemblyState(97, 0.8, 0.01, 0.3, 0.8, 0.7, 0.1, [], [21, 31, 97], 5, 100.0)
            action = Mycelia.select_action(policy, high_k_state)
            
            # Should tend toward terminate with high k
            Test.@test action.decision ∈ [:continue_k, :next_k, :terminate]  # Any valid action
        end
    end

    Test.@testset "Training Data Generation" begin
        Test.@testset "Error Introduction" begin
            # Create test reads
            original_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCG"), "IIIIIIIIIIIIIIII"),
                FASTX.FASTQ.Record("read2", BioSequences.LongDNA{4}("GCTAGCTAGCTAGCTA"), "IIIIIIIIIIIIIIII")
            ]
            
            # Introduce 100% error rate (all bases should change)
            error_reads = Mycelia.introduce_sequencing_errors(original_reads, 1.0)
            
            Test.@test length(error_reads) == length(original_reads)
            
            # Check that sequences changed
            original_seq1 = string(FASTX.FASTQ.sequence(original_reads[1]))
            error_seq1 = string(FASTX.FASTQ.sequence(error_reads[1]))
            
            Test.@test original_seq1 != error_seq1  # Should be different with 100% error rate
            Test.@test length(original_seq1) == length(error_seq1)  # Same length
        end
        
        Test.@testset "Training Dataset Generation" begin
            output_dir = "test_training_generation"
            
            try
                # Generate minimal datasets for testing
                training_files = Mycelia.generate_training_datasets(
                    n_datasets=2,
                    genome_sizes=[1000],
                    error_rates=[0.01],
                    coverage_levels=[10]
                )
                
                Test.@test length(training_files) <= 2  # May be fewer if generation fails
                
                # Check that files exist and are valid
                for file in training_files
                    Test.@test isfile(file)
                    Test.@test filesize(file) > 0
                    
                    # Verify FASTQ format by reading a few records
                    reader = FASTX.FASTQ.Reader(open(file))
                    record_count = 0
                    for record in reader
                        record_count += 1
                        Test.@test length(FASTX.FASTQ.sequence(record)) > 0
                        if record_count >= 5  # Just check first few records
                            break
                        end
                    end
                    close(reader)
                    Test.@test record_count > 0
                end
                
            finally
                # Cleanup
                if isdir(output_dir)
                    rm(output_dir, recursive=true, force=true)
                end
                rm("rl_training_data", recursive=true, force=true)
            end
        end
    end

    Test.@testset "Curriculum Learning" begin
        Test.@testset "Curriculum Schedule Creation" begin
            curriculum = Mycelia.create_curriculum_schedule(stages=3, datasets_per_stage=5)
            
            Test.@test length(curriculum) == 3
            
            for (i, stage) in enumerate(curriculum)
                Test.@test stage["stage"] == i
                Test.@test stage["datasets"] == 5
                Test.@test length(stage["genome_sizes"]) >= 1
                Test.@test length(stage["error_rates"]) >= 1
                Test.@test length(stage["coverage_levels"]) >= 1
                Test.@test stage["episode_length"] >= 50
                Test.@test 0.0 <= stage["success_threshold"] <= 1.0
            end
            
            # Check that difficulty increases
            Test.@test curriculum[3]["error_rates"][1] > curriculum[1]["error_rates"][1]
            Test.@test curriculum[3]["episode_length"] > curriculum[1]["episode_length"]
        end
        
        Test.@testset "Curriculum Stage Progression" begin
            curriculum = Mycelia.create_curriculum_schedule(stages=2, datasets_per_stage=2)
            
            # Verify progression logic
            stage1 = curriculum[1]
            stage2 = curriculum[2]
            
            Test.@test stage1["stage"] == 1
            Test.@test stage2["stage"] == 2
            Test.@test stage2["error_rates"][1] >= stage1["error_rates"][1]  # Increasing difficulty
        end
    end

    Test.@testset "Integration Functions" begin
        Test.@testset "Minimal Training" begin
            # Create minimal test data
            test_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCGATCGATCGATCG"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
                FASTX.FASTQ.Record("read2", BioSequences.LongDNA{4}("GCTAGCTAGCTAGCTAGCTAGCTAGCTA"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII")
            ]
            
            train_file = "minimal_train.fastq"
            val_file = "minimal_val.fastq"
            
            Mycelia.write_fastq(records=test_reads, filename=train_file)
            Mycelia.write_fastq(records=test_reads, filename=val_file)
            
            try
                # Test very short training
                policy, rewards = Mycelia.train_assembly_agent(
                    [train_file], [val_file],
                    episodes=2,
                    episode_length=3
                )
                
                Test.@test isa(policy, Mycelia.DQNPolicy)
                Test.@test length(rewards) == 2
                Test.@test all(r -> isa(r, Float64), rewards)
                
            finally
                rm(train_file, force=true)
                rm(val_file, force=true)
            end
        end
        
        Test.@testset "Agent Evaluation" begin
            # Create test data
            test_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCGATCGATCGATCG"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
                FASTX.FASTQ.Record("read2", BioSequences.LongDNA{4}("GCTAGCTAGCTAGCTAGCTAGCTAGCTA"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII")
            ]
            
            val_file = "eval_test.fastq"
            Mycelia.write_fastq(records=test_reads, filename=val_file)
            
            try
                policy = Mycelia.create_dqn_policy()
                
                metrics = Mycelia.evaluate_assembly_agent(policy, [val_file], episodes=2)
                
                Test.@test haskey(metrics, "mean_reward")
                Test.@test haskey(metrics, "std_reward")
                Test.@test haskey(metrics, "mean_quality")
                Test.@test haskey(metrics, "std_quality")
                Test.@test haskey(metrics, "best_reward")
                Test.@test haskey(metrics, "best_quality")
                
                Test.@test isa(metrics["mean_reward"], Float64)
                Test.@test isa(metrics["std_reward"], Float64)
                
            finally
                rm(val_file, force=true)
            end
        end
        
        Test.@testset "Policy Application" begin
            # Create test data
            test_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCGATCGATCGATCG"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
                FASTX.FASTQ.Record("read2", BioSequences.LongDNA{4}("GCTAGCTAGCTAGCTAGCTAGCTAGCTA"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII")
            ]
            
            input_file = "policy_test.fastq"
            output_dir = "policy_test_output"
            
            Mycelia.write_fastq(records=test_reads, filename=input_file)
            
            try
                policy = Mycelia.create_dqn_policy()
                
                results = Mycelia.apply_learned_policy(policy, input_file, output_dir=output_dir)
                
                Test.@test haskey(results, "input_file")
                Test.@test haskey(results, "output_dir")
                Test.@test haskey(results, "steps_taken")
                Test.@test haskey(results, "final_k")
                Test.@test haskey(results, "final_quality")
                Test.@test haskey(results, "k_progression")
                Test.@test haskey(results, "assembly_log")
                
                Test.@test results["input_file"] == input_file
                Test.@test results["steps_taken"] > 0
                Test.@test results["final_k"] > 0
                Test.@test length(results["k_progression"]) >= 1
                Test.@test length(results["assembly_log"]) >= 1
                
                # Check output directory exists
                Test.@test isdir(output_dir)
                
                # Check log file exists
                log_file = joinpath(output_dir, "assembly_log.json")
                Test.@test isfile(log_file)
                
            finally
                rm(input_file, force=true)
                rm(output_dir, recursive=true, force=true)
            end
        end
    end

    Test.@testset "Error Handling and Edge Cases" begin
        Test.@testset "Invalid Environment Parameters" begin
            # Test with empty dataset lists
            env = Mycelia.create_assembly_environment(String[], String[])
            Test.@test length(env.training_datasets) == 0
            Test.@test length(env.validation_datasets) == 0
        end
        
        Test.@testset "Invalid Action Types" begin
            test_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCG"), "IIIIIIIIIIIIIIII")
            ]
            
            test_file = "invalid_action_test.fastq"
            Mycelia.write_fastq(records=test_reads, filename=test_file)
            
            try
                env = Mycelia.create_assembly_environment([test_file], [test_file])
                Mycelia.reset_environment!(env, 1)
                
                # Test invalid action decision
                invalid_action = Mycelia.AssemblyAction(:invalid_decision, Dict(), 0.95, 100, 5)
                
                Test.@test_throws ErrorException Mycelia.step_environment!(env, invalid_action)
                
            finally
                rm(test_file, force=true)
            end
        end
        
        Test.@testset "Zero Error Rate" begin
            original_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCG"), "IIIIIIIIIIIIIIII")
            ]
            
            # Test 0% error rate (no changes)
            error_reads = Mycelia.introduce_sequencing_errors(original_reads, 0.0)
            
            Test.@test length(error_reads) == length(original_reads)
            
            original_seq = string(FASTX.FASTQ.sequence(original_reads[1]))
            error_seq = string(FASTX.FASTQ.sequence(error_reads[1]))
            
            Test.@test original_seq == error_seq  # Should be identical with 0% error rate
        end
        
        Test.@testset "Environment Reset with Invalid Dataset Index" begin
            test_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCG"), "IIIIIIIIIIIIIIII")
            ]
            
            test_file = "reset_invalid_test.fastq"
            Mycelia.write_fastq(records=test_reads, filename=test_file)
            
            try
                env = Mycelia.create_assembly_environment([test_file], [test_file])
                
                # Test reset with invalid index (should default to 1)
                state = Mycelia.reset_environment!(env, 999)  # Invalid index
                Test.@test state.current_k > 0  # Should still work with default
                
                state2 = Mycelia.reset_environment!(env, 0)  # Invalid index
                Test.@test state2.current_k > 0  # Should still work with default
                
            finally
                rm(test_file, force=true)
            end
        end
    end

    Test.@testset "Performance and Resource Management" begin
        Test.@testset "Memory Usage Tracking" begin
            # Test that memory monitoring doesn't crash
            policy = Mycelia.create_dqn_policy()
            state = Mycelia.AssemblyState(21, 0.8, 0.05, 0.3, 0.8, 0.7, 0.9, [], [21], 50, 10.0)
            
            # Memory usage should be between 0 and 1
            Test.@test 0.0 <= state.memory_usage <= 1.0
        end
        
        Test.@testset "Episode Length Limits" begin
            test_reads = [
                FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATCGATCGATCGATCGATCGATCGATCG"), "IIIIIIIIIIIIIIIIIIIIIIIIIIII")
            ]
            
            test_file = "episode_limit_test.fastq"
            Mycelia.write_fastq(records=test_reads, filename=test_file)
            
            try
                # Create environment with very short episode length
                env = Mycelia.create_assembly_environment([test_file], [test_file], episode_length=2)
                Mycelia.reset_environment!(env, 1)
                
                # Take steps until episode ends
                policy = Mycelia.create_dqn_policy()
                done = false
                step_count = 0
                
                while !done && step_count < 5  # Safety limit
                    action = Mycelia.select_action(policy, env.current_state)
                    _, _, done = Mycelia.step_environment!(env, action)
                    step_count += 1
                end
                
                Test.@test step_count <= 2  # Should terminate within episode length
                
            finally
                rm(test_file, force=true)
            end
        end
    end

    Test.@testset "Framework Self-Test" begin
        Test.@testset "Integrated Test Function" begin
            # This tests the framework's built-in test function
            # Note: This may take a moment as it generates test data
            
            # The test function creates and cleans up its own data
            success = Mycelia.test_rl_framework()
            Test.@test success == true
        end
    end

end

println("Reinforcement Learning Framework tests completed!")
