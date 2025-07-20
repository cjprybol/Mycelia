"""
Test suite for comparing reinforcement learning approaches
"""

# Example: Using the three RL approaches

println("Loading reinforcement learning modules...")

# Test 1: Basic functionality of each approach
Test.@testset "Individual RL Approaches" begin
    # Generate small test dataset
    test_reads = ["ATCGATCG", "TCGATCGA", "CGATCGAT", "GATCGATC", "ATCGATCG"]
    
    Test.@testset "Custom RL" begin
        # Test the existing custom implementation
        agent = Mycelia.DQNPolicy()
        env = Mycelia.AssemblyEnvironment(test_reads)
        state = Mycelia.reset!(env)
        
        Test.@test isa(state, Mycelia.AssemblyState)
        Test.@test state.current_k == 19
        
        # Test action selection
        action = Mycelia.select_action(agent, state)
        Test.@test isa(action, Mycelia.AssemblyAction)
        Test.@test action.decision in [:continue_k, :next_k, :terminate]
    end
    
    Test.@testset "ReinforcementLearning.jl" begin
        # Test RL.jl environment
        env = Mycelia.AssemblyEnvRL()
        state = ReinforcementLearning.state(env)
        
        Test.@test isa(state, Mycelia.AssemblyState)
        Test.@test length(ReinforcementLearning.action_space(env)) == 3
        
        # Test action execution
        ReinforcementLearning.act!(env, :continue_k)
        Test.@test !isempty(env.reward_history)
    end
    
    Test.@testset "POMDPs.jl" begin
        # Test POMDPs.jl MDP
        mdp = Mycelia.AssemblyMDP()
        
        Test.@test POMDPs.discount(mdp) == 0.99
        Test.@test mdp.initial_k == 19
        
        # Test initial state
        init_state = Random.rand(POMDPs.initialstate(mdp))
        Test.@test isa(init_state, Mycelia.AssemblyState)
        Test.@test init_state.current_k == 19
    end
end

# Test 2: Comparison framework
Test.@testset "RL Comparison Framework" begin
    # Small dataset for quick testing
    test_reads = [
        "ATCGATCGATCG",
        "TCGATCGATCGA", 
        "CGATCGATCGAT",
        "GATCGATCGATC",
        "ATCGATCGATCG"
    ]
    
    # Test individual wrapper functions
    Test.@testset "Wrapper Functions" begin
        # Custom RL wrapper
        custom_history = Mycelia.run_custom_rl(test_reads, verbose=false)
        Test.@test isa(custom_history, Mycelia.UnifiedAssemblyHistory)
        Test.@test custom_history.approach == "custom"
        Test.@test length(custom_history.k_values) > 0
        
        # Note: Full RL.jl and POMDPs.jl tests would require more setup
        # and longer runtime, so we test the structure instead
    end
    
    # Test comparison metrics calculation
    Test.@testset "Comparison Metrics" begin
        # Create mock histories for testing
        mock_results = Dict(
            :custom => Mycelia.UnifiedAssemblyHistory(
                "custom",
                [19, 23, 29],          # k_values
                [0.5, 0.7, 0.9],       # quality_scores
                [0.1, 0.2, 0.3],       # correction_rates
                [0.2, 0.3, 0.4],       # memory_usage
                [:continue_k, :next_k, :terminate],  # actions
                [10.0, 20.0, 50.0],    # rewards
                [0.1, 0.2, 0.1],       # time_per_step
                1.5,                    # total_time
                0.9,                    # final_quality
                1000                    # final_assembly_length
            ),
            :rljl => Mycelia.UnifiedAssemblyHistory(
                "ReinforcementLearning.jl (dqn)",
                [19, 23, 29, 31],
                [0.4, 0.6, 0.8, 0.95],
                [0.1, 0.15, 0.25, 0.35],
                [0.2, 0.25, 0.35, 0.45],
                [:continue_k, :next_k, :next_k, :terminate],
                [8.0, 15.0, 30.0, 60.0],
                [0.15, 0.2, 0.25, 0.15],
                2.0,
                0.95,
                1050
            )
        )
        
        metrics = Mycelia.calculate_comparison_metrics(mock_results)
        
        Test.@test haskey(metrics, "final_quality")
        Test.@test haskey(metrics, "average_quality")
        Test.@test haskey(metrics, "total_time")
        Test.@test haskey(metrics, "convergence_speed")
        
        Test.@test metrics["final_quality"]["custom"] == 0.9
        Test.@test metrics["final_quality"]["rljl"] == 0.95
        Test.@test metrics["total_time"]["custom"] == 1.5
        Test.@test metrics["k_values_tried"]["custom"] == 3
        Test.@test metrics["k_values_tried"]["rljl"] == 4
    end
end

# Test 3: Integration test
Test.@testset "Full Comparison Integration" begin
    # Very small dataset for integration test
    mini_reads = ["ATCG", "TCGA", "CGAT", "GATC"]
    
    # Run comparison with only custom approach for speed
    comparison = Mycelia.compare_rl_approaches(
        mini_reads,
        approaches=[:custom],
        dataset_name="test_integration",
        training_episodes=10,
        verbose=false,
        save_results=false
    )
    
    Test.@test isa(comparison, Mycelia.RLComparison)
    Test.@test comparison.dataset_name == "test_integration"
    Test.@test !isnothing(comparison.custom_results)
    Test.@test !isempty(comparison.comparison_metrics)
end

# Test 4: Demonstrate usage patterns
Test.@testset "Usage Examples" begin
    println("\n" * "="^60)
    println("Example Usage of RL Approaches")
    println("="^60)
    
    # Generate example data
    example_reads = [
        "ATCGATCGATCGATCG",
        "TCGATCGATCGATCGA",
        "CGATCGATCGATCGAT",
        "GATCGATCGATCGATC",
        "ATCGATCGATCGATCG"
    ]
    
    println("\n1. Custom RL Approach (Rule-based placeholder):")
    println("   - Fast execution")
    println("   - No training required")
    println("   - Good baseline for comparison")
    
    println("\n2. ReinforcementLearning.jl Approach:")
    println("   - Access to DQN, PPO, A2C, etc.")
    println("   - Experience replay and neural networks")
    println("   - Well-tested implementations")
    
    println("\n3. POMDPs.jl Approach:")
    println("   - Formal MDP/POMDP specification")
    println("   - Value iteration, MCTS, belief tracking")
    println("   - Handles partial observability")
    
    println("\n" * "="^60)
end

# Example: How to use each approach individually
println("\nExample code for using each approach:")
println("""
# Custom approach
agent = Mycelia.DQNPolicy()
env = Mycelia.AssemblyEnvironment(reads)
# ... run episodes ...

# ReinforcementLearning.jl approach
assembly, history = Mycelia.intelligent_assembly_rljl(
    reads,
    algorithm=:dqn,
    train_episodes=1000
)

# POMDPs.jl approach
assembly, history = Mycelia.intelligent_assembly_pomdp(
    reads,
    solver=:value_iteration,
    use_pomdp=false  # Set true for partial observability
)

# Compare all approaches
comparison = Mycelia.compare_rl_approaches(
    reads,
    approaches=[:custom, :rljl, :pomdp],
    rljl_algorithm=:dqn,
    pomdp_solver=:mcts
)
""")

println("\nAll RL comparison tests completed!")