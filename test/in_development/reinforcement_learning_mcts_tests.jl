# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/in_development/reinforcement_learning_mcts_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/in_development/reinforcement_learning_mcts_tests.jl", "test/in_development", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

"""
Tests for Monte Carlo Tree Search (MCTS) implementation
"""

Test.@testset "MCTS Assembly Tests" begin
    
    ## Test basic MCTS node operations
    Test.@testset "MCTS Node Operations" begin
        state = AssemblyState(
            current_k=31,
            assembly_quality=0.85,
            correction_rate=0.05,
            memory_usage=0.3,
            graph_connectivity=0.8,
            coverage_uniformity=0.9,
            error_signal_clarity=0.7,
            reward_history=[0.5, 0.6, 0.7],
            k_history=[19, 23, 31],
            step_count=3,
            total_reward=1.8
        )
        
        ## Test node creation
        node = MCTSNode(state)
        Test.@test node.state == state
        Test.@test node.parent === nothing
        Test.@test isempty(node.children)
        Test.@test node.visits == 0
        Test.@test node.total_reward == 0.0
        Test.@test !isempty(node.untried_actions)
        Test.@test !node.is_terminal
        
        ## Test UCB1 scoring
        node.visits = 10
        node.total_reward = 5.0
        parent_visits = 100
        c = sqrt(2)
        
        score = ucb1_score(node, parent_visits, c)
        expected_exploitation = 5.0 / 10
        expected_exploration = c * sqrt(log(100) / 10)
        Test.@test score ≈ expected_exploitation + expected_exploration
        
        ## Test unvisited node priority
        unvisited_node = MCTSNode(state)
        Test.@test ucb1_score(unvisited_node, parent_visits, c) == Inf
    end
    
    ## Test tree policy and expansion
    Test.@testset "Tree Policy and Expansion" begin
        root_state = AssemblyState(
            current_k=19,
            assembly_quality=0.5,
            correction_rate=0.1,
            memory_usage=0.1,
            graph_connectivity=0.6,
            coverage_uniformity=0.7,
            error_signal_clarity=0.5,
            reward_history=Float64[],
            k_history=[19],
            step_count=0,
            total_reward=0.0
        )
        
        root = MCTSNode(root_state)
        
        ## Test expansion
        initial_untried = length(root.untried_actions)
        child = expand_node(root)
        
        Test.@test child.parent === root
        Test.@test length(root.untried_actions) == initial_untried - 1
        Test.@test length(root.children) == 1
        
        ## Test tree policy selection
        leaf = tree_policy(root, sqrt(2))
        Test.@test leaf !== root  ## Should have expanded
    end
    
    ## Test rollout policy
    Test.@testset "Rollout Policy" begin
        ## High quality state
        high_quality_state = AssemblyState(
            current_k=31,
            assembly_quality=0.9,
            correction_rate=0.02,
            memory_usage=0.3,
            graph_connectivity=0.9,
            coverage_uniformity=0.95,
            error_signal_clarity=0.8,
            reward_history=[0.8, 0.85, 0.9],
            k_history=[19, 23, 31],
            step_count=3,
            total_reward=2.55
        )
        
        action = default_rollout_policy(high_quality_state)
        Test.@test action.decision == :continue_k
        
        ## Low correction rate state
        low_correction_state = AssemblyState(
            current_k=31,
            assembly_quality=0.85,
            correction_rate=0.005,
            memory_usage=0.3,
            graph_connectivity=0.9,
            coverage_uniformity=0.95,
            error_signal_clarity=0.8,
            reward_history=[0.8, 0.85, 0.9],
            k_history=[19, 23, 31],
            step_count=3,
            total_reward=2.55
        )
        
        action = default_rollout_policy(low_correction_state)
        Test.@test action.decision == :next_k
        
        ## Memory pressure state
        memory_pressure_state = AssemblyState(
            current_k=71,
            assembly_quality=0.7,
            correction_rate=0.05,
            memory_usage=0.95,
            graph_connectivity=0.8,
            coverage_uniformity=0.85,
            error_signal_clarity=0.7,
            reward_history=[0.6, 0.65, 0.7],
            k_history=[19, 23, 31, 37, 53, 71],
            step_count=6,
            total_reward=3.0
        )
        
        action = default_rollout_policy(memory_pressure_state)
        Test.@test action.decision == :terminate
    end
    
    ## Test backpropagation
    Test.@testset "Backpropagation" begin
        ## Create a simple tree
        root_state = create_test_assembly_state(k=19, quality=0.5)
        root = MCTSNode(root_state)
        
        child1 = expand_node(root)
        child2 = expand_node(child1)
        
        ## Test backpropagation
        reward = 0.8
        backpropagate!(child2, reward)
        
        Test.@test child2.visits == 1
        Test.@test child2.total_reward == reward
        Test.@test child1.visits == 1
        Test.@test child1.total_reward == reward
        Test.@test root.visits == 1
        Test.@test root.total_reward == reward
    end
    
    ## Test MCTS search
    Test.@testset "MCTS Search" begin
        initial_state = AssemblyState(
            current_k=19,
            assembly_quality=0.5,
            correction_rate=0.1,
            memory_usage=0.1,
            graph_connectivity=0.6,
            coverage_uniformity=0.7,
            error_signal_clarity=0.5,
            reward_history=Float64[],
            k_history=[19],
            step_count=0,
            total_reward=0.0
        )
        
        policy = MCTSPolicy(
            n_simulations=100,
            exploration_constant=sqrt(2),
            max_rollout_depth=10
        )
        
        ## Run MCTS search
        best_action = mcts_search(initial_state, policy)
        
        Test.@test best_action isa AssemblyAction
        Test.@test best_action.decision in [:continue_k, :next_k, :terminate]
    end
    
    ## Test progressive widening
    Test.@testset "Progressive Widening" begin
        state = create_test_assembly_state(k=31, quality=0.7)
        node = MCTSNode(state)
        
        policy = MCTSPolicy(
            use_progressive_widening=true,
            widening_constant=10.0,
            widening_exponent=0.5
        )
        
        ## Simulate visits to trigger widening
        initial_actions = length(node.untried_actions)
        node.visits = 10
        
        apply_progressive_widening!(node, policy)
        
        ## Should have added new actions based on visit count
        max_actions = Int(floor(10.0 * 10^0.5))
        Test.@test length(node.untried_actions) <= max_actions
    end
    
    ## Test MCTS agent training
    Test.@testset "MCTS Agent Training" begin
        ## Generate small training dataset
        training_data = [
            SimulatedAssemblyData(
                "sim_1",
                generate_random_reads(100, 50),
                "ACGT" ^ 25,  ## 100bp reference
                0.01,
                10.0,
                "linear"
            )
        ]
        
        ## Train agent
        agent, rewards = train_mcts_agent(
            training_data;
            n_episodes=5,
            n_simulations=50
        )
        
        Test.@test agent isa MCTSAgent
        Test.@test length(rewards) == 5
        Test.@test !isempty(agent.value_estimates)
    end
    
    ## Test intelligent assembly with MCTS
    Test.@testset "Intelligent Assembly MCTS" begin
        ## Generate test reads
        reference = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        reads = simulate_reads(reference, coverage=5.0, read_length=20, error_rate=0.01)
        
        ## Create agent
        agent = MCTSAgent(n_simulations=50)
        
        ## Run assembly
        assembly, trajectory = intelligent_assembly_mcts(
            reads;
            agent=agent,
            initial_k=11,
            memory_limit_gb=1.0,
            save_intermediates=false
        )
        
        Test.@test !isempty(assembly)
        Test.@test !isempty(trajectory)
        Test.@test all(t -> t[1] isa AssemblyState, trajectory)
        Test.@test all(t -> t[2] isa AssemblyAction, trajectory)
        Test.@test all(t -> t[3] isa Float64, trajectory)
    end
    
    ## Test MCTS comparison performance
    Test.@testset "MCTS Performance Comparison" begin
        ## Create test data
        test_data = [
            SimulatedAssemblyData(
                "test_1",
                generate_random_reads(50, 30),
                "ACGT" ^ 15,
                0.01,
                5.0,
                "linear"
            )
        ]
        
        ## Compare different MCTS configurations
        results = compare_mcts_performance(test_data)
        
        Test.@test haskey(results, "MCTS-100")
        Test.@test haskey(results, "MCTS-500")
        Test.@test haskey(results, "MCTS-1000")
        Test.@test all(r -> length(r) == 1, values(results))
    end
    
    ## Test integration with comparison framework
    Test.@testset "MCTS in Comparison Framework" begin
        reads = generate_random_reads(100, 50)
        
        ## Run comparison including MCTS
        comparison = compare_rl_approaches(
            reads;
            approaches=[:custom, :mcts],
            training_episodes=5,
            mcts_n_simulations=50,
            verbose=false,
            save_results=false
        )
        
        Test.@test comparison.custom_results !== nothing
        Test.@test comparison.mcts_results !== nothing
        Test.@test comparison.mcts_results.approach == "MCTS (n_sim=50)"
    end
end

## Helper function for creating test states
function create_test_assembly_state(; k=31, quality=0.7)
    return AssemblyState(
        current_k=k,
        assembly_quality=quality,
        correction_rate=0.05,
        memory_usage=0.3,
        graph_connectivity=0.8,
        coverage_uniformity=0.9,
        error_signal_clarity=0.7,
        reward_history=Float64[],
        k_history=[k],
        step_count=1,
        total_reward=quality
    )
end
