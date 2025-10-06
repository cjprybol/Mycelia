# """
#     Monte Carlo Tree Search Implementation for Assembly RL
    
# This module implements Monte Carlo Tree Search (MCTS) for genome assembly,
# treating assembly as a sequential decision game. MCTS provides a powerful
# balance between exploration and exploitation through UCB1 selection and
# Monte Carlo rollouts.

# Key advantages of MCTS for assembly:
# - No need for value function approximation
# - Naturally handles large action spaces through sampling
# - Proven success in complex sequential decision problems
# - Can leverage domain knowledge in rollout policy
# """

# ## MCTS-specific types

# """
#     MCTSNode
    
# Node in the Monte Carlo search tree representing an assembly state.
# """
# mutable struct MCTSNode
#     state::AssemblyState
#     parent::Union{Nothing, MCTSNode}
#     children::Dict{AssemblyAction, MCTSNode}
#     visits::Int
#     total_reward::Float64
#     untried_actions::Vector{AssemblyAction}
#     is_terminal::Bool
# end

# function MCTSNode(state::AssemblyState, parent::Union{Nothing, MCTSNode}=nothing)
#     ## Generate available actions for this state
#     available_actions = generate_available_actions(state)
    
#     return MCTSNode(
#         state,
#         parent,
#         Dict{AssemblyAction, MCTSNode}(),
#         0,
#         0.0,
#         available_actions,
#         is_terminal_state(state)
#     )
# end

# """
#     MCTSPolicy
    
# Monte Carlo Tree Search policy for assembly decisions.
# """
# struct MCTSPolicy
#     n_simulations::Int
#     exploration_constant::Float64
#     max_rollout_depth::Int
#     rollout_policy::Function  ## Fast heuristic policy for rollouts
#     use_progressive_widening::Bool
#     widening_constant::Float64
#     widening_exponent::Float64
# end

# function MCTSPolicy(;
#     n_simulations::Int=1000,
#     exploration_constant::Float64=sqrt(2),
#     max_rollout_depth::Int=50,
#     rollout_policy::Function=default_rollout_policy,
#     use_progressive_widening::Bool=true,
#     widening_constant::Float64=10.0,
#     widening_exponent::Float64=0.5
# )
#     return MCTSPolicy(
#         n_simulations,
#         exploration_constant,
#         max_rollout_depth,
#         rollout_policy,
#         use_progressive_widening,
#         widening_constant,
#         widening_exponent
#     )
# end

# ## UCB1 Selection

# """
#     ucb1_score(node::MCTSNode, parent_visits::Int, c::Float64)
    
# Calculate UCB1 score for node selection.
# """
# function ucb1_score(node::MCTSNode, parent_visits::Int, c::Float64)
#     if node.visits == 0
#         return Inf  ## Prioritize unvisited nodes
#     end
    
#     exploitation = node.total_reward / node.visits
#     exploration = c * sqrt(log(parent_visits) / node.visits)
    
#     return exploitation + exploration
# end

# """
#     select_best_child(node::MCTSNode, c::Float64)
    
# Select child with highest UCB1 score.
# """
# function select_best_child(node::MCTSNode, c::Float64)
#     best_score = -Inf
#     best_child = nothing
    
#     for (action, child) in node.children
#         score = ucb1_score(child, node.visits, c)
#         if score > best_score
#             best_score = score
#             best_child = child
#         end
#     end
    
#     return best_child
# end

# ## Tree Policy - Selection and Expansion

# """
#     tree_policy(root::MCTSNode, c::Float64)
    
# Traverse tree using UCB1 until reaching expandable node.
# """
# function tree_policy(root::MCTSNode, c::Float64)
#     node = root
    
#     while !node.is_terminal
#         if !isempty(node.untried_actions)
#             ## Expand with untried action
#             return expand_node(node)
#         else
#             ## Select best child according to UCB1
#             node = select_best_child(node, c)
#             if node === nothing
#                 break
#             end
#         end
#     end
    
#     return node
# end

# """
#     expand_node(node::MCTSNode)
    
# Expand node with one untried action.
# """
# function expand_node(node::MCTSNode)
#     ## Select random untried action
#     action = popfirst!(node.untried_actions)
    
#     ## Generate next state
#     next_state = transition_assembly_state(node.state, action)
    
#     ## Create child node
#     child = MCTSNode(next_state, node)
#     node.children[action] = child
    
#     return child
# end

# ## Default Rollout Policy

# """
#     default_rollout_policy(state::AssemblyState)
    
# Fast heuristic policy for MCTS rollouts.
# Uses simple rules based on assembly metrics.
# """
# function default_rollout_policy(state::AssemblyState)
#     ## Simple heuristic: prefer continuing if quality is improving
#     if state.assembly_quality > 0.8
#         if state.correction_rate > 0.01
#             ## Still making corrections, continue with current k
#             return AssemblyAction(
#                 :continue_k,
#                 Dict(:transition_weight => 1.0, :emission_weight => 1.0, :quality_weight => 0.5),
#                 0.95,
#                 500,
#                 10
#             )
#         else
#             ## Low correction rate, try next k
#             return AssemblyAction(
#                 :next_k,
#                 Dict(:transition_weight => 1.0, :emission_weight => 1.0, :quality_weight => 0.5),
#                 0.99,
#                 1000,
#                 5
#             )
#         end
#     elseif state.memory_usage > 0.9
#         ## Memory pressure, terminate
#         return AssemblyAction(:terminate, Dict(), 0.0, 0, 0)
#     else
#         ## Quality needs improvement, continue with aggressive parameters
#         return AssemblyAction(
#             :continue_k,
#             Dict(:transition_weight => 2.0, :emission_weight => 2.0, :quality_weight => 0.75),
#             0.9,
#             1000,
#             20
#         )
#     end
# end

# ## Rollout Simulation

# """
#     simulate_rollout(state::AssemblyState, policy::Function, max_depth::Int)
    
# Simulate assembly trajectory using rollout policy.
# """
# function simulate_rollout(state::AssemblyState, policy::Function, max_depth::Int)
#     current_state = state
#     total_reward = 0.0
#     depth = 0
    
#     while !is_terminal_state(current_state) && depth < max_depth
#         ## Get action from rollout policy
#         action = policy(current_state)
        
#         ## Transition to next state
#         next_state = transition_assembly_state(current_state, action)
        
#         ## Calculate reward
#         reward = calculate_assembly_reward(next_state)
#         total_reward += reward * (0.99 ^ depth)  ## Discount factor
        
#         current_state = next_state
#         depth += 1
#     end
    
#     return total_reward
# end

# ## Backpropagation

# """
#     backpropagate!(node::MCTSNode, reward::Float64)
    
# Propagate reward up the tree.
# """
# function backpropagate!(node::MCTSNode, reward::Float64)
#     current = node
    
#     while current !== nothing
#         current.visits += 1
#         current.total_reward += reward
#         current = current.parent
#     end
# end

# ## Main MCTS Algorithm

# """
#     mcts_search(root_state::AssemblyState, policy::MCTSPolicy)
    
# Run MCTS to find best action from root state.
# """
# function mcts_search(root_state::AssemblyState, policy::MCTSPolicy)
#     root = MCTSNode(root_state)
    
#     for sim in 1:policy.n_simulations
#         ## Selection - traverse tree to leaf
#         leaf = tree_policy(root, policy.exploration_constant)
        
#         ## Rollout - simulate from leaf
#         reward = simulate_rollout(
#             leaf.state, 
#             policy.rollout_policy, 
#             policy.max_rollout_depth
#         )
        
#         ## Backpropagation - update statistics
#         backpropagate!(leaf, reward)
        
#         ## Progressive widening (optional)
#         if policy.use_progressive_widening
#             apply_progressive_widening!(root, policy)
#         end
#     end
    
#     ## Return action with most visits (robust selection)
#     return select_most_visited_action(root)
# end

# """
#     select_most_visited_action(node::MCTSNode)
    
# Select action with highest visit count (most robust).
# """
# function select_most_visited_action(node::MCTSNode)
#     best_visits = -1
#     best_action = nothing
    
#     for (action, child) in node.children
#         if child.visits > best_visits
#             best_visits = child.visits
#             best_action = action
#         end
#     end
    
#     return best_action
# end

# ## Progressive Widening

# """
#     apply_progressive_widening!(node::MCTSNode, policy::MCTSPolicy)
    
# Add new actions based on visit count (for continuous action spaces).
# """
# function apply_progressive_widening!(node::MCTSNode, policy::MCTSPolicy)
#     ## Calculate maximum allowed actions based on visits
#     max_actions = Int(floor(
#         policy.widening_constant * node.visits ^ policy.widening_exponent
#     ))
    
#     current_actions = length(node.children) + length(node.untried_actions)
    
#     if current_actions < max_actions
#         ## Add new action variant
#         new_action = generate_action_variant(node.state)
#         if new_action !== nothing
#             push!(node.untried_actions, new_action)
#         end
#     end
# end

# """
#     generate_action_variant(state::AssemblyState)
    
# Generate new action variant for progressive widening.
# """
# function generate_action_variant(state::AssemblyState)
#     ## Generate variations in continuous parameters
#     base_decision = rand([:continue_k, :next_k])
    
#     return AssemblyAction(
#         base_decision,
#         Dict(
#             :transition_weight => rand() * 5.0,
#             :emission_weight => rand() * 5.0,
#             :quality_weight => rand()
#         ),
#         0.9 + rand() * 0.09,  ## correction_threshold in [0.9, 0.99]
#         100 + rand(1:900),     ## batch_size in [100, 1000]
#         5 + rand(1:15)         ## max_iterations in [5, 20]
#     )
# end

# ## Helper Functions

# """
#     generate_available_actions(state::AssemblyState)
    
# Generate list of available actions for given state.
# """
# function generate_available_actions(state::AssemblyState)
#     actions = AssemblyAction[]
    
#     ## Always can terminate
#     push!(actions, AssemblyAction(:terminate, Dict(), 0.0, 0, 0))
    
#     ## Can continue with current k if not at memory limit
#     if state.memory_usage < 0.95
#         ## Add several parameter variants
#         for correction_threshold in [0.9, 0.95, 0.99]
#             for quality_weight in [0.0, 0.5, 1.0]
#                 push!(actions, AssemblyAction(
#                     :continue_k,
#                     Dict(
#                         :transition_weight => 1.0,
#                         :emission_weight => 1.0,
#                         :quality_weight => quality_weight
#                     ),
#                     correction_threshold,
#                     500,
#                     10
#                 ))
#             end
#         end
#     end
    
#     ## Can move to next k if available
#     current_k_idx = findfirst(k -> k == state.current_k, state.k_history)
#     if current_k_idx !== nothing && current_k_idx < length(PRIME_KMERS)
#         for quality_weight in [0.25, 0.75]
#             push!(actions, AssemblyAction(
#                 :next_k,
#                 Dict(
#                     :transition_weight => 2.0,
#                     :emission_weight => 2.0,
#                     :quality_weight => quality_weight
#                 ),
#                 0.99,
#                 1000,
#                 5
#             ))
#         end
#     end
    
#     return actions
# end

# """
#     is_terminal_state(state::AssemblyState)
    
# Check if assembly state is terminal.
# """
# function is_terminal_state(state::AssemblyState)
#     return state.memory_usage > 0.95 || 
#            state.correction_rate < 0.001 ||
#            state.current_k >= maximum(PRIME_KMERS)
# end

# ## Training Infrastructure for MCTS

# """
#     MCTSAgent
    
# Agent using MCTS for assembly decisions.
# """
# struct MCTSAgent
#     policy::MCTSPolicy
#     value_estimates::Dict{UInt64, Float64}  ## Cache state value estimates
# end

# function MCTSAgent(; kwargs...)
#     return MCTSAgent(
#         MCTSPolicy(; kwargs...),
#         Dict{UInt64, Float64}()
#     )
# end

# """
#     select_action_mcts(agent::MCTSAgent, state::AssemblyState)
    
# Select action using MCTS.
# """
# function select_action_mcts(agent::MCTSAgent, state::AssemblyState)
#     return mcts_search(state, agent.policy)
# end

# """
#     train_mcts_agent(training_data::Vector{SimulatedAssemblyData}; kwargs...)
    
# Train MCTS agent through self-play on simulated data.
# """
# function train_mcts_agent(
#     training_data::Vector{SimulatedAssemblyData};
#     n_episodes::Int=100,
#     n_simulations::Int=1000,
#     exploration_constant::Float64=sqrt(2),
#     save_path::String="mcts_assembly_agent.jld2"
# )
#     agent = MCTSAgent(
#         n_simulations=n_simulations,
#         exploration_constant=exploration_constant
#     )
    
#     episode_rewards = Float64[]
    
#     for episode in 1:n_episodes
#         ## Select random training sequence
#         sim_data = rand(training_data)
        
#         ## Run episode with MCTS
#         trajectory, total_reward = run_assembly_episode(
#             agent,
#             sim_data.reads,
#             sim_data.reference;
#             method=select_action_mcts
#         )
        
#         push!(episode_rewards, total_reward)
        
#         ## Update value estimates from trajectory
#         for (state, action, reward) in trajectory
#             state_hash = hash(state)
#             current_estimate = get(agent.value_estimates, state_hash, 0.0)
#             ## Incremental mean update
#             n = episode + 1
#             agent.value_estimates[state_hash] = current_estimate + (reward - current_estimate) / n
#         end
        
#         if episode % 10 == 0
#             @info "MCTS Training" episode avg_reward=Statistics.mean(episode_rewards[max(1,episode-9):episode])
#         end
#     end
    
#     ## Save trained agent
#     JLD2.save(save_path, Dict(
#         "policy" => agent.policy,
#         "value_estimates" => agent.value_estimates,
#         "training_rewards" => episode_rewards
#     ))
    
#     return agent, episode_rewards
# end

# ## Integration with Assembly Pipeline

# """
#     intelligent_assembly_mcts(
#         reads::Vector{String};
#         agent::Union{Nothing, MCTSAgent}=nothing,
#         n_simulations::Int=1000,
#         kwargs...
#     )
    
# Run intelligent assembly using MCTS for decision making.
# """
# function intelligent_assembly_mcts(
#     reads::Vector{String};
#     agent::Union{Nothing, MCTSAgent}=nothing,
#     n_simulations::Int=1000,
#     initial_k::Int=19,
#     memory_limit_gb::Float64=32.0,
#     output_dir::String=Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS_mcts_assembly"),
#     save_intermediates::Bool=true,
#     kwargs...
# )
#     ## Create agent if not provided
#     if agent === nothing
#         agent = MCTSAgent(n_simulations=n_simulations)
#     end
    
#     ## Initialize assembly environment
#     env = AssemblyEnvironment(
#         reads,
#         initial_k=initial_k,
#         memory_limit=memory_limit_gb * 1e9
#     )
    
#     ## Create output directory
#     mkpath(output_dir)
    
#     ## Assembly loop
#     step = 0
#     trajectory = []
    
#     while !is_done(env)
#         state = get_state(env)
        
#         ## Use MCTS to select action
#         action = select_action_mcts(agent, state)
        
#         ## Log decision
#         @info "MCTS Decision" step k=state.current_k action=action.decision visits=length(agent.value_estimates)
        
#         ## Take action in environment
#         next_state, reward, done, info = step!(env, action)
        
#         ## Save trajectory
#         push!(trajectory, (state, action, reward, next_state))
        
#         ## Save intermediate results
#         if save_intermediates
#             save_assembly_checkpoint(
#                 joinpath(output_dir, "checkpoint_step_$(step).jld2"),
#                 env,
#                 trajectory
#             )
#         end
        
#         step += 1
#     end
    
#     ## Extract final assembly
#     final_assembly = get_assembly(env)
    
#     ## Save results
#     save_assembly_results(
#         joinpath(output_dir, "final_assembly.fasta"),
#         final_assembly,
#         trajectory,
#         agent
#     )
    
#     return final_assembly, trajectory
# end

# ## Comparison with Other RL Methods

# """
#     compare_mcts_performance(test_data::Vector{SimulatedAssemblyData})
    
# Compare MCTS with other RL methods on test data.
# """
# function compare_mcts_performance(test_data::Vector{SimulatedAssemblyData})
#     results = Dict{String, Vector{Float64}}()
    
#     ## Test MCTS with different configurations
#     for (name, n_sims) in [
#         ("MCTS-100", 100),
#         ("MCTS-500", 500),
#         ("MCTS-1000", 1000),
#         ("MCTS-5000", 5000)
#     ]
#         agent = MCTSAgent(n_simulations=n_sims)
#         rewards = Float64[]
        
#         for sim_data in test_data
#             _, trajectory = intelligent_assembly_mcts(
#                 sim_data.reads,
#                 agent=agent
#             )
            
#             total_reward = sum(r for (_, _, r, _) in trajectory)
#             push!(rewards, total_reward)
#         end
        
#         results[name] = rewards
#     end
    
#     return results
# end

# ## Export main functions
# export MCTSNode, MCTSPolicy, MCTSAgent
# export mcts_search, train_mcts_agent, intelligent_assembly_mcts
# export compare_mcts_performance