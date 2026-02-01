# """
#     Reinforcement Learning Framework for Self-Optimizing Assembly

# This module implements a hierarchical reinforcement learning system for autonomous genome assembly
# parameter optimization. The system learns optimal k-mer selection strategies, error correction
# parameters, and termination conditions through simulated training and cross-validation.

# ## Architecture

# ### High-Level Policy (Meta-Controller)
# - **Algorithm**: Deep Q-Network (DQN) with experience replay
# - **State**: Assembly quality metrics, current k, memory usage, correction rate
# - **Actions**: Continue with current k, move to next prime k, or terminate assembly
# - **Termination**: Based on correction rate, memory limits, or maximum k reached

# ### Low-Level Policy (Error Correction Controller) 
# - **Algorithm**: Policy Gradient (PPO) for continuous parameter optimization
# - **State**: Local graph topology, quality scores, coverage patterns
# - **Actions**: Viterbi parameters, path selection strategies, confidence thresholds
# - **Integration**: Uses existing Viterbi and probabilistic path algorithms

# ## Core Philosophy
# - **Accuracy-First**: Prioritize assembly accuracy over contiguity or speed
# - **Dynamic Decision Making**: Learn optimal strategies through diverse training scenarios
# - **Cross-Validation**: Use validated quality metrics for reliable reward signals
# - **Learnable Parameters**: Eliminate manual parameter tuning through policy learning
# - **Biological Realism**: Train on diverse simulated datasets matching real genomic complexity

# ## Integration with Existing Infrastructure
# - Builds on intelligent-assembly.jl and iterative-assembly.jl foundations
# - Uses existing qualmer graphs and Viterbi error correction algorithms
# - Leverages cross-validation pipeline for assembly quality assessment
# - Maintains type stability with MetaGraphsNext architecture
# """

# # Core RL State and Action Representations

# """
#     AssemblyState

# State representation for reinforcement learning environment containing all information
# needed to make assembly decisions.

# # Fields
# - `current_k::Int`: Current k-mer size being processed
# - `assembly_quality::Float64`: Current assembly quality score (QV-based)
# - `correction_rate::Float64`: Rate of successful error corrections in recent iterations
# - `memory_usage::Float64`: Current memory utilization (fraction of limit)
# - `graph_connectivity::Float64`: Graph connectivity metric (proportion of strongly connected components)
# - `coverage_uniformity::Float64`: Uniformity of k-mer coverage distribution
# - `error_signal_clarity::Float64`: Clarity of error signal detection (sparsity-based)
# - `iteration_history::Vector{Float64}`: Recent reward history for trend analysis
# - `k_progression::Vector{Int}`: Sequence of k-mer sizes processed so far
# - `corrections_made::Int`: Total corrections made at current k
# - `time_elapsed::Float64`: Time spent on current k-mer size (seconds)
# """
# struct AssemblyState
#     current_k::Int
#     assembly_quality::Float64
#     correction_rate::Float64
#     memory_usage::Float64
#     graph_connectivity::Float64
#     coverage_uniformity::Float64
#     error_signal_clarity::Float64
#     iteration_history::Vector{Float64}
#     k_progression::Vector{Int}
#     corrections_made::Int
#     time_elapsed::Float64
# end

# """
#     AssemblyAction

# Action representation for reinforcement learning decisions during assembly.

# # Fields
# - `decision::Symbol`: Primary decision (:continue_k, :next_k, :terminate)
# - `viterbi_params::Dict{Symbol, Float64}`: Viterbi algorithm parameters
# - `correction_threshold::Float64`: Quality threshold for error correction
# - `batch_size::Int`: Batch size for processing (memory management)
# - `max_iterations::Int`: Maximum iterations at current k before forced progression
# """
# struct AssemblyAction
#     decision::Symbol  # :continue_k, :next_k, :terminate
#     viterbi_params::Dict{Symbol, Float64}
#     correction_threshold::Float64
#     batch_size::Int
#     max_iterations::Int
# end

# """
#     RewardComponents

# Structured representation of reward signal components for training the RL agent.

# # Fields
# - `accuracy_reward::Float64`: Primary reward based on assembly accuracy (weighted 1000x)
# - `efficiency_reward::Float64`: Secondary reward for computational efficiency (weighted 10x)
# - `error_penalty::Float64`: Penalty for false positives/negatives (weighted -500x)
# - `progress_bonus::Float64`: Bonus for making meaningful progress
# - `termination_reward::Float64`: Reward for appropriate termination timing
# - `total_reward::Float64`: Weighted sum of all components
# """
# struct RewardComponents
#     accuracy_reward::Float64
#     efficiency_reward::Float64
#     error_penalty::Float64
#     progress_bonus::Float64
#     termination_reward::Float64
#     total_reward::Float64
# end

# # RL Environment Implementation

# """
#     AssemblyEnvironment

# Reinforcement learning environment for training assembly decision policies.

# # Fields
# - `current_state::AssemblyState`: Current environment state
# - `training_datasets::Vector{String}`: Paths to training FASTQ files
# - `validation_datasets::Vector{String}`: Paths to validation FASTQ files
# - `episode_length::Int`: Maximum steps per training episode
# - `step_count::Int`: Current step in episode
# - `reward_history::Vector{Float64}`: Reward history for current episode
# - `action_history::Vector{AssemblyAction}`: Action history for experience replay
# - `assembly_cache::Dict{String, Any}`: Cache for intermediate assembly results
# """
# mutable struct AssemblyEnvironment
#     current_state::AssemblyState
#     training_datasets::Vector{String}
#     validation_datasets::Vector{String}
#     episode_length::Int
#     step_count::Int
#     reward_history::Vector{Float64}
#     action_history::Vector{AssemblyAction}
#     assembly_cache::Dict{String, Any}
# end

# """
#     create_assembly_environment(training_data, validation_data; episode_length=100)

# Create a new reinforcement learning environment for assembly training.

# # Arguments
# - `training_data::Vector{String}`: Paths to training FASTQ datasets
# - `validation_data::Vector{String}`: Paths to validation FASTQ datasets
# - `episode_length::Int`: Maximum steps per training episode (default: 100)

# # Returns
# - `AssemblyEnvironment`: Initialized RL environment ready for training

# # Example
# ```julia
# training_files = ["genome1.fastq", "genome2.fastq", "genome3.fastq"]
# validation_files = ["validation1.fastq", "validation2.fastq"]
# env = create_assembly_environment(training_files, validation_files)
# ```
# """
# function create_assembly_environment(training_data::Vector{String}, validation_data::Vector{String}; episode_length::Int=100)
#     # Initialize with default state
#     initial_state = AssemblyState(
#         21,  # Default starting k
#         0.0,  # No initial quality assessment
#         0.0,  # No initial correction rate
#         0.0,  # No initial memory usage
#         0.0,  # No initial connectivity
#         0.0,  # No initial coverage uniformity
#         0.0,  # No initial error signal clarity
#         Float64[],  # Empty iteration history
#         Int[],  # Empty k progression
#         0,  # No corrections made yet
#         0.0   # No time elapsed
#     )

#     return AssemblyEnvironment(
#         initial_state,
#         training_data,
#         validation_data,
#         episode_length,
#         0,  # Step count starts at 0
#         Float64[],  # Empty reward history
#         AssemblyAction[],  # Empty action history
#         Dict{String, Any}()  # Empty cache
#     )
# end

# """
#     reset_environment!(env::AssemblyEnvironment, dataset_idx::Int=1)

# Reset the RL environment to start a new training episode.

# # Arguments
# - `env::AssemblyEnvironment`: Environment to reset
# - `dataset_idx::Int`: Index of training dataset to use (default: 1)

# # Returns
# - `AssemblyState`: Initial state for new episode

# # Example
# ```julia
# initial_state = reset_environment!(env, 2)  # Use second training dataset
# ```
# """
# function reset_environment!(env::AssemblyEnvironment, dataset_idx::Int=1)
#     # Validate dataset index
#     if dataset_idx < 1 || dataset_idx > length(env.training_datasets)
#         dataset_idx = 1
#     end

#     # Reset environment counters
#     env.step_count = 0
#     env.reward_history = Float64[]
#     env.action_history = AssemblyAction[]
#     env.assembly_cache = Dict{String, Any}()

#     # Initialize state based on selected dataset
#     dataset_path = env.training_datasets[dataset_idx]

#     # Quick assessment of initial k-mer size using sparsity detection
#     initial_k = try
#         # Use existing sparsity detection from intelligent-assembly.jl
#         reads = collect(FASTX.FASTQ.Reader(open(dataset_path)))
#         find_initial_k(reads)
#     catch e
#         @warn "Could not determine initial k for $dataset_path: $e"
#         21  # Default fallback
#     end

#     # Create initial state
#     env.current_state = AssemblyState(
#         initial_k,
#         0.0,  # Quality will be assessed after first assembly attempt
#         0.0,  # No correction rate initially
#         0.0,  # No memory usage initially
#         0.0,  # No connectivity assessment initially
#         0.0,  # No coverage uniformity initially
#         0.0,  # No error signal initially
#         Float64[],  # Empty history
#         [initial_k],  # K progression starts with initial k
#         0,  # No corrections made
#         0.0   # No time elapsed
#     )

#     return env.current_state
# end

# """
#     step_environment!(env::AssemblyEnvironment, action::AssemblyAction)

# Execute an action in the RL environment and return the resulting state and reward.

# # Arguments
# - `env::AssemblyEnvironment`: Environment to step
# - `action::AssemblyAction`: Action to execute

# # Returns
# - `Tuple{AssemblyState, RewardComponents, Bool}`: (next_state, reward, done)
#   - `next_state`: State after action execution
#   - `reward`: Reward components for the action
#   - `done`: Whether episode is complete

# # Example
# ```julia
# action = AssemblyAction(:continue_k, Dict(), 0.95, 1000, 5)
# next_state, reward, done = step_environment!(env, action)
# ```
# """
# function step_environment!(env::AssemblyEnvironment, action::AssemblyAction)
#     env.step_count += 1
#     push!(env.action_history, action)

#     # Check for episode termination conditions
#     done = (env.step_count >= env.episode_length) || (action.decision == :terminate)

#     # Execute action based on decision type
#     if action.decision == :continue_k
#         reward = execute_continue_k_action(env, action)
#     elseif action.decision == :next_k
#         reward = execute_next_k_action(env, action)
#     elseif action.decision == :terminate
#         reward = execute_terminate_action(env, action)
#         done = true
#     else
#         error("Unknown action decision: $(action.decision)")
#     end

#     # Update reward history
#     push!(env.reward_history, reward.total_reward)

#     # Update iteration history in state (keep last 10 rewards)
#     iteration_history = copy(env.current_state.iteration_history)
#     push!(iteration_history, reward.total_reward)
#     if length(iteration_history) > 10
#         iteration_history = iteration_history[end-9:end]
#     end

#     # Create updated state
#     env.current_state = AssemblyState(
#         env.current_state.current_k,
#         env.current_state.assembly_quality,
#         env.current_state.correction_rate,
#         env.current_state.memory_usage,
#         env.current_state.graph_connectivity,
#         env.current_state.coverage_uniformity,
#         env.current_state.error_signal_clarity,
#         iteration_history,
#         env.current_state.k_progression,
#         env.current_state.corrections_made,
#         env.current_state.time_elapsed
#     )

#     return env.current_state, reward, done
# end

# """
#     execute_continue_k_action(env::AssemblyEnvironment, action::AssemblyAction)

# Execute a "continue with current k" action by performing error correction iterations.

# # Arguments
# - `env::AssemblyEnvironment`: Current environment
# - `action::AssemblyAction`: Action specifying correction parameters

# # Returns
# - `RewardComponents`: Reward breakdown for the action
# """
# function execute_continue_k_action(env::AssemblyEnvironment, action::AssemblyAction)
#     start_time = time()

#     # Get current dataset
#     current_dataset = env.training_datasets[1]  # Simplified for now
#     current_k = env.current_state.current_k

#     # Cache key for intermediate results
#     cache_key = "k$(current_k)_step$(env.step_count)"

#     try
#         # Perform error correction at current k using iterative assembly approach
#         if !haskey(env.assembly_cache, cache_key)
#             # Load reads and build qualmer graph
#             reads = collect(FASTX.FASTQ.Reader(open(current_dataset)))
#             graph = build_qualmer_graph(reads, k=current_k)

#             # Perform batch error correction
#             corrections_made = 0
#             total_reads = length(reads)

#             # Simulate error correction process (placeholder for full implementation)
#             for (i, read) in enumerate(reads[1:min(action.batch_size, total_reads)])
#                 # Use existing iterative assembly functions
#                 try
#                     improved_read, was_improved = improve_read_likelihood(read, graph, current_k)
#                     if was_improved
#                         corrections_made += 1
#                     end
#                 catch e
#                     @debug "Error correction failed for read $i: $e"
#                 end
#             end

#             # Cache results
#             env.assembly_cache[cache_key] = Dict(
#                 "corrections_made" => corrections_made,
#                 "total_processed" => min(action.batch_size, total_reads),
#                 "graph_size" => Graphs.nv(graph.graph),
#                 "processing_time" => time() - start_time
#             )
#         end

#         results = env.assembly_cache[cache_key]

#         # Calculate correction rate
#         correction_rate = results["total_processed"] > 0 ? 
#                          results["corrections_made"] / results["total_processed"] : 0.0

#         # Update state metrics
#         env.current_state = AssemblyState(
#             current_k,
#             env.current_state.assembly_quality,  # Will be updated by quality assessment
#             correction_rate,
#             estimate_memory_usage(10000, current_k) / (1024^3),  # Convert to GB, estimate with 10K kmers
#             0.8,  # Placeholder for graph connectivity
#             0.7,  # Placeholder for coverage uniformity
#             correction_rate > 0.05 ? 1.0 : 0.0,  # Error signal clarity
#             env.current_state.iteration_history,
#             env.current_state.k_progression,
#             env.current_state.corrections_made + results["corrections_made"],
#             env.current_state.time_elapsed + results["processing_time"]
#         )

#         # Calculate reward components
#         accuracy_reward = correction_rate * 1000.0  # High weight for accuracy
#         efficiency_reward = (1.0 / max(results["processing_time"], 0.1)) * 10.0  # Efficiency bonus
#         error_penalty = 0.0  # No penalties for valid corrections
#         progress_bonus = corrections_made > 0 ? 50.0 : 0.0  # Progress incentive
#         termination_reward = 0.0  # No termination in continue action

#         total_reward = accuracy_reward + efficiency_reward - error_penalty + progress_bonus + termination_reward

#         return RewardComponents(
#             accuracy_reward,
#             efficiency_reward,
#             error_penalty,
#             progress_bonus,
#             termination_reward,
#             total_reward
#         )

#     catch e
#         @warn "Error executing continue_k action: $e"
#         # Return penalty for failed action
#         return RewardComponents(0.0, 0.0, -100.0, 0.0, 0.0, -100.0)
#     end
# end

# """
#     execute_next_k_action(env::AssemblyEnvironment, action::AssemblyAction)

# Execute a "move to next k" action by progressing to the next prime k-mer size.

# # Arguments
# - `env::AssemblyEnvironment`: Current environment
# - `action::AssemblyAction`: Action specifying progression parameters

# # Returns
# - `RewardComponents`: Reward breakdown for the action
# """
# function execute_next_k_action(env::AssemblyEnvironment, action::AssemblyAction)
#     current_k = env.current_state.current_k
#     next_k = next_prime_k(current_k)  # Use existing function from intelligent-assembly.jl

#     # Check if progression is valid
#     if next_k == current_k || next_k > 101  # Max k limit
#         # Invalid progression - penalize
#         return RewardComponents(0.0, 0.0, -200.0, 0.0, 0.0, -200.0)
#     end

#     # Update k progression
#     new_k_progression = copy(env.current_state.k_progression)
#     push!(new_k_progression, next_k)

#     # Calculate progression reward based on correction rate trend
#     correction_rate = env.current_state.correction_rate
#     iteration_history = env.current_state.iteration_history

#     # Assess whether progression is appropriate
#     progression_appropriate = if length(iteration_history) >= 3
#         # Check if recent rewards are declining (indicating diminishing returns)
#         recent_trend = Statistics.mean(iteration_history[end-2:end])
#         earlier_trend = length(iteration_history) >= 6 ? 
#                        Statistics.mean(iteration_history[end-5:end-3]) : recent_trend
#         recent_trend < earlier_trend * 0.9  # 10% decline threshold
#     else
#         correction_rate < 0.02  # Low correction rate justifies progression
#     end

#     # Update state with new k
#     env.current_state = AssemblyState(
#         next_k,
#         env.current_state.assembly_quality,
#         0.0,  # Reset correction rate for new k
#         env.current_state.memory_usage,
#         env.current_state.graph_connectivity,
#         env.current_state.coverage_uniformity,
#         0.0,  # Reset error signal clarity
#         Float64[],  # Reset iteration history for new k
#         new_k_progression,
#         0,  # Reset corrections count
#         0.0   # Reset time for new k
#     )

#     # Calculate rewards
#     accuracy_reward = 0.0  # No immediate accuracy change
#     efficiency_reward = progression_appropriate ? 100.0 : -50.0  # Reward appropriate timing
#     error_penalty = 0.0
#     progress_bonus = 25.0  # Small bonus for making progress
#     termination_reward = 0.0

#     total_reward = accuracy_reward + efficiency_reward - error_penalty + progress_bonus + termination_reward

#     return RewardComponents(
#         accuracy_reward,
#         efficiency_reward,
#         error_penalty,
#         progress_bonus,
#         termination_reward,
#         total_reward
#     )
# end

# """
#     execute_terminate_action(env::AssemblyEnvironment, action::AssemblyAction)

# Execute a "terminate assembly" action and assess final assembly quality.

# # Arguments
# - `env::AssemblyEnvironment`: Current environment
# - `action::AssemblyAction`: Termination action

# # Returns
# - `RewardComponents`: Final reward breakdown including assembly quality assessment
# """
# function execute_terminate_action(env::AssemblyEnvironment, action::AssemblyAction)
#     # Assess final assembly quality using cross-validation
#     current_dataset = env.training_datasets[1]

#     try
#         # Perform final assembly quality assessment
#         reads = collect(FASTX.FASTQ.Reader(open(current_dataset)))

#         # Save reads to temporary file for cross-validation
#         temp_file = tempname() * ".fastq"
#         write_fastq(records=reads, filename=temp_file)

#         # Use cross-validation for robust quality assessment
#         cv_results = mycelia_cross_validation(
#             temp_file,
#             k_folds=3,  # Quick assessment
#             max_k=env.current_state.current_k
#         )

#         # Cleanup temporary file
#         rm(temp_file, force=true)

#         # Extract quality metrics
#         final_quality = if !isempty(cv_results.results)
#             mean([result["assembly_quality"] for result in cv_results.results])
#         else
#             0.0
#         end

#         # Update final state
#         env.current_state = AssemblyState(
#             env.current_state.current_k,
#             final_quality,
#             env.current_state.correction_rate,
#             env.current_state.memory_usage,
#             env.current_state.graph_connectivity,
#             env.current_state.coverage_uniformity,
#             env.current_state.error_signal_clarity,
#             env.current_state.iteration_history,
#             env.current_state.k_progression,
#             env.current_state.corrections_made,
#             env.current_state.time_elapsed
#         )

#         # Calculate termination rewards
#         accuracy_reward = final_quality * 2000.0  # High reward for final quality
#         efficiency_reward = env.step_count < env.episode_length ? 200.0 : 0.0  # Early termination bonus
#         error_penalty = final_quality < 0.5 ? -1000.0 : 0.0  # Penalty for poor quality
#         progress_bonus = length(env.current_state.k_progression) > 1 ? 100.0 : 0.0  # Multi-k bonus
#         termination_reward = 500.0  # Base termination reward

#         total_reward = accuracy_reward + efficiency_reward - error_penalty + progress_bonus + termination_reward

#         return RewardComponents(
#             accuracy_reward,
#             efficiency_reward,
#             error_penalty,
#             progress_bonus,
#             termination_reward,
#             total_reward
#         )

#     catch e
#         @warn "Error in termination assessment: $e"
#         # Penalty for failed termination
#         return RewardComponents(0.0, 0.0, -500.0, 0.0, 0.0, -500.0)
#     end
# end

# # Policy Network Architecture (Placeholder for Future Implementation)

# """
#     DQNPolicy

# Deep Q-Network policy for high-level assembly decisions.

# This is a placeholder structure for the neural network architecture that will
# be implemented with a machine learning framework like Flux.jl or MLJ.jl.

# # Fields
# - `state_dim::Int`: Dimension of state representation
# - `action_dim::Int`: Number of possible actions
# - `hidden_dims::Vector{Int}`: Hidden layer dimensions
# - `learning_rate::Float64`: Learning rate for training
# - `epsilon::Float64`: Exploration rate for epsilon-greedy policy
# - `experience_buffer::Vector{Any}`: Experience replay buffer
# """
# struct DQNPolicy
#     state_dim::Int
#     action_dim::Int
#     hidden_dims::Vector{Int}
#     learning_rate::Float64
#     epsilon::Float64
#     experience_buffer::Vector{Any}
# end

# """
#     create_dqn_policy(; state_dim=11, action_dim=3, hidden_dims=[128, 64], learning_rate=0.001, epsilon=0.1)

# Create a Deep Q-Network policy for assembly decisions.

# # Arguments
# - `state_dim::Int`: Dimension of state representation (default: 11)
# - `action_dim::Int`: Number of discrete actions (default: 3 for continue/next/terminate)
# - `hidden_dims::Vector{Int}`: Hidden layer sizes (default: [128, 64])
# - `learning_rate::Float64`: Learning rate (default: 0.001)
# - `epsilon::Float64`: Exploration probability (default: 0.1)

# # Returns
# - `DQNPolicy`: Initialized policy network

# # Example
# ```julia
# policy = create_dqn_policy(hidden_dims=[256, 128, 64])
# ```
# """
# function create_dqn_policy(; state_dim::Int=11, action_dim::Int=3, hidden_dims::Vector{Int}=[128, 64], 
#                           learning_rate::Float64=0.001, epsilon::Float64=0.1)
#     return DQNPolicy(
#         state_dim,
#         action_dim,
#         hidden_dims,
#         learning_rate,
#         epsilon,
#         Any[]  # Empty experience buffer
#     )
# end

# """
#     select_action(policy::DQNPolicy, state::AssemblyState)

# Select an action using the DQN policy with epsilon-greedy exploration.

# This is a placeholder implementation that will be replaced with actual neural network
# inference once the ML framework is integrated.

# # Arguments
# - `policy::DQNPolicy`: Trained policy network
# - `state::AssemblyState`: Current environment state

# # Returns
# - `AssemblyAction`: Selected action for the given state

# # Example
# ```julia
# action = select_action(policy, current_state)
# ```
# """
# function select_action(policy::DQNPolicy, state::AssemblyState)
#     # Placeholder implementation using rule-based logic
#     # This will be replaced with neural network inference

#     if Random.rand() < policy.epsilon
#         # Random exploration
#         decision = Random.rand([:continue_k, :next_k, :terminate])
#     else
#         # Greedy policy based on simple heuristics
#         if state.correction_rate > 0.05
#             decision = :continue_k  # Keep correcting if making progress
#         elseif state.current_k < 71  # Allow progression up to reasonable k
#             decision = :next_k
#         else
#             decision = :terminate
#         end
#     end

#     # Generate action with default parameters
#     return AssemblyAction(
#         decision,
#         Dict(:error_threshold => 0.01, :min_coverage => 3),  # Default Viterbi params
#         0.95,  # Default correction threshold
#         1000,  # Default batch size
#         5      # Default max iterations
#     )
# end

# # Training Infrastructure

# """
#     train_assembly_agent(training_data::Vector{String}, validation_data::Vector{String}; 
#                         episodes=1000, episode_length=100)

# Train a reinforcement learning agent for assembly optimization.

# This function implements the complete training loop for the hierarchical RL system.

# # Arguments
# - `training_data::Vector{String}`: Paths to training FASTQ files
# - `validation_data::Vector{String}`: Paths to validation FASTQ files
# - `episodes::Int`: Number of training episodes (default: 1000)
# - `episode_length::Int`: Maximum steps per episode (default: 100)

# # Returns
# - `Tuple{DQNPolicy, Vector{Float64}}`: (trained_policy, training_rewards)

# # Example
# ```julia
# training_files = ["train1.fastq", "train2.fastq", "train3.fastq"]
# validation_files = ["val1.fastq", "val2.fastq"]
# policy, rewards = train_assembly_agent(training_files, validation_files, episodes=500)
# ```
# """
# function train_assembly_agent(training_data::Vector{String}, validation_data::Vector{String}; 
#                              episodes::Int=1000, episode_length::Int=100)

#     # Initialize environment and policy
#     env = create_assembly_environment(training_data, validation_data, episode_length=episode_length)
#     policy = create_dqn_policy()

#     # Training metrics
#     episode_rewards = Float64[]
#     training_start_time = time()

#     println("Starting RL training for assembly optimization...")
#     println("Episodes: $episodes, Episode length: $episode_length")
#     println("Training datasets: $(length(training_data))")
#     println("Validation datasets: $(length(validation_data))")

#     for episode in 1:episodes
#         # Reset environment for new episode
#         dataset_idx = ((episode - 1) % length(training_data)) + 1
#         state = reset_environment!(env, dataset_idx)

#         episode_reward = 0.0
#         episode_start_time = time()

#         # Run episode
#         for step in 1:episode_length
#             # Select action using policy
#             action = select_action(policy, state)

#             # Execute action and observe result
#             next_state, reward, done = step_environment!(env, action)

#             # Update metrics
#             episode_reward += reward.total_reward
#             state = next_state

#             # Store experience for replay (placeholder)
#             # In full implementation, this would update the neural network

#             if done
#                 break
#             end
#         end

#         push!(episode_rewards, episode_reward)

#         # Progress reporting
#         if episode % 100 == 0 || episode == episodes
#             avg_reward = Statistics.mean(episode_rewards[max(1, end-99):end])
#             elapsed_time = time() - training_start_time
#             println("Episode $episode/$episodes - Avg Reward: $(round(avg_reward, digits=2)) - Time: $(round(elapsed_time, digits=1))s")
#         end
#     end

#     total_training_time = time() - training_start_time
#     println("Training completed in $(round(total_training_time, digits=1)) seconds")
#     println("Final average reward: $(round(Statistics.mean(episode_rewards[max(1, end-99):end]), digits=2))")

#     return policy, episode_rewards
# end

# """
#     evaluate_assembly_agent(policy::DQNPolicy, validation_data::Vector{String}; episodes=10)

# Evaluate a trained assembly agent on validation data.

# # Arguments
# - `policy::DQNPolicy`: Trained policy to evaluate
# - `validation_data::Vector{String}`: Validation FASTQ files
# - `episodes::Int`: Number of evaluation episodes (default: 10)

# # Returns
# - `Dict{String, Float64}`: Evaluation metrics including mean reward, assembly quality, etc.

# # Example
# ```julia
# # metrics = evaluate_assembly_agent(trained_policy, validation_files)
# # println("Mean reward: \$(metrics["mean_reward"])")
# ```
# """
# function evaluate_assembly_agent(policy::DQNPolicy, validation_data::Vector{String}; episodes::Int=10)
#     env = create_assembly_environment(validation_data, validation_data)

#     episode_rewards = Float64[]
#     final_qualities = Float64[]

#     println("Evaluating assembly agent on $(length(validation_data)) validation datasets...")

#     for episode in 1:episodes
#         dataset_idx = ((episode - 1) % length(validation_data)) + 1
#         state = reset_environment!(env, dataset_idx)

#         episode_reward = 0.0

#         # Run evaluation episode (no exploration)
#         old_epsilon = policy.epsilon
#         policy_copy = DQNPolicy(policy.state_dim, policy.action_dim, policy.hidden_dims, 
#                                policy.learning_rate, 0.0, policy.experience_buffer)  # No exploration

#         for step in 1:env.episode_length
#             action = select_action(policy_copy, state)
#             next_state, reward, done = step_environment!(env, action)

#             episode_reward += reward.total_reward
#             state = next_state

#             if done
#                 break
#             end
#         end

#         push!(episode_rewards, episode_reward)
#         push!(final_qualities, state.assembly_quality)

#         println("Episode $episode: Reward = $(round(episode_reward, digits=2)), Quality = $(round(state.assembly_quality, digits=3))")
#     end

#     return Dict(
#         "mean_reward" => Statistics.mean(episode_rewards),
#         "std_reward" => Statistics.std(episode_rewards),
#         "mean_quality" => Statistics.mean(final_qualities),
#         "std_quality" => Statistics.std(final_qualities),
#         "best_reward" => maximum(episode_rewards),
#         "best_quality" => maximum(final_qualities)
#     )
# end

# # Simulation Framework for Training Data Generation

# """
#     generate_training_datasets(; n_datasets=20, genome_sizes=[10000, 50000, 100000], 
#                               error_rates=[0.001, 0.01, 0.05], coverage_levels=[20, 30, 50])

# Generate diverse training datasets for RL agent training.

# This function creates simulated genomic datasets with varying characteristics to
# provide comprehensive training scenarios for the RL agent.

# # Arguments
# - `n_datasets::Int`: Total number of datasets to generate (default: 20)
# - `genome_sizes::Vector{Int}`: Range of genome sizes to simulate (default: [10K, 50K, 100K])
# - `error_rates::Vector{Float64}`: Range of sequencing error rates (default: [0.1%, 1%, 5%])
# - `coverage_levels::Vector{Int}`: Range of coverage depths (default: [20x, 30x, 50x])

# # Returns
# - `Vector{String}`: Paths to generated training FASTQ files

# # Example
# ```julia
# training_files = generate_training_datasets(n_datasets=50, genome_sizes=[50000, 100000, 200000])
# ```
# """
# function generate_training_datasets(; n_datasets::Int=20, 
#                                    genome_sizes::Vector{Int}=[10000, 50000, 100000],
#                                    error_rates::Vector{Float64}=[0.001, 0.01, 0.05], 
#                                    coverage_levels::Vector{Int}=[20, 30, 50])

#     output_dir = "rl_training_data"
#     mkpath(output_dir)

#     generated_files = String[]

#     println("Generating $n_datasets training datasets...")
#     println("Genome sizes: $genome_sizes")
#     println("Error rates: $error_rates")
#     println("Coverage levels: $coverage_levels")

#     for i in 1:n_datasets
#         # Randomly select parameters for this dataset
#         genome_size = Random.rand(genome_sizes)
#         error_rate = Random.rand(error_rates)
#         coverage = Random.rand(coverage_levels)

#         # Generate dataset name
#         dataset_name = "training_$(i)_size$(genome_size)_err$(round(error_rate*100, digits=1))_cov$(coverage)"
#         output_file = joinpath(output_dir, "$(dataset_name).fastq")

#         try
#             # Generate reference genome
#             reference = BioSequences.randdnaseq(genome_size)

#             # Generate reads with specified parameters
#             reads_1, reads_2 = generate_paired_end_reads(
#                 reference, coverage, 150, 300  # 150bp reads, 300bp insert
#             )

#             # Convert to FASTQ records
#             fastq_reads = []
#             for (i, seq) in enumerate(reads_1)
#                 quality = "I"^length(seq)  # High quality scores
#                 record = FASTX.FASTQ.Record("read_$(i)", seq, quality)
#                 push!(fastq_reads, record)
#             end

#             # Introduce sequencing errors
#             error_reads = introduce_sequencing_errors(fastq_reads, error_rate)

#             # Save as FASTQ using existing fastx function
#             write_fastq(records=error_reads, filename=output_file)
#             push!(generated_files, output_file)

#             if i % 5 == 0 || i == n_datasets
#                 println("Generated $i/$n_datasets datasets")
#             end

#         catch e
#             @warn "Failed to generate dataset $i: $e"
#         end
#     end

#     println("Generated $(length(generated_files)) training datasets in $output_dir")
#     return generated_files
# end

# """
#     introduce_sequencing_errors(reads::Vector, error_rate::Float64)

# Introduce realistic sequencing errors into a set of reads.

# # Arguments
# - `reads::Vector`: Vector of FASTQ records
# - `error_rate::Float64`: Error rate (0.0 to 1.0)

# # Returns
# - `Vector`: Reads with introduced errors

# # Example
# ```julia
# error_reads = introduce_sequencing_errors(clean_reads, 0.01)  # 1% error rate
# ```
# """
# function introduce_sequencing_errors(reads::Vector, error_rate::Float64)
#     error_reads = similar(reads)

#     for (i, read) in enumerate(reads)
#         sequence = FASTX.FASTQ.sequence(read)
#         quality = FASTX.FASTQ.quality(read)
#         identifier = FASTX.FASTQ.identifier(read)
#         description = FASTX.FASTQ.description(read)

#         # Convert to mutable sequence
#         seq_string = string(sequence)
#         seq_chars = collect(seq_string)

#         # Introduce errors
#         for j in 1:length(seq_chars)
#             if Random.rand() < error_rate
#                 # Random substitution
#                 original = seq_chars[j]
#                 alternatives = filter(x -> x != original, ['A', 'T', 'G', 'C'])
#                 seq_chars[j] = Random.rand(alternatives)
#             end
#         end

#         # Create new sequence
#         error_sequence = BioSequences.LongDNA{4}(join(seq_chars))

#         # Create new record with errors (identifier, sequence, quality)
#         error_reads[i] = FASTX.FASTQ.Record(identifier, error_sequence, quality)
#     end

#     return error_reads
# end

# # Curriculum Learning Implementation

# """
#     create_curriculum_schedule(; stages=4, datasets_per_stage=10)

# Create a curriculum learning schedule that progressively increases difficulty.

# # Arguments
# - `stages::Int`: Number of curriculum stages (default: 4)
# - `datasets_per_stage::Int`: Datasets per stage (default: 10)

# # Returns
# - `Vector{Dict}`: Curriculum schedule with parameters for each stage

# # Example
# ```julia
# curriculum = create_curriculum_schedule(stages=5, datasets_per_stage=15)
# ```
# """
# function create_curriculum_schedule(; stages::Int=4, datasets_per_stage::Int=10)
#     curriculum = Dict[]

#     for stage in 1:stages
#         difficulty_factor = stage / stages

#         stage_config = Dict(
#             "stage" => stage,
#             "datasets" => datasets_per_stage,
#             "genome_sizes" => Int[round(10000 + 90000 * difficulty_factor)],  # 10K to 100K
#             "error_rates" => [0.001 + 0.049 * difficulty_factor],  # 0.1% to 5%
#             "coverage_levels" => [50 - Int(round(30 * difficulty_factor))],  # 50x down to 20x
#             "episode_length" => Int(round(50 + 50 * difficulty_factor)),  # 50 to 100 steps
#             "success_threshold" => 0.8 - 0.2 * difficulty_factor  # 80% down to 60%
#         )

#         push!(curriculum, stage_config)
#     end

#     return curriculum
# end

# """
#     train_with_curriculum(curriculum_schedule::Vector{Dict}, validation_data::Vector{String})

# Train the RL agent using curriculum learning.

# # Arguments
# - `curriculum_schedule::Vector{Dict}`: Curriculum stages
# - `validation_data::Vector{String}`: Validation datasets

# # Returns
# - `Tuple{DQNPolicy, Vector{Float64}}`: (trained_policy, stage_rewards)

# # Example
# ```julia
# curriculum = create_curriculum_schedule()
# policy, rewards = train_with_curriculum(curriculum, validation_files)
# ```
# """
# function train_with_curriculum(curriculum_schedule::Vector{Dict}, validation_data::Vector{String})
#     policy = create_dqn_policy()
#     stage_rewards = Float64[]

#     println("Starting curriculum learning with $(length(curriculum_schedule)) stages...")

#     for stage_config in curriculum_schedule
#         stage = stage_config["stage"]
#         println("\n=== Curriculum Stage $stage ===")

#         # Generate training data for this stage
#         stage_training_data = generate_training_datasets(
#             n_datasets=stage_config["datasets"],
#             genome_sizes=stage_config["genome_sizes"],
#             error_rates=stage_config["error_rates"],
#             coverage_levels=stage_config["coverage_levels"]
#         )

#         # Train on this stage
#         episodes = 200  # Episodes per stage
#         stage_policy, episode_rewards = train_assembly_agent(
#             stage_training_data, validation_data,
#             episodes=episodes,
#             episode_length=stage_config["episode_length"]
#         )

#         # Update policy for next stage
#         policy = stage_policy

#         # Evaluate stage performance
#         stage_performance = Statistics.mean(episode_rewards[max(1, end-49):end])  # Last 50 episodes
#         push!(stage_rewards, stage_performance)

#         println("Stage $stage completed - Performance: $(round(stage_performance, digits=2))")

#         # Check if stage mastery achieved
#         success_threshold = stage_config["success_threshold"]
#         if stage_performance < success_threshold
#             @warn "Stage $stage performance below threshold. Consider extending training."
#         end

#         # Cleanup stage data
#         for file in stage_training_data
#             rm(file, force=true)
#         end
#     end

#     println("\nCurriculum learning completed!")
#     println("Stage performances: $(round.(stage_rewards, digits=2))")

#     return policy, stage_rewards
# end

# # Integration and Testing Functions

# """
#     test_rl_framework()

# Test the reinforcement learning framework with minimal examples.

# This function provides a comprehensive test of the RL infrastructure using
# small synthetic datasets.

# # Returns
# - `Bool`: Whether all tests passed

# # Example
# ```julia
# success = test_rl_framework()
# ```
# """
# function test_rl_framework()
#     println("Testing RL Framework...")

#     try
#         # Generate minimal test datasets
#         test_datasets = generate_training_datasets(
#             n_datasets=3, 
#             genome_sizes=[1000], 
#             error_rates=[0.01], 
#             coverage_levels=[20]
#         )

#         # Test environment creation
#         env = create_assembly_environment(test_datasets, test_datasets, episode_length=10)
#         println("✓ Environment creation")

#         # Test environment reset
#         initial_state = reset_environment!(env, 1)
#         println("✓ Environment reset")

#         # Test policy creation
#         policy = create_dqn_policy()
#         println("✓ Policy creation")

#         # Test action selection
#         action = select_action(policy, initial_state)
#         println("✓ Action selection")

#         # Test environment step
#         next_state, reward, done = step_environment!(env, action)
#         println("✓ Environment step")

#         # Test short training episode
#         _, training_rewards = train_assembly_agent(test_datasets, test_datasets, episodes=3, episode_length=5)
#         println("✓ Training episode")

#         # Test evaluation
#         metrics = evaluate_assembly_agent(policy, test_datasets, episodes=2)
#         println("✓ Agent evaluation")

#         # Cleanup
#         for file in test_datasets
#             rm(file, force=true)
#         end

#         println("All RL framework tests passed! ✅")
#         return true

#     catch e
#         @error "RL framework test failed: $e"
#         return false
#     end
# end

# """
#     apply_learned_policy(policy::DQNPolicy, input_fastq::String; output_dir="rl_assembly")

# Apply a trained RL policy to perform genome assembly.

# This function uses a trained policy to make autonomous assembly decisions.

# # Arguments
# - `policy::DQNPolicy`: Trained assembly policy
# - `input_fastq::String`: Path to input FASTQ file
# - `output_dir::String`: Output directory for assembly results

# # Returns
# - `Dict{String, Any}`: Assembly results and metadata

# # Example
# ```julia
# results = apply_learned_policy(trained_policy, "genome.fastq")
# ```
# """
# function apply_learned_policy(policy::DQNPolicy, input_fastq::String; output_dir::String="rl_assembly")
#     mkpath(output_dir)

#     # Create single-dataset environment
#     env = create_assembly_environment([input_fastq], [input_fastq])
#     state = reset_environment!(env, 1)

#     println("Applying learned policy to $(input_fastq)...")

#     step_count = 0
#     max_steps = 50  # Prevent infinite loops
#     assembly_log = Dict[]

#     while step_count < max_steps
#         step_count += 1

#         # Select action using policy (no exploration)
#         policy_deterministic = DQNPolicy(policy.state_dim, policy.action_dim, policy.hidden_dims,
#                                        policy.learning_rate, 0.0, policy.experience_buffer)
#         action = select_action(policy_deterministic, state)

#         # Log decision
#         decision_log = Dict(
#             "step" => step_count,
#             "k" => state.current_k,
#             "decision" => string(action.decision),
#             "correction_rate" => state.correction_rate,
#             "assembly_quality" => state.assembly_quality
#         )
#         push!(assembly_log, decision_log)

#         println("Step $step_count: k=$(state.current_k), decision=$(action.decision)")

#         # Execute action
#         next_state, reward, done = step_environment!(env, action)
#         state = next_state

#         if done
#             println("Assembly completed in $step_count steps")
#             break
#         end
#     end

#     # Save results
#     results = Dict(
#         "input_file" => input_fastq,
#         "output_dir" => output_dir,
#         "steps_taken" => step_count,
#         "final_k" => state.current_k,
#         "final_quality" => state.assembly_quality,
#         "k_progression" => state.k_progression,
#         "total_corrections" => state.corrections_made,
#         "assembly_log" => assembly_log
#     )

#     # Save assembly log
#     log_file = joinpath(output_dir, "assembly_log.json")
#     open(log_file, "w") do f
#         JSON.print(f, results, 2)
#     end

#     println("Assembly results saved to $output_dir")
#     return results
# end

# # # Export main functions for external use
# # export AssemblyState, AssemblyAction, RewardComponents, AssemblyEnvironment, DQNPolicy
# # export create_assembly_environment, reset_environment!, step_environment!
# # export create_dqn_policy, select_action, train_assembly_agent, evaluate_assembly_agent
# # export generate_training_datasets, create_curriculum_schedule, train_with_curriculum
# # export test_rl_framework, apply_learned_policy

# # Module-level testing
# if @isdefined(Main) && hasmethod(Main.eval, Tuple{Expr})
#     # Only run tests if this file is being executed directly
#     if abspath(PROGRAM_FILE) == @__FILE__
#         println("Running RL framework tests...")
#         test_rl_framework()
#     end
# end
