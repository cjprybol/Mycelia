"""
    Reinforcement Learning Comparison Framework

This module provides a unified interface to compare four RL implementations:
1. Custom implementation (rule-based placeholder)
2. ReinforcementLearning.jl implementation
3. POMDPs.jl implementation
4. Monte Carlo Tree Search (MCTS) implementation

Enables systematic evaluation of different RL approaches for assembly optimization.
"""

# Unified Assembly History for Comparison
struct UnifiedAssemblyHistory
    approach::String
    k_values::Vector{Int}
    quality_scores::Vector{Float64}
    correction_rates::Vector{Float64}
    memory_usage::Vector{Float64}
    actions_taken::Vector{Symbol}
    rewards::Vector{Float64}
    time_per_step::Vector{Float64}
    total_time::Float64
    final_quality::Float64
    final_assembly_length::Int
end

"""
    RLComparison

Results from comparing different RL approaches on the same dataset.
"""
struct RLComparison
    dataset_name::String
    custom_results::UnifiedAssemblyHistory
    rljl_results::Union{UnifiedAssemblyHistory, Nothing}
    pomdp_results::Union{UnifiedAssemblyHistory, Nothing}
    mcts_results::Union{UnifiedAssemblyHistory, Nothing}
    comparison_metrics::Dict{String, Any}
end

"""
    compare_rl_approaches(
        reads::Vector{String};
        approaches::Vector{Symbol}=[:custom, :rljl, :pomdp, :mcts],
        kwargs...
    )

Compare different RL approaches on the same dataset.
"""
function compare_rl_approaches(
    reads::Vector{String};
    approaches::Vector{Symbol}=[:custom, :rljl, :pomdp, :mcts],
    dataset_name::String="assembly_comparison",
    training_episodes::Int=100,
    rljl_algorithm::Symbol=:dqn,
    pomdp_solver::Symbol=:value_iteration,
    mcts_n_simulations::Int=1000,
    verbose::Bool=true,
    save_results::Bool=true,
    output_dir::String="rl_comparison_results"
)
    if save_results && !isdir(output_dir)
        mkdir(output_dir)
    end
    
    results = Dict{Symbol, UnifiedAssemblyHistory}()
    
    # Run each approach
    for approach in approaches
        if verbose
            println("\n" * "="^50)
            println("Running $approach approach...")
            println("="^50)
        end
        
        start_time = Dates.now()
        
        if approach == :custom
            # Use existing custom implementation
            history = run_custom_rl(reads; verbose=verbose)
            results[:custom] = history
            
        elseif approach == :rljl
            # Use ReinforcementLearning.jl implementation
            history = run_rljl_rl(
                reads;
                algorithm=rljl_algorithm,
                training_episodes=training_episodes,
                verbose=verbose
            )
            results[:rljl] = history
            
        elseif approach == :pomdp
            # Use POMDPs.jl implementation
            history = run_pomdp_rl(
                reads;
                solver=pomdp_solver,
                training_episodes=training_episodes,
                verbose=verbose
            )
            results[:pomdp] = history
            
        elseif approach == :mcts
            # Use Monte Carlo Tree Search implementation
            history = run_mcts_rl(
                reads;
                n_simulations=mcts_n_simulations,
                training_episodes=training_episodes,
                verbose=verbose
            )
            results[:mcts] = history
        end
        
        elapsed = Dates.now() - start_time
        if verbose
            println("Completed $approach in $(elapsed)")
        end
    end
    
    # Calculate comparison metrics
    comparison_metrics = calculate_comparison_metrics(results)
    
    # Create comparison object
    comparison = RLComparison(
        dataset_name,
        get(results, :custom, nothing),
        get(results, :rljl, nothing),
        get(results, :pomdp, nothing),
        get(results, :mcts, nothing),
        comparison_metrics
    )
    
    # Save results if requested
    if save_results
        save_comparison_results(comparison, output_dir)
    end
    
    # Generate comparison plots
    if verbose
        plot_comparison_results(comparison, output_dir)
    end
    
    return comparison
end

# Wrapper functions for each approach

"""
    run_custom_rl(reads::Vector{String}; kwargs...)

Run the custom RL implementation and return unified history.
"""
function run_custom_rl(reads::Vector{String}; 
    initial_k::Int=19,
    verbose::Bool=false
)
    start_time = Dates.now()
    time_per_step = Float64[]
    
    # Use existing custom implementation
    agent = DQNPolicy()
    env = AssemblyEnvironment(reads)
    
    k_values = Int[]
    quality_scores = Float64[]
    correction_rates = Float64[]
    memory_usage = Float64[]
    actions_taken = Symbol[]
    rewards = Float64[]
    
    state = reset!(env)
    done = false
    
    while !done
        step_start = Dates.now()
        
        # Get action from custom policy
        action = select_action(agent, state)
        
        # Execute action
        next_state, reward, done, info = step!(env, action)
        
        # Record data
        push!(k_values, state.current_k)
        push!(quality_scores, state.assembly_quality)
        push!(correction_rates, state.correction_rate)
        push!(memory_usage, state.memory_usage)
        push!(actions_taken, action.decision)
        push!(rewards, reward)
        
        step_time = Dates.value(Dates.now() - step_start) / 1000.0
        push!(time_per_step, step_time)
        
        state = next_state
        
        if verbose && length(k_values) % 10 == 0
            println("Step $(length(k_values)): K=$(state.current_k), " *
                   "Quality=$(round(state.assembly_quality, digits=3))")
        end
    end
    
    total_time = Dates.value(Dates.now() - start_time) / 1000.0
    
    return UnifiedAssemblyHistory(
        "custom",
        k_values,
        quality_scores,
        correction_rates,
        memory_usage,
        actions_taken,
        rewards,
        time_per_step,
        total_time,
        length(quality_scores) > 0 ? quality_scores[end] : 0.0,
        0  # Placeholder for assembly length
    )
end

"""
    run_rljl_rl(reads::Vector{String}; kwargs...)

Run the ReinforcementLearning.jl implementation and return unified history.
"""
function run_rljl_rl(reads::Vector{String};
    algorithm::Symbol=:dqn,
    training_episodes::Int=100,
    initial_k::Int=19,
    verbose::Bool=false
)
    start_time = Dates.now()
    
    # Train agent
    training_reads = reads[1:min(1000, length(reads))]
    agent, hook = train_assembly_agent_rljl(
        training_reads,
        episodes=training_episodes,
        algorithm=algorithm
    )
    
    # Apply to full dataset
    env = AssemblyEnvRL(initial_k=initial_k)
    ReinforcementLearning.reset!(env)
    
    k_values = Int[]
    quality_scores = Float64[]
    correction_rates = Float64[]
    memory_usage = Float64[]
    actions_taken = Symbol[]
    rewards = Float64[]
    time_per_step = Float64[]
    
    while !ReinforcementLearning.is_terminated(env)
        step_start = Dates.now()
        
        state = ReinforcementLearning.state(env)
        action = ReinforcementLearning.plan!(agent, state)
        
        ReinforcementLearning.act!(env, action)
        
        # Record data
        push!(k_values, state.current_k)
        push!(quality_scores, state.assembly_quality)
        push!(correction_rates, state.correction_rate)
        push!(memory_usage, state.memory_usage)
        push!(actions_taken, action)
        push!(rewards, ReinforcementLearning.reward(env))
        
        step_time = Dates.value(Dates.now() - step_start) / 1000.0
        push!(time_per_step, step_time)
        
        if verbose && length(k_values) % 10 == 0
            println("Step $(length(k_values)): K=$(state.current_k), " *
                   "Quality=$(round(state.assembly_quality, digits=3))")
        end
    end
    
    total_time = Dates.value(Dates.now() - start_time) / 1000.0
    
    return UnifiedAssemblyHistory(
        "ReinforcementLearning.jl ($algorithm)",
        k_values,
        quality_scores,
        correction_rates,
        memory_usage,
        actions_taken,
        rewards,
        time_per_step,
        total_time,
        length(quality_scores) > 0 ? quality_scores[end] : 0.0,
        0  # Placeholder
    )
end

"""
    run_pomdp_rl(reads::Vector{String}; kwargs...)

Run the POMDPs.jl implementation and return unified history.
"""
function run_pomdp_rl(reads::Vector{String};
    solver::Symbol=:value_iteration,
    use_pomdp::Bool=false,
    training_episodes::Int=100,
    initial_k::Int=19,
    verbose::Bool=false
)
    start_time = Dates.now()
    
    # Train policy
    training_reads = reads[1:min(1000, length(reads))]
    policy, mdp = train_assembly_agent_pomdp(
        training_reads,
        solver=solver,
        use_pomdp=use_pomdp
    )
    
    # Apply to full dataset
    state = AssemblyState(
        initial_k,
        0.0, 0.0, 0.1, 0.0, 0.0, 0.0,
        Float64[],
        Int[initial_k],
        0, 0.0
    )
    
    k_values = Int[]
    quality_scores = Float64[]
    correction_rates = Float64[]
    memory_usage = Float64[]
    actions_taken = Symbol[]
    rewards = Float64[]
    time_per_step = Float64[]
    
    while !POMDPs.isterminal(mdp, state)
        step_start = Dates.now()
        
        action = POMDPs.action(policy, state)
        
        # Execute action
        new_state, reward_components = execute_assembly_action(
            state, action, mdp.max_k, mdp.memory_limit
        )
        
        # Record data
        push!(k_values, state.current_k)
        push!(quality_scores, new_state.assembly_quality)
        push!(correction_rates, new_state.correction_rate)
        push!(memory_usage, new_state.memory_usage)
        push!(actions_taken, action.decision)
        push!(rewards, calculate_total_reward(reward_components))
        
        step_time = Dates.value(Dates.now() - step_start) / 1000.0
        push!(time_per_step, step_time)
        
        state = new_state
        
        if verbose && length(k_values) % 10 == 0
            println("Step $(length(k_values)): K=$(state.current_k), " *
                   "Quality=$(round(state.assembly_quality, digits=3))")
        end
    end
    
    total_time = Dates.value(Dates.now() - start_time) / 1000.0
    
    return UnifiedAssemblyHistory(
        "POMDPs.jl ($solver)",
        k_values,
        quality_scores,
        correction_rates,
        memory_usage,
        actions_taken,
        rewards,
        time_per_step,
        total_time,
        length(quality_scores) > 0 ? quality_scores[end] : 0.0,
        0  # Placeholder
    )
end

"""
    run_mcts_rl(reads::Vector{String}; kwargs...)

Run the MCTS implementation and return unified history.
"""
function run_mcts_rl(reads::Vector{String}; 
    n_simulations::Int=1000,
    exploration_constant::Float64=sqrt(2),
    training_episodes::Int=100,
    verbose::Bool=false
)
    start_time = Dates.now()
    time_per_step = Float64[]
    
    # Train MCTS agent if requested
    if training_episodes > 0
        # Generate training data
        training_data = generate_diverse_training_data(5, 1000)
        
        # Train agent
        agent, _ = train_mcts_agent(
            training_data;
            n_episodes=training_episodes,
            n_simulations=n_simulations,
            exploration_constant=exploration_constant
        )
    else
        # Use untrained agent
        agent = MCTSAgent(
            n_simulations=n_simulations,
            exploration_constant=exploration_constant
        )
    end
    
    # Run assembly with MCTS
    env = AssemblyEnvironment(reads)
    state = get_state(env)
    
    k_values = [state.current_k]
    quality_scores = [state.assembly_quality]
    correction_rates = [state.correction_rate]
    memory_usage = [state.memory_usage]
    actions_taken = Symbol[]
    rewards = Float64[]
    
    while !is_done(env)
        step_start = Dates.now()
        
        # Select action using MCTS
        action = select_action_mcts(agent, state)
        push!(actions_taken, action.decision)
        
        # Take action
        new_state, reward, done, _ = step!(env, action)
        
        # Record metrics
        push!(k_values, new_state.current_k)
        push!(quality_scores, new_state.assembly_quality)
        push!(correction_rates, new_state.correction_rate)
        push!(memory_usage, new_state.memory_usage)
        push!(rewards, reward)
        
        step_time = Dates.value(Dates.now() - step_start) / 1000.0
        push!(time_per_step, step_time)
        
        state = new_state
        
        if verbose && length(k_values) % 10 == 0
            println("Step $(length(k_values)): K=$(state.current_k), " *
                   "Quality=$(round(state.assembly_quality, digits=3)), " *
                   "Simulations=$n_simulations")
        end
    end
    
    total_time = Dates.value(Dates.now() - start_time) / 1000.0
    
    return UnifiedAssemblyHistory(
        "MCTS (n_sim=$n_simulations)",
        k_values,
        quality_scores,
        correction_rates,
        memory_usage,
        actions_taken,
        rewards,
        time_per_step,
        total_time,
        length(quality_scores) > 0 ? quality_scores[end] : 0.0,
        0  # Placeholder
    )
end

# Comparison Metrics

"""
    calculate_comparison_metrics(results::Dict{Symbol, UnifiedAssemblyHistory})

Calculate metrics to compare RL approaches.
"""
function calculate_comparison_metrics(results::Dict{Symbol, UnifiedAssemblyHistory})
    metrics = Dict{String, Any}()
    
    # Final quality comparison
    metrics["final_quality"] = Dict(
        string(k) => v.final_quality for (k, v) in results
    )
    
    # Average quality throughout assembly
    metrics["average_quality"] = Dict(
        string(k) => Statistics.mean(v.quality_scores) for (k, v) in results
    )
    
    # Total time comparison
    metrics["total_time"] = Dict(
        string(k) => v.total_time for (k, v) in results
    )
    
    # Number of k-mer sizes tried
    metrics["k_values_tried"] = Dict(
        string(k) => length(unique(v.k_values)) for (k, v) in results
    )
    
    # Average reward
    metrics["average_reward"] = Dict(
        string(k) => Statistics.mean(v.rewards) for (k, v) in results
    )
    
    # Memory efficiency
    metrics["max_memory_usage"] = Dict(
        string(k) => maximum(v.memory_usage) for (k, v) in results
    )
    
    # Convergence speed (steps to reach 90% of final quality)
    metrics["convergence_speed"] = Dict{String, Any}()
    for (approach, history) in results
        target_quality = 0.9 * history.final_quality
        conv_idx = findfirst(q -> q >= target_quality, history.quality_scores)
        metrics["convergence_speed"][string(approach)] = 
            isnothing(conv_idx) ? length(history.quality_scores) : conv_idx
    end
    
    return metrics
end

# Visualization Functions

"""
    plot_comparison_results(comparison::RLComparison, output_dir::String)

Generate comparison plots for the RL approaches.
"""
function plot_comparison_results(comparison::RLComparison, output_dir::String)
    # Quality over time plot
    fig1 = Makie.Figure(resolution=(800, 600))
    ax1 = Makie.Axis(fig1[1, 1], 
        xlabel="Step", 
        ylabel="Assembly Quality",
        title="Assembly Quality Progression"
    )
    
    if !isnothing(comparison.custom_results)
        Makie.lines!(ax1, comparison.custom_results.quality_scores, 
                    label="Custom", linewidth=2)
    end
    if !isnothing(comparison.rljl_results)
        Makie.lines!(ax1, comparison.rljl_results.quality_scores, 
                    label="RL.jl", linewidth=2)
    end
    if !isnothing(comparison.pomdp_results)
        Makie.lines!(ax1, comparison.pomdp_results.quality_scores, 
                    label="POMDPs.jl", linewidth=2)
    end
    
    Makie.axislegend(ax1)
    Makie.save(joinpath(output_dir, "quality_comparison.png"), fig1)
    
    # K-mer progression plot
    fig2 = Makie.Figure(resolution=(800, 600))
    ax2 = Makie.Axis(fig2[1, 1], 
        xlabel="Step", 
        ylabel="K-mer Size",
        title="K-mer Size Progression"
    )
    
    if !isnothing(comparison.custom_results)
        Makie.stairs!(ax2, comparison.custom_results.k_values, 
                     label="Custom", linewidth=2)
    end
    if !isnothing(comparison.rljl_results)
        Makie.stairs!(ax2, comparison.rljl_results.k_values, 
                     label="RL.jl", linewidth=2)
    end
    if !isnothing(comparison.pomdp_results)
        Makie.stairs!(ax2, comparison.pomdp_results.k_values, 
                     label="POMDPs.jl", linewidth=2)
    end
    
    Makie.axislegend(ax2)
    Makie.save(joinpath(output_dir, "kmer_progression.png"), fig2)
    
    # Performance metrics bar plot
    fig3 = Makie.Figure(resolution=(1000, 600))
    
    # Final quality
    ax3_1 = Makie.Axis(fig3[1, 1], 
        ylabel="Final Quality",
        title="Final Assembly Quality"
    )
    
    approaches = String[]
    final_qualities = Float64[]
    
    for (k, v) in comparison.comparison_metrics["final_quality"]
        push!(approaches, k)
        push!(final_qualities, v)
    end
    
    Makie.barplot!(ax3_1, 1:length(approaches), final_qualities)
    ax3_1.xticks = (1:length(approaches), approaches)
    
    # Time comparison
    ax3_2 = Makie.Axis(fig3[1, 2], 
        ylabel="Time (seconds)",
        title="Total Runtime"
    )
    
    times = [comparison.comparison_metrics["total_time"][k] for k in approaches]
    Makie.barplot!(ax3_2, 1:length(approaches), times)
    ax3_2.xticks = (1:length(approaches), approaches)
    
    Makie.save(joinpath(output_dir, "performance_metrics.png"), fig3)
end

"""
    save_comparison_results(comparison::RLComparison, output_dir::String)

Save comparison results to disk for later analysis.
"""
function save_comparison_results(comparison::RLComparison, output_dir::String)
    # Save as JLD2
    jld2_path = joinpath(output_dir, "$(comparison.dataset_name)_comparison.jld2")
    JLD2.save(jld2_path, "comparison", comparison)
    
    # Save metrics as JSON
    json_path = joinpath(output_dir, "$(comparison.dataset_name)_metrics.json")
    open(json_path, "w") do f
        JSON.print(f, comparison.comparison_metrics, 4)
    end
    
    # Save summary as CSV
    csv_path = joinpath(output_dir, "$(comparison.dataset_name)_summary.csv")
    
    summary_data = []
    for (approach, metrics) in comparison.comparison_metrics["final_quality"]
        push!(summary_data, (
            approach=approach,
            final_quality=metrics,
            avg_quality=comparison.comparison_metrics["average_quality"][approach],
            total_time=comparison.comparison_metrics["total_time"][approach],
            k_values_tried=comparison.comparison_metrics["k_values_tried"][approach],
            avg_reward=comparison.comparison_metrics["average_reward"][approach],
            max_memory=comparison.comparison_metrics["max_memory_usage"][approach],
            convergence_speed=comparison.comparison_metrics["convergence_speed"][approach]
        ))
    end
    
    df = DataFrames.DataFrame(summary_data)
    CSV.write(csv_path, df)
    
    println("Results saved to $output_dir")
end

# Benchmark Functions

"""
    benchmark_rl_approaches(;
        n_datasets::Int=5,
        dataset_sizes::Vector{Int}=[100, 500, 1000],
        kwargs...
    )

Run comprehensive benchmarks comparing RL approaches.
"""
function benchmark_rl_approaches(;
    n_datasets::Int=5,
    dataset_sizes::Vector{Int}=[100, 500, 1000],
    approaches::Vector{Symbol}=[:custom, :rljl, :pomdp],
    output_dir::String="rl_benchmarks",
    kwargs...
)
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    all_results = []
    
    for size in dataset_sizes
        println("\nBenchmarking with dataset size: $size")
        
        for i in 1:n_datasets
            println("Dataset $i/$n_datasets")
            
            # Generate synthetic dataset
            reads = generate_synthetic_reads(size; complexity=:medium)
            
            # Run comparison
            comparison = compare_rl_approaches(
                reads,
                approaches=approaches,
                dataset_name="benchmark_$(size)_$(i)",
                verbose=false,
                save_results=false,
                kwargs...
            )
            
            push!(all_results, (size=size, dataset=i, comparison=comparison))
        end
    end
    
    # Aggregate results
    aggregate_benchmark_results(all_results, output_dir)
    
    return all_results
end

"""
    aggregate_benchmark_results(results::Vector, output_dir::String)

Aggregate and visualize benchmark results.
"""
function aggregate_benchmark_results(results::Vector, output_dir::String)
    # Prepare data for analysis
    benchmark_data = []
    
    for (size, dataset, comparison) in results
        for (approach, final_quality) in comparison.comparison_metrics["final_quality"]
            push!(benchmark_data, (
                dataset_size=size,
                dataset_id=dataset,
                approach=approach,
                final_quality=final_quality,
                total_time=comparison.comparison_metrics["total_time"][approach],
                avg_reward=comparison.comparison_metrics["average_reward"][approach]
            ))
        end
    end
    
    df = DataFrames.DataFrame(benchmark_data)
    
    # Create summary statistics
    summary = DataFrames.combine(
        DataFrames.groupby(df, [:dataset_size, :approach]),
        :final_quality => Statistics.mean => :mean_quality,
        :final_quality => Statistics.std => :std_quality,
        :total_time => Statistics.mean => :mean_time,
        :total_time => Statistics.std => :std_time,
        :avg_reward => Statistics.mean => :mean_reward
    )
    
    # Save results
    CSV.write(joinpath(output_dir, "benchmark_summary.csv"), summary)
    
    # Create visualization
    fig = Makie.Figure(resolution=(1200, 800))
    
    # Quality vs dataset size
    ax1 = Makie.Axis(fig[1, 1], 
        xlabel="Dataset Size", 
        ylabel="Mean Final Quality",
        title="Assembly Quality by Dataset Size"
    )
    
    for approach in unique(summary.approach)
        subset = summary[summary.approach .== approach, :]
        Makie.scatterlines!(ax1, subset.dataset_size, subset.mean_quality,
                           label=approach, markersize=10)
    end
    
    Makie.axislegend(ax1)
    
    # Time vs dataset size
    ax2 = Makie.Axis(fig[2, 1], 
        xlabel="Dataset Size", 
        ylabel="Mean Time (seconds)",
        title="Runtime by Dataset Size"
    )
    
    for approach in unique(summary.approach)
        subset = summary[summary.approach .== approach, :]
        Makie.scatterlines!(ax2, subset.dataset_size, subset.mean_time,
                           label=approach, markersize=10)
    end
    
    Makie.axislegend(ax2)
    
    Makie.save(joinpath(output_dir, "benchmark_results.png"), fig)
    
    println("\nBenchmark Summary:")
    println(summary)
end