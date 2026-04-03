import Pkg
if isinteractive()
    Pkg.activate("..")
end

import CSV
import DataFrames
import Dates
import Distributions
import Mycelia
import Plots
import StableRNGs
import Statistics

include("benchmark_utils.jl")

struct SyntheticRepeatScenario
    repeat_length::Int
    divergence::Float64
    coverage::Int
end

function scenario_name(scenario::SyntheticRepeatScenario)
    divergence_label = lpad(string(round(Int, scenario.divergence * 1000)), 3, '0')
    return "repeat$(scenario.repeat_length)_div$(divergence_label)_cov$(scenario.coverage)"
end

function build_scenarios(config::Dict{String, Any})
    scenarios = SyntheticRepeatScenario[]
    for repeat_length in config["repeat_lengths"]
        for divergence in config["divergences"]
            for coverage in config["coverages"]
                push!(scenarios, SyntheticRepeatScenario(repeat_length, divergence, coverage))
            end
        end
    end
    return scenarios
end

function threshold_label(threshold::Float64)
    rounded = round(threshold; digits = 3)
    return replace(string(rounded), "." => "p")
end

function weight_mode_name(weight_mode::Symbol)
    return weight_mode == :llr ? "llr_mapping_quality" : String(weight_mode)
end

function mode_seed_offset(weight_mode::Symbol)
    if weight_mode == :uniform
        return 0
    elseif weight_mode == :linear
        return 10_000
    elseif weight_mode == :saturating
        return 20_000
    elseif weight_mode == :llr
        return 30_000
    end
    error("Unsupported weight mode: $(weight_mode)")
end

function scenario_flank_length(scenario::SyntheticRepeatScenario)
    return max(80, scenario.repeat_length ÷ 2)
end

function scenario_reference_span(scenario::SyntheticRepeatScenario)
    return scenario.repeat_length + 2 * scenario_flank_length(scenario)
end

function trial_contiguity_span(
        scenario::SyntheticRepeatScenario;
        correct::Bool,
        misassembly::Bool)
    flank_length = scenario_flank_length(scenario)
    if correct
        return scenario_reference_span(scenario)
    elseif misassembly
        return max(40, flank_length ÷ 2)
    end
    return flank_length
end

function synthetic_nga50(length_pairs::Vector{Tuple{Int, Int}})
    isempty(length_pairs) && return 0
    sorted_pairs = sort(copy(length_pairs); by = first, rev = true)
    half_reference = sum(last(pair) for pair in sorted_pairs) / 2
    cumulative_reference = 0.0

    for (aligned_span, reference_span) in sorted_pairs
        cumulative_reference += reference_span
        if cumulative_reference >= half_reference
            return aligned_span
        end
    end

    return sorted_pairs[end][1]
end

function rank_metric(value)
    return ismissing(value) ? -Inf : Float64(value)
end

function plot_metric(value)
    return ismissing(value) ? NaN : Float64(value)
end

function format_metric(value; digits::Int = 3)
    return ismissing(value) ? "n/a" : string(round(Float64(value); digits = digits))
end

function informative_window_count(scenario::SyntheticRepeatScenario, event_count::Int)
    if scenario.divergence >= 0.18
        return min(3, max(1, event_count ÷ 2))
    elseif scenario.divergence >= 0.08
        return min(2, max(1, event_count ÷ 2))
    end
    return 1
end

function simulate_repeat_windows(
        rng::StableRNGs.LehmerRNG,
        scenario::SyntheticRepeatScenario,
        true_branch::Symbol)
    event_count = max(4, cld(scenario.repeat_length, 60) + 2)
    flank_windows = informative_window_count(scenario, event_count)
    base_depth = max(3, ceil(Int, scenario.coverage / 2))
    signal_fraction = clamp(0.5 + 0.45 * sqrt(scenario.divergence), 0.5, 0.95)
    ambiguous_fraction = clamp(0.5 + 0.05 * scenario.divergence, 0.5, 0.6)

    windows = Tuple{Int, Int}[]
    for event_idx in 1:event_count
        event_depth = max(2, base_depth + rand(rng, -2:2))
        is_flank = event_idx <= flank_windows || event_idx > event_count - flank_windows
        favored_fraction = is_flank ? signal_fraction : ambiguous_fraction
        favored_support = rand(rng, Distributions.Binomial(event_depth, favored_fraction))
        competing_support = event_depth - favored_support

        if true_branch == :reference
            push!(windows, (favored_support, competing_support))
        else
            push!(windows, (competing_support, favored_support))
        end
    end

    return windows
end

function observe_weighted_window!(
        state,
        reference_support::Real,
        alternative_support::Real,
        weight_mode::Symbol,
        scenario::SyntheticRepeatScenario;
        alpha::Float64,
        beta::Float64)
    if weight_mode == :llr
        return Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state,
            :reference,
            :alternative,
            reference_support,
            alternative_support;
            alpha = alpha,
            beta = beta
        )
    end

    if weight_mode == :uniform
        weighted_reference = Mycelia.Rhizomorph.MomentumForkResolver.uniform_weight(reference_support)
        weighted_alternative = Mycelia.Rhizomorph.MomentumForkResolver.uniform_weight(alternative_support)
    elseif weight_mode == :linear
        weighted_reference = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight(reference_support)
        weighted_alternative = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight(alternative_support)
    elseif weight_mode == :saturating
        max_weight = max(2.0, scenario.coverage / 3)
        half_saturation = max(1.0, scenario.coverage / 5)
        weighted_reference = Mycelia.Rhizomorph.MomentumForkResolver.saturating_weight(
            reference_support; max_weight = max_weight, half_saturation = half_saturation)
        weighted_alternative = Mycelia.Rhizomorph.MomentumForkResolver.saturating_weight(
            alternative_support; max_weight = max_weight, half_saturation = half_saturation)
    else
        error("Unsupported weight mode: $(weight_mode)")
    end

    state.branch_support[:reference] += Float64(reference_support)
    state.branch_support[:alternative] += Float64(alternative_support)
    state.branch_scores[:reference] += weighted_reference
    state.branch_scores[:alternative] += weighted_alternative
    state.observations += 1
    Mycelia.Rhizomorph.MomentumForkResolver._refresh_state!(state; alpha = alpha, beta = beta)
    return state.decision
end

function run_resolver_trial(
        rng::StableRNGs.LehmerRNG,
        scenario::SyntheticRepeatScenario,
        true_branch::Symbol,
        weight_mode::Symbol;
        alpha::Float64,
        beta::Float64)
    state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
        "synthetic-read",
        [:reference, :alternative];
        reference_branch = :reference
    )

    for (reference_support, alternative_support) in
        simulate_repeat_windows(rng, scenario, true_branch)
        decision = observe_weighted_window!(
            state,
            reference_support,
            alternative_support,
            weight_mode,
            scenario;
            alpha = alpha,
            beta = beta
        )
        if decision != :continue
            break
        end
    end

    resolved = !isnothing(state.resolved_branch)
    correct = resolved && state.resolved_branch == true_branch
    misassembly = resolved && state.resolved_branch != true_branch
    reference_span = scenario_reference_span(scenario)
    contiguity_span = trial_contiguity_span(scenario; correct = correct, misassembly = misassembly)

    return (
        resolved = resolved,
        correct = correct,
        misassembly = misassembly,
        observations = state.observations,
        final_llr = state.log_likelihood_ratio,
        reference_span = reference_span,
        contiguity_span = contiguity_span
    )
end

function scenario_truth_branch(scenario_index::Int, replicate::Int)
    return isodd(scenario_index + replicate) ? :reference : :alternative
end

function summarize_trials(
        scenario::SyntheticRepeatScenario,
        weight_mode::Symbol,
        alpha::Float64,
        beta::Float64,
        trial_rows::Vector)
    total_trials = length(trial_rows)
    resolved_trials = count(row -> row.resolved, trial_rows)
    misassembly_count = count(row -> row.misassembly, trial_rows)
    correct_count = count(row -> row.correct, trial_rows)
    unresolved_trials = total_trials - resolved_trials

    resolved_accuracy = resolved_trials == 0 ? missing : correct_count / resolved_trials
    length_pairs = [(row.contiguity_span, row.reference_span) for row in trial_rows]
    synthetic_nga50_bp = synthetic_nga50(length_pairs)
    mean_reference_span = Statistics.mean(row.reference_span for row in trial_rows)

    return (
        scenario = scenario_name(scenario),
        repeat_length = scenario.repeat_length,
        divergence = scenario.divergence,
        coverage = scenario.coverage,
        weight_mode = weight_mode_name(weight_mode),
        alpha = alpha,
        beta = beta,
        total_trials = total_trials,
        resolved_trials = resolved_trials,
        unresolved_trials = unresolved_trials,
        correct_trials = correct_count,
        misassembly_count = misassembly_count,
        resolution_rate = resolved_trials / total_trials,
        fork_accuracy = resolved_accuracy,
        correct_resolution_rate = correct_count / total_trials,
        misassembly_rate = misassembly_count / total_trials,
        synthetic_nga50_bp = synthetic_nga50_bp,
        synthetic_nga50_fraction = synthetic_nga50_bp / mean_reference_span,
        mean_reference_span = mean_reference_span,
        mean_observations = Statistics.mean(row.observations for row in trial_rows),
        mean_final_llr = Statistics.mean(row.final_llr for row in trial_rows)
    )
end

function aggregate_threshold_rows(
        scenario_rows::Vector,
        length_pairs::Vector{Tuple{Int, Int}},
        weight_mode::Symbol,
        alpha::Float64,
        beta::Float64)
    total_trials = sum(row.total_trials for row in scenario_rows)
    resolved_trials = sum(row.resolved_trials for row in scenario_rows)
    unresolved_trials = sum(row.unresolved_trials for row in scenario_rows)
    correct_trials = sum(row.correct_trials for row in scenario_rows)
    misassembly_count = sum(row.misassembly_count for row in scenario_rows)
    fork_accuracy = resolved_trials == 0 ? missing : correct_trials / resolved_trials
    synthetic_nga50_bp = synthetic_nga50(length_pairs)
    mean_reference_span = Statistics.mean(last(pair) for pair in length_pairs)

    return (
        weight_mode = weight_mode_name(weight_mode),
        alpha = alpha,
        beta = beta,
        scenario_count = length(scenario_rows),
        total_trials = total_trials,
        resolved_trials = resolved_trials,
        unresolved_trials = unresolved_trials,
        correct_trials = correct_trials,
        misassembly_count = misassembly_count,
        resolution_rate = resolved_trials / total_trials,
        fork_accuracy = fork_accuracy,
        correct_resolution_rate = correct_trials / total_trials,
        misassembly_rate = misassembly_count / total_trials,
        synthetic_nga50_bp = synthetic_nga50_bp,
        synthetic_nga50_fraction = synthetic_nga50_bp / mean_reference_span,
        mean_reference_span = mean_reference_span,
        mean_observations = Statistics.mean(row.mean_observations for row in scenario_rows),
        mean_final_llr = Statistics.mean(row.mean_final_llr for row in scenario_rows)
    )
end

function run_threshold_batch(config::Dict{String, Any}, weight_mode::Symbol, alpha::Float64, beta::Float64)
    scenario_rows = NamedTuple[]
    scenarios = build_scenarios(config)
    length_pairs = Tuple{Int, Int}[]
    mode_offset = mode_seed_offset(weight_mode)
    threshold_offset = round(Int, alpha * 100_000)

    for (scenario_index, scenario) in enumerate(scenarios)
        rng = StableRNGs.StableRNG(config["seed"] + scenario_index + mode_offset + threshold_offset)
        trial_rows = NamedTuple[]

        for replicate in 1:config["replicates"]
            true_branch = scenario_truth_branch(scenario_index, replicate)
            push!(
                trial_rows,
                run_resolver_trial(
                    rng,
                    scenario,
                    true_branch,
                    weight_mode;
                    alpha = alpha,
                    beta = beta
                )
            )
        end

        append!(length_pairs, [(row.contiguity_span, row.reference_span) for row in trial_rows])
        push!(scenario_rows, summarize_trials(scenario, weight_mode, alpha, beta, trial_rows))
    end

    threshold_row = aggregate_threshold_rows(scenario_rows, length_pairs, weight_mode, alpha, beta)
    return (
        scenario_rows = scenario_rows,
        threshold_row = threshold_row
    )
end

function comparison_score(row)
    return (
        Float64(row.correct_resolution_rate),
        -Float64(row.misassembly_rate),
        Float64(row.synthetic_nga50_fraction),
        rank_metric(row.fork_accuracy),
        -Float64(row.mean_observations)
    )
end

function select_best_threshold_rows(threshold_df::DataFrames.DataFrame)
    comparison_rows = NamedTuple[]
    for weight_mode in unique(threshold_df[!, :weight_mode])
        mode_df = threshold_df[threshold_df[!, :weight_mode] .== weight_mode, :]
        rows = collect(DataFrames.eachrow(mode_df))
        best_row = NamedTuple(rows[1])
        best_score = comparison_score(best_row)

        for row in rows[2:end]
            candidate = NamedTuple(row)
            candidate_score = comparison_score(candidate)
            if candidate_score > best_score
                best_row = candidate
                best_score = candidate_score
            end
        end

        push!(
            comparison_rows,
            (
                weight_mode = best_row.weight_mode,
                alpha = best_row.alpha,
                beta = best_row.beta,
                resolution_rate = best_row.resolution_rate,
                fork_accuracy = best_row.fork_accuracy,
                correct_resolution_rate = best_row.correct_resolution_rate,
                misassembly_rate = best_row.misassembly_rate,
                synthetic_nga50_bp = best_row.synthetic_nga50_bp,
                synthetic_nga50_fraction = best_row.synthetic_nga50_fraction,
                mean_observations = best_row.mean_observations
            )
        )
    end

    comparison_df = DataFrames.DataFrame(comparison_rows)
    DataFrames.sort!(comparison_df, :weight_mode)
    return comparison_df
end

function write_comparison_markdown(
        comparison_df::DataFrames.DataFrame,
        output_path::AbstractString)
    open(output_path, "w") do io
        println(io, "# Momentum Fork-Resolution Weight Ablation")
        println(io)
        println(io, "Synthetic repeat benchmark comparison for uniform, linear, saturating, and LLR/mapping-quality weighting.")
        println(io)
        println(io, "Synthetic NGA50 is a repeat-resolution contiguity proxy derived from the benchmark's aligned-span model.")
        println(io)
        println(io, "| Weight mode | alpha=beta | Fork accuracy | Correct resolution rate | Misassembly rate | Synthetic NGA50 (bp) | Synthetic NGA50 / reference | Resolution rate | Mean observations |")
        println(io, "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")

        for row in DataFrames.eachrow(comparison_df)
            println(
                io,
                "| $(row.weight_mode) | $(format_metric(row.alpha; digits = 3)) | " *
                "$(format_metric(row.fork_accuracy; digits = 3)) | " *
                "$(format_metric(row.correct_resolution_rate; digits = 3)) | " *
                "$(format_metric(row.misassembly_rate; digits = 3)) | " *
                "$(format_metric(row.synthetic_nga50_bp; digits = 1)) | " *
                "$(format_metric(row.synthetic_nga50_fraction; digits = 3)) | " *
                "$(format_metric(row.resolution_rate; digits = 3)) | " *
                "$(format_metric(row.mean_observations; digits = 2)) |"
            )
        end
    end
    return output_path
end

function render_threshold_comparison_plot(
        threshold_df::DataFrames.DataFrame,
        output_path::AbstractString)
    palette = Dict(
        "uniform" => :gray40,
        "linear" => :steelblue,
        "saturating" => :darkorange,
        "llr_mapping_quality" => :seagreen
    )

    metric_specs = [
        (:fork_accuracy, "Fork accuracy"),
        (:misassembly_rate, "Misassembly rate"),
        (:resolution_rate, "Resolution rate"),
        (:synthetic_nga50_fraction, "Synthetic NGA50 / reference")
    ]

    panels = Any[]
    for (metric, title) in metric_specs
        panel = Plots.plot(
            title = title,
            xlabel = "alpha = beta error rate",
            ylabel = title,
            legend = :bottomright,
            linewidth = 2,
            marker = :circle
        )
        for weight_mode in sort(unique(threshold_df[!, :weight_mode]))
            mode_df = DataFrames.sort(
                threshold_df[threshold_df[!, :weight_mode] .== weight_mode, :],
                :alpha
            )
            Plots.plot!(
                panel,
                Float64.(mode_df[!, :alpha]),
                [plot_metric(value) for value in mode_df[!, metric]];
                label = weight_mode,
                color = get(palette, weight_mode, :black)
            )
        end
        push!(panels, panel)
    end

    comparison_plot = Plots.plot(panels...; layout = (2, 2), size = (1200, 800))
    Plots.savefig(comparison_plot, output_path)
    return output_path
end

function benchmark_config()
    small_config = Dict{String, Any}(
        "repeat_lengths" => [120, 240, 480],
        "divergences" => [0.02, 0.08, 0.16],
        "coverages" => [12, 24, 48],
        "replicates" => 24,
        "error_rates" => [0.2, 0.1, 0.05, 0.01],
        "weight_modes" => [:uniform, :linear, :saturating, :llr],
        "benchmark_samples" => 3,
        "benchmark_seconds" => 3,
        "seed" => 20260403,
        "description" => "Small scale - repeat fork calibration"
    )

    medium_config = Dict{String, Any}(
        "repeat_lengths" => [120, 240, 480, 720],
        "divergences" => [0.01, 0.04, 0.08, 0.16],
        "coverages" => [8, 16, 32, 64],
        "replicates" => 48,
        "error_rates" => [0.2, 0.1, 0.05, 0.01, 0.005],
        "weight_modes" => [:uniform, :linear, :saturating, :llr],
        "benchmark_samples" => 3,
        "benchmark_seconds" => 5,
        "seed" => 20260403,
        "description" => "Medium scale - balanced tradeoff sweep"
    )

    large_config = Dict{String, Any}(
        "repeat_lengths" => [120, 240, 480, 720, 960],
        "divergences" => [0.005, 0.01, 0.04, 0.08, 0.16],
        "coverages" => [8, 16, 32, 64, 96],
        "replicates" => 96,
        "error_rates" => [0.2, 0.1, 0.05, 0.01, 0.005],
        "weight_modes" => [:uniform, :linear, :saturating, :llr],
        "benchmark_samples" => 4,
        "benchmark_seconds" => 8,
        "seed" => 20260403,
        "description" => "Large scale - dense repeat stress sweep"
    )

    scale = get(ENV, "BENCHMARK_SCALE", "small")
    if scale == "medium"
        return medium_config
    elseif scale == "large"
        return large_config
    end
    return small_config
end

println("=== Momentum Fork Resolution Benchmark ===")
println("Start time: $(Dates.now())")

config = benchmark_config()
println("Configuration: $(config["description"])")
println("Repeat lengths: $(config["repeat_lengths"])")
println("Divergences: $(config["divergences"])")
println("Coverages: $(config["coverages"])")
println("Error-rate sweep: $(config["error_rates"])")
println("Weight modes: $(config["weight_modes"])")
println("Replicates per scenario: $(config["replicates"])")

benchmark_suite = BenchmarkSuite("Momentum Fork Resolution Synthetic Repeat Benchmark")
scenario_rows = NamedTuple[]
threshold_rows = NamedTuple[]

for weight_mode in config["weight_modes"]
    println("\n--- Weight mode: $(weight_mode) ---")
    for error_rate in config["error_rates"]
        alpha = Float64(error_rate)
        beta = Float64(error_rate)
        batch = run_threshold_batch(config, weight_mode, alpha, beta)
        append!(scenario_rows, batch.scenario_rows)
        push!(threshold_rows, batch.threshold_row)

        benchmark_result, memory_stats = run_benchmark_with_memory(
            run_threshold_batch,
            config,
            weight_mode,
            alpha,
            beta;
            samples = config["benchmark_samples"],
            seconds = config["benchmark_seconds"]
        )

        result_name = "$(weight_mode_name(weight_mode))_alpha$(threshold_label(alpha))"
        add_benchmark_result!(benchmark_suite, result_name, benchmark_result, memory_stats)
        benchmark_suite.results[result_name]["summary"] = Dict(
            "resolution_rate" => batch.threshold_row.resolution_rate,
            "fork_accuracy" => ismissing(batch.threshold_row.fork_accuracy) ? 0.0 :
                               batch.threshold_row.fork_accuracy,
            "correct_resolution_rate" => batch.threshold_row.correct_resolution_rate,
            "misassembly_count" => batch.threshold_row.misassembly_count,
            "synthetic_nga50_bp" => batch.threshold_row.synthetic_nga50_bp,
            "resolved_trials" => batch.threshold_row.resolved_trials,
            "total_trials" => batch.threshold_row.total_trials
        )

        accuracy_label = ismissing(batch.threshold_row.fork_accuracy) ? "n/a" :
                         string(round(batch.threshold_row.fork_accuracy; digits = 3))
        println(
            "alpha=beta=$(alpha): resolution=$(round(batch.threshold_row.resolution_rate; digits=3)), " *
            "accuracy=$(accuracy_label), misassemblies=$(batch.threshold_row.misassembly_count)"
        )
    end
end

scenario_df = DataFrames.DataFrame(scenario_rows)
threshold_df = DataFrames.DataFrame(threshold_rows)
DataFrames.sort!(scenario_df, [:weight_mode, :alpha, :repeat_length, :divergence, :coverage])
DataFrames.sort!(threshold_df, [:weight_mode, :alpha])
comparison_df = select_best_threshold_rows(threshold_df)

results_dir = mkpath(joinpath(@__DIR__, "results"))
timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
results_prefix = "08_momentum_fork_resolution_benchmark"
scenario_csv = joinpath(results_dir, "$(results_prefix)_scenarios_$(timestamp).csv")
threshold_csv = joinpath(results_dir, "$(results_prefix)_thresholds_$(timestamp).csv")
comparison_csv = joinpath(results_dir, "$(results_prefix)_comparison_$(timestamp).csv")
comparison_md = joinpath(results_dir, "$(results_prefix)_comparison_$(timestamp).md")
comparison_plot = joinpath(results_dir, "$(results_prefix)_comparison_$(timestamp).png")
results_json = joinpath(results_dir, "$(results_prefix)_$(timestamp).json")

CSV.write(scenario_csv, scenario_df)
CSV.write(threshold_csv, threshold_df)
CSV.write(comparison_csv, comparison_df)
write_comparison_markdown(comparison_df, comparison_md)
render_threshold_comparison_plot(threshold_df, comparison_plot)
save_benchmark_results(benchmark_suite, results_json)
format_benchmark_summary(benchmark_suite)

println("\nThreshold tradeoff summary:")
println(threshold_df)
println("\nBest per-mode comparison summary:")
println(comparison_df)
println("\nScenario results saved to: $(scenario_csv)")
println("Threshold results saved to: $(threshold_csv)")
println("Comparison table saved to: $(comparison_csv)")
println("Comparison markdown saved to: $(comparison_md)")
println("Comparison plot saved to: $(comparison_plot)")
println("Benchmark summary saved to: $(results_json)")
println("End time: $(Dates.now())")
