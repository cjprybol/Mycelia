#!/usr/bin/env julia

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import Dates
import DataFrames
import CSV
import JSON
import Statistics

include("rhizomorph_benchmark_suite_common.jl")

println("=== Rhizomorph Reconstruction Suite (Isolate) ===")
println("Start time: $(Dates.now())")

function _parse_int_list_env(var_name::String, default_values::Vector{Int})
    value = strip(get(ENV, var_name, ""))
    if isempty(value)
        return default_values
    end
    return [parse(Int, strip(token)) for token in split(value, ",") if !isempty(strip(token))]
end

function _to_json_row_table(dataframe::DataFrames.DataFrame)
    column_names = DataFrames.names(dataframe)
    rows = Vector{Dict{String, Any}}(undef, DataFrames.nrow(dataframe))
    for (index, row) in enumerate(DataFrames.eachrow(dataframe))
        row_dict = Dict{String, Any}()
        for column_name in column_names
            value = row[column_name]
            row_dict[string(column_name)] = ismissing(value) ? nothing : value
        end
        rows[index] = row_dict
    end
    return rows
end

function _scale_defaults(scale::String)
    if scale == "large"
        return (
            q_values = ISOLATE_Q_VALUES,
            coverage_levels = COVERAGE_LEVELS,
            replicates = 3,
            strategy = :tiered,
            anchor_q_values = ISOLATE_ANCHOR_Q_VALUES,
            anchor_coverage_levels = ISOLATE_ANCHOR_COVERAGE_LEVELS
        )
    elseif scale == "medium"
        return (
            q_values = ISOLATE_Q_VALUES,
            coverage_levels = [10, 25, 50, 100, 250],
            replicates = 2,
            strategy = :tiered,
            anchor_q_values = ISOLATE_ANCHOR_Q_VALUES,
            anchor_coverage_levels = [25, 100, 250]
        )
    end
    return (
        q_values = [10, 30, 50],
        coverage_levels = [10, 50, 250],
        replicates = 1,
        strategy = :tiered,
        anchor_q_values = [30],
        anchor_coverage_levels = [50]
    )
end

function _run_assembly(input_data, condition; k::Int = 31, min_overlap::Int = 31)
    graph_family = Symbol(condition.graph_family)
    memory_profile = Symbol(condition.memory_profile)
    tokenizer = Symbol(condition.tokenizer)

    if graph_family == :olc
        return Mycelia.Rhizomorph.assemble_genome(
            input_data;
            graph_family = graph_family,
            min_overlap = min_overlap,
            memory_profile = memory_profile
        )
    elseif graph_family == :token
        return Mycelia.Rhizomorph.assemble_genome(
            input_data;
            graph_family = graph_family,
            k = k,
            memory_profile = memory_profile,
            tokenizer = tokenizer
        )
    end

    return Mycelia.Rhizomorph.assemble_genome(
        input_data;
        graph_family = graph_family,
        k = k,
        memory_profile = memory_profile
    )
end

scale = lowercase(get(ENV, "BENCHMARK_SCALE", "small"))
defaults = _scale_defaults(scale)
strategy = Symbol(get(ENV, "MYCELIA_BENCHMARK_MATRIX_STRATEGY", string(defaults.strategy)))
default_max_reads = if scale == "large"
    100_000
elseif scale == "medium"
    50_000
else
    20_000
end

q_values = _parse_int_list_env("MYCELIA_BENCHMARK_Q_VALUES", defaults.q_values)
coverage_levels = _parse_int_list_env("MYCELIA_BENCHMARK_COVERAGES", defaults.coverage_levels)
anchor_q_values = _parse_int_list_env("MYCELIA_BENCHMARK_ANCHOR_Q_VALUES",
    defaults.anchor_q_values)
anchor_coverage_levels = _parse_int_list_env("MYCELIA_BENCHMARK_ANCHOR_COVERAGES",
    defaults.anchor_coverage_levels)
replicates = parse(Int, get(ENV, "MYCELIA_BENCHMARK_REPLICATES", string(defaults.replicates)))
include_advanced = lowercase(get(ENV, "MYCELIA_BENCHMARK_ENABLE_ADVANCED", "false")) ==
                   "true"
max_runs = parse(Int, get(ENV, "MYCELIA_BENCHMARK_MAX_RUNS", "0"))
max_reads_per_condition = parse(Int,
    get(ENV, "MYCELIA_BENCHMARK_MAX_READS_PER_CONDITION", string(default_max_reads)))

if strategy == :full_factorial
    anchor_q_values = q_values
    anchor_coverage_levels = coverage_levels
end

genome_specs = default_isolate_genome_specs(scale = scale)
genome_by_id = Dict(spec.genome_id => spec for spec in genome_specs)
genomes = synthesize_isolate_genomes(genome_specs)

conditions = build_isolate_benchmark_matrix(
    genomes = genome_specs,
    q_values = q_values,
    coverage_levels = coverage_levels,
    replicates = replicates,
    strategy = strategy,
    include_advanced = include_advanced,
    anchor_q_values = anchor_q_values,
    anchor_coverage_levels = anchor_coverage_levels
)

if max_runs > 0 && length(conditions) > max_runs
    conditions = conditions[1:max_runs]
end

println("Scale: $(scale)")
println("Strategy: $(strategy)")
println("Conditions: $(length(conditions))")

rows = NamedTuple[]
bundle_cache = Dict{Tuple{String, Int, Int, Int}, Any}()

for (index, condition) in enumerate(conditions)
    println("[$(index)/$(length(conditions))] $(condition.algorithm_id) | $(condition.genome_id) | Q$(condition.q_value) | cov=$(condition.coverage)x | rep=$(condition.replicate)")
    genome_spec = genome_by_id[condition.genome_id]
    genome_sequence = genomes[condition.genome_id]
    bundle_key = (
        condition.genome_id,
        condition.q_value,
        condition.coverage,
        condition.replicate
    )

    if !haskey(bundle_cache, bundle_key)
        simulation_seed = genome_spec.seed +
                          (condition.replicate * 10_000) +
                          (condition.coverage * 10) +
                          condition.q_value
        bundle_cache[bundle_key] = generate_isolate_read_bundle(
            reference_sequence = genome_sequence,
            coverage = condition.coverage,
            q_value = condition.q_value,
            seed = simulation_seed,
            max_records = max_reads_per_condition
        )
    end

    bundle = bundle_cache[bundle_key]
    benchmark_input = select_benchmark_input(bundle, condition.input_mode)
    reference_sequences = [String(genome_sequence)]

    try
        timing = @timed begin
            assembly = _run_assembly(benchmark_input, condition)
            if condition.stage == :basic_correction
                assembly = Mycelia.Rhizomorph.polish_assembly(assembly, bundle.fastq; iterations = 1)
            elseif condition.stage == :advanced_probabilistic
                assembly = Mycelia.Rhizomorph.polish_assembly(assembly, bundle.fastq; iterations = 3)
            end
            assembly
        end

        assembly = timing.value
        metrics = compute_assembly_metrics(reference_sequences, assembly.contigs; k = 31)
        assembly_stats = assembly.assembly_stats
        push!(rows, (
            benchmark = "rhizomorph_reconstruction_suite",
            dataset_type = "isolate",
            scale = scale,
            strategy = string(strategy),
            genome_id = condition.genome_id,
            genome_tier = condition.genome_tier,
            genome_length = condition.genome_length,
            q_value = condition.q_value,
            error_rate = bundle.error_rate,
            coverage = condition.coverage,
            replicate = condition.replicate,
            algorithm_id = condition.algorithm_id,
            graph_family = string(condition.graph_family),
            memory_profile = string(condition.memory_profile),
            tokenizer = string(condition.tokenizer),
            input_mode = string(condition.input_mode),
            stage = string(condition.stage),
            num_reads = bundle.n_reads,
            run_seconds = timing.time,
            peak_memory_bytes = Float64(timing.bytes),
            num_contigs = metrics.num_contigs,
            total_contig_length = metrics.total_contig_length,
            assembly_size_ratio = metrics.assembly_size_ratio,
            n50 = metrics.n50,
            l50 = metrics.l50,
            kmer_precision = metrics.kmer_precision,
            kmer_recall = metrics.kmer_recall,
            kmer_f1 = metrics.kmer_f1,
            effective_graph_family = get(assembly_stats, "effective_graph_family", missing),
            effective_memory_profile = get(assembly_stats, "effective_memory_profile", missing),
            tokenizer_used = get(assembly_stats, "tokenizer_used", missing),
            status = "ok",
            error_message = missing
        ))
    catch error_value
        push!(rows, (
            benchmark = "rhizomorph_reconstruction_suite",
            dataset_type = "isolate",
            scale = scale,
            strategy = string(strategy),
            genome_id = condition.genome_id,
            genome_tier = condition.genome_tier,
            genome_length = condition.genome_length,
            q_value = condition.q_value,
            error_rate = bundle.error_rate,
            coverage = condition.coverage,
            replicate = condition.replicate,
            algorithm_id = condition.algorithm_id,
            graph_family = string(condition.graph_family),
            memory_profile = string(condition.memory_profile),
            tokenizer = string(condition.tokenizer),
            input_mode = string(condition.input_mode),
            stage = string(condition.stage),
            num_reads = bundle.n_reads,
            run_seconds = missing,
            peak_memory_bytes = missing,
            num_contigs = missing,
            total_contig_length = missing,
            assembly_size_ratio = missing,
            n50 = missing,
            l50 = missing,
            kmer_precision = missing,
            kmer_recall = missing,
            kmer_f1 = missing,
            effective_graph_family = missing,
            effective_memory_profile = missing,
            tokenizer_used = missing,
            status = "error",
            error_message = string(typeof(error_value), ": ", error_value)
        ))
    end
end

results = DataFrames.DataFrame(rows)
add_profile_delta_columns!(results)

results_dir = mkpath("results")
timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
csv_file = joinpath(results_dir, "09_rhizomorph_reconstruction_suite_$(timestamp).csv")
json_file = joinpath(results_dir, "09_rhizomorph_reconstruction_suite_$(timestamp).json")
CSV.write(csv_file, results)

ok_rows = results[results.status .== "ok", :]
summary = Dict(
    "num_conditions" => length(conditions),
    "num_rows" => DataFrames.nrow(results),
    "num_ok" => DataFrames.nrow(ok_rows),
    "num_errors" => sum(results.status .== "error"),
    "median_runtime_seconds" => DataFrames.nrow(ok_rows) == 0 ? missing :
                                Statistics.median(ok_rows.run_seconds),
    "median_kmer_f1" => DataFrames.nrow(ok_rows) == 0 ? missing :
                        Statistics.median(ok_rows.kmer_f1)
)

open(json_file, "w") do io
    JSON.print(
        io,
        Dict(
            "benchmark_name" => "rhizomorph_reconstruction_suite",
            "timestamp" => string(Dates.now()),
            "scale" => scale,
            "strategy" => string(strategy),
            "summary" => summary,
            "results" => _to_json_row_table(results)
        ),
        2
    )
end

println("Wrote CSV: $(csv_file)")
println("Wrote JSON: $(json_file)")
println("Done at: $(Dates.now())")
