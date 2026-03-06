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

println("=== Rhizomorph Reconstruction Suite (Metagenome) ===")
println("Start time: $(Dates.now())")

function _parse_int_list_env(var_name::String, default_values::Vector{Int})
    value = strip(get(ENV, var_name, ""))
    if isempty(value)
        return default_values
    end
    return [parse(Int, strip(token)) for token in split(value, ",") if !isempty(strip(token))]
end

function _parse_symbol_list_env(var_name::String, default_values::Vector{Symbol})
    value = strip(get(ENV, var_name, ""))
    if isempty(value)
        return default_values
    end
    return [Symbol(strip(token)) for token in split(value, ",") if !isempty(strip(token))]
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
            complexity_levels = METAGENOME_COMPLEXITY_LEVELS,
            depth_levels = COVERAGE_LEVELS,
            quality_profiles = vcat(
                METAGENOME_PRIMARY_QUALITY_PROFILES,
                [:mixed_balanced, :mixed_extremes]
            ),
            abundance_profiles = [:equal, :log_normal],
            replicates = 3,
            strategy = :tiered,
            anchor_depth_levels = METAGENOME_ANCHOR_DEPTH_LEVELS,
            anchor_quality_profiles = METAGENOME_ANCHOR_QUALITY_PROFILES
        )
    elseif scale == "medium"
        return (
            complexity_levels = METAGENOME_COMPLEXITY_LEVELS,
            depth_levels = [10, 25, 50, 100, 250],
            quality_profiles = vcat(METAGENOME_PRIMARY_QUALITY_PROFILES, [:mixed_balanced]),
            abundance_profiles = [:equal, :log_normal],
            replicates = 2,
            strategy = :tiered,
            anchor_depth_levels = [25, 100, 250],
            anchor_quality_profiles = METAGENOME_ANCHOR_QUALITY_PROFILES
        )
    end
    return (
        complexity_levels = [10, 25],
        depth_levels = [10, 50, 250],
        quality_profiles = [:const_q10, :const_q30, :const_q50, :mixed_balanced],
        abundance_profiles = [:equal, :log_normal],
        replicates = 1,
        strategy = :tiered,
        anchor_depth_levels = [50],
        anchor_quality_profiles = [:const_q30, :mixed_balanced]
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
default_max_records_per_genome = if scale == "large"
    30_000
elseif scale == "medium"
    15_000
else
    8_000
end
default_max_total_records = if scale == "large"
    250_000
elseif scale == "medium"
    120_000
else
    60_000
end

complexity_levels = _parse_int_list_env("MYCELIA_BENCHMARK_METAGENOME_COMPLEXITY",
    defaults.complexity_levels)
depth_levels = _parse_int_list_env("MYCELIA_BENCHMARK_METAGENOME_DEPTHS", defaults.depth_levels)
quality_profiles = _parse_symbol_list_env("MYCELIA_BENCHMARK_METAGENOME_QUALITY_PROFILES",
    defaults.quality_profiles)
abundance_profiles = _parse_symbol_list_env("MYCELIA_BENCHMARK_METAGENOME_ABUNDANCE_PROFILES",
    defaults.abundance_profiles)
anchor_depth_levels = _parse_int_list_env("MYCELIA_BENCHMARK_METAGENOME_ANCHOR_DEPTHS",
    defaults.anchor_depth_levels)
anchor_quality_profiles = _parse_symbol_list_env(
    "MYCELIA_BENCHMARK_METAGENOME_ANCHOR_QUALITY_PROFILES",
    defaults.anchor_quality_profiles
)
replicates = parse(Int, get(ENV, "MYCELIA_BENCHMARK_REPLICATES", string(defaults.replicates)))
include_advanced = lowercase(get(ENV, "MYCELIA_BENCHMARK_ENABLE_ADVANCED", "false")) ==
                   "true"
max_runs = parse(Int, get(ENV, "MYCELIA_BENCHMARK_MAX_RUNS", "0"))
max_records_per_genome = parse(Int, get(ENV, "MYCELIA_BENCHMARK_MAX_RECORDS_PER_GENOME",
    string(default_max_records_per_genome)))
max_total_records = parse(Int, get(ENV, "MYCELIA_BENCHMARK_MAX_TOTAL_RECORDS",
    string(default_max_total_records)))

if strategy == :full_factorial
    anchor_depth_levels = depth_levels
    anchor_quality_profiles = quality_profiles
end

catalog_settings = metagenome_catalog_settings(scale)
catalog = build_synthetic_metagenome_catalog(
    n_references = catalog_settings.n_references,
    min_length = catalog_settings.min_length,
    max_length = catalog_settings.max_length
)

conditions = build_metagenome_benchmark_matrix(
    complexity_levels = complexity_levels,
    depth_levels = depth_levels,
    abundance_profiles = abundance_profiles,
    quality_profiles = quality_profiles,
    replicates = replicates,
    strategy = strategy,
    include_advanced = include_advanced,
    anchor_depth_levels = anchor_depth_levels,
    anchor_quality_profiles = anchor_quality_profiles
)

if max_runs > 0 && length(conditions) > max_runs
    conditions = conditions[1:max_runs]
end

println("Scale: $(scale)")
println("Strategy: $(strategy)")
println("Catalog references: $(length(catalog.sequences))")
println("Conditions: $(length(conditions))")

rows = NamedTuple[]
bundle_cache = Dict{Tuple{Int, Int, Symbol, Symbol, Int}, Any}()
abundance_index = Dict(profile => idx for (idx, profile) in enumerate(abundance_profiles))
quality_index = Dict(profile => idx for (idx, profile) in enumerate(quality_profiles))

for (index, condition) in enumerate(conditions)
    println("[$(index)/$(length(conditions))] $(condition.algorithm_id) | community=$(condition.community_size) | depth=$(condition.depth_target)x | $(condition.quality_profile) | $(condition.abundance_profile) | rep=$(condition.replicate)")

    bundle_key = (
        condition.community_size,
        condition.depth_target,
        Symbol(condition.abundance_profile),
        Symbol(condition.quality_profile),
        condition.replicate
    )

    if !haskey(bundle_cache, bundle_key)
        simulation_seed = 100_000 +
                          (condition.replicate * 10_000) +
                          (condition.depth_target * 100) +
                          (condition.community_size * 7) +
                          (abundance_index[condition.abundance_profile] * 53) +
                          (quality_index[condition.quality_profile] * 97)
        bundle_cache[bundle_key] = generate_metagenome_read_bundle(
            catalog_sequences = catalog.sequences,
            community_size = condition.community_size,
            depth_target = condition.depth_target,
            abundance_profile = condition.abundance_profile,
            quality_profile = condition.quality_profile,
            seed = simulation_seed,
            max_records_per_genome = max_records_per_genome,
            max_total_records = max_total_records
        )
    end

    bundle = bundle_cache[bundle_key]
    benchmark_input = select_benchmark_input(bundle, condition.input_mode)
    reference_sequences = bundle.truth_sequences

    mean_q = DataFrames.nrow(bundle.truth_table) == 0 ? missing :
             Statistics.mean(bundle.truth_table.q_value)
    std_q = DataFrames.nrow(bundle.truth_table) <= 1 ? 0.0 :
            Statistics.std(bundle.truth_table.q_value)
    mean_coverage = DataFrames.nrow(bundle.truth_table) == 0 ? missing :
                    Statistics.mean(bundle.truth_table.coverage)
    std_coverage = DataFrames.nrow(bundle.truth_table) <= 1 ? 0.0 :
                   Statistics.std(bundle.truth_table.coverage)

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
            benchmark = "rhizomorph_metagenome_suite",
            dataset_type = "metagenome",
            scale = scale,
            strategy = string(strategy),
            community_size = condition.community_size,
            depth_target = condition.depth_target,
            abundance_profile = string(condition.abundance_profile),
            quality_profile = string(condition.quality_profile),
            replicate = condition.replicate,
            algorithm_id = condition.algorithm_id,
            graph_family = string(condition.graph_family),
            memory_profile = string(condition.memory_profile),
            tokenizer = string(condition.tokenizer),
            input_mode = string(condition.input_mode),
            stage = string(condition.stage),
            num_reads = length(bundle.fastq),
            num_truth_genomes = length(bundle.selected_ids),
            mean_q_value = mean_q,
            std_q_value = std_q,
            mean_coverage = mean_coverage,
            std_coverage = std_coverage,
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
            benchmark = "rhizomorph_metagenome_suite",
            dataset_type = "metagenome",
            scale = scale,
            strategy = string(strategy),
            community_size = condition.community_size,
            depth_target = condition.depth_target,
            abundance_profile = string(condition.abundance_profile),
            quality_profile = string(condition.quality_profile),
            replicate = condition.replicate,
            algorithm_id = condition.algorithm_id,
            graph_family = string(condition.graph_family),
            memory_profile = string(condition.memory_profile),
            tokenizer = string(condition.tokenizer),
            input_mode = string(condition.input_mode),
            stage = string(condition.stage),
            num_reads = length(bundle.fastq),
            num_truth_genomes = length(bundle.selected_ids),
            mean_q_value = mean_q,
            std_q_value = std_q,
            mean_coverage = mean_coverage,
            std_coverage = std_coverage,
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
csv_file = joinpath(results_dir, "10_rhizomorph_metagenome_suite_$(timestamp).csv")
json_file = joinpath(results_dir, "10_rhizomorph_metagenome_suite_$(timestamp).json")
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
            "benchmark_name" => "rhizomorph_metagenome_suite",
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
