# B8 Viterbi accuracy-vs-error-rate benchmark for real viroid, phage, and text fixtures.
#
# Usage:
#   julia --project=. benchmarking/viterbi_accuracy_benchmark.jl
#   julia --project=. benchmarking/viterbi_accuracy_benchmark.jl --output-dir /tmp/b8 --skip-plots
#   julia --project=. benchmarking/viterbi_accuracy_benchmark.jl --beam-widths exact,1,4,16 --window-scale 4 --repetitions 3 --skip-plots

import BioSequences
if !("--skip-plots" in ARGS)
    @eval import CairoMakie
end
import DataFrames
import Dates
import FASTX
import Kmers

# This benchmark exercises Rhizomorph graph builders plus generalized Viterbi
# correction only. Loading Mycelia in core-benchmark mode avoids unrelated
# plotting/ML imports during HPC smoke runs.
if get(ENV, "MYCELIA_CORE_BENCHMARK", "") == ""
    ENV["MYCELIA_CORE_BENCHMARK"] = "true"
end

import Mycelia
import Test

include(joinpath(@__DIR__, "benchmark_artifacts.jl"))

function _default_viterbi_accuracy_run_id()::String
    timestamp = Dates.format(Dates.now(Dates.UTC), Dates.DateFormat("yyyymmddTHHMMSS"))
    return "b8_viterbi_accuracy_$(timestamp)Z"
end

const VITERBI_ACCURACY_ERROR_RATES = (0.01, 0.05, 0.10)
const VITERBI_ACCURACY_FIXTURE_DIR = joinpath(@__DIR__, "fixtures", "viterbi_accuracy")
const VITERBI_ACCURACY_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "viterbi_accuracy_b8"
)
const VITERBI_ACCURACY_K = 9
const VITERBI_ACCURACY_BASE_VIROID_LENGTH = 360
const VITERBI_ACCURACY_BASE_PHAGE_LENGTH = 432
const VITERBI_ACCURACY_BASE_TEXT_LENGTH = 432

struct ViterbiAccuracyFixture
    dataset_id::String
    dataset_name::String
    domain::Symbol
    alphabet::Vector{Char}
    graph
    truth::Vector{String}
    source::String
end

function main(args::Vector{String} = ARGS)::Nothing
    output_dir = _viterbi_accuracy_arg_value(args, "--output-dir", VITERBI_ACCURACY_DEFAULT_OUTPUT_DIR)
    run_id = _viterbi_accuracy_arg_value(args, "--run-id", _default_viterbi_accuracy_run_id())
    scale = _viterbi_accuracy_arg_value(args, "--scale", "local-smoke")
    beam_widths = _viterbi_accuracy_beam_widths(args)
    error_rates = _viterbi_accuracy_error_rates(args)
    window_scale = _viterbi_accuracy_positive_float(args, "--window-scale", 1.0)
    repetitions = _viterbi_accuracy_positive_int(args, "--repetitions", 1)
    write_plots = !("--skip-plots" in args)
    artifacts = run_viterbi_accuracy_benchmark(
        output_dir;
        run_id = run_id,
        scale = scale,
        command_args = [
            "julia",
            "--project=.",
            "benchmarking/viterbi_accuracy_benchmark.jl",
            args...
        ],
        beam_widths = beam_widths,
        error_rates = error_rates,
        window_scale = window_scale,
        repetitions = repetitions,
        write_plots = write_plots
    )
    println("Wrote B8 Viterbi accuracy benchmark artifacts:")
    println("  root: $(artifacts.root)")
    println("  summary: $(artifacts.summary_csv)")
    println("  figure_png: $(artifacts.figure_png)")
    println("  figure_svg: $(artifacts.figure_svg)")
    return nothing
end

function run_viterbi_accuracy_benchmark(
        output_dir::AbstractString = VITERBI_ACCURACY_DEFAULT_OUTPUT_DIR;
        run_id::AbstractString = _default_viterbi_accuracy_run_id(),
        scale::AbstractString = "local-smoke",
        command_args::Vector{String} = [
            "julia",
            "--project=.",
            "benchmarking/viterbi_accuracy_benchmark.jl"
        ],
        beam_widths::Vector{Union{Nothing, Int}} = Union{Nothing, Int}[nothing],
        error_rates::Vector{Float64} = collect(VITERBI_ACCURACY_ERROR_RATES),
        window_scale::Float64 = 1.0,
        repetitions::Int = 1,
        write_plots::Bool = true
)::NamedTuple
    fixtures = viterbi_accuracy_fixtures(window_scale = window_scale)
    rows = NamedTuple[]

    for fixture in fixtures
        for target_error_rate in error_rates
            for beam_width in beam_widths
                for replicate in 1:repetitions
                    push!(
                        rows,
                        _viterbi_accuracy_row(
                            fixture,
                            target_error_rate,
                            beam_width,
                            replicate,
                            window_scale
                        )
                    )
                end
            end
        end
    end

    summary = DataFrames.DataFrame(rows)
    artifacts = write_benchmark_artifacts(
        ["viterbi_accuracy_summary" => summary];
        output_dir = output_dir,
        run_id = run_id,
        scale = scale,
        dataset_ids = [fixture.dataset_id for fixture in fixtures],
        command_args = command_args,
        metadata = Dict(
            "bead" => "td-he0z.9",
            "benchmark" => "generalized_viterbi_accuracy_vs_error_rate",
            "baseline" => "uncorrected injected-error observations",
            "error_rates" => error_rates,
            "beam_widths" => [_beam_width_label(beam_width) for beam_width in beam_widths],
            "fixture_dir" => relpath(VITERBI_ACCURACY_FIXTURE_DIR, @__DIR__),
            "repetitions" => repetitions,
            "window_scale" => window_scale,
            "unit" => "fixed-length $(VITERBI_ACCURACY_K)-mers/ngrams"
        ),
        table_context_columns = Dict(
            "viterbi_accuracy_summary" => Dict("benchmark_dataset_id" => "dataset_id")
        )
    )

    figure_png = ""
    figure_svg = ""
    if write_plots
        figure_paths = _write_viterbi_accuracy_figure(summary, artifacts.layout.plots)
        figure_png = figure_paths.png
        figure_svg = figure_paths.svg
    end

    return (
        root = artifacts.root,
        summary_csv = artifacts.tables["viterbi_accuracy_summary"].table,
        index = artifacts.index,
        provenance = artifacts.provenance,
        figure_png = figure_png,
        figure_svg = figure_svg,
        rows = DataFrames.nrow(summary)
    )
end

function viterbi_accuracy_fixtures(; window_scale::Float64 = 1.0)::Vector{ViterbiAccuracyFixture}
    pstvd = replace(
        _read_fasta_sequence(joinpath(VITERBI_ACCURACY_FIXTURE_DIR, "pstvd_nc002030.fasta")),
        'T' => 'U'
    )
    phix = _read_fasta_sequence(joinpath(VITERBI_ACCURACY_FIXTURE_DIR, "phix174_nc001422.fasta"))
    text = read(joinpath(VITERBI_ACCURACY_FIXTURE_DIR, "pride_and_prejudice_excerpt.txt"), String)

    viroid_length = _scaled_fixture_length(
        length(pstvd),
        VITERBI_ACCURACY_BASE_VIROID_LENGTH,
        window_scale
    )
    phage_length = _scaled_fixture_length(
        length(phix),
        VITERBI_ACCURACY_BASE_PHAGE_LENGTH,
        window_scale
    )
    text_length = _scaled_fixture_length(
        length(text),
        VITERBI_ACCURACY_BASE_TEXT_LENGTH,
        window_scale
    )

    viroid_window = first(pstvd, viroid_length)
    phage_window = first(phix, phage_length)
    text_window = _normalized_text_window(text, text_length)

    viroid_record = FASTX.FASTA.Record(
        "pstvd_nc002030", BioSequences.LongRNA{4}(viroid_window)
    )
    phage_record = FASTX.FASTA.Record(
        "phix174_nc001422", BioSequences.LongDNA{4}(phage_window)
    )

    return ViterbiAccuracyFixture[
        ViterbiAccuracyFixture(
            "pstvd_nc002030",
            "Potato spindle tuber viroid (NC_002030.1)",
            :viroid,
            collect("ACGU"),
            Mycelia.Rhizomorph.build_kmer_graph(
                [viroid_record], VITERBI_ACCURACY_K;
                dataset_id = "b8_pstvd_nc002030",
                mode = :singlestrand,
                type_hint = :RNA
            ),
            _sequence_kmers(viroid_window, VITERBI_ACCURACY_K),
            "NCBI RefSeq NC_002030.1"
        ),
        ViterbiAccuracyFixture(
            "phix174_nc001422",
            "Escherichia phage phiX174 (NC_001422.1)",
            :phage,
            collect("ACGT"),
            Mycelia.Rhizomorph.build_kmer_graph(
                [phage_record], VITERBI_ACCURACY_K;
                dataset_id = "b8_phix174_nc001422",
                mode = :singlestrand
            ),
            _sequence_kmers(phage_window, VITERBI_ACCURACY_K),
            "NCBI RefSeq NC_001422.1"
        ),
        ViterbiAccuracyFixture(
            "austen_pride_prejudice_excerpt",
            "Pride and Prejudice text excerpt",
            :text,
            collect("abcdefghijklmnopqrstuvwxyz ,.-"),
            Mycelia.Rhizomorph.build_ngram_graph(
                [text_window], VITERBI_ACCURACY_K; dataset_id = "b8_austen_text"
            ),
            _sequence_kmers(text_window, VITERBI_ACCURACY_K),
            "Public-domain Pride and Prejudice excerpt"
        )
    ]
end

function _viterbi_accuracy_row(
        fixture::ViterbiAccuracyFixture,
        target_error_rate::Float64,
        beam_width::Union{Nothing, Int},
        replicate::Int,
        window_scale::Float64
)::NamedTuple
    observed, injected_positions = _inject_fixture_errors(
        fixture.truth, fixture.alphabet, target_error_rate, replicate
    )
    converted_observed = _convert_observations(observed, fixture.domain)
    config = Mycelia.ViterbiCorrectionConfig(
        error_rate = target_error_rate,
        beam_width = beam_width
    )
    correction_start_ns = time_ns()
    result = Mycelia.correct_observations(fixture.graph, [converted_observed]; config = config)
    correction_wall_seconds = (time_ns() - correction_start_ns) / 1.0e9
    corrected_observation = only(result.corrected_observations)
    correction_returned_path = !isnothing(corrected_observation)
    raw_corrected_count = correction_returned_path ? length(corrected_observation) : 0
    raw_corrected = if correction_returned_path
        [string(observation) for observation in corrected_observation]
    else
        String[]
    end
    correction_succeeded = correction_returned_path && length(raw_corrected) == length(fixture.truth)
    corrected = correction_succeeded ? raw_corrected : copy(observed)
    benchmark_reason = if !correction_returned_path
        "no_path"
    elseif !correction_succeeded
        "length_mismatch"
    else
        string(get(result.diagnostics, :reason, ""))
    end

    baseline_edit_distance = _sum_edit_distance(observed, fixture.truth)
    corrected_edit_distance = _sum_edit_distance(corrected, fixture.truth)
    injected_error_count = length(injected_positions)
    recovered_errors = _recovered_error_count(corrected, fixture.truth, injected_positions)
    editable_positions = sum(length, fixture.truth[2:(end - 1)])
    actual_error_rate = injected_error_count / editable_positions
    edit_distance_reduction = if baseline_edit_distance == 0
        0.0
    else
        (baseline_edit_distance - corrected_edit_distance) / baseline_edit_distance
    end
    recall = if injected_error_count == 0
        0.0
    else
        recovered_errors / injected_error_count
    end

    if correction_succeeded
        Test.@test length(corrected) == length(fixture.truth)
    end

    return (
        dataset_id = fixture.dataset_id,
        dataset_name = fixture.dataset_name,
        domain = string(fixture.domain),
        source = fixture.source,
        replicate = replicate,
        window_scale = window_scale,
        target_error_rate = target_error_rate,
        actual_error_rate = actual_error_rate,
        beam_width = _beam_width_label(beam_width),
        beam_width_limit = isnothing(beam_width) ? 0 : beam_width,
        observation_count = length(fixture.truth),
        editable_positions = editable_positions,
        injected_error_count = injected_error_count,
        baseline_method = "uncorrected_observations",
        corrected_method = "generalized_viterbi_correct_observations",
        correction_returned_path = correction_returned_path,
        correction_succeeded = correction_succeeded,
        raw_corrected_observation_count = raw_corrected_count,
        baseline_edit_distance = baseline_edit_distance,
        corrected_edit_distance = corrected_edit_distance,
        edit_distance_reduction = edit_distance_reduction,
        injected_error_recall = recall,
        corrected_all_injected_positions = recovered_errors == injected_error_count,
        correction_wall_seconds = correction_wall_seconds,
        diagnostics_expanded_states = _diagnostic_int(result.diagnostics, :expanded_states),
        diagnostics_generated_states = _diagnostic_int(result.diagnostics, :generated_states),
        diagnostics_max_retained_states = _diagnostic_int(result.diagnostics, :max_retained_states),
        diagnostics_cumulative_retained_states = _diagnostic_int(
            result.diagnostics,
            :cumulative_retained_states
        ),
        diagnostics_alphabet = string(result.diagnostics[:alphabet]),
        diagnostics_strand_mode = string(result.diagnostics[:strand_mode]),
        diagnostics_emission_model = string(result.diagnostics[:emission_model]),
        diagnostics_transition_model = string(result.diagnostics[:transition_model]),
        diagnostics_reason = benchmark_reason
    )
end

function _read_fasta_sequence(path::AbstractString)::String
    lines = readlines(path)
    sequence_lines = [strip(line) for line in lines if !isempty(line) && !startswith(line, ">")]
    return uppercase(join(sequence_lines))
end

function _scaled_fixture_length(total_length::Int, base_length::Int, window_scale::Float64)::Int
    if window_scale <= 0.0
        throw(ArgumentError("window_scale must be positive, got $window_scale"))
    end
    scaled_length = max(VITERBI_ACCURACY_K, round(Int, base_length * window_scale))
    return min(total_length, scaled_length)
end

function _sequence_kmers(sequence::AbstractString, k::Int)::Vector{String}
    return [sequence[index:(index + k - 1)] for index in 1:(lastindex(sequence) - k + 1)]
end

function _normalized_text_window(text::AbstractString, max_length::Int)::String
    normalized = lowercase(replace(text, r"\s+" => " "))
    allowed = collect("abcdefghijklmnopqrstuvwxyz ,.-")
    normalized = filter(character -> character in allowed, normalized)
    return first(normalized, min(length(normalized), max_length))
end

function _inject_fixture_errors(
        truth::Vector{String},
        alphabet::Vector{Char},
        target_error_rate::Float64,
        replicate::Int = 1
)::Tuple{Vector{String}, Vector{Tuple{Int, Int}}}
    candidates = Tuple{Int, Int}[]
    for item_index in 2:(length(truth) - 1)
        for position in eachindex(truth[item_index])
            push!(candidates, (item_index, position))
        end
    end

    target_errors = max(1, round(Int, target_error_rate * length(candidates)))
    target_errors = min(target_errors, length(candidates))
    selected_offsets = _deterministic_error_offsets(length(candidates), target_errors, replicate)
    selected_positions = unique(candidates[selected_offsets])

    observed = copy(truth)
    for (item_index, position) in selected_positions
        characters = collect(observed[item_index])
        characters[position] = _replacement_character(
            characters[position], alphabet, item_index + position
        )
        observed[item_index] = String(characters)
    end

    return observed, selected_positions
end

function _deterministic_error_offsets(candidate_count::Int, target_errors::Int, replicate::Int)::Vector{Int}
    if replicate <= 0
        throw(ArgumentError("replicate must be positive, got $replicate"))
    end
    if target_errors <= 0 || candidate_count <= 0
        return Int[]
    end

    offsets = Int[]
    seen = Set{Int}()
    stride = max(1, floor(Int, candidate_count / target_errors))
    cursor = mod1(replicate, candidate_count)
    while length(offsets) < target_errors
        if !(cursor in seen)
            push!(offsets, cursor)
            push!(seen, cursor)
        end
        cursor = mod1(cursor + stride, candidate_count)
        if cursor in seen
            cursor = mod1(cursor + 1, candidate_count)
        end
    end
    return offsets
end

function _replacement_character(character::Char, alphabet::Vector{Char}, offset::Int)::Char
    current = uppercase(character)
    for shift in 0:(length(alphabet) - 1)
        candidate = alphabet[mod1(offset + shift, length(alphabet))]
        if candidate != current && candidate != character
            return islowercase(character) ? lowercase(candidate) : candidate
        end
    end
    error("no replacement character available for $character")
end

function _convert_observations(observed::Vector{String}, domain::Symbol)::Vector
    if domain == :viroid
        return [Kmers.RNAKmer{VITERBI_ACCURACY_K}(observation) for observation in observed]
    elseif domain == :phage
        return [Kmers.DNAKmer{VITERBI_ACCURACY_K}(observation) for observation in observed]
    elseif domain == :text
        return observed
    end
    throw(ArgumentError("unsupported Viterbi accuracy domain: $domain"))
end

function _sum_edit_distance(left::Vector{String}, right::Vector{String})::Int
    return sum(_edit_distance(left_item, right_item) for (left_item, right_item) in zip(left, right))
end

function _edit_distance(left::AbstractString, right::AbstractString)::Int
    left_chars = collect(left)
    right_chars = collect(right)
    distances = zeros(Int, length(left_chars) + 1, length(right_chars) + 1)

    for index in 0:length(left_chars)
        distances[index + 1, 1] = index
    end
    for index in 0:length(right_chars)
        distances[1, index + 1] = index
    end

    for left_index in 1:length(left_chars)
        for right_index in 1:length(right_chars)
            substitution_cost = left_chars[left_index] == right_chars[right_index] ? 0 : 1
            distances[left_index + 1, right_index + 1] = min(
                distances[left_index, right_index + 1] + 1,
                distances[left_index + 1, right_index] + 1,
                distances[left_index, right_index] + substitution_cost
            )
        end
    end

    return distances[end, end]
end

function _recovered_error_count(
        corrected::Vector{String},
        truth::Vector{String},
        injected_positions::Vector{Tuple{Int, Int}}
)::Int
    recovered = 0
    for (item_index, position) in injected_positions
        recovered += corrected[item_index][position] == truth[item_index][position] ? 1 : 0
    end
    return recovered
end

function _write_viterbi_accuracy_figure(
        summary::DataFrames.DataFrame,
        plots_dir::AbstractString
)::NamedTuple
    if !isdefined(@__MODULE__, :CairoMakie)
        throw(
            ArgumentError(
                "CairoMakie was not loaded; rerun without --skip-plots or " *
                "import CairoMakie before calling this function"
            )
        )
    end

    mkpath(plots_dir)
    fig = CairoMakie.Figure(size = (900, 600), fontsize = 16)
    axis = CairoMakie.Axis(
        fig[1, 1],
        title = "B8 Viterbi correction accuracy vs injected error rate",
        xlabel = "Injected error rate (%)",
        ylabel = "Edit-distance reduction vs uncorrected baseline",
        yticks = 0.0:0.25:1.0
    )

    palette = [:seagreen4, :dodgerblue3, :darkorange3]
    grouped = collect(DataFrames.groupby(summary, :dataset_id))
    elements = CairoMakie.LineElement[]
    labels = String[]
    for (index, group) in enumerate(grouped)
        sorted_group = sort(group, :target_error_rate)
        color = palette[mod1(index, length(palette))]
        x_values = 100 .* sorted_group.target_error_rate
        y_values = sorted_group.edit_distance_reduction
        CairoMakie.lines!(axis, x_values, y_values; color = color, linewidth = 3)
        CairoMakie.scatter!(axis, x_values, y_values; color = color, markersize = 12)
        push!(elements, CairoMakie.LineElement(color = color, linewidth = 3))
        push!(labels, first(sorted_group.dataset_name))
    end

    CairoMakie.ylims!(axis, -0.05, 1.05)
    CairoMakie.axislegend(axis, elements, labels; position = :rb)
    png_path = joinpath(plots_dir, "viterbi_accuracy_vs_error_rate.png")
    svg_path = joinpath(plots_dir, "viterbi_accuracy_vs_error_rate.svg")
    CairoMakie.save(png_path, fig)
    CairoMakie.save(svg_path, fig)
    return (png = png_path, svg = svg_path)
end

function _diagnostic_int(diagnostics::AbstractDict, key::Symbol)::Int
    return Int(get(diagnostics, key, 0))
end

function _beam_width_label(beam_width::Union{Nothing, Int})::String
    return isnothing(beam_width) ? "exact" : string(beam_width)
end

function _viterbi_accuracy_beam_widths(args::Vector{String})::Vector{Union{Nothing, Int}}
    raw = _viterbi_accuracy_arg_value(args, "--beam-widths", "exact")
    beam_widths = Union{Nothing, Int}[]
    for item in split(raw, ",")
        normalized = lowercase(strip(item))
        if normalized in ("", "exact", "none", "nothing", "unbounded")
            push!(beam_widths, nothing)
        else
            beam_width = parse(Int, normalized)
            if beam_width <= 0
                throw(ArgumentError("beam widths must be positive, got $beam_width"))
            end
            push!(beam_widths, beam_width)
        end
    end
    return unique(beam_widths)
end

function _viterbi_accuracy_error_rates(args::Vector{String})::Vector{Float64}
    raw = _viterbi_accuracy_arg_value(
        args,
        "--error-rates",
        join(string.(VITERBI_ACCURACY_ERROR_RATES), ",")
    )
    error_rates = Float64[]
    for item in split(raw, ",")
        normalized = strip(item)
        if isempty(normalized)
            continue
        end
        error_rate = parse(Float64, normalized)
        if error_rate <= 0.0 || error_rate >= 0.5
            throw(ArgumentError("error rates must be in (0, 0.5), got $error_rate"))
        end
        push!(error_rates, error_rate)
    end
    if isempty(error_rates)
        throw(ArgumentError("at least one error rate is required"))
    end
    return unique(error_rates)
end

function _viterbi_accuracy_positive_float(
        args::Vector{String},
        flag::AbstractString,
        default::Float64
)::Float64
    value = parse(Float64, _viterbi_accuracy_arg_value(args, flag, string(default)))
    if value <= 0.0
        throw(ArgumentError("$flag must be positive, got $value"))
    end
    return value
end

function _viterbi_accuracy_positive_int(
        args::Vector{String},
        flag::AbstractString,
        default::Int
)::Int
    value = parse(Int, _viterbi_accuracy_arg_value(args, flag, string(default)))
    if value <= 0
        throw(ArgumentError("$flag must be positive, got $value"))
    end
    return value
end

function _viterbi_accuracy_arg_value(
        args::Vector{String},
        flag::AbstractString,
        default::AbstractString
)::String
    index = findfirst(==(flag), args)
    if isnothing(index) || index == length(args)
        return string(default)
    end
    return args[index + 1]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
