# B8 Viterbi accuracy-vs-error-rate benchmark for real viroid, phage, and text fixtures.
#
# Usage:
#   julia --project=. benchmarking/viterbi_accuracy_benchmark.jl
#   julia --project=. benchmarking/viterbi_accuracy_benchmark.jl --output-dir /tmp/b8 --skip-plots

import BioSequences
import CairoMakie
import DataFrames
import FASTX
import Kmers
import MetaGraphsNext
import Mycelia
import Random
import Test

include(joinpath(@__DIR__, "benchmark_artifacts.jl"))

const VITERBI_ACCURACY_ERROR_RATES = (0.01, 0.05, 0.10)
# Fixed seed for the shuffled-edge-weight null (Control B). Deterministic so the
# null table is byte-reproducible across runs; only the permutation of the real
# edge weights is randomized, never the topology, observations, or emission model.
const VITERBI_NULL_SEED = 20260711
const VITERBI_ACCURACY_FIXTURE_DIR = joinpath(@__DIR__, "fixtures", "viterbi_accuracy")
const VITERBI_ACCURACY_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "viterbi_accuracy_b8"
)
# PRIME k (was 9 = 3x3, the worst period-3 case). A composite k aliases
# period-p tandem repeats onto self-overlapping k-mers, which can make single-k
# correction look either trivially perfect or pathologically wrong; a prime k
# breaks that periodicity. Matches the prime-only k progression the iterative
# corrector uses (find_initial_k draws from Primes.primes; build_k_ladder and
# next_prime_k snap to primes).
const VITERBI_ACCURACY_K = 11

struct ViterbiAccuracyFixture
    dataset_id::String
    dataset_name::String
    domain::Symbol
    alphabet::Vector{Char}
    graph::Any
    truth::Vector{String}
    source::String
end

function main(args::Vector{String} = ARGS)::Nothing
    output_dir = _viterbi_accuracy_arg_value(args, "--output-dir", VITERBI_ACCURACY_DEFAULT_OUTPUT_DIR)
    write_plots = !("--skip-plots" in args)
    artifacts = run_viterbi_accuracy_benchmark(output_dir; write_plots = write_plots)
    println("Wrote B8 Viterbi accuracy benchmark artifacts:")
    println("  root: $(artifacts.root)")
    println("  summary: $(artifacts.summary_csv)")
    println("  figure_png: $(artifacts.figure_png)")
    println("  figure_svg: $(artifacts.figure_svg)")
    return nothing
end

function run_viterbi_accuracy_benchmark(
        output_dir::AbstractString = VITERBI_ACCURACY_DEFAULT_OUTPUT_DIR;
        write_plots::Bool = true
)::NamedTuple
    fixtures = viterbi_accuracy_fixtures()
    rows = NamedTuple[]

    for fixture in fixtures
        for target_error_rate in VITERBI_ACCURACY_ERROR_RATES
            push!(rows, _viterbi_accuracy_row(fixture, target_error_rate))
        end
    end

    summary = DataFrames.DataFrame(rows)

    # Control A (over-correction on un-corrupted input) and Control B (shuffled-
    # edge-weight null) — the two controls the manuscript flags as missing next
    # to the recall table. Same fixtures, same prime k, same corrector.
    overcorrection_rows = NamedTuple[]
    null_rows = NamedTuple[]
    for fixture in fixtures
        for target_error_rate in VITERBI_ACCURACY_ERROR_RATES
            push!(overcorrection_rows, _viterbi_overcorrection_row(fixture, target_error_rate))
            push!(null_rows, _viterbi_null_row(fixture, target_error_rate))
        end
    end
    overcorrection = DataFrames.DataFrame(overcorrection_rows)
    null_control = DataFrames.DataFrame(null_rows)

    artifacts = write_benchmark_artifacts(
        [
            "viterbi_accuracy_summary" => summary,
            "viterbi_overcorrection_summary" => overcorrection,
            "viterbi_null_control_summary" => null_control
        ];
        output_dir = output_dir,
        run_id = "b8_viterbi_accuracy_local_20260625",
        scale = "local-smoke",
        dataset_ids = [fixture.dataset_id for fixture in fixtures],
        command_args = [
            "julia", "--project=.", "benchmarking/viterbi_accuracy_benchmark.jl"],
        metadata = Dict(
            "bead" => "td-he0z.9",
            "benchmark" => "generalized_viterbi_accuracy_vs_error_rate",
            "baseline" => "uncorrected injected-error observations",
            "error_rates" => collect(VITERBI_ACCURACY_ERROR_RATES),
            "fixture_dir" => relpath(VITERBI_ACCURACY_FIXTURE_DIR, @__DIR__),
            "unit" => "fixed-length $(VITERBI_ACCURACY_K)-mers/ngrams",
            "controls" => Dict(
                "over_correction_uncorrupted" => "corrector run on uncorrupted truth; any edit is a false positive",
                "shuffled_weight_null" => "edge weights permuted (seed $(VITERBI_NULL_SEED)); topology/emission held fixed",
                "random_rewire_null" => "edges rewired to random vertex pairs (seed $(VITERBI_NULL_SEED)); vertices/edge-count/weight-multiset held fixed; destroys true adjacency"
            )
        ),
        table_context_columns = Dict(
            "viterbi_accuracy_summary" => Dict("benchmark_dataset_id" => "dataset_id"),
            "viterbi_overcorrection_summary" => Dict("benchmark_dataset_id" => "dataset_id"),
            "viterbi_null_control_summary" => Dict("benchmark_dataset_id" => "dataset_id")
        )
    )

    figure_png = ""
    figure_svg = ""
    control_figure_png = ""
    control_figure_svg = ""
    if write_plots
        figure_paths = _write_viterbi_accuracy_figure(summary, artifacts.layout.plots)
        figure_png = figure_paths.png
        figure_svg = figure_paths.svg
        control_paths = _write_viterbi_control_figure(
            overcorrection, null_control, artifacts.layout.plots)
        control_figure_png = control_paths.png
        control_figure_svg = control_paths.svg
    end

    return (
        root = artifacts.root,
        summary_csv = artifacts.tables["viterbi_accuracy_summary"].table,
        overcorrection_csv = artifacts.tables["viterbi_overcorrection_summary"].table,
        null_control_csv = artifacts.tables["viterbi_null_control_summary"].table,
        index = artifacts.index,
        provenance = artifacts.provenance,
        figure_png = figure_png,
        figure_svg = figure_svg,
        control_figure_png = control_figure_png,
        control_figure_svg = control_figure_svg,
        rows = DataFrames.nrow(summary),
        overcorrection_rows = DataFrames.nrow(overcorrection),
        null_control_rows = DataFrames.nrow(null_control)
    )
end

function viterbi_accuracy_fixtures()::Vector{ViterbiAccuracyFixture}
    pstvd = replace(
        _read_fasta_sequence(joinpath(VITERBI_ACCURACY_FIXTURE_DIR, "pstvd_nc002030.fasta")),
        'T' => 'U'
    )
    phix = _read_fasta_sequence(joinpath(VITERBI_ACCURACY_FIXTURE_DIR, "phix174_nc001422.fasta"))
    text = read(joinpath(VITERBI_ACCURACY_FIXTURE_DIR, "pride_and_prejudice_excerpt.txt"), String)

    viroid_window = first(pstvd, min(length(pstvd), 360))
    phage_window = first(phix, min(length(phix), 432))
    text_window = _normalized_text_window(text, 432)

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
        target_error_rate::Float64
)::NamedTuple
    observed,
    injected_positions = _inject_fixture_errors(
        fixture.truth, fixture.alphabet, target_error_rate
    )
    converted_observed = _convert_observations(observed, fixture.domain)
    config = Mycelia.ViterbiCorrectionConfig(error_rate = target_error_rate)
    result = Mycelia.correct_observations(fixture.graph, [converted_observed]; config = config)
    corrected = [string(observation) for observation in only(result.corrected_observations)]

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

    Test.@test length(corrected) == length(fixture.truth)

    return (
        dataset_id = fixture.dataset_id,
        dataset_name = fixture.dataset_name,
        domain = string(fixture.domain),
        source = fixture.source,
        target_error_rate = target_error_rate,
        actual_error_rate = actual_error_rate,
        observation_count = length(fixture.truth),
        editable_positions = editable_positions,
        injected_error_count = injected_error_count,
        baseline_method = "uncorrected_observations",
        corrected_method = "generalized_viterbi_correct_observations",
        baseline_edit_distance = baseline_edit_distance,
        corrected_edit_distance = corrected_edit_distance,
        edit_distance_reduction = edit_distance_reduction,
        injected_error_recall = recall,
        corrected_all_injected_positions = recovered_errors == injected_error_count,
        diagnostics_alphabet = string(result.diagnostics[:alphabet]),
        diagnostics_strand_mode = string(result.diagnostics[:strand_mode]),
        diagnostics_emission_model = string(result.diagnostics[:emission_model]),
        diagnostics_transition_model = string(result.diagnostics[:transition_model])
    )
end

# --------------------------------------------------------------------------
# Control A — over-correction on UN-CORRUPTED input (specificity / false
# positives). The manuscript's B8 recall table measures only recovery of
# INJECTED errors; it says nothing about whether the corrector also edits
# positions that were already correct. Here we feed the corrector the truth
# itself (zero injected errors) at each assumed error rate and count every edit
# it makes — by construction each such edit is an over-correction, because the
# input had no errors to fix. A corrector that abstains on clean input returns
# an over-correction rate of 0.
# --------------------------------------------------------------------------
function _viterbi_overcorrection_row(
        fixture::ViterbiAccuracyFixture,
        assumed_error_rate::Float64
)::NamedTuple
    # No injection: the observed stream IS the truth.
    converted_observed = _convert_observations(fixture.truth, fixture.domain)
    config = Mycelia.ViterbiCorrectionConfig(error_rate = assumed_error_rate)
    result = Mycelia.correct_observations(
        fixture.graph, [converted_observed]; config = config)
    corrected = [string(observation) for observation in only(result.corrected_observations)]

    Test.@test length(corrected) == length(fixture.truth)

    # Every edit made to already-correct input is a false positive.
    over_correction_edit_distance = _sum_edit_distance(corrected, fixture.truth)
    changed_observations = count(
        index -> corrected[index] != fixture.truth[index], eachindex(fixture.truth))
    correct_positions = sum(length, fixture.truth)
    over_correction_rate = over_correction_edit_distance / correct_positions

    return (
        dataset_id = fixture.dataset_id,
        dataset_name = fixture.dataset_name,
        domain = string(fixture.domain),
        source = fixture.source,
        assumed_error_rate = assumed_error_rate,
        observation_count = length(fixture.truth),
        correct_positions = correct_positions,
        injected_error_count = 0,
        control = "over_correction_uncorrupted",
        corrected_method = "generalized_viterbi_correct_observations",
        changed_observations = changed_observations,
        over_correction_edit_distance = over_correction_edit_distance,
        over_correction_rate = over_correction_rate
    )
end

# --------------------------------------------------------------------------
# Control B — random / shuffled-evidence graph nulls. The manuscript asks how
# much of the perfect recovery is attributable to the graph MODEL rather than to
# the small, low-ambiguity fixtures. We run the SAME injected-error observations
# against two degraded graphs and compare recovery to the real graph:
#
#   (1) shuffled-WEIGHT null — same topology, edge weights permuted. Isolates the
#       learned edge EVIDENCE: if recovery is unchanged, the weights do no work
#       (the decode path is topologically forced on these near-deterministic
#       fixtures); if it collapses, the evidence-weighting drives correction.
#   (2) random-REWIRE null — same vertices + edge count + weight multiset, but
#       edges connect random vertex pairs. Destroys the true k-mer adjacency, so
#       it CAN fail: if recovery survives, the corrector is merely echoing the
#       observations; if it collapses, correction genuinely depends on the graph
#       topology. This is the discriminating null the weight-shuffle cannot be.
# --------------------------------------------------------------------------
function _shuffle_weighted_graph_weights(weighted, seed::Integer)
    edge_labels = collect(MetaGraphsNext.edge_labels(weighted))
    weights = [weighted[src, dst].weight for (src, dst) in edge_labels]
    rng = Random.MersenneTwister(seed)
    shuffled = Random.shuffle(rng, weights)
    for (index, (src, dst)) in enumerate(edge_labels)
        existing = weighted[src, dst]
        weighted[src, dst] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            shuffled[index], existing.src_strand, existing.dst_strand)
    end
    return weighted
end

function _random_rewired_weighted_graph(weighted, seed::Integer)
    labels = collect(MetaGraphsNext.labels(weighted))
    edges = collect(MetaGraphsNext.edge_labels(weighted))
    weights = [weighted[src, dst].weight for (src, dst) in edges]
    rng = Random.MersenneTwister(seed)
    shuffled_weights = Random.shuffle(rng, weights)
    label_type = isempty(labels) ? String : typeof(first(labels))

    rewired = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type = label_type,
        vertex_data_type = Any,
        edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
        weight_function = Mycelia.Rhizomorph.edge_data_weight,
        default_weight = 0.0
    )
    for label in labels
        rewired[label] = nothing
    end

    n_labels = length(labels)
    target_edges = length(edges)
    placed = Set{Tuple{eltype(labels), eltype(labels)}}()
    added = 0
    attempts = 0
    max_attempts = 100 * max(target_edges, 1)
    while added < target_edges && attempts < max_attempts && n_labels > 1
        attempts += 1
        src = labels[rand(rng, 1:n_labels)]
        dst = labels[rand(rng, 1:n_labels)]
        (src == dst || (src, dst) in placed) && continue
        push!(placed, (src, dst))
        rewired[src, dst] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            shuffled_weights[added + 1], Mycelia.Rhizomorph.Forward,
            Mycelia.Rhizomorph.Forward)
        added += 1
    end
    return rewired
end

function _viterbi_null_row(
        fixture::ViterbiAccuracyFixture,
        target_error_rate::Float64
)::NamedTuple
    # Identical injection to the real-graph row (deterministic injector), so the
    # only difference between this row and the recall table is the shuffled graph.
    observed,
    injected_positions = _inject_fixture_errors(
        fixture.truth, fixture.alphabet, target_error_rate
    )
    converted_observed = _convert_observations(observed, fixture.domain)
    config = Mycelia.ViterbiCorrectionConfig(error_rate = target_error_rate)

    injected_error_count = length(injected_positions)
    baseline_edit_distance = _sum_edit_distance(observed, fixture.truth)

    # Decode the SAME observations against three graphs and score each the same
    # way. Only the graph differs across arms.
    function _recovery(weighted_graph)
        result = weighted_graph === nothing ?
                 Mycelia.correct_observations(fixture.graph, [converted_observed]; config = config) :
                 Mycelia.correct_observations(
            fixture.graph, [converted_observed]; config = config,
            weighted_graph = weighted_graph)
        raw = only(result.corrected_observations)
        # A degraded (e.g. randomly-rewired) graph can leave the observation
        # undecodable — the corrector returns `nothing`, or a path of a different
        # length that cannot be positionally scored. That IS the null collapsing:
        # score it as "no correction applied" by falling back to the still-
        # erroneous observed sequence (=> recall 0, edit-distance-reduction 0).
        decoded_seq = raw === nothing ? nothing : [string(o) for o in raw]
        decoded = decoded_seq !== nothing && length(decoded_seq) == length(fixture.truth)
        scored = decoded ? decoded_seq : observed
        recovered = _recovered_error_count(scored, fixture.truth, injected_positions)
        corrected_edit_distance = _sum_edit_distance(scored, fixture.truth)
        recall = injected_error_count == 0 ? 0.0 : recovered / injected_error_count
        reduction = baseline_edit_distance == 0 ? 0.0 :
                    (baseline_edit_distance - corrected_edit_distance) /
                    baseline_edit_distance
        return (recall = recall, reduction = reduction, decoded = decoded)
    end

    real = _recovery(nothing)
    # Shuffled-weight null: same topology, permuted edge weights.
    weight_null = _recovery(
        _shuffle_weighted_graph_weights(
        Mycelia.build_correction_weighted_graph(fixture.graph; config = config),
        VITERBI_NULL_SEED))
    # Random-rewire null: same vertices/edge-count/weight-multiset, random topology.
    rewire_null = _recovery(
        _random_rewired_weighted_graph(
        Mycelia.build_correction_weighted_graph(fixture.graph; config = config),
        VITERBI_NULL_SEED))

    return (
        dataset_id = fixture.dataset_id,
        dataset_name = fixture.dataset_name,
        domain = string(fixture.domain),
        source = fixture.source,
        target_error_rate = target_error_rate,
        injected_error_count = injected_error_count,
        baseline_edit_distance = baseline_edit_distance,
        null_seed = VITERBI_NULL_SEED,
        real_injected_error_recall = real.recall,
        weight_null_injected_error_recall = weight_null.recall,
        rewire_null_injected_error_recall = rewire_null.recall,
        real_edit_distance_reduction = real.reduction,
        weight_null_edit_distance_reduction = weight_null.reduction,
        rewire_null_edit_distance_reduction = rewire_null.reduction,
        real_decoded = real.decoded,
        weight_null_decoded = weight_null.decoded,
        rewire_null_decoded = rewire_null.decoded,
        recall_gap_real_minus_weight_null = real.recall - weight_null.recall,
        recall_gap_real_minus_rewire_null = real.recall - rewire_null.recall
    )
end

function _read_fasta_sequence(path::AbstractString)::String
    lines = readlines(path)
    sequence_lines = [strip(line)
                      for line in lines if !isempty(line) && !startswith(line, ">")]
    return uppercase(join(sequence_lines))
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
        target_error_rate::Float64
)::Tuple{Vector{String}, Vector{Tuple{Int, Int}}}
    candidates = Tuple{Int, Int}[]
    for item_index in 2:(length(truth) - 1)
        for position in eachindex(truth[item_index])
            push!(candidates, (item_index, position))
        end
    end

    target_errors = max(1, round(Int, target_error_rate * length(candidates)))
    target_errors = min(target_errors, length(candidates))
    selected_offsets = round.(Int, range(1, length(candidates); length = target_errors))
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
    return sum(_edit_distance(left_item, right_item)
    for (left_item, right_item) in zip(left, right))
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

function _write_viterbi_control_figure(
        overcorrection::DataFrames.DataFrame,
        null_control::DataFrames.DataFrame,
        plots_dir::AbstractString
)::NamedTuple
    mkpath(plots_dir)
    fig = CairoMakie.Figure(size = (1100, 500), fontsize = 15)

    # Panel 1 — over-correction rate on un-corrupted input (specificity).
    axis1 = CairoMakie.Axis(
        fig[1, 1],
        title = "Control A: over-correction on un-corrupted input",
        xlabel = "Assumed error rate (%)",
        ylabel = "Over-correction rate (edits / correct position)"
    )
    palette = [:seagreen4, :dodgerblue3, :darkorange3]
    grouped_oc = collect(DataFrames.groupby(overcorrection, :dataset_id))
    for (index, group) in enumerate(grouped_oc)
        sorted_group = sort(group, :assumed_error_rate)
        color = palette[mod1(index, length(palette))]
        CairoMakie.scatterlines!(
            axis1,
            100 .* sorted_group.assumed_error_rate,
            sorted_group.over_correction_rate;
            color = color, linewidth = 3, markersize = 12,
            label = first(sorted_group.dataset_name)
        )
    end
    CairoMakie.ylims!(axis1, -0.02, max(0.1, 1.05 *
                                             maximum(overcorrection.over_correction_rate)))

    # Panel 2 — real graph vs the two nulls (Control B). Solid = real graph,
    # dashed = shuffled-weight null, dotted = random-rewire null.
    axis2 = CairoMakie.Axis(
        fig[1, 2],
        title = "Control B: real vs shuffled-weight + random-rewire nulls",
        xlabel = "Injected error rate (%)",
        ylabel = "Injected-error recall",
        yticks = 0.0:0.25:1.0
    )
    grouped_null = collect(DataFrames.groupby(null_control, :dataset_id))
    real_elements = CairoMakie.LineElement[]
    real_labels = String[]
    for (index, group) in enumerate(grouped_null)
        sorted_group = sort(group, :target_error_rate)
        color = palette[mod1(index, length(palette))]
        x_values = 100 .* sorted_group.target_error_rate
        CairoMakie.lines!(
            axis2, x_values, sorted_group.real_injected_error_recall;
            color = color, linewidth = 3)
        CairoMakie.scatter!(
            axis2, x_values, sorted_group.real_injected_error_recall;
            color = color, markersize = 12)
        CairoMakie.lines!(
            axis2, x_values, sorted_group.weight_null_injected_error_recall;
            color = color, linewidth = 2, linestyle = :dash)
        CairoMakie.lines!(
            axis2, x_values, sorted_group.rewire_null_injected_error_recall;
            color = color, linewidth = 2, linestyle = :dot)
        CairoMakie.scatter!(
            axis2, x_values, sorted_group.rewire_null_injected_error_recall;
            color = color, markersize = 10, marker = :xcross)
        push!(real_elements, CairoMakie.LineElement(color = color, linewidth = 3))
        push!(real_labels, first(sorted_group.dataset_name))
    end
    style_elements = [
        CairoMakie.LineElement(color = :gray30, linewidth = 3),
        CairoMakie.LineElement(color = :gray30, linewidth = 2, linestyle = :dash),
        CairoMakie.LineElement(color = :gray30, linewidth = 2, linestyle = :dot)
    ]
    CairoMakie.ylims!(axis2, -0.05, 1.05)
    CairoMakie.axislegend(axis1, position = :lt)
    CairoMakie.axislegend(
        axis2, vcat(real_elements, style_elements),
        vcat(real_labels, ["real graph", "shuffled-weight null", "random-rewire null"]);
        position = :rc)

    png_path = joinpath(plots_dir, "viterbi_control_overcorrection_and_null.png")
    svg_path = joinpath(plots_dir, "viterbi_control_overcorrection_and_null.svg")
    CairoMakie.save(png_path, fig)
    CairoMakie.save(svg_path, fig)
    return (png = png_path, svg = svg_path)
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
