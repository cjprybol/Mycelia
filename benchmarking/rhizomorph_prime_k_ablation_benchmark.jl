# Prime-vs-composite k ablation for the Rhizomorph 2-stage corrector on
# REPEAT-RICH fixtures (bead td-tjym).
#
# Motivation (PR #397): snapping the corrector's k-ladder to primes and rerunning
# the B8 accuracy benchmark at prime k=11 showed that on the existing
# NON-repetitive fixtures (viroid / phiX / text), saturated recall was ROBUST to
# primality — those fixtures cannot exercise the prime-vs-composite distinction.
#
# A composite k that shares a factor with a tandem-repeat period p aliases the
# repeat onto SELF-PERIODIC k-mers: the k-mer ATGATGATG (k=9 over a period-3
# repeat) is itself internally periodic, so its phase within the repeat is
# ambiguous and the path decode can settle on the wrong copy/phase. A prime k
# (11, 19, 23) shares no factor with period 3, so a k-mer covering the same region
# spans a non-integer number of periods and the phase degeneracy is broken.
#
# EXPERIMENT DESIGN (differs from the token-level B8 benchmark on purpose):
#   - Fixtures embed a repeat CORE inside UNIQUE flanks so the k-mer graph
#     genuinely branches (a pure tandem repeat is a non-branching cycle whose
#     correction is trivial at every k — that is why token-level injection
#     saturates regardless of primality).
#   - Errors are injected ONCE at the SEQUENCE level, in the interior with a clean
#     margin >= max(k), so the SAME corrupted sequence is decoded at every k — a
#     fair prime-vs-composite comparison.
#   - The corrupted sequence is windowed into overlapping k-mers (each sequence
#     error therefore spans k consecutive observations), corrected as a PATH, then
#     reconstructed back to a sequence and scored at the sequence level. This
#     surfaces both recovery (recall) AND over-correction (edit-distance-reduction
#     below recall, or negative, when the corrector edits already-correct bases).
#
# SMOKE scale. Usage:
#   julia --project=. benchmarking/rhizomorph_prime_k_ablation_benchmark.jl
#   julia --project=. benchmarking/rhizomorph_prime_k_ablation_benchmark.jl --output-dir /tmp/primek --skip-plots

import BioSequences
import CairoMakie
import DataFrames
import FASTX
import Kmers
import Mycelia
import Statistics
import Test

include(joinpath(@__DIR__, "benchmark_artifacts.jl"))

# Prime k share no factor with periods 2 or 3; composite k are all divisible by 3
# (9 = 3x3, 15 = 3x5, 21 = 3x7) so they alias the period-3 tandem fixture. All k
# are odd, avoiding any even-k edge cases in graph construction.
const PRIME_K = (11, 19, 23)
const COMPOSITE_K = (9, 15, 21)
const ABLATION_K = sort(collect(Iterators.flatten((PRIME_K, COMPOSITE_K))))
# Clean margin (bp) left un-corrupted at each end so every k has clean flanking
# k-mers; must be >= max(ABLATION_K).
const ABLATION_CLEAN_MARGIN = maximum(ABLATION_K)
const ABLATION_ERROR_RATES = (0.02, 0.05, 0.10)
# Error injection modes. :isolated places single evenly-spaced substitutions;
# :burst places adjacent substitution clusters (length ABLATION_BURST_LENGTH),
# which is the regime that can convert a repeat k-mer into a different-phase repeat
# k-mer and thereby expose composite-k phase aliasing.
const ABLATION_ERROR_MODES = (:isolated, :burst)
const ABLATION_BURST_LENGTH = 3
const ABLATION_ALPHABET = collect("ACGT")
const ABLATION_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "rhizomorph_prime_k_ablation"
)

struct AblationFixture
    dataset_id::String
    dataset_name::String
    sequence::String
    period_note::String
end

function main(args::Vector{String} = ARGS)::Nothing
    output_dir = _ablation_arg_value(args, "--output-dir", ABLATION_DEFAULT_OUTPUT_DIR)
    write_plots = !("--skip-plots" in args)
    artifacts = run_prime_k_ablation_benchmark(output_dir; write_plots = write_plots)
    println("Wrote Rhizomorph prime-vs-composite k ablation artifacts:")
    println("  root: $(artifacts.root)")
    println("  per_run_csv: $(artifacts.per_run_csv)")
    println("  delta_csv: $(artifacts.delta_csv)")
    println("  figure_png: $(artifacts.figure_png)")
    println("  figure_svg: $(artifacts.figure_svg)")
    println()
    println("Prime-vs-composite delta per fixture x error-mode (mean over k and error rate):")
    for row in eachrow(artifacts.delta_table)
        println("  $(rpad(row.dataset_id, 24)) $(rpad(row.error_mode, 9)) " *
                "prime_recall=$(round(row.prime_mean_recall; digits = 4)) " *
                "composite_recall=$(round(row.composite_mean_recall; digits = 4)) " *
                "recall_delta=$(round(row.recall_delta_prime_minus_composite; digits = 4)) " *
                "edr_delta=$(round(row.edit_distance_reduction_delta; digits = 4))")
    end
    return nothing
end

function run_prime_k_ablation_benchmark(
        output_dir::AbstractString = ABLATION_DEFAULT_OUTPUT_DIR;
        write_plots::Bool = true
)::NamedTuple
    fixtures = ablation_fixtures()
    rows = NamedTuple[]
    for fixture in fixtures
        for k in ABLATION_K
            for error_rate in ABLATION_ERROR_RATES
                for error_mode in ABLATION_ERROR_MODES
                    push!(rows, _ablation_row(fixture, k, error_rate, error_mode))
                end
            end
        end
    end
    per_run = DataFrames.DataFrame(rows)
    delta = _ablation_delta_table(per_run)

    artifacts = write_benchmark_artifacts(
        [
            "prime_k_ablation_per_run" => per_run,
            "prime_k_ablation_delta" => delta
        ];
        output_dir = output_dir,
        run_id = "prime_k_ablation_local_20260711",
        scale = "local-smoke",
        dataset_ids = [fixture.dataset_id for fixture in fixtures],
        command_args = [
            "julia", "--project=.",
            "benchmarking/rhizomorph_prime_k_ablation_benchmark.jl"],
        metadata = Dict(
            "bead" => "td-tjym",
            "benchmark" => "prime_vs_composite_k_ablation_repeat_rich",
            "baseline" => "uncorrected sequence-level injected-error observations",
            "prime_k" => collect(PRIME_K),
            "composite_k" => collect(COMPOSITE_K),
            "error_rates" => collect(ABLATION_ERROR_RATES),
            "clean_margin_bp" => ABLATION_CLEAN_MARGIN,
            "injection" => "sequence-level, interior, k-independent (same corrupted sequence decoded at every k)",
            "scoring" => "reconstruct corrected overlapping k-mers to a sequence; recall + edit-distance-reduction at sequence level",
            "hypothesis" => "composite k sharing a factor with a tandem period alias the repeat onto self-periodic k-mers, degrading path decode; prime k break the periodicity"
        ),
        table_context_columns = Dict(
            "prime_k_ablation_per_run" => Dict("benchmark_dataset_id" => "dataset_id"),
            "prime_k_ablation_delta" => Dict("benchmark_dataset_id" => "dataset_id")
        )
    )

    figure_png = ""
    figure_svg = ""
    if write_plots
        figure_paths = _write_ablation_figure(per_run, artifacts.layout.plots)
        figure_png = figure_paths.png
        figure_svg = figure_paths.svg
    end

    return (
        root = artifacts.root,
        per_run_csv = artifacts.tables["prime_k_ablation_per_run"].table,
        delta_csv = artifacts.tables["prime_k_ablation_delta"].table,
        index = artifacts.index,
        provenance = artifacts.provenance,
        figure_png = figure_png,
        figure_svg = figure_svg,
        per_run_table = per_run,
        delta_table = delta,
        per_run_rows = DataFrames.nrow(per_run),
        delta_rows = DataFrames.nrow(delta)
    )
end

# Each fixture embeds a repeat core inside two unique 30 bp flanks so the k-mer
# graph branches into and out of the repeat.
function ablation_fixtures()::Vector{AblationFixture}
    prefix = _unique_flank(1, 30)
    suffix = _unique_flank(2, 30)
    return AblationFixture[
        AblationFixture(
            "tandem_period3_atg",
            "Period-3 tandem (ATG) core in unique flanks",
            prefix * repeat("ATG", 100) * suffix,
            "period 3 — shared factor with every composite k (9,15,21); coprime to primes"
        ),
        AblationFixture(
            "tandem_period2_ga",
            "Period-2 tandem (GA) core in unique flanks",
            prefix * repeat("GA", 150) * suffix,
            "period 2 — coprime to both k-sets (odd k); contrast fixture"
        ),
        AblationFixture(
            "palindrome_rich",
            "Palindrome-rich inverted-repeat core in unique flanks",
            prefix * _palindrome_rich_core() * suffix,
            "repeated reverse-complement inverted-repeat unit (period 60)"
        ),
        AblationFixture(
            "control_nonrepetitive",
            "Non-repetitive control (deterministic pseudo-random)",
            _nonrepetitive_sequence(360, 3),
            "no dominant period — primality expected to be neutral"
        )
    ]
end

# A palindrome-rich core: arm + reverse_complement(arm) forms a 60 bp DNA
# palindrome; concatenating copies yields dense reverse-complement self-similarity.
function _palindrome_rich_core()::String
    arm = "ACGGTACCTTGACATGCACGTTGGATCCAT"  # 30 bp
    rc = string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(arm)))
    unit = arm * rc  # 60 bp palindrome
    return repeat(unit, 5)  # 300 bp
end

# Unique (non-repetitive) flank / control via a deterministic linear-congruential
# walk over the alphabet — no RNG dependency, byte-stable across runs.
function _unique_flank(seed_index::Int, length_bp::Int)::String
    return _nonrepetitive_sequence(length_bp, seed_index)
end

function _nonrepetitive_sequence(length_bp::Int, seed_index::Int)::String
    alphabet = ABLATION_ALPHABET
    characters = Vector{Char}(undef, length_bp)
    state = (2246822519 + seed_index * 40503) % 2147483648  # distinct odd-ish seed
    for index in 1:length_bp
        state = (1103515245 * state + 12345) % 2147483648
        characters[index] = alphabet[(state % 4) + 1]
    end
    return String(characters)
end

function _ablation_row(
        fixture::AblationFixture,
        k::Int,
        error_rate::Float64,
        error_mode::Symbol
)::NamedTuple
    clean = fixture.sequence
    corrupted,
    injected_positions = _inject_sequence_errors(
        clean, ABLATION_ALPHABET, error_rate, ABLATION_CLEAN_MARGIN, error_mode)

    record = FASTX.FASTA.Record(fixture.dataset_id, BioSequences.LongDNA{4}(clean))
    graph = Mycelia.Rhizomorph.build_kmer_graph(
        [record], k;
        dataset_id = "prime_k_ablation_$(fixture.dataset_id)_k$(k)",
        mode = :singlestrand
    )

    observed_kmers = _ablation_kmers(corrupted, k)
    converted = [Kmers.DNAKmer{k}(observation) for observation in observed_kmers]
    config = Mycelia.ViterbiCorrectionConfig(error_rate = error_rate)
    result = Mycelia.correct_observations(graph, [converted]; config = config)
    corrected_kmers = [string(o) for o in only(result.corrected_observations)]

    Test.@test length(corrected_kmers) == length(observed_kmers)

    reconstructed = _reconstruct_sequence(corrected_kmers, k)
    # Reconstruction preserves length; scored positionally (Hamming) at seq level.
    baseline_mismatch = _hamming(corrupted, clean)
    corrected_mismatch = _hamming(reconstructed, clean)
    injected_error_count = length(injected_positions)
    recovered = count(pos -> _char_at(reconstructed, pos) == clean[pos], injected_positions)

    recall = injected_error_count == 0 ? 0.0 : recovered / injected_error_count
    edit_distance_reduction = baseline_mismatch == 0 ? 0.0 :
                              (baseline_mismatch - corrected_mismatch) / baseline_mismatch
    # Edits the corrector made at positions that were NOT injected (over-correction).
    over_correction = count(
        pos -> !(pos in injected_positions) && _char_at(reconstructed, pos) != clean[pos],
        eachindex(clean))

    return (
        dataset_id = fixture.dataset_id,
        dataset_name = fixture.dataset_name,
        period_note = fixture.period_note,
        k = k,
        k_class = (k in PRIME_K) ? "prime" : "composite",
        error_mode = string(error_mode),
        target_error_rate = error_rate,
        sequence_length = length(clean),
        kmer_count = length(observed_kmers),
        injected_error_count = injected_error_count,
        baseline_method = "uncorrected_observations",
        corrected_method = "rhizomorph_correct_observations",
        baseline_mismatch = baseline_mismatch,
        corrected_mismatch = corrected_mismatch,
        edit_distance_reduction = edit_distance_reduction,
        injected_error_recall = recall,
        over_correction_positions = over_correction,
        reconstructed_length = length(reconstructed),
        corrected_all_injected_positions = recovered == injected_error_count
    )
end

# Collapse the per-run table to one prime-advantage row per fixture: mean recall
# and mean edit-distance-reduction over prime k vs composite k (across all error
# rates), and their deltas (prime minus composite).
function _ablation_delta_table(per_run::DataFrames.DataFrame)::DataFrames.DataFrame
    rows = NamedTuple[]
    for group in DataFrames.groupby(per_run, [:dataset_id, :error_mode])
        prime = group[group.k_class .== "prime", :]
        composite = group[group.k_class .== "composite", :]
        prime_recall = Statistics.mean(prime.injected_error_recall)
        composite_recall = Statistics.mean(composite.injected_error_recall)
        prime_edr = Statistics.mean(prime.edit_distance_reduction)
        composite_edr = Statistics.mean(composite.edit_distance_reduction)
        prime_oc = Statistics.mean(prime.over_correction_positions)
        composite_oc = Statistics.mean(composite.over_correction_positions)
        push!(rows,
            (
                dataset_id = first(group.dataset_id),
                dataset_name = first(group.dataset_name),
                error_mode = first(group.error_mode),
                period_note = first(group.period_note),
                prime_k = join(string.(PRIME_K), ","),
                composite_k = join(string.(COMPOSITE_K), ","),
                prime_mean_recall = prime_recall,
                composite_mean_recall = composite_recall,
                recall_delta_prime_minus_composite = prime_recall - composite_recall,
                prime_mean_edit_distance_reduction = prime_edr,
                composite_mean_edit_distance_reduction = composite_edr,
                edit_distance_reduction_delta = prime_edr - composite_edr,
                prime_mean_over_correction = prime_oc,
                composite_mean_over_correction = composite_oc
            ))
    end
    return DataFrames.DataFrame(rows)
end

function _ablation_kmers(sequence::AbstractString, k::Int)::Vector{String}
    chars = collect(sequence)
    n = length(chars)
    return [String(chars[index:(index + k - 1)]) for index in 1:(n - k + 1)]
end

function _char_at(sequence::AbstractString, index::Int)::Char
    return collect(sequence)[index]
end

# Reconstruct a sequence from a path of overlapping k-mers: first k-mer in full,
# then the last character of each subsequent k-mer. Length is preserved so long as
# the corrector kept the observation count (asserted by the caller). Any path
# inconsistency surfaces as residual mismatches — exactly what we want to score.
function _reconstruct_sequence(kmers::Vector{String}, k::Int)::String
    isempty(kmers) && return ""
    buffer = collect(first(kmers))
    for kmer in kmers[2:end]
        push!(buffer, last(kmer))
    end
    return String(buffer)
end

function _hamming(left::AbstractString, right::AbstractString)::Int
    left_chars = collect(left)
    right_chars = collect(right)
    n = min(length(left_chars), length(right_chars))
    mismatches = abs(length(left_chars) - length(right_chars))
    for index in 1:n
        mismatches += left_chars[index] == right_chars[index] ? 0 : 1
    end
    return mismatches
end

# Deterministic substitution injector at the SEQUENCE level (no RNG). Errors are
# placed in the interior leaving `clean_margin` bp untouched at each end.
# k-independent: the same corrupted sequence is decoded at every k.
#   :isolated — single substitutions at evenly-spaced positions.
#   :burst    — clusters of ABLATION_BURST_LENGTH adjacent substitutions at
#               evenly-spaced anchors (fewer anchors, same total error budget).
function _inject_sequence_errors(
        clean::AbstractString,
        alphabet::Vector{Char},
        target_error_rate::Float64,
        clean_margin::Int,
        error_mode::Symbol
)::Tuple{String, Vector{Int}}
    chars = collect(clean)
    n = length(chars)
    lo = clean_margin + 1
    hi = n - clean_margin
    candidates = collect(lo:hi)
    isempty(candidates) && return String(chars), Int[]

    target_errors = max(1, round(Int, target_error_rate * length(candidates)))
    target_errors = min(target_errors, length(candidates))

    positions = Int[]
    if error_mode == :isolated
        offsets = round.(Int, range(1, length(candidates); length = target_errors))
        positions = unique(candidates[offsets])
    elseif error_mode == :burst
        n_anchors = max(1, cld(target_errors, ABLATION_BURST_LENGTH))
        anchor_offsets = round.(Int, range(1, length(candidates); length = n_anchors))
        anchors = unique(candidates[anchor_offsets])
        position_set = Int[]
        for anchor in anchors
            for delta in 0:(ABLATION_BURST_LENGTH - 1)
                position = anchor + delta
                position <= hi && push!(position_set, position)
            end
        end
        positions = sort(unique(position_set))
    else
        throw(ArgumentError("unsupported error mode: $error_mode"))
    end

    for position in positions
        chars[position] = _ablation_replacement_char(chars[position], alphabet, position)
    end
    return String(chars), positions
end

function _ablation_replacement_char(character::Char, alphabet::Vector{Char}, offset::Int)::Char
    for shift in 0:(length(alphabet) - 1)
        candidate = alphabet[mod1(offset + shift, length(alphabet))]
        if candidate != character
            return candidate
        end
    end
    error("no replacement character available for $character")
end

function _write_ablation_figure(
        per_run::DataFrames.DataFrame,
        plots_dir::AbstractString
)::NamedTuple
    mkpath(plots_dir)
    fixtures = unique(per_run.dataset_id)
    fig = CairoMakie.Figure(size = (1100, 850), fontsize = 15)
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

    prime_color = :dodgerblue3
    composite_color = :darkorange3

    # Marker shape encodes error mode; color encodes prime vs composite k.
    mode_markers = Dict("isolated" => :circle, "burst" => :xcross)

    for (fixture_index, dataset_id) in enumerate(fixtures)
        row, col = positions[mod1(fixture_index, length(positions))]
        subset = per_run[per_run.dataset_id .== dataset_id, :]
        name = first(subset.dataset_name)
        axis = CairoMakie.Axis(
            fig[row, col],
            title = name,
            xlabel = "k (correction k-mer length)",
            ylabel = "Edit-distance reduction (mean over error rate)",
            xticks = ABLATION_K,
            yticks = 0.0:0.25:1.0
        )
        for mode in string.(ABLATION_ERROR_MODES)
            mode_subset = subset[subset.error_mode .== mode, :]
            ks = sort(unique(mode_subset.k))
            mean_edr = [Statistics.mean(
                            mode_subset[mode_subset.k .== k, :edit_distance_reduction])
                        for k in ks]
            colors = [(k in PRIME_K) ? prime_color : composite_color for k in ks]
            CairoMakie.lines!(axis, ks, mean_edr; color = :gray70, linewidth = 1.0)
            CairoMakie.scatter!(
                axis, ks, mean_edr; color = colors, markersize = 15,
                marker = mode_markers[mode])
        end
        CairoMakie.ylims!(axis, -0.15, 1.05)
    end

    legend_elements = [
        CairoMakie.MarkerElement(color = prime_color, marker = :circle, markersize = 15),
        CairoMakie.MarkerElement(color = composite_color, marker = :circle, markersize = 15),
        CairoMakie.MarkerElement(color = :gray30, marker = :circle, markersize = 15),
        CairoMakie.MarkerElement(color = :gray30, marker = :xcross, markersize = 15)
    ]
    CairoMakie.Legend(
        fig[3, 1:2], legend_elements,
        ["prime k ($(join(PRIME_K, ", ")))", "composite k ($(join(COMPOSITE_K, ", ")))",
            "isolated errors", "burst errors (len $(ABLATION_BURST_LENGTH))"];
        orientation = :horizontal, framevisible = false, nbanks = 1)
    CairoMakie.Label(
        fig[0, 1:2],
        "Rhizomorph corrector accuracy vs k on repeat-rich fixtures (prime vs composite)";
        fontsize = 18, font = :bold)

    png_path = joinpath(plots_dir, "prime_k_ablation_accuracy_vs_k.png")
    svg_path = joinpath(plots_dir, "prime_k_ablation_accuracy_vs_k.svg")
    CairoMakie.save(png_path, fig)
    CairoMakie.save(svg_path, fig)
    return (png = png_path, svg = svg_path)
end

function _ablation_arg_value(
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
