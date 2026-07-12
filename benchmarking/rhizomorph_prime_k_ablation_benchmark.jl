# Prime-vs-composite k ablation for the Rhizomorph assembler on REPEAT-RICH
# fixtures, in the DE-NOVO regime where k-mer aliasing can actually bite
# (bead td-tjym).
#
# BACKGROUND. PR #397 snapped the corrector's k-ladder to primes. The first
# version of this ablation (reviewed on PR #404) could not show any prime-vs-
# composite difference because it was structurally forced to delta=0: it built
# the k-mer graph from the CLEAN reference (so the graph held only the true path)
# and decoded each observation anchored at both true endpoints inside clean
# margins. That measures "recover the true path when the true path is the only
# path and both ends are pinned" — trivially perfect at every k.
#
# THIS VERSION runs the regime where k-mer aliasing is a real phenomenon: DE-NOVO
# assembly from ERROR-CONTAINING READS. For each fixture we simulate a read set
# (coverage x tiling, substitution errors), build+clean the assembly graph FROM
# THE READS at each k (`Mycelia.Rhizomorph.assemble_genome(reads; k=...)`), and
# measure how much of the reference is reconstructed. Nothing is anchored to
# truth; distinct repeat copies and error k-mers genuinely compete in the graph.
#
# WHAT WE FIND (see the committed tables; numbers are the deliverable, not these
# comments):
#   - On the NATURAL tandem / palindrome fixtures, reference recovery rises
#     monotonically with k-SIZE (larger k spans more context). Controlling for
#     size via adjacent prime/composite pairs (9 vs 11, 15 vs 17, 21 vs 23),
#     there is no additional deficit attributable to primality: a composite k
#     divisible by the repeat period recovers the same as its size-matched prime
#     neighbour. The prime-vs-composite gap on natural fixtures is a k-size
#     artifact, not a primality effect.
#   - POSITIVE CONTROL (engineered `shared_repeat_collision` fixture): two
#     otherwise-unique "genes" share an internal period-3 repeat of length
#     exactly 9. A composite k=9 (=3x3) fits ENTIRELY inside the shared repeat,
#     so the two genes collide on the identical 9-mer and the assembly tangles;
#     the size-matched prime k=11 is forced to span the repeat into the two
#     genes' differing flanks, so no shared k-mer forms and both genes are
#     recovered cleanly. This yields a large, reproducible recovery gap
#     (composite << prime) and is the harness's proof that it CAN register a
#     composite-k aliasing deficit when one exists. It is an ENGINEERED
#     demonstration of the shared-repeat mechanism, not a claim that natural
#     tandem repeats exhibit it.
#
# SMOKE scale, deterministic (fixed read-simulation seeds). Usage:
#   julia --project=. benchmarking/rhizomorph_prime_k_ablation_benchmark.jl
#   julia --project=. benchmarking/rhizomorph_prime_k_ablation_benchmark.jl --output-dir /tmp/primek --skip-plots

import BioSequences
import CairoMakie
import DataFrames
import FASTX
import Mycelia
import Random
import Statistics

include(joinpath(@__DIR__, "benchmark_artifacts.jl"))

# Size-matched prime/composite pairs so a prime-minus-composite delta is not
# confounded by k-size: each composite (divisible by 3) is bracketed just below a
# prime of nearly the same length -> (9,11), (15,17), (21,23).
const PRIME_K = (11, 17, 23)
const COMPOSITE_K = (9, 15, 21)
const SIZE_MATCHED_PAIRS = ((9, 11), (15, 17), (21, 23))  # (composite, prime)
const ABLATION_K = sort(collect(Iterators.flatten((PRIME_K, COMPOSITE_K))))

# Read-simulation parameters (SMOKE, deterministic).
const ABLATION_COVERAGE = 30
const ABLATION_READ_LENGTH = 45
const ABLATION_READ_ERROR_RATE = 0.02
const ABLATION_SEEDS = (1, 2, 3)
const ABLATION_ALPHABET = collect("ACGT")
# Fixed reference-recovery window, INDEPENDENT of the assembly k, so recovery is
# comparable across k values.
const ABLATION_RECOVERY_W = 15
# Positive-control firing bar. The COLLISION PENALTY at composite k=9 (recovery of
# the no-collision variant minus the shared-collision variant, at IDENTICAL
# fixture structure so k-size cancels) must be at least this large AND must
# collapse at prime k=11 (the composite-k aliasing must be resolved by the prime k
# spanning the shared repeat). This isolates aliasing from k-size.
const POSITIVE_CONTROL_MIN_GAP = 0.10
const POSITIVE_CONTROL_SHARED_ID = "shared_repeat_collision"
const POSITIVE_CONTROL_DISTINCT_ID = "distinct_repeat_no_collision"
const ABLATION_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "rhizomorph_prime_k_ablation"
)

struct AblationFixture
    dataset_id::String
    dataset_name::String
    category::String            # "natural", "collision_shared", or "collision_distinct"
    references::Vector{String}  # one entry per distinct source sequence
    note::String
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
    println("Size-matched prime-minus-composite recovery delta per fixture x k-pair:")
    for row in eachrow(artifacts.delta_table)
        println("  $(rpad(row.dataset_id, 28)) pair(c=$(row.composite_k),p=$(row.prime_k)) " *
                "composite_recovery=$(round(row.composite_recovery; digits = 3)) " *
                "prime_recovery=$(round(row.prime_recovery; digits = 3)) " *
                "delta=$(round(row.recovery_delta_prime_minus_composite; digits = 3))")
    end
    println()
    println("POSITIVE CONTROL — collision penalty (recovery of no-collision minus " *
            "shared-collision, structure held identical):")
    for row in eachrow(artifacts.positive_control_summary)
        println("  k=$(row.k) ($(row.k_class))  shared=$(round(row.shared_recovery; digits = 3)) " *
                "distinct=$(round(row.distinct_recovery; digits = 3)) " *
                "collision_penalty=$(round(row.collision_penalty; digits = 3))")
    end
    println()
    println("Positive control fires (k=9 collision penalty >= $(POSITIVE_CONTROL_MIN_GAP) " *
            "AND collapses by >= $(POSITIVE_CONTROL_MIN_GAP) at prime k=11): " *
            "$(artifacts.positive_control_fires)")
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
            push!(rows, _ablation_row(fixture, k))
        end
    end
    per_run = DataFrames.DataFrame(rows)
    delta = _ablation_delta_table(per_run, fixtures)
    positive_control = _positive_control_summary(per_run)

    artifacts = write_benchmark_artifacts(
        [
            "prime_k_ablation_per_run" => per_run,
            "prime_k_ablation_delta" => delta,
            "prime_k_ablation_positive_control" => positive_control
        ];
        output_dir = output_dir,
        run_id = "prime_k_ablation_denovo_local_20260711",
        scale = "local-smoke",
        dataset_ids = [fixture.dataset_id for fixture in fixtures],
        command_args = [
            "julia", "--project=.",
            "benchmarking/rhizomorph_prime_k_ablation_benchmark.jl"],
        metadata = Dict(
            "bead" => "td-tjym",
            "benchmark" => "prime_vs_composite_k_ablation_denovo",
            "regime" => "de-novo assembly from error-containing reads (graph built from reads, not the clean reference; no endpoint anchoring)",
            "prime_k" => collect(PRIME_K),
            "composite_k" => collect(COMPOSITE_K),
            "size_matched_pairs" => [collect(pair) for pair in SIZE_MATCHED_PAIRS],
            "coverage" => ABLATION_COVERAGE,
            "read_length" => ABLATION_READ_LENGTH,
            "read_error_rate" => ABLATION_READ_ERROR_RATE,
            "seeds" => collect(ABLATION_SEEDS),
            "recovery_metric" => "fraction of reference $(ABLATION_RECOVERY_W)-mers present in assembled contigs (either strand), mean over seeds",
            "measures" => "whether prime-vs-composite k changes DE-NOVO reference recovery on repeat-rich fixtures; the positive control confirms the harness detects a composite-k aliasing deficit when one exists",
            "positive_control_min_gap" => POSITIVE_CONTROL_MIN_GAP
        ),
        table_context_columns = Dict(
            "prime_k_ablation_per_run" => Dict("benchmark_dataset_id" => "dataset_id"),
            "prime_k_ablation_delta" => Dict("benchmark_dataset_id" => "dataset_id")
        )
    )
    positive_control_csv = artifacts.tables["prime_k_ablation_positive_control"].table

    figure_png = ""
    figure_svg = ""
    if write_plots
        figure_paths = _write_ablation_figure(per_run, fixtures, artifacts.layout.plots)
        figure_png = figure_paths.png
        figure_svg = figure_paths.svg
    end

    return (
        root = artifacts.root,
        per_run_csv = artifacts.tables["prime_k_ablation_per_run"].table,
        delta_csv = artifacts.tables["prime_k_ablation_delta"].table,
        positive_control_csv = positive_control_csv,
        index = artifacts.index,
        provenance = artifacts.provenance,
        figure_png = figure_png,
        figure_svg = figure_svg,
        per_run_table = per_run,
        delta_table = delta,
        positive_control_summary = positive_control,
        per_run_rows = DataFrames.nrow(per_run),
        delta_rows = DataFrames.nrow(delta),
        positive_control_fires = _positive_control_fires(per_run)
    )
end

function ablation_fixtures()::Vector{AblationFixture}
    flank_left = _nonrepetitive_sequence(30, 1)
    flank_right = _nonrepetitive_sequence(30, 2)
    return AblationFixture[
        AblationFixture(
            "tandem_period3_atg",
            "Period-3 tandem (ATG) core in unique flanks",
            "natural",
            [flank_left * repeat("ATG", 60) * flank_right],
            "period 3 — divisible by every composite k (9,15,21)"
        ),
        AblationFixture(
            "tandem_period2_ga",
            "Period-2 tandem (GA) core in unique flanks",
            "natural",
            [flank_left * repeat("GA", 90) * flank_right],
            "period 2 — coprime to both odd k-sets"
        ),
        AblationFixture(
            "palindrome_rich",
            "Palindrome-rich inverted-repeat core in unique flanks",
            "natural",
            [flank_left * _palindrome_rich_core() * flank_right],
            "repeated reverse-complement inverted-repeat unit (period 60)"
        ),
        AblationFixture(
            "control_nonrepetitive",
            "Non-repetitive control (deterministic pseudo-random)",
            "natural",
            [_nonrepetitive_sequence(220, 3)],
            "no dominant period — primality expected to be neutral"
        ),
        # POSITIVE CONTROL (shared). Two unique genes sharing an internal period-3
        # repeat of length exactly 9. Composite k=9 fits inside the shared repeat
        # (genes collide on the identical 9-mer -> tangle); prime k=11 must span
        # the repeat into the genes' DIFFERING flanks (no shared k-mer -> clean).
        AblationFixture(
            POSITIVE_CONTROL_SHARED_ID,
            "Positive control (shared): two genes sharing a 9 bp period-3 repeat",
            "collision_shared",
            [
                _nonrepetitive_sequence(20, 4) * repeat("ATG", 3) *
                _nonrepetitive_sequence(20, 5),
                _nonrepetitive_sequence(20, 6) * repeat("ATG", 3) *
                _nonrepetitive_sequence(20, 7)
            ],
            "engineered 9 bp SHARED repeat: composite k=9 collides, prime k=11 spans into unique flanks"
        ),
        # NEGATIVE CONTROL (distinct). Byte-for-byte the same structure as the
        # shared fixture — same flanks, same 9 bp period-3 insert length — except
        # the two genes carry DIFFERENT period-3 units (ATGx3 vs GCAx3), so there
        # is no shared 9-mer and no collision. Any k=9 recovery gap between this
        # and the shared fixture is attributable to the COLLISION, not to k-size
        # (structure and length are identical across the two).
        AblationFixture(
            POSITIVE_CONTROL_DISTINCT_ID,
            "Negative control (distinct): same structure, non-shared 9 bp repeats",
            "collision_distinct",
            [
                _nonrepetitive_sequence(20, 4) * repeat("ATG", 3) *
                _nonrepetitive_sequence(20, 5),
                _nonrepetitive_sequence(20, 6) * repeat("GCA", 3) *
                _nonrepetitive_sequence(20, 7)
            ],
            "distinct 9 bp inserts (ATGx3 vs GCAx3): no shared k-mer; k-size-matched control for the collision fixture"
        )
    ]
end

# A palindrome-rich core: arm + reverse_complement(arm) forms a 60 bp DNA
# palindrome; concatenating copies yields dense reverse-complement self-similarity.
function _palindrome_rich_core()::String
    arm = "ACGGTACCTTGACATGCACGTTGGATCCAT"  # 30 bp
    rc = string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(arm)))
    return repeat(arm * rc, 5)  # 300 bp, period 60
end

# Deterministic pseudo-random (non-repetitive) DNA via a linear-congruential walk
# — no RNG dependency, byte-stable across runs. Distinct seed_index -> distinct
# sequence, so flanks and control differ.
function _nonrepetitive_sequence(length_bp::Int, seed_index::Int)::String
    characters = Vector{Char}(undef, length_bp)
    state = (2246822519 + seed_index * 40503) % 2147483648
    for index in 1:length_bp
        state = (1103515245 * state + 12345) % 2147483648
        characters[index] = ABLATION_ALPHABET[(state % 4) + 1]
    end
    return String(characters)
end

function _ablation_row(fixture::AblationFixture, k::Int)::NamedTuple
    recoveries = Float64[]
    contig_counts = Int[]
    assembled_flags = Bool[]
    for seed in ABLATION_SEEDS
        reads = _simulate_reads(fixture.references, seed)
        recovery, n_contigs, assembled = _assemble_and_score(reads, fixture.references, k)
        push!(recoveries, recovery)
        push!(contig_counts, n_contigs)
        push!(assembled_flags, assembled)
    end
    return (
        dataset_id = fixture.dataset_id,
        dataset_name = fixture.dataset_name,
        category = fixture.category,
        note = fixture.note,
        k = k,
        k_class = (k in PRIME_K) ? "prime" : "composite",
        n_references = length(fixture.references),
        reference_length = sum(length, fixture.references),
        coverage = ABLATION_COVERAGE,
        read_length = ABLATION_READ_LENGTH,
        read_error_rate = ABLATION_READ_ERROR_RATE,
        n_seeds = length(ABLATION_SEEDS),
        recovery_w = ABLATION_RECOVERY_W,
        mean_reference_recovery = Statistics.mean(recoveries),
        min_reference_recovery = minimum(recoveries),
        max_reference_recovery = maximum(recoveries),
        mean_contig_count = Statistics.mean(contig_counts),
        all_seeds_assembled = all(assembled_flags)
    )
end

# Simulate a read set tiling every reference sequence at ABLATION_COVERAGE, with
# per-base substitution errors. Deterministic given (references, seed).
function _simulate_reads(references::Vector{String}, seed::Int)::Vector{FASTX.FASTA.Record}
    rng = Random.MersenneTwister(seed)
    reads = FASTX.FASTA.Record[]
    read_index = 0
    for reference in references
        n = length(reference)
        read_length = min(ABLATION_READ_LENGTH, n)
        n_reads = max(1, cld(ABLATION_COVERAGE * n, read_length))
        for _ in 1:n_reads
            start = n == read_length ? 1 : rand(rng, 1:(n - read_length + 1))
            characters = collect(reference[start:(start + read_length - 1)])
            for position in eachindex(characters)
                if rand(rng) < ABLATION_READ_ERROR_RATE
                    original = characters[position]
                    characters[position] = rand(
                        rng, filter(candidate -> candidate != original, ABLATION_ALPHABET))
                end
            end
            read_index += 1
            push!(reads, FASTX.FASTA.Record(
                "r$(read_index)", BioSequences.LongDNA{4}(String(characters))))
        end
    end
    return reads
end

function _assemble_and_score(
        reads::Vector{FASTX.FASTA.Record},
        references::Vector{String},
        k::Int
)::Tuple{Float64, Int, Bool}
    try
        result = Mycelia.Rhizomorph.assemble_genome(
            reads; k = k,
            graph_mode = Mycelia.Rhizomorph.SingleStrand,
            error_rate = ABLATION_READ_ERROR_RATE)
        contigs = [string(contig) for contig in result.contigs]
        return _reference_recovery(references, contigs), length(contigs), true
    catch
        # An assembly failure IS a recovery failure at this k; record it as such.
        return 0.0, 0, false
    end
end

# Fraction of reference w-mers (w = ABLATION_RECOVERY_W, fixed) present in the
# assembled contigs on either strand. Sensitive to fragmentation and to
# collision-driven tangles (which drop true w-mers around the junction), and
# independent of the assembly k.
function _reference_recovery(references::Vector{String}, contigs::Vector{String})::Float64
    w = ABLATION_RECOVERY_W
    reference_wmers = Set{String}()
    for reference in references
        for index in 1:(length(reference) - w + 1)
            push!(reference_wmers, reference[index:(index + w - 1)])
        end
    end
    isempty(reference_wmers) && return 0.0

    contig_wmers = Set{String}()
    for contig in contigs
        length(contig) < w && continue
        for strand in (contig, _reverse_complement(contig))
            for index in 1:(length(strand) - w + 1)
                push!(contig_wmers, strand[index:(index + w - 1)])
            end
        end
    end
    return length(intersect(reference_wmers, contig_wmers)) / length(reference_wmers)
end

function _reverse_complement(sequence::AbstractString)::String
    return string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(sequence)))
end

# One delta row per (fixture x size-matched (composite,prime) pair): the
# prime-minus-composite recovery gap at nearly-matched k-size.
function _ablation_delta_table(
        per_run::DataFrames.DataFrame,
        fixtures::Vector{AblationFixture}
)::DataFrames.DataFrame
    rows = NamedTuple[]
    for fixture in fixtures
        subset = per_run[per_run.dataset_id .== fixture.dataset_id, :]
        for (composite_k, prime_k) in SIZE_MATCHED_PAIRS
            composite_recovery = only(
                subset[subset.k .== composite_k, :mean_reference_recovery])
            prime_recovery = only(
                subset[subset.k .== prime_k, :mean_reference_recovery])
            push!(rows,
                (
                    dataset_id = fixture.dataset_id,
                    dataset_name = fixture.dataset_name,
                    category = fixture.category,
                    composite_k = composite_k,
                    prime_k = prime_k,
                    composite_recovery = composite_recovery,
                    prime_recovery = prime_recovery,
                    recovery_delta_prime_minus_composite = prime_recovery -
                                                           composite_recovery
                ))
        end
    end
    return DataFrames.DataFrame(rows)
end

# Size-controlled positive-control summary. At each of k=9 (composite) and k=11
# (prime), the COLLISION PENALTY is recovery(distinct) - recovery(shared): the two
# fixtures are structurally identical (same flanks, same 9 bp insert length), so
# this difference isolates the shared-9-mer collision from any k-size effect. A
# real composite-k aliasing signal shows a large penalty at k=9 that collapses at
# k=11 (where the prime k spans the repeat into the differing flanks).
function _positive_control_summary(per_run::DataFrames.DataFrame)::DataFrames.DataFrame
    recovery_at(id,
        k) = only(per_run[
    (per_run.dataset_id .== id) .& (per_run.k .== k), :mean_reference_recovery])
    rows = NamedTuple[]
    for k in (9, 11)
        shared = recovery_at(POSITIVE_CONTROL_SHARED_ID, k)
        distinct = recovery_at(POSITIVE_CONTROL_DISTINCT_ID, k)
        push!(rows,
            (
                k = k,
                k_class = (k in PRIME_K) ? "prime" : "composite",
                shared_recovery = shared,
                distinct_recovery = distinct,
                collision_penalty = distinct - shared
            ))
    end
    return DataFrames.DataFrame(rows)
end

# The positive control FIRES when the collision penalty at composite k=9 is at
# least POSITIVE_CONTROL_MIN_GAP AND is at least POSITIVE_CONTROL_MIN_GAP larger
# than the penalty at prime k=11 — i.e. the shared-9-mer collision demonstrably
# degrades recovery at the composite k and is resolved at the prime k.
function _positive_control_fires(per_run::DataFrames.DataFrame)::Bool
    summary = _positive_control_summary(per_run)
    penalty_9 = only(summary[summary.k .== 9, :collision_penalty])
    penalty_11 = only(summary[summary.k .== 11, :collision_penalty])
    return penalty_9 >= POSITIVE_CONTROL_MIN_GAP &&
           (penalty_9 - penalty_11) >= POSITIVE_CONTROL_MIN_GAP
end

function _write_ablation_figure(
        per_run::DataFrames.DataFrame,
        fixtures::Vector{AblationFixture},
        plots_dir::AbstractString
)::NamedTuple
    mkpath(plots_dir)
    fig = CairoMakie.Figure(size = (1250, 1150), fontsize = 15)
    positions = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)]
    prime_color = :dodgerblue3
    composite_color = :darkorange3

    for (fixture_index, fixture) in enumerate(fixtures)
        row, col = positions[mod1(fixture_index, length(positions))]
        subset = per_run[per_run.dataset_id .== fixture.dataset_id, :]
        title = fixture.category == "natural" ? fixture.dataset_name :
                "$(fixture.dataset_name)  [CONTROL]"
        axis = CairoMakie.Axis(
            fig[row, col],
            title = title,
            xlabel = "assembly k",
            ylabel = "reference recovery (mean over seeds)",
            xticks = ABLATION_K,
            yticks = 0.0:0.25:1.0
        )
        ks = sort(unique(subset.k))
        recovery = [only(subset[subset.k .== k, :mean_reference_recovery]) for k in ks]
        colors = [(k in PRIME_K) ? prime_color : composite_color for k in ks]
        markers = [(k in PRIME_K) ? :circle : :rect for k in ks]
        CairoMakie.lines!(axis, ks, recovery; color = :gray70, linewidth = 1.5)
        CairoMakie.scatter!(
            axis, ks, recovery; color = colors, markersize = 15, marker = markers)
        CairoMakie.ylims!(axis, -0.05, 1.05)
    end

    legend_elements = [
        CairoMakie.MarkerElement(color = prime_color, marker = :circle, markersize = 15),
        CairoMakie.MarkerElement(color = composite_color, marker = :rect, markersize = 15)
    ]
    CairoMakie.Legend(
        fig[4, 1:2], legend_elements,
        ["prime k ($(join(PRIME_K, ", ")))", "composite k ($(join(COMPOSITE_K, ", ")))"];
        orientation = :horizontal, framevisible = false)
    CairoMakie.Label(
        fig[0, 1:2],
        "Rhizomorph de-novo reference recovery vs k on repeat-rich fixtures";
        fontsize = 18, font = :bold)

    png_path = joinpath(plots_dir, "prime_k_ablation_recovery_vs_k.png")
    svg_path = joinpath(plots_dir, "prime_k_ablation_recovery_vs_k.svg")
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
