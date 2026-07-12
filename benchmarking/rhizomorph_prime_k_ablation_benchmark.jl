# Prime (coprime) vs composite (factor-aligning) k ablation for the Rhizomorph
# assembler, in the DE-NOVO regime, DESIGNED to isolate FACTOR-ALIGNMENT from
# k-SIZE across MANY repeat periods (bead td-tjym).
#
# THE CLAIM UNDER TEST (the real mathematical property of primes). A prime k is
# coprime to EVERY period p < k, so it can never factor-align with a tandem repeat
# of any period below it. A composite k = p*q ALIASES repeats of period p and of
# period q. The purported value of a prime k is therefore NOT "k > repeat length"
# (mere spanning) — it is "k shares no factor with the repeat period", across the
# whole spectrum of periods. To test this you must isolate factor-alignment from
# k-size, which requires MANY periods and, for each period, a factor-sharing
# composite k that is comfortably larger than the repeat unit compared against a
# SIZE-MATCHED prime k coprime to that period.
#
# HISTORY on this branch. v1 (reviewed on #404) was structurally forced to
# delta=0 (clean-reference graph + both-endpoint-anchored decode). v2 moved to the
# valid de-novo regime but tested a single period (p=3) with a single 9 bp shared
# repeat, which conflated factor-alignment with spanning (k=11 beats k=9 only
# because 11 > 9 = the repeat length). v3 (this file) runs the systematic
# (period p x k) sweep that actually isolates factor-alignment.
#
# EXPERIMENT.
#   - Fixtures: tandem repeats of period p in {2,3,4,5,6,7}, a p-bp unit repeated
#     to a fixed-length region, embedded in unique flanks so the de-novo graph
#     must traverse the repeat; plus a non-repetitive control (the k-SIZE
#     reference).
#   - k grid: primes {11,13,17,19,23,29,31} and composites {9,15,21,25,27,33,35},
#     chosen so that for each odd period p there is a composite k with p|k AND a
#     size-matched prime k coprime to p.
#   - Graph is built FROM ERROR-CONTAINING READS at each k (no anchoring); recovery
#     = fraction of reference w-mers (fixed w, independent of assembly k) present
#     in the assembled contigs, mean over deterministic seeds.
#   - THE ISOLATING COMPARISON: for each odd period p, above the k-size threshold,
#     recovery(size-matched coprime PRIME k) - recovery(factor-sharing COMPOSITE k
#     with p|k). If factor-alignment mattered, this delta would be positive
#     (composite aliases, prime does not) DESPITE both k spanning the unit.
#
# OBSERVED RESULT (numbers are the deliverable; see the committed (period x k)
# table). Factor-alignment produces NO recovery deficit isolated from k-size:
#   - k-SIZE dominates: below k ~= 13-15 recovery is low for EVERY fixture,
#     including the non-repetitive control (this is an error-tolerance /
#     graph-fragmentation effect, not repeat-specific). Above it recovery
#     saturates near 1.0 for every fixture.
#   - Above the threshold, a factor-sharing composite k (p|k) recovers the SAME as
#     its size-matched coprime prime k — the isolating deltas are ~0 for every
#     odd period. The (period x k) heatmap shows NO darkening on the p|k cells.
#   - This matches the math: a period-p tandem yields exactly p distinct k-mers (a
#     p-cycle) at ANY k, so single-tandem de-novo recovery is period/primality
#     independent; de-novo graph cleaning washes out any factor-alignment.
#
# HARNESS SENSITIVITY (so the null is not mere insensitivity). A separate
# engineered COLLISION control — two genes sharing a 9 bp period-3 repeat, vs a
# composition-matched no-collision variant — DOES produce a large, isolated
# recovery deficit at k=9 that collapses once k spans the shared repeat. This
# proves the recovery metric CAN detect an aliasing-driven deficit when one truly
# exists; it is a k-vs-repeat-length (collision/spanning) phenomenon, NOT a
# primality mechanism.
#
# INTERPRETATION for the manuscript. In this de-novo w-mer-recovery regime the
# prime k-ladder's value is NOT per-k robustness to factor-alignment (which washes
# out here). Its defensible basis is:
#   (1) ODD k for DNA/RNA: an even-length k-mer can equal its own reverse
#       complement (a palindrome), collapsing RC pairs and creating self-loops;
#       odd k makes palindromic k-mers impossible. (Note in this sweep: among ODD
#       k, only ODD periods 3,5,7 can factor-align at all — even periods 2,4,6 are
#       automatically coprime to every odd k, a built-in control.) Even k are
#       acceptable for amino-acid / natural-language alphabets, which have no
#       reverse complement; this benchmark keeps DNA fixtures on ODD k.
#   (2) Among odd k, primes additionally avoid factor-alignment with odd periods —
#       a real number-theoretic property, but one that does NOT manifest as a
#       de-novo recovery difference at this scale.
#   (3) Geometric ladder spacing (a separate efficiency concern, not tested here).
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

# Tandem-repeat periods. Odd periods (3,5,7) admit a factor-sharing composite k in
# the grid below; even periods (2,4,6) are coprime to EVERY odd k and so are
# built-in "cannot factor-align" controls.
const ABLATION_PERIODS = (2, 3, 4, 5, 6, 7)
# All k are ODD (see header: even k-mers can equal their own reverse complement in
# DNA/RNA). Composites are chosen so each odd period has a factor-sharing k:
# 3 | {9,15,21,27,33}, 5 | {15,25,35}, 7 | {21,35}.
const ABLATION_PRIME_K = (11, 13, 17, 19, 23, 29, 31)
const ABLATION_COMPOSITE_K = (9, 15, 21, 25, 27, 33, 35)
const ABLATION_K = sort(collect(Iterators.flatten((
    ABLATION_PRIME_K, ABLATION_COMPOSITE_K))))

# Tandem region length (bp) and unique flank length (bp). The region is longer
# than the largest k, so recovery is never a whole-region spanning artifact.
const ABLATION_REGION_LENGTH = 45
const ABLATION_FLANK_LENGTH = 30

# Read-simulation parameters (SMOKE, deterministic).
const ABLATION_COVERAGE = 30
const ABLATION_READ_LENGTH = 45
const ABLATION_READ_ERROR_RATE = 0.02
const ABLATION_SEEDS = (1, 2, 3)
const ABLATION_ALPHABET = collect("ACGT")
# Fixed reference-recovery window, INDEPENDENT of the assembly k.
const ABLATION_RECOVERY_W = 15
# The non-repetitive control defines the k-SIZE (error-tolerance) threshold: the
# smallest k at which it recovers at least this fraction. Factor-alignment is
# evaluated only AT OR ABOVE that threshold (so any deficit is not a size effect).
const ABLATION_SIZE_THRESHOLD_RECOVERY = 0.90
# A factor-sharing composite k counts as "isolating a factor-alignment deficit"
# only if a size-matched coprime prime k out-recovers it by at least this margin,
# above the size threshold. Observed deltas are ~0, so this does NOT fire.
const FACTOR_ALIGNMENT_MIN_DELTA = 0.10
# Harness-sensitivity collision control bar (k=9 collision penalty; see below).
const COLLISION_MIN_GAP = 0.10
const COLLISION_SHARED_ID = "collision_shared_9bp"
const COLLISION_DISTINCT_ID = "collision_distinct_9bp"
const CONTROL_NONREPETITIVE_ID = "control_nonrepetitive"

const ABLATION_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "rhizomorph_prime_k_ablation"
)

struct AblationFixture
    dataset_id::String
    dataset_name::String
    category::String            # "tandem", "control", "collision_shared", "collision_distinct"
    period::Int                 # tandem period, or 0 when not a single-period tandem
    references::Vector{String}
    note::String
end

function main(args::Vector{String} = ARGS)::Nothing
    output_dir = _ablation_arg_value(args, "--output-dir", ABLATION_DEFAULT_OUTPUT_DIR)
    write_plots = !("--skip-plots" in args)
    artifacts = run_prime_k_ablation_benchmark(output_dir; write_plots = write_plots)
    println("Wrote Rhizomorph prime(coprime)-vs-composite(factor-aligning) k ablation:")
    println("  root: $(artifacts.root)")
    println("  per_run_csv: $(artifacts.per_run_csv)")
    println("  factor_alignment_csv: $(artifacts.factor_alignment_csv)")
    println("  collision_control_csv: $(artifacts.collision_control_csv)")
    println("  figure_png: $(artifacts.figure_png)")
    println()
    println("k-SIZE threshold (smallest k with non-repetitive control recovery >= " *
            "$(ABLATION_SIZE_THRESHOLD_RECOVERY)): k = $(artifacts.size_threshold_k)")
    println()
    println("ISOLATING COMPARISON — factor-sharing composite (p|k) vs size-matched " *
            "coprime prime, ABOVE the k-size threshold:")
    if DataFrames.nrow(artifacts.factor_alignment_table) == 0
        println("  (no above-threshold factor-sharing pairs found)")
    else
        for row in eachrow(artifacts.factor_alignment_table)
            println("  period p=$(row.period)  composite k=$(row.composite_k) (=p*$(row.composite_k ÷ row.period)) " *
                    "recovery=$(round(row.composite_recovery; digits = 3))  vs  " *
                    "coprime prime k=$(row.prime_k) recovery=$(round(row.prime_recovery; digits = 3))  " *
                    "delta(prime-composite)=$(round(row.recovery_delta_prime_minus_composite; digits = 3))")
        end
    end
    println()
    println("Factor-alignment isolated above threshold (any prime-minus-composite " *
            ">= $(FACTOR_ALIGNMENT_MIN_DELTA)): $(artifacts.factor_alignment_isolated)")
    println()
    println("HARNESS-SENSITIVITY collision control (engineered 9 bp shared repeat, " *
            "composition-matched; NOT a primality test):")
    for row in eachrow(artifacts.collision_control_table)
        println("  k=$(row.k)  shared=$(round(row.shared_recovery; digits = 3))  " *
                "distinct=$(round(row.distinct_recovery; digits = 3))  " *
                "collision_penalty=$(round(row.collision_penalty; digits = 3))")
    end
    println("Collision control fires (k=9 penalty >= $(COLLISION_MIN_GAP) and collapses " *
            "once k spans the repeat): $(artifacts.collision_control_fires)")
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

    size_threshold_k = _size_threshold_k(per_run)
    factor_alignment = _factor_alignment_table(per_run, size_threshold_k)
    collision_control = _collision_control_table(per_run)
    factor_alignment_isolated = _factor_alignment_isolated(factor_alignment)
    collision_fires = _collision_control_fires(collision_control)

    artifacts = write_benchmark_artifacts(
        [
            "prime_k_ablation_per_run" => per_run,
            "prime_k_ablation_factor_alignment" => factor_alignment,
            "prime_k_ablation_collision_control" => collision_control
        ];
        output_dir = output_dir,
        run_id = "prime_k_ablation_denovo_multiperiod_20260712",
        scale = "local-smoke",
        dataset_ids = [fixture.dataset_id for fixture in fixtures],
        command_args = [
            "julia", "--project=.",
            "benchmarking/rhizomorph_prime_k_ablation_benchmark.jl"],
        metadata = Dict(
            "bead" => "td-tjym",
            "benchmark" => "prime_coprime_vs_composite_factor_alignment_denovo_multiperiod",
            "claim_under_test" => "a prime k is coprime to every period p<k so it cannot factor-align with a repeat of any period; a composite k=p*q aliases period-p and period-q repeats. Value of prime k is coprimality, not spanning (k>repeat length).",
            "regime" => "de-novo assembly from error-containing reads (graph built from reads; no endpoint anchoring)",
            "periods" => collect(ABLATION_PERIODS),
            "prime_k" => collect(ABLATION_PRIME_K),
            "composite_k" => collect(ABLATION_COMPOSITE_K),
            "isolating_comparison" => "for each odd period p, ABOVE the k-size threshold, recovery(size-matched coprime prime k) - recovery(factor-sharing composite k with p|k)",
            "coverage" => ABLATION_COVERAGE,
            "read_length" => ABLATION_READ_LENGTH,
            "read_error_rate" => ABLATION_READ_ERROR_RATE,
            "seeds" => collect(ABLATION_SEEDS),
            "recovery_metric" => "fraction of reference $(ABLATION_RECOVERY_W)-mers present in assembled contigs (either strand), mean over seeds",
            "observed_result" => "factor-alignment produces NO recovery deficit isolated from k-size: above the threshold, factor-sharing composite k recovers the same as size-matched coprime prime k (deltas ~0); a period-p tandem is a p-cycle at any k, so single-tandem recovery is primality-independent and de-novo cleaning washes out aliasing.",
            "harness_sensitivity" => "a separate engineered 9 bp shared-repeat collision control (composition-matched) DOES produce a large isolated deficit at k=9 that collapses once k spans the repeat — proving the metric detects aliasing when it exists; this is a k-vs-repeat-length effect, NOT primality.",
            "reverse_complement_odd_k_note" => "the ladder samples ODD k because an even-length DNA/RNA k-mer can equal its own reverse complement (a palindrome), collapsing RC pairs and creating self-loops; odd k forbids palindromic k-mers. Among odd k, primes additionally avoid factor-alignment with odd periods. Even k are acceptable for amino-acid / natural-language alphabets (no reverse complement); DNA fixtures here stay on odd k.",
            "size_threshold_recovery" => ABLATION_SIZE_THRESHOLD_RECOVERY,
            "factor_alignment_min_delta" => FACTOR_ALIGNMENT_MIN_DELTA
        ),
        table_context_columns = Dict(
            "prime_k_ablation_per_run" => Dict("benchmark_dataset_id" => "dataset_id")
        )
    )

    figure_png = ""
    figure_svg = ""
    if write_plots
        paths = _write_heatmap_figure(
            per_run, fixtures, collision_control, size_threshold_k, artifacts.layout.plots)
        figure_png = paths.png
        figure_svg = paths.svg
    end

    return (
        root = artifacts.root,
        per_run_csv = artifacts.tables["prime_k_ablation_per_run"].table,
        factor_alignment_csv = artifacts.tables["prime_k_ablation_factor_alignment"].table,
        collision_control_csv = artifacts.tables["prime_k_ablation_collision_control"].table,
        index = artifacts.index,
        provenance = artifacts.provenance,
        figure_png = figure_png,
        figure_svg = figure_svg,
        per_run_table = per_run,
        factor_alignment_table = factor_alignment,
        collision_control_table = collision_control,
        size_threshold_k = size_threshold_k,
        factor_alignment_isolated = factor_alignment_isolated,
        collision_control_fires = collision_fires,
        per_run_rows = DataFrames.nrow(per_run)
    )
end

function ablation_fixtures()::Vector{AblationFixture}
    flank_left = _nonrepetitive_sequence(ABLATION_FLANK_LENGTH, 1)
    flank_right = _nonrepetitive_sequence(ABLATION_FLANK_LENGTH, 2)
    fixtures = AblationFixture[]
    for period in ABLATION_PERIODS
        region = _tandem_region(period, ABLATION_REGION_LENGTH)
        push!(fixtures,
            AblationFixture(
                "tandem_period$(period)",
                "Period-$(period) tandem in unique flanks",
                "tandem",
                period,
                [flank_left * region * flank_right],
                isodd(period) ?
                "odd period — factor-aligns with composite k divisible by $(period)" :
                "even period — coprime to every odd k (built-in control)"
            ))
    end
    push!(fixtures,
        AblationFixture(
            CONTROL_NONREPETITIVE_ID,
            "Non-repetitive control (defines the k-size threshold)",
            "control",
            0,
            [_nonrepetitive_sequence(ABLATION_REGION_LENGTH + 2 * ABLATION_FLANK_LENGTH, 3)],
            "no dominant period; its low-k recovery collapse is the pure k-size effect"
        ))
    # HARNESS-SENSITIVITY collision controls (NOT a primality test). Two genes
    # share a 9 bp period-3 repeat -> composite k=9 collides on the identical
    # 9-mer; once k spans the 9 bp repeat (k>=10) the collision resolves. The
    # composition-matched distinct variant (gene2 uses TGAx3, a permutation of
    # ATGx3: same GC, no shared 9-mer) isolates the collision from k-size and GC.
    push!(fixtures,
        AblationFixture(
            COLLISION_SHARED_ID,
            "Collision control (shared): two genes sharing a 9 bp period-3 repeat",
            "collision_shared",
            3,
            [
                _nonrepetitive_sequence(20, 4) * repeat("ATG", 3) *
                _nonrepetitive_sequence(20, 5),
                _nonrepetitive_sequence(20, 6) * repeat("ATG", 3) *
                _nonrepetitive_sequence(20, 7)
            ],
            "engineered 9 bp SHARED repeat; harness-sensitivity control, not primality"
        ))
    push!(fixtures,
        AblationFixture(
            COLLISION_DISTINCT_ID,
            "Collision control (distinct): composition-matched, non-shared 9 bp repeat",
            "collision_distinct",
            3,
            [
                _nonrepetitive_sequence(20, 4) * repeat("ATG", 3) *
                _nonrepetitive_sequence(20, 5),
                _nonrepetitive_sequence(20, 6) * repeat("TGA", 3) *
                _nonrepetitive_sequence(20, 7)
            ],
            "gene2 uses TGAx3 (permutation of ATG: same composition, no shared 9-mer)"
        ))
    return fixtures
end

# Tandem region: a p-bp unit repeated to exactly `region_length` bp. The unit is
# deterministic and period-specific so different periods use different sequences.
function _tandem_region(period::Int, region_length::Int)::String
    unit = _nonrepetitive_sequence(period, 200 + period)
    n_copies = cld(region_length, period)
    return first(repeat(unit, n_copies), region_length)
end

# Deterministic pseudo-random (non-repetitive) DNA via a linear-congruential walk;
# no RNG dependency, byte-stable. Distinct seed_index -> distinct sequence.
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
    errors = String[]
    for seed in ABLATION_SEEDS
        reads = _simulate_reads(fixture.references, seed)
        outcome = _assemble_and_score(reads, fixture.references, k)
        push!(recoveries, outcome.recovery)
        push!(contig_counts, outcome.n_contigs)
        push!(assembled_flags, outcome.assembled)
        isempty(outcome.error) || push!(errors, "seed $(seed): $(outcome.error)")
    end
    return (
        dataset_id = fixture.dataset_id,
        dataset_name = fixture.dataset_name,
        category = fixture.category,
        period = fixture.period,
        note = fixture.note,
        k = k,
        k_class = (k in ABLATION_PRIME_K) ? "prime" : "composite",
        factor_aligned = fixture.period > 1 && k % fixture.period == 0,
        n_references = length(fixture.references),
        reference_length = sum(length, fixture.references),
        coverage = ABLATION_COVERAGE,
        read_error_rate = ABLATION_READ_ERROR_RATE,
        n_seeds = length(ABLATION_SEEDS),
        recovery_w = ABLATION_RECOVERY_W,
        mean_reference_recovery = Statistics.mean(recoveries),
        min_reference_recovery = minimum(recoveries),
        max_reference_recovery = maximum(recoveries),
        mean_contig_count = Statistics.mean(contig_counts),
        all_seeds_assembled = all(assembled_flags),
        assembly_errors = isempty(errors) ? "" : join(errors, " | ")
    )
end

# Simulate a read set tiling every reference at ABLATION_COVERAGE with per-base
# substitution errors. Deterministic given (references, seed).
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
)::NamedTuple
    try
        result = Mycelia.Rhizomorph.assemble_genome(
            reads; k = k,
            graph_mode = Mycelia.Rhizomorph.SingleStrand,
            error_rate = ABLATION_READ_ERROR_RATE)
        contigs = [string(contig) for contig in result.contigs]
        return (recovery = _reference_recovery(references, contigs),
            n_contigs = length(contigs), assembled = true, error = "")
    catch exception
        # Never swallow interrupts. A genuine assembly failure IS a recovery
        # failure at this k, but must be VISIBLE: capture the message + warn.
        exception isa InterruptException && rethrow()
        message = sprint(showerror, exception)
        @warn "assemble_genome failed" k=k error=first(message, 200)
        return (recovery = 0.0, n_contigs = 0, assembled = false, error = message)
    end
end

# Fraction of reference w-mers (fixed w) present in the assembled contigs on either
# strand. Independent of the assembly k; sensitive to fragmentation and to
# collision-driven tangles (which drop true w-mers around the junction).
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

function _tandem_recovery(per_run::DataFrames.DataFrame, period::Int, k::Int)::Float64
    return only(per_run[
        (per_run.dataset_id .== "tandem_period$(period)") .& (per_run.k .== k),
        :mean_reference_recovery])
end

# The k-SIZE (error-tolerance) threshold: the smallest k at which the
# non-repetitive control recovers at least ABLATION_SIZE_THRESHOLD_RECOVERY. Below
# this k, recovery is limited by graph fragmentation for EVERY fixture (a k-size
# effect, not repeat-specific); factor-alignment is evaluated only at/above it.
function _size_threshold_k(per_run::DataFrames.DataFrame)::Int
    control = per_run[per_run.dataset_id .== CONTROL_NONREPETITIVE_ID, :]
    passing = sort(control[control.mean_reference_recovery .>= ABLATION_SIZE_THRESHOLD_RECOVERY, :], :k)
    return DataFrames.nrow(passing) == 0 ? maximum(ABLATION_K) : first(passing.k)
end

# THE ISOLATING COMPARISON. For each odd period p and each factor-sharing
# composite k (p|k) at/above the size threshold, pair it with the size-matched
# coprime prime k (nearest prime in the grid that is at/above threshold; every
# grid prime > p so all are coprime to p). Report prime-minus-composite recovery:
# a positive delta would mean the composite aliases where the prime does not.
function _factor_alignment_table(per_run::DataFrames.DataFrame, threshold_k::Int)::DataFrames.DataFrame
    rows = NamedTuple[]
    eligible_primes = [pk for pk in ABLATION_PRIME_K if pk >= threshold_k]
    for period in ABLATION_PERIODS
        isodd(period) || continue  # even periods cannot factor-align with odd k
        for composite_k in ABLATION_COMPOSITE_K
            (composite_k % period == 0 && composite_k >= threshold_k) || continue
            isempty(eligible_primes) && continue
            prime_k = eligible_primes[argmin(abs.(eligible_primes .-
                                                                  composite_k))]
            composite_recovery = _tandem_recovery(per_run, period, composite_k)
            prime_recovery = _tandem_recovery(per_run, period, prime_k)
            push!(rows,
                (
                    period = period,
                    composite_k = composite_k,
                    prime_k = prime_k,
                    size_threshold_k = threshold_k,
                    composite_recovery = composite_recovery,
                    prime_recovery = prime_recovery,
                    recovery_delta_prime_minus_composite = prime_recovery -
                                                           composite_recovery
                ))
        end
    end
    return DataFrames.DataFrame(rows)
end

# Factor-alignment is "isolated" only if SOME above-threshold factor-sharing
# composite is out-recovered by its size-matched coprime prime by at least
# FACTOR_ALIGNMENT_MIN_DELTA. Observed: no — the deltas are ~0.
function _factor_alignment_isolated(factor_alignment::DataFrames.DataFrame)::Bool
    DataFrames.nrow(factor_alignment) == 0 && return false
    return maximum(factor_alignment.recovery_delta_prime_minus_composite) >=
           FACTOR_ALIGNMENT_MIN_DELTA
end

# Harness-sensitivity collision control. At each of k=9 and k=11 the COLLISION
# PENALTY is recovery(distinct) - recovery(shared): the two fixtures are
# structurally identical and base-composition-matched, so this isolates the
# shared-9-mer collision from k-size and GC. The penalty is large at k=9
# (collision) and collapses at k=11 because an 11-mer SPANS the 9 bp repeat — a
# k-vs-repeat-length effect, NOT primality.
function _collision_control_table(per_run::DataFrames.DataFrame)::DataFrames.DataFrame
    recovery_at(id,
        k) = only(per_run[
    (per_run.dataset_id .== id) .& (per_run.k .== k), :mean_reference_recovery])
    rows = NamedTuple[]
    for k in (9, 11)
        shared = recovery_at(COLLISION_SHARED_ID, k)
        distinct = recovery_at(COLLISION_DISTINCT_ID, k)
        push!(rows,
            (
                k = k,
                k_class = (k in ABLATION_PRIME_K) ? "prime" : "composite",
                shared_recovery = shared,
                distinct_recovery = distinct,
                collision_penalty = distinct - shared
            ))
    end
    return DataFrames.DataFrame(rows)
end

function _collision_control_fires(collision_control::DataFrames.DataFrame)::Bool
    penalty_9 = only(collision_control[collision_control.k .== 9, :collision_penalty])
    penalty_11 = only(collision_control[collision_control.k .== 11, :collision_penalty])
    return penalty_9 >= COLLISION_MIN_GAP && (penalty_9 - penalty_11) >= COLLISION_MIN_GAP
end

# THE MONEY FIGURE: a (repeat period x k) recovery heatmap. If factor-alignment
# mattered, the p|k cells (ringed) would be darker than coprime cells at the same
# k; observed, they are not. A bottom row shows the non-repetitive control (the
# k-size reference), and a side panel shows the collision-control sensitivity.
function _write_heatmap_figure(
        per_run::DataFrames.DataFrame,
        fixtures::Vector{AblationFixture},
        collision_control::DataFrames.DataFrame,
        threshold_k::Int,
        plots_dir::AbstractString
)::NamedTuple
    mkpath(plots_dir)
    heatmap_ids = vcat(
        ["tandem_period$(p)" for p in ABLATION_PERIODS], [CONTROL_NONREPETITIVE_ID])
    row_labels = vcat(["p=$(p)" for p in ABLATION_PERIODS], ["non-rep"])
    ks = ABLATION_K
    n_rows = length(heatmap_ids)
    n_cols = length(ks)
    matrix = Array{Float64}(undef, n_rows, n_cols)
    for (r, id) in enumerate(heatmap_ids)
        for (c, k) in enumerate(ks)
            matrix[r, c] = only(per_run[
            (per_run.dataset_id .== id) .& (per_run.k .== k), :mean_reference_recovery])
        end
    end

    fig = CairoMakie.Figure(size = (1350, 780), fontsize = 15)
    axis = CairoMakie.Axis(
        fig[1, 1],
        title = "De-novo reference recovery vs (repeat period, k) — rings mark p | k (factor-aligned)",
        xlabel = "assembly k", ylabel = "repeat period",
        xticks = (1:n_cols, string.(ks)),
        yticks = (1:n_rows, row_labels)
    )
    heat = CairoMakie.heatmap!(
        axis, 1:n_cols, 1:n_rows, permutedims(matrix);
        colormap = :viridis, colorrange = (0.0, 1.0))
    # Ring the factor-aligned (p|k) cells.
    ring_x = Int[]
    ring_y = Int[]
    for (r, id) in enumerate(heatmap_ids)
        startswith(id, "tandem_period") || continue
        period = parse(Int, replace(id, "tandem_period" => ""))
        for (c, k) in enumerate(ks)
            if k % period == 0
                push!(ring_x, c)
                push!(ring_y, r)
            end
        end
    end
    CairoMakie.scatter!(
        axis, ring_x, ring_y; marker = :circle, markersize = 22,
        color = (:white, 0.0), strokecolor = :red, strokewidth = 2.5)
    # k-size threshold marker.
    threshold_col = findfirst(==(threshold_k), ks)
    if threshold_col !== nothing
        CairoMakie.vlines!(
            axis, [threshold_col - 0.5]; color = :white, linestyle = :dash, linewidth = 2)
    end
    CairoMakie.Colorbar(fig[1, 2], heat, label = "reference recovery (mean over seeds)")

    # Collision-control sensitivity panel.
    axis2 = CairoMakie.Axis(
        fig[2, 1],
        title = "Harness-sensitivity control: 9 bp shared-repeat collision (NOT a primality test)",
        xlabel = "assembly k", ylabel = "reference recovery",
        xticks = ([9, 11], ["9", "11"]), yticks = 0.0:0.25:1.0)
    shared = [only(collision_control[collision_control.k .== k, :shared_recovery])
              for k in (9, 11)]
    distinct = [only(collision_control[collision_control.k .== k, :distinct_recovery])
                for k in (9, 11)]
    CairoMakie.scatterlines!(axis2, [9, 11], shared; color = :darkorange3, linewidth = 3,
        markersize = 14, label = "shared 9-mer (collides at k=9)")
    CairoMakie.scatterlines!(axis2, [9, 11], distinct; color = :dodgerblue3, linewidth = 3,
        markersize = 14, marker = :rect, label = "distinct (composition-matched)")
    CairoMakie.ylims!(axis2, -0.05, 1.05)
    CairoMakie.axislegend(axis2; position = :rb)
    CairoMakie.rowsize!(fig.layout, 2, CairoMakie.Relative(0.28))

    CairoMakie.Label(
        fig[0, 1:2],
        "Prime(coprime) vs composite(factor-aligning) k — factor-alignment is NOT isolated from k-size";
        fontsize = 18, font = :bold)

    png_path = joinpath(plots_dir, "prime_k_ablation_period_by_k_heatmap.png")
    svg_path = joinpath(plots_dir, "prime_k_ablation_period_by_k_heatmap.svg")
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
