"""
Iterative Maximum Likelihood Assembly - Phase 5.2a Implementation

This module implements the iterative maximum likelihood assembly system with:
- Complete FASTQ I/O processing per iteration
- Statistical path resampling with likelihood calculations
- Viterbi algorithm integration for optimal path finding
- Timestamped output files for tracking read evolution
- Memory-efficient read set processing

Part of the Mycelia bioinformatics package's iterative assembler framework.
"""

# =============================================================================
# Sparsity Detection and Memory Estimation
# =============================================================================

"""
Calculate k-mer sparsity for a given k-mer size.
Returns fraction of possible k-mers that are NOT observed.
"""
function calculate_sparsity(reads::Vector{<:FASTX.FASTQ.Record}, k::Int)::Float64
    observed_kmers = Set{String}()

    for record in reads
        seq = FASTX.sequence(String, record)
        if length(seq) >= k
            for i in 1:(length(seq) - k + 1)
                kmer = seq[i:(i + k - 1)]
                push!(observed_kmers, kmer)
            end
        end
    end

    unique_observed = length(observed_kmers)
    max_possible = 4.0^k
    return 1.0 - (unique_observed / max_possible)
end

"""
Analyze k-mer coverage distribution to detect if errors are singletons.
Returns true if low-coverage k-mers (likely errors) are well-separated from high-coverage ones.
"""
function errors_are_singletons(reads::Vector{<:FASTX.FASTQ.Record}, k::Int; singleton_threshold::Int = 2)::Bool
    kmer_counts = Dict{String, Int}()

    for record in reads
        seq = FASTX.sequence(String, record)
        if length(seq) >= k
            for i in 1:(length(seq) - k + 1)
                kmer = seq[i:(i + k - 1)]
                kmer_counts[kmer] = get(kmer_counts, kmer, 0) + 1
            end
        end
    end

    coverage_values = collect(values(kmer_counts))
    total_unique = length(coverage_values)
    if total_unique == 0
        return false
    end
    singleton_count = count(x -> x <= singleton_threshold, coverage_values)
    singleton_fraction = singleton_count / total_unique

    if singleton_count > 0 && total_unique > singleton_count
        non_singleton_min = minimum(filter(x -> x > singleton_threshold, coverage_values))
        return singleton_fraction > 0.1 && non_singleton_min > singleton_threshold * 2
    end

    return false
end

"""
Find the optimal starting k-mer size using sparsity detection.
Only considers prime k-mer sizes for optimal performance.
"""
function find_initial_k(reads::Vector{<:FASTX.FASTQ.Record};
        k_range::Vector{Int} = Primes.primes(3, 51),
        sparsity_threshold::Float64 = 0.5)::Int
    for k in k_range
        sparsity = calculate_sparsity(reads, k)
        if sparsity > sparsity_threshold && errors_are_singletons(reads, k)
            return k
        end
    end

    return isempty(k_range) ? 3 : first(k_range)
end

"""
Find the next prime number greater than current_k.
For k-mer progression, we prefer odd numbers and especially primes.
"""
function next_prime_k(current_k::Int; max_k::Int = 1000)::Int
    next_prime = Primes.nextprime(current_k + 1)
    return next_prime <= max_k ? next_prime : current_k
end

"""
    build_k_ladder(initial_k, max_k; k_ladder=nothing, n_k_rungs=nothing)

Compute the ordered set of k-mer sizes the iterative corrector should walk.

- If `k_ladder` is given, use those values that fall in `[initial_k, max_k]`
  (deduplicated and sorted). NOTE: the assembly loop always STARTS at
  `find_initial_k` regardless of the ladder, so if `min(k_ladder) > initial_k` the
  auto-detected initial k is processed first and then the ladder is followed — the
  ladder specifies the k values to visit AFTER the initial k, not a replacement
  for it. (In `n_k_rungs` mode the first rung is pinned to `initial_k`, so they
  coincide.)
- Else if `n_k_rungs` is given, build an ~`n_k_rungs`-rung geometric ladder from
  `initial_k` to `max_k` (LoRMA-style coarse progression). INTERMEDIATE rungs are
  snapped to PRIME values (not merely odd): a composite k like 9=3x3 or 15=3x5
  makes period-3/period-5 tandem repeats collapse to self-overlapping k-mers, the
  worst case for de Bruijn correction — and period-p aliasing is a LOW-k
  phenomenon, so it is exactly the mid-range auto-selected rungs that must avoid
  it. The first rung is pinned to `initial_k` (which `find_initial_k` already
  draws from `Primes.primes`) and the last to the largest ODD `<= max_k`.

  The TOP rung is deliberately `max_k` (largest odd ≤ max_k), NOT `prevprime(max_k)`:
  the final-pass graph-reuse optimization (td-04tb) requires the corrector's final
  k to equal the re-assembly k (`config.k == max_k`), so topping below `max_k`
  would silently disable reuse AND reduce the user's requested assembly resolution.
  At the high k a top rung sits at, k-mers are long enough that a composite k
  rarely spans a short-period repeat, so the aliasing cost there is negligible —
  the win is in the mid-range intermediates, which ARE prime. (The single-k B8
  accuracy benchmark, which has no re-assembly, does use a prime k.)
- Else return `nothing`, signalling the caller to use the original prime-by-prime
  progression (`next_prime_k`) — this preserves legacy behavior.

Returns a sorted `Vector{Int}` (or `nothing`).
"""
function build_k_ladder(initial_k::Int, max_k::Int;
        k_ladder::Union{Nothing, Vector{Int}} = nothing,
        n_k_rungs::Union{Nothing, Int} = nothing)::Union{Nothing, Vector{Int}}
    if k_ladder !== nothing
        rungs = sort(unique(filter(k -> initial_k <= k <= max_k, k_ladder)))
        return isempty(rungs) ? [initial_k] : rungs
    end
    if n_k_rungs === nothing
        return nothing
    end
    if initial_k >= max_k
        return [initial_k]
    end
    n = max(2, n_k_rungs)
    # Top rung stays at max_k (largest odd ≤ max_k) to preserve final-pass graph
    # reuse (td-04tb): the reuse gate needs final_graph_k == reassembly_k == max_k.
    top = isodd(max_k) ? max_k : max_k - 1
    if top <= initial_k
        return [initial_k]
    end
    ratio = (top / initial_k)^(1 / (n - 1))
    ks = Int[]
    for i in 1:n
        kv = round(Int, initial_k * ratio^(i - 1))
        kv = Primes.nextprime(kv)        # snap intermediate rungs UP to a prime
        kv = clamp(kv, initial_k, top)
        push!(ks, kv)
    end
    # First rung pinned to initial_k (prime in practice); last rung to the top
    # (= max_k, kept composite-tolerant so re-assembly graph reuse still fires).
    ks[1] = initial_k
    ks[end] = top
    return sort(unique(ks))
end

"""
    _next_k_in_progression(current_k, max_k, k_schedule)

Return the next k-mer size after `current_k`. When `k_schedule === nothing`, use
the original prime progression (`next_prime_k`). Otherwise advance to the next
scheduled rung strictly greater than `current_k`, or return `current_k` (a
fixed point, which the main loop treats as "stop") when none remain.
"""
function _next_k_in_progression(current_k::Int, max_k::Int,
        k_schedule::Union{Nothing, Vector{Int}})::Int
    if k_schedule === nothing
        return next_prime_k(current_k; max_k = max_k)
    end
    larger = filter(>(current_k), k_schedule)
    return isempty(larger) ? current_k : minimum(larger)
end

# Default ADAPTIVE low-k decode-gate density threshold (td-9h5r): if, at a given
# k-rung, the post-Stage-0 hard-window gate would still send AT LEAST this fraction
# of reads to the expensive per-read graph-Viterbi decode, the gate is judged
# NON-DISCRIMINATING (it is not reducing decode volume) and the whole decode is
# skipped for that pass — the reads fall back on Stage 0 cheap correction +
# skip-solid, and the whole-read Viterbi is deferred to higher, selective rungs
# where the gate actually filters. 0.90 leaves a genuinely selective low-k decode
# (e.g. a gate that skips ~40-60% of reads) ON, and only gates the dense-graph
# "decode almost everything" pathology the #370 profile flagged.
const _DEFAULT_DECODE_GATE_DENSITY = 0.90

"""
Estimate memory usage for a graph with a given number of k-mers.
Provides a rough estimate for memory monitoring.
"""
function estimate_memory_usage(num_kmers::Int, k::Int)::Int
    kmer_size = k * 1
    vertex_overhead = 100
    edge_overhead = 50

    estimated_bytes = num_kmers * (kmer_size + vertex_overhead + edge_overhead * 2)
    return round(Int, estimated_bytes * 1.2)
end

function _label_k_length(label)::Int
    if label isa Qualmer
        return length(label)
    elseif label isa Kmers.Kmer
        return length(label)
    elseif label isa BioSequences.BioSequence
        return length(label)
    elseif label isa AbstractString
        return length(label)
    end
    return length(string(label))
end

"""
Check if memory usage is within acceptable limits.
"""
function check_memory_limits(graph, memory_limit::Int)::Bool
    labels = if applicable(MetaGraphsNext.labels, graph)
        collect(MetaGraphsNext.labels(graph))
    else
        Any[]
    end

    num_kmers = length(labels)
    if num_kmers == 0
        return true
    end
    k = _label_k_length(first(labels))
    estimated_usage = estimate_memory_usage(num_kmers, k)
    return estimated_usage <= memory_limit
end

# =============================================================================
# Sequencing-technology error profiles (indel-aware correction, td-9q84 / 4a)
# =============================================================================
#
# The merged indel pair-HMM decoder (`_viterbi_correct_observation_indel`) is
# driven by an ERROR PROFILE rather than read length or a hardcoded technology:
# correct with the profile you would simulate. Each profile carries BOTH the
# technology's ABSOLUTE per-base error rate (`base_error_rate`, sourced from
# `observe(...)` — illumina 0.005, nanopore 0.10, pacbio 0.11, ultima 1e-6) AND
# the TRUE conditional insertion/deletion fractions (`P(insertion|error)` /
# `P(deletion|error)`, mirroring `observe(...)`'s per-tech error-type split; see
# `simulation.jl`, nanopore 0.30/0.30 and illumina/ultima 0.05/0.05 near line
# 1928). The kernel consumes the fractions as gap-open partitions of the error
# rate (`δ_I = error_rate·f_ins`, `δ_D = error_rate·f_del`; the remainder is the
# substitution mass carried by the existing emission term), so the corrector must
# thread the ABSOLUTE `base_error_rate` into `ViterbiCorrectionConfig.error_rate`
# for the gap masses to be correctly scaled (a config left at the 0.01 default
# would under-weight nanopore/pacbio gaps ~10×).
#
# Illumina/ultima keep their true (nonzero) conditional fractions, but their
# ABSOLUTE per-base indel rate is negligible (Illumina ≈ 0.005 × 0.10 = 5e-4;
# ultima ≈ 1e-6 × 0.05 = 5e-8), so the GATE (`profile_enables_indels`, keyed on
# the absolute rate `base_error_rate × summed fractions`) returns `false` for them
# and the corrector threads NO indel params. The substitution-only decode then
# reproduces the pre-wiring corrector byte-for-byte, so the DEFAULT `:illumina`
# profile is unchanged. Gating on the absolute rate (not the conditional
# fractions) is why illumina/ultima need no hand-tuned zeroing: their arithmetic
# decides `false` on its own, and a future tech's fractions cannot silently
# misfire the gate.

"""
ABSOLUTE per-base indel-rate threshold above which a sequencing-technology profile
enables the indel-aware pair-HMM decode (see [`profile_enables_indels`](@ref)).
The gate compares `base_error_rate × (insertion_fraction + deletion_fraction)`
against this value, so it keys on the absolute rate an indel actually occurs at,
not the conditional (given-an-error) fractions. Chosen (0.02) to sit between the
negligible absolute indel rates of substitution-dominated technologies (Illumina
≈ 5e-4, Ultima ≈ 5e-8) and the substantial rates of single-molecule long reads
(Nanopore ≈ 0.10×0.60 = 0.06, PacBio ≈ 0.11×0.80 = 0.088), with margin on both
sides so all four technologies decide correctly.
"""
const INDEL_PROFILE_THRESHOLD = 0.02

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Per-technology error-profile preset for the indel-aware read corrector.

Returns a `NamedTuple` with the technology's ABSOLUTE per-base `base_error_rate`
(sourced from `observe(...)` in `simulation.jl`: illumina 0.005, nanopore 0.10,
pacbio 0.11, ultima 1e-6), the TRUE CONDITIONAL `insertion_fraction` /
`deletion_fraction` (P(insertion|error) / P(deletion|error), mirroring
`observe(...)`'s per-tech error-type split), and the affine gap-EXTEND
probabilities (`insertion_extend_probability` / `deletion_extend_probability`) —
nanopore/PacBio indels cluster in homopolymer runs, so a one-time gap-open then a
cheaper extend models the geometric run lengths.

All four presets carry their real conditional fractions (no hand-zeroing). The
DECODE gate ([`profile_enables_indels`](@ref)) multiplies `base_error_rate` by the
summed fractions to decide, so substitution-dominated technologies (`:illumina`,
`:ultima`) — whose ABSOLUTE indel rate is negligible — gate OFF and stay on the
substitution-only oracle path, while `:nanopore`/`:pacbio` gate ON. The
`base_error_rate` is threaded into `ViterbiCorrectionConfig.error_rate` so the
kernel's gap-open masses (`δ_I = error_rate·f_ins`, `δ_D = error_rate·f_del`) are
scaled to the real per-base rate.

Supported: `:illumina`, `:nanopore`, `:pacbio`, `:ultima`.
"""
function indel_error_profile(tech::Symbol)
    if tech == :illumina
        # observe(): base 0.005, conditional insertion/deletion 0.05/0.05. Absolute
        # indel rate ≈ 5e-4 (< INDEL_PROFILE_THRESHOLD) ⇒ gate OFF ⇒ the corrector
        # threads no indel params and the decode collapses to the substitution
        # oracle (speed + byte-identity). Extend probabilities unused (gate OFF).
        return (base_error_rate = 0.005,
            insertion_fraction = 0.05, deletion_fraction = 0.05,
            insertion_extend_probability = 0.0, deletion_extend_probability = 0.0)
    elseif tech == :nanopore
        # observe(): base 0.10, split (mismatch 0.40 / insertion 0.30 / deletion
        # 0.30). Absolute indel rate ≈ 0.10 × 0.60 = 0.06 ⇒ gate ON. Indels cluster
        # in homopolymers ⇒ moderate extend.
        return (base_error_rate = 0.10,
            insertion_fraction = 0.30, deletion_fraction = 0.30,
            insertion_extend_probability = 0.30, deletion_extend_probability = 0.30)
    elseif tech == :pacbio
        # observe(): base 0.11, split (mismatch 0.20 / insertion 0.40 / deletion
        # 0.40). Absolute indel rate ≈ 0.11 × 0.80 = 0.088 ⇒ gate ON. Indel-
        # dominated ⇒ stronger extend.
        return (base_error_rate = 0.11,
            insertion_fraction = 0.40, deletion_fraction = 0.40,
            insertion_extend_probability = 0.40, deletion_extend_probability = 0.40)
    elseif tech == :ultima
        # observe(): base 1e-6, conditional insertion/deletion 0.025/0.025. Absolute
        # indel rate ≈ 5e-8 (≪ threshold) ⇒ gate OFF ⇒ substitution-only. Extend
        # probabilities unused (gate OFF).
        return (base_error_rate = 1e-6,
            insertion_fraction = 0.025, deletion_fraction = 0.025,
            insertion_extend_probability = 0.0, deletion_extend_probability = 0.0)
    else
        error("unknown sequencing technology :$(tech); expected one of " *
              ":illumina, :nanopore, :pacbio, :ultima")
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Decide whether a sequencing-technology profile enables the indel-aware decode:
`true` iff the ABSOLUTE per-base indel rate — `base_error_rate ×
(insertion_fraction + deletion_fraction)` — exceeds `threshold` (default
[`INDEL_PROFILE_THRESHOLD`]). Keying on the absolute rate (not the conditional
fractions) means `:illumina` (≈ 5e-4) and `:ultima` (≈ 5e-8) return `false` ⇒
substitution-only, while `:nanopore` (≈ 0.06) and `:pacbio` (≈ 0.088) return
`true` — each profile decides from its own arithmetic, with no hand-tuned zeroing.
"""
function profile_enables_indels(tech::Symbol;
        threshold::Float64 = INDEL_PROFILE_THRESHOLD)
    p = indel_error_profile(tech)
    return p.base_error_rate * (p.insertion_fraction + p.deletion_fraction) > threshold
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Bundle of indel pair-HMM decode parameters threaded from the assembler down to
the per-read `ViterbiCorrectionConfig`. Constructed ONLY when a profile enables
indels (see [`profile_enables_indels`](@ref)); the substitution-only path threads
`nothing`, which reproduces the pre-wiring `ViterbiCorrectionConfig` byte-for-byte
(`indel_moves` stays `false`). `base_error_rate` is the technology's ABSOLUTE
per-base error rate (from the error profile); it is threaded into
`ViterbiCorrectionConfig.error_rate` so the kernel's gap-open masses
(`δ_I = error_rate·f_ins`, `δ_D = error_rate·f_del`) are scaled to the real rate
instead of the 0.01 config default (which would under-weight nanopore/pacbio gaps
~10×). The gap-open fractions + extend probabilities come from the error profile;
the run caps + `band_width` come from the corrector tier knobs (`:scalable`
bounded, `:exhaustive` unbounded).
"""
struct IndelDecodeParams
    base_error_rate::Float64
    insertion_fraction::Float64
    deletion_fraction::Float64
    insertion_extend_probability::Float64
    deletion_extend_probability::Float64
    deletion_max_run::Int
    max_insertion_run::Int
    band_width::Union{Nothing, Int}
end

# Runtime-only frontier budget used by the scalable indel scheduler. The value is
# calibrated against warmed pair-HMM runtimes on production-length windows; graph
# vertex count and correction accuracy are deliberately not decision variables.
# The calibration script in `benchmarking/indel_frontier_runtime.jl` verifies that
# this budget rejects the small maximally-branching graph while admitting large
# linear graphs. Keep the scheduler fail-closed if that separation no longer holds.
# Runtime calibration on warmed 250/500 bp windows places a 500 bp linear walk
# at ~257k conservative probe work while the 64-vertex complete DNA k=3 graph
# exceeds 300k after only a small prefix. The 300k boundary therefore admits the
# affordable large-linear control and fails closed on the small branch explosion.
const _DEFAULT_INDEL_FRONTIER_WORK_LIMIT = 300_000
const _INDEL_FRONTIER_MAX_WINDOW = 500
const _INDEL_FRONTIER_TELEMETRY_SAMPLE_LIMIT = 64

function _empty_indel_rung_telemetry(profile_requested::Bool)::Dict{Symbol, Any}
    return Dict{Symbol, Any}(
        :profile_requested => profile_requested,
        :requested => 0,
        :attempted => 0,
        :completed => 0,
        :truncated => 0,
        :engaged => 0,
        :admitted => false,
        :admitted_windows => 0,
        :rejected_windows => 0,
        :graph_source => :raw,
        :decision_reason => profile_requested ? :not_evaluated : :profile_disabled,
        :frontier_work_limit => _DEFAULT_INDEL_FRONTIER_WORK_LIMIT,
        :frontier_metric_sample_limit => _INDEL_FRONTIER_TELEMETRY_SAMPLE_LIMIT,
        :raw_frontier_evaluated => 0,
        :cleaned_frontier_evaluated => 0,
        :bounded_windowing_forced => false,
        :substitution_decode_memory_gated => false,
        :raw_frontier_metrics => Dict{Symbol, Any}[],
        :cleaned_frontier_metrics => Dict{Symbol, Any}[],
        :graph_cleanup => Dict{String, Any}(),
    )
end

function _symbolize_indel_frontier_metric(metric::AbstractDict)::Dict{Symbol, Any}
    normalized = Dict{Symbol, Any}()
    for (key, value) in metric
        symbol_key = key isa Symbol ? key : Symbol(key)
        normalized[symbol_key] = symbol_key == :reason && value isa AbstractString ?
                                 Symbol(value) : value
    end
    return normalized
end

function _normalize_indel_rung_telemetry(
        row::AbstractDict,
)::Dict{Symbol, Any}
    profile_requested = Bool(get(
        row, :profile_requested, get(row, "profile_requested", false)))
    normalized = _empty_indel_rung_telemetry(profile_requested)
    for (key, value) in row
        symbol_key = key isa Symbol ? key : Symbol(key)
        normalized[symbol_key] = value
    end
    for key in (:graph_source, :decision_reason)
        value = get(normalized, key, nothing)
        value isa AbstractString && (normalized[key] = Symbol(value))
    end
    for key in (:raw_frontier_metrics, :cleaned_frontier_metrics)
        metrics = get(normalized, key, Any[])
        normalized[key] = Dict{Symbol, Any}[
            _symbolize_indel_frontier_metric(metric)
            for metric in metrics if metric isa AbstractDict
        ]
    end
    cleanup = get(normalized, :graph_cleanup, Dict{String, Any}())
    normalized[:graph_cleanup] = if cleanup isa AbstractDict
        Dict{String, Any}(
            string(key) => value for (key, value) in cleanup
        )
    else
        Dict{String, Any}()
    end
    return normalized
end

function _restore_indel_rung_telemetry(
        resume_data::AbstractDict,
)::Vector{Dict{Symbol, Any}}
    rows = get(resume_data, "indel_rung_telemetry", Any[])
    rows isa AbstractVector || return Dict{Symbol, Any}[]
    return Dict{Symbol, Any}[
        _normalize_indel_rung_telemetry(row)
        for row in rows if row isa AbstractDict
    ]
end

function _symbolize_iteration_stats(
        serialized::AbstractDict,
)::Dict{Symbol, Any}
    stats = Dict{Symbol, Any}()
    for (key, value) in serialized
        symbol_key = key isa Symbol ? key : Symbol(String(key))
        stats[symbol_key] = value
    end
    return stats
end

function _restore_iteration_history(
        resume_data::AbstractDict,
)::Dict{Int, Vector{Dict{Symbol, Any}}}
    serialized = get(resume_data, "iteration_history", nothing)
    serialized isa AbstractDict || throw(ArgumentError(
        "checkpoint iteration_history must be an object"))

    history = Dict{Int, Vector{Dict{Symbol, Any}}}()
    for (serialized_k, serialized_rows) in serialized
        serialized_rows isa AbstractVector || throw(ArgumentError(
            "checkpoint history for k=$serialized_k must be an array"))
        k = serialized_k isa Integer ?
            Int(serialized_k) : parse(Int, String(serialized_k))
        all(row isa AbstractDict for row in serialized_rows) ||
            throw(ArgumentError(
                "checkpoint history for k=$serialized_k contains a non-object row"))
        history[k] = Dict{Symbol, Any}[
            _symbolize_iteration_stats(row) for row in serialized_rows
        ]
    end
    return history
end

function _validate_checkpoint_history(
        k_progression::Vector{Int},
        iteration_history::Dict{Int, Vector{Dict{Symbol, Any}}},
        completed_k::Int,
        completed_iteration::Int,
)::Nothing
    isempty(k_progression) && throw(ArgumentError(
        "checkpoint k_progression must contain at least one completed rung"))
    issorted(k_progression) || throw(ArgumentError(
        "checkpoint k_progression must be sorted"))
    length(unique(k_progression)) == length(k_progression) ||
        throw(ArgumentError("checkpoint k_progression must not repeat a rung"))
    last(k_progression) == completed_k || throw(ArgumentError(
        "checkpoint current_k=$completed_k must equal the last completed rung " *
        "$(last(k_progression))"))
    Set(keys(iteration_history)) == Set(k_progression) || throw(ArgumentError(
        "checkpoint iteration_history keys must equal k_progression"))

    for history_k in k_progression
        rows = iteration_history[history_k]
        isempty(rows) && throw(ArgumentError(
            "checkpoint history for k=$history_k must not be empty"))
        observed_iterations = Int[Int(row[:iteration]) for row in rows]
        observed_iterations == collect(1:length(rows)) || throw(ArgumentError(
            "checkpoint history for k=$history_k must contain consecutive " *
            "iterations starting at 1"))
        all(Int(row[:k]) == history_k for row in rows) || throw(ArgumentError(
            "checkpoint history rows for k=$history_k contain a mismatched k"))
    end

    last_history = iteration_history[completed_k]
    Int(last(last_history)[:iteration]) == completed_iteration ||
        throw(ArgumentError(
            "checkpoint cursor k=$completed_k iteration=$completed_iteration " *
            "does not match iteration history"))
    return nothing
end

function _validated_checkpoint_cursor(
        resume_data::AbstractDict,
        completed_k::Int,
        completed_iteration::Int,
        max_k::Int,
        max_iterations_per_k::Int,
        k_schedule::Union{Nothing, Vector{Int}},
)::Tuple{Int, Int, Bool}
    next_k = Int(resume_data["next_k"])
    next_iteration = Int(resume_data["next_iteration"])
    run_complete = Bool(resume_data["run_complete"])
    expected_next_k = _next_k_in_progression(completed_k, max_k, k_schedule)

    if run_complete
        next_k == completed_k || throw(ArgumentError(
            "completed checkpoint next_k must equal current_k=$completed_k"))
        next_iteration == 1 || throw(ArgumentError(
            "completed checkpoint next_iteration must equal 1"))
        expected_next_k == completed_k || throw(ArgumentError(
            "checkpoint cannot be complete because k=$completed_k has a " *
            "remaining scheduled rung k=$expected_next_k"))
    elseif next_k == completed_k
        next_iteration == completed_iteration + 1 || throw(ArgumentError(
            "same-rung checkpoint next_iteration must equal " *
            "$(completed_iteration + 1)"))
        next_iteration <= max_iterations_per_k || throw(ArgumentError(
            "same-rung checkpoint next_iteration=$next_iteration exceeds " *
            "max_iterations_per_k=$max_iterations_per_k"))
    else
        expected_next_k != completed_k || throw(ArgumentError(
            "checkpoint specifies next_k=$next_k after the final scheduled rung"))
        next_k == expected_next_k || throw(ArgumentError(
            "checkpoint next_k=$next_k does not equal the next scheduled rung " *
            "$expected_next_k"))
        next_iteration == 1 || throw(ArgumentError(
            "next-rung checkpoint next_iteration must equal 1"))
    end
    return (next_k, next_iteration, run_complete)
end

function _profile_disabled_checkpoint_telemetry(
        k_progression::Vector{Int},
        iteration_history::Dict{Int, Vector{Dict{Symbol, Any}}},
)::Vector{Dict{Symbol, Any}}
    telemetry = Dict{Symbol, Any}[]
    for (ladder_index, history_k) in enumerate(k_progression)
        for history_row in iteration_history[history_k]
            row = _empty_indel_rung_telemetry(false)
            row[:ladder_index] = ladder_index
            row[:k] = history_k
            row[:iteration] = Int(history_row[:iteration])
            push!(telemetry, row)
        end
    end
    return telemetry
end

# =============================================================================
# Core Iterative Assembly Framework
# =============================================================================

"""
Main iterative maximum likelihood assembly function.
Processes entire read sets per iteration with complete FASTQ I/O tracking.
Enhanced with performance optimizations, caching, and progress tracking.
"""
function mycelia_iterative_assemble(input_fastq::String;
        max_k::Int = 101,
        memory_limit::Int = 32_000_000_000,
        output_dir::String = "iterative_assembly",
        max_iterations_per_k::Int = 10,
        improvement_threshold::Float64 = 0.05,
        stop_on_no_change::Bool = true,
        k_ladder::Union{Nothing, Vector{Int}} = nothing,
        n_k_rungs::Union{Nothing, Int} = nothing,
        graph_mode::Symbol = :canonical,
        verbose::Bool = true,
        enable_parallel::Bool = false,
        batch_size::Int = 10000,
        enable_checkpointing::Bool = true,
        checkpoint_interval::Int = 5,
        skip_solid::Bool = false,
        hard_window::Bool = false,
        windowed_decode::Bool = false,
        soft_em::Bool = false,
        cheap_correct::Bool = false,
        beam_width::Union{Int, Nothing} = nothing,
        calibrated_gap_threshold::Union{Float64, Nothing} = nothing,
        min_decode_k::Union{Int, Nothing} = nothing,
        decode_gate_density::Union{Float64, Nothing} = nothing,
        indel_params::Union{Nothing, IndelDecodeParams} = nothing,
        indel_schedule::Symbol = :unrestricted,
)::Dict{Symbol, Any}
    start_time = time()

    max_k >= 3 || throw(ArgumentError("max_k must be at least 3"))
    memory_limit > 0 || throw(ArgumentError("memory_limit must be positive"))
    max_iterations_per_k >= 1 || throw(ArgumentError(
        "max_iterations_per_k must be at least 1"))
    (!enable_checkpointing || checkpoint_interval >= 1) || throw(ArgumentError(
        "checkpoint_interval must be at least 1 when checkpointing is enabled"))

    indel_schedule in (:unrestricted, :frontier_budgeted) ||
        throw(ArgumentError("indel_schedule must be :unrestricted or " *
                            ":frontier_budgeted; got :$(indel_schedule)"))

    # Accumulates swallowed decode failures (structural / un-k-merizable) across
    # every k and iteration so a systematically-broken corrector that reports
    # "0 improvements" is visible in the result metadata rather than passing as
    # "nothing to fix". See CorrectorDiagnostics.
    corrector_diagnostics = CorrectorDiagnostics()

    # -- Soft-EM registry hygiene (td-e70t, v2 competing-paths + support floor) -
    # Soft-EM v2 ACTIVATES the M-step: within each k's EM loop, iteration N's
    # per-edge path responsibilities are REGISTERED onto iteration N+1's graph
    # (`register_soft_edge_weights!`, floored to each edge's own raw support so
    # supported variation is retained — td-h6w9) so `compute_edge_weight` returns
    # the probability-weighted (soft) evidence and unsupported error edges decay.
    # `improve_read_set_likelihood` owns registration for the whole read-set pass,
    # restores the exact prior registry state in `finally`, and serializes concurrent
    # soft contexts. Do not clear process-global state here: another assembly may be
    # inside its scoped interval.

    # -- Convergence + k-ladder tuning knobs (td-q70n) -------------------------
    #
    # Two literature-backed speedups for the iterative corrector:
    #
    #   1. CONVERGENCE (Musket): `stop_on_no_change=true` breaks the per-k loop the
    #      moment a pass makes 0 changes, independent of `improvement_threshold`
    #      (the old code only bounded the loop when improvement < threshold, so it
    #      never terminated on a 0-change pass when improvement_threshold==0).
    #      `max_iterations_per_k` KEEPS its default of 10 — lowering it to 2
    #      (Musket's ~2-pass heuristic) is a real accuracy/speed tradeoff that must
    #      be measured on error-rich data before becoming a default, so it is
    #      opt-in; the corrector=:iterative route in assemble_genome sets it
    #      explicitly. Note the corrector is not strictly deterministic (the
    #      statistical-resampling arm draws from the global RNG), so a 0-change
    #      pass is a practical, not guaranteed, fixed point.
    #
    #   2. K-LADDER (LoRMA 3-rung): instead of walking every prime
    #      (3,5,7,11,13,...) the caller may request a small, well-spaced ladder.
    #      `k_ladder=[...]` uses exactly those k values; `n_k_rungs=N` builds an
    #      ~N-rung geometric ladder from the initial k to `max_k`. Both default
    #      to `nothing`, which preserves the original prime-by-prime progression
    #      byte-for-byte, so existing callers are unchanged unless they opt in.

    if verbose
        println("Starting Mycelia Iterative Maximum Likelihood Assembly")
        println("Input FASTQ: $input_fastq")
        println("Output directory: $output_dir")
        println("Memory limit: $(memory_limit ÷ 1_000_000_000) GB")
        println("Max k-mer size: $max_k")
        println("Graph mode: $graph_mode")
        parallel_status = enable_parallel ?
                          "enabled ($(Threads.nthreads()) threads)" : "disabled"
        println("Parallel processing: $parallel_status")
        println("Batch size: $batch_size")
        checkpoint_status = enable_checkpointing ?
                            "enabled (every $checkpoint_interval iterations)" :
                            "disabled"
        println("Checkpointing: $checkpoint_status")
    end

    # Create output directory structure
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Create subdirectories for organization
    checkpoints_dir = joinpath(output_dir, "checkpoints")
    graphs_dir = joinpath(output_dir, "graphs")
    progress_dir = joinpath(output_dir, "progress")

    if enable_checkpointing
        mkpath(checkpoints_dir)
        mkpath(graphs_dir)
        mkpath(progress_dir)
    end

    # Check for existing checkpoint to resume from
    checkpoint_file = joinpath(checkpoints_dir, "latest_checkpoint.json")
    resume_data = nothing
    runtime_offset = 0.0

    if enable_checkpointing && isfile(checkpoint_file)
        try
            resume_data = JSON.parsefile(checkpoint_file)
            if verbose
                checkpoint_k = resume_data["current_k"]
                checkpoint_iteration = resume_data["current_iteration"]
                println(
                    "Found existing checkpoint. Resume from k=$checkpoint_k, " *
                    "iteration=$checkpoint_iteration",
                )
            end
        catch e
            if verbose
                println("Warning: Could not load checkpoint file: $e")
            end
        end
    end

    # Initialize or resume from checkpoint
    completed_checkpoint_k = nothing
    completed_checkpoint_iteration = nothing
    resume_iteration = nothing
    resume_existing_rung = false
    resume_run_complete = false
    if resume_data !== nothing
        # Resume from checkpoint
        completed_checkpoint_k = Int(resume_data["current_k"])
        completed_checkpoint_iteration = Int(resume_data["current_iteration"])
        k = completed_checkpoint_k
        k_progression = Int[Int(value) for value in resume_data["k_progression"]]
        iteration_history = _restore_iteration_history(resume_data)
        total_improvements = Int(resume_data["total_improvements"])
        current_fastq_file = String(resume_data["current_fastq_file"])
        runtime_offset = Float64(get(resume_data, "runtime_so_far", 0.0))
        isfinite(runtime_offset) && runtime_offset >= 0.0 || throw(ArgumentError(
            "checkpoint runtime_so_far must be finite and nonnegative"))
        _validate_checkpoint_history(
            k_progression,
            iteration_history,
            completed_checkpoint_k,
            completed_checkpoint_iteration,
        )

        if verbose
            println("Loaded checkpoint after k=$completed_checkpoint_k, " *
                    "iteration=$completed_checkpoint_iteration, " *
                    "total_improvements=$total_improvements")
        end
    else
        # Initialize fresh run
        if verbose
            println("Reading initial FASTQ file...")
        end
        initial_reads = collect(FASTX.FASTQ.Reader(open(input_fastq)))
        k = find_initial_k(initial_reads)  # Reuse from intelligent-assembly.jl
        k <= max_k || throw(ArgumentError(
            "max_k=$max_k is below the detected initial k-mer size k=$k"))

        if verbose
            println("Initial k-mer size: $k (prime: $(Primes.isprime(k)))")
            println("Total reads: $(length(initial_reads))")
        end

        # Track progress across all iterations
        k_progression = Int[]
        iteration_history = Dict{Int, Vector{Dict{Symbol, Any}}}()
        total_improvements = 0
        current_fastq_file = input_fastq
    end

    # Build the k-mer progression schedule. `nothing` => legacy prime-by-prime
    # walk via next_prime_k; a Vector{Int} => explicit / coarse LoRMA-style ladder.
    schedule_initial_k = isempty(k_progression) ? k : first(k_progression)
    k_schedule = build_k_ladder(
        schedule_initial_k, max_k;
        k_ladder = k_ladder,
        n_k_rungs = n_k_rungs,
    )
    if verbose && k_schedule !== nothing
        println("K-mer ladder (coarse progression): $(k_schedule)")
    end

    if resume_data !== nothing
        cursor_keys = ("next_k", "next_iteration", "run_complete")
        cursor_fields_present = Bool[
            haskey(resume_data, key) for key in cursor_keys
        ]
        if all(cursor_fields_present)
            k, resume_iteration, resume_run_complete =
                _validated_checkpoint_cursor(
                    resume_data,
                    something(completed_checkpoint_k),
                    something(completed_checkpoint_iteration),
                    max_k,
                    max_iterations_per_k,
                    k_schedule,
                )
            resume_existing_rung = !resume_run_complete &&
                                   k == completed_checkpoint_k
        elseif any(cursor_fields_present)
            throw(ArgumentError(
                "checkpoint must contain all or none of next_k, next_iteration, " *
                "and run_complete"))
        else
            # Backward-compatible cursor inference for checkpoints written before
            # next-pass fields were persisted. Reproduce the normal convergence
            # decision from the completed history rather than replaying the pass.
            completed_history = iteration_history[completed_checkpoint_k]
            completed_stats = last(completed_history)
            completed_improvements = Int(completed_stats[:improvements_made])
            completed_reads = Int(completed_stats[:total_reads])
            continue_current_k =
                completed_checkpoint_iteration < max_iterations_per_k &&
                !(stop_on_no_change && completed_improvements == 0) &&
                sufficient_improvements(
                    completed_improvements,
                    completed_reads,
                    improvement_threshold;
                    iteration_history = completed_history,
                )
            if continue_current_k
                k = completed_checkpoint_k
                resume_iteration = completed_checkpoint_iteration + 1
                resume_existing_rung = true
            else
                k = _next_k_in_progression(
                    completed_checkpoint_k, max_k, k_schedule)
                resume_iteration = 1
                resume_run_complete = k == completed_checkpoint_k
            end
        end
        resume_existing_rung && soft_em && throw(ArgumentError(
            "cannot resume soft_em=true within k=$k because the prior soft-edge " *
            "accumulator is not checkpointed; resume from the next k rung"))
    end

    # -- Low-k decode gating (td-9h5r) -----------------------------------------
    # The finding (PR #370 profile): the per-read graph-Viterbi decode is 57-78% of
    # every :scalable iteration — the dominant runtime term — yet at LOW k the graph
    # is dense (nearly every k-mer sits on a bubble/repeat vertex), so the
    # hard-window gate cannot discriminate and flags essentially EVERY read for a
    # full whole-read decode. That is correction Stage 0's cheap linear pass already
    # largely did. Two composable gates cut that wasted volume; BOTH still BUILD the
    # graph and run Stage 0 cheap correction + skip-solid, deferring only the
    # expensive whole-read Viterbi to higher, selective rungs:
    #
    #   (1) ADAPTIVE (default): at a rung whose post-Stage-0 hard-window decode
    #       fraction would be >= `decode_gate_density`, the gate is non-
    #       discriminating (it would decode ~everything → no volume reduction), so
    #       the whole decode is skipped for that pass. This fires EXACTLY on the
    #       dense-low-k pathology and is a NO-OP where the gate genuinely filters
    #       (so a load-bearing selective low-k decode is preserved). Measured
    #       post-Stage-0 inside improve_read_set_likelihood (raw-read density is
    #       uninformative — Stage 0 clears most apparent hardness first).
    #   (2) EXPLICIT floor: `min_decode_k` hard-gates every rung `k < min_decode_k`
    #       (the LoRMA-style "start the graph-decode ladder higher" lever). Off by
    #       default; deterministic, for tuning/tests.
    #
    # Both are gated behind the `:scalable` tier (only active when `hard_window` is
    # true), so the `:exhaustive` tier (hard_window=false) is BYTE-IDENTICAL.
    effective_min_decode_k = hard_window ? min_decode_k : nothing
    # Adaptive-gate density threshold (fraction of reads the hard-window gate would
    # still decode after Stage 0 above which the gate counts as non-discriminating).
    # Default 0.90 on :scalable; `nothing` (or the exhaustive tier) disables it.
    effective_decode_gate_density = hard_window ?
                                    (decode_gate_density === nothing ?
                                     _DEFAULT_DECODE_GATE_DENSITY : decode_gate_density) :
                                    nothing
    if verbose && (effective_min_decode_k !== nothing ||
        effective_decode_gate_density !== nothing)
        println("Low-k decode gating (td-9h5r): " *
                (effective_min_decode_k !== nothing ?
                 "explicit floor min_decode_k=$(effective_min_decode_k); " : "") *
                (effective_decode_gate_density !== nothing ?
                 "adaptive gate-off when post-Stage-0 decode fraction >= " *
                 "$(effective_decode_gate_density). " : "") *
                "Gated rungs rely on Stage 0 cheap correction + skip-solid.")
    end
    # Telemetry: the k-rungs whose per-read decode was gated OFF (low-k skip).
    decode_gated_rungs = if resume_data === nothing
        Int[]
    else
        Int[Int(value) for value in
            get(resume_data, "decode_gated_rungs", Any[])]
    end

    # Hard-window skip telemetry (td-nn6l): the fraction of reads the hard-window
    # gate passed through WITHOUT a decode, recorded PER PASS (across every k and
    # iteration) so the run metadata can surface min/mean/max, not just the last
    # value (FIX 6). Empty on the :exhaustive tier (no gate ⇒ no skips).
    skip_fractions = if resume_data === nothing
        Float64[]
    else
        Float64[Float64(value) for value in
            get(resume_data, "skip_fractions", Any[])]
    end

    # Stage 0 cheap-correction telemetry (td-bjnt): bases fixed by the linear
    # k-mer-spectrum pass, recorded PER PASS. Empty/zero when cheap_correct is off
    # (:exhaustive), so the exhaustive tier is unaffected.
    cheap_correction_counts = if resume_data === nothing
        Int[]
    else
        Int[Int(value) for value in
            get(resume_data, "cheap_correction_counts", Any[])]
    end

    # One entry per (ladder index, k, iteration), including profile-disabled and
    # classifier-rejected passes. Repeated iterations at the same k are retained;
    # no k-based deduplication is allowed because requested/attempted/outcome
    # transitions are pass-local runtime evidence.
    indel_rung_telemetry = if resume_data === nothing
        Dict{Symbol, Any}[]
    elseif haskey(resume_data, "indel_rung_telemetry")
        _restore_indel_rung_telemetry(resume_data)
    elseif indel_params === nothing
        _profile_disabled_checkpoint_telemetry(k_progression, iteration_history)
    else
        throw(ArgumentError(
            "checkpoint predates indel rung telemetry; cannot recover historical " *
            "requested/attempted/completed/truncated/engaged counters for an " *
            "indel-enabled resume"))
    end
    if resume_data !== nothing
        expected_passes = sum(length(rows) for rows in values(iteration_history))
        length(indel_rung_telemetry) == expected_passes || throw(ArgumentError(
            "checkpoint indel telemetry has $(length(indel_rung_telemetry)) " *
            "rows; expected $expected_passes completed passes"))
        telemetry_cursor = Tuple{Int, Int, Int}[
            (
                Int(row[:ladder_index]),
                Int(row[:k]),
                Int(row[:iteration]),
            ) for row in indel_rung_telemetry
        ]
        history_cursor = Tuple{Int, Int, Int}[]
        for (ladder_index, history_k) in enumerate(k_progression)
            for row in get(iteration_history, history_k, Dict{Symbol, Any}[])
                push!(history_cursor,
                    (ladder_index, history_k, Int(row[:iteration])))
            end
        end
        telemetry_cursor == history_cursor || throw(ArgumentError(
            "checkpoint indel telemetry cursor does not match iteration history"))
    end
    # The aggregate diagnostic counters cover the whole resumed run, not only the
    # post-resume suffix. Requested work has no atomic diagnostic counterpart and
    # is derived from the preserved rows at finalization.
    Threads.atomic_add!(corrector_diagnostics.indel_attempts,
        sum(get(row, :attempted, 0) for row in indel_rung_telemetry; init = 0))
    Threads.atomic_add!(corrector_diagnostics.indel_decodes,
        sum(get(row, :completed, 0) for row in indel_rung_telemetry; init = 0))
    Threads.atomic_add!(corrector_diagnostics.truncated_decodes,
        sum(get(row, :truncated, 0) for row in indel_rung_telemetry; init = 0))
    Threads.atomic_add!(corrector_diagnostics.indel_engaged,
        sum(get(row, :engaged, 0) for row in indel_rung_telemetry; init = 0))
    if resume_data !== nothing
        saved_diagnostics = get(
            resume_data, "corrector_diagnostics", Dict{String, Any}())
        saved_diagnostics isa AbstractDict || throw(ArgumentError(
            "checkpoint corrector_diagnostics must be an object"))
        Threads.atomic_add!(corrector_diagnostics.structural_errors,
            Int(get(saved_diagnostics, "structural", 0)))
        Threads.atomic_add!(corrector_diagnostics.unkmerizable_reads,
            Int(get(saved_diagnostics, "unkmerizable", 0)))
        Threads.atomic_add!(corrector_diagnostics.gate_skipped,
            Int(get(saved_diagnostics, "gate_skipped", 0)))
        Threads.atomic_add!(corrector_diagnostics.trace_contract_errors,
            Int(get(saved_diagnostics, "trace_contract_errors", 0)))
        Threads.atomic_add!(corrector_diagnostics.window_divergences,
            Int(get(saved_diagnostics, "window_divergences", 0)))
    end

    # --- Final-pass graph reuse (td-04tb) -----------------------------------
    # The corrector builds a qualmer graph FROM SCRATCH each pass; the LAST pass at
    # the largest k builds the graph whose input reads become the final corrected
    # read set — but ONLY when that pass makes 0 improvements (a converged, no-change
    # pass). In that case the graph built from `current_reads` is byte-identical to
    # the graph a downstream re-assembly would build from the corrected reads, so we
    # hand it back and `_assemble_with_iterative_corrector` reuses it instead of
    # rebuilding, saving the redundant build_qualmer_graph (~22-36% of the iterative
    # arm). `final_pass_graph_reusable` is set true ONLY on a 0-improvement pass;
    # when the last pass changed reads the graph is stale and the caller rebuilds.
    # These are overwritten every pass, so after the loop they reflect the LAST
    # pass at the LARGEST k processed (the pass that produced the final FASTQ).
    final_pass_graph = nothing
    final_pass_graph_k = 0
    final_pass_graph_reusable = false
    completed_runtime = nothing

    # Main k-mer progression loop
    while !resume_run_complete && k <= max_k
        if verbose
            println("\n" * "="^60)
            println("PROCESSING K-MER SIZE: $k (prime: $(Primes.isprime(k)))")
            println("="^60)
        end

        if resume_existing_rung
            iteration = something(resume_iteration, 1)
            improvements_this_k = sum(
                Int(row[:improvements_made]) for row in iteration_history[k])
            resume_existing_rung = false
            resume_iteration = nothing
        else
            push!(k_progression, k)
            iteration_history[k] = Dict{Symbol, Any}[]
            iteration = 1
            improvements_this_k = 0
        end
        # Declared in the outer-k scope so it is visible after the while loop
        # (it is assigned inside the loop; Julia loop scope would otherwise leave
        # it undefined at the post-loop "Final improvement rate" line — a latent
        # crash that prevented the pipeline from completing past the first k).
        current_reads = FASTX.FASTQ.Record[]

        # Soft-EM M-step memory (td-e70t v2): the PREVIOUS EM iteration's soft
        # edge-weight accumulator, registered onto THIS iteration's freshly-built
        # graph so the decode consumes the floored, probability-weighted edges.
        # Reset to `nothing` at each new k — accumulator keys are (src,dst) k-mer
        # label tuples, so weights from a different k never match and must not
        # carry over.
        prev_soft_weights = nothing

        # Iterative improvement loop for current k
        while iteration <= max_iterations_per_k
            if verbose
                println("\n--- Iteration $iteration for k=$k ---")
            end

            iteration_start = time()
            timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")

            # Read in entire FASTQ set for this iteration
            if verbose
                println("Reading FASTQ file: $current_fastq_file")
            end
            current_reads = collect(FASTX.FASTQ.Reader(open(current_fastq_file)))

            # Build qualmer graph from current read set
            if verbose
                println("Building qualmer graph with k=$k...")
            end
            graph = Mycelia.Rhizomorph.build_qualmer_graph(current_reads, k; mode = graph_mode)
            num_kmers = length(MetaGraphsNext.labels(graph))

            if verbose
                println("Graph built: $num_kmers unique k-mers")
            end

            # Check memory usage
            if !check_memory_limits(graph, memory_limit)
                if verbose
                    println("Memory limit reached. Stopping at k=$k")
                end
                break
            end

            # --- Soft-EM E-step accumulator (td-e70t, v2 COMPETING-PATHS) -----
            # A fresh accumulator tallies THIS pass's per-edge responsibilities. In
            # v2 the E-step enumerates COMPETING candidate paths per decoded read
            # (the observed read path vs a consensus alternative re-routed through
            # the best-supported sibling) and splits the responsibility across them
            # by a stable softmax over their normalized-transition log-probabilities.
            # This accumulator is then REGISTERED (support-floored) onto the NEXT
            # iteration's graph in the M-step below, closing the loop that lets
            # unsupported error edges decay while supported variation is retained.
            # Only allocated under soft-EM so :exhaustive threads `nothing` and is
            # byte-for-byte unchanged.
            current_soft_weights = soft_em ?
                                   Mycelia.Rhizomorph.SoftEdgeWeightAccumulator() : nothing

            # --- Hard-window gating (td-nn6l) --------------------------------
            # Restrict decoding to reads that touch a "hard" vertex (bubble
            # entry/exit/interior, high-out-degree/repeat-like, or weak/non-solid
            # k-mer). Easy reads pass through untouched, cutting decode volume
            # ~85-95% on clean data. Off unless `hard_window=true` (:scalable).
            # Mode-agnostic (td-nt69): `_hard_vertex_set` operates on graph
            # vertices/edges and `should_decode_read` matches reads in their observed
            # orientation, so the gate works on both :canonical and :doublestrand
            # (:scalable now runs :doublestrand). Only :singlestrand — which has no
            # RC vertices and is not a corrector target — is excluded.
            # Low-k decode gating (td-9h5r) — EXPLICIT floor. Below
            # `effective_min_decode_k` the per-read graph-Viterbi decode is skipped
            # for EVERY read this pass. Stage 0 cheap correction + skip-solid still
            # run below. Only active on the :scalable tier (hard_window=true), so
            # :exhaustive is byte-identical. The ADAPTIVE gate (density-based) is
            # decided inside improve_read_set_likelihood — it needs the post-Stage-0
            # decode fraction — and is reported back via `pass_decode_gated`.
            explicit_floor_gated = effective_min_decode_k !== nothing &&
                                   k < effective_min_decode_k
            # Building the hard-vertex set is pointless when the explicit floor
            # already gates the decode off; the adaptive gate needs it, though, so
            # build it whenever the floor is NOT gating and hard_window is on.
            hard_vertices = (hard_window && !explicit_floor_gated &&
                             graph_mode != :singlestrand) ?
                            _hard_vertex_set(graph, k) : nothing
            if hard_window && !explicit_floor_gated && graph_mode == :singlestrand
                @warn "hard_window gating is not supported for graph_mode=:singlestrand; " *
                      "disabling (decoding all non-solid reads)." graph_mode maxlog = 1
            end

            # Process each read for likelihood improvement with performance optimizations
            if verbose
                println("Processing reads for likelihood improvements...")
            end
            # `beam_width` defaults to the size-aware auto-beam (exact where
            # tractable, bounded on large reads) so the production assembler cannot
            # OOM-crash (td-63qy); :exhaustive passes typemax(Int) to force exact.
            # `corrector_diagnostics` accumulates swallowed decode failures across
            # every k and iteration and is surfaced in the result metadata.
            # --- Soft-EM M-step consumption (td-e70t v2) ---------------------
            # Register the PREVIOUS EM iteration's soft edge weights (support-
            # floored) onto THIS freshly-built graph so `compute_edge_weight` — and
            # thus the Viterbi transition scoring AND the competing-path enumeration
            # in the E-step — consumes the decayed, probability-weighted edges
            # instead of raw counts. The floor holds every >= MIN_SUPPORT edge at
            # its raw coverage (real variation preserved), so only unsupported error
            # edges decay. The read-set helper owns a scoped raw+cleaned registration
            # interval and restores the exact prior state even on error. On
            # :exhaustive, `soft_em` is false: the scope still serializes against a
            # concurrent soft owner, but it registers no weights.
            updated_reads = current_reads
            improvements_made = 0
            pass_skip_fraction = 0.0
            pass_cheap_corrections = 0
            pass_decode_gated = false
            pass_indel_telemetry =
                _empty_indel_rung_telemetry(indel_params !== nothing)
            updated_reads,
            improvements_made,
            pass_skip_fraction,
            pass_cheap_corrections,
            pass_decode_gated = improve_read_set_likelihood(
                current_reads, graph, k,
                verbose = verbose,
                batch_size = batch_size,
                enable_parallel = enable_parallel,
                graph_mode = graph_mode,
                skip_solid = skip_solid,
                cheap_correct = cheap_correct,
                beam_width = beam_width,
                calibrated_gap_threshold = calibrated_gap_threshold,
                soft_weights = current_soft_weights,
                prior_soft_weights = prev_soft_weights,
                hard_vertices = hard_vertices,
                windowed_decode = windowed_decode,
                decode_enabled = !explicit_floor_gated,
                decode_gate_density = effective_decode_gate_density,
                diagnostics = corrector_diagnostics,
                indel_params = indel_params,
                indel_schedule = indel_schedule,
                rung_telemetry = pass_indel_telemetry,
                memory_limit = memory_limit,
            )
            push!(skip_fractions, pass_skip_fraction)
            push!(cheap_correction_counts, pass_cheap_corrections)
            pass_indel_telemetry[:ladder_index] = length(k_progression)
            pass_indel_telemetry[:k] = k
            pass_indel_telemetry[:iteration] = iteration
            push!(indel_rung_telemetry, pass_indel_telemetry)
            # Record a rung whose per-read decode was gated OFF (explicit floor OR
            # adaptive density gate) — telemetry surfaced in the run metadata.
            if (explicit_floor_gated || pass_decode_gated) &&
               (isempty(decode_gated_rungs) || last(decode_gated_rungs) != k)
                push!(decode_gated_rungs, k)
            end
            # Carry this iteration's soft edge memory forward: it becomes the next
            # EM iteration's M-step input (registered onto the next graph). `nothing`
            # on :exhaustive (no soft-EM) so that tier stays byte-identical.
            prev_soft_weights = soft_em ? current_soft_weights : nothing

            # Calculate iteration metrics
            iteration_time = time() - iteration_start
            improvement_rate = improvements_made / length(current_reads)

            iteration_stats = Dict(
                :iteration => iteration,
                :k => k,
                :timestamp => timestamp,
                :improvements_made => improvements_made,
                :total_reads => length(current_reads),
                :improvement_rate => improvement_rate,
                :runtime_seconds => iteration_time,
                :memory_kmers => num_kmers
            )

            push!(iteration_history[k], iteration_stats)
            improvements_this_k += improvements_made
            total_improvements += improvements_made

            if verbose
                println("Improvements made: $improvements_made ($(round(improvement_rate * 100, digits=2))%)")
                println("Iteration runtime: $(round(iteration_time, digits=2)) seconds")
            end

            # Write out complete updated read set to timestamped FASTQ
            output_file = joinpath(output_dir, "reads_k$(k)_iter$(iteration)_$(timestamp).fastq")
            write_fastq(records = updated_reads, filename = output_file)  # Use existing function

            if verbose
                println("Wrote $(length(updated_reads)) reads to $output_file")
            end

            # Capture this pass's graph for potential final-pass reuse (td-04tb).
            # `graph` was built from `current_reads` this pass; `output_file`
            # (== updated_reads) is this pass's corrected output. When the pass made
            # 0 improvements, updated_reads == current_reads content, so `graph` is
            # byte-identical to a rebuild from the corrected reads and is REUSABLE.
            # Overwritten every pass ⇒ after the loop this holds the last pass at the
            # largest k, i.e. the pass whose FASTQ finalize selects as the result.
            final_pass_graph = graph
            final_pass_graph_k = k
            final_pass_graph_reusable = (improvements_made == 0)

            # Persist the NEXT-pass cursor, not merely the pass that just
            # completed. This is the same decision used below during uninterrupted
            # execution, so resume never replays a pass or shifts rung telemetry.
            continue_current_k =
                iteration < max_iterations_per_k &&
                !(stop_on_no_change && improvements_made == 0) &&
                sufficient_improvements(
                    improvements_made,
                    length(current_reads),
                    improvement_threshold;
                    iteration_history = iteration_history[k],
                )
            checkpoint_next_k = continue_current_k ?
                                k : _next_k_in_progression(k, max_k, k_schedule)
            checkpoint_next_iteration = continue_current_k ? iteration + 1 : 1
            checkpoint_run_complete = !continue_current_k && checkpoint_next_k == k
            current_fastq_file = output_file
            pass_runtime_so_far = runtime_offset + (time() - start_time)
            if checkpoint_run_complete
                completed_runtime = pass_runtime_so_far
            end

            # Create checkpoint if enabled and at checkpoint interval
            if enable_checkpointing &&
               (iteration % checkpoint_interval == 0 || checkpoint_run_complete)
                checkpoint_data = Dict(
                    "current_k" => k,
                    "current_iteration" => iteration,
                    "next_k" => checkpoint_next_k,
                    "next_iteration" => checkpoint_next_iteration,
                    "run_complete" => checkpoint_run_complete,
                    "k_progression" => k_progression,
                    "iteration_history" => iteration_history,
                    "total_improvements" => total_improvements,
                    "indel_rung_telemetry" => indel_rung_telemetry,
                    "skip_fractions" => skip_fractions,
                    "cheap_correction_counts" => cheap_correction_counts,
                    "decode_gated_rungs" => decode_gated_rungs,
                    "corrector_diagnostics" => Dict(
                        "structural" => corrector_diagnostics.structural_errors[],
                        "unkmerizable" =>
                            corrector_diagnostics.unkmerizable_reads[],
                        "gate_skipped" => corrector_diagnostics.gate_skipped[],
                        "trace_contract_errors" =>
                            corrector_diagnostics.trace_contract_errors[],
                        "window_divergences" =>
                            corrector_diagnostics.window_divergences[],
                    ),
                    "current_fastq_file" => output_file,
                    "timestamp" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
                    "runtime_so_far" => pass_runtime_so_far,
                )

                try
                    open(checkpoint_file, "w") do f
                        JSON.print(f, checkpoint_data, 2)
                    end

                    # Also save progress summary
                    progress_file = joinpath(progress_dir, "progress_k$(k)_iter$(iteration).json")
                    open(progress_file, "w") do f
                        JSON.print(f,
                            Dict(
                                "k" => k,
                                "iteration" => iteration,
                                "improvements_this_iteration" => improvements_made,
                                "total_improvements" => total_improvements,
                                "improvement_rate" => improvement_rate,
                                "runtime_seconds" => iteration_time,
                                "timestamp" => timestamp
                            ),
                            2)
                    end

                    if verbose
                        println("Checkpoint saved: k=$k, iteration=$iteration")
                    end
                catch e
                    if verbose
                        println("Warning: Could not save checkpoint: $e")
                    end
                end
            end

            if !continue_current_k
                if verbose
                    if stop_on_no_change && improvements_made == 0
                        println("No changes this pass (0 improvements) — " *
                                "converged at k=$k")
                    else
                        println("Insufficient improvements, iteration cap, or " *
                                "convergence detected. Moving to next k-mer size")
                    end
                end
                break
            elseif verbose
                println("Sufficient improvements detected. Continuing with k=$k")
            end
            iteration += 1
        end

        isempty(iteration_history[k]) && throw(ArgumentError(
            "no correction pass completed at k=$k; the graph exceeded the " *
            "memory limit before a checkpointable result was produced"))

        if verbose
            println("\nCompleted k=$k processing:")
            println("  Total iterations: $(length(iteration_history[k]))")
            println("  Total improvements: $improvements_this_k")
            println("  Final improvement rate: $(round(improvements_this_k / length(current_reads), digits=4))")
        end

        # Move to next scheduled k-mer size (coarse ladder or prime progression)
        next_k = _next_k_in_progression(k, max_k, k_schedule)
        if next_k == k
            if verbose
                println("No larger k-mer size available. Stopping at k=$k")
            end
            break
        end
        k = next_k
    end

    # Calculate total runtime
    total_runtime = if resume_run_complete
        runtime_offset
    elseif completed_runtime === nothing
        runtime_offset + (time() - start_time)
    else
        completed_runtime
    end

    if verbose
        println("\n" * "="^60)
        println("ITERATIVE ASSEMBLY COMPLETE")
        println("="^60)
        println("Total runtime: $(round(total_runtime, digits=2)) seconds")
        println("K-mer sizes processed: $(sort(k_progression))")
        println("Total improvements: $total_improvements")
    end

    # Finalize iterative assembly
    return finalize_iterative_assembly(
        output_dir, k_progression, iteration_history, total_runtime,
        verbose = verbose, diagnostics = corrector_diagnostics,
        skip_fractions = skip_fractions,
        cheap_correction_counts = cheap_correction_counts,
        hard_window = hard_window, windowed_decode = windowed_decode,
        soft_em = soft_em, cheap_correct = cheap_correct,
        min_decode_k = effective_min_decode_k,
        decode_gate_density = effective_decode_gate_density,
        decode_gated_rungs = decode_gated_rungs,
        indel_schedule = indel_schedule,
        indel_rung_telemetry = indel_rung_telemetry,
        # Final-pass graph reuse (td-04tb): hand the corrector's last-pass qualmer
        # graph back so a converged run's re-assembly can skip a redundant rebuild.
        final_pass_graph = final_pass_graph,
        final_pass_graph_k = final_pass_graph_k,
        final_pass_graph_mode = graph_mode,
        final_pass_graph_reusable = final_pass_graph_reusable)
end

# =============================================================================
# Read Likelihood Improvement Functions
# =============================================================================

# ------------------------------------------------------------------------------
# Stage 0 skip support (td-1do7): classify k-mers once per graph so that reads
# with no weak k-mers can skip the expensive per-read Viterbi decode. This is the
# volume reduction that motivated the staged-correction architecture — the
# corrector should only decode reads that actually need correction.
# ------------------------------------------------------------------------------

"""
    _solid_kmer_set(graph; classifier) -> Set

Classify every k-mer once and return the set of SOLID (real-genomic) canonical
k-mer labels. Defaults to the auto-threshold `MixtureModelClassifier` (a top
default in the Stage 0 comparison matrix, coverage-based so it needs no per-arm
tuning). Returns the EMPTY set (⇒ `_read_is_all_solid` skips nothing ⇒ correct
every read) when classification cannot run — no labels or no dataset evidence —
so the skip is never applied on evidence it could not compute (review C1). The
conservative default is "decode everything", NOT "skip everything": returning all
labels would mark error k-mers (graph vertices at coverage 1) as solid and
silently skip reads that needed correction.
"""
function _solid_kmer_set(graph;
        classifier = Mycelia.Rhizomorph.MixtureModelClassifier())
    labels = collect(MetaGraphsNext.labels(graph))
    isempty(labels) && return Set{eltype(labels)}()
    dataset_id = _rhizomorph_first_dataset_id(graph[first(labels)])
    if dataset_id === nothing
        @warn "skip_solid: graph has no dataset evidence; classification skipped, correcting all reads (no skip)."
        return Set{eltype(labels)}()   # empty ⇒ skip nothing ⇒ correct everything
    end
    classification = Mycelia.Rhizomorph.classify_kmers(classifier, graph;
        dataset_id = dataset_id)
    solid = Set{eltype(labels)}()
    for (label, is_solid) in classification.verdicts
        is_solid && push!(solid, label)
    end
    return solid
end

# Mode-aware lookup key for classification/gating (td-nt69). `solid_kmers` /
# `hard_vertices` are graph VERTEX LABELS, so the key that a read's observed k-mer
# maps to depends on how the graph was built:
#   :canonical    → labels are canonical k-mers, so canonicalize the read k-mer.
#   :doublestrand → both a k-mer and its RC are SEPARATE vertices, so the read
#                   traverses its observed-orientation k-mer; use it as-is.
#   :singlestrand → no RC vertices; use the observed k-mer as-is.
# Classification is coverage-based and therefore mode-agnostic: under
# :doublestrand each strand-vertex still separates real (≈half-coverage) from
# error (coverage-1) k-mers, so the solid/hard sets remain valid.
@inline function _lookup_key(kmer, graph_mode::Symbol)
    return graph_mode == :canonical ? BioSequences.canonical(kmer) : kmer
end

"""
    _read_is_all_solid(read, k, solid_kmers; graph_mode) -> Bool

True iff every k-mer of `read` (resolved to its graph label via `_lookup_key`) is
in `solid_kmers` (⇒ no weak region ⇒ the read needs no correction and can skip the
decode). A read with no k-mers (shorter than k, or fully ambiguous) returns
`false` so it is never skipped on the basis of absent evidence. Works on both
`:canonical` and `:doublestrand` graphs (td-nt69).
"""
function _read_is_all_solid(read::FASTX.FASTQ.Record, k::Int, solid_kmers::AbstractSet;
        graph_mode::Symbol = :canonical)
    isempty(solid_kmers) && return false
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    saw_kmer = false
    for (kmer, _) in Kmers.UnambiguousDNAMers{k}(sequence)
        saw_kmer = true
        (_lookup_key(kmer, graph_mode) in solid_kmers) || return false
    end
    return saw_kmer
end

# ------------------------------------------------------------------------------
# Stage 0 CHEAP k-mer-spectrum correction (td-bjnt): fix simple single-base
# substitution errors with a LINEAR scan, BFC/Lighter/Bloocoo-style, BEFORE the
# expensive per-read graph Viterbi. A WEAK (low-coverage, non-solid) k-mer flanked
# by SOLID k-mers is the signature of a single-base error; if exactly one base
# substitution turns the whole weak run solid (a UNIQUE solid neighbor), apply it.
# Ambiguous cases (no fix, or >1 candidate — e.g. a balanced heterozygous site
# where both alleles are real) are LEFT untouched for the graph decode. Reserving
# graph Viterbi for genuine ambiguity is the critical-path lever: on err=0.01 data
# ~78% of reads carry an error k-mer, so cheaply clearing the simple ones collapses
# the graph-decode fraction toward the true bubble/repeat ~5–15%.
# Only invoked on the :scalable tier (`cheap_correct=true`, canonical graph); the
# :exhaustive path never calls it and is byte-identical.
# ------------------------------------------------------------------------------

# Graph-label lookup key for the k-mer at 1-based read position `i` (bases
# [i, i+k-1]) of an ACGT-only char vector, comparable against the `_solid_kmer_set`
# labels. Mode-aware (td-nt69): canonicalized under :canonical, observed
# orientation under :doublestrand / :singlestrand. The read is pre-screened for
# ambiguity by the caller so this construction cannot throw.
@inline function _lookup_kmer_at(chars::Vector{Char}, k::Int, i::Int, graph_mode::Symbol)
    return _lookup_key(Kmers.DNAKmer{k}(String(@view chars[i:(i + k - 1)])), graph_mode)
end

"""
    _stage0_correct_read(read, k, solid_kmers) -> (record, n_corrections)

Cheaply correct simple single-substitution errors in ONE read via the k-mer
spectrum. Walk the read's k-mers; for each maximal run of consecutive WEAK
(non-solid) k-mers that is flanked by SOLID k-mers on BOTH sides (an interior
error) and is short enough to be a single substitution (run length ≤ k), search
the base position(s) shared by every k-mer in the run for a substitution that
makes EVERY k-mer overlapping that position solid. Apply the fix only when EXACTLY
ONE (position, base) candidate qualifies — a unique solid neighbor. Zero
candidates (a real but low-coverage allele with no solid neighbor) or multiple
candidates (a balanced variant with two solid alleles, or genuine ambiguity) are
left untouched for the graph decode, so Stage 0 never collapses real variation.

Reads shorter than `k`, reads with ambiguous bases (e.g. `N`), and reads with no
weak run are returned unchanged with a 0 count. Returns the (possibly rewritten)
record and the number of bases corrected. `graph_mode` selects the solid-set
lookup key (canonicalized under :canonical, observed orientation otherwise —
td-nt69).
"""
function _stage0_correct_read(read::FASTX.FASTQ.Record, k::Int, solid_kmers::AbstractSet;
        graph_mode::Symbol = :canonical)
    seq_str = FASTX.sequence(String, read)
    n = length(seq_str)
    n < k && return read, 0
    chars = collect(seq_str)
    # Ambiguity screen: the 2-bit k-mer construction cannot encode non-ACGT bases.
    # Such reads are left for the decode (which tallies un-k-merizable reads).
    for c in chars
        (c == 'A' || c == 'C' || c == 'G' || c == 'T') || return read, 0
    end

    nkmers = n - k + 1
    solid = Vector{Bool}(undef, nkmers)
    @inbounds for i in 1:nkmers
        solid[i] = _lookup_kmer_at(chars, k, i, graph_mode) in solid_kmers
    end

    n_corr = 0
    i = 1
    while i <= nkmers
        if solid[i]
            i += 1
            continue
        end
        # Maximal weak run [a, b].
        a = i
        b = i
        while b + 1 <= nkmers && !solid[b + 1]
            b += 1
        end
        run_len = b - a + 1
        # Interior error: flanked by solid on both sides, run short enough to be a
        # single substitution (the base shared by all k-mers a..b is non-empty).
        flanked = (a - 1 >= 1 && solid[a - 1]) && (b + 1 <= nkmers && solid[b + 1])
        if flanked && run_len <= k
            # Base positions shared by every k-mer in [a, b]: [b, a+k-1] (1-based
            # read coordinates). A single substitution at one of these explains the
            # whole run.
            lo = b
            hi = a + k - 1
            fix_p = 0
            fix_base = ' '
            n_candidates = 0
            for p in lo:hi
                cur = chars[p]
                # k-mers whose window covers base p — the ONLY ones a substitution
                # at p can change. Requiring all of these solid after the edit both
                # fixes the run and guarantees no flanking k-mer is broken.
                lo_j = max(1, p - k + 1)
                hi_j = min(nkmers, p)
                for base in ('A', 'C', 'G', 'T')
                    base == cur && continue
                    chars[p] = base
                    ok = true
                    for j in lo_j:hi_j
                        if !(_lookup_kmer_at(chars, k, j, graph_mode) in solid_kmers)
                            ok = false
                            break
                        end
                    end
                    chars[p] = cur   # restore before testing the next candidate
                    if ok
                        n_candidates += 1
                        n_candidates > 1 && break   # ambiguous: stop early
                        fix_p = p
                        fix_base = base
                    end
                end
                n_candidates > 1 && break
            end
            if n_candidates == 1
                chars[fix_p] = fix_base
                # The corrected k-mers overlapping fix_p are now solid; mark them so
                # the scan does not re-enter this (now resolved) run.
                lo_j = max(1, fix_p - k + 1)
                hi_j = min(nkmers, fix_p)
                for j in lo_j:hi_j
                    solid[j] = true
                end
                n_corr += 1
            end
        end
        i = b + 1
    end

    n_corr == 0 && return read, 0
    corrected_seq = String(chars)
    # Preserve identifier, description, and per-base quality (a substitution does
    # not change the length, so the original quality string still aligns).
    new_record = FASTX.FASTQ.Record(
        FASTX.identifier(read), corrected_seq, FASTX.quality(read))
    return new_record, n_corr
end

"""
    _stage0_cheap_correct(reads, k, solid_kmers) -> (corrected_reads, n_corrections)

Apply `_stage0_correct_read` across a read set, returning a NEW vector of records
(unchanged records are shared, corrected ones replaced) and the total number of
bases corrected. An empty `solid_kmers` (classification could not run) is a no-op:
without a solid reference every k-mer looks weak and there is no trustworthy
neighbor to correct toward, so nothing is changed.
"""
function _stage0_cheap_correct(reads::Vector{<:FASTX.FASTQ.Record}, k::Int,
        solid_kmers::AbstractSet; graph_mode::Symbol = :canonical)
    total = 0
    isempty(solid_kmers) && return collect(reads), 0
    corrected = Vector{FASTX.FASTQ.Record}(undef, length(reads))
    for (idx, read) in enumerate(reads)
        rec, n = _stage0_correct_read(read, k, solid_kmers; graph_mode = graph_mode)
        corrected[idx] = rec
        total += n
    end
    return corrected, total
end

# ------------------------------------------------------------------------------
# Hard-window gating (td-nn6l): decode only reads that touch a "hard" region
# (bubble / repeat) and pass every other read through untouched. This is the
# primary decode-volume reduction on top of skip-solid and Stage 0 cheap
# correction — after Stage 0 clears simple errors, the vast majority of reads
# touch no hard (bubble/repeat) vertex and are skipped.
# ------------------------------------------------------------------------------

"""
    should_decode_read(read, k, hard_vertices; graph_mode) -> Bool

True iff any k-mer of `read` (resolved to its graph label via `_lookup_key`) is a
member of `hard_vertices` (the bubble / repeat-like / weak-k-mer set). Reads that
touch no hard vertex have no region worth decoding and are passed through
untouched. Mirrors `_read_is_all_solid`'s k-mer iteration so the two gates compose
on the same graph. Works on both `:canonical` and `:doublestrand` graphs
(td-nt69). An empty `hard_vertices` means "no gate" ⇒ decode.
"""
function should_decode_read(read::FASTX.FASTQ.Record, k::Int, hard_vertices::AbstractSet;
        graph_mode::Symbol = :canonical)
    isempty(hard_vertices) && return true
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    for (kmer, _) in Kmers.UnambiguousDNAMers{k}(sequence)
        (_lookup_key(kmer, graph_mode) in hard_vertices) && return true
    end
    return false
end

"""
    _hard_vertex_set(graph, k) -> Set

Build the set of "hard" vertices for hard-window gating (td-nn6l/td-bjnt): the
union of
(1) bubble entry/exit/interior vertices (`detect_bubbles_next`), and
(2) high-out-degree (out-degree > 1, repeat-like) vertices.

NARROWED for Stage 0 (td-bjnt): the former clause (3) — "every weak / non-solid
k-mer is hard" — is REMOVED. Under the old clause any read carrying a single
error k-mer (≈78% of reads on err=0.01 data) was classified hard and sent to the
expensive per-read graph Viterbi, which is what times the :scalable corrector out
at 5 kb. Simple single-substitution errors are now fixed CHEAPLY by the Stage 0
k-mer-spectrum pass (`_stage0_cheap_correct`, a linear scan) BEFORE this gate, so
the only reads that must still reach graph Viterbi are those touching GENUINE
ambiguity — a bubble/superbubble (competing balanced alleles) or a repeat-like
high-out-degree vertex. Those are the true ~5–15%. A weak k-mer flanked by solid
neighbors is no longer "hard"; it is an error Stage 0 either already corrected or
left as genuinely ambiguous (in which case it typically also sits in a bubble and
is caught by clause 1). Operates directly on the graph vertices/edges, so it is
mode-agnostic (works on `:canonical` and `:doublestrand` graphs — on a
doublestrand graph a bubble/repeat is flagged on both strands, and reads are
matched in their observed orientation by `should_decode_read`; td-nt69). Only
reached on the :scalable tier (`hard_window=true`), so :exhaustive is unaffected.
"""
function _hard_vertex_set(graph, k::Int)
    labels = collect(MetaGraphsNext.labels(graph))
    T = eltype(labels)
    hard = Set{T}()
    isempty(labels) && return hard

    # (1) Bubble vertices (entry, exit, both alternative interiors) — REAL
    # competing-path ambiguity that a linear cheap-correction cannot resolve.
    for bubble in Mycelia.Rhizomorph.detect_bubbles_next(graph)
        push!(hard, bubble.entry_vertex)
        push!(hard, bubble.exit_vertex)
        for v in bubble.path1
            push!(hard, v)
        end
        for v in bubble.path2
            push!(hard, v)
        end
    end

    # (2) High out-degree (repeat-like) vertices — single O(E) pass over edges.
    # A branch point with >1 successor is a repeat/ambiguity the per-read decode
    # must disambiguate; a cheap linear scan cannot.
    outdeg = Dict{T, Int}()
    for (src, _dst) in MetaGraphsNext.edge_labels(graph)
        outdeg[src] = get(outdeg, src, 0) + 1
    end
    for (v, d) in outdeg
        d > 1 && push!(hard, v)
    end

    return hard
end

# ------------------------------------------------------------------------------
# Stage 3c (td-nn6l): per-hard-region WINDOWED decode — WIRED.
#
# STATUS: ACTIVE under `windowed_decode=true` (the :scalable tier). Instead of
# decoding a hard read end-to-end, `improve_read_likelihood_windowed` extracts the
# k-mer sub-window(s) (<=500 bp) around each hard region (`_hard_window_ranges`)
# and decodes ONLY those windows before splicing each correction back into the
# read. Quality-wrapped observations leave the target endpoint free so the last
# k-mer remains correctable. Windowing bounds the read-length axis to the hard
# neighborhood; graph density remains an independent frontier-cost driver
# (td-2rxh). Whole-read decode (`improve_read_likelihood`) remains the fallback
# when `windowed_decode=false`.
#
# `_hard_window_ranges` (the read-coordinate ranges to decode) is the primitive;
# `improve_read_likelihood_windowed` (below) does the bounded-window decode +
# splice. Both are unit-tested.
# ------------------------------------------------------------------------------

"""
    _hard_window_ranges(read, k, hard_vertices; kwargs...) -> Vector{UnitRange{Int}}

Read-coordinate ranges (1-based, in read bases) covering each hard region of a
read: for every k-mer position whose graph-mode-resolved k-mer is in
`hard_vertices`, take the base span `[pos-pad, pos+k-1+pad]`, clamp to the read,
and merge overlapping spans. With `complete_span=false` (the default), preserve
the pre-indel substitution contract by capping each merged span at its first
`max_window` bases. With `complete_span=true`, cover every merged span using
windows of at most `max_window` bases; consecutive windows overlap by exactly
`k` bases so every later decode has a complete, immutable start k-mer as
context. The windowed splice assigns that overlap to the preceding window.
Empty when the read touches no hard vertex.
"""
function _hard_window_ranges(read::FASTX.FASTQ.Record, k::Int, hard_vertices::AbstractSet;
        pad::Int = 1,
        max_window::Int = 500,
        graph_mode::Symbol = :canonical,
        complete_span::Bool = false)::Vector{UnitRange{Int}}
    max_window >= k || throw(ArgumentError(
        "max_window must be at least k=$k, got $max_window"))
    ranges = UnitRange{Int}[]
    isempty(hard_vertices) && return ranges
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    n = length(sequence)
    n < k && return ranges
    for (km, kpos) in Kmers.UnambiguousDNAMers{k}(sequence)
        if _lookup_key(km, graph_mode) in hard_vertices
            lo = max(1, kpos - pad)
            hi = min(n, kpos + k - 1 + pad)
            push!(ranges, lo:hi)
        end
    end
    isempty(ranges) && return ranges
    # Merge overlapping/adjacent ranges before selecting bounded decode windows.
    sort!(ranges; by = first)
    merged = UnitRange{Int}[]
    cur = ranges[1]
    for r in ranges[2:end]
        if first(r) <= last(cur) + 1
            cur = first(cur):max(last(cur), last(r))
        else
            push!(merged, cur)
            cur = r
        end
    end
    push!(merged, cur)
    if !complete_span
        # This is the historical substitution-only range selection. Keep it exact:
        # changing which tail bases are decoded would violate the byte-identity
        # oracle even though the underlying substitution kernel is unchanged.
        return [
            first(range):min(last(range), first(range) + max_window - 1)
            for range in merged
        ]
    end

    partitioned = UnitRange{Int}[]
    for range in merged
        if length(range) <= max_window
            push!(partitioned, range)
            continue
        end
        owned_capacity = max_window - k
        owned_capacity > 0 || throw(ArgumentError(
            "max_window=$max_window must exceed k=$k to cover a hard span " *
            "longer than one window with anchored overlap"))

        n_windows = 1 + cld(length(range) - max_window, owned_capacity)
        base_owned, longer_owned = divrem(length(range), n_windows)
        ownership_lengths = fill(base_owned, n_windows)
        @inbounds for index in 1:longer_owned
            ownership_lengths[index] += 1
        end
        # A later window can own at most `max_window-k` new bases because its
        # first k decode bases belong to the preceding window. Shift any balanced
        # excess into the larger-capacity first window.
        @inbounds for index in 2:n_windows
            excess = max(0, ownership_lengths[index] - owned_capacity)
            ownership_lengths[index] -= excess
            ownership_lengths[1] += excess
        end
        ownership_lengths[1] <= max_window || error(
            "internal anchored-window partition exceeded first-window capacity")
        ownership_lengths[1] >= k || error(
            "internal anchored-window partition produced a short first window")

        ownership_start = first(range)
        @inbounds for index in eachindex(ownership_lengths)
            ownership_stop = ownership_start + ownership_lengths[index] - 1
            # Later decodes prepend the preceding window's final k bases as an
            # immutable start anchor. Splicing trims exactly this prefix.
            decode_start = index == 1 ? ownership_start : ownership_start - k
            push!(partitioned, decode_start:ownership_stop)
            ownership_start = ownership_stop + 1
        end
    end
    return partitioned
end

function _indel_frontier_metrics_dict(
        metrics::Any,
)::Dict{Symbol, Any}
    return Dict{Symbol, Any}(
        :anchored => metrics.anchored,
        :window_length => metrics.window_length,
        :vertex_count => metrics.vertex_count,
        :edge_count => metrics.edge_count,
        :branch_vertices => metrics.branch_vertices,
        :join_vertices => metrics.join_vertices,
        :branch_fraction => metrics.branch_fraction,
        :join_fraction => metrics.join_fraction,
        :max_out_degree => metrics.max_out_degree,
        :frontier_area => metrics.frontier_area,
        :edge_expansions => metrics.edge_expansions,
        :peak_frontier => metrics.peak_frontier,
        :completed_columns => metrics.completed_columns,
        :frontier_work => metrics.frontier_work,
        :reason => metrics.reason,
    )
end

function _indel_frontier_admitted(
        metrics::Any,
        work_limit::Int,
)::Bool
    return metrics.anchored && metrics.reason == :complete &&
           metrics.completed_columns == metrics.window_length &&
           metrics.frontier_work <= work_limit
end

function _indel_candidate_windows(
        reads::Vector{<:FASTX.FASTQ.Record},
        k::Int,
        hard_vertices::AbstractSet,
        skip_flags::Vector{Bool},
        graph_mode::Symbol,
)::Vector{Tuple{Int, UnitRange{Int}}}
    candidates = Tuple{Int, UnitRange{Int}}[]
    @inbounds for i in eachindex(reads)
        skip_flags[i] && continue
        windows = _hard_window_ranges(
            reads[i], k, hard_vertices;
            pad = k,
            max_window = _INDEL_FRONTIER_MAX_WINDOW,
            graph_mode = graph_mode,
            complete_span = true)
        for window in windows
            length(window) >= k && push!(candidates, (i, window))
        end
    end
    return candidates
end

function _indel_probe_config(
        params::IndelDecodeParams,
        graph_mode::Symbol,
)::Any
    return Mycelia.ViterbiCorrectionConfig(
        alphabet = :DNA,
        strand_mode = graph_mode,
        max_steps = _INDEL_FRONTIER_MAX_WINDOW - 1,
        beam_width = _AUTO_BEAM_BOUNDED_WIDTH,
        error_rate = params.base_error_rate,
        indel_moves = true,
        insertion_fraction = params.insertion_fraction,
        deletion_fraction = params.deletion_fraction,
        insertion_extend_probability = params.insertion_extend_probability,
        deletion_extend_probability = params.deletion_extend_probability,
        deletion_max_run = params.deletion_max_run,
        max_insertion_run = params.max_insertion_run,
        band_width = params.band_width,
    )
end

function _probe_indel_window(
        graph::MetaGraphsNext.MetaGraph,
        read::FASTX.FASTQ.Record,
        window::UnitRange{Int},
        k::Int,
        config::Any,
        graph_mode::Symbol,
        work_limit::Int,
        graph_summary::NamedTuple,
)::Any
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    window_sequence = sequence[window]
    observations = [kmer for (kmer, _) in Kmers.UnambiguousDNAMers{k}(window_sequence)]
    return Mycelia._probe_indel_frontier(
        graph, observations, :DNA;
        config = config,
        strand_mode = graph_mode,
        work_limit = work_limit,
        graph_summary = graph_summary,
    )
end

function _probe_indel_window_metric(
        graph::MetaGraphsNext.MetaGraph,
        read::FASTX.FASTQ.Record,
        window::UnitRange{Int},
        k::Int,
        config::Any,
        graph_mode::Symbol,
        work_limit::Int,
        graph_summary::NamedTuple,
)::Tuple{Dict{Symbol, Any}, Bool}
    try
        metrics = _probe_indel_window(
            graph, read, window, k, config, graph_mode, work_limit,
            graph_summary)
        return (
            _indel_frontier_metrics_dict(metrics),
            _indel_frontier_admitted(metrics, work_limit),
        )
    catch exception
        exception isa BioSequences.EncodeError || rethrow()
        return (
            Dict{Symbol, Any}(
                :anchored => false,
                :reason => :encode_error,
            ),
            false,
        )
    end
end

function _evaluate_indel_frontier_schedule(
        reads::Vector{<:FASTX.FASTQ.Record},
        raw_graph::MetaGraphsNext.MetaGraph,
        k::Int,
        hard_vertices::AbstractSet,
        skip_flags::Vector{Bool},
        params::IndelDecodeParams,
        graph_mode::Symbol;
        work_limit::Int = _DEFAULT_INDEL_FRONTIER_WORK_LIMIT,
        prior_soft_weights::Union{
            Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        memory_limit::Union{Nothing, Int} = nothing,
)::NamedTuple
    return Mycelia.Rhizomorph._with_soft_edge_weight_scope(
        raw_graph, prior_soft_weights) do
        _evaluate_indel_frontier_schedule_impl(
            reads,
            raw_graph,
            k,
            hard_vertices,
            skip_flags,
            params,
            graph_mode;
            work_limit = work_limit,
            prior_soft_weights = prior_soft_weights,
            memory_limit = memory_limit,
        )
    end
end

function _indel_graph_memory_bytes(
        raw_graph::MetaGraphsNext.MetaGraph,
)::Int
    # Measure the live graph rather than inferring it from vertex count. This
    # includes actual edge topology and evidence payloads, so a small branching
    # graph cannot masquerade as cheap merely because `nv` is small.
    return Int(Base.summarysize(raw_graph))
end

function _indel_additional_graph_bytes(
        graph_bytes::Int,
        variants::Int,
)::Int
    additional_variants = variants - 1
    additional_variants <= 0 && return 0
    return graph_bytes > fld(typemax(Int), additional_variants) ?
           typemax(Int) : graph_bytes * additional_variants
end

function _indel_graph_variants_fit_memory(
        raw_graph::MetaGraphsNext.MetaGraph,
        variants::Int,
        memory_limit::Union{Nothing, Int},
)::Bool
    variants > 0 || throw(ArgumentError("variants must be positive, got $variants"))
    memory_limit === nothing && return true
    memory_limit <= 0 && return false
    graph_bytes = _indel_graph_memory_bytes(raw_graph)
    current_live_bytes = Int(Base.gc_live_bytes())
    current_live_bytes > memory_limit && return false
    additional_bytes = _indel_additional_graph_bytes(graph_bytes, variants)
    return additional_bytes <= memory_limit - current_live_bytes
end

function _indel_graph_memory_telemetry(
        raw_graph::MetaGraphsNext.MetaGraph,
        variants::Int,
        memory_limit::Union{Nothing, Int},
        status::String,
)::Dict{String, Any}
    graph_bytes = _indel_graph_memory_bytes(raw_graph)
    current_live_bytes = Int(Base.gc_live_bytes())
    additional_bytes = _indel_additional_graph_bytes(graph_bytes, variants)
    projected_peak = current_live_bytes > typemax(Int) - additional_bytes ?
                     typemax(Int) : current_live_bytes + additional_bytes
    return Dict{String, Any}(
        "graph_cleanup_status" => status,
        "measured_graph_bytes" => graph_bytes,
        "current_live_bytes" => current_live_bytes,
        "projected_graph_variants" => variants,
        "projected_peak_bytes" => projected_peak,
        "memory_limit" => memory_limit,
    )
end

function _try_build_correction_weighted_graph(
        graph::MetaGraphsNext.MetaGraph,
)::Tuple{Any, Bool}
    try
        return (build_correction_weighted_graph(graph), false)
    catch exception
        exception isa OutOfMemoryError || rethrow()
        GC.gc(true)
        return (nothing, true)
    end
end

function _evaluate_indel_frontier_schedule_impl(
        reads::Vector{<:FASTX.FASTQ.Record},
        raw_graph::MetaGraphsNext.MetaGraph,
        k::Int,
        hard_vertices::AbstractSet,
        skip_flags::Vector{Bool},
        params::IndelDecodeParams,
        graph_mode::Symbol;
        work_limit::Int,
        prior_soft_weights::Union{
            Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator},
        memory_limit::Union{Nothing, Int},
)::NamedTuple
    candidates = _indel_candidate_windows(
        reads, k, hard_vertices, skip_flags, graph_mode)
    requested = length(candidates)
    window_sources = Dict{Int, Dict{UnitRange{Int}, Symbol}}()
    for (read_index, window) in candidates
        sources = get!(window_sources, read_index, Dict{UnitRange{Int}, Symbol}())
        sources[window] = :substitution
    end
    if isempty(candidates)
        return (
            graph = raw_graph,
            raw_weighted_graph = nothing,
            cleaned_graph = nothing,
            cleaned_weighted_graph = nothing,
            window_sources,
            admitted = false,
            requested,
            admitted_windows = 0,
            rejected_windows = 0,
            graph_source = :raw,
            decision_reason = :no_candidate_windows,
            raw_evaluated = 0,
            cleaned_evaluated = 0,
            raw_metrics = Dict{Symbol, Any}[],
            cleaned_metrics = Dict{Symbol, Any}[],
            cleanup = Dict{String, Any}(),
        )
    end

    config = _indel_probe_config(params, graph_mode)
    raw_summary = Mycelia._indel_frontier_graph_summary(raw_graph)
    raw_metrics = Dict{Symbol, Any}[]
    raw_evaluated = 0
    raw_admitted = fill(false, requested)
    raw_anchored = fill(false, requested)
    for (candidate_index, (read_index, window)) in enumerate(candidates)
        metric, admitted = _probe_indel_window_metric(
            raw_graph,
            reads[read_index],
            window,
            k,
            config,
            graph_mode,
            work_limit,
            raw_summary,
        )
        metric[:read_index] = read_index
        metric[:window_start] = first(window)
        metric[:window_stop] = last(window)
        metric[:admitted] = admitted
        raw_evaluated += 1
        if length(raw_metrics) < _INDEL_FRONTIER_TELEMETRY_SAMPLE_LIMIT
            push!(raw_metrics, metric)
        end
        raw_admitted[candidate_index] = admitted
        raw_anchored[candidate_index] = get(metric, :anchored, false)
        admitted && (window_sources[read_index][window] = :raw)
    end

    if all(raw_admitted)
        if !_indel_graph_variants_fit_memory(raw_graph, 2, memory_limit)
            for (read_index, window) in candidates
                window_sources[read_index][window] = :substitution
            end
            cleanup = _indel_graph_memory_telemetry(
                raw_graph, 2, memory_limit, "skipped_weighted_memory_limit")
            return (
                graph = raw_graph,
                raw_weighted_graph = nothing,
                cleaned_graph = nothing,
                cleaned_weighted_graph = nothing,
                window_sources,
                admitted = false,
                requested,
                admitted_windows = 0,
                rejected_windows = requested,
                graph_source = :substitution,
                decision_reason = :weighted_graph_memory_limit,
                raw_evaluated,
                cleaned_evaluated = 0,
                raw_metrics,
                cleaned_metrics = Dict{Symbol, Any}[],
                cleanup,
            )
        end
        raw_weighted_graph, weighted_out_of_memory =
            _try_build_correction_weighted_graph(raw_graph)
        if weighted_out_of_memory
            for (read_index, window) in candidates
                window_sources[read_index][window] = :substitution
            end
            return (
                graph = raw_graph,
                raw_weighted_graph = nothing,
                cleaned_graph = nothing,
                cleaned_weighted_graph = nothing,
                window_sources,
                admitted = false,
                requested,
                admitted_windows = 0,
                rejected_windows = requested,
                graph_source = :substitution,
                decision_reason = :weighted_graph_out_of_memory,
                raw_evaluated,
                cleaned_evaluated = 0,
                raw_metrics,
                cleaned_metrics = Dict{Symbol, Any}[],
                cleanup = Dict{String, Any}(
                    "graph_cleanup_status" => "weighted_out_of_memory",
                ),
            )
        end
        return (
            graph = raw_graph,
            raw_weighted_graph,
            cleaned_graph = nothing,
            cleaned_weighted_graph = nothing,
            window_sources,
            admitted = true,
            requested,
            admitted_windows = requested,
            rejected_windows = 0,
            graph_source = :raw,
            decision_reason = :raw_frontier_affordable,
            raw_evaluated,
            cleaned_evaluated = 0,
            raw_metrics,
            cleaned_metrics = Dict{Symbol, Any}[],
            cleanup = Dict{String, Any}(),
        )
    end

    # Cleaning is a rescue attempt on a private copy only. Probe only windows that
    # exceeded the raw frontier budget; raw-affordable windows keep the raw graph,
    # while rejected windows retain the substitution-only kernel/config. Weighted
    # graph copies are built lazily only for sources that will actually decode.
    cleaned_graph = nothing
    cleanup = Dict{String, Any}()
    cleaned_metrics = Dict{Symbol, Any}[]
    cleaned_evaluated = 0
    cleaning_lost_anchor = false
    cleaned_admitted = fill(false, requested)
    cleaned_weighted_graph = nothing
    cleaning_status = :not_attempted
    if _indel_graph_variants_fit_memory(raw_graph, 4, memory_limit)
        try
            cleaned_graph = deepcopy(raw_graph)
            cleanup = Mycelia.Rhizomorph.clean_corrector_graph!(cleaned_graph; k = k)
            # Scope cleaned registrations separately inside the raw-graph scope.
            # If cleaning/probing/building throws (including OOM), the nested
            # finally releases every cleaned edge before recovery allocates or GC
            # runs; the outer scope still preserves the raw registry exactly.
            Mycelia.Rhizomorph._with_soft_edge_weight_scope(
                cleaned_graph, prior_soft_weights) do
                cleaned_summary =
                    Mycelia._indel_frontier_graph_summary(cleaned_graph)
                for (candidate_index, (read_index, window)) in enumerate(candidates)
                    raw_admitted[candidate_index] && continue
                    metric, admitted = _probe_indel_window_metric(
                        cleaned_graph,
                        reads[read_index],
                        window,
                        k,
                        config,
                        graph_mode,
                        work_limit,
                        cleaned_summary,
                    )
                    metric[:read_index] = read_index
                    metric[:window_start] = first(window)
                    metric[:window_stop] = last(window)
                    metric[:admitted] = admitted
                    cleaned_evaluated += 1
                    cleaning_lost_anchor |= raw_anchored[candidate_index] &&
                                            !get(metric, :anchored, false)
                    if length(cleaned_metrics) <
                       _INDEL_FRONTIER_TELEMETRY_SAMPLE_LIMIT
                        push!(cleaned_metrics, metric)
                    end
                    cleaned_admitted[candidate_index] = admitted
                    admitted && (window_sources[read_index][window] = :cleaned)
                end
                if any(cleaned_admitted)
                    cleaned_weighted_graph =
                        build_correction_weighted_graph(cleaned_graph)
                end
            end
            cleaning_status = :complete
        catch exception
            exception isa InterruptException && rethrow()
            out_of_memory = exception isa OutOfMemoryError
            cleaned_graph = nothing
            cleaned_weighted_graph = nothing
            empty!(cleaned_metrics)
            cleaned_evaluated = 0
            fill!(cleaned_admitted, false)
            for (candidate_index, (read_index, window)) in enumerate(candidates)
                raw_admitted[candidate_index] && continue
                window_sources[read_index][window] = :substitution
            end
            if out_of_memory
                # Avoid formatting the exception while memory is exhausted. All
                # cleaned registry references are already gone via the nested scope.
                empty!(cleanup)
                GC.gc(true)
                cleanup = Dict{String, Any}(
                    "graph_cleanup_status" => "out_of_memory",
                    "graph_cleanup_error_type" => "OutOfMemoryError",
                )
                cleaning_status = :out_of_memory
            else
                cleanup = Dict{String, Any}(
                    "graph_cleanup_status" => "error",
                    "graph_cleanup_error_type" => string(typeof(exception)),
                    "graph_cleanup_error" => sprint(showerror, exception),
                )
                cleaning_status = :error
            end
        end
    else
        cleanup = _indel_graph_memory_telemetry(
            raw_graph, 4, memory_limit, "skipped_memory_limit")
        cleaning_status = :memory_limit
    end

    if any(raw_admitted) &&
       !_indel_graph_variants_fit_memory(raw_graph, 2, memory_limit)
        for (candidate_index, (read_index, window)) in enumerate(candidates)
            raw_admitted[candidate_index] || continue
            raw_admitted[candidate_index] = false
            window_sources[read_index][window] = :substitution
        end
        cleanup = _indel_graph_memory_telemetry(
            raw_graph, 2, memory_limit, "skipped_weighted_memory_limit")
        cleaning_status = :weighted_memory_limit
    end

    admitted_mask = raw_admitted .| cleaned_admitted
    admitted_windows = count(identity, admitted_mask)
    rejected_windows = requested - admitted_windows
    used_sources = Set(
        source for sources in values(window_sources) for source in values(sources))
    graph_source = length(used_sources) == 1 ? only(used_sources) : :mixed

    if admitted_windows > 0
        reason = if rejected_windows > 0
            :sparse_frontier_affordable
        elseif graph_source == :mixed
            :mixed_frontier_affordable
        else
            :cleaned_frontier_affordable
        end
        needs_raw_graph = :raw in used_sources || :substitution in used_sources
        needs_cleaned_graph = :cleaned in used_sources
        raw_weighted_graph = nothing
        weighted_out_of_memory = false
        if needs_raw_graph
            raw_weighted_graph, weighted_out_of_memory =
                _try_build_correction_weighted_graph(raw_graph)
        end
        if weighted_out_of_memory
            for (read_index, window) in candidates
                window_sources[read_index][window] = :substitution
            end
            return (
                graph = raw_graph,
                raw_weighted_graph = nothing,
                cleaned_graph = nothing,
                cleaned_weighted_graph = nothing,
                window_sources,
                admitted = false,
                requested,
                admitted_windows = 0,
                rejected_windows = requested,
                graph_source = :substitution,
                decision_reason = :weighted_graph_out_of_memory,
                raw_evaluated,
                cleaned_evaluated,
                raw_metrics,
                cleaned_metrics,
                cleanup = Dict{String, Any}(
                    "graph_cleanup_status" => "weighted_out_of_memory",
                ),
            )
        end
        return (
            graph = raw_graph,
            raw_weighted_graph,
            cleaned_graph = needs_cleaned_graph ? cleaned_graph : nothing,
            cleaned_weighted_graph,
            window_sources,
            admitted = true,
            requested,
            admitted_windows,
            rejected_windows,
            graph_source,
            decision_reason = reason,
            raw_evaluated,
            cleaned_evaluated,
            raw_metrics,
            cleaned_metrics,
            cleanup,
        )
    end

    reason = if cleaning_status == :error
        :cleaning_error
    elseif cleaning_status == :out_of_memory
        :cleaning_out_of_memory
    elseif cleaning_status == :weighted_memory_limit
        :weighted_graph_memory_limit
    elseif cleaning_status == :memory_limit
        :cleaning_memory_limit
    elseif cleaning_lost_anchor
        :cleaning_lost_window_anchor
    else
        :frontier_budget_exceeded
    end
    return (
        graph = raw_graph,
        raw_weighted_graph = nothing,
        cleaned_graph = nothing,
        cleaned_weighted_graph = nothing,
        window_sources,
        admitted = false,
        requested,
        admitted_windows,
        rejected_windows,
        graph_source,
        decision_reason = reason,
        raw_evaluated,
        cleaned_evaluated,
        raw_metrics,
        cleaned_metrics,
        cleanup,
    )
end

"""
    improve_read_likelihood_windowed(read, graph, k, hard_vertices; ...) -> (FASTX.FASTQ.Record, Bool)

Stage 3c WINDOWED decode (td-nn6l). Instead of decoding a hard read end-to-end,
decode ONLY the hard sub-window(s) — the k-mer span(s) around the read's hard
vertices from `_hard_window_ranges` (each `<= max_window` bases) — then splice
each corrected window back into the read. Each window is decoded as its own
`improve_read_likelihood` sub-read. Quality-wrapped observations leave the target
vertex free so the decoder can correct the window's last k-mer; the window bounds
the observation-length axis of the decode, not the full graph frontier. This makes
the read-length contribution scale with the hard-window size rather than the full
read length, but graph density remains a separate cost driver (#375, td-2rxh).

With `indel_params === nothing`, the window Viterbi is length-preserving (the
corrected ML path has one vertex per window k-mer), so each accepted window is
spliced in as a same-length substitution. A length mismatch on that path is
dropped defensively and counted in `divergent_windows`. With non-`nothing`
`indel_params`, the pair-HMM may emit insertions or deletions; those intentional
length changes are accepted and spliced by rebuilding the read from original
window coordinates. Windows shorter than `k` (which cannot be k-merized) are
skipped.

A read with no hard window (empty ranges) returns unchanged (`false`) — this only
occurs for reads the caller's gate already classified as skip, so the whole-read
fallback (`windowed_decode=false`) is the caller's, not this function's.

Each window requests a k-mer-width flank (`pad === nothing` ⇒ `k` bases), keeping
hard vertices away from a boundary when read ends and `max_window` permit. For an
indel schedule, overlong spans use balanced windows whose consecutive decode ranges
overlap by one anchored k-mer; the overlap is trimmed from the later correction
before splice ownership is assigned. Profile-disabled substitution calls retain the
historical first-window cap exactly. Padding is contextual rather than a solidity
guarantee. The quality-aware decoder deliberately leaves the target endpoint free
so the last k-mer remains correctable. The caller decides WHEN to window.

Returns `(record, improved)` where `improved` is `true` iff `>= 1` window was
corrected and spliced. The lower-level `improve_read_likelihood_windowed_detail`
additionally returns the per-read decoded-window count and divergent-window count
for correctness/telemetry checks.
"""
function improve_read_likelihood_windowed(read::FASTX.FASTQ.Record, graph, k::Int,
        hard_vertices::AbstractSet;
        graph_mode::Symbol = :canonical,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        weighted_graph = nothing,
        cleaned_graph = nothing,
        cleaned_weighted_graph = nothing,
        indel_window_sources::Union{
            Nothing, Dict{UnitRange{Int}, Symbol}} = nothing,
        diagnostics = nothing,  # ::Union{Nothing, CorrectorDiagnostics}; struct defined below
        pad::Union{Int, Nothing} = nothing,
        calibrated_gap_threshold::Union{Float64, Nothing} = nothing,
        max_window::Int = 500,
        indel_params::Union{Nothing, IndelDecodeParams} = nothing)::Tuple{
        FASTX.FASTQ.Record, Bool}
    record, improved,
    _decoded,
    _divergent = improve_read_likelihood_windowed_detail(
        read, graph, k, hard_vertices;
        graph_mode = graph_mode, beam_width = beam_width, soft_weights = soft_weights,
        weighted_graph = weighted_graph, cleaned_graph = cleaned_graph,
        cleaned_weighted_graph = cleaned_weighted_graph,
        indel_window_sources = indel_window_sources, diagnostics = diagnostics,
        pad = pad, calibrated_gap_threshold = calibrated_gap_threshold,
        max_window = max_window, indel_params = indel_params)
    return record, improved
end

"""
    improve_read_likelihood_windowed_detail(read, graph, k, hard_vertices; ...)
        -> (FASTX.FASTQ.Record, Bool, Int, Int)

Windowed-decode core (see `improve_read_likelihood_windowed`). Returns
`(record, improved, decoded_windows, divergent_windows)` where `decoded_windows`
is the number of hard windows that reached a Viterbi decode and `divergent_windows`
is the number dropped for a length mismatch.

`indel_params` (td-jt7r): `nothing` ⇒ SUBSTITUTION decode. The window ML path has
one vertex per window k-mer, so the decode is length-preserving and each accepted
window is spliced back IN PLACE; a window whose corrected length differs is DROPPED
(defensive `divergent_windows` guard; expected 0). Non-`nothing` ⇒ the INDEL
pair-HMM decode, which emits I/D moves and thus CHANGES a window's length by design
— the length delta is the correction, not a defect. Those windows are accepted at
their new length and the read is reassembled by ORIGINAL-coordinate segment
rebuild, so an earlier window's length change cannot shift a later window's range.
The substitution path (`indel_params === nothing`) uses the unchanged
substitution-only kernel and parameters. It may correct bases inside its bounded
window; the identifier, qualities paired to unchanged bases, and all untouched
spans remain byte-identical. Exposed for the windowed-decode correctness test.
"""
const _VALID_INDEL_WINDOW_SOURCES = (:raw, :cleaned, :substitution)

function _validate_indel_window_sources(
        sources::Union{Nothing, Dict{UnitRange{Int}, Symbol}},
        cleaned_graph::Any,
        cleaned_weighted_graph::Any,
)::Nothing
    sources === nothing && return nothing
    for (window, source) in sources
        source in _VALID_INDEL_WINDOW_SOURCES || throw(ArgumentError(
            "invalid indel window source :$source for $window; expected one of " *
            join(string.(_VALID_INDEL_WINDOW_SOURCES), ", ")))
        if source == :cleaned &&
           (cleaned_graph === nothing || cleaned_weighted_graph === nothing)
            throw(ArgumentError(
                "indel window source :cleaned for $window requires both a " *
                "cleaned graph and cleaned weighted graph"))
        end
    end
    return nothing
end

function _trim_indel_window_overlap(
        original_sequence::Vector{Char},
        decoded_sequence::Vector{Char},
        decoded_quality::Vector{Char},
        window_start::Int,
        overlap_prefix::Int,
)::Union{Nothing, Tuple{Vector{Char}, Vector{Char}}}
    overlap_prefix == 0 && return (decoded_sequence, decoded_quality)
    overlap_prefix > 0 || throw(ArgumentError(
        "overlap_prefix must be nonnegative"))
    anchor_stop = window_start + overlap_prefix - 1
    1 <= window_start <= anchor_stop <= length(original_sequence) ||
        throw(ArgumentError("overlap anchor lies outside the original sequence"))
    if length(decoded_sequence) < overlap_prefix ||
       length(decoded_quality) < overlap_prefix
        return nothing
    end
    expected_prefix = original_sequence[window_start:anchor_stop]
    decoded_sequence[1:overlap_prefix] == expected_prefix || return nothing
    return (
        decoded_sequence[(overlap_prefix + 1):end],
        decoded_quality[(overlap_prefix + 1):end],
    )
end

function improve_read_likelihood_windowed_detail(read::FASTX.FASTQ.Record, graph, k::Int,
        hard_vertices::AbstractSet;
        graph_mode::Symbol = :canonical,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        weighted_graph = nothing,
        cleaned_graph = nothing,
        cleaned_weighted_graph = nothing,
        indel_window_sources::Union{
            Nothing, Dict{UnitRange{Int}, Symbol}} = nothing,
        diagnostics = nothing,  # ::Union{Nothing, CorrectorDiagnostics}; struct defined below
        pad::Union{Int, Nothing} = nothing,
        calibrated_gap_threshold::Union{Float64, Nothing} = nothing,
        max_window::Int = 500,
        indel_params::Union{Nothing, IndelDecodeParams} = nothing)::Tuple{
        FASTX.FASTQ.Record, Bool, Int, Int}
    _validate_indel_window_sources(
        indel_window_sources, cleaned_graph, cleaned_weighted_graph)
    # Anchor each window with a full solid k-mer flank by default (`pad === nothing`
    # ⇒ `k`): request a k-mer-width flank around the hard span when read boundaries
    # and `max_window` allow it. This is contextual padding, not a guarantee that a
    # boundary k-mer is solid. Quality-aware decoding leaves the target endpoint free
    # so the last k-mer remains correctable.
    effective_pad = pad === nothing ? k : pad
    windows = _hard_window_ranges(
        read, k, hard_vertices;
        pad = effective_pad,
        max_window = max_window,
        graph_mode = graph_mode,
        complete_span = indel_params !== nothing || indel_window_sources !== nothing)
    isempty(windows) && return read, false, 0, 0

    seq_chars = collect(FASTX.sequence(String, read))
    qual_chars = collect(FASTX.quality(read))   # ASCII quality string, one char/base
    id = FASTX.identifier(read)
    decoded_windows = 0
    divergent_windows = 0
    # Accepted (window range, corrected seq, corrected qual) triples, spliced AFTER
    # the loop. Deferring the splice lets the indel path (length-changing) reassemble
    # from original coordinates; the substitution path splices these in place.
    accepted = Tuple{UnitRange{Int}, String, String}[]

    previous_window_stop = 0
    for w in windows
        lo = first(w)
        hi = last(w)
        overlap_prefix = max(0, previous_window_stop - lo + 1)
        previous_window_stop = max(previous_window_stop, hi)
        win_len = hi - lo + 1
        # A window shorter than k cannot be k-merized (no observations) — skip.
        win_len < k && continue
        sub_seq = String(seq_chars[lo:hi])
        sub_qual = String(qual_chars[lo:hi])
        sub_read = FASTX.FASTQ.Record("$(id)_win$(lo)_$(hi)", sub_seq, sub_qual)
        decoded_windows += 1
        # Decode ONLY this window. Quality-wrapped observations keep the target
        # endpoint free; the window bounds observation length, while graph density
        # remains an independent frontier-cost driver (td-2rxh). The pair-HMM runs
        # on the short window rather than the full read.
        window_source = indel_window_sources === nothing ?
                        (indel_params === nothing ? :substitution : :raw) :
                        get(indel_window_sources, w, :substitution)
        window_graph = window_source == :cleaned ? cleaned_graph : graph
        window_weighted_graph = window_source == :cleaned ?
                                cleaned_weighted_graph : weighted_graph
        window_indel_params = window_source == :substitution ? nothing : indel_params
        if window_graph === nothing
            diagnostics === nothing ||
                Threads.atomic_add!(diagnostics.structural_errors, 1)
            continue
        end
        decoded_sub,
        improved = improve_read_likelihood(
            sub_read, window_graph, k; graph_mode = graph_mode,
            beam_width = beam_width, soft_weights = soft_weights,
            weighted_graph = window_weighted_graph, diagnostics = diagnostics,
            calibrated_gap_threshold = calibrated_gap_threshold,
            indel_params = window_indel_params)
        improved || continue
        dseq = FASTX.sequence(String, decoded_sub)
        dqual = String(FASTX.quality(decoded_sub))
        dseq_chars = collect(dseq)
        dqual_chars = collect(dqual)
        if overlap_prefix > 0
            # Complete-span indel windows overlap by one anchored k-mer. The
            # preceding window owns those original coordinates; trim the immutable
            # context before splicing this correction. Fail closed if a lower layer
            # ever rewrites or deletes the anchor instead of assuming coordinates.
            trimmed = _trim_indel_window_overlap(
                seq_chars,
                dseq_chars,
                dqual_chars,
                lo,
                overlap_prefix,
            )
            if trimmed === nothing
                diagnostics === nothing ||
                    Threads.atomic_add!(diagnostics.trace_contract_errors, 1)
                continue
            end
            dseq_chars, dqual_chars = trimmed
        end
        owned_window = (lo + overlap_prefix):hi
        if window_indel_params === nothing
            # Length-preserving splice: a divergent length would corrupt read
            # coordinates for the in-place splice below, so drop it (expected 0).
            if length(dseq) == win_len && length(dqual) == win_len
                push!(accepted,
                    (owned_window, String(dseq_chars), String(dqual_chars)))
            else
                divergent_windows += 1
                diagnostics === nothing ||
                    Threads.atomic_add!(diagnostics.window_divergences, 1)
            end
        else
            # The pair-HMM traceback already rebuilt quality coordinates in
            # `try_viterbi_path_improvement`. Reject a malformed lower-layer record
            # rather than padding/truncating it and silently shifting qualities.
            if length(dqual) != length(dseq)
                diagnostics === nothing ||
                    Threads.atomic_add!(diagnostics.structural_errors, 1)
                continue
            end
            push!(accepted,
                (owned_window, String(dseq_chars), String(dqual_chars)))
        end
    end

    isempty(accepted) &&
        return read, false, decoded_windows, divergent_windows

    if indel_params === nothing
        # Substitution: length-preserving in-place overwrite (byte-identical to the
        # pre-indel path — oracle preservation).
        for (w, dseq, dqual) in accepted
            dseq_chars = collect(dseq)
            dqual_chars = collect(dqual)
            @inbounds for (j, p) in enumerate(first(w):last(w))
                seq_chars[p] = dseq_chars[j]
                qual_chars[p] = dqual_chars[j]
            end
        end
        corrected = FASTX.FASTQ.Record(id, String(seq_chars), String(qual_chars))
        return corrected, true, decoded_windows, divergent_windows
    end

    # Indel: rebuild the read from ORIGINAL-coordinate segments (see
    # `_splice_indel_windows`) so a window's length change cannot shift a later
    # window's range.
    corrected = _splice_indel_windows(id, seq_chars, qual_chars, accepted)
    return corrected, true, decoded_windows, divergent_windows
end

"""
    _splice_indel_windows(id, seq_chars, qual_chars, accepted) -> FASTX.FASTQ.Record

Rebuild a read from ORIGINAL-coordinate segments for the indel (length-changing)
windowed decode. `seq_chars` / `qual_chars` are the ORIGINAL read's characters;
`accepted` is a vector of `(owned_range, corrected_seq, corrected_qual)` triples in
window order. The decode ranges may overlap by an anchored k-mer, but callers trim
that contextual prefix before constructing these sorted, non-overlapping ownership
ranges. Each original ownership span is replaced by its (possibly different-length)
corrected sequence; the solid gaps BETWEEN ranges and the trailing span are copied
verbatim from the original. Because segments are cut at ORIGINAL coordinates (not
a running offset), an earlier window's length change cannot shift a later range.
Exposed for a deterministic unit test of the coordinate math (adjacent / gap /
trailing / length-changing windows) independent of the decoder.
"""
function _splice_indel_windows(id::AbstractString,
        seq_chars::Vector{Char}, qual_chars::Vector{Char},
        accepted::AbstractVector{<:Tuple{
            UnitRange{Int}, <:AbstractString, <:AbstractString}})
    out_seq = IOBuffer()
    out_qual = IOBuffer()
    prev = 0
    for (w, dseq, dqual) in accepted
        lo = first(w)
        hi = last(w)
        print(out_seq, String(seq_chars[(prev + 1):(lo - 1)]))   # solid gap before window
        print(out_qual, String(qual_chars[(prev + 1):(lo - 1)]))
        print(out_seq, dseq)                                     # corrected window
        print(out_qual, dqual)
        prev = hi
    end
    print(out_seq, String(seq_chars[(prev + 1):end]))            # trailing solid span
    print(out_qual, String(qual_chars[(prev + 1):end]))
    return FASTX.FASTQ.Record(id, String(take!(out_seq)), String(take!(out_qual)))
end

"""
    _windowed_decode_read_is_long(read, k) -> Bool

Gate for WHEN to use windowed decode (td-nn6l Stage 3c). Returns `true` only when
the read's observation count exceeds `_AUTO_BEAM_EXACT_THRESHOLD`, bounding the
read-length axis by decoding `<=max_window` sub-windows. Exactness still depends on
graph density and any explicit beam override: a short window on a graph with more
than `_AUTO_BEAM_BOUNDED_WIDTH` vertices remains beam-bounded. Reads below
the length threshold keep whole-read context unless another caller policy (such as
staged indel decoding, td-2rxh) explicitly forces windowing.
"""
function _windowed_decode_read_is_long(read::FASTX.FASTQ.Record, k::Int)::Bool
    n_obs = length(FASTX.sequence(read)) - k + 1
    return n_obs > _AUTO_BEAM_EXACT_THRESHOLD
end

# ------------------------------------------------------------------------------
# Corrector diagnostics + size-aware auto-beam (review: robustness blockers)
# ------------------------------------------------------------------------------

"""
    CorrectorDiagnostics()

Thread-safe tally of decode outcomes the per-read corrector would otherwise
swallow SILENTLY. When `try_viterbi_path_improvement` fails to decode a read it
returns `nothing` and the read passes through uncorrected; without this tally a
SYSTEMATICALLY broken corrector — empty graph, alphabet-inference miss, config
error, or a read that cannot be k-merized — reports `0/N improvements` as if it
had actually run and found nothing to fix. Counters are `Threads.Atomic` so the
parallel (`@threads`) read loop increments them race-free.

- `structural_errors`  : `ArgumentError` from the decoder (empty graph,
  alphabet-inference miss, config error). Means "the decoder could not run",
  NOT "the decoder ran and found no gain".
- `unkmerizable_reads` : `BioSequences.EncodeError` — a read carrying an
  ambiguous base (e.g. `N`) that the 2-bit k-mer iterator cannot encode. The
  read is SKIPPED (passes through uncorrected), never crashed on.
- `indel_attempts`     : calls dispatched to the pair-HMM kernel.
- `indel_decodes`      : full, nontruncated, trace-valid pair-HMM decodes.
- `truncated_decodes`  : pair-HMM frontiers that returned only a decoded prefix;
  rejected before correction or soft-EM side effects.
- `indel_engaged`      : completed traces containing at least one I or D move.
- `trace_contract_errors` : malformed/missing pair-HMM traceback telemetry;
  an internal decoder contract regression rather than a data-dependent miss.
- `window_divergences` : substitution windows rejected for violating the
  length-preserving splice contract.
"""
mutable struct CorrectorDiagnostics
    structural_errors::Threads.Atomic{Int}
    unkmerizable_reads::Threads.Atomic{Int}
    # Calibrated gate requested (threshold set, substitution mode) but a decode
    # contract invariant was violated so gating silently fell open to the ungated
    # decode. On a healthy substitution decode this is always 0; a nonzero value
    # means the opt-in gate disabled itself — a regression signal, not data loss.
    gate_skipped::Threads.Atomic{Int}
    indel_attempts::Threads.Atomic{Int}
    indel_decodes::Threads.Atomic{Int}
    truncated_decodes::Threads.Atomic{Int}
    indel_engaged::Threads.Atomic{Int}
    trace_contract_errors::Threads.Atomic{Int}
    window_divergences::Threads.Atomic{Int}
end
function CorrectorDiagnostics()
    CorrectorDiagnostics(Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
        Threads.Atomic{Int}(0))
end

# Previously-proven-tractable finite beam (td-63qy: beam 256 completed on the
# 48 kb phage that OOM-crashed — ~21B allocations — at the unbounded typemax
# default). A read whose observation count is at/below the threshold keeps the
# EXACT maximum-likelihood guarantee (its frontier is small); above the threshold
# the frontier is bounded so the SHIPPING assembler cannot OOM-crash.
const _AUTO_BEAM_BOUNDED_WIDTH = 256
# Heuristic cutoff on a read's observation (k-mer) count. Typical short reads
# (Illumina, up to a few hundred bp) stay exact; long reads (PacBio/ONT, or the
# 48 kb-scale inputs behind td-63qy) get bounded. This is NOT a correctness
# constant: it only trades exactness for tractability, and only ABOVE the line.
const _AUTO_BEAM_EXACT_THRESHOLD = 1024

# Candidate-GENERATION bounds applied alongside the finite width beam (td-plqi).
# The width beam caps the RETAINED frontier at 256, but on a dense intermediate-k
# graph nearly all 256 retained states are >20 log-prob below the best — improbable
# junk that can never win the ML path yet still generates + emission-scores
# successors at the next depth, so per-depth generation climbs toward the width cap
# as the graph densifies (the empirically-measured residual super-linear term at
# k=9; #386). These bound generation directly:
#
#   * _AUTO_SUCCESSOR_BOUND — top-B outgoing transitions per expanded state, ranked
#     by edge weight before emission. A k-mer de Bruijn vertex has ≤ 4 structural
#     successors, so 16 is a strict no-op on DNA (a robustness guard for pathological
#     high-branching / non-DNA inputs), kept for correctness parity with the design.
#
#   * _AUTO_BEAM_SCORE_MARGIN — Δ log-prob "histogram" beam with an EMISSION
#     exemption: prune a frontier state only when it is >Δ below the depth's best
#     on BOTH the full score AND the cumulative emission (read-consistency). This
#     linearizes the dense-rung decode — the read-INCONSISTENT junk (low on both
#     axes) is discarded so the generating frontier stays O(1) in size — WITHOUT
#     dropping a real-but-rare allele: a supported minor allele in a skewed pool
#     has good emission but a coverage-driven transition penalty, and the emission
#     clause exempts it from pruning (the margin removes WRONG paths, not merely
#     RARE ones; PR #388 variation-safety review). Δ = 30 nats. Both engage ONLY
#     where the width beam is already finite (approximate); exact-ML reads
#     (beam_width == typemax) keep margin = Inf and stay byte-identical.
const _AUTO_SUCCESSOR_BOUND = 16
const _AUTO_BEAM_SCORE_MARGIN = 30.0

"""
    _auto_beam_width(n_observations) -> Int

Size-aware default beam width reconciling the ADR's "exact by default" with the
td-63qy OOM crash. Returns `typemax(Int)` (exact — the unbounded frontier
preserves the ML guarantee) when `n_observations <= _AUTO_BEAM_EXACT_THRESHOLD`,
and the bounded `_AUTO_BEAM_BOUNDED_WIDTH` above it, emitting a one-line `@info`
when it bounds so the exactness→tractability switch is visible in the log.
Callers can still FORCE exact on a large read by passing an explicit
`beam_width = typemax(Int)`.
"""
function _auto_beam_width(n_observations::Integer)::Int
    if n_observations <= _AUTO_BEAM_EXACT_THRESHOLD
        return typemax(Int)
    end
    @info "iterative corrector: read has $(n_observations) observations " *
          "(> $(_AUTO_BEAM_EXACT_THRESHOLD)); bounding Viterbi beam_width=" *
          "$(_AUTO_BEAM_BOUNDED_WIDTH) to stay tractable (exact ML would risk OOM, " *
          "td-63qy). Pass beam_width=typemax(Int) to force exact." maxlog = 3
    return _AUTO_BEAM_BOUNDED_WIDTH
end

"""
    _auto_beam_width(n_observations, n_vertices) -> Int

Graph-density-aware size-aware default beam width (td-35ux). The read-length-only
[`_auto_beam_width`](@ref) keeps the exact/unbounded frontier for any read with
`n_observations <= _AUTO_BEAM_EXACT_THRESHOLD` (typical short reads). That is
correct on a sparse graph, but at a dense intermediate k-rung the EXACT retained
(vertex, strand) frontier grows `O(n_vertices)` per decode depth — i.e.
`O(genome)` per pass, `O(genome^2)` over the read set — even for a single short
150 bp read (~15,723 frontier states at 5 kb / k=9, ~2.15 s per read). Long reads
already dodge this because the read-length rule bounds them; short reads did not,
which is why only short reads blew up (#376).

This two-arg form bounds the beam to `_AUTO_BEAM_BOUNDED_WIDTH` whenever EITHER
the read is large (read-length rule) OR the graph is dense
(`n_vertices > _AUTO_BEAM_BOUNDED_WIDTH`), so a short read on a big graph is
capped regardless of its length. On a small/sparse graph
(`n_vertices <= _AUTO_BEAM_BOUNDED_WIDTH`) the exact frontier is already tiny, so
a small read keeps the exact `typemax(Int)` (ML guarantee preserved). Callers can
still FORCE exact by passing an explicit `beam_width = typemax(Int)`.
"""
function _auto_beam_width(n_observations::Integer, n_vertices::Integer)::Int
    # Read-length rule first: a large read is bounded no matter the graph.
    by_length = _auto_beam_width(n_observations)
    if by_length != typemax(Int)
        return by_length
    end
    # The read is small (exact by read length). If the graph is dense the exact
    # frontier can still grow O(n_vertices) per depth → O(genome^2) per pass, so
    # bound it. On a sparse graph the exact frontier is tiny; keep it exact.
    if n_vertices > _AUTO_BEAM_BOUNDED_WIDTH
        @info "iterative corrector: short read ($(n_observations) observations) on a " *
              "dense graph ($(n_vertices) vertices > $(_AUTO_BEAM_BOUNDED_WIDTH)); " *
              "bounding Viterbi beam_width=$(_AUTO_BEAM_BOUNDED_WIDTH) to keep the " *
              "retained frontier O(1) instead of O(n_vertices) per depth (td-35ux). " *
              "Pass beam_width=typemax(Int) to force exact." maxlog = 3
        return _AUTO_BEAM_BOUNDED_WIDTH
    end
    return typemax(Int)
end

"""
Improve likelihood of entire read set using current graph and k-mer size.
Returns updated reads and count of improvements made.
Uses memory-efficient batch processing for large datasets.

`beam_width` defaults to `nothing` (size-aware auto-beam per read via
`_auto_beam_width`): exact ML where the read is small enough to be tractable,
a bounded frontier on large reads so the SHIPPING assembler cannot OOM-crash
(td-63qy). Pass an explicit `Int` (including `typemax(Int)`) to override — a
caller/test passing `beam_width = typemax(Int)` still gets exact decoding on
every read. Swallowed structural/un-k-merizable decode failures are tallied into
`diagnostics` (auto-created if not supplied) and a high swallowed fraction is
`@warn`ed so a run that silently corrected nothing is visible.

`hard_vertices` (td-nn6l): when non-`nothing`, restricts decoding to reads whose
k-mers overlap the "hard" vertex set (bubbles / repeats / weak k-mers); every
other read passes through untouched. Composes with `skip_solid` — a read is
decoded only when it is NOT all-solid AND (no hard set, or it touches one).

`windowed_decode` (td-nn6l, Stage 3c): when `true` AND a hard-window gate is
active (`hard_vertices !== nothing`), a read selected for decode is corrected
WINDOW-BY-WINDOW (`improve_read_likelihood_windowed`) — only the boundary-
constrained hard sub-window(s) around its hard vertices (each `<=500` bp) are
decoded and spliced back, instead of a whole-read Viterbi. This bounds each hard
read's decode to `O(max_window)` rather than `O(read_length)` (the #375 long-read
super-linear term). `false` (default) keeps whole-read decode as the fallback.

`decode_enabled` (td-9h5r): when `false`, the per-read graph-Viterbi decode is
skipped for EVERY read this pass (skip_fraction == 1.0) — the EXPLICIT low-k decode
floor. Stage 0 cheap correction still runs first (on `work_reads`), so simple
errors are still fixed; only the expensive whole-read Viterbi is deferred to higher
k-rungs. Defaults to `true`; the `:scalable` k-ladder sets it `false` below
`min_decode_k`.

`decode_gate_density` (td-9h5r): the ADAPTIVE low-k decode gate. When non-`nothing`
and a hard-window gate is active, the pass measures the post-Stage-0 fraction of
reads the gate would still decode; if that is `>= decode_gate_density` the gate is
NON-DISCRIMINATING (the dense-low-k "gate skips nothing" pathology) and the whole
decode is skipped for the pass. A genuinely selective low-k decode (gate skips a
meaningful fraction) is left ON. `nothing` (default / `:exhaustive`) disables it.

Returns `(updated_reads, improvements_made, skip_fraction, cheap_corrections,
decode_gated)` where `decode_gated` is `true` iff the per-read decode was gated OFF
for this pass (explicit floor OR adaptive density). Callers destructuring only the
first two/three/four values are unaffected.

`soft_weights` (td-e70t): when non-`nothing`, each decoded read's ML path
accumulates edge responsibilities into it (the soft-EM E-step). Accumulation is
NOT thread-safe, so supplying it forces sequential processing for this pass.

`cheap_correct` (td-bjnt): when `true` (and `graph_mode==:canonical`), run the
Stage 0 k-mer-spectrum correction pass (`_stage0_cheap_correct`) over the read set
BEFORE any gating/decode, cheaply fixing simple single-substitution errors with a
linear scan so the expensive graph Viterbi is reserved for genuine ambiguity
(bubbles/repeats). The decode then operates on the cheaply-corrected reads. Only
enabled on the :scalable tier.

Returns `(updated_reads, improvements_made, skip_fraction, cheap_corrections,
decode_gated)` where `skip_fraction` is the fraction of reads passed through WITHOUT
a decode (solid + hard-window skips, plus the whole set when the low-k decode gate
fires) and `cheap_corrections` is the number of bases fixed by the Stage 0 pass. The
graph-Viterbi decode fraction is `1 - skip_fraction`. `improvements_made` counts
BOTH cheap corrections and decode-accepted corrections so per-k convergence sees
Stage 0 progress. `decode_gated` is `true` iff the per-read decode was gated OFF for
the whole pass (low-k gate, td-9h5r). Callers destructuring only the first
two/three/four values are unaffected.
"""
function improve_read_set_likelihood(reads::Vector{<:FASTX.FASTQ.Record}, graph, k::Int;
        verbose::Bool = false,
        batch_size::Int = 10000,
        enable_parallel::Bool = false,
        graph_mode::Symbol = :canonical,
        skip_solid::Bool = false,
        cheap_correct::Bool = false,
        beam_width::Union{Int, Nothing} = nothing,
        calibrated_gap_threshold::Union{Float64, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        prior_soft_weights::Union{
            Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        hard_vertices::Union{Nothing, AbstractSet} = nothing,
        windowed_decode::Bool = false,
        decode_enabled::Bool = true,
        decode_gate_density::Union{Float64, Nothing} = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing,
        indel_params::Union{Nothing, IndelDecodeParams} = nothing,
        indel_schedule::Symbol = :unrestricted,
        rung_telemetry::Union{Nothing, Dict{Symbol, Any}} = nothing,
        memory_limit::Union{Nothing, Int} = nothing)::Tuple{
        Vector{FASTX.FASTQ.Record}, Int, Float64, Int, Bool}
    return Mycelia.Rhizomorph._with_soft_edge_weight_scope(
        graph, prior_soft_weights) do
        _improve_read_set_likelihood_impl(
            reads,
            graph,
            k;
            verbose = verbose,
            batch_size = batch_size,
            enable_parallel = enable_parallel,
            graph_mode = graph_mode,
            skip_solid = skip_solid,
            cheap_correct = cheap_correct,
            beam_width = beam_width,
            calibrated_gap_threshold = calibrated_gap_threshold,
            soft_weights = soft_weights,
            prior_soft_weights = prior_soft_weights,
            hard_vertices = hard_vertices,
            windowed_decode = windowed_decode,
            decode_enabled = decode_enabled,
            decode_gate_density = decode_gate_density,
            diagnostics = diagnostics,
            indel_params = indel_params,
            indel_schedule = indel_schedule,
            rung_telemetry = rung_telemetry,
            memory_limit = memory_limit,
        )
    end
end

function _improve_read_set_likelihood_impl(
        reads::Vector{<:FASTX.FASTQ.Record}, graph, k::Int;
        verbose::Bool,
        batch_size::Int,
        enable_parallel::Bool,
        graph_mode::Symbol,
        skip_solid::Bool,
        cheap_correct::Bool,
        beam_width::Union{Int, Nothing},
        calibrated_gap_threshold::Union{Float64, Nothing},
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator},
        prior_soft_weights::Union{
            Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator},
        hard_vertices::Union{Nothing, AbstractSet},
        windowed_decode::Bool,
        decode_enabled::Bool,
        decode_gate_density::Union{Float64, Nothing},
        diagnostics::Union{Nothing, CorrectorDiagnostics},
        indel_params::Union{Nothing, IndelDecodeParams},
        indel_schedule::Symbol,
        rung_telemetry::Union{Nothing, Dict{Symbol, Any}},
        memory_limit::Union{Nothing, Int},
)::Tuple{Vector{FASTX.FASTQ.Record}, Int, Float64, Int, Bool}
    indel_schedule in (:unrestricted, :frontier_budgeted) ||
        throw(ArgumentError("indel_schedule must be :unrestricted or " *
                            ":frontier_budgeted; got :$(indel_schedule)"))
    diag = diagnostics === nothing ? CorrectorDiagnostics() : diagnostics
    # Snapshot so the per-call @warn reflects THIS pass's swallowed fraction even
    # when `diag` is a shared accumulator threaded across many passes.
    struct_before = diag.structural_errors[]
    unk_before = diag.unkmerizable_reads[]
    indel_attempts_before = diag.indel_attempts[]
    indel_decodes_before = diag.indel_decodes[]
    truncated_before = diag.truncated_decodes[]
    indel_engaged_before = diag.indel_engaged[]
    trace_contract_before = diag.trace_contract_errors[]
    gate_skipped_before = diag.gate_skipped[]
    window_divergences_before = diag.window_divergences[]
    total_reads = length(reads)
    updated_reads = Vector{FASTX.FASTQ.Record}(undef, total_reads)
    improvements_made = 0

    # Stage 0 skip (td-1do7): classify k-mers once; a read whose every k-mer is
    # solid has no weak region to correct and skips the per-read decode. Opt-in
    # (default off) so existing callers are unchanged. `solid_kmers === nothing`
    # ⇒ correct every read as before.
    # td-nt69: the skip/cheap-correct helpers resolve read k-mers to their graph
    # label via `_lookup_key` (canonicalized under :canonical, observed orientation
    # under :doublestrand), so they are valid on BOTH modes — classification is
    # coverage-based and mode-agnostic. :scalable now runs :doublestrand. Only
    # :singlestrand (no RC vertices, not a corrector target) is excluded.
    if skip_solid && graph_mode == :singlestrand
        @warn "skip_solid is not supported for graph_mode=:singlestrand; disabling skip (correcting all reads)." graph_mode
    end
    if cheap_correct && graph_mode == :singlestrand
        @warn "cheap_correct is not supported for graph_mode=:singlestrand; disabling Stage 0 cheap correction." graph_mode
    end
    # Classify k-mers once if EITHER the skip-solid gate OR the Stage 0 cheap
    # corrector needs the solid set. Compute once, share both consumers.
    need_solid = (skip_solid || cheap_correct) && graph_mode != :singlestrand
    solid_kmers = need_solid ? _solid_kmer_set(graph) : nothing

    # Stage 0 CHEAP correction (td-bjnt): fix simple single-substitution errors
    # with a linear k-mer-spectrum scan BEFORE any gating/decode, so graph Viterbi
    # is reserved for genuine ambiguity. The decode below then runs on the
    # cheaply-corrected reads. Gated on :scalable (cheap_correct=true, canonical).
    cheap_corrections = 0
    work_reads = reads
    if cheap_correct && graph_mode != :singlestrand && solid_kmers !== nothing
        work_reads,
        cheap_corrections = _stage0_cheap_correct(reads, k, solid_kmers;
            graph_mode = graph_mode)
        improvements_made += cheap_corrections
        if verbose && cheap_corrections > 0
            println("  Stage 0 cheap correction: fixed $cheap_corrections base(s) " *
                    "before decode")
        end
    end
    skipped_reads = 0

    # A read is SKIPPED (passed through uncorrected, no decode) when it is either
    # all-solid (skip-solid gate) OR — under hard-window gating — it touches no
    # hard vertex. Both are volume-reduction skips and both count toward the
    # reported skip fraction. `solid_kmers` may be populated for cheap correction
    # even when skip_solid is off, so the all-solid skip is gated on `skip_solid`
    # explicitly. Returns true ⇒ skip. Stage 0 has ALREADY run above (on
    # `work_reads`), so this predicate sees the cheaply-corrected reads.
    _gate_skip = read -> begin
        (skip_solid && solid_kmers !== nothing &&
         _read_is_all_solid(read, k, solid_kmers; graph_mode = graph_mode)) &&
            return true
        (hard_vertices !== nothing &&
         !should_decode_read(read, k, hard_vertices; graph_mode = graph_mode)) &&
            return true
        return false
    end

    # Precompute per-read skip decisions ONCE (cheap k-mer-membership checks) so the
    # adaptive low-k gate can measure the natural decode fraction and the decode loop
    # can reuse the same flags (no double evaluation).
    base_skip_flags = Vector{Bool}(undef, total_reads)
    @inbounds for i in 1:total_reads
        base_skip_flags[i] = _gate_skip(work_reads[i])
    end
    natural_decode_fraction = total_reads > 0 ?
                              count(!, base_skip_flags) / total_reads : 0.0

    # Branching/frontier-aware indel schedule (td-jt7r.2). Profile intent is
    # distinct from runtime admission: `requested` counts the post-Stage-0 hard
    # windows that asked for pair-HMM service, while `attempted` is incremented
    # only immediately before the kernel. Rejected windows execute the existing
    # raw-graph substitution path; admitted windows use either the raw graph or
    # a conservatively cleaned COPY, never mutating the raw/reusable graph.
    pass_indel_telemetry = rung_telemetry === nothing ?
                           _empty_indel_rung_telemetry(indel_params !== nothing) :
                           rung_telemetry
    empty!(pass_indel_telemetry)
    merge!(pass_indel_telemetry, _empty_indel_rung_telemetry(indel_params !== nothing))
    decode_graph = graph
    scheduled_weighted_graph = nothing
    scheduled_cleaned_graph = nothing
    scheduled_cleaned_weighted_graph = nothing
    scheduled_window_sources = Dict{Int, Dict{UnitRange{Int}, Symbol}}()
    effective_indel_params = indel_params
    indel_admitted = indel_params !== nothing && indel_schedule == :unrestricted
    if indel_params === nothing
        pass_indel_telemetry[:decision_reason] = :profile_disabled
    elseif indel_schedule == :unrestricted
        pass_indel_telemetry[:requested] = count(!, base_skip_flags)
        pass_indel_telemetry[:admitted] = true
        pass_indel_telemetry[:decision_reason] = :unrestricted_semantics
    elseif !decode_enabled
        effective_indel_params = nothing
        indel_admitted = false
        pass_indel_telemetry[:decision_reason] = :explicit_decode_floor
    elseif hard_vertices === nothing || !windowed_decode
        effective_indel_params = nothing
        indel_admitted = false
        pass_indel_telemetry[:decision_reason] = :bounded_windows_unavailable
    else
        decision = _evaluate_indel_frontier_schedule(
            work_reads, graph, k, hard_vertices, base_skip_flags,
            indel_params, graph_mode;
            prior_soft_weights = prior_soft_weights,
            memory_limit = memory_limit,
        )
        pass_indel_telemetry[:requested] = decision.requested
        pass_indel_telemetry[:admitted] = decision.admitted
        pass_indel_telemetry[:admitted_windows] = decision.admitted_windows
        pass_indel_telemetry[:rejected_windows] = decision.rejected_windows
        pass_indel_telemetry[:graph_source] = decision.graph_source
        pass_indel_telemetry[:decision_reason] = decision.decision_reason
        pass_indel_telemetry[:raw_frontier_evaluated] = decision.raw_evaluated
        pass_indel_telemetry[:cleaned_frontier_evaluated] = decision.cleaned_evaluated
        pass_indel_telemetry[:raw_frontier_metrics] = decision.raw_metrics
        pass_indel_telemetry[:cleaned_frontier_metrics] = decision.cleaned_metrics
        pass_indel_telemetry[:graph_cleanup] = decision.cleanup
        indel_admitted = decision.admitted
        # Retain the full map even when every pair-HMM window is rejected: those
        # entries explicitly route the same bounded windows through substitution.
        scheduled_window_sources = decision.window_sources
        pass_indel_telemetry[:bounded_windowing_forced] =
            !isempty(scheduled_window_sources)
        if indel_admitted
            scheduled_weighted_graph = decision.raw_weighted_graph
            scheduled_cleaned_graph = decision.cleaned_graph
            scheduled_cleaned_weighted_graph = decision.cleaned_weighted_graph
            if prior_soft_weights !== nothing && scheduled_cleaned_graph !== nothing
                Mycelia.Rhizomorph.register_soft_edge_weights!(
                    scheduled_cleaned_graph, prior_soft_weights)
            end
        else
            effective_indel_params = nothing
        end
    end

    # -- Low-k decode gating (td-9h5r) -----------------------------------------
    # `decode_enabled == false` (EXPLICIT floor `min_decode_k`): skip the per-read
    # graph-Viterbi decode for EVERY read this pass. Stage 0 cheap correction has
    # already run (on `work_reads`), so error correction still happens; only the
    # expensive whole-read Viterbi is deferred to higher, selective k-rungs.
    #
    # ADAPTIVE density gate: even with `decode_enabled`, if the hard-window gate
    # would still send >= `decode_gate_density` of reads to the decode (i.e. it is
    # NON-DISCRIMINATING — the dense-low-k pathology where "the gate skips nothing"),
    # skip the whole decode too: the gate is not reducing volume, so the decode is
    # ~pure waste Stage 0 + the selective higher rungs recover. Requires an active
    # hard-window gate (`hard_vertices !== nothing`) — with no gate there is no
    # discrimination signal to act on, so the decode runs as before.
    adaptive_gated = decode_enabled && !indel_admitted &&
                     isempty(scheduled_window_sources) &&
                     decode_gate_density !== nothing &&
                     hard_vertices !== nothing &&
                     natural_decode_fraction >= decode_gate_density
    low_k_decode_gated = !decode_enabled || adaptive_gated
    pass_decode_off = low_k_decode_gated

    scheduled_sources = Set(
        source for sources in values(scheduled_window_sources)
        for source in values(sources))
    needs_raw_weighted_graph = isempty(scheduled_sources) ||
                               :raw in scheduled_sources ||
                               :substitution in scheduled_sources
    # A fully rejected frontier schedule still promises bounded substitution
    # fallback, which needs one raw weighted copy. Honor the configured memory
    # ceiling before allocating it. Explicit :unrestricted semantics and the
    # profile-disabled legacy path bypass this scheduler-specific gate.
    substitution_memory_gated =
        indel_schedule == :frontier_budgeted &&
        !isempty(scheduled_sources) &&
        needs_raw_weighted_graph &&
        scheduled_weighted_graph === nothing &&
        !_indel_graph_variants_fit_memory(decode_graph, 2, memory_limit)
    if substitution_memory_gated
        pass_decode_off = true
        pass_indel_telemetry[:substitution_decode_memory_gated] = true
    end
    if verbose && adaptive_gated
        println("  Low-k decode gate (td-9h5r): hard-window gate non-discriminating " *
                "at k=$k (would decode $(round(natural_decode_fraction * 100, digits=1))% " *
                ">= $(round(decode_gate_density * 100, digits=1))%); skipping decode " *
                "this pass (Stage 0 + higher-k rungs recover).")
    end

    # Final per-read skip decision: all reads skip when the pass decode is gated off;
    # otherwise use the precomputed gate flags.
    _skip_this_read_at = i -> pass_decode_off || base_skip_flags[i]

    # Soft-EM accumulation into `soft_weights` is not thread-safe (shared Dict),
    # so a soft-EM pass runs sequentially. Otherwise honor the caller's request.
    use_parallel = enable_parallel && Threads.nthreads() > 1 && soft_weights === nothing
    if enable_parallel && soft_weights !== nothing && Threads.nthreads() > 1
        @warn "soft-EM edge accumulation is sequential (race-free); ignoring enable_parallel for this pass." maxlog = 1
    end

    if verbose
        println("  Processing $total_reads reads in batches of $batch_size")
    end

    # Hoist the weighted-decode graph OUT of the per-read loop (td-y4oj). The
    # per-read `correct_observations` previously rebuilt an O(V+E) weighted COPY
    # of the whole graph ONCE PER READ; since n_reads ∝ genome and V ∝ genome,
    # that made the decode O(genome^2) (the dominant #372 super-linear term,
    # alpha≈2.14). Build it once here and thread the shared, read-only instance
    # into every per-read decode so the rebuild cost is amortized across the pass
    # (linear in n_reads). Returns `nothing` for graphs without a weighted-decode
    # path (legacy / non-MetaGraph), which the callee handles by falling back to
    # its unchanged per-read behavior. Skip the build entirely when the pass
    # decode is gated off (no read decodes this pass).
    pass_weighted_graph = if pass_decode_off || !needs_raw_weighted_graph
        nothing
    elseif scheduled_weighted_graph === nothing
        try
            build_correction_weighted_graph(decode_graph)
        catch exception
            if indel_schedule != :frontier_budgeted ||
               !(exception isa OutOfMemoryError)
                rethrow()
            end
            # The preflight estimate is conservative but cannot guarantee an
            # allocation. Fail closed rather than retrying the same raw copy after
            # a cleaning/weighted OOM; the low-k gate flag remains semantically
            # separate and the dedicated telemetry records this memory cause.
            pass_decode_off = true
            substitution_memory_gated = true
            pass_indel_telemetry[:substitution_decode_memory_gated] = true
            GC.gc(true)
            nothing
        end
    else
        scheduled_weighted_graph
    end

    # Stage 3c windowed decode (td-nn6l): when enabled AND a hard-window gate is
    # active, a LONG read selected for decode is corrected WINDOW-BY-WINDOW (only
    # the hard sub-window(s) around its hard vertices, each boundary-constrained and
    # <=500 bp) instead of whole-read. Requires `hard_vertices` — without the gate
    # there are no hard windows to target. The per-read `_windowed_decode_read_is_long`
    # check restricts windowing to reads long enough that whole-read decode is
    # expensive/bounded; short reads keep the cheap, exact whole-read decode
    # (so the :scalable short-read quality gate never regresses).
    use_windowed = windowed_decode && hard_vertices !== nothing
    force_indel_windowing = indel_schedule == :frontier_budgeted &&
                            !isempty(scheduled_window_sources)

    # Process in batches for memory efficiency. `work_reads` is the Stage 0
    # cheaply-corrected read set (== `reads` when cheap_correct is off), so the
    # decode operates on already-simplified reads.
    for batch_start in 1:batch_size:total_reads
        batch_end = min(batch_start + batch_size - 1, total_reads)
        batch_reads = work_reads[batch_start:batch_end]

        batch_improvements = 0

        if use_parallel
            # Parallel processing for large batches
            batch_results = Vector{Tuple{FASTX.FASTQ.Record, Bool}}(undef, length(batch_reads))

            # Vector{Bool} (one byte/elem), NOT falses()/BitVector (64 bits share a
            # word) — the @threads writes below would otherwise race on the shared
            # word and undercount the skip fraction (review I2).
            skip_flags = fill(false, length(batch_reads))
            Threads.@threads for i in eachindex(batch_reads)
                read = batch_reads[i]
                if _skip_this_read_at(batch_start + i - 1)
                    batch_results[i] = (read, false)   # skip the decode
                    skip_flags[i] = true
                else
                    read_index = batch_start + i - 1
                    read_window_sources = get(
                        scheduled_window_sources, read_index, nothing)
                    improved_read,
                    was_improved = (use_windowed &&
                                    (force_indel_windowing ||
                                     _windowed_decode_read_is_long(read, k))) ?
                                   improve_read_likelihood_windowed(
                        read, decode_graph, k, hard_vertices; graph_mode = graph_mode,
                        beam_width = beam_width, weighted_graph = pass_weighted_graph,
                        cleaned_graph = scheduled_cleaned_graph,
                        cleaned_weighted_graph = scheduled_cleaned_weighted_graph,
                        indel_window_sources = read_window_sources,
                        diagnostics = diag,
                        calibrated_gap_threshold = calibrated_gap_threshold,
                        indel_params = effective_indel_params) :
                                   improve_read_likelihood(
                        read, decode_graph, k; graph_mode = graph_mode,
                        beam_width = beam_width, weighted_graph = pass_weighted_graph,
                        diagnostics = diag, indel_params = effective_indel_params,
                        calibrated_gap_threshold = calibrated_gap_threshold)
                    batch_results[i] = (improved_read, was_improved)
                end
            end
            skipped_reads += count(skip_flags)

            # Collect results
            for (i, (improved_read, was_improved)) in enumerate(batch_results)
                updated_reads[batch_start + i - 1] = improved_read
                if was_improved
                    batch_improvements += 1
                end
            end
        else
            # Sequential processing (also the soft-EM path: accumulation into
            # `soft_weights` happens per-read inside improve_read_likelihood).
            for (i, read) in enumerate(batch_reads)
                if _skip_this_read_at(batch_start + i - 1)
                    updated_reads[batch_start + i - 1] = read   # skip the decode
                    skipped_reads += 1
                    continue
                end
                read_index = batch_start + i - 1
                read_window_sources = get(
                    scheduled_window_sources, read_index, nothing)
                improved_read,
                was_improved = (use_windowed &&
                                (force_indel_windowing ||
                                 _windowed_decode_read_is_long(read, k))) ?
                               improve_read_likelihood_windowed(
                    read, decode_graph, k, hard_vertices; graph_mode = graph_mode,
                    beam_width = beam_width, soft_weights = soft_weights,
                    weighted_graph = pass_weighted_graph,
                    cleaned_graph = scheduled_cleaned_graph,
                    cleaned_weighted_graph = scheduled_cleaned_weighted_graph,
                    indel_window_sources = read_window_sources,
                    diagnostics = diag,
                    calibrated_gap_threshold = calibrated_gap_threshold,
                    indel_params = effective_indel_params) :
                               improve_read_likelihood(
                    read, decode_graph, k; graph_mode = graph_mode,
                    beam_width = beam_width, soft_weights = soft_weights,
                    weighted_graph = pass_weighted_graph,
                    diagnostics = diag, indel_params = effective_indel_params,
                    calibrated_gap_threshold = calibrated_gap_threshold)
                updated_reads[batch_start + i - 1] = improved_read

                if was_improved
                    batch_improvements += 1
                end
            end
        end

        improvements_made += batch_improvements

        # Progress reporting and garbage collection
        if verbose
            batch_num = div(batch_start - 1, batch_size) + 1
            total_batches = div(total_reads - 1, batch_size) + 1
            improvement_rate = batch_improvements / length(batch_reads) * 100
            println("    Batch $batch_num/$total_batches: $(batch_improvements)/$(length(batch_reads)) improvements ($(round(improvement_rate, digits=1))%)")
        end

        # Force garbage collection between batches for memory efficiency
        if batch_end < total_reads
            GC.gc()
        end
    end

    skip_fraction = total_reads > 0 ? skipped_reads / total_reads : 0.0
    if verbose
        total_improvement_rate = improvements_made / total_reads * 100
        println("  Total improvements: $improvements_made/$total_reads ($(round(total_improvement_rate, digits=1))%)")
        if skip_solid || hard_vertices !== nothing || pass_decode_off
            println("  Skipped (no decode): $skipped_reads/$total_reads ($(round(skip_fraction * 100, digits=1))%)")
            decode_fraction = 1.0 - skip_fraction
            gated = if low_k_decode_gated
                " [low-k decode gated OFF, td-9h5r]"
            elseif substitution_memory_gated
                " [substitution decode memory gated OFF]"
            else
                ""
            end
            println("  Graph-Viterbi decode fraction: $(round(decode_fraction * 100, digits=1))%$(gated)")
        end
        if cheap_correct
            println("  Stage 0 cheap corrections (this pass): $cheap_corrections")
        end
    end

    # Surface swallowed decode failures for THIS pass (see CorrectorDiagnostics):
    # a systematically broken corrector otherwise reports 0 improvements as if it
    # ran. Distinguish "every read failed structurally" (the corrector is not
    # running at all) from a merely-high swallowed fraction.
    struct_swallowed = diag.structural_errors[] - struct_before
    unk_swallowed = diag.unkmerizable_reads[] - unk_before
    swallowed = struct_swallowed + unk_swallowed
    if total_reads > 0 && swallowed > 0
        frac = swallowed / total_reads
        if struct_swallowed == total_reads
            @warn "iterative corrector: EVERY read failed to decode with a structural " *
                  "error (empty graph / alphabet-inference miss / config error) — the " *
                  "corrector is not running; 0 improvements is NOT 'nothing to fix'." total_reads structural_errors = struct_swallowed
        elseif frac >= 0.5
            @warn "iterative corrector: high swallowed-decode fraction; many reads "*
            "passed through uncorrected (not conflated with 'no likelihood gain')." total_reads structural_errors=struct_swallowed unkmerizable_reads=unk_swallowed fraction=round(
                frac, digits = 3)
        end
    end

    trace_contract_this_pass = diag.trace_contract_errors[] - trace_contract_before
    if trace_contract_this_pass > 0
        @warn "iterative corrector: pair-HMM traceback contract failed; rejected " *
              "decode(s) before correction or soft-EM side effects." total_reads trace_contract_errors = trace_contract_this_pass
    end
    truncated_this_pass = diag.truncated_decodes[] - truncated_before
    successful_indel_this_pass = diag.indel_decodes[] - indel_decodes_before
    completed_indel_outcomes = successful_indel_this_pass + truncated_this_pass +
                               trace_contract_this_pass
    if completed_indel_outcomes > 0 &&
       truncated_this_pass / completed_indel_outcomes >= 0.5
        @warn "iterative corrector: high truncated pair-HMM decode fraction; prefix-only " *
              "paths were rejected before correction or soft-EM side effects." total_reads completed_indel_outcomes truncated_decodes = truncated_this_pass fraction = round(
            truncated_this_pass / completed_indel_outcomes, digits = 3)
    end

    # A requested calibrated gate that silently fell open on a substitution decode
    # is a contract regression (the gate did nothing despite being opted in) — see
    # CorrectorDiagnostics.gate_skipped. Surface it rather than let the frontier look
    # like "the gate doesn't help."
    gate_skipped_this_pass = diag.gate_skipped[] - gate_skipped_before
    if gate_skipped_this_pass > 0
        @warn "iterative corrector: calibrated gate fell open (ungated) on " *
              "$(gate_skipped_this_pass) read(s) despite a substitution decode — the " *
              "gap/path alignment contract was violated; the gate silently did nothing." total_reads gate_skipped = gate_skipped_this_pass
    end

    window_divergences_this_pass =
        diag.window_divergences[] - window_divergences_before
    if window_divergences_this_pass > 0
        @warn "iterative corrector: rejected substitution window(s) whose decoded " *
              "length violated the length-preserving splice contract." total_reads window_divergences = window_divergences_this_pass
    end

    pass_indel_telemetry[:attempted] =
        diag.indel_attempts[] - indel_attempts_before
    pass_indel_telemetry[:completed] =
        diag.indel_decodes[] - indel_decodes_before
    pass_indel_telemetry[:truncated] =
        diag.truncated_decodes[] - truncated_before
    pass_indel_telemetry[:engaged] =
        diag.indel_engaged[] - indel_engaged_before

    return updated_reads, improvements_made, skip_fraction, cheap_corrections,
    low_k_decode_gated
end

"""
Improve likelihood of a single read using maximum likelihood path finding.
Returns improved read and boolean indicating if improvement was made.

`soft_weights` (td-e70t): when supplied, the read's decoded ML path accumulates
edge responsibilities into it (responsibility 1.0 for the single argmax path).
"""
function improve_read_likelihood(read::FASTX.FASTQ.Record, graph, k::Int;
        graph_mode::Symbol = :canonical,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        weighted_graph = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing,
        indel_params::Union{Nothing, IndelDecodeParams} = nothing,
        calibrated_gap_threshold::Union{Float64, Nothing} = nothing)::Tuple{
        FASTX.FASTQ.Record, Bool}
    # Extract sequence and quality
    original_seq = FASTX.sequence(String, read)
    original_qual = FASTX.quality(read)
    read_id = FASTX.identifier(read)

    # Skip reads shorter than k
    if length(original_seq) < k
        return read, false
    end

    # Find optimal path through graph via the trustworthy Viterbi. `beam_width`
    # defaults to `nothing` = size-aware auto-beam (exact where tractable, bounded
    # above the threshold to avoid the td-63qy OOM); an explicit Int overrides it.
    #
    # Two-gate design (intentional, not a narrative/code mismatch): the inner
    # `find_optimal_sequence_path` gates on `improvement > 0.0` (accept the ML path
    # only when it is at least as likely as the original — never return a strictly
    # worse read), while THIS function additionally requires `> 0.01` so a
    # negligible likelihood wiggle is not counted as a real correction. The two
    # thresholds serve different jobs (correctness floor vs. meaningful-change
    # floor) and compose deliberately.
    improved_read,
    likelihood_improvement = find_optimal_sequence_path(
        read, graph, k; graph_mode = graph_mode,
        beam_width = beam_width, soft_weights = soft_weights,
        weighted_graph = weighted_graph,
        diagnostics = diagnostics, indel_params = indel_params,
        calibrated_gap_threshold = calibrated_gap_threshold)

    # Only update if significant improvement
    if likelihood_improvement > 0.01  # Threshold for meaningful improvement
        return improved_read, true
    else
        return read, false
    end
end

"""
Find optimal sequence path through graph using maximum likelihood principles.
Returns improved sequence and likelihood improvement score.
"""
function find_optimal_sequence_path(sequence::AbstractString, quality, graph, k::Int;
        graph_mode::Symbol = :canonical)::Tuple{String, Float64}
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for find_optimal_sequence_path"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    improved_record,
    improvement = find_optimal_sequence_path(record, graph, k; graph_mode = graph_mode)
    return FASTX.sequence(String, improved_record), improvement
end

function find_optimal_sequence_path(read::FASTX.FASTQ.Record, graph, k::Int;
        graph_mode::Symbol = :canonical,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        weighted_graph = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing,
        indel_params::Union{Nothing, IndelDecodeParams} = nothing,
        calibrated_gap_threshold::Union{Float64, Nothing} = nothing)::Tuple{
        FASTX.FASTQ.Record, Float64}
    # Trustworthy Viterbi (graph-as-HMM correction core, td-ak6w). Trust the
    # maximum-likelihood path or report no improvement. The prior heuristic
    # fallback cascade (0.01 abandon-threshold -> statistical resampling -> local
    # edits) is removed from the correctness path: a "corrected" read must lie on
    # an ML graph path, not on an off-path heuristic that may not correspond to
    # any real traversal. `try_statistical_path_resampling` and
    # `try_local_path_improvements` remain DEFINED as opt-in comparison baselines
    # (and for direct callers/tests) but are no longer invoked here.
    #
    # `calculate_read_likelihood` k-merizes the read with the 2-bit iterator and
    # runs BEFORE the decoder, so an un-k-merizable read (ambiguous base, e.g. `N`)
    # throws `BioSequences.EncodeError` HERE, not in `try_viterbi_path_improvement`.
    # Treat it identically: SKIP the read (pass through uncorrected) + count, so a
    # single N-containing read never crashes the (threaded) corrector (td-63qy
    # family). Returning before the decoder means each such read is tallied once.
    original_likelihood = try
        calculate_read_likelihood(read, graph, k; graph_mode = graph_mode)
    catch error
        if error isa BioSequences.EncodeError
            diagnostics === nothing ||
                Threads.atomic_add!(diagnostics.unkmerizable_reads, 1)
            return read, 0.0
        end
        rethrow()
    end

    viterbi_result = try_viterbi_path_improvement(
        read, graph, k; graph_mode = graph_mode,
        beam_width = beam_width, soft_weights = soft_weights,
        weighted_graph = weighted_graph,
        diagnostics = diagnostics, indel_params = indel_params,
        calibrated_gap_threshold = calibrated_gap_threshold)
    if viterbi_result === nothing
        # Viterbi could not decode (empty observation set, structural error, or an
        # un-k-merizable read): report no improvement rather than falling back to a
        # heuristic path. Genuine decode FAILURES (vs. "no observations") are
        # tallied in `diagnostics` inside try_viterbi_path_improvement so they are
        # not silently conflated with "the ML path had no gain".
        return read, 0.0
    end

    viterbi_read, viterbi_likelihood = viterbi_result
    improvement = _comparable_likelihood_improvement(
        original_likelihood,
        viterbi_likelihood,
        length(FASTX.sequence(read)),
        length(FASTX.sequence(viterbi_read)),
        k;
        length_normalized = indel_params !== nothing)
    # Accept the ML path only when it is at least as likely as the original; a
    # non-positive delta means the ML path is not an improvement, so report none
    # (never return a strictly-worse read). This is the ADR's "trust the ML path
    # or report no improvement", NOT the removed 0.01 minimum-magnitude threshold.
    return improvement > 0.0 ? (viterbi_read, improvement) : (read, 0.0)
end

"""
Calculate likelihood of a FASTQ read given the current graph.
"""
function calculate_read_likelihood(
        read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol = :canonical)::Float64
    # Detect sequence type dynamically using existing alphabets.jl functions
    sequence_string = FASTX.sequence(String, read)
    alphabet = Mycelia.detect_alphabet(sequence_string)
    sequence_type = Mycelia.alphabet_to_biosequence_type(alphabet)

    sequence = Mycelia.extract_typed_sequence(read, sequence_type)
    quality_scores = collect(FASTX.quality_scores(read))
    return calculate_sequence_likelihood(
        sequence, quality_scores, graph, k; graph_mode = graph_mode)
end

"""
    _comparable_likelihood_improvement(original_likelihood, corrected_likelihood,
        original_length, corrected_length, k; length_normalized=false)

Compare graph log-likelihoods without giving a mechanical advantage to a shorter
length-changing correction. Substitution mode keeps the historical raw-sum
difference exactly. Indel mode compares mean log-likelihood per emitted k-mer;
the pair-HMM has already applied its affine gap priors while choosing the decoded
path, and this secondary acceptance gate therefore tests graph support on a
length-comparable scale.
"""
function _comparable_likelihood_improvement(
        original_likelihood::Float64,
        corrected_likelihood::Float64,
        original_length::Int,
        corrected_length::Int,
        k::Int;
        length_normalized::Bool = false
)::Float64
    if !length_normalized
        return corrected_likelihood - original_likelihood
    end
    original_terms = max(original_length - k + 1, 1)
    corrected_terms = max(corrected_length - k + 1, 1)
    return corrected_likelihood / corrected_terms -
           original_likelihood / original_terms
end

"""
Calculate quality-aware likelihood of a sequence given the qualmer graph.
Uses both k-mer presence and quality-based confidence from qualmer observations.
"""
function calculate_sequence_likelihood(sequence::AbstractString, quality, graph, k::Int;
        graph_mode::Symbol = :canonical)::Float64
    alphabet = detect_alphabet(sequence)
    typed_sequence = detect_and_extract_sequence(sequence, alphabet)
    quality_scores = if quality isa AbstractString
        Int8.(Int.(collect(codeunits(quality))) .- 33)
    elseif quality isa AbstractVector{<:Integer}
        Int8.(quality)
    else
        throw(ArgumentError("Unsupported quality type for calculate_sequence_likelihood"))
    end
    return calculate_sequence_likelihood(
        typed_sequence, quality_scores, graph, k; graph_mode = graph_mode)
end

function _resolve_kmer_label(graph, kmer; graph_mode::Symbol = :canonical)
    if haskey(graph, kmer)
        return kmer
    end
    if graph_mode == :canonical && kmer isa Union{Kmers.DNAKmer, Kmers.RNAKmer}
        canonical_kmer = BioSequences.canonical(kmer)
        if haskey(graph, canonical_kmer)
            return canonical_kmer
        end
    end
    return nothing
end

function _resolve_qualmer_for_graph(graph, qmer::Qualmer; graph_mode::Symbol = :canonical)
    resolved_kmer = _resolve_kmer_label(graph, qmer.kmer; graph_mode = graph_mode)
    if resolved_kmer === nothing
        return nothing
    end
    if resolved_kmer == qmer.kmer
        return qmer
    end
    if graph_mode == :canonical && qmer.kmer isa Union{Kmers.DNAKmer, Kmers.RNAKmer}
        return canonical(qmer)
    end
    return Qualmer(resolved_kmer, qmer.qualities)
end

"""
Resize a per-base quality vector to `target_len`, preserving the 1:1 position
correspondence a k-mer likelihood walk assumes. Extra tail positions inherit the
last observed quality (a neutral, information-preserving pad); surplus qualities
are dropped. A no-op when the lengths already match. This reconciles the
indel-aware corrector's output (an insertion lengthens the corrected sequence, a
deletion shortens it, breaking the original read's base↔quality correspondence)
with `calculate_sequence_likelihood`, which indexes `quality[pos:pos+k-1]` over
the CORRECTED sequence. Mirrors the length-handling `adjust_quality_scores`
already performs on the same corrected sequence.
"""
function _reconcile_quality_length(quality::AbstractVector{T}, target_len::Int) where {T}
    n = length(quality)
    n == target_len && return quality
    n == 0 && return fill(zero(T), target_len)
    reconciled = Vector{T}(undef, target_len)
    @inbounds for i in 1:target_len
        reconciled[i] = quality[clamp(i, 1, n)]
    end
    return reconciled
end

function calculate_sequence_likelihood(
        sequence::BioSequences.BioSequence, quality::Vector{Int8},
        graph, k::Int; graph_mode::Symbol = :canonical)::Float64
    if length(sequence) < k
        return 0.0
    end

    # Indel-aware correction can hand us a corrected `sequence` whose length differs
    # from the original read's `quality` vector, so the `quality[pos:pos+k-1]` slices
    # below would read out of bounds. Reconcile the quality to the sequence length
    # (no-op on the substitution-only path where they already match → byte-identical).
    quality = _reconcile_quality_length(quality, length(sequence))

    total_log_likelihood = 0.0

    # Use Kmers.jl iterators for efficient k-mer extraction with quality awareness
    if sequence isa BioSequences.LongDNA
        for (pos, kmer) in enumerate(Kmers.FwDNAMers{k}(sequence))
            resolved_kmer = _resolve_kmer_label(graph, kmer; graph_mode = graph_mode)

            if resolved_kmer !== nothing
                vertex_data = graph[resolved_kmer]

                # Enhanced quality-aware likelihood calculation
                kmer_likelihood = calculate_qualmer_likelihood(
                    resolved_kmer,
                    quality[pos:(pos + k - 1)],  # Quality scores for this k-mer
                    vertex_data
                )

                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                # Penalty for unseen k-mers, but consider quality
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:(pos + k - 1)])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    elseif sequence isa BioSequences.LongRNA
        for (pos, kmer) in enumerate(Kmers.FwRNAMers{k}(sequence))
            resolved_kmer = _resolve_kmer_label(graph, kmer; graph_mode = graph_mode)

            if resolved_kmer !== nothing
                vertex_data = graph[resolved_kmer]

                kmer_likelihood = calculate_qualmer_likelihood(
                    resolved_kmer,
                    quality[pos:(pos + k - 1)],
                    vertex_data
                )

                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:(pos + k - 1)])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    else  # LongAA - no canonical form needed
        for (pos, kmer) in enumerate(Kmers.FwAAMers{k}(sequence))
            if haskey(graph, kmer)
                vertex_data = graph[kmer]

                kmer_likelihood = calculate_qualmer_likelihood(
                    kmer,
                    quality[pos:(pos + k - 1)],
                    vertex_data
                )

                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:(pos + k - 1)])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    end

    return total_log_likelihood
end

function _rhizomorph_first_dataset_id(vertex_data)
    dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
    return isempty(dataset_ids) ? nothing : first(dataset_ids)
end

function _rhizomorph_joint_probability(vertex_data)::Float64
    dataset_id = _rhizomorph_first_dataset_id(vertex_data)
    if dataset_id === nothing
        return 0.0
    end

    joint_quality = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, dataset_id)
    if joint_quality === nothing
        return 0.0
    end

    log_prob = 0.0
    for q in joint_quality
        log_prob += log(max(phred_to_probability(q), 1e-10))
    end
    return exp(log_prob)
end

function _vertex_joint_probability(vertex_data)::Float64
    if hasproperty(vertex_data, :joint_probability)
        return vertex_data.joint_probability
    end
    return _rhizomorph_joint_probability(vertex_data)
end

function _vertex_mean_quality(vertex_data)::Float64
    if hasproperty(vertex_data, :mean_quality)
        return vertex_data.mean_quality
    end

    dataset_id = _rhizomorph_first_dataset_id(vertex_data)
    if dataset_id === nothing
        return 0.0
    end

    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, dataset_id)
    if mean_quality === nothing || isempty(mean_quality)
        return 0.0
    end

    return Statistics.mean(mean_quality)
end

function _vertex_coverage(vertex_data)::Int
    if hasproperty(vertex_data, :coverage)
        coverage = vertex_data.coverage
        return coverage isa Integer ? coverage : length(coverage)
    end
    return Int(Mycelia.Rhizomorph.count_evidence(vertex_data))
end

"""
Calculate likelihood of a k-mer given its observed quality scores and qualmer vertex data.
Leverages existing qualmer graph quality calculations.
"""
function calculate_qualmer_likelihood(kmer, quality_scores::AbstractVector{Int8}, vertex_data)::Float64
    # Base likelihood from graph joint probability
    base_likelihood = _vertex_joint_probability(vertex_data)

    # Quality-based confidence adjustment
    # Convert Int8 quality scores to UInt8 for compatibility with existing functions
    quality_uint8 = UInt8.(max.(0, quality_scores))

    # Calculate position-wise quality confidence
    quality_confidence = 1.0
    for q_score in quality_uint8
        prob_correct = phred_to_probability(q_score)
        quality_confidence *= prob_correct
    end

    # Combine graph confidence with observed quality confidence
    # Weight the graph probability by the quality confidence of this observation
    combined_likelihood = base_likelihood * quality_confidence +
                          (1.0 - quality_confidence) * (base_likelihood * 0.1)  # Lower confidence for low quality

    return max(combined_likelihood, 1e-10)
end

"""
Calculate penalty for unseen k-mers based on their quality scores.
High quality unseen k-mers get less penalty than low quality ones.
"""
function calculate_unseen_kmer_penalty(quality_scores::AbstractVector{Int8})::Float64
    quality_uint8 = UInt8.(max.(0, quality_scores))

    # Calculate average quality confidence
    avg_quality_confidence = Statistics.mean(phred_to_probability(q) for q in quality_uint8)

    # High quality unseen k-mers might be real but missing from training data
    # Low quality unseen k-mers are likely errors
    base_penalty = 1e-6
    quality_adjusted_penalty = base_penalty * avg_quality_confidence +
                               base_penalty * 0.01 * (1.0 - avg_quality_confidence)

    return quality_adjusted_penalty
end

"""
Generate alternative k-mers for improvement attempts using proper k-mer objects.
"""
function generate_kmer_alternatives(original_kmer::Kmers.Kmer, graph; graph_mode::Symbol = :canonical)
    sequence_type = if original_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif original_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else
        BioSequences.LongAA
    end
    return generate_kmer_alternatives(original_kmer, graph, sequence_type; graph_mode = graph_mode)
end

function generate_kmer_alternatives(original_kmer::AbstractString, graph; graph_mode::Symbol = :canonical)
    alphabet = detect_alphabet(original_kmer)
    sequence_type = alphabet_to_biosequence_type(alphabet)
    k = length(original_kmer)
    kmer_type = if sequence_type <: BioSequences.LongDNA
        Kmers.DNAKmer{k}
    elseif sequence_type <: BioSequences.LongRNA
        Kmers.RNAKmer{k}
    else
        Kmers.AAKmer{k}
    end
    return generate_kmer_alternatives(kmer_type(original_kmer), graph, sequence_type; graph_mode = graph_mode)
end

function generate_kmer_alternatives(
        original_kmer, graph, sequence_type::Type{<:BioSequences.BioSequence};
        graph_mode::Symbol = :canonical)
    alternatives = []
    k = length(original_kmer)

    # Get alphabet based on sequence type
    if sequence_type <: BioSequences.LongDNA
        alphabet = [
            BioSequences.DNA_A, BioSequences.DNA_T, BioSequences.DNA_G, BioSequences.DNA_C]
        kmer_type = Kmers.DNAKmer{k}
    elseif sequence_type <: BioSequences.LongRNA
        alphabet = [
            BioSequences.RNA_A, BioSequences.RNA_U, BioSequences.RNA_G, BioSequences.RNA_C]
        kmer_type = Kmers.RNAKmer{k}
    else  # LongAA
        # For amino acids, use the standard 20 amino acids
        alphabet = [
            BioSequences.AA_A, BioSequences.AA_R, BioSequences.AA_N, BioSequences.AA_D,
            BioSequences.AA_C, BioSequences.AA_Q, BioSequences.AA_E, BioSequences.AA_G,
            BioSequences.AA_H, BioSequences.AA_I, BioSequences.AA_L, BioSequences.AA_K,
            BioSequences.AA_M, BioSequences.AA_F, BioSequences.AA_P, BioSequences.AA_S,
            BioSequences.AA_T, BioSequences.AA_W, BioSequences.AA_Y, BioSequences.AA_V]
        kmer_type = Kmers.AAKmer{k}
    end

    # Try single position substitutions
    for i in 1:k
        original_symbol = original_kmer[i]
        for symbol in alphabet
            if symbol != original_symbol
                # Build alternative k-mer with substitution at position i
                alt_kmer_seq = sequence_type([original_kmer[j] for j in 1:k])
                alt_kmer_seq[i] = symbol
                alt_kmer = kmer_type(alt_kmer_seq)

                # Resolve k-mer labels against the current graph mode
                resolved_alt = _resolve_kmer_label(graph, alt_kmer; graph_mode = graph_mode)
                resolved_orig = _resolve_kmer_label(graph, original_kmer; graph_mode = graph_mode)

                # Only include if alternative exists in graph with higher probability
                if resolved_alt !== nothing && resolved_orig !== nothing
                    alt_data = graph[resolved_alt]
                    orig_data = graph[resolved_orig]

                    if _vertex_joint_probability(alt_data) >
                       _vertex_joint_probability(orig_data)
                        push!(alternatives, resolved_alt)
                    end
                end
            end
        end
    end

    return alternatives
end

"""
Generate alternative k-mer paths by substituting individual k-mers with plausible alternatives.
Returns a list of alternative paths, each represented as a vector.
"""
function generate_alternative_paths(sequence::AbstractString, graph, k::Int;
        num_samples::Int = 5, graph_mode::Symbol = :canonical)
    base_path = sequence_to_kmer_path(String(sequence), k)
    if isempty(base_path) || num_samples <= 0
        return Vector{Vector{Any}}()
    end

    alternative_paths = Vector{Vector{Any}}()
    for (idx, kmer) in enumerate(base_path)
        alternatives = generate_kmer_alternatives(kmer, graph; graph_mode = graph_mode)
        for alt in alternatives
            if alt != kmer
                alt_path = Any[k for k in base_path]
                alt_path[idx] = alt
                push!(alternative_paths, alt_path)
                if length(alternative_paths) >= num_samples
                    return alternative_paths
                end
                break
            end
        end
    end

    return alternative_paths
end

"""
Adjust quality scores based on likelihood improvement.
"""
function adjust_quality_scores(
        original_quality::AbstractString, improved_sequence::AbstractString,
        likelihood_improvement::Float64)::String
    original_quality_str = String(original_quality)
    improved_sequence_str = String(improved_sequence)

    # For now, return original quality scores
    # This can be enhanced to adjust qualities based on the improvement
    if length(improved_sequence_str) == length(original_quality_str)
        return original_quality_str
    else
        # Handle length changes by extending or truncating quality
        if length(improved_sequence_str) > length(original_quality_str)
            # Extend with median quality
            median_qual = isempty(original_quality_str) ? 'I' :
                          original_quality_str[length(original_quality_str) ÷ 2]
            return original_quality_str *
                   repeat(string(median_qual), length(improved_sequence_str) -
                                               length(original_quality_str))
        else
            # Truncate quality
            return original_quality_str[1:length(improved_sequence_str)]
        end
    end
end

"""
    _quality_from_indel_trace(original_quality, move_trace, read_index_trace, k,
        corrected_length)

Map observed FASTQ qualities onto a pair-HMM corrected sequence from the exact
decoder traceback. The first matched k-mer inherits qualities `1:k`; subsequent
`M` moves inherit the newly consumed observed base, `I` moves drop the consumed
extra base, and `D` moves emit a corrected base with the lower of its adjacent
observed qualities. Returns `nothing` when the trace is malformed or does not
produce `corrected_length`, so callers can fail closed rather than attach shifted
qualities. Runtime and memory are linear in the traceback length.
"""
function _quality_from_indel_trace(
        original_quality::AbstractString,
        move_trace::AbstractVector{Symbol},
        read_index_trace::AbstractVector{Int},
        k::Int,
        corrected_length::Int
)::Union{String, Nothing}
    quality_chars = collect(original_quality)
    n_quality = length(quality_chars)
    if k <= 0 || n_quality < k || isempty(move_trace) ||
       length(move_trace) != length(read_index_trace) || first(move_trace) != :M ||
       first(read_index_trace) != 1
        return nothing
    end

    corrected_quality = copy(quality_chars[1:k])
    n_observations = n_quality - k + 1
    previous_read_index = 1
    @inbounds for trace_index in 2:length(move_trace)
        phase = move_trace[trace_index]
        read_index = read_index_trace[trace_index]
        expected_read_index = phase == :D ? previous_read_index : previous_read_index + 1
        if phase ∉ (:M, :I, :D) || read_index != expected_read_index ||
           !(1 <= read_index <= n_observations)
            return nothing
        end
        observed_position = read_index + k - 1
        if phase == :M
            push!(corrected_quality, quality_chars[observed_position])
        elseif phase == :I
            # The read consumed an extra observed base without advancing the graph;
            # the corrected latent sequence therefore emits no base or quality.
            previous_read_index = read_index
            continue
        elseif phase == :D
            left_position = observed_position
            right_position = min(observed_position + 1, n_quality)
            push!(corrected_quality,
                min(quality_chars[left_position], quality_chars[right_position]))
        end
        previous_read_index = read_index
    end
    length(corrected_quality) == corrected_length || return nothing
    return String(corrected_quality)
end

# =============================================================================
# Decision Making and Termination Conditions
# =============================================================================

"""
Determine if sufficient improvements were made to continue with current k.
Enhanced with convergence detection and adaptive thresholds.
"""
function sufficient_improvements(
        improvements_made::Int, total_reads::Int, threshold::Float64 = 0.05;
        iteration_history::Vector{Dict{Symbol, Any}} = Dict{Symbol, Any}[],
        min_improvement_trend::Int = 3)::Bool
    if total_reads == 0
        return false
    end

    improvement_rate = improvements_made / total_reads

    # Basic threshold check
    if improvement_rate < threshold
        return false
    end

    # Early convergence detection based on improvement trend
    if length(iteration_history) >= min_improvement_trend
        recent_rates = [hist[:improvement_rate]
                        for hist in
                            iteration_history[(end - min_improvement_trend + 1):end]]

        # Check for diminishing returns (decreasing trend)
        if length(recent_rates) >= 2
            trend_decreasing = all(recent_rates[i] >= recent_rates[i + 1]
            for i in 1:(length(recent_rates) - 1))
            avg_recent_rate = Statistics.mean(recent_rates)

            # If trend is consistently decreasing and below adaptive threshold, stop
            if trend_decreasing && avg_recent_rate < threshold * 2.0
                return false
            end
        end

        # Check for plateau (very small improvements)
        if all(rate < threshold * 0.5 for rate in recent_rates)
            return false
        end
    end

    return true
end

"""
Enhanced convergence detection for k-mer progression.
Determines if we should move to the next k-mer size based on multiple criteria.
"""
function should_continue_k_progression(k_history::Vector{Dict{Symbol, Any}},
        current_k::Int,
        max_k::Int;
        convergence_window::Int = 3,
        quality_improvement_threshold::Float64 = 0.001)::Bool
    if current_k >= max_k
        return false
    end

    if length(k_history) < convergence_window
        return true
    end

    # Check recent quality improvements across iterations
    recent_history = k_history[(end - convergence_window + 1):end]

    # Calculate average improvement metrics
    avg_improvement_rate = Statistics.mean([hist[:total_improvements] / hist[:total_reads]
                                            for hist in recent_history])
    avg_runtime_per_read = Statistics.mean([hist[:runtime_seconds] / hist[:total_reads]
                                            for hist in recent_history])

    # Convergence criteria
    improvements_declining = avg_improvement_rate < quality_improvement_threshold
    runtime_increasing = avg_runtime_per_read > 0.1  # More than 0.1 seconds per read indicates inefficiency

    # If improvements are minimal and runtime is high, move to next k
    if improvements_declining && runtime_increasing
        return false
    end

    # If we're making good progress, continue
    return avg_improvement_rate >= quality_improvement_threshold
end

"""
Finalize iterative assembly by combining results from all k-mer sizes and iterations.
"""
function finalize_iterative_assembly(output_dir::String, k_progression::Vector{Int},
        iteration_history::Dict{Int, Vector{Dict{Symbol, Any}}},
        total_runtime::Float64; verbose::Bool = true,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing,
        skip_fractions::Vector{Float64} = Float64[],
        cheap_correction_counts::Vector{Int} = Int[],
        hard_window::Bool = false,
        windowed_decode::Bool = false,
        soft_em::Bool = false,
        cheap_correct::Bool = false,
        min_decode_k::Union{Int, Nothing} = nothing,
        decode_gate_density::Union{Float64, Nothing} = nothing,
        decode_gated_rungs::Vector{Int} = Int[],
        indel_schedule::Symbol = :unrestricted,
        indel_rung_telemetry::Vector{Dict{Symbol, Any}} = Dict{Symbol, Any}[],
        final_pass_graph = nothing,
        final_pass_graph_k::Int = 0,
        final_pass_graph_mode::Symbol = :canonical,
        final_pass_graph_reusable::Bool = false)
    if verbose
        println("Finalizing iterative assembly results...")
    end

    # Only retain the (potentially large) final-pass graph when it is byte-identical-
    # reusable (td-04tb); otherwise drop the reference so a fallback rebuild does not
    # double-hold graph memory.
    if !final_pass_graph_reusable
        final_pass_graph = nothing
    end

    # Find the final output file (last iteration of largest k)
    final_k = maximum(k_progression)
    final_iteration_count = length(iteration_history[final_k])

    # Get the most recent timestamp from final k
    final_timestamp = iteration_history[final_k][end][:timestamp]
    final_fastq = joinpath(
        output_dir, "reads_k$(final_k)_iter$(final_iteration_count)_$(final_timestamp).fastq")

    # Calculate summary statistics
    total_iterations = sum(length(iterations) for iterations in values(iteration_history))
    total_improvements = sum(sum(iter[:improvements_made] for iter in iterations)
    for iterations in values(iteration_history))

    # Read final assembly for k-mer extraction
    if isfile(final_fastq)
        final_reads = collect(FASTX.FASTQ.Reader(open(final_fastq)))
        final_assembly = [FASTX.sequence(String, read) for read in final_reads]
    else
        final_assembly = String[]
    end

    # Swallowed-decode tally (structural / un-k-merizable). Surfaced on the result
    # so a caller can detect a corrector that reported 0 improvements because it
    # never actually ran, distinct from "ran and found nothing to fix".
    corrector_errors = diagnostics === nothing ?
                       Dict(:structural => 0, :unkmerizable => 0,
        :gate_skipped => 0,
        :indel_attempts => 0, :indel_decodes => 0, :truncated_decodes => 0,
        :indel_engaged => 0,
        :trace_contract_errors => 0, :window_divergences => 0) :
                       Dict(:structural => diagnostics.structural_errors[],
        :unkmerizable => diagnostics.unkmerizable_reads[],
        :gate_skipped => diagnostics.gate_skipped[],
        :indel_attempts => diagnostics.indel_attempts[],
        :indel_decodes => diagnostics.indel_decodes[],
        :truncated_decodes => diagnostics.truncated_decodes[],
        :indel_engaged => diagnostics.indel_engaged[],
        :trace_contract_errors => diagnostics.trace_contract_errors[],
        :window_divergences => diagnostics.window_divergences[])

    indel_requested = sum(get(rung, :requested, 0) for rung in indel_rung_telemetry)
    indel_attempted = sum(get(rung, :attempted, 0) for rung in indel_rung_telemetry)
    indel_completed = sum(get(rung, :completed, 0) for rung in indel_rung_telemetry)
    indel_truncated = sum(get(rung, :truncated, 0) for rung in indel_rung_telemetry)
    indel_engaged = sum(get(rung, :engaged, 0) for rung in indel_rung_telemetry)

    # Per-pass skip-fraction summary (FIX 6): min/mean/max across every k+iteration
    # pass, not just the last. `last_skip_fraction` is retained for back-compat.
    last_skip_fraction = isempty(skip_fractions) ? 0.0 : skip_fractions[end]
    skip_min = isempty(skip_fractions) ? 0.0 : minimum(skip_fractions)
    skip_max = isempty(skip_fractions) ? 0.0 : maximum(skip_fractions)
    skip_mean = isempty(skip_fractions) ? 0.0 : sum(skip_fractions) / length(skip_fractions)

    # Graph-Viterbi decode fraction per pass = 1 - skip_fraction (the reads that
    # actually reached the expensive decode). This is the critical-path metric:
    # narrowing the hard-vertex set + Stage 0 cheap correction should drive it down
    # toward the true bubble/repeat ~5-15% (td-bjnt). Zero on tiers without a gate.
    decode_fractions = [1.0 - s for s in skip_fractions]
    decode_min = isempty(decode_fractions) ? 0.0 : minimum(decode_fractions)
    decode_max = isempty(decode_fractions) ? 0.0 : maximum(decode_fractions)
    decode_mean = isempty(decode_fractions) ? 0.0 :
                  sum(decode_fractions) / length(decode_fractions)

    # Stage 0 cheap-correction totals (td-bjnt).
    total_cheap_corrections = isempty(cheap_correction_counts) ? 0 :
                              sum(cheap_correction_counts)

    # Create comprehensive metadata
    metadata = Dict(
        :total_runtime => total_runtime,
        :total_iterations => total_iterations,
        :total_improvements => total_improvements,
        :k_progression => k_progression,
        :final_k => final_k,
        :final_fastq_file => final_fastq,
        :output_directory => output_dir,
        :iteration_history => iteration_history,
        :corrector_errors => corrector_errors,
        :assembly_type => "iterative_maximum_likelihood",
        :version => "Phase_5.2a",
        # Scalable-tier telemetry (td-fuo8/td-nn6l/td-e70t). Zero/false on the
        # exhaustive tier, which threads neither gate.
        #
        # Honest gate provenance (FIX 5): `hard_read_gate` is the SKIP gate
        # (easy reads pass through untouched) — that IS active on :scalable.
        # `windowed_decode` (td-nn6l Stage 3c) is the per-hard-region WINDOWED
        # decode: when `true`, a hard read is decoded window-by-window (only the
        # boundary-constrained hard sub-window(s), <=500 bp) instead of whole-read.
        # Surfaced separately from `hard_window` (the SKIP gate) so the two are not
        # conflated. `hard_window` is kept as a back-compat alias for the skip gate.
        :hard_window => hard_window,
        :hard_read_gate => hard_window,
        :windowed_decode => windowed_decode,
        # Soft-EM v2 ACTIVE: the E-step enumerates competing paths and the M-step
        # registers the support-floored soft edge weights onto the next EM
        # iteration's graph, so unsupported error edges decay while supported
        # variation is retained. `false` on the exhaustive tier (soft-EM off there).
        :soft_em => (soft_em ? "v2-competing-paths-floor" : false),
        # Low-k decode gating (td-9h5r): the explicit floor (if any), the adaptive
        # density threshold, and the specific rungs whose per-read decode was gated
        # OFF (explicit floor OR adaptive). `nothing`/empty on the exhaustive tier.
        :min_decode_k => min_decode_k,
        :decode_gate_density => decode_gate_density,
        :decode_gated_rungs => decode_gated_rungs,
        # Branching/frontier-aware sparse-rung scheduling. The vector retains one
        # entry per actual pass, including repeated iterations at the same k.
        :indel_schedule => indel_schedule,
        :indel_rung_telemetry => indel_rung_telemetry,
        :indel_requested => indel_requested,
        :indel_attempted => indel_attempted,
        :indel_completed => indel_completed,
        :indel_truncated => indel_truncated,
        :indel_engaged => indel_engaged,
        :last_skip_fraction => last_skip_fraction,
        :skip_fraction_per_pass => skip_fractions,
        :skip_fraction_min => skip_min,
        :skip_fraction_mean => skip_mean,
        :skip_fraction_max => skip_max,
        # Graph-Viterbi decode fraction (= 1 - skip_fraction). The critical-path
        # win metric (td-bjnt): after Stage 0 + hard-set narrowing this drops
        # toward the true bubble/repeat ~5-15%.
        :decode_fraction_per_pass => decode_fractions,
        :decode_fraction_min => decode_min,
        :decode_fraction_mean => decode_mean,
        :decode_fraction_max => decode_max,
        # Stage 0 cheap k-mer-spectrum correction (td-bjnt): bases fixed by the
        # linear pass, per pass and in total. Zero/false on :exhaustive.
        :cheap_correct => cheap_correct,
        :cheap_corrections_per_pass => cheap_correction_counts,
        :cheap_corrections_total => total_cheap_corrections,
        # Final-pass graph reuse provenance (td-04tb). `final_graph_reusable` is true
        # only when the final pass at the largest k made 0 improvements, so the
        # returned `:final_graph` is byte-identical to a from-scratch rebuild of the
        # corrected reads at `final_graph_k` under `final_graph_mode`. The caller
        # verifies k + mode match its re-assembly parameters before reusing.
        :final_graph_k => final_pass_graph_k,
        :final_graph_mode => final_pass_graph_mode,
        :final_graph_reusable => final_pass_graph_reusable
    )

    if verbose
        println("Iterative assembly finalized:")
        println("  Total k-mer sizes: $(length(k_progression))")
        println("  Total iterations: $total_iterations")
        println("  Total improvements: $total_improvements")
        println("  Final FASTQ: $final_fastq")
    end

    return Dict(
        :final_assembly => final_assembly,
        :k_progression => k_progression,
        :metadata => metadata,
        # The corrector's final-pass qualmer graph (td-04tb). `nothing` when no pass
        # completed; only byte-identical-reusable when `metadata[:final_graph_reusable]`.
        :final_graph => final_pass_graph
    )
end

# =============================================================================
# Utility and Integration Functions
# =============================================================================

"""
Generate summary report of iterative assembly process.
"""
function iterative_assembly_summary(result::Dict)::String
    metadata = result[:metadata]
    k_prog = metadata[:k_progression]

    report = "Mycelia Iterative Assembly Summary\n"
    report *= "=" * "="^50 * "\n\n"

    report *= "Assembly Type: $(metadata[:assembly_type])\n"
    report *= "Version: $(metadata[:version])\n"
    report *= "Total Runtime: $(round(metadata[:total_runtime], digits=2)) seconds\n"
    report *= "K-mer Progression: $(sort(k_prog))\n"
    report *= "Total Iterations: $(metadata[:total_iterations])\n"
    report *= "Total Improvements: $(metadata[:total_improvements])\n"
    report *= "Final FASTQ: $(metadata[:final_fastq_file])\n\n"

    # Per-k statistics
    report *= "Per-K Statistics:\n"
    for k in sort(k_prog)
        iterations = metadata[:iteration_history][k]
        total_improvements_k = sum(iter[:improvements_made] for iter in iterations)
        report *= "  k=$k: $(length(iterations)) iterations, $total_improvements_k improvements\n"
    end

    return report
end

"""
Run a minimal self-test of the iterative assembly pipeline.
Returns a status dictionary with results or error details.
"""
function test_iterative_assembly(; max_k::Int = 13, max_iterations_per_k::Int = 1, cleanup::Bool = true)::Dict
    temp_dir = mktempdir()
    try
        test_fastq = joinpath(temp_dir, "test_reads.fastq")
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATC",
            "ATCGATCGATCGTTCGATCG"
        ]
        records = [FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
                   for (i, seq) in enumerate(sequences)]
        write_fastq(records = records, filename = test_fastq)

        output_dir = joinpath(temp_dir, "test_output")
        result = mycelia_iterative_assemble(
            test_fastq;
            max_k = max_k,
            output_dir = output_dir,
            max_iterations_per_k = max_iterations_per_k,
            verbose = false
        )

        return Dict(
            :status => :success,
            :final_assembly => result[:final_assembly],
            :k_progression => result[:k_progression],
            :metadata => result[:metadata]
        )
    catch e
        return Dict(:status => :error, :error => sprint(showerror, e))
    finally
        if cleanup
            try
                rm(temp_dir, recursive = true, force = true)
            catch
            end
        end
    end
end

# =============================================================================
# Enhanced Statistical Path Improvement (Phase 5.2b)
# Integration with Viterbi algorithms and statistical path resampling
# =============================================================================

"""
Try to improve FASTQ read using Viterbi algorithm from viterbi-next.jl.
Returns (improved_read, likelihood) or nothing if no improvement.
"""
function try_viterbi_path_improvement(read::FASTX.FASTQ.Record,
        graph,
        k::Int;
        graph_mode::Symbol = :canonical,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        weighted_graph = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing,
        indel_params::Union{Nothing, IndelDecodeParams} = nothing,
        calibrated_gap_threshold::Union{Float64, Nothing} = nothing)::Union{
        Tuple{FASTX.FASTQ.Record, Float64}, Nothing}
    try
        sequence_string = FASTX.sequence(String, read)
        alphabet = Mycelia.detect_alphabet(sequence_string)
        sequence_type = Mycelia.alphabet_to_biosequence_type(alphabet)
        sequence = Mycelia.extract_typed_sequence(read, sequence_type)
        # NOTE: `_record_kmer_iterator` uses a 2-bit (unambiguous) k-mer iterator,
        # so `collect` here throws `BioSequences.EncodeError` on a read carrying an
        # ambiguous base (e.g. a single `N`). That is handled as a SKIP-with-count
        # in the catch below, NOT a crash of the (threaded) batch (td-63qy family).
        kmers = collect(Mycelia._record_kmer_iterator(sequence_type, k, sequence))
        if isempty(kmers)
            return nothing
        end

        # Restore the emission model (graph-as-HMM correction core, td-jqdf): the
        # read's per-base Phred quality is the emission P(observed | hidden). The
        # k-mer at index `i` spans read positions `i:(i+k-1)`, so wrap it with that
        # per-base window (clamped to bounds) instead of discarding the quality and
        # letting the emission fall back to the graph's population-average quality.
        quality_scores = collect(FASTX.quality_scores(read))
        n_qual = length(quality_scores)
        observations = Vector{Mycelia.QualityObservation}(undef, length(kmers))
        for (i, kmer) in enumerate(kmers)
            lo = clamp(i, 1, n_qual)
            hi = clamp(i + k - 1, 1, n_qual)
            window = UInt8.(@view quality_scores[lo:hi])
            observations[i] = Mycelia.QualityObservation(kmer, window)
        end

        # Trustworthy Viterbi (td-ak6w), reconciled with the td-63qy OOM crash and
        # the td-35ux short-read O(genome^2) frontier blowup: the decoder is EXACT
        # where tractable and BOUNDED (with a log line) otherwise. `beam_width ===
        # nothing` (the default) resolves the frontier bound PER READ via the
        # GRAPH-DENSITY-AWARE `_auto_beam_width(length(observations), Graphs.nv(graph))`:
        # `typemax(Int)` (never pruned → global maximum-likelihood path) only when
        # BOTH the read is small AND the graph is sparse; the proven-tractable
        # finite bound when the read is large (read-length rule, td-63qy) OR the
        # graph is dense (density rule, td-35ux — a short read on a big graph
        # otherwise retains an O(n_vertices) frontier at each decode depth). The
        # graph threaded into this call is the one being decoded, so its vertex
        # count is exactly the density the frontier can reach. An explicit
        # `beam_width::Int` (including `typemax(Int)`) is an OPT-IN override honored
        # verbatim — a caller/test can force exact decoding on every read. This
        # trades the ML guarantee for tractability only under the auto-default; it
        # never silently changes correctness of an explicit request.
        effective_beam_width = beam_width === nothing ?
                               _auto_beam_width(length(observations), Graphs.nv(graph)) :
                               beam_width
        # Candidate-generation bounds (td-plqi) engage ONLY where the width beam is
        # already finite (i.e. the decode is already the bounded approximation, not
        # exact ML). When the beam is exact (typemax — small read on a sparse graph,
        # OR an explicit caller override forcing exact), the score margin stays Inf
        # and the successor bound stays unbounded, so exact-ML reads are byte-
        # identical. Where the beam is bounded (dense/large — the k=9-style rungs
        # that carry the #386 super-linear residual) the margin holds the generating
        # frontier O(1) in genome size.
        beam_is_exact = effective_beam_width == typemax(Int)
        effective_margin = beam_is_exact ? Inf : _AUTO_BEAM_SCORE_MARGIN
        effective_successor_bound = beam_is_exact ? typemax(Int) : _AUTO_SUCCESSOR_BOUND
        # Indel-aware wiring (td-9q84 / 4a): when the assembler's error profile
        # enables indels it threads a non-nothing `indel_params`; the config then
        # turns on the pair-HMM gap moves with the profile's gap-open fractions +
        # extend probabilities and the tier's run caps + band. Crucially it also sets
        # `error_rate` to the profile's ABSOLUTE per-base rate so the gap-open masses
        # (`δ_I = error_rate·f_ins`, `δ_D = error_rate·f_del`) are scaled correctly —
        # leaving it at the 0.01 default would under-weight nanopore/pacbio gaps ~10×
        # (e.g. 0.01×0.30 = 0.003 instead of 0.10×0.30 = 0.03). `indel_params ===
        # nothing` (the substitution-only default, e.g. :illumina) builds the config
        # EXACTLY as the pre-wiring corrector did — `error_rate` stays at its 0.01
        # default and `indel_moves` stays false, so the decode is byte-identical
        # (oracle preservation).
        config = if indel_params === nothing
            Mycelia.ViterbiCorrectionConfig(
                alphabet = alphabet,
                strand_mode = graph_mode,
                max_steps = length(observations) - 1,
                beam_width = effective_beam_width,
                max_successors_per_state = effective_successor_bound,
                beam_score_margin = effective_margin,
                record_position_gaps = calibrated_gap_threshold !== nothing
            )
        else
            Mycelia.ViterbiCorrectionConfig(
                alphabet = alphabet,
                strand_mode = graph_mode,
                max_steps = length(observations) - 1,
                beam_width = effective_beam_width,
                max_successors_per_state = effective_successor_bound,
                beam_score_margin = effective_margin,
                error_rate = indel_params.base_error_rate,
                indel_moves = true,
                insertion_fraction = indel_params.insertion_fraction,
                deletion_fraction = indel_params.deletion_fraction,
                insertion_extend_probability = indel_params.insertion_extend_probability,
                deletion_extend_probability = indel_params.deletion_extend_probability,
                deletion_max_run = indel_params.deletion_max_run,
                max_insertion_run = indel_params.max_insertion_run,
                band_width = indel_params.band_width,
                record_position_gaps = calibrated_gap_threshold !== nothing
            )
        end
        if indel_params !== nothing && diagnostics !== nothing
            Threads.atomic_add!(diagnostics.indel_attempts, 1)
        end
        correction = Mycelia.correct_observations(
            graph, [observations]; config = config, weighted_graph = weighted_graph)
        corrected_path = only(correction.corrected_observations)
        if corrected_path === nothing || isempty(corrected_path)
            return nothing
        end
        path_diagnostics = only(correction.paths).diagnostics
        # The substitution kernel is length-preserving by contract: one decoded
        # graph unit per observation. Reject a partial/dead-end path before sequence
        # reconstruction, likelihood comparison, or soft-EM accumulation. Windowed
        # callers then pass through the original bytes rather than recording a
        # post-hoc splice divergence.
        if indel_params === nothing && length(corrected_path) != length(observations)
            diagnostics === nothing ||
                Threads.atomic_add!(diagnostics.structural_errors, 1)
            return nothing
        end
        corrected_sequence = Mycelia.Rhizomorph.path_to_sequence(corrected_path, graph)
        corrected_sequence_string = corrected_sequence isa AbstractString ?
                                    corrected_sequence :
                                    string(corrected_sequence)
        aligned_indel_quality::Union{Nothing, String} = nothing
        if indel_params !== nothing
            if get(path_diagnostics, :algorithm, nothing) != :viterbi_indel_pair_hmm
                diagnostics === nothing ||
                    Threads.atomic_add!(diagnostics.trace_contract_errors, 1)
                return nothing
            end
            if get(path_diagnostics, :truncated, false) === true
                diagnostics === nothing ||
                    Threads.atomic_add!(diagnostics.truncated_decodes, 1)
                return nothing
            end
            aligned_indel_quality = Mycelia._quality_from_indel_trace(
                String(FASTX.quality(read)),
                get(path_diagnostics, :move_trace, Symbol[]),
                get(path_diagnostics, :read_index_trace, Int[]),
                k,
                length(corrected_sequence_string))
            if aligned_indel_quality === nothing
                diagnostics === nothing ||
                    Threads.atomic_add!(diagnostics.trace_contract_errors, 1)
                return nothing
            end
            if diagnostics !== nothing
                Threads.atomic_add!(diagnostics.indel_decodes, 1)
                move_counts = get(
                    path_diagnostics, :move_counts,
                    Dict{Symbol, Int}(:M => 0, :I => 0, :D => 0))
                if get(move_counts, :I, 0) + get(move_counts, :D, 0) > 0
                    Threads.atomic_add!(diagnostics.indel_engaged, 1)
                end
            end
        end

        # Soft-EM E-step (td-e70t v2): build this read's COMPETING candidate paths
        # (the observed read path + a consensus alternative re-routed through the
        # best-supported sibling branch) and split the responsibility across them
        # by a stable softmax over their normalized-transition log-probabilities,
        # accumulating each path's edges weighted by its share. A data-supported
        # branch keeps high mass; an unsupported (error) branch accrues little and
        # decays across EM iterations (the M-step registers this accumulator —
        # support-floored — onto the next graph). Additive + guarded so it cannot
        # perturb the :exhaustive path (soft_weights === nothing there); degenerates
        # to the single observed path at responsibility 1.0 when no distinct
        # alternative exists (a balanced variant, whose branches are retained
        # through their own reads).
        # A substitution decode with a malformed length is rejected by the
        # window splice contract. Do not let that rejected result influence the
        # next soft-EM M-step; only length-valid substitution results (and
        # trace-valid indel results, checked above) contribute responsibilities.
        soft_result_is_valid = indel_params !== nothing ||
                               length(corrected_sequence_string) ==
                               length(sequence_string)
        if soft_weights !== nothing && soft_result_is_valid
            # Bound the competing-path GENERATION with the same discipline as the
            # decode bounds above (td-e70t speed residual C5c): engage the walk band
            # + successor cap ONLY where the width beam is already finite (the
            # dense/large reads that carry the residual), so exact-ML reads keep the
            # unbounded (byte-identical) generation. The rejoin test precedes the
            # band cutoff, so no real variant's read-consistent competing path drops.
            accumulate_competing_paths!(
                soft_weights, read, graph, k; graph_mode = graph_mode,
                walk_band = beam_is_exact ? typemax(Int) : _soft_em_walk_band(k),
                successor_bound = beam_is_exact ? typemax(Int) :
                                  _SOFT_EM_ALT_SUCCESSOR_BOUND)
        end

        # Calibrated gate is a SUBSTITUTION-ONLY tool: the positional revert
        # `corrected_chars[i+k] = original_chars[i+k]` requires a stable base<->column
        # map, which only a length-preserving substitution decode provides. Under
        # `indel_moves` a BALANCED indel decode (one insertion + one deletion) is
        # net-length-neutral — it would pass a length check yet shift the interior
        # map — so we exclude indel mode entirely rather than trust `length` as a
        # proxy for alignment (review: convergent finding, 3 reviewers).
        if calibrated_gap_threshold !== nothing && indel_params === nothing
            path_result = only(correction.paths)
            decoded_path = path_result.path
            gaps = get(path_result.diagnostics, :position_gaps, nothing)
            original_chars = collect(sequence_string)
            corrected_chars = collect(corrected_sequence_string)
            # gaps[i] is the Viterbi margin INTO path.steps[i+1], i.e. read position
            # i+k. Where the margin is below threshold the decode is ambiguous, so
            # revert that base to the observed one — this is what pulls over-correction
            # back down. `Inf` gaps (collapsed frontier ⇒ maximally confident) clear
            # any finite threshold and are kept.
            if decoded_path !== nothing && gaps !== nothing &&
               length(gaps) == length(decoded_path.steps) - 1 &&
               length(corrected_chars) == length(original_chars)
                @inbounds for i in eachindex(gaps)
                    position = i + k
                    position > length(corrected_chars) && break
                    if gaps[i] < calibrated_gap_threshold
                        corrected_chars[position] = original_chars[position]
                    end
                end
                corrected_sequence_string = String(corrected_chars)
            elseif diagnostics !== nothing
                # In substitution mode the decode contract (path present, gaps
                # recorded, len == steps-1, length preserved) is ALWAYS satisfiable;
                # a miss means the gate silently fell open to the ungated decode —
                # count it so an opt-in gate that disabled itself is visible, not
                # invisible (fail-open is the right data-safety choice; silence is not).
                Threads.atomic_add!(diagnostics.gate_skipped, 1)
            end
        end
        improved_quality, improved_likelihood = if indel_params === nothing
            likelihood = Mycelia.calculate_sequence_likelihood(
                corrected_sequence_string, quality_scores, graph, k;
                graph_mode = graph_mode)
            quality = Mycelia.adjust_quality_scores(
                FASTX.quality(read), corrected_sequence_string, likelihood)
            quality, likelihood
        else
            aligned_quality = aligned_indel_quality::String
            scores = Int8.(Int.(collect(codeunits(aligned_quality))) .- 33)
            likelihood = Mycelia.calculate_sequence_likelihood(
                corrected_sequence_string, scores, graph, k; graph_mode = graph_mode)
            aligned_quality, likelihood
        end
        improved_record = FASTX.FASTQ.Record(
            FASTX.identifier(read), corrected_sequence_string, improved_quality)
        return (improved_record, improved_likelihood)
    catch error
        if error isa ArgumentError
            # STRUCTURAL decode failure (empty graph, alphabet-inference miss,
            # config error). Swallow as "no improvement" so one broken read does
            # not abort the (threaded) batch — but COUNT it so a run that silently
            # corrected nothing is visible (not conflated with "decoder ran, no
            # gain"). See CorrectorDiagnostics + the per-pass @warn.
            diagnostics === nothing ||
                Threads.atomic_add!(diagnostics.structural_errors, 1)
            return nothing
        elseif error isa BioSequences.EncodeError
            # UN-K-MERIZABLE read: the 2-bit k-mer iterator throws EncodeError (NOT
            # ArgumentError) on an ambiguous base (e.g. a single `N` in an Illumina
            # read). Treat as a SKIP — the read passes through uncorrected — and
            # count it, rather than rethrowing and crashing the whole threaded
            # corrector on one N-containing read (td-63qy family).
            diagnostics === nothing ||
                Threads.atomic_add!(diagnostics.unkmerizable_reads, 1)
            return nothing
        end
        rethrow()
    end
end

function try_viterbi_path_improvement(
        sequence::AbstractString, quality, graph, k::Int; graph_mode::Symbol = :canonical)
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_viterbi_path_improvement"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    result = try_viterbi_path_improvement(record, graph, k; graph_mode = graph_mode)
    if result === nothing
        return nothing
    end

    improved_record, improvement = result
    return FASTX.sequence(String, improved_record), improvement
end

"""
Try statistical path resampling for alternative high-likelihood paths.
Returns (improved_read, likelihood) or nothing if no improvement.
"""
function try_statistical_path_resampling(read::FASTX.FASTQ.Record,
        graph,
        k::Int;
        graph_mode::Symbol = :canonical)::Union{Tuple{FASTX.FASTQ.Record, Float64}, Nothing}
    try
        # Detect sequence type dynamically
        sequence_string = FASTX.sequence(String, read)
        alphabet = detect_alphabet(sequence_string)
        sequence_type = alphabet_to_biosequence_type(alphabet)

        sequence = extract_typed_sequence(read, sequence_type)
        quality_scores = collect(FASTX.quality_scores(read))
        quality_string = FASTX.quality(read)

        # Generate multiple alternative paths using qualmer sampling
        alternative_paths = generate_alternative_qualmer_paths(
            read, graph, k, num_samples = 5, graph_mode = graph_mode)

        best_sequence = sequence
        best_likelihood = calculate_sequence_likelihood(
            sequence, quality_scores, graph, k; graph_mode = graph_mode)
        best_quality_scores = quality_scores

        # Evaluate each alternative path
        for alt_path in alternative_paths
            if !isempty(alt_path)
                # Convert qualmer path to proper BioSequence + quality
                alt_result = qualmer_path_to_biosequence(alt_path)
                alt_sequence = alt_result.sequence
                alt_quality_scores = alt_result.quality_scores

                alt_likelihood = calculate_sequence_likelihood(
                    alt_sequence, alt_quality_scores, graph, k; graph_mode = graph_mode)

                if alt_likelihood > best_likelihood
                    best_sequence = alt_sequence
                    best_likelihood = alt_likelihood
                    best_quality_scores = alt_quality_scores
                end
            end
        end

        # Return improvement if found
        if best_sequence != sequence
            # Create improved FASTQ record using proper types
            improved_sequence_string = string(best_sequence)
            improved_quality_string = String([Char(Int(q) + 33)
                                              for q in best_quality_scores])  # Convert back to ASCII
            improved_record = FASTX.FASTQ.Record(
                FASTX.identifier(read), improved_sequence_string, improved_quality_string)
            return (improved_record, best_likelihood)
        end

    catch e
        # If statistical resampling fails, return nothing
        return nothing
    end

    return nothing
end

function try_statistical_path_resampling(
        sequence::AbstractString, quality, graph, k::Int; graph_mode::Symbol = :canonical)
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_statistical_path_resampling"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    result = try_statistical_path_resampling(record, graph, k; graph_mode = graph_mode)
    if result === nothing
        return nothing
    end

    improved_record, improvement = result
    return FASTX.sequence(String, improved_record), improvement
end

"""
Local heuristic improvements (fallback method).
Returns (improved_read, likelihood_improvement).
"""
function try_local_path_improvements(
        read::FASTX.FASTQ.Record, graph, k::Int, original_likelihood::Float64;
        graph_mode::Symbol = :canonical)::Tuple{FASTX.FASTQ.Record, Float64}
    sequence_string = FASTX.sequence(String, read)
    alphabet = detect_alphabet(sequence_string)
    sequence_type = alphabet_to_biosequence_type(alphabet)
    quality_scores = collect(FASTX.quality_scores(read))

    improved_seq = sequence_string  # Start with original
    best_likelihood = original_likelihood

    # Try local improvements at positions with low-probability k-mers
    for i in 1:(length(sequence_string) - k + 1)
        kmer_str = sequence_string[i:(i + k - 1)]

        raw_kmer = if sequence_type <: BioSequences.LongDNA
            Kmers.DNAKmer{k}(kmer_str)
        elseif sequence_type <: BioSequences.LongRNA
            Kmers.RNAKmer{k}(kmer_str)
        else
            Kmers.AAKmer{k}(kmer_str)
        end

        target_kmer = _resolve_kmer_label(graph, raw_kmer; graph_mode = graph_mode)
        if target_kmer !== nothing
            vertex_data = graph[target_kmer]
            if _vertex_joint_probability(vertex_data) < 0.7  # Low confidence k-mer
                # Try improving this position
                alternatives = generate_kmer_alternatives(
                    target_kmer, graph, sequence_type; graph_mode = graph_mode)
                for alt_kmer in alternatives
                    alt_kmer_str = string(alt_kmer)
                    if i == 1
                        candidate_seq = alt_kmer_str * sequence_string[(i + k):end]
                    elseif i+k-1 == length(sequence_string)
                        candidate_seq = sequence_string[1:(i - 1)] * alt_kmer_str
                    else
                        candidate_seq = sequence_string[1:(i - 1)] * alt_kmer_str *
                                        sequence_string[(i + k):end]
                    end

                    candidate_bioseq = detect_and_extract_sequence(candidate_seq, alphabet)
                    candidate_likelihood = calculate_sequence_likelihood(
                        candidate_bioseq, quality_scores, graph, k; graph_mode = graph_mode)

                    if candidate_likelihood > best_likelihood
                        improved_seq = candidate_seq
                        best_likelihood = candidate_likelihood
                    end
                end
            end
        end
    end

    likelihood_improvement = best_likelihood - original_likelihood
    improved_quality = adjust_quality_scores(FASTX.quality(read), improved_seq, likelihood_improvement)
    improved_record = FASTX.FASTQ.Record(FASTX.identifier(read), improved_seq, improved_quality)
    return improved_record, likelihood_improvement
end

function try_local_path_improvements(
        sequence::AbstractString, quality, graph, k::Int, original_likelihood::Float64;
        graph_mode::Symbol = :canonical)::Tuple{String, Float64}
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_local_path_improvements"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    improved_record,
    improvement = try_local_path_improvements(
        record, graph, k, original_likelihood; graph_mode = graph_mode)
    return FASTX.sequence(String, improved_record), improvement
end

# ============================================================================
# Soft-EM E-step + M-step accumulation hook (ACTIVE, td-e70t v2 + support floor)
# ============================================================================
#
# See the design note over `SoftEdgeWeightAccumulator` in
# `src/rhizomorph/core/evidence-functions.jl`. This is the correction-side hook
# that turns a read's competing decoded paths into probability-weighted edge
# evidence: each candidate path's RESPONSIBILITY (posterior over the read's own
# candidate set) is accumulated onto the edges it traverses. Summed over reads,
# the accumulator holds soft edge weights in which error edges — traversed only
# by rare, low-probability paths — accrue little weight; once the M-step registers
# them, `compute_edge_weight` (and thus the Viterbi transition scoring) sees the
# decayed weight and the next iteration's decode + rebuild drops the error edge.
# NOTE the decay acts through `compute_edge_weight` / the Viterbi transition
# score, NOT through the vertex `weight > 0.01` gate in
# `generate_alternative_qualmer_paths` (that gate filters vertex candidates).
#
# v2 status — E-step AND M-step both WIRED. Under `soft_em=true`:
#   * E-step: `accumulate_competing_paths!` (below) builds each decoded read's
#     COMPETING candidate paths — the observed read path and a consensus
#     alternative re-routed through the best-supported sibling branch — then
#     softmax-splits the responsibility across them by normalized-transition
#     log-probability, so a real branch keeps high mass and an unsupported (error)
#     branch gets little.
#   * M-step: `mycelia_iterative_assemble` registers the accumulator onto the NEXT
#     EM iteration's graph (`register_soft_edge_weights!`), so `compute_edge_weight`
#     returns the decayed soft weight and the next iteration's transition scoring
#     is biased away from the error edge — a feedback loop that drives the error
#     edge's weight strictly down.
# VARIATION PRESERVATION (td-h6w9): the M-step registration applies a SUPPORT
# FLOOR — an edge backed by >= `SOFT_EM_MIN_SUPPORT` reads is clamped to at least
# its raw coverage, so a real but SKEWED minority allele (e.g. a 10x branch in a
# 20x/10x bubble) NEVER decays toward zero regardless of a heavier sibling. Only
# near-zero-support (error) edges are free to decay below the cleaning gate. This
# fixes the prior v2's collapse of skewed variants (the responsibility split alone
# gave `W_min' = N*W/(W_maj+W_min)`, a geometric decay to zero for any imbalance).

"""
    _decoded_path_edges(result) -> Vector

Ordered `(src_label, dst_label)` edge tuples of a decoded Viterbi path, or an
empty vector when the result carries no path.
"""
function _decoded_path_edges(result::Mycelia.Rhizomorph.ViterbiDecodingResult)
    labels = _decoded_path_labels(result)
    labels === nothing && return Tuple{Any, Any}[]
    return Tuple{Any, Any}[(labels[i], labels[i + 1]) for i in 1:(length(labels) - 1)]
end

"""
    accumulate_soft_em_edge_weights!(accumulator, candidate_paths) -> accumulator

Soft-EM E-step core (Viterbi-result flavor). Given a single read's competing
decoded paths (`ViterbiDecodingResult`s, each with a `.score` log-likelihood and
a `.path`), compute each path's responsibility by softmax-normalizing its score
against the candidate set (`Mycelia.Rhizomorph.path_responsibility`) and
accumulate that responsibility onto the edges it traverses
(`Mycelia.Rhizomorph.accumulate_path_probability!`).

With a single decoded path per read the responsibility is `1.0` (soft weight ==
hard count); the softening comes from `accumulate_competing_paths!`, which
supplies several candidate paths per read. Retained as a primitive (used by the
unit tests) alongside the graph-level `accumulate_competing_paths!` used by the
pipeline E-step.
"""
function accumulate_soft_em_edge_weights!(
        accumulator::Mycelia.Rhizomorph.SoftEdgeWeightAccumulator,
        candidate_paths
)
    scores = Float64[result.score for result in candidate_paths if isfinite(result.score)]
    isempty(scores) && return accumulator
    for result in candidate_paths
        isfinite(result.score) || continue
        responsibility = Mycelia.Rhizomorph.path_responsibility(result.score, scores)
        Mycelia.Rhizomorph.accumulate_path_probability!(
            accumulator, _decoded_path_edges(result), responsibility)
    end
    return accumulator
end

# ----------------------------------------------------------------------------
# Soft-EM v2 competing-paths E-step (td-e70t): graph-level candidate enumeration
# ----------------------------------------------------------------------------

"""
    _graph_oriented_edges(labels, graph) -> Vector{Tuple}

Turn an ordered list of vertex `labels` into the `(src, dst)` edge tuples in the
ORIENTATION the graph actually stores them (checking both directions), skipping
any consecutive pair with no edge in `graph`. Keying the accumulator in the
graph's own orientation is what lets `register_soft_edge_weights!` — which
iterates `MetaGraphsNext.edge_labels(graph)` — find and consume the soft weight.
A candidate substitution that is not a real graph traversal contributes no edge.
"""
function _graph_oriented_edges(labels, graph)
    edges = Tuple{Any, Any}[]
    for i in 1:(length(labels) - 1)
        a, b = labels[i], labels[i + 1]
        if MetaGraphsNext.haskey(graph, a, b)
            push!(edges, (a, b))
        elseif MetaGraphsNext.haskey(graph, b, a)
            push!(edges, (b, a))
        end
    end
    return edges
end

# The observed read path, resolved to graph vertex labels (canonical where the
# graph is canonical). Empty when no k-mer of the read resolves onto the graph.
function _read_resolved_labels(read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol)
    labels = Any[]
    for (qmer, _pos) in qualmers_unambiguous(read, k)
        resolved = _resolve_qualmer_for_graph(graph, qmer; graph_mode = graph_mode)
        resolved === nothing || push!(labels, resolved.kmer)
    end
    return labels
end

# Soft-weight-aware edge weight between two vertex labels in either stored
# orientation, or `nothing` when they are not adjacent. `compute_edge_weight`
# consults the soft-edge registry, so an edge the M-step decayed reads back its
# decayed weight here — the cross-iteration feedback that sharpens the split.
function _edge_weight_between(graph, a, b)
    if MetaGraphsNext.haskey(graph, a, b)
        return Float64(Mycelia.Rhizomorph.compute_edge_weight(graph[a, b]))
    elseif MetaGraphsNext.haskey(graph, b, a)
        return Float64(Mycelia.Rhizomorph.compute_edge_weight(graph[b, a]))
    end
    return nothing
end

# All labels adjacent to `a` (either orientation), deduplicated. Used both as the
# transition normalizer support and as the branch-candidate set.
function _incident_labels(graph, a)
    ns = Any[]
    for n in MetaGraphsNext.outneighbor_labels(graph, a)
        push!(ns, n)
    end
    for n in MetaGraphsNext.inneighbor_labels(graph, a)
        push!(ns, n)
    end
    return unique(ns)
end

# Sum of the soft-weight-aware edge weights on every edge incident to `a` — the
# denominator that turns a raw edge weight into a normalized transition
# probability. Normalizing here counts a single branch decision ONCE (internal
# edges of a linear run normalize to ~0.5 on BOTH competing branches and cancel
# in the responsibility softmax), so the divergence edge carries the coverage
# contrast and a real but skewed variant is not collapsed like a raw
# coverage-ratio prune would.
function _incident_weight_sum(graph, a)
    total = 0.0
    for x in _incident_labels(graph, a)
        w = _edge_weight_between(graph, a, x)
        w === nothing || (total += w)
    end
    return total
end

# Normalized-transition log-probability of a label path: sum of log(w(a,b)/S(a)).
# Returns -Inf if any consecutive pair is non-adjacent or zero-weight, so an
# invalid candidate is dropped by the softmax `isfinite` guard.
function _path_transition_logscore(labels, graph)
    total = 0.0
    for i in 2:length(labels)
        a, b = labels[i - 1], labels[i]
        w = _edge_weight_between(graph, a, b)
        (w === nothing || w <= 0.0) && return -Inf
        s = _incident_weight_sum(graph, a)
        s <= 0.0 && return -Inf
        total += log(w / s)
    end
    return total
end

# ----------------------------------------------------------------------------
# Competing-path generation bounds (td-e70t speed residual C5c)
# ----------------------------------------------------------------------------
#
# After the decode score-margin fix (#388) linearized per-read Viterbi decode,
# the empirical profile's fastest-growing residual is the soft-EM competing-paths
# E-step (C5c, alpha ~1.9 in isolation). The driver is the greedy re-route walk in
# `_consensus_alternative`: it walks highest-weight edges until it rejoins the
# read, bounded only by the read length `n`. In a large genome a k-mer de Bruijn
# graph has many SPURIOUS k-mer repeats (chance collisions whose incident edges
# are themselves well covered), so a greedy re-route increasingly wanders into
# graph bulk and runs the full `n` steps WITHOUT rejoining — the per-read cost
# climbs with graph density even though `n` is fixed, giving the super-linear
# term. These bounds cap the generation with the SAME discipline as the decode
# fix (top-B successors + a near-best band + a read-consistency exemption):
#
#   * `_soft_em_walk_band(k)` — a genome-INDEPENDENT cap on the re-route walk. A
#     read-consistent single-variant bubble at rung k diverges for ~k k-mers
#     before rejoining, so `_SOFT_EM_WALK_BAND_FACTOR * k` gives generous headroom
#     for multi-base variants while never growing with genome. The rejoin test is
#     applied BEFORE the band cutoff, so a real allele that returns to the read
#     within the band is ALWAYS captured; the band prunes only NON-rejoining
#     excursions, which yield no valid spliced alternative anyway. This is the
#     term that removes the super-linearity — the analog of the decode
#     `beam_score_margin`. A support-based test cannot substitute here: spurious-
#     repeat edges in a large genome are themselves high-coverage, so only the
#     rejoin-within-band criterion (read-consistency) separates a real competing
#     allele from graph-bulk wandering. This is the emission-exemption analog:
#     the bound removes WRONG (non-read-consistent) paths, never a real variant's
#     read-consistent competing path.
#
#   * `_SOFT_EM_ALT_SUCCESSOR_BOUND` — caps the incident branches examined per step
#     in EACH direction (out, in). A canonical DNA de Bruijn vertex has <= 4 out and
#     <= 4 in neighbors, so both the default `typemax` and the :scalable value 8 are
#     no-ops on DNA; it is a robustness guard bounding per-step work on pathological
#     high-branching (non-DNA / reduced-alphabet) inputs, matching the intent of the
#     decode `max_successors_per_state`.
#
# Both default to unbounded (`typemax(Int)`), reproducing the pre-td-e70t behavior
# byte-for-byte — the default for direct/unit-test callers and the exact decode
# tier. The `:scalable` corrector opts in ONLY where the width beam is already
# finite (`!beam_is_exact` — the dense/large reads that carry the residual), so
# exact-ML reads stay byte-identical, exactly like the #388 decode bounds.
const _SOFT_EM_WALK_BAND_FACTOR = 4
_soft_em_walk_band(k::Integer) = _SOFT_EM_WALK_BAND_FACTOR * Int(k)
const _SOFT_EM_ALT_SUCCESSOR_BOUND = 8

# Strictly-heavier sibling of the observed branch `a -> b` (excluding `b` and the
# incoming `avoid` vertex), or `nothing` when no strictly-better sibling exists (a
# linear region or a balanced variant — retained, not competed away). Iterates the
# out- then in-neighbors DIRECTLY (no intermediate label vector): finding the
# argmax is idempotent to a bidirectional neighbor appearing in both lists (same
# label, same `_edge_weight_between` value), so no dedup pass is needed — this
# removes the per-position allocation that dominated the C5c profile. `bound` caps
# the neighbors examined in each direction (a per-step work guard; on a canonical
# DNA vertex there are <= 4 each, so the default `typemax` and the :scalable 8 are
# both no-ops), matching the decode `max_successors_per_state`.
function _best_alternative_sibling(graph, a, b, avoid, w_ab::Float64, bound::Int)
    best_x, best_w = nothing, w_ab
    c = 0
    for x in MetaGraphsNext.outneighbor_labels(graph, a)
        c += 1
        c > bound && break
        (x == b || (avoid !== nothing && x == avoid)) && continue
        w = _edge_weight_between(graph, a, x)
        w === nothing && continue
        w > best_w && ((best_w, best_x) = (w, x))
    end
    c = 0
    for x in MetaGraphsNext.inneighbor_labels(graph, a)
        c += 1
        c > bound && break
        (x == b || (avoid !== nothing && x == avoid)) && continue
        w = _edge_weight_between(graph, a, x)
        w === nothing && continue
        w > best_w && ((best_w, best_x) = (w, x))
    end
    return best_x
end

# Highest-weight next vertex from `cur` (excluding the incoming `prev` and any
# vertex already in `seen`), or `nothing` when the walk dead-ends. Same direct
# out-then-in iteration (argmax idempotent to duplicates) and `bound` guard as
# `_best_alternative_sibling`.
function _best_next_step(graph, cur, prev, seen, bound::Int)
    nxt, nxt_w = nothing, 0.0
    c = 0
    for y in MetaGraphsNext.outneighbor_labels(graph, cur)
        c += 1
        c > bound && break
        (y == prev || y in seen) && continue
        w = _edge_weight_between(graph, cur, y)
        w === nothing && continue
        w > nxt_w && ((nxt_w, nxt) = (w, y))
    end
    c = 0
    for y in MetaGraphsNext.inneighbor_labels(graph, cur)
        c += 1
        c > bound && break
        (y == prev || y in seen) && continue
        w = _edge_weight_between(graph, cur, y)
        w === nothing && continue
        w > nxt_w && ((nxt_w, nxt) = (w, y))
    end
    return nxt
end

# Generate ONE competing alternative to `observed` by re-routing its first
# clearly-weak branch through the highest-weight sibling and walking greedily
# until rejoining `observed`. Returns the spliced label path, or `nothing` when
# no weak branch exists (a linear region, or a balanced variant whose sibling is
# not strictly better — which is retained, not competed away).
#
# `walk_band` caps a NON-rejoining re-route to a genome-independent number of steps
# (see `_soft_em_walk_band`); `successor_bound` caps the incident branches examined
# per step (see `_SOFT_EM_ALT_SUCCESSOR_BOUND`). Both default to `typemax(Int)`
# (unbounded — the pre-td-e70t behavior). Always bounded by the observed length so
# it terminates regardless.
function _consensus_alternative(observed, graph;
        walk_band::Int = typemax(Int),
        successor_bound::Int = typemax(Int))
    n = length(observed)
    n < 3 && return nothing
    firstpos = Dict{Any, Int}()
    for i in 1:n
        haskey(firstpos, observed[i]) || (firstpos[observed[i]] = i)
    end
    for i in 1:(n - 1)
        a, b = observed[i], observed[i + 1]
        w_ab = _edge_weight_between(graph, a, b)
        w_ab === nothing && continue
        best_x = _best_alternative_sibling(
            graph, a, b, i > 1 ? observed[i - 1] : nothing, w_ab, successor_bound)
        best_x === nothing && continue   # no strictly-better sibling ⇒ not weak
        # Greedy highest-weight walk from the sibling until we rejoin `observed`
        # at or after position i+1. `seen` is a Set so the membership test is O(1)
        # rather than O(length(mid)) (removes a hidden per-step super-linear term).
        mid = Any[best_x]
        seen = Set{Any}(mid)
        prev, cur = a, best_x
        rejoin = get(firstpos, best_x, 0) >= i + 1 ? firstpos[best_x] : 0
        steps = 0
        while rejoin == 0 && steps < n
            steps += 1
            nxt = _best_next_step(graph, cur, prev, seen, successor_bound)
            nxt === nothing && break
            # Rejoin is tested BEFORE the band cutoff: a read-consistent competing
            # allele that returns to the observed path (at/after i+1) is always
            # captured, so the band never drops a real variant's competing path.
            if get(firstpos, nxt, 0) >= i + 1
                rejoin = firstpos[nxt]
                break
            end
            # Band cutoff (genome-independent): beyond `walk_band` non-rejoining
            # steps the walk is threading spurious k-mer repeats / graph bulk, not a
            # competing allele — abandon it (yields no valid rejoined splice anyway).
            steps >= walk_band && break
            push!(mid, nxt)
            push!(seen, nxt)
            prev, cur = cur, nxt
        end
        if rejoin >= i + 1
            return vcat(observed[1:i], mid, observed[rejoin:end])
        end
    end
    return nothing
end

"""
    accumulate_competing_paths!(accumulator, read, graph, k; graph_mode,
        walk_band, successor_bound) -> accumulator

Soft-EM v2 competing-paths E-step for ONE read (td-e70t). Builds the read's
COMPETING candidate paths — the observed read path and (when the observed path
takes a weak branch) a consensus alternative re-routed through the best-supported
sibling (`_consensus_alternative`) — scores each by its normalized-transition
log-probability on the graph (`_path_transition_logscore`, which reads the
soft-weight-aware `compute_edge_weight`), converts the scores to responsibilities
with a stable softmax (`Mycelia.Rhizomorph.path_responsibility`), and accumulates
each candidate's graph-oriented edges weighted by its responsibility.

The responsibility split alone would decay a real but skewed minority allele
toward zero (`W_min' = N*W/(W_maj+W_min)`); variation is preserved by the SUPPORT
FLOOR applied when the accumulator is registered in the M-step
(`register_soft_edge_weights!`), which clamps every `>= SOFT_EM_MIN_SUPPORT` edge
to at least its raw coverage. So this E-step is free to split responsibility while
supported edges are held at raw and only unsupported (error) edges decay. A
balanced variant produces no strictly-better sibling, so no alternative is
generated and the observed path keeps responsibility 1.0. Never throws (guarded),
so the decode is always safe.

`walk_band` / `successor_bound` bound the alternative GENERATION (td-e70t speed
residual C5c): `walk_band` caps a non-rejoining re-route to a genome-independent
number of steps and `successor_bound` caps the incident branches examined per step
(see `_soft_em_walk_band` / `_SOFT_EM_ALT_SUCCESSOR_BOUND`). Both default to
`typemax(Int)` (unbounded — byte-identical to the pre-bound behavior); the
`:scalable` corrector passes finite bounds only where the width beam is already
finite, so exact-ML reads are unchanged and no real variant's read-consistent
competing path is dropped (the rejoin test precedes the band cutoff).
"""
function accumulate_competing_paths!(
        accumulator::Mycelia.Rhizomorph.SoftEdgeWeightAccumulator,
        read::FASTX.FASTQ.Record,
        graph,
        k::Int;
        graph_mode::Symbol = :canonical,
        walk_band::Int = typemax(Int),
        successor_bound::Int = typemax(Int)
)
    observed = _read_resolved_labels(read, graph, k; graph_mode = graph_mode)
    isempty(observed) && return accumulator

    alternative = try
        _consensus_alternative(
            observed, graph; walk_band = walk_band, successor_bound = successor_bound)
    catch
        nothing
    end

    if alternative === nothing || alternative == observed
        # No genuine competition: the observed path owns all responsibility.
        Mycelia.Rhizomorph.accumulate_path_probability!(
            accumulator, _graph_oriented_edges(observed, graph), 1.0)
        return accumulator
    end

    candidate_paths = (observed, alternative)
    scores = Float64[_path_transition_logscore(labels, graph) for labels in candidate_paths]
    finite = Float64[s for s in scores if isfinite(s)]
    if isempty(finite)
        # Fall back to the observed path if scoring degenerated.
        Mycelia.Rhizomorph.accumulate_path_probability!(
            accumulator, _graph_oriented_edges(observed, graph), 1.0)
        return accumulator
    end
    for (idx, labels) in enumerate(candidate_paths)
        isfinite(scores[idx]) || continue
        responsibility = Mycelia.Rhizomorph.path_responsibility(scores[idx], finite)
        Mycelia.Rhizomorph.accumulate_path_probability!(
            accumulator, _graph_oriented_edges(labels, graph), responsibility)
    end
    return accumulator
end

"""
Generate alternative qualmer paths through the graph using quality-aware probabilistic sampling.
"""
function generate_alternative_qualmer_paths(
        record::FASTX.FASTQ.Record, graph, k::Int; num_samples::Int = 5,
        graph_mode::Symbol = :canonical)::Vector{Vector{<:Qualmer}}
    paths = Vector{Vector{Qualmer}}()

    try
        # Generate qualmers using existing qualmer-analysis functions
        original_qualmers = collect(qualmers_unambiguous(record, k))
        if isempty(original_qualmers)
            return Vector{Vector{Qualmer}}()
        end

        original_path = Qualmer[]
        for (qmer, pos) in original_qualmers
            resolved = _resolve_qualmer_for_graph(graph, qmer; graph_mode = graph_mode)
            if resolved !== nothing
                push!(original_path, resolved)
            end
        end
        if isempty(original_path)
            return Vector{Vector{Qualmer}}()
        end

        # For each sample, create alternative path using quality-aware sampling
        for sample in 1:num_samples
            alternative_path = copy(original_path)

            # Select positions to modify based on quality confidence
            low_confidence_positions = Int[]
            for (i, qmer) in enumerate(original_path)
                current_kmer = qmer.kmer
                if haskey(graph, current_kmer)
                    vertex_data = graph[current_kmer]
                    # Prefer modifying low-confidence positions
                    if _vertex_joint_probability(vertex_data) < 0.8 ||
                       _vertex_mean_quality(vertex_data) < 30.0
                        push!(low_confidence_positions, i)
                    end
                end
            end

            # If no low confidence positions, randomly select some
            modification_positions = if !isempty(low_confidence_positions)
                Random.shuffle(low_confidence_positions)[1:min(
                    2, length(low_confidence_positions))]
            else
                Random.shuffle(1:length(original_path))[1:min(
                    2, length(original_path))]
            end

            for pos in modification_positions
                if pos <= length(alternative_path)
                    current_qualmer = alternative_path[pos]
                    current_kmer = current_qualmer.kmer

                    if haskey(graph, current_kmer)
                        # Get high-quality neighbor k-mers for potential transitions
                        vertex_data = graph[current_kmer]
                        vertex_coverage = _vertex_coverage(vertex_data)
                        all_vertices = collect(MetaGraphsNext.labels(graph))

                        # Filter neighbors that could be valid transitions
                        valid_neighbors = []
                        neighbor_weights = Float64[]

                        for neighbor_kmer in all_vertices
                            neighbor_data = graph[neighbor_kmer]

                            # Quality-aware weighting combining multiple factors
                            weight = _vertex_joint_probability(neighbor_data) *
                                     (_vertex_mean_quality(neighbor_data) / 60.0) *  # Normalize quality (max ~60)
                                     (_vertex_coverage(neighbor_data) /
                                      max(1, vertex_coverage))  # Relative coverage

                            if weight > 0.01  # Minimum threshold for consideration
                                push!(valid_neighbors, neighbor_kmer)
                                push!(neighbor_weights, weight)
                            end
                        end

                        if !isempty(valid_neighbors)
                            # Sample neighbor based on quality-aware weights
                            neighbor_weights ./= sum(neighbor_weights)  # Normalize
                            selected_idx = StatsBase.sample(
                                1:length(valid_neighbors), StatsBase.Weights(neighbor_weights))

                            # Create new qualmer using the selected k-mer and quality from original
                            selected_kmer = valid_neighbors[selected_idx]
                            original_qualities = current_qualmer.qualities

                            new_qualmer = Qualmer(selected_kmer, original_qualities)
                            alternative_path[pos] = new_qualmer
                        end
                    end
                end
            end

            push!(paths, alternative_path)
        end

    catch e
        # If path generation fails, return empty paths
        @debug "Alternative path generation failed: $e"
        return Vector{Vector{Qualmer}}()
    end

    return paths
end

"""
Convert sequence with quality to qualmer path representation using graph vertices.
"""
function sequence_to_qualmer_path(sequence::String, quality::String, graph, k::Int;
        graph_mode::Symbol = :canonical)::Vector{<:Qualmer}
    path = Qualmer[]

    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels) || length(sequence) < k
        return path
    end

    kmer_type = typeof(first(labels))

    for i in 1:(length(sequence) - k + 1)
        kmer_str = sequence[i:(i + k - 1)]
        qual_segment = quality[i:(i + k - 1)]
        quality_scores = UInt8.(codeunits(qual_segment)) .- UInt8(33)

        qmer = Qualmer(kmer_type(kmer_str), quality_scores)
        resolved = _resolve_qualmer_for_graph(graph, qmer; graph_mode = graph_mode)
        if resolved !== nothing
            push!(path, resolved)
        end
    end

    return path
end

"""
Convert sequence to k-mer path representation (simplified version for backward compatibility).
"""
function sequence_to_kmer_path(sequence::String, k::Int)::Vector{String}
    path = String[]

    for i in 1:(length(sequence) - k + 1)
        push!(path, sequence[i:(i + k - 1)])
    end

    return path
end

function _kmer_to_string(kmer)::String
    return kmer isa AbstractString ? String(kmer) : string(kmer)
end

"""
Convert k-mer path back to sequence string.
"""
function kmer_path_to_sequence(path::AbstractVector, k::Int)::String
    if isempty(path)
        return ""
    end

    first_kmer = _kmer_to_string(path[1])
    if isempty(first_kmer)
        return ""
    end

    buffer = IOBuffer()
    write(buffer, first_kmer)

    for i in 2:length(path)
        kmer_str = _kmer_to_string(path[i])
        if !isempty(kmer_str)
            write(buffer, kmer_str[end])
        end
    end

    return String(take!(buffer))
end

"""
Convert a generic path back to sequence string (alias for k-mer paths).
"""
function path_to_sequence(path::AbstractVector, k::Int)::String
    return kmer_path_to_sequence(path, k)
end

"""
Convert qualmer path back to BioSequence and quality vector.
Returns a named tuple with sequence and quality_scores.
"""
function qualmer_path_to_biosequence(path::Vector{<:Qualmer})
    if isempty(path)
        return (sequence = BioSequences.LongDNA{4}(""), quality_scores = Int8[])
    end

    # Determine sequence type from first qualmer
    first_kmer = path[1].kmer
    sequence_type = if first_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif first_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else  # AAKmer
        BioSequences.LongAA
    end

    # Start with first k-mer
    sequence_symbols = collect(path[1].kmer)
    quality_scores = collect(Int8.(path[1].qualities))

    # Add last symbol of each subsequent k-mer
    for i in 2:length(path)
        push!(sequence_symbols, path[i].kmer[end])
        push!(quality_scores, Int8(path[i].qualities[end]))
    end

    # Construct the BioSequence
    final_sequence = sequence_type(sequence_symbols)

    return (sequence = final_sequence, quality_scores = quality_scores)
end

"""
Convert k-mer path back to BioSequence.
"""
function kmer_path_to_biosequence(path::Vector{<:Kmers.Kmer})
    if isempty(path)
        return BioSequences.LongDNA{4}("")
    end

    # Determine sequence type from first k-mer
    first_kmer = path[1]
    sequence_type = if first_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif first_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else  # AAKmer
        BioSequences.LongAA
    end

    # Start with first k-mer
    sequence_symbols = collect(path[1])

    # Add last symbol of each subsequent k-mer
    for i in 2:length(path)
        push!(sequence_symbols, path[i][end])
    end

    # Construct the BioSequence
    return sequence_type(sequence_symbols)
end
