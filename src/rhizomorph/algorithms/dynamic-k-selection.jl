"""
    DynamicKmerSelectionPlan

Scaffold output for dynamic k-mer selection in Rhizomorph.

This plan captures the initial `k` recommendation, the progressive prime
sequence to try next, and the heuristic measurements that informed the choice.
It is intentionally lightweight so future multi-k assembly orchestration can
reuse it without committing to a specific optimization backend.
"""
struct DynamicKmerSelectionPlan
    initial_k::Int
    candidate_ks::Vector{Int}
    search_space::Vector{Int}
    sequence_count::Int
    min_sequence_length::Int
    median_sequence_length::Float64
    max_candidate_k::Int
    sparsity_by_k::Dict{Int, Float64}
    singleton_separation_by_k::Dict{Int, Bool}
end

function _dynamic_k_sequence_string(sequence::BioSequences.BioSequence)
    return String(sequence)
end

function _dynamic_k_sequence_string(record::FASTX.FASTA.Record)
    return String(FASTX.FASTA.sequence(record))
end

function _dynamic_k_sequence_string(record::FASTX.FASTQ.Record)
    return String(FASTX.FASTQ.sequence(record))
end

function _dynamic_k_sequence_string(sequence::AbstractString)
    return String(sequence)
end

function _dynamic_k_sequence_string(observation)
    throw(ArgumentError("Unsupported observation type for dynamic k-mer selection: $(typeof(observation))"))
end

function _collect_dynamic_k_sequences(observations)
    sequences = String[]
    for observation in observations
        sequence = _dynamic_k_sequence_string(observation)
        if !isempty(sequence)
            push!(sequences, sequence)
        end
    end
    return sequences
end

function _dynamic_k_character_sequences(sequences::Vector{String})::Vector{Vector{Char}}
    return [collect(sequence) for sequence in sequences]
end

function _dynamic_k_search_space(min_k::Int, max_k::Int)::Vector{Int}
    if max_k < 1
        return Int[]
    end

    candidate_min_k = max(2, min_k)
    prime_ks = if candidate_min_k <= max_k
        Mycelia.Primes.primes(candidate_min_k, max_k)
    else
        Int[]
    end

    if !isempty(prime_ks)
        return prime_ks
    end

    fallback_primes = Mycelia.Primes.primes(2, max_k)
    if !isempty(fallback_primes)
        return [last(fallback_primes)]
    end

    return [1]
end

function _dynamic_k_alphabet_size(sequences::Vector{String})::Int
    alphabet = Set{Char}()
    for sequence in sequences
        for character in sequence
            push!(alphabet, uppercase(character))
        end
    end
    return max(1, length(alphabet))
end

function _dynamic_k_window_hash(characters::AbstractVector{Char}, start_index::Int, k::Int)::UInt64
    return UInt64(hash(view(characters, start_index:(start_index + k - 1))))
end

function _calculate_dynamic_k_sparsity(
        character_sequences::Vector{Vector{Char}},
        alphabet_size::Int,
        k::Int
)::Float64
    observed_kmers = Set{UInt64}()

    for characters in character_sequences
        if length(characters) >= k
            for start_index in 1:(length(characters) - k + 1)
                push!(observed_kmers, _dynamic_k_window_hash(characters, start_index, k))
            end
        end
    end

    if isempty(observed_kmers)
        return 1.0
    end

    theoretical_space = float(alphabet_size)^k
    sparsity = 1.0 - (length(observed_kmers) / theoretical_space)
    return clamp(sparsity, 0.0, 1.0)
end

function _dynamic_k_errors_are_singletons(
        character_sequences::Vector{Vector{Char}},
        k::Int;
        singleton_threshold::Int = 2
)::Bool
    kmer_counts = Dict{UInt64, Int}()

    for characters in character_sequences
        if length(characters) >= k
            for start_index in 1:(length(characters) - k + 1)
                kmer_hash = _dynamic_k_window_hash(characters, start_index, k)
                kmer_counts[kmer_hash] = get(kmer_counts, kmer_hash, 0) + 1
            end
        end
    end

    if isempty(kmer_counts)
        return false
    end

    coverage_values = collect(values(kmer_counts))
    singleton_count = count(value -> value <= singleton_threshold, coverage_values)
    total_unique = length(coverage_values)

    if singleton_count == 0 || total_unique == singleton_count
        return false
    end

    singleton_fraction = singleton_count / total_unique
    non_singleton_min = minimum(filter(value -> value > singleton_threshold, coverage_values))
    return singleton_fraction > 0.1 && non_singleton_min > singleton_threshold * 2
end

"""
    dynamic_k_prime_pattern(start_k::Int = 11; max_k::Int = 101, initial_step::Int = 2)

Generate a progressively spaced sequence of candidate k-mer sizes.

The scaffold starts at a prime `k`, then increases the step size by two after
each successful prime jump. This mirrors the Phase 2 roadmap's prime
progression idea while keeping the implementation small and deterministic.
"""
function dynamic_k_prime_pattern(start_k::Int = 11; max_k::Int = 101, initial_step::Int = 2)::Vector{Int}
    if initial_step < 1
        throw(ArgumentError("initial_step must be positive, got $initial_step"))
    end

    if start_k < 2
        return start_k > max_k ? Int[] : [start_k]
    end

    if !Mycelia.Primes.isprime(start_k)
        start_k = Mycelia.Primes.nextprime(start_k)
    end

    if start_k > max_k
        return Int[]
    end

    candidate_ks = Int[start_k]
    current_k = start_k
    step = initial_step

    while true
        next_k = current_k + step
        if next_k > max_k
            break
        end
        if !Mycelia.Primes.isprime(next_k)
            break
        end
        push!(candidate_ks, next_k)
        current_k = next_k
        step += 2
    end

    return candidate_ks
end

"""
    select_dynamic_kmer_plan(
        observations;
        min_k::Int = 3,
        max_k::Union{Int, Nothing} = nothing,
        max_search_k::Int = 101,
        initial_step::Int = 2,
        sparsity_threshold::Float64 = 0.5,
        singleton_threshold::Int = 2
    ) -> DynamicKmerSelectionPlan

Create a lightweight dynamic k-mer selection plan from FASTA, FASTQ, string, or
BioSequence observations.

The scaffold searches feasible prime `k` values bounded by the shortest
available observation, chooses the first `k` whose sparsity and singleton
separation exceed configurable thresholds, and then emits a progressive prime
sequence for future multi-k assembly work.
"""
function select_dynamic_kmer_plan(
        observations;
        min_k::Int = 3,
        max_k::Union{Int, Nothing} = nothing,
        max_search_k::Int = 101,
        initial_step::Int = 2,
        sparsity_threshold::Float64 = 0.5,
        singleton_threshold::Int = 2
)::DynamicKmerSelectionPlan
    if min_k < 1
        throw(ArgumentError("min_k must be positive, got $min_k"))
    end
    if max_search_k < 1
        throw(ArgumentError("max_search_k must be positive, got $max_search_k"))
    end
    if max_k !== nothing && max_k < 1
        throw(ArgumentError("max_k must be positive when provided, got $max_k"))
    end

    sequences = _collect_dynamic_k_sequences(observations)
    if isempty(sequences)
        throw(ArgumentError("Dynamic k-mer selection requires at least one non-empty observation"))
    end

    character_sequences = _dynamic_k_character_sequences(sequences)
    sequence_lengths = length.(sequences)
    alphabet_size = _dynamic_k_alphabet_size(sequences)
    bounded_max_k = minimum(sequence_lengths)
    bounded_max_k = min(bounded_max_k, max_search_k)
    if max_k !== nothing
        bounded_max_k = min(bounded_max_k, max_k)
    end

    search_space = _dynamic_k_search_space(min_k, bounded_max_k)
    selected_k = nothing
    sparsity_by_k = Dict{Int, Float64}()
    singleton_separation_by_k = Dict{Int, Bool}()

    for k in search_space
        sparsity = _calculate_dynamic_k_sparsity(character_sequences, alphabet_size, k)
        separated_singletons = _dynamic_k_errors_are_singletons(
            character_sequences,
            k;
            singleton_threshold = singleton_threshold
        )
        sparsity_by_k[k] = sparsity
        singleton_separation_by_k[k] = separated_singletons

        if selected_k === nothing && sparsity > sparsity_threshold && separated_singletons
            selected_k = k
        end
    end

    if selected_k === nothing
        selected_k = first(search_space)
    end

    candidate_ks = dynamic_k_prime_pattern(
        selected_k;
        max_k = bounded_max_k,
        initial_step = initial_step
    )
    if isempty(candidate_ks)
        candidate_ks = [selected_k]
    end

    return DynamicKmerSelectionPlan(
        selected_k,
        candidate_ks,
        search_space,
        length(sequences),
        minimum(sequence_lengths),
        Float64(Statistics.median(sequence_lengths)),
        bounded_max_k,
        sparsity_by_k,
        singleton_separation_by_k
    )
end

# --- Residual-error estimation + survival-based re-assembly k selection ---------
#
# `select_dynamic_kmer_plan` (above) picks a START k for the corrector's k-ladder
# from raw observations; it is NOT the right tool for the *re-assembly* of already
# CORRECTED reads (empirically it returns a flat floor k regardless of error,
# because its singleton-separation heuristic rarely fires on real read sets). The
# re-assembly of corrected reads needs a k that keeps the contig graph connected:
# high k for clean (Illumina) corrected reads, but a LOWER k for high-error long
# reads (nanopore) whose residual errors — substitutions and, dominantly, indels —
# shatter a high-k de Bruijn graph. The functions below estimate the residual error
# reference-free and map it to a prime k via the k-mer survival model.

_read_quality_scores(::Any) = nothing
function _read_quality_scores(record::FASTX.FASTQ.Record)
    quality = FASTX.FASTQ.quality(record)
    return Int[Int(character) - 33 for character in quality]
end

# Mean per-base error probability from Phred quality (e = 10^(-Q/10)), or `nothing`
# when no read carries usable quality (FASTA / string / BioSequence input). A
# placeholder constant-Q40 does no harm: it yields e ≈ 1e-4, which the caller's
# `max` with the k-mer estimate ignores.
function _quality_residual_error(reads)::Union{Float64, Nothing}
    total_error = 0.0
    base_count = 0
    for record in reads
        scores = _read_quality_scores(record)
        scores === nothing && continue
        for score in scores
            total_error += 10.0^(-score / 10.0)
            base_count += 1
        end
    end
    return base_count == 0 ? nothing : clamp(total_error / base_count, 0.0, 0.499)
end

# k-mer spectrum estimate: genomic k-mer positions (error-free) recur at ~coverage
# and are "solid" (count >= solid_min); erroneous positions are singletons. The
# solid FRACTION of k-mer occurrences approximates (1-e)^k_ref, so
# e ≈ 1 - solid_fraction^(1/k_ref). Because an indel disrupts many downstream
# k-mers, this (correctly) reports a higher effective error for indel-heavy reads.
function _kmer_spectrum_residual_error(reads, k_ref::Int, solid_min::Int)::Float64
    sequences = _collect_dynamic_k_sequences(reads)
    isempty(sequences) && return 0.0
    character_sequences = _dynamic_k_character_sequences(sequences)
    counts = Dict{UInt64, Int}()
    total_occurrences = 0
    for characters in character_sequences
        if length(characters) >= k_ref
            for start_index in 1:(length(characters) - k_ref + 1)
                kmer_hash = _dynamic_k_window_hash(characters, start_index, k_ref)
                counts[kmer_hash] = get(counts, kmer_hash, 0) + 1
                total_occurrences += 1
            end
        end
    end
    total_occurrences == 0 && return 0.0
    solid_occurrences = 0
    for occurrence_count in values(counts)
        if occurrence_count >= solid_min
            solid_occurrences += occurrence_count
        end
    end
    solid_fraction = solid_occurrences / total_occurrences
    solid_fraction <= 0.0 && return 0.499
    return clamp(1.0 - solid_fraction^(1.0 / k_ref), 0.0, 0.499)
end

"""
    estimate_residual_error(reads; k_ref = 13, solid_min = 2) -> Float64

Reference-free estimate of the per-base residual error rate of `reads`, used to
choose a re-assembly k that keeps the corrected-read graph connected. Combines two
signals and returns the more conservative (higher-error) one, clamped to `[0, 0.5)`:

- **k-mer spectrum** (always available): solid-fraction of k-mer occurrences at
  `k_ref` inverted through the survival model `(1-e)^k_ref`.
- **per-base Q-values** (FASTQ input only): `mean(10^(-Q/10))`.

`reads` may be FASTQ/FASTA records, strings, or `BioSequence`s.
"""
function estimate_residual_error(reads; k_ref::Int = 13, solid_min::Int = 2)::Float64
    kmer_error = _kmer_spectrum_residual_error(reads, k_ref, solid_min)
    quality_error = _quality_residual_error(reads)
    estimate = quality_error === nothing ? kmer_error : max(kmer_error, quality_error)
    return clamp(estimate, 0.0, 0.499)
end

function _largest_prime_at_most(n::Int)::Int
    candidate = max(n, 2)
    while candidate >= 2 && !Mycelia.Primes.isprime(candidate)
        candidate -= 1
    end
    return max(candidate, 2)
end

"""
    select_reassembly_k(reads, ceiling_k; floor_k = 7) -> Int

Choose a re-assembly k for corrected `reads`, bounded by `[floor_k, ceiling_k]`.
The residual error `e = estimate_residual_error(reads)` is mapped to the largest k
whose expected k-mer survival stays above 0.5 via `k <= log(0.5)/log(1-e)`.

When the requested `ceiling_k` already survives (clean / low-error reads), it is
returned UNCHANGED — this keeps clean/Illumina behavior byte-identical and
preserves the corrector's final-pass graph-reuse eligibility (which requires
`reassembly_k == final_graph_k`). Only when residual error forces a lower k is the
adapted value snapped down to the nearest prime (the `:scalable` k-ladder is prime,
so an adaptively-chosen re-assembly k stays on primes too). High-error long reads
drop toward `floor_k`. Never exceeds `ceiling_k`.
"""
function select_reassembly_k(reads, ceiling_k::Int; floor_k::Int = 7)::Int
    error_rate = estimate_residual_error(reads)
    # Low residual error: the requested ceiling survives — honor it unchanged
    # (identical clean/Illumina behavior; preserves graph-reuse eligibility).
    error_rate <= 1.0e-6 && return ceiling_k
    survival_k = Int(floor(log(0.5) / log(1.0 - error_rate)))
    survival_k >= ceiling_k && return ceiling_k
    # Residual error forces a lower k: drop to the largest prime that keeps the
    # graph connected, floored (but never above the ceiling).
    effective_floor = min(floor_k, ceiling_k)
    target_k = clamp(survival_k, effective_floor, ceiling_k)
    return clamp(_largest_prime_at_most(target_k), min(effective_floor, target_k), ceiling_k)
end
