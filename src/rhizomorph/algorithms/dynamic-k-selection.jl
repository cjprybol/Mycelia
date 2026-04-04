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
    return FASTX.FASTA.sequence(record)
end

function _dynamic_k_sequence_string(record::FASTX.FASTQ.Record)
    return FASTX.FASTQ.sequence(record)
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

function _calculate_dynamic_k_sparsity(sequences::Vector{String}, k::Int)::Float64
    observed_kmers = Set{String}()

    for sequence in sequences
        if length(sequence) >= k
            for start_index in 1:(length(sequence) - k + 1)
                push!(observed_kmers, sequence[start_index:(start_index + k - 1)])
            end
        end
    end

    if isempty(observed_kmers)
        return 1.0
    end

    alphabet_size = _dynamic_k_alphabet_size(sequences)
    theoretical_space = float(alphabet_size)^k
    sparsity = 1.0 - (length(observed_kmers) / theoretical_space)
    return clamp(sparsity, 0.0, 1.0)
end

function _dynamic_k_errors_are_singletons(
        sequences::Vector{String},
        k::Int;
        singleton_threshold::Int = 2
)::Bool
    kmer_counts = Dict{String, Int}()

    for sequence in sequences
        if length(sequence) >= k
            for start_index in 1:(length(sequence) - k + 1)
                kmer = sequence[start_index:(start_index + k - 1)]
                kmer_counts[kmer] = get(kmer_counts, kmer, 0) + 1
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

    sequence_lengths = length.(sequences)
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
        sparsity = _calculate_dynamic_k_sparsity(sequences, k)
        separated_singletons = _dynamic_k_errors_are_singletons(
            sequences,
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
