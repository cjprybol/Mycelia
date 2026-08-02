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
# `select_dynamic_kmer_plan` (above) is a SCAFFOLD with NO production call site, as
# its own docstrings at the top of this file already say ("Scaffold output for
# dynamic k-mer selection", "The scaffold starts at a prime `k`", "The scaffold
# searches feasible prime `k` values"). `grep -rn select_dynamic_kmer_plan src/`
# returns only its docstring, its definition, and this comment — every caller is a
# test. It was DESIGNED to pick a START k for a corrector k-ladder from raw
# observations, but nothing in production consumes that: the iterative corrector's
# ladder comes from `mycelia_iterative_assemble`'s `k_ladder` / `n_k_rungs` kwargs.
# (Earlier wording here asserted the wiring as fact, which made a scaffold read as
# load-bearing.) It is in any case NOT the right tool for the *re-assembly* of
# already CORRECTED reads (empirically it returns a flat floor k regardless of error,
# because its singleton-separation heuristic rarely fires on real read sets). The
# re-assembly of corrected reads needs a k that keeps the contig graph connected:
# high k for clean (Illumina) corrected reads, but a LOWER k for high-error long
# reads (nanopore) whose residual errors — substitutions and, dominantly, indels —
# shatter a high-k de Bruijn graph. `estimate_residual_error` (below) estimates that
# residual error reference-free via the k-mer survival model, but it is a DIAGNOSTIC
# / secondary signal only: the production chooser `select_reassembly_k` keys off the
# COVERAGE-AWARE connectivity criterion (`median_solid_kmer_multiplicity`) instead,
# which measures backbone coverage directly (see its docstring for why).

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

Reference-free estimate of the per-base residual error rate of `reads`. A DIAGNOSTIC
/ secondary signal — the production re-assembly-k chooser `select_reassembly_k` keys
off the coverage-aware connectivity criterion (`median_solid_kmer_multiplicity`), not
this estimate; this remains as an interpretable, reference-free error readout (and for
direct callers/tests). Combines two signals and returns the more conservative
(higher-error) one, clamped to `[0, 0.5)`:

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

# --- k-mer count spectrum over reads (shared by both connectivity measures) -----
_dna_complement_char(c::Char)::Char = c == 'A' ? 'T' :
                                      c == 'T' ? 'A' :
                                      c == 'C' ? 'G' :
                                      c == 'G' ? 'C' :
                                      c == 'a' ? 't' :
                                      c == 't' ? 'a' : c == 'c' ? 'g' : c == 'g' ? 'c' : c

# Canonical hash of a k-length window: min(hash(window), hash(revcomp(window))). The
# corrected-read RE-ASSEMBLY runs on a DOUBLESTRAND (canonical) de Bruijn graph, so
# the connectivity measure must ALSO be canonical — otherwise a forward read and the
# reverse-complement read of the SAME locus are counted as two different k-mers,
# HALVING the measured coverage on both-strand read sets and spuriously triggering
# k-adaptation on well-covered reads (the reuse-oracle regime).
function _dynamic_k_canonical_window_hash(
        characters::AbstractVector{Char}, start_index::Int, k::Int)::UInt64
    forward_hash = _dynamic_k_window_hash(characters, start_index, k)
    revcomp = Vector{Char}(undef, k)
    for offset in 0:(k - 1)
        revcomp[k - offset] = _dna_complement_char(uppercase(characters[start_index + offset]))
    end
    reverse_hash = UInt64(hash(revcomp))
    return min(forward_hash, reverse_hash)
end

function _kmer_count_spectrum(reads, k::Int; canonical::Bool = true)::Dict{UInt64, Int}
    counts = Dict{UInt64, Int}()
    k < 1 && return counts
    sequences = _collect_dynamic_k_sequences(reads)
    isempty(sequences) && return counts
    character_sequences = _dynamic_k_character_sequences(sequences)
    for characters in character_sequences
        if length(characters) >= k
            for start_index in 1:(length(characters) - k + 1)
                kmer_hash = canonical ?
                            _dynamic_k_canonical_window_hash(characters, start_index, k) :
                            _dynamic_k_window_hash(characters, start_index, k)
                counts[kmer_hash] = get(counts, kmer_hash, 0) + 1
            end
        end
    end
    return counts
end

# CAVEAT (td-jt7r, PR #444) — the docstring below contains a REFUTED empirical
# claim. Attributed per claim, because the paragraph bundles several:
#
#   * "solid k-mers ... i.e. the genomic backbone; singleton error k-mers are
#     dropped" and "a robust estimate of the backbone coverage `C·(1-e)^k`":
#     REFUTED. At non-trivial coverage an ERROR k-mer is also observed
#     `>= solid_min` times, so error k-mers ENTER the solid set; being far more
#     numerous than backbone k-mers, they DOMINATE the median. The returned value
#     is then an error-population statistic, not a backbone-coverage estimate.
#     Measured by `benchmarking/indel_toy_overcorrection_diagnostic.jl`'s
#     `solid_composition` audit (`solid_composition.csv`, columns `median_genomic`
#     vs `median_error` vs `median_all`). At 30x the genomic-solid median is 27
#     while this function returns 4.0.
#   * "its absolute scale tracks read coverage": REFUTED as stated — the statistic
#     is NON-MONOTONE in coverage. Raising coverage recruits error k-mers into the
#     solid set faster than it lifts the backbone median, so the returned value can
#     FALL as coverage RISES.
#   * "Canonical (both-strand) counting is REQUIRED ...": NOT affected. That is a
#     separate, still-valid claim about the equivalence class, unrelated to the
#     solid-set composition above.
#   * The two "empirical, 30x corrected reads" bullets (~12–15 clean vs ~2–4
#     shattered): these are the numbers `connectivity_floor = 6.0` was fitted to,
#     and they are two SEPARATE measurements on two different read sets, not one
#     coverage-controlled sweep. They stand as the observations they were, but they
#     do NOT establish that the statistic tracks backbone coverage in general.
#
# The code fix is deferred to its own bead; the docstring is marked here rather than
# left silently asserting a property this PR's own evidence refutes. Read the return
# value as "median of the solid set", nothing stronger.
"""
    median_solid_kmer_multiplicity(reads, k; solid_min = 2) -> Float64

Coverage-aware connectivity measure at size `k`, and the signal `select_reassembly_k`
uses. Counts every length-`k` CANONICAL window over all `reads`, keeps the "solid"
k-mers (occurrence count `>= solid_min`, i.e. the genomic backbone; singleton error
k-mers are dropped), and returns the MEDIAN count over that solid set — a robust
estimate of the backbone coverage `C·(1-e)^k`.

Canonical (both-strand) counting is REQUIRED: the corrected-read re-assembly runs on
a DOUBLESTRAND de Bruijn graph, so a forward read and the reverse-complement read of
the same locus must count as ONE k-mer; a raw-string count would halve the measured
coverage on both-strand read sets and spuriously trigger adaptation on well-covered
reads (the reuse-oracle regime).

Discriminates well-covered from shattered backbones (empirical, 30x corrected reads):

- clean / well-corrected backbone (reuse-oracle, 5% Illumina, k 11–21): ~12–15 —
  the backbone recurs ~coverage, comfortably above the connectivity floor, so the
  ceiling is honored.
- shattered backbone (raw / poorly-corrected high-error long reads): ~2–4 at every
  k — the backbone barely survives, so median-solid sits below the floor and `k`
  drops.

Returns `0.0` when no k-mer is solid. NOTE: this is a per-k backbone-coverage floor;
its absolute scale tracks read coverage (a ~6 floor assumes non-trivial, ~≥10x
coverage — the regime the corrector targets).
"""
function median_solid_kmer_multiplicity(reads, k::Int; solid_min::Int = 2)::Float64
    counts = _kmer_count_spectrum(reads, k; canonical = true)
    solid_counts = Int[c for c in values(counts) if c >= solid_min]
    isempty(solid_counts) && return 0.0
    return Float64(Statistics.median(solid_counts))
end

"""
    effective_kmer_coverage(reads, k) -> Float64

DIAGNOSTIC / SECONDARY signal: the MEAN occurrence count over ALL distinct canonical
k-mers (`total_kmer_occurrences / distinct_kmers`). Unlike `median_solid_...` it
INCLUDES singleton error k-mers, so residual error deflates it toward 1.0; it decays
smoothly with `k` but its dynamic range between a well-covered backbone (reuse-oracle
~1.6 at k21) and a shattered one (~1.1) is too NARROW for a robust byte-identical
oracle. Retained as an interpretable cross-check (it cleanly shows the error-driven
decay) — but `select_reassembly_k` keys off the solid-backbone median, whose range
(~12 vs ~3) separates the regimes with margin. Returns `0.0` on empty input.
"""
function effective_kmer_coverage(reads, k::Int)::Float64
    counts = _kmer_count_spectrum(reads, k; canonical = true)
    isempty(counts) && return 0.0
    total_occurrences = sum(values(counts))
    distinct_kmers = length(counts)
    return total_occurrences / distinct_kmers
end

"""
    _genome_size_floor_k(reads, floor_k, ceiling_k; size_ref_k=17, solid_min=2,
                         occupancy_target=0.01) -> Int

Genome-size-aware UNIQUENESS floor (td-jt7r). The connectivity drop-down can, on a
very shattered read set, fall all the way back to `floor_k` (7); but `4^7 = 16384`
k-mer slots is SMALLER than a multi-kb genome, so its distinct genomic k-mers collide
and the de Bruijn graph tangles into k-mer-sized fragments. Raise the floor to the
smallest `k` whose k-mer space keeps genomic k-mers ~unique — `4^k >= genome_size /
occupancy_target` (default 1% occupancy ⇒ `4^k >= 100·G`). Genome size `G` is estimated
as the number of distinct SOLID canonical k-mers at `size_ref_k` (the genomic backbone;
singleton error k-mers, occurrence `< solid_min`, are excluded). Returns `floor_k` when
the estimate is empty/degenerate, and never exceeds `ceiling_k`.

Only the drop-down / fallback path is affected: on clean reads the connectivity
criterion returns the ceiling directly, so this floor is inert and the illumina path
stays byte-identical (oracle preservation).
"""
function _genome_size_floor_k(reads, floor_k::Int, ceiling_k::Int;
        size_ref_k::Int = 17, solid_min::Int = 2, occupancy_target::Float64 = 0.01)::Int
    ref_k = clamp(size_ref_k, floor_k, ceiling_k)
    counts = _kmer_count_spectrum(reads, ref_k; canonical = true)
    genome_size = count(>=(solid_min), values(counts))
    genome_size <= 0 && return floor_k
    k_unique = ceil(Int, log(genome_size / occupancy_target) / log(4))
    return clamp(k_unique, floor_k, ceiling_k)
end

# CAVEATS on the docstring below (td-jt7r, PR #444), attributed per claim:
#
#   * The CALIBRATION paragraph derives `connectivity_floor = 6.0` from
#     `median_solid_kmer_multiplicity`'s "~12 well-covered vs ~3–4 shattered"
#     separation. That separation rests on the backbone-coverage interpretation this
#     PR REFUTES (see the caveat above `median_solid_kmer_multiplicity`): the
#     statistic is an error-population statistic at higher coverage and is
#     non-monotone in coverage, so 6.0 is calibrated against a signal that does not
#     mean what the docstring says it means. The measurement lives in
#     `benchmarking/indel_toy_overcorrection_diagnostic.jl` (`solid_composition.csv`,
#     `k_ladder.csv`). Re-calibration is deferred to its own bead; the constant is
#     UNCHANGED here.
#   * Within that same paragraph, the sub-claims are of different kinds and should
#     not be read under one stamp: (i) "a floor of 2.0 fails to trip on the
#     shattered regime (3 > 2)" is arithmetic on the measured value and holds
#     whatever the value MEANS; (ii) "a floor near ~12 would clip the well-covered
#     regime" likewise; (iii) "6.0 sits between the two clusters with margin" is the
#     inference that inherits the refuted interpretation; (iv) "the measure MUST be
#     canonical" is an independent, unaffected claim.
#   * SIGNATURE: the `select_reassembly_k(...)` line below omits two real keyword
#     arguments — `genome_size_floor::Bool = true` and `size_ref_k::Int = 17`. Both
#     are ON by default and both materially change the result, because the
#     genome-size uniqueness floor (`_genome_size_floor_k`) can RAISE the effective
#     floor far above `floor_k = 7`.
#   * "If none qualifies, fall back to the smallest prime `>= floor_k`" is wrong in
#     three ways; see the corrected comment at the `return first(candidates)` site
#     at the bottom of this function.
"""
    select_reassembly_k(reads, ceiling_k; floor_k = 7, connectivity_floor = 6.0) -> Int

Choose a re-assembly k for corrected `reads`, bounded by `[floor_k, ceiling_k]`,
using a COVERAGE-AWARE CONNECTIVITY criterion. Iterating candidates DESCENDING from
`ceiling_k`, return the LARGEST candidate whose
`median_solid_kmer_multiplicity(reads, k) >= connectivity_floor` — the highest
specificity whose genomic backbone still recurs enough to keep the corrected-read de
Bruijn graph connected. If none qualifies, fall back to the smallest prime `>= floor_k`.

This measures backbone coverage `C·(1-e)^k` DIRECTLY (median of the solid k-mer peak)
rather than estimating `C` and `e` separately (avoids the circularity in a
residual-error → survival-model chain):

- Clean / well-corrected reads (Illumina): the solid backbone recurs ~coverage even
  at the ceiling (median-solid ~12 ≫ floor), so the ceiling is honored — byte-identical
  to legacy behavior, preserving the corrector's final-pass graph-reuse eligibility
  (which requires `reassembly_k == final_graph_k`).
- Shattered high-error long reads (nanopore): the backbone barely survives at high
  `k` (median-solid ~3 < floor), so `k` DROPS to a lower prime that stays connected.

The ceiling is the caller's explicit k (kept as-is even if non-prime, e.g. 21);
drop-down targets are primes (the `:scalable` k-ladder is prime), so the chosen k is
prime whenever it adapts and always lies in `[floor_k, ceiling_k]`.

CALIBRATION (td-jt7r): the connectivity floor is 6.0, NOT the 2.0 first proposed.
Measured backbone coverage sits at ~12 for a well-covered corrected regime and at
~3–4 for a shattered one; a floor of 2.0 fails to trip on the shattered regime (3 > 2)
so it never adapts, and a floor near ~12 would clip the well-covered regime. 6.0 sits
between the two clusters with margin on both sides (well-covered ÷2, shattered ×2).
Separately, the measure MUST be canonical (see `median_solid_kmer_multiplicity`) or
both-strand read sets spuriously adapt.
"""
function select_reassembly_k(
        reads,
        ceiling_k::Int;
        floor_k::Int = 7,
        connectivity_floor::Float64 = 6.0,
        genome_size_floor::Bool = true,
        size_ref_k::Int = 17
)::Int
    effective_floor = min(floor_k, ceiling_k)
    # Genome-size-aware uniqueness floor (td-jt7r): never let the drop-down fall below
    # the k whose k-mer space keeps genomic k-mers unique, so a multi-kb genome cannot
    # collapse to k-mer-sized fragments at the k=7 fallback. Inert on clean reads (the
    # ceiling is returned before any fallback), so the illumina path is byte-identical.
    if genome_size_floor
        effective_floor = min(
            max(effective_floor,
                _genome_size_floor_k(reads, effective_floor, ceiling_k;
                    size_ref_k = size_ref_k)),
            ceiling_k)
    end
    # The ceiling is the caller's explicit k (may be non-prime, e.g. 21) and is the
    # top candidate — honored UNCHANGED when clean/high-coverage reads keep it
    # connected, so clean/Illumina behavior stays byte-identical and graph-reuse
    # (which needs reassembly_k == final_graph_k) stays eligible. DROP-DOWN targets
    # are primes strictly below the ceiling (the :scalable k-ladder is prime).
    lower_primes = filter(<(ceiling_k), _dynamic_k_search_space(effective_floor, ceiling_k))
    candidates = vcat(lower_primes, ceiling_k)   # ascending; ceiling is the largest
    # Descend from the ceiling: the LARGEST candidate whose solid backbone still meets
    # the connectivity floor is the most specific k that keeps the graph joined.
    for k in Iterators.reverse(candidates)
        if median_solid_kmer_multiplicity(reads, k) >= connectivity_floor
            return k
        end
    end
    # Nothing qualifies (very low coverage / very high error): fall back to the
    # SMALLEST candidate, i.e. the most-connected k still in range. Three
    # corrections to the docstring's "smallest prime >= floor_k" wording (td-jt7r):
    #   1. The floor here is `effective_floor`, which the genome-size uniqueness
    #      floor above may have RAISED well past `floor_k`; the fallback is relative
    #      to that raised floor, not to `floor_k`.
    #   2. `candidates` is `lower_primes ++ [ceiling_k]`, so when no prime lies
    #      strictly below `ceiling_k` at or above `effective_floor`, `lower_primes`
    #      is empty and this returns `ceiling_k` — not a prime below it, and
    #      possibly not prime at all (the ceiling is the caller's explicit k).
    #   3. `_dynamic_k_search_space` has its own degenerate fallbacks: given an
    #      empty prime range it returns the LARGEST prime `<= max_k` (ignoring
    #      `min_k`), and `[1]` when no prime exists at all — so `first(candidates)`
    #      is not guaranteed to be `>= effective_floor` in those corners either.
    return first(candidates)
end
