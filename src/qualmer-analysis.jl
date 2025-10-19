"""
$(DocStringExtensions.TYPEDEF)

A quality-aware k-mer that pairs sequence k-mers with their corresponding quality scores.

# Fields
$(DocStringExtensions.TYPEDFIELDS)

# Examples
```julia
using Mycelia
using Kmers

# Create a DNA k-mer with quality scores
kmer = DNAKmer("ATCG")
qualities = [30, 35, 32, 28]  # Phred quality scores
qmer = Qualmer(kmer, qualities)

# Access individual bases and qualities
base, quality = qmer[1]  # Returns ('A', 30)
```

This structure is essential for quality-aware sequence analysis and assembly algorithms
that consider both sequence content and sequencing quality information.
"""
struct Qualmer{KmerT, K}
    kmer::KmerT
    qualities::NTuple{K,UInt8}
    function Qualmer(kmer::KmerT, qualities::NTuple{K,UInt8}) where {KmerT<:Kmers.Kmer, K}
        @assert K == length(kmer) "Quality length must match kmer length"
        return new{KmerT,K}(kmer, qualities)
    end
end

function Qualmer(kmer::KmerT, qualities::AbstractVector{<:Integer}) where {KmerT<:Kmers.Kmer}
    @assert length(kmer) == length(qualities) "Quality length must match kmer length"
    k = length(kmer)
    Qualmer(kmer, NTuple{k, UInt8}(UInt8.(qualities)))
end

# Basic methods
Base.length(q::Qualmer) = length(q.kmer)
Base.getindex(q::Qualmer, i::Integer) = (q.kmer[i], q.qualities[i])
Base.:(==)(a::T, b::T) where {T<:Qualmer} = a.kmer == b.kmer && a.qualities == b.qualities
Base.hash(q::Qualmer, h::UInt) = hash(q.qualities, hash(q.kmer, h))

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the canonical representation of a DNA qualmer, considering both sequence and quality.
For DNA/RNA, this involves potentially reverse-complementing the k-mer and reversing
the quality scores accordingly.
"""
function canonical(qmer::Qualmer{KmerT}) where {KmerT<:Union{Kmers.DNAKmer, Kmers.RNAKmer}}
    canon_kmer = BioSequences.canonical(qmer.kmer)
    
    # If the canonical form is different from the original, it means we reverse-complemented
    # In that case, we need to reverse the quality scores
    if canon_kmer != qmer.kmer
        K = length(qmer.qualities)
        reversed_quals = ntuple(i -> qmer.qualities[K - i + 1], K)
        return Qualmer(canon_kmer, reversed_quals)
    else
        return Qualmer(canon_kmer, qmer.qualities)
    end
end

# Additional method for the Mer types from FwKmers
function canonical(qmer::Qualmer{KmerT}) where {K, KmerT<:Kmers.Mer{K, BioSequences.DNAAlphabet{N}}} where N
    canon_kmer = BioSequences.canonical(qmer.kmer)
    
    # If the canonical form is different from the original, it means we reverse-complemented
    # In that case, we need to reverse the quality scores
    if canon_kmer != qmer.kmer
        reversed_quals = ntuple(i -> qmer.qualities[K - i + 1], K)
        return Qualmer(canon_kmer, reversed_quals)
    else
        return Qualmer(canon_kmer, qmer.qualities)
    end
end

# Additional method for RNA Mer types
function canonical(qmer::Qualmer{KmerT}) where {K, KmerT<:Kmers.Mer{K, BioSequences.RNAAlphabet{N}}} where N
    canon_kmer = BioSequences.canonical(qmer.kmer)
    
    # If the canonical form is different from the original, it means we reverse-complemented
    # In that case, we need to reverse the quality scores
    if canon_kmer != qmer.kmer
        reversed_quals = ntuple(i -> qmer.qualities[K - i + 1], K)
        return Qualmer(canon_kmer, reversed_quals)
    else
        return Qualmer(canon_kmer, qmer.qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the canonical representation of an amino acid qualmer.
For amino acid sequences, the qualmer is returned unchanged since
amino acids do not have reverse complements like DNA/RNA.
"""
function canonical(qmer::Qualmer{<:Kmers.AAKmer})
    return qmer
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create an iterator that yields DNA qualmers from the given sequence and quality scores.
"""
function qualmers_fw(sequence::BioSequences.LongSequence{A}, quality::AbstractVector{<:Integer}, ::Val{K}) where {A,K}
    return (
        (Qualmer(
            kmer, 
            quality[pos:pos+K-1]
        ), pos) for (pos, kmer) in enumerate(Kmers.FwKmers{A, K}(sequence))
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate unambiguous DNA qualmers from sequence and quality scores.
Only includes k-mers that contain unambiguous nucleotides (A, T, G, C).
Skips k-mers containing ambiguous bases like N.

# Arguments
- `sequence`: DNA sequence as LongDNA
- `quality`: Vector of quality scores corresponding to each base
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples for each unambiguous k-mer
"""
function qualmers_unambiguous(sequence::BioSequences.LongSequence{BioSequences.DNAAlphabet{N}}, quality::AbstractVector{<:Integer}, ::Val{K}) where {N,K}
    return (
        (Qualmer(
            kmer,
            quality[pos:pos+K-1]
        ), pos) for (kmer, pos) in Kmers.UnambiguousDNAMers{K}(sequence)
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate unambiguous RNA qualmers from sequence and quality scores.
Only includes k-mers that contain unambiguous nucleotides (A, U, G, C).
Skips k-mers containing ambiguous bases like N.

# Arguments
- `sequence`: RNA sequence as LongRNA
- `quality`: Vector of quality scores corresponding to each base
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples for each unambiguous k-mer
"""
function qualmers_unambiguous(sequence::BioSequences.LongSequence{BioSequences.RNAAlphabet{N}}, quality::AbstractVector{<:Integer}, ::Val{K}) where {N,K}
    return (
        (Qualmer(
            kmer,
            NTuple{K,UInt8}(quality[pos:pos+K-1])
        ), pos) for (kmer, pos) in Kmers.UnambiguousRNAMers{K}(sequence)
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate amino acid qualmers from protein sequence and quality scores.
For amino acid sequences, all k-mers are considered unambiguous, so this
falls back to the forward qualmer generation.

# Arguments
- `sequence`: Protein sequence as LongAA
- `quality`: Vector of quality scores corresponding to each amino acid
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples for each k-mer
"""
function qualmers_unambiguous(sequence::BioSequences.LongSequence{BioSequences.AminoAcidAlphabet}, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    return qualmers_fw(sequence, quality, Val(K))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate canonical DNA qualmers from sequence and quality scores.
For each k-mer, returns the lexicographically smaller of the k-mer and its
reverse complement, ensuring consistent representation regardless of strand.

# Arguments
- `sequence`: DNA sequence as LongDNA
- `quality`: Vector of quality scores corresponding to each base
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples with canonical k-mer representation
"""
function qualmers_canonical(sequence::BioSequences.LongSequence{BioSequences.DNAAlphabet{N}}, quality::AbstractVector{<:Integer}, ::Val{K}) where {N,K}
    return (
        (canonical(Qualmer(
            kmer,
            quality[pos:pos+K-1]
        )), pos) for (pos, kmer) in enumerate(Kmers.FwKmers{BioSequences.DNAAlphabet{N}, K}(sequence))
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate canonical amino acid qualmers from protein sequence and quality scores.
For amino acid sequences, canonical representation is the same as forward
representation since proteins don't have reverse complements.

# Arguments
- `sequence`: Protein sequence as LongAA
- `quality`: Vector of quality scores corresponding to each amino acid
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples for each k-mer
"""
function qualmers_canonical(sequence::BioSequences.LongSequence{BioSequences.AminoAcidAlphabet}, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    return qualmers_fw(sequence, quality, Val(K))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate forward qualmers from a FASTQ record.
Extracts all k-mers in the forward direction with their associated
quality scores from the FASTQ record.

# Arguments
- `record`: FASTQ record containing sequence and quality data
- `k`: K-mer size

# Returns
- Iterator yielding (Qualmer, position) tuples for each k-mer
"""
function qualmers_fw(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_fw(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate canonical qualmers from a FASTQ record.
Extracts k-mers in canonical representation (lexicographically smaller
of forward and reverse complement) with their associated quality scores.

# Arguments
- `record`: FASTQ record containing sequence and quality data
- `k`: K-mer size

# Returns
- Iterator yielding (Qualmer, position) tuples with canonical k-mer representation
"""
function qualmers_canonical(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.FASTQ.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_canonical(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate unambiguous qualmers from a FASTQ record.
Extracts only k-mers that contain unambiguous nucleotides (no N bases)
with their associated quality scores.

# Arguments
- `record`: FASTQ record containing sequence and quality data
- `k`: K-mer size

# Returns
- Iterator yielding (Qualmer, position) tuples for each unambiguous k-mer
"""
function qualmers_unambiguous(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(UInt8, FASTX.quality_scores(record))
    return qualmers_unambiguous(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate unambiguous canonical qualmers from the given sequence and quality scores.
This function first extracts all unambiguous qualmers and then applies canonical 
representation to each one.
"""
# Also update functions that consume these results, like qualmers_unambiguous_canonical
function qualmers_unambiguous_canonical(sequence, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    # Return both canonical qualmer and position
    return ((canonical(qmer), pos) for (qmer, pos) in qualmers_unambiguous(sequence, quality, Val(K)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate unambiguous canonical qualmers from the given FASTQ record.
"""
function qualmers_unambiguous_canonical(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(UInt8, FASTX.quality_scores(record))
    return qualmers_unambiguous_canonical(sequence, quality, Val(k))
end

"""
Return the k-mer length encoded in a qualmer-supporting k-mer type.
"""
qualmer_kmer_length(::Type{Kmers.Kmer{A, K}}) where {A, K} = K

function _sequence_type_for_kmer(::Type{Kmers.Kmer{A, K}}) where {A, K}
    if A <: BioSequences.DNAAlphabet
        return BioSequences.LongDNA{4}
    elseif A <: BioSequences.RNAAlphabet
        return BioSequences.LongRNA{4}
    elseif A <: BioSequences.AminoAcidAlphabet
        return BioSequences.LongAA
    else
        throw(ArgumentError("Unsupported alphabet for qualmer counting: $(A)"))
    end
end

function _collect_quality_scores(record::FASTX.FASTQ.Record)
    return collect(UInt8, FASTX.quality_scores(record))
end

function _count_qualmers_impl(::Type{KMER_TYPE}, sequence, quality::AbstractVector{<:Integer}; canonical::Bool=false) where {KMER_TYPE<:Kmers.Kmer}
    k = qualmer_kmer_length(KMER_TYPE)
    iterator = canonical ? qualmers_unambiguous_canonical(sequence, quality, Val(k)) :
                           qualmers_unambiguous(sequence, quality, Val(k))
    return sort(StatsBase.countmap([qmer.kmer for (qmer, _) in iterator]))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count qualmers for a single FASTQ record.
"""
function count_qualmers(::Type{KMER_TYPE}, record::FASTX.FASTQ.Record) where {KMER_TYPE<:Kmers.Kmer}
    sequence_type = _sequence_type_for_kmer(KMER_TYPE)
    sequence = FASTX.sequence(sequence_type, record)
    quality = _collect_quality_scores(record)
    return _count_qualmers_impl(KMER_TYPE, sequence, quality; canonical=false)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count canonical qualmers for a single FASTQ record.
"""
function count_canonical_qualmers(::Type{KMER_TYPE}, record::FASTX.FASTQ.Record) where {KMER_TYPE<:Kmers.Kmer}
    sequence_type = _sequence_type_for_kmer(KMER_TYPE)
    sequence = FASTX.sequence(sequence_type, record)
    quality = _collect_quality_scores(record)
    return _count_qualmers_impl(KMER_TYPE, sequence, quality; canonical=true)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count qualmers for a collection of FASTQ records.
"""
function count_qualmers(::Type{KMER_TYPE}, records::AbstractVector{<:FASTX.FASTQ.Record}) where {KMER_TYPE<:Kmers.Kmer}
    combined = Dict{KMER_TYPE, Int}()
    for record in records
        for (kmer, count) in count_qualmers(KMER_TYPE, record)
            combined[kmer] = get(combined, kmer, 0) + count
        end
    end
    return sort(combined)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count canonical qualmers for a collection of FASTQ records.
"""
function count_canonical_qualmers(::Type{KMER_TYPE}, records::AbstractVector{<:FASTX.FASTQ.Record}) where {KMER_TYPE<:Kmers.Kmer}
    combined = Dict{KMER_TYPE, Int}()
    for record in records
        for (kmer, count) in count_canonical_qualmers(KMER_TYPE, record)
            combined[kmer] = get(combined, kmer, 0) + count
        end
    end
    return sort(combined)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count qualmers from a FASTQ reader.
"""
function count_qualmers(::Type{KMER_TYPE}, reader::FASTX.FASTQ.Reader) where {KMER_TYPE<:Kmers.Kmer}
    return count_qualmers(KMER_TYPE, collect(reader))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count canonical qualmers from a FASTQ reader.
"""
function count_canonical_qualmers(::Type{KMER_TYPE}, reader::FASTX.FASTQ.Reader) where {KMER_TYPE<:Kmers.Kmer}
    return count_canonical_qualmers(KMER_TYPE, collect(reader))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count qualmers in a FASTQ file.
"""
function count_qualmers(::Type{KMER_TYPE}, fastq_file::AbstractString) where {KMER_TYPE<:Kmers.Kmer}
    open(fastq_file) do io
        reader = FASTX.FASTQ.Reader(io)
        return count_qualmers(KMER_TYPE, reader)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count canonical qualmers in a FASTQ file.
"""
function count_canonical_qualmers(::Type{KMER_TYPE}, fastq_file::AbstractString) where {KMER_TYPE<:Kmers.Kmer}
    open(fastq_file) do io
        reader = FASTX.FASTQ.Reader(io)
        return count_canonical_qualmers(KMER_TYPE, reader)
    end
end

# ============================================================================
# Joint Probability Calculations for Quality-Aware Assembly
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert PHRED quality score to probability of correctness.
PHRED score Q relates to error probability P by: Q = -10 * log10(P)
Therefore, correctness probability = 1 - P = 1 - 10^(-Q/10)
"""
function phred_to_probability(phred_score::UInt8)
    # Convert PHRED score to probability that base call is correct
    # PHRED score Q = -10 * log10(P_error)
    # Therefore: P_error = 10^(-Q/10)
    # P_correct = 1 - P_error = 1 - 10^(-Q/10)
    error_prob = 10.0^(-Float64(phred_score) / 10.0)
    correct_prob = 1.0 - error_prob
    
    # Ensure we don't return negative probabilities
    return max(0.0, min(1.0, correct_prob))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert PHRED quality score to error probability.
"""
function phred_to_error_probability(phred_score::UInt8)
    return 10.0^(-Float64(phred_score) / 10.0)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate joint probability that a single qualmer is correct.
For a k-mer with quality scores [q1, q2, ..., qk], the joint probability
that all positions are correct is: ∏(1 - 10^(-qi/10))
"""
function qualmer_correctness_probability(qmer::Qualmer)
    log_prob = 0.0
    for q in qmer.qualities
        prob_correct = phred_to_probability(q)
        log_prob += log(prob_correct)
    end
    return exp(log_prob)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate joint probability that multiple observations of the same k-mer are all correct.
For independent observations of the same k-mer sequence, this represents our
confidence that the k-mer truly exists in the data.

# Arguments
- `qualmers`: Vector of Qualmer observations of the same k-mer sequence
- `use_log_space`: Use log-space arithmetic for numerical stability (default: true)

# Returns
- Float64: Joint probability that all observations are correct
"""
function joint_qualmer_probability(qualmers::Vector{<:Qualmer}; use_log_space::Bool=true)
    if isempty(qualmers)
        return 0.0
    end
    
    # Verify all qualmers have the same k-mer sequence
    kmer_seq = qualmers[1].kmer
    for qmer in qualmers[2:end]
        if qmer.kmer != kmer_seq
            throw(ArgumentError("All qualmers must have the same k-mer sequence"))
        end
    end
    
    if use_log_space
        log_joint_prob = 0.0
        for qmer in qualmers
            prob_correct = qualmer_correctness_probability(qmer)
            log_joint_prob += log(prob_correct)
        end
        return exp(log_joint_prob)
    else
        joint_prob = 1.0
        for qmer in qualmers
            joint_prob *= qualmer_correctness_probability(qmer)
        end
        return joint_prob
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate position-wise joint probability for multiple qualmer observations.
This is more sophisticated than `joint_qualmer_probability` as it considers
the quality at each position across all observations.

# Arguments
- `qualmers`: Vector of Qualmer observations of the same k-mer sequence
- `use_log_space`: Use log-space arithmetic for numerical stability (default: true)

# Returns
- Float64: Position-wise joint probability
"""
function position_wise_joint_probability(qualmers::Vector{<:Qualmer}; use_log_space::Bool=true)
    if isempty(qualmers)
        return 0.0
    end
    
    k = length(qualmers[1])
    
    # Verify all qualmers have the same length and k-mer sequence
    kmer_seq = qualmers[1].kmer
    for qmer in qualmers[2:end]
        if length(qmer) != k || qmer.kmer != kmer_seq
            throw(ArgumentError("All qualmers must have the same k-mer sequence and length"))
        end
    end
    
    if use_log_space
        log_joint_prob = 0.0
        for pos in 1:k
            # Calculate joint probability for this position across all observations
            log_pos_prob = 0.0
            for qmer in qualmers
                prob_correct = phred_to_probability(qmer.qualities[pos])
                if prob_correct <= 0.0
                    @warn "Invalid probability $prob_correct for quality $(qmer.qualities[pos]). Using minimum probability."
                    prob_correct = 1e-10  # Use a very small but positive probability
                end
                log_pos_prob += log(prob_correct)
            end
            log_joint_prob += log_pos_prob
        end
        return exp(log_joint_prob)
    else
        joint_prob = 1.0
        for pos in 1:k
            pos_prob = 1.0
            for qmer in qualmers
                pos_prob *= phred_to_probability(qmer.qualities[pos])
            end
            joint_prob *= pos_prob
        end
        return joint_prob
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the quality score for a single base given multiple observations.

This function implements the "Converting to Error Probabilities and Combining" method:
1. Takes error probabilities from multiple reads covering the same base
2. Calculates probability of ALL reads being wrong by multiplying probabilities
3. Calculates final Phred score from this combined probability

To avoid numerical underflow with very small probabilities, the calculation
is performed in log space.

# Arguments
- `error_probabilities::Vector{Float64}`: Vector of error probabilities from 
  multiple reads covering the same base position

# Returns
- `Float64`: Phred quality score representing the combined confidence
"""
function joint_base_quality_score(error_probabilities::Vector{Float64})
    if isempty(error_probabilities)
        return 0.0  # No data available
    end
    
    # Work in log space to avoid underflow
    log_p_all_wrong = sum(log.(error_probabilities))
    
    # Convert back from log space
    p_all_wrong = exp(log_p_all_wrong)
    
    # Prevent underflow/overflow
    if p_all_wrong <= eps(Float64)
        return 999.0  # Cap at a very high quality score
    end
    
    # Convert to Phred score
    return Mycelia.error_rate_to_q_value(p_all_wrong)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate kmer quality score using the specified aggregation method.

Available methods:
- `:min`: Use the minimum base quality (default)
- `:mean`: Use the mean of all base qualities
- `:geometric`: Use the geometric mean (appropriate for probabilities)
- `:harmonic`: Use the harmonic mean (emphasizes lower values)

# Arguments
- `base_qualities::Vector{Float64}`: Vector of quality scores for each base
- `method::Symbol`: Method to use for aggregation

# Returns
- `Float64`: Overall quality score for the kmer
"""
function kmer_quality_score(base_qualities::Vector{Float64}, method::Symbol=:min)
    if isempty(base_qualities)
        return 0.0
    end
    
    if method == :min
        return minimum(base_qualities)
    elseif method == :mean
        return mean(base_qualities)
    elseif method == :geometric
        # Convert to probabilities, compute geometric mean, convert back
        error_probs = [phred_to_error_prob(q) for q in base_qualities]
        geo_mean_prob = exp(sum(log.(error_probs)) / length(error_probs))
        return error_prob_to_phred(geo_mean_prob)
    elseif method == :harmonic
        # Harmonic mean emphasizes lower values
        # which is appropriate for quality scores
        return length(base_qualities) / sum(1.0 ./ base_qualities)
    else
        error("Unknown aggregation method: $method")
    end
end