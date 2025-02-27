"""
# Qualmer Analysis

This module provides functionality for analyzing sequences (DNA, RNA, or amino acids) 
along with their quality scores from FASTQ files. The central concept is a "qualmer" - 
which pairs a k-mer (sequence of length K) with its corresponding quality scores.

Key features:
- Qualmer types for DNA, RNA, and protein sequences
- Iterators for processing all qualmers in a sequence
- Functions for canonical representation of qualmers
- Utilities for working with FASTQ records

## Basic Usage
```julia
import Mycelia
import FASTQ

# Create a DNAQualmer from a FASTQ record
record = FASTQ.Record()
# ... load record from file
qmer = Mycelia.qualmer(record, Mycelia.DNAQualmer{21,4})  # 21-mer with 4-bit encoding

# Process all qualmers in a record
for qmer in Mycelia.EveryQualmer(record, Mycelia.DNAQualmer{21,4})
    # Process each qualmer
end

# Get canonical representations of all qualmers in a record
canonical_qmers = Mycelia.canonical_qualmers(record, Mycelia.DNAQualmer{21,4})
```
"""

"""
    AbstractQualmer{K,N}

Abstract base type for all qualmer types, which combine a k-mer of length `K` 
with its quality scores. The type parameter `N` refers to the number of bits used
to encode each nucleotide or amino acid.
"""
abstract type AbstractQualmer{K,N} end

"""
    DNAQualmer{K,N} <: AbstractQualmer{K,N}

Represents a DNA k-mer of length `K` with corresponding quality scores.
The type parameter `N` specifies the bit encoding used for DNA (typically 4 or 2).

# Fields
- `kmer::Kmers.DNAKmer{K,N}`: The DNA k-mer sequence
- `qualities::NTuple{K,UInt8}`: Quality scores for each position in the k-mer

# Examples
```julia
# Create a DNAQualmer from a DNA sequence and quality scores
seq = BioSequences.LongDNA{4}("ACGTACGT")
quality = UInt8[30, 32, 25, 40, 38, 30, 25, 20]
qmer = qualmer(seq, quality, DNAQualmer{5,4})
```
"""
struct DNAQualmer{K,N} <: AbstractQualmer{K,N}
    kmer::Kmers.DNAKmer{K,N}
    qualities::NTuple{K,UInt8}  # Fixed-size tuple for better performance
    
    # Inner constructor to validate
    function DNAQualmer{K,N}(kmer::Kmers.DNAKmer{K,N}, qualities::NTuple{K,UInt8}) where {K,N}
        @assert length(qualities) == K "Quality length must match kmer length"
        new{K,N}(kmer, qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function DNAQualmer(kmer::Kmers.DNAKmer{K,N}, qualities::NTuple{K,UInt8}) where {K,N}
    DNAQualmer{K,N}(kmer, qualities)
end

# Similar constructors for RNA and AA qualmers
struct RNAQualmer{K,N} <: AbstractQualmer{K,N}
    kmer::Kmers.RNAKmer{K,N}
    qualities::NTuple{K,UInt8}  # Fixed-size tuple for better performance
    
    # Inner constructor to validate
    function RNAQualmer{K,N}(kmer::Kmers.RNAKmer{K,N}, qualities::NTuple{K,UInt8}) where {K,N}
        @assert length(qualities) == K "Quality length must match kmer length"
        new{K,N}(kmer, qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function RNAQualmer(kmer::Kmers.RNAKmer{K,N}, qualities::NTuple{K,UInt8}) where {K,N}
    RNAQualmer{K,N}(kmer, qualities)
end

struct AAQualmer{K,N} <: AbstractQualmer{K,N}
    kmer::Kmers.AAKmer{K,N}
    qualities::NTuple{K,UInt8}  # Fixed-size tuple for better performance
    
    # Inner constructor to validate
    function AAQualmer{K,N}(kmer::Kmers.AAKmer{K,N}, qualities::NTuple{K,UInt8}) where {K,N}
        @assert length(qualities) == K "Quality length must match kmer length"
        new{K,N}(kmer, qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function AAQualmer(kmer::Kmers.AAKmer{K,N}, qualities::NTuple{K,UInt8}) where {K,N}
    AAQualmer{K,N}(kmer, qualities)
end

# Basic methods
Base.length(q::AbstractQualmer) = K
Base.getindex(q::AbstractQualmer, i::Integer) = (q.kmer[i], q.qualities[i])
Base.:(==)(a::T, b::T) where {T<:AbstractQualmer} = a.kmer == b.kmer && a.qualities == b.qualities
Base.hash(q::AbstractQualmer, h::UInt) = hash(q.qualities, hash(q.kmer, h))

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmer(seq::BioSequences.LongDNA{N}, quality::Vector{UInt8}, ::Type{DNAQualmer{K,N}}) where {K,N}
    if length(seq) < K
        throw(ArgumentError("Sequence length must be at least K"))
    end
    if length(seq) != length(quality)
        throw(ArgumentError("Sequence and quality must have same length"))
    end
    
    kmer = Kmers.DNAKmer{K,N}(seq[1:K])
    qual_tuple = NTuple{K,UInt8}(quality[1:K])
    DNAQualmer{K,N}(kmer, qual_tuple)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmer(seq::BioSequences.LongRNA{N}, quality::Vector{UInt8}, ::Type{RNAQualmer{K,N}}) where {K,N}
    if length(seq) < K
        throw(ArgumentError("Sequence length must be at least K"))
    end
    if length(seq) != length(quality)
        throw(ArgumentError("Sequence and quality must have same length"))
    end
    
    kmer = Kmers.RNAKmer{K,N}(seq[1:K])
    qual_tuple = NTuple{K,UInt8}(quality[1:K])
    RNAQualmer{K,N}(kmer, qual_tuple)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmer(seq::BioSequences.LongAA{N}, quality::Vector{UInt8}, ::Type{AAQualmer{K,N}}) where {K,N}
    if length(seq) < K
        throw(ArgumentError("Sequence length must be at least K"))
    end
    if length(seq) != length(quality)
        throw(ArgumentError("Sequence and quality must have same length"))
    end
    
    kmer = Kmers.AAKmer{K,N}(seq[1:K])
    qual_tuple = NTuple{K,UInt8}(quality[1:K])
    AAQualmer{K,N}(kmer, qual_tuple)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the canonical representation of a DNA qualmer, considering both sequence and quality.
For DNA/RNA, this involves potentially reverse-complementing the k-mer and reversing
the quality scores accordingly.

# Arguments
- `qmer::DNAQualmer{K,N}`: The input qualmer

# Returns
- `DNAQualmer{K,N}`: The canonical qualmer

# Notes
For DNA and RNA, if the canonical k-mer is the reverse complement of the original k-mer,
the quality scores are also reversed to maintain position correspondence.
"""
# Get canonical qualmer (considering both sequence and quality)
function canonical(qmer::DNAQualmer{K,N}) where {K,N}
    canon_kmer = BioSequences.canonical(qmer.kmer)
    
    if canon_kmer == BioSequences.reverse_complement(qmer.kmer)
        reversed_quals = NTuple{K,UInt8}(reverse(collect(qmer.qualities)))
        return DNAQualmer{K,N}(canon_kmer, reversed_quals)
    else
        return DNAQualmer{K,N}(canon_kmer, qmer.qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function canonical(qmer::RNAQualmer{K,N}) where {K,N}
    canon_kmer = BioSequences.canonical(qmer.kmer)
    
    if canon_kmer == BioSequences.reverse_complement(qmer.kmer)
        reversed_quals = NTuple{K,UInt8}(reverse(collect(qmer.qualities)))
        return RNAQualmer{K,N}(canon_kmer, reversed_quals)
    else
        return RNAQualmer{K,N}(canon_kmer, qmer.qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function canonical(qmer::AAQualmer{K,N}) where {K,N}
    return qmer
end

# Iterator for qualmers in a sequence
struct EveryQualmer{Q<:AbstractQualmer{K,N}, S<:BioSequence, T<:Integer} where {K,N}
    seq::S
    quality::Vector{T}
    
    function EveryQualmer{Q}(seq::S, quality::Vector{T}) where {Q<:AbstractQualmer{K,N}, S<:BioSequence, T<:Integer, K, N}
        if length(seq) != length(quality)
            throw(ArgumentError("Sequence and quality must have the same length"))
        end
        new{Q,S,T}(seq, quality)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function EveryQualmer(seq::S, quality::Vector{T}, ::Type{Q}) where {S<:BioSequence, T<:Integer, Q<:AbstractQualmer{K,N}} where {K,N}
    EveryQualmer{Q}(seq, quality)
end

# Implement iterator interface
function Base.length(iter::EveryQualmer{Q}) where {Q<:AbstractQualmer{K,N}} where {K,N}
    max(0, length(iter.seq) - K + 1)
end

function Base.eltype(::Type{EveryQualmer{Q}}) where {Q}
    return Q
end

function Base.iterate(iter::EveryQualmer{Q,S,T}) where {Q<:AbstractQualmer{K,N}, S, T} where {K,N}
    if length(iter.seq) < K
        return nothing
    end
    
    first_kmer = _extract_kmer(iter.seq[1:K], Q)
    first_qual = NTuple{K,UInt8}(iter.quality[1:K])
    
    return _construct_qualmer(first_kmer, first_qual, Q), 2
end

function Base.iterate(iter::EveryQualmer{Q,S,T}, state) where {Q<:AbstractQualmer{K,N}, S, T} where {K,N}
    if state + K - 1 > length(iter.seq)
        return nothing
    end
    
    next_kmer = _extract_kmer(iter.seq[state:state+K-1], Q)
    next_qual = NTuple{K,UInt8}(iter.quality[state:state+K-1])
    
    return _construct_qualmer(next_kmer, next_qual, Q), state + 1
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function _extract_kmer(subseq, ::Type{DNAQualmer{K,N}}) where {K,N}
    return Kmers.DNAKmer{K,N}(subseq)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function _extract_kmer(subseq, ::Type{RNAQualmer{K,N}}) where {K,N}
    return Kmers.RNAKmer{K,N}(subseq)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function _extract_kmer(subseq, ::Type{AAQualmer{K,N}}) where {K,N}
    return Kmers.AAKmer{K,N}(subseq)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function _construct_qualmer(kmer::Kmers.DNAKmer{K,N}, qual::NTuple{K,UInt8}, ::Type{DNAQualmer{K,N}}) where {K,N}
    return DNAQualmer{K,N}(kmer, qual)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function _construct_qualmer(kmer::Kmers.RNAKmer{K,N}, qual::NTuple{K,UInt8}, ::Type{RNAQualmer{K,N}}) where {K,N}
    return RNAQualmer{K,N}(kmer, qual)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function _construct_qualmer(kmer::Kmers.AAKmer{K,N}, qual::NTuple{K,UInt8}, ::Type{AAQualmer{K,N}}) where {K,N}
    return AAQualmer{K,N}(kmer, qual)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmer(rec::FASTQ.Record, ::Type{DNAQualmer{K,N}}) where {K,N}
    seq = FASTQ.sequence(BioSequences.LongDNA{N}, rec)
    quality = FASTQ.quality(rec)
    qualmer(seq, quality, DNAQualmer{K,N})
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmer(rec::FASTQ.Record, ::Type{RNAQualmer{K,N}}) where {K,N}
    seq = FASTQ.sequence(BioSequences.LongRNA{N}, rec)
    quality = FASTQ.quality(rec)
    qualmer(seq, quality, RNAQualmer{K,N})
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmer(rec::FASTQ.Record, ::Type{AAQualmer{K,N}}) where {K,N}
    seq = FASTQ.sequence(BioSequences.LongAA{N}, rec)
    quality = FASTQ.quality(rec)
    qualmer(seq, quality, AAQualmer{K,N})
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function EveryQualmer(rec::FASTQ.Record, ::Type{Q}) where {Q<:AbstractQualmer{K,N}} where {K,N}
    seq_type = if Q <: DNAQualmer
        BioSequences.LongDNA{N}
    elseif Q <: RNAQualmer
        BioSequences.LongRNA{N}
    else
        BioSequences.LongAA{N}
    end
    seq = FASTQ.sequence(seq_type, rec)
    quality = FASTQ.quality(rec)
    EveryQualmer{Q}(seq, quality)
end

# function process_qualmers(fastq_file::String, ::Type{Q}, processor::Function) where {Q<:AbstractQualmer{K,N}} where {K,N}
#     open(FASTQ.Reader, fastq_file) do reader
#         for record in reader
#             for qmer in EveryQualmer(record, Q)
#                 processor(qmer)
#             end
#         end
#     end
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract all qualmers from a FASTQ record and convert them to their canonical form.

# Arguments
- `rec::FASTQ.Record`: The FASTQ record to process
- `::Type{Q}`: Type parameter specifying the qualmer type

# Returns
- `Vector{Q}`: Vector of canonical qualmers from the record
"""
function canonical_qualmers(rec::FASTQ.Record, ::Type{Q}) where {Q<:AbstractQualmer{K,N}} where {K,N}
    [canonical(qmer) for qmer in EveryQualmer(rec, Q)]
end