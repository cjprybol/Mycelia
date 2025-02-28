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
    
    if canon_kmer == BioSequences.reverse_complement(qmer.kmer)
        K = length(qmer.qualities)
        reversed_quals = ntuple(i -> qmer.qualities[K - i + 1], K)
        return Qualmer(canon_kmer, reversed_quals)
    else
        return Qualmer(canon_kmer, qmer.qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
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
"""
function qualmers_unambiguous(sequence::BioSequences.LongSequence{BioSequences.AminoAcidAlphabet}, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    return qualmers_fw(sequence, quality, Val(K))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_canonical(sequence::BioSequences.LongSequence{BioSequences.DNAAlphabet{N}}, quality::AbstractVector{<:Integer}, ::Val{K}) where {N,K}
    return (
        (Qualmer(
            kmer,
            quality[pos:pos+K-1]
        ), pos) for (pos, kmer) in enumerate(Kmers.CanonicalKmers{BioSequences.DNAAlphabet{N}, K}(sequence))
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_canonical(sequence::BioSequences.LongSequence{BioSequences.AminoAcidAlphabet}, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    return qualmers_fw(sequence, quality, Val(K))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_fw(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_fw(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_canonical(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.FASTQ.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_canonical(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_unambiguous(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(FASTX.quality_scores(record))
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
    quality = collect(FASTX.quality_scores(record))
    return qualmers_unambiguous_canonical(sequence, quality, Val(k))
end