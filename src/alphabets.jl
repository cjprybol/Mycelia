function get_biosequence_alphabet(s::T) where T<:BioSequences.BioSequence
    return first(T.parameters)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determines the alphabet of a sequence. The function scans through `seq` only once:
- If a 'T' or 't' is found (and no 'U/u'), the sequence is classified as DNA.
- If a 'U' or 'u' is found (and no 'T/t'), it is classified as RNA.
- If both T and U occur, an error is thrown.
- If a character outside the canonical nucleotide and ambiguity codes is encountered,
  the sequence is assumed to be protein.
- If neither T nor U are found, the sequence is assumed to be DNA.
"""
function detect_alphabet(seq::AbstractString)::Symbol
    hasT = false
    hasU = false
    # Define allowed nucleotide characters (both for DNA and RNA, including common ambiguity codes)
    # TODO: define this by merging the alphabets from BioSymbols
    valid_nucleotides = "ACGTacgtACGUacguNRYSWKMBDHnryswkmbdh"
    for c in seq
        if c == 'T' || c == 't'
            hasT = true
        elseif c == 'U' || c == 'u'
            hasU = true
        elseif !(c in valid_nucleotides)
            # If an unexpected character is encountered, assume it's a protein sequence.
            return :AA
        end
        if hasT && hasU
            throw(ArgumentError("Sequence contains both T and U, ambiguous alphabet"))
        end
    end
    if hasT
        return :DNA
    elseif hasU
        return :RNA
    else
        # In the absence of explicit T or U, default to DNA.
        return :DNA
    end
end

# Specific dispatch for DNA sequences
function detect_alphabet(sequence::BioSequences.LongDNA)
    return :DNA
end

# Specific dispatch for RNA sequences
function detect_alphabet(sequence::BioSequences.LongRNA)
    return :RNA
end

# Specific dispatch for protein/amino acid sequences
function detect_alphabet(sequence::BioSequences.LongAA)
    return :AA
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Converts the given sequence (output from FASTX.sequence) into the appropriate BioSequence type:
- DNA sequences are converted using `BioSequences.LongDNA`
- RNA sequences are converted using `BioSequences.LongRNA`
- AA sequences are converted using `BioSequences.LongAA`
"""
function convert_sequence(seq::AbstractString)
    alphabet = detect_alphabet(seq)
    if alphabet == :DNA
        return BioSequences.LongDNA{4}(seq)
    elseif alphabet == :RNA
        return BioSequences.LongRNA{4}(seq)
    elseif alphabet == :AA
        return BioSequences.LongAA(seq)
    else
        throw(ArgumentError("Unrecognized alphabet type"))
    end
end