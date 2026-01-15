"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return the alphabet associated with a `BioSequence` type.

# Arguments
- `s::BioSequences.BioSequence`: A subtype instance.

# Returns
`BioSymbols.Alphabet` of the sequence type.
"""
function get_biosequence_alphabet(s::T) where T<:BioSequences.BioSequence
    return first(T.parameters)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determines the most likely alphabet of a biological sequence.

The function follows a hierarchical approach, from the most specific (unambiguous)
to the most general (ambiguous) alphabets. It checks if the set of characters
in the input sequence is a subset of the character sets defined in `Mycelia.alphabets`.

The order of checking is:
1.  Unambiguous DNA
2.  Unambiguous RNA
3.  DNA + N (ACGTN)
4.  RNA + N (ACGUN)
5.  Unambiguous Amino Acid
6.  Ambiguous DNA (if not a fit for any previous set)
7.  Ambiguous RNA
8.  Ambiguous Amino Acid

If a sequence fits multiple alphabets (e.g., "ACGU" is valid RNA and AA), it is
classified as the first one matched in the hierarchy (RNA in this case). If the
sequence does not fit any of the defined alphabets or is empty, an `ArgumentError` is thrown.

Warning: This is a fast single-sequence heuristic and can misclassify short or
ambiguous reads (for example, short amino-acid sequences composed only of A/C/G/T).
Prefer multi-read inference or explicit hints when the input alphabet is ambiguous.
"""
function detect_alphabet(seq::AbstractString)::Symbol
    if isempty(seq)
        throw(ArgumentError("Input sequence cannot be empty."))
    end

    seq_chars = Set(seq)

    # --- Step 1: Check against Unambiguous Alphabets ---
    # From smallest to largest to find the most specific classification.

    # DNA is the most restrictive unambiguous nucleotide alphabet.
    if issubset(seq_chars, UNAMBIGUOUS_DNA_CHARSET)
        return :DNA
    end
    # RNA is the next most restrictive unambiguous nucleotide alphabet.
    if issubset(seq_chars, UNAMBIGUOUS_RNA_CHARSET)
        return :RNA
    end

    # --- Step 2: Check against N-containing Alphabets ---
    # Check DNA+N and RNA+N before full ambiguous sets.

    if issubset(seq_chars, ACGTN_DNA_CHARSET)
        return :DNA
    end
    if issubset(seq_chars, ACGUN_RNA_CHARSET)
        return :RNA
    end

    # --- Step 3: Check Unambiguous Amino Acid ---
    # Amino Acid is the least restrictive unambiguous alphabet.
    if issubset(seq_chars, UNAMBIGUOUS_AA_CHARSET)
        return :AA
    end

    # --- Step 4: Check against Full Ambiguous Alphabets ---
    # This block is reached only if the sequence contains multiple ambiguous characters.

    if issubset(seq_chars, AMBIGUOUS_DNA_CHARSET)
        return :DNA
    end
    if issubset(seq_chars, AMBIGUOUS_RNA_CHARSET)
        return :RNA
    end
    if issubset(seq_chars, AMBIGUOUS_AA_CHARSET)
        return :AA
    end

    # --- Step 5: No Alphabet Found ---
    # If the sequence characters do not form a subset of any known alphabet.
    unmatched_chars = setdiff(seq_chars, union(AMBIGUOUS_DNA_CHARSET, AMBIGUOUS_RNA_CHARSET, AMBIGUOUS_AA_CHARSET))
    throw(ArgumentError("Sequence contains characters that do not belong to any known alphabet: $(join(unmatched_chars, ", "))"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Detect the alphabet of a `LongDNA` sequence.

Always returns `:DNA`.
"""
function detect_alphabet(sequence::BioSequences.LongDNA)
    return :DNA
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Detect the alphabet of a `LongRNA` sequence.

Always returns `:RNA`.
"""
function detect_alphabet(sequence::BioSequences.LongRNA)
    return :RNA
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Detect the alphabet of a `LongAA` sequence.

Always returns `:AA`.
"""
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Validates that a sequence string can be successfully converted to the specified alphabet type.

Uses BioSequences constructors to validate the sequence without requiring alphabet detection.
Returns `true` if the sequence is valid for the specified alphabet, `false` otherwise.

# Arguments
- `seq::AbstractString`: The sequence string to validate
- `alphabet::Symbol`: The target alphabet (`:DNA`, `:RNA`, or `:AA`)

# Returns
`Bool`: `true` if sequence is valid for the alphabet, `false` otherwise.

# Examples
```julia
validate_alphabet("ACGT", :DNA)   # true
validate_alphabet("ACGU", :RNA)   # true
validate_alphabet("ACGT", :RNA)   # false (contains T)
validate_alphabet("ACDEFGHIKLMNPQRSTVWY", :AA)  # true
```
"""
function validate_alphabet(seq::AbstractString, alphabet::Symbol)::Bool
    try
        if alphabet == :DNA
            BioSequences.LongDNA{4}(seq)
        elseif alphabet == :RNA
            BioSequences.LongRNA{4}(seq)
        elseif alphabet == :AA
            BioSequences.LongAA(seq)
        else
            return false
        end
        return true
    catch
        return false
    end
end
