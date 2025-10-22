"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert an amino acid sequence to a reduced alphabet representation.

Reduced amino acid alphabets group the 20 standard amino acids into fewer classes
based on physicochemical properties (hydrophobicity, charge, size, aromaticity).
This simplifies protein sequence analysis and can improve computational efficiency
while reducing noise in machine learning models.

# Arguments
- `sequence::BioSequences.LongAA`: Input amino acid sequence
- `scheme::Symbol`: Reduced alphabet scheme to use (see available schemes below)

# Returns
- `String`: Sequence converted to reduced alphabet characters

# Available Schemes
- `:HP2` - Binary hydrophobic (H) / polar (P) classification (2 classes)
- `:HYDROPATHY3` - Hydrophobic (H) / neutral (N) / polar (P) (3 classes)
- `:GBMR4` - Glycine (G) / Proline (P) / Hydrophobic (H) / Basic-polar (B) (4 classes)
- `:CHEMICAL5` - Aliphatic (A) / aRomatic (R) / Charged (C) / Tiny (T) / Diverse (D) (5 classes)
- `:CHEMICAL6` - Aliphatic (A) / aRomatic (R) / Positive (+) / Negative (-) / Tiny (T) / Diverse (D) (6 classes)
- `:SDM12` - Structure-dependent model with 12 groups (Murphy et al. 2000)

# Examples
```julia
## Basic usage with different schemes
seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")

## Binary hydrophobic/polar
reduced = reduce_amino_acid_alphabet(seq, :HP2)
## Returns: "HHPPPPHHPPPPPPPPPHPP"

## Three-class hydropathy
reduced = reduce_amino_acid_alphabet(seq, :HYDROPATHY3)
## Returns: "HHPPPPNHPPNNPPPPNHNH"

## Get information about a scheme
info = get_reduced_alphabet_info(:HP2)
println(info[:description])
## Prints: "Simplest reduction grouping amino acids by hydrophobicity (H) vs polarity (P)"

## List all available schemes
schemes = list_reduced_alphabets()
## Returns: [:HP2, :HYDROPATHY3, :GBMR4, :CHEMICAL5, :CHEMICAL6, :SDM12]
```

# References
- Murphy et al. (2000) Protein Eng. 13(3):149-152
- Peterson et al. (2009) BMC Bioinformatics 10:228
- Zheng et al. (2019) Database (Oxford) baz131 (RAACBook)

# See Also
- `list_reduced_alphabets()`: Get list of available schemes
- `get_reduced_alphabet_info()`: Get detailed information about a scheme
"""
function reduce_amino_acid_alphabet(sequence::BioSequences.LongAA, scheme::Symbol)::String
    ## Check if the scheme exists
    if !haskey(Mycelia.REDUCED_ALPHABETS, scheme)
        available = join(sort(collect(keys(Mycelia.REDUCED_ALPHABETS))), ", :")
        throw(ArgumentError("Unknown reduced alphabet scheme: $scheme. Available schemes: :$available"))
    end

    ## Get the mapping dictionary for this scheme
    alphabet_map = Mycelia.REDUCED_ALPHABETS[scheme]

    ## Convert each amino acid to its reduced representation
    reduced_chars = Char[]
    for aa in sequence
        if !haskey(alphabet_map, aa)
            ## Handle unexpected amino acids (ambiguous, gap, termination)
            throw(ArgumentError("Amino acid $aa not found in $scheme reduction scheme. Only standard amino acids are supported."))
        end
        push!(reduced_chars, alphabet_map[aa])
    end

    return String(reduced_chars)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get a list of all available reduced amino acid alphabet schemes.

# Returns
- `Vector{Symbol}`: Vector of scheme names that can be used with `reduce_amino_acid_alphabet()`

# Examples
```julia
schemes = list_reduced_alphabets()
## Returns: [:HP2, :HYDROPATHY3, :GBMR4, :CHEMICAL5, :CHEMICAL6, :SDM12]

## Use with reduce_amino_acid_alphabet
seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
for scheme in schemes
    reduced = reduce_amino_acid_alphabet(seq, scheme)
    println("$scheme: $reduced")
end
```
"""
function list_reduced_alphabets()::Vector{Symbol}
    return sort(collect(keys(Mycelia.REDUCED_ALPHABETS)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get detailed information about a specific reduced alphabet scheme.

# Arguments
- `scheme::Symbol`: The reduced alphabet scheme to query

# Returns
- `Dict{Symbol, Any}`: Dictionary containing:
  - `:name` - Full name of the scheme
  - `:classes` - Number of classes/groups
  - `:description` - Description of the scheme
  - `:groups` - Mapping of reduced alphabet characters to amino acid groups

# Examples
```julia
## Get information about the HP2 scheme
info = get_reduced_alphabet_info(:HP2)
println(info[:name])          ## "Binary Hydrophobic-Polar"
println(info[:classes])       ## 2
println(info[:description])   ## "Simplest reduction grouping amino acids by hydrophobicity (H) vs polarity (P)"
println(info[:groups])        ## Dict('H' => "ACFILMVW", 'P' => "GTSYPNDEQKRH")

## Print all available schemes with their descriptions
for scheme in list_reduced_alphabets()
    info = get_reduced_alphabet_info(scheme)
    println("\$scheme (\$(info[:classes]) classes): \$(info[:description])")
end
```
"""
function get_reduced_alphabet_info(scheme::Symbol)::Dict{Symbol, Any}
    if !haskey(Mycelia.REDUCED_ALPHABET_INFO, scheme)
        available = join(sort(collect(keys(Mycelia.REDUCED_ALPHABET_INFO))), ", :")
        throw(ArgumentError("Unknown reduced alphabet scheme: $scheme. Available schemes: :$available"))
    end

    return Mycelia.REDUCED_ALPHABET_INFO[scheme]
end
