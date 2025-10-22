# const METADATA = joinpath(dirname(dirname(pathof(Mycelia))), "docs", "metadata")

const DEFAULT_BLASTDB_PATH = "$(homedir())/workspace/blastdb"

# fix new error
# ENV["MAMBA_ROOT_PREFIX"] = joinpath(DEPOT_PATH[1], "conda", "3", "x86_64")

# Mycelia.NERSC_MEM * .95
# const NERSC_MEM=512
# const NERSC_MEM=480
const NERSC_MEM=460
# const NERSC_CPU=240
const NERSC_CPU=240

# Phase these out and move to the below more specific options
const DNA_ALPHABET = BioSymbols.ACGT
const RNA_ALPHABET = BioSymbols.ACGU
const AA_ALPHABET = filter(
    x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x) || BioSymbols.isterm(x)),
    BioSymbols.alphabet(BioSymbols.AminoAcid))

# Old version that included termination characters (not valid for FASTQ):
# const AA_ALPHABET = filter(
#     x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x)),
#     BioSymbols.alphabet(BioSymbols.AminoAcid))


# Helper function to convert a tuple of BioSymbols to a case-insensitive Set of Chars
function symbols_to_char_set(symbols)
    chars = Char.(symbols)
    # Return a set containing both upper and lower case versions of the characters
    return Set([chars..., lowercase.(chars)...])
end

# --- ACGTN and ACGUN Alphabets (including N ambiguous character) ---
# These are intermediate alphabets that include the most common ambiguous character (N)
# but avoid the full set of ambiguous characters for more efficient detection.

const ACGTN_DNA_SYMBOLS = BioSymbols.ACGTN
const ACGUN_RNA_SYMBOLS = BioSymbols.ACGUN

const ACGTN_DNA_CHARSET = symbols_to_char_set(ACGTN_DNA_SYMBOLS)
const ACGUN_RNA_CHARSET = symbols_to_char_set(ACGUN_RNA_SYMBOLS)

# --- Unambiguous (Canonical) Alphabets ---
# Filtered to exclude ambiguous symbols and gaps, representing the core characters.

const UNAMBIGUOUS_DNA_SYMBOLS = filter(s -> !BioSymbols.isambiguous(s) && !BioSymbols.isgap(s), BioSymbols.alphabet(BioSequences.DNA))
const UNAMBIGUOUS_RNA_SYMBOLS = filter(s -> !BioSymbols.isambiguous(s) && !BioSymbols.isgap(s), BioSymbols.alphabet(BioSymbols.RNA))
const UNAMBIGUOUS_AA_SYMBOLS  = filter(s -> !BioSymbols.isambiguous(s) && !BioSymbols.isgap(s), BioSymbols.alphabet(BioSymbols.AminoAcid))

const UNAMBIGUOUS_DNA_CHARSET = symbols_to_char_set(UNAMBIGUOUS_DNA_SYMBOLS)
const UNAMBIGUOUS_RNA_CHARSET = symbols_to_char_set(UNAMBIGUOUS_RNA_SYMBOLS)
const UNAMBIGUOUS_AA_CHARSET  = symbols_to_char_set(UNAMBIGUOUS_AA_SYMBOLS)

# --- Full (Ambiguous) Alphabets ---
# These include all symbols: canonical, ambiguous, and gaps.

const ALL_DNA_SYMBOLS = BioSymbols.alphabet(BioSequences.DNA)
const ALL_RNA_SYMBOLS = BioSymbols.alphabet(BioSymbols.RNA)
const ALL_AA_SYMBOLS  = BioSymbols.alphabet(BioSymbols.AminoAcid)

const AMBIGUOUS_DNA_CHARSET = symbols_to_char_set(ALL_DNA_SYMBOLS)
const AMBIGUOUS_RNA_CHARSET = symbols_to_char_set(ALL_RNA_SYMBOLS)
const AMBIGUOUS_AA_CHARSET  = symbols_to_char_set(ALL_AA_SYMBOLS)

# can add support for conda too if needed
# const CONDA_RUNNER = joinpath(Conda.BINDIR, "mamba")
const CONDA_RUNNER = joinpath(Conda.BINDIR, "conda")
const FASTA_REGEX = r"\.(fa|fasta|fna|fas|fsa|ffn|faa|mpfa|frn)(\.gz)?$"
const FASTQ_REGEX = r"\.(fq|fastq)(\.gz)?$"
const XAM_REGEX = r"\.(sam|bam|cram|sam\.gz)$"
const VCF_REGEX = r"\.vcf(\.gz)?$"

# smaller, higher diversity databases do better with >=5 as the denominator - w/ <=4 they run out of memory
# denominator = 5 # produced OOM for NT on NERSC
# denominator = 8 # produced OOM for NT on Lawrencium
# denominator = 10 was only 56% efficient for NT on NERSC
const DEFAULT_MINIMAP_DENOMINATOR=10

# BLOSUM62 diagonal self-scores (others default to 0)
const BLOSUM62_DIAG = Dict(
    'A'=>4, 'R'=>5, 'N'=>6, 'D'=>6, 'C'=>9, 'Q'=>5, 'E'=>5, 'G'=>6,
    'H'=>8, 'I'=>4, 'L'=>4, 'K'=>5, 'M'=>5, 'F'=>6, 'P'=>7, 'S'=>4,
    'T'=>5, 'W'=>11,'Y'=>7, 'V'=>4
)

const BLAST_LAMBDA = 0.318    # Approximate BLAST (BLOSUM62) parameter
const BLAST_K = 0.134

# Reduced amino acid alphabet mappings
# These map each of the 20 standard amino acids to simplified groups based on physicochemical properties
# References:
# - Murphy et al. (2000) Protein Eng. 13(3):149-152
# - Peterson et al. (2009) BMC Bioinformatics 10:228
# - Zheng et al. (2019) Database (Oxford) baz131 (RAACBook)
#
# See Also:
# - [`reduce_amino_acid_alphabet()`](@ref): Main conversion function
# - [`list_reduced_alphabets()`](@ref): Get list of available schemes
# - [`get_reduced_alphabet_info()`](@ref): Get scheme metadata

"""
Binary hydrophobic/polar classification (simplest reduction).

Maps each of the 20 standard amino acids to either:
- 'H' (Hydrophobic): A, C, F, I, L, M, V, W
- 'P' (Polar): G, T, S, Y, P, N, D, E, Q, K, R, H

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:HP2`
- [`REDUCED_ALPHABET_HYDROPATHY3`](@ref): Three-class hydropathy
"""
const REDUCED_ALPHABET_HP2 = Dict(
    BioSymbols.AA_A => 'H', BioSymbols.AA_C => 'H', BioSymbols.AA_F => 'H', BioSymbols.AA_I => 'H',
    BioSymbols.AA_L => 'H', BioSymbols.AA_M => 'H', BioSymbols.AA_V => 'H', BioSymbols.AA_W => 'H',
    BioSymbols.AA_G => 'P', BioSymbols.AA_T => 'P', BioSymbols.AA_S => 'P', BioSymbols.AA_Y => 'P',
    BioSymbols.AA_P => 'P', BioSymbols.AA_N => 'P', BioSymbols.AA_D => 'P', BioSymbols.AA_E => 'P',
    BioSymbols.AA_Q => 'P', BioSymbols.AA_K => 'P', BioSymbols.AA_R => 'P', BioSymbols.AA_H => 'P'
)

"""
Three-class hydropathy (hydrophobic/neutral/hydrophilic) based on IMGT classification.

Maps amino acids to:
- 'H' (Hydrophobic): I, V, L, F, C, M, A, W
- 'N' (Neutral): G, T, S, Y, P, H
- 'P' (Polar/Hydrophilic): D, N, E, Q, K, R

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:HYDROPATHY3`
"""
const REDUCED_ALPHABET_HYDROPATHY3 = Dict(
    BioSymbols.AA_I => 'H', BioSymbols.AA_V => 'H', BioSymbols.AA_L => 'H', BioSymbols.AA_F => 'H',
    BioSymbols.AA_C => 'H', BioSymbols.AA_M => 'H', BioSymbols.AA_A => 'H', BioSymbols.AA_W => 'H',
    BioSymbols.AA_G => 'N', BioSymbols.AA_T => 'N', BioSymbols.AA_S => 'N', BioSymbols.AA_Y => 'N',
    BioSymbols.AA_P => 'N', BioSymbols.AA_H => 'N',
    BioSymbols.AA_D => 'P', BioSymbols.AA_N => 'P', BioSymbols.AA_E => 'P', BioSymbols.AA_Q => 'P',
    BioSymbols.AA_K => 'P', BioSymbols.AA_R => 'P'
)

"""
Four-class alphabet isolating Gly and Pro.

Glycine is the smallest and most flexible, Proline is rigid and breaks secondary structure.

Maps amino acids to:
- 'G' (Glycine): G
- 'P' (Proline): P
- 'H' (Hydrophobic): Y, F, L, I, V, M, C, W, H
- 'B' (Basic/polar): A, D, K, E, R, N, T, S, Q

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:GBMR4`
"""
const REDUCED_ALPHABET_GBMR4 = Dict(
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_Y => 'H', BioSymbols.AA_F => 'H', BioSymbols.AA_L => 'H', BioSymbols.AA_I => 'H',
    BioSymbols.AA_V => 'H', BioSymbols.AA_M => 'H', BioSymbols.AA_C => 'H', BioSymbols.AA_W => 'H',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_A => 'B', BioSymbols.AA_D => 'B', BioSymbols.AA_K => 'B', BioSymbols.AA_E => 'B',
    BioSymbols.AA_R => 'B', BioSymbols.AA_N => 'B', BioSymbols.AA_T => 'B', BioSymbols.AA_S => 'B',
    BioSymbols.AA_Q => 'B'
)

"""
Five classes based on chemical/structural properties.

Maps amino acids to:
- 'A' (Aliphatic): I, V, L
- 'R' (aRomatic): F, Y, W, H
- 'C' (Charged): K, R, D, E
- 'T' (Tiny): G, A, C, S
- 'D' (Diverse): T, M, Q, N, P

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:CHEMICAL5`
- [`REDUCED_ALPHABET_CHEMICAL6`](@ref): Six-class with separated charges
"""
const REDUCED_ALPHABET_CHEMICAL5 = Dict(
    BioSymbols.AA_I => 'A', BioSymbols.AA_V => 'A', BioSymbols.AA_L => 'A',
    BioSymbols.AA_F => 'R', BioSymbols.AA_Y => 'R', BioSymbols.AA_W => 'R', BioSymbols.AA_H => 'R',
    BioSymbols.AA_K => 'C', BioSymbols.AA_R => 'C', BioSymbols.AA_D => 'C', BioSymbols.AA_E => 'C',
    BioSymbols.AA_G => 'T', BioSymbols.AA_A => 'T', BioSymbols.AA_C => 'T', BioSymbols.AA_S => 'T',
    BioSymbols.AA_T => 'D', BioSymbols.AA_M => 'D', BioSymbols.AA_Q => 'D', BioSymbols.AA_N => 'D',
    BioSymbols.AA_P => 'D'
)

"""
Six classes separating positive and negative charges.

Maps amino acids to:
- 'A' (Aliphatic): I, V, L
- 'R' (aRomatic): F, Y, W, H
- '+' (Positive): K, R
- '-' (Negative): D, E
- 'T' (Tiny): G, A, C, S
- 'D' (Diverse): T, M, Q, N, P

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:CHEMICAL6`
- [`REDUCED_ALPHABET_CHEMICAL5`](@ref): Five-class with combined charges
"""
const REDUCED_ALPHABET_CHEMICAL6 = Dict(
    BioSymbols.AA_I => 'A', BioSymbols.AA_V => 'A', BioSymbols.AA_L => 'A',
    BioSymbols.AA_F => 'R', BioSymbols.AA_Y => 'R', BioSymbols.AA_W => 'R', BioSymbols.AA_H => 'R',
    BioSymbols.AA_K => '+', BioSymbols.AA_R => '+',
    BioSymbols.AA_D => '-', BioSymbols.AA_E => '-',
    BioSymbols.AA_G => 'T', BioSymbols.AA_A => 'T', BioSymbols.AA_C => 'T', BioSymbols.AA_S => 'T',
    BioSymbols.AA_T => 'D', BioSymbols.AA_M => 'D', BioSymbols.AA_Q => 'D', BioSymbols.AA_N => 'D',
    BioSymbols.AA_P => 'D'
)

"""
Structure-dependent model with 12 groups (Murphy et al. 2000).

Maintains clusters for acidic/basic, polar, aromatic, and aliphatic groups.

Maps amino acids to:
- 'A': A (Alanine alone)
- 'D': D, E (Acidic)
- 'K': K, R (Basic)
- 'N': N, Q (Amide)
- 'T': T, S (Hydroxyl)
- 'Y': Y, F (Aromatic)
- 'L': L, I, V, M (Aliphatic)
- 'C': C (Cysteine alone)
- 'W': W (Tryptophan alone)
- 'H': H (Histidine alone)
- 'G': G (Glycine alone)
- 'P': P (Proline alone)

# References
Murphy et al. (2000) Protein Eng. 13(3):149-152

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SDM12`
"""
const REDUCED_ALPHABET_SDM12 = Dict(
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_K => 'K', BioSymbols.AA_R => 'K',
    BioSymbols.AA_N => 'N', BioSymbols.AA_Q => 'N',
    BioSymbols.AA_T => 'T', BioSymbols.AA_S => 'T',
    BioSymbols.AA_Y => 'Y', BioSymbols.AA_F => 'Y',
    BioSymbols.AA_L => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_V => 'L', BioSymbols.AA_M => 'L',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_P => 'P'
)

"""
Dictionary mapping scheme names to their reduction dictionaries.

Contains all available reduced amino acid alphabet schemes. Each scheme maps
`BioSymbols.AminoAcid` symbols to reduced alphabet characters.

# Available Schemes
- `:HP2` - Binary hydrophobic/polar (2 classes)
- `:HYDROPATHY3` - Three-class hydropathy (3 classes)
- `:GBMR4` - Four-class with isolated Gly/Pro (4 classes)
- `:CHEMICAL5` - Five-class chemical properties (5 classes)
- `:CHEMICAL6` - Six-class with separated charges (6 classes)
- `:SDM12` - Structure-dependent model (12 classes)

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Main conversion function
- [`list_reduced_alphabets()`](@ref): Get list of available schemes
- [`REDUCED_ALPHABET_INFO`](@ref): Metadata about each scheme
"""
const REDUCED_ALPHABETS = Dict{Symbol, Dict{BioSymbols.AminoAcid, Char}}(
    :HP2 => REDUCED_ALPHABET_HP2,
    :HYDROPATHY3 => REDUCED_ALPHABET_HYDROPATHY3,
    :GBMR4 => REDUCED_ALPHABET_GBMR4,
    :CHEMICAL5 => REDUCED_ALPHABET_CHEMICAL5,
    :CHEMICAL6 => REDUCED_ALPHABET_CHEMICAL6,
    :SDM12 => REDUCED_ALPHABET_SDM12
)

"""
Metadata dictionary for all reduced alphabet schemes.

Contains detailed information about each reduction scheme including:
- `:name` - Full descriptive name
- `:classes` - Number of classes/groups
- `:description` - Purpose and usage description
- `:groups` - Mapping of reduced characters to amino acid groups

# See Also
- [`get_reduced_alphabet_info()`](@ref): Access metadata for a specific scheme
- [`REDUCED_ALPHABETS`](@ref): The actual reduction mappings
"""
const REDUCED_ALPHABET_INFO = Dict{Symbol, Dict{Symbol, Any}}(
    :HP2 => Dict(
        :name => "Binary Hydrophobic-Polar",
        :classes => 2,
        :description => "Simplest reduction grouping amino acids by hydrophobicity (H) vs polarity (P)",
        :groups => Dict('H' => "ACFILMVW", 'P' => "GTSYPNDEQKRH")
    ),
    :HYDROPATHY3 => Dict(
        :name => "Three-class Hydropathy",
        :classes => 3,
        :description => "IMGT hydropathy classification: Hydrophobic (H), Neutral (N), Polar (P)",
        :groups => Dict('H' => "IVLFCMAW", 'N' => "GTSYPH", 'P' => "DNEQKR")
    ),
    :GBMR4 => Dict(
        :name => "GBMR Four-class",
        :classes => 4,
        :description => "Isolates Gly (G) and Pro (P), separates Hydrophobic (H) and Basic/polar (B)",
        :groups => Dict('G' => "G", 'P' => "P", 'H' => "YFILMVCWH", 'B' => "ADKERNTSQ")
    ),
    :CHEMICAL5 => Dict(
        :name => "Five-class Chemical",
        :classes => 5,
        :description => "Aliphatic (A), aRomatic (R), Charged (C), Tiny (T), Diverse (D)",
        :groups => Dict('A' => "IVL", 'R' => "FYWH", 'C' => "KRDE", 'T' => "GACS", 'D' => "TMQNP")
    ),
    :CHEMICAL6 => Dict(
        :name => "Six-class Chemical",
        :classes => 6,
        :description => "Aliphatic (A), aRomatic (R), Positive (+), Negative (-), Tiny (T), Diverse (D)",
        :groups => Dict('A' => "IVL", 'R' => "FYWH", '+' => "KR", '-' => "DE", 'T' => "GACS", 'D' => "TMQNP")
    ),
    :SDM12 => Dict(
        :name => "SDM 12-class",
        :classes => 12,
        :description => "Structure-dependent model preserving acidic/basic, polar, aromatic, and aliphatic clusters",
        :groups => Dict(
            'A' => "A", 'D' => "DE", 'K' => "KR", 'N' => "NQ", 'T' => "TS",
            'Y' => "YF", 'L' => "LIVM", 'C' => "C", 'W' => "W", 'H' => "H",
            'G' => "G", 'P' => "P"
        )
    )
)

ProgressMeter.ijulia_behavior(:clear)