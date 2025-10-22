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

# HP2: Binary hydrophobic/polar classification (simplest reduction)
const REDUCED_ALPHABET_HP2 = Dict(
    BioSymbols.AA_A => 'H', BioSymbols.AA_C => 'H', BioSymbols.AA_F => 'H', BioSymbols.AA_I => 'H',
    BioSymbols.AA_L => 'H', BioSymbols.AA_M => 'H', BioSymbols.AA_V => 'H', BioSymbols.AA_W => 'H',
    BioSymbols.AA_G => 'P', BioSymbols.AA_T => 'P', BioSymbols.AA_S => 'P', BioSymbols.AA_Y => 'P',
    BioSymbols.AA_P => 'P', BioSymbols.AA_N => 'P', BioSymbols.AA_D => 'P', BioSymbols.AA_E => 'P',
    BioSymbols.AA_Q => 'P', BioSymbols.AA_K => 'P', BioSymbols.AA_R => 'P', BioSymbols.AA_H => 'P'
)

# HYDROPATHY3: Three-class hydropathy (hydrophobic/neutral/hydrophilic)
# Based on IMGT hydropathy classification
const REDUCED_ALPHABET_HYDROPATHY3 = Dict(
    BioSymbols.AA_I => 'H', BioSymbols.AA_V => 'H', BioSymbols.AA_L => 'H', BioSymbols.AA_F => 'H',
    BioSymbols.AA_C => 'H', BioSymbols.AA_M => 'H', BioSymbols.AA_A => 'H', BioSymbols.AA_W => 'H',
    BioSymbols.AA_G => 'N', BioSymbols.AA_T => 'N', BioSymbols.AA_S => 'N', BioSymbols.AA_Y => 'N',
    BioSymbols.AA_P => 'N', BioSymbols.AA_H => 'N',
    BioSymbols.AA_D => 'P', BioSymbols.AA_N => 'P', BioSymbols.AA_E => 'P', BioSymbols.AA_Q => 'P',
    BioSymbols.AA_K => 'P', BioSymbols.AA_R => 'P'
)

# GBMR4: Four-class alphabet isolating Gly and Pro
# Glycine is the smallest and most flexible, Proline is rigid and breaks secondary structure
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

# CHEMICAL5: Five classes based on chemical/structural properties
# A=Aliphatic, R=aRomatic, C=Charged, T=Tiny, D=Diverse
const REDUCED_ALPHABET_CHEMICAL5 = Dict(
    BioSymbols.AA_I => 'A', BioSymbols.AA_V => 'A', BioSymbols.AA_L => 'A',
    BioSymbols.AA_F => 'R', BioSymbols.AA_Y => 'R', BioSymbols.AA_W => 'R', BioSymbols.AA_H => 'R',
    BioSymbols.AA_K => 'C', BioSymbols.AA_R => 'C', BioSymbols.AA_D => 'C', BioSymbols.AA_E => 'C',
    BioSymbols.AA_G => 'T', BioSymbols.AA_A => 'T', BioSymbols.AA_C => 'T', BioSymbols.AA_S => 'T',
    BioSymbols.AA_T => 'D', BioSymbols.AA_M => 'D', BioSymbols.AA_Q => 'D', BioSymbols.AA_N => 'D',
    BioSymbols.AA_P => 'D'
)

# CHEMICAL6: Six classes separating positive and negative charges
# A=Aliphatic, R=aRomatic, +=positive, -=negative, T=Tiny, D=Diverse
const REDUCED_ALPHABET_CHEMICAL6 = Dict(
    BioSymbols.AA_I => 'A', BioSymbols.AA_V => 'A', BioSymbols.AA_L => 'A',
    BioSymbols.AA_F => 'R', BioSymbols.AA_Y => 'R', BioSymbols.AA_W => 'R', BioSymbols.AA_H => 'R',
    BioSymbols.AA_K => '+', BioSymbols.AA_R => '+',
    BioSymbols.AA_D => '-', BioSymbols.AA_E => '-',
    BioSymbols.AA_G => 'T', BioSymbols.AA_A => 'T', BioSymbols.AA_C => 'T', BioSymbols.AA_S => 'T',
    BioSymbols.AA_T => 'D', BioSymbols.AA_M => 'D', BioSymbols.AA_Q => 'D', BioSymbols.AA_N => 'D',
    BioSymbols.AA_P => 'D'
)

# SDM12: Structure-dependent model with 12 groups (Murphy et al. 2000)
# Maintains clusters for acidic/basic, polar, aromatic, and aliphatic groups
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

# Map of all available reduced alphabets
const REDUCED_ALPHABETS = Dict{Symbol, Dict{BioSymbols.AminoAcid, Char}}(
    :HP2 => REDUCED_ALPHABET_HP2,
    :HYDROPATHY3 => REDUCED_ALPHABET_HYDROPATHY3,
    :GBMR4 => REDUCED_ALPHABET_GBMR4,
    :CHEMICAL5 => REDUCED_ALPHABET_CHEMICAL5,
    :CHEMICAL6 => REDUCED_ALPHABET_CHEMICAL6,
    :SDM12 => REDUCED_ALPHABET_SDM12
)

# Metadata about each reduced alphabet scheme
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