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
    x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x)),
    BioSymbols.alphabet(BioSymbols.AminoAcid))

# const AA_ALPHABET = filter(
#     x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x) || BioSymbols.isterm(x)),
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

ProgressMeter.ijulia_behavior(:clear)