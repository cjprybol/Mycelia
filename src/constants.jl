# const METADATA = joinpath(dirname(dirname(pathof(Mycelia))), "docs", "metadata")

const DEFAULT_BLASTDB_PATH = "$(homedir())/workspace/blastdb"
const DEFAULT_METABULI_DB_PATH = "$(homedir())/workspace/metabuli"
const DEFAULT_METABULI_DB_NAME = "refseq_release"
const METABULI_DB_URLS = Dict(
    "gtdb+virus+human" => "https://metabuli.steineggerlab.workers.dev/gtdb+virus+human.tar.gz",
    "gtdb" => "https://metabuli.steineggerlab.workers.dev/gtdb.tar.gz",
    "refseq_prokaryote_virus" => "https://metabuli.steineggerlab.workers.dev/refseq_prokaryote_virus.tar.gz",
    "refseq_release" => "https://metabuli.steineggerlab.workers.dev/refseq_release.tar.gz",
    "refseq_virus" => "https://metabuli.steineggerlab.workers.dev/refseq_virus.tar.gz"
)
const DEFAULT_METAPHLAN_DB_PATH = "$(homedir())/workspace/metaphlan"

# fix new error
# ENV["MAMBA_ROOT_PREFIX"] = joinpath(DEPOT_PATH[1], "conda", "3", "x86_64")

# Mycelia.NERSC_MEM * .95
# const NERSC_MEM=512
# const NERSC_MEM=480
const NERSC_MEM=460
# const NERSC_CPU=240
const NERSC_CPU=240

# Default thread count for non-HPC environments or when no HPC allocation detected
# Used as fallback by get_default_threads() function
const DEFAULT_THREADS = 4

# Maximum threads supported by Flye assembler
const FLYE_MAX_THREADS = 128

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
Binary Murphy HP model (hydrophobic/polar) based on BLOSUM50 correlations.

This is the Murphy et al. (2000) version of the classic HP model, capturing ~75-80%
of contact mutual information with the highest statistical significance (Z-score 10.40).

Maps amino acids to:
- 'I': L, V, I, M, C, A, G, S, T, P, F, Y, W (Hydrophobic/Small)
- 'E': E, D, N, Q, K, R, H (Hydrophilic/Polar)

# References
Murphy L.R., Wallqvist A., Levy R.M. (2000) Simplified amino acid alphabets for
protein fold recognition and implications for folding. Protein Eng. 13(3):149-152

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:MURPHY2`
- [`REDUCED_ALPHABET_HP2`](@ref): Alternative 2-letter HP scheme
"""
const REDUCED_ALPHABET_MURPHY2 = Dict(
    BioSymbols.AA_L => 'I', BioSymbols.AA_V => 'I', BioSymbols.AA_I => 'I', BioSymbols.AA_M => 'I',
    BioSymbols.AA_C => 'I', BioSymbols.AA_A => 'I', BioSymbols.AA_G => 'I', BioSymbols.AA_S => 'I',
    BioSymbols.AA_T => 'I', BioSymbols.AA_P => 'I', BioSymbols.AA_F => 'I', BioSymbols.AA_Y => 'I',
    BioSymbols.AA_W => 'I',
    BioSymbols.AA_E => 'E', BioSymbols.AA_D => 'E', BioSymbols.AA_N => 'E', BioSymbols.AA_Q => 'E',
    BioSymbols.AA_K => 'E', BioSymbols.AA_R => 'E', BioSymbols.AA_H => 'E'
)

"""
Three-class Murphy alphabet based on BLOSUM50 correlations.

Captures ~82% of contact MI. Represents minimal folding complexity theorized
as the minimum required for folding heteropolymers.

Maps amino acids to:
- 'L': L, A, S, G, V, T, I, P, M, C
- 'E': E, K, R, D, N, Q, H (Charged/Polar)
- 'F': F, Y, W (Aromatic)

# References
Murphy L.R., Wallqvist A., Levy R.M. (2000) Simplified amino acid alphabets for
protein fold recognition and implications for folding. Protein Eng. 13(3):149-152

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:MURPHY3`
"""
const REDUCED_ALPHABET_MURPHY3 = Dict(
    BioSymbols.AA_L => 'L', BioSymbols.AA_A => 'L', BioSymbols.AA_S => 'L', BioSymbols.AA_G => 'L',
    BioSymbols.AA_V => 'L', BioSymbols.AA_T => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_P => 'L',
    BioSymbols.AA_M => 'L', BioSymbols.AA_C => 'L',
    BioSymbols.AA_E => 'E', BioSymbols.AA_K => 'E', BioSymbols.AA_R => 'E', BioSymbols.AA_D => 'E',
    BioSymbols.AA_N => 'E', BioSymbols.AA_Q => 'E', BioSymbols.AA_H => 'E',
    BioSymbols.AA_F => 'F', BioSymbols.AA_Y => 'F', BioSymbols.AA_W => 'F'
)

"""
Four-class Murphy alphabet (MU4) - minimal practical fold assessment.

Identified as the point where information loss becomes steep during further reduction.
High utility in distant homology detection with accuracy comparable to full alphabets.

Maps amino acids to:
- 'S': A, G, P, S, T (Small/Flexible)
- 'L': C, I, L, M, V (Aliphatic Hydrophobic)
- 'E': D, E, H, K, N, Q, R (Charged/Polar)
- 'F': F, Y, W (Aromatic)

# References
Murphy L.R., Wallqvist A., Levy R.M. (2000) Simplified amino acid alphabets for
protein fold recognition and implications for folding. Protein Eng. 13(3):149-152

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:MU4`
- [`REDUCED_ALPHABET_GBMR4`](@ref): Alternative 4-letter scheme
"""
const REDUCED_ALPHABET_MU4 = Dict(
    BioSymbols.AA_A => 'S', BioSymbols.AA_G => 'S', BioSymbols.AA_P => 'S', BioSymbols.AA_S => 'S',
    BioSymbols.AA_T => 'S',
    BioSymbols.AA_C => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_L => 'L', BioSymbols.AA_M => 'L',
    BioSymbols.AA_V => 'L',
    BioSymbols.AA_D => 'E', BioSymbols.AA_E => 'E', BioSymbols.AA_H => 'E', BioSymbols.AA_K => 'E',
    BioSymbols.AA_N => 'E', BioSymbols.AA_Q => 'E', BioSymbols.AA_R => 'E',
    BioSymbols.AA_F => 'F', BioSymbols.AA_Y => 'F', BioSymbols.AA_W => 'F'
)

"""
Five-class Murphy/BLOSUM alphabet (ML5) - optimal balance.

Maximally balances simplicity and retention (~90% MI retained). Experimentally shown
sufficient to preserve fold/function (e.g., Src SH3 domain).

Maps amino acids to:
- 'L': L, V, I, M, C (Aliphatic)
- 'A': A, S, G, T, P (Small)
- 'F': F, Y, W (Aromatic)
- 'E': E, D, N, Q (Acidic/Polar)
- 'K': K, R, H (Basic)

# References
Murphy L.R., Wallqvist A., Levy R.M. (2000) Simplified amino acid alphabets for
protein fold recognition and implications for folding. Protein Eng. 13(3):149-152

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:ML5`
- [`REDUCED_ALPHABET_WW5`](@ref), [`REDUCED_ALPHABET_MM5`](@ref): Alternative 5-letter schemes
"""
const REDUCED_ALPHABET_ML5 = Dict(
    BioSymbols.AA_L => 'L', BioSymbols.AA_V => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_M => 'L',
    BioSymbols.AA_C => 'L',
    BioSymbols.AA_A => 'A', BioSymbols.AA_S => 'A', BioSymbols.AA_G => 'A', BioSymbols.AA_T => 'A',
    BioSymbols.AA_P => 'A',
    BioSymbols.AA_F => 'F', BioSymbols.AA_Y => 'F', BioSymbols.AA_W => 'F',
    BioSymbols.AA_E => 'E', BioSymbols.AA_D => 'E', BioSymbols.AA_N => 'E', BioSymbols.AA_Q => 'E',
    BioSymbols.AA_K => 'K', BioSymbols.AA_R => 'K', BioSymbols.AA_H => 'K'
)

"""
Five-class Wang & Wang alphabet (WW5).

Alternative 5-letter scheme with different grouping strategy. Performs similarly
to ML5 in fold assessment tasks.

Maps amino acids to:
- 'I': C, F, I, L, M, V, W, Y (Hydrophobic)
- 'A': A, H, T
- 'D': D, E (Acidic)
- 'G': G, P
- 'K': K, N, Q, R, S (Polar/Basic)

# References
Wang J., Wang W. - Citation pending; see docs/src/references.md

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:WW5`
- [`REDUCED_ALPHABET_ML5`](@ref), [`REDUCED_ALPHABET_MM5`](@ref): Alternative 5-letter schemes
"""
const REDUCED_ALPHABET_WW5 = Dict(
    BioSymbols.AA_C => 'I', BioSymbols.AA_F => 'I', BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I',
    BioSymbols.AA_M => 'I', BioSymbols.AA_V => 'I', BioSymbols.AA_W => 'I', BioSymbols.AA_Y => 'I',
    BioSymbols.AA_A => 'A', BioSymbols.AA_H => 'A', BioSymbols.AA_T => 'A',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_G => 'G', BioSymbols.AA_P => 'G',
    BioSymbols.AA_K => 'K', BioSymbols.AA_N => 'K', BioSymbols.AA_Q => 'K', BioSymbols.AA_R => 'K',
    BioSymbols.AA_S => 'K'
)

"""
Five-class Melo & Marti-Renom alphabet (MM5).

Alternative 5-letter scheme. Shows classification performance only slightly lower
than the 20-letter alphabet for large models.

Maps amino acids to:
- 'A': A, G (Tiny)
- 'C': C (Cysteine)
- 'D': D, E, K, N, P, Q, R, S, T (Large mixed group)
- 'I': F, I, L, M, V, W, Y (Hydrophobic)
- 'H': H (Histidine)

# References
Melo F., Marti-Renom M.A. - Citation pending; see docs/src/references.md

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:MM5`
- [`REDUCED_ALPHABET_ML5`](@ref), [`REDUCED_ALPHABET_WW5`](@ref): Alternative 5-letter schemes
"""
const REDUCED_ALPHABET_MM5 = Dict(
    BioSymbols.AA_A => 'A', BioSymbols.AA_G => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D', BioSymbols.AA_K => 'D', BioSymbols.AA_N => 'D',
    BioSymbols.AA_P => 'D', BioSymbols.AA_Q => 'D', BioSymbols.AA_R => 'D', BioSymbols.AA_S => 'D',
    BioSymbols.AA_T => 'D',
    BioSymbols.AA_F => 'I', BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I', BioSymbols.AA_M => 'I',
    BioSymbols.AA_V => 'I', BioSymbols.AA_W => 'I', BioSymbols.AA_Y => 'I',
    BioSymbols.AA_H => 'H'
)

"""
Twelve-class Murphy alphabet - maximum information retention with minimal reduction.

Retains 98.1% of contact MI. Consistently identified as maximizing fold recognition
specificity. Groups: E/Q together, D/N together (differs from SDM12 which groups D/E, N/Q).

Maps amino acids to:
- 'L': L, V, I, M (Aliphatic)
- 'C': C (Cysteine)
- 'A': A (Alanine)
- 'G': G (Glycine)
- 'S': S, T (Hydroxyl)
- 'P': P (Proline)
- 'F': F, Y (Aromatic)
- 'W': W (Tryptophan)
- 'E': E, Q (Glu/Gln - amide-related)
- 'D': D, N (Asp/Asn - amide-related)
- 'K': K, R (Basic)
- 'H': H (Histidine)

# References
Murphy L.R., Wallqvist A., Levy R.M. (2000) Simplified amino acid alphabets for
protein fold recognition and implications for folding. Protein Eng. 13(3):149-152

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:MURPHY12`
- [`REDUCED_ALPHABET_SDM12`](@ref): Alternative 12-letter scheme (legacy)
- [`REDUCED_ALPHABET_SDM12_PRLIC`](@ref): Structure-derived 12-letter scheme
"""
const REDUCED_ALPHABET_MURPHY12 = Dict(
    BioSymbols.AA_L => 'L', BioSymbols.AA_V => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_M => 'L',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_S => 'S', BioSymbols.AA_T => 'S',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_F => 'F', BioSymbols.AA_Y => 'F',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_E => 'E', BioSymbols.AA_Q => 'E',
    BioSymbols.AA_D => 'D', BioSymbols.AA_N => 'D',
    BioSymbols.AA_K => 'K', BioSymbols.AA_R => 'K',
    BioSymbols.AA_H => 'H'
)

"""
Twelve-class structure-derived alphabet (SDM12) by Prlić et al.

Top overall performer in detecting structurally related proteins (best AUC). Derived from
structural alignments of proteins with low sequence identity. Groups K/E/R together (mixed
charge) and T/S/Q together, differing from Murphy12.

Maps amino acids to:
- 'A': A (Alanine)
- 'D': D (Aspartate - kept separate from K/E/R cluster)
- 'K': K, E, R (Mixed charge cluster: Lys, Glu, Arg)
- 'N': N (Asparagine)
- 'T': T, S, Q (Hydroxyl/amide: Thr, Ser, Gln)
- 'Y': Y, F (Aromatic)
- 'L': L, I, V, M (Aliphatic)
- 'C': C (Cysteine)
- 'W': W (Tryptophan)
- 'H': H (Histidine)
- 'G': G (Glycine)
- 'P': P (Proline)

# References
Prlić A., Domingues F.S., Sippl M.J. (2000) Structure-derived substitution matrices
for alignment of distantly related sequences. Protein Eng.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SDM12_PRLIC`
- [`REDUCED_ALPHABET_MURPHY12`](@ref): BLOSUM-derived 12-letter scheme
- [`REDUCED_ALPHABET_SDM12`](@ref): Legacy 12-letter scheme
"""
const REDUCED_ALPHABET_SDM12_PRLIC = Dict(
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_D => 'D',
    BioSymbols.AA_K => 'K', BioSymbols.AA_E => 'K', BioSymbols.AA_R => 'K',
    BioSymbols.AA_N => 'N',
    BioSymbols.AA_T => 'T', BioSymbols.AA_S => 'T', BioSymbols.AA_Q => 'T',
    BioSymbols.AA_Y => 'Y', BioSymbols.AA_F => 'Y',
    BioSymbols.AA_L => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_V => 'L', BioSymbols.AA_M => 'L',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_P => 'P'
)

"""
Seventeen-class structure-derived alphabet (HSDM17) by Prlić et al.

Best selectivity as measured by mean pooled precision (MPP). Only the strongest
associations are maintained: K/E together (acidic/basic) and L/I/V (aliphatic).
Most amino acids remain separate compared to SDM12.

Maps amino acids to:
- 'A': A (Alanine)
- 'D': D (Aspartate)
- 'K': K, E (Lys/Glu cluster)
- 'R': R (Arginine - separate from K/E)
- 'N': N (Asparagine)
- 'T': T (Threonine)
- 'S': S (Serine)
- 'Q': Q (Glutamine)
- 'Y': Y (Tyrosine)
- 'F': F (Phenylalanine)
- 'L': L, I, V (Aliphatic)
- 'M': M (Methionine - separate from L/I/V)
- 'C': C (Cysteine)
- 'W': W (Tryptophan)
- 'H': H (Histidine)
- 'G': G (Glycine)
- 'P': P (Proline)

# References
Prlić A., Domingues F.S., Sippl M.J. (2000) Structure-derived substitution matrices
for alignment of distantly related sequences. Protein Eng.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:HSDM17`
- [`REDUCED_ALPHABET_SDM12_PRLIC`](@ref): More aggressive 12-letter reduction
"""
const REDUCED_ALPHABET_HSDM17 = Dict(
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_D => 'D',
    BioSymbols.AA_K => 'K', BioSymbols.AA_E => 'K',
    BioSymbols.AA_R => 'R',
    BioSymbols.AA_N => 'N',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_S => 'S',
    BioSymbols.AA_Q => 'Q',
    BioSymbols.AA_Y => 'Y',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_L => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_V => 'L',
    BioSymbols.AA_M => 'M',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_P => 'P'
)

"""
    REDUCED_ALPHABET_MURPHY8

Murphy 8-class reduced amino acid alphabet.

# Description
Eight-class reduction from Murphy et al. (2000) providing good balance between
information retention and alphabet simplification.

# Groups (8 classes)
- 'L': L,V,I,M,C (Aliphatic)
- 'A': A,G (Small)
- 'S': S,T (Hydroxyl)
- 'P': P (Proline)
- 'F': F,Y,W (Aromatic)
- 'E': E,D,N,Q (Acidic/Polar)
- 'K': K,R (Basic)
- 'H': H (Histidine)

# Reference
Murphy LR, Wallqvist A, Levy RM. Simplified amino acid alphabets for protein
fold recognition and implications for folding. Protein Eng. 2000;13(3):149-152.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:MURPHY8`
- [`REDUCED_ALPHABET_MURPHY10`](@ref): Less aggressive 10-letter reduction
"""
const REDUCED_ALPHABET_MURPHY8 = Dict(
    BioSymbols.AA_L => 'L', BioSymbols.AA_V => 'L', BioSymbols.AA_I => 'L',
    BioSymbols.AA_M => 'L', BioSymbols.AA_C => 'L',
    BioSymbols.AA_A => 'A', BioSymbols.AA_G => 'A',
    BioSymbols.AA_S => 'S', BioSymbols.AA_T => 'S',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_F => 'F', BioSymbols.AA_Y => 'F', BioSymbols.AA_W => 'F',
    BioSymbols.AA_E => 'E', BioSymbols.AA_D => 'E', BioSymbols.AA_N => 'E', BioSymbols.AA_Q => 'E',
    BioSymbols.AA_K => 'K', BioSymbols.AA_R => 'K',
    BioSymbols.AA_H => 'H'
)

"""
    REDUCED_ALPHABET_MURPHY10

Murphy 10-class reduced amino acid alphabet.

# Description
Ten-class reduction from Murphy et al. (2000) with Cys, Ala, and Gly separated
for better structural representation.

# Groups (10 classes)
- 'L': L,V,I,M (Aliphatic)
- 'C': C (Cysteine)
- 'A': A (Alanine)
- 'G': G (Glycine)
- 'S': S,T (Hydroxyl)
- 'P': P (Proline)
- 'F': F,Y,W (Aromatic)
- 'E': E,D,N,Q (Acidic/Polar)
- 'K': K,R (Basic)
- 'H': H (Histidine)

# Reference
Murphy LR, Wallqvist A, Levy RM. Simplified amino acid alphabets for protein
fold recognition and implications for folding. Protein Eng. 2000;13(3):149-152.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:MURPHY10`
- [`REDUCED_ALPHABET_MURPHY15`](@ref): Less aggressive 15-letter reduction
"""
const REDUCED_ALPHABET_MURPHY10 = Dict(
    BioSymbols.AA_L => 'L', BioSymbols.AA_V => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_M => 'L',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_S => 'S', BioSymbols.AA_T => 'S',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_F => 'F', BioSymbols.AA_Y => 'F', BioSymbols.AA_W => 'F',
    BioSymbols.AA_E => 'E', BioSymbols.AA_D => 'E', BioSymbols.AA_N => 'E', BioSymbols.AA_Q => 'E',
    BioSymbols.AA_K => 'K', BioSymbols.AA_R => 'K',
    BioSymbols.AA_H => 'H'
)

"""
    REDUCED_ALPHABET_MURPHY15

Murphy 15-class reduced amino acid alphabet.

# Description
Fifteen-class reduction from Murphy et al. (2000) with high information retention.
Maintains key structural and chemical distinctions while reducing alphabet complexity.

# Groups (15 classes)
- 'L': L,V,I,M (Aliphatic)
- 'C': C | 'A': A | 'G': G | 'S': S | 'T': T | 'P': P (Individual small/special)
- 'F': F,Y (Aromatic)
- 'W': W (Tryptophan)
- 'E': E | 'D': D | 'N': N | 'Q': Q (Individual polar/charged)
- 'K': K,R (Basic)
- 'H': H (Histidine)

# Reference
Murphy LR, Wallqvist A, Levy RM. Simplified amino acid alphabets for protein
fold recognition and implications for folding. Protein Eng. 2000;13(3):149-152.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:MURPHY15`
- [`REDUCED_ALPHABET_MURPHY12`](@ref): More aggressive 12-letter reduction
"""
const REDUCED_ALPHABET_MURPHY15 = Dict(
    BioSymbols.AA_L => 'L', BioSymbols.AA_V => 'L', BioSymbols.AA_I => 'L', BioSymbols.AA_M => 'L',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_S => 'S',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_F => 'F', BioSymbols.AA_Y => 'F',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_E => 'E',
    BioSymbols.AA_D => 'D',
    BioSymbols.AA_N => 'N',
    BioSymbols.AA_Q => 'Q',
    BioSymbols.AA_K => 'K', BioSymbols.AA_R => 'K',
    BioSymbols.AA_H => 'H'
)

"""
    REDUCED_ALPHABET_SOLIS2

Solis 2-class reduced amino acid alphabet (minimal reduction).

# Description
Binary hydrophobic/hydrophilic classification from Solis (2015). Preserves fold
information from contact interactions.

# Groups (2 classes)
- 'C': C,F,I,L,M,V,W,Y (Hydrophobic)
- 'A': A,D,E,G,H,K,N,P,Q,R,S,T (Hydrophilic)

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.
doi: 10.1002/prot.24936

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS2`
- [`REDUCED_ALPHABET_SOLIS3`](@ref): Next level in Solis hierarchical set
"""
const REDUCED_ALPHABET_SOLIS2 = Dict(
    BioSymbols.AA_C => 'C', BioSymbols.AA_F => 'C', BioSymbols.AA_I => 'C',
    BioSymbols.AA_L => 'C', BioSymbols.AA_M => 'C', BioSymbols.AA_V => 'C',
    BioSymbols.AA_W => 'C', BioSymbols.AA_Y => 'C',
    BioSymbols.AA_A => 'A', BioSymbols.AA_D => 'A', BioSymbols.AA_E => 'A',
    BioSymbols.AA_G => 'A', BioSymbols.AA_H => 'A', BioSymbols.AA_K => 'A',
    BioSymbols.AA_N => 'A', BioSymbols.AA_P => 'A', BioSymbols.AA_Q => 'A',
    BioSymbols.AA_R => 'A', BioSymbols.AA_S => 'A', BioSymbols.AA_T => 'A'
)

"""
    REDUCED_ALPHABET_SOLIS3

Solis 3-class reduced amino acid alphabet.

# Description
Three-class reduction from Solis (2015) hierarchical set.

# Groups (3 classes)
- 'C': C,F,I,L,M,V,W,Y (Hydrophobic)
- 'D': D,E,G,K,N,Q,S (Polar)
- 'A': A,H,P,R,T

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS3`
"""
const REDUCED_ALPHABET_SOLIS3 = Dict(
    BioSymbols.AA_C => 'C', BioSymbols.AA_F => 'C', BioSymbols.AA_I => 'C',
    BioSymbols.AA_L => 'C', BioSymbols.AA_M => 'C', BioSymbols.AA_V => 'C',
    BioSymbols.AA_W => 'C', BioSymbols.AA_Y => 'C',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D', BioSymbols.AA_G => 'D',
    BioSymbols.AA_K => 'D', BioSymbols.AA_N => 'D', BioSymbols.AA_Q => 'D', BioSymbols.AA_S => 'D',
    BioSymbols.AA_A => 'A', BioSymbols.AA_H => 'A', BioSymbols.AA_P => 'A',
    BioSymbols.AA_R => 'A', BioSymbols.AA_T => 'A'
)

"""
    REDUCED_ALPHABET_SOLIS4

Solis 4-class reduced amino acid alphabet.

# Description
Four-class reduction from Solis (2015) with separated aromatic residues.

# Groups (4 classes)
- 'F': F,W,Y (Aromatic)
- 'C': C,I,L,M,V (Aliphatic)
- 'D': D,E,G,K,N,Q,S (Polar)
- 'A': A,H,P,R,T

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS4`
"""
const REDUCED_ALPHABET_SOLIS4 = Dict(
    BioSymbols.AA_F => 'F', BioSymbols.AA_W => 'F', BioSymbols.AA_Y => 'F',
    BioSymbols.AA_C => 'C', BioSymbols.AA_I => 'C', BioSymbols.AA_L => 'C',
    BioSymbols.AA_M => 'C', BioSymbols.AA_V => 'C',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D', BioSymbols.AA_G => 'D',
    BioSymbols.AA_K => 'D', BioSymbols.AA_N => 'D', BioSymbols.AA_Q => 'D', BioSymbols.AA_S => 'D',
    BioSymbols.AA_A => 'A', BioSymbols.AA_H => 'A', BioSymbols.AA_P => 'A',
    BioSymbols.AA_R => 'A', BioSymbols.AA_T => 'A'
)

"""
    REDUCED_ALPHABET_SOLIS5

Solis 5-class reduced amino acid alphabet.

# Description
Five-class reduction from Solis (2015) hierarchical set.

# Groups (5 classes)
- 'F': F,W,Y (Aromatic)
- 'C': C,I,L,M,V (Aliphatic)
- 'D': D,E,G,K,N,S (Charged/Polar)
- 'A': A,P,Q,T
- 'H': H,R (Positive)

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS5`
"""
const REDUCED_ALPHABET_SOLIS5 = Dict(
    BioSymbols.AA_F => 'F', BioSymbols.AA_W => 'F', BioSymbols.AA_Y => 'F',
    BioSymbols.AA_C => 'C', BioSymbols.AA_I => 'C', BioSymbols.AA_L => 'C',
    BioSymbols.AA_M => 'C', BioSymbols.AA_V => 'C',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D', BioSymbols.AA_G => 'D',
    BioSymbols.AA_K => 'D', BioSymbols.AA_N => 'D', BioSymbols.AA_S => 'D',
    BioSymbols.AA_A => 'A', BioSymbols.AA_P => 'A', BioSymbols.AA_Q => 'A', BioSymbols.AA_T => 'A',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H'
)

"""
    REDUCED_ALPHABET_SOLIS6

Solis 6-class reduced amino acid alphabet.

# Description
Six-class reduction from Solis (2015) with separated acidic residues.

# Groups (6 classes)
- 'F': F,W,Y (Aromatic)
- 'C': C,I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'G': G,K,N,Q,S (Mixed polar)
- 'A': A,P,T
- 'H': H,R (Basic)

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS6`
"""
const REDUCED_ALPHABET_SOLIS6 = Dict(
    BioSymbols.AA_F => 'F', BioSymbols.AA_W => 'F', BioSymbols.AA_Y => 'F',
    BioSymbols.AA_C => 'C', BioSymbols.AA_I => 'C', BioSymbols.AA_L => 'C',
    BioSymbols.AA_M => 'C', BioSymbols.AA_V => 'C',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_G => 'G', BioSymbols.AA_K => 'G', BioSymbols.AA_N => 'G',
    BioSymbols.AA_Q => 'G', BioSymbols.AA_S => 'G',
    BioSymbols.AA_A => 'A', BioSymbols.AA_P => 'A', BioSymbols.AA_T => 'A',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H'
)

"""
    REDUCED_ALPHABET_SOLIS7

Solis 7-class reduced amino acid alphabet.

# Description
Seven-class reduction from Solis (2015) with Lysine separated.

# Groups (7 classes)
- 'F': F,W,Y (Aromatic)
- 'C': C,I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'K': K (Lysine)
- 'G': G,N,P,Q,S
- 'A': A,T
- 'H': H,R (Basic)

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS7`
"""
const REDUCED_ALPHABET_SOLIS7 = Dict(
    BioSymbols.AA_F => 'F', BioSymbols.AA_W => 'F', BioSymbols.AA_Y => 'F',
    BioSymbols.AA_C => 'C', BioSymbols.AA_I => 'C', BioSymbols.AA_L => 'C',
    BioSymbols.AA_M => 'C', BioSymbols.AA_V => 'C',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_K => 'K',
    BioSymbols.AA_G => 'G', BioSymbols.AA_N => 'G', BioSymbols.AA_P => 'G',
    BioSymbols.AA_Q => 'G', BioSymbols.AA_S => 'G',
    BioSymbols.AA_A => 'A', BioSymbols.AA_T => 'A',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H'
)

"""
    REDUCED_ALPHABET_SOLIS8

Solis 8-class reduced amino acid alphabet.

# Description
Eight-class reduction from Solis (2015) with Cysteine separated.

# Groups (8 classes)
- 'F': F,W,Y (Aromatic)
- 'I': I,L,M,V (Aliphatic)
- 'C': C (Cysteine)
- 'D': D,E (Acidic)
- 'K': K (Lysine)
- 'G': G,N,P,Q,S
- 'A': A,T
- 'H': H,R (Basic)

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS8`
"""
const REDUCED_ALPHABET_SOLIS8 = Dict(
    BioSymbols.AA_F => 'F', BioSymbols.AA_W => 'F', BioSymbols.AA_Y => 'F',
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I', BioSymbols.AA_M => 'I', BioSymbols.AA_V => 'I',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_K => 'K',
    BioSymbols.AA_G => 'G', BioSymbols.AA_N => 'G', BioSymbols.AA_P => 'G',
    BioSymbols.AA_Q => 'G', BioSymbols.AA_S => 'G',
    BioSymbols.AA_A => 'A', BioSymbols.AA_T => 'A',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H'
)

"""
    REDUCED_ALPHABET_SOLIS9

Solis 9-class reduced amino acid alphabet.

# Description
Nine-class reduction from Solis (2015) hierarchical set.

# Groups (9 classes)
- 'F': F,W,Y (Aromatic)
- 'I': I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'G': G,N,Q,S
- 'P': P,T
- 'A': A (Alanine)
- 'H': H,R (Basic)
- 'C': C (Cysteine)
- 'K': K (Lysine)

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS9`
"""
const REDUCED_ALPHABET_SOLIS9 = Dict(
    BioSymbols.AA_F => 'F', BioSymbols.AA_W => 'F', BioSymbols.AA_Y => 'F',
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I', BioSymbols.AA_M => 'I', BioSymbols.AA_V => 'I',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_G => 'G', BioSymbols.AA_N => 'G', BioSymbols.AA_Q => 'G', BioSymbols.AA_S => 'G',
    BioSymbols.AA_P => 'P', BioSymbols.AA_T => 'P',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS10

Solis 10-class reduced amino acid alphabet.

# Description
Ten-class reduction from Solis (2015) with Phenylalanine separated from aromatics.

# Groups (10 classes)
- 'W': W,Y (Large aromatic)
- 'F': F (Phenylalanine)
- 'I': I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'G': G,N,Q,S
- 'P': P,T
- 'A': A (Alanine)
- 'H': H,R (Basic)
- 'C': C (Cysteine)
- 'K': K (Lysine)

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS10`
"""
const REDUCED_ALPHABET_SOLIS10 = Dict(
    BioSymbols.AA_W => 'W', BioSymbols.AA_Y => 'W',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I', BioSymbols.AA_M => 'I', BioSymbols.AA_V => 'I',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_G => 'G', BioSymbols.AA_N => 'G', BioSymbols.AA_Q => 'G', BioSymbols.AA_S => 'G',
    BioSymbols.AA_P => 'P', BioSymbols.AA_T => 'P',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS11

Solis 11-class reduced amino acid alphabet.

# Description
Eleven-class reduction from Solis (2015) with Glycine and Threonine separated.

# Groups (11 classes)
- 'W': W,Y (Large aromatic)
- 'I': I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'G': G (Glycine)
- 'N': N,P,Q,S
- 'T': T (Threonine)
- 'H': H,R (Basic)
- 'F': F (Phenylalanine)
- 'A': A (Alanine)
- 'C': C (Cysteine)
- 'K': K (Lysine)

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS11`
"""
const REDUCED_ALPHABET_SOLIS11 = Dict(
    BioSymbols.AA_W => 'W', BioSymbols.AA_Y => 'W',
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I', BioSymbols.AA_M => 'I', BioSymbols.AA_V => 'I',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_N => 'N', BioSymbols.AA_P => 'N', BioSymbols.AA_Q => 'N', BioSymbols.AA_S => 'N',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS12

Solis 12-class reduced amino acid alphabet.

# Description
Twelve-class reduction from Solis (2015) with further aliphatic separation.

# Groups (12 classes)
- 'W': W,Y (Large aromatic)
- 'I': I,L (Aliphatic branched)
- 'M': M,V (Aliphatic sulfur/branched)
- 'D': D,E (Acidic)
- 'N': N,P,Q,S
- 'H': H,R (Basic)
- 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS12`
"""
const REDUCED_ALPHABET_SOLIS12 = Dict(
    BioSymbols.AA_W => 'W', BioSymbols.AA_Y => 'W',
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I',
    BioSymbols.AA_M => 'M', BioSymbols.AA_V => 'M',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_N => 'N', BioSymbols.AA_P => 'N', BioSymbols.AA_Q => 'N', BioSymbols.AA_S => 'N',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS13

Solis 13-class reduced amino acid alphabet.

# Description
Thirteen-class reduction from Solis (2015) with Proline separated.

# Groups (13 classes)
- 'W': W,Y (Large aromatic)
- 'I': I,L (Aliphatic branched)
- 'M': M,V (Aliphatic sulfur/branched)
- 'D': D,E (Acidic)
- 'P': P (Proline)
- 'N': N,Q,S
- 'H': H,R (Basic)
- 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS13`
"""
const REDUCED_ALPHABET_SOLIS13 = Dict(
    BioSymbols.AA_W => 'W', BioSymbols.AA_Y => 'W',
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I',
    BioSymbols.AA_M => 'M', BioSymbols.AA_V => 'M',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_N => 'N', BioSymbols.AA_Q => 'N', BioSymbols.AA_S => 'N',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS14

Solis 14-class reduced amino acid alphabet.

# Description
Fourteen-class reduction from Solis (2015) with W and Y separated.

# Groups (14 classes)
- 'W': W | 'Y': Y (Separate aromatics)
- 'I': I,L (Aliphatic branched)
- 'M': M,V (Aliphatic sulfur/branched)
- 'D': D,E (Acidic)
- 'P': P (Proline)
- 'N': N,Q,S
- 'H': H,R (Basic)
- 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS14`
"""
const REDUCED_ALPHABET_SOLIS14 = Dict(
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_Y => 'Y',
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I',
    BioSymbols.AA_M => 'M', BioSymbols.AA_V => 'M',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_N => 'N', BioSymbols.AA_Q => 'N', BioSymbols.AA_S => 'N',
    BioSymbols.AA_H => 'H', BioSymbols.AA_R => 'H',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS15

Solis 15-class reduced amino acid alphabet.

# Description
Fifteen-class reduction from Solis (2015) with H and R separated.

# Groups (15 classes)
- 'I': I,L (Aliphatic branched)
- 'M': M,V (Aliphatic sulfur/branched)
- 'D': D,E (Acidic)
- 'N': N,Q,S
- 'H': H | 'R': R (Separate basic)
- 'W': W | 'Y': Y | 'P': P | 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS15`
"""
const REDUCED_ALPHABET_SOLIS15 = Dict(
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I',
    BioSymbols.AA_M => 'M', BioSymbols.AA_V => 'M',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_N => 'N', BioSymbols.AA_Q => 'N', BioSymbols.AA_S => 'N',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_R => 'R',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_Y => 'Y',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS16

Solis 16-class reduced amino acid alphabet.

# Description
Sixteen-class reduction from Solis (2015) with M and V separated.

# Groups (16 classes)
- 'I': I,L (Aliphatic branched)
- 'M': M | 'V': V (Separate sulfur and branched)
- 'D': D,E (Acidic)
- 'N': N,Q,S
- 'H': H | 'R': R | 'W': W | 'Y': Y | 'P': P | 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS16`
"""
const REDUCED_ALPHABET_SOLIS16 = Dict(
    BioSymbols.AA_I => 'I', BioSymbols.AA_L => 'I',
    BioSymbols.AA_M => 'M',
    BioSymbols.AA_V => 'V',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_N => 'N', BioSymbols.AA_Q => 'N', BioSymbols.AA_S => 'N',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_R => 'R',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_Y => 'Y',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS17

Solis 17-class reduced amino acid alphabet.

# Description
Seventeen-class reduction from Solis (2015) with all aliphatic separated.

# Groups (17 classes)
- 'I': I | 'L': L | 'M': M | 'V': V (All aliphatic separate)
- 'D': D,E (Acidic)
- 'N': N,Q,S
- 'H': H | 'R': R | 'W': W | 'Y': Y | 'P': P | 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS17`
"""
const REDUCED_ALPHABET_SOLIS17 = Dict(
    BioSymbols.AA_I => 'I',
    BioSymbols.AA_L => 'L',
    BioSymbols.AA_M => 'M',
    BioSymbols.AA_V => 'V',
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_N => 'N', BioSymbols.AA_Q => 'N', BioSymbols.AA_S => 'N',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_R => 'R',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_Y => 'Y',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS18

Solis 18-class reduced amino acid alphabet.

# Description
Eighteen-class reduction from Solis (2015) with D,E separated and Q,S paired.

# Groups (18 classes)
- 'D': D,E (Acidic - only remaining pair) → Changed to D only
- 'N': N (Asparagine)
- 'Q': Q,S (Amide/hydroxyl)
- All others separate: I, L, M, V, H, R, W, Y, P, G, T, F, A, C, K, E

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS18`
"""
const REDUCED_ALPHABET_SOLIS18 = Dict(
    BioSymbols.AA_D => 'D', BioSymbols.AA_E => 'D',
    BioSymbols.AA_N => 'N',
    BioSymbols.AA_Q => 'Q', BioSymbols.AA_S => 'Q',
    BioSymbols.AA_I => 'I',
    BioSymbols.AA_L => 'L',
    BioSymbols.AA_M => 'M',
    BioSymbols.AA_V => 'V',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_R => 'R',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_Y => 'Y',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
    REDUCED_ALPHABET_SOLIS19

Solis 19-class reduced amino acid alphabet (nearly complete).

# Description
Nineteen-class reduction from Solis (2015) - only Q,S remain paired.

# Groups (19 classes)
- 'D': D | 'E': E (Separate acidic)
- 'Q': Q,S (Amide/hydroxyl - final pair)
- 'N': N | 'I': I | 'L': L | 'M': M | 'V': V | 'H': H | 'R': R | 'W': W | 'Y': Y | 'P': P | 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K

# Reference
Solis AD. Amino acid alphabet reduction preserves fold information contained in
contact interactions in proteins. Proteins. 2015;83(12):2198-2216.

# See Also
- [`reduce_amino_acid_alphabet()`](@ref): Use this mapping with `:SOLIS19`
"""
const REDUCED_ALPHABET_SOLIS19 = Dict(
    BioSymbols.AA_D => 'D',
    BioSymbols.AA_E => 'E',
    BioSymbols.AA_Q => 'Q', BioSymbols.AA_S => 'Q',
    BioSymbols.AA_N => 'N',
    BioSymbols.AA_I => 'I',
    BioSymbols.AA_L => 'L',
    BioSymbols.AA_M => 'M',
    BioSymbols.AA_V => 'V',
    BioSymbols.AA_H => 'H',
    BioSymbols.AA_R => 'R',
    BioSymbols.AA_W => 'W',
    BioSymbols.AA_Y => 'Y',
    BioSymbols.AA_P => 'P',
    BioSymbols.AA_G => 'G',
    BioSymbols.AA_T => 'T',
    BioSymbols.AA_F => 'F',
    BioSymbols.AA_A => 'A',
    BioSymbols.AA_C => 'C',
    BioSymbols.AA_K => 'K'
)

"""
Dictionary mapping scheme names to their reduction dictionaries.

Contains all available reduced amino acid alphabet schemes. Each scheme maps
`BioSymbols.AminoAcid` symbols to reduced alphabet characters.

# Available Schemes (36 total)
## Legacy Schemes (6)
- `:HP2` - Binary hydrophobic/polar (2 classes)
- `:HYDROPATHY3` - Three-class hydropathy (3 classes)
- `:GBMR4` - Four-class with isolated Gly/Pro (4 classes)
- `:CHEMICAL5` - Five-class chemical properties (5 classes)
- `:CHEMICAL6` - Six-class with separated charges (6 classes)
- `:SDM12` - Structure-dependent model (12 classes)

## Murphy et al. (2000) (8 schemes)
- `:MURPHY2` - Murphy 2-class HP model
- `:MURPHY3` - Murphy 3-class
- `:MU4` - Murphy 4-class
- `:ML5` - Murphy 5-class (optimal)
- `:MURPHY8` - Murphy 8-class
- `:MURPHY10` - Murphy 10-class
- `:MURPHY12` - Murphy 12-class
- `:MURPHY15` - Murphy 15-class

## Prlić et al. (2000) (2 schemes)
- `:SDM12_PRLIC` - Structure-derived 12-class
- `:HSDM17` - Structure-derived 17-class

## Solis (2015) Hierarchical Set (18 schemes)
- `:SOLIS2` through `:SOLIS19` - Complete hierarchical set from 2 to 19 classes

## Other Validated Schemes (2)
- `:WW5` - Wang & Wang 5-class
- `:MM5` - Melo & Marti-Renom 5-class

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
    :SDM12 => REDUCED_ALPHABET_SDM12,
    :MURPHY2 => REDUCED_ALPHABET_MURPHY2,
    :MURPHY3 => REDUCED_ALPHABET_MURPHY3,
    :MU4 => REDUCED_ALPHABET_MU4,
    :ML5 => REDUCED_ALPHABET_ML5,
    :MURPHY8 => REDUCED_ALPHABET_MURPHY8,
    :MURPHY10 => REDUCED_ALPHABET_MURPHY10,
    :MURPHY12 => REDUCED_ALPHABET_MURPHY12,
    :MURPHY15 => REDUCED_ALPHABET_MURPHY15,
    :WW5 => REDUCED_ALPHABET_WW5,
    :MM5 => REDUCED_ALPHABET_MM5,
    :SDM12_PRLIC => REDUCED_ALPHABET_SDM12_PRLIC,
    :HSDM17 => REDUCED_ALPHABET_HSDM17,
    :SOLIS2 => REDUCED_ALPHABET_SOLIS2,
    :SOLIS3 => REDUCED_ALPHABET_SOLIS3,
    :SOLIS4 => REDUCED_ALPHABET_SOLIS4,
    :SOLIS5 => REDUCED_ALPHABET_SOLIS5,
    :SOLIS6 => REDUCED_ALPHABET_SOLIS6,
    :SOLIS7 => REDUCED_ALPHABET_SOLIS7,
    :SOLIS8 => REDUCED_ALPHABET_SOLIS8,
    :SOLIS9 => REDUCED_ALPHABET_SOLIS9,
    :SOLIS10 => REDUCED_ALPHABET_SOLIS10,
    :SOLIS11 => REDUCED_ALPHABET_SOLIS11,
    :SOLIS12 => REDUCED_ALPHABET_SOLIS12,
    :SOLIS13 => REDUCED_ALPHABET_SOLIS13,
    :SOLIS14 => REDUCED_ALPHABET_SOLIS14,
    :SOLIS15 => REDUCED_ALPHABET_SOLIS15,
    :SOLIS16 => REDUCED_ALPHABET_SOLIS16,
    :SOLIS17 => REDUCED_ALPHABET_SOLIS17,
    :SOLIS18 => REDUCED_ALPHABET_SOLIS18,
    :SOLIS19 => REDUCED_ALPHABET_SOLIS19
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
    ),
    :MURPHY2 => Dict(
        :name => "Murphy 2-class",
        :classes => 2,
        :description => "Murphy HP model (BLOSUM50-derived) for fold recognition",
        :groups => Dict('I' => "LVIMC AGSTP FYW", 'E' => "EDNQKRH")
    ),
    :MURPHY3 => Dict(
        :name => "Murphy 3-class",
        :classes => 3,
        :description => "Murphy 3-class BLOSUM50-derived reduction",
        :groups => Dict('L' => "LASGVT IPMC", 'E' => "EKRDNQH", 'F' => "FYW")
    ),
    :MU4 => Dict(
        :name => "Murphy 4-class",
        :classes => 4,
        :description => "Murphy 4-class minimal fold assessment",
        :groups => Dict('S' => "AGPST", 'L' => "CILMV", 'E' => "DEHKNQR", 'F' => "FYW")
    ),
    :ML5 => Dict(
        :name => "Murphy 5-class",
        :classes => 5,
        :description => "Murphy/BLOSUM 5-class optimal balance (~90% MI retained)",
        :groups => Dict('L' => "LVIMC", 'A' => "ASGTP", 'F' => "FYW", 'E' => "EDNQ", 'K' => "KRH")
    ),
    :MURPHY8 => Dict(
        :name => "Murphy 8-class",
        :classes => 8,
        :description => "Murphy 8-class balanced reduction",
        :groups => Dict('L' => "LVIMC", 'A' => "AG", 'S' => "ST", 'P' => "P", 'F' => "FYW", 'E' => "EDNQ", 'K' => "KR", 'H' => "H")
    ),
    :MURPHY10 => Dict(
        :name => "Murphy 10-class",
        :classes => 10,
        :description => "Murphy 10-class with separated Cys, Ala, Gly",
        :groups => Dict('L' => "LVIM", 'C' => "C", 'A' => "A", 'G' => "G", 'S' => "ST", 'P' => "P", 'F' => "FYW", 'E' => "EDNQ", 'K' => "KR", 'H' => "H")
    ),
    :MURPHY12 => Dict(
        :name => "Murphy 12-class",
        :classes => 12,
        :description => "Murphy 12-class high information retention (98.1% MI)",
        :groups => Dict('L' => "LVIM", 'C' => "C", 'A' => "A", 'G' => "G", 'P' => "P", 'W' => "W", 'H' => "H", 'S' => "ST", 'F' => "FY", 'E' => "EQ", 'D' => "DN", 'K' => "KR")
    ),
    :MURPHY15 => Dict(
        :name => "Murphy 15-class",
        :classes => 15,
        :description => "Murphy 15-class high information retention",
        :groups => Dict('L' => "LVIM", 'C' => "C", 'A' => "A", 'G' => "G", 'S' => "S", 'T' => "T", 'P' => "P", 'F' => "FY", 'W' => "W", 'E' => "E", 'D' => "D", 'N' => "N", 'Q' => "Q", 'K' => "KR", 'H' => "H")
    ),
    :WW5 => Dict(
        :name => "Wang & Wang 5-class",
        :classes => 5,
        :description => "Wang & Wang 5-class computational alphabet simplification",
        :groups => Dict('I' => "CFILMVWY", 'A' => "AHT", 'D' => "DE", 'G' => "GP", 'K' => "KNQRS")
    ),
    :MM5 => Dict(
        :name => "Melo & Marti-Renom 5-class",
        :classes => 5,
        :description => "Melo & Marti-Renom 5-class for alignment accuracy",
        :groups => Dict('A' => "AG", 'C' => "C", 'D' => "DEKNPQRST", 'I' => "FILMVWY", 'H' => "H")
    ),
    :SDM12_PRLIC => Dict(
        :name => "Prlić SDM 12-class",
        :classes => 12,
        :description => "Prlić structure-derived 12-class (best AUC for fold recognition)",
        :groups => Dict('A' => "A", 'D' => "D", 'N' => "N", 'C' => "C", 'W' => "W", 'H' => "H", 'G' => "G", 'P' => "P", 'K' => "KER", 'T' => "TSQ", 'Y' => "YF", 'L' => "LIVM")
    ),
    :HSDM17 => Dict(
        :name => "Prlić HSDM 17-class",
        :classes => 17,
        :description => "Prlić structure-derived 17-class (best mean pooled precision)",
        :groups => Dict('A' => "A", 'D' => "D", 'R' => "R", 'N' => "N", 'T' => "T", 'S' => "S", 'Q' => "Q", 'Y' => "Y", 'F' => "F", 'M' => "M", 'C' => "C", 'W' => "W", 'H' => "H", 'G' => "G", 'P' => "P", 'K' => "KE", 'L' => "LIV")
    ),
    :SOLIS2 => Dict(
        :name => "Solis 2-class",
        :classes => 2,
        :description => "Solis 2-class minimal reduction (hydrophobic/hydrophilic)",
        :groups => Dict('C' => "CFILMVWY", 'A' => "ADEGHKNPQRST")
    ),
    :SOLIS3 => Dict(
        :name => "Solis 3-class",
        :classes => 3,
        :description => "Solis 3-class hierarchical reduction",
        :groups => Dict('C' => "CFILMVWY", 'D' => "DEGKNQS", 'A' => "AHPRT")
    ),
    :SOLIS4 => Dict(
        :name => "Solis 4-class",
        :classes => 4,
        :description => "Solis 4-class with separated aromatics",
        :groups => Dict('F' => "FWY", 'C' => "CILMV", 'D' => "DEGKNQS", 'A' => "AHPRT")
    ),
    :SOLIS5 => Dict(
        :name => "Solis 5-class",
        :classes => 5,
        :description => "Solis 5-class hierarchical reduction",
        :groups => Dict('F' => "FWY", 'C' => "CILMV", 'D' => "DEGKNS", 'A' => "APQT", 'H' => "HR")
    ),
    :SOLIS6 => Dict(
        :name => "Solis 6-class",
        :classes => 6,
        :description => "Solis 6-class with separated acidic residues",
        :groups => Dict('F' => "FWY", 'C' => "CILMV", 'D' => "DE", 'G' => "GKNQS", 'A' => "APT", 'H' => "HR")
    ),
    :SOLIS7 => Dict(
        :name => "Solis 7-class",
        :classes => 7,
        :description => "Solis 7-class with Lysine separated",
        :groups => Dict('F' => "FWY", 'C' => "CILMV", 'D' => "DE", 'K' => "K", 'G' => "GNPQS", 'A' => "AT", 'H' => "HR")
    ),
    :SOLIS8 => Dict(
        :name => "Solis 8-class",
        :classes => 8,
        :description => "Solis 8-class with Cysteine separated",
        :groups => Dict('F' => "FWY", 'I' => "ILMV", 'C' => "C", 'D' => "DE", 'K' => "K", 'G' => "GNPQS", 'A' => "AT", 'H' => "HR")
    ),
    :SOLIS9 => Dict(
        :name => "Solis 9-class",
        :classes => 9,
        :description => "Solis 9-class hierarchical reduction",
        :groups => Dict('F' => "FWY", 'I' => "ILMV", 'D' => "DE", 'G' => "GNQS", 'P' => "PT", 'A' => "A", 'H' => "HR", 'C' => "C", 'K' => "K")
    ),
    :SOLIS10 => Dict(
        :name => "Solis 10-class",
        :classes => 10,
        :description => "Solis 10-class with Phe separated",
        :groups => Dict('W' => "WY", 'F' => "F", 'I' => "ILMV", 'D' => "DE", 'G' => "GNQS", 'P' => "PT", 'A' => "A", 'H' => "HR", 'C' => "C", 'K' => "K")
    ),
    :SOLIS11 => Dict(
        :name => "Solis 11-class",
        :classes => 11,
        :description => "Solis 11-class with Gly and Thr separated",
        :groups => Dict('W' => "WY", 'I' => "ILMV", 'D' => "DE", 'G' => "G", 'N' => "NPQS", 'T' => "T", 'H' => "HR", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    ),
    :SOLIS12 => Dict(
        :name => "Solis 12-class",
        :classes => 12,
        :description => "Solis 12-class with aliphatic separation",
        :groups => Dict('W' => "WY", 'I' => "IL", 'M' => "MV", 'D' => "DE", 'N' => "NPQS", 'H' => "HR", 'G' => "G", 'T' => "T", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    ),
    :SOLIS13 => Dict(
        :name => "Solis 13-class",
        :classes => 13,
        :description => "Solis 13-class with Proline separated",
        :groups => Dict('W' => "WY", 'I' => "IL", 'M' => "MV", 'D' => "DE", 'P' => "P", 'N' => "NQS", 'H' => "HR", 'G' => "G", 'T' => "T", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    ),
    :SOLIS14 => Dict(
        :name => "Solis 14-class",
        :classes => 14,
        :description => "Solis 14-class with W and Y separated",
        :groups => Dict('W' => "W", 'Y' => "Y", 'I' => "IL", 'M' => "MV", 'D' => "DE", 'P' => "P", 'N' => "NQS", 'H' => "HR", 'G' => "G", 'T' => "T", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    ),
    :SOLIS15 => Dict(
        :name => "Solis 15-class",
        :classes => 15,
        :description => "Solis 15-class with H and R separated",
        :groups => Dict('I' => "IL", 'M' => "MV", 'D' => "DE", 'N' => "NQS", 'H' => "H", 'R' => "R", 'W' => "W", 'Y' => "Y", 'P' => "P", 'G' => "G", 'T' => "T", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    ),
    :SOLIS16 => Dict(
        :name => "Solis 16-class",
        :classes => 16,
        :description => "Solis 16-class with M and V separated",
        :groups => Dict('I' => "IL", 'M' => "M", 'V' => "V", 'D' => "DE", 'N' => "NQS", 'H' => "H", 'R' => "R", 'W' => "W", 'Y' => "Y", 'P' => "P", 'G' => "G", 'T' => "T", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    ),
    :SOLIS17 => Dict(
        :name => "Solis 17-class",
        :classes => 17,
        :description => "Solis 17-class with all aliphatic separated",
        :groups => Dict('I' => "I", 'L' => "L", 'M' => "M", 'V' => "V", 'D' => "DE", 'N' => "NQS", 'H' => "H", 'R' => "R", 'W' => "W", 'Y' => "Y", 'P' => "P", 'G' => "G", 'T' => "T", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    ),
    :SOLIS18 => Dict(
        :name => "Solis 18-class",
        :classes => 18,
        :description => "Solis 18-class with Q,S paired",
        :groups => Dict('D' => "DE", 'N' => "N", 'Q' => "QS", 'I' => "I", 'L' => "L", 'M' => "M", 'V' => "V", 'H' => "H", 'R' => "R", 'W' => "W", 'Y' => "Y", 'P' => "P", 'G' => "G", 'T' => "T", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    ),
    :SOLIS19 => Dict(
        :name => "Solis 19-class",
        :classes => 19,
        :description => "Solis 19-class nearly complete (only Q,S paired)",
        :groups => Dict('D' => "D", 'E' => "E", 'Q' => "QS", 'N' => "N", 'I' => "I", 'L' => "L", 'M' => "M", 'V' => "V", 'H' => "H", 'R' => "R", 'W' => "W", 'Y' => "Y", 'P' => "P", 'G' => "G", 'T' => "T", 'F' => "F", 'A' => "A", 'C' => "C", 'K' => "K")
    )
)

ProgressMeter.ijulia_behavior(:clear)
