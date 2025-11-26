# # Tutorial 11: Reduced Amino Acid Alphabets
#
# This tutorial demonstrates how to convert amino acid sequences to reduced alphabets
# based on physicochemical properties. Reduced alphabets simplify protein sequences
# by grouping similar amino acids together, which can improve computational efficiency
# and reduce noise in machine learning applications.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - What reduced amino acid alphabets are and why they are useful
# - How to convert protein sequences to different reduction schemes
# - The differences between various reduction strategies
# - Applications of reduced alphabets in sequence analysis
#
# ## Background
#
# The standard genetic code uses 20 amino acids, each with distinct physicochemical
# properties. However, for many computational analyses, this diversity can introduce
# noise and computational complexity. Reduced amino acid alphabets group similar
# amino acids based on properties like:
#
# - **Hydrophobicity**: How water-repelling or water-attracting an amino acid is
# - **Charge**: Positive, negative, or neutral electrical charge
# - **Size**: Physical size of the amino acid side chain
# - **Aromaticity**: Presence of aromatic rings in the structure
# - **Structure**: Special properties like flexibility (Gly) or rigidity (Pro)
#
# Research has shown that reduced alphabets can:
# - Improve protein fold recognition
# - Reduce noise in machine learning models
# - Speed up sequence comparisons
# - Simplify pattern discovery in protein sequences

# ## Setup

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import BioSequences
import Statistics

# ## Part 1: Basic Conversion
#
# Let's start with a simple example using a fragment of human insulin A chain.

## Example protein sequence (fragment of human insulin A chain)
insulin_sequence = BioSequences.LongAA("GIVEQCCTSICSLYQLENYCN")
println("Original sequence: ", insulin_sequence)
println("Length: ", length(insulin_sequence))
println()

# ## Part 2: Exploring Available Schemes
#
# Mycelia provides six well-established reduction schemes from the literature.
# Let's list them all:

println("Available reduction schemes:")
for scheme in Mycelia.list_reduced_alphabets()
    info = Mycelia.get_reduced_alphabet_info(scheme)
    println("  $scheme ($(info[:classes]) classes): $(info[:name])")
end
println()

# ## Part 3: Binary Hydrophobic-Polar (HP2)
#
# The simplest reduction divides amino acids into just two classes:
# - **H** (Hydrophobic): Water-repelling amino acids (A, C, F, I, L, M, V, W)
# - **P** (Polar): Water-attracting amino acids (G, T, S, Y, P, N, D, E, Q, K, R, H)
#
# This is the most aggressive reduction and is useful for identifying
# hydrophobic core regions versus surface-exposed regions in proteins.

println("=" ^ 70)
println("HP2: Binary Hydrophobic-Polar Reduction")
println("=" ^ 70)
reduced_hp2 = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, :HP2)
info = Mycelia.get_reduced_alphabet_info(:HP2)
println("Description: $(info[:description])")
println("Groups:")
for (letter, aas) in sort(collect(info[:groups]))
    println("  $letter: $aas")
end
println("\nOriginal: $insulin_sequence")
println("Reduced:  $reduced_hp2")
println()

# ## Part 4: Three-Class Hydropathy (HYDROPATHY3)
#
# A more nuanced approach uses three classes based on IMGT hydropathy:
# - **H** (Hydrophobic): Strongly water-repelling
# - **N** (Neutral): Neither strongly hydrophobic nor hydrophilic
# - **P** (Polar): Strongly water-attracting

println("=" ^ 70)
println("HYDROPATHY3: Three-class Hydropathy")
println("=" ^ 70)
reduced_hydro3 = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, :HYDROPATHY3)
info = Mycelia.get_reduced_alphabet_info(:HYDROPATHY3)
println("Description: $(info[:description])")
println("Groups:")
for (letter, aas) in sort(collect(info[:groups]))
    println("  $letter: $aas")
end
println("\nOriginal: $insulin_sequence")
println("Reduced:  $reduced_hydro3")
println()

# ## Part 5: GBMR4 - Isolating Special Amino Acids
#
# The GBMR4 scheme recognizes that Glycine and Proline have unique structural roles:
# - **G** (Glycine): Smallest, most flexible
# - **P** (Proline): Rigid, breaks secondary structure
# - **H** (Hydrophobic): Standard hydrophobic amino acids
# - **B** (Basic/polar): Polar and charged amino acids

println("=" ^ 70)
println("GBMR4: Four-class (isolating Gly and Pro)")
println("=" ^ 70)
reduced_gbmr4 = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, :GBMR4)
info = Mycelia.get_reduced_alphabet_info(:GBMR4)
println("Description: $(info[:description])")
println("Groups:")
for (letter, aas) in sort(collect(info[:groups]))
    println("  $letter: $aas")
end
println("\nOriginal: $insulin_sequence")
println("Reduced:  $reduced_gbmr4")
println()

# ## Part 6: Chemical Properties (CHEMICAL6)
#
# The CHEMICAL6 scheme groups by detailed chemical properties:
# - **A** (Aliphatic): Non-aromatic hydrophobic chains
# - **R** (aRomatic): Contains aromatic rings
# - **+** (Positive): Positively charged (basic)
# - **-** (Negative): Negatively charged (acidic)
# - **T** (Tiny): Small amino acids
# - **D** (Diverse): Other properties
#
# This scheme is particularly useful for studying electrostatic interactions.

println("=" ^ 70)
println("CHEMICAL6: Six-class Chemical Properties")
println("=" ^ 70)
reduced_chem6 = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, :CHEMICAL6)
info = Mycelia.get_reduced_alphabet_info(:CHEMICAL6)
println("Description: $(info[:description])")
println("Groups:")
for (letter, aas) in sort(collect(info[:groups]))
    println("  $letter: $aas")
end
println("\nOriginal: $insulin_sequence")
println("Reduced:  $reduced_chem6")
println()

# ## Part 7: Structure-Dependent Model (SDM12)
#
# The SDM12 scheme (Murphy et al. 2000) maintains the most detail with 12 classes.
# It preserves important functional groupings while still reducing complexity.
# This is useful when you need more resolution than HP2 but less than the full 20.

println("=" ^ 70)
println("SDM12: Structure-dependent Model (12 classes)")
println("=" ^ 70)
reduced_sdm12 = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, :SDM12)
info = Mycelia.get_reduced_alphabet_info(:SDM12)
println("Description: $(info[:description])")
println("Groups:")
for (letter, aas) in sort(collect(info[:groups]))
    println("  $letter: $aas")
end
println("\nOriginal: $insulin_sequence")
println("Reduced:  $reduced_sdm12")
println()

# ## Part 8: Comparing All Schemes
#
# Let's see all reductions side-by-side to compare how they simplify the sequence:

println("=" ^ 70)
println("Comparison of All Schemes")
println("=" ^ 70)
println("Original:    $insulin_sequence")
for scheme in Mycelia.list_reduced_alphabets()
    reduced = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, scheme)
    println(rpad("$scheme:", 13), reduced)
end
println()

# ## Part 9: Application - K-mer Pattern Analysis
#
# One practical application of reduced alphabets is simplifying k-mer patterns.
# By reducing the alphabet, we reduce the number of possible k-mers, making
# patterns easier to detect and reducing data sparsity.

println("=" ^ 70)
println("Application: Reduced Alphabet K-mer Patterns")
println("=" ^ 70)

## Extract 3-mers from the original sequence
println("Original 3-mers from insulin sequence:")
original_3mers = [String(insulin_sequence[i:i+2]) for i in 1:(length(insulin_sequence)-2)]
println(join(original_3mers, ", "))
println()

## Extract 3-mers from HP2 reduced sequence
println("HP2 reduced 3-mers:")
reduced_seq = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, :HP2)
reduced_3mers = [reduced_seq[i:i+2] for i in 1:(length(reduced_seq)-2)]
println(join(reduced_3mers, ", "))
println()

## Compare pattern complexity
println("Unique 3-mer patterns:")
println("  Original: ", length(unique(original_3mers)), " unique patterns")
println("  HP2:      ", length(unique(reduced_3mers)), " unique patterns")
reduction_percent = round((1 - length(unique(reduced_3mers))/length(unique(original_3mers))) * 100, digits=1)
println("  Reduction: ", reduction_percent, "% fewer patterns")
println()

# ## Part 10: Practical Considerations
#
# When choosing a reduction scheme, consider:
#
# ### 1. Analysis Goals
# - **HP2**: Maximum simplification, focus on hydrophobicity
# - **HYDROPATHY3**: Balance between simplicity and detail
# - **GBMR4**: Important for structure-related analysis
# - **CHEMICAL5/6**: When charge matters (e.g., binding sites)
# - **SDM12**: When you need more detail but still want reduction
#
# ### 2. Information Loss
# - More aggressive reductions lose more information
# - But they also reduce noise and computational complexity
# - Choose based on signal-to-noise ratio in your data
#
# ### 3. Literature Compatibility
# - Use schemes from published studies for comparison
# - Murphy et al. (2000) SDM12 is widely cited
# - HP2 is the most established for folding studies

# ## Summary
#
# In this tutorial, we learned:
# - How to convert amino acid sequences to reduced alphabets
# - The differences between 6 well-established reduction schemes
# - How reduced alphabets simplify k-mer patterns
# - When to use different reduction strategies
#
# ## References
#
# - Murphy et al. (2000) Protein Eng. 13(3):149-152 - Original SDM12 paper
# - Peterson et al. (2009) BMC Bioinformatics 10:228 - Automated reduction methods
# - Zheng et al. (2019) Database (Oxford) baz131 - RAACBook comprehensive database
#
# ## Further Reading
#
# For more information on protein sequence analysis in Mycelia:
# - Tutorial 2: Quality Control (quality metrics for sequences)
# - Tutorial 3: K-mer Analysis (k-mer counting and analysis)
# - Tutorial 6: Gene Annotation (protein-coding sequence analysis)

println("=== Tutorial Complete ===")
