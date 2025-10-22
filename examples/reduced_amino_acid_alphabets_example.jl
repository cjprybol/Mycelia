## Example: Reduced Amino Acid Alphabet Conversion
##
## This example demonstrates how to convert amino acid sequences to reduced alphabets
## based on physicochemical properties. Reduced alphabets simplify protein sequences
## by grouping similar amino acids together, which can improve computational efficiency
## and reduce noise in machine learning applications.

import Mycelia
import BioSequences

## Example protein sequence (fragment of human insulin A chain)
insulin_sequence = BioSequences.LongAA("GIVEQCCTSICSLYQLENYCN")
println("Original sequence: ", insulin_sequence)
println("Length: ", length(insulin_sequence))
println()

## List all available reduction schemes
println("Available reduction schemes:")
for scheme in Mycelia.list_reduced_alphabets()
    info = Mycelia.get_reduced_alphabet_info(scheme)
    println("  $scheme ($(info[:classes]) classes): $(info[:name])")
end
println()

## Example 1: Binary hydrophobic/polar reduction (HP2)
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

## Example 2: Three-class hydropathy (HYDROPATHY3)
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

## Example 3: GBMR4 - isolating Gly and Pro
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

## Example 4: CHEMICAL6 - separating positive and negative charges
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

## Example 5: SDM12 - high-resolution structure-dependent model
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

## Example 6: Compare all schemes side-by-side
println("=" ^ 70)
println("Comparison of All Schemes")
println("=" ^ 70)
println("Original:    $insulin_sequence")
for scheme in Mycelia.list_reduced_alphabets()
    reduced = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, scheme)
    println(rpad("$scheme:", 13), reduced)
end
println()

## Example 7: Using reduced alphabets for k-mer analysis
println("=" ^ 70)
println("Application: Reduced Alphabet K-mer Patterns")
println("=" ^ 70)
println("Original 3-mers from insulin sequence:")
original_3mers = [String(insulin_sequence[i:i+2]) for i in 1:(length(insulin_sequence)-2)]
println(join(original_3mers, ", "))

println("\nHP2 reduced 3-mers:")
reduced_seq = Mycelia.reduce_amino_acid_alphabet(insulin_sequence, :HP2)
reduced_3mers = [reduced_seq[i:i+2] for i in 1:(length(reduced_seq)-2)]
println(join(reduced_3mers, ", "))

## Count unique patterns
println("\nUnique 3-mer patterns:")
println("  Original: ", length(unique(original_3mers)), " unique patterns")
println("  HP2:      ", length(unique(reduced_3mers)), " unique patterns")
println("  Reduction: ", round((1 - length(unique(reduced_3mers))/length(unique(original_3mers))) * 100, digits=1), "% fewer patterns")
