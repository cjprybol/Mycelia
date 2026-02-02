# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/6_annotation/codon_optimization.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/6_annotation/codon_optimization.jl", "test/6_annotation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import BioSequences
import StatsBase
import DataStructures
import Random
import Kmers

Test.@testset "Codon Optimization Tests" begin
    Test.@testset "reverse_translate function" begin
        # Test simple amino acid sequence
        simple_protein = BioSequences.LongAA("MKT")
        dna_result = Mycelia.reverse_translate(simple_protein)

        Test.@test dna_result isa BioSequences.LongDNA{2}
        Test.@test length(dna_result) == 9  # 3 amino acids * 3 nucleotides
        Test.@test BioSequences.translate(dna_result) == simple_protein

        # Test protein with stop codon
        protein_with_stop = BioSequences.LongAA("MKT*")
        dna_with_stop = Mycelia.reverse_translate(protein_with_stop)
        Test.@test length(dna_with_stop) == 12  # 4 codons * 3 nucleotides
        Test.@test BioSequences.translate(dna_with_stop) == protein_with_stop

        # Test single amino acid
        single_aa = BioSequences.LongAA("M")
        single_dna = Mycelia.reverse_translate(single_aa)
        Test.@test length(single_dna) == 3
        Test.@test BioSequences.translate(single_dna) == single_aa

        # Test empty sequence
        empty_protein = BioSequences.LongAA("")
        empty_dna = Mycelia.reverse_translate(empty_protein)
        Test.@test length(empty_dna) == 0
        Test.@test BioSequences.translate(empty_dna) == empty_protein

        # Test various amino acids
        diverse_protein = BioSequences.LongAA("MKFLILVAGF")
        diverse_dna = Mycelia.reverse_translate(diverse_protein)
        Test.@test BioSequences.translate(diverse_dna) == diverse_protein
    end

    Test.@testset "amino_acids_to_codons function" begin
        aa_to_codon_map = Mycelia.amino_acids_to_codons()

        Test.@test aa_to_codon_map isa Dict
        Test.@test haskey(aa_to_codon_map, BioSequences.AA_M)  # Methionine
        Test.@test haskey(aa_to_codon_map, BioSequences.AA_Term)  # Stop codon

        supported_aas = filter(
            aa -> !(aa in (BioSequences.AA_U, BioSequences.AA_O)),
            vcat(Mycelia.AA_ALPHABET..., [BioSequences.AA_Term])
        )

        # Test that all supported amino acids are represented
        for aa in supported_aas
            Test.@test haskey(aa_to_codon_map, aa)
            Test.@test aa_to_codon_map[aa] isa Kmers.DNACodon
        end

        # Test that codons translate to correct amino acids
        for aa in supported_aas
            codon = aa_to_codon_map[aa]
            translated_aa = first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
            Test.@test translated_aa == aa
        end
    end

    Test.@testset "codons_to_amino_acids function" begin
        codon_to_aa_map = Mycelia.codons_to_amino_acids()

        Test.@test codon_to_aa_map isa Dict
        Test.@test length(codon_to_aa_map) == 64  # All possible codons

        # Test specific known mappings
        expected_mappings = [
            ("ATG", BioSequences.AA_M),  # Start codon
            ("TAA", BioSequences.AA_Term),  # Stop codon
            ("TTT", BioSequences.AA_F),  # Phenylalanine
            ("AAA", BioSequences.AA_K)   # Lysine
        ]

        for (codon_str, expected_aa) in expected_mappings
            codon = Kmers.DNACodon(codon_str)
            Test.@test haskey(codon_to_aa_map, codon)
            Test.@test first(codon_to_aa_map[codon]) == expected_aa
        end
    end

    Test.@testset "normalize_codon_frequencies function" begin
        # Create test codon frequencies
        test_frequencies = Dict(
            BioSequences.AA_M => Dict(Kmers.DNACodon("ATG") => 10),
            BioSequences.AA_K => Dict(
                Kmers.DNACodon("AAA") => 30,
                Kmers.DNACodon("AAG") => 20
            ),
            BioSequences.AA_Term => Dict(
                Kmers.DNACodon("TAA") => 5,
                Kmers.DNACodon("TAG") => 3,
                Kmers.DNACodon("TGA") => 2
            )
        )

        normalized = Mycelia.normalize_codon_frequencies(test_frequencies)

        Test.@test normalized isa Dict

        # Test that frequencies sum to 1.0 for each amino acid
        for (aa, aa_frequencies) in normalized
            freq_sum = sum(values(aa_frequencies))
            Test.@test abs(freq_sum - 1.0) < eps()
        end

        # Test specific values
        Test.@test normalized[BioSequences.AA_M][Kmers.DNACodon("ATG")] ≈ 1.0
        Test.@test normalized[BioSequences.AA_K][Kmers.DNACodon("AAA")] ≈ 0.6
        Test.@test normalized[BioSequences.AA_K][Kmers.DNACodon("AAG")] ≈ 0.4

        # Test with empty frequencies
        empty_freq = Dict(BioSequences.AA_A => Dict{Kmers.DNACodon, Int}())
        normalized_empty = Mycelia.normalize_codon_frequencies(empty_freq)
        Test.@test isempty(normalized_empty[BioSequences.AA_A])
    end

    Test.@testset "normalize_kmer_counts function" begin
        # Test with simple k-mer counts
        test_kmers = Dict(
            "ATG" => 10,
            "GCT" => 20,
            "TAA" => 5
        )

        normalized = Mycelia.normalize_kmer_counts(test_kmers)

        Test.@test normalized isa DataStructures.OrderedDict
        Test.@test abs(sum(values(normalized)) - 1.0) < eps()

        # Test specific values
        total = 35
        Test.@test normalized["ATG"] ≈ 10/total
        Test.@test normalized["GCT"] ≈ 20/total
        Test.@test normalized["TAA"] ≈ 5/total

        # Test with single k-mer
        single_kmer = Dict("ATG" => 100)
        normalized_single = Mycelia.normalize_kmer_counts(single_kmer)
        Test.@test normalized_single["ATG"] ≈ 1.0

        # Test with zero counts (edge case)
        zero_counts = Dict("ATG" => 0, "GCT" => 0)
        normalized_zero = Mycelia.normalize_kmer_counts(zero_counts)
        Test.@test all(iszero, values(normalized_zero))
    end

    Test.@testset "codon_optimize function" begin
        # Create simple normalized codon frequencies for testing
        simple_freq = Dict(
            BioSequences.AA_M => Dict(Kmers.DNACodon("ATG") => 1.0),
            BioSequences.AA_K => Dict(
                Kmers.DNACodon("AAA") => 0.7,
                Kmers.DNACodon("AAG") => 0.3
            ),
            BioSequences.AA_T => Dict(
                Kmers.DNACodon("ACT") => 0.4,
                Kmers.DNACodon("ACC") => 0.3,
                Kmers.DNACodon("ACA") => 0.2,
                Kmers.DNACodon("ACG") => 0.1
            )
        )

        test_protein = BioSequences.LongAA("MKT")
        optimized_dna = Mycelia.codon_optimize(
            normalized_codon_frequencies = simple_freq,
            protein_sequence = test_protein,
            n_iterations = 10
        )

        Test.@test optimized_dna isa BioSequences.LongDNA{2}
        Test.@test length(optimized_dna) == 9  # 3 amino acids * 3 nucleotides
        Test.@test BioSequences.translate(optimized_dna) == test_protein

        # Test with single amino acid
        single_protein = BioSequences.LongAA("M")
        single_optimized = Mycelia.codon_optimize(
            normalized_codon_frequencies = simple_freq,
            protein_sequence = single_protein,
            n_iterations = 5
        )
        Test.@test BioSequences.translate(single_optimized) == single_protein

        # Test that optimization is deterministic given the same random seed
        # (This test might be flaky due to randomness, but useful for debugging)
        Random.seed!(12345)
        opt1 = Mycelia.codon_optimize(
            normalized_codon_frequencies = simple_freq,
            protein_sequence = test_protein,
            n_iterations = 1
        )
        Random.seed!(12345)
        opt2 = Mycelia.codon_optimize(
            normalized_codon_frequencies = simple_freq,
            protein_sequence = test_protein,
            n_iterations = 1
        )
        Test.@test opt1 == opt2
    end

    Test.@testset "Codon optimization integration" begin
        # Test the complete workflow with simulated data
        # Create a more comprehensive frequency table
        comprehensive_freq = Dict()

        # Initialize frequencies for all amino acids
        for aa in vcat(Mycelia.AA_ALPHABET..., [BioSequences.AA_Term])
            comprehensive_freq[aa] = Dict{Kmers.DNACodon, Float64}()
        end

        # Add some realistic frequencies
        comprehensive_freq[BioSequences.AA_M] = Dict(Kmers.DNACodon("ATG") => 1.0)
        comprehensive_freq[BioSequences.AA_K] = Dict(
            Kmers.DNACodon("AAA") => 0.6,
            Kmers.DNACodon("AAG") => 0.4
        )
        comprehensive_freq[BioSequences.AA_F] = Dict(
            Kmers.DNACodon("TTT") => 0.5,
            Kmers.DNACodon("TTC") => 0.5
        )
        comprehensive_freq[BioSequences.AA_L] = Dict(
            Kmers.DNACodon("TTA") => 0.2,
            Kmers.DNACodon("TTG") => 0.2,
            Kmers.DNACodon("CTT") => 0.15,
            Kmers.DNACodon("CTC") => 0.15,
            Kmers.DNACodon("CTA") => 0.15,
            Kmers.DNACodon("CTG") => 0.15
        )

        # Fill in missing amino acids with uniform distributions
        for (aa, freq_dict) in comprehensive_freq
            if isempty(freq_dict)
                # Get all codons for this amino acid
                codons_for_aa = [codon
                                 for codon in
                                     Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
                                 if first(BioSequences.translate(BioSequences.LongDNA{2}(codon))) ==
                                    aa]
                if !isempty(codons_for_aa)
                    uniform_freq = 1.0 / length(codons_for_aa)
                    for codon in codons_for_aa
                        comprehensive_freq[aa][codon] = uniform_freq
                    end
                end
            end
        end

        test_protein = BioSequences.LongAA("MKFL")
        optimized = Mycelia.codon_optimize(
            normalized_codon_frequencies = comprehensive_freq,
            protein_sequence = test_protein,
            n_iterations = 20
        )

        Test.@test BioSequences.translate(optimized) == test_protein
        Test.@test length(optimized) == length(test_protein) * 3
    end

    Test.@testset "Error handling and edge cases" begin
        # Test with protein containing unusual amino acids
        # Note: Some amino acids might not be in the standard alphabet

        # Test empty protein sequence
        empty_protein = BioSequences.LongAA("")
        empty_result = Mycelia.reverse_translate(empty_protein)
        Test.@test length(empty_result) == 0

        # Test frequency normalization with edge cases
        edge_freq = Dict(
            BioSequences.AA_A => Dict{Kmers.DNACodon, Int}(),  # Empty
            BioSequences.AA_R => Dict(Kmers.DNACodon("CGT") => 1)  # Single codon
        )
        normalized_edge = Mycelia.normalize_codon_frequencies(edge_freq)
        Test.@test isempty(normalized_edge[BioSequences.AA_A])
        Test.@test normalized_edge[BioSequences.AA_R][Kmers.DNACodon("CGT")] ≈ 1.0
    end

    Test.@testset "Performance and consistency tests" begin
        # Test with longer sequences
        long_protein = BioSequences.LongAA(repeat("MKFLILVAGFPVK", 10))
        long_result = Mycelia.reverse_translate(long_protein)
        Test.@test BioSequences.translate(long_result) == long_protein
        Test.@test length(long_result) == length(long_protein) * 3

        # Test multiple reverse translations give valid results
        for i in 1:5
            result = Mycelia.reverse_translate(long_protein)
            Test.@test BioSequences.translate(result) == long_protein
        end
    end
end
