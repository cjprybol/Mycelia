# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/sequence-complexity.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/sequence-complexity.jl", "test/2_preprocessing_qc", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import BioSequences
import Mycelia
import Statistics

Test.@testset "Sequence Comparison Functions" begin
    Test.@testset "shannon_entropy" begin
        Test.@testset "DNA sequences" begin
            # Test basic DNA sequence
            dna_seq = BioSequences.LongDNA{4}("ATCGATCGATCG")
            Test.@test Mycelia.shannon_entropy(dna_seq, k = 1) ≈ 2.0  # 4 bases equally distributed
            Test.@test Mycelia.shannon_entropy(dna_seq, k = 2) > 0.0

            # Test uniform sequence (minimum entropy)
            uniform_dna = BioSequences.LongDNA{4}("AAAAAAAAAA")
            Test.@test Mycelia.shannon_entropy(uniform_dna, k = 1) ≈ 0.0

            # Test sequence shorter than k
            short_dna = BioSequences.LongDNA{4}("AT")
            Test.@test Mycelia.shannon_entropy(short_dna, k = 3) == 0.0

            # Test different bases
            Test.@test Mycelia.shannon_entropy(BioSequences.LongDNA{4}("ATCG"), k = 1, base = 10) !=
                       Mycelia.shannon_entropy(BioSequences.LongDNA{4}("ATCG"), k = 1, base = 2)
        end

        Test.@testset "RNA sequences" begin
            rna_seq = BioSequences.LongRNA{4}("AUCGAUCGAUCG")
            Test.@test Mycelia.shannon_entropy(rna_seq, k = 1) ≈ 2.0
            Test.@test Mycelia.shannon_entropy(rna_seq, k = 2) > 0.0

            uniform_rna = BioSequences.LongRNA{4}("UUUUUUUUUU")
            Test.@test Mycelia.shannon_entropy(uniform_rna, k = 1) ≈ 0.0
        end

        Test.@testset "Amino acid sequences" begin
            aa_seq = BioSequences.LongAA("ARNDCQEGHILKMFPSTWYV")  # All 20 standard amino acids
            Test.@test Mycelia.shannon_entropy(aa_seq, k = 1) > 4.0  # log2(20) ≈ 4.32

            uniform_aa = BioSequences.LongAA("AAAAAAAAAA")
            Test.@test Mycelia.shannon_entropy(uniform_aa, k = 1) ≈ 0.0
        end

        Test.@testset "Generic strings" begin
            str = "abcdefgh"
            Test.@test Mycelia.shannon_entropy(str, k = 1) > 2.0
            Test.@test Mycelia.shannon_entropy(str, k = 2) > 0.0

            # Test uniform string
            uniform_str = "aaaaaaa"
            Test.@test Mycelia.shannon_entropy(uniform_str, k = 1) ≈ 0.0

            # Test empty and short strings
            Test.@test Mycelia.shannon_entropy("", k = 1) == 0.0
            Test.@test Mycelia.shannon_entropy("a", k = 2) == 0.0
        end

        Test.@testset "Edge cases" begin
            # Test with different alphabet parameter
            Test.@test Mycelia.shannon_entropy("abc", k = 1, alphabet = 3) > 0.0

            # Test different log bases
            seq = "abcd"
            Test.@test Mycelia.shannon_entropy(seq, k = 1, base = 2) !=
                       Mycelia.shannon_entropy(seq, k = 1, base = 10)
        end
    end

    Test.@testset "renyi_entropy" begin
        Test.@testset "DNA sequences" begin
            dna_seq = BioSequences.LongDNA{4}("ATCGATCGATCG")
            Test.@test Mycelia.renyi_entropy(dna_seq, k = 1, α = 2) > 0.0
            Test.@test Mycelia.renyi_entropy(dna_seq, k = 2, α = 0.5) > 0.0

            # Test that α=1 throws error
            Test.@test_throws AssertionError Mycelia.renyi_entropy(dna_seq, k = 1, α = 1)

            # Test sequence shorter than k
            Test.@test Mycelia.renyi_entropy(BioSequences.LongDNA{4}("AT"), k = 3, α = 2) ==
                       0.0
        end

        Test.@testset "RNA sequences" begin
            rna_seq = BioSequences.LongRNA{4}("AUCGAUCGAUCG")
            Test.@test Mycelia.renyi_entropy(rna_seq, k = 1, α = 2) > 0.0
            Test.@test Mycelia.renyi_entropy(rna_seq, k = 2, α = 0.5) > 0.0
        end

        Test.@testset "Amino acid sequences" begin
            aa_seq = BioSequences.LongAA("ARNDCQEGHILKMFPSTWYV")
            Test.@test Mycelia.renyi_entropy(aa_seq, k = 1, α = 2) > 0.0
            Test.@test Mycelia.renyi_entropy(aa_seq, k = 2, α = 3) > 0.0
        end

        Test.@testset "Generic strings" begin
            str = "abcdefgh"
            Test.@test Mycelia.renyi_entropy(str, k = 1, α = 2) > 0.0
            Test.@test Mycelia.renyi_entropy(str, k = 2, α = 0.5) > 0.0

            # Test uniform string
            uniform_str = "aaaaaaa"
            Test.@test Mycelia.renyi_entropy(uniform_str, k = 1, α = 2) ≈ 0.0
        end

        Test.@testset "Parameter effects" begin
            seq = "aabbbbcccc"  # Non-uniform distribution: a=2, b=4, c=4
            # Different α values should give different results
            Test.@test Mycelia.renyi_entropy(seq, k = 1, α = 0.5) !=
                       Mycelia.renyi_entropy(seq, k = 1, α = 2)
            Test.@test Mycelia.renyi_entropy(seq, k = 1, α = 2) !=
                       Mycelia.renyi_entropy(seq, k = 1, α = 10)

            # Different bases should give different results
            Test.@test Mycelia.renyi_entropy(seq, k = 1, α = 2, base = 2) !=
                       Mycelia.renyi_entropy(seq, k = 1, α = 2, base = 10)
        end
    end

    Test.@testset "kmer_richness" begin
        Test.@testset "DNA sequences" begin
            # Perfect richness (all possible k-mers present)
            dna_seq = BioSequences.LongDNA{4}("ATCGATCGATCGATCG")  # 16 bases
            richness = Mycelia.kmer_richness(dna_seq, 1, normalize = true)
            Test.@test richness ≤ 1.0
            Test.@test richness > 0.0

            # Test unnormalized
            richness_raw = Mycelia.kmer_richness(dna_seq, 1, normalize = false)
            Test.@test richness_raw isa Integer
            Test.@test richness_raw > 0

            # Test sequence shorter than k
            Test.@test Mycelia.kmer_richness(BioSequences.LongDNA{4}("AT"), 3, normalize = true) ==
                       0.0
            Test.@test Mycelia.kmer_richness(BioSequences.LongDNA{4}("AT"), 3, normalize = false) ==
                       0
        end

        Test.@testset "RNA sequences" begin
            rna_seq = BioSequences.LongRNA{4}("AUCGAUCGAUCGAUCG")
            richness = Mycelia.kmer_richness(rna_seq, 1, normalize = true)
            Test.@test richness ≤ 1.0
            Test.@test richness > 0.0
        end

        Test.@testset "Amino acid sequences" begin
            aa_seq = BioSequences.LongAA("ARNDCQEGHILKMFPSTWYV")
            richness = Mycelia.kmer_richness(aa_seq, 1, normalize = true)
            Test.@test richness ≤ 1.0
            Test.@test richness > 0.0

            # Test with custom alphabet - use subset of AA that has < 25 unique residues
            subset_aa = BioSequences.LongAA("ARNDCQEGHILKMF")  # 14 different AAs
            many_positions = repeat(subset_aa, 10)  # 140 positions, 14 unique k-mers for k=1  
            richness_default = Mycelia.kmer_richness(many_positions, 1, normalize = true) # 14/min(140, 20) = 14/20
            richness_custom = Mycelia.kmer_richness(many_positions, 1, alphabet = 25, normalize = true) # 14/min(140, 25) = 14/25
            Test.@test richness_custom != richness_default
        end

        Test.@testset "Generic strings" begin
            str = "abcdefgh"
            richness = Mycelia.kmer_richness(str, 1, normalize = true)
            Test.@test richness ≤ 1.0
            Test.@test richness > 0.0

            # Test repeated characters with explicit alphabet size
            repeated_str = "aaabbbccc"
            richness_repeated = Mycelia.kmer_richness(repeated_str, 1, alphabet = 10, normalize = true)  # Force larger alphabet
            Test.@test richness_repeated < 1.0
        end

        Test.@testset "Edge cases" begin
            # Empty string
            Test.@test Mycelia.kmer_richness("", 1, normalize = true) == 0.0
            Test.@test Mycelia.kmer_richness("", 1, normalize = false) == 0

            # Single character
            Test.@test Mycelia.kmer_richness("a", 1, normalize = true) == 1.0
            Test.@test Mycelia.kmer_richness("a", 2, normalize = true) == 0.0
        end
    end

    Test.@testset "linguistic_complexity" begin
        Test.@testset "DNA sequences" begin
            dna_seq = BioSequences.LongDNA{4}("ATCGATCGATCGATCGATCG")
            profile, summary = Mycelia.linguistic_complexity(dna_seq, kmax = 5)

            Test.@test length(profile) == 5
            Test.@test all(0.0 ≤ x ≤ 1.0 for x in profile)
            Test.@test summary isa Float64
            Test.@test 0.0 ≤ summary ≤ 1.0

            # Test with custom reducer
            _,
            median_summary = Mycelia.linguistic_complexity(dna_seq, kmax = 5, reducer = Statistics.median)
            Test.@test median_summary != summary
        end

        Test.@testset "RNA sequences" begin
            rna_seq = BioSequences.LongRNA{4}("AUCGAUCGAUCGAUCGAUCG")
            profile, summary = Mycelia.linguistic_complexity(rna_seq, kmax = 4)

            Test.@test length(profile) == 4
            Test.@test all(0.0 ≤ x ≤ 1.0 for x in profile)
        end

        Test.@testset "Amino acid sequences" begin
            aa_seq = BioSequences.LongAA("ARNDCQEGHILKMFPSTWYVARNDCQ")
            profile, summary = Mycelia.linguistic_complexity(aa_seq, kmax = 3)

            Test.@test length(profile) == 3
            Test.@test all(0.0 ≤ x ≤ 1.0 for x in profile)

            # Test with custom alphabet
            profile_custom,
            _ = Mycelia.linguistic_complexity(aa_seq, kmax = 3, alphabet = 25)
            Test.@test profile_custom != profile
        end

        Test.@testset "Generic strings" begin
            str = "abcdefghijklmnop"
            profile, summary = Mycelia.linguistic_complexity(str, kmax = 4)

            Test.@test length(profile) == 4
            Test.@test all(0.0 ≤ x ≤ 1.0 for x in profile)

            # Test highly repetitive string
            repetitive = "aaabbbaaabbb"
            profile_rep, summary_rep = Mycelia.linguistic_complexity(repetitive, kmax = 3)
            Test.@test summary_rep < summary  # Should be less complex
        end

        Test.@testset "Parameter effects" begin
            seq = "ababababababab"  # Repetitive pattern to show k-mer complexity differences

            # Default kmax should use full sequence length
            profile_default, _ = Mycelia.linguistic_complexity(seq)
            Test.@test length(profile_default) == length(seq)

            # Custom kmax
            profile_custom, _ = Mycelia.linguistic_complexity(seq, kmax = 5)
            Test.@test length(profile_custom) == 5

            # Different reducers
            _,
            mean_val = Mycelia.linguistic_complexity(seq, kmax = 5, reducer = Statistics.mean)
            _,
            max_val = Mycelia.linguistic_complexity(seq, kmax = 5, reducer = Base.maximum)
            _,
            min_val = Mycelia.linguistic_complexity(seq, kmax = 5, reducer = Base.minimum)

            Test.@test mean_val != max_val
            Test.@test mean_val != min_val
        end

        Test.@testset "Edge cases" begin
            # Very short sequence
            short_seq = "abc"
            profile, summary = Mycelia.linguistic_complexity(short_seq, kmax = 2)
            Test.@test length(profile) == 2

            # Single character sequence
            single_char = "a"
            profile_single, _ = Mycelia.linguistic_complexity(single_char, kmax = 1)
            Test.@test length(profile_single) == 1
            Test.@test profile_single[1] == 1.0

            # Empty sequence (should handle gracefully)
            # This might throw an error or return empty - depends on implementation
            Test.@test_nowarn Mycelia.linguistic_complexity("", kmax = 1)
        end
    end

    Test.@testset "topological_entropy" begin
        Test.@testset "DNA and RNA alphabets" begin
            dna_seq = BioSequences.LongDNA{4}("ACGTACGTACGTACGT")
            dna_default = Mycelia.topological_entropy(dna_seq)
            dna_with_constant = Mycelia.topological_entropy(
                dna_seq, alphabet = Mycelia.DNA_ALPHABET)
            dna_with_symbol = Mycelia.topological_entropy(dna_seq, alphabet = :DNA)

            Test.@test 0.0 ≤ dna_default ≤ 1.0
            Test.@test dna_default ≈ dna_with_constant
            Test.@test dna_default ≈ dna_with_symbol

            rna_seq = BioSequences.LongRNA{4}("AUCGAUCGAUCGAUCG")
            rna_default = Mycelia.topological_entropy(rna_seq)
            rna_with_constant = Mycelia.topological_entropy(
                rna_seq, alphabet = Mycelia.RNA_ALPHABET)
            rna_with_symbol = Mycelia.topological_entropy(rna_seq, alphabet = :RNA)

            Test.@test 0.0 ≤ rna_default ≤ 1.0
            Test.@test rna_default ≈ rna_with_constant
            Test.@test rna_default ≈ rna_with_symbol
        end

        Test.@testset "Amino acids and reduced alphabets" begin
            aa_seq = BioSequences.LongAA("ARNDCEQGHILKMFPSTWYVARND")
            aa_default = Mycelia.topological_entropy(aa_seq)
            aa_with_constant = Mycelia.topological_entropy(
                aa_seq, alphabet = Mycelia.AA_ALPHABET)
            reduced = Mycelia.topological_entropy(
                aa_seq, alphabet = Mycelia.REDUCED_ALPHABET_HP2)

            Test.@test 0.0 ≤ aa_default ≤ 1.0
            Test.@test aa_default ≈ aa_with_constant
            Test.@test 0.0 ≤ reduced ≤ 1.0
        end

        Test.@testset "Unicode and collection alphabets" begin
            unicode_seq = "🙂🙃🙂😉🙂🙃🙂😉"
            unicode_alphabet = ['🙂', '🙃', '😉']
            unicode_entropy = Mycelia.topological_entropy(
                unicode_seq, alphabet = unicode_alphabet)
            unicode_symbol_entropy = Mycelia.topological_entropy(
                unicode_seq, alphabet = :UNICODE)

            Test.@test 0.0 ≤ unicode_entropy ≤ 1.0
            Test.@test 0.0 ≤ unicode_symbol_entropy ≤ 1.0
        end

        Test.@testset "Sliding window behavior" begin
            seq = "abababa"
            entropy_step1 = Mycelia.topological_entropy(
                seq, alphabet = ['a', 'b'], step_size = 1)
            entropy_step2 = Mycelia.topological_entropy(
                seq, alphabet = ['a', 'b'], step_size = 2)

            Test.@test entropy_step1 ≈ 0.5
            Test.@test entropy_step2 ≈ 0.5
        end

        Test.@testset "Parameter validation and edge cases" begin
            Test.@test Mycelia.topological_entropy("") == 0.0
            Test.@test_throws ArgumentError Mycelia.topological_entropy("abc", alphabet = 1)
            Test.@test_throws ArgumentError Mycelia.topological_entropy(
                "abc", alphabet = 3, step_size = 0)

            default_entropy = Mycelia.topological_entropy("aba")
            forced_entropy = Mycelia.topological_entropy("aba", alphabet = 4)
            Test.@test forced_entropy < default_entropy
        end
    end

    Test.@testset "Integration tests" begin
        Test.@testset "Consistency between functions" begin
            # Test that functions give consistent results for same input
            dna_seq = BioSequences.LongDNA{4}("ATCGATCGATCGATCG")

            # Shannon entropy should be related to richness for k=1
            shannon_k1 = Mycelia.shannon_entropy(dna_seq, k = 1)
            richness_k1 = Mycelia.kmer_richness(dna_seq, 1, normalize = false)

            Test.@test shannon_k1 ≥ 0
            Test.@test richness_k1 ≥ 0

            # Renyi entropy should approach Shannon as α approaches 1
            renyi_close_to_1 = Mycelia.renyi_entropy(dna_seq, k = 1, α = 1.001)
            shannon_val = Mycelia.shannon_entropy(dna_seq, k = 1)
            Test.@test abs(renyi_close_to_1 - shannon_val) < 0.1
        end

        Test.@testset "Behavior with increasing k" begin
            seq = "abcdefghijklmnop"

            # Entropy should generally decrease with increasing k (for finite sequences)
            entropies = [Mycelia.shannon_entropy(seq, k = k) for k in 1:5]

            # At least some should be decreasing (though not necessarily monotonic)
            Test.@test any(entropies[i] > entropies[i + 1]
            for i in 1:(length(entropies) - 1))
        end
    end
end
