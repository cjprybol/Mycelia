# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/alphabets.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/alphabets.jl", "test/2_preprocessing_qc", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import BioSequences

Test.@testset "Alphabet Detection and Conversion Tests" begin
    Test.@testset "get_biosequence_alphabet function" begin
        # Test with different BioSequence types
        dna_seq = BioSequences.LongDNA{4}("ATCG")
        rna_seq = BioSequences.LongRNA{4}("AUCG")
        aa_seq = BioSequences.LongAA("MKTL")
        
        Test.@test Mycelia.get_biosequence_alphabet(dna_seq) == BioSequences.DNAAlphabet{4}
        Test.@test Mycelia.get_biosequence_alphabet(rna_seq) == BioSequences.RNAAlphabet{4}
        Test.@test Mycelia.get_biosequence_alphabet(aa_seq) == BioSequences.AminoAcidAlphabet
    end

    Test.@testset "detect_alphabet function - string input" begin
        # Test DNA sequences
        Test.@test Mycelia.detect_alphabet("ATCG") == :DNA
        Test.@test Mycelia.detect_alphabet("atcg") == :DNA
        Test.@test Mycelia.detect_alphabet("ATCGN") == :DNA  # Test ACGTN detection
        Test.@test Mycelia.detect_alphabet("ATCGNatcgn") == :DNA
        Test.@test Mycelia.detect_alphabet("ATCGRYSWKMBDH") == :DNA
        
        # Test RNA sequences
        Test.@test Mycelia.detect_alphabet("AUCG") == :RNA
        Test.@test Mycelia.detect_alphabet("aucg") == :RNA
        Test.@test Mycelia.detect_alphabet("AUCGN") == :RNA  # Test ACGUN detection
        Test.@test Mycelia.detect_alphabet("AUCGNaucgn") == :RNA
        Test.@test Mycelia.detect_alphabet("AUCGRYSWKMBDH") == :RNA
        
        # Test protein sequences
        Test.@test Mycelia.detect_alphabet("MKTLVF") == :AA
        Test.@test Mycelia.detect_alphabet("PROTEIN") == :AA
        Test.@test Mycelia.detect_alphabet("ATCGX") == :AA  # X is not in nucleotide alphabet
        Test.@test_throws ArgumentError Mycelia.detect_alphabet("123")  # Numbers not in any alphabet
        
        # Test edge cases and hierarchy
        Test.@test Mycelia.detect_alphabet("ACGN") == :DNA  # No T or U, fits DNA+N
        Test.@test Mycelia.detect_alphabet("ACG") == :DNA   # Unambiguous DNA comes before RNA
        Test.@test Mycelia.detect_alphabet("ACGU") == :RNA  # Contains U, must be RNA
        Test.@test Mycelia.detect_alphabet("ACGT") == :DNA  # Contains T, must be DNA
        Test.@test Mycelia.detect_alphabet("ACN") == :DNA   # Fits ACGTN before ACGUN
        Test.@test Mycelia.detect_alphabet("ACUN") == :RNA  # Contains U, must be RNA
        
        # Test sequences that contain both T and U should be classified as AA
        Test.@test Mycelia.detect_alphabet("ATUG") == :AA  # Both T and U, not valid nucleotide
        Test.@test Mycelia.detect_alphabet("AtuG") == :AA  # Both t and U, not valid nucleotide
        Test.@test Mycelia.detect_alphabet("ATCGAAU") == :AA  # Mixed T and U, not valid nucleotide
    end

    Test.@testset "detect_alphabet function - BioSequence input" begin
        # Test with BioSequences
        dna_seq = BioSequences.LongDNA{4}("ATCG")
        rna_seq = BioSequences.LongRNA{4}("AUCG")
        aa_seq = BioSequences.LongAA("MKTL")
        
        Test.@test Mycelia.detect_alphabet(dna_seq) == :DNA
        Test.@test Mycelia.detect_alphabet(rna_seq) == :RNA
        Test.@test Mycelia.detect_alphabet(aa_seq) == :AA
    end

    Test.@testset "convert_sequence function" begin
        # Test DNA conversion
        dna_result = Mycelia.convert_sequence("ATCG")
        Test.@test dna_result isa BioSequences.LongDNA{4}
        Test.@test string(dna_result) == "ATCG"
        
        # Test RNA conversion
        rna_result = Mycelia.convert_sequence("AUCG")
        Test.@test rna_result isa BioSequences.LongRNA{4}
        Test.@test string(rna_result) == "AUCG"
        
        # Test protein conversion
        aa_result = Mycelia.convert_sequence("MKTL")
        Test.@test aa_result isa BioSequences.LongAA
        Test.@test string(aa_result) == "MKTL"
        
        # Test case insensitive conversion
        dna_lower = Mycelia.convert_sequence("atcg")
        Test.@test dna_lower isa BioSequences.LongDNA{4}
        
        # Test with ambiguity codes
        dna_ambiguous = Mycelia.convert_sequence("ATCGN")
        Test.@test dna_ambiguous isa BioSequences.LongDNA{4}
        Test.@test string(dna_ambiguous) == "ATCGN"
        
        # Test mixed case
        mixed_dna = Mycelia.convert_sequence("AtCgT")
        Test.@test mixed_dna isa BioSequences.LongDNA{4}
        
        # Test longer sequences
        long_dna = Mycelia.convert_sequence("ATCGATCGATCGATCG")
        Test.@test long_dna isa BioSequences.LongDNA{4}
        Test.@test length(long_dna) == 16
        
        # Test sequences with mixed T and U become AA
        atug_result = Mycelia.convert_sequence("ATUG")  # Mixed T and U becomes AA
        Test.@test atug_result isa BioSequences.LongAA
    end

    Test.@testset "Edge cases and error handling" begin
        # Test empty string should throw error
        Test.@test_throws ArgumentError Mycelia.detect_alphabet("")
        Test.@test_throws ArgumentError Mycelia.convert_sequence("")
        
        # Test single character sequences
        Test.@test Mycelia.detect_alphabet("A") == :DNA  # A is in DNA (comes first)
        Test.@test Mycelia.detect_alphabet("T") == :DNA  # T is DNA-specific
        Test.@test Mycelia.detect_alphabet("U") == :RNA  # U is RNA-specific
        Test.@test Mycelia.detect_alphabet("N") == :DNA  # N fits ACGTN first
        Test.@test Mycelia.detect_alphabet("X") == :AA   # X is AA-specific
        
        # Test sequences with gaps/unknown characters
        Test.@test Mycelia.detect_alphabet("AT-CG") == :DNA  # Gap character is filtered out in alphabet constants
        Test.@test_throws ArgumentError Mycelia.detect_alphabet("ATCG?")  # Unknown character throws error
        
        # Test very long sequences
        long_seq = repeat("ATCG", 1000)
        Test.@test Mycelia.detect_alphabet(long_seq) == :DNA
        long_result = Mycelia.convert_sequence(long_seq)
        Test.@test long_result isa BioSequences.LongDNA{4}
        Test.@test length(long_result) == 4000
    end

    Test.@testset "Consistency tests" begin
        # Test that conversion followed by alphabet detection is consistent
        test_sequences = [
            ("ATCG", :DNA),
            ("AUCG", :RNA), 
            ("MKTL", :AA),
            ("atcgn", :DNA),
            ("aucgn", :RNA),
            ("ATCGN", :DNA),  # Test ACGTN
            ("AUCGN", :RNA),  # Test ACGUN
            ("ACG", :DNA),    # Test hierarchy: DNA comes before RNA
            ("ACGT", :DNA),   # Test T-containing sequence
            ("ACGU", :RNA)    # Test U-containing sequence
        ]
        
        for (seq, expected_alphabet) in test_sequences
            detected = Mycelia.detect_alphabet(seq)
            Test.@test detected == expected_alphabet
            
            converted = Mycelia.convert_sequence(seq)
            detected_after = Mycelia.detect_alphabet(converted)
            Test.@test detected_after == expected_alphabet
        end
    end

    Test.@testset "Performance and memory tests" begin
        # Test with moderately large sequences to ensure no memory issues
        large_dna = repeat("ATCG", 10000)
        large_rna = repeat("AUCG", 10000)
        large_aa = repeat("MKTL", 10000)
        
        Test.@test Mycelia.detect_alphabet(large_dna) == :DNA
        Test.@test Mycelia.detect_alphabet(large_rna) == :RNA
        Test.@test Mycelia.detect_alphabet(large_aa) == :AA
        
        # Conversion should work without memory issues
        converted_dna = Mycelia.convert_sequence(large_dna)
        converted_rna = Mycelia.convert_sequence(large_rna)
        converted_aa = Mycelia.convert_sequence(large_aa)
        
        Test.@test length(converted_dna) == 40000
        Test.@test length(converted_rna) == 40000
        Test.@test length(converted_aa) == 40000
    end

    Test.@testset "ACGTN and ACGUN alphabet detection" begin
        # Test that N-containing sequences are properly detected
        Test.@test Mycelia.detect_alphabet("NNNNN") == :DNA  # N alone should be DNA (ACGTN first)
        Test.@test Mycelia.detect_alphabet("ATCGNNNN") == :DNA
        Test.@test Mycelia.detect_alphabet("AUCGNNNN") == :RNA
        
        # Test sequences that fit ACGTN/ACGUN but not unambiguous alphabets
        Test.@test Mycelia.detect_alphabet("ACGN") == :DNA
        Test.@test Mycelia.detect_alphabet("ACUN") == :RNA
        
        # Test mixed case with N
        Test.@test Mycelia.detect_alphabet("atcgN") == :DNA
        Test.@test Mycelia.detect_alphabet("aucgN") == :RNA
        
        # Verify conversion works for ACGTN/ACGUN sequences
        dna_n_result = Mycelia.convert_sequence("ATCGN")
        Test.@test dna_n_result isa BioSequences.LongDNA{4}
        Test.@test string(dna_n_result) == "ATCGN"
        
        rna_n_result = Mycelia.convert_sequence("AUCGN")
        Test.@test rna_n_result isa BioSequences.LongRNA{4}
        Test.@test string(rna_n_result) == "AUCGN"
    end
end
