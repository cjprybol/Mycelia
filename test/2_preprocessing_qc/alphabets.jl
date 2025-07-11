using Test
using Mycelia
using BioSequences

@testset "Alphabet Detection and Conversion Tests" begin
    @testset "get_biosequence_alphabet function" begin
        # Test with different BioSequence types
        dna_seq = BioSequences.LongDNA{4}("ATCG")
        rna_seq = BioSequences.LongRNA{4}("AUCG")
        aa_seq = BioSequences.LongAA("MKTL")
        
        @test Mycelia.get_biosequence_alphabet(dna_seq) == BioSequences.DNAAlphabet{4}
        @test Mycelia.get_biosequence_alphabet(rna_seq) == BioSequences.RNAAlphabet{4}
        @test Mycelia.get_biosequence_alphabet(aa_seq) == BioSequences.AminoAcidAlphabet
    end

    @testset "detect_alphabet function - string input" begin
        # Test DNA sequences
        @test Mycelia.detect_alphabet("ATCG") == :DNA
        @test Mycelia.detect_alphabet("atcg") == :DNA
        @test Mycelia.detect_alphabet("ATCGNatcgn") == :DNA
        @test Mycelia.detect_alphabet("ATCGRYSWKMBDH") == :DNA
        @test Mycelia.detect_alphabet("") == :DNA  # default case
        
        # Test RNA sequences
        @test Mycelia.detect_alphabet("AUCG") == :RNA
        @test Mycelia.detect_alphabet("aucg") == :RNA
        @test Mycelia.detect_alphabet("AUCGNaucgn") == :RNA
        @test Mycelia.detect_alphabet("AUCGRYSWKMBDH") == :RNA
        
        # Test protein sequences
        @test Mycelia.detect_alphabet("MKTLVF") == :AA
        @test Mycelia.detect_alphabet("PROTEIN") == :AA
        @test Mycelia.detect_alphabet("ATCGX") == :AA  # X is not in nucleotide alphabet
        @test Mycelia.detect_alphabet("123") == :AA  # Numbers not in nucleotide alphabet
        
        # Test edge cases
        @test Mycelia.detect_alphabet("ACGN") == :DNA  # No T or U, defaults to DNA
        
        # Test error cases
        @test_throws ArgumentError Mycelia.detect_alphabet("ATUG")  # Both T and U
        @test_throws ArgumentError Mycelia.detect_alphabet("AtuG")  # Both t and U
        @test_throws ArgumentError Mycelia.detect_alphabet("ATCGAAU")  # Mixed T and U
    end

    @testset "detect_alphabet function - BioSequence input" begin
        # Test with BioSequences
        dna_seq = BioSequences.LongDNA{4}("ATCG")
        rna_seq = BioSequences.LongRNA{4}("AUCG")
        aa_seq = BioSequences.LongAA("MKTL")
        
        @test Mycelia.detect_alphabet(dna_seq) == :DNA
        @test Mycelia.detect_alphabet(rna_seq) == :RNA
        @test Mycelia.detect_alphabet(aa_seq) == :AA
    end

    @testset "convert_sequence function" begin
        # Test DNA conversion
        dna_result = Mycelia.convert_sequence("ATCG")
        @test dna_result isa BioSequences.LongDNA{4}
        @test string(dna_result) == "ATCG"
        
        # Test RNA conversion
        rna_result = Mycelia.convert_sequence("AUCG")
        @test rna_result isa BioSequences.LongRNA{4}
        @test string(rna_result) == "AUCG"
        
        # Test protein conversion
        aa_result = Mycelia.convert_sequence("MKTL")
        @test aa_result isa BioSequences.LongAA
        @test string(aa_result) == "MKTL"
        
        # Test case insensitive conversion
        dna_lower = Mycelia.convert_sequence("atcg")
        @test dna_lower isa BioSequences.LongDNA{4}
        
        # Test with ambiguity codes
        dna_ambiguous = Mycelia.convert_sequence("ATCGN")
        @test dna_ambiguous isa BioSequences.LongDNA{4}
        @test string(dna_ambiguous) == "ATCGN"
        
        # Test mixed case
        mixed_dna = Mycelia.convert_sequence("AtCgT")
        @test mixed_dna isa BioSequences.LongDNA{4}
        
        # Test longer sequences
        long_dna = Mycelia.convert_sequence("ATCGATCGATCGATCG")
        @test long_dna isa BioSequences.LongDNA{4}
        @test length(long_dna) == 16
        
        # Test error conditions
        @test_throws ArgumentError Mycelia.convert_sequence("ATUG")  # Mixed T and U
    end

    @testset "Edge cases and error handling" begin
        # Test empty string
        @test Mycelia.detect_alphabet("") == :DNA
        empty_result = Mycelia.convert_sequence("")
        @test empty_result isa BioSequences.LongDNA{4}
        @test length(empty_result) == 0
        
        # Test single character sequences
        @test Mycelia.detect_alphabet("A") == :DNA
        @test Mycelia.detect_alphabet("U") == :RNA
        @test Mycelia.detect_alphabet("X") == :AA
        
        # Test sequences with gaps/unknown characters
        @test Mycelia.detect_alphabet("AT-CG") == :AA  # Gap character makes it protein
        @test Mycelia.detect_alphabet("ATCG?") == :AA  # Unknown character
        
        # Test very long sequences
        long_seq = repeat("ATCG", 1000)
        @test Mycelia.detect_alphabet(long_seq) == :DNA
        long_result = Mycelia.convert_sequence(long_seq)
        @test long_result isa BioSequences.LongDNA{4}
        @test length(long_result) == 4000
    end

    @testset "Consistency tests" begin
        # Test that conversion followed by alphabet detection is consistent
        test_sequences = [
            ("ATCG", :DNA),
            ("AUCG", :RNA), 
            ("MKTL", :AA),
            ("atcgn", :DNA),
            ("aucgn", :RNA)
        ]
        
        for (seq, expected_alphabet) in test_sequences
            detected = Mycelia.detect_alphabet(seq)
            @test detected == expected_alphabet
            
            converted = Mycelia.convert_sequence(seq)
            detected_after = Mycelia.detect_alphabet(converted)
            @test detected_after == expected_alphabet
        end
    end

    @testset "Performance and memory tests" begin
        # Test with moderately large sequences to ensure no memory issues
        large_dna = repeat("ATCG", 10000)
        large_rna = repeat("AUCG", 10000)
        large_aa = repeat("MKTL", 10000)
        
        @test Mycelia.detect_alphabet(large_dna) == :DNA
        @test Mycelia.detect_alphabet(large_rna) == :RNA
        @test Mycelia.detect_alphabet(large_aa) == :AA
        
        # Conversion should work without memory issues
        converted_dna = Mycelia.convert_sequence(large_dna)
        converted_rna = Mycelia.convert_sequence(large_rna)
        converted_aa = Mycelia.convert_sequence(large_aa)
        
        @test length(converted_dna) == 40000
        @test length(converted_rna) == 40000
        @test length(converted_aa) == 40000
    end
end