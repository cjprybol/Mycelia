# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/reduced_amino_acid_alphabets_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/reduced_amino_acid_alphabets_test.jl", "test/2_preprocessing_qc", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import BioSequences
import BioSymbols

Test.@testset "Reduced Amino Acid Alphabets" begin

    Test.@testset "list_reduced_alphabets" begin
        schemes = Mycelia.list_reduced_alphabets()

        ## Should return a vector of symbols
        Test.@test isa(schemes, Vector{Symbol})

        ## Should include all expected legacy schemes
        Test.@test :HP2 in schemes
        Test.@test :HYDROPATHY3 in schemes
        Test.@test :GBMR4 in schemes
        Test.@test :CHEMICAL5 in schemes
        Test.@test :CHEMICAL6 in schemes
        Test.@test :SDM12 in schemes

        ## Should include Murphy schemes
        Test.@test :MURPHY2 in schemes
        Test.@test :MURPHY3 in schemes
        Test.@test :MU4 in schemes
        Test.@test :ML5 in schemes
        Test.@test :MURPHY8 in schemes
        Test.@test :MURPHY10 in schemes
        Test.@test :MURPHY12 in schemes
        Test.@test :MURPHY15 in schemes

        ## Should include Prlić schemes
        Test.@test :SDM12_PRLIC in schemes
        Test.@test :HSDM17 in schemes

        ## Should include Solis schemes
        Test.@test :SOLIS2 in schemes
        Test.@test :SOLIS19 in schemes

        ## Should include other validated schemes
        Test.@test :WW5 in schemes
        Test.@test :MM5 in schemes

        ## Should be sorted
        Test.@test issorted(schemes)

        ## Should have exactly 36 schemes
        Test.@test length(schemes) == 36
    end

    Test.@testset "get_reduced_alphabet_info" begin
        ## Test HP2 scheme info
        info = Mycelia.get_reduced_alphabet_info(:HP2)
        Test.@test info[:name] == "Binary Hydrophobic-Polar"
        Test.@test info[:classes] == 2
        Test.@test haskey(info, :description)
        Test.@test haskey(info, :groups)
        Test.@test info[:groups]['H'] == "ACFILMVW"
        Test.@test info[:groups]['P'] == "GTSYPNDEQKRH"

        ## Test HYDROPATHY3 scheme info
        info = Mycelia.get_reduced_alphabet_info(:HYDROPATHY3)
        Test.@test info[:classes] == 3
        Test.@test haskey(info[:groups], 'H')
        Test.@test haskey(info[:groups], 'N')
        Test.@test haskey(info[:groups], 'P')

        ## Test SDM12 scheme info
        info = Mycelia.get_reduced_alphabet_info(:SDM12)
        Test.@test info[:classes] == 12

        ## Test error on unknown scheme
        Test.@test_throws ArgumentError Mycelia.get_reduced_alphabet_info(:NONEXISTENT)
    end

    Test.@testset "reduce_amino_acid_alphabet - HP2 scheme" begin
        ## Test all 20 amino acids
        seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :HP2)

        ## Expected mapping:
        ## Hydrophobic (H): A, C, F, I, L, M, V, W
        ## Polar (P): G, T, S, Y, P, N, D, E, Q, K, R, H
        expected = "HHPPHPPHPHHPPPPPPHHP"
        Test.@test reduced == expected

        ## Test shorter sequences
        seq = BioSequences.LongAA("ACFIL")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :HP2)
        Test.@test reduced == "HHHHH"

        seq = BioSequences.LongAA("GTSYP")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :HP2)
        Test.@test reduced == "PPPPP"
    end

    Test.@testset "reduce_amino_acid_alphabet - HYDROPATHY3 scheme" begin
        seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :HYDROPATHY3)

        ## Expected mapping:
        ## Hydrophobic (H): I, V, L, F, C, M, A, W
        ## Neutral (N): G, T, S, Y, P, H
        ## Polar (P): D, N, E, Q, K, R
        expected = "HHPPHNNHPHHPNPPNNHHN"
        Test.@test reduced == expected
    end

    Test.@testset "reduce_amino_acid_alphabet - GBMR4 scheme" begin
        seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :GBMR4)

        ## Expected mapping:
        ## G: G
        ## P: P
        ## Hydrophobic (H): Y, F, L, I, M, V, C, W, H
        ## Basic/polar (B): A, D, K, E, R, N, T, S, Q
        expected = "BHBBHGHHBHHBPBBBBHHH"
        Test.@test reduced == expected

        ## Test glycine and proline specifically
        seq = BioSequences.LongAA("GPGP")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :GBMR4)
        Test.@test reduced == "GPGP"
    end

    Test.@testset "reduce_amino_acid_alphabet - CHEMICAL5 scheme" begin
        seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :CHEMICAL5)

        ## Expected mapping:
        ## Aliphatic (A): I, V, L
        ## aRomatic (R): F, Y, W, H
        ## Charged (C): K, R, D, E
        ## Tiny (T): G, A, C, S
        ## Diverse (D): T, M, Q, N, P
        expected = "TTCCRTRACADDDDCTDARR"
        Test.@test reduced == expected
    end

    Test.@testset "reduce_amino_acid_alphabet - CHEMICAL6 scheme" begin
        seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :CHEMICAL6)

        ## Expected mapping:
        ## Aliphatic (A): I, V, L
        ## aRomatic (R): F, Y, W, H
        ## Positive (+): K, R
        ## Negative (-): D, E
        ## Tiny (T): G, A, C, S
        ## Diverse (D): T, M, Q, N, P
        expected = "TT--RTRA+ADDDD+TDARR"
        Test.@test reduced == expected

        ## Test charge separation
        seq = BioSequences.LongAA("DEKR")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :CHEMICAL6)
        Test.@test reduced == "--++"
    end

    Test.@testset "reduce_amino_acid_alphabet - SDM12 scheme" begin
        seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :SDM12)

        ## Expected mapping:
        ## A: A, D: D/E, K: K/R, N: N/Q, T: T/S, Y: Y/F
        ## L: L/I/V/M, C: C, W: W, H: H, G: G, P: P
        expected = "ACDDYGHLKLLNPNKTTLWY"
        Test.@test reduced == expected

        ## Test that similar amino acids map to same letter
        seq = BioSequences.LongAA("DE")  ## Both acidic
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :SDM12)
        Test.@test reduced == "DD"

        seq = BioSequences.LongAA("KR")  ## Both basic
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :SDM12)
        Test.@test reduced == "KK"

        seq = BioSequences.LongAA("LIVM")  ## All aliphatic
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :SDM12)
        Test.@test reduced == "LLLL"
    end

    Test.@testset "reduce_amino_acid_alphabet - Error handling" begin
        seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")

        ## Test unknown scheme
        Test.@test_throws ArgumentError Mycelia.reduce_amino_acid_alphabet(seq, :UNKNOWN)
        Test.@test_throws ArgumentError Mycelia.reduce_amino_acid_alphabet(seq, :HP3)

        ## Test error message contains available schemes
        try
            Mycelia.reduce_amino_acid_alphabet(seq, :INVALID)
            Test.@test false  ## Should not reach here
        catch e
            Test.@test isa(e, ArgumentError)
            Test.@test occursin("HP2", e.msg)
            Test.@test occursin("Available schemes", e.msg)
        end
    end

    Test.@testset "reduce_amino_acid_alphabet - Empty and single amino acid sequences" begin
        ## Test empty sequence
        seq = BioSequences.LongAA("")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :HP2)
        Test.@test reduced == ""

        ## Test single amino acid
        seq = BioSequences.LongAA("A")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :HP2)
        Test.@test reduced == "H"

        seq = BioSequences.LongAA("K")
        reduced = Mycelia.reduce_amino_acid_alphabet(seq, :HP2)
        Test.@test reduced == "P"
    end

    Test.@testset "reduce_amino_acid_alphabet - All schemes coverage" begin
        ## Ensure all 20 standard amino acids are covered in each scheme
        all_aa = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")

        for scheme in Mycelia.list_reduced_alphabets()
            ## Should not throw an error
            reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, scheme)

            ## Reduced sequence should have same length as input
            Test.@test length(reduced) == length(all_aa)

            ## Should only contain valid reduced alphabet characters
            info = Mycelia.get_reduced_alphabet_info(scheme)
            valid_chars = Set(keys(info[:groups]))
            for char in reduced
                Test.@test char in valid_chars
            end
        end
    end

    Test.@testset "reduce_amino_acid_alphabet - Protein example" begin
        ## Test with a realistic protein sequence fragment
        ## This is part of insulin A chain
        insulin_fragment = BioSequences.LongAA("GIVEQCCTSICSLYQLENYCN")

        ## Test HP2 reduction
        reduced_hp2 = Mycelia.reduce_amino_acid_alphabet(insulin_fragment, :HP2)
        Test.@test length(reduced_hp2) == length(insulin_fragment)
        Test.@test all(c in ['H', 'P'] for c in reduced_hp2)

        ## Test SDM12 reduction (more detailed)
        reduced_sdm12 = Mycelia.reduce_amino_acid_alphabet(insulin_fragment, :SDM12)
        Test.@test length(reduced_sdm12) == length(insulin_fragment)
    end

    Test.@testset "Reduced alphabet constants verification" begin
        ## Verify that all 20 standard amino acids are in each alphabet
        standard_aa = [
            BioSymbols.AA_A, BioSymbols.AA_C, BioSymbols.AA_D, BioSymbols.AA_E,
            BioSymbols.AA_F, BioSymbols.AA_G, BioSymbols.AA_H, BioSymbols.AA_I,
            BioSymbols.AA_K, BioSymbols.AA_L, BioSymbols.AA_M, BioSymbols.AA_N,
            BioSymbols.AA_P, BioSymbols.AA_Q, BioSymbols.AA_R, BioSymbols.AA_S,
            BioSymbols.AA_T, BioSymbols.AA_V, BioSymbols.AA_W, BioSymbols.AA_Y
        ]

        for (scheme_name, alphabet_dict) in Mycelia.REDUCED_ALPHABETS
            for aa in standard_aa
                Test.@test haskey(alphabet_dict, aa)
            end

            ## Verify correct number of amino acids mapped
            Test.@test length(alphabet_dict) == 20
        end
    end

    Test.@testset "New Murphy schemes - Expected outputs" begin
        all_aa = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")

        ## MURPHY2
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :MURPHY2)
        Test.@test reduced == "IIEEIIEIEIIEIEEIIIII"
        Test.@test length(reduced) == 20

        ## MURPHY3
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :MURPHY3)
        Test.@test reduced == "LLEEFLELELLELEELLLFF"
        Test.@test length(reduced) == 20

        ## MU4
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :MU4)
        Test.@test reduced == "SLEEFSELELLESEESSLFF"
        Test.@test length(reduced) == 20

        ## ML5
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :ML5)
        Test.@test reduced == "ALEEFAKLKLLEAEKAALFF"
        Test.@test length(reduced) == 20
    end

    Test.@testset "Other validated schemes - Expected outputs" begin
        all_aa = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")

        ## WW5
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :WW5)
        Test.@test reduced == "AIDDIGAIKIIKGKKKAIII"
        Test.@test length(reduced) == 20

        ## MM5
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :MM5)
        Test.@test reduced == "ACDDIAHIDIIDDDDDDIII"
        Test.@test length(reduced) == 20

        ## MURPHY12
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :MURPHY12)
        Test.@test reduced == "ACDEFGHLKLLDPEKSSLWF"
        Test.@test length(reduced) == 20

        ## SDM12_PRLIC
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :SDM12_PRLIC)
        Test.@test reduced == "ACDKYGHLKLLNPTKTTLWY"
        Test.@test length(reduced) == 20

        ## HSDM17
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :HSDM17)
        Test.@test reduced == "ACDKFGHLKLMNPQRSTLWY"
        Test.@test length(reduced) == 20
    end

    Test.@testset "Solis hierarchical set - Sample tests" begin
        all_aa = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")

        ## Test SOLIS2 (minimal reduction)
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :SOLIS2)
        Test.@test length(reduced) == 20
        ## Should only have 2 characters
        Test.@test length(Set(reduced)) == 2

        ## Test SOLIS10 (middle of range)
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :SOLIS10)
        Test.@test length(reduced) == 20
        ## Should have 10 unique characters
        Test.@test length(Set(reduced)) == 10

        ## Test SOLIS19 (nearly complete)
        reduced = Mycelia.reduce_amino_acid_alphabet(all_aa, :SOLIS19)
        Test.@test length(reduced) == 20
        ## Should have 19 unique characters (only Q,S paired)
        Test.@test length(Set(reduced)) == 19
    end

    Test.@testset "Reduced alphabet metadata verification" begin
        ## Verify metadata is consistent with actual mappings
        for scheme in Mycelia.list_reduced_alphabets()
            info = Mycelia.get_reduced_alphabet_info(scheme)
            mapping = Mycelia.REDUCED_ALPHABETS[scheme]

            ## Count unique values in mapping
            unique_classes = Set(values(mapping))
            Test.@test length(unique_classes) == info[:classes]

            ## Verify groups metadata covers all amino acids
            all_aa_in_groups = join(sort(collect(values(info[:groups]))))
            Test.@test length(all_aa_in_groups) >= 20  ## >= because some might be duplicated in metadata
        end
    end

end
