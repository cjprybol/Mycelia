# Sequence IO tests
import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia

Test.@testset "sequence IO" begin
    Test.@testset "detect alphabet type" begin
        Test.@test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA))) == :DNA
        Test.@test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA))) == :RNA
        Test.@test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA))) == :AA
    end
    Test.@testset "autoconvert sequences" begin
        Test.@test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA)))) <: BioSequences.LongDNA{4}
        Test.@test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA)))) <: BioSequences.LongRNA{4}
        Test.@test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA)))) <: BioSequences.LongAA
    end
    Test.@testset "detect sequence extension" begin
        Test.@test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA))) == ".fna"
        Test.@test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA))) == ".frn"
        Test.@test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA))) == ".faa"
    end
end
