# Sequence IO tests
@testset "sequence IO" begin
    @testset "detect alphabet type" begin
        @test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA))) == :DNA
        @test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA))) == :RNA
        @test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA))) == :AA
    end
    @testset "autoconvert sequences" begin
        @test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA)))) <: BioSequences.LongDNA{4}
        @test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA)))) <: BioSequences.LongRNA{4}
        @test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA)))) <: BioSequences.LongAA
    end
    @testset "detect sequence extension" begin
        @test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA))) == ".fna"
        @test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA))) == ".frn"
        @test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA))) == ".faa"
    end
end
