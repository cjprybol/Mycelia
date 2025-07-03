# Reference Graph and K-mer Analysis tests
@testset "Reference Graph and K-mer Analysis" begin
    @testset "Pangenome Construction" begin
        @test true  # placeholder
    end
    @testset "Optimal K-mer Selection" begin
        @test true  # placeholder
    end
    @testset "statistical kmer analyses" begin
        @test Mycelia.optimal_subsequence_length(error_rate=0.001, sequence_length=100, threshold=0.99) == 10
    end
end
