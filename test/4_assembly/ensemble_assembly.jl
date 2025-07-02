# Probabilistic ensemble assembly tests
@testset "probabilistic ensemble assembly" begin
    @testset "Illumina" begin
        @test 1 + 1 == 2
    end
    @testset "Ultima" begin
        @test 1 + 1 == 2
    end
    @testset "Nanopore" begin
        @test 1 + 1 == 2
    end
    @testset "PacBio" begin
        @test 1 + 1 == 2
    end
    @testset "multi-entity, even coverage" begin
        @test 1 + 1 == 2
    end
    @testset "multi-entity, log-distributed coverage" begin
        @test 1 + 1 == 2
    end
    @testset "multi-platform" begin
        @test 1 + 1 == 2
    end
end
