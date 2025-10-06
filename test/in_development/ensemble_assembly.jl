# Probabilistic ensemble assembly tests
import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
Test.@testset "probabilistic ensemble assembly" begin
    Test.@testset "Illumina" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "Ultima" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "Nanopore" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "PacBio" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "multi-entity, even coverage" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "multi-entity, log-distributed coverage" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "multi-platform" begin
        Test.@test 1 + 1 == 2
    end
end
