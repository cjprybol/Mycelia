using Eisenia
using Test
using Documenter

@testset "Eisenia.jl" begin
    @testset "dummy-tests" begin
        @test 1 + 1 == 2
    end
    
    @testset "documentation" begin
        # https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Doctesting-as-Part-of-Testing
#         doctest(Eisenia; manual = false)
        doctest(Eisenia)
    end
end
