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

# Benchmarking
using PkgBenchmark
# compares this current state against master branch
results = judge(Eisenia, "master")
export_markdown(stdout, results, export_invariants=true)
