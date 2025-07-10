# Strain Resolution tests
@testset "Strain Resolution" begin
    @testset "Strain-aware Reassembly" begin
        @test Mycelia.ks(min=5, max=7) == [5, 7]
    end
end
