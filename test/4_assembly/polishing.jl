# Assembly Polishing tests
@testset "Assembly Polishing" begin
    @testset "Error Correction" begin
        solidity = Bool[true, false, false, true, false, true]
        branchpoints = [1, 4, 6]
        stretches = Mycelia.find_resampling_stretches(
            record_kmer_solidity=solidity,
            solid_branching_kmer_indices=branchpoints,
        )
        @test stretches == [1:4, 4:6]
    end
end
