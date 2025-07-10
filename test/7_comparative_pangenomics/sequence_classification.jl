# Sequence Classification tests for basic utilities

import Mycelia
import BioSequences
import DataFrames

@testset "sequence classification" begin
    seq1 = BioSequences.dna"ATGC"
    seq2 = BioSequences.dna"GCAT"  # reverse complement of seq1
    seq3 = BioSequences.dna"AAAA"

    @test Mycelia.is_equivalent(seq1, seq2)
    @test !Mycelia.is_equivalent(seq1, seq3)

    df = DataFrames.DataFrame(
        top_taxid = [1, 2],
        top_score = [10.0, 9.0],
        ratio_to_next_best_score = [2.1, 2.5],
        additional_taxids = [Dict{Int, Float64}(), Dict(3 => 8.0)],
    )
    classified = Mycelia.apply_conservative_taxonomy(df; ratio_threshold = 2.0)

    @test all(classified.final_assignment .== df.top_taxid)
    @test all(classified.confidence_level .== "high")
end
