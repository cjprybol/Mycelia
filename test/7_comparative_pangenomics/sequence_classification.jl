# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/sequence_classification.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/sequence_classification.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Sequence Classification tests for basic utilities
# import Revise

import Test
import Mycelia
import BioSequences
import DataFrames

Test.@testset "sequence classification" begin
    seq1 = BioSequences.dna"ATGC"
    seq2 = BioSequences.dna"GCAT"  # reverse complement of seq1
    seq3 = BioSequences.dna"AAAA"

    Test.@test Mycelia.is_equivalent(seq1, seq2)
    Test.@test !Mycelia.is_equivalent(seq1, seq3)

    df = DataFrames.DataFrame(
        top_taxid = [1, 2],
        top_score = [10.0, 9.0],
        ratio_to_next_best_score = [2.1, 2.5],
        additional_taxids = [Dict{Int, Float64}(), Dict(3 => 8.0)]
    )
    classified = Mycelia.apply_conservative_taxonomy(df; ratio_threshold = 2.0)

    Test.@test all(classified.final_assignment .== df.top_taxid)
    Test.@test all(classified.confidence_level .== "high")
end
