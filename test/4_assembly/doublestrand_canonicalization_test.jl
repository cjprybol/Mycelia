# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/4_assembly/doublestrand_canonicalization_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/doublestrand_canonicalization_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

# Simple test to verify DoubleStrand mode fix
Test.@testset "DNA DoubleStrand Mode Fix Test (Rhizomorph)" begin
    reference_seq = BioSequences.dna"ATCGATCGATCG"
    reads = [FASTX.FASTA.Record("read1", reference_seq)]

    graph = Mycelia.Rhizomorph.build_kmer_graph(
        reads,
        5;
        dataset_id="ds_fix",
        mode=:doublestrand,
    )

    Test.@test graph isa MetaGraphsNext.MetaGraph
    Test.@test !isempty(MetaGraphsNext.labels(graph))

    has_forward_and_reverse = false
    for label in MetaGraphsNext.labels(graph)
        vertex_data = graph[label]
        Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
        strands = Set(obs.orientation for obs in Iterators.flatten(values(vertex_data.evidence["ds_fix"])))
        if Mycelia.Rhizomorph.Forward in strands && Mycelia.Rhizomorph.Reverse in strands
            has_forward_and_reverse = true
        end
    end

    Test.@test has_forward_and_reverse
end
