# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/singlestrand_canonicalization_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/singlestrand_canonicalization_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

# Simple test to verify SingleStrand mode fix
Test.@testset "SingleStrand Mode Fix Test (Rhizomorph)" begin
    reference_seq = BioSequences.rna"AUCGAUCGAUCG"
    reads = [FASTX.FASTA.Record("read1", reference_seq)]

    graph = Mycelia.Rhizomorph.build_kmer_graph(
        reads,
        5;
        dataset_id = "ss_fix",
        mode = :singlestrand
    )

    Test.@test graph isa MetaGraphsNext.MetaGraph
    Test.@test !isempty(MetaGraphsNext.labels(graph))

    all_forward = true
    for label in MetaGraphsNext.labels(graph)
        vertex_data = graph[label]
        Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
        strands = Set(obs.strand
        for obs in Iterators.flatten(values(vertex_data.evidence["ss_fix"])))
        if isempty(strands)
            all_forward = false
        else
            all_forward &= all(strand == Mycelia.Rhizomorph.Forward for strand in strands)
        end
    end

    Test.@test all_forward
end
