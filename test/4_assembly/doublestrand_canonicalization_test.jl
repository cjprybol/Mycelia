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
import Kmers
import MetaGraphsNext

# Simple test to verify DoubleStrand mode fix
Test.@testset "DNA DoubleStrand Mode Fix Test" begin
    # Create a simple DNA sequence
    reference_seq = BioSequences.dna"ATCGATCGATCG"
    reads = [FASTX.FASTA.Record("read1", reference_seq)]

    # Build k-mer graph in DoubleStrand mode
    kmer_type = Kmers.DNAKmer{5}
    graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.DoubleStrand)

    Test.@test graph isa MetaGraphsNext.MetaGraph
    Test.@test !isempty(MetaGraphsNext.labels(graph))

    # Check if coverage is now populated (this was the failing test)
    has_coverage = false
    for label in MetaGraphsNext.labels(graph)
        vertex_data = graph[label]
        Test.@test vertex_data isa Mycelia.KmerVertexData
        Test.@test vertex_data.canonical_kmer == label

        # Extract strand orientations from coverage
        strand_orientations = [strand for (obs_id, pos, strand) in vertex_data.coverage]
        if !isempty(strand_orientations)
            has_coverage = true
            println("Found coverage: ", strand_orientations)
        end
    end

    Test.@test has_coverage
end