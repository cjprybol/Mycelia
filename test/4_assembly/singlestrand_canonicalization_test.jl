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
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

# Simple test to verify SingleStrand mode fix
Test.@testset "SingleStrand Mode Fix Test" begin
    # Create a simple RNA sequence
    reference_seq = BioSequences.rna"AUCGAUCGAUCG"
    reads = [FASTX.FASTA.Record("read1", reference_seq)]

    # Build k-mer graph in SingleStrand mode
    kmer_type = Kmers.RNAKmer{5}
    graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.SingleStrand)

    Test.@test graph isa MetaGraphsNext.MetaGraph
    Test.@test !isempty(MetaGraphsNext.labels(graph))

    # Check if coverage is now populated (this was failing due to k-mer mismatches)
    has_coverage = false
    for label in MetaGraphsNext.labels(graph)
        vertex_data = graph[label]
        Test.@test vertex_data isa Mycelia.KmerVertexData

        # Extract strand orientations from coverage
        strand_orientations = [strand for (obs_id, pos, strand) in vertex_data.coverage]
        if !isempty(strand_orientations)
            has_coverage = true
            println("Found coverage: ", strand_orientations)
            # In SingleStrand mode, all orientations should be Forward
            Test.@test all(s == Mycelia.Forward for s in strand_orientations)
        end
    end

    Test.@test has_coverage
end