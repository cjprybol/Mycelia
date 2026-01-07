# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/pangenome.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/pangenome.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Pangenome construction tests using simple DNA sequences
# import Revise

import Test
import Mycelia
import FASTX
import Kmers
import Graphs
import MetaGraphsNext
import BioSequences

Test.@testset "pangenome construction" begin
    rec1 = Mycelia.random_fasta_record(moltype = :DNA, seed = 1, L = 15)
    rec2 = Mycelia.random_fasta_record(moltype = :DNA, seed = 2, L = 15)

    graph = Mycelia.Rhizomorph.build_kmer_graph([rec1, rec2], 3; mode=:doublestrand)

    # Collect observed k-mers from both records (strand-specific) and add reverse complements
    function collect_kmers(record)
        seq = FASTX.sequence(record)
        return [kmer for (kmer, _) in Kmers.UnambiguousDNAMers{3}(seq)]
    end
    unique_kmers = unique(vcat(collect_kmers(rec1), collect_kmers(rec2)))
    expected_labels = Set(vcat(unique_kmers, BioSequences.reverse_complement.(unique_kmers)))

    Test.@test graph isa MetaGraphsNext.MetaGraph
    Test.@test Set(MetaGraphsNext.labels(graph)) == expected_labels
    Test.@test Graphs.nv(graph.graph) == length(expected_labels)

    # Verify k-mer size and that evidence was recorded for the dataset
    first_kmer = first(MetaGraphsNext.labels(graph))
    Test.@test Kmers.ksize(typeof(first_kmer)) == 3
    first_vertex = graph[first_kmer]
    Test.@test "dataset_01" in keys(first_vertex.evidence)
end
