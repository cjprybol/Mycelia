# Pangenome construction tests using simple DNA sequences

# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=test -e 'include("test/7_comparative_pangenomics/pangenome.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=test -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/pangenome.jl", "test/7_comparative_pangenomics", execute=false)'
# ````

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Revise
import Test
import Mycelia
import FASTX
import Kmers
import Graphs

Test.@testset "pangenome construction" begin
    rec1 = Mycelia.random_fasta_record(moltype = :DNA, seed = 1, L = 15)
    rec2 = Mycelia.random_fasta_record(moltype = :DNA, seed = 2, L = 15)

    graph = Mycelia.build_stranded_kmer_graph(Kmers.DNAKmer{3}, [rec1, rec2])
    canonical = Mycelia.count_canonical_kmers(Kmers.DNAKmer{3}, [rec1, rec2])

    Test.@test Graphs.nv(graph) == 2 * length(canonical)
    Test.@test graph.gprops[:k] == 3
    Test.@test length(graph.gprops[:observed_paths]) == 2
end
