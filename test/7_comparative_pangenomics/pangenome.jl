# Pangenome construction tests using simple DNA sequences

import Mycelia
import FASTX
import Kmers
import Graphs

@testset "pangenome construction" begin
    rec1 = Mycelia.random_fasta_record(moltype = :DNA, seed = 1, L = 15)
    rec2 = Mycelia.random_fasta_record(moltype = :DNA, seed = 2, L = 15)

    graph = Mycelia.build_stranded_kmer_graph(Kmers.DNAKmer{3}, [rec1, rec2])
    canonical = Mycelia.count_canonical_kmers(Kmers.DNAKmer{3}, [rec1, rec2])

    @test Graphs.nv(graph) == 2 * length(canonical)
    @test graph.gprops[:k] == 3
    @test length(graph.gprops[:observed_paths]) == 2
end
