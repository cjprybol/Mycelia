# Panproteome construction tests on toy protein sequences

import Mycelia
import FASTX
import Kmers

@testset "panproteome construction" begin
    recs = [
        Mycelia.random_fasta_record(moltype = :AA, seed = 1, L = 10),
        Mycelia.random_fasta_record(moltype = :AA, seed = 2, L = 10),
    ]

    counts = Mycelia.count_canonical_kmers(Kmers.AAKmer{2}, recs)
    expected = sum(length(FASTX.sequence(r)) - 1 for r in recs)

    @test sum(values(counts)) == expected
    @test length(counts) <= 400
end
