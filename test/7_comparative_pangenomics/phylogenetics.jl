# Phylogenetic analyses tests using small DNA sequences

import Mycelia
import FASTX
import Kmers
import Distances

@testset "phylogenetic analyses" begin
    recs = [Mycelia.random_fasta_record(moltype = :DNA, seed = i, L = 12) for i in 1:3]
    counts_list = [Mycelia.count_canonical_kmers(Kmers.DNAKmer{3}, [r]) for r in recs]

    all_kmers = sort(collect(union(vcat(map(keys, counts_list)...))))
    matrix = zeros(Int, length(all_kmers), length(recs))
    for (j, counts) in enumerate(counts_list)
        for (kmer, count) in counts
            i = findfirst(isequal(kmer), all_kmers)
            matrix[i, j] = count
        end
    end

    dist = Mycelia.pairwise_distance_matrix(matrix; dist_func = Distances.euclidean, show_progress = false)

    mktemp() do path, io
        newick_file = Mycelia.distance_matrix_to_newick(distance_matrix = dist, labels = 1:3, outfile = path)
        @test isfile(newick_file)
        line = open(readline, newick_file)
        @test endswith(line, ';')
    end
end
