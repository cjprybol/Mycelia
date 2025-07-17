# Pantranscriptome construction tests with RNA sequences

# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=test -e 'include("test/7_comparative_pangenomics/pantranscriptome.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=test -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/pantranscriptome.jl", "test/7_comparative_pangenomics", execute=false)'
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
import Distances

Test.@testset "pantranscriptome construction" begin
    recs = [Mycelia.random_fasta_record(moltype = :RNA, seed = i, L = 12) for i in 1:3]
    counts_list = [Mycelia.count_canonical_kmers(Kmers.RNAKmer{3}, [r]) for r in recs]

    all_kmers = sort(collect(union(vcat(map(keys, counts_list)...))))
    matrix = zeros(Int, length(all_kmers), length(recs))
    for (j, counts) in enumerate(counts_list)
        for (kmer, count) in counts
            i = findfirst(isequal(kmer), all_kmers)
            matrix[i, j] = count
        end
    end

    dist = Mycelia.pairwise_distance_matrix(matrix; dist_func = Distances.euclidean, show_progress = false)

    Test.@test size(dist) == (3, 3)
    Test.@test dist[1, 1] == 0
    Test.@test dist[1, 2] == dist[2, 1]
end
