# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/statistical_analysis.jl")'
# ```

import Test
import Distances
import LinearAlgebra
import Mycelia
import PooledArrays
import Random

Test.@testset "Statistical Analysis Tests" begin
    Test.@testset "PCoA coordinate helper" begin
        sample_space = [
            0.0 0.1 5.0 5.1;
            0.0 0.1 5.0 5.1
        ]
        distance_matrix = Distances.pairwise(Distances.Euclidean(), sample_space; dims = 2)

        pcoa_result = Mycelia.pcoa_from_dist(distance_matrix; maxoutdim = 2)
        row_major_result = Mycelia.pcoa_coordinates(distance_matrix; maxoutdim = 2)

        Test.@test row_major_result.coordinates isa Matrix{Float64}
        Test.@test size(row_major_result.coordinates) == (4, 2)
        Test.@test row_major_result.coordinates ≈ Matrix(transpose(pcoa_result.coordinates))
    end

    Test.@testset "PERMANOVA" begin
        sample_space = [
            0.0 0.1 5.0 5.1;
            0.0 0.1 5.0 5.1
        ]
        distance_matrix = Distances.pairwise(Distances.Euclidean(), sample_space; dims = 2)
        groups = PooledArrays.PooledArray(["A", "A", "B", "B"])

        result = Mycelia.permanova(
            distance_matrix,
            groups;
            n_perm = 199,
            rng = Random.MersenneTwister(1)
        )

        Test.@test result.F > 1.0
        Test.@test 0.0 <= result.p <= 1.0
        Test.@test result.R2 > 0.5
        Test.@test result.n_perm == 199
        Test.@test result.df_between == 1
        Test.@test result.df_within == 2
        Test.@test length(result.groups) == 2
        Test.@test result.group_sizes == [2, 2]
        Test.@test LinearAlgebra.issymmetric(distance_matrix)
    end

    Test.@testset "Adjusted Rand Index" begin
        Test.@test Mycelia.adjusted_rand_index([1, 1, 2, 2], [1, 1, 2, 2]) ≈ 1.0
        Test.@test Mycelia.adjusted_rand_index([1, 1, 2, 2], [1, 2, 1, 2]) <= 0.0
    end
end
