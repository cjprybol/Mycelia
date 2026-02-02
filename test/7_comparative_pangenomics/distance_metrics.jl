# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/distance_metrics.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/distance_metrics.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import Distances
import Clustering
import ProgressMeter
import DataFrames
import SHA
import LinearAlgebra
import BioSequences

Test.@testset "Distance Metrics Tests" begin
    Test.@testset "Jaccard Distance Matrix Functions" begin
        # Test binary matrix Jaccard distance
        binary_matrix = [1 0 1; 0 1 1; 1 1 0]
        jaccard_dist = Mycelia.jaccard_distance(binary_matrix)

        Test.@test jaccard_dist isa Matrix{Float64}
        Test.@test size(jaccard_dist) == (3, 3)
        Test.@test LinearAlgebra.issymmetric(jaccard_dist)
        Test.@test all(LinearAlgebra.diag(jaccard_dist) .== 0.0)  # Diagonal should be zero

        # Manual calculation for verification
        # Column 1: [1,0,1], Column 2: [0,1,1], Column 3: [1,1,0]
        # Jaccard(1,2): intersection=1, union=3, distance=1-1/3=2/3
        Test.@test jaccard_dist[1, 2] ≈ 2/3 atol=1e-10
        Test.@test jaccard_dist[2, 1] ≈ 2/3 atol=1e-10

        # Test Bray-Curtis distance
        count_matrix = [1 2 3; 4 0 2; 0 1 1]
        bc_dist = Mycelia.bray_curtis_distance(count_matrix)

        Test.@test bc_dist isa Matrix{Float64}
        Test.@test size(bc_dist) == (3, 3)
        Test.@test LinearAlgebra.issymmetric(bc_dist)
        Test.@test all(LinearAlgebra.diag(bc_dist) .== 0.0)
    end

    Test.@testset "K-mer Distance Functions" begin
        # Test cosine similarity between k-mer counts
        kmer_counts_1 = Dict("ATG" => 10, "GCT" => 5, "TAA" => 2)
        kmer_counts_2 = Dict("ATG" => 8, "GCT" => 3, "CCC" => 4)

        cosine_dist = Mycelia.kmer_counts_to_cosine_similarity(kmer_counts_1, kmer_counts_2)
        Test.@test cosine_dist isa Float64
        Test.@test 0.0 ≤ cosine_dist ≤ 1.0

        # Test Jensen-Shannon divergence
        js_div = Mycelia.kmer_counts_to_js_divergence(kmer_counts_1, kmer_counts_2)
        Test.@test js_div isa Float64
        Test.@test 0.0 ≤ js_div ≤ 1.0

        # Test identical distributions
        identical_dist = Mycelia.kmer_counts_to_cosine_similarity(kmer_counts_1, kmer_counts_1)
        Test.@test identical_dist ≈ 0.0 atol=1e-10

        identical_js = Mycelia.kmer_counts_to_js_divergence(kmer_counts_1, kmer_counts_1)
        Test.@test identical_js ≈ 0.0 atol=1e-10
    end

    Test.@testset "Set-based Jaccard Functions" begin
        set1 = Set(["A", "B", "C"])
        set2 = Set(["B", "C", "D"])
        set3 = Set(["A", "B", "C"])  # identical to set1

        # Test Jaccard similarity
        jaccard_sim = Mycelia.jaccard_similarity(set1, set2)
        Test.@test jaccard_sim isa Float64
        Test.@test 0.0 ≤ jaccard_sim ≤ 1.0
        Test.@test jaccard_sim ≈ 2/4  # intersection=2, union=4

        # Test identical sets
        identical_sim = Mycelia.jaccard_similarity(set1, set3)
        Test.@test identical_sim ≈ 1.0

        # Test Jaccard distance
        jaccard_dist = Mycelia.jaccard_distance(set1, set2)
        Test.@test jaccard_dist ≈ 1.0 - jaccard_sim
        Test.@test jaccard_dist ≈ 0.5

        # Test disjoint sets
        disjoint_set = Set(["X", "Y", "Z"])
        disjoint_sim = Mycelia.jaccard_similarity(set1, disjoint_set)
        Test.@test disjoint_sim ≈ 0.0

        disjoint_dist = Mycelia.jaccard_distance(set1, disjoint_set)
        Test.@test disjoint_dist ≈ 1.0

        # Test empty sets
        empty_set = Set{String}()
        empty_sim = Mycelia.jaccard_similarity(empty_set, empty_set)
        Test.@test isnan(empty_sim) || empty_sim ≈ 1.0  # 0/0 case handled differently by implementations
    end

    Test.@testset "Matrix Normalization" begin
        # Test distance matrix normalization
        test_matrix = [0.0 0.5 1.0; 0.5 0.0 1.5; 1.0 1.5 0.0]
        normalized = Mycelia.normalize_distance_matrix(test_matrix)

        Test.@test normalized isa Matrix{Float64}
        Test.@test size(normalized) == size(test_matrix)
        Test.@test maximum(normalized) ≈ 1.0

        # Test with matrix containing NaN values
        nan_matrix = [0.0 NaN 1.0; NaN 0.0 0.5; 1.0 0.5 0.0]
        normalized_nan = Mycelia.normalize_distance_matrix(nan_matrix)
        Test.@test maximum(filter(!isnan, normalized_nan)) ≈ 1.0
    end

    Test.@testset "Probability Matrix Conversion" begin
        # Test count matrix to probability matrix conversion
        count_matrix = [10 20 30; 20 10 15; 5 5 10]
        prob_matrix = Mycelia.count_matrix_to_probability_matrix(count_matrix)

        Test.@test prob_matrix isa Matrix{Float64}
        Test.@test size(prob_matrix) == size(count_matrix)

        # Each column should sum to 1.0
        for col in eachcol(prob_matrix)
            Test.@test sum(col) ≈ 1.0 atol=1e-10
        end

        # Test with zero column
        zero_count_matrix = [10 0 30; 20 0 15; 5 0 10]
        zero_prob_matrix = Mycelia.count_matrix_to_probability_matrix(zero_count_matrix)
        Test.@test all(isnan.(zero_prob_matrix[:, 2])) ||
                   all(iszero.(zero_prob_matrix[:, 2]))
    end

    Test.@testset "Pairwise Distance Matrix" begin
        # Test general pairwise distance function
        test_matrix = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

        # Test with Euclidean distance
        euclidean_dist = Mycelia.pairwise_distance_matrix(
            test_matrix,
            dist_func = Distances.euclidean,
            show_progress = false
        )

        Test.@test euclidean_dist isa Matrix{Float64}
        Test.@test size(euclidean_dist) == (3, 3)
        Test.@test LinearAlgebra.issymmetric(euclidean_dist)
        Test.@test all(LinearAlgebra.diag(euclidean_dist) .== 0.0)

        # Test with cosine distance
        cosine_dist = Mycelia.pairwise_distance_matrix(
            test_matrix,
            dist_func = Distances.cosine_dist,
            show_progress = false
        )

        Test.@test cosine_dist isa Matrix{Float64}
        Test.@test size(cosine_dist) == (3, 3)
        Test.@test LinearAlgebra.issymmetric(cosine_dist)
    end

    Test.@testset "Wrapper Functions" begin
        test_matrix = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

        # Test Euclidean distance wrapper
        euclidean_result = Mycelia.frequency_matrix_to_euclidean_distance_matrix(test_matrix)
        Test.@test euclidean_result isa Matrix{Float64}
        Test.@test LinearAlgebra.issymmetric(euclidean_result)

        # Test cosine distance wrapper
        cosine_result = Mycelia.frequency_matrix_to_cosine_distance_matrix(test_matrix)
        Test.@test cosine_result isa Matrix{Float64}
        Test.@test LinearAlgebra.issymmetric(cosine_result)

        # Test Bray-Curtis distance wrapper
        bc_result = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(test_matrix)
        Test.@test bc_result isa Matrix{Float64}
        Test.@test LinearAlgebra.issymmetric(bc_result)

        # Test Jensen-Shannon distance wrapper
        # Normalize matrix first for probability distribution
        prob_matrix = test_matrix ./ sum(test_matrix, dims = 1)
        js_result = Mycelia.frequency_matrix_to_jensen_shannon_distance_matrix(prob_matrix)
        Test.@test js_result isa Matrix{Float64}
        Test.@test LinearAlgebra.issymmetric(js_result)

        # Test binary Jaccard distance wrapper
        binary_matrix = BitMatrix([true false true; false true true; true true false])
        binary_jaccard = Mycelia.binary_matrix_to_jaccard_distance_matrix(binary_matrix)
        Test.@test binary_jaccard isa Matrix{Float64}
        Test.@test LinearAlgebra.issymmetric(binary_jaccard)

        # Test frequency Jaccard distance wrapper
        freq_jaccard = Mycelia.frequency_matrix_to_jaccard_distance_matrix(test_matrix)
        Test.@test freq_jaccard isa Matrix{Float64}
        Test.@test LinearAlgebra.issymmetric(freq_jaccard)
    end

    Test.@testset "Newick Tree Generation" begin
        # Test distance matrix to Newick conversion
        small_dist_matrix = [0.0 0.5 1.0; 0.5 0.0 0.8; 1.0 0.8 0.0]
        labels = ["Sample1", "Sample2", "Sample3"]
        outfile = tempname() * ".newick"

        result_file = Mycelia.distance_matrix_to_newick(
            distance_matrix = small_dist_matrix,
            labels = labels,
            outfile = outfile
        )

        Test.@test result_file == outfile
        Test.@test isfile(outfile)

        # Read and verify Newick format
        newick_content = read(outfile, String)
        Test.@test occursin("(", newick_content)
        Test.@test occursin(")", newick_content)
        Test.@test occursin(";", newick_content)

        # Cleanup
        rm(outfile, force = true)
    end

    Test.@testset "Sequence Hashing" begin
        # Test SHA256 hashing of sequences
        test_seq = "ATCGATCG"
        hash1 = Mycelia.seq2sha256(test_seq)

        Test.@test hash1 isa String
        Test.@test length(hash1) == 64  # SHA256 produces 64-character hex string

        # Test case insensitivity
        hash2 = Mycelia.seq2sha256("atcgatcg")
        Test.@test hash1 == hash2

        # Test with BioSequence
        bio_seq = BioSequences.LongDNA{4}("ATCGATCG")
        hash3 = Mycelia.seq2sha256(bio_seq)
        Test.@test hash1 == hash3

        # Test different sequences produce different hashes
        different_seq = "GCTAGCTA"
        hash4 = Mycelia.seq2sha256(different_seq)
        Test.@test hash1 != hash4
    end

    Test.@testset "Edge Cases and Error Handling" begin
        # Test with single-column matrix
        single_col = reshape([1.0, 2.0, 3.0], 3, 1)
        single_result = Mycelia.pairwise_distance_matrix(single_col, show_progress = false)
        Test.@test size(single_result) == (1, 1)
        Test.@test single_result[1, 1] ≈ 0.0

        # Test with identical columns
        identical_cols = [1.0 1.0 2.0; 2.0 2.0 3.0; 3.0 3.0 4.0]
        identical_result = Mycelia.pairwise_distance_matrix(identical_cols, show_progress = false)
        Test.@test identical_result[1, 2] ≈ 0.0  # Identical columns should have distance 0

        # Test empty k-mer dictionaries
        empty_dict1 = Dict{String, Int}()
        empty_dict2 = Dict{String, Int}()

        # These might throw errors or return NaN, depending on implementation
        Test.@test_nowarn Mycelia.kmer_counts_to_cosine_similarity(empty_dict1, empty_dict2)
        Test.@test_nowarn Mycelia.kmer_counts_to_js_divergence(empty_dict1, empty_dict2)
    end

    Test.@testset "Performance and Memory Tests" begin
        # Test with moderately large matrices to ensure reasonable performance
        large_matrix = rand(100, 50)  # 100 features, 50 samples

        Test.@test_nowarn Mycelia.frequency_matrix_to_euclidean_distance_matrix(large_matrix)
        Test.@test_nowarn Mycelia.frequency_matrix_to_cosine_distance_matrix(large_matrix)

        # Test that results have expected dimensions
        result = Mycelia.frequency_matrix_to_euclidean_distance_matrix(large_matrix)
        Test.@test size(result) == (50, 50)  # Sample × Sample distance matrix
    end

    Test.@testset "Mathematical Properties" begin
        # Test triangle inequality for distance metrics
        test_matrix = [1.0 2.0 3.0 4.0; 5.0 6.0 7.0 8.0; 9.0 10.0 11.0 12.0]
        dist_matrix = Mycelia.frequency_matrix_to_euclidean_distance_matrix(test_matrix)

        # Triangle inequality: d(a,c) ≤ d(a,b) + d(b,c)
        for i in 1:4, j in 1:4, k in 1:4
            Test.@test dist_matrix[i, k] ≤ dist_matrix[i, j] + dist_matrix[j, k] + 1e-10
        end

        # Symmetry
        Test.@test LinearAlgebra.issymmetric(dist_matrix)

        # Non-negativity
        Test.@test all(dist_matrix .≥ 0)

        # Identity: d(x,x) = 0
        Test.@test all(LinearAlgebra.diag(dist_matrix) .≈ 0.0)
    end
end
