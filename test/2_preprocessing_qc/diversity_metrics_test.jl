# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/diversity_metrics_test.jl")'
# ```

import Test
import Mycelia
import Statistics

Test.@testset "Alpha Diversity Metrics" begin
    Test.@testset "shannon_diversity" begin
        # Test with uniform distribution (max entropy)
        uniform = [10, 10, 10, 10]
        shannon_uniform = Mycelia.shannon_diversity(uniform)
        # For 4 equally abundant species, H' = ln(4) ≈ 1.386
        Test.@test isapprox(shannon_uniform, log(4), atol = 1e-10)

        # Test with single species (zero entropy)
        single = [100, 0, 0, 0]
        Test.@test Mycelia.shannon_diversity(single) == 0.0

        # Test with empty/zero vector
        empty = [0, 0, 0, 0]
        Test.@test Mycelia.shannon_diversity(empty) == 0.0

        # Test that more even distributions have higher diversity
        even = [25, 25, 25, 25]
        uneven = [90, 5, 3, 2]
        Test.@test Mycelia.shannon_diversity(even) > Mycelia.shannon_diversity(uneven)
    end

    Test.@testset "simpsons_diversity" begin
        # Test with uniform distribution (high diversity)
        uniform = [10, 10, 10, 10]
        simpson_uniform = Mycelia.simpsons_diversity(uniform)
        # For 4 equally abundant species, 1-D = 1 - 4*(0.25)^2 = 0.75
        Test.@test isapprox(simpson_uniform, 0.75, atol = 1e-10)

        # Test with single species (zero diversity)
        single = [100, 0, 0, 0]
        Test.@test Mycelia.simpsons_diversity(single) == 0.0

        # Test with empty vector
        empty = [0, 0, 0, 0]
        Test.@test Mycelia.simpsons_diversity(empty) == 0.0

        # Test range is [0, 1]
        Test.@test 0.0 <= simpson_uniform <= 1.0
    end

    Test.@testset "species_richness" begin
        # Test basic counting
        Test.@test Mycelia.species_richness([10, 5, 3, 0, 0]) == 3
        Test.@test Mycelia.species_richness([1, 1, 1, 1]) == 4
        Test.@test Mycelia.species_richness([0, 0, 0, 0]) == 0
        Test.@test Mycelia.species_richness([100]) == 1

        # Test that abundance magnitude doesn't matter
        Test.@test Mycelia.species_richness([1, 1, 1]) ==
                   Mycelia.species_richness([100, 50, 25])
    end

    Test.@testset "pielous_evenness" begin
        # Test with perfectly even distribution (J' = 1.0)
        even = [10, 10, 10, 10]
        Test.@test isapprox(Mycelia.pielous_evenness(even), 1.0, atol = 1e-10)

        # Test with single species (undefined, returns 0.0)
        single = [100]
        Test.@test Mycelia.pielous_evenness(single) == 0.0

        # Test with uneven distribution (J' < 1.0)
        uneven = [90, 5, 3, 2]
        Test.@test Mycelia.pielous_evenness(uneven) < 1.0
        Test.@test Mycelia.pielous_evenness(uneven) > 0.0

        # Test range is [0, 1]
        evenness = Mycelia.pielous_evenness([10, 10, 10, 10])
        Test.@test 0.0 <= evenness <= 1.0
    end

    Test.@testset "calculate_alpha_diversity" begin
        # Create a simple abundance matrix (taxa × samples)
        abundance = [10 50 0;   # Taxon 1
                     10 25 100; # Taxon 2
                     10 25 0;   # Taxon 3
                     10 0 0]
        samples = ["S1", "S2", "S3"]

        result = Mycelia.calculate_alpha_diversity(abundance, samples)

        # Check that result is a DataFrame with expected columns
        Test.@test "sample" in names(result)
        Test.@test "shannon" in names(result)
        Test.@test "simpsons" in names(result)
        Test.@test "richness" in names(result)
        Test.@test "evenness" in names(result)

        # Check dimensions
        Test.@test size(result, 1) == 3  # 3 samples

        # Check sample names
        Test.@test result.sample == ["S1", "S2", "S3"]

        # S1 has 4 equally abundant species (highest evenness)
        # S2 has 3 species (uneven)
        # S3 has 1 species (lowest diversity)
        Test.@test result.richness[1] == 4
        Test.@test result.richness[2] == 3
        Test.@test result.richness[3] == 1

        # S1 should have highest diversity (even distribution)
        Test.@test result.shannon[1] > result.shannon[2]
        Test.@test result.shannon[1] > result.shannon[3]
    end
end

Test.@testset "Beta Diversity Metrics" begin
    Test.@testset "bray_curtis_dissimilarity" begin
        # Identical samples (distance = 0)
        a = [10, 20, 30]
        Test.@test Mycelia.bray_curtis_dissimilarity(a, a) == 0.0

        # Completely different samples (distance = 1)
        b = [10, 0, 0]
        c = [0, 20, 30]
        Test.@test Mycelia.bray_curtis_dissimilarity(b, c) == 1.0

        # Symmetric
        x = [10, 20, 5]
        y = [15, 10, 25]
        Test.@test Mycelia.bray_curtis_dissimilarity(x, y) ==
                   Mycelia.bray_curtis_dissimilarity(y, x)

        # Range is [0, 1]
        bc = Mycelia.bray_curtis_dissimilarity(x, y)
        Test.@test 0.0 <= bc <= 1.0

        # Empty vectors
        Test.@test Mycelia.bray_curtis_dissimilarity([0, 0], [0, 0]) == 0.0
    end

    Test.@testset "jaccard_distance_vectors" begin
        # Identical presence/absence (distance = 0)
        a = [10, 20, 0]
        Test.@test Mycelia.jaccard_distance_vectors(a, a) == 0.0

        # Completely different (distance = 1)
        b = [10, 0, 0]
        c = [0, 20, 30]
        Test.@test Mycelia.jaccard_distance_vectors(b, c) == 1.0

        # Partial overlap
        x = [10, 20, 0, 0]  # Present at 1, 2
        y = [0, 15, 30, 0]  # Present at 2, 3
        # Intersection = 1 (position 2), Union = 3 (positions 1, 2, 3)
        # Jaccard = 1 - 1/3 = 2/3
        Test.@test isapprox(Mycelia.jaccard_distance_vectors(x, y), 2/3, atol = 1e-10)

        # Symmetric
        Test.@test Mycelia.jaccard_distance_vectors(x, y) ==
                   Mycelia.jaccard_distance_vectors(y, x)

        # Empty vectors
        Test.@test Mycelia.jaccard_distance_vectors([0, 0], [0, 0]) == 0.0
    end

    Test.@testset "calculate_beta_diversity" begin
        # Create abundance matrix
        abundance = [10 50 10;
                     20 25 20;
                     30 25 30]
        samples = ["S1", "S2", "S3"]

        # Bray-Curtis
        result_bc = Mycelia.calculate_beta_diversity(abundance, samples; metric = :bray_curtis)

        Test.@test size(result_bc.distance_matrix) == (3, 3)
        Test.@test result_bc.sample_names == samples

        # Diagonal should be zero
        for i in 1:3
            Test.@test result_bc.distance_matrix[i, i] == 0.0
        end

        # Should be symmetric
        for i in 1:3, j in 1:3

            Test.@test result_bc.distance_matrix[i, j] == result_bc.distance_matrix[j, i]
        end

        # S1 and S3 are identical, so distance should be 0
        Test.@test result_bc.distance_matrix[1, 3] == 0.0

        # Jaccard
        result_jaccard = Mycelia.calculate_beta_diversity(abundance, samples; metric = :jaccard)
        Test.@test size(result_jaccard.distance_matrix) == (3, 3)

        # Invalid metric should error
        Test.@test_throws ErrorException Mycelia.calculate_beta_diversity(abundance, samples; metric = :invalid)
    end

    Test.@testset "frequency_matrix_to_bray_curtis_distance_matrix" begin
        # Test the wrapper function
        abundance = [10 50 10;
                     20 25 20;
                     30 25 30]

        result = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(abundance)

        Test.@test size(result) == (3, 3)

        # Diagonal should be zero
        for i in 1:3
            Test.@test result[i, i] == 0.0
        end

        # S1 and S3 are identical columns
        Test.@test result[1, 3] == 0.0
    end
end

Test.@testset "PCoA and Visualization Helpers" begin
    Test.@testset "pcoa_to_dataframe" begin
        # Create a simple distance matrix
        dist_matrix = [0.0 0.5 0.8;
                       0.5 0.0 0.3;
                       0.8 0.3 0.0]
        samples = ["S1", "S2", "S3"]

        # Run PCoA
        pcoa_result = Mycelia.pcoa_from_dist(dist_matrix; maxoutdim = 2)

        # Convert to DataFrame
        df = Mycelia.pcoa_to_dataframe(pcoa_result, samples)

        Test.@test "sample" in names(df)
        Test.@test "PC1" in names(df)
        Test.@test "PC2" in names(df)
        Test.@test size(df, 1) == 3
        Test.@test df.sample == ["S1", "S2", "S3"]
    end

    Test.@testset "beta_diversity_pcoa" begin
        # Create abundance matrix
        abundance = [10 50 10 80;
                     20 25 20 10;
                     30 25 30 5;
                     40 0 40 5]
        samples = ["S1", "S2", "S3", "S4"]

        result = Mycelia.beta_diversity_pcoa(abundance, samples; metric = :bray_curtis, maxoutdim = 3)

        # Check all expected fields
        Test.@test haskey(result, :distance_matrix)
        Test.@test haskey(result, :pcoa)
        Test.@test haskey(result, :pcoa_df)
        Test.@test haskey(result, :variance_explained)

        # Check distance matrix
        Test.@test size(result.distance_matrix) == (4, 4)

        # Check PCoA DataFrame
        Test.@test "sample" in names(result.pcoa_df)
        Test.@test "PC1" in names(result.pcoa_df)
        Test.@test "PC2" in names(result.pcoa_df)
        Test.@test size(result.pcoa_df, 1) == 4

        # Check variance explained
        Test.@test length(result.variance_explained) >= 1
        Test.@test all(v -> v >= 0, result.variance_explained)
        # Total variance should sum to ~1 (may be slightly off due to numerical precision)
        Test.@test sum(result.variance_explained) <= 1.0 + 1e-6

        # S1 and S3 are identical, should have same coordinates
        s1_idx = findfirst(result.pcoa_df.sample .== "S1")
        s3_idx = findfirst(result.pcoa_df.sample .== "S3")
        Test.@test isapprox(result.pcoa_df.PC1[s1_idx], result.pcoa_df.PC1[s3_idx], atol = 1e-10)
        Test.@test isapprox(result.pcoa_df.PC2[s3_idx], result.pcoa_df.PC2[s3_idx], atol = 1e-10)
    end
end

println("All diversity metrics tests passed!")
