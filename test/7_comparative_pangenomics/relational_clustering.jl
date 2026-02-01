# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/relational_clustering.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/relational_clustering.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import DataFrames
import Statistics
import LinearAlgebra
import StableRNGs

Test.@testset "Relational Clustering Tests" begin
    Test.@testset "RelationalMatrix Construction" begin
        ## Test long_to_relational_matrix with basic data
        df = DataFrames.DataFrame(
            entity_a = ["A1", "A1", "A2", "A2", "A3"],
            entity_b = ["B1", "B2", "B1", "B2", "B1"],
            value = [1.0, 2.0, 3.0, 4.0, 5.0]
        )

        rm = Mycelia.long_to_relational_matrix(
            df, :entity_a, :entity_b, :value;
            entity_a_name = "test_a",
            entity_b_name = "test_b",
            measurement_name = "score"
        )

        Test.@test rm isa Mycelia.RelationalMatrix{Float64}
        Test.@test size(rm) == (3, 2)
        Test.@test rm.entity_a_name == "test_a"
        Test.@test rm.entity_b_name == "test_b"
        Test.@test rm.measurement_name == "score"
        Test.@test rm["A1", "B1"] == 1.0
        Test.@test rm["A2", "B2"] == 4.0
        Test.@test isnan(rm["A3", "B2"])  ## Missing value

        ## Test accessor functions
        Test.@test Mycelia.entity_a_ids(rm) == ["A1", "A2", "A3"]
        Test.@test Mycelia.entity_b_ids(rm) == ["B1", "B2"]
        Test.@test Mycelia.n_entity_a(rm) == 3
        Test.@test Mycelia.n_entity_b(rm) == 2
        Test.@test Mycelia.n_missing(rm) == 1
        Test.@test Mycelia.n_filled(rm) == 5
        Test.@test Mycelia.coverage(rm) ≈ 5/6
    end

    Test.@testset "RelationalMatrix with Aggregation" begin
        ## Test that duplicate entries are aggregated
        df = DataFrames.DataFrame(
            entity_a = ["A1", "A1", "A1"],
            entity_b = ["B1", "B1", "B1"],
            value = [1.0, 2.0, 3.0]
        )

        rm = Mycelia.long_to_relational_matrix(df, :entity_a, :entity_b, :value)

        ## Default aggregation is median
        Test.@test rm["A1", "B1"] == 2.0  ## median of [1, 2, 3]
    end

    Test.@testset "RelationalMatrix Transpose" begin
        df = DataFrames.DataFrame(
            entity_a = ["A1", "A2"],
            entity_b = ["B1", "B2"],
            value = [1.0, 2.0]
        )

        rm = Mycelia.long_to_relational_matrix(
            df, :entity_a, :entity_b, :value;
            entity_a_name = "row_type",
            entity_b_name = "col_type"
        )
        rm_t = Mycelia.transpose_relational(rm)

        Test.@test size(rm_t) == (2, 2)
        Test.@test Mycelia.entity_a_ids(rm_t) == Mycelia.entity_b_ids(rm)
        Test.@test Mycelia.entity_b_ids(rm_t) == Mycelia.entity_a_ids(rm)
        Test.@test rm_t.entity_a_name == "col_type"
        Test.@test rm_t.entity_b_name == "row_type"
    end

    Test.@testset "RelationalMatrix Summary" begin
        rm = Mycelia.RelationalMatrix{Float64}(
            [1.0 2.0; 3.0 4.0; 5.0 NaN],
            ["A1", "A2", "A3"],
            ["B1", "B2"],
            "entity_a", "entity_b", "value",
            nothing, nothing
        )

        stats = Mycelia.summarize_relational_matrix(rm; verbose = false)

        Test.@test stats.n_entity_a == 3
        Test.@test stats.n_entity_b == 2
        Test.@test stats.n_valid_values == 5
        Test.@test stats.n_missing_values == 1
        Test.@test stats.pct_filled ≈ 83.33 atol=0.01
    end

    Test.@testset "RelationalMatrix to Long Format" begin
        rm = Mycelia.RelationalMatrix{Float64}(
            [1.0 2.0; 3.0 4.0],
            ["A1", "A2"],
            ["B1", "B2"],
            "entity_a", "entity_b", "value",
            nothing, nothing
        )

        long_df = Mycelia.relational_matrix_to_long(rm; include_missing = true)

        Test.@test DataFrames.nrow(long_df) == 4
        Test.@test :entity_a in DataFrames.propertynames(long_df)
        Test.@test :entity_b in DataFrames.propertynames(long_df)
        Test.@test :value in DataFrames.propertynames(long_df)
    end

    Test.@testset "Imputation Methods" begin
        ## Create distance matrix with missing values
        dist = [0.0 0.5 NaN; 0.5 0.0 0.8; NaN 0.8 0.0]

        Test.@testset "IMPUTE_MAX" begin
            imputed = Mycelia.impute_distances(dist; method = Mycelia.IMPUTE_MAX)
            Test.@test imputed[1, 3] == 1.0
            Test.@test imputed[3, 1] == 1.0
            Test.@test LinearAlgebra.issymmetric(imputed)
        end

        Test.@testset "IMPUTE_MAX_OBSERVED" begin
            imputed = Mycelia.impute_distances(dist; method = Mycelia.IMPUTE_MAX_OBSERVED)
            Test.@test imputed[1, 3] == 0.8  ## max of [0.5, 0.8]
            Test.@test imputed[3, 1] == 0.8
        end

        Test.@testset "IMPUTE_MEDIAN" begin
            imputed = Mycelia.impute_distances(dist; method = Mycelia.IMPUTE_MEDIAN)
            expected_median = Statistics.median([0.5, 0.8])
            Test.@test imputed[1, 3] ≈ expected_median
        end

        Test.@testset "IMPUTE_MEAN" begin
            imputed = Mycelia.impute_distances(dist; method = Mycelia.IMPUTE_MEAN)
            expected_mean = Statistics.mean([0.5, 0.8])
            Test.@test imputed[1, 3] ≈ expected_mean
        end

        Test.@testset "Diagonal remains zero" begin
            imputed = Mycelia.impute_distances(dist; method = Mycelia.IMPUTE_MAX)
            Test.@test all(LinearAlgebra.diag(imputed) .== 0.0)
        end
    end

    Test.@testset "Gower Distance" begin
        Test.@testset "Basic numeric features" begin
            ## 3 samples × 2 features (samples as rows)
            matrix = [0.0 0.0; 0.5 0.5; 1.0 1.0]

            gower_dist = Mycelia.gower_distance(matrix)

            Test.@test size(gower_dist) == (3, 3)
            Test.@test LinearAlgebra.issymmetric(gower_dist)
            Test.@test all(LinearAlgebra.diag(gower_dist) .== 0.0)
            ## Sample 1 and Sample 3 should be maximally distant
            Test.@test gower_dist[1, 3] ≈ 1.0
            ## Sample 1 and Sample 2 should be half-distant
            Test.@test gower_dist[1, 2] ≈ 0.5
        end

        Test.@testset "With missing values" begin
            matrix = [1.0 NaN; 2.0 3.0; NaN 4.0]

            gower_dist = Mycelia.gower_distance(matrix; min_shared_features = 1)

            ## Samples 1 and 2 share feature 1
            Test.@test !isnan(gower_dist[1, 2])
            ## Samples 1 and 3 share no features
            Test.@test isnan(gower_dist[1, 3])
            ## Samples 2 and 3 share feature 2
            Test.@test !isnan(gower_dist[2, 3])
        end

        Test.@testset "Mixed feature types" begin
            ## 3 samples with: numeric, binary, categorical
            matrix = Float64[1.0 0.0 1.0; 5.0 1.0 1.0; 3.0 0.0 2.0]

            gower_dist = Mycelia.gower_distance(
                matrix;
                feature_types = [:numeric, :binary, :categorical]
            )

            Test.@test size(gower_dist) == (3, 3)
            Test.@test all(filter(!isnan, vec(gower_dist)) .>= 0.0)
            Test.@test all(filter(!isnan, vec(gower_dist)) .<= 1.0)
        end

        Test.@testset "From RelationalMatrix" begin
            rm = Mycelia.RelationalMatrix{Float64}(
                [1.0 2.0; 3.0 4.0; 5.0 6.0],
                ["A1", "A2", "A3"],
                ["B1", "B2"],
                "entity_a", "entity_b", "value",
                nothing, nothing
            )

            gower_dist = Mycelia.gower_distance(rm)

            Test.@test size(gower_dist) == (3, 3)
        end
    end

    Test.@testset "Rank Normalization" begin
        dist = [0.0 0.3 0.7; 0.3 0.0 0.5; 0.7 0.5 0.0]

        ranked = Mycelia.normalize_ranks(dist)

        Test.@test size(ranked) == size(dist)
        Test.@test all(LinearAlgebra.diag(ranked) .== 0.0)
        ## All values should be in [0, 1]
        Test.@test all(filter(!isnan, vec(ranked)) .>= 0.0)
        Test.@test all(filter(!isnan, vec(ranked)) .<= 1.0)
        ## Smallest distance should have smallest rank
        Test.@test ranked[1, 2] < ranked[1, 3]  ## 0.3 < 0.7
    end

    Test.@testset "Distance Matrix Validation" begin
        Test.@testset "Valid symmetric matrix" begin
            valid_dist = [0.0 0.5 0.8; 0.5 0.0 0.3; 0.8 0.3 0.0]
            validation = Mycelia.validate_distance_matrix(valid_dist)
            Test.@test validation.is_valid
            Test.@test isempty(validation.violations)
            Test.@test validation.n_samples == 3
        end

        Test.@testset "Non-symmetric matrix" begin
            invalid_dist = [0.0 0.5 0.8; 0.4 0.0 0.3; 0.8 0.3 0.0]
            validation = Mycelia.validate_distance_matrix(invalid_dist)
            Test.@test !validation.is_valid
            Test.@test any(occursin.("symmetric", validation.violations))
        end

        Test.@testset "Non-zero diagonal" begin
            invalid_dist = [0.1 0.5; 0.5 0.0]
            validation = Mycelia.validate_distance_matrix(invalid_dist)
            Test.@test !validation.is_valid
            Test.@test any(occursin.("Diagonal", validation.violations))
        end

        Test.@testset "Matrix with NaN" begin
            dist_with_nan = [0.0 NaN; NaN 0.0]
            validation = Mycelia.validate_distance_matrix(dist_with_nan)
            Test.@test validation.has_missing
            Test.@test validation.n_missing_pairs == 1
        end
    end

    Test.@testset "Cluster Ranking" begin
        ## Use seed for reproducibility
        rng = StableRNGs.StableRNG(42)

        ## Create clustered distance matrix with clear structure
        n = 9
        assignments = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        entity_ids = ["E$i" for i in 1:n]

        ## Generate distance matrix with clear cluster structure
        dist = zeros(n, n)
        for i in 1:n, j in 1:n

            if i != j
                if assignments[i] == assignments[j]
                    dist[i, j] = 0.1 + 0.1 * rand(rng)  ## Low intra-cluster
                else
                    dist[i, j] = 0.7 + 0.2 * rand(rng)  ## High inter-cluster
                end
            end
        end
        ## Make symmetric
        dist = (dist + transpose(dist)) / 2

        rankings = Mycelia.rank_cluster_members(dist, assignments, entity_ids)

        Test.@test length(rankings) == 3  ## 3 clusters
        Test.@test all(r -> length(r.entity_indices) == 3, rankings)

        ## Each cluster should have exactly one medoid
        for r in rankings
            medoid_count = count(idx -> idx == r.medoid_index, r.entity_indices)
            Test.@test medoid_count == 1
        end

        ## Medoid should have rank 1 and distance_to_medoid of 0
        for r in rankings
            medoid_local_idx = findfirst(==(r.medoid_index), r.entity_indices)
            Test.@test r.rankings[medoid_local_idx] == 1
            Test.@test r.intra_cluster_distances[medoid_local_idx] == 0.0
        end
    end

    Test.@testset "Rankings to DataFrame" begin
        ## Simple test case
        D = [0.0 0.2 0.8; 0.2 0.0 0.6; 0.8 0.6 0.0]
        assignments = [1, 1, 2]
        entity_ids = ["E1", "E2", "E3"]

        rankings = Mycelia.rank_cluster_members(D, assignments, entity_ids)
        df = Mycelia.rankings_to_dataframe(rankings)

        Test.@test DataFrames.nrow(df) == 3
        Test.@test :cluster_id in DataFrames.propertynames(df)
        Test.@test :entity_id in DataFrames.propertynames(df)
        Test.@test :rank in DataFrames.propertynames(df)
        Test.@test :is_medoid in DataFrames.propertynames(df)
        Test.@test :is_backup in DataFrames.propertynames(df)
        Test.@test :distance_to_medoid in DataFrames.propertynames(df)

        ## Check that medoids are properly marked
        medoid_rows = DataFrames.filter(row -> row.is_medoid, df)
        Test.@test DataFrames.nrow(medoid_rows) == 2  ## One per cluster
    end

    Test.@testset "Cluster Summary" begin
        D = [0.0 0.1 0.5 0.6; 0.1 0.0 0.4 0.5; 0.5 0.4 0.0 0.2; 0.6 0.5 0.2 0.0]
        assignments = [1, 1, 2, 2]
        entity_ids = ["E1", "E2", "E3", "E4"]

        rankings = Mycelia.rank_cluster_members(D, assignments, entity_ids)
        summary = Mycelia.cluster_summary(rankings, D)

        Test.@test DataFrames.nrow(summary) == 2  ## 2 clusters
        Test.@test :cluster_id in DataFrames.propertynames(summary)
        Test.@test :n_members in DataFrames.propertynames(summary)
        Test.@test :medoid_id in DataFrames.propertynames(summary)
        Test.@test :diameter in DataFrames.propertynames(summary)
        Test.@test :avg_intra_distance in DataFrames.propertynames(summary)

        ## Each cluster has 2 members
        Test.@test all(summary.n_members .== 2)
    end

    Test.@testset "Full Pipeline Integration" begin
        ## Create synthetic relational data with clear cluster structure
        rng = StableRNGs.StableRNG(123)

        n_a = 12  ## 12 samples to cluster
        n_b = 8   ## 8 features

        ## Generate data with underlying cluster structure
        ## 3 clusters of 4 samples each
        matrix = zeros(n_a, n_b)
        for i in 1:n_a
            cluster = div(i - 1, 4) + 1  ## Clusters: 1,1,1,1,2,2,2,2,3,3,3,3
            for j in 1:n_b
                if cluster == 1
                    matrix[i, j] = 1.0 + 0.2 * randn(rng)
                elseif cluster == 2
                    matrix[i, j] = 5.0 + 0.2 * randn(rng)
                else
                    matrix[i, j] = 9.0 + 0.2 * randn(rng)
                end
            end
        end

        entity_a_ids = ["A$i" for i in 1:n_a]
        entity_b_ids = ["B$j" for j in 1:n_b]

        rm = Mycelia.RelationalMatrix{Float64}(
            matrix,
            entity_a_ids,
            entity_b_ids,
            "sample", "feature", "score",
            nothing, nothing
        )

        ## Run pipeline
        result = Mycelia.relational_clustering_pipeline(
            rm;
            cluster_by = :entity_a,
            distance_metric = :euclidean,
            imputation = :max_observed,
            ks = 2:5,
            plot_backend = :cairomakie
        )

        ## Check that result has expected fields
        Test.@test haskey(result, :optimal_k)
        Test.@test haskey(result, :assignments)
        Test.@test haskey(result, :rankings)
        Test.@test haskey(result, :rankings_df)
        Test.@test haskey(result, :summary_df)
        Test.@test haskey(result, :distance_matrix)
        Test.@test haskey(result, :hcl)
        Test.@test haskey(result, :figure)

        ## Check that optimal_k is reasonable
        Test.@test result.optimal_k >= 2
        Test.@test result.optimal_k <= 5

        ## Check that all samples are assigned
        Test.@test length(result.assignments) == n_a

        ## Check DataFrame dimensions
        Test.@test DataFrames.nrow(result.rankings_df) == n_a
        Test.@test DataFrames.nrow(result.summary_df) == result.optimal_k
    end

    Test.@testset "Pipeline with Different Metrics" begin
        ## Create simple test data
        rng = StableRNGs.StableRNG(456)

        n_a = 6
        n_b = 4
        matrix = rand(rng, n_a, n_b)

        rm = Mycelia.RelationalMatrix{Float64}(
            matrix,
            ["A$i" for i in 1:n_a],
            ["B$j" for j in 1:n_b],
            "sample", "feature", "value",
            nothing, nothing
        )

        ## Test with cosine distance
        result_cosine = Mycelia.relational_clustering_pipeline(
            rm;
            distance_metric = :cosine,
            ks = 2:3,
            plot_backend = :cairomakie
        )
        Test.@test result_cosine.optimal_k >= 2

        ## Test with Jaccard distance
        result_jaccard = Mycelia.relational_clustering_pipeline(
            rm;
            distance_metric = :jaccard,
            ks = 2:3,
            plot_backend = :cairomakie
        )
        Test.@test result_jaccard.optimal_k >= 2

        ## Test with Gower distance
        result_gower = Mycelia.relational_clustering_pipeline(
            rm;
            distance_metric = :gower,
            ks = 2:3,
            plot_backend = :cairomakie
        )
        Test.@test result_gower.optimal_k >= 2
    end

    Test.@testset "Pipeline Clustering by entity_b" begin
        ## Create test data
        rng = StableRNGs.StableRNG(789)

        n_a = 8
        n_b = 6
        matrix = rand(rng, n_a, n_b)

        rm = Mycelia.RelationalMatrix{Float64}(
            matrix,
            ["A$i" for i in 1:n_a],
            ["B$j" for j in 1:n_b],
            "sample", "feature", "value",
            nothing, nothing
        )

        ## Cluster by entity_b (features) instead of entity_a (samples)
        result = Mycelia.relational_clustering_pipeline(
            rm;
            cluster_by = :entity_b,
            ks = 2:3,
            plot_backend = :cairomakie
        )

        ## Should cluster the 6 features
        Test.@test length(result.assignments) == n_b
        Test.@test DataFrames.nrow(result.rankings_df) == n_b
    end

    Test.@testset "Edge Cases" begin
        Test.@testset "Single entity cluster" begin
            D = [0.0 0.9 0.9; 0.9 0.0 0.1; 0.9 0.1 0.0]
            assignments = [1, 2, 2]  ## Cluster 1 has only one member
            entity_ids = ["E1", "E2", "E3"]

            rankings = Mycelia.rank_cluster_members(D, assignments, entity_ids)

            ## Cluster 1 with single member
            c1 = filter(r -> r.cluster_id == 1, rankings)[1]
            Test.@test length(c1.entity_indices) == 1
            Test.@test c1.medoid_index == 1
            Test.@test isnothing(c1.backup_index)  ## No backup for single member
        end

        Test.@testset "Empty string entity IDs" begin
            df = DataFrames.DataFrame(
                a = ["", "A"],
                b = ["B1", "B1"],
                v = [1.0, 2.0]
            )

            ## Should handle empty strings gracefully
            rm = Mycelia.long_to_relational_matrix(df, :a, :b, :v)
            Test.@test "" in Mycelia.entity_a_ids(rm)
        end

        Test.@testset "All same values" begin
            ## All distances are the same (edge case for ranking)
            D = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
            assignments = [1, 1, 1]
            entity_ids = ["E1", "E2", "E3"]

            ## Should not error
            rankings = Mycelia.rank_cluster_members(D, assignments, entity_ids)
            Test.@test length(rankings) == 1
            Test.@test length(rankings[1].entity_indices) == 3
        end
    end
end
