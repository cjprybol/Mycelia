# # Tutorial 19: Relational Clustering
#
# This tutorial demonstrates how to cluster entities based on their relationships
# with other entities. This is useful for phenotypic clustering where entities are
# characterized by their interaction profiles with a set of features or targets.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - How to construct a RelationalMatrix from tabular data
# - Distance metrics for relational data with missing values (Gower, Euclidean, etc.)
# - How to impute missing pairwise distances
# - How to identify optimal cluster counts using silhouette analysis
# - How to rank entities within clusters by centrality (medoid identification)
# - How to interpret clustering results and generate prioritized entity lists

# ## Setup
#
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/19_relational_clustering.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import DataFrames
import Random
import Statistics
import Distances

Random.seed!(42)

println("=" ^ 60)
println("Tutorial 19: Relational Clustering")
println("=" ^ 60)

# ## Part 1: Understanding Relational Data
#
# Relational data captures measurements between two types of entities.
# This is a common pattern in many biological and scientific domains:
#
# - Interaction assays: Entity A tested against Entity B
# - Host range profiles: Samples characterized by their responses to a panel
# - Drug susceptibility: Compounds tested against pathogen isolates
# - Gene expression: Genes measured across experimental conditions
#
# The key insight is that we can cluster Entity A based on their profiles
# across Entity B (or vice versa).

println("\n" * "=" ^ 60)
println("Part 1: Understanding Relational Data")
println("=" ^ 60)

# ### Create Example Data
#
# We'll simulate an interaction dataset for demonstration.
# This could represent any entity-entity measurement scenario.

## Simulate 20 samples tested against 15 features
n_samples = 20
n_features = 15

## Create long-format data (typical experimental output)
long_data = DataFrames.DataFrame(
    sample_id = String[],
    feature_id = String[],
    measurement = Float64[]
)

Random.seed!(42)

## Generate data with underlying cluster structure
## Samples 1-7: Profile Type A (high response on features 1-5)
## Samples 8-14: Profile Type B (high response on features 6-10)
## Samples 15-20: Profile Type C (high response on features 11-15)

for s in 1:n_samples
    for f in 1:n_features
        ## Simulate some missing data (15% not tested)
        if rand() < 0.15
            continue
        end

        ## Generate measurement based on cluster membership
        if s <= 7
            ## Profile Type A: high on features 1-5
            measurement = f <= 5 ? (0.7 + 0.3 * rand()) : (0.0 + 0.15 * rand())
        elseif s <= 14
            ## Profile Type B: high on features 6-10
            measurement = 6 <= f <= 10 ? (0.6 + 0.4 * rand()) : (0.0 + 0.1 * rand())
        else
            ## Profile Type C: high on features 11-15
            measurement = f >= 11 ? (0.5 + 0.5 * rand()) : (0.05 + 0.1 * rand())
        end

        push!(long_data, (
            sample_id = "Sample_$(lpad(s, 2, '0'))",
            feature_id = "Feature_$(lpad(f, 2, '0'))",
            measurement = round(measurement, digits=3)
        ))
    end
end

println("\nGenerated $(DataFrames.nrow(long_data)) measurements")
println("This simulates $n_samples samples tested against $n_features features")
println("\nPreview of long-format data:")
println(first(long_data, 10))

# ## Part 2: Constructing a RelationalMatrix
#
# Convert long-format data to a matrix representation using `long_to_relational_matrix`.
# This function handles:
# - Automatic detection of unique entities
# - Aggregation of duplicate measurements (median by default)
# - Missing value representation (NaN by default)

println("\n" * "=" ^ 60)
println("Part 2: RelationalMatrix Construction")
println("=" ^ 60)

rm = Mycelia.long_to_relational_matrix(
    long_data,
    :sample_id,      ## Column containing row entity IDs
    :feature_id,     ## Column containing column entity IDs
    :measurement;    ## Column containing values
    entity_a_name = "sample",
    entity_b_name = "feature",
    measurement_name = "response"
)

## Get summary statistics
stats = Mycelia.summarize_relational_matrix(rm)

println("\nMatrix coverage: $(stats.pct_filled)% of cells have values")
println("Missing $(stats.n_missing_values) measurements out of $(stats.n_total_cells) total")

# ## Part 3: Computing Pairwise Distances
#
# To cluster samples, we need to compute pairwise distances between their profiles.
# Mycelia supports several distance metrics:
#
# - **Euclidean**: Standard L2 distance, sensitive to magnitude
# - **Cosine**: Angular distance, invariant to magnitude
# - **Jaccard**: For binary presence/absence data
# - **Gower**: Handles missing values natively, works with mixed types
#
# Gower distance is particularly useful for relational data because it:
# 1. Computes partial distances using only shared (non-missing) features
# 2. Normalizes each feature by its range
# 3. Returns NaN when insufficient features are shared

println("\n" * "=" ^ 60)
println("Part 3: Pairwise Distance Computation")
println("=" ^ 60)

## Extract the matrix (samples as rows, features as columns)
profile_matrix = Mycelia.Matrix(rm)
println("\nProfile matrix shape: $(size(profile_matrix))")

## Compute Gower distance (handles NaN automatically)
gower_dist = Mycelia.gower_distance(profile_matrix; min_shared_features=3)

## Check how many pairs have valid distances
validation = Mycelia.validate_distance_matrix(gower_dist)
println("\nDistance matrix validation:")
println("  Valid: $(validation.is_valid)")
println("  Missing pairs: $(validation.n_missing_pairs) out of $(n_samples * (n_samples-1) ÷ 2)")

# ## Part 4: Handling Missing Distances
#
# When some sample pairs don't have enough shared features, their distance is NaN.
# Before clustering, we need to impute these missing values.
#
# Available imputation methods:
# - `IMPUTE_MAX`: Use theoretical maximum (1.0 for normalized distances)
# - `IMPUTE_MAX_OBSERVED`: Use the maximum observed distance
# - `IMPUTE_MEDIAN`: Use the median of observed distances
# - `IMPUTE_MEAN`: Use the mean of observed distances
#
# `IMPUTE_MAX_OBSERVED` is often a good default: it treats pairs with
# insufficient data as "maximally different within the observed range".

println("\n" * "=" ^ 60)
println("Part 4: Distance Imputation")
println("=" ^ 60)

if validation.has_missing
    println("\nImputing $(validation.n_missing_pairs) missing distance pairs...")

    ## Impute using maximum observed distance
    imputed_dist = Mycelia.impute_distances(gower_dist; method=Mycelia.IMPUTE_MAX_OBSERVED)

    ## Verify imputation worked
    validation_after = Mycelia.validate_distance_matrix(imputed_dist)
    println("After imputation: $(validation_after.n_missing_pairs) missing pairs")
else
    imputed_dist = gower_dist
    println("No missing distances to impute")
end

## Normalize to [0, 1] range
normalized_dist = Mycelia.normalize_distance_matrix(imputed_dist)

# ## Part 5: Identifying Optimal Cluster Count
#
# Use silhouette analysis to find the optimal number of clusters.
# The silhouette score measures how similar each sample is to its own cluster
# compared to other clusters. Higher scores indicate better-defined clusters.

println("\n" * "=" ^ 60)
println("Part 5: Optimal Cluster Identification")
println("=" ^ 60)

cluster_result = Mycelia.identify_optimal_number_of_clusters(
    normalized_dist;
    ks = 2:6,  ## Test 2 to 6 clusters
    plot_backend = :cairomakie
)

println("\nSilhouette scores by k:")
for (k, score) in zip(cluster_result.ks, cluster_result.silhouette_scores)
    marker = k == cluster_result.optimal_number_of_clusters ? " <-- optimal" : ""
    println("  k=$k: $(round(score, digits=3))$marker")
end

optimal_k = cluster_result.optimal_number_of_clusters
assignments = cluster_result.assignments

println("\nOptimal number of clusters: $optimal_k")

## Show cluster sizes
cluster_sizes = [count(==(i), assignments) for i in 1:optimal_k]
println("Cluster sizes: $cluster_sizes")

# ## Part 6: Ranking Cluster Members
#
# Within each cluster, identify the most representative sample (medoid)
# and rank all members by their distance to the medoid.
#
# The **medoid** is the sample that minimizes the sum of distances to all
# other cluster members. This is more robust than the centroid for
# non-Euclidean distances.
#
# Rankings help prioritize samples: rank 1 is the best representative,
# rank 2 is a good backup, etc.

println("\n" * "=" ^ 60)
println("Part 6: Cluster Member Ranking")
println("=" ^ 60)

entity_ids = Mycelia.entity_a_ids(rm)

rankings = Mycelia.rank_cluster_members(
    normalized_dist,
    assignments,
    entity_ids;
    metric_type = :distance
)

## Convert to DataFrame for easy viewing
rankings_df = Mycelia.rankings_to_dataframe(rankings)

println("\nFull rankings table:")
println(rankings_df)

## Show medoids (cluster representatives)
println("\n" * "-" ^ 40)
println("Cluster Representatives (Medoids):")
println("-" ^ 40)
medoids_df = DataFrames.filter(row -> row.is_medoid, rankings_df)
for row in DataFrames.eachrow(medoids_df)
    println("  Cluster $(row.cluster_id): $(row.entity_id)")
end

## Show backups
println("\n" * "-" ^ 40)
println("Backup Representatives:")
println("-" ^ 40)
backups_df = DataFrames.filter(row -> row.is_backup, rankings_df)
for row in DataFrames.eachrow(backups_df)
    println("  Cluster $(row.cluster_id): $(row.entity_id) (dist to medoid: $(round(row.distance_to_medoid, digits=3)))")
end

# ## Part 7: Cluster Summary Statistics
#
# Generate summary statistics for each cluster including:
# - Number of members
# - Cluster diameter (maximum intra-cluster distance)
# - Average intra-cluster distance (cohesion)

println("\n" * "=" ^ 60)
println("Part 7: Cluster Summary Statistics")
println("=" ^ 60)

summary_df = Mycelia.cluster_summary(rankings, normalized_dist)

println("\nCluster Summary:")
println(summary_df)

println("\nInterpretation:")
for row in DataFrames.eachrow(summary_df)
    println("  Cluster $(row.cluster_id):")
    println("    - $(row.n_members) members, medoid: $(row.medoid_id)")
    println("    - Diameter: $(round(row.diameter, digits=3)) (max spread)")
    println("    - Avg cohesion: $(round(row.avg_intra_distance, digits=3))")
end

# ## Part 8: Using the Complete Pipeline
#
# For convenience, all steps can be run with a single function call:
# `relational_clustering_pipeline`
#
# This function:
# 1. Extracts the profile matrix from the RelationalMatrix
# 2. Computes pairwise distances with the specified metric
# 3. Imputes missing distances
# 4. Identifies optimal cluster count via silhouette analysis
# 5. Ranks cluster members by centrality
# 6. Returns all results in a named tuple

println("\n" * "=" ^ 60)
println("Part 8: Full Pipeline")
println("=" ^ 60)

## Run the complete pipeline
## Note: Using Gower distance because our data has missing values.
## Gower handles NaN natively by computing partial distances.
## Euclidean/Cosine would produce NaN for any pair with shared missing features.
pipeline_result = Mycelia.relational_clustering_pipeline(
    rm;
    cluster_by = :entity_a,       ## Cluster samples (rows)
    distance_metric = :gower,     ## Best for data with missing values
    imputation = :max_observed,
    ks = 2:6,
    plot_backend = :cairomakie
)

println("\nPipeline Results Summary:")
println("  Optimal k: $(pipeline_result.optimal_k)")
println("  Rankings DataFrame: $(DataFrames.nrow(pipeline_result.rankings_df)) rows")
println("  Summary DataFrame: $(DataFrames.nrow(pipeline_result.summary_df)) clusters")

# ## Part 9: Clustering the Other Dimension
#
# We can also cluster features (entity_b) based on their sample profiles.
# This is useful when you want to identify groups of features that behave
# similarly across samples.

println("\n" * "=" ^ 60)
println("Part 9: Clustering Features")
println("=" ^ 60)

feature_result = Mycelia.relational_clustering_pipeline(
    rm;
    cluster_by = :entity_b,  ## Cluster features (columns)
    distance_metric = :gower, ## Handles missing values
    ks = 2:5,
    plot_backend = :cairomakie
)

println("\nFeature Clustering Results:")
println("  Optimal k: $(feature_result.optimal_k)")
println(feature_result.summary_df)

# ## Part 10: Distance Metric Comparison Notes
#
# Different distance metrics are appropriate for different data types:
#
# - **Gower**: Best for data with missing values (handles NaN natively)
# - **Euclidean**: Standard distance, but produces NaN if any feature is missing
# - **Cosine**: Angular similarity, also produces NaN with missing values
# - **Jaccard**: For binary presence/absence data
#
# Since our data has ~15% missing values, Gower is the appropriate choice.
# For complete data without missing values, all metrics would work.

println("\n" * "=" ^ 60)
println("Part 10: Distance Metric Notes")
println("=" ^ 60)

println("\nMetric Selection Guide:")
println("  - Gower: Best for data with missing values (used in this tutorial)")
println("  - Euclidean: Standard L2 distance, requires complete data")
println("  - Cosine: Angular similarity, requires complete data")
println("  - Jaccard: For binary presence/absence data")

println("\nWith ~$(round(100 - 100*DataFrames.nrow(long_data)/(n_samples*n_features), digits=1))% missing data,")
println("Gower distance is the recommended choice for this dataset.")

# ## Summary
#
# In this tutorial, we learned how to:
#
# 1. **Construct a RelationalMatrix** from long-format experimental data
# 2. **Compute pairwise distances** between entity profiles using various metrics
# 3. **Handle missing values** with appropriate imputation strategies
# 4. **Identify optimal cluster count** using silhouette analysis
# 5. **Rank cluster members** to identify representatives (medoids) and backups
# 6. **Generate summary statistics** for cluster interpretation
# 7. **Use the complete pipeline** for end-to-end analysis
#
# Key takeaways:
# - Gower distance is robust to missing values
# - Medoids are more robust than centroids for non-Euclidean distances
# - Silhouette analysis helps identify natural cluster structure
# - Rankings help prioritize entities for downstream applications

println("\n" * "=" ^ 60)
println("Tutorial Complete!")
println("=" ^ 60)
println("\nKey functions demonstrated:")
println("  - Mycelia.long_to_relational_matrix()")
println("  - Mycelia.gower_distance()")
println("  - Mycelia.impute_distances()")
println("  - Mycelia.identify_optimal_number_of_clusters()")
println("  - Mycelia.rank_cluster_members()")
println("  - Mycelia.rankings_to_dataframe()")
println("  - Mycelia.cluster_summary()")
println("  - Mycelia.relational_clustering_pipeline()")

nothing
