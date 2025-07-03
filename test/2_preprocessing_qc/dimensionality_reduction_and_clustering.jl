import Random
import Distributions
import MultivariateStats
import UMAP
import Mycelia
import Statistics
import LinearAlgebra
import Test
import Plots
import Clustering

# Set a random seed for reproducibility
Random.seed!(42)

# Function to generate binary (Bernoulli) matrices
function generate_binary_matrix(n_features::Int, n_samples::Int, p::Float64)
    return rand(Distributions.Bernoulli(p), n_features, n_samples)
end

# Function to generate Poisson matrices
function generate_poisson_matrix(n_features::Int, n_samples::Int, λ::Float64)
    return rand(Distributions.Poisson(λ), n_features, n_samples)
end

# Function to compute Jaccard distance for binary matrices
function jaccard_distance(M::AbstractMatrix{<:Integer})
    n_samples = size(M, 2)
    D = zeros(Float64, n_samples, n_samples)
    for i in 1:n_samples
        for j in i+1:n_samples
            intersection = sum(M[:, i] .& M[:, j])
            union = sum(M[:, i] .| M[:, j])
            D[i, j] = 1 - intersection / union
            D[j, i] = D[i, j]
        end
    end
    return D
end

# Function to compute Bray-Curtis distance for count matrices
function bray_curtis_distance(M::AbstractMatrix{<:Integer})
    n_samples = size(M, 2)
    D = zeros(Float64, n_samples, n_samples)
    for i in 1:n_samples
        for j in i+1:n_samples
            sum_abs_diff = sum(abs.(M[:, i] - M[:, j]))
            sum_total = sum(M[:, i]) + sum(M[:, j])
            D[i, j] = sum_abs_diff / sum_total
            D[j, i] = D[i, j]
        end
    end
    return D
end

# Function to plot embeddings with k-means clustering
function plot_embeddings(embeddings; title="", xlabel="", ylabel="", true_labels=nothing, fit_labels=nothing)
    scatter(embeddings[1, :], embeddings[2, :],
           title=title,
           xlabel=xlabel,
           ylabel=ylabel,
           label="",
           legend=:topright)

    if true_labels !== nothing
        for i in unique(true_labels)
            idx = findall(x -> x == i, true_labels)
            scatter!(embeddings[1, idx], embeddings[2, idx],
                     label="True Cluster $i",
                     markershape=:star5,
                     legend=:topright)
        end
    end

    if fit_labels !== nothing
        for i in unique(fit_labels)
            idx = findall(x -> x == i, fit_labels)
            scatter!(embeddings[1, idx], embeddings[2, idx],
                     label="Fit Cluster $i",
                     legend=:bottomright)
        end
    end

    display(plot!)
end

# Test parameters
n_features = 100
n_samples = 50
p_values = [0.1, 0.3, 0.5]
λ_values = [1.0, 3.0, 5.0]

# Generate test matrices
binary_matrices = [generate_binary_matrix(n_features, n_samples, p) for p in p_values]
poisson_matrices = [generate_poisson_matrix(n_features, n_samples, λ) for λ in λ_values]

# Combine all samples into one matrix for clustering
combined_binary = hcat(binary_matrices...)
combined_poisson = hcat(poisson_matrices...)

# Create true labels for clustering
true_labels_binary = repeat(1:length(p_values), inner=n_samples)
true_labels_poisson = repeat(1:length(λ_values), inner=n_samples)
# Test logistic_pca_epca with binary matrices
Test.@testset "logistic_pca_epca" begin
    for (i, M) in enumerate(binary_matrices)
        result = Mycelia.logistic_pca_epca(M, k=5)
        Test.@test size(result.scores) == (5, n_samples)
        Test.@test size(result.loadings) == (5, n_features)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result.scores', 3).assignments

        # Plot embeddings
        plot_embeddings(result.scores;
                       title="Logistic PCA-EPCA - Binary Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=true_labels_binary[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
    end
end

# Test glm_pca_epca with poisson matrices
Test.@testset "glm_pca_epca" begin
    for (i, M) in enumerate(poisson_matrices)
        result = Mycelia.glm_pca_epca(M, k=5)
        Test.@test size(result.scores) == (5, n_samples)
        Test.@test size(result.loadings) == (5, n_features)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result.scores', 3).assignments

        # Plot embeddings
        plot_embeddings(result.scores;
                       title="GLM PCA-EPCA - Poisson Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
    end
end

# Test pca_transform with both types of matrices
Test.@testset "pca_transform" begin
    for (i, M) in enumerate([binary_matrices[1]; poisson_matrices[1]])
        # Test with k specified
        result_k = Mycelia.pca_transform(M, k=5)
        Test.@test size(result_k.scores) == (5, n_samples)
        Test.@test size(result_k.loadings) == (5, n_features)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result_k.scores', 3).assignments

        # Plot embeddings
        plot_embeddings(result_k.scores;
                       title="PCA Transform (k=5) - Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=i <= length(p_values) ? true_labels_binary[1:n_samples] : true_labels_poisson[1:n_samples],
                       fit_labels=fit_labels)

        # Test with var_prop specified
        result_var = Mycelia.pca_transform(M, var_prop=0.95)
        Test.@test size(result_var.scores, 1) <= n_features
        Test.@test size(result_var.loadings, 1) == size(result_var.scores, 1)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result_var.scores', 3).assignments

        # Plot embeddings
        plot_embeddings(result_var.scores;
                       title="PCA Transform (var_prop=0.95) - Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=i <= length(p_values) ? true_labels_binary[1:n_samples] : true_labels_poisson[1:n_samples],
                       fit_labels=fit_labels)
    end
end

# Test negbin_pca_epca with poisson matrices
Test.@testset "negbin_pca_epca" begin
    for (i, M) in enumerate(poisson_matrices)
        result = Mycelia.negbin_pca_epca(M, k=5, r=2)
        Test.@test size(result.scores) == (5, n_samples)
        Test.@test size(result.loadings) == (5, n_features)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result.scores', 3).assignments

        # Plot embeddings
        plot_embeddings(result.scores;
                       title="Negative Binomial PCA-EPCA - Poisson Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
    end
end

# Test distance metrics and PCoA
Test.@testset "distance_metrics_and_pcoa" begin
    for (i, M) in enumerate(binary_matrices)
        D = jaccard_distance(M)
        result = Mycelia.pcoa_from_dist(D, maxoutdim=2)
        Test.@test size(result.coordinates) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result.coordinates', 3).assignments

        # Plot embeddings
        plot_embeddings(result.coordinates;
                       title="PCoA - Jaccard Distance - Binary Matrix $i",
                       xlabel="Coordinate 1",
                       ylabel="Coordinate 2",
                       true_labels=true_labels_binary[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
    end

    for (i, M) in enumerate(poisson_matrices)
        D = bray_curtis_distance(M)
        result = Mycelia.pcoa_from_dist(D, maxoutdim=2)
        Test.@test size(result.coordinates) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result.coordinates', 3).assignments

        # Plot embeddings
        plot_embeddings(result.coordinates;
                       title="PCoA - Bray-Curtis Distance - Poisson Matrix $i",
                       xlabel="Coordinate 1",
                       ylabel="Coordinate 2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
    end
end

# Test UMAP on transformation outputs
Test.@testset "umap_embedding" begin
    for (i, M) in enumerate(binary_matrices)
        # Test UMAP on logistic_pca_epca output
        pca_result = Mycelia.logistic_pca_epca(M, k=5)
        umap_result, umap_model = Mycelia.umap_embed(pca_result.scores)
        Test.@test size(umap_result) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(umap_result', 3).assignments

        # Plot embeddings
        plot_embeddings(umap_result;
                       title="UMAP - Logistic PCA-EPCA - Binary Matrix $i",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=true_labels_binary[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
    end

    for (i, M) in enumerate(poisson_matrices)
        # Test UMAP on glm_pca_epca output
        pca_result = Mycelia.glm_pca_epca(M, k=5)
        umap_result, umap_model = Mycelia.umap_embed(pca_result.scores)
        Test.@test size(umap_result) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(umap_result', 3).assignments

        # Plot embeddings
        plot_embeddings(umap_result;
                       title="UMAP - GLM PCA-EPCA - Poisson Matrix $i",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)

        # Test UMAP on negbin_pca_epca output
        pca_result = Mycelia.negbin_pca_epca(M, k=5, r=2)
        umap_result, umap_model = Mycelia.umap_embed(pca_result.scores)
        Test.@test size(umap_result) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(umap_result', 3).assignments

        # Plot embeddings
        plot_embeddings(umap_result;
                       title="UMAP - Negative Binomial PCA-EPCA - Poisson Matrix $i",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)

        # Test UMAP on pca_transform output
        pca_result = Mycelia.pca_transform(M, k=5)
        umap_result, umap_model = Mycelia.umap_embed(pca_result.scores)
        Test.@test size(umap_result) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(umap_result', 3).assignments

        # Plot embeddings
        plot_embeddings(umap_result;
                       title="UMAP - PCA Transform (k=5) - Poisson Matrix $i",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
    end
end