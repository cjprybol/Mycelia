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

# Parameters
N = 7      # Number of distributions
X = 100      # Number of samples per distribution
L = 1000     # Length of each distribution (number of features)

# Step 1: Generate N distributions (each is a vector of probabilities)
binary_probabilities = [rand(L) for _ in 1:N]  # Each element is a vector of length L with values in [0,1]

# Step 2: For each distribution, sample X binary vectors (each column is a sample)
binary_samples = [hcat([rand.(Distributions.Bernoulli.(p)) for _ in 1:X]...) for p in binary_probabilities]

# Concatenate all samples into one matrix (L x (N*X))
all_samples = hcat(binary_samples...)

# Create a label vector: for each distribution, repeat its index X times
all_labels = repeat(1:N, inner=X)

# Shuffle columns and labels together
perm = Random.shuffle(1:size(all_samples, 2))
shuffled_samples = all_samples[:, perm]
shuffled_labels = all_labels[perm]


# Parameters for Poisson λ distribution
λ_scale = 0.7         # Lower = more bias toward 0
λ_max = 128.0          # User-defined maximum (set to Inf for no max)

# Step 1: Generate N distributions (each is a vector of Poisson means, biased toward 0)
poisson_means = [
    clamp.(rand(Distributions.Exponential(λ_scale), L), 0, λ_max)
    for _ in 1:N
]

# Step 2: For each distribution, sample X count vectors (each column is a sample)
poisson_samples = [hcat([rand.(Distributions.Poisson.(λ)) for _ in 1:X]...) for λ in poisson_means]

# Concatenate all samples into one matrix (L x (N*X))
all_poisson_samples = hcat(poisson_samples...)

# Create a label vector: for each distribution, repeat its index X times
all_poisson_labels = repeat(1:N, inner=X)

# Shuffle columns and labels together
perm_poisson = Random.shuffle(1:size(all_poisson_samples, 2))
shuffled_poisson_samples = all_poisson_samples[:, perm_poisson]
shuffled_poisson_labels = all_poisson_labels[perm_poisson]

# Test logistic_pca_epca with binary matrices
Test.@testset "logistic_pca_epca" begin
    for (i, M) in enumerate(binary_matrices)
        result = Mycelia.logistic_pca_epca(M, k=5)
        Test.@test size(result.scores) == (5, n_samples)
        Test.@test size(result.loadings) == (5, n_features)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result.scores', 3).assignments

        # Plot embeddings
        plt = Mycelia.plot_embeddings(result.scores;
                       title="Logistic PCA-EPCA - Binary Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=true_labels_binary[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)
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
        plt = Mycelia.plot_embeddings(result.scores;
                       title="GLM PCA-EPCA - Poisson Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)
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
        plt = Mycelia.plot_embeddings(result_k.scores;
                       title="PCA Transform (k=5) - Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=i <= length(p_values) ? true_labels_binary[1:n_samples] : true_labels_poisson[1:n_samples],
                       fit_labels=fit_labels)
        display(plt)

        # Test with var_prop specified
        result_var = Mycelia.pca_transform(M, var_prop=0.95)
        Test.@test size(result_var.scores, 1) <= n_features
        Test.@test size(result_var.loadings, 1) == size(result_var.scores, 1)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result_var.scores', 3).assignments

        # Plot embeddings
        plt = Mycelia.plot_embeddings(result_var.scores;
                       title="PCA Transform (var_prop=0.95) - Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=i <= length(p_values) ? true_labels_binary[1:n_samples] : true_labels_poisson[1:n_samples],
                       fit_labels=fit_labels)
        display(plt)
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
        plt = Mycelia.plot_embeddings(result.scores;
                       title="Negative Binomial PCA-EPCA - Poisson Matrix $i",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)
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
        plt = Mycelia.plot_embeddings(result.coordinates;
                       title="PCoA - Jaccard Distance - Binary Matrix $i",
                       xlabel="Coordinate 1",
                       ylabel="Coordinate 2",
                       true_labels=true_labels_binary[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)
    end

    for (i, M) in enumerate(poisson_matrices)
        D = bray_curtis_distance(M)
        result = Mycelia.pcoa_from_dist(D, maxoutdim=2)
        Test.@test size(result.coordinates) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(result.coordinates', 3).assignments

        # Plot embeddings
        plt = Mycelia.plot_embeddings(result.coordinates;
                       title="PCoA - Bray-Curtis Distance - Poisson Matrix $i",
                       xlabel="Coordinate 1",
                       ylabel="Coordinate 2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)
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
        plt = Mycelia.plot_embeddings(umap_result;
                       title="UMAP - Logistic PCA-EPCA - Binary Matrix $i",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=true_labels_binary[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)
    end

    for (i, M) in enumerate(poisson_matrices)
        # Test UMAP on glm_pca_epca output
        pca_result = Mycelia.glm_pca_epca(M, k=5)
        umap_result, umap_model = Mycelia.umap_embed(pca_result.scores)
        Test.@test size(umap_result) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(umap_result', 3).assignments

        # Plot embeddings
        plt = Mycelia.plot_embeddings(umap_result;
                       title="UMAP - GLM PCA-EPCA - Poisson Matrix $i",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)

        # Test UMAP on negbin_pca_epca output
        pca_result = Mycelia.negbin_pca_epca(M, k=5, r=2)
        umap_result, umap_model = Mycelia.umap_embed(pca_result.scores)
        Test.@test size(umap_result) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(umap_result', 3).assignments

        # Plot embeddings
        plt = Mycelia.plot_embeddings(umap_result;
                       title="UMAP - Negative Binomial PCA-EPCA - Poisson Matrix $i",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)

        # Test UMAP on pca_transform output
        pca_result = Mycelia.pca_transform(M, k=5)
        umap_result, umap_model = Mycelia.umap_embed(pca_result.scores)
        Test.@test size(umap_result) == (2, n_samples)

        # Fit k-means clustering
        fit_labels = Clustering.kmeans(umap_result', 3).assignments

        # Plot embeddings
        plt = Mycelia.plot_embeddings(umap_result;
                       title="UMAP - PCA Transform (k=5) - Poisson Matrix $i",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=true_labels_poisson[(i-1)*n_samples+1:i*n_samples],
                       fit_labels=fit_labels)
        display(plt)
    end
end

# Example for logistic_pca_epca on shuffled binary data
Test.@testset "logistic_pca_epca_shuffled" begin
    result = Mycelia.logistic_pca_epca(shuffled_binary, k=5)
    Test.@test size(result.scores) == (5, size(shuffled_binary, 2))
    Test.@test size(result.loadings) == (5, n_features)

    # Fit k-means clustering
    fit_labels = Clustering.kmeans(result.scores', 3).assignments

    # Plot embeddings with color for fit_labels and marker for known group
    plt = Mycelia.plot_embeddings(result.scores;
                   title="Logistic PCA-EPCA - Shuffled Binary Data",
                   xlabel="PC1",
                   ylabel="PC2",
                   true_labels=shuffled_labels_binary,
                   fit_labels=fit_labels,
                   legend_title_shape="Known Group",
                   legend_title_color="K-means Cluster")
    display(plt)
end

# Example for glm_pca_epca on shuffled poisson data
Test.@testset "glm_pca_epca_shuffled" begin
    result = Mycelia.glm_pca_epca(shuffled_poisson, k=5)
    Test.@test size(result.scores) == (5, size(shuffled_poisson, 2))
    Test.@test size(result.loadings) == (5, n_features)

    # Fit k-means clustering
    fit_labels = Clustering.kmeans(result.scores', 3).assignments

    # Plot embeddings with color for fit_labels and marker for known group
    plt = Mycelia.plot_embeddings(result.scores;
                   title="GLM PCA-EPCA - Shuffled Poisson Data",
                   xlabel="PC1",
                   ylabel="PC2",
                   true_labels=shuffled_labels_poisson,
                   fit_labels=fit_labels,
                   legend_title_shape="Known Group",
                   legend_title_color="K-means Cluster")
    display(plt)
end