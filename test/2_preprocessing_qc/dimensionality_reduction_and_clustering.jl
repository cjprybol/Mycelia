import Mycelia
import Random
import Distributions
import MultivariateStats
import UMAP
import Statistics
import LinearAlgebra
import Plots
import Clustering
import Distances

@testset "Confusion Matrix and Metrics - Animal Dataset" begin
    true_labels = ["cat", "dog", "cat", "dog", "cat", "bird"]
    pred_labels = ["cat", "cat", "cat", "dog", "dog", "bird"]

    # Expected confusion matrix:
    #           Pred: bird   cat   dog
    # True: bird   [  1      0     0 ]
    #       cat    [  0      2     1 ]
    #       dog    [  0      1     1 ]
    expected_cm = [1 0 0;
                   0 2 1;
                   0 1 1]
    expected_labels = ["bird", "cat", "dog"]

    cm_out = Mycelia.confusion_matrix(true_labels, pred_labels)
    @test cm_out.cm == expected_cm
    @test cm_out.labels == expected_labels

    metrics = Mycelia.precision_recall_f1(true_labels, pred_labels)
    # Per-label metrics
    @test isapprox(metrics.precisions["bird"], 1.0; atol=1e-8)
    @test isapprox(metrics.recalls["bird"], 1.0; atol=1e-8)
    @test isapprox(metrics.f1s["bird"], 1.0; atol=1e-8)

    @test isapprox(metrics.precisions["cat"], 2/3; atol=1e-8)
    @test isapprox(metrics.recalls["cat"], 2/3; atol=1e-8)
    @test isapprox(metrics.f1s["cat"], 2/3; atol=1e-8)

    @test isapprox(metrics.precisions["dog"], 0.5; atol=1e-8)
    @test isapprox(metrics.recalls["dog"], 0.5; atol=1e-8)
    @test isapprox(metrics.f1s["dog"], 0.5; atol=1e-8)

    # Macro-averaged metrics
    @test isapprox(metrics.macro_precision, (1.0 + 2/3 + 0.5)/3; atol=1e-8)
    @test isapprox(metrics.macro_recall, (1.0 + 2/3 + 0.5)/3; atol=1e-8)
    @test isapprox(metrics.macro_f1, (1.0 + 2/3 + 0.5)/3; atol=1e-8)
end

@testset "Confusion Matrix and Metrics - ABC Dataset" begin
    true_labels = ["A", "A", "A", "B", "B", "C", "C", "C"]
    pred_labels = ["A", "B", "C", "B", "C", "C", "A", "B"]

    # Expected confusion matrix:
    #           Pred: A  B  C
    # True: A   [ 1  1  1 ]
    #       B   [ 0  1  1 ]
    #       C   [ 1  1  1 ]
    expected_cm = [1 1 1;
                   0 1 1;
                   1 1 1]
    expected_labels = ["A", "B", "C"]

    cm_out = Mycelia.confusion_matrix(true_labels, pred_labels)
    @test cm_out.cm == expected_cm
    @test cm_out.labels == expected_labels

    metrics = Mycelia.precision_recall_f1(true_labels, pred_labels)
    # Per-label metrics
    @test isapprox(metrics.precisions["A"], 0.5; atol=1e-8)
    @test isapprox(metrics.recalls["A"], 1/3; atol=1e-8)
    @test isapprox(metrics.f1s["A"], 0.4; atol=1e-8)

    @test isapprox(metrics.precisions["B"], 1/3; atol=1e-8)
    @test isapprox(metrics.recalls["B"], 0.5; atol=1e-8)
    @test isapprox(metrics.f1s["B"], 0.4; atol=1e-8)

    @test isapprox(metrics.precisions["C"], 1/3; atol=1e-8)
    @test isapprox(metrics.recalls["C"], 1/3; atol=1e-8)
    @test isapprox(metrics.f1s["C"], 1/3; atol=1e-8)

    # Macro-averaged metrics
    @test isapprox(metrics.macro_precision, (0.5 + 1/3 + 1/3)/3; atol=1e-8)
    @test isapprox(metrics.macro_recall, (1/3 + 0.5 + 1/3)/3; atol=1e-8)
    @test isapprox(metrics.macro_f1, (0.4 + 0.4 + 1/3)/3; atol=1e-8)
end

@testset "Binary Matrix Processing" begin
    # Set a random seed for reproducibility
    Random.seed!(42)

    # Parameters
    n_distributions = 7      # Number of distributions
    n_samples = 10      # Number of samples per distribution
    n_features = 100     # Length of each distribution (number of features)

    # Bernoulli (binary 0/1)
    binary_probabilities = [rand(n_features) for _ in 1:n_distributions]
    binary_samples = [hcat([rand.(Distributions.Bernoulli.(p)) for _ in 1:n_samples]...) for p in binary_probabilities]
    binary_matrix = hcat(binary_samples...)
    binary_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(binary_labels))
    shuffled_binary_matrix = binary_matrix[:, perm]
    shuffled_binary_labels = binary_labels[perm]
    Mycelia.sanity_check_matrix(shuffled_binary_matrix)

    @testset "Sanity Check - Binary Matrix" begin
        summary = Mycelia.sanity_check_matrix(shuffled_binary_matrix)
        @test summary[:n_features] == n_features
        @test summary[:n_samples] == n_samples * n_distributions
        @test summary[:value_type] == Bool || summary[:value_type] == Int || summary[:value_type] == UInt8
        @test summary[:is_binary] == true
        @test summary[:is_integer] == true
        @test summary[:is_nonnegative] == true
        @test summary[:is_strictly_positive] == false
        @test summary[:is_in_01] == false
        @test summary[:is_overdispersed] == false
        @test summary[:suggested_epca] == :bernoulli_pca_epca
        @test summary[:suggested_distance] == :jaccard_distance
    end

    # Store results for ranking
    binary_method_accuracies = []

    # Distance clustering + Optimal Hierarchical Clustering
    @testset "Jaccard Distance + Optimal Hierarchical Clustering" begin
        println("[Binary] Testing: Distance Clustering (Jaccard Distance + Optimal Hierarchical Clustering)")
        binary_distance_matrix = Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix)
        binary_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(binary_distance_matrix)
        @test binary_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, binary_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(binary_method_accuracies, ("Jaccard Distance + Optimal Hierarchical Clustering", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    @testset "Jaccard Distance + KMeans" begin
        println("[Binary] Testing: Jaccard Distance) + KMeans")
        binary_distance_matrix = Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix)
        # Classical MDS/PCoA to get coordinates from distance matrix
        pcoa_result = Mycelia.pcoa_from_dist(binary_distance_matrix)
        kmeans_labels = Clustering.kmeans(pcoa_result.coordinates, n_distributions).assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, kmeans_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binary_method_accuracies, ("Jaccard Distance + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    @testset "Jaccard Distance + KMedoids" begin
        println("[Binary] Testing: Jaccard Distance + KMedoids")
        binary_distance_matrix = Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix)
        kmedoids_result = Clustering.kmedoids(binary_distance_matrix, n_distributions)
        kmedoids_labels = kmedoids_result.assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, kmedoids_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(binary_method_accuracies, ("Jaccard Distance + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Distance + PCoA + KMeans
    @testset "Jaccard Distance + PCoA + KMeans" begin
        println("[Binary] Testing: Jaccard Distance + PCoA + KMeans")
        pcoa_binary_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix))
        pcoa_fit_binary_labels = Clustering.kmeans(pcoa_binary_result.coordinates, n_distributions).assignments
        pcoa_fit_binary_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, pcoa_fit_binary_labels)
        plt = Mycelia.plot_embeddings(pcoa_binary_result.coordinates;
                       title="Jaccard Distance + PCoA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=pcoa_fit_binary_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, pcoa_fit_binary_labels)
        @test evaluation_result.macro_f1 >= .5
        @test evaluation_result.macro_precision >= .5
        @test evaluation_result.macro_recall >= .5
        @test evaluation_result.accuracy >= .5
        push!(binary_method_accuracies, ("Jaccard Distance + PCoA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    @testset "Jaccard Distance + PCoA + KMedoids" begin
        println("[Binary] Testing: Jaccard Distance + PCoA + KMedoids")
        pcoa_binary_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix))
        # Compute distance matrix from PCoA coordinates (e.g., Euclidean)
        embedding = pcoa_binary_result.coordinates
        dist_matrix = Distances.pairwise(Distances.Euclidean(), embedding; dims=2)
        kmedoids_result = Clustering.kmedoids(dist_matrix, n_distributions)
        pcoa_fit_binary_labels = kmedoids_result.assignments
        pcoa_fit_binary_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, pcoa_fit_binary_labels)
        plt = Mycelia.plot_embeddings(embedding;
                       title="Jaccard Distance + PCoA + KMedoids",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=pcoa_fit_binary_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, pcoa_fit_binary_labels)
        @test evaluation_result.macro_f1 >= 0.5
        @test evaluation_result.macro_precision >= 0.5
        @test evaluation_result.macro_recall >= 0.5
        @test evaluation_result.accuracy >= 0.5
        push!(binary_method_accuracies, ("Jaccard Distance + PCoA + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Jaccard Distance + PCoA + UMAP + KMeans
    @testset "Jaccard Distance + PCoA + UMAP + KMeans" begin
        println("[Binary] Testing: Jaccard Distance + PCoA + UMAP + KMeans")
        pcoa_binary_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix))
        pcoa_binary_umap_model = Mycelia.umap_embed(pcoa_binary_result.coordinates)
        @test size(pcoa_binary_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_binary_umap_fit_labels = Clustering.kmeans(pcoa_binary_umap_model.embedding, n_distributions).assignments
        pcoa_binary_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, pcoa_binary_umap_fit_labels)
        plt = Mycelia.plot_embeddings(pcoa_binary_umap_model.embedding;
                       title="Distance + PCoA (Jaccard Distance) + UMAP + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=pcoa_binary_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, pcoa_binary_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binary_method_accuracies, ("Jaccard Distance + PCoA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    @testset "Jaccard Distance + PCoA + UMAP + KMedoids" begin
        println("[Binary] Testing: Jaccard Distance + PCoA + UMAP + KMedoids")
        pcoa_binary_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix))
        pcoa_binary_umap_model = Mycelia.umap_embed(pcoa_binary_result.coordinates)
        @test size(pcoa_binary_umap_model.embedding) == (2, n_samples * n_distributions)
        # Compute distance matrix from UMAP embedding (Euclidean)
        embedding = pcoa_binary_umap_model.embedding
        dist_matrix = Distances.pairwise(Distances.Euclidean(), embedding; dims=2)
        kmedoids_result = Clustering.kmedoids(dist_matrix, n_distributions)
        pcoa_binary_umap_fit_labels = kmedoids_result.assignments
        pcoa_binary_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, pcoa_binary_umap_fit_labels)
        plt = Mycelia.plot_embeddings(embedding;
                       title="Jaccard Distance + PCoA + UMAP + KMedoids",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=pcoa_binary_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, pcoa_binary_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binary_method_accuracies, ("Jaccard Distance + PCoA + UMAP + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # logisticPCA + KMeans
    @testset "logisticPCA-EPCA + KMeans" begin
        println("[Binary] Testing: logisticPCA-EPCA + KMeans")
        logistic_pca_result = Mycelia.logistic_pca_epca(shuffled_binary_matrix, k=5)
        fit_labels = Clustering.kmeans(logistic_pca_result.scores, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, fit_labels)
        plt = Mycelia.plot_embeddings(logistic_pca_result.scores;
                       title="logisticPCA-EPCA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 0.0
        @test evaluation_result.macro_precision >= 0.0
        @test evaluation_result.macro_recall >= 0.0
        @test evaluation_result.accuracy >= 0.0
        push!(binary_method_accuracies, ("logisticPCA-EPCA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # logisticPCA + KMedoids
    @testset "logisticPCA-EPCA + KMedoids" begin
        println("[Binary] Testing: logisticPCA-EPCA + KMedoids")
        logistic_pca_result = Mycelia.logistic_pca_epca(shuffled_binary_matrix, k=5)
        # Compute distance matrix from logistic PCA scores (Euclidean)
        dist_matrix = Distances.pairwise(Distances.Euclidean(), logistic_pca_result.scores; dims=2)
        kmedoids_result = Clustering.kmedoids(dist_matrix, n_distributions)
        fit_labels = kmedoids_result.assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, fit_labels)
        plt = Mycelia.plot_embeddings(logistic_pca_result.scores;
                       title="logisticPCA-EPCA + KMedoids",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 0.0
        @test evaluation_result.macro_precision >= 0.0
        @test evaluation_result.macro_recall >= 0.0
        @test evaluation_result.accuracy >= 0.0
        push!(binary_method_accuracies, ("logisticPCA-EPCA + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # logisticPCA + UMAP + KMeans
    @testset "logisticPCA-EPCA + UMAP + KMeans" begin
        println("[Binary] Testing: logisticPCA-EPCA + UMAP + KMeans")
        logistic_pca_result = Mycelia.logistic_pca_epca(shuffled_binary_matrix, k=5)
        umap_model = Mycelia.umap_embed(logistic_pca_result.scores)
        fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, fit_labels)
        plt = Mycelia.plot_embeddings(umap_model.embedding;
                       title="logisticPCA-EPCA + UMAP + KMeans",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binary_method_accuracies, ("logisticPCA-EPCA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # logisticPCA + UMAP + KMedoids
    @testset "logisticPCA-EPCA + UMAP + KMedoids" begin
        println("[Binary] Testing: logisticPCA-EPCA + UMAP + KMedoids")
        logistic_pca_result = Mycelia.logistic_pca_epca(shuffled_binary_matrix, k=5)
        umap_model = Mycelia.umap_embed(logistic_pca_result.scores)
        # Compute distance matrix from UMAP embedding (Euclidean)
        embedding = umap_model.embedding
        dist_matrix = Distances.pairwise(Distances.Euclidean(), embedding; dims=2)
        kmedoids_result = Clustering.kmedoids(dist_matrix, n_distributions)
        fit_labels = kmedoids_result.assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, fit_labels)
        plt = Mycelia.plot_embeddings(embedding;
                       title="logisticPCA-EPCA + UMAP + KMedoids",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 1/3
        @test evaluation_result.macro_precision >= 1/3
        @test evaluation_result.macro_recall >= 1/3
        @test evaluation_result.accuracy >= 1/3
        push!(binary_method_accuracies, ("logisticPCA-EPCA + UMAP + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Report ranked list by accuracy
    println("\n[Binary] Accuracy ranking:")
    for (i, (name, acc)) in enumerate(sort(binary_method_accuracies, by=x->-x[2]))
        println("$(i). $(name): $(round(acc, digits=4))")
    end
end

@testset "Poisson (counts) Matrix Processing" begin
    # Set a random seed for reproducibility
    Random.seed!(42)

    # Parameters
    n_distributions = 7      # Number of distributions
    n_samples = 10      # Number of samples per distribution
    n_features = 100     # Length of each distribution (number of features)

    # Poisson (counts)
    poisson_lambdas = [rand(1:10, n_features) for _ in 1:n_distributions]
    poisson_samples = [hcat([rand.(Distributions.Poisson.(λ)) for _ in 1:n_samples]...) for λ in poisson_lambdas]
    poisson_matrix = hcat(poisson_samples...)
    poisson_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(poisson_labels))
    shuffled_poisson_matrix = poisson_matrix[:, perm]
    shuffled_poisson_labels = poisson_labels[perm]

    # Store results for ranking
    poisson_method_accuracies = []

    @testset "Sanity Check - Poisson Matrix" begin
        summary = Mycelia.sanity_check_matrix(shuffled_poisson_matrix)
        @test summary[:n_features] == n_features
        @test summary[:n_samples] == n_samples * n_distributions
        @test summary[:is_integer] == true
        @test summary[:is_nonnegative] == true
        @test summary[:is_binary] == false
        @test summary[:is_strictly_positive] == false
        @test summary[:is_in_01] == false
        @test summary[:is_probability_vector] == false
        @test summary[:suggested_epca] == :poisson_pca_epca || summary[:suggested_epca] == :negbin_pca_epca
        @test summary[:suggested_distance] == :bray_curtis_distance
    end
    # Distance clustering + Optimal Hierarchical Clustering
    @testset "Bray-Curtis Distance + Optimal Hierarchical Clustering" begin
        println("[Poisson] Testing: Bray-Curtis Distance + Optimal Hierarchical Clustering")
        poisson_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_poisson_matrix)
        poisson_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(poisson_distance_matrix)
        @test poisson_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, poisson_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(poisson_method_accuracies, ("Bray-Curtis Distance + Optimal Hierarchical Clustering", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    @testset "Bray-Curtis Distance + KMeans" begin
        println("[Poisson] Testing: Bray-Curtis Distance + KMeans")
        poisson_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_poisson_matrix)
        # Classical MDS/PCoA to get coordinates from distance matrix
        pcoa_result = Mycelia.pcoa_from_dist(poisson_distance_matrix)
        kmeans_labels = Clustering.kmeans(pcoa_result.coordinates, n_distributions).assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, kmeans_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(poisson_method_accuracies, ("Bray-Curtis Distance + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    @testset "Bray-Curtis Distance + KMedoids" begin
        println("[Poisson] Testing: Bray-Curtis Distance + KMedoids")
        poisson_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_poisson_matrix)
        kmedoids_result = Clustering.kmedoids(poisson_distance_matrix, n_distributions)
        kmedoids_labels = kmedoids_result.assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, kmedoids_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(poisson_method_accuracies, ("Bray-Curtis Distance + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # Distance + PCoA + KMeans
    @testset "Bray-Curtis Distance + PCoA + KMeans" begin
        println("[Poisson] Testing: Bray-Curtis Distance + PCoA + KMeans")
        pcoa_poisson_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_poisson_matrix))
        pcoa_fit_poisson_labels = Clustering.kmeans(pcoa_poisson_result.coordinates, n_distributions).assignments
        pcoa_fit_poisson_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, pcoa_fit_poisson_labels)
        plt = Mycelia.plot_embeddings(pcoa_poisson_result.coordinates;
                       title="Bray-Curtis Distance + PCoA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_poisson_labels,
                       fit_labels=pcoa_fit_poisson_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, pcoa_fit_poisson_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(poisson_method_accuracies, ("Bray-Curtis Distance + PCoA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Bray-Curtis Distance + PCoA + UMAP + KMeans
    @testset "Bray-Curtis Distance + PCoA + UMAP + KMeans" begin
        println("[Poisson] Testing: Bray-Curtis Distance + PCoA + UMAP + KMeans")
        pcoa_poisson_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_poisson_matrix))
        pcoa_poisson_umap_model = Mycelia.umap_embed(pcoa_poisson_result.coordinates)
        @test size(pcoa_poisson_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_poisson_umap_fit_labels = Clustering.kmeans(pcoa_poisson_umap_model.embedding, n_distributions).assignments
        pcoa_poisson_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, pcoa_poisson_umap_fit_labels)
        plt = Mycelia.plot_embeddings(pcoa_poisson_umap_model.embedding;
                       title="Bray-Curtis Distance + PCoA + UMAP + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_poisson_labels,
                       fit_labels=pcoa_poisson_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, pcoa_poisson_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(poisson_method_accuracies, ("Bray-Curtis Distance + PCoA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # PoissonPCA-EPCA + KMeans
    @testset "PoissonPCA-EPCA + KMeans" begin
        println("[Poisson] Testing: PoissonPCA-EPCA + KMeans")
        poisson_pca_result = Mycelia.poisson_pca_epca(shuffled_poisson_matrix, k=5)
        fit_labels = Clustering.kmeans(poisson_pca_result.scores, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, fit_labels)
        plt = Mycelia.plot_embeddings(poisson_pca_result.scores;
                       title="PoissonPCA-EPCA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_poisson_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(poisson_method_accuracies, ("PoissonPCA-EPCA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # PoissonPCA-EPCA + UMAP + KMeans
    @testset "PoissonPCA-EPCA + UMAP + KMeans" begin
        println("[Poisson] Testing: PoissonPCA-EPCA + UMAP + KMeans")
        poisson_pca_result = Mycelia.poisson_pca_epca(shuffled_poisson_matrix, k=5)
        umap_model = Mycelia.umap_embed(poisson_pca_result.scores)
        fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, fit_labels)
        plt = Mycelia.plot_embeddings(umap_model.embedding;
                       title="PoissonPCA-EPCA + UMAP + KMeans",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_poisson_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(poisson_method_accuracies, ("PoissonPCA-EPCA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Report ranked list by accuracy
    println("\n[Poisson] Accuracy ranking:")
    for (i, (name, acc)) in enumerate(sort(poisson_method_accuracies, by=x->-x[2]))
        println("$(i). $(name): $(round(acc, digits=4))")
    end
end

@testset "Negative Binomial (overdispersed counts) Matrix Processing" begin
    # Set a random seed for reproducibility
    Random.seed!(42)

    # Parameters
    n_distributions = 7      # Number of distributions
    n_samples = 10      # Number of samples per distribution
    n_features = 100     # Length of each distribution (number of features)

    # Negative Binomial (overdispersed counts)
    nb_r = 5  # dispersion parameter
    nb_ps = [rand(0.2:0.05:0.8, n_features) for _ in 1:n_distributions]
    nb_samples = [hcat([rand.(Distributions.NegativeBinomial.(nb_r, p)) for _ in 1:n_samples]...) for p in nb_ps]
    nb_matrix = hcat(nb_samples...)
    nb_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(nb_labels))
    shuffled_nb_matrix = nb_matrix[:, perm]
    shuffled_nb_labels = nb_labels[perm]

    # Store results for ranking
    nb_method_accuracies = []

    @testset "Sanity Check - Negative Binomial Matrix" begin
        summary = Mycelia.sanity_check_matrix(shuffled_nb_matrix)
        @test summary[:n_features] == n_features
        @test summary[:n_samples] == n_samples * n_distributions
        @test summary[:is_integer] == true
        @test summary[:is_nonnegative] == true
        @test summary[:is_binary] == false
        @test summary[:is_strictly_positive] == false
        @test summary[:is_in_01] == false
        @test summary[:is_probability_vector] == false
        @test summary[:suggested_epca] == :negbin_pca_epca
        @test summary[:suggested_distance] == :bray_curtis_distance
    end
    # Distance clustering + Optimal Hierarchical Clustering
    @testset "Bray-Curtis Distance + Optimal Hierarchical Clustering" begin
        println("[NegBin] Testing: Bray-Curtis Distance + Optimal Hierarchical Clustering")
        nb_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_nb_matrix)
        nb_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(nb_distance_matrix)
        @test nb_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, nb_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(nb_method_accuracies, ("Bray-Curtis Distance + Optimal Hierarchical Clustering", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    @testset "Bray-Curtis Distance + KMeans" begin
        println("[NegBin] Testing: Bray-Curtis Distance + KMeans")
        nb_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_nb_matrix)
        # Classical MDS/PCoA to get coordinates from distance matrix
        pcoa_result = Mycelia.pcoa_from_dist(nb_distance_matrix)
        kmeans_labels = Clustering.kmeans(pcoa_result.coordinates, n_distributions).assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, kmeans_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(nb_method_accuracies, ("Bray-Curtis Distance + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    @testset "Bray-Curtis Distance + KMedoids" begin
        println("[NegBin] Testing: Bray-Curtis Distance + KMedoids")
        nb_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_nb_matrix)
        kmedoids_result = Clustering.kmedoids(nb_distance_matrix, n_distributions)
        kmedoids_labels = kmedoids_result.assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, kmedoids_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(nb_method_accuracies, ("Bray-Curtis Distance + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Bray-Curtis Distance + PCoA + KMeans
    @testset "Bray-Curtis Distance + PCoA + KMeans" begin
        println("[NegBin] Testing: Bray-Curtis Distance + PCoA + KMeans")
        pcoa_nb_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_nb_matrix))
        pcoa_fit_nb_labels = Clustering.kmeans(pcoa_nb_result.coordinates, n_distributions).assignments
        pcoa_fit_nb_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, pcoa_fit_nb_labels)
        plt = Mycelia.plot_embeddings(pcoa_nb_result.coordinates;
                       title="Bray-Curtis Distance + PCoA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_nb_labels,
                       fit_labels=pcoa_fit_nb_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, pcoa_fit_nb_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(nb_method_accuracies, ("Bray-Curtis Distance + PCoA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Bray-Curtis Distance + PCoA + UMAP + KMeans
    @testset "Bray-Curtis Distance + PCoA + UMAP + KMeans" begin
        println("[NegBin] Testing: Bray-Curtis Distance + PCoA + UMAP + KMeans")
        pcoa_nb_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_nb_matrix))
        pcoa_nb_umap_model = Mycelia.umap_embed(pcoa_nb_result.coordinates)
        @test size(pcoa_nb_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_nb_umap_fit_labels = Clustering.kmeans(pcoa_nb_umap_model.embedding, n_distributions).assignments
        pcoa_nb_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, pcoa_nb_umap_fit_labels)
        plt = Mycelia.plot_embeddings(pcoa_nb_umap_model.embedding;
                       title="Bray-Curtis Distance + PCoA + UMAP + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_nb_labels,
                       fit_labels=pcoa_nb_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, pcoa_nb_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(nb_method_accuracies, ("Bray-Curtis Distance + PCoA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # Not working
    # # NegBinPCA-EPCA + KMeans
    # @testset "NegBinPCA-EPCA + KMeans" begin
    #     negbin_pca_result = Mycelia.negbin_pca_epca(shuffled_nb_matrix, k=5)
    #     fit_labels = Clustering.kmeans(negbin_pca_result.scores, n_distributions).assignments
    #     fit_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, fit_labels)
    #     plt = Mycelia.plot_embeddings(negbin_pca_result.scores;
    #                    title="Negative Binomial PCA-EPCA - Negative Binomial Matrix",
    #                    xlabel="PC1",
    #                    ylabel="PC2",
    #                    true_labels=shuffled_nb_labels,
    #                    fit_labels=fit_labels)
    #     display(plt)
    #     evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, fit_labels)
    #     @test evaluation_result.macro_f1 >= 2/3
    #     @test evaluation_result.macro_precision >= 2/3
    #     @test evaluation_result.macro_recall >= 2/3
    #     @test evaluation_result.accuracy >= 2/3
    #     display(evaluation_result.confusion_matrix_plot)
    #     display(evaluation_result.f1_plot)
    #     display(evaluation_result.precision_plot)
    #     display(evaluation_result.recall_plot)
    # end
    
    # # NegBinPCA-EPCA + UMAP + KMeans
    # @testset "NegBinPCA-EPCA + UMAP + KMeans" begin
    #     negbin_pca_result = Mycelia.negbin_pca_epca(shuffled_nb_matrix, k=5)
    #     umap_model = Mycelia.umap_embed(negbin_pca_result.scores)
    #     fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
    #     fit_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, fit_labels)
    #     plt = Mycelia.plot_embeddings(umap_model.embedding;
    #                    title="UMAP - Negative Binomial PCA-EPCA - Negative Binomial Matrix",
    #                    xlabel="UMAP 1",
    #                    ylabel="UMAP 2",
    #                    true_labels=shuffled_nb_labels,
    #                    fit_labels=fit_labels)
    #     display(plt)
    #     evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, fit_labels)
    #     @test evaluation_result.macro_f1 >= 2/3
    #     @test evaluation_result.macro_precision >= 2/3
    #     @test evaluation_result.macro_recall >= 2/3
    #     @test evaluation_result.accuracy >= 2/3
    #     display(evaluation_result.confusion_matrix_plot)
    #     display(evaluation_result.f1_plot)
    #     display(evaluation_result.precision_plot)
    #     display(evaluation_result.recall_plot)
    # end

    # Report ranked list by accuracy
    println("\n[NegBin] Accuracy ranking:")
    for (i, (name, acc)) in enumerate(sort(nb_method_accuracies, by=x->-x[2]))
        println("$(i). $(name): $(round(acc, digits=4))")
    end
end

@testset "Binomial (counts in 0:ntrials) Matrix Processing" begin
    # Set a random seed for reproducibility
    Random.seed!(42)

    # Parameters
    n_distributions = 7      # Number of distributions
    n_samples = 10           # Number of samples per distribution
    n_features = 100         # Length of each distribution (number of features)

    # Binomial (counts in 0:ntrials)
    ntrials = 10
    binom_ps = [rand(n_features) for _ in 1:n_distributions]
    binom_samples = [hcat([rand.(Distributions.Binomial.(ntrials, p)) for _ in 1:n_samples]...) for p in binom_ps]
    binom_matrix = hcat(binom_samples...)
    binom_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(binom_labels))
    shuffled_binom_matrix = binom_matrix[:, perm]
    shuffled_binom_labels = binom_labels[perm]

    # Store results for ranking
    binom_method_accuracies = []

    @testset "Sanity Check - Binomial Matrix" begin
        summary = Mycelia.sanity_check_matrix(shuffled_binom_matrix)
        @test summary[:n_features] == n_features
        @test summary[:n_samples] == n_samples * n_distributions
        @test summary[:is_integer] == true
        @test summary[:is_nonnegative] == true
        @test summary[:is_binary] == false
        @test summary[:is_strictly_positive] == false
        @test summary[:is_in_01] == false
        @test summary[:is_probability_vector] == false
        @test summary[:suggested_epca] == :poisson_pca_epca || summary[:suggested_epca] == :negbin_pca_epca
        @test summary[:suggested_distance] == :bray_curtis_distance
    end

    # Bray-Curtis Distance + Optimal Hierarchical Clustering
    @testset "Bray-Curtis Distance + Optimal Hierarchical Clustering" begin
        println("[Binom] Testing: Bray-Curtis Distance + Optimal Hierarchical Clustering")
        binom_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_binom_matrix)
        binom_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(binom_distance_matrix)
        @test binom_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, binom_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(binom_method_accuracies, ("Bray-Curtis Distance + Optimal Hierarchical Clustering", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    @testset "Bray-Curtis Distance + KMeans" begin
        println("[Binom] Testing: Bray-Curtis Distance + KMeans")
        binom_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_binom_matrix)
        # Classical MDS/PCoA to get coordinates from distance matrix
        pcoa_result = Mycelia.pcoa_from_dist(binom_distance_matrix)
        kmeans_labels = Clustering.kmeans(pcoa_result.coordinates, n_distributions).assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, kmeans_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binom_method_accuracies, ("Bray-Curtis Distance + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    @testset "Bray-Curtis Distance + KMedoids" begin
        println("[Binom] Testing: Bray-Curtis Distance + KMedoids")
        binom_distance_matrix = Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_binom_matrix)
        kmedoids_result = Clustering.kmedoids(binom_distance_matrix, n_distributions)
        kmedoids_labels = kmedoids_result.assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, kmedoids_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binom_method_accuracies, ("Bray-Curtis Distance + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Bray-Curtis Distance + PCoA + KMeans
    @testset "Bray-Curtis Distance + PCoA + KMeans" begin
        println("[Binom] Testing: Bray-Curtis Distance + PCoA + KMeans")
        pcoa_binom_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_binom_matrix))
        pcoa_fit_binom_labels = Clustering.kmeans(pcoa_binom_result.coordinates, n_distributions).assignments
        pcoa_fit_binom_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, pcoa_fit_binom_labels)
        plt = Mycelia.plot_embeddings(pcoa_binom_result.coordinates;
                       title="Bray-Curtis Distance + PCoA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binom_labels,
                       fit_labels=pcoa_fit_binom_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, pcoa_fit_binom_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binom_method_accuracies, ("Bray-Curtis Distance + PCoA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Bray-Curtis Distance + PCoA + UMAP + KMeans
    @testset "Bray-Curtis Distance + PCoA + UMAP + KMeans" begin
        println("[Binom] Testing: Bray-Curtis Distance + PCoA + UMAP + KMeans")
        pcoa_binom_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_bray_curtis_distance_matrix(shuffled_binom_matrix))
        pcoa_binom_umap_model = Mycelia.umap_embed(pcoa_binom_result.coordinates)
        @test size(pcoa_binom_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_binom_umap_fit_labels = Clustering.kmeans(pcoa_binom_umap_model.embedding, n_distributions).assignments
        pcoa_binom_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, pcoa_binom_umap_fit_labels)
        plt = Mycelia.plot_embeddings(pcoa_binom_umap_model.embedding;
                       title="Bray-Curtis Distance + PCoA + UMAP + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binom_labels,
                       fit_labels=pcoa_binom_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, pcoa_binom_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binom_method_accuracies, ("Bray-Curtis Distance + PCoA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # PoissonPCA-EPCA + KMeans
    @testset "PoissonPCA-EPCA + KMeans" begin
        println("[Binom] Testing: PoissonPCA-EPCA + KMeans")
        poisson_pca_result = Mycelia.poisson_pca_epca(shuffled_binom_matrix, k=5)
        fit_labels = Clustering.kmeans(poisson_pca_result.scores, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, fit_labels)
        plt = Mycelia.plot_embeddings(poisson_pca_result.scores;
                       title="PoissonPCA-EPCA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binom_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binom_method_accuracies, ("PoissonPCA-EPCA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # PoissonPCA-EPCA + UMAP + KMeans
    @testset "PoissonPCA-EPCA + UMAP + KMeans" begin
        println("[Binom] Testing: PoissonPCA-EPCA + UMAP + KMeans")
        poisson_pca_result = Mycelia.poisson_pca_epca(shuffled_binom_matrix, k=5)
        umap_model = Mycelia.umap_embed(poisson_pca_result.scores)
        fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, fit_labels)
        plt = Mycelia.plot_embeddings(umap_model.embedding;
                       title="PoissonPCA-EPCA + UMAP + KMeans",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_binom_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(binom_method_accuracies, ("PoissonPCA-EPCA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Report ranked list by accuracy
    println("\n[Binom] Accuracy ranking:")
    for (i, (name, acc)) in enumerate(sort(binom_method_accuracies, by=x->-x[2]))
        println("$(i). $(name): $(round(acc, digits=4))")
    end
end

@testset "Continuous Bernoulli (values in (0,1)) Matrix Processing" begin
    # Set a random seed for reproducibility
    Random.seed!(42)

    # Parameters
    n_distributions = 7      # Number of distributions
    n_samples = 10      # Number of samples per distribution
    n_features = 100     # Length of each distribution (number of features)

    # Continuous Bernoulli (values in (0,1))
    contb_ps = [rand(n_features) for _ in 1:n_distributions]
    contb_samples = [hcat([rand.(Distributions.Beta.(p*0.9 .+ 0.05, (1 .- p)*0.9 .+ 0.05)) for _ in 1:n_samples]...) for p in contb_ps]
    contb_matrix = hcat(contb_samples...)
    contb_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(contb_labels))
    shuffled_contb_matrix = contb_matrix[:, perm]
    shuffled_contb_labels = contb_labels[perm]
    # Ensure all values are strictly in (0, 1) for EPCA
    ϵ = 1e-6
    clipped_contb_matrix = clamp.(shuffled_contb_matrix, ϵ, 1 - ϵ)

    @testset "Sanity Check - Continuous Bernoulli Matrix" begin
        summary = Mycelia.sanity_check_matrix(clipped_contb_matrix)
        @test summary[:n_features] == n_features
        @test summary[:n_samples] == n_samples * n_distributions
        @test summary[:is_integer] == false
        @test summary[:is_nonnegative] == true
        @test summary[:is_binary] == false
        @test summary[:is_strictly_positive] == true
        @test summary[:is_in_01] == true
        @test summary[:is_probability_vector] == false
        @test summary[:suggested_epca] == :contbernoulli_pca_epca
        @test summary[:suggested_distance] == :cosine_distance
    end

    # Store results for ranking
    contb_method_accuracies = []

    # Cosine Distance + Optimal Hierarchical Clustering
    @testset "Cosine Distance + Optimal Hierarchical Clustering" begin
        println("[ContBernoulli] Testing: Cosine Distance + Optimal Hierarchical Clustering")
        contb_distance_matrix = Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_contb_matrix)
        contb_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(contb_distance_matrix)
        @test contb_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_contb_labels, contb_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_contb_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(contb_method_accuracies, ("Cosine Distance + Optimal Hierarchical Clustering", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Cosine Distance + KMeans
    @testset "Cosine Distance + KMeans" begin
        println("[ContBernoulli] Testing: Cosine Distance + KMeans")
        contb_distance_matrix = Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_contb_matrix)
        # Classical MDS/PCoA to get coordinates from distance matrix
        pcoa_result = Mycelia.pcoa_from_dist(contb_distance_matrix)
        kmeans_labels = Clustering.kmeans(pcoa_result.coordinates, n_distributions).assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_contb_labels, kmeans_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_contb_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(contb_method_accuracies, ("Cosine Distance + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # Cosine Distance + KMedoids
    @testset "Cosine Distance + KMedoids" begin
        println("[ContBernoulli] Testing: Cosine Distance + KMedoids")
        contb_distance_matrix = Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_contb_matrix)
        kmedoids_result = Clustering.kmedoids(contb_distance_matrix, n_distributions)
        kmedoids_labels = kmedoids_result.assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_contb_labels, kmedoids_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_contb_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(contb_method_accuracies, ("Cosine Distance + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Cosine Distance + PCoA + KMeans
    @testset "Cosine Distance + PCoA + KMeans" begin
        println("[ContBernoulli] Testing: Cosine Distance + PCoA + KMeans")
        pcoa_contb_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_contb_matrix))
        pcoa_fit_contb_labels = Clustering.kmeans(pcoa_contb_result.coordinates, n_distributions).assignments
        pcoa_fit_contb_labels, mapping = Mycelia.best_label_mapping(shuffled_contb_labels, pcoa_fit_contb_labels)
        plt = Mycelia.plot_embeddings(pcoa_contb_result.coordinates;
                       title="Cosine Distance + PCoA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_contb_labels,
                       fit_labels=pcoa_fit_contb_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_contb_labels, pcoa_fit_contb_labels)
        @test evaluation_result.macro_f1 >= 1/3
        @test evaluation_result.macro_precision >= 1/3
        @test evaluation_result.macro_recall >= 1/3
        @test evaluation_result.accuracy >= 1/3
        push!(contb_method_accuracies, ("Cosine Distance + PCoA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Cosine Distance + PCoA + UMAP + KMeans
    @testset "Cosine Distance + PCoA + UMAP + KMeans" begin
        println("[ContBernoulli] Testing: Cosine Distance + PCoA + UMAP + KMeans")
        pcoa_contb_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_contb_matrix))
        pcoa_contb_umap_model = Mycelia.umap_embed(pcoa_contb_result.coordinates)
        @test size(pcoa_contb_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_contb_umap_fit_labels = Clustering.kmeans(pcoa_contb_umap_model.embedding, n_distributions).assignments
        pcoa_contb_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_contb_labels, pcoa_contb_umap_fit_labels)
        plt = Mycelia.plot_embeddings(pcoa_contb_umap_model.embedding;
                       title="Cosine Distance + PCoA + UMAP + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_contb_labels,
                       fit_labels=pcoa_contb_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_contb_labels, pcoa_contb_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(contb_method_accuracies, ("Cosine Distance + PCoA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # ContBernoulliPCA-EPCA + KMeans
    @testset "ContBernoulliPCA-EPCA + KMeans" begin
        println("[ContBernoulli] Testing: ContBernoulliPCA-EPCA + KMeans")
        contb_pca_result = Mycelia.contbernoulli_pca_epca(clipped_contb_matrix, k=5)
        fit_labels = Clustering.kmeans(contb_pca_result.scores, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_contb_labels, fit_labels)
        plt = Mycelia.plot_embeddings(contb_pca_result.scores;
                       title="ContBernoulliPCA-EPCA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_contb_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_contb_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 1/3
        @test evaluation_result.macro_precision >= 1/3
        @test evaluation_result.macro_recall >= 1/3
        @test evaluation_result.accuracy >= 1/3
        push!(contb_method_accuracies, ("ContBernoulliPCA-EPCA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # ContBernoulliPCA-EPCA + UMAP + KMeans
    @testset "ContBernoulliPCA-EPCA + UMAP + KMeans" begin
        println("[ContBernoulli] Testing: ContBernoulliPCA-EPCA + UMAP + KMeans")
        contb_pca_result = Mycelia.contbernoulli_pca_epca(clipped_contb_matrix, k=5)
        umap_model = Mycelia.umap_embed(contb_pca_result.scores)
        fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_contb_labels, fit_labels)
        plt = Mycelia.plot_embeddings(umap_model.embedding;
                       title="ContBernoulliPCA-EPCA + UMAP + KMeans",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_contb_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_contb_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(contb_method_accuracies, ("ContBernoulliPCA-EPCA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Report ranked list by accuracy
    println("\n[ContBernoulli] Accuracy ranking:")
    for (i, (name, acc)) in enumerate(sort(contb_method_accuracies, by=x->-x[2]))
        println("$(i). $(name): $(round(acc, digits=4))")
    end
end

@testset "Gamma (strictly positive) Matrix Processing" begin
    # Set a random seed for reproducibility
    Random.seed!(42)

    # Parameters
    n_distributions = 7      # Number of distributions
    n_samples = 10      # Number of samples per distribution
    n_features = 100     # Length of each distribution (number of features)

    # Gamma (strictly positive)
    gamma_shapes = [rand(1.0:0.5:5.0, n_features) for _ in 1:n_distributions]
    gamma_scales = [rand(1.0:0.5:3.0, n_features) for _ in 1:n_distributions]
    gamma_samples = [hcat([rand.(Distributions.Gamma.(sh, sc)) for _ in 1:n_samples]...) for (sh, sc) in zip(gamma_shapes, gamma_scales)]
    gamma_matrix = hcat(gamma_samples...)
    gamma_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(gamma_labels))
    shuffled_gamma_matrix = gamma_matrix[:, perm]
    shuffled_gamma_labels = gamma_labels[perm]
    # Ensure all values are strictly positive for GammaPCA-EPCA
    ϵ = 1e-6
    clipped_gamma_matrix = clamp.(shuffled_gamma_matrix, ϵ, Inf)

    @testset "Sanity Check - Gamma Matrix" begin
        summary = Mycelia.sanity_check_matrix(clipped_gamma_matrix)
        @test summary[:n_features] == n_features
        @test summary[:n_samples] == n_samples * n_distributions
        @test summary[:is_integer] == false
        @test summary[:is_nonnegative] == true
        @test summary[:is_binary] == false
        @test summary[:is_strictly_positive] == true
        @test summary[:is_in_01] == false
        @test summary[:is_probability_vector] == false
        @test summary[:suggested_epca] == :gamma_pca_epca
        @test summary[:suggested_distance] == :cosine_distance
    end

    # Store results for ranking
    gamma_method_accuracies = []

    # Cosine Distance + Optimal Hierarchical Clustering
    @testset "Cosine Distance + Optimal Hierarchical Clustering" begin
        println("[Gamma] Testing: Cosine Distance + Optimal Hierarchical Clustering")
        gamma_distance_matrix = Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_gamma_matrix)
        gamma_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(gamma_distance_matrix)
        @test gamma_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_gamma_labels, gamma_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gamma_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(gamma_method_accuracies, ("Cosine Distance + Optimal Hierarchical Clustering", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Cosine Distance + KMeans
    @testset "Cosine Distance + KMeans" begin
        println("[Gamma] Testing: Cosine Distance + KMeans")
        gamma_distance_matrix = Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_gamma_matrix)
        # Classical MDS/PCoA to get coordinates from distance matrix
        pcoa_result = Mycelia.pcoa_from_dist(gamma_distance_matrix)
        kmeans_labels = Clustering.kmeans(pcoa_result.coordinates, n_distributions).assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_gamma_labels, kmeans_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gamma_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(gamma_method_accuracies, ("Cosine Distance + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # Cosine Distance + KMedoids
    @testset "Cosine Distance + KMedoids" begin
        println("[Gamma] Testing: Cosine Distance + KMedoids")
        gamma_distance_matrix = Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_gamma_matrix)
        kmedoids_result = Clustering.kmedoids(gamma_distance_matrix, n_distributions)
        kmedoids_labels = kmedoids_result.assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_gamma_labels, kmedoids_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gamma_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(gamma_method_accuracies, ("Cosine Distance + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Cosine Distance + PCoA + KMeans
    @testset "Cosine Distance + PCoA + KMeans" begin
        println("[Gamma] Testing: Cosine Distance + PCoA + KMeans")
        pcoa_gamma_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_gamma_matrix))
        pcoa_fit_gamma_labels = Clustering.kmeans(pcoa_gamma_result.coordinates, n_distributions).assignments
        pcoa_fit_gamma_labels, mapping = Mycelia.best_label_mapping(shuffled_gamma_labels, pcoa_fit_gamma_labels)
        plt = Mycelia.plot_embeddings(pcoa_gamma_result.coordinates;
                       title="Cosine Distance + PCoA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_gamma_labels,
                       fit_labels=pcoa_fit_gamma_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gamma_labels, pcoa_fit_gamma_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(gamma_method_accuracies, ("Cosine Distance + PCoA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Cosine Distance + PCoA + UMAP + KMeans
    @testset "Cosine Distance + PCoA + UMAP + KMeans" begin
        println("[Gamma] Testing: Cosine Distance + PCoA + UMAP + KMeans")
        pcoa_gamma_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_cosine_distance_matrix(clipped_gamma_matrix))
        pcoa_gamma_umap_model = Mycelia.umap_embed(pcoa_gamma_result.coordinates)
        @test size(pcoa_gamma_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_gamma_umap_fit_labels = Clustering.kmeans(pcoa_gamma_umap_model.embedding, n_distributions).assignments
        pcoa_gamma_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_gamma_labels, pcoa_gamma_umap_fit_labels)
        plt = Mycelia.plot_embeddings(pcoa_gamma_umap_model.embedding;
                       title="Cosine Distance + PCoA + UMAP + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_gamma_labels,
                       fit_labels=pcoa_gamma_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gamma_labels, pcoa_gamma_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(gamma_method_accuracies, ("Cosine Distance + PCoA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # # GammaPCA-EPCA + KMeans
    # @testset "GammaPCA-EPCA + KMeans" begin
    #     println("[Gamma] Testing: GammaPCA-EPCA + KMeans")
    #     gamma_pca_result = Mycelia.gamma_pca_epca(clipped_gamma_matrix, k=5)
    #     fit_labels = Clustering.kmeans(gamma_pca_result.scores, n_distributions).assignments
    #     fit_labels, mapping = Mycelia.best_label_mapping(shuffled_gamma_labels, fit_labels)
    #     plt = Mycelia.plot_embeddings(gamma_pca_result.scores;
    #                    title="GammaPCA-EPCA + KMeans",
    #                    xlabel="PC1",
    #                    ylabel="PC2",
    #                    true_labels=shuffled_gamma_labels,
    #                    fit_labels=fit_labels)
    #     display(plt)
    #     evaluation_result = Mycelia.evaluate_classification(shuffled_gamma_labels, fit_labels)
    #     @test evaluation_result.macro_f1 >= 2/3
    #     @test evaluation_result.macro_precision >= 2/3
    #     @test evaluation_result.macro_recall >= 2/3
    #     @test evaluation_result.accuracy >= 2/3
    #     push!(gamma_method_accuracies, ("GammaPCA-EPCA + KMeans", evaluation_result.accuracy))
    #     display(evaluation_result.confusion_matrix_plot)
    #     display(evaluation_result.f1_plot)
    #     display(evaluation_result.precision_plot)
    #     display(evaluation_result.recall_plot)
    # end

    # # GammaPCA-EPCA + UMAP + KMeans
    # @testset "GammaPCA-EPCA + UMAP + KMeans" begin
    #     println("[Gamma] Testing: GammaPCA-EPCA + UMAP + KMeans")
    #     gamma_pca_result = Mycelia.gamma_pca_epca(clipped_gamma_matrix, k=5)
    #     umap_model = Mycelia.umap_embed(gamma_pca_result.scores)
    #     fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
    #     fit_labels, mapping = Mycelia.best_label_mapping(shuffled_gamma_labels, fit_labels)
    #     plt = Mycelia.plot_embeddings(umap_model.embedding;
    #                    title="GammaPCA-EPCA + UMAP + KMeans",
    #                    xlabel="UMAP 1",
    #                    ylabel="UMAP 2",
    #                    true_labels=shuffled_gamma_labels,
    #                    fit_labels=fit_labels)
    #     display(plt)
    #     evaluation_result = Mycelia.evaluate_classification(shuffled_gamma_labels, fit_labels)
    #     @test evaluation_result.macro_f1 >= 2/3
    #     @test evaluation_result.macro_precision >= 2/3
    #     @test evaluation_result.macro_recall >= 2/3
    #     @test evaluation_result.accuracy >= 2/3
    #     push!(gamma_method_accuracies, ("GammaPCA-EPCA + UMAP + KMeans", evaluation_result.accuracy))
    #     display(evaluation_result.confusion_matrix_plot)
    #     display(evaluation_result.f1_plot)
    #     display(evaluation_result.precision_plot)
    #     display(evaluation_result.recall_plot)
    # end

    # Report ranked list by accuracy
    println("\n[Gamma] Accuracy ranking:")
    for (i, (name, acc)) in enumerate(sort(gamma_method_accuracies, by=x->-x[2]))
        println("$(i). $(name): $(round(acc, digits=4))")
    end
end

@testset "Gaussian (centered, real-valued) Matrix Processing" begin
    # Set a random seed for reproducibility
    Random.seed!(42)

    # Parameters
    n_distributions = 7      # Number of distributions
    n_samples = 10      # Number of samples per distribution
    n_features = 100     # Length of each distribution (number of features)

    # Gaussian (centered, real-valued)
    gauss_means = [randn(n_features) for _ in 1:n_distributions]
    gauss_stds = [rand(0.5:0.1:2.0, n_features) for _ in 1:n_distributions]
    gauss_samples = [hcat([rand.(Distributions.Normal.(μ, σ)) for _ in 1:n_samples]...) for (μ, σ) in zip(gauss_means, gauss_stds)]
    gauss_matrix = hcat(gauss_samples...)
    gauss_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(gauss_labels))
    shuffled_gauss_matrix = gauss_matrix[:, perm]
    shuffled_gauss_labels = gauss_labels[perm]

    @testset "Sanity Check - Gaussian Matrix" begin
        summary = Mycelia.sanity_check_matrix(shuffled_gauss_matrix)
        @test summary[:n_features] == n_features
        @test summary[:n_samples] == n_samples * n_distributions
        @test summary[:is_integer] == false
        @test summary[:is_binary] == false
        # Gaussian can have negative values, so is_nonnegative and is_strictly_positive are likely false
        @test summary[:is_nonnegative] == false
        @test summary[:is_strictly_positive] == false
        @test summary[:is_in_01] == false
        @test summary[:is_probability_vector] == false
        @test summary[:suggested_epca] == :gaussian_pca_epca || summary[:suggested_epca] == :pca_transform
        @test summary[:suggested_distance] == :euclidean_distance
    end

    # Store results for ranking
    gauss_method_accuracies = []

    # Distance clustering + Optimal Hierarchical Clustering
    @testset "Euclidean Distance + Optimal Hierarchical Clustering" begin
        println("[Gaussian] Testing: Euclidean Distance + Optimal Hierarchical Clustering")
        gauss_distance_matrix = Mycelia.frequency_matrix_to_euclidean_distance_matrix(shuffled_gauss_matrix)
        gauss_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(gauss_distance_matrix)
        @test gauss_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_gauss_labels, gauss_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gauss_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(gauss_method_accuracies, ("Euclidean Distance + Optimal Hierarchical Clustering", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # Distance clustering + KMeans
    @testset "Euclidean Distance + KMeans" begin
        println("[Gaussian] Testing: Euclidean Distance + KMeans")
        gauss_distance_matrix = Mycelia.frequency_matrix_to_euclidean_distance_matrix(shuffled_gauss_matrix)
        # Classical MDS/PCoA to get coordinates from distance matrix
        pcoa_result = Mycelia.pcoa_from_dist(gauss_distance_matrix)
        kmeans_labels = Clustering.kmeans(pcoa_result.coordinates, n_distributions).assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_gauss_labels, kmeans_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gauss_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 0.95
        @test evaluation_result.macro_precision >= 0.95
        @test evaluation_result.macro_recall >= 0.95
        @test evaluation_result.accuracy >= 0.95
        push!(gauss_method_accuracies, ("Euclidean Distance + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # Distance clustering + KMedoids
    @testset "Euclidean Distance + KMedoids" begin
        println("[Gaussian] Testing: Euclidean Distance + KMedoids")
        gauss_distance_matrix = Mycelia.frequency_matrix_to_euclidean_distance_matrix(shuffled_gauss_matrix)
        kmedoids_result = Clustering.kmedoids(gauss_distance_matrix, n_distributions)
        kmedoids_labels = kmedoids_result.assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_gauss_labels, kmedoids_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gauss_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 0.95
        @test evaluation_result.macro_precision >= 0.95
        @test evaluation_result.macro_recall >= 0.95
        @test evaluation_result.accuracy >= 0.95
        push!(gauss_method_accuracies, ("Euclidean Distance + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Distance + PCoA + KMeans
    @testset "Euclidean Distance + PCoA + KMeans" begin
        println("[Gaussian] Testing: Euclidean Distance + PCoA + KMeans")
        pcoa_gauss_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_euclidean_distance_matrix(shuffled_gauss_matrix))
        pcoa_fit_gauss_labels = Clustering.kmeans(pcoa_gauss_result.coordinates, n_distributions).assignments
        pcoa_fit_gauss_labels, mapping = Mycelia.best_label_mapping(shuffled_gauss_labels, pcoa_fit_gauss_labels)
        plt = Mycelia.plot_embeddings(pcoa_gauss_result.coordinates;
                       title="Euclidean Distance + PCoA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_gauss_labels,
                       fit_labels=pcoa_fit_gauss_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gauss_labels, pcoa_fit_gauss_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(gauss_method_accuracies, ("Euclidean Distance + PCoA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Euclidean Distance + PCoA + UMAP + KMeans
    @testset "Euclidean Distance + PCoA + UMAP + KMeans" begin
        println("[Gaussian] Testing: Euclidean Distance + PCoA + UMAP + KMeans")
        pcoa_gauss_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_euclidean_distance_matrix(shuffled_gauss_matrix))
        pcoa_gauss_umap_model = Mycelia.umap_embed(pcoa_gauss_result.coordinates)
        @test size(pcoa_gauss_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_gauss_umap_fit_labels = Clustering.kmeans(pcoa_gauss_umap_model.embedding, n_distributions).assignments
        pcoa_gauss_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_gauss_labels, pcoa_gauss_umap_fit_labels)
        plt = Mycelia.plot_embeddings(pcoa_gauss_umap_model.embedding;
                       title="Euclidean Distance + PCoA + UMAP + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_gauss_labels,
                       fit_labels=pcoa_gauss_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gauss_labels, pcoa_gauss_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(gauss_method_accuracies, ("Euclidean Distance + PCoA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # GaussianPCA-EPCA + KMeans
    @testset "GaussianPCA-EPCA + KMeans" begin
        println("[Gaussian] Testing: GaussianPCA-EPCA + KMeans")
        gauss_pca_result = Mycelia.gaussian_pca_epca(shuffled_gauss_matrix, k=5)
        fit_labels = Clustering.kmeans(gauss_pca_result.scores, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_gauss_labels, fit_labels)
        plt = Mycelia.plot_embeddings(gauss_pca_result.scores;
                       title="GaussianPCA-EPCA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_gauss_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gauss_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(gauss_method_accuracies, ("GaussianPCA-EPCA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # GaussianPCA-EPCA + UMAP + KMeans
    @testset "GaussianPCA-EPCA + UMAP + KMeans" begin
        println("[Gaussian] Testing: GaussianPCA-EPCA + UMAP + KMeans")
        gauss_pca_result = Mycelia.gaussian_pca_epca(shuffled_gauss_matrix, k=5)
        umap_model = Mycelia.umap_embed(gauss_pca_result.scores)
        fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_gauss_labels, fit_labels)
        plt = Mycelia.plot_embeddings(umap_model.embedding;
                       title="GaussianPCA-EPCA + UMAP + KMeans",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_gauss_labels,
                       fit_labels=fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_gauss_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(gauss_method_accuracies, ("GaussianPCA-EPCA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Report ranked list by accuracy
    println("\n[Gaussian] Accuracy ranking:")
    for (i, (name, acc)) in enumerate(sort(gauss_method_accuracies, by=x->-x[2]))
        println("$(i). $(name): $(round(acc, digits=4))")
    end
end

@testset "Probability Vector (Compositional) Matrix Processing" begin
    # Set a random seed for reproducibility
    Random.seed!(42)

    # Parameters
    n_distributions = 7      # Number of distributions
    n_samples = 10           # Number of samples per distribution
    n_features = 100         # Length of each distribution (number of features)

    # Dirichlet (probability vectors: non-negative, sum to 1)
    dirichlet_alphas = [rand(0.5:0.1:2.0, n_features) for _ in 1:n_distributions]
    dirichlet_samples = [hcat([rand(Distributions.Dirichlet(alpha)) for _ in 1:n_samples]...) for alpha in dirichlet_alphas]
    dirichlet_matrix = hcat(dirichlet_samples...)
    dirichlet_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(dirichlet_labels))
    shuffled_dirichlet_matrix = dirichlet_matrix[:, perm]
    shuffled_dirichlet_labels = dirichlet_labels[perm]

    @testset "Sanity Check - Probability Vector Matrix" begin
        summary = Mycelia.sanity_check_matrix(shuffled_dirichlet_matrix)
        @test summary[:n_features] == n_features
        @test summary[:n_samples] == n_samples * n_distributions
        @test summary[:is_integer] == false
        @test summary[:is_nonnegative] == true
        @test summary[:is_binary] == false
        @test summary[:is_strictly_positive] == true
        @test summary[:is_in_01] == true
        @test summary[:is_probability_vector] == true
        # Each column should sum to 1 (within tolerance)
        @test all(abs.(sum(shuffled_dirichlet_matrix, dims=1) .- 1) .< 1e-8)
        # No direct EPCA for compositional/probability data
        @test summary[:suggested_epca] === nothing
        @test summary[:suggested_distance] == :jensen_shannon_divergence
    end

    # Store results for ranking
    probvec_method_accuracies = []

    # Distance clustering + Optimal Hierarchical Clustering
    @testset "Jensen-Shannon Divergence + Optimal Hierarchical Clustering" begin
        println("[ProbVec] Testing: Jensen-Shannon Divergence + Optimal Hierarchical Clustering")
        probvec_distance_matrix = Mycelia.frequency_matrix_to_jensen_shannon_distance_matrix(shuffled_dirichlet_matrix)
        probvec_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(probvec_distance_matrix)
        @test probvec_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_dirichlet_labels, probvec_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_dirichlet_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        push!(probvec_method_accuracies, ("Jensen-Shannon Divergence + Optimal Hierarchical Clustering", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    @testset "Jensen-Shannon Divergence + KMeans" begin
        println("[ProbVec] Testing: Jensen-Shannon Divergence + KMeans")
        probvec_distance_matrix = Mycelia.frequency_matrix_to_jensen_shannon_distance_matrix(shuffled_dirichlet_matrix)
        # Classical MDS/PCoA to get coordinates from distance matrix
        pcoa_result = Mycelia.pcoa_from_dist(probvec_distance_matrix)
        kmeans_labels = Clustering.kmeans(pcoa_result.coordinates, n_distributions).assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_dirichlet_labels, kmeans_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_dirichlet_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(probvec_method_accuracies, ("Jensen-Shannon Divergence + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    @testset "Jensen-Shannon Divergence + KMedoids" begin
        println("[ProbVec] Testing: Jensen-Shannon Divergence + KMedoids")
        probvec_distance_matrix = Mycelia.frequency_matrix_to_jensen_shannon_distance_matrix(shuffled_dirichlet_matrix)
        kmedoids_result = Clustering.kmedoids(probvec_distance_matrix, n_distributions)
        kmedoids_labels = kmedoids_result.assignments
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_dirichlet_labels, kmedoids_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_dirichlet_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(probvec_method_accuracies, ("Jensen-Shannon Divergence + KMedoids", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Distance + PCoA + KMeans
    @testset "Jensen-Shannon Divergence + PCoA + KMeans" begin
        println("[ProbVec] Testing: Jensen-Shannon Divergence + PCoA + KMeans")
        pcoa_probvec_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_jensen_shannon_distance_matrix(shuffled_dirichlet_matrix))
        pcoa_fit_probvec_labels = Clustering.kmeans(pcoa_probvec_result.coordinates, n_distributions).assignments
        pcoa_fit_probvec_labels, mapping = Mycelia.best_label_mapping(shuffled_dirichlet_labels, pcoa_fit_probvec_labels)
        plt = Mycelia.plot_embeddings(pcoa_probvec_result.coordinates;
                       title="Jensen-Shannon Divergence + PCoA + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_dirichlet_labels,
                       fit_labels=pcoa_fit_probvec_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_dirichlet_labels, pcoa_fit_probvec_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(probvec_method_accuracies, ("Jensen-Shannon Divergence + PCoA + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Jensen-Shannon Divergence + PCoA + UMAP + KMeans
    @testset "Jensen-Shannon Divergence + PCoA + UMAP + KMeans" begin
        println("[ProbVec] Testing: Jensen-Shannon Divergence + PCoA + UMAP + KMeans")
        pcoa_probvec_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_jensen_shannon_distance_matrix(shuffled_dirichlet_matrix))
        pcoa_probvec_umap_model = Mycelia.umap_embed(pcoa_probvec_result.coordinates)
        @test size(pcoa_probvec_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_probvec_umap_fit_labels = Clustering.kmeans(pcoa_probvec_umap_model.embedding, n_distributions).assignments
        pcoa_probvec_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_dirichlet_labels, pcoa_probvec_umap_fit_labels)
        plt = Mycelia.plot_embeddings(pcoa_probvec_umap_model.embedding;
                       title="Jensen-Shannon Divergence + PCoA + UMAP + KMeans",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_dirichlet_labels,
                       fit_labels=pcoa_probvec_umap_fit_labels)
        display(plt)
        evaluation_result = Mycelia.evaluate_classification(shuffled_dirichlet_labels, pcoa_probvec_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        push!(probvec_method_accuracies, ("Jensen-Shannon Divergence + PCoA + UMAP + KMeans", evaluation_result.accuracy))
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # No direct EPCA for compositional/probability data
    # Report ranked list by accuracy
    println("\n[ProbVec] Accuracy ranking:")
    for (i, (name, acc)) in enumerate(sort(probvec_method_accuracies, by=x->-x[2]))
        println("$(i). $(name): $(round(acc, digits=4))")
    end
end