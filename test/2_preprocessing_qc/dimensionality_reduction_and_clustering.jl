import Mycelia
import Random
import Distributions
import MultivariateStats
import UMAP
import Statistics
import LinearAlgebra
import Plots
import Clustering

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

    # Distance clustering
    @testset "Distance Clustering (Jaccard + Optimal Clusters)" begin
        binary_distance_matrix = Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix)
        binary_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(binary_distance_matrix)
        @test binary_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, binary_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Distance + PCoA
    @testset "Distance + PCoA + KMeans" begin
        pcoa_binary_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix))
        pcoa_fit_binary_labels = Clustering.kmeans(pcoa_binary_result.coordinates, n_distributions).assignments
        pcoa_fit_binary_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, pcoa_fit_binary_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, pcoa_fit_binary_labels)
        @test evaluation_result.macro_f1 >= .5
        @test evaluation_result.macro_precision >= .5
        @test evaluation_result.macro_recall >= .5
        @test evaluation_result.accuracy >= .5
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(pcoa_binary_result.coordinates;
                       title="PCoA - Jaccard Distance - Binary Matrix",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=pcoa_fit_binary_labels)
        display(plt)
    end

    # Distance + PCoA + UMAP
    @testset "Distance + PCoA + UMAP + KMeans" begin
        pcoa_binary_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_jaccard_distance_matrix(shuffled_binary_matrix))
        pcoa_binary_umap_model = Mycelia.umap_embed(pcoa_binary_result.coordinates)
        @test size(pcoa_binary_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_binary_umap_fit_labels = Clustering.kmeans(pcoa_binary_umap_model.embedding, n_distributions).assignments
        pcoa_binary_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, pcoa_binary_umap_fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, pcoa_binary_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(pcoa_binary_umap_model.embedding;
                       title="PCoA - Jaccard Distance - Binary Matrix - UMAP",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=pcoa_binary_umap_fit_labels)
        display(plt)
    end

    # logisticPCA
    @testset "logisticPCA-EPCA + KMeans" begin
        logistic_pca_result = Mycelia.logistic_pca_epca(shuffled_binary_matrix, k=5)
        # @test size(logistic_pca_result.scores) == (5, n_samples * n_distributions)
        # @test size(logistic_pca_result.loadings) == (5, n_features)
        fit_labels = Clustering.kmeans(logistic_pca_result.scores, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 0.0
        @test evaluation_result.macro_precision >= 0.0
        @test evaluation_result.macro_recall >= 0.0
        @test evaluation_result.accuracy >= 0.0
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(logistic_pca_result.scores;
                       title="Logistic PCA-EPCA - Binary Matrix",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=fit_labels)
        display(plt)
    end

    # logisticPCA + UMAP
    @testset "logisticPCA-EPCA + UMAP + KMeans" begin
        logistic_pca_result = Mycelia.logistic_pca_epca(shuffled_binary_matrix, k=5)
        umap_model = Mycelia.umap_embed(logistic_pca_result.scores)
        # @test size(umap_result) == (2, n_samples * n_distributions)
        fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binary_labels, fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binary_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(umap_model.embedding;
                       title="UMAP - Logistic PCA-EPCA - Binary Matrix",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_binary_labels,
                       fit_labels=fit_labels)
        display(plt)
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
    # Distance clustering
    @testset "Distance Clustering (Bray-Curtis + Optimal Clusters)" begin
        poisson_distance_matrix = Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_poisson_matrix)
        poisson_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(poisson_distance_matrix)
        @test poisson_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, poisson_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # Distance + PCoA + KMeans
    @testset "Distance + PCoA + KMeans" begin
        pcoa_poisson_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_poisson_matrix))
        pcoa_fit_poisson_labels = Clustering.kmeans(pcoa_poisson_result.coordinates, n_distributions).assignments
        pcoa_fit_poisson_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, pcoa_fit_poisson_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, pcoa_fit_poisson_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(pcoa_poisson_result.coordinates;
                       title="PCoA - Bray-Curtis Distance - Poisson Matrix",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_poisson_labels,
                       fit_labels=pcoa_fit_poisson_labels)
        display(plt)
    end
    
    # Distance + PCoA + UMAP + KMeans
    @testset "Distance + PCoA + UMAP + KMeans" begin
        pcoa_poisson_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_poisson_matrix))
        pcoa_poisson_umap_model = Mycelia.umap_embed(pcoa_poisson_result.coordinates)
        @test size(pcoa_poisson_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_poisson_umap_fit_labels = Clustering.kmeans(pcoa_poisson_umap_model.embedding, n_distributions).assignments
        pcoa_poisson_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, pcoa_poisson_umap_fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, pcoa_poisson_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(pcoa_poisson_umap_model.embedding;
                       title="PCoA - Bray-Curtis Distance - Poisson Matrix - UMAP",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_poisson_labels,
                       fit_labels=pcoa_poisson_umap_fit_labels)
        display(plt)
    end
    
    # PoissonPCA-EPCA + KMeans
    @testset "PoissonPCA-EPCA + KMeans" begin
        poisson_pca_result = Mycelia.poisson_pca_epca(shuffled_poisson_matrix, k=5)
        fit_labels = Clustering.kmeans(poisson_pca_result.scores, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(poisson_pca_result.scores;
                       title="Poisson PCA-EPCA - Poisson Matrix",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_poisson_labels,
                       fit_labels=fit_labels)
        display(plt)
    end
    
    # PoissonPCA-EPCA + UMAP + KMeans
    @testset "PoissonPCA-EPCA + UMAP + KMeans" begin
        poisson_pca_result = Mycelia.poisson_pca_epca(shuffled_poisson_matrix, k=5)
        umap_model = Mycelia.umap_embed(poisson_pca_result.scores)
        fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_poisson_labels, fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_poisson_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(umap_model.embedding;
                       title="UMAP - Poisson PCA-EPCA - Poisson Matrix",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_poisson_labels,
                       fit_labels=fit_labels)
        display(plt)
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
    # Distance clustering
    @testset "Distance Clustering (Bray-Curtis + Optimal Clusters)" begin
        nb_distance_matrix = Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_nb_matrix)
        nb_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(nb_distance_matrix)
        @test nb_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, nb_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end
    
    # Distance + PCoA + KMeans
    @testset "Distance + PCoA + KMeans" begin
        pcoa_nb_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_nb_matrix))
        pcoa_fit_nb_labels = Clustering.kmeans(pcoa_nb_result.coordinates, n_distributions).assignments
        pcoa_fit_nb_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, pcoa_fit_nb_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, pcoa_fit_nb_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(pcoa_nb_result.coordinates;
                       title="PCoA - Bray-Curtis Distance - Negative Binomial Matrix",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_nb_labels,
                       fit_labels=pcoa_fit_nb_labels)
        display(plt)
    end
    
    # Distance + PCoA + UMAP + KMeans
    @testset "Distance + PCoA + UMAP + KMeans" begin
        pcoa_nb_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_nb_matrix))
        pcoa_nb_umap_model = Mycelia.umap_embed(pcoa_nb_result.coordinates)
        @test size(pcoa_nb_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_nb_umap_fit_labels = Clustering.kmeans(pcoa_nb_umap_model.embedding, n_distributions).assignments
        pcoa_nb_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, pcoa_nb_umap_fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, pcoa_nb_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(pcoa_nb_umap_model.embedding;
                       title="PCoA - Bray-Curtis Distance - Negative Binomial Matrix - UMAP",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_nb_labels,
                       fit_labels=pcoa_nb_umap_fit_labels)
        display(plt)
    end
    
    # Not working
    # # NegBinPCA-EPCA + KMeans
    # @testset "NegBinPCA-EPCA + KMeans" begin
    #     negbin_pca_result = Mycelia.negbin_pca_epca(shuffled_nb_matrix, k=5)
    #     fit_labels = Clustering.kmeans(negbin_pca_result.scores, n_distributions).assignments
    #     fit_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, fit_labels)
    #     evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, fit_labels)
    #     @test evaluation_result.macro_f1 >= 2/3
    #     @test evaluation_result.macro_precision >= 2/3
    #     @test evaluation_result.macro_recall >= 2/3
    #     @test evaluation_result.accuracy >= 2/3
    #     display(evaluation_result.confusion_matrix_plot)
    #     display(evaluation_result.f1_plot)
    #     display(evaluation_result.precision_plot)
    #     display(evaluation_result.recall_plot)
    #     plt = Mycelia.plot_embeddings(negbin_pca_result.scores;
    #                    title="Negative Binomial PCA-EPCA - Negative Binomial Matrix",
    #                    xlabel="PC1",
    #                    ylabel="PC2",
    #                    true_labels=shuffled_nb_labels,
    #                    fit_labels=fit_labels)
    #     display(plt)
    # end
    
    # # NegBinPCA-EPCA + UMAP + KMeans
    # @testset "NegBinPCA-EPCA + UMAP + KMeans" begin
    #     negbin_pca_result = Mycelia.negbin_pca_epca(shuffled_nb_matrix, k=5)
    #     umap_model = Mycelia.umap_embed(negbin_pca_result.scores)
    #     fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
    #     fit_labels, mapping = Mycelia.best_label_mapping(shuffled_nb_labels, fit_labels)
    #     evaluation_result = Mycelia.evaluate_classification(shuffled_nb_labels, fit_labels)
    #     @test evaluation_result.macro_f1 >= 2/3
    #     @test evaluation_result.macro_precision >= 2/3
    #     @test evaluation_result.macro_recall >= 2/3
    #     @test evaluation_result.accuracy >= 2/3
    #     display(evaluation_result.confusion_matrix_plot)
    #     display(evaluation_result.f1_plot)
    #     display(evaluation_result.precision_plot)
    #     display(evaluation_result.recall_plot)
    #     plt = Mycelia.plot_embeddings(umap_model.embedding;
    #                    title="UMAP - Negative Binomial PCA-EPCA - Negative Binomial Matrix",
    #                    xlabel="UMAP 1",
    #                    ylabel="UMAP 2",
    #                    true_labels=shuffled_nb_labels,
    #                    fit_labels=fit_labels)
    #     display(plt)
    # end
end

# need to test and complete from here down!
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
    binom_samples = [hcat([rand.(Distributions.Binomial(ntrials, p)) for _ in 1:n_samples]...) for p in binom_ps]
    binom_matrix = hcat(binom_samples...)
    binom_labels = repeat(1:n_distributions, inner=n_samples)
    perm = Random.shuffle(1:length(binom_labels))
    shuffled_binom_matrix = binom_matrix[:, perm]
    shuffled_binom_labels = binom_labels[perm]

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

    # Distance clustering
    @testset "Distance Clustering (Bray-Curtis + Optimal Clusters)" begin
        binom_distance_matrix = Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_binom_matrix)
        binom_unsupervised_clustering_result = Mycelia.identify_optimal_number_of_clusters(binom_distance_matrix)
        @test binom_unsupervised_clustering_result.optimal_number_of_clusters == n_distributions
        remapped_pred_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, binom_unsupervised_clustering_result.assignments)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, remapped_pred_labels)
        @test evaluation_result.macro_f1 >= .95
        @test evaluation_result.macro_precision >= .95
        @test evaluation_result.macro_recall >= .95
        @test evaluation_result.accuracy >= .95
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
    end

    # Distance + PCoA + KMeans
    @testset "Distance + PCoA + KMeans" begin
        pcoa_binom_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_binom_matrix))
        pcoa_fit_binom_labels = Clustering.kmeans(pcoa_binom_result.coordinates, n_distributions).assignments
        pcoa_fit_binom_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, pcoa_fit_binom_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, pcoa_fit_binom_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(pcoa_binom_result.coordinates;
                       title="PCoA - Bray-Curtis Distance - Binomial Matrix",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binom_labels,
                       fit_labels=pcoa_fit_binom_labels)
        display(plt)
    end

    # Distance + PCoA + UMAP + KMeans
    @testset "Distance + PCoA + UMAP + KMeans" begin
        pcoa_binom_result = Mycelia.pcoa_from_dist(Mycelia.frequency_matrix_to_braycurtis_distance_matrix(shuffled_binom_matrix))
        pcoa_binom_umap_model = Mycelia.umap_embed(pcoa_binom_result.coordinates)
        @test size(pcoa_binom_umap_model.embedding) == (2, n_samples * n_distributions)
        pcoa_binom_umap_fit_labels = Clustering.kmeans(pcoa_binom_umap_model.embedding, n_distributions).assignments
        pcoa_binom_umap_fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, pcoa_binom_umap_fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, pcoa_binom_umap_fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(pcoa_binom_umap_model.embedding;
                       title="PCoA - Bray-Curtis Distance - Binomial Matrix - UMAP",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binom_labels,
                       fit_labels=pcoa_binom_umap_fit_labels)
        display(plt)
    end

    # PoissonPCA-EPCA + KMeans
    @testset "PoissonPCA-EPCA + KMeans" begin
        poisson_pca_result = Mycelia.poisson_pca_epca(shuffled_binom_matrix, k=5)
        fit_labels = Clustering.kmeans(poisson_pca_result.scores, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(poisson_pca_result.scores;
                       title="Poisson PCA-EPCA - Binomial Matrix",
                       xlabel="PC1",
                       ylabel="PC2",
                       true_labels=shuffled_binom_labels,
                       fit_labels=fit_labels)
        display(plt)
    end

    # PoissonPCA-EPCA + UMAP + KMeans
    @testset "PoissonPCA-EPCA + UMAP + KMeans" begin
        poisson_pca_result = Mycelia.poisson_pca_epca(shuffled_binom_matrix, k=5)
        umap_model = Mycelia.umap_embed(poisson_pca_result.scores)
        fit_labels = Clustering.kmeans(umap_model.embedding, n_distributions).assignments
        fit_labels, mapping = Mycelia.best_label_mapping(shuffled_binom_labels, fit_labels)
        evaluation_result = Mycelia.evaluate_classification(shuffled_binom_labels, fit_labels)
        @test evaluation_result.macro_f1 >= 2/3
        @test evaluation_result.macro_precision >= 2/3
        @test evaluation_result.macro_recall >= 2/3
        @test evaluation_result.accuracy >= 2/3
        display(evaluation_result.confusion_matrix_plot)
        display(evaluation_result.f1_plot)
        display(evaluation_result.precision_plot)
        display(evaluation_result.recall_plot)
        plt = Mycelia.plot_embeddings(umap_model.embedding;
                       title="UMAP - Poisson PCA-EPCA - Binomial Matrix",
                       xlabel="UMAP 1",
                       ylabel="UMAP 2",
                       true_labels=shuffled_binom_labels,
                       fit_labels=fit_labels)
        display(plt)
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

    @testset "Sanity Check - Continuous Bernoulli Matrix" begin
        summary = Mycelia.sanity_check_matrix(shuffled_contb_matrix)
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

    @testset "Sanity Check - Gamma Matrix" begin
        summary = Mycelia.sanity_check_matrix(shuffled_gamma_matrix)
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
end