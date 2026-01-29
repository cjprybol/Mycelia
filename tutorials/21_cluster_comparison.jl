# # Tutorial 21: Cluster Comparison Methods
#
# This tutorial demonstrates how to compare two different clusterings of the same
# data using various similarity metrics. This is essential for:
#
# - Evaluating clustering algorithm performance against ground truth
# - Comparing results from different clustering methods
# - Assessing clustering reproducibility across runs
# - Validating clustering consistency under parameter changes
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - How to use multiple cluster comparison metrics (ARI, NMI, V-measure, etc.)
# - When to use each metric and their properties
# - How to interpret metric values and their ranges
# - How to generate comprehensive comparison reports
# - Practical considerations for comparing real-world clusterings

# ## Setup
#
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/21_cluster_comparison.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import DataFrames
import Random
import Statistics

Random.seed!(42)

println("=" ^ 60)
println("Tutorial 21: Cluster Comparison Methods")
println("=" ^ 60)

# ## Part 1: Introduction to Cluster Comparison
#
# When comparing clusterings, we have two vectors of cluster assignments for the
# same set of samples. The key insight is that we care about the *partition structure*,
# not the actual label values. Two clusterings that group samples identically but use
# different labels (e.g., [1,1,2,2] vs [A,A,B,B]) should be considered identical.
#
# ## Available Metrics in Mycelia
#
# | Metric | Range | Properties |
# |--------|-------|------------|
# | Adjusted Rand Index (ARI) | -1 to 1 | Corrected for chance, gold standard |
# | Normalized Mutual Information (NMI) | 0 to 1 | Information-theoretic |
# | Adjusted Mutual Information (AMI) | ~0 to 1 | NMI corrected for chance |
# | V-measure | 0 to 1 | Homogeneity + completeness |
# | Fowlkes-Mallows Index (FMI) | 0 to 1 | Geometric mean of precision/recall |
# | Jaccard Index | 0 to 1 | Pair-based similarity |
# | Cluster Purity | 0 to 1 | Majority class proportion (asymmetric) |

println("\n" * "=" ^ 60)
println("Part 1: Introduction to Cluster Comparison")
println("=" ^ 60)

# ## Part 2: Creating Example Clusterings

println("\n" * "=" ^ 60)
println("Part 2: Creating Example Clusterings")
println("=" ^ 60)

# ### Perfect Agreement

## Two identical clusterings (possibly with different labels)
labels_true = [1, 1, 1, 2, 2, 2, 3, 3, 3]
labels_pred = [1, 1, 1, 2, 2, 2, 3, 3, 3]

println("\nExample 1: Perfect Agreement")
println("  True:      ", labels_true)
println("  Predicted: ", labels_pred)

# ### Same Structure, Different Labels

## Same partition structure but relabeled
labels_relabeled = [3, 3, 3, 1, 1, 1, 2, 2, 2]  ## Same structure, different labels

println("\nExample 2: Same Structure, Different Labels")
println("  True:      ", labels_true)
println("  Relabeled: ", labels_relabeled)

# ### Partial Agreement

## One sample misassigned
labels_partial = [1, 1, 2, 2, 2, 2, 3, 3, 3]

println("\nExample 3: Partial Agreement")
println("  True:    ", labels_true)
println("  Partial: ", labels_partial)

# ### Low Agreement

## Different structure
labels_different = [1, 2, 1, 2, 1, 2, 1, 2, 1]

println("\nExample 4: Low Agreement")
println("  True:      ", labels_true)
println("  Different: ", labels_different)

# ## Part 3: Computing Individual Metrics

println("\n" * "=" ^ 60)
println("Part 3: Computing Individual Metrics")
println("=" ^ 60)

# ### Adjusted Rand Index (ARI)
#
# The ARI is the most commonly used metric for cluster comparison. It measures
# agreement between two clusterings while correcting for chance:
# - 1.0 = perfect agreement
# - 0.0 = agreement no better than random
# - < 0 = agreement worse than random

println("\n--- Adjusted Rand Index (ARI) ---")

ari_perfect = Mycelia.adjusted_rand_index(labels_true, labels_pred)
ari_relabeled = Mycelia.adjusted_rand_index(labels_true, labels_relabeled)
ari_partial = Mycelia.adjusted_rand_index(labels_true, labels_partial)
ari_different = Mycelia.adjusted_rand_index(labels_true, labels_different)

println("Perfect agreement:    ARI = $(round(ari_perfect, digits=4))")
println("Relabeled (same):     ARI = $(round(ari_relabeled, digits=4))")
println("Partial agreement:    ARI = $(round(ari_partial, digits=4))")
println("Different structure:  ARI = $(round(ari_different, digits=4))")

# ### Normalized Mutual Information (NMI)
#
# NMI is an information-theoretic measure based on entropy:
# - 1.0 = perfect agreement (complete mutual information)
# - 0.0 = no mutual information

println("\n--- Normalized Mutual Information (NMI) ---")

nmi_perfect = Mycelia.normalized_mutual_information(labels_true, labels_pred)
nmi_relabeled = Mycelia.normalized_mutual_information(labels_true, labels_relabeled)
nmi_partial = Mycelia.normalized_mutual_information(labels_true, labels_partial)
nmi_different = Mycelia.normalized_mutual_information(labels_true, labels_different)

println("Perfect agreement:    NMI = $(round(nmi_perfect, digits=4))")
println("Relabeled (same):     NMI = $(round(nmi_relabeled, digits=4))")
println("Partial agreement:    NMI = $(round(nmi_partial, digits=4))")
println("Different structure:  NMI = $(round(nmi_different, digits=4))")

# ### V-Measure (Homogeneity and Completeness)
#
# V-measure is the harmonic mean of two complementary metrics:
# - Homogeneity: each predicted cluster contains only members of a single true class
# - Completeness: all members of a true class are assigned to the same predicted cluster

println("\n--- V-Measure ---")

vm_perfect = Mycelia.v_measure(labels_true, labels_pred)
vm_partial = Mycelia.v_measure(labels_true, labels_partial)

println("Perfect agreement:")
println("  Homogeneity:  $(round(vm_perfect.homogeneity, digits=4))")
println("  Completeness: $(round(vm_perfect.completeness, digits=4))")
println("  V-measure:    $(round(vm_perfect.v_measure, digits=4))")

println("\nPartial agreement:")
println("  Homogeneity:  $(round(vm_partial.homogeneity, digits=4))")
println("  Completeness: $(round(vm_partial.completeness, digits=4))")
println("  V-measure:    $(round(vm_partial.v_measure, digits=4))")

# ### Fowlkes-Mallows Index (FMI)
#
# FMI is the geometric mean of precision and recall computed from pair counting:
# - Precision: TP / (TP + FP)
# - Recall: TP / (TP + FN)
# - FMI = sqrt(Precision * Recall)

println("\n--- Fowlkes-Mallows Index (FMI) ---")

fmi_perfect = Mycelia.fowlkes_mallows_index(labels_true, labels_pred)
fmi_partial = Mycelia.fowlkes_mallows_index(labels_true, labels_partial)
fmi_different = Mycelia.fowlkes_mallows_index(labels_true, labels_different)

println("Perfect agreement:    FMI = $(round(fmi_perfect, digits=4))")
println("Partial agreement:    FMI = $(round(fmi_partial, digits=4))")
println("Different structure:  FMI = $(round(fmi_different, digits=4))")

# ### Jaccard Index
#
# The Jaccard Index for clusterings measures agreement based on pairs:
# J = TP / (TP + FP + FN)

println("\n--- Jaccard Index ---")

ji_perfect = Mycelia.jaccard_index(labels_true, labels_pred)
ji_partial = Mycelia.jaccard_index(labels_true, labels_partial)
ji_different = Mycelia.jaccard_index(labels_true, labels_different)

println("Perfect agreement:    Jaccard = $(round(ji_perfect, digits=4))")
println("Partial agreement:    Jaccard = $(round(ji_partial, digits=4))")
println("Different structure:  Jaccard = $(round(ji_different, digits=4))")

# ### Cluster Purity
#
# Purity measures how pure each cluster is with respect to a reference clustering.
# NOTE: Purity is asymmetric - the order of arguments matters!

println("\n--- Cluster Purity (Asymmetric) ---")

purity_result = Mycelia.cluster_purity(labels_true, labels_partial)

println("Purity (true -> partial): $(round(purity_result.overall_purity, digits=4))")
println("Per-cluster purities:")
for (cluster, purity) in purity_result.cluster_purities
    println("  Cluster $cluster: $(round(purity, digits=4))")
end

# ## Part 4: Contingency Matrix

println("\n" * "=" ^ 60)
println("Part 4: Contingency Matrix")
println("=" ^ 60)

# The contingency matrix shows the overlap between clusters from two clusterings.
# It's the foundation for computing many comparison metrics.

cont = Mycelia.contingency_matrix(labels_true, labels_partial)

println("\nContingency Matrix:")
println("  Rows (labels_true): ", cont.labels1)
println("  Cols (labels_partial): ", cont.labels2)
println("  Matrix:")
for i in 1:size(cont.matrix, 1)
    println("    ", cont.matrix[i, :])
end

# ## Part 5: Comprehensive Comparison Summary

println("\n" * "=" ^ 60)
println("Part 5: Comprehensive Comparison Summary")
println("=" ^ 60)

# Mycelia provides convenience functions to compute all metrics at once
# and generate formatted reports.

## Generate comprehensive summary
summary = Mycelia.clustering_comparison_summary(
    labels_true, labels_partial;
    name1 = "Ground Truth",
    name2 = "K-means Result"
)

println("\nSummary fields available:")
for field in keys(summary)
    println("  - $field")
end

# ### Print Formatted Report

Mycelia.print_clustering_comparison(
    labels_true, labels_partial;
    name1 = "Ground Truth",
    name2 = "K-means Result"
)

# ### Export to DataFrame

df = Mycelia.comparison_summary_to_dataframe(summary)
println("\nDataFrame export:")
println(df)

# ## Part 6: Comparing Multiple Clustering Methods

println("\n" * "=" ^ 60)
println("Part 6: Comparing Multiple Clustering Methods")
println("=" ^ 60)

# A common use case is comparing multiple clustering methods against a reference.

## Simulate different clustering "methods"
Random.seed!(42)
n_samples = 50

## Ground truth: 5 clusters of 10 samples each
ground_truth = repeat(1:5, inner=10)

## Method A: Almost perfect (1 error per cluster)
method_a = copy(ground_truth)
method_a[10] = 2  ## One swap
method_a[20] = 3

## Method B: Merges two clusters
method_b = copy(ground_truth)
method_b[method_b .== 2] .= 1  ## Merge cluster 2 into cluster 1

## Method C: Random-ish
method_c = rand(1:5, 50)

## Compare all methods to ground truth
methods = [
    ("Method A (minor errors)", method_a),
    ("Method B (merged clusters)", method_b),
    ("Method C (random)", method_c)
]

println("\nComparing clustering methods to ground truth:")
println("-" ^ 70)
println(rpad("Method", 30), rpad("ARI", 10), rpad("NMI", 10), rpad("V-measure", 10))
println("-" ^ 70)

for (name, labels) in methods
    ari = Mycelia.adjusted_rand_index(ground_truth, labels)
    nmi = Mycelia.normalized_mutual_information(ground_truth, labels)
    vm = Mycelia.v_measure(ground_truth, labels)
    println(
        rpad(name, 30),
        rpad(round(ari, digits=4), 10),
        rpad(round(nmi, digits=4), 10),
        rpad(round(vm.v_measure, digits=4), 10)
    )
end
println("-" ^ 70)

# ## Part 7: Practical Considerations

println("\n" * "=" ^ 60)
println("Part 7: Practical Considerations")
println("=" ^ 60)

# ### When to Use Each Metric
#
# | Metric | Best For |
# |--------|----------|
# | ARI | General-purpose, corrected for chance, gold standard |
# | NMI | When you want information-theoretic interpretation |
# | AMI | When NMI seems inflated (chance correction) |
# | V-measure | When homogeneity vs completeness tradeoff matters |
# | FMI | Alternative to ARI with geometric mean interpretation |
# | Jaccard | Simple, intuitive pair-based similarity |
# | Purity | One-directional quality assessment |

println("\nMetric Selection Guidelines:")
println("  - Use ARI as your primary metric (most widely accepted)")
println("  - Add NMI for information-theoretic perspective")
println("  - Use V-measure when you need homogeneity/completeness breakdown")
println("  - Purity is useful but be aware it's asymmetric")

# ### Interpretation Guidelines
#
# | Score Range | Interpretation |
# |-------------|----------------|
# | > 0.9 | Excellent agreement (nearly identical) |
# | 0.7 - 0.9 | Good agreement |
# | 0.5 - 0.7 | Moderate agreement |
# | 0.3 - 0.5 | Weak agreement |
# | < 0.3 | Poor agreement |

println("\n" * "=" ^ 60)
println("Tutorial 21 Complete!")
println("=" ^ 60)

# ## Summary
#
# This tutorial covered:
# 1. The concept of cluster comparison and why label values don't matter
# 2. Individual metrics: ARI, NMI, AMI, V-measure, FMI, Jaccard, Purity
# 3. Contingency matrices as the foundation for comparison
# 4. Comprehensive comparison reports using `clustering_comparison_summary`
# 5. Practical guidelines for metric selection and interpretation
#
# ## Key Functions Used
#
# - `Mycelia.contingency_matrix(labels1, labels2)` - Build contingency matrix
# - `Mycelia.adjusted_rand_index(labels1, labels2)` - Compute ARI
# - `Mycelia.normalized_mutual_information(labels1, labels2)` - Compute NMI
# - `Mycelia.adjusted_mutual_information(labels1, labels2)` - Compute AMI
# - `Mycelia.v_measure(labels_true, labels_pred)` - Compute V-measure
# - `Mycelia.fowlkes_mallows_index(labels1, labels2)` - Compute FMI
# - `Mycelia.jaccard_index(labels1, labels2)` - Compute Jaccard
# - `Mycelia.cluster_purity(labels_true, labels_pred)` - Compute purity
# - `Mycelia.clustering_comparison_summary(...)` - All metrics at once
# - `Mycelia.print_clustering_comparison(...)` - Formatted report
# - `Mycelia.comparison_summary_to_dataframe(...)` - Export to DataFrame
