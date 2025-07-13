"""
SequenceClustering Submodule for Mycelia.jl

Handles sequence clustering and similarity analysis.
"""
module SequenceClustering

import Clustering
import Distances
import DataFrames
import SparseArrays

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Sequence clustering
export cluster_sequences, hierarchical_clustering
export kmeans_clustering, dbscan_clustering

# Clustering analysis
export assess_clustering_quality, optimal_cluster_number
export cluster_validation, silhouette_analysis

# Distance and similarity
export sequence_distance_matrix, jaccard_similarity
export edit_distance, hamming_distance

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _validate_clustering_params, _compute_cluster_stats
export _optimize_cluster_parameters

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from clustering.jl to here

end # module SequenceClustering