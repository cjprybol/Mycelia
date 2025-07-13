"""
Sequences Submodule for Mycelia.jl

Handles sequence analysis including k-mer counting, quality-aware k-mers, and basic sequence operations.
"""
module Sequences

import DataFrames
import Dates
import JLD2
import Kmers
import ProgressMeter
import SparseArrays
import Statistics

# Import parent module utilities
# (This will be updated when we reorganize the main module)

# ============================================================================
# PUBLIC API EXPORTS  
# ============================================================================

# K-mer analysis
export count_kmers, count_canonical_kmers, kmer_profile, jaccard_similarity
export save_kmer_results, load_kmer_results

# Quality-aware k-mers (from qualmer-analysis.jl)
export qualmers_fw, qualmers_canonical, qualmers_unambiguous

# Sequence comparison and metrics
export sequence_similarity, sequence_distance

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _kmer_generator, _canonical_kmer

# ============================================================================
# TYPE DEFINITIONS
# ============================================================================

# From qualmer-analysis.jl
struct Qualmer
    sequence::String
    quality::Vector{Int}
end

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions from kmer-analysis.jl would be moved here
# Functions from qualmer-analysis.jl would be moved here
# [Truncated for brevity - full implementations would be copied]

end # module Sequences