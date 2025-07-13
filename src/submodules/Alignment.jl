"""
Alignment Submodule for Mycelia.jl

Handles sequence alignment and mapping operations.
"""
module Alignment

import BioAlignments
import BioSequences
import CSV
import DataFrames
import XAM

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Sequence alignment
export minimap_index, minimap_map, run_clustal_omega
export bwa_index, bwa_mem, bowtie2_index, bowtie2_align

# Alignment analysis
export assess_alignment_accuracy, alignment_stats
export parse_sam_bam, extract_mapped_reads

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _parse_alignment_output, _validate_alignment_params
export _setup_alignment_environment

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from alignments-and-mapping.jl to here

end # module Alignment