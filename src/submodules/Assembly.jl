"""
Assembly Submodule for Mycelia.jl

Handles genome assembly tools and quality assessment.
"""
module Assembly

import BioSequences
import CSV
import DataFrames
import Dates
import DelimitedFiles
import JLD2
import ProgressMeter
import Statistics

# Import parent module utilities that these functions depend on
# (This will be updated when we reorganize the main module)

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Core assembly functions
export run_megahit, run_metaspades, run_flye, run_canu, run_unicycler
export run_hifiasm, run_nextdenovo, run_raven

# Assembly quality assessment
export assess_assembly_quality, assess_assembly_completeness_busco
export contig_stats, assembly_stats

# Assembly utilities
export extract_contigs_by_length, merge_assemblies

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _setup_assembly_environment, _parse_assembly_stats
export _validate_assembly_params, _check_assembly_dependencies

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from assembly.jl to here
# Placeholder - actual implementations would be copied from assembly.jl

end # module Assembly