"""
Variants Submodule for Mycelia.jl

Handles variant calling and analysis.
"""
module Variants

import BioSequences
import DataFrames
import XAM
import CSV

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Variant calling
export call_variants, run_bcftools, run_gatk
export identify_snps, identify_indels

# Variant analysis
export annotate_variants, filter_variants
export variant_stats, calculate_variant_frequency

# Population genetics
export calculate_fst, population_structure
export allele_frequency, hardy_weinberg_test

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _parse_vcf_output, _validate_variant_calls
export _setup_variant_calling_environment

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from variant-analysis.jl to here

end # module Variants