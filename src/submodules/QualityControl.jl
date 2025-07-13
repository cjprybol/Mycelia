"""
QualityControl Submodule for Mycelia.jl

Handles quality assessment and benchmarking for biological data.
"""
module QualityControl

import BioSequences
import DataFrames
import Statistics
import CSV
import ProgressMeter

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Quality assessment
export assess_read_quality, quality_metrics
export contamination_check, adapter_detection

# Benchmarking
export benchmark_assembly, benchmark_alignment
export performance_metrics, resource_usage

# Quality filtering
export filter_low_quality_reads, trim_adapters
export remove_duplicates, normalize_coverage

# Validation
export validate_file_format, check_data_integrity
export cross_validation, statistical_validation

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _calculate_quality_scores, _generate_quality_report
export _setup_qc_environment

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from quality-control-and-benchmarking.jl to here

end # module QualityControl