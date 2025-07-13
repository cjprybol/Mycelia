"""
Annotation Submodule for Mycelia.jl

Handles gene prediction and genome annotation.
"""
module Annotation

import GenomicAnnotations
import GFF3
import DataFrames
import BioSequences

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Gene prediction
export run_prodigal, run_augustus, run_genemark
export predict_genes, annotate_genes

# Codon optimization
export codon_optimize, calculate_codon_usage
export reverse_translate, optimize_for_expression

# Annotation analysis
export parse_gff, extract_gene_sequences
export annotation_stats, validate_annotations

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _parse_gene_predictions, _validate_gff_format
export _setup_annotation_environment

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from annotation.jl and codon-optimization.jl to here

end # module Annotation