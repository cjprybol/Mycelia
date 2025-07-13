"""
IO Submodule for Mycelia.jl

Handles file input/output operations for FASTA/FASTQ files and format conversions.
"""
module IO

import BioSequences
import CSV
import CodecZlib
import DataFrames
import FASTX
import HTTP
import OrderedCollections
import ProgressMeter
import SHA
import Statistics
import UUIDs
import uCSV

# Import parent module utilities that these functions depend on
# (This will be updated when we reorganize the main module)

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Core file I/O operations
export open_fastx, write_fasta, write_fastq

# Quality control and assessment
export assess_duplication_rates, fastx_stats
export qc_filter_short_reads_fastp, qc_filter_long_reads_filtlong, qc_filter_long_reads_fastplong
export trim_galore_paired, run_fastqc

# File utilities and conversion
export join_fastqs_with_uuid
export fastx2normalized_table, fastxs2normalized_tables
export write_fastas_from_normalized_fastx_tables

# Sequence manipulation
export equivalent_fasta_sequences, merge_fasta_files, deduplicate_fasta_file
export translate_nucleic_acid_fasta, sort_fastq
export subsample_reads_seqkit, samtools_index_fasta
export genbank_to_fasta

# Data structure conversions  
export fasta_to_table, fasta_table_to_fasta
export fastq_record, fastx_to_contig_lengths

# Analysis functions
export count_records, determine_read_lengths, fasta_genome_size
export q_value_to_error_rate, error_rate_to_q_value

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export detect_sequence_extension, _detect_sequence_extension

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Copy all function implementations from fastx.jl here...
# [Functions would be copied here - truncated for brevity]

end # module IO