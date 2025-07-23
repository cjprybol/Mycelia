"""
Strain-Resolved Assembly Types - Enhanced from Mycelia-Dev

This module provides comprehensive data structures for strain-resolved assembly
with detailed quality metrics, confidence scoring, and statistical validation.

Adapted from Mycelia-Dev/src/core/types.jl with improvements for the current
Mycelia architecture.
"""

# =============================================================================
# Core Enumerations
# =============================================================================

"""
    PolymerType

Enumeration of supported biological polymer types for assembly
"""
@enum PolymerType DNA=1 RNA=2 AA=3

# =============================================================================
# Comprehensive Quality Metrics
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Comprehensive quality metrics for strain-resolved assembly evaluation.

# Fields
- `strain_recall::Float64`: True strains recovered / Total true strains
- `strain_precision::Float64`: Correct strain calls / Total strain calls  
- `strain_f1_score::Float64`: Harmonic mean of strain precision and recall
- `NGA50::Int64`: N50 based on aligned contigs
- `largest_alignment::Int64`: Longest single alignment length
- `total_aligned_length::Int64`: Total aligned sequence length
- `mismatches_per_100kb::Float64`: SNV error rate per 100kb
- `indels_per_100kb::Float64`: Indel error rate per 100kb
- `misassemblies_total::Int64`: Total misassemblies detected
- `misassemblies_local::Int64`: Local misassemblies only
- `completeness::Float64`: Fraction of reference genome covered
- `contamination::Float64`: Fraction of assembly that's contaminant
- `busco_completeness::Float64`: BUSCO gene completeness score
- `read_mapping_rate::Float64`: Fraction of reads that map back to assembly
- `properly_paired_rate::Float64`: Fraction of read pairs properly mapped
- `mean_base_confidence::Float64`: Average confidence score per base
- `uncertain_regions_fraction::Float64`: Fraction of assembly with low confidence
- `peak_memory_gb::Float64`: Maximum memory usage during assembly
- `cpu_time_hours::Float64`: Total computation time in hours
- `homopolymer_accuracy::Float64`: Accuracy in homopolymer regions
- `repeat_resolution_rate::Float64`: Fraction of repeats properly resolved

# Example
```julia
metrics = StrainQualityMetrics()
metrics.strain_recall = 0.95
metrics.NGA50 = 50000
```
"""
struct StrainQualityMetrics
    # Strain Resolution Metrics
    strain_recall::Float64           
    strain_precision::Float64        
    strain_f1_score::Float64         
    
    # Assembly Contiguity Metrics
    NGA50::Int64                     
    largest_alignment::Int64         
    total_aligned_length::Int64      
    
    # Assembly Accuracy Metrics
    mismatches_per_100kb::Float64    
    indels_per_100kb::Float64        
    misassemblies_total::Int64       
    misassemblies_local::Int64       
    
    # Completeness and Contamination
    completeness::Float64            
    contamination::Float64           
    busco_completeness::Float64      
    
    # Read Mapping Metrics
    read_mapping_rate::Float64       
    properly_paired_rate::Float64    
    
    # Confidence and Uncertainty
    mean_base_confidence::Float64    
    uncertain_regions_fraction::Float64  
    
    # Resource Usage
    peak_memory_gb::Float64          
    cpu_time_hours::Float64          
    
    # Technology-Specific Metrics
    homopolymer_accuracy::Float64    
    repeat_resolution_rate::Float64  
    
    # Default constructor with zeros
    StrainQualityMetrics() = new(
        0.0, 0.0, 0.0,  # strain metrics
        0, 0, 0,        # contiguity metrics
        0.0, 0.0, 0, 0, # accuracy metrics
        0.0, 0.0, 0.0,  # completeness metrics
        0.0, 0.0,       # mapping metrics
        0.0, 0.0,       # confidence metrics
        0.0, 0.0,       # resource metrics
        0.0, 0.0        # technology metrics
    )
end

# =============================================================================
# Strain-Resolved Contig Representation
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Represents an assembled contig with strain-specific metadata and confidence scoring.

# Fields
- `id::String`: Unique contig identifier
- `sequence::BioSequences.LongDNA{4}`: Assembled DNA sequence
- `length::Int64`: Sequence length in base pairs
- `coverage::Float64`: Mean coverage depth across contig
- `confidence::Vector{Float64}`: Per-base confidence scores (0.0-1.0)
- `strain_id::Union{String, Nothing}`: Associated strain identifier
- `strain_probability::Float64`: Probability of strain assignment

# Example
```julia
contig = StrainContig("contig_001", seq)
contig.strain_id = "strain_A"
contig.strain_probability = 0.95
```
"""
struct StrainContig
    id::String
    sequence::BioSequences.LongDNA{4}
    length::Int64
    coverage::Float64
    confidence::Vector{Float64}  # Per-base confidence scores
    strain_id::Union{String, Nothing}
    strain_probability::Float64  # Confidence in strain assignment
    
    function StrainContig(id::String, sequence::BioSequences.LongDNA{4})
        length = Base.length(sequence)
        new(id, sequence, length, 0.0, zeros(length), nothing, 0.0)
    end
    
    function StrainContig(id::String, sequence::BioSequences.LongDNA{4}, 
                         strain_id::String, strain_probability::Float64)
        length = Base.length(sequence)
        new(id, sequence, length, 0.0, zeros(length), strain_id, strain_probability)
    end
end

# =============================================================================
# Genetic Variant Representation
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Represents a detected genetic variant with strain-specific context.

# Fields
- `contig_id::String`: Reference contig identifier
- `position::Int64`: 1-based position in contig
- `reference::String`: Reference sequence at position
- `alternative::String`: Alternative sequence
- `variant_type::Symbol`: Variant classification (:SNV, :INDEL, :SV)
- `confidence::Float64`: Confidence in variant call (0.0-1.0)
- `allele_frequency::Float64`: Frequency in strain population
- `supporting_reads::Int64`: Number of supporting reads
- `strain_id::Union{String, Nothing}`: Associated strain if known

# Example
```julia
variant = StrainVariant("contig_001", 1000, "A", "G", :SNV, 0.95, 0.25, 15)
```
"""
struct StrainVariant
    contig_id::String
    position::Int64
    reference::String
    alternative::String
    variant_type::Symbol  # :SNV, :INDEL, :SV
    confidence::Float64
    allele_frequency::Float64
    supporting_reads::Int64
    strain_id::Union{String, Nothing}
    
    function StrainVariant(contig_id::String, position::Int64, reference::String, 
                          alternative::String, variant_type::Symbol, confidence::Float64,
                          allele_frequency::Float64, supporting_reads::Int64)
        new(contig_id, position, reference, alternative, variant_type, 
            confidence, allele_frequency, supporting_reads, nothing)
    end
end

# =============================================================================
# Strain-Resolved Assembly Graph
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Strain-resolved overlap graph for multi-strain assembly with confidence scoring.

# Fields
- `nodes::Dict{T, Int64}`: k-mer nodes with occurrence counts
- `edges::Dict{Tuple{T, T}, Float64}`: Weighted edges between k-mers
- `strain_paths::Dict{String, Vector{T}}`: Strain-specific paths through graph
- `confidence_scores::Dict{T, Float64}`: Per-node confidence scores
- `strain_support::Dict{T, Dict{String, Int64}}`: Strain-specific node support

# Example
```julia
graph = StrainAssemblyGraph{BioSequences.LongDNA{4}}()
graph.strain_paths["strain_A"] = [kmer1, kmer2, kmer3]
```
"""
mutable struct StrainAssemblyGraph{T <: BioSequences.BioSequence}
    nodes::Dict{T, Int64}                    # k-mer -> count
    edges::Dict{Tuple{T, T}, Float64}        # (k-mer1, k-mer2) -> weight
    strain_paths::Dict{String, Vector{T}}    # strain_id -> path through graph
    confidence_scores::Dict{T, Float64}      # k-mer -> confidence
    strain_support::Dict{T, Dict{String, Int64}}  # k-mer -> strain-specific counts
    
    function StrainAssemblyGraph{T}() where {T <: BioSequences.BioSequence}
        new{T}(
            Dict{T, Int64}(),
            Dict{Tuple{T, T}, Float64}(),
            Dict{String, Vector{T}}(),
            Dict{T, Float64}(),
            Dict{T, Dict{String, Int64}}()
        )
    end
end

# =============================================================================
# Complete Assembly Result
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Complete result of strain-resolved assembly process with comprehensive metadata.

# Fields
- `contigs::Vector{StrainContig}`: Assembled contigs with strain assignments
- `variants::Vector{StrainVariant}`: Detected genetic variants
- `quality_metrics::StrainQualityMetrics`: Comprehensive quality assessment
- `assembly_graph::Union{StrainAssemblyGraph, Nothing}`: Optional assembly graph
- `strain_abundances::Dict{String, Float64}`: Estimated strain relative abundances
- `metadata::Dict{String, Any}`: Additional assembly metadata

# Example
```julia
result = StrainAssemblyResult()
result.strain_abundances["strain_A"] = 0.6
result.strain_abundances["strain_B"] = 0.4
```
"""
struct StrainAssemblyResult
    contigs::Vector{StrainContig}
    variants::Vector{StrainVariant}
    quality_metrics::StrainQualityMetrics
    assembly_graph::Union{StrainAssemblyGraph, Nothing}
    strain_abundances::Dict{String, Float64}  # strain_id -> relative abundance
    metadata::Dict{String, Any}
    
    StrainAssemblyResult() = new(
        StrainContig[], 
        StrainVariant[], 
        StrainQualityMetrics(), 
        nothing,
        Dict{String, Float64}(),
        Dict{String, Any}()
    )
end

# =============================================================================
# Statistical Validation Types
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Results from strain-aware cross-validation analysis.

# Fields
- `mean_strain_f1::Float64`: Mean F1 score across strains
- `std_strain_f1::Float64`: Standard deviation of strain F1 scores
- `mean_nga50::Float64`: Mean NGA50 across validation folds
- `strain_consistency::Float64`: Consistency of strain assignments across folds

# Example
```julia
cv_results = StrainCVResults(0.92, 0.05, 45000.0, 0.88)
```
"""
struct StrainCVResults
    mean_strain_f1::Float64
    std_strain_f1::Float64
    mean_nga50::Float64
    strain_consistency::Float64  # Consistency of strain assignments
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Results from bootstrap resampling analysis for strain-resolved assembly.

# Fields
- `mean_performance::Float64`: Mean performance across bootstrap samples
- `confidence_interval_95::Tuple{Float64, Float64}`: 95% confidence interval
- `strain_stability::Float64`: Stability of strain assignments across samples

# Example
```julia
bootstrap_results = StrainBootstrapResults(0.91, (0.87, 0.95), 0.89)
```
"""
struct StrainBootstrapResults
    mean_performance::Float64
    confidence_interval_95::Tuple{Float64, Float64}
    strain_stability::Float64  # Stability of strain detection across samples
end

# =============================================================================
# Benchmarking Dataset Types
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Represents a strain-resolved benchmarking dataset with expected outcomes.

# Fields
- `name::String`: Dataset identifier
- `description::String`: Dataset description
- `file_path::String`: Path to data files
- `expected_metrics::StrainQualityMetrics`: Expected assembly quality
- `known_strains::Vector{String}`: Known strain identifiers
- `strain_abundances::Dict{String, Float64}`: True strain abundances
- `difficulty_level::Symbol`: Difficulty classification (:easy, :medium, :hard, :extreme)

# Example
```julia
dataset = StrainBenchmarkDataset(
    "mixed_ecoli", 
    "E. coli strain mixture", 
    "/data/mixed_ecoli.fastq",
    expected_metrics,
    ["strain_K12", "strain_O157"],
    Dict("strain_K12" => 0.7, "strain_O157" => 0.3),
    :medium
)
```
"""
struct StrainBenchmarkDataset
    name::String
    description::String
    file_path::String
    expected_metrics::StrainQualityMetrics
    known_strains::Vector{String}
    strain_abundances::Dict{String, Float64}  # Expected strain abundances
    difficulty_level::Symbol  # :easy, :medium, :hard, :extreme
end

# =============================================================================
# Helper Functions
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Calculate F1 score from precision and recall values.

# Arguments
- `precision::Float64`: Precision value (0.0-1.0)
- `recall::Float64`: Recall value (0.0-1.0)

# Returns
- `Float64`: F1 score (harmonic mean of precision and recall)

# Example
```julia
f1 = calculate_f1_score(0.90, 0.85)  # Returns ~0.87
```
"""
function calculate_f1_score(precision::Float64, recall::Float64)::Float64
    if precision + recall == 0.0
        return 0.0
    end
    return 2.0 * (precision * recall) / (precision + recall)
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Update strain quality metrics with calculated F1 score.

# Arguments
- `metrics::StrainQualityMetrics`: Quality metrics struct to update

# Returns
- `StrainQualityMetrics`: Updated metrics with F1 score calculated

# Example
```julia
metrics = update_strain_f1_score!(metrics)
```
"""
function update_strain_f1_score!(metrics::StrainQualityMetrics)
    metrics = StrainQualityMetrics(
        metrics.strain_recall,
        metrics.strain_precision,
        calculate_f1_score(metrics.strain_precision, metrics.strain_recall),
        metrics.NGA50,
        metrics.largest_alignment,
        metrics.total_aligned_length,
        metrics.mismatches_per_100kb,
        metrics.indels_per_100kb,
        metrics.misassemblies_total,
        metrics.misassemblies_local,
        metrics.completeness,
        metrics.contamination,
        metrics.busco_completeness,
        metrics.read_mapping_rate,
        metrics.properly_paired_rate,
        metrics.mean_base_confidence,
        metrics.uncertain_regions_fraction,
        metrics.peak_memory_gb,
        metrics.cpu_time_hours,
        metrics.homopolymer_accuracy,
        metrics.repeat_resolution_rate
    )
    return metrics
end