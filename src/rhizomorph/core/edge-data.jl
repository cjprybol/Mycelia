"""
Edge data structures for the Rhizomorph graph ecosystem.

These structures store transitions between vertices and their associated evidence.

# Key Requirements
1. Use nested Dict evidence structure (not flat vectors)
2. NO premature weight calculation - compute on-demand
3. Strand information in evidence entries, not hardcoded in edge
"""

"""
    KmerEdgeData

Edge data for k-mer graphs.

# Fields
- `evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}`: Edge evidence
  - Outer key: dataset_id
  - Inner key: observation_id
  - Value: Set of EdgeEvidenceEntry (from_pos, to_pos, strand)

# Design Notes
- **NO stored weight** - compute on-demand from evidence
- **NO stored strand orientations** - different observations may traverse edge in different orientations
- Evidence tracks all transitions, allowing flexible weight/strand calculations

# Examples
```julia
edge = KmerEdgeData()

add_evidence!(edge, "dataset_01", "read_001",
             EdgeEvidenceEntry(5, 6, Forward))

# Compute weight on-demand
weight = compute_edge_weight(edge)
coverage = compute_edge_coverage(edge)
```
"""
struct KmerEdgeData
    evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}

    function KmerEdgeData()
        evidence = Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}()
        new(evidence)
    end
end

"""
    QualmerEdgeData

Edge data for quality-aware k-mer graphs.

# Fields
- `evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}`: Quality edge evidence

# Examples
```julia
edge = QualmerEdgeData()

add_evidence!(edge, "dataset_01", "read_001",
             EdgeQualityEvidenceEntry(5, 6, Forward, UInt8[30, 35], UInt8[32, 28]))
```
"""
struct QualmerEdgeData
    evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}

    function QualmerEdgeData()
        evidence = Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}()
        new(evidence)
    end
end

"""
    BioSequenceEdgeData

Edge data for variable-length BioSequence graphs.

# Fields
- `overlap_length::Int`: Length of overlap between sequences (structural property)
- `evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}`: Edge evidence

# Details
For variable-length graphs, overlap_length is a structural property determined
by sequence alignment, so it's stored directly rather than computed from evidence.

# Examples
```julia
edge = BioSequenceEdgeData(50)  # 50bp overlap

add_evidence!(edge, "dataset_01", "read_001",
             EdgeEvidenceEntry(1, 150, Forward))
```
"""
struct BioSequenceEdgeData
    overlap_length::Int
    evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}

    function BioSequenceEdgeData(overlap_length::Int)
        evidence = Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}()
        new(overlap_length, evidence)
    end
end

"""
    QualityBioSequenceEdgeData

Edge data for variable-length quality-aware BioSequence graphs.

# Fields
- `overlap_length::Int`: Length of overlap between sequences
- `evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}`: Quality evidence

# Examples
```julia
edge = QualityBioSequenceEdgeData(50)

add_evidence!(edge, "dataset_01", "read_001",
             EdgeQualityEvidenceEntry(1, 150, Forward,
                                     UInt8[30, 35, 32],  # from quality
                                     UInt8[28, 40, 38])) # to quality
```
"""
struct QualityBioSequenceEdgeData
    overlap_length::Int
    evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}

    function QualityBioSequenceEdgeData(overlap_length::Int)
        evidence = Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}()
        new(overlap_length, evidence)
    end
end

"""
    StringEdgeData

Edge data for string graphs.

# Fields
- `overlap_length::Int`: Length of overlap between strings
- `evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}`: Edge evidence

# Examples
```julia
edge = StringEdgeData(2)  # 2-character overlap

add_evidence!(edge, "dataset_01", "text_001",
             EdgeEvidenceEntry(5, 6, Forward))
```
"""
struct StringEdgeData
    overlap_length::Int
    evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}

    function StringEdgeData(overlap_length::Int)
        evidence = Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}()
        new(overlap_length, evidence)
    end
end

"""
    QualityStringEdgeData

Edge data for quality-aware string graphs.

# Fields
- `overlap_length::Int`: Length of overlap between strings
- `evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}`: Quality evidence

# Examples
```julia
edge = QualityStringEdgeData(2)

add_evidence!(edge, "dataset_01", "text_001",
             EdgeQualityEvidenceEntry(5, 6, Forward,
                                     UInt8[30, 35], UInt8[32, 28]))
```
"""
struct QualityStringEdgeData
    overlap_length::Int
    evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}

    function QualityStringEdgeData(overlap_length::Int)
        evidence = Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}()
        new(overlap_length, evidence)
    end
end

# ============================================================================
# Lightweight Edge Data Types
#
# These types store only aggregated counts instead of full evidence sets.
# They pair with the lightweight vertex types to enable low-memory graph
# construction for analysis workflows that only need topology + counts.
# ============================================================================

"""
    LightweightEdgeData

Lightweight edge data that stores only aggregated counts and unique observation
tracking.

Used for both k-mer and n-gram/string edges in lightweight mode. For n-gram
edges the overlap_length is structural and always stored.

# Fields
- `overlap_length::Int`: Length of overlap (0 for k-mer edges, n-1 for n-gram edges)
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `dataset_observations::Dict{String, Set{String}}`: Per-dataset unique
  observation IDs

# Semantic Equivalence with Full Mode
- `count_total_observations` → unique observation IDs (matches full mode)
- `count_evidence_entries` → total entries (matches full mode)
"""
mutable struct LightweightEdgeData
    const overlap_length::Int
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const dataset_observations::Dict{String, Set{String}}

    function LightweightEdgeData(overlap_length::Int = 0)
        new(overlap_length, 0, Dict{String, Int}(), Dict{String, Set{String}}())
    end
end

# ============================================================================
# Ultralight Edge Data Types
#
# These types store only aggregated counts — NO observation ID tracking.
# They pair with the ultralight vertex types.
# ============================================================================

"""
    UltralightEdgeData

Ultralight edge data that stores only aggregated counts. No observation ID
tracking (unlike LightweightEdgeData).

Used for both k-mer and n-gram/string edges in ultralight mode.

# Fields
- `overlap_length::Int`: Length of overlap (0 for k-mer edges, n-1 for n-gram edges)
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts

# Semantic Difference from Lightweight
- `count_total_observations` → returns `total_count` (NOT unique obs count)
- `get_all_observation_ids` → returns `nothing`
"""
mutable struct UltralightEdgeData
    const overlap_length::Int
    total_count::Int
    const dataset_counts::Dict{String, Int}

    function UltralightEdgeData(overlap_length::Int = 0)
        new(overlap_length, 0, Dict{String, Int}())
    end
end

# ============================================================================
# Ultralight Quality Edge Data Types
#
# These types store aggregated counts AND joint quality scores for the overlap
# region, but NO observation ID tracking.
# ============================================================================

"""
    UltralightQualityEdgeData

Ultralight quality-aware edge data. Stores aggregated counts and joint quality
scores for overlap regions, but no observation ID tracking.

# Fields
- `overlap_length::Int`: Length of overlap (0 for k-mer edges)
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `from_joint_quality::Vector{UInt8}`: Joint Phred scores for source k-mer
- `to_joint_quality::Vector{UInt8}`: Joint Phred scores for destination k-mer
"""
mutable struct UltralightQualityEdgeData
    const overlap_length::Int
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const from_joint_quality::Vector{UInt8}
    const to_joint_quality::Vector{UInt8}

    function UltralightQualityEdgeData(overlap_length::Int = 0)
        new(overlap_length, 0, Dict{String, Int}(), Vector{UInt8}(), Vector{UInt8}())
    end
end

# ============================================================================
# Lightweight Quality Edge Data Types
#
# These types store aggregated counts, observation ID tracking, AND joint
# quality scores.
# ============================================================================

"""
    LightweightQualityEdgeData

Lightweight quality-aware edge data. Combines observation tracking with
quality aggregation.

# Fields
- `overlap_length::Int`: Length of overlap (0 for k-mer edges)
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `dataset_observations::Dict{String, Set{String}}`: Per-dataset unique observation IDs
- `from_joint_quality::Vector{UInt8}`: Joint Phred scores for source k-mer
- `to_joint_quality::Vector{UInt8}`: Joint Phred scores for destination k-mer
"""
mutable struct LightweightQualityEdgeData
    const overlap_length::Int
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const dataset_observations::Dict{String, Set{String}}
    const from_joint_quality::Vector{UInt8}
    const to_joint_quality::Vector{UInt8}

    function LightweightQualityEdgeData(overlap_length::Int = 0)
        new(overlap_length, 0, Dict{String, Int}(), Dict{String, Set{String}}(),
            Vector{UInt8}(), Vector{UInt8}())
    end
end
