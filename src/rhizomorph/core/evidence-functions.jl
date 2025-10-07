"""
Evidence manipulation and query helper functions.

These functions provide efficient operations on the nested evidence structure:
Dict{String, Dict{String, Set{EvidenceEntry}}}

# Function Categories
1. Adding evidence to vertices/edges
2. Querying evidence by dataset/observation/position
3. Counting observations and evidence entries
4. Manipulating evidence (merging, filtering, flipping strand)
"""

# ============================================================================
# Adding Evidence
# ============================================================================

"""
    add_evidence!(vertex, dataset_id, observation_id, entry)

Add an evidence entry to a vertex.

Handles initialization of nested dictionary structure automatically.

# Examples
```julia
vertex = KmerVertexData(kmer)
add_evidence!(vertex, "dataset_01", "read_001", EvidenceEntry(5, Forward))
```
"""
function add_evidence!(
    vertex::Union{KmerVertexData, BioSequenceVertexData, StringVertexData},
    dataset_id::String,
    observation_id::String,
    entry::EvidenceEntry
)
    # Initialize dataset if needed
    if !haskey(vertex.evidence, dataset_id)
        vertex.evidence[dataset_id] = Dict{String, Set{EvidenceEntry}}()
    end

    # Initialize observation if needed
    if !haskey(vertex.evidence[dataset_id], observation_id)
        vertex.evidence[dataset_id][observation_id] = Set{EvidenceEntry}()
    end

    # Add evidence entry
    push!(vertex.evidence[dataset_id][observation_id], entry)

    return vertex
end

"""
    add_evidence!(vertex, dataset_id, observation_id, entry::QualityEvidenceEntry)

Add a quality evidence entry to a quality-aware vertex.
"""
function add_evidence!(
    vertex::Union{QualmerVertexData, QualityBioSequenceVertexData, QualityStringVertexData},
    dataset_id::String,
    observation_id::String,
    entry::QualityEvidenceEntry
)
    # Initialize dataset if needed
    if !haskey(vertex.evidence, dataset_id)
        vertex.evidence[dataset_id] = Dict{String, Set{QualityEvidenceEntry}}()
    end

    # Initialize observation if needed
    if !haskey(vertex.evidence[dataset_id], observation_id)
        vertex.evidence[dataset_id][observation_id] = Set{QualityEvidenceEntry}()
    end

    # Add evidence entry
    push!(vertex.evidence[dataset_id][observation_id], entry)

    return vertex
end

"""
    add_evidence!(edge, dataset_id, observation_id, entry::EdgeEvidenceEntry)

Add an edge evidence entry to an edge.
"""
function add_evidence!(
    edge::Union{KmerEdgeData, BioSequenceEdgeData, StringEdgeData},
    dataset_id::String,
    observation_id::String,
    entry::EdgeEvidenceEntry
)
    # Initialize dataset if needed
    if !haskey(edge.evidence, dataset_id)
        edge.evidence[dataset_id] = Dict{String, Set{EdgeEvidenceEntry}}()
    end

    # Initialize observation if needed
    if !haskey(edge.evidence[dataset_id], observation_id)
        edge.evidence[dataset_id][observation_id] = Set{EdgeEvidenceEntry}()
    end

    # Add evidence entry
    push!(edge.evidence[dataset_id][observation_id], entry)

    return edge
end

"""
    add_evidence!(edge, dataset_id, observation_id, entry::EdgeQualityEvidenceEntry)

Add a quality edge evidence entry to a quality-aware edge.
"""
function add_evidence!(
    edge::Union{QualmerEdgeData, QualityBioSequenceEdgeData, QualityStringEdgeData},
    dataset_id::String,
    observation_id::String,
    entry::EdgeQualityEvidenceEntry
)
    # Initialize dataset if needed
    if !haskey(edge.evidence, dataset_id)
        edge.evidence[dataset_id] = Dict{String, Set{EdgeQualityEvidenceEntry}}()
    end

    # Initialize observation if needed
    if !haskey(edge.evidence[dataset_id], observation_id)
        edge.evidence[dataset_id][observation_id] = Set{EdgeQualityEvidenceEntry}()
    end

    # Add evidence entry
    push!(edge.evidence[dataset_id][observation_id], entry)

    return edge
end

# ============================================================================
# Querying Evidence - Dataset Level
# ============================================================================

"""
    get_dataset_evidence(vertex_or_edge, dataset_id)

Get all evidence from a specific dataset.

Returns Dict{String, Set{EvidenceEntry}} mapping observation_id to evidence,
or `nothing` if dataset doesn't exist.

# Examples
```julia
dataset_evidence = get_dataset_evidence(vertex, "dataset_01")
if !isnothing(dataset_evidence)
    for (obs_id, evidence_set) in dataset_evidence
        # Process evidence from this observation
    end
end
```
"""
function get_dataset_evidence(vertex_or_edge, dataset_id::String)
    if haskey(vertex_or_edge.evidence, dataset_id)
        return vertex_or_edge.evidence[dataset_id]
    else
        return nothing
    end
end

# ============================================================================
# Querying Evidence - Observation Level
# ============================================================================

"""
    get_observation_evidence(vertex_or_edge, dataset_id, observation_id)

Get all evidence from a specific observation within a dataset.

Returns Set{EvidenceEntry}, or `nothing` if observation doesn't exist.

# Examples
```julia
obs_evidence = get_observation_evidence(vertex, "dataset_01", "read_001")
if !isnothing(obs_evidence)
    for entry in obs_evidence
        println("Position: \$(entry.position), Strand: \$(entry.strand)")
    end
end
```
"""
function get_observation_evidence(vertex_or_edge, dataset_id::String, observation_id::String)
    dataset_evidence = get_dataset_evidence(vertex_or_edge, dataset_id)
    if isnothing(dataset_evidence)
        return nothing
    end

    if haskey(dataset_evidence, observation_id)
        return dataset_evidence[observation_id]
    else
        return nothing
    end
end

# ============================================================================
# Querying Evidence - Position Level
# ============================================================================

"""
    get_position_evidence(vertex, dataset_id, position)

Get all evidence at a specific position across all observations in a dataset.

Returns vector of evidence entries at that position.

# Examples
```julia
# Get all evidence at position 5
position_5_evidence = get_position_evidence(vertex, "dataset_01", 5)
```
"""
function get_position_evidence(
    vertex::Union{KmerVertexData, BioSequenceVertexData, StringVertexData},
    dataset_id::String,
    position::Int
)
    dataset_evidence = get_dataset_evidence(vertex, dataset_id)
    if isnothing(dataset_evidence)
        return EvidenceEntry[]
    end

    result = EvidenceEntry[]
    for (obs_id, evidence_set) in dataset_evidence
        for entry in evidence_set
            if entry.position == position
                push!(result, entry)
            end
        end
    end

    return result
end

# ============================================================================
# Evidence Manipulation
# ============================================================================

"""
    merge_evidence_sets(set1, set2)

Merge two evidence sets using set union.

# Examples
```julia
merged = merge_evidence_sets(evidence1, evidence2)
```
"""
function merge_evidence_sets(
    set1::Set{T},
    set2::Set{T}
) where {T<:Union{EvidenceEntry, QualityEvidenceEntry}}
    return union(set1, set2)
end

"""
    flip_evidence_strand(entry::EvidenceEntry)

Flip the strand orientation of an evidence entry.

Used for reverse complement operations.

# Examples
```julia
forward_entry = EvidenceEntry(5, Forward)
reverse_entry = flip_evidence_strand(forward_entry)  # EvidenceEntry(5, Reverse)
```
"""
function flip_evidence_strand(entry::EvidenceEntry)
    new_strand = entry.strand == Forward ? Reverse : Forward
    return EvidenceEntry(entry.position, new_strand)
end

"""
    flip_evidence_strand(entry::QualityEvidenceEntry)

Flip the strand orientation of a quality evidence entry.

Also reverses the quality scores (since RC reverses the sequence).

# Examples
```julia
entry = QualityEvidenceEntry(5, Forward, UInt8[30, 35, 40])
flipped = flip_evidence_strand(entry)
# Result: QualityEvidenceEntry(5, Reverse, UInt8[40, 35, 30])
```
"""
function flip_evidence_strand(entry::QualityEvidenceEntry)
    new_strand = entry.strand == Forward ? Reverse : Forward
    reversed_quality = reverse(entry.quality_scores)
    return QualityEvidenceEntry(entry.position, new_strand, reversed_quality)
end

# ============================================================================
# Counting Functions
# ============================================================================

"""
    count_total_observations(vertex_or_edge)

Count total number of unique observations across all datasets.

# Examples
```julia
total_obs = count_total_observations(vertex)  # e.g., 150
```
"""
function count_total_observations(vertex_or_edge)
    total = 0
    for (dataset_id, dataset_evidence) in vertex_or_edge.evidence
        total += length(keys(dataset_evidence))
    end
    return total
end

"""
    count_dataset_observations(vertex_or_edge, dataset_id)

Count number of observations in a specific dataset.

# Examples
```julia
dataset_obs = count_dataset_observations(vertex, "dataset_01")
```
"""
function count_dataset_observations(vertex_or_edge, dataset_id::String)
    dataset_evidence = get_dataset_evidence(vertex_or_edge, dataset_id)
    if isnothing(dataset_evidence)
        return 0
    end
    return length(keys(dataset_evidence))
end

"""
    count_evidence_entries(vertex_or_edge)

Count total number of evidence entries (not unique observations).

An observation may contribute multiple evidence entries if it visits
a vertex/edge multiple times.

# Examples
```julia
total_entries = count_evidence_entries(vertex)  # e.g., 500
```
"""
function count_evidence_entries(vertex_or_edge)
    total = 0
    for (dataset_id, dataset_evidence) in vertex_or_edge.evidence
        for (obs_id, evidence_set) in dataset_evidence
            total += length(evidence_set)
        end
    end
    return total
end

# ============================================================================
# Filtering Functions
# ============================================================================

"""
    filter_evidence_by_strand(vertex, strand)

Create a new vertex containing only evidence from specified strand.

# Examples
```julia
forward_only = filter_evidence_by_strand(vertex, Forward)
```
"""
function filter_evidence_by_strand(
    vertex::KmerVertexData{T},
    strand::StrandOrientation
) where {T}
    filtered = KmerVertexData(vertex.Kmer)

    for (dataset_id, dataset_evidence) in vertex.evidence
        for (obs_id, evidence_set) in dataset_evidence
            for entry in evidence_set
                if entry.strand == strand
                    add_evidence!(filtered, dataset_id, obs_id, entry)
                end
            end
        end
    end

    return filtered
end

# ============================================================================
# Query Functions - IDs
# ============================================================================

"""
    get_all_dataset_ids(vertex_or_edge)

Get list of all dataset IDs that have evidence.

# Examples
```julia
dataset_ids = get_all_dataset_ids(vertex)
# ["dataset_01", "dataset_02", "dataset_03"]
```
"""
function get_all_dataset_ids(vertex_or_edge)
    return collect(keys(vertex_or_edge.evidence))
end

"""
    get_all_observation_ids(vertex_or_edge, dataset_id)

Get list of all observation IDs in a dataset.

# Examples
```julia
obs_ids = get_all_observation_ids(vertex, "dataset_01")
# ["read_001", "read_002", "read_003"]
```
"""
function get_all_observation_ids(vertex_or_edge, dataset_id::String)
    dataset_evidence = get_dataset_evidence(vertex_or_edge, dataset_id)
    if isnothing(dataset_evidence)
        return String[]
    end
    return collect(keys(dataset_evidence))
end

# ============================================================================
# Weight/Coverage Computation (On-Demand, Not Stored)
# ============================================================================

"""
    compute_edge_weight(edge)

Compute edge weight from evidence.

Weight is based on number of observations (coverage).

# Examples
```julia
weight = compute_edge_weight(edge)
```
"""
function compute_edge_weight(edge)
    # Simple weight: count of unique observations
    return Float64(count_total_observations(edge))
end

"""
    compute_edge_coverage(edge)

Compute edge coverage (number of observations).

# Examples
```julia
coverage = compute_edge_coverage(edge)
```
"""
function compute_edge_coverage(edge)
    return count_total_observations(edge)
end
