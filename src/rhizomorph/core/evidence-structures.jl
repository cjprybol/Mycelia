"""
Evidence tracking structures for the Rhizomorph graph ecosystem.

These structures store observations of k-mers, sequences, and transitions
in a memory-efficient, queryable format.

# Key Design Principles
1. **Observation IDs as dictionary keys** - stored once, not repeated in every entry
2. **Nested structure** - dataset_id → observation_id → Set{Evidence}
3. **O(1) queries** - efficient lookup by dataset or observation
4. **Strand awareness** - all evidence tracks strand orientation
"""

"""
    EvidenceEntry

Basic evidence entry for non-quality-aware graphs.

Stores the position and strand orientation where a k-mer or sequence was observed.
The observation_id is stored as a dictionary key, not in the entry itself.

# Fields
- `position::Int`: Position in the observation where this evidence was recorded
- `strand::StrandOrientation`: Strand orientation (Forward or Reverse)

# Examples
```julia
entry = EvidenceEntry(42, Forward)
```
"""
struct EvidenceEntry
    position::Int
    strand::StrandOrientation
end

# Equality and hashing for Set operations
function Base.:(==)(a::EvidenceEntry, b::EvidenceEntry)
    a.position == b.position && a.strand == b.strand
end

Base.hash(e::EvidenceEntry, h::UInt) = hash(e.strand, hash(e.position, h))

"""
    QualityEvidenceEntry

Quality-aware evidence entry that includes per-base quality scores.

# Fields
- `position::Int`: Position in the observation
- `strand::StrandOrientation`: Strand orientation
- `quality_scores::Vector{UInt8}`: Phred quality scores (0-60 typical range)

# Details
Quality scores are stored as Vector{UInt8} to support variable-length sequences.
For fixed-length k-mers in QualmerVertexData, these vectors will all have the same
length (k), but the structure allows flexibility.

# Examples
```julia
entry = QualityEvidenceEntry(10, Forward, UInt8[30, 35, 32, 28, 40])
```
"""
struct QualityEvidenceEntry
    position::Int
    strand::StrandOrientation
    quality_scores::Vector{UInt8}
end

# Equality and hashing
function Base.:(==)(a::QualityEvidenceEntry, b::QualityEvidenceEntry)
    a.position == b.position &&
        a.strand == b.strand &&
        a.quality_scores == b.quality_scores
end

function Base.hash(e::QualityEvidenceEntry, h::UInt)
    hash(e.quality_scores, hash(e.strand, hash(e.position, h)))
end

"""
    EdgeEvidenceEntry

Evidence entry for edges (transitions between vertices).

# Fields
- `from_position::Int`: Position in observation where transition starts
- `to_position::Int`: Position in observation where transition ends
- `strand::StrandOrientation`: Strand orientation for this transition

# Details
For k-mer graphs, edges typically connect adjacent positions (to_position = from_position + 1).
For variable-length graphs, the gap may be larger depending on overlap structure.

# Examples
```julia
# Adjacent positions (typical for k-mer graphs)
entry = EdgeEvidenceEntry(5, 6, Forward)

# Non-adjacent (may occur in variable-length graphs)
entry = EdgeEvidenceEntry(100, 250, Reverse)
```
"""
struct EdgeEvidenceEntry
    from_position::Int
    to_position::Int
    strand::StrandOrientation
end

# Equality and hashing
function Base.:(==)(a::EdgeEvidenceEntry, b::EdgeEvidenceEntry)
    a.from_position == b.from_position &&
        a.to_position == b.to_position &&
        a.strand == b.strand
end

function Base.hash(e::EdgeEvidenceEntry, h::UInt)
    hash(e.strand, hash(e.to_position, hash(e.from_position, h)))
end

"""
    EdgeQualityEvidenceEntry

Quality-aware evidence entry for edges.

# Fields
- `from_position::Int`: Starting position
- `to_position::Int`: Ending position
- `strand::StrandOrientation`: Strand orientation
- `from_quality::Vector{UInt8}`: Quality scores at source vertex
- `to_quality::Vector{UInt8}`: Quality scores at target vertex

# Details
Stores quality information for both ends of the edge transition, enabling
quality-aware path finding and assembly algorithms.

# Examples
```julia
entry = EdgeQualityEvidenceEntry(
    5, 6, Forward,
    UInt8[30, 35, 32],  # Quality at source
    UInt8[28, 40, 38]   # Quality at target
)
```
"""
struct EdgeQualityEvidenceEntry
    from_position::Int
    to_position::Int
    strand::StrandOrientation
    from_quality::Vector{UInt8}
    to_quality::Vector{UInt8}
end

# Equality and hashing
function Base.:(==)(a::EdgeQualityEvidenceEntry, b::EdgeQualityEvidenceEntry)
    a.from_position == b.from_position &&
        a.to_position == b.to_position &&
        a.strand == b.strand &&
        a.from_quality == b.from_quality &&
        a.to_quality == b.to_quality
end

function Base.hash(e::EdgeQualityEvidenceEntry, h::UInt)
    hash(e.to_quality, hash(e.from_quality,
        hash(e.strand, hash(e.to_position, hash(e.from_position, h)))))
end
