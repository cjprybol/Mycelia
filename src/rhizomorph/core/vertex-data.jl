"""
Vertex data structures for the Rhizomorph graph ecosystem.

These structures store k-mers, sequences, and their associated evidence
in a type-stable, memory-efficient format.

# Key Requirements
1. Store OBSERVED k-mer/sequence (NOT canonical)
2. Use nested Dict evidence structure for O(1) queries
3. Generic parameterization for type stability
"""

"""
    KmerVertexData{T}

Vertex data for fixed-length k-mer graphs (DNA, RNA, or amino acid).

# Type Parameter
- `T`: K-mer type (e.g., DNAKmer{31}, RNAKmer{21}, AAKmer{11})

# Fields
- `Kmer::T`: The k-mer AS OBSERVED (not canonical)
- `evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}`: Evidence structure
  - Outer key: dataset_id
  - Inner key: observation_id (e.g., read identifier)
  - Value: Set of EvidenceEntry (position, strand)

# Critical Design Notes
- **NEVER stores canonical k-mer** - stores exactly what was observed
- **Strand-specific by default** - Forward and Reverse k-mers are separate vertices
- **Efficient queries** - O(1) lookup by dataset or observation
- **Memory efficient** - observation_id stored once as key, not repeated

# Examples
```julia
# Create vertex for observed k-mer (NOT canonical)
kmer = DNAKmer{3}("ATG")
vertex = KmerVertexData(kmer)

# Add evidence
add_evidence!(vertex, "dataset_01", "read_001", EvidenceEntry(5, Forward))

# Query by dataset
dataset_evidence = get_dataset_evidence(vertex, "dataset_01")

# Query by observation
obs_evidence = get_observation_evidence(vertex, "dataset_01", "read_001")
```
"""
struct KmerVertexData{T}
    Kmer::T
    evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}

    function KmerVertexData(kmer::T) where {T}
        evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()
        new{T}(kmer, evidence)
    end
end

"""
    QualmerVertexData{T}

Vertex data for quality-aware k-mer graphs.

# Type Parameter
- `T`: K-mer type (DNA or RNA k-mers only - amino acids don't have quality)

# Fields
- `Kmer::T`: The k-mer AS OBSERVED (not canonical)
- `evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}`: Quality evidence
  - Outer key: dataset_id
  - Inner key: observation_id
  - Value: Set of QualityEvidenceEntry (position, strand, quality_scores)

# Examples
```julia
kmer = DNAKmer{3}("ATG")
vertex = QualmerVertexData(kmer)

add_evidence!(vertex, "dataset_01", "read_001",
             QualityEvidenceEntry(5, Forward, UInt8[30, 35, 32]))
```
"""
struct QualmerVertexData{T}
    Kmer::T
    evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}

    function QualmerVertexData(kmer::T) where {T}
        evidence = Dict{String, Dict{String, Set{QualityEvidenceEntry}}}()
        new{T}(kmer, evidence)
    end
end

"""
    BioSequenceVertexData{T}

Vertex data for variable-length BioSequence graphs (FASTA graphs).

# Type Parameter
- `T`: BioSequence type (LongDNA{4}, LongRNA{4}, LongAA)

# Fields
- `sequence::T`: The sequence AS OBSERVED
- `evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}`: Evidence structure

# Examples
```julia
sequence = LongDNA{4}("ATGATGATG")
vertex = BioSequenceVertexData(sequence)

add_evidence!(vertex, "dataset_01", "read_001", EvidenceEntry(1, Forward))
```
"""
struct BioSequenceVertexData{T}
    sequence::T
    evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}

    function BioSequenceVertexData(sequence::T) where {T}
        evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()
        new{T}(sequence, evidence)
    end
end

"""
    QualityBioSequenceVertexData{T}

Vertex data for variable-length quality-aware BioSequence graphs (FASTQ graphs).

# Type Parameter
- `T`: BioSequence type (LongDNA{4}, LongRNA{4})

# Fields
- `sequence::T`: The sequence AS OBSERVED
- `quality_scores::Vector{UInt8}`: Per-base quality scores (Phred+33 encoded)
- `evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}`: Quality evidence

# Examples
```julia
sequence = LongDNA{4}("ATGATGATG")
vertex = QualityBioSequenceVertexData(sequence)

add_evidence!(vertex, "dataset_01", "read_001",
             QualityEvidenceEntry(1, Forward, UInt8[30, 35, 32, 28, 40, 38, 35, 32, 30]))
```
"""
struct QualityBioSequenceVertexData{T}
    sequence::T
    quality_scores::Vector{UInt8}
    evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}

    function QualityBioSequenceVertexData(
            sequence::T, quality_scores::Vector{UInt8} = Vector{UInt8}()) where {T}
        if !isempty(quality_scores) && length(sequence) != length(quality_scores)
            error("Sequence and quality lengths must match")
        end
        evidence = Dict{String, Dict{String, Set{QualityEvidenceEntry}}}()
        new{T}(sequence, quality_scores, evidence)
    end
end

"""
    StringVertexData

Vertex data for string graphs (n-gram graphs or variable-length string graphs).

# Fields
- `string_value::String`: The string AS OBSERVED
- `evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}`: Evidence structure

# Examples
```julia
vertex = StringVertexData("ATG")

add_evidence!(vertex, "dataset_01", "text_001", EvidenceEntry(5, Forward))
```
"""
struct StringVertexData
    string_value::String
    evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}

    function StringVertexData(string_value::String)
        evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()
        new(string_value, evidence)
    end
end

"""
    QualityStringVertexData

Vertex data for quality-aware string graphs.

# Fields
- `string_value::String`: The string AS OBSERVED
- `evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}`: Quality evidence

# Examples
```julia
vertex = QualityStringVertexData("ATG")

add_evidence!(vertex, "dataset_01", "text_001",
             QualityEvidenceEntry(5, Forward, UInt8[30, 35, 32]))
```
"""
struct QualityStringVertexData
    string_value::String
    evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}

    function QualityStringVertexData(string_value::String)
        evidence = Dict{String, Dict{String, Set{QualityEvidenceEntry}}}()
        new(string_value, evidence)
    end
end

# ============================================================================
# Lightweight Vertex Data Types
#
# These types store only topology and aggregated counts instead of full evidence.
# They satisfy the same interface as their full counterparts via method overloads
# in evidence-functions.jl, enabling a ~1875x memory reduction for analysis
# workflows that only need count_total_observations / count_evidence_entries.
#
# Use lightweight=true in public graph construction functions to enable.
# ============================================================================

"""
    LightweightKmerVertexData{T}

Lightweight vertex data for k-mer graphs that stores only aggregated counts
and unique observation tracking.

# Type Parameter
- `T`: K-mer type (e.g., DNAKmer{31}, RNAKmer{21}, AAKmer{11})

# Fields
- `Kmer::T`: The k-mer AS OBSERVED (not canonical)
- `total_count::Int`: Total evidence entries across all datasets (incremented
  per `add_evidence!` call — counts every read+position pair)
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `dataset_observations::Dict{String, Set{String}}`: Per-dataset unique
  observation IDs (e.g., read identifiers)

# Semantic Equivalence with Full Mode
Both counting functions return identical results in full and lightweight modes:
- `count_total_observations` → unique observation IDs across all datasets
- `count_evidence_entries` → total entries (read+position pairs)

# Memory Savings
Full KmerVertexData stores ~300 bytes per observation per vertex due to nested
Dict/Set overhead (positions, strands per entry). LightweightKmerVertexData
stores only the Set of observation ID strings, not their positional details —
a ~1875x reduction for datasets with many positions per observation.

# Limitations
- Cannot retrieve individual EvidenceEntry structs (position, strand)
- `get_dataset_evidence` returns nothing (no positional data)
- Algorithms requiring strand info should use full mode
"""
mutable struct LightweightKmerVertexData{T}
    const Kmer::T
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const dataset_observations::Dict{String, Set{String}}

    function LightweightKmerVertexData(kmer::T) where {T}
        new{T}(kmer, 0, Dict{String, Int}(), Dict{String, Set{String}}())
    end
end

"""
    LightweightBioSequenceVertexData{T}

Lightweight vertex data for variable-length BioSequence graphs.

# Type Parameter
- `T`: BioSequence type (LongDNA{4}, LongRNA{4}, LongAA)

# Fields
- `sequence::T`: The sequence AS OBSERVED
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `dataset_observations::Dict{String, Set{String}}`: Per-dataset unique
  observation IDs
"""
mutable struct LightweightBioSequenceVertexData{T}
    const sequence::T
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const dataset_observations::Dict{String, Set{String}}

    function LightweightBioSequenceVertexData(sequence::T) where {T}
        new{T}(sequence, 0, Dict{String, Int}(), Dict{String, Set{String}}())
    end
end

"""
    LightweightStringVertexData

Lightweight vertex data for string/n-gram graphs.

# Fields
- `string_value::String`: The string AS OBSERVED
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `dataset_observations::Dict{String, Set{String}}`: Per-dataset unique
  observation IDs
"""
mutable struct LightweightStringVertexData
    const string_value::String
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const dataset_observations::Dict{String, Set{String}}

    function LightweightStringVertexData(string_value::String)
        new(string_value, 0, Dict{String, Int}(), Dict{String, Set{String}}())
    end
end

# ============================================================================
# Ultralight Vertex Data Types
#
# These types store only topology and aggregated counts — NO observation ID
# tracking. They satisfy the same interface as lightweight types but return
# total_count instead of unique observation count for count_total_observations.
#
# Use memory_profile=:ultralight in public graph construction functions.
# ============================================================================

"""
    UltralightKmerVertexData{T}

Ultralight vertex data for k-mer graphs. Stores only aggregated counts per
dataset — no observation ID tracking and no positional evidence.

# Type Parameter
- `T`: K-mer type (e.g., DNAKmer{31}, RNAKmer{21}, AAKmer{11})

# Fields
- `Kmer::T`: The k-mer AS OBSERVED (not canonical)
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts

# Semantic Difference from Lightweight
- `count_total_observations` → returns `total_count` (NOT unique obs count)
- `get_all_observation_ids` → returns `nothing` (no obs ID tracking)
- All other count/topology functions work identically

# Memory Savings
Drops both positional evidence AND observation ID sets. For 50K vertices with
5K observations: ~40 MB vs ~100 MB (lightweight) vs ~75 GB (full).
"""
mutable struct UltralightKmerVertexData{T}
    const Kmer::T
    total_count::Int
    const dataset_counts::Dict{String, Int}

    function UltralightKmerVertexData(kmer::T) where {T}
        new{T}(kmer, 0, Dict{String, Int}())
    end
end

"""
    UltralightBioSequenceVertexData{T}

Ultralight vertex data for variable-length BioSequence graphs.

# Type Parameter
- `T`: BioSequence type (LongDNA{4}, LongRNA{4}, LongAA)

# Fields
- `sequence::T`: The sequence AS OBSERVED
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
"""
mutable struct UltralightBioSequenceVertexData{T}
    const sequence::T
    total_count::Int
    const dataset_counts::Dict{String, Int}

    function UltralightBioSequenceVertexData(sequence::T) where {T}
        new{T}(sequence, 0, Dict{String, Int}())
    end
end

"""
    UltralightStringVertexData

Ultralight vertex data for string/n-gram graphs.

# Fields
- `string_value::String`: The string AS OBSERVED
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
"""
mutable struct UltralightStringVertexData
    const string_value::String
    total_count::Int
    const dataset_counts::Dict{String, Int}

    function UltralightStringVertexData(string_value::String)
        new(string_value, 0, Dict{String, Int}())
    end
end

# ============================================================================
# Ultralight Quality Vertex Data Types
#
# These types store aggregated counts AND joint quality scores per position,
# but NO observation ID tracking. Quality scores are aggregated incrementally
# during add_evidence! using Phred score addition (equivalent to multiplying
# error probabilities — associative and commutative, so incremental aggregation
# is mathematically exact).
#
# Use memory_profile=:ultralight_quality in public graph construction functions.
# ============================================================================

"""
    UltralightQualityKmerVertexData{T}

Ultralight quality-aware vertex data for k-mer graphs.

# Type Parameter
- `T`: K-mer type (DNA or RNA k-mers — amino acids don't have quality)

# Fields
- `Kmer::T`: The k-mer AS OBSERVED (not canonical)
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `joint_quality::Vector{UInt8}`: Per-position joint Phred scores across ALL observations
- `dataset_joint_quality::Dict{String, Vector{UInt8}}`: Per-dataset joint quality
"""
mutable struct UltralightQualityKmerVertexData{T}
    const Kmer::T
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const joint_quality::Vector{UInt8}
    const dataset_joint_quality::Dict{String, Vector{UInt8}}

    function UltralightQualityKmerVertexData(kmer::T) where {T}
        new{T}(kmer, 0, Dict{String, Int}(), Vector{UInt8}(), Dict{String, Vector{UInt8}}())
    end
end

"""
    UltralightQualityBioSequenceVertexData{T}

Ultralight quality-aware vertex data for variable-length BioSequence graphs.

# Type Parameter
- `T`: BioSequence type (LongDNA{4}, LongRNA{4})

# Fields
- `sequence::T`: The sequence AS OBSERVED
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `joint_quality::Vector{UInt8}`: Per-position joint Phred scores
- `dataset_joint_quality::Dict{String, Vector{UInt8}}`: Per-dataset joint quality
"""
mutable struct UltralightQualityBioSequenceVertexData{T}
    const sequence::T
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const joint_quality::Vector{UInt8}
    const dataset_joint_quality::Dict{String, Vector{UInt8}}

    function UltralightQualityBioSequenceVertexData(sequence::T) where {T}
        new{T}(sequence, 0, Dict{String, Int}(),
            Vector{UInt8}(), Dict{String, Vector{UInt8}}())
    end
end

# ============================================================================
# Lightweight Quality Vertex Data Types
#
# These types store aggregated counts, observation ID tracking, AND joint
# quality scores. This is the most expensive reduced mode but enables both
# multi-sample comparison (via obs IDs) and quality-aware assembly.
#
# Use memory_profile=:lightweight_quality in public graph construction functions.
# ============================================================================

"""
    LightweightQualityKmerVertexData{T}

Lightweight quality-aware vertex data for k-mer graphs. Combines observation
tracking from lightweight mode with quality aggregation.

# Type Parameter
- `T`: K-mer type (DNA or RNA k-mers)

# Fields
- `Kmer::T`: The k-mer AS OBSERVED (not canonical)
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `dataset_observations::Dict{String, Set{String}}`: Per-dataset unique observation IDs
- `joint_quality::Vector{UInt8}`: Per-position joint Phred scores
- `dataset_joint_quality::Dict{String, Vector{UInt8}}`: Per-dataset joint quality
"""
mutable struct LightweightQualityKmerVertexData{T}
    const Kmer::T
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const dataset_observations::Dict{String, Set{String}}
    const joint_quality::Vector{UInt8}
    const dataset_joint_quality::Dict{String, Vector{UInt8}}

    function LightweightQualityKmerVertexData(kmer::T) where {T}
        new{T}(kmer, 0, Dict{String, Int}(), Dict{String, Set{String}}(),
            Vector{UInt8}(), Dict{String, Vector{UInt8}}())
    end
end

"""
    LightweightQualityBioSequenceVertexData{T}

Lightweight quality-aware vertex data for variable-length BioSequence graphs.

# Type Parameter
- `T`: BioSequence type (LongDNA{4}, LongRNA{4})

# Fields
- `sequence::T`: The sequence AS OBSERVED
- `total_count::Int`: Total evidence entries across all datasets
- `dataset_counts::Dict{String, Int}`: Per-dataset evidence entry counts
- `dataset_observations::Dict{String, Set{String}}`: Per-dataset unique observation IDs
- `joint_quality::Vector{UInt8}`: Per-position joint Phred scores
- `dataset_joint_quality::Dict{String, Vector{UInt8}}`: Per-dataset joint quality
"""
mutable struct LightweightQualityBioSequenceVertexData{T}
    const sequence::T
    total_count::Int
    const dataset_counts::Dict{String, Int}
    const dataset_observations::Dict{String, Set{String}}
    const joint_quality::Vector{UInt8}
    const dataset_joint_quality::Dict{String, Vector{UInt8}}

    function LightweightQualityBioSequenceVertexData(sequence::T) where {T}
        new{T}(sequence, 0, Dict{String, Int}(), Dict{String, Set{String}}(),
            Vector{UInt8}(), Dict{String, Vector{UInt8}}())
    end
end

# Convenience constructors for parameterized vertex types to support explicit instantiation.
KmerVertexData{T}(kmer::T) where {T} = KmerVertexData(kmer)
QualmerVertexData{T}(kmer::T) where {T} = QualmerVertexData(kmer)
BioSequenceVertexData{T}(sequence::T) where {T} = BioSequenceVertexData(sequence)
LightweightKmerVertexData{T}(kmer::T) where {T} = LightweightKmerVertexData(kmer)
function LightweightBioSequenceVertexData{T}(sequence::T) where {T}
    LightweightBioSequenceVertexData(sequence)
end
function QualityBioSequenceVertexData{T}(sequence::T) where {T}
    QualityBioSequenceVertexData(sequence)
end
function QualityBioSequenceVertexData{T}(sequence::T, quality_scores::Vector{UInt8}) where {T}
    QualityBioSequenceVertexData(sequence, quality_scores)
end
UltralightKmerVertexData{T}(kmer::T) where {T} = UltralightKmerVertexData(kmer)
function UltralightBioSequenceVertexData{T}(sequence::T) where {T}
    UltralightBioSequenceVertexData(sequence)
end
function UltralightQualityKmerVertexData{T}(kmer::T) where {T}
    UltralightQualityKmerVertexData(kmer)
end
function UltralightQualityBioSequenceVertexData{T}(sequence::T) where {T}
    UltralightQualityBioSequenceVertexData(sequence)
end
function LightweightQualityKmerVertexData{T}(kmer::T) where {T}
    LightweightQualityKmerVertexData(kmer)
end
function LightweightQualityBioSequenceVertexData{T}(sequence::T) where {T}
    LightweightQualityBioSequenceVertexData(sequence)
end
