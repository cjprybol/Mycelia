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

    function QualityBioSequenceVertexData(sequence::T, quality_scores::Vector{UInt8}=Vector{UInt8}()) where {T}
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

# Convenience constructors for parameterized vertex types to support explicit instantiation.
KmerVertexData{T}(kmer::T) where {T} = KmerVertexData(kmer)
QualmerVertexData{T}(kmer::T) where {T} = QualmerVertexData(kmer)
BioSequenceVertexData{T}(sequence::T) where {T} = BioSequenceVertexData(sequence)
QualityBioSequenceVertexData{T}(sequence::T) where {T} = QualityBioSequenceVertexData(sequence)
QualityBioSequenceVertexData{T}(sequence::T, quality_scores::Vector{UInt8}) where {T} =
    QualityBioSequenceVertexData(sequence, quality_scores)
