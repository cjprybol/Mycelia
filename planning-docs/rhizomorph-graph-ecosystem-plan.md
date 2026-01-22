# Comprehensive Rhizomorph Graph Ecosystem Action Plan

> For a quick status overview of supported graph types and modes, see [RHIZOMORPH_SUPPORT_MATRIX.md](RHIZOMORPH_SUPPORT_MATRIX.md).

A detailed plan for implementing a robust "rhizomorph" graph ecosystem that properly handles all graph types as **directed, strand-aware graphs** where vertices and edges are added only when observed. The design distinguishes between two orthogonal concepts:

## ✅ Status Snapshot (current)
- Implemented in `src/rhizomorph`: enums, evidence structs/functions, quality utilities, graph query helpers, strand-specific singlestrand builders for k-mer/qualmer/n-gram and variable-length string/FASTA/FASTQ graphs, DNA/RNA doublestrand + canonical conversions (fixed and variable-length), Eulerian path + path_to_sequence fixes, GFA export/import, and tip removal.
- Tests now on Rhizomorph: `basic_graph_tests`, `sequence_graphs_next`, `graph_algorithms_next` (bubble helpers, DFS fallback, tip thresholds), `end_to_end_graph_tests` (exact k-mer content/evidence positions, qualmer edge/vertex quality evidence, variable-length overlap evidence, doublestrand/canonical conversion checks, GFA vertex-count round-trip), canonicalization/singlestrand orientation suites, `gfa_io_next`, strand-specific/canonical traversal suites, and GFA parsing basics.
- Missing/planned: error-correction; expand simplification edge-removal coverage (bubble resolution/inner-edge cleanup); add traversal coverage for variable-length/n-gram graphs beyond smoke tests.
- Migration work needed: port remaining high-touch legacy suites (`end_to_end_assembly_tests.jl`, `six_graph_hierarchy_tests.jl`, `comprehensive_*`, `probabilistic_algorithms_next.jl`, older string/assembly end-to-end suites) to `Mycelia.Rhizomorph` APIs, then drop legacy includes (`graph-core.jl`, `kmer-graphs.jl`, `sequence-graphs-next.jl`, `string-graphs.jl`, `qualmer-analysis.jl`, `qualmer-graphs.jl`, `fasta-graphs.jl`, `fastq-graphs.jl`).
- Testing gaps: variable-length and n-gram traversal still light; simplification needs broader edge-removal tests beyond bubble detection/tip pruning smoke tests.

1. **Strand Specificity** (methodological): Whether evidence is preserved per strand orientation or merged across reverse-complement pairs
   - **Strand-specific**: Preserve strand orientation as observed (default construction mode)
   - **Non-strand-specific**: Merge evidence across reverse-complement pairs (post-processing transformation)

2. **Strand Representation** (structural): How many strand orientations exist in the graph
   - **Single-stranded**: Only one strand orientation present (biological reality or filtering choice)
   - **Double-stranded**: Both forward and reverse-complement strands present with bidirectional relationships

---

## 🎯 Core Architecture Principles

### 1. **Directed, Strand-Aware Graphs**
- **All graphs are directed graphs** where vertices have explicit strand orientation
- **Record exactly what we observe**: vertices and edges added only when observed (no inference)
- **Evidence tracking**: `(observation_id, position, strand_orientation)` for all observations
- **Edge construction**: Built from observed transitions (n-grams, k-mers) or alignments/overlaps (variable-length)

### 2. **Strand Specificity (Methodological Choice)**
- **Strand-specific (default)**: Preserve strand orientation exactly as observed
  - Evidence for Forward orientation stored separately from Reverse
  - No canonicalization during construction
  - Biological orientation preserved
- **Non-strand-specific (post-processing)**: Merge evidence across reverse-complement pairs
  - Union operation combining Forward and Reverse evidence
  - Applied after initial strand-specific construction
  - Reversible: filter back to strand-specific by selecting one orientation

### 3. **Strand Representation (Structural Choice)**
- **Single-stranded representation**: Graph contains one strand orientation
  - May reflect biological reality (ssRNA viruses, mRNA)
  - May be filtering choice for double-stranded molecules
  - Simpler graph structure
- **Double-stranded representation**: Graph contains both RC orientations
  - Both strands present with bidirectional relationships
  - May reflect biological reality (dsDNA)
  - Created by replicating vertices/edges in RC orientation
  - Reversible: collapse to single-stranded by filtering

### 4. **Orthogonality of Concepts**
These are **independent** choices:
- Strand-specific + Single-stranded: mRNA with strand-specific library prep
- Strand-specific + Double-stranded: dsDNA genome with strand-specific library prep
- Non-strand-specific + Single-stranded: ssRNA virus with non-specific library prep (shows both orientations)
- Non-strand-specific + Double-stranded: dsDNA genome with non-specific library prep

### 5. **Type System Consistency**
- **Never convert to strings** unless input type is String
- Maintain BioSequence types throughout k-mer/qualmer pipelines
- Use proper Phred scores (0-60) not ASCII characters
- Generic parameterized structures: `VertexData{T}`, `EdgeData{T}`

---

## Graph Representation Modes (Single, Double, Canonical) — Explicit Definitions

All graphs are constructed first as **strand-specific single-strand directed graphs** that record only what is observed. DNA/RNA graphs can be transformed into double-strand (directed) or canonical (undirected) representations; amino acid and string graphs remain single-strand.

- **Single-strand (directed; all alphabets)**: Vertices/edges exactly as observed; evidence is strand-aware but only forward observations exist. This is the default construction.
- **Double-strand (directed; DNA/RNA only)**: Forward and reverse-complement vertices/edges both exist. Evidence is the joint set from forward and reverse, but strand flags are preserved so forward and reverse entries are equal-and-opposite in orientation.
- **Canonical (undirected storage; DNA/RNA only)**: Forward/RC vertices are merged into canonical labels and edges are undirected. Evidence is merged but retains strand flags to recover direction during traversal.

### Interconversion Rules
- **Single → Double**: Add reverse-complement vertices and directed edges; merge evidence; keep strand flags.
- **Single → Canonical**: Merge RC pairs into canonical labels; make edges undirected; keep strand flags for traversal validation.
- **Double → Single**: Filter to forward/+ vertices, edges, and evidence (drop RC copies).
- **Canonical → Single**: Expand canonical vertices into forward/RC directed copies; emit only forward-supported vertices/edges/evidence.
- **Double ↔ Canonical**: Merge RC pairs into canonical labels (Double → Canonical) or split canonical labels back into forward/RC directed pairs (Canonical → Double) while preserving strand-aware evidence.

### Flow Diagram (text)
```
Single-strand (directed, all alphabets)
    |-- add RC vertices/edges + merge evidence --> Double-strand (directed, DNA/RNA)
    |-- canonicalize labels + undirected edges --> Canonical (undirected, DNA/RNA)

Double-strand (directed, DNA/RNA)
    |-- keep only forward vertices/edges/evidence --> Single-strand (directed)
    |-- merge RC pairs into canonical labels ------> Canonical (undirected)

Canonical (undirected, DNA/RNA)
    |-- expand to forward/RC directed copies (keep strand flags) --> Double-strand (directed)
    |-- expand then keep only forward copies --------------------> Single-strand (directed)
```


## 📁 Proposed Module Organization

**File Structure and Organization Clarification:**

The module is organized to separate concerns:
- **core/**: Data structures and fundamental operations
- **fixed-length/**: Graph construction for k-mers, qualmers, n-grams
- **variable-length/**: Graph construction for FASTA, FASTQ, strings
- **algorithms/**: Higher-level algorithms that operate on any graph type

**Note on Strand Conversion Placement:**
- Basic strand operations (flipping evidence, finding RC pairs) → `core/evidence-functions.jl` (Section 1.3.1)
- Graph-level conversions (convert between graph modes) → `variable-length/strand-conversions.jl`
- Rationale: Evidence operations are fundamental building blocks; graph conversions are higher-level transformations

```
src/rhizomorph/
├── core/
│   ├── enums.jl                   # StrandOrientation, GraphMode
│   ├── vertex-data.jl             # Generic VertexData{T} structures
│   ├── edge-data.jl               # Generic EdgeData{T} structures
│   ├── evidence-functions.jl      # Evidence manipulation + basic strand operations
│   ├── graph-construction.jl      # Shared graph building logic
│   └── graph-type-conversions.jl  # Fixed↔Variable, Quality↔Non-quality
├── fixed-length/
│   ├── kmer-graphs.jl             # DNAKmer, RNAKmer, AAKmer graphs
│   ├── qualmer-graphs.jl          # Quality-aware k-mer graphs
│   └── ngram-graphs.jl            # String n-gram graphs
├── variable-length/
│   ├── fasta-graphs.jl            # Variable-length BioSequence graphs
│   ├── fastq-graphs.jl            # Variable-length quality BioSequence graphs
│   ├── strand-conversions.jl      # Graph-level strand mode transformations (variable-length)
│   └── string-graphs.jl           # Variable-length string graphs
├── algorithms/
│   ├── path-finding.jl            # Eulerian paths, reconstruction
│   ├── error-correction.jl        # Read-centric probabilistic correction (planned)
│   ├── simplification.jl          # Bubble popping, tip clipping, variant calling (partial)
│   └── io.jl                      # GFA export/import
└── rhizomorph.jl                  # Main module file
```

### Immediate Migration Plan
1. Continue porting legacy graph suites directly to `Mycelia.Rhizomorph` (no shims): gfa IO, end-to-end assembly, six-graph hierarchy, comprehensive correctness/fixes/type-stable, canonicalization consistency, AA/RNA singlestrand, and probabilistic algorithm tests.
2. Implement or replace legacy-only helpers in Rhizomorph where tests still rely on them (repeat detection, contig/coverage helpers, graph-type conversions) with real code + tests; otherwise re-scope tests to existing APIs.
3. Once the remaining suites run on Rhizomorph, remove legacy includes from `src/Mycelia.jl` and delete the deprecated graph files.
4. Expand traversal coverage for variable-length/n-gram modes alongside simplification edge-case tests (JLD2 and GFA round-trips are covered in end-to-end suites).

---

## 🔧 Detailed Migration Plan

### Phase 1: Core Infrastructure (Week 1)

#### 1.1 Create Rhizomorph Module Structure
```bash
mkdir -p src/rhizomorph/{core,fixed-length,variable-length,algorithms}
```

#### 1.1.1 Module Imports and Organization (`src/rhizomorph/rhizomorph.jl`)

**Import Philosophy: Top-Level Packages Only**

Following the Mycelia package standard:
- **ONLY import top-level packages** using `import PackageName`
- **NEVER use `using` statements** or import specific functions
- **Always use full module namespacing**: `BioSequences.LongDNA{4}`, `FASTX.sequence()`, etc.

**Main Module Structure:**

```julia
module Rhizomorph

# ============================================================================
# External Dependencies - Top-Level Packages Only
# ============================================================================

# Core Julia packages
import Base
import Statistics

# Bioinformatics packages
import BioSequences
import Kmers
import FASTX
import BioAlignments

# Graph packages
import MetaGraphsNext
import Graphs

# I/O and serialization
import JLD2      # Lossless graph serialization
import Arrow     # Alternative serialization option

# Optional: Static arrays for fixed-length quality scores
# import StaticArrays

# ============================================================================
# Include All Source Files
# ============================================================================

include("core/enums.jl")
include("core/vertex-data.jl")
include("core/edge-data.jl")
include("core/evidence-functions.jl")
include("core/graph-construction.jl")
include("core/graph-type-conversions.jl")

include("fixed-length/kmer-graphs.jl")
include("fixed-length/qualmer-graphs.jl")
include("fixed-length/ngram-graphs.jl")

include("variable-length/fasta-graphs.jl")
include("variable-length/fastq-graphs.jl")
include("variable-length/string-graphs.jl")

include("algorithms/path-finding.jl")
include("algorithms/error-correction.jl")
include("algorithms/simplification.jl")
include("algorithms/io.jl")

end # module Rhizomorph
```

**Usage Pattern:**

```julia
# User code - import only top-level Mycelia
import Mycelia

# Use full namespacing
graph = Mycelia.build_kmer_graph(records, k=31)
evidence = Mycelia.get_dataset_evidence(vertex_data, "dataset_01")
qual_phred = Mycelia.probability_to_phred(0.01)
```

**No Explicit Exports:**
- Users access all functions via `Mycelia.function_name()`
- No export declarations needed
- Cleaner namespace, no conflicts
- Clear provenance of all functions

#### 1.2 Move Shared Enums and Types (`src/rhizomorph/core/`)

**`enums.jl`** - Move from `graph-core.jl`:
```julia
@enum StrandOrientation Forward=true Reverse=false
@enum GraphMode SingleStrand DoubleStrand
```

**`vertex-data.jl`** - Consolidate all vertex data structures:
```julia
# Evidence entry structures (observation_id is now a dict key, not stored in entry)
struct EvidenceEntry
    position::Int
    strand::StrandOrientation
end

struct QualityEvidenceEntry
    position::Int
    strand::StrandOrientation
    quality_scores::Vector{UInt8}  # Phred scores 0-60
end

# Fixed-length vertex data
struct KmerVertexData{T}
    Kmer::T
    # dataset_id -> observation_id -> Set{(position, strand)}
    evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}
end

struct QualmerVertexData{T}
    Kmer::T
    # dataset_id -> observation_id -> Set{(position, strand, quality_scores)}
    evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}
end

# Variable-length vertex data
struct BioSequenceVertexData{T}
    sequence::T
    evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}
end

struct QualityBioSequenceVertexData{T}
    sequence::T
    evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}
end

struct StringVertexData
    string_value::String
    evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}
end

struct QualityStringVertexData
    string_value::String
    evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}
end
```

---
**Design Note on Evidence Storage:**

The double-nested dictionary structure `Dict{String, Dict{String, Set{EvidenceEntry}}}` provides:

1. **Efficient dataset-level queries**: O(1) access to all evidence from a specific dataset
2. **Efficient observation-level queries**: O(1) access to all evidence from a specific observation within a dataset
3. **Memory efficiency**: Each observation_id stored once as a dictionary key rather than repeated in every evidence entry
4. **Natural organization**: Matches the hierarchical nature of sequencing data (datasets contain observations)

This design is optimal for:
- Read-centric correction algorithms that need fast per-observation queries
- Comparative analysis across multiple datasets
- High-coverage scenarios where observations contribute multiple positions to vertices
- Long-read sequencing where single reads may revisit vertices multiple times

The structure trades minor iteration complexity for substantial memory savings and query performance, especially important for iterative assembly with varying k-mer sizes and increasingly common long-read technologies.

---

**`edge-data.jl`** - Consolidate edge structures:
```julia
# Edge evidence: observation_id is key in nested structure
struct EdgeEvidenceEntry
    from_position::Int
    to_position::Int
    strand::StrandOrientation
end

struct EdgeQualityEvidenceEntry
    from_position::Int
    to_position::Int
    strand::StrandOrientation
    from_quality::Vector{UInt8}  # Quality scores at source vertex
    to_quality::Vector{UInt8}    # Quality scores at target vertex
end

struct KmerEdgeData
    # dataset_id -> observation_id -> Set{(from_pos, to_pos, strand)}
    evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}
end

struct QualmerEdgeData
    # dataset_id -> observation_id -> Set{(from_pos, to_pos, strand, qualities)}
    evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}
end

struct BioSequenceEdgeData
    overlap_length::Int
    evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}
end

struct StringEdgeData
    overlap_length::Int
    evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}
end

struct QualityBioSequenceEdgeData
    overlap_length::Int
    evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}
end

struct QualityStringEdgeData
    overlap_length::Int
    evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}
end
```

---

#### 1.2.1 Quality Score Specifications and Type Considerations

**Quality Score Vector Lengths and Type Choices:**

The length and type of `quality_scores` depends on context:

1. **Fixed-length k-mers/qualmers** (`QualmerVertexData`):
   - Length = k (the k-mer size)
   - **Type consideration**: Use `NTuple{k, UInt8}` or `SVector{k, UInt8}` for type stability and performance
   - Enforces length at compile time, prevents length mismatches
   - Example: `quality_scores::NTuple{31, UInt8}` for k=31

2. **Fixed-length n-grams** (`QualityStringVertexData` with n-grams):
   - Length = n (the n-gram size)
   - **Type consideration**: Use `NTuple{n, UInt8}` or `SVector{n, UInt8}`
   - Same benefits as k-mers

3. **Variable-length sequences** (`QualityBioSequenceVertexData`):
   - Length = sequence length (varies per vertex)
   - **Type**: Keep as `Vector{UInt8}` (variable length appropriate here)
   - Sequence type (LongDNA, LongRNA, LongAA) is also variable length

4. **Variable-length strings** (`QualityStringVertexData` without fixed n):
   - Length = string length (varies)
   - **Type**: Keep as `Vector{UInt8}`

5. **Edge quality scores** (`EdgeQualityEvidenceEntry`):
   - `from_quality`: Length matches source vertex sequence length
   - `to_quality`: Length matches destination vertex sequence length
   - Use same type strategy as vertices (NTuple for fixed, Vector for variable)

**Recommendation**: Use `NTuple` for fixed-length (simpler, no extra dependency) or `SVector` from StaticArrays.jl (more functionality). This provides type stability and compile-time guarantees.

**Updated Structure Definitions (if using NTuple):**

```julia
# For k-mer graphs with k known at runtime, we still need Vector
# But for specific k values, we could parameterize:

struct QualmerVertexData{T, K}
    Kmer::T
    # K is the k-mer size for type stability
    evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry{K}}}}
end

struct QualityEvidenceEntry{K}
    position::Int
    strand::StrandOrientation
    quality_scores::NTuple{K, UInt8}  # Fixed length K, type-stable
end

# For variable-length, keep Vector:
struct QualityBioSequenceVertexData{T}
    sequence::T
    evidence::Dict{String, Dict{String, Set{VariableLengthQualityEvidenceEntry}}}
end

struct VariableLengthQualityEvidenceEntry
    position::Int
    strand::StrandOrientation
    quality_scores::Vector{UInt8}  # Variable length
end
```

**Phred Score Range:**

- **FASTQ format standard**: Phred scores 0-60 (individual base calls from sequencer)
- **Joint observation scores**: Can exceed 60 up to **typemax(UInt8) = 255**
  - As we observe the same k-mer multiple times, confidence increases
  - Phred 60 ≈ 99.9999% confidence (1 in 1,000,000 error)
  - Phred 100 ≈ 1 in 10 billion error
  - Phred 255 ≈ essentially certainty
  - **This allows differentiation between**: low confidence (10-30), moderate (30-60), high (60-100), extreme (100+)

**Phred Score Mathematics:**

```julia
# ============================================================================
# Phred Score Conversion Functions
# ============================================================================

"""
Convert error probability to Phred score.
Formula: Q = -10 * log10(P_error)

Examples:
- P_error = 0.1 (90% correct)  → Q = 10
- P_error = 0.01 (99% correct) → Q = 20
- P_error = 0.001 (99.9% correct) → Q = 30
- P_error = 1e-10 (99.99999999% correct) → Q = 100
"""
function probability_to_phred(p_error::Float64)
    if p_error <= 0.0
        return 255.0  # Maximum representable quality
    elseif p_error >= 1.0
        return 0.0    # Minimum quality
    else
        q = -10.0 * log10(p_error)
        return min(q, 255.0)  # Clamp to UInt8 range
    end
end

"""
Convert Phred score to error probability.
Formula: P_error = 10^(-Q/10)

Examples:
- Q = 10 → P_error = 0.1 (90% correct)
- Q = 20 → P_error = 0.01 (99% correct)
- Q = 30 → P_error = 0.001 (99.9% correct)
- Q = 60 → P_error = 1e-6 (99.9999% correct)
- Q = 100 → P_error = 1e-10 (99.99999999% correct)
"""
function phred_to_probability(phred::Float64)
    return 10.0^(-phred/10.0)
end

"""
Convert Phred score to correctness probability.
Formula: P_correct = 1 - P_error = 1 - 10^(-Q/10)
"""
function phred_to_correctness(phred::Float64)
    return 1.0 - phred_to_probability(phred)
end

# ============================================================================
# Quality Score Aggregation (Independence Assumption)
# ============================================================================

"""
Aggregate quality scores from multiple observations using independence assumption.

IMPORTANT: This is the primary/default aggregation strategy.

When the same k-mer/sequence is observed multiple times with different quality scores,
combine them assuming independence:

Mathematical derivation:
- Single observation: Q₁ → P₁(error) = 10^(-Q₁/10)
- Second observation: Q₂ → P₂(error) = 10^(-Q₂/10)
- **Independence assumption**: P(both errors) = P₁(error) × P₂(error)
- **In Phred (log) space, this is ADDITIVE**:
  - Log₁₀(P_both_errors) = Log₁₀(P₁ × P₂) = Log₁₀(P₁) + Log₁₀(P₂)
  - Since Q = -10 × Log₁₀(P), this means: **Q_combined = Q₁ + Q₂**

Example:
- Two observations: Q₁=10 (P_error=0.1), Q₂=10 (P_error=0.1)
- P(both errors) = 0.1 × 0.1 = 0.01
- Q_combined = -10 × log₁₀(0.01) = 20 ✓
- **Simple formula**: Q_combined = Q₁ + Q₂ = 10 + 10 = 20 ✓

Multiple observations:
- N observations with qualities [Q₁, Q₂, ..., Qₙ]
- **Q_combined = Q₁ + Q₂ + ... + Qₙ** (can exceed 60, up to 255)

This allows us to differentiate:
- Low confidence: 10-30 (few low-quality observations)
- Moderate: 30-60 (typical single high-quality observation)
- High: 60-100 (multiple high-quality observations)
- Extreme: 100-255 (many high-quality observations = essentially certain)
"""
function aggregate_quality_scores_independence(quality_scores::Vector{Vector{UInt8}})
    if isempty(quality_scores)
        return UInt8[]
    end

    # All quality vectors must have same length
    k = length(quality_scores[1])
    @assert all(length(q) == k for q in quality_scores) "All quality score vectors must have same length"

    aggregated = Vector{UInt8}(undef, k)

    # Aggregate position by position
    for pos in 1:k
        # Sum Phred scores (equivalent to multiplying error probabilities)
        summed_phred = sum(Float64(qs[pos]) for qs in quality_scores)

        # Clamp to UInt8 range [0, 255]
        aggregated[pos] = UInt8(clamp(round(summed_phred), 0, 255))
    end

    return aggregated
end

"""
Conservative quality aggregation (alternative method).

Instead of independence assumption, use correctness probabilities:
- P(all correct) = P₁(correct) × P₂(correct) × ... × Pₙ(correct)
- P(error) = 1 - P(all correct)
- Q_combined = -10 × log₁₀(P_error)

This is MORE CONSERVATIVE than independence assumption for error probabilities.

Example:
- Two observations: Q₁=10, Q₂=10
- P₁(correct) = 0.9, P₂(correct) = 0.9
- P(both correct) = 0.9 × 0.9 = 0.81
- P(error) = 1 - 0.81 = 0.19
- Q_combined = -10 × log₁₀(0.19) ≈ 7.2

Compare to independence: Q_combined = 20

Use this when you want to be conservative and account for possibility
that multiple observations might share systematic errors.
"""
function aggregate_quality_scores_conservative(quality_scores::Vector{Vector{UInt8}})
    if isempty(quality_scores)
        return UInt8[]
    end

    k = length(quality_scores[1])
    @assert all(length(q) == k for q in quality_scores) "All quality score vectors must have same length"

    aggregated = Vector{UInt8}(undef, k)

    # Aggregate position by position
    for pos in 1:k
        # Collect all quality scores at this position
        position_quals = [Float64(qs[pos]) for qs in quality_scores]

        # Convert to correctness probabilities
        p_correct_all = [phred_to_correctness(q) for q in position_quals]

        # Multiply probabilities (conservative assumption)
        p_all_correct = prod(p_correct_all)

        # Convert back to Phred
        p_error = 1.0 - p_all_correct
        q_combined = probability_to_phred(p_error)

        # Clamp to UInt8 range [0, 255]
        aggregated[pos] = UInt8(clamp(round(q_combined), 0, 255))
    end

    return aggregated
end

"""
Simple heuristic aggregations (for comparison/debugging):
"""
function aggregate_quality_scores_max(quality_scores::Vector{Vector{UInt8}})
    k = length(quality_scores[1])
    return [maximum(qs[pos] for qs in quality_scores) for pos in 1:k]
end

function aggregate_quality_scores_mean(quality_scores::Vector{Vector{UInt8}})
    k = length(quality_scores[1])
    return [UInt8(round(sum(qs[pos] for qs in quality_scores) / length(quality_scores))) for pos in 1:k]
end
```

**Default Strategy**: Use `aggregate_quality_scores_independence` (additive in Phred space)
**Alternative**: Use `aggregate_quality_scores_conservative` when concerned about systematic errors

---

#### 1.2.2 Error Rate Model and Dataset/Observation ID Assignment

**Error Rate Model Definition:**

```julia
"""
Inferred or pre-trained error rate model.
Captures insertion, deletion, and substitution rates for a sequencing workflow.
"""
struct ErrorRateModel
    insertion_rate::Float64      # Probability of insertion error
    deletion_rate::Float64        # Probability of deletion error
    substitution_rate::Float64    # Probability of substitution error
    confidence::Float64           # Confidence in these rates [0.0, 1.0]
    source::Symbol                # :inferred_from_graph, :pretrained, :default
end

# Default error models for common platforms
function illumina_error_model()
    ErrorRateModel(0.0003, 0.0003, 0.0004, 0.5, :default)  # ~0.1% total
end

function nanopore_error_model()
    ErrorRateModel(0.03, 0.04, 0.01, 0.5, :default)  # ~8% total
end

function pacbio_hifi_error_model()
    ErrorRateModel(0.0003, 0.0003, 0.0004, 0.5, :default)  # ~0.1% total
end
```

**Dataset and Observation ID Assignment:**

```julia
"""
Extract observation ID from FASTX record.
Uses read identifier (sequence name) up to first whitespace.

Example: @SRR123456.1 1 length=150  →  "SRR123456.1"
"""
function get_observation_id(record::Union{FASTX.FASTA.Record, FASTX.FASTQ.Record})
    identifier = FASTX.identifier(record)
    return String(split(identifier, ' ')[1])
end

"""
Generate observation ID from index (for synthetic data or non-unique read names).
"""
function get_observation_id_from_index(index::Int, prefix::String="obs")
    return "$(prefix)_$(lpad(index, 10, '0'))"
end

"""
Dataset ID assignment strategy:

DEFAULT: Always use input filename (unless user specifies alternative)
- Automatically handles multi-file workflows
- Ensures traceability back to source data
- Only use generic default ("dataset_01") when input is in-memory records
  with no associated file

Examples:
1. Single file: dataset_id = "sample_A"  (from sample_A.fastq)
2. Multiple files: auto-extracted from each filename
3. User override: dataset_id = "patient_001_tumor"  (explicit metadata)
4. In-memory records: dataset_id = "dataset_01"  (no file available)
"""
function get_dataset_id_from_file(filepath::String)
    return splitext(basename(filepath))[1]
end

# Example graph construction with automatic dataset ID:
function build_graph_from_fastq_file(
    filepath::String;
    dataset_id::Union{String,Nothing}=nothing,  # User can override
    k::Int=31
)
    # DEFAULT: Use filename as dataset ID
    if isnothing(dataset_id)
        dataset_id = get_dataset_id_from_file(filepath)
    end

    records = collect(FASTX.FASTQ.Reader(open(filepath)))

    graph = build_graph_from_records(records, dataset_id, k)

    return graph
end

# For in-memory records (no file):
function build_graph_from_records(
    records::Vector{<:FASTX.Record},
    dataset_id::String="dataset_01",  # Generic default only when no file
    k::Int=31
)
    # ... build graph with dataset_id
end
```

---

#### 1.3 Evidence Manipulation Functions (`src/rhizomorph/core/evidence-functions.jl`)

**Helper Functions for Evidence Management:**

These functions abstract the complexity of the double-nested dictionary structure and provide a clean API for working with evidence.

```julia
# ============================================================================
# Adding Evidence
# ============================================================================

"""
Add evidence entry to vertex data structure.
"""
function add_evidence!(
    vertex_data::Union{KmerVertexData, BioSequenceVertexData, StringVertexData},
    dataset_id::String,
    observation_id::String,
    position::Int,
    strand::StrandOrientation
)
    # Ensure dataset exists
    if !haskey(vertex_data.evidence, dataset_id)
        vertex_data.evidence[dataset_id] = Dict{String, Set{EvidenceEntry}}()
    end

    # Ensure observation exists
    if !haskey(vertex_data.evidence[dataset_id], observation_id)
        vertex_data.evidence[dataset_id][observation_id] = Set{EvidenceEntry}()
    end

    # Add evidence entry
    push!(vertex_data.evidence[dataset_id][observation_id],
          EvidenceEntry(position, strand))
end

"""
Add quality-aware evidence entry to vertex data structure.
"""
function add_evidence!(
    vertex_data::Union{QualmerVertexData, QualityBioSequenceVertexData, QualityStringVertexData},
    dataset_id::String,
    observation_id::String,
    position::Int,
    strand::StrandOrientation,
    quality_scores::Vector{UInt8}
)
    if !haskey(vertex_data.evidence, dataset_id)
        vertex_data.evidence[dataset_id] = Dict{String, Set{QualityEvidenceEntry}}()
    end

    if !haskey(vertex_data.evidence[dataset_id], observation_id)
        vertex_data.evidence[dataset_id][observation_id] = Set{QualityEvidenceEntry}()
    end

    push!(vertex_data.evidence[dataset_id][observation_id],
          QualityEvidenceEntry(position, strand, quality_scores))
end

"""
Add edge evidence entry.
"""
function add_edge_evidence!(
    edge_data::Union{KmerEdgeData, BioSequenceEdgeData, StringEdgeData},
    dataset_id::String,
    observation_id::String,
    from_position::Int,
    to_position::Int,
    strand::StrandOrientation
)
    if !haskey(edge_data.evidence, dataset_id)
        edge_data.evidence[dataset_id] = Dict{String, Set{EdgeEvidenceEntry}}()
    end

    if !haskey(edge_data.evidence[dataset_id], observation_id)
        edge_data.evidence[dataset_id][observation_id] = Set{EdgeEvidenceEntry}()
    end

    push!(edge_data.evidence[dataset_id][observation_id],
          EdgeEvidenceEntry(from_position, to_position, strand))
end

# ============================================================================
# Querying Evidence
# ============================================================================

"""
Get all evidence for a specific dataset (coverage-only version).
Returns iterator of (observation_id, position, strand) tuples.
"""
function get_dataset_evidence(
    vertex_data::Union{KmerVertexData, BioSequenceVertexData, StringVertexData},
    dataset_id::String
)
    if !haskey(vertex_data.evidence, dataset_id)
        return Iterators.empty()
    end

    return Iterators.flatten(
        ((obs_id, entry.position, entry.strand)
         for entry in evidence_set)
        for (obs_id, evidence_set) in vertex_data.evidence[dataset_id]
    )
end

"""
Get all evidence for a specific dataset (quality-aware version).
Returns iterator of (observation_id, position, strand, quality_scores) tuples.
"""
function get_dataset_evidence(
    vertex_data::Union{QualmerVertexData, QualityBioSequenceVertexData, QualityStringVertexData},
    dataset_id::String
)
    if !haskey(vertex_data.evidence, dataset_id)
        return Iterators.empty()
    end

    return Iterators.flatten(
        ((obs_id, entry.position, entry.strand, entry.quality_scores)
         for entry in evidence_set)
        for (obs_id, evidence_set) in vertex_data.evidence[dataset_id]
    )
end

"""
Get all evidence for a specific observation (coverage-only).
Returns Set{EvidenceEntry} or empty set if not found.
"""
function get_observation_evidence(
    vertex_data::Union{KmerVertexData, BioSequenceVertexData, StringVertexData},
    dataset_id::String,
    observation_id::String
)
    if !haskey(vertex_data.evidence, dataset_id)
        return Set{EvidenceEntry}()
    end

    if !haskey(vertex_data.evidence[dataset_id], observation_id)
        return Set{EvidenceEntry}()
    end

    return vertex_data.evidence[dataset_id][observation_id]
end

"""
Get all evidence for a specific observation (quality-aware).
Returns Set{QualityEvidenceEntry} or empty set if not found.
"""
function get_observation_evidence(
    vertex_data::Union{QualmerVertexData, QualityBioSequenceVertexData, QualityStringVertexData},
    dataset_id::String,
    observation_id::String
)
    if !haskey(vertex_data.evidence, dataset_id)
        return Set{QualityEvidenceEntry}()
    end

    if !haskey(vertex_data.evidence[dataset_id], observation_id)
        return Set{QualityEvidenceEntry}()
    end

    return vertex_data.evidence[dataset_id][observation_id]
end

"""
Get all evidence across all datasets.
Returns iterator of (dataset_id, observation_id, position, strand) tuples.
"""
function get_all_evidence(vertex_data)
    return Iterators.flatten(
        ((dataset_id, obs_id, entry.position, entry.strand)
         for (obs_id, evidence_set) in dataset_dict
         for entry in evidence_set)
        for (dataset_id, dataset_dict) in vertex_data.evidence
    )
end

"""
Count total number of evidence entries.
"""
function count_evidence(vertex_data)
    return sum(
        length(evidence_set)
        for dataset_dict in values(vertex_data.evidence)
        for evidence_set in values(dataset_dict)
    )
end

"""
Count number of unique observations supporting this vertex.
"""
function count_observations(vertex_data)
    return sum(
        length(dataset_dict)
        for dataset_dict in values(vertex_data.evidence)
    )
end

"""
Count number of datasets contributing evidence.
"""
function count_datasets(vertex_data)
    return length(vertex_data.evidence)
end

"""
Get all observations contributing to this vertex across all datasets.
Returns Set of (dataset_id, observation_id) pairs.
"""
function get_all_observations(vertex_data)
    return Set(
        (dataset_id, obs_id)
        for (dataset_id, dataset_dict) in vertex_data.evidence
        for obs_id in keys(dataset_dict)
    )
end

# ============================================================================
# Filtering Evidence
# ============================================================================

"""
Filter evidence by strand orientation.
Returns new vertex data with only specified strand.
"""
function filter_by_strand(vertex_data::KmerVertexData{T},
                          target_strand::StrandOrientation) where T
    new_evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()

    for (dataset_id, dataset_dict) in vertex_data.evidence
        new_evidence[dataset_id] = Dict{String, Set{EvidenceEntry}}()

        for (obs_id, evidence_set) in dataset_dict
            filtered = filter(e -> e.strand == target_strand, evidence_set)
            if !isempty(filtered)
                new_evidence[dataset_id][obs_id] = filtered
            end
        end

        # Remove empty datasets
        if isempty(new_evidence[dataset_id])
            delete!(new_evidence, dataset_id)
        end
    end

    return KmerVertexData{T}(vertex_data.Kmer, new_evidence)
end

"""
Filter evidence to only include specified datasets.
"""
function filter_by_datasets(vertex_data::KmerVertexData{T},
                            dataset_ids::Set{String}) where T
    new_evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()

    for dataset_id in dataset_ids
        if haskey(vertex_data.evidence, dataset_id)
            new_evidence[dataset_id] = vertex_data.evidence[dataset_id]
        end
    end

    return KmerVertexData{T}(vertex_data.Kmer, new_evidence)
end

"""
Filter evidence to only include specified observations.
"""
function filter_by_observations(vertex_data::KmerVertexData{T},
                                observations::Dict{String, Set{String}}) where T
    # observations format: dataset_id -> Set{observation_id}
    new_evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()

    for (dataset_id, obs_set) in observations
        if haskey(vertex_data.evidence, dataset_id)
            new_evidence[dataset_id] = Dict{String, Set{EvidenceEntry}}()

            for obs_id in obs_set
                if haskey(vertex_data.evidence[dataset_id], obs_id)
                    new_evidence[dataset_id][obs_id] = vertex_data.evidence[dataset_id][obs_id]
                end
            end

            if isempty(new_evidence[dataset_id])
                delete!(new_evidence, dataset_id)
            end
        end
    end

    return KmerVertexData{T}(vertex_data.Kmer, new_evidence)
end

# ============================================================================
# Merging Evidence
# ============================================================================

"""
Merge evidence from two vertices (for non-strand-specific graphs or vertex merging).
Performs union operation on evidence sets.
"""
function merge_evidence!(
    vertex_data1::KmerVertexData{T},
    vertex_data2::KmerVertexData{T}
) where T
    for (dataset_id, dataset_dict) in vertex_data2.evidence
        if !haskey(vertex_data1.evidence, dataset_id)
            vertex_data1.evidence[dataset_id] = Dict{String, Set{EvidenceEntry}}()
        end

        for (obs_id, evidence_set) in dataset_dict
            if !haskey(vertex_data1.evidence[dataset_id], obs_id)
                vertex_data1.evidence[dataset_id][obs_id] = Set{EvidenceEntry}()
            end

            union!(vertex_data1.evidence[dataset_id][obs_id], evidence_set)
        end
    end

    return vertex_data1
end

"""
Merge evidence across reverse-complement strands (for non-strand-specific conversion).
Creates new vertex data with combined evidence from both strand orientations.
"""
function merge_reverse_complement_evidence(
    forward_vertex::KmerVertexData{T},
    reverse_vertex::KmerVertexData{T}
) where T
    # Create deep copy of forward vertex evidence
    merged_evidence = deepcopy(forward_vertex.evidence)

    # Merge in reverse vertex evidence
    temp_vertex = KmerVertexData{T}(forward_vertex.Kmer, merged_evidence)
    merge_evidence!(temp_vertex, reverse_vertex)

    return temp_vertex
end

# ============================================================================
# Comparison and Statistics
# ============================================================================

"""
Calculate evidence depth (total number of observations) at this vertex.
"""
function evidence_depth(vertex_data)
    return count_evidence(vertex_data)
end

"""
Calculate coverage breadth (number of unique observations).
"""
function evidence_breadth(vertex_data)
    return count_observations(vertex_data)
end

"""
Get evidence statistics for a vertex.
"""
function evidence_statistics(vertex_data)
    total_evidence = count_evidence(vertex_data)
    total_observations = count_observations(vertex_data)
    total_datasets = count_datasets(vertex_data)

    # Calculate per-dataset statistics
    dataset_stats = Dict{String, NamedTuple}()
    for (dataset_id, dataset_dict) in vertex_data.evidence
        n_obs = length(dataset_dict)
        n_evidence = sum(length(evidence_set) for evidence_set in values(dataset_dict))
        avg_evidence_per_obs = n_evidence / n_obs

        dataset_stats[dataset_id] = (
            observations = n_obs,
            evidence_entries = n_evidence,
            avg_evidence_per_observation = avg_evidence_per_obs
        )
    end

    return (
        total_evidence = total_evidence,
        total_observations = total_observations,
        total_datasets = total_datasets,
        dataset_statistics = dataset_stats
    )
end

"""
Check if two vertices share any observations (useful for bubble detection).
"""
function share_observations(vertex_data1, vertex_data2)
    obs1 = get_all_observations(vertex_data1)
    obs2 = get_all_observations(vertex_data2)
    return !isempty(intersect(obs1, obs2))
end

# ============================================================================
# Quality-Aware Evidence Functions
# ============================================================================

"""
Filter quality-aware evidence by strand orientation.
Returns new vertex data with only specified strand.
"""
function filter_by_strand(vertex_data::QualmerVertexData{T},
                          target_strand::StrandOrientation) where T
    new_evidence = Dict{String, Dict{String, Set{QualityEvidenceEntry}}}()

    for (dataset_id, dataset_dict) in vertex_data.evidence
        new_evidence[dataset_id] = Dict{String, Set{QualityEvidenceEntry}}()

        for (obs_id, evidence_set) in dataset_dict
            filtered = filter(e -> e.strand == target_strand, evidence_set)
            if !isempty(filtered)
                new_evidence[dataset_id][obs_id] = filtered
            end
        end

        if isempty(new_evidence[dataset_id])
            delete!(new_evidence, dataset_id)
        end
    end

    return QualmerVertexData{T}(vertex_data.Kmer, new_evidence)
end

"""
Filter quality-aware evidence to only include specified datasets.
"""
function filter_by_datasets(vertex_data::QualmerVertexData{T},
                            dataset_ids::Set{String}) where T
    new_evidence = Dict{String, Dict{String, Set{QualityEvidenceEntry}}}()

    for dataset_id in dataset_ids
        if haskey(vertex_data.evidence, dataset_id)
            new_evidence[dataset_id] = vertex_data.evidence[dataset_id]
        end
    end

    return QualmerVertexData{T}(vertex_data.Kmer, new_evidence)
end

"""
Merge quality-aware evidence from two vertices.
Performs union operation on evidence sets.
"""
function merge_evidence!(
    vertex_data1::QualmerVertexData{T},
    vertex_data2::QualmerVertexData{T}
) where T
    for (dataset_id, dataset_dict) in vertex_data2.evidence
        if !haskey(vertex_data1.evidence, dataset_id)
            vertex_data1.evidence[dataset_id] = Dict{String, Set{QualityEvidenceEntry}}()
        end

        for (obs_id, evidence_set) in dataset_dict
            if !haskey(vertex_data1.evidence[dataset_id], obs_id)
                vertex_data1.evidence[dataset_id][obs_id] = Set{QualityEvidenceEntry}()
            end

            union!(vertex_data1.evidence[dataset_id][obs_id], evidence_set)
        end
    end

    return vertex_data1
end

"""
Merge quality-aware evidence across reverse-complement strands.
"""
function merge_reverse_complement_evidence(
    forward_vertex::QualmerVertexData{T},
    reverse_vertex::QualmerVertexData{T}
) where T
    merged_evidence = deepcopy(forward_vertex.evidence)
    temp_vertex = QualmerVertexData{T}(forward_vertex.Kmer, merged_evidence)
    merge_evidence!(temp_vertex, reverse_vertex)
    return temp_vertex
end

# ============================================================================
# Edge Evidence Helper Functions
# ============================================================================

"""
Add quality-aware edge evidence entry.
"""
function add_edge_evidence!(
    edge_data::Union{QualmerEdgeData, QualityBioSequenceEdgeData, QualityStringEdgeData},
    dataset_id::String,
    observation_id::String,
    from_position::Int,
    to_position::Int,
    strand::StrandOrientation,
    from_quality::Vector{UInt8},
    to_quality::Vector{UInt8}
)
    if !haskey(edge_data.evidence, dataset_id)
        edge_data.evidence[dataset_id] = Dict{String, Set{EdgeQualityEvidenceEntry}}()
    end

    if !haskey(edge_data.evidence[dataset_id], observation_id)
        edge_data.evidence[dataset_id][observation_id] = Set{EdgeQualityEvidenceEntry}()
    end

    push!(edge_data.evidence[dataset_id][observation_id],
          EdgeQualityEvidenceEntry(from_position, to_position, strand, from_quality, to_quality))
end

"""
Get all edge evidence for a specific dataset.
"""
function get_dataset_edge_evidence(edge_data, dataset_id::String)
    if !haskey(edge_data.evidence, dataset_id)
        return Iterators.empty()
    end

    return Iterators.flatten(
        ((obs_id, entry.from_position, entry.to_position, entry.strand)
         for entry in evidence_set)
        for (obs_id, evidence_set) in edge_data.evidence[dataset_id]
    )
end

"""
Get all edge evidence for a specific observation.
"""
function get_observation_edge_evidence(edge_data, dataset_id::String, observation_id::String)
    if !haskey(edge_data.evidence, dataset_id)
        return Set{EdgeEvidenceEntry}()
    end

    if !haskey(edge_data.evidence[dataset_id], observation_id)
        return Set{EdgeEvidenceEntry}()
    end

    return edge_data.evidence[dataset_id][observation_id]
end

"""
Count total number of edge evidence entries.
"""
function count_edge_evidence(edge_data)
    return sum(
        length(evidence_set)
        for dataset_dict in values(edge_data.evidence)
        for evidence_set in values(dataset_dict)
    )
end

"""
Filter edge evidence by strand orientation.
"""
function filter_edge_by_strand(edge_data::KmerEdgeData,
                               target_strand::StrandOrientation)
    new_evidence = Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}()

    for (dataset_id, dataset_dict) in edge_data.evidence
        new_evidence[dataset_id] = Dict{String, Set{EdgeEvidenceEntry}}()

        for (obs_id, evidence_set) in dataset_dict
            filtered = filter(e -> e.strand == target_strand, evidence_set)
            if !isempty(filtered)
                new_evidence[dataset_id][obs_id] = filtered
            end
        end

        if isempty(new_evidence[dataset_id])
            delete!(new_evidence, dataset_id)
        end
    end

    return KmerEdgeData(new_evidence)
end

"""
Filter quality-aware edge evidence by strand.
"""
function filter_edge_by_strand(edge_data::QualmerEdgeData,
                               target_strand::StrandOrientation)
    new_evidence = Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}()

    for (dataset_id, dataset_dict) in edge_data.evidence
        new_evidence[dataset_id] = Dict{String, Set{EdgeQualityEvidenceEntry}}()

        for (obs_id, evidence_set) in dataset_dict
            filtered = filter(e -> e.strand == target_strand, evidence_set)
            if !isempty(filtered)
                new_evidence[dataset_id][obs_id] = filtered
            end
        end

        if isempty(new_evidence[dataset_id])
            delete!(new_evidence, dataset_id)
        end
    end

    return QualmerEdgeData(new_evidence)
end

"""
Merge edge evidence (for non-strand-specific graphs).
"""
function merge_edge_evidence!(
    edge_data1::KmerEdgeData,
    edge_data2::KmerEdgeData
)
    for (dataset_id, dataset_dict) in edge_data2.evidence
        if !haskey(edge_data1.evidence, dataset_id)
            edge_data1.evidence[dataset_id] = Dict{String, Set{EdgeEvidenceEntry}}()
        end

        for (obs_id, evidence_set) in dataset_dict
            if !haskey(edge_data1.evidence[dataset_id], obs_id)
                edge_data1.evidence[dataset_id][obs_id] = Set{EdgeEvidenceEntry}()
            end

            union!(edge_data1.evidence[dataset_id][obs_id], evidence_set)
        end
    end

    return edge_data1
end

"""
Merge reverse-complement edge evidence (for non-strand-specific conversion).
"""
function merge_reverse_complement_edge_evidence!(graph::MetaGraphsNext.MetaGraph)
    # Find all reverse-complement edge pairs
    rc_edge_pairs = find_reverse_complement_edge_pairs(graph)

    for (forward_edge, reverse_edge) in rc_edge_pairs
        forward_edge_data = graph[forward_edge...]
        reverse_edge_data = graph[reverse_edge...]

        merge_edge_evidence!(forward_edge_data, reverse_edge_data)
        merge_edge_evidence!(reverse_edge_data, forward_edge_data)
    end

    return graph
end
```

**Design Benefits:**
- **Abstraction**: Complex nested structure hidden behind clean API
- **Type safety**: Multiple dispatch handles different vertex data types
- **Performance**: Functions optimized for common operations
- **Flexibility**: Easy to add new query patterns as needed
- **Read-centric support**: Fast observation-level queries for correction algorithms
- **Complete coverage**: Both vertex and edge evidence, both quality-aware and quality-unaware

---

#### 1.3.0.1 Type Safety and Parameterization Strategy

**Critical for Performance and Correctness:**

Julia's type system and type inference are critical for performance. To maintain type stability:

**1. Return Type Specification:**

```julia
# GOOD: Type-stable with clear return type
function get_observation_evidence(
    vertex_data::KmerVertexData{T},
    dataset_id::String,
    observation_id::String
)::Set{EvidenceEntry} where T
    if !haskey(vertex_data.evidence, dataset_id)
        return Set{EvidenceEntry}()
    end

    if !haskey(vertex_data.evidence[dataset_id], observation_id)
        return Set{EvidenceEntry}()
    end

    return vertex_data.evidence[dataset_id][observation_id]
end

# Quality-aware version with explicit return type
function get_observation_evidence(
    vertex_data::QualmerVertexData{T},
    dataset_id::String,
    observation_id::String
)::Set{QualityEvidenceEntry} where T
    if !haskey(vertex_data.evidence, dataset_id)
        return Set{QualityEvidenceEntry}()
    end

    if !haskey(vertex_data.evidence[dataset_id], observation_id)
        return Set{QualityEvidenceEntry}()
    end

    return vertex_data.evidence[dataset_id][observation_id]
end
```

**2. Iterator Return Types:**

```julia
# Type-stable iterator for coverage-only
function get_dataset_evidence(
    vertex_data::Union{KmerVertexData, BioSequenceVertexData, StringVertexData},
    dataset_id::String
)::Base.Iterators.Flatten{Base.Generator}
    if !haskey(vertex_data.evidence, dataset_id)
        # Return typed empty iterator
        return Iterators.flatten(Iterators.Empty{Tuple{String, Int, StrandOrientation}}())
    end

    return Iterators.flatten(
        ((obs_id, entry.position, entry.strand)
         for entry in evidence_set)
        for (obs_id, evidence_set) in vertex_data.evidence[dataset_id]
    )
end

# Quality-aware version explicitly typed
function get_dataset_evidence(
    vertex_data::Union{QualmerVertexData, QualityBioSequenceVertexData, QualityStringVertexData},
    dataset_id::String
)::Base.Iterators.Flatten{Base.Generator}
    if !haskey(vertex_data.evidence, dataset_id)
        return Iterators.flatten(
            Iterators.Empty{Tuple{String, Int, StrandOrientation, Vector{UInt8}}}()
        )
    end

    return Iterators.flatten(
        ((obs_id, entry.position, entry.strand, entry.quality_scores)
         for entry in evidence_set)
        for (obs_id, evidence_set) in vertex_data.evidence[dataset_id]
    )
end
```

**3. Type Parameterization Consistency:**

Use Union types for functions that handle multiple vertex types with same interface:

```julia
# Consistent across all coverage-only vertex types
const CoverageOnlyVertexData = Union{
    KmerVertexData,
    BioSequenceVertexData,
    StringVertexData
}

# Consistent across all quality-aware vertex types
const QualityAwareVertexData = Union{
    QualmerVertexData,
    QualityBioSequenceVertexData,
    QualityStringVertexData
}

# Then use in function signatures
function count_evidence(vertex_data::Union{CoverageOnlyVertexData, QualityAwareVertexData})
    return sum(
        length(evidence_set)
        for dataset_dict in values(vertex_data.evidence)
        for evidence_set in values(dataset_dict)
    )
end
```

**4. @inferred Testing:**

All evidence functions must pass `@inferred` tests:

```julia
Test.@testset "Type Stability - Evidence Functions" begin
    kmer_data = KmerVertexData(...)
    qualmer_data = QualmerVertexData(...)

    # Coverage-only functions
    Test.@inferred get_observation_evidence(kmer_data, "dataset_01", "obs_001")
    Test.@inferred get_dataset_evidence(kmer_data, "dataset_01")
    Test.@inferred count_evidence(kmer_data)

    # Quality-aware functions
    Test.@inferred get_observation_evidence(qualmer_data, "dataset_01", "obs_001")
    Test.@inferred get_dataset_evidence(qualmer_data, "dataset_01")
    Test.@inferred count_evidence(qualmer_data)
end
```

**Key Principles:**
- **Always specify return types** for public API functions
- **Use Union type aliases** for clarity and consistency
- **Empty containers must be typed** (`Set{EvidenceEntry}()` not `Set()`)
- **Test type stability** with `@inferred` for all critical paths
- **Multiple dispatch** handles different vertex types cleanly

---

#### 1.3.1 Supporting Helper Functions

**Utility Functions Referenced by Evidence Manipulation:**

These functions support the evidence manipulation API and provide core functionality for sequence operations,
reverse complement handling, and graph traversal.

```julia
# ============================================================================
# Sequence Access and Manipulation
# ============================================================================

"""
Get sequence from vertex data structure (works for all vertex types).
"""
function get_sequence(vertex_data::KmerVertexData{T}) where T
    return vertex_data.Kmer
end

function get_sequence(vertex_data::QualmerVertexData{T}) where T
    return vertex_data.Kmer
end

function get_sequence(vertex_data::BioSequenceVertexData{T}) where T
    return vertex_data.sequence
end

function get_sequence(vertex_data::QualityBioSequenceVertexData{T}) where T
    return vertex_data.sequence
end

function get_sequence(vertex_data::StringVertexData)
    return vertex_data.string_value
end

function get_sequence(vertex_data::QualityStringVertexData)
    return vertex_data.string_value
end

# ============================================================================
# Reverse Complement Operations
# ============================================================================

"""
Find all reverse-complement vertex pairs in the graph.
Returns vector of (forward_vertex, reverse_vertex) tuples.
"""
function find_reverse_complement_pairs(graph::MetaGraphsNext.MetaGraph)
    rc_pairs = []
    visited = Set()

    for vertex in vertices(graph)
        if vertex in visited
            continue
        end

        vertex_data = graph[vertex]
        seq = get_sequence(vertex_data)
        rc_seq = BioSequences.reverse_complement(seq)

        # Check if RC vertex exists
        if has_vertex(graph, rc_seq) && rc_seq != seq  # Don't pair palindromes with themselves
            rc_vertex = rc_seq
            push!(rc_pairs, (vertex, rc_vertex))
            push!(visited, vertex)
            push!(visited, rc_vertex)
        end
    end

    return rc_pairs
end

"""
Find all reverse-complement edge pairs in the graph.
Returns vector of ((src1, dst1), (src2, dst2)) tuples.
"""
function find_reverse_complement_edge_pairs(graph::MetaGraphsNext.MetaGraph)
    rc_edge_pairs = []
    visited = Set()

    for edge in edges(graph)
        if edge in visited
            continue
        end

        src, dst = edge
        src_data = graph[src]
        dst_data = graph[dst]

        # Find reverse complement vertices
        src_rc = BioSequences.reverse_complement(get_sequence(src_data))
        dst_rc = BioSequences.reverse_complement(get_sequence(dst_data))

        # RC edge goes from dst_rc to src_rc (reversed direction)
        rc_edge = (dst_rc, src_rc)

        if has_edge(graph, rc_edge...)
            push!(rc_edge_pairs, (edge, rc_edge))
            push!(visited, edge)
            push!(visited, rc_edge)
        end
    end

    return rc_edge_pairs
end

"""
Create reverse-complement vertex data from forward vertex.
Flips strand orientations in evidence.
"""
function create_reverse_complement_vertex(vertex_data::KmerVertexData{T}) where T
    rc_kmer = BioSequences.reverse_complement(vertex_data.Kmer)
    rc_evidence = flip_evidence_strands(vertex_data.evidence)
    return KmerVertexData{T}(rc_kmer, rc_evidence)
end

function create_reverse_complement_vertex(vertex_data::QualmerVertexData{T}) where T
    rc_kmer = BioSequences.reverse_complement(vertex_data.Kmer)
    rc_evidence = flip_quality_evidence_strands(vertex_data.evidence)
    return QualmerVertexData{T}(rc_kmer, rc_evidence)
end

"""
Flip strand orientations in evidence (Forward <-> Reverse).
Used when creating reverse-complement vertices.
"""
function flip_evidence_strands(
    evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}
)
    flipped_evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()

    for (dataset_id, dataset_dict) in evidence
        flipped_evidence[dataset_id] = Dict{String, Set{EvidenceEntry}}()

        for (obs_id, evidence_set) in dataset_dict
            flipped_evidence[dataset_id][obs_id] = Set{EvidenceEntry}()

            for entry in evidence_set
                # Flip strand orientation
                flipped_strand = entry.strand == Forward ? Reverse : Forward
                push!(flipped_evidence[dataset_id][obs_id],
                      EvidenceEntry(entry.position, flipped_strand))
            end
        end
    end

    return flipped_evidence
end

"""
Flip strand orientations in quality-aware evidence.
"""
function flip_quality_evidence_strands(
    evidence::Dict{String, Dict{String, Set{QualityEvidenceEntry}}}
)
    flipped_evidence = Dict{String, Dict{String, Set{QualityEvidenceEntry}}}()

    for (dataset_id, dataset_dict) in evidence
        flipped_evidence[dataset_id] = Dict{String, Set{QualityEvidenceEntry}}()

        for (obs_id, evidence_set) in dataset_dict
            flipped_evidence[dataset_id][obs_id] = Set{QualityEvidenceEntry}()

            for entry in evidence_set
                flipped_strand = entry.strand == Forward ? Reverse : Forward
                # Quality scores also get reversed for RC
                reversed_quality = reverse(entry.quality_scores)
                push!(flipped_evidence[dataset_id][obs_id],
                      QualityEvidenceEntry(entry.position, flipped_strand, reversed_quality))
            end
        end
    end

    return flipped_evidence
end

"""
Check if vertex represents reverse-complement orientation.
"""
function is_reverse_complement_orientation(
    vertex_data,
    target_strand::StrandOrientation
)
    # Check if majority of evidence is from opposite strand
    forward_count = 0
    reverse_count = 0

    for (dataset_id, dataset_dict) in vertex_data.evidence
        for (obs_id, evidence_set) in dataset_dict
            for entry in evidence_set
                if entry.strand == Forward
                    forward_count += 1
                else
                    reverse_count += 1
                end
            end
        end
    end

    # If target is Forward, return true if this is primarily Reverse
    if target_strand == Forward
        return reverse_count > forward_count
    else
        return forward_count > reverse_count
    end
end

"""
Replicate edges in reverse-complement orientation.
For each edge A→B, create RC edge: RC(B)→RC(A).
"""
function replicate_edges_reverse_complement!(graph::MetaGraphsNext.MetaGraph)
    edges_to_add = []

    for edge in edges(graph)
        src, dst = edge
        src_data = graph[src]
        dst_data = graph[dst]

        # Get RC sequences
        src_rc = BioSequences.reverse_complement(get_sequence(src_data))
        dst_rc = BioSequences.reverse_complement(get_sequence(dst_data))

        # RC edge goes in reverse direction: dst_rc → src_rc
        if !has_edge(graph, dst_rc, src_rc)
            # Get edge data and flip it
            edge_data = graph[src, dst]
            rc_edge_data = create_reverse_complement_edge_data(edge_data)
            push!(edges_to_add, (dst_rc, src_rc, rc_edge_data))
        end
    end

    for (src, dst, edge_data) in edges_to_add
        add_edge!(graph, src, dst, edge_data)
    end

    return graph
end

"""
Create reverse-complement edge data from forward edge.
"""
function create_reverse_complement_edge_data(edge_data::KmerEdgeData)
    rc_evidence = flip_edge_evidence_strands(edge_data.evidence)
    return KmerEdgeData(rc_evidence)
end

function create_reverse_complement_edge_data(edge_data::QualmerEdgeData)
    rc_evidence = flip_edge_quality_evidence_strands(edge_data.evidence)
    return QualmerEdgeData(rc_evidence)
end

"""
Flip strand orientations in edge evidence.
"""
function flip_edge_evidence_strands(
    evidence::Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}
)
    flipped = Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}()

    for (dataset_id, dataset_dict) in evidence
        flipped[dataset_id] = Dict{String, Set{EdgeEvidenceEntry}}()

        for (obs_id, evidence_set) in dataset_dict
            flipped[dataset_id][obs_id] = Set{EdgeEvidenceEntry}()

            for entry in evidence_set
                flipped_strand = entry.strand == Forward ? Reverse : Forward
                push!(flipped[dataset_id][obs_id],
                      EdgeEvidenceEntry(entry.from_position, entry.to_position, flipped_strand))
            end
        end
    end

    return flipped
end

"""
Flip strand orientations in quality-aware edge evidence.
"""
function flip_edge_quality_evidence_strands(
    evidence::Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}
)
    flipped = Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}()

    for (dataset_id, dataset_dict) in evidence
        flipped[dataset_id] = Dict{String, Set{EdgeQualityEvidenceEntry}}()

        for (obs_id, evidence_set) in dataset_dict
            flipped[dataset_id][obs_id] = Set{EdgeQualityEvidenceEntry}()

            for entry in evidence_set
                flipped_strand = entry.strand == Forward ? Reverse : Forward
                # Reverse quality scores for RC
                reversed_from_q = reverse(entry.from_quality)
                reversed_to_q = reverse(entry.to_quality)
                push!(flipped[dataset_id][obs_id],
                      EdgeQualityEvidenceEntry(entry.from_position, entry.to_position,
                                              flipped_strand, reversed_from_q, reversed_to_q))
            end
        end
    end

    return flipped
end

# ============================================================================
# Graph Traversal and Path Operations
# ============================================================================

"""
Find all unbranching paths (linear chains) in the graph.
An unbranching path is a maximal path where each internal vertex has exactly one in-edge and one out-edge.
"""
function find_unbranching_paths(graph::MetaGraphsNext.MetaGraph)
    paths = []
    visited = Set()

    for vertex in vertices(graph)
        if vertex in visited
            continue
        end

        # Start new path if this is a branch point or source
        in_degree = length(inneighbors(graph, vertex))
        out_degree = length(outneighbors(graph, vertex))

        if in_degree != 1 || out_degree == 0  # Source or branch point
            # Follow each outgoing edge
            for next_vertex in outneighbors(graph, vertex)
                path = [vertex, next_vertex]
                push!(visited, vertex)
                current = next_vertex

                # Extend path while unbranching
                while true
                    push!(visited, current)
                    out_neighbors = outneighbors(graph, current)
                    in_neighbors = inneighbors(graph, current)

                    # Stop if branching or sink
                    if length(out_neighbors) != 1 || length(in_neighbors) != 1
                        break
                    end

                    current = collect(out_neighbors)[1]
                    push!(path, current)
                end

                if length(path) > 1  # Only add non-trivial paths
                    push!(paths, path)
                end
            end
        end
    end

    return paths
end

"""
Concatenate sequences along a path.
For k-mer graphs: overlap by k-1
For variable-length graphs: use overlap_length from edge data
"""
function concatenate_path_sequences(path::Vector, graph::MetaGraphsNext.MetaGraph)
    if isempty(path)
        return ""
    end

    # Get first sequence
    first_vertex_data = graph[path[1]]
    result = get_sequence(first_vertex_data)

    # For k-mer graphs, detect k from first sequence length
    is_kmer_graph = typeof(first_vertex_data) <: Union{KmerVertexData, QualmerVertexData}

    if is_kmer_graph
        k = length(result)
        # Append last base of each subsequent k-mer
        for i in 2:length(path)
            vertex_data = graph[path[i]]
            kmer = get_sequence(vertex_data)
            result = result * kmer[end]  # Append last character
        end
    else
        # Variable-length: use overlap information from edges
        for i in 2:length(path)
            src = path[i-1]
            dst = path[i]
            edge_data = graph[src, dst]
            overlap_len = edge_data.overlap_length

            dst_vertex_data = graph[dst]
            dst_seq = get_sequence(dst_vertex_data)

            # Append non-overlapping portion
            result = result * dst_seq[(overlap_len+1):end]
        end
    end

    return result
end

"""
Build k-mer edges based on (k-1) overlap.
For each k-mer, add edge to all k-mers that overlap by k-1.
"""
function build_kmer_edges!(graph::MetaGraphsNext.MetaGraph, k::Int)
    # For each vertex, find all vertices that could follow it
    for src_vertex in vertices(graph)
        src_data = graph[src_vertex]
        src_kmer = get_sequence(src_data)

        # Suffix of length k-1
        suffix = src_kmer[2:end]

        # Find all k-mers whose prefix matches this suffix
        for dst_vertex in vertices(graph)
            if src_vertex == dst_vertex
                continue
            end

            dst_data = graph[dst_vertex]
            dst_kmer = get_sequence(dst_data)

            # Prefix of length k-1
            prefix = dst_kmer[1:end-1]

            if suffix == prefix
                # Add edge if it doesn't exist
                if !has_edge(graph, src_vertex, dst_vertex)
                    edge_data = typeof(src_data) <: QualmerVertexData ?
                        QualmerEdgeData(Dict{String, Dict{String, Set{EdgeQualityEvidenceEntry}}}()) :
                        KmerEdgeData(Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}())
                    add_edge!(graph, src_vertex, dst_vertex, edge_data)
                end
            end
        end
    end

    return graph
end

"""
Build edges for variable-length graphs based on sequence overlaps.

Strategy: Three-tiered approach from fast/exact to slow/sensitive

1. **Exact suffix-prefix matching** (fast, for high-quality assemblies)
   - Find exact overlaps between sequence ends
   - Minimum overlap length parameter (e.g., k-1 for consistency with k-mer graphs)

2. **Bounded-error matching** (medium speed, for noisy data)
   - Allow small number of mismatches in overlap
   - Use edit distance or Hamming distance threshold

3. **Full alignment** (slow, for maximum sensitivity - future work)
   - Use BioAlignments.jl for gapped alignment
   - Or integrate external aligners (minimap2, etc.)

Current implementation: Tier 1 (exact matching) with option to upgrade to Tier 2/3 later.
"""
function build_variable_length_edges!(
    graph::MetaGraphsNext.MetaGraph;
    min_overlap::Int=30,
    max_mismatches::Int=0,  # 0 = exact matching only
    dataset_id::String="default"
)
    vertices_list = collect(MetaGraphsNext.vertices(graph))

    # For each pair of vertices, check for suffix-prefix overlap
    for i in 1:length(vertices_list)
        for j in 1:length(vertices_list)
            if i == j
                continue  # Skip self-loops
            end

            src_vertex = vertices_list[i]
            dst_vertex = vertices_list[j]

            src_data = graph[src_vertex]
            dst_data = graph[dst_vertex]

            src_seq = get_sequence(src_data)
            dst_seq = get_sequence(dst_data)

            # Find overlap: suffix of src matches prefix of dst
            overlap_length = find_overlap(src_seq, dst_seq, min_overlap, max_mismatches)

            if overlap_length >= min_overlap
                # Create edge with overlap metadata
                edge_data = create_overlap_edge_data(
                    overlap_length,
                    src_data,
                    dst_data,
                    dataset_id
                )

                # Add edge if not already present
                if !MetaGraphsNext.has_edge(graph, src_vertex, dst_vertex)
                    MetaGraphsNext.add_edge!(graph, src_vertex, dst_vertex, edge_data)
                end
            end
        end
    end

    return graph
end

"""
Find overlap length between suffix of seq1 and prefix of seq2.

Returns: Length of longest overlap meeting criteria, or 0 if none found.
"""
function find_overlap(
    seq1::T,
    seq2::T,
    min_length::Int,
    max_mismatches::Int
) where T <: BioSequences.LongSequence
    max_overlap = min(length(seq1), length(seq2))

    # Try overlaps from longest to shortest
    for overlap_len in max_overlap:-1:min_length
        # Extract suffix of seq1 and prefix of seq2
        suffix_start = length(seq1) - overlap_len + 1
        suffix = seq1[suffix_start:end]
        prefix = seq2[1:overlap_len]

        # Count mismatches
        mismatches = count_mismatches(suffix, prefix)

        if mismatches <= max_mismatches
            return overlap_len
        end
    end

    return 0  # No valid overlap found
end

"""
Find overlap between string sequences (for StringVertexData).
"""
function find_overlap(
    seq1::String,
    seq2::String,
    min_length::Int,
    max_mismatches::Int
)
    max_overlap = min(length(seq1), length(seq2))

    for overlap_len in max_overlap:-1:min_length
        suffix_start = length(seq1) - overlap_len + 1
        suffix = seq1[suffix_start:end]
        prefix = seq2[1:overlap_len]

        mismatches = sum(suffix[i] != prefix[i] for i in 1:overlap_len)

        if mismatches <= max_mismatches
            return overlap_len
        end
    end

    return 0
end

"""
Count mismatches between two BioSequences.
"""
function count_mismatches(seq1::T, seq2::T) where T <: BioSequences.LongSequence
    Base.@assert length(seq1) == length(seq2) "Sequences must have equal length"

    mismatches = 0
    for i in 1:length(seq1)
        if seq1[i] != seq2[i]
            mismatches += 1
        end
    end

    return mismatches
end

"""
Create edge data for overlap-based edge.
"""
function create_overlap_edge_data(
    overlap_length::Int,
    src_vertex_data,
    dst_vertex_data,
    dataset_id::String
)
    # Create evidence from overlapping region
    # This is simplified - real implementation would track which observations
    # support this specific overlap

    evidence = Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}()

    # For now, create minimal edge with overlap length metadata
    # Future: extract actual observation support from vertex evidence

    edge_data = BioSequenceEdgeData(overlap_length, evidence)

    return edge_data
end

"""
Convert edges from quality-aware to quality-unaware.
"""
function convert_edges_remove_quality!(non_quality_graph, quality_graph)
    for edge in edges(quality_graph)
        src, dst = edge
        quality_edge_data = quality_graph[src, dst]

        # Convert quality evidence to non-quality
        non_quality_evidence = Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}()

        for (dataset_id, obs_dict) in quality_edge_data.evidence
            non_quality_evidence[dataset_id] = Dict{String, Set{EdgeEvidenceEntry}}()

            for (obs_id, quality_evidence_set) in obs_dict
                non_quality_evidence[dataset_id][obs_id] = Set{EdgeEvidenceEntry}()

                for quality_entry in quality_evidence_set
                    push!(non_quality_evidence[dataset_id][obs_id],
                          EdgeEvidenceEntry(quality_entry.from_position,
                                          quality_entry.to_position,
                                          quality_entry.strand))
                end
            end
        end

        # Add edge to non-quality graph
        edge_data = KmerEdgeData(non_quality_evidence)  # Adjust type based on graph type
        add_edge!(non_quality_graph, src, dst, edge_data)
    end

    return non_quality_graph
end
```

---

#### 1.4 Extract Graph Construction Logic (`src/rhizomorph/core/graph-construction.jl`)

**SingleStrand Construction Pattern:**
```julia
function build_singlestrand_graph(
    observations::Vector{<:FASTX.Record},
    vertex_data_type::Type{<:VertexData{T}},
    edge_data_type::Type
) where T
    # 1. Create graph with all observed elements
    graph = MetaGraphsNext.MetaGraph(...)

    # 2. Process each observation independently
    for (obs_idx, obs) in enumerate(observations)
        path = extract_path(obs, elements)
        add_observation_to_graph!(graph, path, obs_idx)
    end

    return graph
end
```

**Interconversion Between Graph Modes:**

The architecture supports bidirectional transformation between all four combinations of strand specificity and representation:

```julia
# ============================================================================
# Strand Specificity Transformations
# ============================================================================

"""
Convert strand-specific graph to non-strand-specific by merging evidence
across reverse-complement pairs.
"""
function convert_to_nonstrand_specific!(graph::MetaGraphsNext.MetaGraph)
    # Find all reverse-complement vertex pairs
    rc_pairs = find_reverse_complement_pairs(graph)

    # Merge evidence from RC pairs
    for (forward_v, reverse_v) in rc_pairs
        forward_data = graph[forward_v]
        reverse_data = graph[reverse_v]

        # Use helper function to merge evidence
        merge_reverse_complement_evidence!(forward_data, reverse_data)

        # Both vertices now have identical evidence pools
    end

    # Similarly for edges
    merge_reverse_complement_edge_evidence!(graph)

    return graph
end

"""
Convert non-strand-specific graph to strand-specific by filtering
evidence to retain only specified strand orientation.
"""
function convert_to_strand_specific!(
    graph::MetaGraphsNext.MetaGraph,
    target_strand::StrandOrientation = Forward
)
    # Filter all vertex evidence to target strand
    for vertex in vertices(graph)
        vertex_data = graph[vertex]
        filtered_data = filter_by_strand(vertex_data, target_strand)
        graph[vertex] = filtered_data
    end

    # Filter all edge evidence to target strand
    for edge in edges(graph)
        edge_data = graph[edge]
        filtered_data = filter_edge_by_strand(edge_data, target_strand)
        graph[edge] = filtered_data
    end

    return graph
end

# ============================================================================
# Strand Representation Transformations
# ============================================================================

"""
Convert single-stranded representation to double-stranded by replicating
vertices and edges in reverse-complement orientation.
"""
function convert_to_doublestrand_representation!(graph::MetaGraphsNext.MetaGraph)
    # Collect all vertices to duplicate (can't modify during iteration)
    vertices_to_add = []

    for vertex in vertices(graph)
        vertex_data = graph[vertex]
        rc_sequence = reverse_complement(get_sequence(vertex_data))

        # Create RC vertex if it doesn't exist
        if !has_vertex(graph, rc_sequence)
            rc_vertex_data = create_reverse_complement_vertex(vertex_data)
            push!(vertices_to_add, (rc_sequence, rc_vertex_data))
        end
    end

    # Add RC vertices
    for (rc_seq, rc_data) in vertices_to_add
        add_vertex!(graph, rc_seq, rc_data)
    end

    # Similarly replicate edges in RC orientation
    replicate_edges_reverse_complement!(graph)

    return graph
end

"""
Convert double-stranded representation to single-stranded by collapsing
to one strand orientation (keeping evidence from target strand only).
"""
function convert_to_singlestrand_representation!(
    graph::MetaGraphsNext.MetaGraph,
    target_strand::StrandOrientation = Forward
)
    vertices_to_remove = []

    for vertex in vertices(graph)
        vertex_data = graph[vertex]

        # If this is a reverse-strand vertex, mark for removal
        if is_reverse_complement_orientation(vertex_data, target_strand)
            push!(vertices_to_remove, vertex)
        end
    end

    # Remove RC vertices and their edges
    for vertex in vertices_to_remove
        rem_vertex!(graph, vertex)
    end

    return graph
end

# ============================================================================
# Combined Transformations
# ============================================================================

"""
Transform graph between any two modes.
"""
function transform_graph_mode!(
    graph::MetaGraphsNext.MetaGraph;
    target_specificity::Symbol = :strand_specific,  # :strand_specific or :nonstrand_specific
    target_representation::Symbol = :single_strand  # :single_strand or :double_strand
)
    # Handle specificity transformation
    if target_specificity == :nonstrand_specific
        convert_to_nonstrand_specific!(graph)
    else
        convert_to_strand_specific!(graph)
    end

    # Handle representation transformation
    if target_representation == :double_strand
        convert_to_doublestrand_representation!(graph)
    else
        convert_to_singlestrand_representation!(graph)
    end

    return graph
end
```

**Transformation Matrix:**

| From → To | Strand-Specific Single | Strand-Specific Double | Non-Specific Single | Non-Specific Double |
|-----------|------------------------|------------------------|---------------------|---------------------|
| **Strand-Specific Single** | - | Replicate RC | Merge RC evidence | Merge RC + Replicate |
| **Strand-Specific Double** | Filter to 1 strand | - | Merge RC evidence | Merge RC evidence |
| **Non-Specific Single** | Filter to 1 strand | Replicate RC | - | Replicate RC |
| **Non-Specific Double** | Filter to 1 strand | Filter evidence only | Collapse to 1 strand | - |

**Design Principles:**
- **Lossless where possible**: Strand-specific → Non-specific is **structurally lossless** (can reverse by filtering)
  - **Note on reversibility**: When merging Forward+Reverse evidence, you lose information about which strand each piece of evidence originally came from
  - However, you can filter back to strand-specific by selecting one orientation (e.g., keep only Forward evidence)
  - This is **structurally reversible** (same graph structure) but **informationally lossy** (lose original strand assignment)
  - Example: If kmer ATCG had 5 Forward observations and 3 Reverse, after merging you have 8 total but don't know the 5/3 split
  - After filtering back to Forward-only, you might have 4 Forward (lost 1) or some other subset
- **Lossy transformations explicit**: Single → Double (replication) and filtering are clearly documented as lossy
- **Evidence-first**: All transformations operate on evidence data structures
- **Validation**: Check that transformations maintain graph validity

#### 1.5 Graph Type Interconversion (`src/rhizomorph/core/graph-type-conversions.jl`)

**Additional Transformation Dimensions:**

Beyond strand specificity and representation, the architecture supports transformations between:
1. **Fixed-length ↔ Variable-length graphs**
2. **Quality-aware ↔ Quality-unaware graphs**

These transformations enable flexible analysis workflows where graphs can be converted between representations as needed.

```julia
# ============================================================================
# Fixed-Length ↔ Variable-Length Conversions
# ============================================================================

"""
Convert variable-length graph to fixed-length by fragmenting sequences into
fixed-length k-mers or n-grams.

Note: This is a lossy transformation - path information is lost.
"""
function fragment_to_fixed_length(
    variable_graph::MetaGraphsNext.MetaGraph,
    k::Int;
    graph_type::Symbol = :kmer  # :kmer, :qualmer, :ngram
)
    # Create new fixed-length graph
    fixed_graph = MetaGraphsNext.MetaGraph(...)

    # For each variable-length vertex
    for vertex in vertices(variable_graph)
        vertex_data = variable_graph[vertex]
        sequence = get_sequence(vertex_data)

        # Fragment into k-mers/n-grams
        for (i, kmer) in enumerate(each_kmer(sequence, k))
            # Transfer evidence from variable-length to fixed-length
            # Position in original sequence becomes position in evidence
            for (dataset_id, obs_dict) in vertex_data.evidence
                for (obs_id, evidence_set) in obs_dict
                    for evidence_entry in evidence_set
                        # Map position from variable-length vertex to k-mer position
                        new_position = evidence_entry.position + i - 1
                        add_evidence!(fixed_graph[kmer], dataset_id, obs_id,
                                    new_position, evidence_entry.strand)
                    end
                end
            end
        end
    end

    # Build edges from k-mer overlaps
    build_kmer_edges!(fixed_graph, k)

    return fixed_graph
end

"""
Convert fixed-length graph to variable-length by collapsing unbranching paths
into single vertices.

This is a lossless transformation - all evidence is preserved.
"""
function collapse_to_variable_length(
    fixed_graph::MetaGraphsNext.MetaGraph;
    graph_type::Symbol = :fasta  # :fasta, :fastq, :string
)
    variable_graph = MetaGraphsNext.MetaGraph(...)

    # Find all unbranching paths (linear chains)
    unbranching_paths = find_unbranching_paths(fixed_graph)

    for path in unbranching_paths
        # Concatenate sequences along path
        collapsed_sequence = concatenate_path_sequences(path, fixed_graph)

        # Merge evidence from all vertices in path
        merged_evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()

        for (vertex_idx, vertex) in enumerate(path)
            vertex_data = fixed_graph[vertex]

            # Transfer evidence with adjusted positions
            for (dataset_id, obs_dict) in vertex_data.evidence
                if !haskey(merged_evidence, dataset_id)
                    merged_evidence[dataset_id] = Dict{String, Set{EvidenceEntry}}()
                end

                for (obs_id, evidence_set) in obs_dict
                    if !haskey(merged_evidence[dataset_id], obs_id)
                        merged_evidence[dataset_id][obs_id] = Set{EvidenceEntry}()
                    end

                    for evidence_entry in evidence_set
                        # Adjust position based on vertex position in path
                        # vertex_idx is 1-based index in path
                        # For k-mer graphs: each k-mer advances position by 1
                        # Position in collapsed sequence = position_in_kmer + (vertex_idx - 1)
                        adjusted_position = evidence_entry.position + (vertex_idx - 1)

                        # Edge case validation
                        if adjusted_position < 1
                            @warn "Negative position adjustment: entry.position=$(evidence_entry.position), vertex_idx=$vertex_idx"
                            continue  # Skip invalid evidence
                        end

                        push!(merged_evidence[dataset_id][obs_id],
                              EvidenceEntry(adjusted_position, evidence_entry.strand))
                    end
                end
            end
        end

        # Create variable-length vertex
        if graph_type == :fasta
            new_vertex_data = BioSequenceVertexData(collapsed_sequence, merged_evidence)
        elseif graph_type == :fastq
            new_vertex_data = QualityBioSequenceVertexData(collapsed_sequence, merged_evidence)
        else  # :string
            new_vertex_data = StringVertexData(string(collapsed_sequence), merged_evidence)
        end

        add_vertex!(variable_graph, collapsed_sequence, new_vertex_data)
    end

    # Build edges from remaining branch points
    build_variable_length_edges!(variable_graph)

    return variable_graph
end

# ============================================================================
# Quality-Aware ↔ Quality-Unaware Conversions
# ============================================================================

"""
Remove quality information from quality-aware graph.

Note: This is a lossy transformation - quality scores are discarded.
User must explicitly request this transformation.
"""
function remove_quality_scores(
    quality_graph::MetaGraphsNext.MetaGraph;
    target_type::Symbol = :kmer  # :kmer, :fasta, :string
)
    # Create new graph without quality
    non_quality_graph = MetaGraphsNext.MetaGraph(...)

    for vertex in vertices(quality_graph)
        vertex_data = quality_graph[vertex]

        # Convert quality evidence to non-quality evidence
        non_quality_evidence = Dict{String, Dict{String, Set{EvidenceEntry}}}()

        for (dataset_id, obs_dict) in vertex_data.evidence
            non_quality_evidence[dataset_id] = Dict{String, Set{EvidenceEntry}}()

            for (obs_id, quality_evidence_set) in obs_dict
                non_quality_evidence[dataset_id][obs_id] = Set{EvidenceEntry}()

                for quality_entry in quality_evidence_set
                    # Drop quality_scores field
                    push!(non_quality_evidence[dataset_id][obs_id],
                          EvidenceEntry(quality_entry.position, quality_entry.strand))
                end
            end
        end

        # Create non-quality vertex data
        if target_type == :kmer
            new_vertex_data = KmerVertexData(vertex_data.Kmer, non_quality_evidence)
        elseif target_type == :fasta
            new_vertex_data = BioSequenceVertexData(vertex_data.sequence, non_quality_evidence)
        else  # :string
            new_vertex_data = StringVertexData(vertex_data.string_value, non_quality_evidence)
        end

        add_vertex!(non_quality_graph, get_sequence(vertex_data), new_vertex_data)
    end

    # Similarly convert edges
    convert_edges_remove_quality!(non_quality_graph, quality_graph)

    return non_quality_graph
end

"""
Note: Quality-unaware → Quality-aware conversion is NOT supported.
Quality information cannot be recovered once discarded.
Users should start from raw quality-aware data if quality is needed.
"""
```

**Transformation Summary:**

| From → To | Transformation | Lossless? | Notes |
|-----------|---------------|-----------|-------|
| **Variable → Fixed** | Fragment sequences | ❌ Lossy | Path structure lost; evidence preserved |
| **Fixed → Variable** | Collapse unbranching paths | ✅ Lossless | All evidence preserved with adjusted positions |
| **Quality-aware → Quality-unaware** | Drop quality scores | ❌ Lossy | Requires explicit user request |
| **Quality-unaware → Quality-aware** | ❌ Not supported | - | Quality cannot be recovered; use raw data |

**Combined Transformation Examples:**

```julia
# FASTQ (variable, quality) → Qualmer (fixed, quality)
qualmer_graph = fragment_to_fixed_length(fastq_graph, k=31, graph_type=:qualmer)

# Qualmer (fixed, quality) → Kmer (fixed, no quality)
kmer_graph = remove_quality_scores(qualmer_graph, target_type=:kmer)

# Kmer (fixed, no quality) → FASTA (variable, no quality)
fasta_graph = collapse_to_variable_length(kmer_graph, graph_type=:fasta)

# Full pipeline: FASTQ → Qualmer → Kmer → FASTA
fasta_from_fastq = fastq_graph |>
    g -> fragment_to_fixed_length(g, k=31, graph_type=:qualmer) |>
    g -> remove_quality_scores(g, target_type=:kmer) |>
    g -> collapse_to_variable_length(g, graph_type=:fasta)
```

**Design Principles:**
- **Explicit lossy transformations**: Quality removal and fragmentation require explicit user action
- **Lossless where possible**: Collapsing preserves all evidence with position adjustments
- **No reconstruction of lost data**: Cannot add quality scores or path structure after removal
- **Evidence tracking**: All transformations maintain dataset_id → observation_id hierarchy
- **Composable**: Transformations can be chained for complex workflows

---

### Phase 2: Fixed-Length Graphs (Week 2)

#### 2.1 Migrate K-mer Graphs (`src/rhizomorph/fixed-length/kmer-graphs.jl`)

**From:** `src/kmer-graphs.jl` + `sequence-graphs-next.jl` k-mer sections
**Actions:**
1. Move `build_kmer_graph_next()` implementation
2. **Critical Fix**: Remove canonicalization during construction for SingleStrand
3. Implement proper DoubleStrand as evidence merge across strands
4. Update to use new vertex/edge data structures
5. Add validation for k-mer type consistency

**Key Changes:**
```julia
function build_kmer_graph_next(
    kmer_type::Type{<:Kmers.Kmer},
    observations::Vector{<:FASTX.Record};
    graph_mode::GraphMode = SingleStrand
)
    graph = build_singlestrand_kmer_graph(kmer_type, observations)
    if graph_mode == DoubleStrand
        graph = convert_to_doublestrand!(graph)
    end
    return graph
end

function build_singlestrand_kmer_graph(kmer_type, observations)
    # Extract ALL observed k-mers (no canonicalization)
    all_kmers = extract_all_kmers(kmer_type, observations)
    return build_singlestrand_graph(all_kmers, observations, KmerVertexData{kmer_type}, KmerEdgeData)
end
```

#### 2.2 Migrate Qualmer Graphs (`src/rhizomorph/fixed-length/qualmer-graphs.jl`)

**From:** `src/qualmer-graphs.jl` + `src/qualmer-analysis.jl` bridge functions
**Actions:**
1. **Critical**: Import qualmer data structures from `qualmer-analysis.jl`
2. Create wrapper functions that call core qualmer implementations
3. Implement proper Phred score handling (not ASCII)
4. Add quality aggregation for multiple observations

**Key Implementation:**
```julia
function build_qualmer_graph(
    fastq_records::Vector{FASTX.FASTQ.Record};
    k::Int = 3,
    graph_mode::GraphMode = SingleStrand
)
    # Use existing qualmer counting from qualmer-analysis.jl
    qualmer_counts = if graph_mode == SingleStrand
        count_qualmers(fastq_records; k=k)
    else
        count_canonical_qualmers(fastq_records; k=k)
    end

    # Convert to proper Phred scores
    for (qualmer, count) in qualmer_counts
        qualmer.quality_scores = get_phred_scores(qualmer.quality_scores)
    end

    return build_graph_from_qualmers(qualmer_counts, fastq_records, graph_mode)
end
```

#### 2.3 Create N-gram Graphs (`src/rhizomorph/fixed-length/ngram-graphs.jl`)

**From:** `src/string-graphs.jl` n-gram functions
**Actions:**
1. Move fixed-length string n-gram functionality
2. Separate from variable-length string graphs
3. Implement quality-aware n-gram graphs
4. Implement evidence tracking with double-nested dictionary structure (currently broken)

---

### Phase 3: Variable-Length Graphs (Week 3)

#### 3.1 Consolidate FASTA Graphs (`src/rhizomorph/variable-length/fasta-graphs.jl`)

**From:** `src/fasta-graphs.jl` + BioSequence sections in `sequence-graphs-next.jl`
**Actions:**
1. Move all variable-length BioSequence graph construction
2. Implement proper alignment-based edge detection using exact matching, alignemnts with BioAlignments, and external tools (e.g. minimap2)
3. Implement evidence tracking with double-nested dictionary structure for variable-length vertices
4. Consider extending extension-based alphabet hints (`.fna`, `.frn`, `.faa`) to additional FASTA entry points package-wide using the shared `alphabet_hint_from_path` helper, while keeping `.fa`/`.fasta` ambiguous unless the caller opts in.

#### 3.2 Consolidate FASTQ Graphs (`src/rhizomorph/variable-length/fastq-graphs.jl`)

**From:** `src/fastq-graphs.jl` + quality BioSequence sections
**Actions:**
1. **Critical**: Fix Phred score handling throughout
2. Implement quality aggregation for overlapping sequences
3. Add provenance tracking for multiple observations
4. Implement proper alignment-based edge detection using exact matching, alignemnts with BioAlignments, and external tools (e.g. minimap2)

#### 3.3 Consolidate String Graphs (`src/rhizomorph/variable-length/string-graphs.jl`)

**From:** `src/string-graphs.jl` variable-length functions
**Actions:**
1. Move variable-length string graph functionality
2. Implement text-based alignment detection using exact matching, alignments julia tools (TBD) and external tools (TBD)
3. Add quality-aware variable-length string processing

---

### Phase 4: Algorithms & Utilities (Week 4)

#### 4.1 Path Finding (`src/rhizomorph/algorithms/path-finding.jl`)

**From:** Path reconstruction sections in `sequence-graphs-next.jl`
**Actions:**
1. Move `find_eulerian_paths_next()`
2. **Critical**: Fix `path_to_sequence()` type stability
3. Implement strand-aware path traversal
4. Add path validation for DoubleStrand graphs

#### 4.2 Canonicalization (`src/rhizomorph/algorithms/canonicalization.jl`)

**New Implementation:**
```julia
function canonicalize_doublestrand_graph!(graph::MetaGraphsNext.MetaGraph)
    # Find reverse complement pairs
    rc_pairs = find_reverse_complement_pairs(graph)

    # Merge vertices while preserving strand information
    for (forward_vertex, reverse_vertex) in rc_pairs
        merge_reverse_complement_vertices!(graph, forward_vertex, reverse_vertex)
    end

    return graph
end
```

#### 4.4 Error Correction and Graph Simplification Philosophy

**Philosophical Approach: Distinguishing Errors from True Variation**

A central challenge in sequence assembly is distinguishing sequencing errors from true rare genetic diversity. This challenge exists whether sequencing isolates (which may still contain low-frequency variation) or microbial communities (which definitely contain variation). There is no way to know with certainty whether observed variation represents valid low-frequency diversity or sequencing artifacts.

**Key Insight:** Statistical correction approaches can inadvertently remove valid low-frequency variation along with errors, similar to auto-correct replacing words that algorithms think are wrong based on context or frequency, even when they were exactly what you meant to write. The goal is to lean on quality scores and evidence frequencies to make more productive changes than unproductive ones, while acknowledging that the balance varies with sequencing depth, error rates, and biological context.

**Core Design Principles:**

1. **Probabilistic, Read-Centric Correction (Primary Strategy)**
   - Rather than deterministically modifying the graph structure (bubble-popping, tip-clipping)
   - Evaluate each read against the assembly graph individually
   - Accept updates probabilistically with normalized acceptance rate: P(alt) / (P(alt) + P(current))
     - This ensures proper probability range [0, 1]
     - High-likelihood alternatives accepted frequently but not deterministically
     - Low-likelihood alternatives occasionally accepted (stochastic sampling)
   - Expectation-maximization approach: sometimes make the change, sometimes don't
   - If low-frequency variation appears multiple times, some observations may be incorrectly modified while enough retain original sequence for variant to survive

2. **K-fold Cross-Validation with Random Seeds**
   - Multiple assembly runs with different random seeds
   - Compare whether variation is preserved across runs
   - Variation preserved in multiple assemblies has higher confidence

3. **Graph-Based Operations as Supplementary Tools**
   - Bubble popping, tip clipping, and variant reporting (à la FreeBayes/GATK)
   - Useful for visualization and validation
   - Should not replace probabilistic read-by-read correction as primary strategy

4. **Dual Support: Quality-Aware and Coverage-Only**
   - **Quality-aware**: Primary mode when FASTQ inputs available
   - **Coverage-only**: Fallback for FASTA inputs, pangenome graphs from reference assemblies
   - Both approaches must be robustly supported

#### 4.4.1 Read-Centric Probabilistic Correction Model (Primary Strategy)

**Workflow:**

1. **Initial Graph Construction**
   - Build graph from raw reads (quality-aware or coverage-only)
   - Preserve all observed evidence with strand, position, dataset, observation tracking

2. **Read-by-Read Evaluation**
   - For each read, find its current best path through the graph
   - Identify alternative paths (variants, potential corrections)
   - Calculate likelihood of current path vs alternative paths

3. **Probabilistic Path Selection**
   - For quality-aware graphs:
     - Calculate joint probability from Phred scores and inferred error rates
     - **Acceptance probability** = P(alt) / (P(alt) + P(current))
       - This normalizes to [0, 1] range
       - If P(alt) >> P(current): acceptance ≈ 1.0 (always accept)
       - If P(alt) << P(current): acceptance ≈ 0.0 (rarely accept)
       - If P(alt) = P(current): acceptance = 0.5 (50/50 chance)
     - Uses independence assumption for combining quality scores
   - For coverage-only graphs:
     - Use evidence depth and inferred error rates only
     - Acceptance probability = evidence(alt) / (evidence(alt) + evidence(current))
     - Ensures stochastic updates rather than deterministic

4. **Iterative Refinement**
   - Rebuild graph with updated read paths
   - Repeat correction rounds as needed
   - Monitor convergence and variant preservation

5. **Validation**
   - K-fold cross-validation with different random seeds
   - Track which variants persist across assemblies

6. **Parallelization Strategy**
   - Read-by-read correction is embarrassingly parallel
   - Each read evaluation is independent
   - Use Julia threading (`Threads.@threads`) for shared-memory parallelism
   - Use distributed computing for very large datasets

**Parallelization Implementation:**

```julia
"""
Parallel read-centric correction using Julia threading.

Distributes read evaluation across available CPU cores.
Safe because each iteration only reads from graph (no writes during parallel section).
"""
function probabilistic_read_correction_parallel!(
    graph::MetaGraphsNext.MetaGraph;
    error_rates::ErrorRateModel,
    n_iterations::Int=3,
    random_seed::Int=42
)
    # Get all observations from graph
    observations = collect_all_observations(graph)

    for iteration in 1:n_iterations
        # Thread-safe: each read processed independently
        # No writes to shared graph during parallel section
        corrections = Vector{Any}(undef, length(observations))

        Threads.@threads for i in 1:length(observations)
            dataset_id, obs_id = observations[i]

            # Evaluate this read's path
            current_path = find_read_path(graph, dataset_id, obs_id)
            alternative_paths = find_alternative_paths(graph, current_path)

            # Calculate likelihoods
            current_likelihood = calculate_path_likelihood(current_path, graph)
            best_alt = find_best_alternative(alternative_paths, graph)

            if !isnothing(best_alt)
                alt_likelihood = calculate_path_likelihood(best_alt, graph)
                acceptance_prob = alt_likelihood / (alt_likelihood + current_likelihood)

                # Stochastic acceptance
                if rand() < acceptance_prob
                    corrections[i] = (obs_id, dataset_id, best_alt)
                else
                    corrections[i] = nothing
                end
            else
                corrections[i] = nothing
            end
        end

        # Sequential graph updates (thread-safe)
        for correction in corrections
            if !isnothing(correction)
                obs_id, dataset_id, new_path = correction
                update_read_path!(graph, dataset_id, obs_id, new_path)
            end
        end

        # Rebuild graph with updated paths
        rebuild_graph!(graph)
    end

    return graph
end
```

**Distributed Computing (for very large datasets):**

```julia
# For datasets too large for single machine
using Distributed

"""
Distributed correction across multiple machines/nodes.

Partition observations across workers, collect results, merge graphs.
"""
function probabilistic_read_correction_distributed!(
    graph::MetaGraphsNext.MetaGraph;
    error_rates::ErrorRateModel,
    n_iterations::Int=3
)
    observations = collect_all_observations(graph)

    # Partition observations across workers
    chunks = partition_observations(observations, nworkers())

    for iteration in 1:n_iterations
        # Distribute work across nodes
        results = pmap(chunks) do chunk
            correct_reads_chunk(graph, chunk, error_rates)
        end

        # Merge results
        for result in results
            apply_corrections!(graph, result)
        end

        rebuild_graph!(graph)
    end

    return graph
end
```

**Performance Considerations:**

- **Threading overhead**: Only beneficial for >1000 reads
- **Memory sharing**: Graph is read-only during parallel evaluation phase
- **Lock-free design**: Updates applied sequentially after parallel section
- **Scalability**: Linear speedup expected up to number of available cores
- **Distributed**: For datasets >100GB, partition by chromosome/contig

**Quality Score Integration (Quality-Aware Mode):**

```julia
# Calculate joint probability of observation being correct
# Assuming independence (not Bayesian priors)

# Single observation: 90% confidence = Phred 10
# P(error) = 0.1, P(correct) = 0.9

# Two independent observations at 90% confidence each:
# P(both_errors) = 0.1 × 0.1 = 0.01
# P(correct) = 1 - 0.01 = 0.99 = 99% = Phred 20

# In Phred (log) space: additive rather than multiplicative
# Phred10 + Phred10 = Phred20

function calculate_path_likelihood_quality_aware(
    path_evidence::Vector{QualityEvidenceEntry},
    inferred_error_rates::ErrorRateModel
)
    # Combine Phred scores across observations (additive in log space)
    joint_quality = sum(entry.quality_scores for entry in path_evidence)

    # Weight by inferred insertion/deletion/substitution rates
    weighted_likelihood = apply_error_rate_weights(joint_quality, inferred_error_rates)

    return weighted_likelihood
end

"""
Statistical confidence adjustment for multiple observations.

When combining evidence from multiple independent observations, we need to adjust
our confidence to account for multiple testing. This is NOT traditional Bonferroni
correction (which adjusts significance levels), but rather confidence interval
adjustment for independent observations.

Mathematical approach:
- If each observation has error probability p_error
- For n independent observations, probability ALL are errors = p_error^n
- Therefore confidence increases: P(correct) = 1 - p_error^n

Example:
- Single observation: 99% confidence (p_error = 0.01)
- Two independent observations: P(both errors) = 0.01 × 0.01 = 0.0001
- Combined confidence: 1 - 0.0001 = 99.99%

For Phred scores (quality-aware graphs):
- This is handled by additive combination: Q_combined = Q1 + Q2 + ... + Qn
- No additional adjustment needed (already accounted for in aggregate_quality_scores)

For coverage-only graphs (no quality scores):
- Use conservative estimate based on evidence depth
- Higher coverage = higher confidence, but with diminishing returns
- Cap at reasonable maximum to avoid overconfidence

IMPORTANT: This function is for coverage-only mode. Quality-aware mode uses
direct Phred score addition which inherently accounts for multiple observations.
"""
function coverage_based_confidence(
    evidence_count::Int,
    expected_coverage::Float64,
    base_error_rate::Float64 = 0.01  # Default 1% error rate
)
    # Avoid division by zero
    if expected_coverage <= 0.0
        return 0.0
    end

    # Relative coverage depth
    coverage_ratio = evidence_count / expected_coverage

    # Conservative model: confidence saturates at high coverage
    # Uses logistic function to prevent overconfidence
    # At expected coverage: ~0.95 confidence
    # At 2× coverage: ~0.99 confidence
    # At 10× coverage: ~0.999 confidence (caps at asymptote)

    max_confidence = 1.0 - base_error_rate
    confidence = max_confidence * (1.0 - Base.exp(-coverage_ratio))

    return confidence
end
```

**Coverage-Only Mode (Quality-Unaware Graphs):**

When quality scores are unavailable (FASTA inputs, pangenome graphs), fall back to coverage-based likelihood:

```julia
"""
Calculate path likelihood for coverage-only graphs (no quality scores).

Uses evidence depth relative to expected coverage to estimate confidence.

Parameters:
- path_evidence: Evidence entries supporting this path
- expected_coverage: Expected/average coverage depth for this dataset
- inferred_error_rates: Error model inferred from graph structure

Returns:
- Normalized likelihood score [0.0, 1.0]

IMPORTANT: total_observations specification
- expected_coverage is the mean coverage depth across the genome/assembly
- For Illumina: typically 30-100× coverage
- For Nanopore/PacBio: typically 10-50× coverage
- Calculate from graph: mean(count_evidence(v) for v in vertices(graph))
- Or specify based on known sequencing depth

The ratio evidence_count/expected_coverage gives relative support:
- < 1.0: Below average support (potential error)
- = 1.0: Average support (typical variant)
- > 1.0: Above average support (high confidence)
"""
function calculate_path_likelihood_coverage_only(
    path_evidence::Vector{EvidenceEntry},
    expected_coverage::Float64,
    inferred_error_rates::ErrorRateModel
)
    # Count observations supporting this path
    evidence_count = length(path_evidence)

    # Calculate confidence based on coverage depth
    # Uses logistic model that saturates at high coverage
    base_confidence = coverage_based_confidence(evidence_count, expected_coverage,
                                                 inferred_error_rates.substitution_rate)

    # Weight by inferred error rates from graph structure
    # (Applies error type-specific adjustments: insertions vs deletions vs substitutions)
    weighted_likelihood = apply_error_rate_weights(base_confidence, inferred_error_rates)

    return weighted_likelihood
end
```

#### 4.4.2 Inferred Error Rate Modeling

When gold-standard validation sequences are unavailable, estimate insertion/deletion/substitution rates from graph structure:

**Approach:**

1. **Identify Consensus/Highway Paths**
   - Paths with highest evidence depth and confidence
   - "Main roads" through the graph

2. **Identify Alternative/Back-Road Paths**
   - Lower evidence, lower confidence
   - Variants, potential errors

3. **Calculate Error Rate Estimates**
   - Compare alternate path frequencies to consensus paths
   - Longer alternate paths (insertions) vs shorter (deletions) vs equal-length (substitutions)
   - Estimate rates from relative frequencies

```julia
struct InferredErrorRates
    insertion_rate::Float64
    deletion_rate::Float64
    substitution_rate::Float64
    confidence::Float64  # How confident are we in these estimates?
end

function infer_error_rates_from_graph(graph::MetaGraphsNext.MetaGraph)
    # Find highest-evidence paths (consensus/highways)
    consensus_paths = find_consensus_paths(graph)

    # Find alternative paths (variants/back-roads)
    alternative_paths = find_alternative_paths(graph, consensus_paths)

    # Classify alternatives by type
    insertions = filter(p -> length(p) > consensus_length(p), alternative_paths)
    deletions = filter(p -> length(p) < consensus_length(p), alternative_paths)
    substitutions = filter(p -> length(p) == consensus_length(p), alternative_paths)

    # Calculate rates from frequencies
    total_alts = length(alternative_paths)
    insertion_rate = length(insertions) / total_alts
    deletion_rate = length(deletions) / total_alts
    substitution_rate = length(substitutions) / total_alts

    # Estimate confidence in these rates
    confidence = calculate_rate_confidence(consensus_paths, alternative_paths)

    return InferredErrorRates(insertion_rate, deletion_rate, substitution_rate, confidence)
end

"""
Find consensus paths (highest evidence/confidence paths) in the graph.
Returns paths with highest weighted support.
"""
function find_consensus_paths(graph::MetaGraphsNext.MetaGraph; min_evidence::Int=10)
    consensus_paths = []

    # Find all unbranching paths
    paths = find_unbranching_paths(graph)

    # Score each path by evidence
    for path in paths
        total_evidence = 0
        min_vertex_evidence = typemax(Int)

        for vertex in path
            vertex_data = graph[vertex]
            vertex_evidence = count_evidence(vertex_data)
            total_evidence += vertex_evidence
            min_vertex_evidence = min(min_vertex_evidence, vertex_evidence)
        end

        # Require minimum evidence at each vertex
        if min_vertex_evidence >= min_evidence
            push!(consensus_paths, (path=path, evidence=total_evidence, min_evidence=min_vertex_evidence))
        end
    end

    # Sort by total evidence (highest first)
    sort!(consensus_paths, by=x->x.evidence, rev=true)

    return consensus_paths
end

"""
Find alternative paths (lower evidence, potential variants or errors).
"""
function find_alternative_paths(graph::MetaGraphsNext.MetaGraph, consensus_paths)
    all_paths = find_unbranching_paths(graph)

    # Mark consensus path vertices
    consensus_vertices = Set()
    for cons_path in consensus_paths
        for vertex in cons_path.path
            push!(consensus_vertices, vertex)
        end
    end

    # Find paths not in consensus
    alternative_paths = []
    for path in all_paths
        # Check if this path overlaps with consensus
        is_alternative = false
        for vertex in path
            if vertex ∉ consensus_vertices
                is_alternative = true
                break
            end
        end

        if is_alternative
            # Calculate evidence
            total_evidence = sum(count_evidence(graph[v]) for v in path)
            push!(alternative_paths, (path=path, evidence=total_evidence))
        end
    end

    return alternative_paths
end

"""
Get consensus length for a path (for classification as insertion/deletion/substitution).
"""
function consensus_length(path_info)
    return length(path_info.path)
end

"""
Calculate confidence in inferred error rates.
Based on amount of data available and variance in estimates.
"""
function calculate_rate_confidence(consensus_paths, alternative_paths)
    # More data = higher confidence
    total_paths = length(consensus_paths) + length(alternative_paths)
    total_evidence = sum(p.evidence for p in consensus_paths) + sum(p.evidence for p in alternative_paths)

    # Base confidence on total evidence
    if total_evidence < 100
        return 0.1  # Very low confidence
    elseif total_evidence < 1000
        return 0.3  # Low confidence
    elseif total_evidence < 10000
        return 0.5  # Moderate confidence
    elseif total_evidence < 100000
        return 0.7  # Good confidence
    else
        return 0.9  # High confidence
    end
end

"""
Apply error rate weights to likelihood calculation.
Adjusts path likelihood based on expected insertion/deletion/substitution rates.
"""
function apply_error_rate_weights(
    base_likelihood::Float64,
    error_rates::ErrorRateModel,
    path_type::Symbol  # :insertion, :deletion, :substitution, :match
)
    if path_type == :insertion
        # Penalize by insertion rate
        return base_likelihood * (1.0 - error_rates.insertion_rate)
    elseif path_type == :deletion
        # Penalize by deletion rate
        return base_likelihood * (1.0 - error_rates.deletion_rate)
    elseif path_type == :substitution
        # Penalize by substitution rate
        return base_likelihood * (1.0 - error_rates.substitution_rate)
    else  # :match
        # Reward matches (no error)
        error_rate_total = error_rates.insertion_rate + error_rates.deletion_rate + error_rates.substitution_rate
        return base_likelihood * (1.0 - error_rate_total)
    end
end
```

**Future Enhancement: Pre-trained Error Models**

For high-throughput or clinical workflows:

1. Run gold-standard sequences (known content, wide %GC range, homopolymers, etc.)
2. Through exact library prep and sequencing workflow
3. Observe actual insertion/deletion/substitution rates
4. Train model capturing workflow-specific error patterns
5. Apply trained model weights to unknown samples

This provides validated, context-matched error rates for maximum accuracy.

#### 4.4.3 Graph-Based Simplification Algorithms (Supplementary Tools)

These algorithms support visualization, validation, and variant calling, but are not the primary correction strategy.

##### 4.4.3.1 Tip Removal (Clipping)

Prunes short, dead-end paths (tips) that typically result from sequencing errors or coverage drop-off.

*   **Functionality:** The tip clipping function will accept thresholds based on:
    *   `min_length::Int`: The minimum length of a tip to be retained.
    *   `min_evidence::Float64`: The minimum average evidence depth of a tip to be retained.
    *   `min_confidence::Float64`: For quality-aware graphs, the minimum average vertex confidence to be retained.
*   **Logic:** Users can specify the relationship between these thresholds using `and`, `or`, or `either` logic to determine which tips are pruned.
*   **Default Thresholds:**
    *   **Length:** For fixed-length graphs, the default will be the k-mer/n-gram size. For variable-length graphs, this will not be set by default.
    *   **Evidence/Confidence:** The default threshold will be calculated statistically from the graph properties. It will be defined as a set number of standard deviations below the mean evidence depth/confidence of all vertices. An alternative based on the median and coefficient of variation will also be explored for skewed distributions.

**Important:** Tip clipping can remove valid low-frequency variants. Prefer read-centric correction where possible.

##### 4.4.3.2 Bubble Popping and Variant Calling

Identifies divergent paths (bubbles) representing SNPs, indels, or sequencing errors.

*   **Approach:** Bubbles will be identified as simple alternative paths that start and end at the same vertices.
*   **Context-Awareness:** The algorithm will be aware of the expected bubble size based on the graph type. A SNP in a fixed-length graph creates a bubble of length `k`, whereas in a variable-length graph, the bubble's length corresponds to the variant region.
*   **Two Modes:**
    1. **Bubble Popping:** Remove minor path if evidence suggests sequencing error
    2. **Variant Calling:** Report bubble as potential variant (à la FreeBayes, GATK)

*   **Probabilistic Model:** Drawing inspiration from FreeBayes and GATK:
    - Calculate posterior probability that minor path represents true variant vs artifact
    - Evaluate evidence: depth, quality scores, strand bias, observation support
    - Report variant with confidence scores OR pop bubble if error probability exceeds threshold

**Data Structures:**

```julia
"""
Represents a bubble (alternative paths) in the assembly graph.
"""
struct Bubble
    start_vertex::Any  # Parameterized by vertex type (k-mer, sequence, string)
    end_vertex::Any
    paths::Vector{Vector{Any}}  # Vector of paths, each path is vector of vertices
    evidence_counts::Vector{Int}  # Evidence count for each path
end

"""
Represents a called variant from bubble analysis.
"""
struct VariantCall
    bubble::Bubble
    variant_probability::Float64  # P(true variant | evidence)
    major_path::Vector{Any}  # Reference/consensus path
    minor_path::Vector{Any}  # Alternative/variant path
    variant_type::Symbol  # :snp, :insertion, :deletion, :complex
    quality_score::Float64  # Phred-scaled quality
    strand_bias::Float64  # Evidence of strand bias (0.0 = no bias, 1.0 = complete bias)
    observations_supporting::Set{Tuple{String, String}}  # (dataset_id, observation_id) pairs
end
```

**Bubble Evaluation:**

```julia
function evaluate_bubble(
    bubble::Bubble,
    graph::MetaGraphsNext.MetaGraph;
    mode::Symbol = :variant_calling  # :variant_calling or :bubble_popping
)
    major_path = bubble.paths[1]
    minor_path = bubble.paths[2]

    # Calculate likelihood of each path
    major_likelihood = calculate_path_likelihood(major_path, graph)
    minor_likelihood = calculate_path_likelihood(minor_path, graph)

    # Posterior probability minor path is true variant
    p_variant = minor_likelihood / (major_likelihood + minor_likelihood)

    if mode == :variant_calling
        # Report as variant with confidence
        return VariantCall(bubble, p_variant, major_path, minor_path)
    else  # :bubble_popping
        # Pop bubble if likely error
        if p_variant < error_threshold
            remove_path!(graph, minor_path)
            return :popped
        else
            return :retained
        end
    end
end
```

**Recommendation:** Use `:variant_calling` mode to preserve variation for downstream analysis. Only use `:bubble_popping` when over-assembly is problematic and read-centric correction is insufficient.

##### 4.4.3.3 Cycle and Repeat Resolution

Addresses simple cycles and repeats in the graph.

*   **Primary Strategy:** Leverage read-centric correction model
*   **Graph-Based Supplement:** For simple, short tandem repeats:
    - Assess likelihood of repeat copy number given spanning read evidence
    - Use quality scores when available
    - Report uncertainty rather than forcing resolution

#### 4.4.4 Future Direction: Reinforcement Learning

**Long-term enhancement (post-publication):**

The most promising approach to maximize accuracy vs runtime balance without heuristic thresholds:

1. **Data Collection Phase**
   - Generate large, diverse in-silico datasets (varying error rates, coverage, GC content, etc.)
   - Record all states, actions, and outcomes during assembly
   - States: graph configuration, evidence distribution, quality scores
   - Actions: correct this read, increment k, rebuild graph, run another correction round
   - Outcomes: final assembly accuracy, runtime

2. **Model Training**
   - Train RL model on state-action-outcome dataset
   - Learn optimal reward function balancing accuracy and efficiency
   - Determine which steps to take at each decision point

3. **Deployment**
   - Use trained model to guide assembly decisions
   - Avoid hard-coded heuristic thresholds
   - Adapt to different data types and quality profiles

**Current Status:** Not pursued until core implementation validated and published. Mentioned here as architectural consideration.

#### 4.5 I/O Functions (`src/rhizomorph/algorithms/io.jl`)

**I/O Strategy: Lossless JLD2 + Lossy GFA**

The Rhizomorph graph ecosystem uses a dual I/O strategy that balances complete data preservation with community interoperability:

- **Primary format: JLD2** (lossless) - Preserves all evidence, quality scores, and metadata
- **Secondary format: GFA** (lossy) - Standard format for visualization and tool interoperability

**Primary I/O: JLD2 (Lossless)**

Use JLD2 for full preservation of all graph data:

```julia
"""
Save graph to JLD2 file with complete preservation of all data structures.

Preserves:
- All vertex and edge data
- Double-nested evidence dictionaries
- Dataset and observation IDs
- Quality scores
- Strand orientations
- Graph mode metadata
"""
function save_graph(graph::MetaGraphsNext.MetaGraph, filepath::String)
    # Extract metadata
    metadata = Dict{String, Any}(
        "graph_type" => typeof(first(MetaGraphsNext.labels(graph))),
        "vertex_data_type" => typeof(graph[first(MetaGraphsNext.vertices(graph))]),
        "created" => Dates.now(),
        "rhizomorph_version" => "0.1.0"
    )

    JLD2.jldopen(filepath, "w") do file
        file["graph"] = graph
        file["metadata"] = metadata
    end
end

"""
Load graph from JLD2 file.
Completely lossless - reconstructs exact graph state.
"""
function load_graph(filepath::String)
    JLD2.jldopen(filepath, "r") do file
        graph = file["graph"]
        metadata = file["metadata"]
        return graph, metadata
    end
end

# Usage:
save_graph(my_graph, "assembly.jld2")
loaded_graph, metadata = load_graph("assembly.jld2")
```

**Secondary I/O: GFA (Lossy, Community Standard)**

Use GFA for interoperability with other tools:

```julia
"""
Export graph to GFA format (lossy).

What's preserved:
- Sequences (vertices)
- Connections (edges)
- Total evidence counts
- Basic metadata

What's lost:
- Individual observation IDs
- Position-specific evidence
- Detailed quality scores
- Full dataset provenance
"""
function export_to_gfa(graph::MetaGraphsNext.MetaGraph, filepath::String;
                       include_evidence_summary::Bool=true)
    Base.open(filepath, "w") do io
        # Write header with Rhizomorph-specific tags
        Base.println(io, "H\tVN:Z:1.0")

        # Write segments (vertices)
        for vertex in MetaGraphsNext.vertices(graph)
            vertex_data = graph[vertex]
            seq = get_sequence(vertex_data)

            # Basic GFA: S <name> <sequence>
            # Optional tags: RC:i:<evidence_count>
            if include_evidence_summary
                evidence_count = count_evidence(vertex_data)
                Base.println(io, "S\t$(seq)\t$(seq)\tRC:i:$(evidence_count)")
            else
                Base.println(io, "S\t$(seq)\t$(seq)")
            end
        end

        # Write links (edges)
        for edge in MetaGraphsNext.edges(graph)
            src, dst = edge
            # GFA: L <from> <from_orient> <to> <to_orient> <overlap>
            # For k-mer graphs, overlap is k-1
            Base.println(io, "L\t$(src)\t+\t$(dst)\t+\t*")
        end
    end
end

"""
Import graph from GFA format (lossy).

Reconstructs:
- Graph structure
- Sequences
- Evidence counts (if available in tags)

Cannot reconstruct:
- Individual observations
- Detailed provenance
- Quality scores
- Graph mode (assumes strand-specific, single-strand)
"""
function import_from_gfa(filepath::String; k::Union{Int,Nothing}=nothing)
    # Parse GFA file
    # Reconstruct basic graph structure
    # Return graph with minimal evidence (counts only if available)

    @warn "GFA import is lossy - detailed evidence and quality information not preserved"
    # Implementation...
end
```

**I/O Decision Matrix:**

| **Use Case** | **Format** | **Rationale** |
|--------------|------------|---------------|
| Save working graph | JLD2 | Preserve all data for continued analysis |
| Share with collaborators (Mycelia users) | JLD2 | Full data preservation |
| Archive final assembly | JLD2 | Complete record |
| Interoperability with Bandage, vg, etc. | GFA | Community standard |
| Publish assembly graph | GFA | Widely supported |
| Quick visualization | GFA | Many tools support it |

**Recommended Workflow:**

```julia
# During analysis: use JLD2
graph = build_kmer_graph(reads, k=31)
save_graph(graph, "intermediate_assembly.jld2")

# Later: resume work
graph, metadata = load_graph("intermediate_assembly.jld2")
corrected = probabilistic_read_correction!(graph)
save_graph(corrected, "corrected_assembly.jld2")

# Final: export for publication/sharing
export_to_gfa(corrected, "final_assembly.gfa")  # For community
save_graph(corrected, "final_assembly.jld2")     # For archive
```

**Implementation Actions:**
1. Implement save_graph() and load_graph() for JLD2 format
2. Implement export_to_gfa() with evidence summary tags
3. Implement import_from_gfa() with lossy data warning
4. Add comprehensive documentation explaining lossless vs lossy tradeoffs
5. Add tests verifying JLD2 round-trip preservation
6. Add tests verifying GFA compatibility with standard tools (Bandage, etc.)

---

### Phase 5: Integration & Testing (Week 5)

#### 5.1 Update Main Module Integration

**Update `src/Mycelia.jl`:**
```julia
# Include rhizomorph module
include("rhizomorph/rhizomorph.jl")
```

#### 5.2 Update 24 Test Files

**For each test file:**
1. Update import statements to use rhizomorph functions
2. Add proper Phred score validation
3. Test both SingleStrand and DoubleStrand modes
4. Validate evidence tracking and provenance with double-nested dictionary structure

#### 5.2.1 Error Handling and Input Validation Strategy

**Input Validation Patterns:**

```julia
function build_kmer_graph(records; k::Int=31)
    # Validate k
    if k < 1
        Base.error("k must be positive, got k=$k")
    end

    # Validate records
    if Base.isempty(records)
        Base.error("Cannot build graph from empty record set")
    end

    # Warn about short sequences
    short_reads = Base.count(r -> Base.length(FASTX.sequence(r)) < k, records)
    if short_reads > 0
        @warn "$short_reads sequences shorter than k=$k will be skipped"
    end

    # Continue...
end
```

**Edge Case Handling Matrix:**

| **Condition** | **Handling** | **Rationale** |
|--------------|--------------|---------------|
| k > sequence length | Skip with warning | Common in real data, not fatal |
| Empty datasets | Error | Cannot build graph from nothing |
| Duplicate dataset IDs on merge | Error | Ambiguous provenance |
| Negative position adjustment | Warn and skip | Data corruption or edge case |
| Quality > 255 after aggregation | Clamp to 255 | UInt8 range limit |
| Palindromic k-mers (self-RC) | Handle correctly | Valid biological sequence |
| Missing quality scores in FASTQ | Error | Data format violation |
| Inconsistent quality score length | Error | Data corruption |

**Error vs Warning Criteria:**

- **Error (abort operation):**
  - Invalid input types or parameters
  - Empty or corrupted data files
  - Logical impossibilities (negative k, etc.)
  - Data format violations

- **Warning (continue with caution):**
  - Skippable data (short reads, low quality)
  - Potentially suspicious patterns (high error rates, strand bias)
  - Non-critical information loss (GFA export)
  - Performance concerns (very large graphs)

#### 5.2.2 Memory Budget Guidance and Performance Monitoring

**When to Use Memory Estimation:**

The Mycelia package includes memory estimation utilities. Use these proactively to avoid out-of-memory errors:

```julia
"""
Estimate memory requirements before building graph.

Returns: Estimated memory in bytes
"""
function estimate_graph_memory(
    n_unique_kmers::Int,
    avg_evidence_per_kmer::Float64,
    has_quality_scores::Bool=false
)
    # Base vertex data size
    bytes_per_vertex = has_quality_scores ? 512 : 256  # Rough estimates

    # Evidence overhead (double-nested dictionaries)
    bytes_per_evidence_entry = has_quality_scores ? 64 : 32

    total_evidence_entries = n_unique_kmers * avg_evidence_per_kmer

    estimated_bytes = (n_unique_kmers * bytes_per_vertex) +
                      (total_evidence_entries * bytes_per_evidence_entry)

    return estimated_bytes
end

# Usage before construction
estimated_memory = estimate_graph_memory(10_000_000, 30.0, true)
available_memory = Sys.free_memory()

if estimated_memory > 0.8 * available_memory
    @warn "Graph may exceed available memory. Consider:" *
          "\n  - Reducing k-mer size" *
          "\n  - Filtering low-coverage k-mers" *
          "\n  - Using distributed computing" *
          "\n  - Processing in chunks"
end
```

**Graph Size Thresholds:**

| **Scenario** | **Action** | **Threshold** |
|--------------|------------|---------------|
| Small graph (< 1M k-mers) | Single-threaded, in-memory | < 10 GB RAM |
| Medium graph (1M - 10M k-mers) | Multi-threaded, in-memory | 10-100 GB RAM |
| Large graph (10M - 100M k-mers) | Distributed or chunked | 100 GB - 1 TB RAM |
| Very large graph (> 100M k-mers) | Distributed + disk spilling | > 1 TB RAM |

**Memory Monitoring During Execution:**

```julia
function build_kmer_graph_with_monitoring(records; k::Int=31)
    start_memory = Sys.free_memory()

    @info "Starting graph construction" free_memory_gb=start_memory/1e9

    graph = build_kmer_graph(records, k=k)

    end_memory = Sys.free_memory()
    memory_used = start_memory - end_memory

    @info "Graph construction complete" memory_used_gb=memory_used/1e9 \
          vertices=nv(graph) edges=ne(graph)

    # Warn if memory usage is high
    if memory_used > 0.7 * Sys.total_memory()
        @warn "High memory usage - consider distributed processing for larger datasets"
    end

    return graph
end
```

**Disk Spilling Strategy (Future Work):**

For graphs too large for RAM:
1. **Partition by k-mer prefix** (first 8-10 bases)
2. **Process each partition** independently
3. **Merge partitions** for global analysis
4. **Use memory-mapped I/O** for intermediate results

This is not implemented in initial version but architecture supports future extension.

#### 5.3 Comprehensive Validation and Testing Strategy

**Test Organization:**

```
test/rhizomorph/
├── runtests.jl              # Main test runner
├── test_core.jl             # Core data structures
├── test_evidence.jl         # Evidence manipulation
├── test_quality.jl          # Quality score math
├── test_graphs.jl           # Graph construction (all types)
├── test_transformations.jl  # Graph mode conversions
├── test_algorithms.jl       # Error correction, simplification
├── test_io.jl               # JLD2 and GFA I/O
└── test_integration.jl      # End-to-end workflows
```

**Key Test Categories:**

**1. Type Stability Tests:**
```julia
Test.@testset "Type Stability" begin
    graph = build_test_graph()
    vertex = first(MetaGraphsNext.vertices(graph))
    vertex_data = graph[vertex]

    # All evidence functions should be type-stable
    Test.@inferred get_observation_evidence(vertex_data, "dataset_01", "read_001")
    Test.@inferred count_evidence(vertex_data)
    Test.@inferred filter_by_strand(vertex_data, Forward)
end
```

**2. Quality Score Math Validation:**
```julia
Test.@testset "Phred Score Independence Assumption" begin
    # Q1 + Q2 should equal combined quality
    q1 = UInt8(10)  # 90% correct, 10% error
    q2 = UInt8(10)  # 90% correct, 10% error

    quals = [[q1 for _ in 1:31], [q2 for _ in 1:31]]
    combined = aggregate_quality_scores_independence(quals)

    # Should be 20 (0.1 * 0.1 = 0.01 error = Phred 20)
    Test.@test all(q == UInt8(20) for q in combined)
end
```

**3. Simulated Data Validation:**
```julia
Test.@testset "Error Correction Validation" begin
    # Generate known sequence
    true_seq = generate_random_sequence(1000)

    # Simulate reads with errors
    reads = simulate_reads(true_seq, coverage=50, error_rate=0.01)

    # Build and correct
    graph = build_kmer_graph(reads, k=31)
    corrected = probabilistic_read_correction!(graph)
    assembled = assemble_contigs(corrected)

    # Validate accuracy
    identity = sequence_identity(assembled[1], true_seq)
    Test.@test identity > 0.99
end
```

**4. Evidence Tracking and Provenance:**
```julia
Test.@testset "Double-Nested Dictionary Evidence" begin
    # Test dataset-level queries
    vertex_data = build_test_vertex_data()
    dataset_evidence = collect(get_dataset_evidence(vertex_data, "dataset_01"))
    Test.@test length(dataset_evidence) > 0

    # Test observation-level queries
    obs_evidence = get_observation_evidence(vertex_data, "dataset_01", "read_001")
    Test.@test !isempty(obs_evidence)

    # Test evidence counting
    Test.@test count_evidence(vertex_data) == expected_count
    Test.@test count_observations(vertex_data) == expected_observations
    Test.@test count_datasets(vertex_data) == expected_datasets
end
```

**5. I/O Round-Trip Tests:**
```julia
Test.@testset "JLD2 Lossless Round-Trip" begin
    original_graph = build_test_graph()

    # Save and load
    save_graph(original_graph, "test.jld2")
    loaded_graph, metadata = load_graph("test.jld2")

    # Verify complete preservation
    Test.@test graphs_equal(original_graph, loaded_graph)
    Test.@test all_evidence_preserved(original_graph, loaded_graph)
end

Test.@testset "GFA Lossy Export" begin
    graph = build_test_graph()
    export_to_gfa(graph, "test.gfa", include_evidence_summary=true)

    # Verify GFA format compliance
    Test.@test gfa_is_valid("test.gfa")

    # Verify expected information loss
    imported = import_from_gfa("test.gfa")
    Test.@test structure_preserved(graph, imported)
    Test.@test !detailed_evidence_preserved(graph, imported)  # Expected loss
end
```

**6. Performance Benchmarks:**

```julia
# Memory targets: 2-3× reduction from double-nested dictionaries
BenchmarkTools.@benchmark build_kmer_graph($reads, k=31)

# Speed targets: ≤110% of baseline
BenchmarkTools.@benchmark get_dataset_evidence($vertex_data, "dataset_01")

# Type stability: Zero type instability warnings
Test.@inferred get_observation_evidence(vertex_data, "dataset_id", "obs_id")
```

**Run test suite progression:**
1. Basic construction tests for each graph type (all 6 vertex types)
2. Type stability validation across all operations (@inferred tests)
3. SingleStrand vs DoubleStrand correctness
4. Evidence tracking and provenance with dataset/observation hierarchy
5. Quality score math validation (Phred conversions, aggregation)
6. Path reconstruction and sequence identity
7. I/O round-trip preservation (JLD2 lossless, GFA lossy)
8. Performance benchmarks vs current implementation
9. Integration tests (end-to-end workflows)
10. Simulated data validation with known ground truth

---

## Critical Issues to Address

### 1. **Immediate Phred Score Fix**
```julia
# WRONG (current)
quality_chars = FASTX.quality(record)
quality_ints = [Int(c) for c in quality_chars]

# RIGHT (target)
quality_phred = get_phred_scores(record)  # Returns 0-60 range
# also right and already exists
quality_phred = FASTX.quality_scores(record)  # Returns 0-60 range
```

### 2. **String Conversion Elimination**
```julia
# WRONG (found in multiple places)
sequence_string = string(biosequence)
kmer_string = string(kmer)

# RIGHT (maintain types)
sequence::BioSequences.LongDNA{4} = biosequence
kmer::Kmers.DNAKmer{3} = kmer
```

### 3. **SingleStrand Construction Pattern**
```julia
# WRONG (canonicalizing during construction)
canonical_kmer = BioSequences.canonical(observed_kmer)

# RIGHT (record exactly what we observe)
observed_kmers = [kmer for kmer in Kmers.EveryKmer(sequence)]
# Canonicalize LATER if DoubleStrand mode chosen
```

### 4. **Evidence Tracking Consistency**
```julia
# Standardized evidence tracking across all graph types
# Double-nested dictionary: dataset_id -> observation_id -> Set{EvidenceEntry}

struct EvidenceEntry
    position::Int              # Position within observation
    strand::StrandOrientation  # Forward or Reverse
end

# For quality-aware graphs
struct QualityEvidenceEntry
    position::Int
    strand::StrandOrientation
    quality_scores::Vector{UInt8}  # Phred scores 0-60
end

# Evidence structure in vertex/edge data
evidence::Dict{String, Dict{String, Set{EvidenceEntry}}}
# Key1: dataset_id (enables multi-dataset analysis)
# Key2: observation_id (enables read-centric algorithms)
# Value: Set of evidence entries for that observation
```

---

## 📊 Success Metrics

1. **All 24 test files pass** with new rhizomorph implementation
2. **No string conversions** in k-mer/qualmer/BioSequence pipelines
3. **Proper Phred score ranges** (0-60) throughout quality-aware graphs
4. **Strand-specific graphs maintain exact observations** with proper strand tracking
5. **Non-strand-specific transformations** correctly merge evidence across RC pairs
6. **Double-nested evidence dictionaries** enable O(1) dataset and observation queries
7. **Interconversion functions** work bidirectionally between all four graph modes
8. **Type stability** verified across all path reconstruction operations
9. **Performance maintained or improved** vs current implementation
10. **Clean module organization** with no duplicate functionality
11. **Evidence helper functions** provide clean API abstracting nested structure

---

## Expected Outcomes

After implementing this plan:

- **Clean, maintainable architecture** following biological principles
- **Robust graph ecosystem** handling all 24 graph type combinations
- **Type-safe operations** throughout the assembly pipeline
- **Proper scientific accuracy** in quality score handling
- **Clear conceptual model**:
  - Strand-specific vs non-strand-specific (methodological choice)
  - Single-stranded vs double-stranded representation (structural choice)
  - Four orthogonal combinations with bidirectional interconversion
- **Efficient evidence storage**:
  - Double-nested dictionaries: dataset_id → observation_id → evidence
  - O(1) queries at dataset and observation levels
  - 2-3× memory savings vs flat structure
  - Optimal for read-centric correction algorithms
- **Extensible framework** for future graph algorithms
- **Comprehensive test suite** ensuring correctness across all modes

---

## Implementation Checklist

### Phase 1: Core Infrastructure
- [x] Create rhizomorph module structure
- [x] Move enums from graph-core.jl (StrandOrientation, GraphMode)
- [x] Consolidate vertex data structures with double-nested evidence dictionaries
- [x] Consolidate edge data structures with double-nested evidence dictionaries
- [x] Implement evidence manipulation helper functions (add, query, filter, merge, stats)
- [ ] Implement strand specificity interconversion functions (strand-specific ↔ non-strand-specific) **(missing public API)**
- [x] Implement strand representation interconversion functions (single-stranded ↔ double-stranded)
- [ ] Implement graph type interconversion functions:
  - [ ] Fixed-length → Variable-length (collapse unbranching paths) **(planned)**
  - [ ] Variable-length → Fixed-length (fragment sequences) **(planned)**
  - [ ] Quality-aware → Quality-unaware (remove quality scores) **(planned)**
- [x] Implement shared graph construction patterns

### Phase 2: Fixed-Length Graphs
- [x] Migrate k-mer graphs with SingleStrand-first approach (singlestrand/doublestrand/canonical)
- [x] Fix qualmer graphs with proper Phred handling (singlestrand/doublestrand/canonical)
- [x] Create separate n-gram graphs module (singlestrand only)
- [ ] Update all imports and exports **(needs compatibility shims + legacy removal plan)**

### Phase 3: Variable-Length Graphs
- [x] Consolidate FASTA graph functionality (singlestrand with DNA/RNA double/canonical conversions)
- [x] Fix FASTQ graphs with quality handling (singlestrand with DNA/RNA double/canonical conversions)
- [x] Separate variable-length string graphs (singlestrand only)
- [ ] Implement alignment-based edge detection **(planned)**

### Phase 4: Algorithms & Utilities
- [x] Move path finding algorithms
- [x] Implement canonicalization as post-processing
- [x] Consolidate I/O functionality
- [ ] Implement Read-Centric Probabilistic Correction Model (Primary Strategy):
  - [ ] Quality-aware mode: Joint probability from Phred scores
  - [ ] Coverage-only mode: Evidence depth-based likelihood
  - [ ] Probabilistic path selection: P(alt) / P(current)
  - [ ] Iterative refinement with graph rebuilding
  - [ ] K-fold cross-validation with random seeds
  - [ ] Coverage-based confidence adjustments for high observation counts
- [ ] Implement Inferred Error Rate Modeling:
  - [ ] Identify consensus/highway paths vs alternate/back-road paths
  - [ ] Classify alternatives: insertions, deletions, substitutions
  - [ ] Estimate rates from graph structure
  - [ ] Future: Pre-trained error models from gold-standard sequences
- [ ] Implement Graph-Based Simplification (Supplementary Tools):
  - [ ] Tip Removal (Clipping) with configurable thresholds **(missing)**
  - [ ] Bubble Popping and Variant Calling (FreeBayes/GATK-style) **(partial; detection exists, lacks tests)**
  - [ ] Cycle/Repeat Resolution **(planned)**
- [ ] Future: Reinforcement Learning framework (post-publication)

### Phase 5: Integration & Testing
- [ ] Update main module integration **(needs legacy graph includes removed after shims/tests)**
- [ ] Update all 24 test files **(legacy tests must migrate to rhizomorph APIs)**
- [ ] Run comprehensive validation suite
- [ ] Performance benchmarking
- [ ] Documentation updates

---

## End-to-End Workflow Example

**Complete FASTQ → Assembly Pipeline:**

```julia
import Mycelia

# ============================================================================
# 1. Build quality-aware graph from FASTQ
# ============================================================================

graph = Mycelia.build_graph_from_fastq_file("sample.fastq", k=31)
# Returns: quality-aware k-mer graph with full evidence tracking
# Default mode: strand-specific, single-strand

# ============================================================================
# 2. Save intermediate state (lossless)
# ============================================================================

Mycelia.save_graph(graph, "initial_graph.jld2")
# Preserves: all evidence, quality scores, dataset/observation IDs

# ============================================================================
# 3. Error rate inference from graph structure
# ============================================================================

error_rates = Mycelia.infer_error_rates_from_graph(graph)
# Returns: ErrorRateModel with insertion/deletion/substitution rates
# Estimated from consensus vs alternative path frequencies

# ============================================================================
# 4. Probabilistic read-centric error correction
# ============================================================================

corrected = Mycelia.probabilistic_read_correction!(
    graph,
    error_rates=error_rates,
    n_iterations=3,
    random_seed=42
)
# Primary correction strategy: evaluates each read against graph
# Accepts corrections probabilistically based on quality scores
# Preserves low-frequency variation stochastically

# ============================================================================
# 5. Save corrected graph
# ============================================================================

Mycelia.save_graph(corrected, "corrected_graph.jld2")

# ============================================================================
# 6. Optional: graph-based simplification (supplementary)
# ============================================================================

# Tip clipping (remove short dead-ends)
simplified = Mycelia.clip_tips!(corrected, min_length=31, min_evidence=5)

# Variant calling (identify and report bubbles)
variants = Mycelia.call_variants(simplified, mode=:variant_calling)

# ============================================================================
# 7. Assemble contigs from graph
# ============================================================================

paths = Mycelia.find_eulerian_paths(simplified)
contigs = [Mycelia.path_to_sequence(p, simplified) for p in paths]

# ============================================================================
# 8. Export for publication/visualization
# ============================================================================

# GFA for community tools (Bandage, etc.) - lossy but standard
Mycelia.export_to_gfa(simplified, "assembly.gfa", include_evidence_summary=true)

# JLD2 for archive and future analysis - lossless
Mycelia.save_graph(simplified, "final_assembly.jld2")

# ============================================================================
# 9. Optional: K-fold cross-validation
# ============================================================================

# Run assembly with different random seeds
seeds = [42, 123, 456, 789, 1011]
assemblies = []

for seed in seeds
    g = Mycelia.build_graph_from_fastq_file("sample.fastq", k=31)
    corrected = Mycelia.probabilistic_read_correction!(g, random_seed=seed)
    paths = Mycelia.find_eulerian_paths(corrected)
    push!(assemblies, [Mycelia.path_to_sequence(p, corrected) for p in paths])
end

# Compare: variants preserved across assemblies have higher confidence
consensus_variants = Mycelia.find_consensus_variants(assemblies)
```

**Key Design Principles Illustrated:**

1. **Lossless intermediate storage** (JLD2) allows resuming work without information loss
2. **Probabilistic correction** preserves low-frequency variation stochastically
3. **Error rate inference** from graph structure (no gold standard needed)
4. **Quality-aware** throughout (Phred scores integrated into all likelihood calculations)
5. **Read-centric** primary strategy, graph-based operations supplementary
6. **K-fold validation** for variant confidence
7. **Dual I/O** strategy: JLD2 for preservation, GFA for interoperability

---

*This architecture follows the core insight: SingleStrand graphs are simple and direct, DoubleStrand graphs should canonicalize after construction to maintain traversal validity. The rhizomorph ecosystem will provide a solid foundation for all downstream assembly operations.*
