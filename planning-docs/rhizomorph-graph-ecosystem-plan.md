# Comprehensive Rhizomorph Graph Ecosystem Action Plan

*A detailed plan for implementing a robust "rhizomorph" graph ecosystem that properly handles all graph types with SingleStrand simplicity first and DoubleStrand canonicalization as a post-processing step.*

---

## 🎯 Core Architecture Principles

### 1. **SingleStrand-First Design Philosophy**
- **Always record exactly what we observe** in directed graphs
- Track coverage with `(observation_id, position, strand_orientation)`
- Build edges from k-1/n-1 overlaps or exact alignments
- **No canonicalization during construction** - preserve biological orientation

### 2. **DoubleStrand as Post-Processing**
- Build two parallel directed graphs (forward + reverse)
- **After construction**: identify reverse-complement pairs
- Optionally collapse to canonical representation
- Maintain traversal validity through explicit strand tracking

### 3. **Type System Consistency**
- **Never convert to strings** unless input type is String
- Maintain BioSequence types throughout k-mer/qualmer pipelines
- Use proper Phred scores (0-60) not ASCII characters
- Generic parameterized structures: `VertexData{T}`, `EdgeData{T}`

---

## 📁 Proposed Module Organization

```
src/rhizomorph/
├── core/
│   ├── enums.jl              # StrandOrientation, GraphMode
│   ├── vertex-data.jl        # Generic VertexData{T} structures
│   ├── edge-data.jl          # Generic EdgeData{T} structures
│   └── graph-construction.jl # Shared graph building logic
├── fixed-length/
│   ├── kmer-graphs.jl        # DNAKmer, RNAKmer, AAKmer graphs
│   ├── qualmer-graphs.jl     # Quality-aware k-mer graphs
│   └── ngram-graphs.jl       # String n-gram graphs
├── variable-length/
│   ├── fasta-graphs.jl       # Variable-length BioSequence graphs
│   ├── fastq-graphs.jl       # Variable-length quality BioSequence graphs
│   └── string-graphs.jl      # Variable-length string graphs
├── algorithms/
│   ├── path-finding.jl       # Eulerian paths, reconstruction
│   ├── canonicalization.jl   # DoubleStrand post-processing
│   ├── simplification.jl     # Bubble removal, graph cleanup
│   └── io.jl                # GFA export/import
└── rhizomorph.jl            # Main module file
```

---

## 🔧 Detailed Migration Plan

### Phase 1: Core Infrastructure (Week 1)

#### 1.1 Create Rhizomorph Module Structure
```bash
mkdir -p src/rhizomorph/{core,fixed-length,variable-length,algorithms}
```

#### 1.2 Move Shared Enums and Types (`src/rhizomorph/core/`)

**`enums.jl`** - Move from `graph-core.jl`:
```julia
@enum StrandOrientation Forward=true Reverse=false
@enum GraphMode SingleStrand DoubleStrand
```

**`vertex-data.jl`** - Consolidate all vertex data structures:
```julia
# Generic coverage entry
struct CoverageEntry
    observation_id::String
    position::Int
    strand::StrandOrientation
end

struct QualityCoverageEntry
    observation_id::String
    position::Int
    strand::StrandOrientation
    quality_scores::Vector{UInt8}  # Phred scores 0-60
end

# Fixed-length vertex data
struct KmerVertexData{T}
    canonical_element::T
    coverage::Dict{String, Set{CoverageEntry}}
end

struct QualmerVertexData{T}
    canonical_element::T
    coverage::Dict{String, Set{QualityCoverageEntry}}
    mean_quality::Float64
end

# Variable-length vertex data
struct BioSequenceVertexData{T}
    sequence::T
    coverage::Dict{String, Set{CoverageEntry}}
end

struct QualityBioSequenceVertexData{T}
    sequence::T
    coverage::Dict{String, Set{QualityCoverageEntry}}
    mean_quality::Float64
end

struct StringVertexData
    string_value::String
    coverage::Dict{String, Set{CoverageEntry}}
end

struct QualityStringVertexData
    string_value::String
    coverage::Dict{String, Set{QualityCoverageEntry}}
    mean_quality::Float64
end
```

---
**Design Note on Provenance Storage:**

To ensure both uniqueness of evidence and efficient querying, the `coverage` field will be implemented as a `Dict{Int, Set{CoverageEntry}}` where the key is the `dataset_id` and the value is a `Set` of unique `CoverageEntry` tuples for that dataset. This structure is chosen over a simple `Vector` or `Set` to allow for rapid, dataset-specific queries, which is critical for comparative analysis and resolving complex genomic regions. While this introduces a minor memory overhead compared to a flat list, the performance benefits for downstream algorithms are substantial.

---

**`edge-data.jl`** - Consolidate edge structures:
```julia
struct KmerEdgeData
    coverage::Set{Pair{CoverageEntry, CoverageEntry}}
    weight::Float64
    src_strand::StrandOrientation
    dst_strand::StrandOrientation
end

struct BioSequenceEdgeData
    overlap_length::Int
    coverage::Set{Pair{CoverageEntry, CoverageEntry}}
    weight::Float64
end

struct StringEdgeData
    overlap_length::Int
    coverage::Set{Pair{CoverageEntry, CoverageEntry}}
    weight::Float64
end

struct QualmerEdgeData
    coverage::Set{Pair{QualityCoverageEntry, QualityCoverageEntry}}
    weight::Float64
    src_strand::StrandOrientation
    dst_strand::StrandOrientation
end

struct QualityBioSequenceEdgeData
    overlap_length::Int
    coverage::Set{Pair{QualityCoverageEntry, QualityCoverageEntry}}
    weight::Float64
end

struct QualityStringEdgeData
    overlap_length::Int
    coverage::Set{Pair{QualityCoverageEntry, QualityCoverageEntry}}
    weight::Float64
end
```

#### 1.3 Extract Graph Construction Logic (`src/rhizomorph/core/graph-construction.jl`)

**SingleStrand Construction Pattern:**
```julia
function build_singlestrand_graph(
    elements::Vector{T},
    observations::Vector{<:FASTX.Record},
    vertex_data_type::Type{<:VertexData{T}},
    edge_data_type::Type
) where T
    # 1. Create graph with all observed elements
    graph = MetaGraphsNext.MetaGraph(...)

    # 2. Add vertices for each element
    for element in elements
        graph[element] = vertex_data_type(element)
    end

    # 3. Process each observation independently
    for (obs_idx, obs) in enumerate(observations)
        path = extract_path(obs, elements)
        add_observation_to_graph!(graph, path, obs_idx)
    end

    return graph
end
```

**DoubleStrand Construction Pattern:**
```julia
function build_doublestrand_graph(
    elements::Vector{T},
    observations::Vector{<:FASTX.Record},
    vertex_data_type::Type{<:VertexData{T}},
    edge_data_type::Type
) where T
    # 1. Build SingleStrand graph first
    forward_graph = build_singlestrand_graph(elements, observations, vertex_data_type, edge_data_type)

    # 2. Build reverse complement graph
    rc_elements = [reverse_complement(elem) for elem in elements if can_reverse_complement(elem)]
    reverse_graph = build_singlestrand_graph(rc_elements, observations, vertex_data_type, edge_data_type)

    # 3. Merge graphs maintaining strand information
    merged_graph = merge_strand_graphs(forward_graph, reverse_graph)

    return merged_graph
end
```

---

### Phase 2: Fixed-Length Graphs (Week 2)

#### 2.1 Migrate K-mer Graphs (`src/rhizomorph/fixed-length/kmer-graphs.jl`)

**From:** `src/kmer-graphs.jl` + `sequence-graphs-next.jl` k-mer sections
**Actions:**
1. Move `build_kmer_graph_next()` implementation
2. **Critical Fix**: Remove canonicalization during construction for SingleStrand
3. Implement proper DoubleStrand as two-graph merge
4. Update to use new vertex/edge data structures
5. Add validation for k-mer type consistency

**Key Changes:**
```julia
function build_kmer_graph_next(
    kmer_type::Type{<:Kmers.Kmer},
    observations::Vector{<:FASTX.Record};
    graph_mode::GraphMode = SingleStrand
)
    if graph_mode == SingleStrand
        return build_singlestrand_kmer_graph(kmer_type, observations)
    else
        return build_doublestrand_kmer_graph(kmer_type, observations)
    end
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
4. Fix coverage tracking (currently broken)

---

### Phase 3: Variable-Length Graphs (Week 3)

#### 3.1 Consolidate FASTA Graphs (`src/rhizomorph/variable-length/fasta-graphs.jl`)

**From:** `src/fasta-graphs.jl` + BioSequence sections in `sequence-graphs-next.jl`
**Actions:**
1. Move all variable-length BioSequence graph construction
2. Implement proper alignment-based edge detection using exact matching, alignemnts with BioAlignments, and external tools (e.g. minimap2)
3. Fix coverage tracking for variable-length vertices

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

#### 4.4 Graph Simplification Algorithms (`src/rhizomorph/algorithms/simplification.jl`)

**New Implementation:**
Implement a suite of robust, statistically-grounded graph simplification algorithms. These algorithms will support both graph-centric and read-centric correction models, with a strong emphasis on probabilistic methods inspired by established bioinformatics tools.

##### 4.4.1 Tip Removal (Clipping)
Prunes short, dead-end paths (tips) that typically result from sequencing errors or coverage drop-off.

*   **Functionality:** The tip clipping function will accept thresholds based on:
    *   `min_length::Int`: The minimum length of a tip to be retained.
    *   `min_coverage::Float64`: The minimum average coverage of a tip to be retained.
    *   `min_confidence::Float64`: For quality-aware graphs, the minimum average vertex confidence to be retained.
*   **Logic:** Users can specify the relationship between these thresholds using `and`, `or`, or `either` logic to determine which tips are pruned.
*   **Default Thresholds:**
    *   **Length:** For fixed-length graphs, the default will be the k-mer/n-gram size. For variable-length graphs, this will not be set by default.
    *   **Coverage/Confidence:** The default threshold will be calculated statistically from the graph properties. It will be defined as a set number of standard deviations below the mean coverage/confidence of all vertices. An alternative based on the median and coefficient of variation will also be explored for skewed distributions.

##### 4.4.2 Bubble Popping
Merges simple divergent paths (bubbles) that typically represent SNPs, indels, or sequencing errors.

*   **Approach:** Bubbles will be identified as simple alternative paths that start and end at the same vertices. The decision to "pop" a bubble (i.e., remove the minor path) will be based on a similar set of configurable thresholds as tip clipping (`min_length`, `min_coverage`, `min_confidence`).
*   **Context-Awareness:** The algorithm will be aware of the expected bubble size based on the graph type. A SNP in a fixed-length graph creates a bubble of length `k`, whereas in a variable-length graph, the bubble's length corresponds to the variant region.
*   **Probabilistic Popping:** Drawing inspiration from tools like **FreeBayes** and **GATK**, we will implement a Bayesian or maximum-likelihood model. This model will evaluate the evidence (coverage, quality scores, strand bias) for each path in the bubble and calculate the posterior probability that the minor path represents a true biological variant versus a sequencing artifact. Bubbles will be popped if the probability of the alternative path being an error exceeds a user-defined threshold.

##### 4.4.3 Cycle and Repeat Resolution
Addresses simple cycles and repeats in the graph, which are often the most challenging aspect of genome assembly.

*   **Strategy:** The primary strategy will be to leverage the read-centric correction model (described below). For direct graph-based resolution, simple, short tandem repeats will be identified and resolved based on probabilistic criteria. The algorithm will assess the likelihood of the repeat's copy number given the spanning read evidence and associated quality scores.

##### 4.4.4 Read-Centric Correction Model
As a powerful alternative to direct graph modification, we will implement a read-centric correction framework.

*   **Workflow:**
    1.  An initial graph is constructed from the raw reads.
    2.  Each read is re-aligned to the graph to find its most likely path.
    3.  During alignment, a probabilistic model evaluates whether mismatches between the read and its path are more likely to be sequencing errors in the read or evidence for a true variant missing from the graph path. This model will heavily rely on the read's quality scores.
    4.  Based on this evaluation, reads can be "corrected," or low-confidence regions of the graph can be flagged.
    5.  The graph is then rebuilt using the corrected reads or down-weighted evidence.
*   **Benefit:** This iterative process of graph construction and read re-alignment can resolve complex errors, tips, and bubbles in a more robust and data-driven manner than graph-centric heuristics alone.

#### 4.5 I/O Functions (`src/rhizomorph/algorithms/io.jl`)

**From:** GFA functions scattered across modules
**Actions:**
1. Consolidate all import/export functionality
2. Implement strand-aware GFA output
3. Add metadata preservation for different graph types

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
4. Validate coverage tracking and provenance

#### 5.3 Comprehensive Validation

**Run test suite progression:**
1. Basic construction tests for each graph type
2. Type stability validation across all operations
3. SingleStrand vs DoubleStrand correctness
4. Coverage and provenance tracking
5. Path reconstruction and sequence identity
6. Performance benchmarks vs current implementation

---

## Critical Issues to Address

### 1. **Immediate Phred Score Fix**
```julia
# WRONG (current)
quality_chars = FASTX.quality(record)
quality_ints = [Int(c) for c in quality_chars]

# RIGHT (target)
quality_phred = get_phred_scores(record)  # Returns 0-60 range
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
observed_kmers = [kmer for kmer in extract_kmers(sequence)]
# Canonicalize LATER if DoubleStrand mode chosen
```

### 4. **Coverage Tracking Consistency**
```julia
# Standardize across all graph types
struct CoverageEntry
    observation_id::String      # Which input sequence
    position::Int           # Position within that sequence
    strand::StrandOrientation  # Forward or Reverse
end
```

---

## 📊 Success Metrics

1. **All 24 test files pass** with new rhizomorph implementation
2. **No string conversions** in k-mer/qualmer/BioSequence pipelines
3. **Proper Phred score ranges** (0-60) throughout quality-aware graphs
4. **SingleStrand graphs maintain exact observations**
5. **DoubleStrand graphs show proper canonicalization**
6. **Type stability** verified across all path reconstruction operations
7. **Performance maintained or improved** vs current implementation
8. **Clean module organization** with no duplicate functionality

---

## Expected Outcomes

After implementing this plan:

- **Clean, maintainable architecture** following biological principles
- **Robust graph ecosystem** handling all 24 graph type combinations
- **Type-safe operations** throughout the assembly pipeline
- **Proper scientific accuracy** in quality score handling
- **Clear separation** between observation (SingleStrand) and analysis (DoubleStrand)
- **Extensible framework** for future graph algorithms
- **Comprehensive test coverage** ensuring correctness

---

## Implementation Checklist

### Phase 1: Core Infrastructure
- [ ] Create rhizomorph module structure
- [ ] Move enums from graph-core.jl
- [ ] Consolidate vertex data structures
- [ ] Consolidate edge data structures
- [ ] Implement shared graph construction patterns

### Phase 2: Fixed-Length Graphs
- [ ] Migrate k-mer graphs with SingleStrand-first approach
- [ ] Fix qualmer graphs with proper Phred handling
- [ ] Create separate n-gram graphs module
- [ ] Update all imports and exports

### Phase 3: Variable-Length Graphs
- [ ] Consolidate FASTA graph functionality
- [ ] Fix FASTQ graphs with quality handling
- [ ] Separate variable-length string graphs
- [ ] Implement alignment-based edge detection

### Phase 4: Algorithms & Utilities
- [ ] Move path finding algorithms
- [ ] Implement canonicalization as post-processing
- [ ] Consolidate I/O functionality
- [ ] Implement Tip Removal (Clipping)
- [ ] Implement Bubble Popping with probabilistic models
- [ ] Implement Cycle/Repeat Resolution
- [ ] Implement Read-Centric Correction Model

### Phase 5: Integration & Testing
- [ ] Update main module integration
- [ ] Update all 24 test files
- [ ] Run comprehensive validation suite
- [ ] Performance benchmarking
- [ ] Documentation updates

---

*This architecture follows the core insight: SingleStrand graphs are simple and direct, DoubleStrand graphs should canonicalize after construction to maintain traversal validity. The rhizomorph ecosystem will provide a solid foundation for all downstream assembly operations.*