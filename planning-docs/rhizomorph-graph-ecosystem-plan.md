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
    observation_id::Int
    position::Int
    strand::StrandOrientation
end

struct QualityCoverageEntry
    observation_id::Int
    position::Int
    strand::StrandOrientation
    quality_scores::Vector{Int}  # Phred scores 0-60
end

# Fixed-length vertex data
struct KmerVertexData{T}
    canonical_element::T
    coverage::Vector{CoverageEntry}
end

struct QualmerVertexData{T}
    canonical_element::T
    coverage::Vector{QualityCoverageEntry}
    mean_quality::Float64
end

# Variable-length vertex data
struct BioSequenceVertexData{T}
    sequence::T
    coverage::Vector{CoverageEntry}
end

struct QualityBioSequenceVertexData{T}
    sequence::T
    coverage::Vector{QualityCoverageEntry}
    mean_quality::Float64
end

struct StringVertexData
    string_value::String
    coverage::Vector{CoverageEntry}
end

struct QualityStringVertexData
    string_value::String
    coverage::Vector{QualityCoverageEntry}
    mean_quality::Float64
end
```

**`edge-data.jl`** - Consolidate edge structures:
```julia
struct KmerEdgeData
    coverage::Vector{Tuple{CoverageEntry, CoverageEntry}}
    weight::Float64
    src_strand::StrandOrientation
    dst_strand::StrandOrientation
end

struct BioSequenceEdgeData
    overlap_length::Int
    coverage::Vector{Tuple{CoverageEntry, CoverageEntry}}
    weight::Float64
end

struct StringEdgeData
    overlap_length::Int
    coverage::Vector{Tuple{CoverageEntry, CoverageEntry}}
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
2. Implement proper alignment-based edge detection
3. Fix coverage tracking for variable-length vertices

#### 3.2 Consolidate FASTQ Graphs (`src/rhizomorph/variable-length/fastq-graphs.jl`)

**From:** `src/fastq-graphs.jl` + quality BioSequence sections
**Actions:**
1. **Critical**: Fix Phred score handling throughout
2. Implement quality aggregation for overlapping sequences
3. Add provenance tracking for multiple observations

#### 3.3 Consolidate String Graphs (`src/rhizomorph/variable-length/string-graphs.jl`)

**From:** `src/string-graphs.jl` variable-length functions
**Actions:**
1. Move variable-length string graph functionality
2. Implement text-based alignment detection
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

#### 4.3 I/O Functions (`src/rhizomorph/algorithms/io.jl`)

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

# Export main graph construction functions
export build_kmer_graph_next, build_qualmer_graph, build_fasta_graph, build_fastq_graph
export build_string_graph, build_ngram_graph
export StrandOrientation, GraphMode, Forward, Reverse, SingleStrand, DoubleStrand
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

## 🚨 Critical Issues to Address

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
    observation_id::Int      # Which input sequence
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

## 🎉 Expected Outcomes

After implementing this plan:

- **Clean, maintainable architecture** following biological principles
- **Robust graph ecosystem** handling all 24 graph type combinations
- **Type-safe operations** throughout the assembly pipeline
- **Proper scientific accuracy** in quality score handling
- **Clear separation** between observation (SingleStrand) and analysis (DoubleStrand)
- **Extensible framework** for future graph algorithms
- **Comprehensive test coverage** ensuring correctness

---

## 📋 Implementation Checklist

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
- [ ] Add graph simplification algorithms

### Phase 5: Integration & Testing
- [ ] Update main module integration
- [ ] Update all 24 test files
- [ ] Run comprehensive validation suite
- [ ] Performance benchmarking
- [ ] Documentation updates

---

*This architecture follows the core insight: SingleStrand graphs are simple and direct, DoubleStrand graphs should canonicalize after construction to maintain traversal validity. The rhizomorph ecosystem will provide a solid foundation for all downstream assembly operations.*