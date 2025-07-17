# Strand-Aware K-mer Graph Implementation

## Overview

The next-generation sequence graph implementation introduces a crucial distinction between **vertex representation** and **edge directionality** that was missing in the original design. This addresses the biological reality of strand-specific transitions in genome assembly.

## Key Design Principles

### 1. Canonical Vertices, Strand-Aware Edges

**Previous approach**: Stored both k-mer and reverse complement as separate vertices
- Memory inefficient (2x vertices)
- Complex graph topology
- Difficult to identify canonical paths

**New approach**: Store only canonical k-mers as vertices, track strand in edges
- Memory efficient (canonical representation)
- Cleaner graph structure  
- Strand information preserved in edge metadata

### 2. Biologically Valid Transitions

The graph now enforces that edges represent **biologically valid transitions** where:
- The suffix of the source k-mer overlaps with the prefix of the destination k-mer
- Strand orientations are compatible for the transition
- Both single-strand and double-strand modes are supported

### 3. Two Graph Modes

#### SingleStrand Mode
- For RNA, amino acids, or directional DNA analysis
- All k-mers use forward orientation
- No reverse complement consideration
- Suitable for: transcriptome assembly, RNA-seq, protein/proteome assembly, directional cDNA

#### DoubleStrand Mode (Default)  
- For standard DNA assembly
- K-mers converted to canonical representation
- Strand information tracked in coverage and edges
- Reverse complement transitions allowed
- Suitable for: genome assembly, DNA-seq

## Data Structures

### StrandOrientation Enum
```julia
@enum StrandOrientation Forward=true Reverse=false
```

### KmerVertexData
```julia
struct KmerVertexData
    coverage::Vector{Tuple{Int, Int, StrandOrientation}}  # (obs_id, position, strand)
    canonical_kmer::String  # Always canonical representation
end
```

### KmerEdgeData  
```julia
struct KmerEdgeData
    coverage::Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}}
    weight::Float64
    src_strand::StrandOrientation  # Required source orientation
    dst_strand::StrandOrientation  # Required destination orientation
end
```

## Biological Correctness

### Transition Validation
Each edge represents a transition that is valid only when:
1. **Overlap constraint**: suffix(src_kmer) == prefix(dst_kmer)
2. **Strand compatibility**: orientations allow for proper base pairing
3. **Coverage support**: observed in the sequence data

### Example: Valid Transitions
For k=3, canonical k-mers "ATC" and "TCG":

**Valid transitions**:
- ATC(Forward) → TCG(Forward): 
  ```
  ATC  (suffix: TC)
   TCG (prefix: TC)  ✓ overlap matches
  = ATCG
  ```
- GAT(Reverse) → CGA(Reverse): 
  ```
  ATC  (RC of GAT, suffix: TC)
   TCG (RC of CGA, prefix: TC)  ✓ overlap matches  
  = ATCG
  ```

**Invalid transitions**:  
- ATC(Forward) → GAT(Forward): 
  ```
  ATC  (suffix: TC)
  GAT  (prefix: GA)  ✗ TC ≠ GA, no overlap
  ```
- ATC(Forward) → TCG(Reverse): 
  ```
  ATC  (suffix: TC)
  CGA  (RC of TCG, prefix: CG)  ✗ TC ≠ CG, wrong strand
  ```

## Assembly Implications

### Path Finding
- **Shortest probability paths**: distance ∝ -log(edge_weight)
- **Maximum weight walks**: follow highest coverage edges
- **Strand-consistent paths**: automatically enforced by edge constraints

### Error Correction
- **Quality-aware transitions**: edges weighted by coverage and quality
- **Consensus calling**: multiple observations of same transition
- **Bubble detection**: alternative paths with same src/dst vertices

### Performance Benefits
- **Memory efficiency**: ~50% reduction vs. stranded vertex representation
- **Cleaner algorithms**: canonical vertices simplify path operations
- **Type stability**: all metadata is type-stable for performance

## Migration Path

### Legacy Compatibility
The `legacy_to_next_graph()` function converts old MetaGraphs format:
1. Converts stranded vertices to canonical representation
2. Merges coverage from forward/reverse k-mer pairs
3. Creates strand-aware edges with proper orientation metadata
4. Maintains all original coverage information

### Gradual Adoption
- New code uses `build_kmer_graph_next()`
- Legacy code continues working via compatibility layer
- Tests validate both implementations during transition
- Performance benchmarks guide optimization priorities

## Future Enhancements

### Phase 2 Algorithms
Now that we have proper strand-aware graphs, we can implement:
- **Probabilistic walks** with strand-consistent transitions
- **Viterbi polishing** using canonical representation
- **Bubble resolution** algorithms
- **Repeat detection** and resolution

### Advanced Features
- **Quality score integration** in edge weights
- **Coverage depth thresholds** for edge filtering  
- **Graph simplification** algorithms
- **Interactive visualization** of strand-aware paths

This foundation provides the correct biological semantics needed for sophisticated assembly algorithms while maintaining performance through type-stable, memory-efficient data structures.
