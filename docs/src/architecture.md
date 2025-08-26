# Architecture Overview

This document provides a comprehensive overview of Mycelia's architecture, design decisions, and system organization.

## Table of Contents

1. [System Overview](#system-overview)
2. [Core Components](#core-components)
3. [Module Organization](#module-organization)
4. [Data Flow](#data-flow)
5. [Type System](#type-system)
6. [Design Patterns](#design-patterns)
7. [Key Algorithms](#key-algorithms)
8. [Performance Considerations](#performance-considerations)

## System Overview

Mycelia is designed as a modular bioinformatics toolkit with a focus on probabilistic genome assembly. The architecture emphasizes:

- **Modularity**: Independent components for different bioinformatics tasks
- **Extensibility**: Easy to add new algorithms and methods
- **Performance**: Efficient memory usage and parallel processing
- **Interoperability**: Integration with existing bioinformatics tools

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    User Interface                        │
│  (Julia API, Command Line, Jupyter Notebooks)          │
└─────────────────────────────────────────────────────────┘
                             │
┌─────────────────────────────────────────────────────────┐
│                    Core Modules                          │
├─────────────────────┬───────────────┬──────────────────┤
│   Sequence I/O      │   Assembly    │   Analysis       │
│   - FASTX          │   - Qualmer   │   - K-mer        │
│   - Converters     │   - String    │   - Alignment    │
│                    │   - Intelligent│   - QC           │
└─────────────────────┴───────────────┴──────────────────┘
                             │
┌─────────────────────────────────────────────────────────┐
│                 Foundation Layer                         │
├──────────────┬──────────────────┬──────────────────────┤
│  BioSequences│    Graphs        │     Utilities        │
│  - DNA/RNA   │  - MetaGraphsNext│  - Memory Mgmt      │
│  - Proteins  │  - Algorithms    │  - Parallel         │
└──────────────┴──────────────────┴──────────────────────┘
```

## Core Components

### 1. Sequence Management (`fastx.jl`)

Handles all sequence I/O operations:
- **FASTA/FASTQ parsing** via FASTX.jl
- **Format conversion** between sequence types
- **Quality score handling** for FASTQ data
- **Streaming support** for large files

### 2. Assembly Engines

#### Qualmer Assembly (`qualmer-analysis.jl`)
- **Quality-aware k-mers**: Combines sequence and quality information
- **Probabilistic scoring**: Joint probability calculations
- **Graph construction**: Quality-weighted De Bruijn graphs

#### String Graphs (`string-graphs.jl`)
- **N-gram analysis**: Flexible substring decomposition
- **Graph simplification**: Path collapse algorithms
- **Component assembly**: Connected component traversal

#### Intelligent Assembly (`intelligent-assembly.jl`)
- **Dynamic k-mer selection**: Prime number progression
- **Error correction**: Probabilistic path finding
- **Memory management**: Adaptive processing
- **Reward-based optimization**: Accuracy-prioritized decisions

### 3. Analysis Tools

- **K-mer Analysis** (`kmer-analysis.jl`): Counting, spectrum analysis
- **Alignment** (`alignments-and-mapping.jl`): Sequence comparison
- **Quality Control** (`quality-control-and-benchmarking.jl`): Metrics and validation
- **Clustering** (`clustering.jl`): Sequence grouping algorithms

### 4. Utility Functions (`utility-functions.jl`)

- **Memory estimation**: Predict resource requirements
- **File handling**: Robust I/O operations
- **Progress tracking**: User feedback
- **Error handling**: Graceful failure modes

## Module Organization

### Source File Structure

```
src/
├── Mycelia.jl                          # Main module, imports all dependencies
├── fastx.jl                           # Sequence I/O
├── assembly.jl                        # Assembly coordination
├── intelligent-assembly.jl            # Smart assembly algorithm
├── qualmer-analysis.jl               # Quality-aware k-mers
├── string-graphs.jl                  # String graph methods
├── kmer-analysis.jl                  # K-mer operations
├── alignments-and-mapping.jl         # Alignment algorithms
├── annotation.jl                     # Gene annotation
├── clustering.jl                     # Clustering methods
├── classification.jl                 # Sequence classification
├── quality-control-and-benchmarking.jl # QC and validation
├── plotting-and-visualization.jl      # Data visualization
├── reference-databases.jl            # Database integration
├── taxonomy-and-trees.jl             # Phylogenetic analysis
├── variant-analysis.jl               # Variant calling
└── utility-functions.jl              # Helper functions
```

### Dependency Management

All package dependencies are imported at the top level in `Mycelia.jl`:

```julia
import BioSequences
import FASTX
import Graphs
import MetaGraphsNext
import DataFrames
# ... etc
```

Individual source files **do not** import packages - they access them through the module namespace.

## Data Flow

### Assembly Pipeline Flow

```
Input FASTQ
    │
    ▼
Read Loading & Validation
    │
    ▼
Quality Score Extraction
    │
    ├──► Qualmer Generation
    │         │
    │         ▼
    │    Quality-Aware Graph
    │         │
    │         ▼
    │    Probabilistic Paths
    │
    ├──► K-mer Analysis
    │         │
    │         ▼
    │    Coverage Estimation
    │
    ▼
Intelligent Assembly
    │
    ├──► Dynamic K Selection
    ├──► Error Correction
    ├──► Memory Monitoring
    └──► Reward Optimization
         │
         ▼
    Final Assembly
```

## Type System

### Core Types

#### Qualmer Type
```julia
struct Qualmer{KmerT, K}
    kmer::KmerT
    qualities::NTuple{K,UInt8}
end
```

Combines sequence (k-mer) with quality scores for probabilistic analysis.

#### Graph Types
- **MetaGraphs**: Vertex and edge metadata storage
- **QualmerVertexData**: Coverage, probability, observations
- **QualmerEdgeData**: Transition counts, quality weights

#### Sequence Types
Always use BioSequences types:
- `BioSequences.LongDNA{4}`
- `BioSequences.LongRNA{4}`
- `BioSequences.LongAA`

### Type Hierarchy

```
AbstractSequence
    ├── BioSequence
    │   ├── LongDNA
    │   ├── LongRNA
    │   └── LongAA
    │
    └── Kmer
        ├── DNAKmer
        ├── RNAKmer
        └── AAKmer
```

## Design Patterns

### 1. Functional Core, Imperative Shell
- **Pure functions** for algorithms
- **Mutable state** isolated to graph operations
- **Immutable data** for k-mers and sequences

### 2. Iterator Pattern
- **Lazy evaluation** for k-mer generation
- **Memory efficient** sequence processing
- **Composable** operations

### 3. Strategy Pattern
- **Multiple assembly methods** (qualmer, string, intelligent)
- **Pluggable algorithms** for each step
- **Runtime selection** based on data characteristics

### 4. Builder Pattern
- **Graph construction** with incremental updates
- **Assembly results** with metadata accumulation

## Key Algorithms

### 1. Prime K-mer Progression
```julia
function dynamic_k_prime_pattern(start_prime, max_k, initial_step)
    # Progressive stepping with prime discovery
    # Exploits mathematical properties of prime distribution
end
```

### 2. Quality-Weighted Path Finding
```julia
function find_quality_weighted_path(graph, start_vertex)
    # Greedy selection based on joint probability
    # Prioritizes high-confidence k-mers
end
```

### 3. Error Correction
```julia
function attempt_error_correction(graph, kmer, vertex_data)
    # Identifies low-coverage, low-quality k-mers
    # Suggests corrections based on graph context
end
```

### 4. Reward-Based Optimization
```julia
function calculate_assembly_reward(graph, corrections, k)
    # Multi-metric scoring: coverage, quality, connectivity
    # Guides iterative assembly decisions
end
```

## Performance Considerations

### Memory Management
- **Lazy loading**: Process sequences as streams
- **Graph sparsity**: Only store observed k-mers
- **Memory limits**: Configurable constraints
- **Garbage collection**: Explicit cleanup in loops

### Parallelization
- **Thread safety**: Immutable k-mer operations
- **Parallel k-mer counting**: Partition by sequence
- **Concurrent graph updates**: Lock-free where possible
- **BLAS operations**: Leverage optimized libraries

### Optimization Strategies
1. **K-mer size selection**: Start with optimal k based on error rate
2. **Early termination**: Stop when improvements plateau
3. **Sparse data structures**: Efficient for genomic data
4. **Type stability**: Concrete types for performance

### Scalability
- **Streaming algorithms**: Handle files larger than RAM
- **Hierarchical assembly**: Process in chunks
- **Distributed computing**: Future HPC support
- **Incremental updates**: Resume interrupted assemblies

## Extension Points

### Adding New Assembly Methods
1. Create new file in `src/`
2. Implement assembly function with standard interface
3. Register in main assembly dispatcher
4. Add tests and documentation

### Integrating External Tools
1. Use `run_external_command` wrapper
2. Handle tool installation via bioconda
3. Parse output formats appropriately
4. Add to tool registry

### Custom Quality Metrics
1. Extend `QualmerVertexData` with new fields
2. Implement scoring function
3. Integrate into reward calculation
4. Update graph construction

## Future Architecture Plans

### Planned Enhancements
- **GPU acceleration** for k-mer operations
- **Distributed assembly** across clusters
- **Cloud storage** integration
- **Real-time visualization** of assembly progress
- **Plugin system** for custom algorithms

### API Stability
- **Core API**: Stable for v1.0
- **Internal functions**: Subject to change
- **Deprecation policy**: One version warning
- **Semantic versioning**: Breaking changes in major versions