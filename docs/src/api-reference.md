# API Reference

```@meta
CurrentModule = Mycelia
```

This page provides comprehensive documentation for all public functions in Mycelia, automatically generated from docstrings.

## Table of Contents

```@contents
Pages = ["api-reference.md"]
Depth = 3
```

## Assembly Functions

### Main Assembly Interface

```@docs
Mycelia.Rhizomorph.assemble_genome
```

### Third-Party Assembly Tool Wrappers

See [Assembly Suite](api/workflows/assembly-suite.md) for comprehensive documentation of production-ready assembly wrappers including:

- Short-read assemblers: MEGAHIT, SPAdes, metaSPAdes, SKESA, Velvet
- Long-read assemblers: Flye, metaFlye, Canu, hifiasm, MetaMDBG
- Hybrid assemblers: Unicycler
- Polishing tools: Apollo, HyPo, Strainy

### Experimental Assembly Algorithms

Mycelia includes experimental assembly algorithms in `src/development/` that are not currently part of the public API:

- **Intelligent Assembly**: `mycelia_assemble` - Self-optimizing k-mer progression
- **Assembly Utilities**: `dynamic_k_prime_pattern`, `error_optimized_k_sequence`, `calculate_assembly_reward`

These are research implementations subject to change. For production use, see the third-party assembler wrappers above.

## Qualmer Analysis

### Core Qualmer Types

```@docs
Mycelia.Qualmer
Mycelia.QualmerVertexData
Mycelia.QualmerEdgeData
```

### Qualmer Generation

```@docs
Mycelia.qualmers_fw
Mycelia.qualmers_canonical
Mycelia.qualmers_unambiguous
Mycelia.qualmers_unambiguous_canonical
```

### Quality Analysis

```@docs
Mycelia.phred_to_probability
Mycelia.qualmer_correctness_probability
Mycelia.joint_qualmer_probability
Mycelia.position_wise_joint_probability
```

### Graph Construction

```@docs
Mycelia.build_qualmer_graph
Mycelia.get_qualmer_statistics
```

## String Graphs

### N-gram Analysis

```@docs
Mycelia.ngrams
Mycelia.string_to_ngram_graph
Mycelia.plot_ngram_graph
```

### Graph Simplification

```@docs
Mycelia.collapse_unbranching_paths
Mycelia.find_connected_components
```

**Note**: Additional graph simplification function `assemble_strings` is available but currently lacks documentation.

## Sequence I/O

### File Handling

```@docs
Mycelia.open_fastx
Mycelia.write_fastq
Mycelia.write_fasta
```

### Format Conversion

```@docs
Mycelia.convert_sequence
```

## K-mer Analysis

### K-mer Counting

For k-mer counting functions, see [Sequence Analysis Workflow](api/workflows/sequence-analysis.md).

### Distance Metrics

```@docs
Mycelia.jaccard_distance
Mycelia.kmer_counts_to_js_divergence
```

### Genome Size Estimation

```@docs
Mycelia.estimate_genome_size_from_kmers
```

### K-mer Result Storage

```@docs
Mycelia.save_kmer_results
Mycelia.load_kmer_results
```

## Simulation and Testing

### Sequence Generation

```@docs
Mycelia.random_fasta_record
```

### Read Simulation

See [Data Acquisition Workflow](api/workflows/data-acquisition.md) for read simulation functions.

### Assembly Testing (Experimental)

Experimental assembly testing function `test_intelligent_assembly` is available in `src/development/` but not currently part of the public API.

## Quality Control

### FASTQ Analysis

```@docs
Mycelia.calculate_gc_content
Mycelia.assess_duplication_rates
```

For FASTQ quality analysis functions, see [Quality Control Workflow](api/workflows/quality-control.md).

### External Tool Integration

For quality control and filtering functions, see [Quality Control Workflow](api/workflows/quality-control.md).

## Assembly Validation

### Quality Assessment

```@docs
Mycelia.assess_assembly_quality
Mycelia.Rhizomorph.validate_assembly
```

### External Validators

```@docs
Mycelia.run_quast
Mycelia.run_busco
Mycelia.run_mummer
```

## Coverage and Abundance

```@docs
Mycelia.run_coverm_contig
Mycelia.run_coverm_genome
```

## Annotation

### Gene Prediction

Functions for gene prediction and annotation (using external tools):

- Pyrodigal integration
- BLAST+ searches
- MMSeqs2 searches
- TransTerm terminator prediction
- tRNAscan-SE tRNA detection
- MLST typing

## Comparative Genomics

## Visualization

### Plotting Functions

```@docs
Mycelia.visualize_genome_coverage
Mycelia.plot_embeddings
Mycelia.plot_taxa_abundances
```

For k-mer frequency plotting, see [Sequence Analysis Workflow](api/workflows/sequence-analysis.md).

## Utility Functions

### Memory Management (Experimental)

Experimental memory management functions (`estimate_memory_usage`, `check_memory_limits`) are available in `src/development/` but not currently part of the public API.

### File Operations

Helper functions for robust file handling and data validation.

### Progress Tracking

Functions for user feedback during long-running operations.

## Graph Types and Enums

### Graph Modes

```@docs
Mycelia.GraphMode
Mycelia.StrandOrientation
```

Used to control how graphs are constructed and oriented.

## Data Structures

### Assembly Results

The result of assembly operations typically includes:

- `:final_assembly`: Vector of assembled sequences
- `:k_progression`: K-mer sizes used
- `:metadata`: Assembly statistics and parameters

### Quality Metrics

Standard quality metrics returned by analysis functions:

- Coverage statistics
- Quality score distributions  
- Joint probability calculations
- Assembly completeness measures

## Function Categories

### Core Assembly Pipeline

1. **Input Processing**: `open_fastx`, `convert_sequence`
2. **Quality Analysis**: `analyze_fastq_quality`, qualmer functions
3. **Assembly**: `mycelia_assemble`, graph construction
4. **Validation**: `assess_assembly_quality`, `Mycelia.Rhizomorph.validate_assembly`
5. **Output**: Assembly statistics and final sequences

### Supporting Workflows

- **Data Acquisition**: Download and simulation functions
- **Quality Control**: Filtering and trimming tools
- **Analysis**: K-mer analysis and distance calculations
- **Visualization**: Plotting and reporting functions

## Type System

Mycelia uses a strict type system based on BioSequences.jl:

- **DNA sequences**: `BioSequences.LongDNA{4}`
- **RNA sequences**: `BioSequences.LongRNA{4}`  
- **Protein sequences**: `BioSequences.LongAA`
- **K-mers**: `Kmers.DNAKmer{K}`, `Kmers.RNAKmer{K}`, `Kmers.AAKmer{K}`

### Quality-Aware Types

- **Qualmer**: Combines k-mer with quality scores
- **QualmerVertexData**: Graph vertex with quality metadata
- **QualmerEdgeData**: Graph edge with transition weights

## Performance Characteristics

### Scalability

Most functions scale as follows with input size:

- **Linear**: File I/O, k-mer counting
- **Quadratic**: Distance matrix calculations
- **Variable**: Assembly (depends on complexity)

### Memory Usage

Typical memory requirements:

- **Small genomes** (< 10 Mb): 1-4 GB RAM
- **Bacterial genomes** (1-10 Mb): 4-16 GB RAM  
- **Large genomes** (> 100 Mb): 32-128 GB RAM

Use `memory_limit` parameters to control usage.

### Parallelization

Functions that support parallel execution:

- K-mer counting (automatic)
- Assembly (multi-threaded)
- Quality analysis (vectorized)

Set `JULIA_NUM_THREADS` environment variable for best performance.

## Error Handling

### Common Exceptions

- `ArgumentError`: Invalid parameters
- `MethodError`: Wrong types
- `OutOfMemoryError`: Insufficient RAM
- `SystemError`: File I/O issues

### Debugging Tips

1. **Check inputs**: Use `typeof()` and `length()` 
2. **Monitor memory**: Use `memory_limit` parameters
3. **Enable logging**: Set `verbose=true` in functions where applicable
4. **Test subset**: Use smaller datasets for debugging

## Integration Points

### External Tools

Mycelia integrates with numerous bioinformatics tools via Bioconda:

- **Assemblers**: SPAdes, Megahit, Miniasm, Raven
- **QC Tools**: FastP, Filtlong, Trim Galore
- **Validators**: QUAST, BUSCO, MUMmer
- **Annotators**: Pyrodigal, BLAST+, MMSeqs2

### File Formats

Supported formats:

- **Input**: FASTA, FASTQ (including gzipped)
- **Output**: FASTA, CSV, JSON, JLD2
- **Intermediate**: Graph formats, quality scores

## Version Compatibility

This documentation corresponds to Mycelia version 0.1.0+. 

- **Breaking changes** will be noted in release notes
- **Deprecated functions** include migration guidance
- **Experimental features** are clearly marked

For the latest updates, see the [CHANGELOG](https://github.com/cjprybol/Mycelia/blob/main/CHANGELOG.md).
