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

### Intelligent Assembly

```@docs
mycelia_assemble(::Vector{<:FASTX.FASTQ.Record})
mycelia_assemble(::Union{String, Vector{String}})
```

### Assembly Utilities

```@docs
dynamic_k_prime_pattern
error_optimized_k_sequence
calculate_assembly_reward
```

## Qualmer Analysis

### Core Qualmer Types

```@docs
Qualmer
QualmerVertexData
QualmerEdgeData
```

### Qualmer Generation

```@docs
qualmers_fw
qualmers_canonical
qualmers_unambiguous
qualmers_unambiguous_canonical
```

### Quality Analysis

```@docs
phred_to_probability
qualmer_correctness_probability
joint_qualmer_probability
position_wise_joint_probability
```

### Graph Construction

```@docs
build_qualmer_graph
get_qualmer_statistics
```

## String Graphs

### N-gram Analysis

```@docs
ngrams
string_to_ngram_graph
plot_ngram_graph
```

### Graph Simplification

```@docs
collapse_unbranching_paths
find_connected_components
assemble_strings
```

## Sequence I/O

### File Handling

```@docs
open_fastx
download_genome_by_accession
```

### Format Conversion

```@docs
convert_sequence
```

## K-mer Analysis

### K-mer Counting

```@docs
count_canonical_kmers
```

### Distance Metrics

```@docs
jaccard_distance
kmer_counts_to_js_divergence
```

## Simulation and Testing

### Read Simulation

```@docs
simulate_pacbio_reads
simulate_nanopore_reads
```

### Assembly Testing

```@docs
test_intelligent_assembly
```

## Quality Control

### FASTQ Analysis

```@docs
analyze_fastq_quality
calculate_gc_content
assess_duplication_rates
```

### External Tool Integration

```@docs
qc_filter_short_reads_fastp
qc_filter_long_reads_filtlong
trim_galore_paired
```

## Assembly Validation

### Quality Assessment

```@docs
assess_assembly_quality
validate_assembly
```

### External Validators

```@docs
run_quast
run_busco
run_mummer
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

### Pangenome Analysis

```@docs
analyze_pangenome_kmers
build_genome_distance_matrix
```

## Visualization

### Plotting Functions

```@docs
plot_kmer_frequency_spectra
visualize_genome_coverage
plot_embeddings
plot_taxa_abundances
```

## Utility Functions

### Memory Management

```@docs
estimate_memory_usage
check_memory_limits
```

### File Operations

Helper functions for robust file handling and data validation.

### Progress Tracking

Functions for user feedback during long-running operations.

## Graph Types and Enums

### Graph Modes

```@docs
GraphMode
StrandOrientation
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
4. **Validation**: `assess_assembly_quality`, `validate_assembly`
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
3. **Enable logging**: Set `verbose=true` in functions
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

### Workflow Integration

Functions are designed to work together in pipelines:

```julia
# Typical workflow
data = open_fastx("reads.fastq")
quality_report = analyze_fastq_quality(data)
assembly = mycelia_assemble("reads.fastq")
validation = assess_assembly_quality(assembly)
```

## Version Compatibility

This documentation corresponds to Mycelia version 0.1.0+. 

- **Breaking changes** will be noted in release notes
- **Deprecated functions** include migration guidance
- **Experimental features** are clearly marked

For the latest updates, see the [CHANGELOG](../CHANGELOG.md).