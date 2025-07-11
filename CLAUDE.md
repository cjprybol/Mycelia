# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Mycelia is a Julia package for bioinformatics and computational biology, providing tools for sequence analysis, genome assembly, annotation, and comparative genomics. The package is designed as a comprehensive bioinformatics toolkit with modular functionality.

## Development Commands

### Testing
```bash
# Run all tests from package root
julia --project=. -e "import Pkg; Pkg.test()"

# Run with coverage
julia --project=. --code-coverage=user -e "import Pkg; Pkg.test()"

# Run with color output
julia --project=. --color=yes -e "import Pkg; Pkg.test()"
```

### Static Analysis
```bash
# Run JET static analysis
julia --project=. test/jet.jl
```

### Documentation
```bash
# Build documentation (from docs/ directory)
julia --project=docs make.jl
```

## Architecture

### Core Structure
- **src/Mycelia.jl**: Main module file that dynamically imports all other source files
- **src/**: Modular functionality organized by domain:
  - `fastx.jl`: FASTA/FASTQ sequence handling
  - `assembly.jl`: Genome assembly tools
  - `annotation.jl`: Gene annotation functionality
  - `alignments-and-mapping.jl`: Sequence alignment tools
  - `kmer-analysis.jl`: K-mer counting and analysis
  - `clustering.jl`: Clustering algorithms
  - `classification.jl`: Sequence classification
  - `quality-control-and-benchmarking.jl`: QC and validation
  - `plotting-and-visualization.jl`: Data visualization
  - `utility-functions.jl`: General utility functions
  - `reference-databases.jl`: Database integration
  - `taxonomy-and-trees.jl`: Phylogenetic analysis
  - `variant-analysis.jl`: Variant calling and analysis

### Test Organization
Tests are organized in numbered directories following the bioinformatics workflow:
1. `1_data_acquisition/`: Data download and simulation
2. `2_preprocessing_qc/`: Quality control and preprocessing
3. `3_feature_extraction_kmer/`: K-mer analysis
4. `4_assembly/`: Genome assembly
5. `5_validation/`: Assembly validation
6. `6_annotation/`: Gene annotation
7. `7_comparative_pangenomics/`: Comparative genomics
8. `8_tool_integration/`: External tool integration

### Key Dependencies
- **BioSequences.jl**: Core sequence data structures
- **FASTX.jl**: FASTA/FASTQ parsing
- **BioAlignments.jl**: Sequence alignment
- **Graphs.jl**: Graph algorithms for assembly
- **DataFrames.jl**: Data manipulation
- **Makie.jl/Plots.jl**: Visualization
- **Arrow.jl/JLD2.jl**: Data serialization
- **XAM.jl**: BAM/SAM file handling

## Development Notes

### Module Loading
The main module uses dynamic file inclusion - all `.jl` files in `src/` are automatically included. When adding new functionality, create appropriately named files in `src/` and they will be automatically loaded.

### Memory Management
The package includes utilities for memory estimation and checking (see `utility-functions.jl`). Large-scale genomic analyses should use these tools to avoid memory issues.

### Testing Approach
- Tests follow the bioinformatics workflow order
- Use `include_all_tests()` function to recursively load test files
- Test data is stored in `test/metadata/`

### Documentation
- Functions use Julia docstrings with DocStringExtensions
- Documentation is built using Documenter.jl
- Documentation source is in `docs/src/`

## External Tool Integration

The package integrates with various bioinformatics tools:
- Bioconda for package management
- SLURM for job scheduling
- Rclone for cloud storage
- Neo4j for graph databases