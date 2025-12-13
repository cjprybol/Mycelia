# Mycelia

! Warning: In Development !

**A Julia package for bioinformatics and computational biology**

[![codecov](https://codecov.io/github/cjprybol/Mycelia/graph/badge.svg?token=0ZQSER2FLR)](https://codecov.io/github/cjprybol/Mycelia)
[![Documentation](https://github.com/cjprybol/Mycelia/actions/workflows/documentation.yml/badge.svg)](https://cjprybol.github.io/Mycelia/dev/)

![Banner Logo](banner-logo.jpg)

## Overview

Mycelia is a Julia package for bioinformatics and computational biology that implements graph-based biosequence assembly and quality-aware sequence analysis. The package provides both research-oriented algorithms and practical bioinformatics functionality, including extensive tool integration for genomics workflows. While some components are experimental and in active development, the package includes substantial implemented functionality for data processing, assembly, annotation, and analysis.

## Quick Start

```julia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
import Mycelia
```

## Key Features

### Core Functionality

- **File Format Support**: FASTA/FASTQ/GenBank/GFF/VCF/SAM/BAM with automatic compression handling
- **K-mer Analysis**: k-mer counting, frequency spectra, saturation analysis, and distance metrics
- **Tool Integration**: Multiple assemblers, annotation tools (Pyrodigal, BLAST+, MMSeqs2), and QC utilities
- **Parallel Processing**: Multi-threaded analysis with progress tracking and HPC support (SLURM, rclone)
- **Visualization**: Coverage plots, k-mer spectra, embeddings, taxonomic distributions, and more

### Ongoing Areas of Reseach and Development

- **BioJulia-based Sequence Graph Assembly Framework**: From fixed-length (n-gram, k-mer, qualmer) to variable-length (string, FASTA, FASTQ) graphs
- **Quality-Aware Assembly**: A framework to preserve per-base quality scores throughout assembly process
- **Machine Learning Integration**: Reinforcement learning of optimal assembly workflow for automated parameter optimization

### Bioinformatics Workflows

- **Data Acquisition**: NCBI downloads, read simulation (Illumina, PacBio, Nanopore via ART, Badread)
- **Quality Control**: FastQC integration, native FASTQ analysis, filtering (fastp, filtlong, trim_galore)
- **Assembly**: External assemblers (MEGAHIT, metaSPAdes, SKESA, Flye, Canu, Hifiasm, Unicycler) plus novel Rhizomorph assembly framework
- **Annotation**: Gene prediction (Pyrodigal), homology search (BLAST+, MMSeqs2), specialized tools (tRNAscan-SE, TransTerm)
- **Alignment & Mapping**: Minimap2, Clustal Omega, BAM processing, variant calling
- **Comparative Genomics**: Pangenome analysis, FastANI integration, k-mer based comparisons

## Documentation

- **[Getting Started Guide](https://cjprybol.github.io/Mycelia/dev/getting-started/)** - Install and complete your first analysis
- **[API Reference](https://cjprybol.github.io/Mycelia/dev/)** - Complete function documentation
- **[Tutorials](https://cjprybol.github.io/Mycelia/dev/tutorials/)** - Step-by-step workflows
- **[Workflow & Tool Map](docs/src/workflow-map.md)** - Quick links from inputs to tools, outputs, and tutorials

## Installation

### Prerequisites

- Julia 1.10 or higher (LTS recommended)

### Install

```julia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
```

For detailed installation instructions including HPC setup, see the [Getting Started Guide](https://cjprybol.github.io/Mycelia/dev/getting-started/).

## Development Status

**Status**: Research platform with substantial implemented functionality alongside experimental algorithms

### Working Components

- **File Format Support**: FASTA/FASTQ/GenBank/GFF/VCF/SAM/BAM processing with compression support
- **Data Acquisition**: NCBI genome download, reference database access, read simulation (PacBio, Nanopore, Illumina)
- **Quality Control**: FastQC integration, comprehensive FASTQ analysis, filtering tools (fastp, filtlong, trim_galore)
- **Annotation Pipeline**: Pyrodigal, BLAST+, MMSeqs2, TransTerm, tRNAscan-SE, MLST integration
- **Alignment Tools**: Minimap2, Clustal Omega integration with variant calling support
- **Sequence Analysis**: K-mer counting, canonical k-mer analysis, sequence complexity assessment
- **Visualization**: Coverage plots, k-mer spectra, embeddings, taxonomic analysis, progress tracking

### Experimental/In Development

- **Rhizomorph Assembly Suite**: External assemblers (MEGAHIT, metaSPAdes, SKESA, Flye, Canu, Hifiasm, Unicycler) plus novel quality-aware graph algorithms
- **Graph-Based Assembly**: 6-graph type hierarchy with quality-aware assembly algorithms
- **Reinforcement Learning–Guided Assembly Optimization**
- **Advanced Assembly Validation Metrics**
- **Native Quality Control Implementations** (external tools currently integrated)
- **Pangenome Analysis Workflows**
- **Advanced Phylogenetics Integration**

### Known Limitations

- Some documentation examples may reference experimental features
- Research algorithms may require parameter tuning for optimal results

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- [Documentation](https://cjprybol.github.io/Mycelia/dev/)
- [Issue Tracker](https://github.com/cjprybol/Mycelia/issues)
- [Discussions](https://github.com/cjprybol/Mycelia/discussions)
