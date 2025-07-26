# Mycelia

! Warning: In Development !

**A Julia package for bioinformatics and computational biology**

[![codecov](https://codecov.io/github/cjprybol/Mycelia/graph/badge.svg?token=0ZQSER2FLR)](https://codecov.io/github/cjprybol/Mycelia)
[![Documentation](https://github.com/cjprybol/Mycelia/actions/workflows/documentation.yml/badge.svg)](https://cjprybol.github.io/Mycelia/dev/)

![Banner Logo](banner-logo.jpg)

## Overview

Mycelia is a Julia package for bioinformatics and computational biology that implements graph-based genome assembly and quality-aware sequence analysis. The package provides both research-oriented algorithms and practical bioinformatics functionality, including extensive tool integration for genomics workflows. While some components are experimental and in active development, the package includes substantial implemented functionality for data processing, assembly, annotation, and analysis.

## 🚀 Quick Start

```julia
# Install Mycelia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")

# Run your first analysis
import Mycelia

# Download a reference genome (phiX174 bacteriophage)
genome_file = Mycelia.download_genome_by_accession(accession="NC_001422.1")

# Generate test reads from the reference
reads_file = Mycelia.simulate_pacbio_reads(fasta=genome_file, quantity="50x")

# Assemble the genome
assembly = Mycelia.assemble_genome(reads_file)
```

## ✨ Key Features

### Core Functionality
- **🧬 Comprehensive File Format Support**: FASTA/FASTQ/GenBank/GFF/VCF/SAM/BAM with automatic compression handling
- **📊 Advanced K-mer Analysis**: Canonical k-mer counting, frequency spectra, saturation analysis, and distance metrics
- **🔧 Extensive Tool Integration**: 30+ assemblers, annotation tools (Pyrodigal, BLAST+, MMSeqs2), and QC utilities
- **⚡ Parallel Processing**: Multi-threaded analysis with progress tracking and HPC support (SLURM, rclone)
- **📈 Rich Visualization**: Coverage plots, k-mer spectra, embeddings, taxonomic distributions, and more

### Research Innovations
- **🧪 Novel 6-Graph Assembly Framework**: Unique hierarchy transitioning from fixed-length (k-mer, qualmer) to variable-length (string, FASTQ) graphs
- **🎯 Quality-Aware Assembly**: First framework to preserve per-base quality scores throughout assembly process
- **🤖 Machine Learning Integration**: Four reinforcement learning approaches (DQN, PPO, POMDPs, MCTS) for automated parameter optimization
- **🔬 Zero String Conversion**: Type-safe implementation using native BioSequences types throughout

### Bioinformatics Workflows
- **🧬 Data Acquisition**: NCBI downloads, read simulation (Illumina, PacBio, Nanopore via ART, Badread)
- **🔍 Quality Control**: FastQC integration, native FASTQ analysis, filtering (fastp, filtlong, trim_galore)
- **🧩 Assembly**: External assemblers (MEGAHIT, metaSPAdes, SKESA, Flye, Canu, Hifiasm, Unicycler) plus novel Rhizomorph assembly framework
- **🏷️ Annotation**: Gene prediction (Pyrodigal), homology search (BLAST+, MMSeqs2), specialized tools (tRNAscan-SE, TransTerm)
- **🧮 Alignment & Mapping**: Minimap2, Clustal Omega, BAM processing, variant calling
- **📊 Comparative Genomics**: Pangenome analysis, FastANI integration, k-mer based comparisons

## 📚 Documentation

- **[Getting Started Guide](https://cjprybol.github.io/Mycelia/dev/getting-started/)** - Install and complete your first analysis
- **[API Reference](https://cjprybol.github.io/Mycelia/dev/)** - Complete function documentation
- **[Tutorials](https://cjprybol.github.io/Mycelia/dev/tutorials/)** - Step-by-step workflows

## 🔧 Installation

### Prerequisites
- Julia 1.10 or higher (LTS recommended)

### Install
```julia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
```

For detailed installation instructions including HPC setup, see the [Getting Started Guide](https://cjprybol.github.io/Mycelia/dev/getting-started/).

## 🧪 Development Status

**Status**: Research platform with substantial implemented functionality alongside experimental algorithms

### Working Components
- ✅ **File Format Support**: FASTA/FASTQ/GenBank/GFF/VCF/SAM/BAM processing with compression support
- ✅ **Data Acquisition**: NCBI genome download, reference database access, read simulation (PacBio, Nanopore, Illumina)
- ✅ **Quality Control**: FastQC integration, comprehensive FASTQ analysis, filtering tools (fastp, filtlong, trim_galore)
- ✅ **Rhizomorph Assembly Suite**: External assemblers (MEGAHIT, metaSPAdes, SKESA, Flye, Canu, Hifiasm, Unicycler) plus novel quality-aware graph algorithms
- ✅ **Annotation Pipeline**: Pyrodigal, BLAST+, MMSeqs2, TransTerm, tRNAscan-SE, MLST integration
- ✅ **Alignment Tools**: Minimap2, Clustal Omega integration with variant calling support
- ✅ **Sequence Analysis**: K-mer counting, canonical k-mer analysis, sequence complexity assessment
- ✅ **Graph-Based Assembly**: Novel 6-graph type hierarchy with quality-aware assembly algorithms
- ✅ **Visualization**: Coverage plots, k-mer spectra, embeddings, taxonomic analysis, progress tracking

### Experimental/In Development
- 🧪 Reinforcement learning guided assembly optimization (Custom RL, ReinforcementLearning.jl, POMDPs.jl, Monte Carlo Tree Search)
- 🧪 Advanced assembly validation metrics
- 🚧 Native quality control implementations (external tools currently integrated)
- 🚧 Pangenome analysis workflows
- 🚧 Advanced phylogenetics integration

### Known Limitations
- ⚠️ Testing framework needs fixes (`Pkg.test()` currently fails due to dependency conflicts)
- ⚠️ Function discoverability requires using `Mycelia.function_name()` syntax
- ⚠️ Some documentation examples may reference experimental features
- ⚠️ Research algorithms may require parameter tuning for optimal results

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

- 📖 [Documentation](https://cjprybol.github.io/Mycelia/dev/)
- 🐛 [Issue Tracker](https://github.com/cjprybol/Mycelia/issues)
- 💬 [Discussions](https://github.com/cjprybol/Mycelia/discussions)
