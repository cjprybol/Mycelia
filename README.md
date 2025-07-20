# Mycelia

! Warning: In Development !

**A Julia package for bioinformatics and computational biology**

[![codecov](https://codecov.io/github/cjprybol/Mycelia/graph/badge.svg?token=0ZQSER2FLR)](https://codecov.io/github/cjprybol/Mycelia)
[![Documentation](https://github.com/cjprybol/Mycelia/actions/workflows/documentation.yml/badge.svg)](https://cjprybol.github.io/Mycelia/dev/)

![Banner Logo](banner-logo.jpg)

## Overview

Mycelia is an experimental Julia package exploring novel approaches to bioinformatics, with a focus on graph-based genome assembly and quality-aware sequence analysis. Currently in early development, the package provides a research platform for testing innovative algorithms while building toward a comprehensive toolkit.

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

## ✨ Key Features & Research Areas

### Currently Available
- **🧬 Sequence Processing**: Basic FASTA/FASTQ I/O and read simulation
- **📊 K-mer Analysis**: Canonical k-mer counting and distance metrics
- **🔧 Tool Integration**: Wrappers for MEGAHIT, SPAdes, and other assemblers
- **⚡ HPC Support**: SLURM job submission and rclone integration

### In Active Development
- **🧪 Novel Assembly Algorithms**: Graph-based approaches with quality awareness
- **🌐 Pangenome Analysis**: K-mer based comparative genomics (basic implementation available)
- **📈 Quality Control**: Comprehensive QC pipeline (partially implemented)

### Planned/Early Stage
- **🔍 Annotation**: Gene prediction and functional annotation integration
- **🌳 Phylogenetics**: Tree construction from pangenome data
- **📊 Visualization**: Interactive plots for genomic data
- **🔧 Assembly Polish**: Error correction and consensus calling

## 📚 Documentation

- **[Getting Started Guide](https://cjprybol.github.io/Mycelia/dev/getting-started/)** - Install and complete your first analysis
- **[API Reference](https://cjprybol.github.io/Mycelia/dev/)** - Complete function documentation
- **[Tutorials](https://cjprybol.github.io/Mycelia/dev/tutorials/)** - Step-by-step workflows

## 🔧 Installation

### Prerequisites
- Julia LTS or higher

### Install
```julia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
```

For detailed installation instructions including HPC setup, see the [Getting Started Guide](https://cjprybol.github.io/Mycelia/dev/getting-started/).

## 🧪 Development Status

**Status**: Early development - Research algorithms being implemented and tested

### Working Components
- ✅ Basic FASTA/FASTQ file I/O
- ✅ External tool integration (MEGAHIT, metaSPAdes, badread)
- ✅ Read simulation functions (PacBio, Nanopore)
- ✅ Data download by accession number
- ✅ K-mer counting infrastructure

### Experimental/In Development
- 🧪 Graph-based assembly framework (6-graph hierarchy)
- 🧪 Quality-aware sequence graphs (qualmer graphs)
- 🧪 Machine learning guided assembly (reinforcement learning)
- 🚧 Quality control and filtering functions
- 🚧 Assembly validation and metrics
- 🚧 Pangenome analysis workflows
- 🚧 Visualization capabilities

### Known Limitations
- ⚠️ Testing framework needs fixes (`Pkg.test()` currently fails)
- ⚠️ Many documented functions not yet implemented
- ⚠️ Documentation may reference planned features
- ⚠️ Not recommended for production workflows

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

- 📖 [Documentation](https://cjprybol.github.io/Mycelia/dev/)
- 🐛 [Issue Tracker](https://github.com/cjprybol/Mycelia/issues)
- 💬 [Discussions](https://github.com/cjprybol/Mycelia/discussions)
