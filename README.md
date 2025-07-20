# Mycelia

! Warning: In Development !

**A Julia package for bioinformatics and computational biology**

[![codecov](https://codecov.io/github/cjprybol/Mycelia/graph/badge.svg?token=0ZQSER2FLR)](https://codecov.io/github/cjprybol/Mycelia)
[![Documentation](https://github.com/cjprybol/Mycelia/actions/workflows/documentation.yml/badge.svg)](https://cjprybol.github.io/Mycelia/dev/)

![Banner Logo](banner-logo.jpg)

## Overview

Mycelia aims to provide a complete toolkit for modern bioinformatics workflows, from sequence processing to comparative genomics. Designed for researchers who need powerful, scalable tools for genomic analysis.

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

- **🧬 Sequence Processing**: FASTA/FASTQ handling, simulation, and quality control
- **🔧 Genome Assembly**: short-read, long-read, and hybrid assembly, polishing, and error correction  
- **📊 K-mer Analysis**: Quality-aware k-mer counting and graph construction
- **🌐 Pangenome Analysis**: Multi-genome comparative genomics
- **🔍 Annotation**: Gene prediction and functional annotation
- **🌳 Phylogenetics**: Tree construction and comparative analysis
- **📈 Visualization**: Interactive plots and data exploration
- **⚡ HPC Integration**: SLURM job submission and cloud storage

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

**Status**: Active development - Core functionality stable, documentation expanding

- ✅ Core sequence processing and analysis
- 🚧 K-mer analysis and pangenome tools
- 🚧 Genome assembly and annotation pipelines
- 🚧 Visualization and plotting capabilities
- 🚧 Comprehensive user documentation (in progress)
- 🚧 Tutorial notebooks and workflows
- 🚧 HPC deployment guides

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

- 📖 [Documentation](https://cjprybol.github.io/Mycelia/dev/)
- 🐛 [Issue Tracker](https://github.com/cjprybol/Mycelia/issues)
- 💬 [Discussions](https://github.com/cjprybol/Mycelia/discussions)
