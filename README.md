# Mycelia

**A comprehensive Julia package for bioinformatics and computational biology**

[![codecov](https://codecov.io/github/cjprybol/Mycelia/graph/badge.svg?token=0ZQSER2FLR)](https://codecov.io/github/cjprybol/Mycelia)
[![Documentation](https://github.com/cjprybol/Mycelia/actions/workflows/documentation.yml/badge.svg)](https://cjprybol.github.io/Mycelia/dev/)
[![Julia](https://img.shields.io/badge/julia-1.8+-blue.svg)](https://julialang.org/)

![Banner Logo](banner-logo.jpg)

## Overview

Mycelia provides a complete toolkit for modern bioinformatics workflows, from sequence processing to comparative genomics. Designed for researchers who need powerful, scalable tools for genomic analysis.

## 🚀 Quick Start

```julia
# Install Mycelia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")

# Run your first analysis
using Mycelia
genome = simulate_random_genome(length=10000)
reads = simulate_hifi_reads(genome, coverage=20)
assembly = assemble_genome(reads)
```

## ✨ Key Features

- **🧬 Sequence Processing**: FASTA/FASTQ handling, simulation, and quality control
- **🔧 Genome Assembly**: HiFi assembly with hifiasm, polishing, and error correction  
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
- Julia 1.8 or higher
- Optional: bioconda for external tools

### Install
```julia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
```

For detailed installation instructions including HPC setup, see the [Getting Started Guide](https://cjprybol.github.io/Mycelia/dev/getting-started/).

## 🧪 Development Status

**Status**: Active development - Core functionality stable, documentation expanding

- ✅ Core sequence processing and analysis
- ✅ Genome assembly and annotation pipelines  
- ✅ K-mer analysis and pangenome tools
- ✅ Visualization and plotting capabilities
- 🚧 Comprehensive user documentation (in progress)
- 🚧 Tutorial notebooks and workflows
- 🚧 HPC deployment guides

## 🤝 Contributing

Mycelia is open-source and welcomes contributions! See our [development guide](https://cjprybol.github.io/Mycelia/dev/contributing/) for details.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

- 📖 [Documentation](https://cjprybol.github.io/Mycelia/dev/)
- 🐛 [Issue Tracker](https://github.com/cjprybol/Mycelia/issues)
- 💬 [Discussions](https://github.com/cjprybol/Mycelia/discussions)
