# Mycelia Documentation

**A comprehensive Julia package for bioinformatics and computational biology**

Mycelia provides a complete toolkit for genomic analysis, from sequence processing to comparative genomics. Built for researchers who need powerful, scalable tools for modern bioinformatics workflows.

## Quick Start

New to Mycelia? Start with our [Getting Started Guide](getting-started.md) to install the package and complete your first genomic analysis in minutes.

## Key Features

- **🧬 Sequence Processing**: FASTA/FASTQ handling, simulation, and quality control
- **🔧 Genome Assembly**: HiFi assembly with hifiasm, polishing, and error correction  
- **📊 K-mer Analysis**: Quality-aware k-mer counting and graph construction
- **🌐 Pangenome Analysis**: Multi-genome comparative genomics
- **🔍 Annotation**: Gene prediction and functional annotation
- **🌳 Phylogenetics**: Tree construction and comparative analysis
- **📈 Visualization**: Interactive plots and data exploration
- **⚡ HPC Integration**: SLURM job submission and cloud storage

## Documentation Contents

```@contents
Pages = [
    "getting-started.md",
    "concepts.md",
    "tutorials.md", 
    "api.md",
    "visualization-gallery.md",
    "troubleshooting.md",
    "contributing.md"
]
Depth = 2
```

## Installation

### Quick Install
```julia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
```

### Development Install
```julia
import Pkg
Pkg.develop(url="git@github.com:cjprybol/Mycelia.git")
```

For detailed installation instructions including HPC setup, see the [Getting Started Guide](getting-started.md).

## Function Docstrings

```@autodocs
Modules = [Mycelia]
```
