# Mycelia Documentation

**An experimental Julia package for bioinformatics and computational biology**

Mycelia is a research-oriented package exploring novel approaches to genomic analysis, with a focus on graph-based genome assembly and quality-aware sequence processing. Currently in early development, it provides both experimental algorithms and integrations with established bioinformatics tools.

## Quick Start

New to Mycelia? Start with our [Getting Started Guide](getting-started.md) to install the package and complete your first genomic analysis in minutes.

## Key Features & Research Areas

### Currently Available
- **🧬 Sequence Processing**: Basic FASTA/FASTQ I/O and read simulation
- **📊 K-mer Analysis**: Canonical k-mer counting and distance metrics
- **🔧 Tool Integration**: Wrappers for established assemblers (MEGAHIT, SPAdes, hifiasm)
- **⚡ HPC Support**: SLURM job submission and rclone integration

### In Active Development  
- **🧪 Novel Assembly Algorithms**: Graph-based approaches with quality awareness
- **🌐 Pangenome Analysis**: K-mer based comparative genomics
- **📈 Quality Control**: Integration with QC tools (fastp, filtlong, trim_galore)

### Planned Features
- **🔍 Annotation**: Gene prediction and functional annotation
- **🌳 Phylogenetics**: Tree construction from pangenome data
- **📊 Visualization**: Interactive plots for genomic data

## Documentation Contents

```@contents
Pages = [
    "getting-started.md",
    "concepts.md",
    "tutorials.md", 
    "api.md",
    "visualization-gallery.md",
    "faq.md",
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

For complete API documentation, see the [API Reference](api-reference.md) section.
