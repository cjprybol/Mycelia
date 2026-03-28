# Mycelia Documentation

**An experimental Julia package for bioinformatics and computational biology**

Mycelia is a research-oriented package exploring best-practice methodologies for pan-multi-omics analysis and graph-based assembly with quality-aware sequence processing. Currently in early development, it provides both experimental algorithms and integrations with established bioinformatics tools.

## Quick Start

New to Mycelia? Start with our [Getting Started Guide](getting-started.md) to install the package and complete your first genomic analysis in minutes.

## Key Features & Research Areas

### Production-Ready Features
- **Third-Party Assembly Tools**: Stable wrappers for 15+ established assemblers including MEGAHIT, SPAdes, Flye, hifiasm, Canu, Unicycler, Velvet, and others
- **Sequence Processing**: FASTA/FASTQ I/O and read simulation (Illumina, PacBio, Nanopore)
- **K-mer Analysis**: Canonical k-mer counting and distance metrics
- **Quality Control Tools**: Integration with fastp, filtlong, trim_galore
- **Assembly Validation**: QUAST, BUSCO, CheckM/CheckM2, MUMmer integration
- **Annotation Tools**: Pyrodigal, BLAST+, MMSeqs2, TransTerm, tRNAscan-SE, MLST integration
- **HPC Support**: SLURM job submission and rclone cloud storage integration

### Experimental/Research Features
- **Internal Assembly Algorithms**: Graph-based approaches with quality awareness (6-graph hierarchy)
- **Intelligent Assembly**: Self-optimizing parameter selection and k-mer progression
- **Qualmer Analysis**: Quality-aware k-mer assembly and probabilistic path selection
- **Reinforcement Learning Assembly**: ML-guided assembly optimization (proof-of-concept)
- **Pangenome Analysis**: K-mer based comparative genomics
- **Phylogenetics**: Tree construction from pangenome data

## Documentation Contents

```@contents
Pages = [
    "getting-started.md",
    "installation.md",
    "concepts.md",
    "workflow-map.md",
    "metagenomic-workflow.md",
    "tutorials.md",
    "api/all-functions.md",
    "benchmarks.md",
    "related-projects.md",
    "contributing.md"
]
Depth = 2
```

## Navigating the Docs

- **Workflow & Tool Map**: See [`workflow-map.md`](workflow-map.md) for the end-to-end data-flow map with links to tutorials and third-party tool wrappers.
- **Tutorials**: Start with [`tutorials.md`](tutorials.md) for guided walkthroughs of each workflow stage.

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

For complete API documentation, see the [Complete API Surface](api/all-functions.md) section.
