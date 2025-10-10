# Getting Started with Mycelia

Welcome to Mycelia, an experimental Julia package for bioinformatics and computational biology. This guide will help you install Mycelia and explore its current capabilities.

## What is Mycelia?

Mycelia is a research-oriented bioinformatics package that currently provides:

### Working Features
- **Basic Sequence I/O**: FASTA/FASTQ file reading and writing
- **Read Simulation**: PacBio and Nanopore read simulators  
- **Tool Integration**: Wrappers for established assemblers (MEGAHIT, SPAdes, hifiasm)
- **K-mer Analysis**: Canonical k-mer counting and distance calculations
- **HPC Support**: SLURM job submission and rclone integration

### Experimental Features
- **Novel Assembly Algorithms**: Graph-based approaches incorporating quality scores
- **Pangenome Analysis**: K-mer based comparative genomics
- **Quality Control**: Integration with external QC tools

### Planned Features  
- **Annotation**: Gene prediction integration
- **Phylogenetics**: Tree construction from pangenome data
- **Visualization**: Interactive genomic plots

## Prerequisites

### Julia Installation

Mycelia is tested against Julia LTS and release. Install Julia using [juliaup](https://github.com/JuliaLang/juliaup):

```bash
# Install juliaup
curl -fsSL https://install.julialang.org | sh

# Install latest lts Julia
juliaup add lts
juliaup default lts

# Install latest release Julia
juliaup add release
juliaup default release
```

### System Dependencies

Mycelia integrates with external bioinformatics tools. For full functionality, you'll need conda. If no conda environment is available, Mycelia will install it's own environment using Conda.jl.

### HPC Environment Setup

If you're using Mycelia on HPC systems, you may need to reset the `LD_LIBRARY_PATH` to avoid conflicts with visualization libraries:

```bash
# Launch Julia with clean library path
export LD_LIBRARY_PATH="" && julia
```

For Jupyter kernels, add this to your kernel.json:
```json
{
  "env": {
    "LD_LIBRARY_PATH": ""
  }
}
```

## Installation

### For Users

Install Mycelia directly from GitHub:

```julia
import Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
```

### For Developers

Clone and develop the package:

```julia
import Pkg
Pkg.develop(url="git@github.com:cjprybol/Mycelia.git")
```

## Your First Analysis

Let's walk through a complete workflow using small test datasets included with Mycelia.

### 1. Load Mycelia

```julia
import Mycelia
```

### 2. Download Test Data

Download a reference genome for testing:

```julia
# Download a small viral genome (phiX174)
genome_file = Mycelia.download_genome_by_accession(accession="NC_001422.1")

# Or create a random test sequence
test_genome = Mycelia.random_fasta_record(moltype=:DNA, L=10000)
Mycelia.write_fasta(outfile="test_genome.fasta", records=[test_genome])

# Simulate reads from the genome
reads_file = Mycelia.simulate_pacbio_reads(fasta="test_genome.fasta", quantity="20x")
Mycelia.write_fastq("test_reads.fastq", reads)
```

### 3. Quality Control (Using External Tools)

Filter reads using integrated external tools:

```julia
# Filter long reads with filtlong
Mycelia.qc_filter_long_reads_filtlong(
    input_file="test_reads.fastq",
    output_file="filtered_reads.fastq",
    min_length=1000,
    min_mean_q=7
)

# Note: Native quality assessment functions are planned but not yet implemented
```

### 4. K-mer Analysis

Analyze k-mer composition:

```julia
# Count canonical k-mers
import Kmers
kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{21}, "test_reads.fastq")
println("Unique k-mers: $(length(kmer_counts))")

# Note: K-mer spectrum analysis and plotting functions are planned
```

### 5. Genome Assembly  

Assemble using external tools through Mycelia wrappers:

```julia
# Use MEGAHIT for assembly (works with short reads)
Mycelia.assemble_metagenome_megahit(
    reads=["test_reads.fastq"],
    output_dir="megahit_assembly"
)

# Or use experimental graph-based assembly (research feature)
# Note: This is experimental and may not produce optimal results
# assembly = Mycelia.assemble_with_graph_framework("test_reads.fastq")

# Assembly evaluation functions are planned but not yet implemented
```

### 6. Comparative Analysis (Experimental)

Compare multiple genomes using k-mer analysis:

```julia
# Compare two genomes
genome_files = ["genome1.fasta", "genome2.fasta"]
pangenome_result = Mycelia.analyze_pangenome_kmers(
    genome_files,
    kmer_type=Kmers.DNAKmer{21}
)
println("Core k-mers: $(length(pangenome_result.core_kmers))")
println("Unique k-mers: $(sum(length(v) for v in values(pangenome_result.unique_kmers_by_genome)))")

# Note: Gene prediction and visualization functions are planned
```

## What's Next?

Explore the available features and help improve the package:

### Available Tutorials
- Tutorials are available but currently being reorganized
- See the generated tutorial pages in the documentation sidebar
- Focus on working examples that demonstrate current capabilities

### Research Features
- Explore the experimental graph-based assembly algorithms
- Test the quality-aware k-mer (qualmer) graph implementation
- Try the reinforcement learning guided assembly (very experimental)

### Contributing
- Report issues or feature requests on [GitHub](https://github.com/cjprybol/Mycelia/issues)
- Help implement planned features
- Improve documentation for existing functions

### Note on CLI Tools
Command-line interface tools are planned but not yet implemented.

### API Reference

Browse the [complete API documentation](api-reference.md) for detailed function references and examples.

## Getting Help

If you encounter issues:

1. Check the [FAQ](faq.md) for common issues
2. Browse the API reference for function documentation
3. Report bugs on [GitHub Issues](https://github.com/cjprybol/Mycelia/issues)

## Memory and Performance

For large-scale analyses:

```julia
# Check memory requirements
estimated_memory = Mycelia.estimate_memory_usage(input_file="large_reads.fastq")
println("Estimated memory needed: $(estimated_memory) GB")

# Monitor memory during analysis
Mycelia.with_memory_monitoring() do
    result = Mycelia.assemble_genome("large_reads.fastq")
end
```

## Contributing

Mycelia is open-source and welcomes contributions! Visit the [GitHub repository](https://github.com/cjprybol/Mycelia) for details on how to contribute.

---

*Ready to dive deeper? Explore the documentation sidebar for tutorials and guides on biological analyses with real datasets.*