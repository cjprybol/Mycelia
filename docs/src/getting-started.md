# Getting Started with Mycelia

Welcome to Mycelia, a comprehensive Julia package for bioinformatics and computational biology. This guide will help you install Mycelia and complete your first genomic analysis.

## What is Mycelia?

Mycelia is a modular bioinformatics toolkit that provides:

- **Sequence Processing**: FASTA/FASTQ handling, simulation, and quality control
- **Genome Assembly**: HiFi assembly with hifiasm, polishing, and error correction
- **K-mer Analysis**: Quality-aware k-mer counting and graph construction
- **Pangenome Analysis**: Multi-genome comparative genomics
- **Annotation**: Gene prediction and functional annotation
- **Phylogenetics**: Tree construction and comparative analysis
- **Visualization**: Interactive plots and data exploration
- **HPC Integration**: SLURM job submission and cloud storage

## Prerequisites

### Julia Installation

Mycelia requires Julia 1.8 or higher. Install Julia using [juliaup](https://github.com/JuliaLang/juliaup):

```bash
# Install juliaup
curl -fsSL https://install.julialang.org | sh

# Install latest stable Julia
juliaup add release
juliaup default release
```

### System Dependencies

Mycelia integrates with external bioinformatics tools. For full functionality, you'll need:

```bash
# Install bioconda (recommended for tool management)
conda install -c conda-forge mamba
mamba install -c bioconda hifiasm minimap2 pyrodigal
```

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
using Mycelia
```

### 2. Simulate Test Data

Create synthetic genomic data for testing:

```julia
# Generate random genome sequence
genome = simulate_random_genome(length=10000, gc_content=0.45)

# Simulate HiFi reads
reads = simulate_hifi_reads(genome, coverage=20, error_rate=0.001)

# Write to FASTQ format
write_fastq("test_reads.fastq", reads)
```

### 3. Quality Control

Assess read quality and characteristics:

```julia
# Load reads and calculate statistics
read_stats = analyze_fastq_quality("test_reads.fastq")
println("Read count: $(read_stats.n_reads)")
println("Mean length: $(read_stats.mean_length)")
println("Mean quality: $(read_stats.mean_quality)")
```

### 4. K-mer Analysis

Analyze k-mer composition:

```julia
# Count k-mers
kmer_counts = count_kmers("test_reads.fastq", k=21)

# Analyze k-mer frequency spectrum
spectrum = kmer_frequency_spectrum(kmer_counts)
plot_kmer_spectrum(spectrum)
```

### 5. Genome Assembly

Assemble the genome from reads:

```julia
# Run hifiasm assembly
assembly = assemble_genome("test_reads.fastq", 
                          output_dir="assembly_output",
                          threads=4)

# Assess assembly quality
assembly_stats = evaluate_assembly(assembly)
println("Contigs: $(assembly_stats.n_contigs)")
println("N50: $(assembly_stats.n50)")
println("Total length: $(assembly_stats.total_length)")
```

### 6. Annotation

Predict genes in the assembled genome:

```julia
# Predict genes with Pyrodigal
genes = predict_genes(assembly, genetic_code="standard")
println("Predicted genes: $(length(genes))")

# Write annotations to GFF3
write_gff3("annotations.gff3", genes)
```

### 7. Visualization

Create plots to visualize your results:

```julia
# Plot assembly statistics
plot_assembly_stats(assembly_stats)

# Create k-mer spectrum plot
plot_kmer_spectrum(spectrum)

# Visualize gene features
plot_genome_features(assembly, genes)
```

## What's Next?

Now that you've completed your first analysis, explore these advanced topics:

### Tutorials (Coming Soon)
- **Genome Assembly Tutorial**: Complete HiFi assembly workflow
- **Pangenome Analysis**: Multi-genome comparative analysis
- **Phylogenetic Analysis**: Tree construction and visualization
- **HPC Deployment**: Running large-scale analyses on clusters

### CLI Usage

Mycelia includes command-line tools for common workflows:

```bash
# Run complete assembly pipeline
mycelia assemble --input reads.fastq --output assembly/ --threads 8

# Perform k-mer analysis
mycelia kmer-analysis --input reads.fastq --k 21 --output kmers/

# Annotate genome
mycelia annotate --genome assembly.fasta --output annotations.gff3
```

### API Reference

Browse the [complete API documentation](api.md) for detailed function references and examples.

## Getting Help

If you encounter issues:

1. Check the [troubleshooting guide](troubleshooting.md)
2. Browse [example workflows](examples.md)
3. Report bugs on [GitHub Issues](https://github.com/cjprybol/Mycelia/issues)

## Memory and Performance

For large-scale analyses:

```julia
# Check memory requirements
estimated_memory = estimate_memory_usage(input_file="large_reads.fastq")
println("Estimated memory needed: $(estimated_memory) GB")

# Monitor memory during analysis
with_memory_monitoring() do
    result = assemble_genome("large_reads.fastq")
end
```

## Contributing

Mycelia is open-source and welcomes contributions! See our [development guide](contributing.md) for details.

---

*Ready to dive deeper? Check out our [workflow tutorials](tutorials.md) for complete biological analyses with real datasets.*