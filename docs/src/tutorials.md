# Tutorials

The following tutorials aim to cover the common bioinformatics workflows from data acquisition to comparative genomics.

## Getting Started

- [Assembly in 5 Minutes](generated/tutorials/00_assembly_in_5_minutes.md) - Quick introduction to genome assembly with Mycelia

## Core Workflow Tutorials

1. [Data Acquisition](generated/tutorials/01_data_acquisition.md) - Downloading and simulating sequence data
2. [Quality Control](generated/tutorials/02_quality_control.md) - Quality assessment and filtering
3. [K-mer Analysis](generated/tutorials/03_kmer_analysis.md) - K-mer counting and spectrum analysis
4. [Genome Assembly](generated/tutorials/04_genome_assembly.md) - Assembly algorithms and execution
5. [Assembly Validation](generated/tutorials/05_assembly_validation.md) - Quality metrics and validation
6. [Gene Annotation](generated/tutorials/06_gene_annotation.md) - Gene prediction and annotation
7. [Comparative Genomics](generated/tutorials/07_comparative_genomics.md) - Pangenome and phylogenetic analysis
8. [Tool Integration](generated/tutorials/08_tool_integration.md) - Working with external bioinformatics tools
9. [Binning Workflow](generated/tutorials/14_binning_workflow.md) - Binning and post-binning workflows

## Specialized Topics

### Graph Types and Assembly Methods

- [Graph Type Tutorials](generated/tutorials/04_graph_type_tutorials.md) - Overview of different graph structures

### Round-Trip Examples

9. [String Graphs](generated/tutorials/09_round_trip_01_string_graphs.md) - String graph construction and traversal
10. [N-gram to String](generated/tutorials/09_round_trip_02_ngram_to_string.md) - Converting between representations
11. [FASTA Sequences](generated/tutorials/09_round_trip_03_fasta_sequences.md) - Working with FASTA format
12. [K-mer to Sequence](generated/tutorials/09_round_trip_04_kmer_to_sequence.md) - K-mer graph reconstruction
13. [FASTQ Graphs](generated/tutorials/09_round_trip_05_fastq_graphs.md) - Quality-aware graph construction
14. [Qualmer Graphs](generated/tutorials/09_round_trip_06_qualmer_graphs.md) - Quality-weighted k-mer graphs

## Advanced Topics

- [Advanced Assembly Theory and Practice](generated/tutorials/advanced-assembly-theory-and-practice.md) - Deep dive into assembly algorithms
- [Viroid Assembly Workflow](generated/tutorials/10_viroid_assembly_workflow.md) - Complete workflow for circular viral genomes

## Running Tutorials

All tutorials are executable Julia scripts. To run a tutorial:

```julia
import Mycelia
include("tutorials/01_data_acquisition.jl")
```

Or run all tutorials:

```julia
include("tutorials/run_all_tutorials.jl")
```

## Next Steps

- Review the [API Reference](api-reference.md) for detailed function documentation
- See [Getting Started](getting-started.md) for installation and basic usage
- Check [Examples](examples.md) for common workflow patterns
