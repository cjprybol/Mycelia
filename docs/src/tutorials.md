# Tutorials

The following tutorials aim to cover the common bioinformatics workflows from data acquisition to comparative genomics.

## Getting Started

- [Rhizomorph Assembly](generated/tutorials/13_rhizomorph_assembly.md) - Quick introduction to genome assembly with Mycelia

## Core Workflow Tutorials

1. [01 Data Acquisition](generated/tutorials/01_data_acquisition.md) - Downloading and simulating sequence data
2. [02 Quality Control](generated/tutorials/02_quality_control.md) - Quality assessment and filtering
3. [03 K-mer Analysis](generated/tutorials/03_kmer_analysis.md) - K-mer counting and spectrum analysis
4. [04 Genome Assembly](generated/tutorials/04_genome_assembly.md) - Assembly algorithms and execution
5. [05 Assembly Validation](generated/tutorials/05_assembly_validation.md) - Quality metrics and validation
6. [06 Gene Annotation](generated/tutorials/06_gene_annotation.md) - Gene prediction and annotation
7. [07 Comparative Genomics](generated/tutorials/07_comparative_genomics.md) - Pangenome and phylogenetic analysis
8. [08 Tool Integration](generated/tutorials/08_tool_integration.md) - Working with external bioinformatics tools

## Specialized Topics

### Graph Types and Assembly Methods

4b. [04 Graph Type Tutorials](generated/tutorials/04_graph_type_tutorials.md) - Overview of different graph structures

### Round-Trip Examples

09a. [09 Round Trip 01: String Graphs](generated/tutorials/09_round_trip_01_string_graphs.md) - String graph construction and traversal
09b. [09 Round Trip 02: N-gram to String](generated/tutorials/09_round_trip_02_ngram_to_string.md) - Converting between representations
09c. [09 Round Trip 03: FASTA Sequences](generated/tutorials/09_round_trip_03_fasta_sequences.md) - Working with FASTA format
09d. [09 Round Trip 04: K-mer to Sequence](generated/tutorials/09_round_trip_04_kmer_to_sequence.md) - K-mer graph reconstruction
09e. [09 Round Trip 05: FASTQ Graphs](generated/tutorials/09_round_trip_05_fastq_graphs.md) - Quality-aware graph construction
09f. [09 Round Trip 06: Qualmer Graphs](generated/tutorials/09_round_trip_06_qualmer_graphs.md) - Quality-weighted k-mer graphs

## Advanced Topics

10. [10 Viroid Assembly Workflow](generated/tutorials/10_viroid_assembly_workflow.md) - Complete workflow for circular viral genomes
11. [11 Reduced Amino Acid Alphabets](generated/tutorials/11_reduced_amino_acid_alphabets.md) - Reduced alphabet encodings and analysis
12. [12 CoverM Coverage](generated/tutorials/12_coverm_coverage.md) - Coverage profiling and reporting
13. [13 Rhizomorph Assembly](generated/tutorials/13_rhizomorph_assembly.md) - Rhizomorph assembly workflows
14a. [14 Binning Workflow](generated/tutorials/14_binning_workflow.md) - Binning and post-binning workflows
14b. [14 Mash Classification](generated/tutorials/14_mash_classification.md) - Mash sketching and screening
15. [15 Round Trip Benchmarking](generated/tutorials/15_round_trip_benchmarking.md) - Reference download and simulation harness
16. [16 UN Corpus Ngram vs Token Graphs](generated/tutorials/16_un_corpus_ngram_vs_token_graphs.md) - Graph comparisons on UN corpus data
17. [17 Viroid Sketch Round Trip](generated/tutorials/17_viroid_sketch_round_trip.md) - BLAST DB export, sketching, and validation
18. [18 Advanced Assembly Theory and Practice](generated/tutorials/18_advanced_assembly_theory_and_practice.md) - Deep dive into assembly algorithms

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

- Review the [Complete API Surface](api/all-functions.md) for detailed function documentation
- See [Getting Started](getting-started.md) for installation and basic usage
- Check the [Workflow Map](workflow-map.md) for common workflow patterns
