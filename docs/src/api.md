# API Documentation

Welcome to the Mycelia API documentation! This guide organizes functions by **biological workflows** and **user tasks** rather than source code structure, making it easy to find the right tools for your analysis.

## 🧬 Quick Start

New to Mycelia? Start with our workflow-based guides:

- **[Basic Workflows](api/examples/basic-workflows.md)** - Common analysis patterns
- **[Function Index](api/quick-reference/function-index.md)** - Alphabetical function list
- **[Parameter Guide](api/quick-reference/parameter-guide.md)** - Common parameters explained

## 📋 By Workflow Stage

Follow the typical bioinformatics analysis workflow:

```@contents
Pages = [
    "api/workflows/data-acquisition.md",
    "api/workflows/quality-control.md", 
    "api/workflows/sequence-analysis.md",
    "api/workflows/genome-assembly.md",
    "api/workflows/assembly-validation.md",
    "api/workflows/gene-annotation.md",
    "api/workflows/comparative-genomics.md",
    "api/workflows/visualization.md"
]
Depth = 2
```

### 1. [Data Acquisition & Simulation](api/workflows/data-acquisition.md)
Download genomic data from public databases and simulate synthetic datasets for testing.

**Key Functions:** `download_genome_by_accession`, `simulate_hifi_reads`, `ncbi_genome_download`

### 2. [Quality Control & Preprocessing](api/workflows/quality-control.md)
Assess and improve sequencing data quality before analysis.

**Key Functions:** `analyze_fastq_quality`, `trim_reads`, `filter_sequences`

### 3. [Sequence Analysis & K-mers](api/workflows/sequence-analysis.md)
Analyze sequence composition, count k-mers, and extract genomic features.

**Key Functions:** `count_kmers`, `kmer_frequency_spectrum`, `estimate_genome_size`

### 4. [Genome Assembly](api/workflows/genome-assembly.md)
Assemble genomes from sequencing reads using modern HiFi approaches.

**Key Functions:** `assemble_genome`, `hifiasm_assembly`, `polish_assembly`

### 5. [Assembly Validation](api/workflows/assembly-validation.md)
Validate and assess the quality of genome assemblies.

**Key Functions:** `evaluate_assembly`, `calculate_assembly_stats`, `busco_analysis`

### 6. [Gene Annotation](api/workflows/gene-annotation.md)
Predict genes and assign functional annotations.

**Key Functions:** `predict_genes`, `annotate_functions`, `analyze_go_terms`

### 7. [Comparative Genomics](api/workflows/comparative-genomics.md)
Compare genomes, build pangenomes, and construct phylogenetic trees.

**Key Functions:** `build_pangenome`, `construct_phylogeny`, `calculate_synteny`

### 8. [Visualization & Reporting](api/workflows/visualization.md)
Create plots, figures, and reports for analysis results.

**Key Functions:** `plot_assembly_stats`, `plot_phylogenetic_tree`, `generate_report`

## 📁 By Data Type

Working with specific file formats and data structures:

```@contents
Pages = [
    "api/data-types/fasta-fastq.md",
    "api/data-types/assemblies.md",
    "api/data-types/annotations.md",
    "api/data-types/alignments.md",
    "api/data-types/trees.md"
]
Depth = 2
```

### [FASTA/FASTQ Files](api/data-types/fasta-fastq.md)
Reading, writing, and manipulating sequence files.

### [Assembly Files](api/data-types/assemblies.md)
Working with contigs, scaffolds, and assembly statistics.

### [Annotation Files](api/data-types/annotations.md)
Handling GFF3, GenBank, and other annotation formats.

### [Alignment Files](api/data-types/alignments.md)
Processing BAM/SAM files and alignment results.

### [Phylogenetic Trees](api/data-types/trees.md)
Tree construction, manipulation, and visualization.

## 🎯 By Analysis Goal

Cross-cutting concerns and specific use cases:

```@contents
Pages = [
    "api/examples/basic-workflows.md",
    "api/examples/advanced-usage.md",
    "api/quick-reference/function-index.md",
    "api/quick-reference/parameter-guide.md"
]
Depth = 2
```

### [Basic Workflows](api/examples/basic-workflows.md)
Complete examples for common analysis tasks.

### [Advanced Usage](api/examples/advanced-usage.md)
Complex workflows and optimization techniques.

### [Function Index](api/quick-reference/function-index.md)
Alphabetical listing of all functions with brief descriptions.

### [Parameter Guide](api/quick-reference/parameter-guide.md)
Common parameters and their usage across functions.

## 🔍 Finding What You Need

### By Task
- **"I want to assemble a genome"** → [Genome Assembly](api/workflows/genome-assembly.md)
- **"I need to validate my assembly"** → [Assembly Validation](api/workflows/assembly-validation.md)
- **"I want to compare genomes"** → [Comparative Genomics](api/workflows/comparative-genomics.md)
- **"I need to check data quality"** → [Quality Control](api/workflows/quality-control.md)

### By Data Type
- **"I have FASTQ files"** → [FASTA/FASTQ Files](api/data-types/fasta-fastq.md)
- **"I have assembly contigs"** → [Assembly Files](api/data-types/assemblies.md)
- **"I have gene annotations"** → [Annotation Files](api/data-types/annotations.md)

### By Experience Level
- **Beginner** → [Basic Workflows](api/examples/basic-workflows.md)
- **Intermediate** → Workflow-specific guides
- **Advanced** → [Advanced Usage](api/examples/advanced-usage.md)

## 💡 Usage Patterns

### Function Documentation Format
Each function is documented with:

```julia
"""
    function_name(required_param, optional_param="default")

Brief description of what the function does.

## Purpose
When and why to use this function in your workflow.

## Arguments
- `required_param`: Description and expected data type
- `optional_param`: Description, default value, and alternatives

## Returns
Description of return value and structure.

## Examples
```julia
# Basic usage
result = function_name("input.fasta")

# Advanced usage with options
result = function_name("input.fasta", 
                      optional_param="custom_value",
                      threads=4)
```

## Related Functions
- [`related_function`](@ref) - What it does
- [`workflow_next_step`](@ref) - Next step in workflow

## Performance Notes
- Memory usage: ~X GB for typical datasets
- Runtime: ~X minutes for Y-sized genomes
- Scaling: Linear/quadratic with input size

## See Also
- [Workflow Guide](../workflows/relevant-workflow.md)
- [Data Type Guide](../data-types/relevant-type.md)
"""
```

### Cross-References
Functions are linked to:
- **Workflow context** - Where they fit in analysis pipelines
- **Related functions** - What to use before/after
- **Data types** - What formats they accept/produce
- **Examples** - Real usage scenarios

## 🚀 Integration with Tutorials

This API documentation integrates with the tutorial system:

- **Tutorials** show complete workflows with explanation
- **API docs** provide detailed function reference
- **Examples** bridge the gap with focused use cases

For hands-on learning, see the [Tutorials](../tutorials.md) which use these functions in complete bioinformatics workflows.

---

*This documentation is automatically generated from function docstrings and organized for biological workflows. Functions are tested through the tutorial system to ensure accuracy.*