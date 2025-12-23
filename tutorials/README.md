# Mycelia Tutorials

This directory contains comprehensive bioinformatics tutorials using Literate.jl format. Each tutorial can be run as a Julia script or converted to a Jupyter notebook.

## Tutorial Structure

Tutorials are organized following the standard bioinformatics workflow:

1. **Data Acquisition** - Downloading and simulating genomic data
2. **Quality Control** - Preprocessing and quality assessment
3. **K-mer Analysis** - Sequence feature extraction and analysis
4. **Genome Assembly** - Constructing genomes from sequencing data
5. **Assembly Validation** - Quality assessment and validation
6. **Gene Annotation** - Identifying and annotating genomic features
7. **Comparative Genomics** - Pangenome and phylogenetic analysis
8. **Tool Integration** - Working with external bioinformatics tools
9. **Binning & Post-binning** - MAG binning and dereplication workflows

## Usage

### Running as Julia Scripts

Each tutorial can be run directly as a Julia script:

```bash
julia --project=. tutorials/01_data_acquisition.jl
```

### Converting to Jupyter Notebooks

Convert any tutorial to a Jupyter notebook using Literate.jl:

```bash
julia --project=. -e 'import Literate; Literate.notebook("tutorials/01_data_acquisition.jl", "tutorials/notebooks", execute=false)'
```

### Running All Tutorials

Execute all tutorials (use with caution - may take significant time):

```bash
julia --project=. tutorials/run_all_tutorials.jl
```

## Tutorial Variants

Each tutorial includes:
- **Quick Version**: Small test datasets for rapid learning
- **Realistic Version**: Larger datasets for practical experience
- **Performance Analysis**: Benchmarking and optimization examples

## Requirements

- Julia 1.8+
- Mycelia package installed
- Internet connection for data download
- Optional: External bioinformatics tools (bioconda recommended)

## Target Audience

These tutorials are designed for graduate-level and post-graduate researchers in:
- Bioinformatics
- Computational biology
- Genomics
- Microbiology
- Evolutionary biology

## Getting Help

- Check the main [Getting Started Guide](../docs/src/getting-started.md)
- Browse the [API Documentation](../docs/src/index.md)
- Report issues on [GitHub](https://github.com/cjprybol/Mycelia/issues)
