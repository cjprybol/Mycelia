# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Mycelia is a Julia package for bioinformatics and computational biology, providing tools for sequence analysis, genome assembly, annotation, and comparative genomics. The package is designed as a comprehensive bioinformatics toolkit with modular functionality.

## Development Commands

### Testing
```bash
# Run all tests from package root
julia --project=. -e "import Pkg; Pkg.test()"

# Run with coverage
julia --project=. --code-coverage=user -e "import Pkg; Pkg.test()"

# Run with color output
julia --project=. --color=yes -e "import Pkg; Pkg.test()"
```

### Static Analysis
```bash
# Run JET static analysis
julia --project=. test/jet.jl
```

### Documentation
```bash
# Build documentation (from docs/ directory)
julia --project=docs make.jl
```

## Architecture

### Core Structure
- **src/Mycelia.jl**: Main module file that dynamically imports all other source files
- **src/**: Modular functionality organized by domain:
  - `fastx.jl`: FASTA/FASTQ sequence handling
  - `assembly.jl`: Genome assembly tools
  - `annotation.jl`: Gene annotation functionality
  - `alignments-and-mapping.jl`: Sequence alignment tools
  - `kmer-analysis.jl`: K-mer counting and analysis
  - `clustering.jl`: Clustering algorithms
  - `classification.jl`: Sequence classification
  - `quality-control-and-benchmarking.jl`: QC and validation
  - `plotting-and-visualization.jl`: Data visualization
  - `utility-functions.jl`: General utility functions
  - `reference-databases.jl`: Database integration
  - `taxonomy-and-trees.jl`: Phylogenetic analysis
  - `variant-analysis.jl`: Variant calling and analysis

### Test Organization
Tests are organized in numbered directories following the bioinformatics workflow:
1. `1_data_acquisition/`: Data download and simulation
2. `2_preprocessing_qc/`: Quality control and preprocessing
3. `3_feature_extraction_kmer/`: K-mer analysis
4. `4_assembly/`: Genome assembly
5. `5_validation/`: Assembly validation
6. `6_annotation/`: Gene annotation
7. `7_comparative_pangenomics/`: Comparative genomics
8. `8_tool_integration/`: External tool integration

### Key Dependencies
- **BioSequences.jl**: Core sequence data structures
- **FASTX.jl**: FASTA/FASTQ parsing  
- **BioAlignments.jl**: Sequence alignment
- **Graphs.jl**: Graph algorithms for assembly

### Sequence Type Guidelines
- **ALWAYS use BioSequences types** - `BioSequences.LongDNA{4}`, `BioSequences.LongRNA{4}`, `BioSequences.LongAA`
- **Extract sequences from FASTQ records using proper types**: `FASTX.sequence(BioSequences.LongDNA{4}, record)`
- **String conversions should ONLY be used in string graphs** - everywhere else work with BioSequence objects
- **NO string conversions in k-mer graphs, qualmer graphs, or assembly algorithms** 
- **Use `string()` only when interfacing with external tools or final output**
- **DataFrames.jl**: Data manipulation
- **Makie.jl/Plots.jl**: Visualization
- **Arrow.jl/JLD2.jl**: Data serialization
- **XAM.jl**: BAM/SAM file handling

## Development Notes

### Module Loading
The main module uses dynamic file inclusion - all `.jl` files in `src/` are automatically included. When adding new functionality, create appropriately named files in `src/` and they will be automatically loaded.

### Dependency Management
- **All package dependencies are imported at the top-level** in `src/Mycelia.jl`
- **Individual source files should NOT import dependencies** - they are already available through the main module
- This pattern ensures consistent dependency management and avoids import conflicts

### Memory Management
The package includes utilities for memory estimation and checking (see `utility-functions.jl`). Large-scale genomic analyses should use these tools to avoid memory issues.

### Testing Approach
- Tests follow the bioinformatics workflow order
- Use `include_all_tests()` function to recursively load test files
- Test data is stored in `test/metadata/`

### Documentation
- Functions use Julia docstrings with DocStringExtensions
- Documentation is built using Documenter.jl
- Documentation source is in `docs/src/`

### Tutorial Files and Literate.jl
- Tutorial files in the `tutorials/` directory are processed by Literate.jl
- **Important**: Literate.jl comment conventions:
  - Single `#` at the start of a line becomes markdown text
  - Double `##` at the start of a line remains as a code comment
  - For inline code comments within code blocks, use `##` to prevent breaking the code block
  - Example:
    ```julia
    # This becomes markdown text
    
    open("file.txt") do io
        ## This remains a code comment
        line = readline(io)
        println(line)  ## This also remains a code comment
    end
    ```

## External Tool Integration

The package integrates with various bioinformatics tools:
- Bioconda for package management
- SLURM for job scheduling
- Rclone for cloud storage
<!-- - Neo4j for graph databases -->

## Code Style Guidelines

### Package Imports
- **NEVER use `using` statements or import specific functions** - ONLY import top-level packages with `import`
- **All package dependencies are imported at the top-level** in `src/Mycelia.jl` and are available in all source files
- **Individual source files should NOT import any packages** - they are already available through the main module
- Always use full module namespacing: `Test.@test`, `Dates.now()`, `Statistics.mean()`, etc.
- Example: Use `import Test` then `Test.@test`, NOT `using Test` or `using Test: @test`
- This ensures consistent dependency management and avoids import conflicts

## Communication and Documentation Standards

### Conservative and Understated Claims
- **Be conservative and understated in all statements and claims** made in this repository
- **Do not overpromise and under-deliver** - all claims must be verified facts with backing tests and benchmarks
- **Avoid overselling, overpromising, or unnecessarily over-hyping** what this software can do or what has been accomplished
- **Use measured language** - prefer "may provide", "can help with", "implements" over "revolutionary", "groundbreaking", "best-in-class"
- **Back all performance claims with actual benchmarks** - no theoretical improvements without validation
- **Clearly distinguish between completed features and planned features** in documentation
- **Use precise, technical language** rather than marketing language
- **Focus on functionality and capabilities** rather than superlative claims

## Testing Standards

### Never Disable Tests
- **NEVER disable or skip tests because underlying functionality is broken** - this masks bugs and prevents proper quality assurance
- **Fix the implementation first**, then ensure tests cover complete, expected, and logically correct behavior
- **Tests should define the expected behavior** - if a function fails tests, fix the function to meet the expected behavior
- **Use @test_skip only for features that are explicitly planned but not yet implemented**
- **All implemented functions must have working, comprehensive tests**