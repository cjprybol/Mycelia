# Installation Guide

This comprehensive guide covers installing Mycelia on various platforms and environments.

## Quick Start

For most users, the simplest installation method is:

```julia
using Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
```

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Prerequisites](#prerequisites)
3. [Installation Methods](#installation-methods)
4. [Platform-Specific Instructions](#platform-specific-instructions)
5. [HPC and Cluster Setup](#hpc-and-cluster-setup)
6. [Docker Installation](#docker-installation)
7. [Development Installation](#development-installation)
8. [Verifying Installation](#verifying-installation)
9. [Common Issues](#common-issues)

## System Requirements

### Minimum Requirements
- **Julia**: Version 1.10 or higher
- **RAM**: 8 GB (for small datasets)
- **Storage**: 2 GB free space
- **OS**: Linux

### Recommended Requirements
- **Julia**: Latest stable version
- **RAM**: 32 GB or more (for genomic datasets)
- **Storage**: 50 GB free space
- **CPU**: Multi-core processor for parallel operations

### For Large-Scale Assembly
- **RAM**: 64-256 GB depending on genome size
- **Storage**: SSD with 500 GB+ free space
- **CPU**: 16+ cores recommended

## Prerequisites

### Julia Installation

1. **Download Julia** from [julialang.org](https://julialang.org/downloads/)
2. **Install Julia** following platform-specific instructions
3. **Verify installation**:
   ```bash
   julia --version
   ```

### Optional Dependencies

Some features require additional software:

- **Bioconda tools** (for external tool integration)
- **Git** (for development installation)
- **C compiler** (for some dependencies)

## Installation Methods

### Method 1: Direct Package Installation (Recommended)

```julia
using Pkg
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
```

### Method 2: Development Installation

```bash
git clone https://github.com/cjprybol/Mycelia.git
cd Mycelia
julia --project=.
```

In Julia REPL:
```julia
using Pkg
Pkg.instantiate()
```

## Verifying Installation

### Run Test Suite
```julia
using Pkg
Pkg.test("Mycelia")
```

### Quick Assembly Test
```julia
using Mycelia

# Download test data
ref_file = Mycelia.download_genome_by_accession(accession="NC_001422.1")

# Simulate reads
reads_file = Mycelia.simulate_pacbio_reads(
    fasta=ref_file,
    quantity="10x"
)

# Run assembly
results = Mycelia.mycelia_assemble(reads_file, max_k=31)
```

## Common Issues

### Issue: Package Installation Fails

**Symptom**: Error during `Pkg.add()`

**Solutions**:
1. Update Julia: `juliaup update` (if using juliaup)
2. Clear package cache: `rm -rf ~/.julia/compiled`
3. Update registry: `Pkg.Registry.update()`

### Issue: Missing Dependencies

**Symptom**: Package load errors

**Solutions**:
```julia
using Pkg
Pkg.resolve()
Pkg.instantiate()
```

## Next Steps

- Read the [Getting Started](getting-started.md) guide
- Try the [5-minute assembly tutorial](../tutorials/00_assembly_in_5_minutes.jl)
- Explore the [API Reference](api.md)
- Review [Performance Guide](performance.md) for optimization tips