# Installation Guide

This comprehensive guide covers installing Mycelia on various platforms and environments.

## Quick Start

For most users, the simplest installation method is:

```julia
import Pkg
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

### Recommended Requirements
- **Julia**: 1.10 (LTS) or newer (LTS recommended)
- **RAM**: 32 GB or more (for genomic datasets)
- **Storage**: 50 GB free space
- **CPU**: Multi-core processor for parallel operations

### For Large-Scale Assembly
- **RAM**: 64-256 GB depending on genome size
- **Storage**: 500 GB+ free space
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

## Installation Methods

### Method 1: Direct Package Installation (Recommended)

```julia
import Pkg
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
import Pkg
Pkg.instantiate()
```

## Verifying Installation

### Run Test Suite
```julia
import Pkg
Pkg.test("Mycelia")
```

### Quick Assembly Test
```julia
import Mycelia

# Download test data
ref_file = Mycelia.download_genome_by_accession(accession="NC_001422.1")

# Simulate reads
reads_file = Mycelia.simulate_pacbio_reads(
    fasta=ref_file,
    quantity="10x"
)
```

## HPC and Cluster Setup

Mycelia includes SLURM helpers and an HPC CI runner for extended tests, tutorials, and benchmarks.
Start with the local guide in `ci/hpc/README.md` for scheduler examples and resource sizing.

If you see shared library conflicts (common on clusters), launch Julia with a clean library path:
```bash
export LD_LIBRARY_PATH="" && julia
```

## Next Steps

- Read the [Getting Started](getting-started.md) guide
- Explore the [Complete API Surface](api/all-functions.md)
