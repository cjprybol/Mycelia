# Benchmarks

## Overview

Comprehensive benchmarks comparing Mycelia's various approaches on different datasets.

## Standardized Test Datasets

To ensure rigorous validation across platforms, Mycelia uses the following gold-standard communities:

### Mock Communities (Physical & Sequencing)
| Source | Product | Complexity | Description |
|--------|---------|------------|-------------|
| **Zymo** | [D6331](https://www.zymoresearch.com/products/zymobiomics-gut-microbiome-standard) | Medium | Gut Microbiome Standard (21 strains) |
| **Zymo** | [D6300](https://www.zymoresearch.com/products/zymobiomics-microbial-community-dna-standard) | Low | Microbial Community Standard (8 bacteria, 2 yeast) |
| **ATCC** | [MSA-1002](https://www.atcc.org/products/msa-1002) | Medium | 20 Strain Even Mix |
| **ATCC** | [MSA-1003](https://www.atcc.org/products/msa-1003) | Medium | 20 Strain Staggered Mix |
| **NIST** | [RM 8376](https://www.nist.gov/programs-projects/rm-8376-microbial-pathogen-dna-standards-detection-and-identification) | High | Microbial Pathogen DNA Standard |

### Benchmarking Challenges (Synthetic)
* **CAMI Challenge**: [Toy Datasets (Low/Med/High Complexity)](https://cami-challenge.org/datasets/)
* **Genome in a Bottle**: [HG002 (Ashkenazi Trio)](https://github.com/genome-in-a-bottle) - Standard for variant calling.

### Simulation Targets
For internal testing, we target the following simulation profiles:
* **Depth**: Low (10x), Medium (100x), High (1000x)
* **Diversity**: Isolate, Defined Community (10), Complex Community (100+)
* **Abundance**: Even, Random, Log-normal (staggered)

## Coming Soon

Detailed benchmark results including:
- Runtime comparisons
- Memory usage analysis
- Assembly quality metrics
- Accuracy assessments

For current benchmarking code and data, see the [benchmarking directory](https://github.com/cjprybol/Mycelia/tree/main/benchmarking) in the repository.

## Related Documentation

- [Workflow Map](workflow-map.md)
- [Metagenomic Workflow](metagenomic-workflow.md)
