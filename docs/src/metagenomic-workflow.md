# Metagenomic Analysis Workflow

This document describes the metagenomic analysis workflow implemented in Mycelia, including read-level classification, assembly, taxonomy assignment, and validation strategies.

## Overview

The Mycelia metagenomic workflow implements a triangulated approach to taxonomic classification, combining evidence from multiple sources:

1. **Read-level classification** - Direct mapping of reads to reference databases
2. **Contig-level DNA taxonomy** - BLAST-based classification of assembled contigs
3. **Contig-level protein taxonomy** - MMseqs2-based classification against UniRef databases

```mermaid
flowchart TD
  R[Sequencing reads]

  R --> E1[Read level mapping<br/>to reference databases]

  R --> Asm[Assembly into contigs]
  Asm --> E2[Contig DNA taxonomy<br/>BLAST to NCBI and others]
  Asm --> E3[Contig protein taxonomy<br/>mmseqs2 to UniRef sets]

  R --> MapAsm[Reads mapped to contigs]
  MapAsm --> E2
  MapAsm --> E3

  E1 --> Tri[Triangulated taxonomic assignment]
  E2 --> Tri
  E3 --> Tri

  Tri --> Bench[Benchmark against other tools<br/>kmer marker gene protein based]
```

## Full Workflow

The complete metagenomic workflow includes quality control, multiple assembly strategies, variant calling, and pangenome analysis:

```mermaid
flowchart TD
  R[Raw reads] --> Q[QC]

  Q --> DBs[Reference databases<br/>NCBI Core Nucleotide<br/>IMG VR v5 and others]

  Q --> RL_map
  subgraph RL [Read level classification]
    direction TB
    RL_map[Minimap2 mapping<br/>short and long reads] --> RL_score[Weighted scoring and LCA<br/>top hit vs next hit<br/>coverage threshold]
  end
  DBs --> RL_map

  Q --> RL_val[Validation tools<br/>MetaPhlAn<br/>metabuli<br/>Sylph or Kraken2]

  Q --> SR_asm
  Q --> LR_asm
  subgraph ASM [Metagenomic assembly]
    direction LR
    SR_asm[Short read assembly<br/>MEGAHIT<br/>metaSPAdes<br/>Penguinn]
    LR_asm[Long read assembly<br/>hifiasm meta<br/>metaMDBG<br/>metaFlye]
  end

  SR_asm --> DNA_tax
  LR_asm --> DNA_tax
  DBs --> DNA_tax
  DNA_tax[DNA level taxonomy<br/>BLAST vs CoreNT and others]

  SR_asm --> Prot_tax
  LR_asm --> Prot_tax
  Prot_tax[Protein level taxonomy<br/>mmseqs2 vs UniRef 50 90 100]

  Q --> MapAsm[Reads mapped to contigs<br/>minimap2]
  MapAsm --> Pangraph[Pangenome and graph alignment<br/>metapangenomes]
  Pangraph --> VarPan[Variant calling graph based<br/>PGGB<br/>VG<br/>Cactus]
  Pangraph --> VarAlign[Variant calling alignment based<br/>short reads GATK FreeBayes BCFtools<br/>long reads Clair3]

  MapAsm --> DNA_tax
  MapAsm --> Prot_tax

  SR_asm --> QCcontig
  LR_asm --> QCcontig
  QCcontig[Contig QC<br/>viral CheckV<br/>microbial CheckM<br/>all QUAST BUSCO]

  RL_score --> Tri[Triangulated taxonomic assignment]
  DNA_tax --> Tri
  Prot_tax --> Tri

  Tri --> Bench[Benchmark against other tools<br/>kmer marker gene protein based]
  Bench --> RL_val

  Tri --> C[Classifications]
  VarPan --> V[Variants]
  VarAlign --> V[Variants]
  Pangraph --> P[Pangenomes]
```

## Workflow Components

### Quality Control

Initial read quality assessment and filtering using:
- FastQC for quality metrics
- Adapter trimming and quality filtering

### Read-Level Classification

Direct classification of reads using minimap2 mapping against reference databases, with weighted scoring based on:
- Top hit vs next hit ratio
- Coverage thresholds
- LCA (Lowest Common Ancestor) assignment for ambiguous mappings

### Assembly Strategies

**Short Read Assembly:**
- MEGAHIT - memory-efficient assembly
- metaSPAdes - metagenome-specific assembly
- Penguinn - guided assembly with protein references

**Long Read Assembly:**
- hifiasm-meta - HiFi metagenome assembly
- metaMDBG - minimizer-based assembly
- metaFlye - long-read metagenome assembly

### Taxonomy Assignment

**DNA-level:**
- BLAST against NCBI Core Nucleotide
- Custom reference databases

**Protein-level:**
- MMseqs2 against UniRef50, UniRef90, UniRef100
- Sensitive homology detection for divergent sequences

### Quality Control and Validation

**Viral contigs:**
- CheckV for completeness and contamination

**Microbial contigs:**
- CheckM/CheckM2 for completeness and contamination

**General:**
- QUAST for assembly statistics
- BUSCO for completeness assessment

### Variant Calling

**Graph-based:**
- PGGB (PanGenome Graph Builder)
- VG toolkit
- Cactus

**Alignment-based:**
- Short reads: GATK HaplotypeCaller, FreeBayes, BCFtools
- Long reads: Clair3

## Integration with Rhizomorph

The Mycelia workflow integrates with the Rhizomorph graph-based assembly system to:
- Build k-mer graphs from sequencing reads
- Perform probabilistic assembly
- Validate assemblies through graph traversal
- Compare traditional assembly methods with graph-based approaches

## Tool Status

See the [Tool Wrapper Status](../planning-docs/TOOL_WRAPPER_STATUS.md) document for current implementation status of each tool wrapper.
See the [Tool Wrapper Status](https://github.com/cjprybol/Mycelia/blob/main/planning-docs/TOOL_WRAPPER_STATUS.md) document for current implementation status of each tool wrapper.
