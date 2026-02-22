# Generative Sequence Design with Weighted Biological Graphs

## An Open-Source Framework for AI-Guided Biological Sequence Engineering

---

## Elevator Pitch (60 seconds)

Nature has optimized biological sequences over billions of years, but the space
of functional sequences is vastly larger than what evolution has explored. We
need methods to rationally design new sequences — whether for engineered living
materials, biosynthetic pathways, or precision antimicrobials.

We've developed the **Generative Sequence Design Framework**: a computational
platform built on Mycelia.jl that combines biological sequence graphs with
generative AI to rationally design synthetic genomes with specified functional
properties.

Our approach encodes biological diversity as weighted graphs across DNA, RNA,
and protein layers, where edge weights reflect both evolutionary abundance and
experimentally-derived phenotypic efficacy. A fully Julia-native ML pipeline
predicts phenotype from multi-omic graph features, and genome language models
(Evo2) explore novel sequence spaces beyond what nature has sampled, generating
candidate genomes optimized for target outcomes.

The result: computationally designed biological sequences with predicted
functional profiles, ready for synthesis and validation. We're transforming
sequence engineering from discovery-based serendipity to rational design.

---

## Executive Summary

### The Problem

Traditional biological sequence engineering relies on screening natural
diversity for desired properties — a slow, expensive process limited to the
sequences nature has already produced. Whether designing microbes for DOE
applications (rhizosphere engineering, biosynthetic pathways), engineering
phages for antimicrobial therapy, or building synthetic organisms for
biotechnology, the fundamental bottleneck is the same: we search when we should
be designing.

### Our Solution

The Generative Sequence Design Framework is a computational platform for
rational, AI-guided design of synthetic biological sequences with specified
functional properties. Rather than screening natural isolates for serendipitous
matches, the framework enables:

1. **Knowledge Integration**: Comprehensive encoding of genomic diversity as
   weighted sequence graphs (Mycelia Rhizomorph)
2. **Phenotype-Genotype Mapping**: Julia-native ML pipeline linking multi-omic
   graph features to functional phenotypes
3. **Generative Exploration**: AI-driven exploration of viable sequence spaces
   using state-of-the-art language models
4. **Probabilistic Optimization**: Weighted graph traversal to identify optimal
   genome configurations
5. **Predictive Validation**: In silico evaluation of candidate genomes before
   synthesis

### Key Innovation

Traditional approaches treat sequence design as a search problem through natural
diversity. This framework reframes it as a **generative design problem**: given
target phenotypes, construct optimal genomes by integrating biological
constraints with AI-powered sequence generation — all within a fully
self-contained, Julia-native open-source platform.

The framework operates across three interconnected representational layers:

| Layer       | Representation          | AI Integration                 | Mycelia Module             | Status        |
| ----------- | ----------------------- | ------------------------------ | -------------------------- | ------------- |
| **DNA**     | Pangenome graphs        | Genome Language Models (Evo2)  | `src/rhizomorph/`          | IN VALIDATION |
| **RNA**     | Transcriptome graphs    | Expression optimization        | `src/rhizomorph/`          | IN VALIDATION |
| **Protein** | Protein sequence graphs | Protein LMs + Folding models   | `src/rhizomorph/`          | IN VALIDATION |
| **ML**      | Feature → Phenotype     | MLJ.jl / EvoTrees.jl / Flux.jl | `src/ml-pipeline.jl` (NEW) | PLANNED       |

---

## Technical Overview

### Graph Architecture (Mycelia Rhizomorph)

The foundational data structure is a multi-layer weighted sequence graph
encoding:

**Construction Methods:**

- **k-mer / de Bruijn graphs**: Capture sequence composition and local variation
- **Overlap-Layout-Consensus (OLC) graphs**: Encode longer structural
  relationships
- **Tokenized graphs**: BPE and unigram tokenization adapted from NLP to capture
  functional motifs

**Weighting Schemes:**

| Weighting Type         | Description                                                                              | Status        |
| ---------------------- | ---------------------------------------------------------------------------------------- | ------------- |
| **Genotypic Weights**  | Abundance-based weights reflecting observed allele frequencies across genome collections | IN VALIDATION |
| **Phenotypic Weights** | Julia ML-derived weights linking multi-omic features to application-specific phenotypes  | PLANNED       |

### Julia-Native ML Pipeline

A fully self-contained ML workflow in Julia, from graph features to phenotype
prediction:

| Component         | Julia Package   | Purpose                                 |
| ----------------- | --------------- | --------------------------------------- |
| Framework         | MLJ.jl          | Unified ML training, evaluation, tuning |
| Gradient boosting | EvoTrees.jl     | High-performance ensemble models        |
| Random forests    | DecisionTree.jl | Interpretable classification/regression |
| Neural networks   | Flux.jl         | Custom deep learning architectures      |
| Interpretability  | ShapML.jl       | SHAP values for feature importance      |

Features extracted from existing Mycelia modules: k-mer frequencies, BPE/Unigram
tokens, gene annotations, codon usage metrics, protein domain annotations,
language model embeddings.

### AI Integration Points

**Genome Language Models (Evo2):**

- Generate novel DNA sequences conditioned on graph nodes
- Propose inter-gene connections and regulatory regions
- Enable exploration of unsampled sequence space

**Protein Language Models (ESM-2, ProtTrans):**

- Embed protein sequences in functional space
- Identify functionally equivalent variants
- Guide sequence optimization while preserving function

**Folding Models (ESMFold, AlphaFold):**

- Validate structural integrity of generated proteins
- Enable structure-based filtering of candidates
- Support in-silico saturation mutagenesis

### Optimization Strategy

- Transform edge weights to distances (d = 1/p) for shortest-path algorithms
- Apply A\*/Dijkstra/Beam/MCMC to find maximum-likelihood paths connecting
  essential features
- Biological constraints: essential gene coverage, reading frames, genome
  topology
- Consensus assembly of multiple traversals
- Ranked candidate genomes with confidence scores

---

## Application Domains

### Primary: DOE Applications (LDRD 2026 Aligned)

- **Engineered living materials** — designer microbes for materials science
- **Biosynthetic pathway engineering** — optimizing gene clusters for production
  of bioplastics, polymers, enzymes, energy-storing compounds
- **Rhizosphere engineering** — enhanced energy harvesting from plant-microbe
  interactions
- **Target organisms:** _Pseudomonas syringae_, _Pseudomonas putida_, and others

### Community Applications

- **Bacteriophage therapy** — rational design of synthetic phages targeting
  antibiotic-resistant pathogens (e.g., ESKAPE organisms including _E. coli_,
  _S. aureus_)
- **Microbiome engineering** — designing organisms for environmental or health
  applications
- **Vaccine design** — optimized antigen sequences

---

## Competitive Advantages

| Traditional Approach              | Generative Sequence Design              |
| --------------------------------- | --------------------------------------- |
| Screen natural isolates           | Rational design from first principles   |
| Limited to existing diversity     | Explores AI-generated sequence space    |
| Unpredictable functional outcomes | Designed for specific target phenotypes |
| Slow iteration cycles             | Rapid in silico optimization            |
| Discovery-based                   | Engineering-based                       |
| Python/mixed toolchains           | Fully Julia-native, self-contained      |

---

## Open-Source Ecosystem

- **Mycelia.jl** — Core graph infrastructure, algorithms, and ML pipeline
- **License:** MIT (open source)
- **Affiliation:** Lawrence Berkeley National Laboratory
- **Repository:**
  [github.com/cjprybol/Mycelia](https://github.com/cjprybol/Mycelia)
- **Julia ecosystem:** Fully integrated with Julia package registry
- **Community model:** Open contributions, application-domain templates

---

## Resource Requirements

### Computational Infrastructure

| Resource               | Purpose                                                          |
| ---------------------- | ---------------------------------------------------------------- |
| NERSC (Perlmutter)     | GPU inference (Evo2, ESM-2, AlphaFold2), large-scale ML training |
| Lawrencium             | Development, testing, batch processing                           |
| Lovelace               | CI/CD, continuous integration                                    |
| HPC storage (CFS/HPSS) | Data persistence and archival                                    |

### Expertise

- Computational genomics and bioinformatics
- Machine learning (Julia ecosystem)
- Molecular biology and genetics
- Synthetic biology and genome engineering

---

## Timeline

| Phase                     | Timeline     | Deliverables                                                  |
| ------------------------- | ------------ | ------------------------------------------------------------- |
| **Julia ML Pipeline**     | Q2-Q3 2026   | Feature extraction, model training, evaluation — all in Julia |
| **AI Integration**        | Q3 2026      | Evo2, ESM-2, AlphaFold2 wrappers; phenotypic weighting        |
| **Optimization**          | Q3-Q4 2026   | End-to-end genome generation, candidate ranking               |
| **Validation**            | Q4 2026      | Benchmarks, synthesis of top candidates                       |
| **Publication**           | Q4 2026-2027 | Methods paper (Nature Methods / Genome Research tier)         |
| **Application Templates** | 2027+        | DOE microbes, phage therapy, biosynthetic pathways            |

---

## Contact

**Cameron Prybol** | Lawrence Berkeley National Laboratory GitHub:
[cjprybol/Mycelia](https://github.com/cjprybol/Mycelia)

---

_Document Version: 1.0_ _Last Updated: February 2026_
