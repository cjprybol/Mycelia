# Generative Sequence Design: Open-Source Development Roadmap

## Phase 6 of the Mycelia Comprehensive Roadmap

---

## Vision

Extend Mycelia.jl into the premier open-source platform for AI-guided biological
sequence design, with a fully Julia-native weighted graph + generative AI
framework as a new paradigm for computational genome engineering.

---

## Relationship to Mycelia Roadmap

The Generative Sequence Design Framework is **Phase 6** of the Mycelia
Comprehensive Roadmap, building on all preceding phases:

| Phase | Focus                               | Status            | Relationship to Phase 6                          |
| ----- | ----------------------------------- | ----------------- | ------------------------------------------------ |
| **1** | Rhizomorph graph ecosystem          | In Progress       | **Prerequisite** — graph data structures         |
| **2** | Core algorithms & performance       | In Progress       | **Prerequisite** — path-finding, simplification  |
| **3** | Production readiness & validation   | Planned           | **Prerequisite** — stability, testing            |
| **4** | Documentation & user experience     | Planned           | Enables community adoption of GSD                |
| **5** | Advanced features (RL, TDA, hybrid) | Planned           | **Partial prerequisite** — RL traversal policies |
| **6** | **Generative Sequence Design**      | **This document** | —                                                |

---

## Phase 6 Sub-Phases

### Phase 6.1: Julia-Native ML Pipeline (Q2-Q3 2026)

**New module:** `src/ml-pipeline.jl`

The foundational ML component — a fully self-contained Julia pipeline from
multi-omic graph features to phenotype prediction. This is the Julia-native
equivalent of Python AutoGluon/scikit-learn pipelines, independently developed
as open-source.

#### Components

| Component           | Implementation                                                                                                                    | Dependencies                      |
| ------------------- | --------------------------------------------------------------------------------------------------------------------------------- | --------------------------------- |
| Feature extraction  | Leverage existing `kmer-analysis.jl`, `sentencepiece.jl`, `genome-features.jl`, `codon-optimization.jl`, `amino-acid-analysis.jl` | Existing Mycelia                  |
| Feature engineering | Normalization, selection, pair features (cosine distance, concatenation, element-wise products)                                   | New code                          |
| Model training      | MLJ.jl framework with EvoTrees.jl (gradient boosting), DecisionTree.jl (random forests), Flux.jl (neural networks)                | Julia packages                    |
| Evaluation          | Biological hold-out splits (cluster-representative K-fold via `clustering.jl`), feature importance, SHAP values (ShapML.jl)       | Julia packages + existing Mycelia |
| Prediction          | Inference on new sequences, uncertainty quantification                                                                            | New code                          |
| Active learning     | Identify high-value experiments, rank by model uncertainty                                                                        | New code                          |

#### Milestones

- [ ] Package evaluation: MLJ.jl, EvoTrees.jl, DecisionTree.jl, Flux.jl
      compatibility testing
- [ ] Feature extraction pipeline operational (using existing Mycelia modules)
- [ ] Biological hold-out splits implemented and validated
- [ ] Model training framework with MLJ.jl integration
- [ ] Feature importance / SHAP working end-to-end
- [ ] Benchmarked on publicly available phenotype datasets

### Phase 6.2: Phenotypic Weighting Module (Q3 2026)

**New module:** `src/phenotypic-weighting.jl`

Connects ML pipeline outputs to graph edge/node weights for the joint weighting
system.

#### Components

| Component             | Implementation                                                         |
| --------------------- | ---------------------------------------------------------------------- |
| Weight computation    | Transform ML predictions to graph-compatible weights                   |
| Joint weighting       | W_joint = α × W_genotypic + β × W_phenotypic                           |
| Weight propagation    | Propagate phenotypic weights across graph layers (DNA → RNA → Protein) |
| Hyperparameter tuning | Optimize α/β via cross-validation                                      |

#### Milestones

- [ ] Weight computation from ML model outputs
- [ ] Joint weighting formula implemented
- [ ] Cross-layer weight propagation
- [ ] α/β optimization framework

### Phase 6.3: AI Model Integration Layer (Q3 2026)

**New module:** `src/language-model-integration.jl`

Wrappers for foundation model inference. Uses PythonCall.jl for model inference
while keeping all pre/post-processing in Julia.

#### Components

| Model                | Wrapper                            | Julia Integration            |
| -------------------- | ---------------------------------- | ---------------------------- |
| Evo2 (Arc Institute) | DNA sequence generation            | PythonCall.jl → Julia arrays |
| ESM-2 (Meta AI)      | Protein embeddings                 | PythonCall.jl → Julia arrays |
| ProtTrans            | Protein embeddings (alternative)   | PythonCall.jl → Julia arrays |
| ESMFold              | Fast structure prediction          | PythonCall.jl → Julia arrays |
| AlphaFold2           | Gold-standard structure validation | PythonCall.jl → Julia arrays |

#### HPC Job Templates

- NERSC GPU job template for Evo2 inference
- NERSC GPU job template for ESM-2/ESMFold batch embeddings
- NERSC GPU job template for AlphaFold2 validation
- All via `src/slurm-templates.jl`

#### Milestones

- [ ] Evo2 wrapper with Julia input/output (validated — already working)
- [ ] ESM-2 wrapper with batch embedding extraction
- [ ] ESMFold wrapper with structure prediction
- [ ] AlphaFold2 wrapper for validation
- [ ] NERSC SLURM templates for all models
- [ ] In-silico saturation mutagenesis pipeline

### Phase 6.4: Generative Optimization Pipeline (Q3-Q4 2026)

Extends Rhizomorph path-finding with phenotypic constraints for genome
generation.

#### Components

| Component                 | Implementation                                                     | Extends                                |
| ------------------------- | ------------------------------------------------------------------ | -------------------------------------- |
| Distance transformation   | d(e) = 1/p(e) or -log(p(e))                                        | `src/distance-metrics.jl`              |
| Constrained A\* traversal | Essential genes as required waypoints                              | New code, uses Rhizomorph graphs       |
| Beam search               | Top-k partial paths with biological constraints                    | New code                               |
| MCMC sampling             | Probabilistic walks with phenotypic weighting                      | Extends Rhizomorph probabilistic walks |
| Viterbi decoding          | Maximum likelihood path                                            | `src/viterbi-next.jl`                  |
| Consensus assembly        | Multiple traversals → consensus                                    | `src/assembly.jl`                      |
| Candidate ranking         | Joint probability + functional predictions + synthesis feasibility | New code                               |

#### Milestones

- [ ] Distance transformation integrated with joint weights
- [ ] Constrained path-finding with essential gene requirements
- [ ] Beam search with biological constraint filtering
- [ ] MCMC sampling with phenotypic edge weights
- [ ] Consensus assembly of generated walks
- [ ] Candidate ranking and output formatting (FASTA + GFF3 + predictions)

### Phase 6.5: Validation Framework & Benchmarks (Q4 2026)

#### Components

- Benchmark datasets for generative design evaluation
- Synthetic genome evaluation metrics (essential gene coverage, reading frame
  integrity, GC content, codon usage)
- Community validation protocols
- Comparison to existing tools

#### Milestones

- [ ] Benchmark dataset selection (public genomes with known phenotypes)
- [ ] Evaluation metric suite implemented
- [ ] End-to-end pipeline benchmarked
- [ ] Results documented for publication

### Phase 6.6: Application Templates (Q4 2026+)

Domain-specific configurations using the general framework:

| Template                       | Target Organisms           | Phenotype Model                                 | Status              |
| ------------------------------ | -------------------------- | ----------------------------------------------- | ------------------- |
| **DOE: Designer microbes**     | _P. syringae_, _P. putida_ | Biosynthetic pathway yield, rhizosphere fitness | Primary (LDRD 2026) |
| **DOE: Biosynthetic pathways** | Various                    | Compound production, pathway efficiency         | Primary (LDRD 2026) |
| **Phage therapy**              | ESKAPE pathogens           | Host-range, killing efficacy                    | Community template  |
| **Microbiome engineering**     | Various                    | Community composition, metabolic output         | Future              |

Each template includes:

- Domain-specific phenotype model configuration
- Recommended feature engineering preset
- Validation protocol
- Example notebook (Jupytext .jl format)

---

## Publication Plan

### Paper 1: Generative Sequence Design Methods Paper (Q3-Q4 2026)

- **Target venue:** Nature Methods, Genome Research, or Bioinformatics
- **Content:** Framework description + one application domain validation
- **Key figures:**
  1. Framework architecture diagram (multi-layer graphs + ML + optimization)
  2. Julia ML pipeline performance vs. Python baselines
  3. Feature importance leaderboard across feature categories
  4. Generated genome evaluation metrics
  5. Comparison to existing tools/approaches
- **Depends on:** Phase 6.1-6.5 completion

### Paper 2: Application-Specific Papers (2027+)

- Domain-specific validation papers for each application template
- Co-authored with domain experts

---

## Community Engagement Plan

### Open-Source Infrastructure

- MIT license (consistent with Mycelia)
- GitHub Discussions for design proposals
- Issue templates for bug reports, feature requests, new application templates
- Contributing guide focused on adding new application domains

### Tutorials and Documentation

- Quickstart tutorial: "Design your first genome in 30 minutes"
- Application-specific tutorials (one per template)
- API documentation (auto-generated via Documenter.jl)
- Jupytext notebook examples

### Outreach

- Workshop at relevant conference (ISCB, ASM, Phage Australia, APS)
- Blog post / preprint announcing framework
- Julia community presentation (JuliaCon or similar)

---

## Infrastructure Requirements

| System                 | Purpose                                              | Allocation                       |
| ---------------------- | ---------------------------------------------------- | -------------------------------- |
| **NERSC (Perlmutter)** | GPU inference (Evo2, ESM-2, AlphaFold2), ML training | DOE allocation                   |
| **Lawrencium**         | Development, batch processing                        | LBNL allocation                  |
| **Lovelace**           | CI/CD, Memgraph, large-memory jobs                   | Dedicated node                   |
| **HPC Storage**        | Data persistence                                     | CFS (community), HPSS (archival) |

### Containerization

- Apptainer images for AI model inference on HPC
- Julia depot snapshot for reproducible environments
- CI/CD via Lovelace (GitLab CI pipeline already configured)

---

## Related Beads

| Bead    | Title                                                 | Priority | Relationship                        |
| ------- | ----------------------------------------------------- | -------- | ----------------------------------- |
| td-l4xv | Establish generative design canonical docs in Mycelia | P1       | Parent epic                         |
| td-sozu | Generative genome design model                        | P4       | Research direction                  |
| td-hhr0 | Bacteriophage-specific assembly methods               | P4       | Application template input          |
| td-yb6v | Protein language model conservation scoring           | P3       | Feeds Phase 6.3                     |
| td-zy5t | Multi-omic knowledge graph infrastructure             | P3       | Supports multi-layer graph approach |
| td-y1g7 | Graph embedding stepwise optimization                 | P4       | Phase 6.4 research                  |

---

_Document Version: 1.0_ _Last Updated: February 2026_ _Author: Cameron Prybol,
Lawrence Berkeley National Laboratory_
