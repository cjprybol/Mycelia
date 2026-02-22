# Mycelia Comprehensive Development Roadmap v1.0

## 1. Executive Summary & Core Mission

**Source:** `ASSEMBLY_ROADMAP_UNIFIED.md`

This roadmap consolidates all previous planning documents into a single source
of truth. Our mission is to evolve Mycelia from a research framework into a
production-ready, scientifically rigorous, and user-accessible platform for
sequence graph analysis and assembly. We will achieve this by first implementing
the robust and type-safe **Rhizomorph** graph ecosystem, followed by building
core assembly algorithms, ensuring production stability, and finally, enhancing
user experience and deploying advanced research features.

### Core Mission

Transform Mycelia from a research-oriented assembly framework into an
accessible, production-ready platform that serves three user constituencies:

1.  **Assembly Beginners**: Clear entry points, simple examples, guided
    workflows
2.  **Domain Experts**: Advanced features, performance optimization, comparative
    tools
3.  **Developers**: Extensible architecture, clear APIs, contribution pathways

### Guiding Principles

- **Accuracy Over Everything:** Correctness over speed. Preserve biological
  reality, including quality scores and strand orientation, at every step.
- **Type-Safe by Design:** Eliminate string conversions in biological pipelines
  to ensure correctness and performance.
- **Progressive Disclosure:** Build a system that is simple for beginners to use
  for standard tasks but provides deep, powerful features for experts.
- **Test-Driven Development:** Implement the comprehensive testing framework in
  parallel with the core ecosystem to guarantee correctness from the ground up.
- **Biological Correctness:** Strand-aware edge transitions, canonical k-mer
  awareness, and proper handling of different sequence types.
- **Modular Architecture:** Reusable components, clear conversion pathways, and
  integration with existing bioinformatics tools.
- **Research-First Approach:** Prioritize novel algorithms and experimental
  strategies alongside production features.

## 2. Key Scientific Advances & Innovations

**Source:** `ROADMAP.md`, `ASSEMBLY_ROADMAP.md`

- **Novel 6-Graph Type Hierarchy**: A comprehensive, type-safe graph system that
  maintains biological correctness while enabling efficient computation.
- **Quality-Aware Assembly**: First assembly framework to preserve per-base
  quality scores throughout the entire assembly process, enabling higher
  accuracy.
- **Intelligent Self-Optimizing Assembly**: A learnable, self-optimizing
  assembler that uses reinforcement learning to eliminate manual parameter
  tuning while maximizing assembly accuracy.
- **Topology-Aware Assembly Optimization (TDA)**: Use Betti curves and
  (optionally) persistent homology across coverage/quality/confidence
  filtrations to guide graph cleaning and parameter selection (see
  `planning-docs/TDA_INTEGRATION_PLAN.md`).
- **Iterative Maximum Likelihood Assembly**: A complementary statistical path
  resampling approach for read-level error correction.
- **Zero String Conversions**: Complete elimination of string conversions in
  bioinformatics pipelines for type safety and efficiency.
- **Cross-validation Pipeline**: For robust assembly confidence and accuracy
  assessment.
- **Hierarchical Reinforcement Learning**: For automated, intelligent parameter
  tuning during assembly.
- **Type-Safe BioSequence Architecture**: All graphs use proper `Kmers.*Kmer{K}`
  or `BioSequences.Long*` types, eliminating ambiguity.

## 3. Phased Development Plan

**Source:** `ASSEMBLY_ROADMAP_UNIFIED.md` (for structure), with details from all
documents.

### **Phase 1: Foundational Rhizomorph Ecosystem (CRITICAL PRIORITY)**

_Goal: Implement the core type-safe, biologically accurate graph data
structures._

_Status (current): Core rhizomorph module, fixed-/variable-length builders,
strand conversions (DNA/RNA double/canonical), graph-type conversions
(fixed↔variable, quality↔non-quality), simplification helpers (tip
removal/linear-chain collapse), and path-finding/GFA I/O are implemented. JLD2
round-trip coverage lives in `end_to_end_graph_tests`. Migration of legacy graph
APIs/tests to Rhizomorph is in progress with direct ports (no shims); latest
ports restore coverage for `graph_algorithms_next`, `end_to_end_graph_tests`
(evidence positions, qualmer edge/vertex quality, doublestrand/canonical
conversion, GFA vertex counts), singlestrand k-mer suites, and `gfa_io_next`.
Remaining gaps: broaden variable-length traversal tests and edge-removal
simplification cases; remove legacy includes once coverage lands._

- **1.1: Establish Module Architecture:**
  - Create the `rhizomorph` module structure as defined in the ecosystem plan.
  - Consolidate all enums, vertex/edge data structures, and shared utilities
    into this core module.
- **1.2: Implement Fixed-Length Graphs (TDD Approach):**
  - Implement **N-gram Graphs** for Unicode strings.
  - Implement **K-mer Graphs** for FASTA data (`DNA`, `RNA`, `AA`) using the
    "SingleStrand-First" approach.
  - Implement **Qualmer Graphs** for FASTQ data, ensuring PHRED scores are
    handled correctly as integers (`UInt8`).
- **1.3: Implement Variable-Length Graphs (TDD Approach):**
  - Implement **String Graphs** for variable-length Unicode strings.
  - Implement **FASTA Graphs** using an overlap-layout-consensus approach.
  - Implement **FASTQ Graphs**, ensuring proper quality score aggregation and
    provenance tracking during simplification.
- **1.4: Implement Comprehensive Testing Framework:**
  - Develop and run the full test suite for all 24 graph combinations against
    the standard "pathological" graph cases (linear, bubble, repeat, etc.).
  - Achieve 100% test passage for all graph construction and I/O round-trip
    tests.

### **Phase 2: Core Algorithms & Performance (HIGH PRIORITY)**

_Goal: Build the essential graph manipulation and assembly algorithms on the
Rhizomorph foundation._

- **2.1: Pathfinding & Sequence Reconstruction:**
  - Implement a type-stable `find_eulerian_paths_next()` and
    `path_to_sequence()` for all graph types.
  - Implement strand-aware path traversal logic.
- **2.2: Graph Simplification & Canonicalization:**
  - Implement canonicalization as a post-processing step to collapse
    SingleStrand graphs into a DoubleStrand representation.
  - Implement essential simplification algorithms: tip removal, bubble popping
    (quality-score aware), and basic repeat resolution.
- **2.3: Performance Optimization & Benchmarking:**
  - Profile memory usage and speed for graph construction across all types.
  - Implement multi-threaded graph construction and parallel processing for key
    algorithms.
  - Investigate and implement memory-efficient streaming or on-disk graph
    representations for large datasets.

### **Phase 3: Production Readiness & Validation (HIGH PRIORITY)**

_Goal: Make Mycelia a robust, reliable, and verifiable tool for scientific
research._

- **3.1: Finalize Technical Debt & Migration:**
  - Complete the full migration from `MetaGraphs` to `MetaGraphsNext`.
  - Ensure backward compatibility or provide clear migration guides for I/O
    formats (e.g., GFA).
- **3.2: Create a Validation & CI Suite:**
  - **Test Coverage & Quality**:
    - Fix all remaining `Pkg.test()` failures to achieve 100% test passage.
    - Achieve >90% code coverage (up from 80% goal).
    - Add integration tests for common user workflows.
    - Create performance regression tests to prevent slowdowns.
    - Add fuzzing for robust input validation.
    - Test on diverse biological datasets (viral, bacterial, eukaryotic).
  - **Automation & CI**:
    - Automate regression testing and benchmarking via GitHub Actions.
  - **Validation Datasets**:
    - Curate standard validation datasets (e.g., simulated reads, known
      microbial genomes) with ground-truth assemblies.
    - Create simulated data generators for controlled testing.
    - Benchmark against gold-standard assemblies from other tools.
    - Document expected outcomes for all validation datasets.
- **3.3: Integrate External Validation Tools:**
  - Implement wrappers for industry-standard validation tools: `run_quast()`,
    `run_busco()`, and `run_mummer()` to allow for seamless assembly quality
    assessment.

### **Phase 4: User Accessibility & Documentation (MEDIUM PRIORITY)**

_Goal: Transform Mycelia into a product that is easy to install, learn, and
use._

- **4.1: Develop the "Probabilistic Assembly Hub":**
  - Create a central documentation landing page explaining the core concepts.
  - Write clear "Getting Started" tutorials for common use cases (e.g.,
    microbial genome assembly, text analysis).
  - Provide a Quick Start Kit with pre-configured environments and example
    datasets.
- **4.2: Package Polish & User Experience:**
  - Standardize function naming conventions and provide informative error
    messages.
  - Add progress indicators for long-running operations.
  - Finalize installation guides for different operating systems and HPC/cloud
    environments.
- **4.3: Implement Comprehensive Documentation & Tutorials (Detailed Tasks):**
  - See Section 5 for a full breakdown.

### **Phase 5: Advanced Features & Research (MEDIUM PRIORITY)**

_Goal: Implement the novel scientific breakthroughs that differentiate Mycelia._

- **5.1: Advanced Assembly Strategies:**
  - Implement **Unified Assembly Interface**.
  - Complete the **Hybrid Assembly** methods (e.g., OLC long-reads + Qualmer
    short-reads).
  - Implement and benchmark the **Intelligent Self-Optimizing Assembly**
    strategies (Iterative Maximum Likelihood, Prime k-mer progression).
- **5.2: Visualization & Monitoring:**
  - Develop a real-time assembly dashboard to monitor progress, quality metrics,
    and resource usage.
  - Create visualization tools for graph topology, path selection, and quality
    score heatmaps.
- **5.3: Reinforcement Learning Framework:**
  - Implement and train the hierarchical RL framework for self-optimizing
    parameter selection.
  - Prepare manuscripts documenting the novel quality-aware assembly and
    RL-guided methods.

## 4. Technical Specifications & Implementation Details

### 4.1. Graph Type Hierarchy Specification

**Source:** `ASSEMBLY_ROADMAP.md`

#### Fixed-Length Graphs (Assembly Foundation)

- **N-gram Graphs**: Unicode text analysis foundation based on ngram vertices
  connected by n-1 edge overlaps observed in the dataset
- **K-mer Graphs**: Type-stable `Kmers.DNAKmer{K}` vertices with NO string
  conversions and edges based on K-1 overlaps between kmers observed in the
  dataset
- **Qualmer Graphs**: Quality-aware k-mers preserving PHRED scores and edges
  based on K-1 overlaps between qualmers observed in the dataset

#### Variable-Length Graphs (Simplified Products)

- **String Graphs**: Variable-length unicode strings from N-gram simplification
  or direct creation from Strings via overlap-layout-consensus
- **FASTA Graphs**: `BioSequences.LongDNA/RNA/AA` from k-mer simplification or
  direct creation from FASTA records via overlap-layout-consensus
- **FASTQ Graphs**: Quality-preserving BioSequences with per-base scores from
  qualmer graph simplification or direct creation from FASTQ records via
  overlap-layout-consensus

#### Qualmer Implementation Requirements

- **Core Type**:
  ```julia
  struct Qualmer{K}
      sequence::BioSequences.DNAKmer{K}
      quality_scores::Vector{Int8}  # PHRED scores
      observations::Set{QualmerObservation}
      joint_probability::Float64    # Joint confidence in k-mer existence
  end
  ```

#### Joint Probability Calculation

When a k-mer is observed multiple times with different quality scores:

- Convert PHRED scores to probabilities: `p = 10^(-phred/10)`
- Calculate joint probability of k-mer correctness across all observations
- Use log-space arithmetic for numerical stability
- Example: 3 observations at 99% confidence each → joint confidence > 99%

#### Probabilistic Algorithms

- **Probabilistic walks** with configurable transition probabilities
- **Shortest path algorithms** where distance ∝ (1 - probability)
- **Maximum likelihood/weight walks** for high-confidence path finding
- **Viterbi maximum likelihood** path inference

### 4.2. Unified Assembly Interface

**Source:** `ASSEMBLY_ROADMAP.md`

A high-level API for common assembly tasks.

```julia
function Mycelia.Rhizomorph.assemble_genome(reads; method=:qualmer_graph, k=31, error_rate=0.01)
function Mycelia.Rhizomorph.polish_assembly(assembly, reads; iterations=3)
function Mycelia.Rhizomorph.validate_assembly(assembly, reference=nothing)
```

### 4.3. Intelligent Self-Optimizing Assembly System 🔵

**Source:** `ASSEMBLY_ROADMAP.md`

#### Strategic Vision:

**Create a learnable, self-optimizing assembler** that uses reinforcement
learning to eliminate manual parameter tuning while maximizing assembly
accuracy. The system will learn optimal assembly strategies through
cross-validation and simulated training, making intelligent decisions about
k-mer selection, error correction, and termination conditions.

#### Core Philosophy:

- **Accuracy-First**: Prioritize assembly accuracy over contiguity or speed
- **Dynamic k-mer Selection**: Iteratively process prime k-mer sizes until
  corrections stop
- **Cross-Validation**: Use 5-10 fold validation to create trusted consensus
  assemblies
- **Learnable Parameters**: Train on diverse simulated datasets to generalize to
  real data
- **Sparsity-Based Optimization**: Use sparsity detection to find optimal
  starting k-mer sizes

#### Architecture: Hierarchical Reinforcement Learning with Multiple Approaches

**1. High-Level Policy (Meta-Controller)**:

- **Algorithm**: Deep Q-Network (DQN) with experience replay
- **State**: Overall assembly quality, current k-mer size, memory usage,
  correction rate
- **Actions**: Continue with current k, move to next prime k, or terminate
- **Termination**: Based on correction rate, memory limits (32GB), or max k
  (~101)

**2. Low-Level Policy (Error Correction Controller)**:

- **Algorithm**: Policy Gradient (PPO) for continuous parameter spaces
- **State**: Local graph topology, quality scores, coverage patterns
- **Actions**: Viterbi parameters, path selection strategies, confidence
  thresholds
- **Integration**: Uses existing Viterbi + probabilistic path algorithms

**3. Monte Carlo Tree Search (Game-Based Assembly)** 🆕:

- **Algorithm**: UCB1-based tree search with simulated rollouts
- **Game Formulation**: Assembly as sequential decision game with known target
  (training)
- **State**: Current partial assembly, available reads, quality metrics
- **Actions**: Read selection, overlap decisions, path extensions
- **Reward**: Accuracy vs. known reference (weighted 1000x) + efficiency metrics
- **Rollout Policy**: Fast heuristic assembly for leaf node evaluation
- **Selection**: UCB1 balancing exploration vs exploitation

#### Implementation Phases:

**Foundation (High Priority)**

- a: Iterative prime k-mer progression algorithm
- b: Accuracy-prioritized reward function
- c: Cross-validation pipeline for quality assessment

**Learning System (Medium Priority)**

- a: RL environment for assembly decision making
- b: Simulation framework for training on diverse datasets
- c: Policy networks for parameter optimization

**Visualization & Automation (Medium Priority)**

- a: Decision pathway visualization showing confidence levels
- b: Automated parameter selection based on learned policies
- c: Real-time quality assessment during assembly
- d: Interactive tools for exploring assembly decisions

#### Key Algorithms:

**Dynamic k-mer Selection**:

```julia
function mycelia_assemble(reads; max_k=101, memory_limit=32_000_000_000)
    k = find_initial_k(reads)  # First prime k with sparsity
    while k <= max_k
        graph = build_qualmer_graph(reads, k=k)
        corrections_made = correct_errors_at_k(graph, k)
        if !should_continue_k(graph, corrections_made, k)
            k = next_prime_k(k)
        end
    end
    return finalize_assembly(graphs)
end
```

### 4.4. Iterative Maximum Likelihood Assembly Module

**Source:** `ASSEMBLY_ROADMAP.md`

#### Strategic Vision: Statistical Path Resampling Assembly

**Create an iterative assembly module** that complements the intelligent
assembly system by using statistical path resampling and maximum likelihood
approaches for read-level error correction.

#### Core Philosophy:

- **Sparsity-Based Initialization**: Use same sparsity assessment as intelligent
  assembly for starting k
- **Read-Level Correction**: Consider each read individually for statistical
  path improvement
- **Maximum Likelihood Approach**: Replace reads with higher likelihood paths
  proportionally
- **Iterative Prime Progression**: Like intelligent assembly, iterate through
  prime k-mer sizes
- **Viterbi Integration**: Leverage existing viterbi-next.jl algorithms for path
  finding

#### Algorithm Flow:

1. **Initialize**: Find starting prime k using sparsity assessment.
2. **For each prime k**:
   - Build qualmer graph from current read set.
   - For each read, find its path likelihood, search for higher-likelihood
     alternatives using Viterbi, and replace the read proportionally to the
     likelihood improvement.
   - Write the entire updated read set to a new, timestamped FASTQ file.
   - Continue until corrections saturate, then move to the next prime k.
3. **Output**: A series of timestamped FASTQ files showing read evolution.

### 4.5. Cross-Validation Pipeline for Quality Assessment

- K-fold cross-validation framework (configurable 3-10 folds)
- Holdout validation through read mapping and coverage analysis
- Consensus pangenome generation with statistical comparison

### 4.6. Performance Validation and Optimization

Validate theoretical improvements with real benchmarks:

```julia
# Performance benchmarking suite
function benchmark_graph_construction(legacy_vs_next)
function benchmark_memory_usage(canonical_vs_stranded)
function benchmark_algorithm_performance(phase2_vs_legacy)
function benchmark_quality_aware_assembly(qualmer_vs_kmer)
```

- **Optimization**: Profile and optimize hot paths in core algorithms.
- **Scalability**:
  - [ ] Support for distributed computing frameworks.
  - [ ] Cloud storage integration (S3, GCS).
  - [ ] Integration with batch processing and workflow management systems.
  - [ ] Implement real-time resource usage monitoring.

### 4.7. Strategic Tool Integration & Functionality Gaps

**Source:** `ROADMAP.md`

#### Completed Integrations (2025-12-10)

- [x] **SentencePiece** - Subword tokenization for ML/NLP applications on
      biological sequences
  - Implementation: `src/sentencepiece.jl`
  - Tests: `test/8_tool_integration/sentencepiece.jl`
  - Features: Train models, encode/decode DNA/RNA/AA/text, subword
    regularization

#### Gap Analysis and Tool Integration

- [ ] **Comprehensive functionality audit**: Identify missing functions vs
      available tools.
- [ ] **Tool integration opportunities**: Find external tools for remaining
      gaps.
- [ ] **Wrapper implementation priority**: Focus on high-value integrations.
- [ ] **Custom implementation assessment**: Determine what truly needs native
      implementation.

#### Functionality Gap Assessment

- **Assembly Quality Assessment:**
  - [ ] `detect_assembly_errors()`: For misassembly detection.
- **Quality Control & Visualization:**
  - [ ] `estimate_genome_size_from_kmers()`: For k-mer based size estimation.
  - [ ] `detect_contamination()`: For contamination screening.
- **Sequence Analysis:**
  - [ ] `calculate_gc_content()`: Can be derived from k=1 k-mer counting.
  - [ ] `construct_phylogeny()`: To be handled via FastTree/IQ-TREE integration.
  - [ ] `calculate_evolutionary_distances()`: To be handled via external tools.

#### Assembly Validation Tools

- [x] `QUAST` integration - `run_quast()` (implemented; opt-in extended tests
      and default outdir derivation)
- [x] `BUSCO` integration - `run_busco()` (implemented; auto-lineage default,
      dataset predownload helper, opt-in extended tests)
- [ ] `MUMmer` integration - `run_mummer()`
- [x] `CheckM2` integration - `run_checkm2()`
- [ ] `Mauve` integration - `run_mauve()`

#### Advanced Quality Assessment Tools

- [ ] **ALE**: `run_ale()` for reference-free quality assessment.
- [ ] **FRCbam**: `run_frcbam()` for feature response curves.
- [ ] **4CAC**: `run_4cac()` for ML-based contig classification.

#### Strain-Aware Assembly Methods

- [ ] **HyLight**: `run_hylight()` for hybrid strain-resolved assembly.
- [ ] **STRONG**: `run_strong()` for strain resolution on assembly graphs.
- [ ] **Strainy**: `run_strainy()` for strain phasing from long reads.

#### Assembly Polishing and Merging Tools

- [ ] **Apollo**: `run_apollo()` for HMM-based polishing.
- [ ] **QuickMerge**: `run_quickmerge()` for assembly merging.
- [ ] **Homopolish**: `run_homopolish()` for reference-based correction.

#### Quality Control, Phylogenetics, and Visualization

- [ ] **QC**: `MultiQC`, `NanoPlot`, `PycoQC`.
- [ ] **Phylogenetics**: `FastTree`, `IQ-TREE`, `RAxML`.
- [ ] **Visualization**: `Circos`, IGV integration, interactive assembly graphs.

### 4.8. Multi-Scale K-mer Analysis Framework

**Source:** `ROADMAP.md` (Universal polymer assembler algorithm design)

- **Simultaneous Multi-K Analysis**:
  - [ ] Implement prime k-mer simultaneous analysis (e.g., 3,5,7,11,13,17,19).
  - [ ] Add sliding window instantaneous quality averaging.
  - [ ] Support adaptive k-mer selection based on coverage patterns.
- **Quality-Aware Analysis**:
  - [ ] Incorporate quality scoring into k-mer window analysis.
  - [ ] Implement adaptive k-mer selection algorithms.
  - [ ] Add consensus k-mer identification with confidence scoring.
- **Bootstrap Validation Framework**:
  - [ ] Add `bootstrap_assembly_validation()` function.
  - [ ] Implement robust parameter optimization with confidence intervals.
  - [ ] Support statistical validation of assembly parameters.

### 4.9. Probabilistic Assembly Framework

**Source:** `ROADMAP.md` (Perfect metagenomic workflow design philosophy)

- **Maximum Likelihood Assembly**:
  - [ ] Implement probabilistic path selection through assembly graphs.
  - [ ] Add consensus sequence generation with path probability weighting.
  - [ ] Support multiple path analysis with confidence scoring.
- **Confidence Interval Framework**:
  - [ ] Implement per-base confidence scoring based on path probabilities.
  - [ ] Add uncertainty quantification for ambiguous assembly regions.
  - [ ] Support probabilistic variant calling from graph ambiguities.
- **Statistical Validation**:
  - [ ] Add bootstrap validation for parameter optimization.
  - [ ] Implement confidence interval calculation for assembly metrics.
  - [ ] Support robust statistical assessment of assembly quality.

## 5. Documentation, Tutorials, & Package Polish

### 5.1. Probabilistic Assembly Hub & Tutorials

**Source:** `ASSEMBLY_ROADMAP_UNIFIED.md`

- **1. Create "Probabilistic Assembly Hub" Landing Page**
  - Explain what probabilistic assembly is, why to use it, and how it works
    (with diagrams).
  - Include a quick start guide and a decision tree for method selection.
- **2. Develop "Assembly in 5 Minutes" Tutorial**
  - Provide a simple, complete, copy-pasteable example:
    ```julia
    # Step 1: Load your data
    reads = Mycelia.load_fastq("example_reads.fastq")
    # Step 2: Run probabilistic assembly
    assembly = Mycelia.assemble_probabilistic(reads)
    # Step 3: Evaluate results
    metrics = Mycelia.evaluate_assembly(assembly)
    Mycelia.plot_assembly_quality(metrics)
    # Step 4: Save results
    Mycelia.write_assembly(assembly, "my_assembly.fasta")
    ```
- **3. Create Visual Decision Flow Chart**
  - Interactive flowchart guiding users from data type to recommended approach.

### 5.2. Comprehensive Documentation Plan

**Source:** `ROADMAP.md`, `ASSEMBLY_ROADMAP_UNIFIED.md`

- **User Guides & Tutorials (Beginner Track)**:
  - Step-by-step tutorials with no assumed knowledge.
  - "Getting Started" guide with working, copy-pasteable examples.
  - Installation and troubleshooting guide.
  - FAQ based on common user feedback.
- **Workflows & Cookbooks (Expert Track)**:
  - Create "Assembly Recipe Book" for common scenarios (e.g., "I have Illumina
    reads from E. coli...", "I need to assemble a metagenome...").
  - Performance optimization guide.
  - Develop interactive tutorials using Jupyter/Pluto notebooks with live
    examples.
- **Developer Documentation (Developer Track)**:
  - Comprehensive architecture overview.
  - API reference with function index and parameter guide.
  - Contribution guidelines (`CONTRIBUTING.md`), code style, and testing
    guidelines.
  - Release process documentation.
- **Scientific Documentation**:
  - Detailed algorithm descriptions with citations.
  - Benchmarking methodology and comprehensive comparison with other tools.
  - Scientific use case examples and tutorial notebooks.

### 5.3. Package Polish & User Experience

**Source:** `ROADMAP.md`

- **Function Discovery & Usability**:
  - [ ] Create function categories in the module.
  - [ ] Add `?Mycelia` help documentation.
  - [ ] Generate a function index automatically.
  - [ ] Standardize function naming conventions.
  - [ ] Provide informative error messages and progress indicators.
  - [ ] Implement smart default parameters for common use cases.
  - [ ] Add automatic output directory creation.
- **Robustness & Reliability**:
  - [ ] Add pre-flight checks for memory requirements and data quality.
  - [ ] Implement a checkpoint/resume system for recovery from interrupted
        operations.
  - [ ] Comprehensive input validation and logging options.
- **Fix Documentation-Reality Mismatch**:
  - [ ] Audit all documented functions.
  - [ ] Remove or mark non-existent functions.
  - [ ] Update all examples to use working functions.
  - [ ] Add status badges (stable/experimental/planned).

## 6. Success Metrics & Risk Mitigation

**Source:** `ASSEMBLY_ROADMAP_UNIFIED.md`, `ROADMAP.md`, `ASSEMBLY_ROADMAP.md`

### Success Metrics

- **User Accessibility**:
  - Time to first successful assembly: <10 minutes.
  - Documentation satisfaction: >80%.
- **Technical Excellence**:
  - Test coverage: >90%.
  - All tests passing on LTS Julia versions (e.g., 1.10+).
  - Import time: <10 seconds.
  - Performance vs competitors: Top 3 in speed and accuracy.
  - Memory efficiency: 50% reduction compared to previous versions, with
    optimized parallel processing.
  - Assembly Accuracy: >99% accuracy on simulated datasets.
- **Adoption & Impact**:
  - Growth in GitHub stars, active installations, citations, and contributors.
  - Generalization: Consistent performance across diverse organism types.
  - Automation: Minimal manual parameter tuning required.

### Risk Mitigation

- **Technical Risks**: Mitigate complexity with progressive disclosure; address
  performance issues with profiling.
- **User Adoption Risks**: Address steep learning curves with simplified
  tutorials; build trust with comprehensive benchmarks.
- **Maintenance Risks**: Combat technical debt with refactoring sprints; prevent
  documentation drift with automated doc testing.

### Known Operational Issues (carryover from `ROADMAP.md`)

- **Network/API Reliability**:
  - Completed: `ncbi_genome_download_accession()` in
    `src/reference-databases.jl` now uses the shared `with_retry()` helper.
  - TODO: Wrap remaining network calls with `with_retry()` for consistent
    retries/logging:
    - `download_genome_by_accession()` (~line 1143) and `load_ncbi_metadata()`
      (~lines 1950-2000) in `src/reference-databases.jl`.
    - Taxonomy functions using `ncbi-datasets-cli` (~line 1298) in
      `src/taxonomy-and-trees.jl`.
    - Any other functions invoking `Downloads.download()` or external CLIs that
      can transiently fail.
  - Reference utility: `src/utility-functions.jl` contains `with_retry()` to
    standardize backoff behavior.

## 7. Research Publication Preparation

**Source:** `ROADMAP.md`

### 7.1. Novel Algorithm Documentation & Manuscripts

- [ ] Prepare manuscript documenting the novel 6-graph type hierarchy and its
      performance benchmarks.
- [ ] Prepare manuscript on the quality-aware assembly methodology and its
      impact on accuracy.
- [ ] Prepare manuscript detailing the RL-guided self-optimizing assembly
      framework.

### 7.2. Benchmarking & Validation for Publication

- [ ] Conduct comprehensive benchmarking against standard, state-of-the-art
      assemblers.
- [ ] Perform detailed performance evaluation of all novel algorithms.
- [ ] Generate publication-quality figures, tables, and statistical analyses.

---

## 8. Phase 6: Generative Sequence Design (FUTURE)

_Goal: Extend Mycelia into a platform for AI-guided design of synthetic
biological sequences with specified functional properties._

This phase represents the culmination of the Rhizomorph graph ecosystem: using
the weighted sequence graph infrastructure built in Phases 1-5 as the foundation
for rational, generative genome engineering. Rather than assembling existing
sequences, Phase 6 enables _designing new sequences_ by combining multi-layer
weighted graphs with generative AI models and a Julia-native ML pipeline.

**Key components:**

- **6.1: Julia-Native ML Pipeline** — Feature extraction from graphs → model
  training → phenotype prediction → active learning (MLJ.jl, EvoTrees.jl,
  Flux.jl)
- **6.2: Phenotypic Weighting Module** — Connect ML predictions to graph
  edge/node weights
- **6.3: AI Model Integration Layer** — Evo2, ESM-2, AlphaFold2 wrappers via
  PythonCall.jl
- **6.4: Generative Optimization Pipeline** — Constrained path-finding with
  phenotypic weights
- **6.5: Validation Framework** — Benchmarks and evaluation metrics
- **6.6: Application Templates** — Domain-specific configurations (DOE microbes,
  phage therapy, biosynthetic pathways)

**See:** `planning-docs/GENERATIVE_SEQUENCE_DESIGN_PLAN.md` (technical spec),
`planning-docs/GENERATIVE_SEQUENCE_DESIGN_ROADMAP.md` (detailed Phase 6
roadmap), and `planning-docs/GENERATIVE_SEQUENCE_DESIGN_OVERVIEW.md`
(stakeholder summary).
