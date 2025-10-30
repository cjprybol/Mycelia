# Assembly Suite Overview

Mycelia provides a comprehensive assembly suite with multiple approaches for genome assembly. This document provides an overview of the assembly ecosystem and guidance on choosing the right approach for your data.

> **Important Distinction**: Mycelia provides two categories of assembly functionality:
>
> **Production-Ready**: Stable wrappers for 15+ established third-party assemblers (MEGAHIT, SPAdes, Flye, hifiasm, Canu, Unicycler, Velvet, MetaMDBG, and others). These are tested, documented, and ready for production use.
>
> **Experimental/Research**: Internal graph-based assembly algorithms, intelligent self-optimizing approaches, and reinforcement learning methods. These are research implementations that may require additional configuration and are subject to change.

## Third-Party Assembler Wrappers [Production-Ready]

Mycelia provides production-ready wrappers for established assembly tools:

### Short-Read Assemblers
- **MEGAHIT** (`run_megahit`) - Fast metagenome assembler
- **SPAdes** (`run_spades`) - Isolate genome assembler
- **metaSPAdes** (`run_metaspades`) - Metagenome assembler
- **SKESA** (`run_skesa`) - Strategic K-mer Extension for Scrupulous Assemblies
- **Velvet/MetaVelvet** (`run_velvet`, `run_metavelvet`) - Classic de Bruijn graph assemblers

### Long-Read Assemblers
- **Flye** (`run_flye`) - Long-read assembler for single genomes
- **metaFlye** (`run_metaflye`) - Long-read metagenome assembler
- **Canu** (`run_canu`) - Long-read assembler with error correction
- **hifiasm** (`run_hifiasm`) - Haplotype-resolved assembler for HiFi reads
- **MetaMDBG** (`run_metamdbg`) - Metagenomic de Bruijn graph assembler

### Hybrid Assemblers
- **Unicycler** (`run_unicycler`) - Hybrid short+long read assembler

### Polishing/Refinement Tools
- **Apollo** (`run_apollo`) - Assembly polishing
- **HyPo** (`run_hylight`) - Hybrid polishing
- **Strainy** (`run_strainy`) - Strain phasing and separation

## Internal Assembly Algorithms [Experimental]

The following algorithms are research implementations currently in experimental status:

## Internal Assembly Approaches [Experimental]

The following assembly approaches use Mycelia's internal graph-based algorithms. These are research implementations:

### 1. Standard K-mer Assembly [Experimental]

**Fixed-Length Graphs (Assembly Foundation)**
1. **N-gram Graphs** (`build_ngram_graph`) - Unicode character analysis
2. **K-mer Graphs** (`build_kmer_graph_next`) - Standard k-mer assembly
3. **Qualmer Graphs** (`build_qualmer_graph`) - Quality-aware k-mer assembly

**Variable-Length Graphs (Simplified Products)**
4. **String Graphs** (`string_to_ngram_graph`) - Variable unicode strings
5. **FASTA Graphs** (`build_biosequence_graph`) - Variable BioSequences
6. **FASTQ Graphs** (`build_quality_biosequence_graph`) - Quality-aware BioSequences

### Graph Type Hierarchy [Experimental]

Mycelia implements a 6-level graph hierarchy for assembly research:

For basic k-mer assembly without quality awareness:

```julia
import Mycelia

# Build k-mer graph
graph = Mycelia.build_kmer_graph_next(fasta_records; k=31, graph_mode=DoubleStrand)

# Convert to sequence graph
seq_graph = Mycelia.kmer_graph_to_biosequence_graph(graph, 31)

# Extract contigs
contigs = Mycelia.extract_contigs_from_graph(seq_graph)
```

**Status**: Experimental - Basic graph construction implemented, contig extraction in development

**Use Cases (when fully implemented)**:
- High-quality FASTA data
- Simple genomes without complex repeats
- Research into graph-based assembly

### 2. Quality-Aware Assembly [Experimental]

For assembly utilizing FASTQ quality scores:

#### Option A: Qualmer-Mediated (Recommended for most cases)
```julia
import Mycelia

# Build qualmer graph (k-mer level quality analysis)
qualmer_graph = Mycelia.build_qualmer_graph(fastq_records; k=31)

# Find quality-weighted paths
start_vertex = argmax(v -> qualmer_graph[v].joint_probability, vertices(qualmer_graph))
optimal_path = Mycelia.find_quality_weighted_path(qualmer_graph, start_vertex)

# Convert to quality-aware sequence graph
fastq_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph, 31)

# Extract final assembly
final_records = Mycelia.quality_biosequence_graph_to_fastq(fastq_graph, "assembly")
```

#### Option B: Direct Quality-Aware
```julia
import Mycelia

# Build quality-aware sequence graph directly
fastq_graph = Mycelia.build_quality_biosequence_graph(fastq_records)

# Extract contigs with quality preservation
contigs = Mycelia.quality_biosequence_graph_to_fastq(fastq_graph, "contigs")
```

**Status**: Experimental - Quality-aware graph construction implemented, assembly pipeline in development

**Use Cases (when fully implemented)**:
- Illumina short reads
- Error-prone long reads
- Research into quality-aware assembly

### 3. Intelligent Self-Optimizing Assembly [Experimental]

For automated parameter optimization:

```julia
import Mycelia

# Intelligent assembly with dynamic k-mer progression
results = Mycelia.mycelia_assemble(fastq_file; 
    max_k=101,
    memory_limit=32_000_000_000,
    output_dir="intelligent_assembly"
)

# Cross-validation for quality assessment
validation_results = Mycelia.mycelia_cross_validation(
    fastq_records;
    assembly_methods=[:intelligent, :iterative],
    k_folds=5
)
```

**Status**: Experimental - Core functionality implemented, requires further testing and validation

**Features**:
- Automatic prime k-mer progression
- Sparsity-based k-mer selection
- Memory monitoring and limits
- Cross-validation for quality assessment

**Use Cases (when fully validated)**:
- Research into parameter optimization
- Comparison with established assemblers
- Algorithm development

### 4. Iterative Statistical Assembly [Experimental]

For statistical path improvement:

```julia
import Mycelia

# Iterative assembly with read-level improvements
results = Mycelia.mycelia_iterative_assemble(fastq_file;
    max_k=101,
    memory_limit=32_000_000_000,
    output_dir="iterative_assembly"
)
```

**Status**: Experimental - Core functionality implemented, requires further testing and validation

**Features**:
- Read-level likelihood optimization
- Statistical path resampling
- Timestamped output for tracking evolution
- Viterbi integration (partially implemented)

**Use Cases (when fully validated)**:
- Research into statistical assembly approaches
- Algorithm development and comparison
- Detailed assembly process tracking

### 5. Reinforcement Learning Guided Assembly [Research/Proof-of-Concept]

For machine learning optimization (four experimental implementations available):

#### Custom RL Framework **[Experimental - Under Development]**
```julia
import Mycelia

# Train custom RL agent
agent = Mycelia.DQNPolicy()
env = Mycelia.AssemblyEnvironment(reads)
trained_agent = Mycelia.train_assembly_agent(env, agent; episodes=1000)

# Apply learned policy
results = Mycelia.apply_learned_policy(trained_agent, "genome.fastq")
```

#### ReinforcementLearning.jl Integration **[Experimental - Basic Implementation]**
```julia
import Mycelia

# Use well-tested RL algorithms
assembly, history = Mycelia.intelligent_assembly_rljl(
    reads;
    algorithm=:dqn,
    train_episodes=1000
)
```

#### POMDPs.jl Integration **[Experimental - Basic Implementation]**
```julia
import Mycelia

# Formal MDP/POMDP approach
assembly, history = Mycelia.intelligent_assembly_pomdp(
    reads;
    solver=:value_iteration,
    use_pomdp=false
)
```

#### Monte Carlo Tree Search (MCTS) **[Experimental - Game-Based Approach]**
```julia
import Mycelia

# Assembly as a sequential decision game
assembly, trajectory = Mycelia.intelligent_assembly_mcts(
    reads;
    n_simulations=1000,
    exploration_constant=sqrt(2)
)
```

#### Comparison Framework **[Experimental - Proof of Concept]**
```julia
import Mycelia

# Compare all four RL approaches
comparison = Mycelia.compare_rl_approaches(
    reads;
    approaches=[:custom, :rljl, :pomdp, :mcts],
    rljl_algorithm=:dqn,
    pomdp_solver=:value_iteration,
    mcts_n_simulations=1000
)
```

**Status**: Proof-of-concept - Four experimental RL approaches implemented for research purposes

**Use Cases (research only)**:
- Algorithm development and exploration
- Comparative studies of RL approaches to assembly
- Academic research into ML-guided assembly optimization

**Note**: These are research prototypes not recommended for production use. For production assembly, use the third-party assembler wrappers listed above.

## Quality Assessment and Metrics

### Assembly Quality Metrics
```julia
import Mycelia

# Comprehensive quality assessment [Implemented]
metrics = Mycelia.calculate_assembly_quality_metrics(qualmer_graph)

# Error detection [Partially Implemented]
potential_errors = Mycelia.identify_potential_errors(graph)

# Cross-validation comparison [Implemented]
comparison = Mycelia.compare_assembly_statistics(
    intelligent_results,
    iterative_results
)
```

### Available Metrics
- Mean k-mer coverage and quality (implemented)
- Joint probability confidence scores (implemented)
- Low-confidence k-mer fraction (implemented)
- Assembly accuracy when reference available (basic implementation)
- Read mapping statistics (planned)
- Contig N50 and other contiguity metrics (planned)

## Choosing the Right Approach

### Recommended Assembly Tools by Use Case

**For Production Assembly** - Use third-party assembler wrappers:

| Data Type | Read Length | Organism | Recommended Tool |
|-----------|-------------|----------|------------------|
| Illumina PE | Short | Isolate | SPAdes (`run_spades`) |
| Illumina PE | Short | Metagenome | MEGAHIT (`run_megahit`) or metaSPAdes (`run_metaspades`) |
| PacBio HiFi | Long | Isolate | hifiasm (`run_hifiasm`) or Flye (`run_flye`) |
| PacBio/ONT | Long | Metagenome | metaFlye (`run_metaflye`) or MetaMDBG (`run_metamdbg`) |
| Hybrid PE+Long | Mixed | Isolate | Unicycler (`run_unicycler`) |

**For Assembly Research** - Internal experimental algorithms:

| Research Goal | Recommended Approach |
|---------------|---------------------|
| Quality-aware assembly methods | Qualmer graphs and quality-aware assembly |
| Parameter optimization | Intelligent assembly |
| Statistical assembly approaches | Iterative assembly |
| ML-guided assembly | RL-guided assembly (proof-of-concept) |

### Performance Considerations

| Approach | Speed | Memory | Accuracy | Automation | Status |
|----------|-------|---------|----------|------------|--------|
| Standard | Fastest | Moderate | Good | Limited | Stable |
| Quality-Aware | Fast | Moderate | High | Moderate | Stable |
| Intelligent | Moderate | Moderate | Very High | Very High | Stable |
| Iterative | Slow | High | High | Very High | Stable |
| RL-Guided | Very Slow | Very High | Very High | Very High | Experimental |

## Integration Examples

### Complete Workflow Example
```julia
import Mycelia

# 1. Quality control
qc_results = Mycelia.assess_fastq_quality(fastq_file)

# 2. Choose assembly approach based on quality
if qc_results.mean_quality > 30
    # High quality - use direct approach
    assembly = Mycelia.build_quality_biosequence_graph(records)
else
    # Lower quality - use intelligent approach
    assembly = Mycelia.mycelia_assemble(fastq_file)
end

# 3. Validate assembly
validation = Mycelia.mycelia_cross_validation(records)

# 4. Extract final results
final_contigs = Mycelia.quality_biosequence_graph_to_fastq(assembly, "final")
```

### Benchmarking Workflow
```julia
import Mycelia

# Compare multiple approaches
approaches = [:intelligent, :iterative, :quality_aware]
benchmark_results = Mycelia.benchmark_assembly_approaches(
    test_datasets,
    approaches;
    metrics=[:accuracy, :speed, :memory]
)
```

## Best Practices

### General Guidelines
1. **Start Simple**: Use standard k-mer assembly for initial exploration
2. **Use Quality**: Incorporate quality scores when available
3. **Validate Results**: Always use cross-validation for important assemblies
4. **Monitor Resources**: Use memory limits to prevent system crashes
5. **Document Parameters**: Save assembly parameters for reproducibility

### Performance Optimization
1. **Choose Appropriate k**: Use sparsity detection for optimal k-mer sizes
2. **Memory Management**: Configure memory limits based on available resources
3. **Parallel Processing**: Utilize multi-threading where available
4. **Incremental Assembly**: Use checkpointing for long-running assemblies

### Quality Control
1. **Pre-Assembly QC**: Assess read quality before choosing approach
2. **Assembly Validation**: Use cross-validation and read mapping
3. **Error Detection**: Monitor low-confidence regions
4. **Comparative Analysis**: Compare results from multiple approaches

## Implementation Status and Future Directions

### Current Status
- **Stable**: Core graph-based assembly algorithms, quality-aware assembly, intelligent assembly
- **In Development**: Error correction algorithms, advanced validation metrics
- **Experimental**: Reinforcement learning approaches, POMDP integration
- **Planned**: Enhanced parallel processing, cloud integration, long-read optimizations

### Experimental Features Notice
**Experimental features** are research implementations that may:
- Require additional dependencies
- Have limited testing coverage
- Change significantly between versions
- Be computationally intensive

### Future Directions
The Mycelia assembly suite continues to evolve with:
- Enhanced machine learning integration (in development)
- Improved parallel processing capabilities (planned)
- Additional quality-aware algorithms (in development)
- Extended support for long-read technologies (planned)
- Integration with cloud computing platforms (planned)

For the latest developments, see the [Assembly Roadmap](https://github.com/cjprybol/Mycelia/blob/main/ASSEMBLY_ROADMAP.md).