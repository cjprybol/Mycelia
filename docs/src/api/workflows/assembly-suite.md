# Assembly Suite Overview

Mycelia provides a comprehensive assembly suite with multiple approaches for genome assembly optimization. This document provides an overview of the assembly ecosystem and guidance on choosing the right approach for your data.

> **Development Status**: The assembly suite includes both stable implementations and experimental features. Features are clearly marked with their current status throughout this documentation.

## Assembly Philosophy

Mycelia's assembly suite is built on several core principles:

- **Accuracy First**: Prioritize assembly accuracy over contiguity or speed
- **Quality Awareness**: Utilize quality scores for better assembly decisions
- **Graph Hierarchy**: Support 6 complementary graph types for different use cases
- **Type Safety**: Maintain proper BioSequence types throughout (no string conversions)
- **Modular Design**: Choose components based on your specific needs

## Graph Type Hierarchy

### Fixed-Length Graphs (Assembly Foundation)
1. **N-gram Graphs** (`build_ngram_graph`) - Unicode character analysis
2. **K-mer Graphs** (`build_kmer_graph_next`) - Standard k-mer assembly
3. **Qualmer Graphs** (`build_qualmer_graph`) - Quality-aware k-mer assembly

### Variable-Length Graphs (Simplified Products)
4. **String Graphs** (`string_to_ngram_graph`) - Variable unicode strings
5. **FASTA Graphs** (`build_biosequence_graph`) - Variable BioSequences
6. **FASTQ Graphs** (`build_quality_biosequence_graph`) - Quality-aware BioSequences

## Assembly Approaches

### 1. Standard Assembly (Traditional)

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

**Use Cases:**
- High-quality FASTA data
- Simple genomes without complex repeats
- When computational speed is critical

### 2. Quality-Aware Assembly

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

**Use Cases:**
- Illumina short reads
- Error-prone long reads  
- When maximum accuracy is required
- Complex genomes with repeats

### 3. Intelligent Self-Optimizing Assembly ✅ **Implemented**

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

**Features:**
- ✅ Automatic prime k-mer progression
- ✅ Sparsity-based k-mer selection
- ✅ Memory monitoring and limits
- ✅ Cross-validation for quality assessment

**Use Cases:**
- Unknown optimal parameters
- Production assembly pipelines
- When manual optimization is impractical

### 4. Iterative Statistical Assembly ✅ **Implemented**

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

**Features:**
- ✅ Read-level likelihood optimization
- ✅ Statistical path resampling
- ✅ Timestamped output for tracking evolution
- 🚧 Viterbi integration for optimal paths (partially implemented)

**Use Cases:**
- When individual read optimization is beneficial
- Datasets with systematic errors
- Research applications requiring detailed tracking

### 5. Reinforcement Learning Guided Assembly 🧪 **Experimental**

For machine learning optimization (four experimental implementations available):

#### Custom RL Framework 🧪 **Experimental - Under Development**
```julia
import Mycelia

# Train custom RL agent
agent = Mycelia.DQNPolicy()
env = Mycelia.AssemblyEnvironment(reads)
trained_agent = Mycelia.train_assembly_agent(env, agent; episodes=1000)

# Apply learned policy
results = Mycelia.apply_learned_policy(trained_agent, "genome.fastq")
```

#### ReinforcementLearning.jl Integration 🧪 **Experimental - Basic Implementation**
```julia
import Mycelia

# Use well-tested RL algorithms
assembly, history = Mycelia.intelligent_assembly_rljl(
    reads;
    algorithm=:dqn,
    train_episodes=1000
)
```

#### POMDPs.jl Integration 🧪 **Experimental - Basic Implementation**
```julia
import Mycelia

# Formal MDP/POMDP approach
assembly, history = Mycelia.intelligent_assembly_pomdp(
    reads;
    solver=:value_iteration,
    use_pomdp=false
)
```

#### Monte Carlo Tree Search (MCTS) 🧪 **Experimental - Game-Based Approach**
```julia
import Mycelia

# Assembly as a sequential decision game
assembly, trajectory = Mycelia.intelligent_assembly_mcts(
    reads;
    n_simulations=1000,
    exploration_constant=sqrt(2)
)
```

#### Comparison Framework 🧪 **Experimental - Proof of Concept**
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

**Use Cases:**
- Research into assembly optimization
- When you have training data for optimization
- Automated parameter learning from experience

## Quality Assessment and Metrics

### Assembly Quality Metrics
```julia
import Mycelia

# Comprehensive quality assessment ✅ **Implemented**
metrics = Mycelia.calculate_assembly_quality_metrics(qualmer_graph)

# Error detection 🚧 **Partially Implemented**
potential_errors = Mycelia.identify_potential_errors(graph)

# Cross-validation comparison ✅ **Implemented**
comparison = Mycelia.compare_assembly_statistics(
    intelligent_results,
    iterative_results
)
```

### Available Metrics
- ✅ Mean k-mer coverage and quality
- ✅ Joint probability confidence scores
- ✅ Low-confidence k-mer fraction
- 🚧 Assembly accuracy (when reference available - basic implementation)
- 📋 Read mapping statistics (planned)
- 📋 Contig N50 and other contiguity metrics (planned)

## Choosing the Right Approach

### Decision Matrix

| Data Type | Quality | Complexity | Recommended Approach |
|-----------|---------|------------|---------------------|
| FASTA, High Quality | N/A | Simple | Standard K-mer |
| FASTQ, High Quality | >Q30 | Simple | Direct Quality-Aware |
| FASTQ, Medium Quality | Q20-Q30 | Medium | Qualmer-Mediated |
| FASTQ, Low Quality | <Q20 | Complex | Intelligent + Iterative |
| Unknown Parameters | Any | Any | Intelligent Assembly |
| Research/Optimization | Any | Any | RL-Guided |

### Performance Considerations

| Approach | Speed | Memory | Accuracy | Automation | Status |
|----------|-------|---------|----------|------------|--------|
| Standard | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ✅ Stable |
| Quality-Aware | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ✅ Stable |
| Intelligent | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ✅ Stable |
| Iterative | ⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ✅ Stable |
| RL-Guided | ⭐ | ⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | 🧪 Experimental |

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
- ✅ **Stable**: Core graph-based assembly algorithms, quality-aware assembly, intelligent assembly
- 🚧 **In Development**: Error correction algorithms, advanced validation metrics
- 🧪 **Experimental**: Reinforcement learning approaches, POMDP integration
- 📋 **Planned**: Enhanced parallel processing, cloud integration, long-read optimizations

### Experimental Features Notice
🧪 **Experimental features** are research implementations that may:
- Require additional dependencies
- Have limited testing coverage
- Change significantly between versions
- Be computationally intensive

### Future Directions
The Mycelia assembly suite continues to evolve with:
- 🚧 Enhanced machine learning integration
- 📋 Improved parallel processing capabilities
- 🚧 Additional quality-aware algorithms
- 📋 Extended support for long-read technologies
- 📋 Integration with cloud computing platforms

For the latest developments, see the [Assembly Roadmap](../../ASSEMBLY_ROADMAP.md).