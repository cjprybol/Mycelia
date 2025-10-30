# Theoretical Foundations of Mycelia Assembly Algorithms

Mycelia's assembly algorithms are built on rigorous mathematical foundations derived from extensive research documented in the Mycelia-Dev notebooks. This document provides the theoretical background, design rationales, and mathematical principles that guide the implementation choices in Mycelia.

## Mathematical Framework for Assembly

### Core Assembly Paradigm: Probabilistic Expectation-Maximization

Mycelia treats genome assembly as a probabilistic inference problem, unifying the strengths of both De Bruijn graph and overlap-layout-consensus approaches through a generalized framework:

**Algorithm Overview:**
1. Create De Bruijn assembly graph from sequencing reads
2. Convert graph structure into probabilistic Hidden Markov Model (HMM)
3. Apply Viterbi algorithm for maximum likelihood error correction
4. Iterate until convergence at current k-mer size
5. Increment k-mer size and repeat until final convergence

**Theoretical Foundation:**
The approach treats assembly as an **Expectation-Maximization (EM) problem**:
- **E-step:** Calculate maximum likelihood paths through the assembly graph
- **M-step:** Update graph structure based on error-corrected reads

This probabilistic framework provides mathematical rigor that replaces heuristic graph cleaning decisions with principled, optimal solutions.

## K-mer Size Selection Theory

### Error Rate Mathematical Foundation

The optimal k-mer size selection follows mathematically derived principles:

**Lower Bound Formula:**
```
lower_bound_k = max(3, floor(1/error_rate - 1))
```

**Rationale:** This ensures that error-prone k-mers remain distinguishable from true genomic k-mers based on frequency distributions.

**Examples:**
- 1% error rate → k ≥ 99
- 5% error rate → k ≥ 19  
- 10% error rate → k ≥ 9

### Sequence Length Optimization

**Log₄ Pattern Discovery:**
Empirical analysis revealed that optimal starting k-mer sizes follow a log₄(sequence_length) pattern:

- 100bp sequences → k = 3
- 1,000bp sequences → k = 5
- 10,000bp sequences → k = 7
- 100,000bp sequences → k = 8

**Biological Rationale:** This pattern corresponds to the divergence point where erroneous k-mers begin to dominate true k-mers in the frequency spectrum.

### Coverage vs. Error Rate Trade-offs

**Key Finding:** A 10x increase in sequencing coverage contributes more noise to the k-mer spectrum than a 10x increase in error rate.

**Implication:** Conservative error rate assumptions with moderate coverage are preferable to high coverage with relaxed quality filtering.

## Dynamic Prime Pattern Algorithm

### Mathematical Elegance of Prime Selection

Mycelia employs a dynamic k-mer selection algorithm that exploits the mathematical properties of prime number distribution:

**Algorithm:**
```julia
function dynamic_k_prime_pattern(start_prime=11, max_k=101, initial_step=2)
    k_sequence = [start_prime]
    current_k = start_prime
    step = initial_step
    
    while true
        next_k = current_k + step
        if next_k > max_k || !isprime(next_k)
            break
        end
        push!(k_sequence, next_k)
        current_k = next_k
        step += 2  # Progressive spacing
    end
    
    return k_sequence
end
```

**Advantages:**
1. **Twin Prime Avoidance:** Automatically skips one member of twin prime pairs, reducing redundant analysis
2. **Progressive Spacing:** Increasing gaps minimize computational overlap while maintaining coverage
3. **Hardware Optimization:** k=31 fits optimally in 64-bit integers (2 bits per base × 32 bases)

### Biological Rationale for Prime K-mers

**Theoretical Foundation:** Prime k-mer lengths cannot form perfect repeats, reducing spurious matches and improving assembly specificity. Additionally, odd k-mers cannot be reverse complements of themselves, avoiding palindromic complications.

## Graph Theory Foundations

### Explicit vs. Inferred Edge Storage

**Core Principle:** Mycelia stores only edges that were actually observed in the data, not all theoretically possible k-mer neighbors.

**Mathematical Justification:**
- As k increases, the probability of unobserved neighboring k-mers existing decreases asymptotically toward zero
- For small k and large datasets, false positive edge rates can be extremely high
- Explicit storage provides O(n) space complexity vs. O(4^k) for complete neighbor graphs

### Statistical Tip Clipping Algorithm

**Methodology:** Remove graph tips using statistical thresholds relative to connected component coverage distributions.

**Algorithm:**
```julia
function statistical_tip_clipping(graph; std_dev_multiplier=3.0)
    for component in connected_components(graph)
        coverage_values = get_coverage_values(graph, component)
        median_coverage = median(coverage_values)
        coverage_std = std(coverage_values)
        threshold = median_coverage - std_dev_multiplier * coverage_std
        
        for tip in identify_tips(component)
            if get_coverage(tip) < threshold || get_coverage(tip) == 1
                remove_node!(graph, tip)
            end
        end
    end
end
```

**Statistical Foundation:** Uses 3σ rule to identify coverage outliers, preserving high-confidence tips while removing error artifacts.

## Metagenomic Classification Theory

### MAPQ Score Reinterpretation

**Critical Insight:** In metagenomic contexts, MAPQ=0 indicates taxonomic ambiguity rather than poor alignment quality.

**Theoretical Correction:**
- Traditional approach: Filter out MAPQ=0 alignments
- Mycelia approach: Retain MAPQ=0 alignments, weight by alignment scores

**Mathematical Framework:**
```julia
function weighted_taxonomic_assignment(alignments)
    weighted_assignments = []
    for alignment in alignments
        weight = alignment.score  # NOT alignment.mapq
        taxon = get_taxonomy(alignment.reference)
        push!(weighted_assignments, (taxon, weight))
    end
    return normalize_weights(weighted_assignments)
end
```

**Impact:** Prevents loss of biologically meaningful information from closely related taxa that receive MAPQ=0 due to mapping ambiguity.

### Strain-Level Clustering with FastANI

**Threshold Selection:** 99.5% Average Nucleotide Identity (ANI) threshold for strain-level clustering.

**Biological Rationale:** This threshold corresponds to the conventional species boundary, enabling strain-resolution within species while maintaining taxonomic coherence.

## Assembly Quality Assessment

### Graph Distance Metrics

**Approach:** Use Jaccard similarity between k-mer sets and edge sets to quantify assembly accuracy.

**Metrics:**
```julia
kmer_distance = 1 - jaccard(assembly_kmers, reference_kmers)
edge_distance = 1 - jaccard(assembly_edges, reference_edges)
```

**Theoretical Advantage:** Provides quantitative, parameter-free assessment of assembly quality relative to ground truth.

### Viterbi-Based Error Correction

**Hidden Markov Model Framework:**
- **States:** K-mer positions with associated quality scores
- **Transitions:** Observed k-mer adjacencies in assembly graph
- **Emissions:** Quality score-weighted k-mer observations

**Likelihood Calculation:**
```julia
error_probability = 10^(quality_score / -10)
state_likelihood = (1 - error_probability)^matches × error_probability^mismatches
```

## Genomic Grammar and Zipf's Law

### Linguistic Properties of Genomic Sequences

**Theoretical Discovery:** Genomic sequences follow Zipf's law (power-law distribution) similar to natural languages.

**Implications:**
1. **Frequency Distributions:** K-mer frequencies exhibit log-log linear patterns
2. **Biological Grammar:** Non-coding regions particularly follow linguistic-like structures
3. **Assembly Strategy:** Rare k-mers may represent genuine biological signal rather than noise

**Mathematical Model:**
```
frequency(k-mer) ∝ rank^(-α)
```
where α ≈ 1 for natural languages and genomic sequences.

## Strain-Resolved Assembly Framework

### Comprehensive Quality Metrics

**Strain-Specific Metrics:**
- `strain_recall = true_strains_recovered / total_true_strains`
- `strain_precision = correct_strain_calls / total_strain_calls`
- `strain_f1_score = 2 × (precision × recall) / (precision + recall)`

**Assembly Contiguity:**
- `NGA50`: N50 based on aligned contigs to reference
- `largest_alignment`: Longest single alignment length
- `total_aligned_length`: Total successfully aligned sequence

**Technology-Specific Assessments:**
- `homopolymer_accuracy`: Accuracy in homopolymer-rich regions
- `repeat_resolution_rate`: Fraction of repetitive elements properly resolved

### Confidence-Aware Assembly

**Per-Base Confidence Scoring:**
Each assembled base receives a confidence score based on:
1. Supporting read depth
2. Base quality scores of supporting reads
3. Graph topology confidence (branching vs. linear regions)
4. Strain assignment probability

## Computational Complexity Considerations

### Progressive K-mer Strategy

**Complexity Reduction:** Starting with small k-values and progressively increasing reduces computational complexity from O(4^k) to O(log k) in practice.

**Convergence Criteria:**
1. No corrections made in current iteration
2. Graph structure stabilizes between iterations
3. Quality metrics plateau across k-increments

### Memory-Efficient Graph Representation

**Canonical K-mer Storage:** Store only canonical form of each k-mer, reducing memory by ~50%.

**Sparse Edge Representation:** Use adjacency lists rather than dense matrices, achieving O(E) vs. O(V²) space complexity.

## Integration with Modern Sequencing Technologies

### Long-Read Assembly Considerations

**Overlap-Based Framework:** For error-prone long reads, minimum overlap length follows:
```
minimum_overlap = log₄(genome_size)
```

**Hybrid Strategy:** Combine exact k-mer matching (short reads) with approximate overlaps (long reads) using probabilistic framework.

### Real-Time Assembly

**Streaming Algorithm Design:** Progressive assembly enables real-time processing of sequencing data as it becomes available.

**Hardware Acceleration:** GPU/FPGA implementations can accelerate Viterbi algorithm computations for production-scale applications.

## Validation and Benchmarking

### Cross-Validation Framework

**Statistical Validation:** Bootstrap resampling with confidence interval calculation for robust parameter optimization.

**Multi-Scale Validation:** Simultaneous analysis using prime k-mer sizes (3,5,7,11,13,17,19) for comprehensive assessment.

### Benchmarking Methodology

**Ground Truth Comparison:** Use simulated data with known ground truth for algorithm validation.

**Comparative Analysis:** Quantitative comparison against established assemblers using standardized metrics.

## Research Foundations

This theoretical framework is derived from extensive research documented in the Mycelia-Dev repository, representing several years of algorithm development and validation. The mathematical foundations provide:

1. **Reproducible Results:** Deterministic algorithms with well-defined parameters
2. **Scalable Solutions:** Theoretical framework adapts to different sequencing technologies
3. **Optimal Performance:** Maximum likelihood approaches replace heuristic approximations
4. **Biological Relevance:** Algorithms account for genuine biological sequence properties

The integration of probability theory, graph algorithms, information theory, and computational biology provides a comprehensive foundation for next-generation genome assembly tools that achieve both accuracy and biological relevance.