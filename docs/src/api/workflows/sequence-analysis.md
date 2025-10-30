# Sequence Analysis & K-mers

Functions for analyzing sequence composition, counting k-mers, and extracting genomic features from sequencing data.

## Overview

Sequence analysis forms the foundation of many bioinformatics algorithms. Mycelia provides:

- **K-mer counting and analysis** using efficient algorithms (implemented)
- **Genome size estimation** from k-mer frequency spectra (implemented)
- **Graph-based k-mer analysis** for assembly (implemented)
- **Error detection and correction** using k-mer patterns *(planned)*
- **Sequence composition analysis** and bias detection *(planned)*
- **Memory-efficient processing** of large datasets (partially implemented)

> **Implementation Status**: Core k-mer counting, canonical k-mer analysis, and basic genome size estimation are implemented. Advanced analysis features, error correction, and comprehensive visualization are planned.

## Common Workflows

### 1. Basic K-mer Analysis
```julia
# Count k-mers in sequencing reads
kmer_counts = Mycelia.count_kmers("reads.fastq", k=21)

# Analyze k-mer frequency spectrum
spectrum = Mycelia.kmer_frequency_spectrum(kmer_counts)
```

### 2. Genome Size Estimation (Implemented)

```julia
# Estimate genome size from k-mer spectrum
genome_size = Mycelia.estimate_genome_size_from_kmers(kmer_counts, coverage_peak=30)
```

### 3. Contamination Detection *(Planned)*

```julia
# Detect contamination using k-mer profiles (planned)
contamination = Mycelia.detect_contamination_kmers("reads.fastq", expected_profile)
```

## K-mer Counting

### Basic K-mer Counting

```@docs
Mycelia.count_kmers
Mycelia.count_canonical_kmers
Mycelia.fasta_list_to_dense_kmer_counts
Mycelia.fasta_list_to_sparse_kmer_counts
```

#### Example: K-mer Counting Comparison
```julia
# Compare different k-mer counting approaches
reads_file = "reads.fastq"

# Dense counting (stores all possible k-mers)
dense_counts = Mycelia.count_kmers(reads_file, k=15, method="dense")
println("Dense matrix size: $(size(dense_counts))")

# Sparse counting (only observed k-mers)
sparse_counts = Mycelia.count_kmers(reads_file, k=21, method="sparse")
println("Unique k-mers: $(length(sparse_counts))")

# Canonical k-mers (combines forward and reverse complement)
canonical_counts = Mycelia.count_canonical_kmers(reads_file, k=21)
println("Canonical k-mers: $(length(canonical_counts))")
```

### Advanced K-mer Analysis

<!-- build_kmer_graph, find_kmer_overlaps, extract_kmer_paths, analyze_kmer_connectivity not yet implemented as documented -->

```@docs
Mycelia.build_kmer_graph_next
```

**Note**: Additional k-mer graph construction functions (`build_stranded_kmer_graph`, `build_directed_kmer_graph`) have been deprecated and moved to `sequence-graphs.jl.deprecated`.

#### Example: K-mer Graph Construction
```julia
# Build k-mer overlap graph
kmer_graph = Mycelia.build_kmer_graph(
    sequences,
    k=31,
    min_overlap=20,
    min_coverage=5
)

# Analyze graph properties
graph_stats = Mycelia.analyze_kmer_connectivity(kmer_graph)
println("Graph nodes: $(graph_stats.n_nodes)")
println("Graph edges: $(graph_stats.n_edges)")
println("Connected components: $(graph_stats.n_components)")
```

## K-mer Frequency Analysis

### Frequency Spectra

<!-- kmer_frequency_spectrum, analyze_spectrum_peaks, fit_spectrum_model not yet implemented as documented -->

```@docs
Mycelia.jellyfish_counts_to_kmer_frequency_histogram
Mycelia.plot_kmer_frequency_spectra
```

#### Example: K-mer Spectrum Analysis
```julia
# Generate and analyze k-mer frequency spectrum
kmer_counts = Mycelia.count_kmers("reads.fastq", k=21)
spectrum = Mycelia.kmer_frequency_spectrum(kmer_counts)

# Identify characteristic peaks
peaks = Mycelia.analyze_spectrum_peaks(spectrum)
println("Coverage peak at frequency: $(peaks.main_peak)")
println("Error peak at frequency: $(peaks.error_peak)")

# Fit mathematical model to spectrum
model = Mycelia.fit_spectrum_model(spectrum, model_type="negative_binomial")
println("Model fit R²: $(model.r_squared)")
```

### Genome Characteristics from K-mers

<!-- estimate_genome_size_from_kmers, estimate_coverage_from_kmers, detect_repetitive_kmers, calculate_genome_complexity not yet implemented as documented -->

#### Example: Genome Size Estimation
```julia
# Estimate genome size using k-mer spectrum
kmer_counts = Mycelia.count_kmers("reads.fastq", k=21)
spectrum = Mycelia.kmer_frequency_spectrum(kmer_counts)

# Find coverage peak
coverage_peak = Mycelia.find_coverage_peak(spectrum)
println("Estimated coverage: $(coverage_peak)x")

# Calculate genome size
total_kmers = sum(spectrum.frequencies .* spectrum.counts)
genome_size = total_kmers ÷ coverage_peak
println("Estimated genome size: $(genome_size) bp")

# Alternative method with error correction
corrected_estimate = Mycelia.estimate_genome_size_from_kmers(
    kmer_counts,
    error_correction=true,
    ploidy=1
)
println("Corrected genome size: $(corrected_estimate.size) bp")
println("Confidence interval: $(corrected_estimate.confidence_interval)")
```

## Error Detection and Correction *(Planned)*

### K-mer Based Error Detection *(Planned)*

> **NOTE**: K-mer based error detection and correction functions are planned but not yet implemented.

#### Example: Error Detection *(Planned)*

```julia
# Identify likely error k-mers (planned)
error_kmers = Mycelia.identify_error_kmers(
    kmer_counts,
    min_coverage=3,
    max_coverage=100
)

println("Potential error k-mers: $(length(error_kmers))")
println("Error rate estimate: $(length(error_kmers) / length(kmer_counts) * 100)%")

# Correct errors in reads
corrected_reads = Mycelia.correct_sequencing_errors(
    "reads.fastq",
    error_kmers,
    correction_method="consensus"
)

# Validate correction effectiveness
validation = Mycelia.validate_error_correction(
    "reads.fastq",
    corrected_reads,
    reference_genome="reference.fasta"
)
println("Error reduction: $(validation.error_reduction_percent)%")
```

## Sequence Composition Analysis *(Planned)*

### Nucleotide Composition *(Planned)*

> **NOTE**: Comprehensive composition analysis functions are planned but not yet implemented.

#### Example: Comprehensive Composition Analysis *(Planned)*

```julia
# Analyze nucleotide composition (planned)
composition = Mycelia.calculate_nucleotide_frequencies("sequences.fasta")
println("GC content: $(composition.gc_content)%")
println("AT content: $(composition.at_content)%")

# Dinucleotide analysis
dinuc_freq = Mycelia.analyze_dinucleotide_frequencies("sequences.fasta")
println("CpG frequency: $(dinuc_freq.CG)")
println("Expected CpG: $(dinuc_freq.expected_CG)")
println("CpG O/E ratio: $(dinuc_freq.CG_oe_ratio)")

# Detect composition bias
bias_analysis = Mycelia.detect_sequence_bias(composition)
if bias_analysis.bias_detected
    println("Sequence bias detected:")
    println("  Type: $(bias_analysis.bias_type)")
    println("  Strength: $(bias_analysis.bias_strength)")
end
```

### Codon Usage Analysis

<!-- calculate_codon_usage, analyze_codon_bias, compare_codon_usage, optimize_codon_usage not yet implemented as documented -->

#### Example: Codon Usage Analysis
```julia
# Analyze codon usage in coding sequences
cds_file = "coding_sequences.fasta"
codon_usage = Mycelia.calculate_codon_usage(cds_file, genetic_code="standard")

# Calculate codon bias metrics
bias_metrics = Mycelia.analyze_codon_bias(codon_usage)
println("Effective Number of Codons (ENC): $(bias_metrics.enc)")
println("Codon Adaptation Index (CAI): $(bias_metrics.cai)")
println("Frequency of Optimal Codons (FOP): $(bias_metrics.fop)")

# Compare with reference organism
reference_usage = Mycelia.load_reference_codon_usage("escherichia_coli")
comparison = Mycelia.compare_codon_usage(codon_usage, reference_usage)
println("Similarity to E. coli: $(comparison.similarity_score)")
```

## Memory-Efficient Processing

### Large Dataset Handling

<!-- stream_kmer_counting, process_kmers_in_chunks, parallel_kmer_analysis, memory_efficient_counting not yet implemented as documented -->

#### Example: Large File Processing
```julia
# Process large files with memory constraints
large_file = "large_dataset.fastq"

# Stream processing for memory efficiency
kmer_counts = Mycelia.stream_kmer_counting(
    large_file,
    k=21,
    chunk_size=100000,
    memory_limit_gb=8
)

# Parallel processing for speed
parallel_counts = Mycelia.parallel_kmer_analysis(
    large_file,
    k=21,
    n_workers=8,
    merge_strategy="sum"
)
```

### Optimized Data Structures

<!-- create_sparse_kmer_matrix, build_compressed_kmer_index, use_bloom_filter_counting, implement_count_min_sketch not yet implemented as documented -->

#### Example: Memory Optimization
```julia
# Use probabilistic data structures for very large datasets
bloom_filter = Mycelia.create_kmer_bloom_filter(
    estimated_kmers=1_000_000_000,
    false_positive_rate=0.01
)

# Streaming k-mer presence testing
for read in Mycelia.stream_fastq("huge_dataset.fastq")
    for kmer in Mycelia.extract_kmers(read, k=21)
        if Mycelia.probably_present(bloom_filter, kmer)
            # Process likely present k-mer
            Mycelia.process_kmer(kmer)
        end
    end
end
```

## Comparative K-mer Analysis *(Planned)*

### Multi-Sample Analysis *(Planned)*

> **NOTE**: Multi-sample k-mer comparison functions are planned but not yet implemented.

#### Example: Multi-Sample K-mer Comparison *(Planned)*

```julia
# Compare k-mer profiles across multiple samples (planned)
sample_files = ["sample1.fastq", "sample2.fastq", "sample3.fastq"]

# Build k-mer profiles for each sample
profiles = [Mycelia.count_kmers(file, k=21) for file in sample_files]

# Calculate pairwise distances
distance_matrix = Mycelia.build_kmer_distance_matrix(profiles, metric="jaccard")
println("Pairwise k-mer distances:")
display(distance_matrix)

# Cluster samples by k-mer similarity
clusters = Mycelia.cluster_by_kmer_similarity(profiles, method="hierarchical")
println("Sample clusters: $(clusters)")
```

### Contamination Detection

<!-- detect_contamination_kmers, identify_foreign_kmers, classify_kmer_sources, remove_contaminating_kmers not yet implemented as documented -->

#### Example: Contamination Detection
```julia
# Detect contamination using k-mer profiles
sample_kmers = Mycelia.count_kmers("sample.fastq", k=21)
reference_kmers = Mycelia.count_kmers("reference_genome.fasta", k=21)

# Identify foreign k-mers
foreign_kmers = Mycelia.identify_foreign_kmers(
    sample_kmers,
    reference_kmers,
    min_abundance=5
)

contamination_rate = length(foreign_kmers) / length(sample_kmers)
println("Contamination rate: $(contamination_rate * 100)%")

# Classify contamination sources
contamination_sources = Mycelia.classify_kmer_sources(
    foreign_kmers,
    database_kmers=["human", "bacterial", "viral"]
)

for (source, proportion) in contamination_sources
    println("$(source): $(proportion * 100)%")
end
```

## Visualization and Reporting

### K-mer Plots

<!-- plot_kmer_spectrum, plot_kmer_composition, plot_genome_size_estimation, create_kmer_dashboard not yet implemented as individual functions -->

#### Example: K-mer Visualization
```julia
# Create comprehensive k-mer analysis plots
kmer_counts = Mycelia.count_kmers("reads.fastq", k=21)

# K-mer frequency spectrum
spectrum_plot = Mycelia.plot_kmer_spectrum(kmer_counts, 
                                  title="21-mer Frequency Spectrum",
                                  log_scale=true)

# Genome size estimation plot
size_plot = Mycelia.plot_genome_size_estimation(kmer_counts,
                                       show_confidence_interval=true)

# Combined dashboard
kmer_dashboard = Mycelia.create_kmer_dashboard(kmer_counts)
Mycelia.save_plot(kmer_dashboard, "kmer_analysis.png")
```

## Performance Considerations

### Algorithm Selection
- **Dense matrices**: Use for small k (k ≤ 12) or small genomes
- **Sparse matrices**: Use for large k (k ≥ 15) or large genomes
- **Streaming**: Use for memory-constrained environments
- **Parallel processing**: Use for time-critical applications

### Memory Usage
- **k=15**: ~1 GB for dense counting
- **k=21**: ~17 GB for dense counting (use sparse)
- **k=31**: Sparse counting only
- **Large genomes**: Consider streaming or chunked processing

### Computational Complexity
- **Time**: O(n) for sequence length n
- **Space**: O(4^k) for dense, O(unique k-mers) for sparse
- **Parallel**: Near-linear speedup for independent samples

## Common Issues and Solutions

### Memory Limitations
```julia
# Handle memory constraints
if Mycelia.estimate_memory_usage(file, k=21) > Mycelia.available_memory()
    # Use streaming approach
    kmer_counts = Mycelia.stream_kmer_counting(file, k=21, chunk_size=50000)
else
    # Use standard approach
    kmer_counts = Mycelia.count_kmers(file, k=21)
end
```

### Parameter Selection
```julia
# Choose optimal k-mer size
optimal_k = Mycelia.select_optimal_k(
    genome_size_estimate=5_000_000,
    error_rate=0.01,
    coverage=30
)
println("Recommended k-mer size: $(optimal_k)")
```

## Related Functions

### Data Structures *(planned)*
- `KmerCounts` - K-mer count data structure *(planned)*
- `KmerSpectrum` - Frequency spectrum representation *(planned)*
- `KmerGraph` - K-mer overlap graph *(planned)*

### File I/O
- [`Mycelia.open_fastx`](@ref) - Read sequence files
- [`Mycelia.save_kmer_results`](@ref) - Save k-mer counts
- [`Mycelia.load_kmer_results`](@ref) - Load saved counts

### Related Workflows

### Previous Steps
- [Quality Control](quality-control.md) - Preprocess reads before k-mer analysis

### Next Steps
- Genome Assembly *(planned)* - Use k-mer analysis for assembly
- Comparative Genomics *(planned)* - Compare k-mer profiles

## See Also
- [Tutorial 3: K-mer Analysis](../../generated/tutorials/03_kmer_analysis.md)
- Assembly Workflows *(planned)*
- [Performance Optimization](../examples/advanced-usage.md)