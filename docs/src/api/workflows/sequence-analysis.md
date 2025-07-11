# Sequence Analysis & K-mers

Functions for analyzing sequence composition, counting k-mers, and extracting genomic features from sequencing data.

## Overview

Sequence analysis forms the foundation of many bioinformatics algorithms. Mycelia provides comprehensive tools for:

- **K-mer counting and analysis** using efficient algorithms
- **Genome size estimation** from k-mer frequency spectra
- **Error detection and correction** using k-mer patterns
- **Sequence composition analysis** and bias detection
- **Memory-efficient processing** of large datasets

## Common Workflows

### 1. Basic K-mer Analysis
```julia
# Count k-mers in sequencing reads
kmer_counts = count_kmers("reads.fastq", k=21)

# Analyze k-mer frequency spectrum
spectrum = kmer_frequency_spectrum(kmer_counts)
```

### 2. Genome Size Estimation
```julia
# Estimate genome size from k-mer spectrum
genome_size = estimate_genome_size_from_kmers(kmer_counts, coverage_peak=30)
```

### 3. Contamination Detection
```julia
# Detect contamination using k-mer profiles
contamination = detect_contamination_kmers("reads.fastq", expected_profile)
```

## K-mer Counting

### Basic K-mer Counting

```@docs
count_kmers
count_canonical_kmers
fasta_list_to_dense_kmer_counts
fasta_list_to_sparse_kmer_counts
```

#### Example: K-mer Counting Comparison
```julia
# Compare different k-mer counting approaches
reads_file = "reads.fastq"

# Dense counting (stores all possible k-mers)
dense_counts = count_kmers(reads_file, k=15, method="dense")
println("Dense matrix size: $(size(dense_counts))")

# Sparse counting (only observed k-mers)
sparse_counts = count_kmers(reads_file, k=21, method="sparse")
println("Unique k-mers: $(length(sparse_counts))")

# Canonical k-mers (combines forward and reverse complement)
canonical_counts = count_canonical_kmers(reads_file, k=21)
println("Canonical k-mers: $(length(canonical_counts))")
```

### Advanced K-mer Analysis

```@docs
build_kmer_graph
find_kmer_overlaps
extract_kmer_paths
analyze_kmer_connectivity
```

#### Example: K-mer Graph Construction
```julia
# Build k-mer overlap graph
kmer_graph = build_kmer_graph(
    sequences,
    k=31,
    min_overlap=20,
    min_coverage=5
)

# Analyze graph properties
graph_stats = analyze_kmer_connectivity(kmer_graph)
println("Graph nodes: $(graph_stats.n_nodes)")
println("Graph edges: $(graph_stats.n_edges)")
println("Connected components: $(graph_stats.n_components)")
```

## K-mer Frequency Analysis

### Frequency Spectra

```@docs
kmer_frequency_spectrum
analyze_spectrum_peaks
fit_spectrum_model
plot_kmer_spectrum
```

#### Example: K-mer Spectrum Analysis
```julia
# Generate and analyze k-mer frequency spectrum
kmer_counts = count_kmers("reads.fastq", k=21)
spectrum = kmer_frequency_spectrum(kmer_counts)

# Identify characteristic peaks
peaks = analyze_spectrum_peaks(spectrum)
println("Coverage peak at frequency: $(peaks.main_peak)")
println("Error peak at frequency: $(peaks.error_peak)")

# Fit mathematical model to spectrum
model = fit_spectrum_model(spectrum, model_type="negative_binomial")
println("Model fit R²: $(model.r_squared)")
```

### Genome Characteristics from K-mers

```@docs
estimate_genome_size_from_kmers
estimate_coverage_from_kmers
detect_repetitive_kmers
calculate_genome_complexity
```

#### Example: Genome Size Estimation
```julia
# Estimate genome size using k-mer spectrum
kmer_counts = count_kmers("reads.fastq", k=21)
spectrum = kmer_frequency_spectrum(kmer_counts)

# Find coverage peak
coverage_peak = find_coverage_peak(spectrum)
println("Estimated coverage: $(coverage_peak)x")

# Calculate genome size
total_kmers = sum(spectrum.frequencies .* spectrum.counts)
genome_size = total_kmers ÷ coverage_peak
println("Estimated genome size: $(genome_size) bp")

# Alternative method with error correction
corrected_estimate = estimate_genome_size_from_kmers(
    kmer_counts,
    error_correction=true,
    ploidy=1
)
println("Corrected genome size: $(corrected_estimate.size) bp")
println("Confidence interval: $(corrected_estimate.confidence_interval)")
```

## Error Detection and Correction

### K-mer Based Error Detection

```@docs
identify_error_kmers
correct_sequencing_errors
validate_error_correction
```

#### Example: Error Detection
```julia
# Identify likely error k-mers
error_kmers = identify_error_kmers(
    kmer_counts,
    min_coverage=3,
    max_coverage=100
)

println("Potential error k-mers: $(length(error_kmers))")
println("Error rate estimate: $(length(error_kmers) / length(kmer_counts) * 100)%")

# Correct errors in reads
corrected_reads = correct_sequencing_errors(
    "reads.fastq",
    error_kmers,
    correction_method="consensus"
)

# Validate correction effectiveness
validation = validate_error_correction(
    "reads.fastq",
    corrected_reads,
    reference_genome="reference.fasta"
)
println("Error reduction: $(validation.error_reduction_percent)%")
```

## Sequence Composition Analysis

### Nucleotide Composition

```@docs
calculate_nucleotide_frequencies
analyze_dinucleotide_frequencies
calculate_codon_usage
detect_sequence_bias
```

#### Example: Comprehensive Composition Analysis
```julia
# Analyze nucleotide composition
composition = calculate_nucleotide_frequencies("sequences.fasta")
println("GC content: $(composition.gc_content)%")
println("AT content: $(composition.at_content)%")

# Dinucleotide analysis
dinuc_freq = analyze_dinucleotide_frequencies("sequences.fasta")
println("CpG frequency: $(dinuc_freq.CG)")
println("Expected CpG: $(dinuc_freq.expected_CG)")
println("CpG O/E ratio: $(dinuc_freq.CG_oe_ratio)")

# Detect composition bias
bias_analysis = detect_sequence_bias(composition)
if bias_analysis.bias_detected
    println("Sequence bias detected:")
    println("  Type: $(bias_analysis.bias_type)")
    println("  Strength: $(bias_analysis.bias_strength)")
end
```

### Codon Usage Analysis

```@docs
calculate_codon_usage
analyze_codon_bias
compare_codon_usage
optimize_codon_usage
```

#### Example: Codon Usage Analysis
```julia
# Analyze codon usage in coding sequences
cds_file = "coding_sequences.fasta"
codon_usage = calculate_codon_usage(cds_file, genetic_code="standard")

# Calculate codon bias metrics
bias_metrics = analyze_codon_bias(codon_usage)
println("Effective Number of Codons (ENC): $(bias_metrics.enc)")
println("Codon Adaptation Index (CAI): $(bias_metrics.cai)")
println("Frequency of Optimal Codons (FOP): $(bias_metrics.fop)")

# Compare with reference organism
reference_usage = load_reference_codon_usage("escherichia_coli")
comparison = compare_codon_usage(codon_usage, reference_usage)
println("Similarity to E. coli: $(comparison.similarity_score)")
```

## Memory-Efficient Processing

### Large Dataset Handling

```@docs
stream_kmer_counting
process_kmers_in_chunks
parallel_kmer_analysis
memory_efficient_counting
```

#### Example: Large File Processing
```julia
# Process large files with memory constraints
large_file = "large_dataset.fastq"

# Stream processing for memory efficiency
kmer_counts = stream_kmer_counting(
    large_file,
    k=21,
    chunk_size=100000,
    memory_limit_gb=8
)

# Parallel processing for speed
parallel_counts = parallel_kmer_analysis(
    large_file,
    k=21,
    n_workers=8,
    merge_strategy="sum"
)
```

### Optimized Data Structures

```@docs
create_sparse_kmer_matrix
build_compressed_kmer_index
use_bloom_filter_counting
implement_count_min_sketch
```

#### Example: Memory Optimization
```julia
# Use probabilistic data structures for very large datasets
bloom_filter = create_kmer_bloom_filter(
    estimated_kmers=1_000_000_000,
    false_positive_rate=0.01
)

# Streaming k-mer presence testing
for read in stream_fastq("huge_dataset.fastq")
    for kmer in extract_kmers(read, k=21)
        if probably_present(bloom_filter, kmer)
            # Process likely present k-mer
            process_kmer(kmer)
        end
    end
end
```

## Comparative K-mer Analysis

### Multi-Sample Analysis

```@docs
compare_kmer_profiles
build_kmer_distance_matrix
cluster_by_kmer_similarity
identify_sample_specific_kmers
```

#### Example: Multi-Sample K-mer Comparison
```julia
# Compare k-mer profiles across multiple samples
sample_files = ["sample1.fastq", "sample2.fastq", "sample3.fastq"]

# Build k-mer profiles for each sample
profiles = [count_kmers(file, k=21) for file in sample_files]

# Calculate pairwise distances
distance_matrix = build_kmer_distance_matrix(profiles, metric="jaccard")
println("Pairwise k-mer distances:")
display(distance_matrix)

# Cluster samples by k-mer similarity
clusters = cluster_by_kmer_similarity(profiles, method="hierarchical")
println("Sample clusters: $(clusters)")
```

### Contamination Detection

```@docs
detect_contamination_kmers
identify_foreign_kmers
classify_kmer_sources
remove_contaminating_kmers
```

#### Example: Contamination Detection
```julia
# Detect contamination using k-mer profiles
sample_kmers = count_kmers("sample.fastq", k=21)
reference_kmers = count_kmers("reference_genome.fasta", k=21)

# Identify foreign k-mers
foreign_kmers = identify_foreign_kmers(
    sample_kmers,
    reference_kmers,
    min_abundance=5
)

contamination_rate = length(foreign_kmers) / length(sample_kmers)
println("Contamination rate: $(contamination_rate * 100)%")

# Classify contamination sources
contamination_sources = classify_kmer_sources(
    foreign_kmers,
    database_kmers=["human", "bacterial", "viral"]
)

for (source, proportion) in contamination_sources
    println("$(source): $(proportion * 100)%")
end
```

## Visualization and Reporting

### K-mer Plots

```@docs
plot_kmer_spectrum
plot_kmer_composition
plot_genome_size_estimation
create_kmer_dashboard
```

#### Example: K-mer Visualization
```julia
# Create comprehensive k-mer analysis plots
kmer_counts = count_kmers("reads.fastq", k=21)

# K-mer frequency spectrum
spectrum_plot = plot_kmer_spectrum(kmer_counts, 
                                  title="21-mer Frequency Spectrum",
                                  log_scale=true)

# Genome size estimation plot
size_plot = plot_genome_size_estimation(kmer_counts,
                                       show_confidence_interval=true)

# Combined dashboard
kmer_dashboard = create_kmer_dashboard(kmer_counts)
save_plot(kmer_dashboard, "kmer_analysis.png")
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
if estimate_memory_usage(file, k=21) > available_memory()
    # Use streaming approach
    kmer_counts = stream_kmer_counting(file, k=21, chunk_size=50000)
else
    # Use standard approach
    kmer_counts = count_kmers(file, k=21)
end
```

### Parameter Selection
```julia
# Choose optimal k-mer size
optimal_k = select_optimal_k(
    genome_size_estimate=5_000_000,
    error_rate=0.01,
    coverage=30
)
println("Recommended k-mer size: $(optimal_k)")
```

## Related Functions

### Data Structures
- [`KmerCounts`](@ref) - K-mer count data structure
- [`KmerSpectrum`](@ref) - Frequency spectrum representation
- [`KmerGraph`](@ref) - K-mer overlap graph

### File I/O
- [`read_fastq`](@ref) - Read sequence files
- [`save_kmer_counts`](@ref) - Save k-mer counts
- [`load_kmer_counts`](@ref) - Load saved counts

### Related Workflows

### Previous Steps
- [Quality Control](quality-control.md) - Preprocess reads before k-mer analysis

### Next Steps
- [Genome Assembly](genome-assembly.md) - Use k-mer analysis for assembly
- [Comparative Genomics](comparative-genomics.md) - Compare k-mer profiles

## See Also
- [Tutorial 3: K-mer Analysis](../../tutorials/03_kmer_analysis.md)
- [Assembly Workflows](genome-assembly.md)
- [Performance Optimization](../examples/advanced-usage.md)