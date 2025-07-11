# Quality Control & Preprocessing

Functions for assessing and improving sequencing data quality before downstream analysis.

## Overview

Quality control is essential for reliable bioinformatics results. Mycelia provides comprehensive tools for:

- **Assessing sequencing data quality** using multiple metrics
- **Preprocessing reads** to improve data quality
- **Identifying and removing contaminants**
- **Filtering and trimming** based on quality thresholds
- **Generating quality reports** for documentation

## Common Workflows

### 1. Basic Quality Assessment
```julia
# Analyze FASTQ quality
quality_stats = analyze_fastq_quality("reads.fastq")

# Generate quality report
quality_report = generate_quality_report(quality_stats)
```

### 2. Read Preprocessing
```julia
# Quality trimming and filtering
clean_reads = preprocess_reads(
    "raw_reads.fastq",
    min_quality=20,
    min_length=1000,
    trim_ends=true
)
```

### 3. Contamination Removal
```julia
# Remove host contamination
decontaminated_reads = remove_contamination(
    "reads.fastq",
    reference_genome="host_genome.fasta",
    output="clean_reads.fastq"
)
```

## Quality Assessment

### FASTQ Quality Analysis

```@docs
analyze_fastq_quality
calculate_per_base_quality
calculate_per_sequence_quality
assess_quality_degradation
```

#### Example: Comprehensive Quality Analysis
```julia
# Analyze quality across multiple metrics
quality_data = analyze_fastq_quality("reads.fastq")

# Access quality metrics
println("Total reads: $(quality_data.n_reads)")
println("Mean quality: $(quality_data.mean_quality)")
println("Mean length: $(quality_data.mean_length)")
println("GC content: $(quality_data.gc_content)%")

# Quality score distribution
quality_dist = quality_data.quality_distribution
println("Q20+ reads: $(quality_dist.q20_percent)%")
println("Q30+ reads: $(quality_dist.q30_percent)%")
```

### Sequence Composition Analysis

```@docs
calculate_gc_content
analyze_nucleotide_composition
detect_composition_bias
calculate_complexity_scores
```

#### Example: Composition Analysis
```julia
# Analyze sequence composition
composition = analyze_nucleotide_composition("reads.fastq")

# Check for bias
bias_detected = detect_composition_bias(composition)
if bias_detected.has_bias
    println("Composition bias detected:")
    println("  Type: $(bias_detected.bias_type)")
    println("  Severity: $(bias_detected.severity)")
    println("  Affected positions: $(bias_detected.positions)")
end
```

### Read Length Analysis

```@docs
calculate_read_length_distribution
analyze_length_uniformity
detect_length_artifacts
```

#### Example: Length Distribution Analysis
```julia
# Analyze read length characteristics
length_stats = calculate_read_length_distribution("reads.fastq")

println("Read length statistics:")
println("  Mean: $(length_stats.mean)")
println("  Median: $(length_stats.median)")
println("  Min: $(length_stats.minimum)")
println("  Max: $(length_stats.maximum)")
println("  Std Dev: $(length_stats.std)")

# Check for length artifacts
artifacts = detect_length_artifacts(length_stats)
if !isempty(artifacts)
    println("Length artifacts detected: $(length(artifacts))")
end
```

## Preprocessing and Filtering

### Quality-Based Filtering

```@docs
filter_by_quality
trim_low_quality_ends
remove_low_quality_reads
adaptive_quality_filtering
```

#### Example: Quality Filtering
```julia
# Filter reads by quality thresholds
filtered_reads = filter_by_quality(
    "reads.fastq",
    min_mean_quality=25,
    min_length=1000,
    max_n_percent=5
)

# Save filtered reads
write_fastq("filtered_reads.fastq", filtered_reads)

# Report filtering statistics
println("Original reads: $(count_reads("reads.fastq"))")
println("Filtered reads: $(length(filtered_reads))")
println("Retention rate: $(length(filtered_reads) / count_reads("reads.fastq") * 100)%")
```

### Adapter and Contamination Removal

```@docs
remove_adapters
detect_adapter_contamination
remove_host_contamination
remove_vector_contamination
```

#### Example: Adapter Removal
```julia
# Detect and remove adapters
adapter_results = detect_adapter_contamination("reads.fastq")

if adapter_results.contamination_detected
    println("Adapter contamination detected:")
    println("  Adapter type: $(adapter_results.adapter_type)")
    println("  Contamination rate: $(adapter_results.contamination_rate)%")
    
    # Remove adapters
    clean_reads = remove_adapters(
        "reads.fastq",
        adapter_sequences=adapter_results.adapter_sequences,
        min_overlap=10
    )
    
    write_fastq("adapter_trimmed.fastq", clean_reads)
end
```

### Length-Based Filtering

```@docs
filter_by_length
trim_to_length
remove_short_reads
normalize_read_lengths
```

#### Example: Length Filtering
```julia
# Filter reads by length criteria
length_filtered = filter_by_length(
    "reads.fastq",
    min_length=1000,
    max_length=50000
)

# Normalize length distribution (optional)
normalized_reads = normalize_read_lengths(
    length_filtered,
    target_length=15000,
    tolerance=0.2
)
```

## Contamination Detection

### Host Contamination

```@docs
detect_host_contamination
remove_host_sequences
classify_contamination_sources
```

#### Example: Host Contamination Removal
```julia
# Screen for host contamination
contamination_results = detect_host_contamination(
    "reads.fastq",
    host_genome="human_genome.fasta",
    min_identity=0.9
)

println("Host contamination: $(contamination_results.contamination_rate)%")

# Remove contaminated reads
clean_reads = remove_host_sequences(
    "reads.fastq",
    contamination_results.contaminated_reads
)
```

### Vector and Adapter Contamination

```@docs
screen_vector_contamination
detect_primer_contamination
remove_synthetic_sequences
```

#### Example: Vector Screening
```julia
# Screen for vector contamination
vector_results = screen_vector_contamination(
    "reads.fastq",
    vector_database="vector_db.fasta"
)

if vector_results.contamination_found
    println("Vector contamination detected:")
    for hit in vector_results.contaminated_sequences
        println("  $(hit.read_id): $(hit.vector_name)")
    end
end
```

## Quality Metrics and Reporting

### Standard Quality Metrics

```@docs
calculate_phred_scores
assess_base_call_accuracy
calculate_error_rates
estimate_sequencing_quality
```

#### Example: Quality Metrics Calculation
```julia
# Calculate comprehensive quality metrics
quality_metrics = calculate_comprehensive_metrics("reads.fastq")

# Phred score analysis
phred_analysis = calculate_phred_scores(quality_metrics)
println("Mean Phred score: $(phred_analysis.mean_phred)")
println("Q30+ rate: $(phred_analysis.q30_rate)%")

# Error rate estimation
error_rates = calculate_error_rates(quality_metrics)
println("Estimated error rate: $(error_rates.overall_error_rate)")
```

### Quality Control Reports

```@docs
generate_quality_report
create_quality_dashboard
export_quality_metrics
```

#### Example: Quality Report Generation
```julia
# Generate comprehensive quality report
quality_report = generate_quality_report(
    "reads.fastq",
    output_format="html",
    include_plots=true
)

# Save report
save_quality_report(quality_report, "quality_report.html")

# Export metrics for further analysis
metrics_data = export_quality_metrics(quality_report, format="csv")
```

## Specialized Quality Control

### Platform-Specific QC

```@docs
assess_hifi_quality
assess_nanopore_quality
assess_illumina_quality
```

#### Example: HiFi-Specific Quality Control
```julia
# HiFi-specific quality assessment
hifi_qc = assess_hifi_quality("hifi_reads.fastq")

println("HiFi Quality Assessment:")
println("  Mean accuracy: $(hifi_qc.mean_accuracy)")
println("  Mean length: $(hifi_qc.mean_length)")
println("  Length N50: $(hifi_qc.length_n50)")
println("  Quality distribution:")
println("    Q20+: $(hifi_qc.q20_percent)%")
println("    Q30+: $(hifi_qc.q30_percent)%")
println("    Q40+: $(hifi_qc.q40_percent)%")
```

### Application-Specific QC

```@docs
assess_assembly_readiness
assess_annotation_readiness
assess_variant_calling_readiness
```

#### Example: Assembly Readiness Assessment
```julia
# Check if reads are suitable for assembly
assembly_readiness = assess_assembly_readiness("reads.fastq")

println("Assembly Readiness:")
println("  Overall score: $(assembly_readiness.overall_score)/10")
println("  Coverage estimate: $(assembly_readiness.estimated_coverage)x")
println("  Quality sufficient: $(assembly_readiness.quality_sufficient)")
println("  Length distribution: $(assembly_readiness.length_distribution_score)")

if !assembly_readiness.ready_for_assembly
    println("Issues found:")
    for issue in assembly_readiness.issues
        println("  - $(issue)")
    end
end
```

## Visualization and Plotting

### Quality Plots

```@docs
plot_quality_distribution
plot_length_distribution
plot_gc_content_distribution
plot_base_composition
```

#### Example: Quality Visualization
```julia
# Create quality visualization plots
quality_plots = create_quality_plots("reads.fastq")

# Individual plots
plot_quality_distribution(quality_plots.quality_data, 
                         title="Per-Base Quality Scores")

plot_length_distribution(quality_plots.length_data,
                        title="Read Length Distribution")

# Combined quality dashboard
quality_dashboard = plot_quality_dashboard(quality_plots)
save_plot(quality_dashboard, "quality_dashboard.png")
```

## Performance Optimization

### Memory-Efficient Processing

```@docs
stream_quality_analysis
process_in_chunks
parallel_quality_assessment
```

#### Example: Large File Processing
```julia
# Process large files efficiently
large_file_qc = stream_quality_analysis(
    "large_reads.fastq",
    chunk_size=10000,
    parallel=true,
    threads=8
)

# Memory-efficient filtering
filtered_output = process_in_chunks(
    "large_reads.fastq",
    "filtered_reads.fastq",
    chunk_size=50000,
    filter_function=quality_filter_function
)
```

## Common Issues and Solutions

### Low Quality Data
```julia
# Identify quality issues
quality_issues = diagnose_quality_issues("reads.fastq")

for issue in quality_issues
    println("Issue: $(issue.type)")
    println("  Description: $(issue.description)")
    println("  Severity: $(issue.severity)")
    println("  Recommendation: $(issue.recommendation)")
end
```

### Contamination Problems
```julia
# Comprehensive contamination screening
contamination_screen = comprehensive_contamination_screening(
    "reads.fastq",
    databases=["host", "vector", "adapter", "primer"]
)

# Generate contamination report
contamination_report = generate_contamination_report(contamination_screen)
```

## Related Functions

### Data Input/Output
- [`read_fastq`](@ref) - Read FASTQ files
- [`write_fastq`](@ref) - Write filtered FASTQ files
- [`compress_fastq`](@ref) - Compress output files

### Downstream Analysis
- [`count_kmers`](@ref) - K-mer analysis of cleaned reads
- [`assemble_genome`](@ref) - Genome assembly with quality-controlled reads

### Visualization
- [`plot_quality_metrics`](@ref) - Quality visualization
- [`create_quality_dashboard`](@ref) - Interactive quality dashboard

## Related Workflows

### Previous Steps
- [Data Acquisition](data-acquisition.md) - Obtaining sequencing data

### Next Steps
- [Sequence Analysis](sequence-analysis.md) - K-mer analysis of quality-controlled reads
- [Genome Assembly](genome-assembly.md) - Assembly with preprocessed reads

## See Also
- [Tutorial 2: Quality Control](../../tutorials/02_quality_control.md)
- [FASTA/FASTQ Data Types](../data-types/fasta-fastq.md)
- [Performance Optimization Guide](../examples/advanced-usage.md)