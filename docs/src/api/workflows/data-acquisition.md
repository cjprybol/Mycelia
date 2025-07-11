# Data Acquisition & Simulation

Functions for downloading genomic data from public databases and simulating synthetic datasets for testing and benchmarking.

## Overview

Data acquisition is the first step in any bioinformatics analysis. Mycelia provides tools for:

- **Downloading reference genomes** from NCBI databases
- **Simulating realistic sequencing data** for testing
- **Managing data provenance** and metadata
- **Handling multiple data formats** and sources

## Common Workflows

### 1. Download Reference Genome
```julia
# Download a specific genome by accession
genome_file = download_genome_by_accession("NC_001422.1")

# Download complete assembly with annotations
assembly_data = ncbi_genome_download_accession("GCF_000819615.1")
```

### 2. Simulate Test Data
```julia
# Create synthetic genome
reference = simulate_random_genome(length=100000, gc_content=0.45)

# Generate HiFi reads
reads = simulate_hifi_reads(reference, coverage=30, error_rate=0.001)
```

### 3. Batch Download
```julia
# Download multiple genomes
accessions = ["GCF_000005825.2", "GCF_000009605.1", "GCF_000027325.1"]
genomes = download_genomes_batch(accessions, output_dir="genomes/")
```

## Public Database Downloads

### NCBI Genome Downloads

```@docs
download_genome_by_accession
ncbi_genome_download_accession
download_genomes_batch
```

#### Example: Download phiX174 Genome
```julia
# Download the classic phiX174 bacteriophage genome
phix_file = download_genome_by_accession("NC_001422.1")

# Verify download
@assert isfile(phix_file)
@assert endswith(phix_file, ".fna.gz")

# Read the genome
genome_record = first(read_fasta(phix_file))
sequence = String(FASTX.sequence(genome_record))
println("Downloaded genome: $(length(sequence)) bp")
```

#### Example: Complete Assembly Package
```julia
# Download complete E. coli assembly with all associated files
assembly = ncbi_genome_download_accession("GCF_000005825.2")

# Available files
println("Genome: $(assembly.genome)")
println("Proteins: $(assembly.protein)")
println("Annotations: $(assembly.gff3)")
println("CDS: $(assembly.cds)")
```

### SRA Data Downloads

```@docs
download_sra_data
prefetch_sra_runs
fasterq_dump_parallel
```

#### Example: Download Sequencing Reads
```julia
# Download reads from SRA
sra_run = "SRR1234567"
fastq_files = download_sra_data(sra_run, output_dir="reads/")

# Process paired-end reads
if length(fastq_files) == 2
    read1, read2 = fastq_files
    println("Downloaded paired-end reads:")
    println("  R1: $read1")
    println("  R2: $read2")
end
```

## Data Simulation

### Genome Simulation

```@docs
simulate_random_genome
simulate_genome_with_features
generate_synthetic_chromosome
```

#### Example: Create Test Genome
```julia
# Generate random genome with realistic GC content
test_genome = simulate_random_genome(
    length=50000,
    gc_content=0.42,
    seed=123  # for reproducibility
)

# Add realistic features
genome_with_genes = simulate_genome_with_features(
    test_genome,
    n_genes=50,
    gene_length_dist=(500, 2000),
    intergenic_length_dist=(100, 1000)
)
```

### Sequencing Read Simulation

```@docs
simulate_hifi_reads
simulate_illumina_reads
simulate_nanopore_reads
add_sequencing_errors
```

#### Example: HiFi Read Simulation
```julia
# Simulate high-quality HiFi reads
hifi_reads = simulate_hifi_reads(
    reference_genome,
    coverage=25,
    read_length_mean=15000,
    read_length_std=3000,
    error_rate=0.001
)

# Write to FASTQ
write_fastq("hifi_reads.fastq", hifi_reads)
```

#### Example: Illumina Read Simulation
```julia
# Simulate paired-end Illumina reads
illumina_reads = simulate_illumina_reads(
    reference_genome,
    coverage=50,
    read_length=150,
    fragment_size_mean=300,
    fragment_size_std=50,
    error_rate=0.01
)

# Write paired-end files
write_fastq("illumina_R1.fastq", illumina_reads.read1)
write_fastq("illumina_R2.fastq", illumina_reads.read2)
```

## Data Validation

### Download Validation

```@docs
validate_download_integrity
check_file_format
verify_genome_completeness
```

#### Example: Validate Downloaded Data
```julia
# Check file integrity
integrity_ok = validate_download_integrity("genome.fna.gz")

# Verify file format
format_ok = check_file_format("genome.fna.gz", expected_format="fasta")

# Check genome completeness
completeness = verify_genome_completeness("genome.fna.gz")
println("Genome completeness: $(completeness.complete_sequences)/$(completeness.total_sequences)")
```

### Simulation Validation

```@docs
validate_simulation_parameters
calculate_simulation_statistics
compare_simulated_vs_real
```

#### Example: Validate Simulated Data
```julia
# Check simulation statistics
sim_stats = calculate_simulation_statistics(simulated_reads)
println("Simulated reads: $(sim_stats.n_reads)")
println("Mean length: $(sim_stats.mean_length)")
println("Coverage: $(sim_stats.coverage)x")

# Compare with real data characteristics
comparison = compare_simulated_vs_real(simulated_reads, real_reads)
println("Length distribution similarity: $(comparison.length_similarity)")
println("Quality distribution similarity: $(comparison.quality_similarity)")
```

## Metadata Management

### Data Provenance

```@docs
create_data_manifest
track_data_provenance
generate_metadata_report
```

#### Example: Track Data Sources
```julia
# Create manifest for downloaded data
manifest = create_data_manifest(
    files=["genome.fna.gz", "annotations.gff3"],
    sources=["NCBI:NC_001422.1", "NCBI:GCF_000819615.1"],
    download_date=now(),
    version="1.0"
)

# Save manifest
save_manifest(manifest, "data_manifest.json")
```

### Batch Processing

```@docs
process_accession_list
download_with_retry
parallel_download
```

#### Example: Batch Download with Error Handling
```julia
# Download multiple genomes with retry logic
accession_list = ["GCF_000005825.2", "GCF_000009605.1", "GCF_000027325.1"]

results = parallel_download(
    accession_list,
    max_retries=3,
    delay_between_retries=30,
    max_concurrent=5
)

# Check results
successful = [r for r in results if r.success]
failed = [r for r in results if !r.success]

println("Downloaded: $(length(successful))/$(length(accession_list))")
if !isempty(failed)
    println("Failed downloads: $(length(failed))")
    for failure in failed
        println("  $(failure.accession): $(failure.error)")
    end
end
```

## Performance Considerations

### Memory Usage
- **Genome downloads**: Minimal memory usage (streaming)
- **Large simulations**: Memory scales with genome size
- **Batch operations**: Consider parallel processing limits

### Network Considerations
- **NCBI rate limits**: Respect API rate limits
- **Retry logic**: Implement exponential backoff
- **Parallel downloads**: Limit concurrent connections

### Storage Requirements
- **Compressed files**: Use .gz compression by default
- **Temporary files**: Clean up intermediate files
- **Disk space**: Monitor available space for large datasets

## Common Issues and Solutions

### Download Failures
```julia
# Handle network timeouts
try
    genome = download_genome_by_accession("NC_001422.1", timeout=300)
catch NetworkError
    println("Download failed - retrying with longer timeout")
    genome = download_genome_by_accession("NC_001422.1", timeout=600)
end
```

### Simulation Validation
```julia
# Validate simulation parameters before large runs
validate_simulation_parameters(
    genome_size=1000000,
    coverage=30,
    read_length=15000
)
```

## Related Functions

### File I/O
- [`read_fasta`](@ref) - Read FASTA files
- [`write_fastq`](@ref) - Write FASTQ files
- [`compress_file`](@ref) - Compress downloaded files

### Quality Control
- [`analyze_fastq_quality`](@ref) - Assess downloaded read quality
- [`calculate_genome_stats`](@ref) - Analyze genome characteristics

### Next Steps
- [Quality Control](quality-control.md) - Assess data quality
- [Sequence Analysis](sequence-analysis.md) - Analyze sequence composition

## See Also
- [Tutorial 1: Data Acquisition](../../tutorials/01_data_acquisition.md)
- [FASTA/FASTQ Data Types](../data-types/fasta-fastq.md)
- [Basic Workflows](../examples/basic-workflows.md)