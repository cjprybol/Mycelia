# Parameter Guide

Common parameters used across Mycelia functions with explanations, defaults, and best practices.

## Overview

This guide covers the most frequently used parameters in Mycelia functions. Understanding these parameters helps you:

- **Optimize performance** for your specific use case
- **Choose appropriate values** for different data types
- **Understand parameter interactions** and trade-offs
- **Troubleshoot issues** related to parameter selection

## File Path Parameters

### Input Files
Most functions accept file paths as input. Common patterns:

```julia
# Single file input
result = Mycelia.analyze_function("input.fastq")

# Multiple file input
result = Mycelia.analyze_function(["file1.fastq", "file2.fastq"])

# Directory input (processes all files)
result = Mycelia.analyze_function("input_directory/")
```

**Common Parameters:**
- `input_file::String` - Path to input file
- `input_files::Vector{String}` - Multiple input files
- `input_dir::String` - Directory containing input files

**Best Practices:**
- Use absolute paths for better reliability
- Check file existence before processing
- Handle compressed files (.gz) automatically

### Output Specification
```julia
# Output file specification
result = Mycelia.process_function("input.fastq", output="output.fastq")

# Output directory
result = Mycelia.process_function("input.fastq", output_dir="results/")

# Auto-generated output names
result = Mycelia.process_function("input.fastq", auto_output=true)
```

**Common Parameters:**
- `output::String` - Output file path
- `output_dir::String` - Output directory
- `output_prefix::String` - Prefix for output files
- `overwrite::Bool = false` - Overwrite existing files

## Quality Control Parameters

### Quality Thresholds
```julia
# Quality score thresholds
Mycelia.filter_by_quality("reads.fastq", 
    min_quality=20,        # Minimum average quality score
    min_base_quality=15,   # Minimum per-base quality
    quality_window=10      # Sliding window size
)
```

**Parameter Details:**
- `min_quality::Int = 20` - Minimum average Phred score (Q20 = 1% error)
- `min_base_quality::Int = 10` - Minimum individual base quality
- `quality_window::Int = 4` - Window size for quality assessment

**Quality Score Reference:**
- Q10: 10% error rate (poor)
- Q20: 1% error rate (acceptable)
- Q30: 0.1% error rate (good)
- Q40: 0.01% error rate (excellent)

### Length Filtering
```julia
Mycelia.filter_by_length("reads.fastq",
    min_length=1000,       # Minimum read length
    max_length=50000,      # Maximum read length
    length_tolerance=0.1   # Tolerance for length variation
)
```

**Parameter Details:**
- `min_length::Int = 500` - Minimum acceptable read length
- `max_length::Int = Inf` - Maximum acceptable read length
- `length_tolerance::Float64 = 0.2` - Acceptable length variation

**Platform-Specific Defaults:**
- **Illumina**: min_length=50, max_length=300
- **PacBio HiFi**: min_length=1000, max_length=30000
- **Nanopore**: min_length=500, max_length=100000

## K-mer Analysis Parameters

### K-mer Size Selection
```julia
Mycelia.count_kmers("reads.fastq", 
    k=21,                  # K-mer size
    alphabet=:DNA,         # Sequence alphabet
    canonical=true         # Use canonical k-mers
)
```

**Parameter Details:**
- `k::Int = 21` - K-mer size (length of subsequences)
- `alphabet::Symbol = :DNA` - Sequence alphabet (:DNA, :RNA, :PROTEIN)
- `canonical::Bool = true` - Combine forward and reverse complement

**K-mer Size Guidelines:**
- **k=11-15**: Error correction, small genomes
- **k=19-25**: General analysis, genome size estimation
- **k=31-51**: Large genomes, repeat resolution
- **k>51**: Very large genomes, high specificity

**Memory Considerations:**
- Dense counting: 4^k possible k-mers
- k=15: ~1 GB memory
- k=21: ~17 GB memory (use sparse)
- k≥25: Always use sparse counting

### Counting Methods
```julia
Mycelia.count_kmers("reads.fastq",
    method="sparse",       # Counting method
    min_count=1,          # Minimum count threshold
    max_count=1000        # Maximum count threshold
)
```

**Method Options:**
- `"dense"`: Store all possible k-mers (memory intensive)
- `"sparse"`: Store only observed k-mers (memory efficient)
- `"streaming"`: Process in chunks (very large files)

## Assembly Parameters

### Assembly Configuration
```julia
Mycelia.Rhizomorph.assemble_genome("reads.fastq",
    assembler="hifiasm",   # Assembly software
    k=31,                 # K-mer size for assembly
    min_overlap=1000,     # Minimum overlap length
    threads=8             # Number of CPU threads
)
```

**Assembler Options:**
- `"hifiasm"`: Best for HiFi reads
- `"canu"`: Good for long reads with higher error rates
- `"flye"`: Fast assembly for long reads
- `"spades"`: Best for Illumina reads

**Performance Parameters:**
- `threads::Int = 4` - Number of CPU threads
- `memory_gb::Int = 16` - Maximum memory usage
- `tmp_dir::String = "/tmp"` - Temporary file directory

### Assembly Quality Control
```julia
Mycelia.Rhizomorph.assemble_genome("reads.fastq",
    min_contig_length=1000,  # Minimum contig size
    min_coverage=5,          # Minimum coverage depth
    error_correction=true,   # Enable error correction
    polish=true             # Enable polishing
)
```

## Comparative Genomics Parameters

### Pangenome Construction
```julia
Mycelia.build_pangenome(genomes,
    similarity_threshold=0.95,  # Gene similarity cutoff
    coverage_threshold=0.8,     # Minimum coverage for alignment
    clustering_method="mcl",    # Clustering algorithm
    inflation=2.0              # MCL inflation parameter
)
```

**Similarity Thresholds:**
- `0.95`: Very strict (same species)
- `0.90`: Strict (closely related strains)
- `0.80`: Moderate (related species)
- `0.70`: Permissive (distant relationships)

### Phylogenetic Analysis
```julia
Mycelia.build_phylogenetic_tree(alignment,
    method="ml",              # Tree construction method
    model="GTR+G",           # Evolutionary model
    bootstrap=1000,          # Bootstrap replicates
    outgroup="species_A"     # Outgroup specification
)
```

**Method Options:**
- `"ml"`: Maximum likelihood (most accurate)
- `"nj"`: Neighbor-joining (fast)
- `"mp"`: Maximum parsimony (character-based)

## Performance Parameters

### Parallel Processing
```julia
Mycelia.parallel_function(data,
    threads=8,               # Number of CPU threads
    workers=4,              # Number of worker processes
    chunk_size=1000,        # Data chunk size
    load_balance=true       # Enable load balancing
)
```

**Thread Guidelines:**
- Use `Sys.CPU_THREADS` for maximum threads
- Leave 1-2 threads free for system
- Memory-bound tasks: threads = cores
- I/O-bound tasks: threads = 2x cores

### Memory Management
```julia
Mycelia.memory_intensive_function(data,
    memory_limit_gb=16,      # Maximum memory usage
    chunk_processing=true,   # Process in chunks
    gc_frequency=1000,      # Garbage collection frequency
    tmp_dir="/fast_storage" # Temporary file location
)
```

**Memory Guidelines:**
- Monitor with `memory_usage()` function
- Use streaming for files > available RAM
- Set conservative limits for shared systems

## File Format Parameters

### Compression
```julia
Mycelia.write_output(data,
    compress=true,          # Enable compression
    compression_level=6,    # Compression level (1-9)
    format="auto"          # Output format detection
)
```

**Compression Levels:**
- 1-3: Fast compression, larger files
- 4-6: Balanced compression and speed
- 7-9: Maximum compression, slower

### Format Specification
```julia
Mycelia.read_sequences("input.file",
    format="auto",          # Format detection
    validate=true,         # Validate file format
    encoding="utf-8"       # Text encoding
)
```

**Format Options:**
- `"auto"`: Automatic detection from extension
- `"fasta"`: FASTA format
- `"fastq"`: FASTQ format
- `"gff3"`: GFF3 annotation format

## Error Handling Parameters

### Validation and Checks
```julia
Mycelia.robust_function(input,
    validate_input=true,    # Validate input data
    strict_mode=false,     # Strict error checking
    continue_on_error=false, # Continue despite errors
    max_errors=10          # Maximum allowed errors
)
```

### Retry Logic
```julia
Mycelia.network_function(url,
    max_retries=3,         # Maximum retry attempts
    retry_delay=30,        # Delay between retries (seconds)
    exponential_backoff=true, # Increase delay each retry
    timeout=300           # Operation timeout (seconds)
)
```

## Common Parameter Patterns

### Quality Control Pattern
```julia
standard_qc_params = Dict(
    :min_quality => 20,
    :min_length => 1000,
    :max_n_percent => 5,
    :trim_ends => true,
    :remove_duplicates => false
)
```

### Performance Pattern
```julia
performance_params = Dict(
    :threads => min(8, Sys.CPU_THREADS),
    :memory_gb => 16,
    :chunk_size => 10000,
    :parallel => true
)
```

### Output Pattern
```julia
output_params = Dict(
    :output_dir => "results",
    :compress => true,
    :overwrite => false,
    :create_manifest => true
)
```

## Parameter Validation

### Built-in Validation
```julia
# Most functions automatically validate parameters
try
    result = Mycelia.count_kmers("reads.fastq", k=0)  # Invalid k
catch ArgumentError as e
    println("Parameter error: $e")
end
```

### Manual Validation
```julia
# Validate parameters before expensive operations
if !Mycelia.validate_parameters(k=21, min_quality=20, threads=8)
    error("Invalid parameter combination")
end
```

## Troubleshooting Common Issues

### Memory Problems
```julia
# Reduce memory usage
Mycelia.count_kmers("large_file.fastq", 
    k=21,
    method="sparse",        # Use sparse instead of dense
    chunk_size=50000,      # Process in smaller chunks
    memory_limit_gb=8      # Set memory limit
)
```

### Performance Issues
```julia
# Optimize for speed
Mycelia.process_function(data,
    threads=Sys.CPU_THREADS,  # Use all available cores
    chunk_size=1000,          # Optimize chunk size
    parallel=true,            # Enable parallelization
    cache_results=true        # Cache intermediate results
)
```

### File I/O Problems
```julia
# Handle file I/O robustly
Mycelia.read_function("file.fastq",
    validate=true,            # Validate file format
    buffer_size=8192,        # Optimize buffer size
    encoding="utf-8",        # Specify encoding
    handle_errors="skip"     # Skip problematic records
)
```

## See Also
- [Function Index](function-index.md) - Complete function listing
- [Basic Workflows](../examples/basic-workflows.md) - Parameter usage examples
- [Advanced Usage](../examples/advanced-usage.md) - Complex parameter combinations
- [Performance Guide](../examples/advanced-usage.md#performance-optimization) - Optimization strategies