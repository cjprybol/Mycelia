# Mycelia Benchmarking Suite

This directory contains a comprehensive performance benchmarking infrastructure for Mycelia using realistic datasets and production-ready workflows. The benchmarking suite is designed to:

- **Evaluate computational performance and memory usage** with detailed profiling
- **Test scalability** with configurable dataset sizes and parallel processing
- **Track performance regressions** against baseline results
- **Validate accuracy** against reference implementations and known benchmarks
- **Identify performance bottlenecks** and optimization opportunities
- **Support HPC environments** with SLURM integration and automated resource monitoring

## ⚠️ Important Notes

**These benchmarks are resource-intensive and should NOT be run in CI/CD pipelines.**

- Use realistic dataset sizes (GBs of data)
- Require significant computational time (hours to days)
- Consume substantial memory and storage
- Designed for manual execution on HPC systems

## Benchmark Categories

### 1. Data Processing Benchmarks (`01_data_processing_benchmark.jl`)
- **FASTQ/FASTA Processing**: Large-scale sequence file I/O performance
- **Format Conversion**: Efficiency of format transformations
- **Memory Usage**: Scaling with dataset size and parallel processing
- **Parallel Processing**: Thread scalability and load balancing
- **I/O Performance**: Sequential and concurrent file operations
- **Quality Control**: Duplication rate assessment and read length determination

### 2. K-mer Analysis Benchmarks (`02_kmer_analysis_benchmark.jl`)
- **K-mer Counting**: Dense and sparse counting performance across k values
- **Canonical K-mers**: Memory efficiency with canonicalization
- **Spectrum Analysis**: Frequency histogram generation and genome size estimation
- **Memory Efficiency**: Peak usage, fragmentation, and cache performance
- **Scalability**: Performance vs genome size with efficiency metrics
- **Accuracy**: Genome size estimation validation

### 3. Assembly Benchmarks (`03_assembly_benchmark.jl`)
- **MEGAHIT Performance**: Metagenomic assembly with realistic datasets
- **metaSPAdes Testing**: Comparative assembler performance
- **Quality Metrics**: N50, contiguity, and length recovery assessment
- **Simulated Reads**: Paired-end read generation with error models
- **Resource Usage**: Memory and time scaling with coverage and genome size
- **Thread Scalability**: Multi-core assembly performance

### 4. Annotation Benchmarks (`04_annotation_benchmark.jl`)
- **Pyrodigal Gene Prediction**: Ab initio gene finding performance
- **Parallel Annotation**: Multi-genome processing with threading
- **Accuracy Assessment**: Gene prediction sensitivity and specificity
- **Organism Types**: Bacterial, fungal, and plant genome annotation
- **Throughput Analysis**: Base pairs processed per second
- **Memory Profiling**: Peak usage during annotation pipelines

### 5. Comparative Genomics Benchmarks (Planned)
- Pangenome construction with many genomes
- Phylogenetic tree construction
- Large-scale comparative analysis

### 6. Coverage Profiling Benchmarks (`06_coverm_benchmark.jl`)
- CoverM contig/genome modes on multiple BAMs
- Thread scaling for coverage and abundance summarization
- End-to-end timing from BAM generation through CoverM

### 7. Binning Benchmarks (`07_binning_benchmark.jl`)
- Opt-in timing for binning/post-binning wrappers
- Uses user-supplied contigs/depth/coverage inputs
- Records wall-clock durations for tool runs

## Usage

### Prerequisites

```bash
# Ensure adequate resources
# Recommended: 64GB+ RAM, 1TB+ storage, 16+ cores
# Install required external tools via bioconda
mamba env create -f environment.yml

# Ensure BenchmarkTools.jl is installed
julia --project=. -e "import Pkg; Pkg.add(\"BenchmarkTools\")"
```

### Benchmark Configuration

The benchmarking suite supports three scales via environment variables:

```bash
# Small scale (development/testing) - Default
export BENCHMARK_SCALE=small

# Medium scale (realistic datasets)
export BENCHMARK_SCALE=medium  

# Large scale (scalability testing)
export BENCHMARK_SCALE=large
```

### Running Individual Benchmarks

```bash
# Run specific benchmark with small scale (default)
julia --project=. benchmarking/01_data_processing_benchmark.jl

# Run with medium scale
BENCHMARK_SCALE=medium julia --project=. benchmarking/02_kmer_analysis_benchmark.jl

# Run with performance monitoring and allocation tracking
julia --project=. --track-allocation=user benchmarking/03_assembly_benchmark.jl

# Run annotation benchmark
julia --project=. benchmarking/04_annotation_benchmark.jl
```

### Running Complete Benchmark Suite

```bash
# Run orchestrated benchmark suite with regression checking
julia --project=. benchmarking/benchmark_runner.jl small

# Run medium scale benchmark suite
julia --project=. benchmarking/benchmark_runner.jl medium

# Submit to SLURM cluster (recommended for large scale)
sbatch benchmarking/run_all_benchmarks.sh
```

### Performance Regression Testing

```bash
# Run with regression checking (default)
julia --project=. benchmarking/benchmark_runner.jl small

# Update baselines for future comparisons
UPDATE_BASELINES=true julia --project=. benchmarking/benchmark_runner.jl medium

# Disable regression checking
CHECK_REGRESSION=false julia --project=. benchmarking/benchmark_runner.jl small
```

## Benchmark Results

The benchmarking suite generates comprehensive results in multiple formats:

### Result Files
- **Individual Results**: `results/[benchmark]_[timestamp].json` - Detailed benchmark data
- **Comprehensive Results**: `results/[benchmark]_comprehensive_[timestamp].json` - Enhanced metrics
- **HTML Reports**: `results/benchmark_report_[scale]_[timestamp].html` - Interactive reports
- **Baseline Storage**: `baselines/[benchmark]_baseline.json` - Reference performance data

### Performance Metrics Collected
- **Timing**: Median, mean, min, max execution times
- **Memory**: Peak usage, allocations, garbage collection statistics
- **Throughput**: Data processed per unit time (reads/min, bp/sec, etc.)
- **Scalability**: Performance efficiency vs dataset size
- **Quality**: Assembly N50, gene prediction sensitivity, etc.
- **Resource Utilization**: CPU usage, I/O patterns, memory fragmentation

### Automated Analysis Features
- **Performance Regression Detection**: Compare against baselines with configurable thresholds
- **Scaling Analysis**: Efficiency calculations across different input sizes
- **Quality Assessment**: Accuracy metrics for assembly and annotation results
- **System Information**: Hardware details, Julia version, platform data

## Dataset Requirements

### Small Scale (Development Testing)
- 1-10 MB datasets
- Quick validation of benchmark framework
- Run time: minutes

### Medium Scale (Realistic Testing)
- 100 MB - 1 GB datasets
- Representative of typical research projects
- Run time: hours

### Large Scale (Scalability Testing)
- 10+ GB datasets
- Stress testing and scalability evaluation
- Run time: days

## HPC Integration

### SLURM Job Scripts

The benchmarking suite includes production-ready SLURM integration:

```bash
# Submit complete benchmark suite to cluster
sbatch benchmarking/run_all_benchmarks.sh

# Customize resources
sbatch --partition=compute --time=24:00:00 --mem=64G --cpus-per-task=16 benchmarking/run_all_benchmarks.sh

# Set benchmark scale for cluster execution
BENCHMARK_SCALE=large sbatch benchmarking/run_all_benchmarks.sh
```

### Automated Resource Monitoring

The SLURM integration includes:
- **Timeout Management**: 4-20 hour timeouts with graceful handling
- **Resource Tracking**: Peak memory, CPU time, wall time via `sacct`
- **Output Compression**: Automatic compression of large log files
- **Error Handling**: Robust error reporting and recovery
- **Environment Setup**: Automatic module loading and conda environment activation
- **Cleanup**: Automated temporary file management

### HPC Features
- **Parallel Execution**: Coordinated benchmark runner handles multiple benchmarks
- **Performance Analysis**: Built-in regression checking on cluster results  
- **Report Generation**: Automatic HTML report creation with cluster metadata
- **Baseline Management**: Distributed baseline storage and comparison

## Benchmark Infrastructure

### Core Components

The benchmarking suite is built on several key components:

- **`benchmark_utils.jl`**: Core utilities with BenchmarkTools.jl integration
  - Memory profiling with detailed allocation tracking
  - Scaling benchmark runner for input size analysis  
  - Performance regression checking against baselines
  - Comprehensive result serialization and HTML report generation

- **`benchmark_runner.jl`**: Orchestrated execution framework
  - Coordinated multi-benchmark execution
  - Automated baseline management and regression detection
  - HTML report generation with performance summaries
  - Error handling and graceful failure recovery

- **`run_all_benchmarks.sh`**: Production SLURM integration
  - HPC resource management and monitoring
  - Automated environment setup and cleanup
  - Timeout handling and resource usage reporting

### Validation Framework

#### Accuracy Testing
- **Assembly Quality**: N50, contiguity, length recovery vs reference genomes
- **Gene Prediction**: Sensitivity assessment against simulated gene positions
- **K-mer Analysis**: Genome size estimation accuracy validation
- **Cross-validation**: Compare results with established bioinformatics tools

#### Performance Regression Testing
- **Automated Baseline Comparison**: Configurable performance thresholds (default 10%)
- **Historical Tracking**: Performance metrics stored over time
- **Regression Alerts**: Automatic detection of performance degradations
- **Improvement Detection**: Recognition of performance gains

## Contributing

When adding new benchmarks:

1. **Document Resource Requirements**: Specify minimum system requirements
2. **Include Multiple Scales**: Small, medium, and large dataset variants
3. **Add Reference Results**: Include expected performance metrics
4. **Provide Cleanup Scripts**: Ensure proper resource cleanup
5. **Test on HPC**: Validate SLURM integration

## Safety Guidelines

### Resource Management
- Always specify resource limits in job scripts
- Implement timeout mechanisms
- Monitor disk space usage
- Clean up temporary files

### Data Management
- Use temporary directories for intermediate files
- Implement checkpointing for long-running jobs
- Compress large output files
- Archive results appropriately

## Troubleshooting

### Common Issues
- **Out of Memory**: Increase memory allocation or use disk-based processing
- **Timeout**: Increase job time limits or optimize algorithms
- **Disk Space**: Clean up temporary files or use scratch storage
- **Network Issues**: Implement retry logic for data downloads

### Performance Optimization
- Profile code to identify bottlenecks
- Use parallel processing where appropriate
- Optimize memory usage patterns
- Consider algorithmic improvements

## Results Interpretation

### Performance Metrics
- **Throughput**: Data processed per unit time (reads/min, bp/sec, k-mers/sec)
- **Memory Efficiency**: Peak memory usage relative to data size and scaling patterns
- **Scalability**: Performance vs dataset size relationship with efficiency scores
- **Accuracy**: Quality of results vs reference standards (N50, sensitivity, etc.)
- **Resource Utilization**: CPU usage patterns, I/O efficiency, memory fragmentation

### Comparative Analysis
- **Tool Comparison**: Performance vs MEGAHIT, metaSPAdes, Pyrodigal
- **Hardware Analysis**: Performance across different systems and configurations  
- **Trade-off Evaluation**: Speed vs accuracy vs memory usage relationships
- **Optimization Guidance**: Bottleneck identification and improvement recommendations

### Example Metrics
```
K-mer Analysis Performance Summary:
- Sequences tested: 15
- Total sequence length: 7.5 Mbp
- K-mer throughput estimates:
  - k=15: 2,450,000 k-mers/second
  - k=21: 1,890,000 k-mers/second
- Scaling efficiency:
  - 500000bp: 0.95 (near-linear scaling)
  - 1000000bp: 0.91 (slight overhead)
- Peak memory usage: 245.67 MB
```

## Quick Start

```bash
# 1. Install dependencies
julia --project=. -e "import Pkg; Pkg.add(\"BenchmarkTools\")"

# 2. Run a quick benchmark test
julia --project=. benchmarking/01_data_processing_benchmark.jl

# 3. Run full benchmark suite with regression checking
julia --project=. benchmarking/benchmark_runner.jl small

# 4. View results
firefox results/benchmark_report_small_[timestamp].html
```

---

**Remember**: These benchmarks are designed for performance evaluation and optimization guidance, not routine analysis. Use the regular tutorials for learning and the main package for production analysis.
