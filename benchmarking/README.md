# Mycelia Benchmarking Suite

This directory contains performance benchmarking tests for Mycelia using realistic datasets. These tests are designed to:

- Evaluate computational performance and memory usage
- Test scalability with large datasets
- Validate accuracy against reference implementations
- Identify performance bottlenecks and optimization opportunities

## ⚠️ Important Notes

**These benchmarks are resource-intensive and should NOT be run in CI/CD pipelines.**

- Use realistic dataset sizes (GBs of data)
- Require significant computational time (hours to days)
- Consume substantial memory and storage
- Designed for manual execution on HPC systems

## Benchmark Categories

### 1. Data Processing Benchmarks
- Large-scale FASTQ processing
- Massive FASTA file handling
- Format conversion performance

### 2. K-mer Analysis Benchmarks
- Scalability testing with varying k-mer sizes
- Memory efficiency with sparse vs dense matrices
- Performance on different genome sizes

### 3. Assembly Benchmarks
- HiFi assembly with varying coverage depths
- Different genome sizes and complexities
- Assembly quality vs computational cost

### 4. Annotation Benchmarks
- Gene prediction on large genomes
- Parallel processing efficiency
- Accuracy vs reference annotations

### 5. Comparative Genomics Benchmarks
- Pangenome construction with many genomes
- Phylogenetic tree construction
- Large-scale comparative analysis

## Usage

### Prerequisites

```bash
# Ensure adequate resources
# Recommended: 64GB+ RAM, 1TB+ storage, 16+ cores
# Install required external tools via bioconda
mamba env create -f environment.yml
```

### Running Individual Benchmarks

```bash
# Run specific benchmark
julia --project=. benchmarking/01_data_processing_benchmark.jl

# Run with performance monitoring
julia --project=. --track-allocation=user benchmarking/01_data_processing_benchmark.jl
```

### Running All Benchmarks (HPC)

```bash
# Submit to SLURM (recommended)
sbatch benchmarking/run_all_benchmarks.sh

# Or run locally (use with extreme caution)
julia --project=. benchmarking/run_all_benchmarks.jl
```

## Benchmark Results

Results are automatically saved to:
- `results/performance_metrics.json` - Timing and memory usage
- `results/accuracy_metrics.json` - Quality and accuracy measurements
- `results/benchmark_report.html` - Comprehensive HTML report

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

All benchmarks include SLURM job scripts for HPC execution:

```bash
# Example SLURM submission
sbatch --partition=compute --time=24:00:00 --mem=64G --cpus-per-task=16 benchmarking/assembly_benchmark.sh
```

### Resource Monitoring

Benchmarks automatically collect:
- CPU usage and time
- Memory consumption (peak and average)
- Disk I/O statistics
- Network usage (if applicable)

## Benchmark Validation

### Accuracy Testing
- Compare results against known references
- Validate against published benchmarks
- Cross-validate with other tools

### Performance Regression Testing
- Track performance over time
- Identify performance regressions
- Validate optimizations

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
- **Throughput**: Data processed per unit time
- **Memory Efficiency**: Peak memory usage relative to data size
- **Scalability**: Performance vs dataset size relationship
- **Accuracy**: Quality of results vs reference standards

### Comparative Analysis
- Compare with other bioinformatics tools
- Analyze performance across different hardware
- Evaluate trade-offs between speed and accuracy

---

**Remember**: These benchmarks are designed for performance evaluation, not routine analysis. Use the regular tutorials for learning and the main package for production analysis.