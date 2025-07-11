#!/bin/bash
#SBATCH --job-name=mycelia_benchmarks
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=benchmarks_%j.out
#SBATCH --error=benchmarks_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# Mycelia Benchmarking Suite - HPC Execution Script
# This script runs comprehensive benchmarks on HPC systems using SLURM

echo "=== Mycelia Benchmarking Suite ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo ""

# Set up environment
echo "Setting up environment..."
module purge
module load julia/1.9.0  # Adjust version as needed
module load bioconda
module load gcc/11.2.0    # For compilation if needed

# Activate conda environment with bioinformatics tools
source activate mycelia-bench

# Set Julia environment
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_PROJECT=.

# Create results directory
mkdir -p results
mkdir -p temp_data

# Set temporary directory to node-local storage (if available)
if [ -d "/tmp" ]; then
    export TMPDIR="/tmp/mycelia_$SLURM_JOB_ID"
    mkdir -p $TMPDIR
fi

echo "Environment setup complete"
echo "Julia threads: $JULIA_NUM_THREADS"
echo "Temporary directory: $TMPDIR"
echo ""

# Function to run benchmark with error handling
run_benchmark() {
    local benchmark_name=$1
    local benchmark_file=$2
    
    echo "=== Running $benchmark_name ==="
    echo "Start time: $(date)"
    
    # Run with timeout and resource monitoring
    timeout 4h julia --project=. --track-allocation=user "$benchmark_file" 2>&1 | tee "results/${benchmark_name}_output.log"
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ]; then
        echo "✓ $benchmark_name completed successfully"
    elif [ $exit_code -eq 124 ]; then
        echo "✗ $benchmark_name timed out (4 hours)"
    else
        echo "✗ $benchmark_name failed with exit code $exit_code"
    fi
    
    echo "End time: $(date)"
    echo ""
    
    return $exit_code
}

# Main benchmarking sequence
echo "=== Starting Benchmarking Sequence ==="

# Data Processing Benchmarks
if [ -f "benchmarking/01_data_processing_benchmark.jl" ]; then
    run_benchmark "data_processing" "benchmarking/01_data_processing_benchmark.jl"
fi

# K-mer Analysis Benchmarks
if [ -f "benchmarking/02_kmer_analysis_benchmark.jl" ]; then
    run_benchmark "kmer_analysis" "benchmarking/02_kmer_analysis_benchmark.jl"
fi

# Assembly Benchmarks
if [ -f "benchmarking/03_assembly_benchmark.jl" ]; then
    run_benchmark "assembly" "benchmarking/03_assembly_benchmark.jl"
fi

# Annotation Benchmarks
if [ -f "benchmarking/04_annotation_benchmark.jl" ]; then
    run_benchmark "annotation" "benchmarking/04_annotation_benchmark.jl"
fi

# Comparative Genomics Benchmarks
if [ -f "benchmarking/05_comparative_benchmark.jl" ]; then
    run_benchmark "comparative" "benchmarking/05_comparative_benchmark.jl"
fi

# Generate HTML reports from executed benchmarks
echo "=== Generating HTML Reports ==="
if command -v jupyter &> /dev/null; then
    julia --project=. run_extended_tests.jl reports
    echo "✓ HTML reports generated"
else
    echo "⚠ Jupyter not found - install with: pip install jupyter nbconvert"
    echo "  Reports can be generated later with: julia --project=. run_extended_tests.jl reports"
fi

# System resource usage summary
echo "=== Resource Usage Summary ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Peak memory usage: $(sacct -j $SLURM_JOB_ID --format=MaxRSS --noheader | head -1)"
echo "CPU time: $(sacct -j $SLURM_JOB_ID --format=CPUTime --noheader | head -1)"
echo "Wall time: $(sacct -j $SLURM_JOB_ID --format=Elapsed --noheader | head -1)"

# Cleanup temporary files
echo "=== Cleanup ==="
if [ -d "$TMPDIR" ]; then
    echo "Cleaning up temporary directory: $TMPDIR"
    rm -rf $TMPDIR
fi

# Compress large output files
echo "Compressing large output files..."
find results -name "*.log" -size +10M -exec gzip {} \;
find results -name "*.json" -size +10M -exec gzip {} \;

echo "=== Benchmarking Complete ==="
echo "End time: $(date)"
echo "Results saved to: results/"
echo "Job completed with exit code: $?"

# Exit with success
exit 0