# # Tutorial 3: K-mer Analysis and Feature Extraction
#
# This tutorial explores k-mer analysis, a fundamental technique in bioinformatics
# for sequence analysis, genome assembly, and comparative genomics.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - K-mer theory and biological significance
# - How k-mer size affects analysis sensitivity and specificity
# - Dense vs sparse k-mer counting strategies
# - K-mer frequency spectra and their interpretation
# - Applications in genome size estimation and quality assessment
# - Memory and computational considerations for large-scale analysis

# ## Setup
import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import FASTX
import Random
import Plots
import Statistics

Random.seed!(42)

# ## Part 1: K-mer Theory and Biological Context
#
# K-mers are subsequences of length k extracted from DNA sequences.
# They capture local sequence composition and are fundamental to many algorithms.

println("=== K-mer Analysis Tutorial ===")

# ### K-mer Mathematics
#
# For DNA sequences, there are 4^k possible k-mers.
# Understanding k-mer space helps with parameter selection.

function kmer_space_size(k, alphabet_size=4)
    return alphabet_size^k
end

println("K-mer space sizes:")
for k in [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    space_size = kmer_space_size(k)
    println("k=$k: $(space_size) possible k-mers ($(space_size ÷ 1000)K)")
end

# ### Biological Significance
#
# K-mers capture different biological features depending on their size:
# - Small k-mers (k=3-7): Capture short motifs, sensitive to errors
# - Medium k-mers (k=15-21): Balance sensitivity and specificity
# - Large k-mers (k=25-51): Specific but may miss short overlaps

println("\nK-mer size selection guidelines:")
println("k=3-7:   Short motifs, codon analysis")
println("k=15-21: Error correction, initial assembly")
println("k=25-31: Genome assembly, repeat detection")
println("k=35-51: Specific overlaps, large genome assembly")

# ## Part 2: K-mer Counting Strategies
#
# Different applications require different counting approaches.
# Understanding trade-offs helps optimize performance.

println("\n=== K-mer Counting Strategies ===")

# Generate test sequences for demonstration
test_sequences = [
    Mycelia.random_fasta_record(moltype=:DNA, seed=i, L=1000) 
    for i in 1:10
]

# Write sequences to temporary files
temp_files = String[]
for (i, seq) in enumerate(test_sequences)
    filename = "test_seq_$i.fasta"
    Mycelia.write_fasta(outfile=filename, records=[seq])
    push!(temp_files, filename)
end

println("Generated $(length(temp_files)) test sequences")

# ### Dense K-mer Counting
#
# Dense counting stores all possible k-mers, including those not observed.
# Memory usage: O(4^k) - grows exponentially with k

println("\n--- Dense K-mer Counting ---")

# TODO: Implement comprehensive dense k-mer analysis
# - Count all possible k-mers
# - Analyze frequency distributions
# - Memory usage profiling
# - Comparison across different k values

for k in [3, 5, 7, 9]
    println("Computing dense k-mer counts for k=$k...")
    
    # Memory estimation
    memory_mb = (4^k * 4) / (1024^2)  # Assuming 4 bytes per count
    println("  Estimated memory: $(round(memory_mb, digits=2)) MB")
    
    if memory_mb < 100  # Only run if memory usage is reasonable
        dense_counts = Mycelia.fasta_list_to_dense_kmer_counts(
            fasta_list=temp_files, 
            alphabet=:DNA, 
            k=k
        )
        println("  ✓ Dense counting completed")
        println("  Matrix size: $(size(dense_counts))")
    else
        println("  ⚠ Skipping due to high memory usage")
    end
end

# ### Sparse K-mer Counting
#
# Sparse counting only stores observed k-mers.
# Memory usage: O(n) where n is number of unique k-mers

println("\n--- Sparse K-mer Counting ---")

# TODO: Implement comprehensive sparse k-mer analysis
# - Count only observed k-mers
# - Analyze sparsity patterns
# - Memory efficiency comparison
# - Scalability testing

for k in [11, 13, 15, 17, 19, 21]
    println("Computing sparse k-mer counts for k=$k...")
    
    sparse_counts = Mycelia.fasta_list_to_sparse_kmer_counts(
        fasta_list=temp_files,
        alphabet=:DNA,
        k=k
    )
    println("  ✓ Sparse counting completed")
    println("  Unique k-mers: $(length(sparse_counts))")
end

# ## Part 3: K-mer Frequency Spectra
#
# K-mer frequency spectra reveal genome characteristics and data quality

println("\n=== K-mer Frequency Spectra ===")

# TODO: Implement k-mer spectrum analysis
# - Generate frequency histograms
# - Identify characteristic peaks
# - Estimate genome size and coverage
# - Detect contamination and repeats

# ## Part 4: Applications in Genome Analysis
#
# K-mers have many applications in genomic analysis

println("\n=== Genome Analysis Applications ===")

# ### Genome Size Estimation
#
# Use k-mer frequency spectra to estimate genome size
# Formula: Genome size ≈ Total k-mers / Coverage peak

println("--- Genome Size Estimation ---")

# TODO: Implement genome size estimation
# - Find coverage peak in frequency spectrum
# - Calculate total k-mers
# - Estimate genome size
# - Validate against known genome size

# ### Error Detection and Correction
#
# Low-frequency k-mers often represent sequencing errors

println("--- Error Detection ---")

# TODO: Implement error detection
# - Identify low-frequency k-mers
# - Classify as errors vs rare variants
# - Implement error correction algorithms

# ### Contamination Detection
#
# Foreign DNA creates distinctive k-mer patterns

println("--- Contamination Detection ---")

# TODO: Implement contamination detection
# - Compare k-mer profiles
# - Identify foreign k-mers
# - Quantify contamination levels

# ## Part 5: Performance Optimization
#
# Large-scale k-mer analysis requires optimization

println("\n=== Performance Optimization ===")

# ### Memory Management
#
# Strategies for handling large k-mer datasets

println("--- Memory Management ---")

# TODO: Implement memory optimization
# - Streaming k-mer processing
# - Disk-based storage
# - Memory-mapped files
# - Compression techniques

# ### Parallel Processing
#
# Accelerate k-mer counting with parallelization

println("--- Parallel Processing ---")

# TODO: Implement parallel k-mer counting
# - Multi-threaded counting
# - Distributed processing
# - Load balancing strategies

# ## Part 6: Visualization and Interpretation
#
# Create plots to understand k-mer analysis results

println("\n=== K-mer Visualization ===")

# TODO: Implement k-mer visualization
# - Frequency spectrum plots
# - K-mer composition heatmaps
# - Coverage distribution plots
# - Comparative k-mer analysis

# ## Part 7: Advanced K-mer Techniques
#
# Explore advanced k-mer analysis methods

println("\n=== Advanced Techniques ===")

# ### Minimizers
#
# Reduce k-mer space using minimizer techniques

println("--- Minimizers ---")

# TODO: Implement minimizer techniques
# - Canonical minimizers
# - Syncmers
# - Strobemers

# ### Graph Construction
#
# Build graphs from k-mer overlaps

println("--- Graph Construction ---")

# TODO: Implement k-mer graph construction
# - de Bruijn graphs
# - String graphs
# - Overlap graphs

# ## Summary and Best Practices
println("\n=== K-mer Analysis Summary ===")
println("✓ Understanding k-mer theory and biological significance")
println("✓ Choosing appropriate k-mer sizes for different applications")
println("✓ Implementing dense and sparse counting strategies")
println("✓ Analyzing k-mer frequency spectra")
println("✓ Applying k-mer analysis to genome size estimation")
println("✓ Optimizing performance for large-scale analysis")

println("\nBest Practices:")
println("- Start with k=21 for general analysis")
println("- Use dense counting for small k, sparse for large k")
println("- Monitor memory usage and optimize accordingly")
println("- Validate results with known datasets")
println("- Consider biological context in interpretation")

# Cleanup
for file in temp_files
    rm(file, force=true)
end

println("\nNext: Tutorial 4 - Genome Assembly")

nothing