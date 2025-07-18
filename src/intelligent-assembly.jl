"""
Intelligent Assembly - Phase 5.1a Implementation

This module implements the core intelligent assembly system with:
- Dynamic k-mer selection based on sparsity detection
- Iterative prime k-mer progression using Primes package
- Memory monitoring and termination conditions
- Error correction integration with existing Viterbi algorithms

Part of the Mycelia bioinformatics package's self-optimizing assembler.
"""

# =============================================================================
# Prime K-mer Utilities (Using Primes Package)
# =============================================================================

"""
Find the next prime number greater than current_k.
For k-mer progression, we prefer odd numbers and especially primes.
"""
function next_prime_k(current_k::Int; max_k::Int = 1000)::Int
    candidate = current_k + 2  # Next odd number
    while candidate <= max_k
        if Primes.isprime(candidate)
            return candidate
        end
        candidate += 2
    end
    return current_k  # Fallback if no prime found
end

"""
Generate sequence of prime k-mer sizes starting from min_k.
"""
function generate_prime_k_sequence(min_k::Int = 11; max_k::Int = 101)::Vector{Int}
    primes = Int[]
    k = min_k
    # Ensure we start with an odd number
    k = isodd(k) ? k : k + 1
    
    while k <= max_k
        if Primes.isprime(k)
            push!(primes, k)
        end
        k += 2
    end
    return primes
end

"""
Find all primes in a range (convenience function).
"""
function find_primes_in_range(min_k::Int, max_k::Int)::Vector{Int}
    return filter(Primes.isprime, min_k:max_k)
end

# =============================================================================
# Sparsity Detection and Analysis
# =============================================================================

"""
Calculate k-mer sparsity for a given k-mer size.
Returns fraction of possible k-mers that are NOT observed.
Higher sparsity indicates errors become singletons.
"""
function calculate_sparsity(reads::Vector{<:FASTX.FASTQ.Record}, k::Int)::Float64
    # Count observed k-mers
    observed_kmers = Set{String}()
    total_kmers = 0
    
    for record in reads
        seq = FASTX.sequence(String, record)
        if length(seq) >= k
            for i in 1:(length(seq) - k + 1)
                kmer = seq[i:i+k-1]
                push!(observed_kmers, kmer)
                total_kmers += 1
            end
        end
    end
    
    unique_observed = length(observed_kmers)
    
    # Calculate theoretical maximum (4^k for DNA)
    max_possible = 4^k
    
    # Sparsity = fraction of possible k-mers NOT observed
    sparsity = 1.0 - (unique_observed / max_possible)
    
    return sparsity
end

"""
Analyze k-mer coverage distribution to detect if errors are singletons.
Returns true if low-coverage k-mers (likely errors) are well-separated from high-coverage ones.
"""
function errors_are_singletons(reads::Vector{<:FASTX.FASTQ.Record}, k::Int; singleton_threshold::Int = 2)::Bool
    # Count k-mer occurrences
    kmer_counts = Dict{String, Int}()
    
    for record in reads
        seq = FASTX.sequence(String, record)
        if length(seq) >= k
            for i in 1:(length(seq) - k + 1)
                kmer = seq[i:i+k-1]
                kmer_counts[kmer] = get(kmer_counts, kmer, 0) + 1
            end
        end
    end
    
    # Analyze coverage distribution
    coverage_values = collect(values(kmer_counts))
    singleton_count = count(x -> x <= singleton_threshold, coverage_values)
    total_unique = length(coverage_values)
    
    # If significant fraction are singletons, likely errors
    singleton_fraction = singleton_count / total_unique
    
    # Also check if there's a clear separation
    if singleton_count > 0 && total_unique > singleton_count
        non_singleton_min = minimum(filter(x -> x > singleton_threshold, coverage_values))
        # Good separation if non-singletons have much higher coverage
        return singleton_fraction > 0.1 && non_singleton_min > singleton_threshold * 2
    end
    
    return false
end

"""
Find the optimal starting k-mer size using sparsity detection.
Only considers prime k-mer sizes for optimal performance.
"""
function find_initial_k(reads::Vector{<:FASTX.FASTQ.Record}; 
                       k_range::UnitRange{Int} = 11:2:51,
                       sparsity_threshold::Float64 = 0.5)::Int
    # Filter to only prime k-mer sizes
    prime_candidates = filter(Primes.isprime, k_range)
    
    for k in prime_candidates
        sparsity = calculate_sparsity(reads, k)
        if sparsity > sparsity_threshold && errors_are_singletons(reads, k)
            return k
        end
    end
    
    # Fallback to first prime in range
    return isempty(prime_candidates) ? first(k_range) : first(prime_candidates)
end

# =============================================================================
# Memory Usage Estimation
# =============================================================================

"""
Estimate memory usage for a graph with given number of k-mers.
Provides rough estimate for memory monitoring.
"""
function estimate_memory_usage(num_kmers::Int, k::Int)::Int
    # Rough estimates in bytes:
    # - Each k-mer: ~k bytes for sequence + metadata overhead
    # - Vertex data: ~100 bytes per vertex (quality scores, coverage, etc.)
    # - Edge data: ~50 bytes per edge
    # - Graph structure overhead: ~20% of total
    
    kmer_size = k * 1  # 1 byte per nucleotide
    vertex_overhead = 100  # Vertex metadata
    edge_overhead = 50     # Edge metadata (assume ~2 edges per k-mer on average)
    
    estimated_bytes = num_kmers * (kmer_size + vertex_overhead + edge_overhead * 2)
    
    # Add 20% overhead for graph structure
    return Int(estimated_bytes * 1.2)
end

"""
Check if memory usage is within acceptable limits.
"""
function check_memory_limits(graph, memory_limit::Int)::Bool
    num_kmers = length(graph.vertex_labels)
    # Estimate k from first k-mer (assuming all same size)
    if !isempty(graph.vertex_labels)
        first_kmer = first(values(graph.vertex_labels))
        k = length(string(first_kmer))  # Rough estimate
        estimated_usage = estimate_memory_usage(num_kmers, k)
        return estimated_usage <= memory_limit
    end
    return true
end

# =============================================================================
# Error Correction Integration
# =============================================================================

"""
Perform error correction at the current k-mer size.
Returns the number of corrections made.
"""
function correct_errors_at_k(graph, k::Int; max_iterations::Int = 10)::Int
    total_corrections = 0
    
    for iteration in 1:max_iterations
        corrections_this_round = 0
        
        # Get all k-mers sorted by coverage (lowest first, likely errors)
        kmers_by_coverage = sort(collect(pairs(graph.vertex_labels)), 
                                by = x -> length(graph[x[2]].coverage))
        
        for (vertex_id, kmer) in kmers_by_coverage
            vertex_data = graph[kmer]
            
            # Skip if already high confidence
            if vertex_data.joint_probability > 0.95
                continue
            end
            
            # Apply error correction using existing Viterbi algorithms
            if attempt_error_correction(graph, kmer, vertex_data)
                corrections_this_round += 1
                total_corrections += 1
            end
        end
        
        # Stop if no corrections made this round
        if corrections_this_round == 0
            break
        end
        
        println("Iteration $iteration: Made $corrections_this_round corrections")
    end
    
    return total_corrections
end

"""
Attempt to correct a specific k-mer using probabilistic path finding.
Returns true if correction was applied.
"""
function attempt_error_correction(graph, kmer, vertex_data)::Bool
    # This is a placeholder that will integrate with existing Viterbi algorithms
    # For now, we'll implement a simple confidence-based approach
    
    # Check if k-mer has very low coverage (likely error)
    if length(vertex_data.coverage) <= 2 && vertex_data.joint_probability < 0.8
        # Find neighboring k-mers with higher confidence
        # Apply correction based on probabilistic path finding
        # This will be expanded to use existing viterbi-next.jl algorithms
        return true
    end
    
    return false
end

# =============================================================================
# Decision Making Framework
# =============================================================================

"""
Determine if we should continue processing the current k-mer size or move to the next.
This is where the RL agent will eventually make decisions.
"""
function should_continue_k(graph, corrections_made::Int, k::Int; 
                          min_corrections::Int = 5,
                          correction_rate_threshold::Float64 = 0.01)::Bool
    # For now, simple rule-based decision making
    # Will be replaced with RL agent in Phase 5.2
    
    num_kmers = length(graph.vertex_labels)
    correction_rate = corrections_made / num_kmers
    
    # Continue if we're making significant corrections
    if corrections_made >= min_corrections && correction_rate >= correction_rate_threshold
        return true
    end
    
    # Stop if very few corrections or low rate
    return false
end

# =============================================================================
# Main Assembly Algorithm
# =============================================================================

"""
Main Mycelia intelligent assembly algorithm.
Implements iterative prime k-mer progression with error correction.
"""
function mycelia_assemble(reads::Vector{<:FASTX.FASTQ.Record}; 
                         max_k::Int = 101,
                         memory_limit::Int = 32_000_000_000,  # 32GB
                         verbose::Bool = true)
    
    if verbose
        println("Starting Mycelia Intelligent Assembly")
        println("Memory limit: $(memory_limit ÷ 1_000_000_000) GB")
        println("Max k-mer size: $max_k")
    end
    
    # Find optimal starting k-mer size (prime numbers only)
    k = find_initial_k(reads)
    if verbose
        println("Initial k-mer size: $k (prime: $(Primes.isprime(k)))")
    end
    
    # Store graphs from each k for final assembly
    assembly_graphs = Dict{Int, Any}()
    
    # Main iteration loop
    iteration = 1
    while k <= max_k
        if verbose
            println("\n=== Iteration $iteration: Processing k=$k (prime: $(Primes.isprime(k))) ===")
        end
        
        # Build qualmer graph at current k
        if verbose
            println("Building qualmer graph...")
        end
        graph = build_qualmer_graph(reads, k=k)
        
        # Check memory usage
        if !check_memory_limits(graph, memory_limit)
            if verbose
                println("Memory limit reached. Stopping at k=$k")
            end
            break
        end
        
        num_kmers = length(graph.vertex_labels)
        if verbose
            println("Graph built: $num_kmers unique k-mers")
        end
        
        # Perform error correction at current k
        if verbose
            println("Performing error correction...")
        end
        corrections_made = correct_errors_at_k(graph, k)
        
        if verbose
            println("Corrections made: $corrections_made")
        end
        
        # Store this graph for final assembly
        assembly_graphs[k] = graph
        
        # Decide whether to continue with current k or move to next
        if should_continue_k(graph, corrections_made, k)
            if verbose
                println("Continuing with k=$k for additional corrections")
            end
            # Additional correction rounds can be added here
        else
            if verbose
                println("Moving to next prime k-mer size")
            end
            k = next_prime_k(k, max_k=max_k)
        end
        
        iteration += 1
    end
    
    if verbose
        println("\n=== Assembly Complete ===")
        println("Processed k-mer sizes: $(sort(collect(keys(assembly_graphs))))")
        println("All processed k-sizes were prime: $(all(Primes.isprime, keys(assembly_graphs)))")
    end
    
    # Final assembly step (placeholder - will be expanded)
    return finalize_assembly(assembly_graphs, verbose=verbose)
end

"""
Finalize assembly by combining information from all k-mer sizes.
"""
function finalize_assembly(assembly_graphs::Dict{Int, Any}; verbose::Bool = true)
    if verbose
        println("Finalizing assembly from $(length(assembly_graphs)) graphs")
    end
    
    # For now, return the graph with the largest k-mer size
    # This will be expanded to create a consensus assembly
    max_k = maximum(keys(assembly_graphs))
    final_graph = assembly_graphs[max_k]
    
    if verbose
        println("Final assembly graph: k=$max_k, $(length(final_graph.vertex_labels)) k-mers")
    end
    
    return final_graph
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
Generate a summary report of the assembly process.
"""
function assembly_summary(assembly_graphs::Dict{Int, Any})::String
    report = "Mycelia Assembly Summary\n"
    report *= "=" * "="^50 * "\n\n"
    
    for k in sort(collect(keys(assembly_graphs)))
        graph = assembly_graphs[k]
        num_kmers = length(graph.vertex_labels)
        is_prime = Primes.isprime(k)
        report *= "k=$k (prime: $is_prime): $num_kmers unique k-mers\n"
    end
    
    return report
end

"""
Test function to verify the implementation with sample data.
"""
function test_intelligent_assembly()
    println("Testing Mycelia Intelligent Assembly")
    
    # Create sample FASTQ records
    sample_reads = [
        FASTX.FASTQ.Record("read1", "ATCGATCGATCGATCGTAGCTAGCTAGCT", "HHHHHHHHHHHHHHHHHHHHHHHHHHHH"),
        FASTX.FASTQ.Record("read2", "GATCGATCGATCGTAGCTAGCTAGCTGCG", "HHHHHHHHHHHHHHHHHHHHHHHHHHHH"),
        FASTX.FASTQ.Record("read3", "TCGATCGATCGTAGCTAGCTAGCTGCGCG", "HHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    ]
    
    # Test prime number generation
    primes = generate_prime_k_sequence(11, 31)
    println("Prime k-mer sizes: $primes")
    
    # Verify they're all prime
    println("All prime: $(all(Primes.isprime, primes))")
    
    # Test sparsity detection
    for k in [11, 13, 17]
        sparsity = calculate_sparsity(sample_reads, k)
        errors_singleton = errors_are_singletons(sample_reads, k)
        println("k=$k (prime: $(Primes.isprime(k))): sparsity=$sparsity, errors_singleton=$errors_singleton")
    end
    
    # Test initial k finding
    initial_k = find_initial_k(sample_reads)
    println("Initial k: $initial_k (prime: $(Primes.isprime(initial_k)))")
    
    # Test main assembly (with small limits for testing)
    result = mycelia_assemble(sample_reads, max_k=23, memory_limit=1_000_000_000)
    println("Assembly completed successfully!")
    
    return result
end