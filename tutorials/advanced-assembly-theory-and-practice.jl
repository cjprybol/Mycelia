# # Advanced Assembly Theory and Practice in Mycelia
#
# This tutorial demonstrates the theoretical foundations and practical implementation
# of Mycelia's advanced assembly algorithms, based on research from Mycelia-Dev.
#
# ## Learning Objectives
#
# - Understand the mathematical foundations of k-mer size selection
# - Apply dynamic prime pattern algorithms for optimal k-mer progression
# - Implement strain-resolved assembly with quality metrics
# - Use statistical graph cleanup methods
# - Apply metagenomic classification with MAPQ-aware techniques

using Mycelia

# ## 1. Mathematical Foundations of K-mer Selection
#
# Mycelia's k-mer selection is based on rigorous mathematical principles derived
# from extensive research on error rate relationships and genomic properties.

## ### Error Rate-Based K-mer Selection
##
## The fundamental relationship between error rate and optimal k-mer size:
## lower_bound_k = max(3, floor(1/error_rate - 1))

function demonstrate_error_rate_kmer_selection()
    println("Error Rate-Based K-mer Size Selection")
    println("=====================================")
    
    error_rates = [0.01, 0.05, 0.10, 0.15, 0.20]
    
    for error_rate in error_rates
        lower_bound = max(3, Int(floor(1/error_rate - 1)))
        
        ## Ensure odd k-mer (better for biological sequences)
        if lower_bound % 2 == 0
            lower_bound += 1
        end
        
        println("Error rate: $(error_rate*100)% → Minimum k-mer size: $lower_bound")
    end
    
    return error_rates
end

demonstrate_error_rate_kmer_selection()

# ### Sequence Length Optimization
#
# The log₄(sequence_length) pattern provides optimal starting points for
# k-mer size selection based on the divergence point where erroneous
# k-mers begin to dominate true k-mers.

function demonstrate_log4_optimization()
    println("\nLog₄ Sequence Length Optimization")
    println("=================================")
    
    sequence_lengths = [100, 1_000, 10_000, 100_000, 1_000_000]
    
    for seq_len in sequence_lengths
        optimal_k = Int(round(log(4, seq_len)))
        
        ## Ensure odd and prime when possible
        if optimal_k % 2 == 0
            optimal_k += 1
        end
        
        if !Primes.isprime(optimal_k)
            optimal_k = Primes.nextprime(optimal_k)
        end
        
        println("Sequence length: $(seq_len) bp → Optimal starting k: $optimal_k")
    end
end

demonstrate_log4_optimization()

# ## 2. Dynamic Prime Pattern Algorithm
#
# Mycelia implements a sophisticated dynamic k-mer selection algorithm that
# exploits the mathematical properties of prime number distribution.

## Generate optimal k-mer progression using dynamic prime pattern
function demonstrate_dynamic_prime_pattern()
    println("\nDynamic Prime Pattern K-mer Selection")
    println("====================================")
    
    ## Standard progression for high-quality data
    k_sequence_standard = Mycelia.dynamic_k_prime_pattern(11, 101, 2)
    println("Standard progression (start=11): $k_sequence_standard")
    
    ## Progression for error-prone data
    k_sequence_error_prone = Mycelia.dynamic_k_prime_pattern(7, 51, 2)
    println("Error-prone progression (start=7): $k_sequence_error_prone")
    
    ## Error rate optimized progression
    k_sequence_optimized = Mycelia.error_optimized_k_sequence(0.05, 101, 10000)
    println("Error-optimized (5% error, 10kb): $k_sequence_optimized")
    
    return k_sequence_standard, k_sequence_error_prone, k_sequence_optimized
end

k_sequences = demonstrate_dynamic_prime_pattern()

# ### Theoretical Advantages of Prime K-mers
#
# 1. **Twin Prime Avoidance**: Automatically skips redundant analysis
# 2. **Progressive Spacing**: Reduces computational overlap
# 3. **Hardware Optimization**: k=31 fits optimally in 64-bit integers
# 4. **Biological Relevance**: Primes cannot form perfect repeats

println("\nAdvantages of Prime K-mer Selection:")
println("• Twin prime avoidance reduces redundant analysis")
println("• Progressive spacing minimizes computational overlap")
println("• Hardware-optimized for 64-bit architectures")
println("• Cannot form perfect repeats, improving specificity")

# ## 3. Strain-Resolved Assembly Framework
#
# Mycelia provides comprehensive data structures and algorithms for
# strain-resolved assembly with detailed quality metrics.

## Create a strain-resolved assembly result with quality metrics
function demonstrate_strain_resolved_assembly()
    println("\nStrain-Resolved Assembly Framework")
    println("=================================")
    
    ## Initialize strain quality metrics
    metrics = Mycelia.StrainQualityMetrics()
    
    ## Set example values (in practice, these would be calculated)
    metrics = Mycelia.StrainQualityMetrics(
        0.95,  # strain_recall
        0.92,  # strain_precision  
        0.0,   # strain_f1_score (will be calculated)
        50000, # NGA50
        150000, # largest_alignment
        900000, # total_aligned_length
        2.5,   # mismatches_per_100kb
        0.8,   # indels_per_100kb
        3,     # misassemblies_total
        1,     # misassemblies_local
        0.98,  # completeness
        0.02,  # contamination
        0.94,  # busco_completeness
        0.96,  # read_mapping_rate
        0.93,  # properly_paired_rate
        0.88,  # mean_base_confidence
        0.05,  # uncertain_regions_fraction
        8.5,   # peak_memory_gb
        2.3,   # cpu_time_hours
        0.97,  # homopolymer_accuracy
        0.85   # repeat_resolution_rate
    )
    
    ## Calculate F1 score
    metrics = Mycelia.update_strain_f1_score!(metrics)
    
    println("Strain Quality Metrics:")
    println("• F1 Score: $(round(metrics.strain_f1_score, digits=3))")
    println("• NGA50: $(metrics.NGA50) bp")
    println("• Completeness: $(metrics.completeness*100)%")
    println("• Contamination: $(metrics.contamination*100)%")
    
    return metrics
end

strain_metrics = demonstrate_strain_resolved_assembly()

# ## 4. Statistical Graph Cleanup Methods
#
# Demonstrates the statistical tip clipping algorithm that uses coverage-based
# statistics rather than arbitrary cutoffs.

## Simulate a simple graph cleanup scenario
function demonstrate_statistical_cleanup_theory()
    println("\nStatistical Graph Cleanup Theory")
    println("===============================")
    
    ## Simulate coverage values for a connected component
    ## (mixture of high-coverage true nodes and low-coverage error nodes)
    true_node_coverage = [95, 98, 102, 89, 97, 93, 101, 96, 99, 94]
    error_node_coverage = [1, 2, 1, 3, 1]
    
    all_coverage = vcat(true_node_coverage, error_node_coverage)
    
    ## Calculate statistical thresholds
    median_coverage = Statistics.median(all_coverage)
    coverage_std = Statistics.std(all_coverage)
    
    println("Coverage Statistics:")
    println("• Median: $(round(median_coverage, digits=1))")
    println("• Standard deviation: $(round(coverage_std, digits=1))")
    
    ## Apply 3σ rule for tip removal
    threshold_3sigma = median_coverage - 3.0 * coverage_std
    threshold_1x = 1.0
    
    println("\nCleanup Thresholds:")
    println("• 3σ threshold: $(round(threshold_3sigma, digits=1))")
    println("• Hard 1x threshold: $threshold_1x")
    
    ## Identify nodes for removal
    nodes_for_removal = []
    for (i, coverage) in enumerate(all_coverage)
        if coverage <= threshold_1x || coverage < threshold_3sigma
            push!(nodes_for_removal, i)
        end
    end
    
    println("• Nodes marked for removal: $nodes_for_removal")
    println("• True nodes preserved: $(length(true_node_coverage))")
    println("• Error nodes removed: $(length(error_node_coverage))")
    
    return threshold_3sigma, nodes_for_removal
end

cleanup_results = demonstrate_statistical_cleanup_theory()

# ## 5. Metagenomic MAPQ-Aware Classification
#
# Demonstrates the critical insight about MAPQ score interpretation in
# metagenomic contexts and the correct approach for taxonomic assignment.

function demonstrate_mapq_aware_classification()
    println("\nMAPQ-Aware Metagenomic Classification")
    println("====================================")
    
    ## Simulate alignment records with different MAPQ scores
    simulated_alignments = [
        (read_id="read_001", ref_id="species_A", alignment_score=95.5, mapq=60, taxon="Escherichia coli"),
        (read_id="read_002", ref_id="species_B", alignment_score=88.2, mapq=0, taxon="Escherichia albertii"),  # Ambiguous but informative!
        (read_id="read_003", ref_id="species_A", alignment_score=92.1, mapq=30, taxon="Escherichia coli"),
        (read_id="read_004", ref_id="species_C", alignment_score=85.7, mapq=0, taxon="Shigella flexneri"),     # Ambiguous but informative!
        (read_id="read_005", ref_id="species_A", alignment_score=97.8, mapq=60, taxon="Escherichia coli"),
    ]
    
    println("Traditional Approach (Filter MAPQ=0):")
    traditional_retained = filter(a -> a.mapq > 0, simulated_alignments)
    println("• Reads retained: $(length(traditional_retained))/$(length(simulated_alignments))")
    
    ## Calculate traditional abundances
    traditional_taxa = [a.taxon for a in traditional_retained]
    traditional_counts = Dict{String, Int}()
    for taxon in traditional_taxa
        traditional_counts[taxon] = get(traditional_counts, taxon, 0) + 1
    end
    
    println("• Traditional relative abundances:")
    total_traditional = length(traditional_retained)
    for (taxon, count) in traditional_counts
        abundance = count / total_traditional
        println("  - $taxon: $(round(abundance*100, digits=1))%")
    end
    
    println("\nMycelia Approach (Weight by Alignment Score):")
    println("• Reads retained: $(length(simulated_alignments))/$(length(simulated_alignments))")
    
    ## Calculate weighted abundances using alignment scores
    weighted_abundances = Dict{String, Float64}()
    total_weight = 0.0
    
    for alignment in simulated_alignments
        weight = alignment.alignment_score  # Use alignment score, NOT MAPQ
        taxon = alignment.taxon
        
        weighted_abundances[taxon] = get(weighted_abundances, taxon, 0.0) + weight
        total_weight += weight
    end
    
    ## Normalize to relative abundances
    println("• Weighted relative abundances:")
    for (taxon, weight) in weighted_abundances
        abundance = weight / total_weight
        println("  - $taxon: $(round(abundance*100, digits=1))%")
    end
    
    ## Show the critical insight
    mapq_zero_reads = filter(a -> a.mapq == 0, simulated_alignments)
    println("\nCritical Insight:")
    println("• MAPQ=0 reads retained: $(length(mapq_zero_reads))")
    println("• These reads provide valuable taxonomic information despite ambiguity")
    println("• Traditional filtering would lose $(length(mapq_zero_reads))/$(length(simulated_alignments)) of data")
    
    return weighted_abundances, traditional_counts
end

abundance_comparison = demonstrate_mapq_aware_classification()

# ## 6. Assembly Accuracy Assessment with Graph Metrics
#
# Demonstrates novel graph-based quality metrics using Jaccard similarity
# between k-mer sets and edge sets.

function demonstrate_graph_quality_metrics()
    println("\nGraph-Based Assembly Quality Metrics")
    println("===================================")
    
    ## Simulate k-mer sets from reference and assembly
    reference_kmers = Set(["ATCG", "TCGA", "CGAT", "GATC", "ATCG"])
    assembly_kmers = Set(["ATCG", "TCGA", "CGTT", "GATC", "ATCG"])  # One substitution error
    
    ## Simulate edge sets (k+1 mers representing transitions)
    reference_edges = Set(["ATCGA", "TCGAT", "CGATC"])
    assembly_edges = Set(["ATCGA", "TCGTT", "CGATC"])  # One edge error
    
    ## Calculate Jaccard similarities
    kmer_intersection = length(intersect(reference_kmers, assembly_kmers))
    kmer_union = length(union(reference_kmers, assembly_kmers))
    kmer_jaccard = kmer_intersection / kmer_union
    kmer_distance = 1.0 - kmer_jaccard
    
    edge_intersection = length(intersect(reference_edges, assembly_edges))
    edge_union = length(union(reference_edges, assembly_edges))
    edge_jaccard = edge_intersection / edge_union
    edge_distance = 1.0 - edge_jaccard
    
    println("Assembly Quality Assessment:")
    println("• K-mer Jaccard similarity: $(round(kmer_jaccard, digits=3))")
    println("• K-mer distance: $(round(kmer_distance, digits=3))")
    println("• Edge Jaccard similarity: $(round(edge_jaccard, digits=3))")
    println("• Edge distance: $(round(edge_distance, digits=3))")
    
    ## Interpretation
    println("\nInterpretation:")
    println("• Lower distances indicate higher assembly accuracy")
    println("• Edge distance captures connectivity errors")
    println("• K-mer distance captures sequence content errors")
    
    return kmer_distance, edge_distance
end

quality_metrics = demonstrate_graph_quality_metrics()

# ## 7. Probabilistic Assembly Framework
#
# Demonstrates the theoretical foundation of treating assembly as an
# Expectation-Maximization problem with maximum likelihood inference.

function demonstrate_probabilistic_assembly_theory()
    println("\nProbabilistic Assembly Framework")
    println("===============================")
    
    ## Simulate quality scores and error probabilities
    quality_scores = [30, 35, 25, 40, 20]  # PHRED scores
    
    println("Quality Score to Error Probability Conversion:")
    for q in quality_scores
        error_prob = 10^(q / -10.0)
        accuracy = 1.0 - error_prob
        println("• Q$q → Error rate: $(round(error_prob*100, digits=4))%, Accuracy: $(round(accuracy*100, digits=2))%")
    end
    
    ## Demonstrate likelihood calculation for sequence correction
    println("\nSequence Likelihood Calculation:")
    observed_sequence = "ATCGAT"
    candidate_sequences = ["ATCGAT", "ATCGTT", "ATCCAT"]
    
    ## Assume average quality score of 30 (0.1% error rate)
    avg_error_rate = 0.001
    
    for candidate in candidate_sequences
        matches = sum(observed_sequence[i] == candidate[i] for i in 1:length(observed_sequence))
        mismatches = length(observed_sequence) - matches
        
        ## Calculate likelihood
        likelihood = (1 - avg_error_rate)^matches * avg_error_rate^mismatches
        
        println("• Candidate '$candidate': $matches matches, $mismatches mismatches")
        println("  Likelihood: $(round(likelihood, sigdigits=3))")
    end
    
    println("\nEM Algorithm Interpretation:")
    println("• E-step: Calculate maximum likelihood paths through assembly graph")
    println("• M-step: Update graph structure based on corrected sequences")
    println("• Iterate until convergence or k-mer increment")
    
    return quality_scores
end

probabilistic_demo = demonstrate_probabilistic_assembly_theory()

# ## 8. Genomic Grammar and Zipf's Law
#
# Demonstrates the linguistic properties of genomic sequences and their
# implications for assembly algorithm design.

function demonstrate_zipf_law_genomics()
    println("\nGenomic Grammar and Zipf's Law")
    println("=============================")
    
    ## Simulate k-mer frequency distribution following power law
    ranks = 1:20
    alpha = 1.0  # Zipf exponent (≈1 for natural languages and genomes)
    
    println("Zipf's Law in Genomic Sequences:")
    println("Frequency ∝ Rank^(-α), where α ≈ 1")
    println()
    
    total_frequency = sum(1.0 / rank^alpha for rank in ranks)
    
    for rank in ranks[1:10]  # Show first 10
        frequency = (1.0 / rank^alpha) / total_frequency
        println("• Rank $rank: Relative frequency = $(round(frequency*100, digits=1))%")
    end
    
    println("\nImplications for Assembly:")
    println("• Rare k-mers may represent genuine biological signal")
    println("• Power-law distributions inform graph cleaning strategies")
    println("• Non-coding regions particularly follow linguistic patterns")
    println("• Assembly algorithms should account for natural frequency distributions")
    
    return ranks, alpha
end

zipf_demo = demonstrate_zipf_law_genomics()

# ## 9. Computational Complexity and Optimization
#
# Demonstrates the computational considerations and optimization strategies
# used in Mycelia's assembly algorithms.

function demonstrate_computational_complexity()
    println("\nComputational Complexity Analysis")
    println("================================")
    
    ## K-mer complexity analysis
    k_values = [11, 21, 31, 41, 51]
    
    println("K-mer Space Complexity:")
    for k in k_values
        possible_kmers = 4^k
        if possible_kmers > 1e12
            println("• k=$k: $(round(possible_kmers/1e12, digits=1)) trillion possible k-mers")
        elseif possible_kmers > 1e9
            println("• k=$k: $(round(possible_kmers/1e9, digits=1)) billion possible k-mers")
        elseif possible_kmers > 1e6
            println("• k=$k: $(round(possible_kmers/1e6, digits=1)) million possible k-mers")
        else
            println("• k=$k: $(Int(possible_kmers)) possible k-mers")
        end
    end
    
    println("\nOptimization Strategies:")
    println("• Progressive k-mer approach: O(log k) practical complexity")
    println("• Canonical k-mer storage: ~50% memory reduction")
    println("• Sparse edge representation: O(E) vs O(V²) space")
    println("• Prime k-mer selection: Reduces redundant analysis")
    
    ## Memory estimation for different genome sizes
    genome_sizes = [1e6, 1e7, 1e8, 1e9]  # 1MB to 1GB
    
    println("\nMemory Requirements (estimated):")
    for genome_size in genome_sizes
        ## Rough estimation: ~10 bytes per unique k-mer
        estimated_kmers = genome_size * 0.8  # Assuming 80% unique k-mers
        memory_gb = (estimated_kmers * 10) / 1e9
        
        size_label = if genome_size >= 1e9
            "$(Int(genome_size/1e9))Gb"
        elseif genome_size >= 1e6
            "$(Int(genome_size/1e6))Mb"
        else
            "$(Int(genome_size))bp"
        end
        
        println("• $size_label genome: ~$(round(memory_gb, digits=1)) GB RAM")
    end
    
    return k_values, genome_sizes
end

complexity_analysis = demonstrate_computational_complexity()

# ## Summary and Integration
#
# This tutorial has demonstrated the key theoretical foundations and practical
# implementations that make Mycelia's assembly algorithms both mathematically
# rigorous and biologically relevant.

println("\n" * "="^60)
println("SUMMARY: Theoretical Foundations Integrated into Mycelia")
println("="^60)

println("\n1. Mathematical K-mer Selection:")
println("   • Error rate formula: k ≥ 1/error_rate - 1")
println("   • Log₄ optimization for sequence length")
println("   • Dynamic prime pattern progression")

println("\n2. Graph Theory Applications:")  
println("   • Explicit vs. inferred edge storage")
println("   • Statistical tip clipping (3σ rule)")
println("   • Jaccard similarity quality metrics")

println("\n3. Probabilistic Framework:")
println("   • EM algorithm for assembly optimization")
println("   • Viterbi algorithm for error correction")
println("   • Quality score integration")

println("\n4. Metagenomic Innovations:")
println("   • MAPQ-aware taxonomic assignment")
println("   • Alignment score weighting")
println("   • Strain-level clustering (99.5% ANI)")

println("\n5. Biological Insights:")
println("   • Zipf's law in genomic sequences")
println("   • Prime k-mer advantages")
println("   • Strain-resolved quality metrics")

println("\n6. Computational Optimization:")
println("   • Progressive k-mer complexity reduction")
println("   • Memory-efficient data structures")
println("   • Hardware-optimized algorithms")

println("\nThese theoretical foundations provide Mycelia with:")
println("✓ Mathematically rigorous algorithms")
println("✓ Biologically relevant approaches") 
println("✓ Scalable computational methods")
println("✓ Comprehensive quality assessment")
println("✓ State-of-the-art metagenomic capabilities")

# This tutorial demonstrates how theoretical research has been systematically
# integrated into practical, production-ready algorithms in Mycelia, providing
# users with cutting-edge tools backed by solid mathematical foundations.