"""
Cross-Validation Pipeline for Hybrid Assembly Quality Assessment - Phase 5.1c

This module implements comprehensive cross-validation framework for comparing:
- Intelligent assembly (Phase 5.1a/5.1b) 
- Iterative maximum likelihood assembly (Phase 5.2a)
- Hybrid assembly quality assessment
- K-fold validation through read partitioning
- Holdout validation through read mapping
- Consensus pangenome generation

Part of the Mycelia bioinformatics package's assembly validation framework.
"""

# =============================================================================
# Core Cross-Validation Framework
# =============================================================================

"""
Main cross-validation function for hybrid assembly quality assessment.
Compares intelligent vs iterative assembly approaches across multiple validation folds.
"""
function mycelia_cross_validation(input_fastq::String;
        k_folds::Int = 5,
        max_k::Int = 101,
        memory_limit::Int = 32_000_000_000,
        output_dir::String = "cross_validation",
        validation_split::Float64 = 0.2,
        verbose::Bool = true)
    start_time = time()

    if verbose
        println("Starting Mycelia Cross-Validation Pipeline")
        println("Input FASTQ: $input_fastq")
        println("K-folds: $k_folds")
        println("Validation split: $(validation_split * 100)%")
        println("Output directory: $output_dir")
        println("Max k-mer size: $max_k")
    end

    # Create output directory structure
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    folds_dir = joinpath(output_dir, "folds")
    results_dir = joinpath(output_dir, "results")
    consensus_dir = joinpath(output_dir, "consensus")

    for dir in [folds_dir, results_dir, consensus_dir]
        if !isdir(dir)
            mkpath(dir)
        end
    end

    # Load and partition reads
    if verbose
        println("Loading reads from $input_fastq...")
    end
    all_reads = collect(FASTX.FASTQ.Reader(open(input_fastq)))

    if verbose
        println("Total reads: $(length(all_reads))")
    end

    # Create k-fold partitions
    fold_partitions = create_kfold_partitions(all_reads, k_folds, validation_split)

    # Track results across all folds
    intelligent_results = Dict{Int, Dict}()
    iterative_results = Dict{Int, Dict}()
    validation_metrics = Dict{Int, Dict}()

    # Process each fold
    for fold in 1:k_folds
        if verbose
            println("\n" * "="^60)
            println("PROCESSING FOLD $fold/$k_folds")
            println("="^60)
        end

        train_reads = fold_partitions[fold][:train]
        validation_reads = fold_partitions[fold][:validation]

        if verbose
            println("Training reads: $(length(train_reads))")
            println("Validation reads: $(length(validation_reads))")
        end

        # Create training FASTQ file for this fold
        train_fastq = joinpath(folds_dir, "fold_$(fold)_train.fastq")
        write_fastq(records = train_reads, filename = train_fastq)

        # Create validation FASTQ file for this fold
        validation_fastq = joinpath(folds_dir, "fold_$(fold)_validation.fastq")
        write_fastq(records = validation_reads, filename = validation_fastq)

        # Run intelligent assembly on training data
        if verbose
            println("\n--- Running Intelligent Assembly ---")
        end
        intelligent_output = joinpath(results_dir, "fold_$(fold)_intelligent")
        intelligent_result = mycelia_assemble(train_reads,
            max_k = max_k,
            memory_limit = memory_limit,
            verbose = false)
        intelligent_results[fold] = intelligent_result

        # Run iterative assembly on training data
        if verbose
            println("\n--- Running Iterative Assembly ---")
        end
        iterative_output = joinpath(results_dir, "fold_$(fold)_iterative")
        iterative_result = mycelia_iterative_assemble(train_fastq,
            max_k = max_k,
            memory_limit = memory_limit,
            output_dir = iterative_output,
            verbose = false)
        iterative_results[fold] = iterative_result

        # Validate both assemblies against holdout data
        if verbose
            println("\n--- Validating Assemblies ---")
        end
        validation_metrics[fold] = validate_assemblies_against_holdout(
            intelligent_result, iterative_result, validation_reads, fold, results_dir, verbose = verbose
        )

        if verbose
            println("Fold $fold completed")
            println("Intelligent assembly final k: $(intelligent_result[:metadata][:final_k])")
            println("Iterative assembly final k: $(iterative_result[:metadata][:final_k])")
        end
    end

    # Generate consensus pangenome and final analysis
    if verbose
        println("\n" * "="^60)
        println("GENERATING CONSENSUS ANALYSIS")
        println("="^60)
    end

    consensus_result = generate_consensus_pangenome(
        intelligent_results, iterative_results, validation_metrics,
        consensus_dir, verbose = verbose
    )

    # Calculate total runtime
    total_runtime = time() - start_time

    if verbose
        println("Cross-validation completed in $(round(total_runtime, digits=2)) seconds")
    end

    # Create comprehensive results
    return finalize_cross_validation(
        intelligent_results, iterative_results, validation_metrics,
        consensus_result, total_runtime, output_dir, verbose = verbose
    )
end

# =============================================================================
# K-Fold Partitioning Functions
# =============================================================================

"""
Create k-fold partitions with validation holdout for cross-validation.
"""
function create_kfold_partitions(reads::Vector{<:FASTX.FASTQ.Record}, k_folds::Int, validation_split::Float64)
    n_reads = length(reads)
    shuffled_indices = Random.shuffle(1:n_reads)

    # Calculate fold size
    fold_size = n_reads ÷ k_folds
    validation_size = Int(round(fold_size * validation_split))
    train_size = fold_size - validation_size

    partitions = Dict{Int, Dict{Symbol, Union{Vector{<:FASTX.FASTQ.Record}, Int}}}()

    for fold in 1:k_folds
        start_idx = (fold - 1) * fold_size + 1
        end_idx = fold == k_folds ? n_reads : fold * fold_size

        fold_indices = shuffled_indices[start_idx:end_idx]

        # Split into training and validation
        validation_indices = fold_indices[1:min(validation_size, length(fold_indices))]
        train_indices = fold_indices[(validation_size + 1):end]

        partitions[fold] = Dict{Symbol, Union{Vector{<:FASTX.FASTQ.Record}, Int}}(
            :train => reads[train_indices],
            :validation => reads[validation_indices],
            :fold_size => length(fold_indices),
            :train_size => length(train_indices),
            :validation_size => length(validation_indices)
        )
    end

    return partitions
end

# =============================================================================
# Assembly Validation Functions
# =============================================================================

"""
Validate assemblies against holdout validation data through read mapping.
"""
function validate_assemblies_against_holdout(
        intelligent_result::Dict, iterative_result::Dict,
        validation_reads::Vector{<:FASTX.FASTQ.Record},
        fold::Int, results_dir::String; verbose::Bool = false)

    # Extract final assemblies
    intelligent_assembly = intelligent_result[:final_assembly]
    iterative_assembly = iterative_result[:final_assembly]

    if verbose
        println("Validating against $(length(validation_reads)) holdout reads")
    end

    # Calculate mapping metrics for intelligent assembly
    intelligent_metrics = calculate_assembly_mapping_metrics(
        intelligent_assembly, validation_reads, "intelligent", verbose = verbose
    )

    # Calculate mapping metrics for iterative assembly
    iterative_metrics = calculate_assembly_mapping_metrics(
        iterative_assembly, validation_reads, "iterative", verbose = verbose
    )

    # Compare assembly statistics
    assembly_comparison = compare_assembly_statistics(
        intelligent_result, iterative_result, verbose = verbose
    )

    # Create comprehensive validation report
    validation_report = Dict(
        :fold => fold,
        :validation_reads_count => length(validation_reads),
        :intelligent_metrics => intelligent_metrics,
        :iterative_metrics => iterative_metrics,
        :assembly_comparison => assembly_comparison,
        :timestamp => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    )

    # Save validation report
    report_file = joinpath(results_dir, "fold_$(fold)_validation_report.json")
    save_validation_report(validation_report, report_file)

    if verbose
        println("Validation metrics saved to $report_file")
        println("Intelligent mapping rate: $(round(intelligent_metrics[:mapping_rate] * 100, digits=2))%")
        println("Iterative mapping rate: $(round(iterative_metrics[:mapping_rate] * 100, digits=2))%")
    end

    return validation_report
end

"""
Calculate assembly mapping metrics for validation reads.
"""
function calculate_assembly_mapping_metrics(
        assembly, validation_reads::Vector{<:FASTX.FASTQ.Record},
        assembly_type::String; verbose::Bool = false)

    # Handle different assembly formats
    assembly_sequences = if applicable(MetaGraphsNext.labels, assembly)
        # Graph-based assembly - extract k-mers
        [string(kmer) for kmer in MetaGraphsNext.labels(assembly)]
    elseif isa(assembly, Vector{String})
        # Sequence-based assembly
        assembly
    else
        # Try to extract sequences from various formats
        String[]
    end

    if isempty(assembly_sequences)
        if verbose
            println("Warning: No sequences found in $assembly_type assembly")
        end
        return Dict(
            :assembly_type => assembly_type,
            :total_reads => length(validation_reads),
            :mapped_reads => 0,
            :mapping_rate => 0.0,
            :average_coverage => 0.0,
            :assembly_sequences => 0
        )
    end

    # Create assembly k-mer index for fast lookup
    assembly_kmers = Set{String}()
    for seq in assembly_sequences
        if length(seq) >= 3  # Use minimum k-mer size of 3
            for i in 1:(length(seq) - 3 + 1)
                push!(assembly_kmers, seq[i:(i + 2)])
            end
        end
    end

    # Map validation reads against assembly
    mapped_reads = 0
    total_coverage = 0.0

    for read in validation_reads
        read_seq = FASTX.sequence(String, read)
        if length(read_seq) >= 3
            read_kmers = Set{String}()
            for i in 1:(length(read_seq) - 3 + 1)
                push!(read_kmers, read_seq[i:(i + 2)])
            end

            # Calculate coverage as fraction of read k-mers found in assembly
            overlap = length(intersect(read_kmers, assembly_kmers))
            coverage = overlap / length(read_kmers)
            total_coverage += coverage

            # Consider mapped if >50% of k-mers are found
            if coverage > 0.5
                mapped_reads += 1
            end
        end
    end

    mapping_rate = mapped_reads / length(validation_reads)
    average_coverage = total_coverage / length(validation_reads)

    return Dict(
        :assembly_type => assembly_type,
        :total_reads => length(validation_reads),
        :mapped_reads => mapped_reads,
        :mapping_rate => mapping_rate,
        :average_coverage => average_coverage,
        :assembly_sequences => length(assembly_sequences),
        :assembly_kmers => length(assembly_kmers)
    )
end

"""
Compare statistical properties of intelligent vs iterative assemblies.
"""
function compare_assembly_statistics(intelligent_result::Dict, iterative_result::Dict; verbose::Bool = false)
    intel_meta = get(intelligent_result, :metadata, Dict())
    iter_meta = get(iterative_result, :metadata, Dict())

    comparison = Dict(
        :runtime_comparison => Dict(
            :intelligent_runtime => get(intel_meta, :total_runtime, 0.0),
            :iterative_runtime => get(iter_meta, :total_runtime, 0.0),
            :runtime_ratio => get(iter_meta, :total_runtime, 1.0) /
                              max(get(intel_meta, :total_runtime, 1.0), 1e-10)
        ),
        :k_progression_comparison => Dict(
            :intelligent_k_progression => get(intel_meta, :k_progression, Int[]),
            :iterative_k_progression => get(iter_meta, :k_progression, Int[]),
            :intelligent_final_k => get(intel_meta, :final_k, 0),
            :iterative_final_k => get(iter_meta, :final_k, 0)
        ),
        :assembly_type_comparison => Dict(
            :intelligent_type => get(intel_meta, :assembly_type, "unknown"),
            :iterative_type => get(iter_meta, :assembly_type, "unknown"),
            :intelligent_version => get(intel_meta, :version, "unknown"),
            :iterative_version => get(iter_meta, :version, "unknown")
        )
    )

    # Add improvement metrics if available
    if haskey(iter_meta, :total_improvements)
        comparison[:improvement_metrics] = Dict(
            :iterative_total_improvements => iter_meta[:total_improvements],
            :iterative_total_iterations => iter_meta[:total_iterations]
        )
    end

    if verbose
        println("Runtime comparison:")
        println("  Intelligent: $(round(intel_meta[:total_runtime], digits=2))s")
        println("  Iterative: $(round(iter_meta[:total_runtime], digits=2))s")
        println("  Ratio: $(round(comparison[:runtime_comparison][:runtime_ratio], digits=2))")
    end

    return comparison
end

# =============================================================================
# Consensus Pangenome Generation
# =============================================================================

"""
Generate consensus pangenome from cross-validation results.
"""
function generate_consensus_pangenome(intelligent_results::Dict, iterative_results::Dict,
        validation_metrics::Dict, consensus_dir::String; verbose::Bool = false)
    if verbose
        println("Generating consensus pangenome analysis...")
    end

    n_folds = length(intelligent_results)

    # Aggregate performance metrics across folds
    intel_mapping_rates = Float64[]
    iter_mapping_rates = Float64[]
    intel_runtimes = Float64[]
    iter_runtimes = Float64[]

    for fold in 1:n_folds
        push!(intel_mapping_rates, validation_metrics[fold][:intelligent_metrics][:mapping_rate])
        push!(iter_mapping_rates, validation_metrics[fold][:iterative_metrics][:mapping_rate])
        push!(intel_runtimes, intelligent_results[fold][:metadata][:total_runtime])
        push!(iter_runtimes, iterative_results[fold][:metadata][:total_runtime])
    end

    # Calculate summary statistics
    consensus_stats = Dict(
        :intelligent_performance => Dict(
            :mean_mapping_rate => Statistics.mean(intel_mapping_rates),
            :std_mapping_rate => Statistics.std(intel_mapping_rates),
            :mean_runtime => Statistics.mean(intel_runtimes),
            :std_runtime => Statistics.std(intel_runtimes)
        ),
        :iterative_performance => Dict(
            :mean_mapping_rate => Statistics.mean(iter_mapping_rates),
            :std_mapping_rate => Statistics.std(iter_mapping_rates),
            :mean_runtime => Statistics.mean(iter_runtimes),
            :std_runtime => Statistics.std(iter_runtimes)
        ),
        :statistical_comparison => Dict(
            :mapping_rate_improvement => Statistics.mean(iter_mapping_rates) -
                                         Statistics.mean(intel_mapping_rates),
            :runtime_overhead => Statistics.mean(iter_runtimes) /
                                 Statistics.mean(intel_runtimes),
            :n_folds => n_folds
        )
    )

    # Determine recommended approach
    mapping_improvement = consensus_stats[:statistical_comparison][:mapping_rate_improvement]
    runtime_overhead = consensus_stats[:statistical_comparison][:runtime_overhead]

    recommendation = if mapping_improvement > 0.05 && runtime_overhead < 3.0
        "iterative"  # Significantly better mapping with reasonable runtime
    elseif mapping_improvement < -0.05
        "intelligent"  # Significantly better mapping
    elseif runtime_overhead > 5.0
        "intelligent"  # Much faster
    else
        "hybrid"  # Close performance, recommend hybrid approach
    end

    consensus_stats[:recommendation] = Dict(
        :approach => recommendation,
        :reasoning => generate_recommendation_reasoning(mapping_improvement, runtime_overhead),
        :confidence => calculate_recommendation_confidence(intel_mapping_rates, iter_mapping_rates)
    )

    # Save consensus analysis
    consensus_file = joinpath(consensus_dir, "consensus_analysis.json")
    save_consensus_analysis(consensus_stats, consensus_file)

    if verbose
        println("Consensus analysis completed:")
        println("  Intelligent mean mapping: $(round(consensus_stats[:intelligent_performance][:mean_mapping_rate] * 100, digits=2))%")
        println("  Iterative mean mapping: $(round(consensus_stats[:iterative_performance][:mean_mapping_rate] * 100, digits=2))%")
        println("  Recommended approach: $(recommendation)")
        println("  Consensus saved to: $consensus_file")
    end

    return consensus_stats
end

"""
Generate reasoning for assembly approach recommendation.
"""
function generate_recommendation_reasoning(mapping_improvement::Float64, runtime_overhead::Float64)::String
    if mapping_improvement > 0.05 && runtime_overhead < 3.0
        return "Iterative assembly shows significant mapping improvement (>5%) with acceptable runtime overhead (<3x)"
    elseif mapping_improvement < -0.05
        return "Intelligent assembly shows significantly better mapping performance"
    elseif runtime_overhead > 5.0
        return "Intelligent assembly is much faster (>5x) with comparable mapping performance"
    else
        return "Performance is comparable between approaches, hybrid strategy recommended"
    end
end

"""
Calculate confidence in recommendation based on variance across folds.
"""
function calculate_recommendation_confidence(intel_rates::Vector{Float64}, iter_rates::Vector{Float64})::Float64
    intel_cv = Statistics.std(intel_rates) / Statistics.mean(intel_rates)
    iter_cv = Statistics.std(iter_rates) / Statistics.mean(iter_rates)
    avg_cv = (intel_cv + iter_cv) / 2
    return max(0.0, 1.0 - avg_cv)  # Higher confidence with lower coefficient of variation
end

# =============================================================================
# Utility and Reporting Functions
# =============================================================================

"""
Save validation report to JSON file.
"""
function save_validation_report(report::Dict, filename::String)
    open(filename, "w") do f
        JSON.print(f, report)
    end
end

"""
Save consensus analysis to JSON file.
"""
function save_consensus_analysis(analysis::Dict, filename::String)
    open(filename, "w") do f
        JSON.print(f, analysis)
    end
end

"""
Finalize cross-validation results with comprehensive summary.
"""
function finalize_cross_validation(intelligent_results::Dict, iterative_results::Dict,
        validation_metrics::Dict, consensus_result::Dict,
        total_runtime::Float64, output_dir::String; verbose::Bool = false)
    final_result = Dict(
        :cross_validation_type => "hybrid_assembly_quality_assessment",
        :version => "Phase_5.1c",
        :total_runtime => total_runtime,
        :n_folds => length(intelligent_results),
        :output_directory => output_dir,
        :intelligent_results => intelligent_results,
        :iterative_results => iterative_results,
        :validation_metrics => validation_metrics,
        :consensus_analysis => consensus_result,
        :timestamp => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    )

    # Save final results
    results_file = joinpath(output_dir, "cross_validation_results.json")
    save_consensus_analysis(final_result, results_file)

    if verbose
        println("Cross-validation results saved to: $results_file")
        println("Recommended approach: $(consensus_result[:recommendation][:approach])")
    end

    return final_result
end

"""
Generate cross-validation summary report.
"""
function cross_validation_summary(result::Dict)::String
    consensus = result[:consensus_analysis]

    report = "Mycelia Cross-Validation Summary\n"
    report *= "=" * "="^50 * "\n\n"

    report *= "Cross-validation Type: $(result[:cross_validation_type])\n"
    report *= "Version: $(result[:version])\n"
    report *= "Total Runtime: $(round(result[:total_runtime], digits=2)) seconds\n"
    report *= "Number of Folds: $(result[:n_folds])\n\n"

    report *= "Performance Comparison:\n"
    intel_perf = consensus[:intelligent_performance]
    iter_perf = consensus[:iterative_performance]

    report *= "  Intelligent Assembly:\n"
    report *= "    Mean Mapping Rate: $(round(intel_perf[:mean_mapping_rate] * 100, digits=2))% ± $(round(intel_perf[:std_mapping_rate] * 100, digits=2))%\n"
    report *= "    Mean Runtime: $(round(intel_perf[:mean_runtime], digits=2))s ± $(round(intel_perf[:std_runtime], digits=2))s\n\n"

    report *= "  Iterative Assembly:\n"
    report *= "    Mean Mapping Rate: $(round(iter_perf[:mean_mapping_rate] * 100, digits=2))% ± $(round(iter_perf[:std_mapping_rate] * 100, digits=2))%\n"
    report *= "    Mean Runtime: $(round(iter_perf[:mean_runtime], digits=2))s ± $(round(iter_perf[:std_runtime], digits=2))s\n\n"

    report *= "Statistical Comparison:\n"
    stat_comp = consensus[:statistical_comparison]
    report *= "  Mapping Rate Improvement: $(round(stat_comp[:mapping_rate_improvement] * 100, digits=2))%\n"
    report *= "  Runtime Overhead: $(round(stat_comp[:runtime_overhead], digits=2))x\n\n"

    report *= "Recommendation:\n"
    rec = consensus[:recommendation]
    report *= "  Approach: $(rec[:approach])\n"
    report *= "  Confidence: $(round(rec[:confidence] * 100, digits=1))%\n"
    report *= "  Reasoning: $(rec[:reasoning])\n"

    return report
end

"""
Test function for cross-validation with sample data.
"""
function test_cross_validation()
    try
        println("Testing Mycelia Cross-Validation Pipeline")

        # Create temporary test FASTQ file
        temp_dir = mktempdir()
        test_fastq = joinpath(temp_dir, "test_reads.fastq")

        # Generate test sequences
        sequences = [
            "ATCGATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC",
            "ATCGATCGATCGATCGTTCGATCG",
            "TCGATCGATCGATCGATCGATCGA",
            "GCGATCGATCGATCGATCGATCGT",
            "ACGATCGATCGATCGATCGATCGT",
            "CCGATCGATCGATCGATCGATCGA"
        ]

        test_records = [FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
                        for (i, seq) in enumerate(sequences)]

        # Write test FASTQ
        write_fastq(records = test_records, filename = test_fastq)

        # Test cross-validation with minimal settings
        output_dir = joinpath(temp_dir, "test_cv")
        result = mycelia_cross_validation(test_fastq,
            k_folds = 3,
            max_k = 13,
            output_dir = output_dir,
            verbose = false)

        # Cleanup
        rm(temp_dir, recursive = true)

        println("Cross-validation test completed successfully!")
        return merge(result, Dict(:status => :success))

    catch e
        println("Cross-validation test failed: $e")
        return Dict(:status => :error, :error => string(e))
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform bootstrap validation for assembly parameter optimization.

# Arguments
- `fastq_file::String`: Path to input FASTQ file
- `n_bootstrap::Int`: Number of bootstrap samples (default: 100)
- `sample_ratio::Float64`: Fraction of reads to sample in each bootstrap (default: 0.8)
- `k_values::Vector{Int}`: K-mer sizes to test (default: [21, 31, 41, 51])
- `output_dir::String`: Output directory for results

# Returns
Named tuple containing:
- `optimal_k::Int`: Optimal k-mer size based on bootstrap validation
- `confidence_intervals::Dict`: Confidence intervals for assembly metrics
- `bootstrap_results::Vector{Dict}`: Results from each bootstrap sample

# Details
- Implements bootstrap resampling for robust parameter validation
- Evaluates assembly quality across multiple parameter combinations
- Calculates confidence intervals for assembly statistics
- Supports parameter optimization with uncertainty quantification
"""
function bootstrap_assembly_validation(fastq_file::String;
        n_bootstrap::Int = 100,
        sample_ratio::Float64 = 0.8,
        k_values::Vector{Int} = [21, 31, 41, 51],
        output_dir::String = "bootstrap_validation")
    mkpath(output_dir)
    @info "Starting bootstrap validation with $n_bootstrap samples"

    # Read all records from FASTQ file
    all_records = Vector{FASTX.FASTQ.Record}()
    FASTX.FASTQ.Reader(open(fastq_file, "r")) do reader
        for record in reader
            push!(all_records, record)
        end
    end

    total_reads = length(all_records)
    sample_size = Int(round(total_reads * sample_ratio))

    bootstrap_results = Vector{Dict}()

    for bootstrap_idx in 1:n_bootstrap
        @info "Processing bootstrap sample $bootstrap_idx/$n_bootstrap"

        # Bootstrap sampling with replacement
        Random.seed!(bootstrap_idx)  # For reproducibility
        sampled_indices = rand(1:total_reads, sample_size)
        sampled_records = all_records[sampled_indices]

        # Write bootstrap sample to temporary file
        bootstrap_fastq = joinpath(output_dir, "bootstrap_$(bootstrap_idx).fastq")
        write_fastq(records = sampled_records, filename = bootstrap_fastq)

        bootstrap_result = Dict{String, Any}()
        bootstrap_result["sample_id"] = bootstrap_idx
        bootstrap_result["k_results"] = Dict{Int, Dict}()

        # Test each k-mer size
        for k in k_values
            @info "Testing k=$k for bootstrap $bootstrap_idx"

            try
                # Perform assembly with this k value
                assembly_result = test_assembly_with_k(bootstrap_fastq, k, output_dir, bootstrap_idx)
                bootstrap_result["k_results"][k] = assembly_result

            catch e
                @warn "Failed assembly for k=$k in bootstrap $bootstrap_idx: $e"
                bootstrap_result["k_results"][k] = Dict("error" => string(e))
            end
        end

        push!(bootstrap_results, bootstrap_result)

        # Cleanup temporary FASTQ
        rm(bootstrap_fastq, force = true)
    end

    # Analyze bootstrap results
    confidence_intervals = calculate_bootstrap_confidence_intervals(bootstrap_results, k_values)
    optimal_k = select_optimal_k_bootstrap(bootstrap_results, k_values)

    # Save results
    results_file = joinpath(output_dir, "bootstrap_results.json")
    open(results_file, "w") do io
        JSON3.write(io,
            Dict(
                "optimal_k" => optimal_k,
                "confidence_intervals" => confidence_intervals,
                "bootstrap_results" => bootstrap_results,
                "parameters" => Dict(
                    "n_bootstrap" => n_bootstrap,
                    "sample_ratio" => sample_ratio,
                    "k_values" => k_values
                )
            ))
    end

    return (; optimal_k, confidence_intervals, bootstrap_results)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Test assembly with a specific k-mer size for bootstrap validation.
"""
function test_assembly_with_k(fastq_file::String, k::Int, output_dir::String, bootstrap_idx::Int)
    assembly_dir = joinpath(output_dir, "assembly_k$(k)_b$(bootstrap_idx)")
    mkpath(assembly_dir)

    try
        # Use megahit for consistent assembly testing
        assembly_result = run_megahit(
            fastq1 = fastq_file,
            outdir = assembly_dir,
            k_list = string(k),
            min_contig_len = 200
        )

        # Calculate assembly statistics
        contigs_file = assembly_result.contigs
        if isfile(contigs_file)
            stats = calculate_basic_assembly_stats(contigs_file)
            stats["assembly_success"] = true
            return stats
        else
            return Dict("assembly_success" => false, "error" => "No contigs file produced")
        end

    catch e
        return Dict("assembly_success" => false, "error" => string(e))
    finally
        # Cleanup assembly directory to save space
        rm(assembly_dir, recursive = true, force = true)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate basic assembly statistics for bootstrap validation.
"""
function calculate_basic_assembly_stats(contigs_file::String)
    stats = Dict{String, Any}()

    sequences = Vector{BioSequences.LongDNA{4}}()
    FASTX.FASTA.Reader(open(contigs_file, "r")) do reader
        for record in reader
            seq = FASTX.sequence(BioSequences.LongDNA{4}, record)
            push!(sequences, seq)
        end
    end

    if isempty(sequences)
        return Dict("num_contigs" => 0, "total_length" => 0, "n50" => 0)
    end

    lengths = [length(seq) for seq in sequences]
    total_length = sum(lengths)

    # Calculate N50
    sorted_lengths = sort(lengths, rev = true)
    cumulative_length = 0
    n50 = 0

    for length in sorted_lengths
        cumulative_length += length
        if cumulative_length >= total_length / 2
            n50 = length
            break
        end
    end

    stats["num_contigs"] = length(sequences)
    stats["total_length"] = total_length
    stats["n50"] = n50
    stats["min_length"] = minimum(lengths)
    stats["max_length"] = maximum(lengths)
    stats["mean_length"] = Statistics.mean(lengths)

    return stats
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate confidence intervals from bootstrap results.
"""
function calculate_bootstrap_confidence_intervals(bootstrap_results::Vector{Dict}, k_values::Vector{Int})
    confidence_intervals = Dict{Int, Dict}()

    for k in k_values
        # Collect all successful results for this k
        k_stats = []
        for result in bootstrap_results
            k_result = get(result["k_results"], k, nothing)
            if k_result !== nothing && get(k_result, "assembly_success", false)
                push!(k_stats, k_result)
            end
        end

        if isempty(k_stats)
            confidence_intervals[k] = Dict("error" => "No successful assemblies")
            continue
        end

        # Calculate confidence intervals for key metrics
        ci_dict = Dict{String, Dict}()

        for metric in ["n50", "total_length", "num_contigs", "mean_length"]
            values = [get(stat, metric, 0) for stat in k_stats if haskey(stat, metric)]

            if !isempty(values)
                sorted_values = sort(values)
                n = length(sorted_values)

                # 95% confidence interval
                lower_idx = max(1, Int(ceil(0.025 * n)))
                upper_idx = min(n, Int(floor(0.975 * n)))

                ci_dict[metric] = Dict(
                    "mean" => Statistics.mean(values),
                    "lower_95" => sorted_values[lower_idx],
                    "upper_95" => sorted_values[upper_idx],
                    "median" => Statistics.median(values)
                )
            end
        end

        confidence_intervals[k] = ci_dict
    end

    return confidence_intervals
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Select optimal k-mer size based on bootstrap validation results.
"""
function select_optimal_k_bootstrap(bootstrap_results::Vector{Dict}, k_values::Vector{Int})
    k_scores = Dict{Int, Float64}()

    for k in k_values
        scores = Float64[]

        for result in bootstrap_results
            k_result = get(result["k_results"], k, nothing)
            if k_result !== nothing && get(k_result, "assembly_success", false)
                # Simple scoring based on N50 and total length
                n50 = get(k_result, "n50", 0)
                total_length = get(k_result, "total_length", 0)
                num_contigs = get(k_result, "num_contigs", 1)

                # Score favors higher N50, reasonable total length, fewer contigs
                score = n50 * log(total_length + 1) / sqrt(num_contigs + 1)
                push!(scores, score)
            end
        end

        k_scores[k] = isempty(scores) ? 0.0 : Statistics.mean(scores)
    end

    # Return k with highest average score
    return argmax(k_scores)
end
