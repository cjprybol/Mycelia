"""
Iterative Maximum Likelihood Assembly - Phase 5.2a Implementation

This module implements the iterative maximum likelihood assembly system with:
- Complete FASTQ I/O processing per iteration
- Statistical path resampling with likelihood calculations
- Viterbi algorithm integration for optimal path finding
- Timestamped output files for tracking read evolution
- Memory-efficient read set processing

Part of the Mycelia bioinformatics package's iterative assembler framework.
"""

# =============================================================================
# Core Iterative Assembly Framework
# =============================================================================

"""
Main iterative maximum likelihood assembly function.
Processes entire read sets per iteration with complete FASTQ I/O tracking.
Enhanced with performance optimizations, caching, and progress tracking.
"""
function mycelia_iterative_assemble(input_fastq::String; 
                                   max_k::Int = 101,
                                   memory_limit::Int = 32_000_000_000,
                                   output_dir::String = "iterative_assembly",
                                   max_iterations_per_k::Int = 10,
                                   improvement_threshold::Float64 = 0.05,
                                   graph_mode::Symbol = :canonical,
                                   verbose::Bool = true,
                                   enable_parallel::Bool = false,
                                   batch_size::Int = 10000,
                                   enable_checkpointing::Bool = true,
                                   checkpoint_interval::Int = 5)
    
    start_time = time()
    
    if verbose
        println("Starting Mycelia Iterative Maximum Likelihood Assembly")
        println("Input FASTQ: $input_fastq")
        println("Output directory: $output_dir") 
        println("Memory limit: $(memory_limit ÷ 1_000_000_000) GB")
        println("Max k-mer size: $max_k")
        println("Graph mode: $graph_mode")
        println("Parallel processing: $(enable_parallel ? "enabled ($(Threads.nthreads()) threads)" : "disabled")")
        println("Batch size: $batch_size")
        println("Checkpointing: $(enable_checkpointing ? "enabled (every $checkpoint_interval iterations)" : "disabled")")
    end
    
    # Create output directory structure
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # Create subdirectories for organization
    checkpoints_dir = joinpath(output_dir, "checkpoints")
    graphs_dir = joinpath(output_dir, "graphs")
    progress_dir = joinpath(output_dir, "progress")
    
    if enable_checkpointing
        mkpath(checkpoints_dir)
        mkpath(graphs_dir)
        mkpath(progress_dir)
    end
    
    # Check for existing checkpoint to resume from
    checkpoint_file = joinpath(checkpoints_dir, "latest_checkpoint.json")
    resume_data = nothing
    
    if enable_checkpointing && isfile(checkpoint_file)
        try
            resume_data = JSON.parsefile(checkpoint_file)
            if verbose
                println("Found existing checkpoint. Resume from k=$(resume_data["current_k"]), iteration=$(resume_data["current_iteration"])")
            end
        catch e
            if verbose
                println("Warning: Could not load checkpoint file: $e")
            end
        end
    end
    
    # Initialize or resume from checkpoint
    if resume_data !== nothing
        # Resume from checkpoint
        k = resume_data["current_k"]
        k_progression = resume_data["k_progression"]
        iteration_history = Dict{Int, Vector{Dict{Symbol, Any}}}()
        for (k_str, hist) in resume_data["iteration_history"]
            iteration_history[parse(Int, k_str)] = hist
        end
        total_improvements = resume_data["total_improvements"]
        current_fastq_file = resume_data["current_fastq_file"]
        
        if verbose
            println("Resumed from checkpoint: k=$k, total_improvements=$total_improvements")
        end
    else
        # Initialize fresh run
        if verbose
            println("Reading initial FASTQ file...")
        end
        initial_reads = collect(FASTX.FASTQ.Reader(open(input_fastq)))
        k = find_initial_k(initial_reads)  # Reuse from intelligent-assembly.jl
        
        if verbose
            println("Initial k-mer size: $k (prime: $(Primes.isprime(k)))")
            println("Total reads: $(length(initial_reads))")
        end
        
        # Track progress across all iterations
        k_progression = Int[]
        iteration_history = Dict{Int, Vector{Dict{Symbol, Any}}}()
        total_improvements = 0
        current_fastq_file = input_fastq
    end
    
    # Main k-mer progression loop
    while k <= max_k
        if verbose
            println("\n" * "="^60)
            println("PROCESSING K-MER SIZE: $k (prime: $(Primes.isprime(k)))")
            println("="^60)
        end
        
        push!(k_progression, k)
        iteration_history[k] = Dict{Symbol, Any}[]
        
        iteration = 1
        improvements_this_k = 0
        
        # Iterative improvement loop for current k
        while iteration <= max_iterations_per_k
            if verbose
                println("\n--- Iteration $iteration for k=$k ---")
            end
            
            iteration_start = time()
            timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
            
            # Read in entire FASTQ set for this iteration
            if verbose
                println("Reading FASTQ file: $current_fastq_file")
            end
            current_reads = collect(FASTX.FASTQ.Reader(open(current_fastq_file)))
            
            # Build qualmer graph from current read set
            if verbose
                println("Building qualmer graph with k=$k...")
            end
            graph = Mycelia.Rhizomorph.build_qualmer_graph(current_reads, k; mode=graph_mode)
            num_kmers = length(MetaGraphsNext.labels(graph))
            
            if verbose
                println("Graph built: $num_kmers unique k-mers")
            end
            
            # Check memory usage
            if !check_memory_limits(graph, memory_limit)
                if verbose
                    println("Memory limit reached. Stopping at k=$k")
                end
                break
            end
            
            # Process each read for likelihood improvement with performance optimizations
            if verbose
                println("Processing reads for likelihood improvements...")
            end
            updated_reads, improvements_made = improve_read_set_likelihood(
                current_reads, graph, k, 
                verbose=verbose, 
                batch_size=batch_size,
                enable_parallel=enable_parallel,
                graph_mode=graph_mode
            )
            
            # Calculate iteration metrics
            iteration_time = time() - iteration_start
            improvement_rate = improvements_made / length(current_reads)
            
            iteration_stats = Dict(
                :iteration => iteration,
                :k => k,
                :timestamp => timestamp,
                :improvements_made => improvements_made,
                :total_reads => length(current_reads),
                :improvement_rate => improvement_rate,
                :runtime_seconds => iteration_time,
                :memory_kmers => num_kmers
            )
            
            push!(iteration_history[k], iteration_stats)
            improvements_this_k += improvements_made
            total_improvements += improvements_made
            
            if verbose
                println("Improvements made: $improvements_made ($(round(improvement_rate * 100, digits=2))%)")
                println("Iteration runtime: $(round(iteration_time, digits=2)) seconds")
            end
            
            # Write out complete updated read set to timestamped FASTQ
            output_file = joinpath(output_dir, "reads_k$(k)_iter$(iteration)_$(timestamp).fastq")
            write_fastq(records=updated_reads, filename=output_file)  # Use existing function
            
            if verbose
                println("Wrote $(length(updated_reads)) reads to $output_file")
            end
            
            # Create checkpoint if enabled and at checkpoint interval
            if enable_checkpointing && iteration % checkpoint_interval == 0
                checkpoint_data = Dict(
                    "current_k" => k,
                    "current_iteration" => iteration,
                    "k_progression" => k_progression,
                    "iteration_history" => iteration_history,
                    "total_improvements" => total_improvements,
                    "current_fastq_file" => output_file,
                    "timestamp" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
                    "runtime_so_far" => time() - start_time
                )
                
                try
                    open(checkpoint_file, "w") do f
                        JSON.print(f, checkpoint_data, 2)
                    end
                    
                    # Also save progress summary
                    progress_file = joinpath(progress_dir, "progress_k$(k)_iter$(iteration).json")
                    open(progress_file, "w") do f
                        JSON.print(f, Dict(
                            "k" => k,
                            "iteration" => iteration,
                            "improvements_this_iteration" => improvements_made,
                            "total_improvements" => total_improvements,
                            "improvement_rate" => improvement_rate,
                            "runtime_seconds" => iteration_time,
                            "timestamp" => timestamp
                        ), 2)
                    end
                    
                    if verbose
                        println("Checkpoint saved: k=$k, iteration=$iteration")
                    end
                catch e
                    if verbose
                        println("Warning: Could not save checkpoint: $e")
                    end
                end
            end
            
            # Check if we should continue with this k or move to next
            if sufficient_improvements(improvements_made, length(current_reads), improvement_threshold, 
                                     iteration_history=iteration_history[k])
                if verbose
                    println("Sufficient improvements detected. Continuing with k=$k")
                end
                iteration += 1
                current_fastq_file = output_file  # Use updated reads for next iteration
            else
                if verbose
                    println("Insufficient improvements or convergence detected. Moving to next k-mer size")
                end
                break
            end
        end
        
        if verbose
            println("\nCompleted k=$k processing:")
            println("  Total iterations: $(length(iteration_history[k]))")
            println("  Total improvements: $improvements_this_k")
            println("  Final improvement rate: $(round(improvements_this_k / length(current_reads), digits=4))")
        end
        
        # Move to next prime k-mer size
        next_k = next_prime_k(k, max_k=max_k)
        if next_k == k
            if verbose
                println("No larger prime k-mer size available. Stopping at k=$k")
            end
            break
        end
        k = next_k
    end
    
    # Calculate total runtime
    total_runtime = time() - start_time
    
    if verbose
        println("\n" * "="^60)
        println("ITERATIVE ASSEMBLY COMPLETE")
        println("="^60)
        println("Total runtime: $(round(total_runtime, digits=2)) seconds")
        println("K-mer sizes processed: $(sort(k_progression))")
        println("Total improvements: $total_improvements")
    end
    
    # Finalize iterative assembly
    return finalize_iterative_assembly(output_dir, k_progression, iteration_history, total_runtime, verbose=verbose)
end

# =============================================================================
# Read Likelihood Improvement Functions
# =============================================================================

"""
Improve likelihood of entire read set using current graph and k-mer size.
Returns updated reads and count of improvements made.
Uses memory-efficient batch processing for large datasets.
"""
function improve_read_set_likelihood(reads::Vector{<:FASTX.FASTQ.Record}, graph, k::Int; 
                                   verbose::Bool = false, 
                                   batch_size::Int = 10000,
                                   enable_parallel::Bool = false,
                                   graph_mode::Symbol = :canonical)::Tuple{Vector{FASTX.FASTQ.Record}, Int}
    
    total_reads = length(reads)
    updated_reads = Vector{FASTX.FASTQ.Record}(undef, total_reads)
    improvements_made = 0
    
    if verbose
        println("  Processing $total_reads reads in batches of $batch_size")
    end
    
    # Process in batches for memory efficiency
    for batch_start in 1:batch_size:total_reads
        batch_end = min(batch_start + batch_size - 1, total_reads)
        batch_reads = reads[batch_start:batch_end]
        
        batch_improvements = 0
        
        if enable_parallel && Threads.nthreads() > 1
            # Parallel processing for large batches
            batch_results = Vector{Tuple{FASTX.FASTQ.Record, Bool}}(undef, length(batch_reads))
            
            Threads.@threads for i in eachindex(batch_reads)
                improved_read, was_improved = improve_read_likelihood(batch_reads[i], graph, k; graph_mode=graph_mode)
                batch_results[i] = (improved_read, was_improved)
            end
            
            # Collect results
            for (i, (improved_read, was_improved)) in enumerate(batch_results)
                updated_reads[batch_start + i - 1] = improved_read
                if was_improved
                    batch_improvements += 1
                end
            end
        else
            # Sequential processing
            for (i, read) in enumerate(batch_reads)
                improved_read, was_improved = improve_read_likelihood(read, graph, k; graph_mode=graph_mode)
                updated_reads[batch_start + i - 1] = improved_read
                
                if was_improved
                    batch_improvements += 1
                end
            end
        end
        
        improvements_made += batch_improvements
        
        # Progress reporting and garbage collection
        if verbose
            batch_num = div(batch_start - 1, batch_size) + 1
            total_batches = div(total_reads - 1, batch_size) + 1
            improvement_rate = batch_improvements / length(batch_reads) * 100
            println("    Batch $batch_num/$total_batches: $(batch_improvements)/$(length(batch_reads)) improvements ($(round(improvement_rate, digits=1))%)")
        end
        
        # Force garbage collection between batches for memory efficiency
        if batch_end < total_reads
            GC.gc()
        end
    end
    
    if verbose
        total_improvement_rate = improvements_made / total_reads * 100
        println("  Total improvements: $improvements_made/$total_reads ($(round(total_improvement_rate, digits=1))%)")
    end
    
    return updated_reads, improvements_made
end

"""
Improve likelihood of a single read using maximum likelihood path finding.
Returns improved read and boolean indicating if improvement was made.
"""
function improve_read_likelihood(read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol = :canonical)::Tuple{FASTX.FASTQ.Record, Bool}
    # Extract sequence and quality
    original_seq = FASTX.sequence(String, read)
    original_qual = FASTX.quality(read)
    read_id = FASTX.identifier(read)
    
    # Skip reads shorter than k
    if length(original_seq) < k
        return read, false
    end
    
    # Find optimal path through graph using enhanced statistical path improvement
    improved_read, likelihood_improvement = find_optimal_sequence_path(read, graph, k; graph_mode=graph_mode)
    
    # Only update if significant improvement
    if likelihood_improvement > 0.01  # Threshold for meaningful improvement
        return improved_read, true
    else
        return read, false
    end
end

"""
Find optimal sequence path through graph using maximum likelihood principles.
Returns improved sequence and likelihood improvement score.
"""
function find_optimal_sequence_path(sequence::AbstractString, quality, graph, k::Int; graph_mode::Symbol = :canonical)::Tuple{String, Float64}
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for find_optimal_sequence_path"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    improved_record, improvement = find_optimal_sequence_path(record, graph, k; graph_mode=graph_mode)
    return FASTX.sequence(String, improved_record), improvement
end

function find_optimal_sequence_path(read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol = :canonical)::Tuple{FASTX.FASTQ.Record, Float64}
    # Enhanced Statistical Path Improvement (Phase 5.2b)
    # Integrates iterative assembly with existing Viterbi algorithms from viterbi-next.jl
    
    original_likelihood = calculate_read_likelihood(read, graph, k; graph_mode=graph_mode)
    
    # Option 1: Use Viterbi algorithm for full path optimization
    viterbi_result = try_viterbi_path_improvement(read, graph, k; graph_mode=graph_mode)
    if viterbi_result !== nothing
        viterbi_read, viterbi_likelihood = viterbi_result
        viterbi_improvement = viterbi_likelihood - original_likelihood
        
        # If Viterbi shows significant improvement, use it
        if viterbi_improvement > 0.01
            return viterbi_read, viterbi_improvement
        end
    end
    
    # Option 2: Statistical resampling for alternative paths
    statistical_result = try_statistical_path_resampling(read, graph, k; graph_mode=graph_mode)
    if statistical_result !== nothing
        stat_read, stat_likelihood = statistical_result
        stat_improvement = stat_likelihood - original_likelihood
        
        # If statistical resampling shows improvement, use it
        if stat_improvement > 0.01
            return stat_read, stat_improvement
        end
    end
    
    # Option 3: Local heuristic improvements (fallback)
    return try_local_path_improvements(read, graph, k, original_likelihood; graph_mode=graph_mode)
end

"""
Calculate likelihood of a FASTQ read given the current graph.
"""
function calculate_read_likelihood(read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol = :canonical)::Float64
    # Detect sequence type dynamically using existing alphabets.jl functions
    sequence_string = FASTX.sequence(String, read)
    alphabet = detect_alphabet(sequence_string)
    sequence_type = alphabet_to_biosequence_type(alphabet)
    
    sequence = extract_typed_sequence(read, sequence_type)
    quality_scores = collect(FASTX.quality_scores(read))
    return calculate_sequence_likelihood(sequence, quality_scores, graph, k; graph_mode=graph_mode)
end

"""
Calculate quality-aware likelihood of a sequence given the qualmer graph.
Uses both k-mer presence and quality-based confidence from qualmer observations.
"""
function calculate_sequence_likelihood(sequence::AbstractString, quality, graph, k::Int; graph_mode::Symbol = :canonical)::Float64
    alphabet = detect_alphabet(sequence)
    typed_sequence = detect_and_extract_sequence(sequence, alphabet)
    quality_scores = if quality isa AbstractString
        Int8.(Int.(collect(codeunits(quality))) .- 33)
    elseif quality isa AbstractVector{<:Integer}
        Int8.(quality)
    else
        throw(ArgumentError("Unsupported quality type for calculate_sequence_likelihood"))
    end
    return calculate_sequence_likelihood(typed_sequence, quality_scores, graph, k; graph_mode=graph_mode)
end

function _resolve_kmer_label(graph, kmer; graph_mode::Symbol = :canonical)
    if haskey(graph, kmer)
        return kmer
    end
    if graph_mode == :canonical && kmer isa Union{Kmers.DNAKmer, Kmers.RNAKmer}
        canonical_kmer = BioSequences.canonical(kmer)
        if haskey(graph, canonical_kmer)
            return canonical_kmer
        end
    end
    return nothing
end

function _resolve_qualmer_for_graph(graph, qmer::Qualmer; graph_mode::Symbol = :canonical)
    resolved_kmer = _resolve_kmer_label(graph, qmer.kmer; graph_mode=graph_mode)
    if resolved_kmer === nothing
        return nothing
    end
    if resolved_kmer == qmer.kmer
        return qmer
    end
    if graph_mode == :canonical && qmer.kmer isa Union{Kmers.DNAKmer, Kmers.RNAKmer}
        return canonical(qmer)
    end
    return Qualmer(resolved_kmer, qmer.qualities)
end

function calculate_sequence_likelihood(sequence::BioSequences.BioSequence, quality::Vector{Int8}, graph, k::Int; graph_mode::Symbol = :canonical)::Float64
    if length(sequence) < k
        return 0.0
    end
    
    total_log_likelihood = 0.0
    
    # Use Kmers.jl iterators for efficient k-mer extraction with quality awareness
    if sequence isa BioSequences.LongDNA
        for (pos, kmer) in enumerate(Kmers.FwDNAMers{k}(sequence))
            resolved_kmer = _resolve_kmer_label(graph, kmer; graph_mode=graph_mode)
            
            if resolved_kmer !== nothing
                vertex_data = graph[resolved_kmer]
                
                # Enhanced quality-aware likelihood calculation
                kmer_likelihood = calculate_qualmer_likelihood(
                    resolved_kmer,
                    quality[pos:pos+k-1],  # Quality scores for this k-mer
                    vertex_data
                )
                
                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                # Penalty for unseen k-mers, but consider quality
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:pos+k-1])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    elseif sequence isa BioSequences.LongRNA
        for (pos, kmer) in enumerate(Kmers.FwRNAMers{k}(sequence))
            resolved_kmer = _resolve_kmer_label(graph, kmer; graph_mode=graph_mode)
            
            if resolved_kmer !== nothing
                vertex_data = graph[resolved_kmer]
                
                kmer_likelihood = calculate_qualmer_likelihood(
                    resolved_kmer,
                    quality[pos:pos+k-1],
                    vertex_data
                )
                
                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:pos+k-1])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    else  # LongAA - no canonical form needed
        for (pos, kmer) in enumerate(Kmers.FwAAMers{k}(sequence))
            if haskey(graph, kmer)
                vertex_data = graph[kmer]
                
                kmer_likelihood = calculate_qualmer_likelihood(
                    kmer, 
                    quality[pos:pos+k-1],
                    vertex_data
                )
                
                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:pos+k-1])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    end
    
    return total_log_likelihood
end

function _rhizomorph_first_dataset_id(vertex_data)
    dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
    return isempty(dataset_ids) ? nothing : first(dataset_ids)
end

function _rhizomorph_joint_probability(vertex_data)::Float64
    dataset_id = _rhizomorph_first_dataset_id(vertex_data)
    if dataset_id === nothing
        return 0.0
    end

    joint_quality = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, dataset_id)
    if joint_quality === nothing
        return 0.0
    end

    log_prob = 0.0
    for q in joint_quality
        log_prob += log(max(phred_to_probability(q), 1e-10))
    end
    return exp(log_prob)
end

function _vertex_joint_probability(vertex_data)::Float64
    if hasproperty(vertex_data, :joint_probability)
        return vertex_data.joint_probability
    end
    return _rhizomorph_joint_probability(vertex_data)
end

function _vertex_mean_quality(vertex_data)::Float64
    if hasproperty(vertex_data, :mean_quality)
        return vertex_data.mean_quality
    end

    dataset_id = _rhizomorph_first_dataset_id(vertex_data)
    if dataset_id === nothing
        return 0.0
    end

    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, dataset_id)
    if mean_quality === nothing || isempty(mean_quality)
        return 0.0
    end

    return Statistics.mean(mean_quality)
end

function _vertex_coverage(vertex_data)::Int
    if hasproperty(vertex_data, :coverage)
        coverage = vertex_data.coverage
        return coverage isa Integer ? coverage : length(coverage)
    end
    return Int(Mycelia.Rhizomorph.count_evidence(vertex_data))
end

"""
Calculate likelihood of a k-mer given its observed quality scores and qualmer vertex data.
Leverages existing qualmer graph quality calculations.
"""
function calculate_qualmer_likelihood(kmer, quality_scores::AbstractVector{Int8}, vertex_data)::Float64
    # Base likelihood from graph joint probability
    base_likelihood = _vertex_joint_probability(vertex_data)
    
    # Quality-based confidence adjustment
    # Convert Int8 quality scores to UInt8 for compatibility with existing functions
    quality_uint8 = UInt8.(max.(0, quality_scores))
    
    # Calculate position-wise quality confidence
    quality_confidence = 1.0
    for q_score in quality_uint8
        prob_correct = phred_to_probability(q_score)
        quality_confidence *= prob_correct
    end
    
    # Combine graph confidence with observed quality confidence
    # Weight the graph probability by the quality confidence of this observation
    combined_likelihood = base_likelihood * quality_confidence + 
                         (1.0 - quality_confidence) * (base_likelihood * 0.1)  # Lower confidence for low quality
    
    return max(combined_likelihood, 1e-10)
end

"""
Calculate penalty for unseen k-mers based on their quality scores.
High quality unseen k-mers get less penalty than low quality ones.
"""
function calculate_unseen_kmer_penalty(quality_scores::AbstractVector{Int8})::Float64
    quality_uint8 = UInt8.(max.(0, quality_scores))
    
    # Calculate average quality confidence
    avg_quality_confidence = Statistics.mean(phred_to_probability(q) for q in quality_uint8)
    
    # High quality unseen k-mers might be real but missing from training data
    # Low quality unseen k-mers are likely errors
    base_penalty = 1e-6
    quality_adjusted_penalty = base_penalty * avg_quality_confidence + 
                              base_penalty * 0.01 * (1.0 - avg_quality_confidence)
    
    return quality_adjusted_penalty
end

"""
Generate alternative k-mers for improvement attempts using proper k-mer objects.
"""
function generate_kmer_alternatives(original_kmer::Kmers.Kmer, graph; graph_mode::Symbol = :canonical)
    sequence_type = if original_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif original_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else
        BioSequences.LongAA
    end
    return generate_kmer_alternatives(original_kmer, graph, sequence_type; graph_mode=graph_mode)
end

function generate_kmer_alternatives(original_kmer::AbstractString, graph; graph_mode::Symbol = :canonical)
    alphabet = detect_alphabet(original_kmer)
    sequence_type = alphabet_to_biosequence_type(alphabet)
    k = length(original_kmer)
    kmer_type = if sequence_type <: BioSequences.LongDNA
        Kmers.DNAKmer{k}
    elseif sequence_type <: BioSequences.LongRNA
        Kmers.RNAKmer{k}
    else
        Kmers.AAKmer{k}
    end
    return generate_kmer_alternatives(kmer_type(original_kmer), graph, sequence_type; graph_mode=graph_mode)
end

function generate_kmer_alternatives(original_kmer, graph, sequence_type::Type{<:BioSequences.BioSequence}; graph_mode::Symbol = :canonical)
    alternatives = []
    k = length(original_kmer)
    
    # Get alphabet based on sequence type
    if sequence_type <: BioSequences.LongDNA
        alphabet = [BioSequences.DNA_A, BioSequences.DNA_T, BioSequences.DNA_G, BioSequences.DNA_C]
        kmer_type = Kmers.DNAKmer{k}
    elseif sequence_type <: BioSequences.LongRNA
        alphabet = [BioSequences.RNA_A, BioSequences.RNA_U, BioSequences.RNA_G, BioSequences.RNA_C]
        kmer_type = Kmers.RNAKmer{k}
    else  # LongAA
        # For amino acids, use the standard 20 amino acids
        alphabet = [BioSequences.AA_A, BioSequences.AA_R, BioSequences.AA_N, BioSequences.AA_D, 
                   BioSequences.AA_C, BioSequences.AA_Q, BioSequences.AA_E, BioSequences.AA_G,
                   BioSequences.AA_H, BioSequences.AA_I, BioSequences.AA_L, BioSequences.AA_K,
                   BioSequences.AA_M, BioSequences.AA_F, BioSequences.AA_P, BioSequences.AA_S,
                   BioSequences.AA_T, BioSequences.AA_W, BioSequences.AA_Y, BioSequences.AA_V]
        kmer_type = Kmers.AAKmer{k}
    end
    
    # Try single position substitutions
    for i in 1:k
        original_symbol = original_kmer[i]
        for symbol in alphabet
            if symbol != original_symbol
                # Build alternative k-mer with substitution at position i
                alt_kmer_seq = sequence_type([original_kmer[j] for j in 1:k])
                alt_kmer_seq[i] = symbol
                alt_kmer = kmer_type(alt_kmer_seq)
                
                # Resolve k-mer labels against the current graph mode
                resolved_alt = _resolve_kmer_label(graph, alt_kmer; graph_mode=graph_mode)
                resolved_orig = _resolve_kmer_label(graph, original_kmer; graph_mode=graph_mode)

                # Only include if alternative exists in graph with higher probability
                if resolved_alt !== nothing && resolved_orig !== nothing
                    alt_data = graph[resolved_alt]
                    orig_data = graph[resolved_orig]

                    if _vertex_joint_probability(alt_data) > _vertex_joint_probability(orig_data)
                        push!(alternatives, resolved_alt)
                    end
                end
            end
        end
    end
    
    return alternatives
end

"""
Adjust quality scores based on likelihood improvement.
"""
function adjust_quality_scores(original_quality::String, improved_sequence::String, likelihood_improvement::Float64)::String
    # For now, return original quality scores
    # This can be enhanced to adjust qualities based on the improvement
    if length(improved_sequence) == length(original_quality)
        return original_quality
    else
        # Handle length changes by extending or truncating quality
        if length(improved_sequence) > length(original_quality)
            # Extend with median quality
            median_qual = isempty(original_quality) ? 'I' : original_quality[length(original_quality)÷2]
            return original_quality * repeat(string(median_qual), length(improved_sequence) - length(original_quality))
        else
            # Truncate quality
            return original_quality[1:length(improved_sequence)]
        end
    end
end

# =============================================================================
# Decision Making and Termination Conditions
# =============================================================================

"""
Determine if sufficient improvements were made to continue with current k.
Enhanced with convergence detection and adaptive thresholds.
"""
function sufficient_improvements(improvements_made::Int, total_reads::Int, threshold::Float64 = 0.05;
                                iteration_history::Vector{Dict{Symbol, Any}} = Dict{Symbol, Any}[],
                                min_improvement_trend::Int = 3)::Bool
    if total_reads == 0
        return false
    end
    
    improvement_rate = improvements_made / total_reads
    
    # Basic threshold check
    if improvement_rate < threshold
        return false
    end
    
    # Early convergence detection based on improvement trend
    if length(iteration_history) >= min_improvement_trend
        recent_rates = [hist[:improvement_rate] for hist in iteration_history[end-min_improvement_trend+1:end]]
        
        # Check for diminishing returns (decreasing trend)
        if length(recent_rates) >= 2
            trend_decreasing = all(recent_rates[i] >= recent_rates[i+1] for i in 1:length(recent_rates)-1)
            avg_recent_rate = Statistics.mean(recent_rates)
            
            # If trend is consistently decreasing and below adaptive threshold, stop
            if trend_decreasing && avg_recent_rate < threshold * 2.0
                return false
            end
        end
        
        # Check for plateau (very small improvements)
        if all(rate < threshold * 0.5 for rate in recent_rates)
            return false
        end
    end
    
    return true
end

"""
Enhanced convergence detection for k-mer progression.
Determines if we should move to the next k-mer size based on multiple criteria.
"""
function should_continue_k_progression(k_history::Vector{Dict{Symbol, Any}}, 
                                     current_k::Int, 
                                     max_k::Int;
                                     convergence_window::Int = 3,
                                     quality_improvement_threshold::Float64 = 0.001)::Bool
    
    if current_k >= max_k
        return false
    end
    
    if length(k_history) < convergence_window
        return true
    end
    
    # Check recent quality improvements across iterations
    recent_history = k_history[end-convergence_window+1:end]
    
    # Calculate average improvement metrics
    avg_improvement_rate = Statistics.mean([hist[:total_improvements] / hist[:total_reads] for hist in recent_history])
    avg_runtime_per_read = Statistics.mean([hist[:runtime_seconds] / hist[:total_reads] for hist in recent_history])
    
    # Convergence criteria
    improvements_declining = avg_improvement_rate < quality_improvement_threshold
    runtime_increasing = avg_runtime_per_read > 0.1  # More than 0.1 seconds per read indicates inefficiency
    
    # If improvements are minimal and runtime is high, move to next k
    if improvements_declining && runtime_increasing
        return false
    end
    
    # If we're making good progress, continue
    return avg_improvement_rate >= quality_improvement_threshold
end

"""
Finalize iterative assembly by combining results from all k-mer sizes and iterations.
"""
function finalize_iterative_assembly(output_dir::String, k_progression::Vector{Int}, 
                                    iteration_history::Dict{Int, Vector{Dict{Symbol, Any}}}, 
                                    total_runtime::Float64; verbose::Bool = true)
    if verbose
        println("Finalizing iterative assembly results...")
    end
    
    # Find the final output file (last iteration of largest k)
    final_k = maximum(k_progression)
    final_iteration_count = length(iteration_history[final_k])
    
    # Get the most recent timestamp from final k
    final_timestamp = iteration_history[final_k][end][:timestamp]
    final_fastq = joinpath(output_dir, "reads_k$(final_k)_iter$(final_iteration_count)_$(final_timestamp).fastq")
    
    # Calculate summary statistics
    total_iterations = sum(length(iterations) for iterations in values(iteration_history))
    total_improvements = sum(sum(iter[:improvements_made] for iter in iterations) for iterations in values(iteration_history))
    
    # Read final assembly for k-mer extraction
    if isfile(final_fastq)
        final_reads = collect(FASTX.FASTQ.Reader(open(final_fastq)))
        final_assembly = [FASTX.sequence(String, read) for read in final_reads]
    else
        final_assembly = String[]
    end
    
    # Create comprehensive metadata
    metadata = Dict(
        :total_runtime => total_runtime,
        :total_iterations => total_iterations,
        :total_improvements => total_improvements,
        :k_progression => k_progression,
        :final_k => final_k,
        :final_fastq_file => final_fastq,
        :output_directory => output_dir,
        :iteration_history => iteration_history,
        :assembly_type => "iterative_maximum_likelihood",
        :version => "Phase_5.2a"
    )
    
    if verbose
        println("Iterative assembly finalized:")
        println("  Total k-mer sizes: $(length(k_progression))")
        println("  Total iterations: $total_iterations")
        println("  Total improvements: $total_improvements")
        println("  Final FASTQ: $final_fastq")
    end
    
    return Dict(
        :final_assembly => final_assembly,
        :k_progression => k_progression,
        :metadata => metadata
    )
end

# =============================================================================
# Utility and Integration Functions
# =============================================================================

"""
Generate summary report of iterative assembly process.
"""
function iterative_assembly_summary(result::Dict)::String
    metadata = result[:metadata]
    k_prog = metadata[:k_progression]
    
    report = "Mycelia Iterative Assembly Summary\n"
    report *= "=" * "="^50 * "\n\n"
    
    report *= "Assembly Type: $(metadata[:assembly_type])\n"
    report *= "Version: $(metadata[:version])\n"
    report *= "Total Runtime: $(round(metadata[:total_runtime], digits=2)) seconds\n"
    report *= "K-mer Progression: $(sort(k_prog))\n"
    report *= "Total Iterations: $(metadata[:total_iterations])\n"
    report *= "Total Improvements: $(metadata[:total_improvements])\n"
    report *= "Final FASTQ: $(metadata[:final_fastq_file])\n\n"
    
    # Per-k statistics
    report *= "Per-K Statistics:\n"
    for k in sort(k_prog)
        iterations = metadata[:iteration_history][k]
        total_improvements_k = sum(iter[:improvements_made] for iter in iterations)
        report *= "  k=$k: $(length(iterations)) iterations, $total_improvements_k improvements\n"
    end
    
    return report
end

# =============================================================================
# Enhanced Statistical Path Improvement (Phase 5.2b)
# Integration with Viterbi algorithms and statistical path resampling
# =============================================================================

"""
Try to improve FASTQ read using Viterbi algorithm from viterbi-next.jl.
Returns (improved_read, likelihood) or nothing if no improvement.
"""
function try_viterbi_path_improvement(read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol = :canonical)::Union{Tuple{FASTX.FASTQ.Record, Float64}, Nothing}
    try
        # Detect sequence type dynamically
        sequence_string = FASTX.sequence(String, read)
        alphabet = detect_alphabet(sequence_string)
        sequence_type = alphabet_to_biosequence_type(alphabet)
        
        sequence = extract_typed_sequence(read, sequence_type)
        quality_scores = collect(FASTX.quality_scores(read))
        
        # Create ViterbiConfig for this sequence
        config = ViterbiConfig(
            match_prob = 0.95,
            mismatch_prob = 0.04,
            insertion_prob = 0.005,
            deletion_prob = 0.005,
            use_log_space = true,
            consider_reverse_complement = (alphabet != :AA)  # No reverse complement for amino acids
        )
        
        # Convert sequence to k-mer observations for Viterbi
        observations = String[]
        if length(sequence_string) >= k
            for i in 1:(length(sequence_string) - k + 1)
                push!(observations, sequence_string[i:(i + k - 1)])
            end
        end
        
        if !isempty(observations)
            # Use Viterbi algorithm to find optimal path
            viterbi_result = viterbi_decode_next(graph, observations, config)
            
            if viterbi_result !== nothing && !isempty(viterbi_result.polished_sequence)
                # Convert improved sequence back to proper BioSequence type
                improved_sequence = detect_and_extract_sequence(viterbi_result.polished_sequence, alphabet)
                
                # Calculate likelihood of improved sequence
                improved_likelihood = calculate_sequence_likelihood(improved_sequence, quality_scores, graph, k; graph_mode=graph_mode)
                
                # Create improved FASTQ record
                improved_quality = adjust_quality_scores(FASTX.quality(read), viterbi_result.polished_sequence, improved_likelihood)
                improved_record = FASTX.FASTQ.Record(FASTX.identifier(read), viterbi_result.polished_sequence, improved_quality)
                
                return (improved_record, improved_likelihood)
            end
        end
        
    catch e
        # If Viterbi fails, return nothing to fall back to other methods
        return nothing
    end
    
    return nothing
end

function try_viterbi_path_improvement(sequence::AbstractString, quality, graph, k::Int; graph_mode::Symbol = :canonical)
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_viterbi_path_improvement"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    result = try_viterbi_path_improvement(record, graph, k; graph_mode=graph_mode)
    if result === nothing
        return nothing
    end

    improved_record, improvement = result
    return FASTX.sequence(String, improved_record), improvement
end

"""
Try statistical path resampling for alternative high-likelihood paths.
Returns (improved_read, likelihood) or nothing if no improvement.
"""
function try_statistical_path_resampling(read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol = :canonical)::Union{Tuple{FASTX.FASTQ.Record, Float64}, Nothing}
    try
        # Detect sequence type dynamically
        sequence_string = FASTX.sequence(String, read)
        alphabet = detect_alphabet(sequence_string)
        sequence_type = alphabet_to_biosequence_type(alphabet)
        
        sequence = extract_typed_sequence(read, sequence_type)
        quality_scores = collect(FASTX.quality_scores(read))
        quality_string = FASTX.quality(read)
        
        # Generate multiple alternative paths using qualmer sampling
        alternative_paths = generate_alternative_qualmer_paths(read, graph, k, num_samples=5, graph_mode=graph_mode)
        
        best_sequence = sequence
        best_likelihood = calculate_sequence_likelihood(sequence, quality_scores, graph, k; graph_mode=graph_mode)
        best_quality_scores = quality_scores
        
        # Evaluate each alternative path
        for alt_path in alternative_paths
            if !isempty(alt_path)
                # Convert qualmer path to proper BioSequence + quality
                alt_result = qualmer_path_to_biosequence(alt_path)
                alt_sequence = alt_result.sequence
                alt_quality_scores = alt_result.quality_scores
                
                alt_likelihood = calculate_sequence_likelihood(alt_sequence, alt_quality_scores, graph, k; graph_mode=graph_mode)
                
                if alt_likelihood > best_likelihood
                    best_sequence = alt_sequence
                    best_likelihood = alt_likelihood
                    best_quality_scores = alt_quality_scores
                end
            end
        end
        
        # Return improvement if found
        if best_sequence != sequence
            # Create improved FASTQ record using proper types
            improved_sequence_string = string(best_sequence)
            improved_quality_string = String([Char(Int(q) + 33) for q in best_quality_scores])  # Convert back to ASCII
            improved_record = FASTX.FASTQ.Record(FASTX.identifier(read), improved_sequence_string, improved_quality_string)
            return (improved_record, best_likelihood)
        end
        
    catch e
        # If statistical resampling fails, return nothing
        return nothing
    end
    
    return nothing
end

function try_statistical_path_resampling(sequence::AbstractString, quality, graph, k::Int; graph_mode::Symbol = :canonical)
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_statistical_path_resampling"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    result = try_statistical_path_resampling(record, graph, k; graph_mode=graph_mode)
    if result === nothing
        return nothing
    end

    improved_record, improvement = result
    return FASTX.sequence(String, improved_record), improvement
end

"""
Local heuristic improvements (fallback method).
Returns (improved_read, likelihood_improvement).
"""
function try_local_path_improvements(read::FASTX.FASTQ.Record, graph, k::Int, original_likelihood::Float64; graph_mode::Symbol = :canonical)::Tuple{FASTX.FASTQ.Record, Float64}
    sequence_string = FASTX.sequence(String, read)
    alphabet = detect_alphabet(sequence_string)
    sequence_type = alphabet_to_biosequence_type(alphabet)
    quality_scores = collect(FASTX.quality_scores(read))

    improved_seq = sequence_string  # Start with original
    best_likelihood = original_likelihood
    
    # Try local improvements at positions with low-probability k-mers
    for i in 1:(length(sequence_string) - k + 1)
        kmer_str = sequence_string[i:i+k-1]

        raw_kmer = if sequence_type <: BioSequences.LongDNA
            Kmers.DNAKmer{k}(kmer_str)
        elseif sequence_type <: BioSequences.LongRNA
            Kmers.RNAKmer{k}(kmer_str)
        else
            Kmers.AAKmer{k}(kmer_str)
        end

        target_kmer = _resolve_kmer_label(graph, raw_kmer; graph_mode=graph_mode)
        if target_kmer !== nothing
            vertex_data = graph[target_kmer]
            if _vertex_joint_probability(vertex_data) < 0.7  # Low confidence k-mer
                # Try improving this position
                alternatives = generate_kmer_alternatives(target_kmer, graph, sequence_type; graph_mode=graph_mode)
                for alt_kmer in alternatives
                    alt_kmer_str = string(alt_kmer)
                    if i == 1
                        candidate_seq = alt_kmer_str * sequence_string[i+k:end]
                    elseif i+k-1 == length(sequence_string)
                        candidate_seq = sequence_string[1:i-1] * alt_kmer_str
                    else
                        candidate_seq = sequence_string[1:i-1] * alt_kmer_str * sequence_string[i+k:end]
                    end

                    candidate_bioseq = detect_and_extract_sequence(candidate_seq, alphabet)
                    candidate_likelihood = calculate_sequence_likelihood(candidate_bioseq, quality_scores, graph, k; graph_mode=graph_mode)

                    if candidate_likelihood > best_likelihood
                        improved_seq = candidate_seq
                        best_likelihood = candidate_likelihood
                    end
                end
            end
        end
    end
    
    likelihood_improvement = best_likelihood - original_likelihood
    improved_quality = adjust_quality_scores(FASTX.quality(read), improved_seq, likelihood_improvement)
    improved_record = FASTX.FASTQ.Record(FASTX.identifier(read), improved_seq, improved_quality)
    return improved_record, likelihood_improvement
end

function try_local_path_improvements(sequence::AbstractString, quality, graph, k::Int, original_likelihood::Float64; graph_mode::Symbol = :canonical)::Tuple{String, Float64}
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_local_path_improvements"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    improved_record, improvement = try_local_path_improvements(record, graph, k, original_likelihood; graph_mode=graph_mode)
    return FASTX.sequence(String, improved_record), improvement
end

"""
Generate alternative qualmer paths through the graph using quality-aware probabilistic sampling.
"""
function generate_alternative_qualmer_paths(record::FASTX.FASTQ.Record, graph, k::Int; num_samples::Int = 5, graph_mode::Symbol = :canonical)::Vector{Vector{<:Qualmer}}
    paths = Vector{Vector{Qualmer}}()
    
    try
        # Generate qualmers using existing qualmer-analysis functions
        original_qualmers = collect(qualmers_unambiguous(record, k))
        if isempty(original_qualmers)
            return Vector{Vector{Qualmer}}()
        end

        original_path = Qualmer[]
        for (qmer, pos) in original_qualmers
            resolved = _resolve_qualmer_for_graph(graph, qmer; graph_mode=graph_mode)
            if resolved !== nothing
                push!(original_path, resolved)
            end
        end
        if isempty(original_path)
            return Vector{Vector{Qualmer}}()
        end
        
        # For each sample, create alternative path using quality-aware sampling
        for sample in 1:num_samples
            alternative_path = copy(original_path)
            
            # Select positions to modify based on quality confidence
            low_confidence_positions = Int[]
            for (i, qmer) in enumerate(original_path)
                current_kmer = qmer.kmer
                if haskey(graph, current_kmer)
                    vertex_data = graph[current_kmer]
                    # Prefer modifying low-confidence positions
                    if _vertex_joint_probability(vertex_data) < 0.8 || _vertex_mean_quality(vertex_data) < 30.0
                        push!(low_confidence_positions, i)
                    end
                end
            end
            
            # If no low confidence positions, randomly select some
            modification_positions = if !isempty(low_confidence_positions)
                Random.shuffle(low_confidence_positions)[1:min(2, length(low_confidence_positions))]
            else
                Random.shuffle(1:length(original_path))[1:min(2, length(original_path))]
            end
            
            for pos in modification_positions
                if pos <= length(alternative_path)
                    current_qualmer = alternative_path[pos]
                    current_kmer = current_qualmer.kmer

                    if haskey(graph, current_kmer)
                        # Get high-quality neighbor k-mers for potential transitions
                        vertex_data = graph[current_kmer]
                        vertex_coverage = _vertex_coverage(vertex_data)
                        all_vertices = collect(MetaGraphsNext.labels(graph))
                        
                        # Filter neighbors that could be valid transitions
                        valid_neighbors = []
                        neighbor_weights = Float64[]
                        
                        for neighbor_kmer in all_vertices
                            neighbor_data = graph[neighbor_kmer]
                            
                            # Quality-aware weighting combining multiple factors
                            weight = _vertex_joint_probability(neighbor_data) * 
                                    (_vertex_mean_quality(neighbor_data) / 60.0) *  # Normalize quality (max ~60)
                                    (_vertex_coverage(neighbor_data) / max(1, vertex_coverage))  # Relative coverage
                            
                            if weight > 0.01  # Minimum threshold for consideration
                                push!(valid_neighbors, neighbor_kmer)
                                push!(neighbor_weights, weight)
                            end
                        end
                        
                        if !isempty(valid_neighbors)
                            # Sample neighbor based on quality-aware weights
                            neighbor_weights ./= sum(neighbor_weights)  # Normalize
                            selected_idx = StatsBase.sample(1:length(valid_neighbors), StatsBase.Weights(neighbor_weights))
                            
                            # Create new qualmer using the selected k-mer and quality from original
                            selected_kmer = valid_neighbors[selected_idx]
                            original_qualities = current_qualmer.qualities
                            
                            new_qualmer = Qualmer(selected_kmer, original_qualities)
                            alternative_path[pos] = new_qualmer
                        end
                    end
                end
            end
            
            push!(paths, alternative_path)
        end
        
    catch e
        # If path generation fails, return empty paths
        @debug "Alternative path generation failed: $e"
        return Vector{Vector{Qualmer}}()
    end
    
    return paths
end

"""
Convert sequence with quality to qualmer path representation using graph vertices.
"""
function sequence_to_qualmer_path(sequence::String, quality::String, graph, k::Int; graph_mode::Symbol = :canonical)::Vector{<:Qualmer}
    path = Qualmer[]

    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels) || length(sequence) < k
        return path
    end

    kmer_type = typeof(first(labels))

    for i in 1:(length(sequence) - k + 1)
        kmer_str = sequence[i:i+k-1]
        qual_segment = quality[i:i+k-1]
        quality_scores = UInt8.(codeunits(qual_segment)) .- UInt8(33)

        qmer = Qualmer(kmer_type(kmer_str), quality_scores)
        resolved = _resolve_qualmer_for_graph(graph, qmer; graph_mode=graph_mode)
        if resolved !== nothing
            push!(path, resolved)
        end
    end

    return path
end

"""
Convert sequence to k-mer path representation (simplified version for backward compatibility).
"""
function sequence_to_kmer_path(sequence::String, k::Int)::Vector{String}
    path = String[]
    
    for i in 1:(length(sequence) - k + 1)
        push!(path, sequence[i:i+k-1])
    end
    
    return path
end

"""
Convert qualmer path back to BioSequence and quality vector.
Returns a named tuple with sequence and quality_scores.
"""
function qualmer_path_to_biosequence(path::Vector{<:Qualmer})
    if isempty(path)
        return (sequence=BioSequences.LongDNA{4}(""), quality_scores=Int8[])
    end
    
    # Determine sequence type from first qualmer
    first_kmer = path[1].kmer
    sequence_type = if first_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif first_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else  # AAKmer
        BioSequences.LongAA
    end
    
    # Start with first k-mer
    sequence_symbols = collect(path[1].kmer)
    quality_scores = collect(Int8.(path[1].qualities))
    
    # Add last symbol of each subsequent k-mer
    for i in 2:length(path)
        push!(sequence_symbols, path[i].kmer[end])
        push!(quality_scores, Int8(path[i].qualities[end]))
    end
    
    # Construct the BioSequence
    final_sequence = sequence_type(sequence_symbols)
    
    return (sequence=final_sequence, quality_scores=quality_scores)
end

"""
Convert k-mer path back to BioSequence.
"""
function kmer_path_to_biosequence(path::Vector{<:Kmers.Kmer})
    if isempty(path)
        return BioSequences.LongDNA{4}("")
    end
    
    # Determine sequence type from first k-mer
    first_kmer = path[1]
    sequence_type = if first_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif first_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else  # AAKmer
        BioSequences.LongAA
    end
    
    # Start with first k-mer
    sequence_symbols = collect(path[1])
    
    # Add last symbol of each subsequent k-mer
    for i in 2:length(path)
        push!(sequence_symbols, path[i][end])
    end
    
    # Construct the BioSequence
    return sequence_type(sequence_symbols)
end
