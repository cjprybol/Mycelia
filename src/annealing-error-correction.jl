# # Annealing-based error correction for k-mer graphs
# # Based on analysis from Mycelia-Dev notebooks (2021-09-01 through 2021-09-11)

# """
#     AnnealingCorrectionResult

# Results from annealing-based error correction analysis.

# # Fields
# - `corrected_sequence::String`: Error-corrected sequence
# - `original_sequence::String`: Original input sequence
# - `correction_positions::Vector{Int}`: Positions where corrections were made
# - `confidence_scores::Vector{Float64}`: Confidence scores for each correction
# - `bubbles_resolved::Int`: Number of error bubbles successfully resolved
# - `iterations_performed::Int`: Number of annealing iterations completed
# """
# struct AnnealingCorrectionResult
#     corrected_sequence::String
#     original_sequence::String
#     correction_positions::Vector{Int}
#     confidence_scores::Vector{Float64}
#     bubbles_resolved::Int
#     iterations_performed::Int
# end

# """
#     BubbleDetectionResult

# Results from error bubble detection in k-mer graphs.

# # Fields
# - `bubble_regions::Vector{Tuple{Int, Int}}`: Start and end positions of detected bubbles
# - `bubble_types::Vector{Symbol}`: Types of bubbles (:insertion, :deletion, :substitution, :complex)
# - `solid_anchors::Vector{Tuple{String, String}}`: Solid k-mers bounding each bubble
# - `estimated_errors::Vector{Int}`: Estimated number of errors in each bubble
# """
# struct BubbleDetectionResult
#     bubble_regions::Vector{Tuple{Int, Int}}
#     bubble_types::Vector{Symbol}
#     solid_anchors::Vector{Tuple{String, String}}
#     estimated_errors::Vector{Int}
# end

# """
#     annealing_error_correction(observed_sequence::String,
#                               reference_kmers::Dict{String, Int},
#                               k::Int;
#                               max_iterations::Int=5,
#                               walk_attempts::Int=10,
#                               confidence_threshold::Float64=0.7) -> AnnealingCorrectionResult

# Perform annealing-based error correction using probabilistic k-mer graph traversal.

# This algorithm identifies error regions ("bubbles") in sequences by finding regions bounded 
# by solid k-mers that don't exist in the reference graph, then corrects these regions using 
# probabilistic walks through the reference graph.

# # Arguments
# - `observed_sequence::String`: The sequence to correct
# - `reference_kmers::Dict{String, Int}`: Reference k-mer counts from error-free sequences
# - `k::Int`: K-mer size used for analysis
# - `max_iterations::Int=5`: Maximum number of correction iterations
# - `walk_attempts::Int=10`: Number of probabilistic walks to attempt per bubble
# - `confidence_threshold::Float64=0.7`: Minimum confidence required for correction

# # Returns
# `AnnealingCorrectionResult`: Comprehensive correction results with quality metrics

# # Algorithm Details
# 1. **Bubble Detection**: Identifies regions where consecutive k-mers are missing from reference
# 2. **Edge Probability Calculation**: Computes traversal probabilities based on k-mer counts
# 3. **Probabilistic Walking**: Performs multiple random walks to find correction paths
# 4. **Consensus Selection**: Chooses most frequent path from multiple walk attempts
# 5. **Iterative Refinement**: Repeats process until no more bubbles are found

# # Mathematical Foundation
# Edge probabilities use two-stage weighting:
# ```
# p(edge) = (edge_count / Σ_edge_counts) × (dest_count / Σ_dest_counts)
# ```

# # Examples
# ```julia
# # Basic error correction
# reference = count_kmers(reference_sequences, 7)
# result = annealing_error_correction(observed_seq, reference, 7)

# # High-stringency correction
# result = annealing_error_correction(observed_seq, reference, 7,
#                                    confidence_threshold=0.9,
#                                    walk_attempts=20)

# # Check correction quality
# if result.confidence_scores |> mean > 0.8
#     corrected_seq = result.corrected_sequence
# end
# ```

# # References
# Based on annealing correction algorithms from Mycelia-Dev notebooks:
# - `2021-09-01-annealing-correction.ipynb`
# - `2021-09-05-annealing-correction.ipynb` 
# - `2021-09-06-annealing-correction.ipynb`
# - `2021-09-11-annealing-correction-L10-K5.ipynb`
# - `2021-09-11-annealing-correction-L10-K7.ipynb`
# - `2021-09-11-annealing-correction-L100-K7.ipynb`
# """
# function annealing_error_correction(observed_sequence::String,
#                                    reference_kmers::Dict{String, Int},
#                                    k::Int;
#                                    max_iterations::Int=5,
#                                    walk_attempts::Int=10,
#                                    confidence_threshold::Float64=0.7)
    
#     if length(observed_sequence) < k
#         throw(ArgumentError("Sequence length must be ≥ k-mer size"))
#     end
    
#     if isempty(reference_kmers)
#         throw(ArgumentError("Reference k-mers dictionary cannot be empty"))
#     end
    
#     original_sequence = observed_sequence
#     corrected_sequence = observed_sequence
#     all_correction_positions = Int[]
#     all_confidence_scores = Float64[]
#     total_bubbles_resolved = 0
    
#     @info "Starting annealing correction: $(length(observed_sequence))bp sequence, k=$k"
    
#     ## Iterative correction loop
#     for iteration in 1:max_iterations
        
#         ## Detect error bubbles in current sequence
#         bubble_result = detect_error_bubbles(corrected_sequence, reference_kmers, k)
        
#         if isempty(bubble_result.bubble_regions)
#             @info "No more error bubbles detected after $iteration iterations"
#             break
#         end
        
#         @debug "Iteration $iteration: found $(length(bubble_result.bubble_regions)) bubbles"
        
#         ## Correct each bubble
#         iteration_corrections = 0
        
#         for (i, (start_pos, end_pos)) in enumerate(bubble_result.bubble_regions)
#             bubble_type = bubble_result.bubble_types[i]
#             (opening_kmer, closing_kmer) = bubble_result.solid_anchors[i]
            
#             @debug "Correcting $bubble_type bubble at positions $start_pos-$end_pos"
            
#             ## Perform probabilistic correction
#             correction_result = correct_single_bubble(
#                 corrected_sequence, start_pos, end_pos,
#                 opening_kmer, closing_kmer,
#                 reference_kmers, k, walk_attempts
#             )
            
#             if correction_result.success && correction_result.confidence >= confidence_threshold
#                 ## Apply correction
#                 corrected_sequence = correction_result.corrected_sequence
                
#                 ## Track correction details
#                 push!(all_correction_positions, start_pos)
#                 push!(all_confidence_scores, correction_result.confidence)
#                 iteration_corrections += 1
                
#                 @debug "Applied correction with confidence $(round(correction_result.confidence, digits=3))"
#             else
#                 @debug "Skipped low-confidence correction ($(round(correction_result.confidence, digits=3)) < $confidence_threshold)"
#             end
#         end
        
#         total_bubbles_resolved += iteration_corrections
        
#         if iteration_corrections == 0
#             @info "No corrections applied in iteration $iteration - stopping"
#             break
#         end
#     end
    
#     @info "Annealing correction complete: $total_bubbles_resolved bubbles resolved, " *
#           "$(length(all_correction_positions)) total corrections"
    
#     return AnnealingCorrectionResult(
#         corrected_sequence,
#         original_sequence,
#         all_correction_positions,
#         all_confidence_scores,
#         total_bubbles_resolved,
#         min(max_iterations, total_bubbles_resolved > 0 ? findlast(x -> x > 0, [total_bubbles_resolved]) : 1)
#     )
# end

# """
#     detect_error_bubbles(sequence::String, reference_kmers::Dict{String, Int}, k::Int) -> BubbleDetectionResult

# Detect error bubbles in sequences by identifying regions bounded by solid k-mers.

# An error bubble is a region where consecutive k-mers are missing from the reference
# while being bounded by k-mers that exist in the reference (solid k-mers).

# # Arguments
# - `sequence::String`: Sequence to analyze for bubbles
# - `reference_kmers::Dict{String, Int}`: Reference k-mer counts
# - `k::Int`: K-mer size

# # Returns
# `BubbleDetectionResult`: Detailed information about detected bubbles
# """
# function detect_error_bubbles(sequence::String, reference_kmers::Dict{String, Int}, k::Int)
    
#     if length(sequence) < k
#         return BubbleDetectionResult(Tuple{Int,Int}[], Symbol[], Tuple{String,String}[], Int[])
#     end
    
#     ## Extract k-mers with positions
#     kmers_with_positions = Tuple{String, Int}[]
#     for i in 1:(length(sequence) - k + 1)
#         kmer = sequence[i:(i + k - 1)]
#         canonical_kmer = _get_canonical_kmer(kmer)
#         push!(kmers_with_positions, (canonical_kmer, i))
#     end
    
#     ## Identify solid vs missing k-mers
#     kmer_status = [haskey(reference_kmers, kmer) for (kmer, pos) in kmers_with_positions]
    
#     ## Find bubble regions
#     bubble_regions = Tuple{Int, Int}[]
#     bubble_types = Symbol[]
#     solid_anchors = Tuple{String, String}[]
#     estimated_errors = Int[]
    
#     i = 1
#     while i <= length(kmer_status)
#         if !kmer_status[i]  ## Found missing k-mer - potential bubble start
            
#             ## Find opening solid k-mer (last solid k-mer before this region)
#             opening_pos = i - 1
#             opening_kmer = opening_pos > 0 ? kmers_with_positions[opening_pos][1] : ""
            
#             ## Find extent of missing region
#             bubble_start = i
#             while i <= length(kmer_status) && !kmer_status[i]
#                 i += 1
#             end
#             bubble_end = i - 1
            
#             ## Find closing solid k-mer (first solid k-mer after this region)
#             closing_pos = i
#             closing_kmer = closing_pos <= length(kmers_with_positions) ? kmers_with_positions[closing_pos][1] : ""
            
#             ## Only process bubbles with valid solid anchors
#             if opening_kmer != "" && closing_kmer != ""
#                 sequence_start = kmers_with_positions[bubble_start][2]
#                 sequence_end = bubble_end <= length(kmers_with_positions) ? 
#                               kmers_with_positions[bubble_end][2] + k - 1 : 
#                               length(sequence)
                
#                 push!(bubble_regions, (sequence_start, sequence_end))
#                 push!(solid_anchors, (opening_kmer, closing_kmer))
                
#                 ## Estimate bubble type and error count
#                 bubble_length = sequence_end - sequence_start + 1
#                 missing_kmer_count = bubble_end - bubble_start + 1
                
#                 bubble_type = _classify_bubble_type(bubble_length, missing_kmer_count, k)
#                 push!(bubble_types, bubble_type)
#                 push!(estimated_errors, max(1, missing_kmer_count ÷ 2))  ## Conservative estimate
#             end
#         else
#             i += 1
#         end
#     end
    
#     return BubbleDetectionResult(bubble_regions, bubble_types, solid_anchors, estimated_errors)
# end

# """
#     correct_single_bubble(sequence::String, start_pos::Int, end_pos::Int,
#                          opening_kmer::String, closing_kmer::String,
#                          reference_kmers::Dict{String, Int}, k::Int,
#                          walk_attempts::Int=10)

# Correct a single error bubble using probabilistic graph traversal.
# """
# function correct_single_bubble(sequence::String, start_pos::Int, end_pos::Int,
#                               opening_kmer::String, closing_kmer::String,
#                               reference_kmers::Dict{String, Int}, k::Int,
#                               walk_attempts::Int=10)
    
#     ## Build edge probability graph from reference k-mers
#     edge_probabilities = build_kmer_edge_probabilities(reference_kmers, k)
    
#     ## Perform multiple probabilistic walks
#     correction_attempts = String[]
    
#     for attempt in 1:walk_attempts
#         try
#             path = probabilistic_walk(opening_kmer, closing_kmer, edge_probabilities, k, 20)  ## Max 20 steps
            
#             if !isempty(path)
#                 ## Reconstruct sequence from k-mer path
#                 corrected_region = reconstruct_sequence_from_path(path, k)
#                 push!(correction_attempts, corrected_region)
#             end
#         catch e
#             @debug "Walk attempt $attempt failed: $e"
#             continue
#         end
#     end
    
#     if isempty(correction_attempts)
#         return (success=false, confidence=0.0, corrected_sequence=sequence)
#     end
    
#     ## Select most frequent correction (consensus)
#     correction_counts = Dict{String, Int}()
#     for correction in correction_attempts
#         correction_counts[correction] = get(correction_counts, correction, 0) + 1
#     end
    
#     best_correction = ""
#     best_count = 0
#     for (correction, count) in correction_counts
#         if count > best_count
#             best_count = count
#             best_correction = correction
#         end
#     end
    
#     ## Calculate confidence as fraction of walks agreeing on best correction
#     confidence = best_count / length(correction_attempts)
    
#     ## Apply correction to sequence
#     corrected_sequence = sequence[1:(start_pos-1)] * best_correction * sequence[(end_pos+1):end]
    
#     return (success=true, confidence=confidence, corrected_sequence=corrected_sequence)
# end

# ## Helper function to build edge probabilities from k-mer counts
# function build_kmer_edge_probabilities(reference_kmers::Dict{String, Int}, k::Int)
#     edge_probs = Dict{String, Dict{Char, Float64}}()
    
#     ## For each k-mer, calculate probabilities of extending with each base
#     for (kmer, count) in reference_kmers
#         if length(kmer) != k
#             continue
#         end
        
#         suffix = kmer[2:end]  ## (k-1)-mer suffix
#         base_counts = Dict{Char, Float64}()
        
#         ## Find all k-mers with this suffix as prefix
#         for base in ['A', 'C', 'G', 'T']
#             next_kmer = suffix * base
#             canonical_next = _get_canonical_kmer(next_kmer)
            
#             if haskey(reference_kmers, canonical_next)
#                 base_counts[base] = Float64(reference_kmers[canonical_next])
#             end
#         end
        
#         ## Normalize to probabilities
#         total_count = sum(values(base_counts))
#         if total_count > 0
#             for base in keys(base_counts)
#                 base_counts[base] /= total_count
#             end
#             edge_probs[kmer] = base_counts
#         end
#     end
    
#     return edge_probs
# end

# ## Helper function for probabilistic walk through k-mer graph
# function probabilistic_walk(start_kmer::String, target_kmer::String,
#                            edge_probabilities::Dict{String, Dict{Char, Float64}},
#                            k::Int, max_steps::Int)
    
#     path = [start_kmer]
#     current_kmer = start_kmer
    
#     for step in 1:max_steps
#         if current_kmer == target_kmer
#             return path
#         end
        
#         ## Get possible extensions
#         if !haskey(edge_probabilities, current_kmer)
#             break  ## Dead end
#         end
        
#         base_probs = edge_probabilities[current_kmer]
#         if isempty(base_probs)
#             break
#         end
        
#         ## Sample next base according to probabilities
#         bases = collect(keys(base_probs))
#         probs = [base_probs[base] for base in bases]
        
#         ## Weighted random selection
#         cumulative_probs = cumsum(probs)
#         rand_val = rand()
#         selected_idx = findfirst(x -> x >= rand_val, cumulative_probs)
        
#         if selected_idx === nothing
#             selected_idx = length(bases)
#         end
        
#         selected_base = bases[selected_idx]
        
#         ## Construct next k-mer
#         next_kmer = current_kmer[2:end] * selected_base
#         push!(path, next_kmer)
#         current_kmer = next_kmer
#     end
    
#     return String[]  ## Failed to reach target
# end

# ## Helper function to reconstruct sequence from k-mer path
# function reconstruct_sequence_from_path(path::Vector{String}, k::Int)
#     if isempty(path)
#         return ""
#     end
    
#     if length(path) == 1
#         return path[1]
#     end
    
#     ## Start with first k-mer, then add last base of each subsequent k-mer
#     sequence = path[1]
#     for i in 2:length(path)
#         sequence *= path[i][end]
#     end
    
#     return sequence
# end

# ## Helper function for canonical k-mer representation
# function _get_canonical_kmer(kmer::String)
#     reverse_complement = _reverse_complement(kmer)
#     return kmer <= reverse_complement ? kmer : reverse_complement
# end

# ## Helper function for reverse complement
# function _reverse_complement(seq::String)
#     complement_map = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')
#     return reverse(String([complement_map[c] for c in seq]))
# end

# ## Helper function to classify bubble types
# function _classify_bubble_type(bubble_length::Int, missing_kmer_count::Int, k::Int)
#     if missing_kmer_count == 1
#         return :substitution
#     elseif bubble_length < k
#         return :deletion
#     elseif missing_kmer_count > bubble_length - k + 1
#         return :insertion
#     else
#         return :complex
#     end
# end