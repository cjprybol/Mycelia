# """
# Viterbi Algorithm Migration to Next-Generation Framework

# This module provides an enhanced, strand-aware Viterbi implementation using MetaGraphsNext.jl
# with improvements for batch processing, memory efficiency, and biological accuracy.

# Migrates from the legacy implementation in viterbi-polishing-and-error-correction.jl
# to work with the new strand-aware k-mer graph architecture.
# """

# """
#     ViterbiState

# State information for Viterbi algorithm on k-mer graphs.
# """
# struct ViterbiState
#     vertex_label::String
#     strand::StrandOrientation
#     emission_prob::Float64
#     position::Int
    
#     function ViterbiState(vertex_label::String, strand::StrandOrientation, 
#                          emission_prob::Float64, position::Int)
#         @assert 0.0 <= emission_prob <= 1.0 "Emission probability must be in [0,1]"
#         @assert position >= 0 "Position must be non-negative"
#         new(vertex_label, strand, emission_prob, position)
#     end
# end

# """
#     ViterbiPath

# Complete Viterbi path through the k-mer graph.
# """
# struct ViterbiPath
#     states::Vector{ViterbiState}
#     log_probability::Float64
#     polished_sequence::String
#     corrections_made::Vector{Tuple{Int, String, String}}  # (position, original, corrected)
    
#     function ViterbiPath(states::Vector{ViterbiState}, log_probability::Float64,
#                         polished_sequence::String, corrections_made::Vector{Tuple{Int, String, String}})
#         new(states, log_probability, polished_sequence, corrections_made)
#     end
# end

# """
#     ViterbiConfig

# Configuration for Viterbi algorithm parameters.
# """
# struct ViterbiConfig
#     # Emission probabilities
#     match_prob::Float64
#     mismatch_prob::Float64
#     insertion_prob::Float64
#     deletion_prob::Float64
    
#     # Transition parameters
#     stay_prob::Float64          # Probability of staying in same state
#     move_prob::Float64          # Probability of moving to next state
    
#     # Processing parameters
#     batch_size::Int             # For batch processing
#     memory_limit::Int           # Memory limit in MB
#     use_log_space::Bool         # Use log probabilities for numerical stability
    
#     # Strand handling
#     consider_reverse_complement::Bool
#     strand_switch_penalty::Float64
    
#     function ViterbiConfig(;
#         match_prob::Float64 = 0.95,
#         mismatch_prob::Float64 = 0.04,
#         insertion_prob::Float64 = 0.005,
#         deletion_prob::Float64 = 0.005,
#         stay_prob::Float64 = 0.9,
#         move_prob::Float64 = 0.1,
#         batch_size::Int = 1000,
#         memory_limit::Int = 1024,  # 1GB default
#         use_log_space::Bool = true,
#         consider_reverse_complement::Bool = true,
#         strand_switch_penalty::Float64 = 0.1
#     )
#         # Validate probabilities sum correctly
#         @assert abs((match_prob + mismatch_prob + insertion_prob + deletion_prob) - 1.0) < 1e-10
#         @assert abs((stay_prob + move_prob) - 1.0) < 1e-10
#         @assert all(p -> 0.0 <= p <= 1.0, [match_prob, mismatch_prob, insertion_prob, deletion_prob, stay_prob, move_prob])
        
#         new(match_prob, mismatch_prob, insertion_prob, deletion_prob,
#             stay_prob, move_prob, batch_size, memory_limit, use_log_space,
#             consider_reverse_complement, strand_switch_penalty)
#     end
# end

# """
#     create_hmm_from_graph(graph::MetaGraph, config::ViterbiConfig) -> (states, transitions, emissions)

# Create Hidden Markov Model parameters from a k-mer graph structure.
# """
# function create_hmm_from_graph(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, String, KmerVertexData, KmerEdgeData}, 
#                               config::ViterbiConfig)
#     vertices = collect(MetaGraphsNext.labels(graph))
#     n_states = length(vertices) * (config.consider_reverse_complement ? 2 : 1)
    
#     # State mapping: vertex_label -> (forward_index, reverse_index)
#     state_map = Dict{String, Tuple{Int, Int}}()
#     states = Vector{ViterbiState}()
    
#     idx = 1
#     for vertex_label in vertices
#         forward_idx = idx
#         reverse_idx = config.consider_reverse_complement ? idx + 1 : idx
        
#         # Forward strand state
#         push!(states, ViterbiState(vertex_label, Forward, 1.0, 0))
#         idx += 1
        
#         # Reverse strand state (if enabled)
#         if config.consider_reverse_complement
#             push!(states, ViterbiState(vertex_label, Reverse, 1.0, 0))
#             idx += 1
#         end
        
#         state_map[vertex_label] = (forward_idx, reverse_idx)
#     end
    
#     # Initialize transition matrix
#     transitions = zeros(Float64, n_states, n_states)
    
#     # Fill transition probabilities based on graph edges
#     for edge in MetaGraphsNext.edge_labels(graph)
#         src_label, dst_label = edge
#         edge_data = graph[src_label, dst_label]
        
#         src_forward, src_reverse = state_map[src_label]
#         dst_forward, dst_reverse = state_map[dst_label]
        
#         # Set transitions based on edge strand compatibility
#         if edge_data.src_strand == Forward && edge_data.dst_strand == Forward
#             transitions[src_forward, dst_forward] = edge_data.weight
#         elseif edge_data.src_strand == Forward && edge_data.dst_strand == Reverse
#             transitions[src_forward, dst_reverse] = edge_data.weight * (1.0 - config.strand_switch_penalty)
#         elseif edge_data.src_strand == Reverse && edge_data.dst_strand == Forward
#             transitions[src_reverse, dst_forward] = edge_data.weight * (1.0 - config.strand_switch_penalty)
#         elseif edge_data.src_strand == Reverse && edge_data.dst_strand == Reverse
#             transitions[src_reverse, dst_reverse] = edge_data.weight
#         end
#     end
    
#     # Normalize transition probabilities
#     for i in 1:n_states
#         row_sum = sum(transitions[i, :])
#         if row_sum > 0
#             transitions[i, :] ./= row_sum
#         end
#     end
    
#     # Create emission probabilities (will be updated based on observations)
#     emissions = fill(config.match_prob, n_states)
    
#     return states, transitions, emissions
# end

# """
#     estimate_transition_probabilities(graph::MetaGraph, sequences::Vector) -> Matrix{Float64}

# Estimate transition probabilities from observed sequences in the graph.
# """
# function estimate_transition_probabilities(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
#                                           sequences::Vector)
#     vertices = collect(MetaGraphsNext.labels(graph))
#     n_vertices = length(vertices)
    
#     # Count transitions observed in sequences
#     transition_counts = zeros(Int, n_vertices, n_vertices)
#     vertex_to_idx = Dict(v => i for (i, v) in enumerate(vertices))
    
#     for seq_record in sequences
#         sequence = FASTX.sequence(String, seq_record)
#         k = length(first(vertices))  # Assume all k-mers same length
        
#         if length(sequence) >= k + 1
#             for i in 1:(length(sequence) - k)
#                 kmer1 = sequence[i:(i + k - 1)]
#                 kmer2 = sequence[(i + 1):(i + k)]
                
#                 if haskey(vertex_to_idx, kmer1) && haskey(vertex_to_idx, kmer2)
#                     idx1 = vertex_to_idx[kmer1]
#                     idx2 = vertex_to_idx[kmer2]
#                     transition_counts[idx1, idx2] += 1
#                 end
#             end
#         end
#     end
    
#     # Convert counts to probabilities
#     transition_probs = zeros(Float64, n_vertices, n_vertices)
#     for i in 1:n_vertices
#         row_sum = sum(transition_counts[i, :])
#         if row_sum > 0
#             transition_probs[i, :] = transition_counts[i, :] ./ row_sum
#         end
#     end
    
#     return transition_probs
# end

# """
#     viterbi_decode_next(graph::MetaGraph, observations::Vector, config::ViterbiConfig) -> ViterbiPath

# Enhanced Viterbi decoding with strand awareness and memory efficiency.
# """
# function viterbi_decode_next(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, String, KmerVertexData, KmerEdgeData},
#                             observations::Vector{String},
#                             config::ViterbiConfig = ViterbiConfig())
#     states, transitions, emissions = create_hmm_from_graph(graph, config)
#     n_states = length(states)
#     n_obs = length(observations)
    
#     if n_states == 0 || n_obs == 0
#         return ViterbiPath(ViterbiState[], -Inf, "", Tuple{Int, String, String}[])
#     end
    
#     # Use log probabilities for numerical stability
#     log_transitions = config.use_log_space ? log.(transitions .+ 1e-10) : transitions
#     log_emissions = config.use_log_space ? log.(emissions .+ 1e-10) : emissions
    
#     # Viterbi matrices
#     if config.use_log_space
#         viterbi_probs = fill(-Inf, n_states, n_obs)
#         viterbi_path = zeros(Int, n_states, n_obs)
#     else
#         viterbi_probs = zeros(Float64, n_states, n_obs)
#         viterbi_path = zeros(Int, n_states, n_obs)
#     end
    
#     # Initialize first column
#     for s in 1:n_states
#         emission_prob = calculate_emission_probability(states[s], observations[1], config)
#         if config.use_log_space
#             viterbi_probs[s, 1] = log(1.0 / n_states) + log(emission_prob + 1e-10)
#         else
#             viterbi_probs[s, 1] = (1.0 / n_states) * emission_prob
#         end
#     end
    
#     # Forward pass
#     for t in 2:n_obs
#         for s in 1:n_states
#             emission_prob = calculate_emission_probability(states[s], observations[t], config)
            
#             best_prob = config.use_log_space ? -Inf : 0.0
#             best_prev = 1
            
#             for prev_s in 1:n_states
#                 if config.use_log_space
#                     prob = viterbi_probs[prev_s, t-1] + log_transitions[prev_s, s] + log(emission_prob + 1e-10)
#                 else
#                     prob = viterbi_probs[prev_s, t-1] * transitions[prev_s, s] * emission_prob
#                 end
                
#                 if config.use_log_space ? (prob > best_prob) : (prob > best_prob)
#                     best_prob = prob
#                     best_prev = prev_s
#                 end
#             end
            
#             viterbi_probs[s, t] = best_prob
#             viterbi_path[s, t] = best_prev
#         end
#     end
    
#     # Backward pass to find best path
#     path_states = Vector{ViterbiState}()
    
#     # Find best final state
#     if config.use_log_space
#         best_final_prob, best_final_state = findmax(viterbi_probs[:, end])
#     else
#         best_final_prob, best_final_state = findmax(viterbi_probs[:, end])
#     end
    
#     # Trace back
#     current_state = best_final_state
#     for t in n_obs:-1:1
#         state_copy = ViterbiState(
#             states[current_state].vertex_label,
#             states[current_state].strand,
#             states[current_state].emission_prob,
#             t
#         )
#         pushfirst!(path_states, state_copy)
        
#         if t > 1
#             current_state = viterbi_path[current_state, t]
#         end
#     end
    
#     # Generate polished sequence and corrections
#     polished_sequence, corrections = generate_polished_sequence(path_states, observations, config)
    
#     return ViterbiPath(path_states, best_final_prob, polished_sequence, corrections)
# end

# """
#     calculate_emission_probability(state::ViterbiState, observation::String, config::ViterbiConfig) -> Float64

# Calculate emission probability for a state given an observation.
# """
# function calculate_emission_probability(state::ViterbiState, observation::String, config::ViterbiConfig)
#     kmer_seq = state.strand == Forward ? state.vertex_label : reverse_complement(state.vertex_label)
    
#     if kmer_seq == observation
#         return config.match_prob
#     else
#         # Calculate edit distance for mismatch probability
#         edit_dist = simple_edit_distance(kmer_seq, observation)
#         if edit_dist == 1
#             return config.mismatch_prob
#         else
#             return config.mismatch_prob / (edit_dist + 1)
#         end
#     end
# end

# """
#     simple_edit_distance(s1::String, s2::String) -> Int

# Simple edit distance calculation for k-mers.
# """
# function simple_edit_distance(s1::String, s2::String)
#     if length(s1) != length(s2)
#         return abs(length(s1) - length(s2))
#     end
    
#     mismatches = 0
#     for (c1, c2) in zip(s1, s2)
#         if c1 != c2
#             mismatches += 1
#         end
#     end
#     return mismatches
# end

# """
#     generate_polished_sequence(states::Vector{ViterbiState}, observations::Vector{String}, 
#                               config::ViterbiConfig) -> (String, Vector{Tuple{Int, String, String}})

# Generate polished sequence and track corrections made.
# """
# function generate_polished_sequence(states::Vector{ViterbiState}, 
#                                    observations::Vector{String},
#                                    config::ViterbiConfig)
#     if isempty(states)
#         return "", Tuple{Int, String, String}[]
#     end
    
#     polished_parts = String[]
#     corrections = Tuple{Int, String, String}[]
    
#     for (i, state) in enumerate(states)
#         kmer_seq = state.strand == Forward ? state.vertex_label : reverse_complement(state.vertex_label)
#         observed = observations[i]
        
#         if kmer_seq != observed
#             push!(corrections, (i, observed, kmer_seq))
#         end
        
#         # For overlapping k-mers, only add the new character
#         if i == 1
#             push!(polished_parts, kmer_seq)
#         else
#             # Add only the last character of the k-mer
#             push!(polished_parts, string(kmer_seq[end]))
#         end
#     end
    
#     polished_sequence = join(polished_parts)
#     return polished_sequence, corrections
# end

# """
#     viterbi_batch_process(graph::MetaGraph, sequences::Vector, config::ViterbiConfig) -> Vector{ViterbiPath}

# Process multiple sequences in batches for memory efficiency.
# """
# function viterbi_batch_process(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
#                               sequences::Vector,
#                               config::ViterbiConfig = ViterbiConfig())
#     results = Vector{ViterbiPath}()
#     n_sequences = length(sequences)
    
#     for batch_start in 1:config.batch_size:n_sequences
#         batch_end = min(batch_start + config.batch_size - 1, n_sequences)
#         batch = sequences[batch_start:batch_end]
        
#         println("Processing batch $(div(batch_start - 1, config.batch_size) + 1)/$(ceil(Int, n_sequences / config.batch_size))")
        
#         for seq_record in batch
#             sequence = FASTX.sequence(String, seq_record)
#             k = length(first(MetaGraphsNext.labels(graph)))
            
#             # Convert sequence to k-mer observations
#             observations = String[]
#             if length(sequence) >= k
#                 for i in 1:(length(sequence) - k + 1)
#                     push!(observations, sequence[i:(i + k - 1)])
#                 end
#             end
            
#             if !isempty(observations)
#                 viterbi_result = viterbi_decode_next(graph, observations, config)
#                 push!(results, viterbi_result)
#             end
#         end
        
#         # Force garbage collection between batches
#         GC.gc()
#     end
    
#     return results
# end

# """
#     polish_sequence_next(graph::MetaGraph, sequence::String, config::ViterbiConfig) -> ViterbiPath

# Polish a single sequence using Viterbi algorithm on k-mer graph.
# """
# function polish_sequence_next(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
#                              sequence::String,
#                              config::ViterbiConfig = ViterbiConfig())
#     k = length(first(MetaGraphsNext.labels(graph)))
    
#     # Convert sequence to k-mer observations
#     observations = String[]
#     if length(sequence) >= k
#         for i in 1:(length(sequence) - k + 1)
#             push!(observations, sequence[i:(i + k - 1)])
#         end
#     end
    
#     if isempty(observations)
#         return ViterbiPath(ViterbiState[], -Inf, sequence, Tuple{Int, String, String}[])
#     end
    
#     return viterbi_decode_next(graph, observations, config)
# end

# """
#     correct_errors_next(graph::MetaGraph, sequences::Vector, config::ViterbiConfig) -> Vector{FASTX.FASTA.Record}

# Correct errors in sequences using Viterbi algorithm and return corrected FASTA records.
# """
# function correct_errors_next(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
#                             sequences::Vector,
#                             config::ViterbiConfig = ViterbiConfig())
#     corrected_records = Vector{FASTX.FASTA.Record}()
    
#     for (i, seq_record) in enumerate(sequences)
#         original_sequence = FASTX.sequence(String, seq_record)
#         polished_result = polish_sequence_next(graph, original_sequence, config)
        
#         # Create new record with polished sequence
#         corrected_id = "$(FASTX.identifier(seq_record))_polished"
#         corrected_desc = "$(FASTX.description(seq_record)) | $(length(polished_result.corrections_made)) corrections"
        
#         corrected_record = FASTX.FASTA.Record(corrected_id, polished_result.polished_sequence)
#         push!(corrected_records, corrected_record)
        
#         # Log corrections if any were made
#         if !isempty(polished_result.corrections_made)
#             println("Sequence $i: $(length(polished_result.corrections_made)) corrections made")
#         end
#     end
    
#     return corrected_records
# end