struct Qualmer{KmerT, K}
    kmer::KmerT
    qualities::NTuple{K,UInt8}
    function Qualmer(kmer::KmerT, qualities::NTuple{K,UInt8}) where {KmerT<:Kmers.Kmer, K}
        @assert K == length(kmer) "Quality length must match kmer length"
        return new{KmerT,K}(kmer, qualities)
    end
end

function Qualmer(kmer::KmerT, qualities::AbstractVector{<:Integer}) where {KmerT<:Kmers.Kmer}
    @assert length(kmer) == length(qualities) "Quality length must match kmer length"
    k = length(kmer)
    Qualmer(kmer, NTuple{k, UInt8}(UInt8.(qualities)))
end

# Basic methods
Base.length(q::Qualmer) = length(q.kmer)
Base.getindex(q::Qualmer, i::Integer) = (q.kmer[i], q.qualities[i])
Base.:(==)(a::T, b::T) where {T<:Qualmer} = a.kmer == b.kmer && a.qualities == b.qualities
Base.hash(q::Qualmer, h::UInt) = hash(q.qualities, hash(q.kmer, h))

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the canonical representation of a DNA qualmer, considering both sequence and quality.
For DNA/RNA, this involves potentially reverse-complementing the k-mer and reversing
the quality scores accordingly.
"""
function canonical(qmer::Qualmer{KmerT}) where {KmerT<:Union{Kmers.DNAKmer, Kmers.RNAKmer}}
    canon_kmer = BioSequences.canonical(qmer.kmer)
    
    if canon_kmer == BioSequences.reverse_complement(qmer.kmer)
        K = length(qmer.qualities)
        reversed_quals = ntuple(i -> qmer.qualities[K - i + 1], K)
        return Qualmer(canon_kmer, reversed_quals)
    else
        return Qualmer(canon_kmer, qmer.qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function canonical(qmer::Qualmer{<:Kmers.AAKmer})
    return qmer
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create an iterator that yields DNA qualmers from the given sequence and quality scores.
"""
function qualmers_fw(sequence::BioSequences.LongSequence{A}, quality::AbstractVector{<:Integer}, ::Val{K}) where {A,K}
    return (
        (Qualmer(
            kmer, 
            quality[pos:pos+K-1]
        ), pos) for (pos, kmer) in enumerate(Kmers.FwKmers{A, K}(sequence))
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_unambiguous(sequence::BioSequences.LongSequence{BioSequences.DNAAlphabet{N}}, quality::AbstractVector{<:Integer}, ::Val{K}) where {N,K}
    return (
        (Qualmer(
            kmer,
            quality[pos:pos+K-1]
        ), pos) for (kmer, pos) in Kmers.UnambiguousDNAMers{K}(sequence)
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_unambiguous(sequence::BioSequences.LongSequence{BioSequences.RNAAlphabet{N}}, quality::AbstractVector{<:Integer}, ::Val{K}) where {N,K}
    return (
        (Qualmer(
            kmer,
            NTuple{K,UInt8}(quality[pos:pos+K-1])
        ), pos) for (kmer, pos) in Kmers.UnambiguousRNAMers{K}(sequence)
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_unambiguous(sequence::BioSequences.LongSequence{BioSequences.AminoAcidAlphabet}, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    return qualmers_fw(sequence, quality, Val(K))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_canonical(sequence::BioSequences.LongSequence{BioSequences.DNAAlphabet{N}}, quality::AbstractVector{<:Integer}, ::Val{K}) where {N,K}
    return (
        (Qualmer(
            kmer,
            quality[pos:pos+K-1]
        ), pos) for (pos, kmer) in enumerate(Kmers.CanonicalKmers{BioSequences.DNAAlphabet{N}, K}(sequence))
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_canonical(sequence::BioSequences.LongSequence{BioSequences.AminoAcidAlphabet}, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    return qualmers_fw(sequence, quality, Val(K))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_fw(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_fw(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_canonical(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.FASTQ.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_canonical(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function qualmers_unambiguous(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_unambiguous(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate unambiguous canonical qualmers from the given sequence and quality scores.
This function first extracts all unambiguous qualmers and then applies canonical 
representation to each one.
"""
# Also update functions that consume these results, like qualmers_unambiguous_canonical
function qualmers_unambiguous_canonical(sequence, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    # Return both canonical qualmer and position
    return ((canonical(qmer), pos) for (qmer, pos) in qualmers_unambiguous(sequence, quality, Val(K)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate unambiguous canonical qualmers from the given FASTQ record.
"""
function qualmers_unambiguous_canonical(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_unambiguous_canonical(sequence, quality, Val(k))
end

# ============================================================================
# Joint Probability Calculations for Quality-Aware Assembly
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert PHRED quality score to probability of correctness.
PHRED score Q relates to error probability P by: Q = -10 * log10(P)
Therefore, correctness probability = 1 - P = 1 - 10^(-Q/10)
"""
function phred_to_probability(phred_score::UInt8)
    return 1.0 - 10.0^(-phred_score / 10.0)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert PHRED quality score to error probability.
"""
function phred_to_error_probability(phred_score::UInt8)
    return 10.0^(-phred_score / 10.0)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate joint probability that a single qualmer is correct.
For a k-mer with quality scores [q1, q2, ..., qk], the joint probability
that all positions are correct is: ∏(1 - 10^(-qi/10))
"""
function qualmer_correctness_probability(qmer::Qualmer)
    log_prob = 0.0
    for q in qmer.qualities
        prob_correct = phred_to_probability(q)
        log_prob += log(prob_correct)
    end
    return exp(log_prob)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate joint probability that multiple observations of the same k-mer are all correct.
For independent observations of the same k-mer sequence, this represents our
confidence that the k-mer truly exists in the data.

# Arguments
- `qualmers`: Vector of Qualmer observations of the same k-mer sequence
- `use_log_space`: Use log-space arithmetic for numerical stability (default: true)

# Returns
- Float64: Joint probability that all observations are correct
"""
function joint_qualmer_probability(qualmers::Vector{<:Qualmer}; use_log_space::Bool=true)
    if isempty(qualmers)
        return 0.0
    end
    
    # Verify all qualmers have the same k-mer sequence
    kmer_seq = qualmers[1].kmer
    for qmer in qualmers[2:end]
        if qmer.kmer != kmer_seq
            throw(ArgumentError("All qualmers must have the same k-mer sequence"))
        end
    end
    
    if use_log_space
        log_joint_prob = 0.0
        for qmer in qualmers
            prob_correct = qualmer_correctness_probability(qmer)
            log_joint_prob += log(prob_correct)
        end
        return exp(log_joint_prob)
    else
        joint_prob = 1.0
        for qmer in qualmers
            joint_prob *= qualmer_correctness_probability(qmer)
        end
        return joint_prob
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate position-wise joint probability for multiple qualmer observations.
This is more sophisticated than `joint_qualmer_probability` as it considers
the quality at each position across all observations.

# Arguments
- `qualmers`: Vector of Qualmer observations of the same k-mer sequence
- `use_log_space`: Use log-space arithmetic for numerical stability (default: true)

# Returns
- Float64: Position-wise joint probability
"""
function position_wise_joint_probability(qualmers::Vector{<:Qualmer}; use_log_space::Bool=true)
    if isempty(qualmers)
        return 0.0
    end
    
    k = length(qualmers[1])
    
    # Verify all qualmers have the same length and k-mer sequence
    kmer_seq = qualmers[1].kmer
    for qmer in qualmers[2:end]
        if length(qmer) != k || qmer.kmer != kmer_seq
            throw(ArgumentError("All qualmers must have the same k-mer sequence and length"))
        end
    end
    
    if use_log_space
        log_joint_prob = 0.0
        for pos in 1:k
            # Calculate joint probability for this position across all observations
            log_pos_prob = 0.0
            for qmer in qualmers
                prob_correct = phred_to_probability(qmer.qualities[pos])
                log_pos_prob += log(prob_correct)
            end
            log_joint_prob += log_pos_prob
        end
        return exp(log_joint_prob)
    else
        joint_prob = 1.0
        for pos in 1:k
            pos_prob = 1.0
            for qmer in qualmers
                pos_prob *= phred_to_probability(qmer.qualities[pos])
            end
            joint_prob *= pos_prob
        end
        return joint_prob
    end
end

# ============================================================================
# Graph Types for Quality-Aware Assembly
# ============================================================================

"""
Qualmer observation for tracking k-mer occurrences in sequences.
"""
struct QualmerObservation{QualmerT<:Qualmer}
    qualmer::QualmerT
    sequence_id::Int
    position::Int
    
    function QualmerObservation(qualmer::QualmerT, sequence_id::Int, position::Int) where {QualmerT<:Qualmer}
        @assert sequence_id >= 1 "Sequence ID must be positive"
        @assert position >= 1 "Position must be positive"
        new{QualmerT}(qualmer, sequence_id, position)
    end
end

"""
Vertex data for quality-aware k-mer graphs.
"""
struct QualmerVertexData{QualmerT<:Qualmer}
    canonical_qualmer::QualmerT              # Canonical representation
    observations::Vector{QualmerObservation{QualmerT}}  # All observations
    joint_probability::Float64               # Joint probability of correctness
    coverage::Int                           # Number of observations
    mean_quality::Float64                   # Mean quality across all observations
    
    function QualmerVertexData(observations::Vector{QualmerObservation{QualmerT}}) where {QualmerT<:Qualmer}
        @assert !isempty(observations) "Must have at least one observation"
        
        # Get canonical qualmer (all observations should have the same canonical sequence)
        canonical_qualmer = canonical(observations[1].qualmer)
        
        # Calculate joint probability
        qualmers = [obs.qualmer for obs in observations]
        joint_prob = position_wise_joint_probability(qualmers)
        
        coverage = length(observations)
        
        # Calculate mean quality across all observations
        all_qualities = UInt8[]
        for obs in observations
            append!(all_qualities, obs.qualmer.qualities)
        end
        mean_quality = isempty(all_qualities) ? 0.0 : Statistics.mean(all_qualities)
        
        new{QualmerT}(canonical_qualmer, observations, joint_prob, coverage, mean_quality)
    end
end

"""
Edge data for quality-aware k-mer graphs.
"""
struct QualmerEdgeData
    observations::Vector{Tuple{Int, Int}}   # (sequence_id, position) pairs
    src_strand::StrandOrientation          # Source vertex strand
    dst_strand::StrandOrientation          # Destination vertex strand
    weight::Float64                        # Edge weight based on quality and coverage
    quality_weight::Float64                # Additional weight from quality scores
    
    function QualmerEdgeData(observations::Vector{Tuple{Int, Int}},
                           src_strand::StrandOrientation, dst_strand::StrandOrientation,
                           src_vertex_data::QualmerVertexData, dst_vertex_data::QualmerVertexData)
        @assert !isempty(observations) "Edge must have at least one observation"
        
        # Base weight from coverage (number of observations)
        base_weight = Float64(length(observations))
        
        # Quality weight from joint probabilities of connected k-mers
        quality_weight = src_vertex_data.joint_probability * dst_vertex_data.joint_probability
        
        # Combined weight
        total_weight = base_weight * quality_weight
        
        new(observations, src_strand, dst_strand, total_weight, quality_weight)
    end
end

# ============================================================================
# Quality-Aware Graph Construction
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a quality-aware k-mer graph from FASTQ records using existing Qualmer functionality.
This function leverages the existing qualmer extraction functions and adds graph construction.

# Arguments
- `fastq_records`: Vector of FASTQ records with quality scores
- `k`: K-mer size
- `graph_mode`: SingleStrand or DoubleStrand mode (default: DoubleStrand)
- `min_quality`: Minimum average PHRED quality to include k-mer (default: 10)
- `min_coverage`: Minimum coverage (observations) to include k-mer (default: 1)

# Returns
- MetaGraphsNext.MetaGraph with QualmerVertexData and QualmerEdgeData
"""
function build_qualmer_graph(fastq_records::Vector{FASTX.FASTQ.Record}; 
                            k::Int=31, graph_mode::GraphMode=DoubleStrand,
                            min_quality::UInt8=UInt8(10), min_coverage::Int=1)
    
    # Create the MetaGraph
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=QualmerVertexData,
        edge_data_type=QualmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Track qualmer observations by canonical k-mer sequence
    canonical_observations = Dict{String, Vector{QualmerObservation}}()
    
    # Process each FASTQ record
    for (seq_id, record) in enumerate(fastq_records)
        # Extract qualmers based on graph mode
        qualmer_iterator = if graph_mode == DoubleStrand
            qualmers_unambiguous_canonical(record, k)
        else
            qualmers_unambiguous(record, k)
        end
        
        # Process each qualmer
        for (qmer, pos) in qualmer_iterator
            # Check quality threshold
            mean_qual = Statistics.mean(qmer.qualities)
            if mean_qual >= min_quality
                # Get canonical representation
                canonical_qmer = graph_mode == DoubleStrand ? qmer : canonical(qmer)
                canonical_seq = string(canonical_qmer.kmer)
                
                # Create observation
                observation = QualmerObservation(qmer, seq_id, pos)
                
                # Store observation
                if !haskey(canonical_observations, canonical_seq)
                    canonical_observations[canonical_seq] = QualmerObservation[]
                end
                push!(canonical_observations[canonical_seq], observation)
            end
        end
    end
    
    # Filter by minimum coverage and create vertices
    vertex_data_map = Dict{String, QualmerVertexData}()
    for (canonical_seq, observations) in canonical_observations
        if length(observations) >= min_coverage
            vertex_data = QualmerVertexData(observations)
            vertex_data_map[canonical_seq] = vertex_data
            graph[canonical_seq] = vertex_data
        end
    end
    
    # Create edges based on k-mer adjacency
    _add_qualmer_edges!(graph, vertex_data_map, fastq_records, k, graph_mode)
    
    return graph
end

"""
Add edges to the qualmer graph based on k-mer adjacency in sequences.
"""
function _add_qualmer_edges!(graph::MetaGraphsNext.MetaGraph, 
                            vertex_data_map::Dict{String, QualmerVertexData},
                            fastq_records::Vector{FASTX.FASTQ.Record}, 
                            k::Int, graph_mode::GraphMode)
    
    # Track edge observations
    edge_observations = Dict{Tuple{String, String}, Vector{Tuple{Int, Int}}}()
    edge_strand_info = Dict{Tuple{String, String}, Tuple{StrandOrientation, StrandOrientation}}()
    
    for (seq_id, record) in enumerate(fastq_records)
        # Extract consecutive qualmers
        qualmer_iterator = if graph_mode == DoubleStrand
            qualmers_unambiguous_canonical(record, k)
        else
            qualmers_unambiguous(record, k)
        end
        
        qualmer_pairs = collect(qualmer_iterator)
        
        # Process consecutive pairs
        for i in 1:(length(qualmer_pairs) - 1)
            qmer1, pos1 = qualmer_pairs[i]
            qmer2, pos2 = qualmer_pairs[i + 1]
            
            # Check if positions are adjacent (allowing for gaps in unambiguous extraction)
            if pos2 == pos1 + 1
                # Get canonical representations
                canonical_qmer1 = graph_mode == DoubleStrand ? qmer1 : canonical(qmer1)
                canonical_qmer2 = graph_mode == DoubleStrand ? qmer2 : canonical(qmer2)
                
                canonical_seq1 = string(canonical_qmer1.kmer)
                canonical_seq2 = string(canonical_qmer2.kmer)
                
                # Check if both vertices exist in graph
                if haskey(vertex_data_map, canonical_seq1) && haskey(vertex_data_map, canonical_seq2)
                    # Determine strand orientations
                    src_strand = _determine_strand(qmer1, canonical_qmer1)
                    dst_strand = _determine_strand(qmer2, canonical_qmer2)
                    
                    # Validate biological transition
                    if _is_valid_qualmer_transition(canonical_qmer1.kmer, canonical_qmer2.kmer, src_strand, dst_strand)
                        # Track edge observation
                        edge_key = (canonical_seq1, canonical_seq2)
                        if !haskey(edge_observations, edge_key)
                            edge_observations[edge_key] = Tuple{Int, Int}[]
                            edge_strand_info[edge_key] = (src_strand, dst_strand)
                        end
                        push!(edge_observations[edge_key], (seq_id, pos1))
                    end
                end
            end
        end
    end
    
    # Add edges to graph
    for ((src_seq, dst_seq), observations) in edge_observations
        src_strand, dst_strand = edge_strand_info[(src_seq, dst_seq)]
        src_vertex_data = vertex_data_map[src_seq]
        dst_vertex_data = vertex_data_map[dst_seq]
        
        # Create edge data
        edge_data = QualmerEdgeData(observations, src_strand, dst_strand, src_vertex_data, dst_vertex_data)
        
        # Add edge to graph
        graph[src_seq, dst_seq] = edge_data
    end
end

"""
Determine strand orientation by comparing original and canonical qualmers.
"""
function _determine_strand(original_qmer::Qualmer, canonical_qmer::Qualmer)
    return original_qmer.kmer == canonical_qmer.kmer ? Forward : Reverse
end

"""
Validate that the transition between two k-mers is biologically valid.
"""
function _is_valid_qualmer_transition(kmer1, kmer2, src_strand::StrandOrientation, dst_strand::StrandOrientation)
    # Get string representations based on strand
    str1 = src_strand == Forward ? string(kmer1) : string(BioSequences.reverse_complement(kmer1))
    str2 = dst_strand == Forward ? string(kmer2) : string(BioSequences.reverse_complement(kmer2))
    
    # Check if k-mers overlap properly (k-1 overlap)
    return str1[2:end] == str2[1:end-1]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get comprehensive statistics about a qualmer graph.
"""
function get_qualmer_statistics(graph::MetaGraphsNext.MetaGraph)
    vertices = collect(MetaGraphsNext.labels(graph))
    
    if isempty(vertices)
        return Dict{String, Any}(
            "num_vertices" => 0,
            "num_edges" => 0,
            "mean_coverage" => 0.0,
            "mean_quality" => 0.0,
            "mean_joint_probability" => 0.0
        )
    end
    
    # Collect vertex statistics
    coverages = Float64[]
    qualities = Float64[]
    joint_probs = Float64[]
    
    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        push!(coverages, vertex_data.coverage)
        push!(qualities, vertex_data.mean_quality)
        push!(joint_probs, vertex_data.joint_probability)
    end
    
    return Dict{String, Any}(
        "num_vertices" => length(vertices),
        "num_edges" => length(collect(MetaGraphsNext.edge_labels(graph))),
        "mean_coverage" => Statistics.mean(coverages),
        "mean_quality" => Statistics.mean(qualities),
        "mean_joint_probability" => Statistics.mean(joint_probs),
        "median_coverage" => Statistics.median(coverages),
        "median_quality" => Statistics.median(qualities),
        "median_joint_probability" => Statistics.median(joint_probs),
        "total_observations" => sum(Int(c) for c in coverages)
    )
end