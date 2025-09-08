"""
$(DocStringExtensions.TYPEDEF)

A quality-aware k-mer that pairs sequence k-mers with their corresponding quality scores.

# Fields
$(DocStringExtensions.TYPEDFIELDS)

# Examples
```julia
using Mycelia
using Kmers

# Create a DNA k-mer with quality scores
kmer = DNAKmer("ATCG")
qualities = [30, 35, 32, 28]  # Phred quality scores
qmer = Qualmer(kmer, qualities)

# Access individual bases and qualities
base, quality = qmer[1]  # Returns ('A', 30)
```

This structure is essential for quality-aware sequence analysis and assembly algorithms
that consider both sequence content and sequencing quality information.
"""
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
    
    # If the canonical form is different from the original, it means we reverse-complemented
    # In that case, we need to reverse the quality scores
    if canon_kmer != qmer.kmer
        K = length(qmer.qualities)
        reversed_quals = ntuple(i -> qmer.qualities[K - i + 1], K)
        return Qualmer(canon_kmer, reversed_quals)
    else
        return Qualmer(canon_kmer, qmer.qualities)
    end
end

# Additional method for the Mer types from FwKmers
function canonical(qmer::Qualmer{KmerT}) where {K, KmerT<:Kmers.Mer{K, BioSequences.DNAAlphabet{N}}} where N
    canon_kmer = BioSequences.canonical(qmer.kmer)
    
    # If the canonical form is different from the original, it means we reverse-complemented
    # In that case, we need to reverse the quality scores
    if canon_kmer != qmer.kmer
        reversed_quals = ntuple(i -> qmer.qualities[K - i + 1], K)
        return Qualmer(canon_kmer, reversed_quals)
    else
        return Qualmer(canon_kmer, qmer.qualities)
    end
end

# Additional method for RNA Mer types
function canonical(qmer::Qualmer{KmerT}) where {K, KmerT<:Kmers.Mer{K, BioSequences.RNAAlphabet{N}}} where N
    canon_kmer = BioSequences.canonical(qmer.kmer)
    
    # If the canonical form is different from the original, it means we reverse-complemented
    # In that case, we need to reverse the quality scores
    if canon_kmer != qmer.kmer
        reversed_quals = ntuple(i -> qmer.qualities[K - i + 1], K)
        return Qualmer(canon_kmer, reversed_quals)
    else
        return Qualmer(canon_kmer, qmer.qualities)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the canonical representation of an amino acid qualmer.
For amino acid sequences, the qualmer is returned unchanged since
amino acids do not have reverse complements like DNA/RNA.
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

Generate unambiguous DNA qualmers from sequence and quality scores.
Only includes k-mers that contain unambiguous nucleotides (A, T, G, C).
Skips k-mers containing ambiguous bases like N.

# Arguments
- `sequence`: DNA sequence as LongDNA
- `quality`: Vector of quality scores corresponding to each base
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples for each unambiguous k-mer
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

Generate unambiguous RNA qualmers from sequence and quality scores.
Only includes k-mers that contain unambiguous nucleotides (A, U, G, C).
Skips k-mers containing ambiguous bases like N.

# Arguments
- `sequence`: RNA sequence as LongRNA
- `quality`: Vector of quality scores corresponding to each base
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples for each unambiguous k-mer
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

Generate amino acid qualmers from protein sequence and quality scores.
For amino acid sequences, all k-mers are considered unambiguous, so this
falls back to the forward qualmer generation.

# Arguments
- `sequence`: Protein sequence as LongAA
- `quality`: Vector of quality scores corresponding to each amino acid
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples for each k-mer
"""
function qualmers_unambiguous(sequence::BioSequences.LongSequence{BioSequences.AminoAcidAlphabet}, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    return qualmers_fw(sequence, quality, Val(K))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate canonical DNA qualmers from sequence and quality scores.
For each k-mer, returns the lexicographically smaller of the k-mer and its
reverse complement, ensuring consistent representation regardless of strand.

# Arguments
- `sequence`: DNA sequence as LongDNA
- `quality`: Vector of quality scores corresponding to each base
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples with canonical k-mer representation
"""
function qualmers_canonical(sequence::BioSequences.LongSequence{BioSequences.DNAAlphabet{N}}, quality::AbstractVector{<:Integer}, ::Val{K}) where {N,K}
    return (
        (canonical(Qualmer(
            kmer,
            quality[pos:pos+K-1]
        )), pos) for (pos, kmer) in enumerate(Kmers.FwKmers{BioSequences.DNAAlphabet{N}, K}(sequence))
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate canonical amino acid qualmers from protein sequence and quality scores.
For amino acid sequences, canonical representation is the same as forward
representation since proteins don't have reverse complements.

# Arguments
- `sequence`: Protein sequence as LongAA
- `quality`: Vector of quality scores corresponding to each amino acid
- `::Val{K}`: K-mer size as a type parameter

# Returns
- Iterator yielding (Qualmer, position) tuples for each k-mer
"""
function qualmers_canonical(sequence::BioSequences.LongSequence{BioSequences.AminoAcidAlphabet}, quality::AbstractVector{<:Integer}, ::Val{K}) where {K}
    return qualmers_fw(sequence, quality, Val(K))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate forward qualmers from a FASTQ record.
Extracts all k-mers in the forward direction with their associated
quality scores from the FASTQ record.

# Arguments
- `record`: FASTQ record containing sequence and quality data
- `k`: K-mer size

# Returns
- Iterator yielding (Qualmer, position) tuples for each k-mer
"""
function qualmers_fw(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_fw(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate canonical qualmers from a FASTQ record.
Extracts k-mers in canonical representation (lexicographically smaller
of forward and reverse complement) with their associated quality scores.

# Arguments
- `record`: FASTQ record containing sequence and quality data
- `k`: K-mer size

# Returns
- Iterator yielding (Qualmer, position) tuples with canonical k-mer representation
"""
function qualmers_canonical(record::FASTX.FASTQ.Record, k::Int)
    sequence = Mycelia.convert_sequence(FASTX.FASTQ.sequence(record))
    quality = collect(FASTX.quality_scores(record))
    return qualmers_canonical(sequence, quality, Val(k))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate unambiguous qualmers from a FASTQ record.
Extracts only k-mers that contain unambiguous nucleotides (no N bases)
with their associated quality scores.

# Arguments
- `record`: FASTQ record containing sequence and quality data
- `k`: K-mer size

# Returns
- Iterator yielding (Qualmer, position) tuples for each unambiguous k-mer
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
    # Convert PHRED score to probability that base call is correct
    # PHRED score Q = -10 * log10(P_error)
    # Therefore: P_error = 10^(-Q/10)
    # P_correct = 1 - P_error = 1 - 10^(-Q/10)
    error_prob = 10.0^(-Float64(phred_score) / 10.0)
    correct_prob = 1.0 - error_prob
    
    # Ensure we don't return negative probabilities
    return max(0.0, min(1.0, correct_prob))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert PHRED quality score to error probability.
"""
function phred_to_error_probability(phred_score::UInt8)
    return 10.0^(-Float64(phred_score) / 10.0)
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
                if prob_correct <= 0.0
                    @warn "Invalid probability $prob_correct for quality $(qmer.qualities[pos]). Using minimum probability."
                    prob_correct = 1e-10  # Use a very small but positive probability
                end
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
$(DocStringExtensions.TYPEDEF)

Vertex data for quality-aware k-mer graphs that stores k-mer observations with quality scores.

# Fields
$(DocStringExtensions.TYPEDFIELDS)

This structure aggregates multiple observations of the same k-mer and computes quality-based
statistics for use in quality-aware assembly algorithms.
"""
struct QualmerVertexData{QualmerT<:Qualmer}
    canonical_qualmer::QualmerT              # Canonical representation
    observations::Vector{QualmerObservation{QualmerT}}  # All observations
    joint_probability::Float64               # Joint probability of correctness
    coverage::Int                           # Number of observations
    mean_quality::Float64                   # Mean quality across all observations
    
    function QualmerVertexData(observations::Vector{<:QualmerObservation})
        @assert !isempty(observations) "Must have at least one observation"
        
        # Infer the QualmerT type from the first observation
        QualmerT = typeof(observations[1].qualmer)
        
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
$(DocStringExtensions.TYPEDEF)

Edge data for quality-aware k-mer graphs that connects k-mer vertices with strand information.

# Fields
$(DocStringExtensions.TYPEDFIELDS)

This structure represents connections between k-mers in quality-aware graphs, including
strand orientation and quality-weighted edge weights for assembly algorithms.
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
    
    # Determine k-mer type from first record
    first_seq = Mycelia.convert_sequence(FASTX.sequence(fastq_records[1]))
    kmer_type = if first_seq isa BioSequences.LongDNA
        Kmers.DNAKmer{k}
    elseif first_seq isa BioSequences.LongRNA
        Kmers.RNAKmer{k}
    else
        Kmers.AAKmer{k}
    end
    
    # Track qualmer observations by canonical k-mer
    # We'll create the graph after we know the actual k-mer type
    canonical_observations = Dict{Any, Vector{QualmerObservation}}()
    
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
                canonical_kmer = canonical_qmer.kmer
                
                # Create observation
                observation = QualmerObservation(qmer, seq_id, pos)
                
                # Store observation
                if !haskey(canonical_observations, canonical_kmer)
                    canonical_observations[canonical_kmer] = QualmerObservation[]
                end
                push!(canonical_observations[canonical_kmer], observation)
            end
        end
    end
    
    # Now create the graph with the actual k-mer type (inferred from data)
    if isempty(canonical_observations)
        @warn "No qualmer observations found"
        return MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=kmer_type,
            vertex_data_type=QualmerVertexData,
            edge_data_type=QualmerEdgeData,
            weight_function=edge_data -> edge_data.weight,
            default_weight=0.0
        )
    end
    
    # Infer actual k-mer type from the observations
    first_canonical_kmer = first(keys(canonical_observations))
    actual_kmer_type = typeof(first_canonical_kmer)
    
    # Create the MetaGraph with proper k-mer type labels
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=actual_kmer_type,
        vertex_data_type=QualmerVertexData,
        edge_data_type=QualmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Filter by minimum coverage and create vertices
    vertex_data_map = Dict{actual_kmer_type, QualmerVertexData}()
    for (canonical_kmer, observations) in canonical_observations
        if length(observations) >= min_coverage
            vertex_data = QualmerVertexData(observations)
            vertex_data_map[canonical_kmer] = vertex_data
            graph[canonical_kmer] = vertex_data
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
                            vertex_data_map::Dict{KmerT, QualmerVertexData},
                            fastq_records::Vector{FASTX.FASTQ.Record}, 
                            k::Int, graph_mode::GraphMode) where {KmerT}
    
    # Track edge observations
    edge_observations = Dict{Tuple{KmerT, KmerT}, Vector{Tuple{Int, Int}}}()
    edge_strand_info = Dict{Tuple{KmerT, KmerT}, Tuple{StrandOrientation, StrandOrientation}}()
    
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
                
                canonical_kmer1 = canonical_qmer1.kmer
                canonical_kmer2 = canonical_qmer2.kmer
                
                # Check if both vertices exist in graph
                if haskey(vertex_data_map, canonical_kmer1) && haskey(vertex_data_map, canonical_kmer2)
                    # Determine strand orientations
                    src_strand = _determine_strand(qmer1, canonical_qmer1)
                    dst_strand = _determine_strand(qmer2, canonical_qmer2)
                    
                    # Validate biological transition
                    if _is_valid_qualmer_transition(canonical_qmer1.kmer, canonical_qmer2.kmer, src_strand, dst_strand)
                        # Track edge observation
                        edge_key = (canonical_kmer1, canonical_kmer2)
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
    for ((src_kmer, dst_kmer), observations) in edge_observations
        src_strand, dst_strand = edge_strand_info[(src_kmer, dst_kmer)]
        src_vertex_data = vertex_data_map[src_kmer]
        dst_vertex_data = vertex_data_map[dst_kmer]
        
        # Create edge data
        edge_data = QualmerEdgeData(observations, src_strand, dst_strand, src_vertex_data, dst_vertex_data)
        
        # Add edge to graph
        graph[src_kmer, dst_kmer] = edge_data
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Find a quality-weighted path through a qualmer graph starting from a given vertex.
Uses joint probability as the primary weighting factor for path selection.

# Arguments
- `graph`: Qualmer graph (MetaGraphsNext with QualmerVertexData)
- `start_vertex`: Starting vertex for path traversal
- `max_path_length::Int=1000`: Maximum path length to prevent infinite loops

# Returns
- `Vector{Int}`: Path as sequence of vertex indices

# Details
At each step, selects the unvisited neighbor with the highest joint probability.
Terminates when no unvisited neighbors are available or max length is reached.
"""
function find_quality_weighted_path(graph, start_vertex; max_path_length::Int=1000)
    path = [start_vertex]
    current = start_vertex
    visited = Set([current])
    
    while length(path) < max_path_length
        # Get outgoing edges
        neighbors = Graphs.outneighbors(graph, current)
        
        # Filter out visited vertices
        unvisited = filter(n -> n ∉ visited, neighbors)
        
        if isempty(unvisited)
            break
        end
        
        # Choose next vertex based on joint probability
        next_vertex = unvisited[argmax([graph[v].joint_probability for v in unvisited])]
        
        push!(path, next_vertex)
        push!(visited, next_vertex)
        current = next_vertex
    end
    
    return path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate comprehensive assembly quality metrics for a qualmer graph.

# Arguments
- `graph`: Qualmer graph (MetaGraphsNext with QualmerVertexData)
- `low_confidence_threshold::Float64=0.95`: Threshold for identifying low-confidence k-mers

# Returns
- `NamedTuple`: Assembly quality metrics including coverage, quality, and confidence statistics

# Details
Calculates mean values for coverage, quality scores, and joint probabilities.
Identifies fraction of k-mers below confidence threshold as potential error indicators.
"""
function calculate_assembly_quality_metrics(qualmer_graph; low_confidence_threshold::Float64=0.95)
    # Get all vertex data
    vertices = [qualmer_graph[v] for v in MetaGraphsNext.labels(qualmer_graph)]
    
    if isempty(vertices)
        return (
            mean_coverage = 0.0,
            mean_quality = 0.0,
            mean_confidence = 0.0,
            low_confidence_fraction = 0.0,
            total_kmers = 0
        )
    end
    
    # Calculate metrics
    mean_coverage = Statistics.mean([v.coverage for v in vertices])
    mean_quality = Statistics.mean([v.mean_quality for v in vertices])
    mean_confidence = Statistics.mean([v.joint_probability for v in vertices])
    
    # Find low-confidence k-mers (potential errors)
    low_conf_kmers = filter(v -> v.joint_probability < low_confidence_threshold, vertices)
    
    return (
        mean_coverage = mean_coverage,
        mean_quality = mean_quality,
        mean_confidence = mean_confidence,
        low_confidence_fraction = length(low_conf_kmers) / length(vertices),
        total_kmers = length(vertices)
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identify potential sequencing errors based on quality scores and coverage patterns.

# Arguments
- `graph`: Qualmer graph (MetaGraphsNext with QualmerVertexData)
- `min_coverage::Int=2`: Minimum coverage for reliable k-mers
- `min_quality::Float64=20.0`: Minimum mean quality score
- `min_confidence::Float64=0.95`: Minimum joint probability threshold

# Returns
- `Vector`: K-mer labels of potential error k-mers

# Details
Identifies k-mers that are likely errors based on:
- Low coverage (singleton or few observations)
- Low quality scores
- Low joint probability (low confidence)
"""
function identify_potential_errors(graph; 
                                 min_coverage::Int=2, 
                                 min_quality::Float64=20.0, 
                                 min_confidence::Float64=0.95)
    error_labels = []
    
    for v in MetaGraphsNext.labels(graph)
        vdata = graph[v]
        
        # Check if vertex meets error criteria
        is_low_coverage = vdata.coverage < min_coverage
        is_low_quality = vdata.mean_quality < min_quality
        is_low_confidence = vdata.joint_probability < min_confidence
        
        # Consider it an error if it fails any criterion
        if is_low_coverage || is_low_quality || is_low_confidence
            push!(error_labels, v)
        end
    end
    
    return error_labels
end