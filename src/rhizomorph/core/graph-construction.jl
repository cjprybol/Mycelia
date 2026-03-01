# Graph Construction Functions
#
# Strand-specific graph construction algorithms for the Rhizomorph ecosystem.
#
# CRITICAL DESIGN PRINCIPLES:
# 1. Store k-mers AS OBSERVED (not canonical) - strand-specific by default
# 2. Track evidence with dataset_id and observation_id
# 3. Evidence stored in nested Dict structure for O(1) queries
# 4. Use MetaGraphsNext with proper vertex/edge data structures
# 5. Type-stable core functions with outer wrappers for type detection
#
# Based on rhizomorph-graph-ecosystem-plan.md section 1.4

# ============================================================================
# K-mer Graph Construction (Quality-unaware)
# ============================================================================

function _normalize_alphabet_symbol(value::Symbol)
    return Symbol(uppercase(String(value)))
end

function _normalize_ambiguous_action(value::Symbol)
    return Symbol(lowercase(String(value)))
end

function _possible_alphabets(seq_str::AbstractString)
    if isempty(seq_str)
        throw(ArgumentError("Input sequence cannot be empty."))
    end

    seq_chars = Set(seq_str)
    possible = Set{Symbol}()
    unambiguous = Set{Symbol}()

    if issubset(seq_chars, parentmodule(Rhizomorph).UNAMBIGUOUS_DNA_CHARSET)
        push!(possible, :DNA)
        push!(unambiguous, :DNA)
    end
    if issubset(seq_chars, parentmodule(Rhizomorph).UNAMBIGUOUS_RNA_CHARSET)
        push!(possible, :RNA)
        push!(unambiguous, :RNA)
    end
    if issubset(seq_chars, parentmodule(Rhizomorph).UNAMBIGUOUS_AA_CHARSET)
        push!(possible, :AA)
        push!(unambiguous, :AA)
    end

    if issubset(seq_chars, parentmodule(Rhizomorph).AMBIGUOUS_DNA_CHARSET)
        push!(possible, :DNA)
    end
    if issubset(seq_chars, parentmodule(Rhizomorph).AMBIGUOUS_RNA_CHARSET)
        push!(possible, :RNA)
    end
    if issubset(seq_chars, parentmodule(Rhizomorph).AMBIGUOUS_AA_CHARSET)
        push!(possible, :AA)
    end

    if isempty(possible)
        all_chars = union(
            parentmodule(Rhizomorph).AMBIGUOUS_DNA_CHARSET,
            parentmodule(Rhizomorph).AMBIGUOUS_RNA_CHARSET,
            parentmodule(Rhizomorph).AMBIGUOUS_AA_CHARSET
        )
        unmatched_chars = setdiff(seq_chars, all_chars)
        throw(ArgumentError("Sequence contains characters that do not belong to any known alphabet: $(join(unmatched_chars, ", "))"))
    end

    return (unambiguous = unambiguous, possible = possible)
end

function _resolve_ambiguous_alphabet(possible_alphabets::Set{Symbol}, ambiguous_action::Symbol)
    if isempty(possible_alphabets)
        error("Mixed sequence types detected. Cannot build graph from heterogeneous input.")
    end

    if length(possible_alphabets) == 1
        return first(possible_alphabets)
    end

    action = _normalize_ambiguous_action(ambiguous_action)
    if action == :error
        throw(ArgumentError("Ambiguous sequence type. Candidates: $(sort!(collect(possible_alphabets)))."))
    end

    preference = if action == :dna
        (:DNA, :RNA, :AA)
    elseif action == :rna
        (:RNA, :DNA, :AA)
    elseif action == :aa
        (:AA, :DNA, :RNA)
    else
        error("Invalid ambiguous_action: $ambiguous_action. Expected :dna, :rna, :aa, or :error.")
    end

    for candidate in preference
        if candidate in possible_alphabets
            return candidate
        end
    end

    throw(ArgumentError("Ambiguous sequence type. Candidates: $(sort!(collect(possible_alphabets)))."))
end

function _resolve_sequence_alphabet(
        possible_alphabets::Set{Symbol},
        unambiguous_alphabets::Set{Symbol};
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if !isnothing(type_hint)
        hint = _normalize_alphabet_symbol(type_hint)
        if !(hint in (:DNA, :RNA, :AA))
            error("Invalid type hint: $type_hint. Expected :DNA, :RNA, or :AA.")
        end
        if hint in possible_alphabets
            return hint
        end
        error("Type hint $hint is incompatible with inferred alphabets $(sort!(collect(possible_alphabets))).")
    end

    # If ambiguous characters expand the candidate alphabets, honor the preference order.
    if !isempty(unambiguous_alphabets) && unambiguous_alphabets == possible_alphabets
        return _resolve_ambiguous_alphabet(unambiguous_alphabets, ambiguous_action)
    end

    return _resolve_ambiguous_alphabet(possible_alphabets, ambiguous_action)
end

function _infer_sequence_alphabet(
        records::Vector{<:FASTX.Record};
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot infer sequence type from empty record set"))
    end

    inferred_alphabet = nothing
    unambiguous_alphabet = nothing
    for record in records
        seq_str = String(FASTX.sequence(record))

        record_alphabets = _possible_alphabets(seq_str)
        if isnothing(inferred_alphabet)
            inferred_alphabet = record_alphabets.possible
            unambiguous_alphabet = record_alphabets.unambiguous
        else
            inferred_alphabet = intersect(inferred_alphabet, record_alphabets.possible)
            unambiguous_alphabet = intersect(unambiguous_alphabet, record_alphabets.unambiguous)
        end

        if isempty(inferred_alphabet)
            error("Mixed sequence types detected. Cannot build graph from heterogeneous input.")
        end
    end

    return _resolve_sequence_alphabet(
        inferred_alphabet,
        unambiguous_alphabet;
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )
end

"""
Build strand-specific k-mer de Bruijn graph from FASTA/FASTQ records.

This is the PRIMARY graph construction function. Constructs a graph where:
- Each k-mer is stored AS OBSERVED (not canonicalized)
- Forward and reverse complement k-mers are SEPARATE vertices
- Evidence tracks which reads support which k-mers at which positions
- Edges connect k-mers with (k-1) overlap

**Outer wrapper**: Detects sequence type and dispatches to type-stable core function.

# Arguments
- `records::Vector{<:FASTX.Record}`: Input FASTA or FASTQ records
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `type_hint::Union{Nothing,Symbol}=nothing`: Optional alphabet hint (:DNA, :RNA, :AA)
- `ambiguous_action::Symbol=:dna`: Resolution for ambiguous alphabets (:dna, :rna, :aa, :error)

# Returns
- `MetaGraphsNext.MetaGraph`: Strand-specific k-mer graph with evidence
"""
function build_kmer_graph_singlestrand(
        records::Vector{<:FASTX.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    # Convert FASTQ to FASTA if needed (k-mer graphs don't use quality)
    fasta_records = if records[1] isa FASTX.FASTQ.Record
        [FASTX.FASTA.Record(
             FASTX.identifier(r),
             FASTX.sequence(r)
         ) for r in records]
    else
        records
    end

    sequence_alphabet = _infer_sequence_alphabet(
        fasta_records;
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    # Validate all records are compatible with the dataset alphabet.
    for record in fasta_records[1:min(10, length(fasta_records))]
        seq_str = String(FASTX.sequence(record))
        if !parentmodule(Rhizomorph).validate_alphabet(seq_str, sequence_alphabet)
            error("Mixed sequence types detected. Cannot build graph from heterogeneous input.")
        end
    end

    # Dispatch to type-stable core function
    if sequence_alphabet == :DNA
        return _build_kmer_graph_core(fasta_records, Val(k), Kmers.DNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :RNA
        return _build_kmer_graph_core(fasta_records, Val(k), Kmers.RNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :AA
        return _build_kmer_graph_core(fasta_records, Val(k), Kmers.AAKmer{k}, dataset_id)
    else
        error("Unsupported sequence type for k-mer graph: $sequence_alphabet")
    end
end

"""
Type-stable core function for k-mer graph construction.

# Type Parameters
- `KmerType`: Concrete k-mer type (DNAKmer{k}, RNAKmer{k}, or AAKmer{k})
"""
function _build_kmer_graph_core(
        records::Vector{FASTX.FASTA.Record},
        ::Val{K},
        ::Type{KmerType},
        dataset_id::String
) where {K, KmerType <: Kmers.Kmer}
    # Determine sequence type and iterator based on KmerType
    if KmerType <: Kmers.DNAKmer
        SeqType = BioSequences.LongDNA{4}
        KmerIterator = Kmers.UnambiguousDNAMers{K}
    elseif KmerType <: Kmers.RNAKmer
        SeqType = BioSequences.LongRNA{4}
        KmerIterator = Kmers.UnambiguousRNAMers{K}
    elseif KmerType <: Kmers.AAKmer
        SeqType = BioSequences.LongAA
        # Try UnambiguousAAMers first, fall back to FwAAMers if not available
        KmerIterator = try
            Kmers.UnambiguousAAMers{K}
        catch
            Kmers.FwAAMers{K}
        end
    else
        error("Unsupported k-mer type: $KmerType")
    end

    # Get actual kmer type from the first record that yields a k-mer.
    # Short/ambiguous inputs can yield no k-mers; fall back to the requested type.
    is_aa = KmerType <: Kmers.AAKmer
    actual_kmer_type = nothing
    for record in records
        test_seq = FASTX.sequence(SeqType, record)
        iter = KmerIterator(test_seq)
        first_item = iterate(iter)
        if first_item === nothing
            continue
        end
        value = first_item[1]
        test_kmer = is_aa ? value : value[1]
        actual_kmer_type = typeof(test_kmer)
        break
    end
    ActualKmerType = actual_kmer_type === nothing ? KmerType : actual_kmer_type

    # Create empty directed graph with actual kmer type as vertex labels
    # ALWAYS use directed graphs - MetaGraph with DiGraph backend
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = ActualKmerType,
        vertex_data_type = KmerVertexData{ActualKmerType},
        edge_data_type = KmerEdgeData,
        weight_function = compute_edge_weight
    )

    # Process each record
    for record in records
        observation_id = String(split(FASTX.identifier(record), ' ')[1])
        sequence = FASTX.sequence(SeqType, record)

        # Extract k-mers with positions
        # DNA/RNA iterators return (kmer, position), AA iterators return just kmer
        kmers_with_positions = if is_aa
            # FwAAMers returns just k-mers, need to manually add positions
            [(kmer, i) for (i, kmer) in enumerate(KmerIterator(sequence))]
        else
            # UnambiguousDNAMers/RNAMers return (kmer, position) tuples
            collect(KmerIterator(sequence))
        end

        # Add vertices and evidence
        for (kmer, position) in kmers_with_positions
            # Create vertex if it doesn't exist
            if !haskey(graph, kmer)
                vertex_data = KmerVertexData(kmer)
                graph[kmer] = vertex_data
            end

            # Add evidence to vertex
            vertex_data = graph[kmer]
            add_evidence!(vertex_data, dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        # Add edges between consecutive k-mers
        for i in 1:(length(kmers_with_positions) - 1)
            src_kmer, src_pos = kmers_with_positions[i]
            dst_kmer, dst_pos = kmers_with_positions[i + 1]

            # Only add edge if positions are consecutive (no gap from skipped ambiguous bases)
            if dst_pos == src_pos + 1
                # Create edge if it doesn't exist
                if !haskey(graph, src_kmer, dst_kmer)
                    edge_data = KmerEdgeData()
                    graph[src_kmer, dst_kmer] = edge_data
                end

                # Add evidence to edge
                edge_data = graph[src_kmer, dst_kmer]
                add_evidence!(edge_data, dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

"""
Build doublestrand k-mer de Bruijn graph from FASTA/FASTQ records.

Creates a directed graph containing both forward and reverse complement k-mers
as separate vertices. Each observed k-mer and its RC are added with properly
oriented evidence. Edges are directed with RC edges having reversed direction.

# Arguments
- `records::Vector{<:FASTX.Record}`: Input records
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking

# Returns
- `MetaGraphsNext.MetaGraph`: DoubleStrand k-mer graph (DiGraph with 2x vertices/edges)

# Example
```julia
records = [FASTX.FASTA.Record("seq1", dna"ATGATG")]
graph = build_kmer_graph_doublestrand(records, 3)

# Graph contains: ATG, TGA, TG, CAT (rc ATG), TCA (rc TGA), CA (rc AT)
# With directed edges in both orientations
```
"""
function build_kmer_graph_doublestrand(
        records::Vector{<:FASTX.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    # First build singlestrand graph
    singlestrand_graph = build_kmer_graph_singlestrand(
        records,
        k;
        dataset_id = dataset_id,
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    # Convert to doublestrand representation (replicate vertices/edges)
    return convert_to_doublestrand(singlestrand_graph)
end

# ============================================================================
# Qualmer Graph Construction (Quality-aware)
# ============================================================================

"""
Build strand-specific qualmer de Bruijn graph from FASTQ records.

Quality-aware version of k-mer graph construction. Stores quality scores
alongside k-mer observations.

**Outer wrapper**: Detects sequence type and dispatches to type-stable core function.

# Arguments
- `records::Vector{FASTX.FASTQ.Record}`: Input FASTQ records (MUST be FASTQ for quality)
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `type_hint::Union{Nothing,Symbol}=nothing`: Optional alphabet hint (:DNA, :RNA, :AA)
- `ambiguous_action::Symbol=:dna`: Resolution for ambiguous alphabets (:dna, :rna, :aa, :error)

# Returns
- `MetaGraphsNext.MetaGraph`: Strand-specific qualmer graph with quality-aware evidence
"""
function build_qualmer_graph_singlestrand(
        records::Vector{FASTX.FASTQ.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    sequence_alphabet = _infer_sequence_alphabet(
        records;
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    # Validate all records are compatible with the dataset alphabet.
    for record in records[1:min(10, length(records))]
        seq_str = String(FASTX.sequence(record))
        if !parentmodule(Rhizomorph).validate_alphabet(seq_str, sequence_alphabet)
            error("Mixed sequence types detected. Cannot build graph from heterogeneous input.")
        end
    end

    # Dispatch to type-stable core function
    if sequence_alphabet == :DNA
        return _build_qualmer_graph_core(records, Val(k), Kmers.DNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :RNA
        return _build_qualmer_graph_core(records, Val(k), Kmers.RNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :AA
        return _build_qualmer_graph_core(records, Val(k), Kmers.AAKmer{k}, dataset_id)
    else
        error("Unsupported sequence type for qualmer graph: $sequence_alphabet")
    end
end

"""
Type-stable core function for qualmer graph construction.

# Type Parameters
- `KmerType`: Concrete k-mer type (DNAKmer{k}, RNAKmer{k}, or AAKmer{k})
"""
function _build_qualmer_graph_core(
        records::Vector{FASTX.FASTQ.Record},
        ::Val{K},
        ::Type{KmerType},
        dataset_id::String
) where {K, KmerType <: Kmers.Kmer}
    # Determine sequence type and iterator
    if KmerType <: Kmers.DNAKmer
        SeqType = BioSequences.LongDNA{4}
        KmerIterator = Kmers.UnambiguousDNAMers{K}
    elseif KmerType <: Kmers.RNAKmer
        SeqType = BioSequences.LongRNA{4}
        KmerIterator = Kmers.UnambiguousRNAMers{K}
    elseif KmerType <: Kmers.AAKmer
        SeqType = BioSequences.LongAA
        KmerIterator = try
            Kmers.UnambiguousAAMers{K}
        catch
            Kmers.FwAAMers{K}
        end
    else
        error("Unsupported k-mer type: $KmerType")
    end

    # Get actual kmer type from iterator (includes all type parameters)
    test_seq = FASTX.sequence(SeqType, records[1])

    # DNA/RNA iterators return (kmer, position), AA iterators return just kmer
    is_aa = KmerType <: Kmers.AAKmer
    if is_aa
        test_kmer = first(KmerIterator(test_seq))
    else
        test_kmer, _ = first(KmerIterator(test_seq))
    end
    ActualKmerType = typeof(test_kmer)

    # Create empty directed graph with actual kmer type as vertex labels
    # ALWAYS use directed graphs - MetaGraph with DiGraph backend
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = ActualKmerType,
        vertex_data_type = QualmerVertexData{ActualKmerType},
        edge_data_type = QualmerEdgeData,
        weight_function = compute_edge_weight
    )

    # Process each record
    for record in records
        observation_id = String(split(FASTX.identifier(record), ' ')[1])
        sequence = FASTX.sequence(SeqType, record)
        quality = Vector{UInt8}(FASTX.quality(record))

        # Extract k-mers with positions
        # DNA/RNA iterators return (kmer, position), AA iterators return just kmer
        kmers_with_positions = if is_aa
            # FwAAMers returns just k-mers, need to manually add positions
            [(kmer, i) for (i, kmer) in enumerate(KmerIterator(sequence))]
        else
            # UnambiguousDNAMers/RNAMers return (kmer, position) tuples
            collect(KmerIterator(sequence))
        end

        # Add vertices and evidence
        for (kmer, position) in kmers_with_positions
            # Create vertex if it doesn't exist
            if !haskey(graph, kmer)
                vertex_data = QualmerVertexData(kmer)
                graph[kmer] = vertex_data
            end

            # Extract quality scores for this k-mer
            kmer_quality = quality[position:(position + K - 1)]

            # Add evidence to vertex
            vertex_data = graph[kmer]
            add_evidence!(vertex_data, dataset_id, observation_id,
                QualityEvidenceEntry(position, Forward, kmer_quality))
        end

        # Add edges between consecutive k-mers
        for i in 1:(length(kmers_with_positions) - 1)
            src_kmer, src_pos = kmers_with_positions[i]
            dst_kmer, dst_pos = kmers_with_positions[i + 1]

            # Only add edge if positions are consecutive
            if dst_pos == src_pos + 1
                # Create edge if it doesn't exist
                if !haskey(graph, src_kmer, dst_kmer)
                    edge_data = QualmerEdgeData()
                    graph[src_kmer, dst_kmer] = edge_data
                end

                # Extract quality scores for both k-mers
                src_quality = quality[src_pos:(src_pos + K - 1)]
                dst_quality = quality[dst_pos:(dst_pos + K - 1)]

                # Add evidence to edge
                edge_data = graph[src_kmer, dst_kmer]
                add_evidence!(edge_data, dataset_id, observation_id,
                    EdgeQualityEvidenceEntry(src_pos, dst_pos, Forward,
                        src_quality, dst_quality))
            end
        end
    end

    return graph
end

"""
Build doublestrand qualmer de Bruijn graph from FASTQ records.

Creates a directed graph containing both forward and reverse complement qualmers
as separate vertices, preserving quality information. Each observed qualmer and its
RC are added with properly oriented evidence and quality scores.

# Arguments
- `records::Vector{FASTX.FASTQ.Record}`: Input FASTQ records
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking

# Returns
- `MetaGraphsNext.MetaGraph`: DoubleStrand qualmer graph (DiGraph with 2x vertices/edges)
"""
function build_qualmer_graph_doublestrand(
        records::Vector{FASTX.FASTQ.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    # First build singlestrand qualmer graph
    singlestrand_graph = build_qualmer_graph_singlestrand(
        records,
        k;
        dataset_id = dataset_id,
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    # Convert to doublestrand representation (replicate vertices/edges)
    return convert_to_doublestrand(singlestrand_graph)
end

# ============================================================================
# Incremental Graph Construction
# ============================================================================

"""
Add observations from new records to existing graph.

Allows incremental graph construction from multiple datasets or batches.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Existing graph to add observations to
- `records::Vector{<:FASTX.Record}`: New records to add (FASTA for kmer, FASTQ for qualmer)
- `k::Int`: K-mer size (must match graph)
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking

# Returns
- `MetaGraphsNext.MetaGraph`: Updated graph with new observations
"""
function add_observations_to_graph!(
        graph::MetaGraphsNext.MetaGraph,
        records::Vector{<:FASTX.Record},
        k::Int;
        dataset_id::String = "dataset_01"
)
    if isempty(records)
        return graph
    end

    # Get actual k-mer type and vertex data type from the graph
    # This is safe and doesn't rely on brittle parameter indexing
    if !isempty(MetaGraphsNext.labels(graph))
        KmerType = typeof(first(MetaGraphsNext.labels(graph)))
        first_label = first(MetaGraphsNext.labels(graph))
        vertex_data = graph[first_label]
        is_quality_aware = vertex_data isa QualmerVertexData
    else
        # Empty graph - need to infer from first record
        error("Cannot add observations to empty graph. Use build_*_graph_singlestrand instead.")
    end

    # Validate record type matches graph type
    if is_quality_aware && !(records[1] isa FASTX.FASTQ.Record)
        error("Quality-aware graph requires FASTQ records")
    end

    # Convert FASTQ to FASTA for k-mer graphs
    if !is_quality_aware && records[1] isa FASTX.FASTQ.Record
        records = [FASTX.FASTA.Record(
                       FASTX.identifier(r),
                       FASTX.sequence(r)
                   ) for r in records]
    end

    # Determine sequence type and iterator based on KmerType
    if KmerType <: Kmers.DNAKmer
        SeqType = BioSequences.LongDNA{4}
        KmerIterator = Kmers.UnambiguousDNAMers{k}
    elseif KmerType <: Kmers.RNAKmer
        SeqType = BioSequences.LongRNA{4}
        KmerIterator = Kmers.UnambiguousRNAMers{k}
    elseif KmerType <: Kmers.AAKmer
        SeqType = BioSequences.LongAA
        KmerIterator = try
            Kmers.UnambiguousAAMers{k}
        catch
            Kmers.FwAAMers{k}
        end
    else
        error("Unsupported k-mer type in graph: $KmerType")
    end

    # Process each record
    for record in records
        observation_id = String(split(FASTX.identifier(record), ' ')[1])
        sequence = FASTX.sequence(SeqType, record)

        if is_quality_aware
            quality = Vector{UInt8}(FASTX.quality(record))
        end

        # Extract k-mers with positions
        kmers_with_positions = collect(KmerIterator(sequence))

        # Add vertices and evidence
        for (kmer, position) in kmers_with_positions
            # Create vertex if it doesn't exist
            if !haskey(graph, kmer)
                if is_quality_aware
                    vertex_data = QualmerVertexData(kmer)
                else
                    vertex_data = KmerVertexData(kmer)
                end
                graph[kmer] = vertex_data
            end

            # Add evidence to vertex
            vertex_data = graph[kmer]

            if is_quality_aware
                kmer_quality = quality[position:(position + k - 1)]
                add_evidence!(vertex_data, dataset_id, observation_id,
                    QualityEvidenceEntry(position, Forward, kmer_quality))
            else
                add_evidence!(vertex_data, dataset_id, observation_id,
                    EvidenceEntry(position, Forward))
            end
        end

        # Add edges between consecutive k-mers
        for i in 1:(length(kmers_with_positions) - 1)
            src_kmer, src_pos = kmers_with_positions[i]
            dst_kmer, dst_pos = kmers_with_positions[i + 1]

            # Only add edge if positions are consecutive
            if dst_pos == src_pos + 1
                # Create edge if it doesn't exist
                if !haskey(graph, src_kmer, dst_kmer)
                    if is_quality_aware
                        edge_data = QualmerEdgeData()
                    else
                        edge_data = KmerEdgeData()
                    end
                    graph[src_kmer, dst_kmer] = edge_data
                end

                # Add evidence to edge
                edge_data = graph[src_kmer, dst_kmer]

                if is_quality_aware
                    src_quality = quality[src_pos:(src_pos + k - 1)]
                    dst_quality = quality[dst_pos:(dst_pos + k - 1)]
                    add_evidence!(edge_data, dataset_id, observation_id,
                        EdgeQualityEvidenceEntry(src_pos, dst_pos, Forward,
                            src_quality, dst_quality))
                else
                    add_evidence!(edge_data, dataset_id, observation_id,
                        EdgeEvidenceEntry(src_pos, dst_pos, Forward))
                end
            end
        end
    end

    return graph
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
Get observation ID from FASTX record.

Uses read identifier (sequence name) up to first whitespace.

# Example
- @SRR123456.1 1 length=150  →  "SRR123456.1"
"""
function get_observation_id(record::Union{FASTX.FASTA.Record, FASTX.FASTQ.Record})
    identifier = FASTX.identifier(record)
    return String(split(identifier, ' ')[1])
end

"""
Generate observation ID from index (for synthetic data or non-unique read names).
"""
function get_observation_id_from_index(index::Int, prefix::String = "obs")
    return "$(prefix)_$(lpad(index, 10, '0'))"
end

"""
Extract dataset ID from filepath.

DEFAULT: Use input filename (unless user specifies alternative).

# Examples
1. Single file: dataset_id = "sample_A"  (from sample_A.fastq)
2. User override: dataset_id = "patient_001_tumor"  (explicit metadata)
"""
function get_dataset_id_from_file(filepath::String)
    return splitext(basename(filepath))[1]
end

function _alphabet_hint_from_path(filepath::String)
    return parentmodule(Rhizomorph).alphabet_hint_from_path(filepath)
end

# ============================================================================
# Canonical Graph Construction (Undirected, Merged RC Pairs)
# ============================================================================

"""
Build canonical (doublestrand) k-mer graph from FASTA/FASTQ records.

This constructs a graph where:
- Each canonical k-mer represents BOTH forward and reverse complement
- Evidence from both strands is merged into the canonical k-mer
- Suitable for double-stranded DNA/RNA (not amino acids)

# Arguments
- `records::Vector{<:FASTX.Record}`: Input FASTA or FASTQ records
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking

# Returns
- `MetaGraphsNext.MetaGraph`: Canonical k-mer graph with merged evidence

# Examples
```julia
# DNA reads
records = [FASTX.FASTA.Record("read1", "ATGCAT")]
graph = build_kmer_graph_canonical(records, 3)

# Both ATG and its RC (CAT) map to the canonical form
```
"""
function build_kmer_graph_canonical(
        records::Vector{<:FASTX.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    # First build singlestrand graph
    singlestrand_graph = build_kmer_graph_singlestrand(
        records,
        k;
        dataset_id = dataset_id,
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    # Convert to canonical representation
    return convert_to_canonical(singlestrand_graph)
end

"""
Build canonical (doublestrand) qualmer graph from FASTQ records.

Similar to `build_kmer_graph_canonical` but preserves quality information.

# Arguments
- `records::Vector{FASTX.FASTQ.Record}`: Input FASTQ records
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier

# Returns
- `MetaGraphsNext.MetaGraph`: Canonical qualmer graph with merged evidence
"""
function build_qualmer_graph_canonical(
        records::Vector{FASTX.FASTQ.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    # First build singlestrand qualmer graph
    singlestrand_graph = build_qualmer_graph_singlestrand(
        records,
        k;
        dataset_id = dataset_id,
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    # Convert to canonical representation
    return convert_to_canonical(singlestrand_graph)
end

"""
    convert_to_doublestrand(singlestrand_graph)

Convert a singlestrand graph to doublestrand representation.

Replicates all vertices and edges in both forward and reverse complement orientations.
Creates a DIRECTED graph with 2x vertices and 2x directed edges.

# Arguments
- `singlestrand_graph`: Graph with strand-specific k-mers

# Returns
- New graph with both forward and RC k-mers as separate vertices (directed edges)

# Example
```julia
# If singlestrand has: ATG→TGC
# Doublestrand will have:
#   Vertices: {ATG, TGC, CAT (rc of ATG), GCA (rc of TGC)}
#   Edges: ATG→TGC, GCA→CAT (both directed)
```
"""
function convert_to_doublestrand(singlestrand_graph)
    # Get graph type information
    if isempty(MetaGraphsNext.labels(singlestrand_graph))
        return singlestrand_graph  # Empty graph stays empty
    end

    first_label = first(MetaGraphsNext.labels(singlestrand_graph))
    first_vertex = singlestrand_graph[first_label]
    KmerType = typeof(first_label)
    VertexDataType = typeof(first_vertex)

    # Only DNA/RNA k-mers support doublestrand conversion
    if !(KmerType <: Kmers.DNAKmer || KmerType <: Kmers.RNAKmer)
        error("Cannot create doublestrand graph for non-nucleotide k-mer labels: $(KmerType)")
    end
    # Amino acids don't have reverse complement
    if KmerType <: Kmers.AAKmer
        error("Cannot create doublestrand graph for amino acid sequences (no reverse complement)")
    end

    # Check if this is a reduced (non-full) graph
    is_reduced = !hasfield(VertexDataType, :evidence)

    if is_reduced
        return _convert_to_doublestrand_reduced(singlestrand_graph)
    end

    # Determine edge data type
    EdgeDataType = if !isempty(MetaGraphsNext.edge_labels(singlestrand_graph))
        src, dst = first(MetaGraphsNext.edge_labels(singlestrand_graph))
        typeof(singlestrand_graph[src, dst])
    else
        KmerEdgeData
    end

    # Create new doublestrand graph (still DiGraph!)
    doublestrand_graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = KmerType,
        vertex_data_type = VertexDataType,
        edge_data_type = EdgeDataType,
        weight_function = compute_edge_weight
    )

    # Add all forward vertices
    for kmer in MetaGraphsNext.labels(singlestrand_graph)
        doublestrand_graph[kmer] = singlestrand_graph[kmer]
    end

    # Add all reverse complement vertices (merging evidence if RC already exists)
    for kmer in MetaGraphsNext.labels(singlestrand_graph)
        rc_kmer = BioSequences.reverse_complement(kmer)

        # Skip if RC is same as forward (palindrome) - already added
        if rc_kmer == kmer
            continue
        end

        # Get the forward vertex data to flip
        vertex_data = singlestrand_graph[kmer]

        # Check if RC vertex already exists (it may have been in the original graph)
        if haskey(doublestrand_graph, rc_kmer)
            # RC vertex exists - merge flipped evidence into existing vertex
            existing_rc_vertex = doublestrand_graph[rc_kmer]
            for (dataset_id, observations) in vertex_data.evidence
                for (obs_id, evidence_set) in observations
                    for evidence in evidence_set
                        flipped_evidence = flip_evidence_strand(evidence)
                        add_evidence!(existing_rc_vertex, dataset_id, obs_id, flipped_evidence)
                    end
                end
            end
        else
            # RC vertex doesn't exist - create new one
            rc_vertex_data = if VertexDataType <: KmerVertexData
                KmerVertexData(rc_kmer)
            elseif VertexDataType <: QualmerVertexData
                QualmerVertexData(rc_kmer)
            else
                error("Unsupported vertex data type: $(VertexDataType)")
            end

            # Copy and flip all evidence from forward to RC
            for (dataset_id, observations) in vertex_data.evidence
                for (obs_id, evidence_set) in observations
                    for evidence in evidence_set
                        flipped_evidence = flip_evidence_strand(evidence)
                        add_evidence!(rc_vertex_data, dataset_id, obs_id, flipped_evidence)
                    end
                end
            end

            doublestrand_graph[rc_kmer] = rc_vertex_data
        end
    end

    # Add all forward edges
    for (src_kmer, dst_kmer) in MetaGraphsNext.edge_labels(singlestrand_graph)
        doublestrand_graph[src_kmer, dst_kmer] = singlestrand_graph[src_kmer, dst_kmer]
    end

    # Add all reverse complement edges (with reversed direction!)
    for (src_kmer, dst_kmer) in MetaGraphsNext.edge_labels(singlestrand_graph)
        # RC edge: reverse_complement(dst) → reverse_complement(src)
        # Note the REVERSED direction!
        rc_src = BioSequences.reverse_complement(dst_kmer)
        rc_dst = BioSequences.reverse_complement(src_kmer)

        # Create RC edge with flipped evidence
        edge_data = singlestrand_graph[src_kmer, dst_kmer]

        # Create new edge data for RC
        rc_edge_data = if EdgeDataType <: KmerEdgeData
            KmerEdgeData()
        elseif EdgeDataType <: QualmerEdgeData
            QualmerEdgeData()
        else
            error("Unsupported edge data type: $(EdgeDataType)")
        end

        # Copy and flip all evidence from forward edge to RC edge
        for (dataset_id, observations) in edge_data.evidence
            for (obs_id, evidence_set) in observations
                for evidence in evidence_set
                    flipped_evidence = flip_evidence_strand(evidence)
                    add_evidence!(rc_edge_data, dataset_id, obs_id, flipped_evidence)
                end
            end
        end

        doublestrand_graph[rc_src, rc_dst] = rc_edge_data
    end

    return doublestrand_graph
end

"""
    convert_to_canonical(singlestrand_graph)

Convert a singlestrand graph to canonical representation.

Merges forward and reverse complement k-mers into their canonical forms.
Evidence from both strands is combined. Creates an UNDIRECTED graph where
edges represent bidirectional relationships.

Note: Canonical graphs require on-the-fly validation of traversal directions
using find_valid_canonical_traversal() since edges are undirected.

# Arguments
- `singlestrand_graph`: Graph with strand-specific k-mers

# Returns
- New graph with canonical k-mers and merged evidence (undirected edges)
"""
function convert_to_canonical(singlestrand_graph)
    # Get graph type information from first vertex
    if isempty(MetaGraphsNext.labels(singlestrand_graph))
        return singlestrand_graph  # Empty graph stays empty
    end

    first_label = first(MetaGraphsNext.labels(singlestrand_graph))
    first_vertex = singlestrand_graph[first_label]
    KmerType = typeof(first_label)
    VertexDataType = typeof(first_vertex)

    # Only DNA/RNA k-mers support canonical conversion
    if !(KmerType <: Kmers.DNAKmer || KmerType <: Kmers.RNAKmer)
        error("Cannot create canonical graph for non-nucleotide k-mer labels: $(KmerType)")
    end
    # Amino acids don't have reverse complement
    if KmerType <: Kmers.AAKmer
        error("Cannot create doublestrand graph for amino acid sequences (no reverse complement)")
    end

    # Check if this is a reduced (non-full) graph
    is_reduced = !hasfield(VertexDataType, :evidence)

    if is_reduced
        return _convert_to_canonical_reduced(singlestrand_graph)
    end

    # Determine edge data type
    EdgeDataType = if !isempty(MetaGraphsNext.edge_labels(singlestrand_graph))
        src, dst = first(MetaGraphsNext.edge_labels(singlestrand_graph))
        typeof(singlestrand_graph[src, dst])
    else
        # Default to KmerEdgeData if no edges
        KmerEdgeData
    end

    # Create new canonical graph
    canonical_graph = MetaGraphsNext.MetaGraph(
        Graphs.Graph();
        label_type = KmerType,
        vertex_data_type = VertexDataType,
        edge_data_type = EdgeDataType,
        weight_function = compute_edge_weight
    )

    # Track which k-mers we've already processed
    processed = Set{KmerType}()

    # Process all vertices
    for kmer in MetaGraphsNext.labels(singlestrand_graph)
        if kmer in processed
            continue
        end

        # Get canonical form
        canon_kmer = BioSequences.canonical(kmer)
        rc_kmer = BioSequences.reverse_complement(kmer)

        # Mark both as processed
        push!(processed, kmer)
        push!(processed, rc_kmer)

        # Get vertex data for forward strand
        fwd_data = singlestrand_graph[kmer]

        # Create or get canonical vertex
        if !haskey(canonical_graph, canon_kmer)
            # Create new vertex with canonical k-mer
            if VertexDataType <: KmerVertexData
                canonical_graph[canon_kmer] = KmerVertexData(canon_kmer)
            else  # QualmerVertexData
                canonical_graph[canon_kmer] = QualmerVertexData(canon_kmer)
            end
        end

        canon_vertex = canonical_graph[canon_kmer]

        # Add forward strand evidence
        for (dataset_id, observations) in fwd_data.evidence
            for (obs_id, evidence_set) in observations
                for evidence in evidence_set
                    # If kmer != canon_kmer, we need to flip the strand
                    final_evidence = if kmer != canon_kmer
                        flip_evidence_strand(evidence)
                    else
                        evidence
                    end
                    add_evidence!(canon_vertex, dataset_id, obs_id, final_evidence)
                end
            end
        end

        # Add reverse complement evidence if it exists
        if haskey(singlestrand_graph, rc_kmer)
            rc_data = singlestrand_graph[rc_kmer]
            for (dataset_id, observations) in rc_data.evidence
                for (obs_id, evidence_set) in observations
                    for evidence in evidence_set
                        # Flip strand for RC evidence
                        final_evidence = if rc_kmer != canon_kmer
                            flip_evidence_strand(evidence)
                        else
                            evidence
                        end
                        add_evidence!(canon_vertex, dataset_id, obs_id, final_evidence)
                    end
                end
            end
        end
    end

    # Process edges (similar approach)
    for (src_kmer, dst_kmer) in MetaGraphsNext.edge_labels(singlestrand_graph)
        canon_src = BioSequences.canonical(src_kmer)
        canon_dst = BioSequences.canonical(dst_kmer)

        # Create edge if it doesn't exist
        if !haskey(canonical_graph, canon_src, canon_dst)
            if EdgeDataType == KmerEdgeData
                canonical_graph[canon_src, canon_dst] = KmerEdgeData()
            else
                canonical_graph[canon_src, canon_dst] = QualmerEdgeData()
            end
        end

        # Add edge evidence
        edge_data = singlestrand_graph[src_kmer, dst_kmer]
        canon_edge = canonical_graph[canon_src, canon_dst]

        for (dataset_id, observations) in edge_data.evidence
            for (obs_id, evidence_set) in observations
                for evidence in evidence_set
                    # Flip strand if needed
                    needs_flip = (src_kmer != canon_src) || (dst_kmer != canon_dst)
                    final_evidence = if needs_flip
                        flip_evidence_strand(evidence)
                    else
                        evidence
                    end
                    add_evidence!(canon_edge, dataset_id, obs_id, final_evidence)
                end
            end
        end
    end

    return canonical_graph
end

"""
    convert_to_singlestrand(doublestrand_graph)

Convert a doublestrand (canonical) graph back to singlestrand form.

This "unfolds" canonical k-mers into both forward and reverse complement forms,
creating separate vertices for each orientation with properly oriented evidence.

# Arguments
- `doublestrand_graph`: Graph with canonical k-mers

# Returns
- Singlestrand graph with separate vertices for forward and RC k-mers

# Examples
```julia
# Convert canonical graph back to singlestrand
ss_graph = convert_to_singlestrand(ds_graph)
```
"""
function convert_to_singlestrand(doublestrand_graph)
    # Determine graph type from first vertex
    first_kmer = first(MetaGraphsNext.labels(doublestrand_graph))
    first_vertex = doublestrand_graph[first_kmer]

    # Only DNA/RNA k-mers support doublestrand/canonical modes
    if !(first_kmer isa Kmers.DNAKmer || first_kmer isa Kmers.RNAKmer)
        error("Cannot convert non-nucleotide canonical graph to singlestrand: $(typeof(first_kmer))")
    end

    VertexDataType = typeof(first_vertex)
    EdgeDataType = if VertexDataType <: KmerVertexData
        KmerEdgeData
    else
        QualmerEdgeData
    end

    # Create new singlestrand graph
    singlestrand_graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = typeof(first_kmer),
        vertex_data_type = VertexDataType,
        edge_data_type = EdgeDataType,
        graph_data = "Singlestrand graph (unfolded from canonical)",
        weight_function = compute_edge_weight
    )

    # Process each canonical k-mer
    for canon_kmer in MetaGraphsNext.labels(doublestrand_graph)
        rc_kmer = BioSequences.reverse_complement(canon_kmer)
        canon_vertex = doublestrand_graph[canon_kmer]

        # Create forward k-mer vertex if it doesn't exist
        if !haskey(singlestrand_graph, canon_kmer)
            if VertexDataType <: KmerVertexData
                singlestrand_graph[canon_kmer] = KmerVertexData(canon_kmer)
            else
                singlestrand_graph[canon_kmer] = QualmerVertexData(canon_kmer)
            end
        end

        # Create RC k-mer vertex if it doesn't exist and is different from canonical
        if rc_kmer != canon_kmer && !haskey(singlestrand_graph, rc_kmer)
            if VertexDataType <: KmerVertexData
                singlestrand_graph[rc_kmer] = KmerVertexData(rc_kmer)
            else
                singlestrand_graph[rc_kmer] = QualmerVertexData(rc_kmer)
            end
        end

        # Distribute evidence to appropriate vertices
        for (dataset_id, observations) in canon_vertex.evidence
            for (obs_id, evidence_set) in observations
                for evidence in evidence_set
                    # Determine which vertex to add evidence to based on strand
                    if evidence.strand == Forward
                        # Forward evidence goes to canonical k-mer
                        add_evidence!(singlestrand_graph[canon_kmer], dataset_id, obs_id, evidence)
                    else
                        # Reverse evidence goes to RC k-mer (flipped to Forward)
                        if rc_kmer != canon_kmer
                            flipped_evidence = flip_evidence_strand(evidence)
                            add_evidence!(singlestrand_graph[rc_kmer], dataset_id, obs_id, flipped_evidence)
                        else
                            # Palindrome - keep as is
                            add_evidence!(singlestrand_graph[canon_kmer], dataset_id, obs_id, evidence)
                        end
                    end
                end
            end
        end
    end

    # Process edges
    for (canon_src, canon_dst) in MetaGraphsNext.edge_labels(doublestrand_graph)
        edge_data = doublestrand_graph[canon_src, canon_dst]

        rc_src = BioSequences.reverse_complement(canon_src)
        rc_dst = BioSequences.reverse_complement(canon_dst)

        # Distribute edge evidence based on strand
        for (dataset_id, observations) in edge_data.evidence
            for (obs_id, evidence_set) in observations
                for evidence in evidence_set
                    if evidence.strand == Forward
                        # Forward edge: canon_src → canon_dst
                        if !haskey(singlestrand_graph, canon_src, canon_dst)
                            if EdgeDataType <: KmerEdgeData
                                singlestrand_graph[canon_src, canon_dst] = KmerEdgeData()
                            else
                                singlestrand_graph[canon_src, canon_dst] = QualmerEdgeData()
                            end
                        end
                        add_evidence!(singlestrand_graph[canon_src, canon_dst],
                            dataset_id, obs_id, evidence)
                    else
                        # Reverse edge: rc_dst → rc_src (note the reversal!)
                        # When unfolding, reverse strand means we're on the RC sequence
                        # So the edge goes in the opposite direction on the RC strand
                        if !haskey(singlestrand_graph, rc_dst, rc_src)
                            if EdgeDataType <: KmerEdgeData
                                singlestrand_graph[rc_dst, rc_src] = KmerEdgeData()
                            else
                                singlestrand_graph[rc_dst, rc_src] = QualmerEdgeData()
                            end
                        end
                        # Flip evidence to Forward since it's now on the RC sequence
                        flipped_evidence = flip_evidence_strand(evidence)
                        add_evidence!(singlestrand_graph[rc_dst, rc_src],
                            dataset_id, obs_id, flipped_evidence)
                    end
                end
            end
        end
    end

    return singlestrand_graph
end

# ============================================================================
# N-gram Graph Construction (String-based)
# ============================================================================

"""
    build_ngram_graph_singlestrand(strings, n; dataset_id="dataset_01")

Build strand-specific n-gram de Bruijn graph from Unicode strings.

This is the string equivalent of k-mer graphs, working with fixed-length
character n-grams instead of biological k-mers.

# Arguments
- `strings::Vector{String}`: Input strings (or vector of string records)
- `n::Int`: N-gram size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking

# Returns
- `MetaGraphsNext.MetaGraph`: Strand-specific n-gram graph with evidence

# Key Differences from K-mer Graphs
- No reverse complement (strings don't have biological RC)
- All n-grams are "Forward" strand by convention
- Labels are String type, not k-mer types
- Supports full Unicode, not just biological alphabets

# Examples
```julia
# Simple text analysis
texts = ["hello world", "world hello", "hello there"]
graph = build_ngram_graph_singlestrand(texts, 3)

# Find common 3-character patterns
stats = get_ngram_statistics(graph)
```
"""
function build_ngram_graph_singlestrand(
        strings::Vector{String},
        n::Int;
        dataset_id::String = "dataset_01"
)
    if isempty(strings)
        error("Cannot build graph from empty string vector")
    end

    if n < 1
        error("N-gram size must be positive, got: $n")
    end

    # Create empty directed graph with String labels
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = StringVertexData,
        edge_data_type = StringEdgeData,
        weight_function = compute_edge_weight
    )

    # Process each string
    for (string_idx, input_string) in enumerate(strings)
        observation_id = "string_$(lpad(string_idx, 6, '0'))"

        chars = collect(input_string)
        char_count = length(chars)

        # Skip strings shorter than n (by character count)
        if char_count < n
            continue
        end

        # Extract n-grams with positions (1-indexed)
        ngrams_with_positions = [(String(chars[i:(i + n - 1)]), i)
                                 for i in 1:(char_count - n + 1)]

        # Add vertices and evidence
        for (ngram, position) in ngrams_with_positions
            # Create vertex if it doesn't exist
            if !haskey(graph, ngram)
                vertex_data = StringVertexData(ngram)
                graph[ngram] = vertex_data
            end

            # Add evidence to vertex
            vertex_data = graph[ngram]
            add_evidence!(vertex_data, dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        # Add edges between consecutive n-grams (n-1 character overlap)
        for i in 1:(length(ngrams_with_positions) - 1)
            src_ngram, src_pos = ngrams_with_positions[i]
            dst_ngram, dst_pos = ngrams_with_positions[i + 1]

            # Verify they are actually consecutive (position difference = 1)
            # This ensures we only connect n-grams with proper (n-1) overlap
            if dst_pos == src_pos + 1
                # Create edge if it doesn't exist
                if !haskey(graph, src_ngram, dst_ngram)
                    # N-grams have (n-1) character overlap
                    edge_data = StringEdgeData(n - 1)
                    graph[src_ngram, dst_ngram] = edge_data
                end

                # Add evidence to edge
                edge_data = graph[src_ngram, dst_ngram]
                add_evidence!(edge_data, dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

"""
    add_observations_to_ngram_graph!(graph, strings, n; dataset_id="dataset_01")

Add observations from new strings to existing n-gram graph.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Existing n-gram graph
- `strings::Vector{String}`: New strings to add
- `n::Int`: N-gram size (must match graph)
- `dataset_id::String="dataset_01"`: Dataset identifier

# Returns
- Updated graph with new observations
"""
function add_observations_to_ngram_graph!(
        graph::MetaGraphsNext.MetaGraph,
        strings::Vector{String},
        n::Int;
        dataset_id::String = "dataset_01"
)
    if isempty(strings)
        return graph
    end

    # Process each string
    for (string_idx, input_string) in enumerate(strings)
        observation_id = "string_$(lpad(string_idx, 6, '0'))"

        # Skip strings shorter than n
        if length(input_string) < n
            continue
        end

        # Extract n-grams with positions
        ngrams_with_positions = [(input_string[i:(i + n - 1)], i)
                                 for i in 1:(length(input_string) - n + 1)]

        # Add vertices and evidence
        for (ngram, position) in ngrams_with_positions
            # Create vertex if it doesn't exist
            if !haskey(graph, ngram)
                vertex_data = StringVertexData(ngram)
                graph[ngram] = vertex_data
            end

            # Add evidence to vertex
            vertex_data = graph[ngram]
            add_evidence!(vertex_data, dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        # Add edges between consecutive n-grams
        for i in 1:(length(ngrams_with_positions) - 1)
            src_ngram, src_pos = ngrams_with_positions[i]
            dst_ngram, dst_pos = ngrams_with_positions[i + 1]

            if dst_pos == src_pos + 1
                # Create edge if it doesn't exist
                if !haskey(graph, src_ngram, dst_ngram)
                    # N-grams have (n-1) character overlap
                    edge_data = StringEdgeData(n - 1)
                    graph[src_ngram, dst_ngram] = edge_data
                end

                # Add evidence to edge
                edge_data = graph[src_ngram, dst_ngram]
                add_evidence!(edge_data, dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

# ============================================================================
# Lightweight Graph Construction
#
# These builders create graphs with LightweightVertexData/LightweightEdgeData
# instead of full evidence-tracking types. The graph topology is identical;
# only per-observation evidence storage is replaced by integer counters.
# ============================================================================

"""
Build lightweight k-mer graph that stores only topology + aggregated counts.

Core loop mirrors `_build_kmer_graph_core` but creates `LightweightKmerVertexData`
vertices and `LightweightEdgeData` edges, reducing memory by ~1875x for
high-observation datasets.
"""
function _build_kmer_graph_core_lightweight(
        records::Vector{FASTX.FASTA.Record},
        ::Val{K},
        ::Type{KmerType},
        dataset_id::String
) where {K, KmerType <: Kmers.Kmer}
    if KmerType <: Kmers.DNAKmer
        SeqType = BioSequences.LongDNA{4}
        KmerIterator = Kmers.UnambiguousDNAMers{K}
    elseif KmerType <: Kmers.RNAKmer
        SeqType = BioSequences.LongRNA{4}
        KmerIterator = Kmers.UnambiguousRNAMers{K}
    elseif KmerType <: Kmers.AAKmer
        SeqType = BioSequences.LongAA
        KmerIterator = try
            Kmers.UnambiguousAAMers{K}
        catch
            Kmers.FwAAMers{K}
        end
    else
        error("Unsupported k-mer type: $KmerType")
    end

    is_aa = KmerType <: Kmers.AAKmer
    actual_kmer_type = nothing
    for record in records
        test_seq = FASTX.sequence(SeqType, record)
        iter = KmerIterator(test_seq)
        first_item = iterate(iter)
        if first_item === nothing
            continue
        end
        value = first_item[1]
        test_kmer = is_aa ? value : value[1]
        actual_kmer_type = typeof(test_kmer)
        break
    end
    ActualKmerType = actual_kmer_type === nothing ? KmerType : actual_kmer_type

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = ActualKmerType,
        vertex_data_type = LightweightKmerVertexData{ActualKmerType},
        edge_data_type = LightweightEdgeData,
        weight_function = compute_edge_weight
    )

    for record in records
        observation_id = String(split(FASTX.identifier(record), ' ')[1])
        sequence = FASTX.sequence(SeqType, record)

        kmers_with_positions = if is_aa
            [(kmer, i) for (i, kmer) in enumerate(KmerIterator(sequence))]
        else
            collect(KmerIterator(sequence))
        end

        for (kmer, position) in kmers_with_positions
            if !haskey(graph, kmer)
                graph[kmer] = LightweightKmerVertexData(kmer)
            end
            add_evidence!(graph[kmer], dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        for i in 1:(length(kmers_with_positions) - 1)
            src_kmer, src_pos = kmers_with_positions[i]
            dst_kmer, dst_pos = kmers_with_positions[i + 1]

            if dst_pos == src_pos + 1
                if !haskey(graph, src_kmer, dst_kmer)
                    graph[src_kmer, dst_kmer] = LightweightEdgeData()
                end
                add_evidence!(graph[src_kmer, dst_kmer], dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

"""
Build lightweight strand-specific k-mer graph (outer wrapper with type detection).
"""
function build_kmer_graph_singlestrand_lightweight(
        records::Vector{<:FASTX.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    fasta_records = if records[1] isa FASTX.FASTQ.Record
        [FASTX.FASTA.Record(
             FASTX.identifier(r),
             FASTX.sequence(r)
         ) for r in records]
    else
        records
    end

    sequence_alphabet = _infer_sequence_alphabet(
        fasta_records;
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    if sequence_alphabet == :DNA
        return _build_kmer_graph_core_lightweight(
            fasta_records, Val(k), Kmers.DNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :RNA
        return _build_kmer_graph_core_lightweight(
            fasta_records, Val(k), Kmers.RNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :AA
        return _build_kmer_graph_core_lightweight(
            fasta_records, Val(k), Kmers.AAKmer{k}, dataset_id)
    else
        error("Unsupported sequence type for k-mer graph: $sequence_alphabet")
    end
end

"""
Build lightweight n-gram graph that stores only topology + aggregated counts.

Core loop mirrors `build_ngram_graph_singlestrand` but creates
`LightweightStringVertexData` vertices and `LightweightEdgeData` edges.
"""
function build_ngram_graph_singlestrand_lightweight(
        strings::Vector{String},
        n::Int;
        dataset_id::String = "dataset_01"
)
    if isempty(strings)
        error("Cannot build graph from empty string vector")
    end

    if n < 1
        error("N-gram size must be positive, got: $n")
    end

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = LightweightStringVertexData,
        edge_data_type = LightweightEdgeData,
        weight_function = compute_edge_weight
    )

    for (string_idx, input_string) in enumerate(strings)
        observation_id = "string_$(lpad(string_idx, 6, '0'))"

        chars = collect(input_string)
        char_count = length(chars)

        if char_count < n
            continue
        end

        ngrams_with_positions = [(String(chars[i:(i + n - 1)]), i)
                                 for i in 1:(char_count - n + 1)]

        for (ngram, position) in ngrams_with_positions
            if !haskey(graph, ngram)
                graph[ngram] = LightweightStringVertexData(ngram)
            end
            add_evidence!(graph[ngram], dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        for i in 1:(length(ngrams_with_positions) - 1)
            src_ngram, src_pos = ngrams_with_positions[i]
            dst_ngram, dst_pos = ngrams_with_positions[i + 1]

            if dst_pos == src_pos + 1
                if !haskey(graph, src_ngram, dst_ngram)
                    graph[src_ngram, dst_ngram] = LightweightEdgeData(n - 1)
                end
                add_evidence!(graph[src_ngram, dst_ngram], dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

"""
    add_observations_to_ngram_graph_lightweight!(graph, strings, n; dataset_id)

Add observations from new strings to an existing lightweight n-gram graph.
"""
function add_observations_to_ngram_graph_lightweight!(
        graph::MetaGraphsNext.MetaGraph,
        strings::Vector{String},
        n::Int;
        dataset_id::String = "dataset_01"
)
    if isempty(strings)
        return graph
    end

    for (string_idx, input_string) in enumerate(strings)
        observation_id = "string_$(lpad(string_idx, 6, '0'))"

        if length(input_string) < n
            continue
        end

        ngrams_with_positions = [(input_string[i:(i + n - 1)], i)
                                 for i in 1:(length(input_string) - n + 1)]

        for (ngram, position) in ngrams_with_positions
            if !haskey(graph, ngram)
                graph[ngram] = LightweightStringVertexData(ngram)
            end
            add_evidence!(graph[ngram], dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        for i in 1:(length(ngrams_with_positions) - 1)
            src_ngram, src_pos = ngrams_with_positions[i]
            dst_ngram, dst_pos = ngrams_with_positions[i + 1]

            if dst_pos == src_pos + 1
                if !haskey(graph, src_ngram, dst_ngram)
                    graph[src_ngram, dst_ngram] = LightweightEdgeData(n - 1)
                end
                add_evidence!(graph[src_ngram, dst_ngram], dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

# ============================================================================
# Ultralight Graph Construction
#
# These builders create graphs with UltralightVertexData/UltralightEdgeData
# instead of full or lightweight types. No observation ID tracking — only
# topology + aggregated counts per dataset.
# ============================================================================

"""
Build ultralight k-mer graph that stores only topology + aggregated counts.

Core loop mirrors `_build_kmer_graph_core_lightweight` but creates
`UltralightKmerVertexData` vertices and `UltralightEdgeData` edges.
No observation ID tracking (unlike lightweight).
"""
function _build_kmer_graph_core_ultralight(
        records::Vector{FASTX.FASTA.Record},
        ::Val{K},
        ::Type{KmerType},
        dataset_id::String
) where {K, KmerType <: Kmers.Kmer}
    if KmerType <: Kmers.DNAKmer
        SeqType = BioSequences.LongDNA{4}
        KmerIterator = Kmers.UnambiguousDNAMers{K}
    elseif KmerType <: Kmers.RNAKmer
        SeqType = BioSequences.LongRNA{4}
        KmerIterator = Kmers.UnambiguousRNAMers{K}
    elseif KmerType <: Kmers.AAKmer
        SeqType = BioSequences.LongAA
        KmerIterator = try
            Kmers.UnambiguousAAMers{K}
        catch
            Kmers.FwAAMers{K}
        end
    else
        error("Unsupported k-mer type: $KmerType")
    end

    is_aa = KmerType <: Kmers.AAKmer
    actual_kmer_type = nothing
    for record in records
        test_seq = FASTX.sequence(SeqType, record)
        iter = KmerIterator(test_seq)
        first_item = iterate(iter)
        if first_item === nothing
            continue
        end
        value = first_item[1]
        test_kmer = is_aa ? value : value[1]
        actual_kmer_type = typeof(test_kmer)
        break
    end
    ActualKmerType = actual_kmer_type === nothing ? KmerType : actual_kmer_type

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = ActualKmerType,
        vertex_data_type = UltralightKmerVertexData{ActualKmerType},
        edge_data_type = UltralightEdgeData,
        weight_function = compute_edge_weight
    )

    for record in records
        observation_id = String(split(FASTX.identifier(record), ' ')[1])
        sequence = FASTX.sequence(SeqType, record)

        kmers_with_positions = if is_aa
            [(kmer, i) for (i, kmer) in enumerate(KmerIterator(sequence))]
        else
            collect(KmerIterator(sequence))
        end

        for (kmer, position) in kmers_with_positions
            if !haskey(graph, kmer)
                graph[kmer] = UltralightKmerVertexData(kmer)
            end
            add_evidence!(graph[kmer], dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        for i in 1:(length(kmers_with_positions) - 1)
            src_kmer, src_pos = kmers_with_positions[i]
            dst_kmer, dst_pos = kmers_with_positions[i + 1]

            if dst_pos == src_pos + 1
                if !haskey(graph, src_kmer, dst_kmer)
                    graph[src_kmer, dst_kmer] = UltralightEdgeData()
                end
                add_evidence!(graph[src_kmer, dst_kmer], dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

"""
Build ultralight strand-specific k-mer graph (outer wrapper with type detection).
"""
function build_kmer_graph_singlestrand_ultralight(
        records::Vector{<:FASTX.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    fasta_records = if records[1] isa FASTX.FASTQ.Record
        [FASTX.FASTA.Record(
             FASTX.identifier(r),
             FASTX.sequence(r)
         ) for r in records]
    else
        records
    end

    sequence_alphabet = _infer_sequence_alphabet(
        fasta_records;
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    if sequence_alphabet == :DNA
        return _build_kmer_graph_core_ultralight(
            fasta_records, Val(k), Kmers.DNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :RNA
        return _build_kmer_graph_core_ultralight(
            fasta_records, Val(k), Kmers.RNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :AA
        return _build_kmer_graph_core_ultralight(
            fasta_records, Val(k), Kmers.AAKmer{k}, dataset_id)
    else
        error("Unsupported sequence type for k-mer graph: $sequence_alphabet")
    end
end

"""
Build ultralight n-gram graph that stores only topology + aggregated counts.
No observation ID tracking.
"""
function build_ngram_graph_singlestrand_ultralight(
        strings::Vector{String},
        n::Int;
        dataset_id::String = "dataset_01"
)
    if isempty(strings)
        error("Cannot build graph from empty string vector")
    end

    if n < 1
        error("N-gram size must be positive, got: $n")
    end

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = UltralightStringVertexData,
        edge_data_type = UltralightEdgeData,
        weight_function = compute_edge_weight
    )

    for (string_idx, input_string) in enumerate(strings)
        observation_id = "string_$(lpad(string_idx, 6, '0'))"

        chars = collect(input_string)
        char_count = length(chars)

        if char_count < n
            continue
        end

        ngrams_with_positions = [(String(chars[i:(i + n - 1)]), i)
                                 for i in 1:(char_count - n + 1)]

        for (ngram, position) in ngrams_with_positions
            if !haskey(graph, ngram)
                graph[ngram] = UltralightStringVertexData(ngram)
            end
            add_evidence!(graph[ngram], dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        for i in 1:(length(ngrams_with_positions) - 1)
            src_ngram, src_pos = ngrams_with_positions[i]
            dst_ngram, dst_pos = ngrams_with_positions[i + 1]

            if dst_pos == src_pos + 1
                if !haskey(graph, src_ngram, dst_ngram)
                    graph[src_ngram, dst_ngram] = UltralightEdgeData(n - 1)
                end
                add_evidence!(graph[src_ngram, dst_ngram], dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

"""
    add_observations_to_ngram_graph_ultralight!(graph, strings, n; dataset_id)

Add observations from new strings to an existing ultralight n-gram graph.
"""
function add_observations_to_ngram_graph_ultralight!(
        graph::MetaGraphsNext.MetaGraph,
        strings::Vector{String},
        n::Int;
        dataset_id::String = "dataset_01"
)
    if isempty(strings)
        return graph
    end

    for (string_idx, input_string) in enumerate(strings)
        observation_id = "string_$(lpad(string_idx, 6, '0'))"

        if length(input_string) < n
            continue
        end

        ngrams_with_positions = [(input_string[i:(i + n - 1)], i)
                                 for i in 1:(length(input_string) - n + 1)]

        for (ngram, position) in ngrams_with_positions
            if !haskey(graph, ngram)
                graph[ngram] = UltralightStringVertexData(ngram)
            end
            add_evidence!(graph[ngram], dataset_id, observation_id,
                EvidenceEntry(position, Forward))
        end

        for i in 1:(length(ngrams_with_positions) - 1)
            src_ngram, src_pos = ngrams_with_positions[i]
            dst_ngram, dst_pos = ngrams_with_positions[i + 1]

            if dst_pos == src_pos + 1
                if !haskey(graph, src_ngram, dst_ngram)
                    graph[src_ngram, dst_ngram] = UltralightEdgeData(n - 1)
                end
                add_evidence!(graph[src_ngram, dst_ngram], dataset_id, observation_id,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

# ============================================================================
# Ultralight Quality Graph Construction
#
# These builders create graphs with quality-aware ultralight types.
# Quality scores are aggregated incrementally during construction.
# Requires FASTQ input (quality scores come from reads).
# ============================================================================

"""
Build ultralight quality k-mer graph from FASTQ records.

Aggregates quality scores incrementally during construction using Phred score
addition. Requires FASTQ input — will error on FASTA.
"""
function _build_kmer_graph_core_ultralight_quality(
        records::Vector{FASTX.FASTQ.Record},
        ::Val{K},
        ::Type{KmerType},
        dataset_id::String
) where {K, KmerType <: Kmers.Kmer}
    if KmerType <: Kmers.DNAKmer
        SeqType = BioSequences.LongDNA{4}
        KmerIterator = Kmers.UnambiguousDNAMers{K}
    elseif KmerType <: Kmers.RNAKmer
        SeqType = BioSequences.LongRNA{4}
        KmerIterator = Kmers.UnambiguousRNAMers{K}
    elseif KmerType <: Kmers.AAKmer
        SeqType = BioSequences.LongAA
        KmerIterator = try
            Kmers.UnambiguousAAMers{K}
        catch
            Kmers.FwAAMers{K}
        end
    else
        error("Unsupported k-mer type: $KmerType")
    end

    is_aa = KmerType <: Kmers.AAKmer
    actual_kmer_type = nothing
    for record in records
        test_seq = FASTX.sequence(SeqType, record)
        iter = KmerIterator(test_seq)
        first_item = iterate(iter)
        if first_item === nothing
            continue
        end
        value = first_item[1]
        test_kmer = is_aa ? value : value[1]
        actual_kmer_type = typeof(test_kmer)
        break
    end
    ActualKmerType = actual_kmer_type === nothing ? KmerType : actual_kmer_type

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = ActualKmerType,
        vertex_data_type = UltralightQualityKmerVertexData{ActualKmerType},
        edge_data_type = UltralightQualityEdgeData,
        weight_function = compute_edge_weight
    )

    for record in records
        observation_id = String(split(FASTX.identifier(record), ' ')[1])
        sequence = FASTX.sequence(SeqType, record)
        quality = Vector{UInt8}(FASTX.quality(record))

        kmers_with_positions = if is_aa
            [(kmer, i) for (i, kmer) in enumerate(KmerIterator(sequence))]
        else
            collect(KmerIterator(sequence))
        end

        for (kmer, position) in kmers_with_positions
            if !haskey(graph, kmer)
                graph[kmer] = UltralightQualityKmerVertexData(kmer)
            end

            kmer_quality = quality[position:(position + K - 1)]
            add_evidence!(graph[kmer], dataset_id, observation_id,
                QualityEvidenceEntry(position, Forward, kmer_quality))
        end

        for i in 1:(length(kmers_with_positions) - 1)
            src_kmer, src_pos = kmers_with_positions[i]
            dst_kmer, dst_pos = kmers_with_positions[i + 1]

            if dst_pos == src_pos + 1
                if !haskey(graph, src_kmer, dst_kmer)
                    graph[src_kmer, dst_kmer] = UltralightQualityEdgeData()
                end

                src_quality = quality[src_pos:(src_pos + K - 1)]
                dst_quality = quality[dst_pos:(dst_pos + K - 1)]
                add_evidence!(graph[src_kmer, dst_kmer], dataset_id, observation_id,
                    EdgeQualityEvidenceEntry(src_pos, dst_pos, Forward,
                        src_quality, dst_quality))
            end
        end
    end

    return graph
end

"""
Build ultralight quality strand-specific k-mer graph (outer wrapper).
Requires FASTQ input.
"""
function build_kmer_graph_singlestrand_ultralight_quality(
        records::Vector{FASTX.FASTQ.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    # Infer alphabet from FASTQ records (need FASTA for inference)
    fasta_records = [FASTX.FASTA.Record(
                         FASTX.identifier(r),
                         FASTX.sequence(r)
                     ) for r in records]

    sequence_alphabet = _infer_sequence_alphabet(
        fasta_records;
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    if sequence_alphabet == :DNA
        return _build_kmer_graph_core_ultralight_quality(
            records, Val(k), Kmers.DNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :RNA
        return _build_kmer_graph_core_ultralight_quality(
            records, Val(k), Kmers.RNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :AA
        return _build_kmer_graph_core_ultralight_quality(
            records, Val(k), Kmers.AAKmer{k}, dataset_id)
    else
        error("Unsupported sequence type for k-mer graph: $sequence_alphabet")
    end
end

# ============================================================================
# Lightweight Quality Graph Construction
#
# These builders create graphs with lightweight quality types.
# Combines observation ID tracking with quality aggregation.
# Requires FASTQ input.
# ============================================================================

"""
Build lightweight quality k-mer graph from FASTQ records.

Tracks observation IDs AND aggregates quality scores incrementally.
"""
function _build_kmer_graph_core_lightweight_quality(
        records::Vector{FASTX.FASTQ.Record},
        ::Val{K},
        ::Type{KmerType},
        dataset_id::String
) where {K, KmerType <: Kmers.Kmer}
    if KmerType <: Kmers.DNAKmer
        SeqType = BioSequences.LongDNA{4}
        KmerIterator = Kmers.UnambiguousDNAMers{K}
    elseif KmerType <: Kmers.RNAKmer
        SeqType = BioSequences.LongRNA{4}
        KmerIterator = Kmers.UnambiguousRNAMers{K}
    elseif KmerType <: Kmers.AAKmer
        SeqType = BioSequences.LongAA
        KmerIterator = try
            Kmers.UnambiguousAAMers{K}
        catch
            Kmers.FwAAMers{K}
        end
    else
        error("Unsupported k-mer type: $KmerType")
    end

    is_aa = KmerType <: Kmers.AAKmer
    actual_kmer_type = nothing
    for record in records
        test_seq = FASTX.sequence(SeqType, record)
        iter = KmerIterator(test_seq)
        first_item = iterate(iter)
        if first_item === nothing
            continue
        end
        value = first_item[1]
        test_kmer = is_aa ? value : value[1]
        actual_kmer_type = typeof(test_kmer)
        break
    end
    ActualKmerType = actual_kmer_type === nothing ? KmerType : actual_kmer_type

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = ActualKmerType,
        vertex_data_type = LightweightQualityKmerVertexData{ActualKmerType},
        edge_data_type = LightweightQualityEdgeData,
        weight_function = compute_edge_weight
    )

    for record in records
        observation_id = String(split(FASTX.identifier(record), ' ')[1])
        sequence = FASTX.sequence(SeqType, record)
        quality = Vector{UInt8}(FASTX.quality(record))

        kmers_with_positions = if is_aa
            [(kmer, i) for (i, kmer) in enumerate(KmerIterator(sequence))]
        else
            collect(KmerIterator(sequence))
        end

        for (kmer, position) in kmers_with_positions
            if !haskey(graph, kmer)
                graph[kmer] = LightweightQualityKmerVertexData(kmer)
            end

            kmer_quality = quality[position:(position + K - 1)]
            add_evidence!(graph[kmer], dataset_id, observation_id,
                QualityEvidenceEntry(position, Forward, kmer_quality))
        end

        for i in 1:(length(kmers_with_positions) - 1)
            src_kmer, src_pos = kmers_with_positions[i]
            dst_kmer, dst_pos = kmers_with_positions[i + 1]

            if dst_pos == src_pos + 1
                if !haskey(graph, src_kmer, dst_kmer)
                    graph[src_kmer, dst_kmer] = LightweightQualityEdgeData()
                end

                src_quality = quality[src_pos:(src_pos + K - 1)]
                dst_quality = quality[dst_pos:(dst_pos + K - 1)]
                add_evidence!(graph[src_kmer, dst_kmer], dataset_id, observation_id,
                    EdgeQualityEvidenceEntry(src_pos, dst_pos, Forward,
                        src_quality, dst_quality))
            end
        end
    end

    return graph
end

"""
Build lightweight quality strand-specific k-mer graph (outer wrapper).
Requires FASTQ input.
"""
function build_kmer_graph_singlestrand_lightweight_quality(
        records::Vector{FASTX.FASTQ.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    fasta_records = [FASTX.FASTA.Record(
                         FASTX.identifier(r),
                         FASTX.sequence(r)
                     ) for r in records]

    sequence_alphabet = _infer_sequence_alphabet(
        fasta_records;
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    if sequence_alphabet == :DNA
        return _build_kmer_graph_core_lightweight_quality(
            records, Val(k), Kmers.DNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :RNA
        return _build_kmer_graph_core_lightweight_quality(
            records, Val(k), Kmers.RNAKmer{k}, dataset_id)
    elseif sequence_alphabet == :AA
        return _build_kmer_graph_core_lightweight_quality(
            records, Val(k), Kmers.AAKmer{k}, dataset_id)
    else
        error("Unsupported sequence type for k-mer graph: $sequence_alphabet")
    end
end

# ============================================================================
# Variable-length String Graph Construction (OLC approach)
# ============================================================================

function _normalize_min_overlap(min_overlap::Int)
    if min_overlap < 1
        error("Minimum overlap must be positive, got: $min_overlap")
    end

    # Overlap detection only considers odd lengths.
    return iseven(min_overlap) ? min_overlap + 1 : min_overlap
end

"""
    build_string_graph_olc(strings; dataset_id="dataset_01", min_overlap=3)

Build variable-length string graph using Overlap-Layout-Consensus (OLC) approach.

Unlike de Bruijn graphs (fixed-length k-mers/n-grams), OLC graphs store complete
strings as vertices and connect them based on suffix-prefix overlaps.

**Overlap Requirement**: Only odd-length overlaps are allowed. This ensures
unambiguous overlap resolution and prevents issues with even-length palindromes.

# Arguments
- `strings::Vector{String}`: Input strings (complete sequences)
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `min_overlap::Int=3`: Minimum overlap length (odd-length overlaps only; even values are rounded up)
- `type_hint::Union{Nothing,Symbol}=nothing`: Optional alphabet hint (:DNA, :RNA, :AA)
- `ambiguous_action::Symbol=:dna`: Resolution for ambiguous alphabets (:dna, :rna, :aa, :error)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length string graph with overlap evidence

# Key Differences from N-gram Graphs
- Vertices represent ENTIRE strings, not fixed-length substrings
- Edges represent suffix-prefix overlaps between complete strings
- Variable overlap lengths stored in edge data (all odd-length)
- Suitable for assembly from reads, not pattern analysis

# Examples
```julia
# Assembly-like scenario
reads = ["ATCGATCG", "TCGATCGA", "GATCGATT"]
graph = build_string_graph_olc(reads; min_overlap=5)

# Find overlaps
for (src, dst) in MetaGraphsNext.edge_labels(graph)
    edge_data = graph[src, dst]
    println("Overlap: \$src → \$dst (\$(edge_data.overlap_length) chars)")
end
```
"""
function build_string_graph_olc(
        strings::Vector{String};
        dataset_id::String = "dataset_01",
        min_overlap::Int = 3
)
    if isempty(strings)
        error("Cannot build graph from empty string vector")
    end

    min_overlap = _normalize_min_overlap(min_overlap)

    # Create empty directed graph with String labels
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = StringVertexData,
        edge_data_type = StringEdgeData,
        weight_function = compute_edge_weight
    )

    # Add all strings as vertices
    for (string_idx, input_string) in enumerate(strings)
        observation_id = "string_$(lpad(string_idx, 6, '0'))"

        # Create vertex if it doesn't exist
        if !haskey(graph, input_string)
            vertex_data = StringVertexData(input_string)
            graph[input_string] = vertex_data
        end

        # Add evidence (position=1 since the vertex IS the full string)
        vertex_data = graph[input_string]
        add_evidence!(vertex_data, dataset_id, observation_id,
            EvidenceEntry(1, Forward))
    end

    # Find overlaps between all pairs of strings
    for i in 1:length(strings)
        for j in 1:length(strings)
            if i == j
                continue  # Skip self-loops
            end

            str1 = strings[i]
            str2 = strings[j]

            # Find overlap length where suffix of str1 matches prefix of str2
            overlap_len = find_overlap_length(str1, str2, min_overlap)

            if overlap_len >= min_overlap
                # Create edge if it doesn't exist
                if !haskey(graph, str1, str2)
                    edge_data = StringEdgeData(overlap_len)
                    graph[str1, str2] = edge_data
                end

                # Add evidence to edge
                edge_data = graph[str1, str2]
                observation_id_i = "string_$(lpad(i, 6, '0'))"
                observation_id_j = "string_$(lpad(j, 6, '0'))"

                # Position in str1 is where overlap starts (length - overlap + 1)
                # Position in str2 is always 1 (prefix)
                src_pos = length(str1) - overlap_len + 1
                dst_pos = 1

                add_evidence!(edge_data, dataset_id, observation_id_i,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

"""
    find_overlap_length(str1, str2, min_overlap)

Find the length of odd-length overlap between suffix of str1 and prefix of str2.

Only considers odd-length overlaps to avoid ambiguities with even-length palindromes.
Returns the longest odd-length overlap if >= min_overlap, otherwise 0.

# Examples
```julia
find_overlap_length("ATCGATCG", "TCGATCGA", 3)  # Returns 5 (CGATC, skips even length 6)
find_overlap_length("HELLO", "WORLD", 3)        # Returns 0 (no overlap)
```
"""
function find_overlap_length(str1::String, str2::String, min_overlap::Int)
    max_possible_overlap = min(length(str1), length(str2))

    # Try overlaps from longest to shortest, ONLY ODD LENGTHS
    for overlap_len in max_possible_overlap:-1:min_overlap
        # Skip even-length overlaps
        if iseven(overlap_len)
            continue
        end

        # Check if suffix of str1 matches prefix of str2
        suffix_start = length(str1) - overlap_len + 1
        suffix = str1[suffix_start:end]
        prefix = str2[1:overlap_len]

        if suffix == prefix
            return overlap_len
        end
    end

    return 0  # No sufficient odd-length overlap found
end

# ============================================================================
# Variable-length FASTA Graph Construction (OLC approach with BioSequences)
# ============================================================================

"""
    build_fasta_graph_olc(records; dataset_id="dataset_01", min_overlap=3)

Build variable-length FASTA graph using Overlap-Layout-Consensus (OLC) approach.

Uses complete BioSequences as vertices (NOT strings). This preserves the biological
sequence type and allows proper use of BioSequences operations.

**Overlap Requirement**: Only odd-length overlaps are allowed.

# Arguments
- `records::Vector{FASTX.FASTA.Record}`: Input FASTA records
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `min_overlap::Int=3`: Minimum overlap length (odd-length overlaps only; even values are rounded up)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length FASTA graph with BioSequence vertices

# Key Features
- Vertices are complete BioSequences (LongDNA, LongRNA, or LongAA)
- Edges represent suffix-prefix overlaps between sequences
- Odd-length overlaps only
- Evidence tracking with sequence IDs

# Examples
```julia
# Load FASTA records
records = collect(open_fastx("contigs.fasta"))

# Build OLC graph
graph = build_fasta_graph_olc(records; min_overlap=7)

# Find overlapping contigs
for (src, dst) in MetaGraphsNext.edge_labels(graph)
    edge_data = graph[src, dst]
    println("Overlap: \$(length(src))bp → \$(length(dst))bp (\$(edge_data.overlap_length)bp)")
end
```
"""
function build_fasta_graph_olc(
        records::Vector{FASTX.FASTA.Record};
        dataset_id::String = "dataset_01",
        min_overlap::Int = 3,
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    min_overlap = _normalize_min_overlap(min_overlap)

    sequence_alphabet = _infer_sequence_alphabet(
        records;
        type_hint = type_hint,
        ambiguous_action = ambiguous_action
    )

    for record in records[1:min(10, length(records))]
        seq_str = String(FASTX.sequence(record))
        if !parentmodule(Rhizomorph).validate_alphabet(seq_str, sequence_alphabet)
            error("Mixed sequence types detected. Cannot build graph from heterogeneous input.")
        end
    end

    # Dispatch to type-specific core function
    if sequence_alphabet == :DNA
        return _build_fasta_graph_olc_core(records, BioSequences.LongDNA{4}, dataset_id, min_overlap)
    elseif sequence_alphabet == :RNA
        return _build_fasta_graph_olc_core(records, BioSequences.LongRNA{4}, dataset_id, min_overlap)
    elseif sequence_alphabet == :AA
        return _build_fasta_graph_olc_core(records, BioSequences.LongAA, dataset_id, min_overlap)
    else
        error("Unsupported sequence type for FASTA graph: $sequence_alphabet")
    end
end

"""
Type-stable core function for FASTA OLC graph construction.

Uses BioSequences as vertex labels, not strings.
"""
function _build_fasta_graph_olc_core(
        records::Vector{FASTX.FASTA.Record},
        ::Type{SeqType},
        dataset_id::String,
        min_overlap::Int
) where {SeqType <: BioSequences.LongSequence}
    # Create empty directed graph with BioSequence labels
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = SeqType,
        vertex_data_type = BioSequenceVertexData{SeqType},
        edge_data_type = BioSequenceEdgeData,
        weight_function = compute_edge_weight
    )

    # Convert all sequences to proper BioSequence type
    sequences = [FASTX.sequence(SeqType, record) for record in records]
    identifiers = [String(split(FASTX.identifier(record), ' ')[1]) for record in records]

    # Add all sequences as vertices
    for (seq_idx, (sequence, observation_id)) in enumerate(zip(sequences, identifiers))
        # Create vertex if it doesn't exist
        if !haskey(graph, sequence)
            vertex_data = BioSequenceVertexData(sequence)
            graph[sequence] = vertex_data
        end

        # Add evidence (position=1 since the vertex IS the full sequence)
        vertex_data = graph[sequence]
        add_evidence!(vertex_data, dataset_id, observation_id,
            EvidenceEntry(1, Forward))
    end

    # Find overlaps between all pairs of sequences
    for i in 1:length(sequences)
        for j in 1:length(sequences)
            if i == j
                continue  # Skip self-loops
            end

            seq1 = sequences[i]
            seq2 = sequences[j]

            # Find overlap length where suffix of seq1 matches prefix of seq2
            overlap_len = find_biosequence_overlap_length(seq1, seq2, min_overlap)

            if overlap_len >= min_overlap
                # Create edge if it doesn't exist
                if !haskey(graph, seq1, seq2)
                    edge_data = BioSequenceEdgeData(overlap_len)
                    graph[seq1, seq2] = edge_data
                end

                # Add evidence to edge
                edge_data = graph[seq1, seq2]
                observation_id_i = identifiers[i]

                # Position in seq1 is where overlap starts (length - overlap + 1)
                # Position in seq2 is always 1 (prefix)
                src_pos = length(seq1) - overlap_len + 1
                dst_pos = 1

                add_evidence!(edge_data, dataset_id, observation_id_i,
                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
            end
        end
    end

    return graph
end

"""
    find_biosequence_overlap_length(seq1, seq2, min_overlap)

Find the length of odd-length overlap between suffix of seq1 and prefix of seq2.

Works with BioSequences (LongDNA, LongRNA, LongAA), not strings.
Only considers odd-length overlaps.

# Examples
```julia
seq1 = dna"ATCGATCG"
seq2 = dna"TCGATCGA"
find_biosequence_overlap_length(seq1, seq2, 3)  # Returns 5 (CGATC)
```
"""
function find_biosequence_overlap_length(
        seq1::BioSequences.LongSequence,
        seq2::BioSequences.LongSequence,
        min_overlap::Int
)
    max_possible_overlap = min(length(seq1), length(seq2))

    # Try overlaps from longest to shortest, ONLY ODD LENGTHS
    for overlap_len in max_possible_overlap:-1:min_overlap
        # Skip even-length overlaps
        if iseven(overlap_len)
            continue
        end

        # Check if suffix of seq1 matches prefix of seq2
        suffix_start = length(seq1) - overlap_len + 1
        suffix = seq1[suffix_start:end]
        prefix = seq2[1:overlap_len]

        if suffix == prefix
            return overlap_len
        end
    end

    return 0  # No sufficient odd-length overlap found
end

# ============================================================================
# Variable-length FASTQ Graph Construction (OLC approach with BioSequences and Quality)
# ============================================================================

"""
    build_fastq_graph_olc(records; dataset_id="dataset_01", min_overlap=3)

Build variable-length FASTQ graph using Overlap-Layout-Consensus (OLC) approach.

Uses complete BioSequences as vertices with quality information preserved (NOT strings).
This preserves the biological sequence type and per-base quality scores.

**Overlap Requirement**: Only odd-length overlaps are allowed.

# Arguments
- `records::Vector{FASTX.FASTQ.Record}`: Input FASTQ records
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `min_overlap::Int=3`: Minimum overlap length (odd-length overlaps only; even values are rounded up)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length FASTQ graph with BioSequence vertices and quality data

# Key Features
- Vertices are complete BioSequences (LongDNA, LongRNA, or LongAA)
- Quality scores preserved for each sequence
- Edges represent suffix-prefix overlaps between sequences
- Odd-length overlaps only
- Evidence tracking with sequence IDs and quality

# Examples
```julia
# Load FASTQ records
records = collect(open_fastx("reads.fastq.gz"))

# Build quality-aware OLC graph
graph = build_fastq_graph_olc(records; min_overlap=7)

# Find high-quality overlaps
for (src, dst) in MetaGraphsNext.edge_labels(graph)
    edge_data = graph[src, dst]
    println("Overlap: \$(length(src))bp → \$(length(dst))bp (\$(edge_data.overlap_length)bp)")
end
```
"""
function build_fastq_graph_olc(
        records::Vector{FASTX.FASTQ.Record};
        dataset_id::String = "dataset_01",
        min_overlap::Int = 3
)
    if isempty(records)
        throw(ArgumentError("Cannot build graph from empty record set"))
    end

    min_overlap = _normalize_min_overlap(min_overlap)

    # Detect sequence type from first record
    first_seq_str = String(FASTX.sequence(records[1]))
    biosequence = parentmodule(Rhizomorph).convert_sequence(first_seq_str)

    # FASTQ graphs support DNA, RNA, and amino acid sequences
    if !(biosequence isa
         Union{BioSequences.LongDNA, BioSequences.LongRNA, BioSequences.LongAA})
        error("FASTQ graphs only support DNA, RNA, or AA sequences, got: $(typeof(biosequence))")
    end

    # Validate all records are same type
    for record in records[2:min(10, length(records))]
        seq_str = String(FASTX.sequence(record))
        test_biosequence = parentmodule(Rhizomorph).convert_sequence(seq_str)
        if typeof(test_biosequence) != typeof(biosequence)
            error("Mixed sequence types detected. Cannot build graph from heterogeneous input.")
        end
    end

    # Dispatch to type-specific core function
    if biosequence isa BioSequences.LongDNA
        return _build_fastq_graph_olc_core(records, BioSequences.LongDNA{4}, dataset_id, min_overlap)
    elseif biosequence isa BioSequences.LongRNA
        return _build_fastq_graph_olc_core(records, BioSequences.LongRNA{4}, dataset_id, min_overlap)
    elseif biosequence isa BioSequences.LongAA
        return _build_fastq_graph_olc_core(records, BioSequences.LongAA, dataset_id, min_overlap)
    else
        error("Unsupported sequence type for FASTQ graph: $(typeof(biosequence))")
    end
end

"""
Type-stable core function for FASTQ OLC graph construction.

Uses BioSequences as vertex labels with quality information, not strings.
"""
function _build_fastq_graph_olc_core(
        records::Vector{FASTX.FASTQ.Record},
        ::Type{SeqType},
        dataset_id::String,
        min_overlap::Int
) where {SeqType <: BioSequences.LongSequence}
    # Create empty directed graph with BioSequence labels
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = SeqType,
        vertex_data_type = QualityBioSequenceVertexData{SeqType},
        edge_data_type = QualityBioSequenceEdgeData,
        weight_function = compute_edge_weight
    )

    # Convert all sequences to proper BioSequence type
    sequences = [FASTX.sequence(SeqType, record) for record in records]
    identifiers = [String(split(FASTX.identifier(record), ' ')[1]) for record in records]
    qualities = [Vector{UInt8}(FASTX.quality(record)) for record in records]

    # Add all sequences as vertices
    for (seq_idx, (sequence, observation_id, quality)) in
        enumerate(zip(sequences, identifiers, qualities))
        # Create vertex if it doesn't exist
        if !haskey(graph, sequence)
            vertex_data = QualityBioSequenceVertexData(sequence, quality)
            graph[sequence] = vertex_data
        end

        # Add evidence (position=1 since the vertex IS the full sequence)
        vertex_data = graph[sequence]
        if isempty(vertex_data.quality_scores)
            append!(vertex_data.quality_scores, quality)
        end
        add_evidence!(vertex_data, dataset_id, observation_id,
            QualityEvidenceEntry(1, Forward, quality))
    end

    # Find overlaps between all pairs of sequences
    for i in 1:length(sequences)
        for j in 1:length(sequences)
            if i == j
                continue  # Skip self-loops
            end

            seq1 = sequences[i]
            seq2 = sequences[j]

            # Find overlap length where suffix of seq1 matches prefix of seq2
            overlap_len = find_biosequence_overlap_length(seq1, seq2, min_overlap)

            if overlap_len >= min_overlap
                # Create edge if it doesn't exist
                if !haskey(graph, seq1, seq2)
                    edge_data = QualityBioSequenceEdgeData(overlap_len)
                    graph[seq1, seq2] = edge_data
                end

                # Add evidence to edge with quality scores
                edge_data = graph[seq1, seq2]
                observation_id_i = identifiers[i]

                # Position in seq1 is where overlap starts (length - overlap + 1)
                # Position in seq2 is always 1 (prefix)
                src_pos = length(seq1) - overlap_len + 1
                dst_pos = 1

                # Extract quality scores for overlapping regions
                src_quality = qualities[i][src_pos:end]  # Suffix quality
                dst_quality = qualities[j][1:overlap_len]  # Prefix quality

                add_evidence!(edge_data, dataset_id, observation_id_i,
                    EdgeQualityEvidenceEntry(src_pos, dst_pos, Forward,
                        src_quality, dst_quality))
            end
        end
    end

    return graph
end

# ============================================================================
# Reduced Memory Profile Conversion Helpers
# ============================================================================

"""
    _detect_reduced_memory_profile(vertex_data)

Detect memory profile symbol from a reduced vertex data type.
Returns a Symbol for dispatch to `_create_reduced_vertex`, or `nothing` for full types.
"""
function _detect_reduced_memory_profile(vertex_data)
    if vertex_data isa LightweightKmerVertexData
        return :lightweight
    elseif vertex_data isa UltralightKmerVertexData
        return :ultralight
    elseif vertex_data isa UltralightQualityKmerVertexData
        return :ultralight_quality
    elseif vertex_data isa LightweightQualityKmerVertexData
        return :lightweight_quality
    else
        return nothing
    end
end

"""
    _detect_reduced_memory_profile_edge(edge_data)

Detect memory profile symbol from a reduced edge data type.
"""
function _detect_reduced_memory_profile_edge(edge_data)
    if edge_data isa LightweightEdgeData
        return :lightweight
    elseif edge_data isa UltralightEdgeData
        return :ultralight
    elseif edge_data isa UltralightQualityEdgeData
        return :ultralight_quality
    elseif edge_data isa LightweightQualityEdgeData
        return :lightweight_quality
    else
        return nothing
    end
end

"""
    _merge_reduced_vertex_into!(dest, src; reverse_quality=false)

Merge counts/observations/quality from `src` into `dest`.
When `reverse_quality=true`, reverse quality vectors before merging (for RC strand).
"""
function _merge_reduced_vertex_into!(dest, src; reverse_quality::Bool = false)
    dest.total_count += src.total_count
    for (ds_id, ds_count) in src.dataset_counts
        dest.dataset_counts[ds_id] = get(dest.dataset_counts, ds_id, 0) + ds_count
    end
    if hasfield(typeof(dest), :dataset_observations)
        for (ds_id, obs_set) in src.dataset_observations
            if !haskey(dest.dataset_observations, ds_id)
                dest.dataset_observations[ds_id] = Set{String}()
            end
            union!(dest.dataset_observations[ds_id], obs_set)
        end
    end
    if hasfield(typeof(dest), :joint_quality)
        src_quality = reverse_quality ? reverse(src.joint_quality) : src.joint_quality
        _merge_quality_vector!(dest.joint_quality, src_quality)
        if hasfield(typeof(dest), :dataset_joint_quality)
            for (ds_id, ds_qual) in src.dataset_joint_quality
                src_ds_qual = reverse_quality ? reverse(ds_qual) : ds_qual
                existing = get(dest.dataset_joint_quality, ds_id, nothing)
                if isnothing(existing)
                    dest.dataset_joint_quality[ds_id] = collect(src_ds_qual)
                else
                    _merge_quality_vector!(existing, collect(src_ds_qual))
                end
            end
        end
    end
end

"""
    _merge_reduced_edge_into!(dest, src; swap_quality=false)

Merge counts/observations/quality from `src` edge into `dest` edge.
When `swap_quality=true`, swap and reverse from/to quality vectors (for RC edges).
"""
function _merge_reduced_edge_into!(dest, src; swap_quality::Bool = false)
    dest.total_count += src.total_count
    for (ds_id, ds_count) in src.dataset_counts
        dest.dataset_counts[ds_id] = get(dest.dataset_counts, ds_id, 0) + ds_count
    end
    if hasfield(typeof(dest), :dataset_observations)
        for (ds_id, obs_set) in src.dataset_observations
            if !haskey(dest.dataset_observations, ds_id)
                dest.dataset_observations[ds_id] = Set{String}()
            end
            union!(dest.dataset_observations[ds_id], obs_set)
        end
    end
    if hasfield(typeof(dest), :from_joint_quality)
        if swap_quality
            _merge_quality_vector!(dest.from_joint_quality, reverse(src.to_joint_quality))
            _merge_quality_vector!(dest.to_joint_quality, reverse(src.from_joint_quality))
        else
            _merge_quality_vector!(dest.from_joint_quality, src.from_joint_quality)
            _merge_quality_vector!(dest.to_joint_quality, src.to_joint_quality)
        end
    end
end

"""
    _copy_reduced_vertex_data(src, new_kmer; reverse_quality=false)

Create a new reduced vertex with the same type as `src` but with `new_kmer`,
copying all count/observation/quality data. When `reverse_quality=true`,
reverse quality vectors (for RC conversion).
"""
function _copy_reduced_vertex_data(src, new_kmer; reverse_quality::Bool = false)
    profile = _detect_reduced_memory_profile(src)
    dest = _create_reduced_vertex(new_kmer, profile)
    dest.total_count = src.total_count
    merge!(dest.dataset_counts, src.dataset_counts)
    if hasfield(typeof(dest), :dataset_observations)
        for (ds_id, obs_set) in src.dataset_observations
            dest.dataset_observations[ds_id] = copy(obs_set)
        end
    end
    if hasfield(typeof(dest), :joint_quality)
        q = reverse_quality ? reverse(src.joint_quality) : copy(src.joint_quality)
        append!(dest.joint_quality, q)
        if hasfield(typeof(dest), :dataset_joint_quality)
            for (ds_id, ds_qual) in src.dataset_joint_quality
                dest.dataset_joint_quality[ds_id] = reverse_quality ?
                                                    reverse(ds_qual) : copy(ds_qual)
            end
        end
    end
    return dest
end

"""
    _copy_reduced_edge_data(src; swap_quality=false)

Copy edge data with optional quality swap (for RC edges where src/dst flip).
"""
function _copy_reduced_edge_data(src; swap_quality::Bool = false)
    profile = _detect_reduced_memory_profile_edge(src)
    dest = _create_reduced_edge(src.overlap_length, profile)
    dest.total_count = src.total_count
    merge!(dest.dataset_counts, src.dataset_counts)
    if hasfield(typeof(dest), :dataset_observations)
        for (ds_id, obs_set) in src.dataset_observations
            dest.dataset_observations[ds_id] = copy(obs_set)
        end
    end
    if hasfield(typeof(dest), :from_joint_quality)
        if swap_quality
            append!(dest.from_joint_quality, reverse(src.to_joint_quality))
            append!(dest.to_joint_quality, reverse(src.from_joint_quality))
        else
            append!(dest.from_joint_quality, copy(src.from_joint_quality))
            append!(dest.to_joint_quality, copy(src.to_joint_quality))
        end
    end
    return dest
end

# ============================================================================
# Reduced Memory Profile Doublestrand/Canonical Conversion
# ============================================================================

"""
    _convert_to_doublestrand_reduced(singlestrand_graph)

Convert a singlestrand reduced-type graph to doublestrand representation.
Purely topological: creates RC vertices with copied counts, no evidence iteration.
"""
function _convert_to_doublestrand_reduced(singlestrand_graph)
    first_label = first(MetaGraphsNext.labels(singlestrand_graph))
    first_vertex = singlestrand_graph[first_label]
    KmerType = typeof(first_label)
    VertexDataType = typeof(first_vertex)
    profile = _detect_reduced_memory_profile(first_vertex)

    # Determine edge data type
    EdgeDataType = if !isempty(MetaGraphsNext.edge_labels(singlestrand_graph))
        src, dst = first(MetaGraphsNext.edge_labels(singlestrand_graph))
        typeof(singlestrand_graph[src, dst])
    else
        if profile == :ultralight
            UltralightEdgeData
        elseif profile == :lightweight
            LightweightEdgeData
        elseif profile == :ultralight_quality
            UltralightQualityEdgeData
        elseif profile == :lightweight_quality
            LightweightQualityEdgeData
        end
    end

    doublestrand_graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = KmerType,
        vertex_data_type = VertexDataType,
        edge_data_type = EdgeDataType,
        weight_function = compute_edge_weight
    )

    # Copy all forward vertices (deep copy to avoid aliasing mutable structs)
    for kmer in MetaGraphsNext.labels(singlestrand_graph)
        doublestrand_graph[kmer] = _copy_reduced_vertex_data(
            singlestrand_graph[kmer], kmer
        )
    end

    # Add reverse complement vertices
    for kmer in MetaGraphsNext.labels(singlestrand_graph)
        rc_kmer = BioSequences.reverse_complement(kmer)

        # Skip palindromes
        if rc_kmer == kmer
            continue
        end

        vertex_data = singlestrand_graph[kmer]

        if haskey(doublestrand_graph, rc_kmer)
            # RC already exists (was observed forward) — merge data
            _merge_reduced_vertex_into!(
                doublestrand_graph[rc_kmer], vertex_data;
                reverse_quality = true
            )
        else
            # Create new RC vertex with copied data
            rc_vertex = _copy_reduced_vertex_data(
                vertex_data, rc_kmer;
                reverse_quality = true
            )
            doublestrand_graph[rc_kmer] = rc_vertex
        end
    end

    # Copy all forward edges (deep copy to avoid aliasing mutable structs)
    for (src_kmer, dst_kmer) in MetaGraphsNext.edge_labels(singlestrand_graph)
        doublestrand_graph[src_kmer, dst_kmer] = _copy_reduced_edge_data(
            singlestrand_graph[src_kmer, dst_kmer]
        )
    end

    # Add reverse complement edges (reversed direction: RC(dst) → RC(src))
    for (src_kmer, dst_kmer) in MetaGraphsNext.edge_labels(singlestrand_graph)
        rc_src = BioSequences.reverse_complement(dst_kmer)
        rc_dst = BioSequences.reverse_complement(src_kmer)

        edge_data = singlestrand_graph[src_kmer, dst_kmer]

        if haskey(doublestrand_graph, rc_src, rc_dst)
            _merge_reduced_edge_into!(
                doublestrand_graph[rc_src, rc_dst], edge_data;
                swap_quality = true
            )
        else
            rc_edge = _copy_reduced_edge_data(edge_data; swap_quality = true)
            doublestrand_graph[rc_src, rc_dst] = rc_edge
        end
    end

    return doublestrand_graph
end

"""
    _convert_to_canonical_reduced(singlestrand_graph)

Convert a singlestrand reduced-type graph to canonical representation.
Merges forward and RC k-mers into canonical forms with combined counts.
Creates an UNDIRECTED graph.
"""
function _convert_to_canonical_reduced(singlestrand_graph)
    first_label = first(MetaGraphsNext.labels(singlestrand_graph))
    first_vertex = singlestrand_graph[first_label]
    KmerType = typeof(first_label)
    VertexDataType = typeof(first_vertex)
    profile = _detect_reduced_memory_profile(first_vertex)

    # Determine edge data type
    EdgeDataType = if !isempty(MetaGraphsNext.edge_labels(singlestrand_graph))
        src, dst = first(MetaGraphsNext.edge_labels(singlestrand_graph))
        typeof(singlestrand_graph[src, dst])
    else
        if profile == :ultralight
            UltralightEdgeData
        elseif profile == :lightweight
            LightweightEdgeData
        elseif profile == :ultralight_quality
            UltralightQualityEdgeData
        elseif profile == :lightweight_quality
            LightweightQualityEdgeData
        end
    end

    # Create undirected canonical graph
    canonical_graph = MetaGraphsNext.MetaGraph(
        Graphs.Graph();
        label_type = KmerType,
        vertex_data_type = VertexDataType,
        edge_data_type = EdgeDataType,
        weight_function = compute_edge_weight
    )

    # Track processed kmers
    processed = Set{KmerType}()

    # Process all vertices
    for kmer in MetaGraphsNext.labels(singlestrand_graph)
        if kmer in processed
            continue
        end

        canon_kmer = BioSequences.canonical(kmer)
        rc_kmer = BioSequences.reverse_complement(kmer)

        push!(processed, kmer)
        push!(processed, rc_kmer)

        # Create canonical vertex if it doesn't exist
        if !haskey(canonical_graph, canon_kmer)
            canonical_graph[canon_kmer] = _create_reduced_vertex(canon_kmer, profile)
        end

        # Merge forward data (reverse quality if kmer != canon_kmer)
        fwd_data = singlestrand_graph[kmer]
        _merge_reduced_vertex_into!(
            canonical_graph[canon_kmer], fwd_data;
            reverse_quality = (kmer != canon_kmer)
        )

        # Merge RC data if it exists in source graph
        if rc_kmer != kmer && haskey(singlestrand_graph, rc_kmer)
            rc_data = singlestrand_graph[rc_kmer]
            _merge_reduced_vertex_into!(
                canonical_graph[canon_kmer], rc_data;
                reverse_quality = (rc_kmer != canon_kmer)
            )
        end
    end

    # Process edges
    for (src_kmer, dst_kmer) in MetaGraphsNext.edge_labels(singlestrand_graph)
        canon_src = BioSequences.canonical(src_kmer)
        canon_dst = BioSequences.canonical(dst_kmer)

        edge_data = singlestrand_graph[src_kmer, dst_kmer]

        if !haskey(canonical_graph, canon_src, canon_dst)
            canon_edge = _copy_reduced_edge_data(edge_data)
            canonical_graph[canon_src, canon_dst] = canon_edge
        else
            _merge_reduced_edge_into!(
                canonical_graph[canon_src, canon_dst], edge_data
            )
        end
    end

    return canonical_graph
end
