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

# Returns
- `MetaGraphsNext.MetaGraph`: Strand-specific k-mer graph with evidence
"""
function build_kmer_graph_singlestrand(
    records::Vector{<:FASTX.Record},
    k::Int;
    dataset_id::String="dataset_01"
)
    if isempty(records)
        error("Cannot build graph from empty record set")
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

    # Detect sequence type from first record
    first_seq = FASTX.sequence(fasta_records[1])
    biosequence = parentmodule(Rhizomorph).convert_sequence(String(first_seq))

    # Validate all records are same type
    for record in fasta_records[2:min(10, length(fasta_records))]
        seq = FASTX.sequence(record)
        test_biosequence = parentmodule(Rhizomorph).convert_sequence(String(seq))
        if typeof(test_biosequence) != typeof(biosequence)
            error("Mixed sequence types detected. Cannot build graph from heterogeneous input.")
        end
    end

    # Dispatch to type-stable core function
    if biosequence isa BioSequences.LongDNA
        return _build_kmer_graph_core(fasta_records, Val(k), Kmers.DNAKmer{k}, dataset_id)
    elseif biosequence isa BioSequences.LongRNA
        return _build_kmer_graph_core(fasta_records, Val(k), Kmers.RNAKmer{k}, dataset_id)
    elseif biosequence isa BioSequences.LongAA
        return _build_kmer_graph_core(fasta_records, Val(k), Kmers.AAKmer{k}, dataset_id)
    else
        error("Unsupported sequence type for k-mer graph: $(typeof(biosequence))")
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

    # Get actual kmer type from iterator (includes all type parameters)
    # Need to extract from first sequence to get the concrete type
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
        label_type=ActualKmerType,
        vertex_data_type=KmerVertexData{ActualKmerType},
        edge_data_type=KmerEdgeData
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

# Returns
- `MetaGraphsNext.MetaGraph`: Strand-specific qualmer graph with quality-aware evidence
"""
function build_qualmer_graph_singlestrand(
    records::Vector{FASTX.FASTQ.Record},
    k::Int;
    dataset_id::String="dataset_01"
)
    if isempty(records)
        error("Cannot build graph from empty record set")
    end

    # Detect sequence type from first record
    first_seq = FASTX.sequence(records[1])
    biosequence = parentmodule(Rhizomorph).convert_sequence(String(first_seq))

    # Validate all records are same type
    for record in records[2:min(10, length(records))]
        seq = FASTX.sequence(record)
        test_biosequence = parentmodule(Rhizomorph).convert_sequence(String(seq))
        if typeof(test_biosequence) != typeof(biosequence)
            error("Mixed sequence types detected. Cannot build graph from heterogeneous input.")
        end
    end

    # Dispatch to type-stable core function
    if biosequence isa BioSequences.LongDNA
        return _build_qualmer_graph_core(records, Val(k), Kmers.DNAKmer{k}, dataset_id)
    elseif biosequence isa BioSequences.LongRNA
        return _build_qualmer_graph_core(records, Val(k), Kmers.RNAKmer{k}, dataset_id)
    elseif biosequence isa BioSequences.LongAA
        return _build_qualmer_graph_core(records, Val(k), Kmers.AAKmer{k}, dataset_id)
    else
        error("Unsupported sequence type for qualmer graph: $(typeof(biosequence))")
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
        label_type=ActualKmerType,
        vertex_data_type=QualmerVertexData{ActualKmerType},
        edge_data_type=QualmerEdgeData
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
    dataset_id::String="dataset_01"
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
function get_observation_id_from_index(index::Int, prefix::String="obs")
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

# ============================================================================
# Doublestrand Graph Construction (Canonical K-mers)
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
graph = build_kmer_graph_doublestrand(records, 3)

# Both ATG and its RC (CAT) map to the canonical form
```
"""
function build_kmer_graph_doublestrand(
    records::Vector{<:FASTX.Record},
    k::Int;
    dataset_id::String="dataset_01"
)
    if isempty(records)
        error("Cannot build graph from empty record set")
    end

    # First build singlestrand graph
    singlestrand_graph = build_kmer_graph_singlestrand(records, k; dataset_id=dataset_id)

    # Convert to canonical (doublestrand) representation
    return convert_to_doublestrand(singlestrand_graph)
end

"""
Build canonical (doublestrand) qualmer graph from FASTQ records.

Similar to `build_kmer_graph_doublestrand` but preserves quality information.

# Arguments
- `records::Vector{FASTX.FASTQ.Record}`: Input FASTQ records
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier

# Returns
- `MetaGraphsNext.MetaGraph`: Canonical qualmer graph with merged evidence
"""
function build_qualmer_graph_doublestrand(
    records::Vector{FASTX.FASTQ.Record},
    k::Int;
    dataset_id::String="dataset_01"
)
    if isempty(records)
        error("Cannot build graph from empty record set")
    end

    # First build singlestrand qualmer graph
    singlestrand_graph = build_qualmer_graph_singlestrand(records, k; dataset_id=dataset_id)

    # Convert to canonical (doublestrand) representation
    return convert_to_doublestrand(singlestrand_graph)
end

"""
    convert_to_doublestrand(singlestrand_graph)

Convert a singlestrand graph to doublestrand (canonical) representation.

Merges forward and reverse complement k-mers into their canonical forms.
Evidence from both strands is combined.

# Arguments
- `singlestrand_graph`: Graph with strand-specific k-mers

# Returns
- New graph with canonical k-mers and merged evidence
"""
function convert_to_doublestrand(singlestrand_graph)
    # Get graph type information from first vertex
    if isempty(MetaGraphsNext.labels(singlestrand_graph))
        return singlestrand_graph  # Empty graph stays empty
    end

    first_label = first(MetaGraphsNext.labels(singlestrand_graph))
    first_vertex = singlestrand_graph[first_label]
    KmerType = typeof(first_label)
    VertexDataType = typeof(first_vertex)

    # Amino acids don't have reverse complement
    if KmerType <: Kmers.AAKmer
        error("Cannot create doublestrand graph for amino acid sequences (no reverse complement)")
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
        Graphs.DiGraph();
        label_type=KmerType,
        vertex_data_type=VertexDataType,
        edge_data_type=EdgeDataType
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

    VertexDataType = typeof(first_vertex)
    EdgeDataType = if VertexDataType <: KmerVertexData
        KmerEdgeData
    else
        QualmerEdgeData
    end

    # Create new singlestrand graph
    singlestrand_graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type=typeof(first_kmer),
        vertex_data_type=VertexDataType,
        edge_data_type=EdgeDataType,
        graph_data="Singlestrand graph (unfolded from canonical)"
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
                        add_evidence!(singlestrand_graph[canon_src, canon_dst], dataset_id, obs_id, evidence)
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
                        add_evidence!(singlestrand_graph[rc_dst, rc_src], dataset_id, obs_id, flipped_evidence)
                    end
                end
            end
        end
    end

    return singlestrand_graph
end
