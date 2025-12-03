# Strand conversions for variable-length OLC graphs (FASTA/FASTQ)

import MetaGraphsNext
import Graphs
import BioSequences

"""
    convert_variable_length_to_doublestrand(graph::MetaGraphsNext.MetaGraph)

Convert a singlestrand variable-length graph (DNA/RNA) to a doublestrand directed graph.

Creates reverse-complement vertices and edges; evidence is preserved as-is (assumed strand-aware).

- Supported: `BioSequences.LongDNA{4}` and `BioSequences.LongRNA{4}` labels
- Not supported: amino acids or generic strings (throws)
"""
function convert_variable_length_to_doublestrand(graph::MetaGraphsNext.MetaGraph)
    if isempty(MetaGraphsNext.labels(graph))
        return graph
    end
    first_label = first(MetaGraphsNext.labels(graph))
    if !(first_label isa BioSequences.LongDNA{4} || first_label isa BioSequences.LongRNA{4})
        error("Doublestrand conversion only supported for DNA/RNA variable-length graphs")
    end

    VertexType = typeof(first_label)
    VertexDataType = typeof(graph[first_label])
    EdgeDataType = if !isempty(MetaGraphsNext.edge_labels(graph))
        src, dst = first(MetaGraphsNext.edge_labels(graph))
        typeof(graph[src, dst])
    else
        Any
    end

    ds_graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type=VertexType,
        vertex_data_type=VertexDataType,
        edge_data_type=EdgeDataType
    )

    # Add forward vertices
    for v in MetaGraphsNext.labels(graph)
        ds_graph[v] = graph[v]
    end

    # Add RC vertices
    for v in MetaGraphsNext.labels(graph)
        rc_v = BioSequences.reverse_complement(v)
        if rc_v == v
            continue
        end
        ds_graph[rc_v] = graph[v]  # reuse data; evidence remains strand-aware
    end

    # Add edges
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        ds_graph[src, dst] = graph[src, dst]
    end

    # Add RC edges
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        rc_src = BioSequences.reverse_complement(dst)
        rc_dst = BioSequences.reverse_complement(src)
        ds_graph[rc_src, rc_dst] = graph[src, dst]
    end

    return ds_graph
end

"""
    convert_variable_length_to_canonical(graph::MetaGraphsNext.MetaGraph)

Convert a singlestrand variable-length graph (DNA/RNA) to canonical form.

Canonical vertices use the `BioSequences.canonical` label; evidence from forward/RC is merged
with strand flags preserved via `flip_evidence_strand`.

- Supported: `BioSequences.LongDNA{4}` and `BioSequences.LongRNA{4}` labels
- Not supported: amino acids or generic strings (throws)
"""
function convert_variable_length_to_canonical(graph::MetaGraphsNext.MetaGraph)
    if isempty(MetaGraphsNext.labels(graph))
        return graph
    end
    first_label = first(MetaGraphsNext.labels(graph))
    if !(first_label isa BioSequences.LongDNA{4} || first_label isa BioSequences.LongRNA{4})
        error("Canonical conversion only supported for DNA/RNA variable-length graphs")
    end

    VertexType = typeof(first_label)
    VertexDataType = typeof(graph[first_label])
    EdgeDataType = if !isempty(MetaGraphsNext.edge_labels(graph))
        src, dst = first(MetaGraphsNext.edge_labels(graph))
        typeof(graph[src, dst])
    else
        Any
    end

    canon_graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type=VertexType,
        vertex_data_type=VertexDataType,
        edge_data_type=EdgeDataType
    )

    processed = Set{VertexType}()

    for v in MetaGraphsNext.labels(graph)
        if v in processed
            continue
        end
        rc_v = BioSequences.reverse_complement(v)
        canon_v = BioSequences.canonical(v)
        push!(processed, v)
        push!(processed, rc_v)
        # initialize vertex
        if !haskey(canon_graph, canon_v)
            canon_graph[canon_v] = graph[v]
        end
        # merge evidence from v
        vdata = graph[v]
        if hasproperty(vdata, :evidence)
            for (dsid, obs) in vdata.evidence
                for (obsid, evidset) in obs
                    for ev in evidset
                        Mycelia.Rhizomorph.add_evidence!(canon_graph[canon_v], dsid, obsid, v == canon_v ? ev : Mycelia.Rhizomorph.flip_evidence_strand(ev))
                    end
                end
            end
        end
        # merge evidence from rc_v if present
        if rc_v != v && haskey(graph, rc_v)
            rcdata = graph[rc_v]
            if hasproperty(rcdata, :evidence)
                for (dsid, obs) in rcdata.evidence
                    for (obsid, evidset) in obs
                        for ev in evidset
                            Mycelia.Rhizomorph.add_evidence!(canon_graph[canon_v], dsid, obsid, rc_v == canon_v ? ev : Mycelia.Rhizomorph.flip_evidence_strand(ev))
                        end
                    end
                end
            end
        end
    end

    # edges
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        canon_src = BioSequences.canonical(src)
        canon_dst = BioSequences.canonical(dst)
        if !haskey(canon_graph, canon_src, canon_dst)
            canon_graph[canon_src, canon_dst] = graph[src, dst]
        end
        edata = graph[src, dst]
        if hasproperty(edata, :evidence)
            for (dsid, obs) in edata.evidence
                for (obsid, evidset) in obs
                    for ev in evidset
                        Mycelia.Rhizomorph.add_evidence!(canon_graph[canon_src, canon_dst], dsid, obsid, (src == canon_src && dst == canon_dst) ? ev : Mycelia.Rhizomorph.flip_evidence_strand(ev))
                    end
                end
            end
        end
    end

    return canon_graph
end
