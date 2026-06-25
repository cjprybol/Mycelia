"""
$(DocStringExtensions.TYPEDSIGNATURES)

Default emission model for the graph-agnostic Viterbi correction interface.

The callback shape is intentionally broader than B1 needs so follow-on beads can
swap in alphabet- and quality-aware models without changing the dynamic
programming seam.
"""
function default_viterbi_emission_logp(
        observed_unit,
        node,
        alphabet;
        quality = nothing,
        error_rate::Float64 = 0.01
)::Float64
    if error_rate <= 0.0 || error_rate >= 1.0
        throw(ArgumentError("error_rate must be in (0, 1), got $error_rate"))
    end

    observed = string(observed_unit)
    expected = string(node)
    if observed == expected
        return log1p(-error_rate)
    end

    edit_count = _hamming_like_distance(observed, expected)
    return edit_count == 0 ? log1p(-error_rate) : edit_count * log(error_rate)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Configuration for `correct_observations`.

`emission_logp` is the extension seam for B2-B4: it receives
`(observed_unit, node, alphabet; quality)` and returns a log-probability. B1
keeps the default simple so the legacy B0 oracle remains the behavioral anchor.
"""
struct ViterbiCorrectionConfig{F <: Function}
    error_rate::Float64
    verbosity::String
    emission_logp::F
    alphabet::Symbol
    max_steps::Union{Nothing, Int}
    target_vertex::Any
    start_strand::Rhizomorph.StrandOrientation
    edge_weight::Function

    function ViterbiCorrectionConfig{F}(;
            error_rate::Float64 = 0.01,
            verbosity::String = "dataset",
            emission_logp::F,
            alphabet::Symbol = :auto,
            max_steps::Union{Nothing, Int} = nothing,
            target_vertex = nothing,
            start_strand::Rhizomorph.StrandOrientation = Rhizomorph.Forward,
            edge_weight::Function = Rhizomorph.edge_data_weight
    ) where {F <: Function}
        if error_rate <= 0.0 || error_rate >= 0.5
            throw(ArgumentError("error_rate must be in (0, 0.5), got $error_rate"))
        end
        if !(verbosity in ("debug", "reads", "dataset", "silent"))
            throw(ArgumentError("unsupported verbosity: $verbosity"))
        end
        if max_steps !== nothing && max_steps < 0
            throw(ArgumentError("max_steps must be non-negative, got $max_steps"))
        end
        return new{F}(
            error_rate,
            verbosity,
            emission_logp,
            alphabet,
            max_steps,
            target_vertex,
            start_strand,
            edge_weight
        )
    end
end

function ViterbiCorrectionConfig(;
        error_rate::Float64 = 0.01,
        verbosity::String = "dataset",
        emission_logp::F = default_viterbi_emission_logp,
        alphabet::Symbol = :auto,
        max_steps::Union{Nothing, Int} = nothing,
        target_vertex = nothing,
        start_strand::Rhizomorph.StrandOrientation = Rhizomorph.Forward,
        edge_weight::Function = Rhizomorph.edge_data_weight
)::ViterbiCorrectionConfig{F} where {F <: Function}
    return ViterbiCorrectionConfig{F}(
        error_rate = error_rate,
        verbosity = verbosity,
        emission_logp = emission_logp,
        alphabet = alphabet,
        max_steps = max_steps,
        target_vertex = target_vertex,
        start_strand = start_strand,
        edge_weight = edge_weight
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A dynamic-programming state for graph-agnostic Viterbi correction.
"""
struct ViterbiCorrectionState{V}
    vertex_label::V
    strand::Rhizomorph.StrandOrientation
    logp::Float64
    depth::Int
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Result object returned by the B1 `correct_observations` interface.
"""
struct ViterbiCorrectionResult{C, P}
    corrected_observations::C
    paths::P
    diagnostics::Dict{Symbol, Any}
end

const _VITERBI_CORRECTION_VERTEX_DATA = Union{
    Rhizomorph.KmerVertexData,
    Rhizomorph.QualmerVertexData,
    Rhizomorph.BioSequenceVertexData,
    Rhizomorph.QualityBioSequenceVertexData,
    Rhizomorph.StringVertexData,
    Rhizomorph.QualityStringVertexData,
    Rhizomorph.LightweightKmerVertexData,
    Rhizomorph.LightweightBioSequenceVertexData,
    Rhizomorph.LightweightStringVertexData,
    Rhizomorph.UltralightKmerVertexData,
    Rhizomorph.UltralightBioSequenceVertexData,
    Rhizomorph.UltralightStringVertexData,
    Rhizomorph.UltralightQualityKmerVertexData,
    Rhizomorph.UltralightQualityBioSequenceVertexData,
    Rhizomorph.LightweightQualityKmerVertexData,
    Rhizomorph.LightweightQualityBioSequenceVertexData
}

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Correct observations against a graph with a swappable emission model.

This is the public B1 seam. Legacy stranded-kmer graphs delegate to the B0
oracle-preserving implementation; `MetaGraphsNext.MetaGraph` inputs route
through the Rhizomorph weighted-graph and `viterbi_decode_next` primitives.
"""
function correct_observations(
        graph,
        observations = _default_correction_observations(graph);
        config::ViterbiCorrectionConfig = ViterbiCorrectionConfig()
)
    vertex_data_type = _correction_vertex_data_type(graph)
    return _correct_observations(vertex_data_type, graph, observations; config = config)
end

function _correct_observations(
        ::Type{<:Any},
        graph::MetaGraphs.MetaDiGraph,
        observations;
        config::ViterbiCorrectionConfig
)::ViterbiCorrectionResult
    _assert_legacy_stranded_graph(graph)
    observed_paths = observations === nothing ? graph.gprops[:observed_paths] : observations
    corrected = viterbi_maximum_likelihood_traversals(
        graph;
        error_rate = config.error_rate,
        verbosity = config.verbosity == "silent" ? "dataset" : config.verbosity
    )
    diagnostics = Dict{Symbol, Any}(
        :interface => :legacy_stranded_kmer,
        :vertex_data_type => :legacy_stranded_kmer,
        :algorithm => :viterbi_maximum_likelihood_traversals,
        :observation_count => length(observed_paths),
        :emission_callback => nameof(config.emission_logp)
    )
    return ViterbiCorrectionResult(corrected, observed_paths, diagnostics)
end

function _correct_observations(
        ::Type{<:_VITERBI_CORRECTION_VERTEX_DATA},
        graph::MetaGraphsNext.MetaGraph,
        observations;
        config::ViterbiCorrectionConfig
)::ViterbiCorrectionResult
    return _correct_metagraphs_next_observations(graph, observations; config = config)
end

function _correct_observations(
        ::Type{Any},
        graph::MetaGraphsNext.MetaGraph,
        observations;
        config::ViterbiCorrectionConfig
)::ViterbiCorrectionResult
    return _correct_metagraphs_next_observations(graph, observations; config = config)
end

function _correct_metagraphs_next_observations(
        graph::MetaGraphsNext.MetaGraph,
        observations;
        config::ViterbiCorrectionConfig
)::ViterbiCorrectionResult
    weighted = if _correction_edge_data_type(graph) <: Rhizomorph.StrandWeightedEdgeData
        graph
    else
        Rhizomorph.weighted_graph_from_rhizomorph(graph; edge_weight = config.edge_weight)
    end
    paths = Vector{Rhizomorph.ViterbiDecodingResult}()

    for observation in observations
        start_vertex, target_vertex, max_steps = _decode_observation_bounds(
            observation,
            config
        )
        result = Rhizomorph.viterbi_decode_next(
            weighted,
            start_vertex,
            max_steps;
            target_vertex = target_vertex,
            start_strand = config.start_strand
        )
        push!(paths, result)
    end

    diagnostics = Dict{Symbol, Any}(
        :interface => :metagraphs_next,
        :vertex_data_type => _correction_vertex_data_type(graph),
        :algorithm => :rhizomorph_viterbi_decode_next,
        :observation_count => length(observations),
        :emission_callback => nameof(config.emission_logp)
    )
    return ViterbiCorrectionResult(paths, paths, diagnostics)
end

function _correction_vertex_data_type(
        graph::MetaGraphsNext.MetaGraph{CODE, GRAPH, LABEL, VERTEX_DATA}
)::Type where {CODE, GRAPH, LABEL, VERTEX_DATA}
    return VERTEX_DATA
end

function _correction_edge_data_type(
        graph::MetaGraphsNext.MetaGraph{CODE, GRAPH, LABEL, VERTEX_DATA, EDGE_DATA}
)::Type where {CODE, GRAPH, LABEL, VERTEX_DATA, EDGE_DATA}
    return EDGE_DATA
end

function _correction_vertex_data_type(graph)::Type
    return typeof(graph)
end

function _default_correction_observations(graph::MetaGraphs.MetaDiGraph)
    if haskey(graph.gprops, :observed_paths)
        return graph.gprops[:observed_paths]
    end
    return nothing
end

function _default_correction_observations(graph::MetaGraphsNext.MetaGraph)
    return collect(MetaGraphsNext.labels(graph))
end

function _default_correction_observations(graph)
    throw(ArgumentError("observations are required for $(typeof(graph))"))
end

function _assert_legacy_stranded_graph(graph::MetaGraphs.MetaDiGraph)::Nothing
    required = (:stranded_kmers, :observed_paths, :observation_ids, :k, :K)
    missing = [key for key in required if !haskey(graph.gprops, key)]
    if !isempty(missing)
        throw(ArgumentError("not a legacy stranded k-mer graph; missing graph props: $missing"))
    end
    return nothing
end

function _decode_observation_bounds(
        observation,
        config::ViterbiCorrectionConfig
)
    if observation isa Pair
        start_vertex = first(observation)
        target_vertex = last(observation)
        max_steps = config.max_steps === nothing ? 1 : config.max_steps
        return start_vertex, target_vertex, max_steps
    end

    if observation isa AbstractVector
        if isempty(observation)
            throw(ArgumentError("empty observation path"))
        end
        start_vertex = first(observation)
        target_vertex = config.target_vertex === nothing ? last(observation) : config.target_vertex
        max_steps = config.max_steps === nothing ? length(observation) - 1 : config.max_steps
        return start_vertex, target_vertex, max_steps
    end

    start_vertex = observation
    target_vertex = config.target_vertex
    max_steps = config.max_steps === nothing ? 0 : config.max_steps
    return start_vertex, target_vertex, max_steps
end

function _hamming_like_distance(left::AbstractString, right::AbstractString)::Int
    shared = min(length(left), length(right))
    distance = abs(length(left) - length(right))
    for index in 1:shared
        distance += left[index] == right[index] ? 0 : 1
    end
    return distance
end
