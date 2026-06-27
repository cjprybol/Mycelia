"""
$(DocStringExtensions.TYPEDSIGNATURES)

Default emission model for the graph-agnostic Viterbi correction interface.

The callback shape is intentionally broader than B1 needs so follow-on beads can
swap in alphabet- and quality-aware models without changing the dynamic
programming seam.
"""
function default_viterbi_emission_logp(
        observed_unit::Any,
        node::Any,
        alphabet::Symbol;
        quality = nothing,
        error_rate::Float64 = 0.01
)::Float64
    if error_rate <= 0.0 || error_rate >= 1.0
        throw(ArgumentError("error_rate must be in (0, 1), got $error_rate"))
    end

    resolved_alphabet = _normalize_viterbi_alphabet(alphabet)
    observed = uppercase(_viterbi_unit_string(observed_unit))
    expected = uppercase(_viterbi_unit_string(node))
    _assert_viterbi_unit_matches_alphabet(observed, resolved_alphabet, :observed)
    _assert_viterbi_unit_matches_alphabet(expected, resolved_alphabet, :expected)

    quality_scores = _normalize_viterbi_quality_scores(quality)
    alphabet_size = _viterbi_alphabet_size(resolved_alphabet)
    shared = min(length(observed), length(expected))

    logp = 0.0
    for index in 1:shared
        position_error_rate = _viterbi_position_error_rate(
            quality_scores,
            index,
            error_rate
        )
        if observed[index] == expected[index]
            logp += log1p(-position_error_rate)
        else
            logp += log(position_error_rate / (alphabet_size - 1))
        end
    end

    for index in (shared + 1):max(length(observed), length(expected))
        position_error_rate = _viterbi_position_error_rate(
            quality_scores,
            index,
            error_rate
        )
        logp += log(position_error_rate / alphabet_size)
    end
    return logp
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Configuration for `correct_observations`.

`emission_logp` is the extension seam for B2-B4: it receives
`(observed_unit, node, alphabet; quality)` and returns a log-probability. B1
keeps the default simple so the legacy B0 oracle remains the behavioral anchor.
`strand_mode` controls RC-aware nucleotide correction (`:singlestrand`,
`:doublestrand`, `:canonical`, or `:auto`). BioSequences empirically preserves
RNA type and U-aware complements (`AUGC` reverse-complements to `GCAU`), while
AA/text alphabets remain reverse-complement naive singlestrand paths.
"""
struct ViterbiCorrectionConfig{F <: Function}
    error_rate::Float64
    verbosity::String
    emission_logp::F
    alphabet::Symbol
    strand_mode::Symbol
    max_steps::Union{Nothing, Int}
    beam_width::Union{Nothing, Int}
    target_vertex::Any
    start_strand::Rhizomorph.StrandOrientation
    edge_weight::Function

    function ViterbiCorrectionConfig{F}(;
            error_rate::Float64 = 0.01,
            verbosity::String = "dataset",
            emission_logp::F,
            alphabet::Symbol = :auto,
            strand_mode::Symbol = :auto,
            max_steps::Union{Nothing, Int} = nothing,
            beam_width::Union{Nothing, Int} = nothing,
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
        if beam_width !== nothing && beam_width <= 0
            throw(ArgumentError("beam_width must be positive when set, got $beam_width"))
        end
        return new{F}(
            error_rate,
            verbosity,
            emission_logp,
            alphabet,
            strand_mode,
            max_steps,
            beam_width,
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
        strand_mode::Symbol = :auto,
        max_steps::Union{Nothing, Int} = nothing,
        beam_width::Union{Nothing, Int} = nothing,
        target_vertex = nothing,
        start_strand::Rhizomorph.StrandOrientation = Rhizomorph.Forward,
        edge_weight::Function = Rhizomorph.edge_data_weight
)::ViterbiCorrectionConfig{F} where {F <: Function}
    return ViterbiCorrectionConfig{F}(
        error_rate = error_rate,
        verbosity = verbosity,
        emission_logp = emission_logp,
        alphabet = alphabet,
        strand_mode = strand_mode,
        max_steps = max_steps,
        beam_width = beam_width,
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

const _VITERBI_FIXED_LENGTH_TEXT_VERTEX_DATA = Union{
    Rhizomorph.StringVertexData,
    Rhizomorph.QualityStringVertexData,
    Rhizomorph.LightweightStringVertexData,
    Rhizomorph.UltralightStringVertexData
}

const _VITERBI_VARIABLE_LENGTH_EDGE_DATA = Union{
    Rhizomorph.BioSequenceEdgeData,
    Rhizomorph.QualityBioSequenceEdgeData,
    Rhizomorph.StringEdgeData,
    Rhizomorph.QualityStringEdgeData
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
    alphabet = _resolve_viterbi_alphabet(graph, observations, config.alphabet)
    strand_mode = _resolve_viterbi_strand_mode(graph, alphabet, config.strand_mode)
    transition_edge_weight = _viterbi_transition_edge_weight(graph, config.edge_weight)
    weighted = if _correction_edge_data_type(graph) <: Rhizomorph.StrandWeightedEdgeData
        graph
    else
        Rhizomorph.weighted_graph_from_rhizomorph(graph; edge_weight = transition_edge_weight)
    end
    paths = Vector{Rhizomorph.ViterbiDecodingResult}()

    for observation in observations
        result = if _uses_emission_scored_observation(observation)
            _viterbi_correct_observation(
                weighted,
                observation,
                alphabet;
                config = config,
                strand_mode = strand_mode,
                quality_graph = graph,
                transition_scoring = _viterbi_transition_scoring(graph, transition_edge_weight)
            )
        else
            start_vertex, target_vertex, max_steps = _decode_observation_bounds(
                observation,
                config
            )
            if config.beam_width === nothing
                Rhizomorph.viterbi_decode_next(
                    weighted,
                    start_vertex,
                    max_steps;
                    target_vertex = target_vertex,
                    start_strand = config.start_strand
                )
            else
                Rhizomorph.beam_pruned_viterbi_decode_next(
                    weighted,
                    start_vertex,
                    max_steps;
                    target_vertex = target_vertex,
                    start_strand = config.start_strand,
                    beam_width = config.beam_width
                )
            end
        end
        push!(paths, result)
    end

    corrected = Any[_decoded_path_labels(path_result) for path_result in paths]
    diagnostics = Dict{Symbol, Any}(
        :interface => :metagraphs_next,
        :vertex_data_type => _correction_vertex_data_type(graph),
        :algorithm => :rhizomorph_viterbi_decode_next,
        :observation_count => length(observations),
        :emission_callback => nameof(config.emission_logp),
        :alphabet => alphabet,
        :emission_model => _viterbi_graph_has_quality(graph) ?
                           :quality_aware : :alphabet_parameterized,
        :transition_model => _viterbi_transition_scoring(graph, transition_edge_weight),
        :strand_mode => strand_mode,
        :reverse_complement_support => _viterbi_supports_reverse_complement(alphabet),
        :beam_width => config.beam_width,
        :expanded_states => _sum_path_diagnostic(paths, :expanded_states),
        :generated_states => _sum_path_diagnostic(paths, :generated_states),
        :max_retained_states => _max_path_diagnostic(paths, :max_retained_states),
        :cumulative_retained_states => _sum_path_diagnostic(
            paths,
            :cumulative_retained_states
        )
    )
    return ViterbiCorrectionResult(corrected, paths, diagnostics)
end


function _sum_path_diagnostic(
        paths::AbstractVector{Rhizomorph.ViterbiDecodingResult},
        key::Symbol
)::Int
    return sum(Int(get(path.diagnostics, key, 0)) for path in paths)
end

function _max_path_diagnostic(
        paths::AbstractVector{Rhizomorph.ViterbiDecodingResult},
        key::Symbol
)::Int
    if isempty(paths)
        return 0
    end
    return maximum(Int(get(path.diagnostics, key, 0)) for path in paths)
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

function _viterbi_has_variable_length_edges(graph::MetaGraphsNext.MetaGraph)::Bool
    if _correction_vertex_data_type(graph) <: _VITERBI_FIXED_LENGTH_TEXT_VERTEX_DATA &&
       _is_fixed_length_text_graph(graph)
        return false
    end
    return _correction_edge_data_type(graph) <: _VITERBI_VARIABLE_LENGTH_EDGE_DATA
end

function _viterbi_overlap_edge_weight(edge_data)::Float64
    if hasproperty(edge_data, :overlap_length)
        overlap_length = getproperty(edge_data, :overlap_length)
        if overlap_length > 0
            return Float64(overlap_length)
        end
    end
    return Float64(Rhizomorph.edge_data_weight(edge_data))
end

function _viterbi_transition_edge_weight(
        graph::MetaGraphsNext.MetaGraph,
        configured_edge_weight::Function
)::Function
    if _viterbi_has_variable_length_edges(graph) &&
       configured_edge_weight === Rhizomorph.edge_data_weight
        return _viterbi_overlap_edge_weight
    end
    return configured_edge_weight
end

function _viterbi_transition_scoring(
        graph::MetaGraphsNext.MetaGraph,
        edge_weight::Function
)::Symbol
    if _viterbi_has_variable_length_edges(graph) &&
       edge_weight === _viterbi_overlap_edge_weight
        return :normalized_overlap_length
    end
    return :normalized_edge_weight
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


function _normalize_viterbi_alphabet(alphabet::Symbol)::Symbol
    normalized = Symbol(uppercase(String(alphabet)))
    if !(normalized in (:DNA, :RNA, :AA, :TEXT))
        throw(ArgumentError("unsupported Viterbi correction alphabet: $alphabet"))
    end
    return normalized
end


function _normalize_viterbi_strand_mode(strand_mode::Symbol)::Symbol
    normalized = Symbol(lowercase(String(strand_mode)))
    if normalized in (:auto, :singlestrand, :single)
        return normalized == :auto ? :auto : :singlestrand
    elseif normalized in (:doublestrand, :double)
        return :doublestrand
    elseif normalized == :canonical
        return :canonical
    end
    throw(ArgumentError("unsupported Viterbi correction strand_mode: $strand_mode"))
end

function _viterbi_supports_reverse_complement(alphabet::Symbol)::Bool
    normalized = _normalize_viterbi_alphabet(alphabet)
    return normalized in (:DNA, :RNA)
end

function _resolve_viterbi_strand_mode(
        graph::MetaGraphsNext.MetaGraph,
        alphabet::Symbol,
        configured_mode::Symbol
)::Symbol
    requested = _normalize_viterbi_strand_mode(configured_mode)
    if !_viterbi_supports_reverse_complement(alphabet)
        if requested in (:auto, :singlestrand)
            return :singlestrand
        end
        throw(ArgumentError(
            "strand_mode $configured_mode requires a DNA/RNA alphabet; " *
            "alphabet $alphabet is reverse-complement naive"
        ))
    end

    if requested != :auto
        return requested
    end
    if !Graphs.is_directed(graph.graph)
        return :canonical
    end
    if _viterbi_graph_has_reverse_strand_evidence(graph)
        return :doublestrand
    end
    return :singlestrand
end

function _viterbi_graph_has_reverse_strand_evidence(graph::MetaGraphsNext.MetaGraph)::Bool
    for edge_label in MetaGraphsNext.edge_labels(graph)
        edge_data = graph[edge_label...]
        if edge_data isa Rhizomorph.StrandWeightedEdgeData
            if Rhizomorph._normalize_strand(edge_data.src_strand) == Rhizomorph.Reverse ||
               Rhizomorph._normalize_strand(edge_data.dst_strand) == Rhizomorph.Reverse
                return true
            end
        elseif hasproperty(edge_data, :evidence)
            strand = Rhizomorph.first_evidence_strand(
                getproperty(edge_data, :evidence);
                default = Rhizomorph.Forward
            )
            if Rhizomorph._normalize_strand(strand) == Rhizomorph.Reverse
                return true
            end
        end
    end
    return false
end

function _viterbi_alphabet_size(alphabet::Symbol)::Int
    normalized = _normalize_viterbi_alphabet(alphabet)
    if normalized in (:DNA, :RNA)
        return 4
    elseif normalized == :AA
        return 20
    elseif normalized == :TEXT
        return 256
    end
    throw(ArgumentError("unsupported Viterbi correction alphabet: $alphabet"))
end

function _normalize_viterbi_quality_scores(
        quality
)::Union{Nothing, Vector{Float64}}
    if quality === nothing
        return nothing
    elseif quality isa Number
        return Float64[Float64(quality)]
    elseif quality isa AbstractVector
        return Float64.(quality)
    elseif hasproperty(quality, :quality_scores)
        return _normalize_viterbi_quality_scores(
            getproperty(quality, :quality_scores)
        )
    end
    return nothing
end

function _viterbi_position_error_rate(
        quality_scores::Union{Nothing, Vector{Float64}},
        index::Int,
        fallback_error_rate::Float64
)::Float64
    if quality_scores === nothing || isempty(quality_scores)
        return fallback_error_rate
    end
    quality_index = min(index, length(quality_scores))
    phred = max(0.0, quality_scores[quality_index])
    return clamp(10.0^(-phred / 10.0), eps(Float64), 1.0 - eps(Float64))
end

function _viterbi_graph_has_quality(graph::MetaGraphsNext.MetaGraph)::Bool
    vertex_type = _correction_vertex_data_type(graph)
    return vertex_type <: Union{
        Rhizomorph.QualmerVertexData,
        Rhizomorph.QualityBioSequenceVertexData,
        Rhizomorph.QualityStringVertexData,
        Rhizomorph.UltralightQualityKmerVertexData,
        Rhizomorph.UltralightQualityBioSequenceVertexData,
        Rhizomorph.LightweightQualityKmerVertexData,
        Rhizomorph.LightweightQualityBioSequenceVertexData
    }
end

function _viterbi_emission_quality(
        graph::MetaGraphsNext.MetaGraph,
        observed_unit::Any
)::Union{Nothing, Vector{Float64}}
    direct_quality = _viterbi_direct_quality_scores(observed_unit)
    if direct_quality !== nothing
        return direct_quality
    end

    if !_viterbi_graph_has_quality(graph)
        return nothing
    end

    if haskey(graph, observed_unit)
        return _viterbi_vertex_quality_scores(graph[observed_unit])
    end
    return nothing
end

function _viterbi_direct_quality_scores(unit::Any)::Union{Nothing, Vector{Float64}}
    if hasproperty(unit, :quality_scores)
        return _viterbi_decode_quality_scores(getproperty(unit, :quality_scores))
    end
    return nothing
end

function _viterbi_vertex_quality_scores(
        vertex_data::Any
)::Union{Nothing, Vector{Float64}}
    evidence_quality = _viterbi_evidence_joint_quality_scores(vertex_data)
    if evidence_quality !== nothing
        return evidence_quality
    end

    if hasproperty(vertex_data, :joint_quality)
        joint_quality = getproperty(vertex_data, :joint_quality)
        if !isempty(joint_quality)
            return Float64.(joint_quality)
        end
    end

    if hasproperty(vertex_data, :quality_scores)
        return _viterbi_decode_quality_scores(getproperty(vertex_data, :quality_scores))
    end
    return nothing
end

function _viterbi_evidence_joint_quality_scores(
        vertex_data::Any
)::Union{Nothing, Vector{Float64}}
    quality_vectors = Vector{Vector{Float64}}()
    dataset_ids = _viterbi_dataset_ids(vertex_data)
    for dataset_id in dataset_ids
        observations = _viterbi_observation_ids(vertex_data, dataset_id)
        if observations === nothing
            dataset_quality = _viterbi_dataset_joint_quality(vertex_data, dataset_id)
            dataset_quality === nothing || push!(quality_vectors, dataset_quality)
            continue
        end

        for observation_id in observations
            evidence = Rhizomorph.get_observation_evidence(
                vertex_data,
                dataset_id,
                observation_id
            )
            evidence === nothing && continue
            for entry in evidence
                if entry isa Rhizomorph.QualityEvidenceEntry
                    push!(
                        quality_vectors,
                        _viterbi_decode_quality_scores(entry.quality_scores)
                    )
                end
            end
        end
    end

    if isempty(quality_vectors)
        return nothing
    end
    return _viterbi_combine_quality_vectors(quality_vectors)
end

function _viterbi_dataset_joint_quality(
        vertex_data::Any,
        dataset_id::AbstractString
)::Union{Nothing, Vector{Float64}}
    try
        joint_quality = Rhizomorph.get_vertex_joint_quality(
            vertex_data,
            String(dataset_id)
        )
        return joint_quality === nothing ? nothing : Float64.(joint_quality)
    catch error
        if error isa MethodError
            return nothing
        end
        rethrow()
    end
end

function _viterbi_dataset_ids(vertex_data::Any)::Vector{String}
    try
        return String.(Rhizomorph.get_all_dataset_ids(vertex_data))
    catch error
        if error isa MethodError
            return String[]
        end
        rethrow()
    end
end

function _viterbi_observation_ids(
        vertex_data::Any,
        dataset_id::AbstractString
)::Union{Nothing, Vector{String}}
    try
        observations = Rhizomorph.get_all_observation_ids(
            vertex_data,
            String(dataset_id)
        )
        return observations === nothing ? nothing : String.(observations)
    catch error
        if error isa MethodError
            return nothing
        end
        rethrow()
    end
end

function _viterbi_decode_quality_scores(scores::AbstractVector)::Vector{Float64}
    if isempty(scores)
        return Float64[]
    end
    numeric_scores = Float64.(scores)
    if minimum(numeric_scores) >= 33.0 && maximum(numeric_scores) <= 126.0
        return numeric_scores .- 33.0
    end
    return numeric_scores
end

function _viterbi_combine_quality_vectors(
        quality_vectors::Vector{Vector{Float64}}
)::Vector{Float64}
    max_length = maximum(length.(quality_vectors))
    joint_quality = Vector{Float64}(undef, max_length)
    for index in 1:max_length
        joint_quality[index] = min(
            sum(vector[index] for vector in quality_vectors if index <= length(vector)),
            255.0
        )
    end
    return joint_quality
end

function _viterbi_unit_string(unit::Any)::String
    if hasproperty(unit, :Kmer)
        return string(getproperty(unit, :Kmer))
    elseif hasproperty(unit, :sequence)
        return string(getproperty(unit, :sequence))
    end
    return string(unit)
end

function _assert_viterbi_unit_matches_alphabet(
        unit::AbstractString,
        alphabet::Symbol,
        role::Symbol
)::Nothing
    if _normalize_viterbi_alphabet(alphabet) == :TEXT
        return nothing
    end
    if !validate_alphabet(unit, alphabet)
        throw(ArgumentError("$role unit $unit is not valid for alphabet $alphabet"))
    end
    return nothing
end

function _resolve_viterbi_alphabet(
        graph::MetaGraphsNext.MetaGraph,
        observations::Any,
        configured_alphabet::Symbol
)::Symbol
    if configured_alphabet != :auto
        return _normalize_viterbi_alphabet(configured_alphabet)
    end

    graph_alphabet = _infer_viterbi_graph_alphabet(graph)
    graph_alphabet === nothing || return graph_alphabet

    for label in MetaGraphsNext.labels(graph)
        inferred = _infer_viterbi_alphabet(label)
        inferred === nothing || return inferred
    end
    for observation in observations
        inferred = _infer_viterbi_alphabet(observation)
        inferred === nothing || return inferred
    end
    throw(ArgumentError("could not infer Viterbi correction alphabet; pass config.alphabet"))
end

function _infer_viterbi_graph_alphabet(
        graph::MetaGraphsNext.MetaGraph
)::Union{Nothing, Symbol}
    vertex_type = _correction_vertex_data_type(graph)
    if vertex_type <: _VITERBI_FIXED_LENGTH_TEXT_VERTEX_DATA
        return :TEXT
    end
    return nothing
end

function _is_fixed_length_text_graph(graph::MetaGraphsNext.MetaGraph)::Bool
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels) || !all(label -> label isa AbstractString, labels)
        return false
    end

    label_lengths = unique(length.(labels))
    if length(label_lengths) != 1
        return false
    end

    expected_overlap = only(label_lengths) - 1
    for edge_label in MetaGraphsNext.edge_labels(graph)
        edge_data = graph[edge_label...]
        if !hasproperty(edge_data, :overlap_length) ||
           getproperty(edge_data, :overlap_length) != expected_overlap
            return false
        end
    end
    return true
end

function _infer_viterbi_alphabet(unit::Any)::Union{Nothing, Symbol}
    if unit isa AbstractVector
        for item in unit
            inferred = _infer_viterbi_alphabet(item)
            inferred === nothing || return inferred
        end
        return nothing
    end

    unit_type = typeof(unit)
    if unit_type <: Kmers.DNAKmer
        return :DNA
    elseif unit_type <: Kmers.RNAKmer
        return :RNA
    elseif unit_type <: Kmers.AAKmer
        return :AA
    elseif unit isa BioSequences.LongDNA
        return :DNA
    elseif unit isa BioSequences.LongRNA
        return :RNA
    elseif unit isa BioSequences.LongAA
        return :AA
    elseif hasproperty(unit, :Kmer)
        return _infer_viterbi_alphabet(getproperty(unit, :Kmer))
    elseif hasproperty(unit, :sequence)
        return _infer_viterbi_alphabet(getproperty(unit, :sequence))
    end

    try
        return detect_alphabet(uppercase(_viterbi_unit_string(unit)))
    catch error
        if error isa ArgumentError
            return :TEXT
        end
        rethrow()
    end
end

function _uses_emission_scored_observation(observation::Any)::Bool
    return observation isa AbstractVector && !isempty(observation)
end

function _call_viterbi_emission_logp(
        config::ViterbiCorrectionConfig,
        observed_unit::Any,
        node::Any,
        alphabet::Symbol;
        quality = nothing
)::Float64
    try
        return config.emission_logp(
            observed_unit,
            node,
            alphabet;
            quality = quality,
            error_rate = config.error_rate
        )
    catch error
        if error isa MethodError
            return config.emission_logp(observed_unit, node, alphabet; quality = quality)
        end
        rethrow()
    end
end

function _viterbi_correct_observation(
        graph::MetaGraphsNext.MetaGraph,
        observation::AbstractVector,
        alphabet::Symbol;
        config::ViterbiCorrectionConfig,
        strand_mode::Symbol,
        quality_graph::MetaGraphsNext.MetaGraph = graph,
        transition_scoring::Symbol = :normalized_edge_weight
)::Rhizomorph.ViterbiDecodingResult
    if isempty(observation)
        throw(ArgumentError("empty observation path"))
    end

    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        throw(ArgumentError("cannot correct observations against an empty graph"))
    end

    label_type = eltype(labels)
    start_observed = first(observation)
    target_vertex = _emission_target_vertex(graph, observation, config, alphabet, strand_mode)
    start_candidates = _viterbi_start_candidates(labels, start_observed, alphabet, strand_mode)

    diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_emission_correct_observation,
        :exact => true,
        :alphabet => alphabet,
        :strand_mode => strand_mode,
        :reverse_complement_support => _viterbi_supports_reverse_complement(alphabet),
        :max_steps => length(observation) - 1,
        :beam_width => config.beam_width,
        :target_vertex => target_vertex,
        :start_strand => Rhizomorph._normalize_strand(config.start_strand),
        :score_domain => :log_probability,
        :transition_scoring => transition_scoring,
        :emission_scoring => _viterbi_graph_has_quality(quality_graph) ?
                             :quality_aware : :alphabet_parameterized,
        :expanded_states => 0,
        :generated_states => 0,
        :retained_states => 0,
        :cumulative_retained_states => 0,
        :max_retained_states => 0,
        :skipped_transitions => 0,
        :completed_steps => 0,
        :reached_target => target_vertex === nothing ? nothing : false
    )

    active_scores = Dict{Tuple{label_type, Rhizomorph.StrandOrientation}, Float64}()
    for vertex in start_candidates
        for strand in _viterbi_start_strands(graph, vertex, strand_mode, config.start_strand)
            state = (vertex, strand)
            score = _call_viterbi_state_emission_logp(
                quality_graph,
                config,
                start_observed,
                vertex,
                alphabet,
                strand_mode
            )
            if isfinite(score) && (!haskey(active_scores, state) || score > active_scores[state])
                active_scores[state] = score
            end
        end
    end
    if isempty(active_scores)
        diagnostics[:reason] = :no_finite_start_emission
        return Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diagnostics)
    end

    diagnostics[:retained_states] = length(active_scores)
    diagnostics[:cumulative_retained_states] = length(active_scores)
    diagnostics[:max_retained_states] = length(active_scores)

    best_state, best_score = _best_correction_state(active_scores)
    best_depth = 0
    if target_vertex !== nothing
        if length(observation) == 1
            target_state, target_score = _best_correction_target_state(active_scores, target_vertex)
            if target_state !== nothing
                best_state = target_state
                best_score = target_score
                diagnostics[:reached_target] = true
            else
                best_score = -Inf
            end
        else
            best_score = -Inf
        end
    end

    predecessors_by_depth = Vector{
        Dict{
            Tuple{label_type, Rhizomorph.StrandOrientation},
            Tuple{label_type, Rhizomorph.StrandOrientation}
        }
    }()

    for depth in 1:(length(observation) - 1)
        observed_unit = observation[depth + 1]
        next_scores = Dict{Tuple{label_type, Rhizomorph.StrandOrientation}, Float64}()
        next_predecessors = Dict{
            Tuple{label_type, Rhizomorph.StrandOrientation},
            Tuple{label_type, Rhizomorph.StrandOrientation}
        }()

        for (state, state_score) in active_scores
            current_vertex, current_strand = state
            transitions = Rhizomorph._get_valid_transitions(graph, current_vertex, current_strand)
            diagnostics[:expanded_states] += 1
            if isempty(transitions)
                continue
            end

            total_out = Rhizomorph._total_outgoing_weight(graph, current_vertex, current_strand)
            if !isfinite(total_out) || total_out <= 0.0
                diagnostics[:skipped_transitions] += length(transitions)
                continue
            end

            for transition in transitions
                next_vertex = convert(label_type, transition[:target_vertex])
                next_strand = Rhizomorph._normalize_strand(transition[:target_strand])
                edge_w = Rhizomorph._edge_transition_weight(transition[:edge_data])
                if edge_w <= 0.0
                    diagnostics[:skipped_transitions] += 1
                    continue
                end
                transition_prob = edge_w / total_out
                emission_score = _call_viterbi_state_emission_logp(
                    quality_graph,
                    config,
                    observed_unit,
                    next_vertex,
                    alphabet,
                    strand_mode
                )
                next_score = state_score + log(transition_prob) + emission_score
                if !isfinite(next_score)
                    diagnostics[:skipped_transitions] += 1
                    continue
                end

                next_state = (next_vertex, next_strand)
                diagnostics[:generated_states] += 1
                if !haskey(next_scores, next_state) || next_score > next_scores[next_state]
                    next_scores[next_state] = next_score
                    next_predecessors[next_state] = state
                end
            end
        end

        if isempty(next_scores)
            break
        end

        if config.beam_width !== nothing
            next_scores, next_predecessors = _prune_correction_beam(
                next_scores,
                next_predecessors,
                config.beam_width
            )
        end

        push!(predecessors_by_depth, next_predecessors)
        active_scores = next_scores
        retained_count = length(active_scores)
        diagnostics[:retained_states] = retained_count
        diagnostics[:cumulative_retained_states] += retained_count
        diagnostics[:max_retained_states] = max(diagnostics[:max_retained_states], retained_count)
        diagnostics[:completed_steps] = depth

        if target_vertex === nothing
            best_state, best_score = _best_correction_state(active_scores)
            best_depth = depth
        else
            target_state, target_score = _best_correction_target_state(active_scores, target_vertex)
            if target_state !== nothing
                best_state = target_state
                best_score = target_score
                best_depth = depth
                diagnostics[:reached_target] = true
            end
        end
    end

    if target_vertex !== nothing && !isfinite(best_score)
        diagnostics[:reason] = :target_unreachable
        return Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diagnostics)
    end

    path = _reconstruct_correction_path(graph, best_state, best_depth, predecessors_by_depth)
    diagnostics[:path_length] = length(path.steps)
    return Rhizomorph.ViterbiDecodingResult(path, best_score, diagnostics)
end

function _emission_target_vertex(
        graph::MetaGraphsNext.MetaGraph,
        observation::AbstractVector,
        config::ViterbiCorrectionConfig,
        alphabet::Symbol,
        strand_mode::Symbol
)
    if config.target_vertex !== nothing
        return config.target_vertex
    end
    candidate = last(observation)
    labels = collect(MetaGraphsNext.labels(graph))
    if candidate in labels
        return candidate
    end
    if strand_mode == :canonical && _viterbi_supports_reverse_complement(alphabet)
        canonical = _viterbi_canonical_unit(candidate, alphabet)
        if canonical in labels
            return canonical
        end
    end
    return nothing
end

function _viterbi_start_candidates(
        labels::Vector{T},
        observed_unit::Any,
        alphabet::Symbol,
        strand_mode::Symbol
)::Vector{T} where {T}
    if observed_unit in labels
        return T[convert(T, observed_unit)]
    end
    if strand_mode == :canonical && _viterbi_supports_reverse_complement(alphabet)
        canonical = _viterbi_canonical_unit(observed_unit, alphabet)
        if canonical in labels
            return T[convert(T, canonical)]
        end
    end
    return labels
end

function _viterbi_start_strands(
        graph::MetaGraphsNext.MetaGraph,
        vertex,
        strand_mode::Symbol,
        configured_start_strand
)::Vector{Rhizomorph.StrandOrientation}
    configured = Rhizomorph._normalize_strand(configured_start_strand)
    if strand_mode == :singlestrand
        return Rhizomorph.StrandOrientation[configured]
    end

    outgoing = _viterbi_outgoing_strands(graph, vertex)
    if strand_mode == :canonical
        if configured in outgoing || isempty(outgoing)
            return Rhizomorph.StrandOrientation[configured]
        end
        return Rhizomorph.StrandOrientation[first(outgoing)]
    end

    if isempty(outgoing)
        other = configured == Rhizomorph.Forward ? Rhizomorph.Reverse : Rhizomorph.Forward
        return Rhizomorph.StrandOrientation[configured, other]
    end
    return collect(outgoing)
end

function _viterbi_outgoing_strands(
        graph::MetaGraphsNext.MetaGraph,
        vertex
)::Set{Rhizomorph.StrandOrientation}
    strands = Set{Rhizomorph.StrandOrientation}()
    haskey(graph, vertex) || return strands
    if Graphs.is_directed(graph.graph)
        src_code = MetaGraphsNext.code_for(graph, vertex)
        for dst_code in Graphs.outneighbors(graph.graph, src_code)
            target_vertex = MetaGraphsNext.label_for(graph, dst_code)
            edge_data = graph[vertex, target_vertex]
            push!(strands, Rhizomorph._normalize_strand(edge_data.src_strand))
        end
    else
        for edge_label in MetaGraphsNext.edge_labels(graph)
            if length(edge_label) == 2 && edge_label[1] == vertex
                edge_data = graph[edge_label...]
                push!(strands, Rhizomorph._normalize_strand(edge_data.src_strand))
            end
        end
    end
    return strands
end

function _call_viterbi_state_emission_logp(
        graph::MetaGraphsNext.MetaGraph,
        config::ViterbiCorrectionConfig,
        observed_unit::Any,
        node::Any,
        alphabet::Symbol,
        strand_mode::Symbol
)::Float64
    quality = _viterbi_emission_quality(graph, observed_unit)
    direct = _call_viterbi_emission_logp(
        config,
        observed_unit,
        node,
        alphabet;
        quality = quality
    )
    if strand_mode == :canonical && _viterbi_supports_reverse_complement(alphabet)
        rc_node = _viterbi_reverse_complement_unit(node, alphabet)
        rc = _call_viterbi_emission_logp(
            config,
            observed_unit,
            rc_node,
            alphabet;
            quality = quality
        )
        return max(direct, rc)
    end
    return direct
end

function _viterbi_canonical_unit(unit::Any, alphabet::Symbol)
    if !_viterbi_supports_reverse_complement(alphabet)
        return unit
    end
    try
        return BioSequences.canonical(unit)
    catch error
        if !(error isa MethodError)
            rethrow()
        end
    end
    unit_string = _viterbi_unit_string(unit)
    rc_string = _viterbi_reverse_complement_string(unit_string, alphabet)
    return unit_string <= rc_string ? unit_string : rc_string
end

function _viterbi_reverse_complement_unit(unit::Any, alphabet::Symbol)
    if !_viterbi_supports_reverse_complement(alphabet)
        return unit
    end
    try
        return BioSequences.reverse_complement(unit)
    catch error
        if !(error isa MethodError)
            rethrow()
        end
    end
    return _viterbi_reverse_complement_string(_viterbi_unit_string(unit), alphabet)
end

function _viterbi_reverse_complement_string(
        sequence::AbstractString,
        alphabet::Symbol
)::String
    normalized = _normalize_viterbi_alphabet(alphabet)
    complement = if normalized == :DNA
        Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'N' => 'N')
    elseif normalized == :RNA
        Dict('A' => 'U', 'C' => 'G', 'G' => 'C', 'U' => 'A', 'N' => 'N')
    else
        throw(ArgumentError("alphabet $alphabet has no reverse complement"))
    end
    return join((complement[base] for base in reverse(uppercase(sequence))))
end

function _prune_correction_beam(
        scores::Dict{Tuple{T, Rhizomorph.StrandOrientation}, Float64},
        predecessors::Dict{
            Tuple{T, Rhizomorph.StrandOrientation},
            Tuple{T, Rhizomorph.StrandOrientation}
        },
        beam_width::Int
)::Tuple{
        Dict{Tuple{T, Rhizomorph.StrandOrientation}, Float64},
        Dict{Tuple{T, Rhizomorph.StrandOrientation}, Tuple{T, Rhizomorph.StrandOrientation}}
} where {T}
    if length(scores) <= beam_width
        return scores, predecessors
    end

    retained = first(
        sort(
            collect(scores);
            by = item -> (-item[2], string(item[1][1]), Int(item[1][2]))
        ),
        beam_width
    )
    retained_scores = Dict{Tuple{T, Rhizomorph.StrandOrientation}, Float64}(retained)
    retained_predecessors = Dict{
        Tuple{T, Rhizomorph.StrandOrientation},
        Tuple{T, Rhizomorph.StrandOrientation}
    }()
    for (state, _) in retained
        if haskey(predecessors, state)
            retained_predecessors[state] = predecessors[state]
        end
    end

    return retained_scores, retained_predecessors
end

function _best_correction_state(
        scores::Dict{Tuple{T, Rhizomorph.StrandOrientation}, Float64}
) where {T}
    best_state = nothing
    best_score = -Inf
    for (state, score) in scores
        if best_state === nothing || score > best_score ||
           (isapprox(score, best_score; atol = eps(Float64), rtol = eps(Float64)) &&
            string(state[1]) < string(best_state[1]))
            best_state = state
            best_score = score
        end
    end
    return best_state, best_score
end

function _best_correction_target_state(
        scores::Dict{Tuple{T, Rhizomorph.StrandOrientation}, Float64},
        target_vertex
) where {T}
    target_scores = Dict{Tuple{T, Rhizomorph.StrandOrientation}, Float64}()
    for (state, score) in scores
        if state[1] == target_vertex
            target_scores[state] = score
        end
    end
    if isempty(target_scores)
        return nothing, -Inf
    end
    return _best_correction_state(target_scores)
end

function _reconstruct_correction_path(
        graph::MetaGraphsNext.MetaGraph,
        end_state::Tuple{T, Rhizomorph.StrandOrientation},
        depth::Int,
        predecessors_by_depth::Vector{
            Dict{
                Tuple{T, Rhizomorph.StrandOrientation},
                Tuple{T, Rhizomorph.StrandOrientation}
            }
        }
) where {T}
    path_states = Vector{Tuple{T, Rhizomorph.StrandOrientation}}(undef, depth + 1)
    current_state = end_state
    for path_index in reverse(2:(depth + 1))
        path_states[path_index] = current_state
        current_state = predecessors_by_depth[path_index - 1][current_state]
    end
    path_states[1] = current_state
    return Rhizomorph._build_graph_path_from_vertices(graph, path_states)
end

function _decoded_path_labels(
        result::Rhizomorph.ViterbiDecodingResult
)::Union{Nothing, Vector}
    if result.path === nothing
        return nothing
    end
    return [step.vertex_label for step in something(result.path).steps]
end

function _hamming_like_distance(left::AbstractString, right::AbstractString)::Int
    shared = min(length(left), length(right))
    distance = abs(length(left) - length(right))
    for index in 1:shared
        distance += left[index] == right[index] ? 0 : 1
    end
    return distance
end
