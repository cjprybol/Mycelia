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

    shared = min(length(observed), length(expected))
    substitution_logp = log(error_rate / (_viterbi_alphabet_size(resolved_alphabet) - 1))
    indel_logp = log(error_rate / _viterbi_alphabet_size(resolved_alphabet))
    match_logp = log1p(-error_rate)

    logp = 0.0
    for index in 1:shared
        logp += observed[index] == expected[index] ? match_logp : substitution_logp
    end
    logp += abs(length(observed) - length(expected)) * indel_logp
    return logp
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
    alphabet = _resolve_viterbi_alphabet(graph, observations, config.alphabet)
    paths = Vector{Rhizomorph.ViterbiDecodingResult}()

    for observation in observations
        result = if _uses_emission_scored_observation(observation)
            _viterbi_correct_observation(weighted, observation, alphabet; config = config)
        else
            start_vertex, target_vertex, max_steps = _decode_observation_bounds(
                observation,
                config
            )
            Rhizomorph.viterbi_decode_next(
                weighted,
                start_vertex,
                max_steps;
                target_vertex = target_vertex,
                start_strand = config.start_strand
            )
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
        :emission_model => :alphabet_parameterized
    )
    return ViterbiCorrectionResult(corrected, paths, diagnostics)
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


function _normalize_viterbi_alphabet(alphabet::Symbol)::Symbol
    normalized = Symbol(uppercase(String(alphabet)))
    if !(normalized in (:DNA, :RNA, :AA))
        throw(ArgumentError("unsupported Viterbi correction alphabet: $alphabet"))
    end
    return normalized
end

function _viterbi_alphabet_size(alphabet::Symbol)::Int
    normalized = _normalize_viterbi_alphabet(alphabet)
    if normalized in (:DNA, :RNA)
        return 4
    elseif normalized == :AA
        return 20
    end
    throw(ArgumentError("unsupported Viterbi correction alphabet: $alphabet"))
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

    return detect_alphabet(uppercase(_viterbi_unit_string(unit)))
end

function _uses_emission_scored_observation(observation::Any)::Bool
    return observation isa AbstractVector && !isempty(observation)
end

function _call_viterbi_emission_logp(
        config::ViterbiCorrectionConfig,
        observed_unit::Any,
        node::Any,
        alphabet::Symbol
)::Float64
    try
        return config.emission_logp(
            observed_unit,
            node,
            alphabet;
            quality = nothing,
            error_rate = config.error_rate
        )
    catch error
        if error isa MethodError
            return config.emission_logp(observed_unit, node, alphabet; quality = nothing)
        end
        rethrow()
    end
end

function _viterbi_correct_observation(
        graph::MetaGraphsNext.MetaGraph,
        observation::AbstractVector,
        alphabet::Symbol;
        config::ViterbiCorrectionConfig
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
    target_vertex = _emission_target_vertex(graph, observation, config)
    start_candidates = if start_observed in labels
        label_type[convert(label_type, start_observed)]
    else
        labels
    end

    diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_emission_correct_observation,
        :exact => true,
        :alphabet => alphabet,
        :max_steps => length(observation) - 1,
        :target_vertex => target_vertex,
        :start_strand => Rhizomorph._normalize_strand(config.start_strand),
        :score_domain => :log_probability,
        :transition_scoring => :normalized_edge_weight,
        :emission_scoring => :alphabet_parameterized,
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
        state = (vertex, Rhizomorph._normalize_strand(config.start_strand))
        score = _call_viterbi_emission_logp(config, start_observed, vertex, alphabet)
        if isfinite(score)
            active_scores[state] = score
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
        target_state, target_score = _best_correction_target_state(active_scores, target_vertex)
        if target_state !== nothing
            best_state = target_state
            best_score = target_score
            diagnostics[:reached_target] = true
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
                emission_score = _call_viterbi_emission_logp(
                    config,
                    observed_unit,
                    next_vertex,
                    alphabet
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
            if target_state !== nothing && target_score > best_score
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
        config::ViterbiCorrectionConfig
)
    if config.target_vertex !== nothing
        return config.target_vertex
    end
    candidate = last(observation)
    return candidate in MetaGraphsNext.labels(graph) ? candidate : nothing
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
