# Batched array-frontier Viterbi corrector — CPU proof-of-concept (bead td-qoo3).
#
# This file is the greenfield foundation for the GPU/SIMD corrector-acceleration
# path described in `docs/design/2026-07-06-gpu-simd-corrector-acceleration.md`.
# It implements the ADR's "smallest batched-decode proof-of-concept" (§Prototype
# scope) and NOTHING beyond it: a CPU-only, single-threaded, byte-identical
# reformulation of the scalar per-read decoder (`_viterbi_correct_observation`)
# whose only job is to prove the reformulated recurrence reproduces the scalar
# decoder before any SIMD (Phase A) or GPU kernel (Phase B) work is attempted.
#
# What the ADR asks the PoC to change, and where it happens here:
#
#   1. Array-based frontier (not a Dict). The active frontier is a dense
#      `n_reads x n_states` `Float64` matrix (`-Inf` for unreached states),
#      indexed by a stable contiguous `state_id <- (vertex, strand)` map precomputed
#      once for the immutable graph. Dense rows are SIMD-scannable / GPU-uploadable
#      / coalesced-addressable — the single change that unlocks Phases A and B.
#      (ADR §Decision item 1.)
#
#   2. Depth-outer, read-inner loop nesting (`for depth { for read { relax } }`).
#      The read axis becomes the innermost, fully data-parallel axis; every read in
#      a batch advances one read-position in lock-step against the shared graph.
#      This is the loop-interchange that exposes the zero-write-contention across-
#      read parallelism to the hardware. (ADR §Decision item 2.)
#
#   3. Length-binned batches. Reads are grouped by unit-count so every read in a
#      bin shares a common depth bound; each bin is one lock-step sweep. (ADR
#      §Decision item 3.)
#
#   4. Shared immutable graph as a structure-of-arrays, built ONCE per batch. Per
#      source state we precompute its out-edges as `(next_state_id,
#      log_transition_prob, next_vertex, next_strand)`; `log_transition_prob` is a
#      pure function of the graph (edge weight / total outgoing weight), so it is
#      hoisted out of the hot loop. All lanes read it; none writes it. Diagnostics
#      counters are per-read (per-lane), never a shared mutable scalar in the
#      inner kernel. (ADR §Decision item 4.)
#
# Byte-identical guarantee. Every scoring/topology primitive is the SAME function
# the scalar decoder calls (`_call_viterbi_state_emission_logp`,
# `Rhizomorph._get_valid_transitions`, `Rhizomorph._total_outgoing_weight`,
# `Rhizomorph._edge_transition_weight`), and the arithmetic is evaluated in the
# same associativity: `next = (state_score + log(edge_w/total_out)) + emission`.
# Best-state selection and path reconstruction reuse the scalar helpers verbatim
# (`_best_correction_state`, `_best_correction_target_state`,
# `_reconstruct_correction_path`), so tie-breaks are identical. The
# `batched_equivalence_oracle` below is the acceptance gate AND the regression
# guard for Phases A/B: it asserts the batched decoder returns byte-identical
# corrected paths + scores to `correct_observations` per read.
#
# PoC scope limits (deliberate, per ADR §Prototype scope):
#   * single-strand only (`strand_mode == :singlestrand`) — no strand branching;
#   * `beam_width == typemax(Int)` (exact decode) — no beam pruning is applied,
#     so equivalence is asserted against the exact scalar path;
#   * MetaGraphsNext emission-scored observations only (the corrector's residual
#     hard-window decode path).
#
# This function is NOT wired into the production corrector. `correct_observations`
# / `_correct_metagraphs_next_observations` are untouched, so `strategy=:exhaustive`
# (and every existing caller) stays byte-identical by construction. Wiring the
# batched path in behind `strategy=:scalable` is a follow-on once Phase A (SIMD
# over the read axis) and Phase B (KernelAbstractions/CUDA warp-shuffle wavefront)
# are built on this frontier layout.

"""
$(DocStringExtensions.TYPEDSIGNATURES)

One precomputed out-edge of the shared, immutable decode graph, in the batched
structure-of-arrays layout. `log_transition_prob = log(edge_weight / total_out)`
is graph-only (read-independent) and is therefore hoisted out of the per-read
inner loop; `next_vertex` is retained so the per-read emission term can be scored
against the destination label exactly as the scalar decoder does.
"""
struct BatchedOutEdge{L}
    next_state_id::Int
    log_transition_prob::Float64
    next_vertex::L
    next_strand::Rhizomorph.StrandOrientation
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Immutable, batch-shared decode graph for the array-frontier PoC. `states[i]` is
the `(vertex, strand)` for contiguous `state_id == i`; `state_id[state]` is the
inverse map; `out_edges[i]` is the precomputed out-edge list for `state_id == i`
(empty when the state has no usable outgoing mass). Built ONCE per batch by
[`build_batched_decode_graph`](@ref) and read by every read/lane.
"""
struct BatchedDecodeGraph{L}
    states::Vector{Tuple{L, Rhizomorph.StrandOrientation}}
    state_id::Dict{Tuple{L, Rhizomorph.StrandOrientation}, Int}
    out_edges::Vector{Vector{BatchedOutEdge{L}}}
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build the contiguous `state_id <- (vertex, strand)` index and the precomputed
out-edge structure-of-arrays for `graph`. State ids enumerate every graph label
crossed with both strand orientations (`Forward`, `Reverse`) in label order, so
the map is stable and dense. Every topology/weight query uses the same Rhizomorph
primitives the scalar decoder calls, and `log_transition_prob` is computed with
the identical `log(edge_w / total_out)` expression, so the precomputation is
byte-identical to the scalar per-step arithmetic.
"""
function build_batched_decode_graph(
        graph::MetaGraphsNext.MetaGraph,
        label_type::Type{L}
)::BatchedDecodeGraph{L} where {L}
    labels = collect(MetaGraphsNext.labels(graph))
    strands = (Rhizomorph.Forward, Rhizomorph.Reverse)

    states = Vector{Tuple{L, Rhizomorph.StrandOrientation}}()
    state_id = Dict{Tuple{L, Rhizomorph.StrandOrientation}, Int}()
    for vertex in labels
        for strand in strands
            state = (convert(L, vertex), strand)
            push!(states, state)
            state_id[state] = length(states)
        end
    end

    out_edges = [Vector{BatchedOutEdge{L}}() for _ in 1:length(states)]
    for sid in 1:length(states)
        vertex, strand = states[sid]
        transitions = Rhizomorph._get_valid_transitions(graph, vertex, strand)
        isempty(transitions) && continue
        total_out = Rhizomorph._total_outgoing_weight(graph, vertex, strand)
        (isfinite(total_out) && total_out > 0.0) || continue
        for transition in transitions
            next_vertex = convert(L, transition[:target_vertex])
            next_strand = Rhizomorph._normalize_strand(transition[:target_strand])
            edge_w = Rhizomorph._edge_transition_weight(transition[:edge_data])
            edge_w <= 0.0 && continue
            transition_prob = edge_w / total_out
            log_transition_prob = log(transition_prob)
            next_sid = state_id[(next_vertex, next_strand)]
            push!(
                out_edges[sid],
                BatchedOutEdge{L}(next_sid, log_transition_prob, next_vertex, next_strand)
            )
        end
    end

    return BatchedDecodeGraph{L}(states, state_id, out_edges)
end

# Fresh per-read diagnostics dict, mirroring the scalar decoder's initial
# `diagnostics` (src/viterbi-next.jl:977) field-for-field, but tagging the
# algorithm as the batched variant. The equivalence oracle gates on paths +
# scores + target status, not on the mechanical expand/generate counters, but the
# counters are populated per-read (per-lane) so a lane never mutates a shared
# scalar in the inner kernel (ADR §Decision item 4).
function _batched_new_diagnostics(
        alphabet::Symbol,
        strand_mode::Symbol,
        observation_length::Int,
        target_vertex,
        config::ViterbiCorrectionConfig,
        quality_graph::MetaGraphsNext.MetaGraph,
        transition_scoring::Symbol
)::Dict{Symbol, Any}
    return Dict{Symbol, Any}(
        :algorithm => :viterbi_emission_correct_observation_batched,
        :exact => true,
        :alphabet => alphabet,
        :strand_mode => strand_mode,
        :reverse_complement_support => _viterbi_supports_reverse_complement(alphabet),
        :max_steps => observation_length - 1,
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
end

# Materialize the reached-states score Dict for read `r`'s current frontier row.
# Only used off the hot relaxation path (start init + per-depth best selection),
# so the exact scalar tie-break helpers (`_best_correction_state`,
# `_best_correction_target_state`) can be reused verbatim.
function _batched_reached_scores(
        frontier::Matrix{Float64},
        r::Int,
        states::Vector{S}
)::Dict{S, Float64} where {S}
    reached = Dict{S, Float64}()
    @inbounds for sid in 1:size(frontier, 2)
        score = frontier[r, sid]
        isfinite(score) && (reached[states[sid]] = score)
    end
    return reached
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Decode one length-bin: `bin_observations` are reads that all share the same
unit-count `L`, decoded in lock-step against the shared `decode_graph`. The loop
nesting is depth-outer / read-inner — the ADR's `for depth { for read { relax } }`
— with the frontier held as a dense `n x n_states` matrix. Returns one
`Rhizomorph.ViterbiDecodingResult` per read, byte-identical to what
`_viterbi_correct_observation` returns for that read.

`graph` is the traversal (weighted) graph; `quality_graph` is the emission graph
(the two differ exactly as in `_correct_metagraphs_next_observations`).
"""
function _batched_decode_bin(
        graph::MetaGraphsNext.MetaGraph,
        quality_graph::MetaGraphsNext.MetaGraph,
        config::ViterbiCorrectionConfig,
        alphabet::Symbol,
        strand_mode::Symbol,
        transition_scoring::Symbol,
        labels::Vector{L},
        decode_graph::BatchedDecodeGraph{L},
        bin_observations::Vector
)::Vector{Rhizomorph.ViterbiDecodingResult} where {L}
    states = decode_graph.states
    state_id = decode_graph.state_id
    out_edges = decode_graph.out_edges
    n = length(bin_observations)
    n_states = length(states)
    observation_length = length(bin_observations[1])
    StateT = Tuple{L, Rhizomorph.StrandOrientation}

    # Dense frontier (ADR item 1): current + next score rows, `-Inf` = unreached.
    current_scores = fill(-Inf, n, n_states)
    next_scores = fill(-Inf, n, n_states)
    next_predecessor = zeros(Int, n, n_states)

    finished = falses(n)
    best_state = Vector{Union{Nothing, StateT}}(nothing, n)
    best_score = fill(-Inf, n)
    best_depth = zeros(Int, n)
    targets = Vector{Any}(undef, n)
    predecessors_by_depth = [Vector{Dict{StateT, StateT}}() for _ in 1:n]
    diagnostics = Vector{Dict{Symbol, Any}}(undef, n)
    results = Vector{Rhizomorph.ViterbiDecodingResult}(undef, n)

    # ---- depth 0: per-read start-emission init (mirrors scalar :1001-1041) ----
    for r in 1:n
        observation = bin_observations[r]
        start_observed = first(observation)
        target_vertex = _emission_target_vertex(graph, observation, config, alphabet, strand_mode)
        targets[r] = target_vertex
        diag = _batched_new_diagnostics(
            alphabet, strand_mode, observation_length, target_vertex,
            config, quality_graph, transition_scoring
        )
        diagnostics[r] = diag

        start_candidates = _viterbi_start_candidates(
            labels, _viterbi_label_unit(start_observed), alphabet, strand_mode)
        for vertex in start_candidates
            for strand in _viterbi_start_strands(graph, vertex, strand_mode, config.start_strand)
                sid = state_id[(convert(L, vertex), strand)]
                score = _call_viterbi_state_emission_logp(
                    quality_graph, config, start_observed, vertex, alphabet, strand_mode)
                if isfinite(score) && score > current_scores[r, sid]
                    current_scores[r, sid] = score
                end
            end
        end

        reached = _batched_reached_scores(current_scores, r, states)
        if isempty(reached)
            diag[:reason] = :no_finite_start_emission
            results[r] = Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diag)
            finished[r] = true
            continue
        end
        diag[:retained_states] = length(reached)
        diag[:cumulative_retained_states] = length(reached)
        diag[:max_retained_states] = length(reached)

        bstate, bscore = _best_correction_state(reached)
        best_state[r] = bstate
        best_score[r] = bscore
        best_depth[r] = 0
        if target_vertex !== nothing
            if observation_length == 1
                tstate, tscore = _best_correction_target_state(reached, target_vertex)
                if tstate !== nothing
                    best_state[r] = tstate
                    best_score[r] = tscore
                    diag[:reached_target] = true
                else
                    best_score[r] = -Inf
                end
            else
                best_score[r] = -Inf
            end
        end
    end

    # ---- depth-outer, read-inner relaxation (ADR item 2) ----
    for depth in 1:(observation_length - 1)
        # Reset the next-frontier rows for reads still advancing.
        for r in 1:n
            finished[r] && continue
            @inbounds for sid in 1:n_states
                next_scores[r, sid] = -Inf
                next_predecessor[r, sid] = 0
            end
        end

        # Read axis is the inner, data-parallel axis: every active read relaxes
        # its dense row against the shared, read-only out-edge SoA.
        for r in 1:n
            finished[r] && continue
            observed_unit = bin_observations[r][depth + 1]
            diag = diagnostics[r]
            @inbounds for sid in 1:n_states
                state_score = current_scores[r, sid]
                isfinite(state_score) || continue
                diag[:expanded_states] += 1
                edges = out_edges[sid]
                isempty(edges) && continue
                for edge in edges
                    emission_score = _call_viterbi_state_emission_logp(
                        quality_graph, config, observed_unit, edge.next_vertex,
                        alphabet, strand_mode)
                    candidate = state_score + edge.log_transition_prob + emission_score
                    if !isfinite(candidate)
                        diag[:skipped_transitions] += 1
                        continue
                    end
                    diag[:generated_states] += 1
                    if candidate > next_scores[r, edge.next_state_id]
                        next_scores[r, edge.next_state_id] = candidate
                        next_predecessor[r, edge.next_state_id] = sid
                    end
                end
            end
        end

        # Per-read commit: materialize predecessors, swap frontier, update best.
        for r in 1:n
            finished[r] && continue
            reached_scores = Dict{StateT, Float64}()
            reached_predecessors = Dict{StateT, StateT}()
            @inbounds for sid in 1:n_states
                score = next_scores[r, sid]
                isfinite(score) || continue
                state = states[sid]
                reached_scores[state] = score
                reached_predecessors[state] = states[next_predecessor[r, sid]]
            end

            if isempty(reached_scores)
                finished[r] = true
                continue
            end

            push!(predecessors_by_depth[r], reached_predecessors)
            @inbounds for sid in 1:n_states
                current_scores[r, sid] = next_scores[r, sid]
            end

            diag = diagnostics[r]
            retained = length(reached_scores)
            diag[:retained_states] = retained
            diag[:cumulative_retained_states] += retained
            diag[:max_retained_states] = max(diag[:max_retained_states], retained)
            diag[:completed_steps] = depth

            target_vertex = targets[r]
            if target_vertex === nothing
                bstate, bscore = _best_correction_state(reached_scores)
                best_state[r] = bstate
                best_score[r] = bscore
                best_depth[r] = depth
            else
                tstate, tscore = _best_correction_target_state(reached_scores, target_vertex)
                if tstate !== nothing
                    best_state[r] = tstate
                    best_score[r] = tscore
                    best_depth[r] = depth
                    diag[:reached_target] = true
                end
            end
        end
    end

    # ---- finalize per read (mirrors scalar :1143-1150) ----
    for r in 1:n
        isassigned(results, r) && continue   # already resolved (no finite start)
        diag = diagnostics[r]
        target_vertex = targets[r]
        if target_vertex !== nothing && !isfinite(best_score[r])
            diag[:reason] = :target_unreachable
            results[r] = Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diag)
            continue
        end
        end_state = something(best_state[r])
        path = _reconstruct_correction_path(
            graph, end_state, best_depth[r], predecessors_by_depth[r])
        diag[:path_length] = length(path.steps)
        results[r] = Rhizomorph.ViterbiDecodingResult(path, best_score[r], diag)
    end

    return results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Batched array-frontier corrector — CPU proof-of-concept (bead td-qoo3).

Decodes a batch of `observations` against a shared, immutable Rhizomorph graph
using the ADR reformulation: a dense array frontier, a depth-outer / read-inner
loop, length-binned batches, and a graph structure-of-arrays built once per batch.
Returns a [`ViterbiCorrectionResult`](@ref) with the same shape as
[`correct_observations`](@ref), so per-read corrected paths and scores can be
compared directly against the scalar decoder (see
[`batched_equivalence_oracle`](@ref)).

PoC constraints (throws `ArgumentError` otherwise): single-strand graphs only
(`strand_mode == :singlestrand`), exact decode (`beam_width == typemax(Int)`), and
non-empty emission-scored observations against a `MetaGraphsNext.MetaGraph`. This
is the CPU foundation the future SIMD (Phase A) and GPU (Phase B) kernels build
on; it does not itself add any SIMD/GPU dependency and is not wired into the
production corrector.
"""
function batched_correct_observations(
        graph::MetaGraphsNext.MetaGraph,
        observations = _default_correction_observations(graph);
        config::ViterbiCorrectionConfig = ViterbiCorrectionConfig()
)::ViterbiCorrectionResult
    if config.beam_width != typemax(Int)
        throw(ArgumentError(
            "batched_correct_observations is an exact-decode PoC and requires " *
            "beam_width == typemax(Int) (got $(config.beam_width)); beam pruning " *
            "is a follow-on once Phase A/B are built on the array frontier."))
    end
    if observations === nothing || isempty(observations)
        throw(ArgumentError("no observations to correct"))
    end
    for observation in observations
        _uses_emission_scored_observation(observation) || throw(ArgumentError(
            "batched PoC requires non-empty emission-scored (AbstractVector) " *
            "observations; got $(typeof(observation))"))
    end

    alphabet = _resolve_viterbi_alphabet(graph, observations, config.alphabet)
    strand_mode = _resolve_viterbi_strand_mode(graph, alphabet, config.strand_mode)
    if strand_mode != :singlestrand
        throw(ArgumentError(
            "batched PoC supports single-strand decoding only (resolved " *
            "strand_mode=$strand_mode); multi-strand batching is deferred per the ADR."))
    end

    transition_edge_weight = _viterbi_transition_edge_weight(graph, config.edge_weight)
    transition_scoring = _viterbi_transition_scoring(graph, transition_edge_weight)
    weighted = if _correction_edge_data_type(graph) <: Rhizomorph.StrandWeightedEdgeData
        graph
    else
        Rhizomorph.weighted_graph_from_rhizomorph(graph; edge_weight = transition_edge_weight)
    end

    labels = collect(MetaGraphsNext.labels(weighted))
    if isempty(labels)
        throw(ArgumentError("cannot correct observations against an empty graph"))
    end
    label_type = eltype(labels)

    # Shared, immutable graph built ONCE per batch (ADR item 4).
    decode_graph = build_batched_decode_graph(weighted, label_type)

    # Length-binning (ADR item 3): group observation indices by unit-count.
    bins = Dict{Int, Vector{Int}}()
    for (index, observation) in enumerate(observations)
        push!(get!(bins, length(observation), Int[]), index)
    end

    paths = Vector{Rhizomorph.ViterbiDecodingResult}(undef, length(observations))
    for (_, indices) in bins
        bin_observations = [observations[index] for index in indices]
        bin_results = _batched_decode_bin(
            weighted, graph, config, alphabet, strand_mode, transition_scoring,
            labels, decode_graph, bin_observations)
        for (position, index) in enumerate(indices)
            paths[index] = bin_results[position]
        end
    end

    corrected = Any[_decoded_path_labels(path_result) for path_result in paths]
    diagnostics = Dict{Symbol, Any}(
        :interface => :metagraphs_next_batched,
        :vertex_data_type => _correction_vertex_data_type(graph),
        :algorithm => :rhizomorph_viterbi_decode_next_batched,
        :observation_count => length(observations),
        :emission_callback => nameof(config.emission_logp),
        :alphabet => alphabet,
        :emission_model => _viterbi_graph_has_quality(graph) ?
                           :quality_aware : :alphabet_parameterized,
        :transition_model => transition_scoring,
        :strand_mode => strand_mode,
        :reverse_complement_support => _viterbi_supports_reverse_complement(alphabet),
        :length_bins => sort(collect(keys(bins))),
        :batch_size => length(observations),
        :state_space_size => length(decode_graph.states)
    )
    return ViterbiCorrectionResult(corrected, paths, diagnostics)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Equivalence oracle (the PoC acceptance gate). Runs both the scalar
[`correct_observations`](@ref) and the batched
[`batched_correct_observations`](@ref) on the same graph + observations and checks
that, per read, the corrected label paths and log-probability scores are
byte-identical (`===` on decoded labels, bit-identical `Float64` scores via
`===`), and that path step vertices/strands and `reached_target` status agree.

Returns a named tuple `(; passes, per_read, scalar, batched)` where `per_read` is
a vector of per-read booleans and `passes` is their conjunction. This is both the
acceptance criterion for the reformulation and the regression guard for the
future SIMD (Phase A) and GPU (Phase B) kernels, which must keep passing it.
"""
function batched_equivalence_oracle(
        graph::MetaGraphsNext.MetaGraph,
        observations = _default_correction_observations(graph);
        config::ViterbiCorrectionConfig = ViterbiCorrectionConfig()
)
    scalar = correct_observations(graph, observations; config = config)
    batched = batched_correct_observations(graph, observations; config = config)

    per_read = Bool[]
    for index in 1:length(scalar.paths)
        scalar_path = scalar.paths[index]
        batched_path = batched.paths[index]

        labels_match = _decoded_path_labels(scalar_path) == _decoded_path_labels(batched_path)
        # `===` compares Float64 bit patterns (so -Inf === -Inf and NaN handling
        # are explicit), which is the "byte-identical score" gate the ADR asks for.
        score_match = scalar_path.score === batched_path.score
        target_match = get(scalar_path.diagnostics, :reached_target, nothing) ==
                       get(batched_path.diagnostics, :reached_target, nothing)

        steps_match = if scalar_path.path === nothing || batched_path.path === nothing
            scalar_path.path === batched_path.path
        else
            scalar_steps = something(scalar_path.path).steps
            batched_steps = something(batched_path.path).steps
            length(scalar_steps) == length(batched_steps) &&
                all(
                    scalar_steps[step_index].vertex_label ==
                        batched_steps[step_index].vertex_label &&
                    scalar_steps[step_index].strand == batched_steps[step_index].strand
                    for step_index in 1:length(scalar_steps)
                )
        end

        push!(per_read, labels_match && score_match && target_match && steps_match)
    end

    return (; passes = all(per_read), per_read = per_read, scalar = scalar, batched = batched)
end
