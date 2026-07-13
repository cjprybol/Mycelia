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
        error_rate::Float64 = 0.01,
        count_length_penalty::Bool = true
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
    # Compare by CHARACTER, not byte. `length`/`getindex` on a String count
    # characters vs. address bytes respectively, so `observed[index]` throws a
    # `StringIndexError` on any multibyte UTF-8 unit (e.g. the SentencePiece
    # U+2581 word-boundary marker, accented letters, emoji). Materializing to
    # `Vector{Char}` makes positional indexing character-aware for arbitrary
    # Unicode tokens while keeping the O(n) scan.
    observed_chars = collect(observed)
    expected_chars = collect(expected)
    shared = min(length(observed_chars), length(expected_chars))

    logp = 0.0
    for index in 1:shared
        position_error_rate = _viterbi_position_error_rate(
            quality_scores,
            index,
            error_rate
        )
        if observed_chars[index] == expected_chars[index]
            logp += log1p(-position_error_rate)
        else
            logp += log(position_error_rate / (alphabet_size - 1))
        end
    end

    # Flat trailing length penalty for a length mismatch between observed and
    # expected. This is a FRAMESHIFT charge, not an indel model — with the
    # indel-aware pair-HMM active, insertion/deletion MOVES score the length change
    # explicitly, so keeping this term would DOUBLE-COUNT it and bias the MAP path.
    # The indel kernel passes `count_length_penalty=false` to drop it; the
    # substitution decoder keeps it (default true), where it is score-neutral on the
    # equal-length correct-back-to-reference fixtures anyway.
    if count_length_penalty
        for index in (shared + 1):max(length(observed_chars), length(expected_chars))
            position_error_rate = _viterbi_position_error_rate(
                quality_scores,
                index,
                error_rate
            )
            logp += log(position_error_rate / alphabet_size)
        end
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

Two candidate-generation bounds (both default to exact / no-op; opted into by the
:scalable corrector alongside the size-aware `beam_width`) keep the per-depth
frontier expansion O(1) in graph size instead of growing toward the width beam as
the graph densifies (td-plqi):

  * `max_successors_per_state` — per expanded state, only the top-B outgoing
    transitions by edge weight are materialized before emission scoring, so
    generation per state is O(B) not O(out-degree). On a k-mer de Bruijn graph the
    structural out-degree is ≤ 4 (the four next bases), so any `B ≥ 4` is a strict
    no-op here; it is a robustness guard for pathological high-branching inputs
    (collapsed-repeat multigraphs, non-DNA alphabets). Default `typemax(Int)`.

  * `beam_score_margin` — Δ log-probability threshold ("histogram" beam pruning)
    with an EMISSION exemption. After the width beam, a state is dropped only when
    it is BOTH >Δ below the depth's best FULL score AND >Δ below the depth's best
    cumulative EMISSION (read-consistency, excluding the coverage/transition term).
    This is the term that actually bounds generation at scale: on a dense
    intermediate-k graph the width-256 frontier is ~99% read-INCONSISTENT states
    (low on both axes) that can never win the ML path yet still generate +
    emission-score successors each depth; pruning them holds the generating
    frontier to a few states independent of genome size. The emission clause
    protects variation: a real-but-rare minor allele (skewed-coverage / viral
    quasispecies) has GOOD emission but a coverage-driven transition penalty, so
    the emission exemption keeps it even when its full score trails the dominant
    haplotype — the margin prunes WRONG paths (bad emission), never merely RARE
    ones. It does not weaken error correction (an uncorrected error path has high
    emission → exempt → still available for the full-score ML choice). Default
    `Inf` (no threshold = exact). Only engages where the width beam is already
    finite (approximate), so exact-ML reads (`beam_width == typemax`) stay
    byte-identical.
"""
struct ViterbiCorrectionConfig{F <: Function}
    error_rate::Float64
    verbosity::String
    emission_logp::F
    alphabet::Symbol
    strand_mode::Symbol
    max_steps::Union{Nothing, Int}
    target_vertex::Any
    start_strand::Rhizomorph.StrandOrientation
    edge_weight::Function
    beam_width::Int
    max_successors_per_state::Int
    beam_score_margin::Float64
    # OPT-IN telemetry (td-4osf / Tier-2 calibration): when true, the decoder
    # records a per-position best-vs-2nd-best log-prob GAP (the Viterbi margin at
    # each depth) into its diagnostics. Default false => the decode is byte-
    # identical AND perf-neutral (no top-2 scan on the hot path).
    record_position_gaps::Bool
    # Indel-aware pair-HMM extension (nanopore correction). All fields DEFAULT to
    # the substitution-collapse values, so the constructor default reproduces the
    # substitution-only decoder byte-for-byte: `indel_moves=false` short-circuits
    # before the gap kernel is ever entered, and every gap mass is 0 (⇒ gap
    # transitions score `log(0) = -Inf` ⇒ dropped by the `isfinite` guard even if
    # the kernel is entered with these zero fractions). `insertion_fraction` and
    # `deletion_fraction` partition the per-base error rate into gap-open masses
    # (`δ_I = error_rate·f_ins`, `δ_D = error_rate·f_del`); the remainder is the
    # substitution mass that rides inside the existing emission term. The
    # `*_extend_probability` fields are the affine gap-EXTEND probabilities
    # (`γ_I`, `γ_D`) — one-time open then cheaper extend, because nanopore indels
    # cluster in homopolymer runs (geometric run lengths). `deletion_max_run`
    # (`D_max`) bounds consecutive deletions (the bounded Bellman-Ford relaxation
    # replacing the legacy O(V³) Floyd-Warshall; also guards graph-cycle
    # no-progress loops). `max_insertion_run` bounds consecutive insertions.
    # `band_width` is the adaptive diagonal-band half-width on the net gap
    # (graph-steps − read-index); `nothing` = unbounded (the exact/oracle setting).
    indel_moves::Bool
    insertion_fraction::Float64
    deletion_fraction::Float64
    insertion_extend_probability::Float64
    deletion_extend_probability::Float64
    deletion_max_run::Int
    max_insertion_run::Int
    band_width::Union{Nothing, Int}
    insertion_emission_logp::Union{Nothing, Function}

    function ViterbiCorrectionConfig{F}(;
            error_rate::Float64 = 0.01,
            verbosity::String = "dataset",
            emission_logp::F,
            alphabet::Symbol = :auto,
            strand_mode::Symbol = :auto,
            max_steps::Union{Nothing, Int} = nothing,
            target_vertex = nothing,
            start_strand::Rhizomorph.StrandOrientation = Rhizomorph.Forward,
            edge_weight::Function = Rhizomorph.edge_data_weight,
            beam_width::Int = typemax(Int),
            max_successors_per_state::Int = typemax(Int),
            beam_score_margin::Float64 = Inf,
            record_position_gaps::Bool = false,
            indel_moves::Bool = false,
            insertion_fraction::Float64 = 0.0,
            deletion_fraction::Float64 = 0.0,
            insertion_extend_probability::Float64 = 0.0,
            deletion_extend_probability::Float64 = 0.0,
            deletion_max_run::Int = 0,
            max_insertion_run::Int = 0,
            band_width::Union{Nothing, Int} = nothing,
            insertion_emission_logp::Union{Nothing, Function} = nothing
    ) where {F <: Function}
        if error_rate <= 0.0 || error_rate >= 0.5
            throw(ArgumentError("error_rate must be in (0, 0.5), got $error_rate"))
        end
        if insertion_fraction < 0.0 || deletion_fraction < 0.0 ||
           insertion_fraction + deletion_fraction >= 1.0
            throw(ArgumentError(
                "insertion_fraction/deletion_fraction must be non-negative and sum " *
                "to < 1 (the remainder is substitution mass), got " *
                "$insertion_fraction / $deletion_fraction"))
        end
        for (name, value) in (
            (:insertion_extend_probability, insertion_extend_probability),
            (:deletion_extend_probability, deletion_extend_probability)
        )
            if value < 0.0 || value >= 1.0
                throw(ArgumentError("$name must be in [0, 1), got $value"))
            end
        end
        if deletion_max_run < 0
            throw(ArgumentError("deletion_max_run must be non-negative, got $deletion_max_run"))
        end
        if max_insertion_run < 0
            throw(ArgumentError("max_insertion_run must be non-negative, got $max_insertion_run"))
        end
        if band_width !== nothing && band_width < 0
            throw(ArgumentError("band_width must be non-negative or nothing, got $band_width"))
        end
        if !(verbosity in ("debug", "reads", "dataset", "silent"))
            throw(ArgumentError("unsupported verbosity: $verbosity"))
        end
        if max_steps !== nothing && max_steps < 0
            throw(ArgumentError("max_steps must be non-negative, got $max_steps"))
        end
        if beam_width <= 0
            throw(ArgumentError("beam_width must be positive, got $beam_width"))
        end
        if max_successors_per_state <= 0
            throw(ArgumentError(
                "max_successors_per_state must be positive, got $max_successors_per_state"))
        end
        if isnan(beam_score_margin) || beam_score_margin <= 0.0
            throw(ArgumentError(
                "beam_score_margin must be a positive number or Inf, got $beam_score_margin"))
        end
        # Silent-no-op guard (PR #407 review, silent-failure I2): `indel_moves=true`
        # with NO indel capacity — either both gap masses zero
        # (`insertion_fraction == deletion_fraction == 0`) or both run caps zero
        # (`deletion_max_run == max_insertion_run == 0`) — sends every gap
        # transition to `log(0) = -Inf`, so the "indel-aware" decode silently
        # collapses to substitution-only. This is a legitimate configuration (it is
        # exactly how the collapse/oracle test proves byte-identity), so we WARN
        # rather than error, but we surface it so an operator does not believe they
        # enabled indel correction when they did not.
        if indel_moves && (
            (insertion_fraction == 0.0 && deletion_fraction == 0.0) ||
            (deletion_max_run == 0 && max_insertion_run == 0)
        )
            @warn "indel_moves=true but no indel capacity is configured " *
                  "(insertion_fraction/deletion_fraction both 0, or " *
                  "deletion_max_run/max_insertion_run both 0); the decode is " *
                  "effectively substitution-only."
        end
        return new{F}(
            error_rate,
            verbosity,
            emission_logp,
            alphabet,
            strand_mode,
            max_steps,
            target_vertex,
            start_strand,
            edge_weight,
            beam_width,
            max_successors_per_state,
            beam_score_margin,
            record_position_gaps,
            indel_moves,
            insertion_fraction,
            deletion_fraction,
            insertion_extend_probability,
            deletion_extend_probability,
            deletion_max_run,
            max_insertion_run,
            band_width,
            insertion_emission_logp
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
        target_vertex = nothing,
        start_strand::Rhizomorph.StrandOrientation = Rhizomorph.Forward,
        edge_weight::Function = Rhizomorph.edge_data_weight,
        beam_width::Int = typemax(Int),
        max_successors_per_state::Int = typemax(Int),
        beam_score_margin::Float64 = Inf,
        record_position_gaps::Bool = false,
        indel_moves::Bool = false,
        insertion_fraction::Float64 = 0.0,
        deletion_fraction::Float64 = 0.0,
        insertion_extend_probability::Float64 = 0.0,
        deletion_extend_probability::Float64 = 0.0,
        deletion_max_run::Int = 0,
        max_insertion_run::Int = 0,
        band_width::Union{Nothing, Int} = nothing,
        insertion_emission_logp::Union{Nothing, Function} = nothing
)::ViterbiCorrectionConfig{F} where {F <: Function}
    return ViterbiCorrectionConfig{F}(
        error_rate = error_rate,
        verbosity = verbosity,
        emission_logp = emission_logp,
        alphabet = alphabet,
        strand_mode = strand_mode,
        max_steps = max_steps,
        target_vertex = target_vertex,
        start_strand = start_strand,
        edge_weight = edge_weight,
        beam_width = beam_width,
        max_successors_per_state = max_successors_per_state,
        beam_score_margin = beam_score_margin,
        record_position_gaps = record_position_gaps,
        indel_moves = indel_moves,
        insertion_fraction = insertion_fraction,
        deletion_fraction = deletion_fraction,
        insertion_extend_probability = insertion_extend_probability,
        deletion_extend_probability = deletion_extend_probability,
        deletion_max_run = deletion_max_run,
        max_insertion_run = max_insertion_run,
        band_width = band_width,
        insertion_emission_logp = insertion_emission_logp
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
An observation unit that carries the READ's own per-base Phred quality window
alongside the observed k-mer.

This is the central seam that restores the graph-as-HMM emission model
(`docs/design/2026-07-06-graph-as-hmm-correction-redesign.md`). Bare k-mer
observations force the emission to fall back to the graph vertex's
population-average quality (or a uniform `error_rate`), so `P(observed | hidden)`
loses the *observed* — it measures the graph's confidence in a k-mer rather than
the read's confidence in its base. Wrapping each k-mer with the read's per-base
Phred (`quality_scores[idx-k+1:idx]` for the k-mer ending at read index `idx`)
makes `_viterbi_direct_quality_scores` return the read's real quality, so a
low-quality read base yields a higher tolerance for correcting toward the graph,
and a high-quality base resists it.

!!! warning "RAW Phred, not ASCII Phred+33"
    `quality_scores` stores RAW Phred values (integer quality, as returned by
    `FASTX.quality_scores` — e.g. `0x28` == Q40), NOT ASCII-encoded Phred+33
    bytes (where Q40 is the byte `'I'` == `0x49` == 73). This type deliberately
    bypasses the generic ASCII-vs-Phred decode heuristic used elsewhere and
    consumes the bytes verbatim (see `_viterbi_direct_quality_scores`). A caller
    that passes ASCII Phred+33 bytes here would inflate every quality by ~33
    (Q40 read as Q73), collapsing the emission model. Always construct with the
    integer quality window (`FASTX.quality_scores(read)`), never the ASCII
    quality string.

The inner constructor enforces the one hard invariant — a non-empty quality
window. An empty window makes `_viterbi_direct_quality_scores` return `nothing`,
which silently drops the read's own quality and falls the emission back to the
graph's population-average quality; guarding here surfaces that as an error at
the construction site instead of a silent degradation deep in the decoder. (No
ASCII-detection heuristic is added on purpose: raw-Phred and ASCII-Phred ranges
overlap, so any such heuristic would misfire — the doc invariant above plus this
non-empty guard is the deliberate, non-fragile contract.)
"""
struct QualityObservation{KmerT}
    kmer::KmerT
    quality_scores::Vector{UInt8}

    function QualityObservation(kmer::KmerT,
            quality_scores::AbstractVector{UInt8}) where {KmerT}
        isempty(quality_scores) && throw(ArgumentError(
            "QualityObservation requires a non-empty RAW-Phred quality window; an " *
            "empty window would silently fall the emission back to the graph's " *
            "population-average quality (see the RAW-Phred invariant in the docstring)."))
        return new{KmerT}(kmer, quality_scores)
    end
end

# The observed-unit accessors used by the emission chain. Bare-k-mer behavior is
# preserved by delegating to the wrapped k-mer, so wrapping only ADDS the read's
# per-base quality without changing string/alphabet resolution.
_viterbi_unit_string(unit::QualityObservation)::String = _viterbi_unit_string(unit.kmer)

function _infer_viterbi_alphabet(unit::QualityObservation)::Union{Nothing, Symbol}
    return _infer_viterbi_alphabet(unit.kmer)
end

# QualityObservation stores RAW Phred; return it verbatim (no ASCII decode), which
# is both correct and unambiguous for high-quality reads (Phred >= 33 would
# otherwise be misread as ASCII-offset by the generic `Any` heuristic).
function _viterbi_direct_quality_scores(
        unit::QualityObservation
)::Union{Nothing, Vector{Float64}}
    return isempty(unit.quality_scores) ? nothing : Float64.(unit.quality_scores)
end

# For graph-label comparisons (start-candidate matching, target resolution) the
# emission-carrying wrapper must be transparent: unwrap to the underlying k-mer so
# the decode stays as targeted/efficient as it was for bare k-mers. Emission
# scoring still receives the full `QualityObservation`.
_viterbi_label_unit(unit) = unit
_viterbi_label_unit(unit::QualityObservation) = unit.kmer

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
        config::ViterbiCorrectionConfig = ViterbiCorrectionConfig(),
        weighted_graph = nothing
)
    vertex_data_type = _correction_vertex_data_type(graph)
    return _correct_observations(
        vertex_data_type, graph, observations; config = config,
        weighted_graph = weighted_graph)
end

function _correct_observations(
        ::Type{<:Any},
        graph::MetaGraphs.MetaDiGraph,
        observations;
        config::ViterbiCorrectionConfig,
        weighted_graph = nothing
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
        config::ViterbiCorrectionConfig,
        weighted_graph = nothing
)::ViterbiCorrectionResult
    return _correct_metagraphs_next_observations(
        graph, observations; config = config, weighted_graph = weighted_graph)
end

function _correct_observations(
        ::Type{Any},
        graph::MetaGraphsNext.MetaGraph,
        observations;
        config::ViterbiCorrectionConfig,
        weighted_graph = nothing
)::ViterbiCorrectionResult
    return _correct_metagraphs_next_observations(
        graph, observations; config = config, weighted_graph = weighted_graph)
end

# O(1) label type from the MetaGraph type parameters (used only to type the
# precomputed outgoing-weight table's keys — does NOT touch the decode's own
# state dicts, so it cannot affect decode determinism).
function _viterbi_graph_label_type(
        graph::MetaGraphsNext.MetaGraph{
        CODE, GRAPH, LABEL}
)::Type where {CODE, GRAPH, LABEL}
    return LABEL
end

# Build the weighted decode graph ONCE for a correction pass so callers that
# decode many reads against the same graph (the :scalable corrector's per-read
# loop) can hoist it OUT of the per-read path (td-y4oj). `weighted_graph_from_
# rhizomorph` is O(V+E); with n_reads ∝ genome and V ∝ genome, rebuilding it per
# read is one of two co-dominant O(genome^2) terms (the other is the per-state
# outgoing-weight scan, fixed in `_total_outgoing_weight`). The weighted graph
# depends ONLY on `graph` and the (constant-across-reads) `config.edge_weight`, so
# a single build is bit-identical in DATA to what each per-read
# `_correct_metagraphs_next_observations` call would have produced (build
# determinism verified). Returns `nothing` for graphs that don't use the
# weighted-decode path (legacy stranded / non-MetaGraph inputs). READ-ONLY during
# decode (correction mutates a separate accumulator, never the graph — verified no
# edge mutation across a pass), so one instance is safely shared across reads and
# threads.
function build_correction_weighted_graph(
        graph::MetaGraphsNext.MetaGraph;
        config::ViterbiCorrectionConfig = ViterbiCorrectionConfig()
)
    if _correction_edge_data_type(graph) <: Rhizomorph.StrandWeightedEdgeData
        return graph
    end
    transition_edge_weight = _viterbi_transition_edge_weight(graph, config.edge_weight)
    return Rhizomorph.weighted_graph_from_rhizomorph(
        graph; edge_weight = transition_edge_weight)
end

function build_correction_weighted_graph(graph; config::ViterbiCorrectionConfig = ViterbiCorrectionConfig())
    nothing
end

function _correct_metagraphs_next_observations(
        graph::MetaGraphsNext.MetaGraph,
        observations;
        config::ViterbiCorrectionConfig,
        weighted_graph = nothing
)::ViterbiCorrectionResult
    alphabet = _resolve_viterbi_alphabet(graph, observations, config.alphabet)
    strand_mode = _resolve_viterbi_strand_mode(graph, alphabet, config.strand_mode)
    transition_edge_weight = _viterbi_transition_edge_weight(graph, config.edge_weight)
    weighted = if weighted_graph !== nothing
        # Hoisted precomputed weighted graph (td-y4oj) — reuse across reads.
        weighted_graph
    elseif _correction_edge_data_type(graph) <: Rhizomorph.StrandWeightedEdgeData
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
            start_vertex, target_vertex,
            max_steps = _decode_observation_bounds(
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
        :emission_model => _viterbi_graph_has_quality(graph) ?
                           :quality_aware : :alphabet_parameterized,
        :transition_model => _viterbi_transition_scoring(graph, transition_edge_weight),
        :strand_mode => strand_mode,
        :reverse_complement_support => _viterbi_supports_reverse_complement(alphabet)
    )
    return ViterbiCorrectionResult(corrected, paths, diagnostics)
end

function _correction_vertex_data_type(
        graph::MetaGraphsNext.MetaGraph{
        CODE, GRAPH, LABEL, VERTEX_DATA}
)::Type where {CODE, GRAPH, LABEL, VERTEX_DATA}
    return VERTEX_DATA
end

function _correction_edge_data_type(
        graph::MetaGraphsNext.MetaGraph{CODE, GRAPH, LABEL, VERTEX_DATA,
        EDGE_DATA}
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
        target_vertex = config.target_vertex === nothing ? last(observation) :
                        config.target_vertex
        max_steps = config.max_steps === nothing ? length(observation) - 1 :
                    config.max_steps
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
        quality = nothing,
        count_length_penalty::Bool = true
)::Float64
    # The `count_length_penalty=false` request (indel kernel dropping the flat
    # frameshift charge) is only meaningful for the built-in emission; a custom
    # callback need not accept the keyword, so route it explicitly only for the
    # default and leave every other caller/callback byte-identical.
    if !count_length_penalty && config.emission_logp === default_viterbi_emission_logp
        return default_viterbi_emission_logp(
            observed_unit,
            node,
            alphabet;
            quality = quality,
            error_rate = config.error_rate,
            count_length_penalty = false
        )
    end
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

# Keep only the `beam_width` highest-scoring states, pruning the score and
# predecessor dicts in lockstep so path reconstruction stays consistent. States
# are ranked by score descending; score ties are broken by `hash(state)` for a
# deterministic (run-reproducible) beam. Called once per depth by
# `_viterbi_correct_observation` when the frontier exceeds the beam.
function _prune_correction_beam(
        scores::Dict{S, Float64},
        predecessors::Dict{S, S},
        beam_width::Int
)::Tuple{Dict{S, Float64}, Dict{S, S}} where {S}
    ordered = sort!(collect(scores); by = kv -> (last(kv), hash(first(kv))), rev = true)
    kept = Dict{S, Float64}()
    kept_predecessors = Dict{S, S}()
    for (state, score) in Iterators.take(ordered, beam_width)
        kept[state] = score
        if haskey(predecessors, state)
            kept_predecessors[state] = predecessors[state]
        end
    end
    return kept, kept_predecessors
end

# Bound successor GENERATION per expanded state (td-plqi): keep only the top-B
# outgoing transitions ranked by edge weight BEFORE any emission is scored, so the
# per-state candidate set is O(B) rather than O(out-degree). Ranking is by
# `_edge_transition_weight` descending; ties are broken by target-vertex string for
# a deterministic (run-reproducible) selection. Returns `transitions` unchanged when
# it already fits within `b` (the common case: a k-mer de Bruijn vertex has ≤ 4
# structural successors, so any `b ≥ 4` is a no-op). Only bites on pathological
# high-branching vertices, where the low-weight tail is the least-probable
# (error/spurious) edges — the true, well-supported branch is high-weight and
# always survives.
function _top_b_transitions(transitions, b::Int)
    length(transitions) <= b && return transitions
    ordered = sort(
        transitions;
        by = t -> (Rhizomorph._edge_transition_weight(t[:edge_data]),
            string(t[:target_vertex])),
        rev = true
    )
    return ordered[1:b]
end

# Score-margin ("histogram") beam pruning with an EMISSION EXEMPTION (td-plqi,
# variation-safety hardening from PR #388 review). Applied AFTER the width beam, in
# lockstep on the score / predecessor / emission dicts so path reconstruction stays
# consistent.
#
# A state is pruned only when it is BOTH:
#   (1) more than `margin` log-prob below the frontier's best FULL score
#       (`state + log(transition_prob) + emission`), AND
#   (2) more than `margin` below the frontier's best cumulative EMISSION
#       (read-consistency: `sum(emission)`, EXCLUDING the coverage/transition term).
#
# Rationale — why the AND-gate rather than a full-score-only threshold: the full
# score carries a `log(transition_prob) = log(cov_edge / cov_total)` term that
# penalizes a path for being RARE, not for being WRONG. A real-but-rare minor
# allele (viral quasispecies / skewed pool) has GOOD emission (reads genuinely
# support it) but a coverage-driven transition penalty; a full-score-only margin
# could prune it for rarity. The emission clause EXEMPTS any read-consistent path
# from pruning regardless of how rare its coverage is, so a supported variant is
# never dropped as "improbable." It does NOT weaken error correction: an
# uncorrected error path has HIGH emission (matches the observed read) so it is
# exempt and stays available, but the corrected path still wins the full-score ML
# selection. Genuine junk (read-INCONSISTENT paths threaded through wrong regions)
# is low on BOTH axes and is still pruned, so the O(1)-frontier speed win holds.
# `margin === Inf` (default) is a no-op, preserving exact decoding.
function _prune_correction_beam_by_margin(
        scores::Dict{S, Float64},
        predecessors::Dict{S, S},
        emissions::Dict{S, Float64},
        best_score::Float64,
        best_emission::Float64,
        margin::Float64
)::Tuple{Dict{S, Float64}, Dict{S, S}, Dict{S, Float64}} where {S}
    (isinf(margin) || isempty(scores)) && return scores, predecessors, emissions
    score_threshold = best_score - margin
    emission_threshold = best_emission - margin
    kept = Dict{S, Float64}()
    kept_predecessors = Dict{S, S}()
    kept_emissions = Dict{S, Float64}()
    for (state, score) in scores
        # Keep unless below BOTH thresholds (emission clause exempts read-consistent
        # paths — including rare-but-supported alleles — from coverage-driven pruning).
        emission = get(emissions, state, -Inf)
        if score < score_threshold && emission < emission_threshold
            continue
        end
        kept[state] = score
        kept_emissions[state] = emission
        if haskey(predecessors, state)
            kept_predecessors[state] = predecessors[state]
        end
    end
    return kept, kept_predecessors, kept_emissions
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

    # Indel-aware pair-HMM branch (nanopore correction). When gap MOVES are enabled
    # the decode is a 2D DP over (read-index, graph-state) with an M/I/D phase
    # layer; when disabled (the DEFAULT) we fall straight through to the untouched
    # substitution lock-step decoder below, so it stays byte-identical.
    if config.indel_moves
        return _viterbi_correct_observation_indel(
            graph,
            observation,
            alphabet;
            config = config,
            strand_mode = strand_mode,
            quality_graph = quality_graph,
            transition_scoring = transition_scoring
        )
    end

    # Per-read decode hot path (td-y4oj): resolve the label type and emptiness in
    # O(1) from the graph's type parameters instead of materializing the whole
    # O(V) label vector on EVERY read, and resolve start/target candidates with
    # O(1) `haskey` lookups. The prior `collect(MetaGraphsNext.labels)` + linear
    # `in labels` membership tests were an O(V)-per-read (⇒ O(genome^2)) cost that
    # is CO-dominant with the outgoing-weight scan. `haskey(graph, x)` matches
    # `x in collect(labels(graph))` exactly (verified: 0 disagreements over 53k
    # k-mer and QualityObservation lookups), so this is a pure performance change.
    if Graphs.nv(graph.graph) == 0
        throw(ArgumentError("cannot correct observations against an empty graph"))
    end

    label_type = _viterbi_graph_label_type(graph)
    start_observed = first(observation)
    target_vertex = _emission_target_vertex(
        graph, observation, config, alphabet, strand_mode)
    start_candidates = _viterbi_start_candidates(
        graph, label_type, _viterbi_label_unit(start_observed), alphabet, strand_mode)

    diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_emission_correct_observation,
        :exact => true,
        :alphabet => alphabet,
        :strand_mode => strand_mode,
        :reverse_complement_support => _viterbi_supports_reverse_complement(alphabet),
        :max_steps => length(observation) - 1,
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
    # Cumulative EMISSION (read-consistency) per state, tracked in lockstep with the
    # full score so the score-margin prune can exempt read-consistent paths from
    # coverage-driven pruning (td-plqi variation-safety). At the start there is no
    # transition term, so the start emission == the start score.
    active_emissions = Dict{Tuple{label_type, Rhizomorph.StrandOrientation}, Float64}()
    for vertex in start_candidates
        for strand in
            _viterbi_start_strands(graph, vertex, strand_mode, config.start_strand)
            state = (vertex, strand)
            score = _call_viterbi_state_emission_logp(
                quality_graph,
                config,
                start_observed,
                vertex,
                alphabet,
                strand_mode
            )
            if isfinite(score) &&
               (!haskey(active_scores, state) || score > active_scores[state])
                active_scores[state] = score
                active_emissions[state] = score
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
            target_state,
            target_score = _best_correction_target_state(active_scores, target_vertex)
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

    # Opt-in Tier-2 telemetry: per-depth Viterbi margin (best - 2nd-best surviving
    # log-prob). Empty + untouched unless config.record_position_gaps is set.
    position_gaps = Float64[]
    for depth in 1:(length(observation) - 1)
        observed_unit = observation[depth + 1]
        next_scores = Dict{Tuple{label_type, Rhizomorph.StrandOrientation}, Float64}()
        next_predecessors = Dict{
            Tuple{label_type, Rhizomorph.StrandOrientation},
            Tuple{label_type, Rhizomorph.StrandOrientation}
        }()
        next_emissions = Dict{Tuple{label_type, Rhizomorph.StrandOrientation}, Float64}()

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

            # Bound successor generation to the top-B highest-weight transitions
            # BEFORE emission scoring (td-plqi). No-op when B >= out-degree (the
            # de Bruijn common case, out-degree <= 4). `total_out` is unchanged —
            # transition probabilities stay normalized against the FULL outgoing
            # mass, so a bounded expansion is a strict subset of the exact frontier.
            if length(transitions) > config.max_successors_per_state
                diagnostics[:successor_bounded] = get(diagnostics, :successor_bounded, 0) +
                                                  1
                transitions = _top_b_transitions(transitions, config.max_successors_per_state)
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
                    # Carry the emission (read-consistency) component of THIS path
                    # in lockstep with the Viterbi (max full-score) choice.
                    next_emissions[next_state] = get(active_emissions, state, 0.0) +
                                                 emission_score
                end
            end
        end

        if isempty(next_scores)
            break
        end

        # Beam pruning: cap the frontier to the top `beam_width` states by score
        # before it becomes the next `active_scores`. Without this, the reachable
        # (vertex, strand) set grows ~unboundedly with read-length depth on real
        # branchy graphs (21B allocations / hard crash on a 48 kb phage). Keeping
        # the top-K by score is the standard beam approximation; with the default
        # beam_width = typemax(Int) the guard never fires and the decoder stays
        # exact (B8 correction fixtures are byte-identical).
        if length(next_scores) > config.beam_width
            next_scores,
            next_predecessors = _prune_correction_beam(
                next_scores, next_predecessors, config.beam_width)
            # Keep the emission dict aligned to the width-beam survivors.
            next_emissions = Dict(
                state => next_emissions[state] for state in keys(next_scores))
            diagnostics[:beam_pruned] = get(diagnostics, :beam_pruned, 0) + 1
        end

        # Score-margin ("histogram") prune with emission exemption: keep a state
        # unless it is BOTH >Δ below the best full score AND >Δ below the best
        # cumulative emission (td-plqi). This bounds the GENERATING frontier to O(1)
        # in genome size on dense intermediate-k graphs (the improbable, read-
        # INCONSISTENT junk is pruned), while never dropping a read-consistent path —
        # so a real-but-rare (skewed-coverage) allele is protected from coverage-
        # driven pruning. A no-op under the default Inf margin (exact ML).
        if isfinite(config.beam_score_margin) && !isempty(next_scores)
            depth_best = maximum(values(next_scores))
            depth_best_emission = maximum(values(next_emissions))
            pre_margin = length(next_scores)
            next_scores, next_predecessors,
            next_emissions = _prune_correction_beam_by_margin(
                next_scores, next_predecessors, next_emissions,
                depth_best, depth_best_emission, config.beam_score_margin)
            if length(next_scores) < pre_margin
                diagnostics[:margin_pruned] = get(diagnostics, :margin_pruned, 0) + 1
            end
        end

        # Tier-2 opt-in: record this position's Viterbi margin (best - 2nd-best
        # surviving log-prob). Small margin => an ambiguous decode here; large =>
        # a confident call. Pure telemetry: it reads next_scores but never mutates
        # it or the chosen path, so the decode stays byte-identical when off.
        if config.record_position_gaps
            push!(position_gaps, _top2_score_gap(next_scores))
        end

        push!(predecessors_by_depth, next_predecessors)
        active_scores = next_scores
        active_emissions = next_emissions
        retained_count = length(active_scores)
        diagnostics[:retained_states] = retained_count
        diagnostics[:cumulative_retained_states] += retained_count
        diagnostics[:max_retained_states] = max(diagnostics[:max_retained_states], retained_count)
        diagnostics[:completed_steps] = depth

        if target_vertex === nothing
            best_state, best_score = _best_correction_state(active_scores)
            best_depth = depth
        else
            target_state,
            target_score = _best_correction_target_state(active_scores, target_vertex)
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
        # position_gaps is intentionally NOT stashed on a failed decode: there is no
        # path to align it to. `path === nothing` (+ diagnostics[:reason]) is the
        # failure discriminant — a consumer must not read absence of :position_gaps
        # as "flag was off".
        return Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diagnostics)
    end

    path = _reconstruct_correction_path(graph, best_state, best_depth, predecessors_by_depth)
    diagnostics[:path_length] = length(path.steps)
    if config.record_position_gaps
        # CONSUMPTION CONTRACT (per PR #400 review). `position_gaps[i]` is the
        # Viterbi margin (best − 2nd-best surviving log-prob) at the transition INTO
        # `path.steps[i+1]`; the start position `path.steps[1]` has no gap. Truncate
        # to `best_depth` so `length(position_gaps) == length(path.steps) - 1` in ALL
        # cases — under target anchoring the loop can run PAST the last on-path depth,
        # leaving trailing gaps for abandoned frontiers. TWO caveats for consumers:
        #   • An entry of `Inf` means "no surviving competitor at that depth" — a
        #     ROUTINE value on a collapsed/beam-pruned frontier — so consumers MUST
        #     `filter(isfinite, _)` (or clip) BEFORE any sum-based calibration metric
        #     (reliability/ECE/Brier); a raw `Inf` silently poisons those to Inf/NaN.
        #   • Under a finite `beam_width`/`beam_score_margin` this is a SURVIVOR
        #     margin (2nd-best among survivors), not the exact Viterbi margin.
        diagnostics[:position_gaps] = position_gaps[1:min(best_depth, length(position_gaps))]
    end
    return Rhizomorph.ViterbiDecodingResult(path, best_score, diagnostics)
end

# Keep only the top-`beam_width` states (by score, hash tie-break) of a
# state->score layer dict, IN PLACE. No-op when the dict already fits. Mirrors
# `_prune_correction_beam` but for a single phase layer of the pair-HMM.
function _beam_prune_layer!(scores::Dict{S, Float64}, beam_width::Int) where {S}
    length(scores) <= beam_width && return scores
    ordered = sort!(collect(scores); by = kv -> (last(kv), hash(first(kv))), rev = true)
    keep = Set{S}(first(kv) for kv in Iterators.take(ordered, beam_width))
    for state in collect(keys(scores))
        state in keep || delete!(scores, state)
    end
    return scores
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Indel-aware pair-HMM decode of one observation against a Rhizomorph weighted
graph. This is the gap-move counterpart of `_viterbi_correct_observation`,
selected when `config.indel_moves` is set.

The DP promotes the read index `i` AND an alignment phase `φ ∈ {M, I, D}` into the
state key `(vertex, strand, φ)` (the substitution decoder carries `i` implicitly
in its depth loop and has no phase), so length changes between the read and the
graph walk are scored by TRUE moves rather than a flat frameshift penalty:

  * **M** (match/mismatch) consumes one read unit AND advances one graph edge,
    scoring the existing Phred-aware emission `e_M`.
  * **I** (insertion) consumes one read unit and STAYS on the vertex (`e_I`,
    default `log(1/|alphabet|)`), bounded by `max_insertion_run`.
  * **D** (deletion) advances one graph edge and consumes NO read unit (`e_D = 0`),
    computed as a bounded Bellman-Ford relaxation of up to `deletion_max_run` hops
    within a read column (the scalable replacement for the legacy O(V³)
    Floyd-Warshall).

Affine transition masses come from the error model: `δ_I = ε·f_ins`,
`δ_D = ε·f_del` (substitution mass rides inside `e_M`), with geometric gap-extend
`γ_I`, `γ_D`. Setting the indel fractions to 0 sends every gap transition to
`log(0) = -Inf`, so the reachable frontier collapses to the pure-M substitution
path — Illumina falls out as the special case.

The result diagnostics include `:move_counts`, the ordered `:move_trace`, and the
parallel `:read_index_trace`. The latter two preserve the exact pair-HMM traceback
used to reconstruct length-changing per-base qualities. `:decoded_read_index` and
`:truncated` report whether the frontier consumed the complete observation.
"""
function _viterbi_correct_observation_indel(
        graph::MetaGraphsNext.MetaGraph,
        observation::AbstractVector,
        alphabet::Symbol;
        config::ViterbiCorrectionConfig,
        strand_mode::Symbol,
        quality_graph::MetaGraphsNext.MetaGraph = graph,
        transition_scoring::Symbol = :normalized_edge_weight
)::Rhizomorph.ViterbiDecodingResult
    if Graphs.nv(graph.graph) == 0
        throw(ArgumentError("cannot correct observations against an empty graph"))
    end

    label_type = _viterbi_graph_label_type(graph)
    State = Tuple{label_type, Rhizomorph.StrandOrientation}
    Cell = Tuple{Int, State, Symbol}
    n = length(observation)

    start_observed = first(observation)
    target_vertex = _emission_target_vertex(
        graph, observation, config, alphabet, strand_mode)
    start_candidates = _viterbi_start_candidates(
        graph, label_type, _viterbi_label_unit(start_observed), alphabet, strand_mode)

    alphabet_size = _viterbi_alphabet_size(alphabet)
    insertion_emission = function (observed_unit)
        if config.insertion_emission_logp === nothing
            return -log(alphabet_size)
        end
        return config.insertion_emission_logp(observed_unit, alphabet)
    end

    # Affine transition masses from the error model. δ_I/δ_D = 0 ⇒ log(0) = -Inf,
    # which drops the gap moves and collapses the DP to the substitution path.
    error_rate = config.error_rate
    delta_I = error_rate * config.insertion_fraction
    delta_D = error_rate * config.deletion_fraction
    gamma_I = config.insertion_extend_probability
    gamma_D = config.deletion_extend_probability
    T_MM = log(1.0 - delta_I - delta_D)
    T_MI = log(delta_I)
    T_MD = log(delta_D)
    T_II = log(gamma_I)
    T_IM = log1p(-gamma_I)
    T_DD = log(gamma_D)
    T_DM = log1p(-gamma_D)

    deletion_max_run = config.deletion_max_run
    max_insertion_run = config.max_insertion_run
    band = config.band_width
    finite_beam = config.beam_width != typemax(Int)

    diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_indel_pair_hmm,
        :indel_moves => true,
        :exact => (!finite_beam && band === nothing),
        :alphabet => alphabet,
        :strand_mode => strand_mode,
        :reverse_complement_support => _viterbi_supports_reverse_complement(alphabet),
        :read_length => n,
        :target_vertex => target_vertex,
        :start_strand => Rhizomorph._normalize_strand(config.start_strand),
        :score_domain => :log_probability,
        :transition_scoring => transition_scoring,
        :emission_scoring => _viterbi_graph_has_quality(quality_graph) ?
                             :quality_aware : :alphabet_parameterized,
        :move_counts => Dict{Symbol, Int}(:M => 0, :I => 0, :D => 0),
        :reached_target => target_vertex === nothing ? nothing : false
    )

    # Global backpointer: (read_index, state, phase) -> predecessor cell or nothing.
    backpointers = Dict{Cell, Union{Nothing, Cell}}()

    # Outgoing (target-state, logE) expansion for one state, matching the
    # substitution decoder's normalization (log of edge_weight / total_out).
    # NOTE: these nested helpers are closures, so their signatures/locals must not
    # annotate with the enclosing `State`/`label_type` locals (Julia forbids a local
    # variable in a closure type position). Untyped collections are correctness-
    # equivalent here.
    function _expand(state)
        vertex, strand = state
        transitions = Rhizomorph._get_valid_transitions(graph, vertex, strand)
        out = Tuple{Any, Float64}[]
        isempty(transitions) && return out
        total_out = Rhizomorph._total_outgoing_weight(graph, vertex, strand)
        (!isfinite(total_out) || total_out <= 0.0) && return out
        for transition in transitions
            next_vertex = convert(label_type, transition[:target_vertex])
            next_strand = Rhizomorph._normalize_strand(transition[:target_strand])
            edge_w = Rhizomorph._edge_transition_weight(transition[:edge_data])
            edge_w <= 0.0 && continue
            push!(out, ((next_vertex, next_strand), log(edge_w / total_out)))
        end
        return out
    end

    # Bounded Bellman-Ford deletion relaxation WITHIN one read column: open a
    # deletion from every M-state (run 1) then extend along graph edges up to
    # `deletion_max_run` hops. `deletion_max_run` + negative `log γ_D` guard against
    # unbounded no-progress loops on graph cycles.
    function _relax_deletions!(
            read_index,
            match_scores,
            match_net,
            del_scores,
            del_net,
            del_run
    )
        deletion_max_run <= 0 && return
        for (state, score) in match_scores
            for (next_state, logE) in _expand(state)
                cand = score + T_MD + logE
                isfinite(cand) || continue
                if !haskey(del_scores, next_state) || cand > del_scores[next_state]
                    del_scores[next_state] = cand
                    del_run[next_state] = 1
                    del_net[next_state] = match_net[state] + 1
                    backpointers[(read_index, next_state, :D)] = (read_index, state, :M)
                end
            end
        end
        for hop in 2:deletion_max_run
            frontier = [state for (state, run) in del_run if run == hop - 1]
            for state in frontier
                score = del_scores[state]
                for (next_state, logE) in _expand(state)
                    cand = score + T_DD + logE
                    isfinite(cand) || continue
                    if !haskey(del_scores, next_state) || cand > del_scores[next_state]
                        del_scores[next_state] = cand
                        del_run[next_state] = hop
                        del_net[next_state] = del_net[state] + 1
                        backpointers[(read_index, next_state, :D)] = (read_index, state, :D)
                    end
                end
            end
        end
    end

    # Layer 1: seed the match phase from the start candidates (no gap at the very
    # start — the read is anchored, matching the substitution decoder's seed), then
    # relax deletions within the column.
    match_scores = Dict{State, Float64}()
    ins_scores = Dict{State, Float64}()
    del_scores = Dict{State, Float64}()
    match_net = Dict{State, Int}()
    ins_net = Dict{State, Int}()
    del_net = Dict{State, Int}()
    ins_run = Dict{State, Int}()
    del_run = Dict{State, Int}()

    for vertex in start_candidates
        for strand in
            _viterbi_start_strands(graph, vertex, strand_mode, config.start_strand)
            state = (vertex, strand)
            score = _call_viterbi_state_emission_logp(
                quality_graph, config, start_observed, vertex, alphabet, strand_mode;
                count_length_penalty = false)
            if isfinite(score) &&
               (!haskey(match_scores, state) || score > match_scores[state])
                match_scores[state] = score
                match_net[state] = 0
                backpointers[(1, state, :M)] = nothing
            end
        end
    end
    if isempty(match_scores)
        diagnostics[:reason] = :no_finite_start_emission
        return Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diagnostics)
    end
    _relax_deletions!(1, match_scores, match_net, del_scores, del_net, del_run)

    # The read index of the layer currently held in match/ins/del_scores. If a noisy
    # read kills the whole frontier mid-decode we KEEP the last non-empty layer (its
    # index) and select the best-so-far endpoint from it, rather than mislabeling the
    # surviving cells with the full read length `n`.
    last_index = 1

    for read_index in 2:n
        observed_unit = observation[read_index]
        new_match = Dict{State, Float64}()
        new_ins = Dict{State, Float64}()
        new_del = Dict{State, Float64}()
        new_match_net = Dict{State, Int}()
        new_ins_net = Dict{State, Int}()
        new_del_net = Dict{State, Int}()
        new_ins_run = Dict{State, Int}()
        new_del_run = Dict{State, Int}()

        # M into (i, v): consume a read unit AND advance an edge, from any phase of
        # the previous column.
        for (phase, prev_scores, prev_net, trans) in (
            (:M, match_scores, match_net, T_MM),
            (:I, ins_scores, ins_net, T_IM),
            (:D, del_scores, del_net, T_DM)
        )
            for (state, score) in prev_scores
                base = score + trans
                isfinite(base) || continue
                for (next_state, logE) in _expand(state)
                    emission = _call_viterbi_state_emission_logp(
                        quality_graph, config, observed_unit, next_state[1],
                        alphabet, strand_mode; count_length_penalty = false)
                    cand = base + logE + emission
                    isfinite(cand) || continue
                    if !haskey(new_match, next_state) || cand > new_match[next_state]
                        new_match[next_state] = cand
                        new_match_net[next_state] = prev_net[state]
                        backpointers[(read_index, next_state, :M)] = (
                            read_index - 1, state, phase)
                    end
                end
            end
        end

        # I into (i, v): consume a read unit and STAY on v, from previous M or I,
        # bounded by the insertion-run cap.
        emission_ins = insertion_emission(observed_unit)
        if max_insertion_run >= 1 && isfinite(emission_ins)
            for (state, score) in match_scores
                cand = score + T_MI + emission_ins
                isfinite(cand) || continue
                if !haskey(new_ins, state) || cand > new_ins[state]
                    new_ins[state] = cand
                    new_ins_net[state] = match_net[state] - 1
                    new_ins_run[state] = 1
                    backpointers[(read_index, state, :I)] = (read_index - 1, state, :M)
                end
            end
            for (state, score) in ins_scores
                ins_run[state] >= max_insertion_run && continue
                cand = score + T_II + emission_ins
                isfinite(cand) || continue
                if !haskey(new_ins, state) || cand > new_ins[state]
                    new_ins[state] = cand
                    new_ins_net[state] = ins_net[state] - 1
                    new_ins_run[state] = ins_run[state] + 1
                    backpointers[(read_index, state, :I)] = (read_index - 1, state, :I)
                end
            end
        end

        # D within (i, v): relax deletions from the freshly-built match layer.
        _relax_deletions!(
            read_index, new_match, new_match_net, new_del, new_del_net, new_del_run)

        # Adaptive band: drop states whose net gap (graph-steps − read-index =
        # #deletions − #insertions) exceeds the half-width. `nothing` = unbounded
        # (exact). Applied AFTER the deletion relaxation so a banded state cannot
        # seed a new column.
        if band !== nothing
            for (scores, nets) in (
                (new_match, new_match_net), (new_ins, new_ins_net), (new_del, new_del_net))
                for state in collect(keys(scores))
                    if abs(nets[state]) > band
                        delete!(scores, state)
                    end
                end
            end
        end

        if finite_beam
            _beam_prune_layer!(new_match, config.beam_width)
            _beam_prune_layer!(new_ins, config.beam_width)
            _beam_prune_layer!(new_del, config.beam_width)
        end

        # If the whole frontier died this column, KEEP the previous (non-empty) layer
        # as the decode result and stop — do not swap in the empty layer.
        if isempty(new_match) && isempty(new_ins) && isempty(new_del)
            break
        end

        match_scores, ins_scores, del_scores = new_match, new_ins, new_del
        match_net, ins_net, del_net = new_match_net, new_ins_net, new_del_net
        ins_run, del_run = new_ins_run, new_del_run
        last_index = read_index
    end

    # Endpoint: the best-scoring (state, phase) at the final read index. When a
    # target vertex is pinned, restrict to states on it. Deterministic tie-break by
    # (score, vertex string, phase) so a beam tie is reproducible.
    best_score = -Inf
    best_cell::Union{Nothing, Cell} = nothing
    best_key = (-Inf, "", "")
    for (phase, scores) in ((:M, match_scores), (:I, ins_scores), (:D, del_scores))
        for (state, score) in scores
            (target_vertex !== nothing && state[1] != target_vertex) && continue
            key = (score, string(state[1]), string(phase))
            if best_cell === nothing || score > best_score ||
               (score == best_score && key < best_key)
                best_score = score
                best_cell = (last_index, state, phase)
                best_key = key
            end
        end
    end

    if best_cell === nothing
        diagnostics[:reason] = target_vertex === nothing ? :no_surviving_path :
                               :target_unreachable
        return Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diagnostics)
    end

    # Reconstruct the graph walk: follow backpointers to the start, then keep the
    # vertex of every M and D move (both traverse a real edge); an I move stays on
    # its vertex and contributes NO new vertex.
    chain = Cell[]
    cursor::Union{Nothing, Cell} = best_cell
    while cursor !== nothing
        push!(chain, cursor)
        cursor = backpointers[cursor]
    end
    reverse!(chain)

    path_states = State[]
    move_trace = Symbol[]
    read_index_trace = Int[]
    for (index, (read_index, state, phase)) in enumerate(chain)
        diagnostics[:move_counts][phase] += 1
        push!(move_trace, phase)
        push!(read_index_trace, read_index)
        if index == 1 || phase != :I
            push!(path_states, state)
        end
    end
    diagnostics[:move_trace] = move_trace
    diagnostics[:read_index_trace] = read_index_trace

    path = Rhizomorph._build_graph_path_from_vertices(graph, path_states)
    diagnostics[:path_length] = length(path.steps)
    diagnostics[:decoded_read_index] = last_index
    # Truncation is a first-class diagnostic (PR #407 review, silent-failure I1). On
    # the free-endpoint decode a noisy read can kill the whole frontier before the
    # last read unit; we KEEP the best-so-far prefix as the result, but stamp
    # `:truncated` so callers do not have to INFER truncation by diffing
    # `:decoded_read_index` against the read length. A verbosity-gated `@warn` makes
    # a partial decode visible to an operator watching the reads/debug streams.
    truncated = last_index < n
    diagnostics[:truncated] = truncated
    if truncated && target_vertex === nothing && config.verbosity in ("debug", "reads")
        @warn "indel-aware decode truncated: the read frontier died at unit " *
              "$(last_index) of $(n); the returned path covers only the decoded " *
              "prefix (see diagnostics[:decoded_read_index])."
    end
    if target_vertex !== nothing
        diagnostics[:reached_target] = true
    end
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
    # NOTE: deliberately do NOT unwrap a QualityObservation here. Pinning the
    # decode's endpoint to the read's observed last k-mer would forbid correcting
    # it — if that k-mer is an in-graph error (present from other reads' errors),
    # the ML path would be forced onto the error instead of the higher-likelihood
    # true vertex. Leaving the endpoint free (target === nothing for wrapped
    # observations) lets the ML path + per-base emission choose it, which is the
    # whole point of the correction core. The start (below) is still anchored for
    # efficiency, since a read's first k-mer is a reasonable, usually-correct
    # threading anchor and falls back to all labels when absent from the graph.
    # O(1) `haskey` membership (td-y4oj) replaces the O(V) `collect(labels)` +
    # linear `in labels` scan that ran per read. `haskey(graph, x)` is the graph's
    # vertex-label lookup and matches `x in collect(labels(graph))` exactly (a
    # wrapped QualityObservation candidate is not a vertex label ⇒ `false` ⇒
    # endpoint left free, preserving the documented behavior).
    candidate = last(observation)
    if haskey(graph, candidate)
        return candidate
    end
    if strand_mode == :canonical && _viterbi_supports_reverse_complement(alphabet)
        canonical = _viterbi_canonical_unit(candidate, alphabet)
        if haskey(graph, canonical)
            return canonical
        end
    end
    return nothing
end

# Resolve start candidates with O(1) `haskey` lookups (td-y4oj). The common case
# (the read's first k-mer is a graph vertex) is a single hash probe; only the
# rare fallback (start k-mer absent ⇒ the decode may thread from any vertex)
# materializes the O(V) label vector. Byte-identical to the prior O(V) linear
# `in labels` scan.
function _viterbi_start_candidates(
        graph::MetaGraphsNext.MetaGraph,
        ::Type{T},
        observed_unit::Any,
        alphabet::Symbol,
        strand_mode::Symbol
)::Vector{T} where {T}
    if haskey(graph, observed_unit)
        return T[convert(T, observed_unit)]
    end
    if strand_mode == :canonical && _viterbi_supports_reverse_complement(alphabet)
        canonical = _viterbi_canonical_unit(observed_unit, alphabet)
        if haskey(graph, canonical)
            return T[convert(T, canonical)]
        end
    end
    return collect(MetaGraphsNext.labels(graph))
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
        strand_mode::Symbol;
        count_length_penalty::Bool = true
)::Float64
    quality = _viterbi_emission_quality(graph, observed_unit)
    direct = _call_viterbi_emission_logp(
        config,
        observed_unit,
        node,
        alphabet;
        quality = quality,
        count_length_penalty = count_length_penalty
    )
    if strand_mode == :canonical && _viterbi_supports_reverse_complement(alphabet)
        rc_node = _viterbi_reverse_complement_unit(node, alphabet)
        rc = _call_viterbi_emission_logp(
            config,
            observed_unit,
            rc_node,
            alphabet;
            quality = quality,
            count_length_penalty = count_length_penalty
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

# Best-minus-second-best of a state->score Dict's values (the per-position Viterbi
# margin, td-4osf / Tier-2 calibration). Single pass, no allocation. Returns Inf
# when fewer than two states survive (no competitor => maximally confident). NOTE:
# Inf is a ROUTINE value (any collapsed frontier), so downstream calibration must
# filter/clip it before averaging — see the consumption contract at the
# diagnostics[:position_gaps] stash site.
function _top2_score_gap(scores)::Float64
    best = -Inf
    second = -Inf
    for v in values(scores)
        if v > best
            second = best
            best = v
        elseif v > second
            second = v
        end
    end
    return isfinite(second) ? best - second : Inf
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
    # Character-aware indexing (see `default_viterbi_emission_logp`): byte offsets
    # into a multibyte UTF-8 string throw `StringIndexError`, so iterate the
    # decoded `Char`s to keep the distance well-defined for arbitrary Unicode.
    left_chars = collect(left)
    right_chars = collect(right)
    shared = min(length(left_chars), length(right_chars))
    distance = abs(length(left_chars) - length(right_chars))
    for index in 1:shared
        distance += left_chars[index] == right_chars[index] ? 0 : 1
    end
    return distance
end
