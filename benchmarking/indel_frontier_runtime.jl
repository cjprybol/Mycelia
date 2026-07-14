# Branching/frontier calibration for sparse-rung indel scheduling (td-jt7r.2).
#
# Full calibration (intentionally expensive):
#   LD_LIBRARY_PATH='' julia --project=. benchmarking/indel_frontier_runtime.jl
#
# Fast harness smoke (synthetic graphs only; no fixed-toy decode):
#   LD_LIBRARY_PATH='' julia --project=. benchmarking/indel_frontier_runtime.jl \
#     --smoke --output-dir /tmp/indel-frontier-runtime-smoke
#
# The scheduler is calibrated only against warmed pair-HMM runtime and topology-
# only frontier work. Correction accuracy is neither measured nor available to
# this script, so it cannot leak into threshold selection.

import BioSequences
import CairoMakie
import CSV
import DataFrames
import FASTX
import Graphs
import Mycelia
import Random
import SHA
import Statistics

const INDEL_FRONTIER_K = 9
const INDEL_FRONTIER_TOY_K = 31
const INDEL_FRONTIER_GENOME_LENGTH = 2_000
const INDEL_FRONTIER_SOURCE_READ_LENGTH = 1_200
const INDEL_FRONTIER_COVERAGE = 8
const INDEL_FRONTIER_ERROR_RATE = 0.05
const INDEL_FRONTIER_FIXTURE_SEED = 42
const INDEL_FRONTIER_WINDOW_BASES = (250, 500, 1_200)
const INDEL_FRONTIER_REPEATS = 5
const INDEL_FRONTIER_MISSING_DEPENDENCY_SENTINEL = "MISSING"
# A historical 200 ms value informed early topology calibration and remains
# aspirational context only; it is not an executable or portable wall-clock
# SLA. The executable runtime-only affordability gate requires the 10,001-
# vertex/500 bp warmed p95 to be at most 1 s and at most 3x the same-run 2,001-
# vertex control. The absolute backstop prevents a uniformly pathological host
# from satisfying only the relative check.
const INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_SLOWDOWN = 3.0
const INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_MS = 1_000.0
const INDEL_FRONTIER_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "td-jt7r-2-frontier-runtime"
)
const INDEL_FRONTIER_ARTIFACT_NAMES = (
    "indel_frontier_runtime_summary.csv",
    "indel_frontier_runtime_replicates.csv",
    "indel_frontier_runtime.png",
    "indel_frontier_runtime.svg",
    "indel_frontier_runtime_manifest.csv",
)

struct IndelFrontierCalibrationGraph
    graph_id::String
    graph::Any
    probe_records::Vector{FASTX.FASTQ.Record}
    k::Int
    strand_mode::Symbol
    graph_source::Symbol
    cleanup_stats::Dict{String, Any}
end

function main(args::Vector{String} = ARGS)::Nothing
    options = _indel_frontier_parse_args(args)
    repeats = something(
        options.repeats,
        options.smoke ? 1 : INDEL_FRONTIER_REPEATS,
    )
    _indel_frontier_remove_prior_artifacts(
        options.output_dir, INDEL_FRONTIER_ARTIFACT_NAMES
    )
    run_provenance = _indel_frontier_run_provenance(options.smoke, repeats)

    graph_cases = _indel_frontier_calibration_graphs(options.smoke)
    window_bases = options.smoke ? (50,) : INDEL_FRONTIER_WINDOW_BASES
    summary, replicates = _indel_frontier_run_matrix(
        graph_cases, window_bases, repeats
    )
    if !options.smoke
        _indel_frontier_assert_topology_controls(summary)
    end

    provenance = merge(
        run_provenance,
        _indel_frontier_affordability_provenance(options.smoke, summary),
    )
    staging_dir = Base.Filesystem.mktempdir(
        options.output_dir; prefix = ".frontier-runtime-staging-"
    )
    try
        summary_staging_path = joinpath(
            staging_dir, INDEL_FRONTIER_ARTIFACT_NAMES[1]
        )
        replicates_staging_path = joinpath(
            staging_dir, INDEL_FRONTIER_ARTIFACT_NAMES[2]
        )
        CSV.write(summary_staging_path, summary)
        CSV.write(replicates_staging_path, replicates)
        figure_paths = _indel_frontier_write_figure(summary, staging_dir)
        manifest = DataFrames.DataFrame([
            merge(
                provenance,
                (
                    summary_sha256 = _indel_frontier_file_sha256(
                        summary_staging_path
                    ),
                    replicates_sha256 = _indel_frontier_file_sha256(
                        replicates_staging_path
                    ),
                    png_sha256 = _indel_frontier_file_sha256(figure_paths.png),
                    svg_sha256 = _indel_frontier_file_sha256(figure_paths.svg),
                ),
            ),
        ])
        CSV.write(
            joinpath(staging_dir, INDEL_FRONTIER_ARTIFACT_NAMES[5]), manifest
        )
        _indel_frontier_assert_provenance_unchanged(
            run_provenance, options.smoke, repeats
        )
        _indel_frontier_publish_artifacts(
            staging_dir, options.output_dir, INDEL_FRONTIER_ARTIFACT_NAMES
        )
    finally
        Base.rm(staging_dir; recursive = true, force = true)
    end

    summary_path = joinpath(options.output_dir, INDEL_FRONTIER_ARTIFACT_NAMES[1])
    replicates_path = joinpath(
        options.output_dir, INDEL_FRONTIER_ARTIFACT_NAMES[2]
    )
    png_path = joinpath(options.output_dir, INDEL_FRONTIER_ARTIFACT_NAMES[3])
    svg_path = joinpath(options.output_dir, INDEL_FRONTIER_ARTIFACT_NAMES[4])
    manifest_path = joinpath(
        options.output_dir, INDEL_FRONTIER_ARTIFACT_NAMES[5]
    )

    println("Wrote branching/frontier runtime calibration artifacts:")
    println("  summary:    $(summary_path)")
    println("  replicates: $(replicates_path)")
    println("  png:        $(png_path)")
    println("  svg:        $(svg_path)")
    println("  manifest:   $(manifest_path)")
    println("No correction-accuracy metric was computed or used.")
    return nothing
end

function _indel_frontier_calibration_graphs(
        smoke::Bool
)::Vector{IndelFrontierCalibrationGraph}
    branching = _indel_frontier_branching_graph()
    if smoke
        return IndelFrontierCalibrationGraph[
            branching,
            _indel_frontier_linear_graph(300, "linear_300_smoke"),
        ]
    end

    linear_2k = _indel_frontier_linear_graph(2_001, "linear_2001")
    linear_10k = _indel_frontier_linear_graph(10_001, "linear_10001")
    reads, _ = _indel_frontier_fixed_fixture()
    raw_graph = Mycelia.Rhizomorph.build_qualmer_graph(
        reads, INDEL_FRONTIER_TOY_K; mode = :doublestrand
    )
    cleaned_graph = deepcopy(raw_graph)
    cleanup_stats = Mycelia.Rhizomorph.clean_corrector_graph!(
        cleaned_graph; k = INDEL_FRONTIER_TOY_K
    )

    return IndelFrontierCalibrationGraph[
        branching,
        linear_2k,
        linear_10k,
        IndelFrontierCalibrationGraph(
            "fixed_toy_k31_raw",
            raw_graph,
            reads,
            INDEL_FRONTIER_TOY_K,
            :doublestrand,
            :raw_fixed_toy,
            Dict{String, Any}(),
        ),
        IndelFrontierCalibrationGraph(
            "fixed_toy_k31_cleaned",
            cleaned_graph,
            reads,
            INDEL_FRONTIER_TOY_K,
            :doublestrand,
            :cleaned_fixed_toy,
            cleanup_stats,
        ),
    ]
end

function _indel_frontier_linear_graph(
        target_vertices::Int,
        graph_id::String,
)::IndelFrontierCalibrationGraph
    target_vertices > 0 || throw(ArgumentError("target_vertices must be positive"))
    cycle = _indel_frontier_de_bruijn_sequence(INDEL_FRONTIER_K)
    sequence_length = target_vertices + INDEL_FRONTIER_K - 1
    sequence_length <= length(cycle) || error(
        "de Bruijn prefix is too short for $(target_vertices) vertices"
    )
    sequence = first(cycle, sequence_length)
    record = _indel_frontier_fastq_record(graph_id, sequence)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[record], INDEL_FRONTIER_K; mode = :singlestrand
    )
    n_vertices = Graphs.nv(graph.graph)
    max_outdegree = maximum(
        Graphs.outdegree(graph.graph, vertex) for vertex in Graphs.vertices(graph.graph)
    )
    n_vertices == target_vertices || error(
        "$(graph_id) expected $(target_vertices) unique vertices, observed " *
        "$(n_vertices)"
    )
    max_outdegree <= 1 || error(
        "$(graph_id) is not linear: max outdegree=$(max_outdegree)"
    )
    return IndelFrontierCalibrationGraph(
        graph_id,
        graph,
        FASTX.FASTQ.Record[record],
        INDEL_FRONTIER_K,
        :singlestrand,
        :synthetic_linear,
        Dict{String, Any}(),
    )
end

function _indel_frontier_branching_graph()::IndelFrontierCalibrationGraph
    alphabet = collect("ACGT")
    records = FASTX.FASTQ.Record[]
    record_index = 0
    for a in alphabet, b in alphabet, c in alphabet, d in alphabet
        record_index += 1
        sequence = string(a, b, c, d)
        push!(
            records,
            _indel_frontier_fastq_record("branch_edge_$(record_index)", sequence),
        )
    end
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        records, 3; mode = :singlestrand
    )
    probe_cycle = _indel_frontier_de_bruijn_sequence(3)
    probe_sequence = first(
        repeat(probe_cycle, cld(1_300, length(probe_cycle))), 1_300
    )
    probe_record = _indel_frontier_fastq_record(
        "maximally_branching_probe", probe_sequence
    )

    n_vertices = Graphs.nv(graph.graph)
    outdegrees = [
        Graphs.outdegree(graph.graph, vertex) for vertex in Graphs.vertices(graph.graph)
    ]
    n_vertices == 64 || error(
        "maximally branching k=3 graph expected 64 vertices, observed $(n_vertices)"
    )
    all(==(4), outdegrees) || error(
        "maximally branching k=3 graph expected outdegree 4 at every vertex"
    )
    return IndelFrontierCalibrationGraph(
        "branching_k3_complete",
        graph,
        FASTX.FASTQ.Record[probe_record],
        3,
        :singlestrand,
        :synthetic_branching,
        Dict{String, Any}(),
    )
end

function _indel_frontier_de_bruijn_sequence(order::Int)::String
    order > 0 || throw(ArgumentError("de Bruijn order must be positive"))
    alphabet = collect("ACGT")
    alphabet_size = length(alphabet)
    work = zeros(Int, alphabet_size * order + 1)
    indices = Int[]

    function visit(t::Int, period::Int)::Nothing
        if t > order
            if order % period == 0
                append!(indices, @view work[2:(period + 1)])
            end
            return nothing
        end
        work[t + 1] = work[t - period + 1]
        visit(t + 1, period)
        for index in (work[t - period + 1] + 1):(alphabet_size - 1)
            work[t + 1] = index
            visit(t + 1, t)
        end
        return nothing
    end

    visit(1, 1)
    return String([alphabet[index + 1] for index in indices])
end

function _indel_frontier_fixed_fixture()::Tuple{
        Vector{FASTX.FASTQ.Record}, BioSequences.LongDNA{4}}
    reference_record = Mycelia.random_fasta_record(
        moltype = :DNA,
        seed = INDEL_FRONTIER_FIXTURE_SEED,
        L = INDEL_FRONTIER_GENOME_LENGTH,
    )
    reference = FASTX.sequence(BioSequences.LongDNA{4}, reference_record)
    rng = Random.MersenneTwister(INDEL_FRONTIER_FIXTURE_SEED)
    Random.seed!(INDEL_FRONTIER_FIXTURE_SEED)
    n_reads = ceil(
        Int,
        INDEL_FRONTIER_COVERAGE * INDEL_FRONTIER_GENOME_LENGTH /
        INDEL_FRONTIER_SOURCE_READ_LENGTH,
    )
    reads = FASTX.FASTQ.Record[]
    for read_index in 1:n_reads
        start_position = rand(
            rng,
            1:(INDEL_FRONTIER_GENOME_LENGTH -
               INDEL_FRONTIER_SOURCE_READ_LENGTH + 1),
        )
        fragment = reference[
            start_position:(start_position + INDEL_FRONTIER_SOURCE_READ_LENGTH - 1)
        ]
        if rand(rng, Bool)
            fragment = BioSequences.reverse_complement(fragment)
        end
        observed, qualities = Mycelia.observe(
            fragment; error_rate = INDEL_FRONTIER_ERROR_RATE, tech = :nanopore
        )
        isempty(observed) && continue
        quality_string = String([Char(quality + 33) for quality in qualities])
        push!(
            reads,
            FASTX.FASTQ.Record(
                "nanopore_read_$(read_index)", string(observed), quality_string
            ),
        )
    end
    return reads, reference
end

function _indel_frontier_fastq_record(
        identifier::String,
        sequence::String,
)::FASTX.FASTQ.Record
    return FASTX.FASTQ.Record(identifier, sequence, repeat("I", length(sequence)))
end

function _indel_frontier_run_matrix(
        graph_cases::Vector{IndelFrontierCalibrationGraph},
        window_bases::Tuple{Vararg{Int}},
        repeats::Int,
)::Tuple{DataFrames.DataFrame, DataFrames.DataFrame}
    summary_rows = NamedTuple[]
    replicate_rows = NamedTuple[]
    for graph_case in graph_cases
        structural = _indel_frontier_structural_metrics(graph_case.graph)
        for n_bases in window_bases
            window = _indel_frontier_anchored_window(graph_case, n_bases)
            observations = _indel_frontier_quality_observations(
                window.record, graph_case.k
            )
            config = _indel_frontier_config(
                length(observations), structural.n_vertices, graph_case.strand_mode
            )
            probe = Mycelia._probe_indel_frontier(
                graph_case.graph,
                observations,
                :DNA;
                config = config,
                strand_mode = graph_case.strand_mode,
                work_limit = Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT,
            )
            measurement = _indel_frontier_measure_decode(
                graph_case.graph, observations, config, repeats
            )
            for row in measurement.replicates
                push!(
                    replicate_rows,
                    (
                        graph_id = graph_case.graph_id,
                        graph_source = string(graph_case.graph_source),
                        k = graph_case.k,
                        window_bases = n_bases,
                        n_observations = length(observations),
                        window_start = window.start,
                        phase = row.phase,
                        replicate = row.replicate,
                        elapsed_ms = row.elapsed_ms,
                        algorithm = string(row.algorithm),
                        path_available = row.path_available,
                        pair_hmm_path_valid = row.pair_hmm_path_valid,
                        pair_hmm_valid = row.pair_hmm_valid,
                        truncated = row.truncated,
                        completed_columns = row.completed_columns,
                        decoded_read_index = row.decoded_read_index,
                        terminal_contract_valid = row.terminal_contract_valid,
                        trace_complete = row.trace_complete,
                        complete = row.complete,
                        full_decode = row.full_decode,
                        frontier_area = row.frontier_area,
                        edge_expansions = row.edge_expansions,
                        peak_frontier = row.peak_frontier,
                    ),
                )
            end
            cleanup = graph_case.cleanup_stats
            push!(
                summary_rows,
                (
                    graph_id = graph_case.graph_id,
                    graph_source = string(graph_case.graph_source),
                    k = graph_case.k,
                    strand_mode = string(graph_case.strand_mode),
                    window_bases = n_bases,
                    n_observations = length(observations),
                    window_start = window.start,
                    n_vertices = structural.n_vertices,
                    n_edges = structural.n_edges,
                    branch_vertices = structural.branch_vertices,
                    join_vertices = structural.join_vertices,
                    branch_fraction = structural.branch_fraction,
                    join_fraction = structural.join_fraction,
                    mean_outdegree = structural.mean_outdegree,
                    p95_outdegree = structural.p95_outdegree,
                    max_outdegree = structural.max_outdegree,
                    probe_anchored = probe.anchored,
                    probe_reason = string(probe.reason),
                    probe_completed_columns = probe.completed_columns,
                    probe_frontier_area = probe.frontier_area,
                    probe_edge_expansions = probe.edge_expansions,
                    probe_peak_frontier = probe.peak_frontier,
                    probe_frontier_work = probe.frontier_work,
                    frontier_work_limit =
                        Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT,
                    beam_width = config.beam_width,
                    warmup_ms = measurement.warmup_ms,
                    median_ms = measurement.median_ms,
                    p95_ms = measurement.p95_ms,
                    min_ms = measurement.min_ms,
                    max_ms = measurement.max_ms,
                    repeats = repeats,
                    pair_hmm_samples = measurement.samples,
                    pair_hmm_path_valid = measurement.path_valid,
                    pair_hmm_valid = measurement.valid,
                    pair_hmm_terminal_contract_valid =
                        measurement.terminal_contract_valid,
                    pair_hmm_trace_complete = measurement.trace_complete,
                    pair_hmm_truncated = measurement.truncated,
                    pair_hmm_truncated_samples = measurement.truncated_samples,
                    pair_hmm_complete = measurement.complete,
                    pair_hmm_complete_samples = measurement.complete_samples,
                    pair_hmm_full_decode = measurement.full_decode,
                    cleanup_vertices_before = get(
                        cleanup, "graph_cleanup_vertices_before", missing
                    ),
                    cleanup_vertices_after = get(
                        cleanup, "graph_cleanup_vertices_after", missing
                    ),
                    cleanup_tips_removed = get(
                        cleanup, "graph_cleanup_tips_removed", missing
                    ),
                    cleanup_bubbles_collapsed = get(
                        cleanup, "graph_cleanup_bubbles_collapsed", missing
                    ),
                    accuracy_metric_used = false,
                ),
            )
            println(
                "graph=$(graph_case.graph_id),window=$(n_bases)bp," *
                "nv=$(structural.n_vertices),branch=" *
                "$(round(structural.branch_fraction; digits = 4))," *
                "frontier_work=$(probe.frontier_work),reason=$(probe.reason)," *
                "median_ms=$(round(measurement.median_ms; digits = 3))," *
                "p95_ms=$(round(measurement.p95_ms; digits = 3))"
            )
        end
    end
    return DataFrames.DataFrame(summary_rows), DataFrames.DataFrame(replicate_rows)
end

function _indel_frontier_structural_metrics(graph::Any)::NamedTuple
    internal = graph.graph
    vertices = collect(Graphs.vertices(internal))
    outdegrees = Int[Graphs.outdegree(internal, vertex) for vertex in vertices]
    indegrees = Int[Graphs.indegree(internal, vertex) for vertex in vertices]
    n_vertices = length(vertices)
    branch_vertices = count(>(1), outdegrees)
    join_vertices = count(>(1), indegrees)
    return (
        n_vertices = n_vertices,
        n_edges = Graphs.ne(internal),
        branch_vertices = branch_vertices,
        join_vertices = join_vertices,
        branch_fraction = n_vertices == 0 ? 0.0 : branch_vertices / n_vertices,
        join_fraction = n_vertices == 0 ? 0.0 : join_vertices / n_vertices,
        mean_outdegree = isempty(outdegrees) ? 0.0 : Statistics.mean(outdegrees),
        p95_outdegree = isempty(outdegrees) ?
                        0.0 : Statistics.quantile(outdegrees, 0.95),
        max_outdegree = isempty(outdegrees) ? 0 : maximum(outdegrees),
    )
end

function _indel_frontier_anchored_window(
        graph_case::IndelFrontierCalibrationGraph,
        n_bases::Int,
)::NamedTuple
    n_bases >= graph_case.k || throw(
        ArgumentError("window length must be at least k=$(graph_case.k)")
    )
    for record in graph_case.probe_records
        sequence = FASTX.sequence(String, record)
        length(sequence) < n_bases && continue
        max_start = length(sequence) - n_bases + 1
        for start in 1:max_start
            candidate = _indel_frontier_slice_record(record, start, n_bases)
            observations = _indel_frontier_quality_observations(candidate, graph_case.k)
            if !isempty(observations) && haskey(graph_case.graph, first(observations).kmer)
                return (record = candidate, start = start)
            end
        end
    end
    error(
        "no anchored $(n_bases) bp window found for graph $(graph_case.graph_id)"
    )
end

function _indel_frontier_slice_record(
        record::FASTX.FASTQ.Record,
        start::Int,
        n_bases::Int,
)::FASTX.FASTQ.Record
    sequence = FASTX.sequence(String, record)
    quality = String(FASTX.quality(record))
    stop = start + n_bases - 1
    stop <= length(sequence) || throw(BoundsError(sequence, start:stop))
    return FASTX.FASTQ.Record(
        "$(FASTX.identifier(record))_$(start)_$(stop)",
        sequence[start:stop],
        quality[start:stop],
    )
end

function _indel_frontier_quality_observations(
        read::FASTX.FASTQ.Record,
        k::Int,
)::Vector{Mycelia.QualityObservation}
    sequence_string = FASTX.sequence(String, read)
    alphabet = Mycelia.detect_alphabet(sequence_string)
    sequence_type = Mycelia.alphabet_to_biosequence_type(alphabet)
    sequence = Mycelia.extract_typed_sequence(read, sequence_type)
    kmers = collect(Mycelia._record_kmer_iterator(sequence_type, k, sequence))
    quality_scores = collect(FASTX.quality_scores(read))
    observations = Vector{Mycelia.QualityObservation}(undef, length(kmers))
    for (index, kmer) in enumerate(kmers)
        lo = clamp(index, 1, length(quality_scores))
        hi = clamp(index + k - 1, 1, length(quality_scores))
        observations[index] = Mycelia.QualityObservation(
            kmer, UInt8.(@view quality_scores[lo:hi])
        )
    end
    return observations
end

function _indel_frontier_config(
        n_observations::Int,
        n_vertices::Int,
        strand_mode::Symbol,
)::Mycelia.ViterbiCorrectionConfig
    profile = Mycelia.indel_error_profile(:nanopore)
    knobs = Mycelia.Rhizomorph._corrector_strategy_knobs(:scalable)
    beam_width = Mycelia._auto_beam_width(n_observations, n_vertices)
    beam_is_exact = beam_width == typemax(Int)
    return Mycelia.ViterbiCorrectionConfig(
        alphabet = :DNA,
        strand_mode = strand_mode,
        max_steps = n_observations - 1,
        beam_width = beam_width,
        max_successors_per_state = beam_is_exact ?
                                   typemax(Int) : Mycelia._AUTO_SUCCESSOR_BOUND,
        beam_score_margin = beam_is_exact ?
                            Inf : Mycelia._AUTO_BEAM_SCORE_MARGIN,
        error_rate = profile.base_error_rate,
        indel_moves = true,
        insertion_fraction = profile.insertion_fraction,
        deletion_fraction = profile.deletion_fraction,
        insertion_extend_probability = profile.insertion_extend_probability,
        deletion_extend_probability = profile.deletion_extend_probability,
        deletion_max_run = knobs.deletion_max_run,
        max_insertion_run = knobs.max_insertion_run,
        band_width = knobs.band_width,
    )
end

function _indel_frontier_measure_decode(
        graph::Any,
        observations::Vector{Mycelia.QualityObservation},
        config::Mycelia.ViterbiCorrectionConfig,
        repeats::Int,
)::NamedTuple
    weighted_graph = Mycelia.build_correction_weighted_graph(graph; config = config)
    Base.GC.gc()
    warm_start = time_ns()
    warm = Mycelia.correct_observations(
        graph,
        [observations];
        config = config,
        weighted_graph = weighted_graph,
    )
    warmup_ms = (time_ns() - warm_start) / 1.0e6
    warm_path = only(warm.paths)
    rows = NamedTuple[
        _indel_frontier_replicate_row(
            "warmup", 0, warmup_ms, warm_path, length(observations)
        ),
    ]
    elapsed_ms = Float64[]
    for replicate in 1:repeats
        Base.GC.gc()
        start_ns = time_ns()
        result = Mycelia.correct_observations(
            graph,
            [observations];
            config = config,
            weighted_graph = weighted_graph,
        )
        elapsed = (time_ns() - start_ns) / 1.0e6
        path = only(result.paths)
        push!(elapsed_ms, elapsed)
        push!(
            rows,
            _indel_frontier_replicate_row(
                "measurement",
                replicate,
                elapsed,
                path,
                length(observations),
            ),
        )
    end

    path_valid = all(row.pair_hmm_path_valid for row in rows)
    valid = all(row.pair_hmm_valid for row in rows)
    terminal_contract_valid = all(
        row.terminal_contract_valid for row in rows)
    trace_complete = all(row.trace_complete for row in rows)
    truncated_samples = count(row.truncated for row in rows)
    complete_samples = count(row.complete for row in rows)
    truncated = truncated_samples > 0
    complete = complete_samples == length(rows)

    return (
        warmup_ms = warmup_ms,
        median_ms = Statistics.median(elapsed_ms),
        p95_ms = Statistics.quantile(elapsed_ms, 0.95),
        min_ms = minimum(elapsed_ms),
        max_ms = maximum(elapsed_ms),
        samples = length(rows),
        path_valid = path_valid,
        valid = valid,
        terminal_contract_valid = terminal_contract_valid,
        trace_complete = trace_complete,
        truncated = truncated,
        truncated_samples = truncated_samples,
        complete = complete,
        complete_samples = complete_samples,
        full_decode = valid,
        replicates = rows,
    )
end

function _indel_frontier_trace_complete(
        move_trace::Any,
        read_index_trace::Any,
        expected_columns::Int,
)::Bool
    move_trace isa AbstractVector{Symbol} || return false
    read_index_trace isa AbstractVector{Int} || return false
    isempty(move_trace) && return false
    length(move_trace) == length(read_index_trace) || return false
    first(move_trace) == :M || return false
    first(read_index_trace) == 1 || return false

    previous_read_index = 1
    for trace_index in 2:length(move_trace)
        phase = move_trace[trace_index]
        read_index = read_index_trace[trace_index]
        expected_read_index = phase == :D ?
                              previous_read_index : previous_read_index + 1
        phase in (:M, :I, :D) || return false
        read_index == expected_read_index || return false
        1 <= read_index <= expected_columns || return false
        previous_read_index = read_index
    end
    return previous_read_index == expected_columns
end

function _indel_frontier_replicate_row(
        phase::String,
        replicate::Int,
        elapsed_ms::Float64,
        path::Any,
        expected_columns::Int,
)::NamedTuple
    diagnostics = path.diagnostics
    algorithm = get(diagnostics, :algorithm, :missing)
    path_available = path.path !== nothing
    pair_hmm_path_valid = algorithm == :viterbi_indel_pair_hmm && path_available
    pair_hmm_path_valid || error(
        "pair-HMM $(phase) sample $(replicate) failed validation: " *
        "algorithm=$(algorithm), path_available=$(path_available)"
    )
    truncated_field = get(diagnostics, :truncated, nothing)
    completed_columns_field = get(diagnostics, :completed_columns, nothing)
    decoded_read_index_field = get(
        diagnostics, :decoded_read_index, nothing)
    terminal_contract_valid =
        truncated_field === false &&
        typeof(completed_columns_field) === Int &&
        completed_columns_field == expected_columns &&
        typeof(decoded_read_index_field) === Int &&
        decoded_read_index_field == expected_columns
    trace_complete = _indel_frontier_trace_complete(
        get(diagnostics, :move_trace, nothing),
        get(diagnostics, :read_index_trace, nothing),
        expected_columns,
    )
    truncated = truncated_field === true
    completed_columns = typeof(completed_columns_field) === Int ?
                        completed_columns_field : 0
    decoded_read_index = typeof(decoded_read_index_field) === Int ?
                         decoded_read_index_field : 0
    complete = terminal_contract_valid && trace_complete
    pair_hmm_valid = pair_hmm_path_valid && complete
    return (
        phase = phase,
        replicate = replicate,
        elapsed_ms = elapsed_ms,
        algorithm = algorithm,
        path_available = path_available,
        pair_hmm_path_valid = pair_hmm_path_valid,
        pair_hmm_valid = pair_hmm_valid,
        truncated = truncated,
        completed_columns = completed_columns,
        decoded_read_index = decoded_read_index,
        terminal_contract_valid = terminal_contract_valid,
        trace_complete = trace_complete,
        complete = complete,
        full_decode = pair_hmm_valid && complete,
        frontier_area = get(diagnostics, :frontier_area, 0),
        edge_expansions = get(diagnostics, :edge_expansions, 0),
        peak_frontier = get(diagnostics, :peak_frontier, 0),
    )
end

function _indel_frontier_assert_topology_controls(
        summary::DataFrames.DataFrame
)::Nothing
    branching = summary[
        (summary.graph_id .== "branching_k3_complete") .&
        (summary.window_bases .== 500),
        :,
    ]
    linear = summary[
        (summary.graph_id .== "linear_10001") .&
        (summary.window_bases .== 500),
        :,
    ]
    linear_reference = summary[
        (summary.graph_id .== "linear_2001") .&
        (summary.window_bases .== 500),
        :,
    ]
    DataFrames.nrow(branching) == 1 || error(
        "missing 500 bp highly-branching classifier control"
    )
    DataFrames.nrow(linear) == 1 || error(
        "missing 500 bp large-linear classifier control"
    )
    DataFrames.nrow(linear_reference) == 1 || error(
        "missing 500 bp linear runtime reference control"
    )
    only(branching.probe_reason) == "work_limit" || error(
        "frontier classifier failed to reject the small highly-branching graph"
    )
    only(linear.probe_reason) == "complete" || error(
        "frontier classifier failed to admit the 10,001-vertex linear graph"
    )
    only(linear.pair_hmm_valid) || error(
        "admitted 500 bp large-linear control did not use the pair-HMM path"
    )
    only(linear.pair_hmm_full_decode) || error(
        "admitted 500 bp large-linear control did not complete a full decode"
    )
    !only(linear.pair_hmm_truncated) || error(
        "admitted 500 bp large-linear control was truncated"
    )
    linear_p95_ms = only(linear.p95_ms)
    reference_p95_ms = only(linear_reference.p95_ms)
    isfinite(linear_p95_ms) && linear_p95_ms > 0 || error(
        "10,001-vertex/500 bp warmed p95 is not finite and positive: " *
        "$(linear_p95_ms) ms"
    )
    isfinite(reference_p95_ms) && reference_p95_ms > 0 || error(
        "2,001-vertex/500 bp reference p95 is not finite and positive: " *
        "$(reference_p95_ms) ms"
    )
    linear_p95_ms <= INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_MS || error(
        "10,001-vertex/500 bp warmed p95 $(linear_p95_ms) ms exceeds " *
        "the $(INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_MS) ms " *
        "affordability ceiling"
    )
    slowdown = linear_p95_ms / reference_p95_ms
    slowdown <= INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_SLOWDOWN || error(
        "10,001-vertex/500 bp warmed p95 slowdown $(slowdown)x exceeds " *
        "the same-run $(INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_SLOWDOWN)x " *
        "affordability ceiling"
    )
    return nothing
end

function _indel_frontier_write_figure(
        summary::DataFrames.DataFrame,
        output_dir::String,
)::NamedTuple
    figure = CairoMakie.Figure(
        size = (1_500, 650), fontsize = 14, backgroundcolor = :white
    )
    vertex_axis = CairoMakie.Axis(
        figure[1, 1],
        title = "Vertex count is not a full-decode runtime classifier",
        xlabel = "graph vertices",
        ylabel = "warmed pair-HMM p95 (ms)",
        xscale = log10,
        yscale = log10,
    )
    frontier_axis = CairoMakie.Axis(
        figure[1, 2],
        title = "Full-decode runtime against topology-only frontier work",
        xlabel = "bounded frontier work (area + edge expansions)",
        ylabel = "warmed pair-HMM p95 (ms)",
        xscale = log10,
        yscale = log10,
    )

    colors = Dict(
        "synthetic_branching" => :firebrick3,
        "synthetic_linear" => :seagreen4,
        "raw_fixed_toy" => :darkorange3,
        "cleaned_fixed_toy" => :dodgerblue3,
    )
    markers = Dict(250 => :circle, 500 => :rect, 1_200 => :diamond, 50 => :utriangle)
    vertex_midpoint = sqrt(minimum(summary.n_vertices) * maximum(summary.n_vertices))
    frontier_midpoint = sqrt(
        max(minimum(summary.probe_frontier_work), 1) *
        max(maximum(summary.probe_frontier_work), 1)
    )
    censored_index = 0
    for row in DataFrames.eachrow(summary)
        color = get(colors, row.graph_source, :gray40)
        marker = get(markers, row.window_bases, :circle)
        label = "$(row.graph_id):$(row.window_bases)"
        runtime_evidence_valid =
            row.pair_hmm_valid &&
            row.pair_hmm_full_decode &&
            !row.pair_hmm_truncated
        plotted_color = runtime_evidence_valid ? color : :gray45
        plotted_marker = runtime_evidence_valid ? marker : :xcross
        plotted_size = runtime_evidence_valid ? 13 : 16
        !runtime_evidence_valid && (censored_index += 1)
        label_vertical_offset = runtime_evidence_valid ?
                                5 : 5 + 12 * (censored_index - 1)
        status_label = if runtime_evidence_valid
            label
        elseif row.pair_hmm_truncated
            "$(label) [CENSORED: truncated diagnostic]"
        else
            "$(label) [CENSORED: incomplete diagnostic]"
        end
        vertex_on_left = row.n_vertices <= vertex_midpoint
        frontier_on_left = row.probe_frontier_work <= frontier_midpoint
        CairoMakie.scatter!(
            vertex_axis,
            [row.n_vertices],
            [row.p95_ms];
            color = plotted_color,
            marker = plotted_marker,
            markersize = plotted_size,
        )
        CairoMakie.text!(
            vertex_axis,
            row.n_vertices,
            row.p95_ms;
            text = status_label,
            fontsize = 8,
            align = (vertex_on_left ? :left : :right, :bottom),
            offset = (vertex_on_left ? 5 : -5, label_vertical_offset),
        )
        CairoMakie.scatter!(
            frontier_axis,
            [max(row.probe_frontier_work, 1)],
            [row.p95_ms];
            color = plotted_color,
            marker = plotted_marker,
            markersize = plotted_size,
        )
        CairoMakie.text!(
            frontier_axis,
            max(row.probe_frontier_work, 1),
            row.p95_ms;
            text = "$(status_label) [$(row.probe_reason)]",
            fontsize = 8,
            align = (frontier_on_left ? :left : :right, :bottom),
            offset = (frontier_on_left ? 5 : -5, label_vertical_offset),
        )
    end
    CairoMakie.Label(
        figure[0, 1:2],
        "Indel frontier scheduling calibration — runtime/topology only; no " *
        "correction accuracy";
        fontsize = 19,
        font = :bold,
    )
    CairoMakie.Label(
        figure[2, 1:2],
        "Colored symbols are complete, non-truncated pair-HMM decodes. Gray " *
        "× symbols are explicitly censored timing diagnostics and are not " *
        "full-decode runtime evidence.";
        fontsize = 12,
        color = :gray30,
    )

    png_path = joinpath(output_dir, "indel_frontier_runtime.png")
    svg_path = joinpath(output_dir, "indel_frontier_runtime.svg")
    CairoMakie.save(png_path, figure)
    CairoMakie.save(svg_path, figure)
    return (png = png_path, svg = svg_path)
end

function _indel_frontier_parse_args(args::Vector{String})::NamedTuple
    smoke = false
    output_dir = INDEL_FRONTIER_DEFAULT_OUTPUT_DIR
    repeats::Union{Int, Nothing} = nothing
    seen = Set{String}()
    index = 1
    while index <= length(args)
        flag = args[index]
        flag in seen && throw(ArgumentError("duplicate argument: $(flag)"))
        if flag == "--smoke"
            smoke = true
            push!(seen, flag)
            index += 1
        elseif flag == "--output-dir" || flag == "--repeats"
            push!(seen, flag)
            index == length(args) && throw(
                ArgumentError("$(flag) requires a value")
            )
            value = args[index + 1]
            (isempty(value) || startswith(value, "--")) && throw(
                ArgumentError("$(flag) requires a nonempty value")
            )
            if flag == "--output-dir"
                output_dir = value
            else
                parsed_repeats = tryparse(Int, value)
                isnothing(parsed_repeats) && throw(
                    ArgumentError("--repeats must be an integer, got $(value)")
                )
                parsed_repeats > 0 || throw(
                    ArgumentError("--repeats must be positive")
                )
                repeats = parsed_repeats
            end
            index += 2
        else
            throw(ArgumentError("unknown argument: $(flag)"))
        end
    end
    if !smoke && repeats !== nothing && repeats < INDEL_FRONTIER_REPEATS
        throw(
            ArgumentError(
                "non-smoke --repeats must be at least " *
                "$(INDEL_FRONTIER_REPEATS)"
            ),
        )
    end
    return (smoke = smoke, output_dir = output_dir, repeats = repeats)
end

function _indel_frontier_remove_prior_artifacts(
        output_dir::String,
        artifact_names::Tuple{Vararg{String}},
)::Nothing
    Base.Filesystem.mkpath(output_dir)
    # Invalidate the completion manifest before removing generation members.
    for artifact_index in length(artifact_names):-1:1
        artifact_name = artifact_names[artifact_index]
        artifact_path = joinpath(output_dir, artifact_name)
        isdir(artifact_path) && error(
            "refusing to replace artifact directory: $(artifact_path)"
        )
        Base.rm(artifact_path; force = true)
    end
    return nothing
end

function _indel_frontier_publish_artifacts(
        staging_dir::String,
        output_dir::String,
        artifact_names::Tuple{Vararg{String}},
)::Nothing
    for artifact_name in artifact_names
        staging_path = joinpath(staging_dir, artifact_name)
        isfile(staging_path) || error(
            "staged artifact is missing: $(staging_path)"
        )
    end
    for artifact_name in artifact_names
        Base.Filesystem.rename(
            joinpath(staging_dir, artifact_name),
            joinpath(output_dir, artifact_name),
        )
    end
    return nothing
end

function _indel_frontier_run_provenance(
        smoke::Bool,
        repeats::Int,
)::NamedTuple
    repository_root = normpath(joinpath(@__DIR__, ".."))
    git_head_sha = strip(
        Base.read(
            `git -C $repository_root rev-parse HEAD`,
            String,
        ),
    )
    tracked_diff = Base.read(
        `git -C $repository_root diff --binary --no-ext-diff HEAD --`
    )
    tracked_diff_sha256 = Base.bytes2hex(SHA.sha256(tracked_diff))
    benchmark_source_sha256 = _indel_frontier_file_sha256(@__FILE__)
    dependency = _indel_frontier_dependency_provenance()
    code_environment_components = (
        git_head_sha,
        tracked_diff_sha256,
        benchmark_source_sha256,
        dependency.project_toml_sha256,
        dependency.manifest_toml_sha256,
        string(VERSION),
        string(Threads.nthreads()),
        string(Sys.CPU_NAME),
        string(Sys.ARCH),
        string(Sys.KERNEL),
        string(Sys.CPU_THREADS),
    )
    code_environment_fingerprint = Base.bytes2hex(SHA.sha256(codeunits(join(
        code_environment_components, ":"
    ))))
    run_fingerprint = join(
        (
            code_environment_fingerprint,
            string(INDEL_FRONTIER_FIXTURE_SEED),
            string(smoke),
            string(repeats),
        ),
        ":",
    )
    generation_id = Base.bytes2hex(SHA.sha256(codeunits(run_fingerprint)))
    return (
        manifest_schema_version = 1,
        generation_id = generation_id,
        code_environment_fingerprint = code_environment_fingerprint,
        code_sha = git_head_sha,
        git_tracked_worktree_dirty = !isempty(tracked_diff),
        git_tracked_diff_sha256 = tracked_diff_sha256,
        benchmark_source_sha256 = benchmark_source_sha256,
        active_project_path = dependency.active_project_path,
        project_toml_sha256 = dependency.project_toml_sha256,
        manifest_toml_present = dependency.manifest_toml_present,
        manifest_toml_sha256 = dependency.manifest_toml_sha256,
        julia_version = string(VERSION),
        julia_threads = Threads.nthreads(),
        cpu_name = string(Sys.CPU_NAME),
        architecture = string(Sys.ARCH),
        kernel = string(Sys.KERNEL),
        cpu_threads = Sys.CPU_THREADS,
        smoke = smoke,
        repeats = repeats,
        synthetic_k = INDEL_FRONTIER_K,
        fixed_toy_k = INDEL_FRONTIER_TOY_K,
        genome_length = INDEL_FRONTIER_GENOME_LENGTH,
        source_read_length = INDEL_FRONTIER_SOURCE_READ_LENGTH,
        coverage = INDEL_FRONTIER_COVERAGE,
        error_rate = INDEL_FRONTIER_ERROR_RATE,
        fixture_seed = INDEL_FRONTIER_FIXTURE_SEED,
        window_bases = join(INDEL_FRONTIER_WINDOW_BASES, ";"),
        frontier_work_limit = Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT,
        linear_10k_500_max_p95_ms =
            INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_MS,
        linear_10k_500_max_p95_slowdown =
            INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_SLOWDOWN,
        accuracy_metric_used = false,
    )
end

function _indel_frontier_assert_provenance_unchanged(
        initial::NamedTuple,
        smoke::Bool,
        repeats::Int,
)::Nothing
    current = _indel_frontier_run_provenance(smoke, repeats)
    current.project_toml_sha256 == initial.project_toml_sha256 || error(
        "active Project.toml changed during the frontier runtime run; " *
        "refusing to publish artifacts"
    )
    current.manifest_toml_sha256 == initial.manifest_toml_sha256 || error(
        "active Manifest.toml changed during the frontier runtime run; " *
        "refusing to publish artifacts"
    )
    current.code_environment_fingerprint ==
    initial.code_environment_fingerprint || error(
        "code/worktree/environment fingerprint changed during the frontier " *
        "runtime run; refusing to publish artifacts"
    )
    return nothing
end

function _indel_frontier_dependency_provenance()::NamedTuple
    active_project = Base.active_project()
    active_project isa String || error(
        "an active Project.toml is required for benchmark provenance"
    )
    isfile(active_project) || error(
        "active Project.toml does not exist: $(active_project)"
    )
    manifest_path = joinpath(dirname(active_project), "Manifest.toml")
    manifest_present = isfile(manifest_path)
    return (
        active_project_path = normpath(active_project),
        project_toml_sha256 = _indel_frontier_file_sha256(active_project),
        manifest_toml_present = manifest_present,
        manifest_toml_sha256 = _indel_frontier_optional_dependency_sha256(
            manifest_path
        ),
    )
end

function _indel_frontier_optional_dependency_sha256(path::String)::String
    return isfile(path) ?
           _indel_frontier_file_sha256(path) :
           INDEL_FRONTIER_MISSING_DEPENDENCY_SENTINEL
end

function _indel_frontier_affordability_provenance(
        smoke::Bool,
        summary::DataFrames.DataFrame,
)::NamedTuple
    if smoke
        return (
            affordability_check_enforced = false,
            linear_10k_500_p95_ms = missing,
            linear_2001_500_p95_ms = missing,
            linear_10k_500_p95_slowdown = missing,
            linear_10k_500_affordable = missing,
        )
    end
    linear = summary[
        (summary.graph_id .== "linear_10001") .&
        (summary.window_bases .== 500),
        :,
    ]
    reference = summary[
        (summary.graph_id .== "linear_2001") .&
        (summary.window_bases .== 500),
        :,
    ]
    DataFrames.nrow(linear) == 1 || error(
        "missing 10,001-vertex/500 bp affordability row"
    )
    DataFrames.nrow(reference) == 1 || error(
        "missing 2,001-vertex/500 bp affordability reference row"
    )
    linear_p95_ms = only(linear.p95_ms)
    reference_p95_ms = only(reference.p95_ms)
    slowdown = linear_p95_ms / reference_p95_ms
    affordable =
        linear_p95_ms <= INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_MS &&
        slowdown <= INDEL_FRONTIER_LINEAR_10K_500_MAX_P95_SLOWDOWN
    return (
        affordability_check_enforced = true,
        linear_10k_500_p95_ms = linear_p95_ms,
        linear_2001_500_p95_ms = reference_p95_ms,
        linear_10k_500_p95_slowdown = slowdown,
        linear_10k_500_affordable = affordable,
    )
end

function _indel_frontier_file_sha256(path::String)::String
    return Base.bytes2hex(SHA.sha256(Base.read(path)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
