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
const INDEL_FRONTIER_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "td-jt7r-2-frontier-runtime"
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
    smoke = "--smoke" in args
    output_dir = _indel_frontier_arg_value(
        args, "--output-dir", INDEL_FRONTIER_DEFAULT_OUTPUT_DIR
    )
    repeats = parse(
        Int,
        _indel_frontier_arg_value(
            args, "--repeats", string(smoke ? 1 : INDEL_FRONTIER_REPEATS)
        ),
    )
    repeats > 0 || throw(ArgumentError("--repeats must be positive"))

    graph_cases = _indel_frontier_calibration_graphs(smoke)
    window_bases = smoke ? (50,) : INDEL_FRONTIER_WINDOW_BASES
    summary, replicates = _indel_frontier_run_matrix(
        graph_cases, window_bases, repeats
    )
    if !smoke
        _indel_frontier_assert_topology_controls(summary)
    end

    Base.Filesystem.mkpath(output_dir)
    summary_path = joinpath(output_dir, "indel_frontier_runtime_summary.csv")
    replicates_path = joinpath(output_dir, "indel_frontier_runtime_replicates.csv")
    CSV.write(summary_path, summary)
    CSV.write(replicates_path, replicates)

    figure_paths = _indel_frontier_write_figure(summary, output_dir)

    println("Wrote branching/frontier runtime calibration artifacts:")
    println("  summary:    $(summary_path)")
    println("  replicates: $(replicates_path)")
    println("  png:        $(figure_paths.png)")
    println("  svg:        $(figure_paths.svg)")
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
        truncated = truncated,
        truncated_samples = truncated_samples,
        complete = complete,
        complete_samples = complete_samples,
        full_decode = valid,
        replicates = rows,
    )
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
    truncated = get(diagnostics, :truncated, false)
    completed_columns = get(diagnostics, :completed_columns, 0)
    complete = !truncated && completed_columns == expected_columns
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
    DataFrames.nrow(branching) == 1 || error(
        "missing 500 bp highly-branching classifier control"
    )
    DataFrames.nrow(linear) == 1 || error(
        "missing 500 bp large-linear classifier control"
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

function _indel_frontier_arg_value(
        args::Vector{String},
        flag::String,
        default::String,
)::String
    index = findfirst(==(flag), args)
    if isnothing(index) || index == length(args)
        return default
    end
    return args[index + 1]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
