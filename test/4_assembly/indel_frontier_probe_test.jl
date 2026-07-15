# Score-free topology-frontier probe and pair-HMM runtime telemetry.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/indel_frontier_probe_test.jl")'

import BioSequences
import FASTX
import Graphs
import Kmers
import MetaGraphsNext
import Mycelia
import Test

const INDEL_FRONTIER_TASK_WAIT_SECONDS = 5.0

function indel_frontier_test_throws_message(
        callable::F,
        exception_type::Type{E},
        message_fragment::AbstractString,
)::Nothing where {F, E <: Exception}
    thrown::Union{Exception, Nothing} = nothing
    try
        callable()
    catch exception
        thrown = exception
    end
    Test.@test thrown isa E
    rendered = thrown isa Exception ?
               Base.sprint(Base.showerror, thrown) : ""
    Test.@test Base.occursin(message_fragment, rendered)
    return nothing
end

function indel_frontier_probe_graph()::MetaGraphsNext.MetaGraph
    records = [
        FASTX.FASTA.Record("read_1", BioSequences.dna"ATGCG"),
        FASTX.FASTA.Record("read_2", BioSequences.dna"GCGTA"),
        FASTX.FASTA.Record("read_3", BioSequences.dna"GTACC"),
        FASTX.FASTA.Record("read_4", BioSequences.dna"ACCGT"),
        FASTX.FASTA.Record("read_5", BioSequences.dna"CGTAA")
    ]
    return Mycelia.Rhizomorph.build_fasta_graph_olc(
        records; min_overlap = 3)
end

function indel_frontier_high_degree_graph(
        degree::Int,
)::MetaGraphsNext.MetaGraph
    degree > 0 || throw(ArgumentError("degree must be positive"))
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = Any,
        edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
        weight_function = Mycelia.Rhizomorph.edge_data_weight,
        default_weight = 0.0,
    )
    graph["root"] = nothing
    for index in 1:degree
        target = "target_$index"
        graph[target] = nothing
        graph["root", target] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            1.0,
            Mycelia.Rhizomorph.Forward,
            Mycelia.Rhizomorph.Forward,
        )
    end
    return graph
end


Test.@testset "Checkpoint restores complete per-pass indel telemetry" begin
    serialized = Dict{String, Any}(
        "indel_rung_telemetry" => Any[
            Dict{String, Any}(
                "profile_requested" => true,
                "requested" => 5,
                "attempted" => 3,
                "completed" => 2,
                "truncated" => 1,
                "engaged" => 1,
                "graph_source" => "mixed",
                "decision_reason" => "sparse_frontier_affordable",
                "raw_frontier_metrics" => Any[
                    Dict{String, Any}(
                        "anchored" => true,
                        "reason" => "complete",
                    ),
                ],
            ),
            Dict{String, Any}(
                "profile_requested" => true,
                "requested" => 2,
                "attempted" => 0,
                "completed" => 0,
                "truncated" => 0,
                "engaged" => 0,
                "decision_reason" => "frontier_budget_exceeded",
            ),
        ],
    )
    restored = Mycelia._restore_indel_rung_telemetry(serialized)

    Test.@test length(restored) == 2
    Test.@test restored[1][:requested] == 5
    Test.@test restored[1][:attempted] == 3
    Test.@test restored[1][:completed] == 2
    Test.@test restored[1][:truncated] == 1
    Test.@test restored[1][:engaged] == 1
    Test.@test restored[1][:graph_source] == :mixed
    Test.@test restored[1][:decision_reason] == :sparse_frontier_affordable
    Test.@test only(restored[1][:raw_frontier_metrics])[:reason] == :complete
    Test.@test restored[2][:decision_reason] == :frontier_budget_exceeded
    Test.@test !restored[2][:bounded_windowing_forced]

    missing_outcome = Base.deepcopy(serialized)
    first(missing_outcome["indel_rung_telemetry"])["truncated"] = 0
    indel_frontier_test_throws_message(
        () -> Mycelia._restore_indel_rung_telemetry(missing_outcome),
        ArgumentError,
        "checkpoint indel telemetry completed + truncated must equal attempted",
    )
end

function indel_frontier_probe_config(;
        band_width::Union{Nothing, Int} = 4,
        target_vertex::Any = nothing,
)::Mycelia.ViterbiCorrectionConfig
    return Mycelia.ViterbiCorrectionConfig(
        error_rate = 0.10,
        strand_mode = :singlestrand,
        indel_moves = true,
        insertion_fraction = 0.30,
        deletion_fraction = 0.30,
        insertion_extend_probability = 0.10,
        deletion_extend_probability = 0.10,
        deletion_max_run = 3,
        max_insertion_run = 3,
        band_width = band_width,
        beam_width = typemax(Int),
        target_vertex = target_vertex,
    )
end

Test.@testset "Initial deletion relaxation honors the adaptive band" begin
    graph = indel_frontier_probe_graph()
    observation = [BioSequences.dna"ATGCG"]
    zero_band_config = indel_frontier_probe_config(; band_width = 0)
    unbounded_config = indel_frontier_probe_config(; band_width = nothing)

    zero_band_metrics = Mycelia._probe_indel_frontier(
        graph,
        observation,
        :DNA;
        config = zero_band_config,
        strand_mode = :singlestrand,
    )
    unbounded_metrics = Mycelia._probe_indel_frontier(
        graph,
        observation,
        :DNA;
        config = unbounded_config,
        strand_mode = :singlestrand,
    )

    Test.@test zero_band_metrics.reason == :complete
    Test.@test zero_band_metrics.completed_columns == 1
    Test.@test zero_band_metrics.frontier_area == 1
    Test.@test zero_band_metrics.peak_frontier == 1
    Test.@test unbounded_metrics.reason == :complete
    Test.@test unbounded_metrics.completed_columns == 1
    Test.@test unbounded_metrics.frontier_area > zero_band_metrics.frontier_area
    Test.@test unbounded_metrics.peak_frontier > zero_band_metrics.peak_frontier

    zero_band_result = only(Mycelia.correct_observations(
        graph,
        [observation];
        config = zero_band_config,
    ).paths)
    unbounded_result = only(Mycelia.correct_observations(
        graph,
        [observation];
        config = unbounded_config,
    ).paths)

    Test.@test zero_band_result.diagnostics[:completed_columns] == 1
    Test.@test zero_band_result.diagnostics[:frontier_area] == 1
    Test.@test zero_band_result.diagnostics[:peak_frontier] == 1
    Test.@test unbounded_result.diagnostics[:completed_columns] == 1
    Test.@test unbounded_result.diagnostics[:frontier_area] >
               zero_band_result.diagnostics[:frontier_area]
    Test.@test unbounded_result.diagnostics[:peak_frontier] >
               zero_band_result.diagnostics[:peak_frontier]

    target_vertex = BioSequences.dna"GCGTA"
    zero_band_target = Mycelia.correct_observations(
        graph,
        [observation];
        config = indel_frontier_probe_config(
            band_width = 0,
            target_vertex = target_vertex,
        ),
    )
    one_band_target = Mycelia.correct_observations(
        graph,
        [observation];
        config = indel_frontier_probe_config(
            band_width = 1,
            target_vertex = target_vertex,
        ),
    )
    unbounded_target = Mycelia.correct_observations(
        graph,
        [observation];
        config = indel_frontier_probe_config(
            band_width = nothing,
            target_vertex = target_vertex,
        ),
    )

    Test.@test only(zero_band_target.corrected_observations) === nothing
    for result in (one_band_target, unbounded_target)
        path_result = only(result.paths)
        Test.@test [string(step.vertex_label) for step in path_result.path.steps] ==
                   ["ATGCG", "GCGTA"]
        Test.@test path_result.diagnostics[:move_counts][:D] == 1
    end
end

Test.@testset "Indel topology-frontier probe" begin
    graph = indel_frontier_probe_graph()
    config = indel_frontier_probe_config()
    observation = [
        BioSequences.dna"ATGCG",
        BioSequences.dna"GCGTA",
        BioSequences.dna"GTACC"
    ]

    metrics = Mycelia._probe_indel_frontier(
        graph,
        observation,
        :DNA;
        config = config,
        strand_mode = :singlestrand
    )
    Test.@test metrics isa Mycelia.IndelFrontierMetrics
    Test.@test metrics.anchored
    Test.@test metrics.reason == :complete
    Test.@test metrics.window_length == length(observation)
    Test.@test metrics.vertex_count == 5
    Test.@test metrics.edge_count == 4
    Test.@test metrics.branch_vertices == 0
    Test.@test metrics.join_vertices == 0
    Test.@test metrics.branch_fraction == 0.0
    Test.@test metrics.join_fraction == 0.0
    Test.@test metrics.max_out_degree == 1
    Test.@test metrics.frontier_area >= metrics.peak_frontier >= 1
    Test.@test metrics.edge_expansions > 0
    Test.@test metrics.completed_columns == length(observation)
    Test.@test metrics.frontier_work ==
               metrics.frontier_area + metrics.edge_expansions

    graph_summary = Mycelia._indel_frontier_graph_summary(graph)
    Test.@test graph_summary.successor_index isa
               Mycelia._IndelFrontierSuccessorIndex
    cached_metrics = Mycelia._probe_indel_frontier(
        graph,
        observation,
        :DNA;
        config = config,
        strand_mode = :singlestrand,
        graph_summary = graph_summary,
    )
    for field in fieldnames(Mycelia.IndelFrontierMetrics)
        Test.@test getfield(cached_metrics, field) == getfield(metrics, field)
    end

    unanchored = Mycelia._probe_indel_frontier(
        graph,
        [BioSequences.dna"AAAAA", BioSequences.dna"AAAAC"],
        :DNA;
        config = config,
        strand_mode = :singlestrand
    )
    Test.@test !unanchored.anchored
    Test.@test unanchored.reason == :unanchored_start
    Test.@test unanchored.completed_columns == 0
    Test.@test unanchored.frontier_work == 0

    limited = Mycelia._probe_indel_frontier(
        graph,
        observation,
        :DNA;
        config = config,
        strand_mode = :singlestrand,
        work_limit = 0
    )
    Test.@test limited.anchored
    Test.@test limited.reason == :work_limit
    Test.@test limited.frontier_work > 0
    Test.@test limited.completed_columns < length(observation)

    caught = nothing
    try
        Mycelia._probe_indel_frontier(
            graph,
            observation,
            :DNA;
            config = config,
            strand_mode = :singlestrand,
            work_limit = -1
        )
    catch error
        caught = error
    end
    Test.@test caught isa ArgumentError
    Test.@test occursin("work_limit must be non-negative", sprint(showerror, caught))
end

Test.@testset "Frontier successor topology is indexed once" begin
    graph = indel_frontier_probe_graph()
    resolved_edges = Ref(0)
    strand_resolver = function (edge_data::Any)
        resolved_edges[] += 1
        return Mycelia._indel_frontier_edge_strands(edge_data)
    end
    successor_index = Mycelia._indel_frontier_successor_index(
        graph; strand_resolver = strand_resolver)
    Test.@test resolved_edges[] == Graphs.ne(graph.graph)

    source_vertex = first(first(MetaGraphsNext.edge_labels(graph)))
    cell = Mycelia._IndelFrontierCell(
        source_vertex,
        Mycelia.Rhizomorph.Forward,
        :M,
        0,
        0,
    )
    expected = Mycelia._indel_frontier_successors(
        successor_index, cell, typemax(Int))
    for _ in 1:10
        observed = Mycelia._indel_frontier_successors(
            successor_index, cell, typemax(Int))
        Test.@test observed.successors == expected.successors
        Test.@test observed.overflowed == expected.overflowed
    end
    Mycelia._indel_frontier_start_strands(
        successor_index,
        source_vertex,
        :singlestrand,
        Mycelia.Rhizomorph.Forward,
    )
    Test.@test resolved_edges[] == Graphs.ne(graph.graph)
end

Test.@testset "Undirected frontier topology matches weighted decode" begin
    graph = MetaGraphsNext.MetaGraph(
        Graphs.Graph();
        label_type = String,
        vertex_data_type = Any,
        edge_data_type = Mycelia.Rhizomorph.KmerEdgeData,
    )
    graph["left"] = nothing
    graph["right"] = nothing
    graph["left", "right"] = Mycelia.Rhizomorph.KmerEdgeData()

    stored_source, stored_target = only(MetaGraphsNext.edge_labels(graph))
    resolved_edges = Ref(0)
    strand_resolver = function (edge_data::Any)
        resolved_edges[] += 1
        return Mycelia._indel_frontier_edge_strands(edge_data)
    end
    raw_index = Mycelia._indel_frontier_successor_index(
        graph; strand_resolver = strand_resolver)
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
    weighted_index = Mycelia._indel_frontier_successor_index(weighted)

    for (source, target) in (
            (stored_source, stored_target),
            (stored_target, stored_source),
    )
        cell = Mycelia._IndelFrontierCell(
            source,
            Mycelia.Rhizomorph.Forward,
            :M,
            0,
            0,
        )
        raw_successors = Mycelia._indel_frontier_successors(
            raw_index, cell, typemax(Int))
        weighted_successors = Mycelia._indel_frontier_successors(
            weighted_index, cell, typemax(Int))
        Test.@test raw_successors.successors ==
                   [(target, Mycelia.Rhizomorph.Forward)]
        Test.@test raw_successors.successors == weighted_successors.successors
        Test.@test !raw_successors.overflowed
        Test.@test !weighted_successors.overflowed
    end
    Test.@test resolved_edges[] == Graphs.ne(graph.graph)

    self_loop_graph = MetaGraphsNext.MetaGraph(
        Graphs.Graph();
        label_type = String,
        vertex_data_type = Any,
        edge_data_type = Mycelia.Rhizomorph.KmerEdgeData,
    )
    self_loop_graph["loop"] = nothing
    self_loop_graph["loop", "loop"] = Mycelia.Rhizomorph.KmerEdgeData()
    self_loop_resolutions = Ref(0)
    self_loop_resolver = function (edge_data::Any)
        self_loop_resolutions[] += 1
        return Mycelia._indel_frontier_edge_strands(edge_data)
    end
    self_loop_index = Mycelia._indel_frontier_successor_index(
        self_loop_graph; strand_resolver = self_loop_resolver)
    self_loop_cell = Mycelia._IndelFrontierCell(
        "loop",
        Mycelia.Rhizomorph.Forward,
        :M,
        0,
        0,
    )
    self_loop_successors = Mycelia._indel_frontier_successors(
        self_loop_index, self_loop_cell, typemax(Int))
    self_loop_weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(
        self_loop_graph)
    self_loop_weighted_index = Mycelia._indel_frontier_successor_index(
        self_loop_weighted)
    self_loop_weighted_successors = Mycelia._indel_frontier_successors(
        self_loop_weighted_index, self_loop_cell, typemax(Int))
    Test.@test self_loop_successors.successors ==
               [("loop", Mycelia.Rhizomorph.Forward)]
    Test.@test self_loop_successors.successors ==
               self_loop_weighted_successors.successors
    Test.@test self_loop_resolutions[] == Graphs.ne(self_loop_graph.graph)
end

Test.@testset "Weighted undirected frontier preserves stored direction" begin
    graph = MetaGraphsNext.MetaGraph(
        Graphs.Graph();
        label_type = String,
        vertex_data_type = Any,
        edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
        weight_function = Mycelia.Rhizomorph.edge_data_weight,
        default_weight = 0.0,
    )
    graph["left"] = nothing
    graph["right"] = nothing
    graph["left", "right"] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
        1.0,
        Mycelia.Rhizomorph.Forward,
        Mycelia.Rhizomorph.Reverse,
    )
    stored_source, stored_target = only(MetaGraphsNext.edge_labels(graph))
    Test.@test Mycelia.build_correction_weighted_graph(graph) === graph
    successor_index = Mycelia._indel_frontier_successor_index(graph)
    decode_cache = Dict{
        Tuple{String, Mycelia.Rhizomorph.StrandOrientation},
        Mycelia._IndelDecodeSuccessorBatch{String},
    }()

    for (vertex, strand) in (
            (stored_source, Mycelia.Rhizomorph.Forward),
            (stored_source, Mycelia.Rhizomorph.Reverse),
            (stored_target, Mycelia.Rhizomorph.Forward),
            (stored_target, Mycelia.Rhizomorph.Reverse),
    )
        cell = Mycelia._IndelFrontierCell(vertex, strand, :M, 0, 0)
        frontier_successors = Mycelia._indel_frontier_successors(
            successor_index, cell, typemax(Int))
        decoder_batch = Mycelia._indel_decode_successors!(
            decode_cache, graph, vertex, strand)
        decoder_successors = [
            successor for (successor, _) in decoder_batch.successors
        ]
        Test.@test frontier_successors.successors == decoder_successors
        Test.@test !frontier_successors.overflowed
    end

    source_cell = Mycelia._IndelFrontierCell(
        stored_source,
        Mycelia.Rhizomorph.Forward,
        :M,
        0,
        0,
    )
    source_successors = Mycelia._indel_frontier_successors(
        successor_index, source_cell, typemax(Int))
    Test.@test source_successors.successors ==
               [(stored_target, Mycelia.Rhizomorph.Reverse)]
    Test.@test !haskey(successor_index.outgoing, stored_target)
end

Test.@testset "High-degree frontier expansion is work-bounded" begin
    degree = 10_000
    graph = indel_frontier_high_degree_graph(degree)
    successor_index = Mycelia._indel_frontier_successor_index(graph)
    cell = Mycelia._IndelFrontierCell(
        "root",
        Mycelia.Rhizomorph.Forward,
        :M,
        0,
        0,
    )

    bounded = Mycelia._indel_frontier_successors(successor_index, cell, 3)
    Test.@test bounded.overflowed
    Test.@test length(bounded.successors) == 3

    unrestricted = Mycelia._indel_frontier_successors(
        successor_index, cell, typemax(Int))
    Test.@test !unrestricted.overflowed
    Test.@test length(unrestricted.successors) == degree

    metrics = Mycelia._probe_indel_frontier(
        graph,
        ["root", "target_1"],
        :TOKEN;
        config = indel_frontier_probe_config(),
        strand_mode = :singlestrand,
        work_limit = 3,
    )
    Test.@test metrics.reason == :work_limit
    Test.@test metrics.completed_columns == 0
    Test.@test metrics.edge_expansions == 4
    Test.@test metrics.frontier_work == 4
end

Test.@testset "Pair-HMM frontier telemetry is observational" begin
    graph = indel_frontier_probe_graph()
    config = indel_frontier_probe_config()
    observation = [BioSequences.dna"ATGCG", BioSequences.dna"GTACC"]

    result = Mycelia.correct_observations(
        graph, [observation]; config = config)
    path_result = only(result.paths)
    diagnostics = path_result.diagnostics

    Test.@test diagnostics[:algorithm] == :viterbi_indel_pair_hmm
    Test.@test diagnostics[:frontier_area] >= diagnostics[:peak_frontier] >= 1
    Test.@test diagnostics[:edge_expansions] > 0
    Test.@test diagnostics[:completed_columns] == length(observation)
    Test.@test [string(step.vertex_label) for step in path_result.path.steps] ==
               ["ATGCG", "GCGTA", "GTACC"]
end


function indel_classifier_de_bruijn_sequence(order::Int)::String
    order > 0 || throw(ArgumentError("de Bruijn order must be positive"))
    alphabet = collect("ACGT")
    alphabet_size = length(alphabet)
    work = zeros(Int, alphabet_size * order + 1)
    indices = Int[]

    function visit(t::Int, period::Int)::Nothing
        if t > order
            order % period == 0 && append!(indices, @view work[2:(period + 1)])
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

function indel_classifier_fastq(
        identifier::String,
        sequence::String,
)::FASTX.FASTQ.Record
    return FASTX.FASTQ.Record(identifier, sequence, repeat("I", length(sequence)))
end

function indel_classifier_observations(
        read::FASTX.FASTQ.Record,
        k::Int,
)::Vector{Mycelia.QualityObservation}
    sequence = Mycelia.extract_typed_sequence(
        read, BioSequences.LongDNA{4})
    kmers = collect(Mycelia._record_kmer_iterator(BioSequences.LongDNA{4}, k, sequence))
    quality_scores = collect(FASTX.quality_scores(read))
    observations = Vector{Mycelia.QualityObservation}(undef, length(kmers))
    for (index, kmer) in enumerate(kmers)
        observations[index] = Mycelia.QualityObservation(
            kmer, UInt8.(@view quality_scores[index:(index + k - 1)]))
    end
    return observations
end

function indel_classifier_params()::Mycelia.IndelDecodeParams
    profile = Mycelia.indel_error_profile(:nanopore)
    return Mycelia.IndelDecodeParams(
        profile.base_error_rate,
        profile.insertion_fraction,
        profile.deletion_fraction,
        profile.insertion_extend_probability,
        profile.deletion_extend_probability,
        3,
        3,
        16,
    )
end

mutable struct IndelOOMWeightedGraphBuilder
    calls::Int
end

IndelOOMWeightedGraphBuilder()::IndelOOMWeightedGraphBuilder =
    IndelOOMWeightedGraphBuilder(0)

function (builder::IndelOOMWeightedGraphBuilder)(
        ::MetaGraphsNext.MetaGraph,
)::Any
    builder.calls += 1
    throw(OutOfMemoryError())
end

mutable struct IndelOOMGraphCleaner
    calls::Int
end

IndelOOMGraphCleaner()::IndelOOMGraphCleaner = IndelOOMGraphCleaner(0)

function (cleaner::IndelOOMGraphCleaner)(
        ::MetaGraphsNext.MetaGraph,
        ::Int,
)::Dict{String, Any}
    cleaner.calls += 1
    throw(OutOfMemoryError())
end

Test.@testset "Frontier classifier separates branching from graph size" begin
    alphabet = collect("ACGT")
    branch_records = FASTX.FASTQ.Record[]
    for a in alphabet, b in alphabet, c in alphabet, d in alphabet
        sequence = string(a, b, c, d)
        identifier = "branch_$(length(branch_records) + 1)"
        push!(branch_records, indel_classifier_fastq(identifier, sequence))
    end
    branch_graph = Mycelia.Rhizomorph.build_qualmer_graph(
        branch_records, 3; mode = :singlestrand)
    branch_cycle = indel_classifier_de_bruijn_sequence(3)
    branch_sequence = first(repeat(branch_cycle, cld(500, length(branch_cycle))), 500)
    branch_read = indel_classifier_fastq("branch_probe", branch_sequence)

    linear_cycle = indel_classifier_de_bruijn_sequence(7)
    linear_sequence = first(linear_cycle, 10_007)
    linear_record = indel_classifier_fastq("linear_10001", linear_sequence)
    linear_graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[linear_record], 7; mode = :singlestrand)
    linear_read = indel_classifier_fastq("linear_probe", first(linear_sequence, 500))

    config = Mycelia._indel_probe_config(
        indel_classifier_params(), :singlestrand)
    branch_metrics = Mycelia._probe_indel_frontier(
        branch_graph,
        indel_classifier_observations(branch_read, 3),
        :DNA;
        config = config,
        strand_mode = :singlestrand,
        work_limit = Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT,
    )
    linear_metrics = Mycelia._probe_indel_frontier(
        linear_graph,
        indel_classifier_observations(linear_read, 7),
        :DNA;
        config = config,
        strand_mode = :singlestrand,
        work_limit = Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT,
    )

    Test.@test Graphs.nv(branch_graph.graph) == 64
    Test.@test branch_metrics.max_out_degree == 4
    Test.@test branch_metrics.reason == :work_limit
    Test.@test !Mycelia._indel_frontier_admitted(
        branch_metrics, Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT)
    Test.@test Graphs.nv(linear_graph.graph) == 10_001
    Test.@test linear_metrics.max_out_degree == 1
    Test.@test linear_metrics.reason == :complete
    Test.@test Mycelia._indel_frontier_admitted(
        linear_metrics, Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT)

    # Exercise the production scheduler, including private-copy cleaning. The
    # branch graph must be rejected without mutating its reusable raw graph; the
    # much larger linear graph must be assigned to the raw pair-HMM source.
    branch_labels_before = sort!(string.(collect(MetaGraphsNext.labels(branch_graph))))
    branch_edges_before = sort!(string.(collect(MetaGraphsNext.edge_labels(branch_graph))))
    branch_hard = Set(MetaGraphsNext.labels(branch_graph))
    prior_soft_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
    branch_edge_labels = collect(MetaGraphsNext.edge_labels(branch_graph))
    Mycelia.Rhizomorph.accumulate_path_probability!(
        prior_soft_weights, branch_edge_labels, 0.25)
    Mycelia.Rhizomorph.clear_soft_edge_weights!()
    branch_decision = nothing
    try
        Mycelia.Rhizomorph.register_soft_edge_weights!(
            branch_graph, prior_soft_weights)
        registered_weights = Base.IdDict{Any, Float64}()
        for edge_label in branch_edge_labels
            edge = branch_graph[edge_label...]
            registered_weights[edge] =
                Mycelia.Rhizomorph.compute_edge_weight(edge)
        end
        Test.@test all(isapprox(weight, 0.25) for weight in values(registered_weights))
        branch_decision = Mycelia._evaluate_indel_frontier_schedule(
            FASTX.FASTQ.Record[branch_read],
            branch_graph,
            3,
            branch_hard,
            Bool[false],
            indel_classifier_params(),
            :singlestrand;
            prior_soft_weights = prior_soft_weights,
        )
        Test.@test all(
            Mycelia.Rhizomorph.compute_edge_weight(edge) == weight
            for (edge, weight) in registered_weights
        )
    finally
        Mycelia.Rhizomorph.clear_soft_edge_weights!()
    end
    Test.@test branch_decision.requested == 1
    Test.@test !branch_decision.admitted
    Test.@test branch_decision.admitted_windows == 0
    Test.@test branch_decision.rejected_windows == 1
    Test.@test !isempty(branch_decision.cleanup)
    Test.@test only(values(branch_decision.window_sources[1])) == :substitution
    Test.@test sort!(string.(collect(MetaGraphsNext.labels(branch_graph)))) ==
               branch_labels_before
    Test.@test sort!(string.(collect(MetaGraphsNext.edge_labels(branch_graph)))) ==
               branch_edges_before

    linear_hard = Set(MetaGraphsNext.labels(linear_graph))
    linear_decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[linear_record],
        linear_graph,
        7,
        linear_hard,
        Bool[false],
        indel_classifier_params(),
        :singlestrand,
    )
    max_window = Mycelia._INDEL_FRONTIER_MAX_WINDOW
    expected_linear_windows = 1 + cld(
        length(linear_sequence) - max_window, max_window - 7)
    Test.@test linear_decision.requested == expected_linear_windows
    Test.@test linear_decision.admitted
    Test.@test linear_decision.admitted_windows == expected_linear_windows
    Test.@test linear_decision.rejected_windows == 0
    Test.@test length(linear_decision.window_sources[1]) ==
               expected_linear_windows
    Test.@test all(
        source == :raw
        for source in values(linear_decision.window_sources[1])
    )

    # A 1 kb read remains below the normal whole-read windowing threshold, but
    # the frontier-budgeted schedule must still retain its rejected substitution
    # windows and force bounded windowing. Anchored overlaps require three windows;
    # the fully branching graph makes every one unaffordable without relying on
    # correction accuracy.
    short_sequence = first(
        repeat(branch_cycle, cld(1_000, length(branch_cycle))), 1_000)
    short_read = indel_classifier_fastq("all_rejected_short", short_sequence)
    Test.@test !Mycelia._windowed_decode_read_is_long(short_read, 3)
    short_decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[short_read],
        branch_graph,
        3,
        branch_hard,
        Bool[false],
        indel_classifier_params(),
        :singlestrand,
    )
    expected_short_windows = 1 + cld(
        length(short_sequence) - max_window, max_window - 3)
    Test.@test short_decision.requested == expected_short_windows
    Test.@test short_decision.admitted_windows == 0
    Test.@test short_decision.rejected_windows == expected_short_windows
    Test.@test !short_decision.admitted
    Test.@test Set(keys(short_decision.window_sources)) == Set([1])
    Test.@test length(short_decision.window_sources[1]) == expected_short_windows
    Test.@test all(
        source == :substitution
        for source in values(short_decision.window_sources[1])
    )
    Test.@test all(
        length(window) <= Mycelia._INDEL_FRONTIER_MAX_WINDOW
        for window in keys(short_decision.window_sources[1])
    )

    # The all-rejected schedule builds no weighted graph itself. Force the later
    # substitution pass build to fail and verify that the production recovery
    # branch gates the pass without confusing this allocator failure with a low-k
    # density/floor gate.
    pass_oom_builder = IndelOOMWeightedGraphBuilder()
    pass_oom_telemetry = Dict{Symbol, Any}()
    pass_oom_result = Mycelia.improve_read_set_likelihood(
        FASTX.FASTQ.Record[short_read],
        branch_graph,
        3;
        graph_mode = :singlestrand,
        beam_width = 16,
        hard_vertices = branch_hard,
        windowed_decode = true,
        indel_params = indel_classifier_params(),
        indel_schedule = :frontier_budgeted,
        rung_telemetry = pass_oom_telemetry,
        weighted_graph_builder = pass_oom_builder,
    )
    pass_oom_read = only(pass_oom_result[1])
    Test.@test pass_oom_builder.calls == 1
    Test.@test !pass_oom_result[5]
    Test.@test FASTX.sequence(String, pass_oom_read) == short_sequence
    Test.@test String(FASTX.quality(pass_oom_read)) ==
               String(FASTX.quality(short_read))
    Test.@test pass_oom_telemetry[:decision_reason] ==
               :frontier_budget_exceeded
    Test.@test pass_oom_telemetry[:bounded_windowing_forced]
    Test.@test pass_oom_telemetry[:substitution_decode_memory_gated]
    Test.@test pass_oom_telemetry[:requested] == expected_short_windows
    Test.@test pass_oom_telemetry[:attempted] == 0
    Test.@test pass_oom_telemetry[:completed] == 0
    Test.@test pass_oom_telemetry[:truncated] == 0
    Test.@test pass_oom_telemetry[:engaged] == 0

    # A direct caller owns the complete soft-snapshot lifetime, including the
    # cleaned-graph scheduler and post-schedule errors. Preserve an existing,
    # nonempty task-local snapshot exactly rather than merely clearing it on exit.
    baseline_soft_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
    Mycelia.Rhizomorph.accumulate_path_probability!(
        baseline_soft_weights, branch_edge_labels, 0.125)
    Mycelia.Rhizomorph.clear_soft_edge_weights!()
    Mycelia.Rhizomorph.register_soft_edge_weights!(
        branch_graph,
        baseline_soft_weights;
        min_support = typemax(Int),
    )
    baseline_weights = Base.IdDict{Any, Float64}()
    for edge_label in branch_edge_labels
        edge = branch_graph[edge_label...]
        baseline_weights[edge] =
            Mycelia.Rhizomorph.compute_edge_weight(edge)
    end
    telemetry = Dict{Symbol, Any}()
    try
        result = Mycelia.improve_read_set_likelihood(
            FASTX.FASTQ.Record[short_read],
            branch_graph,
            3;
            graph_mode = :singlestrand,
            beam_width = 16,
            prior_soft_weights = prior_soft_weights,
            hard_vertices = branch_hard,
            windowed_decode = true,
            indel_params = indel_classifier_params(),
            indel_schedule = :frontier_budgeted,
            rung_telemetry = telemetry,
        )
        Test.@test !result[5]
        Test.@test telemetry[:requested] == expected_short_windows
        Test.@test telemetry[:admitted_windows] == 0
        Test.@test telemetry[:rejected_windows] == expected_short_windows
        Test.@test telemetry[:bounded_windowing_forced]
        Test.@test telemetry[:attempted] == 0
        Test.@test all(
            Mycelia.Rhizomorph.compute_edge_weight(edge) == weight
            for (edge, weight) in baseline_weights
        )

        error_weights_before = copy(baseline_weights)
        caught = nothing
        try
            Mycelia.improve_read_set_likelihood(
                FASTX.FASTQ.Record[branch_read],
                branch_graph,
                3;
                batch_size = 0,
                graph_mode = :singlestrand,
                prior_soft_weights = prior_soft_weights,
                hard_vertices = branch_hard,
                windowed_decode = true,
                decode_gate_density = 0.0,
                indel_params = indel_classifier_params(),
                indel_schedule = :frontier_budgeted,
            )
        catch error
            caught = error
        end
        Test.@test caught isa ArgumentError
        Test.@test all(
            Mycelia.Rhizomorph.compute_edge_weight(edge) == weight
            for (edge, weight) in error_weights_before
        )
    finally
        Mycelia.Rhizomorph.clear_soft_edge_weights!()
    end

    # A soft scope is local to its owning task. Hold one open while an unrelated
    # unrestricted task enters concurrently: the owner must retain 0.25 while the
    # unrestricted task and this test task both see raw coverage immediately.
    # This catches both process-global leakage and parent/child lock deadlock.
    concurrent_edge_label = first(branch_edge_labels)
    concurrent_edge = branch_graph[concurrent_edge_label...]
    raw_edge_weight = Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge)
    concurrent_soft_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
    Mycelia.Rhizomorph.accumulate_path_probability!(
        concurrent_soft_weights, [concurrent_edge_label], 0.25)
    soft_scope_entered = Base.Channel{Nothing}(1)
    soft_scope_weight = Base.Channel{Float64}(1)
    release_soft_scope = Base.Channel{Nothing}(1)
    unrestricted_started = Base.Channel{Nothing}(1)
    unrestricted_entered = Base.Channel{Nothing}(1)

    Mycelia.Rhizomorph.clear_soft_edge_weights!()
    soft_task = Base.Threads.@spawn begin
        Mycelia.Rhizomorph._with_soft_edge_weight_scope(
            branch_graph, concurrent_soft_weights) do
            Base.put!(
                soft_scope_weight,
                Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge),
            )
            Base.put!(soft_scope_entered, nothing)
            Base.take!(release_soft_scope)
            return Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge)
        end
    end
    Base.take!(soft_scope_entered)
    Test.@test Base.take!(soft_scope_weight) == 0.25
    Test.@test Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge) ==
               raw_edge_weight

    unrestricted_task = Base.Threads.@spawn begin
        Base.put!(unrestricted_started, nothing)
        return Mycelia.Rhizomorph._with_soft_edge_weight_scope(
            branch_graph, nothing) do
            Base.put!(unrestricted_entered, nothing)
            return Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge)
        end
    end
    Base.take!(unrestricted_started)
    unrestricted_entry_status = try
        Base.timedwait(
            () -> Base.isready(unrestricted_entered),
            INDEL_FRONTIER_TASK_WAIT_SECONDS,
        )
    finally
        Base.put!(release_soft_scope, nothing)
    end
    soft_edge_weight = Base.fetch(soft_task)
    unrestricted_edge_weight = Base.fetch(unrestricted_task)

    Test.@test unrestricted_entry_status == :ok
    Test.@test soft_edge_weight == 0.25
    Test.@test unrestricted_edge_weight == raw_edge_weight
    Test.@test Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge) ==
               raw_edge_weight

    # Two soft owners may overlap without clobbering one another: each task sees
    # its own immutable snapshot while the unrelated test task remains raw.
    alternate_soft_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
    Mycelia.Rhizomorph.accumulate_path_probability!(
        alternate_soft_weights, [concurrent_edge_label], 0.5)
    first_owner_entered = Base.Channel{Nothing}(1)
    second_owner_entered = Base.Channel{Nothing}(1)
    release_first_owner = Base.Channel{Nothing}(1)
    release_second_owner = Base.Channel{Nothing}(1)
    first_owner = Base.Threads.@spawn begin
        Mycelia.Rhizomorph._with_soft_edge_weight_scope(
            branch_graph, concurrent_soft_weights) do
            Base.put!(first_owner_entered, nothing)
            Base.take!(release_first_owner)
            return Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge)
        end
    end
    second_owner = Base.Threads.@spawn begin
        Mycelia.Rhizomorph._with_soft_edge_weight_scope(
            branch_graph, alternate_soft_weights) do
            Base.put!(second_owner_entered, nothing)
            Base.take!(release_second_owner)
            return Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge)
        end
    end
    first_owner_status = nothing
    second_owner_status = nothing
    unrelated_owner_weight = NaN
    try
        first_owner_status =
            Base.timedwait(
                () -> Base.isready(first_owner_entered),
                INDEL_FRONTIER_TASK_WAIT_SECONDS,
            )
        second_owner_status =
            Base.timedwait(
                () -> Base.isready(second_owner_entered),
                INDEL_FRONTIER_TASK_WAIT_SECONDS,
            )
        unrelated_owner_weight =
            Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge)
    finally
        Base.put!(release_first_owner, nothing)
        Base.put!(release_second_owner, nothing)
    end
    Test.@test first_owner_status == :ok
    Test.@test second_owner_status == :ok
    Test.@test Base.fetch(first_owner) == 0.25
    Test.@test Base.fetch(second_owner) == 0.5
    Test.@test unrelated_owner_weight == raw_edge_weight

    # Julia 1.10 does not inherit task-local storage in Threads.@threads tasks.
    # Exercise the production propagation contract explicitly, including the
    # one-thread case where @threads still creates a distinct task.
    propagated_weights = fill(NaN, max(4, Base.Threads.nthreads()))
    unpropagated_weight = Mycelia.Rhizomorph._with_soft_edge_weight_scope(
        branch_graph, concurrent_soft_weights) do
        snapshot = Mycelia.Rhizomorph._current_soft_edge_weight_snapshot()
        unpropagated_task = Base.Threads.@spawn begin
            Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge)
        end
        Base.Threads.@threads for i in eachindex(propagated_weights)
            propagated_weights[i] =
                Mycelia.Rhizomorph._with_soft_edge_weight_snapshot(snapshot) do
                    Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge)
                end
        end
        return Base.fetch(unpropagated_task)
    end
    Test.@test unpropagated_weight == raw_edge_weight
    Test.@test all(==(0.25), propagated_weights)
    Test.@test Mycelia.Rhizomorph.compute_edge_weight(concurrent_edge) ==
               raw_edge_weight
end

Test.@testset "Production parallel decode inherits prior soft weights" begin
    major_sequence = "AAACCCAAA"
    minor_sequence = "AAAGGGAAA"
    query_sequence = "AAATTTAAA"
    graph_reads = FASTX.FASTQ.Record[
        FASTX.FASTQ.Record("major_1", major_sequence, repeat("I", 9)),
        FASTX.FASTQ.Record("major_2", major_sequence, repeat("I", 9)),
        FASTX.FASTQ.Record("minor", minor_sequence, repeat("I", 9)),
    ]
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        graph_reads, 3; mode = :singlestrand)
    queries = FASTX.FASTQ.Record[
        FASTX.FASTQ.Record("query_$(index)", query_sequence, repeat("!", 9))
        for index in 1:4
    ]

    prior_soft_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
    for edge_label in MetaGraphsNext.edge_labels(graph)
        edge_text = string(edge_label[1], edge_label[2])
        weight = occursin('C', edge_text) ?
                 0.001 : occursin('G', edge_text) ? 100.0 : 1.0
        Mycelia.Rhizomorph.accumulate_path_probability!(
            prior_soft_weights, [edge_label], weight)
    end
    no_hoist_builder = _ -> nothing

    raw_result = Mycelia.improve_read_set_likelihood(
        queries,
        graph,
        3;
        graph_mode = :singlestrand,
        beam_width = typemax(Int),
        weighted_graph_builder = no_hoist_builder,
    )
    soft_result = Mycelia.improve_read_set_likelihood(
        queries,
        graph,
        3;
        graph_mode = :singlestrand,
        beam_width = typemax(Int),
        prior_soft_weights = prior_soft_weights,
        weighted_graph_builder = no_hoist_builder,
    )
    parallel_telemetry = Dict{Symbol, Any}()
    parallel_result = Mycelia.improve_read_set_likelihood(
        queries,
        graph,
        3;
        enable_parallel = true,
        graph_mode = :singlestrand,
        beam_width = typemax(Int),
        prior_soft_weights = prior_soft_weights,
        rung_telemetry = parallel_telemetry,
        weighted_graph_builder = no_hoist_builder,
    )

    raw_sequences = [
        FASTX.sequence(String, read) for read in raw_result[1]
    ]
    soft_sequences = [
        FASTX.sequence(String, read) for read in soft_result[1]
    ]
    parallel_sequences = [
        FASTX.sequence(String, read) for read in parallel_result[1]
    ]
    Test.@test raw_sequences == fill(major_sequence, length(queries))
    Test.@test soft_sequences == fill(minor_sequence, length(queries))
    Test.@test parallel_sequences == soft_sequences
    Test.@test parallel_sequences != raw_sequences
    Test.@test parallel_result[2] == length(queries)
    Test.@test parallel_result[3] == 0.0
    Test.@test !parallel_result[5]
    Test.@test parallel_telemetry[:requested] == 0
    Test.@test parallel_telemetry[:attempted] == 0
    Test.@test parallel_telemetry[:completed] == 0
    Test.@test parallel_telemetry[:truncated] == 0
    Test.@test parallel_telemetry[:engaged] == 0

    edge_label = first(MetaGraphsNext.edge_labels(graph))
    edge = graph[edge_label...]
    raw_weight = Mycelia.Rhizomorph.count_total_observations(edge)
    Test.@test Mycelia.Rhizomorph.compute_edge_weight(edge) == raw_weight
    Test.@test Base.fetch(Base.Threads.@spawn(
        Mycelia.Rhizomorph.compute_edge_weight(edge))) == raw_weight
end


Test.@testset "Unrestricted indel schedule preserves exhaustive semantics" begin
    sequence = first(indel_classifier_de_bruijn_sequence(7), 80)
    read = indel_classifier_fastq("unrestricted", sequence)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[read], 7; mode = :singlestrand)
    diagnostics = Mycelia.CorrectorDiagnostics()
    telemetry = Dict{Symbol, Any}()
    indel_frontier_test_throws_message(
        () -> Mycelia.improve_read_set_likelihood(
            FASTX.FASTQ.Record[read],
            graph,
            7;
            indel_schedule = :bogus,
        ),
        ArgumentError,
        "indel_schedule must be :unrestricted or :frontier_budgeted",
    )

    Mycelia.improve_read_set_likelihood(
        FASTX.FASTQ.Record[read],
        graph,
        7;
        graph_mode = :singlestrand,
        diagnostics = diagnostics,
        indel_params = indel_classifier_params(),
        indel_schedule = :unrestricted,
        rung_telemetry = telemetry,
    )

    Test.@test telemetry[:profile_requested]
    Test.@test telemetry[:requested] == 1
    Test.@test telemetry[:attempted] == 1
    Test.@test telemetry[:completed] == 1
    Test.@test telemetry[:truncated] == 0
    Test.@test telemetry[:admitted]
    Test.@test telemetry[:decision_reason] == :unrestricted_semantics
    Test.@test isempty(telemetry[:raw_frontier_metrics])
    Test.@test isempty(telemetry[:cleaned_frontier_metrics])
    Test.@test isnothing(Mycelia._validate_indel_telemetry_counters(telemetry))

    # A legal unrestricted decode still windows genuinely long reads when the
    # caller opts into hard-window decoding. Count the actual pair-HMM dispatch
    # units rather than treating the whole read as one telemetry request.
    long_length = Mycelia._AUTO_BEAM_EXACT_THRESHOLD + 107
    long_sequence = first(indel_classifier_de_bruijn_sequence(7), long_length)
    long_read = indel_classifier_fastq("unrestricted_windowed", long_sequence)
    long_graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[long_read], 7; mode = :singlestrand)
    long_typed_sequence = Mycelia.extract_typed_sequence(
        long_read, BioSequences.LongDNA{4})
    long_kmers = collect(
        Mycelia._record_kmer_iterator(
            BioSequences.LongDNA{4}, 7, long_typed_sequence))
    hard_positions = [100, div(length(long_kmers), 2), length(long_kmers) - 100]
    long_hard = Set(long_kmers[hard_positions])
    long_windows = Mycelia._hard_window_ranges(
        long_read,
        7,
        long_hard;
        pad = 7,
        max_window = Mycelia._INDEL_FRONTIER_MAX_WINDOW,
        graph_mode = :singlestrand,
        complete_span = true,
    )
    Test.@test Mycelia._windowed_decode_read_is_long(long_read, 7)
    Test.@test length(long_windows) == length(hard_positions) > 1

    long_telemetry = Dict{Symbol, Any}()
    Mycelia.improve_read_set_likelihood(
        FASTX.FASTQ.Record[long_read],
        long_graph,
        7;
        graph_mode = :singlestrand,
        hard_vertices = long_hard,
        windowed_decode = true,
        indel_params = indel_classifier_params(),
        indel_schedule = :unrestricted,
        rung_telemetry = long_telemetry,
    )

    expected_units = count(window -> length(window) >= 7, long_windows)
    Test.@test long_telemetry[:requested] == expected_units
    Test.@test long_telemetry[:attempted] == expected_units
    Test.@test long_telemetry[:completed] <= long_telemetry[:attempted]
    Test.@test long_telemetry[:truncated] ==
               long_telemetry[:attempted] - long_telemetry[:completed]
    Test.@test long_telemetry[:engaged] <= long_telemetry[:completed]
    Test.@test all(
        typeof(long_telemetry[key]) === Int && long_telemetry[key] >= 0
        for key in (:requested, :attempted, :completed, :truncated, :engaged)
    )
    Test.@test isnothing(
        Mycelia._validate_indel_telemetry_counters(long_telemetry))

    inexact_telemetry = copy(long_telemetry)
    inexact_telemetry[:requested] = Float64(inexact_telemetry[:requested])
    indel_frontier_test_throws_message(
        () -> Mycelia._validate_indel_telemetry_counters(inexact_telemetry),
        ArgumentError,
        "indel telemetry requested must be an exact Int",
    )
end


Test.@testset "Mixed sparse-window sources preserve substitution bytes" begin
    k = 7
    sequence = first(indel_classifier_de_bruijn_sequence(k), 120)
    read = indel_classifier_fastq("mixed_sources", sequence)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[read], k; mode = :singlestrand)
    typed_sequence = Mycelia.extract_typed_sequence(
        read, BioSequences.LongDNA{4})
    kmers = collect(
        Mycelia._record_kmer_iterator(BioSequences.LongDNA{4}, k, typed_sequence))
    hard = Set(kmers[[10, 50, 90]])
    windows = Mycelia._hard_window_ranges(
        read, k, hard;
        pad = 0,
        max_window = 500,
        graph_mode = :singlestrand,
    )
    Test.@test length(windows) == 3
    sources = Dict{UnitRange{Int}, Symbol}(
        windows[1] => :raw,
        windows[2] => :cleaned,
        windows[3] => :substitution,
    )
    cleaned_graph = deepcopy(graph)
    diagnostics = Mycelia.CorrectorDiagnostics()
    corrected,
    _improved,
    decoded_windows,
    divergent_windows = Mycelia.improve_read_likelihood_windowed_detail(
        read,
        graph,
        k,
        hard;
        graph_mode = :singlestrand,
        weighted_graph = Mycelia.build_correction_weighted_graph(graph),
        cleaned_graph = cleaned_graph,
        cleaned_weighted_graph =
            Mycelia.build_correction_weighted_graph(cleaned_graph),
        indel_window_sources = sources,
        diagnostics = diagnostics,
        pad = 0,
        indel_params = indel_classifier_params(),
    )

    Test.@test decoded_windows == 3
    Test.@test divergent_windows == 0
    Test.@test diagnostics.indel_attempts[] == 2
    Test.@test diagnostics.indel_decodes[] == 2
    Test.@test FASTX.identifier(corrected) == FASTX.identifier(read)
    Test.@test FASTX.sequence(String, corrected) == FASTX.sequence(String, read)
    Test.@test String(FASTX.quality(corrected)) == String(FASTX.quality(read))

    invalid_sources = copy(sources)
    invalid_sources[windows[1]] = :bogus
    invalid_error = nothing
    try
        Mycelia.improve_read_likelihood_windowed_detail(
            read,
            graph,
            k,
            hard;
            graph_mode = :singlestrand,
            weighted_graph = Mycelia.build_correction_weighted_graph(graph),
            cleaned_graph = cleaned_graph,
            cleaned_weighted_graph =
                Mycelia.build_correction_weighted_graph(cleaned_graph),
            indel_window_sources = invalid_sources,
            pad = 0,
            indel_params = indel_classifier_params(),
        )
    catch error
        invalid_error = error
    end
    Test.@test invalid_error isa ArgumentError
    invalid_message = sprint(showerror, invalid_error)
    Test.@test occursin("bogus", invalid_message)
    Test.@test occursin("substitution", invalid_message)

    missing_cleaned_sources = Dict{UnitRange{Int}, Symbol}(
        windows[1] => :cleaned,
        windows[2] => :substitution,
        windows[3] => :substitution,
    )
    missing_cleaned_error = nothing
    try
        Mycelia.improve_read_likelihood_windowed_detail(
            read,
            graph,
            k,
            hard;
            graph_mode = :singlestrand,
            weighted_graph = Mycelia.build_correction_weighted_graph(graph),
            cleaned_graph = nothing,
            cleaned_weighted_graph = nothing,
            indel_window_sources = missing_cleaned_sources,
            pad = 0,
            indel_params = indel_classifier_params(),
        )
    catch error
        missing_cleaned_error = error
    end
    Test.@test missing_cleaned_error isa ArgumentError
    missing_cleaned_message = lowercase(sprint(showerror, missing_cleaned_error))
    Test.@test occursin("cleaned", missing_cleaned_message)
    Test.@test occursin("graph", missing_cleaned_message)
end


Test.@testset "Frontier cleaning respects projected memory" begin
    k = 7
    sequence = first(indel_classifier_de_bruijn_sequence(k), 80)
    read = indel_classifier_fastq("memory_limited", sequence)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[read], k; mode = :singlestrand)
    hard = Set(MetaGraphsNext.labels(graph))
    decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[read],
        graph,
        k,
        hard,
        Bool[false],
        indel_classifier_params(),
        :singlestrand;
        work_limit = 0,
        memory_limit = 1,
    )

    Test.@test !decision.admitted
    Test.@test decision.decision_reason == :cleaning_memory_limit
    Test.@test decision.cleaned_evaluated == 0
    Test.@test decision.graph_source == :substitution
    Test.@test decision.cleanup["graph_cleanup_status"] ==
               "skipped_memory_limit"
    Test.@test decision.cleanup["projected_graph_variants"] == 4
end


Test.@testset "Raw weighted frontier copy respects measured graph memory" begin
    k = 7
    sequence = first(indel_classifier_de_bruijn_sequence(k), 80)
    read = indel_classifier_fastq("raw_weight_memory_limited", sequence)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[read], k; mode = :singlestrand)
    hard = Set(MetaGraphsNext.labels(graph))
    graph_bytes = Mycelia._indel_graph_memory_bytes(graph)
    decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[read],
        graph,
        k,
        hard,
        Bool[false],
        indel_classifier_params(),
        :singlestrand;
        memory_limit = 2 * graph_bytes - 1,
    )

    Test.@test decision.requested == 1
    Test.@test !decision.admitted
    Test.@test decision.admitted_windows == 0
    Test.@test decision.rejected_windows == 1
    Test.@test decision.decision_reason == :weighted_graph_memory_limit
    Test.@test decision.graph_source == :substitution
    Test.@test decision.raw_weighted_graph === nothing
    Test.@test decision.cleanup["graph_cleanup_status"] ==
               "skipped_weighted_memory_limit"
    Test.@test decision.cleanup["measured_graph_bytes"] == graph_bytes
    Test.@test decision.cleanup["projected_graph_variants"] == 2

    raw_oom_builder = IndelOOMWeightedGraphBuilder()
    raw_oom_decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[read],
        graph,
        k,
        hard,
        Bool[false],
        indel_classifier_params(),
        :singlestrand;
        weighted_graph_builder = raw_oom_builder,
    )
    Test.@test raw_oom_builder.calls == 1
    Test.@test !raw_oom_decision.admitted
    Test.@test raw_oom_decision.decision_reason ==
               :weighted_graph_out_of_memory
    Test.@test raw_oom_decision.raw_weighted_graph === nothing
    Test.@test raw_oom_decision.admitted_windows == 0
    Test.@test raw_oom_decision.rejected_windows == raw_oom_decision.requested
    Test.@test all(
        source == :substitution
        for source in values(raw_oom_decision.window_sources[1])
    )
    Test.@test raw_oom_decision.cleanup["graph_cleanup_status"] ==
               "weighted_out_of_memory"

    telemetry = Dict{Symbol, Any}()
    result = Mycelia.improve_read_set_likelihood(
        FASTX.FASTQ.Record[read],
        graph,
        k;
        graph_mode = :singlestrand,
        hard_vertices = hard,
        windowed_decode = true,
        indel_params = indel_classifier_params(),
        indel_schedule = :frontier_budgeted,
        rung_telemetry = telemetry,
        memory_limit = 1,
    )
    # Memory fallback disables this pass operationally but is not a low-k
    # density/floor gate, so it must not populate decode_gated_rungs upstream.
    Test.@test !result[5]
    Test.@test telemetry[:requested] == 1
    Test.@test telemetry[:attempted] == 0
    Test.@test telemetry[:completed] == 0
    Test.@test telemetry[:substitution_decode_memory_gated]
end

Test.@testset "Measured graph variants admit an affordable finite budget" begin
    k = 7
    sequence = first(indel_classifier_de_bruijn_sequence(k), 80)
    read = indel_classifier_fastq("finite_memory_admitted", sequence)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[read], k; mode = :singlestrand)
    graph_bytes = Mycelia._indel_graph_memory_bytes(graph)

    Base.GC.gc()
    live_bytes = Int(Base.gc_live_bytes())
    Test.@test !Mycelia._indel_graph_variants_fit_memory(
        graph, 2, live_bytes + graph_bytes - 1)

    Base.GC.gc()
    live_bytes = Int(Base.gc_live_bytes())
    finite_budget = live_bytes + graph_bytes + 64_000_000
    Test.@test Mycelia._indel_graph_variants_fit_memory(
        graph, 2, finite_budget)

    hard = Set(MetaGraphsNext.labels(graph))
    decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[read],
        graph,
        k,
        hard,
        Bool[false],
        indel_classifier_params(),
        :singlestrand;
        memory_limit = finite_budget,
    )
    Test.@test decision.admitted
    Test.@test decision.raw_weighted_graph !== nothing
end

Test.@testset "Private cleaning rescues a rejected frontier" begin
    k = 7
    backbone =
        "GCAGAAACTCAGTTTCCCTGGGGGTACACACGTAGGCGATTAGGAGTAATGGAGTAACGAAGCGAGCCCATGAGCATGTCGCTAAATTAC"
    error_tips = [
        "GCAGAAACTCAGTTTCCCTGGGACGCAGAGGTGGAGTG",
        "GCAGAAACTCAGTTTCCCTGGGGGTACACACGTAGGCGATAACGACTGCTACTTTC",
        "GCAGAAACTCAGTTTCCCTGGGGGTACACACGTAGGCGATTAGGAGTAATGGAGTAACACGGTTACTTAAGTGA",
    ]

    records = FASTX.FASTQ.Record[
        indel_classifier_fastq("backbone_$index", backbone)
        for index in 1:8
    ]
    for (tip_index, tip) in enumerate(error_tips)
        for replicate in 1:2
            push!(
                records,
                indel_classifier_fastq(
                    "tip_$(tip_index)_$(replicate)",
                    tip,
                ),
            )
        end
    end

    raw_graph = Mycelia.Rhizomorph.build_qualmer_graph(
        records,
        k;
        mode = :singlestrand,
    )
    read = indel_classifier_fastq("cleaning_rescue", backbone)
    observations = indel_classifier_observations(read, k)
    config = Mycelia._indel_probe_config(
        indel_classifier_params(),
        :singlestrand,
    )

    raw_labels_before = Set(MetaGraphsNext.labels(raw_graph))
    raw_edges_before = Set(MetaGraphsNext.edge_labels(raw_graph))
    raw_summary = Mycelia._indel_frontier_graph_summary(raw_graph)
    raw_full = Mycelia._probe_indel_frontier(
        raw_graph,
        observations,
        :DNA;
        config = config,
        strand_mode = :singlestrand,
        work_limit = typemax(Int),
        graph_summary = raw_summary,
    )

    cleaned_preview = deepcopy(raw_graph)
    cleanup_preview = Mycelia.Rhizomorph.clean_corrector_graph!(
        cleaned_preview;
        k = k,
    )
    cleaned_summary = Mycelia._indel_frontier_graph_summary(cleaned_preview)
    cleaned_full = Mycelia._probe_indel_frontier(
        cleaned_preview,
        observations,
        :DNA;
        config = config,
        strand_mode = :singlestrand,
        work_limit = typemax(Int),
        graph_summary = cleaned_summary,
    )

    Test.@test cleanup_preview["graph_cleanup_tips_removed"] > 0
    Test.@test cleaned_summary.vertex_count < raw_summary.vertex_count
    Test.@test cleaned_summary.branch_vertices < raw_summary.branch_vertices
    Test.@test cleaned_full.reason == :complete
    Test.@test cleaned_full.completed_columns == length(observations)
    Test.@test cleaned_full.frontier_work < raw_full.frontier_work

    # The admission boundary depends only on measured topology/frontier work.
    work_limit = cleaned_full.frontier_work
    hard_vertices = Set(MetaGraphsNext.labels(raw_graph))
    decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[read],
        raw_graph,
        k,
        hard_vertices,
        Bool[false],
        indel_classifier_params(),
        :singlestrand;
        work_limit = work_limit,
    )

    raw_metric = only(decision.raw_metrics)
    cleaned_metric = only(decision.cleaned_metrics)
    Test.@test raw_metric[:reason] == :work_limit
    Test.@test !raw_metric[:admitted]
    Test.@test raw_metric[:completed_columns] < raw_metric[:window_length]
    Test.@test cleaned_metric[:reason] == :complete
    Test.@test cleaned_metric[:admitted]
    Test.@test cleaned_metric[:completed_columns] ==
               cleaned_metric[:window_length]
    Test.@test cleaned_metric[:frontier_work] <= work_limit

    Test.@test decision.requested == 1
    Test.@test decision.admitted_windows == 1
    Test.@test decision.rejected_windows == 0
    Test.@test decision.admitted
    Test.@test decision.graph_source == :cleaned
    Test.@test decision.decision_reason == :cleaned_frontier_affordable
    Test.@test only(values(decision.window_sources[1])) == :cleaned
    Test.@test decision.raw_weighted_graph === nothing
    Test.@test decision.cleaned_graph !== nothing
    Test.@test decision.cleaned_weighted_graph !== nothing

    decoded = Mycelia.correct_observations(
        decision.cleaned_graph,
        [observations];
        config = config,
        weighted_graph = decision.cleaned_weighted_graph,
    )
    decoded_path = only(decoded.paths)
    diagnostics = decoded_path.diagnostics
    Test.@test decoded_path.path !== nothing
    Test.@test diagnostics[:algorithm] == :viterbi_indel_pair_hmm
    Test.@test diagnostics[:completed_columns] == length(observations)
    Test.@test !diagnostics[:truncated]
    Test.@test length(only(decoded.corrected_observations)) ==
               length(observations)

    cleaning_oom = IndelOOMGraphCleaner()
    cleaning_builder = IndelOOMWeightedGraphBuilder()
    cleaning_oom_decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[read],
        raw_graph,
        k,
        hard_vertices,
        Bool[false],
        indel_classifier_params(),
        :singlestrand;
        work_limit = work_limit,
        clean_graph! = cleaning_oom,
        weighted_graph_builder = cleaning_builder,
    )
    Test.@test cleaning_oom.calls == 1
    Test.@test cleaning_builder.calls == 0
    Test.@test !cleaning_oom_decision.admitted
    Test.@test cleaning_oom_decision.decision_reason ==
               :cleaning_out_of_memory
    Test.@test cleaning_oom_decision.cleaned_graph === nothing
    Test.@test cleaning_oom_decision.cleaned_weighted_graph === nothing
    Test.@test only(values(cleaning_oom_decision.window_sources[1])) ==
               :substitution
    Test.@test cleaning_oom_decision.cleanup["graph_cleanup_status"] ==
               "out_of_memory"
    Test.@test cleaning_oom_decision.cleanup["graph_cleanup_error_type"] ==
               "OutOfMemoryError"

    # The same fixture is raw-rejected but cleaning-admitted. A throwing builder
    # therefore reaches the cleaned-weighted allocation, not the raw or pass-level
    # allocation branch.
    cleaned_oom_builder = IndelOOMWeightedGraphBuilder()
    cleaned_oom_decision = Mycelia._evaluate_indel_frontier_schedule(
        FASTX.FASTQ.Record[read],
        raw_graph,
        k,
        hard_vertices,
        Bool[false],
        indel_classifier_params(),
        :singlestrand;
        work_limit = work_limit,
        weighted_graph_builder = cleaned_oom_builder,
    )
    Test.@test cleaned_oom_builder.calls == 1
    Test.@test !cleaned_oom_decision.admitted
    Test.@test cleaned_oom_decision.decision_reason ==
               :cleaning_out_of_memory
    Test.@test cleaned_oom_decision.cleaned_graph === nothing
    Test.@test cleaned_oom_decision.cleaned_weighted_graph === nothing
    Test.@test only(values(cleaned_oom_decision.window_sources[1])) ==
               :substitution
    Test.@test cleaned_oom_decision.cleanup["graph_cleanup_status"] ==
               "out_of_memory"
    Test.@test cleaned_oom_decision.cleanup["graph_cleanup_error_type"] ==
               "OutOfMemoryError"

    # Scheduler cleanup operated only on its private graph.
    Test.@test Set(MetaGraphsNext.labels(raw_graph)) == raw_labels_before
    Test.@test Set(MetaGraphsNext.edge_labels(raw_graph)) == raw_edges_before
end

Test.@testset "clean_corrector_graph! includes support-two tips" begin
    k = 7
    backbone =
        "GCAGAAACTCAGTTTCCCTGGGGGTACACACGTAGGCGATTAGGAGTAATGGAGTAACGAAGCGAGCCCATGAGCATGTCGCTAAATTAC"
    support_two_tip =
        "GCAGAAACTCAGTTTCCCTGGGACGCAGAGGTGGAGTG"
    support_three_tip =
        "GCAGAAACTCAGTTTCCCTGGGGGTACACACGTAGGCGATTAGGAGTAATGGAGTAACACGGTTACTTAAGTGA"

    records = FASTX.FASTQ.Record[
        indel_classifier_fastq("backbone_$index", backbone)
        for index in 1:8
    ]
    append!(
        records,
        [
            indel_classifier_fastq("support_two_$index", support_two_tip)
            for index in 1:2
        ],
    )
    append!(
        records,
        [
            indel_classifier_fastq(
                "support_three_$index",
                support_three_tip,
            )
            for index in 1:3
        ],
    )

    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        records,
        k;
        mode = :singlestrand,
    )
    support_two_label = Kmers.DNAKmer{k}("AGAGGTG")
    support_three_label = Kmers.DNAKmer{k}("TTAAGTG")
    Test.@test Mycelia.Rhizomorph.count_evidence(
        graph[support_two_label],
    ) == 2
    Test.@test Mycelia.Rhizomorph.count_evidence(
        graph[support_three_label],
    ) == 3

    below_boundary = deepcopy(graph)
    below_stats = Mycelia.Rhizomorph.clean_corrector_graph!(
        below_boundary;
        k = k,
        max_tip_support = 1,
        collapse_bubbles = false,
        prune_components = false,
    )
    Test.@test below_stats["graph_cleanup_tips_removed"] == 0
    Test.@test haskey(below_boundary, support_two_label)
    Test.@test haskey(below_boundary, support_three_label)

    at_boundary = deepcopy(graph)
    boundary_stats = Mycelia.Rhizomorph.clean_corrector_graph!(
        at_boundary;
        k = k,
        max_tip_support = 2,
        collapse_bubbles = false,
        prune_components = false,
    )
    Test.@test boundary_stats["graph_cleanup_tips_removed"] > 0
    Test.@test !haskey(at_boundary, support_two_label)
    Test.@test haskey(at_boundary, support_three_label)
end
