# Score-free topology-frontier probe and pair-HMM runtime telemetry.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/indel_frontier_probe_test.jl")'

import BioSequences
import FASTX
import Graphs
import MetaGraphsNext
import Mycelia
import Test

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

function indel_frontier_probe_config()::Mycelia.ViterbiCorrectionConfig
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
        band_width = 4,
        beam_width = typemax(Int)
    )
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
        raw_registration_count = length(
            Mycelia.Rhizomorph._SOFT_EDGE_WEIGHT_REGISTRY)
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
        Test.@test length(Mycelia.Rhizomorph._SOFT_EDGE_WEIGHT_REGISTRY) >
                   raw_registration_count
        Test.@test count(
            value -> isapprox(value, 0.25),
            values(Mycelia.Rhizomorph._SOFT_EDGE_WEIGHT_REGISTRY),
        ) > raw_registration_count
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
    Test.@test linear_decision.requested == 1
    Test.@test linear_decision.admitted
    Test.@test linear_decision.admitted_windows == 1
    Test.@test linear_decision.rejected_windows == 0
    Test.@test only(values(linear_decision.window_sources[1])) == :raw
end


Test.@testset "Unrestricted indel schedule preserves exhaustive semantics" begin
    sequence = first(indel_classifier_de_bruijn_sequence(7), 80)
    read = indel_classifier_fastq("unrestricted", sequence)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        FASTX.FASTQ.Record[read], 7; mode = :singlestrand)
    diagnostics = Mycelia.CorrectorDiagnostics()
    telemetry = Dict{Symbol, Any}()
    Test.@test_throws ArgumentError Mycelia.improve_read_set_likelihood(
        FASTX.FASTQ.Record[read],
        graph,
        7;
        indel_schedule = :bogus,
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
end
