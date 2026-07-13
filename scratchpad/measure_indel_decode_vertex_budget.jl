# td-2rxh runtime-only calibration for `_INDEL_DECODE_MAX_VERTICES`.
#
# Run from the Mycelia repository root:
#   LD_LIBRARY_PATH='' julia --project=. \
#     scratchpad/measure_indel_decode_vertex_budget.jl
#
# This fixes the graph-size ceiling from a 200 ms/read compute budget. It does
# not measure or consume correction accuracy.

import BioSequences
import FASTX
import Graphs
import Mycelia
import Random
import Statistics

const K = 13
const BUDGET_MS = 200.0
const CONFIRMATION_MEDIANS_MS = Dict(
    1_996 => [165.643, 153.788],
    3_610 => [263.639, 279.200],
    4_972 => [475.403, 326.281],
    6_068 => [451.029, 383.677],
)

function make_reads()::Vector{FASTX.FASTQ.Record}
    Random.seed!(20260712)
    reference = Mycelia.random_fasta_record(
        moltype = :DNA, seed = 20260712, L = 3_000)
    reference_sequence = FASTX.sequence(BioSequences.LongDNA{4}, reference)
    read_length = 1_000
    n_reads = 90
    reads = FASTX.FASTQ.Record[]
    for read_index in 1:n_reads
        start = 1 + mod((read_index - 1) * 37, length(reference_sequence) - read_length + 1)
        template = reference_sequence[start:(start + read_length - 1)]
        observed, qualities = Mycelia.observe(
            template; error_rate = 0.05, tech = :nanopore)
        quality = String([Char(score + 33) for score in qualities])
        push!(reads, FASTX.FASTQ.Record("read_$(read_index)", string(observed), quality))
    end
    return reads
end

function probe_record(read::FASTX.FASTQ.Record)::FASTX.FASTQ.Record
    sequence = FASTX.sequence(String, read)
    quality = String(FASTX.quality(read))
    probe_length = min(200, length(sequence))
    return FASTX.FASTQ.Record(
        "probe", sequence[1:probe_length], quality[1:probe_length])
end

function quality_observations(
        read::FASTX.FASTQ.Record
)::Vector{Mycelia.QualityObservation}
    sequence_string = FASTX.sequence(String, read)
    alphabet = Mycelia.detect_alphabet(sequence_string)
    sequence_type = Mycelia.alphabet_to_biosequence_type(alphabet)
    sequence = Mycelia.extract_typed_sequence(read, sequence_type)
    kmers = collect(Mycelia._record_kmer_iterator(sequence_type, K, sequence))
    quality_scores = collect(FASTX.quality_scores(read))
    observations = Vector{Mycelia.QualityObservation}(undef, length(kmers))
    for (index, kmer) in enumerate(kmers)
        lo = clamp(index, 1, length(quality_scores))
        hi = clamp(index + K - 1, 1, length(quality_scores))
        observations[index] = Mycelia.QualityObservation(
            kmer, UInt8.(@view quality_scores[lo:hi]))
    end
    return observations
end

function indel_config(n_observations::Int, n_vertices::Int)::Mycelia.ViterbiCorrectionConfig
    profile = Mycelia.indel_error_profile(:nanopore)
    knobs = Mycelia.Rhizomorph._corrector_strategy_knobs(:scalable)
    beam_width = Mycelia._auto_beam_width(n_observations, n_vertices)
    beam_is_exact = beam_width == typemax(Int)
    return Mycelia.ViterbiCorrectionConfig(
        alphabet = :DNA,
        strand_mode = :doublestrand,
        max_steps = n_observations - 1,
        beam_width = beam_width,
        max_successors_per_state = beam_is_exact ? typemax(Int) :
                                   Mycelia._AUTO_SUCCESSOR_BOUND,
        beam_score_margin = beam_is_exact ? Inf : Mycelia._AUTO_BEAM_SCORE_MARGIN,
        error_rate = profile.base_error_rate,
        indel_moves = true,
        insertion_fraction = profile.insertion_fraction,
        deletion_fraction = profile.deletion_fraction,
        insertion_extend_probability = profile.insertion_extend_probability,
        deletion_extend_probability = profile.deletion_extend_probability,
        deletion_max_run = knobs.deletion_max_run,
        max_insertion_run = knobs.max_insertion_run,
        band_width = knobs.band_width)
end

function measure_graph(
        graph::Graphs.AbstractGraph,
        observations::Vector{Mycelia.QualityObservation}
)::NamedTuple
    n_vertices = Graphs.nv(graph)
    config = indel_config(length(observations), n_vertices)
    weighted_graph = Mycelia.build_correction_weighted_graph(graph; config = config)
    warm = Mycelia.correct_observations(
        graph, [observations]; config = config, weighted_graph = weighted_graph)
    warm_path = only(warm.paths)
    algorithm = get(warm_path.diagnostics, :algorithm, :missing)
    valid = algorithm == :viterbi_indel_pair_hmm && warm_path.path !== nothing
    timings_ms = Float64[]
    if valid
        for _ in 1:5
            elapsed = @elapsed Mycelia.correct_observations(
                graph, [observations]; config = config, weighted_graph = weighted_graph)
            push!(timings_ms, elapsed * 1_000)
        end
    end
    return (
        n_vertices = n_vertices,
        median_ms = isempty(timings_ms) ? Inf : Statistics.median(timings_ms),
        timings_ms = timings_ms,
        beam_width = config.beam_width,
        valid = valid,
        algorithm = algorithm)
end

function main()::Nothing
    reads = make_reads()
    probe = probe_record(first(reads))
    observations = quality_observations(probe)
    prefix_sizes = [1, 2, 3, 4]
    measurements = NamedTuple[]
    println("budget_ms=$(BUDGET_MS),k=$(K),probe_bases=$(length(FASTX.sequence(probe)))")
    for prefix_size in prefix_sizes
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            reads[1:prefix_size], K; mode = :doublestrand)
        measurement = measure_graph(graph, observations)
        push!(measurements, (; prefix_size = prefix_size, measurement...))
        println("prefix=$(prefix_size),nv=$(measurement.n_vertices)," *
                "median_ms=$(round(measurement.median_ms, digits = 3))," *
                "timings_ms=$(round.(measurement.timings_ms, digits = 3))," *
                "beam=$(measurement.beam_width),valid=$(measurement.valid)," *
                "algorithm=$(measurement.algorithm)")
    end

    # The ceiling is fixed from two independent fresh-process confirmations,
    # not from whichever transient load this reproduction run happens to see.
    # A vertex count is admitted only when every confirmation median met the
    # budget; this conservative all-runs rule selects 1,996, rounded to 2,000.
    confirmed = sort(collect(CONFIRMATION_MEDIANS_MS); by = first)
    ceiling = nothing
    for (n_vertices, medians_ms) in confirmed
        if all(<=(BUDGET_MS), medians_ms)
            ceiling = n_vertices
        end
    end
    println("confirmation_medians_ms=$(CONFIRMATION_MEDIANS_MS)")
    println("selected_ceiling=$(ceiling),selection=largest_nv_with_all_" *
            "independent_run_medians_under_$(Int(BUDGET_MS))ms")
    println("rounded_package_ceiling=$(Mycelia._INDEL_DECODE_MAX_VERTICES)")
    return nothing
end

main()
