# td-2rxh Step 5 toy proof: staged indel correction on the same long-read fixture.
#
# Run from the Mycelia repository root:
#   LD_LIBRARY_PATH='' julia --project=. scratchpad/step5_toy_proof.jl
#
# The fixture and seeds are fixed before observing any accuracy result. Both arms
# receive byte-identical indel-rich nanopore reads; only the correction profile
# differs. The vertex ceiling remains the package constant derived from the
# independent decode-time calibration and is never tuned here.

import BioSequences
import FASTX
import Mycelia
import Random

const GENOME_LENGTH = 2_000
const SOURCE_READ_LENGTH = 1_200
const COVERAGE = 8
const ERROR_RATE = 0.05
const FIXTURE_SEED = 42
const CORRECTOR_SEED = 1_042
const MAX_K = 31
const MAX_NANOPORE_WALL_SECONDS = 120.0
const DECODE_BUDGET_MS = 200.0
const DECODE_VERTEX_CALIBRATION = [
    (n_vertices = 1_996, confirmation_medians_ms = (165.643, 153.788)),
    (n_vertices = 3_610, confirmation_medians_ms = (263.639, 279.200)),
    (n_vertices = 4_972, confirmation_medians_ms = (475.403, 326.281)),
    (n_vertices = 6_068, confirmation_medians_ms = (451.029, 383.677)),
]

function make_fixture()::Tuple{
        Vector{FASTX.FASTQ.Record}, BioSequences.LongDNA{4}}
    reference_record = Mycelia.random_fasta_record(
        moltype = :DNA,
        seed = FIXTURE_SEED,
        L = GENOME_LENGTH,
    )
    reference = FASTX.sequence(BioSequences.LongDNA{4}, reference_record)
    rng = Random.MersenneTwister(FIXTURE_SEED)
    Random.seed!(FIXTURE_SEED)

    n_reads = ceil(Int, COVERAGE * GENOME_LENGTH / SOURCE_READ_LENGTH)
    reads = FASTX.FASTQ.Record[]
    for read_index in 1:n_reads
        start_position = rand(
            rng,
            1:(GENOME_LENGTH - SOURCE_READ_LENGTH + 1),
        )
        fragment = reference[
            start_position:(start_position + SOURCE_READ_LENGTH - 1)]
        if rand(rng, Bool)
            fragment = BioSequences.reverse_complement(fragment)
        end
        observed, qualities = Mycelia.observe(
            fragment;
            error_rate = ERROR_RATE,
            tech = :nanopore,
        )
        isempty(observed) && continue
        quality_string = String([Char(quality + 33) for quality in qualities])
        push!(
            reads,
            FASTX.FASTQ.Record(
                "nanopore_read_$(read_index)",
                string(observed),
                quality_string,
            ),
        )
    end
    return reads, reference
end

function calculate_n50(contigs::Vector{String})::Int
    isempty(contigs) && return 0
    lengths = sort(length.(contigs); rev = true)
    target = cld(sum(lengths), 2)
    cumulative = 0
    for contig_length in lengths
        cumulative += contig_length
        cumulative >= target && return contig_length
    end
    return 0
end

function best_reference_alignment(
        contigs::Vector{String},
        reference::BioSequences.LongDNA{4},
)::NamedTuple
    reference_forward = string(reference)
    reference_reverse = string(BioSequences.reverse_complement(reference))
    best = (
        identity = 0.0,
        edit_distance = typemax(Int),
        matches = 0,
        aligned_bases = 0,
        contig_length = 0,
        orientation = :none,
    )

    for contig in contigs
        for (orientation, target) in (
                (:forward, reference_forward),
                (:reverse, reference_reverse),
        )
            alignment = Mycelia.assess_alignment(contig, target)
            aligned_bases = alignment.total_matches + alignment.total_edits
            identity = aligned_bases == 0 ?
                       0.0 : alignment.total_matches / aligned_bases
            candidate = (
                identity = identity,
                edit_distance = alignment.total_edits,
                matches = alignment.total_matches,
                aligned_bases = aligned_bases,
                contig_length = length(contig),
                orientation = orientation,
            )
            if candidate.identity > best.identity ||
               (candidate.identity == best.identity &&
                candidate.matches > best.matches) ||
               (candidate.identity == best.identity &&
                candidate.matches == best.matches &&
                candidate.edit_distance < best.edit_distance)
                best = candidate
            end
        end
    end
    return best
end

function run_arm(
        reads::Vector{FASTX.FASTQ.Record},
        reference::BioSequences.LongDNA{4},
        sequencing_tech::Symbol,
)::NamedTuple
    # Reset the corrector RNG so repeated script runs are deterministic and both
    # correction profiles begin from the same random state.
    Random.seed!(CORRECTOR_SEED)
    local assembly::Mycelia.Rhizomorph.AssemblyResult
    wall_seconds = @elapsed begin
        assembly = Mycelia.Rhizomorph.assemble_genome(
            deepcopy(reads);
            k = MAX_K,
            corrector = :iterative,
            strategy = :scalable,
            sequencing_tech = sequencing_tech,
        )
    end

    alignment = best_reference_alignment(assembly.contigs, reference)
    stats = assembly.assembly_stats
    return (
        sequencing_tech = sequencing_tech,
        wall_seconds = wall_seconds,
        identity = alignment.identity,
        edit_distance = alignment.edit_distance,
        matches = alignment.matches,
        aligned_bases = alignment.aligned_bases,
        best_contig_length = alignment.contig_length,
        best_orientation = alignment.orientation,
        n_contigs = length(assembly.contigs),
        total_assembled_bases = sum(length, assembly.contigs; init = 0),
        largest_contig = isempty(assembly.contigs) ? 0 : maximum(length, assembly.contigs),
        n50 = calculate_n50(assembly.contigs),
        k_progression = get(stats, "k_progression", Int[]),
        rung_vertex_counts = get(
            stats,
            "rung_vertex_counts",
            Dict{Int, Vector{Int}}(),
        ),
        staged_indel_rungs = get(stats, "staged_indel_rungs", NamedTuple[]),
        indel_moves = get(stats, "indel_moves", false),
        indel_decodes = get(stats, "indel_decodes", 0),
        truncated_decodes = get(stats, "truncated_decodes", 0),
        trace_contract_errors = get(stats, "trace_contract_errors", 0),
        window_divergences = get(stats, "window_divergences", 0),
        reassembly_k = get(stats, "reassembly_k", nothing),
        reassembly_k_ceiling = get(stats, "reassembly_k_ceiling", nothing),
    )
end

function rung_k(rung::Any)::Union{Int, Nothing}
    if rung isa NamedTuple
        return get(rung, :k, nothing)
    elseif rung isa Pair && first(rung) isa Int
        return first(rung)
    end
    return nothing
end

function rung_vertices(rung::Any)::Union{Int, Nothing}
    if rung isa NamedTuple
        return get(rung, :n_vertices, get(rung, :nv, nothing))
    elseif rung isa Pair && last(rung) isa Int
        return last(rung)
    end
    return nothing
end

function print_arm(result::NamedTuple)::Nothing
    println("\n$(result.sequencing_tech) correction arm")
    println("  wall_seconds:          $(round(result.wall_seconds; digits = 3))")
    println("  identity:              $(round(result.identity; digits = 6))")
    println("  edit_distance:         $(result.edit_distance)")
    println("  matches/aligned:       $(result.matches)/$(result.aligned_bases)")
    println("  best_contig_length:    $(result.best_contig_length)")
    println("  best_orientation:      $(result.best_orientation)")
    println("  contigs/total/largest: $(result.n_contigs)/" *
            "$(result.total_assembled_bases)/$(result.largest_contig)")
    println("  N50:                   $(result.n50)")
    println("  k_progression:         $(result.k_progression)")
    println("  graph vertices/pass:   $(result.rung_vertex_counts)")
    println("  reassembly k/ceiling:  $(result.reassembly_k)/" *
            "$(result.reassembly_k_ceiling)")
    println("  indel_moves:           $(result.indel_moves)")
    println("  indel_decodes:         $(result.indel_decodes)")
    println("  truncated_decodes:     $(result.truncated_decodes)")
    println("  trace_contract_errors: $(result.trace_contract_errors)")
    println("  window_divergences:    $(result.window_divergences)")
    println("  staged indel rungs (actual k/nv):")
    if isempty(result.staged_indel_rungs)
        println("    none")
    else
        for rung in result.staged_indel_rungs
            println("    k=$(rung_k(rung)), nv=$(rung_vertices(rung))")
        end
    end
    return nothing
end

function main()::Nothing
    reads, reference = make_fixture()
    observed_lengths = length.(FASTX.sequence.(String, reads))
    println("td-2rxh fixed Step 5 toy fixture")
    println("  fixture_seed:          $(FIXTURE_SEED)")
    println("  corrector_seed:        $(CORRECTOR_SEED)")
    println("  reference/source read: $(GENOME_LENGTH)/$(SOURCE_READ_LENGTH) bp")
    println("  reads/coverage/error:  $(length(reads))/$(COVERAGE)x/$(ERROR_RATE)")
    println("  observed read lengths: min=$(minimum(observed_lengths)), " *
            "max=$(maximum(observed_lengths))")
    println("  max_k:                 $(MAX_K)")
    println("  nanopore budget:       $(MAX_NANOPORE_WALL_SECONDS) seconds")
    println("  pair-HMM budget:       $(DECODE_BUDGET_MS) ms/read")
    println("  decode calibration (runtime only; not accuracy-fitted):")
    for point in DECODE_VERTEX_CALIBRATION
        println("    nv=$(point.n_vertices): independent medians=" *
                "$(point.confirmation_medians_ms) ms")
    end
    println("  selected ceiling:      $(Mycelia._INDEL_DECODE_MAX_VERTICES) vertices")

    # The nanopore arm runs first, so its wall-time check conservatively includes
    # compilation of the indel kernel rather than benefiting from the control arm.
    nanopore = run_arm(reads, reference, :nanopore)
    illumina = run_arm(reads, reference, :illumina)
    print_arm(nanopore)
    print_arm(illumina)

    staged_ks = filter(!isnothing, rung_k.(nanopore.staged_indel_rungs))
    if !isempty(staged_ks) && all(==(minimum(nanopore.k_progression)), staged_ks)
        println("\nWARNING: indel staging fired only at the initial rung. " *
                "Report this as an initial-rung result, not a sparse high-rung proof.")
    end

    checks = (
        nonempty_nanopore = nanopore.n_contigs > 0,
        nonempty_illumina = illumina.n_contigs > 0,
        max_k_reached = MAX_K in nanopore.k_progression,
        staged_indel_rung_observed = !isempty(nanopore.staged_indel_rungs),
        nanopore_indel_decode_ran = nanopore.indel_decodes > 0,
        nanopore_trace_contract_clean = nanopore.trace_contract_errors == 0,
        illumina_oracle_has_no_indel_rung = isempty(illumina.staged_indel_rungs),
        illumina_oracle_has_no_indel_decode = illumina.indel_decodes == 0,
        nanopore_beats_illumina = nanopore.identity > illumina.identity,
        nanopore_within_budget =
            nanopore.wall_seconds < MAX_NANOPORE_WALL_SECONDS,
    )

    println("\nStep 5 checks")
    for (name, passed) in pairs(checks)
        println("  $(name): $(passed ? "PASS" : "FAIL")")
    end
    failed = [String(name) for (name, passed) in pairs(checks) if !passed]
    if !isempty(failed)
        error("td-2rxh Step 5 toy proof failed: $(join(failed, ", "))")
    end
    println("\ntd-2rxh Step 5 toy proof: PASS")
    return nothing
end

main()
