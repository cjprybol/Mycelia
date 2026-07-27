# Diagnostic for the 2 kb / 1.2 kb-read / 8x / 5%-error correction REGRESSION
# (td-jt7r): at that cell `corrector = :iterative` produces a WORSE assembly than
# `corrector = :none` (largest contig 413 -> 131 bp), while at larger genomes /
# higher coverage correction helps.
#
# The script is a measurement harness, not a fix. For one (genome, read length,
# coverage, error rate) cell it reports, for BOTH arms:
#
#   * assembly quality (n_contigs, largest contig, largest/genome ratio)
#   * the coverage-aware re-assembly k actually chosen (`reassembly_k`) and the
#     ceiling it was chosen from, plus the `median_solid_kmer_multiplicity`
#     ladder that drives the choice, on RAW and CORRECTED reads
#   * DIRECT over-correction evidence: per-read identity of RAW reads and of
#     CORRECTED reads against the known reference (semi-global BioAlignments
#     pairwise alignment, best orientation), plus the fraction of each read's
#     31-mers that occur in the reference.
#
# Corrected reads are obtained by giving `assemble_genome` an `output_dir`, which
# persists the Stage-1 corrector handoff at `<output_dir>/corrected.fastq`
# instead of deleting it (see `_run_stage1_correction`).
#
# Usage:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/indel_toy_overcorrection_diagnostic.jl \
#     --genome=2000 --read-length=1200 --coverage=8 --error=0.05 \
#     --seed=42 --out=/tmp/diag-cell --label=baseline
#
# Emits one CSV row per cell to `<out>/cell_summary.csv` and a k-ladder CSV to
# `<out>/k_ladder.csv`.

import BioAlignments
import BioSequences
import CSV
import DataFrames
import FASTX
import Mycelia
import Random
import Statistics

const DIAG_MAX_K = 31
const DIAG_IDENTITY_KMER = 31

"""
Parse `--key=value` CLI arguments into a `Dict{String, String}`.
"""
function parse_args(args::Vector{String})::Dict{String, String}
    parsed = Dict{String, String}()
    for arg in args
        startswith(arg, "--") || continue
        body = arg[3:end]
        pieces = split(body, '='; limit = 2)
        parsed[String(pieces[1])] = length(pieces) == 2 ? String(pieces[2]) : "true"
    end
    return parsed
end

"""
Build the read fixture exactly as `benchmarking/indel_frontier_fixed_toy_proof.jl`
does, but with genome length / read length / coverage / error rate as parameters.
Reads are simulated through `Mycelia.observe(...; tech = :nanopore)` so they carry
indels as well as substitutions.
"""
function make_fixture(
        genome_length::Int,
        read_length::Int,
        coverage::Int,
        error_rate::Float64,
        seed::Int
)::Tuple{Vector{FASTX.FASTQ.Record}, BioSequences.LongDNA{4}}
    reference_record = Mycelia.random_fasta_record(
        moltype = :DNA, seed = seed, L = genome_length
    )
    reference = FASTX.sequence(BioSequences.LongDNA{4}, reference_record)
    rng = Random.MersenneTwister(seed)
    Random.seed!(seed)
    n_reads = ceil(Int, coverage * genome_length / read_length)
    reads = FASTX.FASTQ.Record[]
    for read_index in 1:n_reads
        start_position = rand(rng, 1:(genome_length - read_length + 1))
        fragment = reference[start_position:(start_position + read_length - 1)]
        if rand(rng, Bool)
            fragment = BioSequences.reverse_complement(fragment)
        end
        observed,
        qualities = Mycelia.observe(
            fragment; error_rate = error_rate, tech = :nanopore
        )
        isempty(observed) && continue
        quality_string = String([Char(quality + 33) for quality in qualities])
        push!(
            reads,
            FASTX.FASTQ.Record(
                "read_$(read_index)", string(observed), quality_string
            )
        )
    end
    return reads, reference
end

"""
Local (Smith-Waterman) identity of `query` against `reference`, taking the better
of the two orientations. Identity is `matches / (matches + mismatches +
insertions + deletions)` INSIDE the aligned block, so indels are penalized while
the unaligned reference flanks are not. `SemiGlobalAlignment` is deliberately NOT
used: BioAlignments still emits the free reference end-gaps as deletion
operations, so a 1.2 kb read against a 2 kb reference scores ~0.57 identity purely
from the 800 bp of untouched flank. Returns the identity and the fraction of the
query consumed by the local block (so a truncated alignment cannot masquerade as
high identity).
"""
function reference_identity(
        query::BioSequences.LongDNA{4},
        reference::BioSequences.LongDNA{4},
        score_model
)::Tuple{Float64, Float64}
    best_identity = 0.0
    best_coverage = 0.0
    query_length = length(query)
    for candidate in (query, BioSequences.reverse_complement(query))
        result = BioAlignments.pairalign(
            BioAlignments.LocalAlignment(), candidate, reference, score_model
        )
        alignment = BioAlignments.alignment(result)
        matches = BioAlignments.count_matches(alignment)
        mismatches = BioAlignments.count_mismatches(alignment)
        insertions = BioAlignments.count_insertions(alignment)
        deletions = BioAlignments.count_deletions(alignment)
        denominator = matches + mismatches + insertions + deletions
        denominator == 0 && continue
        identity = matches / denominator
        if identity > best_identity
            best_identity = identity
            best_coverage = query_length == 0 ? 0.0 :
                            (matches + mismatches + insertions) / query_length
        end
    end
    return (best_identity, best_coverage)
end

"""
Canonical k-mer hash set of `sequence` at size `k`, as `String` windows over the
canonical orientation. Used for the reference-k-mer-recovery statistic.
"""
function canonical_kmer_set(sequence::BioSequences.LongDNA{4}, k::Int)::Set{String}
    kmers = Set{String}()
    length(sequence) < k && return kmers
    for start_index in 1:(length(sequence) - k + 1)
        window = sequence[start_index:(start_index + k - 1)]
        forward = string(window)
        reverse = string(BioSequences.reverse_complement(window))
        push!(kmers, min(forward, reverse))
    end
    return kmers
end

"""
Fraction of `sequence`'s canonical k-mers that occur in `reference_kmers`. A read
that has been corrected TOWARD the reference raises this; a read corrected toward
a wrong consensus lowers it.
"""
function reference_kmer_fraction(
        sequence::BioSequences.LongDNA{4},
        reference_kmers::Set{String},
        k::Int
)::Union{Float64, Nothing}
    length(sequence) < k && return nothing
    total = 0
    hits = 0
    for start_index in 1:(length(sequence) - k + 1)
        window = sequence[start_index:(start_index + k - 1)]
        forward = string(window)
        reverse = string(BioSequences.reverse_complement(window))
        total += 1
        hits += (min(forward, reverse) in reference_kmers) ? 1 : 0
    end
    total == 0 && return nothing
    return hits / total
end

"""
Per-read identity + reference-k-mer-recovery statistics for one read set.
"""
function read_set_stats(
        sequences::Vector{BioSequences.LongDNA{4}},
        reference::BioSequences.LongDNA{4},
        reference_kmers::Set{String},
        score_model
)::NamedTuple
    identities = Float64[]
    coverages = Float64[]
    kmer_fractions = Float64[]
    for sequence in sequences
        identity, coverage = reference_identity(sequence, reference, score_model)
        push!(identities, identity)
        push!(coverages, coverage)
        fraction = reference_kmer_fraction(sequence, reference_kmers, DIAG_IDENTITY_KMER)
        fraction === nothing || push!(kmer_fractions, fraction)
    end
    return (;
        n = length(sequences),
        mean_identity = isempty(identities) ? NaN : Statistics.mean(identities),
        median_identity = isempty(identities) ? NaN : Statistics.median(identities),
        min_identity = isempty(identities) ? NaN : minimum(identities),
        mean_aligned_fraction = isempty(coverages) ? NaN : Statistics.mean(coverages),
        mean_length = isempty(sequences) ? NaN :
                      Statistics.mean(Float64.(length.(sequences))),
        mean_reference_kmer_fraction = isempty(kmer_fractions) ? NaN :
                                       Statistics.mean(kmer_fractions),
        identities = identities
    )
end

"""
Read a FASTQ file into `LongDNA{4}` sequences.
"""
function read_fastq_sequences(path::AbstractString)::Vector{BioSequences.LongDNA{4}}
    sequences = BioSequences.LongDNA{4}[]
    open(FASTX.FASTQ.Reader, String(path)) do reader
        for record in reader
            push!(sequences, FASTX.sequence(BioSequences.LongDNA{4}, record))
        end
    end
    return sequences
end

"""
`median_solid_kmer_multiplicity` at every prime candidate k up to `ceiling_k`,
plus the k `select_reassembly_k` would choose. This is the signal that decides the
corrected-read re-assembly k.
"""
function k_ladder(reads, ceiling_k::Int)::NamedTuple
    candidate_ks = [7, 11, 13, 17, 19, 23, 29, ceiling_k]
    unique!(sort!(candidate_ks))
    filter!(<=(ceiling_k), candidate_ks)
    medians = [Mycelia.Rhizomorph.median_solid_kmer_multiplicity(reads, k)
               for k in candidate_ks]
    chosen = Mycelia.Rhizomorph.select_reassembly_k(reads, ceiling_k)
    return (; candidate_ks, medians, chosen)
end

"""
Assemble one arm and return quality + provenance fields.
"""
function run_arm(
        reads::Vector{FASTX.FASTQ.Record},
        genome_length::Int,
        corrector::Symbol,
        output_dir::Union{String, Nothing}
)::NamedTuple
    Random.seed!(1_042)
    assembly = nothing
    wall_seconds = @elapsed begin
        if corrector == :none
            assembly = Mycelia.Rhizomorph.assemble_genome(
                deepcopy(reads);
                k = DIAG_MAX_K,
                corrector = :none,
                strategy = :scalable,
                sequencing_tech = :illumina
            )
        else
            assembly = Mycelia.Rhizomorph.assemble_genome(
                deepcopy(reads);
                k = DIAG_MAX_K,
                corrector = :iterative,
                strategy = :scalable,
                sequencing_tech = :illumina,
                output_dir = output_dir
            )
        end
    end
    contig_lengths = [length(contig) for contig in assembly.contigs]
    largest = isempty(contig_lengths) ? 0 : maximum(contig_lengths)
    stats = assembly.assembly_stats
    return (;
        corrector = String(corrector),
        n_contigs = length(assembly.contigs),
        largest = largest,
        ratio = largest / genome_length,
        wall_seconds = wall_seconds,
        reassembly_k = get(stats, "reassembly_k", missing),
        reassembly_k_ceiling = get(stats, "reassembly_k_ceiling", missing),
        reassembly_graph_reused = get(stats, "reassembly_graph_reused", missing),
        corrected_read_count = get(stats, "corrected_read_count", missing)
    )
end

function main(args::Vector{String} = ARGS)::Nothing
    parsed = parse_args(args)
    genome_length = parse(Int, get(parsed, "genome", "2000"))
    read_length = parse(Int, get(parsed, "read-length", "1200"))
    coverage = parse(Int, get(parsed, "coverage", "8"))
    error_rate = parse(Float64, get(parsed, "error", "0.05"))
    seed = parse(Int, get(parsed, "seed", "42"))
    label = get(parsed, "label", "cell")
    output_dir = abspath(get(parsed, "out", joinpath(tempdir(), "indel-toy-diag")))
    skip_alignment = get(parsed, "skip-alignment", "false") == "true"
    mkpath(output_dir)

    reads, reference = make_fixture(
        genome_length, read_length, coverage, error_rate, seed
    )
    println("[$(label)] genome=$(genome_length) read_length=$(read_length) " *
            "coverage=$(coverage) error=$(error_rate) reads=$(length(reads))")

    naive = run_arm(reads, genome_length, :none, nothing)
    println("[$(label)] naive: contigs=$(naive.n_contigs) largest=$(naive.largest) " *
            "ratio=$(round(naive.ratio; digits = 4))")

    corrected_dir = joinpath(output_dir, "corrected")
    corrected = run_arm(reads, genome_length, :iterative, corrected_dir)
    println("[$(label)] corrected: contigs=$(corrected.n_contigs) " *
            "largest=$(corrected.largest) ratio=$(round(corrected.ratio; digits = 4)) " *
            "reassembly_k=$(corrected.reassembly_k)/$(corrected.reassembly_k_ceiling)")

    corrected_fastq = joinpath(corrected_dir, "corrected.fastq")
    corrected_sequences = isfile(corrected_fastq) ?
                          read_fastq_sequences(corrected_fastq) :
                          BioSequences.LongDNA{4}[]
    raw_sequences = [FASTX.sequence(BioSequences.LongDNA{4}, record) for record in reads]

    raw_ladder = k_ladder(reads, DIAG_MAX_K)
    corrected_ladder = isempty(corrected_sequences) ? nothing :
                       k_ladder(corrected_sequences, DIAG_MAX_K)

    # Decisive test for the re-assembly-k hypothesis: assemble the SAME raw reads
    # and the SAME corrected reads with corrector=:none at a range of k, so the
    # corrector's contribution and the k choice are separated. If corrected reads
    # at the ceiling k beat raw reads at the ceiling k, correction is fine and the
    # coverage-aware k drop-down is what costs the assembly.
    k_sweep_spec = get(parsed, "k-sweep", "")
    k_sweep = DataFrames.DataFrame(
        label = String[], arm = String[], k = Int[],
        n_contigs = Int[], largest = Int[], ratio = Float64[]
    )
    if !isempty(k_sweep_spec)
        sweep_ks = parse.(Int, split(k_sweep_spec, ','))
        corrected_records = isfile(corrected_fastq) ?
                            collect(open(FASTX.FASTQ.Reader, corrected_fastq) do reader
            collect(reader)
        end) : FASTX.FASTQ.Record[]
        for (arm, arm_reads) in (("raw", reads), ("corrected", corrected_records))
            isempty(arm_reads) && continue
            for k in sweep_ks
                Random.seed!(1_042)
                assembly = Mycelia.Rhizomorph.assemble_genome(
                    deepcopy(arm_reads);
                    k = k, corrector = :none, strategy = :scalable,
                    sequencing_tech = :illumina
                )
                lengths = [length(contig) for contig in assembly.contigs]
                largest = isempty(lengths) ? 0 : maximum(lengths)
                push!(k_sweep,
                    (label, arm, k, length(assembly.contigs), largest,
                        largest / genome_length))
                println("[$(label)] k-sweep $(arm) k=$(k): " *
                        "contigs=$(length(assembly.contigs)) largest=$(largest)")
            end
        end
        CSV.write(joinpath(output_dir, "reassembly_k_sweep.csv"), k_sweep)
    end

    raw_stats = nothing
    corrected_stats = nothing
    if !skip_alignment
        score_model = BioAlignments.AffineGapScoreModel(
            BioAlignments.EDNAFULL; gap_open = -10, gap_extend = -1
        )
        reference_kmers = canonical_kmer_set(reference, DIAG_IDENTITY_KMER)
        raw_stats = read_set_stats(
            raw_sequences, reference, reference_kmers, score_model
        )
        corrected_stats = isempty(corrected_sequences) ? nothing :
                          read_set_stats(
            corrected_sequences, reference, reference_kmers, score_model
        )
        println("[$(label)] raw reads:       n=$(raw_stats.n) " *
                "mean_identity=$(round(raw_stats.mean_identity; digits = 5)) " *
                "mean_ref_kmer_frac=$(round(raw_stats.mean_reference_kmer_fraction; digits = 5))")
        if corrected_stats !== nothing
            println("[$(label)] corrected reads: n=$(corrected_stats.n) " *
                    "mean_identity=$(round(corrected_stats.mean_identity; digits = 5)) " *
                    "mean_ref_kmer_frac=$(round(corrected_stats.mean_reference_kmer_fraction; digits = 5))")
        end
    end

    summary = DataFrames.DataFrame(
        label = [label],
        genome_length = [genome_length],
        read_length = [read_length],
        coverage = [coverage],
        error_rate = [error_rate],
        seed = [seed],
        n_reads = [length(reads)],
        naive_contigs = [naive.n_contigs],
        naive_largest = [naive.largest],
        naive_ratio = [naive.ratio],
        naive_seconds = [naive.wall_seconds],
        corrected_contigs = [corrected.n_contigs],
        corrected_largest = [corrected.largest],
        corrected_ratio = [corrected.ratio],
        corrected_seconds = [corrected.wall_seconds],
        corrected_read_count = [corrected.corrected_read_count],
        reassembly_k = [corrected.reassembly_k],
        reassembly_k_ceiling = [corrected.reassembly_k_ceiling],
        reassembly_graph_reused = [corrected.reassembly_graph_reused],
        select_k_on_raw = [raw_ladder.chosen],
        select_k_on_corrected = [
            corrected_ladder === nothing ? missing : corrected_ladder.chosen
        ],
        raw_mean_identity = [
            raw_stats === nothing ? missing : raw_stats.mean_identity
        ],
        corrected_mean_identity = [
            corrected_stats === nothing ? missing : corrected_stats.mean_identity
        ],
        raw_median_identity = [
            raw_stats === nothing ? missing : raw_stats.median_identity
        ],
        corrected_median_identity = [
            corrected_stats === nothing ? missing : corrected_stats.median_identity
        ],
        raw_mean_ref_kmer_fraction = [
            raw_stats === nothing ? missing : raw_stats.mean_reference_kmer_fraction
        ],
        corrected_mean_ref_kmer_fraction = [
            corrected_stats === nothing ? missing :
            corrected_stats.mean_reference_kmer_fraction
        ],
        raw_mean_aligned_fraction = [
            raw_stats === nothing ? missing : raw_stats.mean_aligned_fraction
        ],
        corrected_mean_aligned_fraction = [
            corrected_stats === nothing ? missing :
            corrected_stats.mean_aligned_fraction
        ],
        raw_mean_length = [raw_stats === nothing ? missing : raw_stats.mean_length],
        corrected_mean_length = [
            corrected_stats === nothing ? missing : corrected_stats.mean_length
        ]
    )
    summary_path = joinpath(output_dir, "cell_summary.csv")
    CSV.write(summary_path, summary)

    ladder = DataFrames.DataFrame(
        label = String[], arm = String[], k = Int[], median_solid = Float64[]
    )
    for (k, median_solid) in zip(raw_ladder.candidate_ks, raw_ladder.medians)
        push!(ladder, (label, "raw", k, median_solid))
    end
    if corrected_ladder !== nothing
        for (k, median_solid) in
            zip(corrected_ladder.candidate_ks, corrected_ladder.medians)
            push!(ladder, (label, "corrected", k, median_solid))
        end
    end
    ladder_path = joinpath(output_dir, "k_ladder.csv")
    CSV.write(ladder_path, ladder)

    println("[$(label)] wrote $(summary_path)")
    println("[$(label)] wrote $(ladder_path)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
