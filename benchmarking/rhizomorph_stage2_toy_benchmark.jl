#!/usr/bin/env julia

import BioSequences
import CairoMakie
import CodecZlib
import CSV
import DataFrames
import Dates
import FASTX
import FileWatching
import JSON
import Mycelia
import Random
import SHA
import StableRNGs
import TOML

include(joinpath(@__DIR__, "quast_report_parsing.jl"))

const FINAL_OUTPUT_ROOT = get(ENV, "RHIZOMORPH_STAGE2_OUTPUT",
    joinpath(@__DIR__, "results", "rhizomorph_stage2_toy"))
const OUTPUT_ROOT = FINAL_OUTPUT_ROOT * ".staging"
const FIXTURE_LENGTH = parse(Int, get(ENV, "RHIZOMORPH_STAGE2_LENGTH", "8000"))
const THREADS = parse(Int, get(ENV, "RHIZOMORPH_STAGE2_THREADS", "4"))
const FIXTURE_SEED = parse(Int, get(ENV, "RHIZOMORPH_STAGE2_SEED", "43001"))
const REFERENCE_SEED = parse(Int,
    get(ENV, "RHIZOMORPH_STAGE2_REFERENCE_SEED", "43003"))
const COVERAGES = String.(
    split(get(ENV, "RHIZOMORPH_STAGE2_COVERAGES", "50x,30x,20x"), ','))
const BADREAD_IDENTITY = get(ENV, "RHIZOMORPH_STAGE2_IDENTITY", "95,99,2.5")
const DIVERGENCES = parse.(Float64,
    split(get(ENV, "RHIZOMORPH_STAGE2_DIVERGENCES", "0.01,0.02"), ','))
const RUN_MODE = Symbol(get(ENV, "RHIZOMORPH_STAGE2_MODE", "development"))
const FIXTURE_COMPLEXITY = Symbol(
    get(ENV, "RHIZOMORPH_STAGE2_COMPLEXITY", "repeat_rich"))
const STAGE2_STAGING_SENTINEL_NAME = ".rhizomorph-stage2-staging.toml"
const STAGE2_STAGING_SENTINEL_SCHEMA = "rhizomorph-stage2-staging-v2"
# The toy gate is bounded at about one hour. A four-hour stale floor avoids
# reclaiming a plausibly live slow attempt while permitting dead-owner recovery.
const STAGE2_STAGING_LOCK_STALE_SECONDS = 4 * 60 * 60
const STAGE2_REPEAT_RICH_MINIMUM_LENGTH = 5_399
const STAGE2_TOY_EVIDENCE_STATUS = "unattested_contract_only"

mutable struct Stage2StagingAttempt{T}
    attempt_id::String
    final_output_root::String
    staging_root::String
    lock_path::String
    lock::T
    active::Bool
end

function stage2_output_paths(final_output_root::AbstractString)::NamedTuple
    isempty(strip(final_output_root)) &&
        error("Stage-2 final output root must be nonempty")
    final = abspath(normpath(String(final_output_root)))
    staging = final * ".staging"
    return (; final, staging, lock = staging * ".lock")
end

function stage2_staging_sentinel_path(staging_root::AbstractString)::String
    return joinpath(staging_root, STAGE2_STAGING_SENTINEL_NAME)
end

function validate_stage2_attempt_id(attempt_id::AbstractString)::String
    value = String(attempt_id)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", value) ||
        error("Stage-2 staging attempt ID is invalid")
    return value
end

function new_stage2_attempt_id()::String
    return bytes2hex(Random.rand(Random.RandomDevice(), UInt8, 16))
end

function validate_stage2_owned_staging(
        staging_root::AbstractString,
        final_output_root::AbstractString;
        expected_attempt_id::Union{Nothing, AbstractString} = nothing,
        expected_lock_state::Union{Nothing, AbstractString} = nothing,
)::Nothing
    paths = stage2_output_paths(final_output_root)
    staging = abspath(normpath(String(staging_root)))
    staging == paths.staging ||
        error("Stage-2 staging path differs from the exact final-path sibling")
    islink(staging) && error("Stage-2 staging root must not be a symbolic link")
    isdir(staging) || error("Stage-2 staging root is not a directory")
    sentinel = stage2_staging_sentinel_path(staging)
    islink(sentinel) && error("Stage-2 staging sentinel must not be a symbolic link")
    isfile(sentinel) || error("refusing to replace unowned Stage-2 staging root")
    table = TOML.parsefile(sentinel)
    get(table, "schema", "") == STAGE2_STAGING_SENTINEL_SCHEMA ||
        error("Stage-2 staging sentinel has an unsupported schema")
    get(table, "final_output_root", "") == paths.final ||
        error("Stage-2 staging sentinel is bound to a different final path")
    get(table, "staging_root", "") == paths.staging ||
        error("Stage-2 staging sentinel is bound to a different staging path")
    get(table, "lock_path", "") == paths.lock ||
        error("Stage-2 staging sentinel is bound to a different lock path")
    attempt_id = get(table, "attempt_id", nothing)
    attempt_id isa AbstractString ||
        error("Stage-2 staging sentinel lacks an attempt ID")
    validated_attempt_id = validate_stage2_attempt_id(attempt_id)
    expected_attempt_id === nothing ||
        validated_attempt_id == validate_stage2_attempt_id(expected_attempt_id) ||
        error("Stage-2 staging sentinel attempt ID differs from the active attempt")
    lock_state = get(table, "lock_state", nothing)
    lock_state in ("active", "promoted") ||
        error("Stage-2 staging sentinel has an invalid lock state")
    expected_lock_state === nothing ||
        lock_state == String(expected_lock_state) ||
        error("Stage-2 staging sentinel lock state differs from the active attempt")
    return nothing
end

function write_stage2_staging_sentinel!(
        attempt::Stage2StagingAttempt;
        lock_state::AbstractString = "active",
)::String
    paths = stage2_output_paths(attempt.final_output_root)
    attempt.staging_root == paths.staging ||
        error("Stage-2 staging path differs from the exact final-path sibling")
    attempt.lock_path == paths.lock ||
        error("Stage-2 lock path differs from the exact staging-path sibling")
    validated_attempt_id = validate_stage2_attempt_id(attempt.attempt_id)
    lock_state in ("active", "promoted") ||
        error("Stage-2 staging lock state must be active or promoted")
    sentinel = stage2_staging_sentinel_path(attempt.staging_root)
    mktemp(attempt.staging_root) do temporary_path, stream
        TOML.print(stream, Dict(
            "schema" => STAGE2_STAGING_SENTINEL_SCHEMA,
            "final_output_root" => paths.final,
            "staging_root" => paths.staging,
            "lock_path" => paths.lock,
            "attempt_id" => validated_attempt_id,
            "lock_state" => String(lock_state),
        ))
        Base.close(stream)
        Base.Filesystem.rename(temporary_path, sentinel)
    end
    return sentinel
end

function release_stage2_staging_attempt!(
        attempt::Stage2StagingAttempt
)::Nothing
    attempt.active || return nothing
    Base.close(attempt.lock)
    attempt.active = false
    return nothing
end

function prepare_stage2_staging!(
        final_output_root::AbstractString;
        attempt_id::AbstractString = new_stage2_attempt_id(),
)::Stage2StagingAttempt
    paths = stage2_output_paths(final_output_root)
    validated_attempt_id = validate_stage2_attempt_id(attempt_id)
    mkpath(dirname(paths.lock))
    lock = FileWatching.Pidfile.trymkpidlock(
        paths.lock; stale_age = STAGE2_STAGING_LOCK_STALE_SECONDS)
    lock === false && error(
        "Stage-2 staging has an active attempt lock: $(paths.lock)",
    )
    attempt = Stage2StagingAttempt(
        validated_attempt_id,
        paths.final,
        paths.staging,
        paths.lock,
        lock,
        true,
    )
    created_staging = false
    try
        islink(paths.final) &&
            error("Stage-2 final output root must not be a symbolic link")
        ispath(paths.final) && error(
            "Stage-2 final output root already exists; refusing to replace it",
        )
        if ispath(paths.staging) || islink(paths.staging)
            validate_stage2_owned_staging(paths.staging, paths.final)
            rm(paths.staging; recursive = true)
        end
        mkdir(paths.staging)
        created_staging = true
        write_stage2_staging_sentinel!(attempt)
        return attempt
    catch
        created_staging && ispath(paths.staging) &&
            rm(paths.staging; force = true, recursive = true)
        release_stage2_staging_attempt!(attempt)
        rethrow()
    end
end

function promote_stage2_staging!(attempt::Stage2StagingAttempt)::String
    attempt.active || error("Stage-2 staging attempt lock is not active")
    try
        paths = stage2_output_paths(attempt.final_output_root)
        attempt.staging_root == paths.staging ||
            error("Stage-2 staging path differs from the active attempt")
        attempt.lock_path == paths.lock ||
            error("Stage-2 lock path differs from the active attempt")
        validate_stage2_owned_staging(
            paths.staging,
            paths.final;
            expected_attempt_id = attempt.attempt_id,
            expected_lock_state = "active",
        )
        islink(paths.final) &&
            error("Stage-2 final output root must not be a symbolic link")
        ispath(paths.final) && error(
            "Stage-2 final output root already exists; refusing to replace it",
        )
        write_stage2_staging_sentinel!(attempt; lock_state = "promoted")
        mv(paths.staging, paths.final)
        return paths.final
    finally
        release_stage2_staging_attempt!(attempt)
    end
end

function make_primary_sequence(
        rng::Random.AbstractRNG;
        fixture_length::Int = FIXTURE_LENGTH,
        fixture_complexity::Symbol = FIXTURE_COMPLEXITY,
)::String
    fixture_length > 0 || error("Stage-2 fixture length must be positive")
    fixture_complexity in (:random, :repeat_rich) ||
        error("RHIZOMORPH_STAGE2_COMPLEXITY must be random or repeat_rich")
    fixture_complexity == :repeat_rich &&
        fixture_length < STAGE2_REPEAT_RICH_MINIMUM_LENGTH && error(
        "repeat-rich Stage-2 fixture length must be at least " *
        "$(STAGE2_REPEAT_RICH_MINIMUM_LENGTH)",
    )
    sequence = collect(String(BioSequences.randdnaseq(rng, fixture_length)))
    fixture_complexity == :random && return String(sequence)
    for (offset, base) in zip(500:500:(fixture_length - 20), "ACGTACGTACGTAC")
        sequence[offset:(offset + 19)] .= base
    end
    repeat_unit = sequence[1000:1099]
    sequence[3000:3399] .= repeat(repeat_unit, 4)
    sequence[5000:5399] .= repeat(reverse(repeat_unit), 4)
    return String(sequence)
end

function mutate_haplotype(
        sequence::String,
        substitution_rate::Float64,
        seed::Int;
        n_indels::Int = 4
)::String
    rng = StableRNGs.StableRNG(seed)
    bases = collect(sequence)
    n_substitutions = round(Int, length(bases) * substitution_rate)
    substitution_positions = Random.randperm(rng, length(bases))[1:n_substitutions]
    alphabet = ['A', 'C', 'G', 'T']
    for position in substitution_positions
        alternatives = filter(base -> base != bases[position], alphabet)
        bases[position] = Random.rand(rng, alternatives)
    end
    for position in sort(Random.randperm(rng, length(bases))[1:n_indels]; rev = true)
        if Random.rand(rng, Bool)
            deleteat!(bases, position)
        else
            insert!(bases, position, Random.rand(rng, alphabet))
        end
    end
    return String(bases)
end

function write_reference(path::String, id::String, sequence::String)::String
    FASTX.FASTA.Writer(open(path, "w")) do writer
        write(writer, FASTX.FASTA.Record(id, sequence))
    end
    return path
end

function append_gzip_fastq(input::String, output::IO)::Nothing
    open(input, "r") do raw
        stream = CodecZlib.GzipDecompressorStream(raw)
        write(output, stream)
        close(stream)
    end
    return nothing
end

function make_fixture(root::String)::NamedTuple
    fixture_dir = mkpath(joinpath(root, "fixture"))
    length(COVERAGES) == 3 || error("RHIZOMORPH_STAGE2_COVERAGES requires 3 values")
    length(DIVERGENCES) == 2 ||
        error("RHIZOMORPH_STAGE2_DIVERGENCES requires 2 values")
    rng = StableRNGs.StableRNG(REFERENCE_SEED)
    primary = make_primary_sequence(rng)
    secondary = mutate_haplotype(primary, DIVERGENCES[1], REFERENCE_SEED + 1)
    tertiary = mutate_haplotype(primary, DIVERGENCES[2], REFERENCE_SEED + 2)
    truth = [
        (id = "primary", sequence = primary, coverage = COVERAGES[1], seed = FIXTURE_SEED + 101),
        (id = "secondary", sequence = secondary, coverage = COVERAGES[2], seed = FIXTURE_SEED + 102),
        (id = "tertiary", sequence = tertiary, coverage = COVERAGES[3], seed = FIXTURE_SEED + 103)
    ]
    reference_paths = String[]
    read_paths = String[]
    for haplotype in truth
        reference = joinpath(fixture_dir, "$(haplotype.id).fasta")
        reads = joinpath(fixture_dir, "$(haplotype.id).fastq.gz")
        write_reference(reference, haplotype.id, haplotype.sequence)
        Mycelia.simulate_badread_reads(;
            fasta = reference,
            quantity = haplotype.coverage,
            outfile = reads,
            identity = BADREAD_IDENTITY,
            quiet = true,
            seed = haplotype.seed)
        push!(reference_paths, reference)
        push!(read_paths, reads)
    end
    combined_reference = joinpath(fixture_dir, "truth.fasta")
    open(combined_reference, "w") do io
        for reference in reference_paths
            open(reference, "r") do input
                write(io, input)
            end
        end
    end
    mixed_fastq = joinpath(fixture_dir, "mixed.fastq")
    open(mixed_fastq, "w") do io
        for reads in read_paths
            append_gzip_fastq(reads, io)
        end
    end
    return (; truth, reference_paths, combined_reference, mixed_fastq)
end

function require_quast_metrics(report::String)::NamedTuple
    base = quast_metrics_for_report(report)
    metrics = (;
        base...,
        quast_mismatches_per_100kbp =
            parse_quast_metric(report, "# mismatches per 100 kbp"),
        quast_indels_per_100kbp =
            parse_quast_metric(report, "# indels per 100 kbp"))
    required = (
        metrics.quast_nga50,
        metrics.quast_num_misassemblies,
        metrics.quast_mismatches_per_100kbp,
        metrics.quast_indels_per_100kbp
    )
    any(ismissing, required) && error("QUAST report lacks a required Stage-2 metric: $(report)")
    return metrics
end

function evaluate_primary_segment_candidate(
        assembly::String,
        reference::String,
        output_dir::String
)::NamedTuple
    quast_dir = Mycelia.run_quast(assembly;
        outdir = output_dir, reference = reference, threads = THREADS, min_contig = 500)
    return require_quast_metrics(joinpath(quast_dir, "report.tsv"))
end

function raw_ranked_segment_candidates(
        strainy_result::NamedTuple,
        fastq::String,
        evaluation_panel_size::Int
)::Vector{<:NamedTuple}
    evaluation_panel_size > 0 || throw(
        ArgumentError("evaluation_panel_size must be positive"),
    )
    segments = Mycelia.Rhizomorph._read_gfa_records(strainy_result.strain_contigs_gfa)
    coverages = Mycelia.Rhizomorph._strainy_coverages(strainy_result.phased_unitig_info)
    support = Mycelia.Rhizomorph._read_support_coverages(segments, fastq, THREADS)
    candidates = [(;
        id,
        sequence = record.sequence,
        coverage = haskey(support, id) ? support[id] :
                   something(record.coverage, get(coverages, id, nothing)))
                  for (id, record) in segments
                  if haskey(support, id) || record.coverage !== nothing ||
                     haskey(coverages, id)]
    sort!(candidates; by = candidate ->
        (-candidate.coverage, -length(candidate.sequence), candidate.id))
    length(candidates) >= evaluation_panel_size ||
        error(
            "raw comparator arm produced only $(length(candidates)) " *
            "experimental segment candidates",
        )
    return candidates[1:evaluation_panel_size]
end

function write_single_fasta(path::String, id::String, sequence::String)::String
    return write_reference(path, id, sequence)
end

function _parse_dnadiff_discrepancy_stream(stream::IO)::NamedTuple
    total_match = nothing
    aligned_match = nothing
    identity_match = nothing
    for line in eachline(stream)
        total_match === nothing &&
            (total_match = match(r"TotalBases\s+(\d+)\s+(\d+)", line))
        aligned_match === nothing &&
            (aligned_match = match(
                r"AlignedBases\s+(\d+)\(([\d.]+)%\)\s+(\d+)\(([\d.]+)%\)",
                line,
            ))
        identity_match === nothing &&
            (identity_match = match(r"AvgIdentity\s+([\d.]+)", line))
        total_match !== nothing && aligned_match !== nothing &&
            identity_match !== nothing && break
    end
    total_match === nothing && error("dnadiff report lacks total base counts")
    aligned_match === nothing && error("dnadiff report lacks reference aligned percentage")
    identity_match === nothing && error("dnadiff report lacks average identity")
    total_ref = parse(Int, total_match.captures[1])
    total_query = parse(Int, total_match.captures[2])
    aligned_ref = parse(Int, aligned_match.captures[1])
    aligned_pct_ref = parse(Float64, aligned_match.captures[2])
    aligned_query = parse(Int, aligned_match.captures[3])
    avg_identity = parse(Float64, identity_match.captures[1])
    discrepancy = (1.0 - avg_identity / 100.0) * aligned_ref +
                  (total_ref - aligned_ref) + (total_query - aligned_query)
    return (;
        aligned_pct_ref,
        avg_identity,
        total_errors_per_100kbp = discrepancy / total_ref * 100_000.0)
end

function parse_dnadiff_discrepancy(report::String)::NamedTuple
    return _parse_dnadiff_discrepancy_stream(IOBuffer(report))
end

function parse_dnadiff_discrepancy_file(path::AbstractString)::NamedTuple
    isfile(path) || error("dnadiff report does not exist: $(path)")
    return open(path, "r") do stream
        _parse_dnadiff_discrepancy_stream(stream)
    end
end

function pair_metrics(
        candidate_id::String,
        candidate_sequence::String,
        truth_id::String,
        truth_path::String,
        output_dir::String
)::NamedTuple
    mkpath(output_dir)
    candidate_path = write_single_fasta(
        joinpath(output_dir, "$(candidate_id).fasta"), candidate_id, candidate_sequence)
    result = Mycelia.run_dnadiff(;
        reference = truth_path,
        query = candidate_path,
        outdir = joinpath(output_dir, "$(candidate_id)_vs_$(truth_id)"),
        force = true)
    return (;
        candidate_id,
        truth_id,
        parse_dnadiff_discrepancy_file(result.report)...,
    )
end

function all_permutations(values::Vector{Int})::Vector{Vector{Int}}
    length(values) == 1 && return [copy(values)]
    permutations = Vector{Vector{Int}}()
    for (index, value) in enumerate(values)
        rest = [values[i] for i in eachindex(values) if i != index]
        for suffix in all_permutations(rest)
            push!(permutations, [value; suffix])
        end
    end
    return permutations
end

function evaluate_ranked_segment_candidates(
        arm::String,
        variants::AbstractVector,
        fixture::NamedTuple,
        output_dir::String
)::DataFrames.DataFrame
    mkpath(output_dir)
    pairwise = Dict{Tuple{Int, Int}, NamedTuple}()
    for (rank, variant) in enumerate(variants)
        for truth_index in eachindex(fixture.truth)
            pairwise[(rank, truth_index)] = pair_metrics(
                "$(arm)_rank$(rank)", variant.sequence,
                fixture.truth[truth_index].id, fixture.reference_paths[truth_index], output_dir)
        end
    end
    assignment = [argmax([
        pairwise[(rank, truth_index)].aligned_pct_ref * 1000.0 +
        pairwise[(rank, truth_index)].avg_identity
        for truth_index in eachindex(fixture.truth)]) for rank in eachindex(variants)]
    return DataFrames.DataFrame([
        (arm = arm, rank = rank, candidate_id = variants[rank].id,
            truth_id = fixture.truth[assignment[rank]].id,
            expected_truth_id = fixture.truth[rank].id,
            aligned_pct_ref = pairwise[(rank, assignment[rank])].aligned_pct_ref,
            avg_identity = pairwise[(rank, assignment[rank])].avg_identity,
            total_errors_per_100kbp =
                pairwise[(rank, assignment[rank])].total_errors_per_100kbp)
        for rank in eachindex(variants)
    ])
end

function create_figure(summary::DataFrames.DataFrame, output_dir::String)::Nothing
    arms = summary.arm
    errors = summary.total_errors_per_100kbp
    nga50 = summary.nga50
    figure = CairoMakie.Figure(size = (1000, 450))
    accuracy_axis = CairoMakie.Axis(figure[1, 1];
        title = "Per-base error (lower is better)", ylabel = "errors / 100 kbp")
    contiguity_axis = CairoMakie.Axis(figure[1, 2];
        title = "Contiguity (higher is better)", ylabel = "NGA50")
    colors = [:gray55, :seagreen3]
    CairoMakie.barplot!(accuracy_axis, 1:length(arms), errors; color = colors)
    CairoMakie.barplot!(contiguity_axis, 1:length(arms), nga50; color = colors)
    for axis in (accuracy_axis, contiguity_axis)
        axis.xticks = (1:length(arms), arms)
    end
    CairoMakie.save(joinpath(output_dir, "stage2_accuracy_contiguity.svg"), figure)
    CairoMakie.save(joinpath(output_dir, "stage2_accuracy_contiguity.png"), figure)
    return nothing
end

function file_sha256(path::AbstractString)::String
    isfile(path) || error("cannot hash a missing file: $(path)")
    return open(path, "r") do stream
        bytes2hex(SHA.sha256(stream))
    end
end

function serialized_stage2_deconvolution_provenance(
        provenance::AbstractDict
)::Dict{String, Any}
    payload = Dict{String, Any}()
    for (key, value) in provenance
        key isa AbstractString || error(
            "Stage-2 deconvolution provenance keys must be strings")
        payload[String(key)] = Base.deepcopy(value)
    end
    attempt_id = get(payload, "stage2_attempt_id", nothing)
    attempt_id isa AbstractString || error(
        "Stage-2 deconvolution provenance lacks an attempt ID")
    validated_attempt_id = validate_stage2_attempt_id(attempt_id)
    get(payload, "stage2_attempt_relative_path", nothing) ==
        validated_attempt_id || error(
        "Stage-2 deconvolution provenance path is not bound to its attempt ID",
    )
    get(payload, "stage2_attempt_path_semantics", nothing) ==
        "relative-to-config-output-dir" || error(
        "Stage-2 deconvolution provenance has invalid path semantics",
    )
    return payload
end

function conda_package_version(environment::String, package::String)::String
    payload = read(`$(Mycelia.CONDA_RUNNER) list -n $(environment) $(package) --json`, String)
    records = JSON.parse(payload)
    isempty(records) && error("package $(package) is absent from conda environment $(environment)")
    return String(only(records)["version"])
end

function run_stage2_toy_benchmark!(
        staging_attempt::Stage2StagingAttempt
)::Nothing
    abspath(normpath(OUTPUT_ROOT)) == staging_attempt.staging_root ||
        error("Stage-2 benchmark output differs from its staging attempt")
    RUN_MODE in (:development, :confirmation) ||
        error("RHIZOMORPH_STAGE2_MODE must be development or confirmation")
    git_dirty = !isempty(strip(read(`git status --porcelain --untracked-files=no`, String)))
    RUN_MODE == :confirmation && git_dirty &&
        error("confirmation runs require a clean tracked worktree")
    reuse_value = get(ENV, "RHIZOMORPH_STAGE2_REUSE", "false")
    reuse_value in ("true", "false") || error(
        "RHIZOMORPH_STAGE2_REUSE must be true or false",
    )
    reuse_value == "true" && error(
        "RHIZOMORPH_STAGE2_REUSE is disabled because stale artifacts are not " *
        "bound to the current configuration and input digests",
    )
    fixture = make_fixture(OUTPUT_ROOT)

    raw_layout = Mycelia.run_metaflye(;
        fastq = fixture.mixed_fastq,
        outdir = joinpath(OUTPUT_ROOT, "raw", "metaflye"),
        genome_size = "$(FIXTURE_LENGTH)",
        read_type = "nano-raw",
        meta = true,
        min_overlap = 1000,
        iterations = 0,
        keep_haplotypes = true,
        no_alt_contigs = true,
        threads = THREADS)
    raw_strainy = Mycelia.run_strainy(;
        gfa_ref = raw_layout.graph,
        fastq = fixture.mixed_fastq,
        outdir = joinpath(OUTPUT_ROOT, "raw", "strainy"),
        read_mode = :nano,
        stage = :e2e,
        min_unitig_coverage = 3,
        unitig_split_length = 0,
        threads = THREADS)

    stage2_config = Mycelia.Rhizomorph.AssemblyConfig(;
        k = 31,
        corrector = :iterative,
        strategy = :scalable,
        sequencing_tech = :nanopore,
        output_dir = joinpath(OUTPUT_ROOT, "stage2"))
    stage2 = Mycelia.Rhizomorph.deconvolve_stage2(
        fixture.mixed_fastq, stage2_config;
        genome_size = "$(FIXTURE_LENGTH)", min_overlap = 1000,
        min_scored_fraction = 0.9, threads = THREADS)

    # The truth count is used only by this evaluator to choose a fixed-size
    # assignment panel. It is not supplied to the production candidate generator.
    raw_candidates = raw_ranked_segment_candidates(
        raw_strainy,
        fixture.mixed_fastq,
        length(fixture.truth),
    )
    raw_primary_fasta = write_single_fasta(
        joinpath(mkpath(joinpath(OUTPUT_ROOT, "raw", "primary")), "primary.fasta"),
        raw_candidates[1].id, raw_candidates[1].sequence)
    raw_metrics = evaluate_primary_segment_candidate(
        raw_primary_fasta, fixture.reference_paths[1],
        joinpath(OUTPUT_ROOT, "raw", "quast"))
    stage2_metrics = evaluate_primary_segment_candidate(
        stage2.primary_segment_candidate_fasta, fixture.reference_paths[1],
        joinpath(OUTPUT_ROOT, "stage2", "quast"))
    raw_primary_metrics = pair_metrics(
        "raw_primary", raw_candidates[1].sequence, fixture.truth[1].id,
        fixture.reference_paths[1], joinpath(OUTPUT_ROOT, "raw", "primary_dnadiff"))
    stage2_primary = only(filter(
        candidate -> candidate.likelihood_rank == 1,
        stage2.segment_candidates,
    ))
    stage2_primary_metrics = pair_metrics(
        "stage2_primary", stage2_primary.sequence, fixture.truth[1].id,
        fixture.reference_paths[1], joinpath(OUTPUT_ROOT, "stage2", "primary_dnadiff"))
    summary = DataFrames.DataFrame([
        (arm = "raw-comparator", nga50 = raw_metrics.quast_nga50,
            misassemblies = raw_metrics.quast_num_misassemblies,
            genome_fraction = raw_metrics.quast_genome_fraction,
            total_errors_per_100kbp = raw_primary_metrics.total_errors_per_100kbp),
        (arm = "rhizomorph-stage2-foundation",
            nga50 = stage2_metrics.quast_nga50,
            misassemblies = stage2_metrics.quast_num_misassemblies,
            genome_fraction = stage2_metrics.quast_genome_fraction,
            total_errors_per_100kbp = stage2_primary_metrics.total_errors_per_100kbp)
    ])

    raw_rank_table = evaluate_ranked_segment_candidates(
        "raw-comparator-abundance", raw_candidates, fixture,
        joinpath(OUTPUT_ROOT, "raw", "dnadiff"))
    n_truth = length(fixture.truth)
    length(stage2.segment_candidates) >= n_truth || error(
        "Stage-2 produced fewer segment candidates than the evaluator truth panel",
    )
    likelihood_candidates = sort(
        stage2.segment_candidates;
        by = candidate -> candidate.likelihood_rank,
    )[1:n_truth]
    abundance_candidates = sort(
        stage2.segment_candidates;
        by = candidate -> candidate.abundance_rank,
    )[1:n_truth]
    likelihood_rank_table = evaluate_ranked_segment_candidates(
        "rhizomorph-segment-candidate-likelihood",
        likelihood_candidates,
        fixture,
        joinpath(OUTPUT_ROOT, "stage2", "dnadiff_likelihood"),
    )
    abundance_rank_table = evaluate_ranked_segment_candidates(
        "rhizomorph-segment-candidate-abundance",
        abundance_candidates,
        fixture,
        joinpath(OUTPUT_ROOT, "stage2", "dnadiff_abundance"),
    )
    ranking = vcat(raw_rank_table, likelihood_rank_table, abundance_rank_table)

    tables_dir = mkpath(joinpath(OUTPUT_ROOT, "tables"))
    plots_dir = mkpath(joinpath(OUTPUT_ROOT, "plots"))
    outputs_dir = mkpath(joinpath(OUTPUT_ROOT, "outputs"))
    CSV.write(joinpath(tables_dir, "stage2_primary_summary.csv"), summary)
    CSV.write(joinpath(tables_dir, "stage2_ranked_candidate_recovery.csv"), ranking)
    cp(
        stage2.primary_segment_candidate_fasta,
        joinpath(outputs_dir, "primary_segment_candidate.fasta");
        force = true,
    )
    cp(
        stage2.likelihood_ranked_segment_candidates_fasta,
        joinpath(outputs_dir, "ranked_segment_candidates_likelihood.fasta");
        force = true,
    )
    cp(
        stage2.abundance_ranked_segment_candidates_fasta,
        joinpath(outputs_dir, "ranked_segment_candidates_abundance.fasta");
        force = true,
    )
    create_figure(summary, plots_dir)

    raw_row = summary[summary.arm .== "raw-comparator", :][1, :]
    stage2_row = summary[summary.arm .== "rhizomorph-stage2-foundation", :][1, :]
    likelihood_ranks = ranking[
        ranking.arm .== "rhizomorph-segment-candidate-likelihood",
        :,
    ]
    abundance_ranks = ranking[
        ranking.arm .== "rhizomorph-segment-candidate-abundance",
        :,
    ]
    gates = Dict(
        "nga50_no_worse" => stage2_row.nga50 >= raw_row.nga50,
        "misassemblies_no_worse" =>
            stage2_row.misassemblies <= raw_row.misassemblies,
        "genome_fraction_at_least_95pct" => stage2_row.genome_fraction >= 95.0,
        "genome_fraction_at_parity" =>
            stage2_row.genome_fraction >= raw_row.genome_fraction - 0.1,
        "per_base_error_strictly_better" =>
            stage2_row.total_errors_per_100kbp < raw_row.total_errors_per_100kbp,
        "rankings_agree_on_primary" =>
            likelihood_ranks.truth_id[1] == abundance_ranks.truth_id[1] == "primary",
        "likelihood_three_distinct_truth_matches" =>
            length(unique(likelihood_ranks.truth_id)) == 3,
        "abundance_three_distinct_truth_matches" =>
            length(unique(abundance_ranks.truth_id)) == 3,
        "truth_coverage_at_least_95pct" =>
            all(likelihood_ranks.aligned_pct_ref .>= 95.0) &&
            all(abundance_ranks.aligned_pct_ref .>= 95.0),
        "truth_identity_at_least_99pct" =>
            all(likelihood_ranks.avg_identity .>= 99.0) &&
            all(abundance_ranks.avg_identity .>= 99.0),
        "abundance_rank_exact" =>
            all(abundance_ranks.truth_id .== abundance_ranks.expected_truth_id)
    )
    contract_passed = all(values(gates))

    provenance_dir = mkpath(joinpath(OUTPUT_ROOT, "provenance"))
    provenance = Dict(
        "scale" => "toy-smoke",
        "claim_boundary" =>
            "unattested toy contract only; not a deconvolution or SOTA claim",
        "evidence_status" => STAGE2_TOY_EVIDENCE_STATUS,
        "generated_at" => string(Dates.now()),
        "staging_attempt_id" => staging_attempt.attempt_id,
        "stage2_deconvolution" =>
            serialized_stage2_deconvolution_provenance(stage2.provenance),
        "git_sha" => strip(read(`git rev-parse HEAD`, String)),
        "fixture_length" => FIXTURE_LENGTH,
        "truth_abundances" => COVERAGES,
        "badread_identity" => BADREAD_IDENTITY,
        "fixture_complexity" => String(FIXTURE_COMPLEXITY),
        "truth_divergences" => DIVERGENCES,
        "fixture_seed" => FIXTURE_SEED,
        "reference_seed" => REFERENCE_SEED,
        "badread_seeds" => FIXTURE_SEED .+ [101, 102, 103],
        "mutation_seeds" => REFERENCE_SEED .+ [1, 2],
        "threads" => THREADS,
        "run_mode" => String(RUN_MODE),
        "git_dirty" => git_dirty,
        "command" => join(Base.ARGS, " "),
        "tool_versions" => Dict(
            "julia" => string(VERSION),
            "badread" => conda_package_version("badread", "badread"),
            "flye" => conda_package_version("flye", "flye"),
            "strainy" => conda_package_version("strainy", "strainy"),
            "quast" => conda_package_version("quast", "quast"),
            "mummer" => conda_package_version("mummer", "mummer")),
        "scientific_passed" => false,
        "contract_passed" => contract_passed,
        "gates" => gates,
        "output_hashes" => Dict(
            "primary_segment_candidate" => file_sha256(
                joinpath(outputs_dir, "primary_segment_candidate.fasta")),
            "ranked_segment_candidates_likelihood" => file_sha256(joinpath(
                outputs_dir,
                "ranked_segment_candidates_likelihood.fasta",
            )),
            "ranked_segment_candidates_abundance" => file_sha256(joinpath(
                outputs_dir,
                "ranked_segment_candidates_abundance.fasta",
            )),
            "primary_summary" => file_sha256(
                joinpath(tables_dir, "stage2_primary_summary.csv")),
            "candidate_recovery" => file_sha256(
                joinpath(tables_dir, "stage2_ranked_candidate_recovery.csv")),
            "mixed_reads" => file_sha256(fixture.mixed_fastq),
            "raw_layout" => file_sha256(raw_layout.assembly),
            "corrected_reads" => file_sha256(stage2.corrected_fastq),
            "primary_reference" => file_sha256(fixture.reference_paths[1]))
    )
    open(joinpath(provenance_dir, "run.provenance.json"), "w") do io
        JSON.print(io, provenance, 2)
    end
    artifact_index = Dict(
        "tables" => [
            "tables/stage2_primary_summary.csv",
            "tables/stage2_ranked_candidate_recovery.csv"],
        "plots" => [
            "plots/stage2_accuracy_contiguity.svg",
            "plots/stage2_accuracy_contiguity.png"],
        "outputs" => [
            "outputs/primary_segment_candidate.fasta",
            "outputs/ranked_segment_candidates_likelihood.fasta",
            "outputs/ranked_segment_candidates_abundance.fasta"],
        "provenance" => "provenance/run.provenance.json")
    open(joinpath(OUTPUT_ROOT, "artifact-index.json"), "w") do io
        JSON.print(io, artifact_index, 2)
    end
    contract_passed || error("Stage-2 toy contract gate failed: $(gates)")
    promote_stage2_staging!(staging_attempt)
    println(
        "Unattested Stage-2 toy contract satisfied: $(FINAL_OUTPUT_ROOT)",
    )
    return nothing
end

function main()::Nothing
    staging_attempt = prepare_stage2_staging!(FINAL_OUTPUT_ROOT)
    try
        return run_stage2_toy_benchmark!(staging_attempt)
    finally
        release_stage2_staging_attempt!(staging_attempt)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
