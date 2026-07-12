#!/usr/bin/env julia

import BioSequences
import CairoMakie
import CodecZlib
import CSV
import DataFrames
import Dates
import FASTX
import JSON
import Mycelia
import Random
import SHA
import StableRNGs

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

function make_primary_sequence(rng::Random.AbstractRNG)::String
    sequence = collect(String(BioSequences.randdnaseq(rng, FIXTURE_LENGTH)))
    FIXTURE_COMPLEXITY == :random && return String(sequence)
    FIXTURE_COMPLEXITY == :repeat_rich ||
        error("RHIZOMORPH_STAGE2_COMPLEXITY must be random or repeat_rich")
    for (offset, base) in zip(500:500:(FIXTURE_LENGTH - 20), "ACGTACGTACGTAC")
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
        write(output, read(stream))
        close(stream)
    end
    return nothing
end

function make_fixture(root::String; reuse::Bool = false)::NamedTuple
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
        if !reuse
            write_reference(reference, haplotype.id, haplotype.sequence)
            Mycelia.simulate_badread_reads(;
                fasta = reference,
                quantity = haplotype.coverage,
                outfile = reads,
                identity = BADREAD_IDENTITY,
                quiet = true,
                seed = haplotype.seed)
        end
        push!(reference_paths, reference)
        push!(read_paths, reads)
    end
    combined_reference = joinpath(fixture_dir, "truth.fasta")
    open(combined_reference, "w") do io
        for reference in reference_paths
            write(io, read(reference))
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

function evaluate_primary(
        assembly::String,
        reference::String,
        output_dir::String
)::NamedTuple
    quast_dir = Mycelia.run_quast(assembly;
        outdir = output_dir, reference = reference, threads = THREADS, min_contig = 500)
    return require_quast_metrics(joinpath(quast_dir, "report.tsv"))
end

function raw_ranked_variants(
        strainy_result,
        fastq::String,
        max_variants::Int
)::Vector{NamedTuple}
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
    length(candidates) >= max_variants ||
        error("raw SOTA arm produced only $(length(candidates)) strain candidates")
    return candidates[1:max_variants]
end

function load_ranked_variants(
        paths::Vector{String},
        table_path::String
)::Vector{NamedTuple}
    sequences = Dict{String, String}()
    for path in paths
        open(FASTX.FASTA.Reader, path) do reader
            for record in reader
                header = FASTX.FASTA.identifier(record)
                id = last(split(header, "id="; limit = 2))
                sequences[id] = FASTX.FASTA.sequence(String, record)
            end
        end
    end
    table = CSV.read(table_path, DataFrames.DataFrame; delim = '\t')
    return [(;
        id = String(row.id),
        sequence = sequences[String(row.id)],
        likelihood_rank = Int(row.likelihood_rank),
        abundance_rank = Int(row.abundance_rank)) for row in DataFrames.eachrow(table)]
end

function write_single_fasta(path::String, id::String, sequence::String)::String
    return write_reference(path, id, sequence)
end

function parse_dnadiff_discrepancy(report::String)::NamedTuple
    total_match = match(r"TotalBases\s+(\d+)\s+(\d+)", report)
    aligned_match = match(
        r"AlignedBases\s+(\d+)\(([\d.]+)%\)\s+(\d+)\(([\d.]+)%\)", report)
    identity_match = match(r"AvgIdentity\s+([\d.]+)", report)
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
    report = read(result.report, String)
    return (; candidate_id, truth_id, parse_dnadiff_discrepancy(report)...)
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

function evaluate_ranked_variants(
        arm::String,
        variants,
        fixture,
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

function file_sha256(path::String)::String
    return bytes2hex(SHA.sha256(read(path)))
end

function conda_package_version(environment::String, package::String)::String
    payload = read(`$(Mycelia.CONDA_RUNNER) list -n $(environment) $(package) --json`, String)
    records = JSON.parse(payload)
    isempty(records) && error("package $(package) is absent from conda environment $(environment)")
    return String(only(records)["version"])
end

function main()::Nothing
    RUN_MODE in (:development, :confirmation) ||
        error("RHIZOMORPH_STAGE2_MODE must be development or confirmation")
    git_dirty = !isempty(strip(read(`git status --porcelain --untracked-files=no`, String)))
    RUN_MODE == :confirmation && git_dirty &&
        error("confirmation runs require a clean tracked worktree")
    reuse = get(ENV, "RHIZOMORPH_STAGE2_REUSE", "false") == "true"
    if !reuse
        isdir(OUTPUT_ROOT) && rm(OUTPUT_ROOT; recursive = true, force = true)
        mkpath(OUTPUT_ROOT)
    end
    fixture = reuse ? make_fixture(OUTPUT_ROOT; reuse = true) : make_fixture(OUTPUT_ROOT)

    raw_layout = reuse ? (;
        assembly = joinpath(OUTPUT_ROOT, "raw", "metaflye", "assembly.fasta"),
        graph = joinpath(OUTPUT_ROOT, "raw", "metaflye", "assembly_graph.gfa")) :
    Mycelia.run_metaflye(;
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
    raw_strainy = reuse ? (;
        strain_contigs_gfa = joinpath(OUTPUT_ROOT, "raw", "strainy", "strain_contigs.gfa"),
        phased_unitig_info = joinpath(
            OUTPUT_ROOT, "raw", "strainy", "phased_unitig_info_table.csv")) :
    Mycelia.run_strainy(;
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
    stage2 = reuse ? (;
        primary_fasta = joinpath(
            OUTPUT_ROOT, "stage2", "ranked", "primary_consensus.fasta"),
        layout_assembly = joinpath(OUTPUT_ROOT, "stage2", "metaflye", "assembly.fasta"),
        corrected_fastq = joinpath(OUTPUT_ROOT, "stage2", "corrected.fastq"),
        likelihood_ranked_variants_fasta = joinpath(
            OUTPUT_ROOT, "stage2", "ranked", "ranked_variants_likelihood.fasta"),
        abundance_ranked_variants_fasta = joinpath(
            OUTPUT_ROOT, "stage2", "ranked", "ranked_variants_abundance.fasta"),
        variants = load_ranked_variants([
            joinpath(OUTPUT_ROOT, "stage2", "ranked", "ranked_variants_likelihood.fasta"),
            joinpath(OUTPUT_ROOT, "stage2", "ranked", "ranked_variants_abundance.fasta")],
            joinpath(OUTPUT_ROOT, "stage2", "ranked", "ranked_variants.tsv"))) :
    Mycelia.Rhizomorph.deconvolve_stage2(
        fixture.mixed_fastq, stage2_config;
        genome_size = "$(FIXTURE_LENGTH)", min_overlap = 1000,
        max_variants = 3, min_scored_fraction = 0.9, threads = THREADS)

    if reuse
        for directory in (
                joinpath(OUTPUT_ROOT, "raw", "quast"),
                joinpath(OUTPUT_ROOT, "stage2", "quast"),
                joinpath(OUTPUT_ROOT, "raw", "dnadiff"),
                joinpath(OUTPUT_ROOT, "stage2", "dnadiff")
        )
            rm(directory; recursive = true, force = true)
        end
    end

    raw_variants = raw_ranked_variants(raw_strainy, fixture.mixed_fastq, 3)
    raw_primary_fasta = write_single_fasta(
        joinpath(mkpath(joinpath(OUTPUT_ROOT, "raw", "primary")), "primary.fasta"),
        raw_variants[1].id, raw_variants[1].sequence)
    raw_metrics = evaluate_primary(raw_primary_fasta, fixture.reference_paths[1],
        joinpath(OUTPUT_ROOT, "raw", "quast"))
    stage2_metrics = evaluate_primary(stage2.primary_fasta, fixture.reference_paths[1],
        joinpath(OUTPUT_ROOT, "stage2", "quast"))
    raw_primary_metrics = pair_metrics(
        "raw_primary", raw_variants[1].sequence, fixture.truth[1].id,
        fixture.reference_paths[1], joinpath(OUTPUT_ROOT, "raw", "primary_dnadiff"))
    stage2_primary = only(filter(variant -> variant.likelihood_rank == 1, stage2.variants))
    stage2_primary_metrics = pair_metrics(
        "stage2_primary", stage2_primary.sequence, fixture.truth[1].id,
        fixture.reference_paths[1], joinpath(OUTPUT_ROOT, "stage2", "primary_dnadiff"))
    summary = DataFrames.DataFrame([
        (arm = "raw-sota", nga50 = raw_metrics.quast_nga50,
            misassemblies = raw_metrics.quast_num_misassemblies,
            genome_fraction = raw_metrics.quast_genome_fraction,
            total_errors_per_100kbp = raw_primary_metrics.total_errors_per_100kbp),
        (arm = "rhizomorph-stage2", nga50 = stage2_metrics.quast_nga50,
            misassemblies = stage2_metrics.quast_num_misassemblies,
            genome_fraction = stage2_metrics.quast_genome_fraction,
            total_errors_per_100kbp = stage2_primary_metrics.total_errors_per_100kbp)
    ])

    raw_rank_table = evaluate_ranked_variants(
        "raw-sota-abundance", raw_variants, fixture,
        joinpath(OUTPUT_ROOT, "raw", "dnadiff"))
    likelihood_variants = sort(stage2.variants; by = variant -> variant.likelihood_rank)[1:3]
    abundance_variants = sort(stage2.variants; by = variant -> variant.abundance_rank)[1:3]
    likelihood_rank_table = evaluate_ranked_variants(
        "rhizomorph-stage2-likelihood", likelihood_variants, fixture,
        joinpath(OUTPUT_ROOT, "stage2", "dnadiff_likelihood"))
    abundance_rank_table = evaluate_ranked_variants(
        "rhizomorph-stage2-abundance", abundance_variants, fixture,
        joinpath(OUTPUT_ROOT, "stage2", "dnadiff_abundance"))
    ranking = vcat(raw_rank_table, likelihood_rank_table, abundance_rank_table)

    tables_dir = mkpath(joinpath(OUTPUT_ROOT, "tables"))
    plots_dir = mkpath(joinpath(OUTPUT_ROOT, "plots"))
    outputs_dir = mkpath(joinpath(OUTPUT_ROOT, "outputs"))
    CSV.write(joinpath(tables_dir, "stage2_primary_summary.csv"), summary)
    CSV.write(joinpath(tables_dir, "stage2_ranked_variant_recovery.csv"), ranking)
    cp(stage2.primary_fasta, joinpath(outputs_dir, "primary_consensus.fasta"); force = true)
    cp(stage2.likelihood_ranked_variants_fasta,
        joinpath(outputs_dir, "ranked_variants_likelihood.fasta"); force = true)
    cp(stage2.abundance_ranked_variants_fasta,
        joinpath(outputs_dir, "ranked_variants_abundance.fasta"); force = true)
    create_figure(summary, plots_dir)

    raw_row = summary[summary.arm .== "raw-sota", :][1, :]
    stage2_row = summary[summary.arm .== "rhizomorph-stage2", :][1, :]
    likelihood_ranks = ranking[ranking.arm .== "rhizomorph-stage2-likelihood", :]
    abundance_ranks = ranking[ranking.arm .== "rhizomorph-stage2-abundance", :]
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
    passed = all(values(gates))

    provenance_dir = mkpath(joinpath(OUTPUT_ROOT, "provenance"))
    provenance = Dict(
        "scale" => "toy-smoke",
        "claim_boundary" => "engineering proof only; not a real-genome SOTA claim",
        "generated_at" => string(Dates.now()),
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
        "passed" => passed,
        "gates" => gates,
        "output_hashes" => Dict(
            "primary_consensus" => file_sha256(
                joinpath(outputs_dir, "primary_consensus.fasta")),
            "ranked_variants_likelihood" => file_sha256(
                joinpath(outputs_dir, "ranked_variants_likelihood.fasta")),
            "ranked_variants_abundance" => file_sha256(
                joinpath(outputs_dir, "ranked_variants_abundance.fasta")),
            "primary_summary" => file_sha256(
                joinpath(tables_dir, "stage2_primary_summary.csv")),
            "variant_recovery" => file_sha256(
                joinpath(tables_dir, "stage2_ranked_variant_recovery.csv")),
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
            "tables/stage2_ranked_variant_recovery.csv"],
        "plots" => [
            "plots/stage2_accuracy_contiguity.svg",
            "plots/stage2_accuracy_contiguity.png"],
        "outputs" => [
            "outputs/primary_consensus.fasta",
            "outputs/ranked_variants_likelihood.fasta",
            "outputs/ranked_variants_abundance.fasta"],
        "provenance" => "provenance/run.provenance.json")
    open(joinpath(OUTPUT_ROOT, "artifact-index.json"), "w") do io
        JSON.print(io, artifact_index, 2)
    end
    passed || error("Stage-2 toy acceptance gate failed: $(gates)")
    isdir(FINAL_OUTPUT_ROOT) && rm(FINAL_OUTPUT_ROOT; recursive = true, force = true)
    mv(OUTPUT_ROOT, FINAL_OUTPUT_ROOT)
    println("Stage-2 toy gate PASS: $(FINAL_OUTPUT_ROOT)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
