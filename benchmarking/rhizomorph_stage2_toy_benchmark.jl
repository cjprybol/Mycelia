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

const OUTPUT_ROOT = get(ENV, "RHIZOMORPH_STAGE2_OUTPUT",
    joinpath(@__DIR__, "results", "rhizomorph_stage2_toy"))
const FIXTURE_LENGTH = parse(Int, get(ENV, "RHIZOMORPH_STAGE2_LENGTH", "8000"))
const THREADS = parse(Int, get(ENV, "RHIZOMORPH_STAGE2_THREADS", "4"))

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
    rng = StableRNGs.StableRNG(42001)
    primary = String(BioSequences.randdnaseq(rng, FIXTURE_LENGTH))
    secondary = mutate_haplotype(primary, 0.01, 42002)
    tertiary = mutate_haplotype(primary, 0.02, 42003)
    truth = [
        (id = "primary", sequence = primary, coverage = "50x", seed = 42101),
        (id = "secondary", sequence = secondary, coverage = "30x", seed = 42102),
        (id = "tertiary", sequence = tertiary, coverage = "20x", seed = 42103)
    ]
    reference_paths = String[]
    read_paths = String[]
    for haplotype in truth
        reference = joinpath(fixture_dir, "$(haplotype.id).fasta")
        reads = joinpath(fixture_dir, "$(haplotype.id).fastq.gz")
        if !reuse
            write_reference(reference, haplotype.id, haplotype.sequence)
            Mycelia.simulate_nanopore_reads(;
                fasta = reference,
                quantity = haplotype.coverage,
                outfile = reads,
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

function raw_ranked_variants(strainy_result, max_variants::Int)::Vector{NamedTuple}
    segments = Mycelia.Rhizomorph._read_gfa_segments(strainy_result.strain_contigs_gfa)
    coverages = Mycelia.Rhizomorph._strainy_coverages(strainy_result.phased_unitig_info)
    candidates = [(; id, sequence, coverage = coverages[id])
                  for (id, sequence) in segments if haskey(coverages, id)]
    sort!(candidates; by = candidate ->
        (-candidate.coverage, -length(candidate.sequence), candidate.id))
    length(candidates) >= max_variants ||
        error("raw SOTA arm produced only $(length(candidates)) strain candidates")
    return candidates[1:max_variants]
end

function load_ranked_variants(path::String)::Vector{NamedTuple}
    variants = NamedTuple[]
    open(FASTX.FASTA.Reader, path) do reader
        for record in reader
            push!(variants, (;
                id = FASTX.FASTA.identifier(record),
                sequence = FASTX.FASTA.sequence(String, record)))
        end
    end
    return variants
end

function write_single_fasta(path::String, id::String, sequence::String)::String
    return write_reference(path, id, sequence)
end

function pair_metrics(
        candidate_id::String,
        candidate_sequence::String,
        truth_id::String,
        truth_path::String,
        output_dir::String
)::NamedTuple
    candidate_path = write_single_fasta(
        joinpath(output_dir, "$(candidate_id).fasta"), candidate_id, candidate_sequence)
    result = Mycelia.run_dnadiff(;
        reference = truth_path,
        query = candidate_path,
        outdir = joinpath(output_dir, "$(candidate_id)_vs_$(truth_id)"),
        force = true)
    report = read(result.report, String)
    aligned_match = match(r"AlignedBases\s+\d+\(([\d.]+)%\)", report)
    identity_match = match(r"AvgIdentity\s+([\d.]+)", report)
    aligned_match === nothing && error("dnadiff report lacks reference aligned percentage")
    identity_match === nothing && error("dnadiff report lacks average identity")
    aligned = parse(Float64, aligned_match.captures[1])
    identity = parse(Float64, identity_match.captures[1])
    return (; candidate_id, truth_id, aligned_pct_ref = Float64(aligned),
        avg_identity = Float64(identity))
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
            avg_identity = pairwise[(rank, assignment[rank])].avg_identity)
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

function main()::Nothing
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
        ranked_variants_fasta = joinpath(
            OUTPUT_ROOT, "stage2", "ranked", "ranked_variants.fasta"),
        variants = load_ranked_variants(joinpath(
            OUTPUT_ROOT, "stage2", "ranked", "ranked_variants.fasta"))) :
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

    raw_metrics = evaluate_primary(raw_layout.assembly, fixture.reference_paths[1],
        joinpath(OUTPUT_ROOT, "raw", "quast"))
    stage2_metrics = evaluate_primary(stage2.primary_fasta, fixture.reference_paths[1],
        joinpath(OUTPUT_ROOT, "stage2", "quast"))
    summary = DataFrames.DataFrame([
        (arm = "raw-sota", nga50 = raw_metrics.quast_nga50,
            misassemblies = raw_metrics.quast_num_misassemblies,
            mismatches_per_100kbp = raw_metrics.quast_mismatches_per_100kbp,
            indels_per_100kbp = raw_metrics.quast_indels_per_100kbp,
            total_errors_per_100kbp = raw_metrics.quast_mismatches_per_100kbp +
                                      raw_metrics.quast_indels_per_100kbp),
        (arm = "rhizomorph-stage2", nga50 = stage2_metrics.quast_nga50,
            misassemblies = stage2_metrics.quast_num_misassemblies,
            mismatches_per_100kbp = stage2_metrics.quast_mismatches_per_100kbp,
            indels_per_100kbp = stage2_metrics.quast_indels_per_100kbp,
            total_errors_per_100kbp = stage2_metrics.quast_mismatches_per_100kbp +
                                      stage2_metrics.quast_indels_per_100kbp)
    ])

    raw_variants = raw_ranked_variants(raw_strainy, 3)
    raw_rank_table = evaluate_ranked_variants(
        "raw-sota", raw_variants, fixture, joinpath(OUTPUT_ROOT, "raw", "dnadiff"))
    stage2_rank_table = evaluate_ranked_variants(
        "rhizomorph-stage2", stage2.variants, fixture,
        joinpath(OUTPUT_ROOT, "stage2", "dnadiff"))
    ranking = vcat(raw_rank_table, stage2_rank_table)

    tables_dir = mkpath(joinpath(OUTPUT_ROOT, "tables"))
    plots_dir = mkpath(joinpath(OUTPUT_ROOT, "plots"))
    outputs_dir = mkpath(joinpath(OUTPUT_ROOT, "outputs"))
    CSV.write(joinpath(tables_dir, "stage2_primary_summary.csv"), summary)
    CSV.write(joinpath(tables_dir, "stage2_ranked_variant_recovery.csv"), ranking)
    cp(stage2.primary_fasta, joinpath(outputs_dir, "primary_consensus.fasta"); force = true)
    cp(stage2.ranked_variants_fasta,
        joinpath(outputs_dir, "ranked_variants.fasta"); force = true)
    create_figure(summary, plots_dir)

    raw_row = summary[summary.arm .== "raw-sota", :][1, :]
    stage2_row = summary[summary.arm .== "rhizomorph-stage2", :][1, :]
    stage2_ranks = ranking[ranking.arm .== "rhizomorph-stage2", :]
    gates = Dict(
        "nga50_no_worse" => stage2_row.nga50 >= raw_row.nga50,
        "misassemblies_no_worse" =>
            stage2_row.misassemblies <= raw_row.misassemblies,
        "per_base_error_strictly_better" =>
            stage2_row.total_errors_per_100kbp < raw_row.total_errors_per_100kbp,
        "three_distinct_truth_matches" => length(unique(stage2_ranks.truth_id)) == 3,
        "truth_coverage_at_least_95pct" => all(stage2_ranks.aligned_pct_ref .>= 95.0),
        "abundance_rank_exact" => all(stage2_ranks.truth_id .== stage2_ranks.expected_truth_id)
    )
    passed = all(values(gates))

    provenance_dir = mkpath(joinpath(OUTPUT_ROOT, "provenance"))
    provenance = Dict(
        "scale" => "toy-smoke",
        "claim_boundary" => "engineering proof only; not a real-genome SOTA claim",
        "generated_at" => string(Dates.now()),
        "git_sha" => strip(read(`git rev-parse HEAD`, String)),
        "fixture_length" => FIXTURE_LENGTH,
        "truth_abundances" => ["50x", "30x", "20x"],
        "badread_seeds" => [42101, 42102, 42103],
        "mutation_seeds" => [42002, 42003],
        "threads" => THREADS,
        "passed" => passed,
        "gates" => gates,
        "output_hashes" => Dict(
            "primary_consensus" => file_sha256(
                joinpath(outputs_dir, "primary_consensus.fasta")),
            "ranked_variants" => file_sha256(
                joinpath(outputs_dir, "ranked_variants.fasta")),
            "primary_summary" => file_sha256(
                joinpath(tables_dir, "stage2_primary_summary.csv")),
            "variant_recovery" => file_sha256(
                joinpath(tables_dir, "stage2_ranked_variant_recovery.csv")))
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
            "outputs/ranked_variants.fasta"],
        "provenance" => "provenance/run.provenance.json")
    open(joinpath(OUTPUT_ROOT, "artifact-index.json"), "w") do io
        JSON.print(io, artifact_index, 2)
    end
    passed || error("Stage-2 toy acceptance gate failed: $(gates)")
    println("Stage-2 toy gate PASS: $(OUTPUT_ROOT)")
    return nothing
end

main()
