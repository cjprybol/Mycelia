import Pkg
if isinteractive()
    Pkg.activate("..")
end

import BioSequences
import CSV
import DataFrames
import Dates
import FASTX
import Mycelia

include("standard_assembler_fixtures.jl")

function _selected_fixture_names(args)
    if "--list-fixtures" in args
        println(list_standard_assembler_fixtures())
        exit()
    end

    fixture_names = filter(arg -> !startswith(arg, "--"), args)
    if isempty(fixture_names)
        return sort(collect(keys(STANDARD_ASSEMBLER_FIXTURE_SPECS)))
    end
    return fixture_names
end

function _contig_metrics(contigs_fasta::String, expected_total_length::Int)
    if !isfile(contigs_fasta)
        return (
            contigs_path = contigs_fasta,
            n_contigs = missing,
            total_length = missing,
            n50 = missing,
            l50 = missing,
            longest_contig = missing,
            length_recovery = missing
        )
    end

    n_contigs, total_length, n50, l50 = Mycelia.assess_assembly_quality(contigs_fasta)
    longest_contig = 0
    for record in Mycelia.open_fastx(contigs_fasta)
        longest_contig = max(longest_contig, length(FASTX.sequence(record)))
    end

    return (
        contigs_path = contigs_fasta,
        n_contigs = n_contigs,
        total_length = total_length,
        n50 = n50,
        l50 = l50,
        longest_contig = longest_contig,
        length_recovery = expected_total_length == 0 ? missing : total_length / expected_total_length
    )
end

function _write_rhizomorph_contigs(result, contigs_fasta::String)
    records = FASTX.FASTA.Record[]
    for (contig_name, contig_sequence) in zip(result.contig_names, result.contigs)
        push!(records, FASTX.FASTA.Record(contig_name, BioSequences.LongDNA{4}(contig_sequence)))
    end
    Mycelia.write_fasta(outfile = contigs_fasta, records = records, gzip = false)
    return contigs_fasta
end

function _run_rhizomorph(fixture, outdir::String)
    records = FASTX.FASTQ.Record[]
    append!(records, collect(Mycelia.open_fastx(fixture.fastq1)))
    append!(records, collect(Mycelia.open_fastx(fixture.fastq2)))

    contigs_fasta = joinpath(outdir, "rhizomorph.contigs.fasta")
    runtime_seconds = @elapsed begin
        result = Mycelia.Rhizomorph.assemble_genome(records; k = 21)
        _write_rhizomorph_contigs(result, contigs_fasta)
    end

    return (
        status = "ok",
        runtime_seconds = runtime_seconds,
        contigs = contigs_fasta
    )
end

function _run_megahit(fixture, outdir::String, threads::Int)
    result = nothing
    runtime_seconds = @elapsed begin
        result = Mycelia.run_megahit(
            fastq1 = fixture.fastq1,
            fastq2 = fixture.fastq2,
            outdir = joinpath(outdir, "megahit"),
            k_list = "21",
            min_contig_len = 100,
            threads = threads
        )
    end
    return (
        status = "ok",
        runtime_seconds = runtime_seconds,
        contigs = result.contigs
    )
end

function _run_metaspades(fixture, outdir::String, threads::Int)
    result = nothing
    runtime_seconds = @elapsed begin
        result = Mycelia.run_metaspades(
            fastq1 = fixture.fastq1,
            fastq2 = fixture.fastq2,
            outdir = joinpath(outdir, "metaspades"),
            k_list = "21",
            threads = threads
        )
    end
    return (
        status = "ok",
        runtime_seconds = runtime_seconds,
        contigs = result.contigs
    )
end

function _run_assembler(assembler_name::String, fixture, outdir::String, threads::Int)
    try
        if assembler_name == "Rhizomorph"
            return _run_rhizomorph(fixture, outdir)
        elseif assembler_name == "MEGAHIT"
            return _run_megahit(fixture, outdir, threads)
        elseif assembler_name == "metaSPAdes"
            return _run_metaspades(fixture, outdir, threads)
        end
        error("Unsupported assembler: $(assembler_name)")
    catch error_value
        return (
            status = "failed",
            runtime_seconds = missing,
            contigs = missing,
            error = sprint(showerror, error_value)
        )
    end
end

function run_standard_fixture_benchmark(; fixture_names = _selected_fixture_names(ARGS),
        results_root = joinpath(@__DIR__, "results", "standard_fixture_assembler_comparison"),
        threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_BENCHMARK_THREADS", "2")), 2), 1, 8))
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    run_root = joinpath(results_root, timestamp)
    input_root = joinpath(run_root, "fixtures")
    mkpath(input_root)

    plan = build_standard_assembler_benchmark_plan(fixture_names = fixture_names)
    results = DataFrames.DataFrame(
        fixture = String[],
        fixture_kind = String[],
        assembler = String[],
        status = String[],
        runtime_seconds = Union{Missing, Float64}[],
        n_contigs = Union{Missing, Int}[],
        total_length = Union{Missing, Int}[],
        n50 = Union{Missing, Int}[],
        l50 = Union{Missing, Int}[],
        longest_contig = Union{Missing, Int}[],
        length_recovery = Union{Missing, Float64}[],
        contigs_path = Union{Missing, String}[],
        error = Union{Missing, String}[]
    )

    println("=== Standard Fixture Assembler Comparison ===")
    println("Fixtures: $(join(fixture_names, ", "))")
    println("Threads: $(threads)")
    println("Run root: $(run_root)")

    for fixture_name in fixture_names
        fixture_outdir = joinpath(input_root, fixture_name)
        fixture = materialize_standard_assembler_fixture(fixture_name; outdir = fixture_outdir, emit_reads = true)
        println("\nFixture: $(fixture_name)")
        println("  Reference bases: $(fixture.total_reference_bases)")

        plan_subset = DataFrames.filter(row -> row.fixture == fixture_name, plan)
        for row in DataFrames.eachrow(plan_subset)
            assembler_outdir = joinpath(run_root, fixture_name, row.assembler)
            mkpath(assembler_outdir)
            println("  Running $(row.assembler)...")
            run_result = _run_assembler(row.assembler, fixture, assembler_outdir, threads)

            metrics = if run_result.status == "ok"
                _contig_metrics(run_result.contigs, fixture.total_reference_bases)
            else
                (
                    contigs_path = missing,
                    n_contigs = missing,
                    total_length = missing,
                    n50 = missing,
                    l50 = missing,
                    longest_contig = missing,
                    length_recovery = missing
                )
            end

            DataFrames.push!(results, (
                fixture = fixture_name,
                fixture_kind = String(fixture.kind),
                assembler = row.assembler,
                status = run_result.status,
                runtime_seconds = run_result.status == "ok" ? run_result.runtime_seconds : missing,
                n_contigs = metrics.n_contigs,
                total_length = metrics.total_length,
                n50 = metrics.n50,
                l50 = metrics.l50,
                longest_contig = metrics.longest_contig,
                length_recovery = metrics.length_recovery,
                contigs_path = metrics.contigs_path,
                error = run_result.status == "ok" ? missing : run_result.error
            ))
        end
    end

    results_csv = joinpath(run_root, "assembler_comparison.csv")
    plan_csv = joinpath(run_root, "benchmark_plan.csv")
    CSV.write(results_csv, results)
    CSV.write(plan_csv, plan)

    println("\nBenchmark results written to:")
    println("  - $(results_csv)")
    println("  - $(plan_csv)")

    return results
end

run_standard_fixture_benchmark()
