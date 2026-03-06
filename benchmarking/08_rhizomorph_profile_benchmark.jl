# # Rhizomorph Evidence Profile Benchmark
#
# Benchmarks `Mycelia.Rhizomorph.assemble_genome` across memory/evidence profiles
# (`:full`, `:lightweight`, `:ultralight`) on deterministic synthetic reads and
# an optional PhiX174 slice.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import BenchmarkTools
import Dates
import JSON

include("benchmark_utils.jl")

println("=== Rhizomorph Evidence Profile Benchmark ===")
println("Start time: $(Dates.now())")

suite = BenchmarkSuite("Rhizomorph Evidence Profiles")
profiles = [:full, :lightweight, :ultralight]

function synthetic_records()
    seqs = [
        "ATGCGATGCAAT",
        "TGCGATGCAATG",
        "GCGATGCAATGC",
        "CGATGCAATGCG"
    ]

    fasta = FASTX.FASTA.Record[]
    fastq = FASTX.FASTQ.Record[]
    for (i, seq) in enumerate(seqs)
        push!(fasta, FASTX.FASTA.Record("syn_$(i)", seq))
        push!(fastq, FASTX.FASTQ.Record("syn_$(i)", seq, repeat("I", length(seq))))
    end

    return fasta, fastq
end

function phix_records()
    fasta = FASTX.FASTA.Record[]
    fastq = FASTX.FASTQ.Record[]

    try
        reference_file = Mycelia.download_genome_by_accession(accession = "NC_001422.1")
        ref = first(Mycelia.open_fastx(reference_file))
        seq = String(FASTX.sequence(ref))

        read_len = 250
        step = 100
        read_id = 1
        for start in 1:step:(length(seq) - read_len + 1)
            fragment = seq[start:(start + read_len - 1)]
            push!(fasta, FASTX.FASTA.Record("phix_$(read_id)", fragment))
            push!(fastq, FASTX.FASTQ.Record("phix_$(read_id)", fragment, repeat("I", length(fragment))))
            read_id += 1
            if read_id > 24
                break
            end
        end
    catch e
        println("Skipping PhiX dataset (download/setup failed): $(e)")
    end

    return fasta, fastq
end

function run_dataset!(suite::BenchmarkSuite, dataset_name::String, reads; k::Int = 21)
    if isempty(reads)
        return
    end

    for profile in profiles
        benchmark_result = BenchmarkTools.@benchmark Mycelia.Rhizomorph.assemble_genome(
            $reads;
            k = $k,
            graph_family = :kmer,
            memory_profile = $profile
        ) samples = 2 seconds = 20

        key = "$(dataset_name)_$(profile)"
        add_benchmark_result!(suite, key, benchmark_result)

        assembled = Mycelia.Rhizomorph.assemble_genome(
            reads;
            k = k,
            graph_family = :kmer,
            memory_profile = profile
        )
        suite.results[key]["num_contigs"] = length(assembled.contigs)
        suite.results[key]["total_contig_length"] = sum(length(c) for c in assembled.contigs)

        println("$(key): median=$(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms, memory=$(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
    end
end

synthetic_fasta, synthetic_fastq = synthetic_records()
run_dataset!(suite, "synthetic_fasta", synthetic_fasta; k = 5)
run_dataset!(suite, "synthetic_fastq", synthetic_fastq; k = 5)

if lowercase(get(ENV, "MYCELIA_INCLUDE_PHIX", "true")) == "true"
    phix_fasta, phix_fastq = phix_records()
    run_dataset!(suite, "phix_fasta", phix_fasta; k = 31)
    run_dataset!(suite, "phix_fastq", phix_fastq; k = 31)
end

results_dir = mkpath("results")
results_file = joinpath(
    results_dir,
    "08_rhizomorph_profile_benchmark_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json"
)
save_benchmark_results(suite, results_file)

format_benchmark_summary(suite)
println("Benchmark complete. End time: $(Dates.now())")
