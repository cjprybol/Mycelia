# # Mash Benchmark
#
# This benchmark measures Mash sketching and screening performance across
# different input sizes and database scales.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import BenchmarkTools
import FASTX
import Random

Random.seed!(42)

println("=== Mash Benchmark ===")

if get(ENV, "MYCELIA_RUN_EXTERNAL", "false") != "true"
    println("Mash benchmarks require external tools. Set MYCELIA_RUN_EXTERNAL=true to run.")
    exit()
end

reference_sizes = [10_000, 50_000, 100_000]
database_sizes = [3, 5, 8]

for ref_size in reference_sizes
    workdir = mktempdir()
    sketch_dir = joinpath(workdir, "sketches")

    reference_record = Mycelia.random_fasta_record(moltype=:DNA, seed=ref_size, L=ref_size)
    reference_fasta = joinpath(workdir, "reference_$(ref_size).fasta")
    Mycelia.write_fasta(outfile=reference_fasta, records=[reference_record])

    reference_sketch_path = joinpath(sketch_dir, "reference_$(ref_size).msh")

    println("\n--- Sketching reference size $(ref_size) ---")
    sketch_trial = BenchmarkTools.@benchmark Mycelia.run_mash_sketch(
        input_files=[$reference_fasta],
        outdir=$sketch_dir,
        k=21,
        s=1000
    ) setup=(rm($reference_sketch_path, force=true)) evals=1
    println(sketch_trial)

    println("Allocated bytes (sketch): $(BenchmarkTools.minimum(sketch_trial).memory)")
end

for db_size in database_sizes
    workdir = mktempdir()
    sketch_dir = joinpath(workdir, "sketches")
    screen_dir = joinpath(workdir, "screen")

    reference_files = String[]
    for i in 1:db_size
        record = Mycelia.random_fasta_record(moltype=:DNA, seed=i, L=50_000)
        ref_file = joinpath(workdir, "reference_db_$(i).fasta")
        Mycelia.write_fasta(outfile=ref_file, records=[record])
        push!(reference_files, ref_file)
    end

    reference_sketches = Mycelia.run_mash_sketch(
        input_files=reference_files,
        outdir=sketch_dir,
        k=21,
        s=1000
    ).sketches

    reference_db = Mycelia.run_mash_paste(
        out_file=joinpath(sketch_dir, "reference_db.msh"),
        in_files=reference_sketches
    )

    reader = FASTX.FASTA.Reader(open(reference_files[1]))
    record = first(reader)
    close(reader)
    reads = Mycelia.create_test_reads(String(FASTX.sequence(record)), 25, 0.01)
    reads_fastq = joinpath(workdir, "reads.fastq")
    Mycelia.write_fastq(records=reads, filename=reads_fastq)

    screen_output = joinpath(screen_dir, "$(basename(reference_db))_vs_reads_mash_screen.tsv")

    println("\n--- Screening against database size $(db_size) ---")
    screen_trial = BenchmarkTools.@benchmark Mycelia.run_mash_screen(
        reference=$reference_db,
        query=$reads_fastq,
        outdir=$screen_dir,
        winner_takes_all=true
    ) setup=(rm($screen_output, force=true)) evals=1
    println(screen_trial)

    println("Allocated bytes (screen): $(BenchmarkTools.minimum(screen_trial).memory)")
end
