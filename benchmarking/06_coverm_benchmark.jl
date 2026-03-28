# # CoverM Coverage Benchmark
#
# This benchmark exercises the CoverM wrappers for contig- and genome-level
# coverage/abundance profiling using small synthetic datasets. It is intended
# for manual runs on a workstation or HPC node (not CI).
#
# ## What it measures
# - End-to-end runtime for `run_coverm_contig` on N BAMs
# - End-to-end runtime for `run_coverm_genome` on the same alignment set
# - Effect of thread count on CoverM throughput
#
# ## Usage
# ```bash
# # Small scale (default)
# julia --project=. benchmarking/06_coverm_benchmark.jl
#
# # Medium / large scales
# BENCHMARK_SCALE=medium julia --project=. benchmarking/06_coverm_benchmark.jl
# BENCHMARK_SCALE=large  julia --project=. benchmarking/06_coverm_benchmark.jl
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import BenchmarkTools
import Dates
import StableRNGs

println("=== CoverM Coverage Benchmark ===")
println("Start time: $(Dates.now())")

# Configuration by scale
small_config = Dict(
    "n_samples" => 2,
    "contig_length" => 5_000,
    "read_count" => 2_500,
    "threads" => min(Sys.CPU_THREADS, 8),
    "description" => "Small scale - quick validation"
)

medium_config = Dict(
    "n_samples" => 4,
    "contig_length" => 25_000,
    "read_count" => 25_000,
    "threads" => min(Sys.CPU_THREADS, 24),
    "description" => "Medium scale - realistic coverage profiling"
)

large_config = Dict(
    "n_samples" => 6,
    "contig_length" => 50_000,
    "read_count" => 50_000,
    "threads" => max(32, Sys.CPU_THREADS),
    "description" => "Large scale - HPC throughput test"
)

config = if get(ENV, "BENCHMARK_SCALE", "small") == "medium"
    medium_config
elseif get(ENV, "BENCHMARK_SCALE", "small") == "large"
    large_config
else
    small_config
end

println("Configuration: $(config["description"]) (samples=$(config["n_samples"]), reads=$(config["read_count"]), threads=$(config["threads"]))")

# Ensure environments are present before timing
Mycelia.add_bioconda_env("coverm")
Mycelia.add_bioconda_env("minimap2")
Mycelia.add_bioconda_env("samtools")

# Workspace and reference
workdir = mkpath(joinpath(pwd(), "benchmark_coverm"))
ref_fasta = joinpath(workdir, "reference.fa")
rng = StableRNGs.StableRNG(2025)
ref_record = Mycelia.random_fasta_record(L=config["contig_length"], seed=rand(rng, 0:typemax(Int)))
Mycelia.write_fasta(outfile=ref_fasta, records=[ref_record])

# Generate BAMs via simulation + mapping
println("\n--- Generating BAMs ---")
bam_paths = String[]
for i in 1:config["n_samples"]
    reads = Mycelia.simulate_illumina_reads(
        fasta=ref_fasta,
        read_count=config["read_count"],
        len=150,
        rndSeed=rand(rng, 0:typemax(Int)),
        quiet=true
    )
    map_result = Mycelia.minimap_map(
        fasta=ref_fasta,
        fastq=reads.forward_reads,
        mapping_type="sr",
        threads=config["threads"]
    )
    run(map_result.cmd)
    sorted_bam = Mycelia.sort_bam(map_result.outfile; threads=config["threads"])
    push!(bam_paths, sorted_bam)
end
println("Generated $(length(bam_paths)) BAMs")

# Build a simple bins directory (reuse the reference as a single bin)
bins_dir = joinpath(workdir, "bins")
mkpath(bins_dir)
bin_fasta = joinpath(bins_dir, "bin1.fa")
cp(ref_fasta, bin_fasta; force=true)

println("\n--- Running CoverM (contig mode) ---")
contig_out = joinpath(workdir, "coverm_contig", "coverm_contig.tsv")
contig_time = @elapsed begin
    Mycelia.run_coverm_contig(
        bam_files=bam_paths,
        reference_fasta=ref_fasta,
        methods=["mean", "covered_fraction"],
        threads=config["threads"],
        output_tsv=contig_out,
        quiet=false
    )
end
println("Contig mode runtime: $(round(contig_time, digits=2))s (output: $(contig_out))")

println("\n--- Running CoverM (genome mode) ---")
genome_out = joinpath(workdir, "coverm_genome", "coverm_genome.tsv")
genome_time = @elapsed begin
    Mycelia.run_coverm_genome(
        bam_files=bam_paths,
        genome_directory=bins_dir,
        genome_extension="fa",
        methods=["relative_abundance", "mean_coverage"],
        threads=config["threads"],
        output_tsv=genome_out,
        quiet=false
    )
end
println("Genome mode runtime: $(round(genome_time, digits=2))s (output: $(genome_out))")

println("\nBenchmark complete at $(Dates.now()). Outputs in $(workdir)")
