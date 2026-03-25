# Real Genome Benchmark — Rhizomorph on real viral/viroid sequences
# Compares Rhizomorph assembly against reference genomes using QUAST metrics
#
# Usage: julia --project=. benchmarking/real_genome_benchmark.jl [--tier N]
#
# Tiers: 1=viroids only (fast), 2=+phages, 3=+RNA viruses (full)

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import BioSequences
import DataFrames
import CSV
import Dates
import Statistics

# === Configuration ===

TIER = let
    idx = findfirst(==("--tier"), ARGS)
    idx !== nothing && idx < length(ARGS) ? parse(Int, ARGS[idx + 1]) : 1
end

# Genome registry: (name, accession, size_bp, tier, molecule_type)
GENOMES = [
    # Tier 1: Viroids
    ("PSTVd", "NC_002030", 359, 1, "RNA"),
    ("CCCVd", "NC_003540", 246, 1, "RNA"),
    ("HSVd", "NC_001351", 297, 1, "RNA"),
    # Tier 2: Bacteriophages
    ("PhiX174", "NC_001422", 5386, 2, "DNA"),
    ("Lambda", "NC_001416", 48502, 2, "DNA"),
    ("T4", "NC_000866", 168903, 2, "DNA"),
    # Tier 3: RNA viruses
    ("SARS-CoV-2", "NC_045512", 29903, 3, "RNA")
]

K_VALUES = [11, 15, 21, 25, 31]

println("=== Real Genome Benchmark (Tier $TIER) ===")
println("Start time: $(Dates.now())")
println("K values: $K_VALUES")

selected_genomes = filter(g -> g[4] <= TIER, GENOMES)
println("Genomes: $(length(selected_genomes)) selected")

# Create working directories
benchmark_dir = mktempdir(prefix = "real_genome_benchmark_")
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
println("Working directory: $benchmark_dir")

# === Phase 1: Download reference genomes ===

function download_reference(name, accession, workdir)
    println("  Downloading $name ($accession)...")
    genome_dir = joinpath(workdir, name)
    mkpath(genome_dir)
    ref_path = Mycelia.download_genome_by_accession(
        accession = accession,
        outdir = genome_dir,
        compressed = false
    )
    if !isfile(ref_path) || filesize(ref_path) == 0
        @warn "Failed to download $name ($accession)"
        return nothing
    end
    println("    Downloaded: $ref_path ($(filesize(ref_path)) bytes)")
    return ref_path
end

println("\n--- Phase 1: Downloading reference genomes ---")
genome_refs = Dict{String, String}()
for (name, accession, size_bp, tier, mol_type) in selected_genomes
    ref = download_reference(name, accession, benchmark_dir)
    if ref !== nothing
        genome_refs[name] = ref
    end
end
println("Successfully downloaded: $(length(genome_refs))/$(length(selected_genomes)) genomes")

# === Phase 2: Rhizomorph assembly ===

function run_rhizomorph_assembly(name, ref_path, k, workdir)
    # Load reference as reads (self-assembly: can we reconstruct from the sequence itself?)
    records = collect(FASTX.FASTA.Reader(open(ref_path)))
    if isempty(records)
        @warn "No records in $ref_path"
        return nothing
    end

    outdir = joinpath(workdir, "rhizomorph_k$(k)")
    mkpath(outdir)

    t0 = time()
    local result
    try
        # use_quality_scores is auto-detected from record type (FASTA -> false)
        result = Mycelia.Rhizomorph.assemble_genome(
            records;
            k = k,
            verbose = false
        )
    catch e
        @warn "Rhizomorph failed on $name k=$k" exception = (e, catch_backtrace())
        return nothing
    end
    runtime = time() - t0

    # Write contigs to FASTA
    contigs_path = joinpath(outdir, "$(name)_k$(k)_contigs.fasta")
    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)  # contigs are String type
        end
    end

    return (
        contigs_path = contigs_path,
        n_contigs = length(result.contigs),
        total_length = sum(length.(result.contigs); init = 0),
        runtime = runtime
    )
end

println("\n--- Phase 2: Rhizomorph assembly (k-mer sweep) ---")

assembly_results = DataFrames.DataFrame(
    genome = String[],
    accession = String[],
    ref_size = Int[],
    k = Int[],
    n_contigs = Int[],
    total_length = Int[],
    runtime_s = Float64[],
    contigs_path = String[],
    ref_path = String[]
)

for (name, accession, size_bp, tier, mol_type) in selected_genomes
    if !haskey(genome_refs, name)
        continue
    end
    ref_path = genome_refs[name]
    for k in K_VALUES
        # Skip k values larger than genome
        if k > size_bp
            println("  Skipping $name k=$k (k > genome size $size_bp)")
            continue
        end
        println("  Assembling $name with k=$k...")
        res = run_rhizomorph_assembly(name, ref_path, k, joinpath(benchmark_dir, name))
        if res !== nothing
            push!(assembly_results,
                (
                    genome = name,
                    accession = accession,
                    ref_size = size_bp,
                    k = k,
                    n_contigs = res.n_contigs,
                    total_length = res.total_length,
                    runtime_s = round(res.runtime; digits = 3),
                    contigs_path = res.contigs_path,
                    ref_path = ref_path
                ))
            println("    -> $(res.n_contigs) contigs, $(res.total_length) bp, $(round(res.runtime; digits=2))s")
        end
    end
end

println("\nAssembly results: $(DataFrames.nrow(assembly_results)) runs completed")
