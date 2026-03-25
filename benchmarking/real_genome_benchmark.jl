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

# === Phase 3: Quality assessment ===

println("\n--- Phase 3: Quality assessment ---")

# Add quality columns
assembly_results[!, :genome_fraction] = fill(0.0, DataFrames.nrow(assembly_results))
assembly_results[!, :largest_contig] = fill(0, DataFrames.nrow(assembly_results))
assembly_results[!, :n50] = fill(0, DataFrames.nrow(assembly_results))
assembly_results[!, :gc_content] = fill(0.0, DataFrames.nrow(assembly_results))

for row_idx in 1:DataFrames.nrow(assembly_results)
    contigs_path = assembly_results[row_idx, :contigs_path]
    ref_size = assembly_results[row_idx, :ref_size]

    if isfile(contigs_path) && filesize(contigs_path) > 0
        try
            metrics = Mycelia.assembly_metrics(contigs_path)
            assembly_results[row_idx, :largest_contig] = metrics.largest_contig
            assembly_results[row_idx, :n50] = metrics.n50
            assembly_results[row_idx, :gc_content] = round(metrics.gc_content; digits = 3)
            # Genome fraction: total assembled / reference size
            assembly_results[row_idx, :genome_fraction] = round(
                assembly_results[row_idx, :total_length] / ref_size * 100; digits = 1)
        catch e
            @warn "Metrics failed for row $row_idx" exception = e
        end
    end
end

# === Phase 4: Results summary ===

println("\n--- Results Summary ---")

# Select display columns
display_cols = [:genome, :ref_size, :k, :n_contigs, :total_length,
    :largest_contig, :n50, :gc_content, :genome_fraction, :runtime_s]
summary_df = assembly_results[:, display_cols]
println(summary_df)

# Save results
timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
csv_path = joinpath(results_dir, "real_genome_benchmark_$(timestamp).csv")
CSV.write(csv_path, summary_df)
println("\nResults saved to: $csv_path")

# Best k per genome
println("\n--- Best k per genome (by genome fraction) ---")
for gname in unique(assembly_results.genome)
    subset = filter(r -> r.genome == gname, assembly_results)
    if DataFrames.nrow(subset) > 0
        best_idx = argmax(subset.genome_fraction)
        best = subset[best_idx, :]
        println("  $gname: k=$(best.k) -> $(best.n_contigs) contigs, $(best.genome_fraction)% coverage, $(best.runtime_s)s")
    end
end

# === Phase 5 (Optional): QUAST validation ===
# Requires QUAST installed via Bioconda. Skip gracefully if unavailable.

run_quast = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") in ["1", "true", "yes"]

if run_quast
    println("\n--- Phase 5: QUAST validation ---")
    for gname in unique(assembly_results.genome)
        subset = filter(r -> r.genome == gname, assembly_results)
        if DataFrames.nrow(subset) == 0
            continue
        end
        ref_path = first(subset.ref_path)
        # Use low min_contig for viroid-scale genomes (246-359 nt)
        ref_size = first(subset.ref_size)
        mc = max(50, ref_size ÷ 10)
        contig_files = String[r.contigs_path
                              for r in eachrow(subset) if isfile(r.contigs_path)]
        if isempty(contig_files)
            continue
        end
        println("  Running QUAST for $gname ($(length(contig_files)) assemblies, min_contig=$mc)...")
        try
            quast_dir = joinpath(benchmark_dir, gname, "quast")
            Mycelia.run_quast(
                contig_files;
                outdir = quast_dir,
                reference = ref_path,
                min_contig = mc
            )
            println("    QUAST output: $quast_dir")
        catch e
            @warn "QUAST failed for $gname" exception = e
        end
    end
else
    println("\n--- Phase 5: QUAST skipped (set MYCELIA_RUN_EXTERNAL=true to enable) ---")
end

println("\n=== Benchmark complete: $(Dates.now()) ===")
