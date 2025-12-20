# # Tutorial 4: Genome Assembly
#
# This tutorial covers comprehensive genome assembly approaches, including short read,
# long read, and hybrid assembly methods, with emphasis on Mycelia's third-party
# assembler orchestration and benchmarking utilities.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - Different assembly algorithms and their applications
# - Short read assembly with MEGAHIT and metaSPAdes
# - Long read assembly with Flye, Canu, and hifiasm
# - Hybrid assembly approaches combining multiple data types (Unicycler)
# - Assembly quality metrics and their interpretation
# - Error correction and polishing techniques
# - Handling repetitive sequences and structural variants
# - Assembly validation and benchmarking approaches
#
# NOTE: Mycelia Rhizomorph assembly has been moved to tutorial 13
# (tutorials/13_rhizomorph_assembly.jl) and is kept as a skeleton while
# Rhizomorph unit tests stabilize.

# ## Setup
import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import FASTX
import Random
import Statistics
import BioSequences
import Kmers

Random.seed!(42)

threads = min(Mycelia.get_default_threads(), 4)
assembly_output = "assembly_output"
mkpath(assembly_output)

# ## Part 1: Assembly Algorithm Overview
#
# Understanding different assembly approaches helps choose the right method
# for your data type and research goals.

println("=== Genome Assembly Tutorial ===")

# ### Assembly Paradigms
#
# Four main approaches to genome assembly:
# 1. de Bruijn Graph - for short reads (MEGAHIT, metaSPAdes)
# 2. Overlap-Layout-Consensus (OLC) - for long reads (Canu)
# 3. String Graph - for long accurate reads (hifiasm, Flye)
# 4. Probabilistic Assembly - Mycelia Rhizomorph (see tutorial 13)

println("Assembly Algorithm Comparison:")
println("de Bruijn Graph:")
println("  - Best for: Short reads (Illumina)")
println("  - Tools: MEGAHIT, metaSPAdes, SPAdes")
println("  - Strengths: Efficient, handles high coverage")
println("  - Weaknesses: Struggles with repeats, requires error correction")
println()
println("OLC (Overlap-Layout-Consensus):")
println("  - Best for: Long reads (PacBio, Nanopore)")
println("  - Tools: Canu, Miniasm")
println("  - Strengths: Handles repeats, intuitive approach")
println("  - Weaknesses: Computationally expensive, error-sensitive")
println()
println("String Graph:")
println("  - Best for: Long accurate reads (HiFi)")
println("  - Tools: hifiasm, Flye")
println("  - Strengths: Efficient, haplotype-aware, handles complexity")
println("  - Weaknesses: Requires high-quality reads")
println()
println("Probabilistic Assembly (Mycelia Rhizomorph):")
println("  - Best for: Any read type with error correction")
println("  - Tools: Rhizomorph string graph + Viterbi")
println("  - Strengths: Handles errors probabilistically, adaptable")
println("  - Weaknesses: Computationally intensive for large genomes")
println("  - See tutorial 13 for Rhizomorph assembly details")

# ## Part 2: Data Preparation for Assembly
#
# Proper data preparation is crucial for successful assembly

println("\n=== Data Preparation ===")

# ### Simulating Multi-Platform Data
#
# Create synthetic data for comprehensive assembly testing

println("--- Generating Test Data ---")

# Create a synthetic genome with known structure
reference_size = 50_000  ## 50 kb for demonstration
reference_genome = Mycelia.random_fasta_record(moltype=:DNA, seed=1, L=reference_size)

println("Reference genome: $(reference_size) bp")

# Simulate different read types
short_read_params = Dict(
    "coverage" => 30,
    "read_length" => 150,
    "error_rate" => 0.001,
    "description" => "Illumina short reads"
)

long_read_params = Dict(
    "coverage" => 20,
    "read_length" => 10_000,
    "error_rate" => 0.05,
    "description" => "Nanopore long reads"
)

hifi_params = Dict(
    "coverage" => 15,
    "read_length" => 15_000,
    "error_rate" => 0.001,
    "description" => "PacBio HiFi reads"
)

println("Simulating read types:")
for (name, params) in [("Short reads", short_read_params),
                      ("Long reads", long_read_params),
                      ("HiFi reads", hifi_params)]
    println("  $(name): $(params["coverage"])x coverage, $(params["read_length"]) bp, $(params["error_rate"] * 100)% error")
end

# Write test data
reference_file = joinpath(assembly_output, "reference_genome.fasta")
Mycelia.write_fasta(outfile=reference_file, records=[reference_genome])

# Generate paired-end short reads (Illumina)
short_reads_base = joinpath(assembly_output, "short_reads")
short_sim = Mycelia.simulate_illumina_reads(
    fasta=reference_file,
    coverage=short_read_params["coverage"],
    outbase=short_reads_base,
    read_length=short_read_params["read_length"],
    rndSeed=42,
    quiet=true
)
short_reads_r1 = short_sim.forward_reads
short_reads_r2 = short_sim.reverse_reads

# Generate long reads for Flye/Canu (Nanopore)
long_reads_gz = Mycelia.simulate_nanopore_reads(
    fasta=reference_file,
    quantity="$(long_read_params["coverage"])x",
    outfile=joinpath(assembly_output, "long_reads.fastq.gz"),
    quiet=true
)
long_reads_file = joinpath(assembly_output, "long_reads.fastq")
if !isfile(long_reads_file)
    run(pipeline(`gzip -dc $(long_reads_gz)`, long_reads_file))
end

# Generate HiFi reads for hifiasm (PacBio HiFi)
hifi_reads_gz = Mycelia.simulate_pacbio_reads(
    fasta=reference_file,
    quantity="$(hifi_params["coverage"])x",
    outfile=joinpath(assembly_output, "hifi_reads.fastq.gz"),
    quiet=true
)
hifi_reads_file = joinpath(assembly_output, "hifi_reads.fastq")
if !isfile(hifi_reads_file)
    run(pipeline(`gzip -dc $(hifi_reads_gz)`, hifi_reads_file))
end

# Generate error-prone reads for polishing demonstrations
error_reads_gz = Mycelia.simulate_very_bad_reads(
    fasta=reference_file,
    quantity="10x",
    outfile=joinpath(assembly_output, "error_reads.fastq.gz"),
    quiet=true
)
error_reads_file = joinpath(assembly_output, "error_reads.fastq")
if !isfile(error_reads_file)
    run(pipeline(`gzip -dc $(error_reads_gz)`, error_reads_file))
end

# ### Read Statistics and Quality Assessment
#
# Analyze read characteristics before assembly

println("--- Read Analysis ---")

short_r1_stats = Mycelia.summarize_fastq("Short reads R1", short_reads_r1, reference_size)
short_r2_stats = Mycelia.summarize_fastq("Short reads R2", short_reads_r2, reference_size)
long_read_stats = Mycelia.summarize_fastq("Long reads", long_reads_file, reference_size)
hifi_stats = Mycelia.summarize_fastq("HiFi reads", hifi_reads_file, reference_size)
error_read_stats = Mycelia.summarize_fastq("Error-prone reads", error_reads_file, reference_size)

# ## Part 3: Multi-Platform Assembly Approaches
#
# Comprehensive coverage of short read, long read, and hybrid assembly

println("\n=== Multi-Platform Assembly Approaches ===")

# ### Short Read Assembly
#
# MEGAHIT and metaSPAdes for short read data

println("--- Short Read Assembly ---")

# Example parameters for short read assembly
short_read_assembly_params = Dict(
    "megahit_k_list" => "21,29,39,59,79,99,119,141",
    "metaspades_k_list" => "21,33,55,77",
    "min_contig_len" => 200,
    "threads" => threads
)

println("Short read assembly parameters:")
for (param, value) in short_read_assembly_params
    println("  $param: $value")
end

short_read_runs = []

megahit_result, megahit_runtime = Mycelia.run_assembler("MEGAHIT") do
    Mycelia.run_megahit(
        fastq1=short_reads_r1,
        fastq2=short_reads_r2,
        outdir=joinpath(assembly_output, "megahit_short"),
        min_contig_len=short_read_assembly_params["min_contig_len"],
        k_list=short_read_assembly_params["megahit_k_list"],
        threads=threads
    )
end
if megahit_result !== nothing
    push!(short_read_runs, (name="MEGAHIT", assembly=megahit_result.contigs, outdir=megahit_result.outdir, runtime=megahit_runtime, reads=short_reads_r1))
end

metaspades_result, metaspades_runtime = Mycelia.run_assembler("metaSPAdes") do
    Mycelia.run_metaspades(
        fastq1=short_reads_r1,
        fastq2=short_reads_r2,
        outdir=joinpath(assembly_output, "metaspades_short"),
        k_list=short_read_assembly_params["metaspades_k_list"],
        threads=threads
    )
end
if metaspades_result !== nothing
    push!(short_read_runs, (name="metaSPAdes", assembly=metaspades_result.contigs, outdir=metaspades_result.outdir, runtime=metaspades_runtime, reads=short_reads_r1))
end

println("Short read assembly metrics:")
for run in short_read_runs
    metrics = Mycelia.assembly_metrics(run.assembly)
    if metrics === nothing
        println("  $(run.name): assembly file missing")
        continue
    end
    output_size = Mycelia.dir_size(run.outdir)
    println("  $(run.name): contigs=$(metrics.n_contigs), N50=$(metrics.n50), L50=$(metrics.l50), total=$(metrics.total_length) bp, output_size=$(round(output_size / 1e6, digits=2)) MB")
end

# ### Long Read Assembly
#
# Flye, Canu, and hifiasm for long read data

println("--- Long Read Assembly ---")

# Example parameters for long read assembly
long_read_assembly_params = Dict(
    "genome_size" => "50k",
    "flye_read_type" => "nano-hq",
    "canu_read_type" => "nanopore",
    "hifiasm_mode" => "primary",
    "threads" => threads
)

println("Long read assembly parameters:")
for (param, value) in long_read_assembly_params
    println("  $param: $value")
end

long_read_runs = []

flye_result, flye_runtime = Mycelia.run_assembler("Flye") do
    Mycelia.run_flye(
        fastq=long_reads_file,
        outdir=joinpath(assembly_output, "flye_long"),
        genome_size=long_read_assembly_params["genome_size"],
        read_type=long_read_assembly_params["flye_read_type"],
        threads=threads
    )
end
if flye_result !== nothing
    push!(long_read_runs, (name="Flye", assembly=flye_result.assembly, outdir=flye_result.outdir, runtime=flye_runtime, reads=long_reads_file))
end

canu_result, canu_runtime = Mycelia.run_assembler("Canu") do
    Mycelia.run_canu(
        fastq=long_reads_file,
        outdir=joinpath(assembly_output, "canu_long"),
        genome_size=long_read_assembly_params["genome_size"],
        read_type=long_read_assembly_params["canu_read_type"],
        stopOnLowCoverage=5,
        threads=threads
    )
end
if canu_result !== nothing
    push!(long_read_runs, (name="Canu", assembly=canu_result.assembly, outdir=canu_result.outdir, runtime=canu_runtime, reads=long_reads_file))
end

hifiasm_result, hifiasm_runtime = Mycelia.run_assembler("hifiasm") do
    Mycelia.run_hifiasm(
        fastq=hifi_reads_file,
        outdir=joinpath(assembly_output, "hifiasm_long"),
        bloom_filter=0,
        threads=threads
    )
end
hifiasm_contigs = Mycelia.hifiasm_primary_contigs(hifiasm_result)
if hifiasm_contigs !== nothing
    push!(long_read_runs, (name="hifiasm", assembly=hifiasm_contigs, outdir=hifiasm_result.outdir, runtime=hifiasm_runtime, reads=hifi_reads_file))
else
    if hifiasm_result !== nothing
        println("  hifiasm produced graph outputs; primary contigs not found at expected .p_ctg.fa path")
    end
end

println("Long read assembly contiguity:")
for run in long_read_runs
    metrics = Mycelia.assembly_metrics(run.assembly)
    if metrics === nothing
        println("  $(run.name): assembly file missing")
        continue
    end
    qv = Mycelia.compute_merqury_qv(run.assembly, run.reads)
    println("  $(run.name): contigs=$(metrics.n_contigs), N50=$(metrics.n50), largest=$(metrics.largest_contig), QV=$(qv)")
end

# ### Hybrid Assembly
#
# Unicycler combining short and long reads

println("--- Hybrid Assembly ---")

# Example parameters for hybrid assembly
hybrid_params = Dict(
    "short_read_accuracy" => 0.99,
    "long_read_accuracy" => 0.90,
    "bridging_mode" => "conservative",
    "threads" => threads
)

println("Hybrid assembly parameters:")
for (param, value) in hybrid_params
    println("  $param: $value")
end

hybrid_runs = []

unicycler_result, unicycler_runtime = Mycelia.run_assembler("Unicycler") do
    Mycelia.run_unicycler(
        short_1=short_reads_r1,
        short_2=short_reads_r2,
        long_reads=long_reads_file,
        outdir=joinpath(assembly_output, "unicycler_hybrid"),
        threads=threads
    )
end
if unicycler_result !== nothing
    push!(hybrid_runs, (name="Unicycler", assembly=unicycler_result.assembly, outdir=unicycler_result.outdir, runtime=unicycler_runtime, reads=short_reads_r1))
end

println("Hybrid assembly comparison:")
for run in hybrid_runs
    metrics = Mycelia.assembly_metrics(run.assembly)
    if metrics === nothing
        println("  $(run.name): assembly file missing")
        continue
    end
    output_size = Mycelia.dir_size(run.outdir)
    println("  $(run.name): contigs=$(metrics.n_contigs), N50=$(metrics.n50), largest=$(metrics.largest_contig), output_size=$(round(output_size / 1e6, digits=2)) MB")
end

if !isempty(short_read_runs) && !isempty(hybrid_runs)
    short_metrics = Mycelia.assembly_metrics(short_read_runs[1].assembly)
    hybrid_metrics = Mycelia.assembly_metrics(hybrid_runs[1].assembly)
    if short_metrics !== nothing && hybrid_metrics !== nothing
        println("  Scaffolding improvement (N50): $(short_metrics.n50) -> $(hybrid_metrics.n50)")
    end
end

println("Assembly approaches comparison completed...")

println("--- Mycelia Rhizomorph Assembly ---")
println("Rhizomorph assembly workflows are covered in tutorial 13 (skeleton) and are intentionally omitted here.")

# ## Part 4: Assembly Quality Assessment
#
# Comprehensive evaluation of assembly quality

println("\n=== Assembly Quality Assessment ===")

# ### Basic Assembly Statistics
#
# Calculate fundamental assembly metrics

println("--- Basic Statistics ---")

primary_assembly = Mycelia.select_first_existing([
    isempty(hybrid_runs) ? nothing : hybrid_runs[1].assembly,
    isempty(long_read_runs) ? nothing : long_read_runs[1].assembly,
    isempty(short_read_runs) ? nothing : short_read_runs[1].assembly
])

assembly_stats = Mycelia.assembly_metrics(primary_assembly)
if assembly_stats === nothing
    println("No assembly file available for basic statistics.")
else
    println("Assembly Statistics:")
    println("  total_length: $(assembly_stats.total_length)")
    println("  n_contigs: $(assembly_stats.n_contigs)")
    println("  n50: $(assembly_stats.n50)")
    println("  n90: $(assembly_stats.n90)")
    println("  l50: $(assembly_stats.l50)")
    println("  l90: $(assembly_stats.l90)")
    println("  largest_contig: $(assembly_stats.largest_contig)")
    println("  gc_content: $(round(assembly_stats.gc_content, digits=3))")
end

# ### Advanced Quality Metrics
#
# More sophisticated quality assessment

println("--- Advanced Quality Metrics ---")

busco_outdir = nothing
if primary_assembly !== nothing
    try
        busco_outdir = Mycelia.run_busco(primary_assembly; threads=threads, outdir=joinpath(assembly_output, "busco_results"))
        println("BUSCO completed: $(busco_outdir)")
    catch e
        @warn "BUSCO failed" exception=e
    end
end

merqury_qv = Mycelia.compute_merqury_qv(primary_assembly, short_reads_r1)
println("Merqury-style QV (k=21): $(merqury_qv)")

println("LAI (LTR Assembly Index): not computed in this tutorial (no built-in wrapper yet)")

if assembly_stats !== nothing && busco_outdir !== nothing
    println("Contiguity vs completeness: compare N50 with BUSCO report in $(busco_outdir)")
end

# ### Comparison with Reference
#
# Validate assembly against known reference

println("--- Reference Comparison ---")

mummer_outdir = nothing
alignment_summary = nothing
if primary_assembly !== nothing
    try
        mummer_outdir = Mycelia.run_mummer(reference_file, primary_assembly; outdir=joinpath(assembly_output, "mummer_compare"), prefix="assembly")
        coords_file = joinpath(mummer_outdir, "assembly.coords")
        alignments = Mycelia.parse_mummer_coords(coords_file)
        if !isempty(alignments)
            total_aligned = sum(a.len_ref for a in alignments)
            weighted_identity = sum(a.len_ref * a.identity for a in alignments) / total_aligned
            coverage = total_aligned / reference_size * 100
            inversions = count(a -> a.start_query > a.end_query, alignments)
            alignment_summary = (;weighted_identity, coverage, inversions, alignment_blocks=length(alignments))
            println("Reference alignment: identity=$(round(weighted_identity, digits=2))%, coverage=$(round(coverage, digits=2))%, inversions=$(inversions)")
        end
    catch e
        @warn "MUMmer comparison failed" exception=e
    end
end

# ## Part 5: Assembly Polishing and Error Correction
#
# Improve assembly accuracy through polishing

println("\n=== Assembly Polishing ===")

# ### Consensus Polishing
#
# Use original reads to polish assembly

println("--- Consensus Polishing ---")

polish_target = Mycelia.select_first_existing([
    isempty(long_read_runs) ? nothing : long_read_runs[1].assembly,
    primary_assembly
])

apollo_result = nothing
if polish_target !== nothing
    try
        apollo_result = Mycelia.run_apollo(polish_target, long_reads_file; outdir=joinpath(assembly_output, "apollo_polish"))
        println("Apollo polishing output: $(apollo_result.polished_assembly)")
    catch e
        @warn "Apollo polishing failed" exception=e
    end
end

homopolish_result = nothing
if apollo_result !== nothing
    try
        homopolish_result = Mycelia.run_homopolish(apollo_result.polished_assembly, long_reads_file; outdir=joinpath(assembly_output, "homopolish"), threads=threads)
        println("Homopolish output: $(homopolish_result.polished_assembly)")
    catch e
        @warn "Homopolish failed" exception=e
    end
end

# ### Structural Error Correction
#
# Fix larger structural errors

println("--- Structural Correction ---")

if alignment_summary !== nothing
    println("Structural signals: alignment_blocks=$(alignment_summary.alignment_blocks), inversions=$(alignment_summary.inversions)")
    if alignment_summary.coverage < 98
        println("  Potential large gaps or structural variants detected (coverage < 98%)")
    end
else
    println("Structural correction requires reference alignments (see MUMmer comparison).")
end

# ## Part 6: Handling Assembly Challenges
#
# Address common assembly difficulties

println("\n=== Assembly Challenges ===")

# ### Repetitive Sequences
#
# Strategies for handling repeats

println("--- Repetitive Sequences ---")

if primary_assembly !== nothing
    repeat_k = 9
    kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{repeat_k}, primary_assembly)
    top_kmers = sort(collect(kmer_counts), by=x -> x[2], rev=true)
    println("Top repeated k-mers (k=$(repeat_k)):")
    for (kmer, count) in Iterators.take(top_kmers, min(5, length(top_kmers)))
        println("  $(kmer): $(count)")
    end
else
    println("Repeat analysis skipped (no assembly available)")
end

# ### Heterozygosity
#
# Handle diploid and polyploid genomes

println("--- Heterozygosity ---")

hap1_seq = FASTX.sequence(reference_genome)
mut_rng = Random.MersenneTwister(123)
hap2_seq = Mycelia.mutate_dna_substitution_fraction(hap1_seq; fraction=0.01, rng=mut_rng)
heterozygous_sites = count(pair -> pair[1] != pair[2], zip(hap1_seq, hap2_seq))
heterozygosity_rate = heterozygous_sites / length(hap1_seq)
println("Simulated heterozygosity rate: $(round(heterozygosity_rate * 100, digits=2))%")

haplotype_fasta = joinpath(assembly_output, "diploid_haplotypes.fasta")
Mycelia.write_fasta(
    outfile=haplotype_fasta,
    records=[
        FASTX.FASTA.Record("haplotype_1", hap1_seq),
        FASTX.FASTA.Record("haplotype_2", hap2_seq)
    ]
)

println("Haplotype-aware assembly recommendation: use hifiasm outputs for phased contigs where available.")

# ### Contamination
#
# Detect and remove contaminating sequences

println("--- Contamination Detection ---")

# Contamination screening options to consider:
# - Mycelia.run_checkm / Mycelia.run_checkm2 for bacterial MAGs
# - Mycelia.run_checkv for viral assemblies
# - NCBI FCS-GX (https://github.com/ncbi/fcs) TODO: add Mycelia.run_ncbi_fcs wrapper and example (bioconda: ncbi-fcs-gx)
# - Read-based taxonomic screens (e.g. Kraken2/Bracken) when wrappers are available

contaminant_candidates = Mycelia.contig_gc_outliers(primary_assembly)
if isempty(contaminant_candidates)
    println("No GC outlier contigs detected (using 2 SD threshold).")
else
    println("Potential contaminant contigs (GC outliers):")
    for candidate in contaminant_candidates
        println("  $(candidate.id): length=$(candidate.length), GC=$(round(candidate.gc, digits=3))")
    end
end

# ## Part 7: Assembly Visualization and Exploration
#
# Create plots and visualizations for assembly analysis

println("\n=== Assembly Visualization ===")

# ### Contig Size Distribution
#
# Visualize assembly contiguity

println("--- Contig Visualization ---")

if primary_assembly !== nothing
    contig_lengths = open(FASTX.FASTA.Reader, primary_assembly) do reader
        [length(FASTX.sequence(record)) for record in reader]
    end
    sorted_lengths = sort(contig_lengths, rev=true)
    cumulative_lengths = cumsum(sorted_lengths)
    println("Top contig lengths: $(sorted_lengths[1:min(5, length(sorted_lengths))])")
    println("Cumulative lengths (first 5): $(cumulative_lengths[1:min(5, length(cumulative_lengths))])")
    println("Use Plots.jl to visualize these distributions as histograms or cumulative curves.")
else
    println("No assembly available for visualization.")
end

# ### Dot Plots
#
# Compare assemblies or validate against reference

println("--- Dot Plot Analysis ---")

if mummer_outdir !== nothing
    delta_file = joinpath(mummer_outdir, "assembly.delta")
    try
        dotplot_path = Mycelia.run_mummer_plot(delta_file; outdir=mummer_outdir, prefix="assembly_plot")
        println("Dot plot generated: $(dotplot_path)")
    catch e
        @warn "Dot plot generation failed" exception=e
    end
else
    println("Dot plot requires MUMmer alignment output.")
end

# ## Part 8: Assembly Benchmarking
#
# Compare different assembly approaches

println("\n=== Assembly Benchmarking ===")

# ### Multi-Assembler Comparison
#
# Compare multiple assembly tools

println("--- Multi-Assembler Comparison ---")

all_runs = vcat(short_read_runs, long_read_runs, hybrid_runs)
benchmark_results = []
for run in all_runs
    metrics = Mycelia.assembly_metrics(run.assembly)
    if metrics === nothing
        continue
    end
    qv = Mycelia.compute_merqury_qv(run.assembly, run.reads)
    push!(benchmark_results, (name=run.name, n50=metrics.n50, total_length=metrics.total_length, runtime=run.runtime, qv=qv))
end

if isempty(benchmark_results)
    println("No benchmark results available.")
else
    sorted_results = sort(benchmark_results, by=x -> x.n50, rev=true)
    println("Assembler ranking by N50:")
    for result in sorted_results
        println("  $(result.name): N50=$(result.n50), total=$(result.total_length), runtime=$(result.runtime), QV=$(result.qv)")
    end
end

# ### Parameter Optimization
#
# Optimize assembly parameters

println("--- Parameter Optimization ---")

megahit_k_options = ["21,33,55", "21,29,39,59"]
megahit_grid_results = []
for k_list in megahit_k_options
    tag = replace(k_list, "," => "_")
    outdir = joinpath(assembly_output, "megahit_k_$(tag)")
    result, runtime = Mycelia.run_assembler("MEGAHIT k=$(k_list)") do
        Mycelia.run_megahit(
            fastq1=short_reads_r1,
            fastq2=short_reads_r2,
            outdir=outdir,
            min_contig_len=short_read_assembly_params["min_contig_len"],
            k_list=k_list,
            threads=threads
        )
    end
    if result !== nothing
        metrics = Mycelia.assembly_metrics(result.contigs)
        if metrics !== nothing
            push!(megahit_grid_results, (k_list=k_list, n50=metrics.n50, runtime=runtime))
        end
    end
end

if isempty(megahit_grid_results)
    println("Parameter optimization skipped (no MEGAHIT results).")
else
    best = sort(megahit_grid_results, by=x -> x.n50, rev=true)[1]
    println("Best MEGAHIT k-list by N50: $(best.k_list) (N50=$(best.n50))")
end

# ## Part 9: Best Practices and Recommendations
#
# Guidelines for successful genome assembly

println("\n=== Best Practices ===")

println("Data Requirements:")
println("- HiFi: 20-30x coverage minimum")
println("- Read N50 > 10 kb preferred")
println("- Low contamination levels")
println("- Balanced coverage distribution")
println()
println("Assembly Strategy:")
println("- Start with hifiasm for HiFi data")
println("- Use haplotype-aware mode for diploids")
println("- Validate with multiple quality metrics")
println("- Polish with original reads")
println()
println("Quality Control:")
println("- Check BUSCO completeness (>90% for eukaryotes)")
println("- Validate N50 vs genome size expectations")
println("- Examine contig count (fewer is better)")
println("- Compare with related genomes")

# ## Summary
println("\n=== Assembly Summary ===")
println("✓ Understanding assembly algorithms and their applications")
println("✓ Short read assembly with MEGAHIT and metaSPAdes")
println("✓ Long read assembly with Flye, Canu, and hifiasm")
println("✓ Hybrid assembly approaches with Unicycler")
println("✓ Comprehensive quality assessment techniques")
println("✓ Assembly polishing and error correction")
println("✓ Handling repetitive sequences and heterozygosity")
println("✓ Visualization and benchmarking approaches")
println("✓ Rhizomorph assembly workflows moved to tutorial 13")

# Cleanup
if isdir(assembly_output)
    rm(assembly_output, recursive=true, force=true)
end

println("\nNext: Tutorial 5 - Assembly Validation")

nothing
