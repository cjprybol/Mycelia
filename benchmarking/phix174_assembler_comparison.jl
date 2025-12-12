# PhiX174 Assembler Comparison Benchmark
# Head-to-head comparison of 8 third-party assemblers on PhiX174 reference genome
# Set `--skip-busco` or `MYCELIA_SKIP_BUSCO=true` to disable BUSCO in constrained environments (CI); BUSCO runs by default for benchmarking.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import BioSequences
import Statistics
import Dates
import DataFrames
import CSV
import Plots

# Optional: skip BUSCO for constrained environments (CI). Default: run BUSCO.
skip_busco = "--skip-busco" in ARGS || lowercase(get(ENV, "MYCELIA_SKIP_BUSCO", "false")) in ["1", "true", "yes"]
is_ci = haskey(ENV, "CI") || haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "TRAVIS") || haskey(ENV, "CIRCLECI")

println("=== PhiX174 Assembler Comparison ===")
println("Start time: $(Dates.now())")

# Create temporary directory for PhiX174 data
phix_dir = mktempdir(prefix="phix174_benchmark_")
    println("PhiX174 benchmark directory: $phix_dir")

    # Create results directory
    results_dir = "results"
    mkpath(results_dir)

    # Simple resource capture helpers
    function _time_and_mem(f)
        t0 = time()
        # Memory is not portable across platforms without extra deps; keep placeholder
        result = f()
        return (result=result, runtime=time() - t0)
    end

try
    # Download PhiX174 reference using Mycelia function
    println("\nDownloading PhiX174 reference...")
    phix_ref = Mycelia.download_genome_by_accession(
        accession="NC_001422",
        outdir=phix_dir,
        compressed=false
    )
    println("  Reference downloaded: $phix_ref")

    # Load reference sequence to get length
    ref_records = collect(FASTX.FASTA.Reader(open(phix_ref)))
    ref_seq = FASTX.sequence(BioSequences.LongDNA{4}, ref_records[1])
    ref_length = length(ref_seq)
    println("  Reference length: $ref_length bp")

    # Generate simulated reads (lighter in CI to reduce memory/time)
    println("\nGenerating simulated reads...")

    short_cov = is_ci ? 30 : 100
    long_cov = is_ci ? "20x" : "50x"

    # Short reads: coverage scaled for CI
    println("  Illumina reads ($(short_cov)x coverage)...")
    Mycelia.simulate_illumina_reads(
        fasta=phix_ref,
        coverage=short_cov,
        outbase=joinpath(phix_dir, "phix174_illumina"),
        read_length=150,
        mflen=300,
        seqSys="HS25",
        paired=true,
        errfree=false
    )
    # ART outputs files as outbase1.fq and outbase2.fq
    short_fastq1 = joinpath(phix_dir, "phix174_illumina1.fq")
    short_fastq2 = joinpath(phix_dir, "phix174_illumina2.fq")
    println("    Generated Illumina reads: $short_fastq1, $short_fastq2")

    # Long reads: coverage scaled for CI, PacBio HiFi
    println("  PacBio reads ($(long_cov) coverage)...")
    long_fastq = Mycelia.simulate_pacbio_reads(
        fasta=phix_ref,
        quantity=long_cov,
        outfile=joinpath(phix_dir, "phix174_pacbio.fastq.gz")
    )
    # Uncompress the file
    run(`gunzip -f $long_fastq`)
    long_fastq = replace(long_fastq, ".gz" => "")
    println("    Generated PacBio reads: $long_fastq")

    # Detect platform
    is_macos = Sys.isapple()
    is_linux = Sys.islinux()

    # Define assemblers to compare
    # Note: Some assemblers are platform-specific
    all_short_read_assemblers = [
        (name="SPAdes", func=Mycelia.run_spades, platforms=[:linux, :macos], is_protein=false),
        (name="SKESA", func=Mycelia.run_skesa, platforms=[:linux, :macos], is_protein=false),
        (name="MEGAHIT", func=Mycelia.run_megahit, platforms=[:linux], is_protein=false),  # Linux only
        (name="metaSPAdes", func=Mycelia.run_metaspades, platforms=[:linux, :macos], is_protein=false),
        (name="PenguiN (guided)", func=Mycelia.run_penguin_guided_nuclassemble, platforms=[:linux, :macos], is_protein=false),
        (name="PenguiN", func=Mycelia.run_penguin_nuclassemble, platforms=[:linux, :macos], is_protein=false),
        (name="PLASS", func=Mycelia.run_plass_assemble, platforms=[:linux, :macos], is_protein=true),
    ]

    all_long_read_assemblers = [
        (name="Flye", func=Mycelia.run_flye, platforms=[:linux, :macos], is_protein=false),
        (name="Canu", func=Mycelia.run_canu, platforms=[:linux, :macos], is_protein=false),
        (name="metaFlye", func=Mycelia.run_metaflye, platforms=[:linux, :macos], is_protein=false),
        (name="metaMDBG", func=Mycelia.run_metamdbg, platforms=[:linux, :macos], is_protein=false),
    ]

    # Filter assemblers based on current platform
    current_platform = is_linux ? :linux : (is_macos ? :macos : :unknown)
    short_read_assemblers = filter(a -> current_platform in a.platforms, all_short_read_assemblers)
    long_read_assemblers = filter(a -> current_platform in a.platforms, all_long_read_assemblers)

    # Trim heavy assemblers in CI to reduce resource pressure
    if is_ci
        short_read_assemblers = filter(a -> !(a.name in ["metaSPAdes"]), short_read_assemblers)
        long_read_assemblers = filter(a -> !(a.name in ["metaMDBG"]), long_read_assemblers)
    end

    println("  Platform: $current_platform (CI=$(is_ci))")
    println("  Testing $(length(short_read_assemblers)) short-read assemblers")
    println("  Testing $(length(long_read_assemblers)) long-read assemblers")

    # Storage for comparison results
comparison_results = DataFrames.DataFrame(
    assembler = String[],
    read_type = String[],
    runtime_s = Float64[],
    qc_quast_runtime_s = Float64[],
    qc_quast_outdir = Vector{Union{String,Missing}}(),
    qc_busco_runtime_s = Float64[],
    qc_busco_outdir = Vector{Union{String,Missing}}(),
    n_contigs = Vector{Union{Int,Missing}}(),
    total_length = Vector{Union{Int,Missing}}(),
    n50 = Vector{Union{Int,Missing}}(),
    longest_contig = Vector{Union{Int,Missing}}(),
    identity_vs_ref = Vector{Union{Float64,Missing}}()
)

    # Helper function to calculate assembly metrics
    function calculate_assembly_metrics(contigs_file, reference_file)
        if !isfile(contigs_file)
            return (n_contigs=0, total_length=0, n50=0, longest=0, identity=0.0)
        end

        # Basic metrics
        n_contigs, total_length, n50 = Mycelia.assess_assembly_quality(contigs_file)

        # Find longest contig
        contig_records = collect(FASTX.FASTA.Reader(open(contigs_file)))
        longest = maximum([length(FASTX.sequence(r)) for r in contig_records]; init=0)

        # Calculate ANI (Average Nucleotide Identity) using FastANI
        identity = 0.0
        try
            # Use FastANI for genome comparison
            ani_outfile = tempname() * "_ani.txt"
            Mycelia.fastani_pair(
                query=contigs_file,
                reference=reference_file,
                outfile=ani_outfile
            )

            # Parse FastANI output (format: query, reference, ANI, count, total_frags)
            if isfile(ani_outfile) && filesize(ani_outfile) > 0
                ani_line = readline(ani_outfile)
                fields = split(ani_line, '\t')
                if length(fields) >= 3
                    identity = parse(Float64, fields[3])
                end
            end

            # Cleanup temporary file
            isfile(ani_outfile) && rm(ani_outfile, force=true)
        catch e
            @warn "ANI calculation failed: $e"
        end

        return (n_contigs=n_contigs, total_length=total_length, n50=n50, longest=longest, identity=identity)
    end

    # Run short-read assemblers
    println("\n--- Running Short-Read Assemblers ---")
    for assembler in short_read_assemblers
        println("\nRunning $(assembler.name)...")
        outdir = joinpath(phix_dir, "assembly_$(assembler.name)")

        try
            is_protein = get(assembler, :is_protein, false)
            # Time the assembly
            asm_call = () -> begin
                if is_protein
                    return assembler.func(
                        reads1=short_fastq1,
                        reads2=short_fastq2,
                        outdir=outdir
                    )
                else
                    return assembler.func(
                        fastq1=short_fastq1,
                        fastq2=short_fastq2,
                        outdir=outdir
                    )
                end
            end
            asm_out = _time_and_mem(asm_call)
            result = asm_out.result
            runtime = asm_out.runtime

            # Calculate metrics (skip for protein-only assemblies)
            contigs_file = get(result, :contigs, nothing)
            metrics = if !is_protein && contigs_file !== nothing
                calculate_assembly_metrics(contigs_file, phix_ref)
            else
                (; n_contigs=missing, total_length=missing, n50=missing, longest=missing, identity=missing)
            end

            # QUAST (with reference) and BUSCO (best-effort, optional skip)
            quast_outdir = missing
            quast_runtime = 0.0
            busco_outdir = missing
            busco_runtime = 0.0

            if !is_protein && contigs_file !== nothing
                try
                    quast_run = _time_and_mem(() -> Mycelia.run_quast(contigs_file; reference=phix_ref))
                    quast_outdir = quast_run.result
                    quast_runtime = quast_run.runtime
                catch e
                    @warn "QUAST failed for $(assembler.name): $e"
                end

                if !skip_busco
                    try
                        busco_run = _time_and_mem(() -> Mycelia.run_busco([contigs_file]; mode="genome"))
                        busco_outdir = busco_run.result
                        busco_runtime = busco_run.runtime
                    catch e
                        @warn "BUSCO failed for $(assembler.name): $e"
                    end
                end
            end

            # Store results
            push!(comparison_results, (
                assembler = assembler.name,
                read_type = "short",
                runtime_s = runtime,
                qc_quast_runtime_s = quast_runtime,
                qc_quast_outdir = quast_outdir,
                qc_busco_runtime_s = busco_runtime,
                qc_busco_outdir = busco_outdir,
                n_contigs = metrics.n_contigs,
                total_length = metrics.total_length,
                n50 = metrics.n50,
                longest_contig = metrics.longest,
                identity_vs_ref = metrics.identity
            ))

            println("  Runtime: $(round(runtime, digits=2)) seconds")
            println("  Contigs: $(something(metrics.n_contigs, \"n/a\"))")
            println("  Total length: $(something(metrics.total_length, \"n/a\")) bp")
            println("  N50: $(something(metrics.n50, \"n/a\")) bp")
            println("  Longest contig: $(something(metrics.longest, \"n/a\")) bp")
            identity_str = isnothing(metrics.identity) || metrics.identity === missing ? "n/a" : "$(round(metrics.identity, digits=2))%"
            println("  Identity vs ref: $identity_str")
            if quast_runtime > 0 && quast_outdir !== missing
                println("  QUAST runtime: $(round(quast_runtime, digits=2)) seconds (outdir=$(quast_outdir))")
            end
            if busco_runtime > 0 && busco_outdir !== missing
                println("  BUSCO runtime: $(round(busco_runtime, digits=2)) seconds (outdir=$(busco_outdir))")
            end

        catch e
            @warn "$(assembler.name) failed: $e"
            # Add failed entry
            push!(comparison_results, (
                assembler = assembler.name,
                read_type = "short",
                runtime_s = 0.0,
                qc_quast_runtime_s = 0.0,
                qc_quast_outdir = missing,
                qc_busco_runtime_s = 0.0,
                qc_busco_outdir = missing,
                n_contigs = missing,
                total_length = missing,
                n50 = missing,
                longest_contig = missing,
                identity_vs_ref = missing
            ))
        end
    end

    # Run long-read assemblers
    println("\n--- Running Long-Read Assemblers ---")
    for assembler in long_read_assemblers
        println("\nRunning $(assembler.name)...")
        outdir = joinpath(phix_dir, "assembly_$(assembler.name)")

        try
            is_protein = get(assembler, :is_protein, false)
            # Time the assembly
            asm_call = () -> begin
                if assembler.name == "Canu"
                    return assembler.func(
                        fastq=long_fastq,
                        outdir=outdir,
                        genome_size="5k"  # PhiX174 is ~5.4kb
                    )
                elseif assembler.name == "metaMDBG"
                    # metaMDBG uses hifi_reads parameter instead of fastq
                    return assembler.func(
                        hifi_reads=long_fastq,
                        outdir=outdir
                    )
                else
                    return assembler.func(
                        fastq=long_fastq,
                        outdir=outdir
                    )
                end
            end
            asm_out = _time_and_mem(asm_call)
            result = asm_out.result
            runtime = asm_out.runtime

            # Calculate metrics
            contigs_file = result.contigs
            metrics = if !is_protein
                calculate_assembly_metrics(contigs_file, phix_ref)
            else
                (; n_contigs=missing, total_length=missing, n50=missing, longest=missing, identity=missing)
            end

            # QUAST (with reference) and BUSCO (best-effort)
            quast_outdir = missing
            quast_runtime = 0.0
            busco_outdir = missing
            busco_runtime = 0.0

            if !is_protein
                try
                    quast_run = _time_and_mem(() -> Mycelia.run_quast(contigs_file; reference=phix_ref))
                    quast_outdir = quast_run.result
                    quast_runtime = quast_run.runtime
                catch e
                    @warn "QUAST failed for $(assembler.name): $e"
                end

                if !skip_busco
                    try
                        busco_run = _time_and_mem(() -> Mycelia.run_busco([contigs_file]; mode="genome"))
                        busco_outdir = busco_run.result
                        busco_runtime = busco_run.runtime
                    catch e
                        @warn "BUSCO failed for $(assembler.name): $e"
                    end
                end
            end

            # Store results
            push!(comparison_results, (
                assembler = assembler.name,
                read_type = "long",
                runtime_s = runtime,
                qc_quast_runtime_s = quast_runtime,
                qc_quast_outdir = quast_outdir,
                qc_busco_runtime_s = busco_runtime,
                qc_busco_outdir = busco_outdir,
                n_contigs = metrics.n_contigs,
                total_length = metrics.total_length,
                n50 = metrics.n50,
                longest_contig = metrics.longest,
                identity_vs_ref = metrics.identity
            ))

            println("  Runtime: $(round(runtime, digits=2)) seconds")
            println("  Contigs: $(something(metrics.n_contigs, \"n/a\"))")
            println("  Total length: $(something(metrics.total_length, \"n/a\")) bp")
            println("  N50: $(something(metrics.n50, \"n/a\")) bp")
            println("  Longest contig: $(something(metrics.longest, \"n/a\")) bp")
            identity_str = isnothing(metrics.identity) || metrics.identity === missing ? "n/a" : "$(round(metrics.identity, digits=2))%"
            println("  Identity vs ref: $identity_str")
            if quast_runtime > 0 && quast_outdir !== missing
                println("  QUAST runtime: $(round(quast_runtime, digits=2)) seconds (outdir=$(quast_outdir))")
            end
            if busco_runtime > 0 && busco_outdir !== missing
                println("  BUSCO runtime: $(round(busco_runtime, digits=2)) seconds (outdir=$(busco_outdir))")
            end

        catch e
            @warn "$(assembler.name) failed: $e"
            # Add failed entry
            push!(comparison_results, (
                assembler = assembler.name,
                read_type = "long",
                runtime_s = 0.0,
                qc_quast_runtime_s = 0.0,
                qc_quast_outdir = missing,
                qc_busco_runtime_s = 0.0,
                qc_busco_outdir = missing,
                n_contigs = missing,
                total_length = missing,
                n50 = missing,
                longest_contig = missing,
                identity_vs_ref = missing
            ))
        end
    end

    # Generate comparison plots
    println("\n--- Generating Comparison Plots ---")

    # Filter successful assemblies
    successful_results = filter(row -> row.n_contigs > 0, comparison_results)

    if !isempty(successful_results)
        # N50 comparison
        n50_plot = Plots.bar(
            successful_results.assembler,
            successful_results.n50,
            title="N50 Comparison (PhiX174)",
            ylabel="N50 (bp)",
            xlabel="Assembler",
            legend=false,
            xrotation=45,
            color=[r.read_type == "short" ? :blue : :green for r in eachrow(successful_results)]
        )
        n50_file = joinpath(results_dir, "phix174_n50_comparison.png")
        Plots.savefig(n50_plot, n50_file)
        println("  Saved N50 plot: $n50_file")

        # Runtime comparison
        runtime_plot = Plots.bar(
            successful_results.assembler,
            successful_results.runtime_s,
            title="Runtime Comparison (PhiX174)",
            ylabel="Runtime (seconds)",
            xlabel="Assembler",
            legend=false,
            xrotation=45,
            color=[r.read_type == "short" ? :blue : :green for r in eachrow(successful_results)]
        )
        runtime_file = joinpath(results_dir, "phix174_runtime_comparison.png")
        Plots.savefig(runtime_plot, runtime_file)
        println("  Saved runtime plot: $runtime_file")

        # Identity comparison
        identity_plot = Plots.bar(
            successful_results.assembler,
            successful_results.identity_vs_ref,
            title="Identity vs Reference (PhiX174)",
            ylabel="Identity (%)",
            xlabel="Assembler",
            legend=false,
            xrotation=45,
            ylim=(90, 100),
            color=[r.read_type == "short" ? :blue : :green for r in eachrow(successful_results)]
        )
        identity_file = joinpath(results_dir, "phix174_identity_comparison.png")
        Plots.savefig(identity_plot, identity_file)
        println("  Saved identity plot: $identity_file")
    end

    # Save results to CSV
    csv_file = joinpath(results_dir, "phix174_comparison.csv")
    CSV.write(csv_file, comparison_results)
    println("  Saved results CSV: $csv_file")

    # Display summary
    println("\n=== PhiX174 Benchmark Summary ===")
    println("Reference: NC_001422.1 ($ref_length bp)")
    println("Assemblers tested: $(DataFrames.nrow(comparison_results))")
    println("Successful assemblies: $(count(row -> row.n_contigs > 0, eachrow(comparison_results)))")
    println("\nResults:")
    for row in eachrow(comparison_results)
        if row.n_contigs > 0
            println("  $(row.assembler) ($(row.read_type)):")
            println("    Runtime: $(round(row.runtime_s, digits=1))s")
            println("    Contigs: $(row.n_contigs)")
            println("    N50: $(row.n50) bp")
            println("    Identity: $(round(row.identity_vs_ref, digits=2))%")
        else
            println("  $(row.assembler) ($(row.read_type)): FAILED")
        end
    end

finally
    # Cleanup
    println("\n--- Cleaning up PhiX174 benchmark data ---")
    rm(phix_dir, force=true, recursive=true)
    println("  Removed: $phix_dir")
end

println("\n=== PhiX174 Assembler Comparison Complete ===")
println("End time: $(Dates.now())")
println("\nResults saved to:")
println("  - $results_dir/phix174_comparison.csv")
println("  - $results_dir/phix174_*.png")

nothing
