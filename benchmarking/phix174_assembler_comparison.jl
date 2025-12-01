# PhiX174 Assembler Comparison Benchmark
# Head-to-head comparison of 8 third-party assemblers on PhiX174 reference genome

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

println("=== PhiX174 Assembler Comparison ===")
println("Start time: $(Dates.now())")

# Create temporary directory for PhiX174 data
phix_dir = mktempdir(prefix="phix174_benchmark_")
println("PhiX174 benchmark directory: $phix_dir")

# Create results directory
results_dir = "results"
mkpath(results_dir)

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

    # Generate simulated reads
    println("\nGenerating simulated reads...")

    # Short reads: 100x coverage, 150bp PE
    println("  Illumina reads (100x coverage)...")
    Mycelia.simulate_illumina_reads(
        fasta=phix_ref,
        coverage=100,
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

    # Long reads: 50x coverage, PacBio HiFi
    println("  PacBio reads (50x coverage)...")
    long_fastq = Mycelia.simulate_pacbio_reads(
        fasta=phix_ref,
        quantity="50x",
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
        (name="SPAdes", func=Mycelia.run_spades, platforms=[:linux, :macos]),
        (name="SKESA", func=Mycelia.run_skesa, platforms=[:linux, :macos]),
        (name="MEGAHIT", func=Mycelia.run_megahit, platforms=[:linux]),  # Linux only
        (name="metaSPAdes", func=Mycelia.run_metaspades, platforms=[:linux, :macos]),
    ]

    all_long_read_assemblers = [
        (name="Flye", func=Mycelia.run_flye, platforms=[:linux, :macos]),
        (name="Canu", func=Mycelia.run_canu, platforms=[:linux, :macos]),
        (name="metaFlye", func=Mycelia.run_metaflye, platforms=[:linux, :macos]),
        (name="metaMDBG", func=Mycelia.run_metamdbg, platforms=[:linux, :macos]),
    ]

    # Filter assemblers based on current platform
    current_platform = is_linux ? :linux : (is_macos ? :macos : :unknown)
    short_read_assemblers = filter(a -> current_platform in a.platforms, all_short_read_assemblers)
    long_read_assemblers = filter(a -> current_platform in a.platforms, all_long_read_assemblers)

    println("  Platform: $current_platform")
    println("  Testing $(length(short_read_assemblers)) short-read assemblers")
    println("  Testing $(length(long_read_assemblers)) long-read assemblers")

    # Storage for comparison results
    comparison_results = DataFrames.DataFrame(
        assembler = String[],
        read_type = String[],
        runtime_s = Float64[],
        n_contigs = Int[],
        total_length = Int[],
        n50 = Int[],
        longest_contig = Int[],
        identity_vs_ref = Float64[]
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
            # Time the assembly
            start_time = time()
            result = assembler.func(
                fastq1=short_fastq1,
                fastq2=short_fastq2,
                outdir=outdir
            )
            runtime = time() - start_time

            # Calculate metrics
            contigs_file = result.contigs
            metrics = calculate_assembly_metrics(contigs_file, phix_ref)

            # Store results
            push!(comparison_results, (
                assembler = assembler.name,
                read_type = "short",
                runtime_s = runtime,
                n_contigs = metrics.n_contigs,
                total_length = metrics.total_length,
                n50 = metrics.n50,
                longest_contig = metrics.longest,
                identity_vs_ref = metrics.identity
            ))

            println("  Runtime: $(round(runtime, digits=2)) seconds")
            println("  Contigs: $(metrics.n_contigs)")
            println("  Total length: $(metrics.total_length) bp")
            println("  N50: $(metrics.n50) bp")
            println("  Longest contig: $(metrics.longest) bp")
            println("  Identity vs ref: $(round(metrics.identity, digits=2))%")

        catch e
            @warn "$(assembler.name) failed: $e"
            # Add failed entry
            push!(comparison_results, (
                assembler = assembler.name,
                read_type = "short",
                runtime_s = 0.0,
                n_contigs = 0,
                total_length = 0,
                n50 = 0,
                longest_contig = 0,
                identity_vs_ref = 0.0
            ))
        end
    end

    # Run long-read assemblers
    println("\n--- Running Long-Read Assemblers ---")
    for assembler in long_read_assemblers
        println("\nRunning $(assembler.name)...")
        outdir = joinpath(phix_dir, "assembly_$(assembler.name)")

        try
            # Time the assembly
            start_time = time()
            result = if assembler.name == "Canu"
                assembler.func(
                    fastq=long_fastq,
                    outdir=outdir,
                    genome_size="5k"  # PhiX174 is ~5.4kb
                )
            elseif assembler.name == "metaMDBG"
                # metaMDBG uses hifi_reads parameter instead of fastq
                assembler.func(
                    hifi_reads=long_fastq,
                    outdir=outdir
                )
            else
                assembler.func(
                    fastq=long_fastq,
                    outdir=outdir
                )
            end
            runtime = time() - start_time

            # Calculate metrics
            contigs_file = result.contigs
            metrics = calculate_assembly_metrics(contigs_file, phix_ref)

            # Store results
            push!(comparison_results, (
                assembler = assembler.name,
                read_type = "long",
                runtime_s = runtime,
                n_contigs = metrics.n_contigs,
                total_length = metrics.total_length,
                n50 = metrics.n50,
                longest_contig = metrics.longest,
                identity_vs_ref = metrics.identity
            ))

            println("  Runtime: $(round(runtime, digits=2)) seconds")
            println("  Contigs: $(metrics.n_contigs)")
            println("  Total length: $(metrics.total_length) bp")
            println("  N50: $(metrics.n50) bp")
            println("  Longest contig: $(metrics.longest) bp")
            println("  Identity vs ref: $(round(metrics.identity, digits=2))%")

        catch e
            @warn "$(assembler.name) failed: $e"
            # Add failed entry
            push!(comparison_results, (
                assembler = assembler.name,
                read_type = "long",
                runtime_s = 0.0,
                n_contigs = 0,
                total_length = 0,
                n50 = 0,
                longest_contig = 0,
                identity_vs_ref = 0.0
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
