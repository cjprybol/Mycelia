# # Genome Assembly Benchmark
#
# This benchmark evaluates genome assembly performance, quality, and
# resource usage across different datasets and assembly parameters.
#
# ## Benchmark Overview
#
# This benchmark tests:
# - Assembly performance with different genome sizes
# - Quality metrics across different assemblers
# - Resource usage optimization
# - Parameter sensitivity analysis
# - Scalability with coverage depth

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import FASTX
import BioSequences
import BenchmarkTools
import Random
import Statistics
import Dates
import DataFrames
import CSV
import Plots

println("=== Genome Assembly Benchmark ===")
println("Start time: $(Dates.now())")

# ## Benchmark Configuration
#
# Configure assembly benchmark parameters

# Small scale configuration
small_config = Dict(
    "genome_sizes" => [50000, 100000],
    "coverage_depths" => [20, 30],
    "assemblers" => ["hifiasm"],
    "n_replicates" => 3,
    "description" => "Small scale - quick validation"
)

# Medium scale configuration
medium_config = Dict(
    "genome_sizes" => [500000, 1000000, 2000000],
    "coverage_depths" => [20, 30, 50],
    "assemblers" => [
        # Short read assemblers (genomic variants)
        "spades", "megahit", "idba_ud", "velvet",
        # Short read assemblers (metagenomic variants)
        "metaspades", "megahit", "metavelvet", "metaidba", "faucet", "ray_meta",
        # Long read assemblers (genomic)
        "flye", "canu", "hifiasm", "nextdenovo", "shasta", "raven", "miniasm", "wtdbg2",
        # Long read assemblers (metagenomic variants)
        "metaflye", "metacanu", "metamdbg",
        # PacBio HiFi specialized
        "hicanu", "lja", "mbg",
        # Hybrid assemblers (genomic)
        "unicycler", "hybridspades", "haslr", "wengan", "masurca", "platanus_allee",
        # Hybrid assemblers (metagenomic)
        "metaplatanus", "dbg2olc",
        # Scaffolding specialists
        "opera_ms", "opera_lg",
        # Mycelia methods
        "mycelia_intelligent", "mycelia_iterative"
    ],
    "n_replicates" => 5,
    "description" => "Medium scale - comprehensive assembler comparison including metagenomic variants (2024-2025)"
)

# Large scale configuration  
large_config = Dict(
    "genome_sizes" => [5000000, 10000000, 50000000],
    "coverage_depths" => [30, 50, 100],
    "assemblers" => [
        # Short read assemblers (genomic variants)
        "spades", "megahit", "idba_ud", "velvet",
        # Short read assemblers (metagenomic variants)
        "metaspades", "megahit", "metavelvet", "metaidba", "faucet", "ray_meta",
        # Long read assemblers (genomic)
        "flye", "canu", "hifiasm", "nextdenovo", "shasta", "raven", "miniasm", "wtdbg2",
        # Long read assemblers (metagenomic variants)
        "metaflye", "metacanu", "metamdbg",
        # PacBio HiFi optimized
        "hicanu", "lja", "mbg", 
        # Hybrid assemblers (genomic)
        "unicycler", "hybridspades", "haslr", "wengan", "masurca", "platanus_allee",
        # Hybrid assemblers (metagenomic)
        "metaplatanus", "dbg2olc",
        # Scaffolding and metagenomic specialists
        "opera_ms", "opera_lg",
        # Mycelia methods
        "mycelia_intelligent", "mycelia_iterative"
    ],
    "n_replicates" => 10,
    "description" => "Large scale - comprehensive benchmarking against all major assemblers including metagenomic variants (2024-2025 state-of-the-art)"
)

# Select configuration
config = if haskey(ENV, "BENCHMARK_SCALE")
    if ENV["BENCHMARK_SCALE"] == "small"
        small_config
    elseif ENV["BENCHMARK_SCALE"] == "medium"
        medium_config
    else
        large_config
    end
else
    small_config
end

println("Benchmark Configuration: $(config["description"])")
println("Genome sizes: $(config["genome_sizes"])")
println("Coverage depths: $(config["coverage_depths"])")
println("Assemblers: $(config["assemblers"])")

# ## Test Data Generation
#
# Generate realistic test datasets for assembly

println("\n--- Generating Test Data ---")

# Include benchmark utilities
include("benchmark_utils.jl")

# Use existing Mycelia simulation functions for paired-end read generation

# Generate test datasets
test_data_dir = "assembly_test_data"
mkpath(test_data_dir)

test_datasets = []

for (idx, genome_size) in enumerate(config["genome_sizes"][1:min(2, length(config["genome_sizes"]))])
    for coverage in config["coverage_depths"][1:min(2, length(config["coverage_depths"]))]
        # Generate reference genome
        reference_genome = BioSequences.randdnaseq(genome_size)
        
        # Generate paired-end reads
        read_length = 150  # Standard Illumina read length
        insert_size = 300  # Standard insert size
        
        reads_1, reads_2 = Mycelia.generate_paired_end_reads(
            reference_genome, coverage, read_length, insert_size
        )
        
        # Save as FASTQ files
        dataset_name = "genome$(idx)_size$(genome_size)_cov$(coverage)"
        fastq1_path = joinpath(test_data_dir, "$(dataset_name)_R1.fastq")
        fastq2_path = joinpath(test_data_dir, "$(dataset_name)_R2.fastq")
        
        # Convert BioSequences to FASTQ records with quality scores
        fastq_records_1 = [FASTX.FASTQ.Record("read_$i", reads_1[i], String([Char(30 + 33) for _ in 1:length(reads_1[i])])) for i in 1:length(reads_1)]
        fastq_records_2 = [FASTX.FASTQ.Record("read_$i", reads_2[i], String([Char(30 + 33) for _ in 1:length(reads_2[i])])) for i in 1:length(reads_2)]
        Mycelia.write_fastq(records=fastq_records_1, filename=fastq1_path)
        Mycelia.write_fastq(records=fastq_records_2, filename=fastq2_path)
        
        push!(test_datasets, Dict(
            "name" => dataset_name,
            "genome_size" => genome_size,
            "coverage" => coverage,
            "fastq1" => fastq1_path,
            "fastq2" => fastq2_path,
            "reference" => reference_genome,
            "n_reads" => length(reads_1)
        ))
    end
end

println("Generated $(length(test_datasets)) test datasets")

# ## Assembly Performance Benchmarks
#
# Test assembly time and resource usage

println("\n--- Assembly Performance Benchmarks ---")

# Initialize benchmark suite
benchmark_suite = BenchmarkSuite("Assembly Performance")

# Test MEGAHIT performance on different dataset sizes
for dataset in test_datasets[1:min(3, length(test_datasets))]  # Limit for benchmarking
    println("\nBenchmarking MEGAHIT assembly: $(dataset["name"])")
    
    # Create output directory for this assembly
    assembly_outdir = joinpath(test_data_dir, "megahit_$(dataset["name"])")
    
    # Benchmark MEGAHIT assembly
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.run_megahit;
        fastq1=dataset["fastq1"],
        fastq2=dataset["fastq2"],
        outdir=assembly_outdir,
        min_contig_len=200,
        samples=1,  # Assembly is slow, so use fewer samples
        seconds=600  # Allow up to 10 minutes per assembly
    )
    
    add_benchmark_result!(
        benchmark_suite,
        "megahit_$(dataset["name"])",
        benchmark_result,
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e9, digits=2)) seconds")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
    
    # Comprehensive assembly quality assessment using existing Mycelia functions
    contigs_file = joinpath(assembly_outdir, "final.contigs.fa")
    if isfile(contigs_file)
        # Primary accuracy metrics using existing k-mer QV function
        reads_files = [dataset["fastq1"], dataset["fastq2"]]
        qv_results = try
            Mycelia.assess_assembly_kmer_quality(
                assembly=contigs_file, 
                observations=reads_files, 
                ks=[21]  # Use k=21 for benchmarking
            )
        catch e
            @warn "K-mer QV assessment failed: $e"
            DataFrames.DataFrame()
        end
        
        # Read mapping quality using existing minimap and qualimap functions
        mapping_results = nothing
        try
            # Map reads to assembly using minimap_map (no index required)
            map_result = Mycelia.minimap_map(
                fasta=contigs_file,
                fastq=dataset["fastq1"],  # Map first read file
                mapping_type="sr",
                threads=4
            )
            
            # Run the mapping command
            run(map_result.cmd)
            
            # Convert compressed SAM to BAM for qualimap
            sam_file = map_result.outfile  # This will be .sam.gz
            if isfile(sam_file)
                # Decompress and convert to BAM
                bam_file = replace(sam_file, ".sam.gz" => ".bam")
                run(`$(Mycelia.CONDA_RUNNER) run -n samtools gunzip -c $(sam_file) | samtools view -bS - > $(bam_file)`)
                
                # Run qualimap on the BAM file
                if isfile(bam_file)
                    Mycelia.add_bioconda_env("qualimap")
                    qualimap_outdir = joinpath(assembly_outdir, "qualimap")
                    run(`$(Mycelia.CONDA_RUNNER) run -n qualimap qualimap bamqc -bam $(bam_file) -outdir $(qualimap_outdir)`)
                    
                    # Parse qualimap results using existing function
                    qualimap_txt = joinpath(qualimap_outdir, "genome_results.txt")
                    if isfile(qualimap_txt)
                        mapping_results = Mycelia.parse_qualimap_contig_coverage(qualimap_txt)
                    end
                end
            end
        catch e
            @warn "Read mapping assessment failed: $e"
        end
        
        # Secondary contiguity metrics
        n_contigs, total_length, n50 = Mycelia.assess_assembly_quality(contigs_file)
        
        println("  Assembly quality:")
        println("    Primary metrics (accuracy-based):")
        if !isempty(qv_results)
            qv_score = qv_results[1, :qv]
            println("      QV score: $(round(qv_score, digits=2))")
            println("      Cosine distance: $(round(qv_results[1, :cosine_distance], digits=4))")
            println("      JS divergence: $(round(qv_results[1, :js_divergence], digits=4))")
        end
        if mapping_results !== nothing
            total_mapped_bases = sum(mapping_results[!, "Mapped bases"])
            println("      Total mapped bases: $total_mapped_bases")
            println("      Mean coverage: $(round(Statistics.mean(mapping_results[!, "Mean coverage"]), digits=2))")
        end
        println("    Secondary metrics (contiguity):")
        println("      Contigs: $n_contigs")
        println("      Total length: $total_length bp")
        println("      N50: $n50 bp")
        
        # Add comprehensive quality metrics to results
        benchmark_suite.results["megahit_$(dataset["name"])"]["assembly_quality"] = Dict(
            "primary_metrics" => Dict(
                "qv_score" => !isempty(qv_results) ? qv_results[1, :qv] : missing,
                "cosine_distance" => !isempty(qv_results) ? qv_results[1, :cosine_distance] : missing,
                "js_divergence" => !isempty(qv_results) ? qv_results[1, :js_divergence] : missing,
                "total_mapped_bases" => mapping_results !== nothing ? sum(mapping_results[!, "Mapped bases"]) : missing,
                "mean_coverage" => mapping_results !== nothing ? Statistics.mean(mapping_results[!, "Mean coverage"]) : missing
            ),
            "secondary_metrics" => Dict(
                "n_contigs" => n_contigs,
                "total_length" => total_length,
                "n50" => n50,
                "expected_length" => dataset["genome_size"]
            )
        )
    end
end

# Test metaSPAdes performance (if smaller datasets)
for dataset in test_datasets[1:min(1, length(test_datasets))]  # Even more limited for metaSPAdes
    println("\nBenchmarking metaSPAdes assembly: $(dataset["name"])")
    
    assembly_outdir = joinpath(test_data_dir, "metaspades_$(dataset["name"])")
    
    try
        benchmark_result, memory_stats = run_benchmark_with_memory(
            Mycelia.run_metaspades;
            fastq1=dataset["fastq1"],
            fastq2=dataset["fastq2"],
            outdir=assembly_outdir,
            samples=1,
            seconds=900  # Allow up to 15 minutes
        )
        
        add_benchmark_result!(
            benchmark_suite,
            "metaspades_$(dataset["name"])",
            benchmark_result,
            memory_stats
        )
        
        println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e9, digits=2)) seconds")
        println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
        
    catch e
        println("  ⚠️  metaSPAdes benchmark failed: $e")
        # Continue with other benchmarks
    end
end

# ## Assembly Quality Benchmarks
#
# Test assembly quality metrics

println("\n--- Assembly Quality Benchmarks ---")

# Assembly quality assessment using existing Mycelia function

# Calculate comprehensive assembly quality metrics for completed assemblies
assembly_quality_summary = Dict()

for (test_name, result) in benchmark_suite.results
    if haskey(result, "assembly_quality")
        quality = result["assembly_quality"]
        
        # Extract primary and secondary metrics
        primary_metrics = get(quality, "primary_metrics", Dict())
        secondary_metrics = get(quality, "secondary_metrics", Dict())
        
        # Calculate derived scores
        length_recovery = secondary_metrics["total_length"] / secondary_metrics["expected_length"]
        contiguity_score = secondary_metrics["n50"] / secondary_metrics["total_length"]  # Normalized N50
        
        assembly_quality_summary[test_name] = Dict(
            # Primary metrics (accuracy-based)
            "qv_score" => get(primary_metrics, "qv_score", missing),
            "cosine_distance" => get(primary_metrics, "cosine_distance", missing),
            "js_divergence" => get(primary_metrics, "js_divergence", missing),
            "total_mapped_bases" => get(primary_metrics, "total_mapped_bases", missing),
            "mean_coverage" => get(primary_metrics, "mean_coverage", missing),
            # Secondary metrics (contiguity-based)
            "length_recovery" => round(length_recovery, digits=3),
            "contiguity_score" => round(contiguity_score, digits=6),
            "n_contigs" => secondary_metrics["n_contigs"],
            "n50" => secondary_metrics["n50"]
        )
        
        println("Assembly quality for $test_name:")
        println("  Primary metrics (accuracy-based):")
        qv_score = get(primary_metrics, "qv_score", missing)
        if !ismissing(qv_score)
            println("    QV score: $(round(qv_score, digits=2))")
            println("    Cosine distance: $(round(get(primary_metrics, "cosine_distance", 0.0), digits=4))")
            println("    JS divergence: $(round(get(primary_metrics, "js_divergence", 0.0), digits=4))")
        end
        total_mapped = get(primary_metrics, "total_mapped_bases", missing)
        if !ismissing(total_mapped)
            println("    Total mapped bases: $total_mapped")
            println("    Mean coverage: $(round(get(primary_metrics, "mean_coverage", 0.0), digits=2))")
        end
        println("  Secondary metrics (contiguity):")
        println("    Length recovery: $(round(length_recovery * 100, digits=1))%")
        println("    Contiguity score: $(round(contiguity_score, digits=6))")
        println("    Contigs: $(secondary_metrics["n_contigs"])")
        println("    N50: $(secondary_metrics["n50"]) bp")
    end
end

# ## Assembler Comparison Benchmarks
#
# Compare different assembly tools

println("\n--- Assembler Comparison Benchmarks ---")

# PhiX174 Reference Genome Benchmarking
# Head-to-head comparison of 8 assemblers on standard reference
println("\n=== PhiX174 Assembler Comparison ===")

# Create temporary directory for PhiX174 data
phix_dir = mktempdir(prefix="phix174_benchmark_")
println("PhiX174 benchmark directory: $phix_dir")

try
    # Download PhiX174 reference
    println("\nDownloading PhiX174 reference...")
    phix_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
    phix_gz = joinpath(phix_dir, "phix174.fna.gz")
    phix_ref = joinpath(phix_dir, "phix174_reference.fasta")

    download(phix_url, phix_gz)
    run(pipeline(`gunzip -c $phix_gz`, phix_ref))
    println("  Reference downloaded: $phix_ref")

    # Load reference sequence
    ref_records = collect(FASTX.FASTA.Reader(open(phix_ref)))
    ref_seq = FASTX.sequence(BioSequences.LongDNA{4}, ref_records[1])
    ref_length = length(ref_seq)
    println("  Reference length: $ref_length bp")

    # Generate simulated reads
    println("\nGenerating simulated reads...")

    # Short reads: 100x coverage, 150bp PE
    println("  Illumina reads (100x coverage)...")
    short_r1, short_r2 = Mycelia.simulate_illumina_reads(
        reference=ref_seq,
        coverage=100,
        read_length=150,
        insert_size=300,
        profile="HS25",
        seed=12345
    )
    short_fastq1 = joinpath(phix_dir, "phix174_illumina_R1.fastq")
    short_fastq2 = joinpath(phix_dir, "phix174_illumina_R2.fastq")
    Mycelia.save_reads_as_fastq(short_r1, short_fastq1)
    Mycelia.save_reads_as_fastq(short_r2, short_fastq2)
    println("    Generated $(length(short_r1)) paired-end reads")

    # Long reads: 50x coverage, PacBio HiFi
    println("  PacBio reads (50x coverage)...")
    long_reads = Mycelia.simulate_pacbio_reads(
        reference=ref_seq,
        coverage=50,
        mean_length=10000,
        length_sd=2000,
        identity=99.9,
        seed=12345
    )
    long_fastq = joinpath(phix_dir, "phix174_pacbio.fastq")
    Mycelia.save_reads_as_fastq(long_reads, long_fastq)
    println("    Generated $(length(long_reads)) long reads")

    # Define assemblers to compare
    short_read_assemblers = [
        (name="SPAdes", func=Mycelia.run_spades),
        (name="SKESA", func=Mycelia.run_skesa),
        (name="MEGAHIT", func=Mycelia.run_megahit),
        (name="metaSPAdes", func=Mycelia.run_metaspades),
    ]

    long_read_assemblers = [
        (name="Flye", func=Mycelia.run_flye),
        (name="Canu", func=Mycelia.run_canu),
        (name="metaFlye", func=Mycelia.run_metaflye),
        (name="metaMDBG", func=Mycelia.run_metamdbg),
    ]

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

        # Align to reference with minimap2 to get identity
        identity = 0.0
        try
            # Create minimap2 index
            idx_result = Mycelia.minimap_index(fasta=reference_file)
            run(idx_result.cmd)

            # Map contigs to reference
            map_result = Mycelia.minimap_map_with_index(
                index=idx_result.outfile,
                fasta=contigs_file,
                mapping_type="asm5"  # Assembly to reference mapping
            )
            run(map_result.cmd)

            # Parse PAF output to get identity
            if isfile(map_result.outfile)
                paf_lines = readlines(map_result.outfile)
                if !isempty(paf_lines)
                    # PAF format: col 10 is number of matches, col 11 is alignment block length
                    identities = []
                    for line in paf_lines
                        fields = split(line, '\t')
                        if length(fields) >= 11
                            matches = parse(Int, fields[10])
                            block_len = parse(Int, fields[11])
                            push!(identities, 100.0 * matches / block_len)
                        end
                    end
                    if !isempty(identities)
                        identity = Statistics.mean(identities)
                    end
                end
            end
        catch e
            @warn "Identity calculation failed: $e"
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
    println("Assemblers tested: $(nrow(comparison_results))")
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

println("\n=== PhiX174 Assembler Comparison Complete ===\n")

# ## Parameter Optimization Benchmarks
#
# Test parameter sensitivity and optimization

println("\n--- Parameter Optimization Benchmarks ---")

# TODO: Implement parameter optimization benchmarks
# - K-mer size optimization
# - Coverage threshold optimization
# - Quality score threshold optimization
# - Multi-objective optimization

# ## Scalability Benchmarks
#
# Test scalability with increasing data sizes

println("\n--- Scalability Benchmarks ---")

# TODO: Implement scalability benchmarks
# - Performance scaling with genome size
# - Memory scaling with coverage
# - Parallel processing efficiency
# - Resource utilization optimization

# ## Error Handling Benchmarks
#
# Test robustness to data quality issues

println("\n--- Error Handling Benchmarks ---")

# TODO: Implement error handling benchmarks
# - Assembly with contaminated data
# - Assembly with low-quality reads
# - Assembly with uneven coverage
# - Assembly with high error rates

# ## Validation Benchmarks
#
# Test assembly validation approaches

println("\n--- Validation Benchmarks ---")

# TODO: Implement validation benchmarks
# - Reference-based validation performance
# - K-mer based validation (Merqury)
# - BUSCO validation performance
# - Read mapping validation

# ## Results Collection
#
# Collect comprehensive benchmark results

println("\n--- Collecting Results ---")

# Calculate assembly performance metrics
total_datasets_processed = length(test_datasets)
total_reads_processed = sum([dataset["n_reads"] for dataset in test_datasets])

# Assembly throughput estimates
assembly_throughput = Dict()
for (test_name, result) in benchmark_suite.results
    if startswith(test_name, "megahit_") && haskey(result, "median_time")
        time_minutes = result["median_time"] / 1e9 / 60
        
        # Find corresponding dataset
        dataset_name = replace(test_name, "megahit_" => "")
        matching_dataset = findfirst(d -> d["name"] == dataset_name, test_datasets)
        
        if matching_dataset !== nothing
            dataset = test_datasets[matching_dataset]
            reads_per_minute = (dataset["n_reads"] * 2) / time_minutes  # Paired reads
            assembly_throughput[test_name] = round(reads_per_minute, digits=0)
        end
    end
end

# Create comprehensive results structure
results_dir = "results"
mkpath(results_dir)

comprehensive_results = Dict(
    "benchmark_name" => "assembly",
    "timestamp" => Dates.now(),
    "configuration" => config,
    "system_info" => Dict(
        "julia_version" => VERSION,
        "cpu_threads" => Sys.CPU_THREADS,
        "total_memory" => Sys.total_memory(),
        "hostname" => gethostname(),
        "platform" => Sys.MACHINE
    ),
    "dataset_summary" => Dict(
        "total_datasets_processed" => total_datasets_processed,
        "total_reads_processed" => total_reads_processed,
        "genome_sizes_tested" => config["genome_sizes"],
        "coverage_depths_tested" => config["coverage_depths"]
    ),
    "performance_metrics" => Dict(
        "assembly_throughput_reads_per_minute" => assembly_throughput,
        "benchmark_duration_minutes" => round((Dates.now() - benchmark_suite.metadata["timestamp"]).value / 60000, digits=2)
    ),
    "quality_metrics" => assembly_quality_summary,
    "resource_utilization" => Dict(
        "peak_memory_mb" => round(maximum([get(result, "memory", 0) for result in values(benchmark_suite.results)]) / 1e6, digits=2),
        "total_allocations" => sum([get(result, "allocations", 0) for result in values(benchmark_suite.results)])
    ),
    "detailed_results" => benchmark_suite.results
)

# Save results
results_file = joinpath(results_dir, "assembly_benchmark_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")
save_benchmark_results(benchmark_suite, results_file)

# Save comprehensive results
comprehensive_file = joinpath(results_dir, "assembly_comprehensive_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")
open(comprehensive_file, "w") do f
    JSON.print(f, comprehensive_results, 2)
end

# Display summary
format_benchmark_summary(benchmark_suite)

println("\n--- Assembly Performance Summary ---")
println("Datasets processed: $total_datasets_processed")
println("Total reads processed: $total_reads_processed")
println("Assembly throughput:")
for (test_name, throughput) in assembly_throughput
    println("  $test_name: $(Int(throughput)) reads/minute")
end

# Cleanup test data
println("\n--- Cleanup ---")
for dataset in test_datasets
    rm(dataset["fastq1"], force=true)
    rm(dataset["fastq2"], force=true)
end
rm(test_data_dir, force=true, recursive=true)

println("\n=== Assembly Benchmark Complete ===")
println("Results saved to: $results_file")
println("Comprehensive results saved to: $comprehensive_file")
println("End time: $(Dates.now())")

nothing