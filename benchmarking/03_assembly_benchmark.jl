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
import BenchmarkTools
import Random
import Statistics
import Dates

println("=== Genome Assembly Benchmark ===")
println("Start time: $(now())")

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
        
        Mycelia.save_reads_as_fastq(reads_1, fastq1_path)
        Mycelia.save_reads_as_fastq(reads_2, fastq2_path)
        
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

# TODO: Implement assembler comparison benchmarks
# - Performance comparison (runtime, memory usage)
# - Quality comparison (N50, completeness, accuracy)
# - Parameter sensitivity analysis
# - Robustness testing
# - Mycelia vs external assemblers:
#   - Short reads: MEGAHIT vs metaSPAdes vs Mycelia
#   - Long reads: Flye vs Canu vs hifiasm vs Mycelia
#   - Hybrid: Unicycler vs Mycelia hybrid approach
# - Error correction effectiveness:
#   - Viterbi polishing vs traditional polishing
#   - Iterative polishing convergence
#   - Probabilistic vs deterministic approaches

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