# # Gene Annotation Benchmark
#
# This benchmark evaluates gene annotation performance, accuracy, and
# computational efficiency across different genome types and annotation methods.
#
# ## Benchmark Overview
#
# This benchmark tests:
# - Gene prediction accuracy across different organisms
# - Functional annotation performance and quality
# - Comparative annotation efficiency
# - Annotation pipeline scalability
# - Database search performance

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

println("=== Gene Annotation Benchmark ===")
println("Start time: $(now())")

# ## Benchmark Configuration
#
# Configure annotation benchmark parameters

# Small scale configuration
small_config = Dict(
    "genome_types" => ["bacterial", "viral"],
    "genome_sizes" => [1000000, 5000000],
    "annotation_methods" => ["ab_initio", "homology"],
    "n_genomes" => 5,
    "description" => "Small scale - quick validation"
)

# Medium scale configuration
medium_config = Dict(
    "genome_types" => ["bacterial", "viral", "fungal"],
    "genome_sizes" => [1000000, 5000000, 10000000],
    "annotation_methods" => ["ab_initio", "homology", "hybrid"],
    "n_genomes" => 20,
    "description" => "Medium scale - realistic datasets"
)

# Large scale configuration
large_config = Dict(
    "genome_types" => ["bacterial", "viral", "fungal", "plant"],
    "genome_sizes" => [1000000, 10000000, 100000000],
    "annotation_methods" => ["ab_initio", "homology", "hybrid", "rnaseq_guided"],
    "n_genomes" => 50,
    "description" => "Large scale - scalability testing"
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
println("Genome types: $(config["genome_types"])")
println("Genome sizes: $(config["genome_sizes"])")
println("Annotation methods: $(config["annotation_methods"])")

# ## Test Data Generation
#
# Download real reference genomes with pre-existing annotations from NCBI

println("\n--- Downloading Reference Genomes from NCBI ---")

# Include benchmark utilities
include("benchmark_utils.jl")

# Define reference genomes from multiple taxa/kingdoms/domains
# These are small, well-annotated genomes suitable for benchmarking
REFERENCE_GENOMES = Dict(
    # Bacteria - E. coli K-12 MG1655 (well-studied model organism)
    "bacterial" => Dict(
        "accession" => "GCF_000005845.2",
        "name" => "Escherichia coli K-12 MG1655",
        "expected_size" => 4_641_652  # ~4.6 Mb
    ),
    # Virus - PhiX174 (classic model phage)
    "viral" => Dict(
        "accession" => "GCF_000819615.1",
        "name" => "Enterobacteria phage phiX174",
        "expected_size" => 5386  # ~5.4 kb
    ),
    # Archaea - Methanocaldococcus jannaschii (first archaeal genome sequenced)
    "archaeal" => Dict(
        "accession" => "GCF_000091665.1",
        "name" => "Methanocaldococcus jannaschii DSM 2661",
        "expected_size" => 1_664_970  # ~1.7 Mb
    ),
    # Fungal - Saccharomyces cerevisiae S288C (yeast model)
    "fungal" => Dict(
        "accession" => "GCF_000146045.2",
        "name" => "Saccharomyces cerevisiae S288C",
        "expected_size" => 12_157_105  # ~12 Mb (only run in medium/large configs)
    ),
    # Plant organelle - Arabidopsis thaliana chloroplast (small plant genome proxy)
    "plant_organelle" => Dict(
        "accession" => "NC_000932.1",  
        "name" => "Arabidopsis thaliana chloroplast",
        "expected_size" => 154_478  # ~154 kb
    )
)

# Select which genomes to use based on configuration
genomes_to_download = if config == large_config
    ["bacterial", "viral", "archaeal", "fungal", "plant_organelle"]
elseif config == medium_config
    ["bacterial", "viral", "archaeal", "plant_organelle"]
else  # small_config
    ["bacterial", "viral"]
end

# Generate test datasets
test_data_dir = "annotation_test_data"
mkpath(test_data_dir)

test_genomes = []

for genome_type in genomes_to_download
    genome_info = REFERENCE_GENOMES[genome_type]
    accession = genome_info["accession"]
    
    println("\nDownloading $(genome_info["name"]) ($(accession))...")
    
    try
        # Use ncbi_genome_download_accession with retry logic
        dataset = Mycelia.with_retry(
            max_attempts=3,
            initial_delay=10.0,
            on_retry=(attempt, ex, delay) -> @warn "Download attempt $attempt/3 for $accession failed, retrying in $(delay)s..." exception=ex
        ) do
            Mycelia.ncbi_genome_download_accession(
                accession=accession,
                outdir=test_data_dir,
                include_string="gff3,cds,protein,genome,seq-report",
                max_attempts=1  # Disable internal retry since we have outer retry
            )
        end
        
        if dataset.genome !== nothing && isfile(dataset.genome)
            # Count genes from existing GFF3 annotation as ground truth
            expected_genes = 0
            if dataset.gff3 !== nothing && isfile(dataset.gff3)
                # Count CDS features in GFF3 as proxy for gene count
                gff_content = read(dataset.gff3, String)
                expected_genes = count(r"\tCDS\t", gff_content)
            end
            
            push!(test_genomes, Dict(
                "name" => "$(genome_type)_$(accession)",
                "organism_type" => genome_type,
                "organism_name" => genome_info["name"],
                "accession" => accession,
                "genome_size" => filesize(dataset.genome),
                "fasta_path" => dataset.genome,
                "gff3_path" => dataset.gff3,
                "cds_path" => dataset.cds,
                "protein_path" => dataset.protein,
                "expected_genes" => expected_genes,
                "has_reference_annotation" => dataset.gff3 !== nothing
            ))
            println("  ✓ Downloaded: $(dataset.genome)")
            println("    Genome size: $(round(filesize(dataset.genome) / 1e6, digits=2)) MB")
            println("    Reference genes (CDS count): $(expected_genes)")
        else
            @warn "Failed to download genome for $accession"
        end
    catch e
        @warn "Failed to download $accession after all retries" exception=e
    end
end

println("\nSuccessfully downloaded $(length(test_genomes)) reference genomes")

# ## Gene Prediction Benchmarks
#
# Test gene prediction performance and accuracy

println("\n--- Gene Prediction Benchmarks ---")

# Initialize benchmark suite
benchmark_suite = BenchmarkSuite("Annotation Performance")

# Test Pyrodigal gene prediction performance
for genome in test_genomes[1:min(3, length(test_genomes))]
    println("\nBenchmarking Pyrodigal gene prediction: $(genome["name"])")
    
    # Benchmark Pyrodigal gene prediction
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.run_pyrodigal;
        fasta_file=genome["fasta_path"],
        samples=3,
        seconds=120  # Allow up to 2 minutes per prediction
    )
    
    add_benchmark_result!(
        benchmark_suite,
        "pyrodigal_$(genome["name"])",
        benchmark_result,
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e9, digits=2)) seconds")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
    
    # Assess gene prediction results if available
    try
        pyrodigal_result = Mycelia.run_pyrodigal(fasta_file=genome["fasta_path"])
        
        if pyrodigal_result isa NamedTuple && haskey(pyrodigal_result, :gff) && isfile(pyrodigal_result.gff)
            predicted_genes = Mycelia.count_predicted_genes(pyrodigal_result.gff)
            
            println("  Gene prediction results:")
            println("    Predicted genes: $predicted_genes")
            println("    Expected genes: $(genome["expected_genes"])")
            
            # Calculate basic accuracy metrics
            sensitivity = min(predicted_genes / genome["expected_genes"], 1.0)
            
            # Add prediction quality metrics to results
            benchmark_suite.results["pyrodigal_$(genome["name"])"]["prediction_quality"] = Dict(
                "predicted_genes" => predicted_genes,
                "expected_genes" => genome["expected_genes"],
                "sensitivity" => round(sensitivity, digits=3)
            )
            
            println("    Sensitivity: $(round(sensitivity * 100, digits=1))%")
        end
    catch e
        println("  ⚠️  Gene prediction quality assessment failed: $e")
    end
end

# Use existing Mycelia function for counting predicted genes

# Test parallel gene prediction performance
if length(test_genomes) >= 2
    println("\nBenchmarking parallel Pyrodigal gene prediction")
    
    # Test parallel processing on multiple genomes
    fasta_files = [genome["fasta_path"] for genome in test_genomes[1:min(2, length(test_genomes))]]
    
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.parallel_pyrodigal, fasta_files;
        samples=2,
        seconds=300  # Allow up to 5 minutes for parallel processing
    )
    
    add_benchmark_result!(
        benchmark_suite,
        "parallel_pyrodigal_$(length(fasta_files))_genomes",
        benchmark_result,
        memory_stats
    )
    
    println("  Parallel processing median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e9, digits=2)) seconds")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
    
    # Calculate speedup vs sequential processing
    sequential_time = 0.0
    for genome in test_genomes[1:length(fasta_files)]
        test_key = "pyrodigal_$(genome["name"])"
        if haskey(benchmark_suite.results, test_key)
            sequential_time += benchmark_suite.results[test_key]["median_time"]
        end
    end
    
    if sequential_time > 0
        parallel_time = BenchmarkTools.median(benchmark_result).time
        speedup = sequential_time / parallel_time
        println("  Estimated speedup: $(round(speedup, digits=2))x")
        
        # Add speedup to results
        benchmark_suite.results["parallel_pyrodigal_$(length(fasta_files))_genomes"]["speedup"] = round(speedup, digits=2)
    end
end

# ## Functional Annotation Benchmarks
#
# Test functional annotation performance

println("\n--- Functional Annotation Benchmarks ---")

# TODO: Implement functional annotation benchmarks
# - Database search performance
# - Annotation transfer accuracy
# - GO term assignment quality
# - Pathway annotation completeness

# ## Comparative Annotation Benchmarks
#
# Test comparative annotation approaches

println("\n--- Comparative Annotation Benchmarks ---")

# TODO: Implement comparative annotation benchmarks
# - Ortholog identification accuracy
# - Synteny-based annotation transfer
# - Multi-genome annotation consistency
# - Phylogenetic annotation validation

# ## Annotation Pipeline Benchmarks
#
# Test complete annotation pipeline performance

println("\n--- Annotation Pipeline Benchmarks ---")

# TODO: Implement annotation pipeline benchmarks
# - End-to-end pipeline performance
# - Memory usage throughout pipeline
# - Intermediate result quality
# - Pipeline robustness

# ## Database Performance Benchmarks
#
# Test database search and annotation performance

println("\n--- Database Performance Benchmarks ---")

# TODO: Implement database performance benchmarks
# - BLAST search performance
# - Database size impact
# - Parallel database searches
# - Memory-efficient database usage

# ## Quality Assessment Benchmarks
#
# Test annotation quality assessment methods

println("\n--- Quality Assessment Benchmarks ---")

# TODO: Implement quality assessment benchmarks
# - Annotation completeness metrics
# - Consistency validation
# - Evidence-based quality scoring
# - Comparative quality assessment

# ## Scalability Benchmarks
#
# Test scalability with increasing genome sizes

println("\n--- Scalability Benchmarks ---")

# TODO: Implement scalability benchmarks
# - Performance vs genome size
# - Memory usage scaling
# - Parallel processing efficiency
# - Batch processing optimization

# ## Accuracy Benchmarks
#
# Test annotation accuracy using reference datasets

println("\n--- Accuracy Benchmarks ---")

# TODO: Implement accuracy benchmarks
# - Gene prediction accuracy (sensitivity, specificity)
# - Functional annotation accuracy
# - Comparative annotation validation
# - Cross-validation studies

# ## Results Collection
#
# Collect comprehensive benchmark results

println("\n--- Collecting Results ---")

# Calculate annotation performance metrics
total_genomes_processed = length(test_genomes)
total_genome_length = sum([genome["genome_size"] for genome in test_genomes])

# Gene prediction throughput estimates
prediction_throughput = Dict()
for (test_name, result) in benchmark_suite.results
    if startswith(test_name, "pyrodigal_") && haskey(result, "median_time")
        time_seconds = result["median_time"] / 1e9
        
        # Find corresponding genome
        genome_name = replace(test_name, "pyrodigal_" => "")
        matching_genome = findfirst(g -> g["name"] == genome_name, test_genomes)
        
        if matching_genome !== nothing
            genome = test_genomes[matching_genome]
            bp_per_second = genome["genome_size"] / time_seconds
            prediction_throughput[test_name] = round(bp_per_second, digits=0)
        end
    end
end

# Analyze prediction accuracy
prediction_accuracy_summary = Dict()
for (test_name, result) in benchmark_suite.results
    if haskey(result, "prediction_quality")
        quality = result["prediction_quality"]
        prediction_accuracy_summary[test_name] = quality
    end
end

# Create comprehensive results structure
results_dir = "results"
mkpath(results_dir)

comprehensive_results = Dict(
    "benchmark_name" => "annotation",
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
        "total_genomes_processed" => total_genomes_processed,
        "total_genome_length" => total_genome_length,
        "genome_sizes_tested" => config["genome_sizes"],
        "genome_types_tested" => config["genome_types"]
    ),
    "performance_metrics" => Dict(
        "prediction_throughput_bp_per_second" => prediction_throughput,
        "benchmark_duration_minutes" => round((Dates.now() - benchmark_suite.metadata["timestamp"]).value / 60000, digits=2)
    ),
    "accuracy_metrics" => prediction_accuracy_summary,
    "resource_utilization" => Dict(
        "peak_memory_mb" => round(maximum([get(result, "memory", 0) for result in values(benchmark_suite.results)]) / 1e6, digits=2),
        "total_allocations" => sum([get(result, "allocations", 0) for result in values(benchmark_suite.results)])
    ),
    "detailed_results" => benchmark_suite.results
)

# Save results
results_file = joinpath(results_dir, "annotation_benchmark_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")
save_benchmark_results(benchmark_suite, results_file)

# Save comprehensive results
comprehensive_file = joinpath(results_dir, "annotation_comprehensive_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")
open(comprehensive_file, "w") do f
    JSON.print(f, comprehensive_results, 2)
end

# Display summary
format_benchmark_summary(benchmark_suite)

println("\n--- Annotation Performance Summary ---")
println("Genomes processed: $total_genomes_processed")
println("Total genome length: $(round(total_genome_length / 1e6, digits=2)) Mbp")
println("Gene prediction throughput:")
for (test_name, throughput) in prediction_throughput
    println("  $test_name: $(Int(throughput)) bp/second")
end
println("Prediction accuracy:")
for (test_name, accuracy) in prediction_accuracy_summary
    println("  $test_name: $(accuracy["sensitivity"] * 100)% sensitivity")
end

# Cleanup test data
println("\n--- Cleanup ---")
for genome in test_genomes
    rm(genome["fasta_path"], force=true)
end
rm(test_data_dir, force=true, recursive=true)

println("\n=== Annotation Benchmark Complete ===")
println("Results saved to: $results_file")
println("Comprehensive results saved to: $comprehensive_file")
println("End time: $(Dates.now())")

nothing