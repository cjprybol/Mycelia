# # Tutorial 3: K-mer Analysis and Feature Extraction
#
# This tutorial explores k-mer analysis, a fundamental technique in bioinformatics
# for sequence analysis, genome assembly, and comparative genomics.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - K-mer theory and biological significance
# - How k-mer size affects analysis sensitivity and specificity
# - Dense vs sparse k-mer counting strategies
# - K-mer frequency spectra and their interpretation
# - Applications in genome size estimation and quality assessment
# - Memory and computational considerations for large-scale analysis

# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/03_kmer_analysis.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import FASTX
import Random
import Statistics

Random.seed!(42)

# ## Part 1: K-mer Theory and Biological Context
#
# K-mers are subsequences of length k extracted from DNA sequences.
# They capture local sequence composition and are fundamental to many algorithms.

println("=== K-mer Analysis Tutorial ===")

# ### K-mer Mathematics
#
# For DNA sequences, there are 4^k possible k-mers.
# Understanding k-mer space helps with parameter selection.

function kmer_space_size(k, alphabet_size=4)
    return alphabet_size^k
end

println("K-mer space sizes:")
for k in [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    space_size = kmer_space_size(k)
    println("k=$k: $(space_size) possible k-mers ($(space_size ÷ 1000)K)")
end

# ### Biological Significance
#
# K-mers capture different biological features depending on their size:
# - Small k-mers (k=3-7): Capture short motifs, sensitive to errors
# - Medium k-mers (k=15-21): Balance sensitivity and specificity
# - Large k-mers (k=25-51): Specific but may miss short overlaps

println("\nK-mer size selection guidelines:")
println("k=3-7:   Short motifs, codon analysis")
println("k=15-21: Error correction, initial assembly")
println("k=25-31: Genome assembly, repeat detection")
println("k=35-51: Specific overlaps, large genome assembly")

# ## Part 2: K-mer Counting Strategies
#
# Different applications require different counting approaches.
# Understanding trade-offs helps optimize performance.

println("\n=== K-mer Counting Strategies ===")

# Generate test sequences for demonstration
test_sequences = [
    Mycelia.random_fasta_record(moltype=:DNA, seed=i, L=1000) 
    for i in 1:10
]

# Write sequences to temporary files
temp_files = String[]
for (i, seq) in enumerate(test_sequences)
    filename = "test_seq_$i.fasta"
    Mycelia.write_fasta(outfile=filename, records=[seq])
    push!(temp_files, filename)
end

println("Generated $(length(temp_files)) test sequences")

# ### Dense K-mer Counting
#
# Dense counting stores all possible k-mers, including those not observed.
# Memory usage: O(4^k) - grows exponentially with k

println("\n--- Dense K-mer Counting ---")

dense_results = Dict{Int, NamedTuple}()
dense_summaries = Dict{Int, NamedTuple}()

for k in [3, 5, 7, 9]
    println("Computing dense k-mer counts for k=$k...")
    
    ## Memory estimation
    memory_mb = (4^k * 4) / (1024^2)  ## Assuming 4 bytes per count
    println("  Estimated memory: $(round(memory_mb, digits=2)) MB")
    
    if memory_mb < 100  ## Only run if memory usage is reasonable
        dense_counts = Mycelia.fasta_list_to_dense_kmer_counts(
            fasta_list=temp_files, 
            alphabet=:DNA, 
            k=k
        )
        println("  ✓ Dense counting completed")
        println("  Matrix size: $(size(dense_counts.counts))")
        dense_results[k] = dense_counts

        kmer_totals = vec(sum(dense_counts.counts; dims=2))
        kmer_hist = Mycelia.kmer_frequency_histogram(kmer_totals)
        zero_kmers = get(kmer_hist, 0, 0)
        peak = Mycelia.coverage_peak_from_hist(kmer_hist; min_coverage=2)
        dense_bytes = Mycelia.estimate_dense_matrix_memory(
            eltype(dense_counts.counts),
            size(dense_counts.counts, 1),
            size(dense_counts.counts, 2)
        )
        dense_summaries[k] = (
            total_kmers=sum(kmer_totals),
            observed_kmers=length(kmer_totals) - zero_kmers,
            zero_fraction=zero_kmers / max(1, length(kmer_totals)),
            peak_coverage=peak.coverage,
            dense_bytes=dense_bytes
        )
        println("  Observed k-mers: $(dense_summaries[k].observed_kmers)")
        println("  Zero fraction: $(round(dense_summaries[k].zero_fraction, digits=3))")
        println("  Peak coverage (>=2): $(dense_summaries[k].peak_coverage)")
        println("  Estimated dense bytes: $(Mycelia.bytes_human_readable(dense_summaries[k].dense_bytes))")
    else
        println("  ⚠ Skipping due to high memory usage")
    end
end

# ### Sparse K-mer Counting
#
# Sparse counting only stores observed k-mers.
# Memory usage: O(n) where n is number of unique k-mers

println("\n--- Sparse K-mer Counting ---")

sparse_example = nothing
sparse_example_k = 15
sparse_summaries = Dict{Int, NamedTuple}()

for k in [11, 13, 15, 17, 19, 21]
    println("Computing sparse k-mer counts for k=$k...")
    
    sparse_counts = Mycelia.fasta_list_to_sparse_kmer_counts(
        fasta_list=temp_files,
        alphabet=:DNA,
        k=k,
        skip_rarefaction_plot=true,
        show_rarefaction_plot=false
    )
    println("  ✓ Sparse counting completed")
    nnz = Mycelia.SparseArrays.nnz(sparse_counts.counts)
    total_entries = prod(size(sparse_counts.counts))
    density = total_entries == 0 ? 0.0 : nnz / total_entries
    sparse_bytes = Mycelia.estimate_sparse_matrix_memory(
        eltype(sparse_counts.counts),
        size(sparse_counts.counts, 1),
        size(sparse_counts.counts, 2);
        nnz=nnz
    )
    sparse_summaries[k] = (
        unique_kmers=length(sparse_counts.kmers),
        nnz=nnz,
        density=density,
        sparse_bytes=sparse_bytes
    )
    println("  Unique k-mers: $(sparse_summaries[k].unique_kmers)")
    println("  Matrix density: $(round(sparse_summaries[k].density, digits=5))")
    println("  Estimated sparse bytes: $(Mycelia.bytes_human_readable(sparse_summaries[k].sparse_bytes))")

    if k == sparse_example_k
        sparse_example = sparse_counts
    end
end

# ## Part 3: K-mer Frequency Spectra
#
# K-mer frequency spectra reveal genome characteristics and data quality

println("\n=== K-mer Frequency Spectra ===")

spectrum_k = 7
canonical_kmer_counts = Mycelia.count_canonical_kmers(
    Mycelia.Kmers.DNAKmer{spectrum_k},
    temp_files
)
kmer_counts_vector = collect(values(canonical_kmer_counts))
kmer_hist = Mycelia.kmer_frequency_histogram(kmer_counts_vector)
coverage_peak = Mycelia.coverage_peak_from_hist(kmer_hist; min_coverage=2)

println("Spectrum histogram bins: $(length(kmer_hist))")
println("Coverage peak (>=2): $(coverage_peak.coverage) with $(coverage_peak.kmers) k-mers")

repeat_threshold = coverage_peak.coverage === missing ? 5 : max(coverage_peak.coverage * 3, 5)
repeat_like = count(c -> c >= repeat_threshold, kmer_counts_vector)
println("Repeat-like k-mers (>= $repeat_threshold): $repeat_like")

spectrum_dir = joinpath(@__DIR__, "..", "results", "tutorial_03_kmer_spectra")
Base.Filesystem.mkpath(spectrum_dir)
spectrum_plot = Mycelia.plot_kmer_frequency_spectra(
    kmer_counts_vector;
    log_scale=log2,
    title="K-mer spectrum (k=$(spectrum_k))"
)
Mycelia.StatsPlots.savefig(
    spectrum_plot,
    joinpath(spectrum_dir, "kmer_spectrum_k$(spectrum_k).png")
)
println("Saved spectrum plot to $(joinpath(spectrum_dir, "kmer_spectrum_k$(spectrum_k).png"))")

# ## Part 4: Applications in Genome Analysis
#
# K-mers have many applications in genomic analysis

println("\n=== Genome Analysis Applications ===")

# ### Genome Size Estimation
#
# Use k-mer frequency spectra to estimate genome size
# Formula: Genome size ≈ Total k-mers / Coverage peak

println("--- Genome Size Estimation ---")

total_kmers = sum(kmer_counts_vector)
estimated_size = coverage_peak.coverage === missing ? missing : Int(round(total_kmers / coverage_peak.coverage))
known_size = sum(length(FASTX.sequence(record)) for record in test_sequences)
basic_estimate = Mycelia.estimate_genome_size_from_kmers(test_sequences, spectrum_k)

println("Total k-mers: $total_kmers")
println("Estimated genome size (coverage peak): $estimated_size")
println("Known total sequence length: $known_size")
if estimated_size !== missing
    size_error = abs(estimated_size - known_size) / max(1, known_size)
    println("Relative error: $(round(size_error * 100, digits=2))%")
end
println("Basic k-mer estimate: $(basic_estimate["estimated_genome_size"]) (k=$(spectrum_k))")

# ### Error Detection and Correction
#
# Low-frequency k-mers often represent sequencing errors

println("--- Error Detection ---")

canonical_counts_dict = Dict(canonical_kmer_counts)
filtered_kmers, clustering_stats, removed_kmers = Mycelia.automatic_error_filtering(canonical_counts_dict)
singleton_kmers = count(c -> c == 1, kmer_counts_vector)
rare_kmers = count(c -> c > 1 && c <= (coverage_peak.coverage === missing ? 3 : max(2, coverage_peak.coverage ÷ 2)), kmer_counts_vector)

println("Singleton k-mers: $singleton_kmers")
println("Rare k-mers: $rare_kmers")
println("Clustering separation quality: $(round(clustering_stats.separation_quality, digits=3))")
println("Suggested coverage threshold: $(clustering_stats.optimal_threshold)")
println("Clustered low-coverage k-mers removed: $(length(removed_kmers))")
println("Clustered high-confidence k-mers retained: $(length(filtered_kmers))")

error_graph = Mycelia.Rhizomorph.build_kmer_graph(
    test_sequences,
    spectrum_k;
    dataset_id="tutorial_error",
    mode=:doublestrand
)
low_coverage_vertices = Mycelia.Rhizomorph.find_low_coverage_kmers(error_graph, 1)
high_coverage_vertices = Mycelia.Rhizomorph.find_high_coverage_kmers(error_graph, 2)
println("Graph low-coverage vertices (<=1): $(length(low_coverage_vertices))")
println("Graph high-coverage vertices (>=2): $(length(high_coverage_vertices))")

# ### Contamination Detection
#
# Foreign DNA creates distinctive k-mer patterns

println("--- Contamination Detection ---")

contaminant_sequences = [
    Mycelia.random_fasta_record(moltype=:DNA, seed=200 + i, L=1000)
    for i in 1:3
]
primary_counts = Mycelia.count_canonical_kmers(
    Mycelia.Kmers.DNAKmer{spectrum_k},
    test_sequences
)
contaminant_counts = Mycelia.count_canonical_kmers(
    Mycelia.Kmers.DNAKmer{spectrum_k},
    contaminant_sequences
)

all_kmers = sort(collect(union(keys(primary_counts), keys(contaminant_counts))))
primary_vec = [get(primary_counts, kmer, 0) for kmer in all_kmers]
contam_vec = [get(contaminant_counts, kmer, 0) for kmer in all_kmers]
cosine_dist = Mycelia.Distances.cosine_dist(primary_vec, contam_vec)
js_div = Mycelia.Distances.js_divergence(
    primary_vec ./ max(1, sum(primary_vec)),
    contam_vec ./ max(1, sum(contam_vec))
)

foreign_kmers = setdiff(keys(contaminant_counts), keys(primary_counts))
foreign_fraction = length(foreign_kmers) / max(1, length(keys(contaminant_counts)))

println("Cosine distance (primary vs contaminant): $(round(cosine_dist, digits=3))")
println("JS divergence (primary vs contaminant): $(round(js_div, digits=3))")
println("Foreign k-mer fraction: $(round(foreign_fraction * 100, digits=2))%")

# ## Part 5: Performance Optimization
#
# Large-scale k-mer analysis requires optimization

println("\n=== Performance Optimization ===")

# ### Memory Management
#
# Strategies for handling large k-mer datasets

println("--- Memory Management ---")

dense_estimated_bytes = Mycelia.estimate_dense_matrix_memory(UInt32, 4^11, length(temp_files))
memory_check = Mycelia.check_matrix_fits_in_memory(dense_estimated_bytes; severity=:warn)
println("Dense k=11 estimate: $(Mycelia.bytes_human_readable(dense_estimated_bytes))")
println("Memory available: $(Mycelia.bytes_human_readable(memory_check.free_memory))")

cache_dir = Base.Filesystem.mktempdir()
cache_file = joinpath(cache_dir, "sparse_counts_k$(sparse_example_k).jld2")
if sparse_example !== nothing
    Mycelia.save_kmer_results(
        filename=cache_file,
        kmers=sparse_example.kmers,
        counts=sparse_example.counts,
        fasta_list=temp_files,
        k=sparse_example_k,
        alphabet=:DNA
    )
    cached = Mycelia.load_kmer_results(cache_file)
    if cached !== nothing
        println("Loaded cached counts: $(size(cached.counts)) from $cache_file")
    end
end

hist_df = Mycelia.DataFrames.DataFrame(
    coverage=collect(keys(kmer_hist)),
    kmers=collect(values(kmer_hist))
)
histogram_path = joinpath(cache_dir, "kmer_histogram.tsv.gz")
Mycelia.write_tsvgz(df=hist_df, filename=histogram_path, force=true)
println("Compressed histogram saved to $histogram_path")

mmap_path = joinpath(cache_dir, "kmer_counts.bin")
open(mmap_path, "w") do io
    write(io, kmer_counts_vector)
end
mapped_counts = Mycelia.Mmap.mmap(mmap_path, Vector{eltype(kmer_counts_vector)}, (length(kmer_counts_vector),))
println("Memory-mapped counts length: $(length(mapped_counts))")

# ### Parallel Processing
#
# Accelerate k-mer counting with parallelization

println("--- Parallel Processing ---")

println("Threads available: $(Threads.nthreads())")
threaded_sparse = Mycelia.fasta_list_to_sparse_kmer_counts(
    fasta_list=temp_files,
    alphabet=:DNA,
    k=spectrum_k,
    force_threading=true,
    skip_rarefaction_plot=true,
    show_rarefaction_plot=false
)
println("Threaded sparse counts: $(size(threaded_sparse.counts))")

chunk_a = temp_files[1:5]
chunk_b = temp_files[6:10]
counts_a = Mycelia.count_canonical_kmers(Mycelia.Kmers.DNAKmer{spectrum_k}, chunk_a)
counts_b = Mycelia.count_canonical_kmers(Mycelia.Kmers.DNAKmer{spectrum_k}, chunk_b)
merged_counts = merge!(+, counts_a, counts_b)
println("Merged counts size (map-reduce style): $(length(merged_counts))")

# ## Part 6: Visualization and Interpretation
#
# Create plots to understand k-mer analysis results

println("\n=== K-mer Visualization ===")

if haskey(dense_results, 3)
    heatmap_plot = Mycelia.StatsPlots.heatmap(
        dense_results[3].counts;
        xlabel="Sequence",
        ylabel="3-mer index",
        title="3-mer composition heatmap"
    )
    Mycelia.StatsPlots.savefig(
        heatmap_plot,
        joinpath(spectrum_dir, "kmer_heatmap_k3.png")
    )
    println("Saved 3-mer heatmap to $(joinpath(spectrum_dir, "kmer_heatmap_k3.png"))")
end

coverage_plot = Mycelia.StatsPlots.histogram(
    kmer_counts_vector;
    bins=:auto,
    xlabel="k-mer count",
    ylabel="Number of k-mers",
    title="Coverage distribution (k=$(spectrum_k))"
)
Mycelia.StatsPlots.savefig(
    coverage_plot,
    joinpath(spectrum_dir, "coverage_distribution_k$(spectrum_k).png")
)

comparative_k = 7
sequence_kmer_counts = [
    Mycelia.count_canonical_kmers(Mycelia.Kmers.DNAKmer{comparative_k}, [record])
    for record in test_sequences
]
similarity_matrix = zeros(length(sequence_kmer_counts), length(sequence_kmer_counts))
for i in eachindex(sequence_kmer_counts)
    for j in eachindex(sequence_kmer_counts)
        all_k = union(keys(sequence_kmer_counts[i]), keys(sequence_kmer_counts[j]))
        a = [get(sequence_kmer_counts[i], kmer, 0) for kmer in all_k]
        b = [get(sequence_kmer_counts[j], kmer, 0) for kmer in all_k]
        similarity_matrix[i, j] = 1 - Mycelia.Distances.cosine_dist(a, b)
    end
end
similarity_plot = Mycelia.StatsPlots.heatmap(
    similarity_matrix;
    xlabel="Sequence index",
    ylabel="Sequence index",
    title="K-mer cosine similarity (k=$(comparative_k))"
)
Mycelia.StatsPlots.savefig(
    similarity_plot,
    joinpath(spectrum_dir, "kmer_similarity_k$(comparative_k).png")
)
println("Saved visualization plots to $spectrum_dir")

# ## Part 7: Advanced K-mer Techniques
#
# Explore advanced k-mer analysis methods

println("\n=== Advanced Techniques ===")

# ### Minimizers
#
# Reduce k-mer space using minimizer techniques

println("--- Minimizers ---")

minimizer_k = 9
minimizer_window = 10
syncmer_s = 3
syncmer_t = 2
strobe_w_min = 1
strobe_w_max = 5

example_sequence = FASTX.sequence(Mycelia.BioSequences.LongDNA{4}, test_sequences[1])
minimizers = Mycelia.canonical_minimizers(example_sequence, minimizer_k, minimizer_window)
syncmers = Mycelia.open_syncmers(example_sequence, minimizer_k, syncmer_s, syncmer_t; canonical=true)
strobes = Mycelia.strobemers(example_sequence, minimizer_k, strobe_w_min, strobe_w_max; canonical=true)

println("Canonical minimizers: $(length(minimizers))")
println("Open syncmers: $(length(syncmers))")
println("Strobemers: $(length(strobes))")

# ### Graph Construction
#
# Build graphs from k-mer overlaps

println("--- Graph Construction ---")

graph_k = 5
kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(
    test_sequences,
    graph_k;
    dataset_id="tutorial",
    mode=:doublestrand
)
println("K-mer graph vertices: $(Mycelia.Graphs.nv(kmer_graph))")
println("K-mer graph edges: $(Mycelia.Graphs.ne(kmer_graph))")
high_coverage_kmers = Mycelia.Rhizomorph.find_high_coverage_kmers(kmer_graph, 2)
println("High-coverage k-mers (>=2): $(length(high_coverage_kmers))")

# ## Summary and Best Practices
println("\n=== K-mer Analysis Summary ===")
println("✓ Understanding k-mer theory and biological significance")
println("✓ Choosing appropriate k-mer sizes for different applications")
println("✓ Implementing dense and sparse counting strategies")
println("✓ Analyzing k-mer frequency spectra")
println("✓ Applying k-mer analysis to genome size estimation")
println("✓ Optimizing performance for large-scale analysis")

println("\nBest Practices:")
println("- Start with k=21 for general analysis")
println("- Use dense counting for small k, sparse for large k")
println("- Monitor memory usage and optimize accordingly")
println("- Validate results with known datasets")
println("- Consider biological context in interpretation")

# Cleanup
for file in temp_files
    rm(file, force=true)
end

println("\nNext: Tutorial 4 - Genome Assembly")

nothing
