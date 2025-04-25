"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the quality score for a single base given multiple observations.

This function implements the "Converting to Error Probabilities and Combining" method:
1. Takes error probabilities from multiple reads covering the same base
2. Calculates probability of ALL reads being wrong by multiplying probabilities
3. Calculates final Phred score from this combined probability

To avoid numerical underflow with very small probabilities, the calculation
is performed in log space.

# Arguments
- `error_probabilities::Vector{Float64}`: Vector of error probabilities from 
  multiple reads covering the same base position

# Returns
- `Float64`: Phred quality score representing the combined confidence
"""
function joint_base_quality_score(error_probabilities::Vector{Float64})
    if isempty(error_probabilities)
        return 0.0  # No data available
    end
    
    # Work in log space to avoid underflow
    log_p_all_wrong = sum(log.(error_probabilities))
    
    # Convert back from log space
    p_all_wrong = exp(log_p_all_wrong)
    
    # Prevent underflow/overflow
    if p_all_wrong <= eps(Float64)
        return 999.0  # Cap at a very high quality score
    end
    
    # Convert to Phred score
    return Mycelia.error_rate_to_q_value(p_all_wrong)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate kmer quality score using the specified aggregation method.

Available methods:
- `:min`: Use the minimum base quality (default)
- `:mean`: Use the mean of all base qualities
- `:geometric`: Use the geometric mean (appropriate for probabilities)
- `:harmonic`: Use the harmonic mean (emphasizes lower values)

# Arguments
- `base_qualities::Vector{Float64}`: Vector of quality scores for each base
- `method::Symbol`: Method to use for aggregation

# Returns
- `Float64`: Overall quality score for the kmer
"""
function kmer_quality_score(base_qualities::Vector{Float64}, method::Symbol=:min)
    if isempty(base_qualities)
        return 0.0
    end
    
    if method == :min
        return minimum(base_qualities)
    elseif method == :mean
        return mean(base_qualities)
    elseif method == :geometric
        # Convert to probabilities, compute geometric mean, convert back
        error_probs = [phred_to_error_prob(q) for q in base_qualities]
        geo_mean_prob = exp(sum(log.(error_probs)) / length(error_probs))
        return error_prob_to_phred(geo_mean_prob)
    elseif method == :harmonic
        # Harmonic mean emphasizes lower values
        # which is appropriate for quality scores
        return length(base_qualities) / sum(1.0 ./ base_qualities)
    else
        error("Unknown aggregation method: $method")
    end
end

# conda install -c bioconda kmer-jellyfish
# count, bc, info, stats, histo, dump, merge, query, cite, mem, jf
# cap at 4 threads, 8Gb per thread by default - this should be plenty fast enough for base usage, but open it up for higher performance!
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers in a FASTA/FASTQ file using Jellyfish.

# Arguments
- `fastx::String`: Path to input FASTA/FASTQ file (can be gzipped)
- `k::Integer`: k-mer length
- `threads::Integer=Sys.CPU_THREADS`: Number of threads to use
- `max_mem::Integer=Int(Sys.free_memory())`: Maximum memory in bytes (defaults to system free memory)
- `canonical::Bool=false`: Whether to count canonical k-mers (both strands combined)
- `outfile::String=auto`: Output filename (auto-generated based on input and parameters)
- `conda_check::Bool=true`: Whether to verify Jellyfish conda installation

# Returns
- `String`: Path to gzipped TSV file containing k-mer counts
"""
function jellyfish_count(;fastx, k, threads=Sys.CPU_THREADS, max_mem=Int(Sys.free_memory()), canonical=false, outfile = ifelse(canonical, "$(fastx).k$(k).canonical.jf", "$(fastx).k$(k).jf"), conda_check=true)
    if conda_check
        Mycelia.add_bioconda_env("kmer-jellyfish")
    end
    mem = Int(floor(max_mem * 0.8))
    
    jellyfish_mem_output = read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish mem --mer-len $(k) --mem $(mem)`, String)
    # sample output = "68719476736 (68G)\n"
    # this version grabs the exact number at the beginning
    jellyfish_buffer_size = first(split(strip(jellyfish_mem_output)))
    
    # this grabs the human readable version in parentheses
    jellyfish_buffer_size = replace(split(strip(jellyfish_mem_output))[2], r"[\(\)]" => "")
    @show jellyfish_buffer_size

    @info "making a temporary copy of the input fastx"
    temp_fastx = copy_to_tempdir(fastx)
    if occursin(r"\.gz$", temp_fastx)
        run(`gzip -d $(temp_fastx)`)
        temp_fastx = replace(temp_fastx, r"\.gz$" => "")
    end
    @assert isfile(temp_fastx)
    if canonical
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish count --canonical --size $(jellyfish_buffer_size) --threads $(threads) --mer-len $(k) --output $(outfile) $(temp_fastx)`
    else
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish count --size $(jellyfish_buffer_size) --threads $(threads) --mer-len $(k) --output $(outfile) $(temp_fastx)`
    end

    temp_tab = outfile * ".tsv"
    tabular_counts = temp_tab * ".gz"

    if !isfile(tabular_counts)
        if !isfile(outfile)
            @info "running kmer counts"
            @show cmd
            run(cmd)
            @info "done counting kmers"
        end
        if !isfile(temp_tab)
            @info "dumping counts to tab-delimited file"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish dump --column --tab --output $(temp_tab) $(outfile)`)
            @info "done dumping counts to tab-delimited file"
        end
        run(`gzip $(temp_tab)`)
    end

    isfile(outfile) && rm(outfile)
    isfile(temp_fastx) && rm(temp_fastx)
    return tabular_counts
end

# 7
# 0.169580 seconds (587.93 k allocations: 39.458 MiB, 76.55% compilation time)
# 11
# 4.420040 seconds (194.49 k allocations: 11.566 MiB)
# 13
# 20.690993 seconds (521.75 k allocations: 37.253 MiB, 0.58% compilation time)
# 17
# 412.670050 seconds (529.86 k allocations: 37.007 MiB, 0.03% compilation time)
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a Jellyfish k-mer count file into a frequency histogram.

# Arguments
- `jellyfish_counts_file::String`: Path to the gzipped TSV file containing Jellyfish k-mer counts
- `outfile::String=replace(jellyfish_counts_file, r"\\.tsv\\.gz\$" => ".count_histogram.tsv")`: Optional output file path

# Returns
- `String`: Path to the generated histogram file

# Description
Processes a Jellyfish k-mer count file to create a frequency histogram where:
- Column 1: Number of k-mers that share the same count
- Column 2: The count they share

Uses system sorting with LC_ALL=C for optimal performance on large files.

Notes
- Requires gzip, sort, uniq, and sed command line tools
- Uses intermediate disk storage for sorting large files
- Skips processing if output file already exists
"""
function jellyfish_counts_to_kmer_frequency_histogram(jellyfish_counts_file, outfile=replace(jellyfish_counts_file, r"\.tsv\.gz$" => ".count_histogram.tsv"))
    # sorting with LC_ALL=C is the biggest speed up here of anything I've found
    if !isfile(outfile)
        io = open(pipeline(
                `gzip -dc $(jellyfish_counts_file)`,
                `cut -f2`,
                Cmd(`sort --temporary-directory . --compress-program gzip --numeric --stable`, env=Dict("LC_ALL" => "C")),
                `uniq --count`,
                `sed 's/^ *//'`,
                `sed 's/ /\t/'`
                ))
        frequency_histogram_table = CSV.read(io, DataFrames.DataFrame, header=["number of kmers", "number of observations"], delim='\t')
        CSV.write(outfile, frequency_histogram_table, delim='\t')
    else
        @info "$(outfile) already exists"
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load k-mer counts from a Jellyfish output file into a DataFrame.

# Arguments
- `jellyfish_counts::String`: Path to a gzipped TSV file (*.jf.tsv.gz) containing Jellyfish k-mer counts

# Returns
- `DataFrame`: Table with columns:
  - `kmer`: Biologically encoded k-mers as `DNAKmer{k}` objects
  - `count`: Integer count of each k-mer's occurrences

# Notes
- Input file must be a gzipped TSV with exactly two columns (k-mer sequences and counts)
- K-mer length is automatically detected from the first entry
- Filename must end with '.jf.tsv.gz'
"""
function load_jellyfish_counts(jellyfish_counts)
    @assert occursin(r"\.jf\.tsv\.gz$", jellyfish_counts)
    open(jellyfish_counts) do io
        table = CSV.read(CodecZlib.GzipDecompressorStream(io), DataFrames.DataFrame, delim="\t", header=["kmer", "count"])
        k = length(table[1, "kmer"])
        table[!, "kmer"] = map(x -> Kmers.DNAKmer{k}(String(x)), collect(table[!, "kmer"]))
        return table
    end
end

# Usage: jellyfish merge [options] input:string+

# Merge jellyfish databases

# Options (default value in (), *required):
#  -o, --output=string                      Output file (mer_counts_merged.jf)
#  -m, --min                                Compute min count instead of sum (false)
#  -M, --max                                Compute max count instead of sum (false)
#  -j, --jaccard                            Compute the jaccard and weighted jaccard similarities (false)
#  -L, --lower-count=uint64                 Don't output k-mer with count < lower-count
#  -U, --upper-count=uint64                 Don't output k-mer with count > upper-count
#      --usage                              Usage
#  -h, --help                               This message
#  -V, --version                            Version

# Usage: jellyfish query [options] file:path mers:string+

# Query a Jellyfish database

# Options (default value in (), *required):
#  -s, --sequence=path                      Output counts for all mers in sequence
#  -o, --output=path                        Output file (stdout)
#  -i, --interactive                        Interactive, queries from stdin (false)
#  -l, --load                               Force pre-loading of database file into memory (false)
#  -L, --no-load                            Disable pre-loading of database file into memory (false)
#  -U, --usage                              Usage
#  -h, --help                               This message
#  -V, --version                            Version

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Counts k-mer occurrences in a FASTA file, considering both forward and reverse complement sequences.

# Arguments
- `kmer_type`: Type specification for k-mers (e.g., `DNAKmer{21}`)
- `fasta`: Path to FASTA file containing reference sequences

# Returns
- `Dict{kmer_type, Int}`: Dictionary mapping each k-mer to its total count across all sequences
"""
function fasta_to_reference_kmer_counts(;kmer_type, fasta)
    kmer_counts = Dict{kmer_type, Int}()
    for record in Mycelia.open_fastx(fasta)
        record_sequence = BioSequences.LongDNA{2}(FASTX.sequence(record))
        forward_counts = StatsBase.countmap(kmer for (i, kmer) in Kmers.EveryKmer{kmer_type}(record_sequence))
        reverse_counts = StatsBase.countmap(kmer for (i, kmer) in Kmers.EveryKmer{kmer_type}(BioSequences.reverse_complement(record_sequence)))
        record_counts = merge(+, forward_counts, reverse_counts)
        merge!(+, kmer_counts, record_counts)
    end
    return kmer_counts
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the cosine similarity between two k-mer count dictionaries.

# Arguments
- `kmer_counts_1::Dict{String,Int}`: First dictionary mapping k-mer sequences to their counts
- `kmer_counts_2::Dict{String,Int}`: Second dictionary mapping k-mer sequences to their counts

# Returns
- `Float64`: Cosine distance between the two k-mer count vectors, in range [0,1]
  where 0 indicates identical distributions and 1 indicates maximum dissimilarity

# Details
Converts k-mer count dictionaries into vectors using a unified set of keys,
then computes cosine distance. Missing k-mers are treated as count 0.
Result is invariant to input order and total counts (normalized internally).
"""
function kmer_counts_to_cosine_similarity(kmer_counts_1, kmer_counts_2)
    sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
    a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    # Distances.cosine_dist(a, b) == Distances.cosine_dist(b, a) == Distances.cosine_dist(a ./ sum(a), b ./ sum(b)) == Distances.cosine_dist(b ./ sum(b), a ./ sum(a))
    return Distances.cosine_dist(a, b)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Jensen-Shannon divergence between two k-mer frequency distributions.

# Arguments
- `kmer_counts_1`: Dictionary mapping k-mers to their counts in first sequence
- `kmer_counts_2`: Dictionary mapping k-mers to their counts in second sequence

# Returns
- Normalized Jensen-Shannon divergence score between 0 and 1, where:
  - 0 indicates identical distributions
  - 1 indicates maximally different distributions

# Notes
- The measure is symmetric: JS(P||Q) = JS(Q||P)
- Counts are automatically normalized to probability distributions
"""
function kmer_counts_to_js_divergence(kmer_counts_1, kmer_counts_2)
    sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
    a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    a_norm = a ./ sum(a)
    b_norm = b ./ sum(b)
    # Distances.js_divergence(a ./ sum(a), b ./ sum(b)) == Distances.js_divergence(b ./ sum(b), a ./ sum(a))
    # Distances.js_divergence(a, b) != Distances.js_divergence(a ./ sum(a), b ./ sum(b))
    return Distances.js_divergence(a_norm, b_norm)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Evaluate genome assembly quality by comparing k-mer distributions between assembled sequences and raw observations.

# Arguments
- `assembly`: Input assembled sequences to evaluate
- `observations`: Raw sequencing data for comparison
- `ks::Vector{Int}`: Vector of k-mer sizes to analyze (default: k=17 to 23)

# Returns
DataFrame containing quality metrics for each k-mer size:
- `k`: K-mer length used
- `cosine_distance`: Cosine similarity between k-mer distributions
- `js_divergence`: Jensen-Shannon divergence between distributions  
- `qv`: MerQury-style quality value score
"""
function assess_assembly_kmer_quality(;assembly, observations, ks::Vector{Int}=filter(x -> 17 <= x <= 23, Mycelia.ks()))
    results = DataFrames.DataFrame()
    @show ks
    ProgressMeter.@showprogress for k in ks
        @show k
        @info "counting assembly kmers..."
        assembled_canonical_kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, assembly)
        # assembled_canonical_kmer_counts_file = Mycelia.jellyfish_count(fastx=assembly, k=k, canonical=true)
        # @info "loading assembly kmer counts..."
        # assembled_canonical_kmer_counts_table = Mycelia.load_jellyfish_counts(assembled_canonical_kmer_counts_file)
        # sort!(assembled_canonical_kmer_counts_table, "kmer")
        # assembled_canonical_kmer_counts = OrderedCollections.OrderedDict(row["kmer"] => row["count"] for row in DataFrames.eachrow(assembled_canonical_kmer_counts_table))
        
        @info "counting observation kmers..."
        observed_canonical_kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, observations)
        # observed_canonical_kmer_counts_file = Mycelia.jellyfish_count(fastx=observations, k=k, canonical=true)
        # @info "loading observation kmer counts..."
        # observed_canonical_kmer_counts_table = Mycelia.load_jellyfish_counts(observed_canonical_kmer_counts_file)
        # sort!(observed_canonical_kmer_counts_table, "kmer")
        # observed_canonical_kmer_counts = OrderedCollections.OrderedDict(row["kmer"] => row["count"] for row in DataFrames.eachrow(observed_canonical_kmer_counts_table))
        cosine_distance = kmer_counts_to_cosine_similarity(observed_canonical_kmer_counts, assembled_canonical_kmer_counts)
        js_divergence = kmer_counts_to_js_divergence(observed_canonical_kmer_counts, assembled_canonical_kmer_counts)
        qv = kmer_counts_to_merqury_qv(raw_data_counts=observed_canonical_kmer_counts, assembly_counts=assembled_canonical_kmer_counts)
        push!(results, (;k, cosine_distance, js_divergence, qv))
    end
    return results
end

# function assess_assembly_quality(;assembled_sequence::BioSequences.LongDNA{2}, fastq::String, k::Int)
#     assess_assembly_quality(assembled_sequence=assembled_sequence, fastq=fastq, ks=[k])
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the Jaccard similarity coefficient between two sets.

The Jaccard similarity is defined as the size of the intersection divided by the size
of the union of two sets:

    J(A,B) = |A ∩ B| / |A ∪ B|

# Arguments
- `set1`: First set for comparison
- `set2`: Second set for comparison

# Returns
- `Float64`: A value between 0.0 and 1.0, where:
  * 1.0 indicates identical sets
  * 0.0 indicates completely disjoint sets
"""
function jaccard_similarity(set1, set2)
    return length(intersect(set1, set2)) / length(union(set1, set2))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Jaccard distance between two sets, which is the complement of the Jaccard similarity.

The Jaccard distance is defined as:
``J_d(A,B) = 1 - J_s(A,B) = 1 - \\frac{|A ∩ B|}{|A ∪ B|}``

# Arguments
- `set1`: First set to compare
- `set2`: Second set to compare

# Returns
- `Float64`: A value in [0,1] where 0 indicates identical sets and 1 indicates disjoint sets
"""
function jaccard_distance(set1, set2)
    return 1.0 - jaccard_similarity(set1, set2)
end

# function kmer_counts_to_jaccard(kmer_counts_1::AbstractDict{Kmers.DNAKmer{k}, Int64}, kmer_counts_2::AbstractDict{Kmers.DNAKmer{k}, Int64}) where k
#     # sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
#     # a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
#     # b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
#     # a_indices = findall(a .> 0)
#     # b_indices = findall(b .> 0)
#     # return Distances.jaccard(a_indices, b_indices)
#     return jaccard(collect(keys(kmer_counts_1)), collect(keys(kmer_counts_2)))
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate assembly Quality Value (QV) score using the Merqury method.

Estimates base-level accuracy by comparing k-mer distributions between raw sequencing
data and assembly. Higher QV scores indicate better assembly quality.

# Arguments
- `raw_data_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}`: K-mer counts from raw sequencing data
- `assembly_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}`: K-mer counts from assembly

# Returns
- `Float64`: Quality Value score in Phred scale (-10log₁₀(error rate))

# Method
QV is calculated using:
1. Ktotal = number of unique kmers in assembly
2. Kshared = number of kmers shared between raw data and assembly
3. P = (Kshared/Ktotal)^(1/k) = estimated base-level accuracy
4. QV = -10log₁₀(1-P)

# Reference
Rhie et al. "Merqury: reference-free quality, completeness, and phasing assessment
for genome assemblies" Genome Biology (2020)
"""
function kmer_counts_to_merqury_qv(;raw_data_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}, assembly_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}) where {k,N}
    # Ktotal = # of kmers found in assembly
    Ktotal = length(keys(assembly_counts))
    # Kshared = # of shared kmers between assembly and readset
    Kshared = length(intersect(keys(raw_data_counts), keys(assembly_counts)))
    # probabilitiy_base_in_assembly_correct
    P = (Kshared/Ktotal)^(1/k)
    # # Error rate
    E = 1-P
    QV = -10log10(E)
    # return (;P, E, QV)
    return QV
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyzes k-mer frequency spectra from sequencing reads to detect sequencing errors and 
coverage patterns.

# Arguments
- `out_directory`: Directory where output files will be saved
- `forward_reads`: Path to forward reads FASTQ file
- `reverse_reads`: Path to reverse reads FASTQ file (optional)
- `k`: Length of k-mers to analyze (default: 17)
- `target_coverage`: Desired coverage level. If non-zero, returns downsampling rate (default: 0)
- `plot_size`: Dimensions of output plot in pixels (default: (600,400))

# Returns
- If `target_coverage > 0`: Returns downsampling rate as Float64
- If `target_coverage = 0`: Returns nothing, saves plots only

# Outputs
- `peak-detected.png`: K-mer spectra plot with detected boundaries
- `peak-detected.svg`: Vector format of the same plot
- `downsampling-rate.txt`: Created only if target_coverage is specified

# Algorithm
1. Counts canonical k-mers from input reads
2. Generates log-scale histogram of k-mer frequencies
3. Detects error threshold using regression or local minimum
4. Identifies coverage peak after error threshold
5. Calculates downsampling rate if target coverage specified
"""
# https://github.com/cjprybol/Mycelia/blob/e7fe50ffe2d18406fb70e0e24ebcfa45e0937596/notebooks/exploratory/2021-08-25-k-medoids-error-cluster-detection-multi-entity-graph-aligner-test.ipynb
function analyze_kmer_spectra(;out_directory, forward_reads="", reverse_reads="", k=17, target_coverage=0, plot_size=(600,400))
    @info "counting $k-mers"
    user_provided_reads = filter(x -> !isempty(x), [forward_reads, reverse_reads])
    canonical_kmer_counts = count_canonical_kmers(Kmers.DNAKmer{k}, user_provided_reads)

    @info "determining max count"
    max_count = maximum(values(canonical_kmer_counts))
    @info "max count = $max_count"

    @info "generating histogram"
    kmer_counts_histogram = sort(collect(StatsBase.countmap(values(canonical_kmer_counts))), by=x->x[1])

    X = log2.(first.(kmer_counts_histogram))
    Y = log2.(last.(kmer_counts_histogram))
    
    @info "plotting kmer spectra"
    p = StatsPlots.scatter(
        X,
        Y,
        xlabel="log2(kmer_frequency)",
        ylabel="log2(# of kmers @ frequency)",
        label="",
        size=plot_size
    )

    earliest_y_min_index = last(findmin(Y))
    lower_boundary = X[earliest_y_min_index]
    lower_boundary_source = "first minimum"

    try
        # take the first 1/denominator datapoints in the set
        # to capture the error line on the left side of the graph
        @info "fitting error curve"
        denominators = [2^i for i in 1:5]
        coeficient_matrix = zeros(length(denominators), 2)
        for (i, denominator) in enumerate(denominators)
            prefix_index = Int(floor(length(X)/denominator))
            _x = X[1:prefix_index]
            _y = Y[1:prefix_index]
            model = GLM.lm(GLM.@formula(Y ~ X), DataFrames.DataFrame(X = _x, Y = _y))
            coeficient_matrix[i, :] = GLM.coef(model)
        end
        median_intercept = Statistics.median(coeficient_matrix[:, 1])
        median_slope = Statistics.median(coeficient_matrix[:, 2])

        X_intercept = (0 - median_intercept) / median_slope

        # some libraries detect the x_intercept being AFTER the end of the data
        # in these instances detect the earliest x-minimum
        if X_intercept < lower_boundary
            lower_boundary = X_intercept
            lower_boundary_source = "detected x-intercept"
        end
    catch
        @info "unable to fit regression"
    end

    p = StatsPlots.vline!(p,
        [lower_boundary],
        label="lower boundary ($(lower_boundary_source))"
    );
    
    is_above_lower_bounds = X .>= lower_boundary
    max_Y_post_error_intercept = first(findmax(Y[is_above_lower_bounds]))
    peak_indices = findall(is_above_lower_bounds .& (Y .== max_Y_post_error_intercept))
    peak_index = Int(round(Statistics.median(peak_indices)))

    p = StatsPlots.vline!([X[peak_index]], label="inferred sample coverage)")
    if isinteractive()
        display(p)
    end
    StatsPlots.savefig(p, "$out_directory/peak-detected.png")
    StatsPlots.savefig(p, "$out_directory/peak-detected.svg")
    
    if target_coverage != 0
        detected_coverage = 2^(X[peak_index])
        downsampling_rate = round(target_coverage/detected_coverage, sigdigits=3)
        downsampling_rate = min(downsampling_rate, 1)
        @info "downsampling rate = $downsampling_rate"

        outfile = "$out_directory/downsampling-rate.txt"
        open(outfile, "w") do io
            @info "writing downsampling rate to $outfile"
            println(io, downsampling_rate)
        end
        return downsampling_rate
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assess k-mer saturation in DNA sequences from FASTX files.

# Arguments
- `fastxs::AbstractVector{<:AbstractString}`: Vector of paths to FASTA/FASTQ files
- `kmer_type`: Type of k-mer to analyze (e.g., DNAKmer{21})
- `kmers_to_assess=Inf`: Maximum number of k-mers to process
- `power=10`: Base for exponential sampling intervals
- `min_count=1`: Minimum count threshold for considering a k-mer

# Returns
Named tuple containing:
- `sampling_points::Vector{Int}`: K-mer counts at which samples were taken
- `unique_kmer_counts::Vector{Int}`: Number of unique canonical k-mers at each sampling point
- `eof::Bool`: Whether the entire input was processed

# Details
Analyzes k-mer saturation by counting unique canonical k-mers at exponentially spaced 
intervals (powers of `power`). Useful for assessing sequence complexity and coverage.
Returns early if all possible k-mers are observed.
"""
function assess_dnamer_saturation(fastxs::AbstractVector{<:AbstractString}, kmer_type; kmers_to_assess=Inf, power=10, min_count = 1)
    # canonical_kmers = Set{kmer_type}()
    canonical_kmer_counts = Dict{kmer_type, Int}()
    
    @show kmer_type
    k = Kmers.ksize(Kmers.kmertype(kmer_type))
    
    max_possible_kmers = determine_max_canonical_kmers(k, DNA_ALPHABET)
    
    if kmers_to_assess == Inf
        # want to read the whole file and predict how long that will take
        # n_records = reduce(sum, map(f -> Mycelia.count_records(f), fastxs))
        kmers_to_assess = max_possible_kmers
        p = ProgressMeter.Progress(kmers_to_assess, 1)
    else
        p = ProgressMeter.Progress(kmers_to_assess, 1)
    end
    
    sampling_points = Int[0]
    i = 0
    while power^i <= kmers_to_assess
        push!(sampling_points, power^i)
        i += 1
    end
    
    unique_kmer_counts = zeros(Int, length(sampling_points))
    
    if length(sampling_points) < 3
        @info "increase the # of reads analyzed or decrease the power to acquire more data points"
        return (;sampling_points, unique_kmer_counts)
    end
    
    kmers_assessed = 0
    for fastx in fastxs
        for record in open_fastx(fastx)      
            record_sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)
            for (index, kmer) in Kmers.EveryKmer{kmer_type}(record_sequence)
                canonical_kmer = BioSequences.canonical(kmer)
                if haskey(canonical_kmer_counts, canonical_kmer)
                    canonical_kmer_counts[canonical_kmer] += 1
                else
                    canonical_kmer_counts[canonical_kmer] = 1
                end
                kmers_assessed += 1
                if (length(canonical_kmer_counts) == max_possible_kmers)                 
                    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
                    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], length(canonical_kmer_counts))
                    return (;sampling_points, unique_kmer_counts, eof = false)
                elseif kmers_assessed in sampling_points
                    i = findfirst(sampling_points .== kmers_assessed)
                    unique_kmer_counts[i] = length(filter(x -> x[2] >= min_count, canonical_kmer_counts))
                    if i == length(sampling_points)
                        return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = false)
                    end
                end
                ProgressMeter.next!(p)
            end
        end
    end
    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], [length(canonical_kmer_counts)])    
    return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = true)
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyze k-mer saturation in DNA sequences to determine optimal k value.

# Arguments
- `fastxs`: Vector of paths to FASTA/FASTQ files to analyze
- `power`: Base of logarithmic sampling points (default: 10)
- `outdir`: Optional output directory for plots and results
- `min_k`: Minimum k-mer size to test (default: 7)
- `max_k`: Maximum k-mer size to test (default: 17)
- `threshold`: Saturation threshold to determine optimal k (default: 0.1)
- `kmers_to_assess`: Maximum number of k-mers to sample (default: 10M)
- `plot`: Whether to generate saturation curves (default: true)

# Returns
Integer representing the first k value that achieves saturation below threshold.
If no k value meets the threshold, returns the k with minimum saturation.

# Details
- Tests only prime k values between min_k and max_k
- Generates saturation curves using logarithmic sampling
- Fits curves to estimate maximum unique k-mers
- If outdir is provided, saves plots as SVG and chosen k value to text file
"""
function assess_dnamer_saturation(fastxs::AbstractVector{<:AbstractString}; power=10, outdir::Union{Missing, String}=missing, min_k=7, max_k=17, threshold=0.1, kmers_to_assess=10_000_000, plot=true)
    ks = Primes.primes(min_k, max_k)
    minimum_saturation = Inf
    midpoint = Inf
    for k in ks
        kmer_type = Kmers.kmertype(Kmers.Kmer{BioSequences.DNAAlphabet{4},k})
        sampling_points, kmer_counts, hit_eof = assess_dnamer_saturation(fastxs, kmer_type, kmers_to_assess=kmers_to_assess, power=power)
        @show sampling_points, kmer_counts, hit_eof
        observed_midpoint_index = findfirst(i -> kmer_counts[i] > last(kmer_counts)/2, 1:length(sampling_points))
        observed_midpoint = sampling_points[observed_midpoint_index]
        initial_parameters = Float64[maximum(kmer_counts), observed_midpoint]
        @time fit = LsqFit.curve_fit(calculate_v, sampling_points, kmer_counts, initial_parameters)
        max_canonical_kmers = determine_max_canonical_kmers(k, DNA_ALPHABET)
        if hit_eof
            inferred_maximum = last(kmer_counts)
        else
            inferred_maximum = max(Int(ceil(fit.param[1])), last(kmer_counts))
            if inferred_maximum > max_canonical_kmers
                inferred_maximum = max_canonical_kmers
            end
        end

        inferred_midpoint = Int(ceil(fit.param[2]))
        predicted_saturation = inferred_maximum / max_canonical_kmers
        @show k, predicted_saturation

        if plot
            scale = 300
            fontsize = 14
            p = StatsPlots.scatter(
                sampling_points,
                kmer_counts,
                label="observed counts",
                ylabel="# unique canonical kmers",
                xlabel="# kmers assessed",
                title = "sequencing saturation @ k = $k",
                # legend=:outertopright,
                # size=(2*scale, 1*scale),
                margins=5Plots.PlotMeasures.mm,
                titlefontsize=fontsize,
                xguidefontsize=fontsize,
                yguidefontsize=fontsize,
                legendfontsize=fontsize-2,
                xtickfontsize=fontsize-6,
                ytickfontsize=fontsize-4,
                # xrotation=45
                )
            StatsPlots.hline!(p, [max_canonical_kmers], label="absolute maximum", line = :solid, linewidth = 2)
            StatsPlots.hline!(p, [inferred_maximum], label="inferred maximum", line = :dash, linewidth = 2)
            StatsPlots.vline!(p, [inferred_midpoint], label="inferred midpoint", line = :dot, linewidth = 2)
            # xs = vcat(sampling_points, [last(sampling_points) * 2^i for i in 1:2])
            xs = sort([sampling_points..., inferred_midpoint])
            ys = calculate_v(xs, fit.param)
            StatsPlots.plot!(
                p,
                xs,
                ys,
                label="fit trendline",
                line=:dashdot,
                linewidth = 2)
            if k != first(ks)
                StatsPlots.plot!(p, legend=false)
            end
            display(p)
            if !ismissing(outdir)
                # StatsPlots.savefig(p, joinpath(outdir, "$k.png"))
                StatsPlots.savefig(p, joinpath(outdir, "$k.svg"))
            end
        end
            

        if predicted_saturation < minimum_saturation
            minimum_saturation = predicted_saturation
            min_k = k
            midpoint = inferred_midpoint 
        end
        if predicted_saturation < threshold
            if !ismissing(outdir)
                chosen_k_file = joinpath(outdir, "chosen_k.txt")
                println("chosen k = $k")
                open(chosen_k_file, "w") do io
                    println(io, k)
                end
            end
            return k
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyzes k-mer saturation in a FASTA/FASTQ file to determine optimal k-mer size.

# Arguments
- `fastx::AbstractString`: Path to input FASTA/FASTQ file
- `power::Int=10`: Exponent for downsampling k-mers (2^power)
- `outdir::String=""`: Output directory for results. Uses current directory if empty
- `min_k::Int=3`: Minimum k-mer size to evaluate
- `max_k::Int=17`: Maximum k-mer size to evaluate
- `threshold::Float64=0.1`: Saturation threshold for k-mer assessment
- `kmers_to_assess::Int=10_000_000`: Maximum number of k-mers to sample

# Returns
`Dict{Int,Float64}`: Dictionary mapping k-mer sizes to their saturation scores
"""
function assess_dnamer_saturation(fastx::AbstractString; power=10, outdir="", min_k=3, max_k=17, threshold=0.1, kmers_to_assess=10_000_000)
    assess_dnamer_saturation([fastx], outdir=outdir, min_k=min_k, max_k=max_k, threshold=threshold, power=power, kmers_to_assess=kmers_to_assess)
end