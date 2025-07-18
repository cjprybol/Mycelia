"""
$(DocStringExtensions.TYPEDSIGNATURES)

Save the kmer counting results (kmers vector, counts sparse matrix)
and the input FASTA file list to a JLD2 file for long-term storage
and reproducibility.

# Arguments
- `filename::AbstractString`: Path to the output JLD2 file.
- `kmers::AbstractVector{<:Kmers.Kmer}`: The sorted vector of unique kmer objects.
- `counts::AbstractMatrix`: The (sparse or dense) matrix of kmer counts.
- `fasta_list::AbstractVector{<:AbstractString}`: The list of FASTA file paths used as input.
- `k::Integer`: The kmer size used.
- `alphabet::Symbol`: The alphabet used (:AA, :DNA, :RNA).

"""
function save_kmer_results(;
    filename::AbstractString,
    kmers::AbstractVector{<:Kmers.Kmer},
    counts::AbstractMatrix,
    fasta_list::AbstractVector{<:AbstractString},
    k::Integer,
    alphabet::Symbol
    )

    if !endswith(filename, ".jld2")
        @warn "Filename '$filename' does not end with '.jld2'. Appending it."
        filename *= ".jld2"
    end

    @info "Saving kmer results to $filename..."
    try
        # Open the file in write mode ('w'). This will create or overwrite the file.
        JLD2.jldopen(filename, "w") do file
            # Write each component as a separate dataset within the JLD2 file
            # The string names ("kmers", "counts", etc.) are how you'll access them upon loading.
            file["kmers"] = kmers
            file["counts"] = counts
            file["fasta_list"] = fasta_list
            # Store metadata as well
            file["metadata/k"] = k
            file["metadata/alphabet"] = string(alphabet) # Store symbol as string for robustness
            file["metadata/creation_timestamp"] = string(Dates.now())
            file["metadata/julia_version"] = string(VERSION)
        end
        @info "Successfully saved results to $filename."
    catch e
        @error "Failed to save results to $filename: $e"
        rethrow(e)
    end
    return nothing
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load kmer counting results previously saved with `save_kmer_results`.

# Arguments
- `filename::AbstractString`: Path to the input JLD2 file.

# Returns
- `NamedTuple`: Contains the loaded `kmers`, `counts`, `fasta_list`, and `metadata`.
  Returns `nothing` if the file cannot be loaded or essential keys are missing.
"""
function load_kmer_results(filename::AbstractString)
    if !Base.Filesystem.isfile(filename)
        @error "File not found: $filename"
        return nothing
    end
    if !endswith(filename, ".jld2")
        @warn "Filename '$filename' does not end with '.jld2'. Attempting to load anyway."
    end

    @info "Loading kmer results from $filename..."
    try
        # JLD2.load returns a dictionary-like object mapping names to loaded data
        loaded_data = JLD2.load(filename)

        # Basic validation (check if essential keys exist)
        required_keys = ["kmers", "counts", "fasta_list"]
        if !all(haskey(loaded_data, k) for k in required_keys)
             @error "File $filename is missing one or more required keys: $required_keys"
             return nothing
        end

        # Reconstruct metadata if present
        metadata = Dict{String, Any}()
        if haskey(loaded_data, "metadata/k")
            metadata["k"] = loaded_data["metadata/k"]
        end
        if haskey(loaded_data, "metadata/alphabet")
             # Attempt to convert back to Symbol, fallback to string if error
             try
                 metadata["alphabet"] = Symbol(loaded_data["metadata/alphabet"])
             catch
                 @warn "Could not convert loaded alphabet '$(loaded_data["metadata/alphabet"])' back to Symbol. Storing as String."
                 metadata["alphabet"] = loaded_data["metadata/alphabet"]
             end
        end
         # Add other metadata fields if they exist
        for key in keys(loaded_data)
            if startswith(key, "metadata/") && !haskey(metadata, replace(key, "metadata/" => ""))
                 metadata[replace(key, "metadata/" => "")] = loaded_data[key]
            end
        end


        @info "Successfully loaded results from $filename."
        # Return as a NamedTuple for convenient access
        return (;
            kmers=loaded_data["kmers"],
            counts=loaded_data["counts"],
            fasta_list=loaded_data["fasta_list"],
            metadata=metadata
        )

    catch e
        @error "Failed to load results from $filename: $e"
        return nothing
    end
end


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
    # probability_base_in_assembly_correct
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

"""
    generate_and_save_kmer_counts(; 
        bioalphabet, 
        fastas, 
        k, 
        output_dir=pwd(),
        filename=nothing
    )

Generates and saves k-mer counts for a list of FASTA files for a single k.

# Keyword Arguments
- `bioalphabet`: Alphabet type (e.g., :DNA).
- `fastas`: List of FASTA file paths.
- `k`: Value of k (e.g., 9).
- `output_dir`: (optional) Directory to write output files (default: current directory).
- `filename`: (optional) Full file name for output (default: 
    `"{Mycelia.normalized_current_date()}.{lowercase(string(bioalphabet))}{k}mers.jld2"`).

# Output
Saves a .jld2 file with the specified file name in `output_dir` if it does not already exist.
"""
function generate_and_save_kmer_counts(; 
    alphabet, 
    fastas, 
    k, 
    output_dir=pwd(),
    filename=nothing
)
    if filename === nothing
        filename = string(
            Mycelia.normalized_current_date(), ".", 
            lowercase(string(alphabet)), 
            k, "mers.jld2"
        )
    end
    kmer_result_file = joinpath(output_dir, filename)
    if !isfile(kmer_result_file)
        if (alphabet == :DNA) || (alphabet == :RNA)
            if k <= 9
                kmer_count_results = Mycelia.fasta_list_to_dense_kmer_counts(
                    fasta_list=fastas,
                    k=k,
                    alphabet=alphabet
                )
            else
                kmer_count_results = Mycelia.fasta_list_to_sparse_kmer_counts(
                    fasta_list=fastas,
                    k=k,
                    alphabet=alphabet
                )
            end
        elseif (alphabet == :AA)
            if k <= 3
                kmer_count_results = Mycelia.fasta_list_to_dense_kmer_counts(
                    fasta_list=fastas,
                    k=k,
                    alphabet=alphabet
                )
            else
                kmer_count_results = Mycelia.fasta_list_to_sparse_kmer_counts(
                    fasta_list=fastas,
                    k=k,
                    alphabet=alphabet
                )
            end
        else
            error("unrecognized alphabet: $(alphabet)")
        end
        Mycelia.save_kmer_results(
            filename = kmer_result_file,
            kmers = kmer_count_results.kmers,
            counts = kmer_count_results.counts,
            fasta_list = fastas,
            k = k,
            alphabet = alphabet
        )
    end
    return kmer_result_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a sparse kmer counts table (SparseMatrixCSC) from a list of FASTA files using a 3-pass approach.
Pass 1 (Parallel): Counts kmers per file and writes to temporary JLD2 files.
Pass 2 (Serial): Aggregates unique kmers, max count, nnz per file, and rarefaction data from temp files.
                 Generates and saves a k-mer rarefaction plot.
Pass 3 (Parallel): Reads temporary counts again to construct the final sparse matrix.

Optionally, a results filename can be provided to save/load the output. If the file exists and `force` is false,
the result is loaded and returned. If `force` is true or the file does not exist, results are computed and saved.

# Output Directory Behavior
- All auxiliary output files (e.g., rarefaction data, plots) are written to a common output directory.
- By default, this is:
    - The value of `out_dir` if provided.
    - Otherwise, the directory containing `result_file` (if provided and has a directory component).
    - Otherwise, the current working directory (`pwd()`).
- If you provide an absolute path for an output file (e.g. `rarefaction_data_filename`), that path is used directly.
- If both `out_dir` and a relative filename are given, the file is written to `out_dir`.

# Arguments
- `fasta_list::AbstractVector{<:AbstractString}`: A list of paths to FASTA files.
- `k::Integer`: The length of the kmer.
- `alphabet::Symbol`: The alphabet type (:AA, :DNA, :RNA).
- `temp_dir_parent::AbstractString`: Parent directory for creating the temporary working directory. Defaults to `Base.tempdir()`.
- `count_element_type::Union{Type{<:Unsigned}, Nothing}`: Optional. Specifies the unsigned integer type for the counts. If `nothing` (default), the smallest `UInt` type capable of holding the maximum observed count is used.
- `rarefaction_data_filename::AbstractString`: Filename for the TSV output of rarefaction data. If a relative path, will be written to `out_dir`.
- `rarefaction_plot_basename::AbstractString`: Basename for the output rarefaction plots. If a relative path, will be written to `out_dir`.
- `show_rarefaction_plot::Bool`: Whether to display the rarefaction plot after generation. Defaults to `true`.
- `result_file::Union{Nothing, AbstractString}`: Optional. If provided, path to a file to save/load the full results (kmers, counts, etc) as a JLD2 file.
- `out_dir::Union{Nothing, AbstractString}`: Optional. Output directory for auxiliary outputs. Defaults as described above.
- `force::Bool`: If true, recompute and overwrite the output file even if it exists. Defaults to `false`.
- `rarefaction_plot_kwargs...`: Keyword arguments to pass to `plot_kmer_rarefaction` for plot customization.

# Returns
- `NamedTuple{(:kmers, :counts, :rarefaction_data_path)}`:
    - `kmers`: A sorted `Vector` of unique kmer objects.
    - `counts`: A `SparseArrays.SparseMatrixCSC{V, Int}` storing kmer counts.
    - `rarefaction_data_path`: Path to the saved TSV file with rarefaction data.

# Raises
- `ErrorException`: If input `fasta_list` is empty, alphabet is invalid, or required Kmer/counting functions are not found.
"""
function fasta_list_to_sparse_kmer_counts(;
    fasta_list::AbstractVector{<:AbstractString},
    k::Integer,
    alphabet::Symbol,
    temp_dir_parent::AbstractString = Base.tempdir(),
    count_element_type::Union{Type{<:Unsigned}, Nothing} = nothing,
    rarefaction_data_filename::AbstractString = "$(normalized_current_datetime()).$(lowercase(string(alphabet)))$(k)mer_rarefaction.tsv",
    rarefaction_plot_basename::AbstractString = replace(rarefaction_data_filename, ".tsv" => ""),
    show_rarefaction_plot::Bool = true,
    result_file::Union{Nothing, AbstractString} = nothing,
    out_dir::Union{Nothing, AbstractString} = nothing,
    force::Bool = false,
    rarefaction_plot_kwargs...
)

    # --- Determine output directory logic ---
    function _get_output_dir()
        if !isnothing(out_dir)
            return out_dir
        elseif !isnothing(result_file) && !isempty(Base.Filesystem.dirname(result_file)) && Base.Filesystem.dirname(result_file) != "."
            return Base.Filesystem.dirname(result_file)
        else
            return Base.pwd()
        end
    end
    output_dir = _get_output_dir()
    if !Base.Filesystem.isdir(output_dir)
        Base.Filesystem.mkpath(output_dir)
    end
    Base.@info "Output directory for auxiliary files: $output_dir"

    # Helper to prepend output_dir only if the path is relative
    function _resolve_outpath(fname::AbstractString)
        (Base.Filesystem.isabspath(fname) || isempty(fname)) ? fname : Base.Filesystem.joinpath(output_dir, fname)
    end

    # --- Results File Short Circuit ---
    if !isnothing(result_file) && Base.Filesystem.isfile(result_file) && !force
        Base.@info "result_file already exists at $result_file. Loading and returning results."
        result = JLD2.load_object(result_file)
        if all(haskey(result, key) for key in (:kmers, :counts, :rarefaction_data_path))
            return result
        else
            Base.@warn "Loaded file does not have required structure; will recompute."
        end
    elseif !isnothing(result_file) && force && Base.Filesystem.isfile(result_file)
        Base.@info "Force is true. Overwriting $result_file with recomputed results."
    end

    # --- 0. Input Validation and Setup ---
    num_files = length(fasta_list)
    if num_files == 0
        error("Input fasta_list is empty.")
    end

    KMER_TYPE = if alphabet == :AA
        isdefined(Kmers, :AAKmer) ? Kmers.AAKmer{k} : error("Kmers.AAKmer not found or Kmers.jl not loaded correctly.")
    elseif alphabet == :DNA
        isdefined(Kmers, :DNAKmer) ? Kmers.DNAKmer{k} : error("Kmers.DNAKmer not found or Kmers.jl not loaded correctly.")
    elseif alphabet == :RNA
        isdefined(Kmers, :RNAKmer) ? Kmers.RNAKmer{k} : error("Kmers.RNAKmer not found or Kmers.jl not loaded correctly.")
    else
        error("Invalid alphabet: $alphabet. Choose from :AA, :DNA, :RNA")
    end

    COUNT_FUNCTION = if alphabet == :DNA
        func_name = :count_canonical_kmers
        isdefined(Main, func_name) ? getfield(Main, func_name) : (isdefined(@__MODULE__, func_name) ? getfield(@__MODULE__, func_name) : error("$func_name not found"))
    else
        func_name = :count_kmers
        isdefined(Main, func_name) ? getfield(Main, func_name) : (isdefined(@__MODULE__, func_name) ? getfield(@__MODULE__, func_name) : error("$func_name not found"))
    end

    temp_dir = Base.Filesystem.mktempdir(temp_dir_parent; prefix="kmer_counts_3pass_")
    Base.@info "Using temporary directory for intermediate counts: $temp_dir"
    temp_file_paths = [Base.Filesystem.joinpath(temp_dir, "counts_$(i).jld2") for i in 1:num_files]

    # --- Pass 1: Count Kmers per file and Write to Temp Files (Parallel) ---
    Base.@info "Pass 1: Counting $(KMER_TYPE) for $num_files files and writing temps using $(Base.Threads.nthreads()) threads..."
    progress_pass1 = ProgressMeter.Progress(num_files; desc="Pass 1 (Counting & Writing Temps): ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:cyan)
    progress_pass1_lock = Base.ReentrantLock()

    Base.Threads.@threads for original_file_idx in 1:num_files
        fasta_file = fasta_list[original_file_idx]
        temp_filename = temp_file_paths[original_file_idx]
        counts_dict = Dict{KMER_TYPE, Int}()
        try
            counts_dict = COUNT_FUNCTION(KMER_TYPE, fasta_file)
            JLD2.save_object(temp_filename, counts_dict)
        catch e
            Base.println("Error processing file $fasta_file (idx $original_file_idx) during Pass 1: $e")
            try
                JLD2.save_object(temp_filename, Dict{KMER_TYPE, Int}())
            catch save_err
                Base.@error "Failed to save empty placeholder for errored file $fasta_file to $temp_filename: $save_err"
            end
        end
        Base.lock(progress_pass1_lock) do
            ProgressMeter.next!(progress_pass1)
        end
    end
    ProgressMeter.finish!(progress_pass1)
    Base.@info "Pass 1 finished."

    # --- Pass 2: Aggregate Unique Kmers, Max Count, NNZ per file, Rarefaction Data (Serial) ---
    Base.@info "Pass 2: Aggregating stats and rarefaction data from temporary files (serially)..."
    all_kmers_set = Set{KMER_TYPE}()
    max_observed_count = 0
    total_non_zero_entries = 0
    nnz_per_file = Vector{Int}(undef, num_files)
    rarefaction_points = Vector{Tuple{Int, Int}}()

    progress_pass2 = ProgressMeter.Progress(num_files; desc="Pass 2 (Aggregating Serially): ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:yellow)

    for i in 1:num_files
        temp_filename = temp_file_paths[i]
        current_file_nnz = 0
        try
            if Base.Filesystem.isfile(temp_filename)
                counts_dict = JLD2.load_object(temp_filename)
                current_file_nnz = length(counts_dict)
                if current_file_nnz > 0
                    for k in keys(counts_dict)
                        push!(all_kmers_set, k)
                    end
                    current_file_max_val = maximum(values(counts_dict))
                    if current_file_max_val > max_observed_count
                        max_observed_count = current_file_max_val
                    end
                end
            else
                Base.@warn "Temporary file not found during Pass 2: $temp_filename. Assuming 0 counts for this file."
            end
        catch e
            Base.@error "Error reading or processing temporary file $temp_filename (for original file $i) during Pass 2: $e. Assuming 0 counts for this file."
        end
        nnz_per_file[i] = current_file_nnz
        total_non_zero_entries += current_file_nnz
        push!(rarefaction_points, (i, length(all_kmers_set)))
        ProgressMeter.next!(progress_pass2; showvalues = [(:unique_kmers, length(all_kmers_set))])
    end
    ProgressMeter.finish!(progress_pass2)
    Base.@info "Pass 2 aggregation finished."

    # --- Process and Save/Plot Rarefaction Data (after Pass 2) ---
    rarefaction_data_path = _resolve_outpath(rarefaction_data_filename)
    Base.@info "Saving rarefaction data to $rarefaction_data_path..."
    try
        data_to_write = [ [pt[1], pt[2]] for pt in rarefaction_points ]
        if !isempty(data_to_write)
            DelimitedFiles.writedlm(rarefaction_data_path, data_to_write, '\t')
        else
            Base.@warn "No rarefaction points recorded. TSV file will be empty or not created."
            DelimitedFiles.writedlm(rarefaction_data_path, Array{Int}(undef,0,2), '\t')
        end
    catch e
        Base.@error "Failed to write rarefaction data to $rarefaction_data_path: $e"
    end

    rarefaction_plot_base = _resolve_outpath(rarefaction_plot_basename)
    if Base.Filesystem.isfile(rarefaction_data_path) && !isempty(rarefaction_points)
        Base.@info "Generating k-mer rarefaction plot..."
        try
            plot_kmer_rarefaction(
                rarefaction_data_path;
                output_dir = Base.Filesystem.dirname(rarefaction_plot_base),
                output_basename = Base.Filesystem.basename(rarefaction_plot_base),
                display_plot = show_rarefaction_plot,
                rarefaction_plot_kwargs...
            )
        catch e
            Base.@error "Failed to generate rarefaction plot: $e. Ensure Makie and a backend are correctly set up. Also ensure 'plot_kmer_rarefaction' function is loaded."
        end
    else
        Base.@warn "Skipping rarefaction plot generation as data file is missing or empty."
    end

    # --- Determine Value Type, Prepare for Pass 3 ---
    ValType = if isnothing(count_element_type)
        if max_observed_count <= typemax(UInt8)
            UInt8
        elseif max_observed_count <= typemax(UInt16)
            UInt16
        elseif max_observed_count <= typemax(UInt32)
            UInt32
        else
            UInt64
        end
    else
        count_element_type
    end
    if !isnothing(count_element_type) && max_observed_count > typemax(count_element_type)
        Base.@warn "User-specified count_element_type ($count_element_type) may be too small for the maximum observed count ($max_observed_count)."
    end
    Base.@info "Using element type $ValType for kmer counts. Max observed count: $max_observed_count."

    if isempty(all_kmers_set) && total_non_zero_entries == 0
        Base.@warn "No kmers found across any files, or errors prevented aggregation."
        result = (; kmers=Vector{KMER_TYPE}(), counts=SparseArrays.spzeros(ValType, Int, 0, num_files), rarefaction_data_path=rarefaction_data_path)
        if !isnothing(result_file)
            JLD2.save_object(result_file, result)
            Base.@info "Saved empty results to $result_file"
        end
        return result
    end

    sorted_kmers = sort(collect(all_kmers_set))
    num_kmers = length(sorted_kmers)
    empty!(all_kmers_set); all_kmers_set = nothing; GC.gc()

    # Estimate memory needed for sparse matrix
    sparse_matrix_bytes_needed = Mycelia.estimate_sparse_matrix_memory(ValType, num_kmers, num_files, nnz=total_non_zero_entries)

    # Check if matrix will fit in memory
    mem_check = check_matrix_fits_in_memory(sparse_matrix_bytes_needed; severity=:error)
    if !mem_check.will_fit_total
        error("Sparse matrix will not fit in available memory. Required: $(bytes_human_readable(sparse_matrix_bytes_needed)), Available: $(bytes_human_readable(mem_check.free_memory))")
    end

    kmer_to_row_map = Dict{KMER_TYPE, Int}(kmer => i for (i, kmer) in enumerate(sorted_kmers))
    Base.@info "Found $num_kmers unique kmers. Total non-zero entries: $total_non_zero_entries."

    # --- Pass 3: Prepare Sparse Matrix Data (In Parallel from Temp Files) ---
    Base.@info "Pass 3: Preparing data for sparse matrix construction using $(Base.Threads.nthreads()) threads..."

    actual_total_nnz_from_files = sum(nnz_per_file)
    if total_non_zero_entries != actual_total_nnz_from_files
        Base.@warn "Sum of nnz_per_file ($actual_total_nnz_from_files) differs from serially accumulated total_non_zero_entries ($total_non_zero_entries). Using sum of nnz_per_file."
        total_non_zero_entries = actual_total_nnz_from_files
    end

    if total_non_zero_entries == 0 && num_kmers > 0
        Base.@warn "Found $num_kmers unique kmers, but $total_non_zero_entries non-zero entries. Matrix will be empty of values."
    elseif total_non_zero_entries == 0 && num_kmers == 0
        Base.@warn "No k-mers and no non-zero entries. Resulting matrix will be empty."
    end

    row_indices = Vector{Int}(undef, total_non_zero_entries)
    col_indices = Vector{Int}(undef, total_non_zero_entries)
    values_vec = Vector{ValType}(undef, total_non_zero_entries)

    write_offsets = Vector{Int}(undef, num_files + 1)
    write_offsets[1] = 0
    for i in 1:num_files
        write_offsets[i+1] = write_offsets[i] + nnz_per_file[i]
    end

    progress_pass3 = ProgressMeter.Progress(num_files; desc="Pass 3 (Filling Sparse Data): ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:blue)
    progress_pass3_lock = Base.ReentrantLock()

    Base.Threads.@threads for original_file_idx in 1:num_files
        if nnz_per_file[original_file_idx] > 0
            temp_filename = temp_file_paths[original_file_idx]
            file_start_offset = write_offsets[original_file_idx]
            current_entry_in_file = 0
            try
                if Base.Filesystem.isfile(temp_filename)
                    counts_dict_pass3 = JLD2.load_object(temp_filename)
                    for (kmer, count_val) in counts_dict_pass3
                        if count_val > 0
                            row_idx_val = Base.get(kmer_to_row_map, kmer, 0)
                            if row_idx_val > 0
                                global_idx = file_start_offset + current_entry_in_file + 1
                                if global_idx <= total_non_zero_entries
                                    row_indices[global_idx] = row_idx_val
                                    col_indices[global_idx] = original_file_idx
                                    values_vec[global_idx] = ValType(count_val)
                                    current_entry_in_file += 1
                                else
                                    Base.@error "Internal error: Exceeded sparse matrix capacity (total_non_zero_entries = $total_non_zero_entries) for file $original_file_idx. Global Idx: $global_idx. Kmer: $kmer."
                                    break
                                end
                            end
                        end
                    end
                    if current_entry_in_file != nnz_per_file[original_file_idx]
                        Base.@warn "Mismatch in NNZ for file $original_file_idx ($(fasta_list[original_file_idx])) during Pass 3. Expected $(nnz_per_file[original_file_idx]), wrote $current_entry_in_file."
                    end
                else
                    Base.@warn "Temp file $temp_filename (for original file $original_file_idx) not found in Pass 3, though nnz_per_file was $(nnz_per_file[original_file_idx])."
                end
            catch e
                Base.@error "Error reading or processing temporary file $temp_filename (for original file $original_file_idx) in Pass 3: $e."
            end
        end
        Base.lock(progress_pass3_lock) do
            ProgressMeter.next!(progress_pass3)
        end
    end
    ProgressMeter.finish!(progress_pass3)
    Base.@info "Pass 3 finished."

    Base.@info "Constructing sparse matrix ($num_kmers rows, $num_files columns, $total_non_zero_entries non-zero entries)..."

    kmer_counts_sparse_matrix = if total_non_zero_entries > 0 && num_kmers > 0
        SparseArrays.sparse(
            row_indices[1:total_non_zero_entries],
            col_indices[1:total_non_zero_entries],
            values_vec[1:total_non_zero_entries],
            num_kmers,
            num_files
        )
    else
        SparseArrays.spzeros(ValType, Int, num_kmers, num_files)
    end

    # --- After constructing the sparse matrix ---
    Base.@info "Done. Returning sorted kmer list, sparse counts matrix, and rarefaction data path."
    # Estimate memory needed for dense matrix
    dense_matrix_bytes_needed = Mycelia.estimate_dense_matrix_memory(ValType, num_kmers, num_files)

    # Compare with sparse matrix memory usage
    if dense_matrix_bytes_needed < sparse_matrix_bytes_needed
        @warn "A dense matrix would be more efficient"
        @warn "Estimated dense matrix memory usage: $(bytes_human_readable(dense_matrix_bytes_needed))"
        @warn "Estimated sparse matrix memory usage: $(bytes_human_readable(sparse_matrix_bytes_needed))"
    end

    final_result = (; kmers=sorted_kmers, counts=kmer_counts_sparse_matrix, rarefaction_data_path=rarefaction_data_path)

    try
        Base.Filesystem.rm(temp_dir; recursive=true, force=true)
        Base.@info "Successfully removed temporary directory: $temp_dir"
    catch e
        Base.@warn "Could not remove temporary directory $temp_dir: $e"
    end

    # Save result if requested
    if !isnothing(result_file)
        try
            JLD2.save_object(result_file, final_result)
            Base.@info "Saved kmer count results to $result_file"
        catch e
            Base.@error "Failed to save kmer count results to $result_file: $e"
        end
    end

    return final_result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the maximum number of possible canonical k-mers for a given alphabet.

# Arguments
- `k::Integer`: Length of k-mer
- `ALPHABET::Vector{Char}`: Character set (nucleotides or amino acids)

# Returns
- `Int`: Maximum number of possible canonical k-mers

# Details
- For amino acids (AA_ALPHABET): returns total possible k-mers
- For nucleotides: returns half of total possible k-mers (canonical form)
- Requires odd k-mer length for nucleotide alphabets
"""
function determine_max_canonical_kmers(k, ALPHABET)
    max_possible_kmers = determine_max_possible_kmers(k, ALPHABET)
    if ALPHABET == AA_ALPHABET
        return max_possible_kmers
    else
        @assert isodd(k) "this calculation is not valid for even length kmers"
        return Int(max_possible_kmers / 2)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the total number of possible unique k-mers that can be generated from a given alphabet.

# Arguments
- `k`: Length of k-mers to consider
- `ALPHABET`: Vector containing the allowed characters/symbols

# Returns
- Integer representing the maximum number of possible unique k-mers (|Σ|ᵏ)
"""
function determine_max_possible_kmers(k, ALPHABET)
    return length(ALPHABET)^k
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Canonicalizes the k-mer counts in the given dictionary.

This function iterates over the provided dictionary `kmer_counts`, which maps k-mers to their respective counts. For each k-mer that is not in its canonical form, it converts the k-mer to its canonical form and updates the count in the dictionary accordingly. If the canonical form of the k-mer already exists in the dictionary, their counts are summed. The original non-canonical k-mer is then removed from the dictionary.

# Arguments
- `kmer_counts::Dict{BioSequences.Kmer, Int}`: A dictionary where keys are k-mers and values are their counts.

# Returns
- The input dictionary `kmer_counts` with all k-mers in their canonical form, sorted by k-mers.
"""
function canonicalize_kmer_counts!(kmer_counts)
    # Only canonicalize nucleic acid k-mers (DNA/RNA), not amino acids
    if !isempty(kmer_counts)
        first_kmer = first(keys(kmer_counts))
        kmer_type = typeof(first_kmer)
        
        # Only canonicalize if it's a nucleic acid k-mer (DNA/RNA)
        # Check if the type parameters indicate a nucleic acid alphabet
        if kmer_type <: Kmers.Kmer{<:BioSequences.NucleicAcidAlphabet}
            for (kmer, count) in kmer_counts
                if !BioSequences.iscanonical(kmer)
                    canonical_kmer = BioSequences.canonical(kmer)
                    if haskey(kmer_counts, canonical_kmer)
                        kmer_counts[canonical_kmer] += count
                    else
                        kmer_counts[canonical_kmer] = count
                    end
                    delete!(kmer_counts, kmer)
                end
            end
        end
        # For amino acids, no canonicalization needed - they are already canonical
    end
    return sort!(kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalize k-mer counts into a canonical form by creating a non-mutating copy.

# Arguments
- `kmer_counts`: Dictionary or collection of k-mer count data

# Returns
- A new normalized k-mer count collection
"""
function canonicalize_kmer_counts(kmer_counts)
    return canonicalize_kmer_counts!(copy(kmer_counts))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count canonical k-mers in biological sequences. A canonical k-mer is the lexicographically 
smaller of a DNA sequence and its reverse complement, ensuring strand-independent counting.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer size and structure
- `sequences`: Iterator of biological sequences to analyze

# Returns
- `Dict{KMER_TYPE,Int}`: Dictionary mapping canonical k-mers to their counts
"""
function count_canonical_kmers(::Type{KMER_TYPE}, sequences) where KMER_TYPE
    kmer_counts = count_kmers(KMER_TYPE, sequences)
    return canonicalize_kmer_counts!(kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the frequency of each k-mer in a DNA sequence.

# Arguments
- `::Type{Kmers.Kmer{A,K}}`: K-mer type with alphabet A and length K
- `sequence::BioSequences.LongSequence`: Input DNA sequence to analyze

# Returns
A sorted dictionary mapping each k-mer to its frequency count in the sequence.

# Type Parameters
- `A <: BioSequences.DNAAlphabet`: DNA alphabet type
- `K`: Length of k-mers
"""
function count_kmers(::Type{Kmers.Kmer{A, K}}, sequence::BioSequences.LongSequence) where {A <: BioSequences.DNAAlphabet, K}
    # return sort(StatsBase.countmap(Kmers.FwDNAMers{K}(sequence)))
    return sort(StatsBase.countmap([kmer for (kmer, index) in Kmers.UnambiguousDNAMers{K}(sequence)]))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the frequency of each k-mer in an RNA sequence.

# Arguments
- `Kmer`: Type parameter specifying the k-mer length K and RNA alphabet
- `sequence`: Input RNA sequence to analyze

# Returns
- `Dict{Kmers.Kmer, Int}`: Sorted dictionary mapping each k-mer to its frequency count
"""
function count_kmers(::Type{Kmers.Kmer{A, K}}, sequence::BioSequences.LongSequence) where {A <: BioSequences.RNAAlphabet, K}
    # return sort(StatsBase.countmap(Kmers.FwRNAMers{K}(sequence)))
    return sort(StatsBase.countmap([kmer for (kmer, index) in Kmers.UnambiguousRNAMers{K}(sequence)]))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the frequency of amino acid k-mers in a biological sequence.

# Arguments
- `Kmers.Kmer{A,K}`: Type parameter specifying amino acid alphabet (A) and k-mer length (K)
- `sequence`: Input biological sequence to analyze

# Returns
A sorted dictionary mapping each k-mer to its frequency count in the sequence.
"""
function count_kmers(::Type{Kmers.Kmer{A, K}}, sequence::BioSequences.LongSequence) where {A <: BioSequences.AminoAcidAlphabet, K}
    return sort(StatsBase.countmap(Kmers.FwAAMers{K}(sequence)))
end

# TODO add a way to handle ambiguity or not
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the frequency of amino acid k-mers in a biological sequence.

# Arguments
- `Kmers.Kmer{A,K}`: Type parameter specifying amino acid alphabet (A) and k-mer length (K)
- `sequence`: Input biological sequence to analyze

# Returns
A sorted dictionary mapping each k-mer to its frequency count in the sequence.
"""
function count_kmers(::Type{KMER_TYPE}, record::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    # TODO: need to figure out how to infer the sequence type
    if eltype(KMER_TYPE) == BioSymbols.DNA
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongDNA{4}, record))
    elseif eltype(KMER_TYPE) == BioSymbols.RNA
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongRNA{4}, record))
    elseif eltype(KMER_TYPE) == BioSymbols.AminoAcid
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongAA, record))
    else
        @error KMER_TYPE
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers across multiple sequence records and return a sorted frequency table.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer length (e.g., `DNAKmer{3}` for 3-mers)
- `records`: Vector of FASTA/FASTQ records to analyze

# Returns
- `Dict{KMER_TYPE, Int}`: Sorted dictionary mapping k-mers to their frequencies
"""
function count_kmers(::Type{KMER_TYPE}, records::AbstractVector{T}) where {KMER_TYPE, T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    kmer_counts = count_kmers(KMER_TYPE, first(records))
    for record in records[2:end]
        _kmer_counts = count_kmers(KMER_TYPE, record)
        merge!(+, kmer_counts, _kmer_counts)
    end
    sort!(kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Counts k-mer occurrences in biological sequences from a FASTA/FASTQ reader.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer length and encoding (e.g., `DNAKmer{4}` for 4-mers)
- `sequences`: A FASTA or FASTQ reader containing the biological sequences to analyze

# Returns
A dictionary mapping k-mers to their counts in the input sequences
"""
function count_kmers(::Type{KMER_TYPE}, sequences::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Reader, FASTX.FASTQ.Reader}}
    return count_kmers(KMER_TYPE, collect(sequences))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers across multiple FASTA/FASTQ files and merge the results.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer length (e.g., `DNAKmer{4}` for 4-mers)
- `fastx_files`: Vector of paths to FASTA/FASTQ files

# Returns
- `Dict{KMER_TYPE, Int}`: Dictionary mapping k-mers to their total counts across all files
"""
function count_kmers(::Type{KMER_TYPE}, fastx_files::AbstractVector{T}) where {KMER_TYPE, T <: AbstractString}
    kmer_counts = count_kmers(KMER_TYPE, first(fastx_files))
    for file in fastx_files[2:end]
        _kmer_counts = count_kmers(KMER_TYPE, file)
        kmer_counts = merge!(+, kmer_counts, _kmer_counts)
    end
    return kmer_counts
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers in a FASTA/FASTQ file and return their frequencies.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer type (e.g., `DNAKmer{K}`)
- `fastx_file`: Path to input FASTA/FASTQ file

# Returns
- `Dict{KMER_TYPE, Int}`: Dictionary mapping each k-mer to its frequency
"""
function count_kmers(::Type{KMER_TYPE}, fastx_file::AbstractString) where {KMER_TYPE}
    fastx_io = open_fastx(fastx_file)
    kmer_counts = count_kmers(KMER_TYPE, fastx_io)
    close(fastx_io)
    return kmer_counts
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a dense k-mer counts table for a set of FASTA files, with disk-backed temporary storage, 
custom element type, robust error handling, and optional output file caching.
"""
function fasta_list_to_dense_kmer_counts(;
    fasta_list::AbstractVector{<:AbstractString},
    k::Integer,
    alphabet::Symbol,
    temp_dir_parent::AbstractString = Base.tempdir(),
    count_element_type::Union{Type{<:Unsigned}, Nothing} = nothing,
    result_file::Union{Nothing, AbstractString} = nothing,
    force::Bool = false,
    cleanup_temp::Bool = true
)

    num_files = length(fasta_list)
    if num_files == 0
        error("Input fasta_list is empty.")
    end
    if !(alphabet in (:AA, :DNA, :RNA))
        error("alphabet must be :AA, :DNA, or :RNA")
    end

    if alphabet == :AA
        if !(isa(k, Integer) && 0 < k <= 5)
            if force
                @warn "k should be a positive integer <= 5 for alphabet :AA.  Sparse counts are recommended for larger k."
            else
                error("k must be a positive integer <= 5 for alphabet :AA. Use sparse counts for larger k")
            end
        end
    else
        if !(isa(k, Integer) && 0 < k <= 11)
            if force
                @warn "k should be a positive integer <= 11 for alphabet :DNA or :RNA. Sparse counts are recommended for larger k."
            else
                error("k must be a positive integer <= 11 for alphabet :DNA or :RNA. Use sparse counts for larger k")
            end
        end
    end
    if any(f -> !Base.Filesystem.isfile(f), fasta_list)
        missing_files = [f for f in fasta_list if !Base.Filesystem.isfile(f)]
        error("Missing FASTA files: $(join(missing_files, ", ")).")
    end

    # Output file short-circuit
    if !isnothing(result_file) && Base.Filesystem.isfile(result_file) && !force
        Base.@info "result_file exists at $result_file. Loading and returning."
        return JLD2.load_object(result_file)
    end

    # KMER TYPE/COUNT/GENERATOR setup
    if alphabet == :AA
        KMER_TYPE = Kmers.AAKmer{k}
        ALPHABET_SEQ = Mycelia.AA_ALPHABET
        COUNT_FUNCTION = Mycelia.count_kmers
        GENERATE_KMERS = Mycelia.generate_all_possible_kmers
        sorted_kmers = sort(GENERATE_KMERS(k, ALPHABET_SEQ))
    elseif alphabet == :DNA
        KMER_TYPE = Kmers.DNAKmer{k}
        ALPHABET_SEQ = Mycelia.DNA_ALPHABET
        COUNT_FUNCTION = Mycelia.count_canonical_kmers
        GENERATE_KMERS = Mycelia.generate_all_possible_canonical_kmers
        sorted_kmers = sort(GENERATE_KMERS(k, ALPHABET_SEQ))
    elseif alphabet == :RNA
        KMER_TYPE = Kmers.RNAKmer{k}
        ALPHABET_SEQ = Mycelia.RNA_ALPHABET
        COUNT_FUNCTION = Mycelia.count_kmers
        GENERATE_KMERS = Mycelia.generate_all_possible_kmers
        sorted_kmers = sort(GENERATE_KMERS(k, ALPHABET_SEQ))
    end

    num_kmers = length(sorted_kmers)
    Base.@info "Counting $num_kmers $(KMER_TYPE) across $num_files files..."

    # Temp files
    temp_dir = Base.Filesystem.mktempdir(temp_dir_parent; prefix="dense_kmer_counts_")
    temp_file_paths = [Base.Filesystem.joinpath(temp_dir, "counts_$(i).jld2") for i in 1:num_files]

    # Pass 1: count kmers per file, record max, save to temp
    progress1 = ProgressMeter.Progress(num_files; desc="Counting: ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:cyan)
    lock = Base.ReentrantLock()
    error_log = Vector{Tuple{Int, String}}()
    successful_indices = Vector{Int}()
    max_observed_count_ref = Ref{Int}(0)

    Threads.@threads for idx in 1:num_files
        fasta_file = fasta_list[idx]
        temp_file = temp_file_paths[idx]
        local_max = 0
        try
            kmer_counts = COUNT_FUNCTION(KMER_TYPE, fasta_file)
            JLD2.save_object(temp_file, kmer_counts)
            if !isempty(kmer_counts)
                local_max = maximum(Base.values(kmer_counts))
            end
            Base.lock(lock)
            try
                push!(successful_indices, idx)
                if local_max > max_observed_count_ref[]
                    max_observed_count_ref[] = local_max
                end
            finally
                Base.unlock(lock)
            end
        catch e
            Base.lock(lock)
            try
                push!(error_log, (idx, string(e)))
            finally
                Base.unlock(lock)
            end
            try JLD2.save_object(temp_file, Dict{KMER_TYPE, Int}()) catch end
        end
        Base.lock(lock)
        try
            ProgressMeter.next!(progress1)
        finally
            Base.unlock(lock)
        end
    end
    ProgressMeter.finish!(progress1)

    sorted_successful_indices = sort(successful_indices)
    successful_fasta_list = fasta_list[sorted_successful_indices]

    num_successful_files = length(successful_fasta_list)
    percent_successful = round((num_successful_files / num_files) * 100, digits=3)
    Base.@info "$num_successful_files of $num_files ($(percent_successful)%) counted successfully"

    max_observed_count = max_observed_count_ref[]

    # Determine element type for counts
    if isnothing(count_element_type)
        ValType = max_observed_count <= typemax(UInt8) ? UInt8 :
                  max_observed_count <= typemax(UInt16) ? UInt16 :
                  max_observed_count <= typemax(UInt32) ? UInt32 : UInt64
    else
        ValType = count_element_type
        if max_observed_count > typemax(ValType)
            Base.@warn "User-specified count_element_type $ValType may be too small for max observed count $max_observed_count"
        end
    end
    Base.@info "Using $ValType for kmer counts: Maximum count observed = $(max_observed_count)"

    dense_matrix_bytes_needed = Mycelia.estimate_dense_matrix_memory(ValType, num_kmers, num_successful_files)
    # Check if matrix will fit in memory
    mem_check = check_matrix_fits_in_memory(dense_matrix_bytes_needed; severity=:error)
    if !mem_check.will_fit_total
        error("Matrix will not fit in available memory. Required: $(bytes_human_readable(bytes_needed)), Available: $(bytes_human_readable(mem_check.free_memory))")
    end

    # Build kmer index for matrix rows
    kmer_index = Dict{KMER_TYPE,Int}(kmer => i for (i, kmer) in enumerate(sorted_kmers))
    kmer_counts_matrix = zeros(ValType, num_kmers, num_successful_files)

    # Pass 2: fill matrix (multi-threaded, only for successful files)
    progress2 = ProgressMeter.Progress(num_successful_files; desc="Filling matrix: ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:green)
    Threads.@threads for col in 1:num_successful_files
        orig_idx = sorted_successful_indices[col]
        try
            kmer_counts = JLD2.load_object(temp_file_paths[orig_idx])
            for (kmer, count) in kmer_counts
                row = kmer_index[kmer]
                kmer_counts_matrix[row, col] = ValType(count)
            end
        catch
            # Optionally log error
        end
        Base.lock(lock)
        try
            ProgressMeter.next!(progress2)
        finally
            Base.unlock(lock)
        end
    end
    ProgressMeter.finish!(progress2)

    if cleanup_temp
        try Base.Filesystem.rm(temp_dir; recursive=true, force=true) catch end
    end

    sparse_matrix_bytes_needed = Mycelia.estimate_sparse_matrix_memory(ValType, num_kmers, num_successful_files, nnz = count(kmer_counts_matrix .!= 0))
    if sparse_matrix_bytes_needed < dense_matrix_bytes_needed
        @warn "A sparse matrix would be more efficient"
        @warn "Estimated dense matrix memory usage: $(bytes_human_readable(dense_matrix_bytes_needed))"
        @warn "Estimated sparse matrix memory usage: $(bytes_human_readable(sparse_matrix_bytes_needed))"
    end

    result = (;kmers=sorted_kmers, counts=kmer_counts_matrix, successful_fasta_list=successful_fasta_list, error_log=error_log)

    if !isnothing(result_file)
        try JLD2.save_object(result_file, result) catch end
    end

    return result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a collection of biological sequences into a dense k-mer count matrix.

# Arguments
- `biosequences`: Collection of DNA, RNA, or amino acid sequences (BioSequence types)
- `k::Integer`: Length of k-mers to count (must be ≤ 13)

# Returns
Named tuple containing:
- `sorted_kmers`: Vector of all possible k-mers in sorted order
- `kmer_counts_matrix`: Dense matrix where rows are k-mers and columns are sequences

# Details
- For DNA sequences, counts canonical k-mers (both strands)
- For RNA and protein sequences, counts exact k-mers
- Uses parallel processing with threads
"""
function biosequences_to_dense_counts_table(;biosequences, k)
    k >= 11 && error("use sparse counts to count k >= 11")    
    if eltype(first(biosequences)) == BioSymbols.AminoAcid
        KMER_TYPE = BioSequences.AminoAcidAlphabet
        sorted_kmers = sort(generate_all_possible_kmers(k, AA_ALPHABET))
        COUNT = count_kmers
    elseif eltype(first(biosequences)) == BioSymbols.DNA
        KMER_TYPE = BioSequences.DNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_canonical_kmers(k, DNA_ALPHABET))
        COUNT = count_canonical_kmers
    elseif eltype(first(biosequences)) == BioSymbols.RNA
        KMER_TYPE = BioSequences.RNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_kmers(k, RNA_ALPHABET))
        COUNT = count_kmers
    else
        error("invalid alphabet, please choose from :AA, :DNA, :RNA")
    end
    kmer_counts_matrix = zeros(length(sorted_kmers), length(biosequences))
    progress = ProgressMeter.Progress(length(biosequences))
    reenrantlock = ReentrantLock()
    Threads.@threads for (entity_index, biosequence) in collect(enumerate(biosequences))
        # Acquire the lock before updating the progress
        lock(reenrantlock) do
            # Update the progress meter
            ProgressMeter.next!(progress)
        end
        entity_mer_counts = COUNT(Kmers.Kmer{KMER_TYPE, k}, biosequence)
        for (i, kmer) in enumerate(sorted_kmers)
            kmer_counts_matrix[i, entity_index] = get(entity_mer_counts, kmer, 0)
        end
    end
    return (;sorted_kmers, kmer_counts_matrix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a collection of biological sequences into a k-mer count matrix.

# Arguments
- `biosequences`: Vector of biological sequences (DNA, RNA, or Amino Acids)
- `k`: Length of k-mers to count

# Returns
Named tuple with:
- `sorted_kmers`: Vector of all unique k-mers found, lexicographically sorted
- `kmer_counts_matrix`: Sparse matrix where rows are k-mers and columns are sequences

# Details
- For DNA sequences, counts canonical k-mers (both strands)
- Uses parallel processing with Thread-safe progress tracking
- Memory efficient sparse matrix representation
- Supports DNA, RNA and Amino Acid sequences
"""
function biosequences_to_counts_table(;biosequences, k)
    if eltype(first(biosequences)) == BioSymbols.AminoAcid
        KMER_TYPE = Kmers.AAKmer{k}
        COUNT = count_kmers
    elseif eltype(first(biosequences)) == BioSymbols.DNA
        KMER_TYPE = Kmers.DNAKmer{k}
        COUNT = count_canonical_kmers
    elseif eltype(first(biosequences)) == BioSymbols.RNA
        KMER_TYPE = Kmers.RNAKmer{k}
        COUNT = count_kmers
    else
        error("invalid alphabet, please choose from :AA, :DNA, :RNA")
    end
    
    kmer_counts = Vector{OrderedCollections.OrderedDict{KMER_TYPE, Int}}(undef, length(biosequences))
    progress = ProgressMeter.Progress(length(biosequences))
    reenrantlock = ReentrantLock()
    Threads.@threads for i in eachindex(biosequences)
        lock(reenrantlock) do
            ProgressMeter.next!(progress)
        end
        kmer_counts[i] = COUNT(KMER_TYPE, biosequences[i])
    end
    sorted_kmers = sort(collect(reduce(union, keys.(kmer_counts))))
    kmer_counts_matrix = SparseArrays.spzeros(Int, length(sorted_kmers), length(biosequences))
    @info "populating sparse counts matrix..."
    for (col, biosequence) in enumerate(biosequences)
        for (row, kmer) in enumerate(sorted_kmers)
            kmer_counts_matrix[row, col] = get(kmer_counts[col], kmer, 0)
        end
    end
    return (;sorted_kmers, kmer_counts_matrix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a dictionary of k-mer counts to a fixed-length numeric vector based on a predefined mapping.

# Arguments
- `kmer_to_index_map`: Dictionary mapping k-mer sequences to their corresponding vector indices
- `kmer_counts`: Dictionary containing k-mer sequences and their occurrence counts

# Returns
- A vector where each position corresponds to a k-mer count, with zeros for absent k-mers
"""
function kmer_counts_dict_to_vector(kmer_to_index_map, kmer_counts)
    kmer_counts_vector = zeros(length(kmer_to_index_map))
    for (kmer, count) in kmer_counts
        kmer_index = kmer_to_index_map[kmer]
        kmer_counts_vector[kmer_index] = count
    end
    return kmer_counts_vector
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

Generate a sorted list of all possible k-mers for a given alphabet.

# Arguments
- `k::Integer`: Length of k-mers to generate
- `alphabet`: Collection of symbols (DNA, RNA, or amino acids) from BioSymbols

# Returns
- Sorted Vector of Kmers of the appropriate type (DNA, RNA, or amino acid)
"""
function generate_all_possible_kmers(k, alphabet)
    kmer_iterator = Iterators.product([alphabet for i in 1:k]...)
    kmer_vectors = collect.(vec(collect(kmer_iterator)))
    if eltype(alphabet) == BioSymbols.AminoAcid
        kmers = [Kmers.AAKmer{k}(BioSequences.LongAA(kv)) for kv in kmer_vectors]
    elseif eltype(alphabet) == BioSymbols.DNA
        kmers = [Kmers.DNAKmer{k}(BioSequences.LongDNA{2}(kv)) for kv in kmer_vectors]
    elseif eltype(alphabet) == BioSymbols.RNA
        kmers = [Kmers.RNAKmer{k}(BioSequences.LongRNA{2}(kv)) for kv in kmer_vectors]
    else
        error()
    end
    return sort!(kmers)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

Generate all possible canonical k-mers of length `k` from the given `alphabet`.

For DNA/RNA sequences, returns unique canonical k-mers where each k-mer is represented by
the lexicographically smaller of itself and its reverse complement.
For amino acid sequences, returns all possible k-mers without canonicalization.

# Arguments
- `k`: Length of k-mers to generate
- `alphabet`: Vector of BioSymbols (DNA, RNA or AminoAcid)

# Returns
- Vector of k-mers, canonicalized for DNA/RNA alphabets
"""
function generate_all_possible_canonical_kmers(k, alphabet)
    kmers = generate_all_possible_kmers(k, alphabet)
    if eltype(alphabet) == BioSymbols.AminoAcid
        return kmers
    elseif eltype(alphabet) in (BioSymbols.DNA, BioSymbols.RNA)
        return unique!(BioSequences.canonical.(kmers))
    else
        error()
    end
end

# CSV is too memory inefficient, the others too slow :(
# # using uCSV
# # k=11
# # 3.444974 seconds (24.58 M allocations: 1.374 GiB, 34.65% gc time, 16.90% compilation time)
# # k=13
# # 362.285866 seconds (357.11 M allocations: 20.550 GiB, 91.60% gc time)

# # using DelimitedFiles.readdlm
# # k=11
# # 2.386620 seconds (16.11 M allocations: 632.732 MiB, 34.16% gc time, 24.25% compilation time)
# # k=13
# # 82.888552 seconds (227.49 M allocations: 8.766 GiB, 82.01% gc time)

# # CSV
# # k=11
# # 12.328422 seconds (7.62 M allocations: 732.639 MiB, 19091.67% compilation time: <1% of which was recompilation)
# # k=13
# # 37.098948 seconds (89.38 k allocations: 2.354 GiB, 93.56% gc time)

# function parse_jellyfish_counts(tabular_counts)
#     # load in the data
#     @assert occursin(r"\.gz$", tabular_counts) "this expects gzipped jellyfish tabular counts"
#     io = CodecZlib.GzipDecompressorStream(open(tabular_counts))
#     canonical_kmer_counts_table = DataFrames.DataFrame(CSV.File(io; delim='\t', header=false))
#     DataFrames.rename!(canonical_kmer_counts_table, [:Column1 => :kmer, :Column2 => :count])
    
#     # recode the kmers from strings to fixed sized kmer types
#     unique_kmer_lengths = unique(length.(canonical_kmer_counts_table[!, "kmer"]))
#     @assert length(unique_kmer_lengths) == 1
#     k = first(unique_kmer_lengths)
#     canonical_kmer_counts_table[!, "kmer"] = Kmers.DNAKmer{k}.(canonical_kmer_counts_table[!, "kmer"])
    
#     return canonical_kmer_counts_table
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Find the closest prime number to the given integer `n`.

Returns the nearest prime number to `n`. If two prime numbers are equally distant 
from `n`, returns the smaller one.

# Arguments
- `n::Int`: The input integer to find the nearest prime for

# Returns
- `Int`: The closest prime number to `n`
"""
function nearest_prime(n::Int)
    if n < 2
        return 2
    end
    next_p = Primes.nextprime(n)
    prev_p = Primes.prevprime(n)
    if n - prev_p <= next_p - n
        return prev_p
    else
        return next_p
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a sequence of Fibonacci numbers strictly less than the input value.

# Arguments
- `n::Int`: Upper bound (exclusive) for the Fibonacci sequence

# Returns
- `Vector{Int}`: Array containing Fibonacci numbers less than n
"""
function fibonacci_numbers_less_than(n::Int)
    if n <= 0
        return []
    elseif n == 1
        return [0]
    else
        fib = [0, 1]
        next_fib = fib[end] + fib[end-1]
        while next_fib < n
            push!(fib, next_fib)
            next_fib = fib[end] + fib[end-1]
        end
        return fib
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generates a specialized sequence of prime numbers combining:
- Odd primes up to 23 (flip_point)
- Primes nearest to Fibonacci numbers above 23 up to max

# Arguments
- `min::Int=0`: Lower bound for the sequence
- `max::Int=10_000`: Upper bound for the sequence

# Returns
Vector of Int containing the specialized prime sequence
"""
function ks(;min=0, max=10_000)
    # flip from all odd primes to only nearest to fibonnaci primes
    flip_point = 23
    # skip 19 because it is so similar to 17
    results = vcat(
        filter(x -> x != 19, filter(isodd, Primes.primes(0, flip_point))),
        filter(x -> x > flip_point, nearest_prime.(fibonacci_numbers_less_than(max*10)))
    )
    return filter(x -> min <= x <= max, results)
end

# function observed_kmer_frequencies(seq::BioSequences.BioSequence{A}, k::Int) where A<:BioSequences.Alphabet
#     kmer_count_dict = BioSequences.kmercounts(seq, k)
#     total_kmers = sum(values(kmer_count_dict))
#     return Dict(kmer => count / total_kmers for (kmer, count) in kmer_count_dict)
# end