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

# conda install -c bioconda kmer-jellyfish
# count, bc, info, stats, histo, dump, merge, query, cite, mem, jf
# cap at 4 threads, 8Gb per thread by default - this should be plenty fast enough for base usage, but open it up for higher performance!
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers in a FASTA/FASTQ file using Jellyfish.

# Arguments
- `fastx::String`: Path to input FASTA/FASTQ file (can be gzipped)
- `k::Integer`: k-mer length
- `threads::Integer=get_default_threads()`: Number of threads to use
- `max_mem::Integer=Int(Sys.free_memory())`: Maximum memory in bytes (defaults to system free memory)
- `canonical::Bool=false`: Whether to count canonical k-mers (both strands combined)
- `outfile::String=auto`: Output filename (auto-generated based on input and parameters)
- `conda_check::Bool=true`: Whether to verify Jellyfish conda installation

# Returns
- `String`: Path to gzipped TSV file containing k-mer counts
"""
function jellyfish_count(;fastx, k, threads=get_default_threads(), max_mem=Int(Sys.free_memory()), canonical=false, outfile = ifelse(canonical, "$(fastx).k$(k).canonical.jf", "$(fastx).k$(k).jf"), conda_check=true)
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

# function assess_assembly_quality(;assembled_sequence::BioSequences.LongDNA{2}, fastq::String, k::Int)
#     assess_assembly_quality(assembled_sequence=assembled_sequence, fastq=fastq, ks=[k])
# end

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
    kmer_frequency_histogram(kmer_counts)

Build a histogram of k-mer coverage counts.

# Arguments
- `kmer_counts::AbstractDict`: Dictionary mapping k-mers to counts
- `kmer_counts::AbstractVector{<:Integer}`: Vector of k-mer counts

# Returns
- `Dict{Int,Int}`: Mapping from coverage -> number of k-mers with that coverage
"""
function kmer_frequency_histogram(kmer_counts::AbstractDict)
    return StatsBase.countmap(values(kmer_counts))
end

function kmer_frequency_histogram(kmer_counts::AbstractVector{<:Integer})
    return StatsBase.countmap(kmer_counts)
end

"""
    coverage_peak_from_hist(kmer_hist; min_coverage=2)

Find the coverage peak in a k-mer frequency histogram.

# Arguments
- `kmer_hist`: Histogram mapping coverage to number of k-mers
- `min_coverage::Int=2`: Minimum coverage to consider for the peak

# Returns
- `NamedTuple`: `(coverage, kmers)` where `coverage` may be `missing` if no peak is found
"""
function coverage_peak_from_hist(kmer_hist; min_coverage::Int=2)
    peak_coverage = missing
    peak_kmers = 0

    for (coverage, kmer_count) in kmer_hist
        if coverage >= min_coverage && kmer_count > peak_kmers
            peak_coverage = coverage
            peak_kmers = kmer_count
        end
    end

    return (coverage=peak_coverage, kmers=peak_kmers)
end

"""
    analyze_kmer_frequency_spectrum(kmer_counts; min_coverage=2)

Analyze a k-mer frequency spectrum and return histogram + peak info.

# Arguments
- `kmer_counts`: Dictionary or vector of k-mer counts
- `min_coverage::Int=2`: Minimum coverage to consider for the peak

# Returns
- `NamedTuple`: `(histogram, peak, total_kmers, unique_kmers)`
"""
function analyze_kmer_frequency_spectrum(kmer_counts; min_coverage::Int=2)
    histogram = kmer_frequency_histogram(kmer_counts)
    peak = coverage_peak_from_hist(histogram; min_coverage=min_coverage)
    total_kmers = if kmer_counts isa AbstractDict
        sum(values(kmer_counts))
    else
        sum(kmer_counts)
    end
    unique_kmers = if kmer_counts isa AbstractDict
        length(kmer_counts)
    else
        length(kmer_counts)
    end

    return (;histogram, peak, total_kmers, unique_kmers)
end

function _sequence_kmer_iterator(sequence, k::Int)
    if sequence isa BioSequences.LongDNA || sequence isa Kmers.DNAKmer
        return Kmers.FwDNAMers{k}(sequence)
    elseif sequence isa BioSequences.LongRNA || sequence isa Kmers.RNAKmer
        return Kmers.FwRNAMers{k}(sequence)
    elseif sequence isa BioSequences.LongAA || sequence isa Kmers.AAKmer
        return Kmers.FwAAMers{k}(sequence)
    else
        error("Unsupported sequence type for k-mer iterator: $(typeof(sequence))")
    end
end

function _collect_kmers(sequence, k::Int)
    return collect(_sequence_kmer_iterator(sequence, k))
end

function _is_nucleotide_sequence(sequence)
    return sequence isa BioSequences.LongDNA ||
           sequence isa BioSequences.LongRNA ||
           sequence isa Kmers.DNAKmer ||
           sequence isa Kmers.RNAKmer
end

"""
    canonical_minimizers(sequence, k::Int, window::Int)

Compute canonical minimizers across a sliding window.

# Arguments
- `sequence`: DNA/RNA/AA sequence or k-mer
- `k::Int`: K-mer size
- `window::Int`: Window size in k-mers

# Returns
- `Vector`: Minimizer k-mers
"""
function canonical_minimizers(sequence, k::Int, window::Int)
    if k < 1
        error("k must be positive, got k=$k")
    end
    if window < 1
        error("window must be positive, got window=$window")
    end

    kmers = _collect_kmers(sequence, k)
    if length(kmers) < window
        return Vector{eltype(kmers)}()
    end

    minimizers = Vector{eltype(kmers)}()
    for i in 1:(length(kmers) - window + 1)
        window_kmers = kmers[i:i + window - 1]
        if _is_nucleotide_sequence(sequence)
            window_kmers = BioSequences.canonical.(window_kmers)
        end
        push!(minimizers, minimum(window_kmers))
    end

    return minimizers
end

"""
    open_syncmers(sequence, k::Int, s::Int, t::Int; canonical::Bool=false)

Compute open syncmers using the minimum s-mer position rule.

# Arguments
- `sequence`: DNA/RNA/AA sequence
- `k::Int`: K-mer size
- `s::Int`: S-mer size within each k-mer
- `t::Int`: 1-based position of the minimum s-mer to select
- `canonical::Bool=false`: Canonicalize s-mers for nucleotide sequences before comparison

# Returns
- `Vector`: Syncmer k-mers
"""
function open_syncmers(sequence, k::Int, s::Int, t::Int; canonical::Bool=false)
    if s < 1 || s > k
        error("s must be between 1 and k (k=$k, s=$s)")
    end
    max_position = k - s + 1
    if t < 1 || t > max_position
        error("t must be between 1 and $(max_position) (t=$t)")
    end

    kmers = _collect_kmers(sequence, k)
    syncmers = Vector{eltype(kmers)}()

    for kmer in kmers
        smers = _collect_kmers(kmer, s)
        if canonical && _is_nucleotide_sequence(kmer)
            smers = BioSequences.canonical.(smers)
        end
        min_smer = minimum(smers)
        min_pos = findfirst(==(min_smer), smers)
        if min_pos == t
            push!(syncmers, kmer)
        end
    end

    return syncmers
end

"""
    strobemers(sequence, k::Int, w_min::Int, w_max::Int; canonical::Bool=false)

Compute strobemers using a minstrobe selection strategy.

# Arguments
- `sequence`: DNA/RNA/AA sequence
- `k::Int`: K-mer size
- `w_min::Int`: Minimum window offset for the second strobe
- `w_max::Int`: Maximum window offset for the second strobe
- `canonical::Bool=false`: Canonicalize nucleotide k-mers before selection

# Returns
- `Vector{Tuple}`: Pairs of k-mers representing strobemers
"""
function strobemers(sequence, k::Int, w_min::Int, w_max::Int; canonical::Bool=false)
    if w_min < 1 || w_max < w_min
        error("w_min must be >= 1 and w_max >= w_min (w_min=$w_min, w_max=$w_max)")
    end

    kmers = _collect_kmers(sequence, k)
    strobes = Vector{Tuple{eltype(kmers), eltype(kmers)}}()

    for i in 1:length(kmers)
        j_start = i + w_min
        j_end = min(i + w_max, length(kmers))
        if j_start <= j_end
            candidates = kmers[j_start:j_end]
            if canonical && _is_nucleotide_sequence(sequence)
                candidates = BioSequences.canonical.(candidates)
            end
            second = minimum(candidates)
            first = canonical && _is_nucleotide_sequence(sequence) ? BioSequences.canonical(kmers[i]) : kmers[i]
            push!(strobes, (first, second))
        end
    end

    return strobes
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
        p = ProgressMeter.Progress(kmers_to_assess; dt=1)
    else
        p = ProgressMeter.Progress(kmers_to_assess; dt=1)
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
- `skip_rarefaction_plot::Bool`: If true, completely skip rarefaction plot generation (including loading plotting libraries). Defaults to `false`.
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
    skip_rarefaction_plot::Bool = false,
    result_file::Union{Nothing, AbstractString} = nothing,
    out_dir::Union{Nothing, AbstractString} = nothing,
    force::Bool = false,
    # Optimization control parameters
    force_threading::Union{Bool, Nothing} = nothing,
    force_temp_files::Union{Bool, Nothing} = nothing,
    force_progress_bars::Union{Bool, Nothing} = nothing,
    skip_rarefaction_data::Bool = false,
    # Heuristic tuning parameters
    threading_file_threshold::Int = 4,
    threading_size_threshold::Int = 10_000_000,
    memory_safety_factor::Float64 = 0.5,
    progress_time_threshold::Float64 = 5.0,
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

    # --- Smart Optimization Heuristics ---
    # Analyze file characteristics for optimization decisions  
    file_sizes = [Base.Filesystem.filesize(f) for f in fasta_list]
    total_file_size = sum(file_sizes)
    
    # Threading decision
    use_threading = if !isnothing(force_threading)
        force_threading
    else
        (num_files >= threading_file_threshold && total_file_size > threading_size_threshold) || 
        num_files >= 2 * threading_file_threshold
    end
    
    # Memory strategy decision
    estimated_avg_kmers_per_file = min(1000, total_file_size ÷ max(1, num_files) ÷ 4)
    estimated_total_kmer_dicts_memory = num_files * estimated_avg_kmers_per_file * 24
    use_temp_files = if !isnothing(force_temp_files)
        force_temp_files
    else
        estimated_total_kmer_dicts_memory > (1_000_000_000 * memory_safety_factor)
    end
    
    # Progress bars decision
    show_progress = if !isnothing(force_progress_bars)
        force_progress_bars
    else
        num_files > 20 || total_file_size > 50_000_000
    end
    
    # Rarefaction decision
    compute_rarefaction = !skip_rarefaction_data && (show_progress || !skip_rarefaction_plot)
    
    Base.@info "Sparse optimization decisions: threading=$use_threading, temp_files=$use_temp_files, progress=$show_progress, rarefaction=$compute_rarefaction ($(num_files) files, $(Base.format_bytes(total_file_size)))"

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
    rarefaction_points = compute_rarefaction ? Vector{Tuple{Int, Int}}() : Tuple{Int, Int}[]

    progress_pass2 = show_progress ? ProgressMeter.Progress(num_files; desc="Pass 2 (Aggregating Serially): ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:yellow) : nothing

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
        if compute_rarefaction
            push!(rarefaction_points, (i, length(all_kmers_set)))
        end
        if show_progress && !isnothing(progress_pass2)
            if compute_rarefaction
                ProgressMeter.next!(progress_pass2; showvalues = [(:unique_kmers, length(all_kmers_set))])
            else
                ProgressMeter.next!(progress_pass2)
            end
        end
    end
    if show_progress && !isnothing(progress_pass2)
        ProgressMeter.finish!(progress_pass2)
    end
    Base.@info "Pass 2 aggregation finished."

    # --- Process and Save/Plot Rarefaction Data (after Pass 2) ---
    rarefaction_data_path = _resolve_outpath(rarefaction_data_filename)
    
    if compute_rarefaction
        Base.@info "Saving rarefaction data to $rarefaction_data_path..."
    else
        Base.@info "Skipping rarefaction data collection and saving (skip_rarefaction_data=true)."
    end
    if compute_rarefaction
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
    end

    if !skip_rarefaction_plot
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
    else
        Base.@info "Skipping rarefaction plot generation (skip_rarefaction_plot=true)."
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
Fast in-memory implementation of dense k-mer counting for small datasets.
Skips temporary files and processes everything in memory for better performance.
"""
function _dense_kmer_counts_in_memory(
    fasta_list, sorted_kmers, KMER_TYPE, COUNT_FUNCTION, 
    use_threading, show_progress, count_element_type, result_file
)
    num_files = length(fasta_list)
    num_kmers = length(sorted_kmers)
    
    # Build kmer index for matrix rows
    kmer_index = Dict{KMER_TYPE,Int}(kmer => i for (i, kmer) in enumerate(sorted_kmers))
    
    # Process files and collect k-mer counts directly in memory
    all_kmer_counts = Vector{Dict{KMER_TYPE,Int}}(undef, num_files)
    max_observed_count = 0
    
    progress = show_progress ? ProgressMeter.Progress(num_files; desc="Counting (in-memory): ", color=:cyan) : nothing
    
    if use_threading
        Threads.@threads for idx in 1:num_files
            try
                kmer_counts = COUNT_FUNCTION(KMER_TYPE, fasta_list[idx])
                all_kmer_counts[idx] = kmer_counts
                if show_progress
                    ProgressMeter.next!(progress)
                end
            catch e
                Base.@error "Error processing file $(fasta_list[idx]): $e"
                all_kmer_counts[idx] = Dict{KMER_TYPE,Int}()
                if show_progress
                    ProgressMeter.next!(progress)
                end
            end
        end
    else
        for idx in 1:num_files
            try
                kmer_counts = COUNT_FUNCTION(KMER_TYPE, fasta_list[idx])
                all_kmer_counts[idx] = kmer_counts
                if show_progress
                    ProgressMeter.next!(progress)
                end
            catch e
                Base.@error "Error processing file $(fasta_list[idx]): $e"
                all_kmer_counts[idx] = Dict{KMER_TYPE,Int}()
                if show_progress
                    ProgressMeter.next!(progress)
                end
            end
        end
    end
    
    if show_progress && !isnothing(progress)
        ProgressMeter.finish!(progress)
    end
    
    # Find maximum count for data type selection
    for kmer_counts in all_kmer_counts
        if !isempty(kmer_counts)
            local_max = maximum(Base.values(kmer_counts))
            if local_max > max_observed_count
                max_observed_count = local_max
            end
        end
    end
    
    # Determine value type
    ValType = if isnothing(count_element_type)
        max_observed_count <= typemax(UInt8) ? UInt8 :
        max_observed_count <= typemax(UInt16) ? UInt16 :
        max_observed_count <= typemax(UInt32) ? UInt32 : UInt64
    else
        count_element_type
        if max_observed_count > typemax(count_element_type)
            Base.@warn "User-specified count_element_type $count_element_type may be too small for max observed count $max_observed_count"
        end
        count_element_type
    end
    
    Base.@info "Using $ValType for kmer counts: Maximum count observed = $max_observed_count"
    
    # Build dense matrix directly
    kmer_counts_matrix = zeros(ValType, num_kmers, num_files)
    
    progress2 = show_progress ? ProgressMeter.Progress(num_files; desc="Building matrix: ", color=:green) : nothing
    
    for (col, kmer_counts) in enumerate(all_kmer_counts)
        for (kmer, count) in kmer_counts
            row = kmer_index[kmer]
            kmer_counts_matrix[row, col] = ValType(count)
        end
        if show_progress && !isnothing(progress2)
            ProgressMeter.next!(progress2)
        end
    end
    
    if show_progress && !isnothing(progress2)
        ProgressMeter.finish!(progress2)
    end
    
    # Create result
    result = (kmers = sorted_kmers, counts = kmer_counts_matrix)
    
    # Save to file if requested
    if !isnothing(result_file)
        try
            JLD2.save_object(result_file, result)
            Base.@info "Results saved to $result_file"
        catch e
            Base.@warn "Failed to save results to $result_file: $e"
        end
    end
    
    return result
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
    cleanup_temp::Bool = true,
    # Optimization control parameters
    force_threading::Union{Bool, Nothing} = nothing,
    force_temp_files::Union{Bool, Nothing} = nothing,
    force_progress_bars::Union{Bool, Nothing} = nothing,
    # Heuristic tuning parameters
    threading_file_threshold::Int = 4,
    threading_size_threshold::Int = 10_000_000,
    memory_safety_factor::Float64 = 0.5,
    progress_time_threshold::Float64 = 5.0
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

    # --- Smart Optimization Heuristics ---
    # Analyze file characteristics for optimization decisions
    file_sizes = [Base.Filesystem.filesize(f) for f in fasta_list]
    total_file_size = sum(file_sizes)
    
    # Threading decision: use threading if we have enough files AND they're large enough,
    # OR if we have many files regardless of size
    use_threading = if !isnothing(force_threading)
        force_threading
    else
        (num_files >= threading_file_threshold && total_file_size > threading_size_threshold) || 
        num_files >= 2 * threading_file_threshold
    end
    
    # Estimate memory requirements for in-memory vs temp file strategy
    # This is a rough estimate - we'll do more precise calculation later
    estimated_avg_kmers_per_file = min(1000, total_file_size ÷ max(1, num_files) ÷ 4)  ## Rough estimate based on file size
    estimated_total_kmer_dicts_memory = num_files * estimated_avg_kmers_per_file * 24  ## Dict overhead estimate
    
    # Memory strategy decision
    use_temp_files = if !isnothing(force_temp_files)
        force_temp_files
    else
        # Use temp files if estimated memory usage exceeds safety factor
        # We'll do a more precise check later, this is just initial heuristic
        estimated_total_kmer_dicts_memory > (1_000_000_000 * memory_safety_factor)  ## 1GB base threshold
    end
    
    # Progress bars decision
    show_progress = if !isnothing(force_progress_bars)
        force_progress_bars
    else
        # Show progress for operations that might take a while
        num_files > 20 || total_file_size > 50_000_000  ## 50MB threshold
    end
    
    Base.@info "Optimization decisions: threading=$use_threading, temp_files=$use_temp_files, progress=$show_progress ($(num_files) files, $(Base.format_bytes(total_file_size)))"

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

    # Branch to optimized execution path based on heuristics
    if !use_temp_files
        return _dense_kmer_counts_in_memory(
            fasta_list, sorted_kmers, KMER_TYPE, COUNT_FUNCTION, 
            use_threading, show_progress, count_element_type, result_file
        )
    end

    # Temp files (original temp file-based implementation)
    temp_dir = Base.Filesystem.mktempdir(temp_dir_parent; prefix="dense_kmer_counts_")
    temp_file_paths = [Base.Filesystem.joinpath(temp_dir, "counts_$(i).jld2") for i in 1:num_files]

    # Pass 1: count kmers per file, record max, save to temp
    progress1 = show_progress ? ProgressMeter.Progress(num_files; desc="Counting: ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:cyan) : nothing
    lock = use_threading ? Base.ReentrantLock() : nothing
    error_log = Vector{Tuple{Int, String}}()
    successful_indices = Vector{Int}()
    max_observed_count_ref = Ref{Int}(0)

    # Define processing loop based on threading decision
    function process_files()
        if use_threading
            Threads.@threads for idx in 1:num_files
                process_single_file(idx)
            end
        else
            for idx in 1:num_files
                process_single_file(idx)
            end
        end
    end
    
    function process_single_file(idx)
        fasta_file = fasta_list[idx]
        temp_file = temp_file_paths[idx]
        local_max = 0
        try
            kmer_counts = COUNT_FUNCTION(KMER_TYPE, fasta_file)
            JLD2.save_object(temp_file, kmer_counts)
            if !isempty(kmer_counts)
                local_max = maximum(Base.values(kmer_counts))
            end
            if use_threading
                Base.lock(lock)
            end
            try
                push!(successful_indices, idx)
                if local_max > max_observed_count_ref[]
                    max_observed_count_ref[] = local_max
                end
            finally
                if use_threading
                    Base.unlock(lock)
                end
            end
        catch e
            if use_threading
                Base.lock(lock)
            end
            try
                push!(error_log, (idx, string(e)))
            finally
                if use_threading
                    Base.unlock(lock)
                end
            end
            try JLD2.save_object(temp_file, Dict{KMER_TYPE, Int}()) catch end
        end
        if show_progress
            if use_threading
                Base.lock(lock)
            end
            try
                ProgressMeter.next!(progress1)
            finally
                if use_threading
                    Base.unlock(lock)
                end
            end
        end
    end
    
    # Execute the file processing
    process_files()
    
    if show_progress && !isnothing(progress1)
        ProgressMeter.finish!(progress1)
    end

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

import Primes
import DocStringExtensions

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a √2-scaled ladder of odd primes for k-mer-based assembly / error-screening.

Starts with a user-supplied list of `seed_primes` (default **[3, 5, 7]**), then
iteratively multiplies the last accepted *k* by `ratio` (default `sqrt(2)`),
rounds **up** to the next odd prime, and appends it **only if** it differs from
the previous accepted prime by at least `min_fractional_gap`.

Keyword arguments
=================
- `max_k::Int = 10_000`          : Absolute upper bound.
- `seed_primes::Vector{Int}`     : Initial primes (e.g. `[3,5,7]` for protein,
                                   `[11,13,17]` for nucleotides).
- `ratio::Float64 = sqrt(2)`     : Target geometric growth factor.
- `min_fractional_gap::Float64 = 0.30` : Minimum (k_new − k_prev)/k_prev to skip
                                         “sister” primes.
- `read_length::Union{Int,Nothing} = nothing` : If set, cap at
                                               `read_length − read_margin`.
- `read_margin::Int = 20`        : Safety margin for short-read data.
- `only_odd::Bool = true`        : Force odd *k* (recommended).
- `return_unique::Bool = true`   : De-duplicate before returning.
- `min_k::Int = 3`               : Drop any *k* below this after generation.

Returns
=======
`Vector{Int}` — ascending prime *k* values suitable for `-k`/`--k-list`.
"""
function k_ladder(; max_k::Int           = 10_000,
                   seed_primes::Vector{Int} = [3, 5, 7],
                   ratio::Float64        = sqrt(2),
                   min_fractional_gap::Float64 = 0.30,
                   read_length::Union{Int,Nothing} = nothing,
                   read_margin::Int      = 20,
                   only_odd::Bool        = true,
                   return_unique::Bool   = true,
                   min_k::Int            = 3)

    # clean & sort the seeds, keep ≥ 2
    seeds = sort(unique(filter(>=(2), seed_primes)))
    seeds = only_odd ? filter(isodd, seeds) : seeds
    ks    = collect(seeds)

    # set effective ceiling
    k_ceiling = max_k
    if read_length !== nothing
        k_ceiling = min(k_ceiling, read_length - read_margin)
    end
    k_ceiling = max(k_ceiling, ks[end])  # never shrink below last seed

    # ladder growth
    lastk = ks[end]
    while true
        cand = ceil(Int, lastk * ratio)
        cand += (only_odd && iseven(cand)) ? 1 : 0       # force odd
        candp = Primes.nextprime(cand)                   # bump to next prime
        candp > k_ceiling && break                       # stop if too large
        if (candp - lastk) / lastk + 1e-12 >= min_fractional_gap
            push!(ks, candp)
            lastk = candp
        else
            cand   = candp + 2        # too close → move forward and retry
        end
    end

    ks = filter(>=(min_k), ks)         # drop anything below min_k
    return return_unique ? unique(ks) : ks
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Estimate genome size from k-mer analysis using total k-mer count.

This function estimates genome size using the basic relationship: 
genome_size ≈ total_kmers - k + 1, where total_kmers is the sum of all k-mer 
counts. This is a simple estimation method; more sophisticated approaches 
accounting for sequencing depth, repeats, and errors may be more accurate.

# Arguments
- `sequence::Union{BioSequences.LongSequence, AbstractString}`: Input sequence or string
- `k::Integer`: K-mer size for analysis

# Returns
- `Dict{String, Any}`: Dictionary containing:
  - "unique_kmers": Number of unique k-mers observed
  - "total_kmers": Total k-mer count (sum of all frequencies)
  - "estimated_genome_size": Estimated genome size
  - "actual_size": Length of input sequence (if provided)

# Examples
```julia
# Estimate genome size from a sequence
sequence = "ATCGATCGATCGATCG"
result = estimate_genome_size_from_kmers(sequence, 5)
```
"""
function estimate_genome_size_from_kmers(sequence::Union{BioSequences.LongSequence, AbstractString}, k::Integer)
    # Convert string to BioSequence if needed
    if isa(sequence, AbstractString)
        bio_sequence = BioSequences.LongDNA{4}(sequence)
    else
        bio_sequence = sequence
    end
    
    # Count k-mers using existing infrastructure
    # Use direct sequence type checks - much more robust than type parameter extraction
    if bio_sequence isa BioSequences.LongDNA
        kmer_counts = count_kmers(Kmers.DNAKmer{k}, bio_sequence)
    elseif bio_sequence isa BioSequences.LongRNA  
        kmer_counts = count_kmers(Kmers.RNAKmer{k}, bio_sequence)
    elseif bio_sequence isa BioSequences.LongAA
        kmer_counts = count_kmers(Kmers.AAKmer{k}, bio_sequence)
    else
        error("Unsupported sequence type: $(typeof(bio_sequence))")
    end
    
    unique_kmers = length(kmer_counts)
    total_kmers = sum(values(kmer_counts))
    
    # Basic genome size estimation: total_kmers - k + 1 ≈ genome_size
    estimated_size = total_kmers - k + 1
    
    return Dict(
        "unique_kmers" => unique_kmers,
        "total_kmers" => total_kmers,
        "estimated_genome_size" => estimated_size,
        "actual_size" => length(bio_sequence)
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Estimate genome size from FASTQ/FASTA records using k-mer analysis.

Overload for processing FASTQ or FASTA records directly.

# Arguments  
- `records::AbstractVector{T}` where T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}`: Input records
- `k::Integer`: K-mer size for analysis

# Returns
- `Dict{String, Any}`: Dictionary with k-mer statistics and genome size estimate
"""
function estimate_genome_size_from_kmers(records::AbstractVector{T}, k::Integer) where {T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    # Determine sequence type from first record using existing robust functions
    first_seq_string = String(FASTX.sequence(records[1]))
    first_seq_biosequence = convert_sequence(first_seq_string)
    
    # Use direct sequence type checks on the converted BioSequence
    if first_seq_biosequence isa BioSequences.LongDNA
        kmer_type = Kmers.DNAKmer{k}
    elseif first_seq_biosequence isa BioSequences.LongRNA
        kmer_type = Kmers.RNAKmer{k}  
    elseif first_seq_biosequence isa BioSequences.LongAA
        kmer_type = Kmers.AAKmer{k}
    else
        error("Unsupported sequence type: $(typeof(first_seq_biosequence))")
    end
    
    # Count k-mers across all records
    kmer_counts = count_kmers(kmer_type, records)
    
    unique_kmers = length(kmer_counts)
    total_kmers = sum(values(kmer_counts))
    
    # Calculate total sequence length
    total_length = sum(length(FASTX.sequence(record)) for record in records)
    
    # Basic genome size estimation
    estimated_size = total_kmers - k + 1
    
    return Dict(
        "unique_kmers" => unique_kmers,
        "total_kmers" => total_kmers,
        "estimated_genome_size" => estimated_size,
        "actual_size" => total_length,
        "num_records" => length(records)
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform simultaneous multi-scale k-mer analysis with prime k-mer sizes.

# Arguments
- `sequences::Vector{BioSequences.LongDNA{4}}`: Vector of DNA sequences to analyze
- `prime_ks::Vector{Int}`: Prime k-mer sizes to analyze simultaneously (default: [3, 5, 7, 11, 13, 17, 19])
- `window_size::Int`: Sliding window size for quality averaging (default: 50)
- `step_size::Int`: Step size for sliding window (default: 10)

# Returns
Named tuple containing:
- `multi_k_counts::Dict{Int, Dict}`: K-mer counts for each k value
- `quality_profiles::Dict{Int, Vector{Float64}}`: Quality profiles for each k
- `consensus_kmers::Dict{Int, Vector}`: Consensus k-mers with highest confidence
- `coverage_estimates::Dict{Int, Float64}`: Coverage estimates for each k

# Details
- Implements universal biological polymer assembly algorithm design
- Uses simultaneous analysis of multiple prime k-mer lengths
- Incorporates sliding window instantaneous quality averaging
- Supports adaptive k-mer selection based on coverage patterns
"""
function multi_scale_kmer_analysis(sequences::Vector{BioSequences.LongDNA{4}}; prime_ks::Vector{Int}=[3, 5, 7, 11, 13, 17, 19], window_size::Int=50, step_size::Int=10)
    multi_k_counts = Dict{Int, Dict}()
    quality_profiles = Dict{Int, Vector{Float64}}()
    consensus_kmers = Dict{Int, Vector}()
    coverage_estimates = Dict{Int, Float64}()
    
    for k in prime_ks
        @info "Analyzing k-mer size: $k"
        
        # Count k-mers for this k value
        kmer_counts = Dict{BioSequences.LongDNA{4}, Int}()
        total_kmers = 0
        
        for seq in sequences
            if length(seq) >= k
                # Extract k-mers with sliding window
                for i in 1:(length(seq) - k + 1)
                    kmer = seq[i:(i + k - 1)]
                    kmer_counts[kmer] = get(kmer_counts, kmer, 0) + 1
                    total_kmers += 1
                end
            end
        end
        
        multi_k_counts[k] = kmer_counts
        
        # Calculate quality profile using sliding window
        quality_profile = Float64[]
        for seq in sequences
            if length(seq) >= window_size
                for start in 1:step_size:(length(seq) - window_size + 1)
                    window_end = min(start + window_size - 1, length(seq))
                    window_seq = seq[start:window_end]
                    
                    # Calculate quality score for this window (simplified)
                    window_quality = calculate_window_quality(window_seq, k)
                    push!(quality_profile, window_quality)
                end
            end
        end
        
        quality_profiles[k] = quality_profile
        
        # Select consensus k-mers (high-frequency, high-quality)
        sorted_kmers = sort(collect(kmer_counts), by=x->x[2], rev=true)
        consensus_threshold = max(3, Int(ceil(0.1 * length(sorted_kmers))))
        consensus_kmers[k] = [kmer for (kmer, count) in sorted_kmers[1:consensus_threshold]]
        
        # Estimate coverage
        if !isempty(kmer_counts)
            coverage_estimates[k] = Statistics.mean(values(kmer_counts))
        else
            coverage_estimates[k] = 0.0
        end
    end
    
    return (;multi_k_counts, quality_profiles, consensus_kmers, coverage_estimates)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate quality score for a sequence window.

# Arguments
- `window_seq::BioSequences.LongDNA{4}`: Sequence window to analyze
- `k::Int`: K-mer size for analysis

# Returns
Float64 quality score for the window
"""
function calculate_window_quality(window_seq::BioSequences.LongDNA{4}, k::Int)
    if length(window_seq) < k
        return 0.0
    end
    
    # Count unique k-mers in window
    unique_kmers = Set{BioSequences.LongDNA{4}}()
    for i in 1:(length(window_seq) - k + 1)
        kmer = window_seq[i:(i + k - 1)]
        push!(unique_kmers, kmer)
    end
    
    # Quality score based on k-mer diversity and GC content
    diversity_score = length(unique_kmers) / max(1, length(window_seq) - k + 1)
    gc_content = BioSequences.gc_content(window_seq)
    gc_balance = 1.0 - abs(gc_content - 0.5) * 2  # Penalize extreme GC content
    
    return diversity_score * gc_balance
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform adaptive k-mer selection based on coverage patterns.

# Arguments
- `multi_k_results`: Results from multi_scale_kmer_analysis
- `target_coverage::Float64`: Target coverage threshold (default: 10.0)

# Returns
Vector{Int} of optimal k-mer sizes for the given data
"""
function adaptive_kmer_selection(multi_k_results; target_coverage::Float64=10.0)
    coverage_estimates = multi_k_results.coverage_estimates
    
    # Select k values with coverage closest to target
    optimal_ks = Int[]
    
    for (k, coverage) in sort(collect(coverage_estimates), by=x->abs(x[2] - target_coverage))
        if coverage >= target_coverage * 0.5  # At least 50% of target coverage
            push!(optimal_ks, k)
        end
        
        # Limit to top 3 k values
        if length(optimal_ks) >= 3
            break
        end
    end
    
    # Ensure we have at least one k value
    if isempty(optimal_ks) && !isempty(coverage_estimates)
        optimal_ks = [argmax(coverage_estimates)]
    end
    
    return sort(optimal_ks)
end

# function observed_kmer_frequencies(seq::BioSequences.BioSequence{A}, k::Int) where A<:BioSequences.Alphabet
#     kmer_count_dict = BioSequences.kmercounts(seq, k)
#     total_kmers = sum(values(kmer_count_dict))
#     return Dict(kmer => count / total_kmers for (kmer, count) in kmer_count_dict)
# end
