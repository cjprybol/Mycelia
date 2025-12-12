"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert SAM/BAM records from a XAM.SAM.Reader into a DataFrame.

Parameters:
- `reader`: A XAM.SAM.Reader object for iterating over records

Returns:
- A DataFrame containing all record data in a structured format
"""
function xam_to_dataframe(reader::XAM.SAM.Reader)::DataFrames.DataFrame
    return _xam_to_dataframe_common(reader, XAM.SAM)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert SAM/BAM records from a XAM.BAM.Reader into a DataFrame.

Parameters:
- `reader`: A XAM.BAM.Reader object for iterating over records

Returns:
- A DataFrame containing all record data in a structured format
"""
function xam_to_dataframe(reader::XAM.BAM.Reader)::DataFrames.DataFrame
    return _xam_to_dataframe_common(reader, XAM.BAM)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Common implementation for converting XAM records to DataFrame.
Handles both SAM and BAM readers using the appropriate module.

Parameters:
- `reader`: A XAM reader object (SAM or BAM)
- `MODULE`: XAM.SAM or XAM.BAM module to use for record parsing

Returns:
- A DataFrame containing all record data in a structured format
"""
function _xam_to_dataframe_common(reader, MODULE)::DataFrames.DataFrame
    # Define empty vectors for each column
    templates = String[]
    ismapped = Bool[]
    isprimary = Bool[]
    flags = UInt16[]
    references = Union{Missing, String}[]
    positions = UnitRange{Int}[]
    mappingqualities = UInt8[]
    cigars = String[]
    rnexts = Union{String, Missing}[]
    pnexts = Union{Int, Missing}[]
    tlens = Union{Int, Missing}[]
    seqs = Union{Nothing, BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}[]
    quals = Union{Missing, Vector{UInt8}}[]
    alignlengths = Union{Int, Missing}[]
    alignment_score = Union{Int, Missing}[]
    mismatches = Union{Int, Missing}[]
    # aux_tags = Vector{Dict{String,Any}}()
    
    # Process each record
    for record in reader
        try
            push!(templates, MODULE.tempname(record))
        catch err
            println("Error extracting tempname from record: ", record)
            rethrow(err)
        end

        try
            push!(ismapped, MODULE.ismapped(record))
        catch err
            println("Error extracting ismapped from record: ", record)
            rethrow(err)
        end

        try
            push!(isprimary, MODULE.isprimary(record))
        catch err
            println("Error extracting isprimary from record: ", record)
            rethrow(err)
        end

        try
            push!(flags, MODULE.flag(record))
        catch err
            println("Error extracting flag from record: ", record)
            rethrow(err)
        end

        try
            push!(references, MODULE.ismapped(record) ? MODULE.refname(record) : missing)
        catch err
            println("Error extracting refname from record: ", record)
            rethrow(err)
        end

        try
            push!(positions, MODULE.position(record):MODULE.rightposition(record))
        catch err
            println("Error extracting positions from record: ", record)
            rethrow(err)
        end

        try
            push!(mappingqualities, MODULE.mappingquality(record))
        catch err
            println("Error extracting mappingquality from record: ", record)
            rethrow(err)
        end

        try
            push!(cigars, MODULE.cigar(record))
        catch err
            println("Error extracting cigar from record: ", record)
            rethrow(err)
        end

        try
            push!(rnexts, MODULE.nextposition(record) == 0 ? missing : MODULE.nextrefname(record))
        catch err
            println("Error extracting rnext from record: ", record)
            rethrow(err)
        end

        try
            push!(pnexts, MODULE.nextposition(record))
        catch err
            println("Error extracting pnext from record: ", record)
            rethrow(err)
        end

        try
            push!(tlens, MODULE.templength(record))
        catch err
            println("Error extracting templength from record: ", record)
            rethrow(err)
        end

        try
            push!(seqs, MODULE.sequence(record))
        catch err
            println("Error extracting sequence from record: ", record)
            rethrow(err)
        end

        try
            push!(quals, MODULE.isprimary(record) ? MODULE.quality(record) : missing)
        catch err
            println("Error extracting quality from record: ", record)
            rethrow(err)
        end

        try
            push!(alignlengths, MODULE.alignlength(record))
        catch err
            println("Error extracting alignlength from record: ", record)
            rethrow(err)
        end

        try
            push!(alignment_score, MODULE.ismapped(record) ? record["AS"] : missing)
        catch err
            println("Error extracting alignment_score from record: ", record)
            rethrow(err)
        end

        try
            push!(mismatches,  MODULE.ismapped(record) ? record["NM"] : missing)
        catch err
            println("Error extracting mismatches from record: ", record)
            rethrow(err)
        end
        
        # not currently working?
        # # Extract auxiliary tags
        # tags = Dict{String,Any}()
        # for tag in XAM.SAM.auxdata(record)
        #     tag_name = String(tag.tag)
        #     tags[tag_name] = tag.value
        # end
        # push!(aux_tags, tags)
    end
    
    # Create and return the DataFrame
    return DataFrames.DataFrame(
        template = templates,
        ismapped = ismapped,
        isprimary = isprimary,
        flag = flags,
        reference = references,
        position = positions,
        mappingquality = mappingqualities,
        cigar = cigars,
        rnext = rnexts,
        pnext = pnexts,
        tlen = tlens,
        seq = seqs,
        qual = quals,
        alignlength = alignlengths,
        alignment_score = alignment_score,
        mismatches = mismatches,
        # aux_tags = aux_tags
    )
end

# """
#     record_to_row(record::XAM.SAM.Record)::NamedTuple

# Helper function to convert a single SAM record to a row for the DataFrame.

# Parameters:
# - `record`: A XAM.SAM.Record object

# Returns:
# - A NamedTuple containing the field values
# """
# function record_to_row(record::XAM.SAM.Record)::NamedTuple
#     # Extract auxiliary tags
#     tags = Dict{String,Any}()
#     for tag in XAM.SAM.auxdata(record)
#         tag_name = String(tag.tag)
#         tags[tag_name] = tag.value
#     end
    
#     # Return as a named tuple
#     return (
#         qname = XAM.SAM.tempname(record),
#         flag = XAM.SAM.flag(record),
#         rname = XAM.SAM.refname(record),
#         pos = XAM.SAM.position(record),
#         mapq = XAM.SAM.mappingquality(record),
#         cigar = XAM.SAM.cigar(record),
#         rnext = XAM.SAM.nextrefname(record),
#         pnext = XAM.SAM.nextposition(record),
#         tlen = XAM.SAM.templength(record),
#         seq = XAM.SAM.sequence(record),
#         qual = XAM.SAM.quality(record),
#         aux_tags = tags
#     )
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a SAM/BAM file to a DataFrame using the open_xam function.

Parameters:
- `xam_path`: Path to the SAM/BAM file
- `header`: Whether to include the header (default: false)

Returns:
- A DataFrame containing the parsed data
"""
function xam_to_dataframe(xam_path::String)::DataFrames.DataFrame
    reader = open_xam(xam_path)
    try
        return xam_to_dataframe(reader)
    finally
        close(reader)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Open a XAM (SAM/BAM/CRAM) file with the specified parser.

# Arguments
- `xam`: Path to the XAM file
- `header::Bool=false`: Whether to include header information
- `parser::Symbol=:auto`: Parser to use (:auto, :xamjl, :samtools)
  - `:auto`: Automatically choose parser based on file format detection
  - `:xamjl`: Use XAM.jl direct parsing (supports SAM, SAM.gz, BAM)
  - `:samtools`: Use samtools view (supports SAM, BAM, CRAM)

# Returns
- XAM.SAM.Reader or XAM.BAM.Reader object for reading records

# Details
- `:auto` detects file format by examining file headers and chooses appropriate parser
- XAM.jl parser is faster for direct file access and supports SAM, SAM.gz, and BAM
- samtools parser handles all formats including CRAM and provides additional validation
- samtools parser always returns XAM.SAM.Reader regardless of input format
"""
function open_xam(xam; header=false, parser=:auto)
    if !isfile(xam)
        error("File not found: ", xam)
    end
    if filesize(xam) == 0
        error("File is empty: ", xam)
    end

    # Determine parser to use - samtools is more robust and handles all formats
    if parser == :auto
        parser = :samtools
    end

    if parser == :xamjl
        return open_xam_xamjl(xam; header=header)
    elseif parser == :samtools
        return open_xam_samtools(xam; header=header)
    else
        error("Invalid parser: $parser. Use :auto, :xamjl, or :samtools")
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Detect the format of a XAM file by examining its header bytes.

# Arguments
- `xam`: Path to the XAM file

# Returns
- Symbol indicating format: :sam, :sam_gz, :bam, :cram, or :unknown
"""
function detect_xam_format(xam)
    # First check file extension for disambiguation of BGZF formats
    ext = lowercase(splitext(xam)[2])
    filename_lower = lowercase(xam)

    # Check for compressed SAM formats by extension first
    if endswith(filename_lower, ".sam.gz") || endswith(filename_lower, ".sam.bgz") || ext == ".bgz"
        return :sam_gz
    elseif ext == ".bam"
        return :bam
    elseif ext == ".cram"
        return :cram
    elseif ext == ".sam"
        return :sam
    end

    # If extension is ambiguous, check magic bytes
    open(xam, "r") do io
        # Read first few bytes to identify format
        header = read(io, min(8, filesize(xam)))

        if length(header) >= 4
            # CRAM format: starts with "CRAM"
            if length(header) >= 4 && String(header[1:4]) == "CRAM"
                return :cram
            # BGZF compressed formats: starts with specific magic bytes
            elseif header[1:4] == UInt8[0x1f, 0x8b, 0x08, 0x04]
                # BGZF format - could be BAM or compressed SAM
                # Without extension context, assume BAM (more common)
                return :bam
            # Standard gzip: starts with gzip magic bytes
            elseif header[1:2] == UInt8[0x1f, 0x8b]
                return :sam_gz
            # SAM format: typically starts with '@' (header) or read name
            elseif header[1] == UInt8('@')
                return :sam
            end
        end

        return :unknown
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Open a XAM file using XAM.jl direct parsing.
Supports SAM, SAM.gz, and BAM files.
Returns appropriate XAM.SAM.Reader or XAM.BAM.Reader.
"""
function open_xam_xamjl(xam; header=false)
    format = detect_xam_format(xam)

    try
        if format == :sam
            return XAM.SAM.Reader(open(xam))
        elseif format == :sam_gz
            ## For BGZF-compressed SAM files, use external samtools for reliable decompression
            ## since concatenated BGZF blocks from samtools are complex to parse directly
            try
                Mycelia.add_bioconda_env("samtools")
                # Create a command that pipes samtools view output to XAM.SAM.Reader
                cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -h $(xam)`
                return XAM.SAM.Reader(cmd)
            catch e
                ## If samtools approach fails, try CodecBGZF as backup
                file_io = open(xam, "r")
                try
                    bgzf_io = CodecBGZF.BGZFDecompressorStream(file_io)
                    return XAM.SAM.Reader(bgzf_io)
                catch e2
                    close(file_io)
                    ## Final fallback to standard gzip
                    file_io = open(xam, "r")
                    gz_io = CodecZlib.GzipDecompressorStream(file_io)
                    return XAM.SAM.Reader(gz_io)
                end
            end
        elseif format == :bam
            return XAM.BAM.Reader(open(xam))
        else
            error("XAM.jl parser does not support format: $format for file: $xam")
        end
    catch e
        error("Invalid XAM file for XAM.jl parser: ", xam, " - ", string(e))
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Open a XAM file using samtools view.
Supports SAM, BAM, and CRAM files with validation.
Always returns XAM.SAM.Reader regardless of input format.
"""
function open_xam_samtools(xam; header=false)
    Mycelia.add_bioconda_env("samtools")

    # First validate the file by testing samtools can read the header
    test_cmd = `$(Mycelia.CONDA_RUNNER) run -n samtools samtools view -H $(xam)`
    try
        # This will fail immediately for malformed files
        run(test_cmd, devnull, devnull, devnull)
    catch e
        error("Invalid or malformed SAM/BAM/CRAM file: $(xam) - samtools cannot parse the file")
    end

    # If validation passed, create the reader
    if header
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -h $(xam)`
    else
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view $(xam)`
    end
    try
        return XAM.SAM.Reader(open(cmd))
    catch e
        error("Invalid SAM/BAM/CRAM file for samtools parser: ", xam, " - ", string(e))
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate alignment statistics for a SAM/BAM/CRAM file using samtools flagstat.

# Arguments
- `xam::AbstractString`: Path to input SAM/BAM/CRAM alignment file
- `samtools_flagstat::AbstractString`: Output path for flagstat results (default: input_path.samtools-flagstat.txt)

# Returns
- `String`: Path to the generated flagstat output file

# Details
Runs samtools flagstat to calculate statistics on the alignment file, including:
- Total reads
- Secondary alignments
- Supplementary alignments  
- Duplicates
- Mapped/unmapped reads
- Proper pairs
- Read 1/2 counts

# Requirements
- Requires samtools to be available via Bioconda
- Input file must be in SAM, BAM or CRAM format
"""
function run_samtools_flagstat(xam, samtools_flagstat=xam * ".samtools-flagstat.txt")
    Mycelia.add_bioconda_env("samtools")
    if !isfile(samtools_flagstat)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools flagstat $(xam)`, samtools_flagstat))
    end
    return samtools_flagstat
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate mapping statistics by comparing sequence alignments (BAM/SAM) to a reference FASTA.

# Arguments
- `fasta::String`: Path to reference FASTA file
- `xam::String`: Path to alignment file (BAM or SAM format)

# Returns
DataFrame with columns:
- `contig`: Reference sequence name
- `contig_length`: Length of reference sequence
- `total_aligned_bases`: Total number of bases aligned to reference
- `mean_depth`: Average depth of coverage (total_aligned_bases/contig_length)
"""
function fasta_xam_mapping_stats(;fasta, xam)
    fastx_contig_lengths = fastx_to_contig_lengths(fasta)
    xam_stats = xam_to_contig_mapping_stats(xam)
    fastx_contig_lengths = fastx_to_contig_lengths(fasta)
    fastx_contig_lengths_table = DataFrames.DataFrame(contig = collect(keys(fastx_contig_lengths)), contig_length = collect(values(fastx_contig_lengths)))
    fastx_contig_mapping_stats_table = DataFrames.innerjoin(fastx_contig_lengths_table, xam_stats, on="contig" => "reference")
    mean_depth = fastx_contig_mapping_stats_table[!, "total_aligned_bases"] ./ fastx_contig_mapping_stats_table[!, "contig_length"]
    DataFrames.insertcols!(fastx_contig_mapping_stats_table, 4, :mean_depth => mean_depth)
    return fastx_contig_mapping_stats_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate detailed mapping statistics for each reference sequence/contig in a XAM (SAM/BAM/CRAM) file.

# Arguments
- `xam`: Path to XAM file or XAM object

# Returns
A DataFrame with per-contig statistics including:
- `n_aligned_reads`: Number of aligned reads
- `total_aligned_bases`: Sum of alignment lengths
- `total_alignment_score`: Sum of alignment scores
- Mapping quality statistics (mean, std, median)
- Alignment length statistics (mean, std, median)
- Alignment score statistics (mean, std, median)
- Percent mismatches statistics (mean, std, median)

Note: Only primary alignments (isprimary=true) and mapped reads (ismapped=true) are considered.
"""
function xam_to_contig_mapping_stats(xam)
    xam_results = Mycelia.parse_xam_to_summary_table(xam)
    xam_results = xam_results[xam_results[!, "isprimary"] .& xam_results[!, "ismapped"], :]
    # Calculate the percentage of mismatches
    xam_results.percent_mismatches = xam_results.mismatches ./ xam_results.alignlength * 100
    
    # Group by the 'reference' column and calculate the summary statistics
    contig_mapping_stats = DataFrames.combine(DataFrames.groupby(xam_results, :reference)) do subdf
        mappingquality_stats = StatsBase.summarystats(subdf.mappingquality)
        alignlength_stats = StatsBase.summarystats(subdf.alignlength)
        alignment_score_stats = StatsBase.summarystats(subdf.alignment_score)
        # mismatches_stats = StatsBase.summarystats(subdf.mismatches)
        percent_mismatches_stats = StatsBase.summarystats(subdf.percent_mismatches)

        (n_aligned_reads = length(subdf[!, "alignlength"]),
         total_aligned_bases = sum(subdf[!, "alignlength"]),
         total_alignment_score = sum(subdf[!, "alignment_score"]),
         mappingquality_mean = mappingquality_stats.mean,
         mappingquality_std = mappingquality_stats.sd,
         # mappingquality_min = mappingquality_stats.min,
         mappingquality_median = mappingquality_stats.median,
         # mappingquality_max = mappingquality_stats.max,

         alignlength_mean = alignlength_stats.mean,
         alignlength_std = alignlength_stats.sd,
         # alignlength_min = alignlength_stats.min,
         alignlength_median = alignlength_stats.median,
         # alignlength_max = alignlength_stats.max,

         alignment_score_mean = alignment_score_stats.mean,
         alignment_score_std = alignment_score_stats.sd,
         # alignment_score_min = alignment_score_stats.min,
         alignment_score_median = alignment_score_stats.median,
         # alignment_score_max = alignment_score_stats.max,

         # mismatches_mean = mismatches_stats.mean,
         # mismatches_std = mismatches_stats.sd,
         # mismatches_min = mismatches_stats.min,
         # mismatches_median = mismatches_stats.median,
         # mismatches_max = mismatches_stats.max,

         percent_mismatches_mean = percent_mismatches_stats.mean,
         percent_mismatches_std = percent_mismatches_stats.sd,
         # percent_mismatches_min = percent_mismatches_stats.min,
         percent_mismatches_median = percent_mismatches_stats.median,
         # percent_mismatches_max = percent_mismatches_stats.max)
        )
    end
    return contig_mapping_stats
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate per-base genomic coverage from a BAM file using bedtools.

# Arguments
- `bam::String`: Path to input BAM file

# Returns
- `DataFrames.DataFrame`: DataFrame with columns: chromosome, position, coverage_depth

# Details
Uses bedtools genomecov to compute per-base coverage. Creates a coverage file
with the format: <chromosome> <position> <coverage_depth> and then reads it into a DataFrame.

# Dependencies
Requires bedtools (automatically installed in conda environment)
"""
function determine_fasta_coverage_from_bam(bam)
    Mycelia.add_bioconda_env("bedtools")
    genome_coverage_file = bam * ".coverage.txt"
    if !isfile(genome_coverage_file)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bedtools bedtools genomecov -d -ibam $(bam)`, genome_coverage_file))
    end
    # Read the coverage file into a DataFrame
    if isfile(genome_coverage_file) && filesize(genome_coverage_file) > 0
        coverage_df = CSV.read(genome_coverage_file, DataFrames.DataFrame;
                              header=[:chromosome, :position, :coverage_depth], delim='\t')
        return coverage_df
    else
        # Return empty DataFrame with expected columns if no coverage data
        return DataFrames.DataFrame(chromosome=String[], position=Int[], coverage_depth=Int[])
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build the CLI arguments for running CoverM in contig mode.

# Internal
Returns a `Vector{String}` suitable for interpolation into a `coverm` command.
"""
function _build_coverm_contig_args(; bam_files::Vector{String},
                                   reference_fasta::Union{Nothing,String},
                                   methods::Vector{String},
                                   threads::Int,
                                   min_read_percent_identity::Union{Nothing,Float64},
                                   min_covered_fraction::Union{Nothing,Float64},
                                   out_path::String,
                                   additional_args::Vector{String})
    args = String[]
    push!(args, "contig")
    push!(args, "--bam-files")
    append!(args, bam_files)

    if reference_fasta !== nothing
        push!(args, "--reference", reference_fasta)
    end

    push!(args, "--methods")
    append!(args, methods)

    push!(args, "--threads", string(threads))

    if min_read_percent_identity !== nothing
        push!(args, "--min-read-percent-identity", string(min_read_percent_identity))
    end
    if min_covered_fraction !== nothing
        push!(args, "--min-covered-fraction", string(min_covered_fraction))
    end

    push!(args, "--output-file", out_path)
    append!(args, additional_args)
    return args
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run CoverM in *contig* mode to compute per-contig coverage statistics.

# Arguments
- `bam_files::Vector{String}`: Paths to BAM files containing read alignments.
- `reference_fasta::Union{Nothing,String}=nothing`: Optional contig FASTA for reference lengths.
- `outdir::Union{Nothing,String}=nothing`: Output directory (defaults to `coverm_contig` next to the first BAM).
- `methods::Vector{String}=["mean", "covered_fraction"]`: Coverage statistics to compute.
- `threads::Int=Sys.CPU_THREADS`: Number of threads for CoverM.
- `min_covered_fraction::Union{Nothing,Float64}=nothing`: Optional `--min-covered-fraction` filter.
- `min_read_percent_identity::Union{Nothing,Float64}=nothing`: Optional `--min-read-percent-identity` filter.
- `output_tsv::Union{Nothing,String}=nothing`: Optional explicit output path.
- `additional_args::Vector{String}=String[]`: Extra CLI arguments passed directly to CoverM.
- `quiet::Bool=true`: Suppress CoverM stdout/stderr.

# Returns
`DataFrames.DataFrame` parsed from the CoverM contig output.
"""
function run_coverm_contig(; bam_files::Vector{String},
                           reference_fasta::Union{Nothing,String}=nothing,
                           outdir::Union{Nothing,String}=nothing,
                           methods::Vector{String}=["mean", "covered_fraction"],
                           threads::Int=Sys.CPU_THREADS,
                           min_covered_fraction::Union{Nothing,Float64}=nothing,
                           min_read_percent_identity::Union{Nothing,Float64}=nothing,
                           output_tsv::Union{Nothing,String}=nothing,
                           additional_args::Vector{String}=String[],
                           quiet::Bool=true)
    @assert !isempty(bam_files) "bam_files must be non-empty"
    for bam in bam_files
        @assert isfile(bam) "BAM file does not exist: $(bam)"
        @assert filesize(bam) > 0 "BAM file is empty: $(bam)"
    end
    if reference_fasta !== nothing
        @assert isfile(reference_fasta) "Reference FASTA file does not exist: $(reference_fasta)"
        @assert filesize(reference_fasta) > 0 "Reference FASTA file is empty: $(reference_fasta)"
    end
    @assert threads > 0 "Thread count must be positive: $(threads)"
    @assert !isempty(methods) "At least one CoverM method must be specified"
    @assert all(!isempty, methods) "CoverM methods cannot contain empty strings"

    resolved_outdir = isnothing(outdir) ? joinpath(dirname(first(bam_files)), "coverm_contig") : outdir
    out_path = isnothing(output_tsv) ? joinpath(resolved_outdir, "coverm_contig.tsv") : output_tsv
    mkpath(dirname(out_path))

    should_run = !isfile(out_path) || filesize(out_path) == 0
    if should_run
        Mycelia.add_bioconda_env("coverm")
        args = _build_coverm_contig_args(; bam_files, reference_fasta, methods, threads,
                                         min_read_percent_identity, min_covered_fraction,
                                         out_path, additional_args)
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n coverm coverm $args`
        if quiet
            run(pipeline(cmd, stdout=devnull, stderr=devnull))
        else
            run(cmd)
        end
        @assert isfile(out_path) "CoverM contig output file was not created: $(out_path)"
        @assert filesize(out_path) > 0 "CoverM contig output file is empty: $(out_path)"
    end

    return CSV.read(out_path, DataFrames.DataFrame; delim='\t', normalizenames=true)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build the CLI arguments for running CoverM in genome mode.

# Internal
Returns a `Vector{String}` suitable for interpolation into a `coverm` command.
"""
function _build_coverm_genome_args(; bam_files::Vector{String},
                                   genome_fasta_files::Union{Nothing,Vector{String}},
                                   genome_directory::Union{Nothing,String},
                                   genome_extension::String,
                                   methods::Vector{String},
                                   threads::Int,
                                   out_path::String,
                                   additional_args::Vector{String})
    args = String[]
    push!(args, "genome")
    push!(args, "--bam-files")
    append!(args, bam_files)

    if genome_fasta_files !== nothing
        push!(args, "--genome-fasta-files")
        append!(args, genome_fasta_files)
    elseif genome_directory !== nothing
        push!(args, "--genome-fasta-directory", genome_directory)
        push!(args, "--genome-fasta-extension", genome_extension)
    end

    push!(args, "--methods")
    append!(args, methods)

    push!(args, "--threads", string(threads))
    push!(args, "--output-file", out_path)
    append!(args, additional_args)
    return args
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run CoverM in *genome* mode to compute per-genome/bin coverage and abundance.

# Arguments
- `bam_files::Vector{String}`: Paths to BAM files containing read alignments.
- `genome_fasta_files::Union{Nothing,Vector{String}}=nothing`: Explicit list of genome/bin FASTA files.
- `genome_directory::Union{Nothing,String}=nothing`: Directory containing genome FASTAs (mutually exclusive with `genome_fasta_files`).
- `genome_extension::String="fa"`: Extension to match within `genome_directory`.
- `outdir::Union{Nothing,String}=nothing`: Output directory (defaults to `coverm_genome` next to the first BAM).
- `methods::Vector{String}=["relative_abundance", "mean_coverage"]`: Coverage/abundance metrics to compute.
- `threads::Int=Sys.CPU_THREADS`: Number of threads for CoverM.
- `output_tsv::Union{Nothing,String}=nothing`: Optional explicit output path.
- `additional_args::Vector{String}=String[]`: Extra CLI arguments passed directly to CoverM.
- `quiet::Bool=true`: Suppress CoverM stdout/stderr.

# Returns
`DataFrames.DataFrame` parsed from the CoverM genome output.
"""
function run_coverm_genome(; bam_files::Vector{String},
                           genome_fasta_files::Union{Nothing,Vector{String}}=nothing,
                           genome_directory::Union{Nothing,String}=nothing,
                           genome_extension::String="fa",
                           outdir::Union{Nothing,String}=nothing,
                           methods::Vector{String}=["relative_abundance", "mean_coverage"],
                           threads::Int=Sys.CPU_THREADS,
                           output_tsv::Union{Nothing,String}=nothing,
                           additional_args::Vector{String}=String[],
                           quiet::Bool=true)
    @assert !isempty(bam_files) "bam_files must be non-empty"
    for bam in bam_files
        @assert isfile(bam) "BAM file does not exist: $(bam)"
        @assert filesize(bam) > 0 "BAM file is empty: $(bam)"
    end

    if genome_fasta_files === nothing && genome_directory === nothing
        error("Specify either genome_fasta_files or genome_directory")
    end
    if genome_fasta_files !== nothing && genome_directory !== nothing
        error("Specify genome_fasta_files OR genome_directory, not both")
    end
    if genome_fasta_files !== nothing
        @assert !isempty(genome_fasta_files) "genome_fasta_files must be non-empty if provided"
        for genome in genome_fasta_files
            @assert isfile(genome) "Genome FASTA file does not exist: $(genome)"
            @assert filesize(genome) > 0 "Genome FASTA file is empty: $(genome)"
        end
    else
        @assert isdir(genome_directory) "Genome directory does not exist: $(genome_directory)"
    end

    @assert threads > 0 "Thread count must be positive: $(threads)"
    @assert !isempty(methods) "At least one CoverM method must be specified"
    @assert all(!isempty, methods) "CoverM methods cannot contain empty strings"
    @assert !isempty(genome_extension) "genome_extension must be non-empty"

    resolved_outdir = isnothing(outdir) ? joinpath(dirname(first(bam_files)), "coverm_genome") : outdir
    out_path = isnothing(output_tsv) ? joinpath(resolved_outdir, "coverm_genome.tsv") : output_tsv
    mkpath(dirname(out_path))

    should_run = !isfile(out_path) || filesize(out_path) == 0
    if should_run
        Mycelia.add_bioconda_env("coverm")
        args = _build_coverm_genome_args(; bam_files, genome_fasta_files, genome_directory,
                                         genome_extension, methods, threads, out_path, additional_args)
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n coverm coverm $args`
        if quiet
            run(pipeline(cmd, stdout=devnull, stderr=devnull))
        else
            run(cmd)
        end
        @assert isfile(out_path) "CoverM genome output file was not created: $(out_path)"
        @assert filesize(out_path) > 0 "CoverM genome output file is empty: $(out_path)"
    end

    return CSV.read(out_path, DataFrames.DataFrame; delim='\t', normalizenames=true)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a BAM file to FASTQ format with gzip compression.

# Arguments
- `bam`: Path to input BAM file
- `fastq`: Optional output path. Defaults to input path with ".fq.gz" extension

# Returns
- Path to the generated FASTQ file

# Details
- Uses samtools through conda environment
- Automatically skips if output file exists
- Output is gzip compressed
- Requires samtools to be available via conda

"""
function bam_to_fastq(;bam, fastq=bam * ".fq.gz")
    Mycelia.add_bioconda_env("samtools")
    bam_to_fastq_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools fastq $(bam)`
    gzip_cmd = `gzip`
    p = pipeline(bam_to_fastq_cmd, gzip_cmd)
    if !isfile(fastq)
        @time run(pipeline(p, fastq))
    else
        @info "$(fastq) already exists"
    end
    return fastq
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function parse_xam(xam; filter_unmapped=false, primary_only=false, min_mapping_quality=0, min_align_length=1)
#     if occursin(r"\.bam$", xam)
#         MODULE = XAM.BAM
#         io = open(xam)
#     elseif occursin(r"\.sam$", xam)
#         MODULE = XAM.SAM
#         io = open(xam)
#     elseif occursin(r"\.sam.gz$", xam)
#         MODULE = XAM.SAM
#         io = CodecZlib.GzipDecompressorStream(open(xam))
#     else
#         error("unrecognized file extension in file: $xam")
#     end
#     # reader = open(MODULE.Reader, io)
#     reader = MODULE.Reader(io)
#     header = reader.header
#     record_iterator = Iterators.filter(record -> true, reader)
#     if filter_unmapped
#         record_iterator = Iterators.filter(record -> MODULE.ismapped(record), record_iterator)
#     end
#     if primary_only
#         record_iterator = Iterators.filter(record -> MODULE.isprimary(record), record_iterator)
#     end
#     record_iterator = Iterators.filter(record -> MODULE.mappingquality(record) >= min_mapping_quality, record_iterator)
#     record_iterator = Iterators.filter(record -> MODULE.alignlength(record) >= min_align_length, record_iterator)
#     records = sort(collect(record_iterator), by=x->[MODULE.refname(x), MODULE.position(x)])
#     # reset header to specify sorted
#     header.metainfo[1] = MODULE.MetaInfo("HD", ["VN" => 1.6, "SO" => "coordinate"])
#     close(io)
#     return (;records, header)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a SAM/BAM file into a summary DataFrame containing alignment metadata.

# Arguments
- `xam::AbstractString`: Path to input SAM (.sam), BAM (.bam), or gzipped SAM (.sam.gz) file

# Returns
DataFrame with columns:
- `template`: Read name
- `flag`: SAM flag
- `reference`: Reference sequence name
- `position`: Alignment position range (start:end)
- `mappingquality`: Mapping quality score
- `alignment_score`: Alignment score (AS tag)
- `isprimary`: Whether alignment is primary
- `alignlength`: Length of the alignment
- `ismapped`: Whether read is mapped
- `mismatches`: Number of mismatches (NM tag)

Note: Only mapped reads are included in the output DataFrame.
"""
function parse_xam_to_summary_table(xam)
    record_table = DataFrames.DataFrame(
        template = String[],
        flag = UInt16[],
        reference = String[],
        position = UnitRange{Int}[],
        mappingquality = UInt8[],
        alignment_score = Int[],
        isprimary = Bool[],
        # cigar = String[],
        # rnext = String[],
        # pnext = Int[],
        # tlen = Int[],
        # sequence = BioSequences.LongDNA{4}[],
        # quality = UInt8[],
        alignlength = Int[],
        ismapped = Bool[],
        # alignment = BioAlignments.Alignment[],
        mismatches = Int[]
    )
    if occursin(r"\.bam$", xam)
        MODULE = XAM.BAM
        io = open(xam)
    elseif occursin(r"\.sam$", xam)
        MODULE = XAM.SAM
        io = open(xam)
    elseif occursin(r"\.sam.gz$", xam)
        MODULE = XAM.SAM
        io = CodecZlib.GzipDecompressorStream(open(xam))
    else
        error("unrecognized file extension in file: $xam")
    end
    # reader = open(MODULE.Reader, io)
    reader = MODULE.Reader(io)
    header = reader.header
    for record in reader
        if XAM.SAM.ismapped(record)
            # @assert !ismissing()
            row = (
                template = XAM.SAM.tempname(record),
                flag = XAM.flag(record),
                reference = XAM.SAM.refname(record),
                position = XAM.SAM.position(record):XAM.SAM.rightposition(record),
                mappingquality = XAM.SAM.mappingquality(record),
                # cigar = XAM.SAM.cigar(record),
                # rnext = XAM.SAM.nextrefname(record),
                # pnext = XAM.SAM.nextposition(record),
                # tlen = XAM.SAM.templength(record),
                # sequence = XAM.SAM.sequence(record),
                # quality = XAM.SAM.quality(record),
                alignlength = XAM.SAM.alignlength(record),
                ismapped = XAM.SAM.ismapped(record),
                isprimary = XAM.SAM.isprimary(record),
                # alignment = XAM.SAM.alignment(record),
                alignment_score = record["AS"],
                mismatches = record["NM"]
                )
            push!(record_table, row, promote=true)
        end
    end
    # records = sort(collect(record_iterator), by=x->[MODULE.refname(x), MODULE.position(x)])
    # reset header to specify sorted
    # header.metainfo[1] = MODULE.MetaInfo("HD", ["VN" => 1.6, "SO" => "coordinate"])
    close(io)
    # return (;records, header)
    return record_table
end

"""
Check if a BAM file is coordinate-sorted.

First attempts a fast header check. If the header is missing or inconclusive,
performs a thorough validation by reading through the file to verify sort order.

Returns true if coordinate-sorted, false otherwise.
"""
function is_bam_coordinate_sorted(bam::String)
    # Try fast header check first
    try
        Mycelia.add_bioconda_env("samtools")
        header_lines = Base.readlines(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -H $(bam)`)
        
        # If we have a header with SO:coordinate, trust it
        if Base.any(line -> Base.occursin(r"@HD.*SO:coordinate", line), header_lines)
            @info "BAM header indicates coordinate sorting" bam
            return true
        end
        
        # If header explicitly says NOT coordinate sorted, trust that too
        if Base.any(line -> Base.occursin(r"@HD.*SO:(queryname|unsorted)", line), header_lines)
            @info "BAM header indicates NOT coordinate sorted" bam
            return false
        end
        
        # Header exists but no SO tag or empty header - need to validate
        @info "BAM header missing or incomplete, validating sort order..." bam
    catch e
        @warn "Failed to read BAM header, validating sort order..." exception=e
    end
    
    # Thorough check: read through file and verify sort order
    try
        reader = Mycelia.open_xam(bam)
        prev_refid = -1
        prev_pos = -1
        record_count = 0
        max_check = 100_000  # Only check first 100k records for performance
        
        for record in reader
            record_count += 1
            if record_count > max_check
                @info "Validated first $max_check records, assuming sorted" bam
                Base.close(reader)
                return true
            end
            
            # Skip unmapped reads
            if !XAM.BAM.ismapped(record)
                continue
            end
            
            refid = XAM.BAM.refid(record)
            pos = XAM.BAM.position(record)
            
            # Check if sort order is violated
            if refid < prev_refid || (refid == prev_refid && pos < prev_pos)
                @info "BAM file is NOT coordinate-sorted (violation at record $record_count)" bam
                Base.close(reader)
                return false
            end
            
            prev_refid = refid
            prev_pos = pos
        end
        
        Base.close(reader)
        @info "BAM file validated as coordinate-sorted" bam
        return true
    catch e
        @warn "Failed to validate BAM sort order, assuming unsorted" exception=e bam
        return false
    end
end

"""
Sort a BAM file by coordinate using samtools sort.

Returns the path to the sorted BAM file. By default, creates a new file with .sorted.bam suffix,
but can be customized with the output_path parameter.

Keyword Arguments:
- threads: Number of threads to use for sorting (default: all available CPUs)
- output_path: Path for the sorted BAM file (default: input with .sorted.bam suffix)
"""
function sort_bam(input_bam::String; threads::Int=get_default_threads(), output_path::Union{String,Nothing}=nothing)
    sorted_bam = Base.isnothing(output_path) ? Base.replace(input_bam, ".bam" => ".sorted.bam") : output_path
    if !isfile(sorted_bam)
        Mycelia.add_bioconda_env("samtools")
        @info "Sorting BAM file..." input_bam sorted_bam threads
        Base.run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(sorted_bam).sort.tmp -o $(sorted_bam) $(input_bam)`)
    else
        @info "Target sorted BAM file already exists at: $(sorted_bam)"
    end    
    return sorted_bam
end

"""
Ensure a BAM file is coordinate-sorted and indexed.
Returns a tuple of (sorted_bam_path, index_path).

If the BAM file is already sorted, it will be indexed in place.
If the BAM file is not sorted, a new sorted BAM file will be created with a .sorted.bam suffix.

This function attempts to index the BAM file directly. If indexing fails due to sort order issues,
it will automatically sort the file and retry indexing.

Keyword Arguments:
- threads: Number of threads to use for sorting (default: all available CPUs)
- skip_sort_check: Skip pre-check for sort order and attempt indexing directly (default: false)
"""
function index_bam(bam_path::String; threads::Int=get_default_threads(), skip_sort_check::Bool=false)

    bai_path = bam_path * ".bai"
    # Check if index already exists
    if Base.isfile(bai_path)
        @info "BAM index already exists at: $bai_path"
        return (bam_path, bai_path)
    end

    current_bam = bam_path
    # Optionally check sort order first (useful for large files where we want to avoid failed index attempts)
    if !skip_sort_check
        if !is_bam_coordinate_sorted(current_bam)
            current_bam = sort_bam(current_bam; threads=threads)
        end
    end
    bai_path = current_bam * ".bai"
    
    # Check if index already exists
    if Base.isfile(bai_path)
        @info "BAM index already exists at: $bai_path"
        return (current_bam, bai_path)
    end
    
    # Try to generate the index
    Mycelia.add_bioconda_env("samtools")
    @info "Generating BAM index for: $current_bam"
    try
        Base.run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools index $(current_bam)`)
        return (current_bam, bai_path)
    catch e
        # Check if error is due to sort order
        error_msg = Base.string(e)
        if Base.occursin(r"NO_COOR.*not in a single block", error_msg) || 
           Base.occursin(r"cannot be indexed", error_msg)
            @warn "Indexing failed due to sort order, sorting BAM file now..." current_bam
            current_bam = sort_bam(current_bam; threads=threads)
            
            # Try indexing the sorted file
            bai_path = current_bam * ".bai"
            @info "Generating BAM index for sorted file: $current_bam"
            Base.run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools index $(current_bam)`)
            return (current_bam, bai_path)
        else
            # Some other error occurred
            Base.rethrow(e)
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate per-base or per-region read depth from BAM files using mosdepth.

# Arguments
- `bam::String`: Path to input BAM or CRAM file
- `prefix::String`: Output prefix for generated files (default: bam path)
- `threads::Int=0`: Number of BAM decompression threads
- `use_median::Bool=false`: Output median coverage instead of mean for regions
- `fast_mode::Bool=false`: Skip internal CIGAR operations and mate overlap correction (faster)
- `no_per_base::Bool=false`: Skip per-base depth output (recommended for speed)
- `by::String=""`: Optional BED file path or window size for regional coverage
- `quantize::String=""`: Quantize coverage into bins (e.g., "0:1:10:50:")
- `thresholds::String=""`: Comma-separated coverage thresholds (e.g., "1,5,10,20")
- `mapq::Int=0`: Minimum mapping quality threshold
- `flag::Int=0`: Exclude reads with these SAM flags set
- `include_flag::Int=0`: Only include reads with these SAM flags set
- `fasta::String=""`: Reference FASTA for CRAM files
- `force::Bool=false`: Rerun even if output files exist

# Common Flag Filtering Strategies

The `flag` parameter excludes reads by summing SAM flag values. Common choices:

1. No filtering (default): `flag=0`
2. Exclude unmapped: `flag=4`
3. Exclude unmapped + failing QC: `flag=516` (4 + 512)
4. Exclude unmapped + failing QC + duplicates: `flag=1540` (4 + 512 + 1024)
5. Exclude unmapped + failing QC + duplicates + secondary alignments: `flag=1796` (4 + 512 + 1024 + 256)
6. Exclude unmapped + failing QC + duplicates + secondary + supplementary: `flag=3844` (4 + 512 + 1024 + 256 + 2048)

SAM Flag Reference:
- 4: Read unmapped
- 256: Secondary alignment (alternative mapping of same read)
- 512: Failed quality checks
- 1024: PCR or optical duplicate
- 2048: Supplementary alignment (chimeric, part of read maps elsewhere)

Note: Options 5 and 6 are most common for variant calling and coverage analysis.
For counting all aligned reads including multi-mappers, use option 4.

# Returns
Named tuple containing paths to output files:
- `prefix::String`: The output prefix used
- `global_dist::String`: Path to global distribution file
- `summary::String`: Path to summary statistics file
- `per_base::String`: Path to per-base coverage (if generated)
- `regions::String`: Path to regional coverage (if --by specified)
- `quantized::String`: Path to quantized coverage (if --quantize specified)
- `thresholds_file::String`: Path to thresholds file (if --thresholds specified)

# Details
mosdepth is a fast BAM/CRAM depth calculation tool that computes:
- Per-base depth across the genome
- Coverage distributions
- Regional depth statistics
- Quantized coverage bins

The function automatically:
- Installs mosdepth via conda if needed
- Indexes the BAM file if index doesn't exist
- Skips computation if output files already exist (unless force=true)

# Output Files
- `{prefix}.mosdepth.global.dist.txt`: Global coverage distribution
- `{prefix}.mosdepth.summary.txt`: Summary statistics per chromosome
- `{prefix}.per-base.bed.gz`: Per-base depth (unless no_per_base=true)
- `{prefix}.regions.bed.gz`: Regional depth (if by is specified)
- `{prefix}.quantized.bed.gz`: Quantized coverage (if quantize is specified)
- `{prefix}.thresholds.bed.gz`: Threshold coverage (if thresholds is specified)
"""
function run_mosdepth(bam::String;
                      prefix::String="",
                      threads::Int=0,
                      use_median::Bool=false,
                      fast_mode::Bool=false,
                      no_per_base::Bool=false,
                      by::String="",
                      quantize::String="",
                      thresholds::String="",
                      mapq::Int=0,
                      flag::Int=0,
                      include_flag::Int=0,
                      fasta::String="",
                      force::Bool=false)


    if isempty(prefix)
        prefix = bam
    end
    # Check if output already exists
    summary_file = prefix * ".mosdepth.summary.txt"
    if isfile(summary_file) && !force
        @info "mosdepth output already exists at: $(summary_file)"
        return (;
            prefix=prefix,
            global_dist=prefix * ".mosdepth.global.dist.txt",
            summary=summary_file,
            per_base=no_per_base ? "" : prefix * ".per-base.bed.gz",
            regions=isempty(by) ? "" : prefix * ".regions.bed.gz",
            quantized=isempty(quantize) ? "" : prefix * ".quantized.bed.gz",
            thresholds_file=isempty(thresholds) ? "" : prefix * ".thresholds.bed.gz"
        )
    end
    
    # Ensure BAM is sorted and indexed
    sorted_bam, bai_path = Mycelia.index_bam(bam, threads=threads)
    if isempty(prefix) || ((prefix == bam) && (sorted_bam != bam))
        prefix = sorted_bam
        # re-check if output already exists
        summary_file = prefix * ".mosdepth.summary.txt"
        if isfile(summary_file) && !force
            @info "mosdepth output already exists at: $(summary_file)"
            return (;
                prefix=prefix,
                global_dist=prefix * ".mosdepth.global.dist.txt",
                summary=summary_file,
                per_base=no_per_base ? "" : prefix * ".per-base.bed.gz",
                regions=isempty(by) ? "" : prefix * ".regions.bed.gz",
                quantized=isempty(quantize) ? "" : prefix * ".quantized.bed.gz",
                thresholds_file=isempty(thresholds) ? "" : prefix * ".thresholds.bed.gz"
            )
        end
    end
    
    # Build command arguments
    cmd_args = String["mosdepth"]
    
    # Add flags
    push!(cmd_args, "--threads", string(threads))
    push!(cmd_args, "--flag", string(flag))
    push!(cmd_args, "--include-flag", string(include_flag))
    push!(cmd_args, "--mapq", string(mapq))
    
    if use_median
        push!(cmd_args, "--use-median")
    end
    
    if fast_mode
        push!(cmd_args, "--fast-mode")
    end
    
    if no_per_base
        push!(cmd_args, "--no-per-base")
    end
    
    if !isempty(by)
        push!(cmd_args, "--by", by)
    end
    
    if !isempty(quantize)
        push!(cmd_args, "--quantize", quantize)
    end
    
    if !isempty(thresholds)
        push!(cmd_args, "--thresholds", thresholds)
    end
    
    if !isempty(fasta)
        push!(cmd_args, "--fasta", fasta)
    end
    
    # Add positional arguments
    push!(cmd_args, prefix, sorted_bam)
    
    # Run mosdepth
    @info "Running mosdepth on $(sorted_bam)"
    # Ensure mosdepth is installed
    Mycelia.add_bioconda_env("mosdepth")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mosdepth $(cmd_args)`)
    
    return (;
        prefix=prefix,
        global_dist=prefix * ".mosdepth.global.dist.txt",
        summary=summary_file,
        per_base=no_per_base ? "" : prefix * ".per-base.bed.gz",
        regions=isempty(by) ? "" : prefix * ".regions.bed.gz",
        quantized=isempty(quantize) ? "" : prefix * ".quantized.bed.gz",
        thresholds_file=isempty(thresholds) ? "" : prefix * ".thresholds.bed.gz"
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse mosdepth global distribution file for coverage QC metrics.

# Arguments
- `dist_file::String`: Path to mosdepth .global.dist.txt or .region.dist.txt file

# Returns
DataFrames.DataFrame with columns:
- `chromosome::String`: Chromosome name or "total" for genome-wide
- `coverage::Int`: Coverage depth threshold
- `proportion::Float64`: Proportion of bases with at least this coverage

# Details
The distribution file contains cumulative coverage data where each row shows
what proportion of bases have at least X coverage. The file format is:
- Column 1: chromosome name (or "total")
- Column 2: coverage level
- Column 3: proportion of bases covered at that level

For QC, focus on "total" rows to get genome-wide metrics like:
- What proportion of bases have ≥10X coverage
- What proportion of bases have ≥30X coverage
- Overall coverage uniformity
"""
function parse_mosdepth_distribution(dist_file::String)
    if !isfile(dist_file)
        error("Distribution file does not exist: $(dist_file)")
    end
    
    df = CSV.read(dist_file, DataFrames.DataFrame;
                  header=[:chromosome, :coverage, :proportion],
                  delim='\t',
                  comment="#")
    
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse mosdepth summary file for per-chromosome coverage statistics.

# Arguments
- `summary_file::String`: Path to mosdepth .summary.txt file

# Returns
DataFrames.DataFrame with columns:
- `chromosome::String`: Chromosome name or "total" for genome-wide
- `length::Int`: Length of chromosome in bases
- `bases::Int`: Number of bases with coverage data
- `mean::Float64`: Mean coverage depth
- `min::Float64`: Minimum coverage depth
- `max::Float64`: Maximum coverage depth

# Details
The summary file provides basic statistics for each chromosome including
mean, min, and max coverage. Use this for quick overview of coverage levels.
"""
function parse_mosdepth_summary(summary_file::String)
    if !isfile(summary_file)
        error("Summary file does not exist: $(summary_file)")
    end
    
    df = CSV.read(summary_file, DataFrames.DataFrame; delim='\t')
    
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse mosdepth thresholds BED file for per-region coverage QC.

# Arguments
- `thresholds_file::String`: Path to mosdepth .thresholds.bed.gz file

# Returns
DataFrames.DataFrame with columns:
- `chrom::String`: Chromosome name
- `start::Int`: Region start position (0-based)
- `end::Int`: Region end position
- `region::String`: Region name from input BED file (or "unknown")
- Plus one column per threshold value (e.g., `1X`, `10X`, `30X`)

# Details
Each threshold column contains the count of bases in that region with
at least that much coverage. Use this to identify regions with poor coverage.

For regions with good coverage, most bases should meet your threshold.
For example, if a 1000bp exon has only 500 bases at ≥10X, that indicates
poor coverage quality for that exon.
"""
function parse_mosdepth_thresholds(thresholds_file::String)
    if !isfile(thresholds_file)
        error("Thresholds file does not exist: $(thresholds_file)")
    end
    
    # Read with automatic header detection
    # mosdepth writes a header line starting with #
    df = CSV.read(thresholds_file, DataFrames.DataFrame; 
                  delim='\t',
                  comment="#")
    
    # If no header was found, mosdepth uses these column names
    if isempty(names(df))
        error("Unable to parse thresholds file - check file format")
    end
    
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract key QC metrics from mosdepth distribution data.

# Arguments
- `dist_df::DataFrames.DataFrame`: DataFrame from parse_mosdepth_distribution
- `thresholds::Vector{Int}=[1, 3, 5, 10, 30, 50, 100, 300, 500, 1000]`: Coverage thresholds to report

# Returns
DataFrames.DataFrame with one row per chromosome containing proportion of bases
meeting each coverage threshold.

# Details
Extracts the proportion of bases with at least X coverage for each specified
threshold. Focuses on interpretable QC metrics rather than full distribution.

For whole genome sequencing QC, typical expectations are:
- ≥1X: >95% (most of genome accessible)
- ≥10X: >90% (sufficient for variant calling)
- ≥30X: >80% (high confidence variant calling)
"""
function summarize_mosdepth_qc(dist_df::DataFrames.DataFrame; 
                               thresholds::Vector{Int}=[1, 3, 5, 10, 30, 50, 100, 300, 500, 1000])
    
    # Get unique chromosomes
    chroms = unique(dist_df.chromosome)
    
    # Initialize result dataframe with columns that can hold missing values
    result = DataFrames.DataFrame(chromosome=String[])
    for threshold in thresholds
        result[!, Symbol("coverage_$(threshold)X")] = Union{Float64, Missing}[]
    end
    
    # For each chromosome, find the proportion at each threshold
    for chrom in chroms
        chrom_data = DataFrames.subset(dist_df, :chromosome => x -> x .== chrom)
        
        row_data = Dict{Symbol, Any}(:chromosome => chrom)
        
        for threshold in thresholds
            # Find the row with this exact coverage value
            matching_rows = DataFrames.subset(chrom_data, :coverage => x -> x .== threshold)
            
            if DataFrames.nrow(matching_rows) > 0
                row_data[Symbol("coverage_$(threshold)X")] = matching_rows[1, :proportion]
            else
                # If exact threshold not found, use missing value
                row_data[Symbol("coverage_$(threshold)X")] = missing
            end
        end
        
        push!(result, row_data)
    end
    
    return result
end
