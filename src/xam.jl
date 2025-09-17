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