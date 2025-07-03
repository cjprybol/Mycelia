"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert SAM/BAM records from a XAM.SAM.Reader into a DataFrame.

Parameters:
- `reader`: A XAM.SAM.Reader object for iterating over records

Returns:
- A DataFrame containing all record data in a structured format
"""
function xam_to_dataframe(reader::XAM.SAM.Reader)::DataFrames.DataFrame
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
    seqs = String[]
    seqs = Union{Nothing, BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}[]
    quals = Union{Missing, Vector{UInt8}}[]
    alignlengths = Union{Int, Missing}[]
    alignment_score = Union{Int, Missing}[]
    mismatches = Union{Int, Missing}[]
    # aux_tags = Vector{Dict{String,Any}}()
    
    # Process each record
    for record in reader
        try
            push!(templates, XAM.SAM.tempname(record))
        catch err
            println("Error extracting tempname from record: ", record)
            rethrow(err)
        end
        
        try
            push!(ismapped, XAM.SAM.ismapped(record))
        catch err
            println("Error extracting ismapped from record: ", record)
            rethrow(err)
        end
        
        try
            push!(isprimary, XAM.SAM.isprimary(record))
        catch err
            println("Error extracting isprimary from record: ", record)
            rethrow(err)
        end

        try
            push!(flags, XAM.SAM.flag(record))
        catch err
            println("Error extracting flag from record: ", record)
            rethrow(err)
        end

        try
            push!(references, XAM.SAM.ismapped(record) ? XAM.SAM.refname(record) : missing)
        catch err
            println("Error extracting refname from record: ", record)
            rethrow(err)
        end

        try
            push!(positions, XAM.SAM.position(record):XAM.SAM.rightposition(record))
        catch err
            println("Error extracting positions from record: ", record)
            rethrow(err)
        end

        try
            push!(mappingqualities, XAM.SAM.mappingquality(record))
        catch err
            println("Error extracting mappingquality from record: ", record)
            rethrow(err)
        end

        try
            push!(cigars, XAM.SAM.cigar(record))
        catch err
            println("Error extracting cigar from record: ", record)
            rethrow(err)
        end

        try
            push!(rnexts, XAM.SAM.nextposition(record) == 0 ? missing : XAM.SAM.nextrefname(record))
        catch err
            println("Error extracting rnext from record: ", record)
            rethrow(err)
        end

        try
            push!(pnexts, XAM.SAM.nextposition(record))
        catch err
            println("Error extracting pnext from record: ", record)
            rethrow(err)
        end

        try
            push!(tlens, XAM.SAM.templength(record))
        catch err
            println("Error extracting templength from record: ", record)
            rethrow(err)
        end

        try
            push!(seqs, XAM.SAM.sequence(record))
        catch err
            println("Error extracting sequence from record: ", record)
            rethrow(err)
        end

        try
            push!(quals, XAM.SAM.isprimary(record) ? XAM.SAM.quality(record) : missing)
        catch err
            println("Error extracting quality from record: ", record)
            rethrow(err)
        end

        try
            push!(alignlengths, XAM.SAM.alignlength(record))
        catch err
            println("Error extracting alignlength from record: ", record)
            rethrow(err)
        end

        try
            push!(alignment_score, XAM.SAM.ismapped(record) ? record["AS"] : missing)
        catch err
            println("Error extracting alignment_score from record: ", record)
            rethrow(err)
        end

        try
            push!(mismatches,  XAM.SAM.ismapped(record) ? record["NM"] : missing)
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
    return xam_to_dataframe(open_xam(xam_path))    
end

function open_xam(xam; header=false)
    Mycelia.add_bioconda_env("samtools")
    if header
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -h $(xam)`
    else
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view $(xam)`
    end
    return XAM.SAM.Reader(open(cmd))
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
- `String`: Path to the generated coverage file (`.coverage.txt`)

# Details
Uses bedtools genomecov to compute per-base coverage. Creates a coverage file 
with the format: <chromosome> <position> <coverage_depth>. 
If the coverage file already exists, returns the existing file path.

# Dependencies
Requires bedtools (automatically installed in conda environment)
"""
function determine_fasta_coverage_from_bam(bam)
    Mycelia.add_bioconda_env("bedtools")
    genome_coverage_file = bam * ".coverage.txt"
    if !isfile(genome_coverage_file)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bedtools bedtools genomecov -d -ibam $(bam)`, genome_coverage_file))
    end
    return genome_coverage_file
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