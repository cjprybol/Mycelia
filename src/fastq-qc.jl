"""
$(DocStringExtensions.TYPEDSIGNATURES)

Trim paired-end FASTQ reads using Trim Galore, a wrapper around Cutadapt and FastQC.

# Arguments
- `outdir::String`: Output directory containing input FASTQ files
- `identifier::String`: Prefix for input/output filenames

# Input files
Expects paired FASTQ files in `outdir` named:
- `{identifier}_1.fastq.gz` (forward reads)
- `{identifier}_2.fastq.gz` (reverse reads)

# Output files
Creates trimmed reads in `outdir/trim_galore/`:
- `{identifier}_1_val_1.fq.gz` (trimmed forward reads)
- `{identifier}_2_val_2.fq.gz` (trimmed reverse reads)

# Dependencies
Requires trim_galore conda environment:
"""
function trim_galore(;outdir="", identifier="")
    
    trim_galore_dir = joinpath(outdir, "trim_galore")
    
    forward_reads = joinpath(outdir, "$(identifier)_1.fastq.gz")
    reverse_reads = joinpath(outdir, "$(identifier)_2.fastq.gz")
    
    trimmed_forward_reads = joinpath(trim_galore_dir, "$(identifier)_1_val_1.fq.gz")
    trimmed_reverse_reads = joinpath(trim_galore_dir, "$(identifier)_2_val_2.fq.gz")
    
    # mamba create -n trim_galore -c bioconda trim_galore
    if !isfile(trimmed_forward_reads) && !isfile(trimmed_reverse_reads)
        cmd = `conda run -n trim_galore trim_galore --suppress_warn --cores $(min(Sys.CPU_THREADS, 4)) --output_dir $(trim_galore_dir) --paired $(forward_reads) $(reverse_reads)`
        run(cmd)
    else
        @info "$(trimmed_forward_reads) & $(trimmed_reverse_reads) already present"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Trim paired-end FASTQ reads using Trim Galore, a wrapper around Cutadapt and FastQC.

# Arguments
- `forward_reads::String`: Path to forward reads FASTQ file
- `reverse_reads::String`: Path to reverse reads FASTQ file
- `outdir::String`: Output directory for trimmed files

# Returns
- `Tuple{String, String}`: Paths to trimmed forward and reverse read files

# Dependencies
Requires trim_galore conda environment:
- `mamba create -n trim_galore -c bioconda trim_galore`
"""
function trim_galore_paired(;forward_reads::String, reverse_reads::String, outdir::String=pwd())
    # Create output directory if it doesn't exist
    trim_galore_dir = mkpath(joinpath(outdir, "trim_galore"))
    
    # Get base filename without path and extension for output naming
    forward_base = basename(forward_reads)
    reverse_base = basename(reverse_reads)
    
    # Construct output filenames according to trim_galore naming convention
    trimmed_forward = joinpath(trim_galore_dir, replace(forward_base, Mycelia.FASTQ_REGEX => "_val_1.fq.gz"))
    trimmed_reverse = joinpath(trim_galore_dir, replace(reverse_base, Mycelia.FASTQ_REGEX => "_val_2.fq.gz"))
    
    if !isfile(trimmed_forward) && !isfile(trimmed_reverse)
        cmd = `conda run -n trim_galore trim_galore --suppress_warn --cores $(min(Sys.CPU_THREADS, 4)) --output_dir $(trim_galore_dir) --paired $(forward_reads) $(reverse_reads)`
        run(cmd)
    else
        @info "$(trimmed_forward) & $(trimmed_reverse) already present"
    end
    
    return (;trimmed_forward, trimmed_reverse, outdir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform quality control (QC) filtering and trimming on short-read FASTQ files using fastp.

# Arguments
- `in_fastq::String`: Path to the input FASTQ file.
- `out_fastq::String`: Path to the output FASTQ file.
- `adapter_seq::String`: Adapter sequence to trim.
- `quality_threshold::Int`: Minimum phred score for trimming (default 20).
- `min_length::Int`: Minimum read length to retain (default 50).

# Returns
- `String`: Path to the filtered and trimmed FASTQ file.

# Details
This function uses fastp to remove adapter contamination, trim low‐quality bases from the 3′ end,
and discard reads shorter than `min_length`. It’s a simple wrapper that executes the external fastp command.
"""
function qc_filter_short_reads_fastp(in_fastq::String, out_fastq::String; adapter_seq::String, quality_threshold::Int=20, min_length::Int=50)
    cmd = `fastp --in1 $(in_fastq) --out1 $(out_fastq) --adapter_sequence $(adapter_seq) --qualified_quality_phred $(quality_threshold) --length_required $(min_length)`
    run(cmd)
    return out_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform QC filtering on long-read FASTQ files using fastplong.

# Arguments
- `in_fastq::String`: Path to the input FASTQ file.
- `out_fastq::String`: Path to the output FASTQ file.
- `quality_threshold::Int`: Minimum average quality to retain a read (default 10).
- `min_length::Int`: Minimum read length (default 1000).
- `max_length::Int=0`: Maximum read length (default 0, no maximum).

# Returns
- `String`: Path to the filtered FASTQ file.

# Details
This function uses fastplong to filter long reads based on quality and length criteria.
It is optimized for Oxford Nanopore, PacBio, or similar long-read datasets.
"""
function qc_filter_long_reads_fastplong(;
                            in_fastq::String,
                            report_title::String=in_fastq * " fastplong report",
                            out_fastq::String=Mycelia.replace(in_fastq, Mycelia.FASTQ_REGEX => ".fastplong.fq.gz"),
                            html_report::String=out_fastq * ".html",
                            json_report::String=out_fastq * ".json",
                            min_length::Int=1000,
                            max_length::Int=0)
    # Build command with required parameters
    Mycelia.add_bioconda_env("fastplong")
    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastplong fastplong
            --in $(in_fastq)
            --out $(out_fastq)
            --report_title $(report_title)
            --html $(html_report)
            --json $(json_report)
            --length_required $(min_length)`
    # Add max length if specified
    if max_length > 0
        push!(cmd, "--length_limit")
        push!(cmd, string(max_length))
    end

    run(`$cmd`)
    return out_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Filter and process long reads from a FASTQ file using Filtlong.

This function filters long sequencing reads based on quality and length criteria, 
then compresses the output using pigz.

# Arguments
- `in_fastq::String`: Path to the input FASTQ file.
- `out_fastq::String`: Path to the output filtered and compressed FASTQ file. 
   Defaults to the input filename with ".filtlong.fq.gz" appended.
- `min_mean_q::Int`: Minimum mean quality score for reads to be kept. Default is 20.
- `keep_percent::Int`: Percentage of reads to keep after filtering. Default is 95.

# Returns
- `out_fastq`

# Details
This function uses Filtlong to filter long reads and pigz for compression. It requires
the Bioconda environment for Filtlong to be set up, which is handled internally.
"""
function qc_filter_long_reads_filtlong(;
        in_fastq,
        out_fastq = replace(in_fastq, r"\.(fq\.gz|fastq\.gz|fastq|fq)$" => ".filtlong.fq.gz"),
        min_mean_q = 20,
        keep_percent = 95
    )
    Mycelia.add_bioconda_env("filtlong")
    p1 = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n filtlong filtlong --min_mean_q $(min_mean_q) --keep_percent $(keep_percent) $(in_fastq)`,
        `pigz`
    )
    p2 = pipeline(p1, out_fastq)
    run(p2)
    return out_fastq
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Perform QC filtering on long-read FASTQ files using chopper.

# # Arguments
# - `in_fastq::String`: Path to the input FASTQ file.
# - `out_fastq::String`: Path to the output FASTQ file.
# - `quality_threshold::Int`: Minimum average quality to retain a read (default 10).
# - `min_length::Int`: Minimum read length (default 1000).

# # Returns
# - `String`: Path to the filtered FASTQ file.

# # Details
# This function uses chopper to discard long reads that do not meet the minimum quality or length thresholds.
# It is intended for Oxford Nanopore or similar long-read datasets.

# # Dependencies
# Requires chopper to be installed via conda
# """
# function qc_filter_long_reads_chopper(in_fastq::String, out_fastq::String; quality_threshold::Int=20, min_length::Int=1000)
#     # chopper reads from STDIN and writes to STDOUT, so we pipe the file
#     Mycelia.add_bioconda_env("chopper")
#     cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n chopper chopper --input (in_fastq) --quality $(quality_threshold) --minlength $(min_length)`
#     p = pipeline(`cat $(in_fastq)`, cmd)
#     , out_fastq)
#     open(out_fastq, "w") do outf
#         run(pipeline(`cat $(in_fastq)`, cmd; stdout=outf))
#     end
#     return out_fastq
# end