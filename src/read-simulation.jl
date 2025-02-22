"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina short reads from a FASTA file using the ART Illumina simulator.

This function wraps ART (installed via Bioconda) to simulate reads from an input
reference FASTA. It supports paired-end (or optionally single-end/mate-pair) simulation,
with options to choose either fold coverage (`--fcov`) or an absolute read count (`--rcount`),
to enable amplicon mode, and to optionally generate a zero-error SAM file.

# Arguments
- `in_fasta::String`: Path to the input FASTA file.
- `coverage::Union{Nothing,Number}`: Desired fold coverage (used with `--fcov`); if `nothing`
  and `read_count` is provided then fold coverage is ignored. (Default: 20)
- `read_count::Union{Nothing,Number}`: Total number of reads (or read pairs) to generate
  (used with `--rcount` instead of fold coverage). (Default: `nothing`)
- `outbase::String`: Output file prefix (default: "\$(in_fasta).art.\$(coverage)x.").
- `read_length::Int`: Length of reads to simulate (default: 150).
- `mflen::Int`: Mean fragment length for paired-end simulations (default: 500).
- `sdev::Int`: Standard deviation of fragment lengths (default: 10).
- `seqSys::String`: Illumina sequencing system ID (e.g. "HS25" for HiSeq 2500) (default: "HS25").
- `paired::Bool`: Whether to simulate paired-end reads (default: true).
- `amplicon::Bool`: Enable amplicon sequencing simulation mode (default: false).
- `errfree::Bool`: Generate an extra SAM file with zero sequencing errors (default: false).
- `noALN::Bool`: Do not output the ALN alignment file (default: true).
- `rndSeed::Union{Nothing,Int}`: Optional seed for reproducibility (default: nothing).

# Outputs
Generates gzipped FASTQ files in the working directory:
- For paired-end: `\$(outbase)1.fq.gz` (forward) and `\$(outbase)2.fq.gz` (reverse).
- For single-end: `\$(outbase)1.fq.gz`.

Additional SAM files may be produced if `--errfree` is enabled and/or if
the ART `--samout` option is specified.

# Details
This function calls ART with the provided options. Note that if `read_count` is supplied,
the function uses the `--rcount` option; otherwise, it uses `--fcov` with the given coverage.
Amplicon mode (via `--amplicon`) restricts the simulation to the amplicon regions, which is
important for targeted sequencing studies.

# Dependencies
Requires ART simulator (installed via Bioconda) and the Mycelia environment helper.

See also: `simulate_nanopore_reads`, `simulate_nearly_perfect_long_reads`, `simulate_pacbio_reads`
"""
function simulate_illumina_paired_reads(;in_fasta::String,
    coverage::Union{Nothing, Number}=20,
    read_count::Union{Nothing, Number}=nothing,
    outbase::String = "$(in_fasta).art.$(coverage)x.",
    read_length::Int = 150,
    mflen::Int = 500,
    sdev::Int = 10,
    seqSys::String = "HS25",
    paired::Bool = true,
    amplicon::Bool = false,
    errfree::Bool = false,
    noALN::Bool = true,
    rndSeed::Union{Nothing,Int} = nothing
    )

    # Ensure ART is available via Bioconda
    Mycelia.add_bioconda_env("art")

    # Precompute option values
    amplicon_flag = amplicon ? "--amplicon" : ""
    errfree_flag = errfree ? "--errfree" : ""
    noALN_flag = noALN ? "--noALN" : ""
    rndSeed_val = (rndSeed !== nothing && rndSeed > -1) ? "$(rndSeed)" : ""

    # Determine read count or coverage
    if read_count !== nothing
        rcount_val = "$(read_count)"
        fcov_val = ""
    elseif coverage !== nothing
        rcount_val = ""
        fcov_val = "$(coverage)"
    else
        error("Either 'coverage' or 'read_count' must be provided.")
    end

    # Input and output values
    in_fasta_val = "$(in_fasta)"
    outbase_val = "$(outbase)"

    # Full command with all options explicitly defined
    full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
        --seqSys $(seqSys) \
        --len $(read_length) \
        $(paired_flag) \
        --mflen $(mflen_val) \
        --sdev $(sdev_val) \
        $(amplicon_flag) \
        $(errfree_flag) \
        $(noALN_flag) \
        --rndSeed $(rndSeed_val) \
        --rcount $(rcount_val) \
        --fcov $(fcov_val) \
        --in $(in_fasta_val) \
        --out $(outbase_val)`

    @info "Running ART with command: $(full_cmd)"
    @time run(full_cmd)

    # Process output FASTQ files: gzip them if not already compressed.
    # For paired-end, output files are expected as outbase1.fq and outbase2.fq.
    # For single-end, only outbase1.fq is produced.
    forward_fq = "$(outbase)1.fq"
    forward_gz = forward_fq * ".gz"
    @assert isfile(forward_fq) "Forward FASTQ file not found: $(forward_fq)"
    run(`gzip $(forward_fq)`)
    @assert isfile(forward_gz) "Gzipped forward FASTQ not found: $(forward_gz)"


    if paired
        reverse_fq = "$(outbase)2.fq"
        reverse_gz = reverse_fq * ".gz"
        @assert isfile(reverse_fq) "Reverse FASTQ file not found: $(reverse_fq)"
        run(`gzip $(reverse_fq)`)
        @assert isfile(reverse_gz) "Gzipped reverse FASTQ not found: $(reverse_gz)"
    else
        reverse_fq = nothing
        reverse_gz = nothing
    end

    # Optionally, you could remove ALN or SAM files here if not needed.
    # return nothing
    return (forward_reads = forward_gz, reverse_reads = reverse_gz)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate PacBio HiFi reads using the Badread error model.

# Arguments
- `fasta::String`: Path to input FASTA file containing reference sequence
- `quantity::String`: Coverage depth (e.g. "50x") or total bases (e.g. "1000000") - NOT TOTAL READS
- `outfile::String`: Output filepath for simulated reads. Defaults to input filename with ".badread.pacbio2021.\${quantity}.fq.gz" suffix

# Returns
- `String`: Path to the generated output file

# Notes
- Requires Badread tool from Bioconda
- Uses PacBio 2021 error and quality score models
- Average read length ~15kb
- Output is gzipped FASTQ format

See also: `simulate_nanopore_reads`, `simulate_nearly_perfect_long_reads`, `simulate_short_reads`
"""
function simulate_pacbio_reads(;fasta, quantity, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.pacbio2021.$(quantity).fq.gz"))
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --error_model pacbio2021 --qscore_model pacbio2021 --identity 30,3 --reference $(fasta) --quantity $(quantity)`, `gzip`)
        run(pipeline(p, outfile))
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Oxford Nanopore sequencing reads using the Badread tool with 2023 error models.

# Arguments
- `fasta::String`: Path to input reference FASTA file
- `quantity::String`: Either fold coverage (e.g. "50x") or total bases to sequence (e.g. "1000000")
- `outfile::String`: Output path for gzipped FASTQ file. Defaults to input filename with modified extension

# Returns
- `String`: Path to the generated output FASTQ file

See also: `simulate_pacbio_reads`, `simulate_nearly_perfect_long_reads`, `simulate_short_reads`
"""
function simulate_nanopore_reads(;fasta, quantity, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.nanopore2023.$(quantity).fq.gz"))
# badread simulate --reference ref.fasta --quantity 50x | gzip > reads.fastq.gz
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --error_model nanopore2023 --qscore_model nanopore2023 --reference $(fasta) --quantity $(quantity)`, `gzip`)
        run(pipeline(p, outfile))
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate high-quality long reads with minimal errors using Badread.

# Arguments
- `reference::String`: Path to reference FASTA file
- `quantity::String`: Coverage depth (e.g. "50x") or total bases (e.g. "1000000")
- `length_mean::Int=40000`: Mean read length
- `length_sd::Int=20000`: Standard deviation of read length

# Returns
Vector of simulated reads in FASTQ format

# Details
Generates nearly perfect long reads by setting error rates and artifacts to minimum values.
Uses ideal quality scores and disables common sequencing artifacts like chimeras and adapters.

See also: `simulate_pacbio_reads`, `simulate_nanopore_reads`, `simulate_short_reads`
"""
function simulate_nearly_perfect_long_reads()
    @error "finish implementing me"
    # badread simulate --reference ref.fasta --quantity 50x --error_model random \
    # --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
    # --identity 30,3 --length 40000,20000 --start_adapter_seq "" --end_adapter_seq "" \
    # | gzip > reads.fastq.gz
end