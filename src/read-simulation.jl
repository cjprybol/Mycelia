"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina paired-end short reads from a FASTA file using ART simulator.

# Arguments
- `in_fasta::String`: Input FASTA file path
- `coverage::Number`: Desired read coverage/depth
- `outbase::String`: Prefix for output files (default: "\${in_fasta}.art.\${coverage}x.")

# Outputs
Generates two gzipped FASTQ files:
- `\${outbase}1.fq.gz`: Forward reads
- `\${outbase}2.fq.gz`: Reverse reads

# Details
Uses ART Illumina with the following parameters:
- Read length: 150bp
- Fragment length: 500bp (SD: 10bp)
- Sequencing system: HiSeq 2500 (HS25)

# Dependencies
Requires ART simulator (automatically installed via Bioconda)

See also: `simulate_nanopore_reads`, `simulate_nearly_perfect_long_reads`, `simulate_pacbio_reads`
"""
function simulate_short_reads(;in_fasta, coverage, outbase = "$(in_fasta).art.$(coverage)x.")
    # -c --rcount
    # total number of reads/read pairs to be generated [per amplicon if for amplicon simulation](not be used together with -f/--fcov)
    # -d --id
    # the prefix identification tag for read ID
    # -ef --errfree
    # indicate to generate the zero sequencing errors SAM file as well the regular one
    # NOTE: the reads in the zero-error SAM file have the same alignment positions as those in the regular SAM file, but have no sequencing errors
    # -f --fcov
    # the fold of read coverage to be simulated or number of reads/read pairs generated for each amplicon
    # --samout
    Mycelia.add_bioconda_env("art")
    p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina --noALN --paired --seqSys HS25 --len 150 --mflen 500 --sdev 10 --in $(in_fasta) --out $(outbase) --fcov $(coverage)`)
    @time run(p)
    out_forward = "$(outbase)1.fq"
    target_forward = out_forward * ".gz"
    @assert isfile(out_forward)
    run(`gzip $(out_forward)`)
    @assert isfile(target_forward)

    out_reverse = "$(outbase)2.fq"
    target_reverse = out_reverse * ".gz"
    @assert isfile(out_reverse)
    run(`gzip $(out_reverse)`)
    @assert isfile(target_reverse)

    # isfile("$(outbase)1.aln") && rm("$(outbase)1.aln")
    # isfile("$(outbase)2.aln") && rm("$(outbase)2.aln")
    # isfile("$(outbase).sam") && rm("$(outbase).sam")
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