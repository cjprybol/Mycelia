const PRINTABLE_ASCII_ALPHABET = filter(isprint, [Char(c) for c in 0x21:0x7E])
const PRINTABLE_GREEK_ALPHABET = vcat(
    filter(isprint, [Char(c) for c in 0x391:0x3A9]),  # Capital Greek letters
    filter(isprint, [Char(c) for c in 0x3B1:0x3C9])   # Lowercase Greek letters
)
const PRINTABLE_LATIN1_ALPHABET = filter(isprint, [Char(c) for c in 0x00:0xFF])
const PRINTABLE_UNICODE_ALPHABET = [Char(c) for c in 0x0000:0x10FFFF if ((c <= 0xD7FF) || (0xE000 <= c <= 0x10FFFF)) && Base.isprint(Char(c))]
const PRINTABLE_BMP_ALPHABET  = [Char(c) for c in 0x0020:0xFFFD if (c < 0xD800 || c > 0xDFFF) && Base.isprint(Char(c))]

# Mutate a string given an alphabet, error rate, and allowed mutation types.
# If `alphabet` is not provided, use the unique set of observed characters in the input string.
function mutate_string(s::String; alphabet::Union{Nothing,AbstractVector{Char}}=nothing, error_rate::Float64=0.01)
    chars = collect(s)
    if alphabet === nothing
        alphabet = unique(chars)
    end
    i = 1
    while i <= length(chars)
        if Random.rand() < error_rate
            mutation = Random.rand(['sub', 'ins', 'del'])
            if mutation == 'sub'
                # Substitute with random char from alphabet (not the same as current)
                choices = setdiff(alphabet, [chars[i]])
                if !isempty(choices)
                    chars[i] = Random.rand(choices)
                end
            elseif mutation == 'ins'
                # Insert random char from alphabet
                insert_char = Random.rand(alphabet)
                insert!(chars, i, insert_char)
                i += 1 # skip inserted char
            elseif mutation == 'del' && length(chars) > 1
                # Delete the char
                deleteat!(chars, i)
                i -= 1 # stay at this index
            end
        end
        i += 1
    end
    return join(chars)
end

"""
    rand_ascii_greek_string(len::Int) -> String

Generate a random string of printable ASCII and Greek characters of length `len`.

The string contains random printable ASCII characters and both uppercase and lowercase Greek letters.
"""
function rand_ascii_greek_string(len::Int)
    alphabet = vcat(PRINTABLE_ASCII_ALPHABET, PRINTABLE_GREEK_ALPHABET)
    return join([Random.rand(alphabet) for _ in 1:len])
end

"""
    rand_latin1_string(len::Int) -> String

Generate a random string of printable Latin-1 characters of length `len`.

The string contains random printable characters from the Latin-1 character set.
"""
function rand_latin1_string(len::Int)
    return join([Random.rand(PRINTABLE_LATIN1_ALPHABET) for _ in 1:len])
end

"""
    rand_printable_unicode_string(len::Int) -> String

Generate a random string of printable Unicode characters of length `len`.

The string contains random printable Unicode characters, excluding surrogate code points.
"""
function rand_printable_unicode_string(len::Int)
    return join([Random.rand(PRINTABLE_UNICODE_ALPHABET) for _ in 1:len])
end

"""
    rand_bmp_printable_string(len::Int) -> String

Generate a random string of printable Basic Multilingual Plane (BMP) characters of length `len`.

The string contains random printable BMP characters, excluding surrogate code points.
"""
function rand_bmp_printable_string(len::Int)
    bmp_chars = [Char(c) for c in 0x0020:0xFFFD if (c < 0xD800 || c > 0xDFFF) && Base.isprint(Char(c))]
    return join([Random.rand(PRINTABLE_BMP_ALPHABET) for _ in 1:len])
end

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
    coverage::Union{Nothing, Number}=nothing,
    read_count::Union{Nothing, Int}=nothing,
    outbase::String = "",
    read_length::Int = 150,
    mflen::Int = 500,
    sdev::Int = 10,
    seqSys::String = "HS25",
    amplicon::Bool = false,
    errfree::Bool = true,
    rndSeed::Int = current_unix_datetime()
    )

    # Ensure ART is available via Bioconda
    Mycelia.add_bioconda_env("art")
    
    # Determine read count or coverage
    if read_count !== nothing
        if isempty(outbase)
            outbase = "$(in_fasta).rcount_$(read_count).art"
        end
        full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
            --paired \
            --errfree \
            --noALN \
            --seqSys $(seqSys) \
            --len $(read_length) \
            --mflen $(mflen) \
            --sdev $(sdev) \
            --rndSeed $(rndSeed) \
            --rcount $(read_count) \
            --in $(in_fasta) \
            --out $(outbase)`
    elseif coverage !== nothing
        if isempty(outbase)
            outbase = "$(in_fasta).fcov_$(coverage)x.art"
        end
        full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
            --paired \
            --errfree \
            --noALN \
            --seqSys $(seqSys) \
            --len $(read_length) \
            --mflen $(mflen) \
            --sdev $(sdev) \
            --rndSeed $(rndSeed) \
            --fcov $(coverage) \
            --in $(in_fasta) \
            --out $(outbase)`
    else
        error("Either 'coverage' or 'read_count' must be provided.")
    end

    forward_fq = "$(outbase)1.fq"
    forward_gz = forward_fq * ".gz"
    
    reverse_fq = "$(outbase)2.fq"
    reverse_gz = reverse_fq * ".gz"
    
    samfile = outbase * ".sam"
    samfile_gz = samfile * ".gz"
    
    error_free_samfile = outbase * "_errFree.sam"
    error_free_samfile_gz = error_free_samfile * ".gz"

    if !(isfile(forward_gz) && isfile(reverse_gz) && isfile(samfile_gz) && isfile(error_free_samfile_gz))
        @info "Running ART with command: $(full_cmd)"
        @time run(full_cmd)
    
        # Process output FASTQ files: gzip them if not already compressed.
        # For paired-end, output files are expected as outbase1.fq and outbase2.fq.
        # For single-end, only outbase1.fq is produced.
    
        @assert isfile(forward_fq) "Forward FASTQ file not found: $(forward_fq)"
        run(`gzip $(forward_fq)`)
        @assert isfile(forward_gz) "Gzipped forward FASTQ not found: $(forward_gz)"
    
        @assert isfile(reverse_fq) "Reverse FASTQ file not found: $(reverse_fq)"
        run(`gzip $(reverse_fq)`)
        @assert isfile(reverse_gz) "Gzipped reverse FASTQ not found: $(reverse_gz)"
    
        #  SAM Alignment File:
    
        @assert isfile(samfile)
        if !isfile(samfile_gz)
            run(`gzip $(samfile)`)
            @assert isfile(samfile_gz)
        end
        
        @assert isfile(error_free_samfile)
        if !isfile(error_free_samfile_gz)
            run(`gzip $(error_free_samfile)`)
            @assert isfile(error_free_samfile_gz)
        end
    else
        @info "All files already present, returning existing paths..."
    end

    return (forward_reads = forward_gz, reverse_reads = reverse_gz, sam = samfile_gz, error_free_sam = error_free_samfile_gz)
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

# src/matrix_generation.jl
"""
Generate a binary (Bernoulli) matrix with given dimensions and probability.

# Arguments
- `n_features::Int`: Number of features (rows)
- `n_samples::Int`: Number of samples (columns)
- `p::Float64`: Probability of 1 in the Bernoulli distribution

# Returns
- `Matrix{Bool}`: Binary matrix with dimensions (n_features, n_samples)
"""
function generate_binary_matrix(n_features::Int, n_samples::Int, p::Float64)
    return rand(Distributions.Bernoulli(p), n_features, n_samples)
end

"""
Generate a Poisson matrix with given dimensions and rate parameter.

# Arguments
- `n_features::Int`: Number of features (rows)
- `n_samples::Int`: Number of samples (columns)
- `λ::Float64`: Rate parameter for the Poisson distribution

# Returns
- `Matrix{Int}`: Poisson matrix with dimensions (n_features, n_samples)
"""
function generate_poisson_matrix(n_features::Int, n_samples::Int, λ::Float64)
    return rand(Distributions.Poisson(λ), n_features, n_samples)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulates genetic variants from sequences in a FASTA file and generates corresponding VCF records.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing sequences to analyze

# Details
1. Processes each record in the input FASTA file
2. Generates simulated variants for each sequence
3. Creates a VCF file with the same base name as input file (.vcf extension)
4. Updates sequences with simulated variants in a new FASTA file (.vcf.fna extension)

# Returns
Path to the modified FASTA file containing sequences with simulated variants
"""
function simulate_variants(fasta_file::String)
    vcf_file = fasta_file * ".vcf"
    modified_fasta_file = vcf_file * ".fna"
    vcf_table = DataFrames.DataFrame()
    for fasta_record in open_fastx(fasta_file)
        append!(vcf_table, simulate_variants(fasta_record))
    end
    Mycelia.write_vcf_table(vcf_file=vcf_file, vcf_table=vcf_table, fasta_file=fasta_file)
    return Mycelia.update_fasta_with_vcf(in_fasta = fasta_file, vcf_file = vcf_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulates genetic variants (substitutions, insertions, deletions, inversions) in a DNA sequence.

# Arguments
- `fasta_record`: Input DNA sequence in FASTA format

# Keywords
- `n_variants=√(sequence_length)`: Number of variants to generate
- `window_size=sequence_length/n_variants`: Size of windows for variant placement
- `variant_size_disbribution=Geometric(1/√window_size)`: Distribution for variant sizes
- `variant_type_likelihoods`: Vector of pairs mapping variant types to probabilities
    - `:substitution => 10⁻¹`
    - `:insertion => 10⁻²` 
    - `:deletion => 10⁻²`
    - `:inversion => 10⁻²`

# Returns
DataFrame in VCF format containing simulated variants with columns:
CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE

# Notes
- Variants are distributed across sequence windows to ensure spread
- Variant sizes are capped by window size
- Equivalent variants are filtered out
- FILTER column indicates variant type
"""
function simulate_variants(fasta_record::FASTX.FASTA.Record;        
        n_variants = Int(floor(sqrt(length(FASTX.sequence(fasta_record))))),
        window_size = Int(ceil(length(FASTX.sequence(fasta_record)) / n_variants)),
        variant_size_disbribution = Distributions.Geometric(1/sqrt(window_size)),
        variant_type_likelihoods = [
            :substitution => 10^-1,
            :insertion => 10^-2,
            :deletion => 10^-2,
            :inversion => 10^-2,
            # special case insertion/deletions, skipping
            # :translocations => 10^-3,
            # :duplication => 10^-3,
        ]
    )
    vcf_table = DataFrames.DataFrame(
        "#CHROM" => String[],
        "POS" => Int[],
        "ID" => String[],
        "REF" => String[],
        "ALT" => String[],
        "QUAL" => Int[],
        "FILTER" => String[],
        "INFO" => String[],
        "FORMAT" => String[],
        "SAMPLE" => String[]
    )
    @assert join(names(vcf_table), '\t') == "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE"
    
    original_sequence = BioSequences.LongDNA{4}(FASTX.sequence(fasta_record))
    original_sequence_id = string(first(split(FASTX.description(fasta_record))))
    modified_sequence_id = "modified__" * original_sequence_id
    
    variant_sizes = rand(variant_size_disbribution, n_variants) .+ 1
    @assert all(variant_sizes .>= 1)
    any(variant_sizes .>= window_size) && @warn "variant size distribution is too large, truncating variant sizes to fit in windows"
    variant_sizes = map(x -> min(x, window_size - 1), variant_sizes)
    # display(StatsPlots.plot(variant_size_disbribution, ylabel = "probability density", xlabel = "variant size", title = "Variant Size Distribution", legend=false))
    # display(StatsPlots.histogram(variant_sizes, ylabel="# of variants", xlabel="variant size", title="Actual samples drawn", nbins=length(unique(variant_sizes)), legend=false))
    variant_types = StatsBase.sample(first.(variant_type_likelihoods), StatsBase.weights(last.(variant_type_likelihoods)), n_variants)
    window_starts = 1:window_size:length(original_sequence)
    window_ends = window_size:window_size:length(original_sequence)
    windows = zip(window_starts, window_ends)
    variant_type_sizes = zip(variant_types, variant_sizes)
        
    for ((variant_type, variant_size), (start, stop)) in collect(zip(variant_type_sizes, windows))
        selected_start = rand(start:stop-variant_size)
        @assert selected_start <= stop
        original_subsequence = original_sequence[selected_start:selected_start+variant_size-1]
        @assert length(original_subsequence) == variant_size
        if variant_type == :substitution
            ref = original_subsequence
            alt = BioSequences.randdnaseq(variant_size)
            while alt == ref
                @info "substitution collision, resampling..."
                alt = BioSequences.randdnaseq(variant_size)
            end
        elseif variant_type == :insertion
            ref = original_subsequence
            alt = original_subsequence * BioSequences.randdnaseq(variant_size)
        elseif variant_type == :deletion
            ref = original_sequence[selected_start:selected_start+variant_size]
            alt = original_sequence[selected_start:selected_start]
        elseif variant_type == :inversion
            ref = original_subsequence
            alt = BioSequences.reverse_complement(original_subsequence)
        end 
        # @show selected_start
        row = Dict(
            "#CHROM" => original_sequence_id,
            "POS" => selected_start,
            "ID" => ".",
            "REF" => string(ref),
            "ALT" => string(alt),
            "QUAL" => 60,
            "FILTER" => string(variant_type),
            "INFO" => ".",
            "FORMAT" => "GT:GQ",
            "SAMPLE" => "1:60"
        )
        push!(vcf_table, row)
        original_prefix = original_sequence[start:selected_start-1]
        if variant_type == :deletion
            original_suffix = original_sequence[selected_start+variant_size+1:stop]
        else
            original_suffix = original_sequence[selected_start+variant_size:stop]
        end
        reconstructed_window = original_prefix * ref * original_suffix
        original_window = original_sequence[start:stop]
        @assert reconstructed_window == original_window
    end
    true_variant = vcf_table[!, "REF"] .!= vcf_table[!, "ALT"]
    if !all(true_variant)
        @warn "filtering equivalent variants"
        vcf_table = vcf_table[true_variant, :]
    end
    return vcf_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate sequencing of a DNA/RNA record by introducing random errors at the specified rate.

# Arguments
- `record`: A FASTA or FASTQ record containing the sequence to be "observed"
- `error_rate`: Probability of error at each position (default: 0.0)

# Returns
A new FASTQ.Record with:
- Random UUID as identifier
- Original record's description 
- Modified sequence with introduced errors
- Generated quality scores
"""
function observe(record::R; error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    converted_sequence = convert_sequence(FASTX.sequence(record))
    new_seq, quality_scores = observe(converted_sequence, error_rate=error_rate)
    new_seq_id = UUIDs.uuid4()
    new_seq_description = FASTX.identifier(record)
    return FASTX.FASTQ.Record(new_seq_id, new_seq_description, new_seq, quality_scores)
end



# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function observe(records::AbstractVector{R};
#                 weights=ones(length(records)),
#                 N = length(records),
#                 outfile = "",
#                 error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
#     if isempty(outfile)
#         error("no file name supplied")
#     end
#     io = open(outfile, "w")
#     fastx_io = FASTX.FASTA.Writer(io)
#     for i in 1:N
#         record = StatsBase.sample(records, StatsBase.weights(weights))
#         new_seq = observe(FASTX.sequence(record), error_rate=error_rate)
#         new_seq_id = Random.randstring(Int(ceil(log(length(new_seq) + 1))))
#         new_seq_description = FASTX.identifier(record)
#         observed_record = FASTX.FASTA.Record(new_seq_id, new_seq_description, new_seq)
#         write(fastx_io, observed_record)
#     end
#     close(fastx_io)
#     close(io)
#     return outfile
# end

# currently this is only for amino acid sequences, expand to include DNA and RNA via multiple dispatch
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a single random mutation in an amino acid sequence.

# Arguments
- `reference_sequence`: Input amino acid sequence to be mutated

# Returns
- `mutant_sequence`: The sequence after applying the mutation
- `haplotype`: A `SequenceVariation.Haplotype` object containing the mutation details

# Details
Performs one of three possible mutation types:
- Substitution: Replace one amino acid with another
- Insertion: Insert 1+ random amino acids at a position
- Deletion: Remove 1+ amino acids from a position

Insertion and deletion sizes follow a truncated Poisson distribution (λ=1, min=1).
"""
function mutate_sequence(reference_sequence)
    i = rand(1:length(reference_sequence))
    mutation_type = rand([SequenceVariation.Substitution, SequenceVariation.Insertion, SequenceVariation.Deletion])
    if mutation_type == SequenceVariation.Substitution
        # new_amino_acid = BioSequences.AA_A
        new_amino_acid = rand(amino_acids)
        mutation_string = "$(reference_sequence[i])$(i)$(new_amino_acid)"
    else
        indel_size = rand(Distributions.truncated(Distributions.Poisson(1), lower=1))
        if mutation_type == SequenceVariation.Insertion
            inserted_sequence = join([rand(amino_acids) for i in 1:indel_size], "")
            mutation_string = "$(i)$(inserted_sequence)"
        else
            @assert mutation_type == SequenceVariation.Deletion
            i_stop = min(i+indel_size, length(reference_sequence))
            mutation_string = "Δ$(i)-$(i_stop)"
        end
    end
    # println(mutation_string)
    mutation = SequenceVariation.Variation(reference_sequence, mutation_string)
    haplotype = SequenceVariation.Haplotype(reference_sequence, [mutation])
    mutant_sequence = SequenceVariation.reconstruct(haplotype)
    return mutant_sequence, haplotype
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function observe(records::AbstractVector{R}; outfile = "", error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
#     if isempty(outfile)
#         error("no file name supplied")
#     end
#     io = open(outfile, 'w')
#     fastx_io = FASTX.FASTQ.Writer(io)
#     for record in records
#         new_seq = observe(FASTX.sequence(record), error_rate=error_rate)
#         new_seq_id = string(hash(new_seq)) * "-" * Random.randstring(32)
#     new_seq_description = FASTX.identifier(record)
#     quality = fill(UInt8(60), length(new_seq))
#     return FASTX.FASTQ.Record(new_seq_id, new_seq_description, new_seq, quality)
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function random_fasta_record(;seed=rand(0:typemax(Int)), L = rand(0:Int(typemax(UInt16))))
#     id = Random.randstring(Int(ceil(log(L + 1))))
#     seq = BioSequences.randdnaseq(Random.seed!(seed), L)
#     return FASTX.FASTA.Record(id, seq)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generates a random FASTA record with a specified molecular type and sequence length.

# Arguments
- `moltype::Symbol=:DNA`: The type of molecule to generate (`:DNA`, `:RNA`, or `:AA` for amino acids).
- `seed`: The random seed used for sequence generation (default: a random integer).
- `L`: The length of the sequence (default: a random integer up to `typemax(UInt16)`).

# Returns
- A `FASTX.FASTA.Record` containing:
  - A randomly generated UUID identifier.
  - A randomly generated sequence of the specified type.

# Errors
- Throws an error if `moltype` is not one of `:DNA`, `:RNA`, or `:AA`.
"""
function random_fasta_record(;moltype::Symbol=:DNA, seed=rand(0:typemax(Int)), L = rand(0:Int(typemax(UInt16))))
    Random.seed!(seed)
    if moltype == :DNA
        seq = BioSequences.randdnaseq(StableRNGs.StableRNG(seed), L)
    elseif moltype == :RNA
        seq = BioSequences.randrnaseq(StableRNGs.StableRNG(seed), L)
    elseif moltype == :AA
        seq = BioSequences.randaaseq(StableRNGs.StableRNG(seed), L)
    else
        error("unrecognized molecule type: $(moltype) ! found in [:DNA, :RNA, :AA]")
    end
    # id = Mycelia.seq2sha256(seq)
    id = string(UUIDs.uuid4())
    return FASTX.FASTA.Record(id, seq)
end

"""
    observe(sequence::BioSequences.LongSequence{T}; error_rate=nothing, tech::Symbol=:illumina) where T

Simulates the “observation” of a biological polymer (DNA, RNA, or protein) by introducing realistic errors along with base‐quality scores.
The simulation takes into account both random and systematic error components. In particular, for technologies:
  
- **illumina**: (mostly substitution errors) the per‐base quality decays along the read (from ~Q40 at the start to ~Q20 at the end);
- **nanopore**: errors are more frequent and include both substitutions and indels (with overall lower quality scores, and an extra “homopolymer” penalty);
- **pacbio**: errors are dominated by indels (with quality scores typical of raw reads);
- **ultima**: (UG 100/ppmSeq™) correct bases are assigned very high quality (~Q60) while errors are extremely rare and, if they occur, are given a modest quality.

An error is introduced at each position with a (possibly position‐dependent) probability. For Illumina, the error probability increases along the read; additionally, if a base is part of a homopolymer run (length ≥ 3) and the chosen technology is one that struggles with homopolymers (nanopore, pacbio, ultima), then the local error probability is multiplied by a constant factor.

Returns a tuple `(new_seq, quality_scores)` where:
- `new_seq` is a `BioSequences.LongSequence{T}` containing the “observed” sequence (which may be longer or shorter than the input if insertions or deletions occur), and 
- `quality_scores` is a vector of integers representing the Phred quality scores (using the Sanger convention) for each base in the output sequence.
"""
function observe(sequence::BioSequences.LongSequence{T}; error_rate=nothing, tech::Symbol=:illumina) where T
    # Determine the appropriate alphabet based on the type T.
    if T <: BioSequences.DNAAlphabet
        alphabet = Mycelia.DNA_ALPHABET
    elseif T <: BioSequences.RNAAlphabet
        alphabet = Mycelia.RNA_ALPHABET
    else
        @assert T <: BioSequences.AminoAcidAlphabet "For amino acid sequences, T must be a subtype of BioSequences.AminoAcidAlphabet."
        alphabet = Mycelia.AA_ALPHABET
    end

    # Set a default baseline error rate if not provided.
    base_error_rate = isnothing(error_rate) ?
         (tech == :illumina  ? 0.005 :
          tech == :nanopore  ? 0.10  :
          tech == :pacbio    ? 0.11  :
          tech == :ultima    ? 1e-6  : error("Unknown technology")) : error_rate

    # Define error type probabilities (mismatch, insertion, deletion) for each technology.
    error_probs = Dict{Symbol, Float64}()
    if tech == :illumina
        error_probs[:mismatch]  = 0.90
        error_probs[:insertion] = 0.05
        error_probs[:deletion]  = 0.05
    elseif tech == :nanopore
        error_probs[:mismatch]  = 0.40
        error_probs[:insertion] = 0.30
        error_probs[:deletion]  = 0.30
    elseif tech == :pacbio
        error_probs[:mismatch]  = 0.20
        error_probs[:insertion] = 0.40
        error_probs[:deletion]  = 0.40
    elseif tech == :ultima
        error_probs[:mismatch]  = 0.95
        error_probs[:insertion] = 0.025
        error_probs[:deletion]  = 0.025
    end

    # Parameters for boosting error probability in homopolymer regions.
    homopolymer_threshold = 3    # If the homopolymer run length is ≥ 3, boost error probability.
    homopolymer_factor    = 2.0

    new_seq = BioSequences.LongSequence{T}()
    quality_scores = Int[]
    n = length(sequence)
    pos = 1

    while pos <= n
        base = sequence[pos]
        # For Illumina, increase error probability along the read; otherwise use a constant baseline.
        pos_error_rate = (tech == :illumina) ?
            base_error_rate * (1 + (pos - 1) / (n - 1)) : base_error_rate

        # Check for a homopolymer run by looking backwards.
        run_length = 1
        if pos > 1 && sequence[pos] == sequence[pos - 1]
            run_length = 2
            j = pos - 2
            while j ≥ 1 && sequence[j] == base
                run_length += 1
                j -= 1
            end
        end
        if run_length ≥ homopolymer_threshold && (tech in (:nanopore, :pacbio, :ultima))
            pos_error_rate *= homopolymer_factor
        end

        # Decide whether to observe the base correctly or introduce an error.
        if rand() > pos_error_rate
            # No error: add the correct base.
            push!(new_seq, base)
            push!(quality_scores, get_correct_quality(tech, pos, n))
            pos += 1
        else
            # An error occurs; choose the error type by sampling.
            r = rand()
            if r < error_probs[:mismatch]
                # Mismatch: choose a random base different from the true base.
                new_base = rand(filter(x -> x != base, alphabet))
                push!(new_seq, new_base)
                push!(quality_scores, get_error_quality(tech))
                pos += 1
            elseif r < (error_probs[:mismatch] + error_probs[:insertion])
                # Insertion: insert one or more random bases (simulate an extra insertion error),
                # then add the correct base.
                num_insertions = 1 + rand(Distributions.Poisson(pos_error_rate))
                for _ in 1:num_insertions
                    ins_base = rand(alphabet)
                    push!(new_seq, ins_base)
                    push!(quality_scores, get_error_quality(tech))
                end
                # Append the original base as a correct base.
                push!(new_seq, base)
                push!(quality_scores, get_correct_quality(tech, pos, n))
                pos += 1
            else
                # Deletion: skip adding the base (simulate a deletion).
                pos += 1
            end
        end
    end
    return new_seq, quality_scores
end

"""
    get_correct_quality(tech::Symbol, pos::Int, read_length::Int) -> Int

Simulates a Phred quality score (using the Sanger convention) for a correctly observed base.
For Illumina, the quality score is modeled to decay linearly from ~40 at the start to ~20 at the end of the read.
For other technologies, the score is sampled from a normal distribution with parameters typical for that platform.

Returns an integer quality score.
"""
function get_correct_quality(tech::Symbol, pos::Int, read_length::Int)
    if tech == :illumina
        q = 40 - 20 * (pos - 1) / (read_length - 1)
        q = clamp(round(Int, rand(Distributions.Normal(q, 2))), 20, 40)
        return q
    elseif tech == :nanopore
        q = clamp(round(Int, rand(Distributions.Normal(12, 2))), 10, 15)
        return q
    elseif tech == :pacbio
        q = clamp(round(Int, rand(Distributions.Normal(15, 2))), 12, 18)
        return q
    elseif tech == :ultima
        q = clamp(round(Int, rand(Distributions.Normal(60, 3))), 55, 65)
        return q
    else
        return 30
    end
end

"""
    get_error_quality(tech::Symbol) -> Int

Simulates a Phred quality score (using the Sanger convention) for a base observed with an error.
Error bases are assigned lower quality scores than correctly observed bases.
For Illumina, scores typically range between 5 and 15; for nanopore and pacbio, slightly lower values are used;
and for ultima, a modest quality score is assigned.

Returns an integer quality score.
"""
function get_error_quality(tech::Symbol)
    if tech == :illumina
        q = clamp(round(Int, rand(Distributions.Normal(10, 2))), 5, 15)
        return q
    elseif tech == :nanopore
        q = clamp(round(Int, rand(Distributions.Normal(7, 2))), 5, 10)
        return q
    elseif tech == :pacbio
        q = clamp(round(Int, rand(Distributions.Normal(7, 2))), 5, 10)
        return q
    elseif tech == :ultima
        q = clamp(round(Int, rand(Distributions.Normal(20, 3))), 15, 25)
        return q
    else
        return 10
    end
end