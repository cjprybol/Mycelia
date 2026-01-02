const PRINTABLE_ASCII_ALPHABET = filter(isprint, [Char(c) for c in 0x21:0x7E])
const PRINTABLE_GREEK_ALPHABET = vcat(
    filter(isprint, [Char(c) for c in 0x391:0x3A9]),  # Capital Greek letters
    filter(isprint, [Char(c) for c in 0x3B1:0x3C9])   # Lowercase Greek letters
)
const PRINTABLE_LATIN1_ALPHABET = filter(isprint, [Char(c) for c in 0x00:0xFF])
const PRINTABLE_UNICODE_ALPHABET = [Char(c) for c in 0x0000:0x10FFFF if ((c <= 0xD7FF) || (0xE000 <= c <= 0x10FFFF)) && Base.isprint(Char(c))]
const PRINTABLE_BMP_ALPHABET  = [Char(c) for c in 0x0020:0xFFFD if (c < 0xD800 || c > 0xDFFF) && Base.isprint(Char(c))]

"""
    mutate_string(s::String; alphabet::Union{Nothing,AbstractVector{Char}}=nothing, error_rate::Float64=0.01)

Introduce random mutations (substitutions, insertions, deletions) into a string.

# Arguments
- `s::String`: The input string to mutate

# Keywords
- `alphabet::Union{Nothing,AbstractVector{Char}}=nothing`: Characters to use for mutations. 
  If `nothing`, uses the unique characters found in the input string.
- `error_rate::Float64=0.01`: Probability of mutation at each position (default: 1%)

# Returns
- `String`: Mutated version of the input string

# Mutation Types
- **Substitution**: Replace character with a different one from the alphabet
- **Insertion**: Insert a random character from the alphabet
- **Deletion**: Remove the character (if string length > 1)

Each mutation type has equal probability when a mutation occurs.

# Example
```julia
# Mutate a DNA sequence with 2% error rate
mutated = mutate_string("ACGTACGT", error_rate=0.02)

# Mutate with custom alphabet
mutated = mutate_string("HELLO", alphabet=['H','E','L','O','W','R','D'], error_rate=0.1)
```

# Usage Context
Used in assembly testing to simulate sequencing errors and evaluate assembly algorithms'
robustness to different error rates and types.
"""
function mutate_string(s::String; alphabet::Union{Nothing,AbstractVector{Char}}=nothing, error_rate::Float64=0.01)
    chars = collect(s)
    if alphabet === nothing
        alphabet = unique(chars)
    end
    i = 1
    while i <= length(chars)
        if Random.rand() < error_rate
            mutation = Random.rand(["sub", "ins", "del"])
            if mutation == "sub"
                # Substitute with random char from alphabet (not the same as current)
                choices = setdiff(alphabet, [chars[i]])
                if !isempty(choices)
                    chars[i] = Random.rand(choices)
                end
            elseif mutation == "ins"
                # Insert random char from alphabet
                insert_char = Random.rand(alphabet)
                insert!(chars, i, insert_char)
                i += 1 # skip inserted char
            elseif mutation == "del" && length(chars) > 1
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
    mutate_dna_substitution_fraction(seq::Union{String,BioSequences.LongDNA{4}}; fraction::Float64, rng::Random.AbstractRNG=Random.default_rng())

Introduce substitutions at approximately `fraction` of positions (no indels) for DNA sequences.
Returns the mutated sequence in the same type as provided.
"""
function mutate_dna_substitution_fraction(seq::Union{String,BioSequences.LongDNA{4}}; fraction::Float64, rng::Random.AbstractRNG=Random.default_rng())
    if !(0.0 ≤ fraction ≤ 1.0)
        error("fraction must be between 0 and 1, got $(fraction)")
    end
    chars = seq isa String ? collect(seq) : collect(String(seq))
    n_mut = ceil(Int, fraction * length(chars))
    idxs = Random.randperm(rng, length(chars))[1:n_mut]
    for idx in idxs
        current = chars[idx]
        choices = setdiff(['A','C','G','T'], [current])
        chars[idx] = rand(rng, choices)
    end
    mutated = String(chars)
    return seq isa String ? mutated : BioSequences.LongDNA{4}(mutated)
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
- `seqSys::String`: Illumina sequencing system ID (default: "HS25"). Available ART profiles:
  - `"GA1"`: Genome Analyzer I (max length: 36bp or 44bp)
  - `"GA2"`: Genome Analyzer II (max length: 50bp or 75bp)  
  - `"HS10"`: HiSeq 1000 (max length: 100bp)
  - `"HS20"`: HiSeq 2000 (max length: 100bp)
  - `"HS25"`: HiSeq 2500 (max length: 125bp or 150bp)
  - `"HSXn"`: HiSeqX v2.5 PCR free (max length: 150bp)
  - `"HSXt"`: HiSeqX v2.5 TruSeq (max length: 150bp)
  - `"MSv1"`: MiSeq v1 (max length: 250bp)
  - `"MSv3"`: MiSeq v3 (max length: 250bp)
  - `"MinS"`: MiniSeq TruSeq (max length: 50bp)
  - `"NS50"`: NextSeq 500 v2 (max length: 75bp)
- `paired::Bool`: Whether to simulate paired-end reads (default: true). Some systems like MinS only support single-end.
- `amplicon::Bool`: Enable amplicon sequencing simulation mode (default: false).
- `errfree::Bool`: Generate an extra SAM file with zero sequencing errors (default: false).
- `rndSeed::Union{Nothing,Int}`: Optional seed for reproducibility (default: nothing).

# Outputs
Generates gzipped FASTQ files in the working directory:
- For paired-end: `\$(outbase)1.fq.gz` (forward) and `\$(outbase)2.fq.gz` (reverse).
- For single-end: `\$(outbase)1.fq.gz`.

SAM files are produced via ART's `--samout` option, and additional error-free
SAM output is created when `--errfree` is enabled.

# Details
This function calls ART with the provided options. Note that if `read_count` is supplied,
the function uses the `--rcount` option; otherwise, it uses `--fcov` with the given coverage.
Amplicon mode (via `--amplicon`) restricts the simulation to the amplicon regions, which is
important for targeted sequencing studies.

# Dependencies
Requires ART simulator (installed via Bioconda) and the Mycelia environment helper.

See also: `simulate_nanopore_reads`, `simulate_nearly_perfect_long_reads`, `simulate_pacbio_reads`
"""
function simulate_illumina_reads(;fasta::String,
    coverage::Union{Nothing, Number}=nothing,
    read_count::Union{Nothing, Int}=nothing,
    outbase::String = "",
    read_length::Int = 150,
    mflen::Int = 500,
    sdev::Int = 10,
    seqSys::String = "HS25",
    amplicon::Bool = false,
    errfree::Bool = false,
    paired::Bool = true,
    rndSeed::Int = current_unix_datetime(),
    quiet::Bool = true
    )

    # Ensure ART is available via Bioconda
    Mycelia.add_bioconda_env("art", quiet=quiet)
    
    errfree_flag = errfree ? ["--errfree"] : String[]
    # Determine read count or coverage
    if read_count !== nothing
        if isempty(outbase)
            outbase = "$(fasta).rcount_$(read_count).art"
        end
        # Build command based on pairing mode
        if paired
            full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
                --paired \
                --samout \
                --seqSys $(seqSys) \
                --len $(read_length) \
                --mflen $(mflen) \
                --sdev $(sdev) \
                --rndSeed $(rndSeed) \
                $(errfree_flag...) \
                --rcount $(read_count) \
                --in $(fasta) \
                --out $(outbase)`
        else
            full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
                --samout \
                --seqSys $(seqSys) \
                --len $(read_length) \
                --rndSeed $(rndSeed) \
                $(errfree_flag...) \
                --rcount $(read_count) \
                --in $(fasta) \
                --out $(outbase)`
        end
    elseif coverage !== nothing
        if isempty(outbase)
            outbase = "$(fasta).fcov_$(coverage)x.art"
        end
        # Build command based on pairing mode
        if paired
            full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
                --paired \
                --samout \
                --seqSys $(seqSys) \
                --len $(read_length) \
                --mflen $(mflen) \
                --sdev $(sdev) \
                --rndSeed $(rndSeed) \
                $(errfree_flag...) \
                --fcov $(coverage) \
                --in $(fasta) \
                --out $(outbase)`
        else
            full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
                --samout \
                --seqSys $(seqSys) \
                --len $(read_length) \
                --rndSeed $(rndSeed) \
                $(errfree_flag...) \
                --fcov $(coverage) \
                --in $(fasta) \
                --out $(outbase)`
        end
    else
        error("Either `coverage` or `read_count` must be provided.")
    end

    # For single-end, ART outputs to outbase.fq, for paired-end it's outbase1.fq  
    forward_fq = paired ? "$(outbase)1.fq" : "$(outbase).fq"
    forward_gz = forward_fq * ".gz"
    
    # For paired-end, we expect a reverse file
    reverse_fq = paired ? "$(outbase)2.fq" : ""
    reverse_gz = paired ? (reverse_fq * ".gz") : ""
    
    samfile = outbase * ".sam"
    samfile_gz = samfile * ".gz"
    
    error_free_samfile = outbase * "_errFree.sam"
    error_free_samfile_gz = error_free_samfile * ".gz"

    # Check if all expected files exist
    expected_files = [forward_gz, samfile_gz]
    if paired
        push!(expected_files, reverse_gz)
    end
    if errfree
        push!(expected_files, error_free_samfile_gz)
    end
    if !all(isfile.(expected_files))
        if !quiet
            @info "Running ART with command: $(full_cmd)"
        end
        if quiet
            run(pipeline(full_cmd, stdout=devnull, stderr=devnull))
        else
            @time run(full_cmd)
        end
    
        # Process output FASTQ files: gzip them if not already compressed.
        # For paired-end, output files are expected as outbase1.fq and outbase2.fq.
        # For single-end, only outbase1.fq is produced.
    
        @assert isfile(forward_fq) "Forward FASTQ file not found: $(forward_fq)"
        run(`gzip $(forward_fq)`)
        @assert isfile(forward_gz) "Gzipped forward FASTQ not found: $(forward_gz)"
    
        # Only process reverse reads if paired-end
        if paired
            @assert isfile(reverse_fq) "Reverse FASTQ file not found: $(reverse_fq)"
            run(`gzip $(reverse_fq)`)
            @assert isfile(reverse_gz) "Gzipped reverse FASTQ not found: $(reverse_gz)"
        end
    
        #  SAM Alignment File:
    
        @assert isfile(samfile)
        if !isfile(samfile_gz)
            run(`gzip $(samfile)`)
            @assert isfile(samfile_gz)
        end
        
        if errfree
            @assert isfile(error_free_samfile)
            if !isfile(error_free_samfile_gz)
                run(`gzip $(error_free_samfile)`)
                @assert isfile(error_free_samfile_gz)
            end
        end
    else
        @info "All files already present, returning existing paths..."
    end

    return (forward_reads = forward_gz, reverse_reads = paired ? reverse_gz : nothing, sam = samfile_gz, error_free_sam = errfree ? error_free_samfile_gz : nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Sample abundance weights for a community.

# Arguments
- `n_organisms::Int`: Number of organisms in the community.
- `balance::Symbol`: `:equal`, `:random`, or `:log_normal`.

# Keywords
- `rng::Random.AbstractRNG=StableRNGs.StableRNG(1)`: RNG for reproducibility.
- `lognorm_mu::Float64=0.0`, `lognorm_sigma::Float64=1.0`: Log-normal parameters.

# Returns
Vector of weights summing to 1.0.
"""
function sample_abundance_weights(;n_organisms::Int, balance::Symbol, rng::Random.AbstractRNG=StableRNGs.StableRNG(1), lognorm_mu::Float64=0.0, lognorm_sigma::Float64=1.0)
    @assert n_organisms > 0 "n_organisms must be positive"
    @assert balance in (:equal, :random, :log_normal)
    if balance == :equal
        return fill(1.0 / n_organisms, n_organisms)
    elseif balance == :random
        return rand(rng, Distributions.Dirichlet(fill(1.0, n_organisms)))
    else
        raw = rand(rng, Distributions.LogNormal(lognorm_mu, lognorm_sigma), n_organisms)
        return raw ./ sum(raw)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate a metagenomic community by sampling reference sequences, applying an
abundance profile, and generating reads per organism before merging.

# Arguments
- `reference_fasta::String`: Reference FASTA with all candidate sequences.
- `reference_table::DataFrames.DataFrame`: Table containing `sequence_id` (or `accession`) and taxonomy columns.
- `n_organisms::Int`: Number of organisms to sample (diversity).
- `depth_target::Real`: Mean per-genome coverage target.
- `abundance_profile::Symbol`: `:equal`, `:random`, `:log_normal`, or `:custom`.
- `readset::Symbol`: `:illumina_pe150`, `:illumina_se250`, `:nanopore`, or `:pacbio_hifi`.
- `outdir::String`: Output directory for reads and truth tables.

# Keywords
- `selected_ids::Union{Nothing,Vector{String}}=nothing`: Optional preselected IDs.
- `weights::Union{Nothing,Vector{Float64}}=nothing`: Optional precomputed abundance weights.
- `rng::Random.AbstractRNG=StableRNGs.StableRNG(1)`: RNG for reproducibility.
- `replicate::Int=1`: Replicate identifier.
- `run_simulation::Bool=true`: Skip external tool execution when false.
- `emit_truth_reads::Bool=true`: For Illumina readsets, also generate error-free reads.
- `lognorm_mu::Float64=0.0`, `lognorm_sigma::Float64=1.0`: Log-normal parameters.
- `mflen::Int=500`, `sdev::Int=10`: ART fragment length parameters.
- `quiet::Bool=true`: Suppress external tool output where supported.
- `cleanup::Bool=true`: Remove temporary per-genome FASTA/FASTQ files.

# Returns
Named tuple with `truth_table`, `reads`, `truth_reads`, `per_genome_fastas`, and `per_genome_reads`.
"""
function simulate_metagenome_community(;
    reference_fasta::String,
    reference_table::DataFrames.DataFrame,
    n_organisms::Int,
    depth_target::Real,
    abundance_profile::Symbol,
    readset::Symbol,
    outdir::String,
    selected_ids::Union{Nothing,Vector{String}}=nothing,
    weights::Union{Nothing,Vector{Float64}}=nothing,
    rng::Random.AbstractRNG=StableRNGs.StableRNG(1),
    replicate::Int=1,
    run_simulation::Bool=true,
    emit_truth_reads::Bool=true,
    lognorm_mu::Float64=0.0,
    lognorm_sigma::Float64=1.0,
    mflen::Int=500,
    sdev::Int=10,
    quiet::Bool=true,
    cleanup::Bool=true
)
    @assert isfile(reference_fasta) "Reference FASTA not found: $(reference_fasta)"
    @assert n_organisms > 0 "n_organisms must be positive"
    @assert depth_target > 0 "depth_target must be positive"
    @assert abundance_profile in (:equal, :random, :log_normal, :custom)
    @assert readset in (:illumina_pe150, :illumina_se250, :nanopore, :pacbio_hifi)
    mkpath(outdir)

    id_col = if "sequence_id" in DataFrames.names(reference_table)
        "sequence_id"
    elseif "accession" in DataFrames.names(reference_table)
        "accession"
    else
        error("reference_table must contain sequence_id or accession column")
    end
    id_col_sym = Symbol(id_col)
    length_col = "length" in DataFrames.names(reference_table) ? "length" : nothing
    length_col_sym = isnothing(length_col) ? nothing : Symbol(length_col)
    valid_rows = DataFrames.filter(row -> !ismissing(row[id_col_sym]), reference_table)
    available_ids = String.(valid_rows[!, id_col_sym])
    @assert length(available_ids) >= n_organisms "Not enough sequences to sample from reference_table"

    if isnothing(selected_ids)
        selected_ids = StatsBase.sample(rng, available_ids, n_organisms; replace=false)
    else
        @assert length(selected_ids) == n_organisms "selected_ids must match n_organisms"
    end
    if isnothing(weights)
        if abundance_profile == :custom
            error("abundance_profile=:custom requires explicit weights")
        else
            weights = Mycelia.sample_abundance_weights(
                n_organisms=n_organisms,
                balance=abundance_profile,
                rng=rng,
                lognorm_mu=lognorm_mu,
                lognorm_sigma=lognorm_sigma
            )
        end
    else
        @assert length(weights) == n_organisms "weights must match n_organisms"
        weights = weights ./ sum(weights)
    end
    coverages = depth_target * n_organisms .* weights

    selected_df = DataFrames.leftjoin(
        DataFrames.DataFrame(sequence_id = selected_ids, selection_order = collect(1:n_organisms)),
        valid_rows,
        on = [:sequence_id => id_col_sym],
        makeunique=true
    )
    DataFrames.sort!(selected_df, :selection_order)
    selected_df[!, :weight] = weights
    selected_df[!, :coverage] = coverages
    if !isnothing(length_col_sym)
        lengths = Float64.(selected_df[!, length_col_sym])
        base_weights = coverages .* lengths
        selected_df[!, :expected_relative_bases] = base_weights ./ sum(base_weights)
    end
    selected_df[!, :replicate] .= replicate
    truth_table = selected_df

    tmpdir = joinpath(outdir, "tmp")
    mkpath(tmpdir)
    fasta_dir = joinpath(tmpdir, "genomes")
    mkpath(fasta_dir)
    reads_dir = joinpath(outdir, "reads")
    mkpath(reads_dir)
    truth_reads_dir = joinpath(outdir, "truth_reads")
    mkpath(truth_reads_dir)

    selected_set = Set(selected_ids)
    per_genome_fastas = Dict{String,String}()
    for record in Mycelia.open_fastx(reference_fasta)
        record_id = String(FASTX.identifier(record))
        if record_id in selected_set
            safe_id = Mycelia.sanitize_fastx_identifier(record_id)
            fasta_path = joinpath(fasta_dir, "$(safe_id).fna")
            Mycelia.write_fasta(outfile=fasta_path, records=[record], gzip=false)
            per_genome_fastas[record_id] = fasta_path
        end
    end
    missing_ids = setdiff(selected_set, Set(keys(per_genome_fastas)))
    @assert isempty(missing_ids) "Missing sequences in reference_fasta: $(collect(missing_ids)[1:min(end, 5)])"

    per_genome_reads = DataFrames.DataFrame(
        sequence_id = String[],
        coverage = Float64[],
        readset = String[],
        forward_reads = Union{String,Missing}[],
        reverse_reads = Union{String,Missing}[],
        truth_forward_reads = Union{String,Missing}[],
        truth_reverse_reads = Union{String,Missing}[],
        sam = Union{String,Missing}[],
        truth_sam = Union{String,Missing}[]
    )

    merged_forward = nothing
    merged_reverse = nothing
    merged_truth_forward = nothing
    merged_truth_reverse = nothing

    if run_simulation
        forward_reads = String[]
        reverse_reads = String[]
        truth_forward_reads = String[]
        truth_reverse_reads = String[]

        for (idx, sequence_id) in enumerate(selected_ids)
            fasta = per_genome_fastas[sequence_id]
            coverage = coverages[idx]
            seed = rand(rng, 1:2^31-1)
            if readset == :illumina_pe150 || readset == :illumina_se250
                paired = readset == :illumina_pe150
                read_length = paired ? 150 : 250
                seq_sys = paired ? "HS25" : "MSv3"
                outbase = joinpath(tmpdir, "$(Mycelia.sanitize_fastx_identifier(sequence_id)).$(readset).r$(replicate)")
                sim = Mycelia.simulate_illumina_reads(
                    fasta=fasta,
                    coverage=coverage,
                    outbase=outbase,
                    read_length=read_length,
                    mflen=mflen,
                    sdev=sdev,
                    seqSys=seq_sys,
                    paired=paired,
                    rndSeed=seed,
                    errfree=false,
                    quiet=quiet
                )
                push!(forward_reads, sim.forward_reads)
                if paired && !isnothing(sim.reverse_reads)
                    push!(reverse_reads, sim.reverse_reads)
                end
                truth_forward = missing
                truth_reverse = missing
                truth_sam = missing
                if emit_truth_reads
                    truth_outbase = outbase * ".errfree"
                    truth_sim = Mycelia.simulate_illumina_reads(
                        fasta=fasta,
                        coverage=coverage,
                        outbase=truth_outbase,
                        read_length=read_length,
                        mflen=mflen,
                        sdev=sdev,
                        seqSys=seq_sys,
                        paired=paired,
                        rndSeed=seed,
                        errfree=true,
                        quiet=quiet
                    )
                    truth_forward = truth_sim.forward_reads
                    if paired && !isnothing(truth_sim.reverse_reads)
                        truth_reverse = truth_sim.reverse_reads
                    end
                    if !isnothing(truth_sim.error_free_sam)
                        truth_sam = truth_sim.error_free_sam
                    end
                    push!(truth_forward_reads, truth_forward)
                    if paired && truth_reverse !== missing
                        push!(truth_reverse_reads, truth_reverse)
                    end
                end
                DataFrames.push!(
                    per_genome_reads,
                    (sequence_id=sequence_id,
                     coverage=coverage,
                     readset=string(readset),
                     forward_reads=sim.forward_reads,
                     reverse_reads=paired ? sim.reverse_reads : missing,
                     truth_forward_reads=truth_forward,
                     truth_reverse_reads=paired ? truth_reverse : missing,
                     sam=sim.sam,
                     truth_sam=truth_sam)
                )
            elseif readset == :nanopore
                outbase = joinpath(tmpdir, "$(Mycelia.sanitize_fastx_identifier(sequence_id)).nanopore.r$(replicate).fq.gz")
                read_path = Mycelia.simulate_nanopore_reads(
                    fasta=fasta,
                    quantity="$(coverage)x",
                    outfile=outbase,
                    quiet=quiet,
                    seed=seed
                )
                push!(forward_reads, read_path)
                DataFrames.push!(
                    per_genome_reads,
                    (sequence_id=sequence_id,
                     coverage=coverage,
                     readset=string(readset),
                     forward_reads=read_path,
                     reverse_reads=missing,
                     truth_forward_reads=missing,
                     truth_reverse_reads=missing,
                     sam=missing,
                     truth_sam=missing)
                )
            else
                outbase = joinpath(tmpdir, "$(Mycelia.sanitize_fastx_identifier(sequence_id)).pacbio_hifi.r$(replicate).fq.gz")
                read_path = Mycelia.simulate_pacbio_reads(
                    fasta=fasta,
                    quantity="$(coverage)x",
                    outfile=outbase,
                    quiet=quiet,
                    seed=seed
                )
                push!(forward_reads, read_path)
                DataFrames.push!(
                    per_genome_reads,
                    (sequence_id=sequence_id,
                     coverage=coverage,
                     readset=string(readset),
                     forward_reads=read_path,
                     reverse_reads=missing,
                     truth_forward_reads=missing,
                     truth_reverse_reads=missing,
                     sam=missing,
                     truth_sam=missing)
                )
            end
        end

        merged_forward = joinpath(reads_dir, "community.$(readset).r$(replicate).R1.fq.gz")
        merged_forward = Mycelia.concatenate_fastq_files(fastq_files=forward_reads, output_fastq=merged_forward, gzip=true, force=true)
        if readset == :illumina_pe150
            merged_reverse = joinpath(reads_dir, "community.$(readset).r$(replicate).R2.fq.gz")
            merged_reverse = Mycelia.concatenate_fastq_files(fastq_files=reverse_reads, output_fastq=merged_reverse, gzip=true, force=true)
        end

        if emit_truth_reads && (readset == :illumina_pe150 || readset == :illumina_se250)
            merged_truth_forward = joinpath(truth_reads_dir, "community.$(readset).r$(replicate).truth.R1.fq.gz")
            merged_truth_forward = Mycelia.concatenate_fastq_files(fastq_files=truth_forward_reads, output_fastq=merged_truth_forward, gzip=true, force=true)
            if readset == :illumina_pe150
                merged_truth_reverse = joinpath(truth_reads_dir, "community.$(readset).r$(replicate).truth.R2.fq.gz")
                merged_truth_reverse = Mycelia.concatenate_fastq_files(fastq_files=truth_reverse_reads, output_fastq=merged_truth_reverse, gzip=true, force=true)
            end
        end
    end

    if cleanup
        rm(tmpdir; force=true, recursive=true)
    end

    return (
        truth_table = truth_table,
        reads = (forward=merged_forward, reverse=merged_reverse),
        truth_reads = (forward=merged_truth_forward, reverse=merged_truth_reverse),
        per_genome_fastas = per_genome_fastas,
        per_genome_reads = per_genome_reads
    )
end

# Individual ART Illumina profile functions with specific parameters
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using GA1 Genome Analyzer I with 36bp read length.
"""
function simulate_illumina_ga1_36bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="GA1", read_length=36, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using GA1 Genome Analyzer I with 44bp read length.
"""
function simulate_illumina_ga1_44bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="GA1", read_length=44, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using GA2 Genome Analyzer II with 50bp read length.
"""
function simulate_illumina_ga2_50bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="GA2", read_length=50, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using GA2 Genome Analyzer II with 75bp read length.
"""
function simulate_illumina_ga2_75bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="GA2", read_length=75, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using HS10 HiSeq 1000 with 100bp read length.
"""
function simulate_illumina_hs10_100bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="HS10", read_length=100, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using HS20 HiSeq 2000 with 100bp read length.
"""
function simulate_illumina_hs20_100bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="HS20", read_length=100, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using HS25 HiSeq 2500 with 125bp read length.
"""
function simulate_illumina_hs25_125bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="HS25", read_length=125, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using HS25 HiSeq 2500 with 150bp read length.
"""
function simulate_illumina_hs25_150bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="HS25", read_length=150, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using HSXn HiSeqX v2.5 PCR free with 150bp read length.
"""
function simulate_illumina_hsxn_150bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="HSXn", read_length=150, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using HSXt HiSeqX v2.5 TruSeq with 150bp read length.
"""
function simulate_illumina_hsxt_150bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="HSXt", read_length=150, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using MSv1 MiSeq v1 with 250bp read length.
"""
function simulate_illumina_msv1_250bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="MSv1", read_length=250, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using MSv3 MiSeq v3 with 250bp read length.
"""
function simulate_illumina_msv3_250bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="MSv3", read_length=250, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using MinS MiniSeq TruSeq with 50bp read length.
"""
function simulate_illumina_mins_50bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="MinS", read_length=50, paired=false, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Illumina reads using NS50 NextSeq 500 v2 with 75bp read length.
"""
function simulate_illumina_ns50_75bp(;fasta::String, coverage::Union{Nothing,Number}=nothing, read_count::Union{Nothing,Int}=nothing, kwargs...)
    return simulate_illumina_reads(;fasta=fasta, coverage=coverage, read_count=read_count, seqSys="NS50", read_length=75, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Ultima Genomics-like reads using ART with custom error profiles.

Uses MiSeq v3 profile as a base with custom insertion/deletion rates to mimic Ultima Genomics sequencing characteristics.
The Ultima platform produces single-end reads with specific error patterns including higher indel rates compared to traditional Illumina systems.

# Arguments
- `fasta::String`: Path to input reference FASTA file
- `coverage::Union{Nothing,Number}=30`: Fold coverage (default: 30x as recommended for Ultima)
- `read_count::Union{Nothing,Int}=nothing`: Alternative to coverage - total number of reads
- `read_length::Int=250`: Read length in base pairs (default: 250bp for Ultima long reads)
- `insertion_rate::Float64=0.004`: Insertion error rate (4 per 1000 bases)
- `deletion_rate::Float64=0.004`: Deletion error rate (4 per 1000 bases)
- `id::String="sim_ultima"`: Read group identifier

# Returns
Named tuple with paths to generated files: (forward_reads, reverse_reads=nothing, sam, error_free_sam)

# Notes
- Produces single-end reads only (reverse_reads will be nothing)
- Uses MSv3 (MiSeq v3) profile as base since it supports longer reads up to 250bp
- Higher indel rates (0.4% each) reflect Ultima's characteristic error profile
- Generates both regular and error-free alignment files for validation

See also: `simulate_illumina_reads`, `simulate_illumina_ns50_75bp`
"""
function simulate_ultima_reads(;fasta::String, 
                              coverage::Union{Nothing,Number}=30,
                              read_count::Union{Nothing,Int}=nothing,
                              read_length::Int=250,
                              insertion_rate::Float64=0.004,
                              deletion_rate::Float64=0.004,
                              id::String="sim_ultima",
                              quiet::Bool=false,
                              kwargs...)
    
    # Ensure ART is available via Bioconda
    Mycelia.add_bioconda_env("art", quiet=quiet)
    
    # Generate output base name
    outbase = isempty(get(kwargs, :outbase, "")) ? "$(fasta).ultima_$(coverage !== nothing ? "$(coverage)x" : "rcount_$(read_count)").art" : kwargs[:outbase]
    
    # Build command for Ultima-like simulation (single-end only)
    if read_count !== nothing
        full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
            --errfree \
            --noALN \
            --seqSys MSv3 \
            --len $(read_length) \
            --rndSeed $(get(kwargs, :rndSeed, Mycelia.current_unix_datetime())) \
            --rcount $(read_count) \
            --id $(id) \
            -ir $(insertion_rate) \
            -dr $(deletion_rate) \
            --in $(fasta) \
            --out $(outbase)`
    elseif coverage !== nothing  
        full_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n art art_illumina \
            --errfree \
            --noALN \
            --seqSys MSv3 \
            --len $(read_length) \
            --rndSeed $(get(kwargs, :rndSeed, Mycelia.current_unix_datetime())) \
            --fcov $(coverage) \
            --id $(id) \
            -ir $(insertion_rate) \
            -dr $(deletion_rate) \
            --in $(fasta) \
            --out $(outbase)`
    else
        error("Either `coverage` or `read_count` must be provided.")
    end

    # Define expected output files (single-end for Ultima)
    forward_fq = "$(outbase).fq"
    forward_gz = forward_fq * ".gz"
    samfile = outbase * ".sam"
    samfile_gz = samfile * ".gz"
    error_free_samfile = outbase * "_errFree.sam"
    error_free_samfile_gz = error_free_samfile * ".gz"

    # Check if all expected files exist
    expected_files = [forward_gz, samfile_gz, error_free_samfile_gz]
    if !all(isfile.(expected_files))
        if !quiet
            @info "Running ART with command: $(full_cmd)"
        end
        if quiet
            run(pipeline(full_cmd, stdout=devnull, stderr=devnull))
        else
            @time run(full_cmd)
        end
    
        # Process output FASTQ file
        @assert isfile(forward_fq) "Forward FASTQ file not found: $(forward_fq)"
        run(`gzip $(forward_fq)`)
        @assert isfile(forward_gz) "Gzipped forward FASTQ not found: $(forward_gz)"
    
        # Process SAM files
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

    return (forward_reads = forward_gz, reverse_reads = nothing, sam = samfile_gz, error_free_sam = error_free_samfile_gz)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate PacBio HiFi reads using optimized Badread settings.

Uses PacBio 2021 error and quality score models with identity settings optimized for HiFi reads.
Follows the recommended PacBio HiFi simulation parameters from the Badread documentation.

# Arguments
- `fasta::String`: Path to input FASTA file containing reference sequence
- `quantity::String`: Coverage depth (e.g. "50x") or total bases (e.g. "1000000")
- `outfile::String`: Output filepath for simulated reads. Defaults to input filename with PacBio HiFi suffix
- `quiet::Bool=false`: If true, suppress badread output to stdout/stderr

# Returns
- `String`: Path to the generated output file

# Notes
- Uses pacbio2021 error and quality score models
- Sets identity to 30,3 for HiFi-appropriate accuracy
- Average read length ~15kb
- Output is gzipped FASTQ format

See also: `simulate_nanopore_reads`, `simulate_nearly_perfect_long_reads`, `simulate_badread_reads`
"""
function simulate_pacbio_reads(;fasta, quantity, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.pacbio_hifi.$(quantity).fq.gz"), quiet=false, seed::Union{Nothing,Int}=nothing)
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        cmd_args = ["badread", "simulate", "--error_model", "pacbio2021", "--qscore_model", "pacbio2021", "--identity", "30,3", "--reference", fasta, "--quantity", quantity]
        if !isnothing(seed)
            push!(cmd_args, "--seed")
            push!(cmd_args, string(seed))
        end
        if quiet
            cmd = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread $(cmd_args)`, stderr=devnull)
            p = pipeline(cmd, `gzip`)
            run(pipeline(p, outfile))
        else
            p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread $(cmd_args)`, `gzip`)
            run(pipeline(p, outfile))
        end
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate Oxford Nanopore R10.4.1 sequencing reads using Badread's default settings.

Badread's default settings correspond to Oxford Nanopore R10.4.1 reads of mediocre quality.
Uses nanopore2023 error and quality models with default identity and length distributions.

# Arguments
- `fasta::String`: Path to input reference FASTA file
- `quantity::String`: Either fold coverage (e.g. "50x") or total bases to sequence (e.g. "1000000")
- `outfile::String`: Output path for gzipped FASTQ file. Defaults to input filename with modified extension

# Returns
- `String`: Path to the generated output FASTQ file

See also: `simulate_pacbio_reads`, `simulate_nanopore_r941_reads`, `simulate_badread_reads`
"""
function simulate_nanopore_reads(;fasta, quantity, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.nanopore_r10.$(quantity).fq.gz"), quiet=false, seed::Union{Nothing,Int}=nothing)
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        cmd_args = ["badread", "simulate", "--reference", fasta, "--quantity", quantity]
        if !isnothing(seed)
            push!(cmd_args, "--seed")
            push!(cmd_args, string(seed))
        end
        if quiet
            cmd = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread $(cmd_args)`, stderr=devnull)
            p = pipeline(cmd, `gzip`)
            run(pipeline(p, outfile))
        else
            p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread $(cmd_args)`, `gzip`)
            run(pipeline(p, outfile))
        end
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate high-quality long reads with minimal errors using Badread.

# Arguments
- `fasta::String`: Path to reference FASTA file
- `quantity::String`: Coverage depth (e.g. "50x") or total bases (e.g. "1000000")
- `length_mean::Int=40000`: Mean read length
- `length_sd::Int=20000`: Standard deviation of read length
- `outfile::String`: Output filepath for simulated reads. Defaults to input filename with ".badread.perfect.\${quantity}.fq.gz" suffix

# Returns
- `String`: Path to the generated output file

# Details
Generates nearly perfect long reads by setting error rates and artifacts to minimum values.
Uses ideal quality scores and disables common sequencing artifacts like chimeras and adapters.

See also: `simulate_pacbio_reads`, `simulate_nanopore_reads`, `simulate_illumina_reads`
"""
function simulate_nearly_perfect_long_reads(;fasta, quantity, length_mean::Int=40000, length_sd::Int=20000, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.perfect.$(quantity).fq.gz"), quiet=false)
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        if quiet
            cmd = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate \
                --reference $(fasta) \
                --quantity $(quantity) \
                --error_model random \
                --qscore_model ideal \
                --glitches 0,0,0 \
                --junk_reads 0 \
                --random_reads 0 \
                --chimeras 0 \
                --identity 30,3 \
                --length $(length_mean),$(length_sd) \
                --start_adapter_seq "" \
                --end_adapter_seq ""`, stderr=devnull)
            p = pipeline(cmd, `gzip`)
            run(pipeline(p, outfile))
        else
            p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate \
                --reference $(fasta) \
                --quantity $(quantity) \
                --error_model random \
                --qscore_model ideal \
                --glitches 0,0,0 \
                --junk_reads 0 \
                --random_reads 0 \
                --chimeras 0 \
                --identity 30,3 \
                --length $(length_mean),$(length_sd) \
                --start_adapter_seq "" \
                --end_adapter_seq ""`, `gzip`)
            run(pipeline(p, outfile))
        end
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate older Oxford Nanopore R9.4.1 sequencing reads with worse basecalling.

Uses nanopore2020 error and quality models with identity settings that reflect the lower
accuracy of older nanopore sequencing technology and basecalling algorithms.

# Arguments
- `fasta::String`: Path to input reference FASTA file
- `quantity::String`: Either fold coverage (e.g. "50x") or total bases to sequence (e.g. "1000000")
- `outfile::String`: Output path for gzipped FASTQ file. Defaults to input filename with R9.4.1 suffix

# Returns
- `String`: Path to the generated output FASTQ file

# Notes
- Uses nanopore2020 error and quality models
- Sets identity to 90,98,5 for R9.4.1-appropriate lower accuracy
- Average read length ~15kb (default Badread length)
- Output is gzipped FASTQ format

See also: `simulate_nanopore_reads`, `simulate_pacbio_reads`, `simulate_badread_reads`
"""
function simulate_nanopore_r941_reads(;fasta, quantity, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.nanopore_r941.$(quantity).fq.gz"), quiet=false)
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        if quiet
            cmd = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --error_model nanopore2020 --qscore_model nanopore2020 --identity 90,98,5 --reference $(fasta) --quantity $(quantity)`, stderr=devnull)
            p = pipeline(cmd, `gzip`)
            run(pipeline(p, outfile))
        else
            p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --error_model nanopore2020 --qscore_model nanopore2020 --identity 90,98,5 --reference $(fasta) --quantity $(quantity)`, `gzip`)
            run(pipeline(p, outfile))
        end
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate very bad quality sequencing reads with high error rates and artifacts.

Simulates reads with high error rates, frequent glitches, junk reads, random reads, 
and chimeras to test assembly algorithms under challenging conditions.

# Arguments
- `fasta::String`: Path to input reference FASTA file
- `quantity::String`: Either fold coverage (e.g. "50x") or total bases to sequence (e.g. "1000000")
- `outfile::String`: Output path for gzipped FASTQ file. Defaults to input filename with "very_bad" suffix

# Returns
- `String`: Path to the generated output FASTQ file

# Notes
- High glitch rates: 1000,100,100
- 5% junk reads, 5% random reads, 10% chimeras
- Low identity: 80,90,6 (80% mean, 90% max, 6% stdev)
- Shorter reads: 4000±2000 bp
- Designed to stress-test assembly algorithms

See also: `simulate_pretty_good_reads`, `simulate_nanopore_reads`, `simulate_badread_reads`
"""
function simulate_very_bad_reads(;fasta, quantity, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.very_bad.$(quantity).fq.gz"), quiet=false)
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        if quiet
            cmd = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --glitches 1000,100,100 --junk_reads 5 --random_reads 5 --chimeras 10 --identity 80,90,6 --length 4000,2000 --reference $(fasta) --quantity $(quantity)`, stderr=devnull)
            p = pipeline(cmd, `gzip`)
            run(pipeline(p, outfile))
        else
            p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --glitches 1000,100,100 --junk_reads 5 --random_reads 5 --chimeras 10 --identity 80,90,6 --length 4000,2000 --reference $(fasta) --quantity $(quantity)`, `gzip`)
            run(pipeline(p, outfile))
        end
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate pretty good quality sequencing reads with low error rates and minimal artifacts.

Simulates high-quality reads with low error rates and minimal sequencing artifacts,
suitable for testing assembly algorithms under favorable conditions.

# Arguments
- `fasta::String`: Path to input reference FASTA file
- `quantity::String`: Either fold coverage (e.g. "50x") or total bases to sequence (e.g. "1000000")
- `outfile::String`: Output path for gzipped FASTQ file. Defaults to input filename with "pretty_good" suffix

# Returns
- `String`: Path to the generated output FASTQ file

# Notes
- Low glitch rates: 10000,10,10
- Minimal artifacts: 0.1% junk reads, 0.1% random reads, 0.1% chimeras
- High identity: 20,3 (mean 20%, stdev 3% - using qscore distribution)
- Default read length: ~15kb
- Designed for high-quality assembly testing

See also: `simulate_very_bad_reads`, `simulate_nearly_perfect_long_reads`, `simulate_badread_reads`
"""
function simulate_pretty_good_reads(;fasta, quantity, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.pretty_good.$(quantity).fq.gz"), quiet=false)
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        if quiet
            cmd = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --glitches 10000,10,10 --junk_reads 0.1 --random_reads 0.1 --chimeras 0.1 --identity 20,3 --reference $(fasta) --quantity $(quantity)`, stderr=devnull)
            p = pipeline(cmd, `gzip`)
            run(pipeline(p, outfile))
        else
            p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --glitches 10000,10,10 --junk_reads 0.1 --random_reads 0.1 --chimeras 0.1 --identity 20,3 --reference $(fasta) --quantity $(quantity)`, `gzip`)
            run(pipeline(p, outfile))
        end
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

General Badread simulator with full parameter control.

Provides access to all Badread simulation parameters for custom read simulation scenarios.
This function exposes the complete Badread parameter set for maximum flexibility.

# Arguments
- `fasta::String`: Path to input reference FASTA file  
- `quantity::String`: Either fold coverage (e.g. "50x") or total bases to sequence (e.g. "1000000")

# Keywords
- `length::String="15000,13000"`: Fragment length distribution (mean,stdev)
- `identity::String="95,99,2.5"`: Sequencing identity distribution (mean,max,stdev for beta or mean,stdev for normal qscore)
- `error_model::String="nanopore2023"`: Error model ("nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016", "pacbio2021", "random", or filename)
- `qscore_model::String="nanopore2023"`: Quality score model ("nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016", "pacbio2021", "random", "ideal", or filename)
- `seed::Union{Nothing,Int}=nothing`: Random seed for deterministic output
- `start_adapter::String="90,60"`: Adapter parameters for read starts (rate,amount)
- `end_adapter::String="50,20"`: Adapter parameters for read ends (rate,amount)  
- `start_adapter_seq::String="AATGTACTTCGTTCAGTTACGTATTGCT"`: Adapter sequence for read starts
- `end_adapter_seq::String="GCAATACGTAACTGAACGAAGT"`: Adapter sequence for read ends
- `junk_reads::Float64=1.0`: Percentage of low-complexity junk reads
- `random_reads::Float64=1.0`: Percentage of random sequence reads
- `chimeras::Float64=1.0`: Percentage at which separate fragments join
- `glitches::String="10000,25,25"`: Read glitch parameters (rate,size,skip)
- `small_plasmid_bias::Bool=false`: Lose small circular plasmids when fragment length is too high
- `outfile::String=""`: Output path for gzipped FASTQ file. Auto-generated if empty

# Returns
- `String`: Path to the generated output FASTQ file

# Examples
```julia
# Custom nanopore simulation with specific error rates
simulate_badread_reads(fasta="ref.fasta", quantity="50x", 
                      error_model="nanopore2020", identity="85,95,8")

# Custom PacBio-like simulation  
simulate_badread_reads(fasta="ref.fasta", quantity="25x",
                      error_model="pacbio2021", qscore_model="pacbio2021", 
                      identity="30,3", length="20000,15000")
```

See also: `simulate_nanopore_reads`, `simulate_pacbio_reads`, `simulate_very_bad_reads`
"""
function simulate_badread_reads(;fasta::String, quantity::String,
                               length::String="15000,13000",
                               identity::String="95,99,2.5", 
                               error_model::String="nanopore2023",
                               qscore_model::String="nanopore2023",
                               seed::Union{Nothing,Int}=nothing,
                               start_adapter::String="90,60",
                               end_adapter::String="50,20",
                               start_adapter_seq::String="AATGTACTTCGTTCAGTTACGTATTGCT",
                               end_adapter_seq::String="GCAATACGTAACTGAACGAAGT",
                               junk_reads::Float64=1.0,
                               random_reads::Float64=1.0,
                               chimeras::Float64=1.0,
                               glitches::String="10000,25,25",
                               small_plasmid_bias::Bool=false,
                               outfile::String="",
                               quiet::Bool=false)
    
    if isempty(outfile)
        outfile = replace(fasta, Mycelia.FASTA_REGEX => ".badread.custom.$(error_model).$(qscore_model).$(quantity).fq.gz")
    end
    
    if !isfile(outfile) || (filesize(outfile) == 0)
        Mycelia.add_bioconda_env("badread")
        
        # Build command with all parameters
        cmd_parts = [
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate`,
            `--reference $(fasta)`,
            `--quantity $(quantity)`,
            `--length $(length)`,
            `--identity $(identity)`,
            `--error_model $(error_model)`,
            `--qscore_model $(qscore_model)`,
            `--start_adapter $(start_adapter)`,
            `--end_adapter $(end_adapter)`,
            `--start_adapter_seq $(start_adapter_seq)`,
            `--end_adapter_seq $(end_adapter_seq)`,
            `--junk_reads $(junk_reads)`,
            `--random_reads $(random_reads)`,
            `--chimeras $(chimeras)`,
            `--glitches $(glitches)`
        ]
        
        # Add optional parameters
        if seed !== nothing
            push!(cmd_parts, `--seed $(seed)`)
        end
        
        if small_plasmid_bias
            push!(cmd_parts, `--small_plasmid_bias`)
        end
        
        # Combine all command parts
        full_cmd = reduce(((a, b) -> `$a $b`), cmd_parts)
        
        if quiet
            cmd = pipeline(full_cmd, stderr=devnull)
            p = pipeline(cmd, `gzip`)
            run(pipeline(p, outfile))
        else
            p = pipeline(full_cmd, `gzip`)
            run(pipeline(p, outfile))
        end
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
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
    new_seq_id = string(UUIDs.uuid4())
    new_seq_description = FASTX.identifier(record)
    # Convert quality scores to string (Phred+33 encoding)
    quality_string = String([Char(q + 33) for q in quality_scores])
    return FASTX.FASTQ.Record(new_seq_id, new_seq, quality_string)
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
    generate_paired_end_reads(reference_seq, coverage, read_length, insert_size; error_rate=0.01)

Generate realistic paired-end sequencing reads from a reference sequence.

Simulates paired-end Illumina sequencing with realistic insert sizes, read lengths, and
optional sequencing errors for assembly benchmarking.

# Arguments
- `reference_seq`: Reference sequence (BioSequences.LongDNA{4})
- `coverage`: Target sequencing coverage depth
- `read_length`: Length of each read in base pairs
- `insert_size`: Insert size between paired reads
- `error_rate`: Sequencing error rate (default: 0.01)

# Returns
- Tuple of (forward_reads, reverse_reads) as vectors of BioSequences.LongDNA{4}

# See Also
- `simulate_illumina_reads`: For more sophisticated read simulation using ART
- `introduce_sequencing_errors`: For adding realistic sequencing errors
"""
function generate_paired_end_reads(reference_seq, coverage, read_length, insert_size; error_rate=0.01)
    ref_length = length(reference_seq)
    n_reads = Int(round((coverage * ref_length) / (2 * read_length)))
    
    reads_1 = BioSequences.LongDNA{4}[]
    reads_2 = BioSequences.LongDNA{4}[]
    
    for i in 1:n_reads
        # Random start position ensuring we can get both reads
        start_pos = rand(1:(ref_length - insert_size + 1))
        
        # Extract forward read
        read1_seq = reference_seq[start_pos:(start_pos + read_length - 1)]
        
        # Extract reverse read (reverse complement from the other end)
        read2_start = start_pos + insert_size - read_length
        read2_seq = BioSequences.reverse_complement(reference_seq[read2_start:(read2_start + read_length - 1)])
        
        # Add sequencing errors
        if error_rate > 0
            read1_seq = introduce_sequencing_errors(read1_seq, error_rate)
            read2_seq = introduce_sequencing_errors(read2_seq, error_rate)
        end
        
        push!(reads_1, read1_seq)
        push!(reads_2, read2_seq)
    end
    
    return reads_1, reads_2
end

"""
    introduce_sequencing_errors(sequence, error_rate)

Introduce realistic sequencing errors into a DNA sequence.

Simulates substitutions (70%), insertions (15%), and deletions (15%) at the specified
error rate for realistic sequencing simulation.

# Arguments
- `sequence`: Input DNA sequence (BioSequences.LongDNA{4})
- `error_rate`: Probability of error per base

# Returns
- Modified sequence with introduced errors (BioSequences.LongDNA{4})

# See Also
- `observe`: For more sophisticated error modeling with quality scores
- `mutate_string`: For string-based mutation operations
"""
function introduce_sequencing_errors(sequence, error_rate)
    seq_array = collect(sequence)
    n_errors = Int(round(length(seq_array) * error_rate))
    
    for _ in 1:n_errors
        pos = rand(1:length(seq_array))
        # 70% substitutions, 15% insertions, 15% deletions
        error_type = rand()
        if error_type < 0.7
            # Substitution
            seq_array[pos] = rand([BioSequences.DNA_A, BioSequences.DNA_C, BioSequences.DNA_G, BioSequences.DNA_T])
        elseif error_type < 0.85 && length(seq_array) < 1000  # Limit insertions to prevent sequence explosion
            # Insertion (insert random base)
            insert!(seq_array, pos, rand([BioSequences.DNA_A, BioSequences.DNA_C, BioSequences.DNA_G, BioSequences.DNA_T]))
        elseif length(seq_array) > 1
            # Deletion
            deleteat!(seq_array, pos)
        end
    end
    
    return BioSequences.LongDNA{4}(seq_array)
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
