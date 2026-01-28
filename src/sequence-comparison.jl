struct SylphProfileResult
    syldb::String
    sample_sketches::Vector{String}
    output_tsv::String
    table::DataFrames.DataFrame
end

"""
    run_sylph_profile(reference_fastas::Vector{String};
                      sample_reads::Vector{String}=String[],
                      first_pairs::Vector{String}=String[],
                      second_pairs::Vector{String}=String[],
                      threads::Int=get_default_threads(),
                      subsampling::Int=200,
                      k::Int=31,
                      min_spacing::Int=30,
                      min_ani::Float64=95.0,
                      min_kmers::Int=50,
                      estimate_unknown::Bool=false,
                      outdir::Union{Nothing,String}=nothing,
                      output_prefix::String="sylph_db",
                      output_tsv::Union{Nothing,String}=nothing,
                      additional_args::Vector{String}=String[],
                      quiet::Bool=true)

Sketch references and samples with Sylph and run `sylph profile`, returning the parsed TSV.

# Arguments
- `reference_fastas`: Reference FASTA files to build the Sylph database.
- `sample_reads`: Single-end read files (FASTQ/FASTA, gz accepted).
- `first_pairs`/`second_pairs`: Paired-end read files. Lengths must match.
- `threads`: Thread count for sketching/profile.
- `subsampling`: Sylph `-c` subsampling rate (default 200).
- `k`: k-mer size (Sylph supports 21 or 31).
- `min_spacing`: Minimum spacing between sampled k-mers.
- `min_ani`: Minimum adjusted ANI threshold for reporting (passed via `-m`).
- `min_kmers`: Minimum sampled k-mers required (`-M`).
- `estimate_unknown`: Pass `-u` to estimate unknown content.
- `outdir`: Output directory (created if missing). Defaults to `mktempdir()`.
- `output_prefix`: Prefix for generated `.syldb`.
- `output_tsv`: Optional explicit path for profile output TSV.
- `additional_args`: Extra CLI args appended to `sylph profile`.
- `quiet`: Suppress Sylph stdout/stderr when true.

# Returns
`SylphProfileResult` containing paths to the database, sample sketches, TSV, and parsed DataFrame.
"""
function run_sylph_profile(reference_fastas::Vector{String};
        sample_reads::Vector{String}=String[],
        first_pairs::Vector{String}=String[],
        second_pairs::Vector{String}=String[],
        threads::Int=get_default_threads(),
        subsampling::Int=200,
        k::Int=31,
        min_spacing::Int=30,
        min_ani::Float64=95.0,
        min_kmers::Int=50,
        estimate_unknown::Bool=false,
        outdir::Union{Nothing,String}=nothing,
        output_prefix::String="sylph_db",
        output_tsv::Union{Nothing,String}=nothing,
        additional_args::Vector{String}=String[],
        quiet::Bool=true)

    if length(first_pairs) != length(second_pairs)
        error("first_pairs and second_pairs must have the same length")
    end
    reference_fastas_abs = abspath.(reference_fastas)
    sample_reads_abs = abspath.(sample_reads)
    first_pairs_abs = abspath.(first_pairs)
    second_pairs_abs = abspath.(second_pairs)

    for file in vcat(reference_fastas_abs, sample_reads_abs, first_pairs_abs, second_pairs_abs)
        if !isfile(file)
            error("File not found: $(file)")
        end
    end

    Mycelia.add_bioconda_env("sylph")

    workdir_is_temp = isnothing(outdir)
    workdir = workdir_is_temp ? mktempdir() : mkpath(outdir)
    workdir = abspath(workdir)
    db_prefix = joinpath(workdir, output_prefix)
    syldb_path = db_prefix * ".syldb"
    sample_dir = workdir

    sketch_args = ["sketch", "-t", string(threads), "-c", string(subsampling), "-k", string(k),
                   "--min-spacing", string(min_spacing), "-o", db_prefix, "-d", sample_dir]
    for fasta in reference_fastas_abs
        push!(sketch_args, "-g")
        push!(sketch_args, fasta)
    end
    if !isempty(first_pairs_abs)
        push!(sketch_args, "-1")
        append!(sketch_args, first_pairs_abs)
    end
    if !isempty(second_pairs_abs)
        push!(sketch_args, "-2")
        append!(sketch_args, second_pairs_abs)
    end
    append!(sketch_args, sample_reads_abs)

    sketch_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n sylph sylph $sketch_args`
    if workdir_is_temp
        sketch_cmd = Cmd(sketch_cmd; dir=workdir)
    end
    if quiet
        run(pipeline(sketch_cmd, stdout=devnull, stderr=devnull))
    else
        run(sketch_cmd)
    end

    sample_sketches = filter(f -> endswith(f, ".sylsp"), readdir(sample_dir; join=true))
    if isempty(sample_sketches)
        error("No sample sketches (*.sylsp) were produced by sylph")
    end
    if !isfile(syldb_path)
        error("Sylph database not found at $(syldb_path)")
    end

    profile_out = isnothing(output_tsv) ? joinpath(workdir, "sylph_profile.tsv") : abspath(output_tsv)
    profile_args = ["profile", "-t", string(threads), "-o", profile_out, "-m", string(min_ani), "-M", string(min_kmers)]
    if estimate_unknown
        push!(profile_args, "-u")
    end
    append!(profile_args, additional_args)
    append!(profile_args, vcat(syldb_path, sample_sketches))

    profile_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n sylph sylph $profile_args`
    if workdir_is_temp
        profile_cmd = Cmd(profile_cmd; dir=workdir)
    end
    if quiet
        run(pipeline(profile_cmd, stdout=devnull, stderr=devnull))
    else
        run(profile_cmd)
    end

    table = CSV.read(profile_out, DataFrames.DataFrame; delim='\t', normalizenames=true)
    return SylphProfileResult(syldb_path, sample_sketches, profile_out, table)
end

"""
    select_sketch_supported_references(reference_scores::AbstractDict{<:AbstractString,<:Real};
                                       min_score::Union{Nothing,Real}=nothing,
                                       max_score::Union{Nothing,Real}=nothing,
                                       max_refs::Union{Nothing,Int}=nothing,
                                       prefer::Symbol=:higher)

Filter and rank sketch scores to identify the pangenome contexts supported by a sample.

Use `min_score` with containment/coverage scores (sourmash/sylph). Use `max_score` and
`prefer=:lower` for distance-style scores (mash). Returns a vector of `Pair{String,Float64}`
sorted by score.
"""
function select_sketch_supported_references(reference_scores::AbstractDict{<:AbstractString,<:Real};
        min_score::Union{Nothing,Real}=nothing,
        max_score::Union{Nothing,Real}=nothing,
        max_refs::Union{Nothing,Int}=nothing,
        prefer::Symbol=:higher)

    prefer in (:higher, :lower) || error("prefer must be :higher or :lower")
    if !isnothing(min_score) && !isnothing(max_score) && min_score > max_score
        error("min_score must be <= max_score")
    end
    if !isnothing(max_refs) && max_refs < 0
        error("max_refs must be non-negative")
    end

    selected = Pair{String,Float64}[]
    for (reference, score) in reference_scores
        score_value = Float64(score)
        if !isnothing(min_score) && score_value < min_score
            continue
        end
        if !isnothing(max_score) && score_value > max_score
            continue
        end
        push!(selected, String(reference) => score_value)
    end

    sort!(selected, by=last, rev=(prefer == :higher))
    if !isnothing(max_refs)
        selected = selected[1:min(max_refs, length(selected))]
    end

    return selected
end

"""
    skani_triangle(fasta_files::Vector{String};
                   small_genomes::Bool=false,
                   sparse::Union{Bool,Nothing}=nothing,
                   threads::Int=3,
                   min_af::Union{Nothing,Float64}=nothing,
                   output_file::Union{Nothing,String}=nothing,
                   additional_args::Vector{String}=String[])

Perform pairwise all-vs-all ANI/AF comparison using skani triangle.

# Arguments
- `fasta_files::Vector{String}`: Vector of paths to FASTA files to compare
- `small_genomes::Bool=false`: Use `--small-genomes` flag for viral/plasmid genomes
- `sparse::Union{Bool,Nothing}=nothing`: Output format. If `nothing` (default), automatically uses sparse
  format (`-E`) when >500 genomes, dense otherwise. Set explicitly to `true` or `false` to override.
- `threads::Int=3`: Number of threads to use
- `min_af::Union{Nothing,Float64}=nothing`: Minimum aligned fraction threshold
- `output_file::Union{Nothing,String}=nothing`: Output file path. If nothing, returns stdout as string
- `additional_args::Vector{String}=String[]`: Additional command-line arguments to pass to skani

# Returns
- If `output_file` is specified: writes to file and returns the file path
- If `output_file` is nothing: returns the output as a String

# Notes
For >500 genomes, sparse output is recommended to avoid large matrix files. The sparse format
produces TSV output similar to `skani dist`, while dense format produces a triangular matrix.
"""
function skani_triangle(fasta_files::Vector{String};
                       small_genomes::Bool=false,
                       sparse::Union{Bool,Nothing}=nothing,
                       threads::Int=3,
                       min_af::Union{Nothing,Float64}=nothing,
                       output_file::Union{Nothing,String}=nothing,
                       additional_args::Vector{String}=String[])

    for file in fasta_files
        if !isfile(file)
            error("File not found: $file")
        end
    end

    # Auto-detect sparse mode based on genome count (skani recommends sparse for >500 genomes)
    n_genomes = length(fasta_files)
    use_sparse = if isnothing(sparse)
        if n_genomes > 500
            @info "Auto-enabling sparse mode (-E) for $n_genomes genomes (>500 threshold)"
            true
        else
            false
        end
    else
        sparse
    end

    Mycelia.add_bioconda_env("skani")

    temp_dir = mktempdir()
    list_file = joinpath(temp_dir, "skani_input_list.txt")

    try
        open(list_file, "w") do io
            for file in fasta_files
                println(io, abspath(file))
            end
        end

        cmd_args = ["triangle", "-l", list_file, "-t", string(threads)]

        if small_genomes
            push!(cmd_args, "--small-genomes")
        end

        if use_sparse
            push!(cmd_args, "-E")
        end
        
        if !isnothing(min_af)
            push!(cmd_args, "--min-af", string(min_af))
        end
        
        if !isnothing(output_file)
            push!(cmd_args, "-o", output_file)
        end
        
        append!(cmd_args, additional_args)
        
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n skani skani $cmd_args`
        
        if isnothing(output_file)
            result = read(cmd, String)
            return result
        else
            run(cmd)
            return output_file
        end
        
    finally
        rm(temp_dir; recursive=true, force=true)
    end
end

"""
    skani_dist(fasta_files::Vector{String};
               threads::Int=3,
               small_genomes::Bool=false,
               min_af::Union{Nothing,Float64}=nothing,
               output_file::Union{Nothing,String}=nothing,
               additional_args::Vector{String}=String[])

Run `skani dist` on a list of FASTA files and return the parsed TSV as a DataFrame.
"""
function skani_dist(fasta_files::Vector{String};
        threads::Int=3,
        small_genomes::Bool=false,
        min_af::Union{Nothing,Float64}=nothing,
        output_file::Union{Nothing,String}=nothing,
        additional_args::Vector{String}=String[])

    for file in fasta_files
        if !isfile(file)
            error("File not found: $file")
        end
    end

    Mycelia.add_bioconda_env("skani")

    temp_dir = mktempdir()
    list_file = joinpath(temp_dir, "skani_input_list.txt")
    out_path = isnothing(output_file) ? joinpath(temp_dir, "skani_dist.tsv") : output_file

    try
        open(list_file, "w") do io
            for file in fasta_files
                println(io, abspath(file))
            end
        end

        cmd_args = ["dist", "--ql", list_file, "--rl", list_file, "-t", string(threads), "-o", out_path]
        if small_genomes
            push!(cmd_args, "--small-genomes")
        end
        if !isnothing(min_af)
            push!(cmd_args, "--min-af", string(min_af))
        end
        append!(cmd_args, additional_args)

        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n skani skani $cmd_args`
        run(cmd)

        df = CSV.read(out_path, DataFrames.DataFrame; delim='\t', normalizenames=true)
        return df
    finally
        rm(temp_dir; recursive=true, force=true)
    end
end

"""
    ngrams(seq::AbstractString, k::Integer)

Return the overlapping k-length substrings of `seq` (character-wise).
Throws `BoundsError` when `seq` is empty or `k < 1`.
"""
function ngrams(seq::AbstractString, k::Integer)
    sstr = String(seq)
    if k < 1
        throw(BoundsError(sstr, k))
    end
    if isempty(sstr)
        throw(BoundsError(sstr, k))
    end

    idxs = collect(eachindex(sstr))
    n = length(idxs)
    if k > n
        return SubString{String}[]
    end

    grams = Vector{SubString{String}}(undef, n - k + 1)
    for i in 1:(n - k + 1)
        start_idx = idxs[i]
        stop_idx = idxs[i + k - 1]
        grams[i] = SubString(sstr, start_idx, stop_idx)
    end
    return grams
end


"""
    shannon_entropy(seq; k::Int=1, base::Real=2, alphabet=nothing)

Compute Shannon entropy of symbols (k=1) or k-mers (k>1) from `seq`
(length-agnostic because it uses relative frequencies). `base` sets the log base.
`alphabet` (optional) lets you force/declare the alphabet size; otherwise inferred.
Returns a Float64.
"""
function shannon_entropy(seq; k::Int=1, base::Real=2, alphabet=nothing)
    # Use existing biosequence infrastructure for biological sequences
    if seq isa BioSequences.BioSequence
        if k > length(seq)
            return 0.0
        end
        
        # Determine k-mer type and counting function based on sequence alphabet
        # TODO: use the constants in this repo instead of hardcoding
        if seq isa BioSequences.LongDNA
            kmer_type = Kmers.DNAKmer{k}
            counts = count_kmers(kmer_type, seq)  # Use non-canonical for entropy
        elseif seq isa BioSequences.LongRNA
            kmer_type = Kmers.RNAKmer{k}
            counts = count_kmers(kmer_type, seq)
        elseif seq isa BioSequences.LongAA
            kmer_type = Kmers.AAKmer{k}
            counts = count_kmers(kmer_type, seq)
        else
            error("Unsupported sequence type: $(typeof(seq))")
        end
    else
        # Use existing ngram infrastructure for generic strings
        sstr = string(seq)
        if k > lastindex(sstr)
            return 0.0
        end
        
        observed_ngrams = ngrams(sstr, k)
        counts = StatsBase.countmap(observed_ngrams)
    end
    
    total = sum(values(counts))
    probs = (v/total for v in values(counts))
    H = -sum(p -> p == 0 ? 0.0 : p * (log(p)/log(base)), probs)
    return H
end

"""
    renyi_entropy(seq; k::Int=1, α::Real=2, base::Real=2)

Rényi entropy of order α (α ≠ 1). For α→1, this tends to Shannon entropy.
"""
function renyi_entropy(seq; k::Int=1, α::Real=2, base::Real=2)
    @assert α != 1 "Use shannon_entropy for α = 1."
    
    # Use existing biosequence infrastructure for biological sequences
    if seq isa BioSequences.BioSequence
        if k > length(seq)
            return 0.0
        end
        
        # Determine k-mer type and counting function based on sequence alphabet
        # TODO: use the constants in this repo instead of hardcoding
        if seq isa BioSequences.LongDNA
            kmer_type = Kmers.DNAKmer{k}
            counts = count_kmers(kmer_type, seq)  # Use non-canonical for entropy
        elseif seq isa BioSequences.LongRNA
            kmer_type = Kmers.RNAKmer{k}
            counts = count_kmers(kmer_type, seq)
        elseif seq isa BioSequences.LongAA
            kmer_type = Kmers.AAKmer{k}
            counts = count_kmers(kmer_type, seq)
        else
            error("Unsupported sequence type: $(typeof(seq))")
        end
    else
        # Use existing ngram infrastructure for generic strings
        sstr = string(seq)
        if k > lastindex(sstr)
            return 0.0
        end
        
        observed_ngrams = ngrams(sstr, k)
        counts = StatsBase.countmap(observed_ngrams)
    end
    
    total = sum(values(counts))
    probsα = sum((v/total)^α for v in values(counts))
    return (1/(1-α)) * (log(probsα)/log(base))
end

"""
    kmer_richness(seq, k; alphabet=nothing, normalize=true)

Return the number of unique k-mers observed.
If `normalize=true` (default), return the linguistic-complexity-style ratio:
unique / min(L-k+1, A^k), where A is alphabet size.
"""
function kmer_richness(seq, k; alphabet=nothing, normalize::Bool=true)
    # Use existing biosequence infrastructure for biological sequences
    if seq isa BioSequences.BioSequence
        L = length(seq)
        if k > L
            return normalize ? 0.0 : 0
        end
        
        # Determine k-mer type and counting function based on sequence alphabet
        # TODO: use the constants in this repo instead of hardcoding
        if seq isa BioSequences.LongDNA
            kmer_type = Kmers.DNAKmer{k}
            kmer_counts = count_canonical_kmers(kmer_type, seq)
            A = isnothing(alphabet) ? 4 : alphabet
        elseif seq isa BioSequences.LongRNA
            kmer_type = Kmers.RNAKmer{k}
            kmer_counts = count_kmers(kmer_type, seq)
            A = isnothing(alphabet) ? 4 : alphabet
        elseif seq isa BioSequences.LongAA
            kmer_type = Kmers.AAKmer{k}
            kmer_counts = count_kmers(kmer_type, seq)
            A = isnothing(alphabet) ? 20 : alphabet
        else
            error("Unsupported sequence type: $(typeof(seq))")
        end
        
        uniq = length(kmer_counts)
    else
        # Use existing ngram infrastructure for generic strings
        sstr = string(seq)
        L = lastindex(sstr)
        if k > L
            return normalize ? 0.0 : 0
        end
        
        observed_ngrams = ngrams(sstr, k)
        uniq = length(unique(observed_ngrams))
        A = isnothing(alphabet) ? length(unique(sstr)) : alphabet
    end
    
    if !normalize
        return uniq
    end
    
    max_possible = min(L - k + 1, A^k)
    return max_possible == 0 ? 0.0 : uniq / max_possible
end

"""
    linguistic_complexity(seq; kmax=nothing, alphabet=nothing, reducer=mean)

Compute richness ratios for k = 1:kmax and reduce them (mean by default).
If `kmax` is omitted, it defaults to the full range 1:floor(Int, log_A(L)) or simply 1:L.
Returns (profile, summary), where profile is a Vector{Float64} of length kmax
and summary is `reducer(profile)`.
"""
function linguistic_complexity(seq; kmax=nothing, alphabet=nothing, reducer=Statistics.mean)
    # Determine sequence type and properties
    if seq isa BioSequences.BioSequence
        L = length(seq)
        # Infer alphabet size from sequence type if not provided
        # TODO: use the constants in this repo instead of hardcoding
        if isnothing(alphabet)
            if seq isa BioSequences.LongDNA
                A = 4
            elseif seq isa BioSequences.LongRNA
                A = 4
            elseif seq isa BioSequences.LongAA
                A = 20
            else
                A = 4  # default fallback
            end
        else
            A = alphabet
        end
    else
        # Generic string case
        sstr = string(seq)
        L = lastindex(sstr)
        A = isnothing(alphabet) ? length(unique(sstr)) : alphabet
    end
    
    kmax = isnothing(kmax) ? L : kmax
    prof = [kmer_richness(seq, k; alphabet=A, normalize=true) for k in 1:kmax]
    return prof, reducer(prof)
end

"""
    merge_and_map_single_end_samples(; 
        fasta_reference::AbstractString, 
        fastq_list::Vector{<:AbstractString}, 
        minimap_index::AbstractString, 
        mapping_type::AbstractString,
        outbase::AbstractString = "results",
        outformats::Vector{<:AbstractString} = [".tsv.gz", ".jld2"]
    ) -> DataFrames.DataFrame

Merge and map single-end sequencing samples, then output results in one or more formats.

# Arguments
- `fasta_reference`: Path to the reference FASTA file.
- `fastq_list`: Vector of paths to input FASTQ files to be merged.
- `minimap_index`: Path to the minimap2 index file (.mmi).
- `mapping_type`: Mapping type string for minimap2 (e.g., "map-ont").
- `outbase`: Base name (optionally including path) for output files (default: `Mycelia.normalized_current_date() * ".joint-minimap-mapping-results"`).
- `outformats`: Vector of output file formats to write results to. Supported: `".tsv.gz"`, `".jld2"`.

# Description
This function merges provided FASTQ files and assigns unique UUIDs to reads, maps the merged FASTQ against the provided reference using minimap2, reads mapping and UUID tables, joins them into a single DataFrame, writes this table to all requested output formats with filenames constructed from the `outbase` and the appropriate extension, and returns the resulting joined DataFrame.

# Output Files
- `.tsv.gz`: Tab-separated, gzip-compressed table of results.
- `".jld2"`: JLD2 file containing results.

# Returns
- The joined results as a `DataFrames.DataFrame`.
"""
function merge_and_map_single_end_samples(; 
    fasta_reference::AbstractString, 
    fastq_list::Vector{<:AbstractString}, 
    minimap_index::AbstractString, 
    mapping_type::AbstractString,
    outbase::AbstractString = Mycelia.normalized_current_date() * ".joint-minimap-mapping-results",
    outformats::Vector{<:AbstractString} = [".tsv.gz", ".jld2"]
)
    # Join FASTQ files
    fastq_out = outbase * ".fq.gz"
    tsv_out = outbase * ".uuid-map.tsv.gz"
    fastq_join_result = Mycelia.join_fastqs_with_uuid(fastq_list, fastq_out=fastq_out, tsv_out=tsv_out)
    
    # Run minimap if needed
    minimap_result = Mycelia.minimap_map_with_index(
        fasta = fasta_reference,
        mapping_type = mapping_type,
        fastq = fastq_out,
        index_file = minimap_index
    )
    if !isfile(minimap_result.outfile)
        @time run(minimap_result.cmd)
    end
    results_table_outfiles = [outbase * fmt for fmt in outformats]
    # Determine file paths for .tsv.gz and .jld2
    tsv_file = ".tsv.gz" in outformats ? outbase * ".tsv.gz" : nothing
    jld2_file = ".jld2" in outformats ? outbase * ".jld2" : nothing
    results_table = nothing
    if !all(isfile, results_table_outfiles) || any(x -> filesize(x) == 0, results_table_outfiles)
        # Read tables
        read_id_mapping_table = CSV.read(
            CodecZlib.GzipDecompressorStream(open(tsv_out)), DataFrames.DataFrame, delim='\t')
        mapping_results_table = Mycelia.xam_to_dataframe(minimap_result.outfile)
        results_table = DataFrames.innerjoin(read_id_mapping_table, mapping_results_table, on="new_uuid" => "template")
        results_table = Mycelia.dataframe_replace_nothing_with_missing(results_table)
        # Write outputs
        for fmt in outformats
            outfile = outbase * fmt
            if fmt == ".tsv.gz"
                @assert outfile == tsv_file
                io = CodecZlib.GzipCompressorStream(open(outfile, "w"))
                CSV.write(io, results_table; delim='\t')
                close(io)
                @assert isfile(tsv_file)
                @show tsv_file
            elseif fmt == ".jld2"
                @assert outfile == jld2_file
                JLD2_write_table(df=results_table, filename=outfile)
                @assert isfile(jld2_file)
            else
                @warn "Unknown output format: $fmt"
            end
        end
    else
        if !isnothing(jld2_file) && isfile(jld2_file)
            results_table = JLD2_read_table(jld2_file)
        elseif !isnothing(tsv_file) && isfile(tsv_file)
            io = CodecZlib.GzipDecompressorStream(open(tsv_file))
            results_table = CSV.read(io, DataFrames.DataFrame; delim='\t')
            close(io)
        end
    end

    return (
        results_table = results_table,
        joint_fastq_file = fastq_out,
        fastq_id_mapping_table = tsv_out,
        bam_file = minimap_result.outfile,
        tsv_file = tsv_file,
        jld2_file = jld2_file
    )
end

"""
    mash_distance_from_jaccard(jaccard_index::Float64, kmer_size::Int)

Calculates the Mash distance (an estimate of Average Nucleotide Identity)
from a given Jaccard Index and k-mer size.

# Arguments
- `jaccard_index::Float64`: The Jaccard similarity between the two k-mer sets. Must be between 0.0 and 1.0.
- `kmer_size::Int`: The length of k-mers used to calculate the Jaccard index.

# Returns
- `Float64`: The estimated Mash distance `D`. The estimated ANI would be `1.0 - D`.
"""
function mash_distance_from_jaccard(jaccard_index::Float64, kmer_size::Int)
    # --- Input Validation ---
    if !(0.0 <= jaccard_index <= 1.0)
    error("Jaccard index must be between 0.0 and 1.0")
    end
    if kmer_size <= 0
    error("k-mer size must be a positive integer")
    end

    # --- Edge Case Handling ---
    # If Jaccard is 0, the genomes share no k-mers. The distance is effectively infinite,
    # conventionally represented as 1.0 (100% divergent). The formula would fail due to log(0).
    if jaccard_index == 0.0
        return 1.0
    end

    # If Jaccard is 1, the genomes are identical. Distance is 0.
    if jaccard_index == 1.0
        return 0.0
    end

    # --- Core Mash Formula ---
    # D = - (1/k) * ln(2J / (1+J))
    # In Julia, log() is the natural logarithm (ln).
    mash_dist = - (1 / kmer_size) * log(2 * jaccard_index / (1 + jaccard_index))

    return mash_dist
end

"""
    run_mash_comparison(fasta1::String, fasta2::String; k::Int=21, s::Int=10000, mash_path::String="mash")

Runs a genome-by-genome comparison using the `mash` command-line tool.

This function first creates sketch files for each FASTA input and then
calculates the distance between them, capturing and parsing the result.

# Arguments
- `fasta1::String`: Path to the first FASTA file.
- `fasta2::String`: Path to the second FASTA file.

# Keyword Arguments
- `k::Int=21`: The k-mer size to use for sketching. Default is 21.
- `s::Int=10000`: The sketch size (number of hashes to keep). Default is 10000.
- `mash_path::String="mash"`: The path to the mash executable if not in the system PATH.

# Returns
- `NamedTuple`: A named tuple containing the parsed results, e.g.,
  `(reference=..., query=..., distance=..., p_value=..., shared_hashes=...)`
- `nothing`: Returns `nothing` if the `mash` command fails.
"""
function run_mash_comparison(fasta1::String, fasta2::String; k::Int=21, s::Int=10000, mash_path::String="mash")
    # --- Step 1: Check if input files exist ---
    if !isfile(fasta1) || !isfile(fasta2)
        error("One or both FASTA files not found.")
    end

    Mycelia.add_bioconda_env("mash")

    # --- Step 2: Create sketch files for each genome ---
    sketch1 = fasta1 * ".msh"
    sketch2 = fasta2 * ".msh"

    println("Sketching $fasta1 (k=$k, s=$s)...")
    sketch_cmd1 = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash sketch -k $k -s $s -o $sketch1 $fasta1`, stdout=devnull, stderr=devnull)

    println("Sketching $fasta2 (k=$k, s=$s)...")
    sketch_cmd2 = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash sketch -k $k -s $s -o $sketch2 $fasta2`, stdout=devnull, stderr=devnull)

    try
        run(sketch_cmd1)
        run(sketch_cmd2)
    catch e
        println("Error: Failed to run 'mash sketch'. Is mash installed and in your PATH?")
        println(e)
        return nothing
    end

    # --- Step 3: Run 'mash dist' on the two sketches and capture output ---
    println("Calculating distance between sketches...")
    dist_cmd = `$mash_path dist $sketch1 $sketch2`

    output = ""
    try
        # read() captures the standard output of the command
        output = read(dist_cmd, String)
    catch e
        println("Error: Failed to run 'mash dist'.")
        println(e)
        return nothing
    finally
        # --- Step 4: Clean up the sketch files ---
        rm(sketch1, force=true)
        rm(sketch2, force=true)
    end

    # --- Step 5: Parse the tab-separated output from mash ---
    if isempty(output)
        println("Warning: Mash command produced no output.")
        return nothing
    end

    # Example output: "genomeA.fasta\tgenomeB.fasta\t0.080539\t0.0\t491/1000"
    parts = split(strip(output), '\t')

    if length(parts) != 5
        println("Error: Unexpected output format from Mash: ", output)
        return nothing
    end

    parsed_result = (
        reference = parts[1],
        query = parts[2],
        distance = parse(Float64, parts[3]),
        p_value = parse(Float64, parts[4]),
        shared_hashes = parts[5]
    )

    return parsed_result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compare two FASTA sequences and calculate alignment statistics.

# Arguments
- `reference_fasta::String`: Path to the reference FASTA file
- `query_fasta::String`: Path to the query FASTA file

# Returns
DataFrame with the following columns:
- `alignment_percent_identity`: Percentage of matching bases in alignment
- `total_equivalent_bases`: Number of equivalent bases between sequences
- `total_alignment_length`: Length of the alignment
- `query_length`: Length of query sequence
- `total_variants`: Total number of variants (SNPs + indels)
- `total_snps`: Number of single nucleotide polymorphisms
- `total_indels`: Number of insertions and deletions
- `alignment_coverage_query`: Percentage of query sequence covered
- `alignment_coverage_reference`: Percentage of reference sequence covered
- `size_equivalence_to_reference`: Size ratio of query to reference (%)

# Notes
- Uses minimap2 with progressively relaxed settings (asm5→asm10→asm20)
- Returns empty string values for alignment statistics if no alignment is found
- Requires minimap2 to be installed and accessible in PATH
"""
# uses minimap
function pairwise_minimap_fasta_comparison(;reference_fasta, query_fasta)
    header = [
        "Query",
        "Query length",
        "Query start",
        "Query end",
        "Query strand",
        "Target",
        "Target length",
        "Target start",
        "Target end",
        "Matches",
        "Alignment length",
        "Mapping quality",
        "Cigar",
        "CS tag"]

    Mycelia.add_bioconda_env("minimap2")
#     asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    results5 = read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -x asm5 --cs -cL $reference_fasta $query_fasta`)
    if !isempty(results5)
        results = results5
    else
        @warn "no hit with asm5, trying asm10"
        results10 = read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -x asm10 --cs -cL $reference_fasta $query_fasta`)
        if !isempty(results10)
            results = results10
        else
            @warn "no hits with asm5 or asm10, trying asm20"
            results20 = read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -x asm20 --cs -cL $reference_fasta $query_fasta`)
            if !isempty(results20)
                results = results20
            end
        end
    end
    if !isempty(results)
        data =  DelimitedFiles.readdlm(IOBuffer(results), '\t')
        data_columns_of_interest = [collect(1:length(header)-2)..., collect(size(data, 2)-1:size(data, 2))...]
        minimap_results = DataFrames.DataFrame(data[:, data_columns_of_interest], header)

        equivalent_matches = reduce(vcat, map(x -> collect(eachmatch(r":([0-9]+)", replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_equivalent_bases = sum(map(match -> parse(Int, first(match.captures)), equivalent_matches))

        insertion_matches = reduce(vcat, map(x -> collect(eachmatch(r"\+([a-z]+)"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_inserted_bases = sum(map(match -> length(first(match.captures)), insertion_matches))
        deletion_matches = reduce(vcat, map(x -> collect(eachmatch(r"\-([a-z]+)"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_deleted_bases = sum(map(match -> length(first(match.captures)), deletion_matches))
        substitution_matches = reduce(vcat, map(x -> collect(eachmatch(r"\*([a-z]{2})"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_substituted_bases = length(substitution_matches)
        total_variants = length(insertion_matches) + length(deletion_matches) + length(substitution_matches)
        total_variable_bases = total_inserted_bases + total_deleted_bases + total_substituted_bases

        total_alignment_length = sum(minimap_results[!, "Alignment length"])
        total_matches = sum(minimap_results[!, "Matches"])
        
        alignment_percent_identity = round(total_matches / total_alignment_length * 100, digits=2)
        size_equivalence_to_reference = round(minimap_results[1, "Query length"]/minimap_results[1, "Target length"] * 100, digits=2)
        alignment_coverage_query = round(total_alignment_length / minimap_results[1, "Query length"] * 100, digits=2)
        alignment_coverage_reference = round(total_alignment_length / minimap_results[1, "Target length"] * 100, digits=2)

        results = DataFrames.DataFrame(
            alignment_percent_identity = alignment_percent_identity,
            total_equivalent_bases = total_equivalent_bases,
            total_alignment_length = total_alignment_length,
            query_length = minimap_results[1, "Query length"],
            total_variants = total_variants,
            total_snps = total_substituted_bases,
            total_indels = length(insertion_matches) + length(deletion_matches),
            alignment_coverage_query = alignment_coverage_query,
            alignment_coverage_reference = alignment_coverage_reference,
            size_equivalence_to_reference = size_equivalence_to_reference,
        )
    else
        query_length = length(FASTX.sequence(first(FASTX.FASTA.Reader(open(query_fasta)))))
        target_length = length(FASTX.sequence(first(FASTX.FASTA.Reader(open(reference_fasta)))))
        size_equivalence_to_reference = round(query_length/target_length * 100, digits=2)

        # unable to find any matches
        results = DataFrames.DataFrame(
            alignment_percent_identity = "",
            total_equivalent_bases = "",
            total_alignment_length = "",
            query_length = query_length,
            total_variants = "",
            total_snps = "",
            total_indels = "",
            alignment_coverage_query = 0,
            alignment_coverage_reference = 0,
            size_equivalence_to_reference = size_equivalence_to_reference
        )
    end
    return results
end

# always interpret as strings to ensure changes in underlying biosequence representation don't change results
# results in 64 character string
# a = "TTANC"
# b = "ttANc"
# c = "ttanc"
# dna_a = BioSequences.LongDNA{4}(a)
# dna_b = BioSequences.LongDNA{4}(b)
# dna_c = BioSequences.LongDNA{4}(c)
# seq2sha256(a) == seq2sha256(dna_a)
# seq2sha256(b) == seq2sha256(dna_b)
# seq2sha256(c) == seq2sha256(dna_c)
# seq2sha256(BioSequences.LongDNA{2}("AAA")) == seq2sha256(BioSequences.LongDNA{4}("AAA"))

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the SHA-256 hash of a sequence string.

# Arguments
- `seq::AbstractString`: Input sequence to be hashed

# Returns
- `String`: Hexadecimal representation of the SHA-256 hash

# Details
The input sequence is converted to uppercase before hashing.
"""
function seq2sha256(seq::AbstractString)
    return SHA.bytes2hex(SHA.sha256(uppercase(seq)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a biological sequence to its SHA256 hash value.

Calculates a cryptographic hash of the sequence by first converting it to a string representation.
This method dispatches to the string version of `seq2sha256`.

# Arguments
- `seq::BioSequences.BioSequence`: The biological sequence to hash

# Returns
- `String`: A 64-character hexadecimal string representing the SHA256 hash
"""
function seq2sha256(seq::BioSequences.BioSequence)
    return seq2sha256(string(seq))
end

const DNA_COMPLEMENT_TABLE = let table = Vector{UInt8}(undef, 256)
    for i in 0:255
        table[i + 1] = UInt8(i)
    end
    complements = Dict(
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'U' => 'A',
        'R' => 'Y',
        'Y' => 'R',
        'S' => 'S',
        'W' => 'W',
        'K' => 'M',
        'M' => 'K',
        'B' => 'V',
        'D' => 'H',
        'H' => 'D',
        'V' => 'B',
        'N' => 'N'
    )
    for (base, comp) in complements
        table[Int(UInt8(base)) + 1] = UInt8(comp)
    end
    table
end

function genome_pair_id(reference::AbstractString, query::AbstractString)
    ref_abs = abspath(reference)
    query_abs = abspath(query)
    ref_stat = stat(ref_abs)
    query_stat = stat(query_abs)
    payload = join([
        ref_abs,
        string(ref_stat.size),
        string(ref_stat.mtime),
        query_abs,
        string(query_stat.size),
        string(query_stat.mtime)
    ], "|")
    return SHA.bytes2hex(SHA.sha256(payload))
end

function count_fasta_records(fasta_path::AbstractString)
    @assert isfile(fasta_path) "FASTA file does not exist: $(fasta_path)"
    count = 0
    open(fasta_path) do io
        reader = FASTX.FASTA.Reader(io)
        for _ in reader
            count += 1
        end
    end
    return count
end

function reverse_complement_ascii(seq::AbstractString)
    bytes = Vector{UInt8}(codeunits(seq))
    n = length(bytes)
    out = Vector{UInt8}(undef, n)
    for i in 1:n
        out[i] = DNA_COMPLEMENT_TABLE[bytes[n - i + 1] + 1]
    end
    return String(out)
end

function booth_min_rotation_index(seq_bytes::Vector{UInt8})
    n = length(seq_bytes)
    if n == 0
        return 1
    end
    doubled = Vector{UInt8}(undef, 2 * n)
    copyto!(doubled, 1, seq_bytes, 1, n)
    copyto!(doubled, n + 1, seq_bytes, 1, n)

    i = 1
    j = 2
    k = 0
    while i <= n && j <= n && k < n
        a = doubled[i + k]
        b = doubled[j + k]
        if a == b
            k += 1
        elseif a < b
            j = j + k + 1
            if j == i
                j += 1
            end
            k = 0
        else
            i = i + k + 1
            if i == j
                i += 1
            end
            k = 0
        end
    end
    return min(i, j)
end

function rotate_bytes(seq_bytes::Vector{UInt8}, start_index::Int)
    n = length(seq_bytes)
    if n == 0
        return ""
    end
    start_index = ((start_index - 1) % n) + 1
    if start_index == 1
        return String(seq_bytes)
    end
    rotated = Vector{UInt8}(undef, n)
    tail_len = n - start_index + 1
    copyto!(rotated, 1, seq_bytes, start_index, tail_len)
    copyto!(rotated, tail_len + 1, seq_bytes, 1, start_index - 1)
    return String(rotated)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return a canonical linear representation of a circular sequence.

Selects the lexicographically minimal rotation across the forward strand and its
reverse complement, returning the canonical sequence and the chosen orientation.
"""
function canonical_circular_sequence(seq::AbstractString)
    seq_upper = uppercase(seq)
    forward_bytes = Vector{UInt8}(codeunits(seq_upper))
    forward_start = booth_min_rotation_index(forward_bytes)
    forward_rot = rotate_bytes(forward_bytes, forward_start)

    reverse_seq = reverse_complement_ascii(seq_upper)
    reverse_bytes = Vector{UInt8}(codeunits(reverse_seq))
    reverse_start = booth_min_rotation_index(reverse_bytes)
    reverse_rot = rotate_bytes(reverse_bytes, reverse_start)

    if reverse_rot < forward_rot
        return reverse_rot, :reverse, reverse_start
    end
    return forward_rot, :forward, forward_start
end

function header_says_circular(header::AbstractString)
    return occursin("circular", lowercase(header))
end

function has_terminal_overlap(seq::AbstractString; overlap_bp::Int=200)
    if overlap_bp <= 0
        return false
    end
    seq_len = length(seq)
    if seq_len < 2 * overlap_bp
        return false
    end
    return seq[1:overlap_bp] == seq[end - overlap_bp + 1:end]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalize a FASTA file by uppercasing sequences and optionally sorting records.
"""
function normalize_fasta(fasta_path::AbstractString; out_path::AbstractString, sort_records::Bool=false, force::Bool=false)
    @assert isfile(fasta_path) "FASTA file does not exist: $(fasta_path)"
    @assert !isempty(out_path) "out_path must be a non-empty string"
    mkpath(dirname(out_path))

    if !force && isfile(out_path) && filesize(out_path) > 0
        return out_path
    end

    records = FASTX.FASTA.Record[]
    open(fasta_path) do io
        reader = FASTX.FASTA.Reader(io)
        for record in reader
            push!(records, record)
        end
    end

    if sort_records
        sort!(records, by=FASTX.identifier)
    end

    open(FASTX.FASTA.Writer, out_path) do writer
        for record in records
            description = FASTX.description(record)
            seq = String(FASTX.sequence(record))
            seq = uppercase(seq)
            header = isnothing(description) || isempty(description) ?
                     String(FASTX.identifier(record)) :
                     String(description)
            fasta_record = FASTX.FASTA.Record(header, seq)
            write(writer, fasta_record)
        end
    end

    return out_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Canonicalize a single-contig circular FASTA by rotating to a canonical start and strand.
"""
function canonicalize_circular_fasta(fasta_path::AbstractString; out_path::AbstractString, force::Bool=false)
    @assert isfile(fasta_path) "FASTA file does not exist: $(fasta_path)"
    @assert !isempty(out_path) "out_path must be a non-empty string"
    mkpath(dirname(out_path))

    if !force && isfile(out_path) && filesize(out_path) > 0
        return out_path
    end

    record = nothing
    record_count = 0
    open(fasta_path) do io
        reader = FASTX.FASTA.Reader(io)
        for rec in reader
            record_count += 1
            if record_count == 1
                record = rec
            else
                break
            end
        end
    end

    @assert record_count == 1 "Circular canonicalization expects a single-contig FASTA: $(fasta_path)"
    description = FASTX.description(record)
    seq = String(FASTX.sequence(record))
    canonical_seq, _, _ = canonical_circular_sequence(seq)
    header = isnothing(description) || isempty(description) ?
             String(FASTX.identifier(record)) :
             String(description)
    fasta_record = FASTX.FASTA.Record(header, canonical_seq)

    open(FASTX.FASTA.Writer, out_path) do writer
        write(writer, fasta_record)
    end

    return out_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Prepare a genome FASTA for alignment-based comparison with optional circular normalization.
"""
function prepare_genome_for_comparison(fasta_path::AbstractString;
        outdir::AbstractString,
        circular::Union{Bool,Symbol}=:auto,
        sort_records::Bool=false,
        circular_heuristic::Bool=false,
        circular_overlap_bp::Int=200,
        force::Bool=false)
    @assert circular in (true, false, :auto) "circular must be true, false, or :auto"
    if circular_heuristic
        @assert circular_overlap_bp > 0 "circular_overlap_bp must be positive"
    end
    mkpath(outdir)
    normalized = joinpath(outdir, "normalized.fna")
    normalize_fasta(fasta_path; out_path=normalized, sort_records=sort_records, force=force)

    if circular == false
        return normalized
    end

    record_count = 0
    header = ""
    seq = ""
    open(normalized) do io
        reader = FASTX.FASTA.Reader(io)
        for record in reader
            record_count += 1
            if record_count == 1
                description = FASTX.description(record)
                header = isnothing(description) || isempty(description) ?
                         String(FASTX.identifier(record)) :
                         string(FASTX.identifier(record), " ", description)
                seq = String(FASTX.sequence(record))
            else
                break
            end
        end
    end

    heuristic_hit = circular_heuristic && has_terminal_overlap(seq; overlap_bp=circular_overlap_bp)
    should_canonicalize = circular == true ||
        (circular == :auto && record_count == 1 && (header_says_circular(header) || heuristic_hit))
    if !should_canonicalize
        return normalized
    end

    canonical = joinpath(outdir, "canonical.fna")
    return canonicalize_circular_fasta(normalized; out_path=canonical, force=force)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Fragment a genome FASTA into consecutive fixed-size chunks.

# Arguments
- `fasta_path::String`: Path to input FASTA file.
- `fragment_size::Int`: Fragment length in bases (default: 1020).
- `out_path::String`: Output FASTA path for fragments.
- `discard_tail::Bool`: If true, drop trailing fragments shorter than `fragment_size`.

# Returns
- `String`: Path to the fragment FASTA file.
"""
function fragment_genome(fasta_path::String, fragment_size::Int=1020;
        out_path::String, discard_tail::Bool=true)
    @assert isfile(fasta_path) "FASTA file does not exist: $(fasta_path)"
    @assert fragment_size > 0 "fragment_size must be positive: $(fragment_size)"
    @assert !isempty(out_path) "out_path must be a non-empty string"
    mkpath(dirname(out_path))

    open(FASTX.FASTA.Writer, out_path) do writer
        open(fasta_path) do in_io
            reader = FASTX.FASTA.Reader(in_io)
            for record in reader
                seq = FASTX.sequence(record)
                seq_id = FASTX.identifier(record)
                seq_len = length(seq)
                num_fragments = discard_tail ? div(seq_len, fragment_size) : cld(seq_len, fragment_size)

                for i in 1:num_fragments
                    start_pos = (i - 1) * fragment_size + 1
                    end_pos = discard_tail ? start_pos + fragment_size - 1 : min(i * fragment_size, seq_len)
                    frag_seq = seq[start_pos:end_pos]
                    frag_id = "$(seq_id)_frag_$(i)"
                    frag_rec = FASTX.FASTA.Record(frag_id, frag_seq)
                    write(writer, frag_rec)
                end
            end
        end
    end

    return out_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Average Nucleotide Identity (ANI) using the standard fragmentation + BLASTN method.

# Arguments
- `query::String`: Path to query genome FASTA.
- `reference::String`: Path to reference genome FASTA (used to build BLAST db).
- `outdir::String`: Output directory (default: temporary directory).
- `threads::Int`: Number of BLAST threads.
- `fragment_size::Int`: Query fragment length (default: 1020).
- `min_coverage_pct::Float64`: Minimum alignment coverage of the query fragment.
- `min_identity_pct::Float64`: Minimum percent identity threshold.
- `task::String`: BLASTN task (default: "blastn").
- `evalue::Float64`: BLAST E-value threshold.
- `discard_tail::Bool`: Drop trailing fragments shorter than `fragment_size`.
- `coverage_mode::Symbol`: Coverage calculation (`:alignment` or `:qspan`).
- `force::Bool`: Force recomputation of fragments and BLAST results.

# Returns
- `Float64`: The ANI value (0.0 - 100.0).
- `DataFrame`: The filtered BLAST results used for the calculation.
"""
function calculate_gold_standard_ani(;
    query::String,
    reference::String,
    outdir::String = mktempdir(),
    threads::Int = get_default_threads(),
    fragment_size::Int = 1020,
    min_coverage_pct::Float64 = 70.0,
    min_identity_pct::Float64 = 30.0,
    task::String = "blastn",
    evalue::Float64 = 1e-15,
    discard_tail::Bool = true,
    coverage_mode::Symbol = :alignment,
    force::Bool = false
)
    @assert isfile(query) "Query file does not exist: $(query)"
    @assert isfile(reference) "Reference file does not exist: $(reference)"
    mkpath(outdir)

    fragmented_query = joinpath(outdir, "query_fragments.fasta")
    if force || !isfile(fragmented_query) || filesize(fragmented_query) == 0
        fragment_genome(query, fragment_size; out_path=fragmented_query, discard_tail=discard_tail)
    end

    db_prefix = Mycelia.ensure_blast_db(fasta=reference, dbtype="nucl")
    blast_report_path = Mycelia.run_blastn(
        outdir=outdir,
        fasta=fragmented_query,
        blastdb=db_prefix,
        threads=threads,
        task=task,
        force=force,
        max_target_seqs=1,
        evalue=evalue
    )

    df = standardize_blast_hits(blast_report_path)
    if DataFrames.nrow(df) == 0
        @warn "No BLAST hits found."
        return 0.0, df
    end

    df[!, :query_coverage] = coverage_fraction(df; coverage_mode=coverage_mode) .* 100
    filtered_df = DataFrames.filter(row ->
        row[:query_coverage] >= min_coverage_pct &&
        row[:pident] >= min_identity_pct,
        df
    )

    if DataFrames.nrow(filtered_df) == 0
        return 0.0, filtered_df
    end

    ani = Statistics.mean(filtered_df.pident)
    return ani, filtered_df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Average Amino Acid Identity (AAI) using BLASTP best-hit filtering.

# Arguments
- `query_proteins::String`: Path to query protein FASTA.
- `reference_proteins::String`: Path to reference protein FASTA.
- `outdir::String`: Output directory (default: temporary directory).
- `threads::Int`: Number of BLAST threads.
- `min_coverage_pct::Float64`: Minimum alignment coverage of the query protein.
- `min_identity_pct::Float64`: Minimum percent identity threshold.
- `evalue::Float64`: BLAST E-value threshold.
- `max_target_seqs::Int`: Maximum target sequences per query (default: 1).
- `force::Bool`: Force recomputation of BLAST results.

# Returns
- `Float64`: The AAI value (0.0 - 100.0).
- `DataFrame`: The filtered BLAST results used for the calculation.
"""
function calculate_gold_standard_aai(;
    query_proteins::String,
    reference_proteins::String,
    outdir::String = mktempdir(),
    threads::Int = get_default_threads(),
    min_coverage_pct::Float64 = 70.0,
    min_identity_pct::Float64 = 30.0,
    evalue::Float64 = 1e-5,
    max_target_seqs::Int = 1,
    force::Bool = false
)
    @assert isfile(query_proteins) "Query protein file does not exist: $(query_proteins)"
    @assert isfile(reference_proteins) "Reference protein file does not exist: $(reference_proteins)"
    mkpath(outdir)

    blast_report_path = Mycelia.run_blastp_search(
        query_fasta=query_proteins,
        reference_fasta=reference_proteins,
        output_dir=outdir,
        threads=threads,
        evalue=evalue,
        max_target_seqs=max_target_seqs
    )

    df = Mycelia.parse_blast_report(blast_report_path)
    if DataFrames.nrow(df) == 0
        @warn "No BLASTP hits found."
        return 0.0, df
    end

    df[!, :query_coverage] = (df[!, Symbol("alignment length")] ./ df[!, Symbol("query length")]) .* 100
    filtered_df = DataFrames.filter(row ->
        row[:query_coverage] >= min_coverage_pct &&
        row[Symbol("% identity")] >= min_identity_pct,
        df
    )

    if DataFrames.nrow(filtered_df) == 0
        return 0.0, filtered_df
    end

    aai = Statistics.mean(filtered_df[!, Symbol("% identity")])
    return aai, filtered_df
end

function empty_hits_table()
    return DataFrames.DataFrame(
        query_id=String[],
        subject_id=String[],
        pident=Float64[],
        alignment_length=Int[],
        query_length=Int[],
        subject_length=Int[],
        evalue=Float64[],
        bitscore=Float64[]
    )
end

function standardize_blast_hits(blast_report::AbstractString)
    df = Mycelia.parse_blast_report(blast_report)
    if DataFrames.nrow(df) == 0
        return empty_hits_table()
    end
    current_names = names(df)
    if any(name -> name isa AbstractString, current_names)
        DataFrames.rename!(df, current_names .=> Symbol.(current_names))
    end
    rename_map = Dict(
        Symbol("query id") => :query_id,
        Symbol("qseqid") => :query_id,
        Symbol("subject id") => :subject_id,
        Symbol("sseqid") => :subject_id,
        Symbol("alignment length") => :alignment_length,
        Symbol("length") => :alignment_length,
        Symbol("% identity") => :pident,
        Symbol("pident") => :pident,
        Symbol("query length") => :query_length,
        Symbol("qlen") => :query_length,
        Symbol("subject length") => :subject_length,
        Symbol("slen") => :subject_length,
        Symbol("q. start") => :query_start,
        Symbol("qstart") => :query_start,
        Symbol("q. end") => :query_end,
        Symbol("qend") => :query_end,
        Symbol("s. start") => :subject_start,
        Symbol("sstart") => :subject_start,
        Symbol("s. end") => :subject_end,
        Symbol("send") => :subject_end,
        Symbol("evalue") => :evalue,
        Symbol("bit score") => :bitscore,
        Symbol("bitscore") => :bitscore
    )
    for (old_name, new_name) in rename_map
        if old_name in names(df)
            DataFrames.rename!(df, old_name => new_name)
        elseif String(old_name) in names(df)
            DataFrames.rename!(df, String(old_name) => new_name)
        end
    end
    return df
end

function read_diamond_hits(results_file::AbstractString)
    if !isfile(results_file) || filesize(results_file) == 0
        return empty_hits_table()
    end
    header = ["query_id", "subject_id", "pident", "alignment_length",
              "query_length", "subject_length", "evalue", "bitscore"]
    return CSV.read(results_file, DataFrames.DataFrame; delim='\t', header=header, normalizenames=true)
end

function best_hits_by_query(df::DataFrames.DataFrame; query_col::Symbol=:query_id,
        bitscore_col::Symbol=:bitscore, evalue_col::Symbol=:evalue)
    if DataFrames.nrow(df) == 0
        return df
    end
    sorted = DataFrames.sort(df, [query_col, bitscore_col, evalue_col], rev=[false, true, false])
    grouped = DataFrames.groupby(sorted, query_col)
    return DataFrames.combine(grouped) do subdf
        subdf[1, :]
    end
end

function filter_hits_by_threshold(df::DataFrames.DataFrame; min_id::Float64,
        min_aln_frac::Float64, use_shorter::Bool=false, evalue_max::Union{Nothing,Float64}=nothing)
    if DataFrames.nrow(df) == 0
        return df
    end
    denom = use_shorter ? min.(df.query_length, df.subject_length) : df.query_length
    coverage = df.alignment_length ./ denom
    coverage_ok = coverage .>= min_aln_frac
    pident_ok = df.pident .>= min_id
    mask = coalesce.(coverage_ok, false) .& coalesce.(pident_ok, false)
    if !isnothing(evalue_max)
        evalue_ok = df.evalue .<= evalue_max
        mask .&= coalesce.(evalue_ok, false)
    end
    return df[mask, :]
end

function coverage_fraction(df::DataFrames.DataFrame; coverage_mode::Symbol=:alignment)
    @assert coverage_mode in (:alignment, :qspan) "coverage_mode must be :alignment or :qspan"
    if DataFrames.nrow(df) == 0
        return Float64[]
    end
    if coverage_mode == :alignment
        return df.alignment_length ./ df.query_length
    end
    required = (:query_start in names(df)) && (:query_end in names(df))
    @assert required "coverage_mode=:qspan requires query_start and query_end columns"
    return (abs.(df.query_end .- df.query_start) .+ 1) ./ df.query_length
end

function compute_directional_ani(hits_df::DataFrames.DataFrame, total_fragments::Int;
        min_id::Float64, min_aln_frac::Float64, coverage_mode::Symbol=:alignment)
    if total_fragments == 0
        return (;ani=missing, af=missing, hits=DataFrames.DataFrame())
    end
    best_hits = best_hits_by_query(hits_df)
    coverage = coalesce.(coverage_fraction(best_hits; coverage_mode=coverage_mode), 0.0)
    pident = coalesce.(best_hits.pident, -Inf)
    mask = (coverage .>= min_aln_frac) .& (pident .>= min_id)
    filtered = best_hits[mask, :]
    if DataFrames.nrow(filtered) == 0
        return (;ani=missing, af=0.0, hits=filtered)
    end
    ani = Statistics.mean(filtered.pident)
    af = DataFrames.nrow(filtered) / total_fragments
    return (;ani, af, hits=filtered)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute ANIm-style ANI using MUMmer `dnadiff`.
"""
function ani_mummer(reference::AbstractString, query::AbstractString;
        outdir::AbstractString="ani_mummer",
        prefix::String="dnadiff",
        force::Bool=false,
        dnadiff_args::Vector{String}=String[])
    mkpath(outdir)
    dnadiff_outputs = Mycelia.run_dnadiff(
        reference=reference,
        query=query,
        outdir=outdir,
        prefix=prefix,
        force=force,
        additional_args=dnadiff_args
    )
    parsed = Mycelia.parse_dnadiff_report(dnadiff_outputs.report)
    sections = parsed.raw_sections
    one_to_one = get(sections, Symbol("1_to_1"), Dict{Symbol, Any}())
    many_to_many = get(sections, :m_to_m, Dict{Symbol, Any}())
    summary = parsed.summary

    function aligned_fraction(metrics::AbstractDict{Symbol, T}, which::Symbol) where {T}
        pct_key = which == :ref ? :aligned_pct_ref : :aligned_pct_query
        if haskey(metrics, pct_key)
            return metrics[pct_key] / 100.0
        end
        aligned_key = which == :ref ? :aligned_bases_ref : :aligned_bases_query
        total_key = which == :ref ? :total_bases_ref : :total_bases_query
        if haskey(metrics, aligned_key) && haskey(metrics, total_key) && metrics[total_key] > 0
            return metrics[aligned_key] / metrics[total_key]
        end
        return missing
    end

    summary_dict = Dict(pairs(summary))
    fallback_ref = aligned_fraction(summary_dict, :ref)
    fallback_query = aligned_fraction(summary_dict, :query)

    return (;
        method=:ANIm,
        ani_1to1=get(one_to_one, :avg_identity, get(summary, :avg_identity, missing)),
        ani_m2m=get(many_to_many, :avg_identity, missing),
        af_ref_1to1=aligned_fraction(one_to_one, :ref) === missing ? fallback_ref : aligned_fraction(one_to_one, :ref),
        af_query_1to1=aligned_fraction(one_to_one, :query) === missing ? fallback_query : aligned_fraction(one_to_one, :query),
        af_ref_m2m=aligned_fraction(many_to_many, :ref),
        af_query_m2m=aligned_fraction(many_to_many, :query),
        aligned_bases_ref=get(one_to_one, :aligned_bases_ref, get(summary, :aligned_bases_ref, missing)),
        aligned_bases_query=get(one_to_one, :aligned_bases_query, get(summary, :aligned_bases_query, missing)),
        total_bases_ref=get(one_to_one, :total_bases_ref, get(summary, :total_bases_ref, missing)),
        total_bases_query=get(one_to_one, :total_bases_query, get(summary, :total_bases_query, missing)),
        report_path=dnadiff_outputs.report,
        parsed_report=parsed
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute ANIb-style ANI using BLASTn on fixed-length fragments in both directions.

Defaults are set by `mode`: `:strict` (1020 bp, 30% id, evalue 1e-15, tail-discard)
or `:modern` (1000 bp, 70% id, evalue 1e-3).
"""
function ani_blast(reference::AbstractString, query::AbstractString;
        outdir::AbstractString="ani_blast",
        mode::Symbol=:strict,
        fragment_size::Union{Nothing,Int}=nothing,
        discard_tail::Union{Nothing,Bool}=nothing,
        min_id::Union{Nothing,Float64}=nothing,
        min_aln_frac::Union{Nothing,Float64}=nothing,
        coverage_mode::Symbol=:alignment,
        threads::Int=get_default_threads(),
        task::String="blastn",
        evalue::Union{Nothing,Float64}=nothing,
        force::Bool=false)
    @assert mode in (:strict, :modern) "mode must be :strict or :modern"
    default_fragment_size = mode == :strict ? 1020 : 1000
    default_min_id = mode == :strict ? 30.0 : 70.0
    default_min_aln_frac = 0.7
    default_evalue = mode == :strict ? 1e-15 : 1e-3

    fragment_size = isnothing(fragment_size) ? default_fragment_size : fragment_size
    min_id = isnothing(min_id) ? default_min_id : min_id
    min_aln_frac = isnothing(min_aln_frac) ? default_min_aln_frac : min_aln_frac
    evalue = isnothing(evalue) ? default_evalue : evalue
    discard_tail = isnothing(discard_tail) ? (mode == :strict) : discard_tail

    @assert fragment_size > 0 "fragment_size must be positive"
    @assert min_aln_frac > 0 "min_aln_frac must be positive"
    mkpath(outdir)

    ref_fragments = joinpath(outdir, "reference_fragments.fna")
    query_fragments = joinpath(outdir, "query_fragments.fna")
    if force || !isfile(ref_fragments) || filesize(ref_fragments) == 0
        fragment_genome(reference, fragment_size; out_path=ref_fragments, discard_tail=discard_tail)
    end
    if force || !isfile(query_fragments) || filesize(query_fragments) == 0
        fragment_genome(query, fragment_size; out_path=query_fragments, discard_tail=discard_tail)
    end

    ref_db = Mycelia.ensure_blast_db(fasta=reference, dbtype="nucl", db_prefix=joinpath(outdir, "reference_db"), force=force)
    query_db = Mycelia.ensure_blast_db(fasta=query, dbtype="nucl", db_prefix=joinpath(outdir, "query_db"), force=force)

    query_vs_ref = Mycelia.run_blastn(
        outdir=outdir,
        fasta=query_fragments,
        blastdb=ref_db,
        threads=threads,
        task=task,
        force=force,
        max_target_seqs=1,
        evalue=evalue
    )
    ref_vs_query = Mycelia.run_blastn(
        outdir=outdir,
        fasta=ref_fragments,
        blastdb=query_db,
        threads=threads,
        task=task,
        force=force,
        max_target_seqs=1,
        evalue=evalue
    )

    query_hits = standardize_blast_hits(query_vs_ref)
    ref_hits = standardize_blast_hits(ref_vs_query)

    query_frag_count = count_fasta_records(query_fragments)
    ref_frag_count = count_fasta_records(ref_fragments)

    query_summary = compute_directional_ani(query_hits, query_frag_count;
        min_id=min_id, min_aln_frac=min_aln_frac, coverage_mode=coverage_mode)
    ref_summary = compute_directional_ani(ref_hits, ref_frag_count;
        min_id=min_id, min_aln_frac=min_aln_frac, coverage_mode=coverage_mode)

    ani_values = skipmissing([query_summary.ani, ref_summary.ani])
    ani = isempty(ani_values) ? missing : Statistics.mean(collect(ani_values))

    return (;
        method=:ANIb,
        ani=ani,
        ani_query_vs_ref=query_summary.ani,
        ani_ref_vs_query=ref_summary.ani,
        af_query_vs_ref=query_summary.af,
        af_ref_vs_query=ref_summary.af,
        n_query_fragments=query_frag_count,
        n_reference_fragments=ref_frag_count,
        params=(; mode, fragment_size, min_id, min_aln_frac, coverage_mode, evalue, discard_tail),
        query_fragments=query_fragments,
        reference_fragments=ref_fragments,
        hits_query_vs_ref=query_summary.hits,
        hits_ref_vs_query=ref_summary.hits,
        raw_query_hits=query_hits,
        raw_ref_hits=ref_hits
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute AAI using reciprocal best hits (RBH) between predicted or provided proteins.

If `tool=:auto`, DIAMOND is used when available, otherwise BLASTP.
"""
function aai_rbh(genome_a::AbstractString, genome_b::AbstractString;
        proteins_a::Union{Nothing,AbstractString}=nothing,
        proteins_b::Union{Nothing,AbstractString}=nothing,
        tool::Symbol=:auto,
        min_id::Float64=30.0,
        min_len_frac::Float64=0.7,
        evalue::Float64=1e-5,
        coverage_denom::Symbol=:shorter,
        threads::Int=get_default_threads(),
        outdir::AbstractString="aai_rbh",
        translation_table::Union{Nothing,Int}=nothing,
        force::Bool=false)
    mkpath(outdir)
    @assert coverage_denom in (:shorter, :query) "coverage_denom must be :shorter or :query"

    prot_a = isnothing(proteins_a) ? Mycelia.run_pyrodigal(
        fasta_file=genome_a,
        out_dir=joinpath(outdir, "proteins_a"),
        translation_table=translation_table
    ).faa : String(proteins_a)
    prot_b = isnothing(proteins_b) ? Mycelia.run_pyrodigal(
        fasta_file=genome_b,
        out_dir=joinpath(outdir, "proteins_b"),
        translation_table=translation_table
    ).faa : String(proteins_b)

    resolved_tool = tool == :auto ? :diamond : tool

    hits_ab_path = ""
    hits_ba_path = ""
    if resolved_tool == :diamond
        try
            hits_ab_path = Mycelia.run_diamond_besthits(
                query_fasta=prot_a,
                reference_fasta=prot_b,
                output_dir=joinpath(outdir, "diamond_ab"),
                threads=threads,
                evalue=evalue,
                force=force
            )
            hits_ba_path = Mycelia.run_diamond_besthits(
                query_fasta=prot_b,
                reference_fasta=prot_a,
                output_dir=joinpath(outdir, "diamond_ba"),
                threads=threads,
                evalue=evalue,
                force=force
            )
        catch e
            if tool == :auto
                @warn "DIAMOND failed for AAI; falling back to BLASTP" exception=e
                resolved_tool = :blastp
            else
                rethrow(e)
            end
        end
    end
    if resolved_tool == :blastp
        hits_ab_path = Mycelia.run_blastp_search(
            query_fasta=prot_a,
            reference_fasta=prot_b,
            output_dir=joinpath(outdir, "blastp_ab"),
            threads=threads,
            evalue=evalue,
            max_target_seqs=1
        )
        hits_ba_path = Mycelia.run_blastp_search(
            query_fasta=prot_b,
            reference_fasta=prot_a,
            output_dir=joinpath(outdir, "blastp_ba"),
            threads=threads,
            evalue=evalue,
            max_target_seqs=1
        )
    elseif resolved_tool != :diamond
        error("tool must be :auto, :diamond, or :blastp")
    end

    hits_ab = resolved_tool == :diamond ? read_diamond_hits(hits_ab_path) : standardize_blast_hits(hits_ab_path)
    hits_ba = resolved_tool == :diamond ? read_diamond_hits(hits_ba_path) : standardize_blast_hits(hits_ba_path)

    best_ab = best_hits_by_query(hits_ab)
    best_ba = best_hits_by_query(hits_ba)

    use_shorter = coverage_denom == :shorter
    filtered_ab = filter_hits_by_threshold(best_ab; min_id=min_id, min_aln_frac=min_len_frac, use_shorter=use_shorter, evalue_max=evalue)
    filtered_ba = filter_hits_by_threshold(best_ba; min_id=min_id, min_aln_frac=min_len_frac, use_shorter=use_shorter, evalue_max=evalue)

    renamed_ab = DataFrames.rename(filtered_ab, Dict(
        :query_id => :query_a,
        :subject_id => :subject_b,
        :pident => :pident_ab
    ))
    renamed_ba = DataFrames.rename(filtered_ba, Dict(
        :query_id => :query_b,
        :subject_id => :subject_a,
        :pident => :pident_ba
    ))
    rbh_table = DataFrames.innerjoin(
        DataFrames.select(renamed_ab, [:query_a, :subject_b, :pident_ab]),
        DataFrames.select(renamed_ba, [:query_b, :subject_a, :pident_ba]),
        on=[:query_a => :subject_a, :subject_b => :query_b]
    )

    n_rbh = DataFrames.nrow(rbh_table)
    n_a = count_fasta_records(prot_a)
    n_b = count_fasta_records(prot_b)
    ortholog_fraction = min(n_a, n_b) == 0 ? missing : n_rbh / min(n_a, n_b)
    aai = n_rbh == 0 ? missing : Statistics.mean((rbh_table.pident_ab .+ rbh_table.pident_ba) ./ 2)

    return (;
        method=:AAI_RBH,
        aai=aai,
        ortholog_fraction=ortholog_fraction,
        n_rbh=n_rbh,
        n_proteins_a=n_a,
        n_proteins_b=n_b,
        params=(; tool=resolved_tool, min_id, min_len_frac, evalue, coverage_denom),
        hits_ab=hits_ab_path,
        hits_ba=hits_ba_path,
        rbh_table=rbh_table
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute Percentage of Conserved Proteins (POCP) between two genomes.

Returns classic POCP (any qualifying hit), plus POCPu variants based on
best-hit uniqueness and reciprocal best hits.
"""
function pocp(genome_a::AbstractString, genome_b::AbstractString;
        proteins_a::Union{Nothing,AbstractString}=nothing,
        proteins_b::Union{Nothing,AbstractString}=nothing,
        tool::Symbol=:auto,
        min_id::Float64=40.0,
        min_len_frac::Float64=0.5,
        evalue::Float64=1e-5,
        coverage_denom::Symbol=:query,
        max_target_seqs::Int=500,
        threads::Int=get_default_threads(),
        outdir::AbstractString="pocp",
        translation_table::Union{Nothing,Int}=nothing,
        force::Bool=false)
    @assert coverage_denom in (:shorter, :query) "coverage_denom must be :shorter or :query"
    mkpath(outdir)

    prot_a = isnothing(proteins_a) ? Mycelia.run_pyrodigal(
        fasta_file=genome_a,
        out_dir=joinpath(outdir, "proteins_a"),
        translation_table=translation_table
    ).faa : String(proteins_a)
    prot_b = isnothing(proteins_b) ? Mycelia.run_pyrodigal(
        fasta_file=genome_b,
        out_dir=joinpath(outdir, "proteins_b"),
        translation_table=translation_table
    ).faa : String(proteins_b)

    resolved_tool = tool == :auto ? :diamond : tool
    hits_ab_path = ""
    hits_ba_path = ""
    if resolved_tool == :diamond
        try
            hits_ab_path = Mycelia.run_diamond_besthits(
                query_fasta=prot_a,
                reference_fasta=prot_b,
                output_dir=joinpath(outdir, "diamond_ab"),
                threads=threads,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                force=force
            )
            hits_ba_path = Mycelia.run_diamond_besthits(
                query_fasta=prot_b,
                reference_fasta=prot_a,
                output_dir=joinpath(outdir, "diamond_ba"),
                threads=threads,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                force=force
            )
        catch e
            if tool == :auto
                @warn "DIAMOND failed for POCP; falling back to BLASTP" exception=e
                resolved_tool = :blastp
            else
                rethrow(e)
            end
        end
    end
    if resolved_tool == :blastp
        hits_ab_path = Mycelia.run_blastp_search(
            query_fasta=prot_a,
            reference_fasta=prot_b,
            output_dir=joinpath(outdir, "blastp_ab"),
            threads=threads,
            evalue=evalue,
            max_target_seqs=max_target_seqs
        )
        hits_ba_path = Mycelia.run_blastp_search(
            query_fasta=prot_b,
            reference_fasta=prot_a,
            output_dir=joinpath(outdir, "blastp_ba"),
            threads=threads,
            evalue=evalue,
            max_target_seqs=max_target_seqs
        )
    elseif resolved_tool != :diamond
        error("tool must be :auto, :diamond, or :blastp")
    end

    hits_ab = resolved_tool == :diamond ? read_diamond_hits(hits_ab_path) : standardize_blast_hits(hits_ab_path)
    hits_ba = resolved_tool == :diamond ? read_diamond_hits(hits_ba_path) : standardize_blast_hits(hits_ba_path)

    n_a = count_fasta_records(prot_a)
    n_b = count_fasta_records(prot_b)
    total_proteins = n_a + n_b
    if total_proteins == 0
        return (;
            method=:POCP,
            pocp=missing,
            pocpu_besthit=missing,
            pocpu_rbh=missing,
            c1=0,
            c2=0,
            c1_besthit=0,
            c2_besthit=0,
            n_rbh=0,
            n_proteins_a=n_a,
            n_proteins_b=n_b,
            params=(; tool=resolved_tool, min_id, min_len_frac, evalue, coverage_denom, max_target_seqs),
            hits_ab=hits_ab_path,
            hits_ba=hits_ba_path,
            rbh_table=DataFrames.DataFrame()
        )
    end

    use_shorter = coverage_denom == :shorter
    filtered_ab_any = filter_hits_by_threshold(hits_ab; min_id=min_id, min_aln_frac=min_len_frac, use_shorter=use_shorter, evalue_max=evalue)
    filtered_ba_any = filter_hits_by_threshold(hits_ba; min_id=min_id, min_aln_frac=min_len_frac, use_shorter=use_shorter, evalue_max=evalue)
    c1 = length(unique(filtered_ab_any.query_id))
    c2 = length(unique(filtered_ba_any.query_id))
    pocp_value = ((c1 + c2) / total_proteins) * 100

    best_ab = best_hits_by_query(hits_ab)
    best_ba = best_hits_by_query(hits_ba)
    filtered_ab_best = filter_hits_by_threshold(best_ab; min_id=min_id, min_aln_frac=min_len_frac, use_shorter=use_shorter, evalue_max=evalue)
    filtered_ba_best = filter_hits_by_threshold(best_ba; min_id=min_id, min_aln_frac=min_len_frac, use_shorter=use_shorter, evalue_max=evalue)
    filtered_ab_best = filtered_ab_best[:, [:query_id, :subject_id, :pident]]
    filtered_ba_best = filtered_ba_best[:, [:query_id, :subject_id, :pident]]
    c1_best = length(unique(filtered_ab_best.query_id))
    c2_best = length(unique(filtered_ba_best.query_id))
    pocpu_besthit = ((c1_best + c2_best) / total_proteins) * 100

    renamed_ab = DataFrames.rename(filtered_ab_best, Dict(
        :query_id => :query_a,
        :subject_id => :subject_b,
        :pident => :pident_ab
    ))
    renamed_ba = DataFrames.rename(filtered_ba_best, Dict(
        :query_id => :query_b,
        :subject_id => :subject_a,
        :pident => :pident_ba
    ))
    rbh_table = DataFrames.innerjoin(
        renamed_ab,
        renamed_ba,
        on=[:query_a => :subject_a, :subject_b => :query_b]
    )
    n_rbh = DataFrames.nrow(rbh_table)
    pocpu_rbh = ((2 * n_rbh) / total_proteins) * 100

    return (;
        method=:POCP,
        pocp=pocp_value,
        pocpu_besthit=pocpu_besthit,
        pocpu_rbh=pocpu_rbh,
        c1=c1,
        c2=c2,
        c1_besthit=c1_best,
        c2_besthit=c2_best,
        n_rbh=n_rbh,
        n_proteins_a=n_a,
        n_proteins_b=n_b,
        params=(; tool=resolved_tool, min_id, min_len_frac, evalue, coverage_denom, max_target_seqs),
        hits_ab=hits_ab_path,
        hits_ba=hits_ba_path,
        rbh_table=rbh_table
    )
end

function flatten_gold_comparison(results::Dict{Symbol, Any})
    fields = Dict{Symbol, Any}()
    if haskey(results, :ANIm)
        anim = results[:ANIm]
        fields[:anim_ani_1to1] = anim.ani_1to1
        fields[:anim_ani_m2m] = anim.ani_m2m
        fields[:anim_af_ref_1to1] = anim.af_ref_1to1
        fields[:anim_af_query_1to1] = anim.af_query_1to1
        fields[:anim_af_ref_m2m] = anim.af_ref_m2m
        fields[:anim_af_query_m2m] = anim.af_query_m2m
    end
    if haskey(results, :ANIb)
        anib = results[:ANIb]
        fields[:anib_ani] = anib.ani
        fields[:anib_ani_query_vs_ref] = anib.ani_query_vs_ref
        fields[:anib_ani_ref_vs_query] = anib.ani_ref_vs_query
        fields[:anib_af_query_vs_ref] = anib.af_query_vs_ref
        fields[:anib_af_ref_vs_query] = anib.af_ref_vs_query
    end
    if haskey(results, :AAI)
        aai = results[:AAI]
        fields[:aai_rbh] = aai.aai
        fields[:ortholog_fraction] = aai.ortholog_fraction
        fields[:aai_rbh_n] = aai.n_rbh
    end
    if haskey(results, :POCP)
        pocp = results[:POCP]
        fields[:pocp] = pocp.pocp
        fields[:pocpu_besthit] = pocp.pocpu_besthit
        fields[:pocpu_rbh] = pocp.pocpu_rbh
        fields[:pocp_c1] = pocp.c1
        fields[:pocp_c2] = pocp.c2
        fields[:pocp_t1] = pocp.n_proteins_a
        fields[:pocp_t2] = pocp.n_proteins_b
        fields[:pocpu_rbh_n] = pocp.n_rbh
    end
    if haskey(results, :PyOrthoANI)
        orthoani = results[:PyOrthoANI]
        fields[:pyorthoani_ani] = orthoani.ani
    end
    return (; (k => fields[k] for k in sort(collect(keys(fields))))...)
end

function gold_comparison_tool_versions(methods::Vector{Symbol};
        aai_tool::Symbol,
        pocp_tool::Symbol)
    versions = Dict{Symbol, Any}()
    if :ANIb in methods || :AAI in methods || :POCP in methods
        Mycelia.add_bioconda_env("blast")
        versions[:blast] = conda_tool_version("blast", ["blastn", "-version"])
    end
    if (:AAI in methods && aai_tool == :diamond) || (:POCP in methods && pocp_tool == :diamond)
        Mycelia.add_bioconda_env("diamond")
        versions[:diamond] = conda_tool_version("diamond", ["diamond", "version"])
    end
    if :ANIm in methods
        Mycelia.add_bioconda_env("mummer")
        versions[:mummer] = conda_tool_version("mummer", ["nucmer", "--version"])
    end
    if :PyOrthoANI in methods
        ensure_pyorthoani_env(quiet=true)
        versions[:pyorthoani] = conda_tool_version(
            "pyorthoani",
            ["python", "-c", "import pyorthoani; print(pyorthoani.__version__)"]
        )
    end
    return versions
end

function gold_comparison_cache_key(pair_id::AbstractString;
        methods::Vector{Symbol},
        params::NamedTuple,
        tool_versions::Dict{Symbol, Any})
    parts = String[]
    push!(parts, "pair_id=$(pair_id)")
    push!(parts, "methods=$(join(sort(string.(methods)), ","))")
    for (key, value) in sort(collect(pairs(params)); by=first)
        push!(parts, "$(key)=$(string(value))")
    end
    for (key, value) in sort(collect(pairs(tool_versions)); by=first)
        push!(parts, "$(key)=$(string(value))")
    end
    return SHA.bytes2hex(SHA.sha256(join(parts, "|")))
end

function write_gold_comparison_report(report_dir::AbstractString;
        pair_id::AbstractString,
        reference_prepped::AbstractString,
        query_prepped::AbstractString,
        methods::Vector{Symbol},
        params::NamedTuple,
        tool_versions::Dict{Symbol, Any},
        cache_key::AbstractString,
        summary::NamedTuple)
    mkpath(report_dir)
    summary_json = joinpath(report_dir, "summary.json")
    summary_tsv = joinpath(report_dir, "summary.tsv")
    metadata_json = joinpath(report_dir, "metadata.json")

    summary_dict = Dict(string(k) => normalize_json_value(v) for (k, v) in pairs(summary))
    open(summary_json, "w") do io
        JSON.print(io, summary_dict, 2)
    end

    open(summary_tsv, "w") do io
        write(io, "metric\tvalue\n")
        for key in sort(collect(keys(summary_dict)))
            value = summary_dict[key]
            value_str = value === nothing ? "NA" : string(value)
            write(io, "$(key)\t$(value_str)\n")
        end
    end

    metadata = Dict(
        "pair_id" => pair_id,
        "cache_key" => cache_key,
        "created_at_unix" => time(),
        "methods" => sort(string.(methods)),
        "params" => normalize_json_value(params),
        "tool_versions" => normalize_json_value(tool_versions),
        "reference_prepped" => reference_prepped,
        "query_prepped" => query_prepped
    )
    open(metadata_json, "w") do io
        JSON.print(io, metadata, 2)
    end

    return (;summary_json, summary_tsv, metadata_json)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run gold-standard genome comparison methods and return a consolidated summary.

Available methods include `:ANIm`, `:ANIb`, `:AAI`, `:POCP`, and `:PyOrthoANI`.
"""
function compare_genomes_gold(reference::AbstractString, query::AbstractString;
        circular::Union{Bool,Symbol}=:auto,
        circular_heuristic::Bool=false,
        circular_overlap_bp::Int=200,
        methods=Symbol[:ANIm, :ANIb, :AAI, :POCP, :PyOrthoANI],
        outdir::AbstractString="genome_compare",
        threads::Int=get_default_threads(),
        use_cache::Bool=true,
        ani_mode::Symbol=:strict,
        fragment_size::Union{Nothing,Int}=nothing,
        min_ani_id::Union{Nothing,Float64}=nothing,
        min_ani_aln_frac::Union{Nothing,Float64}=nothing,
        ani_coverage_mode::Symbol=:alignment,
        ani_discard_tail::Union{Nothing,Bool}=nothing,
        ani_evalue::Union{Nothing,Float64}=nothing,
        aai_tool::Symbol=:auto,
        min_aai_id::Float64=30.0,
        min_aai_len_frac::Float64=0.7,
        aai_coverage_denom::Symbol=:shorter,
        aai_evalue::Float64=1e-5,
        pocp_tool::Symbol=:auto,
        pocp_min_id::Float64=40.0,
        pocp_min_len_frac::Float64=0.5,
        pocp_coverage_denom::Symbol=:query,
        pocp_evalue::Float64=1e-5,
        pocp_max_target_seqs::Int=500,
        pyorthoani_force_env::Bool=false,
        pyorthoani_args::Vector{String}=String[],
        pyorthoani_quiet::Bool=false,
        translation_table::Union{Nothing,Int}=nothing,
        force::Bool=false)
    pair_id = genome_pair_id(reference, query)
    methods = sort(unique(Symbol.(methods)))
    report_dir = joinpath(outdir, "report", pair_id)
    summary_json = joinpath(report_dir, "summary.json")
    metadata_json = joinpath(report_dir, "metadata.json")

    resolved_aai_tool = aai_tool
    if :AAI in methods && aai_tool == :auto
        Mycelia.add_bioconda_env("diamond")
        resolved_aai_tool = _conda_env_exec_ok("diamond", ["diamond", "version"]) ? :diamond : :blastp
    end
    resolved_pocp_tool = pocp_tool
    if :POCP in methods && pocp_tool == :auto
        Mycelia.add_bioconda_env("diamond")
        resolved_pocp_tool = _conda_env_exec_ok("diamond", ["diamond", "version"]) ? :diamond : :blastp
    end
    pre_params = (;
        circular,
        circular_heuristic,
        circular_overlap_bp,
        ani_mode,
        fragment_size,
        min_ani_id,
        min_ani_aln_frac,
        ani_coverage_mode,
        ani_discard_tail,
        ani_evalue,
        aai_tool=resolved_aai_tool,
        min_aai_id,
        min_aai_len_frac,
        aai_coverage_denom,
        aai_evalue,
        pocp_tool=resolved_pocp_tool,
        pocp_min_id,
        pocp_min_len_frac,
        pocp_coverage_denom,
        pocp_evalue,
        pocp_max_target_seqs,
        pyorthoani_force_env,
        pyorthoani_args,
        pyorthoani_quiet,
        translation_table
    )
    pre_versions = gold_comparison_tool_versions(methods;
        aai_tool=resolved_aai_tool,
        pocp_tool=resolved_pocp_tool)
    pre_cache_key = gold_comparison_cache_key(pair_id;
        methods=methods,
        params=pre_params,
        tool_versions=pre_versions)

    if use_cache && !force && isfile(summary_json) && isfile(metadata_json)
        cached_meta = JSON.parsefile(metadata_json)
        if get(cached_meta, "cache_key", "") == pre_cache_key
            cached_summary = JSON.parsefile(summary_json)
            summary = (; (Symbol(k) => v for (k, v) in cached_summary)...)
            reference_prepped = get(cached_meta, "reference_prepped", reference)
            query_prepped = get(cached_meta, "query_prepped", query)
            return (;pair_id, reference=reference_prepped, query=query_prepped, summary,
                details=Dict{Symbol, Any}(), cache_hit=true, report_dir, cache_key=pre_cache_key)
        end
    end
    prep_dir = joinpath(outdir, "prep", pair_id)
    ref_prepped = prepare_genome_for_comparison(
        reference;
        outdir=joinpath(prep_dir, "reference"),
        circular=circular,
        circular_heuristic=circular_heuristic,
        circular_overlap_bp=circular_overlap_bp,
        force=force
    )
    query_prepped = prepare_genome_for_comparison(
        query;
        outdir=joinpath(prep_dir, "query"),
        circular=circular,
        circular_heuristic=circular_heuristic,
        circular_overlap_bp=circular_overlap_bp,
        force=force
    )

    results = Dict{Symbol, Any}()
    if :ANIm in methods
        results[:ANIm] = ani_mummer(
            ref_prepped,
            query_prepped;
            outdir=joinpath(outdir, "anim", pair_id),
            force=force
        )
    end
    if :ANIb in methods
        results[:ANIb] = ani_blast(
            ref_prepped,
            query_prepped;
            outdir=joinpath(outdir, "anib", pair_id),
            mode=ani_mode,
            fragment_size=fragment_size,
            discard_tail=ani_discard_tail,
            min_id=min_ani_id,
            min_aln_frac=min_ani_aln_frac,
            coverage_mode=ani_coverage_mode,
            evalue=ani_evalue,
            threads=threads,
            force=force
        )
    end
    if :AAI in methods
        results[:AAI] = aai_rbh(
            ref_prepped,
            query_prepped;
            tool=aai_tool,
            min_id=min_aai_id,
            min_len_frac=min_aai_len_frac,
            evalue=aai_evalue,
            coverage_denom=aai_coverage_denom,
            threads=threads,
            outdir=joinpath(outdir, "aai", pair_id),
            translation_table=translation_table,
            force=force
        )
    end
    if :POCP in methods
        results[:POCP] = pocp(
            ref_prepped,
            query_prepped;
            tool=pocp_tool,
            min_id=pocp_min_id,
            min_len_frac=pocp_min_len_frac,
            evalue=pocp_evalue,
            coverage_denom=pocp_coverage_denom,
            max_target_seqs=pocp_max_target_seqs,
            threads=threads,
            outdir=joinpath(outdir, "pocp", pair_id),
            translation_table=translation_table,
            force=force
        )
    end
    if :PyOrthoANI in methods
        results[:PyOrthoANI] = run_pyorthoani(
            query=query_prepped,
            reference=ref_prepped,
            outdir=joinpath(outdir, "pyorthoani", pair_id),
            force=force,
            force_env=pyorthoani_force_env,
            additional_args=pyorthoani_args,
            quiet=pyorthoani_quiet
        )
    end

    summary = flatten_gold_comparison(results)
    aai_tool_used = haskey(results, :AAI) ? results[:AAI].params.tool : resolved_aai_tool
    pocp_tool_used = haskey(results, :POCP) ? results[:POCP].params.tool : resolved_pocp_tool
    final_params = (;
        circular,
        circular_heuristic,
        circular_overlap_bp,
        ani_mode,
        fragment_size,
        min_ani_id,
        min_ani_aln_frac,
        ani_coverage_mode,
        ani_discard_tail,
        ani_evalue,
        aai_tool=aai_tool_used,
        min_aai_id,
        min_aai_len_frac,
        aai_coverage_denom,
        aai_evalue,
        pocp_tool=pocp_tool_used,
        pocp_min_id,
        pocp_min_len_frac,
        pocp_coverage_denom,
        pocp_evalue,
        pocp_max_target_seqs,
        pyorthoani_force_env,
        pyorthoani_args,
        pyorthoani_quiet,
        translation_table
    )
    tool_versions = gold_comparison_tool_versions(methods;
        aai_tool=aai_tool_used,
        pocp_tool=pocp_tool_used)
    cache_key = gold_comparison_cache_key(pair_id;
        methods=methods,
        params=final_params,
        tool_versions=tool_versions)
    write_gold_comparison_report(report_dir;
        pair_id=pair_id,
        reference_prepped=ref_prepped,
        query_prepped=query_prepped,
        methods=methods,
        params=final_params,
        tool_versions=tool_versions,
        cache_key=cache_key,
        summary=summary)

    return (;pair_id, reference=ref_prepped, query=query_prepped, summary, details=results,
        cache_hit=false, report_dir, cache_key)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run fastani with a query and reference list

Calculate Average Nucleotide Identity (ANI) between genome sequences using FastANI.

# Arguments
- `query_list::String`: Path to file containing list of query genome paths (one per line)
- `reference_list::String`: Path to file containing list of reference genome paths (one per line)
- `outfile::String`: Path to output file that will contain ANI results
- `threads::Int=get_default_threads()`: Number of parallel threads to use
- `force::Bool=false`: If true, rerun analysis even if output file exists

# Output
Generates a tab-delimited file with columns:
- Query genome
- Reference genome  
- ANI value (%)
- Count of bidirectional fragment mappings
- Total query fragments

# Notes
- Requires FastANI to be available via Bioconda
- Automatically sets up required conda environment
"""
function fastani_list(;query_list="", reference_list="", outfile="", threads=get_default_threads(), force=false)
    Mycelia.add_bioconda_env("fastani")
    if !isfile(outfile) || force
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n fastani fastANI --ql $(query_list) --rl $(reference_list) --threads $(threads) -o $(outfile)`)
        # run(
        # pipeline(
        #     `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastani fastANI --ql $(query_list) --rl $(reference_list) --threads $(threads) -o $(outfile)`,
        #     stdout=outfile * "fastani.stdout.txt",
        #     stderr=outfile * "fastani.stderr.txt"
        #     )
        # )
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Average Nucleotide Identity (ANI) between two genomes using FastANI.

# Arguments
- `query::String`: Path to query genome FASTA file
- `reference::String`: Path to reference genome FASTA file  
- `outfile::String`: Path to save FastANI results
- `force::Bool=false`: If true, overwrite existing output file

# Notes
- Requires FastANI to be available via Bioconda
- Stdout and stderr are captured in separate files with '.stdout.txt' and '.stderr.txt' suffixes
- ANI results are written to the specified outfile
"""
# ./fastANI -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
function fastani_pair(;query="", reference="", outfile="", force=false)
    Mycelia.add_bioconda_env("fastani")
    if !isfile(outfile) || force
        run(
        pipeline(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastani fastANI -q $(query) -r $(reference) -o $(outfile)`,
            stdout=outfile * "fastani.stdout.txt",
            stderr=outfile * "fastani.stderr.txt"
            )
        )
    end
end

function _conda_env_exec_ok(env_name::String, exec_parts::Vector{String})
    cmd = Cmd(vcat([Mycelia.CONDA_RUNNER, "run", "-n", env_name], exec_parts))
    return Base.success(pipeline(cmd, stdout=Base.devnull, stderr=Base.devnull))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Ensure the PyOrthoANI conda environment is available with BLAST+.
"""
function ensure_pyorthoani_env(; force::Bool=false, quiet::Bool=false)
    env_name = "pyorthoani"
    base_packages = ["python=3.11", "pip", "blast"]

    if force && Mycelia.check_bioconda_env_is_installed(env_name)
        run(`$(Mycelia.CONDA_RUNNER) env remove -n $(env_name) -y`)
    end

    if !Mycelia.check_bioconda_env_is_installed(env_name) || force
        cmd_parts = [
            Mycelia.CONDA_RUNNER, "create",
            "-c", "conda-forge", "-c", "bioconda", "-c", "defaults",
            "--strict-channel-priority",
            "-n", env_name
        ]
        append!(cmd_parts, base_packages)
        push!(cmd_parts, "-y")
        if quiet
            push!(cmd_parts, "--quiet")
        end
        run(Cmd(cmd_parts))

        pip_parts = [Mycelia.CONDA_RUNNER, "run", "--live-stream", "-n", env_name, "python", "-m", "pip", "install", "pyorthoani"]
        if quiet
            push!(pip_parts, "-q")
        end
        run(Cmd(pip_parts))

        clean_parts = [Mycelia.CONDA_RUNNER, "clean", "--all", "-y"]
        if quiet
            push!(clean_parts, "--quiet")
        end
        run(Cmd(clean_parts))
        return env_name
    end

    missing = String[]
    if !_conda_env_exec_ok(env_name, ["blastn", "-version"])
        push!(missing, "blast")
    end
    if !_conda_env_exec_ok(env_name, ["python", "-m", "pip", "--version"])
        push!(missing, "pip")
    end

    if !isempty(missing)
        cmd_parts = [
            Mycelia.CONDA_RUNNER, "install",
            "-c", "conda-forge", "-c", "bioconda", "-c", "defaults",
            "--strict-channel-priority",
            "-n", env_name
        ]
        append!(cmd_parts, missing)
        push!(cmd_parts, "-y")
        if quiet
            push!(cmd_parts, "--quiet")
        end
        run(Cmd(cmd_parts))

        clean_parts = [Mycelia.CONDA_RUNNER, "clean", "--all", "-y"]
        if quiet
            push!(clean_parts, "--quiet")
        end
        run(Cmd(clean_parts))
    end

    if !_conda_env_exec_ok(env_name, ["python", "-c", "import pyorthoani"])
        pip_parts = [Mycelia.CONDA_RUNNER, "run", "--live-stream", "-n", env_name, "python", "-m", "pip", "install", "pyorthoani"]
        if quiet
            push!(pip_parts, "-q")
        end
        run(Cmd(pip_parts))
    end

    return env_name
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse the ANI value from PyOrthoANI output.
"""
function parse_pyorthoani_output(output::AbstractString)
    cleaned = strip(output)
    if isempty(cleaned)
        error("PyOrthoANI output was empty.")
    end
    number_match = match(r"[-+]?[0-9]*\.?[0-9]+", cleaned)
    if number_match === nothing
        error("No ANI value found in PyOrthoANI output: $(cleaned)")
    end
    return parse(Float64, number_match.match)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run PyOrthoANI on a query/reference pair and return the ANI value.

# Arguments
- `query::String`: Query FASTA path.
- `reference::String`: Reference FASTA path.
- `outdir::String`: Output directory for stored results.
- `outfile::String`: Optional explicit output file path.
- `force::Bool=false`: Rerun PyOrthoANI if the output file exists.
- `force_env::Bool=false`: Recreate the PyOrthoANI environment.
- `additional_args::Vector{String}`: Extra CLI args passed to `pyorthoani`.
- `quiet::Bool=false`: Reduce conda output during environment setup.

# Returns
`NamedTuple` with `ani`, `output_path`, and `raw_output`.
"""
function run_pyorthoani(;
        query::String,
        reference::String,
        outdir::String=mktempdir(),
        outfile::String="",
        force::Bool=false,
        force_env::Bool=false,
        additional_args::Vector{String}=String[],
        quiet::Bool=false)
    @assert isfile(query) "Query FASTA not found: $(query)"
    @assert isfile(reference) "Reference FASTA not found: $(reference)"

    ensure_pyorthoani_env(force=force_env, quiet=quiet)

    outdir = mkpath(outdir)
    output_path = isempty(outfile) ? joinpath(
        outdir,
        "pyorthoani_$(basename(query))_vs_$(basename(reference)).txt"
    ) : outfile

    if force || !isfile(output_path) || filesize(output_path) == 0
        cmd_parts = ["pyorthoani", "-q", query, "-r", reference]
        append!(cmd_parts, additional_args)
        cmd = Cmd(vcat([Mycelia.CONDA_RUNNER, "run", "--live-stream", "-n", "pyorthoani"], cmd_parts))
        output = read(cmd, String)
        write(output_path, output)
    end

    output = read(output_path, String)
    ani = parse_pyorthoani_output(output)
    return (;ani, output_path, raw_output=output)
end
