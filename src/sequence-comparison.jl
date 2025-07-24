
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
        if BioSequences.alphabet(seq) isa BioSequences.DNAAlphabet
            kmer_type = Kmers.DNAKmer{k}
            counts = count_canonical_kmers(kmer_type, seq)
        elseif BioSequences.alphabet(seq) isa BioSequences.RNAAlphabet
            kmer_type = Kmers.RNAKmer{k}
            counts = count_kmers(kmer_type, seq)
        elseif BioSequences.alphabet(seq) isa BioSequences.AminoAcidAlphabet
            kmer_type = Kmers.AAKmer{k}
            counts = count_kmers(kmer_type, seq)
        else
            error("Unsupported sequence alphabet: $(BioSequences.alphabet(seq))")
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
        if BioSequences.alphabet(seq) isa BioSequences.DNAAlphabet
            kmer_type = Kmers.DNAKmer{k}
            counts = count_canonical_kmers(kmer_type, seq)
        elseif BioSequences.alphabet(seq) isa BioSequences.RNBioSequences.DNAAlphabetlphabet
            kmer_type = Kmers.RNAKmer{k}
            counts = count_kmers(kmer_type, seq)
        elseif BioSequences.alphabet(seq) isa BioSequences.AminoAcidAlphabet
            kmer_type = Kmers.AAKmer{k}
            counts = count_kmers(kmer_type, seq)
        else
            error("Unsupported sequence alphabet: $(BioSequences.alphabet(seq))")
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
        if BioSequences.alphabet(seq) isa BioSequences.DNAAlphabet
            kmer_type = Kmers.DNAKmer{k}
            kmer_counts = count_canonical_kmers(kmer_type, seq)
            A = isnothing(alphabet) ? 4 : alphabet
        elseif BioSequences.alphabet(seq) isa BioSequences.RNAAlphabet
            kmer_type = Kmers.RNAKmer{k}
            kmer_counts = count_kmers(kmer_type, seq)
            A = isnothing(alphabet) ? 4 : alphabet
        elseif BioSequences.alphabet(seq) isa BioSequences.AminoAcidAlphabet
            kmer_type = Kmers.AAKmer{k}
            kmer_counts = count_kmers(kmer_type, seq)
            A = isnothing(alphabet) ? 20 : alphabet
        else
            error("Unsupported sequence alphabet: $(BioSequences.alphabet(seq))")
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
function linguistic_complexity(seq; kmax=nothing, alphabet=nothing, reducer=mean)
    # Determine sequence type and properties
    if seq isa BioSequences.BioSequence
        L = length(seq)
        # Infer alphabet size from sequence type if not provided
        if isnothing(alphabet)
            if BioSequences.alphabet(seq) isa BioSequences.DNAAlphabet
                A = 4
            elseif BioSequences.alphabet(seq) isa BioSequences.RNAAlphabet
                A = 4
            elseif BioSequences.alphabet(seq) isa BioSequences.AminoAcidAlphabet
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
    end

    return (
        # results_table = results_table,
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

    # --- Step 2: Create sketch files for each genome ---
    sketch1 = fasta1 * ".msh"
    sketch2 = fasta2 * ".msh"

    println("Sketching $fasta1 (k=$k, s=$s)...")
    sketch_cmd1 = pipeline(`$mash_path sketch -k $k -s $s -o $sketch1 $fasta1`, stdout=devnull, stderr=devnull)

    println("Sketching $fasta2 (k=$k, s=$s)...")
    sketch_cmd2 = pipeline(`$mash_path sketch -k $k -s $s -o $sketch2 $fasta2`, stdout=devnull, stderr=devnull)

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
    
#     asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    results5 = read(`minimap2 -x asm5 --cs -cL $reference_fasta $query_fasta`)
    if !isempty(results5)
        results = results5
    else
        @warn "no hit with asm5, trying asm10"
        results10 = read(`minimap2 -x asm10 --cs -cL $reference_fasta $query_fasta`)
        if !isempty(results10)
            results = results10
        else
            @warn "no hits with asm5 or asm10, trying asm20"
            results20 = read(`minimap2 -x asm20 --cs -cL $reference_fasta $query_fasta`)
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run fastani with a query and reference list

Calculate Average Nucleotide Identity (ANI) between genome sequences using FastANI.

# Arguments
- `query_list::String`: Path to file containing list of query genome paths (one per line)
- `reference_list::String`: Path to file containing list of reference genome paths (one per line)
- `outfile::String`: Path to output file that will contain ANI results
- `threads::Int=Sys.CPU_THREADS`: Number of parallel threads to use
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
function fastani_list(;query_list="", reference_list="", outfile="", threads=Sys.CPU_THREADS, force=false)
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