i"""
    mash_distance_from_jaccard(jaccard_index::Float64, kmer_size::Int)

Calculates the Mash distance (an estimate of Average Nucleotide Identity)
from a given Jaccard Index and k-mer size.

# Arguments
- `jaccard_index::Float64`: The Jaccard similarity between the two k-mer sets. Must be between 0.0 and 1.0.
- `kmer_size::Int`: The length of k-mers used to calculate the Jaccard index.

# Returns
- `Float64`: The estimated Mash distance `D`. The estimated ANI would be `1.0 - D`.

# Example
```jldoctest
julia> # Example: Two genomes with a Jaccard index of 0.2, using k=21
julia> mash_distance_from_jaccard(0.2, 21)
0.08053896775831388
```
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

"""
    mash_distance_from_jaccard(jaccard_index::Float64, kmer_size::Int)

Calculates the Mash distance (an estimate of Average Nucleotide Identity)
from a given Jaccard Index and k-mer size.

# Arguments
- `jaccard_index::Float64`: The Jaccard similarity between the two k-mer sets. Must be between 0.0 and 1.0.
- `kmer_size::Int`: The length of k-mers used to calculate the Jaccard index.

# Returns
- `Float64`: The estimated Mash distance `D`. The estimated ANI would be `1.0 - D`.

# Example
```jldoctest
julia> # Example: Two genomes with a Jaccard index of 0.2, using k=21
julia> mash_distance_from_jaccard(0.2, 21)
0.08053896775831388
```
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