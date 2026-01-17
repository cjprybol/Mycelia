# Streaming FASTA/FASTQ normalization to a compressed JSON Lines file, record-by-record,
# with an optional progress meter.

"""
    fastx2normalized_jsonl_stream(; fastx_path::AbstractString,
                                   output_path::Union{Nothing,AbstractString}=nothing,
                                   human_readable_id::Union{String,Nothing}=nothing,
                                   force_truncate::Bool=false,
                                   show_progress::Bool=true,
                                   progress_every::Int=1000,
                                   progress_desc::AbstractString="Streaming normalization",
                                   progress_output::IO=stderr)

Stream the input FASTA/FASTQ file and write one JSON object per line (NDJSON) to
a compressed output file. This avoids loading the entire table or very long
sequences into memory and avoids delimiter limitations of TSV/CSV.

When `output_path` is not provided, a default is derived by replacing the file
extension matched by `Mycelia.FAST(A|Q)_REGEX` (including optional .gz) with
`.normalized.jsonl.gz`.

A lightweight progress meter can be shown during streaming. It displays processed
records, elapsed time, and records/second as an indefinite spinner. Control its
frequency with `progress_every` to reduce overhead.
"""
function fastx2normalized_jsonl_stream(; fastx_path::AbstractString,
                                       output_path::Union{Nothing,AbstractString}=nothing,
                                       human_readable_id::Union{String,Nothing}=nothing,
                                       force_truncate::Bool=false,
                                       show_progress::Bool=true,
                                       progress_every::Int=1000,
                                       progress_desc::AbstractString="Streaming normalization",
                                       progress_output::IO=stderr)

    # --- Extract human_readable_id from filename if not provided ---
    if isnothing(human_readable_id)
        human_readable_id = _extract_human_readable_id_from_fastx_path(fastx_path, force_truncate)
    end

    # --- Input Validation ---
    @assert isfile(fastx_path) && filesize(fastx_path) > 0 "File does not exist or is empty."
    if length(human_readable_id) > 16
        if force_truncate
            @warn "Human-readable identifier '$(human_readable_id)' exceeds 16 characters, truncating to '$(human_readable_id[1:16])'"
            human_readable_id = human_readable_id[1:16]
        else
            error("Human-readable identifier '$(human_readable_id)' cannot exceed 16 characters. Use force_truncate=true to allow truncation.")
        end
    end

    # --- File Type Detection (using Mycelia) ---
    file_type = if occursin(Mycelia.FASTA_REGEX, fastx_path)
        :fasta
    elseif occursin(Mycelia.FASTQ_REGEX, fastx_path)
        :fastq
    else
        error("File is not FASTA or FASTQ")
    end

    # --- Derive output path if not provided ---
    out_path = isnothing(output_path) ? _default_out_path(fastx_path) : String(output_path)

    # --- Stream records and write NDJSON.gz ---
    if endswith(out_path, ".gz")
        open(out_path, "w") do raw
            gz = CodecZlib.GzipCompressorStream(raw)
            try
                _write_ndjson_stream(
                    gz, fastx_path, human_readable_id, file_type;
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_desc=progress_desc,
                    progress_output=progress_output,
                )
                flush(gz)
            finally
                close(gz)
            end
        end
    else
        open(out_path, "w") do io
            _write_ndjson_stream(
                io, fastx_path, human_readable_id, file_type;
                show_progress=show_progress,
                progress_every=progress_every,
                progress_desc=progress_desc,
                progress_output=progress_output,
            )
            flush(io)
        end
    end

    return out_path
end

# --- Helpers ---

# Build the default output path by replacing FASTX suffix (including optional .gz) with .normalized.jsonl.gz
function _default_out_path(fastx_path::AbstractString)
    m_fasta = match(Mycelia.FASTA_REGEX, fastx_path)
    m_fastq = match(Mycelia.FASTQ_REGEX, fastx_path)
    if (m_fasta === nothing) && (m_fastq === nothing)
        error("Input path does not match expected regex pattern: $(fastx_path)")
    else
        if !(m_fasta === nothing)
            output_path = replace(fastx_path, Mycelia.FASTA_REGEX => ".normalized.jsonl.gz")
        else
            output_path = replace(fastx_path, Mycelia.FASTQ_REGEX => ".normalized.jsonl.gz")
        end
    end
    return output_path
end

# Core streaming writer: one JSON object per record, with optional progress meter
function _write_ndjson_stream(io::IO, fastx_path::AbstractString, human_readable_id::String, file_type::Symbol;
                              show_progress::Bool,
                              progress_every::Int,
                              progress_desc::AbstractString,
                              progress_output::IO)

    n_processed::Int = 0
    t0 = time()
    p = nothing
    if show_progress
        p = ProgressMeter.ProgressUnknown(desc=progress_desc, dt=0.5, output=progress_output)
    end

    for record in Mycelia.open_fastx(fastx_path)
        record_sequence = FASTX.sequence(String, record)

        # Quality handling
        record_quality = if file_type == :fastq
            Float64.(collect(FASTX.quality_scores(record)))
        else
            nothing
        end

        # Compute derived fields
        record_identifier = FASTX.identifier(record)
        record_description = FASTX.description(record)
        sequence_hash = create_base58_hash(record_sequence, encoded_length=16)
        record_alphabet = join(sort(collect(Set(uppercase(record_sequence)))))
        record_type = string(Mycelia.detect_alphabet(record_sequence))
        mean_record_quality = (record_quality === nothing) ? nothing : Statistics.mean(record_quality)
        median_record_quality = (record_quality === nothing) ? nothing : Statistics.median(record_quality)
        record_length = length(record_sequence)

        # Simplified hierarchical identifier without joint sequence hash
        sequence_identifier = string(human_readable_id, "_", sequence_hash)
        if length(sequence_identifier) > 50
            error("NCBI identifier length limit of 50 characters exceeded: $(sequence_identifier)")
        end

        # Emit JSON in a deterministic key order via NamedTuple literal
        obj = (
            fastx_path = Base.basename(fastx_path),
            human_readable_id = human_readable_id,
            sequence_hash = sequence_hash,
            sequence_identifier = sequence_identifier,
            record_identifier = record_identifier,
            record_description = record_description,
            record_length = record_length,
            record_alphabet = record_alphabet,
            record_type = record_type,
            mean_record_quality = mean_record_quality,
            median_record_quality = median_record_quality,
            record_quality = record_quality,
            record_sequence = record_sequence,
            file_type = String(file_type),
        )

        JSON.print(io, obj)
        write(io, '\n')

        # Progress
        n_processed += 1
        if show_progress && (n_processed % progress_every == 0)
            elapsed = time() - t0
            rps = elapsed > 0 ? (n_processed / elapsed) : 0.0
            ProgressMeter.next!(p; showvalues=[
                (:records, n_processed),
                (:elapsed, _fmt_hms(elapsed)),
                (:rps, round(rps, digits=2)),
            ])
        end
    end

    if show_progress
        elapsed = time() - t0
        rps = elapsed > 0 ? (n_processed / elapsed) : 0.0
        ProgressMeter.finish!(p; showvalues=[
            (:records, n_processed),
            (:elapsed, _fmt_hms(elapsed)),
            (:rps, round(rps, digits=2)),
        ])
    end

    nothing
end

# Format seconds as HH:MM:SS
function _fmt_hms(t::Real)
    t ≤ 0 && return "00:00:00"
    s = floor(Int, t)
    h = s ÷ 3600
    m = (s % 3600) ÷ 60
    ss = s % 60
    return Printf.@sprintf("%02d:%02d:%02d", h, m, ss)
end

"""
    normalized_table2fastx(
        table::DataFrames.DataFrame;
        output_dir::String=".",
        output_basename::Union{String, Nothing}=nothing,
        force::Bool=false,
        verbose::Bool=false,
        gzip::Bool=false,
        deduplicate::Bool=true
    )

Writes a normalized DataFrame to a FASTA or FASTQ file, with an option for GZIP compression.

# Arguments
- `table::DataFrames.DataFrame`: A DataFrame from `fastx2normalized_table`.

# Keyword Arguments
- `output_dir::String="."`: The directory to save the file.
- `output_basename::Union{String, Nothing}=nothing`: The base name for the output file (without extension). Defaults to the normalized fastx identifier in the table.
- `force::Bool=false`: If `true`, overwrite existing files.
- `verbose::Bool=false`: If `true`, print detailed progress information.
- `gzip::Bool=false`: If `true`, the output file will be GZIP compressed.
- `deduplicate::Bool=true`: If `true`, de-duplicate contigs to keep only one copy of each unique sequence.
"""
function normalized_table2fastx(
    table::DataFrames.DataFrame;
    output_dir::AbstractString=".",
    output_basename::Union{AbstractString, Nothing}=nothing,
    force::Bool=false,
    verbose::Bool=false,
    gzip::Bool=false,
    deduplicate::Bool=true
)
    # --- 1. Determine Output Format and Record Type ---
    is_fastq = !all(ismissing, table.record_quality)
    
    # Update file extension logic to handle compression
    file_extension = is_fastq ? ".fq" : ".fna"
    if gzip
        file_extension *= ".gz"
    end
    
    # --- 2. Determine Final Output Path ---
    local basename
    if isnothing(output_basename)
        if hasproperty(table, :fastx_identifier)
            basename = table.fastx_identifier[1]
        elseif hasproperty(table, :genome_identifier)
            basename = table.genome_identifier[1]
        else
            error("Neither 'fastx_identifier' nor 'genome_identifier' columns found in the table.")
        end
    else
        basename = output_basename
    end
    final_outfile = joinpath(output_dir, basename * file_extension)

    if !isfile(final_outfile) || force
        # --- 3. Construct Vector of Records ---
        RecordType = is_fastq ? FASTX.FASTQ.Record : FASTX.FASTA.Record
        records = RecordType[]
        sizehint!(records, DataFrames.nrow(table))
    
        for row in eachrow(table)
            if is_fastq
                quality_integers = round.(Int, row.record_quality)
                record = FASTX.FASTQ.Record(row.sequence_identifier, row.record_sequence, quality_integers)
            else
                record = FASTX.FASTA.Record(row.sequence_identifier, row.record_sequence)
            end
            push!(records, record)
        end

        if deduplicate
            n_initial_records = length(records)
            unique!(records)
            n_removed = n_initial_records - length(records)
            if verbose && n_removed > 0
                percent_duplication = round(n_removed / n_initial_records * 100, digits=3)
                @warn "$(n_removed) duplicate records detected ($(percent_duplication)%); keeping $(length(records)) unique records"
            end
        end
        
        # --- 4. Call Your Existing Writer Function ---
        if is_fastq
            written_outfile = Mycelia.write_fastq(filename=final_outfile, records=records, gzip=gzip)
        else
            written_outfile = Mycelia.write_fasta(outfile=final_outfile, records=records, gzip=gzip)
        end
        @assert written_outfile == final_outfile
    
        verbose && Printf.@printf "Successfully prepared %d records for writing to %s\n" length(records) final_outfile
    else
        verbose && @warn "$(final_outfile) already present, use force=true to overwrite"
    end
    return final_outfile
end

"""
    _extract_human_readable_id_from_fastx_path(fastx_path::String, force_truncate::Bool=false) -> String

Internal function to intelligently extract a human-readable identifier from a FASTX filename.
Handles common bioinformatics naming patterns and attempts to find the longest meaningful prefix.
"""
function _extract_human_readable_id_from_fastx_path(fastx_path::AbstractString, force_truncate::Bool=false)
    filename = basename(fastx_path)
    
    # Remove file extensions using existing regex patterns
    base_name = if occursin(Mycelia.FASTA_REGEX, filename)
        replace(filename, Mycelia.FASTA_REGEX => "")
    elseif occursin(Mycelia.FASTQ_REGEX, filename) 
        replace(filename, Mycelia.FASTQ_REGEX => "")
    else
        error("File '$(filename)' does not match FASTA or FASTQ naming patterns. Supported extensions: .fasta, .fna, .faa, .fa, .frn, .fastq, .fq (optionally .gz compressed)")
    end
    
    # If already <= 16 characters, return as-is
    if length(base_name) <= 16
        return base_name
    end
    
    # Try to find meaningful prefixes by splitting on common delimiters
    delimiters = ['_', '-', ' ', '.']
    viable_prefixes = String[]
    
    for delimiter in delimiters
        if occursin(delimiter, base_name)
            parts = split(base_name, delimiter)
            # Build cumulative prefixes
            current_prefix = ""
            for part in parts
                test_prefix = isempty(current_prefix) ? part : current_prefix * string(delimiter) * part
                if length(test_prefix) <= 16
                    current_prefix = test_prefix
                    push!(viable_prefixes, current_prefix)
                else
                    break
                end
            end
        end
    end
    
    # Find the longest viable prefix
    if !isempty(viable_prefixes)
        longest_prefix = viable_prefixes[argmax(length.(viable_prefixes))]
        
        # Skip very short prefixes that aren't meaningful for common patterns
        if length(longest_prefix) >= 3
            # Check for common bioinformatics prefixes we want to avoid stopping at
            meaningless_prefixes = ["GCA", "GCF", "NC", "NZ", "NW", "NT", "AC", "AE", "AF", "AY", "DQ", "EF", "EU", "FJ", "GQ", "HM", "JF", "JN", "JQ", "JX", "KC", "KF", "KJ", "KM", "KP", "KR", "KT", "KU", "KX", "KY", "MF", "MG", "MH", "MK", "MN", "MT", "MW", "MZ"]
            
            # If we have a short meaningless prefix, try to get more context
            if longest_prefix in meaningless_prefixes && length(viable_prefixes) > 1
                # Look for a longer prefix that includes more meaningful information
                longer_prefixes = filter(p -> length(p) > length(longest_prefix) && length(p) <= 16, viable_prefixes)
                if !isempty(longer_prefixes)
                    longest_prefix = longer_prefixes[argmax(length.(longer_prefixes))]
                end
            end
            
            return longest_prefix
        end
    end
    
    # If no viable prefix found, either truncate with warning or error
    if force_truncate
        truncated = base_name[1:16]
        @warn "Could not find meaningful identifier prefix for '$(base_name)'. Using truncated version: '$(truncated)'"
        return truncated
    else
        error("Could not extract a viable identifier (≤16 chars) from filename '$(base_name)'. Consider using force_truncate=true or providing human_readable_id explicitly.")
    end
end

"""
    fastx2normalized_table(fastx_path::String; human_readable_id::Union{String,Nothing}=nothing, force_truncate::Bool=false) -> DataFrames.DataFrame
    fastx2normalized_table(; fastx_path::String, human_readable_id::Union{String,Nothing}=nothing, force_truncate::Bool=false) -> DataFrames.DataFrame

$(DocStringExtensions.TYPEDSIGNATURES)

Reads a FASTA or FASTQ file and converts its records into a normalized `DataFrames.DataFrame` with stable, hierarchical identifiers.

# Arguments (positional version)
- `fastx_path::String`: Path to a FASTA or FASTQ file.

# Keyword Arguments
- `fastx_path::String`: Path to a FASTA or FASTQ file.
- `human_readable_id::Union{String,Nothing}`: A human-readable identifier for the genome/entity (max 16 characters). 
  If `nothing`, attempts to extract from filename intelligently.
- `force_truncate::Bool`: If true, truncates long identifiers to 16 characters with warning instead of erroring.

# Returns
- `DataFrames.DataFrame`: A data frame with standardized columns, including the new hierarchical identifiers.
"""
function fastx2normalized_table(fastx_path::AbstractString; human_readable_id::Union{String,Nothing}=nothing, force_truncate::Bool=false)
    return fastx2normalized_table(; fastx_path=fastx_path, human_readable_id=human_readable_id, force_truncate=force_truncate)
end

function fastx2normalized_table(; fastx_path::AbstractString, human_readable_id::Union{String,Nothing}=nothing, force_truncate::Bool=false)
    # --- Extract human_readable_id from filename if not provided ---
    if isnothing(human_readable_id)
        human_readable_id = _extract_human_readable_id_from_fastx_path(fastx_path, force_truncate)
    end

    # --- Input Validation ---
    @assert isfile(fastx_path) && filesize(fastx_path) > 0 "File does not exist or is empty."
    if length(human_readable_id) > 16
        if force_truncate
            @warn "Human-readable identifier '$(human_readable_id)' exceeds 16 characters, truncating to '$(human_readable_id[1:16])'"
            human_readable_id = human_readable_id[1:16]
        else
            error("Human-readable identifier '$(human_readable_id)' cannot exceed 16 characters. Use force_truncate=true to allow truncation.")
        end
    end

    # --- DataFrame Initialization ---
    normalized_table = DataFrames.DataFrame(
        record_identifier = String[],
        record_description = String[],
        sequence_hash = String[],
        record_quality = Union{Vector{Float64}, Missing}[],
        record_alphabet = String[],
        record_type = String[],
        mean_record_quality = Union{Float64, Missing}[],
        median_record_quality = Union{Float64, Missing}[],
        record_length = Int[],
        record_sequence = String[],
    )

    # --- File Type Detection (using Mycelia) ---
    file_type = if occursin(Mycelia.FASTA_REGEX, fastx_path)
        :fasta
    elseif occursin(Mycelia.FASTQ_REGEX, fastx_path)
        :fastq
    else
        error("File is not FASTA or FASTQ")
    end
    
    # --- Record Processing (using Mycelia) ---
    for record in Mycelia.open_fastx(fastx_path)
        record_sequence = FASTX.sequence(String, record)
        record_quality = file_type == :fastq ? collect(FASTX.quality_scores(record)) : missing

        push!(normalized_table, (
            record_identifier = FASTX.identifier(record),
            record_description = FASTX.description(record),
            sequence_hash = create_base58_hash(record_sequence, encoded_length=16),
            record_quality = record_quality,
            record_alphabet = join(sort(collect(Set(uppercase(record_sequence))))),
            record_type = string(Mycelia.detect_alphabet(record_sequence)), # Restored Mycelia call
            mean_record_quality = !ismissing(record_quality) ? Statistics.mean(record_quality) : missing,
            median_record_quality = !ismissing(record_quality) ? Statistics.median(record_quality) : missing,
            record_length = length(record_sequence),
            record_sequence = record_sequence,
        ))
    end

    # --- Add New Hierarchical Identifier Columns ---
    current_columns = names(normalized_table)
    fastx_hash = generate_joint_sequence_hash(normalized_table.sequence_hash, encoded_length=16)

    normalized_table[!, "human_readable_id"] .= human_readable_id
    normalized_table[!, "fastx_hash"] .= fastx_hash
    normalized_table[!, "fastx_identifier"] .= human_readable_id .* "_" .* fastx_hash
    normalized_table[!, "sequence_identifier"] .= normalized_table.fastx_identifier .* "_" .* normalized_table.sequence_hash
    @assert all(length.(normalized_table[!, "sequence_identifier"]) .<= 50) "NCBI identifier length limit of 50 characters exceeded."
    normalized_table[!, "fastx_path"] .= Base.basename(fastx_path)

    final_order = [
        "fastx_path", "human_readable_id", "fastx_hash", "sequence_hash",
        "fastx_identifier", "sequence_identifier", "record_identifier",
        "record_description", "record_length", "record_alphabet", "record_type",
        "mean_record_quality", "median_record_quality", "record_quality",
        "record_sequence"
    ]
    return normalized_table[!, final_order]
end

"""
    create_base58_hash(data_to_hash::AbstractString; encoded_length::Int=64, normalize_case::Bool=true) -> String

Hashes a string using BLAKE3 and returns a Base58 encoded string of a
specified *exact* length.

# Arguments
- `data_to_hash::AbstractString`: Input data to hash
- `encoded_length::Int=64`: Length of the Base58 encoded output (default 64 provides ~380 bits entropy for tree-of-life scale applications)
- `normalize_case::Bool=true`: If true, converts input to uppercase before hashing

# Details
By default, sequences are normalized to uppercase before hashing to ensure
case-insensitive hashing for biological sequences where case differences
are typically not meaningful. The default 64-character length provides
collision resistance suitable for massive multi-omics applications spanning
the sequence diversity of the tree of life.

This function now calls the comprehensive create_blake3_hash function.
"""
function create_base58_hash(data_to_hash::AbstractString; encoded_length::Int=64, normalize_case::Bool=true)::String
    return create_blake3_hash(data_to_hash; encoding=:base58, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=true)
end

"""
    generate_joint_sequence_hash(sequences::Vector{String}; hash_function::Symbol=:blake3, encoding::Symbol=:base58, encoded_length::Int=64, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Computes a single, order-independent hash for a collection of sequences.

This function works by hashing each sequence individually, sorting the resulting
hashes, and then hashing the concatenated list of sorted hashes. This ensures
that the final joint sequence hash is stable, even if the sequences are provided in a
different order.

# Arguments
- `sequences::Vector{String}`: A vector of strings, where each string is a biological sequence.
- `hash_function::Symbol=:blake3`: Hash algorithm to use (:blake3, :sha256, :sha512, :md5, :sha1, :sha3_256, :sha3_512, :crc32)
- `encoding::Symbol=:base58`: Output encoding (:hex, :base58, :base64)
- `encoded_length::Int=64`: The desired length for the final encoded hash (default 64 provides ~380 bits entropy)
- `normalize_case::Bool=true`: If true, converts sequences to uppercase before hashing
- `allow_truncation::Bool=false`: Allow truncation if encoded_length < native length

# Returns
- `String`: A single, stable hash representing the joint set of sequences.
"""
function generate_joint_sequence_hash(sequences::Vector{<:AbstractString}; hash_function::Symbol=:blake3, encoding::Symbol=:base58, encoded_length::Int=64, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    if isempty(sequences)
        error("Input sequence vector cannot be empty.")
    end

    # 1. Hash each sequence individually using create_sequence_hash for consistency
    individual_hashes = Vector{String}(undef, length(sequences))
    for i in 1:length(sequences)
        # Use create_sequence_hash with user-specified hash function and encoding
        individual_hashes[i] = create_sequence_hash(sequences[i], hash_function=hash_function, encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
    end

    # 2. Sort the list of hashes. This is the key to making it order-independent!
    sort!(individual_hashes)

    # 3. Combine the sorted hashes into a single string
    combined_data = join(individual_hashes)

    # 4. Hash the combined string to get the final joint sequence hash with same encoding length
    # Use normalize_case=false here because we're hashing encoded strings, not biological sequences
    final_joint_sequence_hash = create_sequence_hash(combined_data, hash_function=hash_function, encoding=encoding, encoded_length=encoded_length, normalize_case=false, allow_truncation=allow_truncation)

    return final_joint_sequence_hash
end

"""
    create_sequence_hash(data_to_hash::AbstractString; hash_function::Symbol=:blake3, encoding::Symbol=:base58, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Unified biological sequence hashing function with selectable algorithms.

# Arguments
- `data_to_hash::AbstractString`: Input sequence data to hash
- `hash_function::Symbol=:blake3`: Hash algorithm (:blake3, :sha256, :sha512, :md5, :sha1, :sha3_256, :sha3_512, :crc32)
- `encoding::Symbol=:base58`: Output encoding (:hex, :base58, :base64)
- `encoded_length::Union{Int,Missing}=missing`: Desired output length (uses function defaults if missing)
- `normalize_case::Bool=true`: If true, converts input to uppercase before hashing
- `allow_truncation::Bool=false`: Allow truncation if encoded_length < native length

# Returns
- `String`: Hash in specified encoding and length

# Details
This is the recommended function for biological sequence hashing, defaulting to
BLAKE3 with Base58 encoding for optimal collision resistance and compactness.
When encoded_length=missing, each hash function uses its optimal default length.
"""
function create_sequence_hash(data_to_hash::AbstractString; hash_function::Symbol=:blake3, encoding::Symbol=:base58, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    # Validate input - empty sequences are not valid in biological context
    if isempty(data_to_hash)
        error("Empty sequences are not valid for biological sequence hashing")
    end
    if hash_function == :blake3
        # For Blake3, pass through encoded_length (will use 64 as default if missing)
        if ismissing(encoded_length)
            # return create_blake3_hash(data_to_hash; encoding=encoding, encoded_length=64, normalize_case=normalize_case, allow_truncation=allow_truncation)
            return create_blake3_hash(data_to_hash; encoding=encoding, normalize_case=normalize_case, allow_truncation=allow_truncation)
        else
            return create_blake3_hash(data_to_hash; encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
        end
    elseif hash_function == :sha256
        return create_sha256_hash(data_to_hash; encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
    elseif hash_function == :sha512
        return create_sha512_hash(data_to_hash; encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
    elseif hash_function == :md5
        return create_md5_hash(data_to_hash; encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
    elseif hash_function == :sha1
        return create_sha1_hash(data_to_hash; encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
    elseif hash_function == :sha3_256
        return create_sha3_256_hash(data_to_hash; encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
    elseif hash_function == :sha3_512
        return create_sha3_512_hash(data_to_hash; encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
    elseif hash_function == :crc32
        return create_crc32_hash(data_to_hash; encoding=encoding, encoded_length=encoded_length, normalize_case=normalize_case, allow_truncation=allow_truncation)
    else
        error("Unsupported hash function: $hash_function")
    end
end

"""
    create_sequence_hash(seq::BioSequences.LongSequence; kwargs...) -> String

Convenience wrapper for hashing BioSequence objects directly.
Converts the BioSequence to string and passes all keyword arguments to create_sequence_hash.

# Arguments
- `seq::BioSequences.LongSequence`: BioSequence object (DNA, RNA, or amino acid)
- `kwargs...`: All keyword arguments supported by create_sequence_hash

# Examples
```julia
dna_seq = BioSequences.dna"ATCGATCGATCG"
hash = create_sequence_hash(dna_seq, encoded_length=32)

rna_seq = BioSequences.rna"AUCGAUCGAUCG" 
hash = create_sequence_hash(rna_seq, hash_function=:sha256)

aa_seq = BioSequences.aa"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAFGVLLWEIACDESGLAGGWEKGQEGAVDASAEVFRHHVMDSSEQGGGGLLAEGWQQQQQQQQQDQDQSQSSQGWQQQQQQGGWAWQGWQSQQAEGQGWAGSAQAEGQAAEGQAAGWAGSAQAEGQAAGSAQWAGGGGWAGAGGWAGSAQWAGGGGWAGAGMGWAGAGGWAGSAQWAGGGGWAGAGMGWAGAGGWAGSQQQWAQGQQMVAEPGFVGQE"
hash = create_sequence_hash(aa_seq, encoded_length=64)
```
"""
function create_sequence_hash(seq::BioSequences.LongSequence; kwargs...)::String
    return create_sequence_hash(string(seq); kwargs...)
end

"""
    create_sequence_hash(kmer::Kmers.Kmer; kwargs...) -> String

Convenience wrapper for hashing k-mer objects directly.
Converts the k-mer to string and passes all keyword arguments to create_sequence_hash.

# Arguments  
- `kmer::Kmers.Kmer`: k-mer object (DNAKmer, RNAKmer, AAKmer)
- `kwargs...`: All keyword arguments supported by create_sequence_hash

# Examples
```julia
import Kmers, BioSequences
dna_seq = BioSequences.LongDNA{4}("ATCGATCG")
kmer = Kmers.DNAKmer{8}(dna_seq)
hash = create_sequence_hash(kmer, encoded_length=16)

# With different parameters
hash = create_sequence_hash(kmer, hash_function=:md5, encoding=:hex)
```
"""
function create_sequence_hash(kmer::Kmers.Kmer; kwargs...)::String
    return create_sequence_hash(string(kmer); kwargs...)
end

"""
    find_fasta_files(input_path::String) -> Vector{String}

Find all FASTA files in a directory or return single file if path is a file.

Uses the existing `FASTA_REGEX` constant to identify FASTA files.

# Arguments
- `input_path`: Path to directory or single FASTA file

# Returns
- Vector of FASTA file paths

# Example
```julia
fasta_files = find_fasta_files("./genomes/")
```
"""
function find_fasta_files(input_path::String)
    if isfile(input_path)
        if occursin(FASTA_REGEX, input_path)
            return [input_path]
        else
            error("Input file does not match FASTA format: $(input_path)")
        end
    elseif isdir(input_path)
        return filter(f -> occursin(FASTA_REGEX, f), readdir(input_path, join=true))
    else
        error("Input path does not exist: $(input_path)")
    end
end

"""
    join_fastqs_with_uuid(
        fastq_files::Vector{String};
        fastq_out::String
        tsv_out::String
    )

Note: does not keep track of paired-end data - assumes single end reads

Designed primarily to allow joint mapping of many long-read samples

Given a collection of fastq files, creates:
- A gzipped TSV mapping original file and read_id to a new UUID per read
- A gzipped joint fastq file with the new UUID as read header

Returns: Tuple of output file paths (tsv_out, fastq_out)
"""
function join_fastqs_with_uuid(
    fastq_files::Vector{String};
    fastq_out::String = Mycelia.normalized_current_datetime() * ".joint_reads.fq.gz",
    tsv_out::String = replace(fastq_out, Mycelia.FASTQ_REGEX => ".tsv.gz")
)
    # Build mapping as a DataFrame in memory
    mapping = DataFrames.DataFrame(
        input_file = String[],
        original_read_id = String[],
        new_uuid = String[]
    )

    if (!isfile(fastq_out) || filesize(fastq_out) == 0) || (!isfile(tsv_out) || filesize(tsv_out) == 0)
        # Prepare gzipped FASTQ writer using CodecZlib
        gz_out = CodecZlib.GzipCompressorStream(open(fastq_out, "w"))
        writer = FASTX.FASTQ.Writer(gz_out)

        for fastq_file in fastq_files
            for record in Mycelia.open_fastx(fastq_file)
                original_id = String(FASTX.identifier(record))
                uuid = string(UUIDs.uuid4())
                # Save mapping to DataFrame
                DataFrames.push!(mapping, (fastq_file, original_id, uuid))
                # Write record with new UUID as id
                new_record = FASTX.FASTQ.Record(uuid, FASTX.sequence(record), FASTX.quality(record))
                FASTX.write(writer, new_record)
            end
        end

        FASTX.close(writer)
        CodecZlib.close(gz_out)

        # Write mapping table as gzipped TSV using CodecZlib
        tsv_io = CodecZlib.GzipCompressorStream(open(tsv_out, "w"))
        CSV.write(tsv_io, mapping; delim='\t')
        CodecZlib.close(tsv_io)
    else
        @warn "Output files already exist and are not empty: $fastq_out, $tsv_out"
    end
    return (;tsv_out, fastq_out)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Concatenate FASTQ files into a single output FASTQ (optionally gzipped).

# Arguments
- `fastq_files::Vector{String}`: Input FASTQ files (gzipped or plain).
- `output_fastq::String`: Output FASTQ path (gzip inferred from extension if `gzip` is nothing).
- `gzip::Union{Nothing,Bool}=nothing`: Override gzip behavior.
- `force::Bool=false`: Overwrite output if it already exists.

# Returns
- `String`: Path to the concatenated FASTQ.
"""
function concatenate_fastq_files(;
    fastq_files::Vector{String},
    output_fastq::String,
    gzip::Union{Nothing,Bool}=nothing,
    force::Bool=false
)
    @assert !isempty(fastq_files) "fastq_files must not be empty"
    for fastq_file in fastq_files
        @assert isfile(fastq_file) "FASTQ file not found: $(fastq_file)"
    end
    gzip_mode = isnothing(gzip) ? endswith(output_fastq, ".gz") : gzip
    if isfile(output_fastq) && filesize(output_fastq) > 0 && !force
        return output_fastq
    end
    mkpath(dirname(output_fastq))
    out_io = gzip_mode ? CodecZlib.GzipCompressorStream(open(output_fastq, "w")) : open(output_fastq, "w")
    buffer = Vector{UInt8}(undef, 1 << 20)
    for fastq_file in fastq_files
        in_io = endswith(fastq_file, ".gz") ? CodecZlib.GzipDecompressorStream(open(fastq_file, "r")) : open(fastq_file, "r")
        while true
            n = readbytes!(in_io, buffer)
            n == 0 && break
            write(out_io, view(buffer, 1:n))
        end
        close(in_io)
    end
    close(out_io)
    return output_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalize a FASTA/FASTQ identifier for filesystem-friendly use.
"""
function sanitize_fastx_identifier(identifier::AbstractString)
    return replace(String(identifier), r"[^\w\.\-]+" => "_")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Split a FASTA into one file per record, returning a mapping of sanitized id to original id.

# Arguments
- `fasta_in::String`: Input FASTA (optionally gzipped).
- `outdir::String`: Output directory for per-record FASTA files.
- `gzip::Bool=false`: Write gzipped outputs if true.
- `force::Bool=false`: Overwrite outputs if they already exist.

# Returns
- `Dict{String,String}` mapping sanitized identifier to original identifier.
"""
function split_fasta_by_record(;
    fasta_in::String,
    outdir::String,
    gzip::Bool=false,
    force::Bool=false
)
    @assert isfile(fasta_in) "Input FASTA not found: $(fasta_in)"
    mkpath(outdir)
    id_map = Dict{String,String}()
    for record in Mycelia.open_fastx(fasta_in)
        record_id = String(FASTX.identifier(record))
        safe_id = sanitize_fastx_identifier(record_id)
        outfile = joinpath(outdir, safe_id * ".fna" * (gzip ? ".gz" : ""))
        if !isfile(outfile) || force
            Mycelia.write_fasta(outfile=outfile, records=[record], gzip=gzip)
        end
        id_map[safe_id] = record_id
    end
    return id_map
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalize a reference label to its original identifier using a sanitized-id map.
"""
function normalize_reference_label(label::AbstractString, id_map::AbstractDict{<:AbstractString,<:AbstractString})
    base = replace(basename(label), r"\.msh$" => "")
    base = replace(base, Mycelia.FASTA_REGEX => "")
    return get(id_map, base, base)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a FASTA file containing only records whose identifiers are in `ids`.

# Arguments
- `fasta_in::String`: Input FASTA (optionally gzipped).
- `ids::AbstractVector{<:AbstractString}`: Identifiers to retain.
- `fasta_out::String`: Output FASTA path.
- `gzip::Union{Nothing,Bool}=nothing`: Override gzip behavior for the output.
- `force::Bool=false`: Overwrite output if it already exists.
- `require_all::Bool=true`: Error if any requested ids are missing.

# Returns
- `String`: Path to the output FASTA.
"""
function subset_fasta_by_ids(;
    fasta_in::String,
    ids::AbstractVector{<:AbstractString},
    fasta_out::String,
    gzip::Union{Nothing,Bool}=nothing,
    force::Bool=false,
    require_all::Bool=true
)
    @assert isfile(fasta_in) "Input FASTA not found: $(fasta_in)"
    ids_set = Set(String.(ids))
    gzip_mode = isnothing(gzip) ? endswith(fasta_out, ".gz") : gzip
    if isfile(fasta_out) && filesize(fasta_out) > 0 && !force
        return fasta_out
    end
    mkpath(dirname(fasta_out))
    out_io = gzip_mode ? CodecZlib.GzipCompressorStream(open(fasta_out, "w")) : open(fasta_out, "w")
    writer = FASTX.FASTA.Writer(out_io)
    found = Set{String}()
    for record in Mycelia.open_fastx(fasta_in)
        record_id = String(FASTX.identifier(record))
        if record_id in ids_set
            FASTX.write(writer, record)
            push!(found, record_id)
        end
    end
    FASTX.close(writer)
    close(out_io)
    if require_all
        missing = setdiff(ids_set, found)
        if !isempty(missing)
            missing_preview = collect(missing)[1:min(length(missing), 5)]
            error("Missing $(length(missing)) FASTA ids (example: $(missing_preview)).")
        end
    end
    return fasta_out
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Prefix FASTQ read identifiers with a sample-specific tag to guarantee uniqueness.

Streaming-friendly: reads are rewritten on the fly; optional mapping output can be
written as TSV (streamed) or Arrow (buffered). If an output path ends with `.gz`,
the file is written uncompressed first and then gzipped with `pigz`/`gzip` via
`Mycelia.gzip_file` to avoid CodecZlib buffering issues.

# Arguments
- `fastq_file::AbstractString`: Input FASTQ (optionally gzipped).
- `sample_tag::AbstractString`: Tag prepended to each read identifier.
- `out_fastq::AbstractString`: Output FASTQ path (gzipped if endswith `.gz`).

# Keywords
- `mapping_out::Union{Nothing,String}=nothing`: Optional path to write read-id map.
- `mapping_format::Symbol=:arrow`: `:arrow` or `:tsv`. TSV is written streaming;
  Arrow buffers in-memory before write.
- `id_delimiter::AbstractString="::"`: Separator between tag and original id.
- `force::Bool=false`: If `true`, overwrite existing outputs.
- `compress_threads::Integer=get_default_threads()`: Threads for gzip compression.

# Returns
Named tuple `(out_fastq, mapping_out, sample_tag)`.
"""
function prefix_fastq_reads(
    fastq_file::AbstractString;
    sample_tag::AbstractString,
    out_fastq::AbstractString,
    mapping_out::Union{Nothing,String}=nothing,
    mapping_format::Symbol=:arrow,
    id_delimiter::AbstractString="::",
    force::Bool=false,
    compress_threads::Integer=get_default_threads()
)
    @assert isfile(fastq_file) "Input FASTQ not found: $fastq_file"
    @assert mapping_format in (:arrow, :tsv) "mapping_format must be :arrow or :tsv"
    mkpath(dirname(out_fastq))
    sanitized_tag = replace(sample_tag, r"[^\w\.\-]+" => "_")

    nonempty_file(path::AbstractString) = isfile(path) && filesize(path) > 0
    input_bytes = filesize(fastq_file)

    # Always write uncompressed first when `.gz` is requested, then gzip externally.
    out_fastq_plain = endswith(out_fastq, ".gz") ? replace(out_fastq, r"\.gz$" => "") : out_fastq
    tmp_fastq_plain = out_fastq_plain * ".tmp"

    map_plain = nothing
    tmp_map_plain = nothing
    if mapping_out !== nothing
        mkpath(dirname(mapping_out))
        if mapping_format == :arrow
            endswith(mapping_out, ".gz") && error("Arrow mapping_out should not be gzipped: $mapping_out")
        else
            map_plain = endswith(mapping_out, ".gz") ? replace(mapping_out, r"\.gz$" => "") : mapping_out
            tmp_map_plain = map_plain * ".tmp"
        end
    end

    if !force && nonempty_file(out_fastq) && (mapping_out === nothing || nonempty_file(mapping_out))
        return (out_fastq=out_fastq, mapping_out=mapping_out, sample_tag=sanitized_tag)
    end

    # Resume: if the uncompressed intermediate exists from a prior run, just gzip it.
    if endswith(out_fastq, ".gz") && !nonempty_file(out_fastq) && nonempty_file(out_fastq_plain)
        mapping_ready = if mapping_out === nothing
            true
        elseif mapping_format == :arrow
            nonempty_file(mapping_out)
        elseif endswith(mapping_out, ".gz")
            nonempty_file(map_plain) || nonempty_file(mapping_out)
        else
            nonempty_file(mapping_out)
        end

        if mapping_ready
            if mapping_out !== nothing && mapping_format == :tsv && endswith(mapping_out, ".gz") && !nonempty_file(mapping_out) && nonempty_file(map_plain)
                Mycelia.gzip_file(map_plain; outfile=mapping_out, threads=compress_threads, force=force, keep_input=false)
            end
            Mycelia.gzip_file(out_fastq_plain; outfile=out_fastq, threads=compress_threads, force=force, keep_input=false)
            return (out_fastq=out_fastq, mapping_out=mapping_out, sample_tag=sanitized_tag)
        end
    end

    isfile(tmp_fastq_plain) && rm(tmp_fastq_plain; force=true)
    if tmp_map_plain !== nothing
        isfile(tmp_map_plain) && rm(tmp_map_plain; force=true)
    end

    # Gracefully handle 0-byte input files (including 0-byte `.gz` placeholders).
    # For a *valid* gzipped FASTQ with zero reads, the gzip stream still has a header
    # and should not be 0 bytes, so this is safe and avoids a CodecZlib truncation error.
    if input_bytes == 0
        open(tmp_fastq_plain, "w") do _
        end
        if mapping_out !== nothing
            if mapping_format == :tsv
                open(tmp_map_plain, "w") do io
                    write(io, "sample_tag\toriginal_read_id\tnew_read_id\n")
                end
            else
                tmp_arrow = mapping_out * ".tmp"
                isfile(tmp_arrow) && rm(tmp_arrow; force=true)
                Arrow.write(tmp_arrow, DataFrames.DataFrame(sample_tag=String[], original_read_id=String[], new_read_id=String[]))
                mv(tmp_arrow, mapping_out; force=true)
            end
        end

        mv(tmp_fastq_plain, out_fastq_plain; force=true)
        if mapping_out !== nothing && mapping_format == :tsv
            mv(tmp_map_plain, map_plain; force=true)
        end
        if endswith(out_fastq, ".gz")
            Mycelia.gzip_file(out_fastq_plain; outfile=out_fastq, threads=compress_threads, force=force, keep_input=false)
        end
        if mapping_out !== nothing && mapping_format == :tsv && endswith(mapping_out, ".gz")
            Mycelia.gzip_file(map_plain; outfile=mapping_out, threads=compress_threads, force=force, keep_input=false)
        end
        return (out_fastq=out_fastq, mapping_out=mapping_out, sample_tag=sanitized_tag)
    end

    tsv_io = nothing
    mapping_rows = mapping_out === nothing ? nothing : Vector{NamedTuple}()
    if mapping_out !== nothing
        if mapping_format == :tsv
            tsv_io = open(tmp_map_plain, "w")
            write(tsv_io, "sample_tag\toriginal_read_id\tnew_read_id\n")
        else
            mapping_rows = Vector{NamedTuple}()
        end
    end

    out_io = open(tmp_fastq_plain, "w")
    writer = FASTX.FASTQ.Writer(out_io)
    reader = nothing
    try
        reader = Mycelia.open_fastx(fastq_file)
        for record in reader
            original_id = String(FASTX.identifier(record))
            new_id = string(sanitized_tag, id_delimiter, original_id)
            new_record = FASTX.FASTQ.Record(new_id, FASTX.sequence(record), FASTX.quality(record))
            FASTX.write(writer, new_record)
            if mapping_out !== nothing
                if mapping_format == :tsv
                    write(tsv_io, string(sanitized_tag, '\t', original_id, '\t', new_id, '\n'))
                else
                    push!(mapping_rows, (sample_tag=sanitized_tag, original_read_id=original_id, new_read_id=new_id))
                end
            end
        end
    catch err
        msg = sprint(showerror, err)
        # Best-effort cleanup of partial outputs
        try
            reader !== nothing && close(reader)
        catch
        end
        try
            FASTX.close(writer)
        catch
        end
        try
            close(out_io)
        catch
        end
        if mapping_out !== nothing && mapping_format == :tsv
            try
                close(tsv_io)
            catch
            end
        end
        try
            isfile(tmp_fastq_plain) && rm(tmp_fastq_plain; force=true)
        catch
        end
        if tmp_map_plain !== nothing
            try
                isfile(tmp_map_plain) && rm(tmp_map_plain; force=true)
            catch
            end
        end

        if occursin("compressed stream may be truncated", msg) || occursin("ZlibError", msg)
            error(
                "Failed to read gzipped FASTQ (gzip stream appears truncated): $(fastq_file)\n" *
                "File size: $(input_bytes) bytes\n" *
                "This usually indicates an incomplete/corrupt `.gz` (e.g., partially copied, still being written, or 0-byte placeholder).\n" *
                "Suggested checks: `gzip -t $(fastq_file)` (or `pigz -t`) and verify upstream file generation/copy.\n" *
                "Original error: $(msg)"
            )
        end
        rethrow()
    finally
        if reader !== nothing
            try
                close(reader)
            catch
            end
        end
    end
    FASTX.close(writer)
    close(out_io)
    if mapping_out !== nothing && mapping_format == :tsv
        close(tsv_io)
    end

    mv(tmp_fastq_plain, out_fastq_plain; force=true)
    if mapping_out !== nothing && mapping_format == :tsv
        mv(tmp_map_plain, map_plain; force=true)
    end

    if mapping_out !== nothing && mapping_format == :arrow
        tmp_arrow = mapping_out * ".tmp"
        isfile(tmp_arrow) && rm(tmp_arrow; force=true)
        Arrow.write(tmp_arrow, DataFrames.DataFrame(mapping_rows))
        mv(tmp_arrow, mapping_out; force=true)
    end

    if endswith(out_fastq, ".gz")
        Mycelia.gzip_file(out_fastq_plain; outfile=out_fastq, threads=compress_threads, force=force, keep_input=false)
    end
    if mapping_out !== nothing && mapping_format == :tsv && endswith(mapping_out, ".gz")
        Mycelia.gzip_file(map_plain; outfile=mapping_out, threads=compress_threads, force=force, keep_input=false)
    end

    return (out_fastq=out_fastq, mapping_out=mapping_out, sample_tag=sanitized_tag)
end

function uuid_fastq_reads(
    fastq_file::AbstractString;
    out_fastq::AbstractString,
    mapping_out::Union{Nothing,String}=nothing,
    mapping_format::Symbol=:arrow,
    source_fastq::Union{Nothing,String}=nothing,
    force::Bool=false,
    compress_threads::Integer=get_default_threads()
)
    @assert isfile(fastq_file) "Input FASTQ not found: $fastq_file"
    @assert mapping_format in (:arrow, :tsv) "mapping_format must be :arrow or :tsv"
    source_fastq = isnothing(source_fastq) ? fastq_file : source_fastq
    mkpath(dirname(out_fastq))

    if mapping_out !== nothing
        mkpath(dirname(mapping_out))
    end

    nonempty_file(path::AbstractString) = isfile(path) && filesize(path) > 0
    input_bytes = filesize(fastq_file)

    out_fastq_plain = endswith(out_fastq, ".gz") ? replace(out_fastq, r"\.gz$" => "") : out_fastq
    tmp_fastq_plain = out_fastq_plain * ".tmp"

    map_plain = nothing
    tmp_map_plain = nothing
    if mapping_out !== nothing
        if mapping_format == :arrow
            endswith(mapping_out, ".gz") && error("Arrow mapping_out should not be gzipped: $mapping_out")
        else
            map_plain = endswith(mapping_out, ".gz") ? replace(mapping_out, r"\.gz$" => "") : mapping_out
            tmp_map_plain = map_plain * ".tmp"
        end
    end

    if !force && nonempty_file(out_fastq) && (mapping_out === nothing || nonempty_file(mapping_out))
        return (out_fastq=out_fastq, mapping_out=mapping_out)
    end

    if endswith(out_fastq, ".gz") && !nonempty_file(out_fastq) && nonempty_file(out_fastq_plain)
        mapping_ready = if mapping_out === nothing
            true
        elseif mapping_format == :arrow
            nonempty_file(mapping_out)
        elseif endswith(mapping_out, ".gz")
            nonempty_file(map_plain) || nonempty_file(mapping_out)
        else
            nonempty_file(mapping_out)
        end

        if mapping_ready
            if mapping_out !== nothing && mapping_format == :tsv && endswith(mapping_out, ".gz") && !nonempty_file(mapping_out) && nonempty_file(map_plain)
                Mycelia.gzip_file(map_plain; outfile=mapping_out, threads=compress_threads, force=force, keep_input=false)
            end
            Mycelia.gzip_file(out_fastq_plain; outfile=out_fastq, threads=compress_threads, force=force, keep_input=false)
            return (out_fastq=out_fastq, mapping_out=mapping_out)
        end
    end

    isfile(tmp_fastq_plain) && rm(tmp_fastq_plain; force=true)
    if tmp_map_plain !== nothing
        isfile(tmp_map_plain) && rm(tmp_map_plain; force=true)
    end

    if input_bytes == 0
        open(tmp_fastq_plain, "w") do _
        end
        if mapping_out !== nothing
            if mapping_format == :tsv
                open(tmp_map_plain, "w") do io
                    write(io, "source_fastq\toriginal_read_id\tuuid\n")
                end
            else
                tmp_arrow = mapping_out * ".tmp"
                isfile(tmp_arrow) && rm(tmp_arrow; force=true)
                Arrow.write(tmp_arrow, DataFrames.DataFrame(source_fastq=String[], original_read_id=String[], uuid=String[]))
                mv(tmp_arrow, mapping_out; force=true)
            end
        end

        mv(tmp_fastq_plain, out_fastq_plain; force=true)
        if mapping_out !== nothing && mapping_format == :tsv
            mv(tmp_map_plain, map_plain; force=true)
        end
        if endswith(out_fastq, ".gz")
            Mycelia.gzip_file(out_fastq_plain; outfile=out_fastq, threads=compress_threads, force=force, keep_input=false)
        end
        if mapping_out !== nothing && mapping_format == :tsv && endswith(mapping_out, ".gz")
            Mycelia.gzip_file(map_plain; outfile=mapping_out, threads=compress_threads, force=force, keep_input=false)
        end
        return (out_fastq=out_fastq, mapping_out=mapping_out)
    end

    tsv_io = nothing
    mapping_rows = mapping_out === nothing ? nothing : Vector{NamedTuple}()
    if mapping_out !== nothing
        if mapping_format == :tsv
            tsv_io = open(tmp_map_plain, "w")
            write(tsv_io, "source_fastq\toriginal_read_id\tuuid\n")
        else
            mapping_rows = Vector{NamedTuple}()
        end
    end

    out_io = open(tmp_fastq_plain, "w")
    writer = FASTX.FASTQ.Writer(out_io)
    reader = nothing
    try
        reader = Mycelia.open_fastx(fastq_file)
        for record in reader
            original_id = String(FASTX.identifier(record))
            new_id = string(UUIDs.uuid4())
            new_record = FASTX.FASTQ.Record(new_id, FASTX.sequence(record), FASTX.quality(record))
            FASTX.write(writer, new_record)
            if mapping_out !== nothing
                if mapping_format == :tsv
                    write(tsv_io, string(source_fastq, '\t', original_id, '\t', new_id, '\n'))
                else
                    push!(mapping_rows, (source_fastq=source_fastq, original_read_id=original_id, uuid=new_id))
                end
            end
        end
    catch err
        msg = sprint(showerror, err)
        try
            reader !== nothing && close(reader)
        catch
        end
        try
            FASTX.close(writer)
        catch
        end
        try
            close(out_io)
        catch
        end
        if mapping_out !== nothing && mapping_format == :tsv
            try
                close(tsv_io)
            catch
            end
        end
        if occursin("compressed stream may be truncated", msg) || occursin("ZlibError", msg)
            error(
                "Failed to read gzipped FASTQ (gzip stream appears truncated): $(fastq_file)\n" *
                "This usually indicates an incomplete/corrupt `.gz` (e.g., partially copied, still being written, or 0-byte placeholder).\n" *
                "Original error: $(msg)"
            )
        end
        rethrow()
    finally
        try
            reader !== nothing && close(reader)
        catch
        end
        try
            FASTX.close(writer)
        catch
        end
        try
            close(out_io)
        catch
        end
        if mapping_out !== nothing && mapping_format == :tsv
            try
                close(tsv_io)
            catch
            end
        end
    end

    mv(tmp_fastq_plain, out_fastq_plain; force=true)
    if mapping_out !== nothing
        if mapping_format == :tsv
            mv(tmp_map_plain, map_plain; force=true)
        else
            tmp_arrow = mapping_out * ".tmp"
            isfile(tmp_arrow) && rm(tmp_arrow; force=true)
            Arrow.write(tmp_arrow, DataFrames.DataFrame(mapping_rows))
            mv(tmp_arrow, mapping_out; force=true)
        end
    end

    if endswith(out_fastq, ".gz")
        Mycelia.gzip_file(out_fastq_plain; outfile=out_fastq, threads=compress_threads, force=force, keep_input=false)
    end
    if mapping_out !== nothing && mapping_format == :tsv && endswith(mapping_out, ".gz")
        Mycelia.gzip_file(map_plain; outfile=mapping_out, threads=compress_threads, force=force, keep_input=false)
    end

    return (out_fastq=out_fastq, mapping_out=mapping_out)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Concatenate FASTQ/FASTA files (gz or plain) into a single output file.
Maintains input order and works for both single-end and paired-end files (one stream).
Uses system `cat`/`zcat` (or `gzip -cd`) and optional `pigz`/`gzip` compression when
the output path ends with `.gz`.
"""
function concatenate_fastx(inputs::AbstractVector{<:AbstractString}; output_path::AbstractString, threads=get_default_threads(), force=false)
    if isempty(inputs)
        error("No input files provided for concatenation")
    end
    for f in inputs
        if !isfile(f)
            error("Input file not found: $(f)")
        end
    end
    mkpath(dirname(output_path))
    out_plain = endswith(output_path, ".gz") ? replace(output_path, r"\.gz$" => "") : output_path
    tmp_out = out_plain * ".tmp"
    cat_cmd = Sys.which("cat")
    zcat_cmd = Sys.which("zcat")
    gzip_cmd = Sys.which("gzip")
    has_gz_inputs = any(endswith.(inputs, ".gz"))
    use_system = cat_cmd !== nothing && (!has_gz_inputs || zcat_cmd !== nothing || gzip_cmd !== nothing)
    if use_system
        open(tmp_out, "w") do out_io
            for f in inputs
                cmd = if endswith(f, ".gz")
                    if zcat_cmd !== nothing
                        Cmd([zcat_cmd, f])
                    else
                        Cmd([gzip_cmd, "-cd", f])
                    end
                else
                    Cmd([cat_cmd, f])
                end
                run(pipeline(cmd, stdout=out_io))
            end
        end
    else
        open(tmp_out, "w") do out_io
            buffer = Vector{UInt8}(undef, 1 << 20)
            for f in inputs
                in_io = endswith(f, ".gz") ? CodecZlib.GzipDecompressorStream(open(f, "r")) : open(f, "r")
                while true
                    n = readbytes!(in_io, buffer)
                    n == 0 && break
                    write(out_io, view(buffer, 1:n))
                end
                close(in_io)
            end
        end
    end
    mv(tmp_out, out_plain; force=true)
    if endswith(output_path, ".gz")
        Mycelia.gzip_file(out_plain; outfile=output_path, threads=threads, force=force, keep_input=false)
    end
    return output_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run FastQC on FASTQ files. Supports both single-end and paired-end reads.

# Arguments
- `fastq::String`: Path to single-end FASTQ file (optional if forward/reverse provided)
- `forward::String`: Path to forward reads FASTQ file (optional if fastq provided)
- `reverse::String`: Path to reverse reads FASTQ file (optional if fastq provided)
- `outdir::String`: Output directory (auto-generated if not provided)

# Examples
```julia
# Single-end reads
run_fastqc(fastq="reads.fastq")

# Paired-end reads
run_fastqc(forward="reads_R1.fastq", reverse="reads_R2.fastq")
```
"""
function run_fastqc(; fastq::Union{String,Nothing}=nothing, 
                      forward::Union{String,Nothing}=nothing, 
                      reverse::Union{String,Nothing}=nothing,
                      outdir::Union{String,Nothing}=nothing)
    
    Mycelia.add_bioconda_env("fastqc")
    
    # Determine input files and output directory
    if fastq !== nothing
        # Single-end mode
        if forward !== nothing || reverse !== nothing
            error("Cannot specify both 'fastq' and 'forward'/'reverse' arguments")
        end
        input_files = [fastq]
        if outdir === nothing
            outdir = replace(fastq, Mycelia.FASTQ_REGEX => "_fastqc")
        end
    elseif forward !== nothing && reverse !== nothing
        # Paired-end mode
        input_files = [forward, reverse]
        if outdir === nothing
            outdir = Mycelia.find_matching_prefix(forward, reverse) * "_fastqc"
        end
    else
        error("Must specify either 'fastq' for single-end or both 'forward' and 'reverse' for paired-end")
    end
    
    # Create output directory and run FastQC
    if !isdir(outdir)
        mkpath(outdir)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n fastqc fastqc --outdir $(outdir) $(input_files...)`)
    else
        @warn "$outdir already exists"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyze sequence duplication rates in a FASTQ file.

This function processes a FASTQ file to quantify both exact sequence duplications and 
canonical duplications (considering sequences and their reverse complements as equivalent).
The function makes two passes through the file: first to count total records, then to
analyze unique sequences.

# Arguments
- `fastq::String`: Path to the input FASTQ file to analyze
- `results_table::String`: Optional. Path where the results will be saved as a tab-separated file.
  Defaults to the same path as the input file but with extension changed to ".duplication_rates.tsv"

# Returns
- `String`: Path to the results table file

# Output
Generates a tab-separated file containing the following metrics:
- `total_records`: Total number of sequence records in the file
- `total_unique_observations`: Count of unique sequence strings
- `total_unique_canonical_observations`: Count of unique canonical sequences 
  (after normalizing for reverse complements)
- `percent_unique_observations`: Percentage of sequences that are unique
- `percent_unique_canonical_observations`: Percentage of sequences that are unique after canonicalization
- `percent_duplication_rate`: Percentage of sequences that are duplicates (100 - percent_unique_observations)
- `percent_canonical_duplication_rate`: Percentage of sequences that are duplicates after canonicalization

# Notes
- If the specified results file already exists and is not empty, the function will
  return early without recomputing.
- Progress is displayed during processing with a progress bar showing speed.

# Example
```julia
# Analyze a FASTQ file and save results to default location
result_path = assess_duplication_rates("data/sample.fastq")

# Specify custom output path
result_path = assess_duplication_rates("data/sample.fastq", results_table="results/duplication_analysis.tsv")
```
"""
function assess_duplication_rates(fastq; results_table=replace(fastq, Mycelia.FASTQ_REGEX => ".duplication_rates.tsv"))
    # @show results_table
    if isfile(results_table) && (filesize(results_table) > 0)
        println("$results_table already exists.")
        display(DataFrames.DataFrame(uCSV.read(results_table, delim='\t', header=1)))
        return results_table
    end
    # First pass: count total records
    println("Counting total records in FASTQ file...")
    total_records = 0
    for _ in Mycelia.open_fastx(fastq)
        total_records += 1
    end
    println("Found $total_records total records.")
    
    # Initialize sets for unique sequences
    unique_observations = Set{String}()
    unique_canonical_observations = Set{String}()
    
    # Create progress meter
    println("Analyzing duplication rates...")
    prog = ProgressMeter.Progress(total_records, desc="Processing: ", showspeed=true)
    
    # Second pass: process records with progress meter
    for record in Mycelia.open_fastx(fastq)
        seq = FASTX.sequence(record)
        # converted_seq = Mycelia.convert_sequence(seq)
        converted_seq = BioSequences.LongDNA{4}(seq)
        push!(unique_observations, string(converted_seq))
        push!(unique_canonical_observations, string(BioSequences.canonical(converted_seq)))
        
        # Update progress meter
        ProgressMeter.next!(prog)
    end
    
    total_unique_observations = length(unique_observations)
    total_unique_canonical_observations = length(unique_canonical_observations)
    percent_unique_observations = (total_unique_observations / total_records) * 100
    percent_unique_canonical_observations = (total_unique_canonical_observations / total_records) * 100
    percent_duplication_rate = 100 - percent_unique_observations
    percent_canonical_duplication_rate = 100 - percent_unique_canonical_observations

    results_df = DataFrames.DataFrame(;total_records,
                total_unique_observations,
                total_unique_canonical_observations,
                percent_unique_observations,
                percent_unique_canonical_observations,
                percent_duplication_rate,
                percent_canonical_duplication_rate
            )
    display(results_df)
    uCSV.write(results_table, results_df, delim='\t')
    return results_table
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
Requires trim_galore conda environment
"""
function trim_galore_paired(;forward_reads::String, reverse_reads::String, outdir::String=pwd(), threads=min(get_default_threads(), 4))
    Mycelia.add_bioconda_env("trim-galore")
    # Create output directory if it doesn't exist
    trim_galore_dir = mkpath(joinpath(outdir, "trim_galore"))
    
    # Get base filename without path and extension for output naming
    forward_base = basename(forward_reads)
    reverse_base = basename(reverse_reads)
    
    # Construct output filenames according to trim_galore naming convention
    trimmed_forward = joinpath(trim_galore_dir, replace(forward_base, Mycelia.FASTQ_REGEX => "_val_1.fq.gz"))
    trimmed_reverse = joinpath(trim_galore_dir, replace(reverse_base, Mycelia.FASTQ_REGEX => "_val_2.fq.gz"))
    
    if !isfile(trimmed_forward) && !isfile(trimmed_reverse)
        cmd = `$(Mycelia.CONDA_RUNNER) run -n trim-galore trim_galore --suppress_warn --cores $(threads) --output_dir $(trim_galore_dir) --paired $(forward_reads) $(reverse_reads)`
        run(cmd)
    else
        @info "$(trimmed_forward) & $(trimmed_reverse) already present"
    end
    
    return (;trimmed_forward, trimmed_reverse, outdir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform quality control (QC) filtering and trimming on paired-end short-read FASTQ files using fastp.

# Arguments
- `forward_reads::String`: Path to the forward (R1) FASTQ file.
- `reverse_reads::String`: Path to the reverse (R2) FASTQ file.
- `out_forward::String`: Output path for filtered forward reads (auto-generated if not specified).
- `out_reverse::String`: Output path for filtered reverse reads (auto-generated if not specified).
- `report_title::String`: Title for the HTML/JSON report (auto-generated if not specified).
- `html::String`: Output path for HTML report (auto-generated if not specified).
- `json::String`: Output path for JSON report (auto-generated if not specified).
- `enable_dedup::Union{Bool,Nothing}`: Control deduplication behavior. If `true`, forces deduplication with memory-aware settings. If `false`, disables deduplication. If `nothing` (default), uses automatic logic based on file size and available memory.

# Returns
- Named tuple containing paths to: `(out_forward, out_reverse, json, html)`

# Details
This function uses fastp to perform quality control, adapter trimming, and optional deduplication on paired-end reads.

## Smart Deduplication Logic
The function includes intelligent memory management for deduplication:

- **User Control**: Set `enable_dedup=true` to force deduplication, or `enable_dedup=false` to disable it.
- **Automatic Mode**: When `enable_dedup=nothing` (default), the function:
  - Skips deduplication for small files (< 100MB total) for efficiency
  - For larger files, enables deduplication if sufficient memory is available
- **Memory-Aware Settings**: Automatically adjusts fastp's `--dup_calc_accuracy` based on available system memory:
  - Default: Level 3 (4GB memory) if sufficient memory available
  - Fallback: Level 2 (2GB memory) or Level 1 (1GB memory) if needed
  - Disables deduplication if < 1GB memory available

This ensures the process won't be killed due to out-of-memory conditions while maintaining deduplication benefits when feasible.
"""
function qc_filter_short_reads_fastp(;
        forward_reads::String,
        reverse_reads::String,
        out_forward::String=replace(forward_reads, Mycelia.FASTQ_REGEX => ".fastp.1.fq.gz"),
        out_reverse::String=replace(reverse_reads, Mycelia.FASTQ_REGEX => ".fastp.2.fq.gz"),
        report_title::String="$(forward_reads) $(reverse_reads) fastp report",
        html::String=Mycelia.find_matching_prefix(out_forward, out_reverse) * ".fastp_report.html",
        json::String=Mycelia.find_matching_prefix(out_forward, out_reverse) * ".fastp_report.json",
        enable_dedup::Union{Bool,Nothing}=nothing
    )
    # usage: fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]
    # options:
    #   # I/O options
    #   -i, --in1                          read1 input file name (string)
    #   -o, --out1                         read1 output file name (string [=])
    #   -I, --in2                          read2 input file name (string [=])
    #   -O, --out2                           read2 output file name (string [=])
    #       --unpaired1                      for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it. (string [=])
    #       --unpaired2                      for PE input, if read2 passed QC but read1 not, it will be written to unpaired2. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])
    #       --failed_out                     specify the file to store reads that cannot pass the filters. (string [=])
    #       --overlapped_out                 for each read pair, output the overlapped region if it has no any mismatched base. (string [=])
    #   -m, --merge                          for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.
    #       --merged_out                     in the merging mode, specify the file name to store merged output, or specify --stdout to stream the merged output (string [=])
    #       --include_unmerged               in the merging mode, write the unmerged or unpaired reads to the file specified by --merge. Disabled by default.
    #   -6, --phred64                      indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
    #   -z, --compression                  compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
    #       --stdin                          input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
    #       --stdout                         output passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end input. Disabled by default.
    #       --interleaved_in                 indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.
    #       --reads_to_process             specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
    #       --dont_overwrite               don't overwrite existing files. Overwritting is allowed by default.
    #       --fix_mgi_id                     the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.
    
    #   # adapter trimming options
    #   -A, --disable_adapter_trimming     adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
    #   -a, --adapter_sequence               the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
    #       --adapter_sequence_r2            the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=])
    #       --adapter_fasta                  specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
    #       --detect_adapter_for_pe          by default, the adapter sequence auto-detection is enabled for SE data only, turn on this option to enable it for PE data.
    
    #   # global trimming options
    #   -f, --trim_front1                    trimming how many bases in front for read1, default is 0 (int [=0])
    #   -t, --trim_tail1                     trimming how many bases in tail for read1, default is 0 (int [=0])
    #   -b, --max_len1                       if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation (int [=0])
    #   -F, --trim_front2                    trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
    #   -T, --trim_tail2                     trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
    #   -B, --max_len2                       if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])
    
    #   # duplication evaluation and deduplication
    #   -D, --dedup                          enable deduplication to drop the duplicated reads/pairs
    #       --dup_calc_accuracy              accuracy level to calculate duplication (1~6), higher level uses more memory (1G, 2G, 4G, 8G, 16G, 24G). Default 1 for no-dedup mode, and 3 for dedup mode. (int [=0])
    #       --dont_eval_duplication          don't evaluate duplication rate to save time and use less memory.
    
    #   # polyG tail trimming, useful for NextSeq/NovaSeq data
    #   -g, --trim_poly_g                  force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
    #       --poly_g_min_len                 the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
    #   -G, --disable_trim_poly_g          disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
    
    #   # polyX tail trimming
    #   -x, --trim_poly_x                    enable polyX trimming in 3' ends.
    #       --poly_x_min_len                 the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
    
    #   # per read cutting by quality options
    #   -5, --cut_front                      move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
    #   -3, --cut_tail                       move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
    #   -r, --cut_right                      move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
    #   -W, --cut_window_size                the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
    #   -M, --cut_mean_quality               the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])
    #       --cut_front_window_size          the window size option of cut_front, default to cut_window_size if not specified (int [=4])
    #       --cut_front_mean_quality         the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
    #       --cut_tail_window_size           the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
    #       --cut_tail_mean_quality          the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
    #       --cut_right_window_size          the window size option of cut_right, default to cut_window_size if not specified (int [=4])
    #       --cut_right_mean_quality         the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])
    
    #   # quality filtering options
    #   -Q, --disable_quality_filtering    quality filtering is enabled by default. If this option is specified, quality filtering is disabled
    #   -q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
    #   -u, --unqualified_percent_limit    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
    #   -n, --n_base_limit                 if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
    #   -e, --average_qual                 if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])
    
    
    #   # length filtering options
    #   -L, --disable_length_filtering     length filtering is enabled by default. If this option is specified, length filtering is disabled
    #   -l, --length_required              reads shorter than length_required will be discarded, default is 15. (int [=15])
    #       --length_limit                 reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
    
    #   # low complexity filtering
    #   -y, --low_complexity_filter          enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
    #   -Y, --complexity_threshold           the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])
    
    #   # filter reads with unwanted indexes (to remove possible contamination)
    #       --filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
    #       --filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
    #       --filter_by_index_threshold      the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])
    
    #   # base correction by overlap analysis options
    #   -c, --correction                   enable base correction in overlapped regions (only for PE data), default is disabled
    #       --overlap_len_require            the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
    #       --overlap_diff_limit             the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
    #       --overlap_diff_percent_limit     the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])
    
    #   # UMI processing
    #   -U, --umi                          enable unique molecular identifier (UMI) preprocessing
    #       --umi_loc                      specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
    #       --umi_len                      if the UMI is in read1/read2, its length should be provided (int [=0])
    #       --umi_prefix                   if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
    #       --umi_skip                       if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])
    
    #   # overrepresented sequence analysis
    #   -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
    #   -P, --overrepresentation_sampling    One in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])
    
    #   # reporting options
    #   -j, --json                         the json format report file name (string [=fastp.json])
    #   -h, --html                         the html format report file name (string [=fastp.html])
    #   -R, --report_title                 should be quoted with ' or ", default is "fastp report" (string [=fastp report])
    
    #   # threading options
    #   -w, --thread                       worker thread number, default is 3 (int [=3])
    
    #   # output splitting options
    #   -s, --split                        split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
    #   -S, --split_by_lines               split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
    #   -d, --split_prefix_digits          the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])
    
    #   # help
    #   -?, --help                         print this message
    if !isfile(out_forward) || !isfile(out_reverse) || !isfile(json) || !isfile(html)
        Mycelia.add_bioconda_env("fastp")
        
        ## Smart deduplication logic based on file size and available memory
        dedup_args = String[]
        
        if enable_dedup === false
            ## User explicitly disabled dedup - don't add any dedup arguments
            @info "Deduplication explicitly disabled by user"
        elseif enable_dedup === true
            ## User explicitly enabled dedup - use default settings but check memory
            file_size_bytes = filesize(forward_reads) + filesize(reverse_reads)
            available_memory = Sys.free_memory()
            
            ## fastp dedup accuracy levels: 1G, 2G, 4G, 8G, 16G, 24G (levels 1-6)
            ## Default is level 3 (4GB) for dedup mode
            dedup_memory_gb = 4
            dedup_memory_bytes = dedup_memory_gb * 1024^3
            
            if dedup_memory_bytes > available_memory
                @warn "Not enough memory for default deduplication (need $(Base.format_bytes(dedup_memory_bytes)), have $(Base.format_bytes(available_memory))). Using lower accuracy level."
                ## Try progressively lower accuracy levels
                if available_memory >= 2 * 1024^3  ## 2GB available
                    push!(dedup_args, "--dedup", "--dup_calc_accuracy", "2")
                    @info "Using dedup accuracy level 2 (2GB memory)"
                elseif available_memory >= 1024^3  ## 1GB available
                    push!(dedup_args, "--dedup", "--dup_calc_accuracy", "1")
                    @info "Using dedup accuracy level 1 (1GB memory)"
                else
                    @warn "Insufficient memory for deduplication (need at least 1GB). Disabling deduplication."
                end
            else
                push!(dedup_args, "--dedup")
                @info "Using default deduplication (accuracy level 3, 4GB memory)"
            end
        else
            ## enable_dedup is nothing - use automatic logic based on file size and memory
            file_size_bytes = filesize(forward_reads) + filesize(reverse_reads)
            available_memory = Sys.free_memory()
            
            ## Heuristic: if files are small (< 100MB total), dedup is probably not worth the memory cost
            if file_size_bytes < 100 * 1024^2
                @info "Small input files ($(Base.format_bytes(file_size_bytes))). Skipping deduplication for efficiency."
            else
                ## For larger files, check if we have enough memory for dedup
                dedup_memory_bytes = 4 * 1024^3  ## Default level 3 needs 4GB
                
                if dedup_memory_bytes <= available_memory
                    push!(dedup_args, "--dedup")
                    @info "Enabling deduplication with default settings (4GB memory, $(Base.format_bytes(available_memory)) available)"
                elseif available_memory >= 2 * 1024^3
                    push!(dedup_args, "--dedup", "--dup_calc_accuracy", "2")
                    @info "Enabling deduplication with reduced accuracy (2GB memory, $(Base.format_bytes(available_memory)) available)"
                elseif available_memory >= 1024^3
                    push!(dedup_args, "--dedup", "--dup_calc_accuracy", "1")
                    @info "Enabling deduplication with minimal accuracy (1GB memory, $(Base.format_bytes(available_memory)) available)"
                else
                    @warn "Insufficient memory for deduplication (need at least 1GB, have $(Base.format_bytes(available_memory))). Disabling deduplication."
                end
            end
        end
        
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastp fastp
                --in1 $(forward_reads)
                --in2 $(reverse_reads)
                --out1 $(out_forward)
                --out2 $(out_reverse)
                --json $(json)
                --html $(html)
                --report_title $(report_title)
                $(dedup_args)`
        run(cmd)
    else
        @show isfile(out_forward)
        @show isfile(out_reverse)
        @show isfile(json)
        @show isfile(html)
    end
    return (
        out_forward = String(out_forward),
        out_reverse = String(out_reverse),
        json = String(json),
        html = String(html)
    )
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
- Named tuple containing:
  - `out_fastq`: Filtered FASTQ file
  - `html_report`: fastplong HTML report
  - `json_report`: fastplong JSON report

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
    outputs_exist = isfile(out_fastq) && isfile(html_report) && isfile(json_report)
    if !outputs_exist
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
    end
    return (
        out_fastq = String(out_fastq),
        html_report = String(html_report),
        json_report = String(json_report)
    )
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
        out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".filtlong.fq.gz"),
        min_mean_q = 20,
        keep_percent = 95
    )
    if !isfile(out_fastq)
        Mycelia.add_bioconda_env("filtlong")
        Mycelia.add_bioconda_env("pigz")
        p1 = pipeline(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n filtlong filtlong --min_mean_q $(min_mean_q) --keep_percent $(keep_percent) $(in_fastq)`,
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz`
        )
        p2 = pipeline(p1, out_fastq)
        run(p2)
    end
    return out_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform QC filtering on long-read FASTQ files using chopper.

# Arguments
- `in_fastq::String`: Path to the input FASTQ file.
- `out_fastq::String`: Path to the output FASTQ file (optional, auto-generated if not provided).
- `quality_threshold::Int`: Minimum average quality to retain a read (default 20).
- `min_length::Int`: Minimum read length (default 1000).

# Returns
- `String`: Path to the filtered FASTQ file.

# Details
This function uses chopper to discard long reads that do not meet the minimum quality or length thresholds.
It is intended for Oxford Nanopore or similar long-read datasets.

# Dependencies
Requires chopper to be installed via conda
"""
function qc_filter_long_reads_chopper(;
        in_fastq::String,
        out_fastq::String = replace(in_fastq, Mycelia.FASTQ_REGEX => ".chopper.fq.gz"),
        quality_threshold::Int = 20,
        min_length::Int = 1000
    )
    if !isfile(out_fastq)
        Mycelia.add_bioconda_env("chopper")
        # chopper reads from STDIN and writes to STDOUT, so we pipe the file
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n chopper chopper --quality $(quality_threshold) --minlength $(min_length)`
        
        # Handle both compressed and uncompressed input files
        input_cmd = if endswith(in_fastq, ".gz")
            `zcat $(in_fastq)`
        else
            `cat $(in_fastq)`
        end
        
        # Create pipeline: input file -> chopper -> gzip -> output
        pipeline_cmd = pipeline(
            input_cmd,
            cmd,
            `gzip`
        )
        
        run(pipeline(pipeline_cmd, out_fastq))
    end
    return out_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate basic statistics for FASTQ/FASTA sequence files using seqkit.

# Arguments
- `fastq::String`: Path to input FASTQ/FASTA file

# Details
Automatically installs and uses seqkit from Bioconda to compute sequence statistics
including number of sequences, total bases, GC content, average length, etc.

# Dependencies
- Requires Conda and Bioconda channel
- Installs seqkit package if not present

# Returns
Returns a DataFrame of the table

https://bioinf.shenwei.me/seqkit/usage/#stats
"""
function fastx_stats(fastx)
    Mycelia.add_bioconda_env("seqkit")
    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit stats --N 90 --all --tabular $(fastx)`
    return DataFrames.DataFrame(uCSV.read(open(cmd), header=1, delim='\t'))
end

"""
    write_fastq(; records, filename=nothing, outfile=nothing, gzip=false)

$(DocStringExtensions.TYPEDSIGNATURES)

Write FASTQ records to file using FASTX.jl.
Validates extension: .fastq, .fq, .fastq.gz, or .fq.gz.
If `gzip` is true or filename endswith .gz, output is gzipped.
`records` must be an iterable of FASTX.FASTQ.Record.
Accepts `outfile` as an alias for `filename` for consistency with write_fasta.
"""
function write_fastq(;records, filename=nothing, outfile=nothing, gzip=false)
    output_path = isnothing(filename) ? outfile : filename
    if isnothing(output_path)
        error("Provide either filename or outfile")
    end
    if !isnothing(filename) && !isnothing(outfile) && filename != outfile
        error("Provide only one of filename or outfile (got $filename and $outfile)")
    end

    function is_valid_fastq_ext(fname)
        any(endswith.(fname, [".fastq", ".fq", ".fastq.gz", ".fq.gz"]))
    end

    if !is_valid_fastq_ext(output_path)
        error("File extension must be .fastq, .fq, .fastq.gz, or .fq.gz")
    end

    gzip_out = gzip || endswith(output_path, ".gz")
    if gzip_out
        io = CodecZlib.GzipCompressorStream(open(output_path, "w"))
    else
        io = open(output_path, "w")
    end

    try
        writer = FASTX.FASTQ.Writer(io)
        for record in records
            FASTX.write(writer, record)
        end
        FASTX.close(writer)
    finally
        close(io)
    end
    return output_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return an alphabet hint from a FASTA/FASTQ filepath extension.

Recognizes `.fna` → `:DNA`, `.frn` → `:RNA`, `.faa` → `:AA`.
Returns `nothing` for ambiguous extensions like `.fa`/`.fasta` or unknown suffixes.
If the path ends with `.gz`, the hint is derived from the underlying filename.
"""
function alphabet_hint_from_path(path::AbstractString)
    path_base = basename(path)
    if endswith(lowercase(path_base), ".gz")
        path_base = replace(path_base, ".gz" => "")
    end

    path_base = lowercase(path_base)
    if endswith(path_base, ".fna")
        return :DNA
    elseif endswith(path_base, ".frn")
        return :RNA
    elseif endswith(path_base, ".faa")
        return :AA
    end

    return nothing
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Detect sequence type from input and suggest appropriate file extension.

# Arguments
- `record`: A FASTA/FASTQ record
- `sequence`: A string or BioSequence containing sequence data

# Returns
- `String`: Suggested file extension:
  - ".fna" for DNA
  - ".frn" for RNA
  - ".faa" for protein
  - ".fa" for unrecognized sequences
"""
function detect_sequence_extension(record::Union{FASTX.FASTA.Record, FASTX.FASTQ.Record})
    return detect_sequence_extension(FASTX.sequence(record))
end
function detect_sequence_extension(sequence::AbstractString)
    sequence_type = detect_alphabet(sequence)
    return _detect_sequence_extension(sequence_type::Symbol)
end
function detect_sequence_extension(sequence::BioSequences.LongSequence)
    sequence_type = detect_alphabet(sequence)
    return _detect_sequence_extension(sequence_type::Symbol)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Internal helper function to convert sequence type to file extension.

Arguments
- sequence_type: Symbol representing sequence type (:DNA, :RNA, or :AA)

Returns
- String: Appropriate file extensions
"""
function _detect_sequence_extension(sequence_type::Symbol)
    @assert sequence_type in [:DNA, :RNA, :AA]
    if sequence_type == :DNA
        return ".fna"
    elseif sequence_type == :RNA
        return ".frn"
    elseif sequence_type == :AA
        return ".faa"
    else
        @warn "unrecognized sequence type: $(seq_type)"
        return ".fa"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compare two FASTA files to determine if they contain the same set of sequences,
regardless of sequence order.

# Arguments
- `fasta_1::String`: Path to first FASTA file
- `fasta_2::String`: Path to second FASTA file

# Returns
- `Bool`: `true` if both files contain exactly the same sequences, `false` otherwise

# Details
Performs a set-based comparison of DNA sequences by hashing each sequence.
Sequence order differences between files do not affect the result.
"""
function equivalent_fasta_sequences(fasta_1, fasta_2)
    fasta_1_hashes = Set(hash(BioSequences.LongDNA{2}(FASTX.sequence(record))) for record in Mycelia.open_fastx(fasta_1))
    fasta_2_hashes = Set(hash(BioSequences.LongDNA{2}(FASTX.sequence(record))) for record in Mycelia.open_fastx(fasta_2))
    @show setdiff(fasta_1_hashes, fasta_2_hashes)
    @show setdiff(fasta_2_hashes, fasta_1_hashes)
    return fasta_1_hashes == fasta_2_hashes
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Join fasta files while adding origin prefixes to the identifiers.

Does not guarantee uniqueness but will warn if conflicts arise
"""
function merge_fasta_files(;fasta_files, fasta_file)
    @info "merging $(length(fasta_files)) files..."
    identifiers = Set{String}()
    open(fasta_file, "w") do io
        fastx_io = FASTX.FASTA.Writer(io)
        ProgressMeter.@showprogress for f in fasta_files
            f_id = replace(basename(f), Mycelia.FASTA_REGEX => "")
            for record in Mycelia.open_fastx(f)
                new_record_id = f_id * "__" * FASTX.identifier(record)
                if new_record_id in identifiers
                    @warn "new identifier $(new_record_id) already in identifiers!!! Skipping"
                    continue
                else
                    push!(identifiers, new_record_id)
                    new_record = FASTX.FASTA.Record(new_record_id, FASTX.sequence(record))
                    write(fastx_io, new_record)
                end
            end
        end
        close(fastx_io)
    end
    @info "$(length(identifiers)) records merged..."
    return fasta_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Remove duplicate sequences from a FASTA file while preserving headers.

# Arguments
- `in_fasta`: Path to input FASTA file
- `out_fasta`: Path where deduplicated FASTA will be written

# Returns
Path to the output FASTA file (same as `out_fasta` parameter)

# Details
- Sequences are considered identical if they match exactly (case-sensitive)
- For duplicate sequences, keeps the first header encountered
- Input sequences are sorted by identifier before deduplication
- Preserves the original sequence formatting
"""
function deduplicate_fasta_file(in_fasta, out_fasta)
    fasta_df = fasta_to_table(collect(open_fastx(in_fasta)))
    sort!(fasta_df, "identifier")
    unique_sequences = DataFrames.combine(DataFrames.groupby(fasta_df, "sequence"), first)
    fasta = fasta_table_to_fasta(unique_sequences)
    open(out_fasta, "w") do io
        writer = FASTX.FASTA.Writer(io)
        for record in fasta
            write(writer, record)
        end
        close(writer)
    end
    return out_fasta
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Translates nucleic acid sequences from a FASTA file into amino acid sequences.

# Arguments
- `fasta_nucleic_acid_file::String`: Path to input FASTA file containing nucleic acid sequences
- `fasta_amino_acid_file::String`: Path where the translated amino acid sequences will be written

# Returns
- `String`: Path to the output FASTA file containing translated amino acid sequences
"""
function translate_nucleic_acid_fasta(fasta_nucleic_acid_file, fasta_amino_acid_file)
    open(fasta_amino_acid_file, "w") do io
        writer = FASTX.FASTA.Writer(io)
        for record in FASTX.FASTA.Reader(open(fasta_nucleic_acid_file))
            try
                raw_seq = FASTX.sequence(record)
                pruned_seq_length = Int(floor(length(raw_seq)/3)) * 3
                truncated_seq = raw_seq[1:pruned_seq_length]
                amino_acid_seq = BioSequences.translate(truncated_seq)
                amino_acid_record = FASTX.FASTA.Record(FASTX.identifier(record), FASTX.description(record), amino_acid_seq)
                write(writer, amino_acid_record)
            catch
                @warn "unable to translate record", record
            end
        end
        close(writer)
    end
    return fasta_amino_acid_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a FASTA file/record iterator to a DataFrame.

# Arguments
- `fasta`: FASTA record iterator from FASTX.jl

# Returns
- `DataFrame` with columns:
  - `identifier`: Sequence identifiers
  - `description`: Full sequence descriptions 
  - `sequence`: Biological sequences as strings
"""
function fasta_to_table(fasta)
    collected_fasta = collect(fasta)
    fasta_df = DataFrames.DataFrame(
        identifier = FASTX.identifier.(collected_fasta),
        description = FASTX.description.(collected_fasta),
        sequence = FASTX.sequence.(collected_fasta)
    )
    return fasta_df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a DataFrame containing FASTA sequence information into a vector of FASTA records.

# Arguments
- `fasta_df::DataFrame`: DataFrame with columns "identifier", "description", and "sequence"

# Returns
- `Vector{FASTX.FASTA.Record}`: Vector of FASTA records
"""
function fasta_table_to_fasta(fasta_df)
    records = Vector{FASTX.FASTA.Record}(undef, DataFrames.nrow(fasta_df))
    for (i, row) in enumerate(DataFrames.eachrow(fasta_df))
        record = FASTX.FASTA.Record(row["identifier"], row["description"], row["sequence"])
        records[i] = record
    end
    return records
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

This turns a 4-line FASTQ entry into a single tab separated line,
adds a column with the length of each read, passes it to Unix sort,
removes the length column, and converts it back into a FASTQ file.

sorts longest to shortest!!

http://thegenomefactory.blogspot.com/2012/11/sorting-fastq-files-by-sequence-length.html
"""
function sort_fastq(input_fastq, output_fastq="")
    
    if endswith(input_fastq, ".gz")
        p = pipeline(
                `gzip -dc $input_fastq`,
                `paste - - - -`,
                `perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;'`,
                `sort -nr`,
                `cut -f2-`,
                `tr "\t" "\n"`,
                `gzip`
                )
    else
        p = pipeline(
                `cat $input_fastq`,
                `paste - - - -`,
                `perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;'`,
                `sort -nr`,
                `cut -f2-`,
                `tr "\t" "\n"`
                )
    end
    run(pipeline(p, output_fastq))
    return output_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Counts the total number of records in a FASTA/FASTQ file.

# Arguments
- `fastx`: Path to a FASTA or FASTQ file (can be gzipped)

# Returns
- Number of records (sequences) in the file
"""
function count_records(fastx)
    n_records = 0
    for record in open_fastx(fastx)
        n_records += 1
    end
    return n_records
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate sequence lengths for reads in a FASTQ file.

# Arguments
- `fastq_file::String`: Path to input FASTQ file
- `total_reads::Integer=Inf`: Number of reads to process (defaults to all reads)

# Returns
- `Vector{Int}`: Array containing the length of each sequence read
"""
function determine_read_lengths(fastq_file; total_reads = Inf)
    if total_reads == Inf
        total_reads = count_records(fastq_file)
    end
    read_lengths = zeros(Int, total_reads)
    @info "determining read lengths"
    p = ProgressMeter.Progress(total_reads; dt=1)
    for (i, record) in enumerate(open_fastx(fastq_file))
#         push!(read_lengths, length(FASTX.sequence(record)))
        read_lengths[i] = length(FASTX.sequence(record))
        ProgressMeter.next!(p)
    end
    return read_lengths
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a Phred quality score (Q-value) to a probability of error.

# Arguments
- `q_value`: Phred quality score, typically ranging from 0 to 40

# Returns
- Error probability in range [0,1], where 0 indicates highest confidence

A Q-value of 10 corresponds to an error rate of 0.1 (10%), while a Q-value of 
30 corresponds to an error rate of 0.001 (0.1%).
"""
function q_value_to_error_rate(q_value)
    error_rate = 10^(q_value/(-10))
    return error_rate
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a sequencing error probability to a Phred quality score (Q-value).

The calculation uses the standard Phred formula: Q = -10 * log₁₀(error_rate)

# Arguments
- `error_rate::Float64`: Probability of error (between 0 and 1)

# Returns
- `q_value::Float64`: Phred quality score
"""
function error_rate_to_q_value(error_rate)
    q_value = -10 * log10(error_rate)
    return q_value
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the total size (in bases) of all sequences in a FASTA file.

# Arguments
- `fasta_file::AbstractString`: Path to the FASTA file

# Returns
- `Int`: Sum of lengths of all sequences in the FASTA file
"""
function total_fasta_size(fasta_file)
    return reduce(sum, map(record -> length(FASTX.sequence(record)), Mycelia.open_fastx(fasta_file)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Subsample reads from a FASTQ file using seqkit.

# Arguments
- `in_fastq::String`: Path to input FASTQ file
- `out_fastq::String=""`: Path to output FASTQ file. If empty, auto-generated based on input filename
- `n_reads::Union{Missing,Int}=missing`: Number of reads to sample
- `proportion_reads::Union{Missing,Float64}=missing`: Proportion of reads to sample (0.0-1.0)

# Returns
- `String`: Path to the output FASTQ file
"""
function subsample_reads_seqkit(;in_fastq::String, out_fastq::String="", n_reads::Union{Missing, Int}=missing, proportion_reads::Union{Missing, Float64}=missing)
    Mycelia.add_bioconda_env("seqkit")
    if ismissing(n_reads) && ismissing(proportion_reads)
        error("please specify the number or proportion of reads")
    elseif !ismissing(n_reads)
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit sample --two-pass --number $(n_reads) $(in_fastq)`, `gzip`)
        if isempty(out_fastq)
            out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqkit.N$(n_reads).fq.gz")
        end
    elseif !ismissing(proportion_reads)
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit sample --proportion $(proportion_reads) $(in_fastq)`, `gzip`)
        if isempty(out_fastq)
            out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqkit.P$(proportion_reads).fq.gz")
        end
    end
    @assert !isempty(out_fastq)
    if !isfile(out_fastq)
        run(pipeline(p, out_fastq))
    else
        @info "$(out_fastq) already present"
    end
    return out_fastq
end

# function subsample_reads_seqtk(;in_fastq::String, out_fastq="", n_reads::Union{Missing, Int}=missing, proportion_reads::Union{Missing, Float64}=missing)
#     Mycelia.add_bioconda_env("seqtk")
#     if ismissing(n_reads) && ismissing(proportion_reads)
#         error("please specify the number or proportion of reads")
#     elseif !ismissing(n_reads)
#         p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqtk seqtk sample -2 $(n_reads) $(in_fastq)`, `gzip`)
#         if isempty(out_fastq)
#             out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqtk.N$(n_reads).fq.gz")
#         end
#     elseif !ismissing(proportion_reads)
#         p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqtk seqtk sample -2 $(proportion_reads) $(in_fastq)`, `gzip`)
#         if isempty(out_fastq)
#             out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqtk.P$(proportion_reads).fq.gz")
#         end
#     end
#     @assert !isempty(out_fastq)
#     run(pipeline(p, out_fastq))
#     return out_fastq
# end

# subsample_reads_seqtk(in_fastq = fastq, n_reads=10)

# function filter_short_reads()
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates an index file (.fai) for a FASTA reference sequence using samtools.

The FASTA index allows efficient random access to the reference sequence. This is 
required by many bioinformatics tools that need to quickly fetch subsequences 
from the reference.

# Arguments
- `fasta`: Path to the input FASTA file

# Side Effects
- Creates a `{fasta}.fai` index file in the same directory as input
- Installs samtools via conda if not already present
"""
function samtools_index_fasta(;fasta)
    Mycelia.add_bioconda_env("samtools")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools faidx $(fasta)`)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Construct a FASTX FASTQ record from its components.

# Arguments
- `identifier::String`: The sequence identifier without the '@' prefix
- `sequence::String`: The nucleotide sequence
- `quality_scores::Vector{Int}`: Quality scores (0-93) as raw integers

# Returns
- `FASTX.FASTQRecord`: A parsed FASTQ record

# Notes
- Quality scores are automatically capped at 93 to ensure FASTQ compatibility
- Quality scores are converted to ASCII by adding 33 (Phred+33 encoding)
- The record is constructed in standard FASTQ format with four lines:
  1. Header line (@ + identifier)
  2. Sequence
  3. Plus line
  4. Quality scores (ASCII encoded)
"""
function fastq_record(;identifier, sequence, quality_scores)
    # Fastx wont parse anything higher than 93
    quality_scores = min.(quality_scores, 93)
    record_string = join(["@" * identifier, sequence, "+", join([Char(x+33) for x in quality_scores])], "\n")
    return FASTX.parse(FASTX.FASTQRecord, record_string)
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
function fastx_to_contig_lengths(fastx)
    OrderedCollections.OrderedDict(String(FASTX.identifier(record)) => length(FASTX.sequence(record)) for record in Mycelia.open_fastx(fastx))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Open and return a reader for FASTA or FASTQ format files.

# Arguments
- `path::AbstractString`: Path to input file. Can be:
    - Local file path
    - HTTP/FTP URL
    - Gzip compressed (.gz extension)

# Supported formats
- FASTA (.fasta, .fna, .faa, .fa)
- FASTQ (.fastq, .fq)

# Returns
- `FASTX.FASTA.Reader` for FASTA files
- `FASTX.FASTQ.Reader` for FASTQ files
"""
function open_fastx(path::AbstractString)
    if isfile(path)
        io = open(path)
    elseif occursin(r"^ftp", path) || occursin(r"^http", path)
        path = replace(path, r"^ftp:" => "http:")
        io = IOBuffer(HTTP.get(path).body)
    else
        error("unable to locate file $path")
    end
    path_base = basename(path)
    if occursin(r"\.gz$", path_base)
        io = CodecZlib.GzipDecompressorStream(io)
        path_base = replace(path_base, ".gz" => "")
    end
    if occursin(r"\.(fasta|fna|faa|fa|frn)$", path_base)
        fastx_io = FASTX.FASTA.Reader(io)
    elseif occursin(r"\.(fastq|fq)$", path_base)
        fastx_io = FASTX.FASTQ.Reader(io)
    else
        error("attempting to open a FASTX file with an unsupported extension: $(path_base)")
    end
    return fastx_io
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Writes FASTA records to a file, optionally gzipped.

# Arguments
- `outfile::AbstractString`: Path to the output FASTA file.  Will append ".gz" if `gzip` is true and ".gz" isn't already the extension.
- `records::Vector{FASTX.FASTA.Record}`: A vector of FASTA records.
- `gzip::Bool`: Optionally force compression of the output with gzip. By default will use the file name to infer.
- `show_progress::Union{Bool,Nothing}`: Whether to show a progress meter. If `nothing` (default), automatically shows progress for ≥1000 records.
- `use_pigz::Union{Bool,Nothing}`: If `true`, uses external pigz for compression (faster for large files). If `nothing` (default), automatically uses pigz for files estimated to be ≥50 MB if pigz is available. Requires pigz to be installed and in PATH.
- `pigz_threads::Int`: Number of threads for pigz to use. Defaults to all available cores.

# Returns
- `outfile::String`: The path to the output FASTA file (including ".gz" if applicable).
"""
function write_fasta(;
    outfile::AbstractString=tempname()*".fna",
    records::Vector{FASTX.FASTA.Record},
    gzip::Bool=false,
    show_progress::Union{Bool,Nothing}=nothing,
    use_pigz::Union{Bool,Nothing}=nothing,
    pigz_threads::Int=get_default_threads()
)
    # Determine if gzip compression should be used based on both the filename and the gzip argument
    gzip = occursin(r"\.gz$", outfile) || gzip

    # Append ".gz" to the filename if gzip is true and the extension isn't already present
    outfile = gzip && !occursin(r"\.gz$", outfile) ? outfile * ".gz" : outfile

    # Determine whether to show progress meter
    display_progress = isnothing(show_progress) ? length(records) >= 1000 : show_progress

    # Estimate file size based on record lengths (approximate)
    # Average FASTA overhead: ~60 bytes per record (identifier line + newlines)
    estimated_size_bytes = sum(length(FASTX.FASTA.sequence(record)) for record in records) + (length(records) * 60)

    # Auto-detect pigz usage: use for files ≥50 MB if pigz is available
    use_pigz_compression = if isnothing(use_pigz) && gzip
        # Check if pigz is available and file is large enough
        estimated_size_bytes >= 50_000_000
    elseif !isnothing(use_pigz)
        use_pigz
    else
        false
    end

    # If using pigz, write to temporary uncompressed file first
    temp_file = if gzip && use_pigz_compression
        # Remove .gz extension for temporary file
        replace(outfile, r"\.gz$" => "")
    else
        outfile
    end

    # Write the FASTA file
    open(temp_file, "w") do io
        if gzip && !use_pigz_compression
            io = CodecZlib.GzipCompressorStream(io)
        end

        FASTX.FASTA.Writer(io) do fastx_io
            if display_progress
                progress = ProgressMeter.Progress(length(records); desc="Writing FASTA records: ")
                for record in records
                    Base.write(fastx_io, record)
                    ProgressMeter.next!(progress)
                end
            else
                for record in records
                    Base.write(fastx_io, record)
                end
            end
        end
    end

    # Compress with pigz if requested
    if gzip && use_pigz_compression
        Mycelia.add_bioconda_env("pigz")
        if display_progress
            println("Compressing with pigz using $pigz_threads threads...")
        end
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $pigz_threads $temp_file`)
        # pigz automatically adds .gz extension and removes original
    end

    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a GenBank format file to FASTA format using EMBOSS seqret.

# Arguments
- `genbank`: Path to input GenBank format file
- `fasta`: Optional output FASTA file path (defaults to input path with .fna extension)
- `force`: If true, overwrites existing output file (defaults to false)

# Returns
Path to the output FASTA file

# Notes
- Requires EMBOSS suite (installed automatically via Conda)
- Will not regenerate output if it already exists unless force=true
"""
function genbank_to_fasta(;genbank, fasta=genbank * ".fna", force=false)
    add_bioconda_env("emboss")
    if !isfile(fasta) || force
        run(`$(Mycelia.CONDA_RUNNER) run -n emboss --live-stream seqret $(genbank) fasta:$(fasta)`)
    end
    return fasta
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get numerical PHRED quality scores from a FASTQ record.

This is a convenience wrapper around FASTX.quality_scores() that returns
the quality scores as a Vector{UInt8} representing PHRED scores.

# Arguments
- `record::FASTX.FASTQ.Record`: FASTQ record to extract quality scores from

# Returns
- `Vector{UInt8}`: PHRED quality scores (0-based, where 0 = lowest quality, 40+ = highest quality)

# Examples
```julia
record = FASTX.FASTQ.Record("read1", "ATCG", "IIII")
scores = get_phred_scores(record)  # Returns [40, 40, 40, 40]
```
"""
function get_phred_scores(record::FASTX.FASTQ.Record)::Vector{UInt8}
    return UInt8.(collect(FASTX.quality_scores(record)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert FASTQ quality string to numerical PHRED scores.

# Arguments
- `quality_string::AbstractString`: Quality string from FASTQ record (e.g., "IIII")

# Returns
- `Vector{UInt8}`: PHRED quality scores

# Examples
```julia
scores = quality_string_to_phred("IIII")  # Returns [40, 40, 40, 40]
scores = quality_string_to_phred("!#%+")  # Returns [0, 2, 4, 10]
```
"""
function quality_string_to_phred(quality_string::AbstractString)::Vector{UInt8}
    return UInt8[Int(c) - 33 for c in quality_string]
end
