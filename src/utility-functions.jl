"""
Apply breadth-first sampling to a DataFrame
"""
function breadth_first_sample_dataframe(df::DataFrames.DataFrame, group_col::Union{Symbol, String}, 
                                       total_sample_size::Int; with_replacement::Bool = false)

    sampled_indices = breadth_first_sample(df[!, group_col], total_sample_size, 
                                         with_replacement=with_replacement)
    return df[sampled_indices, :]
end

"""
    choose_top_n_markers(N)

Return a vector of the top N most visually distinct marker symbols for plotting.

# Arguments
- `N::Int`: Number of distinct markers to return (max 17 for best differentiation).

# Returns
- `Vector{Symbol}`: Vector of marker symbol names.

# Example
    markers = choose_top_n_markers(7)
"""
function choose_top_n_markers(N::Int)
    # Priority order, most differentiable to least
    marker_priority = [
        :circle,
        :rect,
        :diamond,
        :cross,
        :utriangle,
        :dtriangle,
        :ltriangle,
        :rtriangle,
        :star5,
        :pentagon,
        :hexagon,
        :octagon,
        :xcross,
        :x,
        :star4,
        :star6,
        :star7,
        :star8,
        :heptagon,
        :pixel,
        :hline,
        :vline,
        :+,
    ]
    max_n = min(N, length(marker_priority))
    return marker_priority[1:max_n]
end

# """
#     zip_files(output_file::String, input_files::Vector{String})

# Creates a zip archive from input_files.
# Handles many files by writing a file list and using `zip -@`.
# Appends `.zip` if missing.
# """
# function zip_files(output_file::String, input_files::Vector{String})
#     out = endswith(output_file, ".zip") ? output_file : output_file * ".zip"
#     tmpfile = tempname()
#     open(tmpfile, "w") do io
#         for file in input_files
#             println(io, file)
#         end
#     end
#     cmd = `zip $out -@ < $tmpfile`
#     run(cmd)
#     rm(tmpfile)
#     return out
# end

"""
    tar_gz_files(output_file::String, input_files::Vector{String})

Creates a tar.gz archive from input_files.
Handles many files by writing a file list and using `tar -czf ... -T filelist`.
Appends `.tar.gz` if missing.
"""
function tar_gz_files(output_file::String, input_files::Vector{String})
    out = endswith(output_file, ".tar.gz") ? output_file : output_file * ".tar.gz"
    tmpfile = tempname()
    open(tmpfile, "w") do io
        for file in input_files
            println(io, file)
        end
    end
    cmd = `tar -czf $out -T $tmpfile`
    run(cmd)
    rm(tmpfile)
    return out
end

"""
    dictvec_to_dataframe(dictvec::Vector{<:AbstractDict}; symbol_columns::Bool = true)

Convert a vector of dictionaries (with possibly non-uniform keys and any key type) into a DataFrame.
Missing keys in a row are filled with `missing`.

# Arguments
- `dictvec`: Vector of dictionaries.
- `symbol_columns`: If true (default), columns are named as Symbols (when possible), else as raw keys.

# Returns
- `DataFrames.DataFrame` with columns as the union of all keys.
"""
function dictvec_to_dataframe(dictvec::Vector{<:AbstractDict}; symbol_columns::Bool = true)
    # Gather all unique keys from all dictionaries
    all_keys = Set{eltype(keys(dictvec[1]))}()
    for d in dictvec
        union!(all_keys, keys(d))
    end
    all_keys = collect(all_keys)

    # Optionally convert column names to Symbols if possible
    if symbol_columns
        try
            columns = Symbol.(all_keys)
        catch
            error("Some keys cannot be converted to Symbol. Set `symbol_columns=false`.")
        end
    else
        columns = all_keys
    end

    # Build rows as NamedTuples with missing for absent keys
    rows = [
        NamedTuple{Tuple(columns)}(
            map(k -> get(d, k, missing), all_keys)
        )
        for d in dictvec
    ]
    DataFrames.DataFrame(rows)
end

"""
Breadth-first sampling: sample at least one from each group,
then sample remaining proportionally to group frequencies
"""
function breadth_first_sample(group_vector, total_sample_size::Int; 
                            with_replacement::Bool = false)
    if total_sample_size <= 0
        throw(ArgumentError("total_sample_size must be positive"))
    end
    
    # Get unique groups and their counts
    group_counts = StatsBase.countmap(group_vector)
    unique_groups = collect(keys(group_counts))
    n_groups = length(unique_groups)
    
    if total_sample_size < n_groups
        throw(ArgumentError("total_sample_size ($total_sample_size) must be >= number of groups ($n_groups)"))
    end
    
    # Create indices for each group
    group_indices = Dict()
    for (i, group) in enumerate(group_vector)
        if !haskey(group_indices, group)
            group_indices[group] = Int[]
        end
        push!(group_indices[group], i)
    end
    
    # Step 1: Sample one from each group (breadth-first)
    sampled_indices = Int[]
    remaining_indices_by_group = Dict()
    
    for group in unique_groups
        group_idx = group_indices[group]
        # Sample one index from this group
        sampled_idx = StatsBase.sample(group_idx, 1)[1]
        push!(sampled_indices, sampled_idx)
        
        # Store remaining indices for this group
        if with_replacement
            remaining_indices_by_group[group] = group_idx
        else
            remaining_indices_by_group[group] = filter(x -> x != sampled_idx, group_idx)
        end
    end
    
    # Step 2: Calculate remaining sample size
    remaining_sample_size = total_sample_size - n_groups
    
    if remaining_sample_size > 0
        # Calculate total remaining population
        total_remaining = sum(length(indices) for indices in values(remaining_indices_by_group))
        
        if total_remaining == 0 && !with_replacement
            @warn "No remaining indices to sample from. Consider setting with_replacement=true"
            return sampled_indices
        end
        
        # Sample proportionally from remaining indices
        for group in unique_groups
            remaining_group_indices = remaining_indices_by_group[group]
            
            if length(remaining_group_indices) == 0
                continue  # Skip if no remaining indices in this group
            end
            
            # Calculate proportional sample size for this group
            if with_replacement
                group_proportion = group_counts[group] / length(group_vector)
            else
                group_proportion = length(remaining_group_indices) / total_remaining
            end
            
            additional_samples = round(Int, remaining_sample_size * group_proportion)
            
            if additional_samples > 0
                if with_replacement || length(remaining_group_indices) >= additional_samples
                    sampled = StatsBase.sample(remaining_group_indices, additional_samples, 
                                             replace=with_replacement)
                    append!(sampled_indices, sampled)
                else
                    # Take all remaining if not enough for desired sample size
                    append!(sampled_indices, remaining_group_indices)
                end
            end
        end
        
        # Handle any rounding discrepancies by sampling randomly from all remaining
        current_sample_size = length(sampled_indices)
        if current_sample_size < total_sample_size
            shortfall = total_sample_size - current_sample_size
            all_remaining = Int[]
            
            for indices in values(remaining_indices_by_group)
                if with_replacement
                    append!(all_remaining, indices)
                else
                    # Only include indices not already sampled
                    for idx in indices
                        if !(idx in sampled_indices)
                            push!(all_remaining, idx)
                        end
                    end
                end
            end
            
            if length(all_remaining) >= shortfall
                additional = StatsBase.sample(all_remaining, shortfall, replace=with_replacement)
                append!(sampled_indices, additional)
            elseif with_replacement && length(all_remaining) > 0
                additional = StatsBase.sample(all_remaining, shortfall, replace=true)
                append!(sampled_indices, additional)
            end
        end
    end
    
    return sampled_indices
end

function can_downcast_column(col::AbstractVector{T}, ::Type{S}) where {T<:AbstractFloat, S<:AbstractFloat}
    col_vec = collect(col)
    col2 = convert(Vector{S}, col_vec)
    return col_vec == convert(Vector{T}, col2)
end

function downcast_float_columns(df::DataFrames.DataFrame; target_type=Float32)
    for name in DataFrames.names(df)
        col = df[!, name]
        if eltype(col) <: AbstractFloat
            if can_downcast_column(col, target_type)
                # Always collect to Vector before converting types
                df[!, name] = convert(Vector{target_type}, collect(col))
            end
        end
    end
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a DataFrame to a JLD2 file using a standardized internal name.
"""
function JLD2_write_table(;df::DataFrames.DataFrame, filename::String)
    JLD2.jldopen(filename, "w") do file
        file["dataframe"] = df  # Always use the same internal name
    end
    return filename
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Read a DataFrame from a JLD2 file without needing to know the internal name.
If the file contains multiple DataFrames, returns the first one found.
"""
function JLD2_read_table(filename::String)
    df = JLD2.jldopen(filename, "r") do file
        # Try standard name first
        if haskey(file, "dataframe")
            return file["dataframe"]
        end
        
        # Otherwise search for any DataFrame
        for key in keys(file)
            if typeof(file[key]) <: DataFrames.DataFrame
                return file[key]
            end
        end
        
        # No DataFrame found
        error("No DataFrame found in file: $filename")
    end
    
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert all InlineString columns in a DataFrame to standard Strings.
Modifies the dataframe in-place and returns it.
"""
function sanitize_inline_strings!(df::DataFrames.DataFrame)
    for col in names(df)
        if eltype(df[!, col]) <: InlineStrings.InlineString
            df[!, col] = String.(df[!, col])
        end
    end
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a column to standard Strings if it contains InlineStrings,
otherwise return the original column unchanged.
"""
function sanitize_inline_strings(v::AbstractVector)
    if eltype(v) <: InlineStrings.InlineString
        return String.(v)
    else
        return v
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function dataframe_convert_dicts_to_json(df)
    df_copy = DataFrames.copy(df)
    for col in DataFrames.names(df_copy)
        if eltype(df_copy[!, col]) <: AbstractDict || any(x -> isa(x, AbstractDict), df_copy[!, col])
            df_copy[!, col] = [isa(cell, AbstractDict) ? JSON.json(cell) : cell for cell in df_copy[!, col]]
        end
    end
    return df_copy
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return a string representation of the vector `v` with each element on a new line,
mimicking valid Julia syntax. The output encloses the elements in square brackets
and separates them with a comma followed by a newline.
"""
function repr_long(v)
    buf = IOBuffer()
    println(buf, "[")
    for x in v
        println(buf, "    \"$x\",")
    end
    println(buf, "]")
    return String(take!(buf))
end


function rclone_copy_list(;source::String, destination::String, relative_paths::Vector{String})
    # Create a temporary file for storing file paths
    temp_file = joinpath(tempdir(), "rclone_sources_$(Random.randstring(8)).txt")
    
    try
        # Write paths to temp file
        open(temp_file, "w") do file
            for path in relative_paths
                println(file, path)
            end
        end
        
        println("Starting transfer of $(length(relative_paths)) files from $source to $destination...")
        
        # Make sure the destination directory exists if it's local
        if !occursin(":", destination) && !isdir(destination)
            mkdir(destination)
        end
        
        # Download files using rclone with progress reporting
        # --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 
        run(`rclone copy $source $destination --verbose --files-from $temp_file`)
        
        println("Download completed successfully")
        return true
    catch e
        println("Error: $e")
        return false
    finally
        # Clean up temp file
        if isfile(temp_file)
            rm(temp_file)
            println("Temporary file removed")
        end
    end
end
        
"""
    dataframe_to_ndjson(df::DataFrame; outfile::Union{String,Nothing}=nothing)

Converts a DataFrame `df` into a newline-delimited JSON (NDJSON) string.
Each line in the returned string represents one DataFrame row in JSON format,
suitable for upload to Google BigQuery.

# Keyword Arguments
- `outfile::Union{String,Nothing}`: If provided, writes the resulting NDJSON to the file path given.

# Examples
```julia
using DataFrames, Dates

# Sample DataFrame
df = DataFrame(
    id = [1, 2, 3],
    name = ["Alice", "Bob", "Carol"],
    created = [DateTime(2025, 4, 8, 14, 30), DateTime(2025, 4, 8, 15, 0), missing]
)

ndjson_str = dataframe_to_ndjson(df)
println(ndjson_str)

# Optionally, write to a file
dataframe_to_ndjson(df; outfile="output.ndjson")
"""
function dataframe_to_ndjson(df::DataFrames.DataFrame; outfile::Union{String, Nothing}=nothing) ndjson_lines = String[]
    # Iterate over each row in the DataFrame
    for row in eachrow(df)
        row_dict = Dict{String, Any}()
    
        # Build a dictionary for the current row
        for (col, value) in pairs(row)
            # Convert column names to strings
            col_name = string(col)
            if value === missing
                # Convert missing values to `nothing` which JSON prints as null
                row_dict[col_name] = nothing
            elseif isa(value, Dates.DateTime)
                # Format DateTime value in ISO 8601 format with millisecond precision
                # Modify the format string if a different format is needed for BigQuery.
                formatted = Dates.format(value, Dates.dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"
                row_dict[col_name] = formatted
            else
                row_dict[col_name] = value
            end
        end
    
        # Convert the dictionary to a JSON string and push into the list
        push!(ndjson_lines, JSON.json(row_dict))
    end
    
    # Join all JSON strings with newline delimiters (one JSON object per line)
    ndjson_str = join(ndjson_lines, "\n")
    
    # Optionally write the NDJSON content to a file if an outfile path is provided.
    if outfile !== nothing
        open(outfile, "w") do io
            write(io, ndjson_str)
        end
        return outfile
    else
        return ndjson_str
    end
end

function install_cloud_cli()
    # Get the user's home directory
    home_dir = ENV["HOME"]
    
    # Define the target installation directory ($HOME/google-cloud-sdk)
    sdk_install_dir = joinpath(home_dir, "google-cloud-sdk")
    sdk_bin_dir = joinpath(sdk_install_dir, "bin")
    
    # Function to update the PATH in current session and persist to .bashrc if needed.
    function update_path()
        # Update current session PATH
        if !occursin(sdk_bin_dir, ENV["PATH"])
            ENV["PATH"] = "$sdk_bin_dir:" * ENV["PATH"]
        end
        
        bashrc_file = joinpath(home_dir, ".bashrc")
        path_is_set = false
        if isfile(bashrc_file)
            for line in eachline(bashrc_file)
                if occursin("google-cloud-sdk/bin", line)
                    path_is_set = true
                    break
                end
            end
        end
        
        if !path_is_set
            println("Appending SDK bin path to $bashrc_file ...")
            open(bashrc_file, "a") do io
                println(io, "\n# Added by install_cloud_cli() on $(Dates.now())")
                println(io, "export PATH=\"$sdk_bin_dir:\$PATH\"")
            end
            println("Please reload your shell (e.g., run: source ~/.bashrc) to update your PATH.")
        else
            println("PATH already includes google-cloud-sdk/bin")
        end
    end
    
    # Check if the SDK is already installed.
    if isdir(sdk_install_dir)
        println("Google Cloud CLI is already installed.")
        update_path()
        println("No further installation required.")
        return
    end

    # Not installed: proceed with downloading and installing.
    tarball = "google-cloud-cli-linux-x86_64.tar.gz"
    sdk_url = "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/$tarball"
    
    println("Downloading Google Cloud CLI from: $sdk_url ...")
    run(`curl -O $sdk_url`)
    
    println("Extracting $tarball to $home_dir ...")
    run(`tar -xf $tarball -C $home_dir`)
    
    # Run the install script from within the extracted folder.
    install_script = joinpath(sdk_install_dir, "install.sh")
    println("Running the installation script ...")
    run(`bash $install_script --quiet --usage-reporting false --command-completion false --path-update false`)
    
    # Cleanup the downloaded tarball on successful install.
    println("Cleaning up the downloaded tarball...")
    rm(tarball, force=true)
    
    update_path()
    
    println("Installation complete. You now have gsutil and bq installed as part of the Google Cloud CLI.")
    println("Please run `gcloud init` or `gcloud auth login` and follow the interactive prompts to log in.")
end


function upload_dataframe_to_bigquery(;
    ndjson_file::String,
    project_id::String,
    dataset_id::String,
    table_id::String,
    gcs_bucket_name::String,
    bq_location::String, # e.g., "US"
    verbose=true
)
    full_table_id = "$(project_id):$(dataset_id).$(table_id)"
    if verbose
        println("Target BigQuery Table: ", full_table_id)
        println("Target GCS Bucket: ", gcs_bucket_name)
        println("BigQuery Location: ", bq_location)
    end
    
    try
        @assert isfile(ndjson_file) && filesize(ndjson_file) > 0
        verbose && println("NDJSON detected locally at: ", ndjson_file)

        # --- 3. Upload NDJSON to GCS ---
        verbose && println("Uploading to GCS bucket: ", gcs_bucket_name)
        gcs_uri = "gs://$(gcs_bucket_name)/$(basename(ndjson_file))"

        # Using gsutil command (simpler integration)
        upload_cmd = `gsutil cp $(ndjson_file) $(gcs_uri)`
        verbose && println("Running command: ", upload_cmd)
        run(upload_cmd)
        verbose && println("Upload to GCS complete: ", gcs_uri)

        # --- 4. Trigger BigQuery Load Job ---
        verbose && println("Triggering BigQuery load job for: ", full_table_id)

        # Using bq command (simpler integration)
        # Explicitly setting schema source_format and location
        load_cmd = `bq load --source_format=NEWLINE_DELIMITED_JSON --location=$(bq_location) $(full_table_id) $(gcs_uri)`
        # If you have a local schema file (e.g., schema.json):
        # load_cmd = `bq load --source_format=NEWLINE_DELIMITED_JSON --location=$(bq_location) $(full_table_id) $(gcs_uri) ./schema.json`
        # Or rely on auto-detect (use with caution):
        # load_cmd = `bq load --source_format=NEWLINE_DELIMITED_JSON --autodetect --location=$(bq_location) $(full_table_id) $(gcs_uri)`

        verbose && println("Running command: ", load_cmd)
        run(load_cmd)
        verbose && println("BigQuery load job initiated.")

        # --- 5. Cleanup GCS File (Optional) ---
        # Rely on bucket lifecycle rules or explicitly remove the temp file
        cleanup_cmd = `gsutil rm $(gcs_uri)`
        verbose && println("Running command: ", cleanup_cmd)
        run(cleanup_cmd)
        verbose && println("Cleaned up GCS file: ", gcs_uri)

        verbose && println("Process complete for table: ", full_table_id)

    catch e
        println("Error during process: ", e)
        # Add more robust error handling
    end
end

"""
    save_df_jld2(df::DataFrames.DataFrame, filename::String; key::String="dataframe")

Save a DataFrame to a JLD2 file.

# Arguments
- `df`: The DataFrame to save
- `filename`: Path to the JLD2 file (will add .jld2 extension if not present)
- `key`: The name of the dataset within the JLD2 file (defaults to "dataframe")

# Examples
```julia
import DataFrames
df = DataFrames.DataFrame(x = 1:3, y = ["a", "b", "c"])
save_df_jld2(df, "mydata")
```
"""
function save_df_jld2(;df::DataFrames.DataFrame, filename::String, key::String="dataframe")
    # Ensure filename has .jld2 extension
    if !endswith(lowercase(filename), ".jld2")
        filename = filename * ".jld2"
    end
    
    # Save the dataframe to the JLD2 file
    JLD2.jldopen(filename, "w") do file
        file[key] = df
    end
    
    return filename
end

"""
    load_df_jld2(filename::String; key::String="dataframe") -> DataFrames.DataFrame

Load a DataFrame from a JLD2 file.

# Arguments
- `filename`: Path to the JLD2 file (will add .jld2 extension if not present)
- `key`: The name of the dataset within the JLD2 file (defaults to "dataframe")

# Returns
- The loaded DataFrame

# Examples
```julia
df = load_df_jld2("mydata")
```
"""
function load_df_jld2(filename::String; key::String="dataframe")
    # Ensure filename has .jld2 extension
    if !endswith(lowercase(filename), ".jld2")
        filename = filename * ".jld2"
    end
    
    # Check if file exists
    if !isfile(filename)
        error("File not found: $filename")
    end
    
    # Load the dataframe from the JLD2 file
    df = JLD2.jldopen(filename, "r") do file
        if !haskey(file, key)
            error("Key '$key' not found in file $filename")
        end
        return file[key]
    end
    
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the SHA-256 hash of the contents of a file.

# Arguments
- `file::AbstractString`: The path to the file for which the SHA-256 hash is to be computed.

# Returns
- `String`: The SHA-256 hash of the file contents, represented as a hexadecimal string.
"""
function sha256_file(file::AbstractString)
    @assert isfile(file)
    return SHA.bytes2hex(SHA.sha256(read(file)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run a function `f` in parallel over a collection of `items` with a progress meter.

# Arguments
- `f::Function`: The function to be applied to each item in the collection.
- `items::AbstractVector`: A collection of items to be processed.

# Description
This function creates a progress meter to track the progress of processing each item in the `items` collection. 
It uses multithreading to run the function `f` on each item in parallel, updating the progress meter after each item is processed.
"""
function run_parallel_progress(f::Function, items::AbstractVector)
    # Create a progress meter
    p = ProgressMeter.Progress(length(items))
    
    # Create a lock and vector to store errors
    lock = ReentrantLock()
    errors = Vector{Union{Nothing, Tuple{Any, Any}}}(nothing, length(items))
    
    Threads.@threads for (i, item) in enumerate(items)
        try
            f(item)
        catch e
            lock() do
                errors[i] = (e, item)
            end
        end
        lock() do
            ProgressMeter.next!(p)
        end
    end
    return errors
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the current time as a Unix timestamp (seconds since epoch).

# Returns
- `Int`: Current time as an integer Unix timestamp (seconds since January 1, 1970 UTC)

# Examples
```julia
unix_time = current_unix_datetime()
# => 1709071368 (example value, will differ based on current time)
```
"""
function current_unix_datetime()
    return Int(floor(Dates.datetime2unix(Dates.now())))
end