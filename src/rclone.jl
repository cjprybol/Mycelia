"""
    list_gdrive_with_links(remote::AbstractString; link_type::Symbol = :view, link_template::Union{Nothing,String}=nothing, rclone_args::Vector{String}=String[], return_df::Bool=true)

Run `rclone lsjson` on the given Google Drive remote path, parse the returned JSON, add `drive_id` and `drive_link` for each entry, and return either a Vector{Dict{String,Any}} or a DataFrames.DataFrame.

Keyword arguments
- link_type: one of `:view`, `:uc`, `:open`, or `:custom`.
  - :view  -> "https://drive.google.com/file/d/<id>/view?usp=sharing"
  - :uc    -> "https://drive.google.com/uc?id=<id>&export=download"
  - :open  -> "https://drive.google.com/open?id=<id>"
  - :custom -> `link_template` must be provided and contain "{id}"
- link_template: template string for custom links; "{id}" will be replaced with the Drive id.
- rclone_args: extra command-line arguments to pass to `rclone lsjson`.
- return_df: if true, return a DataFrames.DataFrame (requires DataFrames.jl); if false, return Vector{Dict{String,Any}}.

Notes
- When returning a DataFrame, absent values are represented with `missing` (so columns are join-friendly).
"""
function list_gdrive_with_links(remote::AbstractString; link_type::Symbol = :view, link_template::Union{Nothing,String}=nothing, rclone_args::Vector{String}=String[], return_df::Bool=true)
    # Build rclone command
    cmd_parts = ["rclone", "lsjson", remote]
    if !isempty(rclone_args)
        append!(cmd_parts, rclone_args)
    end
    cmd = Base.Cmd(cmd_parts)

    # Run rclone and capture stdout
    output = try
        read(cmd, String)
    catch err
        throw(ErrorException("Failed to run `rclone lsjson`: $(err)"))
    end

    # Parse JSON
    parsed = try
        JSON.parse(output)
    catch err
        throw(ErrorException("Failed to parse `rclone lsjson` output as JSON: $(err)"))
    end

    if !(isa(parsed, Vector))
        throw(ErrorException("Expected JSON array from `rclone lsjson`, got: $(typeof(parsed))"))
    end

    # Helper to build link from id
    function make_link_for_id(id::String)
        if link_type == :view
            return "https://drive.google.com/file/d/$(id)/view?usp=sharing"
        elseif link_type == :uc
            return "https://drive.google.com/uc?id=$(id)&export=download"
        elseif link_type == :open
            return "https://drive.google.com/open?id=$(id)"
        elseif link_type == :custom
            if link_template === nothing
                throw(ArgumentError("link_template must be provided when link_type == :custom"))
            end
            return replace(link_template, "{id}" => id)
        else
            if link_template !== nothing
                return replace(link_template, "{id}" => id)
            end
            throw(ArgumentError("Unknown link_type: $(link_type). Valid options: :view, :uc, :open, :custom"))
        end
    end

    # Common keys that may contain the Drive ID in rclone's lsjson output
    id_keys = ("ID", "Id", "id", "DriveId", "driveId")
    # id_keys = ("ID")

    entries = Vector{Dict{String,Any}}()
    for item in parsed
        if !isa(item, Dict)
            throw(ErrorException("Expected each entry from `rclone lsjson` to be an object/dict, found: $(typeof(item))"))
        end

        # Ensure keys are Strings and we have a mutable Dict
        entry = Dict{String,Any}(item)

        # Find ID-like key
        idval = nothing
        for k in id_keys
            if haskey(entry, k)
                v = entry[k]
                idval = v === nothing ? nothing : string(v)
                break
            end
        end

        # For DataFrame friendliness use `missing` instead of `nothing`
        drive_id_val = idval === nothing ? missing : idval
        drive_link_val = idval === nothing ? missing : make_link_for_id(idval)

        entry["drive_id"] = drive_id_val
        entry["drive_link"] = drive_link_val
        entry["Path"] = joinpath(remote, entry["Path"])

        push!(entries, entry)
    end

    if return_df
        # Convert vector of Dicts into a DataFrame
        df = DataFrames.DataFrame(entries)
        return df
    else
        return entries
    end
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

List all directories at the specified rclone path.

# Arguments
- `path::String`: Remote path to list directories from (e.g. "remote:/path/to/dir")

# Returns
- `Vector{String}`: Full paths to all directories found at the specified location
"""
function rclone_list_directories(path)
    directories = [join(split(line)[5:end], " ") for line in eachline(open(`rclone lsd $(path)`))]
    directories = joinpath.(path, directories)
    return directories
end

"""
    rclone_copy_list(;source::String, destination::String, relative_paths::Vector{String})

Copy a specific list of files from source to destination using rclone.

# Keywords
- `source::String`: Source location (can be remote like "gdrive:folder" or local path)
- `destination::String`: Destination location (can be remote or local path)
- `relative_paths::Vector{String}`: List of relative file paths to copy from source

# Returns
- `Bool`: `true` if transfer succeeded, `false` if an error occurred

# Details
This function creates a temporary file containing the list of paths and uses rclone's
`--files-from` option to copy only the specified files. This is more efficient than
copying files one by one when dealing with many files.

# Example
```julia
# Copy specific files from Google Drive to local directory
success = rclone_copy_list(
    source="gdrive:data",
    destination="/local/data",
    relative_paths=["file1.txt", "subdir/file2.csv", "images/photo.jpg"]
)
```

# Implementation Notes
- Creates destination directory if it's a local path and doesn't exist
- Uses verbose output to show transfer progress
- Automatically cleans up the temporary file list after completion
- Currently not called anywhere in the codebase (appears to be unused utility function)

# TODO
- Add support for additional rclone flags (bandwidth limits, chunk sizes, etc.)
- Add option to preserve directory structure or flatten it
- Add dry-run option for testing
"""
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
$(DocStringExtensions.TYPEDSIGNATURES)

Copy files between local and remote storage using rclone with automated retry logic.

# Arguments
- `source::String`: Source path or remote (e.g. "local/path" or "gdrive:folder")
- `dest::String`: Destination path or remote (e.g. "gdrive:folder" or "local/path")

# Keywords
- `config::String=""`: Optional path to rclone config file
- `max_attempts::Int=3`: Maximum number of retry attempts
- `sleep_timer::Int=60`: Initial sleep duration between retries in seconds (doubles after each attempt)

# Details
Uses optimized rclone settings for large files:
- 2GB chunk size
- 1TB upload cutoff
- Rate limited to 1 transaction per second
"""
function rclone_copy(source, dest; config="", max_attempts=3, sleep_timer=60)
    done = false
    attempts = 0
    while !done && attempts < max_attempts
        attempts += 1
        try
            # https://forum.rclone.org/t/google-drive-uploads-failing-http-429/34147/9
            # --tpslimit                                       Limit HTTP transactions per second to this
            # --drive-chunk-size SizeSuffix                    Upload chunk size (default 8Mi)
            # --drive-upload-cutoff SizeSuffix                 Cutoff for switching to chunked upload (default 8Mi)
            # not currently using these but they may become helpful
            # --drive-pacer-burst int                          Number of API calls to allow without sleeping (default 100)
            # --drive-pacer-min-sleep Duration                 Minimum time to sleep between API calls (default 100ms)
            if isempty(config)
                cmd = `rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $(source) $(dest)`
            else
                cmd = `rclone --config $(config) copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $(source) $(dest)`
            end
            @info "copying $(source) to $(dest) with command: $(cmd)"
            run(cmd)
            done = true
        catch
            @info "copying incomplete, sleeping $(sleep_timer) seconds and trying again..."
            sleep(sleep_timer)
            sleep_timer *= 2
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Copy files between local and remote storage using rclone with automated retry logic.

# Arguments
- `source::String`: Source path or remote (e.g. "local/path" or "gdrive:folder")
- `dest::String`: Destination path or remote (e.g. "gdrive:folder" or "local/path")

# Keywords
- `config::String=""`: Optional path to rclone config file
- `max_attempts::Int=3`: Maximum number of retry attempts
- `sleep_timer::Int=60`: Initial sleep duration between retries in seconds (doubles after each attempt)
- `includes::Vector{String}=[]`: One or more include patterns (each will be passed using `--include`)
- `excludes::Vector{String}=[]`: One or more exclude patterns (each will be passed using `--exclude`)
- `recursive::Bool=false`: If true, adds the flag for recursive traversal
"""
function rclone_copy2(source, dest;
                     config = "",
                     max_attempts = 3, sleep_timer = 60,
                     includes = String[],
                     excludes = String[],
                     recursive = false)
    done = false
    attempts = 0
    while !done && attempts < max_attempts
        attempts += 1
        try
            # Define base flags optimized for large files
            flags = ["--drive-chunk-size", "2G",
                     "--drive-upload-cutoff", "1T",
                     "--tpslimit", "1",
                     "--verbose"]

            # Append each include pattern with its flag
            for pattern in includes
                push!(flags, "--include")
                push!(flags, pattern)
            end

            # Append each exclude pattern with its flag
            for pattern in excludes
                push!(flags, "--exclude")
                push!(flags, pattern)
            end

            # Optionally add the recursive flag
            if recursive
                push!(flags, "--recursive")
            end

            # Build the full argument list as an array of strings.
            args = String[]
            # Add base command and optional config
            push!(args, "rclone")
            if !isempty(config)
                push!(args, "--config")
                push!(args, config)
            end
            push!(args, "copy")
            # Insert all flags (each flag and its parameter are separate elements)
            append!(args, flags)
            # Add source and destination paths
            push!(args, source)
            push!(args, dest)

            # Convert the argument vector into a Cmd object
            cmd = Cmd(args)

            @info "copying $(source) to $(dest) with command: $(cmd)"
            run(cmd)
            done = true
        catch e
            @info "copying incomplete, sleeping $(sleep_timer) seconds and trying again..."
            sleep(sleep_timer)
            sleep_timer *= 2
        end
    end
end