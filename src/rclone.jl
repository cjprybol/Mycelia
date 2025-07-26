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