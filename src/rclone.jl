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