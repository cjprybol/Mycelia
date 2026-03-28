# Checkpointing utilities for long-running genomics workflows
# Provides JLD2-based stage caching to enable iterative development
#
# Usage:
#   import DrWatson
#   const CHECKPOINT_DIR = DrWatson.datadir("checkpoints")
#
#   data = Mycelia.cached_stage("01_loaded", CHECKPOINT_DIR) do
#       expensive_load_operation()
#   end
#
#   # With input-dependent caching (invalidates when files change):
#   data = Mycelia.cached_stage("01_loaded", CHECKPOINT_DIR;
#       input_files=["data/samples.csv", "data/metadata.tsv"]) do
#       expensive_load_operation()
#   end
#
#   # To recompute: Mycelia.clear_stage("01_loaded", CHECKPOINT_DIR)

"""
    _input_hash(input_files::Vector{String})::String

Compute a hash from sorted input file paths, modification times, and sizes.
Returns a 16-character hex string suitable for use as a cache key suffix.
"""
function _input_hash(input_files::Vector{String})::String
    ctx = SHA.SHA2_256_CTX()
    for f in sort(input_files)
        isfile(f) || error("Input file not found: $f")
        SHA.update!(ctx, Vector{UInt8}(f))
        st = Base.stat(f)
        SHA.update!(ctx, Vector{UInt8}(string(st.mtime)))
        SHA.update!(ctx, Vector{UInt8}(string(st.size)))
    end
    return SHA.bytes2hex(SHA.digest!(ctx))[1:16]
end

"""
    cached_stage(compute_fn::Function, name::String, checkpoint_dir::String;
                 input_files::Vector{String}=String[])

Execute `compute_fn()` and cache result to `checkpoint_dir/name.jld2`.
Loads from cache if file exists. Delete checkpoint file to force re-computation.

When `input_files` is provided, the cache key incorporates a hash of the file
paths, modification times, and sizes. If any input file changes, the cache is
automatically invalidated and the computation re-runs.

# Arguments
- `compute_fn::Function`: Zero-argument function that computes the result
- `name::String`: Stage identifier (becomes filename)
- `checkpoint_dir::String`: Directory for checkpoint files
- `input_files::Vector{String}`: Optional list of input file paths to hash into cache key

# Example
```julia
import DrWatson
const CHECKPOINT_DIR = DrWatson.datadir("checkpoints")

# Basic usage (backward-compatible)
data = Mycelia.cached_stage("01_loaded", CHECKPOINT_DIR) do
    expensive_load_operation()
end

# Input-dependent caching
data = Mycelia.cached_stage("01_loaded", CHECKPOINT_DIR;
    input_files=["data/samples.csv"]) do
    expensive_load_operation()
end
```
"""
function cached_stage(compute_fn::Function, name::String, checkpoint_dir::String;
        input_files::Vector{String} = String[])
    mkpath(checkpoint_dir)
    suffix = isempty(input_files) ? "" : "_$(_input_hash(input_files))"
    cache_file = joinpath(checkpoint_dir, "$(name)$(suffix).jld2")

    if isfile(cache_file)
        Logging.@info "Loading cached '$(name)' from $(cache_file)"
        return JLD2.load(cache_file, "data")
    end

    # Clean up stale hash-suffixed cache files for this stage name
    if !isempty(input_files)
        for f in readdir(checkpoint_dir)
            if startswith(f, "$(name)_") && endswith(f, ".jld2")
                old_file = joinpath(checkpoint_dir, f)
                if old_file != cache_file
                    rm(old_file)
                    Logging.@info "Removed stale cache: $(f)"
                end
            end
        end
    end

    Logging.@info "Computing '$(name)'..."
    result = compute_fn()
    JLD2.jldsave(cache_file; data = result)
    Logging.@info "Saved '$(name)' to $(cache_file)"
    return result
end

# Alternative signature with name first (for consistency with do-block syntax)
function cached_stage(name::String, checkpoint_dir::String, compute_fn::Function;
        input_files::Vector{String} = String[])
    cached_stage(compute_fn, name, checkpoint_dir; input_files = input_files)
end

"""
    clear_stage(name::String, checkpoint_dir::String;
                input_files::Vector{String}=String[]) -> Bool

Delete a specific checkpoint to force re-computation.
Returns true if checkpoint was deleted, false if not found.

When `input_files` is provided, clears the hash-suffixed cache file matching
those inputs. Without `input_files`, clears the base cache file.

# Example
```julia
Mycelia.clear_stage("02_processed", CHECKPOINT_DIR)
Mycelia.clear_stage("02_processed", CHECKPOINT_DIR; input_files=["data/samples.csv"])
```
"""
function clear_stage(name::String, checkpoint_dir::String;
        input_files::Vector{String} = String[])::Bool
    suffix = isempty(input_files) ? "" : "_$(_input_hash(input_files))"
    cache_file = joinpath(checkpoint_dir, "$(name)$(suffix).jld2")
    if isfile(cache_file)
        rm(cache_file)
        Logging.@info "Cleared checkpoint: $(name)"
        return true
    end
    Logging.@warn "Checkpoint not found: $(name)"
    return false
end

"""
    clear_all_stages(checkpoint_dir::String) -> Int

Clear all checkpoints in directory (use with caution).
Returns the number of checkpoints deleted.

# Example
```julia
Mycelia.clear_all_stages(CHECKPOINT_DIR)  # Start fresh
```
"""
function clear_all_stages(checkpoint_dir::String)::Int
    if !isdir(checkpoint_dir)
        Logging.@warn "Checkpoint directory not found: $(checkpoint_dir)"
        return 0
    end

    count = 0
    for f in readdir(checkpoint_dir)
        if endswith(f, ".jld2")
            rm(joinpath(checkpoint_dir, f))
            count += 1
        end
    end
    Logging.@info "Cleared $(count) checkpoints from $(checkpoint_dir)"
    return count
end

"""
    list_stages(checkpoint_dir::String) -> Vector{String}

List all cached stages in the checkpoint directory.
Returns stage names (without .jld2 extension) sorted alphabetically.

# Example
```julia
stages = Mycelia.list_stages(CHECKPOINT_DIR)
# ["01_loaded", "02_processed", "03_filtered"]
```
"""
function list_stages(checkpoint_dir::String)::Vector{String}
    if !isdir(checkpoint_dir)
        return String[]
    end

    stages = String[]
    for f in readdir(checkpoint_dir)
        if endswith(f, ".jld2")
            push!(stages, replace(f, ".jld2" => ""))
        end
    end
    return sort(stages)
end

"""
    checkpoint_info(checkpoint_dir::String) -> DataFrame

Get information about all checkpoints: name, size, modification time.
Useful for understanding checkpoint state and cleaning up old stages.

# Example
```julia
info = Mycelia.checkpoint_info(CHECKPOINT_DIR)
# DataFrame with columns: name, size_mb, modified
```
"""
function checkpoint_info(checkpoint_dir::String)
    if !isdir(checkpoint_dir)
        return DataFrames.DataFrame(name = String[], size_mb = Float64[], modified = Dates.DateTime[])
    end

    names = String[]
    sizes = Float64[]
    modified = Dates.DateTime[]

    for f in readdir(checkpoint_dir)
        if endswith(f, ".jld2")
            fpath = joinpath(checkpoint_dir, f)
            push!(names, replace(f, ".jld2" => ""))
            push!(sizes, round(filesize(fpath) / 1e6; digits = 2))
            push!(modified, Dates.unix2datetime(mtime(fpath)))
        end
    end

    return DataFrames.DataFrame(
        name = names,
        size_mb = sizes,
        modified = modified
    )
end
