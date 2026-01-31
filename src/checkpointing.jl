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
#   # To recompute: Mycelia.clear_stage("01_loaded", CHECKPOINT_DIR)

"""
    cached_stage(name::String, compute_fn::Function, checkpoint_dir::String)

Execute `compute_fn()` and cache result to `checkpoint_dir/name.jld2`.
Loads from cache if file exists. Delete checkpoint file to force re-computation.

# Arguments
- `name::String`: Stage identifier (becomes filename)
- `compute_fn::Function`: Zero-argument function that computes the result
- `checkpoint_dir::String`: Directory for checkpoint files

# Example
```julia
import DrWatson
const CHECKPOINT_DIR = DrWatson.datadir("checkpoints")

data = Mycelia.cached_stage("01_loaded", CHECKPOINT_DIR) do
    expensive_load_operation()
end

# With explicit function
function load_data()
    CSV.read("large_file.csv", DataFrames.DataFrame)
end
data = Mycelia.cached_stage("01_loaded", load_data, CHECKPOINT_DIR)
```
"""
function cached_stage(compute_fn::Function, name::String, checkpoint_dir::String)
    mkpath(checkpoint_dir)
    cache_file = joinpath(checkpoint_dir, "$(name).jld2")

    if isfile(cache_file)
        Logging.@info "Loading cached '$(name)' from $(cache_file)"
        return JLD2.load(cache_file, "data")
    end

    Logging.@info "Computing '$(name)'..."
    result = compute_fn()
    JLD2.jldsave(cache_file; data=result)
    Logging.@info "Saved '$(name)' to $(cache_file)"
    return result
end

# Alternative signature with name first (for consistency with do-block syntax)
function cached_stage(name::String, checkpoint_dir::String, compute_fn::Function)
    cached_stage(compute_fn, name, checkpoint_dir)
end

"""
    clear_stage(name::String, checkpoint_dir::String) -> Bool

Delete a specific checkpoint to force re-computation.
Returns true if checkpoint was deleted, false if not found.

# Example
```julia
Mycelia.clear_stage("02_processed", CHECKPOINT_DIR)
# Re-run the cell to recompute
```
"""
function clear_stage(name::String, checkpoint_dir::String)::Bool
    cache_file = joinpath(checkpoint_dir, "$(name).jld2")
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
        return DataFrames.DataFrame(name=String[], size_mb=Float64[], modified=DateTime[])
    end

    names = String[]
    sizes = Float64[]
    modified = DateTime[]

    for f in readdir(checkpoint_dir)
        if endswith(f, ".jld2")
            fpath = joinpath(checkpoint_dir, f)
            push!(names, replace(f, ".jld2" => ""))
            push!(sizes, round(filesize(fpath) / 1e6; digits=2))
            push!(modified, Dates.unix2datetime(mtime(fpath)))
        end
    end

    return DataFrames.DataFrame(
        name=names,
        size_mb=sizes,
        modified=modified
    )
end
