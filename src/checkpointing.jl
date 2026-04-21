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
    _atomic_save_dict(cache_file::String, data::Dict)

Write `data` to `cache_file` atomically via a PID-suffixed temp file and
`mv(...; force=true)`. The temp path includes the process ID to avoid
collisions between concurrent processes sharing a checkpoint directory.

If `jldsave` or `mv` throws (disk full, permission denied, interrupt), the
temp file is removed before the error is rethrown so stray `.tmp.<pid>`
files do not accumulate across failed runs.
"""
function _atomic_save_dict(cache_file::String, data::Dict)
    tmp = "$(cache_file).tmp.$(getpid())"
    try
        JLD2.jldsave(tmp; data = data)
        mv(tmp, cache_file; force = true)
    catch
        isfile(tmp) && rm(tmp; force = true)
        rethrow()
    end
    return nothing
end

"""
    cached_map(fn::Function, name::String, checkpoint_dir::String,
               inputs::AbstractVector;
               keyfn = string, save_every::Int = 50, force::Bool = false)

Per-input memoized map. Stores a single `Dict{String,Any}` at
`<checkpoint_dir>/<name>.jld2` keyed by `keyfn(input)`.

For each input:
- If `keyfn(input)` is already in the cached dict and `!force`, return cached.
- Otherwise, compute `fn(input)` and add to the dict.
- Atomic save every `save_every` newly-computed entries (tmp-file-then-rename).

Returns a `Vector` aligned with `inputs`. The map iteration is parallel-safe
via a `ReentrantLock` around cache writes; thread `fn` accordingly.

Unlike `cached_stage` (which is keyed by filename and treats the whole output
as one opaque blob), `cached_map` memoizes per input, so growing the input set
incrementally adds new entries without invalidating the existing cache. Use
`cached_map` for sample sweeps where the input list changes over time.

`keyfn(input)` must be deterministic and unique across inputs. Two distinct
inputs producing the same key silently overwrite each other.

!!! warning "Single-writer contract"
    `cached_map` is only safe for a single writer per `name`. The cache is
    read into memory once per call and fully overwritten on save, so two
    concurrent processes computing the same `name` will each start from the
    same on-disk state and the last save wins, silently dropping entries
    computed by the other process. If you need concurrent writers, shard by
    assigning each process a distinct `name` (e.g., suffix with worker id)
    and merge the caches after the sweep completes.

# Arguments
- `fn::Function`: Single-argument function mapping one input to one result
- `name::String`: Cache identifier (becomes filename)
- `checkpoint_dir::String`: Directory for the cache file
- `inputs::AbstractVector`: Inputs to map over

# Keywords
- `keyfn`: Function mapping an input to a `String` cache key. Default `string`.
- `save_every::Int`: Flush to disk every N newly-computed entries. Default 50.
- `force::Bool`: If true, ignore any existing cache and recompute everything.

# Example
```julia
import DrWatson
const CHECKPOINT_DIR = DrWatson.datadir("checkpoints")

results = Mycelia.cached_map("phage_mapping", CHECKPOINT_DIR, all_lims_ids) do id
    process_sample(id)
end

# Later, add more samples; only new ones are computed
results = Mycelia.cached_map("phage_mapping", CHECKPOINT_DIR, expanded_ids) do id
    process_sample(id)
end
```
"""
function cached_map(fn::Function, name::String, checkpoint_dir::String,
        inputs::AbstractVector;
        keyfn = string, save_every::Int = 50, force::Bool = false)
    mkpath(checkpoint_dir)
    cache_file = joinpath(checkpoint_dir, "$(name).jld2")

    cached = if isfile(cache_file) && !force
        try
            d = JLD2.load(cache_file, "data")
            if d isa Dict
                # Normalize to Dict{String,Any} so later writes don't fail via
                # `convert` when fn returns a type narrower/wider than what was
                # previously serialized (e.g., cache saved as Dict{String,Int}
                # rejects String values on the next run).
                normalized = Dict{String, Any}()
                for (k, v) in d
                    normalized[string(k)] = v
                end
                normalized
            else
                Logging.@warn "cached_map('$(name)'): cache at $(cache_file) is not a Dict, rebuilding"
                Dict{String, Any}()
            end
        catch err
            Logging.@warn "cached_map('$(name)'): cache at $(cache_file) unreadable ($err), rebuilding"
            Dict{String, Any}()
        end
    else
        Dict{String, Any}()
    end

    # Collect into a 1-based Vector so indexing is safe regardless of whether
    # `inputs` has offset indices (e.g., OffsetArrays).
    inputs_vec = collect(inputs)
    n_inputs = length(inputs_vec)
    keys_in_order = [string(keyfn(x)) for x in inputs_vec]
    missing_idx = [i for i in 1:n_inputs if !haskey(cached, keys_in_order[i])]

    Logging.@info "cached_map('$(name)'): $(n_inputs - length(missing_idx))/$(n_inputs) cached, $(length(missing_idx)) to compute"

    if isempty(missing_idx)
        return Any[cached[k] for k in keys_in_order]
    end

    cache_lock = Base.ReentrantLock()
    n_new = Threads.Atomic{Int}(0)
    # try/finally guarantees partial progress is flushed to disk even if `fn`
    # throws on some input. Without this, an exception in a long-running sweep
    # would lose every entry computed since the last `save_every` flush —
    # exactly the silent-data-loss mode this function was built to prevent.
    try
        Threads.@threads for i in missing_idx
            result = fn(inputs_vec[i])
            Base.lock(cache_lock) do
                cached[keys_in_order[i]] = result
                n = Threads.atomic_add!(n_new, 1) + 1
                if n % save_every == 0
                    _atomic_save_dict(cache_file, cached)
                end
            end
        end
    finally
        _atomic_save_dict(cache_file, cached)
        Logging.@info "cached_map('$(name)'): saved $(length(cached)) total entries to $(cache_file)"
    end
    return Any[cached[k] for k in keys_in_order]
end

# Alternative signature with name first (for consistency with do-block syntax)
function cached_map(name::String, checkpoint_dir::String,
        inputs::AbstractVector, fn::Function;
        keyfn = string, save_every::Int = 50, force::Bool = false)
    cached_map(fn, name, checkpoint_dir, inputs;
        keyfn = keyfn, save_every = save_every, force = force)
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
