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