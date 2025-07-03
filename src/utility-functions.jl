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