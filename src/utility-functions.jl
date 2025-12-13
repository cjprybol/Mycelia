# --- Helper function to generate ellipse points ---
function get_ellipse_points(x_coords, y_coords; scale=3.5)
    if length(x_coords) < 3
        return CairoMakie.Point2f[], Dict()
    end
    
    cov_matrix = Statistics.cov(hcat(x_coords, y_coords))
    mean_vec = [Statistics.mean(x_coords), Statistics.mean(y_coords)]
    
    eigen = LinearAlgebra.eigen(cov_matrix)
    vals, vecs = eigen.values, eigen.vectors
    
    angle = atan(vecs[2, 2], vecs[1, 2])
    major_axis = scale * sqrt(vals[2])
    minor_axis = scale * sqrt(vals[1])
    
    t = range(0, 2π, length=100)
    ellipse_x = major_axis * cos.(t)
    ellipse_y = minor_axis * sin.(t)
    
    rotated_x = ellipse_x * cos(angle) - ellipse_y * sin(angle) .+ mean_vec[1]
    rotated_y = ellipse_x * sin(angle) + ellipse_y * cos(angle) .+ mean_vec[2]
    
    points = [CairoMakie.Point2f(x, y) for (x, y) in zip(rotated_x, rotated_y)]
    
    anchors = Dict(
        :top => points[argmax(rotated_y)],
        :bottom => points[argmin(rotated_y)],
        :left => points[argmin(rotated_x)],
        :right => points[argmax(rotated_x)]
    )
    
    return points, anchors
end


"""
Collapse duplicate rows in a DataFrame by consolidating non-missing values.

For each group defined by `grouping_col`, this function attempts to merge rows
by taking the first non-missing value for each column. If conflicting non-missing
values are found in any column, a warning is issued and all conflicting rows are kept.
"""
function collapse_duplicates(df, grouping_col)
    result_rows = DataFrames.DataFrame[]
    conflicts = []
    
    for g in DataFrames.groupby(df, grouping_col)
        if DataFrames.nrow(g) == 1
            push!(result_rows, g)
            continue
        end
        
        # Check for conflicts across all columns
        has_conflict = false
        conflicting_cols = String[]
        
        for col in DataFrames.names(g)
            if col == grouping_col
                continue
            end
            
            non_missing_values = filter(!ismissing, g[!, col])
            unique_values = unique(non_missing_values)
            
            if length(unique_values) > 1
                has_conflict = true
                push!(conflicting_cols, col)
            end
        end
        
        if has_conflict
            # Keep all rows due to conflict
            group_id = g[1, grouping_col]
            @warn "Conflict found in group $group_id in columns: $(join(conflicting_cols, ", ")). Keeping all rows."
            push!(conflicts, (group_id, conflicting_cols))
            push!(result_rows, g)
        else
            # Merge rows by taking first non-missing value
            merged_row = DataFrames.DataFrame()
            for col in DataFrames.names(g)
                merged_row[!, col] = [Base.coalesce(g[!, col]...)]
            end
            push!(result_rows, merged_row)
        end
    end
    
    collapsed_df = DataFrames.vcat(result_rows...)
    
    if !isempty(conflicts)
        println("\n=== Summary of Conflicts ===")
        println("Total groups with conflicts: $(length(conflicts))")
        for (group_id, cols) in conflicts
            println("  - $group_id: conflicts in $(join(cols, ", "))")
        end
    end
    
    return collapsed_df
end

"""
    with_retry(f; max_attempts=3, initial_delay=5.0, max_delay=120.0, backoff_factor=2.0, on_retry=nothing)

Execute a function with automatic retry logic and exponential backoff.

Useful for wrapping unreliable operations like network requests or API calls
that may fail due to transient issues.

# Arguments
- `f`: Zero-argument function to execute

# Keywords
- `max_attempts::Int=3`: Maximum number of attempts before giving up
- `initial_delay::Float64=5.0`: Initial delay in seconds between retries
- `max_delay::Float64=120.0`: Maximum delay between retries (caps exponential growth)
- `backoff_factor::Float64=2.0`: Multiplier for delay after each failed attempt
- `on_retry::Union{Function,Nothing}=nothing`: Optional callback `(attempt, exception, delay) -> nothing` 
  called before each retry sleep

# Returns
- Result of `f()` on success

# Throws
- Rethrows the last exception if all attempts fail

# Example
```julia
# Retry a download up to 3 times with exponential backoff
result = with_retry(max_attempts=5, initial_delay=10.0) do
    run(`curl -O https://example.com/large_file.zip`)
end

# With custom retry callback for logging
with_retry(on_retry=(att, ex, delay) -> @warn "Attempt \$att failed, retrying in \$(delay)s") do
    unreliable_api_call()
end
```

See also: [`ncbi_genome_download_accession`](@ref) which uses this for robust genome downloads.
"""
function with_retry(f; 
                    max_attempts::Int=3, 
                    initial_delay::Float64=5.0, 
                    max_delay::Float64=120.0, 
                    backoff_factor::Float64=2.0,
                    on_retry::Union{Function,Nothing}=nothing)
    local last_exception
    delay = initial_delay
    
    for attempt in 1:max_attempts
        try
            return f()
        catch e
            last_exception = e
            if attempt < max_attempts
                if on_retry !== nothing
                    on_retry(attempt, e, delay)
                else
                    @warn "Attempt $attempt/$max_attempts failed, retrying in $(delay)s..." exception=(e, catch_backtrace())
                end
                sleep(delay)
                delay = min(delay * backoff_factor, max_delay)
            end
        end
    end
    
    @error "All $max_attempts attempts failed" exception=(last_exception, catch_backtrace())
    throw(last_exception)
end

function recursively_list_directories(dir::AbstractString)
    directories_list = String[]
    for (root, directories, files) in Base.walkdir(dir)
        for d in directories
            directory_path = joinpath(root, d)
            # @assert isdir(directory_path)
            push!(directories_list, directory_path)
        end
    end
    return directories_list
end

function recursively_list_files(dir::AbstractString)
    files_list = String[]
    for (root, directories, files) in Base.walkdir(dir)
        for f in files
            file_path = joinpath(root, f)
            # @assert isfile(file_path)
            push!(files_list, file_path)
        end
    end
    return files_list
end

"""
    read_xlsx(filename::AbstractString) -> NamedTuple

Read all sheets from an XLSX file into a NamedTuple of DataFrames.

Each sheet becomes a field in the NamedTuple where the field name is the sheet name
and the value is the corresponding DataFrame.

# Arguments
- `filename::AbstractString`: Path to the XLSX file

# Returns
- `NamedTuple`: Named tuple with sheet names as keys and DataFrames as values
"""
function read_xlsx(filename::AbstractString)
    
    # Ensure filename ends with .xlsx
    filename = _ensure_xlsx_extension(filename)
    
    # Open the XLSX file
    xf = XLSX.readxlsx(filename)
    
    # Get all sheet names
    sheet_names = XLSX.sheetnames(xf)
    
    # Read each sheet into a DataFrame
    dataframes = []
    symbols = Symbol[]
    
    for sheet_name in sheet_names
        df = DataFrames.DataFrame(XLSX.readtable(filename, sheet_name))
        push!(dataframes, df)
        push!(symbols, Symbol(sheet_name))
    end
    
    # Create and return named tuple
    return NamedTuple{Tuple(symbols)}(dataframes)
end

"""
    write_xlsx(filename::AbstractString, dataframes...)

Write one or more DataFrames to an XLSX file with multiple sheets.

# Arguments
- `filename::AbstractString`: Path for the output XLSX file
- `dataframes...`: Variable number of arguments that can be DataFrames or Pairs of sheet_name => dataframe
"""
function write_xlsx(filename::AbstractString, dataframes...)
    
    # Ensure filename ends with .xlsx
    filename = _ensure_xlsx_extension(filename)
    
    # Process arguments to extract sheet names and dataframes
    sheets = _process_dataframe_args(dataframes...)
    
    # Create new XLSX file
    XLSX.openxlsx(filename, mode="w") do xf
        for (sheet_name, df) in sheets
            sheet = XLSX.addsheet!(xf, sheet_name)
            
            # Write headers
            col_names = DataFrames.names(df)
            for (col_idx, col_name) in enumerate(col_names)
                sheet[1, col_idx] = string(col_name)
            end
            
            # Write data
            for row_idx in 1:DataFrames.nrow(df)
                for (col_idx, col_name) in enumerate(col_names)
                    value = df[row_idx, col_name]
                    # Handle missing values
                    if DataFrames.ismissing(value)
                        sheet[row_idx + 1, col_idx] = ""
                    else
                        sheet[row_idx + 1, col_idx] = value
                    end
                end
            end
        end
    end
    
    println("Successfully wrote $(length(sheets)) sheet(s) to $filename")
end

"""
    write_xlsx_single(filename::AbstractString, dataframe; sheet_name::AbstractString="Sheet1")

Convenience function to write a single DataFrame to an XLSX file.

# Arguments
- `filename::AbstractString`: Path for the output XLSX file
- `dataframe`: DataFrame to write
- `sheet_name::AbstractString`: Name for the sheet (default: "Sheet1")
"""
function write_xlsx_single(filename::AbstractString, dataframe; sheet_name::AbstractString="Sheet1")
    write_xlsx(filename, sheet_name => dataframe)
end

"""
    write_xlsx_single(filename::AbstractString, sheet_pair::Pair)

Convenience function to write a single DataFrame with custom sheet name using pair syntax.

# Arguments
- `filename::AbstractString`: Path for the output XLSX file
- `sheet_pair::Pair`: Pair of sheet_name => dataframe
"""
function write_xlsx_single(filename::AbstractString, sheet_pair::Pair)
    write_xlsx(filename, sheet_pair)
end

# Helper functions

"""
Internal function to ensure filename has .xlsx extension and warn about other extensions.
"""
function _ensure_xlsx_extension(filename::AbstractString)
    # Check for other common tabular data extensions
    other_extensions = [".csv", ".tsv", ".txt", ".tab"]
    filename_lower = lowercase(filename)
    
    for ext in other_extensions
        if endswith(filename_lower, ext)
            @warn "File extension '$ext' detected. This function works with Excel files (.xlsx). " *
                  "Converting extension to .xlsx"
            filename = filename[1:end-length(ext)] * ".xlsx"
            return filename
        end
    end
    
    # Add .xlsx if not present
    if !endswith(filename_lower, ".xlsx")
        filename = filename * ".xlsx"
    end
    
    return filename
end

"""
Internal function to process variable arguments for write_xlsx function.
"""
function _process_dataframe_args(dataframes...)
    
    sheets = Pair{String, DataFrames.DataFrame}[]
    auto_sheet_counter = 1
    
    for arg in dataframes
        if isa(arg, DataFrames.DataFrame)
            # Auto-generate sheet name
            sheet_name = "Sheet$auto_sheet_counter"
            push!(sheets, sheet_name => arg)
            auto_sheet_counter += 1
        elseif isa(arg, Pair)
            # Custom sheet name provided
            sheet_name, df = arg
            if !isa(df, DataFrames.DataFrame)
                error("Value in pair must be a DataFrame, got $(typeof(df))")
            end
            push!(sheets, string(sheet_name) => df)
        else
            error("Arguments must be DataFrames or Pairs of sheet_name => dataframe, got $(typeof(arg))")
        end
    end
    
    if isempty(sheets)
        error("At least one DataFrame must be provided")
    end
    
    return sheets
end


"""
    sanitize_for_arrow(df::DataFrames.DataFrame) -> DataFrames.DataFrame

Creates a new DataFrame with columns sanitized for Arrow compatibility.

It identifies columns with abstract or mixed types (e.g., `Vector{Any}`)
and converts them to the most specific, concrete type possible.

- If types are compatible (e.g., `Int` and `Float64`), they are promoted.
- If types are incompatible (e.g., `Int` and `String`), the column is
  converted to `String` and a warning is issued.
"""
function sanitize_for_arrow(df::DataFrames.DataFrame)
    sanitized_df = copy(df) # Work on a copy
    
    for col_name in names(sanitized_df)
        col = sanitized_df[!, col_name]
        
        # Only process columns with abstract element types
        if !isconcretetype(eltype(col))
            # Get all unique types present in the column, ignoring missings
            present_types = unique(typeof.(skipmissing(col)))
            
            if isempty(present_types)
                # Column is all `missing`, which is fine.
                continue
            end
            
            target_type = nothing
            try
                # Attempt to promote all found types to a common supertype
                target_type = reduce(promote_type, present_types)
            catch
                # Promotion failed (e.g., trying to promote Int and String)
                # Fall back to String as the only safe option
                target_type = String
                @warn "Column '$col_name' has incompatible mixed types. Converting to String."
            end
            
            # Create a new column by converting each element to the target type
            new_col = map(col) do val
                ismissing(val) ? missing : convert(target_type, val)
            end
            
            sanitized_df[!, col_name] = new_col
        end
    end
    
    return sanitized_df
end

"""
    write_arrow(
        df::DataFrames.DataFrame;
        filename::String,
        compress::Symbol=:zstd,
        force::Bool=false,
        sanitize::Bool=true
    ) -> String

Writes a DataFrame to an Apache Arrow file with optional pre-sanitization.

# Keyword Arguments
- `filename::String`: The path for the output Arrow file.
- `compress::Symbol=:zstd`: Compression algorithm (:zstd, :lz4, :gzip, or nothing).
- `force::Bool=false`: If `true`, overwrite an existing file.
- `sanitize::Bool=true`: If `true`, automatically run `sanitize_for_arrow`
  to resolve mixed-type columns before writing. Recommended to leave on.
"""
function write_arrow(
    df::DataFrames.DataFrame;
    filename::String,
    compress::Symbol=:zstd,
    force::Bool=false,
    verbose::Bool=false,
    sanitize::Bool=true # New keyword for safety
)
    if isfile(filename) && !force
        @warn "File '$filename' already exists. Use force=true to overwrite."
        return filename
    end
    
    # Use a temporary variable to hold the DataFrame to be written
    df_to_write = df

    # Sanitize the DataFrame by default or if requested
    if sanitize
        df_to_write = sanitize_for_arrow(df)
    end

    # Arrow.write handles everything else
    Arrow.write(filename, df_to_write, compress=compress)
    
    verbose && println("Successfully wrote DataFrame to $filename")
    return filename
end

"""
    read_arrow(filename::String) -> DataFrames.DataFrame

Reads an Apache Arrow file into a Julia DataFrame.

# Arguments
- `filename::String`: The path to the Arrow file.

# Returns
- `DataFrames.DataFrame`: The loaded DataFrame.
"""
function read_arrow(filename::String)
    if !isfile(filename)
        error("File not found: $filename")
    end
    
    # Reading is a simple two-step process: open a table, then convert to a DataFrame.
    arrow_table = Arrow.Table(filename)
    df = DataFrames.DataFrame(arrow_table)
    
    return df
end


# """
#     @recordtest call_expr

# Evaluate the function call `call_expr`, print a test statement in the form `@test call_expr == result` using the actual result, and return the result.

# # Example

# The printed `@test` statement can be copy-pasted into a test file to verify that the function with the given arguments returns the same result in the future.

# !!! note
# Arguments should be literals or global variables for the printed test to be self-contained and copy-pasteable.
# """
# macro recordtest(call_expr)
#     # Evaluate the function call and capture the result
#     result = eval(call_expr)
#     # Get the source code for the call expression
#     call_str = string(call_expr)
#     # Generate the test statement as a string
#     test_str = "Test.@test $call_str == $(repr(result))"
#     # Print it out (or you could `return test_str` to use elsewhere)
#     println(test_str)
#     return result  # Optionally return the result for interactive feedback
# end

"""
    _get_output_files(output_dir::AbstractString) -> Vector{String}

Get a list of all files in the output directory.
"""
function _get_output_files(output_dir::AbstractString)
    if !Base.isdir(output_dir)
        return String[]
    end
    
    files = String[]
    for (root, dirs, filenames) in Base.walkdir(output_dir)
        for filename in filenames
            push!(files, Base.joinpath(root, filename))
        end
    end
    
    return files
end

"""
    write_tsvgz(;df::DataFrames.DataFrame, filename::String, force::Bool=false, bufsize::Int=2*1024^3, buffer_in_memory::Bool=false)

Write a DataFrame to a gzipped TSV (.tsv.gz) file with robust stream and buffer handling.

# Arguments
- `df`: The DataFrame to write.
- `filename`: Output path, extension will be enforced as `.tsv.gz`.
- `force`: If true, overwrite an existing non-empty file.
- `bufsize`: Buffer size in bytes for both the Gzip stream and CSV write (default: 2GB).
- `buffer_in_memory`: If true, buffer the whole output in memory before writing.

# Returns
- The output filename as a String.
"""
function write_tsvgz(;df::DataFrames.DataFrame, filename::String, force::Bool=false, bufsize::Int=2*1024^3, buffer_in_memory::Bool=false)
    # Enforce .tsv.gz extension
    if !endswith(filename, ".tsv.gz")
        if endswith(filename, ".tsv")
            filename = filename * ".gz"
        else
            filename = filename * ".tsv.gz"
        end
    end

    if isfile(filename) && filesize(filename) > 0 && !force
        @warn "File already exists and is non-empty: $filename. Use force=true to overwrite."
    end

    open(filename, "w") do io
        gzip_stream = CodecZlib.GzipCompressorStream(io, bufsize=bufsize)
        # Use the same bufsize for both the Gzip stream and CSV.write
        CSV.write(gzip_stream, df; delim='\t', bufsize=bufsize, buffer_in_memory=buffer_in_memory)
        close(gzip_stream)
    end
    return filename
end

"""
    read_tsvgz(filename::String; buffer_in_memory::Bool=false, bufsize::Int=2*1024^3) -> DataFrames.DataFrame

Read a DataFrame from a gzipped TSV (.tsv.gz) file with proper decompression and buffering.

# Arguments
- `filename`: Path to the gzipped TSV file (must end with `.tsv.gz`).
- `buffer_in_memory`: If true, decompresses the whole file in memory before parsing (default: false).
- `bufsize`: Buffer size for the decompression stream in bytes (default: 2GB).

# Returns
- The loaded DataFrame.
"""
function read_tsvgz(filename::String; buffer_in_memory::Bool=false, bufsize::Int=2*1024^3)
    # Enforce .tsv.gz extension
    if !endswith(filename, ".tsv.gz")
        error("File must have .tsv.gz extension, got: $filename")
    end
    if !isfile(filename)
        error("File not found: $filename")
    end
    result_df = open(filename, "r") do io
        gzip_stream = CodecZlib.GzipDecompressorStream(io, bufsize=bufsize)
        df = CSV.read(gzip_stream, DataFrames.DataFrame; delim='\t', buffer_in_memory=buffer_in_memory)
        close(gzip_stream)
        return df
    end
    return result_df
end

"""
    dataframe_replace_nothing_with_missing(df::DataFrames.DataFrame) -> DataFrames.DataFrame

Return the DataFrame with all `nothing` values replaced by `missing`.
"""
function dataframe_replace_nothing_with_missing(df::DataFrames.DataFrame)
    for (colname, col) in zip(DataFrames.names(df), DataFrames.eachcol(df))
        # Only process columns that can contain Nothing values
        if Nothing <: eltype(col) && any(isnothing, col)
            df[!, colname] = map(x -> isnothing(x) ? missing : x, col)
        end
    end
    return df
end

"""
    check_matrix_fits_in_memory(bytes_needed::Integer; severity::Symbol=:warn)

Checks whether the specified number of bytes can fit in the computer's memory.

- `bytes_needed`: The number of bytes required (output from `estimate_dense_matrix_memory` or `estimate_sparse_matrix_memory`).
- `severity`: What to do if there is not enough available memory. Can be `:warn` (default) or `:error`.

Returns a named tuple:
    (will_fit_total, will_fit_available, total_memory, free_memory, bytes_needed)
Where:
- `will_fit_total`: `true` if the matrix fits in total system memory.
- `will_fit_available`: `true` if the matrix fits in currently available (free) system memory.
- `total_memory`: Total system RAM in bytes.
- `free_memory`: Currently available system RAM in bytes.
- `bytes_needed`: Bytes requested for the matrix.

If `will_fit_available` is false, either warns or errors depending on `severity`.
"""
function check_matrix_fits_in_memory(bytes_needed::Integer; severity::Symbol=:warn)
    # Get total and available RAM
    total_mem = Sys.total_memory()
    # Sys.free_memory() is available on Linux & macOS, but not in Julia 1.6 on Windows.
    free_mem = try
        Sys.free_memory()
    catch
        # Fallback: If not available, just use total memory as a conservative estimate.
        total_mem
    end

    will_fit_total = bytes_needed ≤ total_mem
    will_fit_available = bytes_needed ≤ free_mem

    
    msg = 
    """
    Requested: $(Mycelia.bytes_human_readable(bytes_needed))
    Free: $(Mycelia.bytes_human_readable(free_mem))
    Total: $(Mycelia.bytes_human_readable(total_mem))"
    """

    if !will_fit_available
        if severity == :error
            error("Matrix will not fit in available memory! $msg")
        elseif severity == :warn
            @warn "Matrix may not fit in available memory! $msg"
        end
    end

    return (
        will_fit_total = will_fit_total,
        will_fit_available = will_fit_available,
        total_memory = total_mem,
        free_memory = free_mem,
        bytes_needed = bytes_needed
    )
end

"""
    estimate_dense_matrix_memory(nrows::Integer, ncols::Integer)
    estimate_dense_matrix_memory(T::DataType, nrows::Integer, ncols::Integer)

Estimate the memory required (in bytes) for a dense matrix.

- If `T` is provided, estimate memory for a matrix with element type `T`.
- If `T` is not provided, defaults to Float64.
"""
function estimate_dense_matrix_memory(args...)
    if length(args) == 2
        nrows, ncols = args
        # T = Int
        T = Float64
    elseif length(args) == 3
        T, nrows, ncols = args
    else
        throw(ArgumentError("Invalid arguments. Provide (nrows, ncols) or (T, nrows, ncols)."))
    end
    elsize = Base.sizeof(T)
    total_bytes = elsize * nrows * ncols
    return total_bytes
end

"""
    estimate_sparse_matrix_memory(nrows::Integer, ncols::Integer; nnz=nothing, density=nothing)
    estimate_sparse_matrix_memory(T::DataType, nrows::Integer, ncols::Integer; nnz=nothing, density=nothing)

Estimate the memory required (in bytes) for a sparse matrix in CSC format.

- If `T` is provided, estimate memory for a matrix with element type `T`.
- If `T` is not provided, defaults to Float64.
- You must specify either `nnz` (number of non-zeros) or `density` (proportion of non-zeros, between 0 and 1).
"""
function estimate_sparse_matrix_memory(args...; nnz::Union{Nothing, Integer}=nothing, density::Union{Nothing, AbstractFloat}=nothing)
    if length(args) == 2
        nrows, ncols = args
        # T = Int
        T = Float64
    elseif length(args) == 3
        T, nrows, ncols = args
    else
        throw(ArgumentError("Invalid arguments. Provide (nrows, ncols) or (T, nrows, ncols)."))
    end

    if nnz !== nothing
        nstored = nnz
    elseif density !== nothing
        nstored = round(Int, density * nrows * ncols)
    else
        throw(ArgumentError("Must specify either nnz or density"))
    end

    val_bytes = Base.sizeof(T) * nstored
    rowind_bytes = Base.sizeof(Int) * nstored
    colptr_bytes = Base.sizeof(Int) * (ncols + 1)

    total_bytes = val_bytes + rowind_bytes + colptr_bytes
    return total_bytes
end

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
        # TODO warn if we find more than one
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
    save_df_jld2(;df::DataFrames.DataFrame, filename::String, key::String="dataframe")

Save a DataFrame to a JLD2 file.

# Arguments
- `df`: The DataFrame to save
- `filename`: Path to the JLD2 file (will add .jld2 extension if not present)
- `key`: The name of the dataset within the JLD2 file (defaults to "dataframe")
"""
function save_df_jld2(;df::DataFrames.DataFrame, filename::String, key::String="dataframe", force::Bool=false)
    # Ensure filename has .jld2 extension
    if !endswith(filename, ".jld2")
        filename = filename * ".jld2"
    end
    
    # Save the dataframe to the JLD2 file
    if !isfile(filename) || force
        JLD2.jldopen(filename, "w") do file
            file[key] = df
        end
    end
    
    return filename
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

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

# seqkit concat $(cat fasta_files.txt) > merged.fasta
# seqtk seq -L $(cat fasta_files.txt) > merged.fasta
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Join fasta files without any regard to record uniqueness.

A cross-platform version of `cat *.fasta > joint.fasta`

See merge_fasta_files

Concatenate multiple FASTA files into a single output file by simple appending.

# Arguments
- `files`: Vector of paths to input FASTA files
- `file`: Path where the concatenated output will be written

# Returns
- Path to the output concatenated file

# Details
Platform-independent implementation of `cat *.fasta > combined.fasta`.
Files are processed sequentially with a progress indicator.
"""
function concatenate_files(;files, file)
    close(open(file, "w"))
    ProgressMeter.@showprogress for f in files
        # stderr=file_path
        run(pipeline(`cat $(f)`, stdout=file, append=true))
    end
    return file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate reaction velocity using Michaelis-Menten enzyme kinetics equation.

# Arguments
- `s::Vector{Float64}`: Substrate concentration(s) [mol/L]
- `p::Vector{Float64}`: Parameters vector where:
    - `p[1]`: vmax - Maximum reaction velocity [mol/(L⋅s)]
    - `p[2]`: km - Michaelis constant [mol/L]

# Returns
- `v::Vector{Float64}`: Reaction velocity [mol/(L⋅s)]

# Description
Implements the Michaelis-Menten equation: v = (vmax * s)/(km + s)

# Reference
Michaelis L., Menten M.L. (1913). Die Kinetik der Invertinwirkung. 
Biochem Z 49:333-369.
"""
# Michaelis–Menten
function calculate_v(s,p)
    vmax = p[1]
    km = p[2]
    v = (vmax .* s) ./ (km .+ s)
    return v
end

struct BandageCompatibilityError <: Exception
    msg::String
end

Base.showerror(io::IO, e::BandageCompatibilityError) = print(io, e.msg)

function _bandage_probe(bin::AbstractString)
    output = IOBuffer()
    cmd = pipeline(ignorestatus(`$(bin) --version`); stdout=output, stderr=output)
    if Sys.islinux()
        return mktempdir() do dir
            Base.withenv(
                "QT_QPA_PLATFORM" => "offscreen",
                "XDG_RUNTIME_DIR" => dir
            ) do
                proc = run(cmd; wait=false)
                wait(proc)
                return success(proc), String(take!(output))
            end
        end
    end

    proc = run(cmd; wait=false)
    wait(proc)
    return success(proc), String(take!(output))
end

function _ensure_bandage_runs(bin::AbstractString)
    ok, output = try
        _bandage_probe(bin)
    catch err
        throw(BandageCompatibilityError("Failed to invoke Bandage binary at $(bin): $(err)"))
    end
    ok && return

    cleaned = strip(output)
    if any(token -> occursin(token, cleaned), ("GLIBC_", "GLIBCXX_", "CXXABI_"))
        throw(BandageCompatibilityError("Bandage binary at $(bin) cannot run on this system (glibc/libstdc++ mismatch).\nOutput:\n$(cleaned)\nProvide a compatible binary via MYCELIA_BANDAGE_CMD or MYCELIA_BANDAGE_URL, or use a system install/module."))
    end

    throw(BandageCompatibilityError("Bandage binary at $(bin) failed to run (--version).\nOutput:\n$(cleaned)"))
end

function _glibc_version_tuple()
    Sys.islinux() || return nothing
    ldd = Sys.which("ldd")
    ldd === nothing && return nothing

    output = IOBuffer()
    proc = run(pipeline(ignorestatus(`$(ldd) --version`); stdout=output, stderr=output); wait=false)
    wait(proc)

    text = String(take!(output))
    lines = split(text, '\n'; keepempty=false)
    isempty(lines) && return nothing

    m = match(r"ldd \(GNU libc\) (\d+)\.(\d+)", strip(lines[1]))
    # m = match(r"ldd \\(GNU libc\\) (\\d+)\\.(\\d+)", strip(lines[1]))
    m === nothing && return nothing

    return (parse(Int, m.captures[1]), parse(Int, m.captures[2]))
end

function _default_bandage_urls()
    if Sys.islinux()
        bandage_ng = "https://github.com/asl/BandageNG/releases/download/v2025.12.1/BandageNG-Linux-e80ad3a.AppImage"
        bandage_legacy = "https://github.com/rrwick/Bandage/releases/download/v0.9.0/Bandage_Ubuntu-x86-64_v0.9.0_AppImage.zip"

        glibc = _glibc_version_tuple()
        if glibc !== nothing && glibc < (2, 29)
            return (bandage_legacy,)
        end

        return (bandage_ng, bandage_legacy)
    end

    if Sys.isapple()
        return ("https://github.com/asl/BandageNG/releases/download/v2025.12.1/BandageNG-macOS-e80ad3a.dmg",)
    end

    error("Bandage download is only automated for Linux and macOS.")
end

function _install_bandage_from_url(url::AbstractString, bandage_executable::AbstractString)
    mktempdir() do tmp
        archive_path = joinpath(tmp, basename(url))
        Downloads.download(url, archive_path)

        lower_url = lowercase(url)
        if endswith(lower_url, ".appimage")
            cp(archive_path, bandage_executable; force=true)
            chmod(bandage_executable, 0o755)
            return
        end

        if endswith(lower_url, ".dmg")
            Sys.isapple() || error("DMG installs are only supported on macOS")
            mountpoint = joinpath(tmp, "mnt")
            mkpath(mountpoint)
            run(`hdiutil attach -nobrowse -mountpoint $(mountpoint) $(archive_path)`)
            candidate = joinpath(mountpoint, "BandageNG.app", "Contents", "MacOS", "BandageNG")
            isfile(candidate) || error("BandageNG executable not found in DMG")
            cp(candidate, bandage_executable; force=true)
            chmod(bandage_executable, 0o755)
            run(`hdiutil detach $(mountpoint)`)
            return
        end

        run(`unzip $(archive_path) -d $(tmp)`)

        best_candidate = nothing
        best_score = 0
        for (root, _, files) in walkdir(tmp)
            for file in files
                path = joinpath(root, file)
                isfile(path) || continue
                name = basename(path)
                lower = lowercase(name)
                score = if name == "BandageNG"
                    2
                elseif name == "Bandage"
                    1
                elseif endswith(lower, ".appimage") && occursin("bandage", lower)
                    3
                else
                    0
                end
                if score > best_score
                    best_score = score
                    best_candidate = path
                end
            end
        end

        best_candidate === nothing && error("Bandage archive did not contain an executable")
        cp(best_candidate, bandage_executable; force=true)
        chmod(bandage_executable, 0o755)
        return
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads and installs Bandage, a bioinformatics visualization tool for genome assembly graphs.

# Arguments
- `outdir`: Target installation directory for the Bandage executable. Defaults to a writable bin directory inside the first depot.

# Returns
- Path to the installed Bandage executable

# Details
- Prefers an existing `BandageNG` or `Bandage` on PATH
- Accepts overrides via `ENV["MYCELIA_BANDAGE_CMD"]` (absolute path) or `ENV["MYCELIA_BANDAGE_URL"]` (archive URL)
- Defaults to BandageNG (and falls back to legacy Bandage on Linux if needed)
- Installs into a user-writable directory without requiring sudo or apt packages
- Skips download if Bandage is already installed at target location
- Validates the selected binary with `--version` and throws `BandageCompatibilityError` if it cannot run

# Dependencies
Requires system commands: unzip (for zip archives), hdiutil on macOS (for DMG)
"""
function download_bandage(outdir=joinpath(first(DEPOT_PATH), "bin"))
    cmd_override = get(ENV, "MYCELIA_BANDAGE_CMD", nothing)
    if cmd_override !== nothing
        _ensure_bandage_runs(cmd_override)
        return cmd_override
    end

    # Prefer a pre-installed BandageNG or Bandage on PATH
    existing = Sys.which("BandageNG")
    existing === nothing && (existing = Sys.which("Bandage"))
    if existing !== nothing
        _ensure_bandage_runs(existing)
        return existing
    end

    mkpath(outdir)
    bandage_executable = joinpath(outdir, "Bandage")
    if isfile(bandage_executable)
        try
            _ensure_bandage_runs(bandage_executable)
            return bandage_executable
        catch e
            (e isa BandageCompatibilityError) || rethrow()
            get(ENV, "MYCELIA_BANDAGE_URL", nothing) === nothing || rethrow()
        end
    end

    urls = if (url_override = get(ENV, "MYCELIA_BANDAGE_URL", nothing)) === nothing
        _default_bandage_urls()
    else
        (url_override,)
    end

    last_error = nothing
    for url in urls
        try
            _install_bandage_from_url(url, bandage_executable)
            _ensure_bandage_runs(bandage_executable)
            return bandage_executable
        catch e
            last_error = e
            if get(ENV, "MYCELIA_BANDAGE_URL", nothing) !== nothing
                rethrow()
            end
        end
    end

    last_error === nothing && error("Bandage download failed.")
    throw(last_error)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a copy of a file in a temporary directory while preserving the original filename.

# Arguments
- `file_path::String`: Path to the source file to be copied

# Returns
- `String`: Path to the newly created temporary file
"""
function copy_to_tempdir(file_path::String)
    # Create a temporary directory
    temp_dir = mktempdir()

    # Get the file name from the original file path
    file_name = basename(file_path)

    # Create the new file path within the temporary directory
    temp_file_path = joinpath(temp_dir, file_name)

    # Copy the original file to the new path
    cp(file_path, temp_file_path)

    # Return the path of the temporary file
    return temp_file_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns the current date and time as a normalized string with all non-word characters removed.

The output format is based on ISO datetime (YYYYMMDDThhmmss) but strips any special characters
like hyphens, colons or dots.
"""
function normalized_current_datetime()
    return replace(Dates.format(Dates.now(), Dates.ISODateTimeFormat), r"[^\w]" => "")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns the current date as a normalized string with all non-word characters removed.

The output format is based on ISO datetime (YYYYMMDD) but strips any special characters
like hyphens, colons or dots.
"""
function normalized_current_date()
    return replace(Dates.format(Dates.today(), Dates.ISODateFormat), r"[^\w]" => "")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns the current git commit hash of the repository.

# Arguments
- `short::Bool=false`: If true, returns abbreviated 8-character hash

# Returns
A string containing the git commit hash (full 40 characters by default)
"""
function githash(;short=false)
    git_hash = rstrip(read(`git rev-parse HEAD`, String))
    if short
        git_hash = git_hash[1:8]
    end
    return git_hash
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identify all columns that have only missing or empty values

Returns as a bit array

See also: drop_empty_columns, drop_empty_columns!
"""
function find_nonempty_columns(df)
    non_empty_columns = Bool[]
    for col in DataFrames.eachcol(df)
        # Skip columns that are entirely Missing or Nothing type
        if eltype(col) == Missing || eltype(col) == Nothing
            push!(non_empty_columns, false)
            continue
        end
        
        # Check if any value in the column is non-empty
        has_nonempty = any(col) do v
            if ismissing(v) || isnothing(v)
                return false
            end
            # For Date types, consider them non-empty
            if isa(v, Dates.Date)
                return true
            end
            # For other types, check if they're not empty
            return !isempty(v)
        end
        
        push!(non_empty_columns, has_nonempty)
    end
    return non_empty_columns
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identify all columns that have only missing or empty values, and remove those columns from the dataframe.

Returns a modified copy of the dataframe.

See also: drop_empty_columns!
"""
function drop_empty_columns(df::DataFrames.AbstractDataFrame)
    # Filter the DataFrame columns by checking if not all values in the column are missing or empty
    non_empty_columns = find_nonempty_columns(df)
    filtered_df = df[:, non_empty_columns]
    return filtered_df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identify all columns that have only missing or empty values, and remove those columns from the dataframe *in-place*.

Returns a modified version of the original dataframe. 

See also: drop_empty_columns
"""
function drop_empty_columns!(df::DataFrames.AbstractDataFrame)
    # Filter the DataFrame columns by checking if not all values in the column are missing or empty
    non_empty_columns = find_nonempty_columns(df)
    # df = df[!, non_empty_columns]
    DataFrames.select!(df, non_empty_columns)
    return df
end

# # Need to add hashdeep & logging
# function tarchive(;directory, tarchive=directory * ".tar.gz")
#     run(`tar --create --gzip --verbose --file=$(tarchive) $(directory)`)
#     return tarchive
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Ensures the hashdeep utility is installed on the system.

Checks if hashdeep is available in PATH and attempts to install it via apt package manager
if not found. Will try with sudo privileges first, then without sudo if that fails.

# Details
- Checks PATH for existing hashdeep executable
- Attempts installation using apt package manager
- Requires a Debian-based Linux distribution

# Returns
- Nothing, but prints status messages during execution
"""
function install_hashdeep()
    if Sys.which("hashdeep") !== nothing
        println("hashdeep executable found")
    else
        println("hashdeep executable not found in PATH, installing")
        try
            run(`sudo apt install hashdeep -y`)
        catch
            run(`apt install hashdeep -y`)
        end
    end
end

# Need to add hashdeep & logging
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a gzipped tar archive of the specified directory along with verification files.

# Arguments
- `directory`: Source directory path to archive
- `tarchive`: Optional output archive path (defaults to directory name with .tar.gz extension)

# Generated Files
- `{tarchive}`: The compressed tar archive
- `{tarchive}.log`: Contents listing of the archive
- `{tarchive}.hashdeep.dfxml`: Cryptographic hashes (MD5, SHA1, SHA256) of the archive

# Returns
- Path to the created tar archive file
"""
function create_tarchive(;directory, tarchive=directory * ".tar.gz")
    directory = normpath(directory)
    tarchive = normpath(tarchive)
    output_dir = mkpath(dirname(tarchive))
    
    install_hashdeep()


    working_dir, source = splitdir(directory)
    if isempty(source)
        source = working_dir
        working_dir = pwd()
    end
    if isempty(working_dir)
        working_dir = pwd()
    end
    println("output_dir: $output_dir\n" *
            "working_dir: $working_dir\n" *
            "source: $source\n" *
            "target_file: $tarchive")
    log_file = tarchive * ".log"
    hashdeep_file = tarchive * ".hashdeep.dfxml"
    
    if !isfile(tarchive)
        run(`tar --create --gzip --verbose --file=$(tarchive) $(directory)`)
    end
    if !isfile(log_file)
        run(pipeline(`tar -tvf $(tarchive)`,log_file))
    end
    if !isfile(hashdeep_file)
        run(pipeline(`hashdeep -c md5,sha1,sha256 -b -d $(tarchive)`,hashdeep_file))
    end
    return tarchive
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract contents of a gzipped tar archive file to a specified directory.

# Arguments
- `tarchive::AbstractString`: Path to the .tar.gz file to extract
- `directory::AbstractString=dirname(tarchive)`: Target directory for extraction (defaults to the archive's directory)

# Returns
- `AbstractString`: Path to the directory where contents were extracted
"""
function tar_extract(;tarchive, directory=dirname(tarchive))
    run(`tar --extract --gzip --verbose --file=$(tarchive) --directory=$(directory)`)
    return directory
end

"""
    gzip_file(
        infile::AbstractString;
        outfile::Union{Nothing,AbstractString}=nothing,
        threads::Integer=get_default_threads(),
        force::Bool=false,
        keep_input::Bool=false
    ) -> String

Compress a file to gzip format using an external compressor (`pigz` preferred, else `gzip`).

If `outfile` is `nothing`, defaults to `infile * ".gz"` (unless `infile` already ends with `.gz`,
in which case it is returned unchanged). Uses an atomic write (`outfile * ".tmp"` then `mv`).

This avoids writing directly to CodecZlib compressor streams, which can hit buffering issues
for large streaming writes.
"""
function gzip_file(
    infile::AbstractString;
    outfile::Union{Nothing,AbstractString}=nothing,
    threads::Integer=get_default_threads(),
    force::Bool=false,
    keep_input::Bool=false
)::String
    @assert isfile(infile) "Input file not found: $infile"

    out = if outfile === nothing
        endswith(infile, ".gz") ? infile : infile * ".gz"
    else
        String(outfile)
    end
    endswith(out, ".gz") || error("gzip_file outfile must end with .gz: $out")
    out == infile && return out

    mkpath(dirname(out))
    if isfile(out) && !force
        error("Refusing to overwrite existing file: $out (use force=true)")
    end

    tmp = out * ".tmp"
    isfile(tmp) && rm(tmp; force=true)

    pigz = Sys.which("pigz")
    gzip = Sys.which("gzip")
    compressor = if pigz !== nothing
        Cmd([pigz, "--processes", string(max(1, threads)), "-c", infile])
    elseif gzip !== nothing
        Cmd([gzip, "-c", infile])
    else
        Mycelia.add_bioconda_env("pigz")
        Cmd([Mycelia.CONDA_RUNNER, "run", "--live-stream", "-n", "pigz", "pigz", "--processes", string(max(1, threads)), "-c", infile])
    end

    open(tmp, "w") do out_io
        run(pipeline(compressor, stdout=out_io))
    end
    mv(tmp, out; force=true)

    if !keep_input
        rm(infile; force=true)
    end

    return out
end

"""
    get_arg_max() -> Union{Int,Nothing}

Return the system `ARG_MAX` (maximum total size of argv+env for `execve`) if it can be determined.
Uses `getconf ARG_MAX` when available; returns `nothing` if unavailable.
"""
function get_arg_max()::Union{Int,Nothing}
    getconf = Sys.which("getconf")
    getconf === nothing && return nothing
    try
        return parse(Int, strip(readchomp(`$getconf ARG_MAX`)))
    catch
        return nothing
    end
end

"""
    estimate_argv_bytes(args::AbstractVector{<:AbstractString}) -> Int

Estimate bytes consumed by argv strings for `execve` (sum of argument lengths plus NULs).
Does not include environment size.
"""
function estimate_argv_bytes(args::AbstractVector{<:AbstractString})::Int
    return sum(length, args) + length(args)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract biosample and barcode information from a PacBio XML metadata file.

# Arguments
- `xml`: Path to PacBio XML metadata file

# Returns
DataFrame with two columns:
- `BioSampleName`: Name of the biological sample
- `BarcodeName`: Associated DNA barcode identifier
"""
function extract_pacbiosample_information(xml)
    xml_dict = XMLDict.parse_xml(read(xml, String))
    wellsample = xml_dict["ExperimentContainer"]["Runs"]["Run"]["Outputs"]["SubreadSets"]["SubreadSet"]["DataSetMetadata"]["Collections"]["CollectionMetadata"]["WellSample"]

    # Initialize empty arrays to store the data
    biosample_names = []
    barcode_names = []
    
    if haskey(wellsample, "BioSamples")
        # display(wellsample)
        for bs in wellsample["BioSamples"]["BioSample"]
            push!(biosample_names, bs[:Name])
            push!(barcode_names, bs["DNABarcodes"]["DNABarcode"][:Name])
        end
    end

    # Create the DataFrame
    df = DataFrames.DataFrame(BioSampleName=biosample_names, BarcodeName=barcode_names)
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalize a dictionary of counts into a probability distribution where values sum to 1.0.

# Arguments
- `countmap::Dict`: Dictionary mapping keys to count values

# Returns
- `Dict`: New dictionary with same keys but values normalized by total sum
"""
function normalize_countmap(countmap)
    sum_total = sum(values(countmap))
    return Dict(k => v/sum_total for (k, v) in countmap)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Copy a file to a new location with a unique identifier prepended to the filename.

# Arguments
- `infile::AbstractString`: Path to the source file to copy
- `out_directory::AbstractString`: Destination directory for the copied file
- `unique_identifier::AbstractString`: String to prepend to the filename
- `force::Bool=true`: If true, overwrite existing files

# Returns
- `String`: Path to the newly created file
"""
function copy_with_unique_identifier(infile, out_directory, unique_identifier; force=true)
    outfile = joinpath(out_directory, unique_identifier * "." * basename(infile))
    cp(infile, outfile, force=force)
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Gets the size of a file and returns it in a human-readable format.

# Arguments
- `f`: The path to the file, either as a `String` or an `AbstractString`.

# Returns
A string representing the file size in a human-readable format (e.g., "3.40 MB").

# Details
This function internally uses `filesize(f)` to get the file size in bytes, then leverages `Base.format_bytes` to convert it into a human-readable format with appropriate units (KB, MB, GB, etc.).

# Examples
```julia
julia> filesize_human_readable("my_image.jpg")
"2.15 MB"
```
See Also
* filesize: Gets the size of a file in bytes.
* Base.format_bytes: Converts a byte count into a human-readable string. 
"""
function filesize_human_readable(f)
    @assert isfile(f) "File does not exist: $f"
    return Base.format_bytes(filesize(f))
end

function bytes_human_readable(num::Real)
    return Base.format_bytes(num)
    # units = ["B", "KB", "MB", "GB", "TB", "PB", "EB"]
    # i = 1
    # while num ≥ 1024 && i < length(units)
    #     num /= 1024
    #     i += 1
    # end
    # return @sprintf("%.2f %s", num, units[i])
end

function scientific_notation(num::Number; precision::Int=2)
    if precision < 0
        error("Number of decimal places (precision) must be non-negative.")
    end
    return Printf.format(Printf.Format("%.$(precision)e"), num)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Find the longest common prefix between two filenames.

# Arguments
- `filename1::String`: First filename to compare
- `filename2::String`: Second filename to compare

# Keywords
- `strip_trailing_delimiters::Bool=true`: If true, removes trailing dots, hyphens, and underscores from the result

# Returns
- `String`: The longest common prefix found between the filenames
"""
function find_matching_prefix(filename1::String, filename2::String; strip_trailing_delimiters=true)
    min_length = min(length(filename1), length(filename2))
    matching_prefix = ""
    
    for i in 1:min_length
        if filename1[i] == filename2[i]
            matching_prefix *= filename1[i]
        else
            break
        end
    end
    if strip_trailing_delimiters
        matching_prefix = replace(matching_prefix, r"[\.\-_]+$" => "")
    end
    
    return matching_prefix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

    parse_jsonl(filepath::String) -> Vector{Dict{String,Any}}

Validate and parse a JSON Lines file (either .ndjson/.jsonl, optionally gzipped)
into a vector of dictionaries, reporting progress in bytes processed.

Validations performed:
  • Extension must be one of: .jsonl, .ndjson, .jsonl.gz, .ndjson.gz
  • File must exist
  • File size must be non-zero

Progress meter shows bytes read from the underlying file (compressed bytes
for .gz). No second full pass is needed.
"""
function parse_jsonl(filepath::String)::Vector{Dict{String,Any}}
    # --- Validate inputs ---
    fname = lowercase(filepath)
    valid_exts = (".jsonl", ".ndjson", ".jsonl.gz", ".ndjson.gz")
    if !any(endswith(fname, ext) for ext in valid_exts)
        error("parse_jsonl: unsupported extension. “$(filepath)” must end in $(join(valid_exts, ", "))")
    end
    if !isfile(filepath)
        error("parse_jsonl: file not found: $filepath")
    end
    file_stat = Base.Filesystem.stat(filepath)
    if file_stat.size == 0
        error("parse_jsonl: file is empty: $filepath")
    end

    # --- Open raw and (optionally) decompress ---
    raw_io = open(filepath, "r")
    io     = endswith(fname, ".gz") ? CodecZlib.GzipDecompressorStream(raw_io) : raw_io

    # --- Set up a byte‐based progress meter ---
    total_bytes = file_stat.size
    p = ProgressMeter.Progress(total_bytes, "Parsing JSONL (bytes): ")
    last_pos = Int64(0)

    # --- Read & parse lines ---
    results = Vector{Dict{String,Any}}()
    for line in eachline(io)
        s = strip(line)
        if !isempty(s)
            push!(results, JSON.parse(s))
        end
        # track compressed‐byte progress even for gzipped files
        curr_pos = position(raw_io)
        delta    = curr_pos - last_pos
        last_pos = curr_pos
        ProgressMeter.next!(p; step = delta)
    end

    close(io)
    return results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

    jsonl_to_dataframe(filepath::String) -> DataFrame

Parse a JSONL (or gzipped JSONL) file and return a DataFrame.
Internally calls `parse_jsonl` for validation and parsing.
Ensures that all rows have the same set of keys by inserting `missing`
for any absent field before constructing the DataFrame.
"""
function jsonl_to_dataframe(filepath::String)::DataFrames.DataFrame
    rows = parse_jsonl(filepath)
    # If no rows, return empty DataFrame
    if isempty(rows)
        return DataFrames.DataFrame()
    end

    # Collect the union of all keys across rows
    all_keys = reduce(union, map(keys, rows))
    # Fill missing columns in each row with `missing`
    ProgressMeter.@showprogress desc="Filling missing values" for row in rows
        for k in setdiff(all_keys, keys(row))
            row[k] = missing
        end
    end

    return DataFrames.DataFrame(rows)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the default number of threads to use for parallel operations.

Checks environment hints in the following order:
1. `JULIA_NUM_THREADS` (explicit user intent)
2. SLURM allocation (`SLURM_JOB_CPUS_PER_NODE` preferred, or `SLURM_CPUS_PER_TASK * SLURM_TASKS_PER_NODE`)
3. PBS allocation (`PBS_NCPUS`)
4. Falls back to half of detected system cores (capped at 16)
5. Uses `DEFAULT_THREADS` (4) if hardware cannot be detected

# Returns
- `Int`: Number of threads to use

# Example
```julia
threads = get_default_threads()  # Returns allocated threads or 4
```
"""
function get_default_threads()::Int
    parse_env_int(var) = begin
        val = get(ENV, var, nothing)
        if val === nothing
            return nothing
        end
        str = String(val)
        # SLURM values can include decorations like "24(x2)"; grab the leading integer.
        m = match(r"^\s*(\d+)", str)
        parsed = m === nothing ? tryparse(Int, str) : tryparse(Int, m.captures[1])
        parsed !== nothing && parsed > 0 ? parsed : nothing
    end

    cpu_threads = Sys.CPU_THREADS > 0 ? Sys.CPU_THREADS : nothing

    # Respect explicit Julia threading request.
    julia_hint = parse_env_int("JULIA_NUM_THREADS")
    if julia_hint !== nothing
        return cpu_threads === nothing ? max(1, julia_hint) : clamp(julia_hint, 1, cpu_threads)
    end

    # SLURM hints: prefer job-level declaration, otherwise derive from per-task × tasks-per-node.
    slurm_job = parse_env_int("SLURM_JOB_CPUS_PER_NODE")
    slurm_cpt = parse_env_int("SLURM_CPUS_PER_TASK")
    slurm_tpn = parse_env_int("SLURM_TASKS_PER_NODE")
    slurm_derived = slurm_cpt !== nothing && slurm_tpn !== nothing ? slurm_cpt * slurm_tpn : nothing
    slurm_on_node = parse_env_int("SLURM_CPUS_ON_NODE")
    slurm_hint =
        slurm_job !== nothing ? slurm_job :
        slurm_derived !== nothing ? slurm_derived :
        slurm_on_node

    if slurm_hint !== nothing
        return cpu_threads === nothing ? slurm_hint : min(slurm_hint, cpu_threads)
    end

    # PBS/Torque hint as a fallback for other schedulers.
    pbs_hint = parse_env_int("PBS_NCPUS")
    if pbs_hint !== nothing
        return cpu_threads === nothing ? pbs_hint : min(pbs_hint, cpu_threads)
    end

    # No scheduler hints: default to half of hardware threads (rounded up), capped at 16.
    if cpu_threads !== nothing
        return min(cld(cpu_threads, 2), 16)
    end
    return DEFAULT_THREADS
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
struct SystemOverview
    system_threads::Int
    julia_threads::Int
    default_threads::Int
    total_memory::Int
    available_memory::Int
    occupied_memory::Int
    memory_occupied_percent::Float64
    memory_running_low::Bool
    total_storage::Int
    available_storage::Int
    occupied_storage::Int
    storage_occupied_percent::Float64
    storage_running_low::Bool
    path::String
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Collect system resource information (raw values) with human-readable display.

Notes:
- `julia_threads` reflects the threads the current session actually started with (`Threads.nthreads()`).
- `default_threads` reflects the configured default based on environment hints (`get_default_threads`; e.g., `SLURM_CPUS_PER_TASK`, `PBS_NCPUS`, `OMP_NUM_THREADS`, `JULIA_NUM_THREADS`, or the conservative fallback) and may differ if the session was launched with different settings.
"""
function system_overview(; path=pwd(), memory_low_threshold=0.90, storage_low_threshold=0.90)
    if !(0.0 <= memory_low_threshold <= 1.0)
        error("memory_low_threshold must be between 0.0 and 1.0 (got $memory_low_threshold)")
    end
    if !(0.0 <= storage_low_threshold <= 1.0)
        error("storage_low_threshold must be between 0.0 and 1.0 (got $storage_low_threshold)")
    end

    total_memory = Sys.total_memory()
    available_memory = Sys.free_memory()
    occupied_memory = total_memory - available_memory
    memory_occupied_percent = total_memory == 0 ? 0.0 : occupied_memory / total_memory

    disk_stats = Base.diskstat(path)
    total_storage = disk_stats.total
    available_storage = disk_stats.available
    occupied_storage = disk_stats.used
    storage_occupied_percent = total_storage == 0 ? 0.0 : occupied_storage / total_storage

    system_threads = Sys.CPU_THREADS
    julia_threads = Threads.nthreads()
    default_threads = get_default_threads()

    if julia_threads > system_threads
        @error "The Julia kernel has $julia_threads allocated but only $system_threads are available on the system"
    elseif julia_threads < system_threads
        @warn "The Julia kernel only has $julia_threads allocated but $system_threads are available on the system"
    else
        @info "The Julia kernel is using all available threads (n=$system_threads) on the system"
    end

    if julia_threads != default_threads
        @info "Julia threads ($julia_threads) differ from configured default ($default_threads)"
    end

    memory_running_low = memory_occupied_percent >= memory_low_threshold
    storage_running_low = storage_occupied_percent >= storage_low_threshold

    if memory_running_low
        @warn "Memory usage high: $(bytes_human_readable(occupied_memory)) / $(bytes_human_readable(total_memory)) ($(round(memory_occupied_percent * 100; digits=1))% used)"
    else
        @info "Memory usage: $(bytes_human_readable(occupied_memory)) / $(bytes_human_readable(total_memory)) ($(round(memory_occupied_percent * 100; digits=1))% used)"
    end

    if storage_running_low
        @warn "Storage usage high: $(bytes_human_readable(occupied_storage)) / $(bytes_human_readable(total_storage)) ($(round(storage_occupied_percent * 100; digits=1))% used)"
    else
        @info "Storage usage: $(bytes_human_readable(occupied_storage)) / $(bytes_human_readable(total_storage)) ($(round(storage_occupied_percent * 100; digits=1))% used)"
    end

    return SystemOverview(
        system_threads,
        julia_threads,
        default_threads,
        total_memory,
        available_memory,
        occupied_memory,
        memory_occupied_percent,
        memory_running_low,
        total_storage,
        available_storage,
        occupied_storage,
        storage_occupied_percent,
        storage_running_low,
        String(path),
    )
end

function Base.show(io::IO, ::MIME"text/plain", overview::SystemOverview)
    memory_pct = round(overview.memory_occupied_percent * 100; digits=1)
    storage_pct = round(overview.storage_occupied_percent * 100; digits=1)

    memory_flag = overview.memory_running_low ? " (LOW)" : ""
    storage_flag = overview.storage_running_low ? " (LOW)" : ""

    println(io, "SystemOverview:")
    println(io, "  Threads: julia=$(overview.julia_threads), system=$(overview.system_threads), default=$(overview.default_threads)")
    println(io, "  Memory: total=$(bytes_human_readable(overview.total_memory)), available=$(bytes_human_readable(overview.available_memory)), occupied=$(bytes_human_readable(overview.occupied_memory)) ($memory_pct% used)$memory_flag")
    println(io, "  Storage: total=$(bytes_human_readable(overview.total_storage)), available=$(bytes_human_readable(overview.available_storage)), occupied=$(bytes_human_readable(overview.occupied_storage)) ($storage_pct% used)$storage_flag")
    println(io, "  Path: $(overview.path)")
end

Base.show(io::IO, overview::SystemOverview) = Base.show(io, MIME"text/plain"(), overview)

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Vertically concatenate DataFrames with different column structures by automatically handling missing values.

# Arguments
- `dfs`: Variable number of DataFrames to concatenate vertically

# Returns
- `DataFrame`: Combined DataFrame containing all rows and columns from input DataFrames, 
  with `missing` values where columns didn't exist in original DataFrames
"""
function vcat_with_missing(dfs::Vararg{DataFrames.AbstractDataFrame})
    # Get all unique column names across all DataFrames
    all_columns = unique(reduce(vcat, [names(df) for df in dfs]))

    # Add missing columns to each DataFrame and fill with missing
    for df in dfs
        for col in all_columns
            if !(col in names(df))
                df[!, col] .= missing
            end
        end
    end

    # Now you can safely vcat the DataFrames
    return vcat(dfs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute a single SHA256 hash from multiple SHA256 hashes.

Takes a vector of hex-encoded SHA256 hashes and produces a new SHA256 hash by:
1. Sorting the input hashes lexicographically
2. Concatenating them in sorted order
3. Computing a new SHA256 hash over the concatenated data

# Arguments
- `vector_of_sha256s`: Vector of hex-encoded SHA256 hash strings

# Returns
- A hex-encoded string representing the computed meta-hash
"""
function metasha256(vector_of_sha256s::Vector{<:AbstractString})
    ctx = SHA.SHA2_256_CTX()
    for sha_hash in sort(vector_of_sha256s)
        SHA.update!(ctx, collect(codeunits(sha_hash)))
    end
    return SHA.bytes2hex(SHA.digest!(ctx))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract the base file extension from a filename, handling compressed files.

For regular files, returns the last extension. For gzipped files, returns the extension
before .gz.
"""
function get_base_extension(filename::String)
  parts = split(basename(filename), "."; limit=3)  # Limit to 3 to handle 2-part extensions
  extension = parts[end]  # Get the last part
  
  if extension == "gz" && length(parts) > 2  # Check for .gz and more parts
    extension = parts[end - 1]  # Get the part before .gz
  end
  
  return "." * extension
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Finds contiguous ranges of `true` values in a boolean vector.

# Arguments
- `bool_vec::AbstractVector{Bool}`: Input boolean vector to analyze
- `min_length=1`: Minimum length requirement for a range to be included

# Returns
Vector of tuples `(start, end)` where each tuple represents the indices of a
contiguous range of `true` values meeting the minimum length requirement.
"""
function find_true_ranges(bool_vec::AbstractVector{Bool}; min_length=1)
    indices = findall(bool_vec)  # Get indices of true values
    if isempty(indices)
    return []  # Handle the case of no true values
    end
    diffs = diff(indices)
    breakpoints = findall(>(1), diffs)  # Find where the difference is greater than 1
    starts = [first(indices); indices[breakpoints .+ 1]]
    ends = [indices[breakpoints]; last(indices)]
    # true_ranges_table = DataFrames.DataFrame(starts = starts, ends = ends, lengths = ends .- starts)
    return collect(Iterators.filter(x -> (x[2] - x[1]) >= min_length, zip(starts, ends)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Sample `n` equally spaced elements from `vector`.

# Arguments
- `vector`: Input vector to sample from
- `n`: Number of samples to return (must be positive)

# Returns
A vector containing `n` equally spaced elements from the input vector.
"""
function equally_spaced_samples(vector, n)
    indices = round.(Int, range(1, length(vector), length=n))
    return vector[indices]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute a centered moving average over a vector using a sliding window.

# Arguments
- `data::AbstractVector{T}`: Input vector to be averaged
- `window_size::Int`: Size of the sliding window (odd number recommended)

# Returns
- `Vector{Float64}`: Vector of same length as input containing moving averages

# Details
- For points near the edges, the window is truncated to available data
- Window is centered on each point, using floor(window_size/2) points on each side
- Result type is always Float64 regardless of input type T
"""
function rolling_centered_avg(data::AbstractVector{T}; window_size::Int) where T
    half_window = Int(floor(window_size / 2))
    result = Vector{Float64}(undef, length(data))
    for i in eachindex(data)
        start_idx = max(1, i - half_window)
        end_idx = min(length(data), i + half_window)
        result[i] = Statistics.mean(data[start_idx:end_idx])
    end
    return result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Select one random row from each group in a grouped DataFrame.

# Arguments
- `gdf::GroupedDataFrame`: A grouped DataFrame created using `groupby`

# Returns
- `DataFrame`: A new DataFrame containing exactly one randomly sampled row from each group
"""
function rand_of_each_group(gdf::DataFrames.GroupedDataFrame{DataFrames.DataFrame})
    result = DataFrames.combine(gdf) do sdf
        sdf[StatsBase.sample(1:DataFrames.nrow(sdf), 1), :]
    end
    return result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a random symmetric distance matrix of size n×n with zeros on the diagonal.

# Arguments
- `n`: Positive integer specifying the matrix dimensions

# Returns
- A symmetric n×n matrix with random values in [0,1), zeros on the diagonal

# Details
- The matrix is symmetric, meaning M[i,j] = M[j,i]
- Diagonal elements M[i,i] are set to 0.0
- Off-diagonal elements are uniformly distributed random values
"""
function random_symmetric_distance_matrix(n)
  # Generate a random matrix
  matrix = rand(n, n)

  # Make the matrix symmetric
  matrix = (matrix + matrix') / 2

  # Ensure the diagonal is zero
  for i in 1:n
    matrix[i, i] = 0.0
  end

  return matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return a DataFrame containing the first row from each group in a GroupedDataFrame.

# Arguments
- `gdf::GroupedDataFrame`: A grouped DataFrame created using `groupby`

# Returns
- `DataFrame`: A new DataFrame containing first row from each group
"""
function first_of_each_group(gdf::DataFrames.GroupedDataFrame{DataFrames.DataFrame})
    return DataFrames.combine(gdf, first)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a type to its string representation, with special handling for Kmer types.

# Arguments
- `T`: The type to convert to string

# Returns
- String representation of the type
  - For Kmer types: Returns "Kmers.DNAKmer{K}" where K is the kmer length
  - For other types: Returns the standard string representation
"""
function type_to_string(T)
    if T <: Kmers.Kmer
        return "Kmers.DNAKmer{$(T.parameters[2])}"
    else
        return string(T)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Converts an AbstractString type to its string representation.

# Arguments
- `T::AbstractString`: The string type to convert

# Returns
A string representation of the input type
"""
function type_to_string(T::AbstractString)
    return string(T)
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Convert between different graph file formats.

# # Arguments
# - `args`: Dictionary with required keys:
#     - `"in"`: Input filepath (supported: .jld2, .gfa, .neo4j)
#     - `"out"`: Output filepath (supported: .jld2, .gfa, .neo4j)

# # Details
# Performs format conversion based on file extensions. For non-JLD2 to non-JLD2
# conversions, uses JLD2 as an intermediate format.
# """
# function convert(args)
#     @show args
#     in_type = missing
#     out_type = missing
#     if occursin(r"\.jld2$", args["in"])
#         in_type = :jld2
#     elseif occursin(r"\.gfa$", args["in"])
#         in_type = :gfa
#     elseif occursin(r"\.neo4j$", args["in"])
#         in_type = :neo4j
#     end

#     if occursin(r"\.jld2$", args["out"])
#         out_type = :jld2
#     elseif occursin(r"\.gfa$", args["out"])
#         out_type = :gfa
#     elseif occursin(r"\.neo4j$", args["out"])
#         out_type = :neo4j
#     end
    
#     if ismissing(in_type) || ismissing(out_type)
#         error("unable to determine in and out types")
#     end
    
#     if (in_type == :jld2) && (out_type == :jld2)
#         # done
#     elseif (in_type == :jld2) && (out_type != :jld2)
#         # convert out of jld2
#         if out_type == :gfa
#             loaded = FileIO.load(args["in"])
#             @assert haskey(loaded, "graph")
#             graph = loaded["graph"]
#             graph_to_gfa(graph, args["out"])
#         end
#     elseif (in_type != :jld2) && (out_type == :jld2)
#         # convert into jld2
#     else
#         # need to convert from input to jld2 as an intermediate first
#     end
# end

# Save the distance matrix to a JLD2 file
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Saves a matrix to a JLD2 file format.

# Arguments
- `matrix`: The matrix to be saved
- `filename`: String path where the file should be saved

# Returns
- The filename string that was used to save the matrix
"""
function save_matrix_jld2(;matrix, filename)
    if !isfile(filename) || (filesize(filename) == 0)
        JLD2.@save filename matrix
    else
        @warn "$(filename) already exists and is non-empty, skipping..."
    end
    return filename
end

# Load the distance matrix from a JLD2 file
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Loads a matrix from a JLD2 file.

# Arguments
- `filename::String`: Path to the JLD2 file containing the matrix under the key "matrix"

# Returns
- `Matrix`: The loaded matrix data
"""
function load_matrix_jld2(filename)
    return JLD2.load(filename, "matrix")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load data stored in a JLD2 file format.

# Arguments
- `filename::String`: Path to the JLD2 file to load

# Returns
- `Dict`: Dictionary containing the loaded data structures
"""
function load_jld2(filename)
    return JLD2.load(filename)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate and display frequency counts for all columns in a DataFrame.

# Arguments
- `table::DataFrame`: Input DataFrame to analyze

# Details
Iterates through each column in the DataFrame and displays:
1. The column name
2. A Dict mapping unique values to their frequencies using StatsBase.countmap
"""
function countmap_columns(table)
    for n in names(table)
        display(n)
        display(StatsBase.countmap(table[!, n]))
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Recursively include all files matching a pattern in a directory and its subdirectories.

# Arguments
- `dir::AbstractString`: Directory path to search recursively
- `pattern::Regex=r"\\.jl\$"`: Regular expression pattern to match files (defaults to .jl files)

# Details
Files are processed in sorted order within each directory. This is useful for 
loading test files, examples, or other Julia modules in a predictable order.

# Examples
```julia
# Include all Julia files in a directory tree
include_all_files("test/modules")

# Include all text files
include_all_files("docs", r"\\.txt\$")
```
"""
function include_all_files(dir::AbstractString; pattern::Regex=r"\.jl$")
    for (root, dirs, files) in walkdir(dir)
        for file in sort(files)
            if occursin(pattern, file)
                include(joinpath(root, file))
            end
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the theoretical k-mer space size for a given k-mer length and alphabet size.

# Arguments
- `k::Integer`: K-mer length
- `alphabet_size::Integer=4`: Size of the alphabet (defaults to 4 for DNA: A,C,G,T)

# Returns
- `Integer`: Total number of possible k-mers (alphabet_size^k)

# Details
For DNA sequences (alphabet_size=4), this computes 4^k. Useful for:
- Memory estimation for k-mer analysis
- Parameter validation and selection
- Understanding computational complexity

# Examples
```julia
# DNA 3-mers: 4^3 = 64 possible k-mers
kmer_space_size(3)

# Protein 5-mers: 20^5 = 3,200,000 possible k-mers  
kmer_space_size(5, 20)
```
"""
function kmer_space_size(k::Integer, alphabet_size::Integer=4)
    return alphabet_size^k
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Clean up a directory by calculating its size and file count, then removing it.

# Arguments
- `directory::AbstractString`: Path to the directory to clean up
- `verbose::Bool=true`: Whether to report cleanup results (default: true)
- `force::Bool=false`: Whether to proceed without confirmation for large directories

# Returns
- Named tuple with fields:
  - `existed`: Whether the directory existed before cleanup
  - `files_deleted`: Number of files that were deleted
  - `bytes_freed`: Total bytes freed up
  - `human_readable_size`: Human-readable representation of bytes freed

# Details
This function will:
1. Check if the directory exists and is non-empty
2. Calculate the total size and number of files recursively
3. Remove the directory and all its contents
4. Report the cleanup results unless verbose=false

For directories larger than 1GB or containing more than 10,000 files,
confirmation is required unless force=true.

# Examples
```julia
# Clean up a temporary directory with reporting
result = cleanup_directory("/tmp/myapp_temp")

# Silent cleanup
cleanup_directory("/tmp/cache", verbose=false)

# Force cleanup of large directory
cleanup_directory("/tmp/large_data", force=true)
```
"""
function cleanup_directory(directory::AbstractString; verbose::Bool=true, force::Bool=false)
    # Check if directory exists
    if !isdir(directory)
        if verbose
            println("Directory does not exist, nothing to clean up: $(directory)")
        end
        return (existed=false, files_deleted=0, bytes_freed=0, human_readable_size="0 B")
    end
    
    # Calculate directory size and file count
    total_bytes = 0
    file_count = 0
    
    for (root, dirs, files) in walkdir(directory)
        for file in files
            file_path = joinpath(root, file)
            if isfile(file_path)  # Additional check to ensure it's a file
                try
                    total_bytes += filesize(file_path)
                    file_count += 1
                catch e
                    # Skip files that can't be read (permissions, etc.)
                    if verbose
                        @warn "Could not read file size for: $(file_path)"
                    end
                end
            end
        end
    end
    
    # Check if directory is empty
    if file_count == 0
        if verbose
            println("Directory is empty, removing: $(directory)")
        end
        try
            rm(directory, recursive=true)
        catch e
            if verbose
                @warn "Failed to remove empty directory: $(directory). Error: $(e)"
            end
            return (existed=true, files_deleted=0, bytes_freed=0, human_readable_size="0 B")
        end
        return (existed=true, files_deleted=0, bytes_freed=0, human_readable_size="0 B")
    end
    
    human_readable_size = Base.format_bytes(total_bytes)
    
    # Safety check for large directories
    if !force && (total_bytes > 1_000_000_000 || file_count > 10_000)  # 1GB or 10k files
        println("Warning: Large directory detected:")
        println("  Path: $(directory)")
        println("  Size: $(human_readable_size)")
        println("  Files: $(file_count)")
        print("Proceed with deletion? (y/N): ")
        response = readline()
        if lowercase(strip(response)) != "y"
            if verbose
                println("Cleanup cancelled by user")
            end
            return (existed=true, files_deleted=0, bytes_freed=0, human_readable_size="0 B")
        end
    end
    
    # Remove the directory
    try
        rm(directory, recursive=true)
        if verbose
            println("Cleanup completed:")
            println("  Directory: $(directory)")
            println("  Files deleted: $(file_count)")
            println("  Storage freed: $(human_readable_size)")
        end
        return (existed=true, files_deleted=file_count, bytes_freed=total_bytes, human_readable_size=human_readable_size)
    catch e
        if verbose
            @warn "Failed to remove directory: $(directory). Error: $(e)"
        end
        return (existed=true, files_deleted=0, bytes_freed=0, human_readable_size="0 B")
    end
end

# =============================================================================
# Sequence Type Utilities
# =============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determine the BioSequence type from an alphabet symbol.

Maps alphabet symbols to the corresponding BioSequences.jl type for 
type-safe sequence operations throughout the codebase.

# Arguments
- `alphabet::Symbol`: The alphabet symbol (`:DNA`, `:RNA`, or `:AA`)

# Returns
- `Type{<:BioSequences.BioSequence}`: The corresponding BioSequence type

# Examples
```julia
alphabet_to_biosequence_type(:DNA)  # Returns BioSequences.LongDNA{4}
alphabet_to_biosequence_type(:RNA)  # Returns BioSequences.LongRNA{4}
alphabet_to_biosequence_type(:AA)   # Returns BioSequences.LongAA
```
"""
function alphabet_to_biosequence_type(alphabet::Symbol)::Type
    if alphabet == :DNA
        return BioSequences.LongDNA{4}
    elseif alphabet == :RNA
        return BioSequences.LongRNA{4}
    elseif alphabet == :AA
        return BioSequences.LongAA
    else
        throw(ArgumentError("Unknown alphabet: $alphabet"))
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract sequence from FASTX record using dynamically determined type.

This function provides type-safe sequence extraction by using the appropriate
BioSequence type, avoiding hardcoded sequence types and string conversions.

# Arguments
- `record::Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}`: Input sequence record
- `sequence_type::Type{<:BioSequences.BioSequence}`: Target BioSequence type

# Returns
- `BioSequences.BioSequence`: Typed sequence from the record

# Examples
```julia
record = FASTX.FASTQ.Record("read1", "ATCG", "IIII")
seq_type = alphabet_to_biosequence_type(:DNA)
sequence = extract_typed_sequence(record, seq_type)
```
"""
function extract_typed_sequence(record::Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}, 
                               sequence_type::Type{<:BioSequences.BioSequence})
    return FASTX.sequence(sequence_type, record)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Detect alphabet and extract typed sequence from FASTX record in one step.

Convenience function that combines alphabet detection with type-safe sequence
extraction, ideal for workflows that need to determine sequence type once
at the beginning.

# Arguments
- `record::Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}`: Input sequence record

# Returns
- `Tuple{Symbol, BioSequences.BioSequence}`: (alphabet_symbol, typed_sequence)

# Examples
```julia
record = FASTX.FASTQ.Record("read1", "ATCG", "IIII")
alphabet, sequence = detect_and_extract_sequence(record)
# alphabet = :DNA, sequence = LongDNA{4} object
```
"""
function detect_and_extract_sequence(record::Union{FASTX.FASTA.Record, FASTX.FASTQ.Record})
    # First extract as string to detect alphabet
    sequence_string = FASTX.sequence(String, record)
    alphabet = detect_alphabet(sequence_string)
    
    # Then extract with proper type
    sequence_type = alphabet_to_biosequence_type(alphabet)
    typed_sequence = extract_typed_sequence(record, sequence_type)
    
    return (alphabet, typed_sequence)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate GC content (percentage of G and C bases) from a BioSequence.

This function calculates the percentage of guanine (G) and cytosine (C) bases
in a nucleotide sequence. Works with both DNA and RNA sequences.

# Arguments
- `sequence::BioSequences.LongSequence`: Input DNA or RNA sequence

# Returns
- `Float64`: GC content as a percentage (0.0-100.0)

# Examples
```julia
# Calculate GC content for DNA
dna_seq = BioSequences.LongDNA{4}("ATCGATCGATCG")
gc_percent = calculate_gc_content(dna_seq)

# Calculate GC content for RNA
rna_seq = BioSequences.LongRNA{4}("AUCGAUCGAUCG") 
gc_percent = calculate_gc_content(rna_seq)
```
"""
function calculate_gc_content(sequence::BioSequences.LongSequence)
    if length(sequence) == 0
        return 0.0
    end
    
    gc_count = 0
    total_bases = 0
    
    for base in sequence
        total_bases += 1
        if BioSequences.alphabet(sequence) isa BioSequences.DNAAlphabet
            if base == BioSequences.DNA_G || base == BioSequences.DNA_C
                gc_count += 1
            end
        elseif BioSequences.alphabet(sequence) isa BioSequences.RNAAlphabet
            if base == BioSequences.RNA_G || base == BioSequences.RNA_C
                gc_count += 1
            end
        else
            error("GC content calculation only supported for DNA and RNA sequences, got $(BioSequences.alphabet(sequence))")
        end
    end
    
    return (gc_count / total_bases) * 100.0
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate GC content from a string sequence.

Convenience function that accepts string input and converts to appropriate BioSequence.
Automatically detects DNA/RNA based on presence of T/U.

# Arguments
- `sequence::AbstractString`: Input DNA or RNA sequence as string

# Returns  
- `Float64`: GC content as a percentage (0.0-100.0)

# Examples
```julia
# Calculate GC content from string
gc_percent = calculate_gc_content("ATCGATCGATCG")
```
"""
function calculate_gc_content(sequence::AbstractString)
    if length(sequence) == 0
        return 0.0
    end
    
    # Detect if DNA (contains T) or RNA (contains U)
    sequence_upper = uppercase(sequence)
    if occursin('T', sequence_upper)
        bio_sequence = BioSequences.LongDNA{4}(sequence_upper)
    elseif occursin('U', sequence_upper)
        bio_sequence = BioSequences.LongRNA{4}(sequence_upper)  
    else
        # Default to DNA if ambiguous
        bio_sequence = BioSequences.LongDNA{4}(sequence_upper)
    end
    
    return calculate_gc_content(bio_sequence)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate GC content from FASTA/FASTQ records.

Processes multiple records and calculates overall GC content across all sequences.

# Arguments
- `records::AbstractVector{T}` where T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}`: Input records

# Returns
- `Float64`: Overall GC content as a percentage (0.0-100.0)

# Examples
```julia
# Calculate GC content from FASTA records
records = collect(FASTX.FASTA.Reader(open("sequences.fasta")))
gc_percent = calculate_gc_content(records)
```
"""
function calculate_gc_content(records::AbstractVector{T}) where {T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    if length(records) == 0
        return 0.0
    end
    
    total_gc = 0
    total_bases = 0
    
    for record in records
        sequence = FASTX.sequence(record)
        gc_count = 0
        
        for base in sequence
            total_bases += 1
            if BioSequences.alphabet(sequence) isa BioSequences.DNAAlphabet
                if base == BioSequences.DNA_G || base == BioSequences.DNA_C
                    gc_count += 1
                end
            elseif BioSequences.alphabet(sequence) isa BioSequences.RNAAlphabet
                if base == BioSequences.RNA_G || base == BioSequences.RNA_C
                    gc_count += 1
                end
            end
        end
        
        total_gc += gc_count
    end
    
    if total_bases == 0
        return 0.0
    end
    
    return (total_gc / total_bases) * 100.0
end

# =============================================================================
# Comprehensive Sequence Hashing Functions
# =============================================================================

"""
    _calculate_required_bytes(encoding::Symbol, encoded_length::Int) -> Int

Calculate minimum bytes needed to produce desired encoded length.
"""
function _calculate_required_bytes(encoding::Symbol, encoded_length::Int)::Int
    if encoding == :hex
        return ceil(Int, encoded_length / 2)
    elseif encoding == :base58
        # Base58 uses log2(58) ≈ 5.86 bits per character
        # Add safety margin to account for encoding variability
        bits_needed = encoded_length * log2(58)
        bytes_needed = ceil(Int, bits_needed / 8)
        # Add 1 extra byte as safety margin for Base58 encoding variability
        return bytes_needed + 1
    elseif encoding == :base64
        # Base64 uses 6 bits per character, but has padding
        return ceil(Int, encoded_length * 3 / 4)
    else
        error("Unsupported encoding: $encoding")
    end
end

"""
    _encode_hash_bytes(hash_bytes::Vector{UInt8}, encoding::Symbol, encoded_length::Union{Int,Missing}, allow_truncation::Bool) -> String

Internal function to encode hash bytes with specified encoding and length handling.
"""
function _encode_hash_bytes(hash_bytes::Vector{UInt8}, encoding::Symbol, encoded_length::Union{Int,Missing}, allow_truncation::Bool)::String
    if encoding == :hex
        encoded = bytes2hex(hash_bytes)
    elseif encoding == :base58
        encoded = String(Base58.base58encode(hash_bytes))
    elseif encoding == :base64
        encoded = Base64.base64encode(hash_bytes)
    else
        error("Unsupported encoding: $encoding. Supported: :hex, :base58, :base64")
    end
    
    # Handle length requirements
    if ismissing(encoded_length)
        return encoded
    elseif length(encoded) == encoded_length
        return encoded
    elseif length(encoded) > encoded_length
        # Auto-truncation for encodings with padding/safety margins
        # Our calculations in _calculate_required_bytes prevent "too short" errors
        # but may produce excess characters. Auto-truncate to requested length.
        if encoding == :base58 || encoding == :base64
            return encoded[1:encoded_length]
        elseif allow_truncation
            return encoded[1:encoded_length]
        else
            error("Encoded hash length ($(length(encoded))) exceeds requested length ($encoded_length). Set allow_truncation=true to truncate.")
        end
    else
        error("Cannot generate encoded hash of length $encoded_length from $(length(hash_bytes)) bytes with $encoding encoding (produces $(length(encoded)) characters)")
    end
end

"""
    crc32_checksum(data_to_hash::AbstractString; normalize_case::Bool=true) -> UInt32

Simple CRC32 checksum function returning raw checksum value.

# Arguments
- `data_to_hash::AbstractString`: Input data to checksum
- `normalize_case::Bool=true`: If true, converts input to uppercase before hashing

# Returns
- `UInt32`: Raw CRC32 checksum value
"""
function crc32_checksum(data_to_hash::AbstractString; normalize_case::Bool=true)::UInt32
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    return CRC32c.crc32c(Vector{UInt8}(normalized_data))
end

"""
    create_md5_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Generate MD5 hash with configurable encoding and length.
"""
function create_md5_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    hash_bytes = Vector{UInt8}(MD5.md5(Vector{UInt8}(normalized_data)))
    return _encode_hash_bytes(hash_bytes, encoding, encoded_length, allow_truncation)
end

"""
    create_sha1_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Generate SHA-1 hash with configurable encoding and length.
"""
function create_sha1_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    hash_bytes = SHA.sha1(Vector{UInt8}(normalized_data))
    return _encode_hash_bytes(hash_bytes, encoding, encoded_length, allow_truncation)
end

"""
    create_sha256_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Generate SHA-256 hash with configurable encoding and length.
"""
function create_sha256_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    hash_bytes = SHA.sha256(Vector{UInt8}(normalized_data))
    return _encode_hash_bytes(hash_bytes, encoding, encoded_length, allow_truncation)
end

"""
    create_sha512_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Generate SHA-512 hash with configurable encoding and length.
"""
function create_sha512_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    hash_bytes = SHA.sha512(Vector{UInt8}(normalized_data))
    return _encode_hash_bytes(hash_bytes, encoding, encoded_length, allow_truncation)
end

"""
    create_sha3_256_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Generate SHA-3 256-bit hash with configurable encoding and length.
"""
function create_sha3_256_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    hash_bytes = SHA.sha3_256(Vector{UInt8}(normalized_data))
    return _encode_hash_bytes(hash_bytes, encoding, encoded_length, allow_truncation)
end

"""
    create_sha3_512_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Generate SHA-3 512-bit hash with configurable encoding and length.
"""
function create_sha3_512_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    hash_bytes = SHA.sha3_512(Vector{UInt8}(normalized_data))
    return _encode_hash_bytes(hash_bytes, encoding, encoded_length, allow_truncation)
end

"""
    create_crc32_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false) -> String

Generate CRC32 checksum with configurable encoding and length.
"""
function create_crc32_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Union{Int,Missing}=missing, normalize_case::Bool=true, allow_truncation::Bool=false)::String
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    crc_val = CRC32c.crc32c(Vector{UInt8}(normalized_data))
    # Convert to 4 bytes (32 bits) in little-endian format
    hash_bytes = Vector{UInt8}(reinterpret(UInt8, [crc_val]))
    return _encode_hash_bytes(hash_bytes, encoding, encoded_length, allow_truncation)
end

"""
    create_blake3_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Int=64, normalize_case::Bool=true, allow_truncation::Bool=false, hash_bytes::Union{Int,Missing}=missing) -> String

Generate BLAKE3 hash with configurable encoding and length.

# Arguments
- `data_to_hash::AbstractString`: Input data to hash
- `encoding::Symbol=:hex`: Output encoding (:hex, :base58, :base64)
- `encoded_length::Int=64`: Desired output length (optimized for tree-of-life scale)
- `normalize_case::Bool=true`: If true, converts input to uppercase before hashing
- `allow_truncation::Bool=false`: Allow truncation if encoded_length < native length
- `hash_bytes::Union{Int,Missing}=missing`: Raw bytes to generate (auto-calculated if missing)
"""
function create_blake3_hash(data_to_hash::AbstractString; encoding::Symbol=:hex, encoded_length::Int=64, normalize_case::Bool=true, allow_truncation::Bool=false, hash_bytes::Union{Int,Missing}=missing)::String
    normalized_data = normalize_case ? uppercase(data_to_hash) : data_to_hash
    
    # Calculate required bytes if not specified
    bytes_needed = ismissing(hash_bytes) ? _calculate_required_bytes(encoding, encoded_length) : hash_bytes
    
    hasher = Blake3Hash.Blake3Ctx()
    Blake3Hash.update!(hasher, Vector{UInt8}(normalized_data))
    
    output_buffer = Vector{UInt8}(undef, bytes_needed)
    Blake3Hash.digest(hasher, output_buffer)
    
    return _encode_hash_bytes(output_buffer, encoding, encoded_length, allow_truncation)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Find the closest prime number to the given integer `n`.

Returns the nearest prime number to `n`. If two prime numbers are equally distant 
from `n`, returns the smaller one.

# Arguments
- `n::Int`: The input integer to find the nearest prime for

# Returns
- `Int`: The closest prime number to `n`
"""
function nearest_prime(n::Int)
    if n < 2
        return 2
    end
    next_p = Primes.nextprime(n)
    prev_p = Primes.prevprime(n)
    if n - prev_p <= next_p - n
        return prev_p
    else
        return next_p
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a sequence of Fibonacci numbers strictly less than the input value.

# Arguments
- `n::Int`: Upper bound (exclusive) for the Fibonacci sequence

# Returns
- `Vector{Int}`: Array containing Fibonacci numbers less than n
"""
function fibonacci_numbers_less_than(n::Int)
    if n <= 0
        return []
    elseif n == 1
        return [0]
    else
        fib = [0, 1]
        next_fib = fib[end] + fib[end-1]
        while next_fib < n
            push!(fib, next_fib)
            next_fib = fib[end] + fib[end-1]
        end
        return fib
    end
end
