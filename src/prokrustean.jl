"""
Prokrustean integration.
Reference: https://github.com/KoslickiLab/prokrustean
"""

const PROKRUSTEAN_REPO_URL = "https://github.com/KoslickiLab/prokrustean.git"
const PROKRUSTEAN_TOOLCHAIN_ENV = "prokrustean-build"
const PROKRUSTEAN_TOOLCHAIN_PACKAGES = [
    "cmake",
    "make",
    "git",
    "gcc_linux-64",
    "gxx_linux-64"
]

function _prokrustean_install_dir()
    return joinpath(dirname(dirname(pathof(Mycelia))), "deps", "prokrustean")
end

function _prokrustean_build_dir(install_dir::AbstractString)
    return joinpath(install_dir, "build")
end

function _prokrustean_executable(name::AbstractString; install_dir::AbstractString = _prokrustean_install_dir())
    return joinpath(_prokrustean_build_dir(install_dir), name)
end

function _prokrustean_missing_tools(tools::Tuple{Vararg{String}})
    missing = String[]
    for tool in tools
        if Sys.which(tool) === nothing
            push!(missing, tool)
        end
    end
    return missing
end

function _prokrustean_compiler()
    if haskey(ENV, "CXX") && !isempty(ENV["CXX"])
        return ENV["CXX"]
    end
    if Sys.which("c++") !== nothing
        return "c++"
    end
    if Sys.which("g++") !== nothing
        return "g++"
    end
    if Sys.which("clang++") !== nothing
        return "clang++"
    end
    return ""
end

function _prokrustean_compiler_supports_ranges(compiler::AbstractString)
    if isempty(compiler) || Sys.which(compiler) === nothing
        return false
    end

    workdir = mktempdir()
    try
        source_path = joinpath(workdir, "ranges_test.cpp")
        open(source_path, "w") do io
            write(io, "#include <ranges>\nint main(){return 0;}\n")
        end
        output_path = joinpath(workdir, "ranges_test")
        try
            run(`$compiler -std=c++20 $source_path -o $output_path`)
        catch
            return false
        end
        return true
    finally
        rm(workdir; recursive = true, force = true)
    end
end

function _prokrustean_module_avail(query::AbstractString)
    try
        return read(`bash -lc "module -t avail $(query) 2>&1"`, String)
    catch
        return ""
    end
end

function _prokrustean_module_version_key(modname::AbstractString)
    parts = split(modname, '/')
    if length(parts) < 2
        return (0, 0, 0)
    end
    version_str = parts[end]
    numbers = [tryparse(Int, part) for part in split(version_str, '.')]
    values = [isnothing(n) ? 0 : n for n in numbers]
    while length(values) < 3
        push!(values, 0)
    end
    return (values[1], values[2], values[3])
end

function _prokrustean_module_candidates()
    candidates = String[]
    for query in ("gcc", "llvm", "clang")
        output = _prokrustean_module_avail(query)
        if isempty(output)
            continue
        end
        for line in split(output, '\n')
            item = strip(line)
            if startswith(item, "$(query)/")
                push!(candidates, item)
            end
        end
    end

    gcc_candidates = filter(name -> startswith(name, "gcc/"), candidates)
    llvm_candidates = filter(
        name -> startswith(name, "llvm/") ||
                startswith(name, "clang/"), candidates)

    sort!(gcc_candidates; by = _prokrustean_module_version_key, rev = true)
    sort!(llvm_candidates; by = _prokrustean_module_version_key, rev = true)

    return vcat(gcc_candidates, llvm_candidates)
end

function _prokrustean_module_supports_ranges(module_name::AbstractString)
    workdir = mktempdir()
    try
        source_path = joinpath(workdir, "ranges_test.cpp")
        open(source_path, "w") do io
            write(io, "#include <ranges>\nint main(){return 0;}\n")
        end
        output_path = joinpath(workdir, "ranges_test")
        cmd = `bash -lc "module load $(module_name) && c++ -std=c++20 $(source_path) -o $(output_path)"`
        run(cmd)
        return true
    catch
        return false
    finally
        rm(workdir; recursive = true, force = true)
    end
end

function _prokrustean_build_with_module(
        module_name::AbstractString, source_dir::AbstractString, build_dir::AbstractString;
        cxx_standard::Int = 20)
    mkpath(build_dir)
    cmd = `bash -lc "module load $(module_name) && cmake -S $(source_dir) -B $(build_dir) -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=$(cxx_standard) -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_CXX_EXTENSIONS=OFF && cmake --build $(build_dir)"`
    run(cmd)
    return nothing
end

function _prokrustean_conda_available()
    return isfile(Mycelia.CONDA_RUNNER)
end

function _prokrustean_ensure_conda_toolchain_env(;
        env_name::AbstractString = PROKRUSTEAN_TOOLCHAIN_ENV,
        force::Bool = false)
    if !_prokrustean_conda_available()
        throw(ErrorException("Conda runner not available at $(Mycelia.CONDA_RUNNER)."))
    end

    if Mycelia.check_bioconda_env_is_installed(env_name)
        if force
            run(`$(Mycelia.CONDA_RUNNER) env remove -n $(env_name) -y`)
        else
            return env_name
        end
    end

    run(`$(Mycelia.CONDA_RUNNER) create -y -n $(env_name) -c conda-forge $(PROKRUSTEAN_TOOLCHAIN_PACKAGES)`)
    run(`$(Mycelia.CONDA_RUNNER) clean --all -y`)
    return env_name
end

function _prokrustean_build_with_conda(
        source_dir::AbstractString, build_dir::AbstractString;
        cxx_standard::Int = 20)
    env_name = _prokrustean_ensure_conda_toolchain_env()
    mkpath(build_dir)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) cmake -S $(source_dir) -B $(build_dir) -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=$(cxx_standard) -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_CXX_EXTENSIONS=OFF`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) cmake --build $(build_dir)`)
    return nothing
end

function _prokrustean_configure_and_build(
        source_dir::AbstractString, build_dir::AbstractString;
        cxx_standard::Int = 20)
    mkpath(build_dir)

    missing_build_tools = _prokrustean_missing_tools(("cmake", "make"))
    compiler = _prokrustean_compiler()
    if isempty(missing_build_tools) && !isempty(compiler) &&
       _prokrustean_compiler_supports_ranges(compiler)
        run(`cmake -S $source_dir -B $build_dir -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=$cxx_standard -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_CXX_EXTENSIONS=OFF`)
        run(`cmake --build $build_dir`)
        return nothing
    end

    if !isempty(missing_build_tools)
        @info "Missing system build tools: $(join(missing_build_tools, ", ")); attempting module or conda toolchain."
    end

    module_candidates = _prokrustean_module_candidates()
    for module_name in module_candidates
        if _prokrustean_module_supports_ranges(module_name)
            @info "Building Prokrustean with module $(module_name)"
            try
                _prokrustean_build_with_module(module_name, source_dir, build_dir; cxx_standard = cxx_standard)
                return nothing
            catch e
                @warn "Module build failed with $(module_name): $(e)"
            end
        end
    end

    if _prokrustean_conda_available()
        @info "Building Prokrustean with conda toolchain $(PROKRUSTEAN_TOOLCHAIN_ENV)"
        _prokrustean_build_with_conda(source_dir, build_dir; cxx_standard = cxx_standard)
        return nothing
    end

    message = "No suitable C++20 toolchain found. " *
              "Install a newer compiler (GCC >= 10 or Clang >= 13), " *
              "load a compiler module, or install the conda toolchain."
    throw(ErrorException(message))
end

function _prokrustean_ensure_executable(
        name::AbstractString; install_dir::AbstractString = _prokrustean_install_dir())
    executable = _prokrustean_executable(name; install_dir = install_dir)
    if !isfile(executable)
        install_prokrustean(dest_dir = install_dir)
    end
    if !isfile(executable)
        error("Prokrustean executable not found at $(executable). Please run install_prokrustean().")
    end
    return executable
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install the Prokrustean tool from source.
Clones the repository, builds with CMake/Make, and sets up executables.

# Arguments
- `dest_dir::String`: Directory to install Prokrustean (default: `deps/prokrustean`).
- `force::Bool=false`: Force re-installation if already present.
"""
function install_prokrustean(;
        dest_dir::String = _prokrustean_install_dir(),
        force::Bool = false
)
    main_executable = _prokrustean_executable("prokrustean"; install_dir = dest_dir)
    if isfile(main_executable) && !force
        @info "Prokrustean already installed at $(dest_dir). Use `force=true` to reinstall."
        return dest_dir
    end

    if isdir(dest_dir) && force
        rm(dest_dir; recursive = true, force = true)
    end

    if !isdir(dest_dir)
        mkpath(dirname(dest_dir))
        if isempty(_prokrustean_missing_tools(("git",)))
            run(`git clone --recursive $(PROKRUSTEAN_REPO_URL) $dest_dir`)
        elseif _prokrustean_conda_available()
            env_name = _prokrustean_ensure_conda_toolchain_env()
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) git clone --recursive $(PROKRUSTEAN_REPO_URL) $dest_dir`)
        else
            message = "Missing git and conda is not available to supply it. " *
                      "Install git or provide a conda runner."
            throw(ErrorException(message))
        end
    end

    build_dir = _prokrustean_build_dir(dest_dir)
    _prokrustean_configure_and_build(dest_dir, build_dir)

    if !isfile(main_executable)
        error("Prokrustean build completed but executable not found at $(main_executable).")
    end

    @info "Prokrustean installed successfully at $(dest_dir)"
    return dest_dir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Construct a Prokrustean graph from an eBWT file.

# Arguments
- `bwt_file::String`: Path to the input eBWT file.
- `output_file::String`: Path for the output binary graph file.

# Keywords
- `kmin::Int=20`: Minimum k-mer size to consider.
- `install_dir::String`: Path to Prokrustean installation (optional).
"""
function prokrustean_build_graph(
        bwt_file::String,
        output_file::String;
        kmin::Int = 20,
        install_dir::String = _prokrustean_install_dir()
)
    if !isfile(bwt_file)
        throw(ArgumentError("BWT file not found: $(bwt_file)"))
    end
    if kmin < 1
        throw(ArgumentError("kmin must be a positive integer, got $(kmin)"))
    end

    executable = _prokrustean_ensure_executable("prokrustean"; install_dir = install_dir)
    mkpath(dirname(output_file))

    cmd = `$executable -i $bwt_file -l $kmin`
    run(cmd)

    default_output = "$(bwt_file).prokrustean"
    if isfile(default_output)
        if default_output != output_file
            mv(default_output, output_file; force = true)
        end
        return output_file
    end

    if isfile(output_file)
        return output_file
    end

    error("Prokrustean did not produce expected output at $(output_file) or $(default_output).")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the number of distinct k-mers for k=kmin...L using a Prokrustean graph.

# Arguments
- `graph_file::String`: Path to the Prokrustean graph file.
- `output_file::String`: Path to save the counts (CSV/TSV format).

# Keywords
- `install_dir::String`: Path to Prokrustean installation.
"""
function prokrustean_kmer_count(
        graph_file::String,
        output_file::String;
        install_dir::String = _prokrustean_install_dir()
)
    if !isfile(graph_file)
        throw(ArgumentError("Graph file not found: $(graph_file)"))
    end
    executable = _prokrustean_ensure_executable("prokrustean_kmer_count"; install_dir = install_dir)
    mkpath(dirname(output_file))

    cmd = pipeline(`$executable -p $graph_file`, stdout = output_file)
    run(cmd)
    return output_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the number of maximal unitigs of de Bruijn graphs for k=kmin...L.

# Arguments
- `graph_file::String`: Path to the Prokrustean graph file.
- `output_file::String`: Path to save the counts.

# Keywords
- `install_dir::String`: Path to Prokrustean installation.
"""
function prokrustean_unitig_count(
        graph_file::String,
        output_file::String;
        install_dir::String = _prokrustean_install_dir()
)
    if !isfile(graph_file)
        throw(ArgumentError("Graph file not found: $(graph_file)"))
    end
    executable = _prokrustean_ensure_executable("prokrustean_unitig_count"; install_dir = install_dir)
    mkpath(dirname(output_file))

    cmd = pipeline(`$executable -p $graph_file`, stdout = output_file)
    run(cmd)
    return output_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute Bray-Curtis dissimilarities between samples for k=kmin...L.

# Arguments
- `graph_file::String`: Path to the Prokrustean graph file (merged).
- `sample_ids_file::String`: Path to the sample IDs file (rows of 0s and 1s).
- `output_file::String`: Path to save the dissimilarity matrix/values.

# Keywords
- `install_dir::String`: Path to Prokrustean installation.
"""
function prokrustean_braycurtis(
        graph_file::String,
        sample_ids_file::String,
        output_file::String;
        install_dir::String = _prokrustean_install_dir()
)
    if !isfile(graph_file)
        throw(ArgumentError("Graph file not found: $(graph_file)"))
    end
    if !isfile(sample_ids_file)
        throw(ArgumentError("Sample IDs file not found: $(sample_ids_file)"))
    end
    executable = _prokrustean_ensure_executable("prokrustean_braycurtis"; install_dir = install_dir)
    mkpath(dirname(output_file))

    cmd = pipeline(`$executable -p $graph_file -s $sample_ids_file`, stdout = output_file)
    run(cmd)
    return output_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count vertex degrees of the overlap graph for k >= kmin.

# Arguments
- `graph_file::String`: Path to the Prokrustean graph file.
- `output_file::String`: Path to save the degree counts.

# Keywords
- `install_dir::String`: Path to Prokrustean installation.
"""
function prokrustean_overlap(
        graph_file::String,
        output_file::String;
        install_dir::String = _prokrustean_install_dir()
)
    if !isfile(graph_file)
        throw(ArgumentError("Graph file not found: $(graph_file)"))
    end
    executable = _prokrustean_ensure_executable("prokrustean_overlap"; install_dir = install_dir)
    mkpath(dirname(output_file))

    cmd = pipeline(`$executable -p $graph_file`, stdout = output_file)
    run(cmd)
    return output_file
end
