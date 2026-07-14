"""
Wrapper for Autocycler and a conservative paired-short polishing pipeline.

References:
- Autocycler: https://github.com/rrwick/Autocycler
- Polypolish: https://github.com/rrwick/Polypolish
- Pypolca: https://github.com/gbouras13/pypolca
"""

const AUTOCYCLER_ENV_NAME = "autocycler"
const AUTOCYCLER_SCRIPT_REVISION =
    "c98b126eb45727584623041db1bfdbdaf7aa0923"
const AUTOCYCLER_SCRIPT_SHA256 =
    "42d9b41385c2095ba05a910511f490cbc97e38b7965e0f8f0978b7a1e1477eaa"
const AUTOCYCLER_SCRIPT_URL =
    "https://raw.githubusercontent.com/rrwick/Autocycler/" *
    AUTOCYCLER_SCRIPT_REVISION *
    "/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/" *
    "autocycler_full.sh"
const AUTOCYCLER_REQUIRED_PACKAGE_SPECS = (
    (name = "autocycler", constraint = :exact, version = "0.5.2"),
    (name = "bwa", constraint = :minimum, version = "0.7.17"),
    (name = "canu", constraint = :minimum, version = "2.3"),
    (name = "flye", constraint = :minimum, version = "2.9.6"),
    (name = "metamdbg", constraint = :minimum, version = "1.0"),
    (name = "miniasm", constraint = :minimum, version = "0.3"),
    (name = "minimap2", constraint = :minimum, version = "2.28"),
    (name = "minipolish", constraint = :minimum, version = "0.2.0"),
    (name = "myloasm", constraint = :minimum, version = "0.1.0"),
    (
        name = "necat",
        constraint = :minimum,
        version = "0.0.1_update20200803",
    ),
    (name = "nextdenovo", constraint = :minimum, version = "2.5.2"),
    (name = "nextpolish", constraint = :minimum, version = "1.4.1"),
    (name = "parallel", constraint = :present, version = ""),
    (name = "plassembler", constraint = :minimum, version = "1.8.0"),
    (name = "polypolish", constraint = :minimum, version = "0.6.0"),
    (name = "pypolca", constraint = :minimum, version = "0.3.1"),
    (name = "racon", constraint = :minimum, version = "1.5.0"),
    (name = "raven-assembler", constraint = :minimum, version = "1.8.3"),
    (name = "sed", constraint = :present, version = ""),
    (name = "wtdbg", constraint = :minimum, version = "2.5"),
)
const AUTOCYCLER_REQUIRED_PACKAGES =
    map(specification -> specification.name, AUTOCYCLER_REQUIRED_PACKAGE_SPECS)
const AUTOCYCLER_READ_TYPES = (
    "ont_r9",
    "ont_r10",
    "pacbio_clr",
    "pacbio_hifi",
)

function _autocycler_paths()::Tuple{String, String, String}
    install_dir = joinpath(dirname(dirname(pathof(Mycelia))), "deps", "autocycler")
    script_path = joinpath(install_dir, "autocycler_full.sh")
    env_file_path = joinpath(install_dir, "environment.yml")
    return install_dir, script_path, env_file_path
end

function _autocycler_sha256(path::AbstractString)::String
    return SHA.bytes2hex(SHA.sha256(read(path)))
end

function _autocycler_script_is_verified(path::AbstractString)::Bool
    return isfile(path) && filesize(path) > 0 &&
           _autocycler_sha256(path) == AUTOCYCLER_SCRIPT_SHA256
end

function _install_verified_autocycler_script!(
        script_path::AbstractString;
        downloader::Function = Downloads.download,
)::String
    normalized_script_path = abspath(script_path)
    mkpath(dirname(normalized_script_path))
    temporary_path, temporary_io = mktemp(dirname(normalized_script_path))
    close(temporary_io)
    try
        downloader(AUTOCYCLER_SCRIPT_URL, temporary_path)
        if !isfile(temporary_path) || filesize(temporary_path) == 0
            throw(ErrorException("Downloaded Autocycler script is empty."))
        end
        actual_sha256 = _autocycler_sha256(temporary_path)
        if actual_sha256 != AUTOCYCLER_SCRIPT_SHA256
            throw(
                ErrorException(
                    "Autocycler script checksum mismatch for revision " *
                    "$(AUTOCYCLER_SCRIPT_REVISION): expected " *
                    "$(AUTOCYCLER_SCRIPT_SHA256), got $(actual_sha256).",
                ),
            )
        end
        mv(temporary_path, normalized_script_path; force = true)
        Base.chmod(normalized_script_path, 0o755)
    finally
        rm(temporary_path; force = true)
    end
    return normalized_script_path
end

function _autocycler_environment_packages(;
        conda_runner::AbstractString = CONDA_RUNNER,
        command_reader::Function = command -> read(command, String),
)::Dict{String, String}
    command = Cmd(
        String[
            String(conda_runner),
            "list",
            "-n",
            AUTOCYCLER_ENV_NAME,
            "--json",
        ],
    )
    package_records = JSON.parse(command_reader(command))
    package_records isa AbstractVector || throw(
        ErrorException("Conda package inventory was not a JSON array."),
    )
    versions = Dict{String, String}()
    for package_record in package_records
        package_record isa AbstractDict || continue
        name = get(package_record, "name", nothing)
        version = get(package_record, "version", nothing)
        if name isa AbstractString && version isa AbstractString
            versions[String(name)] = String(version)
        end
    end
    return versions
end

function _missing_autocycler_packages(
        versions::AbstractDict{<:AbstractString, <:AbstractString},
)::Vector{String}
    return String[
        package for package in AUTOCYCLER_REQUIRED_PACKAGES if
        !haskey(versions, package)
    ]
end

function _autocycler_numeric_version(
        version::AbstractString,
)::Union{Nothing, VersionNumber}
    version_match = match(r"^\d+(?:\.\d+){0,2}", String(version))
    version_match === nothing && return nothing
    return try
        VersionNumber(version_match.match)
    catch
        nothing
    end
end

function _autocycler_version_at_least(
        actual::AbstractString,
        minimum::AbstractString,
)::Bool
    actual_version = _autocycler_numeric_version(actual)
    minimum_version = _autocycler_numeric_version(minimum)
    if actual_version === nothing || minimum_version === nothing
        return false
    elseif actual_version != minimum_version
        return actual_version > minimum_version
    end

    minimum_match = match(r"^\d+(?:\.\d+){0,2}(.*)$", String(minimum))
    actual_match = match(r"^\d+(?:\.\d+){0,2}(.*)$", String(actual))
    minimum_suffix = something(only(minimum_match.captures))
    isempty(minimum_suffix) && return true
    actual_suffix = something(only(actual_match.captures))
    return !isempty(actual_suffix) && !isless(actual_suffix, minimum_suffix)
end

function _autocycler_package_issues(
        versions::AbstractDict{<:AbstractString, <:AbstractString},
)::Vector{String}
    issues = String[]
    for specification in AUTOCYCLER_REQUIRED_PACKAGE_SPECS
        if !haskey(versions, specification.name)
            push!(issues, "$(specification.name) is missing")
            continue
        end
        specification.constraint == :present && continue

        actual = String(versions[specification.name])
        if specification.constraint == :exact
            actual == specification.version || push!(
                issues,
                "$(specification.name) must equal $(specification.version), " *
                "got $(actual)",
            )
            continue
        end

        _autocycler_version_at_least(actual, specification.version) || push!(
            issues,
            "$(specification.name) must be at least $(specification.version), " *
            "got $(actual)",
        )
    end
    return issues
end

function _ensure_autocycler_packages!(
        package_inspector::Function,
        installer::Function,
)::Dict{String, String}
    versions = package_inspector()
    package_issues = _autocycler_package_issues(versions)
    if !isempty(package_issues)
        @warn "Autocycler environment is stale; recreating it before assembly" package_issues
        installer(; force = true)
        versions = package_inspector()
        package_issues = _autocycler_package_issues(versions)
    end
    isempty(package_issues) || throw(
        ErrorException(
            "Autocycler environment has missing or incompatible required " *
            "packages after recreation: $(join(package_issues, "; ")).",
        ),
    )
    return Dict{String, String}(versions)
end

function _autocycler_toolchain_metadata(
        versions::AbstractDict{<:AbstractString, <:AbstractString},
)::Dict{String, Any}
    _, script_path, env_file_path = _autocycler_paths()
    return Dict{String, Any}(
        "autocycler_script_revision" => AUTOCYCLER_SCRIPT_REVISION,
        "autocycler_script_sha256" => _autocycler_sha256(script_path),
        "environment_spec_sha256" => _autocycler_sha256(env_file_path),
        "packages" => Dict{String, String}(versions),
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install Autocycler and the short-read polishing tools in a shared conda
environment.

The package-bundled `environment.yml` is authoritative because it extends the
upstream Autocycler environment with BWA, Polypolish, and Pypolca. The upstream
automation script is pinned to `AUTOCYCLER_SCRIPT_REVISION` and verified against
`AUTOCYCLER_SCRIPT_SHA256` before it can be executed. `force=true` recreates the
environment and refreshes that verified script without overwriting the bundled
environment specification.

# Keywords
- `force::Bool=false`: Recreate the environment and re-download the pinned script.
"""
function install_autocycler(;
        force::Bool = false,
        downloader::Function = Downloads.download,
        paths::Tuple{String, String, String} = _autocycler_paths(),
)::String
    install_dir, script_path, env_file_path = paths
    mkpath(install_dir)

    if !isfile(env_file_path) || filesize(env_file_path) == 0
        throw(
            ErrorException(
                "Bundled Autocycler environment file is missing or empty: " *
                "$(env_file_path). Reinstall Mycelia before installing Autocycler.",
            ),
        )
    end

    create_conda_env_from_yaml(
        env_file_path,
        AUTOCYCLER_ENV_NAME;
        force = force,
    )

    if force || !_autocycler_script_is_verified(script_path)
        @info "Installing pinned, checksum-verified autocycler_full.sh script..."
        _install_verified_autocycler_script!(
            script_path;
            downloader = downloader,
        )
    end

    if !_autocycler_script_is_verified(script_path)
        throw(
            ErrorException(
                "Autocycler script installation failed verification: " *
                "$(script_path)",
            ),
        )
    end

    @info "Autocycler installed successfully."
    return script_path
end

function _require_nonempty_autocycler_file(
        path::AbstractString,
        label::AbstractString,
)::String
    normalized_path = abspath(path)
    if !isfile(normalized_path)
        throw(ArgumentError("$(label) not found: $(normalized_path)"))
    end
    if filesize(normalized_path) == 0
        throw(ArgumentError("$(label) is empty: $(normalized_path)"))
    end
    return normalized_path
end

function _autocycler_pair_identifier(identifier::AbstractString)::String
    first_token = first(split(String(identifier)))
    return replace(first_token, r"/[12]$" => "")
end

function _autocycler_pair_role(identifier::AbstractString)::Union{Nothing, Int}
    first_token = first(split(String(identifier)))
    role_match = match(r"/([12])$", first_token)
    return role_match === nothing ? nothing :
           parse(Int, something(only(role_match.captures)))
end

function _validate_autocycler_fastq(
        path::AbstractString,
        label::AbstractString,
)::Int
    reader = Mycelia.open_fastx(path)
    record_count = 0
    try
        for record in reader
            record isa FASTX.FASTQ.Record || throw(ArgumentError(
                "$(label) must be a FASTQ file: $(abspath(path))",
            ))
            record_count += 1
        end
    finally
        close(reader)
    end
    record_count > 0 || throw(ArgumentError("$(label) must be non-empty."))
    return record_count
end

function _validate_autocycler_paired_fastqs(
        short_reads_1::AbstractString,
        short_reads_2::AbstractString,
)::Int
    !Base.Filesystem.samefile(short_reads_1, short_reads_2) || throw(ArgumentError(
        "Autocycler paired short-read R1 and R2 must be distinct files.",
    ))
    reader_1 = Mycelia.open_fastx(short_reads_1)
    reader_2 = Mycelia.open_fastx(short_reads_2)
    pair_count = 0
    try
        next_1 = iterate(reader_1)
        next_2 = iterate(reader_2)
        while next_1 !== nothing || next_2 !== nothing
            if next_1 === nothing || next_2 === nothing
                throw(
                    ArgumentError(
                        "Autocycler paired short reads have different counts " *
                        "after $(pair_count) complete pairs.",
                    ),
                )
            end
            record_1, state_1 = next_1
            record_2, state_2 = next_2
            pair_count += 1
            if !(record_1 isa FASTX.FASTQ.Record) ||
               !(record_2 isa FASTX.FASTQ.Record)
                throw(
                    ArgumentError(
                        "Autocycler paired short-read inputs must be FASTQ files.",
                    ),
                )
            end
            identifier_1 = String(FASTX.identifier(record_1))
            identifier_2 = String(FASTX.identifier(record_2))
            role_1 = _autocycler_pair_role(identifier_1)
            role_2 = _autocycler_pair_role(identifier_2)
            roles_valid = (role_1 === nothing && role_2 === nothing) ||
                          (role_1 == 1 && role_2 == 2)
            roles_valid || throw(ArgumentError(
                "Autocycler paired short reads have invalid explicit mate " *
                "roles at record $(pair_count): " *
                "R1=$(repr(identifier_1)), R2=$(repr(identifier_2)); " *
                "expected /1 then /2.",
            ))
            if _autocycler_pair_identifier(identifier_1) !=
               _autocycler_pair_identifier(identifier_2)
                throw(
                    ArgumentError(
                        "Autocycler paired short reads are out of sync at " *
                        "record $(pair_count): R1=$(repr(identifier_1)), " *
                        "R2=$(repr(identifier_2)).",
                    ),
                )
            end
            next_1 = iterate(reader_1, state_1)
            next_2 = iterate(reader_2, state_2)
        end
    finally
        close(reader_1)
        close(reader_2)
    end
    pair_count > 0 || throw(
        ArgumentError("Autocycler paired short reads must be non-empty."),
    )
    return pair_count
end

function _validate_autocycler_parameters(
        threads::Integer,
        jobs::Integer,
        read_type::AbstractString,
)::String
    if threads < 1
        throw(ArgumentError("threads must be positive, got $(threads)"))
    end
    if jobs < 1
        throw(ArgumentError("jobs must be positive, got $(jobs)"))
    end

    normalized_read_type = String(read_type)
    if !(normalized_read_type in AUTOCYCLER_READ_TYPES)
        allowed = join(AUTOCYCLER_READ_TYPES, ", ")
        throw(
            ArgumentError(
                "read_type must be one of $(allowed); got $(normalized_read_type)",
            ),
        )
    end
    return normalized_read_type
end

function _validate_autocycler_output_dir(out_dir::AbstractString)::String
    normalized_out_dir = abspath(out_dir)
    if isfile(normalized_out_dir)
        throw(ArgumentError("Autocycler output path is a file: $(normalized_out_dir)"))
    end
    if isdir(normalized_out_dir) && !isempty(readdir(normalized_out_dir))
        throw(
            ArgumentError(
                "Autocycler output directory must be empty: $(normalized_out_dir)",
            ),
        )
    end
    return normalized_out_dir
end

function _prepare_autocycler_output_dir(out_dir::AbstractString)::String
    normalized_out_dir = _validate_autocycler_output_dir(out_dir)
    mkpath(normalized_out_dir)
    return normalized_out_dir
end

function _autocycler_conda_command(
        arguments::Vector{String},
        work_dir::AbstractString;
        conda_runner::AbstractString = CONDA_RUNNER,
)::Cmd
    command_arguments = String[
        String(conda_runner),
        "run",
        "--live-stream",
        "-n",
        AUTOCYCLER_ENV_NAME,
    ]
    append!(command_arguments, arguments)
    command = Cmd(command_arguments)
    return Cmd(command; dir = String(work_dir))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build the exact upstream Autocycler long-read command plan without executing it.
The upstream script accepts only long reads, threads per assembly job, concurrent
jobs, and a long-read type. It writes both
`autocycler_out/consensus_assembly.fasta` and
`autocycler_out/consensus_assembly.gfa` relative to its working directory.
"""
function _autocycler_command_plan(
        long_reads::AbstractString,
        out_dir::AbstractString;
        threads::Integer = max(Sys.CPU_THREADS, 1),
        jobs::Integer = 1,
        read_type::AbstractString = "ont_r10",
        script_path::AbstractString = _autocycler_paths()[2],
        conda_runner::AbstractString = CONDA_RUNNER,
)::NamedTuple
    normalized_read_type = _validate_autocycler_parameters(
        threads,
        jobs,
        read_type,
    )
    normalized_long_reads = abspath(long_reads)
    normalized_out_dir = abspath(out_dir)
    normalized_script_path = abspath(script_path)
    autocycler_out_dir = joinpath(normalized_out_dir, "autocycler_out")
    graph = joinpath(autocycler_out_dir, "consensus_assembly.gfa")
    assembly = joinpath(autocycler_out_dir, "consensus_assembly.fasta")

    command = _autocycler_conda_command(
        String[
            "bash",
            normalized_script_path,
            normalized_long_reads,
            string(threads),
            string(jobs),
            normalized_read_type,
        ],
        normalized_out_dir;
        conda_runner = conda_runner,
    )
    step = (
        name = :autocycler,
        command = command,
        stdout = nothing,
        expected_outputs = String[graph, assembly],
    )
    return (
        steps = (step,),
        outdir = normalized_out_dir,
        assembly = assembly,
        graph = graph,
    )
end

function _autocycler_polishing_command_plan(
        assembly::AbstractString,
        short_reads_1::AbstractString,
        short_reads_2::AbstractString,
        out_dir::AbstractString;
        threads::Integer = max(Sys.CPU_THREADS, 1),
        polypolish_careful::Bool = true,
        conda_runner::AbstractString = CONDA_RUNNER,
)::NamedTuple
    if threads < 1
        throw(ArgumentError("threads must be positive, got $(threads)"))
    end

    normalized_assembly = abspath(assembly)
    normalized_short_reads_1 = abspath(short_reads_1)
    normalized_short_reads_2 = abspath(short_reads_2)
    polishing_dir = joinpath(abspath(out_dir), "short_read_polishing")
    alignments_1 = joinpath(polishing_dir, "alignments_1.sam")
    alignments_2 = joinpath(polishing_dir, "alignments_2.sam")
    filtered_1 = joinpath(polishing_dir, "filtered_1.sam")
    filtered_2 = joinpath(polishing_dir, "filtered_2.sam")
    polypolish_assembly = joinpath(polishing_dir, "polypolish.fasta")
    pypolca_dir = joinpath(polishing_dir, "pypolca")
    pypolca_prefix = "autocycler_polished"
    final_assembly = joinpath(
        pypolca_dir,
        "$(pypolca_prefix)_corrected.fasta",
    )
    pypolca_report = joinpath(pypolca_dir, "$(pypolca_prefix).report")
    bwa_index_files = String[
        "$(normalized_assembly).$(extension)" for
        extension in ("amb", "ann", "bwt", "pac", "sa")
    ]

    bwa_index = (
        name = :bwa_index,
        command = _autocycler_conda_command(
            String["bwa", "index", normalized_assembly],
            polishing_dir;
            conda_runner = conda_runner,
        ),
        stdout = nothing,
        expected_outputs = bwa_index_files,
    )
    bwa_mem_1 = (
        name = :bwa_mem_1,
        command = _autocycler_conda_command(
            String[
                "bwa",
                "mem",
                "-t",
                string(threads),
                "-a",
                normalized_assembly,
                normalized_short_reads_1,
            ],
            polishing_dir;
            conda_runner = conda_runner,
        ),
        stdout = alignments_1,
        expected_outputs = String[alignments_1],
    )
    bwa_mem_2 = (
        name = :bwa_mem_2,
        command = _autocycler_conda_command(
            String[
                "bwa",
                "mem",
                "-t",
                string(threads),
                "-a",
                normalized_assembly,
                normalized_short_reads_2,
            ],
            polishing_dir;
            conda_runner = conda_runner,
        ),
        stdout = alignments_2,
        expected_outputs = String[alignments_2],
    )
    polypolish_filter = (
        name = :polypolish_filter,
        command = _autocycler_conda_command(
            String[
                "polypolish",
                "filter",
                "--in1",
                alignments_1,
                "--in2",
                alignments_2,
                "--out1",
                filtered_1,
                "--out2",
                filtered_2,
            ],
            polishing_dir;
            conda_runner = conda_runner,
        ),
        stdout = nothing,
        expected_outputs = String[filtered_1, filtered_2],
    )
    polypolish_arguments = String["polypolish", "polish"]
    if polypolish_careful
        push!(polypolish_arguments, "--careful")
    end
    append!(
        polypolish_arguments,
        String[normalized_assembly, filtered_1, filtered_2],
    )
    polypolish = (
        name = :polypolish,
        command = _autocycler_conda_command(
            polypolish_arguments,
            polishing_dir;
            conda_runner = conda_runner,
        ),
        stdout = polypolish_assembly,
        expected_outputs = String[polypolish_assembly],
    )
    pypolca = (
        name = :pypolca,
        command = _autocycler_conda_command(
            String[
                "pypolca",
                "run",
                "-a",
                polypolish_assembly,
                "-1",
                normalized_short_reads_1,
                "-2",
                normalized_short_reads_2,
                "-t",
                string(threads),
                "-o",
                pypolca_dir,
                "--careful",
                "-p",
                pypolca_prefix,
            ],
            polishing_dir;
            conda_runner = conda_runner,
        ),
        stdout = nothing,
        expected_outputs = String[final_assembly, pypolca_report],
    )

    return (
        steps = (
            bwa_index,
            bwa_mem_1,
            bwa_mem_2,
            polypolish_filter,
            polypolish,
            pypolca,
        ),
        polishing_dir = polishing_dir,
        polypolish_assembly = polypolish_assembly,
        pypolca_dir = pypolca_dir,
        pypolca_report = pypolca_report,
        bwa_index_files = bwa_index_files,
        intermediate_files = String[
            bwa_index_files...,
            alignments_1,
            alignments_2,
            filtered_1,
            filtered_2,
        ],
        assembly = final_assembly,
    )
end

function _cleanup_autocycler_polishing_intermediates!(
        paths::AbstractVector{<:AbstractString},
)::Nothing
    for path in paths
        rm(path; force = true)
    end
    return nothing
end

function _default_autocycler_step_runner(step::NamedTuple)::Nothing
    if isnothing(step.stdout)
        Base.run(step.command)
    else
        mkpath(dirname(step.stdout))
        Base.open(step.stdout, "w") do io
            Base.run(Base.pipeline(step.command; stdout = io))
        end
    end
    return nothing
end

function _execute_autocycler_steps(
        steps::Tuple;
        runner::Function = _default_autocycler_step_runner,
)::Nothing
    for step in steps
        try
            runner(step)
        catch error
            if error isa InterruptException
                rethrow()
            elseif error isa Base.ProcessFailedException
                message = sprint(showerror, error)
                throw(
                    ErrorException(
                        "Autocycler workflow command $(step.name) failed. " *
                        "Recreate the tool environment with " *
                        "install_autocycler(force=true) if a command is missing. " *
                        "Cause: $(message)",
                    ),
                )
            end
            rethrow()
        end

        for output_path in step.expected_outputs
            if !isfile(output_path) || filesize(output_path) == 0
                throw(
                    ErrorException(
                        "Autocycler workflow step $(step.name) did not create " *
                        "a nonempty artifact: $(output_path)",
                    ),
                )
            end
        end
    end
    return nothing
end

function _ensure_autocycler_installed(;
        package_inspector::Function = _autocycler_environment_packages,
        installer::Function = install_autocycler,
)::Dict{String, Any}
    _, script_path, _ = _autocycler_paths()
    environment_installed = check_bioconda_env_is_installed(
        AUTOCYCLER_ENV_NAME,
    )
    if !environment_installed || !_autocycler_script_is_verified(script_path)
        @warn "Autocycler is not fully installed. Installing now..."
        installer()
    end
    if !_autocycler_script_is_verified(script_path)
        throw(
            ErrorException(
                "Autocycler script is unavailable or failed checksum " *
                "verification: $(script_path)",
            ),
        )
    end
    versions = _ensure_autocycler_packages!(package_inspector, installer)
    return _autocycler_toolchain_metadata(versions)
end

function _run_autocycler(
        long_reads::AbstractString,
        out_dir::AbstractString;
        threads::Integer = max(Sys.CPU_THREADS, 1),
        jobs::Integer = 1,
        read_type::AbstractString = "ont_r10",
        dependency_checker::Function = _ensure_autocycler_installed,
        runner::Function = _default_autocycler_step_runner,
        validate_long_reads::Bool = true,
)::NamedTuple
    normalized_long_reads = _require_nonempty_autocycler_file(
        long_reads,
        "Long-read FASTQ",
    )
    validate_long_reads &&
        _validate_autocycler_fastq(normalized_long_reads, "Long-read input")
    normalized_read_type = _validate_autocycler_parameters(
        threads,
        jobs,
        read_type,
    )
    normalized_out_dir = _validate_autocycler_output_dir(out_dir)
    toolchain = dependency_checker()
    normalized_out_dir = _prepare_autocycler_output_dir(normalized_out_dir)
    _, script_path, _ = _autocycler_paths()
    plan = _autocycler_command_plan(
        normalized_long_reads,
        normalized_out_dir;
        threads = threads,
        jobs = jobs,
        read_type = normalized_read_type,
        script_path = script_path,
    )

    @info "Starting long-read-only Autocycler pipeline..."
    _execute_autocycler_steps(plan.steps; runner = runner)
    assembly = _require_nonempty_autocycler_file(
        plan.assembly,
        "Autocycler consensus FASTA",
    )
    @info "Autocycler pipeline complete" out_dir = normalized_out_dir

    return (
        outdir = normalized_out_dir,
        assembly = assembly,
        graph = plan.graph,
        toolchain = toolchain,
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the upstream long-read-only Autocycler pipeline.

Autocycler consumes one long-read FASTQ and writes its fixed `autocycler_out`
directory relative to the isolated `out_dir` working directory. Its authoritative
consensus FASTA and companion GFA are both validated and returned without
rewriting either artifact. The output directory must be absent or empty.

# Keywords
- `long_reads::AbstractString`: Nonempty long-read FASTQ.
- `out_dir::AbstractString`: Empty isolated working/output directory.
- `threads::Integer`: Threads per input assembler job.
- `jobs::Integer`: Number of simultaneous assembler jobs.
- `read_type::AbstractString`: One of `ont_r9`, `ont_r10`, `pacbio_clr`, or
  `pacbio_hifi`.

# Returns
A named tuple with `outdir`, `assembly`, `graph`, and exact `toolchain`
provenance.
"""
function run_autocycler(;
        long_reads::AbstractString,
        out_dir::AbstractString,
        threads::Integer = max(Sys.CPU_THREADS, 1),
        jobs::Integer = 1,
        read_type::AbstractString = "ont_r10",
)::NamedTuple
    return _run_autocycler(
        long_reads,
        out_dir;
        threads = threads,
        jobs = jobs,
        read_type = read_type,
    )
end

function _run_autocycler_polished(
        long_reads::AbstractString,
        short_reads_1::AbstractString,
        short_reads_2::AbstractString,
        out_dir::AbstractString;
        threads::Integer = max(Sys.CPU_THREADS, 1),
        jobs::Integer = 1,
        read_type::AbstractString = "ont_r10",
        polypolish_careful::Bool = true,
        keep_intermediates::Bool = false,
        dependency_checker::Function = _ensure_autocycler_installed,
        runner::Function = _default_autocycler_step_runner,
)::NamedTuple
    normalized_long_reads = _require_nonempty_autocycler_file(
        long_reads,
        "Long-read FASTQ",
    )
    normalized_short_reads_1 = _require_nonempty_autocycler_file(
        short_reads_1,
        "Paired short-read R1 FASTQ",
    )
    normalized_short_reads_2 = _require_nonempty_autocycler_file(
        short_reads_2,
        "Paired short-read R2 FASTQ",
    )
    _validate_autocycler_fastq(normalized_long_reads, "Long-read input")
    _validate_autocycler_paired_fastqs(
        normalized_short_reads_1,
        normalized_short_reads_2,
    )
    normalized_read_type = _validate_autocycler_parameters(
        threads,
        jobs,
        read_type,
    )
    normalized_out_dir = _validate_autocycler_output_dir(out_dir)
    toolchain = dependency_checker()

    autocycler_result = _run_autocycler(
        normalized_long_reads,
        normalized_out_dir;
        threads = threads,
        jobs = jobs,
        read_type = normalized_read_type,
        dependency_checker = () -> toolchain,
        runner = runner,
        validate_long_reads = false,
    )
    normalized_out_dir = autocycler_result.outdir

    polishing_plan = _autocycler_polishing_command_plan(
        autocycler_result.assembly,
        normalized_short_reads_1,
        normalized_short_reads_2,
        normalized_out_dir;
        threads = threads,
        polypolish_careful = polypolish_careful,
    )
    if isdir(polishing_plan.polishing_dir) &&
       !isempty(readdir(polishing_plan.polishing_dir))
        throw(
            ArgumentError(
                "Autocycler polishing directory must be empty: " *
                "$(polishing_plan.polishing_dir)",
            ),
        )
    end
    mkpath(polishing_plan.polishing_dir)

    @info "Starting paired-short polishing with Polypolish and Pypolca..."
    final_assembly = try
        _execute_autocycler_steps(polishing_plan.steps; runner = runner)
        _require_nonempty_autocycler_file(
            polishing_plan.assembly,
            "Pypolca-polished Autocycler assembly",
        )
    catch
        if !keep_intermediates
            try
                _cleanup_autocycler_polishing_intermediates!(
                    polishing_plan.intermediate_files,
                )
            catch cleanup_error
                @warn "Autocycler failure cleanup could not remove intermediates" cleanup_error
            end
        end
        rethrow()
    end
    retained_intermediates = if keep_intermediates
        polishing_plan.intermediate_files
    else
        _cleanup_autocycler_polishing_intermediates!(
            polishing_plan.intermediate_files,
        )
        String[]
    end
    @info "Autocycler paired-short polishing complete" assembly = final_assembly

    return (
        outdir = normalized_out_dir,
        assembly = final_assembly,
        graph = autocycler_result.graph,
        autocycler_assembly = autocycler_result.assembly,
        polypolish_assembly = polishing_plan.polypolish_assembly,
        pypolca_report = polishing_plan.pypolca_report,
        intermediates = retained_intermediates,
        toolchain = autocycler_result.toolchain,
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assemble corrected long reads with Autocycler, then polish its consensus with a
corrected paired-short read set.

R1 and R2 are aligned separately with `bwa mem -a`, filtered as a pair, and
polished with Polypolish. Pypolca then runs in `--careful` mode on the
Polypolish result. `polypolish_careful=true` is the conservative default and can
be disabled explicitly for sufficiently deep short-read data. The raw
Autocycler GFA and unpolished FASTA are preserved alongside both polishing
stages. Mate count/order/identifiers and all required environment packages are
validated before long-read assembly starts. Large SAM, filtered-SAM, and BWA
index intermediates are removed after success unless `keep_intermediates=true`.
Route-owned BWA/SAM intermediates are also cleaned after a failed polishing run
unless explicit retention was requested; diagnostic assembly artifacts remain.

# Keywords
- `long_reads::AbstractString`: Nonempty long-read FASTQ.
- `short_reads_1::AbstractString`: Corrected paired-short R1 FASTQ.
- `short_reads_2::AbstractString`: Corrected paired-short R2 FASTQ.
- `out_dir::AbstractString`: Empty isolated working/output directory.
- `threads::Integer`: Threads per assembly/alignment job.
- `jobs::Integer`: Number of simultaneous Autocycler assembler jobs.
- `read_type::AbstractString`: Exact Autocycler long-read chemistry.
- `polypolish_careful::Bool`: Enable conservative Polypolish filtering.
- `keep_intermediates::Bool`: Retain BWA indices and SAM intermediates.

# Returns
A named tuple with final `assembly`, raw `graph`, `autocycler_assembly`,
`polypolish_assembly`, `pypolca_report`, `outdir`, and exact `toolchain`
provenance. `intermediates` lists files retained by an explicit
`keep_intermediates=true` request.
"""
function run_autocycler_polished(;
        long_reads::AbstractString,
        short_reads_1::AbstractString,
        short_reads_2::AbstractString,
        out_dir::AbstractString,
        threads::Integer = max(Sys.CPU_THREADS, 1),
        jobs::Integer = 1,
        read_type::AbstractString = "ont_r10",
        polypolish_careful::Bool = true,
        keep_intermediates::Bool = false,
)::NamedTuple
    return _run_autocycler_polished(
        long_reads,
        short_reads_1,
        short_reads_2,
        out_dir;
        threads = threads,
        jobs = jobs,
        read_type = read_type,
        polypolish_careful = polypolish_careful,
        keep_intermediates = keep_intermediates,
    )
end
