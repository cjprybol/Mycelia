"""
Wrapper for Autocycler and a conservative paired-short polishing pipeline.

References:
- Autocycler: https://github.com/rrwick/Autocycler
- Polypolish: https://github.com/rrwick/Polypolish
- Pypolca: https://github.com/gbouras13/pypolca
"""

const AUTOCYCLER_ENV_NAME = "autocycler"
const AUTOCYCLER_SCRIPT_URL = "https://raw.githubusercontent.com/rrwick/Autocycler/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh"
const AUTOCYCLER_ENV_URL = "https://raw.githubusercontent.com/rrwick/Autocycler/main/pipelines/Conda_environment_file_by_Ryan_Wick/environment.yml"
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install Autocycler and the short-read polishing tools in a shared conda
environment.

The package-bundled `environment.yml` is authoritative because it extends the
upstream Autocycler environment with BWA, Polypolish, and Pypolca. `force=true`
recreates the environment and refreshes the upstream automation script without
overwriting that environment specification.

# Keywords
- `force::Bool=false`: Recreate the environment and re-download the script.
"""
function install_autocycler(; force::Bool = false)::String
    install_dir, script_path, env_file_path = _autocycler_paths()
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

    if !isfile(script_path) || filesize(script_path) == 0 || force
        @info "Downloading autocycler_full.sh script..."
        Downloads.download(AUTOCYCLER_SCRIPT_URL, script_path)
        Base.chmod(script_path, 0o755)
    end

    if !isfile(script_path) || filesize(script_path) == 0
        throw(ErrorException("Autocycler script installation failed: $(script_path)"))
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

function _prepare_autocycler_output_dir(out_dir::AbstractString)::String
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

    bwa_index = (
        name = :bwa_index,
        command = _autocycler_conda_command(
            String["bwa", "index", normalized_assembly],
            polishing_dir;
            conda_runner = conda_runner,
        ),
        stdout = nothing,
        expected_outputs = String[],
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
        assembly = final_assembly,
    )
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
            message = sprint(showerror, error)
            throw(
                ErrorException(
                    "Autocycler workflow step $(step.name) failed. " *
                    "Recreate the tool environment with " *
                    "install_autocycler(force=true) if a command is missing. " *
                    "Cause: $(message)",
                ),
            )
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

function _ensure_autocycler_installed()::Nothing
    _, script_path, _ = _autocycler_paths()
    if !isfile(script_path) || filesize(script_path) == 0 ||
       !check_bioconda_env_is_installed(AUTOCYCLER_ENV_NAME)
        @warn "Autocycler is not fully installed. Installing now..."
        install_autocycler()
    end
    if !isfile(script_path) || filesize(script_path) == 0
        throw(ErrorException("Autocycler script is unavailable: $(script_path)"))
    end
    return nothing
end

function _run_autocycler(
        long_reads::AbstractString,
        out_dir::AbstractString;
        threads::Integer = max(Sys.CPU_THREADS, 1),
        jobs::Integer = 1,
        read_type::AbstractString = "ont_r10",
        dependency_checker::Function = _ensure_autocycler_installed,
        runner::Function = _default_autocycler_step_runner,
)::NamedTuple
    normalized_long_reads = _require_nonempty_autocycler_file(
        long_reads,
        "Long-read FASTQ",
    )
    normalized_read_type = _validate_autocycler_parameters(
        threads,
        jobs,
        read_type,
    )
    normalized_out_dir = _prepare_autocycler_output_dir(out_dir)
    dependency_checker()
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
A named tuple with `outdir`, `assembly`, and `graph` paths.
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
    normalized_read_type = _validate_autocycler_parameters(
        threads,
        jobs,
        read_type,
    )
    normalized_out_dir = _prepare_autocycler_output_dir(out_dir)
    dependency_checker()

    autocycler_result = _run_autocycler(
        normalized_long_reads,
        normalized_out_dir;
        threads = threads,
        jobs = jobs,
        read_type = normalized_read_type,
        dependency_checker = () -> nothing,
        runner = runner,
    )

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
    _execute_autocycler_steps(polishing_plan.steps; runner = runner)
    final_assembly = _require_nonempty_autocycler_file(
        polishing_plan.assembly,
        "Pypolca-polished Autocycler assembly",
    )
    @info "Autocycler paired-short polishing complete" assembly = final_assembly

    return (
        outdir = normalized_out_dir,
        assembly = final_assembly,
        graph = autocycler_result.graph,
        autocycler_assembly = autocycler_result.assembly,
        polypolish_assembly = polishing_plan.polypolish_assembly,
        pypolca_report = polishing_plan.pypolca_report,
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
stages.

# Returns
A named tuple with final `assembly`, raw `graph`, `autocycler_assembly`,
`polypolish_assembly`, `pypolca_report`, and `outdir` paths.
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
    )
end
