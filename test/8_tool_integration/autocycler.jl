# From the Mycelia base directory, run these tests with:
#
# ```bash
# julia --project=. --compiled-modules=no -e \
#   'include("test/8_tool_integration/autocycler.jl")'
# ```
#
# A real smoke test is opt-in because Autocycler is compute intensive:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true \
# MYCELIA_RUN_AUTOCYCLER_SMOKE=true \
# MYCELIA_AUTOCYCLER_LONG_READS=/path/to/long.fastq \
# MYCELIA_AUTOCYCLER_SHORT_READS_1=/path/to/R1.fastq \
# MYCELIA_AUTOCYCLER_SHORT_READS_2=/path/to/R2.fastq \
# julia --project=. --compiled-modules=no -e \
#   'include("test/8_tool_integration/autocycler.jl")'
# ```

import Logging
import Mycelia
import Test

function _autocycler_shell_quote(value::AbstractString)::String
    return "'" * replace(String(value), "'" => "'\"'\"'") * "'"
end

function _autocycler_public_test_conda_runner(
        path::AbstractString,
        environment_prefix::AbstractString,
)::String
    inventory = Mycelia.JSON.json(
        _autocycler_test_compatible_packages(),
    )
    quoted_prefix = _autocycler_shell_quote(environment_prefix)
    script = replace(
        raw"""#!/bin/sh
expected_prefix=__EXPECTED_PREFIX__
if [ "$1" = "list" ]; then
    [ "$2" = "-p" ] && [ "$3" = "$expected_prefix" ] || exit 61
    cat <<'MYCELIA_AUTOCYCLER_JSON'
__PACKAGE_INVENTORY__
MYCELIA_AUTOCYCLER_JSON
    exit 0
fi

[ "$1" = "run" ] || exit 62
shift
[ "$1" = "--live-stream" ] || exit 63
shift
[ "$1" = "-p" ] && [ "$2" = "$expected_prefix" ] || exit 64
shift 2

if [ "$1" = "bash" ]; then
    mkdir -p autocycler_out
    printf '>contig_1\nACGTACGT\n' > autocycler_out/consensus_assembly.fasta
    printf 'H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\tdp:f:20.0\n' > autocycler_out/consensus_assembly.gfa
elif [ "$1" = "bwa" ] && [ "$2" = "index" ]; then
    for extension in amb ann bwt pac sa; do
        printf 'index\n' > "$3.$extension"
    done
elif [ "$1" = "bwa" ] && [ "$2" = "mem" ]; then
    printf '@HD\tVN:1.6\n'
elif [ "$1" = "polypolish" ] && [ "$2" = "filter" ]; then
    shift 2
    while [ "$#" -gt 0 ]; do
        case "$1" in
            --out1) out1="$2"; shift 2 ;;
            --out2) out2="$2"; shift 2 ;;
            *) shift ;;
        esac
    done
    printf '@HD\tVN:1.6\n' > "$out1"
    printf '@HD\tVN:1.6\n' > "$out2"
elif [ "$1" = "polypolish" ] && [ "$2" = "polish" ]; then
    printf '>contig_1\nACGTACGT\n'
elif [ "$1" = "pypolca" ] && [ "$2" = "run" ]; then
    shift 2
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -o) output_dir="$2"; shift 2 ;;
            -p) output_prefix="$2"; shift 2 ;;
            *) shift ;;
        esac
    done
    mkdir -p "$output_dir"
    printf '>contig_1\nACGTACGT\n' > "$output_dir/${output_prefix}_corrected.fasta"
    printf 'corrected_bases\t0\n' > "$output_dir/${output_prefix}.report"
else
    exit 65
fi
""",
        "__EXPECTED_PREFIX__" => quoted_prefix,
        "__PACKAGE_INVENTORY__" => inventory,
    )
    normalized_path = abspath(String(path))
    mkpath(dirname(normalized_path))
    write(normalized_path, script)
    chmod(normalized_path, 0o755)
    return normalized_path
end

function _autocycler_test_runner!(step::NamedTuple)::Nothing
    for output_path in step.expected_outputs
        mkpath(dirname(output_path))
        if endswith(output_path, ".gfa")
            write(
                output_path,
                "H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\tdp:f:20.0\n",
            )
        elseif endswith(output_path, ".fasta") ||
               endswith(output_path, ".fa")
            write(output_path, ">contig_1\nACGTACGT\n")
        elseif endswith(output_path, ".sam")
            write(output_path, "@HD\tVN:1.6\n")
        else
            write(output_path, "test artifact\n")
        end
    end
    return nothing
end

function _autocycler_test_graph_only_runner!(step::NamedTuple)::Nothing
    for output_path in step.expected_outputs
        if endswith(output_path, ".gfa")
            mkpath(dirname(output_path))
            write(output_path, "H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\n")
        end
    end
    return nothing
end

function _autocycler_test_missing_output_runner(step::NamedTuple)::Nothing
    return nothing
end

function _autocycler_test_error(function_to_run::Function)::Exception
    try
        function_to_run()
    catch error
        return error
    end
    Base.error("Expected function to throw")
end

function _autocycler_test_compatible_packages(;
        versions::AbstractDict{<:AbstractString, <:AbstractString} =
            Dict{String, String}(),
        builds::AbstractDict{<:AbstractString, <:AbstractString} =
            Dict{String, String}(),
        channels::AbstractDict{<:AbstractString, <:AbstractString} =
            Dict{String, String}(),
        omitted::AbstractSet{<:AbstractString} = Set{String}(),
)::Vector{NamedTuple}
    packages = NamedTuple[
        (
            name = specification.name,
            version = get(
                versions,
                specification.name,
                specification.constraint == :present ?
                "present" : specification.version,
            ),
            build = get(builds, specification.name, "test_0"),
            channel = get(channels, specification.name, "conda-forge"),
        ) for specification in Mycelia.AUTOCYCLER_REQUIRED_PACKAGE_SPECS if
        !(specification.name in omitted)
    ]
    if !("python" in omitted)
        push!(packages, (
            name = "python",
            version = get(versions, "python", "3.11.9"),
            build = get(builds, "python", "h955ad1f_0_cpython"),
            channel = get(channels, "python", "conda-forge"),
        ))
    end
    return packages
end

function _autocycler_test_toolchain(;
        versions::AbstractDict{<:AbstractString, <:AbstractString} =
            Dict{String, String}(),
        builds::AbstractDict{<:AbstractString, <:AbstractString} =
            Dict{String, String}(),
        channels::AbstractDict{<:AbstractString, <:AbstractString} =
            Dict{String, String}(),
        omitted::AbstractSet{<:AbstractString} = Set{String}(),
)::Dict{String, Any}
    return Mycelia._autocycler_toolchain_metadata(
        _autocycler_test_compatible_packages(;
            versions,
            builds,
            channels,
            omitted,
        ),
    )
end

function _autocycler_smoke_env_enabled(
        environment::AbstractDict,
        name::AbstractString,
)::Bool
    value = String(get(environment, name, "false"))
    return lowercase(strip(value)) == "true"
end

function _autocycler_real_smoke_enabled(
        environment::AbstractDict,
)::Bool
    external_enabled =
        _autocycler_smoke_env_enabled(environment, "MYCELIA_RUN_ALL") ||
        _autocycler_smoke_env_enabled(environment, "MYCELIA_RUN_EXTERNAL")
    smoke_enabled = _autocycler_smoke_env_enabled(
        environment,
        "MYCELIA_RUN_AUTOCYCLER_SMOKE",
    )
    if smoke_enabled && !external_enabled
        throw(ArgumentError(
            "MYCELIA_RUN_AUTOCYCLER_SMOKE=true also requires " *
            "MYCELIA_RUN_EXTERNAL=true (or MYCELIA_RUN_ALL=true).",
        ))
    end
    return smoke_enabled
end

function _autocycler_real_smoke_inputs(
        environment::AbstractDict,
)::NamedTuple
    long_reads = String(get(
        environment,
        "MYCELIA_AUTOCYCLER_LONG_READS",
        "",
    ))
    short_reads_1 = String(get(
        environment,
        "MYCELIA_AUTOCYCLER_SHORT_READS_1",
        "",
    ))
    short_reads_2 = String(get(
        environment,
        "MYCELIA_AUTOCYCLER_SHORT_READS_2",
        "",
    ))
    read_type = String(get(
        environment,
        "MYCELIA_AUTOCYCLER_READ_TYPE",
        "ont_r10",
    ))
    isempty(long_reads) && throw(
        ArgumentError(
            "MYCELIA_RUN_AUTOCYCLER_SMOKE=true requires " *
            "MYCELIA_AUTOCYCLER_LONG_READS.",
        ),
    )
    xor(isempty(short_reads_1), isempty(short_reads_2)) && throw(
        ArgumentError(
            "Set both MYCELIA_AUTOCYCLER_SHORT_READS_1 and " *
            "MYCELIA_AUTOCYCLER_SHORT_READS_2, or leave both unset.",
        ),
    )
    return (; long_reads, short_reads_1, short_reads_2, read_type)
end

run_autocycler_smoke = _autocycler_real_smoke_enabled(ENV)
autocycler_smoke_inputs = run_autocycler_smoke ?
                          _autocycler_real_smoke_inputs(ENV) : nothing

Test.@testset "Autocycler wrapper" begin
    Test.@testset "Constants, paths, and bundled environment" begin
        Test.@test Mycelia.AUTOCYCLER_VERSION == "0.5.2"
        Test.@test Mycelia.AUTOCYCLER_ENVIRONMENT_SPEC_SHA256 ==
                   "d6aef758986db23c453203f067a23be6f29b690b536bc24b283b5a6e913624ef"
        Test.@test Mycelia.AUTOCYCLER_ENV_NAME ==
                   "autocycler-0.5.2-d6aef758986db23c"
        Test.@test Mycelia.AUTOCYCLER_MAX_ASSEMBLY_THREADS == 128
        Test.@test Mycelia.AUTOCYCLER_INSTALL_LOCK_STALE_SECONDS > 0
        Test.@test Mycelia.AUTOCYCLER_OUTPUT_LOCK_STALE_SECONDS > 0
        Test.@test Mycelia._effective_autocycler_assembly_threads(64) == 64
        Test.@test Mycelia._effective_autocycler_assembly_threads(256) == 128
        Test.@test Mycelia.AUTOCYCLER_READ_TYPES == (
            "ont_r9",
            "ont_r10",
            "pacbio_clr",
            "pacbio_hifi",
        )
        Test.@test !isempty(Mycelia.AUTOCYCLER_SCRIPT_URL)
        Test.@test occursin(
            Mycelia.AUTOCYCLER_SCRIPT_REVISION,
            Mycelia.AUTOCYCLER_SCRIPT_URL,
        )
        Test.@test length(Mycelia.AUTOCYCLER_SCRIPT_SHA256) == 64
        install_dir, script_path, env_file_path = Mycelia._autocycler_paths()
        Test.@test endswith(script_path, "autocycler_full.sh")
        Test.@test endswith(env_file_path, "environment.yml")
        Test.@test endswith(install_dir, joinpath("deps", "autocycler"))
        Test.@test occursin(
            Mycelia.AUTOCYCLER_ENV_NAME,
            basename(Mycelia._autocycler_install_lock_path()),
        )
        shared_runner = "/opt/shared-conda/bin/conda"
        shared_prefix = joinpath(
            "/opt/shared-conda",
            "envs",
            Mycelia.AUTOCYCLER_ENV_NAME,
        )
        Test.@test Mycelia._autocycler_environment_prefix(shared_runner) ==
                   shared_prefix
        Test.@test Mycelia._autocycler_install_lock_path(shared_runner) ==
                   joinpath(
            "/opt/shared-conda",
            ".mycelia-locks",
            "$(Mycelia.AUTOCYCLER_ENV_NAME).pid",
        )
        Test.@test !occursin(
            first(Base.DEPOT_PATH),
            Mycelia._autocycler_install_lock_path(shared_runner),
        )
        Test.@test Mycelia._require_verified_autocycler_environment_spec(
            env_file_path,
        ) == abspath(env_file_path)

        mktempdir() do temp_dir
            physical_root = joinpath(temp_dir, "physical-conda")
            alias_root = joinpath(temp_dir, "alias-conda")
            physical_runner = joinpath(physical_root, "bin", "conda")
            mkpath(dirname(physical_runner))
            write(physical_runner, "synthetic runner\n")
            symlink(physical_root, alias_root)
            alias_runner = joinpath(alias_root, "bin", "conda")
            Test.@test Mycelia._canonical_autocycler_conda_runner(
                alias_runner,
            ) == realpath(physical_runner)
            Test.@test Mycelia._autocycler_environment_prefix(alias_runner) ==
                       Mycelia._autocycler_environment_prefix(physical_runner)
            Test.@test Mycelia._autocycler_install_lock_path(alias_runner) ==
                       Mycelia._autocycler_install_lock_path(physical_runner)

            environment_file = joinpath(temp_dir, "environment.yml")
            write(environment_file, "name: synthetic\n")
            exact_prefix = joinpath(
                physical_root,
                "envs",
                Mycelia.AUTOCYCLER_ENV_NAME,
            )
            alias_prefix = joinpath(
                alias_root,
                "envs",
                Mycelia.AUTOCYCLER_ENV_NAME,
            )
            Test.@test Mycelia._canonical_autocycler_environment_prefix(
                alias_prefix,
            ) == exact_prefix
            Test.@test Mycelia._autocycler_install_lock_path_from_prefix(
                alias_prefix,
            ) == Mycelia._autocycler_install_lock_path_from_prefix(exact_prefix)
            alias_command = Mycelia._autocycler_conda_command(
                ["autocycler", "--version"],
                temp_dir;
                conda_runner = alias_runner,
                environment_prefix = alias_prefix,
            )
            Test.@test alias_command.exec[1] == realpath(physical_runner)
            Test.@test alias_command.exec[4:5] == ["-p", exact_prefix]
            mkpath(exact_prefix)
            Test.@test Mycelia._canonical_autocycler_environment_prefix(
                alias_prefix,
            ) == realpath(exact_prefix)
            commands = Cmd[]
            created_prefix = Mycelia._create_autocycler_environment_from_yaml(
                environment_file,
                Mycelia.AUTOCYCLER_ENV_NAME;
                conda_runner = alias_runner,
                environment_prefix = exact_prefix,
                command_runner = command -> begin
                    push!(commands, command)
                    return nothing
                end,
            )
            Test.@test created_prefix == exact_prefix
            Test.@test length(commands) == 2
            Test.@test commands[1].exec == String[
                realpath(physical_runner),
                "env",
                "create",
                "-f",
                abspath(environment_file),
                "-p",
                exact_prefix,
            ]
            Test.@test commands[2].exec == String[
                realpath(physical_runner),
                "clean",
                "--all",
                "-y",
            ]
        end

        environment_spec = read(env_file_path, String)
        Test.@test Mycelia._autocycler_sha256(env_file_path) ==
                   Mycelia.AUTOCYCLER_ENVIRONMENT_SPEC_SHA256
        Test.@test occursin("  - autocycler=0.5.2", environment_spec)
        dependency_section = split(environment_spec, "dependencies:"; limit = 2)[2]
        dependency_names = Set{String}()
        for line in split(dependency_section, '\n')
            dependency_match = match(r"^  - ([A-Za-z0-9_-]+)", line)
            dependency_match === nothing || push!(
                dependency_names,
                something(only(dependency_match.captures)),
            )
        end
        Test.@test dependency_names == Set(Mycelia.AUTOCYCLER_REQUIRED_PACKAGES)
        for specification in Mycelia.AUTOCYCLER_REQUIRED_PACKAGE_SPECS
            constraint = if specification.constraint == :exact
                "="
            elseif specification.constraint == :minimum
                ">="
            else
                ""
            end
            Test.@test occursin(
                "  - $(specification.name)$(constraint)$(specification.version)",
                environment_spec,
            )
        end
        Test.@test Mycelia._autocycler_script_is_verified(script_path)
        compatible_packages = _autocycler_test_compatible_packages()
        toolchain = Mycelia._autocycler_toolchain_metadata(compatible_packages)
        Test.@test Set(keys(toolchain)) == Set([
            "toolchain_schema",
            "autocycler_script_revision",
            "autocycler_script_sha256",
            "environment_name",
            "environment_spec_expected_sha256",
            "environment_spec_sha256",
            "inventory_schema",
            "package_inventory_sha256",
            "package_count",
            "packages",
        ])
        Test.@test toolchain["toolchain_schema"] ==
                   Mycelia.AUTOCYCLER_TOOLCHAIN_SCHEMA
        Test.@test toolchain["autocycler_script_revision"] ==
                   Mycelia.AUTOCYCLER_SCRIPT_REVISION
        Test.@test toolchain["autocycler_script_sha256"] ==
                   Mycelia.AUTOCYCLER_SCRIPT_SHA256
        Test.@test toolchain["environment_name"] ==
                   Mycelia.AUTOCYCLER_ENV_NAME
        Test.@test toolchain["environment_spec_expected_sha256"] ==
                   Mycelia.AUTOCYCLER_ENVIRONMENT_SPEC_SHA256
        Test.@test toolchain["environment_spec_sha256"] ==
                   Mycelia.AUTOCYCLER_ENVIRONMENT_SPEC_SHA256
        Test.@test toolchain["inventory_schema"] ==
                   Mycelia.AUTOCYCLER_PACKAGE_INVENTORY_SCHEMA
        Test.@test toolchain["package_inventory_sha256"] ==
                   Mycelia._autocycler_package_inventory_sha256(
            compatible_packages,
        )
        Test.@test toolchain["package_count"] == length(compatible_packages)
        python_packages = filter(
            package -> package["name"] == "python",
            toolchain["packages"],
        )
        Test.@test length(python_packages) == 1
        Test.@test only(python_packages)["version"] == "3.11.9"
        Test.@test only(python_packages)["build"] == "h955ad1f_0_cpython"
        Test.@test only(python_packages)["channel"] == "conda-forge"
        Test.@test Mycelia._require_autocycler_toolchain_provenance(
            toolchain,
        ) == toolchain
        reordered_toolchain = deepcopy(toolchain)
        reverse!(reordered_toolchain["packages"])
        normalized_reordered_snapshot =
            Mycelia._require_autocycler_toolchain_provenance(
                reordered_toolchain,
            )
        Test.@test normalized_reordered_snapshot == toolchain

        mutable_reported_snapshot = deepcopy(toolchain)
        detached_snapshot = Mycelia._require_autocycler_toolchain_provenance(
            mutable_reported_snapshot,
        )
        mutable_reported_snapshot["packages"][1]["build"] = "mutated-after-check"
        mutable_reported_snapshot["package_inventory_sha256"] = repeat("0", 64)
        Test.@test detached_snapshot == toolchain
        Test.@test detached_snapshot !== mutable_reported_snapshot
        Test.@test detached_snapshot["packages"] !==
                   mutable_reported_snapshot["packages"]
        normalized_packages =
            Mycelia._normalize_autocycler_package_inventory(
                toolchain["packages"],
            )
        Test.@test normalized_packages ==
                   Mycelia._normalize_autocycler_package_inventory(
            compatible_packages,
        )

        missing_inventory = copy(toolchain)
        delete!(missing_inventory, "packages")
        missing_inventory_error = _autocycler_test_error() do
            Mycelia._require_autocycler_toolchain_provenance(missing_inventory)
        end
        Test.@test missing_inventory_error isa ErrorException
        Test.@test occursin(
            "no realized package inventory",
            sprint(showerror, missing_inventory_error),
        )

        malformed_inventory = deepcopy(toolchain)
        delete!(first(malformed_inventory["packages"]), "build")
        malformed_inventory_error = _autocycler_test_error() do
            Mycelia._require_autocycler_toolchain_provenance(malformed_inventory)
        end
        Test.@test malformed_inventory_error isa ErrorException
        Test.@test occursin(
            "missing or empty build field",
            sprint(showerror, malformed_inventory_error),
        )

        missing_required_package = deepcopy(toolchain)
        filter!(
            package -> package["name"] != "bwa",
            missing_required_package["packages"],
        )
        missing_required_package["package_count"] =
            length(missing_required_package["packages"])
        missing_required_package["package_inventory_sha256"] =
            Mycelia._autocycler_package_inventory_sha256(
                missing_required_package["packages"],
            )
        missing_required_error = _autocycler_test_error() do
            Mycelia._require_autocycler_toolchain_provenance(
                missing_required_package,
            )
        end
        Test.@test missing_required_error isa ErrorException
        Test.@test occursin(
            "bwa is missing",
            sprint(showerror, missing_required_error),
        )

        for drift_field in ("build", "channel")
            drifted_inventory = deepcopy(toolchain)
            drifted_inventory["packages"][1][drift_field] *= "-drift"
            drift_error = _autocycler_test_error() do
                Mycelia._require_autocycler_toolchain_provenance(
                    drifted_inventory,
                )
            end
            Test.@test drift_error isa ErrorException
            Test.@test occursin(
                "inventory digest does not match",
                sprint(showerror, drift_error),
            )
        end

        for drift_field in ("build", "channel")
            transitive_drift = deepcopy(toolchain)
            python_package = only(filter(
                package -> package["name"] == "python",
                transitive_drift["packages"],
            ))
            python_package[drift_field] *= "-transitive-drift"
            transitive_drift_error = _autocycler_test_error() do
                Mycelia._require_autocycler_toolchain_provenance(
                    transitive_drift,
                )
            end
            Test.@test transitive_drift_error isa ErrorException
            Test.@test occursin(
                "inventory digest does not match",
                sprint(showerror, transitive_drift_error),
            )
        end

        build_drift_toolchain = _autocycler_test_toolchain(
            builds = Dict("bwa" => "different-build_1"),
        )
        build_drift_error = _autocycler_test_error() do
            Mycelia._require_unchanged_autocycler_toolchain(
                toolchain,
                build_drift_toolchain,
                "test boundary",
            )
        end
        Test.@test build_drift_error isa ErrorException
        Test.@test occursin(
            "changed across test boundary",
            sprint(showerror, build_drift_error),
        )

        channel_drift_toolchain = _autocycler_test_toolchain(
            channels = Dict("bwa" => "different-channel"),
        )
        channel_drift_error = _autocycler_test_error() do
            Mycelia._require_unchanged_autocycler_toolchain(
                toolchain,
                channel_drift_toolchain,
                "test boundary",
            )
        end
        Test.@test channel_drift_error isa ErrorException
        Test.@test occursin(
            "changed across test boundary",
            sprint(showerror, channel_drift_error),
        )
    end

    Test.@testset "Pinned script and stale-environment preflight" begin
        mktempdir() do temp_dir
            script_path = joinpath(temp_dir, "autocycler_full.sh")
            downloader = function (
                    url::AbstractString,
                    destination::AbstractString,
            )
                Test.@test url == Mycelia.AUTOCYCLER_SCRIPT_URL
                write(destination, "tampered script\n")
                return String(destination)
            end
            checksum_error = _autocycler_test_error() do
                Mycelia._install_verified_autocycler_script!(
                    script_path;
                    downloader,
                )
            end
            Test.@test checksum_error isa ErrorException
            Test.@test occursin(
                "checksum mismatch",
                sprint(showerror, checksum_error),
            )
            Test.@test !isfile(script_path)
        end

        mktempdir() do temp_dir
            installed_script = joinpath(temp_dir, "autocycler_full.sh")
            bundled_script = Mycelia._autocycler_paths()[2]
            downloader = function (
                    url::AbstractString,
                    destination::AbstractString,
            )
                Test.@test url == Mycelia.AUTOCYCLER_SCRIPT_URL
                cp(bundled_script, destination; force = true)
                return String(destination)
            end
            Test.@test Mycelia._install_verified_autocycler_script!(
                installed_script;
                downloader,
            ) == abspath(installed_script)
            Test.@test Mycelia._autocycler_script_is_verified(installed_script)
        end

        parsed_inventory = Mycelia._autocycler_environment_packages(;
            conda_runner = "/test/conda",
            environment_prefix = "/test/envs/autocycler-exact",
            command_reader = command -> begin
                Test.@test command.exec == String[
                    "/test/conda",
                    "list",
                    "-p",
                    "/test/envs/autocycler-exact",
                    "--json",
                ]
                return "[{\"name\":\"autocycler\",\"version\":\"0.5.2\"," *
                       "\"build_string\":\"test_0\"," *
                       "\"channel\":\"conda-forge\"}]"
            end,
        )
        Test.@test parsed_inventory == NamedTuple[(
            name = "autocycler",
            version = "0.5.2",
            build = "test_0",
            channel = "conda-forge",
        )]

        malformed_inventory_error = _autocycler_test_error() do
            Mycelia._autocycler_environment_packages(;
                command_reader = _command ->
                    "[{\"name\":\"autocycler\",\"version\":\"0.5.2\"," *
                    "\"build_string\":\"test_0\"}]",
            )
        end
        Test.@test malformed_inventory_error isa ErrorException
        Test.@test occursin(
            "missing or empty channel field",
            sprint(showerror, malformed_inventory_error),
        )

        nonarray_inventory_error = _autocycler_test_error() do
            Mycelia._autocycler_environment_packages(;
                command_reader = _command -> "{}",
            )
        end
        Test.@test nonarray_inventory_error isa ErrorException
        Test.@test occursin(
            "was not a JSON array",
            sprint(showerror, nonarray_inventory_error),
        )

        mktempdir() do temp_dir
            paths = (
                temp_dir,
                joinpath(temp_dir, "autocycler_full.sh"),
                joinpath(temp_dir, "environment.yml"),
            )
            lock_path = joinpath(temp_dir, "autocycler-install.pid")
            missing_environment_error = _autocycler_test_error() do
                Mycelia.install_autocycler(; paths, lock_path)
            end
            Test.@test missing_environment_error isa ErrorException
            Test.@test occursin(
                "environment file is missing or empty",
                sprint(showerror, missing_environment_error),
            )
            write(paths[3], "")
            empty_environment_error = _autocycler_test_error() do
                Mycelia.install_autocycler(; paths, lock_path)
            end
            Test.@test empty_environment_error isa ErrorException
            write(paths[3], "tampered environment\n")
            checksum_error = _autocycler_test_error() do
                Mycelia.install_autocycler(; paths, lock_path)
            end
            Test.@test checksum_error isa ErrorException
            Test.@test occursin(
                "environment specification checksum mismatch",
                sprint(showerror, checksum_error),
            )
            Test.@test !ispath(lock_path)
        end

        mktempdir() do temp_dir
            bundled_paths = Mycelia._autocycler_paths()
            paths = (
                temp_dir,
                joinpath(temp_dir, "autocycler_full.sh"),
                joinpath(temp_dir, "environment.yml"),
            )
            cp(bundled_paths[2], paths[2])
            cp(bundled_paths[3], paths[3])
            lock_path = joinpath(temp_dir, "autocycler-install.pid")
            environment_creations = Ref(0)
            script_downloads = Ref(0)
            force_error = _autocycler_test_error() do
                Mycelia.install_autocycler(;
                    force = true,
                    paths,
                    environment_checker = () -> begin
                        Test.@test isfile(lock_path)
                        return true
                    end,
                    environment_creator = function (
                            _environment_file::AbstractString,
                            _environment_name::AbstractString;
                            force::Bool = false,
                    )
                        environment_creations[] += 1
                        return force
                    end,
                    downloader = (_url, _destination) -> begin
                        script_downloads[] += 1
                        return nothing
                    end,
                    lock_path,
                )
            end
            Test.@test force_error isa ErrorException
            Test.@test occursin(
                "Refusing install_autocycler(force=true)",
                sprint(showerror, force_error),
            )
            Test.@test occursin(
                "immutable spec-hash-addressed environment",
                sprint(showerror, force_error),
            )
            Test.@test environment_creations[] == 0
            Test.@test script_downloads[] == 0
            Test.@test !ispath(lock_path)
        end

        mktempdir() do temp_dir
            bundled_paths = Mycelia._autocycler_paths()
            paths = (
                temp_dir,
                joinpath(temp_dir, "autocycler_full.sh"),
                joinpath(temp_dir, "environment.yml"),
            )
            cp(bundled_paths[3], paths[3])
            lock_path = joinpath(temp_dir, "autocycler-install.pid")
            environment_creator_forces = Bool[]
            package_inspections = Ref(0)
            environment_creator = function (
                    environment_file::AbstractString,
                    environment_name::AbstractString;
                    force::Bool = false,
            )
                Test.@test environment_file == abspath(paths[3])
                Test.@test environment_name == Mycelia.AUTOCYCLER_ENV_NAME
                Test.@test isfile(lock_path)
                push!(environment_creator_forces, force)
                return String(environment_name)
            end
            downloader = function (
                    url::AbstractString,
                    destination::AbstractString,
            )
                Test.@test url == Mycelia.AUTOCYCLER_SCRIPT_URL
                cp(bundled_paths[2], destination; force = true)
                return String(destination)
            end
            installed_script = Mycelia.install_autocycler(;
                force = true,
                paths,
                environment_checker = () -> begin
                    Test.@test isfile(lock_path)
                    return false
                end,
                environment_creator,
                package_inspector = () -> begin
                    Test.@test isfile(lock_path)
                    package_inspections[] += 1
                    return _autocycler_test_compatible_packages()
                end,
                downloader,
                lock_path,
            )
            Test.@test installed_script == paths[2]
            Test.@test environment_creator_forces == [false]
            Test.@test package_inspections[] == 1
            Test.@test Mycelia._autocycler_script_is_verified(installed_script)
            Test.@test !ispath(lock_path)
        end

        Test.@testset "Public installer validates packages under its lock" begin
            mktempdir() do temp_dir
                bundled_paths = Mycelia._autocycler_paths()
                paths = (
                    temp_dir,
                    joinpath(temp_dir, "autocycler_full.sh"),
                    joinpath(temp_dir, "environment.yml"),
                )
                cp(bundled_paths[2], paths[2])
                cp(bundled_paths[3], paths[3])
                lock_path = joinpath(temp_dir, "autocycler-install.pid")
                inspection_calls = Ref(0)
                installed_script = Mycelia.install_autocycler(;
                    paths,
                    environment_checker = () -> true,
                    environment_creator = (_args...; _kwargs...) -> error(
                        "existing environment must not be recreated",
                    ),
                    package_inspector = () -> begin
                        Test.@test isfile(lock_path)
                        inspection_calls[] += 1
                        return _autocycler_test_compatible_packages()
                    end,
                    downloader = (_args...) -> error(
                        "verified script must not be downloaded",
                    ),
                    lock_path,
                )
                Test.@test installed_script == paths[2]
                Test.@test inspection_calls[] == 1
                Test.@test !ispath(lock_path)

                incompatible_packages = _autocycler_test_compatible_packages(
                    versions = Dict("autocycler" => "0.5.1"),
                )
                incompatible_error = _autocycler_test_error() do
                    Mycelia.install_autocycler(;
                        paths,
                        environment_checker = () -> true,
                        package_inspector = () -> begin
                            Test.@test isfile(lock_path)
                            return incompatible_packages
                        end,
                        lock_path,
                    )
                end
                Test.@test incompatible_error isa ErrorException
                Test.@test occursin(
                    "autocycler must equal 0.5.2",
                    sprint(showerror, incompatible_error),
                )
                Test.@test !ispath(lock_path)
            end

            mktempdir() do temp_dir
                bundled_paths = Mycelia._autocycler_paths()
                paths = (
                    temp_dir,
                    joinpath(temp_dir, "autocycler_full.sh"),
                    joinpath(temp_dir, "environment.yml"),
                )
                cp(bundled_paths[3], paths[3])
                lock_path = joinpath(temp_dir, "autocycler-install.pid")
                environment_creations = Ref(0)
                package_inspections = Ref(0)
                created_environment_error = _autocycler_test_error() do
                    Mycelia.install_autocycler(;
                        paths,
                        environment_checker = () -> false,
                        environment_creator = function (
                                _environment_file::AbstractString,
                                _environment_name::AbstractString;
                                force::Bool = false,
                        )
                            Test.@test !force
                            Test.@test isfile(lock_path)
                            environment_creations[] += 1
                            return nothing
                        end,
                        package_inspector = () -> begin
                            Test.@test isfile(lock_path)
                            package_inspections[] += 1
                            return _autocycler_test_compatible_packages(
                                omitted = Set(["bwa"]),
                            )
                        end,
                        downloader = function (
                                _url::AbstractString,
                                destination::AbstractString,
                        )
                            cp(bundled_paths[2], destination; force = true)
                            return String(destination)
                        end,
                        lock_path,
                    )
                end
                Test.@test created_environment_error isa ErrorException
                Test.@test occursin(
                    "missing or incompatible required packages",
                    sprint(showerror, created_environment_error),
                )
                Test.@test environment_creations[] == 1
                Test.@test package_inspections[] == 1
                Test.@test !ispath(lock_path)
            end
        end

        inspection_calls = Ref(0)
        package_inspector = function ()
            inspection_calls[] += 1
            return _autocycler_test_compatible_packages(
                versions = Dict("autocycler" => "0.5.1"),
                omitted = Set(["canu"]),
            )
        end
        package_error = _autocycler_test_error() do
            Mycelia._ensure_autocycler_packages!(package_inspector)
        end
        Test.@test package_error isa ErrorException
        Test.@test inspection_calls[] == 1
        Test.@test occursin(
            "immutable spec-hash-addressed environment",
            sprint(showerror, package_error),
        )
        Test.@test occursin(
            "autocycler must equal 0.5.2",
            sprint(showerror, package_error),
        )
        Test.@test occursin(
            "canu is missing",
            sprint(showerror, package_error),
        )

        too_old = _autocycler_test_compatible_packages(
            versions = Dict("flye" => "2.9.5"),
        )
        Test.@test any(
            issue -> occursin("flye must be at least 2.9.6", issue),
            Mycelia._autocycler_package_issues(too_old),
        )
        necat_without_update = _autocycler_test_compatible_packages(
            versions = Dict("necat" => "0.0.1"),
        )
        Test.@test any(
            issue -> occursin("necat must be at least 0.0.1_update20200803", issue),
            Mycelia._autocycler_package_issues(necat_without_update),
        )
        Test.@test Mycelia._autocycler_version_at_least(
            "0.0.1_update20210101",
            "0.0.1_update20200803",
        )

        always_stale = () -> _autocycler_test_compatible_packages(
            omitted = Set([
                specification.name for
                specification in Mycelia.AUTOCYCLER_REQUIRED_PACKAGE_SPECS if
                !(specification.name in ("autocycler", "bwa"))
            ]),
        )
        stale_error = _autocycler_test_error() do
            Mycelia._ensure_autocycler_packages!(always_stale)
        end
        Test.@test stale_error isa ErrorException
        Test.@test occursin(
            "missing or incompatible required packages",
            sprint(showerror, stale_error),
        )
        Test.@test occursin(
            "remove it manually before reinstalling",
            sprint(showerror, stale_error),
        )

        mktempdir() do temp_dir
            lock_path = joinpath(temp_dir, "autocycler-install.pid")
            lock_observed = Ref(false)
            locked_result = Mycelia._with_autocycler_install_lock(
                lock_path,
            ) do
                lock_observed[] = isfile(lock_path)
                return Dict{String, Any}("locked" => true)
            end
            Test.@test lock_observed[]
            Test.@test locked_result == Dict{String, Any}("locked" => true)
            Test.@test !ispath(lock_path)

            observed_stale_age = Ref(0.0)
            pidlock_runner = function (
                    action::Function,
                    observed_path::AbstractString;
                    stale_age::Real,
                    poll_interval::Real,
            )
                Test.@test observed_path == abspath(lock_path)
                Test.@test poll_interval > 0
                observed_stale_age[] = Float64(stale_age)
                write(observed_path, "synthetic bounded lock\n")
                try
                    return action()
                finally
                    rm(observed_path; force = true)
                end
            end
            bounded_result = Mycelia._with_autocycler_install_lock(
                () -> :bounded,
                lock_path;
                pidlock_runner,
            )
            Test.@test bounded_result == :bounded
            Test.@test observed_stale_age[] ==
                       Mycelia.AUTOCYCLER_INSTALL_LOCK_STALE_SECONDS
            Test.@test observed_stale_age[] > 0
            Test.@test !ispath(lock_path)

            lock_active = Ref(false)
            lock_events = Symbol[]
            lock_runner = function (
                    action::Function,
                    observed_path::AbstractString,
            )
                Test.@test observed_path == lock_path
                Test.@test !lock_active[]
                lock_active[] = true
                push!(lock_events, :entered)
                try
                    return action()
                finally
                    push!(lock_events, :exited)
                    lock_active[] = false
                end
            end
            environment_checks = Ref(0)
            environment_checker = function ()
                Test.@test lock_active[]
                environment_checks[] += 1
                return environment_checks[] >= 2
            end
            locked_inspections = Ref(0)
            locked_package_inspector = function ()
                Test.@test lock_active[]
                locked_inspections[] += 1
                return _autocycler_test_compatible_packages()
            end
            locked_installs = Bool[]
            locked_installer = function (; force::Bool = false)
                Test.@test lock_active[]
                push!(locked_installs, force)
                return nothing
            end
            locked_toolchain = Mycelia._ensure_autocycler_installed(;
                package_inspector = locked_package_inspector,
                installer = locked_installer,
                environment_checker,
                lock_path,
                lock_runner,
            )
            Test.@test lock_events == [:entered, :exited]
            Test.@test environment_checks[] == 2
            Test.@test locked_inspections[] == 1
            Test.@test locked_installs == [false]
            Test.@test Mycelia._normalize_autocycler_package_inventory(
                locked_toolchain["packages"],
            ) == Mycelia._normalize_autocycler_package_inventory(
                _autocycler_test_compatible_packages(),
            )
        end
    end

    Test.@testset "Exact long-read-only upstream command" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            out_dir = joinpath(temp_dir, "results")
            script_path = joinpath(temp_dir, "autocycler_full.sh")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(script_path, "#!/usr/bin/env bash\n")
            environment_prefix = joinpath(
                temp_dir,
                "conda",
                "envs",
                Mycelia.AUTOCYCLER_ENV_NAME,
            )

            plan = Mycelia._autocycler_command_plan(
                long_reads,
                out_dir;
                threads = 8,
                jobs = 2,
                read_type = "pacbio_hifi",
                script_path = script_path,
                conda_runner = "/test/conda",
                environment_prefix,
            )
            command = only(plan.steps).command
            Test.@test command.exec == String[
                "/test/conda",
                "run",
                "--live-stream",
                "-p",
                environment_prefix,
                "bash",
                abspath(script_path),
                abspath(long_reads),
                "8",
                "2",
                "pacbio_hifi",
            ]
            Test.@test command.dir == abspath(out_dir)
            Test.@test plan.graph == joinpath(
                abspath(out_dir),
                "autocycler_out",
                "consensus_assembly.gfa",
            )
            Test.@test plan.assembly == joinpath(
                abspath(out_dir),
                "autocycler_out",
                "consensus_assembly.fasta",
            )
            Test.@test only(plan.steps).expected_outputs == String[
                plan.graph,
                plan.assembly,
            ]
            Test.@test plan.requested_threads == 8
            Test.@test plan.autocycler_assembly_threads == 8

            capped_plan = Mycelia._autocycler_command_plan(
                long_reads,
                out_dir;
                threads = 256,
                jobs = 2,
                read_type = "pacbio_hifi",
                script_path = script_path,
                conda_runner = "/test/conda",
            )
            Test.@test capped_plan.requested_threads == 256
            Test.@test capped_plan.autocycler_assembly_threads == 128
            Test.@test only(capped_plan.steps).command.exec[end - 2] == "128"
            Test.@test !ispath(out_dir)
        end
    end

    Test.@testset "Parameter and output validation fails loudly" begin
        invalid_type = _autocycler_test_error() do
            Mycelia._autocycler_command_plan(
                "reads.fastq",
                "out";
                read_type = "illumina",
            )
        end
        Test.@test invalid_type isa ArgumentError
        Test.@test occursin("read_type must be one of", sprint(showerror, invalid_type))
        Test.@test occursin("illumina", sprint(showerror, invalid_type))

        invalid_threads = _autocycler_test_error() do
            Mycelia._autocycler_command_plan(
                "reads.fastq",
                "out";
                threads = 0,
            )
        end
        Test.@test invalid_threads isa ArgumentError
        Test.@test occursin("threads must be positive", sprint(showerror, invalid_threads))

        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            malformed_long_reads = joinpath(temp_dir, "not-fastq.fasta")
            write(malformed_long_reads, ">long\nACGT\n")
            cheap_dependency_checks = Ref(0)
            cheap_runner_calls = Ref(0)
            cheap_runner = function (step::NamedTuple)
                cheap_runner_calls[] += 1
                return nothing
            end

            blank_output_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    malformed_long_reads,
                    "   ";
                    dependency_checker = () -> begin
                        cheap_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = cheap_runner,
                )
            end
            Test.@test blank_output_error isa ArgumentError
            Test.@test occursin(
                "non-blank path",
                sprint(showerror, blank_output_error),
            )

            blocking_ancestor = joinpath(temp_dir, "blocking-ancestor")
            write(blocking_ancestor, "not a directory\n")
            ancestor_output_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    malformed_long_reads,
                    joinpath(blocking_ancestor, "nested", "results");
                    dependency_checker = () -> begin
                        cheap_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = cheap_runner,
                )
            end
            Test.@test ancestor_output_error isa ArgumentError
            Test.@test occursin(
                "non-directory existing ancestor",
                sprint(showerror, ancestor_output_error),
            )

            dangling_output = joinpath(temp_dir, "dangling-output")
            symlink(joinpath(temp_dir, "missing-output-target"), dangling_output)
            dangling_output_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    malformed_long_reads,
                    dangling_output;
                    dependency_checker = () -> begin
                        cheap_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = cheap_runner,
                )
            end
            Test.@test dangling_output_error isa ArgumentError
            Test.@test occursin(
                "must not be a symbolic link",
                sprint(showerror, dangling_output_error),
            )

            dangling_ancestor = joinpath(temp_dir, "dangling-ancestor")
            symlink(joinpath(temp_dir, "missing-ancestor-target"), dangling_ancestor)
            dangling_ancestor_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    malformed_long_reads,
                    joinpath(dangling_ancestor, "nested", "results");
                    dependency_checker = () -> begin
                        cheap_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = cheap_runner,
                )
            end
            Test.@test dangling_ancestor_error isa ArgumentError
            Test.@test occursin(
                "non-directory existing ancestor",
                sprint(showerror, dangling_ancestor_error),
            )

            symlink_target = joinpath(temp_dir, "symlink-target")
            symlink_ancestor = joinpath(temp_dir, "symlink-ancestor")
            mkpath(symlink_target)
            symlink(symlink_target, symlink_ancestor)
            aliased_output = joinpath(symlink_ancestor, "nested", "results")
            Test.@test Mycelia._validate_autocycler_output_dir(aliased_output) ==
                       abspath(aliased_output)
            expected_alias_lock = joinpath(
                realpath(symlink_target),
                "nested",
                ".results.mycelia-output-root.pid",
            )
            Test.@test Mycelia._autocycler_output_lock_path(aliased_output) ==
                       expected_alias_lock

            retarget_target = joinpath(temp_dir, "retarget-target")
            mkpath(retarget_target)
            expected_canonical_output = joinpath(
                realpath(symlink_target),
                "nested",
                "results",
            )
            retargeted_output = joinpath(
                realpath(retarget_target),
                "nested",
                "results",
            )
            retargeting_pidlock_runner = function (
                    action::Function,
                    lock_path::AbstractString;
                    stale_age::Real,
                    poll_interval::Real,
            )
                Test.@test lock_path == expected_alias_lock
                Test.@test stale_age > 0
                Test.@test poll_interval > 0
                rm(symlink_ancestor)
                symlink(retarget_target, symlink_ancestor)
                return action()
            end
            retargeting_output_lock_runner = function (
                    action::Function,
                    output_path::AbstractString,
            )
                return Mycelia._with_autocycler_output_lock(
                    action,
                    output_path;
                    pidlock_runner = retargeting_pidlock_runner,
                )
            end
            retarget_result = Mycelia._run_autocycler(
                long_reads,
                aliased_output;
                dependency_checker = _autocycler_test_toolchain,
                runner = _autocycler_test_runner!,
                output_lock_runner = retargeting_output_lock_runner,
            )
            Test.@test retarget_result.outdir == expected_canonical_output
            Test.@test isfile(retarget_result.assembly)
            Test.@test isfile(retarget_result.graph)
            Test.@test !ispath(retargeted_output)

            cheap_parameter_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    malformed_long_reads,
                    joinpath(temp_dir, "invalid-parameters");
                    threads = 0,
                    dependency_checker = () -> begin
                        cheap_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = cheap_runner,
                )
            end
            Test.@test cheap_parameter_error isa ArgumentError
            Test.@test occursin(
                "threads must be positive",
                sprint(showerror, cheap_parameter_error),
            )
            Test.@test cheap_dependency_checks[] == 0
            Test.@test cheap_runner_calls[] == 0

            nonempty_out_dir = joinpath(temp_dir, "nonempty")
            mkpath(nonempty_out_dir)
            write(joinpath(nonempty_out_dir, "owned.txt"), "keep\n")
            nonempty_dependency_checks = Ref(0)

            nonempty_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    nonempty_out_dir;
                    dependency_checker = () -> begin
                        nonempty_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = _autocycler_test_runner!,
                )
            end
            Test.@test nonempty_error isa ArgumentError
            Test.@test occursin(
                "output directory must be empty",
                sprint(showerror, nonempty_error),
            )
            Test.@test nonempty_dependency_checks[] == 0
            Test.@test isfile(joinpath(nonempty_out_dir, "owned.txt"))

            missing_output_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "missing-output");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = _autocycler_test_missing_output_runner,
                )
            end
            Test.@test missing_output_error isa ErrorException
            Test.@test occursin(
                "did not create a nonempty artifact",
                sprint(showerror, missing_output_error),
            )
            Test.@test occursin(
                "consensus_assembly.gfa",
                sprint(showerror, missing_output_error),
            )

            missing_fasta_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "missing-fasta");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = _autocycler_test_graph_only_runner!,
                )
            end
            Test.@test missing_fasta_error isa ErrorException
            Test.@test occursin(
                "did not create a nonempty artifact",
                sprint(showerror, missing_fasta_error),
            )
            Test.@test occursin(
                "consensus_assembly.fasta",
                sprint(showerror, missing_fasta_error),
            )
        end
    end

    Test.@testset "Step execution preserves exception semantics" begin
        step = (
            name = :synthetic,
            command = `true`,
            stdout = nothing,
            expected_outputs = String[],
        )

        interrupted = _autocycler_test_error() do
            Mycelia._execute_autocycler_steps(
                (step,);
                runner = _ -> throw(InterruptException()),
            )
        end
        Test.@test interrupted isa InterruptException

        runner_bug = ErrorException("synthetic runner bug")
        preserved_bug = _autocycler_test_error() do
            Mycelia._execute_autocycler_steps(
                (step,);
                runner = _ -> throw(runner_bug),
            )
        end
        Test.@test preserved_bug === runner_bug

        process_error = _autocycler_test_error() do
            Mycelia._execute_autocycler_steps(
                (step,);
                runner = _ -> Base.run(`false`),
            )
        end
        Test.@test process_error isa ErrorException
        Test.@test occursin(
            "workflow command synthetic failed",
            sprint(showerror, process_error),
        )
        Test.@test occursin(
            "spec-hash-addressed environment",
            sprint(showerror, process_error),
        )
        Test.@test occursin(
            "remove it manually before reinstalling",
            sprint(showerror, process_error),
        )
        Test.@test !occursin(
            "install_autocycler(force=true)",
            sprint(showerror, process_error),
        )

        mktempdir() do temp_dir
            intermediate_target = joinpath(temp_dir, "intermediate-target.sam")
            intermediate_link = joinpath(temp_dir, "intermediate.sam")
            symlink_step = (
                name = :synthetic_intermediate,
                command = `true`,
                stdout = nothing,
                expected_outputs = String[intermediate_link],
            )
            symlink_error = _autocycler_test_error() do
                Mycelia._execute_autocycler_steps(
                    (symlink_step,);
                    runner = _ -> begin
                        write(intermediate_target, "owned\n")
                        symlink(intermediate_target, intermediate_link)
                        return nothing
                    end,
                )
            end
            Test.@test symlink_error isa ErrorException
            Test.@test occursin(
                "symlinked artifact",
                sprint(showerror, symlink_error),
            )

            intermediate = joinpath(temp_dir, "intermediate.sam")
            rm(intermediate; force = true)
            write(intermediate, "owned\n")
            cleanup_interrupt = _autocycler_test_error() do
                Mycelia._cleanup_autocycler_polishing_intermediates!(
                    [intermediate];
                    workflow_root = temp_dir,
                    remover = (_path::AbstractString) ->
                        throw(InterruptException()),
                )
            end
            Test.@test cleanup_interrupt isa InterruptException
            Test.@test isfile(intermediate)

            dangling_intermediate = joinpath(temp_dir, "dangling.sam")
            symlink(
                joinpath(temp_dir, "missing-cleanup-target"),
                dangling_intermediate,
            )
            dangling_retained = Test.@test_logs (
                :warn,
                r"could not remove a polishing intermediate",
            ) min_level=Logging.Warn begin
                Mycelia._cleanup_autocycler_polishing_intermediates!(
                    [dangling_intermediate];
                    workflow_root = temp_dir,
                    remover = (_path::AbstractString) -> nothing,
                )
            end
            Test.@test dangling_retained == [dangling_intermediate]
            Test.@test islink(dangling_intermediate)
        end
    end

    Test.@testset "Autocycler returns authoritative FASTA and GFA" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            dependency_checks = Ref(0)
            expected_toolchain = _autocycler_test_toolchain()
            out_dir = joinpath(temp_dir, "autocycler-run")
            conda_runner = joinpath(temp_dir, "exact-conda", "bin", "conda")
            mkpath(dirname(conda_runner))
            write(conda_runner, "synthetic runner\n")
            environment_prefix = joinpath(
                temp_dir,
                "exact-conda",
                "envs",
                Mycelia.AUTOCYCLER_ENV_NAME,
            )
            output_lock_path = Mycelia._autocycler_output_lock_path(out_dir)
            environment_lock_path = joinpath(temp_dir, "environment.pid")
            environment_lock_active = Ref(false)
            environment_lock_runner = function (
                    action::Function,
                    observed_path::AbstractString,
            )
                Test.@test observed_path == environment_lock_path
                Test.@test !environment_lock_active[]
                environment_lock_active[] = true
                try
                    return action()
                finally
                    environment_lock_active[] = false
                end
            end
            result = Mycelia._run_autocycler(
                long_reads,
                out_dir;
                threads = 4,
                jobs = 3,
                read_type = "ont_r9",
                conda_runner,
                environment_prefix,
                dependency_checker = () -> begin
                    Test.@test isfile(output_lock_path)
                    Test.@test environment_lock_active[]
                    dependency_checks[] += 1
                    return expected_toolchain
                end,
                runner = (step::NamedTuple) -> begin
                    Test.@test isfile(output_lock_path)
                    Test.@test environment_lock_active[]
                    Test.@test first(step.command.exec) == realpath(conda_runner)
                    prefix_index = findfirst(==("-p"), step.command.exec)
                    Test.@test prefix_index !== nothing
                    Test.@test step.command.exec[prefix_index + 1] ==
                               environment_prefix
                    return _autocycler_test_runner!(step)
                end,
                environment_lock_path,
                environment_lock_runner,
            )
            Test.@test dependency_checks[] == 2
            Test.@test result.toolchain == expected_toolchain
            Test.@test isfile(result.graph)
            Test.@test isfile(result.assembly)
            Test.@test read(result.assembly, String) ==
                       ">contig_1\nACGTACGT\n"
            Test.@test read(result.graph, String) ==
                       "H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\tdp:f:20.0\n"
            Test.@test result.outdir == abspath(joinpath(temp_dir, "autocycler-run"))
            Test.@test result.requested_threads == 4
            Test.@test result.autocycler_assembly_threads == 4
            Test.@test !ispath(output_lock_path)

            mismatched_companion_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    fasta_path = only(filter(
                        path -> endswith(path, ".fasta"),
                        step.expected_outputs,
                    ))
                    write(fasta_path, ">contig_1\nTTAACCGG\n")
                end
                return nothing
            end
            mismatched_companion_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "mismatched-companion");
                    dependency_checker = () -> expected_toolchain,
                    runner = mismatched_companion_runner,
                )
            end
            Test.@test mismatched_companion_error isa ErrorException
            Test.@test occursin(
                "contain different sequences",
                sprint(showerror, mismatched_companion_error),
            )

            mutation_checks = Ref(0)
            mutated_toolchain = _autocycler_test_toolchain(
                builds = Dict("bwa" => "mutated-build_1"),
            )
            mutation_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "mutated-toolchain-run");
                    dependency_checker = () -> begin
                        Test.@test environment_lock_active[]
                        mutation_checks[] += 1
                        return mutation_checks[] == 1 ?
                               expected_toolchain : mutated_toolchain
                    end,
                    runner = (step::NamedTuple) -> begin
                        Test.@test environment_lock_active[]
                        return _autocycler_test_runner!(step)
                    end,
                    environment_lock_path,
                    environment_lock_runner,
                )
            end
            Test.@test mutation_error isa ErrorException
            Test.@test mutation_checks[] == 2
            Test.@test occursin(
                "changed across long-read assembly",
                sprint(showerror, mutation_error),
            )

            malformed_fasta_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    fasta_path = only(filter(
                        path -> endswith(path, ".fasta"),
                        step.expected_outputs,
                    ))
                    write(fasta_path, "not FASTA\n")
                end
                return nothing
            end
            malformed_fasta_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "malformed-fasta");
                    dependency_checker = () -> expected_toolchain,
                    runner = malformed_fasta_runner,
                )
            end
            Test.@test malformed_fasta_error isa ErrorException
            Test.@test occursin(
                "not valid FASTA",
                sprint(showerror, malformed_fasta_error),
            )

            invalid_dna_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    fasta_path = only(filter(
                        path -> endswith(path, ".fasta"),
                        step.expected_outputs,
                    ))
                    write(fasta_path, ">invalid_dna\nACGTZ\n")
                end
                return nothing
            end
            invalid_dna_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "invalid-dna-fasta");
                    dependency_checker = () -> expected_toolchain,
                    runner = invalid_dna_runner,
                )
            end
            Test.@test invalid_dna_error isa ErrorException
            Test.@test occursin(
                "contains invalid DNA",
                sprint(showerror, invalid_dna_error),
            )

            semantic_fasta_cases = (
                (
                    slug = "gapped-dna-fasta",
                    contents = ">gapped\nAC-GT\n",
                    message = "contains invalid DNA",
                ),
                (
                    slug = "duplicate-fasta-identifier",
                    contents =
                        ">duplicate\nACGT\n>duplicate\nTGCA\n",
                    message = "duplicate FASTA identifier",
                ),
            )
            for semantic_case in semantic_fasta_cases
                semantic_runner = function (step::NamedTuple)
                    _autocycler_test_runner!(step)
                    if step.name == :autocycler
                        fasta_path = only(filter(
                            path -> endswith(path, ".fasta"),
                            step.expected_outputs,
                        ))
                        write(fasta_path, semantic_case.contents)
                    end
                    return nothing
                end
                semantic_error = _autocycler_test_error() do
                    Mycelia._run_autocycler(
                        long_reads,
                        joinpath(temp_dir, semantic_case.slug);
                        dependency_checker = () -> expected_toolchain,
                        runner = semantic_runner,
                    )
                end
                Test.@test semantic_error isa ErrorException
                Test.@test occursin(
                    semantic_case.message,
                    sprint(showerror, semantic_error),
                )
            end

            fasta_symlink_target = joinpath(temp_dir, "fasta-symlink-target")
            write(fasta_symlink_target, ">symlinked\nACGT\n")
            symlinked_fasta_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    fasta_path = only(filter(
                        path -> endswith(path, ".fasta"),
                        step.expected_outputs,
                    ))
                    rm(fasta_path)
                    symlink(fasta_symlink_target, fasta_path)
                end
                return nothing
            end
            symlinked_fasta_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "symlinked-fasta");
                    dependency_checker = () -> expected_toolchain,
                    runner = symlinked_fasta_runner,
                )
            end
            Test.@test symlinked_fasta_error isa ErrorException
            Test.@test occursin(
                "symlinked artifact",
                sprint(showerror, symlinked_fasta_error),
            )

            structured_gfa = joinpath(temp_dir, "structured.gfa")
            write(
                structured_gfa,
                "H\tVN:Z:1.0\n" *
                "S\tleft\tACGT\tLN:i:4\n" *
                "S\tright\tTGCA\n" *
                "L\tleft\t+\tright\t-\t1M\n" *
                "P\tprimary\tleft+,right-\t1M\n",
            )
            Test.@test Mycelia._require_valid_autocycler_gfa(
                structured_gfa,
                "Autocycler structured GFA",
            ) == structured_gfa

            comma_identifier_gfa = joinpath(temp_dir, "comma-identifier.gfa")
            write(
                comma_identifier_gfa,
                "S\tleft,half\tACGT\n" *
                "S\tright\tTGCA\n" *
                "P\tprimary\tleft,half+,right-\t1M\n",
            )
            Test.@test Mycelia._require_valid_autocycler_gfa(
                comma_identifier_gfa,
                "Autocycler comma-safe GFA",
            ) == comma_identifier_gfa

            malformed_gfa_cases = (
                (
                    name = "leading-garbage",
                    contents = "not GFA\nS\tcontig_1\tACGT\n",
                    message = "unknown GFA record type",
                ),
                (
                    name = "malformed-link",
                    contents =
                        "S\tleft\tACGT\nS\tright\tTGCA\n" *
                        "L\tleft\t+\tright\t-\n",
                    message = "malformed GFA link",
                ),
                (
                    name = "dangling-link",
                    contents =
                        "S\tleft\tACGT\n" *
                        "L\tleft\t+\tmissing\t-\t1M\n",
                    message = "dangling GFA link segment reference",
                ),
                (
                    name = "dangling-path",
                    contents =
                        "S\tleft\tACGT\n" *
                        "P\tprimary\tleft+,missing-\t1M\n",
                    message = "dangling GFA path segment reference",
                ),
                (
                    name = "shared-segment-path-name",
                    contents =
                        "S\tshared\tACGT\n" *
                        "P\tshared\tshared+\t*\n",
                    message = "duplicate GFA segment/path name",
                ),
                (
                    name = "invalid-segment-name",
                    contents = "S\t*invalid\tACGT\n",
                    message = "invalid GFA segment identifier",
                ),
                (
                    name = "invalid-tag-name",
                    contents = "S\tleft\tACGT\tX:i:4\n",
                    message = "invalid GFA segment optional tag name",
                ),
                (
                    name = "invalid-tag-type",
                    contents = "S\tleft\tACGT\tLN:q:4\n",
                    message =
                        "invalid or unsupported GFA segment optional tag value",
                ),
                (
                    name = "duplicate-tag-name",
                    contents = "S\tleft\tACGT\tLN:i:4\tLN:i:4\n",
                    message = "duplicate GFA segment optional tag name",
                ),
            )
            for malformed_case in malformed_gfa_cases
                malformed_path = joinpath(
                    temp_dir,
                    "$(malformed_case.name).gfa",
                )
                write(malformed_path, malformed_case.contents)
                malformed_error = _autocycler_test_error() do
                    Mycelia._require_valid_autocycler_gfa(
                        malformed_path,
                        "Autocycler malformed GFA",
                    )
                end
                Test.@test malformed_error isa ErrorException
                Test.@test occursin(
                    malformed_case.message,
                    sprint(showerror, malformed_error),
                )
            end

            malformed_gfa_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    gfa_path = only(filter(
                        path -> endswith(path, ".gfa"),
                        step.expected_outputs,
                    ))
                    write(gfa_path, "not GFA\n")
                end
                return nothing
            end
            malformed_gfa_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "malformed-gfa");
                    dependency_checker = () -> expected_toolchain,
                    runner = malformed_gfa_runner,
                )
            end
            Test.@test malformed_gfa_error isa ErrorException
            Test.@test occursin(
                "unknown GFA record type",
                sprint(showerror, malformed_gfa_error),
            )

            gapped_gfa_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    gfa_path = only(filter(
                        path -> endswith(path, ".gfa"),
                        step.expected_outputs,
                    ))
                    write(gfa_path, "S\tgapped\tAC-GT\n")
                end
                return nothing
            end
            gapped_gfa_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "gapped-gfa");
                    dependency_checker = () -> expected_toolchain,
                    runner = gapped_gfa_runner,
                )
            end
            Test.@test gapped_gfa_error isa ErrorException
            Test.@test occursin(
                "invalid DNA for GFA segment",
                sprint(showerror, gapped_gfa_error),
            )

            gfa_symlink_target = joinpath(temp_dir, "gfa-symlink-target")
            write(gfa_symlink_target, "S\tsymlinked\tACGT\n")
            symlinked_gfa_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    gfa_path = only(filter(
                        path -> endswith(path, ".gfa"),
                        step.expected_outputs,
                    ))
                    rm(gfa_path)
                    symlink(gfa_symlink_target, gfa_path)
                end
                return nothing
            end
            symlinked_gfa_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "symlinked-gfa");
                    dependency_checker = () -> expected_toolchain,
                    runner = symlinked_gfa_runner,
                )
            end
            Test.@test symlinked_gfa_error isa ErrorException
            Test.@test occursin(
                "symlinked artifact",
                sprint(showerror, symlinked_gfa_error),
            )
        end
    end

    Test.@testset "Paired-short polishing command plan" begin
        mktempdir() do temp_dir
            assembly = joinpath(temp_dir, "draft.fasta")
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            out_dir = joinpath(temp_dir, "results")
            write(assembly, ">contig\nACGT\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")

            plan = Mycelia._autocycler_polishing_command_plan(
                assembly,
                short_reads_1,
                short_reads_2,
                out_dir;
                threads = 6,
                polypolish_careful = true,
                conda_runner = "/test/conda",
            )
            Test.@test Tuple(step.name for step in plan.steps) == (
                :bwa_index,
                :bwa_mem_1,
                :bwa_mem_2,
                :polypolish_filter,
                :polypolish,
                :pypolca,
            )
            Test.@test plan.steps[1].expected_outputs == String[
                "$(abspath(assembly)).$(extension)" for
                extension in ("amb", "ann", "bwt", "pac", "sa")
            ]
            Test.@test plan.bwa_index_files == plan.steps[1].expected_outputs

            bwa_mem_1 = plan.steps[2].command.exec
            bwa_mem_2 = plan.steps[3].command.exec
            Test.@test "-a" in bwa_mem_1
            Test.@test "-a" in bwa_mem_2
            Test.@test bwa_mem_1[end] == abspath(short_reads_1)
            Test.@test bwa_mem_2[end] == abspath(short_reads_2)
            Test.@test !(abspath(short_reads_2) in bwa_mem_1)
            Test.@test !(abspath(short_reads_1) in bwa_mem_2)

            filter_arguments = plan.steps[4].command.exec
            Test.@test "--in1" in filter_arguments
            Test.@test "--in2" in filter_arguments
            Test.@test "--out1" in filter_arguments
            Test.@test "--out2" in filter_arguments
            Test.@test "--careful" in plan.steps[5].command.exec
            Test.@test "--careful" in plan.steps[6].command.exec
            Test.@test endswith(
                plan.assembly,
                joinpath("pypolca", "autocycler_polished_corrected.fasta"),
            )

            noncareful_plan = Mycelia._autocycler_polishing_command_plan(
                assembly,
                short_reads_1,
                short_reads_2,
                out_dir;
                polypolish_careful = false,
            )
            Test.@test !("--careful" in noncareful_plan.steps[5].command.exec)
            Test.@test "--careful" in noncareful_plan.steps[6].command.exec
            Test.@test !ispath(out_dir)

            requested_polishing_plan = Mycelia._autocycler_polishing_command_plan(
                assembly,
                short_reads_1,
                short_reads_2,
                out_dir;
                threads = 256,
                conda_runner = "/test/conda",
            )
            requested_bwa_arguments = requested_polishing_plan.steps[2].command.exec
            requested_bwa_threads = findfirst(==("-t"), requested_bwa_arguments)
            Test.@test requested_bwa_arguments[requested_bwa_threads + 1] == "256"
            requested_pypolca_arguments =
                requested_polishing_plan.steps[6].command.exec
            requested_pypolca_threads =
                findfirst(==("-t"), requested_pypolca_arguments)
            Test.@test requested_pypolca_arguments[requested_pypolca_threads + 1] ==
                       "256"
        end
    end

    Test.@testset "Full fake polished workflow preserves every stage" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")
            dependency_checks = Ref(0)
            polished_out_dir = joinpath(temp_dir, "polished-run")
            conda_runner = joinpath(
                temp_dir,
                "polished-conda",
                "bin",
                "conda",
            )
            mkpath(dirname(conda_runner))
            write(conda_runner, "synthetic runner\n")
            environment_prefix = joinpath(
                temp_dir,
                "polished-conda",
                "envs",
                Mycelia.AUTOCYCLER_ENV_NAME,
            )
            polished_lock_path =
                Mycelia._autocycler_output_lock_path(polished_out_dir)
            environment_lock_path = joinpath(temp_dir, "polished-environment.pid")
            environment_lock_active = Ref(false)
            environment_lock_runner = function (
                    action::Function,
                    observed_path::AbstractString,
            )
                Test.@test observed_path == environment_lock_path
                Test.@test !environment_lock_active[]
                environment_lock_active[] = true
                try
                    return action()
                finally
                    environment_lock_active[] = false
                end
            end

            result = Mycelia._run_autocycler_polished(
                long_reads,
                short_reads_1,
                short_reads_2,
                polished_out_dir;
                threads = 4,
                jobs = 2,
                read_type = "ont_r10",
                conda_runner,
                environment_prefix,
                dependency_checker = () -> begin
                    Test.@test isfile(polished_lock_path)
                    Test.@test environment_lock_active[]
                    dependency_checks[] += 1
                    return _autocycler_test_toolchain()
                end,
                runner = (step::NamedTuple) -> begin
                    Test.@test isfile(polished_lock_path)
                    Test.@test environment_lock_active[]
                    Test.@test first(step.command.exec) == realpath(conda_runner)
                    prefix_index = findfirst(==("-p"), step.command.exec)
                    Test.@test prefix_index !== nothing
                    Test.@test step.command.exec[prefix_index + 1] ==
                               environment_prefix
                    return _autocycler_test_runner!(step)
                end,
                environment_lock_path,
                environment_lock_runner,
            )
            Test.@test dependency_checks[] == 3
            Test.@test isfile(result.graph)
            Test.@test isfile(result.autocycler_assembly)
            Test.@test isfile(result.polypolish_assembly)
            Test.@test isfile(result.assembly)
            Test.@test isfile(result.pypolca_report)
            Test.@test result.assembly != result.autocycler_assembly
            Test.@test result.assembly != result.polypolish_assembly
            Test.@test result.toolchain == _autocycler_test_toolchain()
            Test.@test isempty(result.intermediates)
            Test.@test result.requested_threads == 4
            Test.@test result.autocycler_assembly_threads == 4
            Test.@test result.polishing_threads == 4
            Test.@test result.input_snapshots.long_reads.path ==
                       abspath(long_reads)
            Test.@test result.input_snapshots.short_reads_1.size_bytes ==
                       filesize(short_reads_1)
            Test.@test length(
                result.input_snapshots.short_reads_2.sha256,
            ) == 64
            Test.@test !ispath(polished_lock_path)
            Test.@test !any(isfile, String[
                "$(result.autocycler_assembly).$(extension)" for
                extension in ("amb", "ann", "bwt", "pac", "sa")
            ])

            external_polishing_dir = joinpath(temp_dir, "external-polishing")
            mkpath(external_polishing_dir)
            external_sentinel = joinpath(external_polishing_dir, "sentinel.txt")
            write(external_sentinel, "preserve\n")
            symlinked_polishing_out = joinpath(
                temp_dir,
                "symlinked-polishing-directory",
            )
            symlinked_polishing_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    symlink(
                        external_polishing_dir,
                        joinpath(
                            symlinked_polishing_out,
                            "short_read_polishing",
                        ),
                    )
                end
                return nothing
            end
            symlinked_polishing_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    symlinked_polishing_out;
                    dependency_checker = _autocycler_test_toolchain,
                    runner = symlinked_polishing_runner,
                )
            end
            Test.@test symlinked_polishing_error isa ErrorException
            Test.@test occursin(
                "symbolic-link component",
                sprint(showerror, symlinked_polishing_error),
            )
            Test.@test read(external_sentinel, String) == "preserve\n"
            Test.@test readdir(external_polishing_dir) == ["sentinel.txt"]

            initial_toolchain = _autocycler_test_toolchain()
            assembly_mutation_toolchain = _autocycler_test_toolchain(
                channels = Dict("pypolca" => "mutated-channel"),
            )
            assembly_mutation_checks = Ref(0)
            assembly_mutation_steps = Symbol[]
            assembly_mutation_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "assembly-toolchain-mutation");
                    dependency_checker = () -> begin
                        Test.@test environment_lock_active[]
                        assembly_mutation_checks[] += 1
                        return assembly_mutation_checks[] == 1 ?
                               initial_toolchain : assembly_mutation_toolchain
                    end,
                    runner = (step::NamedTuple) -> begin
                        Test.@test environment_lock_active[]
                        push!(assembly_mutation_steps, step.name)
                        return _autocycler_test_runner!(step)
                    end,
                    environment_lock_path,
                    environment_lock_runner,
                )
            end
            Test.@test assembly_mutation_error isa ErrorException
            Test.@test assembly_mutation_checks[] == 2
            Test.@test assembly_mutation_steps == [:autocycler]
            Test.@test occursin(
                "changed across long-read assembly within the polished workflow",
                sprint(showerror, assembly_mutation_error),
            )

            lifecycle_mutation_toolchain = _autocycler_test_toolchain(
                builds = Dict("polypolish" => "mutated-build_2"),
            )
            lifecycle_mutation_checks = Ref(0)
            lifecycle_mutation_steps = Symbol[]
            lifecycle_mutation_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "lifecycle-toolchain-mutation");
                    dependency_checker = () -> begin
                        Test.@test environment_lock_active[]
                        lifecycle_mutation_checks[] += 1
                        return lifecycle_mutation_checks[] <= 2 ?
                               initial_toolchain : lifecycle_mutation_toolchain
                    end,
                    runner = (step::NamedTuple) -> begin
                        Test.@test environment_lock_active[]
                        push!(lifecycle_mutation_steps, step.name)
                        return _autocycler_test_runner!(step)
                    end,
                    environment_lock_path,
                    environment_lock_runner,
                )
            end
            Test.@test lifecycle_mutation_error isa ErrorException
            Test.@test lifecycle_mutation_checks[] == 3
            Test.@test last(lifecycle_mutation_steps) == :pypolca
            Test.@test occursin(
                "complete Autocycler and paired-short polishing lifecycle",
                sprint(showerror, lifecycle_mutation_error),
            )

            renamed_polypolish_pypolca_calls = Ref(0)
            renamed_polypolish_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :polypolish
                    write(
                        only(step.expected_outputs),
                        ">renamed_by_polypolish\nACGT\n",
                    )
                elseif step.name == :pypolca
                    renamed_polypolish_pypolca_calls[] += 1
                end
                return nothing
            end
            renamed_polypolish_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "renamed-polypolish-contig");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = renamed_polypolish_runner,
                )
            end
            Test.@test renamed_polypolish_error isa ErrorException
            Test.@test renamed_polypolish_pypolca_calls[] == 0
            Test.@test occursin(
                "changed contig identifiers between Autocycler consensus and Polypolish",
                sprint(showerror, renamed_polypolish_error),
            )

            dropped_polypolish_pypolca_calls = Ref(0)
            dropped_polypolish_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :autocycler
                    consensus_fasta = only(filter(
                        path -> endswith(path, ".fasta"),
                        step.expected_outputs,
                    ))
                    write(
                        consensus_fasta,
                        ">authoritative_consensus\nACGT\n" *
                        ">dropped_contig\nTGCA\n",
                    )
                elseif step.name == :pypolca
                    dropped_polypolish_pypolca_calls[] += 1
                end
                return nothing
            end
            dropped_polypolish_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "dropped-polypolish-contig");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = dropped_polypolish_runner,
                )
            end
            Test.@test dropped_polypolish_error isa ErrorException
            Test.@test dropped_polypolish_pypolca_calls[] == 0
            Test.@test occursin(
                "dropped_contig",
                sprint(showerror, dropped_polypolish_error),
            )

            renamed_pypolca_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    final_fasta = only(filter(
                        path -> endswith(path, "_corrected.fasta"),
                        step.expected_outputs,
                    ))
                    write(final_fasta, ">renamed_by_pypolca\nACGT\n")
                end
                return nothing
            end
            renamed_pypolca_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "renamed-pypolca-contig");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = renamed_pypolca_runner,
                )
            end
            Test.@test renamed_pypolca_error isa ErrorException
            Test.@test occursin(
                "changed contig identifiers between Polypolish and Pypolca",
                sprint(showerror, renamed_pypolca_error),
            )

            mutated_graph_out_dir =
                joinpath(temp_dir, "pypolca-mutated-raw-graph")
            mutated_graph_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    write(
                        joinpath(
                            mutated_graph_out_dir,
                            "autocycler_out",
                            "consensus_assembly.gfa",
                        ),
                        "H\tVN:Z:1.0\nS\tcontig_1\tTGCATGCA\n",
                    )
                end
                return nothing
            end
            mutated_graph_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    mutated_graph_out_dir;
                    dependency_checker = _autocycler_test_toolchain,
                    runner = mutated_graph_runner,
                )
            end
            Test.@test mutated_graph_error isa ErrorException
            Test.@test occursin(
                "Autocycler consensus GFA changed after its validated Autocycler snapshot",
                sprint(showerror, mutated_graph_error),
            )

            malformed_raw_out_dir =
                joinpath(temp_dir, "pypolca-malformed-raw-assembly")
            malformed_raw_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    write(
                        joinpath(
                            malformed_raw_out_dir,
                            "autocycler_out",
                            "consensus_assembly.fasta",
                        ),
                        "not FASTA\n",
                    )
                end
                return nothing
            end
            malformed_raw_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    malformed_raw_out_dir;
                    dependency_checker = _autocycler_test_toolchain,
                    runner = malformed_raw_runner,
                )
            end
            Test.@test malformed_raw_error isa ErrorException
            Test.@test occursin(
                "Autocycler consensus FASTA is not valid FASTA",
                sprint(showerror, malformed_raw_error),
            )

            mutated_polypolish_out_dir =
                joinpath(temp_dir, "pypolca-mutated-polypolish")
            mutated_polypolish_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    write(
                        joinpath(
                            mutated_polypolish_out_dir,
                            "short_read_polishing",
                            "polypolish.fasta",
                        ),
                        ">authoritative_consensus\nCCCCCCCC\n",
                    )
                end
                return nothing
            end
            mutated_polypolish_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    mutated_polypolish_out_dir;
                    dependency_checker = _autocycler_test_toolchain,
                    runner = mutated_polypolish_runner,
                )
            end
            Test.@test mutated_polypolish_error isa ErrorException
            Test.@test occursin(
                "Polypolish Autocycler assembly changed after its validated Autocycler snapshot",
                sprint(showerror, mutated_polypolish_error),
            )

            symlinked_earlier_target =
                joinpath(temp_dir, "symlinked-earlier-target.fasta")
            write(
                symlinked_earlier_target,
                ">authoritative_consensus\nTTAACCGG\n",
            )
            symlinked_earlier_out_dir =
                joinpath(temp_dir, "pypolca-symlinked-polypolish")
            symlinked_earlier_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    polypolish_path = joinpath(
                        symlinked_earlier_out_dir,
                        "short_read_polishing",
                        "polypolish.fasta",
                    )
                    rm(polypolish_path)
                    symlink(symlinked_earlier_target, polypolish_path)
                end
                return nothing
            end
            symlinked_earlier_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    symlinked_earlier_out_dir;
                    dependency_checker = _autocycler_test_toolchain,
                    runner = symlinked_earlier_runner,
                )
            end
            Test.@test symlinked_earlier_error isa ErrorException
            Test.@test occursin(
                "must be a regular, non-symlink artifact",
                sprint(showerror, symlinked_earlier_error),
            )

            escaped_pypolca_out_dir =
                joinpath(temp_dir, "pypolca-parent-symlink-escape")
            escaped_pypolca_target =
                joinpath(temp_dir, "escaped-pypolca-artifacts")
            escaped_pypolca_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    pypolca_dir = dirname(first(step.expected_outputs))
                    cp(pypolca_dir, escaped_pypolca_target; force = true)
                    rm(pypolca_dir; recursive = true)
                    symlink(escaped_pypolca_target, pypolca_dir)
                end
                return nothing
            end
            escaped_pypolca_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    escaped_pypolca_out_dir;
                    dependency_checker = _autocycler_test_toolchain,
                    runner = escaped_pypolca_runner,
                )
            end
            Test.@test escaped_pypolca_error isa ErrorException
            Test.@test occursin(
                "resolves through a symbolic-link component",
                sprint(showerror, escaped_pypolca_error),
            )

            invalid_polypolish_cases = (
                (
                    slug = "malformed-polypolish-fasta",
                    contents = "not FASTA\n",
                    message = "not valid FASTA",
                ),
                (
                    slug = "empty-polypolish-fasta",
                    contents = "",
                    message = "did not create a nonempty artifact",
                ),
                (
                    slug = "duplicate-polypolish-identifier",
                    contents = ">duplicate\nACGT\n>duplicate\nTGCA\n",
                    message = "duplicate FASTA identifier",
                ),
            )
            for invalid_case in invalid_polypolish_cases
                pypolca_calls = Ref(0)
                invalid_polypolish_runner = function (step::NamedTuple)
                    _autocycler_test_runner!(step)
                    if step.name == :polypolish
                        write(only(step.expected_outputs), invalid_case.contents)
                    elseif step.name == :pypolca
                        pypolca_calls[] += 1
                    end
                    return nothing
                end
                invalid_polypolish_error = _autocycler_test_error() do
                    Mycelia._run_autocycler_polished(
                        long_reads,
                        short_reads_1,
                        short_reads_2,
                        joinpath(temp_dir, invalid_case.slug);
                        dependency_checker = _autocycler_test_toolchain,
                        runner = invalid_polypolish_runner,
                    )
                end
                Test.@test invalid_polypolish_error isa ErrorException
                Test.@test occursin(
                    invalid_case.message,
                    sprint(showerror, invalid_polypolish_error),
                )
                Test.@test pypolca_calls[] == 0
            end

            polypolish_symlink_target =
                joinpath(temp_dir, "polypolish-symlink-target")
            write(polypolish_symlink_target, ">symlinked\nACGT\n")
            symlinked_polypolish_pypolca_calls = Ref(0)
            symlinked_polypolish_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :polypolish
                    polypolish_fasta = only(step.expected_outputs)
                    rm(polypolish_fasta)
                    symlink(polypolish_symlink_target, polypolish_fasta)
                elseif step.name == :pypolca
                    symlinked_polypolish_pypolca_calls[] += 1
                end
                return nothing
            end
            symlinked_polypolish_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "symlinked-polypolish-fasta");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = symlinked_polypolish_runner,
                )
            end
            Test.@test symlinked_polypolish_error isa ErrorException
            Test.@test occursin(
                "symlinked artifact",
                sprint(showerror, symlinked_polypolish_error),
            )
            Test.@test symlinked_polypolish_pypolca_calls[] == 0

            retained = Mycelia._run_autocycler_polished(
                long_reads,
                short_reads_1,
                short_reads_2,
                joinpath(temp_dir, "retained-polishing-run");
                dependency_checker = _autocycler_test_toolchain,
                runner = _autocycler_test_runner!,
                keep_intermediates = true,
            )
            Test.@test length(retained.intermediates) == 9
            Test.@test all(isfile, retained.intermediates)

            cleanup_attempts = String[]
            cleanup_failure_path = Ref("")
            cleanup_remover = function (path::AbstractString)
                normalized_path = String(path)
                push!(cleanup_attempts, normalized_path)
                if endswith(normalized_path, "filtered_1.sam")
                    cleanup_failure_path[] = normalized_path
                    throw(ErrorException("synthetic cleanup failure"))
                end
                rm(normalized_path; force = true)
                return nothing
            end
            cleanup_result = Test.@test_logs (
                :warn,
                r"Autocycler could not remove a polishing intermediate",
            ) min_level=Logging.Warn match_mode=:any begin
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "cleanup-failure-run");
                    threads = 256,
                    dependency_checker = _autocycler_test_toolchain,
                    runner = _autocycler_test_runner!,
                    intermediate_remover = cleanup_remover,
                )
            end
            expected_cleanup_plan = Mycelia._autocycler_polishing_command_plan(
                cleanup_result.autocycler_assembly,
                short_reads_1,
                short_reads_2,
                cleanup_result.outdir;
                threads = 256,
            )
            Test.@test cleanup_attempts == expected_cleanup_plan.intermediate_files
            Test.@test cleanup_result.intermediates == [cleanup_failure_path[]]
            Test.@test isfile(cleanup_failure_path[])
            Test.@test all(
                path -> path == cleanup_failure_path[] || !ispath(path),
                expected_cleanup_plan.intermediate_files,
            )
            Test.@test isfile(cleanup_result.assembly)
            Test.@test cleanup_result.requested_threads == 256
            Test.@test cleanup_result.autocycler_assembly_threads == 128
            Test.@test cleanup_result.polishing_threads == 256

            cleanup_mutation_out_dir =
                joinpath(temp_dir, "cleanup-mutated-final-assembly")
            cleanup_mutated_final = Ref(false)
            cleanup_mutating_remover = function (path::AbstractString)
                if !cleanup_mutated_final[]
                    write(
                        joinpath(
                            cleanup_mutation_out_dir,
                            "short_read_polishing",
                            "pypolca",
                            "autocycler_polished_corrected.fasta",
                        ),
                        ">authoritative_consensus\nCCCCCCCC\n",
                    )
                    cleanup_mutated_final[] = true
                end
                rm(path; force = true)
                return nothing
            end
            cleanup_mutation_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    cleanup_mutation_out_dir;
                    dependency_checker = _autocycler_test_toolchain,
                    runner = _autocycler_test_runner!,
                    intermediate_remover = cleanup_mutating_remover,
                )
            end
            Test.@test cleanup_mutation_error isa ErrorException
            Test.@test cleanup_mutated_final[]
            Test.@test occursin(
                "Pypolca-polished Autocycler assembly changed after its validated Autocycler snapshot",
                sprint(showerror, cleanup_mutation_error),
            )

            malformed_final_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    final_fasta = only(filter(
                        path -> endswith(path, "_corrected.fasta"),
                        step.expected_outputs,
                    ))
                    write(final_fasta, "not FASTA\n")
                end
                return nothing
            end
            malformed_final_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "malformed-polished-fasta");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = malformed_final_runner,
                )
            end
            Test.@test malformed_final_error isa ErrorException
            Test.@test occursin(
                "not valid FASTA",
                sprint(showerror, malformed_final_error),
            )

            invalid_final_dna_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    final_fasta = only(filter(
                        path -> endswith(path, "_corrected.fasta"),
                        step.expected_outputs,
                    ))
                    write(final_fasta, ">invalid_dna\nACGTZ\n")
                end
                return nothing
            end
            invalid_final_dna_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "invalid-polished-dna");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = invalid_final_dna_runner,
                )
            end
            Test.@test invalid_final_dna_error isa ErrorException
            Test.@test occursin(
                "contains invalid DNA",
                sprint(showerror, invalid_final_dna_error),
            )

            duplicate_final_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    final_fasta = only(filter(
                        path -> endswith(path, "_corrected.fasta"),
                        step.expected_outputs,
                    ))
                    write(
                        final_fasta,
                        ">duplicate\nACGT\n>duplicate\nTGCA\n",
                    )
                end
                return nothing
            end
            duplicate_final_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "duplicate-polished-identifier");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = duplicate_final_runner,
                )
            end
            Test.@test duplicate_final_error isa ErrorException
            Test.@test occursin(
                "duplicate FASTA identifier",
                sprint(showerror, duplicate_final_error),
            )

            polished_symlink_target =
                joinpath(temp_dir, "polished-symlink-target")
            write(polished_symlink_target, ">symlinked\nACGT\n")
            symlinked_final_runner = function (step::NamedTuple)
                _autocycler_test_runner!(step)
                if step.name == :pypolca
                    final_fasta = only(filter(
                        path -> endswith(path, "_corrected.fasta"),
                        step.expected_outputs,
                    ))
                    rm(final_fasta)
                    symlink(polished_symlink_target, final_fasta)
                end
                return nothing
            end
            symlinked_final_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "symlinked-polished-fasta");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = symlinked_final_runner,
                )
            end
            Test.@test symlinked_final_error isa ErrorException
            Test.@test occursin(
                "symlinked artifact",
                sprint(showerror, symlinked_final_error),
            )
        end
    end

    Test.@testset "Failed polishing cleans only route-owned intermediates" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")

            for keep_intermediates in (false, true)
                out_dir = joinpath(temp_dir, "failure-$(keep_intermediates)")
                marker_error = ErrorException("synthetic pypolca failure")
                runner = function (step::NamedTuple)
                    _autocycler_test_runner!(step)
                    step.name == :pypolca && throw(marker_error)
                    return nothing
                end
                observed_error = _autocycler_test_error() do
                    Mycelia._run_autocycler_polished(
                        long_reads,
                        short_reads_1,
                        short_reads_2,
                        out_dir;
                        keep_intermediates,
                        dependency_checker = _autocycler_test_toolchain,
                        runner,
                    )
                end
                Test.@test observed_error === marker_error

                autocycler_assembly = joinpath(
                    out_dir,
                    "autocycler_out",
                    "consensus_assembly.fasta",
                )
                polishing_plan = Mycelia._autocycler_polishing_command_plan(
                    autocycler_assembly,
                    short_reads_1,
                    short_reads_2,
                    out_dir,
                )
                Test.@test all(
                    path -> isfile(path) == keep_intermediates,
                    polishing_plan.intermediate_files,
                )
                Test.@test isfile(autocycler_assembly)
                Test.@test isfile(joinpath(
                    out_dir,
                    "autocycler_out",
                    "consensus_assembly.gfa",
                ))
                Test.@test isfile(polishing_plan.polypolish_assembly)
                Test.@test isfile(polishing_plan.assembly)
                Test.@test isfile(polishing_plan.pypolca_report)
            end
        end
    end

    Test.@testset "Polished preflight rejects output hazards and input aliases" begin
        mktempdir() do temp_dir
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            write(short_reads_1, "not FASTQ R1\n")
            write(short_reads_2, "not FASTQ R2\n")

            dependency_checks = Ref(0)
            runner_calls = Ref(0)
            dependency_checker = function ()
                dependency_checks[] += 1
                return nothing
            end
            runner = function (step::NamedTuple)
                runner_calls[] += 1
                return nothing
            end

            blank_output_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    short_reads_1,
                    short_reads_1,
                    short_reads_2,
                    "   ";
                    dependency_checker,
                    runner,
                )
            end
            Test.@test blank_output_error isa ArgumentError
            Test.@test occursin(
                "non-blank path",
                sprint(showerror, blank_output_error),
            )

            blocking_ancestor = joinpath(temp_dir, "blocking-output-ancestor")
            write(blocking_ancestor, "not a directory\n")
            ancestor_output_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    short_reads_1,
                    short_reads_1,
                    short_reads_2,
                    joinpath(blocking_ancestor, "nested", "results");
                    dependency_checker,
                    runner,
                )
            end
            Test.@test ancestor_output_error isa ArgumentError
            Test.@test occursin(
                "non-directory existing ancestor",
                sprint(showerror, ancestor_output_error),
            )

            symlink_long_reads = joinpath(temp_dir, "long-symlink.fastq")
            hardlink_long_reads = joinpath(temp_dir, "long-hardlink.fastq")
            symlink(short_reads_2, symlink_long_reads)
            Base.Filesystem.hardlink(short_reads_1, hardlink_long_reads)
            aliased_long_reads = (
                exact = short_reads_1,
                symlink = symlink_long_reads,
                hardlink = hardlink_long_reads,
            )
            for (alias_type, long_reads) in pairs(aliased_long_reads)
                out_dir = joinpath(temp_dir, "$(alias_type)-alias-should-not-run")
                alias_error = _autocycler_test_error() do
                    Mycelia._run_autocycler_polished(
                        long_reads,
                        short_reads_1,
                        short_reads_2,
                        out_dir;
                        dependency_checker,
                        runner,
                    )
                end
                Test.@test alias_error isa ArgumentError
                Test.@test occursin(
                    "physically distinct",
                    sprint(showerror, alias_error),
                )
                Test.@test !ispath(out_dir)
            end
            Test.@test dependency_checks[] == 0
            Test.@test runner_calls[] == 0
        end
    end

    Test.@testset "Paired-short inputs fail before long assembly starts" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            runner_calls = Ref(0)
            runner = (step::NamedTuple) -> begin
                runner_calls[] += 1
                return nothing
            end

            missing_pair_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    joinpath(temp_dir, "missing-R2.fastq"),
                    joinpath(temp_dir, "should-not-run");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = runner,
                )
            end
            Test.@test missing_pair_error isa ArgumentError
            Test.@test occursin(
                "Paired short-read R2 FASTQ not found",
                sprint(showerror, missing_pair_error),
            )
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "should-not-run"))

            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            write(short_reads_2, "@different/2\nACGT\n+\nIIII\n")
            dependency_checks = Ref(0)
            mismatch_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "mismatch-should-not-run");
                    dependency_checker = () -> begin
                        dependency_checks[] += 1
                        return nothing
                    end,
                    runner = runner,
                )
            end
            Test.@test mismatch_error isa ArgumentError
            Test.@test occursin("out of sync", sprint(showerror, mismatch_error))
            Test.@test dependency_checks[] == 0
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "mismatch-should-not-run"))

            same_file_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_1,
                    joinpath(temp_dir, "same-file-should-not-run");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = runner,
                )
            end
            Test.@test same_file_error isa ArgumentError
            Test.@test occursin("must be distinct files", sprint(showerror, same_file_error))
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "same-file-should-not-run"))

            write(short_reads_1, "@pair/2\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/1\nACGT\n+\nIIII\n")
            reversed_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "reversed-should-not-run");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = runner,
                )
            end
            Test.@test reversed_error isa ArgumentError
            Test.@test occursin(
                "invalid explicit mate roles",
                sprint(showerror, reversed_error),
            )
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "reversed-should-not-run"))

            write(short_reads_1, "@pair 1:N:0:ACGT\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair 2:N:0:ACGT\nACGT\n+\nIIII\n")
            Test.@test Mycelia._validate_autocycler_paired_fastqs(
                short_reads_1,
                short_reads_2,
            ) == 1

            write(
                short_reads_1,
                "@pair/1 note 2:N:0:ACGT\nACGT\n+\nIIII\n",
            )
            write(
                short_reads_2,
                "@pair/2 note 1:N:0:ACGT\nACGT\n+\nIIII\n",
            )
            Test.@test Mycelia._validate_autocycler_paired_fastqs(
                short_reads_1,
                short_reads_2,
            ) == 1

            write(
                short_reads_1,
                "@pair/1 2:N:0:ACGT:trailing\nACGT\n+\nIIII\n",
            )
            write(
                short_reads_2,
                "@pair/2 1:N:0:ACGT:trailing\nACGT\n+\nIIII\n",
            )
            Test.@test Mycelia._validate_autocycler_paired_fastqs(
                short_reads_1,
                short_reads_2,
            ) == 1

            write(short_reads_1, "@pair 2:N:0:ACGT\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair 1:N:0:ACGT\nACGT\n+\nIIII\n")
            casava_reversed_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "casava-reversed-should-not-run");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = runner,
                )
            end
            Test.@test casava_reversed_error isa ArgumentError
            Test.@test occursin(
                "invalid explicit mate roles",
                sprint(showerror, casava_reversed_error),
            )
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(
                temp_dir,
                "casava-reversed-should-not-run",
            ))

            write(short_reads_1, "@pair/1 2:N:0:ACGT\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/2 2:N:0:ACGT\nACGT\n+\nIIII\n")
            conflicting_role_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "conflicting-role-should-not-run");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = runner,
                )
            end
            Test.@test conflicting_role_error isa ArgumentError
            Test.@test occursin(
                "conflicting explicit mate roles",
                sprint(showerror, conflicting_role_error),
            )
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(
                temp_dir,
                "conflicting-role-should-not-run",
            ))

            fasta_pair = joinpath(temp_dir, "not-fastq.fasta")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(fasta_pair, ">pair/2\nACGT\n")
            fasta_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    fasta_pair,
                    joinpath(temp_dir, "fasta-should-not-run");
                    dependency_checker = _autocycler_test_toolchain,
                    runner = runner,
                )
            end
            Test.@test fasta_error isa ArgumentError
            Test.@test occursin("must be FASTQ", sprint(showerror, fasta_error))
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "fasta-should-not-run"))

            fasta_long = joinpath(temp_dir, "not-long-fastq.fasta")
            write(fasta_long, ">long\nACGT\n")
            long_dependency_checks = Ref(0)
            long_fasta_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    fasta_long,
                    joinpath(temp_dir, "long-fasta-should-not-run");
                    dependency_checker = () -> begin
                        long_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = runner,
                )
            end
            Test.@test long_fasta_error isa ArgumentError
            Test.@test occursin("must be a FASTQ", sprint(showerror, long_fasta_error))
            Test.@test long_dependency_checks[] == 0
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "long-fasta-should-not-run"))

            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")
            nonempty_out_dir = joinpath(temp_dir, "nonempty-polished")
            mkpath(nonempty_out_dir)
            write(joinpath(nonempty_out_dir, "owned.txt"), "keep\n")
            nonempty_dependency_checks = Ref(0)
            nonempty_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    nonempty_out_dir;
                    dependency_checker = () -> begin
                        nonempty_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = runner,
                )
            end
            Test.@test nonempty_error isa ArgumentError
            Test.@test occursin(
                "output directory must be empty",
                sprint(showerror, nonempty_error),
            )
            Test.@test nonempty_dependency_checks[] == 0
            Test.@test runner_calls[] == 0
            Test.@test isfile(joinpath(nonempty_out_dir, "owned.txt"))

            stale_environment_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "stale-env-should-not-run");
                    dependency_checker = () -> throw(
                        ErrorException("missing required packages: pypolca"),
                    ),
                    runner = runner,
                )
            end
            Test.@test stale_environment_error isa ErrorException
            Test.@test occursin(
                "missing required packages",
                sprint(showerror, stale_environment_error),
            )
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "stale-env-should-not-run"))
        end
    end

    Test.@testset "public custom environment round trip" begin
        mktempdir() do temp_dir
            physical_root = joinpath(temp_dir, "physical-conda")
            alias_root = joinpath(temp_dir, "alias-conda")
            physical_prefix = joinpath(
                physical_root,
                "envs",
                Mycelia.AUTOCYCLER_ENV_NAME,
            )
            alias_prefix = joinpath(
                alias_root,
                "envs",
                Mycelia.AUTOCYCLER_ENV_NAME,
            )
            mkpath(joinpath(physical_prefix, "conda-meta"))
            conda_runner = _autocycler_public_test_conda_runner(
                joinpath(physical_root, "bin", "conda"),
                physical_prefix,
            )
            symlink(physical_root, alias_root)

            installed_script = Mycelia.install_autocycler(;
                conda_runner,
                environment_prefix = alias_prefix,
            )
            Test.@test Mycelia._autocycler_script_is_verified(
                installed_script,
            )

            long_reads = joinpath(temp_dir, "long.fastq")
            short_reads_1 = joinpath(temp_dir, "short_R1.fastq")
            short_reads_2 = joinpath(temp_dir, "short_R2.fastq")
            write(long_reads, "@long\nACGTACGT\n+\nIIIIIIII\n")
            write(short_reads_1, "@pair/1\nACGTACGT\n+\nIIIIIIII\n")
            write(short_reads_2, "@pair/2\nACGTACGT\n+\nIIIIIIII\n")

            long_only = Mycelia.run_autocycler(;
                long_reads,
                out_dir = joinpath(temp_dir, "long-only"),
                threads = 2,
                conda_runner,
                environment_prefix = alias_prefix,
            )
            Test.@test isfile(long_only.assembly)
            Test.@test isfile(long_only.graph)
            Test.@test long_only.toolchain == _autocycler_test_toolchain()

            polished = Mycelia.run_autocycler_polished(;
                long_reads,
                short_reads_1,
                short_reads_2,
                out_dir = joinpath(temp_dir, "polished"),
                threads = 2,
                conda_runner,
                environment_prefix = alias_prefix,
            )
            Test.@test isfile(polished.assembly)
            Test.@test isfile(polished.graph)
            Test.@test isfile(polished.autocycler_assembly)
            Test.@test polished.toolchain == _autocycler_test_toolchain()
        end
    end

    Test.@testset "Real-smoke gate and fixtures fail closed" begin
        Test.@test !_autocycler_real_smoke_enabled(Dict{String, String}())
        for broad_gate in ("MYCELIA_RUN_ALL", "MYCELIA_RUN_EXTERNAL")
            Test.@test !_autocycler_real_smoke_enabled(Dict(
                broad_gate => "true",
            ))
        end
        missing_broad_gate_error = _autocycler_test_error() do
            _autocycler_real_smoke_enabled(Dict(
                "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            ))
        end
        Test.@test missing_broad_gate_error isa ArgumentError
        Test.@test occursin(
            "also requires MYCELIA_RUN_EXTERNAL=true",
            sprint(showerror, missing_broad_gate_error),
        )

        enabled_without_fixtures = Dict(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
        )
        Test.@test _autocycler_real_smoke_enabled(enabled_without_fixtures)
        missing_long_reads_error = _autocycler_test_error() do
            _autocycler_real_smoke_inputs(enabled_without_fixtures)
        end
        Test.@test missing_long_reads_error isa ArgumentError
        Test.@test occursin(
            "requires MYCELIA_AUTOCYCLER_LONG_READS",
            sprint(showerror, missing_long_reads_error),
        )

        incomplete_pair_error = _autocycler_test_error() do
            _autocycler_real_smoke_inputs(Dict(
                "MYCELIA_AUTOCYCLER_LONG_READS" => "long.fastq",
                "MYCELIA_AUTOCYCLER_SHORT_READS_1" => "R1.fastq",
            ))
        end
        Test.@test incomplete_pair_error isa ArgumentError
        Test.@test occursin(
            "Set both MYCELIA_AUTOCYCLER_SHORT_READS_1",
            sprint(showerror, incomplete_pair_error),
        )

        long_only = _autocycler_real_smoke_inputs(Dict(
            "MYCELIA_AUTOCYCLER_LONG_READS" => "long.fastq",
        ))
        Test.@test long_only.long_reads == "long.fastq"
        Test.@test isempty(long_only.short_reads_1)
        Test.@test isempty(long_only.short_reads_2)
        Test.@test long_only.read_type == "ont_r10"
    end
end

if run_autocycler_smoke
    Test.@testset "Autocycler gated real smoke" begin
        smoke_inputs = something(autocycler_smoke_inputs)

        Mycelia.install_autocycler()
        mktempdir() do temp_dir
            if isempty(smoke_inputs.short_reads_1)
                result = Mycelia.run_autocycler(
                    long_reads = smoke_inputs.long_reads,
                    out_dir = joinpath(temp_dir, "autocycler"),
                    read_type = smoke_inputs.read_type,
                )
                Test.@test isfile(result.graph)
                Test.@test filesize(result.graph) > 0
                Test.@test isfile(result.assembly)
                Test.@test filesize(result.assembly) > 0
                Test.@test Mycelia._require_autocycler_toolchain_provenance(
                    result.toolchain,
                ) == result.toolchain
            else
                result = Mycelia.run_autocycler_polished(
                    long_reads = smoke_inputs.long_reads,
                    short_reads_1 = smoke_inputs.short_reads_1,
                    short_reads_2 = smoke_inputs.short_reads_2,
                    out_dir = joinpath(temp_dir, "autocycler-polished"),
                    read_type = smoke_inputs.read_type,
                )
                Test.@test isfile(result.graph)
                Test.@test filesize(result.graph) > 0
                Test.@test isfile(result.assembly)
                Test.@test filesize(result.assembly) > 0
                Test.@test Mycelia._require_autocycler_toolchain_provenance(
                    result.toolchain,
                ) == result.toolchain
            end
        end
    end
else
    @info (
        "Skipping Autocycler real smoke; set MYCELIA_RUN_EXTERNAL=true and " *
        "MYCELIA_RUN_AUTOCYCLER_SMOKE=true to opt in."
    )
end

Test.@testset "Autocycler input content lifecycle binding" begin
    mktempdir() do temp_dir
        toolchain = _autocycler_test_toolchain()
        function write_inputs(name::AbstractString)::NamedTuple
            input_dir = joinpath(temp_dir, String(name))
            mkpath(input_dir)
            long_reads = joinpath(input_dir, "long.fastq")
            short_reads_1 = joinpath(input_dir, "R1.fastq")
            short_reads_2 = joinpath(input_dir, "R2.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")
            return (; long_reads, short_reads_1, short_reads_2)
        end

        semantic_long_inputs = write_inputs("semantic-long")
        semantic_long_dependency_calls = Ref(0)
        semantic_long_runner_calls = Ref(0)
        semantic_long_error = _autocycler_test_error() do
            Mycelia._run_autocycler(
                semantic_long_inputs.long_reads,
                joinpath(temp_dir, "semantic-long-output");
                dependency_checker = () -> begin
                    semantic_long_dependency_calls[] += 1
                    return toolchain
                end,
                runner = step -> begin
                    semantic_long_runner_calls[] += 1
                    _autocycler_test_runner!(step)
                end,
                after_input_semantic_validation_hook = snapshots -> write(
                    snapshots.long_reads.path,
                    "@long\nTGCA\n+\nIIII\n",
                ),
                environment_lock_path = joinpath(
                    temp_dir,
                    "semantic-long-environment.pid",
                ),
            )
        end
        Test.@test semantic_long_error isa ErrorException
        Test.@test occursin(
            "Long-read FASTQ changed after its initial path/size/SHA-256 snapshot",
            sprint(showerror, semantic_long_error),
        )
        Test.@test semantic_long_dependency_calls[] == 0
        Test.@test semantic_long_runner_calls[] == 0

        for role in (:short_reads_1, :short_reads_2)
            semantic_pair_inputs = write_inputs("semantic-$(role)")
            semantic_pair_dependency_calls = Ref(0)
            semantic_pair_runner_calls = Ref(0)
            semantic_pair_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    semantic_pair_inputs.long_reads,
                    semantic_pair_inputs.short_reads_1,
                    semantic_pair_inputs.short_reads_2,
                    joinpath(temp_dir, "semantic-$(role)-output");
                    dependency_checker = () -> begin
                        semantic_pair_dependency_calls[] += 1
                        return toolchain
                    end,
                    runner = step -> begin
                        semantic_pair_runner_calls[] += 1
                        _autocycler_test_runner!(step)
                    end,
                    after_input_semantic_validation_hook = snapshots -> begin
                        selected_snapshot = getproperty(snapshots, role)
                        mate = role == :short_reads_1 ? 1 : 2
                        write(
                            selected_snapshot.path,
                            "@pair/$(mate)\nTGCA\n+\nIIII\n",
                        )
                    end,
                    environment_lock_path = joinpath(
                        temp_dir,
                        "semantic-$(role)-environment.pid",
                    ),
                )
            end
            label = role == :short_reads_1 ?
                    "Paired short-read R1 FASTQ" :
                    "Paired short-read R2 FASTQ"
            Test.@test semantic_pair_error isa ErrorException
            Test.@test occursin(
                "$(label) changed after its initial path/size/SHA-256 snapshot",
                sprint(showerror, semantic_pair_error),
            )
            Test.@test semantic_pair_dependency_calls[] == 0
            Test.@test semantic_pair_runner_calls[] == 0
        end

        before_inputs = write_inputs("before-consumption")
        before_dependency_calls = Ref(0)
        before_runner_calls = Ref(0)
        before_error = _autocycler_test_error() do
            Mycelia._run_autocycler(
                before_inputs.long_reads,
                joinpath(temp_dir, "before-consumption-output");
                dependency_checker = () -> begin
                    before_dependency_calls[] += 1
                    if before_dependency_calls[] == 1
                        write(
                            before_inputs.long_reads,
                            "@long\nTGCA\n+\nIIII\n",
                        )
                    end
                    return toolchain
                end,
                runner = step -> begin
                    before_runner_calls[] += 1
                    _autocycler_test_runner!(step)
                end,
                environment_lock_path = joinpath(
                    temp_dir,
                    "before-consumption-environment.pid",
                ),
            )
        end
        Test.@test before_error isa ErrorException
        Test.@test occursin(
            "Long-read FASTQ changed after its initial path/size/SHA-256 snapshot",
            sprint(showerror, before_error),
        )
        Test.@test before_runner_calls[] == 0

        during_inputs = write_inputs("during-consumption")
        during_runner_calls = Ref(0)
        during_error = _autocycler_test_error() do
            Mycelia._run_autocycler(
                during_inputs.long_reads,
                joinpath(temp_dir, "during-consumption-output");
                dependency_checker = () -> toolchain,
                runner = step -> begin
                    during_runner_calls[] += 1
                    _autocycler_test_runner!(step)
                    write(
                        during_inputs.long_reads,
                        "@long\nTGCA\n+\nIIII\n",
                    )
                end,
                environment_lock_path = joinpath(
                    temp_dir,
                    "during-consumption-environment.pid",
                ),
            )
        end
        Test.@test during_error isa ErrorException
        Test.@test occursin(
            "Long-read FASTQ changed after its initial path/size/SHA-256 snapshot",
            sprint(showerror, during_error),
        )
        Test.@test during_runner_calls[] == 1

        pre_polish_inputs = write_inputs("pre-polish-consumption")
        pre_polish_steps = Symbol[]
        pre_polish_error = _autocycler_test_error() do
            Mycelia._run_autocycler_polished(
                pre_polish_inputs.long_reads,
                pre_polish_inputs.short_reads_1,
                pre_polish_inputs.short_reads_2,
                joinpath(temp_dir, "pre-polish-consumption-output");
                dependency_checker = () -> toolchain,
                runner = step -> begin
                    push!(pre_polish_steps, step.name)
                    _autocycler_test_runner!(step)
                    if step.name == :autocycler
                        write(
                            pre_polish_inputs.short_reads_1,
                            "@pair/1\nTGCA\n+\nIIII\n",
                        )
                    end
                end,
                environment_lock_path = joinpath(
                    temp_dir,
                    "pre-polish-consumption-environment.pid",
                ),
            )
        end
        Test.@test pre_polish_error isa ErrorException
        Test.@test occursin(
            "Paired short-read R1 FASTQ changed after its initial " *
            "path/size/SHA-256 snapshot",
            sprint(showerror, pre_polish_error),
        )
        Test.@test pre_polish_steps == [:autocycler, :bwa_index]

        final_inputs = write_inputs("after-final-consumption")
        final_steps = Symbol[]
        final_error = _autocycler_test_error() do
            Mycelia._run_autocycler_polished(
                final_inputs.long_reads,
                final_inputs.short_reads_1,
                final_inputs.short_reads_2,
                joinpath(temp_dir, "after-final-consumption-output");
                dependency_checker = () -> toolchain,
                runner = step -> begin
                    push!(final_steps, step.name)
                    _autocycler_test_runner!(step)
                    if step.name == :pypolca
                        write(
                            final_inputs.short_reads_2,
                            "@pair/2\nTGCA\n+\nIIII\n",
                        )
                    end
                end,
                environment_lock_path = joinpath(
                    temp_dir,
                    "after-final-consumption-environment.pid",
                ),
            )
        end
        Test.@test final_error isa ErrorException
        Test.@test occursin(
            "Paired short-read R2 FASTQ changed after its initial " *
            "path/size/SHA-256 snapshot",
            sprint(showerror, final_error),
        )
        Test.@test final_steps == [
            :autocycler,
            :bwa_index,
            :bwa_mem_1,
            :bwa_mem_2,
            :polypolish_filter,
            :polypolish,
            :pypolca,
        ]
    end
end
