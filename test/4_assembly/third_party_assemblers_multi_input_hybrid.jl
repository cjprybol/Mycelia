# Real smoke coverage for the Rhizomorph multi-input hybrid workflows.
#
# This file is excluded from default CI by the `third_party_assemblers` filename
# and repeats the external-test gate so direct includes remain safe. Supply a
# small, already paired real-data fixture with:
#
#   MYCELIA_RUN_EXTERNAL=true
#   MYCELIA_HYBRID_SHORT_R1=/path/to/reads_R1.fastq.gz
#   MYCELIA_HYBRID_SHORT_R2=/path/to/reads_R2.fastq.gz
#   MYCELIA_HYBRID_LONG_READS=/path/to/long_reads.fastq.gz
#
# Optional technology selectors are `MYCELIA_HYBRID_SHORT_TECH` (`illumina` or
# `ultima`) and `MYCELIA_HYBRID_LONG_TECH` (`nanopore`, `pacbio_clr`, or
# `pacbio_hifi`). The defaults are `illumina` and `nanopore`.
#
# Autocycler plus Polypolish/Pypolca is substantially more expensive and may
# provision its conda environment. Opt in separately only when that tooling is
# available. The Autocycler fixture must be a bacterial isolate for which most
# or all alternative long-read assemblies are expected to be complete; a
# nonempty result does not validate that biological prerequisite.
#
#   MYCELIA_RUN_AUTOCYCLER_POLISHED=true
#   MYCELIA_AUTOCYCLER_READ_TYPE=ont_r10
#
# Run from the Mycelia repository root:
#
#   julia --project=. -e \
#     'include("test/4_assembly/third_party_assemblers_multi_input_hybrid.jl")'

import Test
import Mycelia

const _HYBRID_INPUT_ENV_VARS = (
    "MYCELIA_HYBRID_SHORT_R1",
    "MYCELIA_HYBRID_SHORT_R2",
    "MYCELIA_HYBRID_LONG_READS",
)
const _HYBRID_AUTOCYCLER_READ_TYPES = (
    :ont_r9,
    :ont_r10,
    :pacbio_clr,
    :pacbio_hifi,
)

function _hybrid_env_enabled(name::AbstractString)::Bool
    return lowercase(strip(get(ENV, name, "false"))) == "true"
end

function _hybrid_env_symbol(
        name::AbstractString,
        default::Symbol,
        allowed::Tuple,
)::Symbol
    value = Symbol(lowercase(strip(get(ENV, name, String(default)))))
    if !(value in allowed)
        throw(ArgumentError(
            "$(name) must be one of $(allowed), got :$(value).",
        ))
    end
    return value
end

function _hybrid_required_env_symbol(
        name::AbstractString,
        allowed::Tuple,
)::Symbol
    text = strip(get(ENV, name, ""))
    isempty(text) && throw(ArgumentError(
        "$(name) is required when MYCELIA_RUN_AUTOCYCLER_POLISHED=true.",
    ))
    value = Symbol(lowercase(text))
    value in allowed || throw(ArgumentError(
        "$(name) must be one of $(allowed), got :$(value).",
    ))
    return value
end

function _hybrid_env_integer(
        name::AbstractString,
        default::Int,
        minimum::Int,
        maximum::Int,
)::Int
    text = strip(get(ENV, name, string(default)))
    value = tryparse(Int, text)
    value === nothing && throw(ArgumentError(
        "$(name) must be an integer, got $(repr(text)).",
    ))
    minimum <= value <= maximum || throw(ArgumentError(
        "$(name) must be between $(minimum) and $(maximum), got $(value).",
    ))
    return value
end

function _hybrid_input_paths()::NamedTuple
    paths = (
        short_r1 = abspath(ENV[_HYBRID_INPUT_ENV_VARS[1]]),
        short_r2 = abspath(ENV[_HYBRID_INPUT_ENV_VARS[2]]),
        long_reads = abspath(ENV[_HYBRID_INPUT_ENV_VARS[3]]),
    )
    for (label, path) in pairs(paths)
        if !isfile(path)
            throw(ArgumentError("Hybrid smoke $(label) FASTQ not found: $(path)"))
        end
        if filesize(path) == 0
            throw(ArgumentError("Hybrid smoke $(label) FASTQ is empty: $(path)"))
        end
    end
    return paths
end

function _assert_persistent_hybrid_result(
        result::Mycelia.Rhizomorph.AssemblyResult,
        workflow::AbstractString,
        output_dir::AbstractString,
        input_paths::NamedTuple,
        short_read_tech::Symbol,
        long_read_tech::Symbol,
        expected_polishers::Vector{String},
)::Nothing
    Test.@test result isa Mycelia.Rhizomorph.AssemblyResult
    Test.@test !isempty(result.contigs)
    Test.@test length(result.contigs) == length(result.contig_names)
    Test.@test result.gfa_compatible == false

    stats = result.assembly_stats
    Test.@test stats["method"] == "HybridAssembly"
    Test.@test stats["workflow"] == workflow
    expected_assembler = workflow == "autocycler_polished" ? "autocycler" : workflow
    Test.@test stats["assembler"] == expected_assembler
    Test.@test stats["input_contract"] == "paired_short_long"
    Test.@test stats["short_read_tech"] == String(short_read_tech)
    Test.@test stats["long_read_tech"] == String(long_read_tech)
    Test.@test stats["polishers"] == expected_polishers
    Test.@test stats["output_dir"] == output_dir
    Test.@test stats["num_contigs"] == length(result.contigs)
    Test.@test stats["input_technologies"] == Dict(
        "short_r1" => String(short_read_tech),
        "short_r2" => String(short_read_tech),
        "long_reads" => String(long_read_tech),
    )
    workflow_settings = stats["workflow_settings"]
    Test.@test workflow_settings["workflow"] == workflow
    Test.@test workflow_settings["threads"] > 0
    tool_artifacts = stats["tool_artifacts"]
    Test.@test tool_artifacts isa Dict{String, String}
    Test.@test isfile(tool_artifacts["final_assembly"])
    Test.@test all(isfile, values(tool_artifacts))

    input_counts = stats["input_read_counts"]
    corrected_counts = stats["corrected_read_counts"]
    for label in ("short_r1", "short_r2", "long_reads")
        Test.@test input_counts[label] > 0
        Test.@test corrected_counts[label] > 0
    end
    Test.@test input_counts["short_r1"] == input_counts["short_r2"]
    Test.@test corrected_counts["short_r1"] == corrected_counts["short_r2"]

    corrected_paths = stats["corrected_fastqs"]
    Test.@test corrected_paths isa Dict{String, String}
    expected_paths = Dict(
        "short_r1" => joinpath(
            output_dir,
            "corrected",
            "short_r1",
            "corrected.fastq",
        ),
        "short_r2" => joinpath(
            output_dir,
            "corrected",
            "short_r2",
            "corrected.fastq",
        ),
        "long_reads" => joinpath(
            output_dir,
            "corrected",
            "long_reads",
            "corrected.fastq",
        ),
    )
    Test.@test corrected_paths == expected_paths
    for path in values(corrected_paths)
        Test.@test isfile(path)
        Test.@test filesize(path) > 0
    end

    content_provenance = stats["read_content_provenance"]
    content_provenance isa AbstractDict || error(
        "Hybrid smoke result is missing read-content provenance.",
    )
    source_contract = content_provenance["source_inputs"]
    corrected_contract = content_provenance["corrected_fastqs"]
    source_contract isa AbstractDict || error(
        "Hybrid smoke result has malformed source-input provenance.",
    )
    corrected_contract isa AbstractDict || error(
        "Hybrid smoke result has malformed corrected-FASTQ provenance.",
    )
    Test.@test source_contract["schema"] ==
               "mycelia-paired-short-long-input-content-v1"
    Test.@test corrected_contract["schema"] ==
               "mycelia-corrected-fastq-content-v1"

    expected_input_paths = Dict(
        "short_r1" => input_paths.short_r1,
        "short_r2" => input_paths.short_r2,
        "long_reads" => input_paths.long_reads,
    )
    for label in ("short_r1", "short_r2", "long_reads")
        expected_input_path = abspath(expected_input_paths[label])
        source_set = source_contract[label]
        source_set isa AbstractDict || error(
            "Hybrid smoke $(label) source provenance is malformed.",
        )
        Test.@test source_set["kind"] == "path_set"
        Test.@test source_set["source_count"] == 1
        _assert_hybrid_sha256(source_set["sha256"], "$(label) source set")
        sources = source_set["sources"]
        sources isa AbstractVector || error(
            "Hybrid smoke $(label) source list is malformed.",
        )
        Test.@test length(sources) == 1
        source = only(sources)
        source isa AbstractDict || error(
            "Hybrid smoke $(label) source identity is malformed.",
        )
        Test.@test source["kind"] == "path"
        Test.@test source["source_index"] == 1
        Test.@test source["path"] == expected_input_path
        Test.@test source["canonical_path"] == realpath(expected_input_path)
        Test.@test source["size_bytes"] == filesize(expected_input_path)
        Test.@test _assert_hybrid_sha256(
            source["sha256"],
            "$(label) source",
        ) == Mycelia.Rhizomorph._multi_input_file_sha256(expected_input_path)

        corrected_path = corrected_paths[label]
        corrected_identity = corrected_contract[label]
        corrected_identity isa AbstractDict || error(
            "Hybrid smoke $(label) corrected identity is malformed.",
        )
        Test.@test corrected_identity["kind"] == "path"
        Test.@test corrected_identity["source_index"] == 1
        Test.@test corrected_identity["path"] == corrected_path
        Test.@test corrected_identity["canonical_path"] == realpath(corrected_path)
        Test.@test corrected_identity["size_bytes"] == filesize(corrected_path)
        Test.@test _assert_hybrid_sha256(
            corrected_identity["sha256"],
            "$(label) corrected FASTQ",
        ) == Mycelia.Rhizomorph._multi_input_file_sha256(corrected_path)
    end
    return nothing
end

function _assert_hybrid_sha256(value::Any, label::AbstractString)::String
    value isa AbstractString || error(
        "Hybrid smoke $(label) SHA-256 is not a string.",
    )
    digest = String(value)
    Test.@test match(r"^[0-9a-f]{64}$", digest) !== nothing
    return digest
end

function _assert_unicycler_toolchain(stats::AbstractDict)::Nothing
    toolchain = stats["toolchain"]
    toolchain isa AbstractDict || error(
        "Unicycler smoke result is missing toolchain provenance.",
    )
    Test.@test toolchain["environment_name"] == "unicycler"
    Test.@test toolchain["inventory_schema"] ==
               "conda-name-version-build-channel-v1"
    inventory_sha256 = _assert_hybrid_sha256(
        toolchain["package_inventory_sha256"],
        "Unicycler package inventory",
    )
    packages = toolchain["packages"]
    packages isa AbstractVector || error(
        "Unicycler smoke result has no realized package inventory.",
    )
    Test.@test !isempty(packages)
    package_names = Set{String}()
    for (package_index, package) in enumerate(packages)
        package isa AbstractDict || error(
            "Unicycler package record $(package_index) is malformed.",
        )
        Test.@test Set(String.(keys(package))) ==
                   Set(["name", "version", "build", "channel"])
        for field in ("name", "version", "build", "channel")
            value = package[field]
            value isa AbstractString || error(
                "Unicycler package $(package_index) $(field) is not a string.",
            )
            Test.@test !isempty(String(value))
        end
        push!(package_names, String(package["name"]))
    end
    Test.@test "unicycler" in package_names
    Test.@test "spades" in package_names
    normalized_packages =
        Mycelia.Rhizomorph._normalize_unicycler_package_inventory(packages)
    expected_inventory_sha256 =
        Mycelia.Rhizomorph._unicycler_package_inventory_sha256(
            normalized_packages,
        )
    Test.@test inventory_sha256 == expected_inventory_sha256
    return nothing
end

function _assert_autocycler_toolchain(
        stats::AbstractDict,
)::Nothing
    toolchain = stats["toolchain"]
    toolchain isa AbstractDict || error(
        "Autocycler smoke result is missing toolchain provenance.",
    )
    Test.@test toolchain["autocycler_script_revision"] ==
               Mycelia.AUTOCYCLER_SCRIPT_REVISION
    Test.@test toolchain["autocycler_script_sha256"] ==
               Mycelia.AUTOCYCLER_SCRIPT_SHA256
    Test.@test match(
        r"^[0-9a-f]{64}$",
        String(toolchain["autocycler_script_sha256"]),
    ) !== nothing

    _, _, environment_path = Mycelia._autocycler_paths()
    expected_environment_sha256 = Mycelia._autocycler_sha256(environment_path)
    Test.@test toolchain["environment_name"] == Mycelia.AUTOCYCLER_ENV_NAME
    Test.@test toolchain["environment_spec_expected_sha256"] ==
               Mycelia.AUTOCYCLER_ENVIRONMENT_SPEC_SHA256
    Test.@test toolchain["environment_spec_sha256"] ==
               expected_environment_sha256
    Test.@test expected_environment_sha256 ==
               Mycelia.AUTOCYCLER_ENVIRONMENT_SPEC_SHA256
    Test.@test match(
        r"^[0-9a-f]{64}$",
        String(toolchain["environment_spec_sha256"]),
    ) !== nothing

    packages = toolchain["packages"]
    packages isa AbstractDict || error(
        "Autocycler smoke result is missing installed package provenance.",
    )
    installed_versions = Dict{String, String}()
    for specification in Mycelia.AUTOCYCLER_REQUIRED_PACKAGE_SPECS
        Test.@test haskey(packages, specification.name)
        haskey(packages, specification.name) || continue
        version = String(packages[specification.name])
        Test.@test !isempty(version)
        installed_versions[specification.name] = version
    end
    Test.@test isempty(Mycelia._autocycler_package_issues(installed_versions))
    return nothing
end

Test.@testset "multi-input hybrid external smoke" begin
    run_external = _hybrid_env_enabled("MYCELIA_RUN_ALL") ||
                   _hybrid_env_enabled("MYCELIA_RUN_EXTERNAL")
    run_autocycler = _hybrid_env_enabled("MYCELIA_RUN_AUTOCYCLER_POLISHED")

    if run_autocycler && !run_external
        throw(ArgumentError(
            "MYCELIA_RUN_AUTOCYCLER_POLISHED=true also requires " *
            "MYCELIA_RUN_EXTERNAL=true (or MYCELIA_RUN_ALL=true).",
        ))
    end

    if !run_external
        @info (
            "Multi-input hybrid smoke skipped: set " *
            "MYCELIA_RUN_EXTERNAL=true to enable external tools."
        )
        Test.@test_skip false
    else
        missing_inputs = [
            name for name in _HYBRID_INPUT_ENV_VARS if isempty(strip(get(ENV, name, "")))
        ]
        if !isempty(missing_inputs)
            throw(ArgumentError(
                "Enabled multi-input hybrid smoke requires all FASTQ " *
                "variables; missing: $(join(missing_inputs, ", ")).",
            ))
        else
            inputs = _hybrid_input_paths()
            short_read_tech = _hybrid_env_symbol(
                "MYCELIA_HYBRID_SHORT_TECH",
                :illumina,
                (:illumina, :ultima),
            )
            long_read_tech = _hybrid_env_symbol(
                "MYCELIA_HYBRID_LONG_TECH",
                :nanopore,
                (:nanopore, :pacbio_clr, :pacbio_hifi),
            )
            threads = clamp(
                something(
                    tryparse(
                        Int,
                        get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2"),
                    ),
                    2,
                ),
                1,
                4,
            )
            autocycler_read_type = nothing
            autocycler_jobs = nothing
            if run_autocycler
                autocycler_read_type = _hybrid_required_env_symbol(
                    "MYCELIA_AUTOCYCLER_READ_TYPE",
                    _HYBRID_AUTOCYCLER_READ_TYPES,
                )
                autocycler_jobs = _hybrid_env_integer(
                    "MYCELIA_AUTOCYCLER_TEST_JOBS",
                    1,
                    1,
                    4,
                )
            end

            mktempdir() do temporary_root
                Test.@testset "corrected paired-short plus long to Unicycler" begin
                    output_dir = joinpath(temporary_root, "unicycler")
                    config = Mycelia.Rhizomorph.UnicyclerHybridConfig(;
                        short_read_tech,
                        long_read_tech,
                        output_dir,
                        threads,
                    )
                    result = Mycelia.Rhizomorph.assemble_hybrid(
                        (inputs.short_r1, inputs.short_r2),
                        inputs.long_reads;
                        config,
                    )
                    _assert_persistent_hybrid_result(
                        result,
                        "unicycler",
                        output_dir,
                        inputs,
                        short_read_tech,
                        long_read_tech,
                        String[],
                    )
                    _assert_unicycler_toolchain(result.assembly_stats)
                end

                if !run_autocycler
                    @info (
                        "Autocycler-polished smoke skipped: set " *
                        "MYCELIA_RUN_AUTOCYCLER_POLISHED=true only when the " *
                        "Autocycler/Polypolish/Pypolca environment is available."
                    )
                    Test.@test_skip false
                else
                    Test.@testset (
                        "Autocycler consensus plus paired-short polishing"
                    ) begin
                        output_dir = joinpath(
                            temporary_root,
                            "autocycler_polished",
                        )
                        config = Mycelia.Rhizomorph.AutocyclerPolishConfig(;
                            short_read_tech,
                            long_read_tech,
                            autocycler_read_type = something(autocycler_read_type),
                            output_dir,
                            threads,
                            jobs = something(autocycler_jobs),
                        )
                        result = Mycelia.Rhizomorph.assemble_hybrid(
                            (inputs.short_r1, inputs.short_r2),
                            inputs.long_reads;
                            config,
                        )
                        resolved_long_read_tech =
                            Mycelia.Rhizomorph._long_read_correction_technology(
                                config,
                            )
                        _assert_persistent_hybrid_result(
                            result,
                            "autocycler_polished",
                            output_dir,
                            inputs,
                            short_read_tech,
                            resolved_long_read_tech,
                            ["polypolish-careful", "pypolca-careful"],
                        )
                        _assert_autocycler_toolchain(result.assembly_stats)
                    end
                end
            end
        end
    end
end
