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
            if run_autocycler
                throw(ArgumentError(
                    "Autocycler-polished smoke requires all FASTQ variables; " *
                    "missing: $(join(missing_inputs, ", ")).",
                ))
            else
                @info (
                    "Multi-input hybrid smoke skipped: provide real FASTQ paths " *
                    "in all required environment variables."
                ) missing_environment_variables = missing_inputs
                Test.@test_skip false
            end
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
                        short_read_tech,
                        long_read_tech,
                        String[],
                    )
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
