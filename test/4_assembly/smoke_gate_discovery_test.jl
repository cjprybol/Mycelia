# Regression coverage for canonical smoke-preflight wiring and its shared parser.

import Test

if !isdefined(@__MODULE__, :_multi_input_hybrid_smoke_prerequisites)
    include(
        Base.joinpath(
            @__DIR__,
            "..",
            "multi_input_hybrid_smoke_support.jl",
        ),
    )
end

const _SMOKE_DISCOVERY_ENV_NAMES = (
    "MYCELIA_RUN_ALL",
    "MYCELIA_RUN_EXTERNAL",
    "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE",
    "MYCELIA_RUN_AUTOCYCLER_POLISHED",
    "MYCELIA_RUN_AUTOCYCLER_SMOKE",
    "MYCELIA_HYBRID_SHORT_R1",
    "MYCELIA_HYBRID_SHORT_R2",
    "MYCELIA_HYBRID_LONG_READS",
    "MYCELIA_HYBRID_SHORT_TECH",
    "MYCELIA_HYBRID_LONG_TECH",
    "MYCELIA_AUTOCYCLER_READ_TYPE",
    "MYCELIA_AUTOCYCLER_TEST_JOBS",
    "MYCELIA_AUTOCYCLER_LONG_READS",
    "MYCELIA_AUTOCYCLER_SHORT_READS_1",
    "MYCELIA_AUTOCYCLER_SHORT_READS_2",
    "MYCELIA_ASSEMBLER_TEST_THREADS",
)

function _run_smoke_discovery(
        overrides::Pair{String, String}...,
)::NamedTuple
    project_root = Base.normpath(Base.joinpath(@__DIR__, "..", ".."))
    preflight_path = Base.joinpath(
        project_root,
        "test",
        "multi_input_hybrid_smoke_preflight.jl",
    )
    runner = `$(Base.julia_cmd()) --startup-file=no
              --project=$(project_root) $(preflight_path)`
    environment = Dict{String, String}(
        String(name) => String(value) for (name, value) in ENV
    )
    for name in _SMOKE_DISCOVERY_ENV_NAMES
        delete!(environment, name)
    end
    environment["LD_LIBRARY_PATH"] = ""
    for override in overrides
        environment[first(override)] = last(override)
    end

    output = IOBuffer()
    process = Base.run(Base.pipeline(
        Base.ignorestatus(Base.setenv(runner, environment));
        stdout = output,
        stderr = output,
    ))
    return (; ok = Base.success(process), output = String(Base.take!(output)))
end

function _smoke_prerequisite_error(action::Function)::Exception
    try
        action()
    catch error
        return error
    end
    Base.error("Expected smoke prerequisite validation to fail")
end

function _hybrid_smoke_environment(
        input_paths::NamedTuple,
        overrides::Pair{String, String}...,
)::Dict{String, String}
    environment = Dict{String, String}(
        "MYCELIA_RUN_EXTERNAL" => "true",
        "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
        "MYCELIA_HYBRID_SHORT_R1" => input_paths.short_r1,
        "MYCELIA_HYBRID_SHORT_R2" => input_paths.short_r2,
        "MYCELIA_HYBRID_LONG_READS" => input_paths.long_reads,
    )
    for override in overrides
        environment[first(override)] = last(override)
    end
    return environment
end

Test.@testset "canonical smoke discovery preflight" begin
    # Keep subprocess coverage only for canonical-driver wiring. Parser and
    # compatibility matrices below call the shared prerequisites directly.
    broad_only = _run_smoke_discovery("MYCELIA_RUN_EXTERNAL" => " true ")
    Test.@test broad_only.ok

    dedicated_only = _run_smoke_discovery(
        "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
    )
    Test.@test !dedicated_only.ok
    Test.@test occursin(
        "also require MYCELIA_RUN_EXTERNAL=true",
        dedicated_only.output,
    )

    standalone_only = _run_smoke_discovery(
        "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
    )
    Test.@test !standalone_only.ok
    Test.@test occursin(
        "MYCELIA_RUN_AUTOCYCLER_SMOKE=true also requires " *
        "MYCELIA_RUN_EXTERNAL=true",
        standalone_only.output,
    )

    Base.mktempdir() do temporary_root
        input_paths = (
            short_r1 = Base.joinpath(temporary_root, "reads_R1.fastq"),
            short_r2 = Base.joinpath(temporary_root, "reads_R2.fastq"),
            long_reads = Base.joinpath(temporary_root, "long.fastq"),
        )
        Base.write(input_paths.short_r1, "@pair/1\nACGT\n+\nIIII\n")
        Base.write(input_paths.short_r2, "@pair/2\nACGT\n+\nIIII\n")
        Base.write(input_paths.long_reads, "@long\nACGT\n+\nIIII\n")

        combined_wiring = _run_smoke_discovery(
            "MYCELIA_RUN_ALL" => " true ",
            "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
            "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_HYBRID_SHORT_R1" => input_paths.short_r1,
            "MYCELIA_HYBRID_SHORT_R2" => input_paths.short_r2,
            "MYCELIA_HYBRID_LONG_READS" => input_paths.long_reads,
            "MYCELIA_AUTOCYCLER_READ_TYPE" => "ont_r10",
            "MYCELIA_AUTOCYCLER_LONG_READS" => input_paths.long_reads,
            "MYCELIA_AUTOCYCLER_SHORT_READS_1" => input_paths.short_r1,
            "MYCELIA_AUTOCYCLER_SHORT_READS_2" => input_paths.short_r2,
            "MYCELIA_ASSEMBLER_TEST_THREADS" => " 4 ",
        )
        Test.@test combined_wiring.ok
        Test.@test isempty(strip(combined_wiring.output))
    end
end

Test.@testset "direct smoke prerequisite parser coverage" begin
    Test.@test _multi_input_hybrid_smoke_prerequisites(
        Dict{String, String}(),
    ) == (; run_smoke = false, run_autocycler = false)
    Test.@test _autocycler_smoke_prerequisites(
        Dict{String, String}(),
    ) == (; run_smoke = false)

    missing_parent = _smoke_prerequisite_error() do
        _multi_input_hybrid_smoke_prerequisites(Dict(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
        ))
    end
    Test.@test missing_parent isa ArgumentError
    Test.@test occursin(
        "also requires MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE=true",
        sprint(showerror, missing_parent),
    )

    missing_fixtures = _smoke_prerequisite_error() do
        _multi_input_hybrid_smoke_prerequisites(Dict(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
        ))
    end
    Test.@test missing_fixtures isa ArgumentError
    Test.@test occursin("missing:", sprint(showerror, missing_fixtures))

    Base.mktempdir() do temporary_root
        input_paths = (
            short_r1 = Base.joinpath(temporary_root, "direct_R1.fastq"),
            short_r2 = Base.joinpath(temporary_root, "direct_R2.fastq"),
            long_reads = Base.joinpath(temporary_root, "direct_long.fastq"),
        )
        Base.write(input_paths.short_r1, "@pair/1\nACGT\n+\nIIII\n")
        Base.write(input_paths.short_r2, "@pair/2\nACGT\n+\nIIII\n")
        Base.write(input_paths.long_reads, "@long\nACGT\n+\nIIII\n")
        environment = _hybrid_smoke_environment(input_paths)
        prerequisites = _multi_input_hybrid_smoke_prerequisites(environment)
        Test.@test prerequisites.run_smoke
        Test.@test !prerequisites.run_autocycler
        Test.@test prerequisites.threads == 2
        Test.@test prerequisites.input_paths == (
            short_r1 = Base.abspath(input_paths.short_r1),
            short_r2 = Base.abspath(input_paths.short_r2),
            long_reads = Base.abspath(input_paths.long_reads),
        )

        threaded_prerequisites = _multi_input_hybrid_smoke_prerequisites(
            _hybrid_smoke_environment(
                input_paths,
                "MYCELIA_ASSEMBLER_TEST_THREADS" => " 4 ",
            ),
        )
        Test.@test threaded_prerequisites.threads == 4

        for (configured_threads, expected_message) in (
                "many" => "must be an integer",
                "0" => "must be between 1 and 4",
                "5" => "must be between 1 and 4",
        )
            hybrid_threads_error = _smoke_prerequisite_error() do
                _multi_input_hybrid_smoke_prerequisites(
                    _hybrid_smoke_environment(
                        input_paths,
                        "MYCELIA_ASSEMBLER_TEST_THREADS" => configured_threads,
                    ),
                )
            end
            Test.@test hybrid_threads_error isa ArgumentError
            Test.@test occursin(
                expected_message,
                sprint(showerror, hybrid_threads_error),
            )

            autocycler_threads_error = _smoke_prerequisite_error() do
                _autocycler_smoke_prerequisites(Dict(
                    "MYCELIA_RUN_EXTERNAL" => "true",
                    "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
                    "MYCELIA_AUTOCYCLER_LONG_READS" => input_paths.long_reads,
                    "MYCELIA_ASSEMBLER_TEST_THREADS" => configured_threads,
                ))
            end
            Test.@test autocycler_threads_error isa ArgumentError
            Test.@test occursin(
                expected_message,
                sprint(showerror, autocycler_threads_error),
            )
        end

        malformed = Base.joinpath(temporary_root, "malformed.fastq")
        Base.write(malformed, ">not_fastq\nACGT\n")
        malformed_error = _smoke_prerequisite_error() do
            _multi_input_hybrid_smoke_prerequisites(
                _hybrid_smoke_environment(
                    merge(input_paths, (; long_reads = malformed)),
                ),
            )
        end
        Test.@test malformed_error isa ArgumentError
        Test.@test occursin(
            "must be a valid FASTQ file",
            sprint(showerror, malformed_error),
        )

        alias_r2 = Base.joinpath(temporary_root, "alias_R2.fastq")
        Base.symlink(input_paths.short_r1, alias_r2)
        alias_error = _smoke_prerequisite_error() do
            _multi_input_hybrid_smoke_prerequisites(
                _hybrid_smoke_environment(
                    merge(input_paths, (; short_r2 = alias_r2)),
                ),
            )
        end
        Test.@test alias_error isa ArgumentError
        Test.@test occursin(
            "must be physically distinct files",
            sprint(showerror, alias_error),
        )

        reversed_r1 = Base.joinpath(temporary_root, "reversed_R1.fastq")
        reversed_r2 = Base.joinpath(temporary_root, "reversed_R2.fastq")
        Base.write(reversed_r1, "@pair/2\nACGT\n+\nIIII\n")
        Base.write(reversed_r2, "@pair/1\nACGT\n+\nIIII\n")
        reversed_error = _smoke_prerequisite_error() do
            _multi_input_hybrid_smoke_prerequisites(
                _hybrid_smoke_environment(merge(
                    input_paths,
                    (; short_r1 = reversed_r1, short_r2 = reversed_r2),
                )),
            )
        end
        Test.@test reversed_error isa ArgumentError
        Test.@test occursin(
            "invalid explicit mate roles",
            sprint(showerror, reversed_error),
        )

        unsynced_r2 = Base.joinpath(temporary_root, "unsynced_R2.fastq")
        Base.write(unsynced_r2, "@different/2\nACGT\n+\nIIII\n")
        unsynced_error = _smoke_prerequisite_error() do
            _multi_input_hybrid_smoke_prerequisites(
                _hybrid_smoke_environment(
                    merge(input_paths, (; short_r2 = unsynced_r2)),
                ),
            )
        end
        Test.@test unsynced_error isa ArgumentError
        Test.@test occursin(
            "out of sync at record 1",
            sprint(showerror, unsynced_error),
        )

        incompatible_read_types = (
            :nanopore => (:pacbio_clr, :pacbio_hifi),
            :pacbio_clr => (:ont_r9, :ont_r10, :pacbio_hifi),
            :pacbio_hifi => (:ont_r9, :ont_r10, :pacbio_clr),
        )
        for (long_read_tech, read_types) in incompatible_read_types
            for read_type in read_types
                incompatible_error = _smoke_prerequisite_error() do
                    _multi_input_hybrid_smoke_prerequisites(
                        _hybrid_smoke_environment(
                            input_paths,
                            "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
                            "MYCELIA_HYBRID_LONG_TECH" =>
                                String(long_read_tech),
                            "MYCELIA_AUTOCYCLER_READ_TYPE" => String(read_type),
                        ),
                    )
                end
                Test.@test incompatible_error isa ArgumentError
                Test.@test occursin(
                    "READ_TYPE=:$(read_type) is incompatible with " *
                    "MYCELIA_HYBRID_LONG_TECH=:$(long_read_tech)",
                    sprint(showerror, incompatible_error),
                )
            end
        end

        compatible_read_types = (
            :nanopore => (:ont_r9, :ont_r10),
            :pacbio_clr => (:pacbio_clr,),
            :pacbio_hifi => (:pacbio_hifi,),
        )
        for (long_read_tech, read_types) in compatible_read_types
            for read_type in read_types
                compatible = _multi_input_hybrid_smoke_prerequisites(
                    _hybrid_smoke_environment(
                        input_paths,
                        "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
                        "MYCELIA_HYBRID_LONG_TECH" => String(long_read_tech),
                        "MYCELIA_AUTOCYCLER_READ_TYPE" => String(read_type),
                    ),
                )
                Test.@test compatible.run_autocycler
                Test.@test compatible.autocycler_read_type == read_type
            end
        end

        half_pair_error = _smoke_prerequisite_error() do
            _autocycler_smoke_prerequisites(Dict(
                "MYCELIA_RUN_EXTERNAL" => "true",
                "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
                "MYCELIA_AUTOCYCLER_LONG_READS" => input_paths.long_reads,
                "MYCELIA_AUTOCYCLER_SHORT_READS_1" => input_paths.short_r1,
                "MYCELIA_AUTOCYCLER_READ_TYPE" => "ont_r10",
            ))
        end
        Test.@test half_pair_error isa ArgumentError
        Test.@test occursin(
            "Set both MYCELIA_AUTOCYCLER_SHORT_READS_1",
            sprint(showerror, half_pair_error),
        )

        for configured_read_type in (nothing, "   ")
            environment = Dict(
                "MYCELIA_RUN_EXTERNAL" => "true",
                "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
                "MYCELIA_AUTOCYCLER_LONG_READS" => input_paths.long_reads,
            )
            configured_read_type === nothing ||
                (environment["MYCELIA_AUTOCYCLER_READ_TYPE"] =
                    configured_read_type)
            missing_read_type = _smoke_prerequisite_error() do
                _autocycler_smoke_prerequisites(environment)
            end
            Test.@test missing_read_type isa ArgumentError
            Test.@test occursin(
                "MYCELIA_AUTOCYCLER_READ_TYPE is required",
                sprint(showerror, missing_read_type),
            )
        end

        invalid_read_type = _smoke_prerequisite_error() do
            _autocycler_smoke_prerequisites(Dict(
                "MYCELIA_RUN_ALL" => "true",
                "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
                "MYCELIA_AUTOCYCLER_LONG_READS" => input_paths.long_reads,
                "MYCELIA_AUTOCYCLER_READ_TYPE" => "illumina",
            ))
        end
        Test.@test invalid_read_type isa ArgumentError
        Test.@test occursin(
            "MYCELIA_AUTOCYCLER_READ_TYPE must be one of",
            sprint(showerror, invalid_read_type),
        )

        standalone = _autocycler_smoke_prerequisites(Dict(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_AUTOCYCLER_LONG_READS" => input_paths.long_reads,
            "MYCELIA_AUTOCYCLER_SHORT_READS_1" => input_paths.short_r1,
            "MYCELIA_AUTOCYCLER_SHORT_READS_2" => input_paths.short_r2,
            "MYCELIA_AUTOCYCLER_READ_TYPE" => "pacbio_hifi",
        ))
        Test.@test standalone.run_smoke
        Test.@test standalone.read_type == "pacbio_hifi"
        Test.@test standalone.threads == 2
    end
end
