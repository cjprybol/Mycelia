# Regression coverage for the preflight driver used by the canonical Pkg.test
# entrypoint.

import Test

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

Test.@testset "canonical smoke discovery preflight" begin
    for broad_gate in ("MYCELIA_RUN_ALL", "MYCELIA_RUN_EXTERNAL")
        for broad_value in ("true", " true ")
            broad_only = _run_smoke_discovery(broad_gate => broad_value)
            Test.@test broad_only.ok
            Test.@test isempty(strip(broad_only.output))
        end
    end

    for dedicated_gate in (
            "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE",
            "MYCELIA_RUN_AUTOCYCLER_POLISHED",
    )
        dedicated_only = _run_smoke_discovery(dedicated_gate => "true")
        Test.@test !dedicated_only.ok
        Test.@test occursin(
            "also require MYCELIA_RUN_EXTERNAL=true",
            dedicated_only.output,
        )
    end

    standalone_without_broad_gate = _run_smoke_discovery(
        "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
    )
    Test.@test !standalone_without_broad_gate.ok
    Test.@test occursin(
        "MYCELIA_RUN_AUTOCYCLER_SMOKE=true also requires " *
        "MYCELIA_RUN_EXTERNAL=true",
        standalone_without_broad_gate.output,
    )

    standalone_without_fixtures = _run_smoke_discovery(
        "MYCELIA_RUN_EXTERNAL" => "true",
        "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
    )
    Test.@test !standalone_without_fixtures.ok
    Test.@test occursin(
        "requires MYCELIA_AUTOCYCLER_LONG_READS",
        standalone_without_fixtures.output,
    )

    missing_fixtures = _run_smoke_discovery(
        "MYCELIA_RUN_EXTERNAL" => "true",
        "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
    )
    Test.@test !missing_fixtures.ok
    Test.@test occursin("missing:", missing_fixtures.output)

    missing_autocycler_fixtures = _run_smoke_discovery(
        "MYCELIA_RUN_EXTERNAL" => "true",
        "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
        "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
    )
    Test.@test !missing_autocycler_fixtures.ok
    Test.@test occursin("missing:", missing_autocycler_fixtures.output)

    missing_parent_gate = _run_smoke_discovery(
        "MYCELIA_RUN_EXTERNAL" => "true",
        "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
    )
    Test.@test !missing_parent_gate.ok
    Test.@test occursin(
        "also requires MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE=true",
        missing_parent_gate.output,
    )

    Base.mktempdir() do temporary_root
        input_paths = (
            short_r1 = Base.joinpath(temporary_root, "reads_R1.fastq"),
            short_r2 = Base.joinpath(temporary_root, "reads_R2.fastq"),
            long_reads = Base.joinpath(temporary_root, "long.fastq"),
        )
        for path in values(input_paths)
            Base.write(path, "@read\nACGT\n+\nIIII\n")
        end
        for broad_gate in ("MYCELIA_RUN_ALL", "MYCELIA_RUN_EXTERNAL")
            valid_opt_in = _run_smoke_discovery(
                broad_gate => " true ",
                "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
                "MYCELIA_HYBRID_SHORT_R1" => input_paths.short_r1,
                "MYCELIA_HYBRID_SHORT_R2" => input_paths.short_r2,
                "MYCELIA_HYBRID_LONG_READS" => input_paths.long_reads,
            )
            Test.@test valid_opt_in.ok
            Test.@test isempty(strip(valid_opt_in.output))
        end

        standalone_paths = (
            long_reads = Base.joinpath(
                temporary_root,
                "standalone_long.fastq",
            ),
            short_r1 = Base.joinpath(
                temporary_root,
                "standalone_R1.fastq",
            ),
            short_r2 = Base.joinpath(
                temporary_root,
                "standalone_R2.fastq",
            ),
        )
        Base.write(
            standalone_paths.long_reads,
            "@long\nACGT\n+\nIIII\n",
        )
        Base.write(
            standalone_paths.short_r1,
            "@pair/1\nACGT\n+\nIIII\n",
        )
        Base.write(
            standalone_paths.short_r2,
            "@pair/2\nACGT\n+\nIIII\n",
        )

        half_pair = _run_smoke_discovery(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_AUTOCYCLER_LONG_READS" =>
                standalone_paths.long_reads,
            "MYCELIA_AUTOCYCLER_SHORT_READS_1" =>
                standalone_paths.short_r1,
        )
        Test.@test !half_pair.ok
        Test.@test occursin(
            "Set both MYCELIA_AUTOCYCLER_SHORT_READS_1",
            half_pair.output,
        )

        malformed_long_reads = Base.joinpath(
            temporary_root,
            "malformed_long.fastq",
        )
        Base.write(malformed_long_reads, ">not_fastq\nACGT\n")
        malformed_fixture = _run_smoke_discovery(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_AUTOCYCLER_LONG_READS" => malformed_long_reads,
        )
        Test.@test !malformed_fixture.ok
        Test.@test occursin("must be a valid FASTQ file", malformed_fixture.output)

        invalid_pair_r1 = Base.joinpath(
            temporary_root,
            "invalid_pair_R1.fastq",
        )
        Base.write(invalid_pair_r1, "@pair/2\nACGT\n+\nIIII\n")
        invalid_pair = _run_smoke_discovery(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_AUTOCYCLER_LONG_READS" =>
                standalone_paths.long_reads,
            "MYCELIA_AUTOCYCLER_SHORT_READS_1" => invalid_pair_r1,
            "MYCELIA_AUTOCYCLER_SHORT_READS_2" =>
                standalone_paths.short_r2,
        )
        Test.@test !invalid_pair.ok
        Test.@test occursin("invalid explicit mate roles", invalid_pair.output)

        missing_pair_fixture = _run_smoke_discovery(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_AUTOCYCLER_LONG_READS" =>
                standalone_paths.long_reads,
            "MYCELIA_AUTOCYCLER_SHORT_READS_1" =>
                standalone_paths.short_r1,
            "MYCELIA_AUTOCYCLER_SHORT_READS_2" => Base.joinpath(
                temporary_root,
                "missing_R2.fastq",
            ),
        )
        Test.@test !missing_pair_fixture.ok
        Test.@test occursin(
            "paired short-read R2 FASTQ not found",
            missing_pair_fixture.output,
        )

        invalid_standalone_read_type = _run_smoke_discovery(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_AUTOCYCLER_LONG_READS" =>
                standalone_paths.long_reads,
            "MYCELIA_AUTOCYCLER_READ_TYPE" => "illumina",
        )
        Test.@test !invalid_standalone_read_type.ok
        Test.@test occursin(
            "MYCELIA_AUTOCYCLER_READ_TYPE must be one of",
            invalid_standalone_read_type.output,
        )

        valid_standalone_long_only = _run_smoke_discovery(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_AUTOCYCLER_LONG_READS" =>
                standalone_paths.long_reads,
        )
        Test.@test valid_standalone_long_only.ok
        Test.@test isempty(strip(valid_standalone_long_only.output))

        valid_standalone_pair = _run_smoke_discovery(
            "MYCELIA_RUN_ALL" => " true ",
            "MYCELIA_RUN_AUTOCYCLER_SMOKE" => "true",
            "MYCELIA_AUTOCYCLER_LONG_READS" =>
                standalone_paths.long_reads,
            "MYCELIA_AUTOCYCLER_SHORT_READS_1" =>
                standalone_paths.short_r1,
            "MYCELIA_AUTOCYCLER_SHORT_READS_2" =>
                standalone_paths.short_r2,
            "MYCELIA_AUTOCYCLER_READ_TYPE" => "pacbio_hifi",
        )
        Test.@test valid_standalone_pair.ok
        Test.@test isempty(strip(valid_standalone_pair.output))

        incompatible_read_types = (
            :nanopore => (:pacbio_clr, :pacbio_hifi),
            :pacbio_clr => (:ont_r9, :ont_r10, :pacbio_hifi),
            :pacbio_hifi => (:ont_r9, :ont_r10, :pacbio_clr),
        )
        for (long_read_tech, read_types) in incompatible_read_types
            for read_type in read_types
                incompatible = _run_smoke_discovery(
                    "MYCELIA_RUN_EXTERNAL" => "true",
                    "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
                    "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
                    "MYCELIA_HYBRID_SHORT_R1" => input_paths.short_r1,
                    "MYCELIA_HYBRID_SHORT_R2" => input_paths.short_r2,
                    "MYCELIA_HYBRID_LONG_READS" => input_paths.long_reads,
                    "MYCELIA_HYBRID_LONG_TECH" => String(long_read_tech),
                    "MYCELIA_AUTOCYCLER_READ_TYPE" => String(read_type),
                )
                Test.@test !incompatible.ok
                Test.@test occursin(
                    "MYCELIA_AUTOCYCLER_READ_TYPE=:$(read_type) is " *
                    "incompatible with MYCELIA_HYBRID_LONG_TECH=" *
                    ":$(long_read_tech)",
                    incompatible.output,
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
                compatible = _run_smoke_discovery(
                    "MYCELIA_RUN_EXTERNAL" => "true",
                    "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
                    "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
                    "MYCELIA_HYBRID_SHORT_R1" => input_paths.short_r1,
                    "MYCELIA_HYBRID_SHORT_R2" => input_paths.short_r2,
                    "MYCELIA_HYBRID_LONG_READS" => input_paths.long_reads,
                    "MYCELIA_HYBRID_LONG_TECH" => String(long_read_tech),
                    "MYCELIA_AUTOCYCLER_READ_TYPE" => String(read_type),
                )
                Test.@test compatible.ok
                Test.@test isempty(strip(compatible.output))
            end
        end
    end
end
