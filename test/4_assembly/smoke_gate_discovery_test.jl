# Regression coverage for the canonical Pkg.test entrypoint's smoke preflight.

import Test

const _SMOKE_DISCOVERY_ENV_NAMES = (
    "MYCELIA_RUN_ALL",
    "MYCELIA_RUN_EXTERNAL",
    "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE",
    "MYCELIA_RUN_AUTOCYCLER_POLISHED",
    "MYCELIA_HYBRID_SHORT_R1",
    "MYCELIA_HYBRID_SHORT_R2",
    "MYCELIA_HYBRID_LONG_READS",
    "MYCELIA_AUTOCYCLER_READ_TYPE",
    "MYCELIA_AUTOCYCLER_TEST_JOBS",
)

function _run_smoke_discovery(
    overrides::Pair{String, String}...,
)::NamedTuple
    project_root = Base.normpath(Base.joinpath(@__DIR__, "..", ".."))
    runtests_path = Base.joinpath(project_root, "test", "runtests.jl")
    runner = `$(Base.julia_cmd()) --startup-file=no
              --project=$(project_root) $(runtests_path)`
    environment = Dict{String, String}(
        String(name) => String(value) for (name, value) in ENV
    )
    for name in _SMOKE_DISCOVERY_ENV_NAMES
        environment[name] = startswith(name, "MYCELIA_RUN_") ? "false" : ""
    end
    environment["MYCELIA_TEST_DISCOVERY_ONLY"] = "true"
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
        broad_only = _run_smoke_discovery(broad_gate => "true")
        Test.@test broad_only.ok
        Test.@test occursin(
            "Mycelia smoke discovery preflight passed.",
            broad_only.output,
        )
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
        valid_opt_in = _run_smoke_discovery(
            "MYCELIA_RUN_EXTERNAL" => "true",
            "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE" => "true",
            "MYCELIA_HYBRID_SHORT_R1" => input_paths.short_r1,
            "MYCELIA_HYBRID_SHORT_R2" => input_paths.short_r2,
            "MYCELIA_HYBRID_LONG_READS" => input_paths.long_reads,
        )
        Test.@test valid_opt_in.ok
        Test.@test occursin(
            "Mycelia smoke discovery preflight passed.",
            valid_opt_in.output,
        )
    end
end
