# Shared private-fixture gate parsing for the multi-input hybrid smoke tests.

const _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS = (
    "MYCELIA_HYBRID_SHORT_R1",
    "MYCELIA_HYBRID_SHORT_R2",
    "MYCELIA_HYBRID_LONG_READS",
)

function _multi_input_hybrid_smoke_env_enabled(
        environment::AbstractDict,
        name::AbstractString,
)::Bool
    return lowercase(strip(String(get(environment, name, "false")))) == "true"
end

function _multi_input_hybrid_smoke_symbol(
        environment::AbstractDict,
        name::AbstractString,
        default::Symbol,
        allowed::Tuple,
)::Symbol
    configured = get(environment, name, String(default))
    value = Symbol(lowercase(strip(String(configured))))
    value in allowed || throw(ArgumentError(
        "$(name) must be one of $(allowed), got :$(value).",
    ))
    return value
end

function _multi_input_hybrid_smoke_integer(
        environment::AbstractDict,
        name::AbstractString,
        default::Int,
        minimum::Int,
        maximum::Int,
)::Int
    text = strip(String(get(environment, name, string(default))))
    value = tryparse(Int, text)
    value === nothing && throw(ArgumentError(
        "$(name) must be an integer, got $(repr(text)).",
    ))
    minimum <= value <= maximum || throw(ArgumentError(
        "$(name) must be between $(minimum) and $(maximum), got $(value).",
    ))
    return value
end

function _multi_input_hybrid_smoke_prerequisites(
        environment::AbstractDict,
)::NamedTuple
    external_enabled =
        _multi_input_hybrid_smoke_env_enabled(environment, "MYCELIA_RUN_ALL") ||
        _multi_input_hybrid_smoke_env_enabled(
            environment,
            "MYCELIA_RUN_EXTERNAL",
        )
    run_smoke = _multi_input_hybrid_smoke_env_enabled(
        environment,
        "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE",
    )
    run_autocycler = _multi_input_hybrid_smoke_env_enabled(
        environment,
        "MYCELIA_RUN_AUTOCYCLER_POLISHED",
    )
    if (run_smoke || run_autocycler) && !external_enabled
        throw(ArgumentError(
            "Dedicated multi-input hybrid smoke gates also require " *
            "MYCELIA_RUN_EXTERNAL=true (or MYCELIA_RUN_ALL=true).",
        ))
    end
    if run_autocycler && !run_smoke
        throw(ArgumentError(
            "MYCELIA_RUN_AUTOCYCLER_POLISHED=true also requires " *
            "MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE=true.",
        ))
    end
    run_smoke || return (;
        run_smoke = false,
        run_autocycler = false,
    )

    missing_inputs = String[
        name for name in _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS
        if isempty(strip(String(get(environment, name, ""))))
    ]
    isempty(missing_inputs) || throw(ArgumentError(
        "Enabled multi-input hybrid smoke requires all FASTQ variables; " *
        "missing: $(join(missing_inputs, ", ")).",
    ))
    input_paths = (
        short_r1 = abspath(String(environment[
            _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS[1]
        ])),
        short_r2 = abspath(String(environment[
            _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS[2]
        ])),
        long_reads = abspath(String(environment[
            _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS[3]
        ])),
    )
    for (label, path) in pairs(input_paths)
        isfile(path) || throw(ArgumentError(
            "Hybrid smoke $(label) FASTQ not found: $(path)",
        ))
        filesize(path) > 0 || throw(ArgumentError(
            "Hybrid smoke $(label) FASTQ is empty: $(path)",
        ))
    end
    short_read_tech = _multi_input_hybrid_smoke_symbol(
        environment,
        "MYCELIA_HYBRID_SHORT_TECH",
        :illumina,
        (:illumina, :ultima),
    )
    long_read_tech = _multi_input_hybrid_smoke_symbol(
        environment,
        "MYCELIA_HYBRID_LONG_TECH",
        :nanopore,
        (:nanopore, :pacbio_clr, :pacbio_hifi),
    )
    autocycler_read_type = nothing
    autocycler_jobs = nothing
    if run_autocycler
        read_type_text = strip(String(get(
            environment,
            "MYCELIA_AUTOCYCLER_READ_TYPE",
            "",
        )))
        isempty(read_type_text) && throw(ArgumentError(
            "MYCELIA_AUTOCYCLER_READ_TYPE is required when " *
            "MYCELIA_RUN_AUTOCYCLER_POLISHED=true.",
        ))
        autocycler_read_type = _multi_input_hybrid_smoke_symbol(
            environment,
            "MYCELIA_AUTOCYCLER_READ_TYPE",
            :ont_r10,
            (:ont_r9, :ont_r10, :pacbio_clr, :pacbio_hifi),
        )
        autocycler_jobs = _multi_input_hybrid_smoke_integer(
            environment,
            "MYCELIA_AUTOCYCLER_TEST_JOBS",
            1,
            1,
            4,
        )
    end
    return (;
        run_smoke,
        run_autocycler,
        input_paths,
        short_read_tech,
        long_read_tech,
        autocycler_read_type,
        autocycler_jobs,
    )
end
