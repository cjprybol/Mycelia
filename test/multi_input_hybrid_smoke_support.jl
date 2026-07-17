# Shared private-fixture gate and FASTQ parsing for hybrid assembly smoke tests.

import CodecZlib
import FASTX

const _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS = (
    "MYCELIA_HYBRID_SHORT_R1",
    "MYCELIA_HYBRID_SHORT_R2",
    "MYCELIA_HYBRID_LONG_READS",
)

const _AUTOCYCLER_SMOKE_READ_TYPES = (
    :ont_r9,
    :ont_r10,
    :pacbio_clr,
    :pacbio_hifi,
)

function _multi_input_hybrid_smoke_env_enabled(
        environment::AbstractDict,
        name::AbstractString,
)::Bool
    return lowercase(strip(String(get(environment, name, "false")))) == "true"
end

function _multi_input_hybrid_external_suite_enabled(
        environment::AbstractDict,
)::Bool
    return _multi_input_hybrid_smoke_env_enabled(
        environment,
        "MYCELIA_RUN_ALL",
    ) || _multi_input_hybrid_smoke_env_enabled(
        environment,
        "MYCELIA_RUN_EXTERNAL",
    )
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

function _smoke_fastq_reader(
        path::AbstractString,
        label::AbstractString,
)::FASTX.FASTQ.Reader
    normalized_path = abspath(String(path))
    path_base = lowercase(basename(normalized_path))
    compressed = endswith(path_base, ".gz")
    uncompressed_base = compressed ? chop(path_base; tail = 3) : path_base
    occursin(r"\.(fastq|fq)$", uncompressed_base) || throw(ArgumentError(
        "$(label) must use a .fastq or .fq extension, optionally followed by " *
        ".gz: $(normalized_path)",
    ))

    raw_input = open(normalized_path, "r")
    input = try
        compressed ? CodecZlib.GzipDecompressorStream(raw_input) : raw_input
    catch
        close(raw_input)
        rethrow()
    end
    return try
        FASTX.FASTQ.Reader(input)
    catch
        close(input)
        rethrow()
    end
end

function _validate_smoke_fastq(
        path::AbstractString,
        label::AbstractString,
)::Int
    reader = try
        _smoke_fastq_reader(path, label)
    catch caught
        caught isa InterruptException && rethrow()
        throw(ArgumentError(
            "$(label) must be a valid FASTQ file. Cause: " *
            sprint(showerror, caught),
        ))
    end
    record_count = 0
    try
        for record in reader
            record isa FASTX.FASTQ.Record || throw(ArgumentError(
                "$(label) must be a FASTQ file: $(abspath(path))",
            ))
            record_count += 1
        end
    catch caught
        caught isa InterruptException && rethrow()
        throw(ArgumentError(
            "$(label) must be a valid FASTQ file. Cause: " *
            sprint(showerror, caught),
        ))
    finally
        close(reader)
    end
    record_count > 0 || throw(ArgumentError("$(label) must be non-empty."))
    return record_count
end

function _smoke_pair_identifier(identifier::AbstractString)::String
    first_token = first(split(String(identifier)))
    return replace(first_token, r"/[12]$" => "")
end

function _smoke_identifier_pair_role(
        identifier::AbstractString,
)::Union{Nothing, Int}
    first_token = first(split(String(identifier)))
    role_match = match(r"/([12])$", first_token)
    return role_match === nothing ? nothing :
           parse(Int, something(only(role_match.captures)))
end

function _smoke_casava_pair_role(
        description::AbstractString,
)::Union{Nothing, Int}
    description_tokens = split(String(description))
    length(description_tokens) >= 2 || return nothing
    role_match = match(
        r"^([12]):[YN]:[0-9]+:[A-Za-z0-9+_-]+$",
        description_tokens[2],
    )
    return role_match === nothing ? nothing :
           parse(Int, something(only(role_match.captures)))
end

function _smoke_pair_role(
        identifier::AbstractString,
        description::AbstractString = "",
)::Union{Nothing, Int}
    identifier_role = _smoke_identifier_pair_role(identifier)
    casava_role = _smoke_casava_pair_role(description)
    if identifier_role !== nothing && casava_role !== nothing &&
       identifier_role != casava_role
        throw(ArgumentError(
            "FASTQ identifier and CASAVA description contain conflicting " *
            "explicit mate roles: identifier=$(repr(String(identifier))), " *
            "description=$(repr(String(description))).",
        ))
    end
    return identifier_role === nothing ? casava_role : identifier_role
end

function _validate_smoke_paired_fastqs(
        short_reads_1::AbstractString,
        short_reads_2::AbstractString;
        smoke_label::AbstractString = "Autocycler smoke",
)::Int
    !Base.Filesystem.samefile(short_reads_1, short_reads_2) || throw(
        ArgumentError(
            "$(smoke_label) paired short-read R1 and R2 must be distinct files.",
        ),
    )
    reader_1 = _smoke_fastq_reader(
        short_reads_1,
        "$(smoke_label) paired short-read R1 FASTQ",
    )
    reader_2 = try
        _smoke_fastq_reader(
            short_reads_2,
            "$(smoke_label) paired short-read R2 FASTQ",
        )
    catch
        close(reader_1)
        rethrow()
    end
    pair_count = 0
    try
        next_1 = iterate(reader_1)
        next_2 = iterate(reader_2)
        while next_1 !== nothing || next_2 !== nothing
            if next_1 === nothing || next_2 === nothing
                throw(ArgumentError(
                    "$(smoke_label) paired short reads have different counts " *
                    "after $(pair_count) complete pairs.",
                ))
            end
            record_1, state_1 = next_1
            record_2, state_2 = next_2
            pair_count += 1
            identifier_1 = String(FASTX.identifier(record_1))
            identifier_2 = String(FASTX.identifier(record_2))
            description_1 = String(FASTX.description(record_1))
            description_2 = String(FASTX.description(record_2))
            role_1 = _smoke_pair_role(identifier_1, description_1)
            role_2 = _smoke_pair_role(identifier_2, description_2)
            roles_valid = (role_1 === nothing && role_2 === nothing) ||
                          (role_1 == 1 && role_2 == 2)
            roles_valid || throw(ArgumentError(
                "$(smoke_label) paired short reads have invalid explicit mate " *
                "roles at record $(pair_count): R1=$(repr(identifier_1)), " *
                "R2=$(repr(identifier_2)).",
            ))
            _smoke_pair_identifier(identifier_1) ==
            _smoke_pair_identifier(identifier_2) || throw(ArgumentError(
                "$(smoke_label) paired short reads are out of sync at record " *
                "$(pair_count): R1=$(repr(identifier_1)), " *
                "R2=$(repr(identifier_2)).",
            ))
            next_1 = iterate(reader_1, state_1)
            next_2 = iterate(reader_2, state_2)
        end
    catch caught
        caught isa InterruptException && rethrow()
        caught isa ArgumentError && rethrow()
        throw(ArgumentError(
            "$(smoke_label) paired short-read inputs must be valid FASTQ " *
            "files. Cause: $(sprint(showerror, caught))",
        ))
    finally
        close(reader_1)
        close(reader_2)
    end
    pair_count > 0 || throw(ArgumentError(
        "$(smoke_label) paired short reads must be non-empty.",
    ))
    return pair_count
end

function _autocycler_smoke_prerequisites(
        environment::AbstractDict,
)::NamedTuple
    run_smoke = _multi_input_hybrid_smoke_env_enabled(
        environment,
        "MYCELIA_RUN_AUTOCYCLER_SMOKE",
    )
    external_enabled = _multi_input_hybrid_external_suite_enabled(environment)
    if run_smoke && !external_enabled
        throw(ArgumentError(
            "MYCELIA_RUN_AUTOCYCLER_SMOKE=true also requires " *
            "MYCELIA_RUN_EXTERNAL=true (or MYCELIA_RUN_ALL=true).",
        ))
    end
    run_smoke || return (; run_smoke = false)

    threads = _multi_input_hybrid_smoke_integer(
        environment,
        "MYCELIA_ASSEMBLER_TEST_THREADS",
        2,
        1,
        4,
    )
    read_type_text = strip(String(get(
        environment,
        "MYCELIA_AUTOCYCLER_READ_TYPE",
        "",
    )))
    isempty(read_type_text) && throw(ArgumentError(
        "MYCELIA_AUTOCYCLER_READ_TYPE is required when " *
        "MYCELIA_RUN_AUTOCYCLER_SMOKE=true.",
    ))
    read_type = _multi_input_hybrid_smoke_symbol(
        environment,
        "MYCELIA_AUTOCYCLER_READ_TYPE",
        :ont_r10,
        _AUTOCYCLER_SMOKE_READ_TYPES,
    )

    long_reads_text = strip(String(get(
        environment,
        "MYCELIA_AUTOCYCLER_LONG_READS",
        "",
    )))
    isempty(long_reads_text) && throw(ArgumentError(
        "MYCELIA_RUN_AUTOCYCLER_SMOKE=true requires " *
        "MYCELIA_AUTOCYCLER_LONG_READS.",
    ))
    short_reads_1_text = strip(String(get(
        environment,
        "MYCELIA_AUTOCYCLER_SHORT_READS_1",
        "",
    )))
    short_reads_2_text = strip(String(get(
        environment,
        "MYCELIA_AUTOCYCLER_SHORT_READS_2",
        "",
    )))
    xor(isempty(short_reads_1_text), isempty(short_reads_2_text)) && throw(
        ArgumentError(
            "Set both MYCELIA_AUTOCYCLER_SHORT_READS_1 and " *
            "MYCELIA_AUTOCYCLER_SHORT_READS_2, or leave both unset.",
        ),
    )
    long_reads = abspath(long_reads_text)
    isfile(long_reads) || throw(ArgumentError(
        "Autocycler smoke long-read FASTQ not found: $(long_reads)",
    ))
    filesize(long_reads) > 0 || throw(ArgumentError(
        "Autocycler smoke long-read FASTQ is empty: $(long_reads)",
    ))
    _validate_smoke_fastq(long_reads, "Autocycler smoke long-read input")

    short_reads_1 = ""
    short_reads_2 = ""
    if !isempty(short_reads_1_text)
        short_reads_1 = abspath(short_reads_1_text)
        short_reads_2 = abspath(short_reads_2_text)
        for (label, path) in (
                ("paired short-read R1", short_reads_1),
                ("paired short-read R2", short_reads_2),
        )
            isfile(path) || throw(ArgumentError(
                "Autocycler smoke $(label) FASTQ not found: $(path)",
            ))
            filesize(path) > 0 || throw(ArgumentError(
                "Autocycler smoke $(label) FASTQ is empty: $(path)",
            ))
        end
        if Base.Filesystem.samefile(long_reads, short_reads_1) ||
           Base.Filesystem.samefile(long_reads, short_reads_2)
            throw(ArgumentError(
                "Autocycler smoke long reads must be physically distinct from " *
                "paired short-read R1 and R2 files.",
            ))
        end
        _validate_smoke_paired_fastqs(short_reads_1, short_reads_2)
    end
    return (;
        run_smoke,
        long_reads,
        short_reads_1,
        short_reads_2,
        read_type = String(read_type),
        threads,
    )
end

function _multi_input_hybrid_smoke_prerequisites(
        environment::AbstractDict,
)::NamedTuple
    external_enabled = _multi_input_hybrid_external_suite_enabled(environment)
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

    threads = _multi_input_hybrid_smoke_integer(
        environment,
        "MYCELIA_ASSEMBLER_TEST_THREADS",
        2,
        1,
        4,
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
        short_r1 = abspath(strip(String(environment[
            _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS[1]
        ]))),
        short_r2 = abspath(strip(String(environment[
            _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS[2]
        ]))),
        long_reads = abspath(strip(String(environment[
            _MULTI_INPUT_HYBRID_SMOKE_INPUT_ENV_VARS[3]
        ]))),
    )
    for (label, path) in pairs(input_paths)
        isfile(path) || throw(ArgumentError(
            "Hybrid smoke $(label) FASTQ not found: $(path)",
        ))
        filesize(path) > 0 || throw(ArgumentError(
            "Hybrid smoke $(label) FASTQ is empty: $(path)",
        ))
    end
    for (first_label, first_path, second_label, second_path) in (
            (:short_r1, input_paths.short_r1, :short_r2, input_paths.short_r2),
            (
                :short_r1,
                input_paths.short_r1,
                :long_reads,
                input_paths.long_reads,
            ),
            (
                :short_r2,
                input_paths.short_r2,
                :long_reads,
                input_paths.long_reads,
            ),
    )
        !Base.Filesystem.samefile(first_path, second_path) || throw(
            ArgumentError(
                "Hybrid smoke $(first_label) and $(second_label) FASTQ inputs " *
                "must be physically distinct files.",
            ),
        )
    end
    _validate_smoke_fastq(
        input_paths.long_reads,
        "Hybrid smoke long-read input",
    )
    _validate_smoke_paired_fastqs(
        input_paths.short_r1,
        input_paths.short_r2;
        smoke_label = "Hybrid smoke",
    )
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
            _AUTOCYCLER_SMOKE_READ_TYPES,
        )
        autocycler_jobs = _multi_input_hybrid_smoke_integer(
            environment,
            "MYCELIA_AUTOCYCLER_TEST_JOBS",
            1,
            1,
            4,
        )
        compatible = if long_read_tech == :nanopore
            autocycler_read_type in (:ont_r9, :ont_r10)
        elseif long_read_tech == :pacbio_clr
            autocycler_read_type == :pacbio_clr
        else
            autocycler_read_type == :pacbio_hifi
        end
        compatible || throw(ArgumentError(
            "MYCELIA_AUTOCYCLER_READ_TYPE=:$(autocycler_read_type) is " *
            "incompatible with MYCELIA_HYBRID_LONG_TECH=:$(long_read_tech).",
        ))
    end
    return (;
        run_smoke,
        run_autocycler,
        input_paths,
        short_read_tech,
        long_read_tech,
        autocycler_read_type,
        autocycler_jobs,
        threads,
    )
end
