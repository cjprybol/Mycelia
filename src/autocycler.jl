"""
Wrapper for Autocycler and a conservative paired-short polishing pipeline.

References:
- Autocycler: https://github.com/rrwick/Autocycler
- Polypolish: https://github.com/rrwick/Polypolish
- Pypolca: https://github.com/gbouras13/pypolca
"""

const AUTOCYCLER_VERSION = "0.5.2"
const AUTOCYCLER_ENVIRONMENT_SPEC_SHA256 =
    "d6aef758986db23c453203f067a23be6f29b690b536bc24b283b5a6e913624ef"
const AUTOCYCLER_ENV_NAME =
    "autocycler-$(AUTOCYCLER_VERSION)-$(first(AUTOCYCLER_ENVIRONMENT_SPEC_SHA256, 16))"
const AUTOCYCLER_MAX_ASSEMBLY_THREADS = 128
const AUTOCYCLER_INSTALL_LOCK_STALE_SECONDS = 6 * 60 * 60
const AUTOCYCLER_OUTPUT_LOCK_STALE_SECONDS = 7 * 24 * 60 * 60
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
    (name = "autocycler", constraint = :exact, version = AUTOCYCLER_VERSION),
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
const AUTOCYCLER_PACKAGE_INVENTORY_SCHEMA =
    "conda-name-version-build-channel-v1"
const AUTOCYCLER_TOOLCHAIN_SCHEMA = "mycelia-autocycler-toolchain-v1"
const AUTOCYCLER_READ_TYPES = (
    "ont_r9",
    "ont_r10",
    "pacbio_clr",
    "pacbio_hifi",
)
const AUTOCYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING = 1_000_000_000_000

function _autocycler_paths()::Tuple{String, String, String}
    install_dir = joinpath(dirname(dirname(pathof(Mycelia))), "deps", "autocycler")
    script_path = joinpath(install_dir, "autocycler_full.sh")
    env_file_path = joinpath(install_dir, "environment.yml")
    return install_dir, script_path, env_file_path
end

function _canonical_autocycler_conda_runner(
        conda_runner::AbstractString = _conda_runner(),
)::String
    executable = Sys.which(String(conda_runner))
    candidate = executable === nothing ?
                abspath(String(conda_runner)) : String(executable)
    return ispath(candidate) ? realpath(candidate) : normpath(candidate)
end

function _canonical_autocycler_environment_prefix(
        environment_prefix::AbstractString,
)::String
    normalized_prefix = normpath(abspath(String(environment_prefix)))
    existing_ancestor = normalized_prefix
    missing_components = String[]
    while !ispath(existing_ancestor) && !islink(existing_ancestor)
        parent = dirname(existing_ancestor)
        parent == existing_ancestor && break
        pushfirst!(missing_components, basename(existing_ancestor))
        existing_ancestor = parent
    end
    isdir(existing_ancestor) || throw(ArgumentError(
        "Autocycler environment prefix has a non-directory existing " *
        "ancestor: $(existing_ancestor)",
    ))
    canonical_ancestor = realpath(existing_ancestor)
    return isempty(missing_components) ? canonical_ancestor :
           normpath(joinpath(canonical_ancestor, missing_components...))
end

function _autocycler_environment_prefix(
        conda_runner::AbstractString = _conda_runner(),
)::String
    canonical_runner = _canonical_autocycler_conda_runner(conda_runner)
    conda_root = normpath(joinpath(dirname(canonical_runner), ".."))
    return _canonical_autocycler_environment_prefix(
        joinpath(conda_root, "envs", AUTOCYCLER_ENV_NAME),
    )
end

function _autocycler_install_lock_path(
        conda_runner::AbstractString = _conda_runner(),
)::String
    environment_prefix = _autocycler_environment_prefix(conda_runner)
    return _autocycler_install_lock_path_from_prefix(environment_prefix)
end

function _autocycler_install_lock_path_from_prefix(
        environment_prefix::AbstractString,
)::String
    normalized_prefix =
        _canonical_autocycler_environment_prefix(environment_prefix)
    conda_root = dirname(dirname(normalized_prefix))
    return joinpath(
        conda_root,
        ".mycelia-locks",
        "$(AUTOCYCLER_ENV_NAME).pid",
    )
end

function _autocycler_sha256(path::AbstractString)::String
    return open(path, "r") do input
        SHA.bytes2hex(SHA.sha256(input))
    end
end

function _autocycler_sha256(bytes::AbstractVector{UInt8})::String
    return SHA.bytes2hex(SHA.sha256(bytes))
end

function _autocycler_spooled_input_snapshot(
        path::AbstractString,
        label::AbstractString,
        ;
        spool_parent::Union{Nothing, AbstractString} = nothing,
        byte_ceiling::Integer = AUTOCYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
        available_bytes_reader::Function = parent -> diskstat(parent).available,
        after_copy_hook::Function = binding -> nothing,
)::NamedTuple
    bindings = _autocycler_spooled_input_snapshots(
        (path,),
        (label,);
        spool_parent,
        byte_ceiling,
        available_bytes_reader,
        after_copy_hook,
    )
    return only(bindings.inputs)
end

function _autocycler_input_source_snapshot(
        path::AbstractString,
        label::AbstractString,
)::NamedTuple
    normalized_path = _require_nonempty_autocycler_file(path, label)
    canonical_path = realpath(normalized_path)
    status = stat(normalized_path)
    return (;
        path = normalized_path,
        canonical_path,
        size_bytes = filesize(normalized_path),
        device = UInt64(status.device),
        inode = UInt64(status.inode),
    )
end

function _autocycler_spool_root_identity(
        spool_root::AbstractString,
)::NamedTuple
    normalized_root = normpath(abspath(String(spool_root)))
    isdir(normalized_root) && !islink(normalized_root) || error(
        "Autocycler private input spool root is missing, replaced, or not " *
        "a directory: $(normalized_root).",
    )
    status = stat(normalized_root)
    return (;
        path = normalized_root,
        device = UInt64(status.device),
        inode = UInt64(status.inode),
    )
end

function _require_unchanged_autocycler_spool_root(
        expected::NamedTuple,
)::Nothing
    _autocycler_spool_root_identity(expected.path) == expected || error(
        "Autocycler private input spool root changed physical identity.",
    )
    return nothing
end

function _autocycler_copy_input_snapshot!(
        source::NamedTuple,
        spool_path::AbstractString,
        label::AbstractString;
        after_copy_hook::Function = binding -> nothing,
)::NamedTuple
    input = Base.Filesystem.open(
        source.canonical_path,
        Base.JL_O_RDONLY |
        Base.JL_O_NOFOLLOW |
        Base.JL_O_CLOEXEC,
    )
    output = nothing
    destination_identity = nothing
    try
        input_status = stat(input)
        UInt64(input_status.device) == source.device &&
            UInt64(input_status.inode) == source.inode &&
            filesize(input) == source.size_bytes || throw(ErrorException(
                "$(label) changed before its private stable input snapshot " *
                "was copied.",
            ))
        output = Base.Filesystem.open(
            spool_path,
            Base.JL_O_WRONLY |
            Base.JL_O_CREAT |
            Base.JL_O_EXCL |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
            0o600,
        )
        output_status = stat(output)
        destination_identity = (;
            device = UInt64(output_status.device),
            inode = UInt64(output_status.inode),
        )
        buffer = Vector{UInt8}(undef, 1024 * 1024)
        remaining = source.size_bytes
        while remaining > 0
            requested = min(length(buffer), remaining)
            count = readbytes!(input, buffer, requested)
            count > 0 || throw(ErrorException(
                "$(label) shrank while its private stable input snapshot " *
                "was materialized.",
            ))
            write(output, @view buffer[1:count])
            remaining -= count
        end
        eof(input) || throw(ErrorException(
            "$(label) grew while its private stable input snapshot was " *
            "materialized.",
        ))
    finally
        output === nothing || close(output)
        close(input)
    end
    try
        after_copy_hook((; source, spool_path = String(spool_path), label))
        observed = _autocycler_input_source_snapshot(source.path, label)
        observed == source || throw(ErrorException(
            "$(label) changed physical identity or size while its private " *
            "stable input snapshot was materialized.",
        ))
        isfile(spool_path) && !islink(spool_path) || throw(ErrorException(
            "$(label) private stable input snapshot was replaced after " *
            "copying.",
        ))
        spool_status = stat(spool_path)
        destination_identity == (;
            device = UInt64(spool_status.device),
            inode = UInt64(spool_status.inode),
        ) || throw(ErrorException(
            "$(label) private stable input snapshot changed physical identity.",
        ))
        chmod(spool_path, 0o400)
        filesize(spool_path) > 0 || throw(ArgumentError(
            "$(label) became empty while its private semantic snapshot was " *
            "materialized: $(source.path)",
        ))
        consumed_snapshot = (;
            path = String(spool_path),
            size_bytes = filesize(spool_path),
            sha256 = _autocycler_sha256(spool_path),
        )
        snapshot = (
            path = source.path,
            canonical_path = source.canonical_path,
            size_bytes = source.size_bytes,
            sha256 = consumed_snapshot.sha256,
            consumed_path = consumed_snapshot.path,
            consumed_size_bytes = consumed_snapshot.size_bytes,
            consumed_sha256 = consumed_snapshot.sha256,
        )
        return (; snapshot, consumed_snapshot, semantic_path = String(spool_path))
    catch
        rethrow()
    end
end

function _autocycler_spooled_input_snapshots(
        paths::Tuple,
        labels::Tuple;
        spool_parent::Union{Nothing, AbstractString} = nothing,
        byte_ceiling::Integer = AUTOCYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
        available_bytes_reader::Function = parent -> diskstat(parent).available,
        after_copy_hook::Function = binding -> nothing,
)::NamedTuple
    length(paths) == length(labels) || throw(ArgumentError(
        "Autocycler input snapshot paths and labels must have equal length.",
    ))
    byte_ceiling > 0 || throw(ArgumentError(
        "Autocycler input spool byte ceiling must be positive, got " *
        "$(byte_ceiling).",
    ))
    sources = NamedTuple[
        _autocycler_input_source_snapshot(path, label) for
        (path, label) in zip(paths, labels)
    ]
    total_bytes_exact = sum(
        BigInt(source.size_bytes) for source in sources;
        init = BigInt(0),
    )
    total_bytes_exact <= BigInt(byte_ceiling) || throw(ArgumentError(
        "Autocycler inputs require $(total_bytes_exact) spool bytes, exceeding the " *
        "configured cumulative ceiling of $(byte_ceiling) bytes.",
    ))
    total_bytes = Int(total_bytes_exact)
    resolved_parent = spool_parent === nothing ? tempdir() :
                      normpath(abspath(String(spool_parent)))
    isdir(resolved_parent) && !islink(resolved_parent) || throw(ArgumentError(
        "Autocycler input spool parent must be an existing non-symlink " *
        "directory: $(resolved_parent).",
    ))
    available_bytes = Int(available_bytes_reader(resolved_parent))
    available_bytes >= total_bytes || throw(ArgumentError(
        "Autocycler input spool preflight requires $(total_bytes) bytes but " *
        "only $(available_bytes) bytes are available under " *
        "$(resolved_parent).",
    ))
    spool_root = mktempdir(
        resolved_parent;
        prefix = "mycelia-autocycler-input-",
    )
    chmod(spool_root, 0o700)
    spool_root_identity = _autocycler_spool_root_identity(spool_root)
    inputs = NamedTuple[]
    try
        for (index, (source, label)) in enumerate(zip(sources, labels))
            _require_unchanged_autocycler_spool_root(spool_root_identity)
            extension = endswith(lowercase(source.path), ".gz") ?
                        ".fastq.gz" : ".fastq"
            spool_path = joinpath(spool_root, "input_$(index)$(extension)")
            copied = _autocycler_copy_input_snapshot!(
                source,
                spool_path,
                label;
                after_copy_hook,
            )
            _require_unchanged_autocycler_spool_root(spool_root_identity)
            push!(
                inputs,
                merge(copied, (; spool_root, spool_root_identity)),
            )
        end
        return (; inputs = Tuple(inputs), spool_root, total_bytes)
    catch caught
        _remove_exact_private_input_spool_root!(
            spool_root_identity,
            "Autocycler",
        )
        caught isa InterruptException && rethrow()
        if caught isa Base.IOError || caught isa SystemError
            throw(ErrorException(
                "Autocycler input spooling failed after cleanup, possibly " *
                "because scratch space was exhausted: " *
                sprint(showerror, caught),
            ))
        end
        rethrow()
    end
end

function _cleanup_autocycler_input_spool!(binding::NamedTuple)::Nothing
    _remove_exact_private_input_spool_root!(
        binding.spool_root_identity,
        "Autocycler",
    )
    return nothing
end

function _with_autocycler_spooled_input_snapshot(
        function_to_run::Function,
        path::AbstractString,
        label::AbstractString,
)::Any
    binding = _autocycler_spooled_input_snapshot(path, label)
    try
        return function_to_run(binding)
    finally
        _cleanup_autocycler_input_spool!(binding)
    end
end

function _with_autocycler_spooled_input_snapshots(
        action::Function,
        paths::Tuple,
        labels::Tuple;
        kwargs...,
)::Any
    bindings = _autocycler_spooled_input_snapshots(paths, labels; kwargs...)
    try
        return action(bindings.inputs)
    finally
        _cleanup_autocycler_input_spool!(first(bindings.inputs))
    end
end

function _autocycler_input_snapshot(
        path::AbstractString,
        label::AbstractString,
)::NamedTuple
    normalized_path = _require_nonempty_autocycler_file(path, label)
    return (
        path = normalized_path,
        canonical_path = realpath(normalized_path),
        size_bytes = filesize(normalized_path),
        sha256 = _autocycler_sha256(normalized_path),
    )
end

function _require_unchanged_autocycler_input(
        expected_snapshot::NamedTuple,
        label::AbstractString,
)::NamedTuple
    observed_snapshot = try
        _autocycler_input_snapshot(expected_snapshot.path, label)
    catch caught
        caught isa InterruptException && rethrow()
        throw(ErrorException(
            "$(label) changed after its initial path/size/SHA-256 snapshot. " *
            "Cause: $(sprint(showerror, caught))",
        ))
    end
    expected_original = (;
        path = expected_snapshot.path,
        canonical_path = expected_snapshot.canonical_path,
        size_bytes = expected_snapshot.size_bytes,
        sha256 = expected_snapshot.sha256,
    )
    observed_snapshot == expected_original || throw(ErrorException(
        "$(label) changed after its initial path/size/SHA-256 snapshot: " *
        "expected canonical_path=$(expected_snapshot.canonical_path), " *
        "size=$(expected_snapshot.size_bytes), " *
        "sha256=$(expected_snapshot.sha256); observed " *
        "canonical_path=$(observed_snapshot.canonical_path), " *
        "size=$(observed_snapshot.size_bytes), " *
        "sha256=$(observed_snapshot.sha256).",
    ))
    return observed_snapshot
end


function _require_unchanged_autocycler_consumed_input(
        expected_snapshot::NamedTuple,
        label::AbstractString,
)::NamedTuple
    path = expected_snapshot.path
    isfile(path) && !islink(path) || throw(ErrorException(
        "$(label) stable consumed snapshot is missing or not a regular " *
        "non-symlink file: $(path).",
    ))
    observed = (;
        path = String(path),
        size_bytes = filesize(path),
        sha256 = _autocycler_sha256(path),
    )
    observed == expected_snapshot || throw(ErrorException(
        "$(label) stable consumed snapshot changed before or during tool " *
        "consumption.",
    ))
    return observed
end

function _require_verified_autocycler_environment_spec(
        path::AbstractString,
)::String
    normalized_path = abspath(path)
    if !isfile(normalized_path) || filesize(normalized_path) == 0
        throw(
            ErrorException(
                "Bundled Autocycler environment file is missing or empty: " *
                "$(normalized_path). Reinstall Mycelia before installing " *
                "Autocycler.",
            ),
        )
    end
    actual_sha256 = _autocycler_sha256(normalized_path)
    if actual_sha256 != AUTOCYCLER_ENVIRONMENT_SPEC_SHA256
        throw(
            ErrorException(
                "Autocycler environment specification checksum mismatch: " *
                "expected $(AUTOCYCLER_ENVIRONMENT_SPEC_SHA256), got " *
                "$(actual_sha256) for $(normalized_path).",
            ),
        )
    end
    return normalized_path
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

function _autocycler_package_record_field(
        package_record::Union{NamedTuple, AbstractDict},
        field::Symbol,
)::Any
    if package_record isa NamedTuple
        return hasproperty(package_record, field) ?
               getproperty(package_record, field) : nothing
    end
    return get(
        package_record,
        String(field),
        get(package_record, field, nothing),
    )
end

function _normalize_autocycler_package_inventory(
        package_records::AbstractVector,
)::Vector{NamedTuple}
    isempty(package_records) && throw(
        ErrorException(
            "Autocycler Conda package inventory must contain at least one package.",
        ),
    )
    inventory = NamedTuple[]
    for (record_index, package_record) in enumerate(package_records)
        package_record isa Union{NamedTuple, AbstractDict} || throw(
            ErrorException(
                "Autocycler Conda package inventory record $(record_index) " *
                "is not an object.",
            ),
        )
        name = _autocycler_package_record_field(package_record, :name)
        version = _autocycler_package_record_field(package_record, :version)
        build = _autocycler_package_record_field(package_record, :build_string)
        build === nothing &&
            (build = _autocycler_package_record_field(package_record, :build))
        channel = _autocycler_package_record_field(package_record, :channel)
        fields = (; name, version, build, channel)
        for (field_name, value) in pairs(fields)
            value isa AbstractString && !isempty(value) || throw(
                ErrorException(
                    "Autocycler Conda package inventory record " *
                    "$(record_index) has a missing or empty $(field_name) field.",
                ),
            )
            occursin(r"[\t\r\n]", value) && throw(
                ErrorException(
                    "Autocycler Conda package inventory record " *
                    "$(record_index) has a noncanonical $(field_name) field.",
                ),
            )
        end
        push!(inventory, (;
            name = String(name),
            version = String(version),
            build = String(build),
            channel = String(channel),
        ))
    end
    sort!(
        inventory;
        by = record -> (
            record.name,
            record.version,
            record.build,
            record.channel,
        ),
    )
    names = getproperty.(inventory, :name)
    allunique(names) || throw(
        ErrorException(
            "Autocycler Conda package inventory contains duplicate package names.",
        ),
    )
    return inventory
end

function _normalize_autocycler_package_inventory(
        package_records::Any,
)::Vector{NamedTuple}
    throw(
        ErrorException(
            "Autocycler Conda package inventory was not a JSON array.",
        ),
    )
end

function _autocycler_package_inventory_contents(
        package_records::AbstractVector,
)::String
    inventory = _normalize_autocycler_package_inventory(package_records)
    return join(
        map(inventory) do record
            join(
                (record.name, record.version, record.build, record.channel),
                '\t',
            )
        end,
        '\n',
    ) * "\n"
end

function _autocycler_package_inventory_sha256(
        package_records::AbstractVector,
)::String
    return SHA.bytes2hex(
        SHA.sha256(_autocycler_package_inventory_contents(package_records)),
    )
end

function _autocycler_environment_packages(;
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::AbstractString =
            _autocycler_environment_prefix(conda_runner),
        command_reader::Function = command -> read(command, String),
)::Vector{NamedTuple}
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix =
        _canonical_autocycler_environment_prefix(environment_prefix)
    command = Cmd(
        String[
            resolved_runner,
            "list",
            "-p",
            resolved_prefix,
            "--json",
        ],
    )
    package_records = JSON.parse(command_reader(command))
    package_records isa AbstractVector || throw(
        ErrorException("Conda package inventory was not a JSON array."),
    )
    return _normalize_autocycler_package_inventory(package_records)
end

function _missing_autocycler_packages(
        package_records::AbstractVector,
)::Vector{String}
    inventory = _normalize_autocycler_package_inventory(package_records)
    names = Set(record.name for record in inventory)
    return String[
        package for package in AUTOCYCLER_REQUIRED_PACKAGES if
        !(package in names)
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
        package_records::AbstractVector,
)::Vector{String}
    inventory = _normalize_autocycler_package_inventory(package_records)
    versions = Dict(record.name => record.version for record in inventory)
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
)::Vector{NamedTuple}
    inventory = _normalize_autocycler_package_inventory(package_inspector())
    package_issues = _autocycler_package_issues(inventory)
    isempty(package_issues) || throw(
        ErrorException(
            "Autocycler immutable spec-hash-addressed environment " *
            "$(repr(AUTOCYCLER_ENV_NAME)) has missing or incompatible " *
            "required packages: $(join(package_issues, "; ")). Refusing " *
            "to recreate or repair it in place. Stop all workflows using " *
            "this environment and remove it manually before reinstalling, " *
            "or use a Mycelia release with a corrected environment " *
            "specification.",
        ),
    )
    return inventory
end

function _autocycler_toolchain_metadata(
        package_records::AbstractVector,
)::Dict{String, Any}
    inventory = _normalize_autocycler_package_inventory(package_records)
    package_issues = _autocycler_package_issues(inventory)
    package_issue_summary = join(package_issues, "; ")
    isempty(package_issues) || throw(
        ErrorException(
            "Autocycler toolchain provenance contains incompatible required " *
            "packages: $(package_issue_summary).",
        ),
    )
    _, script_path, env_file_path = _autocycler_paths()
    verified_env_file_path =
        _require_verified_autocycler_environment_spec(env_file_path)
    return Dict{String, Any}(
        "toolchain_schema" => AUTOCYCLER_TOOLCHAIN_SCHEMA,
        "autocycler_script_revision" => AUTOCYCLER_SCRIPT_REVISION,
        "autocycler_script_sha256" => _autocycler_sha256(script_path),
        "environment_name" => AUTOCYCLER_ENV_NAME,
        "environment_spec_expected_sha256" =>
            AUTOCYCLER_ENVIRONMENT_SPEC_SHA256,
        "environment_spec_sha256" =>
            _autocycler_sha256(verified_env_file_path),
        "inventory_schema" => AUTOCYCLER_PACKAGE_INVENTORY_SCHEMA,
        "package_inventory_sha256" =>
            _autocycler_package_inventory_sha256(inventory),
        "package_count" => length(inventory),
        "packages" => Dict{String, Any}[
            Dict{String, Any}(
                "name" => record.name,
                "version" => record.version,
                "build" => record.build,
                "channel" => record.channel,
            ) for record in inventory
        ],
    )
end

function _require_autocycler_toolchain_provenance(
        toolchain::Any,
)::Dict{String, Any}
    toolchain isa AbstractDict || throw(
        ErrorException(
            "Autocycler workflow did not report realized toolchain provenance.",
        ),
    )
    required_scalars = Dict{String, String}(
        "toolchain_schema" => AUTOCYCLER_TOOLCHAIN_SCHEMA,
        "autocycler_script_revision" => AUTOCYCLER_SCRIPT_REVISION,
        "autocycler_script_sha256" => AUTOCYCLER_SCRIPT_SHA256,
        "environment_name" => AUTOCYCLER_ENV_NAME,
        "environment_spec_expected_sha256" =>
            AUTOCYCLER_ENVIRONMENT_SPEC_SHA256,
        "environment_spec_sha256" => AUTOCYCLER_ENVIRONMENT_SPEC_SHA256,
        "inventory_schema" => AUTOCYCLER_PACKAGE_INVENTORY_SCHEMA,
    )
    for (field, expected_value) in required_scalars
        actual_value = get(toolchain, field, nothing)
        actual_value == expected_value || throw(
            ErrorException(
                "Autocycler workflow toolchain provenance has incompatible " *
                "$(field): expected $(repr(expected_value)), got " *
                "$(repr(actual_value)).",
            ),
        )
    end
    packages = get(toolchain, "packages", nothing)
    packages isa AbstractVector || throw(
        ErrorException(
            "Autocycler workflow toolchain provenance has no realized package " *
            "inventory.",
        ),
    )
    inventory = _normalize_autocycler_package_inventory(packages)
    package_issues = _autocycler_package_issues(inventory)
    package_issue_summary = join(package_issues, "; ")
    isempty(package_issues) || throw(
        ErrorException(
            "Autocycler workflow toolchain provenance has incompatible " *
            "required packages: $(package_issue_summary).",
        ),
    )
    package_count = get(toolchain, "package_count", nothing)
    package_count isa Integer && !(package_count isa Bool) &&
        package_count == length(inventory) || throw(
        ErrorException(
            "Autocycler workflow toolchain package count does not match its " *
            "realized package inventory.",
        ),
    )
    digest = get(toolchain, "package_inventory_sha256", nothing)
    digest isa AbstractString && occursin(r"^[0-9a-f]{64}$", digest) || throw(
        ErrorException(
            "Autocycler workflow toolchain provenance is missing a valid " *
            "package inventory SHA-256 digest.",
        ),
    )
    expected_digest = _autocycler_package_inventory_sha256(inventory)
    digest == expected_digest || throw(
        ErrorException(
            "Autocycler workflow package inventory digest does not match its " *
            "reported name/version/build/channel inventory.",
        ),
    )
    expected_fields = Set([
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
    Set(String.(keys(toolchain))) == expected_fields || throw(
        ErrorException(
            "Autocycler workflow toolchain provenance has unexpected fields.",
        ),
    )
    return Dict{String, Any}(
        "toolchain_schema" => AUTOCYCLER_TOOLCHAIN_SCHEMA,
        "autocycler_script_revision" => AUTOCYCLER_SCRIPT_REVISION,
        "autocycler_script_sha256" => AUTOCYCLER_SCRIPT_SHA256,
        "environment_name" => AUTOCYCLER_ENV_NAME,
        "environment_spec_expected_sha256" =>
            AUTOCYCLER_ENVIRONMENT_SPEC_SHA256,
        "environment_spec_sha256" => AUTOCYCLER_ENVIRONMENT_SPEC_SHA256,
        "inventory_schema" => AUTOCYCLER_PACKAGE_INVENTORY_SCHEMA,
        "package_inventory_sha256" => String(digest),
        "package_count" => Int(package_count),
        "packages" => Dict{String, Any}[
            Dict{String, Any}(
                "name" => record.name,
                "version" => record.version,
                "build" => record.build,
                "channel" => record.channel,
            ) for record in inventory
        ],
    )
end

function _with_autocycler_install_lock(
        action::Function,
        lock_path::AbstractString;
        stale_age::Real = AUTOCYCLER_INSTALL_LOCK_STALE_SECONDS,
        poll_interval::Real = 10,
        pidlock_runner::Function = FileWatching.Pidfile.mkpidlock,
)::Any
    stale_age > 0 || throw(ArgumentError("stale_age must be positive."))
    poll_interval > 0 || throw(ArgumentError("poll_interval must be positive."))
    normalized_lock_path = abspath(lock_path)
    mkpath(dirname(normalized_lock_path))
    return pidlock_runner(
        action,
        normalized_lock_path;
        stale_age,
        poll_interval,
    )
end

function _autocycler_environment_is_installed(
        conda_runner::AbstractString = _conda_runner();
        environment_prefix::AbstractString =
            _autocycler_environment_prefix(conda_runner),
)::Bool
    resolved_prefix =
        _canonical_autocycler_environment_prefix(environment_prefix)
    return isdir(joinpath(resolved_prefix, "conda-meta"))
end

function _create_autocycler_environment_from_yaml(
        environment_file::AbstractString,
        environment_name::AbstractString;
        force::Bool = false,
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::AbstractString =
            _autocycler_environment_prefix(conda_runner),
        command_runner::Function = run,
)::String
    String(environment_name) == AUTOCYCLER_ENV_NAME || throw(ArgumentError(
        "Autocycler environment name must be $(repr(AUTOCYCLER_ENV_NAME)).",
    ))
    force && throw(ArgumentError(
        "Autocycler environments are immutable and cannot be force-created.",
    ))
    normalized_environment_file = abspath(String(environment_file))
    isfile(normalized_environment_file) || throw(ArgumentError(
        "Autocycler environment file does not exist: " *
        normalized_environment_file,
    ))
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix =
        _canonical_autocycler_environment_prefix(environment_prefix)
    _autocycler_environment_is_installed(
        resolved_runner;
        environment_prefix = resolved_prefix,
    ) && return resolved_prefix
    command_runner(Cmd(String[
        resolved_runner,
        "env",
        "create",
        "-f",
        normalized_environment_file,
        "-p",
        resolved_prefix,
    ]))
    command_runner(Cmd(String[resolved_runner, "clean", "--all", "-y"]))
    return resolved_prefix
end

function _install_autocycler_locked(;
        force::Bool = false,
        downloader::Function = Downloads.download,
        paths::Tuple{String, String, String} = _autocycler_paths(),
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        environment_checker::Union{Nothing, Function} = nothing,
        environment_creator::Union{Nothing, Function} = nothing,
        package_inspector::Union{Nothing, Function} = nothing,
)::String
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix = environment_prefix === nothing ?
                      _autocycler_environment_prefix(resolved_runner) :
                      _canonical_autocycler_environment_prefix(
        environment_prefix,
    )
    effective_environment_checker = environment_checker === nothing ?
                                    () -> _autocycler_environment_is_installed(
        resolved_runner;
        environment_prefix = resolved_prefix,
    ) : environment_checker
    effective_environment_creator = environment_creator === nothing ?
                                    (
        environment_file::AbstractString,
        environment_name::AbstractString;
        force::Bool = false,
    ) -> _create_autocycler_environment_from_yaml(
        environment_file,
        environment_name;
        force,
        conda_runner = resolved_runner,
        environment_prefix = resolved_prefix,
    ) : environment_creator
    effective_package_inspector = package_inspector === nothing ?
                                  () -> _autocycler_environment_packages(;
        conda_runner = resolved_runner,
        environment_prefix = resolved_prefix,
    ) : package_inspector
    install_dir, script_path, env_file_path = paths
    mkpath(install_dir)
    verified_env_file_path =
        _require_verified_autocycler_environment_spec(env_file_path)

    environment_installed = Bool(effective_environment_checker())
    if force && environment_installed
        throw(
            ErrorException(
                "Refusing install_autocycler(force=true): immutable " *
                "spec-hash-addressed environment " *
                "$(repr(AUTOCYCLER_ENV_NAME)) already exists and must not " *
                "be recreated in place. Stop all workflows using this " *
                "environment and remove it manually before reinstalling.",
            ),
        )
    end
    if !environment_installed
        effective_environment_creator(
            verified_env_file_path,
            AUTOCYCLER_ENV_NAME;
            force = false,
        )
    end

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

    _ensure_autocycler_packages!(effective_package_inspector)

    @info "Autocycler installed successfully."
    return script_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install Autocycler and the short-read polishing tools in a spec-hash-addressed
conda environment.

The package-bundled `environment.yml` is authoritative because it extends the
upstream Autocycler environment with BWA, Polypolish, and Pypolca. The upstream
automation script is pinned to `AUTOCYCLER_SCRIPT_REVISION` and verified against
`AUTOCYCLER_SCRIPT_SHA256` before it can be executed.

The bundled environment specification is checksum-verified and a 16-character
prefix of its expected SHA-256 is part of `AUTOCYCLER_ENV_NAME`. Consequently,
Mycelia releases with different environment specifications cannot mutate the
environment used by an already-running assembly. An existing environment is
immutable: this function never removes or recreates it. `force=true` is accepted
for a first installation, where it refreshes the verified script, but fails closed
when the spec-hash-addressed environment already exists. A stale or incompatible
environment must not be repaired in place; stop all workflows using it before
manually removing and reinstalling it, or use a Mycelia release with a corrected
environment specification. Direct calls serialize the environment check and
first-use creation under the same per-environment lock held by assembly and
polishing workflows for their complete command lifecycles.

Compatibility is deliberately pinned to Autocycler 0.5.2. This wrapper targets
Autocycler's bacterial-isolate use case, where the alternative input assemblies
are expected to be mostly complete; installing it does not make the workflow a
general metagenome, eukaryotic-genome, or fragmentary-assembly consensus method.

# Keywords
- `force::Bool=false`: Refresh the pinned script during a first installation;
  error if the spec-hash-addressed environment already exists.
- `conda_runner::AbstractString`: Conda-compatible executable that owns the
  selected environment root.
- `environment_prefix::Union{Nothing,AbstractString}`: Exact environment prefix
  to create or inspect; pass the same value to `run_autocycler` or
  `run_autocycler_polished`.
"""
function install_autocycler(;
        force::Bool = false,
        downloader::Function = Downloads.download,
        paths::Tuple{String, String, String} = _autocycler_paths(),
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        environment_checker::Union{Nothing, Function} = nothing,
        environment_creator::Union{Nothing, Function} = nothing,
        package_inspector::Union{Nothing, Function} = nothing,
        lock_path::Union{Nothing, AbstractString} = nothing,
        lock_runner::Function = _with_autocycler_install_lock,
)::String
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix = environment_prefix === nothing ?
                      _autocycler_environment_prefix(resolved_runner) :
                      _canonical_autocycler_environment_prefix(
        environment_prefix,
    )
    resolved_lock_path = lock_path === nothing ?
                         _autocycler_install_lock_path_from_prefix(
        resolved_prefix,
    ) :
                         String(lock_path)
    return lock_runner(resolved_lock_path) do
        _install_autocycler_locked(;
            force,
            downloader,
            paths,
            conda_runner = resolved_runner,
            environment_prefix = resolved_prefix,
            environment_checker,
            environment_creator,
            package_inspector,
        )
    end
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

function _require_contained_regular_autocycler_artifact(
        path::AbstractString,
        workflow_root::AbstractString,
        label::AbstractString,
)::String
    normalized_root = normpath(abspath(workflow_root))
    isdir(normalized_root) || throw(
        ErrorException(
            "Autocycler workflow root is not a directory while validating " *
            "$(label): $(normalized_root).",
        ),
    )
    islink(normalized_root) && throw(
        ErrorException(
            "Autocycler workflow root became a symbolic link while validating " *
            "$(label): $(normalized_root).",
        ),
    )
    canonical_root = realpath(normalized_root)
    canonical_root == normalized_root || throw(
        ErrorException(
            "Autocycler workflow root is no longer canonical while validating " *
            "$(label): expected $(normalized_root), resolved $(canonical_root).",
        ),
    )

    normalized_path = normpath(abspath(path))
    lexical_relative_path = relpath(normalized_path, normalized_root)
    lexical_components = splitpath(lexical_relative_path)
    if isabspath(lexical_relative_path) || isempty(lexical_components) ||
       first(lexical_components) == ".."
        throw(
            ErrorException(
                "$(label) escapes the reserved Autocycler workflow root: " *
                "$(normalized_path) is not contained by $(normalized_root).",
            ),
        )
    end
    normalized_path = _require_nonempty_autocycler_file(path, label)
    islink(normalized_path) && throw(
        ErrorException(
            "$(label) must be a regular, non-symlink artifact: " *
            "$(normalized_path).",
        ),
    )
    canonical_path = realpath(normalized_path)
    canonical_path == normalized_path || throw(
        ErrorException(
            "$(label) resolves through a symbolic-link component: " *
            "$(normalized_path) resolves to $(canonical_path).",
        ),
    )
    canonical_relative_path = relpath(canonical_path, canonical_root)
    canonical_components = splitpath(canonical_relative_path)
    if isabspath(canonical_relative_path) || isempty(canonical_components) ||
       first(canonical_components) == ".."
        throw(
            ErrorException(
                "$(label) resolves outside the reserved Autocycler workflow " *
                "root: $(canonical_path) is not contained by " *
                "$(canonical_root).",
            ),
        )
    end
    return normalized_path
end

function _autocycler_artifact_snapshot(
        path::AbstractString,
        workflow_root::AbstractString,
        label::AbstractString,
)::NamedTuple
    normalized_path = _require_contained_regular_autocycler_artifact(
        path,
        workflow_root,
        label,
    )
    return (
        path = normalized_path,
        size_bytes = filesize(normalized_path),
        sha256 = _autocycler_sha256(normalized_path),
    )
end

function _autocycler_buffered_artifact_snapshot(
        path::AbstractString,
        workflow_root::AbstractString,
        label::AbstractString,
)::NamedTuple
    normalized_path = _require_contained_regular_autocycler_artifact(
        path,
        workflow_root,
        label,
    )
    bytes = open(normalized_path, "r") do input
        read(input)
    end
    isempty(bytes) && throw(ErrorException(
        "$(label) became empty while its immutable byte snapshot was read: " *
        "$(normalized_path)",
    ))
    snapshot = (
        path = normalized_path,
        size_bytes = length(bytes),
        sha256 = _autocycler_sha256(bytes),
    )
    return (; snapshot, bytes)
end

function _require_unchanged_autocycler_artifact(
        expected_snapshot::NamedTuple,
        path::AbstractString,
        workflow_root::AbstractString,
        label::AbstractString,
)::NamedTuple
    observed_snapshot = _autocycler_artifact_snapshot(
        path,
        workflow_root,
        label,
    )
    observed_snapshot == expected_snapshot || throw(
        ErrorException(
            "$(label) changed after its validated Autocycler snapshot: " *
            "expected size=$(expected_snapshot.size_bytes), " *
            "sha256=$(expected_snapshot.sha256); observed " *
            "size=$(observed_snapshot.size_bytes), " *
            "sha256=$(observed_snapshot.sha256). Refusing to return mutated " *
            "assembly artifacts.",
        ),
    )
    return observed_snapshot
end

function _autocycler_fasta_sequence_map_from_reader(
        reader::Any,
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    sequences = Dict{String, String}()
    try
        for (record_count, record) in enumerate(reader)
            record isa FASTX.FASTA.Record || throw(
                ErrorException(
                    "$(label) is not valid FASTA: $(path).",
                ),
            )
            identifier = String(FASTX.identifier(record))
            isempty(identifier) && throw(
                ErrorException(
                    "$(label) contains an empty FASTA identifier at record " *
                    "$(record_count): $(path).",
                ),
            )
            haskey(sequences, identifier) && throw(
                ErrorException(
                    "$(label) contains duplicate FASTA identifier " *
                    "$(repr(identifier)): $(path).",
                ),
            )
            sequence = FASTX.sequence(String, record)
            isempty(sequence) && throw(
                ErrorException(
                    "$(label) contains an empty FASTA sequence at record " *
                    "$(record_count): $(path).",
                ),
            )
            occursin(
                r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$",
                sequence,
            ) || throw(
                ErrorException(
                    "$(label) contains invalid DNA at FASTA record " *
                    "$(record_count): $(path).",
                ),
            )
            try
                BioSequences.LongDNA{4}(sequence)
            catch error
                error isa InterruptException && rethrow()
                throw(
                    ErrorException(
                        "$(label) contains invalid DNA at FASTA record " *
                        "$(record_count): $(path). Cause: " *
                        sprint(showerror, error),
                    ),
                )
            end
            sequences[identifier] = uppercase(sequence)
        end
    catch error
        error isa InterruptException && rethrow()
        if error isa ErrorException && startswith(error.msg, String(label))
            rethrow()
        end
        throw(
            ErrorException(
                "$(label) is not valid FASTA: $(path). Cause: " *
                sprint(showerror, error),
            ),
        )
    finally
        close(reader)
    end
    isempty(sequences) && throw(
        ErrorException("$(label) contains no FASTA records: $(path)."),
    )
    return sequences
end

function _autocycler_fasta_sequence_map(
        bytes::AbstractVector{UInt8},
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    reader = FASTX.FASTA.Reader(IOBuffer(bytes))
    return _autocycler_fasta_sequence_map_from_reader(reader, path, label)
end

function _require_matching_autocycler_contig_identifiers(
        source_identifiers::Set{String},
        target_identifiers::Set{String},
        source_label::AbstractString,
        target_label::AbstractString,
)::Nothing
    source_identifiers == target_identifiers && return nothing
    missing = sort!(collect(setdiff(source_identifiers, target_identifiers)))
    added = sort!(collect(setdiff(target_identifiers, source_identifiers)))
    throw(
        ErrorException(
            "Autocycler polishing changed contig identifiers between " *
            "$(source_label) and $(target_label): missing=$(repr(missing)), " *
            "added=$(repr(added)). Refusing to run or return a downstream " *
            "stage with dropped, added, or renamed contigs.",
        ),
    )
end

function _require_valid_autocycler_gfa(
        path::AbstractString,
        label::AbstractString,
)::String
    return _require_valid_metamdbg_gfa(path, label)
end

function _autocycler_gfa_sequence_map(
        bytes::AbstractVector{UInt8},
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    input = IOBuffer(bytes)
    _require_valid_metamdbg_gfa_input(input, label, path)
    seekstart(input)
    sequences = Dict{String, String}()
    for (line_number, line) in enumerate(eachline(input))
        startswith(line, "S\t") || continue
        fields = split(line, '\t'; keepempty = true)
        length(fields) >= 3 || throw(
            ErrorException(
                "$(label) has a malformed segment at line $(line_number): " *
                "$(path).",
            ),
        )
        identifier = fields[2]
        sequence = uppercase(fields[3])
        (isempty(sequence) || sequence == "*") && throw(
            ErrorException(
                "$(label) segment $(repr(identifier)) has no sequence, so it " *
                "cannot be verified against the companion FASTA.",
            ),
        )
        haskey(sequences, identifier) && throw(
            ErrorException(
                "$(label) contains duplicate segment identifier " *
                "$(repr(identifier)).",
            ),
        )
        sequences[identifier] = sequence
    end
    isempty(sequences) && throw(
        ErrorException("$(label) contains no sequence-bearing GFA segments."),
    )
    return sequences
end

function _require_matching_autocycler_companion_sequence_maps(
        fasta_sequences::AbstractDict,
        gfa_sequences::AbstractDict,
)::Nothing
    fasta_identifiers = Set(keys(fasta_sequences))
    gfa_identifiers = Set(keys(gfa_sequences))
    fasta_identifiers == gfa_identifiers || throw(
        ErrorException(
            "Autocycler consensus FASTA and GFA contain different contig " *
            "identifiers: FASTA-only=" *
            "$(repr(sort!(collect(setdiff(fasta_identifiers, gfa_identifiers))))), " *
            "GFA-only=" *
            "$(repr(sort!(collect(setdiff(gfa_identifiers, fasta_identifiers))))).",
        ),
    )
    mismatched = sort!(String[
        identifier for identifier in fasta_identifiers if
        fasta_sequences[identifier] != gfa_sequences[identifier]
    ])
    isempty(mismatched) || throw(
        ErrorException(
            "Autocycler consensus FASTA and GFA contain different sequences " *
            "for contigs $(repr(mismatched)).",
        ),
    )
    return nothing
end

function _require_matching_autocycler_companion_artifacts(
        assembly_bytes::AbstractVector{UInt8},
        assembly::AbstractString,
        graph_bytes::AbstractVector{UInt8},
        graph::AbstractString,
)::Nothing
    fasta_sequences = _autocycler_fasta_sequence_map(
        assembly_bytes,
        assembly,
        "Autocycler consensus FASTA",
    )
    gfa_sequences = _autocycler_gfa_sequence_map(
        graph_bytes,
        graph,
        "Autocycler consensus GFA",
    )
    return _require_matching_autocycler_companion_sequence_maps(
        fasta_sequences,
        gfa_sequences,
    )
end

function _require_expected_autocycler_artifact_snapshot(
        expected_snapshot::NamedTuple,
        buffered_snapshot::NamedTuple,
        label::AbstractString,
)::Nothing
    buffered_snapshot == expected_snapshot && return nothing
    throw(ErrorException(
        "$(label) changed before immutable semantic validation: expected " *
        "size=$(expected_snapshot.size_bytes), " *
        "sha256=$(expected_snapshot.sha256); observed " *
        "size=$(buffered_snapshot.size_bytes), " *
        "sha256=$(buffered_snapshot.sha256).",
    ))
end

function _autocycler_buffered_semantic_fasta(
        expected_snapshot::NamedTuple,
        path::AbstractString,
        workflow_root::AbstractString,
        label::AbstractString,
)::NamedTuple
    buffered = _autocycler_buffered_artifact_snapshot(
        path,
        workflow_root,
        label,
    )
    _require_expected_autocycler_artifact_snapshot(
        expected_snapshot,
        buffered.snapshot,
        label,
    )
    sequences = _autocycler_fasta_sequence_map(
        buffered.bytes,
        buffered.snapshot.path,
        label,
    )
    _require_unchanged_autocycler_artifact(
        expected_snapshot,
        buffered.snapshot.path,
        workflow_root,
        label,
    )
    return (;
        path = buffered.snapshot.path,
        snapshot = buffered.snapshot,
        sequences,
    )
end

function _autocycler_buffered_semantic_gfa(
        expected_snapshot::NamedTuple,
        path::AbstractString,
        workflow_root::AbstractString,
        label::AbstractString,
)::NamedTuple
    buffered = _autocycler_buffered_artifact_snapshot(
        path,
        workflow_root,
        label,
    )
    _require_expected_autocycler_artifact_snapshot(
        expected_snapshot,
        buffered.snapshot,
        label,
    )
    sequences = _autocycler_gfa_sequence_map(
        buffered.bytes,
        buffered.snapshot.path,
        label,
    )
    _require_unchanged_autocycler_artifact(
        expected_snapshot,
        buffered.snapshot.path,
        workflow_root,
        label,
    )
    return (;
        path = buffered.snapshot.path,
        snapshot = buffered.snapshot,
        sequences,
    )
end

function _autocycler_pair_identifier(identifier::AbstractString)::String
    first_token = first(split(String(identifier)))
    return replace(first_token, r"/[12]$" => "")
end

function _autocycler_identifier_pair_role(
        identifier::AbstractString,
)::Union{Nothing, Int}
    first_token = first(split(String(identifier)))
    role_match = match(r"/([12])$", first_token)
    return role_match === nothing ? nothing :
           parse(Int, something(only(role_match.captures)))
end

function _autocycler_casava_pair_role(
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

function _autocycler_pair_role(
        identifier::AbstractString,
        description::AbstractString = "",
)::Union{Nothing, Int}
    identifier_role = _autocycler_identifier_pair_role(identifier)
    casava_role = _autocycler_casava_pair_role(description)
    if identifier_role !== nothing && casava_role !== nothing &&
       identifier_role != casava_role
        throw(
            ArgumentError(
                "FASTQ identifier and CASAVA description contain conflicting " *
                "explicit mate roles: identifier=$(repr(String(identifier))), " *
                "description=$(repr(String(description))).",
            ),
        )
    end
    return identifier_role === nothing ? casava_role : identifier_role
end

function _validate_autocycler_fastq_reader(
        reader::Any,
        path::AbstractString,
        label::AbstractString,
)::Int
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
        caught isa ArgumentError && rethrow()
        throw(ArgumentError(
            "$(label) must be a FASTQ file: $(abspath(path)). Cause: " *
            sprint(showerror, caught),
        ))
    finally
        close(reader)
    end
    record_count > 0 || throw(ArgumentError("$(label) must be non-empty."))
    return record_count
end

function _validate_autocycler_fastq(
        path::AbstractString,
        label::AbstractString,
)::Int
    return _validate_autocycler_fastq_reader(
        Mycelia.open_fastx(path),
        path,
        label,
    )
end

function _validate_autocycler_paired_fastqs(
        short_reads_1::AbstractString,
        short_reads_2::AbstractString,
)::Int
    !Base.Filesystem.samefile(short_reads_1, short_reads_2) || throw(ArgumentError(
        "Autocycler paired short-read R1 and R2 must be distinct files.",
    ))
    reader_1 = Mycelia.open_fastx(short_reads_1)
    reader_2 = try
        Mycelia.open_fastx(short_reads_2)
    catch
        try
            close(reader_1)
        finally
            rethrow()
        end
    end
    return _validate_autocycler_paired_fastq_readers(
        reader_1,
        reader_2,
    )
end

function _validate_autocycler_paired_fastq_readers(
        reader_1::Any,
        reader_2::Any,
)::Int
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
            description_1 = String(FASTX.description(record_1))
            description_2 = String(FASTX.description(record_2))
            role_1 = _autocycler_pair_role(identifier_1, description_1)
            role_2 = _autocycler_pair_role(identifier_2, description_2)
            roles_valid = (role_1 === nothing && role_2 === nothing) ||
                          (role_1 == 1 && role_2 == 2)
            roles_valid || throw(ArgumentError(
                "Autocycler paired short reads have invalid explicit mate " *
                "roles at record $(pair_count): " *
                "R1=$(repr(identifier_1)), R2=$(repr(identifier_2)); " *
                "expected R1 role 1 then R2 role 2 from /1,/2 suffixes or " *
                "CASAVA descriptions.",
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
    catch caught
        caught isa InterruptException && rethrow()
        caught isa ArgumentError && rethrow()
        throw(ArgumentError(
            "Autocycler paired short-read inputs must be FASTQ files. " *
            "Cause: $(sprint(showerror, caught))",
        ))
    finally
        try
            close(reader_1)
        finally
            close(reader_2)
        end
    end
    pair_count > 0 || throw(
        ArgumentError("Autocycler paired short reads must be non-empty."),
    )
    return pair_count
end

function _validate_autocycler_input_sources(
        long_reads::AbstractString,
        short_reads_1::AbstractString,
        short_reads_2::AbstractString,
)::Nothing
    !Base.Filesystem.samefile(short_reads_1, short_reads_2) || throw(
        ArgumentError(
            "Autocycler paired short-read R1 and R2 must be distinct files.",
        ),
    )
    if Base.Filesystem.samefile(long_reads, short_reads_1) ||
       Base.Filesystem.samefile(long_reads, short_reads_2)
        throw(
            ArgumentError(
                "Autocycler long reads must be physically distinct from paired " *
                "short-read R1 and R2 files.",
            ),
        )
    end
    return nothing
end

function _effective_autocycler_assembly_threads(threads::Integer)::Int
    threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    return Int(min(threads, AUTOCYCLER_MAX_ASSEMBLY_THREADS))
end

function _autocycler_path_entry_exists(path::AbstractString)::Bool
    return ispath(path) || islink(path)
end

function _canonical_autocycler_output_path(out_dir::AbstractString)::String
    normalized_out_dir = abspath(out_dir)
    missing_components = String[]
    existing_ancestor = normalized_out_dir
    while !_autocycler_path_entry_exists(existing_ancestor)
        parent = dirname(existing_ancestor)
        parent == existing_ancestor && break
        pushfirst!(missing_components, basename(existing_ancestor))
        existing_ancestor = parent
    end
    isdir(existing_ancestor) || throw(
        ArgumentError(
            "Autocycler output directory has a non-directory existing " *
            "ancestor: $(existing_ancestor)",
        ),
    )
    canonical_ancestor = realpath(existing_ancestor)
    if isempty(missing_components)
        return canonical_ancestor
    end
    return joinpath(canonical_ancestor, missing_components...)
end

function _autocycler_output_adjacent_spool_parent(
        out_dir::AbstractString,
)::String
    parent = dirname(normpath(abspath(String(out_dir))))
    while !ispath(parent) && !islink(parent)
        next_parent = dirname(parent)
        next_parent == parent && break
        parent = next_parent
    end
    isdir(parent) && !islink(parent) || throw(ArgumentError(
        "Autocycler output-adjacent spool parent must resolve to an existing " *
        "non-symlink directory: $(parent).",
    ))
    return realpath(parent)
end

function _autocycler_output_lock_path(out_dir::AbstractString)::String
    canonical_out_dir = _canonical_autocycler_output_path(out_dir)
    return _autocycler_output_lock_path_from_canonical(canonical_out_dir)
end

function _autocycler_output_lock_path_from_canonical(
        canonical_out_dir::AbstractString,
)::String
    return _output_root_reservation_lock_path_from_canonical(canonical_out_dir)
end

function _with_autocycler_output_lock(
        action::Function,
        out_dir::AbstractString;
        stale_age::Real = AUTOCYCLER_OUTPUT_LOCK_STALE_SECONDS,
        poll_interval::Real = 10,
        pidlock_runner::Function = FileWatching.Pidfile.trymkpidlock,
)::Any
    stale_age > 0 || throw(ArgumentError("stale_age must be positive."))
    poll_interval > 0 || throw(ArgumentError("poll_interval must be positive."))
    canonical_out_dir = _canonical_autocycler_output_path(out_dir)
    lock_path = _autocycler_output_lock_path_from_canonical(canonical_out_dir)
    _require_exclusive_output_root_reservation(
        canonical_out_dir,
        lock_path;
        subject = "Autocycler output directory",
        stale_age,
    )
    mkpath(dirname(lock_path))
    locked_action = function ()
        _require_exclusive_output_root_reservation(
            canonical_out_dir,
            lock_path;
            subject = "Autocycler output directory",
            stale_age,
        )
        return action(canonical_out_dir)
    end
    if pidlock_runner !== FileWatching.Pidfile.trymkpidlock
        return pidlock_runner(
            locked_action,
            lock_path;
            stale_age,
            poll_interval,
        )
    end
    lock_handle = pidlock_runner(
        lock_path;
        stale_age,
        refresh = stale_age / 2,
    )
    lock_handle === false && throw(ArgumentError(
        "Autocycler output directory is already reserved by another " *
        "output-root workflow: $(lock_path)",
    ))
    try
        return locked_action()
    finally
        Base.close(lock_handle)
    end
end

function _validate_autocycler_parameters(
        threads::Integer,
        jobs::Integer,
        read_type::AbstractString,
)::String
    _effective_autocycler_assembly_threads(threads)
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
    isempty(strip(out_dir)) && throw(
        ArgumentError("Autocycler output directory must be a non-blank path."),
    )
    normalized_out_dir = abspath(out_dir)
    if islink(normalized_out_dir)
        throw(
            ArgumentError(
                "Autocycler output path must not be a symbolic link: " *
                "$(normalized_out_dir)",
            ),
        )
    end
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

    existing_ancestor = normalized_out_dir
    while !_autocycler_path_entry_exists(existing_ancestor)
        parent = dirname(existing_ancestor)
        parent == existing_ancestor && break
        existing_ancestor = parent
    end
    if _autocycler_path_entry_exists(existing_ancestor) &&
       !isdir(existing_ancestor)
        throw(
            ArgumentError(
                "Autocycler output directory has a non-directory existing " *
                "ancestor: $(existing_ancestor)",
            ),
        )
    end
    _canonical_autocycler_output_path(normalized_out_dir)
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
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::AbstractString =
            _autocycler_environment_prefix(conda_runner),
)::Cmd
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix =
        _canonical_autocycler_environment_prefix(environment_prefix)
    command_arguments = String[
        resolved_runner,
        "run",
        "--live-stream",
        "-p",
        resolved_prefix,
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
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::AbstractString =
            _autocycler_environment_prefix(conda_runner),
)::NamedTuple
    normalized_read_type = _validate_autocycler_parameters(
        threads,
        jobs,
        read_type,
    )
    requested_threads = Int(threads)
    autocycler_assembly_threads = _effective_autocycler_assembly_threads(threads)
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
            string(autocycler_assembly_threads),
            string(jobs),
            normalized_read_type,
        ],
        normalized_out_dir;
        conda_runner = conda_runner,
        environment_prefix = environment_prefix,
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
        requested_threads = requested_threads,
        autocycler_assembly_threads = autocycler_assembly_threads,
        jobs = Int(jobs),
        read_type = normalized_read_type,
    )
end

function _autocycler_polishing_command_plan(
        assembly::AbstractString,
        short_reads_1::AbstractString,
        short_reads_2::AbstractString,
        out_dir::AbstractString;
        threads::Integer = max(Sys.CPU_THREADS, 1),
        polypolish_careful::Bool = true,
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::AbstractString =
            _autocycler_environment_prefix(conda_runner),
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
            environment_prefix = environment_prefix,
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
            environment_prefix = environment_prefix,
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
            environment_prefix = environment_prefix,
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
            environment_prefix = environment_prefix,
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
            environment_prefix = environment_prefix,
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
            environment_prefix = environment_prefix,
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

function _require_planned_autocycler_path_containment(
        path::AbstractString,
        workflow_root::AbstractString,
        label::AbstractString,
)::String
    normalized_root = normpath(abspath(String(workflow_root)))
    isdir(normalized_root) && !islink(normalized_root) || throw(
        ErrorException(
            "Autocycler workflow root is not a regular directory while " *
            "checking $(label): $(normalized_root).",
        ),
    )
    realpath(normalized_root) == normalized_root || throw(
        ErrorException(
            "Autocycler workflow root resolves through a symbolic-link " *
            "component while checking $(label): $(normalized_root).",
        ),
    )
    normalized_path = normpath(abspath(String(path)))
    relative_path = relpath(normalized_path, normalized_root)
    components = splitpath(relative_path)
    if relative_path == "."
        return normalized_path
    end
    if isabspath(relative_path) || isempty(components) || first(components) == ".."
        throw(
            ErrorException(
                "$(label) escapes the reserved Autocycler workflow root: " *
                "$(normalized_path) is not contained by $(normalized_root).",
            ),
        )
    end
    cursor = normalized_root
    for component in components
        cursor = joinpath(cursor, component)
        islink(cursor) && throw(
            ErrorException(
                "$(label) resolves through a symbolic-link component: " *
                "$(cursor).",
            ),
        )
        if ispath(cursor) && realpath(cursor) != cursor
            throw(
                ErrorException(
                    "$(label) is not bound to its canonical path: " *
                    "$(cursor) resolves to $(realpath(cursor)).",
                ),
            )
        end
    end
    return normalized_path
end

function _autocycler_directory_identity(
        path::AbstractString,
        workflow_root::AbstractString,
        label::AbstractString,
)::NamedTuple
    normalized_path = _require_planned_autocycler_path_containment(
        path,
        workflow_root,
        label,
    )
    isdir(normalized_path) && !islink(normalized_path) || throw(
        ErrorException("$(label) is not a regular directory: $(normalized_path)."),
    )
    metadata = stat(normalized_path)
    return (;
        path = normalized_path,
        device = metadata.device,
        inode = metadata.inode,
        label = String(label),
    )
end

function _require_unchanged_autocycler_directory(
        identity::NamedTuple,
        workflow_root::AbstractString,
)::String
    normalized_path = _require_planned_autocycler_path_containment(
        identity.path,
        workflow_root,
        identity.label,
    )
    isdir(normalized_path) && !islink(normalized_path) || throw(
        ErrorException(
            "$(identity.label) is no longer a regular directory: " *
            "$(normalized_path).",
        ),
    )
    metadata = stat(normalized_path)
    if metadata.device != identity.device || metadata.inode != identity.inode
        throw(
            ErrorException(
                "$(identity.label) changed physical identity during the " *
                "Autocycler workflow: $(normalized_path).",
            ),
        )
    end
    return normalized_path
end

function _require_safe_autocycler_step_paths(
        step::NamedTuple,
        workflow_root::AbstractString,
        directory_identities::Tuple,
)::Nothing
    for identity in directory_identities
        _require_unchanged_autocycler_directory(identity, workflow_root)
    end
    command_directory = String(step.command.dir)
    isempty(command_directory) || begin
        normalized_directory = _require_planned_autocycler_path_containment(
            command_directory,
            workflow_root,
            "Autocycler workflow step $(step.name) working directory",
        )
        isdir(normalized_directory) && !islink(normalized_directory) || throw(
            ErrorException(
                "Autocycler workflow step $(step.name) working directory is " *
                "not a regular directory: $(normalized_directory).",
            ),
        )
    end
    candidate_paths = String[step.expected_outputs...]
    isnothing(step.stdout) || push!(candidate_paths, String(step.stdout))
    for output_path in unique(candidate_paths)
        _require_planned_autocycler_path_containment(
            output_path,
            workflow_root,
            "Autocycler workflow step $(step.name) output",
        )
    end
    return nothing
end

function _cleanup_autocycler_polishing_intermediates!(
        paths::AbstractVector{<:AbstractString};
        workflow_root::AbstractString,
        directory_identities::Tuple = (),
        remover::Function = path -> rm(path; force = true),
)::Vector{String}
    normalized_paths = String[
        normpath(abspath(String(path))) for path in paths
    ]
    cleanup_failures = Set{String}()
    for normalized_path in normalized_paths
        cleanup_error = nothing
        try
            for identity in directory_identities
                _require_unchanged_autocycler_directory(
                    identity,
                    workflow_root,
                )
            end
            _require_planned_autocycler_path_containment(
                normalized_path,
                workflow_root,
                "Autocycler polishing cleanup target",
            )
            remover(normalized_path)
        catch error
            error isa InterruptException && rethrow()
            cleanup_error = error
        end
        if cleanup_error !== nothing
            push!(cleanup_failures, normalized_path)
            cleanup_message =
                "Autocycler could not remove a polishing intermediate"
            @warn cleanup_message path = normalized_path cleanup_error
        end
    end

    # A later remover callback can recreate a path visited earlier. Rescan the
    # complete plan only after every removal attempt so no retained capability
    # can disappear from the returned inventory.
    retained = String[]
    for normalized_path in normalized_paths
        if ispath(normalized_path) || islink(normalized_path)
            push!(retained, normalized_path)
            if !(normalized_path in cleanup_failures)
                retained_message =
                    "Autocycler retained a polishing intermediate after cleanup"
                @warn retained_message path = normalized_path
            end
        end
    end
    return retained
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
        workflow_root::Union{Nothing, AbstractString} = nothing,
        directory_identities::Tuple = (),
        before_step::Function = step -> nothing,
        after_step::Function = step -> nothing,
)::Nothing
    for step in steps
        if workflow_root !== nothing
            _require_safe_autocycler_step_paths(
                step,
                workflow_root,
                directory_identities,
            )
        end
        before_step(step)
        if workflow_root !== nothing
            _require_safe_autocycler_step_paths(
                step,
                workflow_root,
                directory_identities,
            )
        end
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
                        "The spec-hash-addressed environment " *
                        "$(repr(AUTOCYCLER_ENV_NAME)) is immutable. If a " *
                        "command is missing, stop all workflows using this " *
                        "environment and remove it manually before reinstalling. " *
                        "Cause: $(message)",
                    ),
                )
            end
            rethrow()
        end

        for output_path in step.expected_outputs
            if workflow_root !== nothing
                if islink(output_path)
                    throw(
                        ErrorException(
                            "Autocycler workflow step $(step.name) created a " *
                            "symlinked artifact instead of a regular file: " *
                            "$(output_path)",
                        ),
                    )
                end
                if !isfile(output_path) || filesize(output_path) == 0
                    throw(
                        ErrorException(
                            "Autocycler workflow step $(step.name) did not " *
                            "create a nonempty artifact: $(output_path)",
                        ),
                    )
                end
                for identity in directory_identities
                    _require_unchanged_autocycler_directory(
                        identity,
                        workflow_root,
                    )
                end
                _require_planned_autocycler_path_containment(
                    output_path,
                    workflow_root,
                    "Autocycler workflow step $(step.name) output",
                )
                _require_contained_regular_autocycler_artifact(
                    output_path,
                    workflow_root,
                    "Autocycler workflow step $(step.name) output",
                )
                continue
            end
            if islink(output_path)
                throw(
                    ErrorException(
                        "Autocycler workflow step $(step.name) created a " *
                        "symlinked artifact instead of a regular file: " *
                        "$(output_path)",
                    ),
                )
            end
            if !isfile(output_path) || filesize(output_path) == 0
                throw(
                    ErrorException(
                        "Autocycler workflow step $(step.name) did not create " *
                        "a nonempty artifact: $(output_path)",
                    ),
                )
            end
        end
        after_step(step)
        if workflow_root !== nothing
            _require_safe_autocycler_step_paths(
                step,
                workflow_root,
                directory_identities,
            )
        end
    end
    return nothing
end

function _ensure_autocycler_installed_locked(
        package_inspector::Function,
        installer::Function,
        environment_checker::Function,
)::Dict{String, Any}
    _, script_path, env_file_path = _autocycler_paths()
    _require_verified_autocycler_environment_spec(env_file_path)
    environment_installed = Bool(environment_checker())
    if !environment_installed || !_autocycler_script_is_verified(script_path)
        @warn "Autocycler is not fully installed. Installing now..."
        installer()
        Bool(environment_checker()) || throw(
            ErrorException(
                "Autocycler environment was not available after installation.",
            ),
        )
    end
    if !_autocycler_script_is_verified(script_path)
        throw(
            ErrorException(
                "Autocycler script is unavailable or failed checksum " *
                "verification: $(script_path)",
            ),
        )
    end
    inventory = _ensure_autocycler_packages!(package_inspector)
    return _autocycler_toolchain_metadata(inventory)
end

function _ensure_autocycler_installed_locked_default(
        conda_runner::AbstractString = _conda_runner();
        environment_prefix::Union{Nothing, AbstractString} = nothing,
)::Dict{String, Any}
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix = environment_prefix === nothing ?
                      _autocycler_environment_prefix(resolved_runner) :
                      _canonical_autocycler_environment_prefix(
        environment_prefix,
    )
    package_inspector = () -> _autocycler_environment_packages(;
        conda_runner = resolved_runner,
        environment_prefix = resolved_prefix,
    )
    installer = () -> _install_autocycler_locked(;
        conda_runner = resolved_runner,
        environment_prefix = resolved_prefix,
    )
    environment_checker = () -> _autocycler_environment_is_installed(
        resolved_runner;
        environment_prefix = resolved_prefix,
    )
    return _ensure_autocycler_installed_locked(
        package_inspector,
        installer,
        environment_checker,
    )
end

function _inspect_autocycler_toolchain_locked(
        conda_runner::AbstractString = _conda_runner();
        environment_prefix::Union{Nothing, AbstractString} = nothing,
)::Dict{String, Any}
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix = environment_prefix === nothing ?
                      _autocycler_environment_prefix(resolved_runner) :
                      _canonical_autocycler_environment_prefix(
        environment_prefix,
    )
    Bool(_autocycler_environment_is_installed(
        resolved_runner;
        environment_prefix = resolved_prefix,
    )) || throw(
        ErrorException(
            "Autocycler environment disappeared while its shared lock was held.",
        ),
    )
    _, script_path, env_file_path = _autocycler_paths()
    _require_verified_autocycler_environment_spec(env_file_path)
    _autocycler_script_is_verified(script_path) || throw(
        ErrorException(
            "Autocycler script changed or disappeared while its shared lock " *
            "was held.",
        ),
    )
    package_inspector = () -> _autocycler_environment_packages(;
        conda_runner = resolved_runner,
        environment_prefix = resolved_prefix,
    )
    inventory = _ensure_autocycler_packages!(package_inspector)
    return _autocycler_toolchain_metadata(inventory)
end

function _autocycler_toolchain_snapshotter(
        dependency_checker::Function,
        conda_runner::AbstractString = _conda_runner();
        environment_prefix::Union{Nothing, AbstractString} = nothing,
)::Function
    return dependency_checker === _ensure_autocycler_installed_locked_default ?
           () -> _inspect_autocycler_toolchain_locked(
        conda_runner;
        environment_prefix,
    ) :
           dependency_checker
end

function _require_unchanged_autocycler_toolchain(
        expected_toolchain::AbstractDict,
        observed_toolchain::Any,
        boundary::AbstractString,
)::Dict{String, Any}
    normalized_expected =
        _require_autocycler_toolchain_provenance(expected_toolchain)
    normalized_observed =
        _require_autocycler_toolchain_provenance(observed_toolchain)
    expected_digest = normalized_expected["package_inventory_sha256"]
    observed_digest = normalized_observed["package_inventory_sha256"]
    normalized_observed == normalized_expected || throw(
        ErrorException(
            "Autocycler realized toolchain changed across $(boundary) while " *
            "the shared environment lock was held. Expected package inventory " *
            "$(expected_digest), observed $(observed_digest). Refusing " *
            "to return artifacts with ambiguous provenance.",
        ),
    )
    return normalized_observed
end

function _ensure_autocycler_installed(;
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        package_inspector::Union{Nothing, Function} = nothing,
        installer::Union{Nothing, Function} = nothing,
        environment_checker::Union{Nothing, Function} = nothing,
        lock_path::Union{Nothing, AbstractString} = nothing,
        lock_runner::Function = _with_autocycler_install_lock,
)::Dict{String, Any}
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix = environment_prefix === nothing ?
                      _autocycler_environment_prefix(resolved_runner) :
                      _canonical_autocycler_environment_prefix(
        environment_prefix,
    )
    effective_package_inspector = package_inspector === nothing ?
                                  () -> _autocycler_environment_packages(;
        conda_runner = resolved_runner,
        environment_prefix = resolved_prefix,
    ) : package_inspector
    effective_installer = installer === nothing ?
                          () -> _install_autocycler_locked(;
        conda_runner = resolved_runner,
        environment_prefix = resolved_prefix,
    ) : installer
    effective_environment_checker = environment_checker === nothing ?
                                    () -> _autocycler_environment_is_installed(
        resolved_runner;
        environment_prefix = resolved_prefix,
    ) : environment_checker
    resolved_lock_path = lock_path === nothing ?
                         _autocycler_install_lock_path_from_prefix(
        resolved_prefix,
    ) :
                         String(lock_path)
    return lock_runner(resolved_lock_path) do
        _ensure_autocycler_installed_locked(
            effective_package_inspector,
            effective_installer,
            effective_environment_checker,
        )
    end
end

function _run_autocycler_with_reserved_output(
        normalized_long_reads::AbstractString,
        reserved_out_dir::AbstractString;
        threads::Integer,
        jobs::Integer,
        read_type::AbstractString,
        toolchain::AbstractDict,
        runner::Function,
        conda_runner::AbstractString,
        environment_prefix::AbstractString,
        long_read_snapshot::NamedTuple,
        long_read_consumed_snapshot::NamedTuple,
        after_artifact_semantic_validation_hook::Function = artifacts -> nothing,
)::NamedTuple
    validated_out_dir = _validate_autocycler_output_dir(reserved_out_dir)
    normalized_toolchain =
        _require_autocycler_toolchain_provenance(toolchain)
    prepared_out_dir = _prepare_autocycler_output_dir(validated_out_dir)
    _, script_path, _ = _autocycler_paths()
    plan = _autocycler_command_plan(
        normalized_long_reads,
        prepared_out_dir;
        threads = threads,
        jobs = jobs,
        read_type,
        script_path = script_path,
        conda_runner,
        environment_prefix,
    )
    output_root_identity = _autocycler_directory_identity(
        prepared_out_dir,
        prepared_out_dir,
        "Autocycler workflow root",
    )
    function verify_long_read_before_step(step::NamedTuple)::Nothing
        if step.name == :autocycler
            _require_unchanged_autocycler_consumed_input(
                long_read_consumed_snapshot,
                "Long-read FASTQ",
            )
            _require_unchanged_autocycler_input(
                long_read_snapshot,
                "Long-read FASTQ",
            )
        end
        return nothing
    end
    output_snapshots = Ref{Any}(nothing)
    function verify_long_read_after_step(step::NamedTuple)::Nothing
        if step.name == :autocycler
            _require_unchanged_autocycler_consumed_input(
                long_read_consumed_snapshot,
                "Long-read FASTQ",
            )
            _require_unchanged_autocycler_input(
                long_read_snapshot,
                "Long-read FASTQ",
            )
            output_snapshots[] = (
                assembly = _autocycler_artifact_snapshot(
                    plan.assembly,
                    prepared_out_dir,
                    "Autocycler consensus FASTA",
                ),
                graph = _autocycler_artifact_snapshot(
                    plan.graph,
                    prepared_out_dir,
                    "Autocycler consensus GFA",
                ),
            )
        end
        return nothing
    end

    @info "Starting long-read-only Autocycler pipeline..."
    _execute_autocycler_steps(
        plan.steps;
        runner,
        workflow_root = prepared_out_dir,
        directory_identities = (output_root_identity,),
        before_step = verify_long_read_before_step,
        after_step = verify_long_read_after_step,
    )
    autocycler_output_snapshots = something(output_snapshots[])
    buffered_assembly = _autocycler_buffered_artifact_snapshot(
        plan.assembly,
        prepared_out_dir,
        "Autocycler consensus FASTA",
    )
    _require_expected_autocycler_artifact_snapshot(
        autocycler_output_snapshots.assembly,
        buffered_assembly.snapshot,
        "Autocycler consensus FASTA",
    )
    buffered_graph = _autocycler_buffered_artifact_snapshot(
        plan.graph,
        prepared_out_dir,
        "Autocycler consensus GFA",
    )
    _require_expected_autocycler_artifact_snapshot(
        autocycler_output_snapshots.graph,
        buffered_graph.snapshot,
        "Autocycler consensus GFA",
    )
    assembly = buffered_assembly.snapshot.path
    graph = buffered_graph.snapshot.path
    _require_matching_autocycler_companion_artifacts(
        buffered_assembly.bytes,
        assembly,
        buffered_graph.bytes,
        graph,
    )
    after_artifact_semantic_validation_hook((; assembly, graph))
    _require_unchanged_autocycler_artifact(
        autocycler_output_snapshots.assembly,
        assembly,
        prepared_out_dir,
        "Autocycler consensus FASTA",
    )
    _require_unchanged_autocycler_artifact(
        autocycler_output_snapshots.graph,
        graph,
        prepared_out_dir,
        "Autocycler consensus GFA",
    )
    @info "Autocycler pipeline complete" out_dir = prepared_out_dir

    return (
        outdir = prepared_out_dir,
        assembly = assembly,
        graph = graph,
        toolchain = normalized_toolchain,
        output_root_identity,
        assembly_snapshot = autocycler_output_snapshots.assembly,
        graph_snapshot = autocycler_output_snapshots.graph,
        requested_threads = plan.requested_threads,
        autocycler_assembly_threads = plan.autocycler_assembly_threads,
        jobs = plan.jobs,
        read_type = plan.read_type,
    )
end

function _run_autocycler(
        long_reads::AbstractString,
        out_dir::AbstractString;
        threads::Integer = max(Sys.CPU_THREADS, 1),
        jobs::Integer = 1,
        read_type::AbstractString = "ont_r10",
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        dependency_checker::Function =
            _ensure_autocycler_installed_locked_default,
        runner::Function = _default_autocycler_step_runner,
        after_input_semantic_validation_hook::Function = snapshots -> nothing,
        after_artifact_semantic_validation_hook::Function = artifacts -> nothing,
        output_lock_runner::Function = _with_autocycler_output_lock,
        environment_lock_path::Union{Nothing, AbstractString} = nothing,
        environment_lock_runner::Function = _with_autocycler_install_lock,
        input_spool_parent::Union{Nothing, AbstractString} = nothing,
        input_spool_byte_ceiling::Integer =
            AUTOCYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
        input_spool_available_bytes_reader::Function =
            parent -> diskstat(parent).available,
        after_input_spool_copy_hook::Function = binding -> nothing,
)::NamedTuple
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix = environment_prefix === nothing ?
                      _autocycler_environment_prefix(resolved_runner) :
                      _canonical_autocycler_environment_prefix(
        environment_prefix,
    )
    resolved_environment_lock_path = environment_lock_path === nothing ?
                                     _autocycler_install_lock_path_from_prefix(
        resolved_prefix,
    ) : String(environment_lock_path)
    normalized_read_type = _validate_autocycler_parameters(
        threads,
        jobs,
        read_type,
    )
    normalized_out_dir = _validate_autocycler_output_dir(out_dir)
    return output_lock_runner(normalized_out_dir) do reserved_out_dir
        environment_lock_runner(resolved_environment_lock_path) do
            resolved_spool_parent = input_spool_parent === nothing ?
                                    _autocycler_output_adjacent_spool_parent(
                reserved_out_dir,
            ) :
                                    input_spool_parent
            return _with_autocycler_spooled_input_snapshots(
                (long_reads,),
                ("Long-read FASTQ",);
                spool_parent = resolved_spool_parent,
                byte_ceiling = input_spool_byte_ceiling,
                available_bytes_reader = input_spool_available_bytes_reader,
                after_copy_hook = after_input_spool_copy_hook,
            ) do input_bindings
                long_read_binding = only(input_bindings)
                _validate_autocycler_fastq(
                    long_read_binding.semantic_path,
                    "Long-read input",
                )
                long_read_snapshot = long_read_binding.snapshot
                input_snapshots = (; long_reads = long_read_snapshot)
                after_input_semantic_validation_hook(input_snapshots)
                _require_unchanged_autocycler_input(
                    long_read_snapshot,
                    "Long-read FASTQ",
                )
                _require_unchanged_autocycler_consumed_input(
                    long_read_binding.consumed_snapshot,
                    "Long-read FASTQ",
                )
            toolchain_snapshotter = _autocycler_toolchain_snapshotter(
                dependency_checker,
                resolved_runner;
                environment_prefix = resolved_prefix,
            )
            dependency_result =
                dependency_checker === _ensure_autocycler_installed_locked_default ?
                _ensure_autocycler_installed_locked_default(
                    resolved_runner;
                    environment_prefix = resolved_prefix,
                ) :
                dependency_checker()
            initial_toolchain = _require_autocycler_toolchain_provenance(
                dependency_result,
            )
            result = _run_autocycler_with_reserved_output(
                long_read_binding.semantic_path,
                reserved_out_dir;
                threads,
                jobs,
                read_type = normalized_read_type,
                toolchain = initial_toolchain,
                runner,
                conda_runner = resolved_runner,
                environment_prefix = resolved_prefix,
                long_read_snapshot,
                long_read_consumed_snapshot =
                    long_read_binding.consumed_snapshot,
                after_artifact_semantic_validation_hook,
            )
            _require_unchanged_autocycler_toolchain(
                initial_toolchain,
                toolchain_snapshotter(),
                "long-read assembly",
            )
            _require_unchanged_autocycler_directory(
                result.output_root_identity,
                result.outdir,
            )
            final_assembly_binding = _autocycler_buffered_semantic_fasta(
                result.assembly_snapshot,
                result.assembly,
                result.outdir,
                "Autocycler consensus FASTA",
            )
            final_graph_binding = _autocycler_buffered_semantic_gfa(
                result.graph_snapshot,
                result.graph,
                result.outdir,
                "Autocycler consensus GFA",
            )
            _require_matching_autocycler_companion_sequence_maps(
                final_assembly_binding.sequences,
                final_graph_binding.sequences,
            )
            _require_unchanged_autocycler_artifact(
                result.assembly_snapshot,
                result.assembly,
                result.outdir,
                "Autocycler consensus FASTA",
            )
            _require_unchanged_autocycler_artifact(
                result.graph_snapshot,
                result.graph,
                result.outdir,
                "Autocycler consensus GFA",
            )
            _require_unchanged_autocycler_input(
                long_read_snapshot,
                "Long-read FASTQ",
            )
            _require_unchanged_autocycler_consumed_input(
                long_read_binding.consumed_snapshot,
                "Long-read FASTQ",
            )
            return (
                outdir = result.outdir,
                assembly = result.assembly,
                graph = result.graph,
                toolchain = result.toolchain,
                output_root_identity = result.output_root_identity,
                requested_threads = result.requested_threads,
                autocycler_assembly_threads =
                    result.autocycler_assembly_threads,
                jobs = result.jobs,
                read_type = result.read_type,
                artifact_snapshots = (;
                    assembly = result.assembly_snapshot,
                    graph = result.graph_snapshot,
                ),
                input_snapshots,
            )
            end
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the upstream long-read-only Autocycler pipeline.

Autocycler consumes one long-read FASTQ and writes its fixed `autocycler_out`
directory relative to the isolated `out_dir` working directory. Its authoritative
consensus FASTA and companion GFA are both validated and returned without
rewriting either artifact. Their unique contig identifiers and sequences must
agree exactly. The output directory must be absent or empty and is
reserved by an adjacent interprocess lock for the full workflow lifecycle.
Returned artifacts must remain regular non-symlink files contained by that
canonical reserved directory. Its canonical path, device, and inode are bound
across execution, final toolchain inspection, and return validation.
The shared Autocycler environment/install lock is also held from the initial
dependency snapshot through artifact validation. The normalized full Conda
inventory (name, version, build, and channel for every package) is checked before
and after execution; any drift fails rather than returning ambiguous provenance.
The reserved lifecycle first copies exactly the initially sized input into a
read-only stable snapshot. Spooling rejects source growth, shrinkage, or physical
replacement, enforces a cumulative byte ceiling, checks available space before
copying, and removes partial snapshots on failure.

The verified script and environment are compatibility-pinned to Autocycler 0.5.2.
The upstream consensus model is intended for bacterial isolates whose alternative
long-read assemblies are mostly complete. This wrapper does not widen that boundary
to metagenomes, eukaryotic genomes, or highly fragmentary alternative assemblies.

# Keywords
- `long_reads::AbstractString`: Nonempty long-read FASTQ.
- `out_dir::AbstractString`: Empty isolated working/output directory.
- `threads::Integer`: Threads per input assembler job.
- `jobs::Integer`: Number of simultaneous assembler jobs.
- `read_type::AbstractString`: One of `ont_r9`, `ont_r10`, `pacbio_clr`, or
  `pacbio_hifi`.
- `conda_runner::AbstractString`: Conda-compatible executable that owns the
  selected environment root.
- `environment_prefix::Union{Nothing,AbstractString}`: Exact environment prefix
  created by `install_autocycler`; `nothing` derives the spec-hash-addressed
  prefix from `conda_runner`.
- `input_spool_parent::Union{Nothing,AbstractString}`: Existing scratch
  directory for stable inputs; `nothing` uses the output-adjacent directory.
- `input_spool_byte_ceiling::Integer`: Maximum cumulative stable-input bytes.

# Returns
A named tuple with `outdir`, `assembly`, `graph`, and exact realized `toolchain`
provenance, including the deterministic full Conda inventory and its SHA-256.
`input_snapshots.long_reads` records the normalized and canonical input path,
byte size, and SHA-256 value that were verified immediately before consumption
and again after final artifact validation.
`artifact_snapshots` records the exact assembly/GFA byte sizes and SHA-256
values used for semantic validation and final live-path comparison.
`output_root_identity` records the bound canonical path, device, and inode used
for lifecycle swap detection.
`requested_threads` records the caller request and
`autocycler_assembly_threads` records the explicit effective cap used by the
alternative-assembly stage.
"""
function run_autocycler(;
        long_reads::AbstractString,
        out_dir::AbstractString,
        threads::Integer = max(Sys.CPU_THREADS, 1),
        jobs::Integer = 1,
        read_type::AbstractString = "ont_r10",
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        input_spool_parent::Union{Nothing, AbstractString} = nothing,
        input_spool_byte_ceiling::Integer =
            AUTOCYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
)::NamedTuple
    return _run_autocycler(
        long_reads,
        out_dir;
        threads = threads,
        jobs = jobs,
        read_type = read_type,
        conda_runner = conda_runner,
        environment_prefix = environment_prefix,
        input_spool_parent = input_spool_parent,
        input_spool_byte_ceiling = input_spool_byte_ceiling,
    )
end

function _run_autocycler_polished_with_reserved_output(
        normalized_long_reads::AbstractString,
        normalized_short_reads_1::AbstractString,
        normalized_short_reads_2::AbstractString,
        normalized_out_dir::AbstractString;
        threads::Integer,
        jobs::Integer,
        read_type::AbstractString,
        polypolish_careful::Bool,
        keep_intermediates::Bool,
        toolchain::AbstractDict,
        dependency_checker::Function,
        runner::Function,
        intermediate_remover::Function,
        conda_runner::AbstractString,
        environment_prefix::AbstractString,
        input_snapshots::NamedTuple,
        input_consumed_snapshots::NamedTuple,
        after_artifact_semantic_validation_hook::Function = artifacts -> nothing,
)::NamedTuple
    autocycler_result = _run_autocycler_with_reserved_output(
        normalized_long_reads,
        normalized_out_dir;
        threads,
        jobs,
        read_type,
        toolchain,
        runner,
        conda_runner,
        environment_prefix,
        long_read_snapshot = input_snapshots.long_reads,
        long_read_consumed_snapshot = input_consumed_snapshots.long_reads,
        after_artifact_semantic_validation_hook,
    )
    normalized_out_dir = autocycler_result.outdir
    _require_unchanged_autocycler_toolchain(
        autocycler_result.toolchain,
        dependency_checker(),
        "long-read assembly within the polished workflow",
    )
    _require_unchanged_autocycler_artifact(
        autocycler_result.assembly_snapshot,
        autocycler_result.assembly,
        autocycler_result.outdir,
        "Autocycler consensus FASTA",
    )
    _require_unchanged_autocycler_artifact(
        autocycler_result.graph_snapshot,
        autocycler_result.graph,
        autocycler_result.outdir,
        "Autocycler consensus GFA",
    )
    bound_autocycler_assembly = _autocycler_buffered_semantic_fasta(
        autocycler_result.assembly_snapshot,
        autocycler_result.assembly,
        autocycler_result.outdir,
        "Autocycler consensus FASTA",
    )
    autocycler_identifiers = Set(keys(bound_autocycler_assembly.sequences))
    autocycler_assembly_snapshot = autocycler_result.assembly_snapshot
    autocycler_graph_snapshot = autocycler_result.graph_snapshot

    polishing_plan = _autocycler_polishing_command_plan(
        autocycler_result.assembly,
        normalized_short_reads_1,
        normalized_short_reads_2,
        normalized_out_dir;
        threads = threads,
        polypolish_careful = polypolish_careful,
        conda_runner,
        environment_prefix,
    )
    _require_planned_autocycler_path_containment(
        polishing_plan.polishing_dir,
        normalized_out_dir,
        "Autocycler polishing directory",
    )
    if islink(polishing_plan.polishing_dir) ||
       (ispath(polishing_plan.polishing_dir) &&
        !isdir(polishing_plan.polishing_dir))
        throw(
            ArgumentError(
                "Autocycler polishing directory must be a regular " *
                "non-symlink directory: $(polishing_plan.polishing_dir)",
            ),
        )
    end
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
    output_root_identity = autocycler_result.output_root_identity
    _require_unchanged_autocycler_directory(
        output_root_identity,
        normalized_out_dir,
    )
    polishing_directory_identity = _autocycler_directory_identity(
        polishing_plan.polishing_dir,
        normalized_out_dir,
        "Autocycler polishing directory",
    )
    directory_identities = (
        output_root_identity,
        polishing_directory_identity,
    )
    produced_artifact_paths = Dict{Symbol, Vector{String}}(
        step.name => String[step.expected_outputs...] for
        step in polishing_plan.steps
    )
    produced_artifact_snapshots = Dict{Symbol, Vector{NamedTuple}}()
    function polishing_artifact_label(
            producer::Symbol,
            path::AbstractString,
    )::String
        if producer == :bwa_index
            return "BWA index artifact $(basename(path))"
        elseif producer == :bwa_mem_1
            return "BWA R1 alignment"
        elseif producer == :bwa_mem_2
            return "BWA R2 alignment"
        elseif producer == :polypolish_filter
            read_role = endswith(path, "filtered_1.sam") ? "R1" : "R2"
            return "Polypolish filtered $(read_role) alignment"
        elseif producer == :polypolish
            return "Polypolish Autocycler assembly"
        elseif producer == :pypolca && endswith(path, ".report")
            return "Pypolca report"
        elseif producer == :pypolca
            return "Pypolca-polished Autocycler assembly"
        end
        return "Autocycler $(producer) artifact $(basename(path))"
    end
    function captured_polishing_artifact_snapshot(
            producer::Symbol,
            path::AbstractString,
    )::NamedTuple
        snapshots = get(produced_artifact_snapshots, producer, nothing)
        snapshots === nothing && error(
            "Autocycler polishing step $(producer) has no bound output " *
            "snapshots.",
        )
        paths = produced_artifact_paths[producer]
        path_index = findfirst(==(String(path)), paths)
        path_index === nothing && error(
            "Autocycler polishing step $(producer) did not plan output " *
            "$(path).",
        )
        return snapshots[path_index]
    end
    function verify_produced_artifacts(
            producer::Symbol,
            consumer::Symbol,
    )::Nothing
        paths = produced_artifact_paths[producer]
        snapshots = get(produced_artifact_snapshots, producer, nothing)
        snapshots === nothing && error(
            "Autocycler polishing step $(consumer) cannot consume $(producer) " *
            "outputs before they are bound.",
        )
        length(paths) == length(snapshots) || error(
            "Autocycler polishing step $(producer) output snapshot count " *
            "changed before $(consumer).",
        )
        for (path, snapshot) in zip(paths, snapshots)
            _require_unchanged_autocycler_artifact(
                snapshot,
                path,
                normalized_out_dir,
                polishing_artifact_label(producer, path),
            )
        end
        return nothing
    end
    function verify_polishing_intermediate_artifacts(
            paths::AbstractVector{<:AbstractString},
    )::Nothing
        for path in paths
            normalized_path = String(path)
            producer = nothing
            for (candidate_producer, candidate_paths) in produced_artifact_paths
                if normalized_path in candidate_paths
                    producer = candidate_producer
                    break
                end
            end
            producer isa Symbol || error(
                "Autocycler polishing intermediate has no producing step: " *
                "$(normalized_path).",
            )
            snapshot = captured_polishing_artifact_snapshot(
                producer,
                normalized_path,
            )
            _require_unchanged_autocycler_artifact(
                snapshot,
                normalized_path,
                normalized_out_dir,
                polishing_artifact_label(producer, normalized_path),
            )
        end
        return nothing
    end
    function verify_polishing_step_inputs(step::NamedTuple)::Nothing
        if step.name in (
            :bwa_index,
            :bwa_mem_1,
            :bwa_mem_2,
            :polypolish,
        )
            _require_unchanged_autocycler_artifact(
                autocycler_assembly_snapshot,
                autocycler_result.assembly,
                normalized_out_dir,
                "Autocycler consensus FASTA",
            )
        end
        if step.name in (:bwa_mem_1, :bwa_mem_2)
            verify_produced_artifacts(:bwa_index, step.name)
        elseif step.name == :polypolish_filter
            verify_produced_artifacts(:bwa_mem_1, step.name)
            verify_produced_artifacts(:bwa_mem_2, step.name)
        elseif step.name == :polypolish
            verify_produced_artifacts(:polypolish_filter, step.name)
        elseif step.name == :pypolca
            verify_produced_artifacts(:polypolish, step.name)
        end
        if step.name in (:bwa_mem_1, :pypolca)
            _require_unchanged_autocycler_consumed_input(
                input_consumed_snapshots.short_reads_1,
                "Paired short-read R1 FASTQ",
            )
            _require_unchanged_autocycler_input(
                input_snapshots.short_reads_1,
                "Paired short-read R1 FASTQ",
            )
        end
        if step.name in (:bwa_mem_2, :pypolca)
            _require_unchanged_autocycler_consumed_input(
                input_consumed_snapshots.short_reads_2,
                "Paired short-read R2 FASTQ",
            )
            _require_unchanged_autocycler_input(
                input_snapshots.short_reads_2,
                "Paired short-read R2 FASTQ",
            )
        end
        return nothing
    end
    function verify_and_bind_polishing_step_outputs(
            step::NamedTuple,
    )::Nothing
        verify_polishing_step_inputs(step)
        produced_artifact_snapshots[step.name] = NamedTuple[
            _autocycler_artifact_snapshot(
                path,
                normalized_out_dir,
                polishing_artifact_label(step.name, path),
            ) for path in produced_artifact_paths[step.name]
        ]
        return nothing
    end

    @info "Starting paired-short polishing with Polypolish and Pypolca..."
    autocycler_assembly, graph, polypolish_assembly, final_assembly,
        pypolca_report, final_assembly_snapshot, pypolca_report_snapshot,
        polypolish_snapshot = try
        _execute_autocycler_steps(
            polishing_plan.steps[1:(end - 1)];
            runner,
            workflow_root = normalized_out_dir,
            directory_identities,
            before_step = verify_polishing_step_inputs,
            after_step = verify_and_bind_polishing_step_outputs,
        )
        polypolish_snapshot = captured_polishing_artifact_snapshot(
            :polypolish,
            polishing_plan.polypolish_assembly,
        )
        _require_unchanged_autocycler_artifact(
            polypolish_snapshot,
            polishing_plan.polypolish_assembly,
            normalized_out_dir,
            "Polypolish Autocycler assembly",
        )
        bound_polypolish_assembly = _autocycler_buffered_semantic_fasta(
            polypolish_snapshot,
            polishing_plan.polypolish_assembly,
            normalized_out_dir,
            "Polypolish Autocycler assembly",
        )
        validated_polypolish_assembly = bound_polypolish_assembly.path
        _require_unchanged_autocycler_artifact(
            autocycler_assembly_snapshot,
            autocycler_result.assembly,
            normalized_out_dir,
            "Autocycler consensus FASTA",
        )
        _require_unchanged_autocycler_artifact(
            autocycler_graph_snapshot,
            autocycler_result.graph,
            normalized_out_dir,
            "Autocycler consensus GFA",
        )
        polypolish_identifiers = Set(keys(bound_polypolish_assembly.sequences))
        _require_matching_autocycler_contig_identifiers(
            autocycler_identifiers,
            polypolish_identifiers,
            "Autocycler consensus",
            "Polypolish",
        )
        _require_unchanged_autocycler_artifact(
            polypolish_snapshot,
            validated_polypolish_assembly,
            normalized_out_dir,
            "Polypolish Autocycler assembly",
        )
        _execute_autocycler_steps(
            (last(polishing_plan.steps),);
            runner,
            workflow_root = normalized_out_dir,
            directory_identities,
            before_step = verify_polishing_step_inputs,
            after_step = verify_and_bind_polishing_step_outputs,
        )
        final_assembly_snapshot = captured_polishing_artifact_snapshot(
            :pypolca,
            polishing_plan.assembly,
        )
        pypolca_report_snapshot = captured_polishing_artifact_snapshot(
            :pypolca,
            polishing_plan.pypolca_report,
        )
        bound_final_assembly = _autocycler_buffered_semantic_fasta(
            final_assembly_snapshot,
            polishing_plan.assembly,
            normalized_out_dir,
            "Pypolca-polished Autocycler assembly",
        )
        validated_final_assembly = bound_final_assembly.path

        bound_final_autocycler_assembly =
            _autocycler_buffered_semantic_fasta(
                autocycler_assembly_snapshot,
                autocycler_result.assembly,
                normalized_out_dir,
                "Autocycler consensus FASTA",
            )
        validated_autocycler_assembly = bound_final_autocycler_assembly.path
        bound_final_graph = _autocycler_buffered_semantic_gfa(
            autocycler_graph_snapshot,
            autocycler_result.graph,
            normalized_out_dir,
            "Autocycler consensus GFA",
        )
        validated_graph = bound_final_graph.path
        bound_final_polypolish_assembly =
            _autocycler_buffered_semantic_fasta(
                polypolish_snapshot,
                validated_polypolish_assembly,
                normalized_out_dir,
                "Polypolish Autocycler assembly",
            )
        validated_polypolish_assembly =
            bound_final_polypolish_assembly.path
        validated_pypolca_report =
            _require_contained_regular_autocycler_artifact(
                polishing_plan.pypolca_report,
                normalized_out_dir,
                "Pypolca report",
            )
        _require_unchanged_autocycler_artifact(
            pypolca_report_snapshot,
            validated_pypolca_report,
            normalized_out_dir,
            "Pypolca report",
        )

        _require_matching_autocycler_companion_sequence_maps(
            bound_final_autocycler_assembly.sequences,
            bound_final_graph.sequences,
        )

        final_autocycler_identifiers = Set(keys(
            bound_final_autocycler_assembly.sequences,
        ))
        final_polypolish_identifiers = Set(keys(
            bound_final_polypolish_assembly.sequences,
        ))
        final_identifiers = Set(keys(bound_final_assembly.sequences))
        _require_matching_autocycler_contig_identifiers(
            final_autocycler_identifiers,
            final_polypolish_identifiers,
            "Autocycler consensus",
            "Polypolish",
        )
        _require_matching_autocycler_contig_identifiers(
            final_polypolish_identifiers,
            final_identifiers,
            "Polypolish",
            "Pypolca",
        )
        verify_polishing_intermediate_artifacts(
            polishing_plan.intermediate_files,
        )
        (
            validated_autocycler_assembly,
            validated_graph,
            validated_polypolish_assembly,
            validated_final_assembly,
            validated_pypolca_report,
            final_assembly_snapshot,
            pypolca_report_snapshot,
            polypolish_snapshot,
        )
    catch
        if !keep_intermediates
            _cleanup_autocycler_polishing_intermediates!(
                polishing_plan.intermediate_files;
                workflow_root = normalized_out_dir,
                directory_identities,
                remover = intermediate_remover,
            )
        end
        rethrow()
    end
    retained_intermediates = if keep_intermediates
        polishing_plan.intermediate_files
    else
        _cleanup_autocycler_polishing_intermediates!(
            polishing_plan.intermediate_files;
            workflow_root = normalized_out_dir,
            directory_identities,
            remover = intermediate_remover,
        )
    end
    verify_polishing_intermediate_artifacts(retained_intermediates)

    bound_return_autocycler_assembly = _autocycler_buffered_semantic_fasta(
        autocycler_assembly_snapshot,
        autocycler_assembly,
        normalized_out_dir,
        "Autocycler consensus FASTA",
    )
    autocycler_assembly = bound_return_autocycler_assembly.path
    bound_return_graph = _autocycler_buffered_semantic_gfa(
        autocycler_graph_snapshot,
        graph,
        normalized_out_dir,
        "Autocycler consensus GFA",
    )
    graph = bound_return_graph.path
    _require_matching_autocycler_companion_sequence_maps(
        bound_return_autocycler_assembly.sequences,
        bound_return_graph.sequences,
    )
    bound_return_polypolish_assembly = _autocycler_buffered_semantic_fasta(
        polypolish_snapshot,
        polypolish_assembly,
        normalized_out_dir,
        "Polypolish Autocycler assembly",
    )
    polypolish_assembly = bound_return_polypolish_assembly.path
    bound_return_final_assembly = _autocycler_buffered_semantic_fasta(
        final_assembly_snapshot,
        final_assembly,
        normalized_out_dir,
        "Pypolca-polished Autocycler assembly",
    )
    final_assembly = bound_return_final_assembly.path
    pypolca_report = _require_contained_regular_autocycler_artifact(
        pypolca_report,
        normalized_out_dir,
        "Pypolca report",
    )
    _require_unchanged_autocycler_artifact(
        pypolca_report_snapshot,
        pypolca_report,
        normalized_out_dir,
        "Pypolca report",
    )
    final_autocycler_identifiers =
        Set(keys(bound_return_autocycler_assembly.sequences))
    final_polypolish_identifiers =
        Set(keys(bound_return_polypolish_assembly.sequences))
    final_identifiers = Set(keys(bound_return_final_assembly.sequences))
    _require_matching_autocycler_contig_identifiers(
        final_autocycler_identifiers,
        final_polypolish_identifiers,
        "Autocycler consensus",
        "Polypolish",
    )
    _require_matching_autocycler_contig_identifiers(
        final_polypolish_identifiers,
        final_identifiers,
        "Polypolish",
        "Pypolca",
    )
    @info "Autocycler paired-short polishing complete" assembly = final_assembly

    retained_intermediate_snapshots = NamedTuple[]
    for retained_path in retained_intermediates
        producer = only(Symbol[
            candidate_producer for (candidate_producer, candidate_paths) in
            produced_artifact_paths if retained_path in candidate_paths
        ])
        push!(retained_intermediate_snapshots, (;
            path = retained_path,
            snapshot = captured_polishing_artifact_snapshot(
                producer,
                retained_path,
            ),
            label = polishing_artifact_label(producer, retained_path),
        ))
    end
    lifecycle_artifact_snapshots = (;
        autocycler_assembly = autocycler_assembly_snapshot,
        graph = autocycler_graph_snapshot,
        polypolish_assembly = polypolish_snapshot,
        assembly = final_assembly_snapshot,
        pypolca_report = pypolca_report_snapshot,
        intermediates = retained_intermediate_snapshots,
        planned_intermediates = copy(polishing_plan.intermediate_files),
    )
    return (
        outdir = normalized_out_dir,
        assembly = final_assembly,
        graph = graph,
        autocycler_assembly = autocycler_assembly,
        polypolish_assembly = polypolish_assembly,
        pypolca_report = pypolca_report,
        intermediates = retained_intermediates,
        toolchain = autocycler_result.toolchain,
        output_root_identity,
        requested_threads = autocycler_result.requested_threads,
        autocycler_assembly_threads =
            autocycler_result.autocycler_assembly_threads,
        jobs = autocycler_result.jobs,
        read_type = autocycler_result.read_type,
        polishing_threads = Int(threads),
        lifecycle_artifact_snapshots,
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
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        dependency_checker::Function =
            _ensure_autocycler_installed_locked_default,
        runner::Function = _default_autocycler_step_runner,
        after_input_semantic_validation_hook::Function = snapshots -> nothing,
        after_artifact_semantic_validation_hook::Function = artifacts -> nothing,
        intermediate_remover::Function = path -> rm(path; force = true),
        output_lock_runner::Function = _with_autocycler_output_lock,
        environment_lock_path::Union{Nothing, AbstractString} = nothing,
        environment_lock_runner::Function = _with_autocycler_install_lock,
        input_spool_parent::Union{Nothing, AbstractString} = nothing,
        input_spool_byte_ceiling::Integer =
            AUTOCYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
        input_spool_available_bytes_reader::Function =
            parent -> diskstat(parent).available,
        after_input_spool_copy_hook::Function = binding -> nothing,
)::NamedTuple
    resolved_runner = _canonical_autocycler_conda_runner(conda_runner)
    resolved_prefix = environment_prefix === nothing ?
                      _autocycler_environment_prefix(resolved_runner) :
                      _canonical_autocycler_environment_prefix(
        environment_prefix,
    )
    resolved_environment_lock_path = environment_lock_path === nothing ?
                                     _autocycler_install_lock_path_from_prefix(
        resolved_prefix,
    ) : String(environment_lock_path)
    normalized_read_type = _validate_autocycler_parameters(
        threads,
        jobs,
        read_type,
    )
    normalized_out_dir = _validate_autocycler_output_dir(out_dir)
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
    _validate_autocycler_input_sources(
        normalized_long_reads,
        normalized_short_reads_1,
        normalized_short_reads_2,
    )
    return output_lock_runner(normalized_out_dir) do reserved_out_dir
        environment_lock_runner(resolved_environment_lock_path) do
            resolved_spool_parent = input_spool_parent === nothing ?
                                    _autocycler_output_adjacent_spool_parent(
                reserved_out_dir,
            ) :
                                    input_spool_parent
            return _with_autocycler_spooled_input_snapshots(
                (
                    normalized_long_reads,
                    normalized_short_reads_1,
                    normalized_short_reads_2,
                ),
                (
                    "Long-read FASTQ",
                    "Paired short-read R1 FASTQ",
                    "Paired short-read R2 FASTQ",
                );
                spool_parent = resolved_spool_parent,
                byte_ceiling = input_spool_byte_ceiling,
                available_bytes_reader = input_spool_available_bytes_reader,
                after_copy_hook = after_input_spool_copy_hook,
            ) do input_bindings
                long_read_binding, short_read_1_binding,
                    short_read_2_binding = input_bindings
                _validate_autocycler_fastq(
                    long_read_binding.semantic_path,
                    "Long-read input",
                )
                _validate_autocycler_paired_fastqs(
                    short_read_1_binding.semantic_path,
                    short_read_2_binding.semantic_path,
                )
                input_snapshots = (
                    long_reads = long_read_binding.snapshot,
                    short_reads_1 = short_read_1_binding.snapshot,
                    short_reads_2 = short_read_2_binding.snapshot,
                )
                input_consumed_snapshots = (
                    long_reads = long_read_binding.consumed_snapshot,
                    short_reads_1 = short_read_1_binding.consumed_snapshot,
                    short_reads_2 = short_read_2_binding.consumed_snapshot,
                )
                after_input_semantic_validation_hook(input_snapshots)
                for (snapshot, consumed, label) in (
                        (
                            input_snapshots.long_reads,
                            input_consumed_snapshots.long_reads,
                            "Long-read FASTQ",
                        ),
                        (
                            input_snapshots.short_reads_1,
                            input_consumed_snapshots.short_reads_1,
                            "Paired short-read R1 FASTQ",
                        ),
                        (
                            input_snapshots.short_reads_2,
                            input_consumed_snapshots.short_reads_2,
                            "Paired short-read R2 FASTQ",
                        ),
                )
                    _require_unchanged_autocycler_input(snapshot, label)
                    _require_unchanged_autocycler_consumed_input(
                        consumed,
                        label,
                    )
                end
            toolchain_snapshotter = _autocycler_toolchain_snapshotter(
                dependency_checker,
                resolved_runner;
                environment_prefix = resolved_prefix,
            )
            dependency_result =
                dependency_checker === _ensure_autocycler_installed_locked_default ?
                _ensure_autocycler_installed_locked_default(
                    resolved_runner;
                    environment_prefix = resolved_prefix,
                ) :
                dependency_checker()
            initial_toolchain = _require_autocycler_toolchain_provenance(
                dependency_result,
            )
            result = _run_autocycler_polished_with_reserved_output(
                long_read_binding.semantic_path,
                short_read_1_binding.semantic_path,
                short_read_2_binding.semantic_path,
                reserved_out_dir;
                threads,
                jobs,
                read_type = normalized_read_type,
                polypolish_careful,
                keep_intermediates,
                toolchain = initial_toolchain,
                dependency_checker = toolchain_snapshotter,
                runner,
                intermediate_remover,
                conda_runner = resolved_runner,
                environment_prefix = resolved_prefix,
                input_snapshots,
                input_consumed_snapshots,
                after_artifact_semantic_validation_hook,
            )
            _require_unchanged_autocycler_toolchain(
                initial_toolchain,
                toolchain_snapshotter(),
                "the complete Autocycler and paired-short polishing lifecycle",
            )
            _require_unchanged_autocycler_directory(
                result.output_root_identity,
                result.outdir,
            )
            _require_unchanged_autocycler_input(
                input_snapshots.long_reads,
                "Long-read FASTQ",
            )
            _require_unchanged_autocycler_input(
                input_snapshots.short_reads_1,
                "Paired short-read R1 FASTQ",
            )
            _require_unchanged_autocycler_input(
                input_snapshots.short_reads_2,
                "Paired short-read R2 FASTQ",
            )
            for (consumed, label) in (
                    (
                        input_consumed_snapshots.long_reads,
                        "Long-read FASTQ",
                    ),
                    (
                        input_consumed_snapshots.short_reads_1,
                        "Paired short-read R1 FASTQ",
                    ),
                    (
                        input_consumed_snapshots.short_reads_2,
                        "Paired short-read R2 FASTQ",
                    ),
            )
                _require_unchanged_autocycler_consumed_input(
                    consumed,
                    label,
                )
            end
            lifecycle_snapshots = result.lifecycle_artifact_snapshots
            for (snapshot_name, artifact_path, label) in (
                    (
                        :autocycler_assembly,
                        result.autocycler_assembly,
                        "Autocycler consensus FASTA",
                    ),
                    (:graph, result.graph, "Autocycler consensus GFA"),
                    (
                        :polypolish_assembly,
                        result.polypolish_assembly,
                        "Polypolish Autocycler assembly",
                    ),
                    (
                        :assembly,
                        result.assembly,
                        "Pypolca-polished Autocycler assembly",
                    ),
                    (:pypolca_report, result.pypolca_report, "Pypolca report"),
            )
                _require_unchanged_autocycler_artifact(
                    getproperty(lifecycle_snapshots, snapshot_name),
                    artifact_path,
                    result.outdir,
                    label,
                )
            end
            retained_paths = Set(
                retained.path for retained in lifecycle_snapshots.intermediates
            )
            for planned_path in lifecycle_snapshots.planned_intermediates
                normalized_path =
                    normpath(abspath(String(planned_path)))
                _require_planned_autocycler_path_containment(
                    normalized_path,
                    result.outdir,
                    "Autocycler planned polishing intermediate",
                )
                if (ispath(normalized_path) || islink(normalized_path)) &&
                   !(normalized_path in retained_paths)
                    throw(ErrorException(
                        "Autocycler polishing intermediate reappeared after " *
                        "its completed cleanup: $(normalized_path). Refusing " *
                        "to return an incomplete cleanup inventory.",
                    ))
                end
            end
            for retained in lifecycle_snapshots.intermediates
                _require_unchanged_autocycler_artifact(
                    retained.snapshot,
                    retained.path,
                    result.outdir,
                    retained.label,
                )
            end
            return (;
                outdir = result.outdir,
                assembly = result.assembly,
                graph = result.graph,
                autocycler_assembly = result.autocycler_assembly,
                polypolish_assembly = result.polypolish_assembly,
                pypolca_report = result.pypolca_report,
                intermediates = result.intermediates,
                toolchain = result.toolchain,
                output_root_identity = result.output_root_identity,
                requested_threads = result.requested_threads,
                autocycler_assembly_threads =
                    result.autocycler_assembly_threads,
                polishing_threads = result.polishing_threads,
                jobs = result.jobs,
                read_type = result.read_type,
                lifecycle_artifact_snapshots =
                    result.lifecycle_artifact_snapshots,
                input_snapshots,
            )
            end
        end
    end
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
An adjacent interprocess lock reserves the output directory continuously across
both long-read assembly and paired-short polishing. The shared environment/install
lock is held across that same lifecycle. Full normalized Conda inventory snapshots
must match after long-read assembly and again after polishing. Polypolish and
Pypolca must also preserve the exact unique Autocycler contig-identifier set;
renamed, dropped, or added contigs fail before the next stage or before return.
The raw graph, raw consensus, and Polypolish FASTA are content-snapshotted and
must remain byte-identical through every downstream command. Immediately before
return, every reported artifact is revalidated as a contained regular non-symlink
file; all FASTAs/GFA are parsed again, all three identifier sets are compared
again, and the Pypolca report must still be nonempty.
The long-read and both paired-short inputs are bound to normalized and canonical
paths, byte sizes, and SHA-256 values before execution. Each input is rechecked
immediately before every command that consumes it and after final artifact
validation.
The workflow root and polishing directory are bound by canonical path, device,
and inode. Every step output is checked before and after execution, and cleanup
refuses any target whose ancestor or bound directory changed identity.
All three inputs are copied once under the held output and environment
reservations into read-only stable snapshots. Commands consume those snapshots;
source growth, shrinkage, or replacement, insufficient free space, and a
configured cumulative-byte-ceiling violation fail with partial-spool cleanup.

This workflow remains compatibility-pinned to Autocycler 0.5.2 and retains the
same bacterial-isolate/mostly-complete-alternative-assemblies applicability
boundary. Paired-short polishing improves the resulting consensus; it does not
turn Autocycler into a general metagenome, eukaryotic-genome, or fragmentary-input
ensemble method.

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
- `conda_runner::AbstractString`: Conda-compatible executable that owns the
  selected environment root.
- `environment_prefix::Union{Nothing,AbstractString}`: Exact environment prefix
  created by `install_autocycler`; `nothing` derives the spec-hash-addressed
  prefix from `conda_runner`.
- `input_spool_parent::Union{Nothing,AbstractString}`: Existing scratch
  directory for stable inputs; `nothing` uses the output-adjacent directory.
- `input_spool_byte_ceiling::Integer`: Maximum cumulative stable-input bytes.

# Returns
A named tuple with final `assembly`, raw `graph`, `autocycler_assembly`,
`polypolish_assembly`, `pypolca_report`, `outdir`, and exact realized `toolchain`
provenance, including the deterministic full Conda inventory and its SHA-256.
`input_snapshots` records the normalized and canonical paths, byte sizes, and
SHA-256 values bound for the long-read and both paired-short inputs.
`lifecycle_artifact_snapshots` binds the raw consensus/GFA, both polished
assemblies, Pypolca report, and every retained intermediate to their exact byte
sizes and SHA-256 values for the high-level hybrid return boundary.
`output_root_identity` records the bound canonical path, device, and inode used
for lifecycle swap detection.
`intermediates` lists files retained by an explicit
`keep_intermediates=true` request or because best-effort cleanup could not remove
them. `requested_threads` and `autocycler_assembly_threads` distinguish the caller
request from Autocycler's effective assembly cap; `polishing_threads` records the
uncapped caller request used by the threaded BWA and Pypolca polishing steps.
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
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        input_spool_parent::Union{Nothing, AbstractString} = nothing,
        input_spool_byte_ceiling::Integer =
            AUTOCYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
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
        conda_runner = conda_runner,
        environment_prefix = environment_prefix,
        input_spool_parent = input_spool_parent,
        input_spool_byte_ceiling = input_spool_byte_ceiling,
    )
end
