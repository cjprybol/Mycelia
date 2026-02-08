"""
$(DocStringExtensions.TYPEDSIGNATURES)

Canonical description of a portable compute job that can be rendered to
multiple backends (SLURM, Docker, Cloud Build).
"""
struct JobSpec
    job_name::String
    cmd::String
    site::Symbol
    time_limit::String
    workdir::Union{Nothing, String}
    env::Dict{String, String}
    modules::Vector{String}
    nodes::Int
    ntasks::Union{Nothing, Int}
    ntasks_per_node::Union{Nothing, Int}
    cpus_per_task::Int
    mem_gb::Union{Nothing, Float64}
    mem_per_cpu_mb::Union{Nothing, Int}
    gpus_per_node::Union{Nothing, Int}
    gpus_per_task::Union{Nothing, Int}
    gpu_bind::Union{Nothing, String}
    partition::Union{Nothing, String}
    qos::Union{Nothing, String}
    constraint::Union{Nothing, String}
    account::Union{Nothing, String}
    output_path::Union{Nothing, String}
    error_path::Union{Nothing, String}
    mail_user::Union{Nothing, String}
    mail_type::Union{Nothing, Vector{String}}
    requeue::Bool
    container_image::Union{Nothing, String}
    container_args::Vector{String}
    container_mounts::Vector{String}
    container_workdir::Union{Nothing, String}
    container_engine::Union{Nothing, Symbol}
    template_name::Union{Nothing, String}
end

Base.@kwdef struct JobSpecValidation
    errors::Vector{String} = String[]
    warnings::Vector{String} = String[]
end

Base.@kwdef struct SubmitResult
    ok::Bool = false
    dry_run::Bool = true
    site::Symbol = :local
    backend::Symbol = :none
    artifact_path::Union{Nothing, String} = nothing
    artifact_text::Union{Nothing, String} = nothing
    submit_command::Union{Nothing, String} = nothing
    job_id::Union{Nothing, String} = nothing
    stdout::Union{Nothing, String} = nothing
    warnings::Vector{String} = String[]
    errors::Vector{String} = String[]
end

const VALID_JOB_SITES = Set([:nersc, :lawrencium, :scg, :local, :cloudbuild])
const SITE_ALIASES = Dict(
    "nersc" => :nersc,
    "lawrencium" => :lawrencium,
    "scg" => :scg,
    "stanford_scg" => :scg,
    "local" => :local,
    "cloudbuild" => :cloudbuild
)
const VALID_CONTAINER_ENGINES = Set([:docker, :apptainer, :shifter, :podman_hpc])
const CONTAINER_ENGINE_ALIASES = Dict(
    "docker" => :docker,
    "apptainer" => :apptainer,
    "singularity" => :apptainer,
    "shifter" => :shifter,
    "podman-hpc" => :podman_hpc,
    "podman_hpc" => :podman_hpc
)
const TEMPLATE_ROOT = normpath(joinpath(@__DIR__, "..", "templates"))
const DEFAULT_SLURM_LOGDIR = joinpath(homedir(), "workspace", "slurmlogs")
const DEFAULT_SLURM_SCRIPTDIR = joinpath(homedir(), "workspace", "slurm")
const TEMPLATE_TOKEN_REGEX = r"\{\{[A-Z0-9_]+\}\}"

function JobSpec(;
        job_name::AbstractString,
        cmd::AbstractString,
        site::Union{Symbol, AbstractString},
        time_limit::AbstractString,
        workdir = nothing,
        env = nothing,
        modules = nothing,
        nodes::Integer = 1,
        ntasks = nothing,
        ntasks_per_node = nothing,
        cpus_per_task::Integer = 1,
        mem_gb = nothing,
        mem_per_cpu_mb = nothing,
        gpus_per_node = nothing,
        gpus_per_task = nothing,
        gpu_bind = nothing,
        partition = nothing,
        qos = nothing,
        constraint = nothing,
        account = nothing,
        output_path = nothing,
        error_path = nothing,
        mail_user = nothing,
        mail_type = nothing,
        requeue::Bool = false,
        container_image = nothing,
        container_args = nothing,
        container_mounts = nothing,
        container_workdir = nothing,
        container_engine = nothing,
        template_name = nothing
)
    return JobSpec(
        String(job_name),
        String(cmd),
        _normalize_site(site),
        String(time_limit),
        _normalize_optional_string(workdir),
        _normalize_env_dict(env),
        _normalize_string_vector(modules),
        Int(nodes),
        _normalize_optional_int(ntasks),
        _normalize_optional_int(ntasks_per_node),
        Int(cpus_per_task),
        _normalize_optional_float(mem_gb),
        _normalize_optional_int(mem_per_cpu_mb),
        _normalize_optional_int(gpus_per_node),
        _normalize_optional_int(gpus_per_task),
        _normalize_optional_string(gpu_bind),
        _normalize_optional_string(partition),
        _normalize_optional_string(qos),
        _normalize_optional_string(constraint),
        _normalize_optional_string(account),
        _normalize_optional_string(output_path),
        _normalize_optional_string(error_path),
        _normalize_optional_string(mail_user),
        _normalize_mail_type(mail_type),
        requeue,
        _normalize_optional_string(container_image),
        _normalize_string_vector(container_args),
        _normalize_string_vector(container_mounts),
        _normalize_optional_string(container_workdir),
        _normalize_container_engine(container_engine),
        _normalize_optional_string(template_name)
    )
end

function _normalize_site(site::Symbol)::Symbol
    if site in VALID_JOB_SITES
        return site
    end
    site_str = lowercase(String(site))
    return _normalize_site(site_str)
end

function _normalize_site(site::AbstractString)::Symbol
    candidate = lowercase(strip(String(site)))
    if haskey(SITE_ALIASES, candidate)
        return SITE_ALIASES[candidate]
    end
    error("Unsupported site '$site'. Allowed sites: nersc, lawrencium, scg, local, cloudbuild")
end

function _normalize_container_engine(engine::Nothing)
    return nothing
end

function _normalize_container_engine(engine::Symbol)
    if engine in VALID_CONTAINER_ENGINES
        return engine
    end
    return _normalize_container_engine(String(engine))
end

function _normalize_container_engine(engine::AbstractString)
    candidate = lowercase(strip(String(engine)))
    if haskey(CONTAINER_ENGINE_ALIASES, candidate)
        return CONTAINER_ENGINE_ALIASES[candidate]
    end
    error("Unsupported container engine '$engine'. Allowed values: docker, apptainer, shifter, podman-hpc")
end

function _normalize_optional_string(value)
    if value === nothing
        return nothing
    end
    normalized = strip(String(value))
    return isempty(normalized) ? nothing : normalized
end

function _normalize_optional_int(value)
    if value === nothing
        return nothing
    end
    return Int(value)
end

function _normalize_optional_float(value)
    if value === nothing
        return nothing
    end
    return Float64(value)
end

function _normalize_string_vector(value)::Vector{String}
    if value === nothing
        return String[]
    elseif value isa AbstractVector
        return [String(v) for v in value if !isempty(strip(String(v)))]
    else
        return [String(value)]
    end
end

function _normalize_env_dict(value)::Dict{String, String}
    env = Dict{String, String}()
    if value === nothing
        return env
    elseif value isa AbstractDict
        for (k, v) in value
            key = strip(String(k))
            if isempty(key)
                continue
            end
            env[key] = String(v)
        end
        return env
    end
    error("env must be an associative dictionary of KEY => VALUE pairs")
end

function _normalize_mail_type(mail_type)
    if mail_type === nothing
        return nothing
    elseif mail_type isa AbstractString
        normalized = [strip(x) for x in split(String(mail_type), ",") if !isempty(strip(x))]
        return isempty(normalized) ? nothing : normalized
    elseif mail_type isa AbstractVector
        normalized = [strip(String(x)) for x in mail_type if !isempty(strip(String(x)))]
        return isempty(normalized) ? nothing : normalized
    end
    error("mail_type must be a String, Vector{String}, or nothing")
end

function _container_engine_to_string(engine::Union{Nothing, Symbol})
    if engine === nothing
        return nothing
    elseif engine == :podman_hpc
        return "podman-hpc"
    else
        return String(engine)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a `JobSpec` into a plain `Dict` for serialization.
"""
function job_spec_to_dict(job::JobSpec)::Dict{String, Any}
    return Dict(
        "job_name" => job.job_name,
        "cmd" => job.cmd,
        "site" => String(job.site),
        "time_limit" => job.time_limit,
        "workdir" => job.workdir,
        "env" => Dict(job.env),
        "modules" => collect(job.modules),
        "nodes" => job.nodes,
        "ntasks" => job.ntasks,
        "ntasks_per_node" => job.ntasks_per_node,
        "cpus_per_task" => job.cpus_per_task,
        "mem_gb" => job.mem_gb,
        "mem_per_cpu_mb" => job.mem_per_cpu_mb,
        "gpus_per_node" => job.gpus_per_node,
        "gpus_per_task" => job.gpus_per_task,
        "gpu_bind" => job.gpu_bind,
        "partition" => job.partition,
        "qos" => job.qos,
        "constraint" => job.constraint,
        "account" => job.account,
        "output_path" => job.output_path,
        "error_path" => job.error_path,
        "mail_user" => job.mail_user,
        "mail_type" => job.mail_type,
        "requeue" => job.requeue,
        "container_image" => job.container_image,
        "container_args" => collect(job.container_args),
        "container_mounts" => collect(job.container_mounts),
        "container_workdir" => job.container_workdir,
        "container_engine" => _container_engine_to_string(job.container_engine),
        "template_name" => job.template_name
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a `JobSpec` from a dictionary.
"""
function job_spec_from_dict(data::AbstractDict)::JobSpec
    return JobSpec(
        job_name = _dict_required(data, "job_name"),
        cmd = _dict_required(data, "cmd"),
        site = _dict_required(data, "site"),
        time_limit = _dict_required(data, "time_limit"),
        workdir = _dict_get(data, "workdir"),
        env = _dict_get(data, "env"),
        modules = _dict_get(data, "modules"),
        nodes = _dict_get(data, "nodes", 1),
        ntasks = _dict_get(data, "ntasks"),
        ntasks_per_node = _dict_get(data, "ntasks_per_node"),
        cpus_per_task = _dict_get(data, "cpus_per_task", 1),
        mem_gb = _dict_get(data, "mem_gb"),
        mem_per_cpu_mb = _dict_get(data, "mem_per_cpu_mb"),
        gpus_per_node = _dict_get(data, "gpus_per_node"),
        gpus_per_task = _dict_get(data, "gpus_per_task"),
        gpu_bind = _dict_get(data, "gpu_bind"),
        partition = _dict_get(data, "partition"),
        qos = _dict_get(data, "qos"),
        constraint = _dict_get(data, "constraint"),
        account = _dict_get(data, "account"),
        output_path = _dict_get(data, "output_path"),
        error_path = _dict_get(data, "error_path"),
        mail_user = _dict_get(data, "mail_user"),
        mail_type = _dict_get(data, "mail_type"),
        requeue = Bool(_dict_get(data, "requeue", false)),
        container_image = _dict_get(data, "container_image"),
        container_args = _dict_get(data, "container_args"),
        container_mounts = _dict_get(data, "container_mounts"),
        container_workdir = _dict_get(data, "container_workdir"),
        container_engine = _dict_get(data, "container_engine"),
        template_name = _dict_get(data, "template_name")
    )
end

function _dict_get(data::AbstractDict, key::AbstractString, default = nothing)
    if haskey(data, key)
        return data[key]
    end
    symbol_key = Symbol(key)
    if haskey(data, symbol_key)
        return data[symbol_key]
    end
    return default
end

function _dict_required(data::AbstractDict, key::AbstractString)
    value = _dict_get(data, key, nothing)
    if value === nothing
        error("Missing required key '$key' in JobSpec payload")
    end
    return value
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Serialize a `JobSpec` to JSON.
"""
function job_spec_to_json(job::JobSpec; pretty::Bool = true)::String
    payload = job_spec_to_dict(job)
    if pretty
        return JSON.json(payload, 4)
    end
    return JSON.json(payload)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a `JobSpec` to a JSON file.
"""
function write_job_spec_json(job::JobSpec, path::AbstractString; pretty::Bool = true)::String
    mkpath(dirname(String(path)))
    write(String(path), job_spec_to_json(job; pretty = pretty))
    return String(path)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Deserialize a `JobSpec` from JSON text or a JSON file path.
"""
function job_spec_from_json(input::AbstractString)::JobSpec
    input_string = String(input)
    stripped_input = strip(input_string)

    looks_like_json = startswith(stripped_input, "{") || startswith(stripped_input, "[")
    path_exists = if looks_like_json
        false
    else
        try
            isfile(input_string)
        catch
            false
        end
    end

    raw = if path_exists
        read(input_string, String)
    else
        input_string
    end

    payload = JSON.parse(raw)
    if !(payload isa AbstractDict)
        error("JobSpec JSON must decode to an object")
    end
    return job_spec_from_dict(payload)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

YAML serialization is not enabled by default.
Use JSON helpers unless YAML.jl is added to this project.
"""
function job_spec_to_yaml(job::JobSpec)
    _ = job
    error("YAML support is not enabled. Use job_spec_to_json / job_spec_from_json.")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

YAML deserialization is not enabled by default.
Use JSON helpers unless YAML.jl is added to this project.
"""
function job_spec_from_yaml(input::AbstractString)
    _ = input
    error("YAML support is not enabled. Use job_spec_to_json / job_spec_from_json.")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a SLURM walltime string and return total seconds.
Accepts `HH:MM:SS` and `D-HH:MM:SS`.
"""
function parse_time_limit_seconds(time_limit::AbstractString)::Union{Nothing, Int}
    value = strip(String(time_limit))
    if isempty(value)
        return nothing
    end

    day_match = match(r"^(\d+)-(\d{1,2}):(\d{2}):(\d{2})$", value)
    if day_match !== nothing
        days = parse(Int, day_match.captures[1])
        hours = parse(Int, day_match.captures[2])
        minutes = parse(Int, day_match.captures[3])
        seconds = parse(Int, day_match.captures[4])
        if hours > 23 || minutes > 59 || seconds > 59
            return nothing
        end
        return (((days * 24) + hours) * 60 + minutes) * 60 + seconds
    end

    hour_match = match(r"^(\d+):(\d{2}):(\d{2})$", value)
    if hour_match !== nothing
        hours = parse(Int, hour_match.captures[1])
        minutes = parse(Int, hour_match.captures[2])
        seconds = parse(Int, hour_match.captures[3])
        if minutes > 59 || seconds > 59
            return nothing
        end
        return (hours * 60 + minutes) * 60 + seconds
    end

    return nothing
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Validate a `JobSpec`, returning errors and warnings.
Errors should block submission unless an override is explicitly requested.
"""
function validate(job::JobSpec; check_lawrencium_associations::Bool = true)::JobSpecValidation
    normalized_job = normalize_job_spec(job)
    report = JobSpecValidation()

    if isempty(strip(normalized_job.job_name))
        push!(report.errors, "job_name is required")
    end
    if isempty(strip(normalized_job.cmd))
        push!(report.errors, "cmd is required")
    end

    parsed_time = parse_time_limit_seconds(normalized_job.time_limit)
    if parsed_time === nothing
        push!(report.errors, "time_limit must match HH:MM:SS or D-HH:MM:SS")
    end

    if normalized_job.nodes < 1
        push!(report.errors, "nodes must be >= 1")
    end
    if normalized_job.cpus_per_task < 1
        push!(report.errors, "cpus_per_task must be >= 1")
    end
    if normalized_job.ntasks !== nothing && normalized_job.ntasks < 1
        push!(report.errors, "ntasks must be >= 1 when provided")
    end
    if normalized_job.ntasks_per_node !== nothing && normalized_job.ntasks_per_node < 1
        push!(report.errors, "ntasks_per_node must be >= 1 when provided")
    end

    if normalized_job.mem_gb !== nothing && normalized_job.mem_gb <= 0
        push!(report.errors, "mem_gb must be > 0 when provided")
    end
    if normalized_job.mem_per_cpu_mb !== nothing && normalized_job.mem_per_cpu_mb <= 0
        push!(report.errors, "mem_per_cpu_mb must be > 0 when provided")
    end
    if normalized_job.mem_gb !== nothing && normalized_job.mem_per_cpu_mb !== nothing
        push!(report.errors, "set either mem_gb or mem_per_cpu_mb, not both")
    end

    if normalized_job.gpus_per_node !== nothing && normalized_job.gpus_per_node <= 0
        push!(report.errors, "gpus_per_node must be > 0 when provided")
    end
    if normalized_job.gpus_per_task !== nothing && normalized_job.gpus_per_task <= 0
        push!(report.errors, "gpus_per_task must be > 0 when provided")
    end

    for (key, _) in normalized_job.env
        if match(r"^[A-Za-z_][A-Za-z0-9_]*$", key) === nothing
            push!(report.errors, "env key '$key' is invalid; expected shell-style variable names")
        end
    end

    if normalized_job.site in Set([:nersc, :lawrencium, :scg])
        if normalized_job.output_path === nothing
            push!(report.warnings, "output_path is missing; default log path will be used")
        end
        if normalized_job.error_path === nothing
            push!(report.warnings, "error_path is missing; default log path will be used")
        end
    end

    if normalized_job.site == :nersc
        _validate_nersc!(report, normalized_job, parsed_time)
    elseif normalized_job.site == :lawrencium
        _validate_lawrencium!(report, normalized_job; check_associations = check_lawrencium_associations)
    elseif normalized_job.site == :scg
        _validate_scg!(report, normalized_job)
    elseif normalized_job.site == :local
        if normalized_job.container_image === nothing
            push!(report.errors, "local backend requires container_image for docker rendering")
        end
    elseif normalized_job.site == :cloudbuild
        if normalized_job.container_image === nothing && !_has_dockerfile(normalized_job)
            push!(report.errors,
                "cloudbuild requires container_image or a Dockerfile in workdir/current directory")
        end
    end

    if normalized_job.container_engine !== nothing && !(normalized_job.container_engine in VALID_CONTAINER_ENGINES)
        push!(report.errors, "Unsupported container_engine $(normalized_job.container_engine)")
    end

    return report
end

function _validate_nersc!(report::JobSpecValidation, job::JobSpec, parsed_time::Union{Nothing, Int})
    if job.account === nothing
        push!(report.errors, "NERSC jobs require account (-A)")
    end

    qos = lowercase(something(job.qos, "regular"))
    constraint = lowercase(something(job.constraint, "cpu"))

    if qos == "jupyter"
        push!(report.errors, "NERSC jupyter QoS is not supported for arbitrary jobs")
    end

    if parsed_time !== nothing
        if qos == "debug" && parsed_time > 30 * 60
            push!(report.errors, "NERSC debug QoS max walltime is 00:30:00")
        elseif qos == "interactive" && parsed_time > 4 * 60 * 60
            push!(report.errors, "NERSC interactive QoS max walltime is 04:00:00")
        elseif qos in Set(["regular", "shared", "preempt", "premium", "overrun", "shared_overrun"])
            if parsed_time > 48 * 60 * 60
                push!(report.errors, "NERSC $qos QoS max walltime is 48:00:00")
            end
        end

        if qos == "preempt" && parsed_time < 2 * 60 * 60
            push!(report.errors, "NERSC preempt QoS requires minimum walltime of 02:00:00")
        end
    end

    if qos == "debug" && job.nodes > 8
        push!(report.errors, "NERSC debug QoS allows at most 8 nodes")
    end

    if qos == "preempt" && !job.requeue
        push!(report.warnings,
            "NERSC preempt jobs should generally set requeue=true to handle interruptions")
    end

    if !(constraint in Set(["cpu", "gpu", "gpu&hbm80g", "gpu&hbm40g"]))
        push!(report.warnings,
            "NERSC constraint '$constraint' is uncommon; expected cpu, gpu, gpu&hbm80g, or gpu&hbm40g")
    end

    gpu_requested = _job_requests_gpu(job)
    if occursin("gpu", constraint) && !gpu_requested
        push!(report.errors,
            "NERSC GPU constraints require explicit gpus_per_node or gpus_per_task")
    end

    if qos == "shared" && occursin("gpu", constraint)
        gpus_per_node = _effective_gpus_per_node(job)
        if gpus_per_node === nothing
            push!(report.errors, "NERSC shared GPU jobs must request 1 or 2 GPUs")
        elseif !(gpus_per_node in (1, 2))
            push!(report.errors,
                "NERSC shared GPU QoS only supports 1 or 2 GPUs per node")
        elseif gpus_per_node == 1
            _warn_if_unexpected_shared_shape!(report, job; recommended_cpu = 16, recommended_mem = 64)
        elseif gpus_per_node == 2
            _warn_if_unexpected_shared_shape!(report, job; recommended_cpu = 32, recommended_mem = 128)
        end
    end
end

function _validate_lawrencium!(
        report::JobSpecValidation,
        job::JobSpec;
        check_associations::Bool = true
)
    if job.partition === nothing
        push!(report.errors, "Lawrencium jobs require partition")
    end
    if job.account === nothing
        push!(report.errors, "Lawrencium jobs require account")
    end

    if job.partition !== nothing && startswith(lowercase(job.partition), "lr") && job.qos === nothing
        push!(report.warnings,
            "Lawrencium qos not set; defaulting to lr_normal for lr* partitions")
    end

    if check_associations
        associations = list_lawrencium_associations(; execute = true)
        parsed = parse_lawrencium_associations(associations)
        if isempty(parsed)
            push!(report.warnings,
                "Could not parse Lawrencium sacctmgr associations; skipping strict association check")
        elseif job.account !== nothing && job.partition !== nothing && job.qos !== nothing
            if !_association_permits(job, parsed)
                push!(report.warnings,
                    "(account, partition, qos) not observed in sacctmgr associations; verify with your admins")
            end
        end
    end
end

function _validate_scg!(report::JobSpecValidation, job::JobSpec)
    partition = lowercase(something(job.partition, "nih_s10"))

    if job.nodes != 1
        push!(report.errors, "SCG jobs must request nodes=1 by default")
    end

    if partition in Set(["batch", "nih_s10"]) && job.account === nothing
        push!(report.errors, "SCG partition '$partition' requires account=PI_SUNetID")
    end

    if !(partition in Set(["batch", "interactive", "nih_s10"]))
        push!(report.warnings,
            "SCG partition '$partition' is uncommon; expected one of batch, interactive, nih_s10")
    end
end

function _association_permits(job::JobSpec, parsed::Vector{Dict{String, Any}})::Bool
    target_account = lowercase(something(job.account, ""))
    target_partition = lowercase(something(job.partition, ""))
    target_qos = lowercase(something(job.qos, ""))

    for row in parsed
        account = lowercase(get(row, "account", ""))
        if !isempty(target_account) && account != target_account
            continue
        end

        partition = lowercase(something(get(row, "partition", nothing), ""))
        if !isempty(target_partition) && !isempty(partition) && partition != target_partition
            continue
        end

        qos_values = [lowercase(q) for q in get(row, "qos", String[])]
        if isempty(target_qos) || isempty(qos_values) || target_qos in qos_values
            return true
        end
    end

    return false
end

function _warn_if_unexpected_shared_shape!(
        report::JobSpecValidation,
        job::JobSpec;
        recommended_cpu::Int,
        recommended_mem::Int
)
    cpu_per_node = _effective_cpu_per_node(job)
    if cpu_per_node !== nothing && cpu_per_node != recommended_cpu
        push!(report.warnings,
            "NERSC shared GPU job usually maps to $recommended_cpu CPU cores per node")
    end

    if job.mem_gb !== nothing
        rounded_mem = Int(round(job.mem_gb))
        if rounded_mem != recommended_mem
            push!(report.warnings,
                "NERSC shared GPU job usually maps to $(recommended_mem)G memory per node")
        end
    end
end

function _effective_cpu_per_node(job::JobSpec)::Union{Nothing, Int}
    if job.ntasks_per_node !== nothing
        return job.ntasks_per_node * job.cpus_per_task
    end
    if job.ntasks !== nothing
        tasks_per_node = max(1, ceil(Int, job.ntasks / max(job.nodes, 1)))
        return tasks_per_node * job.cpus_per_task
    end
    return nothing
end

function _effective_tasks_per_node(job::JobSpec)::Union{Nothing, Int}
    if job.ntasks_per_node !== nothing
        return job.ntasks_per_node
    end
    if job.ntasks !== nothing
        return max(1, ceil(Int, job.ntasks / max(job.nodes, 1)))
    end
    return nothing
end

function _effective_gpus_per_node(job::JobSpec)::Union{Nothing, Int}
    if job.gpus_per_node !== nothing
        return job.gpus_per_node
    end
    if job.gpus_per_task !== nothing
        tasks_per_node = _effective_tasks_per_node(job)
        if tasks_per_node !== nothing
            return job.gpus_per_task * tasks_per_node
        end
    end
    return nothing
end

function _job_requests_gpu(job::JobSpec)::Bool
    if job.gpus_per_node !== nothing || job.gpus_per_task !== nothing
        return true
    end
    return job.constraint !== nothing && occursin("gpu", lowercase(job.constraint))
end

function normalize_job_spec(job::JobSpec)::JobSpec
    payload = job_spec_to_dict(job)

    if job.site in Set([:nersc, :lawrencium, :scg])
        if payload["output_path"] === nothing
            payload["output_path"] = joinpath(DEFAULT_SLURM_LOGDIR, "%j.%x.out")
        end
        if payload["error_path"] === nothing
            payload["error_path"] = joinpath(DEFAULT_SLURM_LOGDIR, "%j.%x.err")
        end
    end

    if job.site == :nersc
        if payload["qos"] === nothing
            payload["qos"] = "regular"
        end
        if payload["constraint"] === nothing
            payload["constraint"] = "cpu"
        end
    elseif job.site == :lawrencium
        partition = payload["partition"]
        if payload["qos"] === nothing && partition !== nothing &&
           startswith(lowercase(String(partition)), "lr")
            payload["qos"] = "lr_normal"
        end
    elseif job.site == :scg
        if payload["partition"] === nothing
            payload["partition"] = "nih_s10"
        end
    end

    return job_spec_from_dict(payload)
end

function _sanitize_filename(text::AbstractString)
    cleaned = replace(strip(String(text)), r"[^A-Za-z0-9_.-]+" => "_")
    return isempty(cleaned) ? "job" : cleaned
end

function _shell_quote(text::AbstractString)
    return "'" * replace(String(text), "'" => "'\"'\"'") * "'"
end

function _yaml_single_quote(text::AbstractString)
    return "'" * replace(String(text), "'" => "''") * "'"
end

function _render_template_text(template_text::AbstractString, values::Dict{String, String})
    rendered = String(template_text)
    for (key, value) in values
        rendered = replace(rendered, "{{$key}}" => value)
    end

    unresolved = String[]
    for match_entry in eachmatch(TEMPLATE_TOKEN_REGEX, rendered)
        push!(unresolved, match_entry.match)
    end
    if !isempty(unresolved)
        unresolved_unique = join(unique(unresolved), ", ")
        error("Unresolved template variables: $unresolved_unique")
    end

    return rendered
end

function _render_template_file(relative_path::AbstractString, values::Dict{String, String})
    template_path = joinpath(TEMPLATE_ROOT, String(relative_path))
    if !isfile(template_path)
        error("Template not found: $template_path")
    end
    template_text = read(template_path, String)
    return _render_template_text(template_text, values)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

List available template files for a given site/backend.
"""
function list_templates(site::Union{Symbol, AbstractString})::Vector{String}
    site_symbol = _normalize_site(site)
    root = if site_symbol in Set([:nersc, :lawrencium, :scg])
        joinpath(TEMPLATE_ROOT, "slurm", String(site_symbol))
    elseif site_symbol == :local
        joinpath(TEMPLATE_ROOT, "docker")
    elseif site_symbol == :cloudbuild
        joinpath(TEMPLATE_ROOT, "cloudbuild")
    else
        error("Unsupported site '$site'")
    end

    if !isdir(root)
        return String[]
    end

    templates = String[]
    for (current_root, _, files) in walkdir(root)
        for file in sort(files)
            push!(templates, relpath(joinpath(current_root, file), TEMPLATE_ROOT))
        end
    end
    sort!(templates)
    return templates
end

function _select_nersc_template(job::JobSpec)::String
    qos = lowercase(something(job.qos, "regular"))
    constraint = lowercase(something(job.constraint, "cpu"))
    gpu_requested = _job_requests_gpu(job) || occursin("gpu", constraint)

    if occursin("gpu&hbm80g", constraint)
        if _is_single_task_single_gpu(job)
            return joinpath("slurm", "nersc", "gpu_hbm80g_regular_1gpu_1task.sbatch")
        end
        return joinpath("slurm", "nersc", "gpu_hbm80g_regular_4gpus_4tasks.sbatch")
    end

    if gpu_requested
        if qos == "shared"
            gpn = _effective_gpus_per_node(job)
            if gpn == 2
                return joinpath("slurm", "nersc", "gpu_shared_2gpu.sbatch")
            end
            return joinpath("slurm", "nersc", "gpu_shared_1gpu.sbatch")
        elseif qos == "preempt"
            return joinpath("slurm", "nersc", "gpu_preempt.sbatch")
        elseif qos == "premium"
            return joinpath("slurm", "nersc", "gpu_premium.sbatch")
        elseif _is_single_task_single_gpu(job)
            return joinpath("slurm", "nersc", "gpu_regular_1gpu_1task.sbatch")
        else
            return joinpath("slurm", "nersc", "gpu_regular_4gpus_4tasks.sbatch")
        end
    end

    if qos == "debug"
        return joinpath("slurm", "nersc", "cpu_debug.sbatch")
    elseif qos == "shared" || qos == "shared_overrun"
        return joinpath("slurm", "nersc", "cpu_shared.sbatch")
    elseif qos == "preempt"
        return joinpath("slurm", "nersc", "cpu_preempt.sbatch")
    elseif qos == "premium"
        return joinpath("slurm", "nersc", "cpu_premium.sbatch")
    end

    return joinpath("slurm", "nersc", "cpu_regular.sbatch")
end

function _select_lawrencium_template(job::JobSpec)::String
    partition = lowercase(something(job.partition, ""))
    qos = lowercase(something(job.qos, ""))

    if partition == "lr4" && qos == "lr_normal"
        return joinpath("slurm", "lawrencium", "lr4_lr_normal.sbatch")
    elseif partition == "lr6" && qos == "lr_normal"
        return joinpath("slurm", "lawrencium", "lr6_lr_normal.sbatch")
    elseif partition == "lr6" && qos == "lr_debug"
        return joinpath("slurm", "lawrencium", "lr6_lr_debug.sbatch")
    elseif partition == "es1" && qos == "es_normal"
        return joinpath("slurm", "lawrencium", "es1_es_normal.sbatch")
    elseif partition == "es0" && qos in Set(["es_lowprio", "free", "es_free"])
        return joinpath("slurm", "lawrencium", "es0_es_lowprio_or_free.sbatch")
    elseif qos == "lr_lowprio"
        return joinpath("slurm", "lawrencium", "lr_lowprio_generic.sbatch")
    elseif startswith(partition, "lr")
        return joinpath("slurm", "lawrencium", "lr_lowprio_generic.sbatch")
    end

    return joinpath("slurm", "lawrencium", "lr_lowprio_generic.sbatch")
end

function _select_scg_template(job::JobSpec)::String
    partition = lowercase(something(job.partition, "nih_s10"))
    if partition == "batch"
        return joinpath("slurm", "scg", "batch.sbatch")
    elseif partition == "interactive"
        return joinpath("slurm", "scg", "interactive.salloc")
    end
    return joinpath("slurm", "scg", "nih_s10.sbatch")
end

function _selected_template(job::JobSpec)::String
    if job.template_name !== nothing
        template_name = String(job.template_name)
        if startswith(template_name, "slurm/") || startswith(template_name, "docker/") ||
           startswith(template_name, "cloudbuild/")
            return template_name
        end

        if job.site in Set([:nersc, :lawrencium, :scg])
            return joinpath("slurm", String(job.site), template_name)
        elseif job.site == :local
            return joinpath("docker", template_name)
        elseif job.site == :cloudbuild
            return joinpath("cloudbuild", template_name)
        end
    end

    if job.site == :nersc
        return _select_nersc_template(job)
    elseif job.site == :lawrencium
        return _select_lawrencium_template(job)
    elseif job.site == :scg
        return _select_scg_template(job)
    elseif job.site == :local
        return joinpath("docker", "run.sh")
    elseif job.site == :cloudbuild
        return joinpath("cloudbuild", "cloudbuild.yaml")
    end

    error("Unable to resolve template for site $(job.site)")
end

function _is_single_task_single_gpu(job::JobSpec)::Bool
    ntasks = something(job.ntasks, 1)
    gpus_per_node = _effective_gpus_per_node(job)
    return ntasks == 1 && gpus_per_node == 1
end

function _mail_type_string(job::JobSpec)
    if job.mail_type === nothing
        return nothing
    end
    return join(job.mail_type, ",")
end

function _slurm_mem_directive(job::JobSpec)
    if job.mem_per_cpu_mb !== nothing
        return "#SBATCH --mem-per-cpu=$(job.mem_per_cpu_mb)"
    end

    if job.mem_gb !== nothing
        mem_gb = Int(ceil(job.mem_gb))
        return "#SBATCH --mem=$(mem_gb)G"
    end

    return nothing
end

function _build_sbatch_directives(job::JobSpec)::String
    lines = String[]
    push!(lines, "#SBATCH --job-name=$(job.job_name)")

    if job.mail_user !== nothing
        push!(lines, "#SBATCH --mail-user=$(job.mail_user)")
    end

    mail_type = _mail_type_string(job)
    if mail_type !== nothing
        push!(lines, "#SBATCH --mail-type=$(mail_type)")
    end

    if job.error_path !== nothing
        push!(lines, "#SBATCH --error=$(job.error_path)")
    end
    if job.output_path !== nothing
        push!(lines, "#SBATCH --output=$(job.output_path)")
    end

    if job.partition !== nothing
        push!(lines, "#SBATCH --partition=$(job.partition)")
    end
    if job.qos !== nothing
        push!(lines, "#SBATCH --qos=$(job.qos)")
    end
    if job.account !== nothing
        push!(lines, "#SBATCH --account=$(job.account)")
    end

    push!(lines, "#SBATCH --nodes=$(job.nodes)")

    if job.ntasks !== nothing
        push!(lines, "#SBATCH --ntasks=$(job.ntasks)")
    end
    if job.ntasks_per_node !== nothing
        push!(lines, "#SBATCH --ntasks-per-node=$(job.ntasks_per_node)")
    end

    push!(lines, "#SBATCH --time=$(job.time_limit)")
    push!(lines, "#SBATCH --cpus-per-task=$(job.cpus_per_task)")

    memory_line = _slurm_mem_directive(job)
    if memory_line !== nothing
        push!(lines, memory_line)
    end

    if job.constraint !== nothing
        push!(lines, "#SBATCH --constraint=$(job.constraint)")
    end
    if job.gpus_per_node !== nothing
        push!(lines, "#SBATCH --gpus-per-node=$(job.gpus_per_node)")
    end
    if job.gpus_per_task !== nothing
        push!(lines, "#SBATCH --gpus-per-task=$(job.gpus_per_task)")
    end
    if job.gpu_bind !== nothing
        push!(lines, "#SBATCH --gpu-bind=$(job.gpu_bind)")
    end

    if job.requeue
        push!(lines, "#SBATCH --requeue")
    end

    return join(lines, "\n")
end

function _build_prelude_block(job::JobSpec)::String
    lines = String[]

    if job.workdir !== nothing
        push!(lines, "mkdir -p " * _shell_quote(job.workdir))
        push!(lines, "cd " * _shell_quote(job.workdir))
    end

    if !isempty(job.modules)
        for module_name in job.modules
            push!(lines, "module load " * _shell_quote(module_name))
        end
    end

    if !isempty(job.env)
        for key in sort(collect(keys(job.env)))
            push!(lines, "export $key=" * _shell_quote(job.env[key]))
        end
    end

    if isempty(lines)
        return "# No module or environment setup requested"
    end

    return join(lines, "\n")
end

function _apptainer_image_ref(image::AbstractString)
    if occursin("://", String(image))
        return String(image)
    end
    return "docker://" * String(image)
end

function _effective_container_engine(job::JobSpec)::Symbol
    if job.container_engine !== nothing
        return job.container_engine
    end
    if job.site == :nersc
        return :shifter
    elseif job.site in Set([:lawrencium, :scg])
        return :apptainer
    end
    return :docker
end

function _build_mount_flags(mounts::Vector{String}, style::Symbol)::Vector{String}
    flags = String[]
    for mount in mounts
        spec = strip(mount)
        if isempty(spec)
            continue
        end

        if style == :docker
            push!(flags, "--volume " * _shell_quote(spec))
        elseif style == :apptainer
            push!(flags, "--bind " * _shell_quote(spec))
        elseif style == :shifter
            push!(flags, "--volume=" * _shell_quote(spec))
        elseif style == :podman_hpc
            push!(flags, "-v " * _shell_quote(spec))
        end
    end
    return flags
end

function _build_container_exec(job::JobSpec)::String
    if job.container_image === nothing
        return "bash -lc \"\$JOB_CMD\""
    end

    engine = _effective_container_engine(job)
    gpu_requested = _job_requests_gpu(job)

    if engine == :apptainer
        segments = ["apptainer exec"]
        if gpu_requested
            push!(segments, "--nv")
        end
        append!(segments, _build_mount_flags(job.container_mounts, :apptainer))
        if job.container_workdir !== nothing
            push!(segments, "--pwd " * _shell_quote(job.container_workdir))
        end
        append!(segments, [_shell_quote(arg) for arg in job.container_args])
        push!(segments, _shell_quote(_apptainer_image_ref(job.container_image)))
        push!(segments, "bash -lc \"\$JOB_CMD\"")
        return join(segments, " ")
    elseif engine == :shifter
        segments = ["shifter"]
        if gpu_requested
            push!(segments, "--gpu")
        end
        append!(segments, _build_mount_flags(job.container_mounts, :shifter))
        if job.container_workdir !== nothing
            push!(segments, "--workdir=" * _shell_quote(job.container_workdir))
        end
        append!(segments, [_shell_quote(arg) for arg in job.container_args])
        push!(segments, "--image=" * _shell_quote(job.container_image))
        push!(segments, "bash -lc \"\$JOB_CMD\"")
        return join(segments, " ")
    elseif engine == :podman_hpc
        segments = ["podman-hpc run --rm"]
        if gpu_requested
            push!(segments, "--hooks-dir=/usr/share/containers/oci/hooks.d")
        end
        append!(segments, _build_mount_flags(job.container_mounts, :podman_hpc))
        if job.container_workdir !== nothing
            push!(segments, "--workdir " * _shell_quote(job.container_workdir))
        end
        if !isempty(job.env)
            for (key, value) in sort(collect(job.env))
                push!(segments, "--env " * _shell_quote("$key=$value"))
            end
        end
        append!(segments, [_shell_quote(arg) for arg in job.container_args])
        push!(segments, _shell_quote(job.container_image))
        push!(segments, "bash -lc \"\$JOB_CMD\"")
        return join(segments, " ")
    end

    segments = ["docker run --rm"]
    if gpu_requested
        push!(segments, "--gpus all")
    end
    append!(segments, _build_mount_flags(job.container_mounts, :docker))
    if job.container_workdir !== nothing
        push!(segments, "--workdir " * _shell_quote(job.container_workdir))
    end
    if !isempty(job.env)
        for (key, value) in sort(collect(job.env))
            push!(segments, "--env " * _shell_quote("$key=$value"))
        end
    end
    append!(segments, [_shell_quote(arg) for arg in job.container_args])
    push!(segments, _shell_quote(job.container_image))
    push!(segments, "bash -lc \"\$JOB_CMD\"")
    return join(segments, " ")
end

function _build_run_block(job::JobSpec)::String
    lines = String[]
    push!(lines, "JOB_CMD=" * _shell_quote(job.cmd))

    if job.container_image === nothing
        push!(lines, "bash -lc \"\$JOB_CMD\"")
    else
        push!(lines, _build_container_exec(job))
    end

    return join(lines, "\n")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return standardized queue/accounting/resource monitoring hints.
"""
function monitoring_hints(; jobid::AbstractString = "<jobid>")::Vector{String}
    return [
        "squeue -u \$USER",
        "scontrol show job $jobid",
        "sacct -j $jobid --format=JobID,JobName,AllocCPUS,Elapsed,State,ExitCode,MaxRSS",
        "sstat -j $jobid.batch --format=JobID,MaxRSS,MaxVMSize,AvgCPU",
        "seff $jobid  # if seff is available",
        "nvidia-smi",
        "nvidia-smi topo -m  # optional topology check"
    ]
end

function _monitoring_hints_comment_block(job::JobSpec)::String
    lines = ["# Monitoring hints:"]
    for command in monitoring_hints()
        push!(lines, "#   $command")
    end
    if job.site == :nersc && _job_requests_gpu(job)
        push!(lines, "# NERSC note: use --gpu-bind=none when all tasks should see all GPUs")
    end
    return join(lines, "\n")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Render a batch-mode SLURM script from a `JobSpec`.
"""
function render_sbatch(job::JobSpec)::String
    normalized_job = normalize_job_spec(job)
    report = validate(normalized_job)
    if !isempty(report.errors)
        error("JobSpec validation failed:\n" * join(report.errors, "\n"))
    end

    if _is_interactive_job(normalized_job)
        error("Interactive jobs should use render_salloc(job) instead of render_sbatch(job)")
    end

    template_rel = _selected_template(normalized_job)
    if endswith(template_rel, ".salloc")
        error("Template $template_rel is interactive; call render_salloc(job)")
    end

    values = Dict(
        "TEMPLATE_ID" => replace(template_rel, "\\" => "/"),
        "SITE_NAME" => String(normalized_job.site),
        "SBATCH_DIRECTIVES" => _build_sbatch_directives(normalized_job),
        "PRELUDE_BLOCK" => _build_prelude_block(normalized_job),
        "RUN_BLOCK" => _build_run_block(normalized_job),
        "MONITORING_HINTS_BLOCK" => _monitoring_hints_comment_block(normalized_job)
    )

    return _render_template_file(template_rel, values)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Render an interactive allocation command (`salloc ... srun --pty`) where supported.
"""
function render_salloc(job::JobSpec)::String
    normalized_job = normalize_job_spec(job)
    report = validate(normalized_job)
    if !isempty(report.errors)
        error("JobSpec validation failed:\n" * join(report.errors, "\n"))
    end

    if normalized_job.site == :nersc
        qos = lowercase(something(normalized_job.qos, "interactive"))
        if qos != "interactive"
            error("NERSC interactive commands require qos=interactive")
        end

        segments = ["salloc"]
        push!(segments, "-N $(normalized_job.nodes)")
        if normalized_job.ntasks !== nothing
            push!(segments, "-n $(normalized_job.ntasks)")
        end
        if normalized_job.ntasks_per_node !== nothing
            push!(segments, "--ntasks-per-node=$(normalized_job.ntasks_per_node)")
        end
        push!(segments, "-c $(normalized_job.cpus_per_task)")
        push!(segments, "-t $(normalized_job.time_limit)")
        push!(segments, "-q interactive")
        if normalized_job.account !== nothing
            push!(segments, "-A $(normalized_job.account)")
        end
        if normalized_job.constraint !== nothing
            constraint = String(normalized_job.constraint)
            if occursin("&", constraint)
                push!(segments, "-C \"$constraint\"")
            else
                push!(segments, "-C $(constraint)")
            end
        end
        if normalized_job.gpus_per_node !== nothing
            push!(segments, "--gpus-per-node=$(normalized_job.gpus_per_node)")
        end
        if normalized_job.gpus_per_task !== nothing
            push!(segments, "--gpus-per-task=$(normalized_job.gpus_per_task)")
        end
        if normalized_job.mem_gb !== nothing
            push!(segments, "--mem=$(Int(ceil(normalized_job.mem_gb)))G")
        elseif normalized_job.mem_per_cpu_mb !== nothing
            push!(segments, "--mem-per-cpu=$(normalized_job.mem_per_cpu_mb)")
        end
        push!(segments, "srun --pty bash -l")
        return join(segments, " ")
    elseif normalized_job.site == :scg
        partition = lowercase(something(normalized_job.partition, "interactive"))
        if partition != "interactive"
            error("SCG interactive command requires partition=interactive")
        end

        segments = ["salloc"]
        push!(segments, "--partition=interactive")
        push!(segments, "--nodes=1")
        if normalized_job.ntasks !== nothing
            push!(segments, "--ntasks=$(normalized_job.ntasks)")
        end
        push!(segments, "--cpus-per-task=$(normalized_job.cpus_per_task)")
        push!(segments, "--time=$(normalized_job.time_limit)")
        if normalized_job.account !== nothing
            push!(segments, "--account=$(normalized_job.account)")
        end
        if normalized_job.mem_gb !== nothing
            push!(segments, "--mem=$(Int(ceil(normalized_job.mem_gb)))G")
        elseif normalized_job.mem_per_cpu_mb !== nothing
            push!(segments, "--mem-per-cpu=$(normalized_job.mem_per_cpu_mb)")
        end
        push!(segments, "srun --pty bash -l")

        values = Dict(
            "TEMPLATE_ID" => "slurm/scg/interactive.salloc",
            "SALLOC_COMMAND" => join(segments, " ")
        )
        return _render_template_file(joinpath("slurm", "scg", "interactive.salloc"), values)
    end

    error("render_salloc is only supported for site=:nersc or site=:scg")
end

function _default_sbatch_path(job::JobSpec)::String
    mkpath(DEFAULT_SLURM_SCRIPTDIR)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
    filename = "$(timestamp)-$(_sanitize_filename(job.job_name)).sbatch"
    return joinpath(DEFAULT_SLURM_SCRIPTDIR, filename)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a rendered sbatch script to disk and return its path.
"""
function write_sbatch(job::JobSpec; path::Union{Nothing, String} = nothing)::String
    script = render_sbatch(job)
    output_path = something(path, _default_sbatch_path(normalize_job_spec(job)))
    mkpath(dirname(output_path))
    write(output_path, script)
    chmod(output_path, 0o755)
    return output_path
end

function _docker_run_command(job::JobSpec)::String
    if job.container_image === nothing
        error("container_image is required for docker rendering")
    end

    segments = ["docker run --rm"]
    if _job_requests_gpu(job)
        push!(segments, "--gpus all")
    end

    append!(segments, _build_mount_flags(job.container_mounts, :docker))

    if job.container_workdir !== nothing
        push!(segments, "--workdir " * _shell_quote(job.container_workdir))
    end

    if !isempty(job.env)
        for (key, value) in sort(collect(job.env))
            push!(segments, "--env " * _shell_quote("$key=$value"))
        end
    end

    append!(segments, [_shell_quote(arg) for arg in job.container_args])
    push!(segments, _shell_quote(job.container_image))
    push!(segments, "bash -lc " * _shell_quote(job.cmd))

    return join(segments, " ")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Render a Docker execution script from a `JobSpec`.
"""
function render_docker_run(job::JobSpec)::String
    normalized_job = normalize_job_spec(job)
    report = validate(normalized_job)
    if !isempty(report.errors)
        error("JobSpec validation failed:\n" * join(report.errors, "\n"))
    end

    template_rel = _selected_template(normalized_job)
    command_lines = String[]
    if normalized_job.workdir !== nothing
        push!(command_lines, "cd " * _shell_quote(normalized_job.workdir))
    end
    push!(command_lines, _docker_run_command(normalized_job))

    values = Dict(
        "TEMPLATE_ID" => replace(template_rel, "\\" => "/"),
        "DOCKER_COMMAND" => join(command_lines, "\n")
    )

    return _render_template_file(template_rel, values)
end

function _has_dockerfile(job::JobSpec)::Bool
    base = job.workdir === nothing ? pwd() : job.workdir
    return isfile(joinpath(base, "Dockerfile"))
end

function _cloudbuild_substitutions(job::JobSpec)
    image_uri = if job.container_image !== nothing
        job.container_image
    else
        "gcr.io/\$PROJECT_ID/\${_IMAGE_NAME}:\${_IMAGE_TAG}"
    end

    substitutions = [
        "  _IMAGE_NAME: mycelia-job",
        "  _IMAGE_TAG: latest",
        "  _IMAGE_URI: " * _yaml_single_quote(image_uri)
    ]

    return image_uri, substitutions
end

function _cloudbuild_content(job::JobSpec)::String
    has_dockerfile = _has_dockerfile(job)
    _, substitutions = _cloudbuild_substitutions(job)

    lines = String[]
    push!(lines, "steps:")

    if has_dockerfile
        append!(lines, [
            "- id: build-image",
            "  name: gcr.io/cloud-builders/docker",
            "  args: ['build', '-t', '\${_IMAGE_URI}', '.']",
            "- id: push-image",
            "  name: gcr.io/cloud-builders/docker",
            "  args: ['push', '\${_IMAGE_URI}']"
        ])
    elseif job.container_image === nothing
        error("cloudbuild backend requires container_image or Dockerfile")
    end

    push!(lines, "- id: run-job")
    push!(lines, "  name: \${_IMAGE_URI}")
    push!(lines, "  entrypoint: bash")
    push!(lines, "  args: ['-lc', " * _yaml_single_quote(job.cmd) * "]")

    if !isempty(job.env)
        push!(lines, "  env:")
        for (key, value) in sort(collect(job.env))
            push!(lines, "  - " * _yaml_single_quote("$key=$value"))
        end
    end

    run_dir = if job.container_workdir !== nothing
        job.container_workdir
    elseif job.workdir !== nothing
        job.workdir
    else
        nothing
    end
    if run_dir !== nothing
        push!(lines, "  dir: " * _yaml_single_quote(run_dir))
    end

    if has_dockerfile
        push!(lines, "images:")
        push!(lines, "- \${_IMAGE_URI}")
    end

    push!(lines, "options:")
    push!(lines, "  logging: CLOUD_LOGGING_ONLY")
    push!(lines, "substitutions:")
    append!(lines, substitutions)

    return join(lines, "\n")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Render a Google Cloud Build config YAML from a `JobSpec`.
"""
function render_cloudbuild(job::JobSpec)::String
    normalized_job = normalize_job_spec(job)
    report = validate(normalized_job)
    if !isempty(report.errors)
        error("JobSpec validation failed:\n" * join(report.errors, "\n"))
    end

    template_rel = _selected_template(normalized_job)
    values = Dict(
        "TEMPLATE_ID" => replace(template_rel, "\\" => "/"),
        "CLOUDBUILD_CONTENT" => _cloudbuild_content(normalized_job)
    )

    return _render_template_file(template_rel, values)
end

function _is_interactive_job(job::JobSpec)::Bool
    qos = lowercase(something(job.qos, ""))
    partition = lowercase(something(job.partition, ""))
    if job.site == :nersc && qos == "interactive"
        return true
    end
    return job.site == :scg && partition == "interactive"
end

function _extract_sbatch_job_id(stdout_text::AbstractString)
    match_entry = match(r"Submitted batch job (\d+)", String(stdout_text))
    return match_entry === nothing ? nothing : match_entry.captures[1]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a `JobSpec` to its target backend.

- `dry_run=true`: print rendered artifact + exact submit command, no execution.
- `dry_run=false`: write artifact and execute submit command where applicable.
"""
function submit(
        job::JobSpec;
        dry_run::Bool = true,
        path::Union{Nothing, String} = nothing,
        allow_invalid::Bool = false,
        io::IO = stdout,
        execute_interactive::Bool = false,
        execute_cloudbuild::Bool = false
)::SubmitResult
    normalized_job = normalize_job_spec(job)
    report = validate(normalized_job)

    if !isempty(report.errors) && !allow_invalid
        return SubmitResult(
            ok = false,
            dry_run = dry_run,
            site = normalized_job.site,
            backend = :validation,
            warnings = report.warnings,
            errors = report.errors
        )
    end

    if normalized_job.site == :local
        artifact_text = render_docker_run(normalized_job)
        artifact_path = something(path,
            joinpath(pwd(), "docker-run-$(_sanitize_filename(normalized_job.job_name)).sh"))
        submit_command = "bash " * _shell_quote(artifact_path)

        if dry_run
            _print_dry_run(io, artifact_text, submit_command)
            return SubmitResult(
                ok = true,
                dry_run = true,
                site = normalized_job.site,
                backend = :docker,
                artifact_path = artifact_path,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = report.warnings,
                errors = report.errors
            )
        end

        mkpath(dirname(artifact_path))
        write(artifact_path, artifact_text)
        chmod(artifact_path, 0o755)

        try
            output = read(`bash $artifact_path`, String)
            return SubmitResult(
                ok = true,
                dry_run = false,
                site = normalized_job.site,
                backend = :docker,
                artifact_path = artifact_path,
                artifact_text = artifact_text,
                submit_command = submit_command,
                stdout = output,
                warnings = report.warnings,
                errors = report.errors
            )
        catch e
            return SubmitResult(
                ok = false,
                dry_run = false,
                site = normalized_job.site,
                backend = :docker,
                artifact_path = artifact_path,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = report.warnings,
                errors = vcat(report.errors, [string(e)])
            )
        end
    elseif normalized_job.site == :cloudbuild
        artifact_text = render_cloudbuild(normalized_job)
        artifact_path = something(path,
            joinpath(pwd(), "cloudbuild-$(_sanitize_filename(normalized_job.job_name)).yaml"))
        submit_command = "gcloud builds submit --config $(artifact_path) ."

        if dry_run
            _print_dry_run(io, artifact_text, submit_command)
            return SubmitResult(
                ok = true,
                dry_run = true,
                site = normalized_job.site,
                backend = :cloudbuild,
                artifact_path = artifact_path,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = report.warnings,
                errors = report.errors
            )
        end

        mkpath(dirname(artifact_path))
        write(artifact_path, artifact_text)

        if !execute_cloudbuild
            return SubmitResult(
                ok = true,
                dry_run = false,
                site = normalized_job.site,
                backend = :cloudbuild,
                artifact_path = artifact_path,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = vcat(report.warnings, [
                    "Cloud Build config written; set execute_cloudbuild=true to run gcloud submit"
                ]),
                errors = report.errors
            )
        end

        if Sys.which("gcloud") === nothing
            return SubmitResult(
                ok = false,
                dry_run = false,
                site = normalized_job.site,
                backend = :cloudbuild,
                artifact_path = artifact_path,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = report.warnings,
                errors = vcat(report.errors, ["gcloud is not available on PATH"])
            )
        end

        try
            output = read(`gcloud builds submit --config=$artifact_path .`, String)
            return SubmitResult(
                ok = true,
                dry_run = false,
                site = normalized_job.site,
                backend = :cloudbuild,
                artifact_path = artifact_path,
                artifact_text = artifact_text,
                submit_command = submit_command,
                stdout = output,
                warnings = report.warnings,
                errors = report.errors
            )
        catch e
            return SubmitResult(
                ok = false,
                dry_run = false,
                site = normalized_job.site,
                backend = :cloudbuild,
                artifact_path = artifact_path,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = report.warnings,
                errors = vcat(report.errors, [string(e)])
            )
        end
    end

    if _is_interactive_job(normalized_job)
        artifact_text = render_salloc(normalized_job)
        submit_command = artifact_text

        if dry_run
            _print_dry_run(io, artifact_text, submit_command)
            return SubmitResult(
                ok = true,
                dry_run = true,
                site = normalized_job.site,
                backend = :salloc,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = report.warnings,
                errors = report.errors
            )
        end

        if !execute_interactive
            return SubmitResult(
                ok = true,
                dry_run = false,
                site = normalized_job.site,
                backend = :salloc,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = vcat(report.warnings, [
                    "Interactive command generated; set execute_interactive=true to run it"
                ]),
                errors = report.errors
            )
        end

        try
            output = read(`bash -lc $submit_command`, String)
            return SubmitResult(
                ok = true,
                dry_run = false,
                site = normalized_job.site,
                backend = :salloc,
                artifact_text = artifact_text,
                submit_command = submit_command,
                stdout = output,
                warnings = report.warnings,
                errors = report.errors
            )
        catch e
            return SubmitResult(
                ok = false,
                dry_run = false,
                site = normalized_job.site,
                backend = :salloc,
                artifact_text = artifact_text,
                submit_command = submit_command,
                warnings = report.warnings,
                errors = vcat(report.errors, [string(e)])
            )
        end
    end

    artifact_text = render_sbatch(normalized_job)
    artifact_path = something(path, _default_sbatch_path(normalized_job))
    submit_command = "sbatch " * _shell_quote(artifact_path)

    if dry_run
        _print_dry_run(io, artifact_text, submit_command)
        return SubmitResult(
            ok = true,
            dry_run = true,
            site = normalized_job.site,
            backend = :sbatch,
            artifact_path = artifact_path,
            artifact_text = artifact_text,
            submit_command = submit_command,
            warnings = report.warnings,
            errors = report.errors
        )
    end

    mkpath(dirname(artifact_path))
    write(artifact_path, artifact_text)
    chmod(artifact_path, 0o755)

    try
        output = read(`sbatch $artifact_path`, String)
        job_id = _extract_sbatch_job_id(output)
        return SubmitResult(
            ok = true,
            dry_run = false,
            site = normalized_job.site,
            backend = :sbatch,
            artifact_path = artifact_path,
            artifact_text = artifact_text,
            submit_command = submit_command,
            job_id = job_id,
            stdout = output,
            warnings = report.warnings,
            errors = report.errors
        )
    catch e
        return SubmitResult(
            ok = false,
            dry_run = false,
            site = normalized_job.site,
            backend = :sbatch,
            artifact_path = artifact_path,
            artifact_text = artifact_text,
            submit_command = submit_command,
            warnings = report.warnings,
            errors = vcat(report.errors, [string(e)])
        )
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compatibility alias for `submit(job; kwargs...)`.
"""
function submit_sbatch(job::JobSpec; kwargs...)
    return submit(job; kwargs...)
end

function _print_dry_run(io::IO, artifact_text::AbstractString, submit_command::AbstractString)
    println(io, "# --- BEGIN RENDERED ARTIFACT ---")
    println(io, artifact_text)
    println(io, "# --- END RENDERED ARTIFACT ---")
    println(io, "# Submit command")
    println(io, submit_command)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Best-effort helper for querying Lawrencium account/partition/qos associations.
"""
function list_lawrencium_associations(; execute::Bool = true)::String
    cmd = "sacctmgr show association -p user=\$USER"
    if !execute
        return cmd
    end

    if Sys.which("sacctmgr") === nothing
        return cmd
    end

    user = get(ENV, "USER", "")
    if isempty(user)
        return cmd
    end

    try
        return read(`sacctmgr show association -p user=$user`, String)
    catch
        return cmd
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Best-effort helper for listing Lawrencium QoS limits.
"""
function list_lawrencium_qos_limits(; execute::Bool = true)::String
    cmd = "sacctmgr show qos -p format=name,maxtres,maxwall,mintres"
    if !execute
        return cmd
    end

    if Sys.which("sacctmgr") === nothing
        return cmd
    end

    try
        return read(`sacctmgr show qos -p format=name,maxtres,maxwall,mintres`, String)
    catch
        return cmd
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse `sacctmgr ... -p` pipe-delimited output into row dictionaries.
"""
function parse_sacctmgr_pipe_table(raw::AbstractString)::Vector{Dict{String, String}}
    lines = [strip(line) for line in split(String(raw), '\n') if !isempty(strip(line))]
    if isempty(lines)
        return Dict{String, String}[]
    end

    header_tokens = split(first(lines), '|'; keepempty = true)
    if !isempty(header_tokens) && isempty(strip(last(header_tokens)))
        pop!(header_tokens)
    end
    header_parts = [lowercase(strip(x)) for x in header_tokens]
    if isempty(header_parts) || all(isempty, header_parts)
        return Dict{String, String}[]
    end

    rows = Dict{String, String}[]
    for line in Iterators.drop(lines, 1)
        parts = split(line, '|'; keepempty = true)
        if !isempty(parts) && isempty(strip(last(parts)))
            pop!(parts)
        end
        if length(parts) != length(header_parts)
            continue
        end
        row = Dict{String, String}()
        for (header, value) in zip(header_parts, parts)
            row[header] = strip(value)
        end
        push!(rows, row)
    end

    return rows
end

function _row_get(row::Dict{String, String}, keys::Vector{String})
    for key in keys
        if haskey(row, key)
            value = strip(row[key])
            if !isempty(value)
                return value
            end
        end
    end
    return nothing
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse Lawrencium association records into account/partition/qos triples.
"""
function parse_lawrencium_associations(raw::AbstractString)::Vector{Dict{String, Any}}
    rows = parse_sacctmgr_pipe_table(raw)
    parsed = Dict{String, Any}[]

    for row in rows
        account = _row_get(row, ["account", "acct"])
        partition = _row_get(row, ["partition", "partitionname"])
        qos_raw = _row_get(row, ["qos", "defaultqos"])
        qos_values = if qos_raw === nothing
            String[]
        else
            [strip(value) for value in split(qos_raw, ',') if !isempty(strip(value))]
        end

        push!(parsed, Dict(
            "account" => something(account, ""),
            "partition" => partition,
            "qos" => qos_values
        ))
    end

    return parsed
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Summarize scheduler/accounting state for a job id using `scontrol` and `sacct`.
"""
function summarize_job(jobid::Union{Int, AbstractString}; io::IO = stdout)
    jobid_string = string(jobid)
    if match(r"^[0-9]+([.][A-Za-z0-9_]+)?$", jobid_string) === nothing
        error("jobid must be numeric (optionally with a .step suffix)")
    end

    outputs = Dict{String, String}()

    if Sys.which("scontrol") !== nothing
        try
            outputs["scontrol"] = read(`scontrol show job $jobid_string`, String)
        catch e
            outputs["scontrol"] = "ERROR: $(e)"
        end
    else
        outputs["scontrol"] = "scontrol not found"
    end

    if Sys.which("sacct") !== nothing
        try
            outputs["sacct"] = read(`sacct -j $jobid_string --format=JobID,JobName,AllocCPUS,Elapsed,State,ExitCode,MaxRSS`, String)
        catch e
            outputs["sacct"] = "ERROR: $(e)"
        end
    else
        outputs["sacct"] = "sacct not found"
    end

    if Sys.which("sstat") !== nothing
        try
            outputs["sstat"] = read(`sstat -j $(jobid_string).batch --format=JobID,MaxRSS,MaxVMSize,AvgCPU`, String)
        catch e
            outputs["sstat"] = "ERROR: $(e)"
        end
    else
        outputs["sstat"] = "sstat not found"
    end

    if Sys.which("seff") !== nothing
        try
            outputs["seff"] = read(`seff $jobid_string`, String)
        catch e
            outputs["seff"] = "ERROR: $(e)"
        end
    else
        outputs["seff"] = "seff not found"
    end

    println(io, "=== Job Summary: $jobid_string ===")
    println(io, "")
    for label in ["scontrol", "sacct", "sstat", "seff"]
        println(io, "--- $label ---")
        println(io, get(outputs, label, "missing"))
        println(io, "")
    end

    return outputs
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Print NERSC login-node cgroup limits reminder (56 GB memory / 25% CPU throttling).
"""
function nersc_login_node_limits_warning(; io::IO = stderr)
    println(io,
        "NERSC login-node limits: 56 GB RAM per user and aggregate CPU throttling to 25%. Use SLURM allocations for compute-heavy work.")
    return nothing
end
