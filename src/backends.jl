"""
$(DocStringExtensions.TYPEDSIGNATURES)

Abstract backend configuration used to route commands to local or SLURM execution.
"""
abstract type AbstractBackend end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Execute commands directly in the current shell.
"""
struct LocalBackend <: AbstractBackend end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

SLURM backend configuration for SCG or Lawrencium submissions.
"""
Base.@kwdef struct SLURMBackend <: AbstractBackend
    site::Symbol = :scg
    account::Union{Nothing, String} = nothing
    partition::Union{Nothing, String} = nothing
    mail_user::Union{Nothing, String} = nothing
    time_limit::String = "7-00:00:00"
    mem_gb::Int = 96
    cpus_per_task::Int = SCG_THREADS
    qos::Union{Nothing, String} = nothing
    dry_run::Bool = false
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Resolve a backend setting into a concrete backend configuration.
"""
function resolve_backend(::Nothing)::AbstractBackend
    return LocalBackend()
end

function resolve_backend(backend::AbstractBackend)::AbstractBackend
    return backend
end

function resolve_backend(backend::Symbol)::AbstractBackend
    if backend === :local
        return LocalBackend()
    elseif backend in (:scg, :slurm, :stanford_scg)
        return SLURMBackend(site = :scg)
    elseif backend in (:lawrencium, :lrc)
        return SLURMBackend(
            site = :lawrencium,
            partition = "lr6",
            qos = "lr_normal",
            time_limit = "3-00:00:00",
            mem_gb = 64,
            cpus_per_task = LRC_THREADS
        )
    end
    error("Unsupported backend symbol: $(backend)")
end

function resolve_backend(backend)
    error("Unsupported backend: $(typeof(backend))")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return the default backend configuration (`LocalBackend()`).
"""
function default_backend()::AbstractBackend
    return LocalBackend()
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a Julia `Cmd` to a shell command string.
"""
function build_shell_command(cmd::Cmd)::String
    rendered = string(cmd)
    if startswith(rendered, '`') && endswith(rendered, '`')
        return rendered[2:(end - 1)]
    end
    return rendered
end

function build_shell_command(cmd::AbstractString)::String
    return String(cmd)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Execute a command with a local backend.
"""
function execute(::LocalBackend, cmd::Cmd)::Bool
    run(cmd)
    return true
end

function execute(::LocalBackend, cmd::AbstractString)::Bool
    run(`bash -lc $(String(cmd))`)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a command with a configured SLURM backend.
"""
function execute(backend::SLURMBackend, cmd::Union{Cmd, AbstractString}; job_name::AbstractString)
    command_text = build_shell_command(cmd)
    if backend.site in (:scg, :stanford_scg)
        return scg_sbatch(
            job_name = String(job_name),
            mail_user = backend.mail_user,
            partition = something(backend.partition, "nih_s10"),
            account = backend.account,
            time = backend.time_limit,
            cpus_per_task = backend.cpus_per_task,
            mem_gb = backend.mem_gb,
            cmd = command_text,
            dry_run = backend.dry_run
        )
    elseif backend.site in (:lawrencium, :lrc)
        return lawrencium_sbatch(
            job_name = String(job_name),
            mail_user = backend.mail_user,
            partition = something(backend.partition, "lr6"),
            qos = something(backend.qos, "lr_normal"),
            account = backend.account,
            time = backend.time_limit,
            cpus_per_task = backend.cpus_per_task,
            mem_gb = backend.mem_gb,
            cmd = command_text,
            dry_run = backend.dry_run
        )
    end
    error("Unsupported SLURM backend site: $(backend.site)")
end
