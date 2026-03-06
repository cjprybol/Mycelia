"""
$(DocStringExtensions.TYPEDSIGNATURES)

Execution backends for tool wrappers.
"""
abstract type AbstractExecutor end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Execute jobs directly in the current shell session.
"""
struct LocalExecutor <: AbstractExecutor end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit jobs through the SLURM template backend.
"""
Base.@kwdef struct SlurmExecutor <: AbstractExecutor
    dry_run::Bool = false
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

SLURM executor carrying full backend configuration overrides.
"""
Base.@kwdef struct ConfiguredSlurmExecutor <: AbstractExecutor
    backend::SLURMBackend = SLURMBackend()
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Collect jobs without executing them.
"""
mutable struct CollectExecutor <: AbstractExecutor
    jobs::Vector{JobSpec}
end

CollectExecutor() = CollectExecutor(JobSpec[])

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Collect jobs and store `submit(...; dry_run=true)` results.
"""
mutable struct DryRunExecutor <: AbstractExecutor
    jobs::Vector{JobSpec}
    results::Vector{SubmitResult}
end

DryRunExecutor() = DryRunExecutor(JobSpec[], SubmitResult[])

function resolve_executor(executor::AbstractExecutor)::AbstractExecutor
    return executor
end

function resolve_executor(::LocalBackend)::AbstractExecutor
    return LocalExecutor()
end

function resolve_executor(backend::SLURMBackend)::AbstractExecutor
    return ConfiguredSlurmExecutor(backend = backend)
end

function resolve_executor(::Nothing)::AbstractExecutor
    return LocalExecutor()
end

function resolve_executor(executor::Symbol)::AbstractExecutor
    if executor === :local
        return LocalExecutor()
    elseif executor === :slurm
        return SlurmExecutor()
    elseif executor === :collect
        return CollectExecutor()
    elseif executor === :dry_run || executor === :dryrun
        return DryRunExecutor()
    end
    error("Unsupported executor symbol: $(executor)")
end

function resolve_executor(executor)
    error("Unsupported executor: $(typeof(executor))")
end

function command_string(cmd::Cmd)::String
    rendered = string(cmd)
    if startswith(rendered, '`') && endswith(rendered, '`')
        return rendered[2:(end - 1)]
    end
    return rendered
end

function command_string(cmd::AbstractString)::String
    return String(cmd)
end

function build_execution_job(;
        cmd::AbstractString,
        job_name::AbstractString,
        site::Union{Symbol, AbstractString},
        time_limit::AbstractString,
        cpus_per_task::Integer = 1,
        mem_gb = nothing,
        partition = nothing,
        qos = nothing,
        account = nothing,
        mail_user = nothing,
        output_path = nothing,
        error_path = nothing,
        workdir = nothing,
        template_name = nothing
)
    return JobSpec(
        job_name = job_name,
        cmd = cmd,
        site = site,
        time_limit = time_limit,
        cpus_per_task = Int(cpus_per_task),
        mem_gb = mem_gb,
        partition = partition,
        qos = qos,
        account = account,
        mail_user = mail_user,
        output_path = output_path,
        error_path = error_path,
        workdir = workdir,
        template_name = template_name
    )
end

function _local_pipeline(job::JobSpec, cmd::Cmd)
    out_path = job.output_path
    err_path = job.error_path

    if out_path === nothing && err_path === nothing
        return cmd
    end
    if out_path !== nothing
        mkpath(dirname(out_path))
    end
    if err_path !== nothing
        mkpath(dirname(err_path))
    end
    if out_path === nothing
        return pipeline(cmd; stderr = err_path)
    elseif err_path === nothing
        return pipeline(cmd; stdout = out_path)
    end
    return pipeline(cmd; stdout = out_path, stderr = err_path)
end

function execute(job::JobSpec, ::LocalExecutor)
    shell_cmd = `bash -lc $(job.cmd)`
    if job.workdir === nothing
        run(_local_pipeline(job, shell_cmd))
    else
        mkpath(job.workdir)
        cd(job.workdir) do
            run(_local_pipeline(job, shell_cmd))
        end
    end
    return true
end

function execute(job::JobSpec, executor::SlurmExecutor)::SubmitResult
    return submit(job; dry_run = executor.dry_run)
end

function _override_if_set(new_value, current_value)
    return new_value === nothing ? current_value : new_value
end

function _job_with_backend(job::JobSpec, backend::SLURMBackend)::JobSpec
    return JobSpec(
        job_name = job.job_name,
        cmd = job.cmd,
        site = backend.site,
        time_limit = backend.time_limit,
        workdir = job.workdir,
        env = job.env,
        modules = job.modules,
        nodes = job.nodes,
        ntasks = job.ntasks,
        ntasks_per_node = job.ntasks_per_node,
        cpus_per_task = backend.cpus_per_task,
        mem_gb = backend.mem_gb,
        mem_per_cpu_mb = job.mem_per_cpu_mb,
        gpus_per_node = job.gpus_per_node,
        gpus_per_task = job.gpus_per_task,
        gpu_bind = job.gpu_bind,
        partition = _override_if_set(backend.partition, job.partition),
        qos = _override_if_set(backend.qos, job.qos),
        constraint = job.constraint,
        account = _override_if_set(backend.account, job.account),
        output_path = job.output_path,
        error_path = job.error_path,
        mail_user = _override_if_set(backend.mail_user, job.mail_user),
        mail_type = job.mail_type,
        requeue = job.requeue,
        container_image = job.container_image,
        container_args = job.container_args,
        container_mounts = job.container_mounts,
        container_workdir = job.container_workdir,
        container_engine = job.container_engine,
        template_name = job.template_name
    )
end

function execute(job::JobSpec, executor::ConfiguredSlurmExecutor)::SubmitResult
    configured_job = _job_with_backend(job, executor.backend)
    return submit(configured_job; dry_run = executor.backend.dry_run)
end

function execute(job::JobSpec, executor::CollectExecutor)::Int
    push!(executor.jobs, job)
    return length(executor.jobs)
end

function execute(job::JobSpec, executor::DryRunExecutor)::SubmitResult
    push!(executor.jobs, job)
    result = submit(job; dry_run = true)
    push!(executor.results, result)
    return result
end
