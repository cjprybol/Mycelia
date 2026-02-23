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
