module Execution

import Dates

# Remove exports - use fully qualified names only

"""
    JobSpec

A specification for a job to be executed.
"""
struct JobSpec
    cmd::Union{Cmd, String}
    name::String
    time::String        # "HH:MM:SS" or "D-HH:MM:SS"
    cpus::Int
    mem::String         # "8G"
    workdir::String
    stdout::Union{Nothing, String}
    stderr::Union{Nothing, String}
    partition::Union{Nothing, String}
    account::Union{Nothing, String}
    dependency::Union{Nothing, String}
    extra::Dict{Symbol, Any}
end

"""
    JobSpec(cmd; kwargs...)

Construct a JobSpec with default values.
"""
function JobSpec(cmd::Union{Cmd, String};
    name::String="job",
    time::String="01:00:00",
    cpus::Int=1,
    mem::String="4G",
    workdir::String=pwd(),
    stdout::Union{Nothing, String}=nothing,
    stderr::Union{Nothing, String}=nothing,
    partition::Union{Nothing, String}=nothing,
    account::Union{Nothing, String}=nothing,
    dependency::Union{Nothing, String}=nothing,
    extra::Dict{Symbol, Any}=Dict{Symbol, Any}()
)
    return JobSpec(cmd, name, time, cpus, mem, workdir, stdout, stderr, partition, account, dependency, extra)
end

"""
    SlurmTemplate

A template for default Slurm job parameters.
"""
struct SlurmTemplate
    name::Symbol
    partition::Union{Nothing, String}
    account::Union{Nothing, String}
    mem::Union{Nothing, String}
    cpus::Union{Nothing, Int}
    extra::Dict{Symbol, Any}
end

# Define standard templates registry
const SLURM_TEMPLATES = Dict{Symbol, SlurmTemplate}(
    :nersc => SlurmTemplate(:nersc, "regular", nothing, "128G", 256, Dict{Symbol, Any}(:constraint => "cpu")),
    :test => SlurmTemplate(:test, "debug", nothing, "1G", 1, Dict{Symbol, Any}()),
    :lr4 => SlurmTemplate(:lr4, "lr4", nothing, "64G", 32, Dict{Symbol, Any}(:qos => "lr_normal")),
    :nih_s10 => SlurmTemplate(:nih_s10, "norm", "my_account", "128G", 32, Dict{Symbol, Any}()) 
)

"""
    register_slurm_template(name::Symbol, template::SlurmTemplate)

Register a new Slurm template.
"""
function register_slurm_template(name::Symbol, template::SlurmTemplate)
    SLURM_TEMPLATES[name] = template
end

"""
    resolve_job_resources(requested_mem, requested_cpus, template_name, default_mem, default_cpus)

Resolve memory and CPUs strings/ints based on priority:
1. Explicit requested value (if not nothing)
2. Template value (if template exists and has value)
3. Default value
"""
function resolve_job_resources(
    mem::Union{Nothing, String}, 
    cpus::Union{Nothing, Int}, 
    template_name::Union{Nothing, Symbol}, 
    default_mem::String, 
    default_cpus::Int
)
    template = (!isnothing(template_name) && haskey(SLURM_TEMPLATES, template_name)) ? SLURM_TEMPLATES[template_name] : nothing
    
    final_mem = if !isnothing(mem)
        mem
    elseif !isnothing(template) && !isnothing(template.mem)
        template.mem
    else
        default_mem
    end
    
    final_cpus = if !isnothing(cpus)
        cpus
    elseif !isnothing(template) && !isnothing(template.cpus)
        template.cpus
    else
        default_cpus
    end
    
    return final_mem, final_cpus
end

abstract type AbstractExecutor end

"""
    LocalExecutor

Executes jobs locally using `run()`.
"""
struct LocalExecutor <: AbstractExecutor end

"""
    SlurmExecutor

Submits jobs to Slurm using `sbatch`.
"""
struct SlurmExecutor <: AbstractExecutor
    submit_script::Bool # If true, writes a script file. If false, uses --wrap (or pipe stdin).
    template::Union{Nothing, Symbol}
end

SlurmExecutor(;submit_script=false, template=nothing) = SlurmExecutor(submit_script, template)

"""
    CollectExecutor

Collects jobs into a vector for later processing.
"""
struct CollectExecutor <: AbstractExecutor
    jobs::Vector{JobSpec}
end

CollectExecutor() = CollectExecutor(JobSpec[])

"""
    DryRunExecutor

Collects jobs without executing them (alias for CollectExecutor).
"""
struct DryRunExecutor <: AbstractExecutor
    jobs::Vector{JobSpec}
end

DryRunExecutor() = DryRunExecutor(JobSpec[])


# Global executor state
const CURRENT_EXECUTOR = Ref{AbstractExecutor}(LocalExecutor())

"""
    with_executor(f::Function, executor::AbstractExecutor)

Execute the function `f` with the specified `executor` as the current executor.
"""
function with_executor(f::Function, executor::AbstractExecutor)
    old = CURRENT_EXECUTOR[]
    CURRENT_EXECUTOR[] = executor
    try
        return f()
    finally
        CURRENT_EXECUTOR[] = old
    end
end

"""
    submit(job::JobSpec)

Submit a job to the current executor.
"""
function submit(job::JobSpec)
    exec = CURRENT_EXECUTOR[]
    if exec isa LocalExecutor
        return submit_local(job)
    elseif exec isa SlurmExecutor
        # Apply template overrides for partition/account if present
        # Note: mem/cpus are already resolved at JobSpec creation time by the wrapper
        effective_job = apply_slurm_template(job, exec.template)
        return submit_slurm(exec, effective_job)
    elseif (exec isa CollectExecutor) || (exec isa DryRunExecutor)
        push!(exec.jobs, job)
        return length(exec.jobs)
    else
        error("Unknown executor $(typeof(exec))")
    end
end

function apply_slurm_template(job::JobSpec, template_name::Union{Nothing, Symbol})
    if isnothing(template_name) || !haskey(SLURM_TEMPLATES, template_name)
        return job
    end
    
    template = SLURM_TEMPLATES[template_name]
    
    # Only override partition/account/extras if they are missing in the job
    # or if we want template to have precedence?
    # Usually explicit job args > template > defaults.
    # JobSpec doesn't have "unset" for partition/account (they are Union{Nothing,String}).
    
    new_partition = isnothing(job.partition) ? template.partition : job.partition
    new_account = isnothing(job.account) ? template.account : job.account
    new_extra = merge(template.extra, job.extra)
    
    return JobSpec(
        job.cmd, job.name, job.time, job.cpus, job.mem, job.workdir,
        job.stdout, job.stderr, new_partition, new_account, job.dependency, new_extra
    )
end

function submit_local(job::JobSpec)
    cmd_obj = job.cmd isa String ? `bash -c $(job.cmd)` : job.cmd

    if job.workdir != pwd()
        cd(job.workdir) do
            run(cmd_obj)
        end
    else
        run(cmd_obj)
    end
    return true
end

function submit_slurm(exec::SlurmExecutor, job::JobSpec)
    # Generate script and pipe to sbatch
    script = job_to_slurm_script(job)
    return readchomp(pipeline(`sbatch --parsable`, stdin=IOBuffer(script)))
end

"""
    generate_sbatch_script(job)

Create an sbatch script string for a given job.
"""
function generate_sbatch_script(job::JobSpec)
    io = IOBuffer()
    println(io, "#!/bin/bash")
    println(io, "#SBATCH --job-name=$(job.name)")
    println(io, "#SBATCH --time=$(job.time)")
    println(io, "#SBATCH --cpus-per-task=$(job.cpus)")
    println(io, "#SBATCH --mem=$(job.mem)")

    if !isnothing(job.partition)
        println(io, "#SBATCH --partition=$(job.partition)")
    end
    if !isnothing(job.account)
        println(io, "#SBATCH --account=$(job.account)")
    end
    if !isnothing(job.dependency)
        println(io, "#SBATCH --dependency=$(job.dependency)")
    end
    if !isnothing(job.stdout)
        println(io, "#SBATCH --output=$(job.stdout)")
    end
    if !isnothing(job.stderr)
        println(io, "#SBATCH --error=$(job.stderr)")
    end

    for (k, v) in job.extra
        println(io, "#SBATCH --$(k)=$(v)")
    end

    println(io, "")
    println(io, "cd $(job.workdir)")
    println(io, "")

    cmd_str = if job.cmd isa Cmd
        s = string(job.cmd)
        (startswith(s, "`") && endswith(s, "`")) ? s[2:end-1] : s
    else
        job.cmd
    end

    println(io, cmd_str)

    return String(take!(io))
end

function job_to_slurm_script(job::JobSpec)
    return generate_sbatch_script(job)
end


"""
    resolve_executor(x)

Resolve an argument to an executor instance.
"""
function resolve_executor(x::AbstractExecutor)
    return x
end

function resolve_executor(x::Symbol)
    if x === :local
        return LocalExecutor()
    elseif x === :slurm
        return SlurmExecutor()
    elseif x === :dryrun || x === :collect
        return DryRunExecutor()
    else
        return LocalExecutor()
    end
end

function resolve_executor(x)
    return LocalExecutor()
end

end # module
