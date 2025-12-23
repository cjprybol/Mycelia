"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to SLURM scheduler on Lawrence Berkeley Lab's Lawrencium cluster.

# Arguments
- `job_name`: Name identifier for the SLURM job
- `mail_user`: Email address for job notifications
- `mail_type`: Notification type ("ALL", "BEGIN", "END", "FAIL", or "NONE")
- `logdir`: Directory for SLURM output and error logs
- `partition`: Lawrencium compute partition
- `qos`: Quality of Service level
- `account`: Project account for billing
- `nodes`: Number of nodes to allocate
- `ntasks`: Number of tasks to spawn
- `time`: Wall time limit in format "days-hours:minutes:seconds"
- `cpus_per_task`: CPU cores per task
- `mem_gb`: Memory per node in GB
- `cmd`: Shell command to execute

# Returns
- `true` if submission was successful (or job index if collecting)
"""
function lawrencium_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        partition::String="lr4",
        qos::String="lr_normal",
        account::String,
        nodes::Int=1,
        ntasks::Int=1,
        time::String="3-00:00:00",
        cpus_per_task::Int=16,
        mem_gb::Int=64,
        cmd::String
    )
    
    extra = Dict{Symbol, Any}()
    extra[:mail_user] = mail_user
    extra[:mail_type] = mail_type
    extra[:nodes] = nodes
    extra[:ntasks] = ntasks
    extra[:qos] = qos
    
    # Standardize log outputs
    stdout = joinpath(logdir, "%j.%x.out")
    stderr = joinpath(logdir, "%j.%x.err")

    job = Execution.JobSpec(
        cmd;
        name=job_name,
        time=time,
        cpus=cpus_per_task,
        mem="$(mem_gb)G",
        partition=partition,
        account=account,
        stdout=stdout,
        stderr=stderr,
        extra=extra
    )
    
    # We want to force Slurm execution here to match previous behavior if the user calls this directly,
    # BUT the point of the refactor is to allow `submit(job)` to decide.
    # Existing calls expect immediate submission to Slurm.
    # The `Mycelia.Execution.submit` function dispatches based on `CURRENT_EXECUTOR`.
    # If `CURRENT_EXECUTOR` is `LocalExecutor` (default), calling `submit(job)` would run it locally!
    # That is NOT what `lawrencium_sbatch` is supposed to do—it's explicitly a SLURM submitter.
    #
    # However, to allow "collecting" these jobs via `CollectExecutor`, we should use `submit(job)`.
    #
    # So:
    # 1. If we are in `CollectExecutor` mode, `submit(job)` will collect it.
    # 2. If we are in `SlurmExecutor` mode, `submit(job)` will submit it to Slurm.
    # 3. If we are in `LocalExecutor` mode (default), `submit(job)` runs locally.
    #    Wait, `lawrencium_sbatch` assumes we WANT to submit to Slurm.
    #
    # Solution: We should temporarily switch to `SlurmExecutor` if the current executor is `LocalExecutor`,
    # because this function implies an intent to use Slurm.
    # But if the current executor is `CollectExecutor`, we honor it.
    
    exec = Execution.CURRENT_EXECUTOR[]
    if exec isa Execution.LocalExecutor
        # The user called `lawrencium_sbatch` directly in a normal session.
        # They expect it to submit to Slurm, not run locally.
        # So we force SlurmExecutor.
        return Execution.with_executor(Execution.SlurmExecutor(false)) do
             Execution.submit(job)
        end
    else
        # Allow CollectExecutor or pre-configured SlurmExecutor to handle it
        return Execution.submit(job)
    end
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to SLURM using sbatch with specified parameters.

# Arguments
- `job_name::String`: Name identifier for the SLURM job
- `mail_user::String`: Email address for job notifications
- `mail_type::String`: Type of mail notifications (default: "ALL")
- `logdir::String`: Directory for error and output logs
- `partition::String`: SLURM partition to submit job to
- `account::String`: Account to charge for compute resources
- `nodes::Int`: Number of nodes to allocate
- `ntasks::Int`: Number of tasks to run
- `time::String`: Maximum wall time in format "days-hours:minutes:seconds"
- `cpus_per_task::Int`: CPUs per task
- `mem_gb::Int`: Memory in GB
- `cmd::String`: Command to execute

# Returns
- `Bool`: Returns true if submission succeeds
"""
function scg_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=mkpath("$(homedir())/workspace/slurmlogs"),
        partition::String="nih_s10",
        account::String,
        nodes::Int=1,
        ntasks::Int=1,
        time::String="7-00:00:00",
        cpus_per_task::Int=12,
        mem_gb::Int=96,
        cmd::String
    )
    
    extra = Dict{Symbol, Any}()
    extra[:mail_user] = mail_user
    extra[:mail_type] = mail_type
    extra[:nodes] = nodes
    extra[:ntasks] = ntasks
    
    stdout = joinpath(logdir, "%j.%x.out")
    stderr = joinpath(logdir, "%j.%x.err")

    job = Execution.JobSpec(
        cmd;
        name=job_name,
        time=time,
        cpus=cpus_per_task,
        mem="$(mem_gb)G",
        partition=partition,
        account=account,
        stdout=stdout,
        stderr=stderr,
        extra=extra
    )
    
    exec = Execution.CURRENT_EXECUTOR[]
    if exec isa Execution.LocalExecutor
        return Execution.with_executor(Execution.SlurmExecutor(false)) do
             Execution.submit(job)
        end
    else
        return Execution.submit(job)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to NERSC's SLURM scheduler using the shared QOS.

# Arguments
- `job_name`: Identifier for the job
- `mail_user`: Email address for job notifications
- `mail_type`: Notification type
- `logdir`: Directory for storing job output and error logs
- `qos`: Quality of Service level
- `nodes`: Number of nodes to allocate
- `ntasks`: Number of tasks to run
- `time`: Maximum wall time in format "days-hours:minutes:seconds"
- `cpus_per_task`: Number of CPUs per task
- `mem_gb`: Memory per node in GB
- `cmd`: Command to execute
- `constraint`: Node type constraint

# Returns
`true` if job submission succeeds
"""
function nersc_sbatch_shared(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=mkpath("$(homedir())/workspace/slurmlogs"),
        qos::String="shared",
        nodes::Int=1,
        ntasks::Int=1,
        time::String="2-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 2,
        cmd::String,
        constraint::String="cpu"
    )
    
    extra = Dict{Symbol, Any}()
    extra[:mail_user] = mail_user
    extra[:mail_type] = mail_type
    extra[:nodes] = nodes
    extra[:ntasks] = ntasks
    extra[:qos] = qos
    extra[:constraint] = constraint
    
    stdout = joinpath(logdir, "%j.%x.out")
    stderr = joinpath(logdir, "%j.%x.err")

    job = Execution.JobSpec(
        cmd;
        name=job_name,
        time=time,
        cpus=cpus_per_task,
        mem="$(mem_gb)G",
        stdout=stdout,
        stderr=stderr,
        extra=extra
    )
    
    exec = Execution.CURRENT_EXECUTOR[]
    if exec isa Execution.LocalExecutor
        return Execution.with_executor(Execution.SlurmExecutor(false)) do
             Execution.submit(job)
        end
    else
        return Execution.submit(job)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a batch job to NERSC's SLURM workload manager.

# Arguments
- `job_name`: Identifier for the SLURM job
- `mail_user`: Email address for job notifications
- `mail_type`: Notification type
- `logdir`: Directory for storing job output/error logs
- `scriptdir`: Directory for storing generated SLURM scripts
- `qos`: Quality of Service level
- `nodes`: Number of nodes to allocate
- `ntasks`: Number of tasks to run
- `time`: Maximum wall time in format "days-HH:MM:SS"
- `cpus_per_task`: CPU cores per task
- `mem_gb`: Memory per node in GB
- `cmd`: Command(s) to execute (String or Vector{String})
- `constraint`: Node type constraint

# Returns
- `true` if job submission succeeds
"""
function nersc_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        scriptdir::String=mkpath(joinpath(homedir(), "workspace/slurm")),
        qos::String="regular",
        nodes::Int=1,
        ntasks::Int=1,
        time::String="2-00:00:00",
        cpus_per_task::Int=Mycelia.NERSC_CPU,
        mem_gb::Int=Mycelia.NERSC_MEM,
        cmd::Union{String, Vector{String}},
        constraint::String="cpu"
    )
    
    # Process commands
    cmd_block = if isa(cmd, String)
        cmd  # Single command as is
    else
        join(cmd, "\n")  # Multiple commands joined with newlines
    end
    
    extra = Dict{Symbol, Any}()
    extra[:mail_user] = mail_user
    extra[:mail_type] = mail_type
    extra[:nodes] = nodes
    extra[:ntasks] = ntasks
    extra[:qos] = qos
    extra[:constraint] = constraint
    
    stdout = joinpath(logdir, "%j.%x.out")
    stderr = joinpath(logdir, "%j.%x.err")

    job = Execution.JobSpec(
        cmd_block;
        name=job_name,
        time=time,
        cpus=cpus_per_task,
        mem="$(mem_gb)G",
        stdout=stdout,
        stderr=stderr,
        extra=extra
    )
    
    exec = Execution.CURRENT_EXECUTOR[]
    if exec isa Execution.LocalExecutor
        # NERSC sbatch previously created a script file. `SlurmExecutor(true)` creates a script but pipes it.
        # If we really want to preserve the file creation, we'd need to customize or allow SlurmExecutor to write to a path.
        # For now, let's stick to the consistent pipeline unless file artifacts are critical.
        # The previous implementation wrote `script_name` based on timestamp.
        return Execution.with_executor(Execution.SlurmExecutor(false)) do
             Execution.submit(job)
        end
    else
        return Execution.submit(job)
    end
end
# function nersc_sbatch_regular(;
#         job_name::String,
#         mail_user::String,
#         mail_type::String="ALL",
#         logdir::String=mkpath("$(homedir())/workspace/slurmlogs"),
#         qos::String="regular",
#         nodes::Int=1,
#         ntasks::Int=1,
#         time::String="2-00:00:00",
#         cpus_per_task::Int=Mycelia.NERSC_CPU,
#         mem_gb::Int=Mycelia.NERSC_MEM,
#         cmd::String,
#         constraint::String="cpu"
#     )
#     submission = 
#     `sbatch
#     --job-name=$(job_name)
#     --mail-user=$(mail_user)
#     --mail-type=$(mail_type)
#     --error=$(logdir)/%j.%x.err
#     --output=$(logdir)/%j.%x.out
#     --qos=$(qos)
#     --nodes=$(nodes)
#     --ntasks=$(ntasks)
#     --time=$(time)   
#     --cpus-per-task=$(cpus_per_task)
#     --mem=$(mem_gb)G
#     --constraint=cpu
#     --wrap $(cmd)
#     `
#     sleep(5)
#     run(submission)
#     sleep(5)
#     return true
# end
