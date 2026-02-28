"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to SLURM scheduler on Lawrence Berkeley Lab's Lawrencium cluster.

This function preserves the historical public API and now routes through
`JobSpec` + template rendering.
"""
function lawrencium_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String = "ALL",
        logdir::String = mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        partition::String = "lr4",
        qos::String = "lr_normal",
        account::String,
        nodes::Int = 1,
        ntasks::Int = 1,
        time::String = "3-00:00:00",
        cpus_per_task::Int = 16,
        mem_gb::Int = 64,
        cmd::String,
        dry_run::Bool = false,
        executor = nothing
)
    job = JobSpec(
        job_name = job_name,
        cmd = cmd,
        site = :lawrencium,
        time_limit = time,
        partition = partition,
        qos = qos,
        account = account,
        nodes = nodes,
        ntasks = ntasks,
        cpus_per_task = cpus_per_task,
        mem_gb = mem_gb,
        output_path = joinpath(logdir, "%j.%x.out"),
        error_path = joinpath(logdir, "%j.%x.err"),
        mail_user = mail_user,
        mail_type = mail_type
    )

    if executor !== nothing
        resolved_executor = resolve_executor(executor)
        if dry_run && resolved_executor isa SlurmExecutor
            resolved_executor = SlurmExecutor(dry_run = true)
        end
        return execute(job, resolved_executor)
    end

    sleep(5)
    result = submit(job; dry_run = dry_run)
    sleep(5)

    if !result.ok
        @error "Lawrencium submission failed" errors = result.errors warnings = result.warnings
        return false
    end

    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to Stanford SCG via SLURM.

This function preserves the historical public API and now routes through
`JobSpec` + template rendering.
"""
function scg_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String = "ALL",
        logdir::String = mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        partition::String = "nih_s10",
        account::String,
        nodes::Int = 1,
        ntasks::Int = 1,
        time::String = "7-00:00:00",
        cpus_per_task::Int = 12,
        mem_gb::Int = 96,
        cmd::String,
        dry_run::Bool = false,
        executor = nothing
)
    job = JobSpec(
        job_name = job_name,
        cmd = cmd,
        site = :scg,
        time_limit = time,
        partition = partition,
        account = account,
        nodes = nodes,
        ntasks = ntasks,
        cpus_per_task = cpus_per_task,
        mem_gb = mem_gb,
        output_path = joinpath(logdir, "%j.%x.out"),
        error_path = joinpath(logdir, "%j.%x.err"),
        mail_user = mail_user,
        mail_type = mail_type
    )

    if executor !== nothing
        resolved_executor = resolve_executor(executor)
        if dry_run && resolved_executor isa SlurmExecutor
            resolved_executor = SlurmExecutor(dry_run = true)
        end
        return execute(job, resolved_executor)
    end

    sleep(5)
    result = submit(job; dry_run = dry_run)
    sleep(5)

    if !result.ok
        @error "SCG submission failed" errors = result.errors warnings = result.warnings
        return false
    end

    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to NERSC's shared QoS profile.

This compatibility wrapper calls `nersc_sbatch` with `qos="shared"`.
"""
function nersc_sbatch_shared(;
        job_name::String,
        mail_user::String,
        mail_type::String = "ALL",
        logdir::String = mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        qos::String = "shared",
        nodes::Int = 1,
        ntasks::Int = 1,
        time::String = "2-00:00:00",
        cpus_per_task::Int = 1,
        mem_gb::Int = cpus_per_task * 2,
        cmd::String,
        constraint::String = "cpu",
        account::Union{Nothing, String} = nothing,
        dry_run::Bool = false,
        executor = nothing
)
    resolved_account = account === nothing ? get(ENV, "NERSC_ACCOUNT", nothing) : account
    if resolved_account === nothing
        error("NERSC account is required. Pass account=... or set ENV[\"NERSC_ACCOUNT\"].")
    end

    return nersc_sbatch(
        job_name = job_name,
        mail_user = mail_user,
        mail_type = mail_type,
        logdir = logdir,
        qos = qos,
        nodes = nodes,
        ntasks = ntasks,
        time = time,
        cpus_per_task = cpus_per_task,
        mem_gb = mem_gb,
        cmd = cmd,
        constraint = constraint,
        account = resolved_account,
        dry_run = dry_run,
        executor = executor
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a batch job to NERSC's SLURM workload manager.

This function preserves historical defaults and now delegates rendering/
submission to the template system.
"""
function nersc_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String = "ALL",
        logdir::String = mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        scriptdir::String = mkpath(joinpath(homedir(), "workspace/slurm")),
        qos::String = "regular",
        nodes::Int = 1,
        ntasks::Int = 1,
        time::String = "2-00:00:00",
        cpus_per_task::Int = Mycelia.NERSC_CPU,
        mem_gb::Int = Mycelia.NERSC_MEM,
        cmd::Union{String, Vector{String}},
        constraint::String = "cpu",
        account::Union{Nothing, String} = nothing,
        gpus_per_node::Union{Nothing, Int} = nothing,
        gpus_per_task::Union{Nothing, Int} = nothing,
        gpu_bind::Union{Nothing, String} = nothing,
        requeue::Bool = false,
        dry_run::Bool = false,
        executor = nothing
)
    resolved_account = account === nothing ? get(ENV, "NERSC_ACCOUNT", nothing) : account
    if resolved_account === nothing
        error("NERSC account is required. Pass account=... or set ENV[\"NERSC_ACCOUNT\"].")
    end

    cmd_block = if cmd isa String
        cmd
    else
        join(cmd, "\n")
    end

    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
    script_name = "$(timestamp)-$(job_name).sh"
    script_path = joinpath(scriptdir, script_name)

    job = JobSpec(
        job_name = job_name,
        cmd = cmd_block,
        site = :nersc,
        time_limit = time,
        qos = qos,
        account = resolved_account,
        nodes = nodes,
        ntasks = ntasks,
        cpus_per_task = cpus_per_task,
        mem_gb = mem_gb,
        constraint = constraint,
        gpus_per_node = gpus_per_node,
        gpus_per_task = gpus_per_task,
        gpu_bind = gpu_bind,
        requeue = requeue,
        output_path = joinpath(logdir, "%j.%x.out"),
        error_path = joinpath(logdir, "%j.%x.err"),
        mail_user = mail_user,
        mail_type = mail_type
    )

    if executor !== nothing
        resolved_executor = resolve_executor(executor)
        if dry_run && resolved_executor isa SlurmExecutor
            resolved_executor = SlurmExecutor(dry_run = true)
        end
        return execute(job, resolved_executor)
    end

    sleep(5)
    result = submit(job; dry_run = dry_run, path = script_path)
    sleep(5)

    if !result.ok
        @error "NERSC submission failed" errors = result.errors warnings = result.warnings
        return false
    end

    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Execute a command directly on the Lovelace dedicated server (no SLURM).
"""
function lovelace_run(;
        cmd::String,
        logdir::String = mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        job_name::String = "lovelace-job"
)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
    log_out = joinpath(logdir, "$(timestamp).$(job_name).out")
    log_err = joinpath(logdir, "$(timestamp).$(job_name).err")

    try
        run(pipeline(`bash -c $cmd`; stdout = log_out, stderr = log_err))
    catch e
        @error "Lovelace job '$(job_name)' failed: $e"
        return false
    end
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Dispatch a job to the appropriate compute backend based on `SITE_TAG`.

Routes to:
- `lawrencium` -> `lawrencium_sbatch`
- `lovelace` -> `lovelace_run`
- `nersc` -> `nersc_sbatch`
- `stanford_scg` / `scg` -> `scg_sbatch`
"""
function submit_job(; site::String = get(ENV, "SITE_TAG", "default"), kwargs...)
    normalized_site = lowercase(site)
    if normalized_site == "lawrencium"
        return lawrencium_sbatch(; kwargs...)
    elseif normalized_site == "lovelace"
        return lovelace_run(; kwargs...)
    elseif normalized_site == "nersc"
        return nersc_sbatch(; kwargs...)
    elseif normalized_site == "stanford_scg" || normalized_site == "scg"
        return scg_sbatch(; kwargs...)
    else
        error("Unknown SITE_TAG: $site. Expected one of: lawrencium, lovelace, nersc, stanford_scg, scg")
    end
end
