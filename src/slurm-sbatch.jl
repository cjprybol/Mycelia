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
- `true` if submission was successful

# Note
Function includes 5-second delays before and after submission for queue stability.
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
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --partition=$(partition)
    --qos=$(qos)
    --account=$(account)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to SLURM using sbatch with specified parameters.

# Arguments
- `job_name::String`: Name identifier for the SLURM job
- `mail_user::String`: Email address for job notifications
- `mail_type::String`: Type of mail notifications (default: "ALL")
- `logdir::String`: Directory for error and output logs (default: "~/workspace/slurmlogs")
- `partition::String`: SLURM partition to submit job to
- `account::String`: Account to charge for compute resources
- `nodes::Int`: Number of nodes to allocate (default: 1)
- `ntasks::Int`: Number of tasks to run (default: 1)
- `time::String`: Maximum wall time in format "days-hours:minutes:seconds" (default: "1-00:00:00")
- `cpus_per_task::Int`: CPUs per task (default: 12)
- `mem_gb::Int`: Memory in GB, defaults to 96GB
- `cmd::String`: Command to execute

# Returns
- `Bool`: Returns true if submission succeeds

# Notes
- Function includes 5-second delays before and after submission
- Memory is automatically scaled with CPU count
- Log files are named with job ID (%j) and job name (%x)
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
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --partition=$(partition)
    --account=$(account)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to NERSC's SLURM scheduler using the shared QOS (Quality of Service).

# Arguments
- `job_name`: Identifier for the job
- `mail_user`: Email address for job notifications
- `mail_type`: Notification type ("ALL", "BEGIN", "END", "FAIL", "REQUEUE", "STAGE_OUT")
- `logdir`: Directory for storing job output and error logs
- `qos`: Quality of Service level ("shared", "regular", "preempt", "premium")
- `nodes`: Number of nodes to allocate
- `ntasks`: Number of tasks to run
- `time`: Maximum wall time in format "days-hours:minutes:seconds"
- `cpus_per_task`: Number of CPUs per task
- `mem_gb`: Memory per node in GB (default: 2GB per CPU)
- `cmd`: Command to execute
- `constraint`: Node type constraint ("cpu" or "gpu")

# Resource Limits
- Maximum memory per node: 512GB
- Maximum cores per node: 128
- Default memory allocation: 2GB per CPU requested

# QOS Options
- shared: Default QOS for shared node usage
- regular: Standard priority
- preempt: Reduced credit usage but preemptible
- premium: 5x throughput priority (limited usage)

# Returns
`true` if job submission succeeds


https://docs.nersc.gov/jobs/policy/
https://docs.nersc.gov/systems/perlmutter/architecture/#cpu-nodes

default is to use shared qos

use
- regular
- preempt (reduced credit usage but not guaranteed to finish)
- premium (priority runs limited to 5x throughput)

max request is 512Gb memory and 128 cores per node

https://docs.nersc.gov/systems/perlmutter/running-jobs/#tips-and-tricks
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
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --qos=$(qos)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --constraint=cpu
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a batch job to NERSC's SLURM workload manager.

# Arguments
- `job_name`: Identifier for the SLURM job
- `mail_user`: Email address for job notifications
- `mail_type`: Notification type ("ALL", "BEGIN", "END", "FAIL", or "NONE")
- `logdir`: Directory for storing job output/error logs
- `scriptdir`: Directory for storing generated SLURM scripts
- `qos`: Quality of Service level ("regular", "premium", or "preempt")
- `nodes`: Number of nodes to allocate
- `ntasks`: Number of tasks to run
- `time`: Maximum wall time in format "days-HH:MM:SS"
- `cpus_per_task`: CPU cores per task
- `mem_gb`: Memory per node in GB
- `cmd`: Command(s) to execute (String or Vector{String})
- `constraint`: Node type constraint ("cpu" or "gpu")

# Returns
- `true` if job submission succeeds
- `false` if submission fails

# QoS Options
- regular: Standard priority queue
- premium: High priority queue (5x throughput limit)
- preempt: Reduced credit usage but jobs may be interrupted


https://docs.nersc.gov/jobs/policy/
https://docs.nersc.gov/systems/perlmutter/architecture/#cpu-nodes

default is to use shared qos

use
- regular
- preempt (reduced credit usage but not guaranteed to finish)
- premium (priorty runs limited to 5x throughput)

https://docs.nersc.gov/systems/perlmutter/running-jobs/#tips-and-tricks
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
        
    # Generate timestamp
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
    
    # Create script filename
    script_name = "$(timestamp)-$(job_name).sh"
    script_path = joinpath(scriptdir, script_name)
    
    # Process commands
    cmd_block = if isa(cmd, String)
        cmd  # Single command as is
    else
        join(cmd, "\n")  # Multiple commands joined with newlines
    end
    
    # Create script content
    script_content = """
    #!/bin/bash
    #SBATCH --job-name=$(job_name)
    #SBATCH --mail-user=$(mail_user)
    #SBATCH --mail-type=$(mail_type)
    #SBATCH --error=$(logdir)/%j.%x.err
    #SBATCH --output=$(logdir)/%j.%x.out
    #SBATCH --qos=$(qos)
    #SBATCH --nodes=$(nodes)
    #SBATCH --ntasks=$(ntasks)
    #SBATCH --time=$(time)
    #SBATCH --cpus-per-task=$(cpus_per_task)
    #SBATCH --mem=$(mem_gb)G
    #SBATCH --constraint=$(constraint)

    $cmd_block
    """
    
    # Write script to file
    write(script_path, script_content)
    
    # Make script executable
    chmod(script_path, 0o755)
    
    # Submit the job
    sleep(5)
    try
        run(`sbatch $script_path`)
    catch e
        @error "Failed to submit job with sbatch: $e"
        return false
    end
    sleep(5)
    
    return true
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