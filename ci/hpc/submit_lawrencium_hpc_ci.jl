#!/usr/bin/env julia

import Mycelia

function parse_int_env(key::String, default::Int)
    val = get(ENV, key, "")
    if isempty(val)
        return default
    end
    try
        return parse(Int, val)
    catch
        error("Invalid integer for $(key): $(val)")
    end
end

function required_env(key::String)
    val = get(ENV, key, "")
    if isempty(val)
        error("Missing required environment variable: $(key)")
    end
    return val
end

repo_root = normpath(joinpath(@__DIR__, "..", ".."))
default_cmd = "bash " * joinpath(repo_root, "ci", "hpc", "run_hpc_ci.sh")
cmd = get(ENV, "HPC_CMD", default_cmd)

Mycelia.lawrencium_sbatch(
    job_name = get(ENV, "HPC_JOB_NAME", "mycelia-hpc-ci"),
    mail_user = required_env("HPC_MAIL_USER"),
    mail_type = get(ENV, "HPC_MAIL_TYPE", "END,FAIL"),
    logdir = get(ENV, "HPC_LOGDIR", joinpath(homedir(), "workspace", "slurmlogs")),
    partition = get(ENV, "HPC_PARTITION", "lr4"),
    qos = get(ENV, "HPC_QOS", "lr_normal"),
    account = required_env("HPC_ACCOUNT"),
    nodes = parse_int_env("HPC_NODES", 1),
    ntasks = parse_int_env("HPC_NTASKS", 1),
    time = get(ENV, "HPC_TIME", "4-00:00:00"),
    cpus_per_task = parse_int_env("HPC_CPUS_PER_TASK", 8),
    mem_gb = parse_int_env("HPC_MEM_GB", 32),
    cmd = cmd,
)
