import Test
import Mycelia

# td-hhkd: static lint of templates/slurm and of render_sbatch output.
# Reads scripts only. Does not call sbatch or any scheduler.

const SBATCH_PATTERN_LINT = normpath(joinpath(
    @__DIR__, "..", "..", "scripts", "check-sbatch-pattern.sh"))
const SBATCH_TEMPLATE_ROOT = normpath(joinpath(
    @__DIR__, "..", "..", "templates", "slurm"))

function _sbatch_lint(args::AbstractString...)
    stdout_buf = IOBuffer()
    stderr_buf = IOBuffer()
    cmd = pipeline(`bash $SBATCH_PATTERN_LINT $args`; stdout = stdout_buf, stderr = stderr_buf)
    proc = run(cmd; wait = false)
    wait(proc)
    return (
        code = proc.exitcode,
        stdout = String(take!(stdout_buf)),
        stderr = String(take!(stderr_buf)),
    )
end

function _count_sbatch(path::AbstractString)
    if isfile(path)
        return endswith(path, ".sbatch") ? 1 : 0
    end
    n = 0
    for (_, _, files) in walkdir(path)
        for file in files
            if endswith(file, ".sbatch")
                n += 1
            end
        end
    end
    return n
end

function _require_clean(result, label::AbstractString, expected_count::Integer)
    ok_line = "OK: $(expected_count) sbatch scripts clean"
    passed = result.code == 0 && occursin(ok_line, result.stdout)
    if !passed
        println(stderr, "sbatch lint failed: ", label)
        println(stderr, result.stdout)
        println(stderr, result.stderr)
    end
    Test.@test passed
    return nothing
end

function _log_paths(; output_path::AbstractString, error_path::AbstractString)
    return (; output_path, error_path)
end

Test.@testset "sbatch pattern lint (td-hhkd)" begin
    Test.@test isfile(SBATCH_PATTERN_LINT)
    # The lint is a reader. A shell command named sbatch would submit a job.
    # Filename globs such as *.sbatch are not invocations.
    scheduler_invocation = String[]
    for line in eachline(SBATCH_PATTERN_LINT)
        stripped = lstrip(line)
        startswith(stripped, "#") && continue
        if occursin(r"(^|[;&|`(])\s*sbatch(\s|$)", stripped)
            push!(scheduler_invocation, stripped)
        end
    end
    Test.@test isempty(scheduler_invocation)

    Test.@testset "--selftest" begin
        result = _sbatch_lint("--selftest")
        if result.code != 0 || !occursin("SELFTEST: all cases pass.", result.stdout)
            println(stderr, result.stdout)
            println(stderr, result.stderr)
        end
        Test.@test result.code == 0
        Test.@test occursin("SELFTEST: all cases pass.", result.stdout)
    end

    Test.@testset "templates/slurm" begin
        template_count = _count_sbatch(SBATCH_TEMPLATE_ROOT)
        Test.@test template_count > 0
        explicit = _sbatch_lint(SBATCH_TEMPLATE_ROOT)
        default = _sbatch_lint()
        _require_clean(explicit, "templates/slurm", template_count)
        _require_clean(default, "default target", template_count)
        Test.@test explicit.stdout == default.stdout
    end

    Test.@testset "renderer output" begin
        logs = _log_paths(
            output_path = "slurmlogs/%j.%x.out",
            error_path = "slurmlogs/%j.%x.err",
        )
        nersc_julia = Mycelia.JobSpec(;
            job_name = "lint-nersc-julia",
            cmd = "julia --project=. experiments/run.jl",
            site = :nersc,
            account = "m1234",
            qos = "shared",
            constraint = "cpu",
            ntasks = 1,
            cpus_per_task = 16,
            mem_gb = 64,
            time_limit = "06:00:00",
            output_path = logs.output_path,
            error_path = logs.error_path,
        )
        nersc_python = Mycelia.JobSpec(;
            job_name = "lint-nersc-python",
            cmd = "python train.py",
            site = :nersc,
            account = "m1234",
            qos = "regular",
            constraint = "cpu",
            ntasks = 1,
            cpus_per_task = 4,
            mem_gb = 8,
            time_limit = "04:00:00",
            output_path = logs.output_path,
            error_path = logs.error_path,
        )
        lrc_julia = Mycelia.JobSpec(;
            job_name = "lint-lrc-julia",
            cmd = "julia run.jl",
            site = :lawrencium,
            partition = "lr6",
            qos = "lr_normal",
            account = "pc_example",
            ntasks = 1,
            cpus_per_task = 8,
            mem_gb = 16,
            time_limit = "04:00:00",
            output_path = logs.output_path,
            error_path = logs.error_path,
        )
        scg_julia = Mycelia.JobSpec(;
            job_name = "lint-scg-julia",
            cmd = "julia --project=. workflow.jl",
            site = :scg,
            partition = "nih_s10",
            account = "PI_EXAMPLE",
            time_limit = "24:00:00",
            nodes = 1,
            ntasks = 1,
            cpus_per_task = 4,
            mem_gb = 8,
            output_path = logs.output_path,
            error_path = logs.error_path,
        )

        mktempdir() do dir
            rendered = Dict{String, String}(
                joinpath("nersc", "julia.sbatch") => Mycelia.render_sbatch(nersc_julia),
                joinpath("nersc", "python.sbatch") => Mycelia.render_sbatch(nersc_python),
                joinpath("lawrencium", "julia.sbatch") => Mycelia.render_sbatch(lrc_julia),
                joinpath("scg", "julia.sbatch") => Mycelia.render_sbatch(scg_julia),
            )
            for (relative, text) in rendered
                Test.@test !occursin("{{", text)
                path = joinpath(dir, relative)
                mkpath(dirname(path))
                write(path, text)
            end
            # Filename carries no _lrc / lawrencium cue; the site-profile line must.
            lrc_text = rendered[joinpath("lawrencium", "julia.sbatch")]
            Test.@test occursin("# Site profile: lawrencium", lrc_text)
            Test.@test occursin("JULIA_NUM_THREADS", lrc_text)
            Test.@test !occursin("JULIA_NUM_THREADS=auto", lrc_text)
            write(joinpath(dir, "profile-only.sbatch"), lrc_text)

            result = _sbatch_lint(dir)
            _require_clean(result, "renderer output", length(rendered) + 1)

            bad_dir = mktempdir()
            write(
                joinpath(bad_dir, "bad.sbatch"),
                rendered[joinpath("nersc", "julia.sbatch")] * "\njulia -u run.jl\n",
            )
            bad = _sbatch_lint(bad_dir)
            Test.@test bad.code == 1
            Test.@test occursin("julia -u", bad.stdout)
            rm(bad_dir; recursive = true, force = true)
        end
    end
end
