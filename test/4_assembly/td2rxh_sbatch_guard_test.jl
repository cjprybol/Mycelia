module Td2rxhSbatchGuardTest

# Positive controls for the fail-closed guards in
# benchmarking/run_td2rxh_evidence_lrc.sbatch. Each guard is executed against the
# bad input it exists to catch, using a stub `julia` and `sha256sum` (plus failing
# `git` / `sha256sum` stubs where a test needs them) against a real temp git
# repository, so "the guard is present" is replaced by "the guard was seen to
# fire". Covered: row count, stale CSV, modified tracked file, failing rev-parse,
# failing status, failing or garbage sha256sum, missing Manifest, non-zero julia
# exit, unknown cell set, --seeds/--label/--cells passthrough, MYCELIA_PROJECT
# cleared, symlinked TMPDIR base,
# TMPDIR kept on failure, and non-local or unrecognised filesystems refused.
# Not covered: the "julia not runnable", "not in benchmarking/" and mktemp-failure
# branches, and the SIGTERM path. Round-1 review reproduced two
# of these as fail-open on the first draft (a stale CSV passing the row guard, and
# a failing git passing the dirty guard).

import Test

const SBATCH = joinpath(
    @__DIR__, "..", "..", "benchmarking", "run_td2rxh_evidence_lrc.sbatch")

const JULIA_STUB = """
#!/bin/bash
if [[ "\${2:-}" == "--version" ]]; then echo "julia version 1.10.0"; exit 0; fi
# Record the run invocation so tests can check the provenance-carrying arguments.
echo "\$*" > "\$HOME/julia_args.txt"
for v in MYCELIA_PROJECT JULIA_LOAD_PATH JULIA_PROJECT; do
  echo "\$v=\${!v-<unset>}"
done > "\$HOME/julia_env.txt"
out=""
while [[ \$# -gt 0 ]]; do
  if [[ "\$1" == "--out" ]]; then out="\$2"; fi
  shift
done
if [[ -n "\${STUB_ROWS:-}" ]]; then
  echo "header" > "\$out"
  i=0
  while [[ \$i -lt \$STUB_ROWS ]]; do echo "row" >> "\$out"; i=\$((i + 1)); done
fi
exit "\${STUB_RC:-0}"
"""

function write_exe(path::AbstractString, body::AbstractString)
    mkpath(dirname(path))
    write(path, body)
    chmod(path, 0o755)
end

"""Build a clean temp repo + stub HOME; return (repo, home, bindir)."""
function fixture()
    root = mktempdir()
    repo = joinpath(root, "repo")
    mkpath(joinpath(repo, "benchmarking"))
    write(joinpath(repo, "benchmarking", "indel_reroute_evidence.jl"), "# stub\n")
    git(args...) = run(pipeline(`git -C $repo -c user.name=t -c user.email=t@t $args`;
        stdout = devnull, stderr = devnull))
    git("init", "-q")
    git("add", "benchmarking/indel_reroute_evidence.jl")
    git("commit", "-q", "-m", "init")
    write(joinpath(repo, "Manifest.toml"), "# untracked manifest\n")
    home = joinpath(root, "home")
    write_exe(joinpath(home, ".juliaup", "bin", "julia"), JULIA_STUB)
    bindir = joinpath(root, "bin")
    # sha256sum is absent on stock macOS; a stub keeps the test host-independent.
    write_exe(joinpath(bindir, "sha256sum"), "#!/bin/bash\necho \"abcdef0123456789  \$1\"\n")
    # GNU `stat -f -c %T` (filesystem type) has no macOS equivalent; report a local
    # filesystem unless a test overrides it via STUB_FS.
    write_exe(joinpath(bindir, "stat"), "#!/bin/bash\necho \"\${STUB_FS:-ext2/ext3}\"\n")
    return repo, home, bindir
end

"""Run the sbatch script; return (exit code, combined output)."""
function run_sbatch(repo, home, bindir; extra_env = Dict{String, String}())
    env = Dict{String, String}(
        "PATH" => bindir * ":" * ENV["PATH"], "HOME" => home,
        "SLURM_SUBMIT_DIR" => joinpath(repo, "benchmarking"),
        "SLURM_ARRAY_JOB_ID" => "99", "SLURM_ARRAY_TASK_ID" => "1",
        "SLURM_NODELIST" => "n0", "SLURM_CPUS_PER_TASK" => "2",
        # A throwaway base, so failure cases that keep their TMPDIR on purpose do
        # not leave directories in the host's /tmp.
        "TD2RXH_TMP_BASE" => mktempdir())
    merge!(env, extra_env)
    out = IOBuffer()
    proc = run(pipeline(ignorestatus(Cmd(`bash $SBATCH`; env = env));
        stdout = out, stderr = out))
    return proc.exitcode, String(take!(out))
end

Test.@testset "td2rxh sbatch fail-closed guards" begin
    if Sys.which("bash") === nothing || Sys.which("git") === nothing
        Test.@test_skip "bash/git unavailable"
    else
        Test.@testset "clean run with the full row count succeeds" begin
            repo, home, bindir = fixture()
            rc, log = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "12"))
            Test.@test rc == 0
            Test.@test occursin("label=td2rxh@", log) &&
                       occursin("+manifestabcdef012345", log)
        end
        Test.@testset "seed and provenance label reach the harness" begin
            # The analysis rejects malformed or mixed labels, but only if the label
            # is actually passed through; a dropped --label would make every shard
            # "unlabeled" (review round 2, surviving mutants S1 and S4).
            repo, home, bindir = fixture()
            rc,
            _ = run_sbatch(repo, home, bindir;
                extra_env = Dict("STUB_ROWS" => "12", "SLURM_ARRAY_TASK_ID" => "3"))
            args = read(joinpath(home, "julia_args.txt"), String)
            Test.@test rc == 0
            Test.@test occursin(" --seeds 3 ", args)
            Test.@test occursin(r"--label td2rxh@[0-9a-f]{12}\+manifestabcdef012345", args)
            Test.@test occursin("--out results/td2rxh/99/seed3.csv", args)
            # The cell set is passed through, not hardcoded: a hardcoded td2rxh
            # would let a smoke submission run the full grid and still clear the
            # 4-row guard (round-7 review, surviving mutant).
            Test.@test occursin(" --cells td2rxh ", args)
            rc_smoke,
            _ = run_sbatch(repo, home, bindir;
                extra_env = Dict("STUB_ROWS" => "4", "TD2RXH_CELLS" => "smoke"))
            smoke_args = read(joinpath(home, "julia_args.txt"), String)
            Test.@test rc_smoke == 0
            Test.@test occursin(" --cells smoke ", smoke_args)
        end
        Test.@testset "the submit shell's project variables cannot redirect the run" begin
            # Each would run code the provenance label does not name: the harness
            # activates MYCELIA_PROJECT over --project, and a leading
            # JULIA_LOAD_PATH entry is searched before the active project.
            repo, home, bindir = fixture()
            rc,
            _ = run_sbatch(repo, home, bindir;
                extra_env = Dict(
                    "STUB_ROWS" => "12", "MYCELIA_PROJECT" => "/bogus/checkout",
                    "JULIA_LOAD_PATH" => "/bogus/env:", "JULIA_PROJECT" => "/bogus/project"))
            Test.@test rc == 0
            Test.@test readlines(joinpath(home, "julia_env.txt")) ==
                       ["MYCELIA_PROJECT=<unset>", "JULIA_LOAD_PATH=<unset>",
                "JULIA_PROJECT=<unset>"]
        end
        Test.@testset "julia exiting 0 with too few rows is reported as such" begin
            repo, home, bindir = fixture()
            rc, log = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "5"))
            Test.@test rc == 1 && occursin("julia exited 0 but", log)
        end
        Test.@testset "too few rows fails" begin
            repo, home, bindir = fixture()
            rc, log = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "7"))
            Test.@test rc == 1 && occursin("7 of 12", log)
        end
        Test.@testset "stale CSV from a previous attempt cannot pass" begin
            repo, home, bindir = fixture()
            stale = joinpath(repo, "benchmarking", "results", "td2rxh", "99", "seed1.csv")
            mkpath(dirname(stale))
            write(stale, join(["header"; fill("row", 12)], "\n") * "\n")
            rc, _ = run_sbatch(repo, home, bindir)   # stub writes nothing, exits 0
            Test.@test rc == 1
            Test.@test !isfile(stale)
        end
        Test.@testset "modified tracked file fails" begin
            repo, home, bindir = fixture()
            write(joinpath(repo, "benchmarking", "indel_reroute_evidence.jl"), "# edited\n")
            rc, log = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "12"))
            Test.@test rc == 1 && occursin("tracked files are modified", log)
        end
        Test.@testset "failing git fails closed" begin
            repo, home, bindir = fixture()
            write_exe(joinpath(bindir, "git"),
                "#!/bin/bash\necho 'fatal: dubious ownership' >&2\nexit 128\n")
            rc, log = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "12"))
            Test.@test rc == 1 && occursin("cannot resolve HEAD", log)
        end
        Test.@testset "git status failing on its own fails closed" begin
            # rev-parse succeeds, status fails: reaches the second git check, which
            # the rev-parse stub above cannot (review round 2, surviving mutant).
            repo, home, bindir = fixture()
            real_git = Sys.which("git")
            write_exe(
                joinpath(bindir, "git"), """
      #!/bin/bash
      for a in "\$@"; do
        if [[ "\$a" == "status" ]]; then echo 'fatal: dubious ownership' >&2; exit 128; fi
      done
      exec "$(real_git)" "\$@"
      """)
            rc, log = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "12"))
            Test.@test rc == 1 && occursin("git status failed", log)
        end
        Test.@testset "a failing sha256sum cannot yield an empty Manifest hash" begin
            repo, home, bindir = fixture()
            write_exe(joinpath(bindir, "sha256sum"), "#!/bin/bash\nexit 1\n")
            rc, log = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "12"))
            Test.@test rc == 1 && occursin("could not hash", log)
            # The 12-hex half on its own: sha256sum succeeds but prints garbage.
            write_exe(joinpath(bindir, "sha256sum"), "#!/bin/bash\necho \"garbage  \$1\"\n")
            rc2,
            log2 = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "12"))
            Test.@test rc2 == 1 && occursin("could not hash", log2)
        end
        Test.@testset "a per-task TMPDIR is resolved to its physical path and removed" begin
            # Mycelia's Stage-1 corrector rejects a temp root resolving through a
            # symlink; on Lawrencium TMPDIR sits under a symlinked $HOME/workspace.
            repo, home, bindir = fixture()
            physical = mktempdir()
            linked = joinpath(mktempdir(), "tmp_link")
            symlink(physical, linked)
            rc,
            log = run_sbatch(repo, home, bindir;
                extra_env = Dict("STUB_ROWS" => "12", "TD2RXH_TMP_BASE" => linked))
            m = match(r"TMPDIR=(\S+) \(ext2/ext3\)\n", log)
            Test.@test rc == 0
            Test.@test m !== nothing
            task_dir = m === nothing ? "" : m.captures[1]
            # Physical path (no symlink component), fresh mktemp name per task.
            Test.@test startswith(task_dir, joinpath(realpath(physical), "td2rxh_99_1."))
            Test.@test !isdir(task_dir)   # EXIT trap cleaned it up on success
        end
        Test.@testset "a failed run keeps its TMPDIR for debugging" begin
            repo, home, bindir = fixture()
            base = mktempdir()
            rc,
            log = run_sbatch(repo, home, bindir;
                extra_env = Dict("STUB_ROWS" => "12", "STUB_RC" => "3",
                    "TD2RXH_TMP_BASE" => base))
            kept = filter(startswith("td2rxh_99_1."), readdir(base))
            Test.@test rc == 1 && length(kept) == 1
            Test.@test occursin("kept TMPDIR", log)
        end
        Test.@testset "a non-local or unrecognised TMPDIR filesystem is refused" begin
            # renameat2(RENAME_NOREPLACE) returns EINVAL on Lustre (smoke job 27118360);
            # the check is an allowlist so an unknown type cannot slip through either.
            for fs in ("lustre", "nfs", "UNKNOWN (0x19830326)")
                repo, home, bindir = fixture()
                base = mktempdir()
                rc,
                log = run_sbatch(repo, home, bindir;
                    extra_env = Dict("STUB_ROWS" => "12", "TD2RXH_TMP_BASE" => base,
                        "STUB_FS" => fs))
                Test.@test rc == 1 && occursin("is on $(fs)", log)
                Test.@test isempty(readdir(base))   # pre-run refusal leaves nothing behind
            end
        end
        Test.@testset "julia non-zero exit fails even with full rows" begin
            repo, home, bindir = fixture()
            rc,
            log = run_sbatch(repo, home, bindir;
                extra_env = Dict("STUB_ROWS" => "12", "STUB_RC" => "3"))
            Test.@test rc == 1 && occursin("julia exited 3", log)
        end
        Test.@testset "missing Manifest.toml fails" begin
            repo, home, bindir = fixture()
            rm(joinpath(repo, "Manifest.toml"))
            rc, log = run_sbatch(repo, home, bindir; extra_env = Dict("STUB_ROWS" => "12"))
            Test.@test rc == 1 && occursin("Manifest.toml missing", log)
        end
        Test.@testset "smoke cell set expects 4 rows; unknown set fails" begin
            repo, home, bindir = fixture()
            rc,
            _ = run_sbatch(repo, home, bindir;
                extra_env = Dict("STUB_ROWS" => "4", "TD2RXH_CELLS" => "smoke"))
            Test.@test rc == 0
            rc2,
            log2 = run_sbatch(repo, home, bindir;
                extra_env = Dict("STUB_ROWS" => "12", "TD2RXH_CELLS" => "td2rhx"))
            Test.@test rc2 == 1 && occursin("unknown TD2RXH_CELLS=td2rhx", log2)
        end
    end
end

end
