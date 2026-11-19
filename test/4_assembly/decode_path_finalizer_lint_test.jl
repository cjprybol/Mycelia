# DECODE-PATH GC-NEUTRALITY LINT (td-3vge item 5)
# ==============================================
#
# opt5's between-batch GC.gc() is output-neutral only while decode-path modules
# do not observe collection. This test locks the CI lint
# (scripts/check_decode_path_finalizers.py, wired from .github/workflows/ci.yml):
# the current decode-path set is clean, a missing required file fails closed,
# and a new finalizer / WeakRef / WeakKeyDict in code fails the guard.
# Comments and string literals do not count as introductions. A '#' or a
# newline inside a "..." string must not hide a later call. String
# interpolation is out of scope, same as eval of a constructed token.
#
# No Mycelia load. Run directly:
#   julia --project=. test/4_assembly/decode_path_finalizer_lint_test.jl

import Test

const _REPO = normpath(joinpath(@__DIR__, "..", ".."))
const _SCRIPT = joinpath(_REPO, "scripts", "check_decode_path_finalizers.py")

function _lint(args::Vector{String})
    stdout_buf = IOBuffer()
    stderr_buf = IOBuffer()
    cmd = ignorestatus(`python3 $(_SCRIPT) $args`)
    proc = run(pipeline(cmd; stdout = stdout_buf, stderr = stderr_buf))
    return (
        exitcode = proc.exitcode,
        out = String(take!(stdout_buf)),
        err = String(take!(stderr_buf)),
    )
end

Test.@testset "decode-path finalizer lint (td-3vge)" begin
    if Sys.which("python3") === nothing
        Test.@test_skip "python3 is not on PATH; CI lint runs scripts/check_decode_path_finalizers.py"
        return
    end

    Test.@test isfile(_SCRIPT)

    Test.@testset "current decode path is clean and fully listed" begin
        result = _lint(String[])
        Test.@test result.exitcode == 0
        Test.@test isempty(result.err)
        for rel in (
                "src/iterative-assembly.jl",
                "src/viterbi-next.jl",
                "src/viterbi-polishing-and-error-correction.jl",
                "src/rhizomorph/algorithms/path-finding.jl",
                "src/rhizomorph/algorithms/batched-viterbi-poc.jl",
                "src/rhizomorph/algorithms/batched-viterbi-kernel.jl",
                "ext/MyceliaKernelAbstractionsExt.jl",
            )
            Test.@test occursin(rel, result.out)
        end
    end

    Test.@testset "a missing required decode-path file fails closed" begin
        empty = mktempdir()
        result = _lint(["--root", empty])
        Test.@test result.exitcode != 0
        Test.@test occursin("missing decode-path file:", result.err)
        Test.@test occursin("src/iterative-assembly.jl", result.err)
        Test.@test occursin("src/viterbi-next.jl", result.err)
    end

    Test.@testset "code introductions fail; comments and strings do not" begin
        tmp = mktempdir()
        cases = (
            ("call.jl", "finalizer(identity, Ref(1))\n", false, "finalizer"),
            ("weakref.jl", "cache = WeakRef(something)\n", false, "WeakRef"),
            ("weakkey.jl", "cache = WeakKeyDict{Int,Int}()\n", false, "WeakKeyDict"),
            ("qualified.jl", "Base.finalizer(identity, Ref(1))\n", false, "finalizer"),
            ("comment.jl", "# do not add finalizer, WeakRef, or WeakKeyDict\n", true, ""),
            ("block.jl", "#= WeakRef and WeakKeyDict stay out =#\nx = 1\n", true, ""),
            ("nested.jl", "#= outer #= finalizer =# still =#\nx = 1\n", true, ""),
            ("string.jl", "msg = \"finalizer WeakRef WeakKeyDict\"\n", true, ""),
            ("hash_string.jl", "s = \"a#b\"; finalizer(identity, Ref(1))\n", false, "finalizer"),
            ("multiline_string.jl", "s = \"a\nb\"; finalizer(identity, Ref(1))\n", false, "finalizer"),
            ("other_name.jl", "my_finalizer = 1\n", true, ""),
        )
        for (name, body, ok, token) in cases
            path = joinpath(tmp, name)
            write(path, body)
            result = _lint([path])
            if ok
                Test.@test result.exitcode == 0
                Test.@test isempty(result.err)
            else
                Test.@test result.exitcode != 0
                Test.@test occursin(token, result.err)
                Test.@test occursin(name, result.err)
            end
        end
    end
end
