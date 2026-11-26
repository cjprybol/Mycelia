# Guard for the shared Stage-0 harness (td-f581, PR #347 review I1).
#
# benchmarking/stage0_classification_matrix.jl and
# benchmarking/stage0_threshold_calibration.jl must both include
# benchmarking/stage0_common.jl, resolve the same qualmer memory_profile
# (:full), and must not redefine sample_quality, simulate, or
# reference_kmer_set. Comments and string literals do not count.
#
# The file is picked up by test/runtests.jl (include_all_tests on
# test/4_assembly). Run it directly with:
#   julia --project=. test/4_assembly/stage0_common_harness_test.jl

import Test
import Random
import FASTX
import BioSequences
import Kmers

const _STAGE0_HARNESS_HELPERS = ("sample_quality", "simulate", "reference_kmer_set")
const _STAGE0_HARNESS_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _STAGE0_HARNESS_COMMON = joinpath(
    _STAGE0_HARNESS_ROOT, "benchmarking", "stage0_common.jl")
const _STAGE0_HARNESS_SCRIPTS = (
    joinpath(_STAGE0_HARNESS_ROOT, "benchmarking", "stage0_classification_matrix.jl"),
    joinpath(_STAGE0_HARNESS_ROOT, "benchmarking", "stage0_threshold_calibration.jl"),
)

"""Code view of Julia source.

Comments are removed. String and char literals are removed when
`blank_literals` is true, so a docstring cannot look like a definition.
Newlines are kept so `(?m)^` anchors still see statement boundaries.
"""
function _stage0_harness_code(src::String; blank_literals::Bool)::String
    chars = collect(src)
    n = length(chars)
    buf = Char[]
    sizehint!(buf, n)
    i = 1
    while i <= n
        c = chars[i]
        if c == '#'
            if i < n && chars[i + 1] == '='
                depth = 1
                i += 2
                while i <= n && depth > 0
                    if chars[i] == '#' && i < n && chars[i + 1] == '='
                        depth += 1
                        i += 2
                    elseif chars[i] == '=' && i < n && chars[i + 1] == '#'
                        depth -= 1
                        i += 2
                    else
                        chars[i] == '\n' && push!(buf, '\n')
                        i += 1
                    end
                end
            else
                while i <= n && chars[i] != '\n'
                    i += 1
                end
            end
        elseif c == '"'
            triple = i + 2 <= n && chars[i + 1] == '"' && chars[i + 2] == '"'
            if triple
                if !blank_literals
                    push!(buf, '"', '"', '"')
                end
                i += 3
                while i <= n
                    if chars[i] == '"' && i + 2 <= n &&
                       chars[i + 1] == '"' && chars[i + 2] == '"'
                        if !blank_literals
                            push!(buf, '"', '"', '"')
                        end
                        i += 3
                        break
                    end
                    if chars[i] == '\n'
                        push!(buf, '\n')
                    elseif !blank_literals
                        push!(buf, chars[i])
                    end
                    i += 1
                end
            else
                blank_literals || push!(buf, '"')
                i += 1
                while i <= n
                    if chars[i] == '\\' && i < n
                        if !blank_literals
                            push!(buf, chars[i], chars[i + 1])
                        end
                        i += 2
                        continue
                    end
                    if chars[i] == '"'
                        blank_literals || push!(buf, '"')
                        i += 1
                        break
                    end
                    if chars[i] == '\n'
                        push!(buf, '\n')
                    elseif !blank_literals
                        push!(buf, chars[i])
                    end
                    i += 1
                end
            end
        elseif c == '\''
            # Char literal: 'A' or '\n'. Dropped with the string literals so
            # it cannot contribute a fake definition. The closing quote is
            # consumed; surrounding code stays.
            i += 1
            if i <= n && chars[i] == '\\'
                i += 1
            end
            i <= n && (i += 1)
            if i <= n && chars[i] == '\''
                i += 1
            end
        else
            push!(buf, c)
            i += 1
        end
    end
    return String(buf)
end

function _stage0_harness_defined_functions(src::String)::Set{String}
    code = _stage0_harness_code(src; blank_literals = true)
    names = Set{String}()
    for m in eachmatch(r"(?m)^[ \t]*function[ \t]+([A-Za-z_][A-Za-z0-9_]*)\b", code)
        push!(names, m.captures[1])
    end
    short = r"(?m)^[ \t]*([A-Za-z_][A-Za-z0-9_]*)[ \t]*\([^()\n]*\)[ \t]*(?:::[^\n=]+)?="
    for m in eachmatch(short, code)
        push!(names, m.captures[1])
    end
    return names
end

function _stage0_harness_consts(src::String)::Dict{String, String}
    code = _stage0_harness_code(src; blank_literals = true)
    bindings = Dict{String, String}()
    pat = r"(?m)^[ \t]*const[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*=[ \t]*(:[A-Za-z_][A-Za-z0-9_]*)\b"
    for m in eachmatch(pat, code)
        bindings[m.captures[1]] = m.captures[2]
    end
    return bindings
end

function _stage0_harness_resolve(token::AbstractString, consts::Dict{String, String})::String
    startswith(token, ":") && return String(token)
    seen = Set{String}()
    current = String(token)
    while haskey(consts, current)
        current in seen && return "cycle:$current"
        push!(seen, current)
        current = consts[current]
        startswith(current, ":") && return current
    end
    return "unresolved:$current"
end

"""Resolved `memory_profile = ...` values in source order.

A literal `:` symbol is itself. An identifier is resolved through consts in
`src`, then consts in `common_src`.
"""
function _stage0_harness_profiles(src::String, common_src::String)::Vector{String}
    code = _stage0_harness_code(src; blank_literals = true)
    consts = _stage0_harness_consts(common_src)
    for (name, value) in _stage0_harness_consts(src)
        consts[name] = value
    end
    profiles = String[]
    pat = r"(?<![A-Za-z0-9_])memory_profile[ \t]*=[ \t]*((?::[A-Za-z_][A-Za-z0-9_]*)|(?:[A-Za-z_][A-Za-z0-9_]*))"
    for m in eachmatch(pat, code)
        push!(profiles, _stage0_harness_resolve(m.captures[1], consts))
    end
    return profiles
end

function _stage0_harness_walk(f, ex)
    f(ex)
    if ex isa Expr
        for arg in ex.args
            _stage0_harness_walk(f, arg)
        end
    end
    return nothing
end

"""True when `ex` is an `include(...)` call whose arguments name stage0_common.jl.

A docstring or comment that only mentions the filename is not a call, so it
does not count. `Meta.parseall` drops comments before this walk.
"""
function _stage0_harness_include_call(ex)::Bool
    ex isa Expr || return false
    ex.head === :call || return false
    length(ex.args) >= 2 || return false
    ex.args[1] === :include || return false
    return _stage0_harness_expr_mentions_common(ex)
end

function _stage0_harness_expr_mentions_common(ex)::Bool
    if ex isa AbstractString
        return occursin("stage0_common.jl", ex)
    elseif ex isa Expr
        return any(_stage0_harness_expr_mentions_common, ex.args)
    else
        return false
    end
end

function _stage0_harness_includes_common(src::String)::Bool
    ast = try
        Meta.parseall(src)
    catch
        return false
    end
    found = Ref(false)
    _stage0_harness_walk(ast) do ex
        if _stage0_harness_include_call(ex)
            found[] = true
        end
    end
    return found[]
end

"""Default k from `const K = getenv("MYCELIA_S0_K", k)` in a Stage-0 script."""
function _stage0_harness_configured_k(src::String)::Int
    found = Int[]
    ast = Meta.parseall(src)
    _stage0_harness_walk(ast) do ex
        ex isa Expr || return
        ex.head === :const || return
        length(ex.args) == 1 || return
        inner = ex.args[1]
        inner isa Expr && inner.head === :(=) || return
        length(inner.args) == 2 || return
        inner.args[1] === :K || return
        rhs = inner.args[2]
        rhs isa Expr && rhs.head === :call || return
        length(rhs.args) >= 3 || return
        rhs.args[1] === :getenv || return
        rhs.args[2] == "MYCELIA_S0_K" || return
        rhs.args[3] isa Integer || return
        push!(found, Int(rhs.args[3]))
    end
    length(found) == 1 ||
        error("expected one MYCELIA_S0_K default, found $(found)")
    return found[1]
end

"""Problems that mean this script has drifted from the shared harness."""
function _stage0_harness_issues(label::AbstractString, src::String,
                                common_src::String)::Vector{String}
    issues = String[]
    if !_stage0_harness_includes_common(src)
        push!(issues, "$label does not include stage0_common.jl")
    end
    defined = _stage0_harness_defined_functions(common_src)
    for helper in _STAGE0_HARNESS_HELPERS
        if helper ∉ defined
            push!(issues, "stage0_common.jl does not define $helper")
        end
    end
    redefs = [helper for helper in _STAGE0_HARNESS_HELPERS
              if helper in _stage0_harness_defined_functions(src)]
    if !isempty(redefs)
        push!(issues, "$label redefines shared helper(s): $(join(redefs, ", "))")
    end
    profiles = _stage0_harness_profiles(src, common_src)
    if isempty(profiles)
        push!(issues, "$label does not set memory_profile")
    elseif any(profile -> profile != ":full", profiles)
        push!(issues, "$label memory_profile resolved to $profiles; expected [:full]")
    end
    return issues
end

# One smoke module per k configured by the Stage-0 scripts. The module's `K`
# is that script default, and `read_kmers` closes over it, so the k-mer
# iterator is not a literal in the test.
const _STAGE0_HARNESS_SMOKES = Dict{Int, Module}()
for _stage0_k in unique(
    _stage0_harness_configured_k(read(path, String))
    for path in _STAGE0_HARNESS_SCRIPTS
)
    _stage0_mod = Symbol("Stage0CommonHarnessSmoke_k", _stage0_k)
    @eval module $_stage0_mod
        # Error rate 0 makes every read an exact slice of the reference.
        # Read length stays just above this script k so the draw stays small
        # even when a script's default reference is long.
        const K = $_stage0_k
        const READLEN = K + 9
        const REFLEN = 2 * READLEN
        const COVERAGE = 2
        const ERROR_RATE = 0.0
        include($(_STAGE0_HARNESS_COMMON))
        function read_kmers(seq)
            return Kmers.UnambiguousDNAMers{K}(seq)
        end
    end
    _STAGE0_HARNESS_SMOKES[_stage0_k] = getfield(@__MODULE__, _stage0_mod)
end

Test.@testset "Stage-0 shared harness (td-f581)" begin
    common_src = read(_STAGE0_HARNESS_COMMON, String)
    script_srcs = [label => read(path, String)
                   for (label, path) in (
                       ("stage0_classification_matrix.jl", _STAGE0_HARNESS_SCRIPTS[1]),
                       ("stage0_threshold_calibration.jl", _STAGE0_HARNESS_SCRIPTS[2]),
                   )]

    Test.@testset "both scripts use the one harness" begin
        for (label, src) in script_srcs
            issues = _stage0_harness_issues(label, src, common_src)
            Test.@test issues == String[]
        end
        # One resolved value per script. A temporary edit of either
        # memory_profile literal makes these two vectors disagree.
        resolved = [_stage0_harness_profiles(src, common_src) for (_, src) in script_srcs]
        Test.@test resolved == [[":full"], [":full"]]
        Test.@test _stage0_harness_consts(common_src)["STAGE0_MEMORY_PROFILE"] == ":full"
    end

    Test.@testset "detector ignores comments and docstrings" begin
        common = """
        const STAGE0_MEMORY_PROFILE = :full
        function sample_quality(rng, is_error::Bool)::Char
            return 'I'
        end
        function simulate(rng)
            return "AC"
        end
        function reference_kmer_set(ref::AbstractString)
            return Set([ref])
        end
        """
        quiet = """
        include(joinpath(@__DIR__, "stage0_common.jl"))
        # memory_profile = :lightweight_quality
        # function sample_quality(rng, is_error::Bool)::Char
        \"\"\"
        function sample_quality(rng, is_error::Bool)::Char
        memory_profile = :lightweight_quality
        \"\"\"
        memory_profile = STAGE0_MEMORY_PROFILE
        """
        Test.@test _stage0_harness_issues("quiet.jl", quiet, common) == String[]
        Test.@test _stage0_harness_profiles(quiet, common) == [":full"]
    end

    Test.@testset "detector fires when memory_profile values diverge" begin
        common = "const STAGE0_MEMORY_PROFILE = :full\n"
        left = """
        include(joinpath(@__DIR__, "stage0_common.jl"))
        memory_profile = STAGE0_MEMORY_PROFILE
        """
        right = """
        include(joinpath(@__DIR__, "stage0_common.jl"))
        memory_profile = :lightweight_quality
        """
        Test.@test _stage0_harness_profiles(left, common) == [":full"]
        Test.@test _stage0_harness_profiles(right, common) == [":lightweight_quality"]
        Test.@test _stage0_harness_profiles(left, common) != _stage0_harness_profiles(right, common)
        drifted = _stage0_harness_issues("right.jl", right, common)
        Test.@test any(issue -> occursin("memory_profile", issue) &&
                                 occursin(":lightweight_quality", issue), drifted)
    end

    Test.@testset "detector fires when a script redefines a shared helper" begin
        common = """
        const STAGE0_MEMORY_PROFILE = :full
        function sample_quality(rng, is_error::Bool)::Char
            return 'I'
        end
        function simulate(rng)
            return "AC"
        end
        function reference_kmer_set(ref::AbstractString)
            return Set([ref])
        end
        """
        script = """
        include(joinpath(@__DIR__, "stage0_common.jl"))
        memory_profile = STAGE0_MEMORY_PROFILE
        function sample_quality(rng, is_error::Bool)::Char
            return '!'
        end
        """
        issues = _stage0_harness_issues("redef.jl", script, common)
        Test.@test any(issue -> occursin("redefines shared helper", issue) &&
                                 occursin("sample_quality", issue), issues)
    end

    Test.@testset "detector rejects a mention that is not an include call" begin
        common = """
        const STAGE0_MEMORY_PROFILE = :full
        function sample_quality(rng, is_error::Bool)::Char
            return 'I'
        end
        function simulate(rng)
            return "AC"
        end
        function reference_kmer_set(ref::AbstractString)
            return Set([ref])
        end
        """
        # The filename appears only inside a docstring. There is no include call.
        doc_only = """
        \"\"\"
        Shared model: include(joinpath(@__DIR__, "stage0_common.jl"))
        \"\"\"
        memory_profile = :full
        """
        comment_only = """
        # include(joinpath(@__DIR__, "stage0_common.jl"))
        memory_profile = :full
        """
        Test.@test !_stage0_harness_includes_common(doc_only)
        Test.@test !_stage0_harness_includes_common(comment_only)
        for (label, src) in (("doc_only.jl", doc_only), ("comment_only.jl", comment_only))
            issues = _stage0_harness_issues(label, src, common)
            Test.@test issues == ["$label does not include stage0_common.jl"]
        end
    end

    Test.@testset "shared helpers still simulate the Stage-0 model" begin
        checked = Set{Int}()
        for (_, src) in script_srcs
            k = _stage0_harness_configured_k(src)
            smoke = _STAGE0_HARNESS_SMOKES[k]
            # The draw uses the k configured by this script, not a literal in the test.
            Test.@test smoke.K == k
            k in checked && continue
            push!(checked, k)
            Test.@test smoke.STAGE0_MEMORY_PROFILE === :full
            rng = Random.MersenneTwister(7)
            for is_error in (false, true)
                quality = smoke.sample_quality(rng, is_error)
                Test.@test quality isa Char
                Test.@test '#' <= quality <= 'I'   # Phred 2..40, encoded +33
            end
            ref_a, recs_a = smoke.simulate(Random.MersenneTwister(7))
            ref_b, recs_b = smoke.simulate(Random.MersenneTwister(7))
            Test.@test ref_a == ref_b
            Test.@test length(ref_a) == smoke.REFLEN
            n_reads = max(1, round(Int, smoke.COVERAGE * smoke.REFLEN / smoke.READLEN))
            Test.@test length(recs_a) == n_reads
            Test.@test length(recs_a) == length(recs_b)
            truth = smoke.reference_kmer_set(ref_a)
            Test.@test !isempty(truth)
            Test.@test smoke.READLEN >= smoke.K
            for (rec_a, rec_b) in zip(recs_a, recs_b)
                seq_a = FASTX.sequence(String, rec_a)
                Test.@test seq_a == FASTX.sequence(String, rec_b)
                Test.@test length(seq_a) == smoke.READLEN
                Test.@test occursin(seq_a, ref_a)
                read_seq = BioSequences.LongDNA{4}(seq_a)
                for (kmer, _) in smoke.read_kmers(read_seq)
                    Test.@test BioSequences.canonical(kmer) in truth
                end
            end
            for kmer in truth
                Test.@test kmer == BioSequences.canonical(kmer)
            end
        end
        Test.@test !isempty(checked)
    end
end
