# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/workflow_map_section4_nits_test.jl")'
# ```
#
# td-1i5q: PR #461 deferred docs nits in workflow-map section 4 and the
# metaMDBG docstring.

import Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const WORKFLOW_MAP = joinpath(REPO_ROOT, "docs", "src", "workflow-map.md")
const ASSEMBLY_JL = joinpath(REPO_ROOT, "src", "assembly.jl")

function section4(text::AbstractString)
    m = match(r"## 4\. Hybrid assembly wrapper contracts\n"s, text)
    Test.@test m !== nothing
    rest = text[m.offset:end]
    nxt = findnext(r"\n## ", rest, 2)
    return nxt === nothing ? rest : rest[1:nxt]
end

Test.@testset "workflow-map section 4 and metaMDBG docstring (td-1i5q)" begin
    Test.@test isfile(WORKFLOW_MAP)
    Test.@test isfile(ASSEMBLY_JL)
    doc = read(WORKFLOW_MAP, String)
    s4 = section4(doc)
    src = read(ASSEMBLY_JL, String)

    Test.@testset "4.1 states two option-rejection mechanisms" begin
        Test.@test occursin("_reject_route_managed_options", s4) ||
                   occursin(r"denylist", s4)
        Test.@test occursin("_reject_unsupported_hybrid_options", s4) ||
                   occursin(r"allowlist", s4)
        # Must not present both as a single "rejected outright if passed through
        # correction_options or assembler_options" rule with no second mechanism.
        Test.@test occursin(r"different (subsets|mechanisms|exception)", s4) ||
                   occursin(r"two (mechanisms|option bags|rejection)", s4)
    end

    Test.@testset "4.2/4.1 states the 1 TB snapshot ceiling default" begin
        Test.@test occursin("input_snapshot_byte_ceiling", s4)
        Test.@test occursin(r"1\s*TB", s4) || occursin("1_000_000_000_000", s4)
    end

    Test.@testset "Autocycler Polypolish default --careful is stated" begin
        Test.@test occursin(r"(?i)Polypolish", s4)
        # Must name Polypolish's own default --careful, not only careful Pypolca.
        Test.@test occursin(r"(?i)Polypolish.{0,40}--careful", s4)
    end

    Test.@testset "section 4.3 lines are at most 80 characters" begin
        lines = split(s4, '\n')
        # Only the 4.3 subsection.
        idx = findfirst(l -> startswith(l, "### 4.3"), lines)
        Test.@test idx !== nothing
        for line in lines[idx:end]
            startswith(line, "### ") && line != lines[idx] && break
            startswith(line, "## ") && break
            Test.@test length(line) <= 80
        end
    end

    Test.@testset "metaMDBG docstring does not claim unverified upstream rejection" begin
        Test.@test !occursin("metaMDBG v1.4 rejects simultaneous", src)
        Test.@test occursin(r"emits exactly one", src) ||
                   occursin(r"exactly one of `--in-hifi` and `--in-ont`", src)
    end
end
