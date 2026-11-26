# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/ont_identity_selftest_nits_test.jl")'
# ```
#
# td-z2e0: PR #452 deferred selftest/docs nits. This test reads the
# selftest, the sbatch wrapper, and the k-sweep README (it does not
# submit a cluster job).

import Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SELFTEST = joinpath(REPO_ROOT, "benchmarking",
    "run_ont_read_identity_lrc.selftest.sh")
const SBATCH = joinpath(REPO_ROOT, "benchmarking",
    "run_ont_read_identity_lrc.sbatch")
const README = joinpath(REPO_ROOT, "benchmarking", "results",
    "ont_k_sweep", "README.md")

selftest = read(SELFTEST, String)
sbatch = read(SBATCH, String)
readme = read(README, String)

# The wrapper's leading comment block, including the green-exit contract.
header = begin
    idx = findfirst("set -uo pipefail", sbatch)
    idx === nothing ? sbatch : sbatch[1:first(idx)-1]
end

# Provenance intro: from the heading up to the first markdown table.
prov = begin
    m = match(r"## Provenance\n(.*?)(?:\n\|)"s, readme)
    m === nothing ? "" : m.captures[1]
end

Test.@testset "ONT identity selftest and docs nits (td-z2e0)" begin
    Test.@testset "zero-byte empty ladder case" begin
        Test.@test occursin(r"empty_ladder\)", selftest)
        Test.@test occursin(r"STUB_MODE=empty_ladder", selftest)
        Test.@test occursin(r":\s*>\s*\"\$out/kmer_survival_ladder\.tsv\"",
            selftest)
    end

    Test.@testset "reused OUT_DIR uses its own job id" begin
        Test.@test occursin(r"reused \(non-empty\) OUT_DIR", selftest)
        good = match(r"run_case\s+\"good run\"\s+pass\s+(\d+)", selftest)
        reused = match(r"run_case\s+\"reused \(non-empty\) OUT_DIR\"\s+fail\s+(\d+)",
            selftest)
        Test.@test good !== nothing && reused !== nothing
        if good !== nothing && reused !== nothing
            Test.@test good.captures[1] != reused.captures[1]
            jid = reused.captures[1]
            Test.@test occursin(Regex("out/" * jid * "/leftover\\.txt"),
                selftest)
        end
    end

    Test.@testset "sbatch green-exit header names a data row" begin
        Test.@test occursin(r"green exit", header)
        Test.@test occursin(r"at least one data row", header)
    end

    Test.@testset "k-sweep README Provenance pointer" begin
        Test.@test occursin(r"## Provenance", readme)
        Test.@test occursin("Regenerated on a different host", readme)
        Test.@test occursin("Regenerated on a different host", prov)
        Test.@test !occursin("paragraph after the table", prov)
    end
end
