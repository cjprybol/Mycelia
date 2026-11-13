# Self-test for the brute-force strand-mode oracle.
# From the Mycelia base directory, run with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/oracle_0_strand_selftest.jl")'
# ```
#
# The `oracle_0_` prefix makes `include_all_tests` in test/runtests.jl (which
# includes files in `sort(files)` order) run this self-test before the
# `oracle_*_test.jl` files that rely on the oracle.

import Test
import Random

if !isdefined(Main, :StrandOracle)
    include(joinpath(@__DIR__, "strand_oracle.jl"))
end

Test.@testset "StrandOracle self-test (hand-derived expectations)" begin
    Test.@testset "fwd_counts" begin
        v, e = StrandOracle.fwd_counts(["ACGTA"], 3)
        Test.@test v == Dict("ACG" => 1, "CGT" => 1, "GTA" => 1)
        Test.@test e == Dict(("ACG", "CGT") => 1, ("CGT", "GTA") => 1)
    end

    Test.@testset "revcomp and canon" begin
        Test.@test StrandOracle.revcomp("ACCGT") == "ACGGT"
        Test.@test StrandOracle.canon("TAC") == "GTA"
        Test.@test StrandOracle.canon("ACG") == "ACG"
    end

    Test.@testset "doublestrand_expected" begin
        v, e = StrandOracle.doublestrand_expected(["AAC"], 2)
        Test.@test v == Dict("AA" => 1, "AC" => 1, "TT" => 1, "GT" => 1)
        Test.@test e == Dict(("AA", "AC") => 1, ("GT", "TT") => 1)

        v, _ = StrandOracle.doublestrand_expected(["ACGT"], 4)
        Test.@test v == Dict("ACGT" => 2)

        # Self-reverse-complement edge ACG -> CGT (rc(ACG) == CGT): its
        # reverse-complement image is itself, so it counts twice.
        Test.@test StrandOracle.doublestrand_expected(["ACGT"], 3) ==
                   (Dict("ACG" => 2, "CGT" => 2), Dict(("ACG", "CGT") => 2))
    end

    Test.@testset "canonical_expected" begin
        v, e = StrandOracle.canonical_expected(["AAC"], 2)
        Test.@test v == Dict("AA" => 1, "AC" => 1)
        Test.@test e == Dict(("AA", "AC") => 1)

        # k-mers ACG and CGT both canonicalize to ACG; the edge ACG -> CGT is
        # self-reverse-complement, so the self-loop counts on both strands.
        Test.@test StrandOracle.canonical_expected(["ACGT"], 3) ==
                   (Dict("ACG" => 2), Dict(("ACG", "ACG") => 2))

        # Palindromic k-mer ACGT at k = 4 counts on both strands.
        v, _ = StrandOracle.canonical_expected(["ACGT"], 4)
        Test.@test v == Dict("ACGT" => 2)
    end

    Test.@testset "random_non_rc_closed" begin
        s = StrandOracle.random_non_rc_closed(Random.MersenneTwister(1), 30, 5)
        Test.@test length(s) == 30
        kmers, _ = StrandOracle.fwd_counts([s], 5)
        S = Set(keys(kmers))
        Test.@test any(x -> !(StrandOracle.revcomp(x) in S), S)
    end
end
