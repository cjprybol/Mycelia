# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_emission_alloc_test.jl")'
# ```
#
# td-ve02: the corrector's dominant cost was per-emission-call allocation —
# `default_viterbi_emission_logp` materialized `uppercase(string(unit))` for BOTH
# operands on every call (~4 string allocs × billions of calls = 56B allocations
# on Lambda). k-mer/BioSequence operands now take an allocation-free fast path
# that indexes symbols directly. This is a tracked regression: the k-mer path must
# allocate strictly less than the string fallback, and stay correct.

import BioSequences
import Kmers
import Mycelia
import Test

Test.@testset "allocation-free emission fast path (td-ve02)" begin
    k_match_a = Kmers.DNAKmer{5}("ATGCG")
    k_match_b = Kmers.DNAKmer{5}("ATGCG")
    k_mismatch = Kmers.DNAKmer{5}("ATGCA")

    # Warm up compilation for both paths before measuring.
    Mycelia.default_viterbi_emission_logp(k_match_a, k_match_b, :DNA)
    Mycelia.default_viterbi_emission_logp("ATGCG", "ATGCA", :DNA)

    Test.@testset "k-mer fast path allocates less than the string fallback" begin
        kmer_allocs = @allocated Mycelia.default_viterbi_emission_logp(k_match_a, k_mismatch, :DNA)
        str_allocs = @allocated Mycelia.default_viterbi_emission_logp("ATGCG", "ATGCA", :DNA)
        Test.@test kmer_allocs < str_allocs
    end

    Test.@testset "fast path is byte-equivalent to the string path" begin
        # Perfect match scores higher than a single mismatch, and the k-mer and
        # string computations agree exactly.
        m = Mycelia.default_viterbi_emission_logp(k_match_a, k_match_b, :DNA)
        mm = Mycelia.default_viterbi_emission_logp(k_match_a, k_mismatch, :DNA)
        Test.@test m > mm
        Test.@test m ≈ Mycelia.default_viterbi_emission_logp("ATGCG", "ATGCG", :DNA)
        Test.@test mm ≈ Mycelia.default_viterbi_emission_logp("ATGCG", "ATGCA", :DNA)
    end
end
