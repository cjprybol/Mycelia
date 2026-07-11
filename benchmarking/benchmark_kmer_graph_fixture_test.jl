# Unit test for the reusable benchmark_kmer_graph_fixture helper.
#
# Proves the extracted path builds a WORKING correction graph + valid
# observations: it runs correct_observations end-to-end and shows an injected
# substitution is pulled back toward truth.
#
# Run:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/benchmark_kmer_graph_fixture_test.jl

import Test
import Mycelia

include(joinpath(@__DIR__, "benchmark_kmer_graph_fixture.jl"))

# per-position mismatch count over two equal-length k-mer-string vectors
_kmer_mismatches(a, b) = sum(x != y for (x, y) in zip(a, b))

Test.@testset "benchmark_kmer_graph_fixture" begin
    k = 11
    # Non-repetitive DNA truth (> k) so the truth-only graph has an unambiguous path.
    truth = "ACGTTGCAATCGGATCCTAGCATTGACCGATACGGTA"

    fx = benchmark_kmer_graph_fixture(truth, k; moltype = :DNA)

    # truth_kmers is the sliding-window decomposition.
    Test.@test fx.truth_kmers == [truth[i:(i + k - 1)] for i in 1:(length(truth) - k + 1)]
    Test.@test !isempty(fx.truth_kmers)

    # Clean observation has one k-mer per window and runs through the decoder.
    obs_clean = fx.to_observation(truth)
    Test.@test length(obs_clean) == length(fx.truth_kmers)

    # Inject one substitution and correct it via the fixture's graph + converter.
    errored = collect(truth)
    errored[18] = errored[18] == 'A' ? 'C' : 'A'
    obs_err = fx.to_observation(String(errored))
    cfg = Mycelia.ViterbiCorrectionConfig(error_rate = 0.05)
    result = Mycelia.correct_observations(fx.graph, [obs_err]; config = cfg)
    corrected = [string(o) for o in only(result.corrected_observations)]

    Test.@test length(corrected) == length(fx.truth_kmers)
    err_kstr = [string(o) for o in obs_err]
    # Correction moves the observation toward truth (never further from it).
    Test.@test _kmer_mismatches(corrected, fx.truth_kmers) <=
               _kmer_mismatches(err_kstr, fx.truth_kmers)

    # moltype guard.
    Test.@test_throws ArgumentError benchmark_kmer_graph_fixture(truth, k; moltype = :protein)
end
