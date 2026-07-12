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

# Per-position mismatch count over two equal-length k-mer-string vectors.
function _kmer_mismatches(a::AbstractVector, b::AbstractVector)::Int
    length(a) == length(b) || throw(DimensionMismatch("k-mer vectors must have equal length"))
    return sum(x != y for (x, y) in zip(a, b))
end

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
    clean_result = Mycelia.correct_observations(fx.graph, [obs_clean])
    clean_corrected = [string(o) for o in only(clean_result.corrected_observations)]
    Test.@test clean_corrected == fx.truth_kmers

    # Inject one substitution and correct it via the fixture's graph + converter.
    errored = collect(truth)
    errored[18] = errored[18] == 'A' ? 'C' : 'A'
    obs_err = fx.to_observation(String(errored))
    cfg = Mycelia.ViterbiCorrectionConfig(error_rate = 0.05)
    result = Mycelia.correct_observations(fx.graph, [obs_err]; config = cfg)
    corrected = [string(o) for o in only(result.corrected_observations)]

    Test.@test length(corrected) == length(fx.truth_kmers)
    err_kstr = [string(o) for o in obs_err]
    # Correction strictly moves the observation toward truth.
    Test.@test _kmer_mismatches(err_kstr, fx.truth_kmers) > 0
    Test.@test _kmer_mismatches(corrected, fx.truth_kmers) <
               _kmer_mismatches(err_kstr, fx.truth_kmers)

    # Distinct RNA and text graph/observation branches run end-to-end.
    rna_truth = replace(truth, 'T' => 'U')
    for (sequence, moltype) in ((rna_truth, :RNA), ("thequickbrownfox", :text))
        branch_fx = benchmark_kmer_graph_fixture(sequence, 5; moltype = moltype)
        branch_obs = branch_fx.to_observation(sequence)
        branch_result = Mycelia.correct_observations(branch_fx.graph, [branch_obs])
        branch_corrected = [string(o) for o in only(branch_result.corrected_observations)]
        Test.@test branch_corrected == branch_fx.truth_kmers
    end

    # Input guards and AbstractString dataset identifiers.
    dataset_id = SubString("fixture-id", 1, 7)
    Test.@test !isempty(benchmark_kmer_graph_fixture(truth, k; dataset_id = dataset_id).truth_kmers)
    Test.@test_throws ArgumentError benchmark_kmer_graph_fixture(truth, k; moltype = :protein)
    Test.@test_throws ArgumentError benchmark_kmer_graph_fixture(truth, 0)
    Test.@test_throws ArgumentError benchmark_kmer_graph_fixture(truth, length(truth) + 1)
    Test.@test_throws ArgumentError benchmark_kmer_graph_fixture("ACGTβ", 3)
end
