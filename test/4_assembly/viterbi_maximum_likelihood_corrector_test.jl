# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_maximum_likelihood_corrector_test.jl")'
# ```

import Test
import Mycelia
import FASTX
import Kmers
import Random

Test.@testset "Legacy Viterbi maximum-likelihood corrector" begin
    fixture_dir = joinpath(@__DIR__, "..", "fixtures", "viterbi_corrector")
    input_fasta = joinpath(fixture_dir, "dna_error_injection.fasta")
    oracle_tsv = joinpath(fixture_dir, "dna_error_injection.oracle.tsv")

    records = collect(Mycelia.open_fastx(input_fasta))
    graph = Mycelia.build_stranded_kmer_graph(Kmers.DNAKmer{4}, records)

    Random.seed!(1)
    correction = Mycelia.correct_observations(
        graph;
        config = Mycelia.ViterbiCorrectionConfig(
            error_rate = 0.05,
            verbosity = "dataset"
        )
    )
    corrected_records = correction.corrected_observations

    Test.@test correction.diagnostics[:interface] == :legacy_stranded_kmer
    Test.@test correction.diagnostics[:algorithm] == :viterbi_maximum_likelihood_traversals

    observed = Dict(
        String(FASTX.identifier(record)) => string(FASTX.sequence(record))
        for record in corrected_records
    )
    expected = Dict{String, String}()
    for line in eachline(oracle_tsv)
        isempty(strip(line)) && continue
        identifier, sequence = split(line, '\t')
        expected[identifier] = sequence
    end

    Test.@test observed == expected
end
