# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_local_decode_test.jl")'
# ```
#
# Tier 1A (td-ve02): seed-anchored local decode. The corrector optionally decodes
# each read against the sub-dBG reachable within `local_radius` hops of the read's
# exact-match anchor k-mers, instead of the whole graph — bounding the per-read
# state set and eliminating the all-labels start fallback. This is provably
# lossless when the radius contains the read's true path: the induced subgraph
# then contains the optimal path, so the corrected read is IDENTICAL to the
# full-graph decode. The guardrail below asserts that invariance.

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

decoded_strings(result) = [string(l) for l in only(result.corrected_observations)]

Test.@testset "Viterbi seed-anchored local decode (Tier 1A, td-ve02)" begin
    # Linear fixture — the whole graph is the read's neighborhood.
    linear = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
    lgraph = Mycelia.Rhizomorph.build_kmer_graph(
        linear, 3; dataset_id="local_linear", mode=:singlestrand)
    observed = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGA"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
    ]

    Test.@testset "default local_radius=0 is the exact full-graph decode" begin
        r = Mycelia.correct_observations(lgraph, [observed])
        Test.@test decoded_strings(r) == ["ATG", "TGC", "GCG", "CGT"]
    end

    Test.@testset "local decode == full decode (invariance) at large radius" begin
        full = Mycelia.correct_observations(lgraph, [observed])
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, local_radius=100)
        loc = Mycelia.correct_observations(lgraph, [observed]; config=cfg)
        Test.@test decoded_strings(loc) == decoded_strings(full)
        Test.@test only(loc.paths).score ≈ only(full.paths).score
    end

    # Two-component fixture: the read's component plus a distant, unrelated
    # sequence sharing NO 3-mers. Localization must exclude the distant component
    # while returning the identical corrected read — the invariance that makes the
    # localized decode lossless while doing strictly less work.
    two = [
        FASTX.FASTA.Record("read_comp", BioSequences.dna"ATGCGTA"),
        FASTX.FASTA.Record("distant", BioSequences.dna"TTTAAAGGGCCC"),
    ]
    tgraph = Mycelia.Rhizomorph.build_kmer_graph(
        two, 3; dataset_id="local_twocomp", mode=:singlestrand)
    observed2 = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGC"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
        Kmers.DNAKmer{3}("GTA"),
    ]

    Test.@testset "local decode excludes distant component, same result as full" begin
        full = Mycelia.correct_observations(tgraph, [observed2])
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, local_radius=10)
        loc = Mycelia.correct_observations(tgraph, [observed2]; config=cfg)
        Test.@test decoded_strings(loc) == decoded_strings(full)
        Test.@test only(loc.paths).score ≈ only(full.paths).score
    end

    Test.@testset "local_radius must be non-negative" begin
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(local_radius=-1)
    end
end
