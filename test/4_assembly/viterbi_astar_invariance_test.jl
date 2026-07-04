# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_astar_invariance_test.jl")'
# ```
#
# Tier 1C (td-4jdi): exact admissible-bound A* decoder. The monotone-invariance
# guardrail — the returned corrected read + score must be IDENTICAL across
# heuristics h=0 (Dijkstra) -> gap, and equal to the exact unbounded-beam decode.
# This proves each stronger (still-admissible) heuristic only prunes more DP
# states while preserving the optimum (and, for :gap, that the emax bound is
# admissible). `:astar_pops` must be non-increasing as the heuristic strengthens.

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

decoded(result) = [string(l) for l in only(result.corrected_observations)]
pops(result) = get(only(result.paths).diagnostics, :astar_pops, -1)

Test.@testset "Viterbi exact A* decoder invariance (Tier 1C, td-4jdi)" begin
    # B8 substitution-recovery fixture: observed[2] = TGA is an error; the exact
    # decode must recover TGC (the graph path ATG->TGC->GCG->CGT).
    records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
    graph = Mycelia.Rhizomorph.build_kmer_graph(
        records, 3; dataset_id="astar_fix", mode=:singlestrand)
    observed = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGA"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
    ]

    # Exact reference: unbounded-beam decoder (the byte-identical B8 oracle).
    beam = Mycelia.correct_observations(graph, [observed])
    cfg_dij = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, decoder=:astar, heuristic=:dijkstra)
    dij = Mycelia.correct_observations(graph, [observed]; config=cfg_dij)
    cfg_gap = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, decoder=:astar, heuristic=:gap)
    gap = Mycelia.correct_observations(graph, [observed]; config=cfg_gap)

    Test.@testset "A* (Dijkstra) == exact beam, byte-identical" begin
        Test.@test decoded(dij) == ["ATG", "TGC", "GCG", "CGT"]
        Test.@test decoded(dij) == decoded(beam)
        Test.@test only(dij.paths).score ≈ only(beam.paths).score
    end

    Test.@testset "A* (gap) == A* (Dijkstra) — emax heuristic is admissible" begin
        Test.@test decoded(gap) == decoded(dij)
        Test.@test only(gap.paths).score ≈ only(dij.paths).score
    end

    Test.@testset "stronger heuristic only prunes more (pops non-increasing)" begin
        Test.@test pops(gap) <= pops(dij)
    end

    Test.@testset "tangle cap is honored and flagged" begin
        cfg_cap = Mycelia.ViterbiCorrectionConfig(
            alphabet=:DNA, decoder=:astar, heuristic=:gap, max_astar_pops=1)
        capped = Mycelia.correct_observations(graph, [observed]; config=cfg_cap)
        Test.@test get(only(capped.paths).diagnostics, :astar_cap_hit, false) == true
    end

    Test.@testset "decoder/heuristic validation" begin
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(decoder=:foo)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(heuristic=:foo)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(max_astar_pops=0)
    end
end
