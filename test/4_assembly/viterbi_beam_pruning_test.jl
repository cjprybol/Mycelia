# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_beam_pruning_test.jl")'
# ```
#
# Regression for td-63qy: the exact Viterbi corrector (`_viterbi_correct_observation`)
# never beam-capped its (vertex, strand) frontier, so the reachable-state set grew
# ~unboundedly with read-length depth on real branchy graphs — 21 billion
# allocations and a hard crash on a 48 kb phage. `ViterbiCorrectionConfig` now
# carries `beam_width` (default `typemax(Int)` = exact/unbounded); a finite value
# caps the frontier to the top-K states by score each depth.

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

function beam_decoded_label_strings(
    result::Mycelia.ViterbiCorrectionResult
)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

Test.@testset "Viterbi corrector beam pruning (td-63qy)" begin
    # A simple linear graph: the frontier is small, so the exact and bounded
    # decoders must agree, and pruning must never fire.
    linear_records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
    linear_graph = Mycelia.Rhizomorph.build_kmer_graph(
        linear_records, 3; dataset_id="beam_linear", mode=:singlestrand
    )
    linear_observed = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGA"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
    ]

    Test.@testset "default beam_width is unbounded (exact, no regression)" begin
        result = Mycelia.correct_observations(linear_graph, [linear_observed])
        # The exact answer from the existing B8 fixture — must be byte-identical.
        Test.@test beam_decoded_label_strings(result) == ["ATG", "TGC", "GCG", "CGT"]
        # No pruning on the default (unbounded) path.
        Test.@test get(only(result.paths).diagnostics, :beam_pruned, 0) == 0
    end

    Test.@testset "a finite beam still recovers the linear correction" begin
        config = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, beam_width=1)
        result = Mycelia.correct_observations(linear_graph, [linear_observed]; config=config)
        # Anti-crash guarantee: a bounded corrector completes and returns a
        # valid, complete correction on a linear graph.
        Test.@test beam_decoded_label_strings(result) == ["ATG", "TGC", "GCG", "CGT"]
    end

    # A bubble graph (two paths sharing a start/end k-mer) forces the frontier
    # above one state at the branch, so a beam_width=1 corrector MUST prune —
    # and must still complete rather than explode.
    bubble_records = [
        FASTX.FASTA.Record("path_a", BioSequences.dna"ATGCGTA"),
        FASTX.FASTA.Record("path_b", BioSequences.dna"ATGAGTA"),
    ]
    bubble_graph = Mycelia.Rhizomorph.build_kmer_graph(
        bubble_records, 3; dataset_id="beam_bubble", mode=:singlestrand
    )
    bubble_observed = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGC"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
        Kmers.DNAKmer{3}("GTA"),
    ]

    Test.@testset "beam pruning activates on a branch and still completes" begin
        config = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, beam_width=1)
        result = Mycelia.correct_observations(bubble_graph, [bubble_observed]; config=config)
        decoded = only(result.corrected_observations)
        # Completes with a non-empty correction (did not explode / crash).
        Test.@test !isempty(decoded)
        # Pruning fired at least once (the branch pushed the frontier past 1).
        Test.@test get(only(result.paths).diagnostics, :beam_pruned, 0) >= 1
    end

    Test.@testset "unbounded decoder on the bubble is at least as good as beamed" begin
        exact = Mycelia.correct_observations(bubble_graph, [bubble_observed])
        beamed_config = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, beam_width=1)
        beamed = Mycelia.correct_observations(bubble_graph, [bubble_observed]; config=beamed_config)
        # Beam search is an approximation: the exact score is an upper bound.
        Test.@test only(exact.paths).score >= only(beamed.paths).score - 1e-9
    end

    Test.@testset "beam_width must be positive" begin
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_width=0)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_width=-4)
    end

    # Size-aware auto-beam (td-63qy): the SHIPPING corrector default must stay
    # exact on small reads (preserving the ML guarantee) yet bound the frontier on
    # large reads so it cannot OOM-crash (~21B allocations on a 48 kb phage).
    Test.@testset "size-aware auto-beam default (_auto_beam_width)" begin
        # Small read (observation count at/below the threshold): stays EXACT.
        Test.@test Mycelia._auto_beam_width(1) == typemax(Int)
        Test.@test Mycelia._auto_beam_width(300) == typemax(Int)
        Test.@test Mycelia._auto_beam_width(Mycelia._AUTO_BEAM_EXACT_THRESHOLD) ==
                   typemax(Int)
        # Large read (above the threshold): BOUNDED to the proven-tractable width.
        Test.@test Mycelia._auto_beam_width(Mycelia._AUTO_BEAM_EXACT_THRESHOLD + 1) ==
                   Mycelia._AUTO_BEAM_BOUNDED_WIDTH
        Test.@test Mycelia._auto_beam_width(48_000) == Mycelia._AUTO_BEAM_BOUNDED_WIDTH
        # The bounded width is a finite, positive, previously-proven value (256).
        Test.@test Mycelia._AUTO_BEAM_BOUNDED_WIDTH == 256
        Test.@test 0 < Mycelia._AUTO_BEAM_BOUNDED_WIDTH < typemax(Int)
    end

    # Graph-density-aware auto-beam (td-35ux): the read-length-only rule left short
    # reads (n_obs <= threshold) EXACT even on a dense intermediate k-rung, where
    # the exact retained frontier grows O(n_vertices) per decode depth →
    # O(genome^2) per pass (~15,723 frontier states / ~2.15 s for one 150 bp read
    # at 5 kb / k=9). The two-arg form bounds the beam when the graph is dense,
    # regardless of read length.
    Test.@testset "graph-density-aware auto-beam (_auto_beam_width, td-35ux)" begin
        bound = Mycelia._AUTO_BEAM_BOUNDED_WIDTH
        # Short read (130 observations, below the exact threshold) on a DENSE graph
        # (2230 vertices > 256): bounded despite the small read length. This is the
        # #376 short-read blowup case — the exact frontier would grow O(n_vertices).
        Test.@test Mycelia._auto_beam_width(130, 2230) == bound
        # Short read on a SPARSE graph (few vertices): stays exact — the exact
        # frontier is tiny, so the ML guarantee is preserved at no tractability cost.
        Test.@test Mycelia._auto_beam_width(130, 10) == typemax(Int)
        Test.@test Mycelia._auto_beam_width(130, bound) == typemax(Int)
        # Boundary: one vertex past the bounded width flips exact → bounded.
        Test.@test Mycelia._auto_beam_width(130, bound + 1) == bound
        # Large read stays bounded no matter how sparse the graph (read-length rule
        # still dominates — a large read is never exact regardless of density).
        Test.@test Mycelia._auto_beam_width(
            Mycelia._AUTO_BEAM_EXACT_THRESHOLD + 1, 10) == bound
        # A large read on a dense graph is likewise bounded.
        Test.@test Mycelia._auto_beam_width(48_000, 5000) == bound
    end

    # An explicit beam_width override is honored verbatim by the per-read decoder
    # entry point: passing typemax(Int) forces EXACT decoding regardless of the
    # auto-default, and the correction still completes on a small read.
    Test.@testset "explicit beam_width override forces exact on the decoder entry" begin
        records = [FASTX.FASTQ.Record("r", "ATGCGT", repeat("I", 6))]
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; dataset_id="beam_override", mode=:singlestrand)
        read = FASTX.FASTQ.Record("q", "ATGCGT", repeat("I", 6))
        # Explicit exact override.
        forced = Mycelia.try_viterbi_path_improvement(
            read, graph, 3; graph_mode=:singlestrand, beam_width=typemax(Int))
        Test.@test forced isa Union{Tuple{FASTX.FASTQ.Record, Float64}, Nothing}
        # Auto default (nothing) also completes on this small read.
        auto = Mycelia.try_viterbi_path_improvement(read, graph, 3; graph_mode=:singlestrand)
        Test.@test auto isa Union{Tuple{FASTX.FASTQ.Record, Float64}, Nothing}
    end
end

# Candidate-GENERATION bounds (td-plqi): the width beam caps the RETAINED frontier,
# but on a dense intermediate-k graph the generating frontier (the states that
# enumerate successors + emission-score them each depth) still climbs toward the
# width cap as the graph densifies — the empirically-measured residual super-linear
# decode term (#386). Two additive bounds cut generation directly:
#   * `max_successors_per_state` (top-B) — a per-state successor cap (no-op on DNA,
#     out-degree <= 4; a guard for pathological high-branching).
#   * `beam_score_margin` (Δ) — a score-threshold ("histogram") beam that keeps only
#     the near-best band, which does NOT grow with genome, so the generating
#     frontier stays O(1) in size while the improbable tail is discarded.
# Both default to a strict no-op (exact ML) and only engage where the width beam is
# already finite (approximate), so exact-ML reads stay byte-identical.
Test.@testset "Viterbi corrector candidate-generation bounds (td-plqi)" begin
    linear_records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
    linear_graph = Mycelia.Rhizomorph.build_kmer_graph(
        linear_records, 3; dataset_id="candgen_linear", mode=:singlestrand
    )
    linear_observed = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGA"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
    ]

    bubble_records = [
        FASTX.FASTA.Record("path_a", BioSequences.dna"ATGCGTA"),
        FASTX.FASTA.Record("path_b", BioSequences.dna"ATGAGTA"),
    ]
    bubble_graph = Mycelia.Rhizomorph.build_kmer_graph(
        bubble_records, 3; dataset_id="candgen_bubble", mode=:singlestrand
    )
    bubble_observed = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGC"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
        Kmers.DNAKmer{3}("GTA"),
    ]

    Test.@testset "constructor validation" begin
        # A positive successor cap and a positive-or-Inf score margin are required.
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(max_successors_per_state=0)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(max_successors_per_state=-3)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_score_margin=0.0)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_score_margin=-5.0)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_score_margin=NaN)
        # Inf (the default) is valid — it means "no threshold" (exact).
        Test.@test Mycelia.ViterbiCorrectionConfig(beam_score_margin=Inf).beam_score_margin == Inf
    end

    Test.@testset "defaults are exact (both bounds are no-ops)" begin
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA)
        Test.@test cfg.max_successors_per_state == typemax(Int)
        Test.@test cfg.beam_score_margin == Inf
        # Byte-identical to the established B8 exact fixture, and neither bound fires.
        result = Mycelia.correct_observations(linear_graph, [linear_observed])
        Test.@test beam_decoded_label_strings(result) == ["ATG", "TGC", "GCG", "CGT"]
        diag = only(result.paths).diagnostics
        Test.@test get(diag, :successor_bounded, 0) == 0
        Test.@test get(diag, :margin_pruned, 0) == 0
    end

    Test.@testset "top-B successor bound: B >= out-degree is byte-identical (no-op)" begin
        # On a DNA de Bruijn graph the out-degree is <= 4, so B = 16 can never bite:
        # the corrected path must match the unbounded exact answer exactly.
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, max_successors_per_state=16)
        result = Mycelia.correct_observations(bubble_graph, [bubble_observed]; config=cfg)
        exact = Mycelia.correct_observations(bubble_graph, [bubble_observed])
        Test.@test beam_decoded_label_strings(result) == beam_decoded_label_strings(exact)
        Test.@test get(only(result.paths).diagnostics, :successor_bounded, 0) == 0
    end

    Test.@testset "top-B successor bound: a tight B engages on a branch and completes" begin
        # B = 1 forces each expanded state to enumerate only its single top-weight
        # successor. The branch vertex has 2 successors, so the bound MUST fire, and
        # the decode still completes with a valid, non-empty correction.
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, max_successors_per_state=1)
        result = Mycelia.correct_observations(bubble_graph, [bubble_observed]; config=cfg)
        decoded = only(result.corrected_observations)
        Test.@test !isempty(decoded)
        Test.@test get(only(result.paths).diagnostics, :successor_bounded, 0) >= 1
        # Bounded generation is an approximation: its score cannot exceed the exact.
        exact = Mycelia.correct_observations(bubble_graph, [bubble_observed])
        Test.@test only(exact.paths).score >= only(result.paths).score - 1e-9
    end

    Test.@testset "score-margin bound: generous Δ keeps the exact best path" begin
        # A wide margin drops only states far below the best; on a small graph the
        # exact best path is preserved (its states are always within the band).
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, beam_score_margin=1000.0)
        result = Mycelia.correct_observations(bubble_graph, [bubble_observed]; config=cfg)
        exact = Mycelia.correct_observations(bubble_graph, [bubble_observed])
        Test.@test beam_decoded_label_strings(result) == beam_decoded_label_strings(exact)
    end

    Test.@testset "score-margin bound: a tight Δ prunes and never beats exact" begin
        # A tight margin discards the improbable frontier tail (margin_pruned fires),
        # completes, and — as an approximation — never scores above exact.
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, beam_score_margin=1.0)
        result = Mycelia.correct_observations(bubble_graph, [bubble_observed]; config=cfg)
        Test.@test !isempty(only(result.corrected_observations))
        Test.@test get(only(result.paths).diagnostics, :margin_pruned, 0) >= 1
        exact = Mycelia.correct_observations(bubble_graph, [bubble_observed])
        Test.@test only(exact.paths).score >= only(result.paths).score - 1e-9
    end
end
