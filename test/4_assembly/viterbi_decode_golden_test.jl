# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_decode_golden_test.jl")'
# ```
#
# Golden-master (characterization) lock for opt2's fused-typed-transition refactor
# of the substitution hot decode loop in `_viterbi_correct_observation`. Captured
# against PRE-refactor code (commit 047142aa8), before `_get_valid_transitions` +
# `_total_outgoing_weight` were replaced by the fused
# `Rhizomorph._valid_transitions_and_total` and Dict-form `transition[:x]` access
# was converted to `transition.x` field access on `_Transition`. Every literal
# below (path steps, score, diagnostics) is the exact pre-refactor output; any
# byte-level drift after the refactor fails this test.
#
# Fixture: a 3-read k=3 kmer graph where reads "ATGCGTAA" (x2) and "ATGCTTAA" (x1)
# share a prefix (ATG,TGC), branch at TGC into two out-edges (GCG coverage=2, GCT
# coverage=1), then reconverge at TAA. `max_successors_per_state = 1` forces
# `_top_b_transitions` to bound the 2-way branch to a single successor each depth
# (diagnostics[:successor_bounded] == 1 confirms it fires), so this fixture
# exercises both changed call sites: the fused `_valid_transitions_and_total` scan
# and the `_top_b_transitions` field-access sort.

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

Test.@testset "Viterbi decode golden master (opt2 fused-transitions lock)" begin
    records = [
        FASTX.FASTA.Record("r1", BioSequences.dna"ATGCGTAA"),
        FASTX.FASTA.Record("r2", BioSequences.dna"ATGCGTAA"),
        FASTX.FASTA.Record("r3", BioSequences.dna"ATGCTTAA")
    ]
    graph = Mycelia.Rhizomorph.build_kmer_graph(
        records, 3; dataset_id = "golden_fixture", mode = :singlestrand
    )

    observation = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGC"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
        Kmers.DNAKmer{3}("GTA"),
        Kmers.DNAKmer{3}("TAA")
    ]

    config = Mycelia.ViterbiCorrectionConfig(max_successors_per_state = 1)

    alphabet = Mycelia._resolve_viterbi_alphabet(graph, [observation], config.alphabet)
    strand_mode = Mycelia._resolve_viterbi_strand_mode(graph, alphabet, config.strand_mode)
    transition_edge_weight = Mycelia._viterbi_transition_edge_weight(
        graph, config.edge_weight)
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(
        graph; edge_weight = transition_edge_weight)

    result = Mycelia._viterbi_correct_observation(
        weighted,
        observation,
        alphabet;
        config = config,
        strand_mode = strand_mode,
        quality_graph = graph,
        transition_scoring = Mycelia._viterbi_transition_scoring(
            graph, transition_edge_weight)
    )

    Test.@test result.path !== nothing
    path = something(result.path)

    expected_labels = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGC"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
        Kmers.DNAKmer{3}("GTA"),
        Kmers.DNAKmer{3}("TAA")
    ]
    expected_strands = fill(Mycelia.Rhizomorph.Forward, 6)
    expected_probabilities = [1.0, 1.0, 2.0 / 3.0, 1.0, 1.0, 1.0]
    expected_cumulative_probabilities = [
        1.0, 1.0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]

    Test.@test [step.vertex_label for step in path.steps] == expected_labels
    Test.@test [step.strand for step in path.steps] == expected_strands
    Test.@test [step.probability for step in path.steps] == expected_probabilities
    Test.@test [step.cumulative_probability for step in path.steps] ==
               expected_cumulative_probabilities
    Test.@test path.total_probability === 0.6666666666666666

    Test.@test result.score === -0.5863711534711905

    expected_diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_emission_correct_observation,
        :exact => true,
        :alphabet => :DNA,
        :strand_mode => :singlestrand,
        :reverse_complement_support => true,
        :max_steps => 5,
        :target_vertex => Kmers.DNAKmer{3}("TAA"),
        :start_strand => Mycelia.Rhizomorph.Forward,
        :score_domain => :log_probability,
        :transition_scoring => :normalized_edge_weight,
        :emission_scoring => :alphabet_parameterized,
        :expanded_states => 5,
        :generated_states => 5,
        :retained_states => 1,
        :cumulative_retained_states => 6,
        :max_retained_states => 1,
        :skipped_transitions => 0,
        :completed_steps => 5,
        :reached_target => true,
        :path_length => 6,
        :successor_bounded => 1
    )
    Test.@test length(result.diagnostics) == 21
    Test.@test result.diagnostics == expected_diagnostics
end

# ---------------------------------------------------------------------------
# Second golden case (opt2 diagnostics-hoist lock, TRIGGERED side): the fixture
# above never prunes (default beam_width = typemax, margin = Inf), so it only
# tests the ABSENCE of :beam_pruned / :margin_pruned. This case uses a small
# beam (beam_width = 2) and a finite score margin (beam_score_margin = 2.0) on a
# higher-branching 40-read k=5 graph so both counters actually fire, locking the
# hoisted-counter write-back path. Every literal is the exact opt2-branch output
# AND was verified byte-identical against origin/master (0804d0ecc) for this same
# fixture, proving the counter hoist did not change the pruning counts.
Test.@testset "Viterbi decode golden master (opt2 diagnostics-hoist, pruning fires)" begin
    reads = [
        FASTX.FASTA.Record("r1", BioSequences.LongDNA{4}("GCCAACTACGTTTACCATCTTTTG")),
        FASTX.FASTA.Record("r2", BioSequences.LongDNA{4}("AGCTACGTTTACCATCTTTTGAGA")),
        FASTX.FASTA.Record("r3", BioSequences.LongDNA{4}("GAGAGCTGGTAGACTACTCTTCGG")),
        FASTX.FASTA.Record("r4", BioSequences.LongDNA{4}("CCAGGTCCGTTTACCATCTATTGA")),
        FASTX.FASTA.Record("r5", BioSequences.LongDNA{4}("GCTACGATTACCATCTTTTGAGAT")),
        FASTX.FASTA.Record("r6", BioSequences.LongDNA{4}("TGCCAGCTACGTTTAACATCTTTT")),
        FASTX.FASTA.Record("r7", BioSequences.LongDNA{4}("TACGTTTACCATCCTTTGAGAGCT")),
        FASTX.FASTA.Record("r8", BioSequences.LongDNA{4}("GGGCTTGCTGCAAGCTACGTTTAC")),
        FASTX.FASTA.Record("r9", BioSequences.LongDNA{4}("TGAGATCTGGTAGACTACTAATCA")),
        FASTX.FASTA.Record("r10", BioSequences.LongDNA{4}("TCCGTGCTGCCAGCTACGTTTAGC")),
        FASTX.FASTA.Record("r11", BioSequences.LongDNA{4}("CCAGCTACGTGTACCAACTTTTGA")),
        FASTX.FASTA.Record("r12", BioSequences.LongDNA{4}("GCTTTCTGCCAGCTACGTTTACCA")),
        FASTX.FASTA.Record("r13", BioSequences.LongDNA{4}("GCTCGCTGCCAGCTACGTTTACCA")),
        FASTX.FASTA.Record("r14", BioSequences.LongDNA{4}("CTGCCAGCTAAGTTTACCATCTTT")),
        FASTX.FASTA.Record("r15", BioSequences.LongDNA{4}("CTTGCTGAAAGCTACGTTTACCAT")),
        FASTX.FASTA.Record("r16", BioSequences.LongDNA{4}("TATCTTTTGAGATCTGGTAGACTA")),
        FASTX.FASTA.Record("r17", BioSequences.LongDNA{4}("CGTTTACCATCTTTTGAGATCTGG")),
        FASTX.FASTA.Record("r18", BioSequences.LongDNA{4}("GTTTACCATCTTTTGAGCGCTGGT")),
        FASTX.FASTA.Record("r19", BioSequences.LongDNA{4}("ATACGTTTACCATCTTTTGAGATC")),
        FASTX.FASTA.Record("r20", BioSequences.LongDNA{4}("GGGCTAGCTGCCAGCTACGTTTAA")),
        FASTX.FASTA.Record("r21", BioSequences.LongDNA{4}("AGCTACGTTTACCATCTTTTGAGA")),
        FASTX.FASTA.Record("r22", BioSequences.LongDNA{4}("CCTTTTGAGATCTGGTAGACTACT")),
        FASTX.FASTA.Record("r23", BioSequences.LongDNA{4}("GTCTTACTCCCAGCTACGTTGACC")),
        FASTX.FASTA.Record("r24", BioSequences.LongDNA{4}("TTTACCATCTTTTGACATCTGGTA")),
        FASTX.FASTA.Record("r25", BioSequences.LongDNA{4}("TACCATCTTTTTACGTCTGGTAGA")),
        FASTX.FASTA.Record("r26", BioSequences.LongDNA{4}("GGCTTGCTGCCCGCTATGTTTACC")),
        FASTX.FASTA.Record("r27", BioSequences.LongDNA{4}("ACCATCTTTTGAGATCTGGTAAAC")),
        FASTX.FASTA.Record("r28", BioSequences.LongDNA{4}("AGCTACCTTTACCATCTTTTGAGA")),
        FASTX.FASTA.Record("r29", BioSequences.LongDNA{4}("GCTTGCTGCCACCTACGTTTACCA")),
        FASTX.FASTA.Record("r30", BioSequences.LongDNA{4}("CTGCCAGCTATGATTACCATCTTT")),
        FASTX.FASTA.Record("r31", BioSequences.LongDNA{4}("GGGCTTCCTGCCAGCTACCTTTAC")),
        FASTX.FASTA.Record("r32", BioSequences.LongDNA{4}("CTTTTGAGATCTGGTAGACTACTA")),
        FASTX.FASTA.Record("r33", BioSequences.LongDNA{4}("ACCATCTTTCGAGATCTGGAAGAC")),
        FASTX.FASTA.Record("r34", BioSequences.LongDNA{4}("GAGATCTGGTAGACTACTAATCGG")),
        FASTX.FASTA.Record("r35", BioSequences.LongDNA{4}("TTGCTGCCAGATACGTTTACCATC")),
        FASTX.FASTA.Record("r36", BioSequences.LongDNA{4}("CTACGTTTACCAGCTTTTGAGATC")),
        FASTX.FASTA.Record("r37", BioSequences.LongDNA{4}("TGCCAGCTACGTTGACTATCTTTT")),
        FASTX.FASTA.Record("r38", BioSequences.LongDNA{4}("GCCGCGAGCTACGTTTACCTTCTT")),
        FASTX.FASTA.Record("r39", BioSequences.LongDNA{4}("CTTGCCGCCAGCTACGTTTACCAT")),
        FASTX.FASTA.Record("r40", BioSequences.LongDNA{4}("TTGCGGCCAGCTACGTTTACCATC"))
    ]
    graph = Mycelia.Rhizomorph.build_kmer_graph(
        reads, 5; dataset_id = "golden_prune_fixture", mode = :singlestrand
    )

    obsseq = "GCCAACTACGTTTACCATCTTTTG"
    observation = [Kmers.DNAKmer{5}(obsseq[i:(i + 4)]) for i in 1:(length(obsseq) - 4)]

    config = Mycelia.ViterbiCorrectionConfig(
        beam_width = 2, beam_score_margin = 2.0)

    alphabet = Mycelia._resolve_viterbi_alphabet(graph, [observation], config.alphabet)
    strand_mode = Mycelia._resolve_viterbi_strand_mode(graph, alphabet, config.strand_mode)
    transition_edge_weight = Mycelia._viterbi_transition_edge_weight(
        graph, config.edge_weight)
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(
        graph; edge_weight = transition_edge_weight)

    result = Mycelia._viterbi_correct_observation(
        weighted,
        observation,
        alphabet;
        config = config,
        strand_mode = strand_mode,
        quality_graph = graph,
        transition_scoring = Mycelia._viterbi_transition_scoring(
            graph, transition_edge_weight)
    )

    Test.@test result.path !== nothing
    path = something(result.path)

    expected_labels = [
        Kmers.DNAKmer{5}("GCCAA"),
        Kmers.DNAKmer{5}("CCAAC"),
        Kmers.DNAKmer{5}("CAACT"),
        Kmers.DNAKmer{5}("AACTA"),
        Kmers.DNAKmer{5}("ACTAC"),
        Kmers.DNAKmer{5}("CTACG"),
        Kmers.DNAKmer{5}("TACGT"),
        Kmers.DNAKmer{5}("ACGTT"),
        Kmers.DNAKmer{5}("CGTTT"),
        Kmers.DNAKmer{5}("GTTTA"),
        Kmers.DNAKmer{5}("TTTAC"),
        Kmers.DNAKmer{5}("TTACC"),
        Kmers.DNAKmer{5}("TACCA"),
        Kmers.DNAKmer{5}("ACCAT"),
        Kmers.DNAKmer{5}("CCATC"),
        Kmers.DNAKmer{5}("CATCT"),
        Kmers.DNAKmer{5}("ATCTT"),
        Kmers.DNAKmer{5}("TCTTT"),
        Kmers.DNAKmer{5}("CTTTT"),
        Kmers.DNAKmer{5}("TTTTG")
    ]
    expected_strands = fill(Mycelia.Rhizomorph.Forward, 20)
    expected_probabilities = [1.0, 1.0, 1.0, 0.5, 1.0, 0.16666666666666666, 0.9473684210526315, 0.9090909090909091, 0.9, 1.0, 0.8695652173913043, 0.9545454545454546, 0.9545454545454546, 0.9, 1.0, 0.9375, 0.8823529411764706, 1.0, 0.9333333333333333, 0.9375]
    expected_cumulative_probabilities = [1.0, 1.0, 1.0, 0.5, 0.5, 0.08333333333333333, 0.07894736842105263, 0.07177033492822966, 0.0645933014354067, 0.0645933014354067, 0.056168088204701476, 0.05361499328630596, 0.051177948136928414, 0.04606015332323558, 0.04606015332323558, 0.043181393740533355, 0.03810122977105884, 0.03810122977105884, 0.035561147786321586, 0.03333857604967649]

    Test.@test [step.vertex_label for step in path.steps] == expected_labels
    Test.@test [step.strand for step in path.steps] == expected_strands
    Test.@test [step.probability for step in path.steps] == expected_probabilities
    Test.@test [step.cumulative_probability for step in path.steps] ==
               expected_cumulative_probabilities
    Test.@test path.total_probability === 0.03333857604967649

    Test.@test result.score === -4.406073697889444

    expected_diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_emission_correct_observation,
        :exact => true,
        :alphabet => :DNA,
        :strand_mode => :singlestrand,
        :reverse_complement_support => true,
        :max_steps => 19,
        :target_vertex => Kmers.DNAKmer{5}("TTTTG"),
        :start_strand => Mycelia.Rhizomorph.Forward,
        :score_domain => :log_probability,
        :transition_scoring => :normalized_edge_weight,
        :emission_scoring => :alphabet_parameterized,
        :expanded_states => 19,
        :generated_states => 36,
        :retained_states => 1,
        :cumulative_retained_states => 20,
        :max_retained_states => 1,
        :skipped_transitions => 0,
        :completed_steps => 19,
        :reached_target => true,
        :path_length => 20,
        :beam_pruned => 4,
        :margin_pruned => 13
    )
    # 22 keys: the 20 always-present keys plus the two triggered counters
    # (:beam_pruned, :margin_pruned). :successor_bounded stays ABSENT (default
    # max_successors_per_state = typemax never fires).
    Test.@test length(result.diagnostics) == 22
    Test.@test result.diagnostics == expected_diagnostics
end
