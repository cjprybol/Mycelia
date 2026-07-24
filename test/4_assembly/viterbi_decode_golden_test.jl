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
    Test.@test result.diagnostics == expected_diagnostics
end
