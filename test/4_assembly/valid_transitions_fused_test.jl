# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/valid_transitions_fused_test.jl")'
# ```

# FUSED VALID-TRANSITIONS EQUIVALENCE (opt2)
# _valid_transitions_and_total must return, per (vertex,strand): the SAME
# transitions in the SAME order as the Dict-form _get_valid_transitions, and a
# total_out bit-identical to _total_outgoing_weight. This is the byte-identity
# lock for the fused hot-loop function.
import Test
import Mycelia
import FASTX
import Random

Test.@testset "fused valid-transitions equivalence (opt2)" begin
    R = Mycelia.Rhizomorph
    rng = Random.MersenneTwister(31337)
    ref = join(rand(rng, ['A', 'C', 'G', 'T'], 400))
    reads = [FASTX.FASTQ.Record("r$i", ref[s:(s + 79)], String(fill('I', 80)))
             for (i, s) in enumerate(rand(rng, 1:321, 60))]
    k = 11
    raw_graph = R.build_qualmer_graph(reads, k; mode = :canonical)
    # _get_valid_transitions/_total_outgoing_weight/_valid_transitions_and_total
    # all read edge_data.src_strand/dst_strand, which only StrandWeightedEdgeData
    # (the directed weighted-graph form every production walk caller uses) has;
    # the raw qualmer graph's QualmerEdgeData does not carry those fields.
    graph = R.weighted_graph_from_rhizomorph(raw_graph)
    checked = 0
    for label in Mycelia.MetaGraphsNext.labels(graph)
        for strand in (R.Forward, R.Reverse)
            dict_tx = R._get_valid_transitions(graph, label, strand)
            dict_total = R._total_outgoing_weight(graph, label, strand)
            fused_tx, fused_total = R._valid_transitions_and_total(graph, label, strand)
            # fix round 1: the returned vector's static eltype must be concrete
            # (Vector{_Transition{L,E}}), not the abstract Vector{_Transition} that
            # `_Transition[]` produces — checked every iteration (empty and non-empty
            # alike), since eltype concreteness is a static property of the vector.
            Test.@test isconcretetype(eltype(fused_tx))
            Test.@test fused_total === dict_total                       # bit-identical
            Test.@test length(fused_tx) == length(dict_tx)
            for (a, b) in zip(fused_tx, dict_tx)
                Test.@test a.target_vertex == b[:target_vertex]
                Test.@test a.target_strand == b[:target_strand]
                Test.@test R._edge_transition_weight(a.edge_data) === b[:probability]  # bit-identical
                Test.@test a.edge_data === b[:edge_data]
            end
            checked += 1
        end
    end
    Test.@test checked > 0                                             # non-vacuous
end
