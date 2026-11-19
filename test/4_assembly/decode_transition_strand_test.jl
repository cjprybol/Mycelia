# Double-strand edges store a directed observation and its flipped
# reverse-complement image in one evidence set. The decoder must keep that
# summed edge on the orientation the read was observed in. `first_evidence_strand`
# follows set iteration order and flips along a read, which drops the Viterbi
# walk mid-read.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/decode_transition_strand_test.jl")'

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import Kmers
import Random

function _edge_graph(edge)
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = Nothing,
        edge_data_type = typeof(edge),
    )
    graph["A"] = nothing
    graph["B"] = nothing
    graph["A", "B"] = edge
    return graph
end

Test.@testset "decode transition strand on summed double-strand edges" begin
    R = Mycelia.Rhizomorph

    Test.@testset "Forward evidence stays walkable when the set yields Reverse first" begin
        caught_reverse_first = false
        for pos in 1:64
            edge = R.KmerEdgeData()
            R.add_evidence!(
                edge, "d", "o", R.EdgeEvidenceEntry(pos, pos + 1, R.Reverse))
            R.add_evidence!(
                edge, "d", "o", R.EdgeEvidenceEntry(pos + 2, pos + 3, R.Forward))
            R.first_evidence_strand(edge.evidence) == R.Reverse || continue
            weighted = R.weighted_graph_from_rhizomorph(_edge_graph(edge))
            weighted_edge = weighted["A", "B"]
            Test.@test weighted_edge.src_strand == R.Forward
            Test.@test weighted_edge.dst_strand == R.Forward
            Test.@test weighted_edge.weight == 2.0
            transitions, total_out = R._valid_transitions_and_total(
                weighted, "A", R.Forward)
            Test.@test length(transitions) == 1
            Test.@test only(transitions).target_vertex == "B"
            Test.@test only(transitions).target_strand == R.Forward
            Test.@test total_out == 2.0
            reverse_transitions, _ = R._valid_transitions_and_total(
                weighted, "A", R.Reverse)
            Test.@test isempty(reverse_transitions)
            caught_reverse_first = true
            break
        end
        Test.@test caught_reverse_first
    end

    Test.@testset "reverse-complement-only edges stay Reverse" begin
        edge = R.KmerEdgeData()
        R.add_evidence!(edge, "d", "o", R.EdgeEvidenceEntry(1, 2, R.Reverse))
        weighted = R.weighted_graph_from_rhizomorph(_edge_graph(edge))
        weighted_edge = weighted["A", "B"]
        Test.@test weighted_edge.src_strand == R.Reverse
        Test.@test weighted_edge.dst_strand == R.Reverse
        Test.@test weighted_edge.weight == 1.0
        transitions, _ = R._valid_transitions_and_total(weighted, "A", R.Reverse)
        Test.@test length(transitions) == 1
        Test.@test only(transitions).target_strand == R.Reverse
    end

    Test.@testset "partial reads on both strands stay Forward-connected" begin
        rng = Random.MersenneTwister(1)
        ref = BioSequences.randdnaseq(rng, 48)
        k = 7
        readlen = 24
        records = FASTX.FASTA.Record[]
        for i in 1:8
            start = rand(rng, 1:(length(ref) - readlen + 1))
            frag = ref[start:(start + readlen - 1)]
            push!(records, FASTX.FASTA.Record("f$(i)", string(frag)))
            push!(records, FASTX.FASTA.Record(
                "r$(i)", string(BioSequences.reverse_complement(frag))))
        end
        graph = R.build_kmer_graph(records, k; mode = :doublestrand)
        weighted = R.weighted_graph_from_rhizomorph(graph)
        saw_reverse_first = false
        for record in records
            seq = FASTX.sequence(BioSequences.LongDNA{4}, record)
            kmers = collect(Kmers.UnambiguousDNAMers{k}(seq))
            for index in 1:(length(kmers) - 1)
                left, left_pos = kmers[index]
                right, right_pos = kmers[index + 1]
                right_pos == left_pos + 1 || continue
                Test.@test haskey(graph, left, right)
                raw_edge = graph[left, right]
                if R.first_evidence_strand(raw_edge.evidence) == R.Reverse
                    saw_reverse_first = true
                end
                weighted_edge = weighted[left, right]
                Test.@test weighted_edge.src_strand == R.Forward
                Test.@test weighted_edge.weight ==
                           Float64(R.count_evidence(raw_edge))
                transitions, _ = R._valid_transitions_and_total(
                    weighted, left, R.Forward)
                Test.@test any(
                    transition -> transition.target_vertex == right,
                    transitions)
            end
        end
        Test.@test saw_reverse_first
    end
end
