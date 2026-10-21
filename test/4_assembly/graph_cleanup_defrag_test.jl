# ERROR-ONLY GRAPH-CLEANUP DEFRAGMENTATION (td-969e)
# ===================================================
#
# Unit + end-to-end coverage for the error-only defragmentation pass wired into
# the iterative corrector's re-assembly (clean_assembly_graph_next! and its two
# components clip_low_coverage_tips! / collapse_error_bubbles!). The load-bearing
# property is VARIATION PRESERVATION (td-h6w9): the cleanup removes only
# UNSUPPORTED (coverage-1) error structure and must never collapse a real,
# data-supported bubble branch. The variation_preservation_holdout_test.jl
# exercises the same invariant through the full pipeline; these tests pin the
# behavior at the graph-primitive level where the discriminator lives.

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import FASTX
import BioSequences
import Random

const _GCD = Mycelia.Rhizomorph

# Build a String-labeled qualmer-style graph (KmerVertexData/KmerEdgeData) whose
# per-vertex `count_evidence` coverage is exactly `cov` — `cov` distinct evidence
# entries. This mirrors how a real qualmer graph encodes coverage, so the cleanup
# primitives are tested against the same accessor (`count_evidence`) they use in
# production.
function _gcd_add_vertex!(graph, label::String, cov::Int)
    graph[label] = Mycelia.Rhizomorph.KmerVertexData(label)
    for i in 1:cov
        Mycelia.Rhizomorph.add_evidence!(graph[label], "ds", "obs$(i)",
            Mycelia.Rhizomorph.EvidenceEntry(i, Mycelia.Rhizomorph.Forward))
    end
    return graph
end

function _gcd_add_edge!(graph, src::String, dst::String)
    graph[src, dst] = Mycelia.Rhizomorph.KmerEdgeData()
    Mycelia.Rhizomorph.add_evidence!(graph[src, dst], "ds", "obs",
        Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Forward))
    return graph
end

function _gcd_new_graph()
    return MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = Mycelia.Rhizomorph.KmerVertexData,
        edge_data_type = Mycelia.Rhizomorph.KmerEdgeData
    )
end

Test.@testset "graph-cleanup defragmentation (td-969e)" begin
    Test.@testset "clip_low_coverage_tips! peels a coverage-1 dead-end spur (cascade)" begin
        g = _gcd_new_graph()
        # High-coverage backbone A->B->C->D (all 20x). A is a 20x source, D a 20x
        # sink — legitimate contig ends that must NOT be clipped.
        for (v, c) in (("A", 20), ("B", 20), ("C", 20), ("D", 20))
            _gcd_add_vertex!(g, v, c)
        end
        _gcd_add_edge!(g, "A", "B")
        _gcd_add_edge!(g, "B", "C")
        _gcd_add_edge!(g, "C", "D")
        # A two-vertex coverage-1 error spur hanging off C: C->S1->S2 (S2 dead-end).
        _gcd_add_vertex!(g, "S1", 1)
        _gcd_add_vertex!(g, "S2", 1)
        _gcd_add_edge!(g, "C", "S1")
        _gcd_add_edge!(g, "S1", "S2")
        # An isolated coverage-1 fragment (both-ends dead) — also error.
        _gcd_add_vertex!(g, "ISO", 1)

        removed = _GCD.clip_low_coverage_tips!(g; max_coverage = 1)

        Test.@test removed == 3                       # S1, S2, ISO
        Test.@test !MetaGraphsNext.haskey(g, "S1")
        Test.@test !MetaGraphsNext.haskey(g, "S2")
        Test.@test !MetaGraphsNext.haskey(g, "ISO")
        # The high-coverage backbone (incl. its 20x source/sink ends) is untouched.
        for v in ("A", "B", "C", "D")
            Test.@test MetaGraphsNext.haskey(g, v)
        end
    end

    Test.@testset "clip_low_coverage_tips! never clips a high-coverage dead end" begin
        g = _gcd_new_graph()
        _gcd_add_vertex!(g, "A", 20)
        _gcd_add_vertex!(g, "TIP", 5)   # coverage-5 dead-end: real, above threshold
        _gcd_add_edge!(g, "A", "TIP")

        removed = _GCD.clip_low_coverage_tips!(g; max_coverage = 1)
        Test.@test removed == 0
        Test.@test MetaGraphsNext.haskey(g, "TIP")
    end

    Test.@testset "clip_low_coverage_tips! leaves interior (reconverging) branches" begin
        # Interior bubble arms are NOT dead ends, so tip clipping never touches them
        # regardless of coverage — the structural reason a real het bubble is safe.
        g = _gcd_new_graph()
        _gcd_add_vertex!(g, "E", 30)
        _gcd_add_vertex!(g, "P", 15)    # arm interior — coverage-15 but reconverges
        _gcd_add_vertex!(g, "Q", 1)     # arm interior — coverage-1 but reconverges
        _gcd_add_vertex!(g, "X", 30)
        _gcd_add_edge!(g, "E", "P")
        _gcd_add_edge!(g, "E", "Q")
        _gcd_add_edge!(g, "P", "X")
        _gcd_add_edge!(g, "Q", "X")

        removed = _GCD.clip_low_coverage_tips!(g; max_coverage = 1)
        # Even coverage-1 Q survives tip clipping: it has in-degree 1 AND out-degree
        # 1 (not a dead end). Removing it is the BUBBLE collapse's job, under guard.
        Test.@test removed == 0
        Test.@test MetaGraphsNext.haskey(g, "Q")
    end

    Test.@testset "collapse_error_bubbles! removes a coverage-1 error arm" begin
        g = _gcd_new_graph()
        _gcd_add_vertex!(g, "E", 30)
        _gcd_add_vertex!(g, "REAL", 15)   # supported allele
        _gcd_add_vertex!(g, "ERR", 1)     # unsupported error arm
        _gcd_add_vertex!(g, "X", 30)
        _gcd_add_edge!(g, "E", "REAL")
        _gcd_add_edge!(g, "E", "ERR")
        _gcd_add_edge!(g, "REAL", "X")
        _gcd_add_edge!(g, "ERR", "X")

        arms = _GCD.collapse_error_bubbles!(g;
            error_max_coverage = 1, min_supported_coverage = 2)

        Test.@test arms == 1
        Test.@test !MetaGraphsNext.haskey(g, "ERR")   # error arm removed
        Test.@test MetaGraphsNext.haskey(g, "REAL")   # allele retained
        Test.@test MetaGraphsNext.haskey(g, "E")
        Test.@test MetaGraphsNext.haskey(g, "X")
    end

    Test.@testset "collapse_error_bubbles! PRESERVES a balanced (15x/15x) het bubble" begin
        # KEY INVARIANT (td-h6w9): both branches data-supported -> retained.
        g = _gcd_new_graph()
        _gcd_add_vertex!(g, "E", 30)
        _gcd_add_vertex!(g, "A", 15)
        _gcd_add_vertex!(g, "B", 15)
        _gcd_add_vertex!(g, "X", 30)
        _gcd_add_edge!(g, "E", "A")
        _gcd_add_edge!(g, "E", "B")
        _gcd_add_edge!(g, "A", "X")
        _gcd_add_edge!(g, "B", "X")

        arms = _GCD.collapse_error_bubbles!(g;
            error_max_coverage = 1, min_supported_coverage = 2)

        Test.@test arms == 0
        Test.@test MetaGraphsNext.haskey(g, "A")
        Test.@test MetaGraphsNext.haskey(g, "B")
    end

    Test.@testset "collapse_error_bubbles! PRESERVES a skewed-but-supported (20x/10x) bubble" begin
        # The skewed minority (10x) is real and far above the error threshold -> kept.
        g = _gcd_new_graph()
        _gcd_add_vertex!(g, "E", 30)
        _gcd_add_vertex!(g, "MAJ", 20)
        _gcd_add_vertex!(g, "MIN", 10)
        _gcd_add_vertex!(g, "X", 30)
        _gcd_add_edge!(g, "E", "MAJ")
        _gcd_add_edge!(g, "E", "MIN")
        _gcd_add_edge!(g, "MAJ", "X")
        _gcd_add_edge!(g, "MIN", "X")

        arms = _GCD.collapse_error_bubbles!(g;
            error_max_coverage = 1, min_supported_coverage = 2)

        Test.@test arms == 0
        Test.@test MetaGraphsNext.haskey(g, "MAJ")
        Test.@test MetaGraphsNext.haskey(g, "MIN")
    end

    Test.@testset "clean_assembly_graph_next! reports stats and combines both passes" begin
        g = _gcd_new_graph()
        # Backbone with a coverage-1 error tip AND a coverage-1 error bubble arm.
        for (v, c) in (("A", 20), ("B", 20), ("C", 20))
            _gcd_add_vertex!(g, v, c)
        end
        _gcd_add_edge!(g, "A", "B")
        _gcd_add_edge!(g, "B", "C")
        _gcd_add_vertex!(g, "TIP", 1)
        _gcd_add_edge!(g, "C", "TIP")
        # error bubble off B: B->{REAL(20)->X , ERR(1)->X}
        _gcd_add_vertex!(g, "REAL", 20)
        _gcd_add_vertex!(g, "ERR", 1)
        _gcd_add_vertex!(g, "X", 20)
        # Rewire: make B the entry of the bubble by adding B->REAL, B->ERR.
        _gcd_add_edge!(g, "B", "REAL")
        _gcd_add_edge!(g, "B", "ERR")
        _gcd_add_edge!(g, "REAL", "X")
        _gcd_add_edge!(g, "ERR", "X")

        stats = _GCD.clean_assembly_graph_next!(g)

        Test.@test stats.tips_removed >= 1
        Test.@test stats.bubble_arms_removed == 1
        Test.@test stats.vertices_after < stats.vertices_before
        Test.@test !MetaGraphsNext.haskey(g, "TIP")
        Test.@test !MetaGraphsNext.haskey(g, "ERR")
        Test.@test MetaGraphsNext.haskey(g, "REAL")
    end

    Test.@testset "clean_assembly_graph_next! is a no-op on a clean linear graph" begin
        g = _gcd_new_graph()
        for (v, c) in (("A", 20), ("B", 20), ("C", 20))
            _gcd_add_vertex!(g, v, c)
        end
        _gcd_add_edge!(g, "A", "B")
        _gcd_add_edge!(g, "B", "C")

        stats = _GCD.clean_assembly_graph_next!(g)
        Test.@test stats.tips_removed == 0
        Test.@test stats.bubble_arms_removed == 0
        Test.@test stats.vertices_before == stats.vertices_after == 3
    end

    Test.@testset "end-to-end: corrector applies cleanup; plain qualmer path does not" begin
        # Toy 1kb reads (150bp, ~20x, err=0.01) through the iterative corrector.
        rec = Mycelia.random_fasta_record(moltype = :DNA, seed = 7, L = 1000)
        refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
        rng = Random.MersenneTwister(7)
        reads = FASTX.FASTQ.Record[]
        for i in 1:ceil(Int, 20 * 1000 / 150)
            s = rand(rng, 1:(1000 - 150 + 1))
            frag = refseq[s:(s + 149)]
            rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
            obs, q = Mycelia.observe(frag; error_rate = 0.01, tech = :illumina)
            isempty(obs) && continue
            push!(reads, FASTX.FASTQ.Record("r$i", string(obs),
                String([Char(x + 33) for x in q])))
        end

        corrected = _GCD.assemble_genome(reads; k = 21,
            corrector = :iterative, strategy = :scalable)
        Test.@test corrected.assembly_stats["graph_cleanup_applied"] == true
        # frac must stay >= 1.0 (no real sequence lost) even as fragments drop.
        total = sum(length.(corrected.contigs); init = 0)
        Test.@test total >= 1000
        # The cleanup removed unsupported error structure (>=1 vertex).
        removed_total = corrected.assembly_stats["graph_cleanup_tips_removed"] +
                        corrected.assembly_stats["graph_cleanup_bubble_arms_removed"]
        Test.@test removed_total >= 1

        # The plain qualmer path (corrector=:none) is byte-for-byte unchanged: no
        # cleanup applied (graph_cleanup defaults OFF).
        plain = _GCD.assemble_genome(reads; k = 21)
        Test.@test plain.assembly_stats["graph_cleanup_applied"] == false
    end
end
