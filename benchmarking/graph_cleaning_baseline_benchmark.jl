# Heuristic Graph Cleaning — COMPARISON BASELINE benchmark
# ==========================================================
#
# WHAT THIS IS (and, deliberately, is NOT):
#
# This script measures how much assembly-graph fragmentation a classical,
# explicit HEURISTIC cleaning pass removes: dead-end tip pruning
# (`remove_tips!`) plus bubble popping (`detect_bubbles_next` +
# `simplify_graph_next`). It reports contig count and a genome-fraction proxy
# (largest contig / reference length) WITH vs WITHOUT that cleaning on a small,
# deterministic toy assembly graph.
#
# It is a COMPARISON BASELINE ONLY. Per the assembly ADR, heuristic tip/bubble
# removal is NOT the primary cleaning mechanism in Mycelia — the intended
# mechanism is the emergent, soft-EM (probabilistic, evidence-reweighting)
# cleaning that falls out of the iterative assembly loop. The heuristic pass
# exercised here exists so that emergent soft-EM cleaning has a concrete,
# quantified yardstick to be measured against ("how much fragmentation does the
# cheap heuristic remove?"). Nothing in this file is wired into any default
# assembly path, and it must NOT be: the functions it calls
# (`remove_tips!`, `detect_bubbles_next`, `simplify_graph_next`) remain unwired
# from `assemble_genome` on purpose.
#
# The cleaning here is strictly OPT-IN: it only runs when a caller explicitly
# asks for it (the `clean` keyword of `assemble_toy_contigs`, default `false`).
#
# Run from the Mycelia base directory:
#
# ```bash
# LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
#     benchmarking/graph_cleaning_baseline_benchmark.jl
# ```

import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext
import Graphs
import Random

"""
    build_toy_reads(; seed=42) -> (reference, reads)

Construct a small, deterministic toy read set over a single reference. The read
set is engineered to inject the two fragmentation sources that heuristic
cleaning targets:

1. A dead-end **tip** — a read that follows the reference for a while and then
   diverges into a short erroneous suffix that goes nowhere. On the k-mer graph
   this adds an out-branch on the backbone (fragmenting the linear walk) plus a
   low-support dead-end tail.
2. A single-nucleotide **bubble** — a read identical to the reference except for
   one substituted base in the middle, creating an alternate low-support branch
   that rejoins the backbone (an SNP bubble).

Backbone (correct) reads are provided at higher multiplicity so the erroneous
tip/bubble branches are genuinely low-support relative to the truth.
"""
function build_toy_reads(; seed::Int = 42)
    Random.seed!(seed)
    alphabet = ['A', 'C', 'G', 'T']
    reference = join(rand(alphabet, 80))

    # Backbone: the full reference, at multiplicity 4 (high support truth path).
    reads = String[]
    for _ in 1:4
        push!(reads, reference)
    end

    # (1) Tip: reference prefix, then a short erroneous, dead-ending suffix.
    tip_read = reference[1:40] * "TTTTTTTT"
    push!(reads, tip_read)

    # (2) SNP bubble: reference with one interior base flipped (low support: 1).
    mid = 40
    orig = reference[mid]
    alt = orig == 'A' ? 'C' : 'A'
    snp_read = reference[1:(mid - 1)] * alt * reference[(mid + 1):end]
    push!(reads, snp_read)

    return reference, reads
end

"""
    iterated_remove_tips!(graph; min_support, max_passes=10)

OPT-IN helper. Repeatedly apply `remove_tips!` until the graph stops shrinking
(or `max_passes` is hit). A single `remove_tips!` pass only removes the current
terminal dead-end vertices; peeling a multi-k-mer tip back to the backbone
requires iterating. Purely a benchmark convenience — not an assembly default.
"""
function iterated_remove_tips!(graph; min_support::Int, max_passes::Int = 10)
    for _ in 1:max_passes
        before = Graphs.nv(graph.graph)
        Mycelia.Rhizomorph.remove_tips!(graph; min_support = min_support)
        after = Graphs.nv(graph.graph)
        after == before && break
    end
    return graph
end

"""
    assemble_toy_contigs(reference, reads; k, clean=false, min_support=1)

Build a k-mer assembly graph from `reads` and extract contigs.

`clean` is the OPT-IN heuristic-cleaning flag and DEFAULTS TO `false`. When
`false` (the default, and the only behaviour any real assembly path uses),
contigs are extracted from the raw graph. When `true`, an explicit heuristic
cleaning pass — iterated tip removal followed by bubble popping — runs BEFORE
contig extraction. Cleaning is never on by default anywhere.

Returns a NamedTuple of fragmentation metrics.
"""
function assemble_toy_contigs(
        reference::String,
        reads::Vector{String};
        k::Int = 11,
        clean::Bool = false,
        min_support::Int = 1
)
    # Stage reads to a temporary FASTA and build a single-strand k-mer graph.
    metrics = mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        records = [FASTX.FASTA.Record("read_$i", BioSequences.LongDNA{4}(r))
                   for (i, r) in enumerate(reads)]
        Mycelia.write_fasta(outfile = fasta, records = records)

        graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(
            fasta, k; mode = :singlestrand)

        n_vertices_raw = Graphs.nv(graph.graph)

        # --- OPT-IN heuristic cleaning (default OFF) -----------------------
        if clean
            # (a) Dead-end tip pruning: peel low-support dead-ends back to the
            #     high-support backbone.
            iterated_remove_tips!(graph; min_support = min_support)
            # (b) Bubble popping: detect SNP/indel bubbles and drop the
            #     clearly-dominated alternate path.
            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 2 * k)
            graph = Mycelia.Rhizomorph.simplify_graph_next(graph, bubbles)
        end
        # -------------------------------------------------------------------

        contigs = Mycelia.Rhizomorph.find_contigs_next(graph; min_contig_length = 1)

        contig_lengths = [c.length for c in contigs]
        total_len = isempty(contig_lengths) ? 0 : sum(contig_lengths)
        largest = isempty(contig_lengths) ? 0 : maximum(contig_lengths)

        (
            n_vertices = Graphs.nv(graph.graph),
            n_vertices_raw = n_vertices_raw,
            n_contigs = length(contigs),
            total_length = total_len,
            largest_contig = largest,
            genome_fraction_proxy = largest / length(reference)
        )
    end
    return metrics
end

function run_benchmark(; k::Int = 11, seed::Int = 42)
    reference, reads = build_toy_reads(seed = seed)

    println("=" ^ 72)
    println("Heuristic graph cleaning — COMPARISON BASELINE (not primary mechanism)")
    println("=" ^ 72)
    println("Reference length : $(length(reference)) bp")
    println("Reads            : $(length(reads)) (4x backbone + 1 tip + 1 SNP)")
    println("k                : $k")
    println()

    without = assemble_toy_contigs(reference, reads; k = k, clean = false)
    with = assemble_toy_contigs(reference, reads; k = k, clean = true, min_support = 1)

    function report(label, m)
        println(label)
        println("    graph vertices        : $(m.n_vertices)  (raw: $(m.n_vertices_raw))")
        println("    contigs               : $(m.n_contigs)")
        println("    total assembled length: $(m.total_length) bp")
        println("    largest contig        : $(m.largest_contig) bp")
        println("    genome-fraction proxy : $(round(m.genome_fraction_proxy, digits = 3))")
    end

    report("WITHOUT heuristic cleaning (raw graph — the assembly default):", without)
    println()
    report("WITH heuristic cleaning (OPT-IN tip removal + bubble pop):", with)
    println()

    delta_contigs = without.n_contigs - with.n_contigs
    println("-" ^ 72)
    println("Fragmentation removed by heuristic cleaning:")
    println("    Δ contigs             : $(delta_contigs) fewer " *
            "($(without.n_contigs) → $(with.n_contigs))")
    println("    Δ genome-fraction     : " *
            "$(round(without.genome_fraction_proxy, digits = 3)) → " *
            "$(round(with.genome_fraction_proxy, digits = 3))")
    println("-" ^ 72)
    println()
    println("NOTE: this is the BASELINE the emergent soft-EM cleaning is measured")
    println("against, per the assembly ADR — NOT the primary cleaning mechanism.")
    println("=" ^ 72)

    return (; reference, without, with, delta_contigs)
end

# Execute when run as a script (not when included for its function definitions).
if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark()
end
