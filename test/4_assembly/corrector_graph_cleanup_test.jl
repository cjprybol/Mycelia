# CORRECTOR GRAPH CLEANUP (td-969e)
# =================================
#
# Unit tests for the linear-time defragmentation pass that runs on the scalable
# corrector's final graph before contig extraction. The pass removes ONLY
# unambiguous errors — coverage-1 dead-end tips and guarded low-coverage bubble
# branches — and must NEVER collapse a data-supported variant or clip a
# stand-alone linear run / genome terminus (the td-h6w9 variation-preservation
# invariant, tested here at the graph-primitive level; the end-to-end invariant
# lives in variation_preservation_holdout_test.jl).
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_graph_cleanup_test.jl")'

import Test
import Mycelia
import FASTX
import Kmers
import MetaGraphsNext
import Random
import BioSequences

const _CGC = Mycelia.Rhizomorph
const _CGC_BASES = ['A', 'C', 'G', 'T']

# A random DNA string whose k-mers are (with overwhelming probability) all unique,
# so hand-built branch/bubble structures are clean.
_cgc_backbone(rng, L) = join(rand(rng, _CGC_BASES, L))

_cgc_fasta(id, seq) = FASTX.FASTA.Record(id, seq)

# Coverage (evidence count) of a k-mer label in a graph, or 0 if absent.
function _cgc_cov(graph, label)
    haskey(graph, label) || return 0
    return _CGC.count_evidence(graph[label])
end

Test.@testset "corrector graph cleanup primitives (td-969e)" begin
    k = 11

    Test.@testset "clip_error_tips!: coverage-1 dead-end tip removed, backbone kept" begin
        rng = Random.MersenneTwister(11)
        backbone = _cgc_backbone(rng, 60)
        # Error tip: shares the backbone prefix, then diverges into a novel tail
        # that dead-ends (out-degree 0), at coverage 1.
        errtail = _cgc_backbone(Random.MersenneTwister(999), 20)
        err_read = backbone[1:30] * errtail
        reads = FASTX.FASTA.Record[]
        for i in 1:10
            push!(reads, _cgc_fasta("main$i", backbone))   # backbone coverage 10
        end
        push!(reads, _cgc_fasta("err", err_read))          # tip coverage 1

        graph = _CGC.build_kmer_graph(reads, k; dataset_id = "t", mode = :singlestrand)

        # A k-mer fully inside the novel error tail is a coverage-1 dead-end tip.
        tip_kmer = Kmers.DNAKmer{k}(errtail[5:(5 + k - 1)])
        backbone_kmers = [Kmers.DNAKmer{k}(backbone[i:(i + k - 1)])
                          for i in 1:(length(backbone) - k + 1)]

        Test.@test _cgc_cov(graph, tip_kmer) == 1                 # present, coverage-1, pre-clean
        Test.@test all(kmref -> haskey(graph, kmref), backbone_kmers)

        res = _CGC.clip_error_tips!(graph; max_tip_length = 3 * k, max_tip_support = 1)

        Test.@test res.removed >= 1
        Test.@test !haskey(graph, tip_kmer)                       # error tip clipped
        # Every backbone k-mer (coverage 10) survives — no real sequence lost.
        Test.@test all(kmref -> haskey(graph, kmref), backbone_kmers)
    end

    Test.@testset "clip_error_tips!: stand-alone linear run (terminus) NOT clipped" begin
        rng = Random.MersenneTwister(22)
        backbone = _cgc_backbone(rng, 60)
        # A single linear genome with NO branches. Both ends are dead ends but each
        # reaches the other terminus without a junction, so nothing is a clippable
        # tip — even at coverage 1 (mimics a low-coverage genome end).
        graph = _CGC.build_kmer_graph([_cgc_fasta("g", backbone)], k;
            dataset_id = "t", mode = :singlestrand)
        n_before = length(collect(MetaGraphsNext.labels(graph)))

        res = _CGC.clip_error_tips!(graph; max_tip_length = 3 * k, max_tip_support = 1)

        Test.@test res.removed == 0
        Test.@test length(collect(MetaGraphsNext.labels(graph))) == n_before
    end

    Test.@testset "collapse_error_bubbles!: balanced bubble RETAINS both branches" begin
        rng = Random.MersenneTwister(33)
        backbone = collect(_cgc_backbone(rng, 60))
        pvar = 30
        base_a = backbone[pvar]
        base_b = rand(rng, filter(!=(base_a), _CGC_BASES))
        hap_a = copy(backbone)
        hap_b = copy(backbone); hap_b[pvar] = base_b
        reads = FASTX.FASTA.Record[]
        for i in 1:15
            push!(reads, _cgc_fasta("A$i", String(hap_a)))   # allele A 15x
            push!(reads, _cgc_fasta("B$i", String(hap_b)))   # allele B 15x
        end
        graph = _CGC.build_kmer_graph(reads, k; dataset_id = "t", mode = :singlestrand)

        # Allele-distinguishing k-mers (span the variant position).
        kmer_a = Kmers.DNAKmer{k}(String(hap_a[(pvar - k + 1):pvar]))
        kmer_b = Kmers.DNAKmer{k}(String(hap_b[(pvar - k + 1):pvar]))
        Test.@test haskey(graph, kmer_a)
        Test.@test haskey(graph, kmer_b)

        res = _CGC.collapse_error_bubbles!(graph; max_error_support = 2, min_real_support = 3)

        # NEITHER balanced branch is an error -> both retained (variation preserved).
        Test.@test res.collapsed == 0
        Test.@test haskey(graph, kmer_a)
        Test.@test haskey(graph, kmer_b)
    end

    Test.@testset "collapse_error_bubbles!: coverage-1 error bubble collapsed" begin
        rng = Random.MersenneTwister(44)
        backbone = collect(_cgc_backbone(rng, 60))
        pvar = 30
        base_true = backbone[pvar]
        base_err = rand(rng, filter(!=(base_true), _CGC_BASES))
        err_hap = copy(backbone); err_hap[pvar] = base_err
        reads = FASTX.FASTA.Record[]
        for i in 1:30
            push!(reads, _cgc_fasta("C$i", String(backbone)))  # consensus 30x
        end
        push!(reads, _cgc_fasta("err", String(err_hap)))       # error branch 1x

        graph = _CGC.build_kmer_graph(reads, k; dataset_id = "t", mode = :singlestrand)
        kmer_true = Kmers.DNAKmer{k}(String(backbone[(pvar - k + 1):pvar]))
        kmer_err = Kmers.DNAKmer{k}(String(err_hap[(pvar - k + 1):pvar]))
        Test.@test haskey(graph, kmer_true)
        Test.@test haskey(graph, kmer_err)
        Test.@test _cgc_cov(graph, kmer_err) == 1

        res = _CGC.collapse_error_bubbles!(graph; max_error_support = 2, min_real_support = 3)

        Test.@test res.collapsed >= 1
        Test.@test haskey(graph, kmer_true)      # consensus retained
        Test.@test !haskey(graph, kmer_err)      # coverage-1 error branch removed
    end

    Test.@testset "clean_corrector_graph!: composes + reports telemetry" begin
        rng = Random.MersenneTwister(55)
        backbone = collect(_cgc_backbone(rng, 80))
        # Real balanced variant at pvar (both 12x) + coverage-1 error tip.
        pvar = 40
        base_a = backbone[pvar]
        base_b = rand(rng, filter(!=(base_a), _CGC_BASES))
        hap_a = copy(backbone)
        hap_b = copy(backbone); hap_b[pvar] = base_b
        errtail = _cgc_backbone(Random.MersenneTwister(7777), 20)
        reads = FASTX.FASTA.Record[]
        for i in 1:12
            push!(reads, _cgc_fasta("A$i", String(hap_a)))
            push!(reads, _cgc_fasta("B$i", String(hap_b)))
        end
        push!(reads, _cgc_fasta("errtip", String(hap_a[1:20]) * errtail))

        graph = _CGC.build_kmer_graph(reads, k; dataset_id = "t", mode = :singlestrand)
        kmer_a = Kmers.DNAKmer{k}(String(hap_a[(pvar - k + 1):pvar]))
        kmer_b = Kmers.DNAKmer{k}(String(hap_b[(pvar - k + 1):pvar]))

        stats = _CGC.clean_corrector_graph!(graph; k = k)

        Test.@test stats["graph_cleanup_vertices_after"] <= stats["graph_cleanup_vertices_before"]
        Test.@test stats["graph_cleanup_tips_removed"] >= 1
        # Real balanced variant SURVIVES the composed cleanup (both alleles kept).
        Test.@test haskey(graph, kmer_a)
        Test.@test haskey(graph, kmer_b)
    end

    # Reverse complement of a plain DNA string (for the real-sequence guard test).
    _cgc_rc(s) = string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(s)))

    Test.@testset "prune_disconnected_error_components!: coverage-1 island pruned, backbone kept" begin
        rng = Random.MersenneTwister(66)
        backbone = _cgc_backbone(rng, 60)                 # main genome, coverage 10
        # A totally unrelated coverage-1 read: a SEPARATE connected component of
        # error k-mers that shares no k-mer with the backbone.
        island = _cgc_backbone(Random.MersenneTwister(4242), 30)
        reads = FASTX.FASTA.Record[]
        for i in 1:10
            push!(reads, _cgc_fasta("main$i", backbone))
        end
        push!(reads, _cgc_fasta("island", island))

        graph = _CGC.build_kmer_graph(reads, k; dataset_id = "t", mode = :singlestrand)
        island_kmer = Kmers.DNAKmer{k}(island[10:(10 + k - 1)])
        backbone_kmers = [Kmers.DNAKmer{k}(backbone[i:(i + k - 1)])
                          for i in 1:(length(backbone) - k + 1)]
        Test.@test _cgc_cov(graph, island_kmer) == 1
        Test.@test all(kmref -> haskey(graph, kmref), backbone_kmers)

        res = _CGC.prune_disconnected_error_components!(graph; k = k,
            max_component_support = 2, min_real_support = 3)

        Test.@test res.components_pruned >= 1
        Test.@test res.removed >= 1
        Test.@test !haskey(graph, island_kmer)                        # error island removed
        Test.@test all(kmref -> haskey(graph, kmref), backbone_kmers) # backbone untouched
    end

    Test.@testset "prune_disconnected_error_components!: high-coverage disconnected genome RETAINED" begin
        rng = Random.MersenneTwister(77)
        backbone = _cgc_backbone(rng, 80)                 # larger main component
        # A SECOND real genome, disconnected but well-supported (coverage 10). The
        # coverage gate must retain it even though it is not the main component.
        second = _cgc_backbone(Random.MersenneTwister(8484), 50)
        reads = FASTX.FASTA.Record[]
        for i in 1:10
            push!(reads, _cgc_fasta("main$i", backbone))
            push!(reads, _cgc_fasta("second$i", second))
        end

        graph = _CGC.build_kmer_graph(reads, k; dataset_id = "t", mode = :singlestrand)
        second_kmer = Kmers.DNAKmer{k}(second[20:(20 + k - 1)])
        Test.@test _cgc_cov(graph, second_kmer) == 10

        res = _CGC.prune_disconnected_error_components!(graph; k = k,
            max_component_support = 2, min_real_support = 3)

        Test.@test res.components_pruned == 0
        Test.@test haskey(graph, second_kmer)             # well-supported genome kept
    end

    Test.@testset "prune_disconnected_error_components!: real-sequence guard RETAINS overlapping island" begin
        rng = Random.MersenneTwister(88)
        backbone = _cgc_backbone(rng, 60)                 # main genome, coverage 10
        # A coverage-1 read that is the REVERSE COMPLEMENT of a backbone segment.
        # In :singlestrand mode its k-mers are distinct vertices (a separate
        # component) but their CANONICAL forms match the backbone's, so the
        # real-sequence guard must RETAIN it (it could be genuine reverse-strand
        # minor sequence), even though it is small and coverage-1.
        rc_segment = _cgc_rc(backbone[1:30])
        reads = FASTX.FASTA.Record[]
        for i in 1:10
            push!(reads, _cgc_fasta("main$i", backbone))
        end
        push!(reads, _cgc_fasta("rc", rc_segment))

        graph = _CGC.build_kmer_graph(reads, k; dataset_id = "t", mode = :singlestrand)
        rc_kmer = Kmers.DNAKmer{k}(rc_segment[10:(10 + k - 1)])
        Test.@test haskey(graph, rc_kmer)
        Test.@test _cgc_cov(graph, rc_kmer) == 1          # low coverage, but real content

        res = _CGC.prune_disconnected_error_components!(graph; k = k,
            max_component_support = 2, min_real_support = 3)

        # Shares canonical k-mers with the well-supported backbone -> retained.
        Test.@test res.components_pruned == 0
        Test.@test haskey(graph, rc_kmer)
    end

    Test.@testset "prune_disconnected_error_components!: single component -> no-op" begin
        rng = Random.MersenneTwister(99)
        backbone = _cgc_backbone(rng, 60)
        graph = _CGC.build_kmer_graph([_cgc_fasta("g", backbone)], k;
            dataset_id = "t", mode = :singlestrand)
        n_before = length(collect(MetaGraphsNext.labels(graph)))

        res = _CGC.prune_disconnected_error_components!(graph; k = k)

        Test.@test res.removed == 0
        Test.@test res.components_pruned == 0
        Test.@test length(collect(MetaGraphsNext.labels(graph))) == n_before
    end

    Test.@testset "clean_corrector_graph!: prunes error island + reports telemetry" begin
        rng = Random.MersenneTwister(101)
        backbone = _cgc_backbone(rng, 80)
        island = _cgc_backbone(Random.MersenneTwister(1357), 30)
        reads = FASTX.FASTA.Record[]
        for i in 1:10
            push!(reads, _cgc_fasta("main$i", backbone))
        end
        push!(reads, _cgc_fasta("island", island))

        graph = _CGC.build_kmer_graph(reads, k; dataset_id = "t", mode = :singlestrand)
        island_kmer = Kmers.DNAKmer{k}(island[10:(10 + k - 1)])
        backbone_kmers = [Kmers.DNAKmer{k}(backbone[i:(i + k - 1)])
                          for i in 1:(length(backbone) - k + 1)]

        stats = _CGC.clean_corrector_graph!(graph; k = k)

        Test.@test stats["graph_cleanup_components_pruned"] >= 1
        Test.@test stats["graph_cleanup_component_vertices_removed"] >= 1
        Test.@test !haskey(graph, island_kmer)                        # island gone
        Test.@test all(kmref -> haskey(graph, kmref), backbone_kmers) # backbone intact
    end

    # ------------------------------------------------------------------
    # Size-cap opt-out (td-byva): a LARGE (span > default 2000) uniformly
    # coverage-1 disconnected island sharing no real k-mer content is genuine
    # error debris, but the default size cap RETAINS it (it might be an under-
    # sequenced real replicon). `max_component_length = nothing` disables that
    # conservatism knob so the island is pruned — while the coverage floor and
    # real-sequence guard stay binding (proved by the two follow-on cases).
    # ------------------------------------------------------------------
    # A larger k for the large-island cases so a ~2050 bp random sequence has
    # (w.h.p.) all-unique k-mers -> uniform coverage 1, no internal repeat that
    # would spuriously raise a vertex above the coverage floor. span at k=kk is
    # (n_vertices + kk - 1) ~= 2050 > the default 2000 cap.
    kk = 21
    Test.@testset "prune: large coverage-1 island retained by default cap, pruned when disabled" begin
        rng = Random.MersenneTwister(202)
        # Main genome must be LARGER than the island so argmax(component length)
        # keeps the island a prune CANDIDATE (not the retained main component).
        backbone = _cgc_backbone(rng, 3000)               # main genome, coverage 10
        # ~2050 bp error island => span > 2000, coverage 1, shares NO k-mer with
        # the backbone.
        big_island = _cgc_backbone(Random.MersenneTwister(30303), 2050)
        base_reads = FASTX.FASTA.Record[]
        for i in 1:10
            push!(base_reads, _cgc_fasta("main$i", backbone))
        end
        push!(base_reads, _cgc_fasta("big_island", big_island))
        island_kmer = Kmers.DNAKmer{kk}(big_island[100:(100 + kk - 1)])
        backbone_kmers = [Kmers.DNAKmer{kk}(backbone[i:(i + kk - 1)])
                          for i in 1:(length(backbone) - kk + 1)]

        # DEFAULT cap: the large island is above the span cap -> RETAINED.
        graph_default = _CGC.build_kmer_graph(base_reads, kk;
            dataset_id = "t", mode = :singlestrand)
        Test.@test _cgc_cov(graph_default, island_kmer) == 1
        res_default = _CGC.prune_disconnected_error_components!(graph_default; k = kk,
            max_component_support = 2, min_real_support = 3)
        Test.@test res_default.components_pruned == 0
        Test.@test haskey(graph_default, island_kmer)                 # retained by size cap

        # DISABLED cap: coverage floor + guard both pass for this error island,
        # so with the span gate off it is pruned; the backbone is untouched.
        graph_nocap = _CGC.build_kmer_graph(base_reads, kk;
            dataset_id = "t", mode = :singlestrand)
        res_nocap = _CGC.prune_disconnected_error_components!(graph_nocap; k = kk,
            max_component_length = nothing,
            max_component_support = 2, min_real_support = 3)
        Test.@test res_nocap.components_pruned >= 1
        Test.@test res_nocap.removed >= 1
        Test.@test !haskey(graph_nocap, island_kmer)                  # error island removed
        Test.@test all(kmref -> haskey(graph_nocap, kmref), backbone_kmers)
    end

    Test.@testset "prune: cap disabled STILL retains large high-coverage disconnected genome" begin
        rng = Random.MersenneTwister(303)
        # Main genome larger than the second so the second stays a prune CANDIDATE
        # and its retention is genuinely due to the coverage floor.
        backbone = _cgc_backbone(rng, 3000)               # main component
        # A large SECOND real genome (~2050 bp), disconnected but well-supported
        # (coverage 10). The coverage floor must retain it even with the span gate
        # off -> disabling the cap does NOT weaken the coverage-floor safety proof.
        second = _cgc_backbone(Random.MersenneTwister(40404), 2050)
        reads = FASTX.FASTA.Record[]
        for i in 1:10
            push!(reads, _cgc_fasta("main$i", backbone))
            push!(reads, _cgc_fasta("second$i", second))
        end
        graph = _CGC.build_kmer_graph(reads, kk; dataset_id = "t", mode = :singlestrand)
        second_kmer = Kmers.DNAKmer{kk}(second[100:(100 + kk - 1)])
        Test.@test _cgc_cov(graph, second_kmer) == 10

        res = _CGC.prune_disconnected_error_components!(graph; k = kk,
            max_component_length = nothing,
            max_component_support = 2, min_real_support = 3)

        Test.@test res.components_pruned == 0
        Test.@test haskey(graph, second_kmer)             # well-supported genome kept
    end

    Test.@testset "prune: cap disabled STILL retains large real-content-sharing island (guard)" begin
        rng = Random.MersenneTwister(404)
        backbone = _cgc_backbone(rng, 2100)               # large main genome, coverage 10
        # A large coverage-1 read that is the REVERSE COMPLEMENT of a backbone
        # segment: distinct singlestrand vertices (separate component) whose
        # CANONICAL k-mers match the backbone. The real-sequence guard must retain
        # it even with the span gate off -> disabling the cap does NOT weaken the
        # real-sequence guard.
        rc_segment = _cgc_rc(backbone[1:2050])
        reads = FASTX.FASTA.Record[]
        for i in 1:10
            push!(reads, _cgc_fasta("main$i", backbone))
        end
        push!(reads, _cgc_fasta("rc", rc_segment))
        graph = _CGC.build_kmer_graph(reads, kk; dataset_id = "t", mode = :singlestrand)
        rc_kmer = Kmers.DNAKmer{kk}(rc_segment[100:(100 + kk - 1)])
        Test.@test haskey(graph, rc_kmer)
        Test.@test _cgc_cov(graph, rc_kmer) == 1          # low coverage, but real content

        res = _CGC.prune_disconnected_error_components!(graph; k = kk,
            max_component_length = nothing,
            max_component_support = 2, min_real_support = 3)

        Test.@test res.components_pruned == 0             # guard binding despite no size cap
        Test.@test haskey(graph, rc_kmer)
    end
end
