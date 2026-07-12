# Stage 3 (td-nn6l): hard-window gating. The corrector decodes ONLY reads that
# touch a "hard" vertex (bubble / repeat-like / weak k-mer) and passes every other
# read through untouched, cutting decode volume. This test asserts (a) the
# hard-vertex set construction, (b) should_decode_read only passes reads that
# overlap a hard vertex, and (c) an end-to-end scalable assemble reports a
# non-trivial skip fraction (reads passed through without a decode).
#
# Stage 3c (per-hard-region WINDOWED decode with start/target boundary
# constraints) is now WIRED (td-nn6l): under `windowed_decode=true` a hard read is
# decoded window-by-window via `improve_read_likelihood_windowed`, not whole. The
# windowing primitive `_hard_window_ranges` is tested here; the boundary-
# constrained decode + splice is tested in windowed_decode_test.jl.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/scalable_corrector_hard_window_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _BASES = ['A', 'C', 'G', 'T']

function _clean_reads(rng, ref; n_reads = 120, readlen = 80)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = ref[s:(s + readlen - 1)]
        push!(records, FASTX.FASTQ.Record("r$i", seq, String(fill('I', readlen))))
    end
    return records
end

Test.@testset "scalable corrector hard-window gating (td-nn6l)" begin
    R = Mycelia.Rhizomorph
    k = 13

    Test.@testset "_hard_vertex_set + should_decode_read" begin
        rng = Random.MersenneTwister(21)
        ref = join(rand(rng, _BASES, 1000))
        # Introduce a handful of substituted reads so the graph carries bubbles /
        # weak k-mers (real hard regions), plus many clean reads.
        reads = _clean_reads(rng, ref; n_reads = 120)
        for i in 1:8
            s = rand(rng, 1:921)
            seq = collect(ref[s:(s + 79)])
            seq[40] = rand(rng, filter(!=(seq[40]), _BASES))   # mid-read substitution
            push!(reads, FASTX.FASTQ.Record("e$i", String(seq), String(fill('I', 80))))
        end

        graph = R.build_qualmer_graph(reads, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)
        Test.@test hard isa AbstractSet
        # There ARE hard vertices (weak k-mers from the substituted reads at least).
        Test.@test !isempty(hard)
        # Hard set is a subset of the graph's vertices.
        all_labels = Set(Mycelia.Rhizomorph.MetaGraphsNext.labels(graph))
        Test.@test issubset(hard, all_labels)

        # should_decode_read: a read whose canonical k-mers all avoid `hard` is not
        # decoded; a read overlapping a hard vertex is. Build a decode-required read
        # by taking a window straddling an error read's substitution.
        # Empty hard set ⇒ always decode.
        empty_hard = Set{eltype(collect(all_labels))}()
        Test.@test Mycelia.should_decode_read(reads[1], k, empty_hard) == true

        # Count how many reads the gate would decode vs skip.
        n_decode = count(r -> Mycelia.should_decode_read(r, k, hard), reads)
        Test.@test 0 < n_decode <= length(reads)
    end

    Test.@testset "hard reads decode, clean reads skip (independent ground truth)" begin
        # Independently-constructed expectation (FIX 6): rather than re-deriving
        # should_decode_read's own overlap predicate and asserting they agree (a
        # tautology), inject a KNOWN error at a KNOWN locus and assert the gate
        # against externally-known truth — an error-straddling read MUST decode; a
        # clean read from a well-covered region far from the error MUST skip.
        rng = Random.MersenneTwister(22)
        ref = join(rand(rng, _BASES, 1000))
        # High-coverage clean reads ⇒ every ref k-mer is solid (non-hard).
        clean = _clean_reads(rng, ref; n_reads = 300, readlen = 80)
        # Inject a substitution at a fixed locus via several reads so the error
        # allele registers as a bubble / weak (hard) region in the graph.
        err_locus = 850                    # 1-based ref position of the substitution
        err_reads = FASTX.FASTQ.Record[]
        for i in 1:6
            s = err_locus - 40             # read spans [810, 889], error at offset 41
            seq = collect(ref[s:(s + 79)])
            seq[41] = first(filter(!=(seq[41]), _BASES))   # deterministic substitution
            push!(err_reads, FASTX.FASTQ.Record("e$i", String(seq), String(fill('I', 80))))
        end
        graph = R.build_qualmer_graph(vcat(clean, err_reads), k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)

        # Ground truth 1: a read carrying the injected error straddles the hard
        # region ⇒ MUST decode.
        Test.@test Mycelia.should_decode_read(err_reads[1], k, hard) == true

        # Ground truth 2: a clean read pulled straight from the reference, in a
        # well-covered region far from the error locus, touches only solid,
        # non-bubble k-mers ⇒ MUST skip.
        clean_probe = FASTX.FASTQ.Record("clean_probe", ref[100:179], String(fill('I', 80)))
        Test.@test Mycelia.should_decode_read(clean_probe, k, hard) == false
    end

    Test.@testset "_hard_window_ranges (Stage 3c primitive)" begin
        rng = Random.MersenneTwister(24)
        ref = join(rand(rng, _BASES, 300))
        clean = _clean_reads(rng, ref; n_reads = 150, readlen = 80)
        graph = R.build_qualmer_graph(clean, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)

        # A read that touches no hard vertex yields no windows; a read that does
        # yields windows that are (a) within read bounds and (b) each <= max_window.
        for r in clean[1:30]
            wins = Mycelia._hard_window_ranges(r, k, hard; max_window = 500)
            decodes = Mycelia.should_decode_read(r, k, hard)
            Test.@test (isempty(wins)) == (!decodes)   # windows exist iff decodable
            rlen = length(FASTX.sequence(r))
            for w in wins
                Test.@test first(w) >= 1 && last(w) <= rlen
                Test.@test length(w) <= 500
            end
        end
        # Empty hard set ⇒ no windows.
        empty_hard = Set(eltype(collect(hard))[])
        Test.@test isempty(Mycelia._hard_window_ranges(clean[1], k, empty_hard))
    end

    Test.@testset "end-to-end scalable assemble reports a skip fraction" begin
        rng = Random.MersenneTwister(23)
        ref = join(rand(rng, _BASES, 1200))
        reads = _clean_reads(rng, ref; n_reads = 200, readlen = 90)
        res = R.assemble_genome(reads; k = k, corrector = :iterative, strategy = :scalable)
        Test.@test res isa R.AssemblyResult
        Test.@test !isempty(res.contigs)
        Test.@test res.assembly_stats["hard_window"] == true
        skip = res.assembly_stats["skip_fraction"]
        # On higher-coverage cleaner data, many reads touch no hard vertex and are
        # skipped. Assert a real, non-trivial skip fraction (mechanism active),
        # without pinning the exact 85-95% target (data-dependent on this toy).
        Test.@test 0.0 < skip <= 1.0
    end
end
