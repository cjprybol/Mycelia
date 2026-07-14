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

    Test.@testset "overlong indel regions retain anchored tails" begin
        # Substitution-only calls preserve the historical first-window cap. The
        # explicit complete-span indel mode covers the tail with k-base overlaps,
        # giving every later pair-HMM decode a complete immutable start anchor.
        read_length = 1_010
        sequence = first(repeat("ACGT", cld(read_length, 4)), read_length)
        read = FASTX.FASTQ.Record(
            "overlong_hard_region",
            sequence,
            repeat("I", read_length),
        )
        graph = R.build_qualmer_graph(
            FASTX.FASTQ.Record[read], k; mode = :canonical)
        hard = Set(R.MetaGraphsNext.labels(graph))
        legacy_windows = Mycelia._hard_window_ranges(
            read,
            k,
            hard;
            pad = 0,
            max_window = 500,
            graph_mode = :canonical,
        )
        Test.@test legacy_windows == UnitRange{Int}[1:500]

        windows = Mycelia._hard_window_ranges(
            read,
            k,
            hard;
            pad = 0,
            max_window = 500,
            graph_mode = :canonical,
            complete_span = true,
        )

        Test.@test length(windows) == 3
        Test.@test first(first(windows)) == 1
        Test.@test last(last(windows)) == read_length
        Test.@test all(k <= length(window) <= 500 for window in windows)
        Test.@test all(
            last(windows[index]) - first(windows[index + 1]) + 1 == k
            for index in 1:(length(windows) - 1)
        )
        owned_bases = length(first(windows)) +
                      sum(length(window) - k for window in windows[2:end])
        Test.@test owned_bases == read_length

        # Exercise the production anchor-trim helper with true length changes,
        # then splice every balanced ownership range by original coordinates.
        sequence_chars = collect(sequence)
        quality_chars = collect(repeat("I", read_length))
        accepted = Tuple{UnitRange{Int}, String, String}[]

        first_window = first(windows)
        first_anchor_start = last(first_window) - k + 1
        first_prefix = sequence_chars[first(first_window):(first_anchor_start - 1)]
        first_anchor = sequence_chars[first_anchor_start:last(first_window)]
        # Insert immediately BEFORE the terminal anchor. A length change is valid,
        # but the decoded suffix must still name the graph state that anchors the
        # next overlapping window.
        first_decoded = String(vcat(first_prefix, ['A'], first_anchor))
        Test.@test Mycelia._indel_window_terminal_anchor_matches(
            sequence_chars,
            collect(first_decoded),
            last(first_window),
            k,
        )
        push!(accepted,
            (first_window, first_decoded, repeat("I", length(first_decoded))))

        second_window = windows[2]
        second_anchor_stop = first(second_window) + k - 1
        second_anchor = sequence_chars[first(second_window):second_anchor_stop]
        second_owned = sequence_chars[(first(second_window) + k):last(second_window)]
        # Delete the first owned base while retaining the immutable k-base anchor.
        second_decoded = vcat(second_anchor, second_owned[2:end])
        second_quality = fill('I', length(second_decoded))
        second_trimmed = Mycelia._trim_indel_window_overlap(
            sequence_chars,
            second_decoded,
            second_quality,
            first(second_window),
            k,
        )
        Test.@test second_trimmed !== nothing
        second_trimmed_sequence, second_trimmed_quality = something(second_trimmed)
        push!(accepted,
            (
                (first(second_window) + k):last(second_window),
                String(second_trimmed_sequence),
                String(second_trimmed_quality),
            ))

        third_window = windows[3]
        third_decoded = sequence_chars[third_window]
        third_quality = fill('I', length(third_decoded))
        third_trimmed = Mycelia._trim_indel_window_overlap(
            sequence_chars,
            third_decoded,
            third_quality,
            first(third_window),
            k,
        )
        Test.@test third_trimmed !== nothing
        third_trimmed_sequence, third_trimmed_quality = something(third_trimmed)
        push!(accepted,
            (
                (first(third_window) + k):last(third_window),
                String(third_trimmed_sequence),
                String(third_trimmed_quality),
            ))

        corrected = Mycelia._splice_indel_windows(
            "overlap_length_change",
            sequence_chars,
            quality_chars,
            accepted,
        )
        expected_sequence = first_decoded *
                            String(second_owned[2:end]) *
                            String(sequence_chars[(first(third_window) + k):end])
        Test.@test FASTX.sequence(String, corrected) == expected_sequence
        Test.@test length(FASTX.quality(corrected)) == length(expected_sequence)

        mutated_anchor = copy(third_decoded)
        mutated_anchor[1] = mutated_anchor[1] == 'A' ? 'C' : 'A'
        Test.@test Mycelia._trim_indel_window_overlap(
            sequence_chars,
            mutated_anchor,
            third_quality,
            first(third_window),
            k,
        ) === nothing

        mutated_terminal = collect(first_decoded)
        mutated_terminal[end] = mutated_terminal[end] == 'A' ? 'C' : 'A'
        Test.@test !Mycelia._indel_window_terminal_anchor_matches(
            sequence_chars,
            mutated_terminal,
            last(first_window),
            k,
        )
    end

    Test.@testset "doublestrand windows use observed graph orientation" begin
        branch_reads = FASTX.FASTQ.Record[
            FASTX.FASTQ.Record("r1", "TTTA", "IIII"),
            FASTX.FASTQ.Record("r2", "TTTC", "IIII"),
        ]
        graph = R.build_qualmer_graph(branch_reads, 3; mode = :doublestrand)
        hard = Mycelia._hard_vertex_set(graph, 3)
        Test.@test Set(string.(hard)) == Set(["TTT"])
        for read in branch_reads
            Test.@test Mycelia.should_decode_read(
                read, 3, hard; graph_mode = :doublestrand)
            Test.@test Mycelia._hard_window_ranges(
                read, 3, hard;
                pad = 3,
                max_window = 500,
                graph_mode = :doublestrand,
            ) == UnitRange{Int}[1:4]
        end
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
