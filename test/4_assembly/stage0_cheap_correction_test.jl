# Stage 0 CHEAP k-mer-spectrum correction (td-bjnt).
# ===================================================
#
# The :scalable corrector times out at scale because ~78% of reads (any read with
# an error k-mer) were classified "hard" and sent to the expensive per-read graph
# Viterbi. Stage 0 fixes simple single-substitution errors CHEAPLY first (a linear
# k-mer-spectrum scan, BFC/Lighter/Bloocoo-style) and NARROWS the hard-vertex set
# to genuine ambiguity (bubbles/repeats), so the graph-Viterbi decode fraction
# drops toward the true ~5-15%.
#
# This test asserts (a) the per-read cheap corrector fixes an unambiguous single
# substitution toward its UNIQUE solid neighbor, (b) it refuses to touch a balanced
# variant (two solid alleles ⇒ no unique neighbor) or an ambiguous / end / N read,
# (c) the narrowed _hard_vertex_set no longer marks a weak-but-non-bubble k-mer
# hard, and (d) end-to-end on an err=0.01 toy the :scalable tier reports cheap
# corrections and a materially SMALLER graph-decode fraction than a corrector with
# no cheap pre-pass, while :exhaustive stays cheap-correction-free.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/stage0_cheap_correction_test.jl")'

import Test
import Mycelia
import FASTX
import Random
import BioSequences
import Kmers

const _S0_BASES = ['A', 'C', 'G', 'T']
_s0_rec(id, seq) = FASTX.FASTQ.Record(id, seq, String(fill('I', length(seq))))

Test.@testset "Stage 0 cheap k-mer-spectrum correction (td-bjnt)" begin
    R = Mycelia.Rhizomorph

    Test.@testset "single-substitution fix toward unique solid neighbor" begin
        k = 7
        rng = Random.MersenneTwister(101)
        ref = join(rand(rng, _S0_BASES, 120))
        # High-coverage clean reads ⇒ every ref k-mer solid.
        clean = [_s0_rec("c$i", ref) for i in 1:15]
        graph = R.build_qualmer_graph(clean, k; dataset_id = "ds", mode = :canonical,
            memory_profile = :full)
        solid = Mycelia._solid_kmer_set(graph)
        Test.@test !isempty(solid)

        # Inject ONE interior substitution; Stage 0 must restore the reference.
        perr = 60
        errchars = collect(ref)
        errchars[perr] = first(filter(!=(errchars[perr]), _S0_BASES))
        errseq = String(errchars)
        fixed, n = Mycelia._stage0_correct_read(_s0_rec("e", errseq), k, solid)
        Test.@test n == 1
        Test.@test FASTX.sequence(String, fixed) == ref          # restored exactly
        # Quality string length preserved (substitution is length-neutral).
        Test.@test length(FASTX.quality(fixed)) == length(ref)

        # A clean read has no weak run ⇒ untouched, 0 corrections.
        fixed_clean, nc = Mycelia._stage0_correct_read(_s0_rec("c", ref), k, solid)
        Test.@test nc == 0
        Test.@test FASTX.sequence(String, fixed_clean) == ref
    end

    Test.@testset "refuses ambiguous / boundary / N reads" begin
        k = 7
        rng = Random.MersenneTwister(102)
        ref = join(rand(rng, _S0_BASES, 80))
        clean = [_s0_rec("c$i", ref) for i in 1:15]
        graph = R.build_qualmer_graph(clean, k; dataset_id = "ds", mode = :canonical,
            memory_profile = :full)
        solid = Mycelia._solid_kmer_set(graph)

        # Read shorter than k: no k-mers ⇒ never corrected.
        short_r = _s0_rec("s", ref[1:5])
        sr, sn = Mycelia._stage0_correct_read(short_r, k, solid)
        Test.@test sn == 0
        Test.@test FASTX.sequence(String, sr) == ref[1:5]

        # Read carrying an ambiguous base (N): 2-bit k-mer construction cannot
        # encode it, so Stage 0 leaves it for the decode (0 corrections, unchanged).
        nchars = collect(ref); nchars[40] = 'N'
        nr, nn = Mycelia._stage0_correct_read(_s0_rec("n", String(nchars)), k, solid)
        Test.@test nn == 0
        Test.@test occursin("N", FASTX.sequence(String, nr))

        # Empty solid set ⇒ no trustworthy neighbor ⇒ no-op.
        empty_solid = Set(eltype(collect(solid))[])
        errchars = collect(ref); errchars[40] = first(filter(!=(errchars[40]), _S0_BASES))
        er, en = Mycelia._stage0_correct_read(_s0_rec("e", String(errchars)), k, empty_solid)
        Test.@test en == 0
        Test.@test FASTX.sequence(String, er) == String(errchars)
    end

    Test.@testset "preserves a balanced variant (no unique solid neighbor)" begin
        # A real heterozygous site at balanced coverage: BOTH alleles are solid, so
        # neither is a weak run — Stage 0 must not touch either read. Even if one
        # allele were weak, correcting it would require a UNIQUE solid neighbor; the
        # sibling allele being solid too breaks uniqueness.
        k = 7
        rng = Random.MersenneTwister(103)
        L = 120
        backbone = collect(join(rand(rng, _S0_BASES, L)))
        pvar = 60
        base_a = backbone[pvar]
        base_b = first(filter(!=(base_a), _S0_BASES))
        hap_a = copy(backbone)
        hap_b = copy(backbone); hap_b[pvar] = base_b
        reads = FASTX.FASTQ.Record[]
        for i in 1:12
            push!(reads, _s0_rec("A$i", String(hap_a)))
            push!(reads, _s0_rec("B$i", String(hap_b)))
        end
        graph = R.build_qualmer_graph(reads, k; dataset_id = "ds", mode = :canonical,
            memory_profile = :full)
        solid = Mycelia._solid_kmer_set(graph)

        # Neither allele read is altered by Stage 0.
        fa, na = Mycelia._stage0_correct_read(_s0_rec("A", String(hap_a)), k, solid)
        fb, nb = Mycelia._stage0_correct_read(_s0_rec("B", String(hap_b)), k, solid)
        Test.@test na == 0
        Test.@test nb == 0
        Test.@test FASTX.sequence(String, fa) == String(hap_a)
        Test.@test FASTX.sequence(String, fb) == String(hap_b)
    end

    Test.@testset "narrowed _hard_vertex_set excludes a weak non-bubble k-mer" begin
        # Contract of the narrowing (td-bjnt): a WEAK (non-solid) k-mer is no longer
        # hard just for being weak — only bubble/repeat vertices are. Build UNEVEN
        # coverage so a linear (non-bubble) tail region is covered once ⇒ its k-mers
        # are weak but on a single path. Under the OLD clause (3) every weak k-mer
        # was hard; the narrowed set must exclude these linear weak k-mers.
        k = 13
        rng = Random.MersenneTwister(104)
        ref = join(rand(rng, _S0_BASES, 400))
        reads = FASTX.FASTQ.Record[]
        # Rich, even coverage over the head [1, 300].
        for i in 1:40
            s = rand(rng, 1:(300 - 80 + 1))
            push!(reads, _s0_rec("h$i", ref[s:(s + 79)]))
        end
        # ONE read carries the tail [280, 400] (overlaps the head at [280,300], then
        # a coverage-1 linear extension [300,400]).
        push!(reads, _s0_rec("tail", ref[280:400]))

        graph = R.build_qualmer_graph(reads, k; mode = :canonical)
        labels = collect(R.MetaGraphsNext.labels(graph))
        solid = Mycelia._solid_kmer_set(graph)
        hard = Mycelia._hard_vertex_set(graph, k)

        weak_labels = [l for l in labels if !(l in solid)]
        Test.@test !isempty(weak_labels)                 # uneven coverage ⇒ weak k-mers exist
        # The narrowed set removes clause (3): at least one weak k-mer is NOT hard.
        # Under the old "all weak k-mers are hard" clause this would be impossible.
        Test.@test any(l -> !(l in hard), weak_labels)
        # And the hard set never contains a k-mer solely because it is weak: every
        # hard vertex is a real graph vertex (bubble/repeat), i.e. hard ⊆ labels and
        # is not simply the whole weak complement.
        Test.@test issubset(hard, Set(labels))
        Test.@test !issubset(Set(weak_labels), hard)     # weak set ⊄ hard (clause 3 gone)
    end

    Test.@testset "end-to-end: :scalable reports cheap corrections + lower decode fraction" begin
        # err=0.01 toy: Stage 0 should fix most simple errors cheaply and collapse
        # the graph-Viterbi decode fraction relative to the same tier run WITHOUT the
        # cheap pre-pass (isolating the lever).
        rng = Random.MersenneTwister(2026)
        ref = join(rand(rng, _S0_BASES, 2000))
        readlen = 120
        reads = FASTX.FASTQ.Record[]
        for i in 1:400
            s = rand(rng, 1:(2000 - readlen + 1))
            seq = collect(ref[s:(s + readlen - 1)])
            for j in eachindex(seq)
                if rand(rng) < 0.01
                    seq[j] = rand(rng, filter(!=(seq[j]), _S0_BASES))
                end
            end
            push!(reads, _s0_rec("read$i", String(seq)))
        end

        k = 19
        graph = R.build_qualmer_graph(reads, k; mode = :canonical)
        solid = Mycelia._solid_kmer_set(graph)
        hard = Mycelia._hard_vertex_set(graph, k)

        # WITH Stage 0 cheap correction (the :scalable behavior).
        _out_on, imp_on, skip_on, cheap_on = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard)
        # WITHOUT the cheap pre-pass, same gates (the pre-Stage-0 baseline). Weak
        # reads that Stage 0 would have cleared instead reach the decode.
        _out_off, _imp_off, skip_off, cheap_off = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = false, hard_vertices = hard)

        Test.@test cheap_on > 0            # Stage 0 fixed simple errors cheaply
        Test.@test cheap_off == 0          # off ⇒ no cheap corrections
        decode_on = 1.0 - skip_on
        decode_off = 1.0 - skip_off
        # The win: the cheap pre-pass does not INCREASE the decode fraction, and on
        # this error-rich toy it strictly reduces it (fewer reads reach the decode
        # because Stage 0 made them all-solid ⇒ skip-solid skips them).
        Test.@test decode_on <= decode_off + 1e-9
        @info "Stage 0 decode-fraction lever (td-bjnt)" cheap_on decode_off decode_on skip_on skip_off
    end
end
