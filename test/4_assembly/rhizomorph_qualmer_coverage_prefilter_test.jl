# Qualmer coverage-prefilter (min_count) — td-ck03.
#
# The qualmer graph stored a full k-length quality vector per k-mer OCCURRENCE
# (O(reads x read_len)), with NO coverage prefilter, so it did not scale to
# bacterial genomes (E. coli @ 20-50x hit 75-200GB+ and never completed). The
# fix adds an opt-in `min_count` (default 1 = exact no-op) that drops k-mers
# below the floor before creating vertices/evidence, bounding memory to
# O(distinct SOLID k-mers). These testsets fail until E1-E6 land.
#
# Run: julia --project=. -e 'include("test/4_assembly/rhizomorph_qualmer_coverage_prefilter_test.jl")'
import Test
import Mycelia
import FASTX
import BioSequences
import MetaGraphsNext

const PF_REF = "ATCGGCTAATGCCGATTGCACGTACGTTAGCTAGGCATG"  # no exact repeats at k=7

# One overlapping FASTQ tiling of `ref` (all Q40).
function pf_tiling(ref::AbstractString; read_len::Int = 15)
    rs = FASTX.FASTQ.Record[]
    for i in 1:(length(ref) - read_len + 1)
        sub = ref[i:(i + read_len - 1)]
        push!(rs, FASTX.FASTQ.Record("r$(i)", sub, repeat("I", length(sub))))
    end
    return rs
end

# `depth` distinct-observation copies of the tiling (re-IDed so evidence keys differ).
function pf_deep(ref::AbstractString; read_len::Int = 15, depth::Int = 40)
    rs = FASTX.FASTQ.Record[]
    for d in 1:depth, r in pf_tiling(ref; read_len = read_len)

        push!(rs,
            FASTX.FASTQ.Record("r$(d)_$(FASTX.identifier(r))",
                String(FASTX.sequence(r)), String(FASTX.quality(r))))
    end
    return rs
end

# A single read carrying a unique interior substitution => singleton (count-1) k-mers.
function pf_singleton(ref::AbstractString, pos::Int; read_len::Int = 15, id::String = "err")
    sub = collect(ref[pos:(pos + read_len - 1)])
    m = cld(read_len, 2)
    sub[m] = sub[m] == 'A' ? 'C' : 'A'
    return FASTX.FASTQ.Record(id, String(sub), repeat("I", read_len))
end

# Distinct FORWARD k-mer strings (singlestrand labels are forward, not canonical).
function pf_forward_kmers(seqs, k::Int)
    s = Set{String}()
    for seq in seqs
        str = String(seq)
        for i in 1:(length(str) - k + 1)
            push!(s, str[i:(i + k - 1)])
        end
    end
    return s
end

pf_label_strings(g) = Set(string(l) for l in MetaGraphsNext.labels(g))
pf_nv(g) = length(collect(MetaGraphsNext.labels(g)))

Test.@testset "Qualmer coverage prefilter (min_count)" begin
    k = 7

    Test.@testset "A: floor 2 excludes singleton error k-mers; floor 1 no-op" begin
        reads = pf_deep(PF_REF; depth = 40)
        push!(reads, pf_singleton(PF_REF, 8; id = "e1"))
        push!(reads, pf_singleton(PF_REF, 18; id = "e2"))
        solid = pf_forward_kmers([PF_REF], k)

        g1 = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :singlestrand, min_count = 1)
        g2 = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :singlestrand, min_count = 2)

        # floor 1 = exact no-op: keeps the solid set AND the injected singletons
        Test.@test issubset(solid, pf_label_strings(g1))
        Test.@test pf_nv(g1) > length(solid)
        # floor 2: exactly the distinct solid k-mers, singletons gone
        Test.@test pf_label_strings(g2) == solid
        # every survivor has coverage >= 2
        for l in MetaGraphsNext.labels(g2)
            Test.@test Mycelia.Rhizomorph.get_vertex_observation_count(g2, l) >= 2
        end
    end

    Test.@testset "B: floor 2 bounds memory below floor 1 when errors present" begin
        # Inject many unique-error reads: floor 1 stores all of them; floor 2 drops them.
        reads = pf_deep(PF_REF; depth = 30)
        for j in 1:20
            pos = 1 + (j % (length(PF_REF) - 15))
            push!(reads, pf_singleton(PF_REF, pos; id = "err$(j)"))
        end
        g1 = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :singlestrand, min_count = 1)
        g2 = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :singlestrand, min_count = 2)
        Test.@test pf_nv(g2) < pf_nv(g1)                       # fewer vertices (errors pruned)
        # NOTE: Base.summarysize is NOT asserted here — on a ~30-vertex toy graph it
        # is dominated by fixed MetaGraphsNext/Dict overhead + hash-capacity rounding
        # and does not shrink monotonically with vertex count. The real memory win is
        # a SCALE phenomenon (E. coli: 200GB+ -> bounded), validated by the benchmark
        # memory gate, not a toy unit test. Vertex-count reduction is the reliable proxy.
        # All solid (reference) k-mers survive; g2 keeps only kmers seen >= 2x
        # (a few adjacent error reads can share a boundary k-mer >= 2x, so g2 is a
        # superset of the reference set, not exactly equal — exact equality is
        # covered by testset A with clean, distinct singletons).
        Test.@test issubset(pf_forward_kmers([PF_REF], k), pf_label_strings(g2))
        for l in MetaGraphsNext.labels(g2)
            Test.@test Mycelia.Rhizomorph.get_vertex_observation_count(g2, l) >= 2
        end
    end

    Test.@testset "C: doublestrand stays bounded (<= 2x singlestrand-solid)" begin
        reads = pf_deep(PF_REF; depth = 40)
        push!(reads, pf_singleton(PF_REF, 12; id = "e1"))
        ss = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :singlestrand, min_count = 2)
        ds = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :doublestrand, min_count = 2)
        # doublestrand replicates survivors' RC; bounded by 2x, never re-adds pruned errors.
        Test.@test pf_nv(ss) <= pf_nv(ds) <= 2 * pf_nv(ss)
    end

    Test.@testset "E: regression lock — known singleton absent under floor 2" begin
        reads = pf_deep(PF_REF; depth = 40)
        errrec = pf_singleton(PF_REF, 20; id = "lock")
        push!(reads, errrec)
        err_kmers = setdiff(
            pf_forward_kmers([String(FASTX.sequence(errrec))], k),
            pf_forward_kmers([PF_REF], k)
        )
        Test.@test !isempty(err_kmers)
        g2 = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :singlestrand, min_count = 2)
        got = pf_label_strings(g2)
        for ek in err_kmers
            Test.@test !(ek in got)
        end
    end

    Test.@testset "D: config prefilter is fully opt-in (default 1 = no-op everywhere)" begin
        # The prefilter drops low-coverage k-mers => it CHANGES assembly output, so
        # it is NEVER auto-enabled — default 1 on every path incl. the corrector.
        cfg_default = Mycelia.Rhizomorph.AssemblyConfig(k = k)
        cfg_iter = Mycelia.Rhizomorph.AssemblyConfig(k = k, corrector = :iterative)
        Test.@test cfg_default.qualmer_prefilter_min_count == 1
        Test.@test cfg_iter.qualmer_prefilter_min_count == 1  # NOT auto-enabled (opt-in only)
        # explicit opt-in works
        cfg_opt = Mycelia.Rhizomorph.AssemblyConfig(k = k, qualmer_prefilter_min_count = 2)
        Test.@test cfg_opt.qualmer_prefilter_min_count == 2
        # invalid rejected
        # Invalid rejected with a message that names the field (per repo
        # convention: @test_throws must validate message content).
        local threw = false
        try
            Mycelia.Rhizomorph.AssemblyConfig(k = k, qualmer_prefilter_min_count = 0)
        catch err
            threw = true
            Test.@test occursin("qualmer_prefilter_min_count", sprint(showerror, err))
        end
        Test.@test threw
    end

    Test.@testset "F: doublestrand counts merged (canonical) coverage, not per-strand" begin
        # A unique sequence and its reverse complement: every canonical k-mer has
        # MERGED coverage 2 (once forward in A, once forward in revcomp(A)), but
        # each FORWARD orientation has count 1. This is the strand-split case
        # (CodeRabbit + code-review convergent finding). :doublestrand/:canonical
        # is the DEFAULT corrector production mode with auto-2, so the threshold
        # must apply to merged coverage, not per forward orientation.
        A = "ATCGGCTAATGCCGATTGCACGTA"  # 24 bp, no k=7 repeats
        rcA = String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(A)))
        reads = [
            FASTX.FASTQ.Record("fwd", A, repeat("I", length(A))),
            FASTX.FASTQ.Record("rc", rcA, repeat("I", length(rcA))),
        ]
        # doublestrand: threshold on MERGED coverage => canonical count 2 => KEPT.
        # (Without the mode-aware fix this graph would be empty — the regression.)
        ds = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :doublestrand, min_count = 2)
        Test.@test pf_nv(ds) > 0
        # singlestrand: strands are legitimately distinct => each forward k-mer
        # has count 1 => correctly pruned.
        ss = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :singlestrand, min_count = 2)
        Test.@test pf_nv(ss) == 0
    end
end

println("✓ Qualmer coverage prefilter tests completed")
