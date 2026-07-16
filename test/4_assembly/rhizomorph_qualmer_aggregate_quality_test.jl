# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_qualmer_aggregate_quality_test.jl")'
# ```
#
# Tests for corrector-compatible AGGREGATE qualmer quality storage (td-n8ax):
# a coverage counter + running per-position MEAN Phred per DISTINCT k-mer,
# O(distinct k-mers), composed with the opt-in coverage prefilter (PR #425), so
# E. coli-scale assembly stays under a SPAdes-class memory ceiling.
#
# Every behavior is OPT-IN: the default corrector path stays memory_profile=:full
# (per-observation evidence) and byte-identical. These testsets fail until the
# aggregate storage is made corrector-compatible (B1 accessor + B2 exact-mean
# storage + B4 mode/prefilter + wiring).

import Test
import Mycelia
import FASTX
import BioSequences
import MetaGraphsNext
import Graphs

# ---------------------------------------------------------------------------
# Shared fixtures / helpers
# ---------------------------------------------------------------------------

"A short reference with no exact repeats at k=7, so k-mers are distinct."
const AGG_REF = "ATCGGCTAATGCCGATTGCACGTACGTTAGCTAGGCATG"

"""
Build `n` identical FASTQ reads of `seq` with a strictly-ascending per-position
Phred so each k-mer's correct per-position mean is non-uniform (catches position
ordering) AND high enough that a running SUM over `n` observations saturates a
UInt8 accumulator (raw Phred ~40 * n>=7 > 255). This is what discriminates an
exact unclamped mean (correct) from the clamped joint-quality sum (wrong).
"""
function agg_identical_fastq_reads(seq::AbstractString; n::Int = 10)
    qual = join(Char(min(35 + j, 60) + 33) for j in 1:length(seq))  # ~Q36..Q60
    reads = FASTX.FASTQ.Record[]
    for i in 1:n
        push!(reads, FASTX.FASTQ.Record("r$(i)", seq, qual))
    end
    return reads
end

Test.@testset "Rhizomorph aggregate qualmer quality (td-n8ax)" begin

    # -----------------------------------------------------------------------
    # A. EXACT-MEAN EQUIVALENCE: get_vertex_mean_quality on the aggregate
    #    (:lightweight_quality) graph == the :full graph, per position, even at
    #    coverage that saturates a UInt8 running sum. Proves B1 (accessor) + B2
    #    (unclamped storage).
    # -----------------------------------------------------------------------
    Test.@testset "A: exact per-position mean quality matches :full" begin
        k = 7
        reads = agg_identical_fastq_reads(AGG_REF; n = 10)
        fq = Mycelia.Rhizomorph._prepare_fastq_observations(reads)

        full_graph = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :singlestrand, memory_profile = :full)
        agg_graph = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :singlestrand, memory_profile = :lightweight_quality)

        full_labels = Set(MetaGraphsNext.labels(full_graph))
        agg_labels = Set(MetaGraphsNext.labels(agg_graph))
        # Same distinct k-mers (as-observed keys) in both storage modes.
        Test.@test full_labels == agg_labels
        Test.@test !isempty(agg_labels)

        ds = "dataset_01"
        compared = 0
        for label in full_labels
            full_mean = Mycelia.Rhizomorph.get_vertex_mean_quality(full_graph[label], ds)
            agg_mean = Mycelia.Rhizomorph.get_vertex_mean_quality(agg_graph[label], ds)
            # The aggregate MUST expose a real mean, not the nothing/UInt8(2) filler.
            Test.@test agg_mean !== nothing
            Test.@test full_mean !== nothing
            Test.@test length(agg_mean) == k
            # Numerically equal to the per-observation mean (exact, not clamped).
            Test.@test all(isapprox.(agg_mean, full_mean; atol = 1e-9))
            compared += 1
        end
        Test.@test compared == length(full_labels)
    end

    # -----------------------------------------------------------------------
    # B. MODE HONORED: the aggregate profile must build canonical/doublestrand
    #    when asked (today it silently returns singlestrand), AND the exact-mean
    #    accessor must survive the reduced strand conversion (proves the merge
    #    carries dataset_quality_sum with RC reversal). Part of B4.
    # -----------------------------------------------------------------------
    Test.@testset "B: mode honored + mean survives strand conversion" begin
        k = 7
        reads = agg_identical_fastq_reads(AGG_REF; n = 4)
        fq = Mycelia.Rhizomorph._prepare_fastq_observations(reads)

        ss = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :singlestrand, memory_profile = :lightweight_quality)
        dsg = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :doublestrand, memory_profile = :lightweight_quality)
        cn = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :canonical, memory_profile = :lightweight_quality)

        n_ss = length(collect(MetaGraphsNext.labels(ss)))
        n_ds = length(collect(MetaGraphsNext.labels(dsg)))
        n_cn = length(collect(MetaGraphsNext.labels(cn)))

        # Doublestrand adds reverse-complement vertices (NOT a silent singlestrand).
        Test.@test n_ds > n_ss
        Test.@test n_ds <= 2 * n_ss
        Test.@test n_cn <= n_ss

        # The reduced doublestrand vertices still expose an exact mean (the merge
        # must carry dataset_quality_sum through conversion, reversing for RC).
        for label in MetaGraphsNext.labels(dsg)
            m = Mycelia.Rhizomorph.get_vertex_mean_quality(dsg[label], "dataset_01")
            Test.@test m !== nothing
            Test.@test length(m) == k
        end
    end

    # -----------------------------------------------------------------------
    # C. PREFILTER COMPOSES: the PR #425 coverage prefilter (min_count) must
    #    thread into the aggregate builders so singleton error k-mers are pruned
    #    BEFORE aggregate vertex creation — the composition that bounds memory.
    #    Part of B4.
    # -----------------------------------------------------------------------
    Test.@testset "C: coverage prefilter composes with aggregate storage" begin
        k = 7
        # Solid backbone at 8x + one read carrying a unique substitution, whose
        # k-mers spanning the edit appear exactly once (singleton error k-mers).
        solid = agg_identical_fastq_reads(AGG_REF; n = 8)
        # Single-base substitution at pos 20, guaranteed different from the original.
        orig_base = AGG_REF[20]
        new_base = orig_base == 'A' ? 'C' : 'A'
        err_seq = AGG_REF[1:19] * string(new_base) * AGG_REF[21:end]
        Test.@test err_seq != AGG_REF
        qual = join(Char(min(35 + j, 60) + 33) for j in 1:length(err_seq))
        push!(solid, FASTX.FASTQ.Record("err", err_seq, qual))
        fq = Mycelia.Rhizomorph._prepare_fastq_observations(solid)

        g1 = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :singlestrand, memory_profile = :lightweight_quality, min_count = 1)
        g2 = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :singlestrand, memory_profile = :lightweight_quality, min_count = 2)

        n1 = length(collect(MetaGraphsNext.labels(g1)))
        n2 = length(collect(MetaGraphsNext.labels(g2)))
        # min_count=1 is a no-op superset; floor 2 prunes the singleton error k-mers.
        Test.@test n2 < n1
        # Every survivor under floor 2 has coverage >= 2.
        for label in MetaGraphsNext.labels(g2)
            Test.@test g2[label].total_count >= 2
        end
    end
end

println("✓ Rhizomorph aggregate qualmer quality tests completed")
