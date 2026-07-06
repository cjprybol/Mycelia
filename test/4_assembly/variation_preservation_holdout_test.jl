# VARIATION-PRESERVATION HOLDOUT (td-h6w9)
# =========================================
#
# Mandatory anti-overcorrection gate for the scalable read corrector: correction
# must RETAIN data-supported genetic variation and remove only UNSUPPORTED
# errors. It must never "clean" a real, balanced polymorphism down to a single
# allele. This is a holdout/invariant test — it encodes a property the corrector
# must satisfy, not a golden output it happens to produce today.
#
# CONTROLLED FIXTURE (known truth, no randomness in the biology being tested):
#   * A REAL heterozygous site: two haplotypes (allele A, allele B) identical
#     across a 300 bp backbone except for one substitution at the center. Each
#     allele is simulated at 15x (balanced 50/50, 30x total at the site). BOTH
#     are real and MUST survive correction — collapsing either is destruction of
#     variation.
#   * A SEQUENCING ERROR: a single coverage-1 read carrying a substitution at a
#     second, well-covered locus (30x true consensus vs 1x error). The error is
#     unsupported and is the thing correction is allowed to remove.
#
# All reads are error-free except the one injected error, so the truth is exact:
# any allele-distinguishing k-mer that disappears was removed by the corrector,
# not by simulated noise. Reads are Phred 'I' (Q40) so decisions are coverage-
# driven, matching how the corrector treats quality-less input.
#
# WHAT IS ASSERTED (and why at these strengths):
#   1. VARIATION PRESERVED under strategy=:scalable — both allele k-mers survive,
#      simultaneously. HARD assertion. If this fails the corrector is collapsing
#      a real balanced bubble to one branch: a CRITICAL over-correction defect
#      that BLOCKS the merge. This is the load-bearing invariant of the gate.
#   2. NO REGRESSION vs the naive (corrector=:none) baseline — any allele the
#      naive assembly retained must still be present after :scalable correction.
#      HARD. Attributes blame correctly: it isolates variation LOST by
#      correction from variation the naive de Bruijn path never represented.
#   3. ERROR HANDLING — the coverage-1 error IS removable by the machinery:
#      strategy=:exhaustive (hyper-sensitive tier) removes it (HARD positive
#      control, proving the corrector can distinguish error from signal). The
#      conservative strategy=:scalable tier does NOT yet remove it, because
#      skip-solid skips the mostly-solid error read; that is recorded with
#      @test_broken so the suite stays green now and FLIPS to a failure prompt
#      when v2 competing-paths / aggressive cleaning starts removing it (the
#      point at which this should be promoted to a hard @test).
#
# DELIBERATELY SCOPED OUT (kept honest rather than flaky):
#   * Paralog/repeat retention: near-identical repeat copies collapse under any
#     single-k de Bruijn assembly (naive included) — that is expected assembly
#     behavior, not a corrector property, so a contig-content assertion there
#     would test the assembler, not the invariant this gate guards.
#   * Heuristic-popping negative control: no heuristic bubble-popper is reachable
#     from assemble_genome's public surface on this branch, so it is not wired.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/variation_preservation_holdout_test.jl")'

import Test
import Mycelia
import FASTX
import BioSequences
import Random

const _VPH_BASES = ['A', 'C', 'G', 'T']

# Reverse-complement of a plain string (contigs are Vector{String}).
_vph_revcomp(s::AbstractString) =
    string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(s)))

# True iff `sub` (or its reverse complement) occurs in ANY contig. DoubleStrand
# re-assembly can emit either orientation, so both are checked; case-folded so a
# lowercase-masked contig can never spuriously hide a match.
function _vph_present(contigs, sub::AbstractString)
    rc = _vph_revcomp(sub)
    return any(c -> occursin(sub, uppercase(c)) || occursin(rc, uppercase(c)), contigs)
end

"""
Build the controlled diploid-plus-error fixture. Returns the read set and the
four allele-distinguishing / locus-defining k-mers whose presence encodes the
invariant. Deterministic (fixed seed) so the assertions are reproducible.
"""
function _vph_build_fixture(; k::Int = 21, seed::Int = 20260706)
    rng = Random.MersenneTwister(seed)
    L = 300
    backbone = collect(join(rand(rng, _VPH_BASES, L)))

    # Heterozygous site at the center (spanned by every read below).
    pvar = 150
    base_a = backbone[pvar]                                   # allele A == backbone
    base_b = rand(rng, filter(!=(base_a), _VPH_BASES))
    hap_a = copy(backbone)
    hap_b = copy(backbone); hap_b[pvar] = base_b

    # Independent, well-covered locus carrying the injected sequencing error.
    perr = 100
    base_true = backbone[perr]
    base_err = rand(rng, filter(!=(base_true), _VPH_BASES))
    err_hap = copy(backbone); err_hap[perr] = base_err

    half = k ÷ 2
    kmer_a    = uppercase(String(hap_a[(pvar - half):(pvar + half)]))   # allele-A signature
    kmer_b    = uppercase(String(hap_b[(pvar - half):(pvar + half)]))   # allele-B signature
    kmer_err  = uppercase(String(err_hap[(perr - half):(perr + half)])) # error signature (1x)
    kmer_true = uppercase(String(backbone[(perr - half):(perr + half)]))# true consensus (30x)

    readlen = 150
    qual = String(fill('I', readlen))
    reads = FASTX.FASTQ.Record[]
    ncov = 15  # per allele -> 15x each, 30x total at the het site (balanced 50/50)
    for i in 1:ncov
        sa = rand(rng, 1:100)  # start in 1:100 => read spans BOTH perr(100) and pvar(150)
        push!(reads, FASTX.FASTQ.Record("A$i", String(hap_a[sa:(sa + readlen - 1)]), qual))
        sb = rand(rng, 1:100)
        push!(reads, FASTX.FASTQ.Record("B$i", String(hap_b[sb:(sb + readlen - 1)]), qual))
    end
    # Single coverage-1 error read (backbone with the injected error at perr).
    push!(reads, FASTX.FASTQ.Record("ERR", String(err_hap[1:readlen]), qual))

    return (; reads, k, kmer_a, kmer_b, kmer_err, kmer_true, base_a, base_b, base_true, base_err)
end

Test.@testset "variation-preservation holdout (td-h6w9)" begin
    R = Mycelia.Rhizomorph
    fx = _vph_build_fixture()

    # Naive (uncorrected) baseline: the variation the plain de Bruijn path
    # represents, used to attribute any loss specifically to correction.
    naive = R.assemble_genome(fx.reads; k = fx.k, corrector = :none)
    naive_a = _vph_present(naive.contigs, fx.kmer_a)
    naive_b = _vph_present(naive.contigs, fx.kmer_b)

    # The corrector under test.
    scalable = R.assemble_genome(fx.reads; k = fx.k, corrector = :iterative, strategy = :scalable)

    sc_a    = _vph_present(scalable.contigs, fx.kmer_a)
    sc_b    = _vph_present(scalable.contigs, fx.kmer_b)
    sc_err  = _vph_present(scalable.contigs, fx.kmer_err)
    sc_true = _vph_present(scalable.contigs, fx.kmer_true)

    @info "variation-preservation holdout evidence" naive_a naive_b sc_a sc_b sc_err sc_true strategy = scalable.assembly_stats["strategy"] skip_fraction = get(scalable.assembly_stats, "skip_fraction", nothing)

    Test.@testset "sanity: naive baseline represents both real alleles" begin
        # If the naive path already dropped an allele the fixture is degenerate and
        # testset 2's differential guard would be vacuous; assert the premise.
        Test.@test naive_a
        Test.@test naive_b
    end

    Test.@testset "KEY INVARIANT: :scalable correction preserves both real alleles" begin
        # Both balanced alleles survive correction, simultaneously — the bubble is
        # NOT collapsed to a single branch. Failure => CRITICAL over-correction
        # that destroys real variation and BLOCKS the merge.
        Test.@test sc_a
        Test.@test sc_b
        Test.@test sc_a && sc_b            # both at once (not "one branch won")
        # The real consensus at the error locus is also retained.
        Test.@test sc_true
    end

    Test.@testset "no regression vs naive baseline (differential)" begin
        # Correction must not lose an allele the naive assembly kept.
        Test.@test !(naive_a && !sc_a)     # allele A not lost by correction
        Test.@test !(naive_b && !sc_b)     # allele B not lost by correction
    end

    Test.@testset "error handling: removable by machinery; :scalable tradeoff documented" begin
        # Positive control: the hyper-sensitive tier DOES remove the coverage-1
        # error (proves the corrector can tell error from balanced signal), while
        # still keeping both real alleles.
        exhaustive = R.assemble_genome(fx.reads; k = fx.k, corrector = :iterative, strategy = :exhaustive)
        ex_a   = _vph_present(exhaustive.contigs, fx.kmer_a)
        ex_b   = _vph_present(exhaustive.contigs, fx.kmer_b)
        ex_err = _vph_present(exhaustive.contigs, fx.kmer_err)
        Test.@test ex_a
        Test.@test ex_b
        Test.@test !ex_err                 # exhaustive removes the unsupported error

        # Conservative :scalable tier does NOT yet remove the error, because
        # skip-solid skips the mostly-solid error read. Recorded as broken so the
        # suite is green now and FLIPS when v2 aggressive cleaning starts removing
        # it — at which point promote this to `Test.@test !sc_err`.
        Test.@test_broken !sc_err
    end
end
