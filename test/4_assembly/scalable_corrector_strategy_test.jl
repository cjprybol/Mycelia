# Stage 1 (td-fuo8): the corrector `strategy` fork. The :scalable tier (DEFAULT)
# routes to the coarse/skip engine settings; the :exhaustive tier is the
# maximum-sensitivity EXACT-ML engine (prime-walk / 10 iters / exact UNBOUNDED
# beam / no skip / no hard-window / no soft-EM). NOTE :exhaustive is NOT a
# byte-identical reproduction of the prior corrector default — master's route used
# the bounded auto-beam, so :exhaustive's exact unbounded beam can OOM above that
# threshold; it is for small-scale / high-sensitivity use. The knob mapping is a
# pure function so it can be asserted without running the (slow) corrector, and a
# full end-to-end assemble on a toy proves both tiers reach a real AssemblyResult.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/scalable_corrector_strategy_test.jl")'

import Test
import Mycelia
import FASTX
import BioSequences
import Random

const _BASES = ['A', 'C', 'G', 'T']

function _toy_fastq_records(rng; reflen = 1000, n_reads = 80, readlen = 80, err = 0.01)
    ref = join(rand(rng, _BASES, reflen))
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        for j in 1:readlen
            rand(rng) < err && (seq[j] = rand(rng, filter(!=(seq[j]), _BASES)))
        end
        push!(records, FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
    end
    return records
end

# 1kb quality-gap fixture (mirrors benchmarking/quality_gap_diagnostic.jl): reads
# sampled from BOTH strands and passed through Mycelia.observe at err=0.01, so the
# corrector runs on a realistic double-stranded read set. Used by the td-nt69
# doublestrand-recovery testset below.
function _qgd_reads(; genome_len = 1000, readlen = 100, coverage = 20,
        err = 0.01, seed = 42)
    rec = Mycelia.random_fasta_record(moltype = :DNA, seed = seed, L = genome_len)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
    rng = Random.MersenneTwister(seed)
    glen = length(refseq)
    n_reads = max(1, ceil(Int, coverage * glen / readlen))
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        start = rand(rng, 1:(glen - readlen + 1))
        frag = refseq[start:(start + readlen - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        obs_seq, quals = Mycelia.observe(frag; error_rate = err, tech = :illumina)
        isempty(obs_seq) && continue
        qstr = String([Char(q + 33) for q in quals])
        push!(records, FASTX.FASTQ.Record("read_$(i)", string(obs_seq), qstr))
    end
    return records
end

Test.@testset "scalable corrector strategy fork (td-fuo8)" begin
    R = Mycelia.Rhizomorph

    Test.@testset "strategy knob routing (pure)" begin
        ex = R._corrector_strategy_knobs(:exhaustive)
        # Exhaustive = maximum-sensitivity exact-ML engine (prime-by-prime walk,
        # 10 iterations/k, exact UNBOUNDED beam, none of the new gates). This is
        # the intended hyper-sensitive tier, NOT a reproduction of the prior
        # corrector default (which used the bounded auto-beam).
        Test.@test ex.n_k_rungs === nothing
        Test.@test ex.max_iterations_per_k == 10
        Test.@test ex.skip_solid == false
        Test.@test ex.hard_window == false
        Test.@test ex.windowed_decode == false
        Test.@test ex.soft_em == false
        Test.@test ex.beam_width == typemax(Int)
        # Exhaustive is the explicit exact/oracle tier: it must never inherit the
        # scalable branching/frontier classifier.
        Test.@test ex.indel_schedule == :unrestricted
        # td-nt69: :exhaustive derives its corrector graph_mode from
        # config.graph_mode (nothing ⇒ derive), a byte-identical passthrough.
        Test.@test ex.graph_mode === nothing

        sc = R._corrector_strategy_knobs(:scalable)
        # Scalable = coarse 3-rung ladder, low iteration cap, all volume/quality
        # gates on, size-aware auto-beam (nothing). `soft_em=true` is the ENGINE
        # switch that runs the v2 competing-paths E-step AND the support-floored
        # M-step consumption, which the surfaced provenance flag
        # ("v2-competing-paths-floor") reflects — asserted end-to-end below.
        Test.@test sc.n_k_rungs == 3
        Test.@test sc.max_iterations_per_k == 2
        Test.@test sc.skip_solid == true
        Test.@test sc.hard_window == true
        # Stage 3c (td-nn6l): :scalable now decodes hard reads window-by-window.
        Test.@test sc.windowed_decode == true
        Test.@test sc.soft_em == true
        Test.@test sc.beam_width === nothing
        Test.@test sc.indel_schedule == :frontier_budgeted
        # td-nt69: :scalable runs the corrector on a :doublestrand graph. Forcing
        # :canonical was THE cause of the quality gap (skip machinery was over-
        # constrained to require it); the classification is coverage-based and
        # mode-agnostic, so skip stays active under :doublestrand.
        Test.@test sc.graph_mode == :doublestrand

        Test.@test_throws Exception R._corrector_strategy_knobs(:bogus)
    end

    Test.@testset "AssemblyConfig threads + validates strategy" begin
        # SCALABLE is the default tier.
        Test.@test R.AssemblyConfig(k = 13).strategy == :scalable
        Test.@test R.AssemblyConfig(k = 13, strategy = :exhaustive).strategy == :exhaustive
        Test.@test_throws Exception R.AssemblyConfig(k = 13, strategy = :bogus)
    end

    Test.@testset "both tiers reach a real AssemblyResult" begin
        reads = _toy_fastq_records(Random.MersenneTwister(11))

        sc = R.assemble_genome(reads; k = 13, corrector = :iterative, strategy = :scalable)
        Test.@test sc isa R.AssemblyResult
        Test.@test !isempty(sc.contigs)
        Test.@test sc.assembly_stats["strategy"] == "scalable"
        Test.@test sc.assembly_stats["hard_window"] == true
        Test.@test sc.assembly_stats["hard_read_gate"] == true
        # Per-hard-region windowed decode (td-nn6l Stage 3c) is now ACTIVE on
        # :scalable — hard reads are decoded window-by-window (bounded), not whole.
        Test.@test sc.assembly_stats["windowed_decode"] == true
        # soft-EM v2 is ACTIVE (E-step enumerates competing paths, M-step registers
        # the support-floored soft weights), so the surfaced provenance is the v2
        # marker, never a bare `true`.
        Test.@test sc.assembly_stats["soft_em"] == "v2-competing-paths-floor"
        Test.@test sc.assembly_stats["skip_solid"] == true
        # Skip fraction is a real fraction in [0, 1].
        skip = sc.assembly_stats["skip_fraction"]
        Test.@test 0.0 <= skip <= 1.0

        ex = R.assemble_genome(reads; k = 13, corrector = :iterative, strategy = :exhaustive)
        Test.@test ex isa R.AssemblyResult
        Test.@test !isempty(ex.contigs)
        Test.@test ex.assembly_stats["strategy"] == "exhaustive"
        # Exhaustive threads neither gate nor soft-EM.
        Test.@test ex.assembly_stats["hard_window"] == false
        Test.@test ex.assembly_stats["windowed_decode"] == false
        Test.@test ex.assembly_stats["soft_em"] == false
        Test.@test ex.assembly_stats["skip_solid"] == false
    end

    Test.@testset "default corrector tier is scalable" begin
        reads = _toy_fastq_records(Random.MersenneTwister(11))
        # No explicit strategy ⇒ scalable.
        res = R.assemble_genome(reads; k = 13, corrector = :iterative)
        Test.@test res.assembly_stats["strategy"] == "scalable"
    end

    Test.@testset ":scalable doublestrand recovers low-contig regime (td-nt69)" begin
        # THE quality fix. On the same 1kb fixture the quality-gap diagnostic used
        # (err=0.01, ~20x, 100bp), the previously-forced :canonical corrector
        # shattered the re-assembly to ~197 contigs / N50 43. Running the corrector
        # on :doublestrand (all other :scalable knobs held) recovers the historical
        # near-complete regime (~14-21 contigs / N50 ~900), WITH the skip gate still
        # active. This is an invariant (regime recovery + skip-active), not a golden
        # contig count, so the thresholds are deliberately loose.
        reads = _qgd_reads()
        asm = R.assemble_genome(reads; k = 21, corrector = :iterative, strategy = :scalable)

        n_contigs = length(asm.contigs)
        lens = sort(length.(asm.contigs); rev = true)
        total = sum(lens; init = 0)
        # N50 (contig length at which cumulative coverage reaches half the total).
        n50 = 0
        acc = 0
        for l in lens
            acc += l
            if acc >= total / 2
                n50 = l
                break
            end
        end
        skip_fraction = get(asm.assembly_stats, "skip_fraction", 0.0)

        Test.@info "td-nt69 doublestrand recovery" n_contigs n50 largest = (isempty(lens) ? 0 : lens[1]) skip_fraction skip_solid = asm.assembly_stats["skip_solid"] hard_window = asm.assembly_stats["hard_window"]

        # RECOVERED low-contig regime: nowhere near the ~197 canonical shatter.
        # Loose upper bound (observed ~14) that still fails hard on a regression to
        # the fragmented canonical behavior.
        Test.@test n_contigs <= 40
        # High N50 — the historical near-complete regime was ~891-903; require a
        # large jump over the canonical N50 of 43.
        Test.@test n50 >= 400
        # Skip MUST stay ACTIVE under :doublestrand (the whole point: the fix does
        # not disable the volume-reduction gate, it just runs it on the right graph).
        Test.@test skip_fraction > 0.0
        Test.@test asm.assembly_stats["skip_solid"] == true
        Test.@test asm.assembly_stats["hard_window"] == true
    end
end
