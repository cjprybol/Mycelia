# Stage 1 (td-fuo8): the corrector `strategy` fork. The :scalable tier (DEFAULT)
# routes to the coarse/skip engine settings; the :exhaustive tier reproduces
# today's hyper-sensitive engine byte-for-byte (prime-walk / 10 iters / exact
# beam / no skip / no hard-window / no soft-EM). The knob mapping is a pure
# function so it can be asserted without running the (slow) corrector, and a
# full end-to-end assemble on a toy proves both tiers reach a real AssemblyResult.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/scalable_corrector_strategy_test.jl")'

import Test
import Mycelia
import FASTX
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

Test.@testset "scalable corrector strategy fork (td-fuo8)" begin
    R = Mycelia.Rhizomorph

    Test.@testset "strategy knob routing (pure)" begin
        ex = R._corrector_strategy_knobs(:exhaustive)
        # Exhaustive = today's hyper-sensitive engine (byte-identical passthrough):
        # prime-by-prime walk, 10 iterations/k, exact beam, none of the new gates.
        Test.@test ex.n_k_rungs === nothing
        Test.@test ex.max_iterations_per_k == 10
        Test.@test ex.skip_solid == false
        Test.@test ex.hard_window == false
        Test.@test ex.soft_em == false
        Test.@test ex.beam_width == typemax(Int)

        sc = R._corrector_strategy_knobs(:scalable)
        # Scalable = coarse 3-rung ladder, low iteration cap, all volume/quality
        # gates on, size-aware auto-beam (nothing).
        Test.@test sc.n_k_rungs == 3
        Test.@test sc.max_iterations_per_k == 2
        Test.@test sc.skip_solid == true
        Test.@test sc.hard_window == true
        Test.@test sc.soft_em == true
        Test.@test sc.beam_width === nothing

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
        Test.@test sc.assembly_stats["soft_em"] == true
        Test.@test sc.assembly_stats["skip_solid"] == true
        # Skip fraction is a real fraction in [0, 1].
        skip = sc.assembly_stats["skip_fraction"]
        Test.@test 0.0 <= skip <= 1.0

        ex = R.assemble_genome(reads; k = 13, corrector = :iterative, strategy = :exhaustive)
        Test.@test ex isa R.AssemblyResult
        Test.@test !isempty(ex.contigs)
        Test.@test ex.assembly_stats["strategy"] == "exhaustive"
        # Exhaustive threads neither gate (byte-identical passthrough).
        Test.@test ex.assembly_stats["hard_window"] == false
        Test.@test ex.assembly_stats["soft_em"] == false
        Test.@test ex.assembly_stats["skip_solid"] == false
    end

    Test.@testset "default corrector tier is scalable" begin
        reads = _toy_fastq_records(Random.MersenneTwister(11))
        # No explicit strategy ⇒ scalable.
        res = R.assemble_genome(reads; k = 13, corrector = :iterative)
        Test.@test res.assembly_stats["strategy"] == "scalable"
    end
end
