# Convergence + k-ladder tuning for mycelia_iterative_assemble (td-q70n).
#
# Guards two literature-backed speedups for the iterative corrector:
#
#   1. No-change convergence stop (Musket): a pass that makes 0 improvements
#      breaks the per-k loop immediately, INDEPENDENT of improvement_threshold.
#      Isolated here by driving threshold to 0.0 (so the threshold-based early
#      stop can never fire) on perfectly clean reads (0 improvements guaranteed):
#        - stop_on_no_change=true  => exactly 1 iteration at the first k
#        - stop_on_no_change=false => runs the full max_iterations_per_k
#
#   2. Coarse k-ladder (LoRMA): build_k_ladder / _next_k_in_progression drive a
#      small, well-spaced progression instead of every prime, while nothing =>
#      legacy prime-by-prime behavior.

import Test
import Mycelia
import FASTX
import Random

const CBASES = ['A', 'C', 'G', 'T']

# Clean reads (err=0) so the corrector makes exactly 0 improvements every pass.
function write_clean_fastq(path, rng; reflen = 800, n_reads = 60, readlen = 80)
    ref = join(rand(rng, CBASES, reflen))
    open(path, "w") do io
        for i in 1:n_reads
            s = rand(rng, 1:(reflen - readlen + 1))
            seq = ref[s:(s + readlen - 1)]
            println(io, "@r$i"); println(io, seq); println(io, "+")
            println(io, String(fill('I', readlen)))
        end
    end
    return path
end

Test.@testset "iterative corrector: no-change convergence + k-ladder" begin

    Test.@testset "build_k_ladder" begin
        # Coarse geometric ladder: 3 rungs from initial to max.
        Test.@test Mycelia.build_k_ladder(3, 31; n_k_rungs = 3) == [3, 11, 31]
        # Explicit ladder is filtered to [initial_k, max_k] and sorted/deduped.
        Test.@test Mycelia.build_k_ladder(3, 31; k_ladder = [15, 5, 9, 100, 5]) == [5, 9, 15]
        # No knob set => nothing (legacy prime progression preserved).
        Test.@test Mycelia.build_k_ladder(3, 31) === nothing
        # First rung pinned to initial_k, last rung <= max_k and odd.
        ladder = Mycelia.build_k_ladder(5, 40; n_k_rungs = 4)
        Test.@test first(ladder) == 5
        Test.@test last(ladder) <= 40
        Test.@test isodd(last(ladder))
        Test.@test issorted(ladder)
    end

    Test.@testset "_next_k_in_progression" begin
        # Scheduled ladder: advance to next strictly-larger rung.
        Test.@test Mycelia._next_k_in_progression(3, 31, [3, 11, 31]) == 11
        Test.@test Mycelia._next_k_in_progression(11, 31, [3, 11, 31]) == 31
        # Top of ladder is a fixed point (main loop treats == current as stop).
        Test.@test Mycelia._next_k_in_progression(31, 31, [3, 11, 31]) == 31
        # nothing => legacy prime progression.
        Test.@test Mycelia._next_k_in_progression(5, 10, nothing) == 7
    end

    Test.@testset "no-change stop is independent of improvement_threshold" begin
        dir = mktempdir()
        fastq = write_clean_fastq(joinpath(dir, "clean.fastq"), Random.MersenneTwister(7))

        # threshold=0.0 disarms the threshold-based early stop entirely, so any
        # bound on per-k iterations must come from the no-change stop.
        with_stop = Mycelia.mycelia_iterative_assemble(fastq;
            max_k = 7, max_iterations_per_k = 4, improvement_threshold = 0.0,
            stop_on_no_change = true, verbose = false,
            enable_checkpointing = false, output_dir = joinpath(dir, "with_stop"))

        without_stop = Mycelia.mycelia_iterative_assemble(fastq;
            max_k = 7, max_iterations_per_k = 4, improvement_threshold = 0.0,
            stop_on_no_change = false, verbose = false,
            enable_checkpointing = false, output_dir = joinpath(dir, "without_stop"))

        first_k = first(with_stop[:k_progression])

        # With the no-change stop, a 0-improvement first pass ends that k at once.
        Test.@test length(with_stop[:metadata][:iteration_history][first_k]) == 1
        # Without it (and with threshold=0.0), the same clean data burns the full
        # per-k iteration budget doing no useful work.
        Test.@test length(without_stop[:metadata][:iteration_history][first_k]) == 4
        # The no-change stop therefore does strictly fewer total iterations while
        # producing the same assembly (clean reads => no corrections either way).
        Test.@test with_stop[:metadata][:total_iterations] <
                   without_stop[:metadata][:total_iterations]
        Test.@test with_stop[:final_assembly] == without_stop[:final_assembly]
    end

    Test.@testset "coarse ladder completes with fewer k-steps than prime walk" begin
        dir = mktempdir()
        fastq = write_clean_fastq(joinpath(dir, "clean2.fastq"), Random.MersenneTwister(11))

        legacy = Mycelia.mycelia_iterative_assemble(fastq;
            max_k = 31, max_iterations_per_k = 2, n_k_rungs = nothing,
            verbose = false, enable_checkpointing = false,
            output_dir = joinpath(dir, "legacy"))

        coarse = Mycelia.mycelia_iterative_assemble(fastq;
            max_k = 31, max_iterations_per_k = 2, n_k_rungs = 3,
            verbose = false, enable_checkpointing = false,
            output_dir = joinpath(dir, "coarse"))

        # Both complete past the first k-transition.
        Test.@test length(legacy[:k_progression]) >= 2
        Test.@test length(coarse[:k_progression]) >= 2
        # The coarse ladder walks strictly fewer k-mer sizes.
        Test.@test length(coarse[:k_progression]) < length(legacy[:k_progression])
        # Assembly is not degraded (identical corrected read set on clean input).
        Test.@test coarse[:final_assembly] == legacy[:final_assembly]
    end
end
