# Scale-model probe: does k_ref sequence-space SATURATION break the error estimator?
#
# The hypothesis is that k_ref=13 (space 4^13 = 67.1M) saturates on real-sized inputs.
# At E. coli scale (5 Mb, 30x, 10% ONT error) the reads need ~200M distinct 13-mers,
# ~3x the whole space, so unrelated loci must collide, counts inflate, solid_fraction
# rises toward 1, and the estimate collapses toward 0.
#
# Simulating 5 Mb at 30x on a laptop is out of scope. Instead this is a SCALE MODEL:
# hold the genome small and shrink the k-mer SPACE by lowering k_ref, which reproduces
# the same saturation ratio at negligible cost. If the estimate collapses as
# distinct_kmers/4^k_ref approaches 1, the mechanism is confirmed and the k_ref=13
# default is unsafe at real genome scale.
#
# This is a probe, not a deliverable. Ground truth is exact (uniform substitutions).

import Mycelia
import Random
import Printf

const GENOME_LENGTH = 20_000
const READ_LENGTH = 5_000
const BASES = ['A', 'C', 'G', 'T']
const TRUE_ERROR = 0.05
const COVERAGE = 30
const SEED = 42
const K_REFS = [5, 6, 7, 8, 9, 11, 13, 17, 21]

function apply_substitutions(rng, sequence::String, error_rate::Float64)::String
    out = Char[]
    sizehint!(out, length(sequence))
    for character in sequence
        if rand(rng) < error_rate
            push!(out, rand(rng, filter(b -> b != character, BASES)))
        else
            push!(out, character)
        end
    end
    return String(out)
end

function main()
    rng = Random.MersenneTwister(SEED)
    reference = String([rand(rng, BASES) for _ = 1:GENOME_LENGTH])
    read_count = cld(COVERAGE * GENOME_LENGTH, READ_LENGTH)
    reads = String[]
    for _ = 1:read_count
        start_index = rand(rng, 1:(GENOME_LENGTH-READ_LENGTH+1))
        fragment = reference[start_index:(start_index+READ_LENGTH-1)]
        push!(reads, apply_substitutions(rng, fragment, TRUE_ERROR))
    end

    Printf.@printf(
        "true error = %.3f, coverage = %d, genome = %d bp, %d reads\n\n",
        TRUE_ERROR,
        COVERAGE,
        GENOME_LENGTH,
        length(reads)
    )
    Printf.@printf(
        "%6s %14s %14s %10s %10s %10s\n",
        "k_ref",
        "space 4^k",
        "distinct",
        "saturation",
        "estimate",
        "bias"
    )

    for k_ref in K_REFS
        # distinct k-mers actually observed, using the same hashing the estimator uses
        seen = Set{UInt64}()
        for read in reads
            characters = collect(read)
            for start_index = 1:(length(characters)-k_ref+1)
                push!(
                    seen,
                    UInt64(hash(view(characters, start_index:(start_index+k_ref-1)))),
                )
            end
        end
        space = 4.0^k_ref
        estimate = Mycelia.Rhizomorph._kmer_spectrum_residual_error(reads, k_ref, 2)
        Printf.@printf(
            "%6d %14.0f %14d %10.4f %10.5f %+10.5f\n",
            k_ref,
            space,
            length(seen),
            length(seen) / space,
            estimate,
            estimate - TRUE_ERROR
        )
    end

    println("\nInterpretation: saturation -> 1 means the k-mer space is exhausted, so")
    println("distinct loci necessarily share k-mers, counts inflate, and the estimate")
    println("collapses toward 0. Watch whether the estimate degrades as saturation rises.")
end

main()
