# Wiring test: the iterative+skip corrector must be reachable through the public
# assemble_genome interface via the opt-in `corrector=:iterative` keyword, while
# the default path (`corrector=:none`) stays on the existing single-k pipeline.
#
# Guards td-k6oc: previously mycelia_iterative_assemble was only callable from
# benchmarks; AssemblyConfig now threads `corrector`/`skip_solid` and
# assemble_genome routes to it. This test proves (a) the opt-in path runs to
# completion and returns non-empty contigs, and (b) the default path is
# untouched (no corrector metadata, still produces contigs).
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/rhizomorph_iterative_corrector_wiring_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _BASES = ['A', 'C', 'G', 'T']

# ~50-100 fragmented FASTQ reads from a ~1 kb reference (mirrors the known-good
# toy input used by iterative_assemble_completion_test.jl, but returns records so
# they flow through the assemble_genome front door rather than a file path).
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

Test.@testset "assemble_genome iterative corrector wiring (td-k6oc)" begin
    reads = _toy_fastq_records(Random.MersenneTwister(42))

    Test.@testset "config threads corrector/skip_solid + validates" begin
        cfg = Mycelia.Rhizomorph.AssemblyConfig(k = 13, corrector = :iterative, skip_solid = true)
        Test.@test cfg.corrector == :iterative
        Test.@test cfg.skip_solid == true
        # default preserves today's behavior
        Test.@test Mycelia.Rhizomorph.AssemblyConfig(k = 13).corrector == :none
        Test.@test Mycelia.Rhizomorph.AssemblyConfig(k = 13).skip_solid == false
        # invalid corrector is rejected
        Test.@test_throws Exception Mycelia.Rhizomorph.AssemblyConfig(k = 13, corrector = :bogus)
    end

    Test.@testset "corrector=:iterative runs end-to-end, non-empty contigs" begin
        result = Mycelia.Rhizomorph.assemble_genome(reads; k = 13, corrector = :iterative)
        Test.@test result isa Mycelia.Rhizomorph.AssemblyResult
        Test.@test !isempty(result.contigs)
        Test.@test all(!isempty, result.contigs)
        # corrector provenance is recorded
        Test.@test get(result.assembly_stats, "corrector", nothing) == "iterative"
        Test.@test haskey(result.assembly_stats, "k_progression")
    end

    Test.@testset "default (corrector=:none) path is unchanged" begin
        result = Mycelia.Rhizomorph.assemble_genome(reads; k = 13)
        Test.@test result isa Mycelia.Rhizomorph.AssemblyResult
        Test.@test !isempty(result.contigs)
        # default path never touches the corrector -> no corrector metadata
        Test.@test !haskey(result.assembly_stats, "corrector")
    end
end
