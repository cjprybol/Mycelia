# Regression tests: Rhizomorph must recover a genome from a branchy read graph
# efficiently, on BOTH decoder arms. Reads here are error-laden full-length
# copies of the genome (via create_test_reads/observe) — not fragmented shotgun
# reads — which is enough to produce the error-induced branch vertices both
# defects depend on, without modeling fragment length or indels.
#
# Two distinct defects motivated these tests, neither caught by the prior suite
# (which only assembled 16-28 bp references — graphs too small to exhibit either):
#
#   * k-mer arm (quality-off): recovered the genome but ran the unused, O(n^2)
#     detect_bubbles_next / resolve_repeats_next analysis on every assembly,
#     making it minutes-slow on real read graphs (>22 min on a 48.5 kb phage).
#   * qualmer arm (quality-on): primary path-finding found no substantial paths
#     and fell to a placeholder walk that emits ~one short fragment per vertex,
#     so the largest contig collapsed well below the genome length as size grew.
#
# Both are fixed by extracting maximal unitigs (find_contigs_next) as the primary
# contig source and dropping the dead graph-cleaning calls.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/rhizomorph_read_assembly_recovery_test.jl")'

import Test
import Mycelia
import FASTX
import BioSequences
import Random
import StableRNGs

function _synthetic_reads(genome_len, coverage, error_rate; seed = 42)
    # StableRNGs gives a version-stable genome (the stage-4 suite convention).
    genome = BioSequences.LongDNA{4}(join(
        rand(StableRNGs.StableRNG(seed), ['A', 'C', 'G', 'T'], genome_len)))
    # create_test_reads -> observe injects substitution errors using the GLOBAL
    # RNG (it takes no rng argument), so seed it for a reproducible error topology
    # — the branch/bubble structure both defects depend on. Snapshot and restore
    # the global RNG around the call so this does not leak state into sibling
    # tests sharing the process.
    saved_rng = copy(Random.default_rng())
    Random.seed!(seed)
    reads = Mycelia.create_test_reads(genome, coverage, error_rate)
    copy!(Random.default_rng(), saved_rng)
    return genome, reads
end

Test.@testset "Rhizomorph read assembly recovery" begin
    k = 31

    Test.@testset "qualmer arm recovers genome (quality-on / FASTQ input)" begin
        genome_len = 2000
        _, fastq_reads = _synthetic_reads(genome_len, 20, 0.01)
        result = Mycelia.Rhizomorph.assemble_genome(fastq_reads; k = k, verbose = false)

        Test.@test !isempty(result.contigs)
        largest = maximum(length.(result.contigs); init = 0)
        # Degenerate output collapses largest contig far below genome length.
        Test.@test largest >= 0.8 * genome_len
    end

    Test.@testset "k-mer arm recovers genome without pathological slowdown" begin
        genome_len = 3000
        _, fastq_reads = _synthetic_reads(genome_len, 20, 0.01)
        fasta_reads = [FASTX.FASTA.Record(FASTX.identifier(r), FASTX.sequence(r))
                       for r in fastq_reads]

        elapsed = @elapsed result = Mycelia.Rhizomorph.assemble_genome(
            fasta_reads; k = k, verbose = false)

        largest = maximum(length.(result.contigs); init = 0)
        Test.@test largest >= 0.8 * genome_len
        # The dead O(n^2) cleaning calls make a 3 kb assembly take ~1 min; the
        # unitig path finishes in seconds. Generous bound to avoid CI flake.
        Test.@test elapsed < 30
    end
end
