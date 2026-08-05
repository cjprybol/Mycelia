# Quality-channel tests for the qualmer assembly arm (bead td-4e19d.2).
#
# Two DIFFERENT things are pinned here, and conflating them is the mistake this file
# exists to prevent:
#
#   1. CHARACTERIZATION. Under the `:evidence` default, per-base quality reaches contig
#      EMISSION only and never a traversal decision, so the qualmer arm returns
#      byte-identical contigs no matter what the quality scores are. That is a deliberate
#      design choice (PR #335 moved the arm onto the k-mer arm's `find_contigs_next` to
#      fix degenerate assemblies), so it is pinned as CURRENT BEHAVIOUR — if it ever
#      changes, every committed benchmark result silently changes with it.
#
#   2. THE OPT-IN FIX. Under `traversal_weighting = :quality`, quality DOES change the
#      assembly. Written against the default this assertion fails; it passes only with
#      the option wired through.
#
# Method: hold the read SEQUENCES fixed and vary ONLY the quality vector. Any output
# difference is then attributable to quality alone.
#
# The full multi-chemistry, multi-coverage evidence lives in
# `benchmarking/qualmer_quality_channel_probe.jl`; this file is the fast in-CI guard.

import Test
import StableRNGs
import BioSequences
import FASTX
import Mycelia

Test.@testset "Qualmer quality channel (td-4e19d.2)" begin
    rng = StableRNGs.StableRNG(4219)

    # A small genome with enough complexity that k=11 k-mers are mostly unique, so the
    # assembly is structured rather than a single degenerate blob.
    genome = BioSequences.LongDNA{4}(rand(rng, collect("ACGT"), 600))
    k = 11
    read_length = 120
    coverage = 12

    # Reads carry substitution errors. The error POSITIONS are recorded so the two
    # quality conditions below can be built from the same ground truth.
    reads = NamedTuple{
        (:sequence, :corrupted), Tuple{BioSequences.LongDNA{4}, Vector{Bool}}}[]
    n_reads = (length(genome) * coverage) ÷ read_length
    for _ in 1:n_reads
        start = rand(rng, 1:(length(genome) - read_length + 1))
        bases = collect(genome[start:(start + read_length - 1)])
        corrupted = falses(read_length)
        for position in 1:read_length
            if rand(rng) < 0.02
                alternatives = filter(!=(bases[position]),
                    [BioSequences.DNA_A, BioSequences.DNA_C,
                        BioSequences.DNA_G, BioSequences.DNA_T])
                bases[position] = rand(rng, alternatives)
                corrupted[position] = true
            end
        end
        push!(reads,
            (sequence = BioSequences.LongDNA{4}(bases), corrupted = collect(corrupted)))
    end

    # `oracle`: quality is CORRECT by construction (errors get Q2, clean bases Q40).
    # `flat`:   every base Q40 regardless of whether it is an error.
    # Identical sequences, maximally different quality.
    oracle_records = [FASTX.FASTQ.Record("read_$(i)", r.sequence,
                          UInt8[c ? 0x02 : 0x28 for c in r.corrupted])
                      for (i, r) in enumerate(reads)]
    flat_records = [FASTX.FASTQ.Record("read_$(i)", r.sequence,
                        fill(UInt8(0x28), length(r.sequence)))
                    for (i, r) in enumerate(reads)]

    sorted_contigs(result) = sort(String.(result.contigs))

    Test.@testset "sanity: the two conditions really do differ only in quality" begin
        Test.@test length(oracle_records) == length(flat_records)
        Test.@test all(FASTX.sequence(String, oracle_records[i]) ==
                       FASTX.sequence(String, flat_records[i])
        for i in eachindex(oracle_records))
        # If the quality strings were equal the comparisons below would be vacuous.
        Test.@test any(collect(FASTX.quality_scores(oracle_records[i])) !=
                       collect(FASTX.quality_scores(flat_records[i]))
        for i in eachindex(oracle_records))
    end

    Test.@testset "characterization: :evidence default is quality-invariant" begin
        oracle_assembly = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false)
        flat_assembly = Mycelia.Rhizomorph.assemble_genome(
            flat_records; k = k, verbose = false)

        # The documented, deliberate behaviour: quality does not reach traversal.
        Test.@test sorted_contigs(oracle_assembly) == sorted_contigs(flat_assembly)
        Test.@test oracle_assembly.assembly_stats["traversal_weighting"] == "evidence"
        Test.@test oracle_assembly.assembly_stats["quality_influences_traversal"] == false
        Test.@test oracle_assembly.assembly_stats["traversal_quality_pruned_vertices"] == 0

        # The k-mer arm (same reads, quality stripped to FASTA) reaches the SAME
        # contigs — this is why the Track-A `qualmer` and `kmer` decoder arms agree on
        # every assembly metric.
        fasta_records = [FASTX.FASTA.Record("read_$(i)", r.sequence)
                         for (i, r) in enumerate(reads)]
        kmer_assembly = Mycelia.Rhizomorph.assemble_genome(
            fasta_records; k = k, verbose = false)
        Test.@test sorted_contigs(kmer_assembly) == sorted_contigs(oracle_assembly)
    end

    Test.@testset "opt-in: traversal_weighting=:quality makes quality load-bearing" begin
        oracle_assembly = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false,
            traversal_weighting = :quality, traversal_min_quality = 20.0)
        flat_assembly = Mycelia.Rhizomorph.assemble_genome(
            flat_records; k = k, verbose = false,
            traversal_weighting = :quality, traversal_min_quality = 20.0)

        Test.@test oracle_assembly.assembly_stats["traversal_weighting"] == "quality"
        Test.@test oracle_assembly.assembly_stats["quality_influences_traversal"] == true

        # THE ASSERTION THAT FAILS WITHOUT THE FIX: identical reads, different quality,
        # different assembly. Under the `:evidence` default these are byte-identical.
        Test.@test sorted_contigs(oracle_assembly) != sorted_contigs(flat_assembly)

        # And the mechanism is the quality gate, not incidental churn: the oracle
        # condition marks error-bearing k-mers low-quality, so it must prune vertices
        # that the uniformly-Q40 condition keeps.
        Test.@test oracle_assembly.assembly_stats["traversal_quality_pruned_vertices"] > 0
        Test.@test flat_assembly.assembly_stats["traversal_quality_pruned_vertices"] == 0
    end

    Test.@testset "option validation" begin
        bad_weighting = try
            Mycelia.Rhizomorph.AssemblyConfig(; k = 11, traversal_weighting = :phred)
            nothing
        catch error
            error
        end
        Test.@test bad_weighting !== nothing
        Test.@test occursin("traversal_weighting must be :evidence or :quality",
            sprint(showerror, bad_weighting))

        bad_floor = try
            Mycelia.Rhizomorph.AssemblyConfig(;
                k = 11, traversal_weighting = :quality, traversal_min_quality = 99.0)
            nothing
        catch error
            error
        end
        Test.@test bad_floor !== nothing
        Test.@test occursin("traversal_min_quality must be between 0.0 and 60.0",
            sprint(showerror, bad_floor))
    end
end
