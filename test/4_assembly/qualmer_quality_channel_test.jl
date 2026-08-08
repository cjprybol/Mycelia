# Quality-channel tests for the qualmer assembly arm (bead td-4e19d.2).
#
# Three things are pinned here, and conflating them is the mistake this file exists to
# prevent.
#
#   1. CHARACTERIZATION. Under the `:evidence` default, per-base quality reaches contig
#      EMISSION only and never a traversal decision, so the qualmer arm returns
#      byte-identical contigs no matter what the quality scores are. That is a deliberate
#      design choice (PR #335 moved the arm onto the k-mer arm's `find_contigs_next` to
#      fix degenerate assemblies), so it is pinned as CURRENT BEHAVIOUR — if it ever
#      changes, every committed benchmark result silently changes with it.
#
#   2. THE OPT-IN FIX. Under `traversal_weighting = :quality`, quality DOES change the
#      assembly. Written against the default this assertion fails.
#
#   3. THE DESIGN'S ACTUAL CLAIM — that compounding observations of HIGH-QUALITY k-mers
#      is worth more than compounding observations alone
#      (`planning-docs/rhizomorph-graph-ecosystem-plan.md:498-550`). Tests 1 and 2 cannot
#      see this, because they inject only RANDOM errors, and random errors do not recur.
#      The distinction the framework rests on appears only against a RECURRENT artifact:
#      something observed many times, at low quality. An evidence count cannot tell that
#      from real sequence; summed Phred can. The third testset builds exactly that case.
#
# Method throughout: hold the read SEQUENCES fixed and vary ONLY the quality vector, so
# any output difference is attributable to quality alone.
#
# Broader multi-chemistry evidence lives in
# `benchmarking/qualmer_quality_channel_probe.jl`; this file is the fast in-CI guard.

import Test
import StableRNGs
import BioSequences
import FASTX
import Mycelia

const QC_BASES = [BioSequences.DNA_A, BioSequences.DNA_C,
    BioSequences.DNA_G, BioSequences.DNA_T]

Test.@testset "Qualmer quality channel (td-4e19d.2)" begin
    rng = StableRNGs.StableRNG(4219)

    # A genome complex enough that k=11 k-mers are mostly unique, so the assembly is
    # structured rather than one degenerate blob.
    genome = BioSequences.LongDNA{4}(rand(rng, collect("ACGT"), 600))
    k = 11
    read_length = 120
    coverage = 12

    # Reads carry random substitution errors. Error POSITIONS are recorded so the
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
                bases[position] = rand(rng, filter(!=(bases[position]), QC_BASES))
                corrupted[position] = true
            end
        end
        push!(reads,
            (sequence = BioSequences.LongDNA{4}(bases), corrupted = collect(corrupted)))
    end

    # `oracle`: quality is CORRECT by construction (errors Q2, clean bases Q40).
    # `flat`:   every base Q40 regardless of whether it is an error.
    oracle_records = [FASTX.FASTQ.Record("read_$(i)", r.sequence,
                          UInt8[c ? 0x02 : 0x28 for c in r.corrupted])
                      for (i, r) in enumerate(reads)]
    flat_records = [FASTX.FASTQ.Record("read_$(i)", r.sequence,
                        fill(UInt8(0x28), length(r.sequence)))
                    for (i, r) in enumerate(reads)]

    sorted_contigs(result) = sort(String.(result.contigs))

    Test.@testset "sanity: the conditions differ ONLY in quality" begin
        Test.@test all(FASTX.sequence(String, oracle_records[i]) ==
                       FASTX.sequence(String, flat_records[i])
        for i in eachindex(oracle_records))
        # If the quality strings were equal, every comparison below would be vacuous.
        Test.@test any(collect(FASTX.quality_scores(oracle_records[i])) !=
                       collect(FASTX.quality_scores(flat_records[i]))
        for i in eachindex(oracle_records))
    end

    Test.@testset "characterization: :evidence default is quality-invariant" begin
        oracle_assembly = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false)
        flat_assembly = Mycelia.Rhizomorph.assemble_genome(
            flat_records; k = k, verbose = false)

        Test.@test sorted_contigs(oracle_assembly) == sorted_contigs(flat_assembly)
        Test.@test oracle_assembly.assembly_stats["traversal_weighting"] == "evidence"
        Test.@test oracle_assembly.assembly_stats["quality_influences_traversal"] == false
        Test.@test oracle_assembly.assembly_stats["traversal_quality_pruned_vertices"] == 0

        # The k-mer arm (same reads, quality stripped to FASTA) reaches the SAME contigs.
        # This is why the Track-A `qualmer` and `kmer` decoder arms agree on every
        # assembly metric.
        fasta_records = [FASTX.FASTA.Record("read_$(i)", r.sequence)
                         for (i, r) in enumerate(reads)]
        kmer_assembly = Mycelia.Rhizomorph.assemble_genome(
            fasta_records; k = k, verbose = false)
        Test.@test sorted_contigs(kmer_assembly) == sorted_contigs(oracle_assembly)
    end

    Test.@testset "opt-in: traversal_weighting=:quality makes quality load-bearing" begin
        oracle_assembly = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false, traversal_weighting = :quality)
        flat_assembly = Mycelia.Rhizomorph.assemble_genome(
            flat_records; k = k, verbose = false, traversal_weighting = :quality)

        Test.@test oracle_assembly.assembly_stats["traversal_weighting"] == "quality"
        Test.@test oracle_assembly.assembly_stats["quality_influences_traversal"] == true

        # THE ASSERTION THAT FAILS WITHOUT THE FIX: identical reads, different quality,
        # different assembly. Under `:evidence` these are byte-identical (above).
        Test.@test sorted_contigs(oracle_assembly) != sorted_contigs(flat_assembly)

        # The mechanism is the quality gate, not incidental churn.
        Test.@test oracle_assembly.assembly_stats["traversal_quality_pruned_vertices"] >
                   flat_assembly.assembly_stats["traversal_quality_pruned_vertices"]

        # QUALITY, not coverage. Under uniform Q40 the joint score is exactly
        # 40 x n_observations, so a 60.0 floor IS "seen at least twice" — meaning the
        # assertion above is satisfiable by a mechanism containing no quality
        # information at all. Pin the distinction: the ORACLE condition (where quality
        # varies per base) must NOT be reproducible by a pure coverage prefilter.
        coverage_only = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false, qualmer_prefilter_min_count = 2)
        Test.@test sorted_contigs(oracle_assembly) != sorted_contigs(coverage_only)

        # And the opt-in must not be achieving its difference by collapsing the
        # assembly: `!=` alone would pass if `:quality` pruned the graph to nothing.
        Test.@test !isempty(oracle_assembly.contigs)
        Test.@test maximum(length.(oracle_assembly.contigs)) > k
    end

    Test.@testset "opt-in survives the corrector route" begin
        # The corrector rebuilds its re-assembly config from a hand-listed kwarg
        # subset, which silently dropped both new options: requesting `:quality`
        # yielded `evidence` with no error and a stats dict that confirmed the wrong
        # thing. This assertion fails against that behaviour.
        corrected = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false, corrector = :iterative,
            traversal_weighting = :quality)
        Test.@test corrected.assembly_stats["traversal_weighting"] == "quality"
        Test.@test corrected.assembly_stats["quality_influences_traversal"] == true
    end

    Test.@testset "joint quality: recurrent low-quality artifact is rejected" begin
        # THE DESIGN'S CENTRAL CLAIM. A technical artifact recurs — so evidence COUNT
        # cannot distinguish it from real sequence — but it recurs at LOW quality, so
        # summed Phred can. Random errors cannot test this because they do not recur.
        #
        # Construct a systematic artifact: one fixed locus, one fixed wrong base, in
        # EVERY read covering it, always low quality. Its evidence count is therefore as
        # high as any real k-mer's.
        artifact_rng = StableRNGs.StableRNG(99)
        artifact_locus = 300
        artifact_base = genome[artifact_locus] == BioSequences.DNA_A ?
                        BioSequences.DNA_C : BioSequences.DNA_A

        artifact_reads = NamedTuple{
            (:sequence, :artifact_positions),
            Tuple{BioSequences.LongDNA{4}, Vector{Int}}}[]
        for _ in 1:n_reads
            start = rand(artifact_rng, 1:(length(genome) - read_length + 1))
            bases = collect(genome[start:(start + read_length - 1)])
            artifact_positions = Int[]
            offset = artifact_locus - start + 1
            if 1 <= offset <= read_length
                bases[offset] = artifact_base
                push!(artifact_positions, offset)
            end
            push!(artifact_reads,
                (sequence = BioSequences.LongDNA{4}(bases),
                    artifact_positions = artifact_positions))
        end

        # Same sequences, two quality stories: the artifact is either honestly reported
        # as low quality, or indistinguishable from real signal.
        honest = [FASTX.FASTQ.Record("art_$(i)", r.sequence,
                      UInt8[p in r.artifact_positions ? 0x02 : 0x28
                            for p in 1:length(r.sequence)])
                  for (i, r) in enumerate(artifact_reads)]
        credulous = [FASTX.FASTQ.Record("art_$(i)", r.sequence,
                         fill(UInt8(0x28), length(r.sequence)))
                     for (i, r) in enumerate(artifact_reads)]

        # Every k-mer spanning the corrupted locus, in both orientations (contigs may be
        # emitted on either strand). Full k-mers, not a short window: a 5-mer probe
        # matches elsewhere in a 600 bp genome by chance and would report the artifact
        # as present in every assembly.
        artifact_window = copy(genome[(artifact_locus - k + 1):(artifact_locus + k - 1)])
        artifact_window[k] = artifact_base
        artifact_kmers = String[]
        for i in 1:k
            kmer = BioSequences.LongDNA{4}(artifact_window[i:(i + k - 1)])
            push!(artifact_kmers, String(kmer))
            push!(artifact_kmers, String(BioSequences.reverse_complement(kmer)))
        end
        contains_artifact(result) = any(
            contig -> any(kmer -> occursin(kmer, contig), artifact_kmers),
            String.(result.contigs))

        honest_assembly = Mycelia.Rhizomorph.assemble_genome(
            honest; k = k, verbose = false, traversal_weighting = :quality)
        credulous_assembly = Mycelia.Rhizomorph.assemble_genome(
            credulous; k = k, verbose = false, traversal_weighting = :quality)
        evidence_assembly = Mycelia.Rhizomorph.assemble_genome(
            honest; k = k, verbose = false)

        # Joint quality compounds across observations, so honest low-quality evidence
        # stays below the floor no matter how many times the artifact is seen...
        Test.@test !contains_artifact(honest_assembly)
        # ...while the SAME recurring sequence survives when its quality is not believed.
        Test.@test contains_artifact(credulous_assembly)
        # ...and survives under count-based weighting, which is exactly the failure mode
        # the quality channel is supposed to fix.
        Test.@test contains_artifact(evidence_assembly)
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

        # The floor lives on the JOINT Phred scale [0, 255], not the single-observation
        # 0-60 scale, so 99.0 must be ACCEPTED and only >255 rejected.
        Test.@test Mycelia.Rhizomorph.AssemblyConfig(;
            k = 11, traversal_weighting = :quality,
            traversal_min_quality = 99.0).traversal_min_quality == 99.0

        bad_floor = try
            Mycelia.Rhizomorph.AssemblyConfig(;
                k = 11, traversal_weighting = :quality, traversal_min_quality = 300.0)
            nothing
        catch error
            error
        end
        Test.@test bad_floor !== nothing
        Test.@test occursin("traversal_min_quality must be between 0.0 and 255.0",
            sprint(showerror, bad_floor))
    end
end
