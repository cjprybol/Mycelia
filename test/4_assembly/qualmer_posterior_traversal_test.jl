# H-11 P2: traversal_weighting=:posterior (prune + fallback only).
#
# Design note: docs/h11-p1-p3-design-note.md §§E1, E3, P2.
#
# - Default remains `:evidence`; `:evidence` / `:quality` paths must stay
#   bit-identical to pre-P2 behaviour (this file pins non-identity for
#   `:posterior` and quality-invariance under `:evidence`).
# - Constant-Q control (E1): same sequences rewritten to Q40 / Q2; oracle vs
#   constant vs antioracle change contig digests or prune sets under `:posterior`.
# - Non-identity (E3): on a fixture where qualities discriminate, `:posterior`
#   output differs from `:evidence`.
# - v1 scope: quality influences prune + fallback walk only — not Viterbi /
#   decoder transitions (documented on AssemblyConfig / prune helpers).

import Test
import StableRNGs
import BioSequences
import FASTX
import Mycelia
import SHA

const POST_BASES = [BioSequences.DNA_A, BioSequences.DNA_C,
    BioSequences.DNA_G, BioSequences.DNA_T]

"""Rewrite every base quality to a constant raw Phred (FASTX UInt8 API)."""
function _constant_q_records(records, raw_phred::Integer)
    q = UInt8(raw_phred)
    return [FASTX.FASTQ.Record(
                String(FASTX.identifier(r)),
                FASTX.sequence(BioSequences.LongDNA{4}, r),
                fill(q, length(FASTX.sequence(r))))
            for r in records]
end

"""Oracle Q: errors Q2, clean Q40. Antioracle: swapped."""
function _oracle_antioracle_records(reads)
    oracle = [FASTX.FASTQ.Record("read_$(i)", r.sequence,
                  UInt8[c ? 0x02 : 0x28 for c in r.corrupted])
              for (i, r) in enumerate(reads)]
    antioracle = [FASTX.FASTQ.Record("read_$(i)", r.sequence,
                      UInt8[c ? 0x28 : 0x02 for c in r.corrupted])
                  for (i, r) in enumerate(reads)]
    return oracle, antioracle
end

sorted_contigs(result) = sort(String.(result.contigs))

function contig_digest(result)
    return bytes2hex(SHA.sha256(join(sorted_contigs(result), '\n')))
end

Test.@testset "H-11 P2 posterior traversal (E1/E3)" begin
    rng = StableRNGs.StableRNG(1122)
    genome = BioSequences.LongDNA{4}(rand(rng, POST_BASES, 600))
    k = 11
    read_length = 120
    coverage = 12

    reads = NamedTuple{
        (:sequence, :corrupted), Tuple{BioSequences.LongDNA{4}, Vector{Bool}}}[]
    n_reads = (length(genome) * coverage) ÷ read_length
    for _ in 1:n_reads
        start = rand(rng, 1:(length(genome) - read_length + 1))
        bases = collect(genome[start:(start + read_length - 1)])
        corrupted = falses(read_length)
        for position in 1:read_length
            if rand(rng) < 0.02
                bases[position] = rand(rng, filter(!=(bases[position]), POST_BASES))
                corrupted[position] = true
            end
        end
        push!(reads,
            (sequence = BioSequences.LongDNA{4}(bases), corrupted = collect(corrupted)))
    end

    oracle_records, antioracle_records = _oracle_antioracle_records(reads)
    constant40 = _constant_q_records(oracle_records, 40)
    constant02 = _constant_q_records(oracle_records, 2)

    Test.@testset "sanity: sequences match; qualities differ across conditions" begin
        Test.@test all(FASTX.sequence(String, oracle_records[i]) ==
                       FASTX.sequence(String, constant40[i]) ==
                       FASTX.sequence(String, constant02[i]) ==
                       FASTX.sequence(String, antioracle_records[i])
        for i in eachindex(oracle_records))
        Test.@test any(collect(FASTX.quality_scores(oracle_records[i])) !=
                       collect(FASTX.quality_scores(constant40[i]))
        for i in eachindex(oracle_records))
        Test.@test any(collect(FASTX.quality_scores(constant40[i])) !=
                       collect(FASTX.quality_scores(constant02[i]))
        for i in eachindex(oracle_records))
    end

    Test.@testset "AssemblyConfig accepts :posterior; defaults unchanged" begin
        cfg = Mycelia.Rhizomorph.AssemblyConfig(; k = k, traversal_weighting = :posterior)
        Test.@test cfg.traversal_weighting == :posterior
        Test.@test cfg.traversal_min_log_odds == 0.0
        Test.@test cfg.traversal_min_quality == 60.0

        default = Mycelia.Rhizomorph.AssemblyConfig(; k = k)
        Test.@test default.traversal_weighting == :evidence
        Test.@test default.traversal_min_log_odds == 0.0

        cfg2 = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, traversal_weighting = :posterior, traversal_min_log_odds = 5.0)
        Test.@test cfg2.traversal_min_log_odds == 5.0

        bad = try
            Mycelia.Rhizomorph.AssemblyConfig(;
                k = k, traversal_weighting = :posterior,
                traversal_min_log_odds = NaN)
            nothing
        catch error
            error
        end
        Test.@test bad !== nothing
        Test.@test occursin("traversal_min_log_odds must be finite",
            sprint(showerror, bad))
    end

    Test.@testset "E1: constant-Q / oracle / antioracle change under :posterior" begin
        kw = (; k = k, verbose = false, traversal_weighting = :posterior,
            traversal_min_log_odds = 0.0)
        oracle = Mycelia.Rhizomorph.assemble_genome(oracle_records; kw...)
        c40 = Mycelia.Rhizomorph.assemble_genome(constant40; kw...)
        c02 = Mycelia.Rhizomorph.assemble_genome(constant02; kw...)
        anti = Mycelia.Rhizomorph.assemble_genome(antioracle_records; kw...)

        for result in (oracle, c40, c02, anti)
            Test.@test result.assembly_stats["traversal_weighting"] == "posterior"
            Test.@test result.assembly_stats["quality_influences_traversal"] == true
        end

        digests = (
            contig_digest(oracle), contig_digest(c40), contig_digest(c02),
            contig_digest(anti))
        prunes = (
            oracle.assembly_stats["traversal_quality_pruned_vertices"],
            c40.assembly_stats["traversal_quality_pruned_vertices"],
            c02.assembly_stats["traversal_quality_pruned_vertices"],
            anti.assembly_stats["traversal_quality_pruned_vertices"])

        # Quality channel is load-bearing: not all four conditions agree.
        Test.@test length(unique(digests)) > 1 || length(unique(prunes)) > 1

        # Constant Q40 vs Q2 must differ (positive control that quality, not
        # sequence, is the variable) — digests or prune counts.
        Test.@test contig_digest(c40) != contig_digest(c02) ||
                   c40.assembly_stats["traversal_quality_pruned_vertices"] !=
                   c02.assembly_stats["traversal_quality_pruned_vertices"]

        # Oracle vs antioracle: swapped quality stories change the assembly.
        Test.@test contig_digest(oracle) != contig_digest(anti) ||
                   oracle.assembly_stats["traversal_quality_pruned_vertices"] !=
                   anti.assembly_stats["traversal_quality_pruned_vertices"]

        Test.@test !isempty(oracle.contigs)
        Test.@test maximum(length.(oracle.contigs)) > k
    end

    Test.@testset "E3: :posterior differs from :evidence when qualities discriminate" begin
        # Pinned recurrent low-Q artifact (mirrors qualmer_quality_channel_test.jl):
        # random substitution errors alone often leave Λ ≥ 0 for every vertex at the
        # default floor, so prune is a no-op and arms are identical by design. A
        # fixed locus observed many times at Q2 forces Λ below the floor.
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
        honest = [FASTX.FASTQ.Record("art_$(i)", r.sequence,
                      UInt8[p in r.artifact_positions ? 0x02 : 0x28
                            for p in 1:length(r.sequence)])
                  for (i, r) in enumerate(artifact_reads)]

        evidence = Mycelia.Rhizomorph.assemble_genome(
            honest; k = k, verbose = false, traversal_weighting = :evidence)
        posterior = Mycelia.Rhizomorph.assemble_genome(
            honest; k = k, verbose = false, traversal_weighting = :posterior,
            traversal_min_log_odds = 0.0)

        Test.@test evidence.assembly_stats["quality_influences_traversal"] == false
        Test.@test evidence.assembly_stats["traversal_quality_pruned_vertices"] == 0
        Test.@test posterior.assembly_stats["quality_influences_traversal"] == true
        Test.@test posterior.assembly_stats["traversal_quality_pruned_vertices"] > 0
        Test.@test sorted_contigs(evidence) != sorted_contigs(posterior)
        Test.@test !isempty(posterior.contigs)
    end

    Test.@testset ":evidence stays quality-invariant (bit-identical across Q rewrite)" begin
        e_oracle = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false)
        e_c02 = Mycelia.Rhizomorph.assemble_genome(
            constant02; k = k, verbose = false)
        Test.@test sorted_contigs(e_oracle) == sorted_contigs(e_c02)
        Test.@test e_oracle.assembly_stats["traversal_weighting"] == "evidence"
    end

    Test.@testset ":quality path still accepted and distinct from :posterior API" begin
        q = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false, traversal_weighting = :quality)
        p = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false, traversal_weighting = :posterior)
        Test.@test q.assembly_stats["traversal_weighting"] == "quality"
        Test.@test p.assembly_stats["traversal_weighting"] == "posterior"
        # Different floors / combiners — not required to differ on every fixture,
        # but both must report quality influence.
        Test.@test q.assembly_stats["quality_influences_traversal"]
        Test.@test p.assembly_stats["quality_influences_traversal"]
    end

    Test.@testset "opt-in survives the corrector route" begin
        corrected = Mycelia.Rhizomorph.assemble_genome(
            oracle_records; k = k, verbose = false, corrector = :iterative,
            traversal_weighting = :posterior, traversal_min_log_odds = 0.0)
        Test.@test corrected.assembly_stats["traversal_weighting"] == "posterior"
        Test.@test corrected.assembly_stats["quality_influences_traversal"] == true
        Test.@test corrected.assembly_stats["traversal_min_log_odds"] == 0.0
    end
end
