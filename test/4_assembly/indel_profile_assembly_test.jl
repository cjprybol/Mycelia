# Step 4a (indel-wiring): drive the merged indel pair-HMM decoder from an ERROR
# PROFILE (sequencing technology), not read length or a hardcoded tech. The
# indel kernel (`_viterbi_correct_observation_indel`) is reachable only via
# `correct_observations`; this suite proves `assemble_genome` now enables it for
# indel-prone technologies (nanopore/pacbio) while keeping the DEFAULT
# (:illumina) on the substitution-only oracle path byte-for-byte.
#
# The profile -> fractions -> indel_moves mapping and the threshold gate are pure
# functions (asserted without running the slow corrector); two end-to-end
# assemblies prove the nanopore path enables indels and the illumina path stays
# byte-identical to a run with no new parameter.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/indel_profile_assembly_test.jl")'

import Test
import Mycelia
import FASTX
import BioSequences
import Random

# Nanopore/illumina read set sampled from a random reference through
# `Mycelia.observe(...)` at the requested tech (returns a (seq, quals) TUPLE).
function _profile_reads(; tech::Symbol = :nanopore, genome_len = 1000, readlen = 100,
        coverage = 20, err = 0.05, seed = 42)
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
        obs = Mycelia.observe(frag; error_rate = err, tech = tech)
        obs_seq = obs[1]
        quals = obs[2]
        isempty(obs_seq) && continue
        qstr = String([Char(q + 33) for q in quals])
        push!(records, FASTX.FASTQ.Record("read_$(i)", string(obs_seq), qstr))
    end
    return records
end

Test.@testset "indel-aware correction wired via sequencing-tech error profile" begin
    R = Mycelia.Rhizomorph

    Test.@testset "error-profile presets mirror simulation fractions (pure)" begin
        np = Mycelia.indel_error_profile(:nanopore)
        # Nanopore conditional insertion/deletion fractions mirror observe()'s
        # per-tech split (0.30 / 0.30 near simulation.jl:1931).
        Test.@test np.insertion_fraction ≈ 0.30
        Test.@test np.deletion_fraction ≈ 0.30
        Test.@test np.insertion_extend_probability > 0.0
        Test.@test np.deletion_extend_probability > 0.0

        il = Mycelia.indel_error_profile(:illumina)
        # Illumina's ABSOLUTE indel rate is negligible; the correction-side preset
        # carries ~0 indel fractions so the decoder collapses to substitution-only.
        Test.@test il.insertion_fraction + il.deletion_fraction < 0.02

        pb = Mycelia.indel_error_profile(:pacbio)
        Test.@test pb.insertion_fraction + pb.deletion_fraction > 0.02

        Test.@test_throws Exception Mycelia.indel_error_profile(:bogus)
    end

    Test.@testset "threshold gate: profile -> indel_moves" begin
        # Above threshold => indel-aware; negligible => substitution-only.
        Test.@test Mycelia.profile_enables_indels(:nanopore) == true
        Test.@test Mycelia.profile_enables_indels(:pacbio) == true
        Test.@test Mycelia.profile_enables_indels(:illumina) == false
        Test.@test Mycelia.profile_enables_indels(:ultima) == false
        # The gate is a strict threshold on the summed indel fractions.
        Test.@test Mycelia.profile_enables_indels(:nanopore; threshold = 0.7) == false
        Test.@test Mycelia.profile_enables_indels(:illumina; threshold = 0.0) == false
    end

    Test.@testset "AssemblyConfig threads + validates sequencing_tech" begin
        # Default is :illumina (unchanged behavior).
        Test.@test R.AssemblyConfig(k = 13).sequencing_tech == :illumina
        Test.@test R.AssemblyConfig(k = 13, sequencing_tech = :nanopore).sequencing_tech ==
                   :nanopore
        Test.@test_throws Exception R.AssemblyConfig(k = 13, sequencing_tech = :bogus)
    end

    Test.@testset "tier knobs carry indel bounds" begin
        sc = R._corrector_strategy_knobs(:scalable)
        # :scalable => bounded band + small run caps.
        Test.@test sc.band_width isa Int
        Test.@test sc.band_width > 0
        Test.@test sc.deletion_max_run == 3
        Test.@test sc.max_insertion_run >= 1

        ex = R._corrector_strategy_knobs(:exhaustive)
        # :exhaustive => unbounded band + larger caps.
        Test.@test ex.band_width === nothing
        Test.@test ex.deletion_max_run > sc.deletion_max_run
    end

    Test.@testset "illumina default is the substitution-only oracle (byte-identical)" begin
        reads = _profile_reads(tech = :illumina, err = 0.01, seed = 7)
        base = R.assemble_genome(reads; k = 13, corrector = :iterative)
        il = R.assemble_genome(reads; k = 13, corrector = :iterative,
            sequencing_tech = :illumina)
        # Contigs byte-identical: the illumina profile threads NO indel params, so
        # the corrector runs the substitution decode unchanged.
        Test.@test il.contigs == base.contigs
        Test.@test base.assembly_stats["indel_moves"] == false
        Test.@test il.assembly_stats["indel_moves"] == false
        Test.@test il.assembly_stats["sequencing_tech"] == "illumina"
    end

    Test.@testset "nanopore enables indel-aware correction" begin
        reads = _profile_reads(tech = :nanopore, err = 0.05, seed = 7)
        np = R.assemble_genome(reads; k = 13, corrector = :iterative,
            sequencing_tech = :nanopore)
        Test.@test np isa R.AssemblyResult
        Test.@test !isempty(np.contigs)
        # The provenance stamp proves the indel decode path was engaged.
        Test.@test np.assembly_stats["indel_moves"] == true
        Test.@test np.assembly_stats["sequencing_tech"] == "nanopore"
    end
end
