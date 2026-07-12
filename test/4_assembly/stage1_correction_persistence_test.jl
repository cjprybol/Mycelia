# td-ohob: Stage-1 corrected-FASTQ persistence (hybrid-OLC plumbing).
#
# `_assemble_with_iterative_corrector` used to run the corrector inside an
# `mktempdir` that its `finally` block deleted, so the corrected reads never
# outlived the call — a public accessor for just the corrected reads did not
# exist. The hybrid-OLC route (Stage-2 route (a)) needs those corrected reads to
# survive so an external OLC assembler can consume them. The reusable helper
# `_run_stage1_correction(reads, config)` now materializes the corrected reads AND
# persists the corrected FASTQ to a caller-owned location.
#
# These tests prove:
#   (1) With config.output_dir set, the corrected FASTQ exists AFTER the call
#       returns (it outlived the corrector's temp-dir cleanup) at
#       joinpath(output_dir, "corrected.fastq"), and round-trips (re-reading it
#       yields the same reads the helper returned in memory).
#   (2) With config.output_dir unset (nothing), the helper still returns an
#       existing FASTQ, but at an ephemeral tempname() path the CALLER owns — this
#       is the path the native re-assembly deletes, preserving the historical
#       no-stray-file behavior.
#   (3) The end-to-end assemble_genome(...; corrector=:iterative) path still
#       returns a real AssemblyResult with stamped corrector provenance (the
#       thin-consumer wiring from helper -> re-assembly tail is intact).
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/stage1_correction_persistence_test.jl")'

import Test
import Mycelia
import FASTX
import BioSequences
import Random

const R = Mycelia.Rhizomorph

# Small synthetic Illumina-style reads from a random reference (the recipe the
# corrector fixtures use). Kept intentionally small — the iterative corrector
# costs 5-10x naive, so the genome/coverage are sized for a fast unit test.
function _stage1_reads(; genome_len = 600, readlen = 100, coverage = 8,
        err = 0.02, seed = 7)
    rec = Mycelia.random_fasta_record(moltype = :DNA, seed = seed, L = genome_len)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
    rng = Random.MersenneTwister(seed)
    glen = length(refseq)
    erl = min(readlen, glen)
    n = max(1, ceil(Int, coverage * glen / erl))
    recs = FASTX.FASTQ.Record[]
    for i in 1:n
        s = glen == erl ? 1 : rand(rng, 1:(glen - erl + 1))
        frag = refseq[s:(s + erl - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        obs, q = Mycelia.observe(frag; error_rate = err, tech = :illumina)
        isempty(obs) && continue
        push!(recs, FASTX.FASTQ.Record("read_$i", string(obs),
            String([Char(x + 33) for x in q])))
    end
    return recs
end

# Re-read a FASTQ into memory (mirrors the eager collect the helper uses).
function _read_fastq(path)
    open(FASTX.FASTQ.Reader, path) do reader
        collect(reader)
    end
end

Test.@testset "Stage-1 correction persistence (td-ohob)" begin
    # Mycelia.observe draws error/quality from the GLOBAL RNG (bare rand()), so
    # seed it for reproducibility — matching reassembly_graph_reuse_test.jl's
    # convention. The local MersenneTwister in _stage1_reads only governs fragment
    # start + revcomp; the seed here pins the rest.
    Random.seed!(1234)
    reads = _stage1_reads()
    Test.@test !isempty(reads)

    Test.@testset "persist to caller-owned output_dir + round-trip" begin
        outdir = mktempdir()
        config = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, output_dir = outdir)
        res = R._run_stage1_correction(reads, config)

        # Persistence: the corrected FASTQ outlived the corrector's temp-dir cleanup.
        Test.@test isfile(res.corrected_fastq)
        Test.@test dirname(res.corrected_fastq) == outdir
        Test.@test basename(res.corrected_fastq) == "corrected.fastq"
        # A caller-supplied output_dir is NOT ephemeral — the caller keeps the file.
        Test.@test res.ephemeral == false
        # The returned metadata path is repointed at the persisted file, not the
        # now-deleted original inside the corrector temp dir.
        Test.@test res.result_dict[:metadata][:final_fastq_file] == res.corrected_fastq

        # Non-empty corrected read set (the 0-read guard would have thrown).
        Test.@test !isempty(res.corrected_reads)

        # Round-trip: re-reading the persisted FASTQ yields the same reads the
        # helper returned in memory (same count + same sequences, in order).
        reread = _read_fastq(res.corrected_fastq)
        Test.@test length(reread) == length(res.corrected_reads)
        Test.@test [FASTX.sequence(String, r) for r in reread] ==
                   [FASTX.sequence(String, r) for r in res.corrected_reads]

        # The tail-consumed pieces are present so the native path needs no recompute.
        Test.@test haskey(res.result_dict, :metadata)
        Test.@test res.max_k >= 13
        Test.@test hasproperty(res.knobs, :skip_solid)
    end

    Test.@testset "ephemeral tempfile when output_dir unset" begin
        config = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable)  # output_dir defaults to nothing
        Test.@test config.output_dir === nothing

        res = R._run_stage1_correction(reads, config)
        # Still an existing FASTQ on return (the helper does NOT delete it — the
        # native caller is what cleans an ephemeral tempfile).
        Test.@test isfile(res.corrected_fastq)
        Test.@test endswith(res.corrected_fastq, ".fastq")
        Test.@test res.ephemeral == true
        Test.@test !isempty(res.corrected_reads)

        # We own it here (as the native caller would); clean up.
        rm(res.corrected_fastq; force = true)
        Test.@test !isfile(res.corrected_fastq)
    end

    Test.@testset "end-to-end corrector wiring intact" begin
        # The thin consumer (_assemble_with_iterative_corrector) still routes helper
        # output through the re-assembly tail and returns a real AssemblyResult with
        # stamped corrector provenance.
        config = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable)
        assembly = R.assemble_genome(reads, config)
        Test.@test assembly isa R.AssemblyResult
        Test.@test assembly.assembly_stats["corrector"] == "iterative"
        Test.@test assembly.assembly_stats["corrected_read_count"] > 0
    end

    Test.@testset "corrected FASTQ survives a full native assembly (output_dir set)" begin
        # The headline behavior + the caller's conditional cleanup: a full
        # assemble_genome(...; corrector=:iterative, output_dir=set) must LEAVE
        # corrected.fastq in place for the Stage-2 OLC handoff. This guards the
        # `if stage1.ephemeral` branch in _assemble_with_iterative_corrector — a
        # regression that deleted unconditionally would strand the handoff, and the
        # helper-level tests above would not catch it.
        outdir = mktempdir()
        config = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, output_dir = outdir)
        assembly = R.assemble_genome(reads, config)
        Test.@test assembly isa R.AssemblyResult
        # Persisted output survived the full native re-assembly + its finally block.
        persisted = joinpath(outdir, "corrected.fastq")
        Test.@test isfile(persisted)
        # And it round-trips to a non-empty read set the external assembler can use.
        Test.@test !isempty(_read_fastq(persisted))
    end

    Test.@testset "output_dir constructor validation" begin
        # Empty-string output_dir is rejected at construction (it would otherwise
        # take the persist branch and silently write corrected.fastq into cwd).
        Test.@test_throws Exception R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, output_dir = "")
        # A non-empty output_dir (dir need not exist yet) constructs fine.
        Test.@test R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, output_dir = joinpath(mktempdir(), "nested")) isa
                   R.AssemblyConfig
    end
end
