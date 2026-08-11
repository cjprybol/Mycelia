# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/8_tool_integration/minimap_split_temp_cleanup.jl")'
# ```
#
# minimap2 writes one `PREFIX.NNNN.tmp` chunk per index part when run against a
# multi-part index, and unlinks them itself on clean exit. They survive only when
# the process is killed first (walltime, OOM, scancel) -- and nothing can consume
# them afterwards, since minimap2 cannot resume from them. A 2026-07 CAMI_I_LOW
# run aborted this way and stranded 783 GiB on NERSC scratch.
#
# The first three testsets need no external tools -- they exercise the prefix
# derivation and the reclaim helper against synthetic files -- so they run on the
# default CI path. The fourth calls `minimap_map`, which invokes
# `add_bioconda_env` before it builds any command string, and so self-gates
# behind MYCELIA_RUN_EXTERNAL=true.
#
# Note what a `finally` alone CANNOT cover: SLURM sends SIGTERM then SIGKILL at
# walltime, and Julia runs no `finally` on either -- it dies in the signal
# handler (verified on 1.10.10). That is why `minimap_merge_map_and_split`
# sweeps BEFORE the run as well as after; the next attempt is the only moment a
# live Julia process exists after a job-level kill.

import Test
import Mycelia

Test.@testset "minimap2 split-index temp cleanup" begin
    Test.@testset "split prefix derives from the output path" begin
        Test.@test Mycelia.minimap_split_prefix("/x/y.sorted.bam") ==
                   "/x/y.sorted.bam.tmp"
    end

    Test.@testset "reclaims chunks, spares siblings" begin
        mktempdir() do dir
            outfile = joinpath(dir, "sample.ref.mmi.minimap2.sorted.bam")
            split_prefix = Mycelia.minimap_split_prefix(outfile)

            # Four minimap2 split chunks, 1 KiB each. The last one carries a
            # FIVE-digit index on purpose: minimap2 formats the chunk number
            # with %04d, which stops being four digits once a run exceeds 9999
            # index parts. Without this fixture a literal `\d{4}` predicate
            # passes the whole testset, and the choice of `\d+` is untested.
            chunks = [string(split_prefix, ".", lpad(i, 4, '0'), ".tmp") for i in 0:2]
            push!(chunks, string(split_prefix, ".10000.tmp"))
            for c in chunks
                write(c, rand(UInt8, 1024))
            end

            # Siblings that MUST survive. Two groups, and the second is the one
            # that gives this test teeth:
            #
            #   (a) files that share the OUTPUT stem -- the BAM, its index, and
            #       samtools' own sort temps. A loose `startswith(name, base)`
            #       predicate already spares these, so on their own they prove
            #       nothing about the shape anchor.
            #
            #   (b) files that DO start with `split_prefix` but are not chunks.
            #       These are the discriminating cases: a loose prefix sweep
            #       deletes them, the shape-anchored predicate does not. Without
            #       them this testset passes under either implementation and so
            #       tests nothing about the choice between them.
            #
            #   (c) a chunk belonging to a DIFFERENT job, correctly shaped. Every
            #       fixture in (a) and (b) fails the right-hand regex on its own
            #       merits, so deleting the `startswith(name, base)` guard
            #       entirely leaves them all untouched and the testset green.
            #       This one is the only fixture that pins the LEFT anchor --
            #       the guard that scopes the sweep to this run's output rather
            #       than everything in a shared SLURM-array output directory.
            #
            #       Its stem must be the SAME BYTE LENGTH as this run's, which
            #       is why it is derived rather than written out. Without the
            #       guard the predicate slices at `ncodeunits(base) + 1`
            #       regardless of what the name is, so a sibling of any OTHER
            #       length lands mid-name and fails the regex for the wrong
            #       reason -- sparing the file by accident and letting the
            #       mutant survive. Equal length forces the slice to land
            #       exactly on `.0000.tmp`, so only the left anchor can reject
            #       it. (A hand-written `othersample....` fixture was tried
            #       first and did NOT kill the mutant, for exactly this reason.)
            sibling_base = replace(basename(split_prefix), "sample." => "s4mple.")
            Test.@test ncodeunits(sibling_base) ==
                       ncodeunits(basename(split_prefix))
            Test.@test sibling_base != basename(split_prefix)
            sibling_chunk = joinpath(dir, sibling_base * ".0000.tmp")
            survivors = [
                outfile,
                outfile * ".bai",
                string(outfile, ".sort.tmp.0000.bam"),
                joinpath(dir, "unrelated.sorted.bam"),
                split_prefix,                                # bare prefix, no chunk index
                split_prefix * ".log",                        # non-numeric segment
                string(split_prefix, ".0000.tmp.bak"),        # trailing suffix past .tmp
                string(split_prefix, ".notanumber.tmp"),      # right shape, wrong segment
                sibling_chunk                                 # another job's LIVE chunk
            ]
            for s in survivors
                write(s, "keep")
            end

            result = Mycelia.cleanup_minimap_split_temps(split_prefix; verbose = false)

            Test.@test result.removed == length(chunks)
            Test.@test result.bytes == length(chunks) * 1024
            for c in chunks
                Test.@test !isfile(c)
            end
            for s in survivors
                Test.@test isfile(s)
            end
        end
    end

    Test.@testset "is idempotent and safe on a missing directory" begin
        # Idempotence needs a SECOND call after a real sweep. Both cases below
        # start with zero matching chunks, so on their own they hold for a
        # function that does nothing at all -- they test the empty case, not
        # the repeat case. The pre-run sweep added to every caller means a
        # second call on an already-swept directory is now the common path.
        mktempdir() do dir
            outfile = joinpath(dir, "idem.ref.mmi.minimap2.sorted.bam")
            prefix = Mycelia.minimap_split_prefix(outfile)
            for i in 0:1
                write(string(prefix, ".", lpad(i, 4, '0'), ".tmp"), rand(UInt8, 512))
            end
            keep = string(prefix, ".notanumber.tmp")
            write(keep, "keep")

            first_pass = Mycelia.cleanup_minimap_split_temps(prefix; verbose = false)
            Test.@test first_pass.removed == 2
            Test.@test first_pass.bytes == 2 * 512

            second_pass = Mycelia.cleanup_minimap_split_temps(prefix; verbose = false)
            Test.@test second_pass.removed == 0
            Test.@test second_pass.bytes == 0
            Test.@test isfile(keep)
        end
        mktempdir() do dir
            prefix = joinpath(dir, "nothing-here.bam.tmp")
            r = Mycelia.cleanup_minimap_split_temps(prefix; verbose = false)
            Test.@test r.removed == 0
            Test.@test r.bytes == 0
        end
        r = Mycelia.cleanup_minimap_split_temps(
            joinpath("/nonexistent-dir-for-mycelia-test", "x.bam.tmp"); verbose = false)
        Test.@test r.removed == 0
    end

    # ------------------------------------------------------------------
    # The PRE-RUN sweep at each call site.
    #
    # This is the half of the fix that actually covers the incident: a
    # `finally` never runs on SIGTERM or SIGKILL, so the reclaim has to happen
    # either in an `atexit` hook (SIGTERM only, inside SLURM's KillWait window)
    # or at the START of the next attempt. Until these tests existed, deleting
    # all three pre-run sweep calls left the entire suite green -- the primary
    # fix was reverted and nothing noticed.
    #
    # Structural rather than behavioural, and deliberately so: reaching the
    # sweep through any of the three entry points requires a command builder,
    # and every builder calls `add_bioconda_env` before it assembles a string.
    # A behavioural version lives below behind the external gate. These run on
    # the default CI path so that DELETING a sweep call cannot pass unnoticed,
    # which is the specific regression worth guarding.
    Test.@testset "pre-run sweep is wired at every call site" begin
        srcdir = joinpath(dirname(dirname(@__DIR__)), "src")

        function body_of(path, marker)
            text = read(joinpath(srcdir, path), String)
            idx = findfirst(marker, text)
            Test.@test idx !== nothing
            return text[first(idx):end]
        end

        # minimap_merge_map_and_split: swept unconditionally, BEFORE the
        # resume/caching branch. Order is the whole point -- after the branch
        # it would never run on the cached path, which is where a run killed
        # post-output strands its chunks.
        merge_body = body_of("alignments-and-mapping.jl", "    if run_mapping")
        sweep = findfirst("cleanup_minimap_split_temps(split_prefix)", merge_body)
        branch = findfirst("if nonempty_file(merged_bam)", merge_body)
        Test.@test sweep !== nothing
        Test.@test branch !== nothing
        Test.@test first(sweep) < first(branch)

        # merge_and_map_single_end_samples is the deliberate ASYMMETRY: its
        # sweep sits INSIDE the branch, not before it. `outbase` defaults to a
        # date-only string, so two same-day runs in one CWD share a split
        # prefix, and an entry-time sweep would delete a live peer's chunks --
        # trading a recoverable leak for unrecoverable truncated output. If
        # someone "fixes the inconsistency" by hoisting this call, that is a
        # correctness regression, and this assertion is what catches it.
        seq_body = body_of("sequence-comparison.jl",
            "    if !isfile(minimap_result.outfile)")
        seq_sweep = findfirst(
            "Mycelia.cleanup_minimap_split_temps(minimap_result.split_prefix)", seq_body)
        Test.@test seq_sweep !== nothing
        Test.@test first(seq_sweep) < first(findfirst("run(minimap_result.cmd)", seq_body))

        # prepare_binning_test_inputs: unconditional, like the merge path --
        # its inputs_dir comes from mktempdir(), so it cannot collide.
        util = read(joinpath(srcdir, "testing-utilities.jl"), String)
        u_sweep = findfirst(
            "Mycelia.cleanup_minimap_split_temps(mapping.split_prefix)", util)
        u_branch = findfirst("if !isfile(mapping.outfile)", util)
        Test.@test u_sweep !== nothing
        Test.@test u_branch !== nothing
        Test.@test first(u_sweep) < first(u_branch)
    end


    # Behavioural counterpart to the structural testset above. Gated for the
    # same reason the builder testset is: the entry point needs a command
    # builder, and every builder calls `add_bioconda_env` first.
    #
    # `run_mapping = true` is what reaches the sweep; a pre-existing nonempty
    # `merged_bam` then sends execution down the resume branch, which never
    # invokes minimap2. So this observes the real sweep, at the real call site,
    # without needing minimap2 to run.
    if get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
        Test.@testset "pre-run sweep reclaims on the resume path" begin
            mktempdir() do dir
                ref = joinpath(dir, "ref.fa")
                write(ref, ">ref\n" * "ACGT"^500 * "\n")
                fq = joinpath(dir, "reads.fq")
                write(fq, join(["@r$i\nACGT\n+\nIIII\n" for i in 1:100]))

                merged = joinpath(dir, "merged.sorted.bam")
                write(merged, "nonempty, so the resume branch is taken")

                prefix = Mycelia.minimap_split_prefix(merged)
                chunks = [string(prefix, ".", lpad(i, 4, '0'), ".tmp") for i in 0:1]
                for c in chunks
                    write(c, rand(UInt8, 2048))
                end
                # Must survive: right stem, wrong shape.
                survivor = string(prefix, ".notanumber.tmp")
                write(survivor, "keep")

                Mycelia.minimap_merge_map_and_split(
                    reference_fasta = ref,
                    mapping_type = "sr",
                    single_end_fastqs = [fq],
                    outdir = dir,
                    tmpdir = dir,
                    merged_bam = merged,
                    run_mapping = true,
                    run_splitting = false,
                    gzip_prefixed_fastqs = false
                )

                for c in chunks
                    Test.@test !isfile(c)
                end
                Test.@test isfile(survivor)
                # The cached BAM must be left exactly as found -- the sweep
                # reclaims orphans, it does not invalidate the resume.
                Test.@test isfile(merged)
                Test.@test read(merged, String) ==
                          "nonempty, so the resume branch is taken"
            end
        end
    end

    # Gated, because `minimap_map` calls `add_bioconda_env("minimap2")` and
    # `add_bioconda_env("samtools")` unconditionally BEFORE it builds any
    # command string -- so merely asking for a command needs a working conda.
    # The three testsets above genuinely need no external tools, which is what
    # lets this file run on the default (MYCELIA_RUN_EXTERNAL=false) CI path
    # where the cleanup predicate is the part that actually matters.
    if get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
        Test.@testset "builders return a split_prefix the caller can clean up" begin
            mktempdir() do dir
                ref = joinpath(dir, "ref.fa")
                write(ref, ">ref\n" * "ACGT"^100 * "\n")
                fq = joinpath(dir, "reads.fq")
                write(fq, "@r1\nACGT\n+\nIIII\n")
                outfile = joinpath(dir, "out.sorted.bam")

                res = Mycelia.minimap_map(;
                    fasta = ref, fastq = fq, mapping_type = "sr",
                    outfile = outfile, as_string = true)

                Test.@test haskey(res, :split_prefix)
                Test.@test res.split_prefix == Mycelia.minimap_split_prefix(res.outfile)
                # The prefix handed back must be the one minimap2 is actually told
                # to use -- otherwise cleanup would target the wrong files.
                Test.@test occursin("--split-prefix=$(res.split_prefix)", res.cmd)
            end
        end
    end
end
