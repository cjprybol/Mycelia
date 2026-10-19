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
# Four testsets run on the default CI path: they exercise the prefix
# derivation, the reclaim helper against synthetic files, its idempotence, and
# the structural wiring of the pre-run sweep at each call site. Two self-gate
# behind MYCELIA_RUN_EXTERNAL=true, because both reach a command builder and
# every builder calls `add_bioconda_env` before it assembles a string.
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

    Test.@testset "cleanup cannot throw, including on an unreadable parent" begin
        # The no-throw contract is what lets all three call sites invoke this
        # from a `finally` with no wrapper. It was false until now: `isdir` sat
        # ABOVE the try, and Julia's `stat` re-raises on every errno except
        # ENOENT/ENOTDIR/EINVAL, so an unreadable parent threw an IOError that
        # would REPLACE the mapping exception -- and, from the atexit hook,
        # would fail an otherwise successful job at shutdown.
        mktempdir() do root
            locked = joinpath(root, "locked")
            mkpath(locked)
            chmod(locked, 0o000)
            try
                # NB: calling `Base.isdir` on a path under `locked` would itself
                # throw here -- that IS the defect under test, so it must not
                # appear in the test's own setup.
                r = Mycelia.cleanup_minimap_split_temps(
                    joinpath(locked, "sub", "out.bam.tmp"); verbose = false)
                Test.@test r.removed == 0
                Test.@test r.bytes == 0
            catch err
                Test.@test false  # any throw at all violates the contract
            finally
                chmod(locked, 0o755)
            end
        end
    end

    Test.@testset "the delete contract is documented" begin
        # Julia binds a docstring to the NEXT expression and does NOT skip
        # comments. Twelve comment lines and three consts once sat between this
        # docstring and its function, so it was dropped with no warning and
        # `?cleanup_minimap_split_temps` printed "No documentation found" --
        # leaving the function that decides which files get DELETED with no
        # written contract at all. `checkdocs = :none` means the docs build
        # cannot catch it either.
        doc = string(Base.Docs.doc(Mycelia.cleanup_minimap_split_temps))
        Test.@test !occursin("No documentation found", doc)
        Test.@test occursin("does not throw", doc)
        Test.@test occursin("split_prefix", doc)

        # The ownership contract is what decides deletions now, so it has to be
        # written down too.
        owner_doc = string(Base.Docs.doc(Mycelia.minimap_split_owner_state))
        Test.@test !occursin("No documentation found", owner_doc)
        Test.@test occursin("absent", owner_doc)
    end

    # Behavioural replacement for the source-text assertions that used to live
    # here. Those greped `src/alignments-and-mapping.jl` for the literal
    # `owns_output_paths = isnothing(tmpdir) && isnothing(merged_bam)` -- which
    # tested the codebase's VOCABULARY, not its behaviour, and would have gone
    # on passing for as long as the string survived. It did survive, and the
    # gate it described was dead code: when it held, `tmpdir` had just been set
    # to a fresh `mktempdir()`, so the swept directory was empty by
    # construction; when it did not hold, the sweep was skipped. A test that
    # asserts its own literal cannot notice that.
    Test.@testset "ownership sidecar classifies who owns the chunks" begin
        mktempdir() do dir
            prefix = joinpath(dir, "x.sorted.bam.tmp")
            owner = Mycelia.minimap_split_owner_path(prefix)

            # No sidecar. Every path that runs minimap2 writes one first, so
            # its absence means no live owner. This is the load-bearing case:
            # it is what reclaims the chunks a SIGKILL stranded.
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :absent

            # Our own pid on our own host.
            Mycelia.write_minimap_split_owner(prefix)
            Test.@test isfile(owner)
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :live

            # A pid that existed and no longer does. Captured while the process
            # is alive because `Base.getpid` throws once it has been reaped.
            proc = run(`sleep 30`; wait = false)
            dead_pid = Base.getpid(proc)
            kill(proc)
            wait(proc)
            write(owner, "pid=$(dead_pid)\nhost=$(Base.gethostname())\n")
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :dead

            # Another host: our pid space says nothing about it, and under a
            # SLURM array sharing one tmpdir that is the COMMON case. Removing
            # this check is the mutant that deletes a live peer's chunks across
            # nodes, so it must not be reclaimable.
            write(owner, "pid=$(Base.Libc.getpid())\nhost=definitely-not-this-host\n")
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :unverifiable

            # Garbage is not evidence of absence.
            write(owner, "this is not a sidecar\n")
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :unverifiable

            # A pid with no host is equally unusable.
            write(owner, "pid=1\n")
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :unverifiable

            Mycelia.remove_minimap_split_owner(prefix)
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :absent
        end
    end

    Test.@testset "pre-run sweep reclaims orphans and spares a live peer" begin
        function seed(dir)
            prefix = joinpath(dir, "x.sorted.bam.tmp")
            chunks = [string(prefix, ".", lpad(i, 4, '0'), ".tmp") for i in 0:2]
            for c in chunks
                write(c, rand(UInt8, 256))
            end
            return prefix, chunks
        end

        # :absent -> reclaim. Mutating this to "skip" leaves the SIGKILL case
        # uncovered, which is the entire reason the sweep exists.
        mktempdir() do dir
            prefix, chunks = seed(dir)
            r = Mycelia.cleanup_minimap_split_temps(
                prefix; skip_if_owner_live = true, verbose = false)
            Test.@test r.skipped == false
            Test.@test r.removed == length(chunks)
            Test.@test all(c -> !isfile(c), chunks)
        end

        # :live -> skip. This is the mutant that unlinked a live peer's chunks;
        # minimap2 reopens every chunk by name at merge time, so the victim
        # aborts rather than silently degrading.
        mktempdir() do dir
            prefix, chunks = seed(dir)
            Mycelia.write_minimap_split_owner(prefix)   # us: alive by construction
            r = Mycelia.cleanup_minimap_split_temps(
                prefix; skip_if_owner_live = true, verbose = false)
            Test.@test r.skipped == true
            Test.@test r.removed == 0
            Test.@test all(isfile, chunks)
            # The sidecar must SURVIVE a skip. If a skip deleted it, the very
            # next sweep would read `:absent` and delete the chunks this one
            # just spared -- a two-pass version of the same bug.
            Test.@test isfile(Mycelia.minimap_split_owner_path(prefix))
        end

        # :dead -> reclaim, and the stale sidecar goes with it. `removed`
        # counting the chunks and NOT the sidecar also proves `.owner` is never
        # matched as a chunk.
        mktempdir() do dir
            prefix, chunks = seed(dir)
            proc = run(`sleep 30`; wait = false)
            dead_pid = Base.getpid(proc)
            kill(proc)
            wait(proc)
            write(Mycelia.minimap_split_owner_path(prefix),
                "pid=$(dead_pid)\nhost=$(Base.gethostname())\n")
            r = Mycelia.cleanup_minimap_split_temps(
                prefix; skip_if_owner_live = true, verbose = false)
            Test.@test r.skipped == false
            Test.@test r.removed == length(chunks)
            Test.@test !isfile(Mycelia.minimap_split_owner_path(prefix))
        end

        # Foreign host -> skip.
        mktempdir() do dir
            prefix, chunks = seed(dir)
            write(Mycelia.minimap_split_owner_path(prefix),
                "pid=$(Base.Libc.getpid())\nhost=some-other-node\n")
            r = Mycelia.cleanup_minimap_split_temps(
                prefix; skip_if_owner_live = true, verbose = false)
            Test.@test r.skipped == true
            Test.@test all(isfile, chunks)
        end

        # The DEFAULT (skip_if_owner_live = false) must ignore a live sidecar.
        # That is the `finally`/atexit path, where the caller IS the owner --
        # consulting its own record there would strand the chunks it is in the
        # middle of abandoning.
        mktempdir() do dir
            prefix, chunks = seed(dir)
            Mycelia.write_minimap_split_owner(prefix)
            r = Mycelia.cleanup_minimap_split_temps(prefix; verbose = false)
            Test.@test r.skipped == false
            Test.@test r.removed == length(chunks)
        end
    end

    Test.@testset "track writes the sidecar, untrack removes it" begin
        # The sweep's `:absent -> reclaim` rule is only sound because every
        # site that runs minimap2 tracks first. If tracking stopped writing the
        # sidecar, a concurrent peer would read as `:absent` and be reclaimed
        # mid-flight -- so this pairing is the mechanism's keystone.
        mktempdir() do dir
            prefix = joinpath(dir, "y.sorted.bam.tmp")
            Mycelia.track_minimap_split_prefix(prefix)
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :live
            Mycelia.untrack_minimap_split_prefix(prefix)
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :absent
        end
    end

    Test.@testset "write refuses to steal a live peer's sidecar" begin
        # Without this guard, a NEW run's track_minimap_split_prefix would
        # silently overwrite a live peer's sidecar with its own pid/host. The
        # new caller's later `finally`-block cleanup would then delete both
        # the peer's chunks AND the peer's (now-stolen) sidecar, and a third
        # reader would see `:absent` and reclaim chunks that are still in
        # flight -- defeating `skip_if_owner_live` entirely rather than merely
        # losing disk. This is the mutant that removing the guard reintroduces.
        mktempdir() do dir
            prefix = joinpath(dir, "z.sorted.bam.tmp")
            owner = Mycelia.minimap_split_owner_path(prefix)

            # A foreign PID, alive by construction (our own test process' pid,
            # but recorded under a different host so the state resolves to
            # `:unverifiable` -- still not something a write may steal).
            write(owner, "pid=$(Base.Libc.getpid())\nhost=some-other-node\n")
            Test.@test_throws Exception Mycelia.write_minimap_split_owner(prefix)
            # The foreign sidecar must survive the refused write untouched.
            Test.@test read(owner, String) ==
                       "pid=$(Base.Libc.getpid())\nhost=some-other-node\n"

            # A genuinely live, same-host, different-pid owner (`:live`).
            proc = run(`sleep 30`; wait = false)
            live_pid = Base.getpid(proc)
            try
                write(owner, "pid=$(live_pid)\nhost=$(Base.gethostname())\n")
                Test.@test Mycelia.minimap_split_owner_state(prefix) == :live
                Test.@test_throws Exception Mycelia.write_minimap_split_owner(prefix)
                Test.@test Mycelia.minimap_split_owner_state(prefix) == :live
            finally
                kill(proc)
                wait(proc)
            end

            # A dead owner is NOT a collision -- overwrite must succeed and
            # hand ownership to us.
            write(owner, "pid=$(live_pid)\nhost=$(Base.gethostname())\n")
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :dead
            Mycelia.write_minimap_split_owner(prefix)
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :live

            # Re-tracking the SAME process/host (e.g. calling track twice for
            # one prefix without an intervening untrack) is not a collision.
            Mycelia.write_minimap_split_owner(prefix)
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :live

            # No sidecar at all (`:absent`) is the ordinary, non-colliding case.
            Mycelia.remove_minimap_split_owner(prefix)
            Mycelia.write_minimap_split_owner(prefix)
            Test.@test Mycelia.minimap_split_owner_state(prefix) == :live
        end
    end

    Test.@testset "pre-run sweep is wired at every call site" begin
        srcdir = joinpath(dirname(dirname(@__DIR__)), "src")

        # Line-anchored AND asserted unique. A bare substring marker is not
        # enough: "    if run_mapping" also occurs inside
        # "            if run_mapping && !isfile(index_file)" 180 lines earlier,
        # so `findfirst` sliced from the wrong block and a sweep hoisted ABOVE
        # `if run_mapping` still satisfied `sweep < branch`. That mutant --
        # which would make the sweep fire on the command-generation-only path
        # and unlink chunks belonging to a process running that same command --
        # survived the whole suite until this was fixed.
        # Plain substring counting, NOT a regex: `escape_string` escapes for a
        # Julia string literal, so `Regex(escape_string("\n    if x\n"))` looks
        # for a literal backslash-n and matches nothing. That mistake made this
        # very assertion vacuous on its first run.
        function count_occurrences(needle, haystack)
            n = 0
            i = firstindex(haystack)
            while true
                r = findnext(needle, haystack, i)
                r === nothing && break
                n += 1
                i = first(r) + 1
            end
            return n
        end

        function body_of(path, marker)
            text = read(joinpath(srcdir, path), String)
            Test.@test count_occurrences(marker, text) == 1
            idx = findfirst(marker, text)
            if idx === nothing
                Test.@test false  # clean failure, not a MethodError on first(nothing)
                return ""
            end
            return text[first(idx):end]
        end

        # minimap_merge_map_and_split: swept BEFORE the resume/caching branch.
        # Order is the whole point -- after the branch it would never run on
        # the cached path, which is where a run killed post-output strands its
        # chunks. The ownership check is what makes running it here safe, so
        # both facts are asserted together.
        merge_body = body_of("alignments-and-mapping.jl", "\n    if run_mapping\n")
        sweep = findfirst(
            "cleanup_minimap_split_temps(split_prefix; skip_if_owner_live = true)",
            merge_body)
        branch = findfirst("if nonempty_file(merged_bam)", merge_body)
        Test.@test sweep !== nothing
        Test.@test branch !== nothing
        Test.@test first(sweep) < first(branch)

        # merge_and_map_single_end_samples keeps its sweep INSIDE the branch.
        # `outbase` defaults to a date-only string, so two same-day runs in one
        # CWD share a split prefix; the ownership check now makes that safe,
        # but the position is left alone because hoisting it is a behaviour
        # change this PR does not need. The remaining hole (a run killed AFTER
        # writing `outfile` keeps its chunks, since the branch is then skipped
        # forever) is a leak, not a deletion, and is now CLOSEABLE by hoisting
        # -- tracked separately rather than done here.
        seq_body = body_of("sequence-comparison.jl",
            "    if !isfile(minimap_result.outfile)")
        seq_sweep = findfirst("skip_if_owner_live = true", seq_body)
        Test.@test seq_sweep !== nothing
        Test.@test first(seq_sweep) < first(findfirst("run(minimap_result.cmd)", seq_body))

        # prepare_binning_test_inputs: swept on BOTH paths now. Its previous
        # `outdir === nothing` gate had the same defect as the merge entry
        # point's -- when it held, `inputs_dir` was a fresh `mktempdir()`, so
        # the swept directory was empty by construction. Non-collision is still
        # the safety condition; it is now established from the sidecar rather
        # than inferred from which arguments the caller passed.
        util = read(joinpath(srcdir, "testing-utilities.jl"), String)
        u_sweep = findfirst(
            "Mycelia.cleanup_minimap_split_temps(mapping.split_prefix; skip_if_owner_live = true)",
            util)
        u_branch = findfirst("if !isfile(mapping.outfile)", util)
        Test.@test u_sweep !== nothing
        Test.@test u_branch !== nothing
        Test.@test first(u_sweep) < first(u_branch)
        # The dead gate must be GONE, not merely bypassed. Leaving it in place
        # would restore the condition under which the sweep can never fire.
        Test.@test !occursin(
            "if outdir === nothing\n        Mycelia.cleanup_minimap_split_temps",
            util)
    end

    Test.@testset "every site that maps also tracks first" begin
        # The keystone invariant. `:absent -> reclaim` is sound ONLY because a
        # live run always has a sidecar on disk, so a call site that invokes
        # minimap2 without calling `track_minimap_split_prefix` first would
        # leave its chunks classified as orphans and reclaimable mid-flight by
        # any concurrent peer. Structural because the alternative needs conda.
        srcdir = joinpath(dirname(dirname(@__DIR__)), "src")
        for (path, runcall) in (
            ("alignments-and-mapping.jl", "run(minimap_cmd)"),
            ("sequence-comparison.jl", "run(minimap_result.cmd)"),
            ("testing-utilities.jl", "run(mapping.cmd)")
        )
            text = read(joinpath(srcdir, path), String)
            r = findfirst(runcall, text)
            Test.@test r !== nothing
            r === nothing && continue
            before = text[1:first(r)]
            track = findlast("track_minimap_split_prefix(", before)
            Test.@test track !== nothing
            # ...and it must be the TRACK call, not the UNtrack one that
            # follows in the `finally`. `findlast` on the prefix would happily
            # match `untrack_minimap_split_prefix(` otherwise.
            if track !== nothing
                # `prevind(before, first(track), 2)`, not `first(track) - 2`:
                # the source comments around these call sites use em-dashes
                # (multi-byte UTF-8), so raw byte arithmetic can land mid-
                # character and throw a StringIndexError. `prevind` steps back
                # by codepoints and is clamped to `firstindex` when the match
                # is near the start of `before`.
                window_start = max(firstindex(before), prevind(before, first(track), 2))
                Test.@test !occursin("untrack_minimap_split_prefix(",
                    before[window_start:end])
            end
        end
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
    # The four ungated testsets above need no external tools, which is what
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
