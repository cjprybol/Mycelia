# Per-cell peak-RSS sampler shared by the Track A harness and the D5 grid
# driver (td-yddpr.5). Moved verbatim from track_a_baseline_benchmark.jl, whose
# top-level ARGS parsing makes it unsafe to include() from another CLI.
# No ARGS, no side effects. Requires nothing beyond Base.

# === Per-cell memory measurement ===

# `Sys.maxrss()` is the process-lifetime high-water mark (ru_maxrss) and is
# monotonically non-decreasing, so the old `max(0, Sys.maxrss() - rss0)` could only
# be nonzero when a cell beat the all-time PROCESS peak. It measured "excess over
# historical peak", not per-cell peak.
#
# Measured on the 2026-07-25 Lovelace run (288 cells): 271/288 rows (94%) were
# exactly 0 — 144/144 kmer AND 127/144 qualmer. The effect is a global ratchet, not
# a per-pair alternation: all 17 nonzero rows are qualmer only because qualmer is
# the inner loop's first arm, they sit in the first 48% of the run, and the ratchet
# closes for good at T4/ont/100x/seed42 = 51.3 GiB — after which nothing beats it.
# The 6% that survived understate true usage by orders of magnitude (e.g.
# T4/illumina/100x/seed456 = 11 MiB for a multi-GiB cell), because each is a
# delta over the previous all-time peak rather than a measurement. LRC's 144-cell
# run agrees: 131/144 zero (72/72 kmer, 59/72 qualmer).
#
# Read VmRSS from /proc/self/status instead: it reports explicit kB, so unlike
# /proc/self/statm it needs no Sys.PAGESIZE (which does not exist in Julia 1.10).
function current_rss_bytes()
    status_path = "/proc/self/status"
    isfile(status_path) || return nothing
    # `open(...) do` rather than bare `eachline(path)`: eachline closes its handle at
    # end-of-iteration, and VmRSS sits near the top of the file so the early return
    # below ALWAYS exits first. Measured: 2000 polls leaked 2000 descriptors, released
    # only by GC. At 20 Hz that is ~72k/hour, and EMFILE would surface inside the
    # sampler task where it is hardest to diagnose.
    return open(status_path) do io
        for line in eachline(io)
            startswith(line, "VmRSS:") || continue
            fields = split(line)
            length(fields) >= 2 || return nothing
            kilobytes = tryparse(Int, fields[2])
            return kilobytes === nothing ? nothing : kilobytes * 1024
        end
        return nothing
    end
end

"""
Run `f`, returning `(; value, wall_seconds, peak_rss_bytes, rss_baseline_bytes, method)`.

`peak_rss_bytes` is the highest WHOLE-PROCESS resident set observed while `f` ran. It is
deliberately not called a per-cell figure: it includes the Julia runtime, the loaded
package images, every other thread, and memory retained (not returned to the OS) by
earlier cells. `rss_baseline_bytes` is the process RSS immediately before `f`, so a caller
wanting a cell-attributable increment can take the difference — but the absolute level is
what an HPC memory request actually needs, so that is what the primary column holds.

The sampler runs via `Threads.@spawn`, which requires a second thread: the assembly is
CPU-bound Julia that never yields, so a same-thread `@async` sampler would never be
scheduled. Both committed sbatch wrappers export JULIA_NUM_THREADS=\$SLURM_CPUS_PER_TASK —
run_track_a_baseline_lrc.sbatch (16) and run_track_a_baseline_nersc.sbatch (4) — so
sampling is the production path under SLURM. Julia defaults to ONE thread, so a bare
`julia track_a_baseline_benchmark.jl` reverts to the broken high-water semantics; launch
with `julia -t N` or the column is meaningless. That revert is no longer silent — a
one-time `@warn` fires from the fallback branch below whenever it happens.

A second thread is necessary but not sufficient: the sampler's inter-poll wait must be
event-loop-free AND must yield. It uses `Libc.systemsleep` followed by `yield()` —
`sleep` parks on a libuv timer that a non-yielding CPU-bound `f` prevents from ever being
serviced, while an unyielding wait can starve the thread-1-sticky root task instead. See
the comment at the call site for the measured poll counts and hang rates.

`method` records how the number was obtained, so a reader never has to infer it and rows
measured different ways are never averaged together:
  "sampled"          - process-RSS high-water, polled while the cell ran, with at least
                       TWO polls and EVERY poll succeeding. A flat cell whose RSS never
                       exceeds its baseline is still "sampled" — that is a real
                       observation.
  "sampled-partial"  - at least two polls succeeded and at least one FAILED (threw, or
                       returned `nothing`). Sampling continued, but the peak may have been
                       missed while the probe was failing, so the number is a lower bound,
                       not a measurement.
  "sampled-degraded" - the sampler completed at most ONE successful read (never scheduled,
                       starved, or nearly every read failed). Below two samples there is
                       no high-water mark to speak of: the number is the entry baseline or
                       a single arbitrary instant, NOT a measurement.
  "highwater-delta"  - single-threaded or no /proc: the old, known-broken semantics
  "unknown"          - checkpoint predates this column (see OPTIONAL_KEY_DEFAULTS)
Values under different methods are different quantities. Always filter on
`peak_rss_method` before aggregating.

`probe` is the RSS reader, injectable ONLY so the failure labelling above can be tested:
`current_rss_bytes` returns `nothing` on macOS and there is no portable way to make a
real /proc read fail on the third poll.
"""
function timed_with_peak_rss(f; interval_seconds::Float64 = 0.05,
        probe = current_rss_bytes)
    # Guard the baseline read too: everything here is careful never to let telemetry fail
    # science, and an unguarded /proc read error would abort the cell before f() even ran.
    baseline = try
        probe()
    catch e
        @warn "baseline RSS read failed; falling back to high-water semantics" exception = e
        nothing
    end

    can_sample = baseline !== nothing && Threads.nthreads() > 1
    if !can_sample
        # The docstring says a bare `julia` invocation "silently reverts to the broken
        # high-water semantics". `can_sample == false` knows that at runtime, so say it
        # instead of leaving the operator to infer it from a column of near-zeros after a
        # multi-day run. `maxlog` keeps it to one line, not one per cell.
        if Threads.nthreads() == 1
            @warn "peak_rss_bytes is falling back to the known-broken high-water " *
                  "semantics: Julia has ONE thread, so the RSS sampler cannot be " *
                  "scheduled against a CPU-bound assembly. Relaunch with `julia -t N` " *
                  "(both sbatch wrappers export JULIA_NUM_THREADS) or treat the column " *
                  "as unusable." maxlog = 1
        elseif baseline === nothing
            @warn "peak_rss_bytes is falling back to the known-broken high-water " *
                  "semantics: /proc/self/status is unreadable on this host." maxlog = 1
        end
        rss0 = Sys.maxrss()
        timed = @timed f()
        return (value = timed.value, wall_seconds = timed.time,
            peak_rss_bytes = max(0, Sys.maxrss() - rss0),
            rss_baseline_bytes = baseline === nothing ? -1 : baseline,
            method = "highwater-delta")
    end

    peak = Threads.Atomic{Int}(baseline)
    keep_sampling = Threads.Atomic{Bool}(true)
    # Count SUCCESSFUL READS, not increases over the baseline. The failure worth flagging
    # is a sampler that never observed anything — starved, or every /proc read failed —
    # because it then returns the baseline as though it were a measurement: a large,
    # stable, plausible byte count, strictly harder to spot than the exactly-0 signature
    # of the bug being replaced. A cell whose RSS never rises above its baseline is NOT
    # degraded; that happens whenever the allocation fits inside already-resident heap.
    reads = Threads.Atomic{Int}(0)
    # Count FAILED reads separately, because "never observed anything" and "stopped
    # observing partway" are different lies. With the try OUTSIDE the loop, the first
    # throwing poll ended sampling for the rest of the cell, yet `reads > 0` still
    # labelled the frozen partial maximum "sampled" — a plausible number presented as a
    # measurement, which is the exact failure mode the read count was added to catch.
    failed_reads = Threads.Atomic{Int}(0)
    sampler = Threads.@spawn begin
        # The sampler must never be able to fail a good cell: a transient /proc read error
        # is telemetry, not science. The try is INSIDE the loop, so one bad poll costs one
        # sample rather than the rest of the cell, and the counts below report it.
        try
            while keep_sampling[]
                observed = try
                    probe()
                catch
                    nothing
                end
                if observed === nothing
                    # A `nothing` is a FAILED read, not a quiet no-op. current_rss_bytes
                    # never throws — it signals failure by RETURNING nothing on all four
                    # of its paths (no /proc file, no VmRSS line, a short line, an
                    # unparseable field) — so incrementing only inside the `catch` left a
                    # sampler that had gone blind mid-cell reporting failed_reads == 0 and
                    # being labelled "sampled". Counting here and NOT in the catch keeps a
                    # thrown read and a nothing read worth exactly one failure each.
                    Threads.atomic_add!(failed_reads, 1)
                else
                    Threads.atomic_add!(reads, 1)
                    observed > peak[] && Threads.atomic_max!(peak, observed)
                end
                # Libc.systemsleep, NOT sleep. `sleep` parks the task on a libuv timer, and
                # that timer only fires when the event loop is serviced — which does not
                # happen while the main task is inside a non-yielding CPU-bound call, which
                # is exactly what Rhizomorph.assemble_genome is. A second thread is
                # necessary but NOT sufficient; the wait has to be event-loop-free too.
                #
                # Measured on Linux, julia 1.10.11, `-t 4`, over a 3.0 s non-yielding
                # workload at interval 0.05: sleep completed 1 poll of the ~60 due,
                # systemsleep 55-59. On a workload whose RSS ramps ~1.5 GiB and is then
                # released, that is the difference between reporting peak - baseline =
                # 0.3 MiB and 1465.1 MiB. The single poll lands at an arbitrary instant —
                # near the baseline in practice, since the cell has barely started — so the
                # column understates true usage by orders of magnitude while carrying the
                # "sampled" label. That is the same class of defect as the Sys.maxrss()
                # ratchet this sampler replaced, and just as hard to spot: a plausible
                # byte count, not an obvious zero.
                #
                # The trailing `yield()` is NOT decoration — without it this hangs. Julia's
                # ROOT task is sticky to thread 1, and the scheduler is cooperative, so a
                # spawned task that occupies thread 1 and never yields can make the root
                # task permanently unresumable. systemsleep alone does exactly that: it
                # blocks the thread and the surrounding loop has no yield point, so once
                # the sampler lands on thread 1 the root task starves. The old `sleep` hid
                # this because yielding is essentially all it did. Measured over this
                # harness's own test sequence, which parks the root task on `sleep()` while
                # a sampler runs: systemsleep alone completed 0/4 runs and a GC-safepointed
                # busy-wait 2/4, both hanging until killed, while systemsleep + yield
                # completed 30/30 across -t 2, -t 4 and -t 16 and still polled 40/40.
                # `yield()` re-queues the sampler; a CPU-bound root task holds thread 1, so
                # an idle thread picks the sampler back up within microseconds and the
                # interval is unaffected.
                Libc.systemsleep(interval_seconds)
                yield()
            end
        catch
            # Last-resort net for a failure OUTSIDE the probe (an interrupted wait, say).
            # Deliberately swallowed; the two read counts report the degradation.
        end
    end

    timed = try
        @timed f()
    finally
        # Stop AND reap the sampler here, not after the block. With `wait` outside, a
        # sampler exception on the success path would raise TaskFailedException and
        # DISCARD a completed assembly — telemetry destroying science — while on the throw
        # path `wait` would never be reached at all and the sampler would outlive the cell.
        keep_sampling[] = false
        try
            wait(sampler)
        catch
            # sampler already reported via `reads`; never let it mask f()'s own outcome
        end
    end
    # Read the counters only after the `finally` above has reaped the sampler, so the
    # label describes the whole cell rather than a racing snapshot of it.
    method = if reads[] <= 1
        # ONE successful poll is not a peak. A high-water mark needs something to take a
        # maximum OVER; with a single sample the returned value is one arbitrary instant
        # of the cell (and with zero, the entry baseline), yet in the output column it is
        # a large, stable, plausible byte count indistinguishable from a real measurement
        # — which is exactly what the "sampled" label would promise. `== 0` alone let the
        # under-polling failure through: a sampler that manages exactly one poll per cell
        # returns ~the baseline for every cell and was labelled clean.
        "sampled-degraded"
    elseif failed_reads[] > 0
        "sampled-partial"
    else
        "sampled"
    end
    return (value = timed.value, wall_seconds = timed.time,
        peak_rss_bytes = peak[], rss_baseline_bytes = baseline, method = method)
end
