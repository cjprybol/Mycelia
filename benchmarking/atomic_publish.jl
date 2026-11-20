# Unique-temp atomic publish for Track A artifacts.
#
# A FIXED temp path (`<basename>.tmp`, or `target * ".tmp"`) is safe for a single
# writer and unsafe wherever several processes share one results tree. The harness
# header advertises that an HPC array job can split the matrix and share one tree,
# and `write_aggregate` runs once per cell, so two writers otherwise open the SAME
# temp inode and both rename it. The rename then publishes interleaved content as
# a complete table, which is worse than a short table: a short one is detectable,
# a corrupt one is not.
#
# `mktemp(dir)` creates the file (O_EXCL) in the same directory as the target, so
# two writers cannot be handed the same path before either creates it. A
# name-only temporary path does not: it returns a path and leaves creation to
# the caller, and two callers can observe the same free name. The rename stays
# atomic because the temp is on the same filesystem as the target.
#
# This used to live only as `_publish_atomically` in `track_a_baseline_benchmark.jl`.
# That file imports Mycelia and reads `ARGS` at load, so `track_a_merge_hosts.jl`
# cannot include it. Both call this copy.

"""
    _publish_atomically(write!, target) -> String

Write via `write!(tmp_path)` to a unique temp file in `dirname(target)`, then
rename it onto `target`.

The temp file is created with `mktemp` before `write!` runs, and that handle is
closed so `write!` can open the path itself. An exception removes the temp and
leaves `target` unchanged. `InterruptException` is rethrown after that cleanup;
it is not turned into a successful publish.
"""
function _publish_atomically(write!, target::AbstractString)
    dir = dirname(target)
    mkpath(dir)
    tmp, io = mktemp(dir)
    try
        close(io)
        write!(tmp)
        mv(tmp, target; force = true)
    catch e
        isopen(io) && close(io)
        isfile(tmp) && rm(tmp; force = true)
        # Ctrl-C must propagate. The cleanup above is the only work this catch
        # does; every exception, interrupt included, is rethrown.
        e isa InterruptException && rethrow()
        rethrow()
    end
    return target
end

"""
    _atomic_write(f, path) -> String

IO form of [`_publish_atomically`](@ref). `f(io)` writes a unique temp file in
the same directory, then the rename publishes it. A reader never observes a
partially-written file, and an interrupted run leaves the previous contents.
"""
function _atomic_write(f, path::AbstractString)
    _publish_atomically(path) do tmp
        open(f, tmp, "w")
    end
    return path
end
