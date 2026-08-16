# Gate DECISION for the qualmer profile quality comparison (td-n8ax).
#
# Deliberately dependency-free (no `import Mycelia`) and side-effect-free so the
# gate's failure path can be exercised on synthetic results without running a
# real assembly. `qualmer_profile_quality_compare.jl` includes this file and
# `error()`s on a failing verdict, which is what makes the process exit non-zero.
#
# WHY THIS IS A SEPARATE FILE
# ---------------------------
# The decision logic used to be inline in the benchmark script, where it only
# `println`ed — so the "authoritative" gate exited 0 on an assembly error and on
# a DEGRADED verdict alike. A gate that cannot go red reads as coverage and stops
# anyone from looking. Reading the fixed code is not evidence that it goes red;
# the failure path has to be RUN. Splitting the decision out is what makes that
# possible without a multi-hour assembly, and it is covered by
# `test/4_assembly/qualmer_profile_quality_verdict_test.jl`.

"""
    qualmer_qc_gate_verdict(results, requested; tol = 0.01)

Pure gate decision for the aggregate-vs-`:full` qualmer quality comparison.

# Arguments
- `results`: vector of `(profile, n_contigs, genome_fraction, longest, elapsed)`
  tuples, one per profile that assembled SUCCESSFULLY. Profiles that errored are
  simply absent (that absence is itself a failure — see below).
- `requested`: every profile that was asked for, so a silently-missing result is
  detectable.
- `tol`: how far below `:full`'s genome fraction an aggregate profile may sit and
  still count as PARITY.

# Returns
`(verdicts, failures)` — `verdicts` are human-readable `VERDICT …` lines for the
transcript; `failures` is empty iff the gate PASSES.

Three independent ways to fail, all of which previously exited 0:
1. a requested profile produced no result (assembly errored or was skipped),
2. the `:full` baseline is absent, so parity cannot be established at all,
3. an aggregate profile's genome fraction is below `full_gf - tol` (DEGRADED).
"""
function qualmer_qc_gate_verdict(results, requested; tol::Float64 = 0.01)
    verdicts = String[]
    failures = String[]

    completed = Set(r[1] for r in results)
    for profile in requested
        profile in completed || push!(failures,
            "profile $profile produced NO result (assembly errored or was skipped)")
    end

    if !(:full in completed)
        push!(failures,
            "baseline profile :full is absent from results — parity cannot be established")
        return verdicts, failures
    end

    full_gf = first(r[3] for r in results if r[1] == :full)
    for (profile, _, g, _, _) in results
        profile == :full && continue
        degraded = g < full_gf - tol
        push!(verdicts,
            "VERDICT $profile vs :full genome_fraction: $g vs $full_gf -> " *
            (degraded ? "DEGRADED" : "PARITY"))
        degraded && push!(failures,
            "profile $profile DEGRADED vs :full: genome_fraction $g < $(full_gf - tol)")
    end

    return verdicts, failures
end

"""
    qualmer_qc_gate_check(results, requested; tol = 0.01, io = stdout) -> Bool

Print the verdict lines and a single `GATE PASS` / `GATE FAIL …` summary, then
report whether the gate passed. Callers must turn `false` into a non-zero process
exit — printing alone is what made this gate unfalsifiable in the first place.
"""
function qualmer_qc_gate_check(results, requested; tol::Float64 = 0.01, io::IO = stdout)
    verdicts, failures = qualmer_qc_gate_verdict(results, requested; tol = tol)
    for v in verdicts
        println(io, v)
    end
    if isempty(failures)
        println(io, "GATE PASS: every requested profile completed and is at parity with :full")
        return true
    end
    for f in failures
        println(io, "GATE FAIL: $f")
    end
    return false
end
