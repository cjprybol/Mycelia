# Pure VERDICT gate for the correction-accuracy benchmark (no Mycelia, no I/O).
#
# The harness may print VERDICT only when three gates all pass:
#   (1) scale floor (computed by the caller via rhizomorph_scale_guard.jl)
#   (2) real-arm aggregate scored fraction
#   (3) Control B null-arm health on the cells that carry injected errors
#
# Gate (3) is what td-womt closes. `null_ok` and `null_scored_fraction` used to
# be computed per cell and then ignored, so a null arm that crashed, mis-joined
# read ids, or mass-dropped rows could still sit under a VERDICT header while
# Control B had silently done nothing. A failure of (3) downgrades the run to
# SMOKE-ONLY. Cells with no injected errors (Control A, err=0) do not enter
# gate (3): the null has no errors to fail to recover there.
#
# This file is dependency-free so fixture tests can call the same function the
# harness uses, without loading Mycelia or launching a benchmark.

"""
    null_arm_misjoined(null_reads_scored, null_corrected_unjoined) -> Bool

Return `true` when the null corrector emitted reads and none of them joined to
truth. That is a read-ID mismatch, not a genuine recall collapse (`null_ok`
stays true because the corrector did not throw).
"""
function null_arm_misjoined(null_reads_scored::Integer, null_corrected_unjoined::Integer)
    return null_reads_scored == 0 && null_corrected_unjoined > 0
end

"""
One sweep cell's Control B bookkeeping, as the VERDICT gate consumes it.

`injected_errors` is the number of substitutions the simulator injected (0 for
Control A). `null_ok` is false when that cell's null corrector threw.
`null_reads_scored` / `null_corrected_unjoined` are the null arm's join counts
from `per_base_metrics`.
"""
Base.@kwdef struct NullArmCell
    injected_errors::Int
    n_reads::Int
    null_ok::Bool
    null_reads_scored::Int
    null_corrected_unjoined::Int
end

"""
    null_arm_health(cells) -> NamedTuple

Pool Control B health over cells with `injected_errors > 0`.

`null_ok` is true only when every such cell's null corrector ran and joined at
least one read whenever it emitted any (`null_arm_misjoined`).
`null_scored_fraction` is pooled `null_reads_scored / n_reads` over those cells
(a mass drop drives it down). Both are vacuous when no cell carries injected
errors: `null_ok=true` and `null_scored_fraction=NaN`.
"""
function null_arm_health(cells::AbstractVector{NullArmCell})
    n_reads = 0
    n_scored = 0
    n_injected = 0
    null_ok = true
    for cell in cells
        cell.injected_errors > 0 || continue
        n_injected += 1
        n_reads += cell.n_reads
        n_scored += cell.null_reads_scored
        if !cell.null_ok || null_arm_misjoined(cell.null_reads_scored, cell.null_corrected_unjoined)
            null_ok = false
        end
    end
    frac = n_injected == 0 ? NaN : (n_reads == 0 ? 0.0 : n_scored / n_reads)
    return (null_ok = null_ok, null_scored_fraction = Float64(frac), n_injected = n_injected)
end

"""
    null_scored_fraction_ok(null_scored_fraction, min_scored_fraction) -> Bool

`NaN` means no injected-error cell existed, so the null scored floor is
vacuously met. A non-finite or below-floor fraction is not.
"""
function null_scored_fraction_ok(null_scored_fraction::Real, min_scored_fraction::Real)
    frac = Float64(null_scored_fraction)
    isnan(frac) && return true
    return isfinite(frac) && frac >= Float64(min_scored_fraction)
end

"""
    correction_accuracy_verdict_allowed(scale_ok, scored_ok, null_ok,
        null_scored_fraction; min_scored_fraction=0.5) -> Bool

Single VERDICT decision. All of `scale_ok`, `scored_ok`, `null_ok`, and the
null scored floor must hold. `null_ok` / `null_scored_fraction` describe
injected-error cells only (see `null_arm_health`).
"""
function correction_accuracy_verdict_allowed(scale_ok::Bool, scored_ok::Bool, null_ok::Bool,
        null_scored_fraction::Real; min_scored_fraction::Real = 0.5)
    return scale_ok && scored_ok && null_ok &&
           null_scored_fraction_ok(null_scored_fraction, min_scored_fraction)
end

function _format_gate_fraction(x::Real)
    x = Float64(x)
    isfinite(x) || return "NA"
    return string(round(x; digits = 3))
end

"""
    correction_accuracy_summary_line(; verdict_allowed, scale_ok, scored_ok,
        null_ok, metric, scale_floor, agg_scored_fraction, null_scored_fraction,
        min_scored_fraction) -> String

SMOKE-ONLY / VERDICT header. Both labels carry `null_scored_fraction` so a
Control B failure is visible on the same line as the decision.
"""
function correction_accuracy_summary_line(; verdict_allowed::Bool, scale_ok::Bool, scored_ok::Bool,
        null_ok::Bool, metric::Real, scale_floor::Real, agg_scored_fraction::Real,
        null_scored_fraction::Real, min_scored_fraction::Real)
    null_fmt = _format_gate_fraction(null_scored_fraction)
    null_bit = "null_scored_fraction=$(null_fmt)"
    scored_fmt = _format_gate_fraction(agg_scored_fraction)
    scale_txt = string(round(Float64(metric); digits = 0))
    if verdict_allowed
        null_clause = isfinite(Float64(null_scored_fraction)) ?
                      "$(null_bit) >= $(min_scored_fraction)" : null_bit
        return "=== VERDICT (scale $(scale_txt) bases >= $(scale_floor); " *
               "scored $(scored_fmt) >= $(min_scored_fraction); $(null_clause)) ==="
    end
    reasons = String[]
    scale_ok || push!(reasons, "scale metric $(scale_txt) < floor $(scale_floor)")
    scored_ok || push!(reasons, "scored fraction $(scored_fmt) < floor $(min_scored_fraction)")
    null_ok || push!(reasons, "null arm unhealthy on injected-error cells")
    if !null_scored_fraction_ok(null_scored_fraction, min_scored_fraction)
        push!(reasons, "null scored fraction $(null_fmt) < floor $(min_scored_fraction)")
    end
    isempty(reasons) && push!(reasons, "guard failed")
    return "=== SMOKE-ONLY ($(join(reasons, "; ")); $(null_bit)) ==="
end

"""
    correction_accuracy_verdict(; scale_ok, agg_scored_fraction, cells,
        min_scored_fraction=0.5, metric=0, scale_floor=0) -> NamedTuple

Apply the real-arm scored floor and the null-arm health of `cells`, then
combine them with `scale_ok`. Returns the mode (`"VERDICT"` or `"SMOKE-ONLY"`),
the summary line, and the gate bits the harness prints and records.
"""
function correction_accuracy_verdict(; scale_ok::Bool, agg_scored_fraction::Real,
        cells::AbstractVector{NullArmCell}, min_scored_fraction::Real = 0.5,
        metric::Real = 0.0, scale_floor::Real = 0.0)
    scored_ok = agg_scored_fraction >= min_scored_fraction
    health = null_arm_health(cells)
    null_scored_ok = null_scored_fraction_ok(health.null_scored_fraction, min_scored_fraction)
    verdict_allowed = correction_accuracy_verdict_allowed(
        scale_ok, scored_ok, health.null_ok, health.null_scored_fraction;
        min_scored_fraction = min_scored_fraction)
    mode = verdict_allowed ? "VERDICT" : "SMOKE-ONLY"
    line = correction_accuracy_summary_line(;
        verdict_allowed = verdict_allowed,
        scale_ok = scale_ok,
        scored_ok = scored_ok,
        null_ok = health.null_ok,
        metric = metric,
        scale_floor = scale_floor,
        agg_scored_fraction = agg_scored_fraction,
        null_scored_fraction = health.null_scored_fraction,
        min_scored_fraction = min_scored_fraction)
    return (;
        verdict_allowed = verdict_allowed,
        mode = mode,
        line = line,
        scale_ok = scale_ok,
        scored_ok = scored_ok,
        null_ok = health.null_ok,
        null_scored_fraction = health.null_scored_fraction,
        null_scored_ok = null_scored_ok,
        null_arm_ok = health.null_ok && null_scored_ok,
        n_injected = health.n_injected)
end

"""
    print_correction_accuracy_summary(io, decision)

Print the harness SMOKE-ONLY / VERDICT block for `decision` to `io`. The header
is `decision.line`, which includes `null_scored_fraction`. This is the print
the accuracy benchmark uses after `correction_accuracy_verdict`, so a synthetic
cell table can drive the same text without running the sweep.
"""
function print_correction_accuracy_summary(io::IO, decision)
    println(io, "\n$(decision.line)")
    if decision.verdict_allowed
        println(io, "Interpretation per cell:")
        println(io, "  recall high + precision high + over_correction_rate ~0 => corrector fixes errors cleanly")
        println(io, "  precision < recall via mis_fixes => failed corrections; via over_corrections => collateral damage")
        println(io, "  correction_rate should track error_rate (edits ~ true errors); >> error_rate => over-editing")
    else
        println(io, "This run MUST NOT be quoted as validation. " *
                    (decision.scale_ok ? "" : "Raise coverage/genome size. ") *
                    (decision.scored_ok ? "" :
                     "Too few reads scored (mass drops/exclusions) — investigate before trusting metrics. ") *
                    (decision.null_arm_ok ? "" :
                     "Control B null arm is unhealthy on injected-error cells (crash, mis-join, or mass drop). "))
        println(io, "Run the full Lambda-phage config on Lovelace for a VERDICT.")
    end
    return nothing
end
