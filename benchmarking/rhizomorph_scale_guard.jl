# Rhizomorph correction-validation SCALE GUARD (pure helpers, no dependencies).
#
# WHY this file is import-free: the scale guard is the one piece of the
# validation sweep that must be unit-testable without loading Mycelia,
# downloading a genome, or running an assembly. Keeping it dependency-free lets
# `rhizomorph_correction_validation_sweep_test.jl` exercise the SMOKE-ONLY vs
# VERDICT boundary in milliseconds, and lets the full harness `include` the same
# definitions so the guard used in tests is byte-identical to the guard used in
# production.
#
# The guard exists to make a toy run UNQUOTABLE as validation: a VERDICT is
# only emitted when there is enough sequencing (coverage x effective genome
# length) to trust an assembly-vs-assembly comparison. Below the floor the
# harness prints a SMOKE-ONLY notice instead of a VERDICT.

# Minimum total sequenced bases (effective coverage x effective genome length)
# required before a VERDICT may be printed. 1 Mbase is comfortably above any
# toy/self-test configuration (e.g. a 2 kb synthetic genome at 10x = 20 kb) and
# comfortably below a realistic Lambda-phage validation run
# (48,502 bp x 30x ~ 1.46 Mbase). Override with MYCELIA_RGV_SCALE_FLOOR.
const SCALE_FLOOR_BASES = 1_000_000

"""
    scale_verdict_allowed(effective_coverage, effective_genome_length; floor=SCALE_FLOOR_BASES) -> Bool

Return `true` iff `effective_coverage * effective_genome_length` (the total
number of sequenced bases backing the sweep) meets or exceeds `floor`.

This is the single decision point that gates whether the harness may print a
VERDICT. A toy-scale configuration (tiny genome and/or low coverage) falls
below the floor and yields `false`, forcing a SMOKE-ONLY notice so the run can
never be misquoted as validation.
"""
function scale_verdict_allowed(effective_coverage::Real, effective_genome_length::Real;
        floor::Real = SCALE_FLOOR_BASES)
    return effective_coverage * effective_genome_length >= floor
end

"""
    scale_metric_bases(effective_coverage, effective_genome_length) -> Float64

Total sequenced bases backing the sweep — the quantity compared against the
scale floor. Surfaced separately so the harness can report it in both the
SMOKE-ONLY notice and the VERDICT header.
"""
function scale_metric_bases(effective_coverage::Real, effective_genome_length::Real)
    return Float64(effective_coverage) * Float64(effective_genome_length)
end
