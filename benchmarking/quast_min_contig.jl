# QUAST `--min-contig` threshold policy.
# ======================================
#
# Pure, dependency-free (like `rhizomorph_scale_guard.jl` and
# `quast_report_parsing.jl`): no Mycelia, no network, no external tools. Included
# by the benchmark harnesses AND by the unit test, so the threshold policy has
# exactly one definition and a regression test can pin it.
#
# WHY THIS FILE EXISTS (bead td-28o0)
# -----------------------------------
# The harnesses computed the threshold inline as `max(50, glen ÷ 10)`. That
# expression was introduced to LOWER QUAST's threshold for viroid-scale
# references, where QUAST's 500 bp default excludes the entire assembly — the
# surviving comment at its sibling call sites still reads "Use low min_contig for
# viroid-scale genomes (246-359 nt)".
#
# But `glen ÷ 10` is monotone in genome length, so above 5 kb it does the
# OPPOSITE of its stated purpose: it RAISES the threshold above QUAST's default.
# For T4 (`NC_000866`, 168,903 bp) it demands a 16,890 bp contig. The naive arm's
# largest contig is 11,622 bp, so QUAST exits with
#
#   "None of the assembly files contains correct contigs. Please, provide
#    different files or decrease --min-contig threshold."
#
# and the harness degrades to internal size-ratio metrics
# (`metric_source="internal:quast-failed"`) while still reporting `ok=true`.
# Observed live on Lawrencium job 24247925 (RGV T4 @ 30x).
#
# THE FIX: a CEILING, not a different scaling
# -------------------------------------------
# `quast_min_contig` keeps the down-scaling for tiny references (the original,
# correct intent) and clamps it from above so the threshold can never exceed a
# length that a legitimately fragmented assembly is unable to reach.
#
#   min_contig = clamp(genome_len ÷ 10, 50, 5_000)
#
# WHY 5,000 bp as the ceiling. It is the longest read length on the sweep's own
# read-regime axis (`MYCELIA_RGV_READLEN` defaults to `150,5000`). A contig
# shorter than a single read carries no assembly evidence, so "at least as long
# as the longest read" is the strongest filter that is still a filter. Past that
# point the threshold stops excluding noise and starts requiring a particular
# DEGREE of assembly success — which is the quantity under measurement, not a
# precondition for measuring it.
#
# WHY NOT the two rejected alternatives (both weighed explicitly):
#
#   * Clamp to QUAST's own default of 500 bp. Most principled reading of the
#     original intent, but it MOVES every reference above 5 kb: Lambda
#     4,850 -> 500, phi29 1,928 -> 500, SARS-CoV-2 2,990 -> 500. That silently
#     redefines the metric for genomes that were already working and breaks
#     comparability with results already committed and with the Lambda jobs in
#     flight (Lawrencium 24247923 / 24247924). Rejected on comparability grounds,
#     not on principle.
#
#   * Derive the threshold from the observed contig-length distribution
#     (e.g. `min(glen ÷ 10, largest_observed_contig)`). This makes the metric
#     definition a function of the assembly being scored: the naive and iterative
#     arms of the SAME cell would be filtered at different thresholds, so any
#     naive-vs-iterative delta would confound a real effect with a definition
#     change. That is precisely the mixed-definition failure `td-9p91` exists to
#     reject, hidden inside a single cell. Rejected.
#
# The chosen policy is deterministic, depends only on the reference, and is
# identical for every arm of every cell — so arms stay comparable by
# construction. Callers should also RECORD the value they used (the RGV sweep
# writes it to the `quast_min_contig` CSV column) so the metric definition
# travels with the data and any future change to this policy is detectable by
# `metric_source_guard.jl` instead of silent.
#
# PRESERVED VALUES (regression-pinned in the unit test):
#
#   reference     genome_len   before      after
#   viroid               300       50         50   unchanged
#   phi29             19,282    1,928      1,928   unchanged
#   SARS-CoV-2        29,903    2,990      2,990   unchanged
#   Lambda            48,502    4,850      4,850   unchanged  <- comparability
#   T4               168,903   16,890      5,000   was FAILING
#   E. coli        4,641,652  464,165      5,000   would have failed

"""
Divisor applied to the reference length: the historical `genome_len ÷ 10`
down-scaling that exists to bring the threshold BELOW QUAST's 500 bp default for
viroid-scale references.
"""
const MIN_CONTIG_GENOME_DIVISOR = 10

"""
Absolute lower bound in bp. Keeps the threshold positive (and QUAST happy) for
references so short that `genome_len ÷ 10` would round toward zero.
"""
const MIN_CONTIG_FLOOR_BP = 50

"""
Absolute upper bound in bp, anchored to the longest read length on the sweep's
read-regime axis (`MYCELIA_RGV_READLEN` default `150,5000`). Above this the
threshold would demand a particular degree of assembly success rather than
filtering sub-read-length noise. See the file header for the full rationale.
"""
const MIN_CONTIG_CEILING_BP = 5_000

"""
    quast_min_contig(genome_len; divisor, floor_bp, ceiling_bp) -> Int

QUAST `--min-contig` threshold for a reference of `genome_len` bp.

Scales the threshold down for very short references and clamps it into
`[floor_bp, ceiling_bp]` so it can never exceed a contig length that a
legitimately fragmented assembly cannot produce. The result depends ONLY on the
reference, never on the assembly being scored, so every arm of a cell is
filtered identically.

# Arguments

- `genome_len`: reference length in bp (must be non-negative).

# Keywords

- `divisor`: down-scaling divisor (default `MIN_CONTIG_GENOME_DIVISOR` = 10).
- `floor_bp`: absolute lower bound (default `MIN_CONTIG_FLOOR_BP` = 50).
- `ceiling_bp`: absolute upper bound (default `MIN_CONTIG_CEILING_BP` = 5,000).

# Returns

- `Int` threshold in bp, in `[floor_bp, ceiling_bp]`.

# Example

```julia
quast_min_contig(48_502)   # Lambda      -> 4850  (unchanged from `max(50, glen ÷ 10)`)
quast_min_contig(168_903)  # T4          -> 5000  (was 16890, which made QUAST fail)
quast_min_contig(300)      # viroid      -> 50    (floor)
```
"""
function quast_min_contig(genome_len::Integer;
        divisor::Integer = MIN_CONTIG_GENOME_DIVISOR,
        floor_bp::Integer = MIN_CONTIG_FLOOR_BP,
        ceiling_bp::Integer = MIN_CONTIG_CEILING_BP)
    if genome_len < 0
        throw(ArgumentError("genome_len must be non-negative, got $genome_len"))
    end
    if divisor < 1
        throw(ArgumentError("divisor must be >= 1, got $divisor"))
    end
    if floor_bp < 1
        throw(ArgumentError("floor_bp must be >= 1, got $floor_bp"))
    end
    if ceiling_bp < floor_bp
        throw(ArgumentError(
            "ceiling_bp ($ceiling_bp) must be >= floor_bp ($floor_bp)"))
    end
    return Int(clamp(genome_len ÷ divisor, floor_bp, ceiling_bp))
end
