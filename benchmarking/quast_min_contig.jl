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
#
# SOURCE OF THE 11,622 bp FIGURE: the dispatch note on bead td-28o0, recording
# Lawrencium job 24247925 (RGV T4 NC_000866 @ 30x). The cluster log is not
# reachable from a checkout, so that number is REPORTED, not reproduced here. The
# defect does not rest on it: the threshold arithmetic (168,903 / 10 = 16,890) is
# checkable in-repo, and `t4_ksweep.jl`'s own committed results contain T4 rows
# with largest contigs of 174 bp and 5,800 bp — local proof that assemblies of
# this organism land far below a 16,890 bp floor.
#
# THE FIX: clamp at QUAST's own default, which is what the term was for
# ---------------------------------------------------------------------
# `quast_min_contig` keeps the down-scaling for tiny references (the original,
# correct intent) and clamps it from above at QUAST's own default:
#
#   min_contig = clamp(genome_len ÷ 10, 50, 500)
#
# The scaling term exists ONLY to go below 500 for references shorter than QUAST's
# default assumes. Letting it go above 500 was never the intent, and every value
# it produced above 500 was a threshold nobody chose.
#
# WHY NOT a 5,000 bp ceiling (this file's first attempt, corrected 2026-07-27)
# ---------------------------------------------------------------------------
# The first version clamped at 5,000 bp, anchored to the sweep's longest read
# length, specifically to hold Lambda at 4,850 and so "preserve comparability with
# results already committed". Both halves of that rationale were checked against
# primary sources during review and neither survived:
#
#   * THERE ARE NO COMMITTED LAMBDA NUMBERS TO PRESERVE. `git log --all
#     --diff-filter=A` over `benchmarking/results/rhizomorph_correction_validation_sweep_*.csv`
#     returns exactly one file (commit 548dc984d): two rows of `synthetic_2000bp`
#     with no QUAST columns at all, whose threshold is 200 under every candidate
#     policy. The comparability being defended did not exist.
#
#   * LAMBDA WAS NOT "ALREADY WORKING" AT 4,850. Bead `td-4e19d.1` records the
#     completed Lovelace Lambda sweep: QUAST failed on 8 of 12 cells
#     (`metric_source=internal:quast-failed`) because the largest contigs at
#     err>=0.05 were 163-1,939 bp against a 4,850 bp threshold. Holding Lambda at
#     4,850 preserved the failure, not the measurement. At 500 the long-read
#     err=0.05 cells (largest 1,100 and 1,939 bp) become scorable; the short-read
#     err=0.10 cells (163-277 bp) remain unscored, which is correct — those
#     assemblies are genuinely below any usable threshold.
#
# A third consequence, invisible until the guard and the policy were considered
# together: a ceiling above 500 leaves `quast_min_contig` VARYING BY REFERENCE
# (Lambda 4,850 / phi29 1,928 / SARS-CoV-2 2,990 / T4 5,000). The
# pre-registration pools across organisms, so `metric_source_guard.jl` would
# refuse the pooled analysis as mixed-definition, and the only escape would be the
# `allow_mixed_src` override. Clamping at 500 makes the threshold UNIFORM across
# every reference above 5 kb, so the pooled analysis is legal by construction
# rather than by override. See beads td-28o0 / td-9p91.
#
# WHY NOT derive the threshold from the observed contig-length distribution
# ------------------------------------------------------------------------
# (e.g. `min(glen ÷ 10, largest_observed_contig)`). This makes the metric
# definition a function of the assembly being scored: the naive and iterative arms
# of the SAME cell would be filtered at different thresholds, so any
# naive-vs-iterative delta would confound a real effect with a definition change.
# That is precisely the mixed-definition failure `td-9p91` exists to reject,
# hidden inside a single cell. Rejected.
#
# The chosen policy is deterministic, depends only on the reference, and is
# identical for every arm of every cell — so arms stay comparable by construction.
# Callers should also RECORD the value they used (the RGV sweep writes it to the
# `quast_min_contig` CSV column) so the metric definition travels with the data and
# any future change to this policy is detectable by `metric_source_guard.jl`
# instead of silent.
#
# RESULTING THRESHOLDS (regression-pinned in the unit test):
#
#   reference     genome_len   old inline   now
#   viroid               300           50    50   unchanged (floor; intent preserved)
#   phi29             19,282        1,928   500
#   SARS-CoV-2        29,903        2,990   500
#   Lambda            48,502        4,850   500   recovers cells 4,850 zeroed out
#   T4               168,903       16,890   500   was FAILING outright
#   E. coli        4,641,652      464,165   500   would have failed
#
# Every reference above 5 kb now shares one threshold, which is the property the
# pooled pre-registered analysis needs.

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
Absolute upper bound in bp: QUAST's own `--min-contig` default. The scaling term
exists only to go BELOW this for viroid-scale references; it must never go above
it. See the file header for why an earlier 5,000 bp ceiling was wrong.
"""
const MIN_CONTIG_CEILING_BP = 500

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
- `ceiling_bp`: absolute upper bound (default `MIN_CONTIG_CEILING_BP` = 500, QUAST's own default).

# Returns

- `Int` threshold in bp, in `[floor_bp, ceiling_bp]`.

# Example

```julia
quast_min_contig(48_502)   # Lambda -> 500  (ceiling; was 4850 inline)
quast_min_contig(168_903)  # T4     -> 500  (ceiling; was 16890, which made QUAST fail)
quast_min_contig(300)      # viroid -> 50   (floor; the down-scaling term's real purpose)
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
