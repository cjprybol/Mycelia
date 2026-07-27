# opt4 frozen-read-skip — representative-scale accuracy sign-off

td-jbjd, Rule-of-5 pass 5 ("Verification"). Design doc:
`docs/design/2026-07-26-opt4-frozen-read-skip-design.md` ("Accuracy sign-off
methodology"). This is the load-bearing measurement for the maintainer's
enable/keep-off decision on `skip_frozen_reads` (opt-in, default off). It does
**not** make that decision — it reports the accuracy/speed tradeoff at
representative scale so the maintainer can apply their own tolerance.

Script: `benchmarking/opt4_frozen_read_accuracy_signoff.jl` (new, does not
mutate `rhizomorph_correction_accuracy_benchmark.jl`; reuses its
`acquire_reference` / `simulate_substitution_reads` / `per_base_metrics`
building blocks and adds the `skip_frozen_reads` passthrough that harness does
not expose). Config-tuning probes (kept for reproducibility, not part of the
sign-off itself): `benchmarking/opt4_signoff_probe.jl`,
`benchmarking/opt4_signoff_supplementary_probe.jl`.

## Why representative scale, not the pass-4 toy fixtures

`test/4_assembly/corrector_opt4_positive_control_test.jl` (pass 4) found that
at toy scale (900bp–2kb synthetic genomes, tens to ~150 reads), freezing is
**mechanically inert**: Stage 0 cheap correction (`_stage0_cheap_correct`,
`src/iterative-assembly.jl:4688`) is freeze-immune and runs unconditionally
before the freeze gate, and at toy scale it accounts for effectively all
correction volume — the per-read Viterbi decode (the only thing freezing
actually gates) made **zero edits** at any k tested. Separately, the adaptive
low-k decode gate (`decode_gate_density`, default 0.90,
`src/iterative-assembly.jl:4831`) additionally disables decode entirely on toy
graphs because a small/dense graph cannot achieve genuine (<90%)
decode-selectivity. Freezing an operation that never runs cannot be measured to
degrade anything — that is a finding about mechanism, not a clean bill of
health for real data.

This pass therefore targets a scale where decode genuinely fires, confirmed
directly against `result[:metadata][:decode_fraction_mean]` /
`:decode_fraction_max` (the fraction of reads that actually reached the
per-read Viterbi decode each pass, `= 1 - skip_fraction`) and
`total_improvements - cheap_corrections_total` (edits attributable to decode,
not Stage 0).

## Configuration reached

Reference: **NC_001416 (Escherichia phage Lambda, 48,502 bp)**, downloaded via
`Mycelia.download_genome_by_accession`. Substitution-only injected errors
(`Mycelia.mutate_dna_substitution_fraction`, length-preserving) at 2% rate,
read length 150bp, uniform Phred 20. Corrector knobs mirror the wired
`:scalable` production tier exactly (`skip_solid=true`,
`graph_mode=:doublestrand`, `n_k_rungs=3`, `max_iterations_per_k=2`,
`hard_window=true`, `soft_em=true`, `cheap_correct=true`, `k=21`,
`batch_size=500`) — i.e. the sign-off measures the actual knobs
`assemble_genome(corrector=:iterative, strategy=:scalable)` would use in
production, not a bespoke configuration built to make decode fire artificially.

Two coverage tiers were run (both reuse the identical simulated read set per
tier, seed 42, across all four arms):

| Tier | Coverage | Reads | Injected errors | Sequenced bases (cov × genome) | vs. harness's own `SCALE_FLOOR_BASES=1,000,000` |
|---|---|---|---|---|---|
| A (smoke) | 8.0x | 2,587 | 7,761 | 388,016 | below floor — SMOKE-ONLY by the harness's own guard |
| B (verdict) | 21.0x | 6,791 | 20,373 | 1,018,542 | **at/above floor — VERDICT scale** |

**Evidence decode is active** (not gated to zero, unlike the pass-4 toy
fixtures):

| Tier | decode_fraction_mean | decode_fraction_max | decode-attributable edits (total_improvements − cheap_corrections_total) |
|---|---|---|---|
| A | 0.385 | 0.837 | 1,356 |
| B | 0.256 | 0.623 | 2,981 |

Both tiers show a substantial, nonzero fraction of reads reaching the
expensive per-read Viterbi decode every pass, and a real block of edits
attributable to decode rather than Stage 0. This is the prerequisite the pass-4
header called out as missing at toy scale.

**Scale caveat:** Tier B clears the harness's own 1 Mbase `SCALE_FLOOR_BASES`
guard (`benchmarking/rhizomorph_scale_guard.jl`) — the same floor
`rhizomorph_correction_accuracy_benchmark.jl` uses to gate a VERDICT vs.
SMOKE-ONLY label — but this is still a single genome (Lambda phage, no
repeats/structural complexity beyond what a 48.5 kb phage genome has), a single
seed, and a single machine with no repeated-run averaging. Do not extrapolate
past 21x coverage / this genome without re-running.

## Three/four-arm results

Four arms on the **identical** read set per tier (same simulated reads/errors,
same seed — only `skip_frozen_reads` / `freeze_streak_threshold` /
`freeze_across_rungs` vary):

- **exact** — `skip_frozen_reads=false` (the always-on exact path)
- **freeze_within** — `skip_frozen_reads=true, threshold=2, across=false` (the
  **default** scope per the design doc)
- **freeze_across** — `skip_frozen_reads=true, threshold=2, across=true`
  (aggressive scope, conservative threshold)
- **freeze_max_aggressive** — `skip_frozen_reads=true, threshold=1,
  across=true` (the positive-control candidate: most aggressive setting tried)

### Tier A (8x coverage, 2,587 reads, SMOKE-ONLY by the scale guard)

| Arm | Runtime (s) | Speedup vs exact | frozen_reads_skipped | decode_fraction_mean | Recall | Precision | Over-corr. rate | True fixes | Byte-identical to exact? |
|---|---|---|---|---|---|---|---|---|---|
| exact | 46.29 | — | 0 | 0.385 | 0.8382 | 0.8441 | 0.00308 | 6,505 | — |
| freeze_within (default) | 32.20 | +30%* | **0** | 0.385 | 0.8382 | 0.8441 | 0.00308 | 6,505 | **yes** (0/2587 differ) |
| freeze_across | 30.63 | +34% | 724 | 0.327 | 0.8356 (Δ −0.0026) | 0.8445 (Δ +0.0004) | 0.00307 | 6,485 (Δ −20) | no (14/2587 differ) |
| freeze_max_aggressive | 28.29 | +39% | 2,258 | 0.240 | 0.7977 (Δ **−0.0405**) | 0.8706 (Δ +0.0265) | 0.00236 | 6,191 (Δ −314) | no (244/2587 differ) |

\* freeze_within skipped **zero** reads at this config (see mechanism note
below) — its reported "+30% speedup" is single-run wall-clock noise, not a real
effect of the freeze mechanism. It is reported for transparency, not as
evidence of speedup.

### Tier B (21x coverage, 6,791 reads, VERDICT scale)

| Arm | Runtime (s) | Speedup vs exact | frozen_reads_skipped | decode_fraction_mean | Recall | Precision | Over-corr. rate | True fixes | Byte-identical to exact? |
|---|---|---|---|---|---|---|---|---|---|
| exact | 123.59 | — | 0 | 0.256 | 0.9001 | 0.9779 | 0.000407 | 18,337 | — |
| freeze_within (default) | 111.11 | +10%* | **0** | 0.256 | 0.9001 | 0.9779 | 0.000407 | 18,337 | **yes** (0/6791 differ) |
| freeze_across | 95.84 | +22% | 1,276 | 0.213 | 0.8999 (Δ −0.0002) | 0.9779 (Δ ≈0) | 0.000407 | 18,333 (Δ −4) | no (4/6791 differ) |
| freeze_max_aggressive | 98.55 | +20% | 2,922 | 0.159 | 0.8907 (Δ **−0.0094**) | 0.9994 (Δ +0.0215) | 0.000004 | 18,146 (Δ −191) | no (150/6791 differ) |

\* same caveat as Tier A — freeze_within skipped zero reads here too.

Full per-arm CSVs:
`benchmarking/results/opt4_frozen_read_skip_signoff_20260727_134304.csv` (Tier
A), `benchmarking/results/opt4_frozen_read_skip_signoff_20260727_135800.csv`
(Tier B).

## Mechanism finding: the default (within-rung) scope never fired at either scale

`freeze_within` — the DEFAULT, "safe" scope per the design doc — skipped
**zero** reads at both coverage tiers, and produced byte-identical output to
`exact` in both. This is not an accuracy result; it means the freeze condition
had **no opportunity to trigger** under the production `:scalable` knob
`max_iterations_per_k=2`: within-rung scope resets each read's no-improvement
streak to 0 at every k-change, and a `threshold=2` freeze requires **2
consecutive** no-improvement passes at the *same* k before a 3rd pass can be
skipped — but with only 2 iterations allowed per k, the earliest a streak can
reach 2 is on the 2nd (last) iteration of that rung, by which point there is no
further same-k pass left to skip.

A supplementary probe (`benchmarking/opt4_signoff_supplementary_probe.jl`,
Tier-A read set, `max_iterations_per_k=4` instead of 2) confirms this
mechanism: with more per-rung headroom, within-rung freezing **did** engage
(724 reads skipped, matching Tier A's `freeze_across` count exactly — the 3rd
and 4th same-k iterations are where it fires), runtime dropped 41.3s → 31.5s
(−24%), and the accuracy cost was negligible: recall 0.8402 → 0.8400 (Δ
−0.0002), precision 0.8366 → 0.8372 (Δ +0.0006), true fixes 6,521 → 6,519 (Δ
−2).

**Implication:** at the actual production config (`max_iterations_per_k=2`),
enabling `skip_frozen_reads` with the documented default scope
(`freeze_across_rungs=false`) changes **nothing** — it is a true no-op, not
merely a low-risk option. Any speedup from opt4 under the current `:scalable`
production knobs requires the **aggressive** `freeze_across_rungs=true` scope,
which the design doc calls out explicitly as higher-risk and non-default.

## Positive control (this pass's requirement)

`freeze_max_aggressive` (`threshold=1, across=true`) is this pass's
positive-control candidate — the most aggressive setting the pass-4 toy fixture
could not distinguish from exact. At representative scale it **does**
measurably degrade recall:

- Tier A: recall −0.0405 absolute (−4.8% relative), 244/2,587 reads (9.4%)
  differ from exact.
- Tier B (VERDICT scale): recall −0.0094 absolute (−1.0% relative), 150/6,791
  reads (2.2%) differ from exact.

Both are asserted with a real, reproducible (fixed-seed) delta — this confirms
the per-base metric **can** detect opt4-induced degradation at representative
scale, resolving the pass-4 "cannot demonstrate degradation" finding, which was
scale-limited, not a property of opt4 itself.

**Notable asymmetry:** in both tiers, `freeze_max_aggressive` *also improves*
precision and collapses the over-correction rate (Tier A: over-corrections
1,173 → 898; Tier B: 406 → 4, over_correction_rate 0.000407 → 0.000004). The
mechanism: freezing a read early removes it from later, potentially
over-eager, decode passes as much as it removes it from passes that would have
made a genuine fix — the net effect observed here is fewer total edits
(lower `correction_rate` in both tiers) with a *shift* in error type from
"borderline over-correction" toward "missed correction" (lower recall), not a
uniform accuracy loss. This is a real tradeoff the maintainer should weigh, not
simply "worse."

`freeze_across` (the non-default but non-maximal `threshold=2, across=true`
scope) sits between the two: a small, real recall cost (Tier A −0.31%
relative, Tier B −0.02% relative) shrinking as scale/coverage grows, for a
~22–34% single-run speedup and a genuine reduction in decode volume
(`decode_fraction_mean` drops 15–17% relative to exact in both tiers,
independent of wall-clock noise).

## Assembly-level dnadiff/ANI check (design's second accuracy check)

The per-base check above measures the corrector's own read output. This second
check re-assembles the corrected reads into contigs and runs MUMmer `dnadiff`
against the Lambda reference — confirming (or challenging) the per-base −0.02%
result at the assembly level. Run on the **identical** Tier B fixture (21x, seed
42), for **exact** vs **freeze_across** (threshold 2, across-rung). Script:
`benchmarking/opt4_frozen_read_dnadiff_signoff.jl`. **Critical method note:** the
corrected reads are re-assembled with `corrector=:none`
(`Rhizomorph.assemble_genome`, no iterative correction), so the exact-vs-freeze
difference in the input reads propagates into the contigs rather than being
masked by re-correction.

| Arm (Tier B, 21x) | n_contigs | contig bases | largest contig | AvgIdentity | Ref aligned | TotalSNPs | TotalIndels |
| --- | --- | --- | --- | --- | --- | --- | --- |
| exact | 2,983 | 184,523 | 37,565 | 100.00% | 48,469 (99.93%) | 0 | 0 |
| freeze_across | 2,985 | 184,567 | 37,565 | 100.00% | 48,469 (99.93%) | 0 | 0 |
| _(raw, uncorrected — baseline)_ | 30,664 | 1,331,734 | 2,071 | 99.96% | — | — | — |

**Verdict — the assembly-level check confirms the per-base −0.02% finding.**
`exact` and `freeze_across` are **indistinguishable at the assembly level**:
identical 100.00% average identity, identical 99.93% reference coverage
(48,469/48,502 bp), and **zero SNPs and zero indels** in both. The only
difference is 44 query bases (0.02%) of redundant overlapping-contig sequence,
which carries no accuracy cost (identity, SNP, and indel counts are unchanged).
The raw-reads baseline (30,664 fragmented contigs, largest 2,071 bp) confirms
correction is doing substantial work and the assembler is functioning — so the
identical exact/freeze contig quality is a real result, not an artifact of a
no-op assembler. The `freeze_across` recall cost measured at the read level does
**not** propagate to any assembly-level accuracy loss on this fixture.

Caveat: the `corrector=:none` k=21 assembly is intentionally fragmented (no
scaffolding), so query aligned-fraction is ~52% (redundant overlapping contigs);
the load-bearing signal is the exact-vs-freeze equivalence (identical
identity/SNP/indel), not the absolute contiguity. Single genome/seed/machine, as
with the per-base check.

## Caveats

- **Single machine, single run per arm — no repeated-run variance estimate.**
  `freeze_within`'s reported "+10% to +30% speedup" occurred with **zero**
  freeze events, i.e. it is pure wall-clock noise on this box between
  otherwise-identical runs. Treat every raw runtime number in this doc as
  directional, not a precise measurement; `frozen_reads_skipped` and
  `decode_fraction_mean` (both deterministic, seed-driven counts, not
  wall-clock) are the more trustworthy speed proxies for `freeze_across` and
  `freeze_max_aggressive`, where real (nonzero) skip volume was observed.
- **Scale reached:** Tier B (21x coverage, Lambda phage, 48,502 bp) clears the
  harness's own 1 Mbase VERDICT floor, but is still one genome / one seed / one
  machine. Do not extrapolate to larger, more repetitive, or long-read genomes
  without re-running.
- **Assembly-level (dnadiff/ANI) check: now included** — see the
  "Assembly-level dnadiff/ANI check" section above. It confirms the per-base
  result (exact and freeze_across both reach 100.00% AvgIdentity, 0 SNPs, 0
  indels vs the Lambda reference; assemblies differ by 0.02% redundant bases
  only). Per-base accuracy against injected-error ground truth remains the more
  direct signal (it measures the corrector's own output, unmediated by a
  downstream assembler); the dnadiff check is confirmatory.
- **`path_to_sequence` warnings** emitted during every run are pre-existing
  synthetic-genome noise (per the task brief) and were filtered from all logs;
  they do not affect any of the numbers above.
- **`substitution_length_divergences`** (decode reconstructions that failed
  open to the original read to preserve the length contract) were nonzero in
  every arm/tier (roughly 0.2–3% of reads per pass) — this is pre-existing
  corrector behavior unrelated to opt4 (present identically in the `exact`
  arm) and does not bias the exact-vs-freeze comparison, since it affects both
  arms' shared decode path identically.

## Regression checks

All five existing opt4 test files still pass on this branch after the sign-off
run (source untouched — this pass added only new benchmarking scripts and this
results doc):

| Test file | Result |
|---|---|
| `corrector_opt4_frozen_read_skip_test.jl` | 2/2 pass |
| `corrector_opt4_freeze_state_test.jl` | 24/24 pass |
| `corrector_opt4_golden_master_lock_test.jl` | 2/2 pass |
| `corrector_opt4_positive_control_test.jl` | 6/6 pass |
| `corrector_opt4_edge_cases_test.jl` | 16/16 pass |

## Summary for the maintainer's ship decision

At this scale (Lambda phage, 8x–21x coverage, production `:scalable` knobs):

- The **shipped default** (`skip_frozen_reads=true` with `freeze_across_rungs=false`,
  the documented safe scope) is a **byte-identical no-op** under the current
  production `max_iterations_per_k=2` — it neither speeds anything up nor
  changes accuracy, because the freeze condition never has a chance to fire.
  Recommending it as-is provides zero benefit at the current config (though it
  remains harmless).
- The **aggressive scope** (`freeze_across_rungs=true`, still `threshold=2`)
  gives a real, nonzero decode-volume reduction (~15–17% fewer per-read
  decodes) for a recall cost that is small and **shrinks with scale**
  (−0.31% relative at 8x coverage, −0.02% relative at 21x coverage).
- The **positive control** (`threshold=1, across=true`) confirms the metric
  detects real degradation at this scale (recall −1.0% to −4.8% relative,
  precision improves), giving confidence the "no accuracy loss" verdicts above
  are not simply an insensitive measurement.

No enable/ship recommendation is made here — the tolerance for the
aggressive-scope recall cost vs. its decode-volume savings, and whether the
default scope's config-dependent inertness is acceptable to ship as-is or
warrants raising `max_iterations_per_k` for the default tier, are the
maintainer's calls per the design doc's ship gate.
