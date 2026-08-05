# ADR: joint quality aggregation in Rhizomorph (td-2tg8)

| Field          | Value                                                                                                                                                                       |
| -------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Status**     | **Proposed** — acceptance WITHDRAWN 2026-08-05 pending re-derivation of §4 (see Errata)                                                                                     |
| **Decision**   | Stop treating summed Phred as a calibrated probability. Keep it as a named traversal _heuristic_. Move production toward molecule-aware, dependence-discounted aggregation. |
| **Interim**    | **WITHDRAWN pending re-derivation — see Errata.** The empirical basis for preferring the conservative model was a censored statistic.                                       |
| **Scope**      | The opt-in `traversal_weighting = :quality` path, and every manuscript claim interpreting joint quality as a probability                                                    |
| **Supersedes** | The `100-255 = essentially certain` confidence band in `planning-docs/rhizomorph-graph-ecosystem-plan.md:498-550`                                                           |

## Errata — acceptance withdrawn 2026-08-05

This ADR was accepted, then un-accepted the same day after a nine-reviewer pass
on PR #453 falsified two of its supporting claims. **The central decision is
unaffected** — it rests on the calibration experiment, whose arithmetic was
independently audited against 512-bit reference values and confirmed correct to
15 digits across 115 orders of magnitude. What failed is §4, the section that
departs from the commissioned review.

**E1 — §4.2's comparison is computed on a censored subset, and the censoring
removes exactly the disconfirming evidence.** (3 independent reviewers.) The
calibration script emits `-Inf` when a bin has zero observed errors, so **49 of
125 bins were silently dropped** from every §4.2 statistic. The censoring is not
symmetric in consequence: a zero-error bin cannot falsify independence, which
predicts ~0 errors there, but it is precisely where the conservative model's
substantial predicted error rate _can_ be falsified. Scored against this ADR's
own rule of three (§6), six dropped bins exceed the stated ±1.8-order ceiling,
up to **2.38 orders** — and five of the six are random-regime, high-support
bins, which is verbatim the case the commissioned reviewer predicted
conservative would overpredict. **§4.2's rebuttal of that prediction is
therefore unsupported; the prediction is confirmed in the removed subset.** The
claim "conservative is never off by more than ~1.8 orders in either direction"
must be restated as holding only over bins with ≥1 observed error, or re-derived
with a one-sided bound for zero-error bins.

**E2 — the Stage-C "the decoder never ran" mechanism is false for Illumina.** (3
independent reviewers.) The `decoder_ever_ran` flag regex-matched a
_stringified_ vector with a `^`-anchor, so it only ever inspected pass 1 — which
is `0.0` by construction in every row, making the flag structurally incapable of
returning true. Illumina's actual `decode_fraction_per_pass` is
`[0.0, 0.058, 0.017, 0.014]`. The corrected finding is **stronger** than the one
it replaces: the decoder ran, ran on a _quality-dependent_ number of reads (3.6%
under flat quality vs 1.7% under oracle), and the assemblies were still
byte-identical. PacBio and ONT are genuinely all-zero.

**What this episode is an instance of.** §2.3 of this document reports catching
a fallback whose value coincided with the healthy answer, and generalises the
lesson. Both errata above are further instances of the same class that §2.3
claims to have learned: E2's flag failed to the "clean" value, and E1's `-Inf`
sentinel reads as the _safe_ end of a scale documented as "positive =
overconfident." Recognising a defect class is not the same as auditing for it.

## Context

Rhizomorph aggregates per-position Phred across read observations by summing
(`Q_combined = ΣQᵢ`, clamped to 255) — `combine_phred_scores` /
`get_vertex_joint_quality`. Summing in log space multiplies error probabilities,
which is valid only under **independence** of observations. Bead td-4e19d made
this load-bearing by adding an opt-in quality-weighted traversal; the default
greedy path remains quality-blind.

Two external reviews were commissioned (a decision-oriented ADR draft and a
citation-backed methods survey). Their claims were **verified before use**, and
both our own prior claims and theirs were revised as a result. What follows
separates _verified_, _measured_, and _taken on trust_.

## Decision

1. **Summed Phred is not a calibrated probability that a k-mer is correct, and
   must not be reported as one.** It may remain as an explicitly-named
   evidence/traversal score.
2. **The `100-255 = essentially certain` band is withdrawn.** It is not
   supportable — see §3 and §6.
3. **Production direction is dependence-aware aggregation**: collapse evidence
   to physical molecules, then discount residual within-source correlation, with
   independence as the fitted zero-dependence limit rather than an assumption.
4. **Interim, if a probability-like number is required before (3) lands:** use
   the conservative aggregation. It is theoretically estimating the wrong event
   (§4.1) but is empirically within ~1.8 orders of magnitude everywhere we
   measured, against up to 39 for independence.
5. **Joint confidence stops being a `UInt8`.** Accumulate in floating-point log
   space; derive a Phred-like value only at the reporting boundary, capped at
   the empirically supported range.
6. **Quality-weighted traversal stays experimental** until validated on real
   data with sample-specific truth, stratified by chemistry.

## 1. What we verified against primary sources

Both production precedents for dependence discounting are real. This matters
because it moves the recommendation from "statistically principled" to "what
shipping tools do."

| Claim                                                                  | Verified how                           | Result                                                                                                                                                                                                  |
| ---------------------------------------------------------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FreeBayes discounts successive observations                            | Fetched `freebayes/src/Parameters.cpp` | `RDF = 0.9; // -D --read-dependence-factor`, help text: _"Incorporate non-independence of reads by scaling successive observations by this factor during data likelihood calculations."_ **Confirmed.** |
| bcftools/htslib attenuates repeated same-base/same-strand observations | Fetched `htslib/errmod.c`              | `em->fk[n] = pow(1. - depcorr, n) * (1.0 - eta) + eta`, with `eta = 0.03`. Weight decays toward 3% as identical same-strand observations accumulate. **Confirmed.**                                     |

Reported by the surveys and **not independently verified here** (treat as
secondary): GATK's overlapping-mate handling and BQSR covariates (the GATK docs
returned HTTP 403); BayesHammer's independence assumption; DeepVariant/medaka
architecture details; Racon/SPOA and Pilon evidence sums; Merqury's QV
derivation. The two surveys disagree slightly on GATK's overlapping-mate
mechanism (Q20 cap vs. cap at half the PCR error rate), which is one reason it
is listed as unverified.

**Relevant precedent for our actual usage:** Racon/SPOA reportedly sums decoded
Phred as partial-order-graph _edge weights_, and Pilon sums `Q × (MQ+1)` as
allele evidence. Neither is a posterior. That is a good description of what our
traversal gate legitimately is — quality-weighted voting — and a poor
description of a confidence band.

## 2. Corrections to our own prior claims

Both surveys found real errors in our briefing. Both are accepted.

**2.1 We over-generalised across chemistries.** The briefing stated that
random-regime bins with ≥3 observations contain no errors. False: six bins do,
reaching **0.465** for pacbio at 40x and 0.186 for ont at 40x in the 3-4 bin.
The claim was drawn from the Illumina rows and generalised — the precise pooling
error this project repeatedly warns against. Corrected in place. Consequence:
for long reads the "random" regime is **not** a clean independent-error control,
because random long-read error recurs measurably.

**2.2 `n_observations` counts evidence entries, not molecules.** Confirmed by
reading: both the calibration script and production `get_vertex_joint_quality`
iterate `observation_id => [entries]` and push **every** `QualityEvidenceEntry`.
So the model's `n` is an entry count. _However_, we then measured
`entries per observation_id = 1.0` in every bin of this simulator, so this did
**not** inflate the results reported here. It remains a live risk for
paired-end, duplicate-bearing, or multi-entry real data, where one molecule
could contribute several entries and be counted as independent evidence.

**2.3 Two bugs we found ourselves while fixing the above.** Both surfaced only
because a result looked too clean, which is worth recording as a pattern.

_Catastrophic cancellation._ The literal `1 - ∏(1-p)` collapses: once
`p < ~1e-16`, `(1-p)` rounds to exactly `1.0`, the product is exactly `1.0`, and
the union probability becomes a hard `0.0` — manufacturing apparent _infinite_
overconfidence (we briefly recorded fake values near 10^297). Verified at
Q300×31: naive returns `0.0`, stable `-expm1(Σ log1p(-p))` returns `3.1e-29`.
All numbers below use the stable form.

_A silent fallback that reported a broken counter as a clean result._ The
entries-per-observation counter in §2.2 was first computed as
`n_observation_ids == 0 ? 1.0 : ratio`. A dropped assignment left the id count
at zero everywhere, so every bin took the `1.0` fallback — and `1.0` is
precisely the "no problem here" value. The measurement appeared to confirm that
entries equal observations while in fact measuring nothing. Fixed to emit `-1.0`
on a broken count so the failure is loud, then re-measured. **The 1.0 in §2.2 is
now a real measurement**; it was not when first observed. A fallback whose value
coincides with the expected healthy answer is indistinguishable from success —
the same class of defect as the manuscript guard that could never fire.

**2.4 Estimand mismatch — valid, but empirically negligible.** The reviews
correctly noted that we scored the _weakest-position_ value against
_whole-k-mer_ truth, whereas `P(k-mer wrong) = 1 - ∏ⱼ(1-pⱼ)`, bounded above by a
factor of k (log₁₀31 = 1.49). We now compute both. Measured difference: **median
0.05, max 1.09 log₁₀.** So the estimand error accounts for at most ~1 of up to
39 orders of magnitude. **The headline finding survives this objection intact**
— the reviews' strongest methodological criticism does not explain the result.

## 3. What the corrected measurement shows

Simulated phiX174, k=31, three chemistries, never pooled. `random` = independent
substitutions; `artifact` = one fixed reference locus mis-called identically in
every read covering it, at low reported quality.

**Independence fails in proportion to correlation, and the failure grows with
depth.**

| chemistry | coverage | obs bin | observed error | independence off by | conservative off by |
| --------- | -------- | ------- | -------------- | ------------------- | ------------------- |
| illumina  | 40x      | 17+     | 5.9e-3         | **10^39.3**         | 10^-0.28            |
| illumina  | 20x      | 17+     | 1.3e-2         | 10^21.1             | 10^0.00             |
| illumina  | 10x      | 9-16    | 8.3e-3         | 10^11.3             | 10^0.13             |
| pacbio    | 40x      | 17+     | 1.3e-2         | 10^12.2             | 10^-0.95            |
| ont       | 40x      | 17+     | 0.143          | 10^7.4              | 10^-0.44            |

**Saturation is a separate, orthogonal failure.** Fraction of _true_ k-mers
pinned at the 255 clamp: illumina 62% @10x, **97.6% @20x**, 99% @40x (median
unclamped joint 1147 at 40x); pacbio 66% @40x; ont 3.4% @40x. At ordinary
Illumina depth the gate is effectively a coverage gate. This is a representation
limit, not a modelling error, and would persist under any unbounded aggregation.

## 4. Where we depart from the external ADR

The commissioned ADR rejected the conservative model outright. We accept its
_theoretical_ argument and reject its _practical_ conclusion, on our own
measurements.

**4.1 Its theoretical objection is correct.** `1 - ∏(1-pᵢ)` is the probability
that **at least one contributing observation is erroneous**, not that the
**consensus label** is wrong. It also degrades with coverage: 100 agreeing Q40
observations score ~Q20, worse than one. As a general-purpose confidence measure
that is incoherent.

**4.2 But it is dramatically better calibrated than independence, in both
regimes:**

| regime   | independence (median / max) | conservative (median / max / min) |
| -------- | --------------------------- | --------------------------------- |
| random   | 1.14 / **3.80**             | 0.49 / 1.20 / -1.83               |
| artifact | 1.92 / **39.25**            | 0.41 / 1.20 / -1.35               |

Conservative is never off by more than ~1.8 orders in either direction;
independence reaches 39. The ADR predicted conservative would "substantially
overpredict error in many random high-support bins" — measured, its worst
overprediction is 1.83 orders, against independence's 39-order underprediction.

**4.3 Why it probably tracks so well here, stated as a hypothesis not a
finding.** When errors are _systematic_, "at least one observation is wrong" and
"the consensus is wrong" largely coincide — if the cause is shared, one wrong
implies all wrong. The conservative formula may be accidentally well-matched to
the correlated regime for that reason. That is a hypothesis this experiment does
not test, and it would not transfer to regimes where errors are independent and
the consensus is right despite individual errors.

**4.4 Practical consequence.** Neither documented model is the destination, and
the ADR's recommended ESS-tempered model is **untested here**. But between the
two options that exist today, independence is the worse one by a wide margin,
and it is the one currently shipping.

## 5. What we are NOT concluding

- **Not** that these magnitudes transfer to real data. The error model is
  synthetic, the artifact is idealised (uniform low quality), and real
  systematic artifacts frequently are **not** flagged by low instrument Phred —
  which would make the failure harder to detect, not easier.
- **Not** that quality distinguishes repeat collapse. Our truth test asks only
  whether a k-mer occurs _anywhere_ in the reference, so a real sequence in the
  wrong repeat copy or at the wrong copy number scores as true. The design's
  claim that quality separates "repeat-collapsed k-mers" from real sequence is
  **untested** and should be narrowed in the manuscript to recurrent
  _sequencing-error_ k-mers.
- **Not** that chemistry ordering (illumina ≫ pacbio ≫ ont) reflects real
  platforms. It primarily reflects the per-observation Phred assigned in our
  profiles.
- **Not** that the conservative model is correct — only that it is less wrong
  (§4).

## 6. Empirically supportable ceilings

With zero observed errors in a bin of N independent trials, the one-sided 95%
bound is ~3/N. So Q60 requires ~3×10⁶ independent error-free opportunities; Q100
requires ~3×10¹⁰. Overlapping k-mers are **not** independent trials. Scores
above the validated range should be reported as `Q ≥ x`, never extrapolated.
`Q255` is not a meaningful claim.

## 7. Manuscript obligations

Must state: the default greedy path is quality-blind and produced every
committed benchmark row; quality-weighted traversal is opt-in and unvalidated
for assembly improvement; the aggregation model and its independence assumption;
that the score is a traversal heuristic, not a calibrated posterior; that
results are chemistry-stratified and never pooled; that the evidence unit is an
entry, not a verified molecule.

Remove or qualify: "essentially certain"; and the claim that quality
distinguishes **repeat-collapsed** k-mers, which our evidence does not support
(§5).

## 8. What would change this decision

Toward simple independence: after molecule collapse and recalibration, fitted
within-source dependence ≈ 0; calibration holds across observation count and
coverage; duplicating same-source observations does not inflate confidence.

Toward a richer hierarchical model: dependence varies by locus/context beyond
source-level shrinkage; recurrent high-quality false vertices survive molecule
collapse; independent libraries share artifacts.

## 9. Next actions

1. Real-data calibration on GIAB HG002 with sample-specific (not GRCh38) truth,
   stratified by chemistry and difficult-region class. **Blocking for any
   production default.**
2. Falsification tests: duplicate invariance, coverage invariance,
   source-diversity contrast, count-dependent overconfidence trend.
3. Implement molecule collapse so `n` counts molecules, not entries (§2.2).
4. Move joint confidence off `UInt8` (§3 saturation).
5. Implement `aggregate_quality_scores_conservative` — it has never existed in
   `src/`.
6. Rare-variant retention guard: any dependence discount must be checked against
   loss of genuine low-frequency variation, per the design's own §4.4 warning.
