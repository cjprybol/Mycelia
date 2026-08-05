# Deep-research context: is the independent-observation assumption in Rhizomorph's joint quality defensible?

**Prepared** 2026-08-05 · **Tracking** td-2tg8 (decision), td-8zle (supporting
audit) **Repo** Mycelia (Julia), branch `cp/td-4e19d-2-qualmer` · **Project**
Rhizomorph, a probabilistic graph-based genome assembler **Deliverable wanted**
an ADR-style decision: keep independence / adopt the conservative model / adopt
a middle model — justified by evidence, with an empirical calibration check that
would detect the assumption failing on real data.

---

## 1. Project summary

Rhizomorph is the assembly system inside Mycelia. Its distinguishing claim
against ordinary de Bruijn assembly is that it is **quality-aware**: it builds
"qualmer" graphs whose vertices (k-mers) carry per-base Phred quality from every
read observation, not just an occurrence count.

The philosophical case for this, stated in the project's design document, is
worth quoting because the whole research question turns on it:

> Technical artifacts and repeat-collapsed k-mers **also** recur, so an evidence
> count cannot distinguish them from biologically real sequence — but they recur
> at **low quality**. Compounding observations of high-quality k-mers should
> therefore be worth strictly more than compounding observations alone.

The mechanism implementing that claim is **additive Phred aggregation**
(`planning-docs/rhizomorph-graph-ecosystem-plan.md:498-550`):

```
Q_combined = Q₁ + Q₂ + … + Qₙ      (clamped to 255)

Documented confidence bands:
  10-30    low        a few low-quality observations
  30-60    moderate   a typical SINGLE high-quality observation
  60-100   high       multiple high-quality observations
  100-255  extreme    many high-quality observations = essentially certain
```

Implemented as `combine_phred_scores` / `get_vertex_joint_quality`
(`src/rhizomorph/core/quality-functions.jl`).

**This is the assumption under review.** Phred adds in log space, so summing
Phred is _multiplying error probabilities_ — valid only if the observations are
**independent**. The design names this explicitly as the
`aggregate_quality_scores_independence` model.

---

## 2. Why the assumption is suspect

The concern is not abstract, and it is self-undermining in a specific way:

Technical artifacts recur across reads _by construction_ — that is what makes
them artifacts rather than random errors. Their errors are therefore
**correlated**, not independent. Summing over correlated observations overstates
confidence, and it overstates fastest exactly where artifacts recur — **the very
case the quality channel exists to catch.**

Known correlated error sources in real sequencing: systematic/context-specific
miscalls (homopolymers, GC extremes), strand bias, tile/flowcell effects, PCR
duplicates, and instrument Phred that is itself miscalibrated and
context-dependent.

The design **anticipates this** and specifies an alternative,
`aggregate_quality_scores_conservative`:

```
P(all correct) = ∏ Pᵢ(correct)
Q_combined     = -10·log₁₀(1 - P(all correct))
```

with the note: _"Use this when you want to be conservative and account for the
possibility that multiple observations might share systematic errors."_

**That alternative has no implementation anywhere in `src/`** — it exists only
in the design document (verified by audit td-8zle). So today the codebase can
express only the independence model; whichever way this decision goes, the
alternative has to be written.

---

## 3. What we measured, and what it showed

Two experiments were run in this repo. Both are committed and reproducible
without network or Conda access. Reference is phiX174 (5,386 bp), k=31, three
chemistries simulated with distinct error/quality profiles, **never pooled**.

Script: `benchmarking/qualmer_joint_quality_calibration.jl` Data:
`benchmarking/results/qualmer_quality_channel_probe/qualmer_joint_quality_calibration_bins.tsv`,
`…_saturation.tsv`

### 3.1 Method

A k-mer's joint quality is a falsifiable **prediction**: `Q` implies
`P(this k-mer is wrong) = 10^(-Q/10)`. On simulated data ground truth is
knowable — a k-mer either does or does not occur in the reference — so:

```
observed error rate  = fraction of k-mers in a bin absent from the reference
predicted error rate = mean of 10^(-Q_joint/10) over that bin
```

Calibrated ⇒ observed ≈ predicted. Overconfident ⇒ observed ≫ predicted.

Binned by **observation count**, because the entire question is how confidence
should grow with repeated observation. Three confounds were deliberately
separated:

1. **Clamping vs model.** Joint quality is recomputed _unclamped_ from raw
   evidence, so "the 255 clamp saturated" is never mistaken for "the model is
   overconfident."
2. **Random vs correlated error.** Two regimes: random substitutions only; and
   the same reads plus a **systematic artifact** — one fixed reference locus
   mis-called the same way in every read covering it, at low reported quality.
   Its evidence count matches real sequence at the same depth; only its quality
   differs.
3. **Independence vs conservative.** Both aggregations computed for comparison.

### 3.2 Result 1 — independence fails exactly where the framework's claim lives

Under **random (independent) errors**, error k-mers essentially never recur.
Bins with ≥3 observations contain none, so there is no measurable miscalibration
— the assumption is untestable there because the failure mode does not arise.

Under a **correlated artifact**, the model becomes overconfident by margins that
**grow with observation count** — the signature of the independence assumption
failing:

| chemistry | coverage | obs bin | observed error rate | predicted (independence) | overconfident by |
| --------- | -------- | ------- | ------------------- | ------------------------ | ---------------- |
| illumina  | 40x      | 17+     | 5.9e-3              | 3.4e-42                  | **10^39**        |
| illumina  | 20x      | 17+     | 1.3e-2              | 9.6e-24                  | 10^21            |
| illumina  | 10x      | 9-16    | 8.3e-3              | 4.1e-14                  | 10^11            |
| pacbio    | 40x      | 17+     | 1.3e-2              | —                        | 10^12            |
| pacbio    | 20x      | 9-16    | 1.3e-2              | —                        | 10^5.7           |
| ont       | 20x      | 9-16    | 3.6e-2              | —                        | 10^3.8           |

Two structural observations:

- **Overconfidence grows with observation count** under correlation (illumina
  10x: 10^11 → 20x: 10^21 → 40x: 10^39), whereas it is flat (~10^1.2) in the
  random regime and confined to singletons. More observations of a correlated
  artifact make the model _more_ wrong, not less.
- **Severity is chemistry-ordered**: illumina ≫ pacbio ≫ ont, tracking
  per-observation Phred. The higher the reported quality, the faster summing
  overstates confidence. This must never be reported as a pooled number.

### 3.3 Result 2 — the 255 clamp saturates at ordinary Illumina depth

Fraction of **true** k-mers pinned at the clamp, where the score can no longer
discriminate at all:

| chemistry | 5x    | 10x   | 20x       | 40x       | median unclamped joint @40x  |
| --------- | ----- | ----- | --------- | --------- | ---------------------------- |
| illumina  | 13.9% | 62.0% | **97.6%** | **99.0%** | 1147 (4.5x over the ceiling) |
| pacbio    | 0%    | 3.6%  | 37.5%     | 66.3%     | 400                          |
| ont       | 0%    | 0%    | 0%        | 3.4%      | 126                          |

At ≥20x Illumina — an entirely ordinary depth — a joint-quality gate is
effectively a **coverage gate**, because almost every real k-mer sits at the
ceiling. Note this is an orthogonal problem to §3.2: it is a representation
limit (UInt8), not a modelling error, and it would persist under any aggregation
model that grows without bound.

### 3.4 What these results do NOT establish

- The error model is **synthetic**. The artifact regime is an idealised
  correlated error (one locus, every covering read, uniform low quality). Real
  correlated error is messier and its quality signature may be weaker or absent
  — instruments often do _not_ flag systematic artifacts with low Phred, which
  would make the picture worse, not better.
- No real sequencing data was used (no Conda root on the host, so
  ART/pbsim/badread were unavailable).
- The **conservative** model was computed but its calibration was not analysed
  in depth; whether it fixes the correlated case or merely dampens it is open.
- Real instrument Phred miscalibration is not modelled at all here —
  oracle-style quality was used, which is _more_ informative than reality.

---

## 4. Where this sits in the codebase

Relevant, all verified by reading (not inferred from names):

| Location                                                                  | What it does                                                                                                |
| ------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `core/quality-functions.jl:43` `combine_phred_scores`                     | The additive model. Sums, clamps to 255. Carries an assumption notice as of commit `d2886234`.              |
| `core/quality-functions.jl` `get_vertex_joint_quality`                    | Per-position joint quality; the production entry point.                                                     |
| `core/quality-functions.jl` `get_vertex_mean_quality`                     | Mean, not joint. **Cannot express compounding** — a k-mer seen 10x at Q40 scores the same as one seen once. |
| `core/quality-functions.jl` `filter_by_quality`, `get_vertex_min_quality` | The conservative per-position-minimum posture. **Zero production callers.**                                 |
| `rhizomorph/assembly.jl` `_qualmer_joint_confidence`                      | New opt-in traversal gate: joint quality across observations, weakest position across the k-mer.            |
| `viterbi-next.jl:10,847`                                                  | Viterbi emission model; genuinely consumes per-position Phred.                                              |
| `iterative-assembly.jl:6190`                                              | Corrector likelihood; accumulates Phred-derived log-probabilities.                                          |

**Immediate context for why this matters now.** A separate investigation (bead
td-4e19d family) established that the _default_ greedy qualmer assembly path is
entirely quality-blind — its traversal is topological and shared with the plain
k-mer arm, so the two produce byte-identical assemblies on every metric across
66 committed benchmark runs and 18/18 independent probe cells. An opt-in
`traversal_weighting = :quality` was added that makes quality load-bearing, and
it uses joint quality. So the independence assumption has just moved from
dormant to load-bearing, and the Rhizomorph manuscript credits Phred-derived
weights as a core construction.

---

## 5. The actual research questions

Ordered by how much they change the decision.

1. **How do established tools aggregate per-base quality across observations,
   and do any assume independence?** Specifically: GATK HaplotypeCaller,
   FreeBayes, DeepVariant, bcftools; the k-mer correctors BayesHammer, Lighter,
   Musket, Bloocoo; Merqury/Meryl QV estimation; polishers medaka, racon, Pilon.
   What is the _standard_ treatment, and where do they document departing from
   independence?

2. **How does the field handle known correlated error sources** — strand bias,
   tile/flowcell/cycle effects, PCR duplicates, homopolymer/context-specific
   error? Duplicate marking, strand-bias filters and BQSR all look like
   admissions that raw quality summation is unsafe; is that the right reading?

3. **Is there a principled middle model?** Candidates to evaluate: capped or
   discounted independence; effective-sample-size / variance-inflation shrinkage
   (n_eff < n); beta-binomial or Dirichlet-multinomial per-position posteriors;
   hierarchical models with a per-locus random effect; explicit
   strand/source-stratified aggregation requiring support from ≥2 independent
   sources. Which are actually used in production bioinformatics rather than
   merely available in theory?

4. **What is the right empirical calibration check on REAL data?** We have a
   working check on simulated data (§3.1). What is the field-standard analogue —
   e.g. BQSR-style empirical-vs-reported quality, or QV concordance against a
   truth set (GIAB, HG002) — and what would we need to run it?

5. **Representation.** Is a UInt8 joint-quality field defensible given §3.3, or
   should joint confidence be stored as a float / log-probability? What do
   comparable tools store?

6. **Reporting obligation.** If independence is retained, what must the
   manuscript state so the claim is defensible to a reviewer who knows artifacts
   are correlated?

---

## 6. What the deep-research AI should and should not do

**Should:**

- Answer §5 with **citations to primary sources** — tool documentation, methods
  papers, source code — not general reasoning about statistics.
- Distinguish sharply between what tools _do_ and what tools _claim_.
- Recommend a specific model with an explicit rationale, and say what evidence
  would change the recommendation.
- Flag anything in §3 that looks like an artifact of our synthetic setup rather
  than a real property of the estimator. Adversarial reading of our own evidence
  is wanted.

**Should not:**

- Propose changes to Julia code — the implementation is not the bottleneck; the
  decision is.
- Treat the §3 numbers as validated on real data. They are simulation.
- Assume the assembler's default path uses quality at all — it does not (§4).
- Recommend simply switching to the conservative model without addressing
  whether it actually fixes the correlated case, and at what cost in sensitivity
  to genuine low-frequency variation. The design's own §4.4 warns that
  aggressive correction removes real rare variants along with errors; a model
  that is too conservative is a real failure mode, not a safe default.

---

## 7. Files to upload with this document

| #   | File                                                                                                | Why                                                                                                                      |
| --- | --------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| 1   | this briefing                                                                                       | Start here; self-contained.                                                                                              |
| 2   | `planning-docs/rhizomorph-graph-ecosystem-plan.md` (§1.2.1 aggregation, §4.4 correction philosophy) | The primary design source; contains both models and the correction philosophy. Large — the two sections are what matter. |
| 3   | `benchmarking/results/qualmer_quality_channel_probe/qualmer_joint_quality_calibration_bins.tsv`     | Raw calibration numbers behind §3.2.                                                                                     |
| 4   | `benchmarking/results/qualmer_quality_channel_probe/qualmer_joint_quality_saturation.tsv`           | Raw saturation numbers behind §3.3.                                                                                      |
| 5   | `src/rhizomorph/core/quality-functions.jl`                                                          | The implementation, with the assumption notice in place.                                                                 |
| 6   | `benchmarking/qualmer_joint_quality_calibration.jl`                                                 | The experiment, if methodology is questioned.                                                                            |
| 7   | `benchmarking/results/qualmer_quality_channel_probe/README.md`                                      | Context on why the quality channel was dormant until now.                                                                |
