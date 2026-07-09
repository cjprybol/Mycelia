# The graph as a hidden Markov model: a probabilistic read corrector

- **Type:** Methods / design draft (manuscript-pipeline seed, not a final paper)
- **Date:** 2026-07-09
- **Scope:** `src/viterbi-next.jl`, `src/rhizomorph/algorithms/error-correction.jl`,
  the corrector tier routing in `src/rhizomorph/assembly.jl`, and the iterative
  corrector in `src/iterative-assembly.jl`.
- **Status:** Draft. Numbers below are drawn only from merged benchmarks and
  committed test fixtures; each is cited to its source. Several are SMOKE-scale
  toy results and are labelled as such.

## Abstract

We treat a de Bruijn / Rhizomorph assembly graph as a hidden Markov model in
which vertices are hidden states, coverage-weighted edges are state transitions,
and each read is an observation sequence emitted with a per-base, quality-aware
error model; read correction is then the maximum-likelihood decode (Viterbi) of
the read against the graph, replacing heuristic tip-clipping with a single
principled probabilistic operation. Because the objective is likelihood over the
graph rather than structural pruning, a supported-but-rare variant is retained by
construction — an emission-exempt beam and a coverage-tied support floor keep a
minority allele down to the most extreme skew our holdout fixtures assert
(100:3 ≈ 33:1) while error k-mers still decay — and because the recurrence is
alphabet-, strand-, and quality-generic, one corrector runs unchanged over DNA,
RNA, protein, and tokenized text. A sequence of loop-invariant / bound-the-work
fixes linearized a pipeline that previously exhausted memory on a 48 kb input,
cutting the dominant per-read decode term by roughly an order of magnitude at
the 20 kb scale we have profiled end-to-end.

## 1. Thesis: correction as maximum-likelihood decoding on a graph-HMM

### 1.1 The model

A k-mer (or qualmer) assembly graph is reinterpreted as a hidden Markov model:

- **States = vertices.** Each `(vertex, strand orientation)` is a hidden state.
  The true sequence that generated a read is a *path* through these states.
- **Transitions = evidence-weighted edges.** The transition probability from one
  state to a successor is the successor's normalized edge weight,
  `log(cov_edge / cov_total_outgoing)` — the coverage the data provide for taking
  that step (`_total_outgoing_weight`, `src/viterbi-next.jl`). This is the
  "how well does the graph support this continuation" term.
- **Observations = the read.** The read is the emitted symbol sequence. Decoding
  it means finding the state path most likely to have produced it.
- **Emission = a quality-aware, per-base error model.** For an observed unit
  against a candidate node, the emission log-probability sums per-position terms:
  a match contributes `log1p(-p_err)` and a mismatch `log(p_err / (|Σ| - 1))`,
  where `p_err` is derived from the read's **per-base Phred quality** at that
  position, not a population average (`default_viterbi_emission_logp`,
  `src/rhizomorph/algorithms/error-correction.jl`). Quality is threaded to the
  decoder by wrapping each k-mer with its read-relative Phred window
  (`QualityObservation`, restored in PR #359). A checkpoint test proves the
  decision is quality-driven: the same k-mer/node pair yields a different
  emission and a flipped correction (`TGA`→`TGC` under low quality, kept under
  high quality) on quality alone (PR #359).

### 1.2 Correction is the decode

Given the model, correction of a read is the **Viterbi maximum-likelihood path**:
the state sequence maximizing `Σ (emission + transition)` over the read. The
corrected read is the sequence spelled by that path. Where multiple candidate
paths compete (a genuine bubble), an EM-style soft assignment splits each read's
responsibility across the competing decoded paths and accumulates it onto edges
(`accumulate_competing_paths!` / `register_soft_edge_weights!`, PR #369), so that
error edges — traversed only by rare, low-probability paths — lose weight across
iterations and fall below the emergent-cleaning gate without any structural
tip-clipping step.

### 1.3 Why à-la-carte / principled, not heuristic

The repository also carries a greedy baseline
(`correct_sequence_greedy`, `src/rhizomorph/algorithms/error-correction.jl`):
for each position, pick the highest-coverage k-mer within Hamming distance 1 that
extends the previous choice. That is a local, first-order heuristic with no global
optimality guarantee and no quality model. The graph-HMM decoder replaces it with
a global objective:

- **Exact by default, no heuristic fallback.** `find_optimal_sequence_path`
  returns the Viterbi ML result or reports "no improvement" — it never falls
  through to a statistical-resampling / local-edit cascade that could land a
  "corrected" read on a path corresponding to no real graph traversal (removed in
  PR #359). The beam is `typemax(Int)` (frontier never pruned) unless a caller
  opts into a finite width as an explicit speed/accuracy trade.
- **À-la-carte emission.** The emission callback is a documented extension seam:
  the same decode machinery accepts a swap-in alphabet- or quality-aware model
  without touching the dynamic-programming core (`ViterbiCorrectionConfig`).
- The result is that every correction is a defensible statement — "this is the
  most likely sequence given the read's own base qualities and the coverage
  evidence in the graph" — rather than a structural edit justified only by local
  topology.

## 2. Variation preservation: the scientific differentiator

### 2.1 The invariant

The corrector must **remove error while retaining real variation**: a true
heterozygous site, paralog, or skewed minor allele stays; only unsupported error
is removed. This is the property that separates a probabilistic corrector from a
pruner — a tip-clipper cannot tell a low-coverage error from a low-coverage real
allele, so it risks erasing biology.

Two mechanisms make retention structural rather than incidental.

### 2.2 The support floor (soft-EM)

The soft-EM M-step, left unguarded, mathematically **collapses** a skewed minor
allele: when a strictly heavier sibling exists, the minority read shares
responsibility with the majority alternative, and the recurrence
`W_min' = N·W / (W_maj + W_min)` contracts the minority geometrically — the same
form as an error, with balanced 50/50 as the only stable fixed point (PR #369,
bead td-h6w9). The fix clamps each registered edge weight to
`max(responsibility_weighted_value, floor)`, where `floor == raw_coverage` for an
edge backed by at least `SOFT_EM_MIN_SUPPORT = 3` reads and 0 otherwise. A
well-supported edge therefore **never decays below its own coverage** regardless
of a heavier sibling; a coverage-1 error edge keeps a floor of 0 and is free to
decay. Decay is *unsupported*-based, not *less-than-sibling*-based.

### 2.3 The emission exemption (beam pruning)

The scaling work introduced a score-margin ("histogram") beam that drops states
more than Δ below the depth's best. Naively this compares the **full** cumulative
score, whose transition term `log(cov_edge/cov_total)` penalizes a path for being
*rare*, not *wrong* — so in a skewed pool a real minor allele could in principle
be pruned (PR #388 review flagged this as ship-blocking for the viral-
quasispecies domain). The fix is an **emission-exemption AND-gate**: a frontier
state is pruned only when it is more than Δ below the best on **both** the full
score **and** the cumulative emission (read-consistency, excluding the transition
term) (`beam_score_margin`, `src/viterbi-next.jl`; PR #388). Consequences:

- A read-consistent path (good emission) is never dropped for coverage rarity —
  supported variation is protected. A single-base mismatch costs roughly 5.7 nats
  of emission per spanning k-mer (~120 nats over k = 21), which dwarfs the
  coverage penalty (bounded by `log(skew)`, ~2–6 nats), so the margin prunes
  *wrong* paths, never merely *rare* ones.
- Error correction is not weakened: an uncorrected error path also has high
  emission → exempt → still available for the full-score ML choice, which selects
  the corrected path. (An emission-*only* margin, the reviewer's literal
  suggestion, would instead keep the higher-emission uncorrected path — the
  AND-gate is the form that avoids that regression.)

### 2.4 Structural cleanup respects the invariant

The linear-time graph-cleanup passes that run before contig extraction
(`clip_error_tips!`, `collapse_error_bubbles!`, `prune_disconnected_error_components!`,
PRs #384/#385) are all gated to remove only *unambiguous* error: a tip is clipped
only if it is a coverage-1 dead-end off a junction (a real variant rejoins as a
bubble, never a dead end); a bubble branch is collapsed only when its mean
coverage is ~1 *and* the sibling is well supported; a disconnected component is
pruned only when it is separate, small, uniformly low-coverage, *and* shares no
canonical k-mer with any well-supported vertex. Each is deliberately designed to
**under-prune** rather than risk a real low-coverage variant (the td-h6w9
invariant).

### 2.5 Quasispecies safety (evidence and honest bound)

The committed variation-preservation holdout asserts retention of both balanced
and skewed alleles alongside removal of a coverage-1 error, across the pipeline
(24/24 assertions green as of PR #388) and in a dense-graph regime where the beam
and score margin are actually engaged (16/16). The most extreme skew the
fixtures assert is **100× / 3× ≈ 33:1** (`_vph_build_dense_skewed_fixture`,
`test/4_assembly/variation_preservation_holdout_test.jl`, cases
`(15,15)`, `(20,4)`, `(100,3)` at k = 21 and `(20,4)` at k = 9). At 100:3 the
minority allele is retained with 100% recall in the recorded run (PR #388). We
note explicitly: a higher skew (e.g. 67:1) has **not** been validated in a
committed fixture; the defensible claim today is retention to ~33:1, with the
mechanism (emission dominance over the `log(skew)` coverage penalty) predicting
that the safe skew scales with k.

## 3. Universality: one corrector across modalities

The decode is generic along three independent axes, so a single implementation
serves every modality:

- **Alphabet-generic.** The emission model resolves and validates the alphabet
  (`:dna`, `:rna`, `:protein`, `:text`) and uses `|Σ|` in the mismatch term; the
  graph carries opaque labels. DNA/RNA/protein k-mers and SentencePiece text
  tokens all flow through the same `correct_observations` (PR #382).
- **Strand-generic.** `strand_mode` selects `:singlestrand`, `:doublestrand`,
  `:canonical`, or `:auto`. Reverse-complement-defined modes are validated out
  for alphabets where RC is undefined (amino acids and general strings are
  single-strand only; `AssemblyConfig` throws otherwise).
- **Quality-optional.** Emission consumes per-base Phred when present (FASTQ) and
  falls back to a uniform `error_rate` otherwise (FASTA), with no change to the
  recurrence.

Committed modality validation (PR #382):

| Modality | Fixture | Result |
| --- | --- | --- |
| dsDNA | primary path, all assembly benchmarks | see §5 |
| ssRNA | ~300 nt single-strand molecule, simulated substitutions, `mode=:singlestrand`, RNA k-mers | substitutions recovered |
| Protein (AA) | RC-undefined guard | `:doublestrand`/`:canonical` correctly throw; `:singlestrand` path exercised |
| Text | SentencePiece → token graph (`StringVertexData`/`StringEdgeData`) → Viterbi | one corrupted token restored |

Honest caveat on text: the TEXT emission indexes token strings by byte and
throws on multibyte UTF-8 pieces, including the U+2581 word-boundary marker real
SentencePiece emits; the gated test strips U+2581 as a workaround, and a
character-aware emission fix is deferred to the correction-core track (PR #382).

## 4. Scaling: bound-the-work, applied repeatedly

### 4.1 The tiered strategy

The corrector exposes two tiers (`_corrector_strategy_knobs`,
`src/rhizomorph/assembly.jl`):

- **`:scalable` (default)** — a coarse LoRMA-style 3-rung k-ladder, a low (2)
  iteration cap, skip-solid volume reduction, a Stage-0 linear k-mer-spectrum
  pre-correction, hard-read gating, soft-EM edge memory, the size-aware auto-beam,
  and a `:doublestrand` graph. Built for real-scale inputs.
- **`:exhaustive`** — maximum-sensitivity exact-ML: prime-by-prime k-walk, 10
  iterations/k, an unbounded (`typemax(Int)`) Viterbi beam, no skip / no
  hard-window / no soft-EM. Byte-identical to the prior exact decoder, and able
  to OOM on very large reads — intended for small, high-sensitivity inputs.

The exactness knob is the beam width: with `beam_width = typemax(Int)` the
accelerated `:scalable` path is byte-identical to the exact decoder; the tiered
knobs (`beam_score_margin`, `max_successors_per_state`) default to no-ops and
engage only where the width beam is already finite, so exact-ML reads never
change.

### 4.2 The O(genome²) → O(genome) story

The same anti-pattern — an unbounded or linearly-scanned inner structure inside a
per-read/per-position loop, giving a per-pass cost that grows with genome size —
recurred at several sites and was fixed the same way each time (bound the work /
restore a loop invariant):

| Fix | Anti-pattern removed | Effect | Source |
| --- | --- | --- | --- |
| Beam cap on the decode frontier | frontier grew ~unboundedly with read depth; 21-billion-allocation OOM crash on a 48 kb phage | made iterative correction runnable at all | PR #344 |
| Convergence stop + coarse 3-rung k-ladder | 10 k-steps × 10 iters walked every prime | 4.84× faster, byte-identical (100-read/1 kb toy: 126.1 s → 26.0 s) | PR #352 |
| Bubble/out-adjacency index | O(V·E) bubble detection | O(V+E) | PR #373 |
| Decode-hoist: build weighted graph / out-weight / candidates once per pass | rebuilt per read; O(E) edge scan per state; O(V) candidate scan per read | decode α 2.15 → 1.04; ~26× at 5 kb (210 s → 1.1 s/pass); byte-identical | PR #374 |
| Density-aware auto-beam | short reads decoded with an unbounded frontier on dense mid-k graphs | mid-rung frontier α 1.878 → 0 (capped at 256) | PR #379 |
| Set-membership in `find_linear_path` | O(L) vector scan per step → O(L²) ≈ O(genome²) contig walk | walk slope ~2.0 → ~1.1, ~100× at 32 k-vertex contigs; byte-identical | PR #383 |
| Re-assembly graph reuse | rebuilt the qualmer graph from scratch after correction | iterative-arm wall −16–26%, byte-identical | PR #377 |
| Score-margin (emission-exempt) decode frontier | ~99% read-inconsistent states still generated successors each depth | decode 391.7 s → 18.3 s (21×), wall 472 s → 54.7 s (8.6×) at 20 kb; decode α 1.62 → 1.33 | PR #388 |

The recurring signature: the profiler localizes the super-linear term to one
component, and the fix bounds that component's per-item work so the per-pass cost
returns to linear in genome size. The 48 kb input that once OOM-crashed
(21-billion allocations, before any beam; PR #344) is the anchor of the arc;
after the decode-bounding work the dominant per-read decode term at the 20 kb
scale we profile end-to-end drops from 472 s to 54.7 s (PR #388).

**Honest caveat on the "48 kb → 4.7 min" framing.** The 48 kb OOM crash is
verified (PR #344). A *post-fix* 48 kb wall-clock of ~4.7 min is **not** present
in a committed benchmark: the merged end-to-end profile
(`benchmarking/results/empirical_48k_component_profile.txt`) tops out at 20 kb
(455 s pre-#388), and the #388 harness reports the 20 kb wall at 54.7 s post-fix.
The closest defensible statements are the two just given; a direct 48 kb
end-to-end re-measurement after the full scaling campaign is pending (see §6).

## 5. Results (verified numbers only)

All numbers below are toy / SMOKE-scale simulated data unless noted; each is
cited. They are *engineering* validation of the pipeline, not a benchmarking
study against external assemblers.

### 5.1 Naive vs. iterative correction (assembly quality)

Simulated Illumina, k = 21, 20×, 150 bp reads, err = 0.01 (SMOKE toy):

| Genome | Config | Contigs | N50 | Frac | Source |
| --- | --- | --- | --- | --- | --- |
| 1 kb | naive (no correction) | 458 | 41 | — | `quality_gap_diagnostic.csv` |
| 1 kb | iterative, canonical (mis-configured) | 197 | 43 | — | `quality_gap_diagnostic.csv` |
| 1 kb | iterative, doublestrand (fixed) | 16 | 891 | — | `quality_gap_diagnostic.csv` / PR #368 |
| 1 kb | iterative + RC-dedup | 6 | 904 | 1.07 | PR #380 |
| 1 kb | iterative + defrag + component prune | 1 | 974 | 0.97 | PR #385 |
| 3 kb | iterative + RC-dedup | 16 | 2926 | 1.15 | PR #380 |
| 3 kb | iterative + defrag + component prune | 1 | 2957 | 0.99 | PR #385 |

The N50 and largest-contig are identical across the defrag/prune steps (the main
genome is untouched); the contig-count drop and the fraction settling from >1.0
toward ~0.98 reflect removal of spurious error-debris islands, not real sequence
(PR #385). The single dominant configuration lever was `graph_mode`: flipping
only `:canonical → :doublestrand` on the 1 kb fixture collapsed 197 → 16 contigs
and lifted N50 43 → 891 (PR #367/#368).

### 5.2 Scaling (end-to-end, real branchy corrected graph)

`assemble_genome(corrector=:iterative, strategy=:scalable, k=21)`, err = 0.01,
20×, 150 bp, k-ladder [3, 9, 21]
(`benchmarking/results/empirical_48k_component_profile.txt`, pre-#388):

| Genome | Wall | Contigs | Max contig | Per-read decode share |
| --- | --- | --- | --- | --- |
| 5 kb | 59.5 s | 15 | 2807 bp | ~80% |
| 10 kb | 151.0 s | 36 | 5355 bp | ~80% |
| 20 kb | 455.1 s | 71 | 13529 bp | 82% |

Total-wall α ≈ 1.47 pre-#388, with the residual localized to the low/mid-k dense
decode. Post-#388 (decode score-margin), the 20 kb wall drops 472 s → 54.7 s
(8.6×) and the decode term 391.7 s → 18.3 s (21×) on the #388 harness.

### 5.3 Recall / accuracy (toy)

Post-#388, on toy fixtures: 1 kb k-mer recall 92.04% (identical to pre-fix), 3 kb
98.36% vs 97.52% (higher, with fewer spurious k-mers); the 20 kb assembly is
cleaner (25 vs 44 contigs, longer max contig) (PR #388).

### 5.4 Variation holdout

Balanced + skewed pipeline holdout 24/24; dense-regime skewed rare-allele holdout
16/16 across balanced (15/15), mild (20/4) and extreme (100/3 ≈ 33:1) skews at
k = 21, plus 20/4 at k = 9; the 100:3 minority is retained with 100% recall
(PR #388, `variation_preservation_holdout_test.jl`).

### 5.5 Stage-0 classifier (k-mer solidity)

Ground-truth grid (3 seeds × 3 error rates, reflen 10 kb;
`benchmarking/results/scale_calibrate_findings.md`): coverage-based arms saturate
at AUC ≈ 1.0; balanced accuracy MixtureModel 0.998 ± 0.001, BayesianMixture
0.997 ± 0.001; quality-only is honestly weaker (AUC ≈ 0.81) — per-base quality
alone separates genomic from error k-mers far less cleanly than coverage. Stage-0
cheap correction cuts the graph-Viterbi decode fraction from ~80% to ~25% per
pass (PR #364).

### 5.6 Modality coverage

See §3 table. dsDNA is the primary, benchmarked path; ssRNA / protein-guard /
text are functional-validation fixtures (PR #382).

## 6. Limitations and future work

- **Soft-EM speed.** The competing-paths E/M step (`C5c`) is a small share of the
  wall today (≈2% at 20 kb) but has the steepest large-end slope (α_lg ≈ 1.9,
  `empirical_48k_component_profile.txt`); it is the next super-linear term to
  bound if scale increases.
- **Pruner debris.** After defrag + component pruning, a residual of disconnected
  coverage-1/2 error islands can survive when they re-use real sequence content
  (the real-sequence guard deliberately retains them). The design under-prunes;
  a tighter but still variant-safe island test is open.
- **GPU/SIMD acceleration is proposed, not implemented.** The array-frontier
  reformulation (depth-outer, read-inner, dense state array) and the two-phase
  acceleration (CPU SIMD-POA on hard windows; a GPU warp-shuffle wavefront behind
  a KernelAbstractions extension) are specified in the merged ADR
  (`docs/design/2026-07-06-gpu-simd-corrector-acceleration.md`) with a
  byte-identical equivalence oracle as the acceptance gate. A CPU batched-decode
  proof-of-concept and a Phase-B GPU frontier kernel exist (PRs #362, #365) but
  are not on the default path.
- **Real-data validation pending.** All quality/scaling numbers above are
  simulated toy/SMOKE fixtures. The committed real-genome CSVs
  (`real_genome_benchmark_*.csv`) exercise the *plain* (`corrector=:none`)
  assembler, not the iterative corrector, and are not presented here as corrector
  output. A real-genome iterative-corrector benchmark (and a direct 48 kb
  post-scaling re-measurement, §4.2) is the primary open validation task before
  any manuscript claim.
- **Quasispecies skew bound.** Committed retention is validated to ~33:1
  (§2.5); higher skews (viral quasispecies routinely exceed this) need explicit
  fixtures before the safe-skew claim can be widened.
- **Text emission UTF-8.** Multibyte pieces (U+2581) throw in the byte-indexed
  TEXT emission; a character-aware fix is deferred (§3).

## Provenance

Every quantitative claim traces to one of: a merged PR commit body (cited inline
as PR #NNN), a committed benchmark result under `benchmarking/results/`, a
committed test fixture under `test/4_assembly/`, or the merged GPU/SIMD ADR. This
draft touches no `src/` or test code; it is a synthesis of already-merged work
intended to seed a Stanford/LBNL methods manuscript. The next manuscript step is
the real-data validation in §6, after which §5 can be re-cast from engineering
SMOKE fixtures to a reportable benchmark against external correctors.
