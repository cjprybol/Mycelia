# ADR: Restore the graph-as-HMM probabilistic correction core

- Status: proposed
- Date: 2026-07-06
- Scope: Rhizomorph iterative assembly / read correction (`src/iterative-assembly.jl`,
  `src/viterbi-next.jl`, `src/rhizomorph/assembly.jl`, `src/rhizomorph/core/graph-construction.jl`)

## Context — the design thesis

Rhizomorph's assembler is meant to be a **principled probabilistic method**, not a
heuristic bag of tricks. The foundational thesis:

> The assembly graph **is** a Hidden Markov Model. Graph vertices (k-mers /
> qualmers / tokens) are HMM states; edges are transitions whose probabilities come
> from evidence (coverage, quality). A read is an observation sequence. Correcting a
> read = threading it through the graph and taking the **maximum-likelihood state
> path (Viterbi)**, where the **emission probability is the per-read per-base Phred
> quality** — "how much do I trust this read's base here" vs "how much does the graph
> support the alternative." The probabilistic variant resamples paths and updates the
> graph's edge/node weights from those path probabilities (**soft EM**): reads
> re-thread against the updated graph, iterating until the graph **coalesces to its
> most-likely form**. Error edges lose support and disappear on their own — no
> heuristic tip-clipping required, because cleaning is an *emergent* property of the
> EM, not a bolted-on threshold.

This ADR records that the current implementation has the correct skeleton but has
been undermined by (a) a dead emission model, (b) accuracy-sacrificing speedups, and
(c) contig extraction from an uncleaned graph in a partly-broken graph mode. It sets
the corrected design and a restructured work plan.

## Evidence — what the audits found (2026-07-06)

### Correct (the bones are real)
- The graph is the state space: states are `(vertex_label, StrandOrientation)`
  tuples; reads are threaded through graph edges (`viterbi-next.jl:925, 985`).
- The DP is a real Viterbi, not greedy — full frontier + backtracking over
  `log(transition) + log(emission)` (`viterbi-next.jl:975-1066`).
- Transitions come from edge evidence, normalized over out-edges
  (`viterbi-next.jl:1000-1005`, `path-finding.jl:725-730`).
- An EM-like loop exists (rebuild graph from corrected reads, re-thread).

### Broken / undermining the thesis
1. **Dead emission model (central defect).** Per-read per-base Phred quality is
   `collect`ed at `iterative-assembly.jl:1332` and then **never passed to the
   Viterbi**; observations are plain `DNAKmer` objects with no quality
   (`iterative-assembly.jl:1327`). Emission falls back to the graph vertex's
   **population-average** quality for in-graph k-mers, or a uniform
   `error_rate=0.01` for error k-mers (`viterbi-next.jl:558-575`). The emission
   therefore measures the graph's confidence in a k-mer, not the read's confidence
   in its base — `P(observed | hidden)` has lost the *observed*.
2. **Quality stripping.** FASTA inputs get a uniform Q40 placeholder
   (`assembly.jl:599`, `_placeholder_qual(n) = repeat("I", n)`), and even for FASTQ
   inputs the real quality is discarded before the emission (defect #1). The
   "quality-aware" qualmer path runs on fake or ignored quality.
3. **Beam pruning drops the ML guarantee.** `beam_width=256` is hardcoded in the
   iterative path (`iterative-assembly.jl:1342`); the exact default (`typemax`) is
   overridden, so the returned path may not be the global ML path.
4. **Heuristic fallback cascade.** `find_optimal_sequence_path`
   (`iterative-assembly.jl:673-698`) tries Viterbi, then falls back to statistical
   resampling, then to local heuristics — so a "corrected" read may not lie on any
   ML graph path.
5. **Hard-assignment EM, not soft.** The M-step rebuilds the graph from corrected
   *sequences* rather than updating edge weights from path *probabilities*; the
   graph keeps no path-probability memory.
6. **No graph cleaning + broken reconstruction mode.** Contigs are extracted from an
   *uncleaned* graph (`bubble_resolution`/`repeat_resolution` are no-op log stubs;
   `min_coverage` is a dead field), so every sequencing error breaks a unitig — a
   1 kb region yields ~240 contigs vs the ~1 expected from Lander–Waterman. And
   Canonical-mode contig reconstruction is broken for both k-mer and qualmer graphs
   (undirected traversal unimplemented) — yet `skip_solid` forced the corrector into
   Canonical.
7. **Matrix gaps.** Only {k-mer, qualmer} × {SingleStrand, DoubleStrand} + n-gram
   SingleStrand produce valid assemblies. Canonical (2 cells) is broken; the
   token/SentencePiece node type is not wired into `assemble_genome`.

## Empirical consequence

A corrected-then-reassembled toy validation (assembly-vs-assembly, err=0.01) showed
correction does **not** improve the assembly (genome fraction 63.3% naive → 61.8%
iterative) and slightly hurts it — consistent with a correction that runs on a dead
emission model against an uncleaned/partly-broken graph. The earlier "63→97.5%"
result was retracted (it compared naive assembly vs raw corrected reads).

## Decision

Restore the thesis rather than paper over it with more heuristics.

1. **Emission = per-read per-base quality.** Wrap each observation as
   `kmer + quality_scores[pos:pos+k-1]` (a `QualityObservation`) so
   `_viterbi_direct_quality_scores` returns the read's real per-base Phred, not the
   graph-vertex population average.
2. **Preserve quality end-to-end.** Remove the FASTQ→FASTA→Q40 stripping; real
   quality must reach the emission. FASTA inputs are explicitly flagged as
   quality-absent (not silently Q40).
3. **Trustworthy Viterbi.** Exact by default (or a principled, wide, *documented*
   beam only as an opt-in scaling knob); remove the heuristic fallback cascade from
   the correctness path — trust the ML path or report no improvement.
4. **Soft EM.** Update edge/node weights from path probabilities (M-step) so the
   graph coalesces probabilistically; cleaning **emerges** (low-support error edges
   die). Heuristic tip/bubble removal is retained only as a *comparison baseline*,
   never the primary mechanism.
5. **Graph mode.** Run correction + reconstruction on DoubleStrand (Canonical is
   broken); fix Canonical traversal as a separate, optional track if the
   strand-collapsed model is wanted.
6. **à la carte, honestly compared.** Once the core is correct, compare across the
   *working* matrix cells ({k-mer, qualmer} × {Single, Double}) with real metrics;
   wire the token graph and fix Canonical as explicit follow-ons.

Speedups (`skip_solid`, beam pruning, iteration caps, coarse k-ladder) are **opt-in
scaling knobs measured against the exact baseline**, never defaults that silently
change correctness. This reverses the prior stance where they were on by default.

## Consequences

- The correction core is re-architected; `try_viterbi_path_improvement` and the
  emission plumbing change materially. Existing "corrector improves things" claims
  are unproven until re-validated on the corrected core.
- `skip_solid`, beam pruning, and the k-ladder/iteration tuning (td-q70n) are
  demoted from defaults to measured opt-ins.
- Validation must be assembly-vs-assembly with graph cleaning emergent, swept over
  error rate (correction can only help where there is error to correct) and read
  regime (short-low-error vs long-high-error), since all-solid-read skips vanish for
  long reads.

## Restructured epic (replaces the per-read-decode framing under td-9q84)

1. **[core] Emission model = per-read per-base quality** (the central fix).
2. **[core] Quality flow**: stop FASTQ→FASTA→Q40; preserve real quality to emission.
3. **[core] Trustworthy Viterbi**: exact default; remove heuristic fallback from the
   correctness path; beam becomes an opt-in scaling knob.
4. **[core] Soft EM**: edge-weight updates from path probabilities; validate that
   cleaning emerges (error edges die) without heuristic tip-clipping.
5. **[validate] Assembly-vs-assembly, error-rate × regime sweep** on the corrected
   core (DoubleStrand), with a scale-assertion guard so toy runs can't be quoted as
   validation.
6. **[matrix] Fix Canonical reconstruction** (orientation-aware undirected
   traversal) — unblocks 2 cells and the strand-collapsed model.
7. **[matrix] Wire the token/SentencePiece node type** into `assemble_genome`.
8. **[baseline] Heuristic graph cleaning** (`remove_tips!` / bubble pop) as a
   *comparison arm only*, to quantify how close emergent-EM cleaning gets to it.
9. **[manuscript] Document** the graph-as-HMM/EM method and the corrected results.

Deprecated framing: "per-read decode is the bottleneck, optimize it" and "heuristic
graph cleaning is the fix" — both were shortcuts around the real defect (dead
emission + untrustworthy Viterbi).
