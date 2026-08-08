# opt4: Opt-in approximate frozen-read skip (with dnadiff accuracy sign-off)

Date: 2026-07-26 Status: Implemented (Rule-of-5 passes 1–5 landed; see
`benchmarking/results/opt4_frozen_read_skip_signoff.md` for measured results)
Scope: skip re-decoding converged ("frozen") reads across corrector iterations —
an OPT-IN, APPROXIMATE tier, default off, gated by an accuracy sign-off.
**NOT byte-identical.**

> **Implementation update (pass 4/5):** the toy-genome positive control below
> (methodology §3) was found **mechanically inert** at pass 4 — Stage 0 cheap
> correction runs before the freeze gate and dominates toy-scale correction, and
> the Viterbi decode (the only thing freezing gates) makes no edits at toy scale.
> The working positive control therefore moved to **representative scale**
> (Lambda phage 21x) in pass 5, where `freeze_max_aggressive` degrades recall
> −1.0% to −4.8%. See the sign-off doc; the §3 toy methodology below is retained
> for provenance but is superseded by the representative-scale control.

## Context

This is step 5 of the corrector performance campaign. Steps opt5 (GC gating),
opt1 (parallel decode, ~2x measured), and opt2 (typed transitions) all shipped
and were **byte-identical** — they transformed the code that produces the
output without changing the output. opt3 (batched/SIMD kernel) is deferred (its
kernel is single-strand-only while production is doublestrand/canonical).

opt4 is categorically different: it is the campaign's first **approximate**
optimization. It trades a bounded, measured amount of accuracy for skipping
redundant work, and therefore **cannot** be byte-identical. The reason is
structural: the corrector rebuilds the qualmer graph from the *entire* current
read set on **every iteration** (`iterative-assembly.jl:2317`), so a read that
has "converged" (stopped changing) still influences the graph that every *other*
read decodes against on subsequent passes. Skipping a converged read's re-decode
keeps its last correction fixed even though the graph it would decode against has
changed — the changed graph might have corrected it differently. That divergence
from the exact path is the approximation.

Because it is approximate, the deliverable's foundation is not an identity lock
but an **accuracy sign-off**: a measurement showing that, on representative data,
the approximate tier's accuracy loss is within an acceptable tolerance for the
speedup gained — plus a **positive control** proving the measurement can detect
degradation when it occurs. The tier ships **opt-in, default off**.

## Key facts (verified against current master, incl. opt5+opt1+opt2)

File: `src/iterative-assembly.jl` unless noted. Line numbers verified 2026-07-26
against merge `b032fab3`.

1. **The loop nesting is k-rung → iteration → graph rebuild → batch → per-read
   decode.** Outer k-rung loop `while ... && k <= max_k` (`:2254`); inner
   per-k iteration loop `while iteration <= max_iterations_per_k` (`:2290`);
   the graph is rebuilt inside the iteration loop via
   `build_qualmer_graph(current_reads, k; mode = graph_mode, ...)` (`:2317`);
   the per-read decode runs in `_improve_read_set_likelihood_impl` (defined
   `:4555`, reached via the public wrapper `improve_read_set_likelihood` called
   at `:2408` which forwards to the impl at `:4513`); batch loop
   `for batch_start in 1:batch_size:total_reads` (`:4879`),
   per-read decode parallel `Threads.@threads for i in eachindex(batch_reads)`
   (`:4902`) / serial (`:4984`).

2. **Graph rebuilt EVERY iteration, not just every rung** (`:2317`, inside the
   iteration `while`). This is the crux that makes opt4 approximate: a frozen
   read still contributes its (fixed) corrected sequence to every subsequent
   rebuild, and every other read's decode target changes pass-to-pass.

3. **A per-read change signal already exists, but only per-pass.**
   `improve_read_likelihood` returns `was_improved::Bool` per read
   (`:4924`, `:4954`, `:4961`, `:4994`, `:5017`), summed into the pass counter
   `improvements_made`. It is **never persisted per-read across iterations** —
   there is no `converged`/`frozen`/`stable` per-read state today. opt4 must add
   a per-read consecutive-no-improvement streak counter.

4. **Existing convergence detection is aggregate/pass-level only.**
   `stop_on_no_change` breaks the k-loop on a zero-improvement pass (`:2503`,
   `:2585`); `sufficient_improvements` (`:5824`) does aggregate plateau
   detection. None of these is per-read.

5. **The skip mechanism is already proven and separable from decode.** The
   `skip_solid`/`hard_window` path uses `_gate_skip` (`:4647`), `base_skip_flags`
   (`:4660`), and `_skip_this_read_at` (`:4804`) to bypass a read's decode for a
   pass. A skipped read flows into `updated_reads[...] = read` unchanged
   (`:4908`/`:4986`), which is written to the next pass's FASTQ (`:2483`) and
   consumed by the next `build_qualmer_graph` (`:2317`). So a frozen read's last
   correction is *mechanically guaranteed* to reach the graph rebuild even with
   its decode skipped. opt4 reuses this exact path; only the "should I skip"
   decision is new.

6. **Opt-in kwarg precedent (opt5, shallow).** opt5's `gc_between_batches`
   kwarg threads: entry on `mycelia_iterative_assemble` (`:1739`) → forwarded at
   the `improve_read_set_likelihood` call (`:2412`) → kwarg on
   `improve_read_set_likelihood` (`:4489`) → forwarded to
   `_improve_read_set_likelihood_impl` (`:4519`) → kwarg on impl (`:4559`). It
   does NOT thread up into `AssemblyConfig`/`assemble_genome`. opt4 follows this
   same shallow path (kept minimal + directly benchmarkable).

7. **Accuracy tooling already exists.** `run_dnadiff` wraps MUMmer `dnadiff` via
   a bioconda env (`src/alignments-and-mapping.jl:548`). The per-base accuracy
   harness `benchmarking/rhizomorph_correction_accuracy_benchmark.jl` calls
   `mycelia_iterative_assemble` directly, reads `final_fastq_file`, and computes
   recall/precision/over-correction against injected-error ground truth. It has a
   read-scramble null (Control B) that destroys the mechanism and confirms the
   metric collapses — the same positive-control pattern opt4 needs.

## Design

### The freeze criterion

A read **freezes** when its `was_improved` has been `false` for **N consecutive
passes** (N = `freeze_streak_threshold`, a tunable; N≈2 is the proposed
default). This requires one new piece of cross-pass state: a per-read
`Vector{Int}` streak counter, threaded across the `improve_read_set_likelihood`
call boundary the same way `prev_soft_weights` already is (`:2287`/`:2447`).
Each pass: if `was_improved`, reset the read's streak to 0; else increment; a
read with `streak >= N` is added to the skip set for the *next* pass.

### Configurable scope (default within-rung) — the risk control

The dominant accuracy risk is freezing a read across a **k-change**: a read that
looks converged at low k (dense, ambiguous graph, lowest per-k-mer coverage) may
still need correction once the graph sharpens at a higher k. Two scopes, behind
a `freeze_across_rungs::Bool = false` flag:

- **Within-rung (default, safe):** reset every read's streak to 0 whenever k
  advances (`_next_k_in_progression`, `:2612`). A read is only ever skipped
  across same-k iterations; at every new rung — the biggest graph change — every
  read decodes fresh. This bounds the approximation to intra-rung refinements.
- **Across-rung (opt-in, aggressive):** persist the streak across k-changes for
  maximum decodes skipped. Larger speedup, higher accuracy risk. Only for users
  who have run the sign-off on their data and accepted the tradeoff.

### The skip hook

A frozen read is added to the `_skip_this_read_at` predicate (`:4804`), reusing
the `skip_solid` path (fact 5). No new decode-bypass plumbing is needed — the
frozen read's last correction already flows to `updated_reads` → next FASTQ →
next graph. The only additions are (a) the streak state and (b) OR-ing the
frozen predicate into the existing skip decision.

### Opt-in kwarg

Follow opt5's shallow pattern. Add `skip_frozen_reads::Bool = false` (and
`freeze_streak_threshold::Int = 2`, `freeze_across_rungs::Bool = false`) to
`mycelia_iterative_assemble` (`:1739`) and thread through the 3 hops
(`:2412` → `:4489`/`:4519` → `:4559`). Default off ⇒ zero behavior change for
every existing caller (the exact path is preserved bit-for-bit when the flag is
off — this IS byte-identical-when-disabled, and should be locked by a test).

### Interaction with existing skip/soft-EM/diagnostics

The frozen-skip must compose with, not replace, the existing per-pass
`skip_solid`/`hard_window` gating; soft-EM (`soft_weights`/`soft_weights_sink`)
accumulation; and `CorrectorDiagnostics`. A frozen read that is skipped
contributes **nothing** to that pass's soft-EM E-step — the skipped read never
reaches the `soft_weights` accumulator (a fresh per-pass accumulator), so it
contributes no weight rather than a stale/replayed one, matching pre-existing
`skip_solid`/`hard_window` skip behavior. A new diagnostics counter
`frozen_reads_skipped` records the volume of skipped decodes (mirrors opt1's
`parallel_decode_batches` counter, incl. the `fieldcount(CorrectorDiagnostics)`
completeness guard + its `corrector_errors` export site — a struct field has
three sites, per the opt1 incident).

## Accuracy sign-off methodology (replaces the byte-identity lock)

Three complementary checks, run skip-ON vs skip-OFF:

1. **Per-base read accuracy** (primary): extend
   `rhizomorph_correction_accuracy_benchmark.jl` to run with
   `skip_frozen_reads=false` (exact) and `=true` (approximate, both scopes), and
   report the delta in recall / precision / over-correction. This measures the
   corrector's own read output directly.
2. **Assembly-level identity**: re-assemble the corrected reads
   (`assemble_genome(corrector=:iterative)` → `_assemble_with_iterative_corrector`,
   `assembly.jl:1591`) and run `run_dnadiff` on the contigs vs the reference
   genome (skip-ON vs skip-OFF vs reference), reporting ANI / aligned-fraction /
   indel+SNP counts. This catches degradation that per-base read metrics might
   miss after assembly.
3. **Positive control** (the load-bearing gate): a toy genome (~2 kb, reusing
   the `MYCELIA_RCA_SMOKE` path) with a **repeat region carrying a variant that
   only resolves at a higher k rung** — a read genuinely wrong/unchanged at low
   k that flips correct at a later k. Run exact vs freeze (N=1, across-rung, to
   force early freezing). The freeze arm's accuracy metrics **must measurably
   degrade** on this engineered case. If the metric does NOT move, the
   measurement is not trustworthy and the sign-off is invalid — mirrors the
   existing Control B read-scramble null.

**Ship gate:** the tier ships opt-in **only if** (a) the positive control
demonstrates the metric detects degradation, and (b) on representative data the
default (within-rung) scope's accuracy loss is within tolerance for the measured
speedup. The tolerance is a research judgment set by the maintainer, recorded in
the results doc, not hardcoded.

## Testing

- **Disabled-path identity lock:** with `skip_frozen_reads=false`, corrected
  output is bit-identical to pre-opt4 master (a golden/lock test — opt4 must be a
  true no-op when off). This composes with the existing corrector lock tests.
- **Freeze-state unit test:** the streak counter increments/resets correctly on
  a synthetic `was_improved` sequence; within-rung scope resets at k-change;
  across-rung persists.
- **Positive control** (above) as an automated test asserting degradation.
- **Diagnostics completeness:** `frozen_reads_skipped` added to
  `CorrectorDiagnostics` + its constructor + the `corrector_errors` export (the
  three-site guard) so `fieldcount`-based completeness tests stay green.

## Risks

- **Silent accuracy loss on real data** — the whole reason for the opt-in default
  and the sign-off. Mitigated by within-rung default + the positive control +
  representative-data sign-off before any recommendation to enable it.
- **Early low-k freezing** (fact/risk 7) — the highest-divergence regime.
  Mitigated by within-rung default (reset at each k) and by N≥2 (a single
  no-improvement pass at low k does not freeze).
- **State-threading bug** — the streak counter is new cross-pass state; a
  threading error could freeze the wrong reads. Mitigated by the freeze-state
  unit test and the disabled-path identity lock (off ⇒ no state effect).
- **Interaction with `stop_on_no_change`** — if all reads freeze, the pass makes
  zero improvements and the existing zero-improvement early-stop fires, which is
  correct behavior (nothing left to do), but the test suite should confirm the
  frozen-skip doesn't cause premature k-advance that the exact path wouldn't.
- **Interaction with `sufficient_improvements`** (a *second*, distinct early-stop
  from `stop_on_no_change`) — `continue_current_k` is also gated on
  `sufficient_improvements(improvements_made, ...)`, an aggregate improvement-rate
  plateau check. Because frozen reads are structurally excluded from
  `improvements_made`, a rung with many frozen reads shows a lower improvement
  rate and can trip this plateau check **sooner** than the exact path. This is
  directionally safe (frozen reads can only *lower* the rate → early-stop fires
  sooner, never later, so the corrector never runs *more* than exact) but it is a
  distinct vector from the zero-improvement case and is named here for
  completeness.

## Out of scope

- Threading the flag up into `AssemblyConfig`/`assemble_genome` (a follow-on once
  the sign-off passes; opt4 uses the shallow direct-call path).
- opt3 (batched kernel, deferred). Any GPU/SIMD acceleration.
- Changing the default: opt4 ships off; enabling it by default would be a
  separate, evidence-gated decision.

## Acceptance criteria

- `skip_frozen_reads=false` (default) is a bit-identical no-op vs pre-opt4 master
  (locked by test).
- Per-read freeze streak + within-rung/across-rung scope implemented and
  unit-tested; `freeze_streak_threshold` and `freeze_across_rungs` honored.
- Frozen reads skipped via the existing `_skip_this_read_at` path; their frozen
  correction reaches the graph rebuild; `frozen_reads_skipped` diagnostic added
  (three-site guard satisfied).
- Accuracy sign-off run and written up: per-base + dnadiff/ANI deltas,
  skip-ON vs OFF, both scopes; the positive control demonstrates the metric
  detects degradation; speedup measured.
- Ship decision (enable-recommend or keep-off) recorded with the tolerance
  rationale in the results doc — no default change without that evidence.

## Rule-of-5 bead chain (proposed)

Per the campaign's Rule-of-5 pattern, opt4 decomposes into a dependency chain
(all children of `td-jbjd`, [B]):

1. **Draft** — implement the freeze streak state + skip hook + opt-in kwargs
   (shallow thread), default off; disabled-path identity lock test.
2. **Self-review** — "what did I miss?": diagnostics three-site guard, soft-EM
   interaction, `stop_on_no_change` interaction; freeze-state unit test.
3. **Architecture** — within-rung vs across-rung scope correctness; confirm the
   skip path composes with `skip_solid`/`hard_window` rather than conflicting.
4. **Edge cases** — all-reads-freeze, single-read, N=1 vs N≥2, k-boundary reset;
   the positive-control toy fixture.
5. **Verification** — the accuracy sign-off (per-base + dnadiff/ANI, both
   scopes), speedup measurement, and the ship-decision writeup.
