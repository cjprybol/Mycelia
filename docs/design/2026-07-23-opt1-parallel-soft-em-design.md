# opt1: Thread-safe soft-EM via deferred read-order reduce

Date: 2026-07-23 Status: Approved (design) Scope: `:scalable` iterative
corrector — parallelize the per-read Viterbi decode while soft-EM is enabled,
byte-identically to the current serial path.

## Context

The per-read Viterbi decode is ~80–90% of the corrector's wall time and is
embarrassingly parallel, but it runs serially whenever soft-EM is on — which is
always, for the production `:scalable` strategy. This is the second step of the
corrector performance campaign and its largest single lever (targeting a 3→16
core scaling). It follows the opt5 GC change (merged), which removed the forced
between-batch `GC.gc()` that would otherwise park decode workers.

The serial path is forced by one guard in `src/iterative-assembly.jl`:

```julia
use_parallel = enable_parallel && soft_weights === nothing
```

with a `@warn` that soft-EM accumulation "is sequential (race-free)". The
rationale is correct: the soft-EM accumulator is a single shared
`SoftEdgeWeightAccumulator` (a `Dict{Any, Float64}`) mutated per read, so a
naive `Threads.@threads` decode would race on it.

## Key facts (from a code survey of the post-opt5 baseline)

1. **Corrected-reads output is already order-deterministic under `@threads`.**
   Each thread writes `batch_results[i]` and `skip_flags[i]` by its own index
   (indexed writes into preallocated `Vector`s — `skip_flags` is a
   `Vector{Bool}` specifically to avoid the `BitVector` shared-word race), and
   results are collected into `updated_reads` in read order after the loop. No
   change needed.
2. **The soft-EM read (consumption) side is already parallel-safe.** It uses a
   read-only `_SoftEdgeWeightSnapshot` (an `IdDict`) rebound per task via
   `_with_soft_edge_weight_snapshot`, already wired inside the `@threads` loop.
   The accumulator is consumed only _between_ passes (this pass's accumulator
   becomes next pass's snapshot), never mid-pass.
3. **The per-read decode already stages into a private accumulator and merges at
   the end.** `improve_read_likelihood` allocates a per-read
   `staged_soft_weights`, accumulates into it during the E-step, and only on
   passing all decode contracts calls
   `_merge_soft_edge_weights!(soft_weights, staged)` into the shared
   accumulator. That merge is the sole shared mutation.
4. **A reduce primitive already exists.**
   `_merge_soft_edge_weights!(dest, staged)` is
   `dest.weights[e] += staged.weights[e]` over the Dict — exactly the fold
   needed to combine per-thread/per-read accumulators.

## The core problem: float non-associativity

The serial path folds per-read contributions staged₁, staged₂, …, staged_N into
the shared accumulator in ascending read order. Floating-point addition is not
associative, so any parallel scheme that changes the summation grouping produces
last-ULP-different weights. Because this pass's accumulator becomes the decode
snapshot for the next pass, a ULP difference in a weight can eventually change a
decode decision and thus a corrected read — breaking the campaign's
byte-identity requirement (the corrector lock tests assert byte-identical
corrected reads across passes). **The parallel design must replay the exact
serial fold order.**

## Approaches considered

|                | Approach                                                                                                                                                                                                          | Byte-identical weights?                                  | Memory                                                                  |
| -------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- | ----------------------------------------------------------------------- |
| **A (chosen)** | Deferred read-order reduce: parallelize the decode; each read writes its own staged accumulator to `batch_local[i]`; after the loop, serially fold `batch_local[1..N]` into the shared accumulator in read order. | Yes — exactly reproduces the serial left-fold.           | up to `batch_size` small accumulators per batch (freed after the fold). |
| B              | Per-thread accumulators, reduce in thread order.                                                                                                                                                                  | No — chunk-wise grouping ≠ serial left-fold (ULP drift). | T accumulators (T = threads).                                           |
| C              | Lock-striped / atomic shared Dict.                                                                                                                                                                                | No — adds in completion order (nondeterministic).        | 1 shared Dict.                                                          |

Only **A** meets the hard byte-identity requirement. The parallelized portion in
A is the expensive decode; the serial fold is only Dict additions, so A retains
essentially the full multi-core speedup. B and C are lower-memory but change the
fold grouping — they would require an approximate-accuracy sign-off, which
belongs to a separate (approximate) campaign step, not here.

## Design

### 1. Relax the guard

In `src/iterative-assembly.jl`, change the gate to
`use_parallel = enable_parallel` and remove the now-inaccurate "soft-EM is
sequential" `@warn`. Parallel + soft-EM now coexist.

### 2. Parallel branch: per-window ordered capture + flat read-order replay

The serial path folds soft-EM contributions onto a running accumulator in
**read×window order**: for an edge `e`, `0, +r1w1, +r1w2, +r2w1, …` (flat
sequential). Two contribution granularities exist:

- **Non-windowed decode** contributes **once per read** — it stages into a
  private `staged_soft_weights` and does a single terminal
  `_merge_soft_edge_weights!(soft_weights, staged)`. A per-read accumulator
  reproduces this exactly.
- **Windowed decode** (`windowed_decode=true`, which the production `:scalable`
  strategy sets) contributes **once per window** — it stages a fresh accumulator
  per window and merges each into the passed accumulator in window order.

A per-read accumulator therefore is NOT sufficient: it would pre-group a read's
windows as `(w1+w2)` before folding onto the running total `R`, giving
`R+(w1+w2)` where the serial path computes `(R+w1)+w2`. Float non-associativity
makes these differ whenever an edge is hit by ≥2 windows of one read that a
prior read also touched (a within-read cross-window repeat). To be bit-identical
on all paths, the parallel branch must replay the **per-window** contributions
in the exact serial flat order.

Mechanism — a capture "sink" instead of an internal merge:

- Add an optional
  `soft_weights_sink::Union{Nothing, Vector{SoftEdgeWeightAccumulator}}` kwarg
  threaded `improve_read_likelihood` → `try_viterbi_path_improvement`
  (non-windowed) and `improve_read_likelihood_windowed` (windowed). At each
  existing merge site, when `sink !== nothing`, `push!(sink, staged)` (append
  the staged accumulator in decode order) instead of
  `_merge_soft_edge_weights!(soft_weights, staged)`. Staging guards become
  `soft_weights !== nothing || sink !== nothing`. Serial behavior (sink
  `nothing`) is untouched.
- In the `if use_parallel` branch, allocate
  `batch_local = Vector{Vector{SoftEdgeWeightAccumulator}}(undef, n)` when
  `soft_weights !== nothing`. Each `@threads` iteration `i` passes
  `soft_weights_sink = batch_local[i]` (a fresh empty vector; the decode appends
  its per-read/per-window staged accumulators in order). Indexed write, no
  races.
- **After** the loop, replay flat in read×window order:
  `for i in 1:n; isassigned(batch_local, i) || continue; for staged in batch_local[i]; _merge_soft_edge_weights!(soft_weights, staged); end; end`.

For an edge `e`, this reproduces the serial running-total sequence
`0, +r1w1, +r1w2, +r2w1, …` exactly (non-windowed reads simply contribute a
length-1 list), so weights are bit-identical to the serial path on both windowed
and non-windowed decodes. The read-side snapshot rebinding and the
reads/skip-flags collection are unchanged.

### 3. Caller wiring — default parallel on when `nthreads > 1`

In `src/rhizomorph/assembly.jl`:

- Add `enable_parallel = Threads.nthreads() > 1` to the `:scalable` knobs in
  `_corrector_strategy_knobs`, and pass it in the `mycelia_iterative_assemble`
  call (which currently omits it, so it defaults to `false`).
- `:exhaustive` remains explicitly serial/exact (unchanged: `soft_em = false`,
  exact beam). The knob route keeps `:exhaustive` byte-identical while
  defaulting `:scalable` on.

The `gc_between_batches` keyword and `MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES` env
fallback from opt5 remain available; parallel decode is the scenario opt5's GC
change was ultimately for.

## Testing

Byte-identity under real parallelism is the load-bearing claim, and it is the
easiest to leave silently unverified: CI currently runs `Pkg.test` with **one**
thread, so a `Threads.@threads` test runs sequentially and proves nothing (and
an `nthreads() > 1` guard would auto-skip it). The test plan forces real
threads.

1. **Real-threads byte-identity (core lock).** A test that launches a child
   `julia -t 4` subprocess running the corrector twice on the same input —
   `enable_parallel = true` vs `false` — and asserts **identical corrected reads
   AND identical soft-EM weights** (full accumulator Dict, key set and every
   `Float64` value). One test case MUST enable `windowed_decode=true` on a
   repeat-heavy input (a short tandem repeat so an edge is hit by ≥2 windows of
   one read that a prior read also touched) — this is the case the per-window
   ordered capture exists for, and it would diverge under a naive per-read fold.
   Subprocess-based so it exercises genuine concurrency regardless of the parent
   harness's thread count. This is the primary guard against any hidden
   shared-state race beyond `soft_weights`.
2. **Equal-weights exposure.** The accumulator must be observable to the test —
   via a test hook or the corrector's returned metadata — so the weights can be
   compared, not just the reads.
3. **Parallel × GC on/off byte-identity.** The `enable_parallel = true` variant
   deferred from the opt5 review: assert byte-identical corrected reads with the
   between-batch GC on vs off, now that the parallel path is exercisable.
4. **Threaded CI job (in scope).** Add a CI job that runs the test suite with
   `JULIA_NUM_THREADS = 4` (e.g. `Pkg.test(julia_args = ["--threads=4"])` or a
   matrix entry with the env set), so the whole `test/4_assembly` suite runs
   multi-threaded, not only the dedicated subprocess test.
5. **Regression.** Re-verify the three corrector lock tests
   (`batched_viterbi_kernel`, `low_k_decode_gating`, `reassembly_graph_reuse`)
   and the opt5 GC test stay green.

## Risks

- **Peak memory:** up to `batch_size` per-read accumulators held per batch,
  freed after each fold. Bounded and tunable via `batch_size`; each accumulator
  holds only the edges one read touches.
- **Hidden shared mutable state** touched during decode beyond `soft_weights`.
  The survey found none (snapshot is read-only; reads/skips are indexed writes),
  and the subprocess byte-identity test is the guard that would catch any.
- **`@threads` load imbalance** across reads of differing length — acceptable;
  the decode dominates and scheduling does not affect correctness under Approach
  A (the fold is in read order regardless of which thread ran which read).

## Out of scope

Later campaign steps, each landing sequentially with byte-identity re-verified:
inner-DP allocation reduction, cross-pass memoization (approximate; separate
accuracy sign-off), and wiring the batched/SIMD Viterbi kernel into production.

## Acceptance criteria

- `:scalable` corrector runs the decode under `Threads.@threads` with soft-EM on
  and defaults parallel when `nthreads > 1`.
- Real-threads subprocess test: corrected reads and soft-EM weights
  byte-identical between `enable_parallel = true` and `false`.
- Parallel × GC on/off byte-identity test passes.
- Threaded CI job added and green.
- The three corrector lock tests and the opt5 GC test remain green.
- Measured multi-core speedup on a representative genome (recorded via the
  campaign benchmark), reported honestly against the serial baseline.
