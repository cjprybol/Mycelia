# opt2: Kill the inner-DP allocation storm (typed transitions, fused edge-scan)

Date: 2026-07-24
Status: Approved (design)
Scope: `_viterbi_correct_observation` hot decode loop — remove per-transition
allocations, byte-identically.

## Context

The per-read Viterbi decode is the corrector's dominant runtime term. Its inner
loop allocates heavily: an untyped `Vector{Any}` of per-edge `Dict{Symbol,Any}`
per active state per depth, plus boxed-`Int` diagnostics increments, plus a
redundant second full scan of every state's out-edges. This is the third step of
the corrector performance campaign (after opt5's GC gating and opt1's
parallelism). It matters more now that opt1 runs the decode across threads:
allocation pressure is multiplied across cores, so cutting per-transition
allocation reduces GC contention on every worker.

The load-bearing constraint is **byte-identical decode output** — the corrector
lock tests assert byte-identical corrected reads across passes, and the opt1
`parallel_soft_em_byte_identity` test asserts byte-identical soft-EM weights. Any
refactor must preserve exact iteration order, exact floating-point accumulation
order, and the weight floor.

## Key facts (post-opt1 baseline)

Files: `src/viterbi-next.jl` (hot loop) and
`src/rhizomorph/algorithms/path-finding.jl` (transition builder + weight scan).

1. **The transition container.** The hot loop (`_viterbi_correct_observation`,
   depth loop ~1411-1552) calls `_get_valid_transitions(graph, vertex, strand)`
   (path-finding.jl ~531), which builds `transitions = []` (a `Vector{Any}`) and
   `push!`es a `Dict{Symbol,Any}` per edge (`_maybe_push_transition!`,
   path-finding.jl ~515-526) with four keys: `:target_vertex` (dst label),
   `:target_strand` (`StrandOrientation`), `:probability` (`Float64`),
   `:edge_data` (edge-data struct).
2. **Enumeration order + filters.** Directed/production path iterates
   `Graphs.outneighbors(graph.graph, src_code)` in ascending dst-code order, with
   a strand filter `_normalize_strand(edge_data.src_strand) == strand` and a
   positive-weight filter (`_edge_transition_weight(edge_data) <= 0.0` skips).
3. **The double scan.** The loop separately calls
   `_total_outgoing_weight(graph, vertex, strand)` (path-finding.jl ~760) which
   scans the **same** `outneighbors` set in the **same** order with the **same**
   strand filter, summing `total += _edge_transition_weight(edge_data)` and
   returning `max(total, _KSP_MIN_WEIGHT)` (`_KSP_MIN_WEIGHT = 1e-10`).
4. **Diagnostics.** `diagnostics::Dict{Symbol,Any}` counters are bumped inside
   the loop (`:expanded_states`, `:skipped_transitions`, `:generated_states`,
   `:successor_bounded`, `:beam_pruned`, `:margin_pruned`, and the per-depth
   `:retained_states` / `:cumulative_retained_states` / `:max_retained_states` /
   `:completed_steps`). Each `+= 1` on an `Any`-typed Dict boxes an `Int`. All
   read-backs are in-loop self-reads; external consumers read only final values.
5. **Consumers of `_get_valid_transitions`.** Besides the hot loop, the return is
   consumed by `_calculate_transition_probabilities`, `_sample_transition`,
   `_top_b_transitions` (sorts by `(_edge_transition_weight(t[:edge_data]),
   string(t[:target_vertex]))`), and a few other sites — all read `:edge_data`
   and/or `:probability`/`:target_vertex`.
6. **Byte-identity surfaces.** Float accumulations: the `total_out` running sum;
   `transition_prob = edge_w / total_out`; `next_score = state_score +
   log(transition_prob) + emission_score`; the cumulative emission
   `get(active_emissions, state, 0.0) + emission_score`. Ordered iterations that
   drive tie-breaks: the `outneighbors` order, the `active_scores` Dict iteration
   (line ~1420), the `transitions` iteration, and the beam-prune sorts (which
   tie-break on `hash(state_tuple)` — so the **state-tuple key type must not
   change**).

## Design (changes 1-3; change 4 deferred)

### 1. Typed transition struct

Define a concrete, parametric struct and replace the per-edge Dict:

```julia
struct _Transition{L, E}
    target_vertex::L
    target_strand::StrandOrientation
    probability::Float64
    edge_data::E
end
```

`_get_valid_transitions` returns a typed `Vector{_Transition{L,E}}` (element type
inferred from the first push, or built via a typed accumulator). Update every
consumer from `t[:key]` to `t.key` field access
(`_maybe_push_transition!`, the hot-loop consumer at ~1445-1478,
`_calculate_transition_probabilities`, `_sample_transition`,
`_top_b_transitions`, and the other `t[:edge_data]` sites). Byte-identical: the
push order (outneighbors) is unchanged, the four values are unchanged, and the
state-tuple key type (`(next_vertex, next_strand)`) is untouched — so
`hash`-based tie-breaks are unaffected.

### 2. Fuse `_total_outgoing_weight` into the transition pass

Replace the two separate `outneighbors` scans (transition build at ~1422 +
weight sum at ~1428) with one pass that returns both:

```julia
# returns (transitions::Vector{_Transition}, total_out::Float64)
function _valid_transitions_and_total(graph, vertex, strand) ... end
```

The single pass iterates `outneighbors` in ascending dst-code order; for each
strand-matched edge it (a) adds `_edge_transition_weight(edge_data)` to `total`
(including 0.0-weight edges — they contribute exactly `0.0`, so the sum is
bit-identical to the current `_total_outgoing_weight`), and (b) pushes a
`_Transition` iff the weight is `> 0.0`. Return `max(total, _KSP_MIN_WEIGHT)` for
`total_out`. **Order invariants:** `total_out` is summed over the full
strand-matched set in `outneighbors` order and fully computed before any
`edge_w / total_out` divide and before the top-B truncation — exactly as today.
Keep `_total_outgoing_weight` and `_get_valid_transitions` as thin wrappers (or
in place) for their other call sites, so only the hot loop switches to the fused
path.

### 3. Hoist diagnostics to local `Int` accumulators

Replace each in-loop `diagnostics[:key] += 1` (and the `get(diagnostics, :key,
0) + 1` forms) with unboxed `Int` locals, initialized from the pre-loop values
where a running max / cumulative sum is seeded, updated in the loop, and written
back into `diagnostics` once after the loop. Final values are identical (the
counters are only self-read in-loop); the win is removing per-transition boxed
allocation.

### 4. (Deferred) Reuse per-depth Dicts

**Out of scope — byte-identity risk.** `next_predecessors` is retained for
traceback (can't be reused); reusing `next_scores`/`next_emissions` via `empty!`
retains hash-table capacity and changes iteration order versus a fresh Dict,
which can flip an equal-score tie-break in the next depth's `active_scores`
iteration — a byte-identity break for a marginal allocation gain. Recorded as a
possible approximate-tier item, not part of opt2.

## Testing

- **Golden-master decode test:** capture `_viterbi_correct_observation`'s full
  output (decode path, best score, and the diagnostics Dict) on a fixed
  graph+read+config BEFORE the refactor; assert bit-identical AFTER. This
  directly locks the hot loop independent of the full corrector.
- **Existing lock tests** (`low_k_decode_gating`, `reassembly_graph_reuse`,
  `batched_viterbi_kernel`) assert byte-identical corrected reads through the
  full corrector — must stay green.
- **opt1 `parallel_soft_em_byte_identity`** asserts byte-identical soft-EM
  weights (parallel and serial) — must stay green (opt2 is inside the decode
  both paths share).
- **Allocation check (optional):** a `@allocated` micro-benchmark on the hot loop
  before/after to confirm the allocation reduction (directional; not a pass/fail
  gate).

## Risks

- **Type instability from `_Transition{L,E}`** if `L`/`E` aren't concrete at the
  call site — verify the element type is inferred concretely (the win depends on
  it). Mitigate by constructing the vector with the known element type.
- **Missed consumer** of the old Dict interface — grep every `t[:` / `[:target_`
  / `[:edge_data` / `[:probability` / `[:target_strand` site and convert. A
  missed site is a `MethodError` (loud), not a silent byte-identity break.
- **Fusion reordering** the float sum — mitigated by summing in `outneighbors`
  order over the full strand-matched set (identical to today).

## Out of scope

Change 4 (Dict reuse). Later campaign steps: opt4 (cross-pass memo, approximate)
and opt3 (batched/SIMD kernel wiring).

## Acceptance criteria

- Per-edge `Dict{Symbol,Any}` and `Vector{Any}` transitions replaced by the typed
  struct/vector; all consumers converted to field access.
- The two out-edge scans fused into one; `total_out` bit-identical.
- Diagnostics counters hoisted to `Int` locals; final diagnostics identical.
- Golden-master decode test passes (path + score + diagnostics bit-identical);
  all corrector lock tests + the opt1 soft-EM byte-identity test stay green.
- Allocation reduction demonstrated on the hot loop (reported as measured).
