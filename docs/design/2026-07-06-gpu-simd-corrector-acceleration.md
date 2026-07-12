# ADR: GPU / SIMD acceleration of the corrector's residual per-read decode

- **Status:** Proposed
- **Date:** 2026-07-06
- **Scope:** `src/viterbi-next.jl` decode core (`_viterbi_correct_observation`,
  `_correct_metagraphs_next_observations`)
- **Decision type:** Architecture direction + phased plan. This ADR proposes a
  reformulation and a two-phase acceleration path; it does not implement kernels.
  No `src/` change is made by this document.

## Context

The graph-as-HMM corrector decodes each read against a de Bruijn / Rhizomorph
assembly graph with a pair-HMM / Viterbi recurrence. For a single read the
recurrence is `O(read_length x graph_states_reachable)` and is **strictly
sequential in read position**: the score frontier at read position `d+1` depends
on the frontier at position `d` (see `src/viterbi-next.jl:1050-1141`). There is
no intra-read parallelism to exploit without changing the algorithm.

The available parallelism is **across reads**. Every read is decoded against the
*same* graph, and decoding only *reads* graph structure (topology, edge weights,
emission parameters) — it never mutates it. A batch of `N` reads decoded against
one shared, immutable graph therefore has **zero write contention**: the graph is
a read-only oracle queried by `Rhizomorph._get_valid_transitions`,
`Rhizomorph._total_outgoing_weight`, and `_call_viterbi_state_emission_logp`
(`src/viterbi-next.jl:1060,1066,1081`). This is the structural fact that makes
batched acceleration tractable.

### Why acceleration is needed

The current decoder is a CPU scalar Viterbi with a `Dict`-based frontier. Cell-
updates-per-second (CUPS) throughput for scalar pair-HMM/Viterbi on CPU is on the
order of **0.02-0.6 GCUPS**. On real branchy graphs the reachable
`(vertex, strand)` frontier grows with read-length depth until it is capped by
the beam guard (`src/viterbi-next.jl:1108-1119` documents a 21-billion-allocation
hard crash on a 48 kb phage before beam pruning was added). Even with beam
pruning, the residual **hard windows** — regions where the graph is genuinely
branchy and the frontier stays wide — dominate wall-clock and are the reason
`strategy=:exhaustive` decoding times out at scale.

Two external reference points bound the opportunity:

- **gpuPairHMM** (length-binned read batches + a warp-shuffle wavefront that
  propagates the antidiagonal of the DP matrix through register shuffles rather
  than shared memory) reaches **1.82 TCUPS on an L40S** — roughly **44x** a
  128-thread CPU baseline and **174x** the prior best GPU pair-HMM. That is a
  ~3-to-4-order-of-magnitude headroom over the scalar-CPU ceiling above.
- **NVIDIA Parabricks** established the production precedent that this class of
  genomics DP is GPU-appropriate: HaplotypeCaller's pair-HMM step drops from
  ~36 h to ~33 min when moved to GPU.

A **CPU middle ground** also exists and is lower-lift: **abPOA-style windowed
SIMD partial-order alignment** vectorizes the DP across the alphabet / adjacent
cells and reports **2.7-9.5x** over scalar POA. Crucially, the hard-window work
already in flight (`<=500 bp` bounded windows) sets up exactly the fixed,
length-bounded work units that SIMD-POA needs — the window bound *is* the SIMD
tile bound.

### The current structure blocks batching

`_viterbi_correct_observation` decodes **one** read and is called in a per-read
Julia loop (`src/viterbi-next.jl:350-374`). Three concrete properties block a
batched formulation:

1. **`Dict`-based frontier.** The active frontier is
   `Dict{Tuple{label_type, StrandOrientation}, Float64}`
   (`src/viterbi-next.jl:1000`, `:1052`). A hash map is pointer-chasing,
   non-contiguous, and per-read-allocated — it cannot be uploaded to a GPU as a
   dense buffer, cannot be SIMD-scanned, and defeats coalesced memory access.
   The `next_scores` / `next_predecessors` dicts are rebuilt every depth step
   (`src/viterbi-next.jl:1052-1056`) with per-key `haskey`/insert
   (`src/viterbi-next.jl:1097-1100`).

2. **Per-read wrapping.** State — `active_scores`, `predecessors_by_depth`,
   `diagnostics`, `target_vertex` — is scoped inside the single-read function
   (`src/viterbi-next.jl:1000,1043,977,973`) and returned as one
   `ViterbiDecodingResult`. There is no batch object: `N` reads means `N`
   independent function invocations each allocating its own dicts, so there is no
   shared buffer to vectorize over.

3. **Loop nesting `for-read { for-depth }`.** The read loop is the outer loop
   (`src/viterbi-next.jl:350`) and the depth loop is inner
   (`src/viterbi-next.jl:1050`). Within a read the depth loop is inherently
   sequential, so the outer read loop is the only parallel axis — but it is
   expressed as sequential Julia iteration over `observations`, with each
   iteration touching its own dicts. The parallel axis is on the *outside*,
   invisible to any vectorizer or GPU launch.

## Decision

Reformulate the decode core so the **read axis is the innermost, contiguous,
data-parallel axis**, then accelerate that axis in two phases. Specifically:

### Required reformulation

1. **Array-based frontier.** Replace the per-read
   `Dict{(vertex,strand), Float64}` with a dense score array indexed by a
   contiguous state id. Precompute a stable `state_id <- (vertex, strand)` map
   once for the graph (it is immutable during a decode). The frontier becomes a
   `Float32`/`Float64` vector (or `N_reads x N_states` matrix for a batch) with
   `-Inf` for unreached states. Predecessors become an integer array of
   `state_id`s, not a dict. This is the single change that unlocks everything
   downstream: dense arrays are SIMD-scannable, GPU-uploadable, and
   coalesced-addressable.

2. **Depth-outer, read-inner loop nesting: `for-depth { for-read }`.** Invert the
   nesting relative to today. The depth (read-position) step is the outer,
   sequential axis — every read in the batch advances one position in lock-step.
   The read axis becomes the inner, fully parallel axis (SIMD lanes on CPU, warp
   lanes / thread blocks on GPU). At each depth, all `N` reads perform the same
   frontier-relaxation against the shared graph. This is the loop-interchange that
   exposes the zero-contention read parallelism to the hardware.

3. **Length-binned batches.** Group reads by length into bins so a batch shares a
   common depth bound and a common (or padded) frontier width. This is the
   gpuPairHMM batching discipline: equal-length work units keep every lane busy
   for the same number of depth steps and avoid divergence / wasted padding.
   The hard-window bound (`<=500 bp`) already produces naturally length-bounded
   units for this binning.

4. **Shared immutable graph.** Pass the graph (topology, edge weights, emission
   tables) as a single read-only structure-of-arrays uploaded once per batch, not
   per read. All lanes/threads read it; none writes it. Diagnostics that are
   currently mutated in the hot loop (`diagnostics[:expanded_states] += 1`, etc.,
   `src/viterbi-next.jl:1061,1068,1096`) must move out of the inner kernel to
   per-lane counters reduced after the batch, so the kernel body has no shared
   mutable scalar.

### Phased plan

**Phase A — CPU windowed SIMD-POA on residual hard windows (near-term, lower
lift).** Apply abPOA-style SIMD partial-order alignment to the bounded
(`<=500 bp`) hard windows produced by the hard-window work. This requires no new
dependency (SIMD is native via `@simd` / `SIMD.jl` / LLVM auto-vectorization on
the array frontier), targets the exact regions that dominate `:scalable`'s
residual cost, and delivers the reference 2.7-9.5x. Phase A is gated only on the
array-frontier reformulation (items 1-2 above) and the hard-window bound, both of
which are already in motion. It is the risk-reducing first step: it validates the
array frontier and the depth-outer loop on CPU before any GPU commitment.

**Phase B — GPU warp-shuffle wavefront kernel (GPU-gated).** Implement the
gpuPairHMM-style kernel: length-binned batches, one warp (or block) per read,
antidiagonal wavefront propagated by warp-shuffle so the DP frontier lives in
registers rather than shared memory. Target **KernelAbstractions.jl** for a
backend-agnostic kernel (CUDA / ROCm / Metal / oneAPI) with **CUDA.jl** as the
first concrete backend. This is **GPU-gated**: it is only entered when a GPU is
present and the batch is large enough to amortize upload; on CPU-only hosts the
corrector transparently falls back to Phase A / the scalar path. `Project.toml`
currently has **no GPU or accelerator dependency** (confirmed greenfield — no
CUDA / KernelAbstractions / Metal / oneAPI / AMDGPU), so the GPU stack is added
behind a package extension (weak dependency) to keep CPU-only installs lean.

### How it composes with the tiered strategy

The corrector's tiered strategy selects decode effort per read:

- **`strategy=:scalable`** decodes the easy majority cheaply and isolates a
  residual of **hard windows** for expensive exact decode. Phases A and B
  accelerate *exactly that residual*: the hard-window decode is where the wide
  frontier lives, so SIMD-POA (A) and the GPU wavefront (B) attack the part of
  `:scalable` that currently dominates wall-clock. `:scalable` gets faster
  without changing its accuracy contract.
- **`strategy=:exhaustive`** decodes every read exactly and currently **times out
  at scale** for the same reason single reads blow up the frontier. With a
  batched, length-binned GPU wavefront the per-read cost drops by the
  gpuPairHMM-class factor (tens-to-hundreds x), which is what makes exhaustive
  decoding *viable at scale* rather than a small-input-only mode. The beam guard
  (`config.beam_width`, `src/viterbi-next.jl:1108-1119`) remains the exactness
  knob: with `beam_width = typemax(Int)` the accelerated path must remain
  byte-identical to the scalar decoder on the B8 correction fixtures — the
  acceleration changes *where* the arithmetic runs, not the recurrence.

Both tiers keep the shared read-only graph; acceleration is orthogonal to
strategy selection and slots under it.

## Prototype scope (smallest batched-decode proof-of-concept)

The minimal PoC that de-risks the reformulation, before any kernel work:

1. **Fixed emission model, no strand branching.** Single-strand
   (`strand_mode = :single`), quality-aware emission held constant, so the PoC
   isolates the frontier data structure, not emission edge cases.
2. **Array frontier + depth-outer loop, on CPU, single-threaded.** Reimplement
   `_viterbi_correct_observation` for a *batch* of reads as
   `for depth { for read { relax against dense frontier } }` using dense score /
   predecessor arrays and a precomputed `state_id` map. No SIMD, no GPU yet —
   just prove the reformulated recurrence reproduces the scalar decoder.
3. **Equivalence oracle.** On a small fixed graph and a handful of reads
   (including one hard-window read), assert the batched array decoder returns
   **byte-identical paths and scores** to the current `Dict`-based decoder with
   `beam_width = typemax(Int)`. This equivalence test is the acceptance gate and
   the regression guard for Phases A and B.
4. **Then, and only then**, add `@simd`/`SIMD.jl` over the read axis (Phase A) and
   measure speedup on the hard-window batch; defer the GPU kernel (Phase B) to a
   separate PR behind a KernelAbstractions extension.

Success criterion for the PoC: identical output to the scalar decoder on the
equivalence oracle, with the read axis expressed as the inner loop over a dense
frontier — i.e. the data layout a vectorizer/GPU can consume, proven correct at
toy scale before optimizing.

## Consequences

- **Positive:** unlocks the across-read parallelism that the current per-read
  `Dict` structure hides; makes `:exhaustive` viable at scale; accelerates
  `:scalable`'s residual; adds no GPU dependency to CPU-only installs (extension-
  gated); Phase A ships value with no new dependency.
- **Negative / risk:** the array-frontier reformulation touches the hot decode
  core, so it must be landed behind the byte-identical equivalence oracle before
  any optimization. A dense `N_states` frontier can waste memory when the true
  frontier is sparse relative to total graph states; length-binning and beam
  pruning bound this, but very large graphs may need a compressed/active-set
  representation as a follow-up. GPU work is gated on hardware availability and on
  Phase A having validated the reformulation.
- **Deferred:** multi-strand batching, sparse/compressed frontier for very large
  graphs, and multi-GPU sharding are out of scope for this ADR.

## Evidence (current structure)

| Claim | File:line |
| --- | --- |
| Per-read decode entry point | `src/viterbi-next.jl:953` |
| Per-read loop over observations (`for-read` outer) | `src/viterbi-next.jl:350-374` |
| `Dict`-based active frontier | `src/viterbi-next.jl:1000`, `:1052` |
| Depth loop inner + sequential dependence | `src/viterbi-next.jl:1050-1141` |
| Graph queried read-only (transitions/weights/emission) | `src/viterbi-next.jl:1060,1066,1081` |
| Per-state `haskey`/insert into next-frontier dict | `src/viterbi-next.jl:1097-1100` |
| In-hot-loop mutable diagnostics counters | `src/viterbi-next.jl:1061,1068,1096` |
| Beam guard = exactness knob (`typemax(Int)` -> exact) | `src/viterbi-next.jl:1108-1119` |
| No GPU/accelerator dep (greenfield) | `Project.toml` (`[deps]`: no CUDA/KernelAbstractions/Metal/oneAPI/AMDGPU) |
