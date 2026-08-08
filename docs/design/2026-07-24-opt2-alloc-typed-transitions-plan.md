# opt2 Typed Transitions + Fused Edge-Scan Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Kill the inner-DP allocation storm in the hot Viterbi decode by
replacing the per-edge `Vector{Any}`+`Dict{Symbol,Any}` transitions with a typed
struct built in a single fused out-edge pass, and hoisting boxed diagnostics
counters to `Int` locals — byte-identically.

**Architecture:** Add a concrete `struct _Transition{L,E}` and one fused
function
`_valid_transitions_and_total(graph, vertex, strand) -> (Vector{_Transition{L,E}}, Float64)`
that does a single `outneighbors` pass (sum + build). Wire it into the two
decode call sites (substitution hot loop + memoized indel successor batch),
converting their consumers from `t[:key]` to `t.key`. Leave the Dict-form
`_get_valid_transitions`/`_total_outgoing_weight` untouched for the ~10 cold
consumers. Hoist in-loop diagnostics-Dict mutations to `Int` locals merged once
after the loop.

**Tech Stack:** Julia; MetaGraphsNext; Graphs; `Test` stdlib.

## Global Constraints

- **BYTE-IDENTITY is the hard requirement.** Corrector lock tests assert
  byte-identical corrected reads; the opt1 `parallel_soft_em_byte_identity` test
  asserts byte-identical soft-EM weights. Re-verify after every task. Preserve:
  `Graphs.outneighbors` ascending dst-code order; the
  `total += _edge_transition_weight(edge_data)` summation order (fuse without
  reordering; `total_out` fully summed before any `edge_w / total_out` divide
  and before top-B truncation); the `max(total, _KSP_MIN_WEIGHT)` floor
  (`_KSP_MIN_WEIGHT = 1e-10`); the strand filter
  `_normalize_strand(edge_data.src_strand) == strand`; and the state-tuple key
  type (unchanged → `hash`-based beam-prune tie-breaks unaffected).
- **`.jl` edits MUST bypass the formatter hook.** The PostToolUse Edit/Write
  hook whole-file-churns `.jl` to non-SciML. Apply every `src/`/`test/` `.jl`
  change via a Bash-invoked Python/heredoc in-place script (the hook fires on
  Edit/Write tools, not Bash). After each edit, `git -C <wt> diff --stat <file>`
  to confirm a minimal diff; if churned, `git checkout -- <file>` and redo via
  Python.
- **SciML style** (`.JuliaFormatter.toml`: `style="sciml"`), trailing commas,
  4-space indent.
- **Julia:** `LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. ...`.
  Parse-check each edited `.jl` with `Meta.parseall`.
- **No `td-*`** in committed non-comment artifacts (test names, docs prose);
  allowed in commit messages + code comments.
- **Worktree:** `/Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc`,
  branch `cp/opt2-alloc` (off `origin/master` @ 0804d0ec). Commit per task; **do
  NOT merge** (human gate).
- **`git -C <wt> …`** for all git (never `cd <wt> && git` — a guard blocks it).

## Reference: verified baseline locations (post-opt1)

- Hot substitution loop: `src/viterbi-next.jl` `_viterbi_correct_observation`
  (depth loop ~1411-1552). Transition build @~1422 `_get_valid_transitions`,
  total @~1428 `_total_outgoing_weight`, top-B @~1442 `_top_b_transitions`,
  consumer @~1445-1478 (`transition[:target_vertex/:target_strand/:edge_data]`).
- `_top_b_transitions`: `src/viterbi-next.jl:1213` — sort key
  `(_edge_transition_weight(t[:edge_data]), string(t[:target_vertex]))`,
  `rev=true`. ONLY caller is @1442 (verified) → safe to convert to field access.
- Memoized indel site: `src/viterbi-next.jl` `_indel_decode_successors!` (def
  ~1809), `_get_valid_transitions` @~1820, `_total_outgoing_weight` @~1830,
  consumer @~1837-1841.
- Dict-form builders (LEAVE UNCHANGED):
  `src/rhizomorph/algorithms/path-finding.jl` `_get_valid_transitions` (~531),
  `_maybe_push_transition!` (~515), `_total_outgoing_weight` (~760),
  `_edge_transition_weight` (~725), `_KSP_MIN_WEIGHT = 1e-10` (~721). Cold
  consumers of the Dict form (LEAVE UNCHANGED): `information-theory.jl:30/57`,
  `generation.jl:147/233`, `batched-viterbi-poc.jl:126/128`,
  `path-finding.jl:407/467/665/1110/1400/1838`.
- Diagnostics (hoist): `src/viterbi-next.jl` init ~1326-1347 + pre-seed
  ~1379-1381; in-loop `:expanded_states` ~1423, `:skipped_transitions`
  ~1430/1450/1464, `:successor_bounded` ~1440, `:generated_states` ~1469,
  `:beam_pruned` ~1499, `:margin_pruned` ~1518, per-depth
  `:retained_states`/`:cumulative_retained_states`/`:max_retained_states`/`:completed_steps`
  ~1534-1537.

---

### Task 1: `_Transition` struct + fused `_valid_transitions_and_total` + equivalence lock

**Files:**

- Modify: `src/rhizomorph/algorithms/path-finding.jl` (add struct + fused fn
  near `_get_valid_transitions`/`_total_outgoing_weight`)
- Create: `test/4_assembly/valid_transitions_fused_test.jl`

**Interfaces:**

- Produces: `Mycelia.Rhizomorph._Transition{L,E}` (fields `target_vertex::L`,
  `target_strand::StrandOrientation`, `probability::Float64`, `edge_data::E`);
  `Mycelia.Rhizomorph._valid_transitions_and_total(graph, vertex, strand)::Tuple{Vector, Float64}`.
- Consumes: existing `_get_valid_transitions`, `_total_outgoing_weight`,
  `_edge_transition_weight`, `_normalize_strand`, `_KSP_MIN_WEIGHT`,
  `_maybe_push_transition!` semantics.

- [ ] **Step 1: Write the failing equivalence test** (Python-bypass writer). The
      test builds a small qualmer graph, and for each vertex/strand asserts the
      fused function returns transitions whose
      `(target_vertex, target_strand, probability, edge_data)` equal — in the
      SAME order — the Dict-form `_get_valid_transitions`, and `total_out`
      bit-identical (`===`) to `_total_outgoing_weight`.

```julia
# FUSED VALID-TRANSITIONS EQUIVALENCE (opt2)
# _valid_transitions_and_total must return, per (vertex,strand): the SAME
# transitions in the SAME order as the Dict-form _get_valid_transitions, and a
# total_out bit-identical to _total_outgoing_weight. This is the byte-identity
# lock for the fused hot-loop function.
import Test
import Mycelia
import FASTX
import Random

Test.@testset "fused valid-transitions equivalence (opt2)" begin
    R = Mycelia.Rhizomorph
    rng = Random.MersenneTwister(31337)
    ref = join(rand(rng, ['A', 'C', 'G', 'T'], 400))
    reads = [FASTX.FASTQ.Record("r$i", ref[s:(s + 79)], String(fill('I', 80)))
             for (i, s) in enumerate(rand(rng, 1:321, 60))]
    k = 11
    graph = R.build_qualmer_graph(reads, k; mode = :canonical)
    checked = 0
    for label in Mycelia.MetaGraphsNext.labels(graph)
        for strand in (R.Forward, R.Reverse)
            dict_tx = R._get_valid_transitions(graph, label, strand)
            dict_total = R._total_outgoing_weight(graph, label, strand)
            fused_tx, fused_total = R._valid_transitions_and_total(graph, label, strand)
            Test.@test fused_total === dict_total                       # bit-identical
            Test.@test length(fused_tx) == length(dict_tx)
            for (a, b) in zip(fused_tx, dict_tx)
                Test.@test a.target_vertex == b[:target_vertex]
                Test.@test a.target_strand == b[:target_strand]
                Test.@test a.probability === b[:probability]            # bit-identical
                Test.@test a.edge_data === b[:edge_data]
            end
            checked += 1
        end
    end
    Test.@test checked > 0                                             # non-vacuous
end
```

- [ ] **Step 2: Run to verify it fails.**

Run:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e 'include("test/4_assembly/valid_transitions_fused_test.jl")'`
Expected: FAIL — `UndefVarError: _valid_transitions_and_total` (function/struct
not defined).

- [ ] **Step 3: Add the struct + fused function** (Python-bypass) in
      `src/rhizomorph/algorithms/path-finding.jl`, immediately after
      `_total_outgoing_weight` (or after `_maybe_push_transition!`). The fused
      pass MUST mirror `_get_valid_transitions` (directed + undirected branches)
      and `_total_outgoing_weight` exactly. Read both current functions verbatim
      first and mirror their `outneighbors`/`edge_labels` iteration + strand
      filter.

```julia
struct _Transition{L, E}
    target_vertex::L
    target_strand::StrandOrientation
    probability::Float64
    edge_data::E
end

# opt2: one outneighbors pass building the typed transition list AND summing the
# outgoing weight — fuses _get_valid_transitions + _total_outgoing_weight. Sum is
# over ALL strand-matched edges in outneighbors order (0.0-weight edges add 0.0,
# so bit-identical to _total_outgoing_weight); a _Transition is pushed only for
# weight > 0.0 (identical to the positive-weight skip in _maybe_push_transition!).
function _valid_transitions_and_total(graph, vertex_label, strand)
    total = 0.0
    transitions = _Transition[]
    haskey(graph, vertex_label) || return transitions, max(total, _KSP_MIN_WEIGHT)
    if Graphs.is_directed(graph.graph)
        src_code = MetaGraphsNext.code_for(graph, vertex_label)
        for dst_code in Graphs.outneighbors(graph.graph, src_code)
            target_vertex = MetaGraphsNext.label_for(graph, dst_code)
            edge_data = graph[vertex_label, target_vertex]
            _normalize_strand(edge_data.src_strand) == strand || continue
            w = _edge_transition_weight(edge_data)
            total += w
            w <= 0.0 && continue
            push!(transitions,
                _Transition(target_vertex, _normalize_strand(edge_data.dst_strand),
                    w, edge_data))
        end
    else
        for edge_labels in MetaGraphsNext.edge_labels(graph)
            if length(edge_labels) == 2 && edge_labels[1] == vertex_label
                target_vertex = edge_labels[2]
                edge_data = graph[vertex_label, target_vertex]
                _normalize_strand(edge_data.src_strand) == strand || continue
                w = _edge_transition_weight(edge_data)
                total += w
                w <= 0.0 && continue
                push!(transitions,
                    _Transition(target_vertex,
                        _normalize_strand(edge_data.dst_strand), w, edge_data))
            end
        end
    end
    return transitions, max(total, _KSP_MIN_WEIGHT)
end
```

> NOTE (shipped form differs from this sketch): the shipped `_Transition` list
> is allocated as `_Transition{L, E}[]` (concrete eltype, pulled from the
> MetaGraph type parameters `L`/`E` via dispatch — fix round 1), not the
> abstract `_Transition[]` shown above; and the shipped struct is the 3-field
> `_Transition{L, E}(target_vertex, target_strand, edge_data)` after the
> vestigial `probability::Float64` field was dropped (no src reads it — the hot
> consumers recompute the weight via `_edge_transition_weight(edge_data)`), so
> the two `push!` sites take no `w` argument.

NOTE on order: `_maybe_push_transition!` computes `probability` and applies the
`<= 0.0` skip BEFORE reading `dst_strand`; here
`w = _edge_transition_weight(...)` is the same `probability`. The sum adds `w`
for every strand-matched edge (matching `_total_outgoing_weight`, which adds
`_edge_transition_weight` for every strand-matched edge, including those with
`w == 0.0`). Verify by the equivalence test — do not reorder the `total += w`
relative to the `w <= 0.0 && continue`.

- [ ] **Step 4: Run to verify it passes.**

Run:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e 'include("test/4_assembly/valid_transitions_fused_test.jl")'`
Expected: PASS — fused transitions + total bit-identical to the Dict form across
all vertices/strands.

- [ ] **Step 5: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc add src/rhizomorph/algorithms/path-finding.jl test/4_assembly/valid_transitions_fused_test.jl
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc commit -m "feat: typed _Transition + fused _valid_transitions_and_total (opt2)

Single outneighbors pass building a typed transition vector AND summing total_out,
equivalent bit-for-bit to _get_valid_transitions + _total_outgoing_weight
(equivalence-locked by test). Not yet wired into the hot loop.

td-cppm / td-jbjd opt2"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc push
```

---

### Task 2: Wire the fused typed path into the substitution hot loop

**Files:**

- Modify: `src/viterbi-next.jl` (`_viterbi_correct_observation`
  ~1422/1428/1442/1445-1478; `_top_b_transitions` ~1213)

**Interfaces:**

- Consumes: `_valid_transitions_and_total` (Task 1), `_Transition` field access.

- [ ] **Step 1: Capture the golden baseline** — before editing, record the
      current decode output as a committed fixture so the refactor is locked
      even beyond the full lock tests. Write
      `test/4_assembly/viterbi_decode_golden_test.jl` (Python-bypass) that
      builds a FIXED graph+read+config, runs
      `Mycelia.Rhizomorph._viterbi_correct_observation(...)`, and asserts its
      `.path`, `.best_score` (`===`), and full `.diagnostics` Dict equal
      hardcoded expected values. To get the expected values: run the builder
      once on CURRENT code, print the outputs, and paste them in as the
      literals. (This is a characterization/golden test — it must be written
      against pre-refactor output so it stays RED-able.)

Run (to capture):
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e '<build fixture; @show result.path result.best_score result.diagnostics>'`
Then encode those exact values as `Test.@test` literals in the golden test and
confirm it PASSES on current code.

- [ ] **Step 2: Convert `_top_b_transitions` to field access** (Python-bypass).
      At `src/viterbi-next.jl:1213-1220`, change the sort key from
      `t[:edge_data]`/`t[:target_vertex]` to `t.edge_data`/`t.target_vertex`.
      (Only caller is the hot loop @1442, verified — no cold caller.)

```julia
function _top_b_transitions(transitions, b::Int)
    length(transitions) <= b && return transitions
    return partialsort(transitions,
        1:b;
        by = t -> (Rhizomorph._edge_transition_weight(t.edge_data),
            string(t.target_vertex)),
        rev = true)
end
```

(Read the current body first; preserve its exact structure — only the two
`t[:x]` → `t.x` edits.)

- [ ] **Step 3: Wire the fused call + convert the consumer** (Python-bypass).
      Read the current ~1420-1478 region verbatim. Replace the paired calls:
  - `transitions = Rhizomorph._get_valid_transitions(graph, current_vertex, current_strand)`
    (@~1422) AND
    `total_out = Rhizomorph._total_outgoing_weight(graph, current_vertex, current_strand)`
    (@~1428) → a single
    `transitions, total_out = Rhizomorph._valid_transitions_and_total(graph, current_vertex, current_strand)`.
    Preserve the relative position so `total_out` is available where used
    (@~1453) and `transitions` before `_top_b_transitions` (@~1442).
  - In the consumer loop (@~1446-1448): `transition[:target_vertex]` →
    `transition.target_vertex`, `transition[:target_strand]` →
    `transition.target_strand`, `transition[:edge_data]` →
    `transition.edge_data`. Leave all scoring math (the `edge_w / total_out`
    divide, `log`, `next_score` add) byte-identical.

- [ ] **Step 4: Parse-check + run golden + lock tests.**

```bash
LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --startup-file=no -e 'Meta.parseall(read("src/viterbi-next.jl",String)); println("ok")'
LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e 'include("test/4_assembly/viterbi_decode_golden_test.jl")'   # path+score+diagnostics unchanged
for t in low_k_decode_gating_test reassembly_graph_reuse_test corrector_gc_between_batches_test parallel_soft_em_byte_identity_test; do
  LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e "include(\"test/4_assembly/$t.jl\")"
done
```

Expected: golden PASS (bit-identical decode); all lock/soft-EM tests PASS.
(`parallel_soft_em_byte_identity` self-hoists to `-t4`.) For
`batched_viterbi_kernel_test` use the stacked KernelAbstractions env (temp env
with `Pkg.add("KernelAbstractions")` + `JULIA_LOAD_PATH="@:<tmp>:@stdlib"`).

- [ ] **Step 5: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc add src/viterbi-next.jl test/4_assembly/viterbi_decode_golden_test.jl
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc commit -m "perf: fused typed transitions in the substitution hot loop (opt2, byte-identical)

td-cppm / td-jbjd opt2"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc push
```

---

### Task 3: Wire the fused typed path into the memoized indel successor batch

**Files:**

- Modify: `src/viterbi-next.jl` (`_indel_decode_successors!` ~1809-1845)

**Interfaces:** Consumes `_valid_transitions_and_total`, `_Transition` field
access.

- [ ] **Step 1: Read the current ~1809-1845 region verbatim.** It calls
      `_get_valid_transitions` (@~1820) then `_total_outgoing_weight` (@~1830)
      then consumes `transition[:target_vertex/:target_strand/:edge_data]`
      (@~1837-1841) while pushing to a `successors` batch.

- [ ] **Step 2: Wire the fused call + convert the consumer** (Python-bypass).
      Replace the paired `_get_valid_transitions` + `_total_outgoing_weight`
      with
      `transitions, total_out = Rhizomorph._valid_transitions_and_total(graph, vertex, strand)`.
      Preserve the `isempty(transitions)` early return and the
      `!isfinite(total_out) || total_out <= 0.0` guard (compute `total_out` from
      the fused return — note the fused `total_out` is already
      `max(total,_KSP_MIN_WEIGHT)`, so `total_out <= 0.0` is only reachable via
      the non-finite check; this matches the current `_total_outgoing_weight`
      return semantics — verify the guard still behaves identically). Convert
      `transition[:target_vertex]` → `transition.target_vertex`,
      `[:target_strand]` → `.target_strand`, `[:edge_data]` → `.edge_data`.

- [ ] **Step 3: Parse-check + re-run golden + lock tests + any indel-decode
      test.**

```bash
LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --startup-file=no -e 'Meta.parseall(read("src/viterbi-next.jl",String)); println("ok")'
# indel path tests + the standard locks
for t in viterbi_indel_decode_test indel_frontier_probe_test low_k_decode_gating_test reassembly_graph_reuse_test parallel_soft_em_byte_identity_test; do
  LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e "include(\"test/4_assembly/$t.jl\")" 2>&1 | tail -3 || echo "MISSING/FAIL $t"
done
```

Expected: all PASS. (If a named indel test file doesn't exist, locate the
indel-decode tests via `ls test/4_assembly/ | grep -i indel` and run those
instead — report which ran.) Byte-identity: the memoized batch's successor
order + weights are unchanged.

- [ ] **Step 4: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc add src/viterbi-next.jl
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc commit -m "perf: fused typed transitions in the memoized indel successor batch (opt2)

td-cppm / td-jbjd opt2"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc push
```

---

### Task 4: Hoist in-loop diagnostics to `Int` locals

**Files:**

- Modify: `src/viterbi-next.jl` (`_viterbi_correct_observation` diagnostics:
  init ~1326-1347, pre-seed ~1379-1381, in-loop mutations ~1423-1537)

**Interfaces:** none new.

- [ ] **Step 1: Read the diagnostics init + all in-loop mutation sites
      verbatim.** Enumerate every `diagnostics[:key]` write in the depth loop.

- [ ] **Step 2: Hoist to locals** (Python-bypass). Before the depth loop,
      declare `Int` locals seeded from the pre-loop values (e.g.
      `expanded_states = 0`, `skipped_transitions = 0`, `successor_bounded = 0`,
      `generated_states = 0`, `beam_pruned = 0`, `margin_pruned = 0`; and for
      the per-depth ones seed from the pre-seed @1379-1381:
      `cumulative_retained_states = diagnostics[:cumulative_retained_states]`,
      `max_retained_states = diagnostics[:max_retained_states]`, plus
      `retained_states = diagnostics[:retained_states]`, `completed_steps = 0`).
      Replace each in-loop `diagnostics[:k] += 1` / `get(diagnostics,:k,0)+1` /
      `max(...)` / overwrite with the local. After the loop, write each local
      back: `diagnostics[:expanded_states] = expanded_states`, etc. Preserve
      semantics exactly (running max, cumulative sum, last-write-wins for
      `retained_states`/`completed_steps`).

- [ ] **Step 3: Parse-check + golden + lock tests.**

```bash
LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --startup-file=no -e 'Meta.parseall(read("src/viterbi-next.jl",String)); println("ok")'
LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e 'include("test/4_assembly/viterbi_decode_golden_test.jl")'   # diagnostics Dict still bit-identical
for t in low_k_decode_gating_test reassembly_graph_reuse_test parallel_soft_em_byte_identity_test; do
  LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e "include(\"test/4_assembly/$t.jl\")"
done
```

Expected: golden PASS (the diagnostics Dict values are identical — the golden
test's diagnostics assertions are the lock here); locks PASS.

- [ ] **Step 4: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc add src/viterbi-next.jl
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc commit -m "perf: hoist hot-loop diagnostics counters to Int locals (opt2)

td-cppm / td-jbjd opt2"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc push
```

---

### Task 5: Allocation-reduction benchmark

**Files:**

- Create: `benchmarking/results/opt2_alloc_hotloop.md`

**Interfaces:** none.

- [ ] **Step 1: Measure `@allocated` on the decode** before/after — since
      "before" is already merged on `origin/master`, measure the CURRENT
      (post-opt2) allocation and compare to a run on `origin/master`. Practical
      form: run a fixed corrector decode under `@allocated` (or the accuracy
      benchmark with a bytes report) on this branch AND on a checkout of
      `origin/master`, record both.

```bash
# On this branch:
LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e '<build fixed graph+reads; @allocated Mycelia.Rhizomorph._viterbi_correct_observation(...) ; print bytes>'
```

- [ ] **Step 2: Record honestly** in
      `benchmarking/results/opt2_alloc_hotloop.md`: allocation bytes before
      (master) vs after (opt2), the reduction, the fixture, and the caveat that
      this is a micro-measurement (directional). Confirm decode output
      (path/score) identical across the two — the allocation win must not have
      changed output.

- [ ] **Step 3: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc add benchmarking/results/opt2_alloc_hotloop.md
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc commit -m "bench: record opt2 hot-loop allocation reduction

td-cppm / td-jbjd opt2"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt2-alloc push
```

---

### Task 6: PR + review gate

- [ ] **Step 1: Open PR** `cp/opt2-alloc → master` (title
      `perf: typed transitions + fused edge-scan in the Viterbi hot loop (byte-identical)`;
      body: summary + byte-identity evidence + allocation numbers + test plan;
      no internal IDs).
- [ ] **Step 2: Run `/comprehensive-pr-review`** on the head SHA; fold in remote
      CodeRabbit + Codex; resolve all Critical + Convergent findings; re-run the
      local review after any fix commit.
- [ ] **Step 3: Watch CI** (full suite incl. 4-thread job) to green. NOTE:
      adding the two test files must not trip a `fieldcount`/completeness guard
      (opt2 adds no struct field to a guarded struct, but re-verify the full
      suite).
- [ ] **Step 4: STOP at the merge gate** — merge is a human action (`--admin`
      bypass, user-run, as opt5/opt1). Do not self-merge.
- [ ] **Step 5: On merge:** close `td-cppm`, update epic `td-jbjd` (opt2 done;
      next opt4 `td-cvd8` approximate / opt3 `td-h4o8`), remove the `opt2-alloc`
      worktree + branch, run `/sync`.

## Notes for the implementer

- **Every `src/`/`test/` `.jl` edit goes through a Bash/Python in-place script**
  — never the Edit/Write tools (churning formatter hook).
  `git -C <wt> diff --stat` after each to confirm a minimal diff.
- **Read the current call region verbatim immediately before editing** each hot
  site — line numbers drift as tasks land; locate by content
  (`_get_valid_transitions`, `_total_outgoing_weight`, `_top_b_transitions`,
  `transition[:`).
- **Do NOT touch** `_get_valid_transitions` / `_maybe_push_transition!` /
  `_total_outgoing_weight` in `path-finding.jl` (Dict form retained for cold
  consumers) — only ADD the struct + fused fn.
- **The golden test (Task 2) is the primary byte-identity lock for the decode**
  (path + score + diagnostics); the equivalence test (Task 1) locks the fused
  function; the corrector lock tests + opt1 soft-EM test lock the end-to-end
  corrector.
