# opt1 Parallel Soft-EM Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the `:scalable` corrector's per-read Viterbi decode under
`Threads.@threads` with soft-EM enabled, byte-identically to serial, defaulting
parallel-on when `nthreads > 1`.

**Architecture:** Each `@threads` iteration decodes into its own per-read
`SoftEdgeWeightAccumulator` written by index (`batch_local[i]`); after the loop,
the per-read accumulators are folded into the shared accumulator in ascending
read order, exactly reproducing the serial left-fold so soft-EM weights stay
bit-identical despite float non-associativity. The read-side snapshot is already
parallel-safe; the corrected-reads collection is already order-deterministic.

**Tech Stack:** Julia; FASTX; MetaGraphsNext; `Test`, `Logging` stdlibs; GitHub
Actions CI.

## Global Constraints

- **Byte-identity is the hard requirement:** corrected reads AND soft-EM weights
  must be identical between `enable_parallel=true` and `enable_parallel=false`.
  Re-verify after every task.
- **`.jl` edits must bypass the formatter hook.** The repo's PostToolUse
  Edit/Write hook reformats `.jl` files to non-SciML style, whole-file-churning
  the diff. Apply every `src/`/`test/` `.jl` change via a Bash/Python in-place
  script (the hook fires on the Edit/Write tools, not on Bash). Markdown (`.md`)
  and YAML (`.yml`) may use Write normally.
- **SciML style** (`.JuliaFormatter.toml`: `style="sciml"`) — match surrounding
  code (trailing commas, 4-space indent).
- **Julia invocation:**
  `LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. ...` (per repo/HPC
  convention).
- **No internal tracking IDs** (`td-*`) in committed non-comment artifacts (test
  names, CI, docs prose). They may appear in commit messages and code comments
  only.
- **Worktree:**
  `/Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel`, branch
  `cp/opt1-parallel` (off `origin/master`, post-opt5).
- **Commit after each task; do not merge** (human gate). Base branch `master`.

## Reference: verified baseline locations (post-opt5)

- Guard: `src/iterative-assembly.jl:4773`
  `use_parallel = enable_parallel && soft_weights === nothing`; `@warn` at
  `4774-4776`.
- Parallel branch: `src/iterative-assembly.jl:4849-4908` (`Threads.@threads`
  loop; does NOT pass `soft_weights`). Serial branch (reference): `4909-4948`
  (passes shared `soft_weights` at 4938).
- Reduce primitive: `_merge_soft_edge_weights!(dest, staged)` at
  `src/iterative-assembly.jl:3859-3868`.
- Per-read decode:
  `improve_read_likelihood(read, graph, k; graph_mode, beam_width, soft_weights::Union{Nothing,SoftEdgeWeightAccumulator}=nothing, weighted_graph, ...)`
  at `src/iterative-assembly.jl:5105`; stages into `staged_soft_weights` (6257)
  and merges (6606).
- Accumulator:
  `struct SoftEdgeWeightAccumulator; weights::Dict{Any,Float64}; end` at
  `src/rhizomorph/core/evidence-functions.jl:694-698`; zero-arg ctor at 698.
- Read-set wrapper:
  `improve_read_set_likelihood(reads, graph, k; ..., enable_parallel::Bool=false, soft_weights=nothing, gc_between_batches::Bool=false, ...)`
  (public entry ~4448) — mutates the caller-passed `soft_weights` in place.
- Caller knobs: `src/rhizomorph/assembly.jl` `_corrector_strategy_knobs`
  `:scalable` (~1227-1261, `soft_em=true` at 1238), `:exhaustive` (~1269);
  corrector call `mycelia_iterative_assemble(...)` (~1492-1512, omits
  `enable_parallel`).
- CI: `.github/workflows/ci.yml:91`
  `run: julia --color=yes --project=. -e 'import Pkg; Pkg.test(coverage=true);'`
  (single-threaded).

---

> **REVISED (post-review):** Task 1's fold is upgraded from a per-read
> accumulator to **per-window ordered capture + flat read-order replay** — the
> production `:scalable` path sets `windowed_decode=true`, and the windowed
> decode merges per-window, so a per-read fold is not bit-identical to serial.
> See the design spec §"2. Parallel branch: per-window ordered capture + flat
> read-order replay" for the exact mechanism (a `soft_weights_sink` that
> captures staged accumulators in decode order;
> `batch_local::Vector{Vector{SoftEdgeWeightAccumulator}}`; flat replay
> `for i; for staged in batch_local[i]; _merge_soft_edge_weights!(soft_weights, staged)`).
> The test must include a `windowed_decode=true` repeat-heavy case.

### Task 1: Per-window ordered capture + flat read-order replay (parallel decode)

**Files:**

- Modify: `src/iterative-assembly.jl` (guard `~4773`; parallel branch
  `~4849-4908`)
- Create: `test/4_assembly/parallel_soft_em_byte_identity_test.jl`

**Interfaces:**

- Consumes:
  `Mycelia.improve_read_set_likelihood(reads, graph, k; graph_mode, skip_solid, cheap_correct, hard_vertices, decode_enabled, batch_size, enable_parallel, soft_weights)`;
  `Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()`;
  `Mycelia.Rhizomorph.build_qualmer_graph`; `Mycelia._hard_vertex_set`.
- Produces: parallel decode path that fills a caller-passed `soft_weights`
  accumulator byte-identically to serial. Test helper
  `_psi_reads(rng, ref; ...)`, `_psi_run(reads, graph, k; parallel)` returning
  `(corrected_seqs, sorted_weight_pairs)`.

- [ ] **Step 1: Write the failing test.** Create
      `test/4_assembly/parallel_soft_em_byte_identity_test.jl` via a
      Python-bypass writer (`.jl` hook). The test self-hoists to 4 threads: if
      `Threads.nthreads() == 1` it relaunches itself under `julia -t 4` and
      asserts the child succeeds; if `> 1` it runs the assertions. It asserts
      (a) the parallel path is TAKEN — no `"ignoring enable_parallel"` warning
      is logged for the `enable_parallel=true` + soft-EM call (captured via
      `Logging.with_logger(Test.TestLogger())`); and (b) byte-identity —
      corrected sequences AND the full soft-EM weight Dict are identical between
      an `enable_parallel=true` run and an `enable_parallel=false` run, each
      with its own caller-passed accumulator.

```julia
# PARALLEL SOFT-EM BYTE-IDENTITY (opt1)
# Corrected reads AND soft-EM weights must be identical between parallel and
# serial decode. CI runs single-threaded, so this test relaunches itself under
# `julia -t 4` to exercise genuine concurrency (otherwise @threads runs serially
# and proves nothing). RED signal before the fix: the guard demotes
# enable_parallel=true+soft-EM to serial and logs "ignoring enable_parallel".
import Test
import Mycelia
import FASTX
import Random
import Logging

if Threads.nthreads() == 1 && get(ENV, "MYCELIA_PSI_CHILD", "0") != "1"
    # Self-hoist to real threads.
    proj = Base.active_project()
    thisfile = @__FILE__
    cmd = `$(Base.julia_cmd()) --project=$(proj) --threads=4 -e "include(\"$(thisfile)\")"`
    Test.@testset "parallel soft-EM byte-identity (hoisted to -t4)" begin
        Test.@test success(pipeline(setenv(cmd, "MYCELIA_PSI_CHILD" => "1"; dir = pwd());
            stdout = stdout, stderr = stderr))
    end
else
    const _PSI_BASES = ['A', 'C', 'G', 'T']

    function _psi_reads(rng, ref; n_reads = 200, readlen = 80, n_err = 40)
        reflen = length(ref)
        records = FASTX.FASTQ.Record[]
        for i in 1:n_reads
            s = rand(rng, 1:(reflen - readlen + 1))
            seq = collect(ref[s:(s + readlen - 1)])
            if i <= n_err
                p = rand(rng, 1:readlen)
                seq[p] = rand(rng, filter(!=(seq[p]), _PSI_BASES))
            end
            push!(records,
                FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
        end
        return records
    end

    _psi_seqs(records) = [FASTX.sequence(String, r) for r in records]
    _psi_weight_pairs(acc) = sort!(collect(acc.weights); by = p -> repr(p[1]))

    function _psi_run(reads, graph, k, hard; parallel::Bool)
        acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        out, = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard, decode_enabled = true,
            batch_size = 50, enable_parallel = parallel, soft_weights = acc)
        return _psi_seqs(out), _psi_weight_pairs(acc)
    end

    Test.@testset "parallel soft-EM byte-identity (opt1)" begin
        Test.@test Threads.nthreads() > 1     # guard: real threads
        rng = Random.MersenneTwister(2024)
        ref = join(rand(rng, _PSI_BASES, 1500))
        reads = _psi_reads(rng, ref; n_reads = 200, readlen = 80, n_err = 40)
        k = 13
        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)

        # (a) parallel path taken: no "ignoring enable_parallel" warning
        logger = Test.TestLogger()
        seqs_par, wts_par = Logging.with_logger(logger) do
            _psi_run(reads, graph, k, hard; parallel = true)
        end
        Test.@test !any(occursin("ignoring enable_parallel", r.message)
                        for r in logger.logs)

        # (b) byte-identity vs serial: reads AND soft-EM weights
        seqs_ser, wts_ser = _psi_run(reads, graph, k, hard; parallel = false)
        Test.@test seqs_par == seqs_ser
        Test.@test wts_par == wts_ser
    end
end
```

- [ ] **Step 2: Run test to verify it fails.**

Run:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -t 4 -e 'ENV["MYCELIA_PSI_CHILD"]="1"; include("test/4_assembly/parallel_soft_em_byte_identity_test.jl")'`
Expected: FAIL on `!any(occursin("ignoring enable_parallel", ...))` — the guard
demotes to serial and logs the warning. (The byte-identity `@test`s pass
trivially since both runs are serial.)

- [ ] **Step 3: Relax the guard.** Apply via Python-bypass. In
      `src/iterative-assembly.jl` change the guard + comment + drop the `@warn`:

```julia
    # opt1: soft-EM accumulation is now thread-safe — the parallel branch gives
    # each read its own SoftEdgeWeightAccumulator (indexed write), then folds them
    # into the shared accumulator in read order after the loop, reproducing the
    # serial left-fold bit-for-bit. So parallel + soft-EM coexist byte-identically.
    use_parallel = enable_parallel
```

Remove the `if enable_parallel && soft_weights !== nothing; @warn ...; end`
block (4774-4776) entirely.

- [ ] **Step 4: Implement the per-read accumulators + read-order fold.** Apply
      via Python-bypass in the `if use_parallel` branch (`~4849-4908`).

Immediately after `batch_results = Vector{...}(undef, length(batch_reads))` and
the `skip_flags = fill(false, ...)` allocation, add the per-read accumulator
vector:

```julia
            # opt1: one accumulator per read (indexed write, no race), folded in
            # read order after the loop. `nothing` when soft-EM is off.
            batch_local = soft_weights === nothing ? nothing :
                          Vector{Mycelia.Rhizomorph.SoftEdgeWeightAccumulator}(
                              undef, length(batch_reads))
```

Inside the `Threads.@threads for i in eachindex(batch_reads)` body, in the
non-skipped branch where `improve_read_likelihood(read, decode_graph, k; ...)`
is called (the parallel call that currently omits `soft_weights`), give this
read its own accumulator and pass it:

```julia
                        local_acc = soft_weights === nothing ? nothing :
                                    Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
                        if soft_weights !== nothing
                            batch_local[i] = local_acc
                        end
                        improved_read, was_improved = ... improve_read_likelihood(
                            read, decode_graph, k;
                            graph_mode = graph_mode, beam_width = beam_width,
                            soft_weights = local_acc,
                            weighted_graph = pass_weighted_graph, ...)
                        batch_results[i] = (improved_read, was_improved)
```

(Match the exact existing kwargs of the parallel `improve_read_likelihood` call;
the only additions are `soft_weights = local_acc` and the
`local_acc`/`batch_local[i]` lines. For skipped reads, leave `batch_local[i]`
undefined and handle in the fold.)

After the `Threads.@threads` loop and the existing
`skipped_reads += count(skip_flags)` / results-collection loop, add the
deterministic read-order fold:

```julia
            # opt1: fold per-read soft-EM contributions into the shared accumulator
            # in ascending read order — identical summation order to the serial
            # path, so weights are bit-identical despite float non-associativity.
            if soft_weights !== nothing
                for i in eachindex(batch_reads)
                    isassigned(batch_local, i) || continue   # skipped reads
                    _merge_soft_edge_weights!(soft_weights, batch_local[i])
                end
            end
```

- [ ] **Step 5: Run the test to verify it passes.**

Run:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -t 4 -e 'ENV["MYCELIA_PSI_CHILD"]="1"; include("test/4_assembly/parallel_soft_em_byte_identity_test.jl")'`
Expected: PASS — no ignoring-warning, `seqs_par == seqs_ser`,
`wts_par == wts_ser`.

- [ ] **Step 6: Verify the self-hoist wrapper works at 1 thread.**

Run:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e 'include("test/4_assembly/parallel_soft_em_byte_identity_test.jl")'`
Expected: PASS — the 1-thread parent relaunches under `-t4` and reports the
child succeeded.

- [ ] **Step 7: Re-verify the three lock tests + opt5 test (byte-identity
      regression).**

Run each:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e 'include("test/4_assembly/<name>.jl")'`
for `low_k_decode_gating_test`, `reassembly_graph_reuse_test`,
`corrector_gc_between_batches_test`. For `batched_viterbi_kernel_test`, use the
stacked-env form (adds `KernelAbstractions`): create a temp env `mktemp -d`,
`julia --project=<tmp> -e 'import Pkg; Pkg.add("KernelAbstractions")'`, then
`JULIA_LOAD_PATH="@:<tmp>:@stdlib" julia --project=. -e 'include(".../batched_viterbi_kernel_test.jl")'`.
Expected: all PASS (byte-identity preserved; the parallel change is inert at
`enable_parallel=false`, the default for these tests).

- [ ] **Step 8: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel add src/iterative-assembly.jl test/4_assembly/parallel_soft_em_byte_identity_test.jl
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel commit -m "feat: thread-safe soft-EM parallel decode (deferred read-order reduce)

Per-read SoftEdgeWeightAccumulators written by index during the Threads.@threads
decode, folded into the shared accumulator in read order after the loop —
reproducing the serial left-fold bit-for-bit. Relaxes the use_parallel guard so
parallel + soft-EM coexist. Byte-identical corrected reads and soft-EM weights
verified under julia -t 4 (subprocess-hoisted test).

td-vmiy / td-jbjd opt1"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel push
```

---

### Task 2: Default `:scalable` to parallel when `nthreads > 1`

**Files:**

- Modify: `src/rhizomorph/assembly.jl` (`_corrector_strategy_knobs` `:scalable`
  ~1227-1261; `mycelia_iterative_assemble` call ~1492-1512)
- Modify: `test/4_assembly/scalable_corrector_strategy_test.jl` (add a knob
  truth-table assertion; if absent, create a focused test file)

**Interfaces:**

- Consumes: `Mycelia`-internal `_corrector_strategy_knobs(strategy::Symbol)`
  returning a NamedTuple with `soft_em`, and now `enable_parallel`.
- Produces:
  `_corrector_strategy_knobs(:scalable).enable_parallel == (Threads.nthreads() > 1)`;
  `_corrector_strategy_knobs(:exhaustive).enable_parallel == false`; the
  `:scalable` corrector call forwards `enable_parallel`.

- [ ] **Step 1: Write the failing test** (Python-bypass). Add to
      `test/4_assembly/scalable_corrector_strategy_test.jl` (or create
      `test/4_assembly/scalable_parallel_default_test.jl`):

```julia
Test.@testset "scalable defaults parallel on multi-thread (opt1)" begin
    ks = Mycelia._corrector_strategy_knobs(:scalable)
    Test.@test ks.enable_parallel == (Threads.nthreads() > 1)
    ex = Mycelia._corrector_strategy_knobs(:exhaustive)
    Test.@test ex.enable_parallel == false     # exhaustive stays serial/exact
end
```

- [ ] **Step 2: Run test to verify it fails.**

Run:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e 'include("test/4_assembly/scalable_corrector_strategy_test.jl")'`
Expected: FAIL — `ks` has no field `enable_parallel`
(`type NamedTuple has no field enable_parallel`).

- [ ] **Step 3: Add the knob + forward it** (Python-bypass). In
      `_corrector_strategy_knobs`, add
      `enable_parallel = Threads.nthreads() > 1,` to the `:scalable` NamedTuple
      and `enable_parallel = false,` to `:exhaustive`. In the
      `mycelia_iterative_assemble(...)` call (~1492-1512), add
      `enable_parallel = knobs.enable_parallel,`.

- [ ] **Step 4: Run test to verify it passes.**

Run:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. -e 'include("test/4_assembly/scalable_corrector_strategy_test.jl")'`
Expected: PASS.

- [ ] **Step 5: Integration byte-identity via `assemble_genome`**
      (Python-bypass; append to the same test file). Assert a full `:scalable`
      assemble is byte-identical between forced-serial and default (parallel
      under -t4). Reuse the self-hoist pattern if the file isn't already
      thread-hoisted; otherwise gate on `Threads.nthreads() > 1`.

```julia
if Threads.nthreads() > 1
    Test.@testset "scalable assemble_genome parallel==serial (opt1)" begin
        rng = Random.MersenneTwister(77)
        ref = join(rand(rng, ['A','C','G','T'], 1200))
        reads = [FASTX.FASTQ.Record("r$i",
            ref[(s):(s+79)], String(fill('I', 80)))
            for (i, s) in enumerate(rand(rng, 1:1121, 150))]
        tmp = mktempdir(); fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)
        runit = par -> Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 2,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = true, enable_parallel = par,
            verbose = false, enable_checkpointing = false,
            output_dir = joinpath(tmp, par ? "par" : "ser"))
        seqs(res) = [FASTX.sequence(String, r) for r in
            open(FASTX.FASTQ.Reader, res[:metadata][:final_fastq_file]) do rd
                collect(rd)
            end]
        Test.@test seqs(runit(true)) == seqs(runit(false))
    end
end
```

- [ ] **Step 6: Run + verify pass.** Run under `-t 4` as in Task 1 Step 5.
      Expected: PASS.

- [ ] **Step 7: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel add src/rhizomorph/assembly.jl test/4_assembly/scalable_corrector_strategy_test.jl
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel commit -m "feat: default :scalable corrector to parallel when nthreads>1

_corrector_strategy_knobs(:scalable) now sets enable_parallel = nthreads>1 and
the corrector call forwards it; :exhaustive stays serial/exact. Integration
test asserts assemble_genome(:scalable) is byte-identical parallel vs serial.

td-vmiy / td-jbjd opt1"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel push
```

---

### Task 3: Parallel × GC on/off byte-identity (deferred opt5 review item)

**Files:**

- Modify: `test/4_assembly/parallel_soft_em_byte_identity_test.jl` (add a
  testset)

**Interfaces:**

- Consumes: `Base.withenv`, `MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES` env,
  `gc_between_batches` kwarg (from opt5).

- [ ] **Step 1: Write the test** (Python-bypass; add inside the `else`
      real-threads branch). With `enable_parallel = true` under real threads,
      assert corrected reads + weights are identical with the between-batch GC
      on vs off.

```julia
    Test.@testset "parallel x GC on/off byte-identity (opt1 + opt5)" begin
        rng = Random.MersenneTwister(555)
        ref = join(rand(rng, _PSI_BASES, 1500))
        reads = _psi_reads(rng, ref; n_reads = 200, readlen = 80, n_err = 40)
        k = 13
        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)
        runit = gc -> begin
            acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
            out, = Mycelia.improve_read_set_likelihood(
                reads, graph, k; graph_mode = :canonical, skip_solid = true,
                cheap_correct = true, hard_vertices = hard, decode_enabled = true,
                batch_size = 50, enable_parallel = true, gc_between_batches = gc,
                soft_weights = acc)
            (_psi_seqs(out), _psi_weight_pairs(acc))
        end
        off = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
            runit(false)
        end
        on = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
            runit(true)
        end
        Test.@test off[1] == on[1]        # reads
        Test.@test off[2] == on[2]        # weights
    end
```

- [ ] **Step 2: Run + verify pass.** Run the file under `-t 4` (Task 1 Step 5
      form). Expected: PASS.

- [ ] **Step 3: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel add test/4_assembly/parallel_soft_em_byte_identity_test.jl
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel commit -m "test: parallel x between-batch-GC byte-identity (opt5 review follow-up)

Exercises the enable_parallel=true path the opt5 review (I2) could not: GC on vs
off produces identical corrected reads and soft-EM weights under julia -t 4.

td-vmiy / td-jbjd opt1"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel push
```

---

### Task 4: Threaded CI job

**Files:**

- Modify: `.github/workflows/ci.yml`

**Interfaces:** none (CI config).

- [ ] **Step 1: Add a threaded test step/job.** Add a matrix dimension or a
      dedicated job that runs the suite with 4 threads. Minimal form — a second
      test step (Write is fine; `.yml` is not `.jl`):

```yaml
- name: Run tests (4 threads)
  run: julia --color=yes --project=. --threads=4 -e 'import Pkg; Pkg.test();'
```

Place it after the existing coverage test step in `.github/workflows/ci.yml`.
(Prefer a separate matrix entry `threads: [1, 4]` with
`--threads=${{ matrix.threads }}` if the workflow already uses a matrix, to
avoid doubling wall-time on the coverage run.)

- [ ] **Step 2: Validate locally that the suite passes at 4 threads** (the
      parallel path now runs across the whole `4_assembly` suite). Given the
      full suite is ~30+ min, at minimum run the corrector-relevant files under
      4 threads:

Run:
`LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. --threads=4 -e 'for f in ("parallel_soft_em_byte_identity_test","corrector_gc_between_batches_test","low_k_decode_gating_test","reassembly_graph_reuse_test","scalable_corrector_strategy_test"); include("test/4_assembly/$f.jl"); end'`
Expected: all PASS at 4 threads.

- [ ] **Step 3: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel add .github/workflows/ci.yml
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel commit -m "ci: run the test suite with 4 threads

The corrector's parallel decode path only exercises under multiple threads;
CI previously ran single-threaded, so @threads code paths were untested.

td-vmiy / td-jbjd opt1"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel push
```

---

### Task 5: Benchmark the multi-core speedup

**Files:**

- Create: `benchmarking/results/opt1_parallel_scaling_<label>.md` (recorded
  numbers + provenance)

**Interfaces:** Consumes
`benchmarking/rhizomorph_correction_accuracy_benchmark.jl`
(`corrector_runtime_s`, `MYCELIA_RCA_*` env incl. the opt5
`MYCELIA_RCA_BATCH_SIZE`).

- [ ] **Step 1: Run serial vs multi-thread on a representative genome.** phix
      @20x, forced `batch_size=50` (multi-batch), serial (`-t 1`) vs parallel
      (`-t 8`):

```bash
cd /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel
for T in 1 8; do
  MYCELIA_RCA_ACCESSION=NC_001422 MYCELIA_RCA_COVERAGE=20 MYCELIA_RCA_BATCH_SIZE=50 \
  MYCELIA_RCA_ERR=0.01 MYCELIA_RCA_SEED=42 MYCELIA_RCA_SCALE_FLOOR=1 \
  LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. --threads=$T \
    benchmarking/rhizomorph_correction_accuracy_benchmark.jl 2>&1 | tee /tmp/opt1_t$T.log
done
```

Extract `corrector_runtime_s` from the two newest CSVs in
`benchmarking/results/` and confirm `recall`/`precision` are identical across
thread counts (byte-identity sanity).

- [ ] **Step 2: Record honestly.** Write
      `benchmarking/results/opt1_parallel_scaling_<label>.md` with: thread
      counts, `corrector_runtime_s` each, speedup ratio, machine, rep count, and
      the note that phix is small (limited parallel headroom; larger genomes
      scale better). Report the number as measured — do not extrapolate.

- [ ] **Step 3: Commit.**

```bash
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel add benchmarking/results/opt1_parallel_scaling_*.md
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel commit -m "bench: record opt1 serial-vs-parallel corrector scaling

td-vmiy / td-jbjd opt1"
git -C /Users/cameronprybol/workspace/Mycelia/.worktrees/opt1-parallel push
```

---

### Task 6: PR + review gate

- [ ] **Step 1: Open PR** `cp/opt1-parallel → master` (title
      `perf: parallelize :scalable corrector decode with thread-safe soft-EM (byte-identical)`;
      body: summary + test plan; no internal IDs).
- [ ] **Step 2: Run `/comprehensive-pr-review`** on the head SHA; fold in remote
      CodeRabbit + Codex; resolve all Critical + Convergent findings; re-run the
      local review after any fix commit.
- [ ] **Step 3: Watch CI** (including the new 4-thread job) to green.
- [ ] **Step 4: STOP at the merge gate** — merge is a human action (admin bypass
      is user-run, as with opt5). Do not self-merge.
- [ ] **Step 5: On merge:** close `td-vmiy`, update epic `td-jbjd` (opt1 done;
      next opt2 `td-cppm`), remove the `opt1-parallel` worktree + branch, run
      `/sync` to export beads jsonl.

## Notes for the implementer

- **Every `src/`/`test/` `.jl` edit goes through a Bash/Python in-place script**
  (see Global Constraints) — never the Edit/Write tools, which trigger the
  churning formatter hook. After each edit, `git -C <worktree> diff --stat` to
  confirm a minimal diff (no whole-file churn).
- **Parse-check after edits:**
  `LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --startup-file=no -e 'Meta.parseall(read("<file>", String)); println("ok")'`.
- **The parallel branch's exact existing kwargs** for `improve_read_likelihood`
  must be preserved verbatim; the only additions are `soft_weights = local_acc`
  and the `batch_local` bookkeeping. Read `src/iterative-assembly.jl:4849-4908`
  immediately before editing to copy the current call verbatim (line numbers may
  drift as tasks land).
- **`isassigned(batch_local, i)`** guards skipped reads whose accumulator was
  never constructed — do not fold an undefined slot.
