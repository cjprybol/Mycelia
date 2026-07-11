# Hybrid OLC: handing the accurized graph to an established layout assembler

- **Type:** Design draft (Stage-2 deconvolution route 1). Implementation is
  handed off to the `src/`-owning session; this document is the plan, not the
  code.
- **Date:** 2026-07-11
- **Scope:** `src/rhizomorph/assembly.jl` (entry point + config + iterative-corrector
  bridge), the external-assembler wrappers in `src/assembly.jl`, GFA export in
  `src/rhizomorph/algorithms/io.jl`, and the benchmark harness under
  `benchmarking/`.
- **Status:** Proposed. Grounded in the current `master` foundation (verified
  read-only). No `src/` changes are made by this document.
- **Companion:** `docs/design/2026-07-graph-as-hmm-corrector-methods.md` (Stage 1,
  the corrector this route consumes). The Rhizomorph manuscript's Methods §"Stage 2
  deconvolution" describes both native traversal and this reuse route.

## Thesis

Rhizomorph's Stage 1 (the graph-as-HMM corrector) produces **accurized reads**.
Stage 2 must turn those into a contiguous assembly. The native route
(Viterbi/`K`-shortest-path traversal) is one option; this document specifies the
**reuse route** — hand the accurized output to an established
overlap-layout-consensus (OLC) / long-read assembler so contiguity is inherited
from a tool already tuned for it, while the per-base accuracy gain is contributed
by Stage 1. The design goal from the manuscript is that accuracy and contiguity
need not trade off: correction happens first, layout second.

The key observation is that **most of this route already exists as foundation**.
The corrector already materializes corrected reads; the external assemblers are
already wrapped and tested. The route is therefore a thin adapter plus one small
plumbing change, not new machinery.

## Foundation already in place (verified on `master`)

- **Corrected-reads bridge.** `_assemble_with_iterative_corrector(reads, config)`
  (`src/rhizomorph/assembly.jl:842`) writes reads to a temp FASTQ, runs
  `mycelia_iterative_assemble` (`src/iterative-assembly.jl:233`, which returns
  `:final_fastq_file` and the in-memory corrected `Vector{FASTX.FASTQ.Record}`),
  reads the corrected reads back (`~:900`/`:912`), then **re-assembles** them via
  `assemble_genome(corrected_reads; corrector=:none)` (`~:979`). The hybrid route
  diverges exactly at this last step.
- **OLC / long-read wrappers** (`src/assembly.jl`), uniform convention
  (keyword-only, `add_bioconda_env("<tool>")` at entry, take a `fastq` path,
  return a NamedTuple of output paths): `run_hifiasm` (`:948`, +
  `hifiasm_primary_contigs` `:1021`), `run_flye` (`:647`, emits
  `assembly_graph.gfa`), `run_metaflye` (`:744`), `run_canu` (`:854`),
  `run_metamdbg` (`:2713`), `run_autocycler` (`src/autocycler.jl:58`). Tested
  under `test/4_assembly/` (gated by `MYCELIA_RUN_EXTERNAL=true`).
- **GFA export** (for route b): `Mycelia.write_gfa(result, path)`
  (`src/rhizomorph/assembly.jl:416`) / `Rhizomorph.write_gfa_next`
  (`src/rhizomorph/algorithms/io.jl:51`). `AssemblyResult` carries `graph`,
  `simplified_graph`, `gfa_compatible::Bool`.
- **Entry point + config.** `assemble_genome(reads, config)`
  (`src/rhizomorph/assembly.jl:648`) dispatches on `corrector`; `AssemblyConfig`
  keyword constructor (`:182`) with validation (`~:256`).
- **Tool provisioning.** `add_bioconda_env("<tool>")` (`src/bioconda.jl:780`) —
  adding a new assembler is a bioconda recipe name.
- **Benchmark harness.** `benchmarking/phix174_assembler_comparison.jl`,
  `mode_comparison.jl`, `real_data_corrector_validation.jl`, and the manifest
  `rhizomorph_benchmark_manifest.toml`.

## The one real gap

The corrected FASTQ produced inside `_assemble_with_iterative_corrector` lives in
an `mktempdir` that is deleted in the function's `finally` block, and no public
accessor returns just the corrected reads (post-`td-zru6`, the function always
re-assembles internally). The hybrid route needs the corrected reads to **outlive
that cleanup** so an external wrapper can consume them. This is the single genuine
code change; everything else is composition of existing functions.

## Route (a): corrected-reads bridge (primary — "possible with our foundation")

### Control flow

1. Run Stage-1 correction exactly as today, up to the point where
   `corrected_reads` / `corrected_fastq` are in hand
   (`src/rhizomorph/assembly.jl:~900–912`).
2. **Persist** the corrected reads to a caller-visible path (see plumbing below)
   instead of letting the temp dir delete them.
3. Call the selected external wrapper on that FASTQ, e.g.
   `run_hifiasm(fastq = corrected_fastq)` or `run_flye(fastq = corrected_fastq,
   read_type = ...)`.
4. Read the wrapper's primary contigs (`hifiasm_primary_contigs(result)` /
   the wrapper's `graph`/contig path) and wrap them into an `AssemblyResult` (or a
   thin hybrid result type) so downstream code and metrics are uniform.

### Config surface

Add an explicit layout selector rather than overloading `corrector`:

```julia
AssemblyConfig(;
    ...,
    layout::Symbol   = :native,    # :native (current behavior) | :olc
    olc_tool::Symbol = :hifiasm,   # :hifiasm | :flye | :canu | :metaflye | ...
    olc_options      = (;),        # pass-through NamedTuple to the wrapper
)
```

Validate `layout ∈ (:native, :olc)` and `olc_tool ∈ (:hifiasm, :flye, :canu,
:metaflye, :metamdbg, :autocycler)` near the existing validation block
(`~:256`). `layout = :olc` requires `corrector = :iterative` (an OLC layout with
no correction is just the plain external assembler, already reachable via the
wrapper directly). Branch in `assemble_genome(reads, config)` (`:648`): when
`layout == :olc`, dispatch to a new `_assemble_with_hybrid_olc(reads, config)`
that reuses the Stage-1 half of `_assemble_with_iterative_corrector` and replaces
the internal re-assembly with the wrapper call.

### Plumbing change (the gap)

Refactor `_assemble_with_iterative_corrector` so the "materialize corrected
reads" half is a reusable helper (e.g. `_run_stage1_correction(reads, config) ->
(corrected_reads, corrected_fastq_path)`) that writes the corrected FASTQ to a
persistent location owned by the caller (config `output_dir`, or a returned
tempfile the caller is responsible for). Both the existing native re-assembly and
the new hybrid route call this helper; only their tail differs. This keeps the
corrected-FASTQ lifetime explicit and testable, and avoids duplicating the
correction call.

### Read-type mapping

`run_flye`/`run_canu` need a `read_type`/chemistry flag; `run_hifiasm` assumes
HiFi. Derive the wrapper's read-type argument from the same input-type detection
`_auto_configure_assembly` already performs, or require the caller to pass it via
`olc_options`. For the first validation (Illumina short reads → hifiasm/flye) the
mapping is trivial and can be hard-wired; generalize later.

## Route (b): GFA-substrate bridge (stretch)

Instead of round-tripping through corrected reads, emit the accurized graph as
GFA (`Mycelia.write_gfa(result, path)`) and hand it to an engine that accepts a
graph as its layout substrate (Flye-family `assembly_graph.gfa` consumers). This
avoids re-deriving overlaps from reads but is higher-risk: (1) external tools vary
in whether they accept an externally supplied GFA as input vs. only emit one;
(2) preserving the Stage-1 accuracy gain through the tool's own graph-cleaning is
unproven. Stage this behind Route (a); it is not required for a first hybrid arm.

## Validation plan (handed to the `benchmarking/`-owning session)

- Extend `benchmarking/phix174_assembler_comparison.jl` (or
  `real_data_corrector_validation.jl`) with a **`hybrid`/`olc` arm** alongside the
  existing `naive` and `scalable` arms, on phiX174 + lambda (same ART HS25 50×
  simulated reads, MUMmer `dnadiff` metrics), so the three arms are directly
  comparable in one committed CSV.
- Register the arm in `rhizomorph_benchmark_manifest.toml`.
- The manuscript's pre-registered **H5** (state-of-the-art contiguity) is the
  eventual home for the cross-assembler comparison; the interim CSV is engineering
  validation only and must be labelled as such, exactly like the existing
  `real_data_corrector_validation` result.

## Missing wrappers (optional, Conda-only)

`run_verkko`, `run_raven`, `run_miniasm`, `run_wtdbg2`, `run_shasta`,
`run_nextdenovo` do not exist yet (names appear only as strings in
`benchmarking/03_assembly_benchmark.jl`). None is required for a first hybrid arm
(hifiasm + flye suffice). Each is a small addition following the
`add_bioconda_env` + `CONDA_RUNNER run -n <env>` pattern; file separately if a
benchmark tier needs them.

## Risks and open questions

- **Correction→OLC accuracy transfer.** Whether feeding corrected reads to an OLC
  assembler retains the Stage-1 per-base gain (vs. the assembler's own correction
  undoing or duplicating it) is exactly what the benchmark must measure; do not
  assert parity in advance (mirrors the manuscript's H5 hedge).
- **Double correction.** hifiasm/canu perform their own error correction; running
  them on already-corrected reads may be redundant or counterproductive. The
  validation arm should include a "corrected-reads → assembler-with-its-own-
  correction-disabled" variant where the tool supports it.
- **Cost.** Route (a) pays for correction + a full external assembly. Report
  wall-clock alongside quality (the corrector already costs 5–10× naive on the
  committed phiX174/lambda run).

## Work breakdown (handoff beads)

Filed as `[B]` beads in the todo control plane (src/+benchmarking work routes to
the `src/`-owning session; this doc is the spec):

1. **Plumbing:** refactor `_assemble_with_iterative_corrector` to expose a
   persistent corrected-FASTQ helper (`_run_stage1_correction`), with a unit test
   asserting the corrected FASTQ exists and round-trips.
2. **Route (a) wiring:** add `layout`/`olc_tool` config + validation +
   `_assemble_with_hybrid_olc`, dispatching to `run_hifiasm`/`run_flye`; gated
   real-run test under `test/4_assembly/`.
3. **Validation arm:** add the `hybrid` arm to the phiX174/lambda harness; commit
   the comparison CSV.
4. **Manuscript update (docs, this session's domain):** once (2)–(3) land, tighten
   the Methods §"Stage 2 deconvolution" wording from "designed" to
   "implemented, benchmark pending," and add the hybrid arm to the interim
   validation table in Results.
5. **Stretch:** Route (b) GFA-substrate bridge; missing long-read wrappers as
   needed.
