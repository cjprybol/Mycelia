# Hybrid OLC: handing the accurized graph to an established layout assembler

- **Type:** Design and implementation note (Stage-2 deconvolution route 1).
- **Date:** 2026-07-11
- **Scope:** `src/rhizomorph/assembly.jl` (entry points, configs, and
  iterative-corrector bridge), `src/autocycler.jl`, the external-assembler
  wrappers in `src/assembly.jl`, GFA export in
  `src/rhizomorph/algorithms/io.jl`, and the benchmark harness under
  `benchmarking/`.
- **Status:** Route (a) and the td-06er multi-read-set workflows are
  implemented, with external-tool smoke tests gated.
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
  `:final_fastq_file` under `:metadata`, plus the in-memory corrected reads as
  `:final_assembly`), reads the corrected reads back (`~:900`/`:912`), then
  **re-assembles** them via `assemble_genome(corrected_reads; corrector=:none)`
  (`~:979`). The hybrid route diverges exactly at this last step.
- **External-assembler wrappers** (`src/assembly.jl`). The **long-read** wrappers
  `run_hifiasm` (`:948`, + `hifiasm_primary_contigs` `:1021`), `run_flye`
  (`:647`, emits `assembly_graph.gfa`), `run_metaflye` (`:744`), and `run_canu`
  (`:854`) share a convention: keyword-only, `add_bioconda_env("<tool>")` at
  entry, take a `fastq` path, return a NamedTuple of output paths. **This
  convention is not universal** — `run_metamdbg` (`:2713`) takes
  `hifi_reads`/`ont_reads`, and `run_autocycler` (`src/autocycler.jl:58`) takes
  `long_reads` and provisions via its own installer (not `add_bioconda_env`), so
  both need per-tool argument adapters rather than a `fastq=` call. **Short-read**
  wrappers for Illumina data exist separately: `run_megahit` (`:27`),
  `run_metaspades` (`:325`), `run_spades` (`:421`), `run_skesa` (`:514`). Tested
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
accessor returns just the corrected reads — the function always re-assembles
internally. The hybrid route needs the corrected reads to **outlive that cleanup**
so an external wrapper can consume them. This is the single genuine code change;
everything else is composition of existing functions.

There is also a **dead-code stub to reconcile**: `_assemble_hybrid_olc(observations,
config)` already exists at `src/rhizomorph/assembly.jl:1585` as an unreachable
placeholder (`@warn "Hybrid OLC not fully implemented"` → falls back to
`_assemble_kmer_graph`; never dispatched, since `strategy` validation admits only
`:scalable`/`:exhaustive`). The implementer should **repurpose or remove this
stub** rather than add a confusingly-named sibling — prefer wiring the real route
into the existing `_assemble_hybrid_olc` name.

## Route (a): corrected-reads bridge (primary — "possible with our foundation")

### Control flow

1. Run Stage-1 correction exactly as today, up to the point where
   `corrected_reads` / `corrected_fastq` are in hand
   (`src/rhizomorph/assembly.jl:~900–912`).
2. **Persist** the corrected reads to a caller-visible path (see plumbing below)
   instead of letting the temp dir delete them.
3. Call the **read-type-appropriate** wrapper on that FASTQ through a per-tool
   argument adapter, e.g. `run_megahit(fastq = corrected_fastq)` /
   `run_metaspades(...)` for short-read-corrected data, or
   `run_hifiasm(fastq = corrected_fastq)` / `run_flye(fastq = corrected_fastq,
   read_type = ...)` for long-read-corrected data. (Do not feed short reads to a
   long-read assembler — see Read-type mapping.)
4. Read the wrapper's primary contigs (`hifiasm_primary_contigs(result)` /
   the wrapper's `graph`/contig path) and wrap them into an `AssemblyResult` (or a
   thin hybrid result type) so downstream code and metrics are uniform.

### Config surface

Add an explicit layout selector rather than overloading `corrector`:

```julia
AssemblyConfig(;
    ...,
    layout::Symbol   = :native,    # :native (current behavior) | :olc
    olc_tool::Symbol = :auto,      # :auto (pick by read type) | short-read
                                   # :megahit|:metaspades | long-read
                                   # :hifiasm|:flye|:canu|:metaflye
    olc_options      = (;),        # pass-through NamedTuple to the wrapper
)
```

Validate `layout ∈ (:native, :olc)`, and validate `olc_tool` **against the
corrected read type** (see Read-type mapping): short-read layout assemblers
(`:megahit`, `:metaspades`) for Illumina-corrected reads, long-read assemblers
(`:hifiasm`, `:flye`, `:canu`, `:metaflye`) for long-read-corrected reads;
`:auto` picks by read type. `layout = :olc` requires `corrector = :iterative` (an
OLC layout with no correction is just the plain external assembler, already
reachable via the wrapper directly). Branch in `assemble_genome(reads, config)`
(`:648`): when `layout == :olc`, dispatch by **repurposing the existing dead
`_assemble_hybrid_olc` stub** (`:1585`) — reuse the Stage-1 half of
`_assemble_with_iterative_corrector` and replace the internal re-assembly with a
read-type-appropriate wrapper call routed through a per-tool argument adapter
(`fastq` / `hifi_reads` / `long_reads`).

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

The OLC tool must match the type of the corrected reads. `run_hifiasm` /
`run_flye` / `run_canu` / `run_metaflye` are **long-read** assemblers (hifiasm
assumes HiFi; flye/canu need a chemistry flag) and must **not** be handed Illumina
short reads — doing so errors or yields garbage. Short-read-corrected output must
route to `run_megahit` / `run_metaspades` instead. Derive the read type from the
same input-type detection `_auto_configure_assembly` already performs, and reject
(fail fast) an `olc_tool` incompatible with it rather than silently feeding a HiFi
assembler 150 bp reads. Because the committed validation corrects ART HS25
**Illumina** reads (see Validation plan), the first hybrid arm pairs the corrector
with a short-read layout assembler, not hifiasm/flye; the long-read tools come in
with a long-read validation dataset.

## Multi-read-set workflows (td-06er)

The single-FASTQ route above remains intentionally unchanged. Multi-read-set
workflows use a sibling adapter and correct each input independently before the
combined-input assembler runs:

```julia
Mycelia.Rhizomorph.assemble_hybrid((short_r1, short_r2), long_reads;
    config = Mycelia.Rhizomorph.UnicyclerHybridConfig())

Mycelia.Rhizomorph.assemble_hybrid((short_r1, short_r2), long_reads;
    config = Mycelia.Rhizomorph.AutocyclerPolishConfig())

Mycelia.Rhizomorph.assemble_dual_long(hifi_reads, ont_reads;
    config = Mycelia.Rhizomorph.MetaMDBGHybridConfig())
```

The paired-short inputs use their short-read error profile and the long input
uses its Nanopore or PacBio profile. Mate count, order, and normalized
identifiers are checked before and after correction so independent correction
cannot silently drop or reorder one mate. Persistent runs store each corrected
set in its own subdirectory; ephemeral runs clean every corrected FASTQ and tool
directory on both success and failure.

The contracts are deliberately non-uniform:

- **Unicycler** consumes corrected R1, corrected R2, and corrected long reads in
  one hybrid assembly call.
- **Autocycler is long-read-only.** Its upstream full-pipeline script receives
  only the corrected long FASTQ (`threads`, `jobs`, and an explicit
  `ont_r9`/`ont_r10`/`pacbio_clr`/`pacbio_hifi` type). The resulting consensus
  FASTA is then polished with the corrected paired-short reads: R1 and R2 are
  aligned separately with `bwa mem -a`, followed by Polypolish and careful
  Pypolca. The raw Autocycler GFA is retained as provenance and is never
  mislabeled as a sequence-polished graph.
- **metaMDBG** is a separate HiFi-plus-ONT dual-long workflow, not an
  Illumina-plus-long assembler. Its compressed contigs output is read through
  the normal gzip-aware FASTX path.

Every result records the workflow, technologies, input and corrected read
counts, correction settings, polishers, and caller-owned artifacts. This common
provenance surface is the assembly-side foundation for benchmark matrices and
ensemble comparison (method agreement, disagreement, and confidence strata).
Pangenome/variation consensus logic remains a downstream analysis layer; this
route supplies comparable, explicitly attributed assemblies without asserting
in advance which method is empirically best.

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
  comparable in one committed CSV. Because these are **Illumina** reads, the
  hybrid arm pairs the `:scalable` corrector with a short-read layout assembler
  (`run_megahit` / `run_metaspades`), not a long-read tool; add a long-read hybrid
  arm only once a long-read validation dataset is introduced.
- Register the arm in `rhizomorph_benchmark_manifest.toml`.
- The manuscript's pre-registered **H5** (state-of-the-art contiguity) is the
  eventual home for the cross-assembler comparison; the interim CSV is engineering
  validation only and must be labelled as such, exactly like the existing
  `real_data_corrector_validation` result.

## Missing wrappers (optional, Conda-only)

`run_verkko`, `run_raven`, `run_miniasm`, `run_wtdbg2`, `run_shasta`,
`run_nextdenovo` do not exist yet (`raven`/`miniasm`/`wtdbg2`/`shasta`/`nextdenovo`
appear only as strings in `benchmarking/03_assembly_benchmark.jl`; `verkko`
appears nowhere). None is required for a first hybrid arm — the Illumina arm uses
the existing short-read wrappers, and hifiasm/flye cover the long-read case. Each
is a small addition following the `add_bioconda_env` + `CONDA_RUNNER run -n <env>`
pattern; file separately if a benchmark tier needs them.

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
2. **Route (a) wiring:** add `layout`/`olc_tool` config + read-type-aware
   validation, and repurpose the dead `_assemble_hybrid_olc` stub (`:1585`),
   dispatching to a read-type-appropriate wrapper (short-read
   `run_megahit`/`run_metaspades` for the Illumina arm; hifiasm/flye for
   long-read) through a per-tool argument adapter; gated real-run test under
   `test/4_assembly/`.
3. **Validation arm:** add the `hybrid` arm to the phiX174/lambda harness; commit
   the comparison CSV.
4. **Manuscript update (docs, this session's domain):** once (2)–(3) land, tighten
   the Methods §"Stage 2 deconvolution" wording from "designed" to
   "implemented, benchmark pending," and add the hybrid arm to the interim
   validation table in Results.
5. **Stretch:** Route (b) GFA-substrate bridge; missing long-read wrappers as
   needed.
