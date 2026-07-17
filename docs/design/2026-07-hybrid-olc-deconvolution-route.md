# Corrected multi-input hybrid assembly

**Status:** implemented for Unicycler and Autocycler-polished workflows

## Purpose

Mycelia needs comparable assembly workflows for benchmarking Rhizomorph against
external methods and for building later ensemble, pangenome, and variation
analyses. This route supplies corrected inputs and explicit provenance without
asserting that any assembler is universally best.

The existing single-input `AssemblyConfig(layout = :olc)` route remains a
one-corrected-FASTQ contract. Multi-input workflows use the separate
`_run_multi_input_assembler` adapter and do not widen `_run_olc_tool`.

## Supported contracts

| Workflow | Corrected inputs | External contract | Final sequence |
| --- | --- | --- | --- |
| Unicycler hybrid | paired short R1/R2 plus long reads | one combined `-1`, `-2`, `-l` run | Unicycler assembly FASTA |
| Autocycler-polished | long reads, then paired short R1/R2 | long-only Autocycler; separate R1/R2 polishing | careful Pypolca FASTA |

Unicycler receives corrected R1, corrected R2, and corrected long reads in the
same hybrid assembly call.

Autocycler itself remains long-read-only. Mycelia intentionally pins the
upstream 0.5.2 script revision for reproducibility; this is a compatibility pin,
not a claim that the environment tracks the current upstream release. The script
receives one corrected long FASTQ, threads per job, concurrent jobs, and an
explicit `ont_r9`, `ont_r10`, `pacbio_clr`, or `pacbio_hifi` chemistry. Its
consensus is then polished by aligning corrected R1 and R2 separately with
`bwa mem -a`, filtering the paired alignments with Polypolish, running
Polypolish, and finally running Pypolca in `--careful` mode. The raw Autocycler
GFA is never described as a sequence-polished graph. It is retained only when
the caller supplies a persistent output directory; ephemeral high-level runs
delete route-owned tool outputs after wrapping the final assembly.

Autocycler is intended for bacterial isolate genomes where most or all input
assemblies can be complete. It is not an appropriate benchmark arm when complete
alternative assemblies are infeasible, including genomes with repeats longer
than the reads. A nonempty result is not evidence that these biological
preconditions were met.

## Explicit metaMDBG exclusion

metaMDBG v1.4 accepts exactly one platform option, `--in-hifi` or `--in-ont`,
with one or more FASTQs for that technology. Upstream rejects a command that
supplies both flags and selects one platform parameter set internally. Therefore
metaMDBG is not a HiFi-plus-ONT combined assembler and is excluded from the
multi-input adapter. `Mycelia.run_metamdbg` remains available for one or more
HiFi *or* one or more ONT inputs and fails before provisioning when both,
neither, an empty path set, or a missing/empty input file is supplied. The
wrapper uses a spec-hash-addressed metaMDBG 1.4 environment and
records a digest of the complete normalized resolved Conda inventory (package
name, version, build, and channel, including transitive packages). Environment
discovery, creation, inventory, local commands, and generated runtime scripts
all use one canonical Conda executable. The install lock is anchored beneath
that executable's Conda root rather than a Julia depot, so two wrappers cannot
concurrently mutate the same environment through different depots and two
independent Conda roots do not block each other. The wrapper reuses only
structurally valid, sequence-bearing contigs plus exactly one dynamic graph in
the entire output root, and that graph must have the requested k. A requested-k
graph alongside any foreign-k graph is rejected. Reuse also requires the durable
input contract to match the normalized input technology, ordered paths, file
sizes, input SHA-256 digests, `abundance_min`, and pinned environment
specification. Modification time is a transient same-invocation race snapshot
only and is intentionally absent from the durable schema-v4 contract. A
separate atomic completion manifest binds `graph_k`, the workflow signature,
the realized package-inventory digest, and canonical artifact identities,
sizes, and SHA-256 digests. Legacy plain contigs are normalized to
`contigs.fasta.gz`. Each tool-execution lifecycle stages mode-0400 private input
copies, verifies them against the captured SHA-256 contract, and runs metaMDBG
only against those staged paths. Local and submitted lifecycles capture the
full normalized Conda inventory before and after all tool commands and reject
any drift. Submitted jobs also recompute every source input digest after
acquiring the runtime output lock and immediately before publication. Any
nonempty output directory without the input contract, or any complete artifact
set without its matching completion manifest, is rejected rather than adopted.
Partial contracted contigs are likewise never resumed because they lack
realized-stage provenance. One output root owns exactly one `graph_k` lifecycle;
a different requested graph requires a fresh output root.

This exclusion is intentional: modeling metaMDBG as Illumina-plus-long or as a
false dual-long workflow would make benchmark comparisons invalid.

## Correction and pairing model

Each read set is corrected independently before assembly:

1. R1 uses the configured short-read correction profile.
2. R2 uses the same short-read profile in a separate correction call.
3. Long reads use their own profile.

Every source must be FASTQ, whether supplied as a path or as in-memory records.
FASTA input is rejected before correction; the route never converts FASTA to
FASTQ or synthesizes Q40 qualities.

PacBio CLR and HiFi are exact correction boundaries. Autocycler's read type
selects `:pacbio_clr` or `:pacbio_hifi`; HiFi does not inherit the CLR-like
11% indel-heavy profile. The high-accuracy HiFi profile remains on the
substitution-only path, with its nominal 0.001 fallback error rate, until a
validated HiFi insertion/deletion composition is modeled. The legacy
`long_read_tech = :pacbio` alias remains CLR-like and resolves to
`:pacbio_clr`; HiFi callers must select `:pacbio_hifi` explicitly. Likewise,
new single-input hifiasm callers should use the exact `:pacbio_hifi` symbol.
Its historical `:pacbio` contract remains accepted only as a deprecated alias
that is normalized to `:pacbio_hifi` before correction.

Every corrected read set must preserve its input count and identifier order.
Pair synchronization and normalized identifiers are also checked before and
after correction. Explicit `/1` and `/2` suffixes and CASAVA `1:N:...` /
`2:N:...` descriptions are checked when present, so reversed libraries and
identifier/description role conflicts fail loudly. R1, R2, and long reads must
also be distinct in-memory and physical sources. The direct
`run_autocycler_polished` wrapper performs its input mate and output-directory
checks before dependency provisioning or long-read assembly.

## Public entry points

```julia
unicycler_config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
    output_dir = "results/unicycler",
)
unicycler_result = Mycelia.Rhizomorph.assemble_unicycler_hybrid(
    "reads_R1.fastq.gz",
    "reads_R2.fastq.gz",
    "reads_ont.fastq.gz";
    config = unicycler_config,
)

autocycler_config = Mycelia.Rhizomorph.AutocyclerPolishConfig(
    long_read_tech = :pacbio_hifi,
    autocycler_read_type = :pacbio_hifi,
    output_dir = "results/autocycler-polished",
)
autocycler_result = Mycelia.Rhizomorph.assemble_autocycler_polished(
    "reads_R1.fastq.gz",
    "reads_R2.fastq.gz",
    "reads_hifi.fastq.gz";
    config = autocycler_config,
)
```

These documented convenience aliases accept R1, R2, and long reads as separate
arguments. The typed `assemble_hybrid((r1, r2), long_reads; config)` entry point
remains available when selecting workflows dynamically.

## Ownership and failure behavior

With `output_dir = nothing`, the route owns a temporary workflow root and every
corrected FASTQ. Cleanup is best effort for ordinary filesystem failures, with
retained paths recorded in result provenance when a result can be returned;
cancellation is always rethrown. With a caller-supplied output directory,
corrected FASTQs and final tool artifacts are
preserved for auditing. A non-empty persistent output directory is rejected so
stale assemblies cannot be mistaken for current results. Large alignment SAMs,
filtered SAMs, and BWA index files are removed after successful polishing by
default and after polishing failures. `keep_intermediates = true` retains and
records them explicitly and therefore requires a persistent `output_dir`.

Missing, empty, malformed, reordered, or count-changing inputs and corrected
read sets fail before the combined assembler. Missing/empty corrected FASTQs,
zero corrected read counts, absent assemblies, empty or invalid-DNA FASTA
sequences, and malformed, duplicate-segment, or sequence-free reported GFAs
also fail loudly before result provenance is built. Persistent workflow and tool
output roots are protected by adjacent interprocess locks held from output
validation through artifact finalization and cleanup, so concurrent invocations
cannot share a lifecycle.

## Provenance and reproducibility

Every `AssemblyResult` records the workflow, exact correction technologies,
input and corrected read counts, requested options, effective per-read-set
correction knobs, `max_k`, indel parameters, substitution-error rate, polishing
order, cleanup retention, and caller-owned artifacts. Injected correction
runners that do not report an effective field record it as unavailable rather
than copying the requested value. `read_content_provenance` binds each original
source to a deterministic content identity and each corrected FASTQ to its
canonical path, size, and SHA-256. Path-backed sources record per-file canonical
paths, sizes, and SHA-256 digests; in-memory sources record their read count,
identifier digest, and full record-content digest.

The common `run_unicycler` wrapper, and therefore the Unicycler hybrid arm, uses
the realized `unicycler` Conda environment rather than claiming a spec-pinned
solve. It owns one interprocess environment lock across environment preparation,
the complete normalized package inventory (name, version, build, and channel),
and fresh local assembly. It requires both Unicycler and SPAdes and verifies
that the inventory and its SHA-256 digest are identical immediately before and
after assembly. Reused historical outputs do not claim current realized-
execution provenance. Collectors and dry-run executors return an unrealized
plan with expected FASTA/GFA paths but no durable reuse-contract path; real
nonlocal execution is rejected because it would outlive the mutable-environment
lock. The higher-level exact hybrid route therefore receives only a fresh local
exact result. Incomplete non-empty legacy output directories are preserved and
rejected rather than deleted. Its public result treats `outdir`, `assembly`, and
`graph` as canonical absolute safety fields; `requested_outdir` preserves caller
spelling for audit only and carries no artifact-access or cleanup authority.
Autocycler additionally records:

- the immutable upstream script revision and verified SHA-256;
- the bundled environment specification SHA-256; and
- the complete realized name/version/build/channel package inventory and its
  SHA-256 digest.

The installer uses the immutable environment name
`autocycler-0.5.2-d6aef758986db23c`, derived from the compatibility version and
the first 16 hexadecimal characters of the pinned bundled-environment SHA-256.
A same-version interprocess lock serializes direct installation, first-use
creation, revalidation, assembly, polishing, cleanup, and final artifact audit.
The final audit rechecks artifact bytes, FASTA/GFA semantics, and contig-ID
identity across the Autocycler, Polypolish, and Pypolca stages. An existing
spec-hash-addressed environment is never recreated in place. Missing or
incompatible packages fail closed with instructions to stop active workflows
before manually removing and reinstalling the environment, while a future
spec/version receives a different name and cannot replace an environment that
an older workflow is actively using.

`install_autocycler`, `run_autocycler`, and `run_autocycler_polished` expose the
same `conda_runner` and optional `environment_prefix` contract. Explicit
prefixes are canonicalized through their nearest existing physical ancestor,
so an install performed through one stable alias can be consumed through
another alias without selecting a different environment or lock domain.

The single-technology metaMDBG wrapper applies the same fail-closed principles with
the separate `metamdbg-1.4-3b51b282e8aa768d` environment specification. Local
execution returns
`status = :complete` only after semantic FASTA/GFA validation and durable
contract finalization. A nonlocal executor returns `status = :planned` for a
collected or dry-run job and `status = :submitted` for a real submission,
together with the executor result and `expected_artifacts`; it does not present
planned paths as completed artifacts. The submitted script independently locks
the output root and performs the same validation before atomically committing
completion provenance. Real nonlocal execution is limited to the verifiable
Slurm backend. Durable pre-submit reservations never expire automatically. A
new reservation contains only a mode-0600 `contract.json` owner record and is
reported as `submission_state = :reserved`. The wrapper submits the job held,
atomically adds a mode-0600 `job.json` sidecar binding the exact scheduler job
ID, workflow signature, and owner token, and only then releases it; inspection
then reports `submission_state = :submitted`. A runtime job that starts after
release must match the sidecar against `SLURM_JOB_ID` before atomically consuming
the reservation. If that released job starts and consumes the reservation before
the caller returns, the returned capability is reported as
`submission_state = :consumed` rather than pretending the sidecar still exists.

A caller that dies before submission can inspect the on-disk reserved
capability and reclaim it only with the exact owner token plus explicit
confirmation that submission never occurred. In the narrower process-death
window after scheduler acceptance but before normal sidecar publication, an
operator must independently prove which exact scheduler job belongs to the
reservation, then call
`Mycelia.bind_metamdbg_submission_reservation_job!` with the inspected record,
exact owner token and job ID, and
`confirm_submitted = true`. Submitted jobs can be reclaimed only with the exact
recorded scheduler job ID after explicit cancellation or terminal-failed
confirmation. Reservation directories, owner records, and job sidecars must
remain current-user-owned with modes 0700, 0600, and 0600, respectively;
noncanonical, modified, or replacement records fail closed.

The Slurm executor now submits this lifecycle with `sbatch --hold --parsable`,
accepts exactly one numeric job ID with an optional federation-cluster suffix,
durably binds that numeric ID in `job.json` under the output lock, and only then
releases that exact job with `scontrol release`. Runtime execution requires the
bound sidecar and an exact `SLURM_JOB_ID` match before consuming the reservation.
A crash before binding therefore leaves a non-runnable held job, while a crash
after binding leaves an inspectable job ID that can be released or cancelled.
Scheduler acceptance is recorded as not attempted, accepted, or unknown. Only
the not-attempted state permits automatic reservation cleanup; malformed output,
transport interruption, a thrown post-attempt runner, or any other unknown
outcome preserves the unbound reservation for inspection, explicit recovery
binding, or confirmed-not-submitted reclamation.

Updating the 0.5.2 compatibility pin is a deliberate maintenance change, not an
automatic environment refresh. An upgrade must re-pin and verify the upstream
script checksum, reconcile package constraints and command/output contracts,
and rerun the wrapper unit tests plus the real Autocycler smoke test before the
new version is used in a benchmark matrix.

## Verification

Default tests cover typed contracts, independent correction, pair preservation,
sibling-adapter argument mapping, provenance, persistent and ephemeral cleanup,
stale-output rejection, and fail-loud artifacts. The existing single-input OLC
suite remains part of the regression gate.

The real **multi-input Autocycler-polished** smoke is opt-in because it is
compute intensive. The fixture must be a bacterial isolate for which mostly
complete alternative long-read assemblies are expected. Set the broad
external-suite gate and both workflow-specific gates, provide all three FASTQs,
and select the exact Autocycler chemistry:

```bash
MYCELIA_RUN_EXTERNAL=true \
MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE=true \
MYCELIA_RUN_AUTOCYCLER_POLISHED=true \
MYCELIA_HYBRID_SHORT_R1=/path/to/reads_R1.fastq.gz \
MYCELIA_HYBRID_SHORT_R2=/path/to/reads_R2.fastq.gz \
MYCELIA_HYBRID_LONG_READS=/path/to/long_reads.fastq.gz \
MYCELIA_AUTOCYCLER_READ_TYPE=ont_r10 \
julia --project=. -e \
  'include("test/4_assembly/third_party_assemblers_multi_input_hybrid.jl")'
```

The broad external-suite gates alone leave private-fixture smokes disabled.
`MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE=true` opts into the Unicycler
hybrid smoke, and `MYCELIA_RUN_AUTOCYCLER_POLISHED=true` additionally opts into
the Autocycler-polished arm. A dedicated gate without the broad external gate,
or a dedicated-plus-external run with missing prerequisites, fails loudly.

The lower-level direct-wrapper smoke in
`test/8_tool_integration/autocycler.jl` is a separate gate: it uses
`MYCELIA_RUN_EXTERNAL=true` plus `MYCELIA_RUN_AUTOCYCLER_SMOKE=true`, with
`MYCELIA_AUTOCYCLER_LONG_READS` and optional paired
`MYCELIA_AUTOCYCLER_SHORT_READS_1` / `_2` fixtures.

## Downstream benchmark and ensemble boundary

This route is an assembly-side foundation. A downstream benchmark matrix can
run Rhizomorph, Unicycler, Autocycler-polished, and other compatible assemblers,
then score assemblies empirically. Pangenome/variation logic can compare
supported sequence or graph features across methods to identify agreement,
high-confidence consensus, and divergent low-confidence regions. That scoring
and ensemble inference is not claimed as part of this adapter.

## Upstream references

- [Autocycler automated script](https://github.com/rrwick/Autocycler/blob/c98b126eb45727584623041db1bfdbdaf7aa0923/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh)
- [Autocycler bacterial-genome scope](https://github.com/rrwick/Autocycler/blob/c98b126eb45727584623041db1bfdbdaf7aa0923/README.md)
- [Autocycler completeness requirement](https://github.com/rrwick/Autocycler/wiki#an-important-requirement)
- [Autocycler combine outputs](https://github.com/rrwick/Autocycler/wiki/Autocycler-combine)
- [Polypolish workflow](https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish)
- [Pypolca v0.3.1](https://github.com/gbouras13/pypolca/tree/v0.3.1)
- [Unicycler quick usage](https://github.com/rrwick/Unicycler#quick-usage)
- [metaMDBG v1.4 input validation](https://github.com/GaetanBenoitDev/metaMDBG/blob/62951a74ca03d044ba520d118b9de4f65eba85fd/src/Commons.hpp#L1805-L1834)
