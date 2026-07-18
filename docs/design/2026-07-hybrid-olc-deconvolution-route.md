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
GFA is never described as a sequence-polished graph. Persistent runs retain it.
Ephemeral high-level runs attempt exact identity-bound cleanup of route-owned
tool outputs after wrapping the final assembly. Cleanup is never silent: a
failure propagates when no result can be returned or retains and reports exact
evidence with the result.

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
HiFi *or* one or more explicitly attested Nanopore R10.4-or-later inputs.
`ont_reads` requires `ont_r10_4_plus = true`; generic, R9, and unknown ONT
inputs are rejected before provisioning. HiFi calls retain their existing
ergonomics and need no attestation keyword. The wrapper also fails before
provisioning when both technologies, neither technology, an empty path set, a
missing/empty input file, or a nonpositive `graph_k` is supplied. Completed
local or generated-runtime execution uses a spec-hash-addressed metaMDBG 1.4
environment and records a digest of the complete normalized resolved Conda
inventory (package name, version, build, and channel, including transitive
packages). Contract-verified reuse reports the inventory digest and package
count bound by its completion manifest. Planned and submitted results carry
expected environment-specification provenance only and do not claim that
installation or a resolved inventory has been realized. Environment discovery,
creation, inventory, local commands, and generated runtime scripts all use one
canonical Conda executable. The install lock is anchored beneath
that executable's Conda root rather than a Julia depot, so two wrappers cannot
concurrently mutate the same environment through different depots and two
independent Conda roots do not block each other. The wrapper reuses only
structurally valid, sequence-bearing contigs plus exactly one dynamic graph in
the entire output root, and that graph must have the requested k. A requested-k
graph alongside any foreign-k graph is rejected. Reuse also requires the durable
input contract to match the normalized input technology, ordered paths, file
sizes, input SHA-256 digests, `abundance_min`, and pinned environment
specification. The schema-v5 contract and its signature bind the explicit ONT
R10.4-or-later attestation. Modification time is a transient same-invocation
race snapshot only and is intentionally absent from that durable contract. A
separate atomic completion manifest binds `graph_k`, the workflow signature,
the realized package-inventory digest, and canonical artifact identities,
sizes, and SHA-256 digests. Legacy plain contigs are normalized to
`contigs.fasta.gz`. Each tool-execution lifecycle stages mode-0400 private input
copies, verifies them against the captured SHA-256 contract, and runs metaMDBG
only against those staged paths. Local execution and a submitted job's
generated runtime capture the full normalized Conda inventory before and after
all tool commands and reject any drift. Submitted jobs also recompute every
source input digest after acquiring the runtime output lock and immediately
before publication. Any
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
    # Cumulative stable source plus corrected-copy budget for this workflow.
    input_snapshot_byte_ceiling = 500_000_000_000,
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
    input_snapshot_byte_ceiling = 500_000_000_000,
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
corrected FASTQ. Exact identity-bound cleanup is attempted after success or
failure. A cleanup failure fails loudly or reports retained evidence alongside
an existing primary failure rather than deleting a replacement or silently
losing failure context; cancellation remains fail-loud. When a result can still
be returned, its provenance records retained paths. With a caller-supplied
output directory, corrected FASTQs and final tool artifacts are preserved for
auditing. A non-empty persistent output directory is rejected so
stale assemblies cannot be mistaken for current results. Exact identity-bound
removal of large alignment SAMs, filtered SAMs, and BWA index files is attempted
after successful polishing by default and after polishing failures.
`keep_intermediates = true` retains and records them explicitly and therefore
requires a persistent `output_dir`.

Both high-level configs apply `input_snapshot_byte_ceiling` cumulatively to
workflow-owned stable copies of the three original inputs and all three
correction outputs. Before each copy, the route checks both the remaining
configured budget and currently available bytes under the reserved workflow
root. A ceiling, free-space, write, or content-integrity failure attempts to
remove only the exact inode created for the partial snapshot and fails before
the next correction or combined assembler call. Cleanup failure is reported
and may retain the exact evidence rather than remove a replacement. Path-backed
copies bind the streamed source/copy digest, a source rehash, and a
consumed-snapshot rehash; all three must agree. In-memory FASTQ identity is
streamed rather than duplicating the full record set. The reserved high-level
`output_dir` (or its ephemeral equivalent) owns these snapshots, and the combined-input child
revalidates and consumes the already-bound corrected snapshots without making
a second scratch copy. Standalone `run_unicycler` and
`run_autocycler_polished` calls retain their public `input_spool_parent` and
`input_spool_byte_ceiling` controls for direct-wrapper scratch; those keys are
intentionally rejected from the high-level configs' `assembler_options`.

Missing, empty, malformed, reordered, or count-changing inputs and corrected
read sets fail before the combined assembler. Missing/empty corrected FASTQs,
zero corrected read counts, absent assemblies, empty or invalid-DNA FASTA
sequences, and malformed, duplicate-segment, or sequence-free reported GFAs
also fail loudly before result provenance is built. Persistent workflow and tool
output roots are protected by adjacent interprocess locks held from output
validation through artifact finalization and cleanup, so concurrent invocations
cannot share a lifecycle.

## Provenance and reproducibility

Every multi-input workflow `AssemblyResult` records the workflow, exact
correction technologies, input and corrected read counts, requested options,
effective per-read-set correction knobs, `max_k`, indel parameters,
substitution-error rate, polishing order, cleanup retention, and caller-owned
artifacts. Injected correction runners that do not report an effective field
record it as unavailable rather than copying the requested value.
`read_content_provenance` binds each original source to a deterministic content
identity and each corrected FASTQ to its canonical path, size, and SHA-256.
Path-backed sources record per-file canonical paths, sizes, and SHA-256 digests;
in-memory sources record their read count, identifier digest, and full
record-content digest.

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

The single-technology metaMDBG wrapper applies the same fail-closed principles
with the separate `metamdbg-1.4-3b51b282e8aa768d` environment specification.
Local execution returns
`status = :complete` only after semantic FASTA/GFA validation and durable
contract finalization. Local and generated-runtime contract and completion
publication fsync both the complete mode-0600 file and its containing directory
before reporting success. Completed execution and contract-verified reuse
report resolved package-inventory digest/count provenance. A nonlocal executor
returns `status = :planned` for a collected or dry-run job. It returns
`status = :submitted` only after a real held Slurm job's exact normalized
`(job_id, job_cluster)` reference is durably bound, its pending publication name
is durably removed, and that exact job is released. Planned and submitted
results carry expected environment-specification provenance only; neither
claims a realized package inventory. Unscoped submissions record
`job_cluster = nothing`; federation submissions retain the cluster returned by
`sbatch --parsable`. The result includes the executor result and
`expected_artifacts`; those are paths the command is expected to produce, not
completed artifacts. The submitted script independently locks the output root
and performs the same validation before atomically committing completion
provenance. Real nonlocal execution is limited to the verifiable Slurm backend.
Durable pre-submit reservations never expire automatically. A new reservation
contains only a mode-0600 `contract.json` owner record and is reported as
`submission_state = :reserved`. It is durably published through an atomic
same-parent rename; provisional owner/shared-marker crash windows remain
inspectable rather than being reported as a completed claim. The wrapper
submits the job held, then exclusively creates an adjacent empty mode-0600
pending record whose name binds the output identity, owner capability,
normalized scheduler job ID, and optional federation cluster. It fsyncs that
descriptor and then its parent directory, and only then writes and fsyncs its
canonical contents. It hardlinks that record to
`job.json`, fsyncs the owner directory, and only then releases the job. The
`job.json` directory fsync is the irreversible binding commit point; later
ordinary errors retain the exact committed binding and clean only the pending
name. Ordinary errors before that point also retain the accepted-ID pending
record; exact public binding recovery resumes only the inspected pending job
reference instead of treating it as an unsubmitted reservation. Inspection
also accepts canonical pre-upgrade schema-1 records as unscoped jobs, while all
new records use schema 2 with an explicit nullable `job_cluster`. Inspection reports
`submission_state = :submission_pending` while the accepted-ID record exists
without `job.json`, and `:submission_commit_cleanup_pending` when `job.json` has
committed but the exact pending name still needs durable removal. Only a
successfully bound owner with no pending name is `:submitted`. Operators must
resume exact pending promotion or committed-pending cleanup before manually
releasing the held Slurm job; if a hard-killed binder left lifecycle ownership,
confirmed-dead inspection must recover that ownership first. The generated
runtime rejects its exact output-, capability-, job-ID-, and cluster-bound
pending sidecar before it can publish a runtime marker, and rejects any other
same-output pending evidence before it can consume the queued owner. A runtime
job that starts after release must match the sidecar against `SLURM_JOB_ID` and,
for a federation record, `SLURM_CLUSTER_NAME` before moving through a durable,
capability-keyed state machine. It first creates and
fsyncs its shared runtime marker, atomically renames the same-parent queued
owner directory to `runtime.<capability>`, fsyncs the parent, removes the queued
shared marker, and fsyncs the parent again. At final cleanup it releases the
exact private and shared runtime locks before renaming the owner record to the
persistent,
nonblocking `consumed.<capability>` audit state and fsyncing that transition.
The owner record always retains `contract.json` and `job.json`; it is never
hidden under an ephemeral runtime temporary directory. The caller reports
`:runtime` or `:consumed` only after reconstructing the corresponding exact
on-disk state, never from absence of the original queued path.

Before attempting the private runtime lock, the generated job publishes and
fsyncs its exact runtime marker and checks the complete output domain. That
marker blocks competitors and makes a hard kill during private-lock acquisition
inspectable as scheduler-owned `runtime_claiming` state. After acquiring the
private lock, the job rechecks domain exclusivity before renaming its owner.
Each output-domain check uses a descriptor-rooted, no-follow enumerator rather
than a pathname-recursive `find`. Owner-record and output-content scans use a
strict one-level mode that captures every byte-preserving relative name, entry
type, device/inode identity, and the directory's exact name set. Output-domain
ancestor scans instead manifest only the exact PID reservation name and the
durable and pending reservation-name prefixes that are relevant at that
ancestor. Recursive descendant scans traverse every encountered directory with
an explicit iterative descriptor stack, but manifest only entries matching the
PID or durable reservation prefixes. Each traversal retains descriptor and
parent-entry identity checks and fails if a race prevents complete traversal;
unrelated file or directory churn does not invalidate an otherwise complete
reservation snapshot. The enumerator writes its private NUL-delimited inventory
only after the first complete manifest, then captures a second complete manifest
and returns the inventory only when the two reservation-relevant manifests are
identical. Persistent relevant insertions, removals, and same-name type or
identity replacements observed between those captures fail before the runtime
consumes any inventory. Directory descriptors remain bound while their entries
are inspected, so a temporary pathname replacement cannot redirect traversal
into the replacement tree.

The durable, fsynced runtime output-root reservation is the synchronization
fence for supported concurrency. A cooperating Mycelia output-root writer must
observe an already-established fence and refuse an overlapping reservation
before creating an absent parent or publishing its own marker.
metaMDBG, Unicycler, Autocycler, and the persistent multi-input workflow all run
a read-only hierarchy preflight on their canonical planned root before creating
an absent parent or PID record, then repeat the exclusivity scan after acquiring
their own PID reservation. If the runtime fence appears during that race, both
claims may briefly have symmetric reservation markers. Symmetric
post-acquisition scans ensure at least one claimant yields before overlapping
output actions; the losing claimant cleans its own reservation.
The two manifests establish a complete snapshot of reservation-relevant
pre-existing state and add fail-loud relevant drift detection around that fence;
they are not a lock-free atomic snapshot against a same-UID process that
deliberately bypasses the reservation protocol and continuously mutates the
tree. Such noncooperating filesystem tampering after the runtime fence is
unsupported and receives best-effort detection, not an atomic-exclusion
guarantee.

Enumeration is intentionally bounded and fails loudly above depth 256, above
1,000,000 entries, or above 268,435,456 bytes of manifest accounting. These
bounds prevent an adversarial or unexpectedly large output tree from exhausting
the runtime job while preserving ample space for ordinary assembly layouts. A
failed or unstable capture is retried at most three times. Each attempt uses
isolated inventory and diagnostic paths under the secure temporary directory.
Unpublished partial inventories and per-attempt diagnostics stay isolated until
the secure temporary directory's lifecycle cleanup and may be retained on
fail-closed exits. A successful helper result is accepted only when the
inventory is empty or its final byte is NUL, so an exit-zero helper that emits
a truncated record is also rejected before Bash consumption.

A caller that dies before submission can inspect the on-disk reserved
capability and reclaim it only with the exact owner token plus explicit
confirmation that submission never occurred. The complete temporary owner
record and its parent entry are durable before the shared marker, and that shared
marker is durable before the private atomic rename. A hard termination before
that rename, or before its following parent-directory fsync, therefore leaves
either the final reserved record or an exact shared-marker-paired temporary owner
record. Inspection reports the latter as `publication_state = :provisional` and
binds both private and shared filesystem identities so explicit recovery cannot
remove replacements. A hard-killed production caller also leaves its private
lifecycle lock, cleanup sentinel, and PID file. After independently proving the
caller dead, inspection with `confirm_process_dead = true` removes only those
exact same-user submitter remnants and requires the pre-existing canonical local
PID record as its evidence. Runtime owner or marker evidence, or a private lock
without a PID record, is treated as scheduler ownership and is never mutated by
that path; a live, remote, empty, malformed, noncanonical, or replaced PID file
also fails closed. The explicit flag also rejects no-op recovery when no
pre-existing PID record exists. Normal inspection, including queued inspection,
is fully read-only and lockless: phase one captures every owner, pending record,
and exact filesystem identity across the collection; only after that complete
snapshot does phase two reparse and revalidate every record plus the complete
before-and-after owner and pending-path inventories. A concurrent runtime
therefore completes without encountering an
inspector-created PID, sentinel, or private lock, while the inspector returns a
coherent stable state or fails loudly on the transition. Job-ID binding holds
the PID, sentinel, and private-lock domain. A hard kill after the pending name is
fsynced but before `job.json` commits retains the exact accepted scheduler ID;
confirmed-dead inspection validates and completes an empty pending record if
necessary, promotes it to `job.json`, or verifies the already committed exact
hardlink, then durably removes the pending name. Pending-only evidence is never
eligible for `confirm_not_submitted`. Confirmed-dead recovery removes the exact
private lock and cleanup sentinel while its replacement PID record is held,
then attempts pending promotion; if promotion fails, it releases that PID and
leaves the pending accepted-ID evidence available for exact public binding
resume rather than recreating ambiguous lifecycle locks. Explicit pre-submit
cleanup atomically
renames the owner to a same-parent `reclaiming.<capability>` record and fsyncs
the parent before it releases the queued marker. The transition holds the
canonical PID record,
cleanup sentinel, and private lock throughout. A hard kill therefore leaves an
inspectable `reclaiming` or `reclaim_release_pending` record; a second reclaim
while the first PID is live fails without changing any inode. Takeover requires
independently confirming process death and reinspecting with
`confirm_process_dead = true`. Inspection reports durable owner, queued-marker,
runtime-marker, private-lock, and cleanup-sentinel identities. An incomplete
temporary remnant is nonblocking
only when no durable same-root marker could pair with it; otherwise inspection
fails loudly instead of hiding a possibly corrupted capability. In the narrower
process-death
window after scheduler acceptance but before normal sidecar publication, an
operator must independently prove which exact scheduler job belongs to the
reservation, then call
`Mycelia.bind_metamdbg_submission_reservation_job!` with the inspected record,
exact owner token, job ID, and optional federation cluster, and
`confirm_submitted = true`. Binding and pre-submit reclaim require both exact
filesystem identities returned by inspection and refuse same-content
replacements. Runtime and consumed states can be reclaimed only from freshly
inspected identity-bound metadata with the exact recorded scheduler job
reference after explicit cancellation or terminal-failed/terminal-completed
confirmation.
Reservation directories, owner records, and job sidecars must remain
current-user-owned with modes 0700, 0600, and 0600, respectively; noncanonical,
modified, or replacement records fail closed.

The Slurm executor submits this lifecycle with `sbatch --hold --parsable`,
accepts exactly one numeric job ID with an optional federation-cluster suffix,
and durably records the exact normalized pair in the pending name and contents.
It binds the same inode as `job.json` under the output lock, durably removes the
pending name, and only then releases an unscoped job with `scontrol release
<job_id>` or a federation job with `scontrol -M <job_cluster> release <job_id>`.
Runtime execution requires the bound sidecar and an exact `SLURM_JOB_ID` match;
federation records additionally require an exact `SLURM_CLUSTER_NAME` match
before consuming the reservation. A crash after scheduler acceptance but before any pending record
leaves an unbound reservation and requires independent scheduler inspection. A
crash after pending publication is inspectable as `:submission_pending`; a crash
after `job.json` commits but before pending cleanup is
`:submission_commit_cleanup_pending`. Neither state is reported as
`:submitted`, released, or eligible for confirmed-not-submitted reclamation.
Scheduler acceptance is classified as not attempted, accepted, or unknown. Only
the not-attempted case permits automatic reservation cleanup; malformed output,
transport interruption, a thrown post-attempt runner, or any other unknown
outcome preserves the unbound reservation for inspection, explicit recovery
binding, or independently confirmed-not-submitted reclamation.

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
MYCELIA_AUTOCYCLER_TEST_JOBS=1 \
MYCELIA_ASSEMBLER_TEST_THREADS=2 \
julia --project=. -e \
  'include("test/4_assembly/third_party_assemblers_multi_input_hybrid.jl")'
```

The broad external-suite gates alone leave private-fixture smokes disabled.
`MYCELIA_RUN_MULTI_INPUT_HYBRID_SMOKE=true` opts into the Unicycler
hybrid smoke, and `MYCELIA_RUN_AUTOCYCLER_POLISHED=true` additionally opts into
the Autocycler-polished arm. A dedicated gate without the broad external gate,
or a dedicated-plus-external run with missing prerequisites, fails loudly.
`MYCELIA_ASSEMBLER_TEST_THREADS` defaults to `2` for these private-fixture smoke
paths and must parse as an integer from `1` through `4`; malformed and
out-of-range values fail before either assembler runs.
`MYCELIA_AUTOCYCLER_TEST_JOBS` controls concurrent alternative-assembler jobs
for the Autocycler-polished arm, defaults to `1`, and likewise must be an integer
from `1` through `4`. `MYCELIA_AUTOCYCLER_READ_TYPE` selects its exact long-read
chemistry and must agree with `MYCELIA_HYBRID_LONG_TECH`.

The lower-level direct-wrapper smoke in
`test/8_tool_integration/autocycler.jl` is a separate gate: it uses
`MYCELIA_RUN_EXTERNAL=true` plus `MYCELIA_RUN_AUTOCYCLER_SMOKE=true`, with
nonblank `MYCELIA_AUTOCYCLER_LONG_READS` and
`MYCELIA_AUTOCYCLER_READ_TYPE`, plus optional paired
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
