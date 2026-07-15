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

metaMDBG v1.4 accepts exactly one platform input: `--in-hifi` or `--in-ont`.
Upstream rejects a command that supplies both flags and selects one platform
parameter set internally. Therefore metaMDBG is not a HiFi-plus-ONT combined
assembler and is excluded from the multi-input adapter. `Mycelia.run_metamdbg`
remains available for a single HiFi *or* single ONT input and fails before
provisioning when both, neither, an empty path set, or a missing/empty input file
is supplied. The wrapper uses an exact metaMDBG 1.4, spec-hash-addressed
environment and validates its package inventory before execution. It reuses only
structurally valid, sequence-bearing contigs plus the exact requested-k dynamic
graph when the durable provenance contract matches the normalized input
technology, ordered paths, file size/mtime metadata, `abundance_min`, and pinned
toolchain identity; legacy plain contigs are normalized to `contigs.fasta.gz`.
Any nonempty output directory without that contract is rejected rather than
adopted as partial state.

This exclusion is intentional: modeling metaMDBG as Illumina-plus-long or as a
false dual-long workflow would make benchmark comparisons invalid.

## Correction and pairing model

Each read set is corrected independently before assembly:

1. R1 uses the configured short-read correction profile.
2. R2 uses the same short-read profile in a separate correction call.
3. Long reads use their own profile.

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
after correction. Explicit `/1` and `/2` mate roles are checked when present,
so reversed libraries fail loudly. R1, R2, and long reads must also be distinct
in-memory and physical sources. The direct `run_autocycler_polished` wrapper
performs its input mate and output-directory checks before dependency
provisioning or long-read assembly.

## Public entry points

```julia
unicycler_config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
    output_dir = "results/unicycler",
)
unicycler_result = Mycelia.Rhizomorph.assemble_hybrid(
    ("reads_R1.fastq.gz", "reads_R2.fastq.gz"),
    "reads_ont.fastq.gz";
    config = unicycler_config,
)

autocycler_config = Mycelia.Rhizomorph.AutocyclerPolishConfig(
    long_read_tech = :pacbio_hifi,
    autocycler_read_type = :pacbio_hifi,
    output_dir = "results/autocycler-polished",
)
autocycler_result = Mycelia.Rhizomorph.assemble_hybrid(
    ("reads_R1.fastq.gz", "reads_R2.fastq.gz"),
    "reads_hifi.fastq.gz";
    config = autocycler_config,
)
```

Convenience aliases `assemble_unicycler_hybrid` and
`assemble_autocycler_polished` accept R1, R2, and long reads as separate
arguments.

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
zero corrected read counts, absent assemblies, and absent Autocycler graphs
also fail loudly. Persistent workflow and tool output roots are protected by
adjacent interprocess locks held from output validation through artifact
finalization and cleanup, so concurrent invocations cannot share a lifecycle.

## Provenance and reproducibility

Every `AssemblyResult` records the workflow, exact correction technologies,
input and corrected read counts, requested options, effective per-read-set
correction knobs, `max_k`, indel parameters, substitution-error rate, polishing
order, cleanup retention, and caller-owned artifacts. Injected correction
runners that do not report an effective field record it as unavailable rather
than copying the requested value. Autocycler additionally records:

- the immutable upstream script revision and verified SHA-256;
- the bundled environment specification SHA-256; and
- the installed versions of required conda packages.

The installer uses the immutable environment name
`autocycler-0.5.2-d6aef758`, derived from the compatibility version and the
pinned bundled-environment SHA-256. A same-version interprocess lock serializes
direct installation, first-use creation, and revalidation. An existing
spec-hash-addressed environment is never recreated in place. Missing or
incompatible packages fail closed with instructions to stop active workflows
before manually removing and reinstalling the environment, while a future
spec/version receives a different name and cannot replace an environment that
an older workflow is actively using.

The single-input metaMDBG wrapper applies the same fail-closed principles with
the separate `metamdbg-1.4-c6fbbb3c2e85ffae` environment specification. Local
execution returns
`status = :complete` only after semantic FASTA/GFA validation and durable
contract finalization. A nonlocal executor returns `status = :planned` for a
collected or dry-run job and `status = :submitted` for a real submission,
together with the executor result and `expected_artifacts`; it does not present
planned paths as completed artifacts. The submitted script independently locks
the output root and performs the same validation before committing provenance.

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

Real Autocycler execution is opt-in because it is compute intensive. The fixture
must be a bacterial isolate for which mostly complete alternative long-read
assemblies are expected. Set both external-tool gates, provide all three FASTQs,
and select the exact Autocycler chemistry:

```bash
MYCELIA_RUN_EXTERNAL=true \
MYCELIA_RUN_AUTOCYCLER_POLISHED=true \
MYCELIA_HYBRID_SHORT_R1=/path/to/reads_R1.fastq.gz \
MYCELIA_HYBRID_SHORT_R2=/path/to/reads_R2.fastq.gz \
MYCELIA_HYBRID_LONG_READS=/path/to/long_reads.fastq.gz \
MYCELIA_AUTOCYCLER_READ_TYPE=ont_r10 \
julia --project=. -e \
  'include("test/4_assembly/third_party_assemblers_multi_input_hybrid.jl")'
```

Default CI leaves the dedicated Autocycler gate disabled. Once it is enabled,
missing or invalid gate prerequisites fail loudly rather than producing a green
skip.

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
