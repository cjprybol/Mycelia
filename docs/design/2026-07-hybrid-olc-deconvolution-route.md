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

Autocycler itself remains long-read-only. The pinned upstream script receives
one corrected long FASTQ, threads per job, concurrent jobs, and an explicit
`ont_r9`, `ont_r10`, `pacbio_clr`, or `pacbio_hifi` chemistry. Its consensus is
then polished by aligning corrected R1 and R2 separately with `bwa mem -a`,
filtering the paired alignments with Polypolish, running Polypolish, and finally
running Pypolca in `--careful` mode. The raw Autocycler GFA is never described as
a sequence-polished graph. It is retained only when the caller supplies a
persistent output directory; ephemeral high-level runs delete route-owned tool
outputs after wrapping the final assembly.

## Explicit metaMDBG exclusion

metaMDBG v1.4 accepts exactly one platform input: `--in-hifi` or `--in-ont`.
Upstream rejects a command that supplies both flags and selects one platform
parameter set internally. Therefore metaMDBG is not a HiFi-plus-ONT combined
assembler and is excluded from the multi-input adapter. `Mycelia.run_metamdbg`
remains available for a single HiFi *or* single ONT input and fails before
provisioning when both or neither are supplied.

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
substitution-only path until a validated HiFi insertion/deletion composition is
modeled.

Pair count, order, and normalized identifiers are checked before and after
correction. The direct `run_autocycler_polished` wrapper performs the same input
mate check before dependency provisioning or long-read assembly.

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
corrected FASTQ. It deletes those artifacts on success and failure. With a
caller-supplied output directory, corrected FASTQs and final tool artifacts are
preserved for auditing. A non-empty persistent output directory is rejected so
stale assemblies cannot be mistaken for current results. Large alignment SAMs,
filtered SAMs, and BWA index files are removed after successful polishing by
default; `keep_intermediates = true` retains and records them explicitly.

Missing, empty, malformed, reordered, or count-changing paired inputs fail
before the combined assembler. Missing/empty corrected FASTQs, zero corrected
read counts, absent assemblies, and absent Autocycler graphs also fail loudly.

## Provenance and reproducibility

Every `AssemblyResult` records the workflow, exact correction technologies,
input and corrected read counts, correction settings, polishing order, and
caller-owned artifacts. Autocycler additionally records:

- the immutable upstream script revision and verified SHA-256;
- the bundled environment specification SHA-256; and
- the installed versions of required conda packages.

The installer recreates a stale pre-existing `autocycler` environment before
assembly when any newly required polishing package is absent.

## Verification

Default tests cover typed contracts, independent correction, pair preservation,
sibling-adapter argument mapping, provenance, persistent and ephemeral cleanup,
stale-output rejection, and fail-loud artifacts. The existing single-input OLC
suite remains part of the regression gate.

Real Autocycler execution is opt-in because it is compute intensive. Set
`MYCELIA_RUN_EXTERNAL=true` and provide the documented FASTQ environment
variables to run the gated smoke test.

## Downstream benchmark and ensemble boundary

This route is an assembly-side foundation. A downstream benchmark matrix can
run Rhizomorph, Unicycler, Autocycler-polished, and other compatible assemblers,
then score assemblies empirically. Pangenome/variation logic can compare
supported sequence or graph features across methods to identify agreement,
high-confidence consensus, and divergent low-confidence regions. That scoring
and ensemble inference is not claimed as part of this adapter.

## Upstream references

- [Autocycler automated script](https://github.com/rrwick/Autocycler/blob/c98b126eb45727584623041db1bfdbdaf7aa0923/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh)
- [Autocycler combine outputs](https://github.com/rrwick/Autocycler/wiki/Autocycler-combine)
- [Polypolish workflow](https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish)
- [Pypolca v0.3.1](https://github.com/gbouras13/pypolca/tree/v0.3.1)
- [Unicycler quick usage](https://github.com/rrwick/Unicycler#quick-usage)
- [metaMDBG v1.4 input validation](https://github.com/GaetanBenoitDev/metaMDBG/blob/62951a74ca03d044ba520d118b9de4f65eba85fd/src/Commons.hpp#L1805-L1834)
