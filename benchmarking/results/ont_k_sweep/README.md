# ONT k-selection sweep — is the pilot's ONT degeneracy a k artifact?

Date: 2026-08-05 Scripts: `benchmarking/ont_k_sweep.jl`,
`benchmarking/ont_read_identity.jl`,
`benchmarking/ont_alignment_threshold_diagnostic.jl` Data:
`ont_k_sweep_results.tsv` (one row per organism × technology × k × coverage ×
seed), `ont_k_sweep_summary.tsv` (per-stratum medians), `../ont_read_identity/`,
`../ont_alignment_threshold/`

## Question

The Track-A greedy-baseline pilot (2026-07-24) returned NGA50 = 0 in every
replicate of four of its six complete ONT cells, at 10x and 30x. Two facts made
the obvious reading — "ONT is too noisy for naked k-mer assembly" — insufficient
to act on:

1. All 132 pilot rows used **k = 31**, the Illumina-tuned primary k, including
   all 36 ONT rows. k was a hardcoded constant in that harness.
2. ONT reads came from `Mycelia.simulate_nanopore_reads`, which invokes Badread
   with **no** `--error_model`, `--qscore_model`, or `--identity`
   (`src/simulation.jl:1091`), inheriting whatever the installed binary defaults
   to.

Three hypotheses were under test:

- **H-a** — a k-selection artifact that dissolves into a per-chemistry k choice.
- **H-b** — degeneracy persists at every appropriate k; a genuine property of
  naked k-mer assembly on long reads.
- **H-c** — an assembler defect manifesting at low error rates.

## Provenance

Every number below was measured on **one host with one toolchain**. Earlier
cells computed on a macOS laptop (Badread 0.4.1, QUAST forced to 1 thread by a
known Python-3.8+ bug) were discarded rather than merged, because an
unpinned-toolchain difference inside a single table is precisely the hazard this
work documents.

| component | version                                                                           |
| --------- | --------------------------------------------------------------------------------- |
| host      | Lovelace (`foundry-lovelace`), 112-core Linux x86_64, 2 TB RAM                    |
| Badread   | **0.4.2**                                                                         |
| QUAST     | **5.3.0**                                                                         |
| Julia     | 1.10.11                                                                           |
| assembler | `Mycelia.Rhizomorph.assemble_genome`, `corrector = :none` (single-k, uncorrected) |

Badread's defaults, read from the installed binary rather than from the
wrapper's docstring:

| setting                                          | default                                       |
| ------------------------------------------------ | --------------------------------------------- |
| `--identity`                                     | `95,99,2.5` (beta: mean 95%, max 99%, sd 2.5) |
| `--error_model` / `--qscore_model`               | `nanopore2023`                                |
| `--length`                                       | `15000,13000`                                 |
| `--junk_reads` / `--random_reads` / `--chimeras` | 1% each                                       |

**The wrapper's docstring claim is accurate but unpinned.** It asserts "Oxford
Nanopore R10.4.1, nanopore2023", and for the installed binary that is correct.
But it is documentation _about_ Badread's defaults, not a parameter the wrapper
enforces — the sibling `simulate_nanopore_r941_reads` (`src/simulation.jl:1212`)
does pin
`--error_model nanopore2020 --qscore_model nanopore2020 --identity 90,98,5`. A
future Badread release that changed a default would silently move every ONT
benchmark in this repo without touching a line of Mycelia.

Badread 0.4.1 and 0.4.2 were checked against each other and are
**output-equivalent for this workload**: the per-read identity table and the
k-mer survival ladder regenerated under 0.4.2 are byte-identical to the 0.4.1
versions, and the two produce the same first read for the same seed.

## Measured read identity

Reads were mapped back to the reference with `minimap2 -ax map-ont` and identity
recomputed from CIGAR + NM, rather than inherited from a model name.

| estimate                   | mean       | median | q05    | q95    | n   |
| -------------------------- | ---------- | ------ | ------ | ------ | --- |
| BLAST identity (alignment) | **0.9440** | 0.9469 | 0.8973 | 0.9769 | 136 |
| gap-compressed (alignment) | 0.9502     | 0.9528 | 0.9104 | 0.9806 | 136 |
| Badread-declared (header)  | 0.9461     | 0.9496 | 0.9007 | 0.9796 | 139 |

3 of 139 reads (2.2%) did not align, consistent with Badread's 1% junk + 1%
random defaults. The alignment-measured and Badread-declared values agree within
~0.3 pp, so the identity distribution behaves as documented.

BLAST identity is the operative estimate because a k-mer is destroyed by every
erroneous **base** it covers, not by every error event. So **e = 0.056**:

| k   | P(error-free k-mer) | error-free k-mer coverage at 10x / 30x / 50x / 100x |
| --- | ------------------- | --------------------------------------------------- |
| 11  | 0.530               | 5.30 / 15.91 / 26.52 / 53.04                        |
| 15  | 0.421               | 4.21 / 12.64 / 21.06 / 42.12                        |
| 21  | 0.298               | 2.98 / 8.94 / 14.90 / 29.80                         |
| 31  | 0.168               | 1.67 / **5.02** / 8.37 / 16.75                      |

This is what made the pilot's result worth resolving: at k=31 and 30x there is
still ~5x error-free 31-mer coverage, so total degeneracy there is **not**
arithmetically forced.

## Validation against the pilot

All 24 cells this sweep shares with the pilot (Lambda, k=31, both technologies,
all coverages and seeds) reproduce **exactly** — contig counts, genome
fractions, and NGA50 all identical. Badread's `--seed` and this host's ART are
both reproducible against the pilot's environment. The sweep is therefore a
controlled extension of the pilot, not an independent re-measurement.

## NGA50 = 0 is a censored floor, in two distinct ways

The pilot's parser coerced QUAST's `-` to `0.0`
(`something(tryparse(Float64, "-"), 0.0)`, `track_a_baseline_benchmark.jl:227`).
This sweep keeps `NGA50` missing and records **why** in `nga50_status`. The
pilot's zeros map onto two different facts:

- `censored_no_alignment` — nothing aligned at all (pilot ONT 10x).
- `censored_gf_below_50` — contigs **did** align, but **NGA50 is undefined below
  50% genome fraction** by construction: it is the aligned-block length at which
  blocks that size or larger cover half the reference, so under half coverage no
  such length exists. Pilot ONT 30x sits here, with 5–14% of the genome
  genuinely recovered.

Reporting those as "NGA50 = 0" understates recovery in the second case. It also
makes NGA50 a **step function** of genome fraction at the 50% boundary, so
replicate seeds of one condition can land on opposite sides: ONT/k=15/30x gave
genome fraction 50.744% / 36.601% / 55.371% across seeds 42 / 123 / 456,
reporting NGA50 as 545 / absent / 592 for what is nearly the same assembly.

QUAST also reports `# misassemblies = 0` on unaligned contigs, so a 0 there on a
censored cell means non-alignment, not correctness.

## Outcome by (k, coverage) — stratified by chemistry

ONT and Illumina are different error processes and are never pooled. Medians
over 3 seeds; NGA50 shown only where defined.

### ONT

| k      | 10x                            | 30x                       | 50x                      | 100x                     |
| ------ | ------------------------------ | ------------------------- | ------------------------ | ------------------------ |
| 11     | degenerate — max contig 239 bp | degenerate — 241 bp       | degenerate — 244 bp      | degenerate — 321 bp      |
| **15** | degenerate                     | **GF 50.7%, NGA50 568.5**   | **GF 86.0%, NGA50 1509** | **GF 98.8%, NGA50 5391** |
| 21     | degenerate                     | GF 31.7%, NGA50 undefined | GF 77.5%, NGA50 1014     | GF 98.8%, NGA50 4541     |
| 31     | degenerate                     | GF 9.8%, NGA50 undefined  | GF 52.0%, NGA50 552      | GF 96.8%, NGA50 2451     |

### Illumina (control, same k ladder and coverages)

| k      | 10x                            | 30x                       | 50x                      | 100x                     |
| ------ | ------------------------------ | ------------------------- | ------------------------ | ------------------------ |
| 11     | degenerate — max contig 301 bp | degenerate — 268 bp       | degenerate — 268 bp      | degenerate — 268 bp      |
| 15     | GF 97.0%, NGA50 2135           | GF 99.7%, NGA50 11024     | GF 99.7%, NGA50 11024    | GF 99.7%, NGA50 11024    |
| **21** | GF 96.2%, NGA50 2161           | **GF 99.9%, NGA50 48125** | **GF 99.96%, NGA50 48482** | **GF 99.97%, NGA50 48488** |
| 31     | GF 94.1%, NGA50 1468           | GF 99.9%, NGA50 48464     | GF 99.96%, NGA50 48479     | GF 99.97%, NGA50 48488     |

Lambda is 48,502 bp, so Illumina at k=21/31 and ≥30x is producing an essentially
complete assembly.

### Is there a k at which ONT stops being degenerate?

**Yes at ≥30x, and it is k=15 — but "not degenerate" is not "good".**

- At **10x**, no k works. Error-free k-mer coverage is 1.7–5.3x across the
  ladder and nothing survives.
- At **30x**, k=15 moves the cell from the pilot's censored state (k=31: GF
  9.8%) to GF 50.7% with a defined NGA50 of 568. One of three seeds still falls
  below the line, so this is the boundary rather than a clean escape.
- At **50x and 100x**, k=15 is best at every coverage (GF 86.0% and 98.8%).
- **k=15 is the best ONT k at every coverage tested**, and the ordering k=15 >
  k=21 > k=31 is monotone in exactly the direction (1−e)^k predicts.

The mechanism is a two-sided squeeze, and the Illumina control is what separates
its two sides:

- **Long k fails on ONT because k-mers do not survive the error rate.** At
  k=31/30x only 16.8% of 31-mers are error-free. The same k on clean Illumina
  reads gives a complete assembly, so this is chemistry, not the assembler.
- **Short k fails on both chemistries because unitigs cannot grow.** k=11 is
  degenerate for ONT _and_ Illumina at every coverage. Across all 24 k=11 cells
  the single longest contig anywhere is 350 bp, and the per-stratum median
  longest contig stays inside 239–321 bp (ONT) and 268–301 bp (Illumina) while
  coverage rises 10x → 100x — a tenfold increase in data moves it by tens of
  bases, not past the 500 bp scoring threshold. That is a repeat-resolution
  floor, independent of error rate.

The optimum is therefore chemistry-dependent — **k=15 for ONT, k=21/31 for
Illumina** — which is what H-a predicted. But the two technologies are limited
by different things: Illumina at k=15 saturates hard at NGA50 11,024 from 30x
onward (a repeat-resolution ceiling that coverage cannot lift), while ONT at
k=15 keeps climbing 568 → 1509 → 5391 and never reaches that ceiling, because it
is still error-limited well before it becomes repeat-limited.

## Alignment-threshold diagnostic

QUAST's default minimum alignment identity is 95.0%; measured read identity is
94.4%. Censored cells were rescored at 95 / 90 / 85 / 80.5%.

Genome fraction for ONT at 30x:

| k   | seed | 95.0% | 90.0% | 85.0%     | 80.5% |
| --- | ---- | ----- | ----- | --------- | ----- |
| 31  | 42   | 9.8%  | 21.4% | **35.0%** | 35.0% |
| 31  | 123  | 5.0%  | 23.8% | **44.9%** | 45.0% |
| 31  | 456  | 14.5% | 35.9% | **53.7%** | 53.7% |
| 21  | 42   | 31.7% | 31.7% | 31.7%     | 31.7% |
| 21  | 123  | 24.3% | 24.3% | 24.3%     | 24.3% |
| 15  | 123  | 36.6% | 36.6% | 36.6%     | 36.6% |

**The threshold effect is specific to k=31.** There, relaxing the cut multiplies
genome fraction several-fold and makes NGA50 computable (521 / 610 / 697 at
85%), so those contigs _are_ approximate reconstructions scoring just under the
default. At k=15 and k=21 the same relaxation changes nothing — contigs either
align well or not at all.

Two things follow. The pilot's k=31 censoring is **partly a scoring artifact**,
but that artifact is specific to the long k it happened to use. And since the
contigs align at 85–90% while the reads measure 94.4%, contigs are **less
accurate than the reads they are built from** — consistent with chimeric joins
at graph branch points rather than simple inheritance of read error.

## Replicate variance — relevant to the pre-registration

NGA50 coefficient of variation across the 3 seeds, computed only where all three
seeds have a defined NGA50, stratified by chemistry:

| technology | evaluable strata | median CV  | max CV     |
| ---------- | ---------------- | ---------- | ---------- |
| Illumina   | 12 of 16         | **0.0047** | 0.1401     |
| ONT        | **6 of 16**      | **0.1480** | **0.3298** |

Illumina sits far inside the pre-registered CV ≈ 0.15 assumption, reproducing
the pilot's conclusion. ONT does not: its median CV is at the assumption and its
max is more than double it — and, more importantly, **10 of 16 ONT strata cannot
be evaluated on NGA50 at all** because NGA50 is undefined in at least one seed.

Genome fraction is far better behaved on ONT (CV 0.005–0.05 at 50x/100x) and is
defined wherever anything aligns. **A pre-registration using NGA50 as the ONT
endpoint would be measuring a discontinuous statistic on a stratum where it is
frequently undefined; genome fraction does not have that pathology.**

## Verdict

The data support a **compound answer**, not one hypothesis cleanly.

- **H-a — supported, and it explains the pilot's headline result.** The pilot's
  total degeneracy at 30x is substantially an artifact of k=31. Holding
  everything else fixed and moving to k=15 takes ONT/30x from genome fraction
  9.8% (NGA50 undefined) to 50.7% (NGA50 568.5). The optimal k is genuinely
  chemistry-dependent, and the pilot applied the Illumina-tuned k to ONT.
- **H-b — also supported, in the strong form.** No k rescues ONT to anything
  resembling the Illumina result. The best ONT cell measured (k=15, 100x) gives
  NGA50 5,391 against Illumina's 48,125 at k=21/30x — a ~9x gap at 3.3x the
  coverage — and at 10x no k in the ladder produces a scorable assembly at all.
  Naked k-mer assembly on these reads is poor at every k tested.
- **H-c — not supported.** The same assembler at the same k values on clean
  reads produces essentially complete assemblies (NGA50 ≈ genome length from
  30x). ONT behaviour is monotone in k in the direction (1−e)^k predicts, with
  no anomaly that requires a defect to explain. Nothing here is evidence of a
  bug.

Practically: the pilot's ONT cells should be re-run at a per-chemistry k before
they are used as a baseline, and the ONT endpoint should be genome fraction
rather than NGA50.

## What this does not determine

- **One genome.** Only Lambda (48,502 bp). Repeat structure sets the short-k
  floor and the k=15 saturation ceiling, so both are genome-specific. The pilot
  saw ONT degeneracy on T4 (168,903 bp) as well; that was not re-tested here.
- **One error regime.** Only Badread's `nanopore2023` defaults (~94.4% measured
  identity). Nothing here speaks to R9-era (~90%) or duplex/HiFi-grade (~99%)
  long reads, where the (1−e)^k arithmetic moves substantially.
- **The optimum is bracketed, not located.** k=15 is the best of {11, 15, 21,
  31}; k=11 fails, so the ONT optimum lies somewhere in (11, 21). The grid is
  too coarse to place it, and k=13 was not tested.
- **One decoder arm.** Only the `kmer` arm. This is justified by the pilot, in
  which all 66 complete qualmer/kmer pairs produced identical metrics, but arm
  equivalence was not re-verified in this sweep.
- **One assembler configuration.** `corrector = :none` throughout, deliberately,
  so k is never confounded with the iterative corrector's k-ladder. Whether the
  corrector closes the ONT gap is a separate question this sweep cannot answer.
- **Contig accuracy was not measured directly.** The claim that contigs are less
  accurate than their reads is inferred from alignment behaviour at relaxed
  identity thresholds, not from a per-contig alignment identity measurement.
- **Wall times are indicative only.** Cells ran as up to 32 concurrent shards on
  a shared host; `wall_seconds` covers assembly only (not QUAST) and is not
  benchmark-grade. No scientific column is affected.
- **One cell failed transiently and was recomputed.** ONT/k=31/100x/seed456
  first failed with `ZlibError: the compressed stream may be truncated` — a
  partially-written read file from an interrupted earlier run being reused,
  since the simulator regenerates only when the output is absent or zero-length.
  The cell was recomputed from clean inputs and the final grid is 96/96 `ok`.
