# ONT k-selection sweep — is the pilot's ONT degeneracy a k artifact?

Date: 2026-08-05

Scripts: `benchmarking/ont_k_sweep.jl`, `benchmarking/ont_read_identity.jl`,
`benchmarking/ont_alignment_threshold_diagnostic.jl`

Data: `ont_k_sweep_results.tsv` (one row per organism x technology x k x
coverage x seed), `ont_k_sweep_summary.tsv` (per-stratum medians),
`../ont_read_identity/`, `../ont_alignment_threshold/`

Throughout, **GF** means QUAST genome fraction (percent of the reference covered
by at least one aligned contig) and **NGA50** is QUAST's aligned-block NG50.
Contigs are scored at `--min-contig 500`, matching the Track-A pilot; that
threshold is load-bearing for the k=11 result below.

## Question

The Track-A greedy-baseline pilot (2026-07-24,
`rhizomorph-paper/benchmarking/track-a-pilot/track_a_pilot_results_20260724.tsv`)
returned NGA50 = 0 in every replicate of four of its six complete ONT cells, at
10x and 30x. Two facts made the obvious reading — "ONT is too noisy for naked
k-mer assembly" — insufficient to act on:

1. All 132 pilot rows used **k = 31**, the Illumina-tuned primary k, including
   all 36 ONT rows. k was a hardcoded constant in that harness.
2. ONT reads came from `Mycelia.simulate_nanopore_reads`, which at that time
   passed Badread no `--error_model`, `--qscore_model`, or `--identity`, and so
   inherited whatever the installed binary defaulted to. (Those settings are now
   pinned explicitly, to the same values previously inherited, so the reads are
   unchanged.)

Three hypotheses were under test:

- **H-a** — a k-selection artifact that dissolves into a per-chemistry k choice.
- **H-b** — degeneracy persists at every appropriate k; a genuine property of
  naked k-mer assembly on long reads.
- **H-c** — an assembler defect manifesting at low error rates.

## Provenance

Every number below was measured on **one host with one toolchain**. Earlier
cells computed on a macOS laptop (Badread 0.4.1, QUAST forced to 1 thread by a
known Python-3.8+ bug) were discarded rather than merged.

| component | version                                                        |
| --------- | -------------------------------------------------------------- |
| host      | Lovelace (`foundry-lovelace`), 112-core Linux x86_64, 2 TB RAM |
| Badread   | **0.4.2**                                                      |
| QUAST     | **5.3.0** (`--min-contig 500`, default `--min-identity 95.0`)  |
| Julia     | 1.10.11                                                        |
| assembler | `Mycelia.Rhizomorph.assemble_genome`, `corrector = :none`      |

Badread's defaults, read from the installed binary rather than from any
docstring: `--identity 95,99,2.5`, `--error_model`/`--qscore_model`
`nanopore2023`, `--length 15000,13000`, junk/random/chimera 1% each.

Badread 0.4.1 and 0.4.2 are **output-equivalent for this workload**, and that is
checkable in this repo's history rather than asserted: `git diff` across the
commit that switched hosts shows `per_read_identity.tsv` and
`kmer_survival_ladder.tsv` completely unchanged, with only the recorded
`badread_version` field differing.

## Measured read identity

Reads were mapped back with `minimap2 -ax map-ont` and identity recomputed from
CIGAR + NM, rather than inherited from a model name. Lambda, 30x, seed 42.

| estimate                   | mean       | median | q05    | q95    | n   |
| -------------------------- | ---------- | ------ | ------ | ------ | --- |
| BLAST identity (alignment) | **0.9440** | 0.9469 | 0.8973 | 0.9769 | 136 |
| gap-compressed (alignment) | 0.9502     | 0.9528 | 0.9104 | 0.9806 | 136 |
| Badread-declared (header)  | 0.9461     | 0.9496 | 0.9007 | 0.9796 | 139 |

3 of 139 reads (2.2%) did not align. `e` is the **unweighted** mean of per-read
BLAST identities; weighting each read by its aligned length gives 0.9456, which
moves the ladder below by 1.9% (k=11) to 5.4% (k=31) and changes no conclusion.

So **e = 0.056**, and P(error-free k-mer) = (1-e)^k:

| k   | P(clean k-mer) | error-free k-mer coverage at 10x / 30x / 50x / 100x |
| --- | -------------- | --------------------------------------------------- |
| 11  | 0.530          | 5.30 / 15.91 / 26.52 / 53.04                        |
| 15  | 0.421          | 4.21 / 12.64 / 21.06 / 42.12                        |
| 21  | 0.298          | 2.98 / 8.94 / 14.90 / 29.80                         |
| 31  | 0.167          | 1.67 / **5.02** / 8.37 / 16.75                      |

At k=31 and 30x there is still ~5x error-free 31-mer coverage, so total
degeneracy there is **not** arithmetically forced. That is what made the pilot's
result worth resolving.

## Validation against the pilot

All 24 cells this sweep shares with the pilot (Lambda, k=31, both technologies,
all coverages and seeds) reproduce it: contig counts and largest contigs match
in 24/24, and genome fraction and NGA50 match wherever the pilot's parser did
not coerce a censored value to 0.0. In 6 ONT cells the pilot records `0.0` where
this sweep records a censored marker — that difference is the coercion described
next, not a disagreement in measurement.

## What a censored NGA50 means

The pilot's parser coerced QUAST's `-` to `0.0`
(`track_a_baseline_benchmark.jl:227-229`). This sweep keeps `NGA50` missing and
records **why** in `nga50_status`:

- `no_contigs_ge_min` — contigs exist but none reaches 500 bp, so QUAST declines
  to score. This is a result, not a tool failure.
- `censored_no_alignment` — scorable contigs exist, nothing aligned.
- `censored_partial_alignment` — contigs **did** align (genome fraction is a
  real, nonzero measurement) but NGA50 is still absent.

**On that last case, note what the threshold actually is.** NGA50 is an
NG-family statistic: it exists only once aligned blocks of a given length or
longer _total_ at least half the reference length. That total counts duplicated
and overlapping alignments, so it is **not** genome fraction, which measures
unique coverage. The two diverge whenever duplication ratio > 1.

An earlier version of this document claimed NGA50 was undefined below 50% genome
fraction. That is false, and this PR's own diagnostic refutes it:
Lambda/ONT/k=31/30x/seed123 rescored at 90% identity has genome fraction
**23.762%** and a defined NGA50 of **508**. Six such rows exist. The status
label was renamed accordingly, because the old name asserted a condition the
code never tests.

Two consequences follow. `misassemblies` is only interpretable where NGA50 was
measured — where nothing aligned QUAST omits the metric and this sweep records
`NA` (the pilot's `0`s there are its own coercion, not QUAST output). And
replicate seeds of one ONT condition can still land on opposite sides of the
censoring boundary: ONT/k=15/30x gave genome fraction 50.744% / 36.601% /
55.371% across seeds 42/123/456, reporting NGA50 as 545 / absent / 592.

## How "degenerate" is defined

The verdict below turns on this word, so the rule is stated here rather than
only in the source. A cell is classified by these author-chosen thresholds
(`classify_outcome` in `ont_k_sweep.jl`, which flags them as the harness's one
judgment call):

| tier            | rule                                                        |
| --------------- | ----------------------------------------------------------- |
| `degenerate`    | NGA50 censored for any reason, **or** genome fraction < 25% |
| `partial`       | recovered, but below the substantial bar                    |
| `substantial`   | genome fraction >= 90% **and** NGA50 >= 10% of the genome   |
| `near_complete` | genome fraction >= 95% **and** NGA50 >= 50% of the genome   |

**Only one of these choices is doing any work.** `classify_outcome`
short-circuits on censoring before it ever evaluates genome fraction, and
`verdict_stats.tsv` records the consequence: **0 of 240 cells** are degenerate
because of the 25% floor. Every degenerate cell in this grid is degenerate
because NGA50 was censored. So the verdict rests entirely on the choice to treat
**any** censored NGA50 as degenerate — which classifies ONT/k=15/30x/seed123
that way despite 36.6% genome fraction and a 1,303 bp largest alignment. A
reader who prefers a different rule can re-derive the whole table:
`nga50_status` and `outcome` are derived columns, recomputed from stored
measurements on every read.

## Outcome by (k, coverage) — stratified by chemistry and organism

ONT and Illumina are different error processes and are never pooled. Medians
over 3 seeds. NGA50 shown only where defined.

### Lambda (48,502 bp) — ONT

| k      | 10x                            | 30x                | 50x               | 100x               |
| ------ | ------------------------------ | ------------------ | ----------------- | ------------------ |
| 11     | degenerate (max contig 239 bp) | degenerate (241)   | degenerate (244)  | degenerate (321)   |
| 13     | degenerate                     | 35.9% / N517       | 66.2% / N779      | 81.0% / N1170      |
| **15** | degenerate                     | **50.7% / N568.5** | **86.0% / N1509** | 98.8% / N5391      |
| 17     | degenerate                     | 39.3% / N608       | 85.3% / N1307     | 99.8% / N4498      |
| 19     | degenerate                     | 39.1% / N560       | 81.9% / N1174     | **100.0% / N4385** |
| 21     | degenerate                     | 31.7% / —          | 77.5% / N1014     | 98.8% / N4541      |
| 31     | degenerate                     | 9.8% / —           | 52.0% / N552      | 96.8% / N2451      |

### Lambda — Illumina control (k=13 and k=19 omitted for width; both in the TSV)

| k   | 10x                            | 30x              | 50x              | 100x             |
| --- | ------------------------------ | ---------------- | ---------------- | ---------------- |
| 11  | degenerate (max contig 301 bp) | degenerate (268) | degenerate (268) | degenerate (268) |
| 15  | 97.0% / N2135                  | 99.7% / N11024   | 99.7% / N11024   | 99.7% / N11024   |
| 17  | 97.0% / N2259                  | 99.9% / N48464   | 99.96% / N48479  | 99.97% / N48488  |
| 21  | 96.2% / N2161                  | 99.9% / N48125   | 99.96% / N48482  | 99.97% / N48488  |
| 31  | 94.1% / N1468                  | 99.9% / N48464   | 99.96% / N48479  | 99.97% / N48488  |

### T4 (168,903 bp) — the generality check

| k      | ONT 10x    | ONT 30x      | ONT 50x           | Illumina 30x   | Illumina 50x       |
| ------ | ---------- | ------------ | ----------------- | -------------- | ------------------ |
| 11     | degenerate | degenerate   | degenerate        | degenerate     | degenerate         |
| 15     | 2.3% / —   | 48.0% / N552 | 75.0% / N1003     | 80.0% / N1304  | 78.4% / N1304      |
| **21** | 0.7% / —   | 46.6% / N544 | **89.0% / N1601** | 99.1% / N18587 | 98.9% / N18587     |
| **31** | 0.4% / —   | 21.8% / —    | 72.0% / N1017     | 99.8% / N48457 | **99.8% / N52429** |

T4 at 100x was not run: its ONT cells at that depth take hours each, and the k
ordering the question turns on is already unambiguous at 30x and 50x.

## Is there a k at which ONT stops being degenerate?

**Yes at >=30x — but which k depends on the genome, and "not degenerate" is not
"good".**

- At **10x**, no k reaches even the degenerate/partial boundary. On **Lambda**
  nothing aligns at all at QUAST's default identity (0 of 21 cells). On **T4**,
  8 of 12 cells do align, but recover only 0.4–2.6% of the genome with largest
  alignments of 508–910 bp. The mechanism is not uniform across that row:
  k>=13 fails because too few error-free k-mers survive (1.67–4.73x clean
  coverage at those k), while **k=11 fails for a different reason** — see below.
- On **Lambda**, k=15 is a genuine interior optimum at 30x and 50x: k=13 (35.9%)
  and k=17 (39.3%) are both worse at 30x. At 100x the ordering by genome
  fraction shifts to k=19 (100.0%) while NGA50 still favours k=15 — the two
  endpoints disagree.
- On **T4**, the optimum is **k=21**, not k=15: at 50x, k=21 reaches 89.0% /
  N1601 against k=15's 75.0% / N1003.

So **optimal k is not a single per-chemistry constant.** It rises with genome
size, and the Illumina control shows the same shift and explains it: Lambda is
fully resolved by k=17, while T4 needs k=31 (k=21 reaches only N18587 there).
Larger genome, more repeat structure, longer k required — on both chemistries.
Coverage pushes in the same direction, since error-free k-mer coverage is
C·(1−e)^k and raising C buys back the exponential penalty of a longer k.

The k=11 floor is separate and chemistry-independent: k=11 is degenerate for
**both** technologies and **both** organisms at every coverage. Across all 24
Lambda k=11 cells the single longest contig anywhere is 350 bp, and per-stratum
median longest contig stays inside 239–321 bp (ONT) and 268–301 bp (Illumina)
while coverage rises tenfold. That is a repeat-resolution floor meeting the 500
bp scoring threshold, not an error-rate effect.

## Alignment-threshold diagnostic

QUAST's default minimum alignment identity is 95.0%; measured read identity is
94.4%. **All 38 censored cells** were rescored at 95 / 90 / 85 / 80.5%, each
against its own reference. Genome fraction for the Lambda k=31 cells, which are
among the only ones that move:

| k   | coverage | seed | 95.0% | 90.0% | 85.0%     | 80.5% |
| --- | -------- | ---- | ----- | ----- | --------- | ----- |
| 31  | 30x      | 42   | 9.8%  | 21.4% | **35.0%** | 35.0% |
| 31  | 30x      | 123  | 5.0%  | 23.8% | **44.9%** | 45.0% |
| 31  | 30x      | 456  | 14.5% | 35.9% | **53.7%** | 53.7% |
| 31  | 10x      | 42   | —     | 13.7% | 17.9%     | 17.9% |
| 31  | 10x      | 123  | —     | 18.7% | 20.8%     | 20.8% |
| 31  | 10x      | 456  | —     | 14.8% | 22.1%     | 22.1% |
| 21  | 30x      | 42   | 31.7% | 31.7% | 31.7%     | 31.7% |
| 21  | 30x      | 123  | 24.3% | 24.3% | 24.3%     | 24.3% |
| 21  | 30x      | 456  | 39.4% | 39.4% | 39.4%     | 39.4% |
| 15  | 30x      | 123  | 36.6% | 36.6% | 36.6%     | 36.6% |

**The threshold effect is specific to k=31, and it replicates across both
genomes.** Median genome-fraction gain from relaxing the cut 95% → 85%:

| organism | k=13 | k=15 | k=17 | k=19 | k=21 | **k=31** |
| --- | --- | --- | --- | --- | --- | --- |
| Lambda | +0.0 pp | +0.0 pp | +0.0 pp | +0.0 pp | +0.0 pp | **+25.2 pp** |
| T4 | — | +0.0 pp | — | — | +0.0 pp | **+25.5 pp** |

Every k below 31 is completely insensitive to the identity threshold on both
organisms — contigs either align well or not at all. At k=31 the gain is
+25.2 pp on Lambda and +25.5 pp on T4, and NGA50 becomes computable on Lambda
(521 / 610 / 697 at 85%). Two independent genomes agreeing to within 0.3 pp is
considerably stronger evidence than the single-organism version of this finding. Note the 10x rows qualify the "nothing survives at
10x" statement above: at k=31/10x, 17.9–22.1% of the genome does align once the
identity cut is relaxed; it simply does not at QUAST's default.

Since these contigs align at 85–90% while the reads measure 94.4%, contigs are
**less accurate than the reads they are built from** — consistent with chimeric
joins at graph branch points rather than simple inheritance of read error.

## Replicate variance, and why the endpoint question is not settled

NGA50 coefficient of variation across 3 seeds, computed only where all three
seeds have a defined NGA50, stratified by chemistry (Lambda):

| technology | NGA50-evaluable strata | median CV  | max CV     |
| ---------- | ---------------------- | ---------- | ---------- |
| Illumina   | 24 of 28               | **0.0047** | 0.1401     |
| ONT        | **12 of 28**           | 0.1359     | **0.4443** |

That table is the case for _not_ using NGA50 as the ONT endpoint: at least one
seed has an undefined NGA50 in 16 of 28 ONT strata. But the obvious replacement does not survive
its own variance check. **ONT genome-fraction CV, by coverage:**

| organism | 30x               | 50x           | 100x          |
| -------- | ----------------- | ------------- | ------------- |
| Lambda   | 0.206 – **0.486** | 0.016 – 0.067 | 0.002 – 0.040 |
| T4       | 0.045 – 0.230     | 0.040 – 0.106 | not run       |

At 50x and 100x genome fraction is very well behaved. **At 30x it is not** —
Lambda/k=31 reaches CV 0.486, which is 3.2x the pre-registered 0.15 assumption
and slightly worse than NGA50's own worst stratum (0.444). The two endpoints are
therefore comparably unstable at 30x, which is precisely the stratum that
prompted this investigation. Genome fraction's advantage over NGA50 is that it
is *defined* where NGA50 is not; that advantage does not extend to being
low-variance at 30x.

Genome fraction also cannot see fragmentation. Lambda/ONT/k=15/100x has genome
fraction 97.9–98.9% — and 246,914–261,181 contigs, total scored length 2.2–5.0x
the genome, and duplication ratio 1.66–1.76. A genome-fraction-only endpoint
reads that as near-complete success.

This matters because the harness already records two censoring-immune contiguity
statistics that were not considered: **`largest_alignment` is defined on every
censored cell that aligned at all** (the k=21/30x cells report 1,094–1,576 bp
with NGA50 absent), and `NA50` is defined on several. A paired endpoint — genome
fraction plus `largest_alignment` or `NA50` — keeps the contiguity signal that
genome fraction alone discards. **This document does not establish that genome
fraction alone is the right ONT endpoint.**

## Verdict

- **H-a — supported, and it explains the pilot's headline result.** The pilot's
  degeneracy at 30x is substantially an artifact of k=31. Holding everything
  else fixed and moving to k=15 takes Lambda/ONT/30x from genome fraction 9.8%
  to 50.7%. But the corollary one would want — "use k=15 for ONT" — does **not**
  hold: T4 prefers k=21, and optimal k rises with genome size on both
  chemistries.
- **H-b — also supported, in the strong form.** No k rescues ONT to anything
  resembling the Illumina result. Across the **84** Lambda ONT cells, `outcome`
  reaches `substantial` in **4** and `near_complete` in **zero**; the 84 Lambda
  Illumina cells reach `near_complete` **36** times. On T4 the gap is starker
  still — ONT reaches `substantial` **0** times. The best single ONT cell
  anywhere gives NGA50 **8,103** (Lambda, k=19, 100x) against **48,500** for the
  best Lambda Illumina cell. At 10x no k in either organism recovers more than
  2.6% of the genome. (All counts from `verdict_stats.tsv`, which is emitted by
  the harness rather than transcribed.)
- **H-c — unrefuted, not refuted.** The same assembler at the same k values on
  clean reads produces essentially complete assemblies, and ONT degrades
  monotonically in k as (1−e)^k predicts, so nothing here _requires_ a defect.
  But the Illumina control differs from ONT in read length (~150 bp vs a 15 kb
  mean) as well as error rate, so a defect specific to long-read handling would
  be invisible to it and would produce this same monotone pattern. This design
  has low power to detect that class of bug.

Practically: the pilot's ONT cells should be re-run at a per-**organism**,
per-chemistry k before being used as a baseline, and the ONT endpoint question
should be settled with the paired-metric analysis above rather than by
substituting genome fraction for NGA50.

## What this does not determine

- **Two genomes, one of them partial.** Lambda (48,502 bp) across seven k
  values; T4 (168,903 bp) across four, with no 100x. Two points cannot
  characterise how optimal k scales with genome size — they establish only that
  it is not constant.
- **n = 3 seeds.** Every CV above is a 3-point estimate with roughly 40%
  sampling error, and the endpoint recommendation rests on those. No confidence
  intervals or tests are reported anywhere in this document; differences of a
  few tenths of a percent between adjacent k values (e.g. Lambda k=15 vs k=21 at
  100x, 98.78% vs 98.76%) are not resolvable at this n.
- **The outcome tiers are author-chosen.** The `degenerate` verdict depends on
  the 25% genome-fraction floor and on treating any censored NGA50 as
  degenerate. Both are disclosed above; neither is derived from anything.
- **One error regime.** Only Badread's `nanopore2023` defaults (~94.4% measured
  identity). Nothing here speaks to R9-era (~90%) or duplex/HiFi-grade (~99%)
  reads, where the (1−e)^k arithmetic moves substantially.
- **One decoder arm** (`kmer`), justified by the pilot's 66 identical
  qualmer/kmer pairs but not re-verified here.
- **One assembler configuration** — `corrector = :none` throughout, so k is
  never confounded with the corrector's k-ladder. Whether correction closes the
  ONT gap is a separate question, and it may make the per-k choice moot.
- **Contig accuracy was not measured directly.** The claim that contigs are less
  accurate than their reads is inferred from alignment behaviour at relaxed
  identity thresholds, not from per-contig alignment identity.
- **Both endpoints are single-aligner, reference-based.** There is no
  assembly-free check (e.g. k-mer completeness against the reads).
- **Wall times are indicative only** — up to 32 concurrent shards on a shared
  host, and `wall_seconds` covers assembly only, not read simulation or QUAST.
- **Not verifiable from this repo:** the host/QUAST/Julia versions above, and
  the discarded macOS cells. The Badread 0.4.1↔0.4.2 equivalence _is_
  verifiable, from git history.
