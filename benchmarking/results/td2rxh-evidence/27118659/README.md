# td-2rxh evidence: :nanopore vs :illumina corrector profile at 1000bp+

## What this is

This is the pre-registered paired comparison of `assemble_genome(...; corrector=:iterative, strategy=:scalable)` run with `sequencing_tech = :nanopore` against `:illumina`.
Both arms receive byte-identical reads simulated with `Mycelia.observe(...; tech=:nanopore)`.

- Grid: 6 kb synthetic genome × read length {1000, 2000} × 20x × error {0.01, 0.02, 0.05}.
- Seeds: 1–16, so 96 pairs.
- Decision rule: the original pre-registration plus amendments 1–3, all recorded on the tracking issue before any data existed.
  `benchmarking/td2rxh_paired_analysis.jl` implements the rule.
- Tier: **PRELIMINARY**. This is a small synthetic genome with simulated errors, not the 4-tier benchmark.

## Provenance

- Data: Lawrencium array job 27118659, 16 tasks.
  Every task finished `COMPLETED` with exit 0 and 12/12 rows.
  The single run label is `td2rxh@e984af449ac3+manifest6893ba277d71`, which is Mycelia HEAD plus the `Manifest.toml` hash.
  Each task used a node-local TMPDIR.
- Harness: `benchmarking/indel_reroute_evidence.jl --cells td2rxh`, launched by `benchmarking/run_td2rxh_evidence_lrc.sbatch`.
- Analysis: `benchmarking/td2rxh_paired_analysis.jl` at Mycelia 498b2822, run over `shards/seed*.csv`.
  Its outputs are in `analysis/`: `report.md`, `td2rxh_cell_verdicts.csv`, `td2rxh_paired.csv`, `td2rxh_unpaired.csv` (empty), and `td2rxh_paired_deltas.{png,svg}`.

## Result

| cell (genome/read/cov/err) | illumina coverage | nanopore coverage | median Δ | p (Holm) | H2 | incremental s/read |
| --- | --- | --- | --- | --- | --- | --- |
| 6000/1000/20x/e0.01 | 0.930 | 0.501 | −0.426 | 0.00018 | no (worse) | 5.14 |
| 6000/1000/20x/e0.02 | 0.493 | 0.388 | −0.117 | 0.067 | no | 4.99 |
| 6000/1000/20x/e0.05 | 0.072 | 0.109 | +0.020 | 0.0013 | yes | 5.28 |
| 6000/2000/20x/e0.01 | 0.848 | 0.595 | −0.194 | 0.00061 | no (worse) | 9.21 |
| 6000/2000/20x/e0.02 | 0.624 | 0.596 | −0.0001 | 0.89 | no | 9.13 |
| 6000/2000/20x/e0.05 | 0.070 | 0.154 | +0.066 | 0.00018 | yes | 10.15 |

Coverage is the median best-contig reference coverage, a contiguity measure.
Δ is nanopore minus illumina, paired by seed.

- **H1′ (decode cost ≤ 0.30 s/read): REFUTED.** Every cell is 17–34× over the threshold.
- **Acceptance (pre-registered rule): MET.** This is at error 0.05 only: nanopore is better at both read lengths, and the median fit-identity Δ is ≥ 0.

## How to read it

- The only cells where nanopore wins are at 5% error, where both arms are severely fragmented: the best contig covers 7–15% of the genome.
  The relative gain there is real, but the absolute assemblies are poor in both arms.
- At 1% error, nanopore is significantly and substantially **worse** than illumina: about 0.43 coverage lost at 1000 bp and about 0.19 at 2000 bp.
- At 2% error the comparison is null at 2000 bp and non-significant at 1000 bp.
- In the figure, each `p_Holm` label is two-sided and carries no direction.
  The 0.01 cells have small p-values because nanopore is worse there; the point clouds sit below zero.
- The indel profile adds 5–10 s of wall time per read relative to illumina on identical reads.

The run therefore meets the pre-registered acceptance rule without showing a usable nanopore advantage.
The low-error regression is the actionable finding.
