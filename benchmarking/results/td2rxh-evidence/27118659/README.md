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
- Analysis: `benchmarking/td2rxh_paired_analysis.jl`, first run at Mycelia 498b2822 over `shards/seed*.csv`.
  `report.md` and `td2rxh_cell_verdicts.csv` were regenerated after review round 7, which rounds the median deltas to 1e-6 before the median gates.
  That rounding is an implementation detail, not part of the pre-registered rule (which rounds only inside the signed-rank test).
  The inputs carry 6 decimals, so no real difference is lost and float residue can no longer decide a gate's sign: no verdict, status or p-value changed, and one displayed fit-identity median moved from 0.001333 to 0.001334.
  Tracked outputs in `analysis/`: `report.md`, `td2rxh_cell_verdicts.csv`, `td2rxh_paired.csv`, `td2rxh_unpaired.csv` (empty).
  Regenerate figures locally (not committed). From the repo root, with all
  paths relative to that root (also regenerates `report.md` and the CSVs —
  `td2rxh_paired_analysis.jl` pre-deletes them):

  ```bash
  D=benchmarking/results/td2rxh-evidence/27118659
  julia --project=. benchmarking/td2rxh_paired_analysis.jl \
    --csv "$(ls "$D"/shards/seed*.csv | paste -sd,)" \
    --output-dir "$D/analysis"
  ```

  That writes `td2rxh_paired_deltas.{png,svg}` under `$D/analysis`.

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

- **H1′ (decode cost ≤ 0.30 s/read): REFUTED.** Every cell is over the threshold, by 16.6–33.8× (4.99–10.15 s/read).
- **Acceptance (pre-registered rule): MET.** This is at error 0.05 only: nanopore is better at both read lengths, and the median fit-identity Δ is ≥ 0.

Coverage values are medians across the 16 seeds; individual seeds vary, and nanopore reaches 0.39 at 2000/e0.05 on one seed.
"Coverage lost" below is the median **paired** Δ.
At 2000/e0.01 the difference of the arm medians is larger, −0.25.
In 1000/e0.02, one of 16 nanopore runs (seed 6) did not engage the indel path.
The rule needs only ≥1 engaged run, so the cell is still COMPLETE.
`sacct_27118659.txt` is the SLURM accounting dump that backs the "16/16 COMPLETED, exit 0" statement.

## How to read it

- The only cells where nanopore wins are at 5% error, where both arms are severely fragmented: the median best contig covers 7–15% of the genome.
  The gain is in the **single best contig**.
  At that error rate nanopore produces about twice as many contigs (median about 1,490 → 2,870 at 1000 bp and about 1,550 → 3,260 at 2000 bp) and about half the total contig bases.
  In both arms the total contig bases are about 11–22× the genome length at this error rate (nanopore 10.7–13.5×, illumina 19.1–22.1×), so this is not a better assembly overall.
- The fit-identity gate is a median.
  At 2000/e0.05 the median Δ is +0.0013, but two seeds lose 0.13–0.17 fit identity.
- At 1% error, nanopore is significantly and substantially **worse** than illumina: the median paired coverage loss is 0.43 at 1000 bp and 0.19 at 2000 bp.
  That is 16/16 seeds negative at 1000 bp and 15/16 at 2000 bp.
- At 2% error the comparison is null at 2000 bp (p_Holm 0.89).
  At 1000 bp it is worse in direction (11/16 seeds negative, raw p 0.034) but not significant after Holm (0.067).
- Earlier evidence (td-4mbg: 3 seeds, 1200 bp reads, largest-contig ratio) points the same way at 1% error.
  This run is consistent with it; it is not a replication of the same measurement.
- In the figure, each `p_Holm` label is two-sided and carries no direction.
  The 0.01 cells have small p-values because nanopore is worse there; the point clouds sit below zero.
- The indel profile adds 5–10 s of wall time per read relative to illumina on identical reads.

The run therefore meets the pre-registered acceptance rule without showing a usable nanopore advantage.
The low-error regression is the actionable finding.
