# docs/plots

Figure-generator scripts for the documentation site. Each script writes
publication-quality PNG **and** SVG into `docs/src/assets/`, which the Markdown
pages embed. Run them from the Mycelia base directory, e.g.:

```bash
julia --project=. docs/plots/rhizomorph_corrector_before_after.jl
```

## Convention: figures for a committed claim come from the committed data

**A figure that backs a committed benchmark/results claim MUST be generated from
the committed results file (CSV/parquet), never from a fresh live run.** Table,
figure, and prose then stay consistent *by construction* — they all read the
same source of truth.

Why this is a rule and not a preference: assembly (and many pipelines here) is
**non-deterministic** in details like the contig-length distribution. A figure
built from a fresh live run can therefore disagree with the committed results
table for the same inputs. This actually happened while authoring the corrector
docs — a live-run contiguity figure reported 588 naive contigs next to a table
that (from the committed CSV) reported 958. The fix was to make every figure
read the committed CSV.

Practical guidance:

- Read the committed results file (e.g. `benchmarking/results/*.csv`); do not
  re-run the analysis inside the plot script.
- If the committed file lacks a value you want to draw (e.g. per-item detail it
  never recorded), plot only the aggregates it *does* record and label anything
  derived as schematic — do not fabricate the missing detail.
- Cite the exact source file in the figure's caption/source line.
- Validate categorical palettes for colorblind-safety and keep light/dark
  legibility in mind (see the `dataviz` guidance).

## Examples

- `rhizomorph_corrector_before_after.jl` — contigs + N50, naive vs corrected,
  from `benchmarking/results/real_data_corrector_validation_*.csv`.
- `rhizomorph_corrector_graph.jl` — assembly redundancy (total length ÷ genome
  length), same CSV.
