# gqc v0.1 specification

Mycelia graph QC packet (`src/gqc.jl`). Semver **0.1.0**.

## Input layout

### Single-graph folder

```text
<input_dir>/
  graph.gfa
  metadata.json          # optional but recommended
  reference/             # optional; preferred for repeat nulls
    *.fa|*.fas|*.fna|*.ffn|*.fasta[+ .gz]
```

### Graph-set folder

No top-level `graph.gfa`. Each *immediate* child that contains `graph.gfa` is one graph.
Top-level `graph.gfa` always wins (single-graph mode).

## `run_gqc` keywords (highlights)

| Keyword | Default | Notes |
|---------|---------|-------|
| `null_source` | `:auto` | `:auto` / `:reference` / `:segments` |
| `bh_ks` | `[11]` | Within-graph BH family only |
| `strand_mode` | `nothing` → metadata → `:single` | `:single` / `:canonical` / `:double` |
| `canonical` | `nothing` → derived | Explicit bool is honored |
| `overwrite` | `:error` | `:skip` resume / `:replace` / `:error` |
| `n_shuffles` | `100` | Monte Carlo null replicates |

## `null_source`

- **`:auto`**: use loaded reference records when `n_records > 0`; else GFA segments.
- **`:reference`**: require loaded records; error if none.
- **`:segments`**: force segments.

**K-mer-node guard:** if every segment has the same length `L` and
(`L >= max(ks)` **or** all GFA overlaps equal `L-1` **or** `L == metadata["k"]`),
segment-sourced nulls set `null_applicable=false` (no p/q). Does **not** require
`metadata["k"]` to match a sweep k.

Fallback reasons are stored in `params.null_source_reason` and
`manifest.null_source_reason`.

Reference loader accepts `.fa`, `.fas`, `.fna`, `.ffn`, `.fasta`, optional `.gz`,
one sequence per FASTA record (no chimeric joins). Bytes before the first `>`
are ignored. Load count is logged and stored under `manifest.reference_load`.

## Length metrics (overlap-corrected)

Unitig and assembly spans subtract link overlaps along each path:

`span = sum(segment lengths) − sum(L-line overlaps on consecutive edges)`.

- Overlap from the GFA `L` CIGAR when present (leading integer, e.g. `10M` → 10).
- Else, on uniform-length (k-mer-node) graphs, default overlap is `L − 1`.

| Field | Meaning |
|-------|---------|
| `total_bp` | Sum of overlap-corrected unitig spans |
| `total_segment_bp` | Raw sum of segment lengths (no overlap correction) |
| `unitig_n50` / `unitig_l50` / `unitig_max` | From overlap-corrected unitig lengths |

A path of `k`-mers that spells an `n`-bp sequence has unitig length `n`
(not `k × n_nodes`).

## Outputs

Packet directory: `output_root/<graph_id_slug>/` (slug refuses absolute/`..` ids).

```text
manifest.json
CHECKSUMS.sha256
summary.csv / summary.json / report.md
nodes/<slug>.tsv  edges/<slug>.tsv  unitigs/<slug>.tsv
nulls/<slug>__k<k>__<kind>.tsv
plots/...
formats/
  graph.graphml
  nodes.tsv  edges.tsv
  adjacency.mtx      # Matrix Market coordinate pattern *symmetric* (lower triangle)
  node_link.json     # optional; directed GFA L lines
```

Graph-set failures append `output_root/gqc_graph_set_errors.jsonl` and continue.

## Node TSV schema (unified, BREAKING)

`nodes/<slug>.tsv` and `formats/nodes.tsv` share one frozen header
(`Mycelia.GQC_NODE_TSV_COLUMNS`):

```text
graph_id	node_id	index	length	depth	degree	in_degree	out_degree	sequence
```

| Column | Meaning |
|--------|---------|
| `graph_id` | Packet graph id (see below) |
| `node_id` | GFA segment id |
| `index` | 1-based order matching MTX row/col (`formats/nodes.tsv` / adjacency) |
| `length` | Segment length (bp) |
| `depth` | GFA `DP:f` when present; else empty |
| `degree` | Undirected simple-graph degree |
| `in_degree` / `out_degree` | Directed GFA link degrees |
| `sequence` | Segment sequence (`*` → empty) |

Both files contain the same rows for a packet (identical schema and values).

**`graph_id` filling:** `run_gqc` always passes the run’s graph id (metadata
`graph_id`, else `basename(input_dir)`). Standalone `gqc_write_formats` /
`_gqc_write_nodes_tsv` without an explicit `graph_id` fall back to
`basename(dirname(g.path))` (the folder containing `graph.gfa`), matching that
default. The column is empty only when neither an explicit id nor a usable path
parent is available.

## `summary.csv` columns

Frozen order: `Mycelia.GQC_SUMMARY_COLUMNS` (includes `null_source`,
`null_applicable`, dinucleotide/codon null fields, within-graph BH q, coverage,
concordance, `total_segment_bp`).

Missing values: empty CSV field / JSON `null`.

### Null sentinels

| Situation | Fields |
|-----------|--------|
| `null_applicable=false` | observed / mean / p / q → missing; degenerate → false |
| Degenerate shuffle | `*_null_degenerate=true`; p still conservative Monte Carlo |
| BH k not in `bh_ks` | `dinucleotide_bh_q_within_graph` missing |

### Q-value scope (important)

`dinucleotide_bh_q_within_graph` is Benjamini–Hochberg over `bh_ks` **inside one
graph only**. With the default `bh_ks=[11]` it equals that graph’s k=11 p-value.

**Cross-genus / cross-graph FDR is the caller’s job.** Use:

```julia
Mycelia.gqc_cross_bh(["pkt_a/summary.json", "pkt_b/summary.json"]; k=11)
```

which returns `(graph_id, k, p, q)` rows with BH across those graphs.

## Concordance limitation (I6)

`ref_kmer_fraction` / `ref_walkability` always RC-fold k-mers and treat
walkability as consecutive reference k-mers both present as nodes — they do
**not** read GFA `L` lines. Under `strand_mode=:single` the repeat null does not
fold RC; concordance still does. Prefer `strand_mode=:canonical` for HVP
k-mer-node graphs.

## MTX convention

`%%MatrixMarket matrix coordinate pattern symmetric` with the **lower triangle**
(`row >= col`), 1-based indices matching `formats/nodes.tsv`. Directed multi-edges
remain in `edges.tsv` / GraphML / `node_link.json`.

## Reproducibility

- RNG: `Xoshiro(shuffle_seed + 1000*salt + genus_index)`
- Shuffles generated once per null kind, then counted at each k
- `params_hash` sorts integer vectors numerically; `ks`/`bh_ks` in the manifest
  match analysis order
- `input_sha256` covers `graph.gfa`, `metadata.json`, and every matched reference file
