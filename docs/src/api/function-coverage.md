# Function Coverage Audit

This audit maps source modules in `src/` to documentation anchors so newly added capabilities stay discoverable. Counts come from a regex sweep of the Julia sources and should be refreshed after major changes.

## Module-Level Coverage

| Source Module | Function Count | Primary Documentation | Coverage Notes |
| --- | --- | --- | --- |
| `assembly.jl` | 68 | [Assembly Suite](../api/workflows/assembly-suite.md) | Wrapper coverage for MEGAHIT/metaSPAdes/Flye/hifiasm is documented; link polishing and long-read retries into tutorials. |
| `quality-control-and-benchmarking.jl` | 30 | [Quality Control](../api/workflows/quality-control.md) | External filtering and QUAST/BUSCO/MUMmer paths are covered; native metrics (`robust_cv`, `filter_genome_outliers`) referenced in the [Workflow & Tool Map](../workflow-map.md). |
| `kmer-analysis.jl` | 39 | [Sequence Analysis](../api/workflows/sequence-analysis.md) | K-mer spectrum utilities are referenced; add examples for saved-result loaders. |
| `distance-metrics.jl` | 15 | [Sequence Analysis](../api/workflows/sequence-analysis.md) | Mash/Jaccard/JS divergence functions are linked from tutorials. |
| `pangenome-analysis.jl` | 8 | [Comparative Genomics](../api/workflows/comparative-genomics.md) | PGGB/cactus wrappers present; tutorials pending. |
| `taxonomy-and-trees.jl` | 38 | [Workflow & Tool Map](../workflow-map.md) | Classification + plotting functions surfaced in the map; needs dedicated taxonomy tutorial. |
| `variant-analysis.jl` | 15 | [Workflow & Tool Map](../workflow-map.md) | Functions exist but lack tutorial/API detail; planned variant calling tutorial will fill the gap. |
| `alignments-and-mapping.jl` | 25 | [Tool Integration](../generated/tutorials/08_tool_integration.md) | Alignment wrappers (BLAST/DIAMOND/MMSeqs2/minimap2) are referenced from the capability matrix. |
| `reference-databases.jl` | 34 | [Data Acquisition](../api/workflows/data-acquisition.md) | NCBI/SRA helpers documented; ensure coverage of metadata loaders. |
| `simulation.jl` | 41 | [Data Acquisition](../api/workflows/data-acquisition.md) | Read simulators surfaced in tutorials and workflow map. |
| `plotting-and-visualization.jl` | 27 | [Visualization Gallery](../visualization-gallery.md) | Embedding/taxa plots linked; gallery currently has some broken refs. |
| `xam.jl` | 22 | [Tool Integration](../generated/tutorials/08_tool_integration.md) | Mapping stats and coverage plotting referenced; align with BAM/CRAM primer. |
| `alignments-and-mapping.jl` | 25 | [Tool Integration](../../tutorials/08_tool_integration.md) | Alignment wrappers (BLAST/DIAMOND/MMSeqs2/minimap2) are referenced from the capability matrix. |
| `reference-databases.jl` | 34 | [Data Acquisition](../api/workflows/data-acquisition.md) | NCBI/SRA helpers documented; ensure coverage of metadata loaders. |
| `simulation.jl` | 41 | [Data Acquisition](../api/workflows/data-acquisition.md) | Read simulators surfaced in tutorials and workflow map. |
| `plotting-and-visualization.jl` | 27 | [Visualization Gallery](../visualization-gallery.md) | Embedding/taxa plots linked; gallery currently has some broken refs. |
| `xam.jl` | 22 | [Tool Integration](../../tutorials/08_tool_integration.md) | Mapping stats and coverage plotting referenced; align with BAM/CRAM primer. |
| `fastx.jl` | 50 | [Workflow & Tool Map](../workflow-map.md) | FASTA/FASTQ normalization utilities used by tutorials; add quick-start snippet to API reference. |
| `utility-functions.jl` | 102 | Internal | Shared helpers; document only when promoted to user-facing APIs. |

## How to Keep This Audit Fresh

1. When adding a new function, link it from:
   - the nearest workflow page and/or tutorial, **and**
   - the [Workflow & Tool Map](../workflow-map.md) if it creates a new input/output path.
2. For wrappers around external tools, record expected inputs/outputs in the workflow map.
3. Periodically update the counts with a small helper script, for example:

```bash
rg '^function ' src | awk -F: '{print $1}' | sort | uniq -c
```

and edit the table above.

## Cross-check with Complete API

Use the [`Complete API Surface`](all-functions.md) page as a ground truth index of documented symbols:

- If a function appears in the **code** but not on the **Complete API Surface**, it is missing a docstring.
- If it appears in the **API**, but not:

  - the [Workflow & Tool Map](../workflow-map.md), **and**
  - any tutorial,

  then it is likely undiscoverable in real workflows and should be linked.

## Immediate Documentation Tasks

- Add API snippets and at least one tutorial for `variant-analysis.jl` (normalize VCF → evaluate ROC).
- Surface taxonomy abundance plotting (`plot_taxa_abundances`) in the Comparative Genomics or Tool Integration tutorial.
- Expand the Visualization Gallery once broken references are resolved and ensure all plotting helpers appear there.
