# Gold-Standard Genome Comparison Plan

## Goal
Build a rigorous whole-genome comparison suite that:
- reports gold-standard ANI, AAI, and POCP with coverage and orthology fractions,
- supports rearrangements and complex eukaryotic WGA,
- handles circular genomes reproducibly,
- produces standardized, cacheable outputs.

## Status legend
- [ ] not started
- [~] in progress
- [x] done
- [!] blocked

## Current baseline (audited)
- `src/sequence-comparison.jl` includes `fragment_genome`, `calculate_gold_standard_ani`,
  `calculate_gold_standard_aai`, `ani_blast` (bidirectional ANIb), `ani_mummer` (ANIm),
  `aai_rbh`, `pocp`, `compare_genomes_gold` (ANIm/ANIb/AAI/POCP/PyOrthoANI),
  `normalize_fasta`, `prepare_genome_for_comparison` (optional circular heuristic), and
  circular canonicalization (`canonical_circular_sequence`, `canonicalize_circular_fasta`).
- `src/sequence-comparison.jl` also provides `run_pyorthoani` and `fastani_list` wrappers.
- `src/alignments-and-mapping.jl` wraps MUMmer (`run_nucmer`, `run_dnadiff`,
  `run_show_coords`, `parse_dnadiff_report`, `parse_mummer_coords_table`,
  `summarize_mummer_coords`, `calculate_mummer_genome_distance`) and DIAMOND
  (`run_diamond_search`, `run_diamond_besthits`).
- External-tool tests exist in `test/7_comparative_pangenomics/gold_standard_genome_comparison.jl`
  for ANIb, ANIm, AAI, AAI-RBH, compare_genomes_gold, PyOrthoANI, and MUMmer wrappers.

## Gaps to close
- Align ANIb defaults to strict definition: 1020 bp fragments, discard tail fragments,
  min identity 30, min coverage 70, evalue 1e-15, and explicit directional reporting.
- Decide whether `calculate_gold_standard_ani`/`calculate_gold_standard_aai` remain or
  delegate to `ani_blast`/`aai_rbh` for a single canonical implementation.
- POCP and POCPu-style unique counting (no current implementation).
- Extend ANIm and ANIb outputs to include explicit parameter provenance and aligned counts.
- Eukaryotic WGA (AnchorWave, LAST) and SV calling (SyRI, Assemblytics).
- Synteny visualization (plotsr and dotplot).
- Orchestration, caching, standardized reports, and additional tests.

## Definitions and outputs (contract)
- ANIb: fragment query into 1020 bp non-overlapping pieces, discard tail; BLASTN task
  blastn, evalue 1e-15; keep best hit per fragment with identity >= 30 and coverage >= 70
  (using alignment length or qend-qstart+1). Report ANI mean and aligned-fragment fraction.
- ANIm: use MUMmer dnadiff; report 1-to-1 and M-to-M identity plus aligned fractions and
  SNP and indel counts.
- AAI: RBH from BLASTP or DIAMOND; identity >= 30, coverage >= 70, evalue <= 1e-5; report
  mean identity and orthologous fraction.
- POCP: BLASTP or DIAMOND hits; identity >= 40, coverage >= 50, evalue <= 1e-5; count
  unique conserved proteins in each direction and compute POCP; expose POCPu option.
- PyOrthoANI: run `pyorthoani` with default settings; report ANI and output path.
- All metrics: report directionality, symmetric summary, tool versions, and parameters.
- Use BioSequences.LongDNA and LongAA internally; avoid string conversions except for
  external tool IO.

## Roadmap and progress
### Phase 0 - Decisions and alignment on definitions
- [x] ANIb defaults: strict mode uses 1020 bp, tail-discard, min_id=30, min_aln_frac=0.7,
  evalue=1e-15; support `mode=:strict|:modern` with 1000 bp override.
- [x] Coverage calculation: default `coverage_mode=:alignment`, allow `:qspan`.
- [x] AAI engine: auto-select DIAMOND when available, fallback to BLASTP; allow override.
- [x] Coverage denominator: AAI default `coverage_denom=:shorter`; POCP default
  `coverage_denom=:query`, with overrides.
- [x] Circular `:auto`: keep header keyword + single-contig; add optional heuristic flag
  (off by default) and preserve explicit override.

### Phase 1 - Metrics foundation (file-level checklist)
- [x] `src/sequence-comparison.jl`: align `ani_blast` defaults to strict ANIb
  (1020 bp, tail-discard, min_id=30, min_aln_frac=0.7, evalue=1e-15); add `mode` and
  `coverage_mode`; expected: ANIb results include `ani`, `af_*`, `n_fragments_*`, and `params`.
- [x] `src/sequence-comparison.jl`: `aai_rbh` exists with RBH filtering and ortholog fraction.
- [x] `src/sequence-comparison.jl`: add POCP + POCPu (best-hit + RBH) with counts; expected:
  `pocp`, `pocpu_besthit`, `pocpu_rbh`, `c1`, `c2`, `t1`, `t2`, `params`, and `coverage_denom`.
- [x] `src/sequence-comparison.jl`: expose `compare_genomes_gold` support for POCP and
  PyOrthoANI; expected: summary fields `pocp`, `pocpu_besthit`, `pocpu_rbh`, `pyorthoani_ani`.
- [x] `src/sequence-comparison.jl`: add stable report writer (TSV/JSON) with tool versions,
  params, and directional metrics.
- [ ] `test/7_comparative_pangenomics/gold_standard_genome_comparison.jl`: add POCP tests
  and strict ANIb default regression checks; expected: gated external tests pass.

### Phase 2 - MUMmer core (ANIm and alignment parsing)
- [x] Wrap `nucmer`, `dnadiff`, and `show-coords`.
- [x] Parse `show-coords` to DataFrame and summarize alignment metrics.
- [~] Add `promer` and `delta-filter` wrappers; expected: `.delta` and filtered `.delta` paths.
- [~] Implement ANIm report extraction with aligned fractions and QC fields.

### Phase 3 - Eukaryotic WGA and SV
- [ ] Add AnchorWave wrapper (proali and genoali) and anchor prep helper (GFF + FASTA).
- [ ] Add LAST wrapper (`lastdb`, `lastal`, `last-split`) for sensitive distant homology.
- [ ] Add SyRI wrapper (consume MUMmer delta).
- [ ] Add Assemblytics wrapper (delta input) for variant size summaries.

### Phase 4 - Visualization
- [ ] Wrap plotsr for SyRI outputs.
- [ ] Provide dotplot CLI option (dotPlotly or dnadotplot); keep D-GENIES as optional
  heavy dependency.
- [ ] Add small Julia-native plot helper for quick diagnostics (CairoMakie).

### Phase 5 - Performance and scale
- [ ] Add skani wrapper and set as default fast ANI (keep ANIb and ANIm as gold standard).
- [ ] Add MMseqs2 wrappers for protein clustering at pangenome scale.

### Phase 6 - Orchestration and caching
- [~] `compare_genomes_gold` orchestrator supports ANIm/ANIb/AAI/POCP/PyOrthoANI;
  still needs versioning and cache keys.
- [x] Deterministic per-pair output directories and cache keys (inputs + params + versions).
- [x] Emit provenance and runtime metadata in reports.

### Phase 7 - Tests and documentation
- [x] Unit tests for circular canonicalization and report parsing (existing in gold standard tests).
- [x] Small fixture tests for ANIb and AAI thresholds exist; POCP tests and strict ANIb
  default regression checks added.
- [x] External tool tests gated by `MYCELIA_RUN_EXTERNAL=true` exist for ANIb/ANIm/AAI/RBH/POCP.
- [ ] Tutorial and docs updates (WGA, SV, and metrics workflow).

## Notes
- Do not re-import dependencies in leaf files; use full qualification as per repo guidelines.
- Prefer `joinpath` and stable temp dirs for external tool IO.
- Use deterministic RNGs for any synthetic test data (StableRNGs).

## Open questions
- Where to store anchor and alignment intermediates for large eukaryotic runs.
- Whether to add a native parser for .delta later (optional).
- Which dotplot tool best fits the Conda.jl footprint for CI.
