# Tool Wrapper Status

**Verified**: 2026-04-03
**Supersedes**: the archived 2026-01-17 snapshot that still referenced deleted `planning-docs/TODO.md` files

## Scope

This document is a current wrapper-status snapshot, not an archive. Claims here
were cross-checked against:

- `src/Mycelia.jl` include list
- active wrapper sources in `src/`
- `test/runtests.jl` test gating
- wrapper tests in `test/`
- workflow docs in `docs/src/workflow-map.md` and `docs/src/metagenomic-workflow.md`

## Status Legend

- `implemented`: wrapper entry point exists in active `src/` code included by `src/Mycelia.jl`
- `default coverage`: exercised by `Pkg.test()` without `MYCELIA_RUN_EXTERNAL=true`
- `opt-in integration`: covered only when external/tool-heavy suites are enabled
- `validation only`: tests cover input validation, parsing, or executor/job generation, not a real tool run
- `legacy/in-development`: coverage exists only in legacy or `test/in_development`
- `broken`: wrapper exists, but the current smoke path is explicitly marked broken and should not be counted as supported execution coverage

## Current Snapshot

- The main gap is no longer missing wrappers in these workflow areas; it is
  uneven test depth and documentation depth.
- Stage 8 tool-integration tests are skipped unless `MYCELIA_RUN_EXTERNAL=true`.
- `GenomeFace` is intentionally disabled in code.
- The old archived claims that variant-calling and pangenome wrappers were
  missing are now false: those wrappers exist in `src/variant-analysis.jl` and
  `src/pangenome-analysis.jl`.
- The archived PenguiN "missing from metagenomic workflow docs" claim is also
  stale: `docs/src/metagenomic-workflow.md` already lists PenguiN in the
  short-read assembly workflow and tool table.
- Recently expanded comparative/pangenome wrapper coverage is now green in the
  default suite: `test/7_comparative_pangenomics/pangenome_wrappers.jl` passed
  25/25 assertions, `test/7_comparative_pangenomics/pangenome_analysis.jl`
  passed 80/80 assertions overall with 9/9 wrapper-specific validation checks,
  and `test/7_comparative_pangenomics/sequence_comparison_helpers_test.jl`
  passed 27/27 assertions on 2026-04-03.
- HyLight should be treated as broken rather than merely undercovered: the
  active suite only proves a missing-input guard (1/1 passing assertion in
  `test/4_assembly/third_party_assemblers_hybrid.jl`), while the legacy smoke
  path remains `@test_broken` in
  `test/in_development/third_party_assemblers_legacy_hylight.jl`.

## Recent Verification Counts

- `test/7_comparative_pangenomics/pangenome_wrappers.jl`: 25/25 passing on
  2026-04-03 for PGGB/Cactus validation, executor collection, and index guards.
- `test/7_comparative_pangenomics/pangenome_analysis.jl`: 80/80 passing on
  2026-04-03 overall; wrapper-specific PGGB/Cactus/vg conversion/index blocks
  contributed 9/9 passing assertions.
- `test/7_comparative_pangenomics/sequence_comparison_helpers_test.jl`: 27/27
  passing on 2026-04-03, including default validation for `skani_triangle`,
  `skani_dist`, and `run_sylph_profile`.
- `test/4_assembly/third_party_assemblers_hybrid.jl`: 8/8 passing on
  2026-04-03 with external tool runs disabled; HyLight contributes only a 1/1
  missing-input validation check in this file.
- `test/in_development/third_party_assemblers_legacy_hylight.jl`: HyLight
  legacy smoke remains explicitly `@test_broken`, so it is not counted as a
  passing supported wrapper path.

## Status Matrix

| Domain | Tools / entry points | Source | Coverage | Notes |
| --- | --- | --- | --- | --- |
| Read profiling | `run_sourmash_sketch`, `run_sourmash_search`, `run_sourmash_gather` | `src/classification.jl` | opt-in integration via `test/8_tool_integration/classification_tools.jl` | Documented in metagenomic workflow and workflow map. |
| Read profiling | `run_metaphlan`, `parse_metaphlan_profile` | `src/classification.jl` | opt-in integration via `test/8_tool_integration/classification_tools.jl` | Parser and external execution both covered when enabled. |
| Read profiling | `run_metabuli_classify`, `run_metabuli_build_db`, parsers | `src/classification.jl` | opt-in integration via `test/8_tool_integration/classification_tools.jl` | Documented in metagenomic workflow and workflow map. |
| Sketch comparison | `run_sylph_profile`, `skani_dist`, `skani_triangle` | `src/sequence-comparison.jl` | opt-in integration via `test/7_comparative_pangenomics/sequence_comparison.jl`; helper validation in `test/7_comparative_pangenomics/sequence_comparison_helpers_test.jl` | Implemented and documented; real tool runs are external-gated. |
| Mapping | `minimap_index`, `minimap_map`, `minimap_map_with_index`, `minimap_merge_map_and_split`, paired-end helpers | `src/alignments-and-mapping.jl` | opt-in integration via Stage 8 minimap tests; additional comparative coverage is external-gated | Strong wrapper surface area; no default execution coverage. |
| Coverage QC | `run_mosdepth`, `parse_mosdepth_*`, `summarize_mosdepth_qc` | `src/xam.jl` | opt-in integration via `test/5_validation/mosdepth_coverage_qc.jl` | Wrapper exists and is exercised with external validation runs. |
| Assembly validation | `run_quast` | `src/quality-control-and-benchmarking.jl` | opt-in integration via `test/5_validation/quast_busco_wrappers_test.jl` | Real smoke tests exist; external-gated. |
| Assembly validation | `run_busco`, `list_busco_datasets` | `src/quality-control-and-benchmarking.jl`, `src/busco-datasets.jl` | opt-in integration via `test/5_validation/quast_busco_wrappers_test.jl` | Real smoke tests exist; external-gated. |
| Assembly validation | `run_mummer`, `run_mummer_plot` | `src/quality-control-and-benchmarking.jl` | implemented, no direct wrapper test found in active suite | Wrapper exists; coverage gap remains. |
| Short-read assembly | `run_megahit`, `run_metaspades`, `run_skesa` | `src/assembly.jl` | opt-in integration via Stage 4 third-party assembler suites | Covered when external assembler tests are enabled. |
| Long-read assembly | `run_flye`, `run_metaflye`, `run_canu`, `run_hifiasm`, `run_metamdbg` | `src/assembly.jl` | opt-in integration via Stage 4 long-read/metagenomic suites | `run_metaflye` has an active metagenomic suite; others remain external-gated. |
| Hybrid / specialized assembly | `run_unicycler`, `run_plass_assemble`, `run_penguin_guided_nuclassemble`, `run_penguin_nuclassemble` | `src/assembly.jl` | opt-in integration for Unicycler/PLASS/PenguiN; some older coverage is legacy/in-development | Wrappers exist; not part of default suite. |
| Legacy/edge assembly | `run_hifiasm_meta` | `src/assembly.jl` | legacy/in-development via `test/4_assembly/third_party_assemblers_legacy_metagenome.jl` | Implemented but not promoted into the main active test sweep. |
| Strain-aware assembly | `run_hylight` | `src/assembly.jl` | validation only in `test/4_assembly/third_party_assemblers_hybrid.jl`; broken legacy smoke in `test/in_development/third_party_assemblers_legacy_hylight.jl` | Active coverage is limited to a missing-input guard; the legacy smoke path is explicitly `@test_broken`, so HyLight should currently be treated as broken. |
| Strain-aware assembly | `run_strong`, `run_strainy` | `src/assembly.jl` | legacy/in-development via `test/in_development/third_party_assemblers_legacy_hylight.jl` | Implemented, but active-suite support still lags behind the core assemblers. |
| Binning | `run_vamb`, `run_metabat2`, `run_metacoag`, `run_comebin`, `run_taxometer`, `run_taxvamb` | `src/binning.jl` | opt-in integration via `test/8_tool_integration/binning_tools.jl` | This file includes both validation and real tool-run blocks under external gating. |
| Post-binning | `run_drep_dereplicate`, `run_magmax_merge` | `src/binning.jl` | opt-in integration via `test/8_tool_integration/binning_tools.jl` | Implemented and externally exercised. |
| Binning gap | `run_genomeface` | `src/binning.jl` | validation only | Wrapper intentionally disabled: code raises `GenomeFace wrapper disabled`. |
| Annotation gap | `run_egapx` | `src/annotation.jl` | validation only | Intentional placeholder wrapper: `run_egapx` always throws `ArgumentError` documenting that EGAPx is not currently supported in Mycelia. |
| Pangenome | `construct_pangenome_pggb`, `construct_pangenome_cactus` | `src/pangenome-analysis.jl` | default coverage via `test/7_comparative_pangenomics/pangenome_wrappers.jl` and validation coverage in `test/7_comparative_pangenomics/pangenome_analysis.jl` | Contrary to the archived snapshot, these wrappers do exist. |
| Graph/pangenome utilities | `call_variants_from_pggb_graph`, indexing/conversion helpers | `src/pangenome-analysis.jl` | default validation via `test/7_comparative_pangenomics/pangenome_wrappers.jl` | Current coverage is validation/executor oriented. |
| Variant evaluation | `run_vcfeval`, parsing and summary helpers | `src/variant-analysis.jl` | default helper coverage in `test/5_validation/variant_analysis_parsing_test.jl`; no direct wrapper execution test found | Wrapper exists; execution coverage gap remains. |
| Variant calling | `run_gatk_haplotypecaller`, `run_freebayes`, `run_clair3`, `run_bcftools_call` | `src/variant-analysis.jl` | implemented, no direct wrapper tests found in active suite | Archived “missing wrapper” claim is obsolete; test/doc depth is the real gap. |
| Sequence encoding | `train_sentencepiece_model*`, `encode_sentencepiece`, `decode_sentencepiece`, model helpers | `src/sentencepiece.jl` | Stage 8 validation plus opt-in roundtrip coverage in `test/8_tool_integration/sentencepiece.jl` | Pure checks and external roundtrips live in an external-gated suite. |

## What Changed Relative to the Archived Snapshot

The January 17, 2026 document is no longer trustworthy as a current planning
artifact. The main corrections are:

- `planning-docs/TODO.md` is not a valid source of truth because it was removed.
- `run_mummer` exists in `src/quality-control-and-benchmarking.jl`.
- `construct_pangenome_pggb` and `construct_pangenome_cactus` exist in
  `src/pangenome-analysis.jl`.
- `run_vcfeval`, `run_gatk_haplotypecaller`, `run_freebayes`, `run_clair3`, and
  `run_bcftools_call` exist in `src/variant-analysis.jl`.
- `docs/src/metagenomic-workflow.md` already documents PenguiN support; the
  stale claim lived in this planning snapshot, not in the workflow docs.
- The dominant problem is coverage quality, not outright absence of wrappers.

## Highest-Priority Gaps

1. Add direct active-suite coverage for `run_mummer`.
2. Add direct wrapper tests for `run_vcfeval`, `run_gatk_haplotypecaller`,
   `run_freebayes`, `run_clair3`, and `run_bcftools_call`.
3. Decide whether to re-enable `run_genomeface` or remove it from workflow docs.
4. Repair or retire `run_hylight`, then decide whether `run_hifiasm_meta`,
   `run_strong`, and `run_strainy` should be promoted out of
   legacy/in-development coverage.
5. Audit stale planning-text carryovers from archived wrapper snapshots; the
   PenguiN follow-up here outlived the actual docs gap and is the pattern to
   avoid when wrapper support changes.
