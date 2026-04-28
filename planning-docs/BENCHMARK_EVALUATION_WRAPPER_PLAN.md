# Benchmark Evaluation Wrapper Plan

## Goal

Define a reproducible assembly benchmark evaluation stack for manuscript-scale
experiments. The stack should report complementary measures of:

- assembly contiguity and reference agreement,
- conserved-gene completeness,
- k-mer based base accuracy and completeness,
- long-read supported small-scale and structural assembly errors,
- tool versions, parameters, and output provenance.

## Status

Audited: 2026-04-28

Status legend:

- implemented: entry point exists in active `src/` code
- partial: Mycelia has a related native helper but not the full external tool
  wrapper
- planned: wrapper, parser, or orchestrator contract is defined here but not
  implemented

## Current Inventory

| Tool | Mycelia status | Entry points | Active tests | Gap |
| --- | --- | --- | --- | --- |
| QUAST | implemented | `Mycelia.run_quast` | `test/5_validation/quast_busco_wrappers_test.jl` under `MYCELIA_RUN_EXTERNAL=true` | Add parsers that normalize `report.tsv` into benchmark tables and capture version/provenance. |
| BUSCO | implemented | `Mycelia.run_busco`, `Mycelia.list_busco_datasets`, `Mycelia.run_busco_download_datasets` | `test/5_validation/quast_busco_wrappers_test.jl` under `MYCELIA_RUN_EXTERNAL=true`; parser helper tested without a tool run | Add summary/full-table parsers and normalized benchmark tables. |
| Merqury | partial | `Mycelia.kmer_counts_to_merqury_qv`, `Mycelia.compute_merqury_qv`, `Mycelia.assess_assembly_quality` | native QV path in `test/5_validation/validation.jl` | Add external `run_merqury` or explicitly treat native QV as the no-CLI interim path. |
| Inspector | planned | none | none | Add wrapper, dependency guard, output parsers, and external smoke tests. |
| Alternatives | partial | `Mycelia.run_mummer`, `Mycelia.run_mummer_plot`, `Mycelia.parse_mummer_coords`, CheckM/CheckM2, CheckV, CoverM, mosdepth | uneven, mostly external-gated | Use as documented alternatives for reference structural comparison, MAG/viral completeness, and coverage QC. |

## Metric Stack by Hypothesis

| Question | Primary metric family | Tool | Required inputs | Notes |
| --- | --- | --- | --- | --- |
| Are contigs contiguous and broadly reference-consistent? | N50, L50, total length, genome fraction, misassemblies, mismatches and indels per 100 kbp | QUAST | assembly FASTA; optional reference FASTA | Reference-free mode covers contiguity only; reference mode is needed for correctness metrics. |
| Is expected conserved gene content present? | Complete, single-copy, duplicated, fragmented, and missing BUSCO fractions | BUSCO | assembly FASTA; lineage dataset or auto-lineage | Completeness is lineage-dependent and does not measure non-conserved sequence or base accuracy. |
| Are assembled bases supported by high-accuracy reads? | QV, error rate, k-mer completeness, assembly-only k-mers, spectra plots | Merqury or native QV fallback | assembly FASTA plus high-accuracy reads or meryl databases | Requires high-accuracy reads; not a structural-error detector. |
| Are long-read-supported structural errors present? | structural errors, small-scale errors, read mapping rate, Inspector QV | Inspector | assembly FASTA plus long reads; optional reference | Best for long-read assemblies with enough long-read depth. |
| Are benchmark claims robust across tool assumptions? | cross-tool metric matrix and provenance | Mycelia orchestrator | all configured tool inputs | Do not collapse metrics into a single score; retain per-tool applicability flags. |

## Common Wrapper Contract

All follow-on wrappers and parsers should use this shared contract:

- Validate required input files before invoking Bioconda or external commands
  when possible, so default tests can exercise missing-input errors without
  requiring tools.
- Use `add_bioconda_env("<tool>")` for supported Bioconda tools and run through
  `Mycelia.CONDA_RUNNER`.
- Return output paths from `run_*` wrappers and provide separate parser helpers
  that return `DataFrames.DataFrame` rows.
- Emit clear errors for missing files, invalid argument combinations, missing
  expected output files, and failed command execution.
- External execution tests must be gated by `MYCELIA_RUN_EXTERNAL=true`; default
  tests should cover argument validation and parser fixtures.
- Record `tool`, `tool_version`, `assembly_id`, `input_paths`, `outdir`,
  `params_json`, `started_at`, and `completed_at` in normalized summaries.
- Keep tool citations in docstrings or docs when adding wrapper code.

## QUAST Contract

Status: implemented wrapper, planned parser.

Citation: Gurevich et al. QUAST: quality assessment tool for genome assemblies.
Bioinformatics 29(8), 1072-1075 (2013). DOI:
`10.1093/bioinformatics/btt086`.

Applicability boundary:

- Use for assembly contiguity summaries with or without a reference.
- Use reference mode for genome fraction, duplication ratio, misassemblies,
  NGA/NA metrics, mismatch rates, and indel rates.
- Use `metaquast.py` as a separate follow-on for metagenomes; do not overload
  `run_quast` with metagenome-specific defaults.
- Avoid interpreting reference-based mismatches as assembly errors when the
  sample is expected to differ from the reference.

Current interface:

```julia
Mycelia.run_quast(assembly_files::Vector{String};
    outdir::Union{String, Nothing}=nothing,
    reference::Union{String, Nothing}=nothing,
    threads::Int=Mycelia.get_default_threads(),
    min_contig::Int=500,
    gene_finding::Bool=false)
```

Planned parser interface:

```julia
Mycelia.parse_quast_report(report_tsv::String; assembly_id::Union{String, Nothing}=nothing)
Mycelia.summarize_quast(outdir::String; reference::Union{String, Nothing}=nothing)
```

Expected summary schema: `quast_summary`.

| Column | Type | Description |
| --- | --- | --- |
| `assembly_id` | `String` | Stable assembly label from the report column or caller override. |
| `tool` | `String` | Always `QUAST`. |
| `tool_version` | `Union{String, Missing}` | Parsed from QUAST logs or command probe. |
| `reference_id` | `Union{String, Missing}` | Reference FASTA basename or explicit identifier. |
| `min_contig` | `Int` | Minimum contig length used. |
| `n_contigs` | `Union{Int, Missing}` | Number of contigs at the configured threshold. |
| `largest_contig_bp` | `Union{Int, Missing}` | Largest contig length. |
| `total_length_bp` | `Union{Int, Missing}` | Total assembly length. |
| `n50_bp` | `Union{Int, Missing}` | N50 contig length. |
| `l50` | `Union{Int, Missing}` | L50 contig count. |
| `gc_percent` | `Union{Float64, Missing}` | GC percentage. |
| `genome_fraction_percent` | `Union{Float64, Missing}` | Reference-covered fraction, reference mode only. |
| `duplication_ratio` | `Union{Float64, Missing}` | Aligned assembly bases divided by aligned reference bases. |
| `misassemblies` | `Union{Int, Missing}` | QUAST global misassembly count. |
| `local_misassemblies` | `Union{Int, Missing}` | QUAST local misassembly count. |
| `mismatches_per_100kbp` | `Union{Float64, Missing}` | Reference mismatch rate. |
| `indels_per_100kbp` | `Union{Float64, Missing}` | Reference indel rate. |
| `unaligned_length_bp` | `Union{Int, Missing}` | Total unaligned assembly length. |
| `report_tsv` | `String` | Source `report.tsv`. |
| `report_html` | `Union{String, Missing}` | Source `report.html`, if present. |
| `params_json` | `String` | JSON-encoded parameters. |

Missing dependency behavior:

- If the Bioconda environment cannot be created or `quast.py` cannot be run,
  throw an error that includes `QUAST unavailable` and the attempted environment
  name `quast`.
- Tests without `MYCELIA_RUN_EXTERNAL=true` should only assert argument
  validation and parser behavior from committed fixture reports.

## BUSCO Contract

Status: implemented wrapper, planned parsers.

Citation: Simao et al. BUSCO: assessing genome assembly and annotation
completeness with single-copy orthologs. Bioinformatics 31(19), 3210-3212
(2015). DOI: `10.1093/bioinformatics/btv351`.

Applicability boundary:

- Use for conserved single-copy ortholog completeness in genomes,
  transcriptomes, proteins, and bins when the lineage is biologically
  appropriate.
- Report lineage and auto-lineage mode with every result; BUSCO scores are not
  comparable when lineage datasets differ.
- Treat duplicated BUSCOs as possible biological duplication, haplotig
  retention, bin contamination, or assembly duplication rather than as a single
  universal failure mode.
- Do not use BUSCO as a base-level accuracy or structural correctness metric.

Current interface:

```julia
Mycelia.run_busco(assembly_files::Vector{String};
    outdir::Union{String, Nothing}=nothing,
    lineage::Union{String, Nothing}=nothing,
    mode::String="genome",
    threads::Int=Mycelia.get_default_threads(),
    force::Bool=false,
    auto_lineage::Bool=true,
    auto_lineage_euk::Bool=false,
    auto_lineage_prok::Bool=false)
```

Planned parser interface:

```julia
Mycelia.parse_busco_short_summary(summary_file::String)
Mycelia.parse_busco_full_table(full_table::String)
Mycelia.summarize_busco(outdir::String)
```

Expected summary schema: `busco_summary`.

| Column | Type | Description |
| --- | --- | --- |
| `assembly_id` | `String` | Assembly basename or caller-provided label. |
| `tool` | `String` | Always `BUSCO`. |
| `tool_version` | `Union{String, Missing}` | BUSCO version. |
| `mode` | `String` | `genome`, `transcriptome`, or `proteins`. |
| `lineage_dataset` | `String` | Specific lineage used after auto-lineage resolution. |
| `auto_lineage_mode` | `Union{String, Missing}` | `auto`, `auto_euk`, `auto_prok`, or missing. |
| `complete_percent` | `Float64` | Complete BUSCO percentage. |
| `single_copy_percent` | `Float64` | Single-copy complete percentage. |
| `duplicated_percent` | `Float64` | Duplicated complete percentage. |
| `fragmented_percent` | `Float64` | Fragmented percentage. |
| `missing_percent` | `Float64` | Missing percentage. |
| `n_markers` | `Int` | Total BUSCO groups searched. |
| `n_complete` | `Int` | Complete BUSCO count. |
| `n_single_copy` | `Int` | Single-copy BUSCO count. |
| `n_duplicated` | `Int` | Duplicated BUSCO count. |
| `n_fragmented` | `Int` | Fragmented BUSCO count. |
| `n_missing` | `Int` | Missing BUSCO count. |
| `short_summary` | `String` | Source summary path. |
| `full_table` | `Union{String, Missing}` | Full BUSCO table path, if present. |
| `params_json` | `String` | JSON-encoded parameters. |

Missing dependency behavior:

- If BUSCO is unavailable, throw an error that includes `BUSCO unavailable` and
  the attempted environment name `busco`.
- If lineage data must be downloaded and network access fails, report the
  lineage and BUSCO download directory in the error.
- Default tests should use short-summary and full-table fixtures; external
  smoke tests may use a small bacterial genome with a pinned lineage.

## Merqury Contract

Status: partial native QV support, planned external CLI wrapper.

Citation: Rhie et al. Merqury: reference-free quality, completeness, and phasing
assessment for genome assemblies. Genome Biology 21, 245 (2020). DOI:
`10.1186/s13059-020-02134-9`.

Applicability boundary:

- Use for reference-free base accuracy and completeness when high-accuracy reads
  or trusted meryl databases are available.
- Prefer the external Merqury path for manuscript figures because it includes
  spectra plots, completeness, and optional haplotype-specific reports.
- Use Mycelia native QV helpers as a lightweight fallback for internal benchmark
  sweeps, and label those rows as `method = "mycelia_native_qv"`.
- Do not use Merqury alone to claim structural correctness; pair with QUAST,
  MUMmer, or Inspector depending on available evidence.

Planned wrapper interface:

```julia
Mycelia.run_merqury(assembly_file::String, reads::Vector{String};
    outdir::Union{String, Nothing}=nothing,
    k::Int=21,
    threads::Int=Mycelia.get_default_threads(),
    meryl_db::Union{String, Nothing}=nothing,
    haplotype_meryl_dbs::Vector{String}=String[],
    label::Union{String, Nothing}=nothing)
```

Native fallback interface:

```julia
Mycelia.summarize_native_merqury_qv(assembly_file::String, reads_file::String; k::Int=21)
```

Expected summary schema: `merqury_summary`.

| Column | Type | Description |
| --- | --- | --- |
| `assembly_id` | `String` | Assembly label. |
| `tool` | `String` | `Merqury` or `MyceliaNativeQV`. |
| `tool_version` | `Union{String, Missing}` | Merqury and meryl version, or Mycelia version. |
| `method` | `String` | `merqury_cli` or `mycelia_native_qv`. |
| `k` | `Int` | K-mer size. |
| `reads_id` | `String` | Read set label or meryl database label. |
| `qv` | `Union{Float64, Missing}` | Assembly QV. |
| `error_rate` | `Union{Float64, Missing}` | Error rate implied by QV, when available. |
| `completeness_percent` | `Union{Float64, Missing}` | Merqury completeness estimate. |
| `shared_kmers` | `Union{Int, Missing}` | Assembly k-mers found in reads. |
| `assembly_only_kmers` | `Union{Int, Missing}` | Assembly k-mers absent from reads. |
| `total_assembly_kmers` | `Union{Int, Missing}` | Total unique assembly k-mers. |
| `spectra_cn_plot` | `Union{String, Missing}` | Copy-number spectra plot path. |
| `spectra_asm_plot` | `Union{String, Missing}` | Assembly spectra plot path. |
| `haplotype_block_n50` | `Union{Float64, Missing}` | Optional phasing metric. |
| `outdir` | `String` | Merqury output directory. |
| `params_json` | `String` | JSON-encoded parameters. |

Missing dependency behavior:

- If the external path is requested and `meryl` or `merqury.sh` is unavailable,
  throw an error that includes `Merqury unavailable` and names the missing
  executable.
- Do not silently fall back to native QV unless the caller requested
  `method = :auto`; if fallback happens, record it in `method`.
- Default tests should cover missing input files, invalid `k`, and parsers for
  representative Merqury output snippets.

## Inspector Contract

Status: planned.

Citation: Chen et al. Accurate long-read de novo assembly evaluation with
Inspector. Genome Biology 22, 312 (2021). DOI:
`10.1186/s13059-021-02527-4`.

Applicability boundary:

- Use for long-read assembly evaluation when long reads from the assembly
  sample are available at sufficient depth.
- Use reference-free mode to identify read-supported small-scale and structural
  assembly errors without requiring a closely matched reference.
- Use reference-guided outputs only when the reference is appropriate for the
  biological sample.
- Do not run on short-read-only assemblies without long reads; use Merqury,
  QUAST, and BUSCO instead.

Planned wrapper interface:

```julia
Mycelia.run_inspector(assembly_file::String, reads::Vector{String};
    outdir::Union{String, Nothing}=nothing,
    datatype::String="pacbio-hifi",
    genome_size::Union{String, Int, Nothing}=nothing,
    reference::Union{String, Nothing}=nothing,
    threads::Int=Mycelia.get_default_threads(),
    run_correction::Bool=false)
```

Planned parser interface:

```julia
Mycelia.parse_inspector_summary(summary_statistics::String)
Mycelia.parse_inspector_errors(error_file::String)
Mycelia.summarize_inspector(outdir::String)
```

Expected summary schema: `inspector_summary`.

| Column | Type | Description |
| --- | --- | --- |
| `assembly_id` | `String` | Assembly label. |
| `tool` | `String` | Always `Inspector`. |
| `tool_version` | `Union{String, Missing}` | Inspector version. |
| `datatype` | `String` | Long-read datatype passed to Inspector. |
| `reads_id` | `String` | Long-read set label. |
| `reference_id` | `Union{String, Missing}` | Reference label when reference-guided mode is used. |
| `total_length_bp` | `Union{Int, Missing}` | Assembly total length. |
| `contig_n50_bp` | `Union{Int, Missing}` | Contig N50 from Inspector summary. |
| `read_mapping_rate` | `Union{Float64, Missing}` | Fraction or percentage of reads mapped to contigs. |
| `effective_depth` | `Union{Float64, Missing}` | Estimated usable read depth. |
| `qv` | `Union{Float64, Missing}` | Inspector QV. |
| `small_scale_errors` | `Union{Int, Missing}` | Total small-scale errors. |
| `structural_errors` | `Union{Int, Missing}` | Total structural errors. |
| `substitution_errors` | `Union{Int, Missing}` | Substitution-like errors. |
| `insertion_errors` | `Union{Int, Missing}` | Small expansion or insertion errors. |
| `deletion_errors` | `Union{Int, Missing}` | Small collapse or deletion errors. |
| `collapse_errors` | `Union{Int, Missing}` | Structural collapse errors. |
| `expansion_errors` | `Union{Int, Missing}` | Structural expansion errors. |
| `inversion_errors` | `Union{Int, Missing}` | Structural inversion errors. |
| `haplotype_switch_errors` | `Union{Int, Missing}` | Haplotype switch errors. |
| `summary_statistics` | `String` | Inspector `summary_statistics` path. |
| `structural_error_file` | `Union{String, Missing}` | Structural error report path. |
| `small_scale_error_file` | `Union{String, Missing}` | Small-scale error report path. |
| `corrected_assembly` | `Union{String, Missing}` | Corrected assembly path when correction is requested. |
| `params_json` | `String` | JSON-encoded parameters. |

Missing dependency behavior:

- If Inspector is unavailable, throw an error that includes
  `Inspector unavailable` and the executable that failed.
- Validate `datatype` against supported long-read values before execution.
- Default tests should cover invalid datatype, missing inputs, and parser
  fixtures. External tests should use tiny synthetic long-read fixtures and be
  gated by `MYCELIA_RUN_EXTERNAL=true`.

## Benchmark Orchestrator Contract

Planned high-level entry point:

```julia
Mycelia.run_assembly_benchmark_evaluation(assemblies::Vector{String};
    references::Vector{String}=String[],
    short_reads::Vector{String}=String[],
    long_reads::Vector{String}=String[],
    busco_lineage::Union{String, Nothing}=nothing,
    outdir::String="assembly_benchmark_evaluation",
    tools::Vector{Symbol}=[:quast, :busco, :merqury, :inspector],
    threads::Int=Mycelia.get_default_threads())
```

Planned output tables:

| Table | One row per | Source tools | Purpose |
| --- | --- | --- | --- |
| `benchmark_runs.tsv` | tool execution | all | Provenance, versions, parameters, status, output directories, runtime. |
| `assembly_metrics.tsv` | assembly | native metrics, QUAST | Contiguity and length metrics. |
| `reference_agreement.tsv` | assembly-reference pair | QUAST, MUMmer alternatives | Genome fraction, mismatch/indel rates, misassemblies, alignment coverage. |
| `completeness.tsv` | assembly-lineage pair | BUSCO, CheckM/CheckM2 alternatives | Conserved-gene or marker completeness and contamination where applicable. |
| `kmer_accuracy.tsv` | assembly-read-k tuple | Merqury, native QV | QV, error rate, k-mer recovery, completeness. |
| `long_read_error_support.tsv` | assembly-long-read pair | Inspector | Structural and small-scale error counts, read mapping support. |
| `benchmark_summary.tsv` | assembly | all configured tools | Wide summary for manuscript tables, with missing values and applicability flags retained. |

The orchestrator should not coerce unavailable metrics to zero. Missing values
must mean "not run" or "not applicable" and should be accompanied by a
`status`/`reason` column in `benchmark_runs.tsv`.

## Missing Dependency and Test Policy

- Default CI: parser fixtures, input validation, argument conflict validation,
  and native QV tests only.
- External smoke CI: `MYCELIA_RUN_EXTERNAL=true` and small fixtures for QUAST,
  BUSCO, Merqury, and Inspector as wrappers are implemented.
- Long-running benchmark CI or HPC: full manuscript datasets, larger references,
  BUSCO downloads, and long-read Inspector runs.
- A skipped external test must state the missing environment variable or missing
  fixture. A failed enabled external run must throw an error with the tool name,
  executable/environment, and expected output path.

## Implementation Slices

1. Add parser fixtures and default tests for QUAST and BUSCO summaries.
2. Add `parse_quast_report`, `summarize_quast`,
   `parse_busco_short_summary`, `parse_busco_full_table`, and
   `summarize_busco`.
3. Add a shared dependency probe helper for external wrappers that reports
   tool-specific unavailable messages while preserving existing Bioconda
   behavior.
4. Add `run_merqury` plus Merqury output parsers; keep native QV fallback
   explicit in the output schema.
5. Add `run_inspector` plus Inspector summary/error parsers.
6. Add `run_assembly_benchmark_evaluation` only after individual normalized
   summary tables are stable.
7. Add docs examples and external smoke tests, split so CI can enable one tool
   family at a time.

## References

- QUAST: https://doi.org/10.1093/bioinformatics/btt086
- BUSCO: https://doi.org/10.1093/bioinformatics/btv351
- Merqury: https://doi.org/10.1186/s13059-020-02134-9
- Inspector: https://doi.org/10.1186/s13059-021-02527-4
