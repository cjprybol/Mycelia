# Workflow coverage and documentation gaps

See `planning-docs/FUNCTION_COVERAGE_AUDIT.md` for per-module coverage counts and doc link audits.

## Auto-included surface area (from `src/Mycelia.jl`)

- **Core graph + IO**: `utility-functions`, `alphabets`, `constants`, `fastx`, `graph-core`, `sequence-graphs-next`, `string-graphs`, `qualmer-analysis/graphs`, `fasta-graphs`, `fastq-graphs`.
  - Owner: graph core team (Rhizomorph maintainers).
- **Assembly**: `assembly.jl` (external assemblers), `iterative-assembly.jl` (native likelihood improver), `viterbi-next.jl`.
  - Owner: assembly pipeline.
- **Analytics/QC**: `quality-control-and-benchmarking.jl`, `performance-benchmarks.jl`, `kmer-analysis.jl`, `distance-metrics.jl`.
  - Owner: QC + benchmarking.
- **Taxonomy/annotation**: `taxonomy-and-trees.jl`, `classification.jl`, `reference-databases.jl`, `annotation.jl`.
  - Owner: taxonomy/annotation.
- **Wrappers/orchestration**: `bioconda.jl`, `rclone.jl`, `slurm-sbatch.jl`, `neo4jl.jl`, `xam.jl`.
  - Owner: tooling/integration.

## Documented vs. undocumented entry points

- Documented: external assemblers (`run_megahit`, `run_metaspades`, `run_spades`, `run_skesa`) and QC helpers (`robust_cv`, `robust_threshold`, `filter_genome_outliers`) already had docstrings; now referenced on the workflow map.
- Newly surfaced: iterative assembly helpers (`mycelia_iterative_assemble`, `improve_read_set_likelihood`, `find_optimal_sequence_path`) were missing from docs/tutorials and are now summarized on the map with proposed tutorial coverage.
- Still sparse: Rhizomorph path-finding/simplification utilities have docstrings in-file but no tutorial; taxonomy assignment utilities lack an end-to-end example.

## Tutorial coverage findings

- Present: acquisition (01), QC (02), k-mer analysis (03), assembly (04), validation (05), annotation/comparative genomics (06–08), graph round-trips (09 series), advanced assembly theory.
- Missing:
  - Rhizomorph graph experiments (k-mer/qualmer/string graph builds, path-finding, simplification).
  - Taxonomy/reference database workflow end-to-end.
  - QC outlier filtering example using `filter_genome_outliers` and `robust_cv`.
  - Iterative likelihood-improvement handoff after external assembly.
  - Variant calling example (VCF normalization + evaluation).
  - Metagenomic binning + profiling paths.
  - Sketch-guided pangenome context selection (sourmash/sylph/mash prefilter -> subset reference paths -> minimap2 mapping).

## Proposed actions

- Add `tutorials/12_rhizomorph_graphs.jl` with small FASTA/FASTQ fixtures from `assembly_test_data/` to exercise graph building, simplification, Eulerian paths, and GFA I/O.
- Extend `tutorials/02_quality_control.jl` with a short QC-outlier filtering example (`filter_genome_outliers`, `robust_cv`, `robust_threshold`).
- Add `tutorials/13_taxonomy_assignment.jl` using a tiny mock reference FASTA + TSV taxonomy to showcase database fetch + abundance summarization.
- Add a variant-calling tutorial demonstrating `run_gatk_haplotypecaller` → `normalize_vcf` → `evaluate_variant_calling_accuracy`.
- Add a metagenomics tutorial tying together metagenomic assembly, binning, and taxonomic profiling.
