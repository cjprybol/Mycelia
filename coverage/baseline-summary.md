# Mycelia Coverage Baseline

Generated: 2026-04-07T08:51:02.847

Scope: 57 top-level modules spanning 85 loaded source files

## Overall: 35.1% (9266/26431 executable lines)

## Tier Summary

| Tier | Name | Coverage | Target | Gap | Lines | Modules |
|------|------|----------|--------|-----|-------|---------|
| 1 | Data Acquisition | 26.9% | 90% | 63.1pp | 711/2641 | 5 |
| 2 | Preprocessing/QC | 39.6% | 90% | 50.4pp | 634/1599 | 5 |
| 3 | K-mer Analysis | 54.7% | 95% | 40.3pp | 1073/1962 | 5 |
| 4 | Assembly Core | 67.0% | 95% | 28.0pp | 4962/7403 | 8 |
| 5 | Validation | 10.0% | 80% | 70.0pp | 154/1533 | 3 |
| 6 | Annotation | 0.0% | 80% | 80.0pp | 0/1268 | 3 |
| 7 | Comparative | 6.2% | 85% | 78.8pp | 153/2483 | 5 |
| 8 | Tool Integration | 0.0% | 75% | 75.0pp | 0/3518 | 9 |
| 9 | Analysis/Viz | 35.1% | 85% | 49.9pp | 865/2463 | 5 |
| 10 | Infrastructure | 45.7% | 90% | 44.3pp | 714/1561 | 3 |

## Module Detail

| Module | Tier | Coverage | Target | Lines Hit | Lines Total | Source Files | Excluded |
|--------|------|----------|--------|-----------|-------------|--------------|----------|
| src/busco-datasets.jl | Data Acquisition | 0.0% | 90% | 0 | 28 | 1 |  |
| src/ncbi-datasets-cli.jl | Data Acquisition | 60.5% | 90% | 147 | 243 | 1 |  |
| src/protein-databases.jl | Data Acquisition | 0.0% | 90% | 0 | 708 | 1 |  |
| src/reference-databases.jl | Data Acquisition | 21.0% | 90% | 221 | 1053 | 1 |  |
| src/simulation.jl | Data Acquisition | 56.3% | 90% | 343 | 609 | 1 |  |
| src/alphabets.jl | Preprocessing/QC | 98.0% | 90% | 49 | 50 | 1 |  |
| src/constants.jl | Preprocessing/QC | 100.0% | 90% | 7 | 7 | 1 |  |
| src/fastx.jl | Preprocessing/QC | 57.6% | 90% | 568 | 986 | 1 |  |
| src/read-quality-control.jl | Preprocessing/QC | 15.9% | 90% | 10 | 63 | 1 |  |
| src/xam.jl | Preprocessing/QC | 0.0% | 90% | 0 | 493 | 1 |  |
| src/coverage-clustering.jl | K-mer Analysis | 96.1% | 95% | 147 | 153 | 1 |  |
| src/distance-metrics.jl | K-mer Analysis | 28.5% | 95% | 108 | 379 | 1 |  |
| src/kmer-analysis.jl | K-mer Analysis | 53.1% | 95% | 543 | 1023 | 1 |  |
| src/kmer-saturation-analysis.jl | K-mer Analysis | 87.3% | 95% | 178 | 204 | 1 |  |
| src/qualmer-analysis.jl | K-mer Analysis | 47.8% | 95% | 97 | 203 | 1 |  |
| src/assembly.jl | Assembly Core | 6.8% | 95% | 57 | 840 | 1 |  |
| src/autocycler.jl | Assembly Core | 0.0% | 95% | 0 | 42 | 1 |  |
| src/bcalm.jl | Assembly Core | 0.0% | 95% | 0 | 55 | 1 |  |
| src/ggcat.jl | Assembly Core | 0.0% | 95% | 0 | 81 | 1 |  |
| src/graph-cleanup.jl | Assembly Core | 54.7% | 95% | 129 | 236 | 1 |  |
| src/iterative-assembly.jl | Assembly Core | 66.8% | 95% | 534 | 799 | 1 |  |
| src/rhizomorph/rhizomorph.jl | Assembly Core | 79.3% | 95% | 4242 | 5350 | 29 |  |
| src/viterbi-next.jl | Assembly Core | 0.0% | 95% | 0 | 0 | 1 |  |
| src/quality-control-and-benchmarking.jl | Validation | 19.1% | 80% | 154 | 805 | 1 |  |
| src/variant-analysis.jl | Validation | 0.0% | 80% | 0 | 427 | 1 |  |
| src/viterbi-polishing-and-error-correction.jl | Validation | 0.0% | 80% | 0 | 301 | 1 |  |
| src/annotation.jl | Annotation | 0.0% | 80% | 0 | 1011 | 1 |  |
| src/codon-optimization.jl | Annotation | 0.0% | 80% | 0 | 111 | 1 |  |
| src/genome-features.jl | Annotation | 0.0% | 80% | 0 | 146 | 1 |  |
| src/metagraph.jl | Comparative | 0.0% | 85% | 0 | 377 | 1 |  |
| src/pangenome-analysis.jl | Comparative | 0.0% | 85% | 0 | 246 | 1 |  |
| src/relational-matrices.jl | Comparative | 0.0% | 85% | 0 | 189 | 1 |  |
| src/sequence-comparison.jl | Comparative | 14.4% | 85% | 153 | 1065 | 1 |  |
| src/taxonomy-and-trees.jl | Comparative | 0.0% | 85% | 0 | 606 | 1 |  |
| src/alignments-and-mapping.jl | Tool Integration | 0.0% | 75% | 0 | 1424 | 1 |  |
| src/binning.jl | Tool Integration | 0.0% | 75% | 0 | 260 | 1 |  |
| src/bioconda.jl | Tool Integration | 20.8% | 75% | 27 | 130 | 1 | yes |
| src/classification.jl | Tool Integration | 0.0% | 75% | 0 | 615 | 1 |  |
| src/foldseek.jl | Tool Integration | 0.0% | 75% | 0 | 121 | 1 |  |
| src/metagenomic-classification.jl | Tool Integration | 0.0% | 75% | 0 | 192 | 1 |  |
| src/pantools.jl | Tool Integration | 0.0% | 75% | 0 | 346 | 1 |  |
| src/prokrustean.jl | Tool Integration | 0.0% | 75% | 0 | 210 | 1 |  |
| src/rclone.jl | Tool Integration | 0.0% | 75% | 0 | 124 | 1 |  |
| src/sentencepiece.jl | Tool Integration | 0.0% | 75% | 0 | 226 | 1 |  |
| src/amino-acid-analysis.jl | Analysis/Viz | 95.2% | 85% | 20 | 21 | 1 |  |
| src/clustering.jl | Analysis/Viz | 60.9% | 85% | 403 | 662 | 1 |  |
| src/dimensionality-reduction.jl | Analysis/Viz | 66.4% | 85% | 85 | 128 | 1 |  |
| src/plotting-and-visualization.jl | Analysis/Viz | 18.6% | 85% | 296 | 1588 | 1 |  |
| src/tda.jl | Analysis/Viz | 95.3% | 85% | 61 | 64 | 1 |  |
| src/Mycelia.jl | Infrastructure | 100.0% | 90% | 8 | 8 | 1 | yes |
| src/checkpointing.jl | Infrastructure | 100.0% | 90% | 76 | 76 | 1 |  |
| src/execution.jl | Infrastructure | 0.0% | 90% | 0 | 59 | 1 | yes |
| src/precompile_workload.jl | Infrastructure | 3.4% | 90% | 1 | 29 | 1 | yes |
| src/slurm-sbatch.jl | Infrastructure | 0.0% | 90% | 0 | 98 | 1 |  |
| src/slurm-templates.jl | Infrastructure | 0.0% | 90% | 0 | 953 | 1 | yes |
| src/testing-utilities.jl | Infrastructure | 28.0% | 90% | 47 | 168 | 1 | yes |
| src/utility-functions.jl | Infrastructure | 46.0% | 90% | 638 | 1387 | 1 |  |
