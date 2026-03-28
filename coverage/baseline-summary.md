# Mycelia Coverage Baseline

Generated: 2026-03-28T13:55:30.140

## Overall: 82.3% (9219/11199 executable lines)

## Tier Summary

| Tier | Name | Coverage | Target | Gap | Lines | Modules |
|------|------|----------|--------|-----|-------|---------|
| 1 | Data Acquisition | 66.0% | 90% | 24.0pp | 167/253 | 2 |
| 2 | Preprocessing/QC | 75.0% | 90% | 15.0pp | 300/400 | 3 |
| 3 | K-mer Analysis | 75.9% | 95% | 19.1pp | 1081/1425 | 4 |
| 4 | Assembly Core | 83.7% | 95% | 11.3pp | 4883/5832 | 29 |
| 5 | Validation | 95.2% | 80% | -15.2pp | 99/104 | 1 |
| 6 | Annotation | 92.9% | 80% | -12.9pp | 221/238 | 2 |
| 7 | Comparative | 76.7% | 85% | 8.3pp | 353/460 | 2 |
| 8 | Tool Integration | 87.1% | 75% | -12.1pp | 27/31 | 1 |
| 9 | Analysis/Viz | 85.2% | 85% | -0.2pp | 861/1010 | 3 |
| 10 | Infrastructure | 82.1% | 90% | 7.9pp | 546/665 | 2 |

## Module Detail

| Module | Tier | Coverage | Target | Lines Hit | Lines Total | Excluded |
|--------|------|----------|--------|-----------|-------------|----------|
| Mycelia.jl | Unassigned | 87.5% | 0% | 7 | 8 | yes |
| amino-acid-analysis.jl | Unassigned | 95.2% | 0% | 20 | 21 |  |
| bioconda.jl | Unassigned | 0.0% | 0% | 0 | 44 | yes |
| codon-optimization.jl | Unassigned | 98.6% | 0% | 73 | 74 |  |
| coverage-clustering.jl | Unassigned | 96.1% | 0% | 147 | 153 |  |
| precompile_workload.jl | Unassigned | 3.4% | 0% | 1 | 29 | yes |
| relational-matrices.jl | Unassigned | 77.5% | 0% | 141 | 182 |  |
| taxonomy-and-trees.jl | Unassigned | 77.3% | 0% | 116 | 150 |  |
| tda.jl | Unassigned | 95.3% | 0% | 61 | 64 |  |
| testing-utilities.jl | Unassigned | 97.0% | 0% | 32 | 33 | yes |
| variant-analysis.jl | Unassigned | 100.0% | 0% | 42 | 42 |  |
| viterbi-polishing-and-error-correction.jl | Unassigned | 85.3% | 0% | 81 | 95 |  |
| reference-databases.jl | Data Acquisition | 56.6% | 90% | 47 | 83 |  |
| simulation.jl | Data Acquisition | 70.6% | 90% | 120 | 170 |  |
| alphabets.jl | Preprocessing/QC | 98.0% | 90% | 49 | 50 |  |
| constants.jl | Preprocessing/QC | 100.0% | 90% | 7 | 7 |  |
| fastx.jl | Preprocessing/QC | 71.1% | 90% | 244 | 343 |  |
| distance-metrics.jl | K-mer Analysis | 89.1% | 95% | 310 | 348 |  |
| kmer-analysis.jl | K-mer Analysis | 65.5% | 95% | 496 | 757 |  |
| kmer-saturation-analysis.jl | K-mer Analysis | 87.3% | 95% | 178 | 204 |  |
| qualmer-analysis.jl | K-mer Analysis | 83.6% | 95% | 97 | 116 |  |
| assembly.jl | Assembly Core | 67.7% | 95% | 327 | 483 |  |
| contigs.jl | Assembly Core | 83.5% | 95% | 111 | 133 |  |
| edge-data.jl | Assembly Core | 100.0% | 95% | 23 | 23 |  |
| error-correction.jl | Assembly Core | 95.6% | 95% | 43 | 45 |  |
| evidence-functions.jl | Assembly Core | 96.3% | 95% | 287 | 298 |  |
| evidence-structures.jl | Assembly Core | 100.0% | 95% | 19 | 19 |  |
| fasta-graphs.jl | Assembly Core | 87.5% | 95% | 63 | 72 |  |
| fastq-graphs.jl | Assembly Core | 78.5% | 95% | 62 | 79 |  |
| generation.jl | Assembly Core | 88.5% | 95% | 116 | 131 |  |
| graph-cleanup.jl | Assembly Core | 83.2% | 95% | 129 | 155 |  |
| graph-construction.jl | Assembly Core | 87.3% | 95% | 1110 | 1272 |  |
| graph-query.jl | Assembly Core | 89.2% | 95% | 165 | 185 |  |
| graph-type-conversions.jl | Assembly Core | 81.8% | 95% | 126 | 154 |  |
| information-theory.jl | Assembly Core | 89.5% | 95% | 68 | 76 |  |
| io.jl | Assembly Core | 84.9% | 95% | 186 | 219 |  |
| iterative-assembly.jl | Assembly Core | 71.0% | 95% | 532 | 749 |  |
| kmer-graphs.jl | Assembly Core | 80.5% | 95% | 153 | 190 |  |
| metrics.jl | Assembly Core | 95.1% | 95% | 58 | 61 |  |
| ngram-graphs.jl | Assembly Core | 94.4% | 95% | 85 | 90 |  |
| path-finding.jl | Assembly Core | 85.8% | 95% | 319 | 372 |  |
| quality-functions.jl | Assembly Core | 96.1% | 95% | 74 | 77 |  |
| qualmer-graphs.jl | Assembly Core | 73.7% | 95% | 56 | 76 |  |
| repeats.jl | Assembly Core | 95.6% | 95% | 151 | 158 |  |
| rhizomorph.jl | Assembly Core | 100.0% | 95% | 1 | 1 |  |
| sequence-quality.jl | Assembly Core | 95.1% | 95% | 97 | 102 |  |
| simplification.jl | Assembly Core | 85.2% | 95% | 283 | 332 |  |
| strand-conversions.jl | Assembly Core | 61.7% | 95% | 58 | 94 |  |
| string-graphs.jl | Assembly Core | 97.1% | 95% | 133 | 137 |  |
| vertex-data.jl | Assembly Core | 98.0% | 95% | 48 | 49 |  |
| quality-control-and-benchmarking.jl | Validation | 95.2% | 80% | 99 | 104 |  |
| annotation.jl | Annotation | 95.7% | 80% | 154 | 161 |  |
| genome-features.jl | Annotation | 87.0% | 80% | 67 | 77 |  |
| pangenome-analysis.jl | Comparative | 56.7% | 85% | 115 | 203 |  |
| sequence-comparison.jl | Comparative | 92.6% | 85% | 238 | 257 |  |
| alignments-and-mapping.jl | Tool Integration | 87.1% | 75% | 27 | 31 |  |
| clustering.jl | Analysis/Viz | 86.5% | 85% | 467 | 540 |  |
| dimensionality-reduction.jl | Analysis/Viz | 89.5% | 85% | 85 | 95 |  |
| plotting-and-visualization.jl | Analysis/Viz | 82.4% | 85% | 309 | 375 |  |
| checkpointing.jl | Infrastructure | 100.0% | 90% | 76 | 76 |  |
| utility-functions.jl | Infrastructure | 79.8% | 90% | 470 | 589 |  |
