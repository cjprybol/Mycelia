# Caveat: `genome_fraction` in committed real_genome_benchmark_*.csv files

The `genome_fraction` column in the committed
`benchmarking/results/real_genome_benchmark_*.csv` files is
`total_length / ref_size * 100`. That is an assembly **size_ratio**, not
QUAST's aligned reference coverage. A fragmented or over-assembled
result can exceed 100% (viroids in these tables reach ~200% when two
full-length contigs are emitted).

`benchmarking/t4_ksweep.jl` now names the live column `size_ratio` and
picks best-k by `n50 / n_contigs`, not `argmax(genome_fraction)`. The
committed CSVs are historical evidence and were not rewritten; do not
cite their `genome_fraction` column as aligned coverage, and do not
treat an argmax of that column as the most accurate k.
