# Viterbi accuracy benchmark fixtures

Small committed real-data fixtures for `benchmarking/viterbi_accuracy_benchmark.jl`.

- `pstvd_nc002030.fasta`: Potato spindle tuber viroid, complete genome, NCBI RefSeq `NC_002030.1`.
- `phix174_nc001422.fasta`: Escherichia phage phiX174, complete genome, NCBI RefSeq `NC_001422.1`.
- `pride_and_prejudice_excerpt.txt`: Public-domain Project Gutenberg-style excerpt from Jane Austen's *Pride and Prejudice* used as a text-corpus fixture.

The benchmark uses short deterministic windows/chunks from these fixtures so local smoke runs produce committed B8 artifacts without requiring network or HPC access.
