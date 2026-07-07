# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/batched_viterbi_kernel_test.jl")'
# ```
#
# (`import KernelAbstractions` below triggers the MyceliaKernelAbstractionsExt
# extension; KernelAbstractions is in the package's test `[extras]`.)
#
# Equivalence-oracle acceptance gate for the Phase B backend-agnostic frontier
# KERNEL of the batched array-frontier corrector (bead td-qoo3; ADR
# docs/design/2026-07-06-gpu-simd-corrector-acceleration.md).
#
# The kernel advances the dense `n_reads x n_states` frontier one depth-step for
# every read in a length-bin in parallel (one thread per read) over a shared,
# read-only CSR transition SoA + precomputed emission table. This test runs the
# kernel on the CPU KernelAbstractions backend (so CI needs no GPU) and asserts
# its corrected paths + log-probability scores are BYTE-IDENTICAL, per read, to
# BOTH the scalar `correct_observations` AND the merged CPU PoC
# `batched_correct_observations`. It is also the regression guard for the GPU
# path: the same kernel runs on a CUDA backend on GPU hosts.

import BioSequences
import FASTX
import KernelAbstractions
import Kmers
import Mycelia
import Test

function poc_decoded_label_strings(
    result::Mycelia.ViterbiCorrectionResult, index::Integer
)::Vector{String}
    return [string(label) for label in result.corrected_observations[index]]
end

Test.@testset "Batched Viterbi corrector GPU kernel — Phase B (td-qoo3)" begin
    cpu_backend = KernelAbstractions.CPU()

    Test.@testset "n-gram (TEXT, singlestrand) batch — kernel equivalence oracle" begin
        graph = Mycelia.Rhizomorph.build_ngram_graph(
            ["ABCDEFGH"], 3; dataset_id="kernel_poc_ngram"
        )
        observations = [
            ["ABC", "BCX", "CDE", "DEF"],
            ["DEF", "EFG", "FGH"],
            ["ABC", "BCD", "CDX", "DEF"],
            ["DEF", "EFZ", "FGH"],
            ["CDE", "DEF", "EFG", "FGH"],
        ]

        oracle = Mycelia.kernel_batched_equivalence_oracle(
            graph, observations; backend=cpu_backend)

        # THE acceptance gate: kernel == scalar AND kernel == CPU PoC, byte-identical.
        Test.@test oracle.passes
        Test.@test all(oracle.per_read)
        Test.@test all(oracle.per_read_vs_scalar)
        Test.@test all(oracle.per_read_vs_poc)
        Test.@test length(oracle.per_read) == length(observations)

        # The error reads were corrected back onto the reference chain.
        Test.@test poc_decoded_label_strings(oracle.kernel, 1) ==
            ["ABC", "BCD", "CDE", "DEF"]
        Test.@test poc_decoded_label_strings(oracle.kernel, 3) ==
            ["ABC", "BCD", "CDE", "DEF"]
        Test.@test poc_decoded_label_strings(oracle.kernel, 4) ==
            ["DEF", "EFG", "FGH"]

        for index in 1:length(observations)
            Test.@test poc_decoded_label_strings(oracle.kernel, index) ==
                poc_decoded_label_strings(oracle.scalar, index)
            Test.@test oracle.kernel.paths[index].score ===
                oracle.scalar.paths[index].score
            Test.@test oracle.kernel.paths[index].score ===
                oracle.poc.paths[index].score
        end

        Test.@test oracle.kernel.diagnostics[:interface] == :metagraphs_next_kernel
        Test.@test oracle.kernel.diagnostics[:length_bins] == [3, 4]
        Test.@test oracle.kernel.diagnostics[:strand_mode] == :singlestrand
    end

    Test.@testset "DNA k-mer (singlestrand) batch — kernel equivalence oracle" begin
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGTAC")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id="kernel_poc_kmer", mode=:singlestrand
        )
        observations = [
            [
                Kmers.DNAKmer{3}("ATG"),
                Kmers.DNAKmer{3}("TGA"),
                Kmers.DNAKmer{3}("GCG"),
                Kmers.DNAKmer{3}("CGT"),
            ],
            [
                Kmers.DNAKmer{3}("GCG"),
                Kmers.DNAKmer{3}("CGT"),
                Kmers.DNAKmer{3}("GTA"),
            ],
        ]

        oracle = Mycelia.kernel_batched_equivalence_oracle(
            graph, observations; backend=cpu_backend)

        Test.@test oracle.passes
        Test.@test all(oracle.per_read)
        Test.@test poc_decoded_label_strings(oracle.kernel, 1) ==
            ["ATG", "TGC", "GCG", "CGT"]
        for index in 1:length(observations)
            Test.@test poc_decoded_label_strings(oracle.kernel, index) ==
                poc_decoded_label_strings(oracle.scalar, index)
            Test.@test oracle.kernel.paths[index].score ===
                oracle.scalar.paths[index].score
        end
    end

    Test.@testset "single-read batch reproduces the scalar decoder" begin
        graph = Mycelia.Rhizomorph.build_ngram_graph(
            ["ABCDEF"], 3; dataset_id="kernel_poc_single"
        )
        observations = [["ABC", "BCX", "CDE", "DEF"]]
        oracle = Mycelia.kernel_batched_equivalence_oracle(
            graph, observations; backend=cpu_backend)
        Test.@test oracle.passes
        Test.@test poc_decoded_label_strings(oracle.kernel, 1) ==
            ["ABC", "BCD", "CDE", "DEF"]
    end

    Test.@testset "kernel guards reject out-of-scope configs" begin
        graph = Mycelia.Rhizomorph.build_ngram_graph(
            ["ABCDEF"], 3; dataset_id="kernel_poc_guard"
        )
        observations = [["ABC", "BCD", "CDE", "DEF"]]
        beam_config = Mycelia.ViterbiCorrectionConfig(beam_width=2)
        Test.@test_throws ArgumentError Mycelia.kernel_batched_correct_observations(
            graph, observations; config=beam_config, backend=cpu_backend
        )
    end

    Test.@testset "benchmark harness runs (CPU backend; GPU if available)" begin
        # Small sizes keep the CI run fast; the harness reports GPU throughput
        # only when a functional CUDA backend is detected (e.g. on Lovelace).
        bench = Mycelia.kernel_batched_benchmark(
            n_reads=256, n_states=64, avg_degree=3, n_depths=8, iterations=1)
        Test.@test bench.cpu_edge_relaxations_per_sec > 0
        Test.@test haskey(bench, :cuda_available)
        if bench.cuda_available
            Test.@test bench.results_match
            Test.@test bench.gpu_edge_relaxations_per_sec > 0
        end
    end
end
