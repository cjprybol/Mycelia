# Reusable decoder-benchmark k-mer/n-gram graph fixture.
# ======================================================
#
# ONE vetted graph-build + observation-conversion path, extracted from
# viterbi_accuracy_benchmark.jl (B8), so the Tier-2 gap calibration (td-4osf) and
# future decoder benchmarks stop each re-deriving how to build a correction graph
# from a truth sequence and how to convert reads into the label-observation form
# `Mycelia.correct_observations` expects.

import Mycelia
import FASTX
import Kmers

# Sliding-window k-mers of a sequence (byte-indexed — ASCII alphabets only, which
# is what the DNA/RNA/text decoder fixtures use). One string per window position.
_bench_sequence_kmers(sequence::AbstractString,
    k::Int)::Vector{String} = [String(sequence[i:(i + k - 1)])
                               for i in 1:(lastindex(sequence) - k + 1)]

# Convert a k-mer-string sequence into the label-observation form the decoder
# expects for this molecule type (DNA/RNA -> typed Kmers; text -> raw strings).
function _bench_convert_kmers(kmer_strings::Vector{String}, k::Int, moltype::Symbol)
    moltype === :DNA && return [Kmers.DNAKmer{k}(s) for s in kmer_strings]
    moltype === :RNA && return [Kmers.RNAKmer{k}(s) for s in kmer_strings]
    moltype === :text && return kmer_strings
    throw(ArgumentError("unsupported moltype $moltype (expected :DNA, :RNA, :text)"))
end

"""
    benchmark_kmer_graph_fixture(sequence, k; moltype=:DNA, dataset_id="benchmark")

Build a decoder-benchmark correction graph from a single truth `sequence`, plus a
matching sequence -> observation converter. Returns a NamedTuple:

  * `graph`             — the k-mer (DNA/RNA, singlestrand) or n-gram (text)
                          correction graph, ready for `Mycelia.correct_observations`.
  * `truth_kmers`       — `Vector{String}` of the truth's sliding-window k-mers.
  * `to_observation(seq)` — `seq::AbstractString` -> the label-observation vector
                          such that `correct_observations(graph, [to_observation(seq)])`
                          runs. Pass the truth to get the clean observation, or an
                          errored read to get one for correction.

This is the extraction of B8's `build_kmer_graph`/`build_ngram_graph` +
`_sequence_kmers` + `_convert_observations` into a single vetted entry point.
"""
function benchmark_kmer_graph_fixture(sequence::AbstractString, k::Int;
        moltype::Symbol = :DNA, dataset_id::AbstractString = "benchmark")
    graph = if moltype === :DNA
        Mycelia.Rhizomorph.build_kmer_graph(
            [FASTX.FASTA.Record(dataset_id, sequence)], k;
            dataset_id = dataset_id, mode = :singlestrand)
    elseif moltype === :RNA
        Mycelia.Rhizomorph.build_kmer_graph(
            [FASTX.FASTA.Record(dataset_id, sequence)], k;
            dataset_id = dataset_id, mode = :singlestrand, type_hint = :RNA)
    elseif moltype === :text
        Mycelia.Rhizomorph.build_ngram_graph([String(sequence)], k; dataset_id = dataset_id)
    else
        throw(ArgumentError("unsupported moltype $moltype (expected :DNA, :RNA, :text)"))
    end
    to_observation = seq -> _bench_convert_kmers(_bench_sequence_kmers(seq, k), k, moltype)
    return (graph = graph,
        truth_kmers = _bench_sequence_kmers(sequence, k),
        to_observation = to_observation)
end
