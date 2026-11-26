# Pin: the indel pair-HMM match move emits one new base, not the whole k-mer.
# Written against the decoder-repair series (landed in #497, per-base emission).
# The whole-k-mer helper `_call_viterbi_state_emission_logp` and its
# `count_length_penalty` keyword are still on master, so the old-vs-new
# comparison below is kept.
import Test
import Mycelia
import FASTX
import BioSequences
import Random
import Kmers

Test.@testset "indel match move emits one base, not k" begin
    R = Mycelia.Rhizomorph
    k = 17
    Q = 15
    rng = Random.MersenneTwister(11)
    G = BioSequences.randdnaseq(rng, 200)
    gs = string(G)
    qs(n) = repeat(Char(Q + 33), n)
    reads = FASTX.FASTQ.Record[]
    for i in 1:60
        s = rand(rng, 1:141)
        push!(reads, FASTX.FASTQ.Record("c$i", gs[s:s+59], qs(60)))
    end
    truth = gs[71:130]
    mutate(t, kind) = kind == :sub ? t[1:30] * (t[31] == 'A' ? "C" : "A") * t[32:end] :
                      kind == :ins ? t[1:30] * (t[31] == 'A' ? "C" : "A") * t[31:end] :
                      kind == :del ? t[1:30] * t[32:end] : t
    profile = Mycelia.indel_error_profile(:nanopore)
    ip = Mycelia.IndelDecodeParams(profile.base_error_rate, profile.insertion_fraction,
        profile.deletion_fraction, profile.insertion_extend_probability,
        profile.deletion_extend_probability, 3, 3, nothing)
    single_mismatch = log1p(-10^(-Q / 10)) - log(10^(-Q / 10) / 3)
    qobs(s, i) = Mycelia.QualityObservation(Kmers.DNAKmer{k}(s[i:i+k-1]), UInt8.(fill(Q, k)))
    function path_emissions(graph, cfg, observed, hidden)
        n = length(observed) - k + 1
        first = Mycelia._call_viterbi_state_emission_logp(graph, cfg, qobs(observed, 1),
            Kmers.DNAKmer{k}(hidden[1:k]), :DNA, :singlestrand; count_length_penalty = false)
        rest = 0.0
        charged = 0
        for i in 2:n
            e = Mycelia._call_viterbi_new_base_emission_logp(graph, cfg, qobs(observed, i),
                Kmers.DNAKmer{k}(hidden[i:i+k-1]), :DNA, :singlestrand)
            rest += e
            e_own = Mycelia._call_viterbi_new_base_emission_logp(graph, cfg, qobs(observed, i),
                Kmers.DNAKmer{k}(observed[i:i+k-1]), :DNA, :singlestrand)
            e != e_own && (charged += 1)
        end
        return first + rest, charged
    end
    o = mutate(truth, :sub)
    test = FASTX.FASTQ.Record("test", o, qs(length(o)))
    graph = R.build_qualmer_graph(vcat(reads, [test]), k; mode = :singlestrand)
    cfg = Mycelia.ViterbiCorrectionConfig(error_rate = profile.base_error_rate, indel_moves = true)
    own, _ = path_emissions(graph, cfg, o, o)
    tru, n_charged = path_emissions(graph, cfg, o, truth)
    gap = own - tru
    n = length(o) - k + 1
    old_own = sum(Mycelia._call_viterbi_state_emission_logp(graph, cfg, qobs(o, i),
        Kmers.DNAKmer{k}(o[i:i+k-1]), :DNA, :singlestrand; count_length_penalty = false) for i in 1:n)
    old_tru = sum(Mycelia._call_viterbi_state_emission_logp(graph, cfg, qobs(o, i),
        Kmers.DNAKmer{k}(truth[i:i+k-1]), :DNA, :singlestrand; count_length_penalty = false) for i in 1:n)
    old_gap = old_own - old_tru
    Test.@test n_charged == 1
    Test.@test isapprox(gap, single_mismatch; atol = 0.05)
    Test.@test abs(old_gap) > 10 * abs(single_mismatch)
    Test.@test abs(gap) < abs(old_gap) / 5
    for kind in (:sub, :ins, :del)
        obs = mutate(truth, kind)
        rec = FASTX.FASTQ.Record("t_$(kind)", obs, qs(length(obs)))
        g = R.build_qualmer_graph(vcat(reads, [rec]), k; mode = :singlestrand)
        vr = Mycelia.try_viterbi_path_improvement(rec, g, k; graph_mode = :singlestrand, indel_params = ip)
        Test.@test vr !== nothing
        got = FASTX.sequence(String, vr[1])
        Test.@test got == truth
    end
end
