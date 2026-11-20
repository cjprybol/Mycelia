# ==============================================================================
# Shared Stage-0 harness (td-f581, PR #347 review I1).
#
# stage0_classification_matrix.jl and stage0_threshold_calibration.jl both
# include this file. It is the one copy of the read simulator, the genomic
# k-mer truth set, and the qualmer memory_profile, so those cannot drift apart.
#
# Include it only after the caller has defined the run configuration:
#   REFLEN, COVERAGE, READLEN, ERROR_RATE, K
# Those stay in each script (the matrix defaults to a 2 kb reference;
# calibration defaults to 10 kb). The simulation model itself is here.
# ==============================================================================

import FASTX
import BioSequences
import Kmers

# :full keeps every per-observation quality value, which is what
# get_vertex_mean_quality needs for the depth-independent per-base mean the
# fused arms use (review C1/F2). :lightweight_quality stores only accumulated
# joint quality and would zero that signal.
const STAGE0_MEMORY_PROFILE = :full

const BASES = ['A', 'C', 'G', 'T']

"""Phred quality model: correct bases ~Q35, error bases ~Q22 (lower but
overlapping), so quality carries PARTIAL signal — not a perfect oracle. This is
what lets the matrix show where fusing quality with coverage actually helps."""
function sample_quality(rng, is_error::Bool)::Char
    q = is_error ? round(Int, 22 + 8 * randn(rng)) : round(Int, 35 + 3 * randn(rng))
    q = clamp(q, 2, 40)
    return Char(q + 33)  # Phred+33
end

"""Draw an error-free reference and substitution-noised FASTQ reads.

The draw uses the caller's `REFLEN`, `COVERAGE`, `READLEN`, and `ERROR_RATE`.
`rng` is the only argument so both Stage-0 scripts share this model.
"""
function simulate(rng)
    ref = join(rand(rng, BASES, REFLEN))
    n_reads = max(1, round(Int, COVERAGE * REFLEN / READLEN))
    recs = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        start = rand(rng, 1:(REFLEN - READLEN + 1))
        seq = collect(ref[start:(start + READLEN - 1)])
        quals = Vector{Char}(undef, READLEN)
        for j in 1:READLEN
            is_err = rand(rng) < ERROR_RATE
            if is_err
                alt = rand(rng, filter(!=(seq[j]), BASES))
                seq[j] = alt
            end
            quals[j] = sample_quality(rng, is_err)
        end
        push!(recs, FASTX.FASTQ.Record("r$i", join(seq), join(quals)))
    end
    return ref, recs
end

"""Truth set: canonical k-mers of the ERROR-FREE reference = genomic k-mers."""
function reference_kmer_set(ref::AbstractString)
    truth = Set{Kmers.DNAKmer{K}}()
    seq = BioSequences.LongDNA{4}(ref)
    for (kmer, _) in Kmers.UnambiguousDNAMers{K}(seq)
        push!(truth, BioSequences.canonical(kmer))
    end
    return truth
end
