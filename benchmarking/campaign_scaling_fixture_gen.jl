#!/usr/bin/env julia
# gen_fixture.jl — deterministic synthetic genome + short-read FASTQ generator.
#
# Self-contained (Base + Random + FASTX only, no Mycelia-internal functions) so
# it behaves identically regardless of which Mycelia checkout is later loaded
# for correction. Read IDs are deterministic sequential strings (NOT UUIDs) so
# repeated invocations with the same ENV produce a byte-identical FASTQ.
#
# ENV:
#   FX_GENOME_LEN   genome length in bp           (default 2000)
#   FX_READLEN      read length in bp             (default 150)
#   FX_COVERAGE     target fold coverage          (default 50)
#   FX_ERR          per-base substitution rate    (default 0.01)
#   FX_SEED         RNG seed                      (default 42)
#   FX_OUT          output FASTQ path              (required)

import Random
import FASTX

genome_len = parse(Int, get(ENV, "FX_GENOME_LEN", "2000"))
readlen = parse(Int, get(ENV, "FX_READLEN", "150"))
coverage = parse(Float64, get(ENV, "FX_COVERAGE", "50"))
err_rate = parse(Float64, get(ENV, "FX_ERR", "0.01"))
seed = parse(Int, get(ENV, "FX_SEED", "42"))
outpath = ENV["FX_OUT"]

genome_len > 0 || error("FX_GENOME_LEN must be positive, got $genome_len")
readlen > 0 || error("FX_READLEN must be positive, got $readlen")
coverage > 0 || error("FX_COVERAGE must be positive, got $coverage")
(0.0 <= err_rate <= 1.0) || error("FX_ERR must be in [0, 1], got $err_rate")

rng = Random.MersenneTwister(seed)
bases = ('A', 'C', 'G', 'T')

# 1. Deterministic random genome.
genome = join(rand(rng, collect(bases), genome_len))

function revcomp(s::AbstractString)
    comp = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')
    return String(reverse([comp[c] for c in s]))
end

function inject_errors(s::AbstractString, rate::Float64, rng::Random.AbstractRNG)
    chars = collect(s)
    for i in eachindex(chars)
        if rand(rng) < rate
            alts = filter(!=(chars[i]), collect(bases))
            chars[i] = alts[rand(rng, 1:length(alts))]
        end
    end
    return String(chars)
end

effective_readlen = min(readlen, genome_len)
effective_readlen < readlen &&
    @warn "FX_READLEN ($readlen) > FX_GENOME_LEN ($genome_len); reads clipped to $effective_readlen bp"
n_reads = max(1, ceil(Int, coverage * genome_len / effective_readlen))

open(outpath, "w") do io
    writer = FASTX.FASTQ.Writer(io)
    for i = 1:n_reads
        start = rand(rng, 1:(genome_len-effective_readlen+1))
        frag = genome[start:(start+effective_readlen-1)]
        strand = rand(rng, (:fwd, :rev))
        frag = strand == :fwd ? frag : revcomp(frag)
        observed = inject_errors(frag, err_rate, rng)
        qual = String(fill(Char(20 + 33), length(observed)))  # uniform Q20
        rec = FASTX.FASTQ.Record("read_$(i)", observed, qual)
        write(writer, rec)
    end
    close(writer)
end

println(
    "genome_len=$genome_len readlen=$readlen effective_readlen=$effective_readlen coverage=$coverage err_rate=$err_rate n_reads=$n_reads seed=$seed -> $outpath",
)
