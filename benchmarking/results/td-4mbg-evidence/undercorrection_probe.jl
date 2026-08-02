# code-audit I1: does the truncated partition UNDER-CORRECT?
#
# Contiguity is a downstream proxy. The direct question is per-base: after
# correction, how close is each read to the reference it came from? If truncating
# the partition costs real corrections, the fix's corrected reads will sit at
# LOWER identity than master's on the same input.
#
# Both code versions get byte-identical reads. `assemble_genome` with an
# `output_dir` persists the Stage-1 corrector handoff at
# `<output_dir>/corrected.fastq`, which is what we align.
#
# Identity is local (Smith-Waterman) over both orientations with
# AffineGapScoreModel(EDNAFULL), denominator matches+mismatches+insertions+
# deletions -- NOT a degenerate match=0 model, which collapses SW to the empty
# alignment and reports impossible identities.

import Pkg
const PROJECT = get(ENV, "MYCELIA_PROJECT", "")
Pkg.activate(PROJECT; io = devnull)

import Mycelia
import FASTX
import BioSequences
import BioAlignments
import Random
import Statistics

const LABEL = get(ENV, "PROBE_LABEL", "unlabeled")
const CELL_GENOME = parse(Int, get(ENV, "PROBE_GENOME", "6000"))
const CELL_READLEN = parse(Int, get(ENV, "PROBE_READLEN", "1200"))
const CELL_COV = parse(Int, get(ENV, "PROBE_COV", "20"))
const CELL_ERR = parse(Float64, get(ENV, "PROBE_ERR", "0.01"))
const SEEDS = [parse(Int, s) for s in split(get(ENV, "PROBE_SEEDS", "42,43,44"), ",")]

function simulate(seed::Int)
    reference_record =
        Mycelia.random_fasta_record(moltype = :DNA, seed = seed, L = CELL_GENOME)
    reference = FASTX.sequence(BioSequences.LongDNA{4}, reference_record)
    rng = Random.MersenneTwister(seed)
    Random.seed!(seed)
    n_reads = ceil(Int, CELL_COV * CELL_GENOME / CELL_READLEN)
    reads = FASTX.FASTQ.Record[]
    for i = 1:n_reads
        start = rand(rng, 1:(CELL_GENOME-CELL_READLEN+1))
        fragment = reference[start:(start+CELL_READLEN-1)]
        rand(rng, Bool) && (fragment = BioSequences.reverse_complement(fragment))
        observed, quals = Mycelia.observe(fragment; error_rate = CELL_ERR, tech = :nanopore)
        isempty(observed) && continue
        push!(
            reads,
            FASTX.FASTQ.Record(
                "read_$i",
                string(observed),
                String([Char(q + 33) for q in quals]),
            ),
        )
    end
    return reads, reference
end

const SCORE = BioAlignments.AffineGapScoreModel(
    BioAlignments.EDNAFULL;
    gap_open = -10,
    gap_extend = -1,
)

"""Local identity of `query` against `reference`, best of both orientations."""
function local_identity(
    query::BioSequences.LongDNA{4},
    reference::BioSequences.LongDNA{4},
)::Float64
    best = 0.0
    for candidate in (query, BioSequences.reverse_complement(query))
        result = BioAlignments.pairalign(
            BioAlignments.LocalAlignment(),
            candidate,
            reference,
            SCORE,
        )
        aln = BioAlignments.alignment(result)
        matches = BioAlignments.count_matches(aln)
        denom =
            matches +
            BioAlignments.count_mismatches(aln) +
            BioAlignments.count_insertions(aln) +
            BioAlignments.count_deletions(aln)
        denom == 0 && continue
        best = max(best, matches / denom)
    end
    return best
end

function corrected_identity(reads, reference, seed::Int)
    outdir = mktempdir()
    Random.seed!(1_042)
    Mycelia.Rhizomorph.assemble_genome(
        deepcopy(reads);
        k = 31,
        corrector = :iterative,
        strategy = :scalable,
        sequencing_tech = :nanopore,
        output_dir = outdir,
    )
    path = joinpath(outdir, "corrected.fastq")
    isfile(path) || return (n = 0, mean_identity = NaN)
    identities = Float64[]
    reader = FASTX.FASTQ.Reader(open(path))
    for record in reader
        seq = FASTX.sequence(BioSequences.LongDNA{4}, record)
        length(seq) < 50 && continue
        push!(identities, local_identity(seq, reference))
    end
    close(reader)
    rm(outdir; recursive = true, force = true)
    return (
        n = length(identities),
        mean_identity = isempty(identities) ? NaN : Statistics.mean(identities),
    )
end

println("label,seed,n_raw,raw_mean_identity,n_corrected,corrected_mean_identity")
for seed in SEEDS
    reads, reference = simulate(seed)
    raw = [
        local_identity(FASTX.sequence(BioSequences.LongDNA{4}, r), reference) for r in reads
    ]
    corrected = corrected_identity(reads, reference, seed)
    println(
        LABEL,
        ",",
        seed,
        ",",
        length(raw),
        ",",
        round(Statistics.mean(raw), digits = 5),
        ",",
        corrected.n,
        ",",
        round(corrected.mean_identity, digits = 5),
    )
    flush(stdout)
end
