# Diagnostic harness for the td-jt7r QUESTION: at the 2 kb / 1.2 kb-read / 8x /
# 5%-error cell, WHY does `corrector = :iterative` report a smaller largest contig
# (131 bp) than `corrector = :none` (413 bp), when at larger genomes / higher
# coverage correction helps?
#
# The harness was built around the working hypothesis that this is a CORRECTION
# REGRESSION — that the corrector damages the reads. That hypothesis is REFUTED by
# this harness's own output. Verdict as of this commit:
#
#   * Correction is NEUTRAL at the CEILING k. Re-assembling the RAW reads and the
#     CORRECTED reads under ONE fixed, production-matching config at k = 31 gives
#     the SAME largest contig, 492 bp (217 vs 219 contigs) — no correction cost at
#     all (`--k-sweep`; see `reassembly_k_sweep.csv`).
#   * The 413 -> 131 bp drop between the two arms is a k DROP-DOWN, not corrupted
#     reads. `select_reassembly_k` picks k = 11 (ceiling 31) for the corrected
#     reads, and under that same fixed config BOTH read sets collapse at k = 11:
#     raw 128 bp / 610 contigs, corrected 131 bp / 526 contigs. The corrected k = 11
#     sweep row REPRODUCES the corrected arm's production result exactly (526
#     contigs, 131 bp), which is the check that the sweep config now matches
#     production.
#   * The 413 bp naive baseline is itself a CONFIG ARTIFACT: the same raw reads at
#     the same k = 31 assemble to 492 bp once `dedup_revcomp` and `graph_cleanup`
#     are set the way the corrector's re-assembly sets them. The naive arm runs the
#     `corrector = :none` defaults, where both resolve FALSE, so the two arms differ
#     in three config values, not one (see `run_arm`).
#   * Provenance note: an earlier revision of this sweep left `dedup_revcomp` and
#     `graph_cleanup` at the `corrector = :none` defaults and reported 413 bp (raw)
#     vs 343 bp (corrected) at k = 31, and 166 bp for both at k = 11. Every one of
#     those numbers was carrying the config confound described above; the figures in
#     this block are the production-matched replacements.
#
# The script is a measurement harness, not a fix. For one (genome, read length,
# coverage, error rate) cell it reports:
#
#   * assembly quality (n_contigs, largest contig, largest/genome ratio) — BOTH arms
#   * the coverage-aware re-assembly k actually chosen (`reassembly_k`) and the
#     ceiling it was chosen from — CORRECTED ARM ONLY. `select_reassembly_k` lives on
#     the corrector path, so `run_arm` returns `missing` for these on the `:none` arm
#     and the summary columns are populated from the corrected arm alone. The
#     `select_k_on_raw` column is a COUNTERFACTUAL computed here by `k_ladder` (what
#     the chooser WOULD return on raw reads); production never evaluates it there.
#   * the `median_solid_kmer_multiplicity` ladder that drives the choice, on RAW and
#     on CORRECTED reads — both arms, but on a fixed DIAGNOSTIC ladder that is not
#     production's candidate set (see `k_ladder`)
#   * DIRECT over-correction evidence: per-read identity of RAW reads and of
#     CORRECTED reads against the known reference — LOCAL (Smith-Waterman)
#     BioAlignments pairwise alignment, scored over the ALIGNED BLOCK ONLY, best of
#     the two orientations. Semi-global alignment is deliberately REJECTED; see
#     `reference_identity` for the arithmetic. Plus the fraction of each read's
#     31-mers that occur in the reference.
#   * a solid-k-mer COMPOSITION audit (`solid_composition`) splitting the solid set
#     into genomic and error k-mers — this is what refutes
#     `median_solid_kmer_multiplicity`'s backbone-coverage claim.
#
# Corrected reads are obtained by giving `assemble_genome` an `output_dir`, which
# persists the Stage-1 corrector handoff at `<output_dir>/corrected.fastq`
# instead of deleting it (see `_run_stage1_correction`).
#
# Usage:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/indel_toy_overcorrection_diagnostic.jl \
#     --genome=2000 --read-length=1200 --coverage=8 --error=0.05 \
#     --seed=42 --out=/tmp/diag-cell --label=baseline \
#     [--k-sweep=11,31] [--skip-alignment=true]
#
# `--k-sweep=<comma-separated ks>` is OFF by default and runs the decisive
# raw-vs-corrected fixed-config re-assembly sweep (the expensive part).
# `--skip-alignment=true` skips the per-read pairwise-alignment identity statistics.
#
# Emits, under `<out>/`:
#   * `cell_summary.csv`        one row per cell
#   * `k_ladder.csv`            median-solid-k-mer ladder, raw + corrected arms
#   * `solid_composition.csv`   genomic-vs-error solid-k-mer audit
#   * `reassembly_k_sweep.csv`  ONLY when `--k-sweep` is given

import BioAlignments
import BioSequences
import CSV
import DataFrames
import FASTX
import Mycelia
import Random
import Statistics

const DIAG_MAX_K = 31
const DIAG_IDENTITY_KMER = 31

"""
Parse `--key=value` CLI arguments into a `Dict{String, String}`.
"""
function parse_args(args::Vector{String})::Dict{String, String}
    parsed = Dict{String, String}()
    for arg in args
        startswith(arg, "--") || continue
        body = arg[3:end]
        pieces = split(body, '='; limit = 2)
        parsed[String(pieces[1])] = length(pieces) == 2 ? String(pieces[2]) : "true"
    end
    return parsed
end

"""
Build the read fixture the same way `benchmarking/indel_frontier_fixed_toy_proof.jl`
does, but with genome length / read length / coverage / error rate as parameters.
Reads are simulated through `Mycelia.observe(...; tech = :nanopore)` so they carry
INDELS as well as substitutions.

NOT byte-for-byte the sibling's procedure: the sibling hardcodes its cell and then
ASSERTS a minimum observed read length (`INDEL_TOY_MIN_OBSERVED_READ_LENGTH`) before
assembling. This harness has no such guard, because it is swept across cells where
short observed reads are a legitimate outcome to measure rather than a fixture bug.
Reads that `observe` returns empty are dropped, so `n_reads` requested is an upper
bound on the reads actually emitted.

Note the pairing this creates with `run_arm`: the reads carry INDELS (60% of the
simulated errors at `tech = :nanopore`), but both of `run_arm`'s `assemble_genome`
calls hardcode `sequencing_tech = :illumina`, which selects the SUBSTITUTION-ONLY
corrector (Mycelia logs a warning to this effect on every run). That pairing is
deliberate — it reproduces the production configuration under investigation — but it
means "correction" in this harness never repairs an indel, only a substitution.
"""
function make_fixture(
        genome_length::Int,
        read_length::Int,
        coverage::Int,
        error_rate::Float64,
        seed::Int
)::Tuple{Vector{FASTX.FASTQ.Record}, BioSequences.LongDNA{4}}
    reference_record = Mycelia.random_fasta_record(
        moltype = :DNA, seed = seed, L = genome_length
    )
    reference = FASTX.sequence(BioSequences.LongDNA{4}, reference_record)
    rng = Random.MersenneTwister(seed)
    Random.seed!(seed)
    n_reads = ceil(Int, coverage * genome_length / read_length)
    reads = FASTX.FASTQ.Record[]
    for read_index in 1:n_reads
        start_position = rand(rng, 1:(genome_length - read_length + 1))
        fragment = reference[start_position:(start_position + read_length - 1)]
        if rand(rng, Bool)
            fragment = BioSequences.reverse_complement(fragment)
        end
        observed,
        qualities = Mycelia.observe(
            fragment; error_rate = error_rate, tech = :nanopore
        )
        isempty(observed) && continue
        quality_string = String([Char(quality + 33) for quality in qualities])
        push!(
            reads,
            FASTX.FASTQ.Record(
                "read_$(read_index)", string(observed), quality_string
            )
        )
    end
    return reads, reference
end

"""
Local (Smith-Waterman) identity of `query` against `reference`, taking the better
of the two orientations. Identity is `matches / (matches + mismatches +
insertions + deletions)` INSIDE the aligned block, so indels are penalized while
the unaligned reference flanks are not. `SemiGlobalAlignment` is deliberately NOT
used: BioAlignments still emits the free reference end-gaps as deletion
operations, so a 1.2 kb read against a 2 kb reference scores ~0.57 identity purely
from the 800 bp of untouched flank.

Returns `(identity, aligned_fraction)`, where `aligned_fraction` is the share of the
query consumed by the local block. That fraction is REPORTED, not ENFORCED: the
better orientation is selected on identity alone, and neither this function nor
`read_set_stats` applies any coverage gate or threshold. A truncated alignment can
therefore still post a high identity — it is DETECTABLE by reading
`mean_aligned_fraction` alongside `mean_identity`, but nothing rejects or flags it.

`score_model` MUST be MATCH-POSITIVE (positive match scores, negative mismatch and
gap penalties). `LocalAlignment` MAXIMIZES score, so a DISTANCE model — matches
costing 0 and every edit costing more — makes the empty alignment optimal and this
function degenerates SILENTLY to the `(0.0, 0.0)` sentinel below on every read.
`main` supplies `AffineGapScoreModel(EDNAFULL; gap_open = -10, gap_extend = -1)`:
EDNAFULL is the standard match-positive nucleotide substitution matrix (+5 match /
-4 mismatch), and `-10 / -1` is the conventional BLASTN-like affine setting, chosen
so that a single indel — the dominant error class in this fixture — costs less than
abandoning the block and does not split an otherwise-good local alignment in two.

SENTINEL: when BOTH orientations yield an empty alignment (denominator 0, so both
`continue`), the function returns `(0.0, 0.0)`. That is INDISTINGUISHABLE from a
genuine zero-identity / zero-coverage alignment, and it enters the means in
`read_set_stats` as a real 0.0, depressing `mean_identity` and
`mean_aligned_fraction`. It should be rare for kilobase reads against a kilobase
reference, but it is not detected or counted, so read those means accordingly.
"""
function reference_identity(
        query::BioSequences.LongDNA{4},
        reference::BioSequences.LongDNA{4},
        score_model
)::Tuple{Float64, Float64}
    best_identity = 0.0
    best_coverage = 0.0
    query_length = length(query)
    for candidate in (query, BioSequences.reverse_complement(query))
        result = BioAlignments.pairalign(
            BioAlignments.LocalAlignment(), candidate, reference, score_model
        )
        alignment = BioAlignments.alignment(result)
        matches = BioAlignments.count_matches(alignment)
        mismatches = BioAlignments.count_mismatches(alignment)
        insertions = BioAlignments.count_insertions(alignment)
        deletions = BioAlignments.count_deletions(alignment)
        denominator = matches + mismatches + insertions + deletions
        denominator == 0 && continue
        identity = matches / denominator
        if identity > best_identity
            best_identity = identity
            best_coverage = query_length == 0 ? 0.0 :
                            (matches + mismatches + insertions) / query_length
        end
    end
    return (best_identity, best_coverage)
end

"""
The set of `sequence`'s canonical k-mers at size `k`. Each k-mer is represented by
the `String` that is lexicographically smaller of the window and its reverse
complement — a literal string, NOT a numeric hash, so representatives can be
compared directly across read sets. Used for the reference-k-mer-recovery statistic.

Returns an EMPTY set when `length(sequence) < k` (no window exists). Callers cannot
distinguish that from a genuinely k-mer-free sequence; `reference_kmer_fraction`
against an empty reference set scores every read 0.0.
"""
function canonical_kmer_set(sequence::BioSequences.LongDNA{4}, k::Int)::Set{String}
    kmers = Set{String}()
    length(sequence) < k && return kmers
    for start_index in 1:(length(sequence) - k + 1)
        window = sequence[start_index:(start_index + k - 1)]
        forward = string(window)
        reverse = string(BioSequences.reverse_complement(window))
        push!(kmers, min(forward, reverse))
    end
    return kmers
end

"""
Fraction of `sequence`'s canonical k-mers that occur in `reference_kmers`. A read
that has been corrected TOWARD the reference raises this; a read corrected toward
a wrong consensus lowers it.

Returns `nothing` — NOT `0.0` — when `length(sequence) < k` (no window exists), and
`read_set_stats` DROPS those reads instead of scoring them. Consequence:
`mean_reference_kmer_fraction` is averaged over only the reads at least
`DIAG_IDENTITY_KMER` long, while `mean_identity` is averaged over ALL reads. The two
headline read statistics therefore have DIFFERENT denominators whenever any read is
shorter than `DIAG_IDENTITY_KMER`, and the harness does not report how many reads
were dropped. Compare them with that in mind.
"""
function reference_kmer_fraction(
        sequence::BioSequences.LongDNA{4},
        reference_kmers::Set{String},
        k::Int
)::Union{Float64, Nothing}
    length(sequence) < k && return nothing
    total = 0
    hits = 0
    for start_index in 1:(length(sequence) - k + 1)
        window = sequence[start_index:(start_index + k - 1)]
        forward = string(window)
        reverse = string(BioSequences.reverse_complement(window))
        total += 1
        hits += (min(forward, reverse) in reference_kmers) ? 1 : 0
    end
    total == 0 && return nothing
    return hits / total
end

"""
Per-read identity + reference-k-mer-recovery statistics for one read set.

Purely descriptive: every read is scored and every score is kept. There is NO
coverage gate — a read whose local alignment consumed only a fraction of the query
contributes its (possibly high) identity to `mean_identity` unchanged; the only
signal that this happened is `mean_aligned_fraction`. See `reference_identity`.

`mean_reference_kmer_fraction` is averaged over a SMALLER read set than the identity
statistics whenever any read is shorter than `DIAG_IDENTITY_KMER` — see
`reference_kmer_fraction`.
"""
function read_set_stats(
        sequences::Vector{BioSequences.LongDNA{4}},
        reference::BioSequences.LongDNA{4},
        reference_kmers::Set{String},
        score_model
)::NamedTuple
    identities = Float64[]
    coverages = Float64[]
    kmer_fractions = Float64[]
    for sequence in sequences
        identity, coverage = reference_identity(sequence, reference, score_model)
        push!(identities, identity)
        push!(coverages, coverage)
        fraction = reference_kmer_fraction(sequence, reference_kmers, DIAG_IDENTITY_KMER)
        fraction === nothing || push!(kmer_fractions, fraction)
    end
    return (;
        n = length(sequences),
        mean_identity = isempty(identities) ? NaN : Statistics.mean(identities),
        median_identity = isempty(identities) ? NaN : Statistics.median(identities),
        min_identity = isempty(identities) ? NaN : minimum(identities),
        mean_aligned_fraction = isempty(coverages) ? NaN : Statistics.mean(coverages),
        mean_length = isempty(sequences) ? NaN :
                      Statistics.mean(Float64.(length.(sequences))),
        mean_reference_kmer_fraction = isempty(kmer_fractions) ? NaN :
                                       Statistics.mean(kmer_fractions),
        identities = identities
    )
end

"""
Canonical k-mer occurrence counts over `sequences`, keyed by the canonical window
string.

Uses the SAME equivalence class as `_kmer_count_spectrum` (a window and its reverse
complement are one k-mer) but a DIFFERENT representative: this keys on the
lexicographically smaller window STRING, whereas `_kmer_count_spectrum` keys on the
smaller of the two window HASHES. The resulting count multisets agree; the KEYS are
not interchangeable and must not be joined across the two functions. The literal
string is kept here precisely so the solid set can be split against the reference.
"""
function canonical_counts(sequences, k::Int)::Dict{String, Int}
    counts = Dict{String, Int}()
    for sequence in sequences
        length(sequence) < k && continue
        for start_index in 1:(length(sequence) - k + 1)
            window = sequence[start_index:(start_index + k - 1)]
            key = min(
                string(window), string(BioSequences.reverse_complement(window))
            )
            counts[key] = get(counts, key, 0) + 1
        end
    end
    return counts
end

"""
Split the "solid" k-mer set that `median_solid_kmer_multiplicity` medians over into
GENOMIC k-mers (present in the reference) and ERROR k-mers (absent), and report the
median occurrence of each population separately.

Also reports `genomic_recall` = (number of SOLID genomic k-mers) / (number of
DISTINCT canonical k-mers in the REFERENCE at `k`) — the share of the true backbone
that survives the solidity cut. The denominator is the reference's own canonical
k-mer set, not the read set's, so `genomic_recall` is a backbone-completeness
measure and is `NaN` when the reference is shorter than `k`.

The solidity threshold is HARDCODED here at occurrence `>= 2`. That matches the
DEFAULT `solid_min` of `median_solid_kmer_multiplicity` and `select_reassembly_k`,
but production PARAMETERIZES it, so this audit tracks production only while callers
leave that default in place; it is not wired to the production value.

This is the audit that shows whether the statistic is measuring the genomic
backbone it claims to measure. When the error population outnumbers the genomic
one, the OVERALL median is an error-population statistic and no longer tracks
read coverage.
"""
function solid_composition(
        sequences,
        reference::BioSequences.LongDNA{4},
        k::Int
)::NamedTuple
    counts = canonical_counts(sequences, k)
    genomic = Set(keys(canonical_counts([reference], k)))
    solid = [(key, count) for (key, count) in counts if count >= 2]
    genomic_solid = [count for (key, count) in solid if key in genomic]
    error_solid = [count for (key, count) in solid if !(key in genomic)]
    safe_median(values) = isempty(values) ? NaN : Statistics.median(values)
    return (;
        k = k,
        n_solid = length(solid),
        n_solid_genomic = length(genomic_solid),
        n_solid_error = length(error_solid),
        median_all = safe_median([count for (_, count) in solid]),
        median_genomic = safe_median(genomic_solid),
        median_error = safe_median(error_solid),
        genomic_recall = isempty(genomic) ? NaN :
                         length(genomic_solid) / length(genomic)
    )
end

"""
Read a FASTQ file into `LongDNA{4}` sequences.
"""
function read_fastq_sequences(path::AbstractString)::Vector{BioSequences.LongDNA{4}}
    sequences = BioSequences.LongDNA{4}[]
    open(FASTX.FASTQ.Reader, String(path)) do reader
        for record in reader
            push!(sequences, FASTX.sequence(BioSequences.LongDNA{4}, record))
        end
    end
    return sequences
end

"""
`median_solid_kmer_multiplicity` at each k on a FIXED DIAGNOSTIC ladder, plus the k
`select_reassembly_k` would choose. This is the signal that decides the
corrected-read re-assembly k.

The ladder is NOT "every prime candidate k up to `ceiling_k`". It is hardcoded as
`[7, 11, 13, 17, 19, 23, 29, ceiling_k]`, then deduplicated and filtered to
`<= ceiling_k`. So it EXCLUDES the primes 2, 3 and 5, and it APPENDS `ceiling_k`
whether or not that value is prime.

More importantly, the ladder is a READOUT, not production's candidate set. `chosen`
calls `select_reassembly_k` with its default kwargs — including
`genome_size_floor = true` and `size_ref_k = 17` — and that genome-size uniqueness
floor (`_genome_size_floor_k`) can raise the EFFECTIVE floor well above 7. The
emitted `k_ladder.csv` will therefore routinely show low rungs (7, 11, 13, ...) that
the chooser can NEVER return for that read set. Any statement that the connectivity
floor is unmet "at every candidate k" is a statement about THIS measured ladder, not
about production's genome-size-floored candidate set.

Corollary for `cell_summary.csv`: `select_k_on_raw` is a COUNTERFACTUAL. Production
invokes `select_reassembly_k` only on the corrector path, so no production run ever
evaluates it on raw reads; the column name reads like an observation but records a
what-if computed here. `select_k_on_corrected` does correspond to the corrected
arm's actual `reassembly_k`.
"""
function k_ladder(reads, ceiling_k::Int)::NamedTuple
    candidate_ks = [7, 11, 13, 17, 19, 23, 29, ceiling_k]
    unique!(sort!(candidate_ks))
    filter!(<=(ceiling_k), candidate_ks)
    medians = [Mycelia.Rhizomorph.median_solid_kmer_multiplicity(reads, k)
               for k in candidate_ks]
    chosen = Mycelia.Rhizomorph.select_reassembly_k(reads, ceiling_k)
    return (; candidate_ks, medians, chosen)
end

"""
Assemble one arm and return quality + provenance fields.

Each arm calls production `assemble_genome` with the corrector setting under test
and NO other config overrides, so each arm inherits THAT corrector's production
defaults. This is faithful to production, but it means the two arms differ in THREE
config values, not one:

  * `corrector`     — `:none` vs `:iterative` (the variable under test)
  * `dedup_revcomp` — resolved to `corrector == :iterative`
                      (`src/rhizomorph/assembly.jl:271-272`), so FALSE on the naive
                      arm and TRUE on the corrected arm's re-assembly, which is
                      handed the outer resolved value at `assembly.jl:1654`
  * `graph_cleanup` — the plain (`:none`) route cleans only on an EXPLICIT `true`
                      (`assembly.jl:2163`), so FALSE here; the corrector's
                      re-assembly defaults it ON (`assembly.jl:1651`, passed at
                      `assembly.jl:1655`)

So the naive-vs-corrected largest-contig delta reported by this function is a
THREE-variable difference and cannot on its own attribute the change to correction.
The `--k-sweep` path is the arm that isolates it: it re-assembles both read sets
under ONE fixed config that matches the corrector's production re-assembly.

`reassembly_k`, `reassembly_k_ceiling`, `reassembly_graph_reused` and
`corrected_read_count` are emitted only on the corrector path, so the `:none` arm
returns `missing` for all four.
"""
function run_arm(
        reads::Vector{FASTX.FASTQ.Record},
        genome_length::Int,
        corrector::Symbol,
        output_dir::Union{String, Nothing}
)::NamedTuple
    Random.seed!(1_042)
    assembly = nothing
    wall_seconds = @elapsed begin
        if corrector == :none
            assembly = Mycelia.Rhizomorph.assemble_genome(
                deepcopy(reads);
                k = DIAG_MAX_K,
                corrector = :none,
                strategy = :scalable,
                sequencing_tech = :illumina
            )
        else
            assembly = Mycelia.Rhizomorph.assemble_genome(
                deepcopy(reads);
                k = DIAG_MAX_K,
                corrector = :iterative,
                strategy = :scalable,
                sequencing_tech = :illumina,
                output_dir = output_dir
            )
        end
    end
    contig_lengths = [length(contig) for contig in assembly.contigs]
    largest = isempty(contig_lengths) ? 0 : maximum(contig_lengths)
    stats = assembly.assembly_stats
    return (;
        corrector = String(corrector),
        n_contigs = length(assembly.contigs),
        largest = largest,
        ratio = largest / genome_length,
        wall_seconds = wall_seconds,
        reassembly_k = get(stats, "reassembly_k", missing),
        reassembly_k_ceiling = get(stats, "reassembly_k_ceiling", missing),
        reassembly_graph_reused = get(stats, "reassembly_graph_reused", missing),
        corrected_read_count = get(stats, "corrected_read_count", missing)
    )
end

function main(args::Vector{String} = ARGS)::Nothing
    parsed = parse_args(args)
    genome_length = parse(Int, get(parsed, "genome", "2000"))
    read_length = parse(Int, get(parsed, "read-length", "1200"))
    coverage = parse(Int, get(parsed, "coverage", "8"))
    error_rate = parse(Float64, get(parsed, "error", "0.05"))
    seed = parse(Int, get(parsed, "seed", "42"))
    label = get(parsed, "label", "cell")
    output_dir = abspath(get(parsed, "out", joinpath(tempdir(), "indel-toy-diag")))
    skip_alignment = get(parsed, "skip-alignment", "false") == "true"
    mkpath(output_dir)

    reads, reference = make_fixture(
        genome_length, read_length, coverage, error_rate, seed
    )
    println("[$(label)] genome=$(genome_length) read_length=$(read_length) " *
            "coverage=$(coverage) error=$(error_rate) reads=$(length(reads))")

    naive = run_arm(reads, genome_length, :none, nothing)
    println("[$(label)] naive: contigs=$(naive.n_contigs) largest=$(naive.largest) " *
            "ratio=$(round(naive.ratio; digits = 4))")

    corrected_dir = joinpath(output_dir, "corrected")
    corrected = run_arm(reads, genome_length, :iterative, corrected_dir)
    println("[$(label)] corrected: contigs=$(corrected.n_contigs) " *
            "largest=$(corrected.largest) ratio=$(round(corrected.ratio; digits = 4)) " *
            "reassembly_k=$(corrected.reassembly_k)/$(corrected.reassembly_k_ceiling)")

    corrected_fastq = joinpath(corrected_dir, "corrected.fastq")
    corrected_sequences = isfile(corrected_fastq) ?
                          read_fastq_sequences(corrected_fastq) :
                          BioSequences.LongDNA{4}[]
    raw_sequences = [FASTX.sequence(BioSequences.LongDNA{4}, record) for record in reads]

    raw_ladder = k_ladder(reads, DIAG_MAX_K)
    corrected_ladder = isempty(corrected_sequences) ? nothing :
                       k_ladder(corrected_sequences, DIAG_MAX_K)

    # Decisive test for the re-assembly-k hypothesis: assemble the SAME raw reads
    # and the SAME corrected reads with corrector=:none at a range of k, so the
    # corrector's contribution and the k choice are separated. If corrected reads
    # at the ceiling k beat raw reads at the ceiling k, correction is fine and the
    # coverage-aware k drop-down is what costs the assembly.
    #
    # CONFIG: `corrector = :none` here is the SWEEP MECHANISM (it pins k instead of
    # letting the corrector choose it), NOT the production corrector setting under
    # test. Left at its defaults, a `corrector = :none` call resolves BOTH
    # `dedup_revcomp` and `graph_cleanup` to FALSE (assembly.jl:271-272 and
    # assembly.jl:2163) — which would make the sweep matched arm-to-arm but NOT
    # matched to production, since the corrector's real re-assembly runs both TRUE
    # (assembly.jl:1651, threaded at assembly.jl:1654-1655). That is the exact
    # confound this sweep exists to remove, so both are passed EXPLICITLY TRUE
    # below. Both sweep arms share that one config, so the only variables across
    # rows are (read set, k).
    k_sweep_spec = get(parsed, "k-sweep", "")
    k_sweep = DataFrames.DataFrame(
        label = String[], arm = String[], k = Int[],
        n_contigs = Int[], largest = Int[], ratio = Float64[]
    )
    if !isempty(k_sweep_spec)
        sweep_ks = parse.(Int, split(k_sweep_spec, ','))
        corrected_records = isfile(corrected_fastq) ?
                            collect(open(FASTX.FASTQ.Reader, corrected_fastq) do reader
            collect(reader)
        end) : FASTX.FASTQ.Record[]
        for (arm, arm_reads) in (("raw", reads), ("corrected", corrected_records))
            isempty(arm_reads) && continue
            for k in sweep_ks
                Random.seed!(1_042)
                assembly = Mycelia.Rhizomorph.assemble_genome(
                    deepcopy(arm_reads);
                    k = k, corrector = :none, strategy = :scalable,
                    sequencing_tech = :illumina,
                    dedup_revcomp = true, graph_cleanup = true
                )
                lengths = [length(contig) for contig in assembly.contigs]
                largest = isempty(lengths) ? 0 : maximum(lengths)
                push!(k_sweep,
                    (label, arm, k, length(assembly.contigs), largest,
                        largest / genome_length))
                println("[$(label)] k-sweep $(arm) k=$(k): " *
                        "contigs=$(length(assembly.contigs)) largest=$(largest)")
            end
        end
        CSV.write(joinpath(output_dir, "reassembly_k_sweep.csv"), k_sweep)
    end

    raw_stats = nothing
    corrected_stats = nothing
    if !skip_alignment
        score_model = BioAlignments.AffineGapScoreModel(
            BioAlignments.EDNAFULL; gap_open = -10, gap_extend = -1
        )
        reference_kmers = canonical_kmer_set(reference, DIAG_IDENTITY_KMER)
        raw_stats = read_set_stats(
            raw_sequences, reference, reference_kmers, score_model
        )
        corrected_stats = isempty(corrected_sequences) ? nothing :
                          read_set_stats(
            corrected_sequences, reference, reference_kmers, score_model
        )
        println("[$(label)] raw reads:       n=$(raw_stats.n) " *
                "mean_identity=$(round(raw_stats.mean_identity; digits = 5)) " *
                "mean_ref_kmer_frac=$(round(raw_stats.mean_reference_kmer_fraction; digits = 5))")
        if corrected_stats !== nothing
            println("[$(label)] corrected reads: n=$(corrected_stats.n) " *
                    "mean_identity=$(round(corrected_stats.mean_identity; digits = 5)) " *
                    "mean_ref_kmer_frac=$(round(corrected_stats.mean_reference_kmer_fraction; digits = 5))")
        end
    end

    summary = DataFrames.DataFrame(
        label = [label],
        genome_length = [genome_length],
        read_length = [read_length],
        coverage = [coverage],
        error_rate = [error_rate],
        seed = [seed],
        n_reads = [length(reads)],
        naive_contigs = [naive.n_contigs],
        naive_largest = [naive.largest],
        naive_ratio = [naive.ratio],
        naive_seconds = [naive.wall_seconds],
        corrected_contigs = [corrected.n_contigs],
        corrected_largest = [corrected.largest],
        corrected_ratio = [corrected.ratio],
        corrected_seconds = [corrected.wall_seconds],
        corrected_read_count = [corrected.corrected_read_count],
        reassembly_k = [corrected.reassembly_k],
        reassembly_k_ceiling = [corrected.reassembly_k_ceiling],
        reassembly_graph_reused = [corrected.reassembly_graph_reused],
        select_k_on_raw = [raw_ladder.chosen],
        select_k_on_corrected = [
            corrected_ladder === nothing ? missing : corrected_ladder.chosen
        ],
        raw_mean_identity = [
            raw_stats === nothing ? missing : raw_stats.mean_identity
        ],
        corrected_mean_identity = [
            corrected_stats === nothing ? missing : corrected_stats.mean_identity
        ],
        raw_median_identity = [
            raw_stats === nothing ? missing : raw_stats.median_identity
        ],
        corrected_median_identity = [
            corrected_stats === nothing ? missing : corrected_stats.median_identity
        ],
        raw_mean_ref_kmer_fraction = [
            raw_stats === nothing ? missing : raw_stats.mean_reference_kmer_fraction
        ],
        corrected_mean_ref_kmer_fraction = [
            corrected_stats === nothing ? missing :
            corrected_stats.mean_reference_kmer_fraction
        ],
        raw_mean_aligned_fraction = [
            raw_stats === nothing ? missing : raw_stats.mean_aligned_fraction
        ],
        corrected_mean_aligned_fraction = [
            corrected_stats === nothing ? missing :
            corrected_stats.mean_aligned_fraction
        ],
        raw_mean_length = [raw_stats === nothing ? missing : raw_stats.mean_length],
        corrected_mean_length = [
            corrected_stats === nothing ? missing : corrected_stats.mean_length
        ]
    )
    summary_path = joinpath(output_dir, "cell_summary.csv")
    CSV.write(summary_path, summary)

    ladder = DataFrames.DataFrame(
        label = String[], arm = String[], k = Int[], median_solid = Float64[]
    )
    for (k, median_solid) in zip(raw_ladder.candidate_ks, raw_ladder.medians)
        push!(ladder, (label, "raw", k, median_solid))
    end
    if corrected_ladder !== nothing
        for (k, median_solid) in
            zip(corrected_ladder.candidate_ks, corrected_ladder.medians)
            push!(ladder, (label, "corrected", k, median_solid))
        end
    end
    ladder_path = joinpath(output_dir, "k_ladder.csv")
    CSV.write(ladder_path, ladder)

    # Audit the statistic itself: is the solid-k-mer median tracking the genomic
    # backbone, or the error population?
    composition = DataFrames.DataFrame(
        label = String[], arm = String[], k = Int[], n_solid = Int[],
        n_solid_genomic = Int[], n_solid_error = Int[], median_all = Float64[],
        median_genomic = Float64[], median_error = Float64[],
        genomic_recall = Float64[]
    )
    for (arm, arm_sequences) in
        (("raw", raw_sequences), ("corrected", corrected_sequences))
        isempty(arm_sequences) && continue
        for k in (7, 11, 17, DIAG_MAX_K)
            row = solid_composition(arm_sequences, reference, k)
            push!(composition,
                (label, arm, row.k, row.n_solid, row.n_solid_genomic,
                    row.n_solid_error, row.median_all, row.median_genomic,
                    row.median_error, row.genomic_recall))
        end
    end
    composition_path = joinpath(output_dir, "solid_composition.csv")
    CSV.write(composition_path, composition)
    println("[$(label)] wrote $(composition_path)")

    println("[$(label)] wrote $(summary_path)")
    println("[$(label)] wrote $(ladder_path)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
