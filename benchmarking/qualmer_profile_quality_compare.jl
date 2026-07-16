# Measure the corrector QUALITY gap between qualmer_memory_profile=:full and an
# aggregate profile (:ultralight_quality / :lightweight_quality) at real genome
# scale (td-n8ax). Canonical-kmer genome fraction is computed in Julia — no QUAST
# needed; the only external dep is ART for read simulation (Mycelia bootstraps it
# via Conda.jl). Answers "does aggregate storage preserve correction quality?" —
# the decision-2 authoritative gate, complementing the E. coli MEMORY gate.
#
# The 140bp unit-test toy is BELOW the corrector's regime and shows a spurious
# quality gap; this runs on real small genomes (phix 5.4kb, lambda 48kb) where the
# corrector operates as designed. Empirically (2026-07-16, phix @30x): :full and
# :ultralight_quality both reach genome fraction 0.987 (identical single 5316bp
# contig) — aggregate storage is quality-preserving at real scale.
#
# Env knobs:
#   MYCELIA_QC_ACCESSION  default "NC_001422.1" (phix); lambda = "J02459.1"
#   MYCELIA_QC_COVERAGE   default "30"
#   MYCELIA_QC_K          default "19"
#   MYCELIA_QC_SEED       default "42"
#   MYCELIA_QC_PROFILES   default "full,ultralight_quality" (comma-separated)
#
# Usage (bare metal; ART bootstraps on first use):
#   julia --project=. benchmarking/qualmer_profile_quality_compare.jl
#   MYCELIA_QC_ACCESSION=J02459.1 julia --project=. benchmarking/qualmer_profile_quality_compare.jl

import Mycelia, FASTX, BioSequences, Dates

const ACCESSION = get(ENV, "MYCELIA_QC_ACCESSION", "NC_001422.1")
const COV = parse(Int, get(ENV, "MYCELIA_QC_COVERAGE", "30"))
const K = parse(Int, get(ENV, "MYCELIA_QC_K", "19"))
const SEED = parse(Int, get(ENV, "MYCELIA_QC_SEED", "42"))
const PROFILES = Symbol.(split(get(ENV, "MYCELIA_QC_PROFILES", "full,ultralight_quality"), ","))

println("=== qualmer profile quality compare (td-n8ax) START (UTC): $(Dates.now(Dates.UTC)) ===")
println("accession=$ACCESSION coverage=$(COV)x k=$K seed=$SEED profiles=$(PROFILES)")

tmp = mktempdir(prefix = "qualmer_qc_")
refdir = joinpath(tmp, "ref")
mkpath(refdir)
ref = Mycelia.download_genome_by_accession(accession = ACCESSION, outdir = refdir, compressed = false)
refstr = uppercase(String(FASTX.sequence(first(collect(Mycelia.open_fastx(ref))))))
println("reference: $(length(refstr)) bp")

sim = Mycelia.simulate_illumina_reads(fasta = ref, coverage = COV, rndSeed = SEED,
    read_length = 150, paired = true, quiet = true, outbase = joinpath(tmp, "reads_$(COV)x"))
reads = FASTX.FASTQ.Record[]
for p in filter(!isnothing, [sim.forward_reads, sim.reverse_reads])
    for r in Mycelia.open_fastx(p)
        push!(reads, r)
    end
end
println("simulated $(length(reads)) reads")

rc(s) = String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(String(s))))
function ckmers(seqs)
    s = Set{String}()
    for x in String.(seqs), i in 1:(length(x) - K + 1)

        sub = x[i:(i + K - 1)]
        occursin(r"^[ACGT]+$", sub) || continue
        push!(s, min(sub, rc(sub)))
    end
    return s
end
refk = ckmers([refstr])
function gf(cs)
    isempty(cs) ? 0.0 :
    round(length(intersect(refk, ckmers(cs))) / length(refk); digits = 4)
end

results = Tuple{Symbol, Int, Float64, Int, Float64}[]
for profile in PROFILES
    cfg = Mycelia.Rhizomorph.AssemblyConfig(; k = K, corrector = :iterative,
        qualmer_memory_profile = profile, qualmer_prefilter_min_count = 2, verbose = false)
    t0 = time()
    res = try
        Mycelia.Rhizomorph.assemble_genome(reads, cfg)
    catch e
        println("  profile=$profile ERRORED: $e")
        continue
    end
    el = round(time() - t0; digits = 1)
    lens = sort(length.(String.(res.contigs)); rev = true)
    g = gf(res.contigs)
    push!(results, (profile, length(res.contigs), g, isempty(lens) ? 0 : lens[1], el))
    println("RESULT profile=$profile : n_contigs=$(length(res.contigs)) GF=$g " *
            "longest=$(isempty(lens) ? 0 : lens[1]) top5=$(lens[1:min(5, end)]) elapsed=$(el)s")
end

# Verdict: aggregate profiles must not materially underperform :full at real scale.
if any(r -> r[1] == :full, results)
    full_gf = first(r[3] for r in results if r[1] == :full)
    for (prof, _, g, _, _) in results
        prof == :full && continue
        verdict = g >= full_gf - 0.01 ? "PARITY" : "DEGRADED"
        println("VERDICT $prof vs :full genome_fraction: $g vs $full_gf -> $verdict")
    end
end
println("=== qualmer profile quality compare DONE (UTC): $(Dates.now(Dates.UTC)) ===")
