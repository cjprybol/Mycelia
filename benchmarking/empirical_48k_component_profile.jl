# empirical_48k_component_profile.jl
#
# EMPIRICALLY localize the residual super-linear (alpha ~1.47) term in the
# :scalable corrector at large genome scale (bead td-firp).
#
# WHY THIS HARNESS EXISTS (what the two prior attempts got wrong)
# --------------------------------------------------------------
# #381 (benchmarking/largescale_component_profile.jl) localized the alpha~1.47
# term to CONTIG-EXTRACTION (find_linear_path's O(L^2) neighbor-membership walk)
# by SOURCE ANALYSIS on a LINEAR PROXY graph: its Part B timed the reassembly
# sub-components on an ERROR-FREE read set, arguing the corrected-read graph is
# "near-linear". The #383 fix that followed did NOT move the 48 kb wall time
# (2109 -> 2163 s). The reason: the REAL corrected graph at 48 kb is BRANCHY
# (~98 separate contigs, none genome-spanning), so find_linear_path walks ~98
# SHORT paths and never enters its quadratic regime. The linear proxy was the
# WORST case for find_linear_path but the WRONG case for the real corrector.
#
# LESSON (and this harness's discipline): profile the ACTUAL corrector run
# empirically. Do NOT reconstruct the reassembly graph from an error-free proxy;
# do NOT reason about O(.) from source on a graph the corrector never produces.
#
# WHAT THIS HARNESS MEASURES (fully empirical, no proxy)
# -----------------------------------------------------
# For each genome size, it runs the SINGLE real public entry point that performs
# the WHOLE corrector pipeline in one call --
#
#     Mycelia.Rhizomorph.assemble_genome(reads;
#         corrector = :iterative, strategy = :scalable, k = 21, verbose = false)
#
# -- which internally does: iterative correction (mycelia_iterative_assemble,
# :scalable knobs) THEN reassembly of the corrected reads WITH graph_cleanup ON
# (clean_corrector_graph!, defaults ON only on this corrector reassembly path --
# NOT on a standalone corrector=:none reassembly, which is why the prior Part A
# split missed the new cleanup cost). The run is executed under Julia's stdlib
# Profile sampler, and every sample is bucketed to the component whose entry
# function sits on its call stack. Per-component seconds = wall * sample-share,
# aggregated across ALL k-rungs, ALL iterations, AND the final reassembly. This
# attributes the real end-to-end wall time -- including the freshly-added
# clean_corrector_graph! / _dedup_reverse_complements passes -- to each
# component with NO isolated proxy anywhere.
#
# It then fits alpha per component (log-log regression of seconds vs genome
# length, plus the large-end pairwise alpha) and names the component whose alpha
# exceeds ~1.3 at the large end as the TRUE alpha~1.47 driver.
#
# Read-only w.r.t. src: no library function is modified or monkey-patched. Only
# the stdlib sampler observes the unmodified corrector.
#
# RUN (from the worktree):
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
#     benchmarking/empirical_48k_component_profile.jl
#
# Optional env:
#   MYCELIA_E48_SIZES     comma genome sizes bp   (default 5000,10000,20000)
#   MYCELIA_E48_COVERAGE  fold coverage           (default 20)
#   MYCELIA_E48_READLEN   read length bp          (default 150)
#   MYCELIA_E48_ERR       per-base error rate     (default 0.01)
#   MYCELIA_E48_K         corrector max_k         (default 21)
#   MYCELIA_E48_SEED      RNG seed                (default 42)
#   MYCELIA_E48_DELAY     profile sampler delay s (default 0.002)
#   MYCELIA_E48_OUT       results file (append)   (default benchmarking/results/empirical_48k_component_profile.txt)

import Mycelia
import FASTX
import BioSequences
import Random
import Profile
import Printf
import Dates
import Statistics

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
_parse_sizes(s) = Int[parse(Int, strip(x)) for x in split(s, ",") if !isempty(strip(x))]

const SIZES    = _parse_sizes(get(ENV, "MYCELIA_E48_SIZES", "5000,10000,20000"))
const COVERAGE = parse(Float64, get(ENV, "MYCELIA_E48_COVERAGE", "20"))
const READLEN  = parse(Int, get(ENV, "MYCELIA_E48_READLEN", "150"))
const ERR_RATE = parse(Float64, get(ENV, "MYCELIA_E48_ERR", "0.01"))
const MAX_K    = parse(Int, get(ENV, "MYCELIA_E48_K", "21"))
const SEED     = parse(Int, get(ENV, "MYCELIA_E48_SEED", "42"))
const DELAY    = parse(Float64, get(ENV, "MYCELIA_E48_DELAY", "0.002"))
const TECH     = :illumina
const OUT_PATH = get(ENV, "MYCELIA_E48_OUT",
    joinpath(@__DIR__, "results", "empirical_48k_component_profile.txt"))

# ---------------------------------------------------------------------------
# Read simulation -- mirrors the validated illumina cell used by the prior
# scaling profilers (random_fasta_record + Mycelia.observe).
# ---------------------------------------------------------------------------
function simulate_reads(refseq::BioSequences.LongDNA{4}, readlen::Int, coverage::Real,
        error_rate::Float64, tech::Symbol, rng::Random.AbstractRNG)
    glen = length(refseq)
    effective_readlen = min(readlen, glen)
    n_reads = max(1, ceil(Int, coverage * glen / effective_readlen))
    records = Vector{FASTX.FASTQ.Record}()
    for i in 1:n_reads
        start = glen == effective_readlen ? 1 : rand(rng, 1:(glen - effective_readlen + 1))
        frag = refseq[start:(start + effective_readlen - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        obs_seq, quals = Mycelia.observe(frag; error_rate = error_rate, tech = tech)
        isempty(obs_seq) && continue
        qstr = String([Char(q + 33) for q in quals])
        push!(records, FASTX.FASTQ.Record("read_$(i)", string(obs_seq), qstr))
    end
    return records
end

function build_fixture(glen::Int)
    rng = Random.MersenneTwister(SEED)
    ref_rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = glen)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, ref_rec)
    reads = simulate_reads(refseq, READLEN, COVERAGE, ERR_RATE, TECH, rng)
    return refseq, reads
end

# ---------------------------------------------------------------------------
# Component buckets. Ordered leaf->root FIRST-MATCH, so put phase-specific,
# deeper-in-the-stack entry functions FIRST. The reassembly-side cleanup/dedup/
# extract buckets use function names that appear ONLY on the reassembly path, so
# they never steal correction-phase samples even though both phases build graphs
# and touch bubbles.
#   * clean_corrector_graph! -> clip_error_tips! / collapse_error_bubbles!
#     are the reassembly-only defrag pass (td-969e / #384) -- the PRIME new
#     suspect; distinct from the correction gate's _hard_vertex_set/
#     detect_bubbles_next.
# ---------------------------------------------------------------------------
const BUCKETS = [
    # --- reassembly-only components (matched first; names are phase-unique) ---
    ("R1. clean_corrector_graph! (NEW defrag)",
        ["clean_corrector_graph!", "clip_error_tips!", "collapse_error_bubbles!"]),
    ("R2. contig extraction (find_contigs)",
        ["find_contigs_next", "find_linear_path", "find_eulerian_paths_next",
         "_select_linear_neighbor", "_dominant_neighbor"]),
    ("R3. RC-dedup (_dedup_reverse_complements)",
        ["_dedup_reverse_complements", "_canonical_string"]),
    ("R4. consensus/quality (path->fastq)",
        ["_qualmer_path_to_consensus_fastq", "path_to_sequence",
         "_reconstruct_oriented_kmer_path", "_generate_fastq_contigs_from_qualmer_graph"]),
    # --- correction-phase components ---
    ("C5. per-read Viterbi decode",
        ["improve_read_likelihood", "try_viterbi_path_improvement",
         "find_optimal_sequence_path", "calculate_read_likelihood",
         "calculate_sequence_likelihood", "build_correction_weighted_graph",
         "weighted_graph_from_rhizomorph", "viterbi_decode_next",
         "beam_pruned_viterbi_decode_next", "_beam_pruned_viterbi_decode_next",
         "_get_valid_transitions", "_prune_viterbi_beam",
         "accumulate_competing_paths", "_enumerate_competing_paths",
         "correct_observations", "adjust_quality_scores"]),
    ("C4. stage0 cheap_correct", ["_stage0_cheap_correct", "_stage0_correct_read"]),
    ("C3. solid_kmer classification",
        ["_solid_kmer_set", "classify_kmers", "MixtureModelClassifier"]),
    ("C2. hard_vertex_set (corr bubbles)", ["_hard_vertex_set", "detect_bubbles_next"]),
    # --- shared build + IO (matched last so phase-specific work wins) ---
    ("B1. build_qualmer_graph", ["build_qualmer_graph"]),
    ("IO. write/read fastq", ["write_fastq", "_write_reads_to_fastq"]),
]

function classify_sample(frame_names::Vector{String})
    for name in frame_names            # leaf -> root
        for (label, keys) in BUCKETS
            for kw in keys
                occursin(kw, name) && return label
            end
        end
    end
    return "0. other (corrector/reassembly/GC/overhead)"
end

function split_samples(data::Vector{UInt64})
    samples = Vector{Vector{UInt64}}()
    cur = UInt64[]
    for ip in data
        if ip == 0
            isempty(cur) || (push!(samples, copy(cur)); empty!(cur))
        else
            push!(cur, ip)
        end
    end
    isempty(cur) || push!(samples, cur)
    return samples
end

function bucket_profile()
    data = Profile.fetch(include_meta = false)
    lidict = Profile.getdict(data)
    samples = split_samples(data)
    counts = Dict{String, Int}()
    for s in samples
        names = String[]
        for ip in s
            for sf in get(lidict, ip, Base.StackTraces.StackFrame[])
                push!(names, string(sf.func))
            end
        end
        label = classify_sample(names)
        counts[label] = get(counts, label, 0) + 1
    end
    total = sum(values(counts); init = 0)
    return counts, total
end

# ---------------------------------------------------------------------------
# One real end-to-end profiled run at a genome size.
# ---------------------------------------------------------------------------
struct SizeResult
    glen::Int
    nreads::Int
    wall::Float64
    total_samples::Int
    ncontigs::Int
    max_contig::Int
    k_progression::Vector{Int}
    n_passes::Int            # total correction iterations across all k-rungs
    decode_fracs::Vector{Float64}
    verts_before::Int        # graph vertices entering clean_corrector_graph!
    verts_after::Int
    tips_removed::Int
    bubbles_collapsed::Int
    seconds::Dict{String, Float64}   # component label -> attributed seconds
end

function run_size(io, glen::Int)
    _, reads = build_fixture(glen)
    println(io, "\n[size $(glen)bp] simulated $(length(reads)) reads " *
                "(readlen=$(READLEN), $(COVERAGE)x, err=$(ERR_RATE))")
    flush(io)

    Profile.clear()
    Profile.init(n = 10^8, delay = DELAY)
    local result
    wall = @elapsed begin
        redirect_stdout(devnull) do
            Profile.@profile begin
                result = Mycelia.Rhizomorph.assemble_genome(reads;
                    corrector = :iterative, strategy = :scalable,
                    k = MAX_K, verbose = false)
            end
        end
    end

    counts, total = bucket_profile()
    seconds = Dict{String, Float64}()
    for (label, n) in counts
        seconds[label] = total == 0 ? 0.0 : wall * n / total
    end

    # AssemblyResult introspection. contigs::Vector{String} (length = contig bp);
    # assembly_stats::Dict{String,Any} (STRING keys). Best-effort but with the
    # exact keys verified against the AssemblyResult shape.
    ncontigs = 0
    max_contig = 0
    kprog = Int[]
    n_passes = 0
    decode_fracs = Float64[]
    verts_before = 0
    verts_after = 0
    tips_removed = 0
    bubbles_collapsed = 0
    try
        contigs = result.contigs            # Vector{String}
        ncontigs = length(contigs)
        ncontigs > 0 && (max_contig = maximum(length(c) for c in contigs))
    catch
    end
    try
        stats = result.assembly_stats
        _iv(key) = (v = get(stats, key, nothing); v === nothing ? 0 : Int(v))
        haskey(stats, "k_progression") && (kprog = collect(Int, stats["k_progression"]))
        if haskey(stats, "decode_fraction_per_pass")
            decode_fracs = collect(Float64, stats["decode_fraction_per_pass"])
            n_passes = length(decode_fracs)
        end
        verts_before = _iv("graph_cleanup_vertices_before")
        verts_after = _iv("graph_cleanup_vertices_after")
        tips_removed = _iv("graph_cleanup_tips_removed")
        bubbles_collapsed = _iv("graph_cleanup_bubbles_collapsed")
    catch
    end

    println(io, Printf.@sprintf(
        "[size %dbp] wall=%.1fs samples=%d contigs=%d max_contig=%dbp k_prog=%s passes=%d cleanup=%d->%dV(-%dtips,-%dbub)",
        glen, wall, total, ncontigs, max_contig, isempty(kprog) ? "?" : string(kprog),
        n_passes, verts_before, verts_after, tips_removed, bubbles_collapsed))
    flush(io)

    return SizeResult(glen, length(reads), wall, total, ncontigs, max_contig, kprog,
        n_passes, decode_fracs, verts_before, verts_after, tips_removed, bubbles_collapsed,
        seconds)
end

# ---------------------------------------------------------------------------
# alpha fit: slope of log(seconds) vs log(genome_len).
# ---------------------------------------------------------------------------
function loglog_alpha(xs::Vector{Float64}, ys::Vector{Float64})
    # keep only strictly-positive y (component actually observed)
    idx = findall(y -> y > 0, ys)
    length(idx) < 2 && return NaN
    lx = log.(xs[idx]); ly = log.(ys[idx])
    mx = Statistics.mean(lx); my = Statistics.mean(ly)
    denom = sum((lx .- mx) .^ 2)
    denom == 0 && return NaN
    return sum((lx .- mx) .* (ly .- my)) / denom
end

pairwise_alpha(g1, t1, g2, t2) =
    (t1 > 0 && t2 > 0) ? log(t2 / t1) / log(g2 / g1) : NaN

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    mkpath(dirname(OUT_PATH))
    io = open(OUT_PATH, "a")
    header(sink) = begin
        println(sink, "="^92)
        println(sink, "EMPIRICAL :scalable per-component scaling on the REAL corrector run (td-firp)")
        println(sink, "  started : $(Dates.now())")
        println(sink, "  entry   : assemble_genome(reads; corrector=:iterative, strategy=:scalable, k=$(MAX_K))")
        println(sink, "  sizes   : $(SIZES) bp   readlen=$(READLEN)  cov=$(COVERAGE)x  err=$(ERR_RATE)")
        println(sink, "  sampler : Profile delay=$(DELAY)s   (per-component sec = wall * sample-share)")
        println(sink, "  NOTE    : full end-to-end real run (correction + reassembly WITH clean_corrector_graph!).")
        println(sink, "            NOT a proxy graph, NOT source O(.) reasoning, NOT an isolated component.")
        println(sink, "="^92)
    end
    header(stdout); header(io)

    # --- Warmup: compile the WHOLE assemble_genome(corrector=:iterative) path on a
    # tiny input so the profiled runs are ~pure runtime, not first-call compilation.
    # Warm up at the SAME k=MAX_K and a non-trivial genome so EVERY k-rung on the
    # ladder (3..21) plus the reassembly / clean_corrector_graph! / find_contigs /
    # dedup paths are compiled before the first measured run -- otherwise the first
    # real size absorbs k=21-rung + reassembly compilation into "0. other".
    println(stdout, "\n[warmup] compiling the full corrector+reassembly path (k=$(MAX_K)) ...")
    flush(stdout)
    _, wreads = build_fixture(1500)
    redirect_stdout(devnull) do
        Mycelia.Rhizomorph.assemble_genome(wreads;
            corrector = :iterative, strategy = :scalable,
            k = MAX_K, verbose = false)
    end
    println(stdout, "[warmup] done")
    flush(stdout)

    results = SizeResult[]
    for glen in SIZES
        push!(results, run_size(stdout, glen))
        # also log to file as we go so a long run is crash-durable
        rl = results[end]
        println(io, Printf.@sprintf(
            "[size %dbp] wall=%.1fs samples=%d contigs=%d max_contig=%dbp k_prog=%s passes=%d cleanup=%d->%dV(-%dtips,-%dbub)",
            rl.glen, rl.wall, rl.total_samples, rl.ncontigs, rl.max_contig,
            isempty(rl.k_progression) ? "?" : string(rl.k_progression),
            rl.n_passes, rl.verts_before, rl.verts_after, rl.tips_removed, rl.bubbles_collapsed))
        flush(io)
    end

    # --- Per-component scaling table + alpha fit ---
    all_labels = sort(collect(Set(Iterators.flatten(keys(r.seconds) for r in results))))
    gfloat = Float64[Float64(r.glen) for r in results]

    function emit_table(sink)
        println(sink, "\n" * "="^92)
        println(sink, "REAL-RUN GRAPH DIAGNOSTICS (evidence the corrected graph is branchy, not linear)")
        println(sink, Printf.@sprintf("  %-8s %8s %8s %10s %8s %14s %10s %10s",
            "size", "nreads", "contigs", "maxcontig", "passes", "k_prog", "cleanupV", "tips/bub"))
        for r in results
            println(sink, Printf.@sprintf("  %-8s %8d %8d %9dbp %8d %14s %5d->%-5d %5d/%-5d",
                "$(r.glen)bp", r.nreads, r.ncontigs, r.max_contig, r.n_passes,
                isempty(r.k_progression) ? "?" : string(r.k_progression),
                r.verts_before, r.verts_after, r.tips_removed, r.bubbles_collapsed))
        end
        println(sink, "  (max_contig << genome and contigs>>1 => branchy; find_linear_path never")
        println(sink, "   enters its quadratic regime -- confirming the #381 linear-proxy premise was wrong.)")

        println(sink, "\n" * "="^92)
        println(sink, "PER-COMPONENT SCALING TABLE (seconds attributed from the REAL run)")
        hdr = Printf.@sprintf("%-42s", "component")
        for r in results
            hdr *= Printf.@sprintf(" %10s", "$(r.glen)bp")
        end
        hdr *= Printf.@sprintf(" %8s %10s", "alpha", "alpha_lg")
        println(sink, hdr)
        println(sink, "-"^length(hdr))
        # sort components by seconds at the largest size (descending)
        largest = results[end]
        order = sort(all_labels; by = l -> -get(largest.seconds, l, 0.0))
        for label in order
            ys = Float64[get(r.seconds, label, 0.0) for r in results]
            a = loglog_alpha(gfloat, ys)
            # large-end pairwise alpha: last two sizes
            alg = length(results) >= 2 ?
                pairwise_alpha(gfloat[end-1], ys[end-1], gfloat[end], ys[end]) : NaN
            row = Printf.@sprintf("%-42s", label)
            for y in ys
                row *= Printf.@sprintf(" %10.1f", y)
            end
            row *= Printf.@sprintf(" %8s %10s",
                isnan(a) ? "-" : Printf.@sprintf("%.2f", a),
                isnan(alg) ? "-" : Printf.@sprintf("%.2f", alg))
            println(sink, row)
        end
        println(sink, "-"^length(hdr))
        totrow = Printf.@sprintf("%-42s", "TOTAL (wall)")
        walls = Float64[r.wall for r in results]
        for w in walls
            totrow *= Printf.@sprintf(" %10.1f", w)
        end
        atot = loglog_alpha(gfloat, walls)
        atotlg = length(results) >= 2 ?
            pairwise_alpha(gfloat[end-1], walls[end-1], gfloat[end], walls[end]) : NaN
        totrow *= Printf.@sprintf(" %8s %10s",
            isnan(atot) ? "-" : Printf.@sprintf("%.2f", atot),
            isnan(atotlg) ? "-" : Printf.@sprintf("%.2f", atotlg))
        println(sink, totrow)
        println(sink, "\n  alpha    = log-log regression slope of seconds vs genome length (full range)")
        println(sink, "  alpha_lg = large-end pairwise slope over the last two sizes")
        println(sink, "  Interpretation: the component with alpha_lg > ~1.3 AND a large second-count")
        println(sink, "  at the biggest size is the TRUE residual super-linear (alpha~1.47) driver.")

        # --- auto-verdict: dominant super-linear component at the large end ---
        best_label = ""
        best_score = -Inf
        for label in all_labels
            ys = Float64[get(r.seconds, label, 0.0) for r in results]
            alg = length(results) >= 2 ?
                pairwise_alpha(gfloat[end-1], ys[end-1], gfloat[end], ys[end]) : NaN
            (isnan(alg) || alg <= 1.3) && continue
            # score = large-end seconds weighted by how super-linear it is
            score = ys[end] * (alg - 1.0)
            if score > best_score
                best_score = score
                best_label = label
            end
        end
        println(sink, "\n" * "="^92)
        if isempty(best_label)
            println(sink, "AUTO-VERDICT: no single component exceeds alpha_lg 1.3 at the large end;")
            println(sink, "  the residual is diffuse across components (re-examine the table above).")
        else
            ys = Float64[get(r.seconds, label, 0.0) for label in [best_label] for r in results]
            println(sink, "AUTO-VERDICT: dominant super-linear component at the large end =")
            println(sink, "  >>> $(best_label) <<<")
            println(sink, "  ($(Printf.@sprintf("%.1f", ys[end]))s at $(results[end].glen)bp, " *
                          "$(Printf.@sprintf("%.0f%%", 100*ys[end]/results[end].wall)) of wall).")
            println(sink, "  This is the EMPIRICALLY MEASURED driver on the real corrector run.")
        end
        println(sink, "  finished: $(Dates.now())")
        println(sink, "="^92)
    end

    emit_table(stdout)
    emit_table(io)
    close(io)
    println(stdout, "\n[results appended to] $(OUT_PATH)")
end

main()
