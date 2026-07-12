# shortread_scalable_scaling_profile.jl
#
# Localize the SHORT-READ super-linear term in the :scalable corrector (td-9q84).
#
# WHY THIS EXISTS
# ---------------
# After two O(V^2) fixes (per-read weighted-graph rebuild hoisted out of the
# decode loop, #374; O(V+E) bubble detection, #373), a 5 kb :scalable run
# COMPLETES but the SWEEP end-to-end iterative time is still super-linear for
# SHORT (readlen=150, realistic Illumina) reads:
#     1kb 18s / 2kb 36s / 3kb 138s / 5kb 1346s   (alpha ~2.7, ACCELERATING).
# The long-read profiler (readlen=5000) blamed decode VOLUME and recommended
# windowed decode — but 150 bp reads have NO windows, so a DIFFERENT term
# dominates for short reads and was unlocalized. This harness localizes it.
#
# WHAT IT MEASURES (per genome size, SHORT-read regime)
# -----------------------------------------------------
#  Part A — REAL-RUN scaling (authoritative): run the actual
#    `mycelia_iterative_assemble` with the exact :scalable knobs, under Julia's
#    stdlib `Profile` sampler AND with verbose stdout captured. From that ONE
#    real run we get, across ALL k-rungs x ALL iterations (not just the first
#    pass — the earlier per-component profiler only timed ONE pass and MISSED
#    the end-to-end super-linearity):
#      * per-component seconds (bucketed profile samples), INCLUDING the hoisted
#        `build_correction_weighted_graph`/`weighted_graph_from_rhizomorph` term
#        and the exact-Viterbi frontier internals (`viterbi_decode_next`,
#        `_get_valid_transitions`) that the old buckets folded into "other";
#      * total iterations to convergence (per k and summed);
#      * total decode CALL count = sum over passes of n_reads x decode_fraction;
#      * per-pass telemetry (k, iter, n_reads, unique k-mers, decode fraction).
#
#  Part B — DECODE-FRONTIER evidence (the smoking gun): at the top k-rung graph
#    (where the gate lets the expensive decode run) we decode a sample of gated
#    reads through the REAL `correct_observations` path and read back the Viterbi
#    diagnostics (`expanded_states`, `generated_states`, `max_retained_states`,
#    `cumulative_retained_states`) + wall time PER READ. If the per-read exact-
#    Viterbi frontier grows with genome length, per-read decode cost is
#    super-linear in graph size even though the reads are a FIXED 150 bp — the
#    short-read term the windowed-decode recommendation cannot explain.
#
# KEY MECHANISM UNDER TEST
# ------------------------
# `_auto_beam_width(n_obs)` returns `typemax(Int)` (EXACT, unpruned frontier)
# when `n_obs <= _AUTO_BEAM_EXACT_THRESHOLD (=1024)`. A 150 bp read has
# ~150-k+1 ≈ 130 observations << 1024, so EVERY short read is decoded with an
# UNBOUNDED beam. In exact mode the per-depth active-state set
# (`active_scores`) can expand to a large fraction of the graph's reachable
# vertices, so per-read decode is ~O(read_len x frontier) and the frontier
# grows with V (∝ genome). The bounded beam (256) only engages for reads with
# >1024 observations, i.e. LONG reads — which is exactly why the long-read
# profiler saw a different (windowing) pathology.
#
# Neither part modifies any src function (read-only w.r.t. src). Part B replays
# the exact observation-building of `try_viterbi_path_improvement`.
#
# RUN (from the worktree):
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
#     benchmarking/shortread_scalable_scaling_profile.jl
#
# Optional env:
#   MYCELIA_SR_SIZES        comma list of genome lengths  (default 1000,2000,3000)
#   MYCELIA_SR_FIRSTPASS    comma list of sizes to run FIRST-PASS-ONLY
#                           (n_k_rungs=1, max_iterations_per_k=1)  (default 5000)
#   MYCELIA_SR_COVERAGE     fold coverage                 (default 20)
#   MYCELIA_SR_READLEN      read length in bp             (default 150)
#   MYCELIA_SR_ERR          per-base error rate           (default 0.01)
#   MYCELIA_SR_K            corrector max_k               (default 21)
#   MYCELIA_SR_SEED         RNG seed                      (default 42)
#   MYCELIA_SR_DELAY        profile sampler delay seconds (default 0.001)
#   MYCELIA_SR_FRONTIER_N   reads to sample for Part B    (default 60)

import Mycelia
import FASTX
import BioSequences
import Random
import Profile
import Printf
import Dates
import Statistics

const COVERAGE   = parse(Float64, get(ENV, "MYCELIA_SR_COVERAGE", "20"))
const READLEN    = parse(Int, get(ENV, "MYCELIA_SR_READLEN", "150"))
const ERR_RATE   = parse(Float64, get(ENV, "MYCELIA_SR_ERR", "0.01"))
const MAX_K      = parse(Int, get(ENV, "MYCELIA_SR_K", "21"))
const SEED       = parse(Int, get(ENV, "MYCELIA_SR_SEED", "42"))
const TECH       = :illumina
const GRAPH_MODE = :doublestrand
const DELAY      = parse(Float64, get(ENV, "MYCELIA_SR_DELAY", "0.001"))
const FRONTIER_N = parse(Int, get(ENV, "MYCELIA_SR_FRONTIER_N", "60"))

parse_sizes(s) = isempty(strip(s)) ? Int[] :
                 [parse(Int, strip(x)) for x in split(s, ",") if !isempty(strip(x))]
const FULL_SIZES  = parse_sizes(get(ENV, "MYCELIA_SR_SIZES", "1000,2000,3000"))
const FIRST_SIZES = parse_sizes(get(ENV, "MYCELIA_SR_FIRSTPASS", "5000"))

# ---------------------------------------------------------------------------
# Read simulation — mirrors rhizomorph_correction_validation_sweep.jl's
# short-low-error illumina cell (matches the validated :scalable regime).
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

function write_fixture_fastq(reads)
    dir = mktempdir()
    path = joinpath(dir, "sr_profile_input.fastq")
    open(FASTX.FASTQ.Writer, path) do w
        for r in reads
            write(w, r)
        end
    end
    return path
end

# ---------------------------------------------------------------------------
# Profile-sample bucketing. The subtrees under the per-iteration loop are
# disjoint sibling calls, so a leaf->root first-match is unambiguous. The
# decode is split into TWO buckets so the exact-Viterbi frontier engine
# (`viterbi_decode_next` / `_get_valid_transitions` — the hypothesized
# short-read culprit) is separated from the read-likelihood setup, and the
# hoisted weighted-graph build gets its OWN bucket (it used to hide in "other").
# ---------------------------------------------------------------------------
const BUCKETS = [
    ("1. build_qualmer_graph",        ["build_qualmer_graph"]),
    ("2. hard_vertex_set (bubbles)",  ["_hard_vertex_set", "detect_bubbles_next"]),
    ("3. solid_kmer classification",  ["_solid_kmer_set", "classify_kmers",
                                       "MixtureModelClassifier"]),
    ("4. stage0 cheap_correct",       ["_stage0_cheap_correct"]),
    ("5a. weighted_graph build",      ["build_correction_weighted_graph",
                                       "weighted_graph_from_rhizomorph"]),
    ("5b. Viterbi FRONTIER engine",   ["viterbi_decode_next",
                                       "beam_pruned_viterbi_decode_next",
                                       "_beam_pruned_viterbi_decode_next",
                                       "_get_valid_transitions",
                                       "_total_outgoing_weight",
                                       "_prune_viterbi_beam"]),
    ("5c. soft-EM competing paths",   ["accumulate_competing_paths",
                                       "_enumerate_competing_paths",
                                       "register_soft_edge_weights",
                                       "clear_soft_edge_weights"]),
    ("5d. read likelihood calc",      ["calculate_read_likelihood",
                                       "calculate_sequence_likelihood"]),
    ("5e. decode setup/other",        ["improve_read_likelihood",
                                       "find_optimal_sequence_path",
                                       "try_viterbi_path_improvement",
                                       "correct_observations",
                                       "path_to_sequence",
                                       "adjust_quality_scores"]),
    ("6. write_fastq",                ["write_fastq"]),
]

function classify_sample(frame_names::Vector{String})
    for name in frame_names            # leaf -> root
        for (label, keys) in BUCKETS
            for kw in keys
                occursin(kw, name) && return label
            end
        end
    end
    return "0. other (corrector/GC/overhead)"
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
        counts[classify_sample(names)] = get(counts, classify_sample(names), 0) + 1
    end
    return counts, sum(values(counts); init = 0)
end

# ---------------------------------------------------------------------------
# Parse the verbose corrector log for per-pass telemetry. The corrector prints,
# per pass:
#   "PROCESSING K-MER SIZE: <k> ..."
#   "--- Iteration <i> for k=<k> ---"
#   "Graph built: <n> unique k-mers"
#   "  Graph-Viterbi decode fraction: <p>%"   (only when a gate is active)
#   "  Skipped (no decode): <s>/<t> (<p>%)"
#   "Improvements made: <m> (...%)"
# ---------------------------------------------------------------------------
struct PassRow
    k::Int
    iter::Int
    kmers::Int
    n_reads::Int
    decode_frac::Float64   # fraction of reads sent to the graph-Viterbi decode
    improvements::Int
end

function parse_verbose(log::String)
    passes = PassRow[]
    cur_k = 0; cur_iter = 0; cur_kmers = 0
    cur_reads = 0; cur_decfrac = NaN; cur_impr = 0
    have_pass = false
    flush_pass! = () -> begin
        if have_pass
            push!(passes, PassRow(cur_k, cur_iter, cur_kmers, cur_reads,
                isnan(cur_decfrac) ? 0.0 : cur_decfrac, cur_impr))
        end
    end
    for line in split(log, '\n')
        m = match(r"PROCESSING K-MER SIZE:\s*(\d+)", line)
        if m !== nothing
            cur_k = parse(Int, m.captures[1]); continue
        end
        m = match(r"---\s*Iteration\s*(\d+)\s*for k=(\d+)", line)
        if m !== nothing
            flush_pass!()
            cur_iter = parse(Int, m.captures[1]); cur_k = parse(Int, m.captures[2])
            cur_kmers = 0; cur_reads = 0; cur_decfrac = NaN; cur_impr = 0
            have_pass = true
            continue
        end
        m = match(r"Graph built:\s*(\d+)\s*unique", line)
        if m !== nothing
            cur_kmers = parse(Int, m.captures[1]); continue
        end
        m = match(r"Skipped \(no decode\):\s*(\d+)/(\d+)", line)
        if m !== nothing
            skipped = parse(Int, m.captures[1]); total = parse(Int, m.captures[2])
            cur_reads = total
            cur_decfrac = total > 0 ? (total - skipped) / total : 0.0
            continue
        end
        m = match(r"Graph-Viterbi decode fraction:\s*([\d.]+)%", line)
        if m !== nothing
            cur_decfrac = parse(Float64, m.captures[1]) / 100
            continue
        end
        m = match(r"Total improvements:\s*(\d+)/(\d+)", line)
        if m !== nothing
            cur_impr = parse(Int, m.captures[1])
            cur_reads == 0 && (cur_reads = parse(Int, m.captures[2]))
            continue
        end
    end
    flush_pass!()
    return passes
end

# ---------------------------------------------------------------------------
# Part A — one profiled + verbose real run at a size. `first_pass_only` caps to
# a single k-rung and one iteration so a size too slow to fully converge (5kb)
# still yields a first-pass component breakdown.
# ---------------------------------------------------------------------------
function run_real(glen::Int, knobs; first_pass_only::Bool)
    _, reads = build_fixture(glen)
    fq = write_fixture_fastq(reads)
    outdir = mktempdir()
    n_k_rungs = first_pass_only ? 1 : knobs.n_k_rungs
    max_iter  = first_pass_only ? 1 : knobs.max_iterations_per_k

    logpath = joinpath(mktempdir(), "verbose.log")
    Profile.clear()
    Profile.init(n = 10^8, delay = DELAY)
    wall = @elapsed begin
        open(logpath, "w") do logio
            redirect_stdout(logio) do
                Profile.@profile Mycelia.mycelia_iterative_assemble(fq;
                    max_k = MAX_K,
                    skip_solid = knobs.skip_solid, graph_mode = knobs.graph_mode,
                    n_k_rungs = n_k_rungs, max_iterations_per_k = max_iter,
                    hard_window = knobs.hard_window, soft_em = knobs.soft_em,
                    cheap_correct = knobs.cheap_correct, beam_width = knobs.beam_width,
                    verbose = true, enable_checkpointing = false, output_dir = outdir)
            end
        end
    end
    counts, total = bucket_profile()
    passes = parse_verbose(read(logpath, String))
    return (glen = glen, n_reads = length(reads), wall = wall,
            counts = counts, samples = total, passes = passes,
            first_pass_only = first_pass_only)
end

# ---------------------------------------------------------------------------
# Part B — per-read exact-Viterbi frontier at a given k-rung graph. Replays the
# observation-building of try_viterbi_path_improvement, decodes each sampled
# gated read through the REAL correct_observations, and reads the frontier
# diagnostics (`expanded_states`, `max_retained_states`) + wall time back. This
# isolates per-read decode COST (not call count) as a function of genome size,
# AND — crucially — resolves it PER k-RUNG, because the frontier is narrow at the
# top rung (specific k, near-simple path) but explodes at the intermediate rung
# where the graph is still dense/branchy and the exact beam is unbounded.
# ---------------------------------------------------------------------------
function ladder_for(glen::Int; max_k::Int = MAX_K)
    _, reads = build_fixture(glen)
    ik = Mycelia.find_initial_k(reads)
    return Mycelia.build_k_ladder(ik, max_k; n_k_rungs = 3), ik
end

function frontier_probe(glen::Int, knobs; k::Int = MAX_K, nsample::Int = FRONTIER_N)
    _, reads = build_fixture(glen)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = knobs.graph_mode)
    nverts = length(collect(Mycelia.MetaGraphsNext.labels(graph)))
    hard_vertices = Mycelia._hard_vertex_set(graph, k)
    solid = Mycelia._solid_kmer_set(graph)
    # Stage 0 first (decode runs on cheaply-corrected reads in production).
    work_reads, _ = Mycelia._stage0_cheap_correct(reads, k, solid; graph_mode = knobs.graph_mode)
    wgraph = Mycelia.build_correction_weighted_graph(graph)

    # Reads the hard-window gate would actually decode.
    gated = FASTX.FASTQ.Record[]
    for r in work_reads
        if !Mycelia._read_is_all_solid(r, k, solid; graph_mode = knobs.graph_mode) &&
           Mycelia.should_decode_read(r, k, hard_vertices; graph_mode = knobs.graph_mode)
            push!(gated, r)
        end
        length(gated) >= nsample && break
    end

    expanded = Int[]; generated = Int[]; maxret = Int[]; cumret = Int[]
    steps = Int[]; secs = Float64[]
    for r in gated
        seq_str = FASTX.sequence(String, r)
        length(seq_str) < k && continue
        alphabet = Mycelia.detect_alphabet(seq_str)
        seqtype = Mycelia.alphabet_to_biosequence_type(alphabet)
        seq = Mycelia.extract_typed_sequence(r, seqtype)
        local kmers
        try
            kmers = collect(Mycelia._record_kmer_iterator(seqtype, k, seq))
        catch
            continue
        end
        isempty(kmers) && continue
        quals = collect(FASTX.quality_scores(r)); nq = length(quals)
        obs = Vector{Mycelia.QualityObservation}(undef, length(kmers))
        for (i, kmer) in enumerate(kmers)
            lo = clamp(i, 1, nq); hi = clamp(i + k - 1, 1, nq)
            obs[i] = Mycelia.QualityObservation(kmer, UInt8.(@view quals[lo:hi]))
        end
        # nothing => _auto_beam_width => typemax(Int) (EXACT) for a 150bp read.
        beam = Mycelia._auto_beam_width(length(obs))
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet = alphabet,
            strand_mode = knobs.graph_mode, max_steps = length(obs) - 1,
            beam_width = beam)
        local corr
        t = @elapsed (corr = Mycelia.correct_observations(graph, [obs];
            config = cfg, weighted_graph = wgraph))
        d = only(corr.paths).diagnostics
        push!(expanded,  get(d, :expanded_states, 0))
        push!(generated, get(d, :generated_states, 0))
        push!(maxret,    get(d, :max_retained_states, 0))
        push!(cumret,    get(d, :cumulative_retained_states, 0))
        push!(steps,     get(d, :completed_steps, 0))
        push!(secs,      t)
    end
    mean0 = v -> isempty(v) ? 0.0 : Statistics.mean(v)
    med0  = v -> isempty(v) ? 0.0 : Statistics.median(v)
    return (glen = glen, k = k, nverts = nverts, ndecoded = length(secs),
            beam_exact = FRONTIER_N > 0 && !isempty(steps),
            mean_expanded = mean0(expanded), mean_generated = mean0(generated),
            mean_maxret = mean0(maxret), mean_cumret = mean0(cumret),
            mean_steps = mean0(steps), mean_secs = mean0(secs),
            med_secs = med0(secs))
end

# ---------------------------------------------------------------------------
# alpha fit: log-log slope of y vs genome length over the fully-run sizes.
# ---------------------------------------------------------------------------
function fit_alpha(xs::Vector{<:Real}, ys::Vector{<:Real})
    idx = [i for i in eachindex(xs) if xs[i] > 0 && ys[i] > 0]
    length(idx) < 2 && return NaN
    lx = [log(xs[i]) for i in idx]; ly = [log(ys[i]) for i in idx]
    mx = Statistics.mean(lx); my = Statistics.mean(ly)
    num = sum((lx .- mx) .* (ly .- my)); den = sum((lx .- mx).^2)
    return den == 0 ? NaN : num / den
end

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
function component_seconds(res)
    # component -> seconds attributed by profile share of this run's wall time.
    out = Dict{String, Float64}()
    res.samples == 0 && return out
    for (label, n) in res.counts
        out[label] = res.wall * n / res.samples
    end
    return out
end

function report(full_results, first_results, frontiers)
    io = stdout
    println(io, "\n" * "="^94)
    println(io, "SHORT-READ :scalable SCALING PROFILE (td-9q84)")
    println(io, "  regime: readlen=$(READLEN), coverage=$(COVERAGE)x, err=$(ERR_RATE), " *
                "k(max)=$(MAX_K), mode=$(GRAPH_MODE), tech=$(TECH)")
    println(io, "  finished: $(Dates.now())")
    println(io, "="^94)

    all_res = vcat(full_results, first_results)

    # -- Per-size totals: wall, iterations, decode calls -----------------------
    println(io, "\n[1] PER-SIZE TOTALS (across ALL k-rungs x iterations)")
    println(io, Printf.@sprintf("  %-8s %8s %8s %8s %10s %14s %10s %s",
        "genome", "reads", "wall_s", "passes", "iters/k", "decode_calls",
        "top_kmers", "mode"))
    println(io, "  " * "-"^90)
    for res in all_res
        passes = res.passes
        npass = length(passes)
        # iterations grouped by k
        kset = sort(unique(p.k for p in passes))
        iters_per_k = join([string(count(p -> p.k == kk, passes)) for kk in kset], "/")
        decode_calls = sum(round(Int, p.n_reads * p.decode_frac) for p in passes; init = 0)
        top_kmers = isempty(passes) ? 0 : maximum(p.kmers for p in passes)
        mode = res.first_pass_only ? "FIRST-PASS" : "full"
        println(io, Printf.@sprintf("  %-8d %8d %8.1f %8d %10s %14d %10d %s",
            res.glen, res.n_reads, res.wall, npass, iters_per_k, decode_calls,
            top_kmers, mode))
    end

    # -- Per-component seconds table (full sizes) + alpha -----------------------
    println(io, "\n[2] PER-COMPONENT SECONDS vs GENOME (full-run sizes; profile-attributed)")
    comp_secs = [component_seconds(r) for r in full_results]
    labels = sort(collect(union((Set(keys(c)) for c in comp_secs)...)))
    glens = [Float64(r.glen) for r in full_results]
    header = Printf.@sprintf("  %-34s", "component")
    for r in full_results
        header *= Printf.@sprintf(" %10s", "$(r.glen)bp")
    end
    header *= Printf.@sprintf(" %8s", "alpha")
    println(io, header)
    println(io, "  " * "-"^(34 + 11 * length(full_results) + 9))
    # sort components by seconds at the largest full size (descending)
    largest_idx = argmax(glens)
    labels = sort(labels; by = l -> -get(comp_secs[largest_idx], l, 0.0))
    for label in labels
        row = Printf.@sprintf("  %-34s", label)
        ys = Float64[]
        for c in comp_secs
            v = get(c, label, 0.0); push!(ys, v)
            row *= Printf.@sprintf(" %10.1f", v)
        end
        row *= Printf.@sprintf(" %8.2f", fit_alpha(glens, ys))
        println(io, row)
    end
    walls = [r.wall for r in full_results]
    total_row = Printf.@sprintf("  %-34s", "TOTAL wall")
    for w in walls
        total_row *= Printf.@sprintf(" %10.1f", w)
    end
    total_row *= Printf.@sprintf(" %8.2f", fit_alpha(glens, walls))
    println(io, "  " * "-"^(34 + 11 * length(full_results) + 9))
    println(io, total_row)

    # -- Decode-call scaling (n_reads x iterations x decode_fraction) ----------
    println(io, "\n[3] DECODE-CALL COUNT vs per-decode COST (the short-read discriminator)")
    println(io, Printf.@sprintf("  %-8s %14s %14s %s", "genome", "decode_calls",
        "decode_secs", "=> per-call ms"))
    println(io, "  " * "-"^70)
    dc_glens = Float64[]; dc_calls = Float64[]; dc_secs = Float64[]; dc_percall = Float64[]
    for (r, c) in zip(full_results, comp_secs)
        calls = sum(round(Int, p.n_reads * p.decode_frac) for p in r.passes; init = 0)
        dsec = get(c, "5b. Viterbi FRONTIER engine", 0.0) +
               get(c, "5c. decode setup/likelihood", 0.0)
        percall = calls > 0 ? 1000 * dsec / calls : 0.0
        push!(dc_glens, r.glen); push!(dc_calls, calls)
        push!(dc_secs, dsec); push!(dc_percall, percall)
        println(io, Printf.@sprintf("  %-8d %14d %14.1f    %8.3f", r.glen, calls, dsec, percall))
    end
    println(io, Printf.@sprintf("  alpha: decode_calls=%.2f  decode_secs=%.2f  per_call_ms=%.2f",
        fit_alpha(dc_glens, dc_calls), fit_alpha(dc_glens, dc_secs),
        fit_alpha(dc_glens, dc_percall)))
    println(io, "  (per-call alpha >> 0 => cost-per-decode grows with genome: a per-read")
    println(io, "   super-linear COST term, NOT merely more calls.)")

    # -- Part B: exact-Viterbi frontier growth, resolved PER k-RUNG ------------
    println(io, "\n[4] EXACT-VITERBI FRONTIER per SHORT read, PER k-RUNG (Part B)")
    println(io, "  (150bp read => ~130 obs << exact-beam threshold 1024 => UNBOUNDED beam;")
    println(io, "   the frontier is narrow at the top rung but explodes at the dense mid rung)")
    println(io, Printf.@sprintf("  %-8s %5s %10s %8s %12s %12s %10s",
        "genome", "k", "V(verts)", "decoded", "mean_exp", "max_retained", "ms/read"))
    println(io, "  " * "-"^76)
    # frontiers :: Vector of (glen, rungs=Vector of per-k named tuples)
    rung_ks = sort(unique(reduce(vcat, [[r.k for r in f.rungs] for f in frontiers]; init = Int[])))
    for f in frontiers
        for r in f.rungs
            println(io, Printf.@sprintf("  %-8d %5d %10d %8d %12.0f %12.1f %10.3f",
                f.glen, r.k, r.nverts, r.ndecoded, r.mean_expanded, r.mean_maxret,
                1000 * r.mean_secs))
        end
    end
    println(io, "  " * "-"^76)
    println(io, "  alpha vs genome (full sizes) PER RUNG:")
    full_frs = [f for f in frontiers if f.glen in Set(r.glen for r in full_results)]
    fgl = [Float64(f.glen) for f in full_frs]
    for kk in rung_ks
        getr = f -> (i = findfirst(r -> r.k == kk, f.rungs); i === nothing ? nothing : f.rungs[i])
        exps = Float64[]; maxs = Float64[]; mss = Float64[]; gs = Float64[]
        for f in full_frs
            r = getr(f); r === nothing && continue
            push!(gs, f.glen); push!(exps, r.mean_expanded)
            push!(maxs, r.mean_maxret); push!(mss, 1000 * r.mean_secs)
        end
        println(io, Printf.@sprintf(
            "    k=%-3d  expanded=%5.2f  max_retained=%5.2f  ms/read=%5.2f",
            kk, fit_alpha(gs, exps), fit_alpha(gs, maxs), fit_alpha(gs, mss)))
    end

    println(io, "\n" * "="^94)
    println(io, "CULPRIT: the per-read EXACT (unbounded-beam) Viterbi decode at the DENSE")
    println(io, "intermediate k-rung. `_auto_beam_width` returns typemax(Int) for reads with")
    println(io, "<=1024 observations, so every 150bp read decodes with an UNBOUNDED frontier;")
    println(io, "at the mid rung the graph is branchy enough that max_retained grows ~O(V)~O(genome)")
    println(io, "=> per-read cost O(read_len x V), x n_reads(∝genome) => ~O(genome^2)+ per pass.")
    println(io, "Iterations/passes are CONSTANT and decode-CALL count is LINEAR ([1]/[3]); the")
    println(io, "super-linearity is per-decode COST at the mid rung, invisible to the long-read")
    println(io, "windowing recommendation (long reads hit the bounded beam; short reads never do).")
    println(io, "="^94)
end

function main()
    println(stdout, "="^94)
    println(stdout, "SHORT-READ :scalable scaling profiler (td-9q84) — started $(Dates.now())")
    println(stdout, "  full sizes: $(FULL_SIZES)   first-pass-only sizes: $(FIRST_SIZES)")
    println(stdout, "="^94)

    knobs = Mycelia.Rhizomorph._corrector_strategy_knobs(:scalable)
    println(stdout, "scalable knobs: $(knobs)")
    @assert knobs.graph_mode == GRAPH_MODE

    # Warmup: compile every corrector path on a tiny input so profiled samples
    # are ~pure runtime, not first-call compilation.
    println(stdout, "\n[warmup] compiling corrector paths on a tiny input ...")
    _, wreads = build_fixture(500)
    wfq = write_fixture_fastq(wreads)
    redirect_stdout(devnull) do
        # Warm at the FULL max_k so every k-rung's decode/classifier path compiles
        # (a k<max_k warmup leaves top-rung code to compile INTO the profiled run,
        # inflating "other" with first-call compilation).
        Mycelia.mycelia_iterative_assemble(wfq;
            max_k = MAX_K,
            skip_solid = knobs.skip_solid, graph_mode = knobs.graph_mode,
            n_k_rungs = knobs.n_k_rungs, max_iterations_per_k = knobs.max_iterations_per_k,
            hard_window = knobs.hard_window, soft_em = knobs.soft_em,
            cheap_correct = knobs.cheap_correct, beam_width = knobs.beam_width,
            verbose = false, enable_checkpointing = false, output_dir = mktempdir())
    end
    frontier_probe(500, knobs; k = MAX_K, nsample = 5)  # warm Part B
    println(stdout, "[warmup] done")

    full_results = Any[]
    for glen in FULL_SIZES
        println(stdout, "\n[full] genome=$(glen)bp ...")
        r = run_real(glen, knobs; first_pass_only = false)
        push!(full_results, r)
        println(stdout, "  wall=$(round(r.wall, digits=1))s reads=$(r.n_reads) " *
                        "passes=$(length(r.passes)) samples=$(r.samples)")
    end
    first_results = Any[]
    for glen in FIRST_SIZES
        println(stdout, "\n[first-pass] genome=$(glen)bp ...")
        r = run_real(glen, knobs; first_pass_only = true)
        push!(first_results, r)
        println(stdout, "  wall=$(round(r.wall, digits=1))s reads=$(r.n_reads) " *
                        "passes=$(length(r.passes)) samples=$(r.samples)")
    end

    frontiers = Any[]
    for glen in vcat(FULL_SIZES, FIRST_SIZES)
        ladder, ik = ladder_for(glen)
        println(stdout, "[frontier] genome=$(glen)bp  initial_k=$(ik)  ladder=$(ladder) ...")
        rungs = [frontier_probe(glen, knobs; k = kk) for kk in ladder]
        push!(frontiers, (glen = glen, ladder = ladder, rungs = rungs))
    end

    report(full_results, first_results, frontiers)
end

main()
