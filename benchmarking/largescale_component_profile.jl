# largescale_component_profile.jl
#
# Localize the MILD RESIDUAL super-linearity in the :scalable corrector at LARGE
# genome scale (bead td-firp).
#
# WHY THIS EXISTS
# ---------------
# Four scaling fixes (#371 low-k decode gating, #373 O(V+E) bubbles, #374
# linearized per-read decode, #377 re-assembly graph reuse) made the SHORT-read
# :scalable sweep LINEAR over 1-5 kb (alpha ~0.89). But stretching to real scale
# exposes a slowly-growing term only visible at large genomes:
#     5 kb  ~75 s   ->   48 kb  ~2109 s (35 min)   => alpha ~1.47 over 5k->48k.
# The run COMPLETES (tractable) but we want to LOCALIZE the alpha~1.47 term and
# push 48 kb back toward minutes.
#
# The prior short-read profiler (benchmarking/shortread_scalable_scaling_profile.jl)
# localized the 1-5 kb behavior to the per-read EXACT Viterbi frontier during
# CORRECTION. This harness asks a DIFFERENT question at a DIFFERENT scale: with
# per-read decode now linearized, WHICH end-to-end component carries the residual
# alpha~1.47 at 48 kb? Candidates (from the bead):
#   (1) CONTIG-EXTRACTION  — find_contigs_next / find_linear_path over the large
#       (~380k-vertex) reassembly graph. find_linear_path filters neighbors with
#       `!(n in path)` where `path` is a Vector => O(L) membership per step =>
#       O(L^2) per contig of length L. On a near-linear corrected-read graph L is
#       ~genome-length, so this is the prime O(V^2) suspect.
#   (2) RC-DEDUP           — the post-hoc _dedup_reverse_complements pass. Uses a
#       Set (O(1) membership) => O(n) in contig count; measured here to CONFIRM it
#       is NOT quadratic (and, on the default qualmer reassembly path, not even
#       invoked — dedup_revcomp defaults false).
#   (3) RESIDUAL PER-READ DECODE — a decode COST term that still grows with graph
#       size at large scale (frontier probe, Part D).
#   (4) RE-ASSEMBLY GRAPH BUILD  — build_qualmer_graph over the larger read set.
#
# METHOD (read-only w.r.t. src; no library function modified or monkey-patched)
# ----------------------------------------------------------------------------
# Part B (CORE, cheap, ALL sizes incl. 48 kb): directly time the reassembly
#   sub-components on an ERROR-FREE read set (the near-perfect CORRECTED-read
#   graph the reassembly actually runs on is well-approximated by error-free
#   reads; that graph is near-linear => the WORST case for find_linear_path's
#   O(L^2) walk). For each size we @elapsed:
#     t_graph   = build_qualmer_graph(reads, k; mode=:doublestrand)
#     t_find    = find_contigs_next(graph; min_contig_length=k+1)   [CONTIG-EXTRACT]
#     t_dedup   = _dedup_reverse_complements(contig_seqs)           [RC-DEDUP]
#     t_consens = per-path _qualmer_path_to_consensus_fastq         [consensus]
#   and record V (vertices), n_contigs, max/median contig length. A raw
#   (err=0.01) read set is timed alongside for contrast. This part fits alpha per
#   sub-component across the FULL size range — the decisive measurement, since it
#   reaches 48 kb without paying the 35-min correction cost.
#
# Part A (full sizes only, 3/5/10 kb): run the REAL iterative arm end to end
#   (correction = mycelia_iterative_assemble with the :scalable knobs; reassembly
#   = assemble_genome(corrected; corrector=:none)), each @elapsed, to show the
#   end-to-end correction-vs-reassembly split and confirm which grows.
#
# Part C (large sizes, 20/48 kb): FIRST-ITERATION-ONLY correction (n_k_rungs=1,
#   max_iterations_per_k=1) under the Profile sampler — a SINGLE instrumented
#   pass, NOT the 35-min full run — to read the correction-side per-component
#   breakdown (incl. per-read decode) at large graph size.
#
# Part D (all sizes): per-read exact-Viterbi frontier probe (mean ms/read,
#   mean_retained) over the large graph — checks for a residual decode-COST term.
#
# RUN (from the worktree):
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
#     benchmarking/largescale_component_profile.jl
#
# Optional env (defaults chosen so the whole sweep finishes in ~10-25 min):
#   MYCELIA_LC_FULL_SIZES   full e2e correction+reassembly sizes (default 3000,5000,10000)
#   MYCELIA_LC_LARGE_SIZES  first-pass-only correction sizes     (default 20000,48000)
#   MYCELIA_LC_REASM_SIZES  reassembly-component probe sizes      (default 3000,5000,10000,20000,48000)
#   MYCELIA_LC_COVERAGE     fold coverage                         (default 20)
#   MYCELIA_LC_READLEN      read length, bp                       (default 150)
#   MYCELIA_LC_ERR          per-base error rate (raw reads)       (default 0.01)
#   MYCELIA_LC_K            corrector/assembler max_k             (default 21)
#   MYCELIA_LC_SEED         RNG seed                              (default 42)
#   MYCELIA_LC_DELAY        profile sampler delay seconds         (default 0.001)
#   MYCELIA_LC_FRONTIER_N   reads to sample for Part D            (default 40)
#   MYCELIA_LC_SKIP_A       truthy -> skip Part A (full e2e)      (default false)
#   MYCELIA_LC_SKIP_C       truthy -> skip Part C (first-pass)    (default false)

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
_parse_sizes(s) = isempty(strip(s)) ? Int[] :
    Int[parse(Int, strip(x)) for x in split(s, ",") if !isempty(strip(x))]
_truthy(s) = lowercase(strip(s)) in ("1", "true", "yes", "on")

const FULL_SIZES  = _parse_sizes(get(ENV, "MYCELIA_LC_FULL_SIZES", "3000,5000,10000"))
const LARGE_SIZES = _parse_sizes(get(ENV, "MYCELIA_LC_LARGE_SIZES", "20000,48000"))
const REASM_SIZES = _parse_sizes(get(ENV, "MYCELIA_LC_REASM_SIZES", "3000,5000,10000,20000,48000"))
const COVERAGE    = parse(Float64, get(ENV, "MYCELIA_LC_COVERAGE", "20"))
const READLEN     = parse(Int, get(ENV, "MYCELIA_LC_READLEN", "150"))
const ERR_RATE    = parse(Float64, get(ENV, "MYCELIA_LC_ERR", "0.01"))
const MAX_K       = parse(Int, get(ENV, "MYCELIA_LC_K", "21"))
const SEED        = parse(Int, get(ENV, "MYCELIA_LC_SEED", "42"))
const DELAY       = parse(Float64, get(ENV, "MYCELIA_LC_DELAY", "0.001"))
const FRONTIER_N  = parse(Int, get(ENV, "MYCELIA_LC_FRONTIER_N", "40"))
const TECH        = :illumina
const GRAPH_MODE  = :doublestrand
const SKIP_A      = _truthy(get(ENV, "MYCELIA_LC_SKIP_A", "false"))
const SKIP_C      = _truthy(get(ENV, "MYCELIA_LC_SKIP_C", "false"))

# ---------------------------------------------------------------------------
# Read simulation — mirrors the validated short-low-error illumina cell used by
# rhizomorph_correction_validation_sweep.jl / the prior scaling profilers.
# `error_rate=0.0` yields the error-free "corrected-read proxy" read set.
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

function build_fixture(glen::Int; error_rate::Float64 = ERR_RATE)
    rng = Random.MersenneTwister(SEED)
    ref_rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = glen)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, ref_rec)
    reads = simulate_reads(refseq, READLEN, COVERAGE, error_rate, TECH, rng)
    return refseq, reads
end

function write_fixture_fastq(reads)
    dir = mktempdir()
    path = joinpath(dir, "lc_profile_input.fastq")
    open(FASTX.FASTQ.Writer, path) do w
        for r in reads
            write(w, r)
        end
    end
    return path
end

# ---------------------------------------------------------------------------
# alpha fit: log-log slope of y vs genome length over the given (x,y) points.
# ---------------------------------------------------------------------------
function fit_alpha(xs::Vector{<:Real}, ys::Vector{<:Real})
    idx = [i for i in eachindex(xs) if xs[i] > 0 && ys[i] > 0]
    length(idx) < 2 && return NaN
    lx = [log(xs[i]) for i in idx]; ly = [log(ys[i]) for i in idx]
    mx = Statistics.mean(lx); my = Statistics.mean(ly)
    num = sum((lx .- mx) .* (ly .- my)); den = sum((lx .- mx) .^ 2)
    return den == 0 ? NaN : num / den
end

# alpha over only the LARGE end (the two largest sizes) — the residual term the
# whole-range alpha can under-report if the small sizes are noise-dominated.
function fit_alpha_tail(xs::Vector{<:Real}, ys::Vector{<:Real})
    length(xs) < 2 && return NaN
    o = sortperm(xs)
    return fit_alpha(xs[o[end-1:end]], ys[o[end-1:end]])
end

# ===========================================================================
# Part B — reassembly sub-component direct timing (CORE, all sizes incl. 48 kb)
# ===========================================================================
struct ReasmRow
    glen::Int
    label::String        # "err-free" | "raw"
    nreads::Int
    nverts::Int
    ncontigs::Int
    max_contig::Int
    med_contig::Int
    t_graph::Float64
    t_find::Float64      # find_contigs_next  [CONTIG-EXTRACTION]
    t_dedup::Float64     # _dedup_reverse_complements  [RC-DEDUP]
    t_consensus::Float64 # per-path consensus fastq
end

function reasm_probe(glen::Int, label::String, reads)
    k = MAX_K
    local graph
    t_graph = @elapsed (graph = Mycelia.Rhizomorph.build_qualmer_graph(
        reads, k; mode = GRAPH_MODE))
    nverts = length(collect(Mycelia.MetaGraphsNext.labels(graph)))

    local contig_paths
    t_find = @elapsed (contig_paths = Mycelia.Rhizomorph.find_contigs_next(
        graph; min_contig_length = k + 1))
    ncontigs = length(contig_paths)
    lens = Int[length(cp.sequence) for cp in contig_paths]
    max_contig = isempty(lens) ? 0 : maximum(lens)
    med_contig = isempty(lens) ? 0 : round(Int, Statistics.median(lens))

    contig_seqs = String[string(cp.sequence) for cp in contig_paths]
    local _deduped
    t_dedup = @elapsed (_deduped = Mycelia.Rhizomorph._dedup_reverse_complements(contig_seqs))

    t_consensus = @elapsed begin
        for (i, cp) in enumerate(contig_paths)
            isempty(cp.vertices) && continue
            Mycelia.Rhizomorph._qualmer_path_to_consensus_fastq(
                cp.vertices, graph, "contig_$i")
        end
    end

    return ReasmRow(glen, label, length(reads), nverts, ncontigs, max_contig,
        med_contig, t_graph, t_find, t_dedup, t_consensus)
end

# ===========================================================================
# Part A — full end-to-end correction + reassembly at a size (@elapsed each).
# ===========================================================================
struct E2ERow
    glen::Int
    nreads::Int
    ncorrected::Int
    t_correct::Float64
    t_reassemble::Float64
end

function e2e_cell(glen::Int, knobs)
    _, reads = build_fixture(glen; error_rate = ERR_RATE)
    fq = write_fixture_fastq(reads)
    outdir = mktempdir()
    corrector_graph_mode = knobs.graph_mode === nothing ? GRAPH_MODE : knobs.graph_mode

    local result_dict
    t_correct = @elapsed redirect_stdout(devnull) do
        result_dict = Mycelia.mycelia_iterative_assemble(fq;
            max_k = MAX_K, skip_solid = knobs.skip_solid, graph_mode = corrector_graph_mode,
            n_k_rungs = knobs.n_k_rungs, max_iterations_per_k = knobs.max_iterations_per_k,
            hard_window = knobs.hard_window, soft_em = knobs.soft_em,
            cheap_correct = knobs.cheap_correct, beam_width = knobs.beam_width,
            verbose = false, enable_checkpointing = false, output_dir = outdir)
    end

    corrected_fastq = get(result_dict[:metadata], :final_fastq_file, nothing)
    corrected_reads = open(FASTX.FASTQ.Reader, corrected_fastq) do reader
        collect(reader)
    end

    local iter_asm
    t_reassemble = @elapsed redirect_stdout(devnull) do
        iter_asm = Mycelia.Rhizomorph.assemble_genome(
            corrected_reads; k = MAX_K, corrector = :none, verbose = false)
    end

    return E2ERow(glen, length(reads), length(corrected_reads), t_correct, t_reassemble)
end

# ===========================================================================
# Part C — first-iteration-only correction under the Profile sampler.
# Buckets mirror the short-read profiler so the correction-side breakdown
# (incl. per-read decode frontier) is directly comparable.
# ===========================================================================
const BUCKETS = [
    ("1. build_qualmer_graph",       ["build_qualmer_graph"]),
    ("2. hard_vertex_set (bubbles)", ["_hard_vertex_set", "detect_bubbles_next"]),
    ("3. solid_kmer classification", ["_solid_kmer_set", "classify_kmers",
                                      "MixtureModelClassifier"]),
    ("4. stage0 cheap_correct",      ["_stage0_cheap_correct"]),
    ("5a. weighted_graph build",     ["build_correction_weighted_graph",
                                      "weighted_graph_from_rhizomorph"]),
    ("5b. Viterbi FRONTIER engine",  ["viterbi_decode_next",
                                      "beam_pruned_viterbi_decode_next",
                                      "_beam_pruned_viterbi_decode_next",
                                      "_get_valid_transitions",
                                      "_total_outgoing_weight",
                                      "_prune_viterbi_beam"]),
    ("5c. soft-EM competing paths",  ["accumulate_competing_paths",
                                      "_enumerate_competing_paths",
                                      "register_soft_edge_weights",
                                      "clear_soft_edge_weights"]),
    ("5d. read likelihood calc",     ["calculate_read_likelihood",
                                      "calculate_sequence_likelihood"]),
    ("5e. decode setup/other",       ["improve_read_likelihood",
                                      "find_optimal_sequence_path",
                                      "try_viterbi_path_improvement",
                                      "correct_observations",
                                      "path_to_sequence",
                                      "adjust_quality_scores"]),
    ("6. find_contigs (extract)",    ["find_contigs_next", "find_linear_path",
                                      "_select_linear_neighbor", "_dominant_neighbor"]),
    ("7. reverse-complement dedup",  ["_dedup_reverse_complements", "_canonical_string"]),
    ("8. write_fastq",               ["write_fastq"]),
]

function classify_sample(frame_names::Vector{String})
    for name in frame_names
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
        lbl = classify_sample(names)
        counts[lbl] = get(counts, lbl, 0) + 1
    end
    return counts, sum(values(counts); init = 0)
end

function firstpass_correction(glen::Int, knobs)
    _, reads = build_fixture(glen; error_rate = ERR_RATE)
    fq = write_fixture_fastq(reads)
    outdir = mktempdir()
    corrector_graph_mode = knobs.graph_mode === nothing ? GRAPH_MODE : knobs.graph_mode

    Profile.clear()
    Profile.init(n = 10^8, delay = DELAY)
    wall = @elapsed redirect_stdout(devnull) do
        Profile.@profile Mycelia.mycelia_iterative_assemble(fq;
            max_k = MAX_K, skip_solid = knobs.skip_solid, graph_mode = corrector_graph_mode,
            n_k_rungs = 1, max_iterations_per_k = 1,
            hard_window = knobs.hard_window, soft_em = knobs.soft_em,
            cheap_correct = knobs.cheap_correct, beam_width = knobs.beam_width,
            verbose = false, enable_checkpointing = false, output_dir = outdir)
    end
    counts, total = bucket_profile()
    return (glen = glen, nreads = length(reads), wall = wall, counts = counts, samples = total)
end

# ===========================================================================
# Part D — per-read exact-Viterbi frontier probe over the large graph.
# (Replays try_viterbi_path_improvement's observation build; read-only.)
# ===========================================================================
function frontier_probe(glen::Int, knobs; k::Int = MAX_K, nsample::Int = FRONTIER_N)
    _, reads = build_fixture(glen; error_rate = ERR_RATE)
    graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = knobs.graph_mode)
    nverts = length(collect(Mycelia.MetaGraphsNext.labels(graph)))
    hard_vertices = Mycelia._hard_vertex_set(graph, k)
    solid = Mycelia._solid_kmer_set(graph)
    work_reads, _ = Mycelia._stage0_cheap_correct(reads, k, solid; graph_mode = knobs.graph_mode)
    wgraph = Mycelia.build_correction_weighted_graph(graph)

    gated = FASTX.FASTQ.Record[]
    for r in work_reads
        if !Mycelia._read_is_all_solid(r, k, solid; graph_mode = knobs.graph_mode) &&
           Mycelia.should_decode_read(r, k, hard_vertices; graph_mode = knobs.graph_mode)
            push!(gated, r)
        end
        length(gated) >= nsample && break
    end

    maxret = Int[]; steps = Int[]; secs = Float64[]
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
        beam = Mycelia._auto_beam_width(length(obs))
        cfg = Mycelia.ViterbiCorrectionConfig(alphabet = alphabet,
            strand_mode = knobs.graph_mode, max_steps = length(obs) - 1, beam_width = beam)
        local corr
        t = @elapsed (corr = Mycelia.correct_observations(graph, [obs];
            config = cfg, weighted_graph = wgraph))
        d = only(corr.paths).diagnostics
        push!(maxret, get(d, :max_retained_states, 0))
        push!(steps, get(d, :completed_steps, 0))
        push!(secs, t)
    end
    mean0 = v -> isempty(v) ? 0.0 : Statistics.mean(v)
    return (glen = glen, k = k, nverts = nverts, ndecoded = length(secs),
        mean_maxret = mean0(maxret), mean_steps = mean0(steps),
        mean_ms = 1000 * mean0(secs))
end

# ===========================================================================
# Reporting
# ===========================================================================
function report(reasm_rows, e2e_rows, firstpass_rows, frontier_rows)
    io = stdout
    println(io, "\n" * "="^96)
    println(io, "LARGE-SCALE COMPONENT PROFILE — :scalable residual super-linearity (td-firp)")
    println(io, "  regime: readlen=$(READLEN), coverage=$(COVERAGE)x, raw_err=$(ERR_RATE), " *
                "k=$(MAX_K), mode=$(GRAPH_MODE)")
    println(io, "  finished: $(Dates.now())")
    println(io, "="^96)

    # -- Part B: reassembly sub-component scaling (the decisive table) ----------
    println(io, "\n[B] REASSEMBLY SUB-COMPONENT SECONDS vs GENOME (direct @elapsed)")
    println(io, "    err-free reads = corrected-read proxy (near-linear graph = worst case for")
    println(io, "    find_linear_path's O(L^2) walk); raw reads shown for contrast.")
    for label in ("err-free", "raw")
        rows = sort([r for r in reasm_rows if r.label == label]; by = r -> r.glen)
        isempty(rows) && continue
        println(io, "\n  -- $(label) reads --")
        println(io, Printf.@sprintf("  %-8s %8s %9s %8s %9s %9s %9s %9s %11s",
            "glen", "V", "ncontig", "maxLen", "t_graph", "t_find", "t_dedup", "t_consen", "t_reasm_tot"))
        println(io, "  " * "-"^92)
        for r in rows
            tot = r.t_graph + r.t_find + r.t_dedup + r.t_consensus
            println(io, Printf.@sprintf("  %-8d %8d %9d %8d %9.2f %9.2f %9.3f %9.2f %11.2f",
                r.glen, r.nverts, r.ncontigs, r.max_contig,
                r.t_graph, r.t_find, r.t_dedup, r.t_consensus, tot))
        end
        gl = Float64[r.glen for r in rows]
        for (name, sel) in (("t_graph (reassembly build)", r -> r.t_graph),
                            ("t_find  (CONTIG-EXTRACTION)", r -> r.t_find),
                            ("t_dedup (RC-DEDUP)", r -> r.t_dedup),
                            ("t_consensus", r -> r.t_consensus),
                            ("V (vertices)", r -> Float64(r.nverts)),
                            ("maxLen (longest contig)", r -> Float64(r.max_contig)))
            ys = Float64[sel(r) for r in rows]
            println(io, Printf.@sprintf("    alpha[%-28s] full-range=%5.2f  large-tail=%5.2f",
                name, fit_alpha(gl, ys), fit_alpha_tail(gl, ys)))
        end
    end

    # -- Part A: full e2e correction vs reassembly ------------------------------
    if !isempty(e2e_rows)
        println(io, "\n[A] FULL END-TO-END correction vs reassembly (@elapsed; full sizes)")
        println(io, Printf.@sprintf("  %-8s %8s %10s %12s %12s %10s",
            "glen", "reads", "correct_s", "reassemble_s", "e2e_s", "reasm%"))
        println(io, "  " * "-"^70)
        rows = sort(e2e_rows; by = r -> r.glen)
        for r in rows
            tot = r.t_correct + r.t_reassemble
            pct = tot == 0 ? 0.0 : 100 * r.t_reassemble / tot
            println(io, Printf.@sprintf("  %-8d %8d %10.1f %12.1f %12.1f %9.1f%%",
                r.glen, r.nreads, r.t_correct, r.t_reassemble, tot, pct))
        end
        gl = Float64[r.glen for r in rows]
        println(io, Printf.@sprintf("  alpha: correct=%.2f  reassemble=%.2f  e2e=%.2f",
            fit_alpha(gl, Float64[r.t_correct for r in rows]),
            fit_alpha(gl, Float64[r.t_reassemble for r in rows]),
            fit_alpha(gl, Float64[r.t_correct + r.t_reassemble for r in rows])))
    end

    # -- Part C: first-pass correction component breakdown at large scale -------
    if !isempty(firstpass_rows)
        println(io, "\n[C] FIRST-ITERATION-ONLY correction, per-component seconds (large sizes)")
        println(io, "    (single instrumented pass; NOT the full 35-min run — decode-side check)")
        labels = sort(collect(union((Set(keys(r.counts)) for r in firstpass_rows)...)))
        comp = [Dict(l => (r.samples == 0 ? 0.0 : r.wall * get(r.counts, l, 0) / r.samples)
                     for l in labels) for r in firstpass_rows]
        rows = sort(collect(zip(firstpass_rows, comp)); by = t -> t[1].glen)
        hdr = Printf.@sprintf("  %-34s", "component")
        for (r, _) in rows
            hdr *= Printf.@sprintf(" %10s", "$(r.glen)bp")
        end
        hdr *= Printf.@sprintf(" %8s", "alpha")
        println(io, hdr)
        println(io, "  " * "-"^(34 + 11 * length(rows) + 9))
        gl = Float64[r.glen for (r, _) in rows]
        labels = sort(labels; by = l -> -get(rows[end][2], l, 0.0))
        for l in labels
            line = Printf.@sprintf("  %-34s", l)
            ys = Float64[]
            for (_, c) in rows
                v = get(c, l, 0.0); push!(ys, v); line *= Printf.@sprintf(" %10.1f", v)
            end
            line *= Printf.@sprintf(" %8.2f", fit_alpha(gl, ys))
            println(io, line)
        end
        wl = Printf.@sprintf("  %-34s", "TOTAL wall (1 pass)")
        for (r, _) in rows
            wl *= Printf.@sprintf(" %10.1f", r.wall)
        end
        wl *= Printf.@sprintf(" %8.2f", fit_alpha(gl, Float64[r.wall for (r, _) in rows]))
        println(io, "  " * "-"^(34 + 11 * length(rows) + 9))
        println(io, wl)
    end

    # -- Part D: per-read decode frontier over the large graph ------------------
    if !isempty(frontier_rows)
        println(io, "\n[D] PER-READ exact-Viterbi decode over the large graph (residual-decode check)")
        println(io, Printf.@sprintf("  %-8s %10s %9s %12s %12s",
            "glen", "V", "decoded", "max_retain", "ms/read"))
        println(io, "  " * "-"^58)
        rows = sort(frontier_rows; by = r -> r.glen)
        for r in rows
            println(io, Printf.@sprintf("  %-8d %10d %9d %12.1f %12.3f",
                r.glen, r.nverts, r.ndecoded, r.mean_maxret, r.mean_ms))
        end
        gl = Float64[r.glen for r in rows]
        println(io, Printf.@sprintf("  alpha: max_retained=%.2f  ms/read=%.2f",
            fit_alpha(gl, Float64[r.mean_maxret for r in rows]),
            fit_alpha(gl, Float64[r.mean_ms for r in rows])))
    end

    println(io, "\n" * "="^96)
    println(io, "See the CULPRIT readout printed after this table.")
    println(io, "="^96)
end

# ===========================================================================
# Main
# ===========================================================================
function main()
    println(stdout, "="^96)
    println(stdout, "large-scale component profiler (td-firp) — started $(Dates.now())")
    println(stdout, "  full(e2e) sizes: $(FULL_SIZES)")
    println(stdout, "  large(first-pass) sizes: $(LARGE_SIZES)")
    println(stdout, "  reassembly-probe sizes: $(REASM_SIZES)")
    println(stdout, "="^96)

    knobs = Mycelia.Rhizomorph._corrector_strategy_knobs(:scalable)
    println(stdout, "scalable knobs: $(knobs)")

    # --- Warmup: compile every path we time on a tiny input ---
    println(stdout, "\n[warmup] compiling all timed paths on a 500 bp fixture ...")
    _, wreads = build_fixture(500; error_rate = ERR_RATE)
    _, wreads0 = build_fixture(500; error_rate = 0.0)
    reasm_probe(500, "warm", wreads0)
    redirect_stdout(devnull) do
        Mycelia.mycelia_iterative_assemble(write_fixture_fastq(wreads);
            max_k = MAX_K, skip_solid = knobs.skip_solid, graph_mode = knobs.graph_mode,
            n_k_rungs = 1, max_iterations_per_k = 1, hard_window = knobs.hard_window,
            soft_em = knobs.soft_em, cheap_correct = knobs.cheap_correct,
            beam_width = knobs.beam_width, verbose = false,
            enable_checkpointing = false, output_dir = mktempdir())
    end
    frontier_probe(500, knobs; nsample = 5)
    println(stdout, "[warmup] done")

    # --- Part B: reassembly-component probe (CORE, cheap, all sizes) ---
    reasm_rows = ReasmRow[]
    for glen in REASM_SIZES
        for (label, er) in (("err-free", 0.0), ("raw", ERR_RATE))
            _, reads = build_fixture(glen; error_rate = er)
            println(stdout, "[B reasm] glen=$(glen) $(label) reads=$(length(reads)) ...")
            r = reasm_probe(glen, label, reads)
            push!(reasm_rows, r)
            println(stdout, Printf.@sprintf(
                "  V=%d ncontig=%d maxLen=%d | graph=%.2fs find=%.2fs dedup=%.3fs consen=%.2fs",
                r.nverts, r.ncontigs, r.max_contig, r.t_graph, r.t_find, r.t_dedup, r.t_consensus))
        end
    end

    # --- Part A: full e2e correction+reassembly (full sizes) ---
    e2e_rows = E2ERow[]
    if !SKIP_A
        for glen in FULL_SIZES
            println(stdout, "[A e2e] glen=$(glen) ...")
            r = e2e_cell(glen, knobs)
            push!(e2e_rows, r)
            println(stdout, Printf.@sprintf("  correct=%.1fs reassemble=%.1fs",
                r.t_correct, r.t_reassemble))
        end
    end

    # --- Part C: first-pass-only correction (large sizes) ---
    firstpass_rows = Any[]
    if !SKIP_C
        for glen in LARGE_SIZES
            println(stdout, "[C first-pass] glen=$(glen) ...")
            r = firstpass_correction(glen, knobs)
            push!(firstpass_rows, r)
            println(stdout, Printf.@sprintf("  wall(1 pass)=%.1fs samples=%d",
                r.wall, r.samples))
        end
    end

    # --- Part D: per-read decode frontier (all reassembly sizes) ---
    frontier_rows = Any[]
    for glen in REASM_SIZES
        println(stdout, "[D frontier] glen=$(glen) ...")
        push!(frontier_rows, frontier_probe(glen, knobs))
    end

    report(reasm_rows, e2e_rows, firstpass_rows, frontier_rows)

    # --- CULPRIT readout: identify the reassembly sub-component with alpha>1.3 ---
    println(stdout, "\n" * "="^96)
    println(stdout, "CULPRIT READOUT")
    println(stdout, "="^96)
    ef = sort([r for r in reasm_rows if r.label == "err-free"]; by = r -> r.glen)
    if length(ef) >= 2
        gl = Float64[r.glen for r in ef]
        comps = (("build_qualmer_graph (reassembly build)", Float64[r.t_graph for r in ef]),
                 ("find_contigs_next (CONTIG-EXTRACTION)",   Float64[r.t_find for r in ef]),
                 ("_dedup_reverse_complements (RC-DEDUP)",   Float64[r.t_dedup for r in ef]),
                 ("consensus fastq",                          Float64[r.t_consensus for r in ef]))
        worst_name = ""; worst_alpha = -Inf; worst_secs = 0.0
        for (name, ys) in comps
            a = fit_alpha(gl, ys)
            isnan(a) && continue
            if a > worst_alpha
                worst_alpha = a; worst_name = name; worst_secs = ys[end]
            end
        end
        println(stdout, Printf.@sprintf(
            "  Reassembly sub-component with the HIGHEST scaling exponent:\n    %s\n    alpha=%.2f, %.2fs at glen=%d (err-free/corrected-read proxy).",
            worst_name, worst_alpha, worst_secs, Int(gl[end])))
        println(stdout, "  (Compare against Part B alpha rows + Part A/C/D above for the end-to-end picture.)")
    end
    println(stdout, "  finished: $(Dates.now())")
    println(stdout, "="^96)
end

main()
