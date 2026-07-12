# scalable_corrector_profile.jl
#
# Per-iteration time breakdown for the :scalable corrector (corrector=:iterative,
# strategy=:scalable) — see bead td-9q84.
#
# WHY: :scalable produces near-complete assemblies but is slow (~552s for a 1kb
# genome, times out at 5kb). The cost is per-iteration overhead x k-ladder x
# iterations, NOT exponential blow-up (the Viterbi beam is bounded). We need to
# know WHICH per-iteration component dominates before optimizing.
#
# The per-iteration loop (src/iterative-assembly.jl::mycelia_iterative_assemble)
# does, per k per iteration (roughly):
#   1. build_qualmer_graph        — rebuild the qualmer graph from current reads
#   2. _hard_vertex_set           — hard-window vertex set (detect_bubbles_next + ...)
#   3. _solid_kmer_set            — coverage-based solid/weak k-mer classification
#   4. _stage0_cheap_correct      — linear k-mer-spectrum scan over ALL reads
#   5. per-read Viterbi decode    — improve_read_likelihood on the ~hard reads
#   6. write_fastq                — write corrected FASTQ (graph rebuild = #1 next iter)
#
# METHOD (non-invasive — corrector behavior UNCHANGED):
#   * Method 1 (primary): run the REAL mycelia_iterative_assemble with the exact
#     :scalable knobs under Julia's stdlib `Profile` sampler, then bucket the
#     samples by the component entry-function on each sample's call stack. This
#     aggregates seconds+% across ALL k-rungs and iterations automatically.
#   * Method 2 (cross-check): directly @elapsed each component for ONE
#     representative iteration (build graph at the initial k, then time each
#     stage on the real graph/reads). Independent seconds, no shared state with
#     Method 1.
#
# Neither method modifies the functions under test.
#
# RUN (from repo root, worktree, or a checkout with this branch):
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. benchmarking/scalable_corrector_profile.jl
#
# Optional env:
#   MYCELIA_SCP_GENOME_LEN   reference length in bp   (default 1000)
#   MYCELIA_SCP_COVERAGE     fold coverage            (default 20)
#   MYCELIA_SCP_READLEN      read length in bp        (default 100)
#   MYCELIA_SCP_ERR          per-base error rate      (default 0.01)
#   MYCELIA_SCP_K            corrector max_k          (default 21, matches the sweep)
#   MYCELIA_SCP_SEED         RNG seed                 (default 42)
#   MYCELIA_SCP_PROFILE_DELAY sampler delay seconds   (default 0.002)

import Mycelia
import FASTX
import BioSequences
import Random
import Profile
import Printf
import Dates

# ---------------------------------------------------------------------------
# Config (1kb fixture: err=0.01, 20x, 100bp reads, Illumina — reuses the sweep's
# simulation recipe via Mycelia.random_fasta_record + Mycelia.observe).
# ---------------------------------------------------------------------------
const GENOME_LEN = parse(Int, get(ENV, "MYCELIA_SCP_GENOME_LEN", "1000"))
const COVERAGE   = parse(Float64, get(ENV, "MYCELIA_SCP_COVERAGE", "20"))
const READLEN    = parse(Int, get(ENV, "MYCELIA_SCP_READLEN", "100"))
const ERR_RATE   = parse(Float64, get(ENV, "MYCELIA_SCP_ERR", "0.01"))
const MAX_K      = parse(Int, get(ENV, "MYCELIA_SCP_K", "21"))
const SEED       = parse(Int, get(ENV, "MYCELIA_SCP_SEED", "42"))
const TECH       = :illumina
const PROFILE_DELAY = parse(Float64, get(ENV, "MYCELIA_SCP_PROFILE_DELAY", "0.002"))

# ---------------------------------------------------------------------------
# Read simulation — mirrors benchmarking/rhizomorph_correction_validation_sweep.jl
# ::simulate_regime_reads so the fixture matches the validated :scalable regime.
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

function build_fixture()
    rng = Random.MersenneTwister(SEED)
    ref_rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = GENOME_LEN)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, ref_rec)
    reads = simulate_reads(refseq, READLEN, COVERAGE, ERR_RATE, TECH, rng)
    return refseq, reads
end

# Build a fixture at an arbitrary genome length (for the scaling sweep). Uses the
# same simulation recipe as build_fixture but with a per-size seed derived from
# the base SEED so each size is deterministic and independent.
function build_fixture_at(genome_len::Int)
    rng = Random.MersenneTwister(SEED)
    ref_rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = genome_len)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, ref_rec)
    reads = simulate_reads(refseq, READLEN, COVERAGE, ERR_RATE, TECH, rng)
    return refseq, reads
end

function write_fixture_fastq(reads)
    dir = mktempdir()
    path = joinpath(dir, "scalable_profile_input.fastq")
    open(FASTX.FASTQ.Writer, path) do w
        for r in reads
            write(w, r)
        end
    end
    return path
end

# ---------------------------------------------------------------------------
# Component buckets. Each bucket lists the entry-function name(s) whose presence
# on a sample's call stack attributes that sample to the component. The subtrees
# are disjoint under the per-iteration loop (build_qualmer_graph, _hard_vertex_set,
# _solid_kmer_set, _stage0_cheap_correct, and the improve_read_likelihood decode
# are sibling calls), so a leaf->root first-match is unambiguous.
# ---------------------------------------------------------------------------
const BUCKETS = [
    ("1. build_qualmer_graph", ["build_qualmer_graph"]),
    ("2. hard_vertex_set (bubbles)", ["_hard_vertex_set", "detect_bubbles_next"]),
    ("3. solid_kmer classification", ["_solid_kmer_set", "classify_kmers"]),
    ("4. stage0 cheap_correct", ["_stage0_cheap_correct"]),
    ("5. per-read Viterbi decode", ["improve_read_likelihood", "find_optimal_sequence_path",
        "try_viterbi_path_improvement", "calculate_read_likelihood"]),
    ("6. write_fastq", ["write_fastq"]),
]

# Match a sample's ordered (leaf-first) frame names to a bucket. Returns the
# bucket label, or "0. other (corrector/GC/overhead)" if no entry function matched.
function classify_sample(frame_names::Vector{String})
    for name in frame_names            # leaf -> root
        for (label, keys) in BUCKETS
            for kw in keys
                if occursin(kw, name)
                    return label
                end
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

function bucket_profile(wall_seconds::Float64)
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

function print_table(io, counts::Dict{String, Int}, total::Int, wall_seconds::Float64)
    println(io, "\nMethod 1 — Profile-sampled breakdown (aggregated over ALL k-rungs x iterations)")
    println(io, "  total wall time (mycelia_iterative_assemble): $(Printf.@sprintf("%.1f", wall_seconds)) s")
    println(io, "  total profile samples: $total  (delay=$(PROFILE_DELAY)s)")
    println(io, "")
    println(io, Printf.@sprintf("  %-38s %12s %10s %9s", "component", "seconds", "% runtime", "samples"))
    println(io, "  " * "-"^72)
    labels = sort(collect(keys(counts)))
    # Put "0. other" last visually while keeping numeric-prefix ordering
    for label in labels
        n = counts[label]
        pct = total == 0 ? 0.0 : 100 * n / total
        secs = wall_seconds * (total == 0 ? 0.0 : n / total)
        println(io, Printf.@sprintf("  %-38s %12.1f %9.1f%% %9d", label, secs, pct, n))
    end
    println(io, "  " * "-"^72)
    println(io, Printf.@sprintf("  %-38s %12.1f %9.1f%% %9d", "TOTAL (attributed)", wall_seconds, 100.0, total))
end

# ---------------------------------------------------------------------------
# Method 2 — direct @elapsed of each component for ONE representative iteration.
# Independent of the sampler; calls the real functions on the real graph/reads.
# ---------------------------------------------------------------------------
function method2_single_iteration(reads)
    graph_mode = :doublestrand
    # initial k the corrector would pick, clamped into the corrector's floor/ceiling.
    k0 = try
        Mycelia.find_initial_k(reads)
    catch
        min(MAX_K, 13)
    end
    k = min(max(k0, 13), MAX_K)

    results = Tuple{String, Float64, String}[]  # (component, seconds, note)

    local graph
    t = @elapsed (graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = graph_mode))
    nverts = length(collect(Mycelia.MetaGraphsNext.labels(graph)))
    nedges = length(collect(Mycelia.MetaGraphsNext.edge_labels(graph)))
    push!(results, ("1. build_qualmer_graph", t, "$(nverts) vertices, $(nedges) edges, k=$k"))

    local hard_vertices
    t = @elapsed (hard_vertices = Mycelia._hard_vertex_set(graph, k))
    push!(results, ("2. hard_vertex_set (bubbles)", t, "$(length(hard_vertices)) hard vertices"))

    local solid
    t = @elapsed (solid = Mycelia._solid_kmer_set(graph))
    push!(results, ("3. solid_kmer classification", t, "$(length(solid)) solid k-mers"))

    local cheap_out
    t = @elapsed (cheap_out = Mycelia._stage0_cheap_correct(reads, k, solid; graph_mode = graph_mode))
    work_reads = cheap_out[1]
    push!(results, ("4. stage0 cheap_correct", t, "$(cheap_out[2]) bases fixed over $(length(reads)) reads"))

    # Decode ONLY the reads the hard-window gate would decode (matches production).
    diag = Mycelia.CorrectorDiagnostics()
    decoded = 0
    t = @elapsed begin
        for r in work_reads
            if Mycelia.should_decode_read(r, k, hard_vertices; graph_mode = graph_mode)
                Mycelia.improve_read_likelihood(r, graph, k; graph_mode = graph_mode,
                    beam_width = nothing, diagnostics = diag)
                decoded += 1
            end
        end
    end
    push!(results, ("5. per-read Viterbi decode", t,
        "$(decoded)/$(length(work_reads)) reads decoded (rest skipped by gate)"))

    dir = mktempdir()
    outp = joinpath(dir, "m2_out.fastq")
    t = @elapsed Mycelia.write_fastq(records = work_reads, filename = outp)
    push!(results, ("6. write_fastq", t, "$(length(work_reads)) reads"))

    meta = (k = k, n_reads = length(reads), n_vertices = nverts, n_edges = nedges,
        n_decoded = decoded, n_work_reads = length(work_reads))
    return k, results, meta
end

function print_method2(io, k, results)
    total = sum(r[2] for r in results)
    println(io, "\nMethod 2 — direct @elapsed, ONE representative iteration (initial k=$k)")
    println(io, "  (single-iteration cross-check; production runs this ~k_rungs x iters times)")
    println(io, "")
    println(io, Printf.@sprintf("  %-38s %12s %10s  %s", "component", "seconds", "% of iter", "note"))
    println(io, "  " * "-"^96)
    for (name, secs, note) in results
        pct = total == 0 ? 0.0 : 100 * secs / total
        println(io, Printf.@sprintf("  %-38s %12.3f %9.1f%%  %s", name, secs, pct, note))
    end
    println(io, "  " * "-"^96)
    println(io, Printf.@sprintf("  %-38s %12.3f %9.1f%%", "TOTAL (single iteration)", total, 100.0))
end

# ---------------------------------------------------------------------------
# FINDINGS (2026-07-07 sweep, 1/2/3/5 kb, err=0.01, 20x, 100bp reads; see PR):
#
#   component                        1kb      2kb      3kb      5kb    alpha  verdict
#   1. build_qualmer_graph          0.06s    0.19s    0.28s    0.39s   1.15   ~linear
#   2. hard_vertex_set (bubbles)    1.42s    6.11s   13.81s   47.29s   2.16   O(V^2)!!
#   3. solid_kmer classification    0.01s    0.01s    0.02s    0.04s   1.00   linear
#   4. stage0 cheap_correct         0.01s    0.02s    0.03s    0.04s   1.02   linear
#   5. per-read Viterbi decode      6.69s   26.63s   63.77s  210.26s   2.14   O(V^2)!!  <-- DOMINANT
#   6. write_fastq                  0.00s    0.00s    0.00s    0.00s     -    negligible
#   TOTAL (single pass)             8.19s   32.96s   77.90s  258.03s   2.14
#   (n_reads = 200/400/600/1000; n_vertices = 8278/16320/24422/39912 — both ~linear in genome)
#
# TWO super-linear culprits, BOTH quadratic (alpha ~2.1), SAME anti-pattern —
# a per-item step that internally does a full O(V+E) graph scan:
#
#   (5) DECODE (dominant, 81% of a pass): the bounded beam DOES bound the Viterbi
#   DP (auto-beam is sized by observation-count ~= read length, a constant 100bp
#   at every genome size). The quadratic is NOT the DP — it is that
#   `correct_observations` -> `_correct_metagraphs_next_observations` calls
#   `weighted_graph_from_rhizomorph(graph)` (src/rhizomorph/algorithms/
#   path-finding.jl:324), which REBUILDS a whole-graph weighted COPY in O(V+E),
#   ONCE PER READ (try_viterbi_path_improvement decodes each read separately).
#   n_reads ∝ genome and V ∝ genome, so per-read O(V) x O(V) reads = O(genome^2).
#   FIX: build the StrandWeightedEdgeData graph ONCE per pass and reuse it across
#   all reads (hoist the weighting out of the per-read loop), or have
#   build_qualmer_graph emit StrandWeightedEdgeData directly so the `weighted =
#   graph` fast-path (path-finding.jl) is taken and NO per-read copy occurs.
#   Expected: decode drops from O(genome^2) to O(genome) (linear in n_reads).
#
#   (2) HARD_VERTEX_SET / detect_bubbles_next: `detect_bubbles_next`
#   (src/rhizomorph/algorithms/simplification.jl:90) loops over ALL V vertices and
#   per vertex calls `get_out_neighbors`, which scans ALL E edges -> O(V*E) =
#   O(V^2). FIX: build an out-adjacency index once in O(E) (Dict{vertex =>
#   Vector{neighbor}}) or use native Graphs.outneighbors, reducing to O(V+E).
#   (`_hard_vertex_set` already builds an outdeg Dict in O(E) for its high-degree
#   pass — the same index can feed bubble detection.)
#
# NOT culprits: build_qualmer_graph, solid_kmer, stage0, write_fastq are all
# ~linear; the iteration count is fixed at <=6 passes (not size-driven).
# ---------------------------------------------------------------------------
# SCALING SWEEP (td-9q84 follow-up): measure how each per-iteration component's
# per-PASS cost grows with genome size, to localize the super-linear term.
#
# WHY single-pass (not full runs): the :scalable knobs cap iterations at
# n_k_rungs(3) x max_iterations_per_k(2) = 6 passes REGARDLESS of genome size, so
# the total blow-up is NOT from a growing iteration count — it must be a
# per-pass component whose cost scales super-linearly with graph size. Method 2
# times each component ONCE on the real first-rung graph/reads, so it isolates
# per-pass scaling directly and (critically) works at 5kb WITHOUT running the
# full assembly to completion (which does not finish in 3h).
#
# The 6-pass cap is verified structurally from _corrector_strategy_knobs(:scalable)
# (n_k_rungs=3, max_iterations_per_k=2) and printed in the report so the reader
# can confirm iteration count is bounded, not the culprit.
# ---------------------------------------------------------------------------

# Component keys, in loop order, shared by the per-size rows and the fit table.
const SCALING_COMPONENTS = [
    "1. build_qualmer_graph",
    "2. hard_vertex_set (bubbles)",
    "3. solid_kmer classification",
    "4. stage0 cheap_correct",
    "5. per-read Viterbi decode",
    "6. write_fastq",
]

# One scaling row = first-pass per-component times + graph/read sizes at a size.
function scaling_row(genome_len::Int)
    _refseq, reads = build_fixture_at(genome_len)
    k, results, meta = method2_single_iteration(reads)
    times = Dict(name => secs for (name, secs, _note) in results)
    total = sum(secs for (_n, secs, _note) in results)
    return (genome_len = genome_len, k = k, n_reads = meta.n_reads,
        n_vertices = meta.n_vertices, n_edges = meta.n_edges,
        n_decoded = meta.n_decoded, n_work_reads = meta.n_work_reads,
        times = times, total = total)
end

# Ordinary-least-squares slope of log(y) on log(x): the empirical scaling
# exponent alpha in y ~ x^alpha. Needs >=2 finite, positive (x,y) pairs; returns
# NaN otherwise (e.g. a component too fast to time reliably -> 0.0 seconds).
function loglog_alpha(xs::Vector{<:Real}, ys::Vector{<:Real})
    lx = Float64[]
    ly = Float64[]
    for (x, y) in zip(xs, ys)
        if x > 0 && y > 0 && isfinite(x) && isfinite(y)
            push!(lx, log(x))
            push!(ly, log(y))
        end
    end
    length(lx) < 2 && return NaN
    mx = sum(lx) / length(lx)
    my = sum(ly) / length(ly)
    num = sum((lx[i] - mx) * (ly[i] - my) for i in eachindex(lx))
    den = sum((lx[i] - mx)^2 for i in eachindex(lx))
    den == 0 && return NaN
    return num / den
end

function print_scaling(io, rows, fit_sizes::Vector{Int})
    knobs = Mycelia.Rhizomorph._corrector_strategy_knobs(:scalable)
    max_passes = knobs.n_k_rungs * knobs.max_iterations_per_k
    println(io, "\n" * "="^96)
    println(io, "SCALING SWEEP — first-pass per-component cost vs genome size (td-9q84)")
    println(io, "  method: Method 2 (single first-rung pass) at each size; behavior UNCHANGED.")
    println(io, "  iteration budget (fixed, size-independent): n_k_rungs=$(knobs.n_k_rungs) x " *
                "max_iterations_per_k=$(knobs.max_iterations_per_k) = $(max_passes) passes max.")
    println(io, "  => total runtime = (per-pass cost) x (<= $(max_passes) passes); a super-linear")
    println(io, "     per-pass component is therefore the whole-run super-linear culprit.")
    println(io, "="^96)

    # Graph/read size table
    println(io, "\nGraph & read sizes per genome length:")
    println(io, Printf.@sprintf("  %8s %8s %8s %10s %10s %10s %10s",
        "genome", "k", "n_reads", "n_verts", "n_edges", "decoded", "work_rds"))
    println(io, "  " * "-"^74)
    for r in rows
        println(io, Printf.@sprintf("  %8d %8d %8d %10d %10d %10d %10d",
            r.genome_len, r.k, r.n_reads, r.n_vertices, r.n_edges, r.n_decoded, r.n_work_reads))
    end

    # Per-component seconds at each size
    sizes = [r.genome_len for r in rows]
    println(io, "\nPer-component FIRST-PASS seconds (columns = genome length in bp):")
    hdr = Printf.@sprintf("  %-32s", "component")
    for s in sizes
        hdr *= Printf.@sprintf(" %12s", "$(s)bp")
    end
    println(io, hdr)
    println(io, "  " * "-"^(34 + 13 * length(sizes)))
    for comp in SCALING_COMPONENTS
        line = Printf.@sprintf("  %-32s", comp)
        for r in rows
            line *= Printf.@sprintf(" %12.4f", get(r.times, comp, NaN))
        end
        println(io, line)
    end
    total_line = Printf.@sprintf("  %-32s", "TOTAL (single pass)")
    for r in rows
        total_line *= Printf.@sprintf(" %12.4f", r.total)
    end
    println(io, "  " * "-"^(34 + 13 * length(sizes)))
    println(io, total_line)

    # Scaling-exponent fits. alpha vs genome length AND vs n_vertices, using only
    # the `fit_sizes` rows (the sizes that complete cleanly; 5kb is shown but is a
    # single first-pass probe, included in the fit too since Method 2 never runs
    # it to completion).
    fit_rows = [r for r in rows if r.genome_len in fit_sizes]
    glens = Float64[r.genome_len for r in fit_rows]
    nverts = Float64[r.n_vertices for r in fit_rows]
    println(io, "\nEmpirical scaling exponent alpha  (time ~ x^alpha), fit over sizes " *
                "$(sort(collect(fit_sizes))):")
    println(io, Printf.@sprintf("  %-32s %14s %14s %10s",
        "component", "alpha(genome)", "alpha(n_verts)", "verdict"))
    println(io, "  " * "-"^74)
    for comp in vcat(SCALING_COMPONENTS, ["TOTAL (single pass)"])
        ys = comp == "TOTAL (single pass)" ?
            Float64[r.total for r in fit_rows] :
            Float64[get(r.times, comp, NaN) for r in fit_rows]
        a_g = loglog_alpha(glens, ys)
        a_v = loglog_alpha(nverts, ys)
        # Verdict keys off the genome-length exponent (falls back to n_verts).
        a_ref = isnan(a_g) ? a_v : a_g
        verdict = isnan(a_ref) ? "n/a (too fast)" :
            a_ref > 1.7 ? "SUPER-LINEAR!!" :
            a_ref > 1.3 ? "super-linear" :
            a_ref > 1.1 ? "~linear+" :
            "~linear/sub"
        gs = isnan(a_g) ? "   n/a" : Printf.@sprintf("%14.2f", a_g)
        vs = isnan(a_v) ? "   n/a" : Printf.@sprintf("%14.2f", a_v)
        println(io, Printf.@sprintf("  %-32s %14s %14s %10s", comp, gs, vs, verdict))
    end
    println(io, "  " * "-"^74)
    println(io, "\nReading the table: alpha ~ 1.0 => linear (fine). alpha > ~1.3 => super-linear")
    println(io, "  (the blow-up term). alpha ~ 2.0 => quadratic in that size axis.")
end

function main_scaling()
    io = stdout
    sizes_str = get(ENV, "MYCELIA_SCP_SIZES", "1000,2000,3000,5000")
    sizes = [parse(Int, strip(s)) for s in split(sizes_str, ",") if !isempty(strip(s))]
    # Sizes to include in the alpha fit (default: all of them — every size is a
    # single first-pass Method 2 probe, so even 5kb is cheap and fittable).
    fit_str = get(ENV, "MYCELIA_SCP_FIT_SIZES", sizes_str)
    fit_sizes = [parse(Int, strip(s)) for s in split(fit_str, ",") if !isempty(strip(s))]

    println(io, "="^96)
    println(io, ":scalable corrector SCALING sweep (td-9q84)")
    println(io, "  started: $(Dates.now())")
    println(io, "  sizes: $(sizes) bp | $(COVERAGE)x | $(READLEN)bp reads | err=$(ERR_RATE) | " *
                "tech=$(TECH) | max_k=$(MAX_K)")
    println(io, "="^96)

    # Warmup: compile every timed component once on a tiny fixture so the first
    # measured size is not polluted by first-call compilation time.
    println(io, "\n[warmup] compiling timed components on a tiny fixture ...")
    _wref, wreads = build_fixture_at(300)
    method2_single_iteration(wreads)
    println(io, "[warmup] done")

    rows = NamedTuple[]
    for gl in sizes
        println(io, "\n[size $(gl)bp] building fixture + timing one first-rung pass ...")
        t0 = time()
        row = scaling_row(gl)
        println(io, Printf.@sprintf("[size %dbp] done in %.1fs  (n_reads=%d, n_verts=%d, n_edges=%d)",
            gl, time() - t0, row.n_reads, row.n_vertices, row.n_edges))
        push!(rows, row)
    end

    print_scaling(io, rows, fit_sizes)

    println(io, "\n" * "="^96)
    println(io, "DONE. The 'verdict' column flags components with alpha > ~1.3 as the")
    println(io, "super-linear culprit(s). Total runtime = per-pass cost x <=6 fixed passes,")
    println(io, "so a super-linear per-pass component drives the whole-run blow-up.")
    println(io, "  finished: $(Dates.now())")
    println(io, "="^96)
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    io = stdout
    println(io, "="^80)
    println(io, ":scalable corrector per-iteration profile (td-9q84)")
    println(io, "  started: $(Dates.now())")
    println(io, "  fixture: $(GENOME_LEN)bp genome, $(COVERAGE)x, $(READLEN)bp reads, " *
                "err=$(ERR_RATE), tech=$(TECH), max_k=$(MAX_K)")
    println(io, "="^80)

    refseq, reads = build_fixture()
    println(io, "simulated $(length(reads)) reads over a $(length(refseq))bp reference")

    knobs = Mycelia.Rhizomorph._corrector_strategy_knobs(:scalable)
    println(io, "scalable knobs: $(knobs)")

    fastq_path = write_fixture_fastq(reads)

    # --- Warmup: force compilation of every corrector code path on a TINY input so
    # the profiled run's samples are ~pure runtime, not first-call compilation. ---
    println(io, "\n[warmup] compiling corrector paths on a tiny input ...")
    warm_rng = Random.MersenneTwister(SEED + 1)
    warm_ref = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED + 1, L = 300)
    warm_seq = FASTX.sequence(BioSequences.LongDNA{4}, warm_ref)
    warm_reads = simulate_reads(warm_seq, READLEN, 10, ERR_RATE, TECH, warm_rng)
    warm_path = write_fixture_fastq(warm_reads)
    warm_out = mktempdir()
    Mycelia.mycelia_iterative_assemble(warm_path;
        max_k = max(min(MAX_K, 15), 13),
        skip_solid = knobs.skip_solid, graph_mode = knobs.graph_mode,
        n_k_rungs = knobs.n_k_rungs, max_iterations_per_k = knobs.max_iterations_per_k,
        hard_window = knobs.hard_window, soft_em = knobs.soft_em,
        cheap_correct = knobs.cheap_correct, beam_width = knobs.beam_width,
        verbose = false, enable_checkpointing = false, output_dir = warm_out)
    println(io, "[warmup] done")

    # --- Method 1: profile the real per-iteration loop on the 1kb fixture. ---
    println(io, "\n[method 1] profiling mycelia_iterative_assemble (:scalable) on the fixture ...")
    Profile.clear()
    Profile.init(n = 10^8, delay = PROFILE_DELAY)
    prof_out = mktempdir()
    wall = @elapsed Profile.@profile Mycelia.mycelia_iterative_assemble(fastq_path;
        max_k = MAX_K,
        skip_solid = knobs.skip_solid, graph_mode = knobs.graph_mode,
        n_k_rungs = knobs.n_k_rungs, max_iterations_per_k = knobs.max_iterations_per_k,
        hard_window = knobs.hard_window, soft_em = knobs.soft_em,
        cheap_correct = knobs.cheap_correct, beam_width = knobs.beam_width,
        verbose = false, enable_checkpointing = false, output_dir = prof_out)

    counts, total = bucket_profile(wall)
    print_table(io, counts, total, wall)

    # --- Method 2: single-iteration direct timing cross-check. ---
    println(io, "\n[method 2] direct component timing, one representative iteration ...")
    k, m2, _meta = method2_single_iteration(reads)
    print_method2(io, k, m2)

    println(io, "\n" * "="^80)
    println(io, "DOMINANT TERM: read the '% runtime' column above (Method 1 is authoritative;")
    println(io, "Method 2 corroborates on a single iteration). See the harness header for the")
    println(io, "loop structure each component sits in.")
    println(io, "  finished: $(Dates.now())")
    println(io, "="^80)
end

# Mode dispatch (td-9q84 scaling follow-up):
#   MYCELIA_SCP_MODE=scaling (default) -> multi-size per-component scaling sweep,
#                                          emits the SCALING TABLE + alpha fits.
#   MYCELIA_SCP_MODE=single            -> original single-size Method1+Method2 run.
const RUN_MODE = get(ENV, "MYCELIA_SCP_MODE", "scaling")
if RUN_MODE == "single"
    main()
else
    main_scaling()
end
