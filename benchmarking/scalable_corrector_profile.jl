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
    push!(results, ("1. build_qualmer_graph", t, "$(nverts) vertices, k=$k"))

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

    return k, results
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
    k, m2 = method2_single_iteration(reads)
    print_method2(io, k, m2)

    println(io, "\n" * "="^80)
    println(io, "DOMINANT TERM: read the '% runtime' column above (Method 1 is authoritative;")
    println(io, "Method 2 corroborates on a single iteration). See the harness header for the")
    println(io, "loop structure each component sits in.")
    println(io, "  finished: $(Dates.now())")
    println(io, "="^80)
end

main()
