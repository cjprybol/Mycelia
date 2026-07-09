"""
Iterative Maximum Likelihood Assembly - Phase 5.2a Implementation

This module implements the iterative maximum likelihood assembly system with:
- Complete FASTQ I/O processing per iteration
- Statistical path resampling with likelihood calculations
- Viterbi algorithm integration for optimal path finding
- Timestamped output files for tracking read evolution
- Memory-efficient read set processing

Part of the Mycelia bioinformatics package's iterative assembler framework.
"""

# =============================================================================
# Sparsity Detection and Memory Estimation
# =============================================================================

"""
Calculate k-mer sparsity for a given k-mer size.
Returns fraction of possible k-mers that are NOT observed.
"""
function calculate_sparsity(reads::Vector{<:FASTX.FASTQ.Record}, k::Int)::Float64
    observed_kmers = Set{String}()

    for record in reads
        seq = FASTX.sequence(String, record)
        if length(seq) >= k
            for i in 1:(length(seq) - k + 1)
                kmer = seq[i:(i + k - 1)]
                push!(observed_kmers, kmer)
            end
        end
    end

    unique_observed = length(observed_kmers)
    max_possible = 4.0^k
    return 1.0 - (unique_observed / max_possible)
end

"""
Analyze k-mer coverage distribution to detect if errors are singletons.
Returns true if low-coverage k-mers (likely errors) are well-separated from high-coverage ones.
"""
function errors_are_singletons(reads::Vector{<:FASTX.FASTQ.Record}, k::Int; singleton_threshold::Int = 2)::Bool
    kmer_counts = Dict{String, Int}()

    for record in reads
        seq = FASTX.sequence(String, record)
        if length(seq) >= k
            for i in 1:(length(seq) - k + 1)
                kmer = seq[i:(i + k - 1)]
                kmer_counts[kmer] = get(kmer_counts, kmer, 0) + 1
            end
        end
    end

    coverage_values = collect(values(kmer_counts))
    total_unique = length(coverage_values)
    if total_unique == 0
        return false
    end
    singleton_count = count(x -> x <= singleton_threshold, coverage_values)
    singleton_fraction = singleton_count / total_unique

    if singleton_count > 0 && total_unique > singleton_count
        non_singleton_min = minimum(filter(x -> x > singleton_threshold, coverage_values))
        return singleton_fraction > 0.1 && non_singleton_min > singleton_threshold * 2
    end

    return false
end

"""
Find the optimal starting k-mer size using sparsity detection.
Only considers prime k-mer sizes for optimal performance.
"""
function find_initial_k(reads::Vector{<:FASTX.FASTQ.Record};
        k_range::Vector{Int} = Primes.primes(3, 51),
        sparsity_threshold::Float64 = 0.5)::Int
    for k in k_range
        sparsity = calculate_sparsity(reads, k)
        if sparsity > sparsity_threshold && errors_are_singletons(reads, k)
            return k
        end
    end

    return isempty(k_range) ? 3 : first(k_range)
end

"""
Find the next prime number greater than current_k.
For k-mer progression, we prefer odd numbers and especially primes.
"""
function next_prime_k(current_k::Int; max_k::Int = 1000)::Int
    next_prime = Primes.nextprime(current_k + 1)
    return next_prime <= max_k ? next_prime : current_k
end

"""
    build_k_ladder(initial_k, max_k; k_ladder=nothing, n_k_rungs=nothing)

Compute the ordered set of k-mer sizes the iterative corrector should walk.

- If `k_ladder` is given, use those values that fall in `[initial_k, max_k]`
  (deduplicated and sorted). NOTE: the assembly loop always STARTS at
  `find_initial_k` regardless of the ladder, so if `min(k_ladder) > initial_k` the
  auto-detected initial k is processed first and then the ladder is followed — the
  ladder specifies the k values to visit AFTER the initial k, not a replacement
  for it. (In `n_k_rungs` mode the first rung is pinned to `initial_k`, so they
  coincide.)
- Else if `n_k_rungs` is given, build an ~`n_k_rungs`-rung geometric ladder from
  `initial_k` to `max_k` (LoRMA-style coarse progression), snapping intermediate
  rungs to odd values and pinning the first rung to `initial_k` and the last to
  the largest odd `<= max_k`.
- Else return `nothing`, signalling the caller to use the original prime-by-prime
  progression (`next_prime_k`) — this preserves legacy behavior.

Returns a sorted `Vector{Int}` (or `nothing`).
"""
function build_k_ladder(initial_k::Int, max_k::Int;
        k_ladder::Union{Nothing, Vector{Int}} = nothing,
        n_k_rungs::Union{Nothing, Int} = nothing)::Union{Nothing, Vector{Int}}
    if k_ladder !== nothing
        rungs = sort(unique(filter(k -> initial_k <= k <= max_k, k_ladder)))
        return isempty(rungs) ? [initial_k] : rungs
    end
    if n_k_rungs === nothing
        return nothing
    end
    if initial_k >= max_k
        return [initial_k]
    end
    n = max(2, n_k_rungs)
    top = isodd(max_k) ? max_k : max_k - 1
    if top <= initial_k
        return [initial_k]
    end
    ratio = (top / initial_k)^(1 / (n - 1))
    ks = Int[]
    for i in 1:n
        kv = round(Int, initial_k * ratio^(i - 1))
        iseven(kv) && (kv += 1)          # prefer odd k
        kv = clamp(kv, initial_k, top)
        push!(ks, kv)
    end
    ks[1] = initial_k
    ks[end] = top
    return sort(unique(ks))
end

"""
    _next_k_in_progression(current_k, max_k, k_schedule)

Return the next k-mer size after `current_k`. When `k_schedule === nothing`, use
the original prime progression (`next_prime_k`). Otherwise advance to the next
scheduled rung strictly greater than `current_k`, or return `current_k` (a
fixed point, which the main loop treats as "stop") when none remain.
"""
function _next_k_in_progression(current_k::Int, max_k::Int,
        k_schedule::Union{Nothing, Vector{Int}})::Int
    if k_schedule === nothing
        return next_prime_k(current_k; max_k = max_k)
    end
    larger = filter(>(current_k), k_schedule)
    return isempty(larger) ? current_k : minimum(larger)
end

# Default ADAPTIVE low-k decode-gate density threshold (td-9h5r): if, at a given
# k-rung, the post-Stage-0 hard-window gate would still send AT LEAST this fraction
# of reads to the expensive per-read graph-Viterbi decode, the gate is judged
# NON-DISCRIMINATING (it is not reducing decode volume) and the whole decode is
# skipped for that pass — the reads fall back on Stage 0 cheap correction +
# skip-solid, and the whole-read Viterbi is deferred to higher, selective rungs
# where the gate actually filters. 0.90 leaves a genuinely selective low-k decode
# (e.g. a gate that skips ~40-60% of reads) ON, and only gates the dense-graph
# "decode almost everything" pathology the #370 profile flagged.
const _DEFAULT_DECODE_GATE_DENSITY = 0.90

"""
Estimate memory usage for a graph with a given number of k-mers.
Provides a rough estimate for memory monitoring.
"""
function estimate_memory_usage(num_kmers::Int, k::Int)::Int
    kmer_size = k * 1
    vertex_overhead = 100
    edge_overhead = 50

    estimated_bytes = num_kmers * (kmer_size + vertex_overhead + edge_overhead * 2)
    return round(Int, estimated_bytes * 1.2)
end

function _label_k_length(label)::Int
    if label isa Qualmer
        return length(label)
    elseif label isa Kmers.Kmer
        return length(label)
    elseif label isa BioSequences.BioSequence
        return length(label)
    elseif label isa AbstractString
        return length(label)
    end
    return length(string(label))
end

"""
Check if memory usage is within acceptable limits.
"""
function check_memory_limits(graph, memory_limit::Int)::Bool
    labels = if applicable(MetaGraphsNext.labels, graph)
        collect(MetaGraphsNext.labels(graph))
    else
        Any[]
    end

    num_kmers = length(labels)
    if num_kmers == 0
        return true
    end
    k = _label_k_length(first(labels))
    estimated_usage = estimate_memory_usage(num_kmers, k)
    return estimated_usage <= memory_limit
end

# =============================================================================
# Core Iterative Assembly Framework
# =============================================================================

"""
Main iterative maximum likelihood assembly function.
Processes entire read sets per iteration with complete FASTQ I/O tracking.
Enhanced with performance optimizations, caching, and progress tracking.
"""
function mycelia_iterative_assemble(input_fastq::String;
        max_k::Int = 101,
        memory_limit::Int = 32_000_000_000,
        output_dir::String = "iterative_assembly",
        max_iterations_per_k::Int = 10,
        improvement_threshold::Float64 = 0.05,
        stop_on_no_change::Bool = true,
        k_ladder::Union{Nothing, Vector{Int}} = nothing,
        n_k_rungs::Union{Nothing, Int} = nothing,
        graph_mode::Symbol = :canonical,
        verbose::Bool = true,
        enable_parallel::Bool = false,
        batch_size::Int = 10000,
        enable_checkpointing::Bool = true,
        checkpoint_interval::Int = 5,
        skip_solid::Bool = false,
        hard_window::Bool = false,
        soft_em::Bool = false,
        cheap_correct::Bool = false,
        beam_width::Union{Int, Nothing} = nothing,
        min_decode_k::Union{Int, Nothing} = nothing,
        decode_gate_density::Union{Float64, Nothing} = nothing)
    start_time = time()

    # Accumulates swallowed decode failures (structural / un-k-merizable) across
    # every k and iteration so a systematically-broken corrector that reports
    # "0 improvements" is visible in the result metadata rather than passing as
    # "nothing to fix". See CorrectorDiagnostics.
    corrector_diagnostics = CorrectorDiagnostics()

    # -- Soft-EM registry hygiene (td-e70t, v2 competing-paths + support floor) -
    # Soft-EM v2 ACTIVATES the M-step: within each k's EM loop, iteration N's
    # per-edge path responsibilities are REGISTERED onto iteration N+1's graph
    # (`register_soft_edge_weights!`, floored to each edge's own raw support so
    # supported variation is retained — td-h6w9) so `compute_edge_weight` returns
    # the probability-weighted (soft) evidence and unsupported error edges decay.
    # Registration is ALWAYS paired with `clear!` in a try/finally INSIDE the loop,
    # so the process-global registry is EMPTY at corrector entry (and between
    # iterations). Defensively clear it here so no prior crash mid-registration (a
    # direct unit-test primitive call, an aborted run) can leak soft weights into
    # this correction, then assert the invariant (belt-and-suspenders).
    Mycelia.Rhizomorph.clear_soft_edge_weights!()
    @assert isempty(Mycelia.Rhizomorph._SOFT_EDGE_WEIGHT_REGISTRY) (
        "soft-EM registry must be empty at corrector entry (registration is scoped " *
        "to each EM iteration and always cleared in a try/finally)")

    # -- Convergence + k-ladder tuning knobs (td-q70n) -------------------------
    #
    # Two literature-backed speedups for the iterative corrector:
    #
    #   1. CONVERGENCE (Musket): `stop_on_no_change=true` breaks the per-k loop the
    #      moment a pass makes 0 changes, independent of `improvement_threshold`
    #      (the old code only bounded the loop when improvement < threshold, so it
    #      never terminated on a 0-change pass when improvement_threshold==0).
    #      `max_iterations_per_k` KEEPS its default of 10 — lowering it to 2
    #      (Musket's ~2-pass heuristic) is a real accuracy/speed tradeoff that must
    #      be measured on error-rich data before becoming a default, so it is
    #      opt-in; the corrector=:iterative route in assemble_genome sets it
    #      explicitly. Note the corrector is not strictly deterministic (the
    #      statistical-resampling arm draws from the global RNG), so a 0-change
    #      pass is a practical, not guaranteed, fixed point.
    #
    #   2. K-LADDER (LoRMA 3-rung): instead of walking every prime
    #      (3,5,7,11,13,...) the caller may request a small, well-spaced ladder.
    #      `k_ladder=[...]` uses exactly those k values; `n_k_rungs=N` builds an
    #      ~N-rung geometric ladder from the initial k to `max_k`. Both default
    #      to `nothing`, which preserves the original prime-by-prime progression
    #      byte-for-byte, so existing callers are unchanged unless they opt in.

    if verbose
        println("Starting Mycelia Iterative Maximum Likelihood Assembly")
        println("Input FASTQ: $input_fastq")
        println("Output directory: $output_dir")
        println("Memory limit: $(memory_limit ÷ 1_000_000_000) GB")
        println("Max k-mer size: $max_k")
        println("Graph mode: $graph_mode")
        println("Parallel processing: $(enable_parallel ? "enabled ($(Threads.nthreads()) threads)" : "disabled")")
        println("Batch size: $batch_size")
        println("Checkpointing: $(enable_checkpointing ? "enabled (every $checkpoint_interval iterations)" : "disabled")")
    end

    # Create output directory structure
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Create subdirectories for organization
    checkpoints_dir = joinpath(output_dir, "checkpoints")
    graphs_dir = joinpath(output_dir, "graphs")
    progress_dir = joinpath(output_dir, "progress")

    if enable_checkpointing
        mkpath(checkpoints_dir)
        mkpath(graphs_dir)
        mkpath(progress_dir)
    end

    # Check for existing checkpoint to resume from
    checkpoint_file = joinpath(checkpoints_dir, "latest_checkpoint.json")
    resume_data = nothing

    if enable_checkpointing && isfile(checkpoint_file)
        try
            resume_data = JSON.parsefile(checkpoint_file)
            if verbose
                println("Found existing checkpoint. Resume from k=$(resume_data["current_k"]), iteration=$(resume_data["current_iteration"])")
            end
        catch e
            if verbose
                println("Warning: Could not load checkpoint file: $e")
            end
        end
    end

    # Initialize or resume from checkpoint
    if resume_data !== nothing
        # Resume from checkpoint
        k = resume_data["current_k"]
        k_progression = resume_data["k_progression"]
        iteration_history = Dict{Int, Vector{Dict{Symbol, Any}}}()
        for (k_str, hist) in resume_data["iteration_history"]
            iteration_history[parse(Int, k_str)] = hist
        end
        total_improvements = resume_data["total_improvements"]
        current_fastq_file = resume_data["current_fastq_file"]

        if verbose
            println("Resumed from checkpoint: k=$k, total_improvements=$total_improvements")
        end
    else
        # Initialize fresh run
        if verbose
            println("Reading initial FASTQ file...")
        end
        initial_reads = collect(FASTX.FASTQ.Reader(open(input_fastq)))
        k = find_initial_k(initial_reads)  # Reuse from intelligent-assembly.jl

        if verbose
            println("Initial k-mer size: $k (prime: $(Primes.isprime(k)))")
            println("Total reads: $(length(initial_reads))")
        end

        # Track progress across all iterations
        k_progression = Int[]
        iteration_history = Dict{Int, Vector{Dict{Symbol, Any}}}()
        total_improvements = 0
        current_fastq_file = input_fastq
    end

    # Build the k-mer progression schedule. `nothing` => legacy prime-by-prime
    # walk via next_prime_k; a Vector{Int} => explicit / coarse LoRMA-style ladder.
    k_schedule = build_k_ladder(k, max_k; k_ladder = k_ladder, n_k_rungs = n_k_rungs)
    if verbose && k_schedule !== nothing
        println("K-mer ladder (coarse progression): $(k_schedule)")
    end

    # -- Low-k decode gating (td-9h5r) -----------------------------------------
    # The finding (PR #370 profile): the per-read graph-Viterbi decode is 57-78% of
    # every :scalable iteration — the dominant runtime term — yet at LOW k the graph
    # is dense (nearly every k-mer sits on a bubble/repeat vertex), so the
    # hard-window gate cannot discriminate and flags essentially EVERY read for a
    # full whole-read decode. That is correction Stage 0's cheap linear pass already
    # largely did. Two composable gates cut that wasted volume; BOTH still BUILD the
    # graph and run Stage 0 cheap correction + skip-solid, deferring only the
    # expensive whole-read Viterbi to higher, selective rungs:
    #
    #   (1) ADAPTIVE (default): at a rung whose post-Stage-0 hard-window decode
    #       fraction would be >= `decode_gate_density`, the gate is non-
    #       discriminating (it would decode ~everything → no volume reduction), so
    #       the whole decode is skipped for that pass. This fires EXACTLY on the
    #       dense-low-k pathology and is a NO-OP where the gate genuinely filters
    #       (so a load-bearing selective low-k decode is preserved). Measured
    #       post-Stage-0 inside improve_read_set_likelihood (raw-read density is
    #       uninformative — Stage 0 clears most apparent hardness first).
    #   (2) EXPLICIT floor: `min_decode_k` hard-gates every rung `k < min_decode_k`
    #       (the LoRMA-style "start the graph-decode ladder higher" lever). Off by
    #       default; deterministic, for tuning/tests.
    #
    # Both are gated behind the `:scalable` tier (only active when `hard_window` is
    # true), so the `:exhaustive` tier (hard_window=false) is BYTE-IDENTICAL.
    effective_min_decode_k = hard_window ? min_decode_k : nothing
    # Adaptive-gate density threshold (fraction of reads the hard-window gate would
    # still decode after Stage 0 above which the gate counts as non-discriminating).
    # Default 0.90 on :scalable; `nothing` (or the exhaustive tier) disables it.
    effective_decode_gate_density = hard_window ?
                                    (decode_gate_density === nothing ?
                                     _DEFAULT_DECODE_GATE_DENSITY : decode_gate_density) :
                                    nothing
    if verbose && (effective_min_decode_k !== nothing ||
                   effective_decode_gate_density !== nothing)
        println("Low-k decode gating (td-9h5r): " *
                (effective_min_decode_k !== nothing ?
                 "explicit floor min_decode_k=$(effective_min_decode_k); " : "") *
                (effective_decode_gate_density !== nothing ?
                 "adaptive gate-off when post-Stage-0 decode fraction >= " *
                 "$(effective_decode_gate_density). " : "") *
                "Gated rungs rely on Stage 0 cheap correction + skip-solid.")
    end
    # Telemetry: the k-rungs whose per-read decode was gated OFF (low-k skip).
    decode_gated_rungs = Int[]

    # Hard-window skip telemetry (td-nn6l): the fraction of reads the hard-window
    # gate passed through WITHOUT a decode, recorded PER PASS (across every k and
    # iteration) so the run metadata can surface min/mean/max, not just the last
    # value (FIX 6). Empty on the :exhaustive tier (no gate ⇒ no skips).
    skip_fractions = Float64[]

    # Stage 0 cheap-correction telemetry (td-bjnt): bases fixed by the linear
    # k-mer-spectrum pass, recorded PER PASS. Empty/zero when cheap_correct is off
    # (:exhaustive), so the exhaustive tier is unaffected.
    cheap_correction_counts = Int[]

    # --- Final-pass graph reuse (td-04tb) -----------------------------------
    # The corrector builds a qualmer graph FROM SCRATCH each pass; the LAST pass at
    # the largest k builds the graph whose input reads become the final corrected
    # read set — but ONLY when that pass makes 0 improvements (a converged, no-change
    # pass). In that case the graph built from `current_reads` is byte-identical to
    # the graph a downstream re-assembly would build from the corrected reads, so we
    # hand it back and `_assemble_with_iterative_corrector` reuses it instead of
    # rebuilding, saving the redundant build_qualmer_graph (~22-36% of the iterative
    # arm). `final_pass_graph_reusable` is set true ONLY on a 0-improvement pass;
    # when the last pass changed reads the graph is stale and the caller rebuilds.
    # These are overwritten every pass, so after the loop they reflect the LAST
    # pass at the LARGEST k processed (the pass that produced the final FASTQ).
    final_pass_graph = nothing
    final_pass_graph_k = 0
    final_pass_graph_reusable = false

    # Main k-mer progression loop
    while k <= max_k
        if verbose
            println("\n" * "="^60)
            println("PROCESSING K-MER SIZE: $k (prime: $(Primes.isprime(k)))")
            println("="^60)
        end

        push!(k_progression, k)
        iteration_history[k] = Dict{Symbol, Any}[]

        iteration = 1
        improvements_this_k = 0
        # Declared in the outer-k scope so it is visible after the while loop
        # (it is assigned inside the loop; Julia loop scope would otherwise leave
        # it undefined at the post-loop "Final improvement rate" line — a latent
        # crash that prevented the pipeline from completing past the first k).
        current_reads = FASTX.FASTQ.Record[]

        # Soft-EM M-step memory (td-e70t v2): the PREVIOUS EM iteration's soft
        # edge-weight accumulator, registered onto THIS iteration's freshly-built
        # graph so the decode consumes the floored, probability-weighted edges.
        # Reset to `nothing` at each new k — accumulator keys are (src,dst) k-mer
        # label tuples, so weights from a different k never match and must not
        # carry over.
        prev_soft_weights = nothing

        # Iterative improvement loop for current k
        while iteration <= max_iterations_per_k
            if verbose
                println("\n--- Iteration $iteration for k=$k ---")
            end

            iteration_start = time()
            timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")

            # Read in entire FASTQ set for this iteration
            if verbose
                println("Reading FASTQ file: $current_fastq_file")
            end
            current_reads = collect(FASTX.FASTQ.Reader(open(current_fastq_file)))

            # Build qualmer graph from current read set
            if verbose
                println("Building qualmer graph with k=$k...")
            end
            graph = Mycelia.Rhizomorph.build_qualmer_graph(current_reads, k; mode = graph_mode)
            num_kmers = length(MetaGraphsNext.labels(graph))

            if verbose
                println("Graph built: $num_kmers unique k-mers")
            end

            # Check memory usage
            if !check_memory_limits(graph, memory_limit)
                if verbose
                    println("Memory limit reached. Stopping at k=$k")
                end
                break
            end

            # --- Soft-EM E-step accumulator (td-e70t, v2 COMPETING-PATHS) -----
            # A fresh accumulator tallies THIS pass's per-edge responsibilities. In
            # v2 the E-step enumerates COMPETING candidate paths per decoded read
            # (the observed read path vs a consensus alternative re-routed through
            # the best-supported sibling) and splits the responsibility across them
            # by a stable softmax over their normalized-transition log-probabilities.
            # This accumulator is then REGISTERED (support-floored) onto the NEXT
            # iteration's graph in the M-step below, closing the loop that lets
            # unsupported error edges decay while supported variation is retained.
            # Only allocated under soft-EM so :exhaustive threads `nothing` and is
            # byte-for-byte unchanged.
            current_soft_weights = soft_em ?
                                   Mycelia.Rhizomorph.SoftEdgeWeightAccumulator() : nothing

            # --- Hard-window gating (td-nn6l) --------------------------------
            # Restrict decoding to reads that touch a "hard" vertex (bubble
            # entry/exit/interior, high-out-degree/repeat-like, or weak/non-solid
            # k-mer). Easy reads pass through untouched, cutting decode volume
            # ~85-95% on clean data. Off unless `hard_window=true` (:scalable).
            # Mode-agnostic (td-nt69): `_hard_vertex_set` operates on graph
            # vertices/edges and `should_decode_read` matches reads in their observed
            # orientation, so the gate works on both :canonical and :doublestrand
            # (:scalable now runs :doublestrand). Only :singlestrand — which has no
            # RC vertices and is not a corrector target — is excluded.
            # Low-k decode gating (td-9h5r) — EXPLICIT floor. Below
            # `effective_min_decode_k` the per-read graph-Viterbi decode is skipped
            # for EVERY read this pass. Stage 0 cheap correction + skip-solid still
            # run below. Only active on the :scalable tier (hard_window=true), so
            # :exhaustive is byte-identical. The ADAPTIVE gate (density-based) is
            # decided inside improve_read_set_likelihood — it needs the post-Stage-0
            # decode fraction — and is reported back via `pass_decode_gated`.
            explicit_floor_gated = effective_min_decode_k !== nothing &&
                                   k < effective_min_decode_k
            # Building the hard-vertex set is pointless when the explicit floor
            # already gates the decode off; the adaptive gate needs it, though, so
            # build it whenever the floor is NOT gating and hard_window is on.
            hard_vertices = (hard_window && !explicit_floor_gated &&
                             graph_mode != :singlestrand) ?
                            _hard_vertex_set(graph, k) : nothing
            if hard_window && !explicit_floor_gated && graph_mode == :singlestrand
                @warn "hard_window gating is not supported for graph_mode=:singlestrand; " *
                      "disabling (decoding all non-solid reads)." graph_mode maxlog = 1
            end

            # Process each read for likelihood improvement with performance optimizations
            if verbose
                println("Processing reads for likelihood improvements...")
            end
            # `beam_width` defaults to the size-aware auto-beam (exact where
            # tractable, bounded on large reads) so the production assembler cannot
            # OOM-crash (td-63qy); :exhaustive passes typemax(Int) to force exact.
            # `corrector_diagnostics` accumulates swallowed decode failures across
            # every k and iteration and is surfaced in the result metadata.
            # --- Soft-EM M-step consumption (td-e70t v2) ---------------------
            # Register the PREVIOUS EM iteration's soft edge weights (support-
            # floored) onto THIS freshly-built graph so `compute_edge_weight` — and
            # thus the Viterbi transition scoring AND the competing-path enumeration
            # in the E-step — consumes the decayed, probability-weighted edges
            # instead of raw counts. The floor holds every >= MIN_SUPPORT edge at
            # its raw coverage (real variation preserved), so only unsupported error
            # edges decay. Paired with `clear!` in a `finally` so the process-global
            # registry never leaks past this iteration — including on :exhaustive,
            # where `soft_em` is false, nothing is registered, and the decode is
            # byte-identical. Assignments in the `try` share the enclosing scope.
            updated_reads = current_reads
            improvements_made = 0
            pass_skip_fraction = 0.0
            pass_cheap_corrections = 0
            pass_decode_gated = false
            try
                if soft_em && prev_soft_weights !== nothing
                    Mycelia.Rhizomorph.register_soft_edge_weights!(graph, prev_soft_weights)
                end
                updated_reads,
                improvements_made,
                pass_skip_fraction,
                pass_cheap_corrections,
                pass_decode_gated = improve_read_set_likelihood(
                    current_reads, graph, k,
                    verbose = verbose,
                    batch_size = batch_size,
                    enable_parallel = enable_parallel,
                    graph_mode = graph_mode,
                    skip_solid = skip_solid,
                    cheap_correct = cheap_correct,
                    beam_width = beam_width,
                    soft_weights = current_soft_weights,
                    hard_vertices = hard_vertices,
                    decode_enabled = !explicit_floor_gated,
                    decode_gate_density = effective_decode_gate_density,
                    diagnostics = corrector_diagnostics
                )
            finally
                Mycelia.Rhizomorph.clear_soft_edge_weights!()
            end
            push!(skip_fractions, pass_skip_fraction)
            push!(cheap_correction_counts, pass_cheap_corrections)
            # Record a rung whose per-read decode was gated OFF (explicit floor OR
            # adaptive density gate) — telemetry surfaced in the run metadata.
            if (explicit_floor_gated || pass_decode_gated) &&
               (isempty(decode_gated_rungs) || last(decode_gated_rungs) != k)
                push!(decode_gated_rungs, k)
            end
            # Carry this iteration's soft edge memory forward: it becomes the next
            # EM iteration's M-step input (registered onto the next graph). `nothing`
            # on :exhaustive (no soft-EM) so that tier stays byte-identical.
            prev_soft_weights = soft_em ? current_soft_weights : nothing

            # Calculate iteration metrics
            iteration_time = time() - iteration_start
            improvement_rate = improvements_made / length(current_reads)

            iteration_stats = Dict(
                :iteration => iteration,
                :k => k,
                :timestamp => timestamp,
                :improvements_made => improvements_made,
                :total_reads => length(current_reads),
                :improvement_rate => improvement_rate,
                :runtime_seconds => iteration_time,
                :memory_kmers => num_kmers
            )

            push!(iteration_history[k], iteration_stats)
            improvements_this_k += improvements_made
            total_improvements += improvements_made

            if verbose
                println("Improvements made: $improvements_made ($(round(improvement_rate * 100, digits=2))%)")
                println("Iteration runtime: $(round(iteration_time, digits=2)) seconds")
            end

            # Write out complete updated read set to timestamped FASTQ
            output_file = joinpath(output_dir, "reads_k$(k)_iter$(iteration)_$(timestamp).fastq")
            write_fastq(records = updated_reads, filename = output_file)  # Use existing function

            if verbose
                println("Wrote $(length(updated_reads)) reads to $output_file")
            end

            # Capture this pass's graph for potential final-pass reuse (td-04tb).
            # `graph` was built from `current_reads` this pass; `output_file`
            # (== updated_reads) is this pass's corrected output. When the pass made
            # 0 improvements, updated_reads == current_reads content, so `graph` is
            # byte-identical to a rebuild from the corrected reads and is REUSABLE.
            # Overwritten every pass ⇒ after the loop this holds the last pass at the
            # largest k, i.e. the pass whose FASTQ finalize selects as the result.
            final_pass_graph = graph
            final_pass_graph_k = k
            final_pass_graph_reusable = (improvements_made == 0)

            # Create checkpoint if enabled and at checkpoint interval
            if enable_checkpointing && iteration % checkpoint_interval == 0
                checkpoint_data = Dict(
                    "current_k" => k,
                    "current_iteration" => iteration,
                    "k_progression" => k_progression,
                    "iteration_history" => iteration_history,
                    "total_improvements" => total_improvements,
                    "current_fastq_file" => output_file,
                    "timestamp" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
                    "runtime_so_far" => time() - start_time
                )

                try
                    open(checkpoint_file, "w") do f
                        JSON.print(f, checkpoint_data, 2)
                    end

                    # Also save progress summary
                    progress_file = joinpath(progress_dir, "progress_k$(k)_iter$(iteration).json")
                    open(progress_file, "w") do f
                        JSON.print(f,
                            Dict(
                                "k" => k,
                                "iteration" => iteration,
                                "improvements_this_iteration" => improvements_made,
                                "total_improvements" => total_improvements,
                                "improvement_rate" => improvement_rate,
                                "runtime_seconds" => iteration_time,
                                "timestamp" => timestamp
                            ),
                            2)
                    end

                    if verbose
                        println("Checkpoint saved: k=$k, iteration=$iteration")
                    end
                catch e
                    if verbose
                        println("Warning: Could not save checkpoint: $e")
                    end
                end
            end

            # No-change convergence stop (Musket): a pass that made 0 changes is
            # converged for this k — break immediately, independent of
            # `improvement_threshold`. This is a strict superset of the
            # threshold-based early stop below (which would also fire for 0
            # improvements only when threshold > 0), and it is the guarantee that
            # keeps the per-k loop bounded even when `improvement_threshold == 0`.
            if stop_on_no_change && improvements_made == 0
                if verbose
                    println("No changes this pass (0 improvements) — converged at k=$k")
                end
                break
            end

            # Check if we should continue with this k or move to next
            if sufficient_improvements(
                improvements_made, length(current_reads), improvement_threshold,
                iteration_history = iteration_history[k])
                if verbose
                    println("Sufficient improvements detected. Continuing with k=$k")
                end
                iteration += 1
                current_fastq_file = output_file  # Use updated reads for next iteration
            else
                if verbose
                    println("Insufficient improvements or convergence detected. Moving to next k-mer size")
                end
                break
            end
        end

        if verbose
            println("\nCompleted k=$k processing:")
            println("  Total iterations: $(length(iteration_history[k]))")
            println("  Total improvements: $improvements_this_k")
            println("  Final improvement rate: $(round(improvements_this_k / length(current_reads), digits=4))")
        end

        # Move to next scheduled k-mer size (coarse ladder or prime progression)
        next_k = _next_k_in_progression(k, max_k, k_schedule)
        if next_k == k
            if verbose
                println("No larger k-mer size available. Stopping at k=$k")
            end
            break
        end
        k = next_k
    end

    # Calculate total runtime
    total_runtime = time() - start_time

    if verbose
        println("\n" * "="^60)
        println("ITERATIVE ASSEMBLY COMPLETE")
        println("="^60)
        println("Total runtime: $(round(total_runtime, digits=2)) seconds")
        println("K-mer sizes processed: $(sort(k_progression))")
        println("Total improvements: $total_improvements")
    end

    # Finalize iterative assembly
    return finalize_iterative_assembly(
        output_dir, k_progression, iteration_history, total_runtime,
        verbose = verbose, diagnostics = corrector_diagnostics,
        skip_fractions = skip_fractions,
        cheap_correction_counts = cheap_correction_counts,
        hard_window = hard_window, soft_em = soft_em, cheap_correct = cheap_correct,
        min_decode_k = effective_min_decode_k,
        decode_gate_density = effective_decode_gate_density,
        decode_gated_rungs = decode_gated_rungs,
        # Final-pass graph reuse (td-04tb): hand the corrector's last-pass qualmer
        # graph back so a converged run's re-assembly can skip a redundant rebuild.
        final_pass_graph = final_pass_graph,
        final_pass_graph_k = final_pass_graph_k,
        final_pass_graph_mode = graph_mode,
        final_pass_graph_reusable = final_pass_graph_reusable)
end

# =============================================================================
# Read Likelihood Improvement Functions
# =============================================================================

# ------------------------------------------------------------------------------
# Stage 0 skip support (td-1do7): classify k-mers once per graph so that reads
# with no weak k-mers can skip the expensive per-read Viterbi decode. This is the
# volume reduction that motivated the staged-correction architecture — the
# corrector should only decode reads that actually need correction.
# ------------------------------------------------------------------------------

"""
    _solid_kmer_set(graph; classifier) -> Set

Classify every k-mer once and return the set of SOLID (real-genomic) canonical
k-mer labels. Defaults to the auto-threshold `MixtureModelClassifier` (a top
default in the Stage 0 comparison matrix, coverage-based so it needs no per-arm
tuning). Returns the EMPTY set (⇒ `_read_is_all_solid` skips nothing ⇒ correct
every read) when classification cannot run — no labels or no dataset evidence —
so the skip is never applied on evidence it could not compute (review C1). The
conservative default is "decode everything", NOT "skip everything": returning all
labels would mark error k-mers (graph vertices at coverage 1) as solid and
silently skip reads that needed correction.
"""
function _solid_kmer_set(graph;
        classifier = Mycelia.Rhizomorph.MixtureModelClassifier())
    labels = collect(MetaGraphsNext.labels(graph))
    isempty(labels) && return Set{eltype(labels)}()
    dataset_id = _rhizomorph_first_dataset_id(graph[first(labels)])
    if dataset_id === nothing
        @warn "skip_solid: graph has no dataset evidence; classification skipped, correcting all reads (no skip)."
        return Set{eltype(labels)}()   # empty ⇒ skip nothing ⇒ correct everything
    end
    classification = Mycelia.Rhizomorph.classify_kmers(classifier, graph;
        dataset_id = dataset_id)
    solid = Set{eltype(labels)}()
    for (label, is_solid) in classification.verdicts
        is_solid && push!(solid, label)
    end
    return solid
end

# Mode-aware lookup key for classification/gating (td-nt69). `solid_kmers` /
# `hard_vertices` are graph VERTEX LABELS, so the key that a read's observed k-mer
# maps to depends on how the graph was built:
#   :canonical    → labels are canonical k-mers, so canonicalize the read k-mer.
#   :doublestrand → both a k-mer and its RC are SEPARATE vertices, so the read
#                   traverses its observed-orientation k-mer; use it as-is.
#   :singlestrand → no RC vertices; use the observed k-mer as-is.
# Classification is coverage-based and therefore mode-agnostic: under
# :doublestrand each strand-vertex still separates real (≈half-coverage) from
# error (coverage-1) k-mers, so the solid/hard sets remain valid.
@inline function _lookup_key(kmer, graph_mode::Symbol)
    return graph_mode == :canonical ? BioSequences.canonical(kmer) : kmer
end

"""
    _read_is_all_solid(read, k, solid_kmers; graph_mode) -> Bool

True iff every k-mer of `read` (resolved to its graph label via `_lookup_key`) is
in `solid_kmers` (⇒ no weak region ⇒ the read needs no correction and can skip the
decode). A read with no k-mers (shorter than k, or fully ambiguous) returns
`false` so it is never skipped on the basis of absent evidence. Works on both
`:canonical` and `:doublestrand` graphs (td-nt69).
"""
function _read_is_all_solid(read::FASTX.FASTQ.Record, k::Int, solid_kmers::AbstractSet;
        graph_mode::Symbol = :canonical)
    isempty(solid_kmers) && return false
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    saw_kmer = false
    for (kmer, _) in Kmers.UnambiguousDNAMers{k}(sequence)
        saw_kmer = true
        (_lookup_key(kmer, graph_mode) in solid_kmers) || return false
    end
    return saw_kmer
end

# ------------------------------------------------------------------------------
# Stage 0 CHEAP k-mer-spectrum correction (td-bjnt): fix simple single-base
# substitution errors with a LINEAR scan, BFC/Lighter/Bloocoo-style, BEFORE the
# expensive per-read graph Viterbi. A WEAK (low-coverage, non-solid) k-mer flanked
# by SOLID k-mers is the signature of a single-base error; if exactly one base
# substitution turns the whole weak run solid (a UNIQUE solid neighbor), apply it.
# Ambiguous cases (no fix, or >1 candidate — e.g. a balanced heterozygous site
# where both alleles are real) are LEFT untouched for the graph decode. Reserving
# graph Viterbi for genuine ambiguity is the critical-path lever: on err=0.01 data
# ~78% of reads carry an error k-mer, so cheaply clearing the simple ones collapses
# the graph-decode fraction toward the true bubble/repeat ~5–15%.
# Only invoked on the :scalable tier (`cheap_correct=true`, canonical graph); the
# :exhaustive path never calls it and is byte-identical.
# ------------------------------------------------------------------------------

# Graph-label lookup key for the k-mer at 1-based read position `i` (bases
# [i, i+k-1]) of an ACGT-only char vector, comparable against the `_solid_kmer_set`
# labels. Mode-aware (td-nt69): canonicalized under :canonical, observed
# orientation under :doublestrand / :singlestrand. The read is pre-screened for
# ambiguity by the caller so this construction cannot throw.
@inline function _lookup_kmer_at(chars::Vector{Char}, k::Int, i::Int, graph_mode::Symbol)
    return _lookup_key(Kmers.DNAKmer{k}(String(@view chars[i:(i + k - 1)])), graph_mode)
end

"""
    _stage0_correct_read(read, k, solid_kmers) -> (record, n_corrections)

Cheaply correct simple single-substitution errors in ONE read via the k-mer
spectrum. Walk the read's k-mers; for each maximal run of consecutive WEAK
(non-solid) k-mers that is flanked by SOLID k-mers on BOTH sides (an interior
error) and is short enough to be a single substitution (run length ≤ k), search
the base position(s) shared by every k-mer in the run for a substitution that
makes EVERY k-mer overlapping that position solid. Apply the fix only when EXACTLY
ONE (position, base) candidate qualifies — a unique solid neighbor. Zero
candidates (a real but low-coverage allele with no solid neighbor) or multiple
candidates (a balanced variant with two solid alleles, or genuine ambiguity) are
left untouched for the graph decode, so Stage 0 never collapses real variation.

Reads shorter than `k`, reads with ambiguous bases (e.g. `N`), and reads with no
weak run are returned unchanged with a 0 count. Returns the (possibly rewritten)
record and the number of bases corrected. `graph_mode` selects the solid-set
lookup key (canonicalized under :canonical, observed orientation otherwise —
td-nt69).
"""
function _stage0_correct_read(read::FASTX.FASTQ.Record, k::Int, solid_kmers::AbstractSet;
        graph_mode::Symbol = :canonical)
    seq_str = FASTX.sequence(String, read)
    n = length(seq_str)
    n < k && return read, 0
    chars = collect(seq_str)
    # Ambiguity screen: the 2-bit k-mer construction cannot encode non-ACGT bases.
    # Such reads are left for the decode (which tallies un-k-merizable reads).
    for c in chars
        (c == 'A' || c == 'C' || c == 'G' || c == 'T') || return read, 0
    end

    nkmers = n - k + 1
    solid = Vector{Bool}(undef, nkmers)
    @inbounds for i in 1:nkmers
        solid[i] = _lookup_kmer_at(chars, k, i, graph_mode) in solid_kmers
    end

    n_corr = 0
    i = 1
    while i <= nkmers
        if solid[i]
            i += 1
            continue
        end
        # Maximal weak run [a, b].
        a = i
        b = i
        while b + 1 <= nkmers && !solid[b + 1]
            b += 1
        end
        run_len = b - a + 1
        # Interior error: flanked by solid on both sides, run short enough to be a
        # single substitution (the base shared by all k-mers a..b is non-empty).
        flanked = (a - 1 >= 1 && solid[a - 1]) && (b + 1 <= nkmers && solid[b + 1])
        if flanked && run_len <= k
            # Base positions shared by every k-mer in [a, b]: [b, a+k-1] (1-based
            # read coordinates). A single substitution at one of these explains the
            # whole run.
            lo = b
            hi = a + k - 1
            fix_p = 0
            fix_base = ' '
            n_candidates = 0
            for p in lo:hi
                cur = chars[p]
                # k-mers whose window covers base p — the ONLY ones a substitution
                # at p can change. Requiring all of these solid after the edit both
                # fixes the run and guarantees no flanking k-mer is broken.
                lo_j = max(1, p - k + 1)
                hi_j = min(nkmers, p)
                for base in ('A', 'C', 'G', 'T')
                    base == cur && continue
                    chars[p] = base
                    ok = true
                    for j in lo_j:hi_j
                        if !(_lookup_kmer_at(chars, k, j, graph_mode) in solid_kmers)
                            ok = false
                            break
                        end
                    end
                    chars[p] = cur   # restore before testing the next candidate
                    if ok
                        n_candidates += 1
                        n_candidates > 1 && break   # ambiguous: stop early
                        fix_p = p
                        fix_base = base
                    end
                end
                n_candidates > 1 && break
            end
            if n_candidates == 1
                chars[fix_p] = fix_base
                # The corrected k-mers overlapping fix_p are now solid; mark them so
                # the scan does not re-enter this (now resolved) run.
                lo_j = max(1, fix_p - k + 1)
                hi_j = min(nkmers, fix_p)
                for j in lo_j:hi_j
                    solid[j] = true
                end
                n_corr += 1
            end
        end
        i = b + 1
    end

    n_corr == 0 && return read, 0
    corrected_seq = String(chars)
    # Preserve identifier, description, and per-base quality (a substitution does
    # not change the length, so the original quality string still aligns).
    new_record = FASTX.FASTQ.Record(
        FASTX.identifier(read), corrected_seq, FASTX.quality(read))
    return new_record, n_corr
end

"""
    _stage0_cheap_correct(reads, k, solid_kmers) -> (corrected_reads, n_corrections)

Apply `_stage0_correct_read` across a read set, returning a NEW vector of records
(unchanged records are shared, corrected ones replaced) and the total number of
bases corrected. An empty `solid_kmers` (classification could not run) is a no-op:
without a solid reference every k-mer looks weak and there is no trustworthy
neighbor to correct toward, so nothing is changed.
"""
function _stage0_cheap_correct(reads::Vector{<:FASTX.FASTQ.Record}, k::Int,
        solid_kmers::AbstractSet; graph_mode::Symbol = :canonical)
    total = 0
    isempty(solid_kmers) && return collect(reads), 0
    corrected = Vector{FASTX.FASTQ.Record}(undef, length(reads))
    for (idx, read) in enumerate(reads)
        rec, n = _stage0_correct_read(read, k, solid_kmers; graph_mode = graph_mode)
        corrected[idx] = rec
        total += n
    end
    return corrected, total
end

# ------------------------------------------------------------------------------
# Hard-window gating (td-nn6l): decode only reads that touch a "hard" region
# (bubble / repeat) and pass every other read through untouched. This is the
# primary decode-volume reduction on top of skip-solid and Stage 0 cheap
# correction — after Stage 0 clears simple errors, the vast majority of reads
# touch no hard (bubble/repeat) vertex and are skipped.
# ------------------------------------------------------------------------------

"""
    should_decode_read(read, k, hard_vertices; graph_mode) -> Bool

True iff any k-mer of `read` (resolved to its graph label via `_lookup_key`) is a
member of `hard_vertices` (the bubble / repeat-like / weak-k-mer set). Reads that
touch no hard vertex have no region worth decoding and are passed through
untouched. Mirrors `_read_is_all_solid`'s k-mer iteration so the two gates compose
on the same graph. Works on both `:canonical` and `:doublestrand` graphs
(td-nt69). An empty `hard_vertices` means "no gate" ⇒ decode.
"""
function should_decode_read(read::FASTX.FASTQ.Record, k::Int, hard_vertices::AbstractSet;
        graph_mode::Symbol = :canonical)
    isempty(hard_vertices) && return true
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    for (kmer, _) in Kmers.UnambiguousDNAMers{k}(sequence)
        (_lookup_key(kmer, graph_mode) in hard_vertices) && return true
    end
    return false
end

"""
    _hard_vertex_set(graph, k) -> Set

Build the set of "hard" vertices for hard-window gating (td-nn6l/td-bjnt): the
union of
(1) bubble entry/exit/interior vertices (`detect_bubbles_next`), and
(2) high-out-degree (out-degree > 1, repeat-like) vertices.

NARROWED for Stage 0 (td-bjnt): the former clause (3) — "every weak / non-solid
k-mer is hard" — is REMOVED. Under the old clause any read carrying a single
error k-mer (≈78% of reads on err=0.01 data) was classified hard and sent to the
expensive per-read graph Viterbi, which is what times the :scalable corrector out
at 5 kb. Simple single-substitution errors are now fixed CHEAPLY by the Stage 0
k-mer-spectrum pass (`_stage0_cheap_correct`, a linear scan) BEFORE this gate, so
the only reads that must still reach graph Viterbi are those touching GENUINE
ambiguity — a bubble/superbubble (competing balanced alleles) or a repeat-like
high-out-degree vertex. Those are the true ~5–15%. A weak k-mer flanked by solid
neighbors is no longer "hard"; it is an error Stage 0 either already corrected or
left as genuinely ambiguous (in which case it typically also sits in a bubble and
is caught by clause 1). Operates directly on the graph vertices/edges, so it is
mode-agnostic (works on `:canonical` and `:doublestrand` graphs — on a
doublestrand graph a bubble/repeat is flagged on both strands, and reads are
matched in their observed orientation by `should_decode_read`; td-nt69). Only
reached on the :scalable tier (`hard_window=true`), so :exhaustive is unaffected.
"""
function _hard_vertex_set(graph, k::Int)
    labels = collect(MetaGraphsNext.labels(graph))
    T = eltype(labels)
    hard = Set{T}()
    isempty(labels) && return hard

    # (1) Bubble vertices (entry, exit, both alternative interiors) — REAL
    # competing-path ambiguity that a linear cheap-correction cannot resolve.
    for bubble in Mycelia.Rhizomorph.detect_bubbles_next(graph)
        push!(hard, bubble.entry_vertex)
        push!(hard, bubble.exit_vertex)
        for v in bubble.path1
            push!(hard, v)
        end
        for v in bubble.path2
            push!(hard, v)
        end
    end

    # (2) High out-degree (repeat-like) vertices — single O(E) pass over edges.
    # A branch point with >1 successor is a repeat/ambiguity the per-read decode
    # must disambiguate; a cheap linear scan cannot.
    outdeg = Dict{T, Int}()
    for (src, _dst) in MetaGraphsNext.edge_labels(graph)
        outdeg[src] = get(outdeg, src, 0) + 1
    end
    for (v, d) in outdeg
        d > 1 && push!(hard, v)
    end

    return hard
end

# ------------------------------------------------------------------------------
# Stage 3c SCAFFOLD (td-nn6l): per-hard-region WINDOWED decode.
#
# STATUS: NOT wired into the gate. Stage 3 currently decodes each hard read
# WHOLE (the skip-gate part of Stage 3 is fully delivered — easy reads are
# skipped, hard reads are decoded). The windowed variant below is the primitive
# for the remaining 3c work: instead of decoding a hard read end-to-end, extract
# the k-mer sub-window (<=500 bp) around each hard region and decode ONLY that
# window with start_vertex/target_vertex boundary constraints (which
# `Mycelia.correct_observations` supports via the observation's first/last vertex
# and `ViterbiCorrectionConfig.target_vertex`), then splice the corrected window
# back into the read. Windowing bounds decode cost to the hard neighborhood
# rather than the whole read.
#
# `_hard_window_ranges` (the read-coordinate ranges to decode) is implemented and
# unit-tested; the correct_observations-with-boundaries call + splice is the
# deferred wiring. See the PR description for the 3c caveat.
# ------------------------------------------------------------------------------

"""
    _hard_window_ranges(read, k, hard_vertices; pad=1, max_window=500) -> Vector{UnitRange{Int}}

Read-coordinate ranges (1-based, in read bases) covering each hard region of a
read: for every k-mer position whose canonical k-mer is in `hard_vertices`, take
the base span `[pos-pad, pos+k-1+pad]`, clamp to the read, merge overlapping
spans, and cap each merged span at `max_window` bases. These are the windows a
Stage 3c windowed decode would correct in isolation (whole-read decode is the
union of all bases; windowing decodes only these). Empty when the read touches
no hard vertex.
"""
function _hard_window_ranges(read::FASTX.FASTQ.Record, k::Int, hard_vertices::AbstractSet;
        pad::Int = 1, max_window::Int = 500)
    ranges = UnitRange{Int}[]
    isempty(hard_vertices) && return ranges
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    n = length(sequence)
    n < k && return ranges
    for (km, kpos) in Kmers.UnambiguousDNAMers{k}(sequence)
        if BioSequences.canonical(km) in hard_vertices
            lo = max(1, kpos - pad)
            hi = min(n, kpos + k - 1 + pad)
            push!(ranges, lo:hi)
        end
    end
    isempty(ranges) && return ranges
    # Merge overlapping/adjacent ranges, capping each merged span at max_window.
    sort!(ranges; by = first)
    merged = UnitRange{Int}[]
    cur = ranges[1]
    for r in ranges[2:end]
        if first(r) <= last(cur) + 1
            cur = first(cur):max(last(cur), last(r))
        else
            push!(merged, cur)
            cur = r
        end
    end
    push!(merged, cur)
    return [first(r):min(last(r), first(r) + max_window - 1) for r in merged]
end

# ------------------------------------------------------------------------------
# Corrector diagnostics + size-aware auto-beam (review: robustness blockers)
# ------------------------------------------------------------------------------

"""
    CorrectorDiagnostics()

Thread-safe tally of decode outcomes the per-read corrector would otherwise
swallow SILENTLY. When `try_viterbi_path_improvement` fails to decode a read it
returns `nothing` and the read passes through uncorrected; without this tally a
SYSTEMATICALLY broken corrector — empty graph, alphabet-inference miss, config
error, or a read that cannot be k-merized — reports `0/N improvements` as if it
had actually run and found nothing to fix. Counters are `Threads.Atomic` so the
parallel (`@threads`) read loop increments them race-free.

- `structural_errors`  : `ArgumentError` from the decoder (empty graph,
  alphabet-inference miss, config error). Means "the decoder could not run",
  NOT "the decoder ran and found no gain".
- `unkmerizable_reads` : `BioSequences.EncodeError` — a read carrying an
  ambiguous base (e.g. `N`) that the 2-bit k-mer iterator cannot encode. The
  read is SKIPPED (passes through uncorrected), never crashed on.
"""
mutable struct CorrectorDiagnostics
    structural_errors::Threads.Atomic{Int}
    unkmerizable_reads::Threads.Atomic{Int}
end
function CorrectorDiagnostics()
    CorrectorDiagnostics(Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
end

# Previously-proven-tractable finite beam (td-63qy: beam 256 completed on the
# 48 kb phage that OOM-crashed — ~21B allocations — at the unbounded typemax
# default). A read whose observation count is at/below the threshold keeps the
# EXACT maximum-likelihood guarantee (its frontier is small); above the threshold
# the frontier is bounded so the SHIPPING assembler cannot OOM-crash.
const _AUTO_BEAM_BOUNDED_WIDTH = 256
# Heuristic cutoff on a read's observation (k-mer) count. Typical short reads
# (Illumina, up to a few hundred bp) stay exact; long reads (PacBio/ONT, or the
# 48 kb-scale inputs behind td-63qy) get bounded. This is NOT a correctness
# constant: it only trades exactness for tractability, and only ABOVE the line.
const _AUTO_BEAM_EXACT_THRESHOLD = 1024

# Candidate-GENERATION bounds applied alongside the finite width beam (td-plqi).
# The width beam caps the RETAINED frontier at 256, but on a dense intermediate-k
# graph nearly all 256 retained states are >20 log-prob below the best — improbable
# junk that can never win the ML path yet still generates + emission-scores
# successors at the next depth, so per-depth generation climbs toward the width cap
# as the graph densifies (the empirically-measured residual super-linear term at
# k=9; #386). These bound generation directly:
#
#   * _AUTO_SUCCESSOR_BOUND — top-B outgoing transitions per expanded state, ranked
#     by edge weight before emission. A k-mer de Bruijn vertex has ≤ 4 structural
#     successors, so 16 is a strict no-op on DNA (a robustness guard for pathological
#     high-branching / non-DNA inputs), kept for correctness parity with the design.
#
#   * _AUTO_BEAM_SCORE_MARGIN — Δ log-prob "histogram" beam with an EMISSION
#     exemption: prune a frontier state only when it is >Δ below the depth's best
#     on BOTH the full score AND the cumulative emission (read-consistency). This
#     linearizes the dense-rung decode — the read-INCONSISTENT junk (low on both
#     axes) is discarded so the generating frontier stays O(1) in size — WITHOUT
#     dropping a real-but-rare allele: a supported minor allele in a skewed pool
#     has good emission but a coverage-driven transition penalty, and the emission
#     clause exempts it from pruning (the margin removes WRONG paths, not merely
#     RARE ones; PR #388 variation-safety review). Δ = 30 nats. Both engage ONLY
#     where the width beam is already finite (approximate); exact-ML reads
#     (beam_width == typemax) keep margin = Inf and stay byte-identical.
const _AUTO_SUCCESSOR_BOUND = 16
const _AUTO_BEAM_SCORE_MARGIN = 30.0

"""
    _auto_beam_width(n_observations) -> Int

Size-aware default beam width reconciling the ADR's "exact by default" with the
td-63qy OOM crash. Returns `typemax(Int)` (exact — the unbounded frontier
preserves the ML guarantee) when `n_observations <= _AUTO_BEAM_EXACT_THRESHOLD`,
and the bounded `_AUTO_BEAM_BOUNDED_WIDTH` above it, emitting a one-line `@info`
when it bounds so the exactness→tractability switch is visible in the log.
Callers can still FORCE exact on a large read by passing an explicit
`beam_width = typemax(Int)`.
"""
function _auto_beam_width(n_observations::Integer)::Int
    if n_observations <= _AUTO_BEAM_EXACT_THRESHOLD
        return typemax(Int)
    end
    @info "iterative corrector: read has $(n_observations) observations " *
          "(> $(_AUTO_BEAM_EXACT_THRESHOLD)); bounding Viterbi beam_width=" *
          "$(_AUTO_BEAM_BOUNDED_WIDTH) to stay tractable (exact ML would risk OOM, " *
          "td-63qy). Pass beam_width=typemax(Int) to force exact." maxlog = 3
    return _AUTO_BEAM_BOUNDED_WIDTH
end

"""
    _auto_beam_width(n_observations, n_vertices) -> Int

Graph-density-aware size-aware default beam width (td-35ux). The read-length-only
[`_auto_beam_width`](@ref) keeps the exact/unbounded frontier for any read with
`n_observations <= _AUTO_BEAM_EXACT_THRESHOLD` (typical short reads). That is
correct on a sparse graph, but at a dense intermediate k-rung the EXACT retained
(vertex, strand) frontier grows `O(n_vertices)` per decode depth — i.e.
`O(genome)` per pass, `O(genome^2)` over the read set — even for a single short
150 bp read (~15,723 frontier states at 5 kb / k=9, ~2.15 s per read). Long reads
already dodge this because the read-length rule bounds them; short reads did not,
which is why only short reads blew up (#376).

This two-arg form bounds the beam to `_AUTO_BEAM_BOUNDED_WIDTH` whenever EITHER
the read is large (read-length rule) OR the graph is dense
(`n_vertices > _AUTO_BEAM_BOUNDED_WIDTH`), so a short read on a big graph is
capped regardless of its length. On a small/sparse graph
(`n_vertices <= _AUTO_BEAM_BOUNDED_WIDTH`) the exact frontier is already tiny, so
a small read keeps the exact `typemax(Int)` (ML guarantee preserved). Callers can
still FORCE exact by passing an explicit `beam_width = typemax(Int)`.
"""
function _auto_beam_width(n_observations::Integer, n_vertices::Integer)::Int
    # Read-length rule first: a large read is bounded no matter the graph.
    by_length = _auto_beam_width(n_observations)
    if by_length != typemax(Int)
        return by_length
    end
    # The read is small (exact by read length). If the graph is dense the exact
    # frontier can still grow O(n_vertices) per depth → O(genome^2) per pass, so
    # bound it. On a sparse graph the exact frontier is tiny; keep it exact.
    if n_vertices > _AUTO_BEAM_BOUNDED_WIDTH
        @info "iterative corrector: short read ($(n_observations) observations) on a " *
              "dense graph ($(n_vertices) vertices > $(_AUTO_BEAM_BOUNDED_WIDTH)); " *
              "bounding Viterbi beam_width=$(_AUTO_BEAM_BOUNDED_WIDTH) to keep the " *
              "retained frontier O(1) instead of O(n_vertices) per depth (td-35ux). " *
              "Pass beam_width=typemax(Int) to force exact." maxlog = 3
        return _AUTO_BEAM_BOUNDED_WIDTH
    end
    return typemax(Int)
end

"""
Improve likelihood of entire read set using current graph and k-mer size.
Returns updated reads and count of improvements made.
Uses memory-efficient batch processing for large datasets.

`beam_width` defaults to `nothing` (size-aware auto-beam per read via
`_auto_beam_width`): exact ML where the read is small enough to be tractable,
a bounded frontier on large reads so the SHIPPING assembler cannot OOM-crash
(td-63qy). Pass an explicit `Int` (including `typemax(Int)`) to override — a
caller/test passing `beam_width = typemax(Int)` still gets exact decoding on
every read. Swallowed structural/un-k-merizable decode failures are tallied into
`diagnostics` (auto-created if not supplied) and a high swallowed fraction is
`@warn`ed so a run that silently corrected nothing is visible.

`hard_vertices` (td-nn6l): when non-`nothing`, restricts decoding to reads whose
k-mers overlap the "hard" vertex set (bubbles / repeats / weak k-mers); every
other read passes through untouched. Composes with `skip_solid` — a read is
decoded only when it is NOT all-solid AND (no hard set, or it touches one).

`decode_enabled` (td-9h5r): when `false`, the per-read graph-Viterbi decode is
skipped for EVERY read this pass (skip_fraction == 1.0) — the EXPLICIT low-k decode
floor. Stage 0 cheap correction still runs first (on `work_reads`), so simple
errors are still fixed; only the expensive whole-read Viterbi is deferred to higher
k-rungs. Defaults to `true`; the `:scalable` k-ladder sets it `false` below
`min_decode_k`.

`decode_gate_density` (td-9h5r): the ADAPTIVE low-k decode gate. When non-`nothing`
and a hard-window gate is active, the pass measures the post-Stage-0 fraction of
reads the gate would still decode; if that is `>= decode_gate_density` the gate is
NON-DISCRIMINATING (the dense-low-k "gate skips nothing" pathology) and the whole
decode is skipped for the pass. A genuinely selective low-k decode (gate skips a
meaningful fraction) is left ON. `nothing` (default / `:exhaustive`) disables it.

Returns `(updated_reads, improvements_made, skip_fraction, cheap_corrections,
decode_gated)` where `decode_gated` is `true` iff the per-read decode was gated OFF
for this pass (explicit floor OR adaptive density). Callers destructuring only the
first two/three/four values are unaffected.

`soft_weights` (td-e70t): when non-`nothing`, each decoded read's ML path
accumulates edge responsibilities into it (the soft-EM E-step). Accumulation is
NOT thread-safe, so supplying it forces sequential processing for this pass.

`cheap_correct` (td-bjnt): when `true` (and `graph_mode==:canonical`), run the
Stage 0 k-mer-spectrum correction pass (`_stage0_cheap_correct`) over the read set
BEFORE any gating/decode, cheaply fixing simple single-substitution errors with a
linear scan so the expensive graph Viterbi is reserved for genuine ambiguity
(bubbles/repeats). The decode then operates on the cheaply-corrected reads. Only
enabled on the :scalable tier.

Returns `(updated_reads, improvements_made, skip_fraction, cheap_corrections,
decode_gated)` where `skip_fraction` is the fraction of reads passed through WITHOUT
a decode (solid + hard-window skips, plus the whole set when the low-k decode gate
fires) and `cheap_corrections` is the number of bases fixed by the Stage 0 pass. The
graph-Viterbi decode fraction is `1 - skip_fraction`. `improvements_made` counts
BOTH cheap corrections and decode-accepted corrections so per-k convergence sees
Stage 0 progress. `decode_gated` is `true` iff the per-read decode was gated OFF for
the whole pass (low-k gate, td-9h5r). Callers destructuring only the first
two/three/four values are unaffected.
"""
function improve_read_set_likelihood(reads::Vector{<:FASTX.FASTQ.Record}, graph, k::Int;
        verbose::Bool = false,
        batch_size::Int = 10000,
        enable_parallel::Bool = false,
        graph_mode::Symbol = :canonical,
        skip_solid::Bool = false,
        cheap_correct::Bool = false,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        hard_vertices::Union{Nothing, AbstractSet} = nothing,
        decode_enabled::Bool = true,
        decode_gate_density::Union{Float64, Nothing} = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing)::Tuple{
        Vector{FASTX.FASTQ.Record}, Int, Float64, Int, Bool}
    diag = diagnostics === nothing ? CorrectorDiagnostics() : diagnostics
    # Snapshot so the per-call @warn reflects THIS pass's swallowed fraction even
    # when `diag` is a shared accumulator threaded across many passes.
    struct_before = diag.structural_errors[]
    unk_before = diag.unkmerizable_reads[]
    total_reads = length(reads)
    updated_reads = Vector{FASTX.FASTQ.Record}(undef, total_reads)
    improvements_made = 0

    # Stage 0 skip (td-1do7): classify k-mers once; a read whose every k-mer is
    # solid has no weak region to correct and skips the per-read decode. Opt-in
    # (default off) so existing callers are unchanged. `solid_kmers === nothing`
    # ⇒ correct every read as before.
    # td-nt69: the skip/cheap-correct helpers resolve read k-mers to their graph
    # label via `_lookup_key` (canonicalized under :canonical, observed orientation
    # under :doublestrand), so they are valid on BOTH modes — classification is
    # coverage-based and mode-agnostic. :scalable now runs :doublestrand. Only
    # :singlestrand (no RC vertices, not a corrector target) is excluded.
    if skip_solid && graph_mode == :singlestrand
        @warn "skip_solid is not supported for graph_mode=:singlestrand; disabling skip (correcting all reads)." graph_mode
    end
    if cheap_correct && graph_mode == :singlestrand
        @warn "cheap_correct is not supported for graph_mode=:singlestrand; disabling Stage 0 cheap correction." graph_mode
    end
    # Classify k-mers once if EITHER the skip-solid gate OR the Stage 0 cheap
    # corrector needs the solid set. Compute once, share both consumers.
    need_solid = (skip_solid || cheap_correct) && graph_mode != :singlestrand
    solid_kmers = need_solid ? _solid_kmer_set(graph) : nothing

    # Stage 0 CHEAP correction (td-bjnt): fix simple single-substitution errors
    # with a linear k-mer-spectrum scan BEFORE any gating/decode, so graph Viterbi
    # is reserved for genuine ambiguity. The decode below then runs on the
    # cheaply-corrected reads. Gated on :scalable (cheap_correct=true, canonical).
    cheap_corrections = 0
    work_reads = reads
    if cheap_correct && graph_mode != :singlestrand && solid_kmers !== nothing
        work_reads, cheap_corrections = _stage0_cheap_correct(reads, k, solid_kmers;
            graph_mode = graph_mode)
        improvements_made += cheap_corrections
        if verbose && cheap_corrections > 0
            println("  Stage 0 cheap correction: fixed $cheap_corrections base(s) " *
                    "before decode")
        end
    end
    skipped_reads = 0

    # A read is SKIPPED (passed through uncorrected, no decode) when it is either
    # all-solid (skip-solid gate) OR — under hard-window gating — it touches no
    # hard vertex. Both are volume-reduction skips and both count toward the
    # reported skip fraction. `solid_kmers` may be populated for cheap correction
    # even when skip_solid is off, so the all-solid skip is gated on `skip_solid`
    # explicitly. Returns true ⇒ skip. Stage 0 has ALREADY run above (on
    # `work_reads`), so this predicate sees the cheaply-corrected reads.
    _gate_skip = read -> begin
        (skip_solid && solid_kmers !== nothing &&
             _read_is_all_solid(read, k, solid_kmers; graph_mode = graph_mode)) &&
            return true
        (hard_vertices !== nothing &&
             !should_decode_read(read, k, hard_vertices; graph_mode = graph_mode)) &&
            return true
        return false
    end

    # Precompute per-read skip decisions ONCE (cheap k-mer-membership checks) so the
    # adaptive low-k gate can measure the natural decode fraction and the decode loop
    # can reuse the same flags (no double evaluation).
    base_skip_flags = Vector{Bool}(undef, total_reads)
    @inbounds for i in 1:total_reads
        base_skip_flags[i] = _gate_skip(work_reads[i])
    end
    natural_decode_fraction = total_reads > 0 ?
                              count(!, base_skip_flags) / total_reads : 0.0

    # -- Low-k decode gating (td-9h5r) -----------------------------------------
    # `decode_enabled == false` (EXPLICIT floor `min_decode_k`): skip the per-read
    # graph-Viterbi decode for EVERY read this pass. Stage 0 cheap correction has
    # already run (on `work_reads`), so error correction still happens; only the
    # expensive whole-read Viterbi is deferred to higher, selective k-rungs.
    #
    # ADAPTIVE density gate: even with `decode_enabled`, if the hard-window gate
    # would still send >= `decode_gate_density` of reads to the decode (i.e. it is
    # NON-DISCRIMINATING — the dense-low-k pathology where "the gate skips nothing"),
    # skip the whole decode too: the gate is not reducing volume, so the decode is
    # ~pure waste Stage 0 + the selective higher rungs recover. Requires an active
    # hard-window gate (`hard_vertices !== nothing`) — with no gate there is no
    # discrimination signal to act on, so the decode runs as before.
    adaptive_gated = decode_enabled && decode_gate_density !== nothing &&
                     hard_vertices !== nothing &&
                     natural_decode_fraction >= decode_gate_density
    pass_decode_off = !decode_enabled || adaptive_gated
    if verbose && adaptive_gated
        println("  Low-k decode gate (td-9h5r): hard-window gate non-discriminating " *
                "at k=$k (would decode $(round(natural_decode_fraction * 100, digits=1))% " *
                ">= $(round(decode_gate_density * 100, digits=1))%); skipping decode " *
                "this pass (Stage 0 + higher-k rungs recover).")
    end

    # Final per-read skip decision: all reads skip when the pass decode is gated off;
    # otherwise use the precomputed gate flags.
    _skip_this_read_at = i -> pass_decode_off || base_skip_flags[i]

    # Soft-EM accumulation into `soft_weights` is not thread-safe (shared Dict),
    # so a soft-EM pass runs sequentially. Otherwise honor the caller's request.
    use_parallel = enable_parallel && Threads.nthreads() > 1 && soft_weights === nothing
    if enable_parallel && soft_weights !== nothing && Threads.nthreads() > 1
        @warn "soft-EM edge accumulation is sequential (race-free); ignoring enable_parallel for this pass." maxlog = 1
    end

    if verbose
        println("  Processing $total_reads reads in batches of $batch_size")
    end

    # Hoist the weighted-decode graph OUT of the per-read loop (td-y4oj). The
    # per-read `correct_observations` previously rebuilt an O(V+E) weighted COPY
    # of the whole graph ONCE PER READ; since n_reads ∝ genome and V ∝ genome,
    # that made the decode O(genome^2) (the dominant #372 super-linear term,
    # alpha≈2.14). Build it once here and thread the shared, read-only instance
    # into every per-read decode so the rebuild cost is amortized across the pass
    # (linear in n_reads). Returns `nothing` for graphs without a weighted-decode
    # path (legacy / non-MetaGraph), which the callee handles by falling back to
    # its unchanged per-read behavior. Skip the build entirely when the pass
    # decode is gated off (no read decodes this pass).
    pass_weighted_graph = pass_decode_off ? nothing : build_correction_weighted_graph(graph)

    # Process in batches for memory efficiency. `work_reads` is the Stage 0
    # cheaply-corrected read set (== `reads` when cheap_correct is off), so the
    # decode operates on already-simplified reads.
    for batch_start in 1:batch_size:total_reads
        batch_end = min(batch_start + batch_size - 1, total_reads)
        batch_reads = work_reads[batch_start:batch_end]

        batch_improvements = 0

        if use_parallel
            # Parallel processing for large batches
            batch_results = Vector{Tuple{FASTX.FASTQ.Record, Bool}}(undef, length(batch_reads))

            # Vector{Bool} (one byte/elem), NOT falses()/BitVector (64 bits share a
            # word) — the @threads writes below would otherwise race on the shared
            # word and undercount the skip fraction (review I2).
            skip_flags = fill(false, length(batch_reads))
            Threads.@threads for i in eachindex(batch_reads)
                read = batch_reads[i]
                if _skip_this_read_at(batch_start + i - 1)
                    batch_results[i] = (read, false)   # skip the decode
                    skip_flags[i] = true
                else
                    improved_read,
                    was_improved = improve_read_likelihood(
                        read, graph, k; graph_mode = graph_mode,
                        beam_width = beam_width, weighted_graph = pass_weighted_graph,
                        diagnostics = diag)
                    batch_results[i] = (improved_read, was_improved)
                end
            end
            skipped_reads += count(skip_flags)

            # Collect results
            for (i, (improved_read, was_improved)) in enumerate(batch_results)
                updated_reads[batch_start + i - 1] = improved_read
                if was_improved
                    batch_improvements += 1
                end
            end
        else
            # Sequential processing (also the soft-EM path: accumulation into
            # `soft_weights` happens per-read inside improve_read_likelihood).
            for (i, read) in enumerate(batch_reads)
                if _skip_this_read_at(batch_start + i - 1)
                    updated_reads[batch_start + i - 1] = read   # skip the decode
                    skipped_reads += 1
                    continue
                end
                improved_read,
                was_improved = improve_read_likelihood(
                    read, graph, k; graph_mode = graph_mode,
                    beam_width = beam_width, soft_weights = soft_weights,
                    weighted_graph = pass_weighted_graph,
                    diagnostics = diag)
                updated_reads[batch_start + i - 1] = improved_read

                if was_improved
                    batch_improvements += 1
                end
            end
        end

        improvements_made += batch_improvements

        # Progress reporting and garbage collection
        if verbose
            batch_num = div(batch_start - 1, batch_size) + 1
            total_batches = div(total_reads - 1, batch_size) + 1
            improvement_rate = batch_improvements / length(batch_reads) * 100
            println("    Batch $batch_num/$total_batches: $(batch_improvements)/$(length(batch_reads)) improvements ($(round(improvement_rate, digits=1))%)")
        end

        # Force garbage collection between batches for memory efficiency
        if batch_end < total_reads
            GC.gc()
        end
    end

    skip_fraction = total_reads > 0 ? skipped_reads / total_reads : 0.0
    if verbose
        total_improvement_rate = improvements_made / total_reads * 100
        println("  Total improvements: $improvements_made/$total_reads ($(round(total_improvement_rate, digits=1))%)")
        if skip_solid || hard_vertices !== nothing || pass_decode_off
            println("  Skipped (no decode): $skipped_reads/$total_reads ($(round(skip_fraction * 100, digits=1))%)")
            decode_fraction = 1.0 - skip_fraction
            gated = pass_decode_off ? " [low-k decode gated OFF, td-9h5r]" : ""
            println("  Graph-Viterbi decode fraction: $(round(decode_fraction * 100, digits=1))%$(gated)")
        end
        if cheap_correct
            println("  Stage 0 cheap corrections (this pass): $cheap_corrections")
        end
    end

    # Surface swallowed decode failures for THIS pass (see CorrectorDiagnostics):
    # a systematically broken corrector otherwise reports 0 improvements as if it
    # ran. Distinguish "every read failed structurally" (the corrector is not
    # running at all) from a merely-high swallowed fraction.
    struct_swallowed = diag.structural_errors[] - struct_before
    unk_swallowed = diag.unkmerizable_reads[] - unk_before
    swallowed = struct_swallowed + unk_swallowed
    if total_reads > 0 && swallowed > 0
        frac = swallowed / total_reads
        if struct_swallowed == total_reads
            @warn "iterative corrector: EVERY read failed to decode with a structural " *
                  "error (empty graph / alphabet-inference miss / config error) — the " *
                  "corrector is not running; 0 improvements is NOT 'nothing to fix'." total_reads structural_errors = struct_swallowed
        elseif frac >= 0.5
            @warn "iterative corrector: high swallowed-decode fraction; many reads "*
            "passed through uncorrected (not conflated with 'no likelihood gain')." total_reads structural_errors=struct_swallowed unkmerizable_reads=unk_swallowed fraction=round(
                frac, digits = 3)
        end
    end

    return updated_reads, improvements_made, skip_fraction, cheap_corrections,
    pass_decode_off
end

"""
Improve likelihood of a single read using maximum likelihood path finding.
Returns improved read and boolean indicating if improvement was made.

`soft_weights` (td-e70t): when supplied, the read's decoded ML path accumulates
edge responsibilities into it (responsibility 1.0 for the single argmax path).
"""
function improve_read_likelihood(read::FASTX.FASTQ.Record, graph, k::Int;
        graph_mode::Symbol = :canonical,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        weighted_graph = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing)::Tuple{
        FASTX.FASTQ.Record, Bool}
    # Extract sequence and quality
    original_seq = FASTX.sequence(String, read)
    original_qual = FASTX.quality(read)
    read_id = FASTX.identifier(read)

    # Skip reads shorter than k
    if length(original_seq) < k
        return read, false
    end

    # Find optimal path through graph via the trustworthy Viterbi. `beam_width`
    # defaults to `nothing` = size-aware auto-beam (exact where tractable, bounded
    # above the threshold to avoid the td-63qy OOM); an explicit Int overrides it.
    #
    # Two-gate design (intentional, not a narrative/code mismatch): the inner
    # `find_optimal_sequence_path` gates on `improvement > 0.0` (accept the ML path
    # only when it is at least as likely as the original — never return a strictly
    # worse read), while THIS function additionally requires `> 0.01` so a
    # negligible likelihood wiggle is not counted as a real correction. The two
    # thresholds serve different jobs (correctness floor vs. meaningful-change
    # floor) and compose deliberately.
    improved_read,
    likelihood_improvement = find_optimal_sequence_path(
        read, graph, k; graph_mode = graph_mode,
        beam_width = beam_width, soft_weights = soft_weights,
        weighted_graph = weighted_graph,
        diagnostics = diagnostics)

    # Only update if significant improvement
    if likelihood_improvement > 0.01  # Threshold for meaningful improvement
        return improved_read, true
    else
        return read, false
    end
end

"""
Find optimal sequence path through graph using maximum likelihood principles.
Returns improved sequence and likelihood improvement score.
"""
function find_optimal_sequence_path(sequence::AbstractString, quality, graph, k::Int;
        graph_mode::Symbol = :canonical)::Tuple{String, Float64}
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for find_optimal_sequence_path"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    improved_record,
    improvement = find_optimal_sequence_path(record, graph, k; graph_mode = graph_mode)
    return FASTX.sequence(String, improved_record), improvement
end

function find_optimal_sequence_path(read::FASTX.FASTQ.Record, graph, k::Int;
        graph_mode::Symbol = :canonical,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        weighted_graph = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing)::Tuple{
        FASTX.FASTQ.Record, Float64}
    # Trustworthy Viterbi (graph-as-HMM correction core, td-ak6w). Trust the
    # maximum-likelihood path or report no improvement. The prior heuristic
    # fallback cascade (0.01 abandon-threshold -> statistical resampling -> local
    # edits) is removed from the correctness path: a "corrected" read must lie on
    # an ML graph path, not on an off-path heuristic that may not correspond to
    # any real traversal. `try_statistical_path_resampling` and
    # `try_local_path_improvements` remain DEFINED as opt-in comparison baselines
    # (and for direct callers/tests) but are no longer invoked here.
    #
    # `calculate_read_likelihood` k-merizes the read with the 2-bit iterator and
    # runs BEFORE the decoder, so an un-k-merizable read (ambiguous base, e.g. `N`)
    # throws `BioSequences.EncodeError` HERE, not in `try_viterbi_path_improvement`.
    # Treat it identically: SKIP the read (pass through uncorrected) + count, so a
    # single N-containing read never crashes the (threaded) corrector (td-63qy
    # family). Returning before the decoder means each such read is tallied once.
    original_likelihood = try
        calculate_read_likelihood(read, graph, k; graph_mode = graph_mode)
    catch error
        if error isa BioSequences.EncodeError
            diagnostics === nothing ||
                Threads.atomic_add!(diagnostics.unkmerizable_reads, 1)
            return read, 0.0
        end
        rethrow()
    end

    viterbi_result = try_viterbi_path_improvement(
        read, graph, k; graph_mode = graph_mode,
        beam_width = beam_width, soft_weights = soft_weights,
        weighted_graph = weighted_graph,
        diagnostics = diagnostics)
    if viterbi_result === nothing
        # Viterbi could not decode (empty observation set, structural error, or an
        # un-k-merizable read): report no improvement rather than falling back to a
        # heuristic path. Genuine decode FAILURES (vs. "no observations") are
        # tallied in `diagnostics` inside try_viterbi_path_improvement so they are
        # not silently conflated with "the ML path had no gain".
        return read, 0.0
    end

    viterbi_read, viterbi_likelihood = viterbi_result
    improvement = viterbi_likelihood - original_likelihood
    # Accept the ML path only when it is at least as likely as the original; a
    # non-positive delta means the ML path is not an improvement, so report none
    # (never return a strictly-worse read). This is the ADR's "trust the ML path
    # or report no improvement", NOT the removed 0.01 minimum-magnitude threshold.
    return improvement > 0.0 ? (viterbi_read, improvement) : (read, 0.0)
end

"""
Calculate likelihood of a FASTQ read given the current graph.
"""
function calculate_read_likelihood(
        read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol = :canonical)::Float64
    # Detect sequence type dynamically using existing alphabets.jl functions
    sequence_string = FASTX.sequence(String, read)
    alphabet = Mycelia.detect_alphabet(sequence_string)
    sequence_type = Mycelia.alphabet_to_biosequence_type(alphabet)

    sequence = Mycelia.extract_typed_sequence(read, sequence_type)
    quality_scores = collect(FASTX.quality_scores(read))
    return calculate_sequence_likelihood(
        sequence, quality_scores, graph, k; graph_mode = graph_mode)
end

"""
Calculate quality-aware likelihood of a sequence given the qualmer graph.
Uses both k-mer presence and quality-based confidence from qualmer observations.
"""
function calculate_sequence_likelihood(sequence::AbstractString, quality, graph, k::Int;
        graph_mode::Symbol = :canonical)::Float64
    alphabet = detect_alphabet(sequence)
    typed_sequence = detect_and_extract_sequence(sequence, alphabet)
    quality_scores = if quality isa AbstractString
        Int8.(Int.(collect(codeunits(quality))) .- 33)
    elseif quality isa AbstractVector{<:Integer}
        Int8.(quality)
    else
        throw(ArgumentError("Unsupported quality type for calculate_sequence_likelihood"))
    end
    return calculate_sequence_likelihood(
        typed_sequence, quality_scores, graph, k; graph_mode = graph_mode)
end

function _resolve_kmer_label(graph, kmer; graph_mode::Symbol = :canonical)
    if haskey(graph, kmer)
        return kmer
    end
    if graph_mode == :canonical && kmer isa Union{Kmers.DNAKmer, Kmers.RNAKmer}
        canonical_kmer = BioSequences.canonical(kmer)
        if haskey(graph, canonical_kmer)
            return canonical_kmer
        end
    end
    return nothing
end

function _resolve_qualmer_for_graph(graph, qmer::Qualmer; graph_mode::Symbol = :canonical)
    resolved_kmer = _resolve_kmer_label(graph, qmer.kmer; graph_mode = graph_mode)
    if resolved_kmer === nothing
        return nothing
    end
    if resolved_kmer == qmer.kmer
        return qmer
    end
    if graph_mode == :canonical && qmer.kmer isa Union{Kmers.DNAKmer, Kmers.RNAKmer}
        return canonical(qmer)
    end
    return Qualmer(resolved_kmer, qmer.qualities)
end

function calculate_sequence_likelihood(
        sequence::BioSequences.BioSequence, quality::Vector{Int8},
        graph, k::Int; graph_mode::Symbol = :canonical)::Float64
    if length(sequence) < k
        return 0.0
    end

    total_log_likelihood = 0.0

    # Use Kmers.jl iterators for efficient k-mer extraction with quality awareness
    if sequence isa BioSequences.LongDNA
        for (pos, kmer) in enumerate(Kmers.FwDNAMers{k}(sequence))
            resolved_kmer = _resolve_kmer_label(graph, kmer; graph_mode = graph_mode)

            if resolved_kmer !== nothing
                vertex_data = graph[resolved_kmer]

                # Enhanced quality-aware likelihood calculation
                kmer_likelihood = calculate_qualmer_likelihood(
                    resolved_kmer,
                    quality[pos:(pos + k - 1)],  # Quality scores for this k-mer
                    vertex_data
                )

                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                # Penalty for unseen k-mers, but consider quality
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:(pos + k - 1)])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    elseif sequence isa BioSequences.LongRNA
        for (pos, kmer) in enumerate(Kmers.FwRNAMers{k}(sequence))
            resolved_kmer = _resolve_kmer_label(graph, kmer; graph_mode = graph_mode)

            if resolved_kmer !== nothing
                vertex_data = graph[resolved_kmer]

                kmer_likelihood = calculate_qualmer_likelihood(
                    resolved_kmer,
                    quality[pos:(pos + k - 1)],
                    vertex_data
                )

                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:(pos + k - 1)])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    else  # LongAA - no canonical form needed
        for (pos, kmer) in enumerate(Kmers.FwAAMers{k}(sequence))
            if haskey(graph, kmer)
                vertex_data = graph[kmer]

                kmer_likelihood = calculate_qualmer_likelihood(
                    kmer,
                    quality[pos:(pos + k - 1)],
                    vertex_data
                )

                total_log_likelihood += log(max(kmer_likelihood, 1e-10))
            else
                unseen_penalty = calculate_unseen_kmer_penalty(quality[pos:(pos + k - 1)])
                total_log_likelihood += log(max(unseen_penalty, 1e-10))
            end
        end
    end

    return total_log_likelihood
end

function _rhizomorph_first_dataset_id(vertex_data)
    dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
    return isempty(dataset_ids) ? nothing : first(dataset_ids)
end

function _rhizomorph_joint_probability(vertex_data)::Float64
    dataset_id = _rhizomorph_first_dataset_id(vertex_data)
    if dataset_id === nothing
        return 0.0
    end

    joint_quality = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, dataset_id)
    if joint_quality === nothing
        return 0.0
    end

    log_prob = 0.0
    for q in joint_quality
        log_prob += log(max(phred_to_probability(q), 1e-10))
    end
    return exp(log_prob)
end

function _vertex_joint_probability(vertex_data)::Float64
    if hasproperty(vertex_data, :joint_probability)
        return vertex_data.joint_probability
    end
    return _rhizomorph_joint_probability(vertex_data)
end

function _vertex_mean_quality(vertex_data)::Float64
    if hasproperty(vertex_data, :mean_quality)
        return vertex_data.mean_quality
    end

    dataset_id = _rhizomorph_first_dataset_id(vertex_data)
    if dataset_id === nothing
        return 0.0
    end

    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, dataset_id)
    if mean_quality === nothing || isempty(mean_quality)
        return 0.0
    end

    return Statistics.mean(mean_quality)
end

function _vertex_coverage(vertex_data)::Int
    if hasproperty(vertex_data, :coverage)
        coverage = vertex_data.coverage
        return coverage isa Integer ? coverage : length(coverage)
    end
    return Int(Mycelia.Rhizomorph.count_evidence(vertex_data))
end

"""
Calculate likelihood of a k-mer given its observed quality scores and qualmer vertex data.
Leverages existing qualmer graph quality calculations.
"""
function calculate_qualmer_likelihood(kmer, quality_scores::AbstractVector{Int8}, vertex_data)::Float64
    # Base likelihood from graph joint probability
    base_likelihood = _vertex_joint_probability(vertex_data)

    # Quality-based confidence adjustment
    # Convert Int8 quality scores to UInt8 for compatibility with existing functions
    quality_uint8 = UInt8.(max.(0, quality_scores))

    # Calculate position-wise quality confidence
    quality_confidence = 1.0
    for q_score in quality_uint8
        prob_correct = phred_to_probability(q_score)
        quality_confidence *= prob_correct
    end

    # Combine graph confidence with observed quality confidence
    # Weight the graph probability by the quality confidence of this observation
    combined_likelihood = base_likelihood * quality_confidence +
                          (1.0 - quality_confidence) * (base_likelihood * 0.1)  # Lower confidence for low quality

    return max(combined_likelihood, 1e-10)
end

"""
Calculate penalty for unseen k-mers based on their quality scores.
High quality unseen k-mers get less penalty than low quality ones.
"""
function calculate_unseen_kmer_penalty(quality_scores::AbstractVector{Int8})::Float64
    quality_uint8 = UInt8.(max.(0, quality_scores))

    # Calculate average quality confidence
    avg_quality_confidence = Statistics.mean(phred_to_probability(q) for q in quality_uint8)

    # High quality unseen k-mers might be real but missing from training data
    # Low quality unseen k-mers are likely errors
    base_penalty = 1e-6
    quality_adjusted_penalty = base_penalty * avg_quality_confidence +
                               base_penalty * 0.01 * (1.0 - avg_quality_confidence)

    return quality_adjusted_penalty
end

"""
Generate alternative k-mers for improvement attempts using proper k-mer objects.
"""
function generate_kmer_alternatives(original_kmer::Kmers.Kmer, graph; graph_mode::Symbol = :canonical)
    sequence_type = if original_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif original_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else
        BioSequences.LongAA
    end
    return generate_kmer_alternatives(original_kmer, graph, sequence_type; graph_mode = graph_mode)
end

function generate_kmer_alternatives(original_kmer::AbstractString, graph; graph_mode::Symbol = :canonical)
    alphabet = detect_alphabet(original_kmer)
    sequence_type = alphabet_to_biosequence_type(alphabet)
    k = length(original_kmer)
    kmer_type = if sequence_type <: BioSequences.LongDNA
        Kmers.DNAKmer{k}
    elseif sequence_type <: BioSequences.LongRNA
        Kmers.RNAKmer{k}
    else
        Kmers.AAKmer{k}
    end
    return generate_kmer_alternatives(kmer_type(original_kmer), graph, sequence_type; graph_mode = graph_mode)
end

function generate_kmer_alternatives(
        original_kmer, graph, sequence_type::Type{<:BioSequences.BioSequence};
        graph_mode::Symbol = :canonical)
    alternatives = []
    k = length(original_kmer)

    # Get alphabet based on sequence type
    if sequence_type <: BioSequences.LongDNA
        alphabet = [
            BioSequences.DNA_A, BioSequences.DNA_T, BioSequences.DNA_G, BioSequences.DNA_C]
        kmer_type = Kmers.DNAKmer{k}
    elseif sequence_type <: BioSequences.LongRNA
        alphabet = [
            BioSequences.RNA_A, BioSequences.RNA_U, BioSequences.RNA_G, BioSequences.RNA_C]
        kmer_type = Kmers.RNAKmer{k}
    else  # LongAA
        # For amino acids, use the standard 20 amino acids
        alphabet = [
            BioSequences.AA_A, BioSequences.AA_R, BioSequences.AA_N, BioSequences.AA_D,
            BioSequences.AA_C, BioSequences.AA_Q, BioSequences.AA_E, BioSequences.AA_G,
            BioSequences.AA_H, BioSequences.AA_I, BioSequences.AA_L, BioSequences.AA_K,
            BioSequences.AA_M, BioSequences.AA_F, BioSequences.AA_P, BioSequences.AA_S,
            BioSequences.AA_T, BioSequences.AA_W, BioSequences.AA_Y, BioSequences.AA_V]
        kmer_type = Kmers.AAKmer{k}
    end

    # Try single position substitutions
    for i in 1:k
        original_symbol = original_kmer[i]
        for symbol in alphabet
            if symbol != original_symbol
                # Build alternative k-mer with substitution at position i
                alt_kmer_seq = sequence_type([original_kmer[j] for j in 1:k])
                alt_kmer_seq[i] = symbol
                alt_kmer = kmer_type(alt_kmer_seq)

                # Resolve k-mer labels against the current graph mode
                resolved_alt = _resolve_kmer_label(graph, alt_kmer; graph_mode = graph_mode)
                resolved_orig = _resolve_kmer_label(graph, original_kmer; graph_mode = graph_mode)

                # Only include if alternative exists in graph with higher probability
                if resolved_alt !== nothing && resolved_orig !== nothing
                    alt_data = graph[resolved_alt]
                    orig_data = graph[resolved_orig]

                    if _vertex_joint_probability(alt_data) >
                       _vertex_joint_probability(orig_data)
                        push!(alternatives, resolved_alt)
                    end
                end
            end
        end
    end

    return alternatives
end

"""
Generate alternative k-mer paths by substituting individual k-mers with plausible alternatives.
Returns a list of alternative paths, each represented as a vector.
"""
function generate_alternative_paths(sequence::AbstractString, graph, k::Int;
        num_samples::Int = 5, graph_mode::Symbol = :canonical)
    base_path = sequence_to_kmer_path(String(sequence), k)
    if isempty(base_path) || num_samples <= 0
        return Vector{Vector{Any}}()
    end

    alternative_paths = Vector{Vector{Any}}()
    for (idx, kmer) in enumerate(base_path)
        alternatives = generate_kmer_alternatives(kmer, graph; graph_mode = graph_mode)
        for alt in alternatives
            if alt != kmer
                alt_path = Any[k for k in base_path]
                alt_path[idx] = alt
                push!(alternative_paths, alt_path)
                if length(alternative_paths) >= num_samples
                    return alternative_paths
                end
                break
            end
        end
    end

    return alternative_paths
end

"""
Adjust quality scores based on likelihood improvement.
"""
function adjust_quality_scores(
        original_quality::AbstractString, improved_sequence::AbstractString,
        likelihood_improvement::Float64)::String
    original_quality_str = String(original_quality)
    improved_sequence_str = String(improved_sequence)

    # For now, return original quality scores
    # This can be enhanced to adjust qualities based on the improvement
    if length(improved_sequence_str) == length(original_quality_str)
        return original_quality_str
    else
        # Handle length changes by extending or truncating quality
        if length(improved_sequence_str) > length(original_quality_str)
            # Extend with median quality
            median_qual = isempty(original_quality_str) ? 'I' :
                          original_quality_str[length(original_quality_str) ÷ 2]
            return original_quality_str *
                   repeat(string(median_qual), length(improved_sequence_str) -
                                               length(original_quality_str))
        else
            # Truncate quality
            return original_quality_str[1:length(improved_sequence_str)]
        end
    end
end

# =============================================================================
# Decision Making and Termination Conditions
# =============================================================================

"""
Determine if sufficient improvements were made to continue with current k.
Enhanced with convergence detection and adaptive thresholds.
"""
function sufficient_improvements(
        improvements_made::Int, total_reads::Int, threshold::Float64 = 0.05;
        iteration_history::Vector{Dict{Symbol, Any}} = Dict{Symbol, Any}[],
        min_improvement_trend::Int = 3)::Bool
    if total_reads == 0
        return false
    end

    improvement_rate = improvements_made / total_reads

    # Basic threshold check
    if improvement_rate < threshold
        return false
    end

    # Early convergence detection based on improvement trend
    if length(iteration_history) >= min_improvement_trend
        recent_rates = [hist[:improvement_rate]
                        for hist in
                            iteration_history[(end - min_improvement_trend + 1):end]]

        # Check for diminishing returns (decreasing trend)
        if length(recent_rates) >= 2
            trend_decreasing = all(recent_rates[i] >= recent_rates[i + 1]
            for i in 1:(length(recent_rates) - 1))
            avg_recent_rate = Statistics.mean(recent_rates)

            # If trend is consistently decreasing and below adaptive threshold, stop
            if trend_decreasing && avg_recent_rate < threshold * 2.0
                return false
            end
        end

        # Check for plateau (very small improvements)
        if all(rate < threshold * 0.5 for rate in recent_rates)
            return false
        end
    end

    return true
end

"""
Enhanced convergence detection for k-mer progression.
Determines if we should move to the next k-mer size based on multiple criteria.
"""
function should_continue_k_progression(k_history::Vector{Dict{Symbol, Any}},
        current_k::Int,
        max_k::Int;
        convergence_window::Int = 3,
        quality_improvement_threshold::Float64 = 0.001)::Bool
    if current_k >= max_k
        return false
    end

    if length(k_history) < convergence_window
        return true
    end

    # Check recent quality improvements across iterations
    recent_history = k_history[(end - convergence_window + 1):end]

    # Calculate average improvement metrics
    avg_improvement_rate = Statistics.mean([hist[:total_improvements] / hist[:total_reads]
                                            for hist in recent_history])
    avg_runtime_per_read = Statistics.mean([hist[:runtime_seconds] / hist[:total_reads]
                                            for hist in recent_history])

    # Convergence criteria
    improvements_declining = avg_improvement_rate < quality_improvement_threshold
    runtime_increasing = avg_runtime_per_read > 0.1  # More than 0.1 seconds per read indicates inefficiency

    # If improvements are minimal and runtime is high, move to next k
    if improvements_declining && runtime_increasing
        return false
    end

    # If we're making good progress, continue
    return avg_improvement_rate >= quality_improvement_threshold
end

"""
Finalize iterative assembly by combining results from all k-mer sizes and iterations.
"""
function finalize_iterative_assembly(output_dir::String, k_progression::Vector{Int},
        iteration_history::Dict{Int, Vector{Dict{Symbol, Any}}},
        total_runtime::Float64; verbose::Bool = true,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing,
        skip_fractions::Vector{Float64} = Float64[],
        cheap_correction_counts::Vector{Int} = Int[],
        hard_window::Bool = false,
        soft_em::Bool = false,
        cheap_correct::Bool = false,
        min_decode_k::Union{Int, Nothing} = nothing,
        decode_gate_density::Union{Float64, Nothing} = nothing,
        decode_gated_rungs::Vector{Int} = Int[],
        final_pass_graph = nothing,
        final_pass_graph_k::Int = 0,
        final_pass_graph_mode::Symbol = :canonical,
        final_pass_graph_reusable::Bool = false)
    if verbose
        println("Finalizing iterative assembly results...")
    end

    # Only retain the (potentially large) final-pass graph when it is byte-identical-
    # reusable (td-04tb); otherwise drop the reference so a fallback rebuild does not
    # double-hold graph memory.
    if !final_pass_graph_reusable
        final_pass_graph = nothing
    end

    # Find the final output file (last iteration of largest k)
    final_k = maximum(k_progression)
    final_iteration_count = length(iteration_history[final_k])

    # Get the most recent timestamp from final k
    final_timestamp = iteration_history[final_k][end][:timestamp]
    final_fastq = joinpath(
        output_dir, "reads_k$(final_k)_iter$(final_iteration_count)_$(final_timestamp).fastq")

    # Calculate summary statistics
    total_iterations = sum(length(iterations) for iterations in values(iteration_history))
    total_improvements = sum(sum(iter[:improvements_made] for iter in iterations)
    for iterations in values(iteration_history))

    # Read final assembly for k-mer extraction
    if isfile(final_fastq)
        final_reads = collect(FASTX.FASTQ.Reader(open(final_fastq)))
        final_assembly = [FASTX.sequence(String, read) for read in final_reads]
    else
        final_assembly = String[]
    end

    # Swallowed-decode tally (structural / un-k-merizable). Surfaced on the result
    # so a caller can detect a corrector that reported 0 improvements because it
    # never actually ran, distinct from "ran and found nothing to fix".
    corrector_errors = diagnostics === nothing ?
                       Dict(:structural => 0, :unkmerizable => 0) :
                       Dict(:structural => diagnostics.structural_errors[],
        :unkmerizable => diagnostics.unkmerizable_reads[])

    # Per-pass skip-fraction summary (FIX 6): min/mean/max across every k+iteration
    # pass, not just the last. `last_skip_fraction` is retained for back-compat.
    last_skip_fraction = isempty(skip_fractions) ? 0.0 : skip_fractions[end]
    skip_min = isempty(skip_fractions) ? 0.0 : minimum(skip_fractions)
    skip_max = isempty(skip_fractions) ? 0.0 : maximum(skip_fractions)
    skip_mean = isempty(skip_fractions) ? 0.0 : sum(skip_fractions) / length(skip_fractions)

    # Graph-Viterbi decode fraction per pass = 1 - skip_fraction (the reads that
    # actually reached the expensive decode). This is the critical-path metric:
    # narrowing the hard-vertex set + Stage 0 cheap correction should drive it down
    # toward the true bubble/repeat ~5-15% (td-bjnt). Zero on tiers without a gate.
    decode_fractions = [1.0 - s for s in skip_fractions]
    decode_min = isempty(decode_fractions) ? 0.0 : minimum(decode_fractions)
    decode_max = isempty(decode_fractions) ? 0.0 : maximum(decode_fractions)
    decode_mean = isempty(decode_fractions) ? 0.0 :
                  sum(decode_fractions) / length(decode_fractions)

    # Stage 0 cheap-correction totals (td-bjnt).
    total_cheap_corrections = isempty(cheap_correction_counts) ? 0 :
                              sum(cheap_correction_counts)

    # Create comprehensive metadata
    metadata = Dict(
        :total_runtime => total_runtime,
        :total_iterations => total_iterations,
        :total_improvements => total_improvements,
        :k_progression => k_progression,
        :final_k => final_k,
        :final_fastq_file => final_fastq,
        :output_directory => output_dir,
        :iteration_history => iteration_history,
        :corrector_errors => corrector_errors,
        :assembly_type => "iterative_maximum_likelihood",
        :version => "Phase_5.2a",
        # Scalable-tier telemetry (td-fuo8/td-nn6l/td-e70t). Zero/false on the
        # exhaustive tier, which threads neither gate.
        #
        # Honest gate provenance (FIX 5): `hard_read_gate` is the SKIP gate
        # (easy reads pass through untouched) — that IS active on :scalable.
        # `windowed_decode` is the per-hard-region WINDOWED decode, which is
        # scaffolded (hard reads are decoded WHOLE), so it is always false in v1;
        # surfaced separately so `hard_window` can't be misread as "windowed
        # decode active". `hard_window` is kept as a back-compat alias.
        :hard_window => hard_window,
        :hard_read_gate => hard_window,
        :windowed_decode => false,
        # Soft-EM v2 ACTIVE: the E-step enumerates competing paths and the M-step
        # registers the support-floored soft edge weights onto the next EM
        # iteration's graph, so unsupported error edges decay while supported
        # variation is retained. `false` on the exhaustive tier (soft-EM off there).
        :soft_em => (soft_em ? "v2-competing-paths-floor" : false),
        # Low-k decode gating (td-9h5r): the explicit floor (if any), the adaptive
        # density threshold, and the specific rungs whose per-read decode was gated
        # OFF (explicit floor OR adaptive). `nothing`/empty on the exhaustive tier.
        :min_decode_k => min_decode_k,
        :decode_gate_density => decode_gate_density,
        :decode_gated_rungs => decode_gated_rungs,
        :last_skip_fraction => last_skip_fraction,
        :skip_fraction_per_pass => skip_fractions,
        :skip_fraction_min => skip_min,
        :skip_fraction_mean => skip_mean,
        :skip_fraction_max => skip_max,
        # Graph-Viterbi decode fraction (= 1 - skip_fraction). The critical-path
        # win metric (td-bjnt): after Stage 0 + hard-set narrowing this drops
        # toward the true bubble/repeat ~5-15%.
        :decode_fraction_per_pass => decode_fractions,
        :decode_fraction_min => decode_min,
        :decode_fraction_mean => decode_mean,
        :decode_fraction_max => decode_max,
        # Stage 0 cheap k-mer-spectrum correction (td-bjnt): bases fixed by the
        # linear pass, per pass and in total. Zero/false on :exhaustive.
        :cheap_correct => cheap_correct,
        :cheap_corrections_per_pass => cheap_correction_counts,
        :cheap_corrections_total => total_cheap_corrections,
        # Final-pass graph reuse provenance (td-04tb). `final_graph_reusable` is true
        # only when the final pass at the largest k made 0 improvements, so the
        # returned `:final_graph` is byte-identical to a from-scratch rebuild of the
        # corrected reads at `final_graph_k` under `final_graph_mode`. The caller
        # verifies k + mode match its re-assembly parameters before reusing.
        :final_graph_k => final_pass_graph_k,
        :final_graph_mode => final_pass_graph_mode,
        :final_graph_reusable => final_pass_graph_reusable
    )

    if verbose
        println("Iterative assembly finalized:")
        println("  Total k-mer sizes: $(length(k_progression))")
        println("  Total iterations: $total_iterations")
        println("  Total improvements: $total_improvements")
        println("  Final FASTQ: $final_fastq")
    end

    return Dict(
        :final_assembly => final_assembly,
        :k_progression => k_progression,
        :metadata => metadata,
        # The corrector's final-pass qualmer graph (td-04tb). `nothing` when no pass
        # completed; only byte-identical-reusable when `metadata[:final_graph_reusable]`.
        :final_graph => final_pass_graph
    )
end

# =============================================================================
# Utility and Integration Functions
# =============================================================================

"""
Generate summary report of iterative assembly process.
"""
function iterative_assembly_summary(result::Dict)::String
    metadata = result[:metadata]
    k_prog = metadata[:k_progression]

    report = "Mycelia Iterative Assembly Summary\n"
    report *= "=" * "="^50 * "\n\n"

    report *= "Assembly Type: $(metadata[:assembly_type])\n"
    report *= "Version: $(metadata[:version])\n"
    report *= "Total Runtime: $(round(metadata[:total_runtime], digits=2)) seconds\n"
    report *= "K-mer Progression: $(sort(k_prog))\n"
    report *= "Total Iterations: $(metadata[:total_iterations])\n"
    report *= "Total Improvements: $(metadata[:total_improvements])\n"
    report *= "Final FASTQ: $(metadata[:final_fastq_file])\n\n"

    # Per-k statistics
    report *= "Per-K Statistics:\n"
    for k in sort(k_prog)
        iterations = metadata[:iteration_history][k]
        total_improvements_k = sum(iter[:improvements_made] for iter in iterations)
        report *= "  k=$k: $(length(iterations)) iterations, $total_improvements_k improvements\n"
    end

    return report
end

"""
Run a minimal self-test of the iterative assembly pipeline.
Returns a status dictionary with results or error details.
"""
function test_iterative_assembly(; max_k::Int = 13, max_iterations_per_k::Int = 1, cleanup::Bool = true)::Dict
    temp_dir = mktempdir()
    try
        test_fastq = joinpath(temp_dir, "test_reads.fastq")
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "CGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATC",
            "ATCGATCGATCGTTCGATCG"
        ]
        records = [FASTX.FASTQ.Record("read$i", seq, repeat("I", length(seq)))
                   for (i, seq) in enumerate(sequences)]
        write_fastq(records = records, filename = test_fastq)

        output_dir = joinpath(temp_dir, "test_output")
        result = mycelia_iterative_assemble(
            test_fastq;
            max_k = max_k,
            output_dir = output_dir,
            max_iterations_per_k = max_iterations_per_k,
            verbose = false
        )

        return Dict(
            :status => :success,
            :final_assembly => result[:final_assembly],
            :k_progression => result[:k_progression],
            :metadata => result[:metadata]
        )
    catch e
        return Dict(:status => :error, :error => sprint(showerror, e))
    finally
        if cleanup
            try
                rm(temp_dir, recursive = true, force = true)
            catch
            end
        end
    end
end

# =============================================================================
# Enhanced Statistical Path Improvement (Phase 5.2b)
# Integration with Viterbi algorithms and statistical path resampling
# =============================================================================

"""
Try to improve FASTQ read using Viterbi algorithm from viterbi-next.jl.
Returns (improved_read, likelihood) or nothing if no improvement.
"""
function try_viterbi_path_improvement(read::FASTX.FASTQ.Record,
        graph,
        k::Int;
        graph_mode::Symbol = :canonical,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        weighted_graph = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing)::Union{
        Tuple{FASTX.FASTQ.Record, Float64}, Nothing}
    try
        sequence_string = FASTX.sequence(String, read)
        alphabet = Mycelia.detect_alphabet(sequence_string)
        sequence_type = Mycelia.alphabet_to_biosequence_type(alphabet)
        sequence = Mycelia.extract_typed_sequence(read, sequence_type)
        # NOTE: `_record_kmer_iterator` uses a 2-bit (unambiguous) k-mer iterator,
        # so `collect` here throws `BioSequences.EncodeError` on a read carrying an
        # ambiguous base (e.g. a single `N`). That is handled as a SKIP-with-count
        # in the catch below, NOT a crash of the (threaded) batch (td-63qy family).
        kmers = collect(Mycelia._record_kmer_iterator(sequence_type, k, sequence))
        if isempty(kmers)
            return nothing
        end

        # Restore the emission model (graph-as-HMM correction core, td-jqdf): the
        # read's per-base Phred quality is the emission P(observed | hidden). The
        # k-mer at index `i` spans read positions `i:(i+k-1)`, so wrap it with that
        # per-base window (clamped to bounds) instead of discarding the quality and
        # letting the emission fall back to the graph's population-average quality.
        quality_scores = collect(FASTX.quality_scores(read))
        n_qual = length(quality_scores)
        observations = Vector{Mycelia.QualityObservation}(undef, length(kmers))
        for (i, kmer) in enumerate(kmers)
            lo = clamp(i, 1, n_qual)
            hi = clamp(i + k - 1, 1, n_qual)
            window = UInt8.(@view quality_scores[lo:hi])
            observations[i] = Mycelia.QualityObservation(kmer, window)
        end

        # Trustworthy Viterbi (td-ak6w), reconciled with the td-63qy OOM crash and
        # the td-35ux short-read O(genome^2) frontier blowup: the decoder is EXACT
        # where tractable and BOUNDED (with a log line) otherwise. `beam_width ===
        # nothing` (the default) resolves the frontier bound PER READ via the
        # GRAPH-DENSITY-AWARE `_auto_beam_width(length(observations), Graphs.nv(graph))`:
        # `typemax(Int)` (never pruned → global maximum-likelihood path) only when
        # BOTH the read is small AND the graph is sparse; the proven-tractable
        # finite bound when the read is large (read-length rule, td-63qy) OR the
        # graph is dense (density rule, td-35ux — a short read on a big graph
        # otherwise retains an O(n_vertices) frontier at each decode depth). The
        # graph threaded into this call is the one being decoded, so its vertex
        # count is exactly the density the frontier can reach. An explicit
        # `beam_width::Int` (including `typemax(Int)`) is an OPT-IN override honored
        # verbatim — a caller/test can force exact decoding on every read. This
        # trades the ML guarantee for tractability only under the auto-default; it
        # never silently changes correctness of an explicit request.
        effective_beam_width = beam_width === nothing ?
                               _auto_beam_width(length(observations), Graphs.nv(graph)) :
                               beam_width
        # Candidate-generation bounds (td-plqi) engage ONLY where the width beam is
        # already finite (i.e. the decode is already the bounded approximation, not
        # exact ML). When the beam is exact (typemax — small read on a sparse graph,
        # OR an explicit caller override forcing exact), the score margin stays Inf
        # and the successor bound stays unbounded, so exact-ML reads are byte-
        # identical. Where the beam is bounded (dense/large — the k=9-style rungs
        # that carry the #386 super-linear residual) the margin holds the generating
        # frontier O(1) in genome size.
        beam_is_exact = effective_beam_width == typemax(Int)
        effective_margin = beam_is_exact ? Inf : _AUTO_BEAM_SCORE_MARGIN
        effective_successor_bound = beam_is_exact ? typemax(Int) : _AUTO_SUCCESSOR_BOUND
        config = Mycelia.ViterbiCorrectionConfig(
            alphabet = alphabet,
            strand_mode = graph_mode,
            max_steps = length(observations) - 1,
            beam_width = effective_beam_width,
            max_successors_per_state = effective_successor_bound,
            beam_score_margin = effective_margin
        )
        correction = Mycelia.correct_observations(
            graph, [observations]; config = config, weighted_graph = weighted_graph)
        corrected_path = only(correction.corrected_observations)
        if corrected_path === nothing || isempty(corrected_path)
            return nothing
        end

        # Soft-EM E-step (td-e70t v2): build this read's COMPETING candidate paths
        # (the observed read path + a consensus alternative re-routed through the
        # best-supported sibling branch) and split the responsibility across them
        # by a stable softmax over their normalized-transition log-probabilities,
        # accumulating each path's edges weighted by its share. A data-supported
        # branch keeps high mass; an unsupported (error) branch accrues little and
        # decays across EM iterations (the M-step registers this accumulator —
        # support-floored — onto the next graph). Additive + guarded so it cannot
        # perturb the :exhaustive path (soft_weights === nothing there); degenerates
        # to the single observed path at responsibility 1.0 when no distinct
        # alternative exists (a balanced variant, whose branches are retained
        # through their own reads).
        if soft_weights !== nothing
            accumulate_competing_paths!(
                soft_weights, read, graph, k; graph_mode = graph_mode)
        end

        corrected_sequence = Mycelia.Rhizomorph.path_to_sequence(corrected_path, graph)
        corrected_sequence_string = corrected_sequence isa AbstractString ?
                                    corrected_sequence :
                                    string(corrected_sequence)
        improved_likelihood = Mycelia.calculate_sequence_likelihood(
            corrected_sequence_string, quality_scores, graph, k; graph_mode = graph_mode)
        improved_quality = Mycelia.adjust_quality_scores(
            FASTX.quality(read), corrected_sequence_string, improved_likelihood)
        improved_record = FASTX.FASTQ.Record(
            FASTX.identifier(read), corrected_sequence_string, improved_quality)
        return (improved_record, improved_likelihood)
    catch error
        if error isa ArgumentError
            # STRUCTURAL decode failure (empty graph, alphabet-inference miss,
            # config error). Swallow as "no improvement" so one broken read does
            # not abort the (threaded) batch — but COUNT it so a run that silently
            # corrected nothing is visible (not conflated with "decoder ran, no
            # gain"). See CorrectorDiagnostics + the per-pass @warn.
            diagnostics === nothing ||
                Threads.atomic_add!(diagnostics.structural_errors, 1)
            return nothing
        elseif error isa BioSequences.EncodeError
            # UN-K-MERIZABLE read: the 2-bit k-mer iterator throws EncodeError (NOT
            # ArgumentError) on an ambiguous base (e.g. a single `N` in an Illumina
            # read). Treat as a SKIP — the read passes through uncorrected — and
            # count it, rather than rethrowing and crashing the whole threaded
            # corrector on one N-containing read (td-63qy family).
            diagnostics === nothing ||
                Threads.atomic_add!(diagnostics.unkmerizable_reads, 1)
            return nothing
        end
        rethrow()
    end
end

function try_viterbi_path_improvement(
        sequence::AbstractString, quality, graph, k::Int; graph_mode::Symbol = :canonical)
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_viterbi_path_improvement"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    result = try_viterbi_path_improvement(record, graph, k; graph_mode = graph_mode)
    if result === nothing
        return nothing
    end

    improved_record, improvement = result
    return FASTX.sequence(String, improved_record), improvement
end

"""
Try statistical path resampling for alternative high-likelihood paths.
Returns (improved_read, likelihood) or nothing if no improvement.
"""
function try_statistical_path_resampling(read::FASTX.FASTQ.Record,
        graph,
        k::Int;
        graph_mode::Symbol = :canonical)::Union{Tuple{FASTX.FASTQ.Record, Float64}, Nothing}
    try
        # Detect sequence type dynamically
        sequence_string = FASTX.sequence(String, read)
        alphabet = detect_alphabet(sequence_string)
        sequence_type = alphabet_to_biosequence_type(alphabet)

        sequence = extract_typed_sequence(read, sequence_type)
        quality_scores = collect(FASTX.quality_scores(read))
        quality_string = FASTX.quality(read)

        # Generate multiple alternative paths using qualmer sampling
        alternative_paths = generate_alternative_qualmer_paths(
            read, graph, k, num_samples = 5, graph_mode = graph_mode)

        best_sequence = sequence
        best_likelihood = calculate_sequence_likelihood(
            sequence, quality_scores, graph, k; graph_mode = graph_mode)
        best_quality_scores = quality_scores

        # Evaluate each alternative path
        for alt_path in alternative_paths
            if !isempty(alt_path)
                # Convert qualmer path to proper BioSequence + quality
                alt_result = qualmer_path_to_biosequence(alt_path)
                alt_sequence = alt_result.sequence
                alt_quality_scores = alt_result.quality_scores

                alt_likelihood = calculate_sequence_likelihood(
                    alt_sequence, alt_quality_scores, graph, k; graph_mode = graph_mode)

                if alt_likelihood > best_likelihood
                    best_sequence = alt_sequence
                    best_likelihood = alt_likelihood
                    best_quality_scores = alt_quality_scores
                end
            end
        end

        # Return improvement if found
        if best_sequence != sequence
            # Create improved FASTQ record using proper types
            improved_sequence_string = string(best_sequence)
            improved_quality_string = String([Char(Int(q) + 33)
                                              for q in best_quality_scores])  # Convert back to ASCII
            improved_record = FASTX.FASTQ.Record(
                FASTX.identifier(read), improved_sequence_string, improved_quality_string)
            return (improved_record, best_likelihood)
        end

    catch e
        # If statistical resampling fails, return nothing
        return nothing
    end

    return nothing
end

function try_statistical_path_resampling(
        sequence::AbstractString, quality, graph, k::Int; graph_mode::Symbol = :canonical)
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_statistical_path_resampling"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    result = try_statistical_path_resampling(record, graph, k; graph_mode = graph_mode)
    if result === nothing
        return nothing
    end

    improved_record, improvement = result
    return FASTX.sequence(String, improved_record), improvement
end

"""
Local heuristic improvements (fallback method).
Returns (improved_read, likelihood_improvement).
"""
function try_local_path_improvements(
        read::FASTX.FASTQ.Record, graph, k::Int, original_likelihood::Float64;
        graph_mode::Symbol = :canonical)::Tuple{FASTX.FASTQ.Record, Float64}
    sequence_string = FASTX.sequence(String, read)
    alphabet = detect_alphabet(sequence_string)
    sequence_type = alphabet_to_biosequence_type(alphabet)
    quality_scores = collect(FASTX.quality_scores(read))

    improved_seq = sequence_string  # Start with original
    best_likelihood = original_likelihood

    # Try local improvements at positions with low-probability k-mers
    for i in 1:(length(sequence_string) - k + 1)
        kmer_str = sequence_string[i:(i + k - 1)]

        raw_kmer = if sequence_type <: BioSequences.LongDNA
            Kmers.DNAKmer{k}(kmer_str)
        elseif sequence_type <: BioSequences.LongRNA
            Kmers.RNAKmer{k}(kmer_str)
        else
            Kmers.AAKmer{k}(kmer_str)
        end

        target_kmer = _resolve_kmer_label(graph, raw_kmer; graph_mode = graph_mode)
        if target_kmer !== nothing
            vertex_data = graph[target_kmer]
            if _vertex_joint_probability(vertex_data) < 0.7  # Low confidence k-mer
                # Try improving this position
                alternatives = generate_kmer_alternatives(
                    target_kmer, graph, sequence_type; graph_mode = graph_mode)
                for alt_kmer in alternatives
                    alt_kmer_str = string(alt_kmer)
                    if i == 1
                        candidate_seq = alt_kmer_str * sequence_string[(i + k):end]
                    elseif i+k-1 == length(sequence_string)
                        candidate_seq = sequence_string[1:(i - 1)] * alt_kmer_str
                    else
                        candidate_seq = sequence_string[1:(i - 1)] * alt_kmer_str *
                                        sequence_string[(i + k):end]
                    end

                    candidate_bioseq = detect_and_extract_sequence(candidate_seq, alphabet)
                    candidate_likelihood = calculate_sequence_likelihood(
                        candidate_bioseq, quality_scores, graph, k; graph_mode = graph_mode)

                    if candidate_likelihood > best_likelihood
                        improved_seq = candidate_seq
                        best_likelihood = candidate_likelihood
                    end
                end
            end
        end
    end

    likelihood_improvement = best_likelihood - original_likelihood
    improved_quality = adjust_quality_scores(FASTX.quality(read), improved_seq, likelihood_improvement)
    improved_record = FASTX.FASTQ.Record(FASTX.identifier(read), improved_seq, improved_quality)
    return improved_record, likelihood_improvement
end

function try_local_path_improvements(
        sequence::AbstractString, quality, graph, k::Int, original_likelihood::Float64;
        graph_mode::Symbol = :canonical)::Tuple{String, Float64}
    quality_string = if quality isa AbstractString
        String(quality)
    elseif quality isa AbstractVector{<:Integer}
        String([Char(Int(q) + 33) for q in quality])
    else
        throw(ArgumentError("Unsupported quality type for try_local_path_improvements"))
    end

    record = FASTX.FASTQ.Record("input_sequence", sequence, quality_string)
    improved_record,
    improvement = try_local_path_improvements(
        record, graph, k, original_likelihood; graph_mode = graph_mode)
    return FASTX.sequence(String, improved_record), improvement
end

# ============================================================================
# Soft-EM E-step + M-step accumulation hook (ACTIVE, td-e70t v2 + support floor)
# ============================================================================
#
# See the design note over `SoftEdgeWeightAccumulator` in
# `src/rhizomorph/core/evidence-functions.jl`. This is the correction-side hook
# that turns a read's competing decoded paths into probability-weighted edge
# evidence: each candidate path's RESPONSIBILITY (posterior over the read's own
# candidate set) is accumulated onto the edges it traverses. Summed over reads,
# the accumulator holds soft edge weights in which error edges — traversed only
# by rare, low-probability paths — accrue little weight; once the M-step registers
# them, `compute_edge_weight` (and thus the Viterbi transition scoring) sees the
# decayed weight and the next iteration's decode + rebuild drops the error edge.
# NOTE the decay acts through `compute_edge_weight` / the Viterbi transition
# score, NOT through the vertex `weight > 0.01` gate in
# `generate_alternative_qualmer_paths` (that gate filters vertex candidates).
#
# v2 status — E-step AND M-step both WIRED. Under `soft_em=true`:
#   * E-step: `accumulate_competing_paths!` (below) builds each decoded read's
#     COMPETING candidate paths — the observed read path and a consensus
#     alternative re-routed through the best-supported sibling branch — then
#     softmax-splits the responsibility across them by normalized-transition
#     log-probability, so a real branch keeps high mass and an unsupported (error)
#     branch gets little.
#   * M-step: `mycelia_iterative_assemble` registers the accumulator onto the NEXT
#     EM iteration's graph (`register_soft_edge_weights!`), so `compute_edge_weight`
#     returns the decayed soft weight and the next iteration's transition scoring
#     is biased away from the error edge — a feedback loop that drives the error
#     edge's weight strictly down.
# VARIATION PRESERVATION (td-h6w9): the M-step registration applies a SUPPORT
# FLOOR — an edge backed by >= `SOFT_EM_MIN_SUPPORT` reads is clamped to at least
# its raw coverage, so a real but SKEWED minority allele (e.g. a 10x branch in a
# 20x/10x bubble) NEVER decays toward zero regardless of a heavier sibling. Only
# near-zero-support (error) edges are free to decay below the cleaning gate. This
# fixes the prior v2's collapse of skewed variants (the responsibility split alone
# gave `W_min' = N*W/(W_maj+W_min)`, a geometric decay to zero for any imbalance).

"""
    _decoded_path_edges(result) -> Vector

Ordered `(src_label, dst_label)` edge tuples of a decoded Viterbi path, or an
empty vector when the result carries no path.
"""
function _decoded_path_edges(result::Mycelia.Rhizomorph.ViterbiDecodingResult)
    labels = _decoded_path_labels(result)
    labels === nothing && return Tuple{Any, Any}[]
    return Tuple{Any, Any}[(labels[i], labels[i + 1]) for i in 1:(length(labels) - 1)]
end

"""
    accumulate_soft_em_edge_weights!(accumulator, candidate_paths) -> accumulator

Soft-EM E-step core (Viterbi-result flavor). Given a single read's competing
decoded paths (`ViterbiDecodingResult`s, each with a `.score` log-likelihood and
a `.path`), compute each path's responsibility by softmax-normalizing its score
against the candidate set (`Mycelia.Rhizomorph.path_responsibility`) and
accumulate that responsibility onto the edges it traverses
(`Mycelia.Rhizomorph.accumulate_path_probability!`).

With a single decoded path per read the responsibility is `1.0` (soft weight ==
hard count); the softening comes from `accumulate_competing_paths!`, which
supplies several candidate paths per read. Retained as a primitive (used by the
unit tests) alongside the graph-level `accumulate_competing_paths!` used by the
pipeline E-step.
"""
function accumulate_soft_em_edge_weights!(
        accumulator::Mycelia.Rhizomorph.SoftEdgeWeightAccumulator,
        candidate_paths
)
    scores = Float64[result.score for result in candidate_paths if isfinite(result.score)]
    isempty(scores) && return accumulator
    for result in candidate_paths
        isfinite(result.score) || continue
        responsibility = Mycelia.Rhizomorph.path_responsibility(result.score, scores)
        Mycelia.Rhizomorph.accumulate_path_probability!(
            accumulator, _decoded_path_edges(result), responsibility)
    end
    return accumulator
end

# ----------------------------------------------------------------------------
# Soft-EM v2 competing-paths E-step (td-e70t): graph-level candidate enumeration
# ----------------------------------------------------------------------------

"""
    _graph_oriented_edges(labels, graph) -> Vector{Tuple}

Turn an ordered list of vertex `labels` into the `(src, dst)` edge tuples in the
ORIENTATION the graph actually stores them (checking both directions), skipping
any consecutive pair with no edge in `graph`. Keying the accumulator in the
graph's own orientation is what lets `register_soft_edge_weights!` — which
iterates `MetaGraphsNext.edge_labels(graph)` — find and consume the soft weight.
A candidate substitution that is not a real graph traversal contributes no edge.
"""
function _graph_oriented_edges(labels, graph)
    edges = Tuple{Any, Any}[]
    for i in 1:(length(labels) - 1)
        a, b = labels[i], labels[i + 1]
        if MetaGraphsNext.haskey(graph, a, b)
            push!(edges, (a, b))
        elseif MetaGraphsNext.haskey(graph, b, a)
            push!(edges, (b, a))
        end
    end
    return edges
end

# The observed read path, resolved to graph vertex labels (canonical where the
# graph is canonical). Empty when no k-mer of the read resolves onto the graph.
function _read_resolved_labels(read::FASTX.FASTQ.Record, graph, k::Int; graph_mode::Symbol)
    labels = Any[]
    for (qmer, _pos) in qualmers_unambiguous(read, k)
        resolved = _resolve_qualmer_for_graph(graph, qmer; graph_mode = graph_mode)
        resolved === nothing || push!(labels, resolved.kmer)
    end
    return labels
end

# Soft-weight-aware edge weight between two vertex labels in either stored
# orientation, or `nothing` when they are not adjacent. `compute_edge_weight`
# consults the soft-edge registry, so an edge the M-step decayed reads back its
# decayed weight here — the cross-iteration feedback that sharpens the split.
function _edge_weight_between(graph, a, b)
    if MetaGraphsNext.haskey(graph, a, b)
        return Float64(Mycelia.Rhizomorph.compute_edge_weight(graph[a, b]))
    elseif MetaGraphsNext.haskey(graph, b, a)
        return Float64(Mycelia.Rhizomorph.compute_edge_weight(graph[b, a]))
    end
    return nothing
end

# All labels adjacent to `a` (either orientation), deduplicated. Used both as the
# transition normalizer support and as the branch-candidate set.
function _incident_labels(graph, a)
    ns = Any[]
    for n in MetaGraphsNext.outneighbor_labels(graph, a)
        push!(ns, n)
    end
    for n in MetaGraphsNext.inneighbor_labels(graph, a)
        push!(ns, n)
    end
    return unique(ns)
end

# Sum of the soft-weight-aware edge weights on every edge incident to `a` — the
# denominator that turns a raw edge weight into a normalized transition
# probability. Normalizing here counts a single branch decision ONCE (internal
# edges of a linear run normalize to ~0.5 on BOTH competing branches and cancel
# in the responsibility softmax), so the divergence edge carries the coverage
# contrast and a real but skewed variant is not collapsed like a raw
# coverage-ratio prune would.
function _incident_weight_sum(graph, a)
    total = 0.0
    for x in _incident_labels(graph, a)
        w = _edge_weight_between(graph, a, x)
        w === nothing || (total += w)
    end
    return total
end

# Normalized-transition log-probability of a label path: sum of log(w(a,b)/S(a)).
# Returns -Inf if any consecutive pair is non-adjacent or zero-weight, so an
# invalid candidate is dropped by the softmax `isfinite` guard.
function _path_transition_logscore(labels, graph)
    total = 0.0
    for i in 2:length(labels)
        a, b = labels[i - 1], labels[i]
        w = _edge_weight_between(graph, a, b)
        (w === nothing || w <= 0.0) && return -Inf
        s = _incident_weight_sum(graph, a)
        s <= 0.0 && return -Inf
        total += log(w / s)
    end
    return total
end

# Generate ONE competing alternative to `observed` by re-routing its first
# clearly-weak branch through the highest-weight sibling and walking greedily
# until rejoining `observed`. Returns the spliced label path, or `nothing` when
# no weak branch exists (a linear region, or a balanced variant whose sibling is
# not strictly better — which is retained, not competed away). Bounded by the
# observed length so it always terminates.
function _consensus_alternative(observed, graph)
    n = length(observed)
    n < 3 && return nothing
    firstpos = Dict{Any, Int}()
    for i in 1:n
        haskey(firstpos, observed[i]) || (firstpos[observed[i]] = i)
    end
    for i in 1:(n - 1)
        a, b = observed[i], observed[i + 1]
        w_ab = _edge_weight_between(graph, a, b)
        w_ab === nothing && continue
        best_x, best_w = nothing, w_ab
        for x in _incident_labels(graph, a)
            x == b && continue
            (i > 1 && x == observed[i - 1]) && continue
            wx = _edge_weight_between(graph, a, x)
            wx === nothing && continue
            if wx > best_w
                best_w, best_x = wx, x
            end
        end
        best_x === nothing && continue   # no strictly-better sibling ⇒ not weak
        # Greedy highest-weight walk from the sibling until we rejoin `observed`
        # at or after position i+1.
        mid = Any[best_x]
        prev, cur = a, best_x
        rejoin = get(firstpos, best_x, 0) >= i + 1 ? firstpos[best_x] : 0
        steps = 0
        while rejoin == 0 && steps < n
            steps += 1
            nxt, nxt_w = nothing, 0.0
            for y in _incident_labels(graph, cur)
                y == prev && continue
                y in mid && continue
                wy = _edge_weight_between(graph, cur, y)
                wy === nothing && continue
                if wy > nxt_w
                    nxt_w, nxt = wy, y
                end
            end
            nxt === nothing && break
            if get(firstpos, nxt, 0) >= i + 1
                rejoin = firstpos[nxt]
                break
            end
            push!(mid, nxt)
            prev, cur = cur, nxt
        end
        if rejoin >= i + 1
            return vcat(observed[1:i], mid, observed[rejoin:end])
        end
    end
    return nothing
end

"""
    accumulate_competing_paths!(accumulator, read, graph, k; graph_mode) -> accumulator

Soft-EM v2 competing-paths E-step for ONE read (td-e70t). Builds the read's
COMPETING candidate paths — the observed read path and (when the observed path
takes a weak branch) a consensus alternative re-routed through the best-supported
sibling (`_consensus_alternative`) — scores each by its normalized-transition
log-probability on the graph (`_path_transition_logscore`, which reads the
soft-weight-aware `compute_edge_weight`), converts the scores to responsibilities
with a stable softmax (`Mycelia.Rhizomorph.path_responsibility`), and accumulates
each candidate's graph-oriented edges weighted by its responsibility.

The responsibility split alone would decay a real but skewed minority allele
toward zero (`W_min' = N*W/(W_maj+W_min)`); variation is preserved by the SUPPORT
FLOOR applied when the accumulator is registered in the M-step
(`register_soft_edge_weights!`), which clamps every `>= SOFT_EM_MIN_SUPPORT` edge
to at least its raw coverage. So this E-step is free to split responsibility while
supported edges are held at raw and only unsupported (error) edges decay. A
balanced variant produces no strictly-better sibling, so no alternative is
generated and the observed path keeps responsibility 1.0. Never throws (guarded),
so the decode is always safe.
"""
function accumulate_competing_paths!(
        accumulator::Mycelia.Rhizomorph.SoftEdgeWeightAccumulator,
        read::FASTX.FASTQ.Record,
        graph,
        k::Int;
        graph_mode::Symbol = :canonical
)
    observed = _read_resolved_labels(read, graph, k; graph_mode = graph_mode)
    isempty(observed) && return accumulator

    alternative = try
        _consensus_alternative(observed, graph)
    catch
        nothing
    end

    if alternative === nothing || alternative == observed
        # No genuine competition: the observed path owns all responsibility.
        Mycelia.Rhizomorph.accumulate_path_probability!(
            accumulator, _graph_oriented_edges(observed, graph), 1.0)
        return accumulator
    end

    candidate_paths = (observed, alternative)
    scores = Float64[_path_transition_logscore(labels, graph) for labels in candidate_paths]
    finite = Float64[s for s in scores if isfinite(s)]
    if isempty(finite)
        # Fall back to the observed path if scoring degenerated.
        Mycelia.Rhizomorph.accumulate_path_probability!(
            accumulator, _graph_oriented_edges(observed, graph), 1.0)
        return accumulator
    end
    for (idx, labels) in enumerate(candidate_paths)
        isfinite(scores[idx]) || continue
        responsibility = Mycelia.Rhizomorph.path_responsibility(scores[idx], finite)
        Mycelia.Rhizomorph.accumulate_path_probability!(
            accumulator, _graph_oriented_edges(labels, graph), responsibility)
    end
    return accumulator
end

"""
Generate alternative qualmer paths through the graph using quality-aware probabilistic sampling.
"""
function generate_alternative_qualmer_paths(
        record::FASTX.FASTQ.Record, graph, k::Int; num_samples::Int = 5,
        graph_mode::Symbol = :canonical)::Vector{Vector{<:Qualmer}}
    paths = Vector{Vector{Qualmer}}()

    try
        # Generate qualmers using existing qualmer-analysis functions
        original_qualmers = collect(qualmers_unambiguous(record, k))
        if isempty(original_qualmers)
            return Vector{Vector{Qualmer}}()
        end

        original_path = Qualmer[]
        for (qmer, pos) in original_qualmers
            resolved = _resolve_qualmer_for_graph(graph, qmer; graph_mode = graph_mode)
            if resolved !== nothing
                push!(original_path, resolved)
            end
        end
        if isempty(original_path)
            return Vector{Vector{Qualmer}}()
        end

        # For each sample, create alternative path using quality-aware sampling
        for sample in 1:num_samples
            alternative_path = copy(original_path)

            # Select positions to modify based on quality confidence
            low_confidence_positions = Int[]
            for (i, qmer) in enumerate(original_path)
                current_kmer = qmer.kmer
                if haskey(graph, current_kmer)
                    vertex_data = graph[current_kmer]
                    # Prefer modifying low-confidence positions
                    if _vertex_joint_probability(vertex_data) < 0.8 ||
                       _vertex_mean_quality(vertex_data) < 30.0
                        push!(low_confidence_positions, i)
                    end
                end
            end

            # If no low confidence positions, randomly select some
            modification_positions = if !isempty(low_confidence_positions)
                Random.shuffle(low_confidence_positions)[1:min(
                    2, length(low_confidence_positions))]
            else
                Random.shuffle(1:length(original_path))[1:min(
                    2, length(original_path))]
            end

            for pos in modification_positions
                if pos <= length(alternative_path)
                    current_qualmer = alternative_path[pos]
                    current_kmer = current_qualmer.kmer

                    if haskey(graph, current_kmer)
                        # Get high-quality neighbor k-mers for potential transitions
                        vertex_data = graph[current_kmer]
                        vertex_coverage = _vertex_coverage(vertex_data)
                        all_vertices = collect(MetaGraphsNext.labels(graph))

                        # Filter neighbors that could be valid transitions
                        valid_neighbors = []
                        neighbor_weights = Float64[]

                        for neighbor_kmer in all_vertices
                            neighbor_data = graph[neighbor_kmer]

                            # Quality-aware weighting combining multiple factors
                            weight = _vertex_joint_probability(neighbor_data) *
                                     (_vertex_mean_quality(neighbor_data) / 60.0) *  # Normalize quality (max ~60)
                                     (_vertex_coverage(neighbor_data) /
                                      max(1, vertex_coverage))  # Relative coverage

                            if weight > 0.01  # Minimum threshold for consideration
                                push!(valid_neighbors, neighbor_kmer)
                                push!(neighbor_weights, weight)
                            end
                        end

                        if !isempty(valid_neighbors)
                            # Sample neighbor based on quality-aware weights
                            neighbor_weights ./= sum(neighbor_weights)  # Normalize
                            selected_idx = StatsBase.sample(
                                1:length(valid_neighbors), StatsBase.Weights(neighbor_weights))

                            # Create new qualmer using the selected k-mer and quality from original
                            selected_kmer = valid_neighbors[selected_idx]
                            original_qualities = current_qualmer.qualities

                            new_qualmer = Qualmer(selected_kmer, original_qualities)
                            alternative_path[pos] = new_qualmer
                        end
                    end
                end
            end

            push!(paths, alternative_path)
        end

    catch e
        # If path generation fails, return empty paths
        @debug "Alternative path generation failed: $e"
        return Vector{Vector{Qualmer}}()
    end

    return paths
end

"""
Convert sequence with quality to qualmer path representation using graph vertices.
"""
function sequence_to_qualmer_path(sequence::String, quality::String, graph, k::Int;
        graph_mode::Symbol = :canonical)::Vector{<:Qualmer}
    path = Qualmer[]

    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels) || length(sequence) < k
        return path
    end

    kmer_type = typeof(first(labels))

    for i in 1:(length(sequence) - k + 1)
        kmer_str = sequence[i:(i + k - 1)]
        qual_segment = quality[i:(i + k - 1)]
        quality_scores = UInt8.(codeunits(qual_segment)) .- UInt8(33)

        qmer = Qualmer(kmer_type(kmer_str), quality_scores)
        resolved = _resolve_qualmer_for_graph(graph, qmer; graph_mode = graph_mode)
        if resolved !== nothing
            push!(path, resolved)
        end
    end

    return path
end

"""
Convert sequence to k-mer path representation (simplified version for backward compatibility).
"""
function sequence_to_kmer_path(sequence::String, k::Int)::Vector{String}
    path = String[]

    for i in 1:(length(sequence) - k + 1)
        push!(path, sequence[i:(i + k - 1)])
    end

    return path
end

function _kmer_to_string(kmer)::String
    return kmer isa AbstractString ? String(kmer) : string(kmer)
end

"""
Convert k-mer path back to sequence string.
"""
function kmer_path_to_sequence(path::AbstractVector, k::Int)::String
    if isempty(path)
        return ""
    end

    first_kmer = _kmer_to_string(path[1])
    if isempty(first_kmer)
        return ""
    end

    buffer = IOBuffer()
    write(buffer, first_kmer)

    for i in 2:length(path)
        kmer_str = _kmer_to_string(path[i])
        if !isempty(kmer_str)
            write(buffer, kmer_str[end])
        end
    end

    return String(take!(buffer))
end

"""
Convert a generic path back to sequence string (alias for k-mer paths).
"""
function path_to_sequence(path::AbstractVector, k::Int)::String
    return kmer_path_to_sequence(path, k)
end

"""
Convert qualmer path back to BioSequence and quality vector.
Returns a named tuple with sequence and quality_scores.
"""
function qualmer_path_to_biosequence(path::Vector{<:Qualmer})
    if isempty(path)
        return (sequence = BioSequences.LongDNA{4}(""), quality_scores = Int8[])
    end

    # Determine sequence type from first qualmer
    first_kmer = path[1].kmer
    sequence_type = if first_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif first_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else  # AAKmer
        BioSequences.LongAA
    end

    # Start with first k-mer
    sequence_symbols = collect(path[1].kmer)
    quality_scores = collect(Int8.(path[1].qualities))

    # Add last symbol of each subsequent k-mer
    for i in 2:length(path)
        push!(sequence_symbols, path[i].kmer[end])
        push!(quality_scores, Int8(path[i].qualities[end]))
    end

    # Construct the BioSequence
    final_sequence = sequence_type(sequence_symbols)

    return (sequence = final_sequence, quality_scores = quality_scores)
end

"""
Convert k-mer path back to BioSequence.
"""
function kmer_path_to_biosequence(path::Vector{<:Kmers.Kmer})
    if isempty(path)
        return BioSequences.LongDNA{4}("")
    end

    # Determine sequence type from first k-mer
    first_kmer = path[1]
    sequence_type = if first_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif first_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else  # AAKmer
        BioSequences.LongAA
    end

    # Start with first k-mer
    sequence_symbols = collect(path[1])

    # Add last symbol of each subsequent k-mer
    for i in 2:length(path)
        push!(sequence_symbols, path[i][end])
    end

    # Construct the BioSequence
    return sequence_type(sequence_symbols)
end
