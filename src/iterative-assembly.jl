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
        beam_width::Union{Int, Nothing} = nothing)
    start_time = time()

    # Accumulates swallowed decode failures (structural / un-k-merizable) across
    # every k and iteration so a systematically-broken corrector that reports
    # "0 improvements" is visible in the result metadata rather than passing as
    # "nothing to fix". See CorrectorDiagnostics.
    corrector_diagnostics = CorrectorDiagnostics()

    # -- Soft-EM registry hygiene (td-e70t, v1 record-only) --------------------
    # Soft-EM v1 is a RECORD-ONLY scaffold: the E-step accumulates per-edge path
    # responsibilities (for the future v2 competing-paths M-step) but the pipeline
    # NEVER calls `register_soft_edge_weights!`, so `compute_edge_weight` always
    # returns raw coverage counts and the identity-keyed soft-weight registry MUST
    # stay empty for this whole run. Defensively clear it at entry so no prior
    # crash mid-registration (a direct unit-test primitive call, or a future v2
    # path) can leak soft weights into this correction, then assert the invariant
    # (belt-and-suspenders). The registry is reserved for v2.
    Mycelia.Rhizomorph.clear_soft_edge_weights!()
    @assert isempty(Mycelia.Rhizomorph._SOFT_EDGE_WEIGHT_REGISTRY) (
        "soft-EM registry must be empty at corrector entry (v1 is record-only; " *
        "the pipeline never registers soft weights)")

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

    # Hard-window skip telemetry (td-nn6l): the fraction of reads the hard-window
    # gate passed through WITHOUT a decode, recorded PER PASS (across every k and
    # iteration) so the run metadata can surface min/mean/max, not just the last
    # value (FIX 6). Empty on the :exhaustive tier (no gate ⇒ no skips).
    skip_fractions = Float64[]

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

            # --- Soft-EM E-step accumulation (td-e70t, v1 RECORD-ONLY) --------
            # v1 is a DORMANT scaffold: a fresh accumulator tallies THIS pass's
            # decoded-path responsibilities (the E-step), but the pipeline does NOT
            # consume it — `register_soft_edge_weights!` is deliberately NOT called,
            # so `compute_edge_weight` always returns raw coverage counts and the
            # graph is never soft-reweighted. This keeps the accumulation machinery
            # live and exercised (recording path probabilities for the future v2
            # competing-paths M-step) while removing v1's uncontrolled perturbation:
            # a single-path decode makes every responsibility 1.0 and hard-skipped
            # reads never accumulate, so consuming the accumulator would mix
            # DECODED-SUBSET coverage with raw counts on unvisited edges — an
            # uncontrolled reweighting, not principled soft-EM. Only allocated under
            # soft-EM so :exhaustive threads `nothing`.
            current_soft_weights = soft_em ?
                                   Mycelia.Rhizomorph.SoftEdgeWeightAccumulator() : nothing

            # --- Hard-window gating (td-nn6l) --------------------------------
            # Restrict decoding to reads that touch a "hard" vertex (bubble
            # entry/exit/interior, high-out-degree/repeat-like, or weak/non-solid
            # k-mer). Easy reads pass through untouched, cutting decode volume
            # ~85-95% on clean data. Off unless `hard_window=true` (:scalable).
            hard_vertices = (hard_window && graph_mode == :canonical) ?
                            _hard_vertex_set(graph, k) : nothing
            if hard_window && graph_mode != :canonical
                @warn "hard_window gating is only supported for graph_mode=:canonical; " *
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
            updated_reads,
            improvements_made,
            pass_skip_fraction = improve_read_set_likelihood(
                current_reads, graph, k,
                verbose = verbose,
                batch_size = batch_size,
                enable_parallel = enable_parallel,
                graph_mode = graph_mode,
                skip_solid = skip_solid,
                beam_width = beam_width,
                soft_weights = current_soft_weights,
                hard_vertices = hard_vertices,
                diagnostics = corrector_diagnostics
            )
            push!(skip_fractions, pass_skip_fraction)
            # Soft-EM v1 is record-only: `current_soft_weights` was populated by the
            # E-step but is deliberately NOT registered/consumed, so there is nothing
            # to clear from the registry (it stayed empty) and no cross-iteration
            # carry-forward — the accumulator is dropped once recorded. The v2
            # competing-paths M-step will consume it.

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
        hard_window = hard_window, soft_em = soft_em)
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

"""
    _read_is_all_solid(read, k, solid_kmers) -> Bool

True iff every canonical k-mer of `read` is in `solid_kmers` (⇒ no weak region ⇒
the read needs no correction and can skip the decode). A read with no k-mers
(shorter than k, or fully ambiguous) returns `false` so it is never skipped on the
basis of absent evidence.
"""
function _read_is_all_solid(read::FASTX.FASTQ.Record, k::Int, solid_kmers::AbstractSet)
    isempty(solid_kmers) && return false
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    saw_kmer = false
    for (kmer, _) in Kmers.UnambiguousDNAMers{k}(sequence)
        saw_kmer = true
        (BioSequences.canonical(kmer) in solid_kmers) || return false
    end
    return saw_kmer
end

# ------------------------------------------------------------------------------
# Hard-window gating (td-nn6l): decode only reads that touch a "hard" region
# (bubble / repeat / weak k-mer) and pass every other read through untouched.
# This is the primary decode-volume reduction on top of skip-solid — on clean
# data the vast majority of reads touch no hard vertex and are skipped.
# ------------------------------------------------------------------------------

"""
    should_decode_read(read, k, hard_vertices) -> Bool

True iff any canonical k-mer of `read` is a member of `hard_vertices` (the
bubble / repeat-like / weak-k-mer set). Reads that touch no hard vertex have no
region worth decoding and are passed through untouched. Mirrors
`_read_is_all_solid`'s canonical k-mer iteration so the two gates compose on the
same `:canonical` graph. An empty `hard_vertices` means "no gate" ⇒ decode.
"""
function should_decode_read(read::FASTX.FASTQ.Record, k::Int, hard_vertices::AbstractSet)
    isempty(hard_vertices) && return true
    sequence = FASTX.sequence(BioSequences.LongDNA{4}, read)
    for (kmer, _) in Kmers.UnambiguousDNAMers{k}(sequence)
        (BioSequences.canonical(kmer) in hard_vertices) && return true
    end
    return false
end

"""
    _hard_vertex_set(graph, k) -> Set

Build the set of "hard" vertices for hard-window gating (td-nn6l): the union of
(1) bubble entry/exit/interior vertices (`detect_bubbles_next`),
(2) high-out-degree (out-degree > 1, repeat-like) vertices, and
(3) weak / non-solid k-mers (the complement of `_solid_kmer_set`).
When k-mer classification cannot run (empty solid set) every vertex is treated as
hard — the conservative "decode everything" fallback that never skips on absent
evidence. Intended for `:canonical` graphs (matches the skip-solid gate).
"""
function _hard_vertex_set(graph, k::Int)
    labels = collect(MetaGraphsNext.labels(graph))
    T = eltype(labels)
    hard = Set{T}()
    isempty(labels) && return hard

    # (1) Bubble vertices (entry, exit, both alternative interiors).
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
    outdeg = Dict{T, Int}()
    for (src, _dst) in MetaGraphsNext.edge_labels(graph)
        outdeg[src] = get(outdeg, src, 0) + 1
    end
    for (v, d) in outdeg
        d > 1 && push!(hard, v)
    end

    # (3) Weak (non-solid) k-mers = complement of the solid set. Empty solid set
    # (classification could not run) ⇒ treat all vertices as hard (decode all).
    solid = _solid_kmer_set(graph)
    if isempty(solid)
        union!(hard, labels)
    else
        for v in labels
            (v in solid) || push!(hard, v)
        end
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

`soft_weights` (td-e70t): when non-`nothing`, each decoded read's ML path
accumulates edge responsibilities into it (the soft-EM E-step). Accumulation is
NOT thread-safe, so supplying it forces sequential processing for this pass.

Returns `(updated_reads, improvements_made, skip_fraction)` where `skip_fraction`
is the fraction of reads passed through WITHOUT a decode (solid + hard-window
skips). Callers destructuring only the first two values are unaffected.
"""
function improve_read_set_likelihood(reads::Vector{<:FASTX.FASTQ.Record}, graph, k::Int;
        verbose::Bool = false,
        batch_size::Int = 10000,
        enable_parallel::Bool = false,
        graph_mode::Symbol = :canonical,
        skip_solid::Bool = false,
        beam_width::Union{Int, Nothing} = nothing,
        soft_weights::Union{Nothing, Mycelia.Rhizomorph.SoftEdgeWeightAccumulator} = nothing,
        hard_vertices::Union{Nothing, AbstractSet} = nothing,
        diagnostics::Union{Nothing, CorrectorDiagnostics} = nothing)::Tuple{
        Vector{FASTX.FASTQ.Record}, Int, Float64}
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
    # The skip helpers canonicalize read k-mers, so they are only valid on a
    # :canonical graph; on any other mode disable the skip (correct all) rather
    # than mis-match strands and wrongly skip (review I1/I2 graph_mode).
    if skip_solid && graph_mode != :canonical
        @warn "skip_solid is only supported for graph_mode=:canonical; disabling skip (correcting all reads)." graph_mode
    end
    solid_kmers = (skip_solid && graph_mode == :canonical) ? _solid_kmer_set(graph) :
                  nothing
    skipped_reads = 0

    # A read is SKIPPED (passed through uncorrected, no decode) when it is either
    # all-solid (Stage 0) OR — under hard-window gating — it touches no hard
    # vertex. Both are volume-reduction skips and both count toward the reported
    # skip fraction. Returns true ⇒ skip.
    _skip_this_read = read -> begin
        (solid_kmers !== nothing && _read_is_all_solid(read, k, solid_kmers)) && return true
        (hard_vertices !== nothing && !should_decode_read(read, k, hard_vertices)) &&
            return true
        return false
    end

    # Soft-EM accumulation into `soft_weights` is not thread-safe (shared Dict),
    # so a soft-EM pass runs sequentially. Otherwise honor the caller's request.
    use_parallel = enable_parallel && Threads.nthreads() > 1 && soft_weights === nothing
    if enable_parallel && soft_weights !== nothing && Threads.nthreads() > 1
        @warn "soft-EM edge accumulation is sequential (race-free); ignoring enable_parallel for this pass." maxlog = 1
    end

    if verbose
        println("  Processing $total_reads reads in batches of $batch_size")
    end

    # Process in batches for memory efficiency
    for batch_start in 1:batch_size:total_reads
        batch_end = min(batch_start + batch_size - 1, total_reads)
        batch_reads = reads[batch_start:batch_end]

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
                if _skip_this_read(read)
                    batch_results[i] = (read, false)   # skip the decode
                    skip_flags[i] = true
                else
                    improved_read,
                    was_improved = improve_read_likelihood(
                        read, graph, k; graph_mode = graph_mode,
                        beam_width = beam_width, diagnostics = diag)
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
                if _skip_this_read(read)
                    updated_reads[batch_start + i - 1] = read   # skip the decode
                    skipped_reads += 1
                    continue
                end
                improved_read,
                was_improved = improve_read_likelihood(
                    read, graph, k; graph_mode = graph_mode,
                    beam_width = beam_width, soft_weights = soft_weights,
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
        if skip_solid || hard_vertices !== nothing
            println("  Skipped (no decode): $skipped_reads/$total_reads ($(round(skip_fraction * 100, digits=1))%)")
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

    return updated_reads, improvements_made, skip_fraction
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
        hard_window::Bool = false,
        soft_em::Bool = false)
    if verbose
        println("Finalizing iterative assembly results...")
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
        # Soft-EM v1 is RECORD-ONLY (E-step accumulates, M-step does NOT consume),
        # so report a scaffold marker rather than a bare `true` that would imply
        # active soft reweighting (FIX 1/5). `false` on the exhaustive tier.
        :soft_em => (soft_em ? "scaffold-v1-record-only" : false),
        :last_skip_fraction => last_skip_fraction,
        :skip_fraction_per_pass => skip_fractions,
        :skip_fraction_min => skip_min,
        :skip_fraction_mean => skip_mean,
        :skip_fraction_max => skip_max
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
        :metadata => metadata
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

        # Trustworthy Viterbi (td-ak6w), reconciled with the td-63qy OOM crash:
        # the decoder is EXACT where tractable and BOUNDED (with a log line) above
        # a size threshold. `beam_width === nothing` (the default) resolves the
        # frontier bound PER READ via `_auto_beam_width(length(observations))`:
        # `typemax(Int)` (never pruned → global maximum-likelihood path) for small
        # reads, the proven-tractable finite bound for large reads (the unbounded
        # frontier grows ~unboundedly with read-length depth — 21B allocations on
        # a 48 kb phage). An explicit `beam_width::Int` (including `typemax(Int)`)
        # is an OPT-IN override that is honored verbatim — a caller/test can force
        # exact decoding on every read. This trades the ML guarantee for
        # tractability ONLY above the threshold, and only under the auto-default;
        # it never silently changes correctness of an explicit request.
        effective_beam_width = beam_width === nothing ?
                               _auto_beam_width(length(observations)) : beam_width
        config = Mycelia.ViterbiCorrectionConfig(
            alphabet = alphabet,
            strand_mode = graph_mode,
            max_steps = length(observations) - 1,
            beam_width = effective_beam_width
        )
        correction = Mycelia.correct_observations(graph, [observations]; config = config)
        corrected_path = only(correction.corrected_observations)
        if corrected_path === nothing || isempty(corrected_path)
            return nothing
        end

        # Soft-EM E-step (td-e70t): accumulate this read's decoded ML path edges
        # onto the shared accumulator. With a single argmax path per read the
        # responsibility is 1.0 (soft weight == hard count); the softening emerges
        # once competing paths per read are wired (v2). Additive + guarded so it
        # cannot perturb the :exhaustive path (soft_weights === nothing there).
        if soft_weights !== nothing
            accumulate_soft_em_edge_weights!(soft_weights, correction.paths)
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
# Soft-EM M-step accumulation hook (SCAFFOLD, td-e70t)
# ============================================================================
#
# See the design note over `SoftEdgeWeightAccumulator` in
# `src/rhizomorph/core/evidence-functions.jl`. This is the correction-side hook
# that turns a read's competing decoded paths into probability-weighted edge
# evidence: each candidate path's RESPONSIBILITY (posterior over the read's own
# candidate set) is accumulated onto the edges it traverses. Summed over reads,
# the accumulator holds soft edge weights in which error edges — traversed only
# by rare, low-probability paths — accrue little weight and decay below the
# `generate_alternative_qualmer_paths` `weight > 0.01` gate WITHOUT tip-clipping.
#
# v1 status — E-step WIRED (record-only), M-step NOT wired. Under `soft_em=true`
# this hook IS called per decode from `try_viterbi_path_improvement`, so the
# accumulator records decoded-path responsibilities each pass. But the pipeline
# deliberately does NOT consume them: `register_soft_edge_weights!` is never
# called from `mycelia_iterative_assemble`, so `compute_edge_weight` always
# returns raw coverage counts and the graph is never soft-reweighted. Emergent
# coalescence needs the M-step to consume these weights (and, before that, the
# v2 competing-paths E-step so responsibilities stop being a degenerate 1.0),
# which is the tracked td-e70t follow-on. Keeping the E-step live but the M-step
# dormant avoids v1's uncontrolled reweighting (single-path 1.0 responsibilities
# + hard-skipped reads never accumulating would mix decoded-subset coverage with
# raw counts) while exercising the recording machinery for v2.

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

Soft-EM M-step seam (td-e70t scaffold). Given a single read's competing decoded
paths (`ViterbiDecodingResult`s, each with a `.score` log-likelihood and a
`.path`), compute each path's responsibility by softmax-normalizing its score
against the candidate set (`Mycelia.Rhizomorph.path_responsibility`) and
accumulate that responsibility onto the edges it traverses
(`Mycelia.Rhizomorph.accumulate_path_probability!`).

With a single argmax path per read the responsibility is `1.0` (soft weight ==
hard count); the softening emerges once several candidate paths per read compete
(e.g. from `generate_alternative_qualmer_paths`), which is the regime the
follow-on wires into the iteration loop.
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
