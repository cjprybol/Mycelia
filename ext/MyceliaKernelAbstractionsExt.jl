# Mycelia package extension: backend-agnostic (CPU + CUDA/GPU) frontier kernel for
# the batched array-frontier Viterbi corrector — Phase B of bead td-qoo3, built on
# the merged CPU PoC (`src/rhizomorph/algorithms/batched-viterbi-poc.jl`).
#
# This module is loaded lazily by Julia's package-extension mechanism when (and
# only when) the user does `import KernelAbstractions`. It adds real methods to
# the `Mycelia.kernel_batched_*` stubs declared in
# `src/rhizomorph/algorithms/batched-viterbi-kernel.jl`.
#
# Byte-identity contract (the whole point): the kernel produces `next_scores` /
# `next_predecessor` that are BIT-identical to the CPU PoC's inner relaxation loop
# (`_batched_decode_bin`), because (a) one GPU thread handles one read and runs
# the SAME sequential scan over states/out-edges in the SAME order, (b) the
# arithmetic keeps the SAME `(state_score + log_transition) + emission`
# associativity, (c) the relax test is the SAME strict `>` (first-seen wins on a
# tie, so predecessors match), and (d) every value the kernel reads — the CSR
# transition arrays and the `emission[read, state, depth]` table — is precomputed
# on the host with the SAME `_call_viterbi_state_emission_logp` /
# `log(edge_w/total_out)` the scalar decoder uses. Start-init, per-depth
# predecessor materialization + best-state selection, and path reconstruction stay
# host-side and reuse the scalar tie-break helpers verbatim, exactly as the PoC.

module MyceliaKernelAbstractionsExt

import Mycelia
import KernelAbstractions
import MetaGraphsNext
import Random

# KernelAbstractions is a kernel DSL: the `@kernel` macro rewrites the kernel
# body and only recognizes the index/const DSL macros when they are in scope
# unqualified, so these two DSL macro names are brought in by name (the rest of
# the API — including `@kernel` itself, invoked below as `KA.@kernel`, plus
# allocate/synchronize/CPU/etc. — stays fully qualified as `KA.`).
using KernelAbstractions: @index, @Const

const KA = KernelAbstractions

# ---------------------------------------------------------------------------
# The wavefront kernel: advance the dense `n_reads x n_states` frontier ONE
# depth-step for every read in a length-bin in parallel (one thread per read —
# the ADR's zero-write-contention, fully-data-parallel read axis). Runs unchanged
# on the CPU backend (CI, no GPU) and on a CUDA backend (GPU hosts).
#
# `current_scores`, the CSR arrays, and `emission` are read-only (@Const). Each
# thread writes ONLY its own row of `next_scores` / `next_pred`, so there is zero
# cross-thread write contention and no atomics are needed — the parallelism is
# purely across reads, exactly as the CPU PoC's depth-outer/read-inner loop.
# ---------------------------------------------------------------------------
KA.@kernel function _frontier_step_kernel!(
        next_scores, next_pred,
        @Const(current_scores), @Const(edge_offsets),
        @Const(edge_dst), @Const(edge_logprob), @Const(emission),
        depth, n_states)
    r = @index(Global)
    # `ndrange` may be padded up to a multiple of the workgroup size; mask the
    # out-of-range lanes (KernelAbstractions does not do this automatically).
    if r <= size(current_scores, 1)
        # Reset this read's next-frontier row (mirrors PoC lines 309-315).
        for sid in 1:n_states
            @inbounds next_scores[r, sid] = -Inf
            @inbounds next_pred[r, sid] = Int32(0)
        end
        # Relax every finite state against the shared read-only out-edge SoA,
        # in state order then edge order — identical to the PoC inner loop.
        for sid in 1:n_states
            @inbounds state_score = current_scores[r, sid]
            if isfinite(state_score)
                @inbounds lo = edge_offsets[sid]
                @inbounds hi = edge_offsets[sid + 1]
                for e in (lo + 1):hi
                    @inbounds dst = edge_dst[e]
                    @inbounds emission_score = emission[r, dst, depth]
                    # Same associativity as the PoC: (state + log_transition) + emission.
                    @inbounds candidate = state_score + edge_logprob[e] + emission_score
                    if isfinite(candidate)
                        @inbounds if candidate > next_scores[r, dst]
                            next_scores[r, dst] = candidate
                            next_pred[r, dst] = Int32(sid)
                        end
                    end
                end
            end
        end
    end
end

# Allocate a backend array and copy a host array into it (host->device when the
# backend is a GPU; a plain copy on the CPU backend). `Base.copyto!` is overloaded
# by CUDA.jl for host<->device transfers, so this stays backend-agnostic.
function _to_device(backend, host_array::AbstractArray{T}) where {T}
    device_array = KA.allocate(backend, T, size(host_array))
    Base.copyto!(device_array, host_array)
    return device_array
end

# Flatten a BatchedDecodeGraph's per-state out-edge lists into a read-only CSR
# structure-of-arrays. Edge order within a state is preserved from the PoC's
# `out_edges[sid]`, so the kernel's tie-break order matches the PoC exactly.
function _build_csr(decode_graph::Mycelia.BatchedDecodeGraph)
    n_states = length(decode_graph.states)
    edge_offsets = Vector{Int32}(undef, n_states + 1)
    edge_offsets[1] = 0
    for sid in 1:n_states
        edge_offsets[sid + 1] = edge_offsets[sid] + length(decode_graph.out_edges[sid])
    end
    n_edges = Int(edge_offsets[end])
    edge_dst = Vector{Int32}(undef, n_edges)
    edge_logprob = Vector{Float64}(undef, n_edges)
    pos = 1
    for sid in 1:n_states
        for edge in decode_graph.out_edges[sid]
            edge_dst[pos] = Int32(edge.next_state_id)
            edge_logprob[pos] = edge.log_transition_prob
            pos += 1
        end
    end
    return edge_offsets, edge_dst, edge_logprob
end

# Precompute the read x state x depth emission table on the host using the SAME
# emission function the scalar/PoC decoder calls, so every value the kernel reads
# is byte-identical to the scalar per-step arithmetic. `emission[r, sid, depth]`
# is the emission log-prob of read `r`'s observed unit at position `depth+1`
# against the destination vertex of state `sid`.
function _build_emission_table(
        quality_graph, config, alphabet, strand_mode,
        states, bin_observations, observation_length)
    n = length(bin_observations)
    n_states = length(states)
    n_depths = max(observation_length - 1, 0)
    emission = fill(-Inf, n, n_states, n_depths)
    for depth in 1:n_depths
        for r in 1:n
            observed_unit = bin_observations[r][depth + 1]
            for sid in 1:n_states
                vertex = states[sid][1]
                emission[r, sid, depth] = Mycelia._call_viterbi_state_emission_logp(
                    quality_graph, config, observed_unit, vertex, alphabet, strand_mode)
            end
        end
    end
    return emission
end

# Kernel-driven length-bin decode. Faithful copy of the CPU PoC's
# `_batched_decode_bin` host scaffolding (start-init, per-depth commit, finalize),
# with ONLY the inner frontier relaxation replaced by a `@kernel` launch. Every
# other line reuses the exact scalar/PoC helpers, so the result is byte-identical.
function _kernel_decode_bin(
        backend,
        graph, quality_graph, config, alphabet, strand_mode,
        transition_scoring, labels, decode_graph, bin_observations)
    states = decode_graph.states
    state_id = decode_graph.state_id
    n = length(bin_observations)
    n_states = length(states)
    observation_length = length(bin_observations[1])
    L = eltype(labels)
    StateT = Tuple{L, Mycelia.Rhizomorph.StrandOrientation}

    current_scores = fill(-Inf, n, n_states)
    next_scores = fill(-Inf, n, n_states)
    next_predecessor = zeros(Int, n, n_states)

    finished = falses(n)
    best_state = Vector{Union{Nothing, StateT}}(nothing, n)
    best_score = fill(-Inf, n)
    best_depth = zeros(Int, n)
    targets = Vector{Any}(undef, n)
    predecessors_by_depth = [Vector{Dict{StateT, StateT}}() for _ in 1:n]
    diagnostics = Vector{Dict{Symbol, Any}}(undef, n)
    results = Vector{Mycelia.Rhizomorph.ViterbiDecodingResult}(undef, n)

    # ---- depth 0: per-read start-emission init (identical to PoC) ----
    for r in 1:n
        observation = bin_observations[r]
        start_observed = first(observation)
        target_vertex = Mycelia._emission_target_vertex(
            graph, observation, config, alphabet, strand_mode)
        targets[r] = target_vertex
        diag = Mycelia._batched_new_diagnostics(
            alphabet, strand_mode, observation_length, target_vertex,
            config, quality_graph, transition_scoring)
        diagnostics[r] = diag

        start_candidates = Mycelia._viterbi_start_candidates(
            labels, Mycelia._viterbi_label_unit(start_observed), alphabet, strand_mode)
        for vertex in start_candidates
            for strand in Mycelia._viterbi_start_strands(
                    graph, vertex, strand_mode, config.start_strand)
                sid = state_id[(convert(L, vertex), strand)]
                score = Mycelia._call_viterbi_state_emission_logp(
                    quality_graph, config, start_observed, vertex, alphabet, strand_mode)
                if isfinite(score) && score > current_scores[r, sid]
                    current_scores[r, sid] = score
                end
            end
        end

        reached = Mycelia._batched_reached_scores(current_scores, r, states)
        if isempty(reached)
            diag[:reason] = :no_finite_start_emission
            results[r] = Mycelia.Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diag)
            finished[r] = true
            continue
        end
        diag[:retained_states] = length(reached)
        diag[:cumulative_retained_states] = length(reached)
        diag[:max_retained_states] = length(reached)

        bstate, bscore = Mycelia._best_correction_state(reached)
        best_state[r] = bstate
        best_score[r] = bscore
        best_depth[r] = 0
        if target_vertex !== nothing
            if observation_length == 1
                tstate, tscore = Mycelia._best_correction_target_state(reached, target_vertex)
                if tstate !== nothing
                    best_state[r] = tstate
                    best_score[r] = tscore
                    diag[:reached_target] = true
                else
                    best_score[r] = -Inf
                end
            else
                best_score[r] = -Inf
            end
        end
    end

    # ---- device-resident, read-only shared graph + emission (uploaded once) ----
    edge_offsets_h, edge_dst_h, edge_logprob_h = _build_csr(decode_graph)
    emission_h = _build_emission_table(
        quality_graph, config, alphabet, strand_mode,
        states, bin_observations, observation_length)

    if observation_length > 1
        d_edge_offsets = _to_device(backend, edge_offsets_h)
        d_edge_dst = _to_device(backend, edge_dst_h)
        d_edge_logprob = _to_device(backend, edge_logprob_h)
        d_emission = _to_device(backend, emission_h)
        d_current = _to_device(backend, current_scores)
        d_next = KA.allocate(backend, Float64, (n, n_states))
        d_next_pred = KA.allocate(backend, Int32, (n, n_states))
        kernel = _frontier_step_kernel!(backend)

        # ---- depth-outer, read-inner relaxation ON THE KERNEL (ADR item 2) ----
        for depth in 1:(observation_length - 1)
            # Keep the device current frontier in sync with the host commit.
            Base.copyto!(d_current, current_scores)
            kernel(
                d_next, d_next_pred, d_current,
                d_edge_offsets, d_edge_dst, d_edge_logprob, d_emission,
                depth, n_states; ndrange = n)
            KA.synchronize(backend)
            Base.copyto!(next_scores, d_next)
            Base.copyto!(next_predecessor, d_next_pred)

            # Per-read commit: materialize predecessors, swap frontier, update
            # best. Byte-for-byte the PoC commit (lines 348-392).
            for r in 1:n
                finished[r] && continue
                reached_scores = Dict{StateT, Float64}()
                reached_predecessors = Dict{StateT, StateT}()
                @inbounds for sid in 1:n_states
                    score = next_scores[r, sid]
                    isfinite(score) || continue
                    state = states[sid]
                    reached_scores[state] = score
                    reached_predecessors[state] = states[next_predecessor[r, sid]]
                end

                if isempty(reached_scores)
                    finished[r] = true
                    continue
                end

                push!(predecessors_by_depth[r], reached_predecessors)
                @inbounds for sid in 1:n_states
                    current_scores[r, sid] = next_scores[r, sid]
                end

                diag = diagnostics[r]
                retained = length(reached_scores)
                diag[:retained_states] = retained
                diag[:cumulative_retained_states] += retained
                diag[:max_retained_states] = max(diag[:max_retained_states], retained)
                diag[:completed_steps] = depth

                target_vertex = targets[r]
                if target_vertex === nothing
                    bstate, bscore = Mycelia._best_correction_state(reached_scores)
                    best_state[r] = bstate
                    best_score[r] = bscore
                    best_depth[r] = depth
                else
                    tstate, tscore = Mycelia._best_correction_target_state(
                        reached_scores, target_vertex)
                    if tstate !== nothing
                        best_state[r] = tstate
                        best_score[r] = tscore
                        best_depth[r] = depth
                        diag[:reached_target] = true
                    end
                end
            end
        end
    end

    # ---- finalize per read (identical to PoC lines 396-410) ----
    for r in 1:n
        isassigned(results, r) && continue
        diag = diagnostics[r]
        target_vertex = targets[r]
        if target_vertex !== nothing && !isfinite(best_score[r])
            diag[:reason] = :target_unreachable
            results[r] = Mycelia.Rhizomorph.ViterbiDecodingResult(nothing, -Inf, diag)
            continue
        end
        end_state = something(best_state[r])
        path = Mycelia._reconstruct_correction_path(
            graph, end_state, best_depth[r], predecessors_by_depth[r])
        diag[:path_length] = length(path.steps)
        results[r] = Mycelia.Rhizomorph.ViterbiDecodingResult(path, best_score[r], diag)
    end

    return results
end

# ---------------------------------------------------------------------------
# Public method: kernel-driven batched decode. Mirrors
# `batched_correct_observations`'s setup/guards exactly, then dispatches each
# length-bin to `_kernel_decode_bin` on the selected backend.
# ---------------------------------------------------------------------------
function Mycelia.kernel_batched_correct_observations(
        graph::MetaGraphsNext.MetaGraph,
        observations = Mycelia._default_correction_observations(graph);
        config::Mycelia.ViterbiCorrectionConfig = Mycelia.ViterbiCorrectionConfig(),
        backend = KA.CPU())
    if config.beam_width != typemax(Int)
        throw(ArgumentError(
            "kernel_batched_correct_observations is an exact-decode kernel and " *
            "requires beam_width == typemax(Int) (got $(config.beam_width)); beam " *
            "pruning is a follow-on."))
    end
    if observations === nothing || isempty(observations)
        throw(ArgumentError("no observations to correct"))
    end
    for observation in observations
        Mycelia._uses_emission_scored_observation(observation) || throw(ArgumentError(
            "kernel PoC requires non-empty emission-scored (AbstractVector) " *
            "observations; got $(typeof(observation))"))
    end

    alphabet = Mycelia._resolve_viterbi_alphabet(graph, observations, config.alphabet)
    strand_mode = Mycelia._resolve_viterbi_strand_mode(graph, alphabet, config.strand_mode)
    if strand_mode != :singlestrand
        throw(ArgumentError(
            "kernel PoC supports single-strand decoding only (resolved " *
            "strand_mode=$strand_mode); multi-strand batching is deferred."))
    end

    transition_edge_weight = Mycelia._viterbi_transition_edge_weight(graph, config.edge_weight)
    transition_scoring = Mycelia._viterbi_transition_scoring(graph, transition_edge_weight)
    weighted = if Mycelia._correction_edge_data_type(graph) <: Mycelia.Rhizomorph.StrandWeightedEdgeData
        graph
    else
        Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(
            graph; edge_weight = transition_edge_weight)
    end

    labels = collect(MetaGraphsNext.labels(weighted))
    if isempty(labels)
        throw(ArgumentError("cannot correct observations against an empty graph"))
    end
    label_type = eltype(labels)

    decode_graph = Mycelia.build_batched_decode_graph(weighted, label_type)

    bins = Dict{Int, Vector{Int}}()
    for (index, observation) in enumerate(observations)
        push!(get!(bins, length(observation), Int[]), index)
    end

    paths = Vector{Mycelia.Rhizomorph.ViterbiDecodingResult}(undef, length(observations))
    for (_, indices) in bins
        bin_observations = [observations[index] for index in indices]
        bin_results = _kernel_decode_bin(
            backend, weighted, graph, config, alphabet, strand_mode,
            transition_scoring, labels, decode_graph, bin_observations)
        for (position, index) in enumerate(indices)
            paths[index] = bin_results[position]
        end
    end

    corrected = Any[Mycelia._decoded_path_labels(path_result) for path_result in paths]
    diagnostics = Dict{Symbol, Any}(
        :interface => :metagraphs_next_kernel,
        :vertex_data_type => Mycelia._correction_vertex_data_type(graph),
        :algorithm => :rhizomorph_viterbi_decode_next_kernel,
        :observation_count => length(observations),
        :emission_callback => nameof(config.emission_logp),
        :alphabet => alphabet,
        :emission_model => Mycelia._viterbi_graph_has_quality(graph) ?
                           :quality_aware : :alphabet_parameterized,
        :transition_model => transition_scoring,
        :strand_mode => strand_mode,
        :reverse_complement_support => Mycelia._viterbi_supports_reverse_complement(alphabet),
        :length_bins => sort(collect(keys(bins))),
        :batch_size => length(observations),
        :state_space_size => length(decode_graph.states),
        :backend => string(nameof(typeof(backend))))
    return Mycelia.ViterbiCorrectionResult(corrected, paths, diagnostics)
end

# ---------------------------------------------------------------------------
# Public method: equivalence oracle. Asserts the kernel decode is byte-identical
# to BOTH the scalar `correct_observations` and the CPU PoC
# `batched_correct_observations` (the merged CPU reference).
# ---------------------------------------------------------------------------
function Mycelia.kernel_batched_equivalence_oracle(
        graph::MetaGraphsNext.MetaGraph,
        observations = Mycelia._default_correction_observations(graph);
        config::Mycelia.ViterbiCorrectionConfig = Mycelia.ViterbiCorrectionConfig(),
        backend = KA.CPU())
    scalar = Mycelia.correct_observations(graph, observations; config = config)
    poc = Mycelia.batched_correct_observations(graph, observations; config = config)
    kernel = Mycelia.kernel_batched_correct_observations(
        graph, observations; config = config, backend = backend)

    function _paths_match(reference, candidate, index)
        rp = reference.paths[index]
        cp = candidate.paths[index]
        labels_match = Mycelia._decoded_path_labels(rp) == Mycelia._decoded_path_labels(cp)
        score_match = rp.score === cp.score
        target_match = get(rp.diagnostics, :reached_target, nothing) ==
                       get(cp.diagnostics, :reached_target, nothing)
        steps_match = if rp.path === nothing || cp.path === nothing
            rp.path === cp.path
        else
            rs = something(rp.path).steps
            cs = something(cp.path).steps
            length(rs) == length(cs) && all(
                rs[i].vertex_label == cs[i].vertex_label && rs[i].strand == cs[i].strand
                for i in 1:length(rs))
        end
        return labels_match && score_match && target_match && steps_match
    end

    n = length(scalar.paths)
    per_read_vs_scalar = Bool[_paths_match(scalar, kernel, i) for i in 1:n]
    per_read_vs_poc = Bool[_paths_match(poc, kernel, i) for i in 1:n]
    per_read = per_read_vs_scalar .& per_read_vs_poc

    return (;
        passes = all(per_read),
        per_read = per_read,
        per_read_vs_scalar = per_read_vs_scalar,
        per_read_vs_poc = per_read_vs_poc,
        scalar = scalar, poc = poc, kernel = kernel,
        backend = string(nameof(typeof(backend))))
end

# ---------------------------------------------------------------------------
# Benchmark harness. Isolates the frontier-kernel throughput (the thing Phase B
# accelerates) with a synthetic, device-resident depth sweep — no per-depth host
# round-trips — so the GPU number reflects the kernel, not host sync. Always runs
# the CPU backend; if a functional CUDA backend is available (passed explicitly
# or auto-detected from an already-imported CUDA) it also runs the GPU and reports
# the speedup. Degrades gracefully to CPU-only.
# ---------------------------------------------------------------------------

# Device-resident sweep of `n_depths` frontier steps (no host commit). Returns the
# final `current` frontier copied back to the host so results can be cross-checked
# across backends.
function _frontier_sweep!(
        backend, current_h, edge_offsets_h, edge_dst_h, edge_logprob_h, emission_h)
    n, n_states = size(current_h)
    n_depths = size(emission_h, 3)
    d_current = _to_device(backend, current_h)
    d_next = KA.allocate(backend, Float64, (n, n_states))
    d_next_pred = KA.allocate(backend, Int32, (n, n_states))
    d_edge_offsets = _to_device(backend, edge_offsets_h)
    d_edge_dst = _to_device(backend, edge_dst_h)
    d_edge_logprob = _to_device(backend, edge_logprob_h)
    d_emission = _to_device(backend, emission_h)
    kernel = _frontier_step_kernel!(backend)
    for depth in 1:n_depths
        kernel(d_next, d_next_pred, d_current,
            d_edge_offsets, d_edge_dst, d_edge_logprob, d_emission,
            depth, n_states; ndrange = n)
        KA.synchronize(backend)
        d_current, d_next = d_next, d_current
    end
    out = Array{Float64}(undef, n, n_states)
    Base.copyto!(out, d_current)
    return out
end

# Best-effort detection of an already-loaded, functional CUDA backend WITHOUT
# adding CUDA as a dependency. Returns a `CUDABackend()` or `nothing`.
function _detect_cuda_backend()
    cuda_id = Base.PkgId(
        Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA")
    if !Base.root_module_exists(cuda_id)
        return nothing  # user has not imported CUDA; do not force-load a heavy dep
    end
    cuda = Base.root_module(cuda_id)
    try
        if cuda.functional()
            return cuda.CUDABackend()
        end
    catch
        return nothing
    end
    return nothing
end

function Mycelia.kernel_batched_benchmark(;
        n_reads::Integer = 4096,
        n_states::Integer = 256,
        avg_degree::Integer = 4,
        n_depths::Integer = 64,
        iterations::Integer = 3,
        cuda_backend = nothing,
        seed::Integer = 20260707)
    rng = Random.MersenneTwister(seed)

    # Synthetic read-only shared graph (CSR) + emission table + start frontier.
    edge_offsets = Vector{Int32}(undef, n_states + 1)
    edge_offsets[1] = 0
    edge_dst_list = Int32[]
    edge_logprob_list = Float64[]
    for sid in 1:n_states
        deg = max(1, round(Int, avg_degree))
        edge_offsets[sid + 1] = edge_offsets[sid] + deg
        for _ in 1:deg
            push!(edge_dst_list, Int32(rand(rng, 1:n_states)))
            push!(edge_logprob_list, log(rand(rng)))
        end
    end
    edge_dst = collect(edge_dst_list)
    edge_logprob = collect(edge_logprob_list)
    emission = log.(rand(rng, n_reads, n_states, n_depths))
    current0 = fill(-Inf, n_reads, n_states)
    for r in 1:n_reads
        current0[r, rand(rng, 1:n_states)] = 0.0  # one live start state per read
    end

    edge_relaxations = Int(edge_offsets[end]) * n_reads * n_depths

    function _time_backend(backend)
        # Warmup (compile + allocate), then time `iterations` full sweeps.
        _frontier_sweep!(backend, current0, edge_offsets, edge_dst, edge_logprob, emission)
        best = Inf
        result = nothing
        for _ in 1:iterations
            t0 = Base.time_ns()
            result = _frontier_sweep!(
                backend, current0, edge_offsets, edge_dst, edge_logprob, emission)
            elapsed = (Base.time_ns() - t0) / 1e9
            best = min(best, elapsed)
        end
        throughput = edge_relaxations / best
        return (; seconds = best, edge_relaxations_per_sec = throughput, result = result)
    end

    cpu = _time_backend(KA.CPU())

    gpu_backend = cuda_backend === nothing ? _detect_cuda_backend() : cuda_backend
    if gpu_backend === nothing
        return (;
            cuda_available = false,
            n_reads = Int(n_reads), n_states = Int(n_states),
            avg_degree = Int(avg_degree), n_depths = Int(n_depths),
            edge_relaxations = edge_relaxations,
            cpu_seconds = cpu.seconds,
            cpu_edge_relaxations_per_sec = cpu.edge_relaxations_per_sec,
            gpu_seconds = nothing,
            gpu_edge_relaxations_per_sec = nothing,
            gpu_speedup = nothing,
            results_match = nothing,
            note = "No functional CUDA backend detected; CPU-only. `import CUDA` " *
                   "on a GPU host (e.g. Lovelace) to enable the GPU benchmark.")
    end

    gpu = _time_backend(gpu_backend)
    # Cross-backend sanity: the GPU sweep must reproduce the CPU sweep bit-for-bit
    # (the kernel is deterministic and identical across backends).
    results_match = cpu.result == gpu.result
    return (;
        cuda_available = true,
        n_reads = Int(n_reads), n_states = Int(n_states),
        avg_degree = Int(avg_degree), n_depths = Int(n_depths),
        edge_relaxations = edge_relaxations,
        cpu_seconds = cpu.seconds,
        cpu_edge_relaxations_per_sec = cpu.edge_relaxations_per_sec,
        gpu_seconds = gpu.seconds,
        gpu_edge_relaxations_per_sec = gpu.edge_relaxations_per_sec,
        gpu_speedup = cpu.seconds / gpu.seconds,
        results_match = results_match,
        backend = string(nameof(typeof(gpu_backend))))
end

end # module MyceliaKernelAbstractionsExt
