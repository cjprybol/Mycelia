# Batched array-frontier Viterbi corrector — backend-agnostic GPU kernel (Phase B,
# bead td-qoo3). This file lives in the BASE package and declares only the public
# entry points; the actual KernelAbstractions kernel + backend-agnostic decoder,
# equivalence oracle, and GPU benchmark live in the package EXTENSION
# `ext/MyceliaKernelAbstractionsExt.jl`, which is loaded lazily when (and only
# when) the user does `import KernelAbstractions`.
#
# Why an extension (not a hard dep)?  KernelAbstractions.jl is a weakdep in
# `Project.toml`, so the base Mycelia package carries NO GPU dependency: the CPU
# array-frontier PoC (`batched-viterbi-poc.jl`) and every existing caller stay
# exactly as they are. The same `@kernel` runs on the CPU backend (used by the
# equivalence-oracle test so CI needs no GPU) AND on a CUDA/GPU backend (used by
# the benchmark on GPU hosts like Lovelace). CUDA itself is never a dependency —
# the GPU path is discovered at runtime only if the user has already imported a
# functional CUDA and passes its backend.
#
# Phase B design (built on the merged CPU PoC, byte-identical by construction):
#
#   * The wavefront is a `@kernel` that advances the dense `n_reads x n_states`
#     frontier ONE depth-step for every read in a length-bin in parallel — one
#     thread per read (the ADR's zero-write-contention, fully-data-parallel read
#     axis). Each thread runs the SAME sequential inner scan over states/out-edges
#     as the CPU PoC (same order, same `(state + log_transition) + emission`
#     associativity, same strict-`>` relax/tie-break), so the produced
#     `next_scores`/`next_predecessor` are bit-identical to the PoC inner loop.
#
#   * The graph is uploaded ONCE as read-only structure-of-arrays: a CSR
#     transition layout (`edge_offsets`, `edge_dst`, `edge_logprob`) plus a
#     precomputed emission table (`emission[read, state, depth]`) whose entries
#     come from the SAME `_call_viterbi_state_emission_logp` the scalar/PoC
#     decoder calls — so every value the kernel reads is byte-identical to the
#     scalar path.
#
#   * Start-init, per-depth predecessor materialization + best-state selection,
#     and path reconstruction stay host-side and reuse the scalar tie-break
#     helpers verbatim (`_best_correction_state`, `_best_correction_target_state`,
#     `_reconstruct_correction_path`), exactly as the CPU PoC does. Only the hot
#     frontier relaxation moves into the kernel.
#
# The three entry points below are declared here (with a helpful fallback that
# fires when the extension is not loaded) and given real methods by the
# extension:
#
#   * `kernel_batched_correct_observations` — decode a batch through the kernel.
#   * `kernel_batched_equivalence_oracle`   — assert byte-identity vs the scalar
#     `correct_observations` AND the CPU `batched_correct_observations` PoC.
#   * `kernel_batched_benchmark`            — throughput microbenchmark of the
#     frontier kernel on the CPU backend and, if a functional CUDA backend is
#     supplied/detected, on the GPU; degrades gracefully to CPU-only.

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Decode a batch of `observations` through the backend-agnostic frontier kernel
(Phase B of bead td-qoo3). Requires the KernelAbstractions extension — do
`import KernelAbstractions` first. See `ext/MyceliaKernelAbstractionsExt.jl`.

Keyword `backend` selects the KernelAbstractions backend (defaults to the CPU
backend); pass a `CUDA.CUDABackend()` to run on GPU. Returns a
[`ViterbiCorrectionResult`](@ref) identical in shape to
[`batched_correct_observations`](@ref).
"""
function kernel_batched_correct_observations end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Equivalence oracle for the Phase B kernel decoder: runs the scalar
[`correct_observations`](@ref), the CPU PoC [`batched_correct_observations`](@ref),
and the kernel [`kernel_batched_correct_observations`](@ref) on the same inputs
and checks that the kernel's corrected paths + log-probability scores are
byte-identical to BOTH (`===` on scores, label/step/target equality per read).
Runs on the CPU KernelAbstractions backend by default so it works with no GPU
(CI-safe). Requires `import KernelAbstractions`.
"""
function kernel_batched_equivalence_oracle end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Throughput microbenchmark of the frontier kernel. Always benchmarks the CPU
KernelAbstractions backend; if a functional CUDA backend is available (either
passed explicitly or auto-detected from an already-imported, `CUDA.functional()`
CUDA) it also benchmarks the GPU and reports the speedup. Degrades gracefully to
CPU-only when no GPU is present. Requires `import KernelAbstractions`.
"""
function kernel_batched_benchmark end

# Friendly fallback: these fire only when the KernelAbstractions extension has
# NOT been loaded (the extension installs strictly-more-specific typed methods
# that take precedence once `import KernelAbstractions` has run).
const _KERNEL_EXT_HINT = "requires the KernelAbstractions extension; run " *
    "`import KernelAbstractions` to load `MyceliaKernelAbstractionsExt`."

kernel_batched_correct_observations(args...; kwargs...) =
    error("kernel_batched_correct_observations " * _KERNEL_EXT_HINT)
kernel_batched_equivalence_oracle(args...; kwargs...) =
    error("kernel_batched_equivalence_oracle " * _KERNEL_EXT_HINT)
kernel_batched_benchmark(args...; kwargs...) =
    error("kernel_batched_benchmark " * _KERNEL_EXT_HINT)
