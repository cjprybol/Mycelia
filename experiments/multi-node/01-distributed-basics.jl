# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Julia 1.10
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Distributed.jl Basics — Parallel K-mer Counting
#
# **Learning objectives:**
# - `addprocs()`, `@everywhere`, `@distributed`, `pmap`
# - `ClusterManagers.jl` + `SlurmManager(N)` for SLURM integration
# - SharedArrays for shared-memory parallelism
#
# **Target system:** Lawrencium (2-4 nodes)
# **Week 1-2 of learning path**

# %% Setup
import Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=joinpath(homedir(), "workspace", "Mycelia"))
Pkg.add(["Distributed", "ClusterManagers", "SharedArrays", "BenchmarkTools"])

# %%
import Distributed
import ClusterManagers
import SharedArrays
import BenchmarkTools

# %% [markdown]
# ## 1. Basic `addprocs` — Local Workers
#
# Start with local workers to understand the API before scaling to SLURM nodes.

# %%
# Add local workers (adjust nprocs based on available cores)
Distributed.addprocs(4)
println("Workers: ", Distributed.workers())
println("Number of workers: ", Distributed.nworkers())

# %% [markdown]
# ## 2. `@everywhere` — Loading Code on All Workers

# %%
Distributed.@everywhere begin
    import Mycelia
end

# %% [markdown]
# ## 3. `pmap` — Parallel Map for K-mer Counting
#
# K-mer counting across multiple genomes is embarrassingly parallel —
# each genome can be processed independently.

# %%
# Example: parallel k-mer counting across genome files
# genome_files = readdir(joinpath(homedir(), "workspace/data/genomes"); join=true)
# results = Distributed.pmap(genome_files) do genome_file
#     Mycelia.count_kmers(genome_file; k=21)
# end

# %% [markdown]
# ## 4. `@distributed` — Reduction Operations

# %%
# Example: distributed reduction to count total k-mers
# total = Distributed.@distributed (+) for f in genome_files
#     length(Mycelia.count_kmers(f; k=21))
# end

# %% [markdown]
# ## 5. SLURM Integration with ClusterManagers.jl
#
# On Lawrencium, use `SlurmManager` to add workers across SLURM-allocated nodes.

# %%
# NOTE: Only run this on Lawrencium within a SLURM allocation
# salloc --nodes=2 --ntasks-per-node=16 --time=01:00:00 --partition=lr4 --qos=lr_normal --account=<account>
#
# using ClusterManagers
# addprocs(SlurmManager(32))  # 2 nodes × 16 tasks

# %% [markdown]
# ## 6. Benchmarking: Serial vs Parallel

# %%
# BenchmarkTools.@btime begin
#     # Serial
#     results_serial = map(genome_files) do f
#         Mycelia.count_kmers(f; k=21)
#     end
# end

# %%
# BenchmarkTools.@btime begin
#     # Parallel
#     results_parallel = Distributed.pmap(genome_files) do f
#         Mycelia.count_kmers(f; k=21)
#     end
# end

# %% [markdown]
# ## Cleanup

# %%
Distributed.rmprocs(Distributed.workers())
