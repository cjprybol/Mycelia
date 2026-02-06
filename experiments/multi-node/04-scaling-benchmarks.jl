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
# # Multi-Node Scaling Benchmarks
#
# **Benchmark design:**
# - Scaling study: 1, 2, 4, 8, 16, 32 nodes
# - Workloads: k-mer counting, graph construction, PageRank
# - Metrics: wall time, parallel efficiency, communication overhead
# - Datasets: small (1 genome), medium (100), large (10,000+ RefSeq)
#
# **Target system:** NERSC Perlmutter (Week 5-8)

# %% Setup
import Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=joinpath(homedir(), "workspace", "Mycelia"))
Pkg.add(["MPI", "DataFrames", "CSV", "CairoMakie", "BenchmarkTools", "Dates"])

# %%
import MPI
import DataFrames
import CSV
import CairoMakie
import BenchmarkTools
import Dates
import Mycelia

# %% [markdown]
# ## Benchmark Configuration

# %%
const BENCHMARK_CONFIG = (
    node_counts = [1, 2, 4, 8, 16, 32],
    workloads = [:kmer_counting, :graph_construction, :pagerank],
    datasets = Dict(
        :small => "1 genome",
        :medium => "100 genomes",
        :large => "10,000+ RefSeq genomes"
    ),
    k = 21,
    pagerank_iterations = 100,
    repeats = 3
)

# %% [markdown]
# ## 1. K-mer Counting Scaling

# %%
# function benchmark_kmer_counting(genome_files, n_nodes; k=21)
#     MPI.Init()
#     comm = MPI.COMM_WORLD
#     rank = MPI.Comm_rank(comm)
#
#     t_start = MPI.Wtime()
#
#     # Partition genomes across ranks
#     my_files = partition_files(genome_files, rank, n_nodes)
#
#     # Count k-mers locally
#     local_counts = Dict{String,Int}()
#     for f in my_files
#         merge!(+, local_counts, count_kmers(f; k=k))
#     end
#
#     # Reduce to rank 0
#     # (In practice, use custom MPI reduction for dictionaries)
#
#     t_end = MPI.Wtime()
#     wall_time = t_end - t_start
#
#     if rank == 0
#         return (n_nodes=n_nodes, wall_time=wall_time, n_genomes=length(genome_files))
#     end
#
#     MPI.Finalize()
# end

# %% [markdown]
# ## 2. Results Analysis Template

# %%
# After running benchmarks, analyze results:
#
# results = DataFrames.DataFrame(
#     n_nodes = Int[],
#     workload = Symbol[],
#     dataset = Symbol[],
#     wall_time = Float64[],
#     speedup = Float64[],
#     efficiency = Float64[]
# )
#
# # Calculate parallel efficiency
# for row in DataFrames.eachrow(results)
#     baseline = results[results.n_nodes .== 1 .& results.workload .== row.workload .& results.dataset .== row.dataset, :wall_time][1]
#     row.speedup = baseline / row.wall_time
#     row.efficiency = row.speedup / row.n_nodes
# end

# %% [markdown]
# ## 3. Scaling Plots

# %%
# function plot_scaling(results; workload=:kmer_counting)
#     df = DataFrames.filter(r -> r.workload == workload, results)
#
#     fig = CairoMakie.Figure(; size=(800, 600))
#
#     # Strong scaling plot
#     ax1 = CairoMakie.Axis(fig[1, 1];
#         xlabel="Number of Nodes",
#         ylabel="Speedup",
#         title="Strong Scaling: $workload"
#     )
#
#     for dataset in unique(df.dataset)
#         sub = DataFrames.filter(r -> r.dataset == dataset, df)
#         CairoMakie.lines!(ax1, sub.n_nodes, sub.speedup; label=String(dataset))
#     end
#
#     # Ideal scaling line
#     CairoMakie.lines!(ax1, [1, 32], [1, 32]; linestyle=:dash, color=:gray, label="Ideal")
#     CairoMakie.axislegend(ax1)
#
#     # Efficiency plot
#     ax2 = CairoMakie.Axis(fig[1, 2];
#         xlabel="Number of Nodes",
#         ylabel="Parallel Efficiency",
#         title="Efficiency: $workload"
#     )
#
#     for dataset in unique(df.dataset)
#         sub = DataFrames.filter(r -> r.dataset == dataset, df)
#         CairoMakie.lines!(ax2, sub.n_nodes, sub.efficiency; label=String(dataset))
#     end
#
#     CairoMakie.hlines!(ax2, [1.0]; linestyle=:dash, color=:gray)
#     CairoMakie.axislegend(ax2)
#
#     # Save both formats
#     CairoMakie.save("scaling_$(workload).svg", fig)
#     CairoMakie.save("scaling_$(workload).png", fig)
#
#     return fig
# end

# %% [markdown]
# ## 4. SLURM Submission Template for Scaling Study
#
# ```bash
# #!/bin/bash
# # scaling-study.sh — submit benchmarks at each node count
#
# for NODES in 1 2 4 8 16 32; do
#   sbatch --nodes=$NODES \
#     --ntasks-per-node=4 \
#     --cpus-per-task=32 \
#     --constraint=cpu \
#     --qos=regular \
#     --time=02:00:00 \
#     --job-name="bench-${NODES}nodes" \
#     --output="bench-${NODES}nodes-%j.out" \
#     --wrap="srun julia --project 04-scaling-benchmarks.jl --nodes=$NODES"
# done
# ```

# %% [markdown]
# ## Learning Progression
#
# | Week | Topic                        | System      |
# | ---- | ---------------------------- | ----------- |
# | 1    | Distributed.jl + pmap        | Lawrencium  |
# | 2    | ClusterManagers.jl + SLURM   | Lawrencium  |
# | 3    | MPI.jl basics                | Lawrencium  |
# | 4    | Distributed graph algorithms | Lawrencium  |
# | 5-6  | Scale to NERSC, benchmark    | NERSC       |
# | 7-8  | Integrate into Mycelia       | All systems |
