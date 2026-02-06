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
# # MPI.jl Basics — Inter-Node Communication
#
# **Learning objectives:**
# - `MPI.Init()`, `Comm_rank`, `Comm_size`
# - Point-to-point: `Send`, `Recv`
# - Collective: `Bcast`, `Reduce`, `Allgather`
# - NERSC specifics: `srun` instead of `mpirun`, module loading
#
# **Exercise:** Distributed graph construction — each rank builds a subgraph
# from a genome partition, then merge via `Allgather`.
#
# **Target system:** Lawrencium (Week 3-4)
#
# **NOTE:** MPI scripts must be run as standalone scripts, not interactively.
# Use: `srun -n 4 julia 02-mpi-basics.jl`

# %% Setup (run interactively for package installation only)
import Pkg
Pkg.activate(; temp=true)
Pkg.add(["MPI"])

# %%
import MPI

# %% [markdown]
# ## 1. Hello World — Rank and Size
#
# Each MPI process has a unique rank (0 to N-1) and knows the total size.

# %%
# This block should be saved to a standalone script and run with srun/mpirun
#
# MPI.Init()
# comm = MPI.COMM_WORLD
# rank = MPI.Comm_rank(comm)
# size = MPI.Comm_size(comm)
# println("Hello from rank $rank of $size")
# MPI.Finalize()

# %% [markdown]
# ## 2. Point-to-Point: Send and Recv
#
# Rank 0 sends data to rank 1.

# %%
# MPI.Init()
# comm = MPI.COMM_WORLD
# rank = MPI.Comm_rank(comm)
#
# if rank == 0
#     data = [1.0, 2.0, 3.0, 4.0]
#     MPI.Send(data, comm; dest=1, tag=0)
#     println("Rank 0 sent: $data")
# elseif rank == 1
#     data = zeros(4)
#     MPI.Recv!(data, comm; source=0, tag=0)
#     println("Rank 1 received: $data")
# end
#
# MPI.Finalize()

# %% [markdown]
# ## 3. Collective: Broadcast
#
# Rank 0 broadcasts a value to all ranks.

# %%
# MPI.Init()
# comm = MPI.COMM_WORLD
# rank = MPI.Comm_rank(comm)
#
# root = 0
# if rank == root
#     data = Float64[42.0]
# else
#     data = Float64[0.0]
# end
#
# MPI.Bcast!(data, comm; root=root)
# println("Rank $rank has data: $data")
#
# MPI.Finalize()

# %% [markdown]
# ## 4. Collective: Reduce
#
# Sum values from all ranks to rank 0.

# %%
# MPI.Init()
# comm = MPI.COMM_WORLD
# rank = MPI.Comm_rank(comm)
#
# local_value = Float64[rank + 1.0]  # rank 0 → 1.0, rank 1 → 2.0, etc.
# global_sum = Float64[0.0]
#
# MPI.Reduce!(local_value, global_sum, MPI.SUM, comm; root=0)
# if rank == 0
#     println("Global sum: $(global_sum[1])")
# end
#
# MPI.Finalize()

# %% [markdown]
# ## 5. Allgather — Distributed Graph Construction
#
# Each rank builds a local subgraph from its genome partition,
# then all ranks gather the complete graph.

# %%
# MPI.Init()
# comm = MPI.COMM_WORLD
# rank = MPI.Comm_rank(comm)
# size = MPI.Comm_size(comm)
#
# # Each rank processes a subset of genomes
# all_genomes = readdir(genome_dir; join=true)
# chunk_size = ceil(Int, length(all_genomes) / size)
# my_start = rank * chunk_size + 1
# my_end = min((rank + 1) * chunk_size, length(all_genomes))
# my_genomes = all_genomes[my_start:my_end]
#
# # Build local subgraph
# local_edges = build_subgraph(my_genomes)
#
# # Gather all edges to all ranks
# all_edges = MPI.Allgather(local_edges, comm)
#
# # Merge into complete graph
# complete_graph = merge_edges(all_edges)
#
# MPI.Finalize()

# %% [markdown]
# ## NERSC-Specific Notes
#
# On Perlmutter, use `srun` instead of `mpirun`:
# ```bash
# module load julia/1.10.10
# srun -n 4 -c 32 julia --project 02-mpi-basics.jl
# ```
#
# SLURM script:
# ```bash
# #!/bin/bash
# #SBATCH --nodes=2
# #SBATCH --ntasks-per-node=2
# #SBATCH --cpus-per-task=32
# #SBATCH --constraint=cpu
# #SBATCH --qos=regular
# #SBATCH --time=00:30:00
#
# module load julia/1.10.10
# srun julia --project 02-mpi-basics.jl
# ```
