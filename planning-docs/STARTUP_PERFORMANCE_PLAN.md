# Startup & Precompilation Performance Plan

Goal: cut time-to-first-analysis for Mycelia by enabling precompilation, avoiding eager loading of heavy stacks, and providing reproducible warm-start workflows for temporary environments.

## Current Observations
- Precompilation explicitly disabled in `src/Mycelia.jl`.
- Eager top-level imports pull in Makie/Plots/GraphMakie/Luxor/SankeyPlots, Conda/HTTP, POMDPs/ReinforcementLearning, etc., even when not used.
- Conda/tool bootstrap helpers exist; must stay off the import path to avoid side effects during precompile.

## Prioritized Work (impact/effort)
1) Re-enable package precompilation and add a minimal precompile workload to cover the “hello world” flow (load, tiny graph build, one plot save).  
2) Lazily load heavy/optional stacks (Makie/Plots and RL ecosystems, Conda-bound helpers) via extensions or runtime `Base.require` so base load avoids them.  
3) Add a warm-cache script/CI step to run `Pkg.instantiate` + `Pkg.precompile` and stash the depot tarball for new machines/temp envs.  
4) Reduce duplicate plotting backends (prefer Makie or Plots, not both) and trim unused deps after confirming usage.  
5) Add docs for temp-env workflows (shared depot, toggles, and how to skip optional stacks).

## Execution Checklist
- [ ] Enable `__precompile__()` in `src/Mycelia.jl` (with note about disabling if specific tests require).  
- [ ] Add `PrecompileTools` to `[deps]` and create `src/precompile_workload.jl` with small fixtures-driven workload; include from `Mycelia.jl`.  
- [ ] Convert plotting/RL/Conda-heavy code to optional extensions; move imports inside functions or extension modules.  
- [ ] Add `bin/warm_cache.jl` (instantiate + precompile) and document usage; wire optional CI artifact publishing.  
- [ ] Audit and prune unused deps (or guard them) to shrink precompile surface.  
- [ ] Update README/docs with startup tips (shared depot, env toggles, how to disable optional stacks).  

## Notes / Risks
- Watch for the previous “testing framework fix” that led to disabling precompile; rerun test suite after re-enabling.  
- Keep Conda/tool downloads behind explicit user calls; precompilation must remain side-effect free.  
- Extensions require Julia ≥1.9; ensure compat remains satisfied.
