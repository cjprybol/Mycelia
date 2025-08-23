#!/bin/bash
export LD_LIBRARY_PATH="" && julia -e 'import Pkg; Pkg.add(["IJulia"]); using IJulia; n=Sys.CPU_THREADS; ks=[2^i for i in 0:floor(Int,log2(n-1)) if 2^i < n]; ks = push!(ks, n); foreach(k -> IJulia.installkernel("julia-$k", env=Dict("JULIA_NUM_THREADS"=>"$k", "LD_LIBRARY_PATH"=>"")), ks)'
npm install -g @anthropic-ai/claude-code
npm install -g @google/gemini-cli