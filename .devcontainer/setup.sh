#!/bin/bash
export LD_LIBRARY_PATH="" && julia -e 'import Pkg; Pkg.add(["IJulia"]); using IJulia; n=Sys.CPU_THREADS; ks=[2^i for i in 0:floor(Int,log2(n-1)) if 2^i < n]; ks = push!(ks, n); foreach(k -> IJulia.installkernel("julia-$k", env=Dict("JULIA_NUM_THREADS"=>"$k", "LD_LIBRARY_PATH"=>"")), ks)'
npm install -g @anthropic-ai/claude-code@latest
npm install -g @openai/codex@latest
# alternate one-liner installation
# npm install -g @anthropic-ai/claude-code @github/copilot @openai/codex @google/gemini-cli

# if installing in an environment without sudo, use:
# install nvm (node version manager)
# curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
# nvm install --lts
# nvm use --lts
# install non-globally for just my personal user
# npm install @anthropic-ai/claude-code@latest
# npm install @openai/codex@latest
# need to run with npx prefix when not installed globally
# npx codex
# npx claude

# Note - running codex on extra high (xhigh) is painfully slow, prefer high or default (medium)
# If xhigh is equivalent to pro level in main app, use pro in app to plan and then switch to high in codex for faster response times.