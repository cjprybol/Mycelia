julia --code-coverage test/runtests.jl
julia assess-code-coverage.jl
julia assess-benchmarks.jl