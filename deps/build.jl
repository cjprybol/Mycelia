import Pkg

# Force rebuild of Conda to ensure proper setup
println("Rebuilding Conda.jl to ensure proper conda environment setup...")
try
    Pkg.build("Conda")
    println("Conda.jl rebuild successful")
catch e
    @warn "Conda.jl rebuild failed, but continuing installation" exception=e
end