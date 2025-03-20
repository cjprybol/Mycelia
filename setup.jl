#!/usr/bin/env julia

"""
    setup.jl

Setup script for the Mycelia CLI tool.
This script should be run from the project root directory.

Usage:
    julia setup.jl install        # Install dependencies and prepare the CLI
    julia setup.jl build          # Build all CLI optimized versions
    julia setup.jl test           # Run CLI tests
"""

using Pkg

function install_dependencies()
    @info "Installing Mycelia CLI dependencies..."
    Pkg.activate(".")
    Pkg.add(["ArgParse", "PackageCompiler"])
    
    # Make sure the bin directory exists
    mkpath("bin")
    
    # Ensure scripts are executable
    if Sys.isunix()
        try
            run(`chmod +x bin/build_mycelia.jl bin/mycelia.jl`)
        catch
            @warn "Could not set execute permissions. This is fine if the files don't exist yet."
        end
    end
    
    @info "Dependencies installed successfully."
end

function build_cli()
    @info "Building Mycelia CLI..."
    
    # Check if build script exists
    if !isfile("bin/build_mycelia.jl")
        @error "build_mycelia.jl not found in bin directory. Run setup.jl install first."
        return 1
    end
    
    # Build both system image and app
    cmd = `julia bin/build_mycelia.jl all`
    run(cmd)
    
    @info "Mycelia CLI built successfully."
end

function run_tests()
    @info "Running Mycelia CLI tests..."
    
    # Check if mycelia script exists
    if !isfile("bin/mycelia")
        @error "mycelia executable not found in bin directory. Run setup.jl build first."
        return 1
    end
    
    # Run the test command
    cmd = `bin/mycelia test`
    run(cmd)
    
    @info "Tests completed."
end

function main()
    if length(ARGS) < 1
        println("""
        Usage:
            julia setup.jl install  # Install dependencies and prepare the CLI
            julia setup.jl build    # Build all CLI optimized versions
            julia setup.jl test     # Run CLI tests
        """)
        return 1
    end
    
    cmd = lowercase(ARGS[1])
    
    if cmd == "install"
        install_dependencies()
    elseif cmd == "build"
        build_cli()
    elseif cmd == "test"
        run_tests()
    else
        println("Unknown command: $cmd")
        return 1
    end
    
    return 0
end

# Run the main function only when executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end