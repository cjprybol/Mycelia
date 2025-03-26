#!/usr/bin/env julia

"""
    build_mycelia.jl

Build script for Mycelia CLI that can create both a system image and a standalone application.
This script is designed to work when placed in a bin/ directory of a project.

Run with:
    julia bin/build_mycelia.jl sysimage  # To build just the system image
    julia bin/build_mycelia.jl app       # To build just the standalone app
    julia bin/build_mycelia.jl all       # To build both
"""

using Pkg

# Determine the project root and bin directory
const BIN_DIR = dirname(abspath(PROGRAM_FILE))
const PROJECT_ROOT = dirname(BIN_DIR)

# Ensure we're activating the project at the root level
if isfile(joinpath(PROJECT_ROOT, "Project.toml"))
    @info "Activating project at $(PROJECT_ROOT)"
    Pkg.activate(PROJECT_ROOT)
else
    @warn "No Project.toml found in $(PROJECT_ROOT). Creating a new project."
    Pkg.activate(PROJECT_ROOT; shared=false)
end

# Make sure PackageCompiler is installed
if !haskey(Pkg.project().dependencies, "PackageCompiler")
    @info "Installing PackageCompiler..."
    Pkg.add("PackageCompiler")
end

using PackageCompiler

# Common configuration
const MYCELIA_SCRIPT = joinpath(BIN_DIR, "mycelia.jl")
const PRECOMPILE_FILE = joinpath(BIN_DIR, "precompile_statements.jl")
const PACKAGES = ["Mycelia", "ArgParse"]
const APP_NAME = "mycelia"
const CPU_TARGET = "native"

"""
    build_sysimage()

Build a system image for Mycelia CLI.
"""
function build_sysimage()
    @info "Building system image for Mycelia..."
    
    sysimage_path = joinpath(BIN_DIR, "$(APP_NAME).so")
    
    create_sysimage(
        PACKAGES,
        sysimage_path = sysimage_path,
        precompile_execution_file = MYCELIA_SCRIPT,
        precompile_statements_file = PRECOMPILE_FILE,
        cpu_target = CPU_TARGET
    )
    
    # Create a convenience script in the bin directory
    wrapper_path = joinpath(BIN_DIR, APP_NAME)
    open(wrapper_path, "w") do io
        script = """
        #!/bin/bash
        SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
        julia --project="\$(dirname \$SCRIPT_DIR)" --sysimage "\${SCRIPT_DIR}/$(APP_NAME).so" "\${SCRIPT_DIR}/mycelia.jl" "\$@"
        """
        write(io, script)
    end
    chmod(wrapper_path, 0o755)
    
    @info "System image created as $(sysimage_path)"
    @info "Wrapper script created as $(wrapper_path)"
    @info "You can now run the application with: $(wrapper_path)"
    
    return sysimage_path
end

"""
    build_app()

Build a standalone application for Mycelia CLI.
"""
function build_app()
    @info "Building standalone application for Mycelia..."
    
    app_dir = joinpath(BIN_DIR, "$(APP_NAME)_app")
    
    # We need to create the app from the project root
    create_app(
        PROJECT_ROOT,   # Source directory (project root)
        app_dir,        # Destination directory (in bin/)
        executables = [
            ExecutableConfig(
                MYCELIA_SCRIPT, 
                APP_NAME; 
                precompile_execution_file = MYCELIA_SCRIPT,
                precompile_statements_file = PRECOMPILE_FILE
            )
        ],
        force = true,
        include_lazy_artifacts = true,
        cpu_target = CPU_TARGET
    )
    
    # Create a symlink for easier access
    symlink_path = joinpath(BIN_DIR, "$(APP_NAME)-app")
    islink(symlink_path) && rm(symlink_path)
    symlink(joinpath(app_dir, "bin", APP_NAME), symlink_path)
    
    @info "Standalone application created at: $(app_dir)"
    @info "Symlink created at: $(symlink_path)"
    @info "You can now run the standalone app with: $(symlink_path)"
    
    return app_dir
end

"""
    main()

Main entry point for the build script.
"""
function main()
    if length(ARGS) < 1
        println("""
        Usage:
            julia bin/build_mycelia.jl sysimage  # Build system image only
            julia bin/build_mycelia.jl app       # Build standalone app only
            julia bin/build_mycelia.jl all       # Build both
        """)
        return 1
    end
    
    cmd = lowercase(ARGS[1])
    
    if cmd == "sysimage"
        build_sysimage()
    elseif cmd == "app"
        build_app()
    elseif cmd == "all"
        build_sysimage()
        build_app()
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