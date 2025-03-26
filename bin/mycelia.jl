#!/usr/bin/env julia

"""
    mycelia.jl

Command line interface for the Mycelia package.
This script is designed to work when placed in a bin/ directory of a project.
"""

# Ensure we have the correct project activated
const BIN_DIR = dirname(abspath(PROGRAM_FILE))
const PROJECT_ROOT = dirname(BIN_DIR)

# Add the project directory to the load path and activate it
if isfile(joinpath(PROJECT_ROOT, "Project.toml"))
    import Pkg
    Pkg.activate(PROJECT_ROOT)
end

using ArgParse
import Pkg

# Import the package. The commented line shows how to add the src directory to the load path
# if needed during development
# push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using Mycelia

"""
    parse_commandline()

Parse command line arguments for the Mycelia CLI tool.
"""
function parse_commandline()
    settings = ArgParseSettings(
        description = "Mycelia: A tool for genome graph construction and manipulation",
        version = Pkg.TOML.parsefile(joinpath(pkgdir(Mycelia), "Project.toml"))["version"],
        add_version = true
    )

    @add_arg_table settings begin
        "construct"
            help = "construct a reference graph"
            action = :command
        "assemble"
            help = "assemble a sequence graph"
            action = :command
        "convert"
            help = "convert between graph formats"
            action = :command
        "analyze"
            help = "analyze graphs or sequences"
            action = :command
        "visualize"
            help = "visualize genome graphs"
            action = :command
        "test"
            help = "run tests"
            action = :command
    end
    
    # Construct command
    settings["construct"].description = "Construct a reference kmer graph from observations"
    @add_arg_table settings["construct"] begin
        "--k"
            help = "kmer size that will be used. Must be prime and < 64"
            arg_type = Int
            required = true
        "--out"
            help = "outpath for the genome graph that will be constructed"
            arg_type = String
            required = true
        "--fastx"
            help = "one or more fasta or fastq file(s). Can be gzipped as long as files end in .gz"
            arg_type = String
            required = true
            nargs = '*'
        "--min-count"
            help = "minimum kmer count threshold"
            arg_type = Int
            default = 3
        "--threads"
            help = "number of threads to use"
            arg_type = Int
            default = Threads.nthreads()
    end
    
    # Assemble command
    settings["assemble"].description = "Assemble graph genomes from observations"
    @add_arg_table settings["assemble"] begin
        "--kmers"
            help = "sorted file of trusted kmers"
            arg_type = String
            required = true
        "--sequences"
            help = """
                  One or more paths to fasta or fastq files.

                  Sequences will be used for determining edges between kmers
                  and depth of coverage for kmers and resulting segments in the
                  GFA graph
                  """
            arg_type = String
            nargs = '*'
        "--colors"
            help = """
                  Paths to files that represent the color groups for the colored debruijn graphs.

                  Each file contains a list of fastq and/or fasta files. Segments will be colored
                  using the tag CL:Z:# in the output GFA file, where # is an integer representing
                  the order of the color groups as they were presented in the command line.
                  """
            arg_type = String
            nargs = '*'
        "--out"
            help = "output file path (GFA format by default)"
            arg_type = String
            required = true
        "--threads"
            help = "number of threads to use"
            arg_type = Int
            default = Threads.nthreads()
    end
   
    # Convert command
    settings["convert"].description = "Interconvert between different graph formats (.jld2, .gfa, .neo4j)"
    @add_arg_table settings["convert"] begin
        "--in"
            help = "input graph"
            arg_type = String
            required = true
        "--out"
            help = "outpath of the reformatted genome graph"
            arg_type = String
            required = true
        "--username"
            help = "neo4j username (required if outformat is neo4j)"
            arg_type = String
        "--password"
            help = "neo4j password (required if outformat is neo4j)"
    end

    # Analyze command (placeholder)
    settings["analyze"].description = "Analyze genome graphs or sequences"
    @add_arg_table settings["analyze"] begin
        "--in"
            help = "input graph or sequence file"
            arg_type = String
            required = true
        "--out"
            help = "output file for analysis results"
            arg_type = String
            required = true
        "--type"
            help = "type of analysis to run (metrics, stats, etc.)"
            arg_type = String
            required = true
    end

    # Visualize command (placeholder)
    settings["visualize"].description = "Visualize genome graphs in various formats"
    @add_arg_table settings["visualize"] begin
        "--in"
            help = "input graph file"
            arg_type = String
            required = true
        "--out"
            help = "output file for visualization"
            arg_type = String
            required = true
        "--format"
            help = "output visualization format (svg, png, etc.)"
            arg_type = String
            default = "svg"
        "--layout"
            help = "graph layout algorithm to use"
            arg_type = String
            default = "force"
    end

    # Test command
    settings["test"].description = "Run test suite for Mycelia"
    @add_arg_table settings["test"] begin
        "--coverage"
            help = "generate test coverage report"
            action = :store_true
    end

    return parse_args(settings)
end

"""
    dispatch_command(command, args)

Dispatch the appropriate function based on the command.
"""
function dispatch_command(command, args)
    command_functions = Dict(
        "construct" => Mycelia.construct,
        "assemble" => Mycelia.assemble,
        "convert" => Mycelia.convert_graph,
        "analyze" => Mycelia.analyze,
        "visualize" => Mycelia.visualize,
        "test" => Mycelia.run_tests
    )
    
    if haskey(command_functions, command)
        # Call the appropriate function with args
        return command_functions[command](args)
    else
        error("Unknown command: $command")
    end
end

"""
    main()

Main entry point for the Mycelia CLI tool.
"""
function main()
    parsed_args = parse_commandline()
    
    # Get the command (first positional argument)
    command = parsed_args["%COMMAND%"]
    if command === nothing
        println("No command specified. Use --help for usage information.")
        return 1
    end
    
    # Dispatch to the appropriate command handler
    try
        return dispatch_command(command, parsed_args[command])
    catch e
        @error "Error executing command" exception=(e, catch_backtrace())
        return 1
    end
end

# Only run main when script is executed directly, not when included
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end