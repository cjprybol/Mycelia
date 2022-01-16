#!/usr/bin/env julia
# push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using ArgParse
using Mycelia
import Pkg

function parse_commandline()
    settings = ArgParseSettings()
    settings.description = "description"
    settings.version = Pkg.TOML.parsefile(joinpath(pkgdir(Mycelia), "Project.toml"))["version"]
    settings.add_version = true

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
       "test"
           help = "run tests"
           action = :command
    end
    
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
    end
    
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
    

#     settings["plot"].description = ""
#     @add_arg_table settings["plot"] begin
#         "histogram"
#             help = "plot a histogram of kmer counts"
#             action = :command
#         "rank-frequency"
#             help = "plot the rank-frequency zipfs law relationship of kmer counts"
#             action = :command
#         "gfa"
#             help = "render a sequence graph as an SVG"
#             action = :command
#     end

#     settings["plot"]["histogram"].description = ""
#     @add_arg_table settings["plot"]["histogram"] begin
#         "--histogram"
#             help = "kmer counts histogram file to plot"
#             arg_type = String
#             required = true
#     end

#     settings["plot"]["rank-frequency"].description = ""
#     @add_arg_table settings["plot"]["rank-frequency"] begin
#         "--kmer-counts"
#             help = "kmer count file to plot"
#             arg_type = String
#             required = true
#     end

#     settings["plot"]["gfa"].description = ""
#     @add_arg_table settings["plot"]["gfa"] begin
#         "--gfa"
#             help = "gfa file to plot"
#             arg_type = String
#             required = true
#     end

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
                  using the tag CL:Z:# in the output GFA file, where # is an interger representing
                  the order of the color groups as they were presented in the command line.
                  """
           arg_type = String
           nargs = '*'
    end

    return parse_args(settings)
end

function main()
    parsed_args = parse_commandline()
    command = parsed_args["%COMMAND%"]
    f_string = "Mycelia.$(command)"
    f = eval(Meta.parse(f_string))
    f(parsed_args[command])
    @show parsed_args
    @show f
end

main()