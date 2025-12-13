"""
PanTools (pangenome graph toolkit) integration.

This file provides thin Julia wrappers around the PanTools CLI installed via Bioconda.
All executions go through Mycelia's Conda runner and a dedicated `pantools` environment.
"""

"""
    PantoolsDB(path::AbstractString)

Lightweight handle for a PanTools database directory.
"""
struct PantoolsDB
    path::String
end

pantools_db_path(db::PantoolsDB) = db.path
pantools_db_path(db::AbstractString) = String(db)

const PANTOOLS_CONDA_ENV = "pantools"

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Construct Java memory options for PanTools.

# Keywords
- `xms`: e.g. `"2g"` to emit `-Xms2g`
- `xmx`: e.g. `"8g"` to emit `-Xmx8g`
"""
function pantools_java_opts(; xms::Union{Nothing,AbstractString}=nothing,
                              xmx::Union{Nothing,AbstractString}=nothing)
    opts = String[]
    if !isnothing(xms)
        push!(opts, "-Xms$(xms)")
    end
    if !isnothing(xmx)
        push!(opts, "-Xmx$(xmx)")
    end
    return opts
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a `Cmd` that runs PanTools inside the dedicated conda environment.

This function does not install anything; use [`run_pantools`](@ref) to auto-install.
"""
function pantools_cmd(args::Vector{String};
                      env::AbstractString=PANTOOLS_CONDA_ENV,
                      live_stream::Bool=true)
    if live_stream
        return `$(CONDA_RUNNER) run --live-stream -n $(env) pantools $(args)`
    end
    return `$(CONDA_RUNNER) run -n $(env) pantools $(args)`
end

pantools_cmd(args...; kwargs...) = pantools_cmd(string.(collect(args)); kwargs...)

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run PanTools inside the dedicated conda environment, creating it if needed.

# Keywords
- `env`: Conda environment name (default: `"pantools"`)
- `force_env`: Recreate the environment
- `quiet_env`: Suppress conda output during installation
- `workdir`: Run in this working directory
"""
function run_pantools(args::Vector{String};
                     env::AbstractString=PANTOOLS_CONDA_ENV,
                     live_stream::Bool=true,
                     force_env::Bool=false,
                     quiet_env::Bool=false,
                     workdir::Union{Nothing,AbstractString}=nothing)
    add_bioconda_env("pantools"; force=force_env, quiet=quiet_env)
    cmd = pantools_cmd(args; env=env, live_stream=live_stream)
    if isnothing(workdir)
        run(cmd)
    else
        run(Cmd(cmd; dir=workdir))
    end
    return cmd
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a PanTools genome locations file (one FASTA path per line).
"""
function write_pantools_genome_locations_file(genomes::AbstractVector{<:AbstractString};
                                              outfile::Union{Nothing,AbstractString}=nothing)
    @assert !isempty(genomes) "No genomes provided"
    for genome in genomes
        @assert isfile(genome) "Genome FASTA does not exist: $(genome)"
    end
    if isnothing(outfile)
        outfile = tempname() * ".genome_locations.txt"
    end
    mkpath(dirname(outfile))
    open(outfile, "w") do io
        for genome in genomes
            println(io, abspath(String(genome)))
        end
    end
    return String(outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a PanTools annotation locations file (`genome_number` + `GFF`/`GFF3` path).
"""
function write_pantools_annotation_locations_file(annotation_by_genome::AbstractDict{<:Integer,<:AbstractString};
                                                  outfile::Union{Nothing,AbstractString}=nothing)
    @assert !isempty(annotation_by_genome) "No annotations provided"
    for (genome_number, gff) in annotation_by_genome
        @assert genome_number > 0 "Genome numbers must be positive, got $(genome_number)"
        @assert isfile(gff) "Annotation file does not exist: $(gff)"
    end
    if isnothing(outfile)
        outfile = tempname() * ".annotation_locations.txt"
    end
    mkpath(dirname(outfile))
    open(outfile, "w") do io
        for genome_number in sort!(collect(keys(annotation_by_genome)))
            println(io, "$(genome_number) $(abspath(String(annotation_by_genome[genome_number])))")
        end
    end
    return String(outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a PanTools "genome numbers" file (one integer genome id per line).
"""
function write_pantools_genome_numbers_file(genome_numbers::AbstractVector{<:Integer};
                                            outfile::Union{Nothing,AbstractString}=nothing)
    @assert !isempty(genome_numbers) "No genome numbers provided"
    for genome_number in genome_numbers
        @assert genome_number > 0 "Genome numbers must be positive, got $(genome_number)"
    end
    if isnothing(outfile)
        outfile = tempname() * ".genome_numbers.txt"
    end
    mkpath(dirname(outfile))
    open(outfile, "w") do io
        for genome_number in genome_numbers
            println(io, genome_number)
        end
    end
    return String(outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a PanTools regions file.

Each entry should be a tuple:
- `(genome_number, contig, start, stop)` or
- `(genome_number, contig, start, stop, strand)`
"""
function write_pantools_regions_file(regions::AbstractVector{<:Tuple};
                                     outfile::Union{Nothing,AbstractString}=nothing)
    @assert !isempty(regions) "No regions provided"
    if isnothing(outfile)
        outfile = tempname() * ".regions.txt"
    end
    mkpath(dirname(outfile))
    open(outfile, "w") do io
        for region in regions
            @assert length(region) == 4 || length(region) == 5 "Region tuples must be length 4 or 5, got $(length(region))"
            genome_number, contig, start_pos, stop_pos = region[1:4]
            @assert genome_number isa Integer && genome_number > 0 "Region genome_number must be a positive integer"
            @assert contig isa AbstractString && !isempty(contig) "Region contig must be a non-empty string"
            @assert start_pos isa Integer && start_pos > 0 "Region start must be a positive integer"
            @assert stop_pos isa Integer && stop_pos > 0 "Region stop must be a positive integer"
            if length(region) == 4
                println(io, "$(genome_number) $(contig) $(start_pos) $(stop_pos)")
            else
                strand = region[5]
                @assert strand isa AbstractString || strand isa Char "Region strand must be a Char or String"
                println(io, "$(genome_number) $(contig) $(start_pos) $(stop_pos) $(strand)")
            end
        end
    end
    return String(outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a PanTools pangenome database from genome FASTA files.

# Keywords
- `db_dir`: Output directory for the PanTools database
- `genomes`: Vector of genome FASTA paths (writes a temporary locations file)
- `genome_locations_file`: Pre-existing locations file (one FASTA path per line)
- `threads`: Number of threads (passed as `-tN`)
- `xms`, `xmx`: Java heap settings (e.g. `"2g"`, `"8g"`)
- `force`: Remove `db_dir` if it exists
- `additional_args`: Extra CLI args inserted after the subcommand
"""
function pantools_build_pangenome(; db_dir::AbstractString,
                                  genomes::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
                                  genome_locations_file::Union{Nothing,AbstractString}=nothing,
                                  threads::Union{Nothing,Integer}=get_default_threads(),
                                  xms::Union{Nothing,AbstractString}=nothing,
                                  xmx::Union{Nothing,AbstractString}="5g",
                                  force::Bool=false,
                                  additional_args::Vector{String}=String[],
                                  live_stream::Bool=true,
                                  force_env::Bool=false,
                                  quiet_env::Bool=false,
                                  workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(String(db_dir))
    mkpath(dirname(db_dir))

    if isdir(db_dir)
        if force
            rm(db_dir; recursive=true, force=true)
        else
            if isempty(readdir(db_dir))
                rm(db_dir; recursive=true, force=true)
            else
                @info "PanTools DB already exists; skipping build" db_dir
                return PantoolsDB(db_dir)
            end
        end
    end

    if !isnothing(genomes)
        genome_locations_file = write_pantools_genome_locations_file(
            genomes;
            outfile=genome_locations_file
        )
    end
    @assert !isnothing(genome_locations_file) "Provide `genomes` or `genome_locations_file`"
    @assert isfile(genome_locations_file) "Genome locations file does not exist: $(genome_locations_file)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "build_pangenome")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)
    push!(args, abspath(String(genome_locations_file)))

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    return PantoolsDB(db_dir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add genomes to an existing PanTools pangenome database.
"""
function pantools_add_genomes(db::Union{PantoolsDB,AbstractString};
                              genomes::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
                              genome_locations_file::Union{Nothing,AbstractString}=nothing,
                              threads::Union{Nothing,Integer}=get_default_threads(),
                              xms::Union{Nothing,AbstractString}=nothing,
                              xmx::Union{Nothing,AbstractString}="5g",
                              additional_args::Vector{String}=String[],
                              live_stream::Bool=true,
                              force_env::Bool=false,
                              quiet_env::Bool=false,
                              workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    if !isnothing(genomes)
        genome_locations_file = write_pantools_genome_locations_file(
            genomes;
            outfile=genome_locations_file
        )
    end
    @assert !isnothing(genome_locations_file) "Provide `genomes` or `genome_locations_file`"
    @assert isfile(genome_locations_file) "Genome locations file does not exist: $(genome_locations_file)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "add_genomes")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)
    push!(args, abspath(String(genome_locations_file)))

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    return PantoolsDB(db_dir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add genome annotations (GFF/GFF3) to an existing PanTools pangenome database.

# Keywords
- `annotation_by_genome`: `Dict(genome_number => gff_path)`; writes a locations file
- `annotation_locations_file`: Pre-existing locations file
- `connect`: Pass `--connect` to connect features to nucleotide nodes
"""
function pantools_add_annotations(db::Union{PantoolsDB,AbstractString};
                                 annotation_by_genome::Union{Nothing,AbstractDict{<:Integer,<:AbstractString}}=nothing,
                                 annotation_locations_file::Union{Nothing,AbstractString}=nothing,
                                 connect::Bool=true,
                                 threads::Union{Nothing,Integer}=get_default_threads(),
                                 xms::Union{Nothing,AbstractString}=nothing,
                                 xmx::Union{Nothing,AbstractString}="5g",
                                 additional_args::Vector{String}=String[],
                                 live_stream::Bool=true,
                                 force_env::Bool=false,
                                 quiet_env::Bool=false,
                                 workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    if !isnothing(annotation_by_genome)
        annotation_locations_file = write_pantools_annotation_locations_file(
            annotation_by_genome;
            outfile=annotation_locations_file
        )
    end
    @assert !isnothing(annotation_locations_file) "Provide `annotation_by_genome` or `annotation_locations_file`"
    @assert isfile(annotation_locations_file) "Annotation locations file does not exist: $(annotation_locations_file)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "add_annotations")
    if connect
        push!(args, "--connect")
    end
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)
    push!(args, abspath(String(annotation_locations_file)))

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    overview = joinpath(db_dir, "annotation_overview.txt")
    return (; db=PantoolsDB(db_dir), annotation_overview=isfile(overview) ? overview : nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add phenotypes to an existing PanTools pangenome database.
"""
function pantools_add_phenotypes(db::Union{PantoolsDB,AbstractString},
                                 phenotypes_file::AbstractString;
                                 xms::Union{Nothing,AbstractString}=nothing,
                                 xmx::Union{Nothing,AbstractString}="2g",
                                 additional_args::Vector{String}=String[],
                                 live_stream::Bool=true,
                                 force_env::Bool=false,
                                 quiet_env::Bool=false,
                                 workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"
    @assert isfile(phenotypes_file) "Phenotypes file does not exist: $(phenotypes_file)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "add_phenotypes")
    append!(args, additional_args)
    push!(args, db_dir)
    push!(args, abspath(String(phenotypes_file)))

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    return PantoolsDB(db_dir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Cluster proteins into homology groups.
"""
function pantools_group(db::Union{PantoolsDB,AbstractString};
                        relaxation::Integer=4,
                        threads::Union{Nothing,Integer}=get_default_threads(),
                        xms::Union{Nothing,AbstractString}=nothing,
                        xmx::Union{Nothing,AbstractString}="5g",
                        additional_args::Vector{String}=String[],
                        live_stream::Bool=true,
                        force_env::Bool=false,
                        quiet_env::Bool=false,
                        workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"
    @assert relaxation > 0 "relaxation must be positive, got $(relaxation)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "group")
    push!(args, "--relaxation=$(relaxation)")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    return PantoolsDB(db_dir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run PanTools metrics and return common output paths when present.
"""
function pantools_metrics(db::Union{PantoolsDB,AbstractString};
                          threads::Union{Nothing,Integer}=nothing,
                          xms::Union{Nothing,AbstractString}=nothing,
                          xmx::Union{Nothing,AbstractString}="2g",
                          additional_args::Vector{String}=String[],
                          live_stream::Bool=true,
                          force_env::Bool=false,
                          quiet_env::Bool=false,
                          workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "metrics")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    metrics_per_genome = joinpath(db_dir, "metrics_per_genome.csv")
    return (; db=PantoolsDB(db_dir), metrics_per_genome=isfile(metrics_per_genome) ? metrics_per_genome : nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run PanTools gene classification and return common output paths when present.
"""
function pantools_gene_classification(db::Union{PantoolsDB,AbstractString};
                                      core_threshold::Union{Nothing,Integer}=nothing,
                                      include::Union{Nothing,Integer}=nothing,
                                      exclude::Union{Nothing,Integer}=nothing,
                                      threads::Union{Nothing,Integer}=nothing,
                                      xms::Union{Nothing,AbstractString}=nothing,
                                      xmx::Union{Nothing,AbstractString}="2g",
                                      additional_args::Vector{String}=String[],
                                      live_stream::Bool=true,
                                      force_env::Bool=false,
                                      quiet_env::Bool=false,
                                      workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "gene_classification")
    if !isnothing(core_threshold)
        @assert core_threshold > 0 "core_threshold must be positive, got $(core_threshold)"
        push!(args, "--core-threshold=$(core_threshold)")
    end
    if !isnothing(include)
        @assert include > 0 "include must be positive, got $(include)"
        push!(args, "--include=$(include)")
    end
    if !isnothing(exclude)
        @assert exclude > 0 "exclude must be positive, got $(exclude)"
        push!(args, "--exclude=$(exclude)")
    end
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    overview = joinpath(db_dir, "gene_classification_overview.txt")
    additional_copies = joinpath(db_dir, "additional_copies.csv")
    return (;
        db=PantoolsDB(db_dir),
        gene_classification_overview=isfile(overview) ? overview : nothing,
        additional_copies=isfile(additional_copies) ? additional_copies : nothing
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run PanTools pangenome structure analysis.
"""
function pantools_pangenome_structure(db::Union{PantoolsDB,AbstractString};
                                      threads::Union{Nothing,Integer}=nothing,
                                      xms::Union{Nothing,AbstractString}=nothing,
                                      xmx::Union{Nothing,AbstractString}="2g",
                                      additional_args::Vector{String}=String[],
                                      live_stream::Bool=true,
                                      force_env::Bool=false,
                                      quiet_env::Bool=false,
                                      workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "pangenome_structure")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    outdir = joinpath(db_dir, "pangenome_structure")
    return (; db=PantoolsDB(db_dir), pangenome_structure_dir=isdir(outdir) ? outdir : nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Retrieve regions or whole genomes from a PanTools database.

Provide exactly one of:
- `regions`: vector of tuples as accepted by [`write_pantools_regions_file`](@ref)
- `genome_numbers`: vector of genome numbers (one per line)
- `infile`: a pre-existing file in either format
"""
function pantools_retrieve_regions(db::Union{PantoolsDB,AbstractString};
                                  regions::Union{Nothing,AbstractVector{<:Tuple}}=nothing,
                                  genome_numbers::Union{Nothing,AbstractVector{<:Integer}}=nothing,
                                  infile::Union{Nothing,AbstractString}=nothing,
                                  threads::Union{Nothing,Integer}=nothing,
                                  xms::Union{Nothing,AbstractString}=nothing,
                                  xmx::Union{Nothing,AbstractString}="2g",
                                  additional_args::Vector{String}=String[],
                                  live_stream::Bool=true,
                                  force_env::Bool=false,
                                  quiet_env::Bool=false,
                                  workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    provided = count(x -> !isnothing(x), (regions, genome_numbers, infile))
    @assert provided == 1 "Provide exactly one of `regions`, `genome_numbers`, or `infile`"

    if !isnothing(regions)
        infile = write_pantools_regions_file(regions)
    elseif !isnothing(genome_numbers)
        infile = write_pantools_genome_numbers_file(genome_numbers)
    end
    @assert !isnothing(infile)
    @assert isfile(infile) "Input file does not exist: $(infile)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "retrieve_regions")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)
    push!(args, abspath(String(infile)))

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    retrieval_dir = joinpath(db_dir, "retrieval")
    return (; db=PantoolsDB(db_dir), retrieval_dir=isdir(retrieval_dir) ? retrieval_dir : nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a PanTools functions locations file (`genome_number` + functions file path).

PanTools expects a per-genome mapping similar to the annotations locations file.
"""
function write_pantools_functions_locations_file(functions_by_genome::AbstractDict{<:Integer,<:AbstractString};
                                                 outfile::Union{Nothing,AbstractString}=nothing)
    return write_pantools_annotation_locations_file(functions_by_genome; outfile=outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add functional annotations to an existing PanTools database.
"""
function pantools_add_functions(db::Union{PantoolsDB,AbstractString};
                                functions_by_genome::Union{Nothing,AbstractDict{<:Integer,<:AbstractString}}=nothing,
                                functions_locations_file::Union{Nothing,AbstractString}=nothing,
                                threads::Union{Nothing,Integer}=get_default_threads(),
                                xms::Union{Nothing,AbstractString}=nothing,
                                xmx::Union{Nothing,AbstractString}="2g",
                                additional_args::Vector{String}=String[],
                                live_stream::Bool=true,
                                force_env::Bool=false,
                                quiet_env::Bool=false,
                                workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    if !isnothing(functions_by_genome)
        functions_locations_file = write_pantools_functions_locations_file(
            functions_by_genome;
            outfile=functions_locations_file
        )
    end
    @assert !isnothing(functions_locations_file) "Provide `functions_by_genome` or `functions_locations_file`"
    @assert isfile(functions_locations_file) "Functions locations file does not exist: $(functions_locations_file)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "add_functions")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)
    push!(args, abspath(String(functions_locations_file)))

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    return PantoolsDB(db_dir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate functional overviews from previously-added functional annotations.
"""
function pantools_function_overview(db::Union{PantoolsDB,AbstractString};
                                    threads::Union{Nothing,Integer}=nothing,
                                    xms::Union{Nothing,AbstractString}=nothing,
                                    xmx::Union{Nothing,AbstractString}="2g",
                                    additional_args::Vector{String}=String[],
                                    live_stream::Bool=true,
                                    force_env::Bool=false,
                                    quiet_env::Bool=false,
                                    workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "function_overview")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    overview = joinpath(db_dir, "function_overview_per_group.csv")
    return (; db=PantoolsDB(db_dir), function_overview_per_group=isfile(overview) ? overview : nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run `group_info` for a provided homology groups CSV file.
"""
function pantools_group_info(db::Union{PantoolsDB,AbstractString},
                             homology_groups_file::AbstractString;
                             threads::Union{Nothing,Integer}=nothing,
                             xms::Union{Nothing,AbstractString}=nothing,
                             xmx::Union{Nothing,AbstractString}="2g",
                             additional_args::Vector{String}=String[],
                             live_stream::Bool=true,
                             force_env::Bool=false,
                             quiet_env::Bool=false,
                             workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"
    @assert isfile(homology_groups_file) "Homology groups file does not exist: $(homology_groups_file)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "group_info")
    push!(args, "-H")
    push!(args, abspath(String(homology_groups_file)))
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    return PantoolsDB(db_dir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run GO enrichment for a provided homology groups CSV file.
"""
function pantools_go_enrichment(db::Union{PantoolsDB,AbstractString},
                                homology_groups_file::AbstractString;
                                threads::Union{Nothing,Integer}=nothing,
                                xms::Union{Nothing,AbstractString}=nothing,
                                xmx::Union{Nothing,AbstractString}="2g",
                                additional_args::Vector{String}=String[],
                                live_stream::Bool=true,
                                force_env::Bool=false,
                                quiet_env::Bool=false,
                                workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"
    @assert isfile(homology_groups_file) "Homology groups file does not exist: $(homology_groups_file)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "go_enrichment")
    push!(args, "-H")
    push!(args, abspath(String(homology_groups_file)))
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    go_dir = joinpath(db_dir, "go_enrichment")
    go_csv = joinpath(go_dir, "go_enrichment.csv")
    return (; db=PantoolsDB(db_dir), go_enrichment_csv=isfile(go_csv) ? go_csv : nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run functional classification.
"""
function pantools_functional_classification(db::Union{PantoolsDB,AbstractString};
                                            threads::Union{Nothing,Integer}=nothing,
                                            xms::Union{Nothing,AbstractString}=nothing,
                                            xmx::Union{Nothing,AbstractString}="2g",
                                            additional_args::Vector{String}=String[],
                                            live_stream::Bool=true,
                                            force_env::Bool=false,
                                            quiet_env::Bool=false,
                                            workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "functional_classification")
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    outdir = joinpath(db_dir, "functional_classification")
    return (; db=PantoolsDB(db_dir), functional_classification_dir=isdir(outdir) ? outdir : nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run multiple sequence alignments from homology groups.
"""
function pantools_msa(db::Union{PantoolsDB,AbstractString},
                      homology_groups_file::AbstractString;
                      threads::Union{Nothing,Integer}=nothing,
                      xms::Union{Nothing,AbstractString}=nothing,
                      xmx::Union{Nothing,AbstractString}="2g",
                      additional_args::Vector{String}=String[],
                      live_stream::Bool=true,
                      force_env::Bool=false,
                      quiet_env::Bool=false,
                      workdir::Union{Nothing,AbstractString}=nothing)
    db_dir = abspath(pantools_db_path(db))
    @assert isdir(db_dir) "PanTools DB directory not found: $(db_dir)"
    @assert isfile(homology_groups_file) "Homology groups file does not exist: $(homology_groups_file)"

    args = String[]
    append!(args, pantools_java_opts(; xms=xms, xmx=xmx))
    push!(args, "msa")
    push!(args, "-H")
    push!(args, abspath(String(homology_groups_file)))
    if !isnothing(threads)
        @assert threads > 0 "threads must be positive, got $(threads)"
        push!(args, "-t$(threads)")
    end
    append!(args, additional_args)
    push!(args, db_dir)

    run_pantools(args;
        live_stream=live_stream,
        force_env=force_env,
        quiet_env=quiet_env,
        workdir=workdir
    )

    alignments_dir = joinpath(db_dir, "alignments")
    return (; db=PantoolsDB(db_dir), alignments_dir=isdir(alignments_dir) ? alignments_dir : nothing)
end
