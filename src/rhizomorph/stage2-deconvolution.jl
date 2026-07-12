"""
Compact transition-probability index extracted from the Stage-1 accurized graph.

The index deliberately excludes read-level evidence and graph topology so the
full corrector graph can be released before the external layout/deconvolution
tools run.
"""
struct TransitionLikelihoodIndex
    k::Int
    log2_probability::Dict{Tuple{String, String}, Float64}
end

"""
One ranked Stage-2 haplotype candidate.
"""
struct RankedVariant
    rank::Int
    role::Symbol
    id::String
    sequence::String
    graph_mean_log2_probability::Float64
    scored_transition_fraction::Float64
    strainy_coverage::Float64
end

"""
Persistent outputs from the Stage-2 metaFlye + Strainy route.
"""
struct Stage2DeconvolutionResult
    primary_fasta::String
    layout_gfa::String
    ranked_variants_fasta::String
    ranking_tsv::String
    corrected_fastq::String
    variants::Vector{RankedVariant}
    provenance::Dict{String, Any}
end

function _transition_likelihood_index(
        graph::MetaGraphsNext.MetaGraph,
        k::Int
)::TransitionLikelihoodIndex
    outgoing = Dict{String, Float64}()
    weights = Dict{Tuple{String, String}, Float64}()
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        src_string = string(src)
        dst_string = string(dst)
        weight = _edge_transition_weight(graph[src, dst])
        weights[(src_string, dst_string)] = weight
        outgoing[src_string] = get(outgoing, src_string, 0.0) + weight
    end

    probabilities = Dict{Tuple{String, String}, Float64}()
    for (edge, weight) in weights
        total = get(outgoing, edge[1], 0.0)
        total > 0.0 || continue
        probabilities[edge] = log2(weight / total)
    end
    return TransitionLikelihoodIndex(k, probabilities)
end

function _score_variant_sequence(
        sequence::AbstractString,
        index::TransitionLikelihoodIndex
)::NamedTuple{(:mean_log2_probability, :scored_fraction), Tuple{Float64, Float64}}
    n_transitions = length(sequence) - index.k
    n_transitions > 0 || return (mean_log2_probability = -Inf, scored_fraction = 0.0)

    score = 0.0
    n_scored = 0
    for i in 1:n_transitions
        src = String(sequence[i:(i + index.k - 1)])
        dst = String(sequence[(i + 1):(i + index.k)])
        edge_score = get(index.log2_probability, (src, dst), nothing)
        if edge_score !== nothing
            score += edge_score
            n_scored += 1
        end
    end
    fraction = n_scored / n_transitions
    mean_score = n_scored == 0 ? -Inf : score / n_scored
    return (mean_log2_probability = mean_score, scored_fraction = fraction)
end

function _read_gfa_segments(path::String)::Dict{String, String}
    isfile(path) || error("GFA file does not exist: $(path)")
    segments = Dict{String, String}()
    open(path, "r") do io
        for line in eachline(io)
            fields = split(line, '\t')
            length(fields) >= 3 || continue
            fields[1] == "S" || continue
            fields[3] == "*" && continue
            segments[String(fields[2])] = String(fields[3])
        end
    end
    isempty(segments) && error("GFA contains no sequence-bearing segments: $(path)")
    return segments
end

function _find_column(
        table::DataFrames.DataFrame,
        candidates::Vector{String},
        purpose::String
)::String
    lookup = Dict(lowercase(String(name)) => String(name) for name in DataFrames.names(table))
    for candidate in candidates
        haskey(lookup, candidate) && return lookup[candidate]
    end
    error("Strainy phased-unitig table lacks a $(purpose) column; " *
          "columns=$(join(DataFrames.names(table), ", "))")
end

function _strainy_coverages(path::String)::Dict{String, Float64}
    isfile(path) || error("Strainy phased-unitig table does not exist: $(path)")
    table = CSV.read(path, DataFrames.DataFrame)
    id_column = _find_column(table,
        ["unitig", "unitig_id", "strain_unitig", "id", "name"], "unitig identifier")
    coverage_column = _find_column(table,
        ["coverage", "mean_coverage", "median_coverage", "cov"], "coverage")
    coverages = Dict{String, Float64}()
    for row in DataFrames.eachrow(table)
        id = String(row[id_column])
        coverage = tryparse(Float64, string(row[coverage_column]))
        coverage === nothing && error("non-numeric Strainy coverage for $(id)")
        coverages[id] = coverage
    end
    return coverages
end

function _variant_role(rank::Int)::Symbol
    rank == 1 && return :primary
    rank == 2 && return :secondary
    rank == 3 && return :tertiary
    return :alternate
end

function _rank_stage2_variants(
        strain_unitigs_gfa::String,
        phased_unitig_info::String,
        index::TransitionLikelihoodIndex;
        max_variants::Int = 3,
        min_scored_fraction::Float64 = 0.9
)::Vector{RankedVariant}
    max_variants > 0 || throw(ArgumentError("max_variants must be positive"))
    0.0 <= min_scored_fraction <= 1.0 ||
        throw(ArgumentError("min_scored_fraction must be between 0 and 1"))

    segments = _read_gfa_segments(strain_unitigs_gfa)
    coverages = _strainy_coverages(phased_unitig_info)
    candidates = NamedTuple[]
    for (id, sequence) in segments
        haskey(coverages, id) || continue
        score = _score_variant_sequence(sequence, index)
        score.scored_fraction >= min_scored_fraction || continue
        isfinite(score.mean_log2_probability) || continue
        push!(candidates, (;
            id,
            sequence,
            graph_mean_log2_probability = score.mean_log2_probability,
            scored_transition_fraction = score.scored_fraction,
            strainy_coverage = coverages[id]
        ))
    end
    sort!(candidates;
        by = candidate -> (-candidate.graph_mean_log2_probability,
            -candidate.strainy_coverage, -length(candidate.sequence), candidate.id))
    length(candidates) >= max_variants ||
        error("Stage-2 produced only $(length(candidates)) graph-supported variants; " *
              "$(max_variants) required")

    return [RankedVariant(rank, _variant_role(rank), candidate.id,
                candidate.sequence, candidate.graph_mean_log2_probability,
                candidate.scored_transition_fraction, candidate.strainy_coverage)
            for (rank, candidate) in enumerate(candidates[1:max_variants])]
end

function _write_ranked_variants(
        variants::Vector{RankedVariant},
        output_dir::String
)::NamedTuple{(:primary, :fasta, :tsv), Tuple{String, String, String}}
    mkpath(output_dir)
    primary = joinpath(output_dir, "primary_consensus.fasta")
    fasta = joinpath(output_dir, "ranked_variants.fasta")
    tsv = joinpath(output_dir, "ranked_variants.tsv")
    open(fasta, "w") do io
        for variant in variants
            println(io, ">rank=$(variant.rank)|role=$(variant.role)|id=$(variant.id)")
            println(io, variant.sequence)
        end
    end
    open(primary, "w") do io
        variant = first(variants)
        println(io, ">primary|id=$(variant.id)")
        println(io, variant.sequence)
    end
    table = DataFrames.DataFrame(
        rank = [variant.rank for variant in variants],
        role = [String(variant.role) for variant in variants],
        id = [variant.id for variant in variants],
        length = [length(variant.sequence) for variant in variants],
        graph_mean_log2_probability =
            [variant.graph_mean_log2_probability for variant in variants],
        scored_transition_fraction =
            [variant.scored_transition_fraction for variant in variants],
        strainy_coverage = [variant.strainy_coverage for variant in variants]
    )
    CSV.write(tsv, table; delim = '\t')
    return (; primary, fasta, tsv)
end

"""
    deconvolve_stage2(reads, config::AssemblyConfig; kwargs...)

Run Stage-1 graph accurization once, reuse metaFlye for long-read repeat-graph
layout, and reuse Strainy for strain-haplotype deconvolution. Candidate
haplotypes are ranked by normalized likelihood under the accurized graph.
"""
function deconvolve_stage2(
        reads::R,
        config::AssemblyConfig;
        genome_size::Union{Nothing, String} = nothing,
        min_overlap::Union{Nothing, Int} = nothing,
        max_variants::Int = 3,
        min_scored_fraction::Float64 = 0.9,
        threads::Int = Mycelia.get_default_threads()
)::Stage2DeconvolutionResult where {R}
    config.corrector == :iterative ||
        throw(ArgumentError("Stage-2 requires corrector=:iterative"))
    config.sequencing_tech == :nanopore ||
        throw(ArgumentError("the first Stage-2 slice requires sequencing_tech=:nanopore"))
    config.output_dir === nothing &&
        throw(ArgumentError("Stage-2 requires a persistent AssemblyConfig.output_dir"))

    output_dir = config.output_dir
    layout_dir = joinpath(output_dir, "metaflye")
    strainy_dir = joinpath(output_dir, "strainy")
    ranked_dir = joinpath(output_dir, "ranked")
    any(isdir, (layout_dir, strainy_dir, ranked_dir)) &&
        error("Stage-2 output directory already contains a prior run: $(output_dir)")

    stage1 = _run_stage1_correction(reads, config)
    corrected_fastq = stage1.corrected_fastq
    graph = get(stage1.result_dict, :final_graph, nothing)
    graph === nothing && error("Stage-1 did not return a final accurized graph")
    graph_k = get(stage1.result_dict[:metadata], :final_graph_k, nothing)
    graph_k isa Int || error("Stage-1 metadata lacks an integer :final_graph_k")
    index = _transition_likelihood_index(graph, graph_k)

    # The external OLC/deconvolution tail needs only corrected reads plus the
    # compact likelihood index. Drop the full evidence graph before launching it.
    stage1.result_dict[:final_graph] = nothing
    graph = nothing

    layout = Mycelia.run_metaflye(;
        fastq = corrected_fastq,
        outdir = layout_dir,
        genome_size = genome_size,
        read_type = "nano-corr",
        meta = true,
        min_overlap = min_overlap,
        iterations = 0,
        keep_haplotypes = true,
        no_alt_contigs = true,
        threads = threads)
    isfile(layout.assembly) || error("metaFlye did not produce assembly.fasta")
    isfile(layout.graph) || error("metaFlye did not produce assembly_graph.gfa")

    phased = Mycelia.run_strainy(;
        gfa_ref = layout.graph,
        fastq = corrected_fastq,
        outdir = strainy_dir,
        read_mode = :nano,
        stage = :e2e,
        min_unitig_coverage = 3,
        unitig_split_length = 0,
        threads = threads)
    variants = _rank_stage2_variants(
        phased.strain_contigs_gfa, phased.phased_unitig_info, index;
        max_variants = max_variants,
        min_scored_fraction = min_scored_fraction)
    ranked = _write_ranked_variants(variants, ranked_dir)

    provenance = Dict{String, Any}(
        "stage1_corrector" => "iterative",
        "stage1_strategy" => String(config.strategy),
        "sequencing_tech" => "nanopore",
        "layout_tool" => "metaflye",
        "layout_read_type" => "nano-corr",
        "deconvolution_tool" => "strainy",
        "graph_k" => graph_k,
        "max_variants" => max_variants,
        "min_scored_fraction" => min_scored_fraction
    )
    return Stage2DeconvolutionResult(
        ranked.primary,
        layout.graph,
        ranked.fasta,
        ranked.tsv,
        corrected_fastq,
        variants,
        provenance)
end
