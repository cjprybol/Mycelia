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
    likelihood_rank::Int
    abundance_rank::Int
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
    layout_assembly::String
    layout_gfa::String
    likelihood_ranked_variants_fasta::String
    abundance_ranked_variants_fasta::String
    ranking_tsv::String
    corrected_fastq::String
    variants::Vector{RankedVariant}
    provenance::Dict{String, Any}
end

function _prepare_stage2_inputs(
        reads::R,
        config::AssemblyConfig
)::NamedTuple where {R}
    stage1 = _run_stage1_correction(
        reads, config; materialize_corrected_reads = false)
    graph = get(stage1.result_dict, :final_graph, nothing)
    graph_k = get(stage1.result_dict[:metadata], :final_graph_k, stage1.max_k)
    graph_k isa Int || (graph_k = stage1.max_k)
    graph_source = "stage1-final"
    if graph === nothing
        graph_source = "corrected-fastq-rebuild"
        graph = open(FASTX.FASTQ.Reader, stage1.corrected_fastq) do reader
            build_kmer_graph(collect(reader), graph_k; mode = :doublestrand)
        end
    end
    index = _transition_likelihood_index(graph, graph_k)
    return (;
        corrected_fastq = stage1.corrected_fastq,
        corrected_read_count = stage1.corrected_read_count,
        graph_source,
        graph_k,
        index)
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

function _read_gfa_records(path::String)::Dict{String, NamedTuple}
    isfile(path) || error("GFA file does not exist: $(path)")
    segments = Dict{String, NamedTuple}()
    open(path, "r") do io
        for line in eachline(io)
            fields = split(line, '\t')
            length(fields) >= 3 || continue
            fields[1] == "S" || continue
            fields[3] == "*" && continue
            coverage = nothing
            for tag in fields[4:end]
                startswith(tag, "dp:") || continue
                parsed = tryparse(Float64, split(tag, ':'; limit = 3)[3])
                parsed === nothing || (coverage = parsed)
            end
            segments[String(fields[2])] = (;
                sequence = String(fields[3]), coverage)
        end
    end
    isempty(segments) && error("GFA contains no sequence-bearing segments: $(path)")
    return segments
end

function _read_gfa_segments(path::String)::Dict{String, String}
    return Dict(id => record.sequence for (id, record) in _read_gfa_records(path))
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

function _read_support_coverages(
        records::Dict{String, NamedTuple},
        fastq::String,
        threads::Int
)::Dict{String, Float64}
    isfile(fastq) || error("support FASTQ does not exist: $(fastq)")
    mktempdir() do directory
        candidates = joinpath(directory, "candidates.fasta")
        open(candidates, "w") do io
            for (id, record) in sort(collect(records); by = first)
                println(io, ">$(id)")
                println(io, record.sequence)
            end
        end
        command = `$(Mycelia.CONDA_RUNNER) run -n strainy minimap2 -x map-ont -t $(threads) $(candidates) $(fastq)`
        paf = read(command, String)
        best = Dict{String, Tuple{Int, String, Int}}()
        for line in eachline(IOBuffer(paf))
            fields = split(line, '\t')
            length(fields) >= 12 || continue
            query = String(fields[1])
            target = String(fields[6])
            matches = something(tryparse(Int, fields[10]), 0)
            aligned = something(tryparse(Int, fields[11]), 0)
            current = get(best, query, (-1, "", 0))
            (matches, aligned, target) > (current[1], current[3], current[2]) || continue
            best[query] = (matches, target, aligned)
        end
        aligned_bases = Dict(id => 0.0 for id in keys(records))
        for (_, target, aligned) in values(best)
            aligned_bases[target] = get(aligned_bases, target, 0.0) + aligned
        end
        return Dict(id => aligned_bases[id] / length(record.sequence)
            for (id, record) in records)
    end
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
        min_scored_fraction::Float64 = 0.9,
        support_fastq::Union{Nothing, String} = nothing,
        threads::Int = Mycelia.get_default_threads()
)::Vector{RankedVariant}
    max_variants > 0 || throw(ArgumentError("max_variants must be positive"))
    0.0 <= min_scored_fraction <= 1.0 ||
        throw(ArgumentError("min_scored_fraction must be between 0 and 1"))

    segments = _read_gfa_records(strain_unitigs_gfa)
    coverages = _strainy_coverages(phased_unitig_info)
    support_coverages = support_fastq === nothing ? Dict{String, Float64}() :
                        _read_support_coverages(segments, support_fastq, threads)
    candidates = NamedTuple[]
    for (id, record) in segments
        sequence = record.sequence
        coverage = haskey(support_coverages, id) ? support_coverages[id] :
                   something(record.coverage, get(coverages, id, nothing))
        coverage === nothing && continue
        score = _score_variant_sequence(sequence, index)
        score.scored_fraction >= min_scored_fraction || continue
        isfinite(score.mean_log2_probability) || continue
        push!(candidates, (;
            id,
            sequence,
            graph_mean_log2_probability = score.mean_log2_probability,
            scored_transition_fraction = score.scored_fraction,
            strainy_coverage = Float64(coverage)
        ))
    end
    likelihood_order = sort(copy(candidates);
        by = candidate -> (-candidate.graph_mean_log2_probability,
            -candidate.strainy_coverage, -length(candidate.sequence), candidate.id))
    abundance_order = sort(copy(candidates);
        by = candidate -> (-candidate.strainy_coverage,
            -candidate.graph_mean_log2_probability, -length(candidate.sequence), candidate.id))
    length(likelihood_order) >= max_variants ||
        error("Stage-2 produced only $(length(likelihood_order)) graph-supported variants; " *
              "$(max_variants) required")
    likelihood_rank = Dict(candidate.id => rank
        for (rank, candidate) in enumerate(likelihood_order))
    abundance_rank = Dict(candidate.id => rank
        for (rank, candidate) in enumerate(abundance_order))
    selected_ids = union(
        Set(candidate.id for candidate in likelihood_order[1:max_variants]),
        Set(candidate.id for candidate in abundance_order[1:max_variants]))
    selected = filter(candidate -> candidate.id in selected_ids, likelihood_order)
    return [RankedVariant(likelihood_rank[candidate.id], abundance_rank[candidate.id],
                _variant_role(likelihood_rank[candidate.id]), candidate.id,
                candidate.sequence, candidate.graph_mean_log2_probability,
                candidate.scored_transition_fraction, candidate.strainy_coverage)
            for candidate in selected]
end

function _write_ranked_variants(
        variants::Vector{RankedVariant},
        output_dir::String,
        max_variants::Int
)::NamedTuple
    mkpath(output_dir)
    likelihood_primary = only(filter(variant -> variant.likelihood_rank == 1, variants))
    abundance_primary = only(filter(variant -> variant.abundance_rank == 1, variants))
    likelihood_primary.id == abundance_primary.id ||
        error("likelihood and abundance rankings disagree on primary: " *
              "$(likelihood_primary.id) != $(abundance_primary.id)")
    primary = joinpath(output_dir, "primary_consensus.fasta")
    likelihood_fasta = joinpath(output_dir, "ranked_variants_likelihood.fasta")
    abundance_fasta = joinpath(output_dir, "ranked_variants_abundance.fasta")
    tsv = joinpath(output_dir, "ranked_variants.tsv")
    open(likelihood_fasta, "w") do io
        for variant in filter(variant -> variant.likelihood_rank <= max_variants,
                sort(variants; by = variant -> variant.likelihood_rank))
            println(io, ">likelihood_rank=$(variant.likelihood_rank)|id=$(variant.id)")
            println(io, variant.sequence)
        end
    end
    open(abundance_fasta, "w") do io
        for variant in filter(variant -> variant.abundance_rank <= max_variants,
                sort(variants; by = variant -> variant.abundance_rank))
            println(io, ">abundance_rank=$(variant.abundance_rank)|id=$(variant.id)")
            println(io, variant.sequence)
        end
    end
    open(primary, "w") do io
        variant = likelihood_primary
        println(io, ">primary|id=$(variant.id)")
        println(io, variant.sequence)
    end
    table = DataFrames.DataFrame(
        likelihood_rank = [variant.likelihood_rank for variant in variants],
        abundance_rank = [variant.abundance_rank for variant in variants],
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
    return (; primary, likelihood_fasta, abundance_fasta, tsv)
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
        genome_size::Union{Nothing, Integer, AbstractString} = nothing,
        min_overlap::Union{Nothing, Int} = nothing,
        layout_iterations::Int = 0,
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
    layout_iterations >= 0 ||
        throw(ArgumentError("layout_iterations must be nonnegative"))

    output_dir = config.output_dir
    layout_dir = joinpath(output_dir, "metaflye")
    strainy_dir = joinpath(output_dir, "strainy")
    ranked_dir = joinpath(output_dir, "ranked")
    any(isdir, (layout_dir, strainy_dir, ranked_dir)) &&
        error("Stage-2 output directory already contains a prior run: $(output_dir)")

    prepared = _prepare_stage2_inputs(reads, config)
    (; corrected_fastq, graph_k, index) = prepared

    layout = Mycelia.run_metaflye(;
        fastq = corrected_fastq,
        outdir = layout_dir,
        genome_size = genome_size,
        read_type = "nano-corr",
        meta = true,
        min_overlap = min_overlap,
        iterations = layout_iterations,
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
        min_scored_fraction = min_scored_fraction,
        support_fastq = corrected_fastq,
        threads = threads)
    ranked = _write_ranked_variants(variants, ranked_dir, max_variants)

    provenance = Dict{String, Any}(
        "stage1_corrector" => "iterative",
        "stage1_strategy" => String(config.strategy),
        "sequencing_tech" => "nanopore",
        "layout_tool" => "metaflye",
        "layout_read_type" => "nano-corr",
        "layout_iterations" => layout_iterations,
        "deconvolution_tool" => "strainy",
        "graph_k" => graph_k,
        "graph_source" => prepared.graph_source,
        "corrected_read_count" => prepared.corrected_read_count,
        "max_variants" => max_variants,
        "min_scored_fraction" => min_scored_fraction
    )
    return Stage2DeconvolutionResult(
        ranked.primary,
        layout.assembly,
        layout.graph,
        ranked.likelihood_fasta,
        ranked.abundance_fasta,
        ranked.tsv,
        corrected_fastq,
        variants,
        provenance)
end
