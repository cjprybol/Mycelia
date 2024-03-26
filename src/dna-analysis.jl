function parse_xam(xam; filter_unmapped=false, primary_only=false, min_mapping_quality=0, min_align_length=1)
    if occursin(r"\.bam$", xam)
        MODULE = XAM.BAM
    elseif occursin(r"\.sam$", xam)
        MODULE = XAM.SAM
    else
        error("unrecognized file extension in file: $xam")
    end
    reader = open(MODULE.Reader, xam)
    header = reader.header
    record_iterator = Iterators.filter(record -> true, reader)
    if filter_unmapped
        record_iterator = Iterators.filter(record -> MODULE.ismapped(record), record_iterator)
    end
    if primary_only
        record_iterator = Iterators.filter(record -> MODULE.isprimary(record), record_iterator)
    end
    record_iterator = Iterators.filter(record -> MODULE.mappingquality(record) >= min_mapping_quality, record_iterator)
    record_iterator = Iterators.filter(record -> MODULE.alignlength(record) >= min_align_length, record_iterator)
    records = sort(collect(record_iterator), by=x->[MODULE.refname(x), MODULE.position(x)])
    # reset header to specify sorted
    header.metainfo[1] = MODULE.MetaInfo("HD", ["VN" => 1.6, "SO" => "coordinate"])
    return (;records, header)
end


# # map reads to the assembly and run qualimap QC
# bwt_index = "$(assembled_fasta).bwt"
# if !isfile(bwt_index)
#     run(`bwa index $(assembled_fasta)`)
# end

# mapped_reads_bam = "$(assembled_fasta).bwa.bam"
# if !isfile(mapped_reads_bam)
#     run(pipeline(
#         `bwa mem -R "@RG\tID:$(config["sample identifier"])\tSM:bar" -t $(Sys.CPU_THREADS) $(assembled_fasta) $(TRIMMED_FORWARD) $(TRIMMED_REVERSE)`,
#         `samtools collate -O - -`,
#         `samtools fixmate -m - -`,
#         `samtools sort`,
#         `samtools markdup - -`,
#         `samtools view -buh`,
#         mapped_reads_bam))
# end

# if !isfile("$(mapped_reads_bam).bai")
#     run(`samtools index $(mapped_reads_bam)`)
# end

# qualimap_report_pdf = "$(assembly_dir)/qualimap/report.pdf"
# qualimap_report_txt = "$(assembly_dir)/qualimap/genome_results.txt"

# if !isfile(qualimap_report_pdf) || !isfile(qualimap_report_txt)
#     run(`
#         qualimap bamqc
#         -nt $(Sys.CPU_THREADS)
#         -bam $(mapped_reads_bam)
#         -outdir $(assembly_dir)/qualimap
#         -outformat PDF:HTML
#         --output-genome-coverage $(mapped_reads_bam).genome_coverage.txt
#         `)
# end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run fastani with a query and reference list

```jldoctest
julia> 1 + 1
2
```
"""
function fastani_list(;query_list="", reference_list="", outfile="", force=false)
    Mycelia.add_bioconda_env("fastani")
    if !isfile(outfile) || force
        # run(`fastANI --ql $(query_list) --rl $(reference_list) -o $(outfile)`)
        run(
        pipeline(
            `$(Mycelia.MAMBA) run --live-stream -n fastani fastANI --ql $(query_list) --rl $(reference_list) --threads $(Sys.CPU_THREADS) -o $(outfile)`,
            stdout=outfile * "fastani.stdout.txt",
            stderr=outfile * "fastani.stderr.txt"
            )
        )
    end
end

# ./fastANI -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
function fastani_pair(;query="", reference="", outfile="", force=false)
    Mycelia.add_bioconda_env("fastani")
    if !isfile(outfile) || force
        run(
        pipeline(
            `$(Mycelia.MAMBA) run --live-stream -n fastani fastANI -q $(query) -r $(reference) -o $(outfile)`,
            stdout=outfile * "fastani.stdout.txt",
            stderr=outfile * "fastani.stderr.txt"
            )
        )
    end
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determines the contig with the greatest number of total bases mapping to it

```jldoctest
julia> 1 + 1
2
```
"""
function determine_primary_contig(qualimap_results)
    primary_contig_index = last(findmax(qualimap_results[!, "Mapped bases"]))
    primary_contig = qualimap_results[primary_contig_index, "Contig"]
    return primary_contig
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Primary contig is defined as the contig with the most bases mapped to it

In the context of picking out phage from metagenomic assemblies
the longest contig is often bacteria whereas the highest coverage contigs are often primer-dimers or other PCR amplification artifacts.

Taking the contig that has the most bases mapped to it as a product of length * depth is cherry picked as our phage
"""
function isolate_normalized_primary_contig(assembled_fasta, assembled_gfa, qualimap_report_txt, identifier, k::Int; primary_contig_fasta = "$(identifier).primary_contig.fna")
    
    qualimap_results = parse_qualimap_contig_coverage(qualimap_report_txt)
    primary_contig = determine_primary_contig(qualimap_results)

    # Find primary contig from scaffolds, then export as primary_contig.fasta
    for record in FASTX.FASTA.Reader(open(assembled_fasta))
        record_id = FASTX.identifier(record)
        if record_id == primary_contig
            primary_contig_sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)

            # If the primary contig is circular, need to trim to remove closure scar
            if Mycelia.contig_is_circular(assembled_gfa, primary_contig)
                # trim k-length from end before writing if it matches the first k of the contig
                if primary_contig_sequence[1:k] == primary_contig_sequence[end-k+1:end]
                    for i in 1:k pop!(primary_contig_sequence) end
                end
            end

        w = FASTX.FASTA.Writer(open(primary_contig_fasta, "w")) 
            write(w, FASTX.FASTA.Record(identifier, primary_contig_sequence))
            close(w)
        end
    end
    return primary_contig_fasta
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns bool indicating whether the contig is cleanly assembled

graph_file = path to assembly graph.gfa file
contig_name = name of the contig

```jldoctest
julia> 1 + 1
2
```
"""
function contig_is_cleanly_assembled(graph_file::String, contig_name::String)
    gfa_graph = Mycelia.parse_gfa(graph_file)
    if !isempty(gfa_graph.gprops[:paths])
        # probably spades
        # gfa segment identifiers have _1 appended to fasta sequence identifiers for spades
        segment_identifier = contig_name * "_1"
        segment_node_string = gfa_graph.gprops[:paths][segment_identifier]["segments"]
        nodes_in_segment = replace.(split(segment_node_string, ','), r"[^\d]" => "")
        node_indices = [gfa_graph[n, :identifier] for n in nodes_in_segment]
    else
        # megahit
        node_indices = [gfa_graph[contig_name, :identifier]]
    end
    connected_components = Graphs.connected_components(gfa_graph)
    component_of_interest = first(filter(cc -> all(n -> n in cc, node_indices), connected_components))
    if (1 <= length(component_of_interest) <= 2)
        return true
    else
        return false
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns bool indicating whether the contig is a circle

graph_file = path to assembly graph.gfa file
contig_name = name of the contig

```jldoctest
julia> 1 + 1
2
```
"""
function contig_is_circular(graph_file::String, contig_name::String)
    gfa_graph = Mycelia.parse_gfa(graph_file)
    if !isempty(gfa_graph.gprops[:paths])
        # probably spades
        # gfa segment identifiers have _1 appended to fasta sequence identifiers for spades
        segment_identifier = contig_name * "_1"
        segment_node_string = gfa_graph.gprops[:paths][segment_identifier]["segments"]
        nodes_in_segment = replace.(split(segment_node_string, ','), r"[^\d]" => "")
        node_indices = [gfa_graph[n, :identifier] for n in nodes_in_segment]
    else
        # megahit
        node_indices = [gfa_graph[contig_name, :identifier]]
    end
    connected_components = Graphs.connected_components(gfa_graph)
    component_of_interest = first(filter(cc -> all(n -> n in cc, node_indices), connected_components))
    subgraph, vertex_map = Graphs.induced_subgraph(gfa_graph, component_of_interest)
    if Graphs.is_cyclic(subgraph)
        return true
    else
        return false
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
# uses minimap
function determine_percent_identity(reference_fasta, query_fasta)
    header = [
        "Query",
        "Query length",
        "Query start",
        "Query end",
        "Query strand",
        "Target",
        "Target length",
        "Target start",
        "Target end",
        "Matches",
        "Alignment length",
        "Mapping quality",
        "Cigar",
        "CS tag"]
    
#     asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    results5 = read(`minimap2 -x asm5 --cs -cL $reference_fasta $query_fasta`)
    if !isempty(results5)
        results = results5
    else
        @warn "no hit with asm5, trying asm10"
        results10 = read(`minimap2 -x asm10 --cs -cL $reference_fasta $query_fasta`)
        if !isempty(results10)
            results = results10
        else
            @warn "no hits with asm5 or asm10, trying asm20"
            results20 = read(`minimap2 -x asm20 --cs -cL $reference_fasta $query_fasta`)
            if !isempty(results20)
                results = results20
            end
        end
    end
    if !isempty(results)
        data =  DelimitedFiles.readdlm(IOBuffer(results), '\t')
        data_columns_of_interest = [collect(1:length(header)-2)..., collect(size(data, 2)-1:size(data, 2))...]
        minimap_results = DataFrames.DataFrame(data[:, data_columns_of_interest], header)

        equivalent_matches = reduce(vcat, map(x -> collect(eachmatch(r":([0-9]+)", replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_equivalent_bases = sum(map(match -> parse(Int, first(match.captures)), equivalent_matches))

        insertion_matches = reduce(vcat, map(x -> collect(eachmatch(r"\+([a-z]+)"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_inserted_bases = sum(map(match -> length(first(match.captures)), insertion_matches))
        deletion_matches = reduce(vcat, map(x -> collect(eachmatch(r"\-([a-z]+)"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_deleted_bases = sum(map(match -> length(first(match.captures)), deletion_matches))
        substitution_matches = reduce(vcat, map(x -> collect(eachmatch(r"\*([a-z]{2})"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_substituted_bases = length(substitution_matches)
        total_variants = length(insertion_matches) + length(deletion_matches) + length(substitution_matches)
        total_variable_bases = total_inserted_bases + total_deleted_bases + total_substituted_bases

        total_alignment_length = sum(minimap_results[!, "Alignment length"])
        total_matches = sum(minimap_results[!, "Matches"])
        
        alignment_percent_identity = round(total_matches / total_alignment_length * 100, digits=2)
        size_equivalence_to_reference = round(minimap_results[1, "Query length"]/minimap_results[1, "Target length"] * 100, digits=2)
        alignment_coverage_query = round(total_alignment_length / minimap_results[1, "Query length"] * 100, digits=2)
        alignment_coverage_reference = round(total_alignment_length / minimap_results[1, "Target length"] * 100, digits=2)

        results = DataFrames.DataFrame(
            alignment_percent_identity = alignment_percent_identity,
            total_equivalent_bases = total_equivalent_bases,
            total_alignment_length = total_alignment_length,
            query_length = minimap_results[1, "Query length"],
            total_variants = total_variants,
            total_snps = total_substituted_bases,
            total_indels = length(insertion_matches) + length(deletion_matches),
            alignment_coverage_query = alignment_coverage_query,
            alignment_coverage_reference = alignment_coverage_reference,
            size_equivalence_to_reference = size_equivalence_to_reference,
        )
    else
        query_length = length(FASTX.sequence(first(FASTX.FASTA.Reader(open(query_fasta)))))
        target_length = length(FASTX.sequence(first(FASTX.FASTA.Reader(open(reference_fasta)))))
        size_equivalence_to_reference = round(query_length/target_length * 100, digits=2)

        # unable to find any matches
        results = DataFrames.DataFrame(
            alignment_percent_identity = "",
            total_equivalent_bases = "",
            total_alignment_length = "",
            query_length = query_length,
            total_variants = "",
            total_snps = "",
            total_indels = "",
            alignment_coverage_query = 0,
            alignment_coverage_reference = 0,
            size_equivalence_to_reference = size_equivalence_to_reference
        )
    end
    return results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
# https://github.com/cjprybol/Mycelia/blob/e7fe50ffe2d18406fb70e0e24ebcfa45e0937596/notebooks/exploratory/2021-08-25-k-medoids-error-cluster-detection-multi-entity-graph-aligner-test.ipynb
function analyze_kmer_spectra(;out_directory, forward_reads="", reverse_reads="", k=17, target_coverage=0, plot_size=(600,400))
    @info "counting $k-mers"
    user_provided_reads = filter(x -> !isempty(x), [forward_reads, reverse_reads])
    canonical_kmer_counts = count_canonical_kmers(Kmers.Kmer{BioSequences.DNAAlphabet{4},k}, user_provided_reads)

    @info "determining max count"
    max_count = maximum(values(canonical_kmer_counts))
    @info "max count = $max_count"

    @info "generating histogram"
    kmer_counts_histogram = sort(collect(StatsBase.countmap(values(canonical_kmer_counts))), by=x->x[1])

    X = log2.(first.(kmer_counts_histogram))
    Y = log2.(last.(kmer_counts_histogram))
    
    @info "plotting kmer spectra"
    p = StatsPlots.scatter(
        X,
        Y,
        xlabel="log2(kmer_frequency)",
        ylabel="log2(# of kmers @ frequency)",
        label="",
        size=plot_size
    )

    earliest_y_min_index = last(findmin(Y))
    lower_boundary = X[earliest_y_min_index]
    lower_boundary_source = "first minimum"

    try
        # take the first 1/denominator datapoints in the set
        # to capture the error line on the left side of the graph
        @info "fitting error curve"
        denominators = [2^i for i in 1:5]
        coeficient_matrix = zeros(length(denominators), 2)
        for (i, denominator) in enumerate(denominators)
            prefix_index = Int(floor(length(X)/denominator))
            _x = X[1:prefix_index]
            _y = Y[1:prefix_index]
            model = GLM.lm(GLM.@formula(Y ~ X), DataFrames.DataFrame(X = _x, Y = _y))
            coeficient_matrix[i, :] = GLM.coef(model)
        end
        median_intercept = Statistics.median(coeficient_matrix[:, 1])
        median_slope = Statistics.median(coeficient_matrix[:, 2])

        X_intercept = (0 - median_intercept) / median_slope

        # some libraries detect the x_intercept being AFTER the end of the data
        # in these instances detect the earliest x-minimum
        if X_intercept < lower_boundary
            lower_boundary = X_intercept
            lower_boundary_source = "detected x-intercept"
        end
    catch
        @info "unable to fit regression"
    end

    p = StatsPlots.vline!(p,
        [lower_boundary],
        label="lower boundary ($(lower_boundary_source))"
    );
    
    is_above_lower_bounds = X .>= lower_boundary
    max_Y_post_error_intercept = first(findmax(Y[is_above_lower_bounds]))
    peak_indices = findall(is_above_lower_bounds .& (Y .== max_Y_post_error_intercept))
    peak_index = Int(round(Statistics.median(peak_indices)))

    p = StatsPlots.vline!([X[peak_index]], label="inferred sample coverage)")
    if isinteractive()
        display(p)
    end
    StatsPlots.savefig(p, "$out_directory/peak-detected.png")
    StatsPlots.savefig(p, "$out_directory/peak-detected.svg")
    
    if target_coverage != 0
        detected_coverage = 2^(X[peak_index])
        downsampling_rate = round(target_coverage/detected_coverage, sigdigits=3)
        downsampling_rate = min(downsampling_rate, 1)
        @info "downsampling rate = $downsampling_rate"

        outfile = "$out_directory/downsampling-rate.txt"
        open(outfile, "w") do io
            @info "writing downsampling rate to $outfile"
            println(io, downsampling_rate)
        end
        return downsampling_rate
    end
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_dnamer_saturation(fastxs::AbstractVector{<:AbstractString}, kmer_type; kmers_to_assess=Inf, power=10, min_count = 1)
    # canonical_kmers = Set{kmer_type}()
    canonical_kmer_counts = Dict{kmer_type, Int}()
    
    @show kmer_type
    k = Kmers.ksize(Kmers.kmertype(kmer_type))
    
    max_possible_kmers = determine_max_canonical_kmers(k, DNA_ALPHABET)
    
    if kmers_to_assess == Inf
        # want to read the whole file and predict how long that will take
        # n_records = reduce(sum, map(f -> Mycelia.count_records(f), fastxs))
        kmers_to_assess = max_possible_kmers
        p = ProgressMeter.Progress(kmers_to_assess, 1)
    else
        p = ProgressMeter.Progress(kmers_to_assess, 1)
    end
    
    sampling_points = Int[0]
    i = 0
    while power^i <= kmers_to_assess
        push!(sampling_points, power^i)
        i += 1
    end
    
    unique_kmer_counts = zeros(Int, length(sampling_points))
    
    if length(sampling_points) < 3
        @info "increase the # of reads analyzed or decrease the power to acquire more data points"
        return (;sampling_points, unique_kmer_counts)
    end
    
    kmers_assessed = 0
    for fastx in fastxs
        for record in open_fastx(fastx)      
            record_sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)
            for (index, kmer) in Kmers.EveryKmer{kmer_type}(record_sequence)
                canonical_kmer = BioSequences.canonical(kmer)
                if haskey(canonical_kmer_counts, canonical_kmer)
                    canonical_kmer_counts[canonical_kmer] += 1
                else
                    canonical_kmer_counts[canonical_kmer] = 1
                end
                kmers_assessed += 1
                if (length(canonical_kmer_counts) == max_possible_kmers)                 
                    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
                    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], length(canonical_kmer_counts))
                    return (;sampling_points, unique_kmer_counts, eof = false)
                elseif kmers_assessed in sampling_points
                    i = findfirst(sampling_points .== kmers_assessed)
                    unique_kmer_counts[i] = length(filter(x -> x[2] >= min_count, canonical_kmer_counts))
                    if i == length(sampling_points)
                        return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = false)
                    end
                end
                ProgressMeter.next!(p)
            end
        end
    end
    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], [length(canonical_kmer_counts)])    
    return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = true)
end

function assess_dnamer_saturation(fastxs::AbstractVector{<:AbstractString}; power=10, outdir::Union{Missing, String}=missing, min_k=3, max_k=31, threshold=0.1, kmers_to_assess=10_000_000, plot=true)
    ks = Primes.primes(min_k, max_k)
    minimum_saturation = Inf
    midpoint = Inf
    for k in ks
        kmer_type = Kmers.kmertype(Kmers.Kmer{BioSequences.DNAAlphabet{4},k})
        sampling_points, kmer_counts, hit_eof = assess_dnamer_saturation(fastxs, kmer_type, kmers_to_assess=kmers_to_assess, power=power)
        @show sampling_points, kmer_counts, hit_eof
        observed_midpoint_index = findfirst(i -> kmer_counts[i] > last(kmer_counts)/2, 1:length(sampling_points))
        observed_midpoint = sampling_points[observed_midpoint_index]
        initial_parameters = Float64[maximum(kmer_counts), observed_midpoint]
        @time fit = LsqFit.curve_fit(calculate_v, sampling_points, kmer_counts, initial_parameters)
        max_canonical_kmers = determine_max_canonical_kmers(k, DNA_ALPHABET)
        if hit_eof
            inferred_maximum = last(kmer_counts)
        else
            inferred_maximum = max(Int(ceil(fit.param[1])), last(kmer_counts))
            if inferred_maximum > max_canonical_kmers
                inferred_maximum = max_canonical_kmers
            end
        end

        inferred_midpoint = Int(ceil(fit.param[2]))
        predicted_saturation = inferred_maximum / max_canonical_kmers
        @show k, predicted_saturation

        if plot
            scale = 300
            p = StatsPlots.scatter(
                sampling_points,
                kmer_counts,
                label="observed kmer counts",
                ylabel="# unique kmers",
                xlabel="# kmers assessed",
                title = "sequencing saturation @ k = $k",
                legend=:outertopright,
                size=(2*scale, 1*scale),
                margins=3Plots.PlotMeasures.mm
                )
            StatsPlots.hline!(p, [max_canonical_kmers], label="absolute maximum")
            StatsPlots.hline!(p, [inferred_maximum], label="inferred maximum")
            StatsPlots.vline!(p, [inferred_midpoint], label="inferred midpoint")
            # xs = vcat(sampling_points, [last(sampling_points) * 2^i for i in 1:2])
            xs = sort([sampling_points..., inferred_midpoint])
            ys = calculate_v(xs, fit.param)
            StatsPlots.plot!(
                p,
                xs,
                ys,
                label="fit trendline")
            display(p)
            if !ismissing(outdir)
                StatsPlots.savefig(p, joinpath(outdir, "$k.png"))
                StatsPlots.savefig(p, joinpath(outdir, "$k.svg"))
            end
        end
            

        if predicted_saturation < minimum_saturation
            minimum_saturation = predicted_saturation
            min_k = k
            midpoint = inferred_midpoint 
        end
        if predicted_saturation < threshold
            if !ismissing(outdir)
                chosen_k_file = joinpath(outdir, "chosen_k.txt")
                println("chosen k = $k")
                open(chosen_k_file, "w") do io
                    println(io, k)
                end
            end
            return k
        end
    end
end

function assess_dnamer_saturation(fastx::AbstractString; power=10, outdir="", min_k=3, max_k=31, threshold=0.1, kmers_to_assess=10_000_000)
    assess_dnamer_saturation([fastx], outdir=outdir, min_k=min_k, max_k=max_k, threshold=threshold, power=power, kmers_to_assess=kmers_to_assess)
end