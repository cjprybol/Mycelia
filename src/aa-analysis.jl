function make_diamond_db(fasta_file, db_file=fasta_file)
    @time run(`diamond makedb --in $(fasta_file) -d $(db_file)`)
end

# in order to change this to be a standard blast where we don't need all pairwise hits
# just drop the parameters id, min-score, max-target-seqs
function pairwise_diamond(joint_fasta_file)
    if !isfile("$(joint_fasta_file).dmnd")
        make_diamond_db(joint_fasta_file)
    end
    n_records = count_records(joint_fasta_file)
    # max_target_seqs = Int(ceil(sqrt(n_records)))
    @show "here!"
    sensitivity = "--iterate"
    # --block-size/-b
    # https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options
    # set block size to total memory / 8
    available_gigabytes = floor(Sys.free_memory() / 1e9)
    block_size = floor(available_gigabytes / 8)
    
    @time run(`diamond blastp $(sensitivity) --block-size $(block_size) --id 0 --min-score 0 --max-target-seqs $(n_records) --unal 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore -d $(joint_fasta_file).dmnd -q $(joint_fasta_file) -o $(joint_fasta_file).dmnd.tsv`)
    # # pairwise output is all of the alignments, super helpful!
    # # @time run(`diamond blastp $(sensitivity) --id 0 --min-score 0 --max-target-seqs $(N_RECORDS) --unal 1 --outfmt 0  -d $(joint_fasta_outfile).dmnd -q $(joint_fasta_outfile) -o $(joint_fasta_outfile).diamond.pairwise.txt`)
end

function diamond_line_to_named_tuple(diamond_line)
    sline = split(line)
    values_named_tuple = (
        qseqid = sline[1],
        sseqid = sline[2],
        pident = parse(Float64, sline[3]),
        length = parse(Int, sline[4]),
        mismatch = parse(Int, sline[5]),
        gapopen = parse(Int, sline[6]),
        qlen = parse(Int, sline[7]),
        qstart = parse(Int, sline[8]),
        qend = parse(Int, sline[9]),
        slen = parse(Int, sline[10]),
        sstart = parse(Int, sline[11]),
        send = parse(Int, sline[12]),
        evalue = parse(Float64, sline[13]),
        bitscore = parse(Float64, sline[14])
        )
    return values_named_tuple
end

function read_diamond_alignments_file(diamond_file)
    column_names_to_types = [
        "qseqid" => String,
        "sseqid" => String,
        "pident" => Float64,
        "length" => Int,
        "mismatch" => Int,
        "gapopen" => Int,
        "qlen" => Int,
        "qstart" => Int,
        "qend" => Int,
        "slen" => Int,
        "sstart" => Int,
        "send" => Int,
        "evalue" => Float64,
        "bitscore" => Float64,
    ]
    types = Dict(i => t for (i, t) in enumerate(last.(column_names_to_types)))
    data, header = uCSV.read(diamond_file, header=1, delim='\t', types = types)
    # header = first.(column_names_to_types)
    @assert header == first.(column_names_to_types)
    table = DataFrames.DataFrame(data, header)
    return table
end

# function diamond_alignments_file_to_distance_matrix(diamond_file)
#     column_names_to_types = [
#         "qseqid" => String,
#         "sseqid" => String,
#         "pident" => Float64,
#         "length" => Int,
#         "mismatch" => Int,
#         "gapopen" => Int,
#         "qlen" => Int,
#         "qstart" => Int,
#         "qend" => Int,
#         "slen" => Int,
#         "sstart" => Int,
#         "send" => Int,
#         "evalue" => Float64,
#         "bitscore" => Float64,
#     ]
#     types = Dict(i => t for (i, t) in enumerate(last.(column_names_to_types)))
#     data, header = uCSV.read(diamond_file, header=1, delim='\t', types = types)
#     # header = first.(column_names_to_types)
#     @assert header == first.(column_names_to_types)
#     table = DataFrames.DataFrame(data, header)
#     return table
# end

function add_header_to_diamond_file(infile, outfile=replace(infile, ".tsv" => ".with-header.tsv"))
    column_names = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qlen",
        "qstart",
        "qend",
        "slen",
        "sstart",
        "send",
        "evalue",
        "bitscore"
    ]
    # dangerous but fast
    # try
    #     inserted_text = join(columns_names, '\t') * '\n'
    #     sed_cmd = "1s/^/$(inserted_text)/"
    #     full_cmd = `sed -i $sed_cmd $infile`
    # catch
    open(outfile, "w") do io
        println(io, join(column_names, "\t"))
        for line in eachline(infile)
            println(io, line)
        end
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return an iterator of each aamer for a given record and k length

```jldoctest
julia> 1 + 1
2
```
"""
function each_aamer(k, sequence::T) where T <: BioSequences.BioSequence
    return (sequence[i:i+k-1] for i in 1:length(sequence)-k+1)
end

function each_aamer(k, record::FASTX.FASTA.Record)
    return each_aamer(k, FASTX.sequence(record))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the aamers for a given record

```jldoctest
julia> 1 + 1
2
```
"""
function count_aamers(k, fasta_protein::FASTX.FASTA.Record)
    return sort(StatsBase.countmap(each_aamer(k, fasta_protein)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count aamers for the entire file

```jldoctest
julia> 1 + 1
2
```
"""
# function count_aamers_by_file(k, fasta_proteins::FP) where {FP <: Union{AbstractVector{FASTX.FASTA.Record}, FASTX.FASTA.Reader}}
#     aamer_counts = OrderedCollections.OrderedDict{BioSequences.LongAminoAcidSeq, Int64}()
#     for protein in fasta_proteins
#         if !(FASTX.sequence(protein) isa BioSequences.LongAminoAcidSeq)
#             # @warn "record $(protein) is not encoded as a protein sequence, skipping..."
#             continue
#         end
#         these_counts = count_aamers(k, protein)
#         merge!(+, aamer_counts, these_counts)
#     end
#     return sort(aamer_counts)
# end

function count_aamers_by_file(k, fastx)
    counts_by_record = collect(values(count_aamers_by_record(k, fastx_file)))
    counts_by_file = first(counts_by_record)
    for these_counts in counts_by_record[2:end]
        merge!(+, counts_by_file, these_counts)
    end
    return counts_by_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count aamers by record for a given file

```jldoctest
julia> 1 + 1
2
```
"""
function count_aamers_by_record(k, fastx_file)
    return Dict(FASTX.identifier(record) => count_aamers(k, record) for record in open_fastx(fastx_file))
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Count the aamers across all records of a given file

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function count_aamers_by_file(k, fastx_file)
#     kmer_counts = StatsBase.countmap(each_aamer(k, record) for record in open_fastx(fastx_file))
#     kmer_counts = sort(kmer_counts)
#     kmer_counts_table = 
#     DataFrames.DataFrame(
#         kmer = collect(keys(kmer_counts)),
#         count = collect(values(kmer_counts))
#     )
#     return kmer_counts_table
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Count the aamers of each record in a file, returning the counts as a table

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function count_aamers_by_record(k, fastx_file)
#     kmer_counts_table = 
#     DataFrames.DataFrame(
#         record_identifier = String[],
#         kmer = BioSequences.LongAminoAcidSeq[],
#         count = Int[]
#         )
#     for record in open_fastx(fastx_file)
#         kmer_counts = StatsBase.countmap(each_aamer(k, record))
#         kmer_counts = sort(kmer_counts)
#         for (kmer, count) in sort(kmer_counts)
#             row = (
#                 record_identifier = FASTX.identifier(record),
#                 kmer = kmer,
#                 count = count
#                 )
#             push!(kmer_counts_table, row)
#         end
#     end
#     return kmer_counts_table
# end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_aamer_saturation(fastxs, k; kmers_to_assess=Inf, power=10)
    kmers = Set{BioSequences.LongAASeq}()
    
    max_possible_kmers = determine_max_possible_kmers(k, AA_ALPHABET)
    
    if kmers_to_assess == Inf
        kmers_to_assess = max_possible_kmers
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
    
    p = ProgressMeter.Progress(kmers_to_assess, 1)
    
    kmers_assessed = 0
    for fastx in fastxs
        for record in open_fastx(fastx)
            
            for kmer in each_aamer(k, FASTX.sequence(record))
                push!(kmers, canonical_kmer)
                kmers_assessed += 1
                if (length(kmers) == max_possible_kmers)                 
                    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
                    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], length(kmers))
                    return (;sampling_points, unique_kmer_counts, eof = false)
                elseif kmers_assessed in sampling_points
                    i = findfirst(sampling_points .== kmers_assessed)
                    unique_kmer_counts[i] = length(kmers)
                    if i == length(sampling_points)
                        return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = false)
                    end
                end
                ProgressMeter.next!(p)
            end
        end
    end
    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], [length(kmers)])    
    return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = true)
end

function assess_aamer_saturation(fastxs; outdir="", min_k=1, max_k=15, threshold=0.1)
    
    if isempty(outdir)
        outdir = joinpath(pwd(), "kmer-saturation")
    end
    mkpath(outdir)
    
    ks = Primes.primes(min_k, max_k)
    minimum_saturation = Inf
    midpoint = Inf
    for k in ks
        kmers_to_assess = 10_000_000
        max_possible_kmers = determine_max_possible_kmers(k, AA_ALPHABET)
        sampling_points, kmer_counts, hit_eof = assess_aamer_saturation(fastxs, k, kmers_to_assess=kmers_to_assess)
        @show sampling_points, kmer_counts, hit_eof
        observed_midpoint_index = findfirst(i -> kmer_counts[i] > last(kmer_counts)/2, 1:length(sampling_points))
        observed_midpoint = sampling_points[observed_midpoint_index]
        initial_parameters = Float64[maximum(kmer_counts), observed_midpoint]
        @time fit = LsqFit.curve_fit(calculate_v, sampling_points, kmer_counts, initial_parameters)
        if hit_eof
            inferred_maximum = last(kmer_counts)
        else
            inferred_maximum = max(Int(ceil(fit.param[1])), last(kmer_counts))
        end

        inferred_midpoint = Int(ceil(fit.param[2]))
        predicted_saturation = inferred_maximum / max_possible_kmers
        @show k, predicted_saturation

        p = StatsPlots.scatter(
            sampling_points,
            kmer_counts,
            label="observed kmer counts",
            ylabel="# unique kmers",
            xlabel="# kmers assessed",
            title = "sequencing saturation @ k = $k",
            legend=:outertopright,
            size=(800, 400),
            margins=3Plots.PlotMeasures.mm
            )
        StatsPlots.hline!(p, [max_possible_kmers], label="absolute maximum")
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
        StatsPlots.savefig(p, joinpath(outdir, "$k.png"))
        StatsPlots.savefig(p, joinpath(outdir, "$k.svg"))

        if predicted_saturation < minimum_saturation
            minimum_saturation = predicted_saturation
            min_k = k
            midpoint = inferred_midpoint 
        end
        if predicted_saturation < threshold
            chosen_k_file = joinpath(outdir, "chosen_k.txt")
            println("chosen k = $k")
            open(chosen_k_file, "w") do io
                println(io, k)
            end
            return k
        end
    end
end