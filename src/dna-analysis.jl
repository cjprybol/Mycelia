
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
    # kmer_type = Kmers.kmertype(Kmers.Kmer{BioSequences.DNAAlphabet{2},31})
    
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

function assess_dnamer_saturation(fastxs::AbstractVector{<:AbstractString}; power=10, outdir="", min_k=3, max_k=31, threshold=0.1)
    if isempty(outdir)
        outdir = joinpath(pwd(), "kmer-saturation")
    end
    mkpath(outdir)
    
    ks = Primes.primes(min_k, max_k)
    minimum_saturation = Inf
    midpoint = Inf
    for k in ks
        # kmer_type = Kmers.Kmers.DNAKmer{k}
        kmer_type = Kmers.kmertype(Kmers.Kmer{BioSequences.DNAAlphabet{2},k})
        kmers_to_assess = 10_000_000
        sampling_points, kmer_counts, hit_eof = assess_dnamer_saturation(fastxs, kmer_type, kmers_to_assess=kmers_to_assess, power=power)
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
        predicted_saturation = inferred_maximum / determine_max_canonical_kmers(k, DNA_ALPHABET)
        @show k, predicted_saturation

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
        StatsPlots.hline!(p, [determine_max_canonical_kmers(k, DNA_ALPHABET)], label="absolute maximum")
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

function assess_dnamer_saturation(fastx::AbstractString; power=10, outdir="", min_k=3, max_k=31, threshold=0.1)
    assess_dnamer_saturation([fastx], outdir=outdir, min_k=min_k, max_k=max_k, threshold=threshold, power=power)
end