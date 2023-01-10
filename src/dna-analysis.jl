# https://github.com/cjprybol/Mycelia/blob/e7fe50ffe2d18406fb70e0e24ebcfa45e0937596/notebooks/exploratory/2021-08-25-k-medoids-error-cluster-detection-multi-entity-graph-aligner-test.ipynb

# TODO, replace jellyfish functions with internal code
function analyze_kmer_spectra(;kmer_directory, forward_reads, reverse_reads, k=17, target_coverage=0)
    @info "counting $k-mers"
    jf_file = "$kmer_directory/$k.downsampled.jf"
    p = pipeline(
        `gzip -dc $forward_reads $reverse_reads`,
        `jellyfish count -m $k -s 100M -t $(Sys.CPU_THREADS) -C --output $jf_file /dev/fd/0`
    )
    run(p)

    @info "determining max count"
    p = pipeline(
        `jellyfish dump -ct $jf_file`,
        `awk '{print $2}'`,
        `sort --numeric-sort --reverse`,
        `head -n1`)
    max_count = parse(Int, read(p, String))
    @info "max count = $max_count"

    @info "generating histogram"
    histogram_matrix = Int.(DelimitedFiles.readdlm(open(`jellyfish histo --high $max_count $jf_file`)))

    X = log2.(histogram_matrix[:, 1])
    Y = log2.(histogram_matrix[:, 2])
    
    @info "plotting kmer spectra"
    p = StatsPlots.scatter(
        X,
        Y,
        xlabel="log2(kmer_frequency)",
        ylabel="log2(# of kmers @ frequency)",
        label=""
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
    StatsPlots.savefig(p, "$jf_file.peak-detected.png")
    StatsPlots.savefig(p, "$jf_file.peak-detected.svg")
    
    if target_coverage != 0
        detected_coverage = 2^(X[peak_index])
        downsampling_rate = round(target_coverage/detected_coverage, sigdigits=3)
        downsampling_rate = min(downsampling_rate, 1)
        @info "downsampling rate = $downsampling_rate"

        outfile = "$kmer_directory/downsampling-rate.txt"
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

function assess_dnamer_saturation(fastxs::AbstractVector{<:AbstractString}; power=10, outdir="", min_k=3, max_k=31, threshold=0.1, kmers_to_assess=10_000_000)
    if isempty(outdir)
        outdir = joinpath(pwd(), "kmer-saturation")
    end
    mkpath(outdir)
    
    ks = Primes.primes(min_k, max_k)
    minimum_saturation = Inf
    midpoint = Inf
    for k in ks
        # kmer_type = Kmers.DNAKmer{k}
        kmer_type = Kmers.kmertype(Kmers.Kmer{BioSequences.DNAAlphabet{2},k})
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

function assess_dnamer_saturation(fastx::AbstractString; power=10, outdir="", min_k=3, max_k=31, threshold=0.1, kmers_to_assess=10_000_000)
    assess_dnamer_saturation([fastx], outdir=outdir, min_k=min_k, max_k=max_k, threshold=threshold, power=power, kmers_to_assess=kmers_to_assess)
end