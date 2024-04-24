function simulate_variants(fasta_record::FASTX.FASTA.Record;        
        n_variants = Int(floor(sqrt(length(FASTX.sequence(fasta_record))))),
        window_size = Int(ceil(length(FASTX.sequence(fasta_record)) / n_variants)),
        variant_size_disbribution = Distributions.Geometric(1/sqrt(window_size)),
        variant_type_likelihoods = [
            :substitution => 10^-1,
            :insertion => 10^-2,
            :deletion => 10^-2,
            :inversion => 10^-2,
            # special case insertion/deletions, skipping
            # :translocations => 10^-3,
            # :duplication => 10^-3,
        ]
    )
    vcf_table = DataFrames.DataFrame(
        "#CHROM" => String[],
        "POS" => Int[],
        "ID" => String[],
        "REF" => String[],
        "ALT" => String[],
        "QUAL" => Int[],
        "FILTER" => String[],
        "INFO" => String[],
        "FORMAT" => String[],
        "SAMPLE" => String[]
    )
    @assert join(names(vcf_table), '\t') == "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE"
    
    original_sequence = BioSequences.LongDNA{4}(FASTX.sequence(fasta_record))
    original_sequence_id = string(first(split(FASTX.description(fasta_record))))
    modified_sequence_id = "modified__" * original_sequence_id
    
    variant_sizes = rand(variant_size_disbribution, n_variants) .+ 1
    @assert all(variant_sizes .>= 1)
    any(variant_sizes .>= window_size) && @warn "variant size distribution is too large, truncating variant sizes to fit in windows"
    variant_sizes = map(x -> min(x, window_size - 1), variant_sizes)
    # display(StatsPlots.plot(variant_size_disbribution, ylabel = "probability density", xlabel = "variant size", title = "Variant Size Distribution", legend=false))
    # display(StatsPlots.histogram(variant_sizes, ylabel="# of variants", xlabel="variant size", title="Actual samples drawn", nbins=length(unique(variant_sizes)), legend=false))
    variant_types = StatsBase.sample(first.(variant_type_likelihoods), StatsBase.weights(last.(variant_type_likelihoods)), n_variants)
    window_starts = 1:window_size:length(original_sequence)
    window_ends = window_size:window_size:length(original_sequence)
    windows = zip(window_starts, window_ends)
    variant_type_sizes = zip(variant_types, variant_sizes)
        
    for ((variant_type, variant_size), (start, stop)) in collect(zip(variant_type_sizes, windows))
        selected_start = rand(start:stop-variant_size)
        @assert selected_start <= stop
        original_subsequence = original_sequence[selected_start:selected_start+variant_size-1]
        @assert length(original_subsequence) == variant_size
        if variant_type == :substitution
            ref = original_subsequence
            alt = BioSequences.randdnaseq(variant_size)
            while alt == ref
                @info "substitution collision, resampling..."
                alt = BioSequences.randdnaseq(variant_size)
            end
        elseif variant_type == :insertion
            ref = original_subsequence
            alt = original_subsequence * BioSequences.randdnaseq(variant_size)
        elseif variant_type == :deletion
            ref = original_sequence[selected_start:selected_start+variant_size]
            alt = original_sequence[selected_start:selected_start]
        elseif variant_type == :inversion
            ref = original_subsequence
            alt = BioSequences.reverse_complement(original_subsequence)
        end 
        # @show selected_start
        row = Dict(
            "#CHROM" => original_sequence_id,
            "POS" => selected_start,
            "ID" => ".",
            "REF" => string(ref),
            "ALT" => string(alt),
            "QUAL" => 60,
            "FILTER" => string(variant_type),
            "INFO" => ".",
            "FORMAT" => "GT:GQ",
            "SAMPLE" => "1:60"
        )
        push!(vcf_table, row)
        original_prefix = original_sequence[start:selected_start-1]
        if variant_type == :deletion
            original_suffix = original_sequence[selected_start+variant_size+1:stop]
        else
            original_suffix = original_sequence[selected_start+variant_size:stop]
        end
        reconstructed_window = original_prefix * ref * original_suffix
        original_window = original_sequence[start:stop]
        @assert reconstructed_window == original_window
    end
    true_variant = vcf_table[!, "REF"] .!= vcf_table[!, "ALT"]
    if !all(true_variant)
        @warn "filtering equivalent variants"
        vcf_table = vcf_table[true_variant, :]
    end
    return vcf_table
end

function simulate_variants(fasta_file::String)
    vcf_file = fasta_file * ".vcf"
    modified_fasta_file = vcf_file * ".fna"
    vcf_table = DataFrames.DataFrame()
    for fasta_record in open_fastx(fasta_file)
        append!(vcf_table, simulate_variants(fasta_record))
    end
    Mycelia.write_vcf_table(vcf_file=vcf_file, vcf_table=vcf_table, fasta_file=fasta_file)
    return Mycelia.update_fasta_with_vcf(in_fasta = fasta_file, vcf_file = vcf_file)
end
    
function write_vcf_table(;vcf_file, vcf_table, fasta_file)
    true_variant = vcf_table[!, "REF"] .!= vcf_table[!, "ALT"]
    if !all(true_variant)
        @warn "filtering equivalent variants"
        vcf_table = vcf_table[true_variant, :]
    end
    open(vcf_file, "w") do io
        VCF_HEADER = 
        """
        ##fileformat=VCFv4.3
        ##fileDate=$(Dates.today())
        ##source=simulated-variants
        ##reference=$(fasta_file)
        ##FILTER=<ID=substitution,Description="substitution variant">
        ##FILTER=<ID=insertion,Description="insertion variant">
        ##FILTER=<ID=deletion,Description="deletion variant">
        ##FILTER=<ID=inversion,Description="inversion variant">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        """
        print(io, VCF_HEADER)
        println(io, join(names(vcf_table), '\t'))
        for row in DataFrames.eachrow(vcf_table)
            println(io, join([row[col] for col in names(vcf_table)], '\t'))
        end
    end
end
    
function update_fasta_with_vcf(;in_fasta, vcf_file, out_fasta=replace(vcf_file, ".vcf" => ".normalized.vcf.fna"))
    isfile("$(vcf_file).gz") && rm("$(vcf_file).gz")
    isfile("$(vcf_file).gz.tbi") && rm("$(vcf_file).gz.tbi")
    normalized_vcf_file = replace(vcf_file, ".vcf" => ".normalized.vcf")
    run(`$(Conda.conda) run --live-stream -n htslib bgzip $(vcf_file)`)
    run(`$(Conda.conda) run --live-stream -n tabix tabix -f -p vcf $(vcf_file).gz`)
    run(pipeline(`$(Conda.conda) run --live-stream -n bcftools bcftools norm -cs --fasta-ref $(in_fasta) $(vcf_file).gz`, normalized_vcf_file))
    rm("$(vcf_file).gz")
    rm("$(vcf_file).gz.tbi")
    isfile("$(normalized_vcf_file).gz") && rm("$(normalized_vcf_file).gz")
    isfile("$(normalized_vcf_file).gz.tbi") && rm("$(normalized_vcf_file).gz.tbi")
    isfile("$(normalized_vcf_file).fna") && rm("$(normalized_vcf_file).fna")
    run(`$(Conda.conda) run --live-stream -n htslib bgzip $(normalized_vcf_file)`)
    run(`$(Conda.conda) run --live-stream -n tabix tabix -p vcf $(normalized_vcf_file).gz`)
    run(`$(Conda.conda) run --live-stream -n bcftools bcftools consensus -f $(in_fasta) $(normalized_vcf_file).gz -o $(out_fasta)`)
    return out_fasta
end

function equivalent_fasta_sequences(fasta_1, fasta_2)
    fasta_1_hashes = Set(hash(BioSequences.LongDNA{2}(FASTX.sequence(record))) for record in Mycelia.open_fastx(fasta_1))
    fasta_2_hashes = Set(hash(BioSequences.LongDNA{2}(FASTX.sequence(record))) for record in Mycelia.open_fastx(fasta_2))
    @show setdiff(fasta_1_hashes, fasta_2_hashes)
    @show setdiff(fasta_2_hashes, fasta_1_hashes)
    return fasta_1_hashes == fasta_2_hashes
end

# function normalize_vcf(;reference_fasta, vcf, normalized_vcf=)
function normalize_vcf(;reference_fasta, vcf_file)
    normalized_vcf=replace(vcf_file, Mycelia.VCF_REGEX => ".sorted.normalized.vcf")
    out_vcf = normalized_vcf * ".gz"
    if !isfile(out_vcf)
        if occursin(r"\.gz$", x)
            gzipped_vcf = vcf_file
        else
            gzipped_vcf = "$(vcf).gz"
            run(`$(Conda.conda) run --live-stream -n htslib bgzip $(vcf_file)`)
        end
        tabix_index = "$(gzipped_vcf).tbi"
        if !isfile(tabix_index)
            run(`$(Conda.conda) run --live-stream -n tabix tabix -f -p vcf $(vcf_file).gz`)
        end
        run(pipeline(`$(Conda.conda) run --live-stream -n bcftools bcftools norm -cs --fasta-ref $(in_fasta) $(vcf_file).gz`, normalized_vcf_file))
        run(`$(Conda.conda) run --live-stream -n htslib bgzip $(normalized_vcf_file)`)
        run(`$(Conda.conda) run --live-stream -n tabix tabix -p vcf $(out_vcf)`)
    end
    return out_vcf
end