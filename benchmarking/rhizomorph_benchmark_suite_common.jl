import Mycelia
import FASTX
import BioSequences
import Kmers
import Random
import StableRNGs
import DataFrames
import Statistics

const ISOLATE_Q_VALUES = [10, 20, 30, 40, 50, 60]
const COVERAGE_LEVELS = [10, 25, 50, 100, 250, 500, 1000]
const ISOLATE_ANCHOR_Q_VALUES = [20, 40, 60]
const ISOLATE_ANCHOR_COVERAGE_LEVELS = [25, 100, 500]
const METAGENOME_COMPLEXITY_LEVELS = [10, 25, 50]
const METAGENOME_ANCHOR_DEPTH_LEVELS = [25, 100, 500]
const METAGENOME_PRIMARY_QUALITY_PROFILES = [
    :const_q10, :const_q20, :const_q30, :const_q40, :const_q50, :const_q60
]
const METAGENOME_ANCHOR_QUALITY_PROFILES = [
    :const_q20, :const_q40, :const_q60, :mixed_balanced, :mixed_extremes
]

"""
Return quality profile definitions as `(fraction, q_value)` vectors.
"""
function metagenome_quality_profiles()
    return Dict{Symbol, Vector{Tuple{Float64, Int}}}(
        :const_q10 => [(1.0, 10)],
        :const_q20 => [(1.0, 20)],
        :const_q30 => [(1.0, 30)],
        :const_q40 => [(1.0, 40)],
        :const_q50 => [(1.0, 50)],
        :const_q60 => [(1.0, 60)],
        :mixed_balanced => [(0.5, 20), (0.3, 30), (0.2, 40)],
        :mixed_extremes => [(0.4, 10), (0.2, 30), (0.4, 60)]
    )
end

"""
Default isolate panel from viral to mid-size bacterial genomes.
"""
function default_isolate_genome_specs(; scale::String = "small")
    if scale == "large"
        return [
            (genome_id = "viral_5kb", genome_tier = "viral", genome_length = 5386, seed = 11),
            (genome_id = "bacterial_small_1mb", genome_tier = "bacterial_small",
                genome_length = 1_000_000, seed = 22),
            (genome_id = "bacterial_mid_5mb", genome_tier = "bacterial_mid",
                genome_length = 5_000_000, seed = 33)
        ]
    elseif scale == "medium"
        return [
            (genome_id = "viral_5kb", genome_tier = "viral", genome_length = 5386, seed = 11),
            (genome_id = "bacterial_small_750kb", genome_tier = "bacterial_small",
                genome_length = 750_000, seed = 22),
            (genome_id = "bacterial_mid_3mb", genome_tier = "bacterial_mid",
                genome_length = 3_000_000, seed = 33)
        ]
    end
    return [
        (genome_id = "viral_5kb", genome_tier = "viral", genome_length = 5386, seed = 11),
        (genome_id = "bacterial_small_250kb", genome_tier = "bacterial_small",
            genome_length = 250_000, seed = 22),
        (genome_id = "bacterial_mid_1200kb", genome_tier = "bacterial_mid",
            genome_length = 1_200_000, seed = 33)
    ]
end

function isolate_algorithm_specs(; include_advanced::Bool = false)
    specs = [
        # Full matrix core
        (algorithm_id = "kmer_full_fastq", graph_family = :kmer,
            memory_profile = :full, tokenizer = :none, input_mode = :fastq,
            stage = :baseline, anchor_only = false),
        (algorithm_id = "kmer_lightweight_fastq", graph_family = :kmer,
            memory_profile = :lightweight, tokenizer = :none, input_mode = :fastq,
            stage = :baseline, anchor_only = false),
        (algorithm_id = "kmer_ultralight_fastq", graph_family = :kmer,
            memory_profile = :ultralight, tokenizer = :none, input_mode = :fastq,
            stage = :baseline, anchor_only = false),

        # Anchor cross-family checks
        (algorithm_id = "kmer_full_fasta", graph_family = :kmer,
            memory_profile = :full, tokenizer = :none, input_mode = :fasta,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "kmer_lightweight_fasta", graph_family = :kmer,
            memory_profile = :lightweight, tokenizer = :none, input_mode = :fasta,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "kmer_ultralight_fasta", graph_family = :kmer,
            memory_profile = :ultralight, tokenizer = :none, input_mode = :fasta,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "ngram_full_string", graph_family = :ngram,
            memory_profile = :full, tokenizer = :none, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "ngram_lightweight_string", graph_family = :ngram,
            memory_profile = :lightweight, tokenizer = :none, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "ngram_ultralight_string", graph_family = :ngram,
            memory_profile = :ultralight, tokenizer = :none, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "olc_full_fastq", graph_family = :olc,
            memory_profile = :full, tokenizer = :none, input_mode = :fastq,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "olc_full_fasta", graph_family = :olc,
            memory_profile = :full, tokenizer = :none, input_mode = :fasta,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "olc_full_string", graph_family = :olc,
            memory_profile = :full, tokenizer = :none, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "token_bpe_string", graph_family = :token,
            memory_profile = :full, tokenizer = :bpe, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "token_unigram_string", graph_family = :token,
            memory_profile = :full, tokenizer = :unigram, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "token_char_string", graph_family = :token,
            memory_profile = :full, tokenizer = :char_fallback, input_mode = :string,
            stage = :baseline, anchor_only = true)
    ]

    if include_advanced
        push!(specs,
            (algorithm_id = "kmer_full_fastq_basic_correction", graph_family = :kmer,
                memory_profile = :full, tokenizer = :none, input_mode = :fastq,
                stage = :basic_correction, anchor_only = true))
        push!(specs,
            (algorithm_id = "kmer_full_fastq_advanced_probabilistic", graph_family = :kmer,
                memory_profile = :full, tokenizer = :none, input_mode = :fastq,
                stage = :advanced_probabilistic, anchor_only = true))
    end

    return specs
end

function metagenome_algorithm_specs(; include_advanced::Bool = false)
    specs = [
        # Full matrix core
        (algorithm_id = "kmer_full_fastq", graph_family = :kmer,
            memory_profile = :full, tokenizer = :none, input_mode = :fastq,
            stage = :baseline, anchor_only = false),
        (algorithm_id = "kmer_lightweight_fastq", graph_family = :kmer,
            memory_profile = :lightweight, tokenizer = :none, input_mode = :fastq,
            stage = :baseline, anchor_only = false),
        (algorithm_id = "kmer_ultralight_fastq", graph_family = :kmer,
            memory_profile = :ultralight, tokenizer = :none, input_mode = :fastq,
            stage = :baseline, anchor_only = false),

        # Anchor cross-family checks
        (algorithm_id = "kmer_full_fasta", graph_family = :kmer,
            memory_profile = :full, tokenizer = :none, input_mode = :fasta,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "kmer_lightweight_fasta", graph_family = :kmer,
            memory_profile = :lightweight, tokenizer = :none, input_mode = :fasta,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "kmer_ultralight_fasta", graph_family = :kmer,
            memory_profile = :ultralight, tokenizer = :none, input_mode = :fasta,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "olc_full_fastq", graph_family = :olc,
            memory_profile = :full, tokenizer = :none, input_mode = :fastq,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "olc_full_fasta", graph_family = :olc,
            memory_profile = :full, tokenizer = :none, input_mode = :fasta,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "olc_full_string", graph_family = :olc,
            memory_profile = :full, tokenizer = :none, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "ngram_full_string", graph_family = :ngram,
            memory_profile = :full, tokenizer = :none, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "token_bpe_string", graph_family = :token,
            memory_profile = :full, tokenizer = :bpe, input_mode = :string,
            stage = :baseline, anchor_only = true),
        (algorithm_id = "token_unigram_string", graph_family = :token,
            memory_profile = :full, tokenizer = :unigram, input_mode = :string,
            stage = :baseline, anchor_only = true)
    ]

    if include_advanced
        push!(specs,
            (algorithm_id = "kmer_full_fastq_basic_correction", graph_family = :kmer,
                memory_profile = :full, tokenizer = :none, input_mode = :fastq,
                stage = :basic_correction, anchor_only = true))
        push!(specs,
            (algorithm_id = "kmer_full_fastq_advanced_probabilistic", graph_family = :kmer,
                memory_profile = :full, tokenizer = :none, input_mode = :fastq,
                stage = :advanced_probabilistic, anchor_only = true))
    end

    return specs
end

function build_isolate_benchmark_matrix(;
        genomes = default_isolate_genome_specs(),
        q_values::Vector{Int} = ISOLATE_Q_VALUES,
        coverage_levels::Vector{Int} = COVERAGE_LEVELS,
        replicates::Int = 3,
        strategy::Symbol = :tiered,
        include_advanced::Bool = false,
        anchor_q_values::Vector{Int} = ISOLATE_ANCHOR_Q_VALUES,
        anchor_coverage_levels::Vector{Int} = ISOLATE_ANCHOR_COVERAGE_LEVELS)
    if !(strategy in (:tiered, :full_factorial))
        error("Unknown isolate matrix strategy: $(strategy)")
    end

    specs = isolate_algorithm_specs(include_advanced = include_advanced)
    matrix = NamedTuple[]

    for genome in genomes
        for replicate in 1:replicates
            for spec in specs
                local_q_values = strategy == :full_factorial ? q_values :
                                 (spec.anchor_only ? anchor_q_values : q_values)
                local_coverages = strategy == :full_factorial ? coverage_levels :
                                  (spec.anchor_only ? anchor_coverage_levels : coverage_levels)
                for q in local_q_values
                    for coverage in local_coverages
                        push!(matrix, (
                            dataset_type = :isolate,
                            genome_id = genome.genome_id,
                            genome_tier = genome.genome_tier,
                            genome_length = genome.genome_length,
                            q_value = q,
                            coverage = coverage,
                            replicate = replicate,
                            algorithm_id = spec.algorithm_id,
                            graph_family = spec.graph_family,
                            memory_profile = spec.memory_profile,
                            tokenizer = spec.tokenizer,
                            input_mode = spec.input_mode,
                            stage = spec.stage
                        ))
                    end
                end
            end
        end
    end

    return matrix
end

function build_metagenome_benchmark_matrix(;
        complexity_levels::Vector{Int} = METAGENOME_COMPLEXITY_LEVELS,
        depth_levels::Vector{Int} = COVERAGE_LEVELS,
        abundance_profiles::Vector{Symbol} = [:equal, :log_normal],
        quality_profiles::Vector{Symbol} = METAGENOME_PRIMARY_QUALITY_PROFILES,
        replicates::Int = 3,
        strategy::Symbol = :tiered,
        include_advanced::Bool = false,
        anchor_depth_levels::Vector{Int} = METAGENOME_ANCHOR_DEPTH_LEVELS,
        anchor_quality_profiles::Vector{Symbol} = METAGENOME_ANCHOR_QUALITY_PROFILES)
    if !(strategy in (:tiered, :full_factorial))
        error("Unknown metagenome matrix strategy: $(strategy)")
    end

    valid_quality_profiles = Set(keys(metagenome_quality_profiles()))
    for profile in quality_profiles
        if !(profile in valid_quality_profiles)
            error("Unknown metagenome quality profile: $(profile)")
        end
    end

    specs = metagenome_algorithm_specs(include_advanced = include_advanced)
    matrix = NamedTuple[]

    for complexity in complexity_levels
        for abundance_profile in abundance_profiles
            for replicate in 1:replicates
                for spec in specs
                    local_depths = strategy == :full_factorial ? depth_levels :
                                   (spec.anchor_only ? anchor_depth_levels : depth_levels)
                    local_qualities = strategy == :full_factorial ? quality_profiles :
                                      (spec.anchor_only ? anchor_quality_profiles : quality_profiles)
                    for depth_target in local_depths
                        for quality_profile in local_qualities
                            push!(matrix, (
                                dataset_type = :metagenome,
                                community_size = complexity,
                                depth_target = depth_target,
                                abundance_profile = abundance_profile,
                                quality_profile = quality_profile,
                                replicate = replicate,
                                algorithm_id = spec.algorithm_id,
                                graph_family = spec.graph_family,
                                memory_profile = spec.memory_profile,
                                tokenizer = spec.tokenizer,
                                input_mode = spec.input_mode,
                                stage = spec.stage
                            ))
                        end
                    end
                end
            end
        end
    end

    return matrix
end

function synthesize_isolate_genomes(genome_specs)
    genomes = Dict{String, BioSequences.LongDNA{4}}()
    for spec in genome_specs
        rng = StableRNGs.StableRNG(spec.seed)
        genomes[spec.genome_id] = BioSequences.randdnaseq(rng, spec.genome_length)
    end
    return genomes
end

function _quality_string(length_value::Int, q_value::Int)
    q_clamped = clamp(q_value, 0, 60)
    return repeat(string(Char(q_clamped + 33)), length_value)
end

function _record_cap_indices(total_records::Int, max_records::Union{Nothing, Int})
    if isnothing(max_records) || total_records <= max_records
        return collect(1:total_records)
    end
    return unique(round.(Int, range(1, total_records, length = max_records)))
end

function generate_isolate_read_bundle(;
        reference_sequence::BioSequences.LongDNA{4},
        coverage::Real,
        q_value::Int,
        seed::Int,
        read_length::Int = 150,
        insert_size::Int = 300,
        max_records::Union{Nothing, Int} = nothing)
    Random.seed!(seed)
    error_rate = Mycelia.q_value_to_error_rate(q_value)
    reads_1, reads_2 = Mycelia.generate_paired_end_reads(
        reference_sequence, coverage, read_length, insert_size; error_rate = error_rate)

    fastq_records = FASTX.FASTQ.Record[]
    fasta_records = FASTX.FASTA.Record[]

    for (idx, read) in enumerate(reads_1)
        seq_string = String(read)
        quality = _quality_string(length(seq_string), q_value)
        push!(fastq_records, FASTX.FASTQ.Record("read_$(idx)_R1", seq_string, quality))
        push!(fasta_records, FASTX.FASTA.Record("read_$(idx)_R1", seq_string))
    end
    for (idx, read) in enumerate(reads_2)
        seq_string = String(read)
        quality = _quality_string(length(seq_string), q_value)
        push!(fastq_records, FASTX.FASTQ.Record("read_$(idx)_R2", seq_string, quality))
        push!(fasta_records, FASTX.FASTA.Record("read_$(idx)_R2", seq_string))
    end

    selected_indices = _record_cap_indices(length(fastq_records), max_records)
    fastq_records = fastq_records[selected_indices]
    fasta_records = fasta_records[selected_indices]
    read_strings = [String(FASTX.sequence(record)) for record in fastq_records]
    return (
        fastq = fastq_records,
        fasta = fasta_records,
        strings = read_strings,
        n_reads = length(fastq_records),
        error_rate = error_rate
    )
end

function metagenome_catalog_settings(scale::String)
    if scale == "large"
        return (n_references = 400, min_length = 8000, max_length = 80_000)
    elseif scale == "medium"
        return (n_references = 180, min_length = 6000, max_length = 40_000)
    end
    return (n_references = 100, min_length = 4000, max_length = 20_000)
end

function build_synthetic_metagenome_catalog(;
        n_references::Int,
        min_length::Int = 4000,
        max_length::Int = 20_000,
        seed::Int = 777)
    rng = StableRNGs.StableRNG(seed)
    sequences = Dict{String, BioSequences.LongDNA{4}}()
    rows = NamedTuple[]

    for idx in 1:n_references
        sequence_id = "ref_$(lpad(string(idx), 4, '0'))"
        sequence_length = Random.rand(rng, min_length:max_length)
        sequence = BioSequences.randdnaseq(rng, sequence_length)
        sequences[sequence_id] = sequence
        push!(rows, (
            sequence_id = sequence_id,
            length = sequence_length,
            species = "species_$(idx)",
            genus = "genus_$(Int(ceil(idx / 3)))",
            family = "family_$(Int(ceil(idx / 10)))"
        ))
    end

    return (sequences = sequences, table = DataFrames.DataFrame(rows))
end

function _sample_quality_values(profile::Vector{Tuple{Float64, Int}},
        n::Int, rng::Random.AbstractRNG)
    cumulative = cumsum(first.(profile))
    q_values = Int[]
    for _ in 1:n
        sample_value = Random.rand(rng)
        index = findfirst(x -> sample_value <= x, cumulative)
        if isnothing(index)
            index = length(cumulative)
        end
        push!(q_values, profile[index][2])
    end
    return q_values
end

function generate_metagenome_read_bundle(;
        catalog_sequences::Dict{String, BioSequences.LongDNA{4}},
        community_size::Int,
        depth_target::Real,
        abundance_profile::Symbol,
        quality_profile::Symbol,
        seed::Int,
        read_length::Int = 150,
        insert_size::Int = 300,
        max_records_per_genome::Union{Nothing, Int} = nothing,
        max_total_records::Union{Nothing, Int} = nothing)
    rng = StableRNGs.StableRNG(seed)
    all_ids = sort(collect(keys(catalog_sequences)))
    if community_size > length(all_ids)
        error("community_size=$(community_size) exceeds catalog size=$(length(all_ids))")
    end

    shuffled_ids = Random.shuffle(rng, all_ids)
    selected_ids = shuffled_ids[1:community_size]
    weights = Mycelia.sample_abundance_weights(
        n_organisms = community_size,
        balance = abundance_profile,
        rng = rng
    )
    coverages = depth_target * community_size .* weights

    profile_definitions = metagenome_quality_profiles()
    if !haskey(profile_definitions, quality_profile)
        error("Unknown quality profile: $(quality_profile)")
    end
    q_values = _sample_quality_values(profile_definitions[quality_profile], community_size, rng)

    pooled_fastq = FASTX.FASTQ.Record[]
    pooled_fasta = FASTX.FASTA.Record[]
    pooled_strings = String[]
    truth_sequences = String[]
    truth_rows = NamedTuple[]

    for idx in eachindex(selected_ids)
        sequence_id = selected_ids[idx]
        sequence = catalog_sequences[sequence_id]
        this_coverage = coverages[idx]
        this_q = q_values[idx]
        sequence_seed = Random.rand(rng, 1:(2 ^ 31 - 1))
        bundle = generate_isolate_read_bundle(
            reference_sequence = sequence,
            coverage = this_coverage,
            q_value = this_q,
            seed = sequence_seed,
            read_length = read_length,
            insert_size = insert_size,
            max_records = max_records_per_genome
        )

        for record in bundle.fastq
            new_id = "$(sequence_id)__$(FASTX.identifier(record))"
            push!(pooled_fastq, FASTX.FASTQ.Record(
                new_id, FASTX.sequence(record), FASTX.quality(record)))
        end
        for record in bundle.fasta
            new_id = "$(sequence_id)__$(FASTX.identifier(record))"
            push!(pooled_fasta, FASTX.FASTA.Record(new_id, FASTX.sequence(record)))
        end
        append!(pooled_strings, bundle.strings)
        push!(truth_sequences, String(sequence))
        push!(truth_rows, (
            sequence_id = sequence_id,
            coverage = this_coverage,
            q_value = this_q,
            read_count = bundle.n_reads
        ))
    end

    selected_indices = _record_cap_indices(length(pooled_fastq), max_total_records)
    pooled_fastq = pooled_fastq[selected_indices]
    pooled_fasta = pooled_fasta[selected_indices]
    pooled_strings = pooled_strings[selected_indices]

    return (
        fastq = pooled_fastq,
        fasta = pooled_fasta,
        strings = pooled_strings,
        truth_sequences = truth_sequences,
        truth_table = DataFrames.DataFrame(truth_rows),
        selected_ids = selected_ids
    )
end

function select_benchmark_input(read_bundle, input_mode::Symbol)
    if input_mode == :fastq
        return read_bundle.fastq
    elseif input_mode == :fasta
        return read_bundle.fasta
    elseif input_mode == :string
        return read_bundle.strings
    end
    error("Unsupported input_mode: $(input_mode)")
end

function n50_l50(lengths::Vector{Int})
    if isempty(lengths)
        return (n50 = 0, l50 = 0)
    end
    sorted_lengths = sort(lengths; rev = true)
    total_length = sum(sorted_lengths)
    half_length = total_length / 2
    cumulative = 0
    for (index, length_value) in enumerate(sorted_lengths)
        cumulative += length_value
        if cumulative >= half_length
            return (n50 = length_value, l50 = index)
        end
    end
    return (n50 = sorted_lengths[end], l50 = length(sorted_lengths))
end

function compute_kmer_presence_metrics(reference_sequences::Vector{String},
        assembled_contigs::Vector{String}; k::Int = 31)
    reference_records = FASTX.FASTA.Record[]
    for (idx, sequence) in enumerate(reference_sequences)
        if length(sequence) >= k
            push!(reference_records, FASTX.FASTA.Record("ref_$(idx)", sequence))
        end
    end
    if isempty(reference_records)
        return (tp = 0, fp = 0, fn = 0, precision = 1.0, recall = 1.0, f1 = 1.0)
    end

    assembly_records = FASTX.FASTA.Record[]
    for (idx, sequence) in enumerate(assembled_contigs)
        if length(sequence) >= k
            push!(assembly_records, FASTX.FASTA.Record("asm_$(idx)", sequence))
        end
    end

    reference_kmers = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, reference_records)
    assembly_kmer_keys = if isempty(assembly_records)
        typeof(first(keys(reference_kmers)))[]
    else
        collect(keys(Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, assembly_records)))
    end

    return Mycelia.presence_precision_recall_f1(
        collect(keys(reference_kmers)),
        assembly_kmer_keys
    )
end

function compute_assembly_metrics(reference_sequences::Vector{String},
        assembled_contigs::Vector{String}; k::Int = 31)
    contig_lengths = length.(assembled_contigs)
    total_contig_length = sum(contig_lengths)
    total_truth_length = sum(length.(reference_sequences))
    n_stats = n50_l50(contig_lengths)
    kmer_scores = compute_kmer_presence_metrics(
        reference_sequences,
        assembled_contigs;
        k = k
    )
    size_ratio = total_truth_length == 0 ? 0.0 : total_contig_length / total_truth_length

    return (
        num_contigs = length(assembled_contigs),
        total_contig_length = total_contig_length,
        n50 = n_stats.n50,
        l50 = n_stats.l50,
        assembly_size_ratio = size_ratio,
        kmer_precision = kmer_scores.precision,
        kmer_recall = kmer_scores.recall,
        kmer_f1 = kmer_scores.f1
    )
end

function add_profile_delta_columns!(results::DataFrames.DataFrame)
    if DataFrames.nrow(results) == 0
        return results
    end

    results[!, :f1_delta_vs_full] = Union{Missing, Float64}[
        missing for _ in 1:DataFrames.nrow(results)]
    results[!, :runtime_ratio_vs_full] = Union{Missing, Float64}[
        missing for _ in 1:DataFrames.nrow(results)]
    results[!, :memory_ratio_vs_full] = Union{Missing, Float64}[
        missing for _ in 1:DataFrames.nrow(results)]

    key_columns = filter(
        column_name -> column_name in DataFrames.names(results),
        [:dataset_type, :genome_id, :community_size, :q_value, :quality_profile,
            :coverage, :depth_target, :abundance_profile, :replicate, :graph_family,
            :input_mode, :tokenizer]
    )

    baseline_by_key = Dict{Tuple, NamedTuple{(:kmer_f1, :run_seconds, :peak_memory_bytes),
        Tuple{Float64, Float64, Float64}}}()

    for row in DataFrames.eachrow(results)
        if row.memory_profile == "full" && row.status == "ok"
            key = Tuple(row[c] for c in key_columns)
            baseline_by_key[key] = (
                kmer_f1 = Float64(row.kmer_f1),
                run_seconds = Float64(row.run_seconds),
                peak_memory_bytes = Float64(row.peak_memory_bytes)
            )
        end
    end

    for (row_index, row) in enumerate(DataFrames.eachrow(results))
        key = Tuple(row[c] for c in key_columns)
        if haskey(baseline_by_key, key) && row.status == "ok"
            baseline = baseline_by_key[key]
            results[row_index, :f1_delta_vs_full] = row.kmer_f1 - baseline.kmer_f1
            if baseline.run_seconds > 0
                results[row_index, :runtime_ratio_vs_full] =
                    row.run_seconds / baseline.run_seconds
            end
            if baseline.peak_memory_bytes > 0
                results[row_index, :memory_ratio_vs_full] =
                    row.peak_memory_bytes / baseline.peak_memory_bytes
            end
        end
    end

    return results
end
