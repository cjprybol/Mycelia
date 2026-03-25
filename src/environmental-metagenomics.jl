const ENVIRONMENTAL_METAGENOME_METADATA_PATH = joinpath(
    @__DIR__, "..", "test", "metadata", "20250702.sra-public-fastqs.csv.gz")
const ENVIRONMENTAL_METAGENOME_CASE_STUDY_ID = :yukon_river_pilot_station

_normalized_metadata_value(value) = value === missing ? "" : strip(string(value))

function _parse_metadata_bytes(value)
    normalized = _normalized_metadata_value(value)
    isempty(normalized) && return 0
    try
        return parse(Int, normalized)
    catch
        return 0
    end
end

function _default_environmental_sample_name(path::AbstractString)
    basename_no_gz = replace(basename(path), r"\.gz$" => "")
    return replace(basename_no_gz, r"\.(fastq|fq)$" => "")
end

function _read_environmental_metadata_table(
        metadata_path::AbstractString = ENVIRONMENTAL_METAGENOME_METADATA_PATH)
    isfile(metadata_path) || error("Environmental metadata file not found: $(metadata_path)")

    if endswith(metadata_path, ".gz")
        return open(metadata_path) do io
            CSV.read(CodecZlib.GzipDecompressorStream(io), DataFrames.DataFrame)
        end
    end

    return CSV.read(metadata_path, DataFrames.DataFrame)
end

"""
    environmental_metagenome_case_study(; metadata_path=ENVIRONMENTAL_METAGENOME_METADATA_PATH, max_runs=3)

Select a lightweight public environmental metagenome cohort from the bundled SRA
metadata table.

The current case study targets paired-end WGS runs from Yukon River water at
Pilot Station, Alaska, USA. Runs are sorted by size so the smallest samples can
be explored first on local hardware.
"""
function environmental_metagenome_case_study(;
        metadata_path::AbstractString = ENVIRONMENTAL_METAGENOME_METADATA_PATH,
        max_runs::Int = 3
)
    max_runs > 0 || throw(ArgumentError("max_runs must be positive, got $(max_runs)"))

    metadata = _read_environmental_metadata_table(metadata_path)

    library_source = lowercase.(_normalized_metadata_value.(metadata[!, "LibrarySource"]))
    assay_type = lowercase.(_normalized_metadata_value.(metadata[!, "Assay Type"]))
    library_selection = lowercase.(_normalized_metadata_value.(metadata[!, "LibrarySelection"]))
    isolation_source = lowercase.(_normalized_metadata_value.(metadata[!, "isolation_source"]))
    geo_loc_name = lowercase.(_normalized_metadata_value.(metadata[!, "geo_loc_name"]))

    mask = (library_source .== "metagenomic") .&
           (assay_type .== "wgs") .&
           (library_selection .== "random") .&
           occursin.("river water", isolation_source) .&
           occursin.("yukon river", geo_loc_name)

    selected = metadata[mask, :]
    DataFrames.nrow(selected) > 0 || error(
        "No environmental metagenome runs matched the Yukon River case-study filter.")

    selected = DataFrames.select(
        selected,
        "Run" => DataFrames.ByRow(_normalized_metadata_value) => :run,
        "LibraryLayout" => DataFrames.ByRow(_normalized_metadata_value) => :library_layout,
        "LibrarySelection" => DataFrames.ByRow(_normalized_metadata_value) => :library_selection,
        "Instrument" => DataFrames.ByRow(_normalized_metadata_value) => :instrument,
        "geo_loc_name" => DataFrames.ByRow(_normalized_metadata_value) => :geo_loc_name,
        "isolation_source" => DataFrames.ByRow(_normalized_metadata_value) => :isolation_source,
        "Bytes" => DataFrames.ByRow(_parse_metadata_bytes) => :bytes
    )
    selected[!, :bytes_gb] = round.(selected.bytes ./ 1_000_000_000; digits = 3)
    DataFrames.sort!(selected, :bytes)
    selected = first(selected, min(max_runs, DataFrames.nrow(selected)))

    return (
        case_study_id = ENVIRONMENTAL_METAGENOME_CASE_STUDY_ID,
        title = "Yukon River Water Environmental Metagenome Cohort",
        habitat = "freshwater river water",
        location = "Yukon River at Pilot Station, Alaska, USA",
        rationale = "Bundled public SRA metadata already contains a small WGS river-water cohort suitable for local QC and k-mer ordination.",
        metadata_path = metadata_path,
        samples = selected
    )
end

"""
    environmental_fastq_quality_summary(fastq_files; sample_names=nothing)

Summarize per-sample FASTQ quality metrics for an environmental metagenome
cohort.
"""
function environmental_fastq_quality_summary(
        fastq_files::AbstractVector{<:AbstractString};
        sample_names = nothing
)
    isempty(fastq_files) && throw(ArgumentError("fastq_files must not be empty"))
    all(isfile, fastq_files) || error("All FASTQ files must exist")

    if isnothing(sample_names)
        sample_names = map(_default_environmental_sample_name, fastq_files)
    end

    length(sample_names) == length(fastq_files) ||
        throw(ArgumentError("sample_names length must match fastq_files length"))

    summary = DataFrames.DataFrame(
        sample = String[],
        fastq = String[],
        n_reads = Int[],
        mean_quality = Float64[],
        mean_length = Float64[],
        gc_content = Float64[],
        q20_percent = Float64[],
        q30_percent = Float64[],
        q40_percent = Float64[]
    )

    for (sample, fastq_file) in zip(sample_names, fastq_files)
        quality = Mycelia.analyze_fastq_quality(fastq_file)
        push!(summary, (
            sample = string(sample),
            fastq = string(fastq_file),
            n_reads = quality.n_reads,
            mean_quality = quality.mean_quality,
            mean_length = quality.mean_length,
            gc_content = quality.gc_content,
            q20_percent = quality.quality_distribution.q20_percent,
            q30_percent = quality.quality_distribution.q30_percent,
            q40_percent = quality.quality_distribution.q40_percent
        ))
    end

    return summary
end

function _count_canonical_kmers_from_fastq(
        fastq_file::AbstractString,
        k::Int;
        max_reads::Union{Nothing, Int} = nothing
)
    k > 0 || throw(ArgumentError("k must be positive, got $(k)"))
    if !isnothing(max_reads)
        max_reads > 0 || throw(ArgumentError("max_reads must be positive, got $(max_reads)"))
    end

    KMER_TYPE = Core.apply_type(Kmers.DNAKmer, k)
    counts = OrderedCollections.OrderedDict{KMER_TYPE, Int}()
    reads_used = 0

    reader = Mycelia.open_fastx(fastq_file)
    reader isa FASTX.FASTQ.Reader ||
        error("Expected FASTQ input for environmental analysis: $(fastq_file)")

    try
        for record in reader
            record_counts = Mycelia.count_kmers(KMER_TYPE, record)
            Mycelia.canonicalize_kmer_counts!(record_counts)
            merge!(+, counts, record_counts)

            reads_used += 1
            if !isnothing(max_reads) && reads_used >= max_reads
                break
            end
        end
    finally
        close(reader)
    end

    reads_used > 0 || error("No reads found in FASTQ file: $(fastq_file)")
    return (counts = sort!(counts), reads_used = reads_used)
end

"""
    analyze_environmental_kmer_profiles(fastq_files; sample_names=nothing, k=21, max_reads=5000)

Build lightweight k-mer profiles from a cohort of FASTQ files, compute a
Jaccard distance matrix, and return a PCoA embedding for sample-level
comparison.
"""
function analyze_environmental_kmer_profiles(
        fastq_files::AbstractVector{<:AbstractString};
        sample_names = nothing,
        k::Int = 21,
        max_reads::Union{Nothing, Int} = 5_000
)
    length(fastq_files) >= 2 ||
        throw(ArgumentError("Need at least two FASTQ files for cohort comparison"))
    all(isfile, fastq_files) || error("All FASTQ files must exist")

    if isnothing(sample_names)
        sample_names = map(_default_environmental_sample_name, fastq_files)
    end

    length(sample_names) == length(fastq_files) ||
        throw(ArgumentError("sample_names length must match fastq_files length"))

    counts_by_sample = OrderedCollections.OrderedDict{String, Any}()
    spectrum_summary = DataFrames.DataFrame(
        sample = String[],
        reads_used = Int[],
        total_kmers = Int[],
        unique_kmers = Int[],
        peak_coverage = Union{Missing, Int}[],
        peak_kmers = Int[]
    )

    for (sample, fastq_file) in zip(sample_names, fastq_files)
        result = _count_canonical_kmers_from_fastq(fastq_file, k; max_reads = max_reads)
        spectrum = Mycelia.analyze_kmer_frequency_spectrum(result.counts)
        counts_by_sample[string(sample)] = result.counts
        push!(spectrum_summary, (
            sample = string(sample),
            reads_used = result.reads_used,
            total_kmers = spectrum.total_kmers,
            unique_kmers = spectrum.unique_kmers,
            peak_coverage = spectrum.peak.coverage,
            peak_kmers = spectrum.peak.kmers
        ))
    end

    all_kmers = Set{Any}()
    for counts in values(counts_by_sample)
        union!(all_kmers, keys(counts))
    end
    !isempty(all_kmers) || error("No k-mers were observed across the provided FASTQ files")

    sorted_kmers = sort!(collect(all_kmers), by = string)
    kmer_index = Dict{Any, Int}(kmer => i for (i, kmer) in enumerate(sorted_kmers))
    counts_matrix = zeros(Int, length(sorted_kmers), length(sample_names))

    for (sample_index, sample_name) in enumerate(sample_names)
        for (kmer, count) in counts_by_sample[string(sample_name)]
            counts_matrix[kmer_index[kmer], sample_index] = count
        end
    end

    distance_matrix = Mycelia.frequency_matrix_to_jaccard_distance_matrix(counts_matrix)
    pcoa_result = Mycelia.pcoa_from_dist(distance_matrix; maxoutdim = 2)
    pcoa_df = Mycelia.pcoa_to_dataframe(pcoa_result, collect(string.(sample_names)))

    return (
        k = k,
        sample_names = collect(string.(sample_names)),
        counts_by_sample = counts_by_sample,
        counts_matrix = counts_matrix,
        distance_matrix = distance_matrix,
        pcoa = pcoa_result,
        pcoa_df = pcoa_df,
        spectrum_summary = spectrum_summary
    )
end
