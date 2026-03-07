import BioSequences
import CSV
import DataFrames
import FASTX
import Mycelia
import StableRNGs

const STANDARD_ASSEMBLER_FIXTURE_SPECS = Dict(
    "synthetic_isolate_5386" => (
        kind = :isolate,
        description = "Single synthetic short-read isolate with fixed 5.4 kb genome length.",
        genomes = [
            (name = "synthetic_isolate_5386", length = 5_386, seed = 17_401, coverage = 30.0)
        ],
        read_length = 150,
        insert_size = 300
    ),
    "synthetic_metagenome_pair" => (
        kind = :metagenome,
        description = "Two-genome low-complexity metagenome with deterministic uneven coverage.",
        genomes = [
            (name = "synthetic_metagenome_genome_1", length = 4_000, seed = 17_402, coverage = 24.0),
            (name = "synthetic_metagenome_genome_2", length = 3_000, seed = 17_403, coverage = 12.0)
        ],
        read_length = 150,
        insert_size = 320
    )
)

function list_standard_assembler_fixtures()
    rows = NamedTuple[]
    for fixture_name in sort(collect(keys(STANDARD_ASSEMBLER_FIXTURE_SPECS)))
        spec = STANDARD_ASSEMBLER_FIXTURE_SPECS[fixture_name]
        push!(rows, (
            name = fixture_name,
            kind = String(spec.kind),
            description = spec.description,
            genomes = length(spec.genomes),
            total_reference_bases = sum(genome.length for genome in spec.genomes),
            read_length = spec.read_length,
            insert_size = spec.insert_size
        ))
    end
    return DataFrames.DataFrame(rows)
end

function build_standard_assembler_benchmark_plan(;
        fixture_names = sort(collect(keys(STANDARD_ASSEMBLER_FIXTURE_SPECS))),
        assemblers = ["Rhizomorph", "MEGAHIT", "metaSPAdes"])
    rows = NamedTuple[]
    for fixture_name in fixture_names
        @assert haskey(STANDARD_ASSEMBLER_FIXTURE_SPECS, fixture_name) "Unknown fixture: $(fixture_name)"
        spec = STANDARD_ASSEMBLER_FIXTURE_SPECS[fixture_name]
        for assembler in assemblers
            push!(rows, (
                fixture = fixture_name,
                fixture_kind = String(spec.kind),
                assembler = assembler,
                total_reference_bases = sum(genome.length for genome in spec.genomes)
            ))
        end
    end
    return DataFrames.DataFrame(rows)
end

function _fixture_reference_path(outdir::String, fixture_name::String)
    return joinpath(outdir, "$(fixture_name).reference.fasta")
end

function _fixture_fastq_paths(outdir::String, fixture_name::String)
    return (
        forward = joinpath(outdir, "$(fixture_name).R1.fastq"),
        reverse = joinpath(outdir, "$(fixture_name).R2.fastq")
    )
end

function _fixture_truth_table_path(outdir::String, fixture_name::String)
    return joinpath(outdir, "$(fixture_name).truth.csv")
end

function _quality_string(length_value::Int)
    return repeat("I", length_value)
end

function _generate_paired_fastq_records(genome_name::String, sequence, coverage::Real, read_length::Int,
        insert_size::Int)
    reads_1, reads_2 = Mycelia.generate_paired_end_reads(
        sequence,
        coverage,
        read_length,
        insert_size;
        error_rate = 0.0
    )

    forward_records = FASTX.FASTQ.Record[]
    reverse_records = FASTX.FASTQ.Record[]
    for (read_index, (read_1, read_2)) in enumerate(zip(reads_1, reads_2))
        read_id = "$(genome_name)_read_$(read_index)"
        push!(forward_records, FASTX.FASTQ.Record(
            read_id * "/1",
            read_1,
            _quality_string(length(read_1))
        ))
        push!(reverse_records, FASTX.FASTQ.Record(
            read_id * "/2",
            read_2,
            _quality_string(length(read_2))
        ))
    end
    return forward_records, reverse_records
end

function materialize_standard_assembler_fixture(fixture_name::String; outdir::String, emit_reads::Bool = true)
    @assert haskey(STANDARD_ASSEMBLER_FIXTURE_SPECS, fixture_name) "Unknown fixture: $(fixture_name)"
    spec = STANDARD_ASSEMBLER_FIXTURE_SPECS[fixture_name]

    mkpath(outdir)
    reference_fasta = _fixture_reference_path(outdir, fixture_name)
    fastq_paths = _fixture_fastq_paths(outdir, fixture_name)
    truth_table_path = _fixture_truth_table_path(outdir, fixture_name)

    fasta_records = FASTX.FASTA.Record[]
    forward_records = FASTX.FASTQ.Record[]
    reverse_records = FASTX.FASTQ.Record[]
    truth_rows = NamedTuple[]

    for genome in spec.genomes
        rng = StableRNGs.StableRNG(genome.seed)
        sequence = BioSequences.randdnaseq(rng, genome.length)
        push!(fasta_records, FASTX.FASTA.Record(genome.name, sequence))
        push!(truth_rows, (
            sequence_id = genome.name,
            length = genome.length,
            coverage = genome.coverage,
            seed = genome.seed
        ))

        if emit_reads
            generated_forward, generated_reverse = _generate_paired_fastq_records(
                genome.name,
                sequence,
                genome.coverage,
                spec.read_length,
                spec.insert_size
            )
            append!(forward_records, generated_forward)
            append!(reverse_records, generated_reverse)
        end
    end

    Mycelia.write_fasta(outfile = reference_fasta, records = fasta_records, gzip = false)
    CSV.write(truth_table_path, DataFrames.DataFrame(truth_rows))

    if emit_reads
        Mycelia.write_fastq(outfile = fastq_paths.forward, records = forward_records, gzip = false)
        Mycelia.write_fastq(outfile = fastq_paths.reverse, records = reverse_records, gzip = false)
    end

    return (
        name = fixture_name,
        kind = spec.kind,
        description = spec.description,
        reference_fasta = reference_fasta,
        fastq1 = emit_reads ? fastq_paths.forward : nothing,
        fastq2 = emit_reads ? fastq_paths.reverse : nothing,
        read_length = spec.read_length,
        insert_size = spec.insert_size,
        total_reference_bases = sum(genome.length for genome in spec.genomes),
        truth_table = DataFrames.DataFrame(truth_rows),
        truth_table_path = truth_table_path
    )
end
