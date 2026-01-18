import Test
import Mycelia

Test.@testset "FASTX Utilities" begin
    Test.@testset "Identifier extraction" begin
        path = "Sample_001_extra_long.fasta"
        extracted = Mycelia._extract_human_readable_id_from_fastx_path(path)
        Test.@test extracted == "Sample_001_extra"
    end

    Test.@testset "Hash helpers" begin
        sequences = ["ATGC", "AACC"]
        hash1 = Mycelia.generate_joint_sequence_hash(sequences; encoded_length=16)
        hash2 = Mycelia.generate_joint_sequence_hash(reverse(sequences); encoded_length=16)
        Test.@test hash1 == hash2

        hash_single = Mycelia.create_sequence_hash("atgc"; encoded_length=16)
        Test.@test length(hash_single) == 16
    end

    Test.@testset "Identifier normalization" begin
        Test.@test Mycelia.sanitize_fastx_identifier("id with spaces/") == "id_with_spaces_"

        id_map = Dict("foo" => "original")
        Test.@test Mycelia.normalize_reference_label("foo.fna", id_map) == "original"
        Test.@test Mycelia.normalize_reference_label("bar.msh", id_map) == "bar"
    end

    Test.@testset "FASTA and FASTQ IO" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "test.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGC"),
            Mycelia.FASTX.FASTA.Record("seq2", "GCTA")
        ]
        Mycelia.write_fasta(outfile=fasta_path, records=records, gzip=false, show_progress=false)

        reader = Mycelia.open_fastx(fasta_path)
        read_records = collect(reader)
        close(reader)
        Test.@test length(read_records) == 2
        Test.@test Mycelia.count_records(fasta_path) == 2

        fastq_path = joinpath(temp_dir, "reads.fastq")
        fastq_records = [
            Mycelia.fastq_record(identifier="r1", sequence="ATGC", quality_scores=[30, 30, 30, 30]),
            Mycelia.fastq_record(identifier="r2", sequence="GCTA", quality_scores=[20, 20, 20, 20])
        ]
        Mycelia.write_fastq(records=fastq_records, filename=fastq_path)
        Test.@test Mycelia.count_records(fastq_path) == 2
    end

    Test.@testset "FASTQ quality helpers" begin
        record = Mycelia.FASTX.FASTQ.Record("read1", "ATGC", "IIII")
        Test.@test Mycelia.get_phred_scores(record) == UInt8[40, 40, 40, 40]
        Test.@test Mycelia.quality_string_to_phred("!#%+") == UInt8[0, 2, 4, 10]

        error_rate = Mycelia.q_value_to_error_rate(20)
        Test.@test isapprox(Mycelia.error_rate_to_q_value(error_rate), 20.0; atol=1e-12)
    end

    Test.@testset "Sequence extension helpers" begin
        Test.@test Mycelia.alphabet_hint_from_path("sample.fna.gz") == :DNA
        Test.@test Mycelia.alphabet_hint_from_path("sample.fa") === nothing
        Test.@test Mycelia.detect_sequence_extension("ATGC") == ".fna"
        Test.@test Mycelia.detect_sequence_extension("AUGC") == ".frn"
    end

    Test.@testset "FASTA comparison" begin
        temp_dir = mktempdir()
        fasta_a = joinpath(temp_dir, "a.fna")
        fasta_b = joinpath(temp_dir, "b.fna")
        records_a = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGC"),
            Mycelia.FASTX.FASTA.Record("seq2", "GCTA")
        ]
        records_b = [
            Mycelia.FASTX.FASTA.Record("seq2", "GCTA"),
            Mycelia.FASTX.FASTA.Record("seq1", "ATGC")
        ]
        Mycelia.write_fasta(outfile=fasta_a, records=records_a, gzip=false, show_progress=false)
        Mycelia.write_fasta(outfile=fasta_b, records=records_b, gzip=false, show_progress=false)
        Test.@test Mycelia.equivalent_fasta_sequences(fasta_a, fasta_b)
    end

    Test.@testset "Find FASTA files" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "seqs.fna")
        other_path = joinpath(temp_dir, "notes.txt")
        Mycelia.write_fasta(outfile=fasta_path, records=[Mycelia.FASTX.FASTA.Record("seq", "ATGC")], gzip=false, show_progress=false)
        open(other_path, "w") do io
            write(io, "ignore")
        end
        found = Mycelia.find_fasta_files(temp_dir)
        Test.@test fasta_path in found
        Test.@test !(other_path in found)
    end

    Test.@testset "Join FASTQs with UUID" begin
        temp_dir = mktempdir()
        fastq_a = joinpath(temp_dir, "a.fastq")
        fastq_b = joinpath(temp_dir, "b.fastq")
        records_a = [
            Mycelia.fastq_record(identifier="r1", sequence="ATGC", quality_scores=[30, 30, 30, 30])
        ]
        records_b = [
            Mycelia.fastq_record(identifier="r2", sequence="GCTA", quality_scores=[20, 20, 20, 20])
        ]
        Mycelia.write_fastq(records=records_a, filename=fastq_a)
        Mycelia.write_fastq(records=records_b, filename=fastq_b)

        fastq_out = joinpath(temp_dir, "joint_reads.fq.gz")
        tsv_out = joinpath(temp_dir, "joint_reads.tsv.gz")
        outputs = Mycelia.join_fastqs_with_uuid([fastq_a, fastq_b]; fastq_out=fastq_out, tsv_out=tsv_out)
        Test.@test Mycelia.nonempty_file(outputs.fastq_out)
        Test.@test Mycelia.nonempty_file(outputs.tsv_out)

        mapping = Mycelia.read_tsvgz(outputs.tsv_out)
        Test.@test Mycelia.DataFrames.nrow(mapping) == 2
        Test.@test Mycelia.count_records(outputs.fastq_out) == 2
    end

    Test.@testset "Normalized tables and streams" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "normalized.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGC"),
            Mycelia.FASTX.FASTA.Record("seq2", "GCTA")
        ]
        Mycelia.write_fasta(outfile=fasta_path, records=records, gzip=false, show_progress=false)

        table = Mycelia.fastx2normalized_table(fasta_path; human_readable_id="sample")
        Test.@test Mycelia.DataFrames.nrow(table) == 2
        Test.@test all(table[!, "human_readable_id"] .== "sample")
        Test.@test "fastx_identifier" in String.(Mycelia.DataFrames.names(table))

        out_path = Mycelia.fastx2normalized_jsonl_stream(fastx_path=fasta_path, human_readable_id="sample", show_progress=false)
        Test.@test Mycelia.nonempty_file(out_path)

        io = Mycelia.CodecZlib.GzipDecompressorStream(open(out_path, "r"))
        lines = readlines(io)
        close(io)
        Test.@test length(lines) == 2
        parsed = Mycelia.JSON.parse(lines[1])
        Test.@test parsed["file_type"] == "fasta"
        Test.@test parsed["human_readable_id"] == "sample"
    end
end
