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
        hash1 = Mycelia.generate_joint_sequence_hash(sequences; encoded_length = 16)
        hash2 = Mycelia.generate_joint_sequence_hash(reverse(sequences); encoded_length = 16)
        Test.@test hash1 == hash2

        hash_single = Mycelia.create_sequence_hash("atgc"; encoded_length = 16)
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
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        reader = Mycelia.open_fastx(fasta_path)
        read_records = collect(reader)
        close(reader)
        Test.@test length(read_records) == 2
        Test.@test Mycelia.count_records(fasta_path) == 2

        fastq_path = joinpath(temp_dir, "reads.fastq")
        fastq_records = [
            Mycelia.fastq_record(identifier = "r1", sequence = "ATGC", quality_scores = [
                30, 30, 30, 30]),
            Mycelia.fastq_record(identifier = "r2", sequence = "GCTA", quality_scores = [
                20, 20, 20, 20])
        ]
        Mycelia.write_fastq(records = fastq_records, filename = fastq_path)
        Test.@test Mycelia.count_records(fastq_path) == 2
    end

    Test.@testset "FASTQ quality helpers" begin
        record = Mycelia.FASTX.FASTQ.Record("read1", "ATGC", "IIII")
        Test.@test Mycelia.get_phred_scores(record) == UInt8[40, 40, 40, 40]
        Test.@test Mycelia.quality_string_to_phred("!#%+") == UInt8[0, 2, 4, 10]

        error_rate = Mycelia.q_value_to_error_rate(20)
        Test.@test isapprox(Mycelia.error_rate_to_q_value(error_rate), 20.0; atol = 1e-12)
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
        Mycelia.write_fasta(outfile = fasta_a, records = records_a, gzip = false, show_progress = false)
        Mycelia.write_fasta(outfile = fasta_b, records = records_b, gzip = false, show_progress = false)
        Test.@test Mycelia.equivalent_fasta_sequences(fasta_a, fasta_b)
    end

    Test.@testset "Find FASTA files" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "seqs.fna")
        other_path = joinpath(temp_dir, "notes.txt")
        Mycelia.write_fasta(
            outfile = fasta_path, records = [Mycelia.FASTX.FASTA.Record("seq", "ATGC")],
            gzip = false, show_progress = false)
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
            Mycelia.fastq_record(identifier = "r1", sequence = "ATGC", quality_scores = [
            30, 30, 30, 30])
        ]
        records_b = [
            Mycelia.fastq_record(identifier = "r2", sequence = "GCTA", quality_scores = [
            20, 20, 20, 20])
        ]
        Mycelia.write_fastq(records = records_a, filename = fastq_a)
        Mycelia.write_fastq(records = records_b, filename = fastq_b)

        fastq_out = joinpath(temp_dir, "joint_reads.fq.gz")
        tsv_out = joinpath(temp_dir, "joint_reads.tsv.gz")
        outputs = Mycelia.join_fastqs_with_uuid([fastq_a, fastq_b]; fastq_out = fastq_out, tsv_out = tsv_out)
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
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        table = Mycelia.fastx2normalized_table(fasta_path; human_readable_id = "sample")
        Test.@test Mycelia.DataFrames.nrow(table) == 2
        Test.@test all(table[!, "human_readable_id"] .== "sample")
        Test.@test "fastx_identifier" in String.(Mycelia.DataFrames.names(table))

        out_path = Mycelia.fastx2normalized_jsonl_stream(
            fastx_path = fasta_path, human_readable_id = "sample", show_progress = false)
        Test.@test Mycelia.nonempty_file(out_path)

        io = Mycelia.CodecZlib.GzipDecompressorStream(open(out_path, "r"))
        lines = readlines(io)
        close(io)
        Test.@test length(lines) == 2
        parsed = Mycelia.JSON.parse(lines[1])
        Test.@test parsed["file_type"] == "fasta"
        Test.@test parsed["human_readable_id"] == "sample"
    end

    Test.@testset "_fmt_hms helper" begin
        Test.@test Mycelia._fmt_hms(0) == "00:00:00"
        Test.@test Mycelia._fmt_hms(-5) == "00:00:00"
        Test.@test Mycelia._fmt_hms(61) == "00:01:01"
        Test.@test Mycelia._fmt_hms(3661) == "01:01:01"
    end

    Test.@testset "_default_out_path helper" begin
        Test.@test Mycelia._default_out_path("reads.fna") == "reads.normalized.jsonl.gz"
        Test.@test Mycelia._default_out_path("reads.fna.gz") == "reads.normalized.jsonl.gz"
        Test.@test Mycelia._default_out_path("reads.fastq") == "reads.normalized.jsonl.gz"
        Test.@test Mycelia._default_out_path("reads.fq.gz") == "reads.normalized.jsonl.gz"
        Test.@test_throws ErrorException Mycelia._default_out_path("reads.txt")
    end

    Test.@testset "_extract_human_readable_id edge cases" begin
        # Short name - returned as-is
        Test.@test Mycelia._extract_human_readable_id_from_fastx_path("short.fna") ==
                   "short"

        # FASTQ extension
        Test.@test Mycelia._extract_human_readable_id_from_fastx_path("reads.fastq") ==
                   "reads"
        Test.@test Mycelia._extract_human_readable_id_from_fastx_path("reads.fq.gz") ==
                   "reads"

        # Long name with delimiters - finds longest prefix ≤16
        Test.@test length(Mycelia._extract_human_readable_id_from_fastx_path(
            "Sample_001_extra_long_name.fna")) <= 16

        # Meaningless prefix bypass (e.g. GCA accession)
        id = Mycelia._extract_human_readable_id_from_fastx_path("GCA_000001405.fna")
        Test.@test length(id) <= 16
        Test.@test length(id) >= 3

        # force_truncate on very long name with no delimiters
        long_name = "abcdefghijklmnopqrstuvwxyz.fna"
        Test.@test length(Mycelia._extract_human_readable_id_from_fastx_path(long_name, true)) ==
                   16

        # Error without force_truncate on very long name with no delimiters
        Test.@test_throws ErrorException Mycelia._extract_human_readable_id_from_fastx_path(
            long_name, false)

        # Non-FASTX extension errors
        Test.@test_throws ErrorException Mycelia._extract_human_readable_id_from_fastx_path(
            "file.txt")
    end

    Test.@testset "normalized_table2fastx round-trip" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "input.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("s1", "ATGCATGC"),
            Mycelia.FASTX.FASTA.Record("s2", "GCTAGCTA")
        ]
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        table = Mycelia.fastx2normalized_table(fasta_path; human_readable_id = "rt")
        out_fasta = Mycelia.normalized_table2fastx(table; output_dir = temp_dir)
        Test.@test isfile(out_fasta)
        Test.@test Mycelia.count_records(out_fasta) == 2

        # Test force=false skips existing
        out_fasta2 = Mycelia.normalized_table2fastx(table; output_dir = temp_dir)
        Test.@test out_fasta2 == out_fasta

        # Test with explicit output_basename
        out_fasta3 = Mycelia.normalized_table2fastx(
            table; output_dir = temp_dir, output_basename = "custom", force = true)
        Test.@test endswith(out_fasta3, "custom.fna")
        Test.@test isfile(out_fasta3)
    end

    Test.@testset "split_fasta_by_record" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "multi.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("contig1", "ATGCATGC"),
            Mycelia.FASTX.FASTA.Record("contig2", "GCTAGCTA")
        ]
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        outdir = joinpath(temp_dir, "split")
        id_map = Mycelia.split_fasta_by_record(fasta_in = fasta_path, outdir = outdir)
        Test.@test length(id_map) == 2
        Test.@test haskey(id_map, "contig1")
        Test.@test haskey(id_map, "contig2")
        Test.@test id_map["contig1"] == "contig1"

        for (safe_id, _) in id_map
            Test.@test isfile(joinpath(outdir, safe_id * ".fna"))
        end

        # Re-run should skip existing files (force=false)
        id_map2 = Mycelia.split_fasta_by_record(fasta_in = fasta_path, outdir = outdir)
        Test.@test id_map2 == id_map
    end

    Test.@testset "subset_fasta_by_ids" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "all.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("s1", "AAAA"),
            Mycelia.FASTX.FASTA.Record("s2", "CCCC"),
            Mycelia.FASTX.FASTA.Record("s3", "GGGG")
        ]
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        out_path = joinpath(temp_dir, "subset.fna")
        result = Mycelia.subset_fasta_by_ids(
            fasta_in = fasta_path, ids = ["s1", "s3"], fasta_out = out_path)
        Test.@test isfile(result)
        Test.@test Mycelia.count_records(result) == 2

        # require_all=true should error on missing ids
        Test.@test_throws ErrorException Mycelia.subset_fasta_by_ids(
            fasta_in = fasta_path, ids = ["s1", "missing"], fasta_out = joinpath(temp_dir, "fail.fna"),
            force = true)

        # require_all=false should not error
        out2 = joinpath(temp_dir, "partial.fna")
        result2 = Mycelia.subset_fasta_by_ids(
            fasta_in = fasta_path, ids = ["s1", "missing"], fasta_out = out2,
            require_all = false)
        Test.@test Mycelia.count_records(result2) == 1

        # Existing output is returned without re-processing (force=false)
        result3 = Mycelia.subset_fasta_by_ids(
            fasta_in = fasta_path, ids = ["s1", "s3"], fasta_out = out_path)
        Test.@test result3 == out_path
    end

    Test.@testset "concatenate_fastq_files" begin
        temp_dir = mktempdir()
        fq1 = joinpath(temp_dir, "a.fastq")
        fq2 = joinpath(temp_dir, "b.fastq")
        Mycelia.write_fastq(
            records = [Mycelia.fastq_record(
                identifier = "r1", sequence = "ATGC", quality_scores = [30, 30, 30, 30])],
            filename = fq1)
        Mycelia.write_fastq(
            records = [Mycelia.fastq_record(
                identifier = "r2", sequence = "GCTA", quality_scores = [20, 20, 20, 20])],
            filename = fq2)

        concat_out = joinpath(temp_dir, "concat.fastq")
        result = Mycelia.concatenate_fastq_files(
            fastq_files = [fq1, fq2], output_fastq = concat_out)
        Test.@test result == concat_out
        Test.@test Mycelia.count_records(concat_out) == 2

        # force=false returns existing
        result2 = Mycelia.concatenate_fastq_files(
            fastq_files = [fq1, fq2], output_fastq = concat_out)
        Test.@test result2 == concat_out
    end

    Test.@testset "create_sequence_hash dispatches" begin
        # Test BioSequences.LongSequence dispatch
        dna_seq = Mycelia.BioSequences.LongDNA{4}("ATCGATCG")
        hash_bio = Mycelia.create_sequence_hash(dna_seq; encoded_length = 16)
        hash_str = Mycelia.create_sequence_hash("ATCGATCG"; encoded_length = 16)
        Test.@test hash_bio == hash_str
        Test.@test length(hash_bio) == 16

        # Test different hash functions (all except crc32 can produce 16 chars)
        for hf in [:sha256, :sha512, :md5, :sha1, :sha3_256, :sha3_512]
            h = Mycelia.create_sequence_hash(
                "ATGC"; hash_function = hf, encoded_length = 16,
                allow_truncation = true)
            Test.@test length(h) == 16
        end
        # CRC32 produces only 4 bytes (6 base58 chars)
        h_crc = Mycelia.create_sequence_hash(
            "ATGC"; hash_function = :crc32, encoded_length = 6,
            allow_truncation = true)
        Test.@test length(h_crc) == 6

        # Unsupported hash function
        Test.@test_throws ErrorException Mycelia.create_sequence_hash(
            "ATGC"; hash_function = :bogus)

        # Empty sequence
        Test.@test_throws ErrorException Mycelia.create_sequence_hash("")
    end

    Test.@testset "create_base58_hash" begin
        h = Mycelia.create_base58_hash("ATGCATGC"; encoded_length = 32)
        Test.@test length(h) == 32
        # Case insensitive by default
        Test.@test Mycelia.create_base58_hash("atgc") == Mycelia.create_base58_hash("ATGC")
    end

    Test.@testset "generate_joint_sequence_hash edge cases" begin
        # Empty vector should error
        Test.@test_throws ErrorException Mycelia.generate_joint_sequence_hash(String[])

        # Order independence
        h1 = Mycelia.generate_joint_sequence_hash(["AAAA", "CCCC"]; encoded_length = 16)
        h2 = Mycelia.generate_joint_sequence_hash(["CCCC", "AAAA"]; encoded_length = 16)
        Test.@test h1 == h2
    end

    Test.@testset "find_fasta_files error paths" begin
        # Non-FASTA file input
        temp_dir = mktempdir()
        txt_file = joinpath(temp_dir, "notes.txt")
        open(txt_file, "w") do io
            write(io, "not fasta")
        end
        Test.@test_throws ErrorException Mycelia.find_fasta_files(txt_file)

        # Non-existent path
        Test.@test_throws ErrorException Mycelia.find_fasta_files("/nonexistent/path/xyz")

        # Single FASTA file input
        fasta_file = joinpath(temp_dir, "single.fna")
        Mycelia.write_fasta(outfile = fasta_file,
            records = [Mycelia.FASTX.FASTA.Record("s", "ATGC")],
            gzip = false, show_progress = false)
        found = Mycelia.find_fasta_files(fasta_file)
        Test.@test found == [fasta_file]
    end

    Test.@testset "write_fastq validation" begin
        # No filename provided
        Test.@test_throws ErrorException Mycelia.write_fastq(
            records = [Mycelia.fastq_record(
            identifier = "r1", sequence = "ATGC", quality_scores = [30, 30, 30, 30])])

        # Invalid extension
        Test.@test_throws ErrorException Mycelia.write_fastq(
            records = [Mycelia.fastq_record(
                identifier = "r1", sequence = "ATGC", quality_scores = [30, 30, 30, 30])],
            filename = "/tmp/bad.txt")

        # Conflicting filename and outfile
        Test.@test_throws ErrorException Mycelia.write_fastq(
            records = [Mycelia.fastq_record(
                identifier = "r1", sequence = "ATGC", quality_scores = [30, 30, 30, 30])],
            filename = "/tmp/a.fastq", outfile = "/tmp/b.fastq")
    end

    Test.@testset "FASTQ gzipped IO" begin
        temp_dir = mktempdir()
        gz_path = joinpath(temp_dir, "reads.fq.gz")
        recs = [
            Mycelia.fastq_record(identifier = "r1", sequence = "ATGC", quality_scores = [
            30, 30, 30, 30])
        ]
        Mycelia.write_fastq(records = recs, filename = gz_path)
        Test.@test Mycelia.count_records(gz_path) == 1
    end

    Test.@testset "total_fasta_size" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "size.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("s1", "AAAA"),
            Mycelia.FASTX.FASTA.Record("s2", "CCCCCC")
        ]
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)
        Test.@test Mycelia.total_fasta_size(fasta_path) == 10
    end

    Test.@testset "fasta_to_table and fasta_table_to_fasta round-trip" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "table_rt.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("id1", "ATGC"),
            Mycelia.FASTX.FASTA.Record("id2", "GCTA")
        ]
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        reader = Mycelia.open_fastx(fasta_path)
        df = Mycelia.fasta_to_table(reader)
        close(reader)
        Test.@test Mycelia.DataFrames.nrow(df) == 2
        Test.@test "identifier" in Mycelia.DataFrames.names(df)

        rt_records = Mycelia.fasta_table_to_fasta(df)
        Test.@test length(rt_records) == 2
    end

    Test.@testset "deduplicate_fasta_file" begin
        temp_dir = mktempdir()
        in_fasta = joinpath(temp_dir, "dups.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("s1", "ATGC"),
            Mycelia.FASTX.FASTA.Record("s2", "ATGC"),
            Mycelia.FASTX.FASTA.Record("s3", "GCTA")
        ]
        Mycelia.write_fasta(outfile = in_fasta, records = records, gzip = false, show_progress = false)

        out_fasta = joinpath(temp_dir, "unique.fna")
        result = Mycelia.deduplicate_fasta_file(in_fasta, out_fasta)
        Test.@test result == out_fasta
        Test.@test Mycelia.count_records(out_fasta) == 2
    end

    Test.@testset "fastx_to_contig_lengths" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "contigs.fna")
        records = [
            Mycelia.FASTX.FASTA.Record("c1", "AAAA"),
            Mycelia.FASTX.FASTA.Record("c2", "CCCCCC")
        ]
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        lengths = Mycelia.fastx_to_contig_lengths(fasta_path)
        Test.@test lengths["c1"] == 4
        Test.@test lengths["c2"] == 6
    end

    Test.@testset "fastx2normalized_table with FASTQ" begin
        temp_dir = mktempdir()
        fastq_path = joinpath(temp_dir, "norm.fastq")
        recs = [
            Mycelia.fastq_record(identifier = "r1", sequence = "ATGC", quality_scores = [
                30, 30, 30, 30]),
            Mycelia.fastq_record(identifier = "r2", sequence = "GCTA", quality_scores = [
                20, 20, 20, 20])
        ]
        Mycelia.write_fastq(records = recs, filename = fastq_path)

        table = Mycelia.fastx2normalized_table(fastq_path; human_readable_id = "fq")
        Test.@test Mycelia.DataFrames.nrow(table) == 2
        Test.@test !any(ismissing, table[!, "record_quality"])
        Test.@test !any(ismissing, table[!, "mean_record_quality"])
    end

    Test.@testset "fastx2normalized_jsonl_stream with FASTQ" begin
        temp_dir = mktempdir()
        fastq_path = joinpath(temp_dir, "stream.fastq")
        recs = [
            Mycelia.fastq_record(identifier = "r1", sequence = "ATGC", quality_scores = [
            30, 30, 30, 30])
        ]
        Mycelia.write_fastq(records = recs, filename = fastq_path)

        out_path = Mycelia.fastx2normalized_jsonl_stream(
            fastx_path = fastq_path, human_readable_id = "fq", show_progress = false)
        Test.@test Mycelia.nonempty_file(out_path)
        io = Mycelia.CodecZlib.GzipDecompressorStream(open(out_path, "r"))
        lines = readlines(io)
        close(io)
        parsed = Mycelia.JSON.parse(lines[1])
        Test.@test parsed["file_type"] == "fastq"
        Test.@test parsed["record_quality"] !== nothing
    end

    Test.@testset "fastx2normalized_jsonl_stream non-gz output" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "plain.fna")
        records = [Mycelia.FASTX.FASTA.Record("s1", "ATGC")]
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        out_path = joinpath(temp_dir, "output.jsonl")
        result = Mycelia.fastx2normalized_jsonl_stream(
            fastx_path = fasta_path, human_readable_id = "plain",
            output_path = out_path, show_progress = false)
        Test.@test result == out_path
        Test.@test isfile(out_path)
        lines = readlines(out_path)
        Test.@test length(lines) == 1
    end

    Test.@testset "fastx2normalized_table force_truncate" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "trunc.fna")
        records = [Mycelia.FASTX.FASTA.Record("s1", "ATGC")]
        Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false, show_progress = false)

        # Too-long explicit ID without force_truncate should error
        Test.@test_throws ErrorException Mycelia.fastx2normalized_table(
            fasta_path; human_readable_id = "this_is_way_too_long_id")

        # With force_truncate it should work
        table = Mycelia.fastx2normalized_table(
            fasta_path; human_readable_id = "this_is_way_too_long_id", force_truncate = true)
        Test.@test length(table[1, "human_readable_id"]) == 16
    end

    Test.@testset "alphabet_hint_from_path extended" begin
        Test.@test Mycelia.alphabet_hint_from_path("sample.faa") == :AA
        Test.@test Mycelia.alphabet_hint_from_path("sample.frn.gz") == :RNA
        Test.@test Mycelia.alphabet_hint_from_path("sample.fna.gz") == :DNA
        Test.@test Mycelia.alphabet_hint_from_path("sample.fastq") === nothing
        Test.@test Mycelia.alphabet_hint_from_path("sample.txt") === nothing
    end

    Test.@testset "detect_sequence_extension extended" begin
        # Protein sequence
        Test.@test Mycelia.detect_sequence_extension("MKTAYIAK") == ".faa"
        # From FASTA record
        rec = Mycelia.FASTX.FASTA.Record("s1", "ATGCATGC")
        Test.@test Mycelia.detect_sequence_extension(rec) == ".fna"
    end

    Test.@testset "equivalent_fasta_sequences non-equal" begin
        temp_dir = mktempdir()
        fa = joinpath(temp_dir, "a.fna")
        fb = joinpath(temp_dir, "b.fna")
        Mycelia.write_fasta(outfile = fa,
            records = [Mycelia.FASTX.FASTA.Record("s1", "ATGC")],
            gzip = false, show_progress = false)
        Mycelia.write_fasta(outfile = fb,
            records = [Mycelia.FASTX.FASTA.Record("s1", "GGGG")],
            gzip = false, show_progress = false)
        Test.@test !Mycelia.equivalent_fasta_sequences(fa, fb)
    end

    Test.@testset "translate_nucleic_acid_fasta" begin
        temp_dir = mktempdir()
        in_fasta = joinpath(temp_dir, "dna.fna")
        records = [Mycelia.FASTX.FASTA.Record("gene1", "ATGATGATG")]
        Mycelia.write_fasta(outfile = in_fasta, records = records, gzip = false, show_progress = false)

        out_fasta = joinpath(temp_dir, "protein.faa")
        result = Mycelia.translate_nucleic_acid_fasta(in_fasta, out_fasta)
        Test.@test result == out_fasta
        Test.@test isfile(out_fasta)
        Test.@test Mycelia.count_records(out_fasta) == 1
    end
end
