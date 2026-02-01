import Test
import Mycelia

Test.@testset "Utility Functions" begin
    Test.@testset "File helpers" begin
        Test.@test Mycelia.strip_gz_extension("sample.tsv.gz") == "sample.tsv"
        Test.@test Mycelia.strip_gz_extension("sample.tsv") == "sample.tsv"

        temp_dir = mktempdir()
        nonempty_path = joinpath(temp_dir, "nonempty.txt")
        empty_path = joinpath(temp_dir, "empty.txt")

        open(nonempty_path, "w") do io
            write(io, "data")
        end
        close(open(empty_path, "w"))

        Test.@test Mycelia.nonempty_file(nonempty_path)
        Test.@test !Mycelia.nonempty_file(empty_path)
        Test.@test !Mycelia.nonempty_file(joinpath(temp_dir, "missing.txt"))
    end

    Test.@testset "Collapse duplicates" begin
        df = Mycelia.DataFrames.DataFrame(id = [1, 1, 2], a = [1, missing, 2], b = [
            "x", "x", "y"])
        collapsed = Mycelia.collapse_duplicates(df, :id)
        Test.@test Mycelia.DataFrames.nrow(collapsed) == 2
        Test.@test collapsed[collapsed.id .== 1, :a][1] == 1
        Test.@test collapsed[collapsed.id .== 1, :b][1] == "x"

        conflict_df = Mycelia.DataFrames.DataFrame(id = [1, 1], a = [1, 2])
        conflict_collapsed = Mycelia.collapse_duplicates(
            conflict_df, :id; warn_on_conflict = false, summarize_conflicts = false)
        Test.@test Mycelia.DataFrames.nrow(conflict_collapsed) == 2
    end

    Test.@testset "Retry helper" begin
        attempts = Ref(0)
        result = Mycelia.with_retry(max_attempts = 3, initial_delay = 0.0,
            backoff_factor = 1.0, log_on_retry = false) do
            attempts[] += 1
            attempts[] < 2 && error("fail")
            "ok"
        end
        Test.@test result == "ok"
        Test.@test attempts[] == 2

        Test.@test_throws ErrorException Mycelia.with_retry(
            max_attempts = 2, initial_delay = 0.0,
            log_on_retry = false, log_on_failure = false) do
            error("always")
        end
    end

    Test.@testset "Directory listing" begin
        temp_dir = mktempdir()
        nested_dir = joinpath(temp_dir, "a", "b")
        mkpath(nested_dir)
        file_path = joinpath(nested_dir, "file.txt")
        open(file_path, "w") do io
            write(io, "content")
        end

        dirs = Mycelia.recursively_list_directories(temp_dir)
        files = Mycelia.recursively_list_files(temp_dir)
        Test.@test nested_dir in dirs
        Test.@test file_path in files

        Test.@test isempty(Mycelia._get_output_files(joinpath(temp_dir, "missing")))
        Test.@test file_path in Mycelia._get_output_files(temp_dir)
    end

    Test.@testset "TSV GZ roundtrip" begin
        temp_dir = mktempdir()
        df = Mycelia.DataFrames.DataFrame(a = [1, 2], b = ["x", "y"])
        out_path = joinpath(temp_dir, "data.tsv.gz")
        written = Mycelia.write_tsvgz(df = df, filename = out_path)
        Test.@test written == out_path

        loaded = Mycelia.read_tsvgz(out_path)
        Test.@test Mycelia.DataFrames.nrow(loaded) == 2
        Test.@test loaded.a == df.a
        Test.@test loaded.b == df.b
    end

    Test.@testset "DataFrame utilities" begin
        df = Mycelia.DataFrames.DataFrame(a = [1, nothing], b = ["x", "y"])
        Mycelia.dataframe_replace_nothing_with_missing(df)
        Test.@test df.a[2] === missing

        dict_df = Mycelia.DataFrames.DataFrame(meta = [Dict("a" => 1), Dict("b" => 2)])
        json_df = Mycelia.dataframe_convert_dicts_to_json(dict_df)
        Test.@test json_df.meta[1] == Mycelia.JSON.json(Dict("a" => 1))

        downcast_df = Mycelia.DataFrames.DataFrame(a = [1.0, 2.0], b = [1.5, 2.5])
        Mycelia.downcast_float_columns(downcast_df; target_type = Float32)
        Test.@test eltype(downcast_df.a) == Float32
        Test.@test eltype(downcast_df.b) == Float32

        dictvec = [Dict("a" => 1, "b" => 2), Dict("b" => 3)]
        dict_df = Mycelia.dictvec_to_dataframe(dictvec)
        Test.@test :a in Symbol.(Mycelia.DataFrames.names(dict_df))
        Test.@test dict_df[2, :a] === missing
    end

    Test.@testset "JSON helpers" begin
        Test.@test Mycelia.normalize_json_value(missing) === nothing
        Test.@test Mycelia.normalize_json_value(:alpha) == "alpha"
        Test.@test Mycelia.normalize_json_value((a = 1, b = :two)) ==
                   Dict("a" => 1, "b" => "two")

        df = Mycelia.DataFrames.DataFrame(
            id = [1, 2],
            name = ["a", "b"],
            created = [Mycelia.Dates.DateTime(2024, 1, 2, 3, 4), missing]
        )
        ndjson = Mycelia.dataframe_to_ndjson(df)
        lines = split(ndjson, '\n')
        Test.@test length(lines) == 2
        parsed = Mycelia.JSON.parse(lines[1])
        Test.@test parsed["id"] == 1
        Test.@test parsed["name"] == "a"
    end

    Test.@testset "Math and memory helpers" begin
        dense_bytes = Mycelia.estimate_dense_matrix_memory(2, 3)
        Test.@test dense_bytes == Base.sizeof(Float64) * 2 * 3

        sparse_bytes = Mycelia.estimate_sparse_matrix_memory(2, 3; nnz = 4)
        expected_sparse = Base.sizeof(Float64) * 4 + Base.sizeof(Int) * 4 +
                          Base.sizeof(Int) * 4
        Test.@test sparse_bytes == expected_sparse

        mem_check = Mycelia.check_matrix_fits_in_memory(1)
        Test.@test mem_check.bytes_needed == 1
    end

    Test.@testset "Sampling and markers" begin
        markers = Mycelia.choose_top_n_markers(3)
        Test.@test length(markers) == 3
        Test.@test markers[1] == :circle
    end

    Test.@testset "Hashing and formatting" begin
        temp_dir = mktempdir()
        hash_path = joinpath(temp_dir, "hash.txt")
        open(hash_path, "w") do io
            write(io, "abc")
        end
        expected = Mycelia.SHA.bytes2hex(Mycelia.SHA.sha256(codeunits("abc")))
        Test.@test Mycelia.sha256_file(hash_path) == expected

        normalized = Mycelia.normalized_current_datetime()
        Test.@test !isempty(normalized)
        Test.@test occursin(r"^[A-Za-z0-9_]+$", normalized)
    end

    Test.@testset "Sequence utilities" begin
        Test.@test Mycelia.alphabet_to_biosequence_type(:DNA) ==
                   Mycelia.BioSequences.LongDNA{4}
        Test.@test Mycelia.alphabet_to_biosequence_type(:RNA) ==
                   Mycelia.BioSequences.LongRNA{4}
        Test.@test Mycelia.alphabet_to_biosequence_type(:AA) == Mycelia.BioSequences.LongAA
        Test.@test_throws ArgumentError Mycelia.alphabet_to_biosequence_type(:XYZ)

        record = Mycelia.FASTX.FASTA.Record("id", "ATGC")
        alphabet, typed_seq = Mycelia.detect_and_extract_sequence(record)
        Test.@test alphabet == :DNA
        Test.@test typed_seq isa Mycelia.BioSequences.LongDNA

        typed = Mycelia.detect_and_extract_sequence("AUGC", :RNA)
        Test.@test typed isa Mycelia.BioSequences.LongRNA

        Test.@test Mycelia.calculate_gc_content("ATGC") == 50.0
        records = [
            Mycelia.FASTX.FASTA.Record("r1", "ATGC"),
            Mycelia.FASTX.FASTA.Record("r2", "GG")
        ]
        gc_pct = Mycelia.calculate_gc_content(records)
        Test.@test isapprox(gc_pct, 66.6666667; atol = 1e-6)
    end

    Test.@testset "Ranges and sampling" begin
        Test.@test Mycelia.get_base_extension("sample.fna.gz") == ".fna"
        Test.@test Mycelia.get_base_extension("sample.txt") == ".txt"

        ranges = Mycelia.find_true_ranges([false, true, true, false, true]; min_length = 2)
        Test.@test ranges == [(2, 3)]

        samples = Mycelia.equally_spaced_samples(collect(1:10), 3)
        Test.@test samples == [1, 6, 10]

        avg = Mycelia.rolling_centered_avg([1, 2, 3, 4, 5]; window_size = 3)
        Test.@test avg[1] == 1.5
        Test.@test avg[3] == 3.0
    end

    Test.@testset "Distance matrix and types" begin
        matrix = Mycelia.random_symmetric_distance_matrix(4)
        Test.@test all(matrix .== matrix')
        Test.@test all(matrix[i, i] == 0.0 for i in 1:size(matrix, 1))

        Test.@test Mycelia.type_to_string(Int) == string(Int)
        Test.@test Mycelia.type_to_string(Mycelia.Kmers.DNAKmer{3}) == "Kmers.DNAKmer{3}"
        Test.@test Mycelia.kmer_space_size(3, 4) == 64
    end

    Test.@testset "DataFrame structure helpers" begin
        df = Mycelia.DataFrames.DataFrame(a = [missing, missing], b = ["", ""], c = [
            "x", ""])
        nonempty = Mycelia.find_nonempty_columns(df)
        Test.@test nonempty == [false, false, true]

        dropped = Mycelia.drop_empty_columns(df)
        Test.@test Symbol.(Mycelia.DataFrames.names(dropped)) == [:c]

        df_copy = Mycelia.DataFrames.copy(df)
        Mycelia.drop_empty_columns!(df_copy)
        Test.@test Symbol.(Mycelia.DataFrames.names(df_copy)) == [:c]
    end

    Test.@testset "JSONL parsing" begin
        temp_dir = mktempdir()
        jsonl_path = joinpath(temp_dir, "data.jsonl")
        open(jsonl_path, "w") do io
            println(io, "{\"a\":1,\"b\":\"x\"}")
            println(io, "{\"a\":2}")
        end

        rows = Mycelia.parse_jsonl(jsonl_path)
        Test.@test length(rows) == 2
        Test.@test rows[1]["a"] == 1

        df = Mycelia.jsonl_to_dataframe(jsonl_path)
        Test.@test Mycelia.DataFrames.nrow(df) == 2
        Test.@test df[2, "b"] === missing
    end

    Test.@testset "Distribution helpers" begin
        normalized = Mycelia.normalize_countmap(Dict("a" => 2, "b" => 1))
        Test.@test isapprox(normalized["a"], 2 / 3; atol = 1e-8)
        Test.@test isapprox(normalized["b"], 1 / 3; atol = 1e-8)

        sha_values = [repeat("a", 64), repeat("b", 64)]
        meta1 = Mycelia.metasha256(sha_values)
        meta2 = Mycelia.metasha256(reverse(sha_values))
        Test.@test meta1 == meta2
        Test.@test length(meta1) == 64
    end

    Test.@testset "Filesystem helpers" begin
        temp_dir = mktempdir()
        source_path = joinpath(temp_dir, "source.txt")
        open(source_path, "w") do io
            write(io, "payload")
        end
        target_dir = joinpath(temp_dir, "out")
        mkpath(target_dir)
        copied = Mycelia.copy_with_unique_identifier(source_path, target_dir, "id")
        Test.@test isfile(copied)
        Test.@test read(copied, String) == "payload"

        size_str = Mycelia.filesize_human_readable(copied)
        Test.@test !isempty(size_str)

        Test.@test Mycelia.dir_size(temp_dir) > 0
    end

    Test.@testset "String helpers" begin
        Test.@test Mycelia.scientific_notation(12.34; precision = 2) == "1.23e+01"
        Test.@test Mycelia.find_matching_prefix("sample-01.fastq", "sample-02.fastq") ==
                   "sample-0"
    end

    Test.@testset "Math utilities" begin
        Test.@test Mycelia.nearest_prime(10) == 11
        Test.@test Mycelia.nearest_prime(12) == 11
        Test.@test Mycelia.fibonacci_numbers_less_than(10) == [0, 1, 1, 2, 3, 5, 8]
        Test.@test Mycelia.fibonacci_numbers_less_than(1) == [0]
        Test.@test Mycelia.fibonacci_numbers_less_than(0) == Int[]
    end

    Test.@testset "Path selection" begin
        temp_dir = mktempdir()
        missing_path = joinpath(temp_dir, "missing.txt")
        existing_path = joinpath(temp_dir, "exists.txt")
        open(existing_path, "w") do io
            write(io, "data")
        end
        Test.@test Mycelia.select_first_existing([missing_path, existing_path]) ==
                   existing_path
        Test.@test Mycelia.select_first_existing([missing_path]) === nothing
    end
end
