import Test
import Mycelia
import CodecZlib
import Tar

const UN_REQUEST_FIXTURE_ROOT = Ref("")

@eval Mycelia begin
    function _un_request(url::AbstractString, method::AbstractString;
            headers::Dict{String, String} = Dict{String, String}(),
            max_redirects::Int = 5
    )
        filename = split(String(url), '/')[end]
        if startswith(filename, "unexpected")
            return (status = 500,)
        end

        path = joinpath(Main.UN_REQUEST_FIXTURE_ROOT[], filename)
        exists = isfile(path)
        if method == "HEAD"
            if !exists
                return (status = 404,)
            elseif endswith(filename, ".00")
                return (status = 405,)
            elseif endswith(filename, ".01")
                return (status = 403,)
            else
                return (status = 200,)
            end
        elseif method == "GET"
            if !exists
                return (status = 404,)
            elseif get(headers, "Range", "") == "bytes=0-0"
                if endswith(filename, ".00")
                    return (status = 206,)
                elseif endswith(filename, ".01")
                    return (status = 403,)
                end
            end
            return (status = 200,)
        end

        return (status = 500,)
    end
end

function write_tar_gz(outfile::AbstractString, files::Dict{String, String})
    mktempdir() do dir
        for (relative_path, content) in files
            path = joinpath(dir, relative_path)
            mkpath(dirname(path))
            write(path, content)
        end
        open(outfile, "w") do raw
            gzip_stream = CodecZlib.GzipCompressorStream(raw)
            try
                Tar.create(dir, gzip_stream)
            finally
                close(gzip_stream)
            end
        end
    end
    return outfile
end

function split_archive_file(infile::AbstractString, first_part::AbstractString, second_part::AbstractString)
    bytes = read(infile)
    midpoint = cld(length(bytes), 2)
    write(first_part, bytes[1:midpoint])
    write(second_part, bytes[(midpoint + 1):end])
    return (first_part, second_part)
end

function file_url(path::AbstractString)
    return "file://" * abspath(path)
end

Test.@testset "Reference Database Archive Helpers" begin
    mktempdir() do fixture_dir
        UN_REQUEST_FIXTURE_ROOT[] = fixture_dir

        direct_archive = joinpath(fixture_dir, "UNv1.0.en-fr.tar.gz")
        write_tar_gz(direct_archive, Dict("direct/en-fr.en" => "hello world\n"))

        split_source = joinpath(fixture_dir, "split-source.tar.gz")
        write_tar_gz(split_source, Dict("split/en-fr.fr" => "bonjour monde\n"))
        split_archive_file(
            split_source,
            joinpath(fixture_dir, "UNv1.0.fr-ru.tar.gz.00"),
            joinpath(fixture_dir, "UNv1.0.fr-ru.tar.gz.01")
        )
        rm(split_source; force = true)

        testset_archive = joinpath(fixture_dir, "UNv1.0.testsets.tar.gz")
        write_tar_gz(testset_archive, Dict("testsets/dev.en" => "dev line\n"))

        Test.@test Mycelia._un_part_exists(file_url(joinpath(fixture_dir, "UNv1.0.fr-ru.tar.gz.00")))
        Test.@test Mycelia._un_part_exists(file_url(joinpath(fixture_dir, "UNv1.0.fr-ru.tar.gz.01")))
        Test.@test Mycelia._un_part_exists(file_url(direct_archive))
        Test.@test !Mycelia._un_part_exists(file_url(joinpath(fixture_dir, "missing.tar.gz")))
        Test.@test_throws ErrorException Mycelia._un_part_exists(file_url(joinpath(fixture_dir, "unexpected.tar.gz")))

        mktempdir() do outdir
            split_result = Mycelia._download_and_extract_split_archive(
                "UNv1.0.fr-ru.tar.gz",
                file_url(fixture_dir),
                outdir
            )
            Test.@test split_result == outdir
            Test.@test isfile(joinpath(outdir, "split", "en-fr.fr"))
            Test.@test !isfile(joinpath(outdir, "UNv1.0.fr-ru.tar.gz"))
            Test.@test !isfile(joinpath(outdir, "UNv1.0.fr-ru.tar.gz.00"))
            Test.@test !isfile(joinpath(outdir, "UNv1.0.fr-ru.tar.gz.01"))
        end

        mktempdir() do outdir
            direct_result = Mycelia._download_and_extract_split_archive(
                "UNv1.0.en-fr.tar.gz",
                file_url(fixture_dir),
                outdir
            )
            Test.@test direct_result == outdir
            Test.@test isfile(joinpath(outdir, "direct", "en-fr.en"))
            Test.@test !isfile(joinpath(outdir, "UNv1.0.en-fr.tar.gz"))
        end

        mktempdir() do outdir
            existing_archive = joinpath(outdir, "UNv1.0.existing.tar.gz")
            write_tar_gz(existing_archive, Dict("existing/cache.txt" => "cached\n"))
            reused = Mycelia._download_and_extract_split_archive(
                "UNv1.0.existing.tar.gz",
                file_url(fixture_dir),
                outdir
            )
            Test.@test reused == outdir
            Test.@test isfile(joinpath(outdir, "existing", "cache.txt"))
            Test.@test isfile(existing_archive)
        end

        mktempdir() do outdir
            Test.@test_throws ErrorException Mycelia._download_and_extract_split_archive(
                "UNv1.0.missing.tar.gz",
                file_url(fixture_dir),
                outdir
            )
        end

        mktempdir() do outdir
            corpus_dir = Mycelia.download_un_parallel_corpus(
                outdir = outdir,
                subsets = ["en-fr"],
                formats = ["txt"],
                base_url = file_url(fixture_dir)
            )
            Test.@test corpus_dir == joinpath(outdir, "un_corpus")
            Test.@test isfile(joinpath(corpus_dir, "direct", "en-fr.en"))

            testset_dir = Mycelia.download_un_parallel_corpus_testset(
                outdir = outdir,
                base_url = file_url(fixture_dir)
            )
            Test.@test testset_dir == joinpath(outdir, "un_corpus_testsets")
            Test.@test isfile(joinpath(testset_dir, "testsets", "dev.en"))
        end
    end
end
