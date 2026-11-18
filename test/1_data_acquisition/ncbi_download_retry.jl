# Offline retry coverage for NCBI download helpers.
#
# ```bash
# julia --project=. -e 'include("test/1_data_acquisition/ncbi_download_retry.jl")'
# ```

import Dates
import Logging
import Test
import Mycelia

if !isdefined(Main, :test_throws_message)
    include(joinpath(dirname(@__DIR__), "test_helpers.jl"))
end

function with_temp_home(f::Function)
    original_home = get(ENV, "HOME", nothing)
    mktempdir() do dir
        ENV["HOME"] = dir
        try
            return f(dir)
        finally
            if isnothing(original_home)
                delete!(ENV, "HOME")
            else
                ENV["HOME"] = original_home
            end
        end
    end
end

function fasta_records(header::AbstractString, sequence::AbstractString)
    reader = Mycelia.FASTX.FASTA.Reader(IOBuffer(">" * header * "\n" * sequence * "\n"))
    try
        return collect(reader)
    finally
        close(reader)
    end
end

function ncbi_summary_fixture(;
        assembly_accession::AbstractString = "GCF_000001",
        ftp_path::AbstractString = "ftp://example/GCF_000001",
        taxid::AbstractString = "562",
        gc_percent::AbstractString = "50.0")
    header = [
        "assembly_accession",
        "bioproject",
        "biosample",
        "wgs_master",
        "refseq_category",
        "organism_name",
        "infraspecific_name",
        "isolate",
        "version_status",
        "assembly_level",
        "release_type",
        "genome_rep",
        "seq_rel_date",
        "asm_name",
        "asm_submitter",
        "gbrs_paired_asm",
        "paired_asm_comp",
        "ftp_path",
        "excluded_from_refseq",
        "relation_to_type_material",
        "asm_not_live_date",
        "assembly_type",
        "group",
        "annotation_provider",
        "annotation_name",
        "annotation_date",
        "pubmed_id",
        "taxid",
        "species_taxid",
        "genome_size",
        "genome_size_ungapped",
        "replicon_count",
        "scaffold_count",
        "contig_count",
        "total_gene_count",
        "protein_coding_gene_count",
        "non_coding_gene_count",
        "gc_percent"
    ]
    row = [
        assembly_accession,
        "PRJNA1",
        "SAMN1",
        "WGS1",
        "reference genome",
        "Escherichia coli",
        "strain=K12",
        "K12",
        "latest",
        "Complete Genome",
        "Major",
        "Full",
        "2023-01-01",
        "ASM1",
        "Submitter",
        "none",
        "none",
        ftp_path,
        "none",
        "none",
        "none",
        "haploid",
        "bacteria",
        "Provider",
        "Annot1",
        "2023-01-02",
        "12345",
        taxid,
        taxid,
        "5000",
        "4800",
        "1",
        "1",
        "1",
        "4500",
        "4300",
        "200",
        gc_percent
    ]
    return join([
        "## Assembly summary test",
        "# " * join(header, '\t'),
        join(row, '\t'),
        ""
    ], '\n')
end

function write_fake_datasets_runner(dir::AbstractString; succeed_on::Int)
    logdir = joinpath(dir, "logs")
    mkpath(logdir)
    runner = joinpath(dir, "conda")
    template = raw"""
    #!/bin/sh
    logdir="__LOGDIR__"
    succeed_on=__SUCCEED_ON__
    printf '%s\n' "$*" >> "$logdir/commands.log"
    case " $* " in
      *" datasets "*)
        count_file="$logdir/datasets_attempts"
        n=0
        if [ -f "$count_file" ]; then
          n=$(cat "$count_file")
        fi
        n=$((n + 1))
        printf '%s\n' "$n" > "$count_file"
        prev=""
        for arg in "$@"; do
          if [ "$prev" = "--inputfile" ]; then
            cp "$arg" "$logdir/last_input.txt"
            printf '%s\n' "$arg" > "$logdir/last_input_path.txt"
          fi
          prev=$arg
        done
        if [ "$n" -lt "$succeed_on" ]; then
          echo "transient taxonomy failure" >&2
          exit 1
        fi
        printf '%s\n' '{"tax_id":562}'
        ;;
      *" dataformat "*)
        cat >/dev/null
        printf 'tax_id\trank\n562\tspecies\n'
        ;;
      *)
        echo "unexpected conda command: $*" >&2
        exit 1
        ;;
    esac
    """
    script = replace(template, "__LOGDIR__" => logdir)
    script = replace(script, "__SUCCEED_ON__" => string(succeed_on))
    script = join(strip.(split(script, '\n'; keepempty = false)), '\n') * "\n"
    write(runner, script)
    chmod(runner, 0o755)
    return (runner = runner, logdir = logdir)
end

Test.@testset "NCBI download retry" begin
    Test.@testset "download_genome_by_accession retries the fetch seam" begin
        mktempdir() do dir
            attempts = Ref(0)
            fetcher = function ()
                attempts[] += 1
                attempts[] < 3 && error("transient nucleotide failure")
                return fasta_records("NC_TEST", "ACGTACGT")
            end

            outfile = Test.@test_logs (
                :warn, r"NCBI nucleotide download attempt 1/3",
            ) min_level=Logging.Warn match_mode=:any begin
                Mycelia.download_genome_by_accession(
                    accession = "NC_TEST",
                    outdir = dir,
                    compressed = false,
                    max_attempts = 3,
                    initial_retry_delay = 0.0,
                    record_fetcher = fetcher
                )
            end

            Test.@test outfile == joinpath(dir, "NC_TEST.fna")
            Test.@test attempts[] == 3
            Test.@test occursin("ACGTACGT", read(outfile, String))
            Test.@test occursin(">NC_TEST", read(outfile, String))

            Mycelia.download_genome_by_accession(
                accession = "NC_TEST",
                outdir = dir,
                compressed = false,
                max_attempts = 3,
                initial_retry_delay = 0.0,
                record_fetcher = fetcher
            )
            Test.@test attempts[] == 3

            write(outfile, "")
            Mycelia.download_genome_by_accession(
                accession = "NC_TEST",
                outdir = dir,
                compressed = false,
                max_attempts = 3,
                initial_retry_delay = 0.0,
                record_fetcher = () -> fasta_records("NC_TEST", "GGGG")
            )
            Test.@test occursin("GGGG", read(outfile, String))
        end
    end

    Test.@testset "download_genome_by_accession writes gzip and rethrows" begin
        mktempdir() do dir
            gz = Mycelia.download_genome_by_accession(
                accession = "NC_GZIP",
                outdir = dir,
                compressed = true,
                max_attempts = 1,
                initial_retry_delay = 0.0,
                record_fetcher = () -> fasta_records("NC_GZIP", "TTTT")
            )
            Test.@test gz == joinpath(dir, "NC_GZIP.fna.gz")
            Test.@test filesize(gz) > 0
            open(gz) do io
                Test.@test read(io, 2) == UInt8[0x1f, 0x8b]
            end

            attempts = Ref(0)
            test_throws_message(ErrorException, "nucleotide download failed hard") do
                Mycelia.download_genome_by_accession(
                    accession = "NC_FAIL",
                    outdir = dir,
                    compressed = false,
                    max_attempts = 2,
                    initial_retry_delay = 0.0,
                    record_fetcher = function ()
                        attempts[] += 1
                        error("nucleotide download failed hard")
                    end
                )
            end
            Test.@test attempts[] == 2
            Test.@test !isfile(joinpath(dir, "NC_FAIL.fna"))
        end
    end

    Test.@testset "load_ncbi_metadata retries the downloader" begin
        test_throws_message(ArgumentError, "Invalid database") do
            Mycelia.load_ncbi_metadata("invalid-db"; initial_retry_delay = 0.0)
        end

        fixture = ncbi_summary_fixture(
            assembly_accession = "GCF_RETRY",
            ftp_path = "ftp://example/GCF_RETRY",
            taxid = "111"
        )

        with_temp_home() do home
            attempts = Ref(0)
            seen_urls = String[]
            downloader = function (url, dest)
                attempts[] += 1
                push!(seen_urls, String(url))
                if attempts[] == 1
                    write(dest, "partial")
                    error("transient metadata failure")
                elseif attempts[] == 2
                    write(dest, "")
                    return dest
                end
                write(dest, fixture)
                return dest
            end

            df = Test.@test_logs (
                :warn, r"NCBI metadata download attempt 1/3",
            ) min_level=Logging.Warn match_mode=:any begin
                Mycelia.load_ncbi_metadata(
                    "refseq";
                    max_attempts = 3,
                    initial_retry_delay = 0.0,
                    downloader = downloader
                )
            end
            Test.@test attempts[] == 3
            Test.@test df[1, "assembly_accession"] == "GCF_RETRY"
            Test.@test df[1, "taxid"] == 111
            Test.@test df[1, "ftp_path"] == "ftp://example/GCF_RETRY"
            Test.@test length(seen_urls) == 3
            Test.@test all(url -> occursin("ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt", url), seen_urls)
            cache_path = joinpath(
                home, "workspace", ".ncbi",
                Dates.format(Dates.today(), "yyyy-mm-dd") * ".assembly_summary_refseq.txt"
            )
            Test.@test occursin("GCF_RETRY", read(cache_path, String))
            Test.@test !occursin("partial", read(cache_path, String))

            cached_attempts = Ref(0)
            cached = Mycelia.load_ncbi_metadata(
                "refseq";
                initial_retry_delay = 0.0,
                downloader = function (url, dest)
                    cached_attempts[] += 1
                    error("cache hit should not download")
                end
            )
            Test.@test cached_attempts[] == 0
            Test.@test cached[1, "assembly_accession"] == "GCF_RETRY"
        end

        with_temp_home() do home
            cache_dir = joinpath(home, "workspace", ".ncbi")
            mkpath(cache_dir)
            empty_cache = joinpath(
                cache_dir,
                Dates.format(Dates.today(), "yyyy-mm-dd") * ".assembly_summary_genbank.txt"
            )
            write(empty_cache, "")
            attempts = Ref(0)
            seen_url = Ref("")
            df = Test.@test_logs (
                :warn, r"Cached file exists but is empty",
            ) min_level=Logging.Warn match_mode=:any begin
                Mycelia.load_ncbi_metadata(
                    "genbank";
                    initial_retry_delay = 0.0,
                    summary_url = "http://127.0.0.1/assembly_summary_genbank.txt",
                    downloader = function (url, dest)
                        attempts[] += 1
                        seen_url[] = String(url)
                        write(dest, ncbi_summary_fixture(
                            assembly_accession = "GCA_EMPTY_REPLACED",
                            taxid = "222"
                        ))
                        return dest
                    end
                )
            end
            Test.@test attempts[] == 1
            Test.@test seen_url[] == "http://127.0.0.1/assembly_summary_genbank.txt"
            Test.@test df[1, "assembly_accession"] == "GCA_EMPTY_REPLACED"
            Test.@test df[1, "taxid"] == 222
        end

        with_temp_home() do home
            attempts = Ref(0)
            test_throws_message(ErrorException, "Could not load NCBI metadata") do
                Mycelia.load_ncbi_metadata(
                    "refseq";
                    max_attempts = 2,
                    initial_retry_delay = 0.0,
                    downloader = function (url, dest)
                        attempts[] += 1
                        write(dest, "partial")
                        error("metadata downloader failed")
                    end
                )
            end
            Test.@test attempts[] == 2
            cache_path = joinpath(
                home, "workspace", ".ncbi",
                Dates.format(Dates.today(), "yyyy-mm-dd") * ".assembly_summary_refseq.txt"
            )
            Test.@test !isfile(cache_path)
        end

        with_temp_home() do home
            attempts = Ref(0)
            test_throws_message(ErrorException, "No valid cache file") do
                Mycelia.load_ncbi_metadata(
                    "refseq";
                    date = Dates.Date(2024, 2, 4),
                    initial_retry_delay = 0.0,
                    downloader = function (url, dest)
                        attempts[] += 1
                        error("dated load must not download")
                    end
                )
            end
            Test.@test attempts[] == 0
        end
    end

    Test.@testset "taxids2ncbi_taxonomy_table retries the datasets CLI" begin
        mktempdir() do dir
            fake = write_fake_datasets_runner(dir; succeed_on = 3)
            table = Test.@test_logs (
                :warn, r"NCBI taxonomy summary attempt 1/3",
            ) min_level=Logging.Warn match_mode=:any begin
                Mycelia.taxids2ncbi_taxonomy_table(
                    [10, 2, 10];
                    max_attempts = 3,
                    initial_retry_delay = 0.0,
                    conda_runner = fake.runner
                )
            end
            Test.@test Mycelia.DataFrames.nrow(table) == 1
            Test.@test string(table[1, "tax_id"]) == "562"
            Test.@test table[1, "rank"] == "species"
            Test.@test read(joinpath(fake.logdir, "datasets_attempts"), String) == "3\n"
            Test.@test read(joinpath(fake.logdir, "last_input.txt"), String) == "2\n10\n"
            input_path = strip(read(joinpath(fake.logdir, "last_input_path.txt"), String))
            Test.@test !isfile(input_path)
            commands = read(joinpath(fake.logdir, "commands.log"), String)
            Test.@test occursin("datasets summary taxonomy taxon --inputfile", commands)
            Test.@test occursin("dataformat tsv taxonomy --template tax-summary", commands)
            Test.@test occursin("ncbi-datasets-cli", commands)
            Test.@test length(filter(line -> occursin(" datasets ", " " * line * " "), split(commands, '\n'))) == 3
        end

        mktempdir() do dir
            fake = write_fake_datasets_runner(dir; succeed_on = 100)
            test_throws_message(ProcessFailedException, "failed") do
                Mycelia.taxids2ncbi_taxonomy_table(
                    [562];
                    max_attempts = 2,
                    initial_retry_delay = 0.0,
                    conda_runner = fake.runner
                )
            end
            Test.@test read(joinpath(fake.logdir, "datasets_attempts"), String) == "2\n"
            commands = read(joinpath(fake.logdir, "commands.log"), String)
            Test.@test !occursin("dataformat", commands)
            input_path = strip(read(joinpath(fake.logdir, "last_input_path.txt"), String))
            Test.@test !isfile(input_path)
        end
    end
end
