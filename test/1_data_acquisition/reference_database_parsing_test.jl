import Test
import Dates
import Logging
import Sockets
import Mycelia

const RUN_ALL = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true"
const RUN_EXTERNAL = RUN_ALL || lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"

function install_noop_add_bioconda_env!()
    @eval Mycelia begin
        function add_bioconda_env(pkg; force = false, quiet = false)
            return nothing
        end
    end
    return nothing
end

function restore_add_bioconda_env!()
    @eval Mycelia begin
        function add_bioconda_env(pkg; force = false, quiet = false)
            _ensure_conda_env_vars!()
            channel = nothing
            if occursin("::", pkg)
                if !quiet
                    println("splitting $(pkg)")
                end
                channel, pkg = split(pkg, "::")
                if !quiet
                    println("into channel:$(channel) pkg:$(pkg)")
                end
            end
            already_installed = check_bioconda_env_is_installed(pkg)
            if !already_installed || force
                if !quiet
                    @info "installing conda environment $(pkg)"
                end
                if isnothing(channel)
                    if quiet
                        run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y --quiet`)
                    else
                        run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
                    end
                else
                    if quiet
                        run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(channel)::$(pkg) -y --quiet`)
                    else
                        run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(channel)::$(pkg) -y`)
                    end
                end
                if quiet
                    run(`$(CONDA_RUNNER) clean --all -y --quiet`)
                else
                    run(`$(CONDA_RUNNER) clean --all -y`)
                end
            end
        end
    end
    return nothing
end

function write_executable(path::AbstractString, contents::AbstractString)
    write(path, contents)
    chmod(path, 0o755)
    return path
end

function install_fake_sra_runner!(runner_path::AbstractString)
    script = raw"""#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "run" ]]; then
    shift
fi
if [[ "${1:-}" == "--live-stream" ]]; then
    shift
fi
if [[ "${1:-}" == "-n" ]]; then
    shift 2
fi

tool="${1:-}"
if [[ -z "$tool" ]]; then
    echo "missing tool" >&2
    exit 1
fi
shift

case "$tool" in
    prefetch)
        srr=""
        outdir="."
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -O)
                    outdir="$2"
                    shift 2
                    ;;
                *)
                    srr="$1"
                    shift
                    ;;
            esac
        done
        if [[ "$srr" == *FAIL ]]; then
            exit 1
        fi
        mkdir -p "$outdir/$srr"
        printf 'stub archive' > "$outdir/$srr/$srr.sra"
        ;;
    fasterq-dump)
        outdir="."
        target=""
        while [[ $# -gt 0 ]]; do
            case "$1" in
                --outdir)
                    outdir="$2"
                    shift 2
                    ;;
                --mem|--split-3|--skip-technical)
                    shift
                    ;;
                --threads)
                    shift 2
                    ;;
                *)
                    target="$1"
                    shift
                    ;;
            esac
        done
        srr="$(basename "$target")"
        if [[ "$srr" == *FAIL ]]; then
            exit 1
        fi
        mkdir -p "$outdir"
        if [[ "$srr" == *PAIRED* ]]; then
            printf 'forward reads' > "$outdir/$srr"_1.fastq
            printf 'reverse reads' > "$outdir/$srr"_2.fastq
        elif [[ "$srr" == *SINGLE* ]]; then
            printf 'single reads' > "$outdir/$srr.fastq"
        fi
        ;;
    *)
        echo "unexpected tool: $tool" >&2
        exit 1
        ;;
esac
"""
    return write_executable(runner_path, script)
end

function install_fake_gzip!(gzip_path::AbstractString)
    script = raw"""#!/usr/bin/env bash
set -euo pipefail

input="${1:?missing input}"
cp "$input" "$input.gz"
rm -f "$input"
"""
    return write_executable(gzip_path, script)
end

function with_fake_sra_tooling(f::Function)
    original_path = get(ENV, "PATH", "")
    conda_runner = Mycelia.CONDA_RUNNER
    @assert isfile(conda_runner)

    mktempdir() do dir
        fake_runner = joinpath(dir, "conda")
        fake_gzip_dir = joinpath(dir, "bin")
        backup_runner = joinpath(dir, "conda.backup")

        mkpath(fake_gzip_dir)
        install_fake_sra_runner!(fake_runner)
        install_fake_gzip!(joinpath(fake_gzip_dir, "gzip"))

        cp(conda_runner, backup_runner; force = true)
        cp(fake_runner, conda_runner; force = true)
        chmod(conda_runner, 0o755)
        ENV["PATH"] = fake_gzip_dir * ":" * original_path

        try
            return f()
        finally
            ENV["PATH"] = original_path
            cp(backup_runner, conda_runner; force = true)
            chmod(conda_runner, 0o755)
        end
    end
end

function install_synthetic_sra_helpers!()
    @eval Mycelia begin
        function prefetch(; SRR, outdir = pwd())
            if endswith(SRR, "FAIL")
                error("synthetic prefetch failure for $(SRR)")
            end
            output_dir = joinpath(outdir, SRR)
            archive = joinpath(output_dir, "$(SRR).sra")
            mkpath(output_dir)
            write(archive, "stub archive")
            return (directory = output_dir, archive = archive)
        end

        function fasterq_dump(; outdir = pwd(), srr_identifier = "")
            if endswith(srr_identifier, "FAIL")
                error("synthetic fasterq failure for $(srr_identifier)")
            end
            mkpath(outdir)
            if occursin("PAIRED", srr_identifier)
                forward_reads = joinpath(outdir, "$(srr_identifier)_1.fastq.gz")
                reverse_reads = joinpath(outdir, "$(srr_identifier)_2.fastq.gz")
                write(forward_reads, "forward reads")
                write(reverse_reads, "reverse reads")
                return (
                    forward_reads = forward_reads,
                    reverse_reads = reverse_reads,
                    unpaired_reads = missing
                )
            elseif occursin("SINGLE", srr_identifier)
                unpaired_reads = joinpath(outdir, "$(srr_identifier).fastq.gz")
                write(unpaired_reads, "single reads")
                return (
                    forward_reads = missing,
                    reverse_reads = missing,
                    unpaired_reads = unpaired_reads
                )
            else
                return (
                    forward_reads = missing,
                    reverse_reads = missing,
                    unpaired_reads = missing
                )
            end
        end
    end
    return nothing
end

function restore_sra_helpers!()
    @eval Mycelia begin
        function prefetch(; SRR, outdir = pwd())
            Mycelia.add_bioconda_env("sra-tools")
            final_dir = joinpath(outdir, SRR)
            sra_archive = joinpath(final_dir, "$(SRR).sra")
            if !isfile(sra_archive)
                run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n sra-tools prefetch $(SRR) -O $(outdir)`)
            else
                @info "SRA archive already present: $(sra_archive)"
            end
            return (directory = final_dir, archive = sra_archive)
        end

        function fasterq_dump(; outdir = pwd(), srr_identifier = "")
            Mycelia.add_bioconda_env("sra-tools")
            prefetch_results = Mycelia.prefetch(SRR = srr_identifier, outdir = outdir)

            final_outdir = prefetch_results.directory

            forward_reads = joinpath(final_outdir, "$(srr_identifier)_1.fastq")
            reverse_reads = joinpath(final_outdir, "$(srr_identifier)_2.fastq")
            unpaired_reads = joinpath(final_outdir, "$(srr_identifier).fastq")

            forward_reads_gz = forward_reads * ".gz"
            reverse_reads_gz = reverse_reads * ".gz"
            unpaired_reads_gz = unpaired_reads * ".gz"

            forward_and_reverse_present = isfile(forward_reads_gz) && isfile(reverse_reads_gz)
            unpaired_present = isfile(unpaired_reads_gz)

            if !(forward_and_reverse_present || unpaired_present)
                fasterq_dump_cmd = `
                    $(Mycelia.CONDA_RUNNER) run --live-stream -n sra-tools fasterq-dump
                        --outdir $(final_outdir)
                        --mem 1G
                        --split-3
                        --threads $(min(get_default_threads(), 4))
                        --skip-technical
                        $(final_outdir)`
                @time run(fasterq_dump_cmd)
                isfile(forward_reads) && run(`gzip $(forward_reads)`)
                isfile(reverse_reads) && run(`gzip $(reverse_reads)`)
                isfile(unpaired_reads) && run(`gzip $(unpaired_reads)`)
            else
                @info "$(forward_reads_gz) & $(reverse_reads_gz) already present"
            end
            return (
                forward_reads = isfile(forward_reads_gz) ? forward_reads_gz : missing,
                reverse_reads = isfile(reverse_reads_gz) ? reverse_reads_gz : missing,
                unpaired_reads = isfile(unpaired_reads_gz) ? unpaired_reads_gz : missing
            )
        end
    end
    return nothing
end

function blast_env_available()
    if !isfile(Mycelia.CONDA_RUNNER)
        return false
    end
    try
        read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd -version`, String)
        read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --version`, String)
        return true
    catch
        return false
    end
end

function make_test_blastdb(dir::AbstractString; dbtype::AbstractString)
    source_path = if dbtype == "nucl"
        joinpath(dir, "test.fna")
    elseif dbtype == "prot"
        joinpath(dir, "test.faa")
    else
        error("Unsupported BLAST dbtype fixture: $dbtype")
    end

    if dbtype == "nucl"
        write(source_path, ">seq1\nACGTACGT\n>seq2\nTTTTAAAA\n")
        db_path = joinpath(dir, "nucl_db")
    else
        write(source_path, ">prot1\nMKT\n>prot2\nGGA\n")
        db_path = joinpath(dir, "prot_db")
    end

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast makeblastdb -in $(source_path) -dbtype $(dbtype) -parse_seqids -out $(db_path)`)
    return db_path
end

function reserve_local_port()
    server = Sockets.listen(Sockets.InetAddr(parse(Sockets.IPAddr, "127.0.0.1"), 0))
    port = last(Sockets.getsockname(server))
    close(server)
    return port
end

function make_test_tarball_bytes(files::Dict{String, String})
    mktempdir() do dir
        payload_dir = joinpath(dir, "payload")
        mkpath(payload_dir)
        for (relative_path, contents) in files
            output_path = joinpath(payload_dir, relative_path)
            mkpath(dirname(output_path))
            write(output_path, contents)
        end
        archive_path = joinpath(dir, "payload.tar.gz")
        run(`tar -czf $(archive_path) -C $(payload_dir) .`)
        return read(archive_path)
    end
end

function split_bytes(payload::Vector{UInt8}; chunks::Int = 2)
    chunk_count = max(1, chunks)
    chunk_size = cld(length(payload), chunk_count)
    return [
        payload[start_index:min(start_index + chunk_size - 1, length(payload))]
        for start_index in 1:chunk_size:length(payload)
    ]
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

function write_ncbi_metadata_fixture(
        path::AbstractString;
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

    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "## Assembly summary test")
        println(io, "# " * join(header, '\t'))
        println(io, join(row, '\t'))
    end

    return path
end

function start_un_test_server(file_map::Dict{String, Vector{UInt8}})
    port = reserve_local_port()
    base_url = "http://127.0.0.1:$(port)"

    function make_response(status::Integer, body = UInt8[]; headers = Pair{String, String}[])
        response_headers = ["Connection" => "close"]
        append!(response_headers, headers)
        return Mycelia.HTTP.Response(status, response_headers, body)
    end

    function handler(request)
        path = String(request.target)

        if path == "/redirect-start"
            return make_response(302; headers = ["Location" => "$(base_url)/redirect-target"])
        elseif path == "/redirect-target"
            return make_response(200, "redirected")
        elseif path == "/redirect-loop"
            return make_response(302; headers = ["Location" => "$(base_url)/redirect-loop"])
        elseif path == "/redirect-no-location"
            return make_response(302)
        elseif path == "/unexpected-status"
            return make_response(500, "server error")
        elseif path == "/range-fallback"
            if request.method == "HEAD"
                return make_response(405)
            end
            return make_response(206, file_map[path][1:min(length(file_map[path]), 1)])
        elseif path == "/range-missing-fallback"
            if request.method == "HEAD"
                return make_response(405)
            end
            return make_response(404, "missing after fallback")
        elseif path == "/forbidden-fallback"
            return make_response(403)
        end

        if !haskey(file_map, path)
            return make_response(404, "missing")
        end

        payload = file_map[path]
        if request.method == "HEAD"
            return make_response(200; headers = ["Content-Length" => string(length(payload))])
        end

        return make_response(200, payload)
    end

    server = Mycelia.HTTP.serve!(handler, "127.0.0.1", port)
    return (server = server, base_url = base_url)
end

Test.@testset "Reference Database Parsing" begin
    Test.@testset "NCBI metadata parsing" begin
        mktempdir() do dir
            path = joinpath(dir, "assembly_summary.txt")
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
                "GCF_000001",
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
                "ftp://example/GCF_000001",
                "none",
                "none",
                "none",
                "haploid",
                "bacteria",
                "Provider",
                "Annot1",
                "2023-01-02",
                "12345",
                "562",
                "562",
                "5000",
                "4800",
                "1",
                "1",
                "1",
                "4500",
                "4300",
                "200",
                "50.0"
            ]

            open(path, "w") do io
                println(io, "## Assembly summary test")
                println(io, "# " * join(header, '\t'))
                println(io, join(row, '\t'))
            end

            df = Mycelia.parse_ncbi_metadata_file(path)
            Test.@test Mycelia.DataFrames.nrow(df) == 1
            Test.@test df[1, "assembly_accession"] == "GCF_000001"
            Test.@test df[1, "taxid"] == 562
            Test.@test df[1, "gc_percent"] == 50.0
        end
    end

    Test.@testset "parse_ncbi_metadata_file error paths" begin
        Test.@test_throws ErrorException Mycelia.parse_ncbi_metadata_file("/nonexistent/path/file.txt")

        mktempdir() do dir
            empty_file = joinpath(dir, "empty.txt")
            touch(empty_file)
            Test.@test_throws ErrorException Mycelia.parse_ncbi_metadata_file(empty_file)
        end
    end

    Test.@testset "load_ncbi_metadata cache coverage" begin
        fixture_date = Dates.Date(2024, 2, 3)
        date_prefix = Dates.format(fixture_date, "yyyy-mm-dd")
        today_prefix = Dates.format(Dates.today(), "yyyy-mm-dd")

        Test.@test_throws ArgumentError Mycelia.load_ncbi_metadata("invalid-db"; date = fixture_date)

        with_temp_home() do home
            cache_dir = joinpath(home, "workspace", ".ncbi")

            requested_refseq_cache = write_ncbi_metadata_fixture(
                joinpath(cache_dir, "$(date_prefix).assembly_summary_refseq.txt");
                assembly_accession = "GCF_REQUESTED_REFSEQ",
                ftp_path = "ftp://example/GCF_REQUESTED_REFSEQ",
                taxid = "111"
            )

            requested_df = Test.@test_logs (:info, r"Found cached file for requested date") min_level=Logging.Info match_mode=:any begin
                Mycelia.load_ncbi_metadata("refseq"; date = fixture_date)
            end
            Test.@test requested_df[1, "assembly_accession"] == "GCF_REQUESTED_REFSEQ"
            Test.@test requested_df[1, "ftp_path"] == "ftp://example/GCF_REQUESTED_REFSEQ"
            Test.@test requested_refseq_cache == joinpath(cache_dir, "$(date_prefix).assembly_summary_refseq.txt")

            today_genbank_cache = write_ncbi_metadata_fixture(
                joinpath(cache_dir, "$(today_prefix).assembly_summary_genbank.txt");
                assembly_accession = "GCA_TODAY_GENBANK",
                ftp_path = "ftp://example/GCA_TODAY_GENBANK",
                taxid = "222"
            )

            today_df = Test.@test_logs (:info, r"Found valid cached file for today") min_level=Logging.Info match_mode=:any begin
                Mycelia.load_ncbi_metadata("genbank")
            end
            Test.@test today_df[1, "assembly_accession"] == "GCA_TODAY_GENBANK"
            Test.@test today_df[1, "taxid"] == 222
            Test.@test today_genbank_cache == joinpath(cache_dir, "$(today_prefix).assembly_summary_genbank.txt")

            refseq_wrapper = Mycelia.load_refseq_metadata(date = fixture_date)
            Test.@test refseq_wrapper[1, "assembly_accession"] == "GCF_REQUESTED_REFSEQ"

            requested_genbank_cache = write_ncbi_metadata_fixture(
                joinpath(cache_dir, "$(date_prefix).assembly_summary_genbank.txt");
                assembly_accession = "GCA_REQUESTED_GENBANK",
                ftp_path = "ftp://example/GCA_REQUESTED_GENBANK",
                taxid = "333"
            )

            genbank_wrapper = Mycelia.load_genbank_metadata(date = fixture_date)
            Test.@test genbank_wrapper[1, "assembly_accession"] == "GCA_REQUESTED_GENBANK"
            Test.@test requested_genbank_cache == joinpath(cache_dir, "$(date_prefix).assembly_summary_genbank.txt")

            Test.@test_throws ErrorException Mycelia.load_ncbi_metadata(
                "refseq";
                date = Dates.Date(2024, 2, 4)
            )
        end

        with_temp_home() do home
            write(joinpath(home, "workspace"), "occupied by file")
            Test.@test_logs (:error, r"Failed to create cache directory") min_level=Logging.Error match_mode=:any begin
                Test.@test_throws Exception Mycelia.load_ncbi_metadata("refseq"; date = fixture_date)
            end
        end
    end

    Test.@testset "parse_ncbi_metadata_file warning and logged failure paths" begin
        mktempdir() do dir
            warned_path = joinpath(dir, "assembly_summary_extra_header.txt")
            warned_header = [
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
                "gc_percent",
                "assembly_notes"
            ]
            warned_row = [
                "GCF_000002",
                "PRJNA2",
                "SAMN2",
                "WGS2",
                "representative genome",
                "Bacillus subtilis",
                "strain=168",
                "168",
                "latest",
                "Complete Genome",
                "Major",
                "Full",
                "2024-01-01",
                "ASM2",
                "Submitter",
                "none",
                "none",
                "ftp://example/GCF_000002",
                "none",
                "none",
                "none",
                "haploid",
                "bacteria",
                "Provider",
                "Annot2",
                "2024-01-02",
                "67890",
                "1423",
                "1423",
                "4200",
                "4000",
                "1",
                "1",
                "1",
                "3900",
                "3700",
                "200",
                "43.5",
                "manually-curated"
            ]

            open(warned_path, "w") do io
                println(io, "## Assembly summary warning test")
                println(io, "# " * join(warned_header, '\t'))
                println(io, join(warned_row, '\t'))
            end

            warned_df = Test.@test_logs (:warn, r"Headers found in data file") min_level=Logging.Warn match_mode=:any begin
                Mycelia.parse_ncbi_metadata_file(warned_path)
            end
            Test.@test warned_df[1, "assembly_notes"] == "manually-curated"

            invalid_path = joinpath(dir, "assembly_summary_invalid_value.txt")
            invalid_row = copy(warned_row)
            invalid_row[28] = "not_an_int"
            open(invalid_path, "w") do io
                println(io, "## Assembly summary invalid value test")
                println(io, "# " * join(warned_header, '\t'))
                println(io, join(invalid_row, '\t'))
            end

            Test.@test_logs (:error, r"Failed to parse NCBI metadata file") min_level=Logging.Error match_mode=:any begin
                Test.@test_throws Exception Mycelia.parse_ncbi_metadata_file(invalid_path)
            end
        end
    end

    Test.@testset "NCBI FTP path helper" begin
        ftp_path = "ftp://example/GCF_000001"
        url = Mycelia.ncbi_ftp_path_to_url(ftp_path = ftp_path, extension = "genomic.fna.gz")
        Test.@test url == "ftp://example/GCF_000001/GCF_000001_genomic.fna.gz"
    end

    Test.@testset "ncbi_ftp_path_to_url extensions" begin
        ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"
        Test.@test endswith(
            Mycelia.ncbi_ftp_path_to_url(ftp_path = ftp_path, extension = "genomic.gff.gz"),
            "_genomic.gff.gz"
        )
        Test.@test endswith(
            Mycelia.ncbi_ftp_path_to_url(ftp_path = ftp_path, extension = "protein.faa.gz"),
            "_protein.faa.gz"
        )
        Test.@test endswith(
            Mycelia.ncbi_ftp_path_to_url(ftp_path = ftp_path, extension = "assembly_report.txt"),
            "_assembly_report.txt"
        )
    end

    Test.@testset "Constants" begin
        Test.@test Mycelia.NCBI_SUPERKINGDOMS isa Dict
        Test.@test length(Mycelia.NCBI_SUPERKINGDOMS) == 6
        Test.@test Mycelia.NCBI_SUPERKINGDOMS["Bacteria"] == 2
        Test.@test Mycelia.NCBI_SUPERKINGDOMS["Archaea"] == 2157
        Test.@test Mycelia.NCBI_SUPERKINGDOMS["Eukaryota"] == 2759
        Test.@test Mycelia.NCBI_SUPERKINGDOMS["Viruses"] == 10239

        Test.@test Mycelia.UN_LANGUAGES isa Vector{String}
        Test.@test length(Mycelia.UN_LANGUAGES) == 6
        Test.@test "en" in Mycelia.UN_LANGUAGES
        Test.@test "zh" in Mycelia.UN_LANGUAGES

        Test.@test Mycelia.UN_PAIRS isa Vector{String}
        Test.@test length(Mycelia.UN_PAIRS) == 15
        Test.@test "en-fr" in Mycelia.UN_PAIRS
        Test.@test "ar-zh" in Mycelia.UN_PAIRS

        Test.@test Mycelia.UN_CORPUS_BASE_URL isa String
        Test.@test occursin("un.org", Mycelia.UN_CORPUS_BASE_URL)

        Test.@test Mycelia.UN_CORPUS_USER_AGENT == "Mycelia/UNCorpus"
    end

    Test.@testset "_resolve_un_archives" begin
        # Single bitext pair
        archives = Mycelia._resolve_un_archives(["en-fr"], ["txt"])
        Test.@test archives == ["UNv1.0.en-fr.tar.gz"]

        # Reversed pair order should be normalized
        archives = Mycelia._resolve_un_archives(["fr-en"], ["txt"])
        Test.@test archives == ["UNv1.0.en-fr.tar.gz"]

        # 6way subset
        archives = Mycelia._resolve_un_archives(["6way"], ["txt"])
        Test.@test archives == ["UNv1.0.6way.tar.gz"]

        # 6-way alias
        archives = Mycelia._resolve_un_archives(["6-way"], ["txt"])
        Test.@test archives == ["UNv1.0.6way.tar.gz"]

        # TEI subset
        archives = Mycelia._resolve_un_archives(["tei"], ["xml"])
        Test.@test length(archives) == 6
        Test.@test all(a -> occursin("TEI", a), archives)

        # bitext subset generates all 15 pairs
        archives = Mycelia._resolve_un_archives(["bitext"], ["txt"])
        Test.@test length(archives) == 15

        # "all" gets everything
        archives = Mycelia._resolve_un_archives(["all"], ["txt"])
        Test.@test length(archives) > 15  # bitext + 6way + TEI

        # Error: empty subsets
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(String[], ["txt"])

        # Error: empty formats
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["bitext"], String[])

        # Error: effectively empty formats after stripping
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["bitext"], [" ", "\t"])

        # Error: effectively empty subsets after stripping
        Test.@test_throws ErrorException Mycelia._resolve_un_archives([" ", "\t"], ["txt"])

        # Error: unsupported format
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["bitext"], ["pdf"])

        # Error: unsupported subset
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["invalid"], ["txt"])

        # Error: unsupported language
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["en-de"], ["txt"])

        # Error: language codes are valid but pair itself is unsupported
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["en-en"], ["txt"])

        # Error: bitext without txt format
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["bitext"], ["xml"])

        # Error: tei without xml format
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["tei"], ["txt"])
    end

    Test.@testset "_merge_un_archive_parts" begin
        mktempdir() do dir
            # Create small test parts
            part1 = joinpath(dir, "part.00")
            part2 = joinpath(dir, "part.01")
            merged = joinpath(dir, "merged.dat")

            write(part1, "Hello, ")
            write(part2, "World!")

            Mycelia._merge_un_archive_parts([part1, part2], merged)
            Test.@test isfile(merged)
            Test.@test read(merged, String) == "Hello, World!"

            # Single part
            merged2 = joinpath(dir, "merged2.dat")
            Mycelia._merge_un_archive_parts([part1], merged2)
            Test.@test read(merged2, String) == "Hello, "

            # Empty parts list
            merged3 = joinpath(dir, "merged3.dat")
            Mycelia._merge_un_archive_parts(String[], merged3)
            Test.@test isfile(merged3)
            Test.@test filesize(merged3) == 0
        end
    end

    Test.@testset "collect_un_language_files" begin
        mktempdir() do dir
            # Create mock UN corpus structure
            mkpath(joinpath(dir, "sub"))
            write(joinpath(dir, "UNv1.0.en-fr.en"), "english text")
            write(joinpath(dir, "UNv1.0.en-fr.fr"), "french text")
            write(joinpath(dir, "UNv1.0.6way.en"), "6way english")
            write(joinpath(dir, "sub", "UNv1.0.ar-en.ar"), "arabic text")
            # Non-matching files
            write(joinpath(dir, "README.md"), "readme")
            write(joinpath(dir, "data.tar.gz"), "archive")
            write(joinpath(dir, "noversion.en"), "no UNv1.0 prefix")

            # all subsets
            result = Mycelia.collect_un_language_files(dir; subsets = ["all"])
            Test.@test haskey(result, "en")
            Test.@test haskey(result, "fr")
            Test.@test haskey(result, "ar")
            Test.@test length(result["en"]) == 2  # en-fr.en, 6way.en

            # specific subset filter
            result2 = Mycelia.collect_un_language_files(dir; subsets = ["en-fr"])
            Test.@test haskey(result2, "en")
            Test.@test haskey(result2, "fr")
            Test.@test !haskey(result2, "ar")

            # strict_naming=false allows files without UNv1.0 prefix
            write(joinpath(dir, "testset.en"), "test english")
            result3 = Mycelia.collect_un_language_files(dir; subsets = ["all"], strict_naming = false)
            # Should include testset.en and noversion.en now
            en_files = result3["en"]
            Test.@test length(en_files) > length(result["en"])
        end
    end

    Test.@testset "collect_un_language_files skips malformed filenames" begin
        mktempdir() do dir
            candidate_path = joinpath(dir, "UNv1.0.testset.en")
            valid_path = joinpath(dir, "UNv1.0.6way.en")
            write(joinpath(dir, "README"), "missing language separator")
            write(candidate_path, "missing pair marker")
            write(joinpath(dir, "UNv1.0.6way"), "missing language suffix")
            write(joinpath(dir, "UNv1.0.6way.engl"), "language suffix too long")
            write(joinpath(dir, "UNv1.0.6way.e1"), "non-letter language suffix")
            write(valid_path, "valid english file")

            strict_result = Mycelia.collect_un_language_files(dir; subsets = ["all"])
            Test.@test haskey(strict_result, "en")
            Test.@test strict_result["en"] == [valid_path]

            loose_result = Mycelia.collect_un_language_files(dir; subsets = ["all"], strict_naming = false)
            Test.@test haskey(loose_result, "en")
            Test.@test Set(loose_result["en"]) == Set([candidate_path, valid_path])
            Test.@test !haskey(loose_result, "engl")
            Test.@test !haskey(loose_result, "e1")

            filtered_result = Mycelia.collect_un_language_files(dir; subsets = ["6way"], strict_naming = false)
            Test.@test haskey(filtered_result, "en")
            Test.@test filtered_result["en"] == [valid_path]
        end
    end

    Test.@testset "download_sra_data error paths" begin
        Test.@test_throws ErrorException Mycelia.download_sra_data("")
    end

    Test.@testset "download_sra_data path resolution" begin
        requested_outdir = "/tmp/sra-downloads"
        resolved_outdir = joinpath(requested_outdir, "SRR1234567")

        paired_result = Mycelia._download_sra_data_result(
            "SRR1234567",
            (
                forward_reads = joinpath(resolved_outdir, "SRR1234567_1.fastq.gz"),
                reverse_reads = joinpath(resolved_outdir, "SRR1234567_2.fastq.gz"),
                unpaired_reads = missing
            )
        )
        Test.@test paired_result.outdir == resolved_outdir
        Test.@test paired_result.files == [
            joinpath(resolved_outdir, "SRR1234567_1.fastq.gz"),
            joinpath(resolved_outdir, "SRR1234567_2.fastq.gz")
        ]
        Test.@test paired_result.is_paired

        single_result = Mycelia._download_sra_data_result(
            "SRR7654321",
            (
                forward_reads = missing,
                reverse_reads = missing,
                unpaired_reads = joinpath(requested_outdir, "SRR7654321", "SRR7654321.fastq.gz")
            )
        )
        Test.@test single_result.outdir == joinpath(requested_outdir, "SRR7654321")
        Test.@test single_result.files == [joinpath(requested_outdir, "SRR7654321", "SRR7654321.fastq.gz")]
        Test.@test !single_result.is_paired

        Test.@test_throws ErrorException Mycelia._download_sra_data_result(
            "SRR000001",
            (
                forward_reads = joinpath(resolved_outdir, "SRR000001_1.fastq.gz"),
                reverse_reads = missing,
                unpaired_reads = missing
            )
        )
    end

    install_noop_add_bioconda_env!()
    try

    if RUN_EXTERNAL
        Test.@testset "Cached SRA helper coverage is core-only" begin
            Test.@test_skip "Cached prefetch/fasterq coverage is skipped when external tooling is enabled"
        end

        Test.@testset "prefetch and fasterq_dump cached outputs" begin
            mktempdir() do dir
                paired_srr = "SRR_PAIRED_CACHED"
                paired_dir = joinpath(dir, paired_srr)
                mkpath(paired_dir)
                paired_archive = joinpath(paired_dir, "$(paired_srr).sra")
                write(paired_archive, "cached archive")

                paired_prefetch = Test.@test_logs (:info, r"SRA archive already present") begin
                    Mycelia.prefetch(SRR = paired_srr, outdir = dir)
                end
                Test.@test paired_prefetch.directory == paired_dir
                Test.@test paired_prefetch.archive == paired_archive

                paired_forward = joinpath(paired_dir, "$(paired_srr)_1.fastq.gz")
                paired_reverse = joinpath(paired_dir, "$(paired_srr)_2.fastq.gz")
                write(paired_forward, "forward")
                write(paired_reverse, "reverse")

                paired_result = Mycelia.fasterq_dump(outdir = dir, srr_identifier = paired_srr)
                Test.@test paired_result.forward_reads == paired_forward
                Test.@test paired_result.reverse_reads == paired_reverse
                Test.@test ismissing(paired_result.unpaired_reads)

                single_srr = "SRR_SINGLE_CACHED"
                single_dir = joinpath(dir, single_srr)
                mkpath(single_dir)
                write(joinpath(single_dir, "$(single_srr).sra"), "cached archive")
                single_unpaired = joinpath(single_dir, "$(single_srr).fastq.gz")
                write(single_unpaired, "single")

                single_result = Mycelia.fasterq_dump(outdir = dir, srr_identifier = single_srr)
                Test.@test ismissing(single_result.forward_reads)
                Test.@test ismissing(single_result.reverse_reads)
                Test.@test single_result.unpaired_reads == single_unpaired
            end
        end

        Test.@testset "prefetch and fasterq_dump command branches" begin
            with_fake_sra_tooling() do
                mktempdir() do dir
                    prefetched = Mycelia.prefetch(SRR = "SRR_PREFETCH_LIVE", outdir = dir)
                    Test.@test prefetched.directory == joinpath(dir, "SRR_PREFETCH_LIVE")
                    Test.@test isfile(prefetched.archive)

                    paired_srr = "SRR_PAIRED_LIVE"
                    paired_result = Mycelia.fasterq_dump(outdir = dir, srr_identifier = paired_srr)
                    Test.@test paired_result.forward_reads == joinpath(dir, paired_srr, "$(paired_srr)_1.fastq.gz")
                    Test.@test paired_result.reverse_reads == joinpath(dir, paired_srr, "$(paired_srr)_2.fastq.gz")
                    Test.@test ismissing(paired_result.unpaired_reads)
                    Test.@test isfile(paired_result.forward_reads)
                    Test.@test isfile(paired_result.reverse_reads)

                    single_srr = "SRR_SINGLE_LIVE"
                    single_result = Mycelia.fasterq_dump(outdir = dir, srr_identifier = single_srr)
                    Test.@test ismissing(single_result.forward_reads)
                    Test.@test ismissing(single_result.reverse_reads)
                    Test.@test single_result.unpaired_reads == joinpath(dir, single_srr, "$(single_srr).fastq.gz")
                    Test.@test isfile(single_result.unpaired_reads)
                end
            end
        end

        install_synthetic_sra_helpers!()
        try

            Test.@testset "download_sra_data wrapper branches" begin
                mktempdir() do dir
                    paired = Mycelia.download_sra_data("SRR_PAIRED_WRAPPER"; outdir = dir)
                    Test.@test paired.srr_id == "SRR_PAIRED_WRAPPER"
                    Test.@test paired.is_paired
                    Test.@test length(paired.files) == 2
                    Test.@test all(isfile, paired.files)

                    single = Mycelia.download_sra_data("SRR_SINGLE_WRAPPER"; outdir = dir)
                    Test.@test single.srr_id == "SRR_SINGLE_WRAPPER"
                    Test.@test !single.is_paired
                    Test.@test length(single.files) == 1
                    Test.@test isfile(only(single.files))

                    Test.@test_throws ErrorException Mycelia.download_sra_data("SRR_MISSING_WRAPPER"; outdir = dir)
                end
            end

            Test.@testset "prefetch_sra_runs error path" begin
                Test.@test_throws ErrorException Mycelia.prefetch_sra_runs(String[])
            end

            Test.@testset "prefetch_sra_runs collects success and failure results" begin
                mktempdir() do dir
                    results = Mycelia.prefetch_sra_runs(
                        ["SRR_BATCH_A", "SRR_BATCH_FAIL", "SRR_BATCH_B"];
                        outdir = dir,
                        max_parallel = 2
                    )
                    Test.@test length(results) == 3
                    successes = Dict(result.srr_id => result for result in results)
                    Test.@test successes["SRR_BATCH_A"].success
                    Test.@test successes["SRR_BATCH_B"].success
                    Test.@test !successes["SRR_BATCH_FAIL"].success
                    Test.@test occursin("synthetic prefetch failure", successes["SRR_BATCH_FAIL"].error)
                end
            end

            Test.@testset "fasterq_dump_parallel error path" begin
                Test.@test_throws ErrorException Mycelia.fasterq_dump_parallel(String[])
            end

            Test.@testset "fasterq_dump_parallel collects success and failure results" begin
                mktempdir() do dir
                    results = Mycelia.fasterq_dump_parallel(
                        ["SRR_PAIRED_A", "SRR_PAIRED_FAIL", "SRR_SINGLE_B"];
                        outdir = dir,
                        max_parallel = 2
                    )
                    Test.@test length(results) == 3
                    outcomes = Dict(result.srr_id => result for result in results)
                    Test.@test outcomes["SRR_PAIRED_A"].success
                    Test.@test outcomes["SRR_SINGLE_B"].success
                    Test.@test !outcomes["SRR_PAIRED_FAIL"].success
                    Test.@test occursin("synthetic fasterq failure", outcomes["SRR_PAIRED_FAIL"].error)
                end
            end
        finally
            restore_sra_helpers!()
        end
    finally
        restore_add_bioconda_env!()
    end

    Test.@testset "BLAST database helper coverage" begin
        if !RUN_EXTERNAL
            @info "Skipping BLAST helper coverage tests; external tool execution is opt-in via MYCELIA_RUN_EXTERNAL=true"
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true to run BLAST helper coverage tests"
        elseif !blast_env_available()
            @info "Skipping BLAST helper coverage tests; blast env unavailable"
            Test.@test_skip "BLAST conda env unavailable"
        else
            mktempdir() do dir
                nucleotide_db = make_test_blastdb(dir; dbtype = "nucl")
                protein_db = make_test_blastdb(dir; dbtype = "prot")

                nucleotide_metadata = Mycelia.get_blastdb_metadata(blastdb = nucleotide_db)
                Test.@test nucleotide_metadata["dbtype"] == "Nucleotide"
                Test.@test nucleotide_metadata["number-of-sequences"] == 2
                Test.@test Dates.DateTime(nucleotide_metadata["last-updated"]) isa Dates.DateTime

                protein_metadata = Mycelia.get_blastdb_metadata(blastdb = protein_db)
                Test.@test protein_metadata["dbtype"] == "Protein"
                Test.@test protein_metadata["number-of-sequences"] == 2

                Test.@test begin
                    Mycelia.get_blastdb_info(blastdb = nucleotide_db)
                    true
                end

                tax_info = Mycelia.get_blastdb_tax_info(blastdb = nucleotide_db)
                Test.@test Mycelia.DataFrames.nrow(tax_info) == 1
                Test.@test tax_info[1, "num of seqs"] == 2

                tax_info_from_entries = Mycelia.get_blastdb_tax_info(
                    blastdb = nucleotide_db,
                    entries = ["seq1"]
                )
                Test.@test Mycelia.DataFrames.nrow(tax_info_from_entries) == 1

                tax_info_from_taxids = Mycelia.get_blastdb_tax_info(
                    blastdb = nucleotide_db,
                    taxids = [0]
                )
                Test.@test Mycelia.DataFrames.nrow(tax_info_from_taxids) == 1

                Test.@test_throws ErrorException Mycelia.get_blastdb_tax_info(
                    blastdb = nucleotide_db,
                    entries = ["seq1"],
                    taxids = [0]
                )

                nucleotide_date = Dates.format(
                    Dates.DateTime(nucleotide_metadata["last-updated"]),
                    "yyyy-mm-dd"
                )
                nucleotide_outfile = "$(nucleotide_db).$(nucleotide_date).fna.gz"
                write(nucleotide_outfile, "existing nucleotide archive")
                Test.@test Mycelia.blastdb_to_fasta(
                    blastdb = nucleotide_db,
                    force = false
                ) == nucleotide_outfile

                protein_date = Dates.format(
                    Dates.DateTime(protein_metadata["last-updated"]),
                    "yyyy-mm-dd"
                )
                protein_outfile = "$(protein_db).$(protein_date).faa.gz"
                write(protein_outfile, "existing protein archive")
                Test.@test Mycelia.blastdb_to_fasta(
                    blastdb = protein_db,
                    force = false
                ) == protein_outfile

                explicit_outfile = joinpath(dir, "already-there.fna.gz")
                write(explicit_outfile, "no-op")
                Test.@test Mycelia.blastdb_to_fasta(
                    blastdb = nucleotide_db,
                    outfile = explicit_outfile,
                    force = false
                ) == explicit_outfile

                exported_outfile = joinpath(dir, "exports", "full-db.fna.gz")
                exported_path = Mycelia.blastdb_to_fasta(
                    blastdb = nucleotide_db,
                    outfile = exported_outfile,
                    force = true,
                    max_cores = 1
                )
                Test.@test exported_path == exported_outfile
                Test.@test isfile(exported_path)
                exported_contents = read(`gzip -dc $(exported_path)`, String)
                Test.@test occursin(">seq1", exported_contents)
                Test.@test occursin("ACGTACGT", exported_contents)

                entry_subset_outfile = joinpath(dir, "exports", "entry-subset.fna.gz")
                entry_subset_path = Mycelia.blastdb_to_fasta(
                    blastdb = nucleotide_db,
                    entries = ["seq1"],
                    outfile = entry_subset_outfile,
                    force = true,
                    max_cores = 1
                )
                Test.@test entry_subset_path == entry_subset_outfile
                Test.@test isfile(entry_subset_path)
                entry_subset_contents = read(`gzip -dc $(entry_subset_path)`, String)
                Test.@test occursin(">seq1", entry_subset_contents)
                Test.@test !occursin(">seq2", entry_subset_contents)

                taxid_subset_outfile = joinpath(dir, "exports", "taxid-subset.fna.gz")
                taxid_subset_path = Mycelia.blastdb_to_fasta(
                    blastdb = nucleotide_db,
                    taxids = [0],
                    outfile = taxid_subset_outfile,
                    force = true,
                    max_cores = 1
                )
                Test.@test taxid_subset_path == taxid_subset_outfile
                Test.@test isfile(taxid_subset_path)
                Test.@test filesize(taxid_subset_path) > 0
                taxid_subset_contents = read(`gzip -dc $(taxid_subset_path)`, String)
                Test.@test occursin(">seq1", taxid_subset_contents)

                Test.@test_throws ErrorException Mycelia.blastdb_to_fasta(
                    blastdb = nucleotide_db,
                    entries = ["seq1"],
                    taxids = [0]
                )
            end
        end
    end

    Test.@testset "UN download helper coverage" begin
        split_archive = make_test_tarball_bytes(Dict("split/archive.txt" => "split archive contents"))
        direct_archive = make_test_tarball_bytes(Dict("direct/archive.txt" => "direct archive contents"))
        testset_archive = make_test_tarball_bytes(Dict("testsets/sample.en" => "testset contents"))
        split_archive_parts = split_bytes(split_archive)

        server_files = Dict(
            "/download.txt" => Vector{UInt8}(codeunits("downloaded via helper")),
            "/range-fallback" => Vector{UInt8}(codeunits("range fallback")),
            "/forbidden-fallback" => Vector{UInt8}(codeunits("forbidden fallback")),
            "/UNv1.0.en-fr.tar.gz.00" => split_archive_parts[1],
            "/UNv1.0.en-fr.tar.gz.01" => split_archive_parts[2],
            "/direct.tar.gz" => direct_archive,
            "/UNv1.0.testsets.tar.gz" => testset_archive
        )

        server = start_un_test_server(server_files)
        try
            redirected_response = Mycelia._un_request("$(server.base_url)/redirect-start", "GET")
            Test.@test redirected_response.status == 200
            Test.@test String(redirected_response.body) == "redirected"

            Test.@test_throws ArgumentError Mycelia._un_request(
                "$(server.base_url)/redirect-no-location",
                "GET"
            )
            Test.@test_throws ErrorException Mycelia._un_request(
                "$(server.base_url)/redirect-loop",
                "GET";
                max_redirects = 1
            )

            Test.@test Mycelia._un_part_exists("$(server.base_url)/download.txt")
            Test.@test !Mycelia._un_part_exists("$(server.base_url)/missing.txt")
            Test.@test Mycelia._un_part_exists("$(server.base_url)/range-fallback")
            Test.@test !Mycelia._un_part_exists("$(server.base_url)/range-missing-fallback")
            Test.@test Mycelia._un_part_exists("$(server.base_url)/forbidden-fallback")
            Test.@test_throws ErrorException Mycelia._un_part_exists("$(server.base_url)/unexpected-status")

            mktempdir() do dir
                downloaded_path = joinpath(dir, "download.txt")
                Test.@test Mycelia._un_download(
                    "$(server.base_url)/download.txt",
                    downloaded_path
                ) == downloaded_path
                Test.@test read(downloaded_path, String) == "downloaded via helper"

                cached_archive_dir = joinpath(dir, "cached")
                mkpath(cached_archive_dir)
                cached_archive_path = joinpath(cached_archive_dir, "cached.tar.gz")
                write(cached_archive_path, direct_archive)
                cached_result = Mycelia._download_and_extract_split_archive(
                    "cached.tar.gz",
                    server.base_url,
                    cached_archive_dir
                )
                Test.@test cached_result == cached_archive_dir
                Test.@test read(
                    joinpath(cached_archive_dir, "direct", "archive.txt"),
                    String
                ) == "direct archive contents"
                Test.@test isfile(cached_archive_path)

                split_outdir = joinpath(dir, "split")
                split_result = Mycelia._download_and_extract_split_archive(
                    "UNv1.0.en-fr.tar.gz",
                    server.base_url,
                    split_outdir
                )
                Test.@test split_result == split_outdir
                Test.@test read(joinpath(split_outdir, "split", "archive.txt"), String) ==
                    "split archive contents"
                Test.@test !isfile(joinpath(split_outdir, "UNv1.0.en-fr.tar.gz"))
                Test.@test !isfile(joinpath(split_outdir, "UNv1.0.en-fr.tar.gz.00"))
                Test.@test !isfile(joinpath(split_outdir, "UNv1.0.en-fr.tar.gz.01"))

                cached_split_outdir = joinpath(dir, "cached-split")
                mkpath(cached_split_outdir)
                write(joinpath(cached_split_outdir, "UNv1.0.en-fr.tar.gz.00"), split_archive_parts[1])
                write(joinpath(cached_split_outdir, "UNv1.0.en-fr.tar.gz.01"), split_archive_parts[2])
                cached_split_result = Test.@test_logs (:info, r"Using existing part") min_level=Logging.Info match_mode=:any begin
                    Mycelia._download_and_extract_split_archive(
                        "UNv1.0.en-fr.tar.gz",
                        server.base_url,
                        cached_split_outdir
                    )
                end
                Test.@test cached_split_result == cached_split_outdir
                Test.@test read(joinpath(cached_split_outdir, "split", "archive.txt"), String) ==
                    "split archive contents"
                Test.@test !isfile(joinpath(cached_split_outdir, "UNv1.0.en-fr.tar.gz.00"))
                Test.@test !isfile(joinpath(cached_split_outdir, "UNv1.0.en-fr.tar.gz.01"))

                direct_outdir = joinpath(dir, "direct")
                direct_result = Mycelia._download_and_extract_split_archive(
                    "direct.tar.gz",
                    server.base_url,
                    direct_outdir
                )
                Test.@test direct_result == direct_outdir
                Test.@test read(joinpath(direct_outdir, "direct", "archive.txt"), String) ==
                    "direct archive contents"
                Test.@test !isfile(joinpath(direct_outdir, "direct.tar.gz"))

                Test.@test_throws ErrorException Mycelia._download_and_extract_split_archive(
                    "missing.tar.gz",
                    server.base_url,
                    joinpath(dir, "missing")
                )

                corpus_dir = Mycelia.download_un_parallel_corpus(
                    outdir = dir,
                    subsets = ["en-fr"],
                    base_url = "$(server.base_url)/"
                )
                Test.@test corpus_dir == joinpath(dir, "un_corpus")
                Test.@test read(joinpath(corpus_dir, "split", "archive.txt"), String) ==
                    "split archive contents"

                testset_dir = Mycelia.download_un_parallel_corpus_testset(
                    outdir = dir,
                    base_url = "$(server.base_url)/"
                )
                Test.@test testset_dir == joinpath(dir, "un_corpus_testsets")
                Test.@test read(joinpath(testset_dir, "testsets", "sample.en"), String) ==
                    "testset contents"
            end
        finally
            Mycelia.HTTP.forceclose(server.server)
        end
    end

    Test.@testset "_resolve_un_archives empty-selection error path" begin
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["6way"], ["xml"])
    end
end
