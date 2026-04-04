import Test
import Mycelia
import Dates

function write_ncbi_assembly_summary(path::AbstractString; accession::AbstractString, taxid::Int)
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
        accession,
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
        "ftp://example/$(accession)",
        "none",
        "none",
        "none",
        "haploid",
        "bacteria",
        "Provider",
        "Annot1",
        "2023-01-02",
        "12345",
        string(taxid),
        string(taxid),
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

    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "## Assembly summary test")
        println(io, "# " * join(header, '\t'))
        println(io, join(row, '\t'))
    end
    return path
end

function file_url(path::AbstractString)
    return "file://" * abspath(path)
end

Test.@testset "Reference Database Cache Wrappers" begin
    Test.@test_throws ArgumentError Mycelia.load_ncbi_metadata("invalid-db")

    mktempdir() do home_dir
        withenv("HOME" => home_dir) do
            cache_dir = joinpath(home_dir, "workspace", ".ncbi")
            mkpath(cache_dir)

            refseq_date = Dates.Date(2026, 4, 2)
            refseq_cache = joinpath(
                cache_dir,
                "$(Dates.format(refseq_date, "yyyy-mm-dd")).assembly_summary_refseq.txt"
            )
            write_ncbi_assembly_summary(refseq_cache; accession = "GCF_000001", taxid = 562)

            refseq_df = Mycelia.load_ncbi_metadata("refseq"; date = refseq_date)
            Test.@test Mycelia.DataFrames.nrow(refseq_df) == 1
            Test.@test refseq_df[1, "assembly_accession"] == "GCF_000001"

            refseq_wrapper_df = Mycelia.load_refseq_metadata(date = refseq_date)
            Test.@test refseq_wrapper_df[1, "taxid"] == 562

            genbank_date = Dates.Date(2026, 4, 3)
            genbank_cache = joinpath(
                cache_dir,
                "$(Dates.format(genbank_date, "yyyy-mm-dd")).assembly_summary_genbank.txt"
            )
            write_ncbi_assembly_summary(genbank_cache; accession = "GCA_000002", taxid = 9606)

            genbank_df = Mycelia.load_genbank_metadata(date = genbank_date)
            Test.@test Mycelia.DataFrames.nrow(genbank_df) == 1
            Test.@test genbank_df[1, "assembly_accession"] == "GCA_000002"

            today_cache = joinpath(
                cache_dir,
                "$(Dates.format(Dates.today(), "yyyy-mm-dd")).assembly_summary_refseq.txt"
            )
            write_ncbi_assembly_summary(today_cache; accession = "GCF_000003", taxid = 7227)

            today_df = Mycelia.load_ncbi_metadata("refseq")
            Test.@test Mycelia.DataFrames.nrow(today_df) == 1
            Test.@test today_df[1, "assembly_accession"] == "GCF_000003"

            Test.@test_throws ErrorException Mycelia.load_ncbi_metadata(
                "refseq";
                date = Dates.Date(2026, 4, 1)
            )
        end
    end
end

Test.@testset "Reference Database Local Download Wrappers" begin
    mktempdir() do dir
        source_file = joinpath(dir, "source.txt")
        write(source_file, "downloaded fixture")
        copied_file = joinpath(dir, "copied.txt")

        downloaded = Mycelia._un_download(file_url(source_file), copied_file)
        Test.@test downloaded == copied_file
        Test.@test read(downloaded, String) == "downloaded fixture"

        ftp_paths = String[]
        expected_payloads = Set{String}()
        for accession in ["GCF_100001", "GCF_100002"]
            remote_dir = joinpath(dir, accession)
            mkpath(remote_dir)
            remote_file = joinpath(remote_dir, accession * "_genomic.fna.gz")
            write(remote_file, accession)
            push!(ftp_paths, file_url(remote_dir))
            push!(expected_payloads, accession)
        end

        outdir = joinpath(dir, "downloads")
        first_download = Mycelia.download_genome_by_ftp(ftp = ftp_paths[1], outdir = outdir)
        Test.@test basename(first_download) == "GCF_100001_genomic.fna.gz"
        Test.@test read(first_download, String) == "GCF_100001"

        cached_download = Mycelia.download_genome_by_ftp(ftp = ftp_paths[1], outdir = outdir)
        Test.@test cached_download == first_download

        downloads_df = Mycelia.download_genomes_by_ftp(ftp_paths = ftp_paths, outdir = outdir)
        Test.@test Mycelia.DataFrames.nrow(downloads_df) == 2
        Test.@test Set(String.(downloads_df.ftp_path)) == Set(ftp_paths)
        Test.@test Set(read(path, String) for path in downloads_df.fna_path) == expected_payloads
    end
end
