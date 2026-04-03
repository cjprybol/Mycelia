import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function write_fake_conda_runner(path::AbstractString, capture_dir::AbstractString)
    script = """
#!/usr/bin/env bash
set -euo pipefail

capture_dir="$capture_dir"

if [[ "\${1:-}" == "run" ]]; then
    shift
fi

if [[ "\${1:-}" == "--live-stream" ]]; then
    shift
fi

if [[ "\${1:-}" == "-n" ]]; then
    env_name="\${2:-}"
    shift 2
fi

tool="\${1:-}"
if [[ -n "\$tool" ]]; then
    shift
fi

if [[ "\$tool" == "blastdbcmd" ]]; then
    db=""
    batch_file=""
    taxid_file=""
    metadata=0
    info=0
    tax_info=0

    while [[ \$# -gt 0 ]]; do
        case "\$1" in
            -db)
                db="\${2:-}"
                shift 2
                ;;
            -entry)
                shift 2
                ;;
            -entry_batch)
                batch_file="\${2:-}"
                shift 2
                ;;
            -taxidlist)
                taxid_file="\${2:-}"
                shift 2
                ;;
            -metadata)
                metadata=1
                shift
                ;;
            -info)
                info=1
                shift
                ;;
            -tax_info)
                tax_info=1
                shift
                ;;
            -outfmt)
                shift 2
                ;;
            *)
                shift
                ;;
        esac
    done

    if [[ \$metadata -eq 1 ]]; then
        case "\$db" in
            nt)
                printf '{"dbtype":"Nucleotide","last-updated":"2024-01-02T03:04:05"}'
                ;;
            nr)
                printf '{"dbtype":"Protein","last-updated":"2024-02-03T04:05:06"}'
                ;;
            weird)
                printf '{"dbtype":"Unsupported","last-updated":"2024-05-06T07:08:09"}'
                ;;
            *)
                printf '{"dbtype":"Protein","last-updated":"2024-03-05T12:00:00"}'
                ;;
        esac
        exit 0
    fi

    if [[ \$info -eq 1 ]]; then
        printf 'Database: %s\\n' "\$db"
        exit 0
    fi

    if [[ \$tax_info -eq 1 ]]; then
        printf '# comment\\n'
        printf '2\\tBacillus subtilis\\tbacterium\\tBacteria\\tBacteria\\t1\\n'
        exit 0
    fi

    if [[ -n "\$batch_file" ]]; then
        printf '%s\\n' "\$batch_file" > "\$capture_dir/entries_path.txt"
        cp "\$batch_file" "\$capture_dir/entries_capture.txt"
        while IFS= read -r entry; do
            printf '>%s\\nSEQUENCE\\n' "\$entry"
        done < "\$batch_file"
        exit 0
    fi

    if [[ -n "\$taxid_file" ]]; then
        printf '%s\\n' "\$taxid_file" > "\$capture_dir/taxids_path.txt"
        cp "\$taxid_file" "\$capture_dir/taxids_capture.txt"
        while IFS= read -r taxid; do
            printf '>taxid:%s\\nSEQUENCE\\n' "\$taxid"
        done < "\$taxid_file"
        exit 0
    fi

    printf '>all\\nSEQUENCE\\n'
    exit 0
fi

if [[ "\$tool" == "pigz" ]]; then
    cat
    exit 0
fi

printf 'unexpected tool: %s\\n' "\$tool" >&2
exit 1
"""

    write(path, script)
    chmod(path, 0o755)
    return path
end

function with_fake_conda_runner(f::Function)
    mktempdir() do dir
        # Use PATH-based isolation instead of overwriting Mycelia.CONDA_RUNNER.
        # Create a fake conda script in the temp dir and prepend to PATH so
        # the module finds our fake before the real one.
        fake_conda = write_fake_conda_runner(joinpath(dir, "conda"), dir)
        original_path = get(ENV, "PATH", "")
        ENV["PATH"] = dir * ":" * original_path
        try
            return f(dir)
        finally
            ENV["PATH"] = original_path
        end
    end
end

Test.@testset "Reference Database Parsing" begin
    Test.@testset "BLAST database helpers" begin
        with_fake_conda_runner() do capture_dir
            Test.@test_throws ErrorException Mycelia.blastdb_to_fasta(
                blastdb = "nr",
                entries = ["seq1"],
                taxids = [2]
            )

            mktempdir() do dir
                cd(dir) do
                    write("nt.2024-01-02.fna.gz", "cached nucleotide")
                    write("nr.2024-02-03.faa.gz", "cached protein")

                    Test.@test Mycelia.blastdb_to_fasta(
                        blastdb = "nt",
                        force = false
                    ) == "nt.2024-01-02.fna.gz"

                    Test.@test Mycelia.blastdb_to_fasta(
                        blastdb = "nr",
                        force = false
                    ) == "nr.2024-02-03.faa.gz"

                    Test.@test_throws ErrorException Mycelia.blastdb_to_fasta(
                        blastdb = "weird",
                        force = false
                    )
                end

                entries_outfile = joinpath(dir, "nested", "entries.faa.gz")
                Test.@test Mycelia.blastdb_to_fasta(
                    blastdb = "nr",
                    entries = ["seqA", "seqB"],
                    outfile = entries_outfile,
                    max_cores = 0
                ) == entries_outfile
                Test.@test read(entries_outfile, String) ==
                           ">seqA\nSEQUENCE\n>seqB\nSEQUENCE\n"
                Test.@test read(joinpath(capture_dir, "entries_capture.txt"), String) ==
                           "seqA\nseqB\n"
                Test.@test !isfile(strip(read(joinpath(capture_dir, "entries_path.txt"), String)))

                taxids_outfile = joinpath(dir, "taxids.faa.gz")
                Test.@test Mycelia.blastdb_to_fasta(
                    blastdb = "nr",
                    taxids = [2, 2157],
                    outfile = taxids_outfile,
                    max_cores = typemax(Int)
                ) == taxids_outfile
                Test.@test read(taxids_outfile, String) ==
                           ">taxid:2\nSEQUENCE\n>taxid:2157\nSEQUENCE\n"
                Test.@test read(joinpath(capture_dir, "taxids_capture.txt"), String) ==
                           "2\n2157\n"
                Test.@test !isfile(strip(read(joinpath(capture_dir, "taxids_path.txt"), String)))
            end

            blast_metadata = Mycelia.get_blastdb_metadata(blastdb = "nr")
            Test.@test blast_metadata["dbtype"] == "Protein"
            Test.@test blast_metadata["last-updated"] == "2024-02-03T04:05:06"

            Test.@test success(Mycelia.get_blastdb_info(blastdb = "nr"))

            tax_info = Mycelia.get_blastdb_tax_info(blastdb = "nr")
            Test.@test Mycelia.DataFrames.nrow(tax_info) == 1
            Test.@test tax_info[1, "taxid"] == 2
            Test.@test tax_info[1, "scientific name"] == "Bacillus subtilis"
            Test.@test tax_info[1, "taxonomic super kingdom"] == "Bacteria"
        end
    end

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

        # Error: unsupported format
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["bitext"], ["pdf"])

        # Error: unsupported subset
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["invalid"], ["txt"])

        # Error: unsupported language
        Test.@test_throws ErrorException Mycelia._resolve_un_archives(["en-de"], ["txt"])

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

    Test.@testset "download_sra_data error paths" begin
        Test.@test_throws ErrorException Mycelia.download_sra_data("")
    end

    Test.@testset "prefetch_sra_runs error path" begin
        Test.@test_throws ErrorException Mycelia.prefetch_sra_runs(String[])
    end

    Test.@testset "fasterq_dump_parallel error path" begin
        Test.@test_throws ErrorException Mycelia.fasterq_dump_parallel(String[])
    end
end
