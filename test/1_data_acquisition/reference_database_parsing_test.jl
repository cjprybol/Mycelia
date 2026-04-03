import Test
import Dates
import Logging
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

const RUN_ALL = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true"
const RUN_EXTERNAL = RUN_ALL || lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"

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
end
