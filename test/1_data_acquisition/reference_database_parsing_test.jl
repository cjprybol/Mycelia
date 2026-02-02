import Test
import Mycelia

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

    Test.@testset "NCBI FTP path helper" begin
        ftp_path = "ftp://example/GCF_000001"
        url = Mycelia.ncbi_ftp_path_to_url(ftp_path = ftp_path, extension = "genomic.fna.gz")
        Test.@test url == "ftp://example/GCF_000001/GCF_000001_genomic.fna.gz"
    end
end
