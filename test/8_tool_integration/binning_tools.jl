"""
Binning tool wrappers and parsers.

Lightweight parser and validation tests run by default. Integration tests that
invoke external tools are opt-in via `MYCELIA_RUN_EXTERNAL=true`.

# Example:
#   julia --project=. -e 'include("test/8_tool_integration/binning_tools.jl")'
#   MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/binning_tools.jl")'
"""

import DataFrames
import Test
import Mycelia

const RUN_ALL = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
const RUN_EXTERNAL = RUN_ALL || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

Test.@testset "Binning Tools Integration" begin

    Test.@testset "Parser utilities" begin
        # Generic contig/bin table
        tmp = tempname()
        open(tmp, "w") do io
            write(io, "contig\tbin\tlength\n")
            write(io, "ctg1\tbin1\t1500\n")
            write(io, "ctg2\tbin2\t2100\n")
            write(io, "ctg3\tbin1\t900\n")
        end

        df = Mycelia.parse_bin_assignments(tmp)
        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) == 3
        Test.@test df.bin[1] == "bin1"
        Test.@test df.contig[3] == "ctg3"

        # Custom column names
        tmp_custom = tempname()
        open(tmp_custom, "w") do io
            write(io, "contig_id\tcluster\tcov\n")
            write(io, "aaa\tBinA\t3.2\n")
            write(io, "bbb\tBinB\t4.1\n")
        end

        df_custom = Mycelia.parse_bin_assignments(tmp_custom; contig_col="contig_id", bin_col="cluster")
        Test.@test DataFrames.nrow(df_custom) == 2
        Test.@test df_custom.bin[2] == "BinB"
        Test.@test df_custom.contig[1] == "aaa"

        # dRep cluster table
        tmp_drep = tempname()
        open(tmp_drep, "w") do io
            write(io, "genome,secondary_cluster,representative,ani\n")
            write(io, "g1,1,g1,0.99\n")
            write(io, "g2,1,g1,0.99\n")
            write(io, "g3,2,g3,0.98\n")
        end

        df_drep = Mycelia.parse_drep_clusters(tmp_drep)
        Test.@test DataFrames.nrow(df_drep) == 3
        Test.@test df_drep.secondary_cluster[1] == "1"
        Test.@test df_drep.representative[2] == "g1"
    end

    Test.@testset "Input validation" begin
        outdir = mktempdir()
        Test.@test_throws ErrorException Mycelia.run_vamb(
            contigs_fasta="missing_contigs.fna",
            depth_file="missing_depth.tsv",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_metabat2(
            contigs_fasta="missing_contigs.fna",
            depth_file="missing_depth.tsv",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_metacoag(
            contigs_fasta="missing_contigs.fna",
            assembly_graph="missing.gfa",
            mapping_file="missing.tsv",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_comebin(
            contigs_fasta="missing_contigs.fna",
            coverage_table="missing_cov.tsv",
            marker_file="missing_markers.tsv",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_drep_dereplicate(
            genomes=String[],
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_taxometer(
            contigs_fasta="missing_contigs.fna",
            depth_file="missing_depth.tsv",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_taxvamb(
            contigs_fasta="missing_contigs.fna",
            depth_file="missing_depth.tsv",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_genomeface(
            contigs_fasta="missing_contigs.fna",
            coverage_table="missing_cov.tsv",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_magmax_merge(
            bins_dirs=String[],
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_magmax_merge(
            bins_dirs=["missing_bins_dir"],
            outdir=outdir
        )
    end

    if RUN_EXTERNAL
        @info "No external binning tools executed in this suite; add commands when tools are available."
    else
        @info "Skipping binning external execution; set MYCELIA_RUN_EXTERNAL=true to enable when tools are installed."
    end
end
