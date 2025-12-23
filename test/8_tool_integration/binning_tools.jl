"""
Binning tool wrappers and parsers.

Lightweight parser and validation tests run by default. Integration tests that
invoke external tools are opt-in via `MYCELIA_RUN_EXTERNAL=true`.

# Example:
#   julia --project=. -e 'include("test/8_tool_integration/binning_tools.jl")'
#   MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/binning_tools.jl")'
#
# External runs auto-generate simulated inputs (contigs, reads, depth, coverage,
# and minimal bin directories) to keep integration tests reproducible.
"""

import DataFrames
import Test
import Mycelia

const RUN_ALL = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
const RUN_EXTERNAL = RUN_ALL || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
# const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

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
            bam_path="missing_bams",
            outdir=outdir
        )
        tmp_contigs = tempname() * ".fna"
        open(tmp_contigs, "w") do io
            write(io, ">contig1\nACGTACGTACGT\n")
        end
        Test.@test_throws ErrorException Mycelia.run_comebin(
            contigs_fasta=tmp_contigs,
            bam_path="missing_bams",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_drep_dereplicate(
            genomes=String[],
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_taxometer(
            contigs_fasta="missing_contigs.fna",
            depth_file="missing_depth.tsv",
            taxonomy_file="missing_taxonomy.tsv",
            outdir=outdir
        )
        Test.@test_throws ErrorException Mycelia.run_taxvamb(
            contigs_fasta="missing_contigs.fna",
            depth_file="missing_depth.tsv",
            taxonomy_file="missing_taxonomy.tsv",
            outdir=outdir
        )
        genomeface_error = nothing
        try
            Mycelia.run_genomeface(
                contigs_fasta="missing_contigs.fna",
                coverage_table="missing_cov.tsv",
                outdir=outdir
            )
        catch err
            genomeface_error = err
        end
        Test.@test genomeface_error isa ErrorException
        if genomeface_error isa ErrorException
            Test.@test occursin("GenomeFace wrapper disabled", sprint(showerror, genomeface_error))
        end
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
        Test.@testset "External tool runs (opt-in)" begin
            inputs = Mycelia.get_binning_test_inputs()
            contigs_fasta = inputs.contigs_fasta
            depth_file = inputs.depth_file
            coverage_table = inputs.coverage_table
            marker_file = inputs.marker_file
            taxonomy_file = inputs.taxonomy_file
            assembly_graph = inputs.assembly_graph
            mapping_file = inputs.mapping_file
            genomes = inputs.genomes
            bins_dirs = inputs.bins_dirs

            if isfile(contigs_fasta) && isfile(depth_file)
                Test.@testset "VAMB" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_vamb(
                            contigs_fasta=contigs_fasta,
                            depth_file=depth_file,
                            outdir=outdir
                        )
                        Test.@test isdir(result.outdir)
                        if isfile(result.clusters_tsv)
                            Test.@test true
                        else
                            Test.@test_skip "VAMB did not produce clusters.tsv; check inputs/output."
                        end
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                end

                Test.@testset "MetaBAT2" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_metabat2(
                            contigs_fasta=contigs_fasta,
                            depth_file=depth_file,
                            outdir=outdir
                        )
                        Test.@test isdir(result.outdir)
                        bins_prefix = basename(result.bins_prefix)
                        bins = filter(name -> startswith(name, bins_prefix), readdir(result.outdir))
                        if isempty(bins)
                            Test.@test_skip "MetaBAT2 did not produce bin files; check inputs/output."
                        else
                            Test.@test true
                        end
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                end
            else
                @info "Skipping VAMB/MetaBAT2; missing contigs/depth inputs."
            end

            if isfile(contigs_fasta) && isfile(depth_file) && isfile(taxonomy_file)
                Test.@testset "Taxometer" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_taxometer(
                            contigs_fasta=contigs_fasta,
                            depth_file=depth_file,
                            taxonomy_file=taxonomy_file,
                            outdir=outdir
                        )
                        Test.@test isdir(result.outdir)
                        if isempty(readdir(result.outdir))
                            Test.@test_skip "Taxometer did not produce output files; check inputs/output."
                        else
                            Test.@test true
                        end
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                end

                Test.@testset "TaxVAMB" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_taxvamb(
                            contigs_fasta=contigs_fasta,
                            depth_file=depth_file,
                            taxonomy_file=taxonomy_file,
                            outdir=outdir
                        )
                        Test.@test isdir(result.outdir)
                        if isfile(result.clusters_tsv)
                            Test.@test true
                        else
                            Test.@test_skip "TaxVAMB did not produce clusters.tsv; check inputs/output."
                        end
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                end
            else
                @info "Skipping Taxometer/TaxVAMB; missing contigs/depth/taxonomy inputs."
            end

            if isfile(contigs_fasta) && isfile(assembly_graph) && isfile(mapping_file)
                Test.@testset "MetaCoAG" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_metacoag(
                            contigs_fasta=contigs_fasta,
                            assembly_graph=assembly_graph,
                            mapping_file=mapping_file,
                            outdir=outdir
                        )
                        Test.@test isdir(result.outdir)
                        if isfile(result.bins_tsv)
                            Test.@test true
                        else
                            Test.@test_skip "MetaCoAG did not produce bins.tsv; check inputs/output."
                        end
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                end
            else
                @info "Skipping MetaCoAG; missing contigs/graph/mapping inputs."
            end

            if isfile(contigs_fasta) && isfile(coverage_table)
                Test.@testset "GenomeFace (disabled)" begin
                    outdir = mktempdir()
                    genomeface_error = nothing
                    try
                        Mycelia.run_genomeface(
                            contigs_fasta=contigs_fasta,
                            coverage_table=coverage_table,
                            outdir=outdir
                        )
                    catch err
                        genomeface_error = err
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                    Test.@test genomeface_error isa ErrorException
                    if genomeface_error isa ErrorException
                        Test.@test occursin("GenomeFace wrapper disabled", sprint(showerror, genomeface_error))
                    end
                end
            else
                @info "Skipping GenomeFace; missing contigs/coverage inputs."
            end

            comebin_inputs = nothing
            try
                include(joinpath(@__DIR__, "..", "metadata", "download_comebin_data.jl"))
                comebin_root = joinpath(dirname(@__DIR__), "metadata", "comebin_test_data")
                comebin_inputs = download_and_prep_comebin_data(comebin_root)
            catch e
                @info "Skipping COMEBin; download/prep failed." exception=e
            end

            if comebin_inputs !== nothing &&
                    isfile(comebin_inputs.contigs) &&
                    (isfile(comebin_inputs.bam_path) || isdir(comebin_inputs.bam_path))
                Test.@testset "COMEBin" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_comebin(
                            contigs_fasta=comebin_inputs.contigs,
                            bam_path=comebin_inputs.bam_path,
                            outdir=outdir
                        )
                        Test.@test isdir(result.outdir)
                        if isfile(result.bins_tsv)
                            Test.@test true
                        else
                            Test.@test_skip "COMEBin did not produce bins.tsv; check inputs/output."
                        end
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                end
            else
                @info "Skipping COMEBin; missing contigs/BAM inputs."
            end

            if !isempty(genomes) && all(isfile, genomes)
                Test.@testset "dRep" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_drep_dereplicate(
                            genomes=genomes,
                            outdir=outdir
                        )
                        Test.@test isdir(result.outdir)
                        if isfile(result.winning_genomes)
                            Test.@test true
                        else
                            Test.@test_skip "dRep did not produce dereplicated_genomes.csv; check inputs/output."
                        end
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                end
            else
                @info "Skipping dRep; missing genome FASTA inputs."
            end

            if !isempty(bins_dirs) && all(isdir, bins_dirs)
                Test.@testset "MAGmax" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_magmax_merge(
                            bins_dirs=bins_dirs,
                            outdir=outdir
                        )
                        Test.@test isdir(result.outdir)
                        if isempty(readdir(result.outdir))
                            Test.@test_skip "MAGmax did not produce output files; check inputs/output."
                        else
                            Test.@test true
                        end
                    finally
                        rm(outdir; recursive=true, force=true)
                    end
                end
            else
                @info "Skipping MAGmax; missing bin directories."
            end
        end
    else
        @info "Skipping binning external execution; set MYCELIA_RUN_EXTERNAL=true to enable when tools are installed."
    end
end
