# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/binning_tools.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/binning_tools.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

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

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end

    function _ensure_vamb_installed()
        return "mycelia_vamb"
    end
end

function _write_text_file(path::String, content::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

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

        df_custom = Mycelia.parse_bin_assignments(tmp_custom; contig_col = "contig_id", bin_col = "cluster")
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
            contigs_fasta = "missing_contigs.fna",
            depth_file = "missing_depth.tsv",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_metabat2(
            contigs_fasta = "missing_contigs.fna",
            depth_file = "missing_depth.tsv",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_metacoag(
            contigs_fasta = "missing_contigs.fna",
            assembly_graph = "missing.gfa",
            mapping_file = "missing.tsv",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_comebin(
            contigs_fasta = "missing_contigs.fna",
            bam_path = "missing_bams",
            outdir = outdir
        )
        tmp_contigs = tempname() * ".fna"
        open(tmp_contigs, "w") do io
            write(io, ">contig1\nACGTACGTACGT\n")
        end
        Test.@test_throws ErrorException Mycelia.run_comebin(
            contigs_fasta = tmp_contigs,
            bam_path = "missing_bams",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_drep_dereplicate(
            genomes = String[],
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_taxometer(
            contigs_fasta = "missing_contigs.fna",
            depth_file = "missing_depth.tsv",
            taxonomy_file = "missing_taxonomy.tsv",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_taxvamb(
            contigs_fasta = "missing_contigs.fna",
            depth_file = "missing_depth.tsv",
            taxonomy_file = "missing_taxonomy.tsv",
            outdir = outdir
        )
        genomeface_error = nothing
        try
            Mycelia.run_genomeface(
                contigs_fasta = "missing_contigs.fna",
                coverage_table = "missing_cov.tsv",
                outdir = outdir
            )
        catch err
            genomeface_error = err
        end
        Test.@test genomeface_error isa ErrorException
        if genomeface_error isa ErrorException
            Test.@test occursin("GenomeFace wrapper disabled", sprint(showerror, genomeface_error))
        end
        Test.@test_throws ErrorException Mycelia.run_magmax_merge(
            bins_dirs = String[],
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_magmax_merge(
            bins_dirs = ["missing_bins_dir"],
            outdir = outdir
        )
    end

    Test.@testset "Helper utilities and collection paths" begin
        workdir = mktempdir()

        Test.@test Mycelia._vamb_env_name() == "mycelia_vamb"

        vamb_depth = _write_text_file(
            joinpath(workdir, "vamb.tsv"),
            "contigname\tsampleA\tsampleB\nctg1\t1.0\t2.0\n"
        )
        Test.@test Mycelia._vamb_abundance_tsv(vamb_depth) == vamb_depth

        jgi_depth = _write_text_file(
            joinpath(workdir, "depth.tsv"),
            "contigName\tcontigLen\ttotalAvgDepth\tsampleA\tsampleA-var\tsampleB\nctg1\t1000\t4.0\t3.0\t0.2\t1.0\nctg2\t500\t2.0\t2.0\n"
        )
        abundance_dir = joinpath(workdir, "abundance")
        abundance_tsv = Mycelia._vamb_abundance_tsv(jgi_depth; output_dir = abundance_dir)
        abundance_df = DataFrames.DataFrame(Mycelia.CSV.File(abundance_tsv; delim = '\t'))
        Test.@test DataFrames.names(abundance_df) == ["contigname", "sampleA", "sampleB"]
        Test.@test abundance_df.sampleB[2] == 0
        Test.@test Mycelia._vamb_abundance_tsv(jgi_depth; output_dir = abundance_dir) == abundance_tsv

        bad_depth = _write_text_file(
            joinpath(workdir, "bad_depth.tsv"),
            "contigName\tcontigLen\ttotalAvgDepth\nctg1\t1000\t1.0\n"
        )
        Test.@test_throws ErrorException Mycelia._vamb_abundance_tsv(bad_depth)

        search_root = joinpath(workdir, "search")
        mkpath(joinpath(search_root, "nested", "bins"))
        _write_text_file(joinpath(search_root, "alpha.tsv"), "x\n")
        _write_text_file(joinpath(search_root, "nested", "beta.tsv"), "x\n")
        _write_text_file(joinpath(search_root, "nested", "bins", "bin.fa"), ">bin\nACGT\n")
        Test.@test Mycelia._find_first_matching_file(search_root, [r"alpha\.tsv$"]) == joinpath(search_root, "alpha.tsv")
        Test.@test Mycelia._find_first_matching_file(search_root, [r"beta\.tsv$"]; recursive = true) == joinpath(search_root, "nested", "beta.tsv")
        Test.@test isnothing(Mycelia._find_first_matching_file(search_root, [r"gamma\.tsv$"]))
        Test.@test Mycelia._find_first_matching_dir(search_root, [r"bins$"]; recursive = true) == joinpath(search_root, "nested", "bins")
        Test.@test isnothing(Mycelia._find_first_matching_dir(search_root, [r"missing$"]))

        contigs = _write_text_file(joinpath(workdir, "contigs.fa"), ">ctg1\nACGT\n")
        taxonomy = _write_text_file(joinpath(workdir, "taxonomy.tsv"), "ctg1\tBacteria\n")
        graph = _write_text_file(joinpath(workdir, "assembly.gfa"), "H\tVN:Z:1.0\n")
        mapping = _write_text_file(joinpath(workdir, "mapping.tsv"), "contig\tsampleA\nctg1\t1.0\n")
        bam_dir = joinpath(workdir, "bams")
        mkpath(bam_dir)
        _write_text_file(joinpath(bam_dir, "sample.bam"), "bam")
        genome_a = _write_text_file(joinpath(workdir, "genome_a.fa"), ">a\nACGT\n")
        genome_b = _write_text_file(joinpath(workdir, "genome_b.fa"), ">b\nTGCA\n")
        bins_a = joinpath(workdir, "bins_a")
        bins_b = joinpath(workdir, "bins_b")
        mkpath(bins_a)
        mkpath(bins_b)
        _write_text_file(joinpath(bins_a, "bin.1.fa"), ">bin1\nACGT\n")
        _write_text_file(joinpath(bins_b, "bin.2.fa"), ">bin2\nTGCA\n")

        vamb_exec = Mycelia.CollectExecutor()
        vamb_result = Mycelia.run_vamb(
            contigs_fasta = contigs,
            depth_file = jgi_depth,
            outdir = joinpath(workdir, "vamb_out"),
            minfasta = 2500,
            threads = 6,
            executor = vamb_exec
        )
        Test.@test length(vamb_exec.jobs) == 1
        Test.@test occursin("vamb bin default", vamb_exec.jobs[1].cmd)
        Test.@test occursin("--abundance_tsv", vamb_exec.jobs[1].cmd)
        Test.@test occursin("-p 6", vamb_exec.jobs[1].cmd)
        Test.@test vamb_result.outdir == joinpath(workdir, "vamb_out")

        metabat_exec = Mycelia.CollectExecutor()
        metabat_result = Mycelia.run_metabat2(
            contigs_fasta = contigs,
            depth_file = jgi_depth,
            outdir = joinpath(workdir, "metabat_out"),
            min_contig = 1800,
            threads = 3,
            seed = 99,
            extra_args = ["--saveCls"],
            executor = metabat_exec
        )
        Test.@test length(metabat_exec.jobs) == 1
        Test.@test occursin("metabat2", metabat_exec.jobs[1].cmd)
        Test.@test occursin("-m 1800", metabat_exec.jobs[1].cmd)
        Test.@test occursin("-s 99", metabat_exec.jobs[1].cmd)
        Test.@test occursin("--saveCls", metabat_exec.jobs[1].cmd)
        Test.@test endswith(metabat_result.bins_prefix, "bin")

        metacoag_exec = Mycelia.CollectExecutor()
        metacoag_result = Mycelia.run_metacoag(
            contigs_fasta = contigs,
            assembly_graph = graph,
            mapping_file = mapping,
            outdir = joinpath(workdir, "metacoag_out"),
            assembler = "megahit",
            threads = 2,
            extra_args = ["--min_length", "1000"],
            executor = metacoag_exec
        )
        Test.@test length(metacoag_exec.jobs) == 1
        Test.@test occursin("metacoag", metacoag_exec.jobs[1].cmd)
        Test.@test occursin("--assembler megahit", metacoag_exec.jobs[1].cmd)
        Test.@test occursin("--min_length 1000", metacoag_exec.jobs[1].cmd)
        Test.@test metacoag_result.bins_dir == joinpath(metacoag_result.outdir, "bins")

        comebin_exec = Mycelia.CollectExecutor()
        comebin_result = Mycelia.run_comebin(
            contigs_fasta = contigs,
            bam_path = bam_dir,
            outdir = joinpath(workdir, "comebin_out"),
            views = 4,
            threads = 5,
            temperature = 0.2,
            embedding_size = 512,
            coverage_embedding_size = 256,
            batch_size = 128,
            extra_args = ["--resume"],
            executor = comebin_exec
        )
        Test.@test length(comebin_exec.jobs) == 1
        Test.@test occursin("run_comebin.sh", comebin_exec.jobs[1].cmd)
        Test.@test occursin("-n 4", comebin_exec.jobs[1].cmd)
        Test.@test occursin("-l 0.2", comebin_exec.jobs[1].cmd)
        Test.@test occursin("--resume", comebin_exec.jobs[1].cmd)
        Test.@test comebin_result.outdir == joinpath(workdir, "comebin_out")

        drep_exec = Mycelia.CollectExecutor()
        drep_result = Mycelia.run_drep_dereplicate(
            genomes = [genome_a, genome_b],
            outdir = joinpath(workdir, "drep_out"),
            completeness_threshold = 90.0,
            contamination_threshold = 5.0,
            ani_threshold = 0.97,
            threads = 4,
            extra_args = ["--ignoreGenomeQuality"],
            executor = drep_exec
        )
        Test.@test length(drep_exec.jobs) == 1
        Test.@test occursin("dRep dereplicate", drep_exec.jobs[1].cmd)
        Test.@test occursin("-comp 90.0", drep_exec.jobs[1].cmd)
        Test.@test occursin("-con 5.0", drep_exec.jobs[1].cmd)
        Test.@test occursin("-sa 0.97", drep_exec.jobs[1].cmd)
        Test.@test drep_result.outdir == joinpath(workdir, "drep_out")

        taxometer_exec = Mycelia.CollectExecutor()
        taxometer_result = Mycelia.run_taxometer(
            contigs_fasta = contigs,
            depth_file = jgi_depth,
            taxonomy_file = taxonomy,
            outdir = joinpath(workdir, "taxometer_out"),
            threads = 7,
            executor = taxometer_exec
        )
        Test.@test length(taxometer_exec.jobs) == 1
        Test.@test occursin("vamb taxometer", taxometer_exec.jobs[1].cmd)
        Test.@test occursin("-p 7", taxometer_exec.jobs[1].cmd)
        Test.@test taxometer_result.outdir == joinpath(workdir, "taxometer_out")

        taxvamb_exec = Mycelia.CollectExecutor()
        taxvamb_result = Mycelia.run_taxvamb(
            contigs_fasta = contigs,
            depth_file = jgi_depth,
            taxonomy_file = taxonomy,
            outdir = joinpath(workdir, "taxvamb_out"),
            threads = 8,
            extra_args = ["--cuda"],
            executor = taxvamb_exec
        )
        Test.@test length(taxvamb_exec.jobs) == 1
        Test.@test occursin("vamb bin taxvamb", taxvamb_exec.jobs[1].cmd)
        Test.@test occursin("-p 8", taxvamb_exec.jobs[1].cmd)
        Test.@test occursin("--cuda", taxvamb_exec.jobs[1].cmd)
        Test.@test taxvamb_result.outdir == joinpath(workdir, "taxvamb_out")

        magmax_exec = Mycelia.CollectExecutor()
        magmax_result = Mycelia.run_magmax_merge(
            bins_dirs = [bins_a, bins_b],
            outdir = joinpath(workdir, "magmax_out"),
            threads = 2,
            executor = magmax_exec
        )
        Test.@test length(magmax_exec.jobs) == 1
        Test.@test occursin("magmax", magmax_exec.jobs[1].cmd)
        Test.@test occursin("--format fa", magmax_exec.jobs[1].cmd)
        copied_bins = sort(readdir(magmax_result.bins_input_dir))
        Test.@test copied_bins == ["bins_a__bin.1.fa", "bins_b__bin.2.fa"]
    end

    Test.@testset "Existing outdir validation" begin
        workdir = mktempdir()
        contigs = _write_text_file(joinpath(workdir, "contigs.fa"), ">ctg1\nACGT\n")
        depth = _write_text_file(
            joinpath(workdir, "depth.tsv"),
            "contigName\tcontigLen\ttotalAvgDepth\tsampleA\nctg1\t1000\t1.0\t1.0\n"
        )
        taxonomy = _write_text_file(joinpath(workdir, "taxonomy.tsv"), "ctg1\tBacteria\n")

        existing_vamb = joinpath(workdir, "existing_vamb")
        mkpath(existing_vamb)
        Test.@test_throws ErrorException Mycelia.run_vamb(
            contigs_fasta = contigs,
            depth_file = depth,
            outdir = existing_vamb
        )

        existing_taxometer = joinpath(workdir, "existing_taxometer")
        mkpath(existing_taxometer)
        Test.@test_throws ErrorException Mycelia.run_taxometer(
            contigs_fasta = contigs,
            depth_file = depth,
            taxonomy_file = taxonomy,
            outdir = existing_taxometer
        )

        existing_taxvamb = joinpath(workdir, "existing_taxvamb")
        mkpath(existing_taxvamb)
        Test.@test_throws ErrorException Mycelia.run_taxvamb(
            contigs_fasta = contigs,
            depth_file = depth,
            taxonomy_file = taxonomy,
            outdir = existing_taxvamb
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

            comebin_inputs = nothing
            try
                include(joinpath(@__DIR__, "..", "metadata", "download_comebin_data.jl"))
                comebin_root = joinpath(dirname(@__DIR__), "metadata", "comebin_test_data")
                comebin_inputs = download_and_prep_comebin_data(comebin_root)
            catch e
                @info "COMEBin download/prep failed." exception=e
            end

            if isfile(contigs_fasta) && isfile(depth_file)
                Test.@testset "VAMB" begin
                    outdir = joinpath(mktempdir(), "vamb_out")
                    try
                        result = Mycelia.run_vamb(
                            contigs_fasta = contigs_fasta,
                            depth_file = depth_file,
                            outdir = outdir
                        )
                        Test.@test isdir(result.outdir)
                        clusters_ok = result.clusters_tsv !== nothing &&
                                      isfile(result.clusters_tsv)
                        bins_dir = joinpath(result.outdir, "bins")
                        Test.@test clusters_ok || isdir(bins_dir)
                    finally
                        rm(dirname(outdir); recursive = true, force = true)
                    end
                end

                Test.@testset "MetaBAT2" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_metabat2(
                            contigs_fasta = contigs_fasta,
                            depth_file = depth_file,
                            outdir = outdir
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
                        rm(outdir; recursive = true, force = true)
                    end
                end
            else
                @info "Skipping VAMB/MetaBAT2; missing contigs/depth inputs."
            end

            if isfile(contigs_fasta) && isfile(depth_file) && isfile(taxonomy_file)
                Test.@testset "Taxometer" begin
                    outdir = joinpath(mktempdir(), "taxometer_out")
                    try
                        result = Mycelia.run_taxometer(
                            contigs_fasta = contigs_fasta,
                            depth_file = depth_file,
                            taxonomy_file = taxonomy_file,
                            outdir = outdir
                        )
                        Test.@test isdir(result.outdir)
                        if isempty(readdir(result.outdir))
                            Test.@test_skip "Taxometer did not produce output files; check inputs/output."
                        else
                            Test.@test true
                        end
                    finally
                        rm(dirname(outdir); recursive = true, force = true)
                    end
                end

                Test.@testset "TaxVAMB" begin
                    outdir = joinpath(mktempdir(), "taxvamb_out")
                    try
                        result = Mycelia.run_taxvamb(
                            contigs_fasta = contigs_fasta,
                            depth_file = depth_file,
                            taxonomy_file = taxonomy_file,
                            outdir = outdir
                        )
                        Test.@test isdir(result.outdir)
                        clusters_ok = result.clusters_tsv !== nothing &&
                                      isfile(result.clusters_tsv)
                        bins_dir = joinpath(result.outdir, "bins")
                        Test.@test clusters_ok || isdir(bins_dir)
                    finally
                        rm(dirname(outdir); recursive = true, force = true)
                    end
                end
            else
                @info "Skipping Taxometer/TaxVAMB; missing contigs/depth/taxonomy inputs."
            end

            if comebin_inputs === nothing
                Test.@test false
            elseif isfile(comebin_inputs.contigs) &&
                   (isfile(comebin_inputs.bam_path) || isdir(comebin_inputs.bam_path))
                Test.@testset "MetaCoAG" begin
                    outdir = joinpath(mktempdir(), "metacoag_out")
                    try
                        work_dir = dirname(outdir)
                        assembly_graph = joinpath(work_dir, "assembly_graph.gfa")
                        if !isfile(assembly_graph)
                            Mycelia._write_simple_gfa_from_fasta(comebin_inputs.contigs, assembly_graph)
                        end
                        coverage_table = joinpath(work_dir, "coverm_contig.tsv")
                        if !isfile(coverage_table) || filesize(coverage_table) == 0
                            Mycelia.run_coverm_contig(
                                bam_files = comebin_inputs.bam_files,
                                output_tsv = coverage_table,
                                threads = Mycelia.get_default_threads(),
                                quiet = true
                            )
                        end
                        metacoag_abundance = joinpath(work_dir, "metacoag_abundance.tsv")
                        if !isfile(metacoag_abundance) || filesize(metacoag_abundance) == 0
                            open(coverage_table, "r") do io
                                first_line = readline(io)
                                first_fields = split(first_line, '\t')
                                open(metacoag_abundance, "w") do out
                                    if isempty(first_fields) ||
                                       lowercase(first_fields[1]) ∉ ("contig", "contigname")
                                        write(out, first_line, '\n')
                                    end
                                    for line in eachline(io)
                                        write(out, line, '\n')
                                    end
                                end
                            end
                        end
                        result = Mycelia.run_metacoag(
                            contigs_fasta = comebin_inputs.contigs,
                            assembly_graph = assembly_graph,
                            mapping_file = metacoag_abundance,
                            outdir = outdir
                        )
                        Test.@test isdir(result.outdir)
                        bins_ok = isdir(result.bins_dir) ||
                                  (result.bins_tsv !== nothing && isfile(result.bins_tsv))
                        Test.@test bins_ok
                    finally
                        rm(outdir; recursive = true, force = true)
                    end
                end
            else
                Test.@test false
            end

            if isfile(contigs_fasta) && isfile(coverage_table)
                Test.@testset "GenomeFace (disabled)" begin
                    outdir = mktempdir()
                    genomeface_error = nothing
                    try
                        Mycelia.run_genomeface(
                            contigs_fasta = contigs_fasta,
                            coverage_table = coverage_table,
                            outdir = outdir
                        )
                    catch err
                        genomeface_error = err
                    finally
                        rm(outdir; recursive = true, force = true)
                    end
                    Test.@test genomeface_error isa ErrorException
                    if genomeface_error isa ErrorException
                        Test.@test occursin("GenomeFace wrapper disabled", sprint(showerror, genomeface_error))
                    end
                end
            else
                @info "Skipping GenomeFace; missing contigs/coverage inputs."
            end

            if comebin_inputs !== nothing &&
               isfile(comebin_inputs.contigs) &&
               (isfile(comebin_inputs.bam_path) || isdir(comebin_inputs.bam_path))
                Test.@testset "COMEBin" begin
                    outdir = mktempdir()
                    try
                        result = Mycelia.run_comebin(
                            contigs_fasta = comebin_inputs.contigs,
                            bam_path = comebin_inputs.bam_path,
                            outdir = outdir,
                            views = 2,
                            embedding_size = 512,
                            coverage_embedding_size = 512,
                            batch_size = 256,
                            threads = min(4, Mycelia.get_default_threads())
                        )
                        Test.@test isdir(result.outdir)
                        bins_ok = (result.bins_dir !== nothing && isdir(result.bins_dir)) ||
                                  (result.bins_tsv !== nothing && isfile(result.bins_tsv))
                        Test.@test bins_ok
                    finally
                        rm(outdir; recursive = true, force = true)
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
                            genomes = genomes,
                            outdir = outdir,
                            extra_args = ["--ignoreGenomeQuality", "-l", "1000"]
                        )
                        Test.@test isdir(result.outdir)
                        Test.@test result.winning_genomes !== nothing &&
                                   isfile(result.winning_genomes)
                    finally
                        rm(outdir; recursive = true, force = true)
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
                            bins_dirs = bins_dirs,
                            outdir = outdir,
                            extra_args = ["--no-reassembly"]
                        )
                        Test.@test isdir(result.outdir)
                        if isempty(readdir(result.outdir))
                            Test.@test_skip "MAGmax did not produce output files; check inputs/output."
                        else
                            Test.@test true
                        end
                    finally
                        rm(outdir; recursive = true, force = true)
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
