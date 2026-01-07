# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/metabuli_metaphlan_strainphlan.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/metabuli_metaphlan_strainphlan.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import DataFrames
import StableRNGs

const RUN_ALL = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
const RUN_EXTERNAL = RUN_ALL || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

Test.@testset "Metabuli / MetaPhlAn / StrainPhlAn integration" begin

    if !RUN_EXTERNAL
        @info "Skipping integration tests; set MYCELIA_RUN_EXTERNAL=true to run external tools"
        return
    end

    workdir = mktempdir()
    virus_accession = "GCF_000819615.1"
    bacteria_accession = "GCF_000005845.2"
    virus_genome = Mycelia.get_test_genome_fasta(; use_ncbi=true, accession=virus_accession)
    bacteria_genome = Mycelia.get_test_genome_fasta(; use_ncbi=true, accession=bacteria_accession)

    virus_ref = virus_genome.fasta
    bacteria_ref = bacteria_genome.fasta
    if !isfile(virus_ref) || filesize(virus_ref) == 0 || !isfile(bacteria_ref) || filesize(bacteria_ref) == 0
        @warn "Could not download reference genomes; skipping integration suite"
        rm(workdir; recursive=true, force=true)
        return
    end

    rng = StableRNGs.StableRNG(42)

    virus_reads = Mycelia.simulate_illumina_reads(
        fasta=virus_ref;
        read_count=12000,
        read_length=150,
        paired=true,
        quiet=true,
        outbase=joinpath(workdir, "virus_illumina"),
        rndSeed=rand(rng, 0:typemax(Int))
    )
    virus_r1 = virus_reads.forward_reads
    virus_r2 = virus_reads.reverse_reads

    bacteria_reads = Mycelia.simulate_illumina_reads(
        fasta=bacteria_ref;
        read_count=20000,
        read_length=150,
        paired=true,
        quiet=true,
        outbase=joinpath(workdir, "bacteria_illumina"),
        rndSeed=rand(rng, 0:typemax(Int))
    )
    mp_r1 = bacteria_reads.forward_reads
    mp_r2 = bacteria_reads.reverse_reads

    mp_long_reads = Mycelia.simulate_nearly_perfect_long_reads(
        fasta=bacteria_ref,
        quantity="20000000",
        length_mean=2000,
        length_sd=500,
        outfile=joinpath(workdir, "metaphlan_long_reads.fq.gz"),
        quiet=true
    )

    virus_long_reads = Mycelia.simulate_nearly_perfect_long_reads(
        fasta=virus_ref,
        quantity="5000000",
        length_mean=2000,
        length_sd=500,
        outfile=joinpath(workdir, "metabuli_long_reads.fq.gz"),
        quiet=true
    )

    try
        # ------------------------------------------------------------------
        # Metabuli: short, long, and contig modes
        # ------------------------------------------------------------------
        metabuli_db = Mycelia.get_metabuli_db_path(db_name="refseq_virus")
        Test.@testset "Metabuli classify (short/long/contig)" begin
            short_out = mkpath(joinpath(workdir, "metabuli_short"))
            short_res = Mycelia.run_metabuli_classify(
                virus_r1;
                reads2=virus_r2,
                dbdir=metabuli_db,
                outdir=short_out,
                jobid="short_sample",
                read_platform=:illumina,
                threads=2,
                force=true
            )
            Test.@test isfile(short_res.classifications_tsv)
            Test.@test isfile(short_res.report_tsv)
            Test.@test DataFrames.nrow(short_res.report) >= 0

            long_out = mkpath(joinpath(workdir, "metabuli_long"))
            long_res = Mycelia.run_metabuli_classify(
                virus_long_reads;
                reads2=nothing,
                dbdir=metabuli_db,
                outdir=long_out,
                jobid="long_sample",
                read_platform=:ont,
                threads=2,
                force=true
            )
            Test.@test isfile(long_res.classifications_tsv)
            Test.@test isfile(long_res.report_tsv)

            contig_out = mkpath(joinpath(workdir, "metabuli_contig"))
            contig_res = Mycelia.run_metabuli_classify(
                virus_ref;
                reads2=nothing,
                dbdir=metabuli_db,
                outdir=contig_out,
                jobid="contig_sample",
                read_platform=:contig,
                threads=2,
                force=true
            )
            Test.@test isfile(contig_res.classifications_tsv)
            Test.@test isfile(contig_res.report_tsv)
        end

        # ------------------------------------------------------------------
        # MetaPhlAn: short and long modes with mapout for StrainPhlAn
        # ------------------------------------------------------------------
        metaphlan_db_dir = Mycelia.get_metaphlan_db_path()
        Test.@testset "MetaPhlAn profiling (short/long)" begin
            mp_out = mkpath(joinpath(workdir, "metaphlan_short"))
            mp_res = Mycelia.run_metaphlan(
                mp_r1;
                reads2=mp_r2,
                sample_name="short_sample",
                outdir=mp_out,
                db_dir=metaphlan_db_dir,
                input_type=:fastq,
                long_reads=false,
                mapout=true,
                threads=2,
                force=true
            )
            Test.@test isfile(mp_res.profile_tsv)
            Test.@test !isempty(mp_res.table)
            Test.@test mp_res.map_output !== nothing && isfile(mp_res.map_output)

            mp_long_out = mkpath(joinpath(workdir, "metaphlan_long"))
            mp_long = Mycelia.run_metaphlan(
                mp_long_reads;
                sample_name="long_sample",
                outdir=mp_long_out,
                db_dir=metaphlan_db_dir,
                input_type=:fastq,
                long_reads=true,
                mapout=true,
                threads=2,
                force=true
            )
            Test.@test isfile(mp_long.profile_tsv)
            Test.@test mp_long.map_output !== nothing && isfile(mp_long.map_output)
        end

        # --------------------------------------------------------------
        # StrainPhlAn: run sample2markers and strainphlan on MetaPhlAn outputs
        # --------------------------------------------------------------
        Test.@testset "StrainPhlAn end-to-end" begin
            map_files = String[]
            push!(map_files, Mycelia.run_metaphlan(
                mp_r1;
                reads2=mp_r2,
                sample_name="strain_sample1",
                outdir=joinpath(workdir, "strain_mp1"),
                db_dir=metaphlan_db_dir,
                input_type=:fastq,
                mapout=true,
                threads=2,
                force=true
            ).map_output)
            push!(map_files, Mycelia.run_metaphlan(
                mp_long_reads;
                sample_name="strain_sample2",
                outdir=joinpath(workdir, "strain_mp2"),
                db_dir=metaphlan_db_dir,
                input_type=:fastq,
                long_reads=true,
                mapout=true,
                threads=2,
                force=true
            ).map_output)

            marker_dir = mkpath(joinpath(workdir, "markers"))
            s2m = Mycelia.run_sample2markers(
                map_files;
                outdir=marker_dir,
                threads=2,
                force=true
            )
            Test.@test !isempty(s2m.marker_fastas)
            Test.@test all(isfile.(s2m.marker_fastas))

            strain_out = mkpath(joinpath(workdir, "strainphlan"))
            sp_res = Mycelia.run_strainphlan(
                s2m.marker_fastas;
                outdir=strain_out,
                clade="synthetic_clade",
                threads=2,
                force=true
            )
            Test.@test sp_res.tree_file === nothing || isfile(sp_res.tree_file)
            Test.@test sp_res.marker_alignment === nothing || isfile(sp_res.marker_alignment)
        end
    finally
        rm(workdir; recursive=true, force=true)
        for g in (virus_genome, bacteria_genome)
            if hasproperty(g, :cleanup) && g.cleanup !== nothing
                g.cleanup()
            end
        end
    end
end
