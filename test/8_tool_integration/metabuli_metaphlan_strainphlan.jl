import Test
import Mycelia
import DataFrames
import StableRNGs

const RUN_ALL = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
const RUN_EXTERNAL = RUN_ALL || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

Test.@testset "Metabuli / MetaPhlAn / StrainPhlAn integration" begin

    if !RUN_EXTERNAL
        @info "Skipping integration tests; set MYCELIA_RUN_EXTERNAL=true or MYCELIA_RUN_ALL=true to run external tools"
        return
    end

    workdir = mktempdir()
    genomes = Mycelia.get_test_genome_fasta.(;
        use_ncbi=true,
        accession=["NC_001422.1", "GCF_000819615.1"],
    )

    refs = [g.fasta for g in genomes if isfile(g.fasta) && filesize(g.fasta) > 0]
    if length(refs) < 2
        @warn "Could not download reference genomes; skipping integration suite"
        rm(workdir; recursive=true, force=true)
        return
    end

    rng = StableRNGs.StableRNG(42)

    illumina_a = Mycelia.simulate_illumina_reads(
        fasta=refs[1]; coverage=10, read_length=150, paired=true, quiet=true, rng=rng
    )
    illumina_b = Mycelia.simulate_illumina_reads(
        fasta=refs[2]; coverage=5, read_length=150, paired=true, quiet=true, rng=rng
    )

    r1 = Mycelia.concatenate_fastx(
        [illumina_a.forward_reads, illumina_b.forward_reads];
        output_path=joinpath(workdir, "illumina_R1.fq.gz")
    )
    r2 = Mycelia.concatenate_fastx(
        [illumina_a.reverse_reads, illumina_b.reverse_reads];
        output_path=joinpath(workdir, "illumina_R2.fq.gz")
    )

    long_reads = Mycelia.simulate_nanopore_reads(
        fasta=refs[1], quantity="10x", outfile=joinpath(workdir, "long_reads.fq.gz"), quiet=true
    )

    try
        # ------------------------------------------------------------------
        # Metabuli: short, long, and contig modes
        # ------------------------------------------------------------------
        metabuli_db = get(ENV, "METABULI_DB", "")
        if isempty(metabuli_db) || !isdir(metabuli_db)
            @info "Skipping Metabuli tests; METABULI_DB not set"
        else
            Test.@testset "Metabuli classify (short/long/contig)" begin
                short_out = mkpath(joinpath(workdir, "metabuli_short"))
                short_res = Mycelia.run_metabuli_classify(
                    r1;
                    reads2=r2,
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
                    long_reads;
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
                    refs[1];
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
        end

        # ------------------------------------------------------------------
        # MetaPhlAn: short and long modes with mapout for StrainPhlAn
        # ------------------------------------------------------------------
        metaphlan_db_dir = get(ENV, "METAPHLAN_DB_DIR", nothing)
        bowtie2db = get(ENV, "METAPHLAN_BOWTIE2DB", nothing)

        if metaphlan_db_dir === nothing && bowtie2db === nothing
            @info "Skipping MetaPhlAn/StrainPhlAn tests; METAPHLAN_DB_DIR or METAPHLAN_BOWTIE2DB not set"
        else
            Test.@testset "MetaPhlAn profiling (short/long)" begin
                mp_out = mkpath(joinpath(workdir, "metaphlan_short"))
                mp_res = Mycelia.run_metaphlan(
                    r1;
                    reads2=r2,
                    sample_name="short_sample",
                    outdir=mp_out,
                    db_dir=metaphlan_db_dir,
                    bowtie2db=bowtie2db,
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
                    long_reads;
                    sample_name="long_sample",
                    outdir=mp_long_out,
                    db_dir=metaphlan_db_dir,
                    bowtie2db=bowtie2db,
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
                    r1;
                    reads2=r2,
                    sample_name="strain_sample1",
                    outdir=joinpath(workdir, "strain_mp1"),
                    db_dir=metaphlan_db_dir,
                    bowtie2db=bowtie2db,
                    input_type=:fastq,
                    mapout=true,
                    threads=2,
                    force=true
                ).map_output)
                push!(map_files, Mycelia.run_metaphlan(
                    long_reads;
                    sample_name="strain_sample2",
                    outdir=joinpath(workdir, "strain_mp2"),
                    db_dir=metaphlan_db_dir,
                    bowtie2db=bowtie2db,
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
        end
    finally
        rm(workdir; recursive=true, force=true)
        for g in genomes
            if hasproperty(g, :cleanup) && g.cleanup !== nothing
                g.cleanup()
            end
        end
    end
end
