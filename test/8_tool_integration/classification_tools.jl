# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/classification_tools.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/classification_tools.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Classification Tools Integration Tests
# To run with external tool execution (requires conda/network):

import Test
import Mycelia
import DataFrames
import BioSequences
import FASTX
import CodecZlib
import Random

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function _write_text_file(path::String, content::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

function _make_metaphlan_db(dir::String; index::String = "mpa_vTest")
    mkpath(dir)
    _write_text_file(joinpath(dir, "$(index).pkl"), "db")
    return dir
end

function _make_metabuli_db(dir::String)
    mkpath(dir)
    _write_text_file(joinpath(dir, "db.info"), "db")
    return dir
end

function _make_kraken2_db(dir::String)
    mkpath(dir)
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d")
        _write_text_file(joinpath(dir, name), "db")
    end
    return dir
end

# Check if external tool tests should run
const RUN_ALL = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
const RUN_EXTERNAL = RUN_ALL || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
# const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

Test.@testset "Classification Tools Integration Tests" begin

    # ========================================================================
    # Sourmash Tests
    # ========================================================================
    Test.@testset "Sourmash Integration" begin
        if RUN_EXTERNAL
            Test.@testset "Sourmash sketch, search, and gather" begin
                workdir = mktempdir()

                # Download PhiX genome as test data
                phix_fasta = Mycelia.download_genome_by_accession(
                    accession = "NC_001422.1",
                    outdir = workdir,
                    compressed = false
                )
                Test.@test isfile(phix_fasta)

                # Create a mutated version to have two different sequences
                fasta_records = collect(Mycelia.open_fastx(phix_fasta))
                seq_record = first(fasta_records)
                seq_str = String(FASTX.sequence(seq_record))

                mutated_seq = Mycelia.mutate_dna_substitution_fraction(seq_str; fraction = 0.05)
                mutated_fasta = joinpath(workdir, "phix_mutated.fasta")
                Mycelia.write_fasta(
                    outfile = mutated_fasta,
                    records = [FASTX.FASTA.Record("phix_mutated", BioSequences.LongDNA{4}(mutated_seq))]
                )
                Test.@test isfile(mutated_fasta)

                # Test 1: Create signatures (sketches)
                sketch_result = Mycelia.run_sourmash_sketch(
                    input_files = [phix_fasta, mutated_fasta],
                    outdir = joinpath(workdir, "sketches"),
                    k_sizes = [21, 31],
                    scaled = 100,
                    molecule = "dna"
                )

                Test.@test isdir(sketch_result.outdir)
                Test.@test length(sketch_result.signatures) == 2
                for sig in sketch_result.signatures
                    Test.@test isfile(sig)
                    Test.@test endswith(sig, ".sig")
                end

                # Test 2: Search one signature against the other
                search_result = Mycelia.run_sourmash_search(
                    query_sig = sketch_result.signatures[1],
                    database_sig = sketch_result.signatures[2],
                    outdir = joinpath(workdir, "search_results"),
                    threshold = 0.01,
                    k_size = 31
                )

                Test.@test isdir(search_result.outdir)
                Test.@test isfile(search_result.results_csv)

                # Read and verify search results
                search_df = DataFrames.DataFrame(Mycelia.CSV.File(search_result.results_csv))
                Test.@test DataFrames.nrow(search_df) >= 0  # May have matches or not

                # Test 3: Create a combined signature for gather test
                combined_fasta = joinpath(workdir, "combined.fasta")
                open(combined_fasta, "w") do io
                    for record in collect(Mycelia.open_fastx(phix_fasta))
                        write(io, ">$(FASTX.identifier(record))\n$(FASTX.sequence(record))\n")
                    end
                    for record in collect(Mycelia.open_fastx(mutated_fasta))
                        write(io, ">$(FASTX.identifier(record))\n$(FASTX.sequence(record))\n")
                    end
                end

                combined_sketch = Mycelia.run_sourmash_sketch(
                    input_files = [combined_fasta],
                    outdir = joinpath(workdir, "combined_sketch"),
                    k_sizes = [31],
                    scaled = 100
                )

                # Create a database from the reference signatures
                db_dir = joinpath(workdir, "sig_db")
                mkpath(db_dir)
                for sig in sketch_result.signatures
                    cp(sig, joinpath(db_dir, basename(sig)))
                end

                gather_result = Mycelia.run_sourmash_gather(
                    query_sig = combined_sketch.signatures[1],
                    database_sig = db_dir,
                    outdir = joinpath(workdir, "gather_results"),
                    k_size = 31,
                    threshold_bp = 100
                )

                Test.@test isdir(gather_result.outdir)
                Test.@test isfile(gather_result.results_csv)

                # Cleanup
                rm(workdir; recursive = true, force = true)
            end
        else
            @info "Skipping sourmash integration test; external tool execution is opt-in via MYCELIA_RUN_EXTERNAL=true"
        end
    end

    # ========================================================================
    # MetaPhlAn Tests
    # ========================================================================
    Test.@testset "MetaPhlAn Integration" begin
        Test.@testset "MetaPhlAn profile parsing" begin
            # Test parsing with mock profile data (no external tools needed)
            temp_profile = tempname()

            # Create a realistic MetaPhlAn profile
            open(temp_profile, "w") do io
                write(io, "#mpa_vJan21_CHOCOPhlAnSGB_202103\n")
                write(io, "#/usr/local/bin/metaphlan reads.fq --input_type fastq\n")
                write(io, "#SampleID\tMetaphlan_Analysis\n")
                write(io, "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n")
                write(io, "k__Bacteria\t2\t100.0\t\n")
                write(io, "k__Bacteria|p__Proteobacteria\t1224\t65.5\t\n")
                write(io, "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria\t1236\t45.2\t\n")
                write(io, "k__Bacteria|p__Firmicutes\t1239\t34.5\t\n")
                write(io, "k__Bacteria|p__Firmicutes|c__Bacilli\t91061\t20.1\t\n")
            end

            try
                df = Mycelia.parse_metaphlan_profile(temp_profile)

                Test.@test df isa DataFrames.DataFrame
                Test.@test DataFrames.nrow(df) == 5
                Test.@test :clade_name in DataFrames.propertynames(df)
                Test.@test :relative_abundance in DataFrames.propertynames(df)

                # Verify specific values
                Test.@test df.clade_name[1] == "k__Bacteria"
                Test.@test df.relative_abundance[1] == 100.0
                Test.@test df.relative_abundance[2] == 65.5

                # Verify taxonomic hierarchy is preserved
                Test.@test occursin("Proteobacteria", df.clade_name[2])
                Test.@test occursin("Firmicutes", df.clade_name[4])
            finally
                rm(temp_profile, force = true)
            end
        end

        if RUN_EXTERNAL
            Test.@testset "MetaPhlAn full execution" begin
                workdir = mktempdir()

                try
                    # Download real genome (with fallback to simulated if NCBI unavailable)
                    test_genome_info = Mycelia.get_test_genome_fasta(
                        use_ncbi = true,
                        accession = "GCF_000005845.2"  # E. coli K-12 MG1655
                    )
                    test_fasta = test_genome_info.fasta

                    if !isfile(test_fasta) || filesize(test_fasta) == 0
                        @warn "Failed to obtain test genome for MetaPhlAn test"
                        Test.@test_broken false
                    else
                        # Simulate realistic paired-end Illumina reads (100bp, >70bp requirement)
                        sim_result = Mycelia.simulate_illumina_hs20_100bp(
                            fasta = test_fasta,
                            read_count = 10000,
                            quiet = true
                        )
                        forward_reads = sim_result.forward_reads
                        reverse_reads = sim_result.reverse_reads

                        # Create synthetic reads with random sequences to verify
                        # tool doesn't misclassify fake data (using randdnaseq for valid DNA)
                        Random.seed!(42)  # For reproducibility
                        synthetic_records = [FASTX.FASTQ.Record(
                                                 "synthetic_read$(i)",
                                                 BioSequences.randdnaseq(Random.default_rng(), 100),
                                                 fill(Int8(40), 100)  # Quality score 40 ('I')
                                             )
                                             for i in 1:50]
                        synthetic_fastq = joinpath(workdir, "synthetic_reads.fq")
                        Mycelia.write_fastq(records = synthetic_records, filename = synthetic_fastq)

                        # Test with single-end reads (forward + synthetic)
                        combined_fastq = joinpath(workdir, "combined_reads.fq")
                        Mycelia.concatenate_fastx(
                            [forward_reads, synthetic_fastq], output_path = combined_fastq)

                        result = Mycelia.run_metaphlan(
                            input_file = combined_fastq,
                            outdir = joinpath(workdir, "metaphlan_output_single"),
                            input_type = "fastq",
                            nprocs = 2
                        )

                        Test.@test isdir(result.outdir)
                        Test.@test isfile(result.profile_txt)

                        # Parse the output
                        if isfile(result.profile_txt) && filesize(result.profile_txt) > 0
                            df = Mycelia.parse_metaphlan_profile(result.profile_txt)
                            Test.@test df isa DataFrames.DataFrame
                        end
                    end

                    # Cleanup genome if downloaded
                    if test_genome_info.cleanup !== nothing
                        test_genome_info.cleanup()
                    end
                finally
                    rm(workdir; recursive = true, force = true)
                end
            end
        else
            @info "Skipping MetaPhlAn execution test; external tool execution is opt-in via MYCELIA_RUN_EXTERNAL=true"
        end
    end

    # ========================================================================
    # Metabuli Tests
    # ========================================================================
    Test.@testset "Metabuli Integration" begin
        Test.@testset "Metabuli report parsing" begin
            # Test parsing with actual Metabuli report format (no header row)
            # Format: percentage, num_reads, num_direct_reads, rank, taxid, name
            temp_report = tempname()

            # Create a realistic Metabuli report (no header, tab-separated)
            open(temp_report, "w") do io
                write(io, "0.6970\t789\t789\tno rank\t0\tunclassified\n")
                write(io, "99.3030\t112410\t1\tno rank\t1\troot\n")
                write(io, "45.50\t1000\t455\tspecies\t562\tEscherichia coli\n")
                write(io, "30.20\t800\t302\tspecies\t632\tYersinia pestis\n")
            end

            try
                df = Mycelia.parse_metabuli_report(temp_report)

                Test.@test df isa DataFrames.DataFrame
                Test.@test DataFrames.nrow(df) == 4

                # Verify all expected columns exist
                Test.@test :percentage in DataFrames.propertynames(df)
                Test.@test :num_reads in DataFrames.propertynames(df)
                Test.@test :num_direct_reads in DataFrames.propertynames(df)
                Test.@test :rank in DataFrames.propertynames(df)
                Test.@test :taxid in DataFrames.propertynames(df)
                Test.@test :name in DataFrames.propertynames(df)

                # Verify data types
                Test.@test eltype(df.percentage) <: AbstractFloat
                Test.@test eltype(df.num_reads) <: Integer
                Test.@test eltype(df.num_direct_reads) <: Integer
                Test.@test eltype(df.taxid) <: Integer
                Test.@test eltype(df.rank) <: AbstractString
                Test.@test eltype(df.name) <: AbstractString

                # Verify specific values
                Test.@test df.percentage[1] ≈ 0.6970
                Test.@test df.num_reads[1] == 789
                Test.@test df.taxid[3] == 562
                Test.@test df.name[3] == "Escherichia coli"
            finally
                rm(temp_report, force = true)
            end
        end

        Test.@testset "Metabuli report vcat (combining multiple samples)" begin
            # Test that multiple parsed reports can be combined with vcat
            # This was the original bug: different samples had different "column names"
            # because the parser was treating data rows as headers
            temp_report1 = tempname()
            temp_report2 = tempname()

            # Create two different sample reports
            open(temp_report1, "w") do io
                write(io, "50.00\t500\t250\tspecies\t562\tEscherichia coli\n")
                write(io, "30.00\t300\t150\tspecies\t632\tYersinia pestis\n")
                write(io, "20.00\t200\t100\tno rank\t0\tunclassified\n")
            end

            open(temp_report2, "w") do io
                write(io, "60.00\t600\t300\tspecies\t287\tPseudomonas aeruginosa\n")
                write(io, "25.00\t250\t125\tspecies\t562\tEscherichia coli\n")
                write(io, "15.00\t150\t75\tno rank\t0\tunclassified\n")
            end

            try
                df1 = Mycelia.parse_metabuli_report(temp_report1)
                df2 = Mycelia.parse_metabuli_report(temp_report2)

                # Add sample identifiers (as notebooks do)
                df1[!, :sample] .= "sample1"
                df2[!, :sample] .= "sample2"

                # This should NOT throw an error - the original bug caused
                # ArgumentError due to mismatched column names
                combined = DataFrames.vcat(df1, df2)

                Test.@test DataFrames.nrow(combined) == 6
                Test.@test length(unique(combined.sample)) == 2
                Test.@test "sample1" in combined.sample
                Test.@test "sample2" in combined.sample
            finally
                rm(temp_report1, force = true)
                rm(temp_report2, force = true)
            end
        end

        Test.@testset "Metabuli classifications parsing" begin
            # Test parsing per-read classifications
            temp_classifications = tempname()

            open(temp_classifications, "w") do io
                write(io, "C\tread_001\t562\n")
                write(io, "C\tread_002\t562\n")
                write(io, "C\tread_003\t632\n")
                write(io, "U\tread_004\t0\n")
                write(io, "C\tread_005\t287\n")
            end

            try
                df = Mycelia.parse_metabuli_classifications(temp_classifications)

                Test.@test df isa DataFrames.DataFrame
                Test.@test DataFrames.nrow(df) == 5
                Test.@test :classified in DataFrames.propertynames(df)
                Test.@test :read_id in DataFrames.propertynames(df)
                Test.@test :taxid in DataFrames.propertynames(df)

                # Verify classification status
                Test.@test df.classified[1] == "C"  # Classified
                Test.@test df.classified[4] == "U"  # Unclassified

                # Verify read IDs
                Test.@test df.read_id[1] == "read_001"
                Test.@test df.read_id[5] == "read_005"
            finally
                rm(temp_classifications, force = true)
            end
        end

        if RUN_EXTERNAL
            Test.@testset "Metabuli full execution" begin
                # Use the virus database for unit testing to avoid larger downloads.
                db_path = Mycelia.get_metabuli_db_path(db_name = "refseq_virus")
                workdir = mktempdir()

                try
                    # Download real genome (with fallback to simulated if NCBI unavailable)
                    test_genome_info = Mycelia.get_test_genome_fasta(
                        use_ncbi = true,
                        accession = "GCF_000819615.1"  # phiX174 (virus)
                    )
                    test_fasta = test_genome_info.fasta

                    if !isfile(test_fasta) || filesize(test_fasta) == 0
                        @warn "Failed to obtain test genome for Metabuli test"
                        Test.@test_broken false
                    else
                        # Simulate realistic paired-end Illumina reads
                        sim_result = Mycelia.simulate_illumina_hs20_100bp(
                            fasta = test_fasta,
                            read_count = 10000,
                            quiet = true
                        )
                        forward_reads = sim_result.forward_reads
                        reverse_reads = sim_result.reverse_reads

                        # Create synthetic reads with random sequences
                        Random.seed!(42)
                        synthetic_records = [FASTX.FASTQ.Record(
                                                 "synthetic_read$(i)",
                                                 BioSequences.randdnaseq(Random.default_rng(), 100),
                                                 fill(Int8(40), 100)
                                             )
                                             for i in 1:50]
                        synthetic_fastq = joinpath(workdir, "synthetic_reads.fq")
                        Mycelia.write_fastq(records = synthetic_records, filename = synthetic_fastq)

                        # Test with single-end reads (forward + synthetic)
                        combined_fastq = joinpath(workdir, "combined_reads.fq.gz")
                        Mycelia.concatenate_fastx(
                            [forward_reads, synthetic_fastq], output_path = combined_fastq)

                        result = Mycelia.run_metabuli_classify(
                            input_files = [combined_fastq],
                            database_path = db_path,
                            outdir = joinpath(workdir, "metabuli_output"),
                            seq_mode = "1",
                            threads = 1
                        )

                        Test.@test isdir(result.outdir)
                        Test.@test isfile(result.report_file)
                        Test.@test isfile(result.classifications_file)

                        # Parse and verify outputs
                        if isfile(result.report_file) && filesize(result.report_file) > 0
                            report_df = Mycelia.parse_metabuli_report(result.report_file)
                            Test.@test report_df isa DataFrames.DataFrame
                        end

                        if isfile(result.classifications_file) &&
                           filesize(result.classifications_file) > 0
                            class_df = Mycelia.parse_metabuli_classifications(result.classifications_file)
                            Test.@test class_df isa DataFrames.DataFrame
                        end
                    end

                    # Cleanup genome if downloaded
                    if test_genome_info.cleanup !== nothing
                        test_genome_info.cleanup()
                    end
                finally
                    rm(workdir; recursive = true, force = true)
                end
            end
        else
            @info "Skipping Metabuli execution test; external tool execution is opt-in via MYCELIA_RUN_EXTERNAL=true"
        end
    end

    # ========================================================================
    # Fast Wrapper Coverage Tests
    # ========================================================================
    Test.@testset "Sourmash wrapper collection and cache paths" begin
        workdir = mktempdir()

        fasta_a = _write_text_file(joinpath(workdir, "a.fasta"), ">a\nACGT\n")
        fasta_b = _write_text_file(joinpath(workdir, "b.fastq"), "@b\nACGT\n+\nIIII\n")
        query_sig = _write_text_file(joinpath(workdir, "query.sig"), "sig")
        db_sig = _write_text_file(joinpath(workdir, "db.sig"), "sig")

        sketch_exec = Mycelia.CollectExecutor()
        sketch_result = Mycelia.run_sourmash_sketch(
            input_files = [fasta_a, fasta_b],
            outdir = joinpath(workdir, "sourmash_sketch"),
            k_sizes = [21, 31],
            scaled = 500,
            singleton = true,
            name = "demo",
            executor = sketch_exec,
            site = :scg
        )
        Test.@test length(sketch_exec.jobs) == 1
        Test.@test sketch_exec.jobs[1].site == :scg
        Test.@test occursin("sourmash sketch dna", sketch_exec.jobs[1].cmd)
        Test.@test occursin("--singleton", sketch_exec.jobs[1].cmd)
        Test.@test occursin("--name demo", sketch_exec.jobs[1].cmd)
        Test.@test occursin("k=21,k=31,scaled=500", sketch_exec.jobs[1].cmd)
        Test.@test length(sketch_result.signatures) == 2

        search_exec = Mycelia.CollectExecutor()
        search_result = Mycelia.run_sourmash_search(
            query_sig = query_sig,
            database_sig = db_sig,
            outdir = joinpath(workdir, "sourmash_search"),
            threshold = 0.2,
            k_size = 51,
            best_only = true,
            num_results = 3,
            executor = search_exec
        )
        Test.@test length(search_exec.jobs) == 1
        Test.@test occursin("sourmash search", search_exec.jobs[1].cmd)
        Test.@test occursin("--best-only", search_exec.jobs[1].cmd)
        Test.@test occursin("-n 3", search_exec.jobs[1].cmd)
        Test.@test occursin("--threshold 0.2", search_exec.jobs[1].cmd)
        Test.@test endswith(search_result.results_csv, "_search_results.csv")

        gather_exec = Mycelia.CollectExecutor()
        gather_result = Mycelia.run_sourmash_gather(
            query_sig = query_sig,
            database_sig = db_sig,
            outdir = joinpath(workdir, "sourmash_gather"),
            k_size = 21,
            threshold_bp = 250,
            executor = gather_exec
        )
        Test.@test length(gather_exec.jobs) == 1
        Test.@test occursin("sourmash gather", gather_exec.jobs[1].cmd)
        Test.@test occursin("--threshold-bp 250", gather_exec.jobs[1].cmd)
        Test.@test endswith(gather_result.results_csv, "_gather.csv")
        Test.@test endswith(gather_result.results_matches, "_gather_matches.sig")

        cached_search_dir = joinpath(workdir, "cached_search")
        mkpath(cached_search_dir)
        cached_search = joinpath(cached_search_dir, "query_search_results.csv")
        _write_text_file(cached_search, "similarity,name\n0.9,ref\n")
        cached_result = Mycelia.run_sourmash_search(
            query_sig = query_sig,
            database_sig = db_sig,
            outdir = cached_search_dir
        )
        Test.@test cached_result.results_csv == cached_search

        cached_sketch_dir = joinpath(workdir, "cached_sketch")
        mkpath(cached_sketch_dir)
        _write_text_file(joinpath(cached_sketch_dir, "a.sig"), "sig")
        cached_sketch = Mycelia.run_sourmash_sketch(
            input_files = [fasta_a],
            outdir = cached_sketch_dir
        )
        Test.@test cached_sketch.signatures == [joinpath(cached_sketch_dir, "a.sig")]

        gather_csv = _write_text_file(joinpath(workdir, "gather.csv"), "intersect_bp,name\n10,ref\n")
        gather_df = Mycelia.parse_sourmash_gather_output(gather_csv)
        Test.@test DataFrames.nrow(gather_df) == 1

        search_csv = _write_text_file(joinpath(workdir, "search.csv"), "similarity,name\n0.2,ref\n")
        search_df = Mycelia.parse_sourmash_search_output(search_csv)
        Test.@test DataFrames.nrow(search_df) == 1
    end

    Test.@testset "Mash wrapper collection and cache paths" begin
        workdir = mktempdir()

        fasta = _write_text_file(joinpath(workdir, "reads.fa"), ">r1\nACGT\n")
        fastq = _write_text_file(joinpath(workdir, "reads.fastq"), "@r1\nACGT\n+\nIIII\n")
        sketch_a = _write_text_file(joinpath(workdir, "sketch_a.msh"), "msh")
        sketch_b = _write_text_file(joinpath(workdir, "sketch_b.msh"), "msh")

        sketch_exec = Mycelia.CollectExecutor()
        sketch_result = Mycelia.run_mash_sketch(
            input_files = [fasta, fastq],
            outdir = joinpath(workdir, "mash_sketch"),
            r = true,
            min_copies = 2,
            output_prefix = "panel",
            additional_args = ["-I"],
            executor = sketch_exec
        )
        Test.@test length(sketch_exec.jobs) == 2
        Test.@test all(job -> occursin("mash sketch", job.cmd), sketch_exec.jobs)
        Test.@test any(job -> occursin("-r", job.cmd), sketch_exec.jobs)
        Test.@test any(job -> occursin("-m 2", job.cmd), sketch_exec.jobs)
        Test.@test any(job -> occursin("-I", job.cmd), sketch_exec.jobs)
        Test.@test length(sketch_result.sketches) == 2

        paste_exec = Mycelia.CollectExecutor()
        pasted = Mycelia.run_mash_paste(
            out_file = joinpath(workdir, "combined.msh"),
            in_files = [sketch_a, sketch_b],
            executor = paste_exec
        )
        Test.@test length(paste_exec.jobs) == 1
        Test.@test occursin("mash paste", paste_exec.jobs[1].cmd)
        Test.@test endswith(pasted, ".msh")

        dist_exec = Mycelia.CollectExecutor()
        dist_result = Mycelia.run_mash_dist(
            reference = sketch_a,
            query = sketch_b,
            outdir = joinpath(workdir, "mash_dist"),
            threads = 4,
            additional_args = ["-v"],
            executor = dist_exec
        )
        Test.@test length(dist_exec.jobs) == 1
        Test.@test occursin("mash dist", dist_exec.jobs[1].cmd)
        Test.@test occursin("-p 4", dist_exec.jobs[1].cmd)
        Test.@test occursin("-v", dist_exec.jobs[1].cmd)
        Test.@test endswith(dist_result.results_tsv, "_mash_dist.tsv")

        screen_exec = Mycelia.CollectExecutor()
        screen_result = Mycelia.run_mash_screen(
            reference = sketch_a,
            query = [fasta, fastq],
            outdir = joinpath(workdir, "mash_screen"),
            winner_takes_all = false,
            min_identity = 0.95,
            additional_args = ["-v"],
            executor = screen_exec
        )
        Test.@test length(screen_exec.jobs) == 1
        Test.@test occursin("mash screen", screen_exec.jobs[1].cmd)
        Test.@test !occursin(" -w", screen_exec.jobs[1].cmd)
        Test.@test occursin("-i 0.95", screen_exec.jobs[1].cmd)
        Test.@test occursin("-v", screen_exec.jobs[1].cmd)
        Test.@test endswith(screen_result.results_tsv, "_mash_screen.tsv")

        cached_dist_dir = joinpath(workdir, "cached_dist")
        mkpath(cached_dist_dir)
        cached_dist = joinpath(cached_dist_dir, "existing.tsv")
        _write_text_file(cached_dist, "ref\tqry\t0.0\t1\t1/1\n")
        cached_dist_result = Mycelia.run_mash_dist(
            reference = sketch_a,
            query = sketch_b,
            outdir = cached_dist_dir,
            output_tsv = cached_dist
        )
        Test.@test cached_dist_result.results_tsv == cached_dist

        cached_default_sketch_dir = joinpath(workdir, "cached_default_sketch")
        mkpath(cached_default_sketch_dir)
        _write_text_file(joinpath(cached_default_sketch_dir, "reads.msh"), "msh")
        cached_default_sketch = Mycelia.run_mash_sketch(
            input_files = [fasta],
            outdir = cached_default_sketch_dir
        )
        Test.@test cached_default_sketch.sketches == [joinpath(cached_default_sketch_dir, "reads.msh")]

        cached_custom_sketch_dir = joinpath(workdir, "cached_custom_sketch")
        mkpath(cached_custom_sketch_dir)
        _write_text_file(joinpath(cached_custom_sketch_dir, "custom.msh"), "msh")
        cached_custom_sketch = Mycelia.run_mash_sketch(
            input_files = [fasta],
            outdir = cached_custom_sketch_dir,
            output_prefix = "custom"
        )
        Test.@test cached_custom_sketch.sketches == [joinpath(cached_custom_sketch_dir, "custom.msh")]

        winner_exec = Mycelia.CollectExecutor()
        winner_screen = Mycelia.run_mash_screen(
            reference = sketch_a,
            query = fasta,
            outdir = joinpath(workdir, "winner_screen"),
            winner_takes_all = true,
            executor = winner_exec
        )
        Test.@test length(winner_exec.jobs) == 1
        Test.@test occursin("-w", winner_exec.jobs[1].cmd)
        Test.@test endswith(winner_screen.results_tsv, "_mash_screen.tsv")

        cached_screen_dir = joinpath(workdir, "cached_screen")
        mkpath(cached_screen_dir)
        cached_screen = joinpath(cached_screen_dir, "existing.tsv")
        _write_text_file(cached_screen, "0.99\t1/1\t1\t1e-5\tq\tr\n")
        cached_screen_result = Mycelia.run_mash_screen(
            reference = sketch_a,
            query = fasta,
            outdir = cached_screen_dir,
            output_tsv = cached_screen
        )
        Test.@test cached_screen_result.results_tsv == cached_screen
        Test.@test_throws ErrorException Mycelia.run_mash_screen(
            reference = sketch_a,
            query = fasta,
            outdir = joinpath(workdir, "bad_screen"),
            min_identity = 1.5
        )
    end

    Test.@testset "MetaPhlAn wrapper database and parser paths" begin
        workdir = mktempdir()
        db_dir = _make_metaphlan_db(joinpath(workdir, "metaphlan_db"))
        reads = _write_text_file(joinpath(workdir, "reads.sam"), "@HD\tVN:1.0\n")

        default_index, default_explicit = Mycelia.resolve_metaphlan_db_index()
        Test.@test isnothing(default_index)
        Test.@test !default_explicit

        resolved_index, index_explicit = Mycelia.resolve_metaphlan_db_index(db_index = "mpa_vTest")
        Test.@test resolved_index == "mpa_vTest"
        Test.@test index_explicit

        withenv("METAPHLAN_DB_INDEX" => "env_index") do
            env_index, env_explicit = Mycelia.resolve_metaphlan_db_index()
            Test.@test env_index == "env_index"
            Test.@test env_explicit
        end

        Test.@test Mycelia._metaphlan_db_present(db_dir, "mpa_vTest")
        Test.@test !Mycelia._metaphlan_db_present(joinpath(workdir, "missing_db"), "mpa_vTest")
        Test.@test Mycelia._metaphlan_db_present(db_dir)
        Test.@test Mycelia.get_metaphlan_db_path(download = false, db_dir = db_dir, db_index = "mpa_vTest") == db_dir
        Test.@test isnothing(
            Mycelia.get_metaphlan_db_path(
                require = false,
                download = false,
                db_dir = joinpath(workdir, "missing_meta_db"),
                db_index = "mpa_vTest"
            )
        )

        withenv(
            "METAPHLAN_DB_DIR" => "",
            "METAPHLAN_DB_PATH" => "",
            "METAPHLAN_BOWTIE2DB" => db_dir
        ) do
            Test.@test_logs (:warn, r"deprecated") begin
                Test.@test Mycelia.get_metaphlan_db_path(download = false, db_index = "mpa_vTest") == db_dir
            end
        end

        metaphlan_exec = Mycelia.CollectExecutor()
        result = Mycelia.run_metaphlan(
            input_file = reads,
            outdir = joinpath(workdir, "metaphlan_out"),
            input_type = "sam",
            nprocs = 4,
            db_dir = db_dir,
            db_index = "mpa_vTest",
            unknown_estimation = false,
            long_reads = true,
            executor = metaphlan_exec
        )
        Test.@test length(metaphlan_exec.jobs) == 1
        Test.@test occursin("metaphlan", metaphlan_exec.jobs[1].cmd)
        Test.@test occursin("--input_type sam", metaphlan_exec.jobs[1].cmd)
        Test.@test occursin("--index mpa_vTest", metaphlan_exec.jobs[1].cmd)
        Test.@test occursin("--skip_unclassified_estimation", metaphlan_exec.jobs[1].cmd)
        Test.@test occursin("--long_reads", metaphlan_exec.jobs[1].cmd)
        Test.@test result.mapout == joinpath(result.outdir, "reads_mapout.bz2")

        cached_dir = joinpath(workdir, "metaphlan_cached")
        mkpath(cached_dir)
        _write_text_file(joinpath(cached_dir, "reads_profile.txt"), "#profile\n")
        cached_result = Mycelia.run_metaphlan(
            input_file = reads,
            outdir = cached_dir,
            db_dir = db_dir,
            db_index = "mpa_vTest"
        )
        Test.@test cached_result.profile_txt == joinpath(cached_dir, "reads_profile.txt")

        empty_profile = _write_text_file(joinpath(workdir, "empty_profile.txt"), "# comment\n")
        empty_df = Mycelia.parse_metaphlan_profile(empty_profile)
        Test.@test DataFrames.nrow(empty_df) == 0

        mixed_profile = _write_text_file(
            joinpath(workdir, "mixed_profile.txt"),
            "# header\nk__Bacteria\t2\t100.0\t  \ninvalid\nk__Archaea\t\t25.0\tMethanogens\nk__Skip\t1\tnot-a-number\tignored\n"
        )
        mixed_df = Mycelia.parse_metaphlan_profile(mixed_profile)
        Test.@test DataFrames.nrow(mixed_df) == 2
        Test.@test mixed_df.taxid[2] === missing
        Test.@test mixed_df.additional_species[1] === missing
        Test.@test mixed_df.additional_species[2] == "Methanogens"
    end

    Test.@testset "Metabuli wrapper database and parser paths" begin
        workdir = mktempdir()
        reads_1 = _write_text_file(joinpath(workdir, "reads_1.fastq"), "@r1\nACGT\n+\nIIII\n")
        reads_2 = _write_text_file(joinpath(workdir, "reads_2.fastq"), "@r2\nTGCA\n+\nIIII\n")
        db_root = joinpath(workdir, "metabuli_root")
        named_db = _make_metabuli_db(joinpath(db_root, "db_named"))

        Test.@test Mycelia._metabuli_db_present(named_db)
        Test.@test !Mycelia._metabuli_db_present(joinpath(workdir, "missing_metabuli"))
        Test.@test Mycelia._resolve_metabuli_db_dir(db_root, "db_named") == named_db

        withenv("METABULI_DB" => named_db) do
            Test.@test Mycelia.get_metabuli_db_path(download = false) == named_db
        end
        withenv("METABULI_DB" => joinpath(workdir, "bad_metabuli")) do
            Test.@test isnothing(Mycelia.get_metabuli_db_path(require = false, download = false))
        end
        Test.@test Mycelia.get_metabuli_db_path(download = false, db_root = db_root, db_name = "db_named") == named_db
        Test.@test isnothing(
            Mycelia.get_metabuli_db_path(
                require = false,
                download = false,
                db_root = joinpath(workdir, "missing_metabuli_root"),
                db_name = "missing_db"
            )
        )
        Test.@test 1 <= Mycelia._default_metabuli_max_ram_gb() <= 128

        metabuli_exec = Mycelia.CollectExecutor()
        metabuli_result = Mycelia.run_metabuli_classify(
            input_files = [reads_1, reads_2],
            database_path = named_db,
            outdir = joinpath(workdir, "metabuli_out"),
            seq_mode = "3",
            threads = 2,
            min_score = 0.4,
            min_sp_score = 0.6,
            max_ram_gb = 7,
            executor = metabuli_exec
        )
        Test.@test length(metabuli_exec.jobs) == 1
        Test.@test occursin("metabuli classify", metabuli_exec.jobs[1].cmd)
        Test.@test occursin("--seq-mode 3", metabuli_exec.jobs[1].cmd)
        Test.@test occursin("--min-score 0.4", metabuli_exec.jobs[1].cmd)
        Test.@test occursin("--min-sp-score 0.6", metabuli_exec.jobs[1].cmd)
        Test.@test occursin("--max-ram 7", metabuli_exec.jobs[1].cmd)
        Test.@test endswith(metabuli_result.report_file, "_report.tsv")

        taxonomy_dir = joinpath(workdir, "taxonomy")
        mkpath(taxonomy_dir)
        _write_text_file(joinpath(taxonomy_dir, "names.dmp"), "1\t|\troot\t|\n")
        _write_text_file(joinpath(taxonomy_dir, "nodes.dmp"), "1\t|\t1\t|\n")
        reference_fasta = _write_text_file(joinpath(workdir, "reference.fa"), ">ref\nACGT\n")

        build_exec = Mycelia.CollectExecutor()
        build_result = Mycelia.run_metabuli_build_db(
            reference_fasta = reference_fasta,
            taxonomy_dir = taxonomy_dir,
            outdir = joinpath(workdir, "metabuli_build"),
            threads = 3,
            split_num = 64,
            executor = build_exec
        )
        Test.@test length(build_exec.jobs) == 1
        Test.@test occursin("metabuli build", build_exec.jobs[1].cmd)
        Test.@test occursin("--threads 3", build_exec.jobs[1].cmd)
        Test.@test occursin("--split-num 64", build_exec.jobs[1].cmd)
        Test.@test build_result.database_path == joinpath(workdir, "metabuli_build")

        cached_build_dir = joinpath(workdir, "cached_metabuli_build")
        _make_metabuli_db(cached_build_dir)
        cached_build = Mycelia.run_metabuli_build_db(
            reference_fasta = reference_fasta,
            taxonomy_dir = taxonomy_dir,
            outdir = cached_build_dir
        )
        Test.@test cached_build.database_path == cached_build_dir

        partial_report = _write_text_file(
            joinpath(workdir, "partial_report.tsv"),
            "20.0\t4\t2\t  species  \n"
        )
        partial_df = Mycelia.parse_metabuli_report(partial_report)
        Test.@test DataFrames.names(partial_df) == ["percentage", "num_reads", "num_direct_reads", "rank"]

        partial_classifications = _write_text_file(
            joinpath(workdir, "partial_classifications.tsv"),
            "C\tread_1\n"
        )
        partial_class_df = Mycelia.parse_metabuli_classifications(partial_classifications)
        Test.@test DataFrames.ncol(partial_class_df) == 2
    end

    Test.@testset "Legacy and Kraken wrapper fast paths" begin
        workdir = mktempdir()
        assembly = _write_text_file(joinpath(workdir, "assembly.fna"), ">ctg\nACGT\n")
        assembly_gz = joinpath(workdir, "assembly.fna.gz")
        open(assembly_gz, "w") do io
            write(io, "placeholder")
        end

        clamlst_out = _write_text_file(joinpath(workdir, "assembly_clamlst.tsv"), "scheme\tst\n")
        Test.@test Mycelia.run_clamlst(assembly_gz; outdir = workdir) == clamlst_out

        ectyper_dir = joinpath(workdir, "ectyper_out")
        mkpath(ectyper_dir)
        ectyper_out = _write_text_file(joinpath(ectyper_dir, "output.tsv"), "serotype\n")
        Test.@test Mycelia.run_ectyper(assembly; outdir = ectyper_dir) == ectyper_out
        Test.@test_throws ErrorException Mycelia.run_ectyper(
            joinpath(workdir, "missing_assembly.fna");
            outdir = joinpath(workdir, "ectyper_missing")
        )

        ezclermont_out = _write_text_file(joinpath(workdir, "assembly_ezclermont.csv"), "sample,group\n")
        Test.@test Mycelia.run_ezclermont(assembly_gz; outdir = workdir) == ezclermont_out

        kleborate_outdir = joinpath(workdir, "kleborate")
        mkpath(kleborate_outdir)
        _write_text_file(joinpath(kleborate_outdir, "existing_output.txt"), "done\n")
        Test.@test Mycelia.run_kleborate([assembly]; outdir = kleborate_outdir) == kleborate_outdir
        Test.@test Mycelia.run_kleborate(assembly; outdir = kleborate_outdir) == kleborate_outdir
        Test.@test_throws ErrorException Mycelia.run_kleborate(String[]; outdir = joinpath(workdir, "kleborate_empty"))
        Test.@test_throws ErrorException Mycelia.run_kleborate(
            [assembly];
            outdir = joinpath(workdir, "kleborate_invalid"),
            preset = nothing,
            modules = nothing
        )

        kraken_db = _make_kraken2_db(joinpath(workdir, "kraken_db"))
        Test.@test Mycelia._kraken2_db_present(kraken_db)
        Test.@test !Mycelia._kraken2_db_present(joinpath(workdir, "missing_kraken"))
        withenv("KRAKEN2_DB" => kraken_db) do
            Test.@test Mycelia.get_kraken2_db_path(download = false) == kraken_db
        end
        withenv("KRAKEN2_DB" => joinpath(workdir, "bad_kraken")) do
            Test.@test isnothing(Mycelia.get_kraken2_db_path(require = false, download = false))
        end

        kraken_exec = Mycelia.CollectExecutor()
        gz_reads = _write_text_file(joinpath(workdir, "reads.fastq.gz"), "@r1\nACGT\n+\nIIII\n")
        kraken_result = Mycelia.run_kraken2_classify(
            input_files = [gz_reads, gz_reads],
            outdir = joinpath(workdir, "kraken_out"),
            database_path = kraken_db,
            threads = 8,
            confidence = 0.2,
            minimum_hit_groups = 3,
            paired = true,
            gzip_compressed = true,
            report_minimizer_data = true,
            additional_args = ["--quick"],
            executor = kraken_exec
        )
        Test.@test length(kraken_exec.jobs) == 1
        Test.@test occursin("kraken2", kraken_exec.jobs[1].cmd)
        Test.@test occursin("--paired", kraken_exec.jobs[1].cmd)
        Test.@test occursin("--gzip-compressed", kraken_exec.jobs[1].cmd)
        Test.@test occursin("--report-minimizer-data", kraken_exec.jobs[1].cmd)
        Test.@test occursin("--quick", kraken_exec.jobs[1].cmd)
        Test.@test endswith(kraken_result.report_file, "_kraken2_report.txt")

        cached_kraken_dir = joinpath(workdir, "kraken_cached")
        mkpath(cached_kraken_dir)
        _write_text_file(joinpath(cached_kraken_dir, "reads_kraken2_output.txt"), "C\tread\t1\n")
        _write_text_file(joinpath(cached_kraken_dir, "reads_kraken2_report.txt"), "100.0\t1\t1\tS\t1\t name\n")
        cached_kraken = Mycelia.run_kraken2_classify(
            input_files = [gz_reads],
            outdir = cached_kraken_dir,
            database_path = kraken_db
        )
        Test.@test cached_kraken.report_file == joinpath(cached_kraken_dir, "reads_kraken2_report.txt")

        kraken_report = _write_text_file(
            joinpath(workdir, "kraken_report.txt"),
            "100.0\t5\t5\tS\t562\t  Escherichia coli  \n"
        )
        kraken_df = Mycelia.parse_kraken2_report(kraken_report)
        Test.@test DataFrames.nrow(kraken_df) == 1
        Test.@test kraken_df.name[1] == "Escherichia coli"

        Test.@test haskey(Mycelia.list_kraken2_databases(), Mycelia.DEFAULT_KRAKEN2_DB_NAME)
    end

    # ========================================================================
    # Input Validation Tests (no external tools needed)
    # ========================================================================
    Test.@testset "Input Validation" begin
        Test.@testset "Sourmash input validation" begin
            # Test nonexistent file
            Test.@test_throws ErrorException Mycelia.run_sourmash_sketch(
                input_files = ["nonexistent_$(rand(1000:9999)).fq"],
                outdir = tempdir()
            )

            # Test invalid molecule type
            temp_file = tempname()
            touch(temp_file)
            try
                Test.@test_throws ErrorException Mycelia.run_sourmash_sketch(
                    input_files = [temp_file],
                    outdir = tempdir(),
                    molecule = "invalid"
                )
            finally
                rm(temp_file, force = true)
            end
        end

        Test.@testset "MetaPhlAn input validation" begin
            # Test nonexistent file
            Test.@test_throws ErrorException Mycelia.run_metaphlan(
                input_file = "nonexistent_$(rand(1000:9999)).fq",
                outdir = tempdir()
            )

            # Test invalid input type
            temp_file = tempname()
            touch(temp_file)
            try
                Test.@test_throws ErrorException Mycelia.run_metaphlan(
                    input_file = temp_file,
                    outdir = tempdir(),
                    input_type = "invalid_type"
                )
            finally
                rm(temp_file, force = true)
            end
        end

        Test.@testset "Metabuli input validation" begin
            # Test nonexistent input file
            Test.@test_throws ErrorException Mycelia.run_metabuli_classify(
                input_files = ["nonexistent_$(rand(1000:9999)).fq"],
                database_path = tempdir(),
                outdir = tempdir()
            )

            # Test nonexistent database
            temp_file = tempname()
            touch(temp_file)
            try
                Test.@test_throws ErrorException Mycelia.run_metabuli_classify(
                    input_files = [temp_file],
                    database_path = "nonexistent_db_$(rand(1000:9999))",
                    outdir = tempdir()
                )
            finally
                rm(temp_file, force = true)
            end

            # Test invalid seq_mode
            temp_file = tempname()
            touch(temp_file)
            try
                Test.@test_throws ErrorException Mycelia.run_metabuli_classify(
                    input_files = [temp_file],
                    database_path = tempdir(),
                    outdir = tempdir(),
                    seq_mode = "invalid"
                )
            finally
                rm(temp_file, force = true)
            end
        end

        Test.@testset "Parser file validation" begin
            # Test parsing nonexistent files
            Test.@test_throws ErrorException Mycelia.parse_metaphlan_profile(
                "nonexistent_$(rand(1000:9999)).txt"
            )
            Test.@test_throws ErrorException Mycelia.parse_metabuli_report(
                "nonexistent_$(rand(1000:9999)).tsv"
            )
            Test.@test_throws ErrorException Mycelia.parse_metabuli_classifications(
                "nonexistent_$(rand(1000:9999)).tsv"
            )
        end
    end
end
