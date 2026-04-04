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

# Check if external tool tests should run
const RUN_ALL = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
const RUN_EXTERNAL = RUN_ALL || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
# const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

function write_text_file(path::AbstractString, contents::AbstractString)
    open(path, "w") do io
        write(io, contents)
    end
    return path
end

function write_gzip_text_file(path::AbstractString, contents::AbstractString)
    open(path, "w") do io
        gzip_io = CodecZlib.GzipCompressorStream(io)
        try
            write(gzip_io, contents)
        finally
            close(gzip_io)
        end
    end
    return path
end

function create_metabuli_db_dir(path::AbstractString; marker::AbstractString = "db.info")
    mkpath(path)
    touch(joinpath(path, marker))
    return path
end

function create_kraken2_db_dir(path::AbstractString)
    mkpath(path)
    for marker in ("hash.k2d", "opts.k2d", "taxo.k2d")
        touch(joinpath(path, marker))
    end
    return path
end

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
    # Wrapper Helper Coverage Tests (no external tools needed)
    # ========================================================================
    Test.@testset "Wrapper Helper Coverage" begin
        Test.@testset "MetaPhlAn database resolution" begin
            Test.@test Mycelia.resolve_metaphlan_db_index(db_index = "mpa_vTest") ==
                        ("mpa_vTest", true)

            Base.withenv("METAPHLAN_DB_INDEX" => "mpa_env_index") do
                Test.@test Mycelia.resolve_metaphlan_db_index() == ("mpa_env_index", true)
            end

            mktempdir() do dir
                write_text_file(joinpath(dir, "mpa_vTest.pkl"), "pickle-placeholder")

                Test.@test Mycelia._metaphlan_db_present(dir)
                Test.@test Mycelia._metaphlan_db_present(dir, "mpa_vTest")
                Test.@test !Mycelia._metaphlan_db_present(dir, "missing")

                Base.withenv("METAPHLAN_DB_DIR" => dir, "METAPHLAN_DB_INDEX" => "mpa_vTest") do
                    Test.@test Mycelia.get_metaphlan_db_path(download = false) == dir
                end

                missing_dir = joinpath(dir, "missing")
                Base.withenv("METAPHLAN_DB_DIR" => missing_dir) do
                    Test.@test Mycelia.get_metaphlan_db_path(
                        require = false,
                        download = false
                    ) === nothing
                    Test.@test_throws ErrorException Mycelia.get_metaphlan_db_path(
                        download = false
                    )
                end
            end
        end

        Test.@testset "Metabuli database resolution" begin
            mktempdir() do dir
                direct_db = create_metabuli_db_dir(joinpath(dir, "refseq_virus"))
                fallback_db = create_metabuli_db_dir(joinpath(dir, "metabuli_refseq_bacteria"); marker = "info")

                Test.@test Mycelia._metabuli_db_present(direct_db)
                Test.@test !Mycelia._metabuli_db_present(joinpath(dir, "missing"))
                Test.@test Mycelia._resolve_metabuli_db_dir(dir, "refseq_virus") == direct_db
                Test.@test Mycelia._resolve_metabuli_db_dir(dir, "bacteria") == fallback_db

                Base.withenv("METABULI_DB" => direct_db) do
                    Test.@test Mycelia.get_metabuli_db_path(download = false) == direct_db
                end

                Base.withenv(
                    "METABULI_DB" => "",
                    "METABULI_DB_PATH" => "",
                    "METABULI_DB_ROOT" => dir,
                    "METABULI_DB_NAME" => "refseq_virus"
                ) do
                    Test.@test Mycelia.get_metabuli_db_path(download = false) == direct_db
                end

                missing_db = joinpath(dir, "missing_db")
                Base.withenv("METABULI_DB" => missing_db) do
                    Test.@test Mycelia.get_metabuli_db_path(
                        require = false,
                        download = false
                    ) === nothing
                    Test.@test_throws ErrorException Mycelia.get_metabuli_db_path(
                        download = false
                    )
                end
            end

            mktempdir() do dir
                Test.@test_throws ErrorException Mycelia.download_metabuli_db(
                    db_name = "not_a_real_metabuli_db",
                    db_root = dir
                )
            end
        end

        Test.@testset "Cached wrapper outputs" begin
            mktempdir() do dir
                genome_file = joinpath(dir, "sample.fasta.gz")
                write_gzip_text_file(genome_file, ">contig_1\nACGTACGT\n")
                outdir = joinpath(dir, "clamlst")
                mkpath(outdir)
                expected_output = write_text_file(
                    joinpath(outdir, "sample_clamlst.tsv"),
                    "sample\t42\n"
                )

                Test.@test Mycelia.run_clamlst(
                    genome_file;
                    outdir = outdir,
                    db_path = joinpath(dir, "claMLSTDB")
                ) == expected_output
            end

            mktempdir() do dir
                genome_file = write_text_file(
                    joinpath(dir, "sample.fasta"),
                    ">contig_1\nACGTACGT\n"
                )
                cached_outdir = joinpath(dir, "ectyper_cached")
                mkpath(cached_outdir)
                expected_output = write_text_file(
                    joinpath(cached_outdir, "output.tsv"),
                    "sample\tO1:H7\n"
                )

                Test.@test Mycelia.run_ectyper(genome_file; outdir = cached_outdir) ==
                            expected_output

                Test.@test_throws ErrorException Mycelia.run_ectyper(
                    joinpath(dir, "missing.fasta");
                    outdir = joinpath(dir, "ectyper_missing")
                )
            end

            mktempdir() do dir
                genome_file = joinpath(dir, "isolate.fna.gz")
                write_gzip_text_file(genome_file, ">contig_1\nACGTACGT\n")
                outdir = joinpath(dir, "ezclermont")
                mkpath(outdir)
                expected_output = write_text_file(
                    joinpath(outdir, "isolate_ezclermont.csv"),
                    "isolate,A\n"
                )

                Test.@test Mycelia.run_ezclermont(genome_file; outdir = outdir) ==
                            expected_output
            end

            mktempdir() do dir
                assembly_file = write_text_file(
                    joinpath(dir, "assembly.fasta"),
                    ">contig_1\nACGTACGT\n"
                )
                cached_outdir = joinpath(dir, "kleborate_cached")
                mkpath(cached_outdir)
                write_text_file(joinpath(cached_outdir, "kpsc_output.txt"), "complete\n")

                Test.@test Mycelia.run_kleborate([assembly_file]; outdir = cached_outdir) ==
                            cached_outdir
                Test.@test Mycelia.run_kleborate(assembly_file; outdir = cached_outdir) ==
                            cached_outdir
                Test.@test_throws ErrorException Mycelia.run_kleborate(
                    String[];
                    outdir = joinpath(dir, "kleborate_empty")
                )
                Test.@test_throws ErrorException Mycelia.run_kleborate(
                    [assembly_file];
                    outdir = joinpath(dir, "kleborate_invalid"),
                    preset = nothing,
                    modules = String[]
                )
            end
        end

        Test.@testset "Kraken2 helper resolution and cached execution" begin
            databases = Mycelia.list_kraken2_databases()
            Test.@test databases isa Dict
            Test.@test !isempty(databases)

            mktempdir() do dir
                db_dir = create_kraken2_db_dir(joinpath(dir, "standard"))
                input_file = joinpath(dir, "reads.fastq.gz")
                outdir = joinpath(dir, "kraken2")
                mkpath(outdir)

                write_gzip_text_file(input_file, "@read1\nACGT\n+\nIIII\n")
                output_file = write_text_file(
                    joinpath(outdir, "reads_kraken2_output.txt"),
                    "C\tread1\t562\t4\tA:1\n"
                )
                report_file = write_text_file(
                    joinpath(outdir, "reads_kraken2_report.txt"),
                    "100.00\t1\t1\tS\t562\t  Escherichia coli\n"
                )

                Test.@test Mycelia._kraken2_db_present(db_dir)

                Base.withenv("KRAKEN2_DB" => db_dir) do
                    Test.@test Mycelia.get_kraken2_db_path(download = false) == db_dir
                end

                Base.withenv(
                    "KRAKEN2_DB" => "",
                    "KRAKEN2_DB_PATH" => "",
                    "KRAKEN2_DB_ROOT" => dir,
                    "KRAKEN2_DB_NAME" => "standard"
                ) do
                    Test.@test Mycelia.get_kraken2_db_path(download = false) == db_dir
                end

                Base.withenv("KRAKEN2_DB_PATH" => joinpath(dir, "missing")) do
                    Test.@test Mycelia.get_kraken2_db_path(
                        require = false,
                        download = false
                    ) === nothing
                    Test.@test_throws ErrorException Mycelia.get_kraken2_db_path(
                        download = false
                    )
                end

                result = Mycelia.run_kraken2_classify(
                    input_files = [input_file],
                    outdir = outdir,
                    database_path = db_dir
                )

                Test.@test result.output_file == output_file
                Test.@test result.report_file == report_file

                df = Mycelia.parse_kraken2_report(report_file)
                Test.@test df isa DataFrames.DataFrame
                Test.@test DataFrames.nrow(df) == 1
                Test.@test df.name[1] == "Escherichia coli"
                Test.@test df.taxid[1] == 562
            end

            mktempdir() do dir
                Test.@test_throws ErrorException Mycelia.download_kraken2_db(
                    db_name = "not_a_real_kraken2_db",
                    db_root = dir
                )
            end
        end
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

        Test.@testset "Kraken2 input validation" begin
            Test.@test_throws ErrorException Mycelia.run_kraken2_classify(
                input_files = ["nonexistent_$(rand(1000:9999)).fq"],
                database_path = tempdir(),
                outdir = tempdir()
            )

            temp_file = tempname() * ".fq"
            touch(temp_file)
            try
                Test.@test_throws ErrorException Mycelia.run_kraken2_classify(
                    input_files = [temp_file],
                    database_path = "nonexistent_db_$(rand(1000:9999))",
                    outdir = tempdir()
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
            Test.@test_throws ErrorException Mycelia.parse_kraken2_report(
                "nonexistent_$(rand(1000:9999)).tsv"
            )
        end
    end
end
