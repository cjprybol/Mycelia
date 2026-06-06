import Test
import DataFrames
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function _write_file(path::String, content::String = "")
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

function _make_kraken_db(path::String)
    mkpath(path)
    for filename in ("hash.k2d", "opts.k2d", "taxo.k2d")
        _write_file(joinpath(path, filename), filename)
    end
    return path
end

function _make_tar_archive(source_dir::String, archive_path::String)
    mkpath(dirname(archive_path))
    run(`tar -czf $(archive_path) -C $(source_dir) .`)
    return archive_path
end

function _with_dict_override(f::Function, dict::AbstractDict, key, value)
    had_key = haskey(dict, key)
    old_value = had_key ? dict[key] : nothing
    dict[key] = value

    try
        return f()
    finally
        if had_key
            dict[key] = old_value
        else
            delete!(dict, key)
        end
    end
end

function _with_env(f::Function, bindings::AbstractVector{<:Pair{String, <:Any}})
    saved = Dict{String, Union{Nothing, String}}()
    for (key, value) in bindings
        saved[key] = get(ENV, key, nothing)
        if value === nothing
            if haskey(ENV, key)
                delete!(ENV, key)
            end
        else
            ENV[key] = value
        end
    end

    try
        return f()
    finally
        for (key, value) in saved
            if value === nothing
                if haskey(ENV, key)
                    delete!(ENV, key)
                end
            else
                ENV[key] = value
            end
        end
    end
end

Test.@testset "classification.jl wrapper coverage" begin
    Test.@testset "Sourmash wrappers and parsers" begin
        mktempdir() do tmpdir
            fasta = _write_file(joinpath(tmpdir, "reads.fasta"), ">r1\nACGT\n")
            query_sig = _write_file(joinpath(tmpdir, "query.sig"), "signature\n")
            db_sig = _write_file(joinpath(tmpdir, "db.sig"), "signature\n")
            db_dir = mkpath(joinpath(tmpdir, "sig_db"))
            _write_file(joinpath(db_dir, "ref.sig"), "signature\n")

            collector = Mycelia.CollectExecutor()
            sketch_result = Mycelia.run_sourmash_sketch(
                input_files = [fasta],
                outdir = joinpath(tmpdir, "sketch"),
                k_sizes = [21, 31],
                scaled = 50,
                molecule = "protein",
                singleton = true,
                name = "sample",
                executor = collector,
                site = :scg
            )
            Test.@test sketch_result.signatures == [joinpath(tmpdir, "sketch", "reads.sig")]
            Test.@test length(collector.jobs) == 1
            Test.@test collector.jobs[1].site == :scg
            Test.@test occursin("sourmash sketch protein", collector.jobs[1].cmd)
            Test.@test occursin("--singleton", collector.jobs[1].cmd)
            Test.@test occursin("--name sample", collector.jobs[1].cmd)

            preexisting_sig = _write_file(joinpath(tmpdir, "existing_sketch", "reads.sig"), "done\n")
            existing_collector = Mycelia.CollectExecutor()
            existing_result = Mycelia.run_sourmash_sketch(
                input_files = [fasta],
                outdir = dirname(preexisting_sig),
                executor = existing_collector
            )
            Test.@test existing_result.signatures == [preexisting_sig]
            Test.@test length(existing_collector.jobs) == 1
            Test.@test occursin("if [ ! -f", existing_collector.jobs[1].cmd)
            Test.@test_throws ErrorException Mycelia.run_sourmash_sketch(
                input_files = [joinpath(tmpdir, "missing.fasta")],
                outdir = joinpath(tmpdir, "missing")
            )
            Test.@test_throws ErrorException Mycelia.run_sourmash_sketch(
                input_files = [fasta],
                outdir = joinpath(tmpdir, "invalid_molecule"),
                molecule = "rna"
            )

            search_collector = Mycelia.CollectExecutor()
            search_result = Mycelia.run_sourmash_search(
                query_sig = query_sig,
                database_sig = db_sig,
                outdir = joinpath(tmpdir, "search"),
                threshold = 0.2,
                k_size = 51,
                best_only = true,
                num_results = 7,
                executor = search_collector
            )
            Test.@test endswith(search_result.results_csv, "_search_results.csv")
            Test.@test length(search_collector.jobs) == 1
            Test.@test occursin("sourmash search", search_collector.jobs[1].cmd)
            Test.@test occursin("--best-only", search_collector.jobs[1].cmd)
            Test.@test occursin("-n 7", search_collector.jobs[1].cmd)

            existing_search = _write_file(
                joinpath(tmpdir, "search_done", "query_search_results.csv"),
                "similarity,name,filename,md5\n0.9,ref,ref.sig,abcd\n"
            )
            search_skip = Mycelia.run_sourmash_search(
                query_sig = query_sig,
                database_sig = db_dir,
                outdir = dirname(existing_search),
                executor = Mycelia.CollectExecutor()
            )
            Test.@test search_skip.results_csv == existing_search

            gather_collector = Mycelia.CollectExecutor()
            gather_result = Mycelia.run_sourmash_gather(
                query_sig = query_sig,
                database_sig = db_dir,
                outdir = joinpath(tmpdir, "gather"),
                k_size = 31,
                threshold_bp = 1000,
                executor = gather_collector
            )
            Test.@test endswith(gather_result.results_csv, "_gather.csv")
            Test.@test endswith(gather_result.results_matches, "_gather_matches.sig")
            Test.@test length(gather_collector.jobs) == 1
            Test.@test occursin("sourmash gather", gather_collector.jobs[1].cmd)
            Test.@test occursin("--save-matches", gather_collector.jobs[1].cmd)

            gather_csv = _write_file(
                joinpath(tmpdir, "gather.csv"),
                "intersect_bp,f_orig_query,f_match,name\n100,0.4,0.5,ref_a\n"
            )
            gather_df = Mycelia.parse_sourmash_gather_output(gather_csv)
            Test.@test DataFrames.nrow(gather_df) == 1
            Test.@test gather_df.name[1] == "ref_a"

            search_csv = _write_file(
                joinpath(tmpdir, "search.csv"),
                "similarity,name,filename,md5\n0.95,ref_b,ref.sig,deadbeef\n"
            )
            search_df = Mycelia.parse_sourmash_search_output(search_csv)
            Test.@test DataFrames.nrow(search_df) == 1
            Test.@test search_df.similarity[1] == 0.95
            Test.@test_throws ErrorException Mycelia.parse_sourmash_gather_output(joinpath(tmpdir, "missing.csv"))
            Test.@test_throws ErrorException Mycelia.parse_sourmash_search_output(joinpath(tmpdir, "missing_search.csv"))
        end
    end

    Test.@testset "Mash wrappers and parsers" begin
        mktempdir() do tmpdir
            fasta = _write_file(joinpath(tmpdir, "reads.fasta"), ">r1\nACGT\n")
            fastq = _write_file(joinpath(tmpdir, "reads.fastq"), "@r1\nACGT\n+\n!!!!\n")
            ref_sketch = _write_file(joinpath(tmpdir, "ref.msh"), "sketch\n")
            query_sketch = _write_file(joinpath(tmpdir, "query.msh"), "sketch\n")

            single_sketch_path = _write_file(joinpath(tmpdir, "mash_single", "custom.msh"), "msh\n")
            sketch_result = Mycelia.run_mash_sketch(
                input_files = [fasta],
                outdir = dirname(single_sketch_path),
                output_prefix = "custom"
            )
            Test.@test sketch_result.sketches == [single_sketch_path]

            multi_outdir = joinpath(tmpdir, "mash_multi")
            first_multi = _write_file(joinpath(multi_outdir, "combo_reads.msh"), "msh\n")
            second_multi = _write_file(joinpath(multi_outdir, "combo_reads_2.msh"), "msh\n")
            fastq2 = _write_file(joinpath(tmpdir, "reads_2.fastq"), "@r2\nTGCA\n+\n!!!!\n")
            multi_result = Mycelia.run_mash_sketch(
                input_files = [fastq, fastq2],
                outdir = multi_outdir,
                output_prefix = "combo"
            )
            Test.@test multi_result.sketches == [first_multi, second_multi]
            Test.@test_throws ErrorException Mycelia.run_mash_sketch(
                input_files = [fasta],
                outdir = joinpath(tmpdir, "bad_k"),
                k = 0
            )
            Test.@test_throws ErrorException Mycelia.run_mash_sketch(
                input_files = [fastq],
                outdir = joinpath(tmpdir, "bad_min_copies"),
                min_copies = 2
            )
            Test.@test_throws ErrorException Mycelia.run_mash_sketch(
                input_files = [fastq],
                outdir = joinpath(tmpdir, "bad_threads"),
                r = true,
                min_copies = 1,
                threads = 0
            )

            pasted = _write_file(joinpath(tmpdir, "db", "combined.msh"), "msh\n")
            paste_result = Mycelia.run_mash_paste(
                out_file = joinpath(tmpdir, "db", "combined"),
                in_files = [ref_sketch, query_sketch]
            )
            Test.@test paste_result == pasted
            Test.@test_throws ErrorException Mycelia.run_mash_paste(
                out_file = joinpath(tmpdir, "db", "missing"),
                in_files = [joinpath(tmpdir, "absent.msh")]
            )

            dist_path = _write_file(
                joinpath(tmpdir, "dist", "custom.tsv"),
                "ref.msh\tquery.msh\t0.01\t0\t100/100\n"
            )
            dist_result = Mycelia.run_mash_dist(
                reference = ref_sketch,
                query = query_sketch,
                outdir = dirname(dist_path),
                output_tsv = dist_path
            )
            Test.@test dist_result.results_tsv == dist_path
            dist_df = Mycelia.parse_mash_dist_output(dist_path)
            Test.@test DataFrames.nrow(dist_df) == 1
            Test.@test dist_df.distance[1] == 0.01
            Test.@test_throws ErrorException Mycelia.run_mash_dist(
                reference = ref_sketch,
                query = query_sketch,
                outdir = joinpath(tmpdir, "bad_dist"),
                threads = 0
            )

            screen_path = _write_file(
                joinpath(tmpdir, "screen", "custom.tsv"),
                "#comment\n0.99\t100/100\t1\t0\treads.fastq\tref.msh\n"
            )
            screen_result = Mycelia.run_mash_screen(
                reference = ref_sketch,
                query = [fastq],
                outdir = dirname(screen_path),
                output_tsv = screen_path,
                winner_takes_all = false
            )
            Test.@test screen_result.results_tsv == screen_path
            screen_df = Mycelia.parse_mash_screen_output(screen_path)
            Test.@test DataFrames.nrow(screen_df) == 1
            Test.@test screen_df.identity[1] == 0.99
            Test.@test_throws ErrorException Mycelia.run_mash_screen(
                reference = ref_sketch,
                query = fastq,
                outdir = joinpath(tmpdir, "bad_screen"),
                min_identity = 1.5
            )
        end
    end

    Test.@testset "MetaPhlAn helpers and wrapper" begin
        mktempdir() do tmpdir
            db_dir = mkpath(joinpath(tmpdir, "metaphlan_db"))
            _write_file(joinpath(db_dir, "mpa_vJan21.pkl"), "db\n")
            input_fastq = _write_file(joinpath(tmpdir, "reads.fastq"), "@r1\nACGT\n+\n!!!!\n")

            Test.@test Mycelia.resolve_metaphlan_db_index(db_index = "custom") == ("custom", true)
            _with_env(["METAPHLAN_DB_INDEX" => "from_env"]) do
                Test.@test Mycelia.resolve_metaphlan_db_index() == ("from_env", true)
            end
            _with_env(["METAPHLAN_DB_INDEX" => nothing]) do
                Test.@test Mycelia.resolve_metaphlan_db_index() == (nothing, false)
            end

            Test.@test Mycelia._metaphlan_db_present(db_dir)
            Test.@test Mycelia._metaphlan_db_present(db_dir, "mpa_vJan21")
            Test.@test !Mycelia._metaphlan_db_present(db_dir, "other")
            Test.@test !Mycelia._metaphlan_db_present(joinpath(tmpdir, "missing_db"))

            Test.@test Mycelia.get_metaphlan_db_path(
                require = true,
                db_dir = db_dir,
                db_index = "mpa_vJan21",
                download = false
            ) == db_dir
            _with_env([
                "METAPHLAN_DB_DIR" => db_dir,
                "METAPHLAN_DB_PATH" => nothing,
                "METAPHLAN_BOWTIE2DB" => nothing
            ]) do
                Test.@test Mycelia.get_metaphlan_db_path(download = false) == db_dir
            end
            Test.@test isnothing(Mycelia.get_metaphlan_db_path(
                require = false,
                db_dir = joinpath(tmpdir, "absent_db"),
                download = false
            ))

            collector = Mycelia.CollectExecutor()
            run_result = Mycelia.run_metaphlan(
                input_file = input_fastq,
                outdir = joinpath(tmpdir, "metaphlan_run"),
                input_type = "fastq",
                nprocs = 3,
                db_dir = db_dir,
                db_index = "mpa_vJan21",
                unknown_estimation = false,
                stat_q = 0.3,
                long_reads = true,
                executor = collector,
                site = :scg
            )
            Test.@test endswith(run_result.profile_txt, "_profile.txt")
            Test.@test length(collector.jobs) == 1
            Test.@test collector.jobs[1].cpus_per_task == 3
            Test.@test occursin("--db_dir", collector.jobs[1].cmd)
            Test.@test occursin("--index mpa_vJan21", collector.jobs[1].cmd)
            Test.@test occursin("--skip_unclassified_estimation", collector.jobs[1].cmd)
            Test.@test occursin("--long_reads", collector.jobs[1].cmd)

            existing_profile = _write_file(
                joinpath(tmpdir, "metaphlan_done", "reads_profile.txt"),
                "#header\n"
            )
            _write_file(joinpath(tmpdir, "metaphlan_done", "reads_mapout.bz2"), "")
            skip_collector = Mycelia.CollectExecutor()
            skip_result = Mycelia.run_metaphlan(
                input_file = input_fastq,
                outdir = dirname(existing_profile),
                db_dir = db_dir,
                executor = skip_collector
            )
            Test.@test skip_result.profile_txt == existing_profile
            Test.@test isempty(skip_collector.jobs)
            Test.@test_throws ErrorException Mycelia.run_metaphlan(
                input_file = input_fastq,
                outdir = joinpath(tmpdir, "bad_type"),
                db_dir = db_dir,
                input_type = "bamish"
            )

            profile_file = _write_file(
                joinpath(tmpdir, "profile.txt"),
                "#comment\nk__Bacteria\t2\t100.0\t\ninvalid\trow\nk__Bacteria|p__Firmicutes\t\t33.5\tBacillus\n"
            )
            parsed_profile = Mycelia.parse_metaphlan_profile(profile_file)
            Test.@test DataFrames.nrow(parsed_profile) == 2
            Test.@test isequal(parsed_profile.taxid[2], missing)
            Test.@test parsed_profile.additional_species[2] == "Bacillus"

            empty_profile = _write_file(joinpath(tmpdir, "empty_profile.txt"), "#comment only\n")
            empty_df = Mycelia.parse_metaphlan_profile(empty_profile)
            Test.@test DataFrames.nrow(empty_df) == 0
        end
    end

    Test.@testset "Metabuli helpers and wrappers" begin
        mktempdir() do tmpdir
            direct_db = mkpath(joinpath(tmpdir, "metabuli_root", "refseq"))
            _write_file(joinpath(direct_db, "db.info"), "db\n")
            alt_db = mkpath(joinpath(tmpdir, "metabuli_root", "refseq_alt"))
            _write_file(joinpath(alt_db, "info"), "db\n")
            input_fastq = _write_file(joinpath(tmpdir, "reads.fastq"), "@r1\nACGT\n+\n!!!!\n")

            Test.@test Mycelia._metabuli_db_present(direct_db)
            Test.@test Mycelia._metabuli_db_present(alt_db)
            Test.@test !Mycelia._metabuli_db_present(joinpath(tmpdir, "missing_metabuli"))
            Test.@test !Mycelia._metabuli_db_present(mkpath(joinpath(tmpdir, "empty_metabuli")))
            Test.@test Mycelia._resolve_metabuli_db_dir(dirname(direct_db), "refseq") == direct_db

            named_root = mkpath(joinpath(tmpdir, "named_root"))
            named_match = mkpath(joinpath(named_root, "virus_refseq"))
            _write_file(joinpath(named_match, "db.info"), "db\n")
            other_match = mkpath(joinpath(named_root, "other"))
            _write_file(joinpath(other_match, "db.info"), "db\n")
            Test.@test Mycelia._resolve_metabuli_db_dir(named_root, "virus") == named_match

            singleton_root = mkpath(joinpath(tmpdir, "singleton_root"))
            lone_db = mkpath(joinpath(singleton_root, "only_db"))
            _write_file(joinpath(lone_db, "info"), "db\n")
            Test.@test Mycelia._resolve_metabuli_db_dir(singleton_root, "missing_name") == lone_db
            Test.@test_throws ErrorException Mycelia._resolve_metabuli_db_dir(mkpath(joinpath(tmpdir, "empty_root")), "nothing")

            _with_env([
                "METABULI_DB" => direct_db,
                "METABULI_DB_PATH" => nothing,
                "METABULI_DB_ROOT" => nothing,
                "METABULI_DB_NAME" => nothing
            ]) do
                Test.@test Mycelia.get_metabuli_db_path(download = false) == direct_db
            end
            _with_env([
                "METABULI_DB" => nothing,
                "METABULI_DB_PATH" => nothing,
                "METABULI_DB_ROOT" => dirname(direct_db),
                "METABULI_DB_NAME" => basename(direct_db)
            ]) do
                Test.@test Mycelia.get_metabuli_db_path(download = false) == direct_db
            end
            Test.@test isnothing(Mycelia.get_metabuli_db_path(
                require = false,
                db_root = joinpath(tmpdir, "absent_root"),
                db_name = "missing",
                download = false
            ))
            _with_env([
                "METABULI_DB" => joinpath(tmpdir, "missing_env_db"),
                "METABULI_DB_PATH" => nothing
            ]) do
                Test.@test isnothing(Mycelia.get_metabuli_db_path(require = false, download = false))
            end
            Test.@test_throws ErrorException Mycelia.get_metabuli_db_path(
                require = true,
                db_root = joinpath(tmpdir, "absent_root"),
                db_name = "missing",
                download = false
            )

            staged_metabuli = mkpath(joinpath(tmpdir, "metabuli_archive_src", "refseq_release"))
            _write_file(joinpath(staged_metabuli, "db.info"), "db\n")
            metabuli_archive = _make_tar_archive(
                dirname(staged_metabuli),
                joinpath(tmpdir, "metabuli_refseq_release.tar.gz")
            )
            _with_dict_override(
                Mycelia.METABULI_DB_URLS,
                "refseq_release",
                "file://" * metabuli_archive
            ) do
                downloaded_db = Mycelia.download_metabuli_db(
                    db_name = "refseq_release",
                    db_root = joinpath(tmpdir, "downloaded_metabuli"),
                    force = true
                )
                Test.@test downloaded_db == joinpath(tmpdir, "downloaded_metabuli", "refseq_release")
                Test.@test isfile(joinpath(downloaded_db, "db.info"))
                Test.@test !isfile(joinpath(tmpdir, "downloaded_metabuli", basename(metabuli_archive)))
            end
            Test.@test_throws ErrorException Mycelia.download_metabuli_db(
                db_name = "not_a_real_metabuli_db",
                db_root = joinpath(tmpdir, "bad_metabuli_download")
            )

            staged_metabuli_get = mkpath(joinpath(tmpdir, "metabuli_get_src", "refseq_release"))
            _write_file(joinpath(staged_metabuli_get, "db.info"), "db\n")
            metabuli_get_archive = _make_tar_archive(
                dirname(staged_metabuli_get),
                joinpath(tmpdir, "metabuli_get_refseq_release.tar.gz")
            )
            _with_dict_override(
                Mycelia.METABULI_DB_URLS,
                "refseq_release",
                "file://" * metabuli_get_archive
            ) do
                resolved_download = Mycelia.get_metabuli_db_path(
                    db_root = joinpath(tmpdir, "metabuli_get_download"),
                    db_name = "refseq_release",
                    download = true
                )
                Test.@test resolved_download == joinpath(tmpdir, "metabuli_get_download", "refseq_release")
                Test.@test isfile(joinpath(resolved_download, "db.info"))
            end
            Test.@test Mycelia._default_metabuli_max_ram_gb() >= 1
            Test.@test Mycelia._default_metabuli_max_ram_gb() <= 128

            classify_collector = Mycelia.CollectExecutor()
            classify_result = Mycelia.run_metabuli_classify(
                input_files = [input_fastq],
                outdir = joinpath(tmpdir, "metabuli_run"),
                database_path = direct_db,
                seq_mode = "2",
                threads = 4,
                min_score = 0.5,
                min_sp_score = 0.7,
                max_ram_gb = 8,
                executor = classify_collector,
                site = :scg
            )
            Test.@test endswith(classify_result.report_file, "_report.tsv")
            Test.@test endswith(classify_result.classifications_file, "_classifications.tsv")
            Test.@test length(classify_collector.jobs) == 1
            Test.@test classify_collector.jobs[1].cpus_per_task == 4
            Test.@test occursin("--seq-mode 2", classify_collector.jobs[1].cmd)
            Test.@test occursin("--min-score 0.5", classify_collector.jobs[1].cmd)
            Test.@test occursin("--min-sp-score 0.7", classify_collector.jobs[1].cmd)
            Test.@test occursin("--max-ram 8", classify_collector.jobs[1].cmd)

            existing_report = _write_file(
                joinpath(tmpdir, "metabuli_done", "metabuli_classification_report.tsv"),
                "100\t10\t10\tS\t1\tSpecies\n"
            )
            _write_file(
                joinpath(tmpdir, "metabuli_done", "metabuli_classification_classifications.tsv"),
                "C\tread1\t1\n"
            )
            skip_collect = Mycelia.CollectExecutor()
            skip_classify = Mycelia.run_metabuli_classify(
                input_files = [input_fastq],
                outdir = dirname(existing_report),
                database_path = direct_db,
                executor = skip_collect
            )
            Test.@test skip_classify.report_file == existing_report
            Test.@test isempty(skip_collect.jobs)
            Test.@test_throws ErrorException Mycelia.run_metabuli_classify(
                input_files = [input_fastq],
                outdir = joinpath(tmpdir, "bad_seq_mode"),
                database_path = direct_db,
                seq_mode = "9"
            )
            Test.@test_throws ErrorException Mycelia.run_metabuli_classify(
                input_files = [input_fastq],
                outdir = joinpath(tmpdir, "bad_ram"),
                database_path = direct_db,
                max_ram_gb = 0
            )

            taxonomy_dir = mkpath(joinpath(tmpdir, "taxonomy"))
            _write_file(joinpath(taxonomy_dir, "names.dmp"), "1\t|\troot\t|\n")
            _write_file(joinpath(taxonomy_dir, "nodes.dmp"), "1\t|\t1\t|\tno rank\t|\n")
            reference_fasta = _write_file(joinpath(tmpdir, "reference.fna"), ">ref\nACGT\n")
            build_collector = Mycelia.CollectExecutor()
            build_result = Mycelia.run_metabuli_build_db(
                reference_fasta = reference_fasta,
                taxonomy_dir = taxonomy_dir,
                outdir = joinpath(tmpdir, "metabuli_build"),
                threads = 2,
                split_num = 128,
                executor = build_collector,
                site = :scg
            )
            Test.@test build_result.database_path == joinpath(tmpdir, "metabuli_build")
            Test.@test length(build_collector.jobs) == 1
            Test.@test build_collector.jobs[1].cpus_per_task == 2
            Test.@test occursin("--split-num 128", build_collector.jobs[1].cmd)

            report_file = _write_file(
                joinpath(tmpdir, "metabuli_report.tsv"),
                "100.0\t10\t4\tS\t562\t  Escherichia coli  \n"
            )
            report_df = Mycelia.parse_metabuli_report(report_file)
            Test.@test DataFrames.nrow(report_df) == 1
            Test.@test report_df.name[1] == "Escherichia coli"

            partial_report = _write_file(
                joinpath(tmpdir, "metabuli_partial_report.tsv"),
                "100.0\t10\t4\n"
            )
            partial_df = Mycelia.parse_metabuli_report(partial_report)
            Test.@test :percentage in DataFrames.propertynames(partial_df)
            Test.@test :num_direct_reads in DataFrames.propertynames(partial_df)

            classifications_file = _write_file(
                joinpath(tmpdir, "metabuli_classifications.tsv"),
                "C\tread1\t562\tfoo\nU\tread2\t0\tbar\n"
            )
            classifications_df = Mycelia.parse_metabuli_classifications(classifications_file)
            Test.@test DataFrames.nrow(classifications_df) == 2
            Test.@test :classified in DataFrames.propertynames(classifications_df)
            Test.@test classifications_df.read_id[1] == "read1"
        end
    end

    Test.@testset "Kleborate, EzClermont, ECTyper, claMLST, and Kraken2 helpers" begin
        mktempdir() do tmpdir
            assembly = _write_file(joinpath(tmpdir, "assembly.fna"), ">contig\nACGT\n")

            kleborate_outdir = mkpath(joinpath(tmpdir, "kleborate"))
            _write_file(joinpath(kleborate_outdir, "sample_output.txt"), "done\n")
            Test.@test Mycelia.run_kleborate(String[]; outdir = kleborate_outdir) == kleborate_outdir
            Test.@test Mycelia.run_kleborate(assembly; outdir = kleborate_outdir) == kleborate_outdir
            Test.@test_throws ErrorException Mycelia.run_kleborate(String[]; outdir = joinpath(tmpdir, "empty_kleborate"))
            redirect_stdout(devnull) do
                Test.@test length(Mycelia.list_kraken2_databases()) > 0
            end

            ectyper_dir = mkpath(joinpath(tmpdir, "ectyper_output"))
            ectyper_output = _write_file(joinpath(ectyper_dir, "output.tsv"), "done\n")
            Test.@test Mycelia.run_ectyper(assembly; outdir = ectyper_dir) == ectyper_output

            ez_outdir = mkpath(joinpath(tmpdir, "ez"))
            ez_output = _write_file(joinpath(ez_outdir, "assembly_ezclermont.csv"), "sample,group\n")
            Test.@test Mycelia.run_ezclermont(assembly; outdir = ez_outdir) == ez_output

            clamlst_outdir = mkpath(joinpath(tmpdir, "clamlst"))
            clamlst_output = _write_file(joinpath(clamlst_outdir, "assembly_clamlst.tsv"), "sample\tst\n")
            Test.@test Mycelia.run_clamlst(assembly; outdir = clamlst_outdir) == clamlst_output

            kraken_db = _make_kraken_db(joinpath(tmpdir, "kraken_db"))
            Test.@test Mycelia._kraken2_db_present(kraken_db)
            Test.@test !Mycelia._kraken2_db_present(joinpath(tmpdir, "missing_kraken"))

            _with_env([
                "KRAKEN2_DB" => kraken_db,
                "KRAKEN2_DB_PATH" => nothing,
                "KRAKEN2_DB_ROOT" => nothing,
                "KRAKEN2_DB_NAME" => nothing
            ]) do
                Test.@test Mycelia.get_kraken2_db_path(download = false) == kraken_db
            end
            _with_env([
                "KRAKEN2_DB" => nothing,
                "KRAKEN2_DB_PATH" => nothing,
                "KRAKEN2_DB_ROOT" => tmpdir,
                "KRAKEN2_DB_NAME" => "kraken_db"
            ]) do
                Test.@test Mycelia.get_kraken2_db_path(download = false) == kraken_db
            end
            Test.@test isnothing(Mycelia.get_kraken2_db_path(
                require = false,
                db_root = joinpath(tmpdir, "absent_root"),
                db_name = "missing",
                download = false
            ))
            _with_env([
                "KRAKEN2_DB" => joinpath(tmpdir, "missing_kraken_env"),
                "KRAKEN2_DB_PATH" => nothing
            ]) do
                Test.@test isnothing(Mycelia.get_kraken2_db_path(require = false, download = false))
            end
            Test.@test_throws ErrorException Mycelia.get_kraken2_db_path(
                require = true,
                db_root = joinpath(tmpdir, "absent_root"),
                db_name = "missing",
                download = false
            )

            staged_kraken = mkpath(joinpath(tmpdir, "kraken_archive_src"))
            _make_kraken_db(staged_kraken)
            kraken_archive = _make_tar_archive(
                staged_kraken,
                joinpath(tmpdir, "kraken_standard.tar.gz")
            )
            _with_dict_override(
                Mycelia.KRAKEN2_DB_URLS,
                "standard",
                "file://" * kraken_archive
            ) do
                downloaded_kraken = Mycelia.download_kraken2_db(
                    db_name = "standard",
                    db_root = joinpath(tmpdir, "downloaded_kraken"),
                    force = true
                )
                Test.@test downloaded_kraken == joinpath(tmpdir, "downloaded_kraken", "standard")
                Test.@test Mycelia._kraken2_db_present(downloaded_kraken)
                Test.@test !isfile(joinpath(tmpdir, "downloaded_kraken", basename(kraken_archive)))
            end
            Test.@test_throws ErrorException Mycelia.download_kraken2_db(
                db_name = "not_a_real_kraken_db",
                db_root = joinpath(tmpdir, "bad_kraken_download")
            )

            staged_kraken_get = mkpath(joinpath(tmpdir, "kraken_get_src"))
            _make_kraken_db(staged_kraken_get)
            kraken_get_archive = _make_tar_archive(
                staged_kraken_get,
                joinpath(tmpdir, "kraken_get_standard.tar.gz")
            )
            _with_dict_override(
                Mycelia.KRAKEN2_DB_URLS,
                "standard",
                "file://" * kraken_get_archive
            ) do
                resolved_kraken = Mycelia.get_kraken2_db_path(
                    db_root = joinpath(tmpdir, "kraken_get_download"),
                    db_name = "standard",
                    download = true
                )
                Test.@test resolved_kraken == joinpath(tmpdir, "kraken_get_download", "standard")
                Test.@test Mycelia._kraken2_db_present(resolved_kraken)
            end

            reads_fastq = _write_file(joinpath(tmpdir, "reads.fastq"), "@r1\nACGT\n+\n!!!!\n")
            kraken_outdir = mkpath(joinpath(tmpdir, "kraken_run"))
            output_file = _write_file(joinpath(kraken_outdir, "reads_kraken2_output.txt"), "C\tread1\t562\tA:1\n")
            report_file = _write_file(
                joinpath(kraken_outdir, "reads_kraken2_report.txt"),
                "100.00\t10\t10\tS\t562\t  Escherichia coli\n"
            )
            kraken_result = Mycelia.run_kraken2_classify(
                input_files = [reads_fastq],
                outdir = kraken_outdir,
                database_path = kraken_db
            )
            Test.@test kraken_result.output_file == output_file
            Test.@test kraken_result.report_file == report_file
            Test.@test_throws ErrorException Mycelia.run_kraken2_classify(
                input_files = String[],
                outdir = joinpath(tmpdir, "kraken_empty"),
                database_path = kraken_db
            )
            Test.@test_throws ErrorException Mycelia.run_kraken2_classify(
                input_files = [reads_fastq],
                outdir = joinpath(tmpdir, "kraken_bad_db"),
                database_path = joinpath(tmpdir, "bad_db")
            )

            kraken_df = Mycelia.parse_kraken2_report(report_file)
            Test.@test DataFrames.nrow(kraken_df) == 1
            Test.@test kraken_df[1, "name"] == "Escherichia coli"
        end
    end
end
