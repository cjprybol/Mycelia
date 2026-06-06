import Test
import Mycelia
import DataFrames
import CodecZlib

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function _write_file(path::String, content::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

function _write_gzip_file(path::String, content::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        gzip = CodecZlib.GzipCompressorStream(io)
        write(gzip, content)
        close(gzip)
    end
    return path
end

function _mock_tool_script()::String
    return """#!/usr/bin/env bash
set -euo pipefail

tool="\$(basename "\$0")"

if [ -n "\${MYCELIA_MOCK_TOOL_LOG:-}" ]; then
  printf '%s %s\n' "\${tool}" "\$*" >> "\${MYCELIA_MOCK_TOOL_LOG}"
fi

case "\${tool}" in
  mash)
    sub="\${1:-}"
    shift || true
    case "\${sub}" in
      sketch)
        prefix=""
        while [ "\$#" -gt 0 ]; do
          case "\$1" in
            -o)
              prefix="\$2"
              shift 2
              ;;
            *)
              shift
              ;;
          esac
        done
        touch "\${prefix}.msh"
        ;;
      paste)
        output_prefix="\${1:-}"
        touch "\${output_prefix}.msh"
        ;;
      dist)
        printf 'ref.msh\tquery.msh\t0.01\t1e-06\t10/1000\n'
        ;;
      screen)
        printf '0.95\t15/1000\t2\t1e-06\treads.fastq\treference.msh\n'
        ;;
      *)
        echo "unsupported mash subcommand: \${sub}" >&2
        exit 98
        ;;
    esac
    ;;
  metaphlan)
    if [ "\${1:-}" = "--install" ]; then
        shift
        db_dir=""
        db_index="mock_metaphlan"
        while [ "\$#" -gt 0 ]; do
          case "\$1" in
            --db_dir)
              db_dir="\$2"
              shift 2
              ;;
            --index)
              db_index="\$2"
              shift 2
              ;;
            *)
              shift
              ;;
          esac
        done
        mkdir -p "\${db_dir}"
        touch "\${db_dir}/\${db_index}.pkl"
    else
        profile=""
        mapout=""
        while [ "\$#" -gt 0 ]; do
          case "\$1" in
            -o)
              profile="\$2"
              shift 2
              ;;
            --mapout)
              mapout="\$2"
              shift 2
              ;;
            *)
              shift
              ;;
          esac
        done
        mkdir -p "\$(dirname "\${profile}")"
        cat > "\${profile}" <<'EOF'
#mpa_vMock_CHOCOPhlAnSGB_202604
#SampleID	Metaphlan_Analysis
k__Bacteria	2	100.0
EOF
        printf 'mock mapout\n' > "\${mapout}"
    fi
    ;;
  metabuli)
    sub="\${1:-}"
    shift || true
    case "\${sub}" in
      build)
        reference_fasta="\${1:-}"
        taxonomy_dir="\${2:-}"
        outdir="\${3:-}"
        mkdir -p "\${outdir}"
        printf 'mock db\n' > "\${outdir}/db.info"
        ;;
      classify)
        input_file="\${1:-}"
        database_path="\${2:-}"
        outdir="\${3:-}"
        job_id="\${4:-metabuli_classification}"
        mkdir -p "\${outdir}"
        cat > "\${outdir}/\${job_id}_report.tsv" <<'EOF'
99.0	99	99	species	562	Escherichia coli
EOF
        cat > "\${outdir}/\${job_id}_classifications.tsv" <<'EOF'
C	read_001	562
EOF
        ;;
      *)
        echo "unsupported metabuli subcommand: \${sub}" >&2
        exit 99
        ;;
    esac
    ;;
  claMLST)
    sub="\${1:-}"
    shift || true
    case "\${sub}" in
      import)
        db_path="\${1:-}"
        mkdir -p "\$(dirname "\${db_path}")"
        printf 'mock clamlst db\n' > "\${db_path}"
        ;;
      search)
        db_path="\${1:-}"
        input_path="\${2:-}"
        printf 'sample\tst\n%s\t42\n' "\$(basename "\${input_path}")"
        ;;
      *)
        echo "unsupported claMLST subcommand: \${sub}" >&2
        exit 100
        ;;
    esac
    ;;
  ectyper)
    outdir=""
    while [ "\$#" -gt 0 ]; do
      case "\$1" in
        -o)
          outdir="\$2"
          shift 2
          ;;
        *)
          shift
          ;;
      esac
    done
    mkdir -p "\${outdir}"
    cat > "\${outdir}/output.tsv" <<'EOF'
Name	Species
sample	Escherichia coli
EOF
    ;;
  ezclermont)
    input_path="\${1:-}"
    printf 'sample,phylogroup\n%s,A\n' "\$(basename "\${input_path}")"
    ;;
  kleborate)
    if [ "\${1:-}" = "--list_modules" ]; then
      printf 'module_a\nmodule_b\n'
    else
      outdir=""
      while [ "\$#" -gt 0 ]; do
        case "\$1" in
          --outdir)
            outdir="\$2"
            shift 2
            ;;
          *)
            shift
            ;;
        esac
      done
      mkdir -p "\${outdir}"
      printf 'mock kleborate results\n' > "\${outdir}/mock_output.txt"
    fi
    ;;
  kraken2)
    output=""
    report=""
    while [ "\$#" -gt 0 ]; do
      case "\$1" in
        --output)
          output="\$2"
          shift 2
          ;;
        --report)
          report="\$2"
          shift 2
          ;;
        *)
          shift
          ;;
      esac
    done
    mkdir -p "\$(dirname "\${output}")"
    printf 'C\tread_001\t562\t100\n' > "\${output}"
    printf '100.0\t1\t1\tS\t562\tEscherichia coli\n' > "\${report}"
    ;;
  *)
    echo "unsupported mock tool: \${tool}" >&2
    exit 101
    ;;
esac
"""
end

function _read_log(path::String)::String
    return isfile(path) ? read(path, String) : ""
end

function _conda_envs_root()::String
    return joinpath(dirname(dirname(Mycelia.CONDA_RUNNER)), "envs")
end

function _install_mock_env(
        workspace::String,
        backups::Vector{Tuple{String, Union{Nothing, String}}},
        env_name::String,
        executable_name::String,
        script::String)
    env_dir = joinpath(_conda_envs_root(), env_name)
    backup_dir = nothing
    if isdir(env_dir)
        backup_dir = joinpath(workspace, "backup_" * env_name)
        mv(env_dir, backup_dir; force = true)
    end
    push!(backups, (env_dir, backup_dir))

    mkpath(joinpath(env_dir, "bin"))
    mkpath(joinpath(env_dir, "conda-meta"))
    touch(joinpath(env_dir, "conda-meta", "history"))
    wrapper_path = _write_file(joinpath(env_dir, "bin", executable_name), script)
    run(`chmod +x $wrapper_path`)
    return env_dir
end

function _with_mock_conda_envs(f::Function)
    workspace = mktempdir()
    script = _mock_tool_script()
    log_path = joinpath(workspace, "mock_tools.log")
    backups = Tuple{String, Union{Nothing, String}}[]

    for (env_name, executable_name) in (
            ("mash", "mash"),
            ("metaphlan", "metaphlan"),
            ("metabuli", "metabuli"),
            ("pymlst", "claMLST"),
            ("ectyper", "ectyper"),
            ("ezclermont", "ezclermont"),
            ("kleborate", "kleborate"),
            ("kraken2", "kraken2"))
        _install_mock_env(workspace, backups, env_name, executable_name, script)
    end

    previous_log = get(ENV, "MYCELIA_MOCK_TOOL_LOG", nothing)
    ENV["MYCELIA_MOCK_TOOL_LOG"] = log_path

    try
        return f(workspace, log_path)
    finally
        if previous_log === nothing
            delete!(ENV, "MYCELIA_MOCK_TOOL_LOG")
        else
            ENV["MYCELIA_MOCK_TOOL_LOG"] = previous_log
        end

        for (env_dir, backup_dir) in reverse(backups)
            rm(env_dir; recursive = true, force = true)
            if backup_dir !== nothing
                mv(backup_dir, env_dir; force = true)
            end
        end
    end
end

Test.@testset "classification wrapper unit tests" begin
    _with_mock_conda_envs() do workspace, log_path
        Test.@testset "mash wrappers" begin
            ref_fasta = _write_file(joinpath(workspace, "reference.fasta"), ">ref\nACGTACGT\n")
            reads_fastq = _write_file(
                joinpath(workspace, "reads.fastq"),
                "@read1\nACGT\n+\n!!!!\n"
            )

            sketch_result = Mycelia.run_mash_sketch(
                input_files = [ref_fasta, reads_fastq],
                outdir = joinpath(workspace, "mash_sketch"),
                output_prefix = "combo",
                k = 31,
                s = 500,
                r = true,
                min_copies = 2,
                threads = 3,
                additional_args = ["--seed", "7"]
            )

            Test.@test length(sketch_result.sketches) == 2
            Test.@test all(isfile, sketch_result.sketches)
            Test.@test occursin("mash sketch", _read_log(log_path))
            Test.@test occursin("-m 2", _read_log(log_path))

            pasted = Mycelia.run_mash_paste(
                out_file = joinpath(workspace, "mash_db", "combined.msh"),
                in_files = sketch_result.sketches
            )
            Test.@test isfile(pasted)

            dist_result = Mycelia.run_mash_dist(
                reference = sketch_result.sketches[1],
                query = sketch_result.sketches[2],
                outdir = joinpath(workspace, "mash_dist"),
                threads = 2,
                additional_args = ["-v"]
            )
            Test.@test isfile(dist_result.results_tsv)
            dist_table = Mycelia.parse_mash_dist_output(dist_result.results_tsv)
            Test.@test DataFrames.nrow(dist_table) == 1
            Test.@test dist_table.reference[1] == "ref.msh"

            screen_result = Mycelia.run_mash_screen(
                reference = sketch_result.sketches[1],
                query = [reads_fastq],
                outdir = joinpath(workspace, "mash_screen"),
                winner_takes_all = false,
                threads = 4,
                min_identity = 0.91,
                additional_args = ["-v"]
            )
            Test.@test isfile(screen_result.results_tsv)
            screen_table = Mycelia.parse_mash_screen_output(screen_result.results_tsv)
            Test.@test DataFrames.nrow(screen_table) == 1
            Test.@test screen_table.reference[1] == "reference.msh"

            Test.@test_throws ErrorException Mycelia.run_mash_sketch(
                input_files = [ref_fasta],
                outdir = joinpath(workspace, "invalid_sketch"),
                min_copies = 2
            )
            Test.@test_throws ErrorException Mycelia.run_mash_screen(
                reference = sketch_result.sketches[1],
                query = reads_fastq,
                outdir = joinpath(workspace, "invalid_screen"),
                min_identity = 1.5
            )
        end

        Test.@testset "metaphlan and metabuli helpers" begin
            withenv("METAPHLAN_DB_INDEX" => nothing) do
                Test.@test Mycelia.resolve_metaphlan_db_index() == (nothing, false)
            end
            withenv("METAPHLAN_DB_INDEX" => "env_index") do
                Test.@test Mycelia.resolve_metaphlan_db_index() == ("env_index", true)
            end
            Test.@test Mycelia.resolve_metaphlan_db_index(db_index = "explicit") ==
                        ("explicit", true)

            metaphlan_dir = joinpath(workspace, "metaphlan_db")
            downloaded_dir = Mycelia.download_metaphlan_db(
                db_dir = metaphlan_dir,
                db_index = "mock_v1",
                force = true
            )
            Test.@test downloaded_dir == metaphlan_dir
            Test.@test Mycelia._metaphlan_db_present(metaphlan_dir, "mock_v1")
            Test.@test Mycelia.get_metaphlan_db_path(
                db_dir = metaphlan_dir,
                db_index = "mock_v1",
                download = false
            ) == metaphlan_dir
            Test.@test isnothing(Mycelia.get_metaphlan_db_path(
                require = false,
                db_dir = joinpath(workspace, "missing_metaphlan"),
                download = false
            ))

            input_fastq = _write_file(joinpath(workspace, "metaphlan.fastq"), "@r1\nACGT\n+\n!!!!\n")
            metaphlan_result = Mycelia.run_metaphlan(
                input_file = input_fastq,
                outdir = joinpath(workspace, "metaphlan_run"),
                db_dir = metaphlan_dir,
                db_index = "mock_v1",
                input_type = "fastq",
                nprocs = 4,
                unknown_estimation = false,
                long_reads = true
            )
            Test.@test isfile(metaphlan_result.profile_txt)
            Test.@test isfile(metaphlan_result.mapout)
            metaphlan_df = Mycelia.parse_metaphlan_profile(metaphlan_result.profile_txt)
            Test.@test DataFrames.nrow(metaphlan_df) == 1
            Test.@test metaphlan_df.relative_abundance[1] == 100.0

            collect_executor = Mycelia.CollectExecutor()
            queued_metaphlan = Mycelia.run_metaphlan(
                input_file = input_fastq,
                outdir = joinpath(workspace, "metaphlan_collect"),
                db_dir = metaphlan_dir,
                db_index = "mock_v1",
                executor = collect_executor,
                site = :scg
            )
            Test.@test queued_metaphlan.profile_txt ==
                        joinpath(workspace, "metaphlan_collect", "metaphlan_profile.txt")
            Test.@test length(collect_executor.jobs) == 1
            Test.@test occursin("metaphlan", collect_executor.jobs[1].cmd)

            taxonomy_dir = joinpath(workspace, "taxonomy")
            _write_file(joinpath(taxonomy_dir, "names.dmp"), "1\t|\troot\t|\n")
            _write_file(joinpath(taxonomy_dir, "nodes.dmp"), "1\t|\t1\t|\n")
            reference_fasta = _write_file(joinpath(workspace, "reference.fa"), ">ref\nACGT\n")
            metabuli_db = Mycelia.run_metabuli_build_db(
                reference_fasta = reference_fasta,
                taxonomy_dir = taxonomy_dir,
                outdir = joinpath(workspace, "metabuli_db"),
                threads = 2,
                split_num = 128
            )
            Test.@test Mycelia._metabuli_db_present(metabuli_db.database_path)
            Test.@test Mycelia._resolve_metabuli_db_dir(dirname(metabuli_db.database_path), basename(metabuli_db.database_path)) ==
                        metabuli_db.database_path
            Test.@test Mycelia.get_metabuli_db_path(
                db_root = dirname(metabuli_db.database_path),
                db_name = basename(metabuli_db.database_path),
                download = false
            ) == metabuli_db.database_path
            Test.@test isnothing(Mycelia.get_metabuli_db_path(
                require = false,
                db_root = joinpath(workspace, "missing_metabuli_root"),
                db_name = "missing",
                download = false
            ))

            metabuli_input = _write_file(
                joinpath(workspace, "metabuli.fastq"),
                "@read1\nACGT\n+\n!!!!\n"
            )
            metabuli_result = Mycelia.run_metabuli_classify(
                input_files = [metabuli_input],
                database_path = metabuli_db.database_path,
                outdir = joinpath(workspace, "metabuli_run"),
                seq_mode = "1",
                threads = 3,
                min_score = 0.5,
                min_sp_score = 0.7,
                max_ram_gb = 4
            )
            Test.@test isfile(metabuli_result.report_file)
            Test.@test isfile(metabuli_result.classifications_file)
            metabuli_report = Mycelia.parse_metabuli_report(metabuli_result.report_file)
            metabuli_classifications = Mycelia.parse_metabuli_classifications(
                metabuli_result.classifications_file
            )
            Test.@test metabuli_report.taxid[1] == 562
            Test.@test metabuli_classifications.read_id[1] == "read_001"

            collect_executor = Mycelia.CollectExecutor()
            queued = Mycelia.run_metabuli_build_db(
                reference_fasta = reference_fasta,
                taxonomy_dir = taxonomy_dir,
                outdir = joinpath(workspace, "metabuli_db_collect"),
                threads = 5,
                executor = collect_executor,
                site = :scg
            )
            Test.@test queued.database_path == joinpath(workspace, "metabuli_db_collect")
            Test.@test length(collect_executor.jobs) == 1
            Test.@test occursin("metabuli build", collect_executor.jobs[1].cmd)
        end

        Test.@testset "single-genome classification wrappers" begin
            genome_fasta = _write_file(joinpath(workspace, "sample.fasta"), ">chr1\nACGTACGT\n")
            genome_gz = _write_gzip_file(joinpath(workspace, "sample.fasta.gz"), ">chr1\nACGTACGT\n")
            clamlst_db = _write_file(
                joinpath(workspace, "clamlst", "claMLSTDB"),
                "preexisting db\n"
            )

            clamlst_output = Mycelia.run_clamlst(
                genome_gz;
                db_path = clamlst_db,
                outdir = joinpath(workspace, "clamlst_out"),
                species = "Escherichia coli",
                force_db_update = true
            )
            Test.@test isfile(clamlst_output)
            Test.@test occursin("sample.fasta", read(clamlst_output, String))
            Test.@test occursin("claMLST import", _read_log(log_path))
            Test.@test occursin("--force", _read_log(log_path))

            ectyper_result = Mycelia.run_ectyper(
                genome_fasta;
                outdir = joinpath(workspace, "ectyper_out"),
                verify_species = false,
                threads = 2
            )
            Test.@test ectyper_result == joinpath(workspace, "ectyper_out")
            Test.@test isfile(joinpath(ectyper_result, "output.tsv"))
            ectyper_cached = Mycelia.run_ectyper(
                genome_fasta;
                outdir = joinpath(workspace, "ectyper_out")
            )
            Test.@test ectyper_cached == joinpath(workspace, "ectyper_out", "output.tsv")

            ezclermont_result = Mycelia.run_ezclermont(
                genome_gz;
                outdir = joinpath(workspace, "ezclermont_out")
            )
            Test.@test isfile(ezclermont_result)
            Test.@test occursin("phylogroup", read(ezclermont_result, String))
            Test.@test !isfile(joinpath(workspace, "ezclermont_out", "tmp_ez_sample.fasta"))

            Mycelia.list_kleborate_modules()
            Test.@test occursin("kleborate --list_modules", _read_log(log_path))

            kleborate_outdir = Mycelia.run_kleborate(
                [genome_fasta];
                outdir = joinpath(workspace, "kleborate_out"),
                modules = ["amr", "virulence"],
                trim_headers = true
            )
            Test.@test kleborate_outdir == joinpath(workspace, "kleborate_out")
            Test.@test isfile(joinpath(kleborate_outdir, "mock_output.txt"))
            Test.@test Mycelia.run_kleborate(genome_fasta; outdir = kleborate_outdir) ==
                        kleborate_outdir
            Test.@test_throws ErrorException Mycelia.run_kleborate(
                String[];
                outdir = joinpath(workspace, "kleborate_empty")
            )
        end

        Test.@testset "kraken2 helpers and wrapper" begin
            kraken_db = joinpath(workspace, "kraken_db")
            mkpath(kraken_db)
            Test.@test !Mycelia._kraken2_db_present(kraken_db)
            for marker in ("hash.k2d", "opts.k2d", "taxo.k2d")
                _write_file(joinpath(kraken_db, marker), marker)
            end
            Test.@test Mycelia._kraken2_db_present(kraken_db)

            withenv("KRAKEN2_DB" => kraken_db) do
                Test.@test Mycelia.get_kraken2_db_path(download = false) == kraken_db
            end
            Test.@test Mycelia.get_kraken2_db_path(
                db_root = workspace,
                db_name = "kraken_db",
                download = false
            ) == kraken_db
            Test.@test isnothing(Mycelia.get_kraken2_db_path(
                require = false,
                db_root = joinpath(workspace, "missing_kraken_root"),
                db_name = "missing",
                download = false
            ))

            reads_1 = _write_gzip_file(joinpath(workspace, "reads_1.fastq.gz"), "@r1\nACGT\n+\n!!!!\n")
            reads_2 = _write_gzip_file(joinpath(workspace, "reads_2.fastq.gz"), "@r2\nTGCA\n+\n!!!!\n")
            kraken_result = Mycelia.run_kraken2_classify(
                input_files = [reads_1, reads_2],
                outdir = joinpath(workspace, "kraken_out"),
                database_path = kraken_db,
                threads = 2,
                confidence = 0.25,
                minimum_hit_groups = 3,
                paired = true,
                gzip_compressed = nothing,
                report_minimizer_data = true,
                additional_args = ["--quick"]
            )
            Test.@test isfile(kraken_result.output_file)
            Test.@test isfile(kraken_result.report_file)
            kraken_report = Mycelia.parse_kraken2_report(kraken_result.report_file)
            Test.@test DataFrames.nrow(kraken_report) == 1
            Test.@test kraken_report.name[1] == "Escherichia coli"
        end
    end
end
