import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function write_fake_sra_conda_runner(path::AbstractString)
    script = """
#!/usr/bin/env bash
set -euo pipefail

if [[ "\${1:-}" == "run" ]]; then
    shift
fi

if [[ "\${1:-}" == "--live-stream" ]]; then
    shift
fi

if [[ "\${1:-}" == "-n" ]]; then
    shift 2
fi

tool="\${1:-}"
if [[ -n "\$tool" ]]; then
    shift
fi

if [[ "\$tool" == "prefetch" ]]; then
    srr=""
    outdir=""

    while [[ \$# -gt 0 ]]; do
        case "\$1" in
            -O)
                outdir="\${2:-}"
                shift 2
                ;;
            *)
                if [[ -z "\$srr" ]]; then
                    srr="\$1"
                fi
                shift
                ;;
        esac
    done

    if [[ "\$srr" == *FAIL* ]]; then
        printf 'prefetch failure for %s\\n' "\$srr" >&2
        exit 1
    fi

    mkdir -p "\$outdir/\$srr"
    printf 'sra:%s\\n' "\$srr" > "\$outdir/\$srr/\$srr.sra"
    exit 0
fi

if [[ "\$tool" == "fasterq-dump" ]]; then
    outdir=""
    target=""

    while [[ \$# -gt 0 ]]; do
        case "\$1" in
            --outdir)
                outdir="\${2:-}"
                shift 2
                ;;
            --mem|--threads)
                shift 2
                ;;
            --split-3|--skip-technical)
                shift
                ;;
            *)
                target="\$1"
                shift
                ;;
        esac
    done

    srr="\$(basename "\$target")"

    if [[ "\$srr" == *FAIL* ]]; then
        printf 'fasterq failure for %s\\n' "\$srr" >&2
        exit 1
    fi

    mkdir -p "\$outdir"

    if [[ "\$srr" == *EMPTY* ]]; then
        exit 0
    fi

    if [[ "\$srr" == *PAIRED* ]]; then
        printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\$outdir/\${srr}_1.fastq"
        printf '@r2\\nTGCA\\n+\\n!!!!\\n' > "\$outdir/\${srr}_2.fastq"
    else
        printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\$outdir/\${srr}.fastq"
    fi
    exit 0
fi

printf 'unexpected tool: %s\\n' "\$tool" >&2
exit 1
"""

    write(path, script)
    chmod(path, 0o755)
    return path
end

function with_fake_sra_conda_runner(f::Function)
    mktempdir() do dir
        runner_path = Mycelia.CONDA_RUNNER
        backup_path = joinpath(dir, "conda-runner-backup")
        runner_existed = isfile(runner_path)
        if runner_existed
            cp(runner_path, backup_path; force = true)
        else
            mkpath(dirname(runner_path))
        end
        write_fake_sra_conda_runner(runner_path)
        try
            return f()
        finally
            if runner_existed
                mv(backup_path, runner_path; force = true)
            else
                rm(runner_path; force = true)
            end
        end
    end
end

Test.@testset "Reference Database SRA Helpers" begin
    with_fake_sra_conda_runner() do
        mktempdir() do dir
            paired_prefetch = Mycelia.prefetch(SRR = "SRR_PAIRED", outdir = dir)
            Test.@test paired_prefetch.directory == joinpath(dir, "SRR_PAIRED")
            Test.@test paired_prefetch.archive == joinpath(dir, "SRR_PAIRED", "SRR_PAIRED.sra")
            Test.@test isfile(paired_prefetch.archive)

            cached_prefetch = Mycelia.prefetch(SRR = "SRR_PAIRED", outdir = dir)
            Test.@test cached_prefetch == paired_prefetch

            paired_fastq = Mycelia.fasterq_dump(outdir = dir, srr_identifier = "SRR_PAIRED")
            Test.@test paired_fastq.forward_reads == joinpath(dir, "SRR_PAIRED", "SRR_PAIRED_1.fastq.gz")
            Test.@test paired_fastq.reverse_reads == joinpath(dir, "SRR_PAIRED", "SRR_PAIRED_2.fastq.gz")
            Test.@test ismissing(paired_fastq.unpaired_reads)
            Test.@test isfile(paired_fastq.forward_reads)
            Test.@test isfile(paired_fastq.reverse_reads)

            cached_paired_fastq = Mycelia.fasterq_dump(outdir = dir, srr_identifier = "SRR_PAIRED")
            Test.@test isequal(cached_paired_fastq, paired_fastq)

            single_fastq = Mycelia.fasterq_dump(outdir = dir, srr_identifier = "SRR_SINGLE")
            Test.@test ismissing(single_fastq.forward_reads)
            Test.@test ismissing(single_fastq.reverse_reads)
            Test.@test single_fastq.unpaired_reads == joinpath(dir, "SRR_SINGLE", "SRR_SINGLE.fastq.gz")
            Test.@test isfile(single_fastq.unpaired_reads)
        end
    end

    with_fake_sra_conda_runner() do
        mktempdir() do dir
            paired_download = Mycelia.download_sra_data("SRR_PAIRED"; outdir = dir)
            Test.@test paired_download.srr_id == "SRR_PAIRED"
            Test.@test paired_download.outdir == joinpath(dir, "SRR_PAIRED")
            Test.@test paired_download.is_paired
            Test.@test paired_download.files == [
                joinpath(dir, "SRR_PAIRED", "SRR_PAIRED_1.fastq.gz"),
                joinpath(dir, "SRR_PAIRED", "SRR_PAIRED_2.fastq.gz")
            ]

            single_download = Mycelia.download_sra_data("SRR_SINGLE"; outdir = dir)
            Test.@test single_download.srr_id == "SRR_SINGLE"
            Test.@test single_download.outdir == joinpath(dir, "SRR_SINGLE")
            Test.@test !single_download.is_paired
            Test.@test single_download.files == [joinpath(dir, "SRR_SINGLE", "SRR_SINGLE.fastq.gz")]

            Test.@test_throws ErrorException Mycelia.download_sra_data("SRR_EMPTY"; outdir = dir)
        end
    end

    with_fake_sra_conda_runner() do
        mktempdir() do dir
            prefetch_results = Mycelia.prefetch_sra_runs(
                ["SRR_PAIRED", "SRR_FAIL"];
                outdir = dir,
                max_parallel = 1
            )
            Test.@test length(prefetch_results) == 2
            Test.@test prefetch_results[1].srr_id == "SRR_PAIRED"
            Test.@test prefetch_results[1].success
            Test.@test prefetch_results[1].result.archive ==
                       joinpath(dir, "SRR_PAIRED", "SRR_PAIRED.sra")
            Test.@test prefetch_results[2].srr_id == "SRR_FAIL"
            Test.@test !prefetch_results[2].success
            Test.@test !isnothing(prefetch_results[2].error)

            fasterq_results = Mycelia.fasterq_dump_parallel(
                ["SRR_SINGLE", "SRR_FAIL"];
                outdir = dir,
                max_parallel = 1
            )
            Test.@test length(fasterq_results) == 2
            Test.@test fasterq_results[1].srr_id == "SRR_SINGLE"
            Test.@test fasterq_results[1].success
            Test.@test fasterq_results[2].srr_id == "SRR_FAIL"
            Test.@test !fasterq_results[2].success
            Test.@test !isnothing(fasterq_results[2].error)
        end
    end
end
