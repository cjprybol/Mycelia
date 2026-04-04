import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function write_fake_art_conda_runner_default_outbase(path::AbstractString)
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

if [[ "\$tool" == "art_illumina" ]]; then
    outbase=""
    paired=0

    while [[ \$# -gt 0 ]]; do
        case "\$1" in
            --out)
                outbase="\${2:-}"
                shift 2
                ;;
            --paired)
                paired=1
                shift
                ;;
            --samout|--errfree)
                shift
                ;;
            --fcov|--rcount|--seqSys|--len|--mflen|--sdev|--rndSeed|--in)
                shift 2
                ;;
            *)
                shift
                ;;
        esac
    done

    mkdir -p "\$(dirname \"\$outbase\")"
    if [[ \$paired -eq 1 ]]; then
        printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\${outbase}1.fq"
        printf '@r2\\nTGCA\\n+\\n!!!!\\n' > "\${outbase}2.fq"
    else
        printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\${outbase}.fq"
    fi
    printf '@SQ\\tSN:ref\\tLN:4\\n' > "\${outbase}.sam"
    exit 0
fi

printf 'unexpected tool: %s\\n' "\$tool" >&2
exit 1
"""

    write(path, script)
    chmod(path, 0o755)
    return path
end

function with_fake_art_conda_runner_default_outbase(f::Function)
    mktempdir() do dir
        runner_path = Mycelia.CONDA_RUNNER
        backup_path = joinpath(dir, "conda-runner-backup")
        runner_existed = isfile(runner_path)
        if runner_existed
            cp(runner_path, backup_path; force = true)
        else
            mkpath(dirname(runner_path))
        end
        write_fake_art_conda_runner_default_outbase(runner_path)
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

Test.@testset "Simulation Default Outbase Coverage" begin
    with_fake_art_conda_runner_default_outbase() do
        mktempdir() do dir
            fasta = joinpath(dir, "reference.fna")
            write(fasta, ">ref\nACGT\n")

            paired_default = Mycelia.simulate_illumina_reads(
                fasta = fasta,
                read_count = 6,
                quiet = true
            )

            Test.@test paired_default.forward_reads == fasta * ".rcount_6.art1.fq.gz"
            Test.@test paired_default.reverse_reads == fasta * ".rcount_6.art2.fq.gz"
            Test.@test paired_default.sam == fasta * ".rcount_6.art.sam.gz"
            Test.@test paired_default.error_free_sam === nothing
            Test.@test isfile(paired_default.forward_reads)
            Test.@test isfile(paired_default.reverse_reads)
            Test.@test isfile(paired_default.sam)

            single_default = Mycelia.simulate_illumina_reads(
                fasta = fasta,
                coverage = 2,
                paired = false,
                quiet = true
            )

            Test.@test single_default.forward_reads == fasta * ".fcov_2x.art.fq.gz"
            Test.@test single_default.reverse_reads === nothing
            Test.@test single_default.sam == fasta * ".fcov_2x.art.sam.gz"
            Test.@test single_default.error_free_sam === nothing
            Test.@test isfile(single_default.forward_reads)
            Test.@test isfile(single_default.sam)
        end
    end
end
