import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function write_profile_logging_conda_runner(path::AbstractString)
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

if [[ "\$tool" != "art_illumina" ]]; then
    printf 'unexpected tool: %s\\n' "\$tool" >&2
    exit 1
fi

outbase=""
paired=0
seq_sys=""
read_len=""
coverage=""
read_count=""
logfile=""

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
        --seqSys)
            seq_sys="\${2:-}"
            shift 2
            ;;
        --len)
            read_len="\${2:-}"
            shift 2
            ;;
        --fcov)
            coverage="\${2:-}"
            shift 2
            ;;
        --rcount)
            read_count="\${2:-}"
            shift 2
            ;;
        --samout|--errfree)
            shift
            ;;
        --mflen|--sdev|--rndSeed|--in)
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

mkdir -p "\$(dirname \"\$outbase\")"
logfile="\$(dirname \"\$outbase\")/wrapper_calls.tsv"
printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "\$seq_sys" "\$read_len" "\$paired" "\$coverage" "\$read_count" >> "\$logfile"

if [[ \$paired -eq 1 ]]; then
    printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\${outbase}1.fq"
    printf '@r2\\nTGCA\\n+\\n!!!!\\n' > "\${outbase}2.fq"
else
    printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\${outbase}.fq"
fi

printf '@SQ\\tSN:ref\\tLN:4\\n' > "\${outbase}.sam"
exit 0
"""

    write(path, script)
    chmod(path, 0o755)
    return path
end

function with_profile_logging_conda_runner(f::Function)
    mktempdir() do dir
        runner_path = Mycelia.CONDA_RUNNER
        backup_path = joinpath(dir, "conda-runner-backup")
        runner_existed = isfile(runner_path)
        if runner_existed
            cp(runner_path, backup_path; force = true)
        else
            mkpath(dirname(runner_path))
        end
        write_profile_logging_conda_runner(runner_path)
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

function read_wrapper_call(logfile::AbstractString)
    fields = split(chomp(read(logfile, String)), '\t'; keepempty = true, limit = 0)
    return (
        seqSys = fields[1],
        read_length = parse(Int, fields[2]),
        paired = fields[3] == "1",
        coverage = isempty(fields[4]) ? nothing : parse(Float64, fields[4]),
        read_count = isempty(fields[5]) ? nothing : parse(Int, fields[5])
    )
end

Test.@testset "Simulation Profile Wrappers" begin
    wrapper_specs = [
        (fn = :simulate_illumina_ga1_36bp, seqSys = "GA1", read_length = 36, paired = true, use_read_count = false),
        (fn = :simulate_illumina_ga1_44bp, seqSys = "GA1", read_length = 44, paired = true, use_read_count = true),
        (fn = :simulate_illumina_ga2_50bp, seqSys = "GA2", read_length = 50, paired = true, use_read_count = false),
        (fn = :simulate_illumina_ga2_75bp, seqSys = "GA2", read_length = 75, paired = true, use_read_count = true),
        (fn = :simulate_illumina_hs10_100bp, seqSys = "HS10", read_length = 100, paired = true, use_read_count = false),
        (fn = :simulate_illumina_hs20_100bp, seqSys = "HS20", read_length = 100, paired = true, use_read_count = true),
        (fn = :simulate_illumina_hs25_125bp, seqSys = "HS25", read_length = 125, paired = true, use_read_count = false),
        (fn = :simulate_illumina_hs25_150bp, seqSys = "HS25", read_length = 150, paired = true, use_read_count = true),
        (fn = :simulate_illumina_hsxn_150bp, seqSys = "HSXn", read_length = 150, paired = true, use_read_count = false),
        (fn = :simulate_illumina_hsxt_150bp, seqSys = "HSXt", read_length = 150, paired = true, use_read_count = true),
        (fn = :simulate_illumina_msv1_250bp, seqSys = "MSv1", read_length = 250, paired = true, use_read_count = false),
        (fn = :simulate_illumina_msv3_250bp, seqSys = "MSv3", read_length = 250, paired = true, use_read_count = true),
        (fn = :simulate_illumina_mins_50bp, seqSys = "MinS", read_length = 50, paired = false, use_read_count = false),
        (fn = :simulate_illumina_ns50_75bp, seqSys = "NS50", read_length = 75, paired = true, use_read_count = true)
    ]

    with_profile_logging_conda_runner() do
        mktempdir() do dir
            cd(dir) do
                fasta = joinpath(dir, "reference.fna")
                write(fasta, ">ref\nACGTACGT\n")

                for (index, spec) in enumerate(wrapper_specs)
                    logfile = joinpath(dir, "wrapper_calls.tsv")
                    rm(logfile; force = true)
                    outbase = joinpath(dir, "wrapper_$(index)")
                    wrapper = getfield(Mycelia, spec.fn)
                    result = if spec.use_read_count
                        wrapper(fasta = fasta, read_count = 7, outbase = outbase, quiet = true)
                    else
                        wrapper(fasta = fasta, coverage = 3.5, outbase = outbase, quiet = true)
                    end

                    call = read_wrapper_call(logfile)
                    Test.@test call.seqSys == spec.seqSys
                    Test.@test call.read_length == spec.read_length
                    Test.@test call.paired == spec.paired
                    if spec.use_read_count
                        Test.@test call.read_count == 7
                        Test.@test call.coverage === nothing
                    else
                        Test.@test call.coverage == 3.5
                        Test.@test call.read_count === nothing
                    end

                    Test.@test isfile(result.forward_reads)
                    Test.@test isfile(result.sam)
                    if spec.paired
                        Test.@test !isnothing(result.reverse_reads)
                        Test.@test isfile(result.reverse_reads)
                    else
                        Test.@test result.reverse_reads === nothing
                    end
                    Test.@test result.error_free_sam === nothing
                end
            end
        end
    end
end
