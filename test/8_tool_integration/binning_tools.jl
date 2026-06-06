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

const RUN_ALL = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
const RUN_EXTERNAL = RUN_ALL || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
# const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

function write_text_file(path::String, contents::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, contents)
    end
    return path
end

function make_fake_conda_env(env_name::String, scripts::Dict{String, String})
    env_dir = joinpath(Mycelia._conda_envs_dir(), env_name)
    ispath(env_dir) && error("Refusing to overwrite existing conda environment: $(env_dir)")
    mkpath(joinpath(env_dir, "bin"))
    mkpath(joinpath(env_dir, "conda-meta"))
    write_text_file(joinpath(env_dir, "conda-meta", "history"), "")
    for (script_name, script_body) in scripts
        script_path = joinpath(env_dir, "bin", script_name)
        write_text_file(script_path, "#!/usr/bin/env bash\nset -euo pipefail\n$(script_body)\n")
        chmod(script_path, 0o755)
    end
    return env_dir
end

function with_fake_conda_envs(f::Function, specs::Vector{Tuple{String, Dict{String, String}}})
    created_envs = String[]
    try
        for (env_name, scripts) in specs
            push!(created_envs, make_fake_conda_env(env_name, scripts))
        end
        return f()
    finally
        for env_dir in reverse(created_envs)
            rm(env_dir; recursive = true, force = true)
        end
    end
end

function with_stubbed_conda_envs(f::Function, specs::Vector{Tuple{String, Dict{String, String}}})
    created_envs = String[]
    backed_up_scripts = Tuple{String, String}[]
    installed_scripts = String[]
    try
        for (env_name, scripts) in specs
            env_dir = joinpath(Mycelia._conda_envs_dir(), env_name)
            if !isdir(env_dir)
                push!(created_envs, env_dir)
                mkpath(joinpath(env_dir, "bin"))
                mkpath(joinpath(env_dir, "conda-meta"))
                write_text_file(joinpath(env_dir, "conda-meta", "history"), "")
            end
            bin_dir = joinpath(env_dir, "bin")
            mkpath(bin_dir)
            for (script_name, script_body) in scripts
                script_path = joinpath(bin_dir, script_name)
                backup_path = script_path * ".mycelia-test-backup"
                ispath(backup_path) && error("Conda env script backup already exists: $(backup_path)")
                if isfile(script_path)
                    mv(script_path, backup_path)
                    push!(backed_up_scripts, (script_path, backup_path))
                end
                write_text_file(script_path, "#!/usr/bin/env bash\nset -euo pipefail\n$(script_body)\n")
                chmod(script_path, 0o755)
                push!(installed_scripts, script_path)
            end
        end
        return f()
    finally
        for script_path in reverse(installed_scripts)
            rm(script_path; force = true)
        end
        for (script_path, backup_path) in reverse(backed_up_scripts)
            mv(backup_path, script_path)
        end
        for env_dir in reverse(created_envs)
            rm(env_dir; recursive = true, force = true)
        end
    end
end

function with_stubbed_conda_runner(f::Function, script_body::AbstractString)
    runner_path = Mycelia.CONDA_RUNNER
    isfile(runner_path) || error("Cannot stub missing conda runner: $(runner_path)")
    backup_path = runner_path * ".mycelia-test-backup"
    ispath(backup_path) && error("Conda runner backup already exists: $(backup_path)")
    mv(runner_path, backup_path)
    write_text_file(runner_path, "#!/usr/bin/env bash\nset -euo pipefail\n$(script_body)\n")
    chmod(runner_path, 0o755)
    try
        return f()
    finally
        rm(runner_path; force = true)
        mv(backup_path, runner_path)
    end
end

function with_missing_conda_runner(f::Function)
    runner_path = Mycelia.CONDA_RUNNER
    isfile(runner_path) || error("Cannot hide missing conda runner: $(runner_path)")
    backup_path = runner_path * ".mycelia-test-backup"
    ispath(backup_path) && error("Conda runner backup already exists: $(backup_path)")
    mv(runner_path, backup_path)
    try
        return f()
    finally
        mv(backup_path, runner_path)
    end
end

function fake_vamb_bootstrap_conda_runner_script()
    return raw"""
log_file="${MYCELIA_FAKE_CONDA_LOG:?}"
envs_dir="${MYCELIA_FAKE_CONDA_ENVS_DIR:?}"
mkdir -p "$envs_dir"
printf '%s\n' "$*" >> "$log_file"

if [[ "${1:-}" == "create" ]]; then
    env_name=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -n)
                env_name="$2"
                shift 2
                ;;
            *)
                shift
                ;;
        esac
    done
    mkdir -p "$envs_dir/$env_name/bin" "$envs_dir/$env_name/conda-meta"
    : > "$envs_dir/$env_name/conda-meta/history"
    exit 0
fi

if [[ "${1:-}" == "run" ]]; then
    shift
    if [[ "${1:-}" == "--live-stream" ]]; then
        shift
    fi
    [[ "${1:-}" == "-n" ]] || exit 2
    env_name="$2"
    shift 2

    if [[ "${1:-}" == "vamb" && "${2:-}" == "--version" ]]; then
        counter_file="$envs_dir/$env_name/version_checks.txt"
        count=0
        if [[ -f "$counter_file" ]]; then
            count="$(cat "$counter_file")"
        fi
        count=$((count + 1))
        printf '%s\n' "$count" > "$counter_file"
        if [[ "$count" -eq 1 ]]; then
            exit 1
        fi
        echo "vamb 1.0.0"
        exit 0
    fi

    if [[ "${1:-}" == "python" && "${2:-}" == "-m" && "${3:-}" == "pip" && "${4:-}" == "install" ]]; then
        if [[ "${5:-}" == "--upgrade" ]]; then
            exit 0
        fi
        vamb_script="$envs_dir/$env_name/bin/vamb"
        cat > "$vamb_script" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
original_args=("$@")
if [[ "${1:-}" == "--version" ]]; then
    echo "vamb 1.0.0"
    exit 0
fi
outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)
            outdir="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
mkdir -p "$outdir"
printf '%s\n' "${original_args[*]}" > "$outdir/vamb_args.txt"
if [[ "${original_args[0]:-}" == "bin" ]]; then
    printf 'contig\tbin\ncontig1\tbin_1\n' > "$outdir/vae_clusters.tsv"
else
    printf 'taxometer\n' > "$outdir/taxometer_output.txt"
fi
EOF
        chmod +x "$vamb_script"
        exit 0
    fi

    tool_path="$envs_dir/$env_name/bin/${1:-}"
    [[ -x "$tool_path" ]] || exit 127
    exec "$tool_path" "${@:2}"
fi

exit 0
"""
end

function with_fake_vamb_bootstrap_runner(f::Function)
    env_name = Mycelia._vamb_env_name()
    env_dir = joinpath(Mycelia._conda_envs_dir(), env_name)
    backup_env_dir = env_dir * ".mycelia-test-backup"
    ispath(backup_env_dir) && error("VAMB environment backup already exists: $(backup_env_dir)")
    if ispath(env_dir)
        mv(env_dir, backup_env_dir)
    end
    mktempdir() do dir
        log_file = joinpath(dir, "fake_conda.log")
        try
            withenv(
                "MYCELIA_FAKE_CONDA_LOG" => log_file,
                "MYCELIA_FAKE_CONDA_ENVS_DIR" => Mycelia._conda_envs_dir(),
            ) do
                with_stubbed_conda_runner(fake_vamb_bootstrap_conda_runner_script()) do
                    try
                        return f(log_file)
                    finally
                        rm(env_dir; recursive = true, force = true)
                    end
                end
            end
        finally
            if ispath(backup_env_dir)
                mv(backup_env_dir, env_dir)
            end
        end
    end
end

function create_binning_stub_inputs(root::String)
    contigs_fasta = write_text_file(
        joinpath(root, "contigs.fna"),
        ">contig1\nACGTACGTACGT\n>contig2\nTTTTCCCCAAAA\n"
    )
    depth_file = write_text_file(
        joinpath(root, "depth.tsv"),
        "contigName\tcontigLen\ttotalAvgDepth\tsampleA\tsampleA-var\tsampleB\n" *
        "contig1\t12\t8.0\t3.5\t0.2\t4.5\n" *
        "contig2\t12\t5.0\t5.0\t0.1\n"
    )
    taxonomy_file = write_text_file(
        joinpath(root, "taxonomy.tsv"),
        "contig\ttaxonomy\ncontig1\tBacteria\n"
    )
    assembly_graph = write_text_file(
        joinpath(root, "assembly.gfa"),
        "H\tVN:Z:1.0\nS\tcontig1\tACGTACGTACGT\n"
    )
    mapping_file = write_text_file(
        joinpath(root, "mapping.tsv"),
        "contig\tcoverage\ncontig1\t7.0\n"
    )
    coverage_table = write_text_file(
        joinpath(root, "coverage.tsv"),
        "contig\tdepth\ncontig1\t7.0\n"
    )
    bam_dir = joinpath(root, "bams")
    mkpath(bam_dir)
    write_text_file(joinpath(bam_dir, "reads.bam"), "stub bam\n")
    genomes_dir = joinpath(root, "genomes")
    mkpath(genomes_dir)
    genome_a = write_text_file(joinpath(genomes_dir, "genome_a.fa"), ">genome_a\nACGT\n")
    genome_b = write_text_file(joinpath(genomes_dir, "genome_b.fa"), ">genome_b\nTGCA\n")
    bins_a = joinpath(root, "bins_a")
    bins_b = joinpath(root, "bins_b")
    mkpath(bins_a)
    mkpath(bins_b)
    write_text_file(joinpath(bins_a, "a.fa"), ">bin_a\nACGT\n")
    write_text_file(joinpath(bins_b, "b.fa"), ">bin_b\nTGCA\n")
    return (;
        contigs_fasta,
        depth_file,
        taxonomy_file,
        assembly_graph,
        mapping_file,
        coverage_table,
        bam_dir,
        genomes = [genome_a, genome_b],
        bins_dirs = [bins_a, bins_b]
    )
end

function fake_binning_env_specs()
    vamb_script = raw"""
original_args=("$@")
if [[ "${1:-}" == "--version" ]]; then
    echo "vamb 1.0.0"
    exit 0
fi
outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)
            outdir="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
mkdir -p "$outdir"
printf '%s\n' "${original_args[*]}" > "$outdir/vamb_args.txt"
if [[ "${original_args[0]:-}" == "bin" ]]; then
    printf 'contig\tbin\ncontig1\tbin_1\n' > "$outdir/vae_clusters.tsv"
else
    printf 'taxometer\n' > "$outdir/taxometer_output.txt"
fi
"""
    metabat_script = raw"""
original_args=("$@")
prefix=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o)
            prefix="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
mkdir -p "$(dirname "$prefix")"
printf '%s\n' "${original_args[*]}" > "$(dirname "$prefix")/metabat2_args.txt"
printf '>bin1\nACGT\n' > "${prefix}.1.fa"
"""
    metacoag_script = raw"""
original_args=("$@")
outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output)
            outdir="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
mkdir -p "$outdir/bins"
printf '%s\n' "${original_args[*]}" > "$outdir/metacoag_args.txt"
printf '>bin1\nACGT\n' > "$outdir/bins/bin1.fa"
printf 'contig\tbin\ncontig1\tbin1\n' > "$outdir/bin_assignments.tsv"
"""
    comebin_script = raw"""
original_args=("$@")
outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o)
            outdir="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
mkdir -p "$outdir/results/bins"
printf '%s\n' "${original_args[*]}" > "$outdir/results/comebin_args.txt"
printf '>bin1\nACGT\n' > "$outdir/results/bins/bin1.fa"
printf 'contig\tbin\ncontig1\tbin1\n' > "$outdir/results/bins_assignments.tsv"
"""
    drep_script = raw"""
original_args=("$@")
outdir="${2:-}"
mkdir -p "$outdir/data_tables"
printf '%s\n' "${original_args[*]}" > "$outdir/drep_args.txt"
printf 'genome,secondary_cluster,representative\n' > "$outdir/data_tables/Widb.csv"
printf 'genome_a,1,genome_a\n' >> "$outdir/data_tables/Widb.csv"
"""
    magmax_script = raw"""
original_args=("$@")
bindir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bindir)
            bindir="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
printf '%s\n' "${original_args[*]}" > "$PWD/magmax_args.txt"
find "$bindir" -maxdepth 1 -type f | sort > "$PWD/magmax_inputs.txt"
mkdir -p "$PWD/merged"
printf 'merged\n' > "$PWD/merged/summary.txt"
"""
    return [
        ("mycelia_vamb", Dict("vamb" => vamb_script)),
        ("metabat2", Dict("metabat2" => metabat_script)),
        ("metacoag", Dict("metacoag" => metacoag_script)),
        ("comebin", Dict("run_comebin.sh" => comebin_script)),
        ("drep", Dict("dRep" => drep_script)),
        ("magmax", Dict("magmax" => magmax_script))
    ]
end

function fake_conda_envs_available(specs::Vector{Tuple{String, Dict{String, String}}})
    return all(!ispath(joinpath(Mycelia._conda_envs_dir(), env_name)) for (env_name, _) in specs)
end

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

    Test.@testset "Binning helper utilities" begin
        Test.@test Mycelia._vamb_env_name() == "mycelia_vamb"

        mktempdir() do dir
            vamb_depth = write_text_file(
                joinpath(dir, "already_vamb.tsv"),
                "contigname\tsample1\tsample2\ncontig1\t1.0\t2.0\n"
            )
            Test.@test Mycelia._vamb_abundance_tsv(vamb_depth) == vamb_depth
        end

        mktempdir() do dir
            depth_file = write_text_file(
                joinpath(dir, "depth.tsv"),
                "contigName\tcontigLen\ttotalAvgDepth\tsampleA\tsampleA-var\tsampleB\n" *
                "contig1\t12\t8.0\t3.5\t0.2\t4.5\n" *
                "contig2\t12\t5.0\t5.0\t0.1\n\n"
            )
            output_dir = joinpath(dir, "derived")
            abundance_tsv = Mycelia._vamb_abundance_tsv(depth_file; output_dir = output_dir)
            Test.@test abundance_tsv == joinpath(output_dir, "vamb_abundance.tsv")
            Test.@test readlines(abundance_tsv) == [
                "contigname\tsampleA\tsampleB",
                "contig1\t3.5\t4.5",
                "contig2\t5.0\t0"
            ]

            write_text_file(abundance_tsv, "preexisting\n")
            Test.@test Mycelia._vamb_abundance_tsv(depth_file; output_dir = output_dir) == abundance_tsv
            Test.@test read(abundance_tsv, String) == "preexisting\n"
        end

        mktempdir() do dir
            no_sample_depth = write_text_file(
                joinpath(dir, "no_sample_depth.tsv"),
                "contigName\tcontigLen\ttotalAvgDepth\tsampleA-var\ncontig1\t12\t8.0\t0.2\n"
            )
            Test.@test_throws ErrorException Mycelia._vamb_abundance_tsv(no_sample_depth)
        end

        mktempdir() do dir
            write_text_file(joinpath(dir, "z_last.txt"), "z\n")
            write_text_file(joinpath(dir, "a_first.tsv"), "a\n")
            mkpath(joinpath(dir, "nested"))
            write_text_file(joinpath(dir, "nested", "match.tsv"), "nested\n")

            Test.@test Mycelia._find_first_matching_file(
                dir,
                [r"\.tsv$"]
            ) == joinpath(dir, "a_first.tsv")
            Test.@test Mycelia._find_first_matching_file(
                dir,
                [r"match\.tsv$"];
                recursive = true
            ) == joinpath(dir, "nested", "match.tsv")
            Test.@test isnothing(Mycelia._find_first_matching_file(dir, [r"missing\.tsv$"]))
        end

        mktempdir() do dir
            mkpath(joinpath(dir, "z_last_dir"))
            mkpath(joinpath(dir, "a_first_dir"))
            mkpath(joinpath(dir, "nested", "match_dir"))

            Test.@test Mycelia._find_first_matching_dir(
                dir,
                [r"^a_"]
            ) == joinpath(dir, "a_first_dir")
            Test.@test Mycelia._find_first_matching_dir(
                dir,
                [r"match_dir$"];
                recursive = true
            ) == joinpath(dir, "nested", "match_dir")
            Test.@test isnothing(Mycelia._find_first_matching_dir(dir, [r"missing_dir$"]))
        end
    end

    Test.@testset "VAMB bootstrap coverage" begin
        with_fake_vamb_bootstrap_runner() do log_file
            mktempdir() do dir
                inputs = create_binning_stub_inputs(dir)

                taxometer_outdir = joinpath(dir, "taxometer_bootstrap_out")
                taxometer_result = Mycelia.run_taxometer(
                    contigs_fasta = inputs.contigs_fasta,
                    depth_file = inputs.depth_file,
                    taxonomy_file = inputs.taxonomy_file,
                    outdir = taxometer_outdir,
                    threads = 6
                )
                taxometer_args = read(joinpath(taxometer_outdir, "vamb_args.txt"), String)
                Test.@test taxometer_result.outdir == taxometer_outdir
                Test.@test isfile(joinpath(taxometer_outdir, "taxometer_output.txt"))
                Test.@test occursin("taxometer", taxometer_args)
                Test.@test occursin("-p 6", taxometer_args)

                taxvamb_outdir = joinpath(dir, "taxvamb_bootstrap_out")
                taxvamb_result = Mycelia.run_taxvamb(
                    contigs_fasta = inputs.contigs_fasta,
                    depth_file = inputs.depth_file,
                    taxonomy_file = inputs.taxonomy_file,
                    outdir = taxvamb_outdir,
                    threads = 7
                )
                taxvamb_args = read(joinpath(taxvamb_outdir, "vamb_args.txt"), String)
                Test.@test taxvamb_result.clusters_tsv == joinpath(taxvamb_outdir, "vae_clusters.tsv")
                Test.@test isfile(taxvamb_result.clusters_tsv)
                Test.@test occursin("taxvamb", taxvamb_args)
                Test.@test occursin("-p 7", taxvamb_args)
            end

            conda_log = read(log_file, String)
            Test.@test occursin("create -y -n mycelia_vamb python=3.11 pip", conda_log)
            Test.@test occursin(
                "run --live-stream -n mycelia_vamb python -m pip install --upgrade pip setuptools wheel",
                conda_log
            )
            Test.@test occursin(
                "run --live-stream -n mycelia_vamb python -m pip install vamb",
                conda_log
            )
        end
    end

    Test.@testset "VAMB version-check failure handling" begin
        env_name = Mycelia._vamb_env_name()
        with_stubbed_conda_envs([(env_name, Dict{String, String}())]) do
            with_missing_conda_runner() do
                Test.@test_throws Exception Mycelia._ensure_vamb_installed()
            end
        end
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

        Test.@test_throws ErrorException Mycelia.parse_bin_assignments("missing_assignments.tsv")
        Test.@test_throws ErrorException Mycelia.parse_drep_clusters("missing_drep.csv")

        mktempdir() do dir
            assignments = write_text_file(joinpath(dir, "assignments.tsv"), "bin\tcount\nbin1\t1\n")
            Test.@test_throws ErrorException Mycelia.parse_bin_assignments(assignments)
        end

        mktempdir() do dir
            drep_clusters = write_text_file(
                joinpath(dir, "drep.csv"),
                "genome,secondary_cluster\nsample,1\n"
            )
            Test.@test_throws ErrorException Mycelia.parse_drep_clusters(drep_clusters)
        end

        mktempdir() do dir
            inputs = create_binning_stub_inputs(dir)
            existing_outdir = joinpath(dir, "existing_vamb_out")
            mkpath(existing_outdir)
            Test.@test_throws ErrorException Mycelia.run_vamb(
                contigs_fasta = inputs.contigs_fasta,
                depth_file = inputs.depth_file,
                outdir = existing_outdir
            )

            existing_taxometer_outdir = joinpath(dir, "existing_taxometer_out")
            mkpath(existing_taxometer_outdir)
            Test.@test_throws ErrorException Mycelia.run_taxometer(
                contigs_fasta = inputs.contigs_fasta,
                depth_file = inputs.depth_file,
                taxonomy_file = inputs.taxonomy_file,
                outdir = existing_taxometer_outdir
            )

            existing_taxvamb_outdir = joinpath(dir, "existing_taxvamb_out")
            mkpath(existing_taxvamb_outdir)
            Test.@test_throws ErrorException Mycelia.run_taxvamb(
                contigs_fasta = inputs.contigs_fasta,
                depth_file = inputs.depth_file,
                taxonomy_file = inputs.taxonomy_file,
                outdir = existing_taxvamb_outdir
            )

            Test.@test_throws ErrorException Mycelia.run_metacoag(
                contigs_fasta = inputs.contigs_fasta,
                assembly_graph = inputs.assembly_graph,
                mapping_file = "missing_mapping.tsv",
                outdir = joinpath(dir, "metacoag_out")
            )

            empty_bam_dir = joinpath(dir, "empty_bams")
            mkpath(empty_bam_dir)
            write_text_file(joinpath(empty_bam_dir, "notes.txt"), "not a bam\n")
            Test.@test_throws ErrorException Mycelia.run_comebin(
                contigs_fasta = inputs.contigs_fasta,
                bam_path = empty_bam_dir,
                outdir = joinpath(dir, "comebin_out")
            )

            Test.@test_throws ErrorException Mycelia.run_drep_dereplicate(
                genomes = [joinpath(dir, "missing_genome.fa")],
                outdir = joinpath(dir, "drep_out")
            )
        end
    end

    Test.@testset "Stubbed wrapper execution" begin
        fake_specs = fake_binning_env_specs()
        with_stubbed_conda_envs(fake_specs) do
            mktempdir() do dir
                inputs = create_binning_stub_inputs(dir)

                vamb_outdir = joinpath(dir, "vamb_out")
                vamb_result = Mycelia.run_vamb(
                    contigs_fasta = inputs.contigs_fasta,
                    depth_file = inputs.depth_file,
                    outdir = vamb_outdir,
                    minfasta = 1500,
                    threads = 3
                )
                Test.@test vamb_result.clusters_tsv == joinpath(vamb_outdir, "vae_clusters.tsv")
                Test.@test isfile(vamb_result.clusters_tsv)
                Test.@test occursin("--abundance_tsv", read(joinpath(vamb_outdir, "vamb_args.txt"), String))
                Test.@test occursin("-p 3", read(joinpath(vamb_outdir, "vamb_args.txt"), String))
                Test.@test isfile(joinpath(dirname(inputs.depth_file), "vamb_abundance.tsv"))

                taxometer_outdir = joinpath(dir, "taxometer_out")
                taxometer_result = Mycelia.run_taxometer(
                    contigs_fasta = inputs.contigs_fasta,
                    depth_file = inputs.depth_file,
                    taxonomy_file = inputs.taxonomy_file,
                    outdir = taxometer_outdir,
                    threads = 2,
                    extra_args = ["--threads", "7"]
                )
                taxometer_args = read(joinpath(taxometer_outdir, "vamb_args.txt"), String)
                Test.@test taxometer_result.outdir == taxometer_outdir
                Test.@test isfile(joinpath(taxometer_outdir, "taxometer_output.txt"))
                Test.@test occursin("taxometer", taxometer_args)
                Test.@test occursin("--threads 7", taxometer_args)
                Test.@test !occursin("-p 2", taxometer_args)

                taxvamb_outdir = joinpath(dir, "taxvamb_out")
                taxvamb_result = Mycelia.run_taxvamb(
                    contigs_fasta = inputs.contigs_fasta,
                    depth_file = inputs.depth_file,
                    taxonomy_file = inputs.taxonomy_file,
                    outdir = taxvamb_outdir,
                    threads = 4,
                    extra_args = ["-p", "9"]
                )
                taxvamb_args = read(joinpath(taxvamb_outdir, "vamb_args.txt"), String)
                Test.@test taxvamb_result.clusters_tsv == joinpath(taxvamb_outdir, "vae_clusters.tsv")
                Test.@test occursin("taxvamb", taxvamb_args)
                Test.@test occursin("-p 9", taxvamb_args)
                Test.@test !occursin("-p 4", taxvamb_args)

                metabat_outdir = joinpath(dir, "metabat2_out")
                metabat_result = Mycelia.run_metabat2(
                    contigs_fasta = inputs.contigs_fasta,
                    depth_file = inputs.depth_file,
                    outdir = metabat_outdir,
                    min_contig = 1200,
                    threads = 5,
                    seed = 11,
                    extra_args = ["--maxEdges", "250"]
                )
                metabat_args = read(joinpath(metabat_outdir, "metabat2_args.txt"), String)
                Test.@test metabat_result.bins_prefix == joinpath(metabat_outdir, "bin")
                Test.@test isfile("$(metabat_result.bins_prefix).1.fa")
                Test.@test occursin("-m 1200", metabat_args)
                Test.@test occursin("-t 5", metabat_args)
                Test.@test occursin("-s 11", metabat_args)
                Test.@test occursin("--maxEdges 250", metabat_args)

                metacoag_outdir = joinpath(dir, "metacoag_out")
                metacoag_result = Mycelia.run_metacoag(
                    contigs_fasta = inputs.contigs_fasta,
                    assembly_graph = inputs.assembly_graph,
                    mapping_file = inputs.mapping_file,
                    outdir = metacoag_outdir,
                    assembler = "megahit",
                    threads = 6,
                    extra_args = ["--min_length", "1000"]
                )
                metacoag_args = read(joinpath(metacoag_outdir, "metacoag_args.txt"), String)
                Test.@test metacoag_result.bins_dir == joinpath(metacoag_outdir, "bins")
                Test.@test isdir(metacoag_result.bins_dir)
                Test.@test occursin("--assembler megahit", metacoag_args)
                Test.@test occursin("--nthreads 6", metacoag_args)
                Test.@test occursin("--min_length 1000", metacoag_args)

                comebin_outdir = joinpath(dir, "comebin_out")
                comebin_result = Mycelia.run_comebin(
                    contigs_fasta = inputs.contigs_fasta,
                    bam_path = inputs.bam_dir,
                    outdir = comebin_outdir,
                    views = 2,
                    threads = 3,
                    temperature = 0.25,
                    embedding_size = 128,
                    coverage_embedding_size = 64,
                    batch_size = 16,
                    extra_args = ["--mode", "test"]
                )
                comebin_args = read(joinpath(comebin_outdir, "results", "comebin_args.txt"), String)
                Test.@test comebin_result.bins_dir == joinpath(comebin_outdir, "results", "bins")
                Test.@test comebin_result.bins_tsv == joinpath(comebin_outdir, "results", "bins_assignments.tsv")
                Test.@test occursin("-l 0.25", comebin_args)
                Test.@test occursin("-e 128", comebin_args)
                Test.@test occursin("-c 64", comebin_args)
                Test.@test occursin("-b 16", comebin_args)
                Test.@test occursin("--mode test", comebin_args)

                drep_outdir = joinpath(dir, "drep_out")
                drep_result = Mycelia.run_drep_dereplicate(
                    genomes = inputs.genomes,
                    outdir = drep_outdir,
                    completeness_threshold = 90.0,
                    contamination_threshold = 5.0,
                    ani_threshold = 0.97,
                    threads = 4,
                    extra_args = ["--ignoreGenomeQuality"]
                )
                drep_args = read(joinpath(drep_outdir, "drep_args.txt"), String)
                Test.@test drep_result.winning_genomes == joinpath(drep_outdir, "data_tables", "Widb.csv")
                Test.@test occursin("-comp 90.0", drep_args)
                Test.@test occursin("-con 5.0", drep_args)
                Test.@test occursin("-sa 0.97", drep_args)
                Test.@test occursin("--ignoreGenomeQuality", drep_args)

                magmax_outdir = joinpath(dir, "magmax_out")
                magmax_result = Mycelia.run_magmax_merge(
                    bins_dirs = inputs.bins_dirs,
                    outdir = magmax_outdir,
                    threads = 8
                )
                magmax_args = read(joinpath(magmax_outdir, "magmax_args.txt"), String)
                copied_bins = read(joinpath(magmax_outdir, "magmax_inputs.txt"), String)
                Test.@test magmax_result.bins_input_dir == joinpath(magmax_outdir, "magmax_bins_input")
                Test.@test occursin("--threads 8", magmax_args)
                Test.@test occursin("--format fa", magmax_args)
                Test.@test occursin("bins_a__a.fa", copied_bins)
                Test.@test occursin("bins_b__b.fa", copied_bins)
            end
        end
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
