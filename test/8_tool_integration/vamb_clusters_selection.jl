# Pure-Julia coverage for selecting VAMB's contig-to-bin clusters table.
#
# VAMB writes `<prefix>_clusters_metadata.tsv` (per-cluster statistics) next to
# `<prefix>_clusters_unsplit.tsv` and, when binsplitting, `<prefix>_clusters_split.tsv`
# (both `clustername\tcontigname`). The metadata table sorts first by name, so a
# name-sorted pattern match returned it instead of an assignments table.
#
# The wrapper testset reaches run_vamb/run_taxvamb through a fake conda runner
# selected with MYCELIA_CONDA_RUNNER, inside a temp dir; no real conda env or
# runner is touched.
#
# Runs in default CI (no conda, no external tools):
#   julia --project=. -e 'include("test/8_tool_integration/vamb_clusters_selection.jl")'

import Test
import Mycelia

function _touch_all(dir::String, names::Vector{String})
    for name in names
        touch(joinpath(dir, name))
    end
    return dir
end

function _write_executable(path::String, body::String)
    mkpath(dirname(path))
    write(path, body)
    chmod(path, 0o755)
    return path
end

# Fake conda runner at <root>/bin/conda; `_conda_envs_dir()` resolves to <root>/envs.
const _FAKE_VAMB_CONDA_RUNNER = raw"""
#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
case "${1:-}" in
    env)
        exit 0
        ;;
    run)
        shift
        [[ "${1:-}" == "--live-stream" ]] && shift
        [[ "${1:-}" == "-n" ]] || exit 2
        env_name="$2"
        shift 2
        exec "$root/envs/$env_name/bin/$1" "${@:2}"
        ;;
esac
echo "fake conda: unexpected command: $*" >&2
exit 3
"""

# Fake `vamb` writing the file names current VAMB writes for `bin default` and
# `bin taxvamb`, including the metadata table that sorts first.
const _FAKE_VAMB = raw"""
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "--version" ]]; then
    echo "vamb 5.0.4"
    exit 0
fi
mode="${2:-}"
outdir=""
args=("$@")
while [[ $# -gt 0 ]]; do
    if [[ "$1" == "--outdir" ]]; then
        outdir="$2"
        shift 2
    else
        shift
    fi
done
mkdir "$outdir"
printf '%s\n' "${args[*]}" > "$outdir/vamb_args.txt"
model="vae"
[[ "$mode" == "taxvamb" ]] && model="vaevae"
printf 'name\tradius\n' > "$outdir/${model}_clusters_metadata.tsv"
printf 'clustername\tcontigname\nbin_1\tcontig1\n' > "$outdir/${model}_clusters_unsplit.tsv"
printf 'clustername\tcontigname\nbin_1\tcontig1\n' > "$outdir/${model}_clusters_split.tsv"
"""

Test.@testset "VAMB clusters table selection" begin
    Test.@testset "prefers split over unsplit, never metadata" begin
        mktempdir() do dir
            _touch_all(dir,
                ["vae_clusters_metadata.tsv", "vae_clusters_split.tsv",
                    "vae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_split.tsv")
        end
    end

    Test.@testset "falls back to unsplit when binsplitting is off" begin
        mktempdir() do dir
            _touch_all(dir, ["vae_clusters_metadata.tsv", "vae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_unsplit.tsv")
        end
    end

    Test.@testset "metadata alone is not an assignments table" begin
        mktempdir() do dir
            _touch_all(dir, ["vae_clusters_metadata.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") === nothing
        end
    end

    Test.@testset "legacy single-file names" begin
        mktempdir() do dir
            _touch_all(dir, ["vae_clusters.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters.tsv")
        end
        mktempdir() do dir
            _touch_all(dir, ["clusters.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "clusters.tsv")
        end
    end

    Test.@testset "current names outrank the defensive fallback names" begin
        mktempdir() do dir
            _touch_all(dir, [
                "vae_clusters_unsplit.tsv", "vae_clusters.tsv", "clusters.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_unsplit.tsv")
        end
        mktempdir() do dir
            _touch_all(dir, ["vae_clusters.tsv", "clusters.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters.tsv")
        end
    end

    Test.@testset "prefix selects the model's own tables" begin
        mktempdir() do dir
            # TaxVAMB names contain "vae_clusters" as a substring; the prefix must
            # match exactly so a `vae` lookup does not pick up `vaevae` tables.
            _touch_all(dir,
                ["vaevae_clusters_metadata.tsv", "vaevae_clusters_split.tsv",
                    "vaevae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") === nothing
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vaevae") ==
                       joinpath(dir, "vaevae_clusters_split.tsv")
        end
        mktempdir() do dir
            # AVAMB writes both VAE and adversarial tables in one directory.
            _touch_all(dir, ["aae_z_clusters_split.tsv", "vae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_unsplit.tsv")
        end
    end

    Test.@testset "ignores directories and missing outdir" begin
        mktempdir() do dir
            mkdir(joinpath(dir, "vae_clusters_split.tsv"))
            _touch_all(dir, ["vae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_unsplit.tsv")
        end
        mktempdir() do dir
            Test.@test Mycelia._vamb_clusters_tsv(joinpath(dir, "absent"); prefix = "vae") ===
                       nothing
        end
    end

    Test.@testset "run_vamb and run_taxvamb return the split assignments" begin
        mktempdir() do root
            runner = _write_executable(joinpath(root, "bin", "conda"), _FAKE_VAMB_CONDA_RUNNER)
            env_dir = joinpath(root, "envs", Mycelia._vamb_env_name())
            mkpath(joinpath(env_dir, "conda-meta"))
            _write_executable(joinpath(env_dir, "bin", "vamb"), _FAKE_VAMB)

            inputs = joinpath(root, "inputs")
            mkpath(inputs)
            contigs = joinpath(inputs, "contigs.fna")
            write(contigs, ">contig1\nACGTACGTACGT\n")
            # VAMB-format abundance table, used as-is (no JGI conversion).
            depth = joinpath(inputs, "abundance.tsv")
            write(depth, "contigname\tS1\ncontig1\t2.0\n")
            taxonomy = joinpath(inputs, "taxonomy.tsv")
            write(taxonomy, "contigs\tpredictions\ncontig1\td__Bacteria\n")

            # _check_conda_env_exists exports CONDA_EXE for the active runner;
            # restore the caller's values so the fake path does not leak.
            withenv("MYCELIA_CONDA_RUNNER" => runner,
                "CONDA_EXE" => get(ENV, "CONDA_EXE", nothing),
                "CONDA_PYTHON_EXE" => get(ENV, "CONDA_PYTHON_EXE", nothing)) do
                vamb_outdir = joinpath(root, "out", "vamb")
                vamb_result = Mycelia.run_vamb(
                    contigs_fasta = contigs, depth_file = depth, outdir = vamb_outdir,
                    threads = 2)
                Test.@test vamb_result.clusters_tsv ==
                           joinpath(vamb_outdir, "vae_clusters_split.tsv")
                Test.@test occursin("bin default",
                    read(joinpath(vamb_outdir, "vamb_args.txt"), String))

                taxvamb_outdir = joinpath(root, "out", "taxvamb")
                taxvamb_result = Mycelia.run_taxvamb(
                    contigs_fasta = contigs, depth_file = depth,
                    taxonomy_file = taxonomy, outdir = taxvamb_outdir, threads = 2)
                Test.@test taxvamb_result.clusters_tsv ==
                           joinpath(taxvamb_outdir, "vaevae_clusters_split.tsv")
                Test.@test occursin("bin taxvamb",
                    read(joinpath(taxvamb_outdir, "vamb_args.txt"), String))
            end
        end
    end
end
