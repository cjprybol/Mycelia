# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/sbatch_depot_before_version_guard_test.jl")'
# ```
#
# td-67fj: the free-standing run_*.sbatch wrappers must export JULIA_DEPOT_PATH
# and LD_LIBRARY_PATH before the live `julia --version` guard, and must keep
# their existing version pin (`+lts` on Lawrencium, `julia/1.10.10` on NERSC).
# The two wrappers that create a node-local TMPDIR also pin a free-space floor
# and a tmpfs opt-in. This file only reads the scripts. It does not run sbatch.
#
# Not in the ordered set:
# - benchmarking/run_ont_read_identity_lrc.sbatch belongs to PR #509.
# - benchmarking/run_d5_grid_*.sbatch is the Lawrencium pilot; its version
#   guard lives in d5_grid_task.sh.

import Test

const BENCH = normpath(joinpath(@__DIR__, "..", "..", "benchmarking"))

# Lawrencium `+lts` launchers, including the td-2rxh evidence wrapper.
const LRC_PLUS_LTS = (
    "run_ecoli_qualmer_memory_gate_lrc.sbatch",
    "run_rgv_sweep_lrc.sbatch",
    "run_talk_subset_benchmark_lrc.sbatch",
    "run_td2rxh_evidence_lrc.sbatch",
    "run_track_a_baseline_lrc.sbatch",
)

# NERSC `module load julia/1.10.10`, including the non-LRC talk wrapper.
const NERSC_EXACT = (
    "run_rgv_sweep_nersc.sbatch",
    "run_talk_subset_benchmark.sbatch",
    "run_track_a_baseline_nersc.sbatch",
)

const ORDERED = (LRC_PLUS_LTS..., NERSC_EXACT...)

# PR #509 owns this file. It has the same depot/version pattern and is listed
# so a new wrapper cannot appear unnoticed, but this test does not require the
# export-before-guard order there.
const OWNED_BY_PR_509 = "run_ont_read_identity_lrc.sbatch"

const TMPDIR_GUARDS = (
    ("run_ecoli_qualmer_memory_gate_lrc.sbatch",
        "QUALMER_GATE_ALLOW_TMPFS", "QUALMER_GATE_TMP_MIN_FREE_KB"),
    ("run_td2rxh_evidence_lrc.sbatch",
        "TD2RXH_ALLOW_TMPFS", "TD2RXH_TMP_MIN_FREE_KB"),
)

# True when the line is shell code rather than a comment or a blank.
function is_shell_code(line::AbstractString)
    stripped = lstrip(line)
    return !isempty(stripped) && !startswith(stripped, '#')
end

# 1-based line number of the first code line for which `pred` is true.
function first_code_line(lines, pred)::Union{Int, Nothing}
    for (i, line) in enumerate(lines)
        is_shell_code(line) || continue
        pred(line) && return i
    end
    return nothing
end

function code_lines_matching(lines, pred)
    return [i for (i, line) in enumerate(lines) if is_shell_code(line) && pred(line)]
end

function read_sbatch(name::AbstractString)
    path = joinpath(BENCH, name)
    return split(read(path, String), '\n'; keepempty = true)
end

Test.@testset "sbatch depot exports precede the version guard (td-67fj)" begin
    Test.@testset "the free-standing wrappers are exactly the known set" begin
        found = String[]
        for name in readdir(BENCH)
            startswith(name, "run_") && endswith(name, ".sbatch") || continue
            lines = read_sbatch(name)
            has_version = first_code_line(lines, line -> occursin("--version", line)) !== nothing
            has_depot = first_code_line(lines,
                line -> occursin(r"^export JULIA_DEPOT_PATH=", lstrip(line))) !== nothing
            if has_version && has_depot
                push!(found, name)
            end
        end
        Test.@test sort(found) == sort([ORDERED..., OWNED_BY_PR_509])
    end

    Test.@testset "exports come before the live julia --version guard" begin
        for name in ORDERED
            lines = read_sbatch(name)
            version_at = first_code_line(lines, line -> occursin("--version", line))
            depot_lines = code_lines_matching(lines,
                line -> occursin(r"^export JULIA_DEPOT_PATH=", lstrip(line)))
            ld_lines = code_lines_matching(lines,
                line -> occursin(r"^export LD_LIBRARY_PATH=", lstrip(line)))
            Test.@test version_at !== nothing
            Test.@test length(depot_lines) == 1
            Test.@test length(ld_lines) == 1
            if version_at !== nothing && length(depot_lines) == 1 && length(ld_lines) == 1
                Test.@test depot_lines[1] < version_at
                Test.@test ld_lines[1] < version_at
                Test.@test occursin("LD_LIBRARY_PATH=\"\"", lines[ld_lines[1]])
            end
        end
    end

    Test.@testset "the version pin is unchanged" begin
        for name in LRC_PLUS_LTS
            text = read(joinpath(BENCH, name), String)
            Test.@test occursin("JULIA_ARGS=(+lts)", text)
            Test.@test !occursin("JULIA_ARGS=(+1.", text)
        end
        for name in NERSC_EXACT
            text = read(joinpath(BENCH, name), String)
            Test.@test occursin("module load julia/1.10.10", text)
            Test.@test !occursin("JULIA_ARGS=(+", text)
        end
    end

    Test.@testset "TMPDIR free-space floor and tmpfs opt-in" begin
        for (name, allow_var, floor_var) in TMPDIR_GUARDS
            lines = read_sbatch(name)
            text = join(lines, '\n')
            allow_arm = filter(is_shell_code, lines)
            silent = filter(line -> occursin("ext2/ext3|", line), allow_arm)
            Test.@test length(silent) == 1
            Test.@test !occursin("tmpfs", only(silent))
            Test.@test any(line -> lstrip(line) == "tmpfs)", allow_arm)
            Test.@test occursin("$(allow_var)=1", text)
            Test.@test occursin("WARN:", text)
            Test.@test occursin("counts against --mem", text)
            Test.@test occursin("$(floor_var):-2097152", text)
            Test.@test occursin("below the", text)
            Test.@test occursin("df -P -k", text)
            Test.@test occursin("non-negative integer", text)
        end
    end

    Test.@testset "resolve-mode hashing note on the qualmer gate" begin
        text = read(joinpath(BENCH, "run_ecoli_qualmer_memory_gate_lrc.sbatch"), String)
        Test.@test occursin("Resolve mode leaves manifest_src_sha256 empty", text)
        Test.@test occursin(
            "manifest_src_sha256=\${QUALMER_GATE_MANIFEST_SRC_SHA:-<resolve mode>}",
            text)
        # The note records the empty source hash; it does not start hashing one.
        Test.@test occursin("there is no source Manifest to", text)
    end
end
