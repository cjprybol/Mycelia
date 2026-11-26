# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/enum_flag_fail_loud_guard_test.jl")'
# ```
#
# td-70kz: fail-loud guard so a known-broken-but-exposed GraphMode or
# AssemblyConfig mode flag cannot ship silently. Every public GraphMode value
# and every closed mode-flag value must either have a functional assembly
# test or emit a @warn (or a louder error) when used.
#
# memory_profile, qualmer_memory_profile, corrector, strategy, layout, and
# traversal_weighting are read from the constants AssemblyConfig validates
# against. graph_mode, sequencing_tech, and olc_tool are discovered from
# GraphMode, _correction_profile_technologies, and _olc_taxonomy.

import Test
import Mycelia
import FASTX
import BioSequences
import Logging

const R = Mycelia.Rhizomorph
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# A small reference with no exact repeats at k=7 (same fixture idea as the
# efficiency-modes tests), so a working mode produces a non-empty assembly.
const GUARD_REF = "ATCGGCTAATGCCGATTGCACGTACGTTAGCTAGGCATG"

function guard_tiling_reads(ref::AbstractString; read_len::Int = 15)
    reads = FASTX.FASTA.Record[]
    for i in 1:(length(ref) - read_len + 1)
        push!(reads, FASTX.FASTA.Record("r$(i)", ref[i:(i + read_len - 1)]))
    end
    return reads
end

function guard_tiling_fastq(ref::AbstractString; read_len::Int = 15)
    reads = FASTX.FASTQ.Record[]
    for i in 1:(length(ref) - read_len + 1)
        sub = ref[i:(i + read_len - 1)]
        push!(reads, FASTX.FASTQ.Record("r$(i)", sub, repeat("I", length(sub))))
    end
    return reads
end

struct WarnCollector <: Logging.AbstractLogger
    msgs::Vector{String}
end
Logging.min_enabled_level(::WarnCollector) = Logging.Warn
Logging.shouldlog(::WarnCollector, args...) = true
Logging.catch_exceptions(::WarnCollector) = false
function Logging.handle_message(logger::WarnCollector, level, message, _module, group, id,
        file, line; kwargs...)
    push!(logger.msgs, string(message))
    return nothing
end

function run_capturing_warnings(f)
    collector = WarnCollector(String[])
    result = nothing
    err = nothing
    try
        Logging.with_logger(collector) do
            result = f()
        end
    catch e
        err = e
    end
    return result, collector.msgs, err
end

function is_fail_loud(result, warnings, err)
    ok = result !== nothing && !isempty(result.contigs)
    loud = !isempty(warnings) || err !== nothing
    return ok || loud
end

# Relative test paths that already exercise a mode. The guard also runs its
# own smokes; these entries keep the inventory honest when a new value appears
# in the constructor but is not added here.
const MODE_FLAG_COVERAGE = Dict{Tuple{Symbol, Any}, Tuple{String, String}}(
    # GraphMode (enum values, not Symbols)
    (:graph_mode, R.SingleStrand) =>
        ("test/4_assembly/dna_kmer_singlestrand_test.jl", "SingleStrand"),
    (:graph_mode, R.DoubleStrand) =>
        ("test/4_assembly/rhizomorph_efficiency_modes_test.jl", "graph_mode = Mycelia.Rhizomorph.DoubleStrand"),
    (:graph_mode, R.Canonical) =>
        ("test/4_assembly/rhizomorph_efficiency_modes_test.jl", "graph_mode = Mycelia.Rhizomorph.Canonical"),
    # memory_profile
    (:memory_profile, :full) =>
        ("test/4_assembly/rhizomorph_efficiency_modes_test.jl", "memory_profile = :full"),
    (:memory_profile, :lightweight) =>
        ("test/4_assembly/rhizomorph_efficiency_modes_test.jl", "memory_profile = :lightweight"),
    (:memory_profile, :ultralight) =>
        ("test/4_assembly/rhizomorph_efficiency_modes_test.jl", "memory_profile = :ultralight"),
    (:memory_profile, :lightweight_quality) =>
        ("test/4_assembly/rhizomorph_qualmer_aggregate_quality_test.jl", "memory_profile = :lightweight_quality"),
    (:memory_profile, :ultralight_quality) =>
        ("test/4_assembly/rhizomorph_qualmer_aggregate_quality_test.jl", ":ultralight_quality"),
    # qualmer_memory_profile — :lightweight_quality is exercised in this file
    (:qualmer_memory_profile, :full) =>
        ("test/4_assembly/rhizomorph_qualmer_aggregate_quality_test.jl", "qualmer_memory_profile = :full"),
    (:qualmer_memory_profile, :ultralight_quality) =>
        ("test/4_assembly/rhizomorph_qualmer_aggregate_quality_test.jl", "qualmer_memory_profile = :ultralight_quality"),
    (:qualmer_memory_profile, :lightweight_quality) =>
        ("test/4_assembly/enum_flag_fail_loud_guard_test.jl", "qualmer_memory_profile = :lightweight_quality"),
    # corrector / strategy / layout / traversal_weighting
    (:corrector, :none) =>
        ("test/4_assembly/variation_preservation_holdout_test.jl", "corrector = :none"),
    (:corrector, :iterative) =>
        ("test/4_assembly/rhizomorph_iterative_corrector_wiring_test.jl", "corrector = :iterative"),
    (:strategy, :scalable) =>
        ("test/4_assembly/scalable_corrector_strategy_test.jl", "strategy = :scalable"),
    (:strategy, :exhaustive) =>
        ("test/4_assembly/scalable_corrector_strategy_test.jl", "strategy = :exhaustive"),
    (:layout, :native) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "default_cfg.layout == :native"),
    (:layout, :olc) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "layout = :olc"),
    (:traversal_weighting, :evidence) =>
        ("test/4_assembly/qualmer_quality_channel_test.jl", "traversal_weighting"),
    (:traversal_weighting, :quality) =>
        ("test/4_assembly/qualmer_quality_channel_test.jl", "traversal_weighting = :quality"),
    # sequencing_tech
    (:sequencing_tech, :illumina) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "sequencing_tech = :illumina"),
    (:sequencing_tech, :ultima) =>
        ("test/4_assembly/indel_profile_assembly_test.jl", "sequencing_tech = :ultima"),
    (:sequencing_tech, :nanopore) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "sequencing_tech = :nanopore"),
    (:sequencing_tech, :pacbio) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "sequencing_tech = :pacbio"),
    (:sequencing_tech, :pacbio_clr) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "sequencing_tech = :pacbio_clr"),
    (:sequencing_tech, :pacbio_hifi) =>
        ("test/4_assembly/indel_profile_assembly_test.jl", "sequencing_tech = :pacbio_hifi"),
    # olc_tool
    (:olc_tool, :auto) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "olc_tool = :auto"),
    (:olc_tool, :megahit) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "olc_tool = :megahit"),
    (:olc_tool, :metaspades) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "olc_tool = :metaspades"),
    (:olc_tool, :flye) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "olc_tool = :flye"),
    (:olc_tool, :metaflye) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", ":metaflye"),
    (:olc_tool, :canu) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "olc_tool = :canu"),
    (:olc_tool, :hifiasm) =>
        ("test/4_assembly/hybrid_olc_config_test.jl", "olc_tool = :hifiasm"),
)

function discovered_mode_values()
    tax = R._olc_taxonomy()
    return Dict{Symbol, Vector}(
        :graph_mode => collect(instances(R.GraphMode)),
        :memory_profile => collect(R._VALID_MEMORY_PROFILES),
        :qualmer_memory_profile => collect(R._VALID_QUALMER_MEMORY_PROFILES),
        :corrector => collect(R._VALID_CORRECTORS),
        :strategy => collect(R._VALID_STRATEGIES),
        :layout => collect(R._VALID_LAYOUTS),
        :traversal_weighting => collect(R._VALID_TRAVERSAL_WEIGHTINGS),
        :sequencing_tech => collect(R._correction_profile_technologies()),
        :olc_tool => collect((:auto, tax.short_read_tools..., tax.long_read_tools...)),
    )
end

function repo_file_has(relpath::AbstractString, needle::AbstractString)
    path = joinpath(REPO_ROOT, relpath)
    isfile(path) || return false
    return occursin(needle, read(path, String))
end

Test.@testset "enum/flag fail-loud guard (td-70kz)" begin
    Test.@testset "every GraphMode assembles or is fail-loud" begin
        reads = guard_tiling_reads(GUARD_REF)
        for mode in instances(R.GraphMode)
            Test.@testset "GraphMode $(mode)" begin
                cfg = R.AssemblyConfig(; k = 7, graph_mode = mode, use_quality_scores = false)
                result, warnings, err = run_capturing_warnings() do
                    R.assemble_genome(reads, cfg)
                end
                Test.@test is_fail_loud(result, warnings, err)
                # Working modes must actually produce contigs, not just warn.
                # A brand-new unimplemented value may warn or error instead.
                if err === nothing && isempty(warnings)
                    Test.@test result !== nothing && !isempty(result.contigs)
                end
            end
        end
    end

    Test.@testset "qualmer_memory_profile = :lightweight_quality assembles" begin
        reads = guard_tiling_fastq(GUARD_REF)
        cfg = R.AssemblyConfig(;
            k = 7,
            use_quality_scores = true,
            qualmer_memory_profile = :lightweight_quality,
        )
        result, warnings, err = run_capturing_warnings() do
            R.assemble_genome(reads, cfg)
        end
        Test.@test is_fail_loud(result, warnings, err)
        if err === nothing && isempty(warnings)
            Test.@test result !== nothing && !isempty(result.contigs)
        end
    end

    Test.@testset "inventory covers every discovered mode flag" begin
        discovered = discovered_mode_values()
        for (flag, values) in discovered
            Test.@testset "$(flag)" begin
                for value in values
                    entry = get(MODE_FLAG_COVERAGE, (flag, value), nothing)
                    Test.@test entry !== nothing  # missing coverage for $(flag)=$(value)
                    if entry !== nothing
                        relpath, needle = entry
                        Test.@test repo_file_has(relpath, needle)  # $(flag)=$(value) in $(relpath)
                    end
                end
            end
        end
        # The table must not list a stale value the constructor no longer accepts.
        for (flag, value) in keys(MODE_FLAG_COVERAGE)
            Test.@test haskey(discovered, flag)
            Test.@test value in discovered[flag]
        end
    end
end
