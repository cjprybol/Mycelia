# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/traversal_min_quality_guard_test.jl")'
# ```
#
# td-zhua: in-range traversal_min_quality must not silently empty a graph;
# prune counts actual deletions; get_vertex_joint_quality validates
# observation lengths; the withdrawn 60.0 band rationale is not stated
# as fact.

import Test
import Mycelia
import FASTX
import BioSequences
import Kmers
import MetaGraphsNext
import Logging

const R = Mycelia.Rhizomorph
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const HIGH_REF = "ATCGGCTAATGCCGATTGCACGTACGTTAGCTAGGCATG"  # 39 bp, no k=7 repeats
const LOW_REF  = "GGTTAACCGGAATTCCGGTTAACCATGCGTTAGCCTAAGGTTAACCGGAATTCCGGTTAACC"

# One observation per k-mer so joint quality equals the raw per-base Phred
# (overlapping tiles would compound and clamp, hiding the floor).
function singleton_fastq(ref::AbstractString, qual_char::Char; k::Int = 7, prefix::String = "r")
    reads = FASTX.FASTQ.Record[]
    for i in 1:(length(ref) - k + 1)
        sub = ref[i:(i + k - 1)]
        push!(reads, FASTX.FASTQ.Record("$(prefix)$(i)", sub, repeat(qual_char, k)))
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

Test.@testset "traversal_min_quality guards (td-zhua)" begin
    Test.@testset "prune counts actual deletions" begin
        reads = singleton_fastq(HIGH_REF, 'I')  # Q40
        graph = R.build_qualmer_graph(R._prepare_fastq_observations(reads), 7;
            mode = :singlestrand)
        before = length(collect(MetaGraphsNext.labels(graph)))
        Test.@test before > 0
        pruned, n = R._prune_qualmer_graph_by_quality(graph, 50.0)  # all joint Q40 < 50
        after = length(collect(MetaGraphsNext.labels(pruned)))
        Test.@test n == before - after
        Test.@test after == 0
        Test.@test n == before
    end

    Test.@testset "in-range floor that empties a non-empty graph errors" begin
        reads = singleton_fastq(HIGH_REF, 'I')
        cfg = R.AssemblyConfig(;
            k = 7,
            use_quality_scores = true,
            graph_mode = R.SingleStrand,
            traversal_weighting = :quality,
            traversal_min_quality = 200.0,
        )
        result, warnings, err = run_capturing_warnings() do
            R.assemble_genome(reads, cfg)
        end
        Test.@test err !== nothing
        Test.@test result === nothing
        msg = sprint(showerror, err)
        Test.@test occursin("traversal_min_quality", msg)
        Test.@test occursin("non-empty", msg)
    end

    Test.@testset "pruning a large fraction warns and records vertices_before" begin
        # Q40 k-mers (joint 40) plus a longer stretch of Q2 k-mers (joint 2) so
        # a floor of 20 prunes a majority but not all.
        reads = vcat(singleton_fastq(HIGH_REF, 'I'; prefix = "h"),
                     singleton_fastq(LOW_REF, '#'; prefix = "l"))  # '#' is Phred+33 Q2
        cfg = R.AssemblyConfig(;
            k = 7,
            use_quality_scores = true,
            graph_mode = R.SingleStrand,
            traversal_weighting = :quality,
            traversal_min_quality = 20.0,
        )
        result, warnings, err = run_capturing_warnings() do
            R.assemble_genome(reads, cfg)
        end
        Test.@test err === nothing
        Test.@test result !== nothing
        # Singleton k-mer reads have no transitions, so contigs may be empty;
        # the guard is that a partial prune must not error and must warn.
        Test.@test result !== nothing
        Test.@test any(occursin("pruned", w) for w in warnings)
        stats = result === nothing ? Dict{String,Any}() : result.assembly_stats
        Test.@test haskey(stats, "traversal_quality_vertices_before")
        if haskey(stats, "traversal_quality_vertices_before") &&
           haskey(stats, "traversal_quality_pruned_vertices")
            Test.@test stats["traversal_quality_vertices_before"] >
                       stats["traversal_quality_pruned_vertices"]
            Test.@test stats["traversal_quality_pruned_vertices"] > 0
            frac = stats["traversal_quality_pruned_vertices"] /
                   stats["traversal_quality_vertices_before"]
            Test.@test frac >= 0.5
        end
    end

    Test.@testset "get_vertex_joint_quality rejects unequal observation lengths" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = R.QualmerVertexData(kmer)
        long = R.QualityEvidenceEntry(1, R.Forward, UInt8[63, 63, 63])  # Q30 + 33
        short = R.QualityEvidenceEntry(1, R.Forward, UInt8[63, 63])     # length 2
        vertex.evidence["test"] = Dict(
            "read_a" => Set([long]),
            "read_b" => Set([short]),
        )
        err = Test.@test_throws ArgumentError R.get_vertex_joint_quality(vertex, "test")
        Test.@test occursin("length", sprint(showerror, err.value))
    end

    Test.@testset "60.0 default's withdrawn band rationale is not stated as fact" begin
        assembly_src = read(joinpath(REPO_ROOT, "src", "rhizomorph", "assembly.jl"), String)
        # Isolate the traversal_min_quality field comment (before the constructor).
        m = match(r"traversal_min_quality::Float64\n"s, assembly_src)
        Test.@test m !== nothing
        # The comment block immediately above the field.
        idx = findfirst("JOINT-quality floor consulted ONLY when traversal_weighting=:quality", assembly_src)
        Test.@test idx !== nothing
        block = assembly_src[idx.start:m.offset]
        Test.@test occursin(r"(?i)withdrawn", block)
        Test.@test !occursin(
            r"100-255\s+extreme\s+\(many high-quality observations; essentially certain\)",
            block)
        # The probe README justified the same default with the withdrawn bands.
        readme = read(joinpath(REPO_ROOT, "benchmarking", "results",
            "qualmer_quality_channel_probe", "README.md"), String)
        Test.@test !occursin(r"60-100 high \(multiple\), 100-255 extreme\. The 60\.0 default", readme)
    end
end
