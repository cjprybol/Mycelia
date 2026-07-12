# td-yymj: hybrid-OLC route (a) CONFIG + routing tests (no external tools).
#
# The hybrid-OLC route hands Stage-1-corrected reads to an external OLC assembler
# (layout=:olc). This file covers the parts that need NO external assembler and so
# run in default CI: the AssemblyConfig validation surface (layout / olc_tool /
# olc_options), the :auto tool resolution, and the per-tool adapter's guard. The
# real end-to-end run (corrected reads -> megahit/metaspades -> contigs) is the
# gated companion test in third_party_assemblers_hybrid_olc_route.jl.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/hybrid_olc_config_test.jl")'

import Test
import Mycelia

const R = Mycelia.Rhizomorph

Test.@testset "hybrid-OLC route (a) config + routing (td-yymj)" begin
    Test.@testset "layout selector validation" begin
        # Default: native layout, auto tool.
        default_cfg = R.AssemblyConfig(; k = 13)
        Test.@test default_cfg.layout == :native
        Test.@test default_cfg.olc_tool == :auto
        Test.@test default_cfg.olc_options == (;)

        # Valid :olc config (short-read tech + iterative corrector).
        olc_cfg = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, layout = :olc, sequencing_tech = :illumina)
        Test.@test olc_cfg.layout == :olc
        Test.@test olc_cfg.olc_tool == :auto

        # Unknown layout is rejected.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13, layout = :bogus)

        # :olc requires the corrector — an OLC layout with no correction is just the
        # plain external assembler.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13, layout = :olc,
            corrector = :none)
    end

    Test.@testset "olc_tool <-> sequencing_tech compatibility" begin
        # Unknown tool rejected under :olc.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
            corrector = :iterative, strategy = :scalable, layout = :olc,
            olc_tool = :not_a_tool)

        # Explicit short-read tools are accepted with a short-read tech.
        for tool in (:megahit, :metaspades)
            cfg = R.AssemblyConfig(; k = 13, corrector = :iterative,
                strategy = :scalable, layout = :olc, olc_tool = tool,
                sequencing_tech = :illumina)
            Test.@test cfg.olc_tool == tool
        end

        # Short-read tool with a long-read tech is a mismatch.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
            corrector = :iterative, strategy = :scalable, layout = :olc,
            olc_tool = :megahit, sequencing_tech = :nanopore)

        # Long-read tools (td-wvto) are accepted with a compatible long-read tech.
        for (tool, tech, opts) in ((:flye, :nanopore, (;)), (:flye, :pacbio, (;)),
            (:metaflye, :nanopore, (;)), (:hifiasm, :pacbio, (;)),
            (:canu, :pacbio, (; genome_size = "5m")))
            cfg = R.AssemblyConfig(; k = 13, corrector = :iterative,
                strategy = :scalable, layout = :olc, olc_tool = tool,
                sequencing_tech = tech, olc_options = opts)
            Test.@test cfg.olc_tool == tool
        end

        # Long-read tool with a short-read tech is a mismatch.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
            corrector = :iterative, strategy = :scalable, layout = :olc,
            olc_tool = :flye, sequencing_tech = :illumina)

        # hifiasm is PacBio-HiFi-specific: a Nanopore tech is rejected.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
            corrector = :iterative, strategy = :scalable, layout = :olc,
            olc_tool = :hifiasm, sequencing_tech = :nanopore)

        # canu requires an estimated genome_size in olc_options.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
            corrector = :iterative, strategy = :scalable, layout = :olc,
            olc_tool = :canu, sequencing_tech = :pacbio)  # no genome_size → error
    end

    Test.@testset ":auto tool resolution" begin
        # :auto resolves to :megahit for a short-read tech, :flye for a long-read tech.
        cfg_short = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, layout = :olc, olc_tool = :auto,
            sequencing_tech = :illumina)
        Test.@test R._resolve_olc_tool(cfg_short) == :megahit
        for tech in (:nanopore, :pacbio)
            cfg_long = R.AssemblyConfig(; k = 13, corrector = :iterative,
                strategy = :scalable, layout = :olc, olc_tool = :auto,
                sequencing_tech = tech)
            Test.@test R._resolve_olc_tool(cfg_long) == :flye
        end

        # An explicit tool passes through unchanged.
        cfg_ms = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, layout = :olc, olc_tool = :metaspades,
            sequencing_tech = :illumina)
        Test.@test R._resolve_olc_tool(cfg_ms) == :metaspades
    end

    Test.@testset "taxonomy single-source + adapter guard" begin
        tax = R._olc_taxonomy()
        wired = (tax.short_read_tools..., tax.long_read_tools...)
        # This testset verifies what it can WITHOUT external tools: (1) :auto
        # resolves to a member of the wired-tool set for every valid tech, and
        # (2) _run_olc_tool rejects a tool not in the taxonomy. (A stronger check
        # that every wired tool has a live _run_olc_tool branch needs a wrapper
        # stub — a follow-up.)
        for tech in (tax.short_read_techs..., tax.long_read_techs...)
            cfg = R.AssemblyConfig(; k = 13, corrector = :iterative,
                strategy = :scalable, layout = :olc, olc_tool = :auto,
                sequencing_tech = tech)
            Test.@test R._resolve_olc_tool(cfg) in wired
        end
        # An unknown tool fails loud in the adapter (guards a routing bug reaching
        # _run_olc_tool with a tool that has no wrapper branch).
        cfg = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, layout = :olc, sequencing_tech = :illumina)
        Test.@test_throws ErrorException R._run_olc_tool(:bogus_tool, "unused.fastq",
            mktempdir(), cfg)

        # read-type mappers: pure functions, deterministic, cover the pacbio path
        # the gated flye run never exercises. Fail loud on an unexpected tech.
        Test.@test R._flye_read_type(:nanopore) == "nano-corr"
        Test.@test R._flye_read_type(:pacbio) == "pacbio-corr"
        Test.@test_throws ErrorException R._flye_read_type(:illumina)
        Test.@test R._canu_read_type(:nanopore) == "nanopore"
        Test.@test R._canu_read_type(:pacbio) == "pacbio"
        Test.@test_throws ErrorException R._canu_read_type(:illumina)

        # Partition invariant (single source of truth): the corrector's accepted
        # techs are exactly the taxonomy's short ∪ long, so :auto can never see a
        # tech the resolver/read-type maps don't classify.
        Test.@test Set((tax.short_read_techs..., tax.long_read_techs...)) ==
                   Set((:illumina, :ultima, :nanopore, :pacbio))
    end

    Test.@testset "read_type is a reserved olc_options key" begin
        # read_type is route-managed (derived from sequencing_tech) — a caller
        # cannot override it via olc_options.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
            corrector = :iterative, strategy = :scalable, layout = :olc,
            olc_tool = :flye, sequencing_tech = :nanopore,
            olc_options = (; read_type = "nano-hq"))
    end

    Test.@testset "discoverability + incompatibility warnings" begin
        # olc_tool/olc_options are inert unless layout=:olc → warn (mirrors the
        # output_dir/skip_solid inert-field warnings). match_mode=:any tolerates
        # other construction-time logs.
        Test.@test_logs (:warn, r"only used when layout=:olc") R.AssemblyConfig(;
            k = 13, olc_tool = :megahit)  # layout defaults to :native → inert → warn
        Test.@test_logs (:warn, r"only used when layout=:olc") R.AssemblyConfig(;
            k = 13, olc_options = (; threads = 4))  # non-empty olc_options under :native → warn
        # metaspades needs paired-end; the single-end corrected handoff will fail →
        # warn at construction (kept selectable, per follow-on paired-end support).
        Test.@test_logs (:warn, r"expects paired-end") R.AssemblyConfig(;
            k = 13, corrector = :iterative, strategy = :scalable, layout = :olc,
            olc_tool = :metaspades, sequencing_tech = :illumina)
    end

    Test.@testset "olc_options reserved-key rejection" begin
        # Route-managed keys cannot be overridden via olc_options (they would
        # redirect the assembler / break cleanup / bypass the corrected handoff).
        for k in (:fastq1, :fastq2, :fastq, :outdir, :out_dir, :executor)
            Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
                corrector = :iterative, strategy = :scalable, layout = :olc,
                sequencing_tech = :illumina,
                olc_options = NamedTuple{(k,)}(("x",)))
        end
        # A non-reserved key is accepted.
        Test.@test R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, layout = :olc, sequencing_tech = :illumina,
            olc_options = (; k_list = "21")) isa R.AssemblyConfig
    end
end
