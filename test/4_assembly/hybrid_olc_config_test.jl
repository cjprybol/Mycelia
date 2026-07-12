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

        # Long-read tools are not wired in this PR (td-wvto): rejected at construction.
        for tool in (:hifiasm, :flye, :canu, :metaflye)
            Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
                corrector = :iterative, strategy = :scalable, layout = :olc,
                olc_tool = tool, sequencing_tech = :pacbio)
        end

        # :auto with a long-read tech resolves to a not-yet-wired long-read assembler.
        Test.@test_throws ErrorException R.AssemblyConfig(; k = 13,
            corrector = :iterative, strategy = :scalable, layout = :olc,
            olc_tool = :auto, sequencing_tech = :nanopore)
    end

    Test.@testset ":auto tool resolution" begin
        # :auto resolves to :megahit for a short-read tech.
        cfg_auto = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, layout = :olc, olc_tool = :auto,
            sequencing_tech = :illumina)
        Test.@test R._resolve_olc_tool(cfg_auto) == :megahit

        # An explicit tool passes through unchanged.
        cfg_ms = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, layout = :olc, olc_tool = :metaspades,
            sequencing_tech = :illumina)
        Test.@test R._resolve_olc_tool(cfg_ms) == :metaspades
    end

    Test.@testset "adapter guard for unwired tools" begin
        cfg = R.AssemblyConfig(; k = 13, corrector = :iterative,
            strategy = :scalable, layout = :olc, sequencing_tech = :illumina)
        # _run_olc_tool only wires :megahit / :metaspades; anything else fails loud
        # rather than silently doing nothing (guards against a future routing bug
        # reaching the adapter with a long-read tool before td-wvto lands).
        Test.@test_throws ErrorException R._run_olc_tool(:flye, "unused.fastq",
            mktempdir(), cfg)
    end
end
