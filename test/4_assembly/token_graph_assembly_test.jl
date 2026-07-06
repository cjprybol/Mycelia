# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/token_graph_assembly_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/token_graph_assembly_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Token Graph Assembly - end-to-end through assemble_genome
#
# Exercises the TokenGraph route: pre-tokenized SentencePiece/word-token
# sequences are supplied via AssemblyConfig.token_sequences and assembled
# directly (reads are ignored). SingleStrand only.

import Test
import Mycelia
import MetaGraphsNext

Test.@testset "Token Graph Assembly (assemble_genome route)" begin
    Test.@testset "TokenGraph enum member exists" begin
        Test.@test Mycelia.Rhizomorph.TokenGraph isa Mycelia.Rhizomorph.AssemblyMethod
    end

    Test.@testset "config validation: token_sequences requires SingleStrand" begin
        toks = [["the", "cat", "sat"]]
        # DoubleStrand (the default) must be rejected for tokens.
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(;
            token_sequences = toks,
            graph_mode = Mycelia.Rhizomorph.DoubleStrand)
        # SingleStrand constructs cleanly and stores the tokens.
        cfg_ok = Mycelia.Rhizomorph.AssemblyConfig(;
            token_sequences = toks,
            graph_mode = Mycelia.Rhizomorph.SingleStrand)
        Test.@test cfg_ok.token_sequences == toks
    end

    Test.@testset "default config has token_sequences === nothing" begin
        cfg = Mycelia.Rhizomorph.AssemblyConfig()
        Test.@test cfg.token_sequences === nothing
    end

    Test.@testset "end-to-end linear token assembly" begin
        # A single linear chain: the -> quick -> brown -> fox -> jumps.
        token_sequences = [["the", "quick", "brown", "fox", "jumps"]]

        config = Mycelia.Rhizomorph.AssemblyConfig(;
            token_sequences = token_sequences,
            graph_mode = Mycelia.Rhizomorph.SingleStrand)

        # reads are ignored when token_sequences is supplied; pass an empty vector.
        result = Mycelia.Rhizomorph.assemble_genome(String[], config)

        Test.@test result isa Mycelia.Rhizomorph.AssemblyResult
        Test.@test Mycelia.Rhizomorph.has_graph_structure(result)
        Test.@test result.graph isa MetaGraphsNext.MetaGraph

        # 5 distinct tokens -> 5 vertices, 4 edges.
        Test.@test length(collect(MetaGraphsNext.labels(result.graph))) == 5
        Test.@test length(collect(MetaGraphsNext.edge_labels(result.graph))) == 4

        # At least one contig, names aligned, and stats stamped.
        Test.@test !isempty(result.contigs)
        Test.@test length(result.contigs) == length(result.contig_names)
        Test.@test result.assembly_stats["method"] == "TokenGraph"
        Test.@test result.assembly_stats["num_input_sequences"] == 1

        # A valid linear assembly recovers the full ordered token string via the
        # Eulerian fast path, joined by the default separator.
        Test.@test "the quick brown fox jumps" in result.contigs

        # Structural consistency check passes.
        report = Mycelia.Rhizomorph.validate_assembly_structure(result)
        Test.@test isempty(report["issues"])
    end

    Test.@testset "shared tokens across multiple sequences" begin
        # Two sequences sharing a prefix token exercise the multi-sequence path
        # and the branching-graph contig fallback.
        token_sequences = [
            ["the", "cat", "sat"],
            ["the", "dog", "ran"]
        ]

        config = Mycelia.Rhizomorph.AssemblyConfig(;
            token_sequences = token_sequences,
            graph_mode = Mycelia.Rhizomorph.SingleStrand)

        result = Mycelia.Rhizomorph.assemble_genome(String[], config)

        Test.@test result isa Mycelia.Rhizomorph.AssemblyResult
        Test.@test Mycelia.Rhizomorph.has_graph_structure(result)
        # 5 distinct tokens: the, cat, sat, dog, ran.
        Test.@test length(collect(MetaGraphsNext.labels(result.graph))) == 5
        Test.@test !isempty(result.contigs)
        Test.@test length(result.contigs) == length(result.contig_names)
        # Contigs are token strings; every emitted token appears somewhere.
        joined = join(result.contigs, " ")
        for tok in ["the", "cat", "sat", "dog", "ran"]
            Test.@test occursin(tok, joined)
        end
    end
end
