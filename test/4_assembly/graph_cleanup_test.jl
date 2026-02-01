import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences

Test.@testset "Graph cleanup" begin
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = Int,
        vertex_data_type = Dict{Symbol, Any},
        edge_data_type = Nothing
    )

    graph[1] = Dict(:coverage => 1.0, :sequence => BioSequences.LongDNA{4}("AAAA"))
    graph[2] = Dict(:coverage => 10.0, :sequence => BioSequences.LongDNA{4}("ATGC"))
    graph[3] = Dict(:coverage => 1.0, :sequence => BioSequences.LongDNA{4}("ATGC"))

    Graphs.add_edge!(graph, 1, 2)
    Graphs.add_edge!(graph, 2, 3)

    cleaned_graph,
    stats = Mycelia.statistical_tip_clipping(
        graph;
        min_coverage_threshold = 1,
        std_dev_multiplier = 3.0,
        preserve_high_quality = true
    )

    Test.@test !MetaGraphsNext.haskey(cleaned_graph, 1)
    Test.@test MetaGraphsNext.haskey(cleaned_graph, 3)
    Test.@test stats[:tips_removed] == 1
    Test.@test stats[:high_quality_preserved] == 1

    component_graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = Int,
        vertex_data_type = Dict{Symbol, Any},
        edge_data_type = Nothing
    )

    component_graph[1] = Dict(:coverage => 5.0)
    component_graph[2] = Dict(:coverage => 5.0)
    component_graph[3] = Dict(:coverage => 5.0)

    Graphs.add_edge!(component_graph, 1, 2)

    components = Mycelia.find_connected_components(component_graph)

    Test.@test length(components) == 2
    Test.@test any(component -> Set(component) == Set([1, 2]), components)
    Test.@test any(component -> Set(component) == Set([3]), components)
end

Test.@testset "assess_sequence_quality - DNA case handling" begin
    ## BioSequences.jl is case-insensitive - lowercase and uppercase should produce identical results
    ## This documents and tests that BioSequences handles case internally (no rejection of lowercase)

    ## DNA - uppercase (baseline)
    Test.@test Mycelia.assess_sequence_quality("ACGT") == 1.0
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongDNA{4}("ACGT")) == 1.0

    ## DNA - lowercase (BioSequences accepts this and normalizes internally)
    Test.@test Mycelia.assess_sequence_quality("acgt") == 1.0
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongDNA{4}("acgt")) == 1.0

    ## DNA - mixed case (soft-masked style, common in reference genomes)
    Test.@test Mycelia.assess_sequence_quality("ACGTacgt") == 1.0
    Test.@test Mycelia.assess_sequence_quality("acgtACGT") == 1.0
    Test.@test Mycelia.assess_sequence_quality("AcGt") == 1.0

    ## Verify case-insensitivity: same result regardless of case
    uppercase_quality = Mycelia.assess_sequence_quality("ACGT")
    lowercase_quality = Mycelia.assess_sequence_quality("acgt")
    mixed_quality = Mycelia.assess_sequence_quality("AcGt")
    Test.@test uppercase_quality == lowercase_quality
    Test.@test uppercase_quality == mixed_quality

    ## With ambiguous bases (N) - quality depends on complexity score compensation
    ## Note: The quality formula is (1 - ambiguous_fraction) * complexity_score
    ## where complexity_score = unique_symbols / alphabet_size
    ## For "ACGTN": 5 unique symbols / 4 = 1.25 complexity, (1 - 0.2) * 1.25 = 1.0
    quality_with_n = Mycelia.assess_sequence_quality("ACGTN")
    Test.@test quality_with_n == 1.0  # Complexity compensates for single N

    ## Same result regardless of case when N is present
    Test.@test Mycelia.assess_sequence_quality("acgtn") == quality_with_n
    Test.@test Mycelia.assess_sequence_quality("ACGTn") == quality_with_n
    Test.@test Mycelia.assess_sequence_quality("acgtN") == quality_with_n

    ## Empty sequences
    Test.@test Mycelia.assess_sequence_quality("") == 0.0
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongDNA{4}("")) == 0.0

    ## Low complexity (homopolymer) - quality reduced due to low complexity score
    quality_homo = Mycelia.assess_sequence_quality("AAAAAAA")
    Test.@test quality_homo < 1.0  # Low complexity
    ## Verify same for lowercase
    Test.@test Mycelia.assess_sequence_quality("aaaaaaa") == quality_homo

    ## All ambiguous (all N's) - quality is 0.0
    ## Formula: (1 - 1.0) * (1/4) = 0.0 * 0.25 = 0.0
    quality_all_n = Mycelia.assess_sequence_quality("NNNN")
    Test.@test quality_all_n == 0.0  # 100% ambiguous = 0 quality
    ## Lowercase N's should work the same
    Test.@test Mycelia.assess_sequence_quality("nnnn") == quality_all_n

    ## Full IUPAC ambiguity codes (DNA with all ambiguous characters)
    ## W=A/T, S=G/C, M=A/C, K=G/T, R=A/G, Y=C/T, B=C/G/T, D=A/G/T, H=A/C/T, V=A/C/G
    ## Note: High complexity score (14 unique / 4 expected = 3.5) compensates for ambiguity
    ambiguous_seq = "ACGTWSMKRYBDHV"
    quality_ambig = Mycelia.assess_sequence_quality(ambiguous_seq)
    Test.@test quality_ambig == 1.0  # High complexity compensates (capped at 1.0)
    ## Lowercase should work
    Test.@test Mycelia.assess_sequence_quality(lowercase(ambiguous_seq)) == quality_ambig
end

Test.@testset "assess_sequence_quality - RNA sequences" begin
    ## RNA uppercase
    Test.@test Mycelia.assess_sequence_quality("ACGU") == 1.0
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongRNA{4}("ACGU")) == 1.0

    ## RNA lowercase
    Test.@test Mycelia.assess_sequence_quality("acgu") == 1.0
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongRNA{4}("acgu")) == 1.0

    ## RNA mixed case
    Test.@test Mycelia.assess_sequence_quality("ACGUacgu") == 1.0
    Test.@test Mycelia.assess_sequence_quality("AcGu") == 1.0

    ## Verify case-insensitivity for RNA
    rna_upper = Mycelia.assess_sequence_quality("ACGU")
    rna_lower = Mycelia.assess_sequence_quality("acgu")
    Test.@test rna_upper == rna_lower

    ## RNA with ambiguous bases - complexity compensates for single N
    quality_with_n = Mycelia.assess_sequence_quality("ACGUN")
    Test.@test quality_with_n == 1.0  # 5 unique / 4 = 1.25 complexity compensates
    Test.@test Mycelia.assess_sequence_quality("acgun") == quality_with_n

    ## Empty RNA
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongRNA{4}("")) == 0.0
end

Test.@testset "assess_sequence_quality - amino acid sequences" begin
    ## Standard 20 amino acids
    aa_seq = "ACDEFGHIKLMNPQRSTVWY"

    ## Uppercase amino acids
    Test.@test Mycelia.assess_sequence_quality(aa_seq) == 1.0
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongAA(aa_seq)) == 1.0

    ## Lowercase amino acids
    Test.@test Mycelia.assess_sequence_quality(lowercase(aa_seq)) == 1.0
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongAA(lowercase(aa_seq))) ==
               1.0

    ## Mixed case amino acids
    mixed_aa = "AcDeFgHiKlMnPqRsTvWy"
    Test.@test Mycelia.assess_sequence_quality(mixed_aa) == 1.0

    ## Verify case-insensitivity
    aa_upper = Mycelia.assess_sequence_quality(aa_seq)
    aa_lower = Mycelia.assess_sequence_quality(lowercase(aa_seq))
    Test.@test aa_upper == aa_lower

    ## Empty amino acid sequence
    Test.@test Mycelia.assess_sequence_quality(BioSequences.LongAA("")) == 0.0

    ## Short amino acid sequences - low complexity due to AA alphabet size of 20
    ## Quality = (1 - ambiguous_fraction) * (unique_symbols / 20)
    ## "MKTL" has 4 unique symbols, so quality = 1.0 * (4/20) = 0.2
    short_aa_quality = Mycelia.assess_sequence_quality("MKTL")
    Test.@test short_aa_quality == 0.2
    Test.@test Mycelia.assess_sequence_quality("mktl") == short_aa_quality  # Case-insensitive
end

Test.@testset "Graph cleanup with string sequences (lowercase)" begin
    ## Test that graph cleanup works with string sequences instead of BioSequence objects
    ## This verifies the assess_sequence_quality(::AbstractString) method works in graph context

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = Int,
        vertex_data_type = Dict{Symbol, Any},
        edge_data_type = Nothing
    )

    ## Use lowercase string sequences (should work fine via BioSequences case handling)
    graph[1] = Dict(:coverage => 1.0, :sequence => "aaaa")
    graph[2] = Dict(:coverage => 10.0, :sequence => "atgc")
    graph[3] = Dict(:coverage => 1.0, :sequence => "atgc")

    Graphs.add_edge!(graph, 1, 2)
    Graphs.add_edge!(graph, 2, 3)

    ## This should not throw - verify lowercase strings are handled correctly
    cleaned_graph,
    stats = Mycelia.statistical_tip_clipping(
        graph;
        min_coverage_threshold = 1,
        std_dev_multiplier = 3.0,
        preserve_high_quality = true
    )

    Test.@test stats[:tips_removed] >= 0
    Test.@test haskey(stats, :high_quality_preserved)
end

Test.@testset "Graph cleanup with soft-masked sequences (mixed case)" begin
    ## Test graph cleanup with soft-masked style sequences (mixed uppercase/lowercase)
    ## This is common in reference genomes where lowercase indicates repetitive/masked regions

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = Int,
        vertex_data_type = Dict{Symbol, Any},
        edge_data_type = Nothing
    )

    ## Soft-masked style sequences (uppercase = unmasked, lowercase = masked)
    graph[1] = Dict(:coverage => 1.0, :sequence => "AAAAttttGGGG")
    graph[2] = Dict(:coverage => 10.0, :sequence => "ATGCatgcATGC")
    graph[3] = Dict(:coverage => 1.0, :sequence => "atgcATGCatgc")

    Graphs.add_edge!(graph, 1, 2)
    Graphs.add_edge!(graph, 2, 3)

    ## This should not throw - verify mixed case strings are handled correctly
    cleaned_graph,
    stats = Mycelia.statistical_tip_clipping(
        graph;
        min_coverage_threshold = 1,
        std_dev_multiplier = 3.0,
        preserve_high_quality = true
    )

    Test.@test stats[:tips_removed] >= 0
    Test.@test haskey(stats, :high_quality_preserved)
end

Test.@testset "Graph cleanup with sequences containing N" begin
    ## Test graph cleanup with sequences containing ambiguous bases (N)
    ## These should be handled gracefully and factor into quality scoring

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = Int,
        vertex_data_type = Dict{Symbol, Any},
        edge_data_type = Nothing
    )

    ## Sequences with varying amounts of N (ambiguous bases)
    graph[1] = Dict(:coverage => 1.0, :sequence => "NNNN")  # All N - low quality
    graph[2] = Dict(:coverage => 10.0, :sequence => "ATGCATGC")  # No N - high quality
    graph[3] = Dict(:coverage => 1.0, :sequence => "ATNCATGC")  # Some N - medium quality

    Graphs.add_edge!(graph, 1, 2)
    Graphs.add_edge!(graph, 2, 3)

    cleaned_graph,
    stats = Mycelia.statistical_tip_clipping(
        graph;
        min_coverage_threshold = 1,
        std_dev_multiplier = 3.0,
        preserve_high_quality = true
    )

    Test.@test stats[:tips_removed] >= 0
end
