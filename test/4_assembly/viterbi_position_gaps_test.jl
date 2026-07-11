# Tier-2 opt-in per-position Viterbi margin (td-4osf / calibration).
#
# Verifies that ViterbiCorrectionConfig(record_position_gaps=true):
#   1. populates result.paths[i].diagnostics[:position_gaps] (a Vector{Float64}),
#   2. is ABSENT by default (flag off), and
#   3. does NOT perturb the decode — the path + score are byte-identical with the
#      flag on vs off, so the telemetry is provably side-effect-free.
#
# Run:
#   julia --project=. test/4_assembly/viterbi_position_gaps_test.jl

import Test
import Mycelia
import MetaGraphsNext

Test.@testset "Viterbi per-position gap telemetry (record_position_gaps)" begin
    function create_weighted_graph(labels)
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type = eltype(collect(labels)),
            vertex_data_type = Any,
            edge_data_type = Mycelia.Rhizomorph.StrandWeightedEdgeData,
            weight_function = Mycelia.Rhizomorph.edge_data_weight,
            default_weight = 0.0
        )
        for label in labels
            graph[label] = nothing
        end
        return graph
    end
    function add_edge!(graph, s, t, w)
        graph[s, t] = Mycelia.Rhizomorph.StrandWeightedEdgeData(
            w, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward)
        return graph
    end
    path_labels(p) = [step.vertex_label for step in p.steps]

    # A graph with genuine competing paths (S->A weight 9 vs S->B weight 4; both
    # reach T), so the frontier has >=2 states per depth and the margin is finite.
    function ambiguous_graph()
        g = create_weighted_graph(["S", "A", "B", "T", "X"])
        add_edge!(g, "S", "A", 9.0)
        add_edge!(g, "S", "B", 4.0)
        add_edge!(g, "A", "T", 1.0)
        add_edge!(g, "A", "X", 9.0)
        add_edge!(g, "B", "T", 1.0)
        return g
    end

    obs = [["S", "A", "X"]]

    # --- Default: no gaps recorded ------------------------------------------
    base = Mycelia.correct_observations(ambiguous_graph(), obs)
    Test.@test !haskey(base.paths[1].diagnostics, :position_gaps)

    # --- Opt-in: gaps recorded ----------------------------------------------
    cfg = Mycelia.ViterbiCorrectionConfig(record_position_gaps = true)
    withgaps = Mycelia.correct_observations(ambiguous_graph(), obs; config = cfg)
    gaps = withgaps.paths[1].diagnostics[:position_gaps]
    Test.@test gaps isa Vector{Float64}
    Test.@test !isempty(gaps)                       # >=1 decode depth for a 3-unit observation
    Test.@test all(g -> g >= 0.0, gaps)             # margin = best - 2nd-best >= 0 (or Inf)
    Test.@test all(g -> !isnan(g), gaps)

    # --- Byte-identity: the telemetry must not change the decode ------------
    Test.@test path_labels(something(withgaps.paths[1].path)) ==
               path_labels(something(base.paths[1].path))
    Test.@test withgaps.paths[1].score == base.paths[1].score

    # --- _top2_score_gap unit behavior --------------------------------------
    Test.@test Mycelia._top2_score_gap(Dict(:a => 5.0, :b => 2.0, :c => -1.0)) == 3.0
    Test.@test Mycelia._top2_score_gap(Dict(:a => 5.0)) == Inf   # no competitor => maximally confident
    Test.@test Mycelia._top2_score_gap(Dict{Symbol, Float64}()) == Inf
end
