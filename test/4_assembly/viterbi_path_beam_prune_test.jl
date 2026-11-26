# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_path_beam_prune_test.jl")'
# ```
#
# Direct coverage of `Rhizomorph._prune_viterbi_beam` (td-63vh).

import Mycelia
import Test

Test.@testset "_prune_viterbi_beam (td-63vh)" begin
    Forward = Mycelia.Rhizomorph.Forward
    Reverse = Mycelia.Rhizomorph.Reverse

    # Rows are (vertex, strand, score, predecessor vertex, predecessor strand,
    # transition probability). Insertion order is the row order, so a second
    # call with the reversed rows checks that ties do not follow dict order.
    function frontier(rows)
        state_type = Tuple{String, Mycelia.Rhizomorph.StrandOrientation}
        scores = Dict{state_type, Float64}()
        predecessors = Dict{state_type, Tuple{state_type, Float64}}()
        for (vertex, strand, score, pred_vertex, pred_strand, probability) in rows
            state = (vertex, strand)
            scores[state] = score
            predecessors[state] = ((pred_vertex, pred_strand), probability)
        end
        return scores, predecessors
    end

    function assert_retained(kept_scores, kept_predecessors, scores, predecessors, expected)
        Test.@test Set(keys(kept_scores)) == Set(expected)
        Test.@test Set(keys(kept_predecessors)) == Set(expected)
        for state in expected
            Test.@test kept_scores[state] == scores[state]
            Test.@test kept_predecessors[state] == predecessors[state]
        end
        return nothing
    end

    function prune(scores, predecessors, beam_width, beam_log_threshold)
        return Mycelia.Rhizomorph._prune_viterbi_beam(
            scores,
            predecessors,
            beam_width,
            beam_log_threshold,
        )
    end

    # "A" is the lexicographically best vertex and the worst score, so a
    # name-first or insertion-order ranking would keep it and fail these tests.
    ranked_rows = [
        ("A", Forward, -5.0, "S", Forward, 0.1),
        ("B", Forward, -1.0, "S", Reverse, 0.2),
        ("C", Forward, 0.0, "S", Forward, 0.3),
        ("D", Forward, -2.0, "S", Reverse, 0.4),
    ]

    Test.@testset "top-k keeps the highest scores and aligned predecessors" begin
        scores, predecessors = frontier(ranked_rows)
        kept_scores, kept_predecessors = prune(scores, predecessors, 2, Inf)
        assert_retained(
            kept_scores,
            kept_predecessors,
            scores,
            predecessors,
            [("C", Forward), ("B", Forward)],
        )
        Test.@test !haskey(kept_scores, ("D", Forward))
        Test.@test !haskey(kept_predecessors, ("A", Forward))

        kept_one, kept_one_pred = prune(scores, predecessors, 1, Inf)
        assert_retained(
            kept_one,
            kept_one_pred,
            scores,
            predecessors,
            [("C", Forward)],
        )

        # Width above the frontier is a no-op when the threshold is disabled.
        kept_all, kept_all_pred = prune(scores, predecessors, 8, Inf)
        Test.@test kept_all == scores
        Test.@test kept_all_pred == predecessors

        # A runner-up inside any score margin is still dropped when the width
        # is already full. Threshold Inf disables the margin.
        Test.@test haskey(scores, ("B", Forward))
        Test.@test !haskey(kept_one, ("B", Forward))
    end

    Test.@testset "log-score threshold drops states beyond the margin" begin
        # Differences from the best score: 0, 1, 2, 2.5, 5.
        # The comparison is strict (`best - score > threshold`), so a state
        # sitting exactly on the margin is kept.
        rows = [
            ("best", Forward, 0.0, "S", Forward, 0.9),
            ("near", Forward, -1.0, "S", Forward, 0.8),
            ("margin", Forward, -2.0, "S", Reverse, 0.7),
            ("past", Forward, -2.5, "S", Forward, 0.6),
            ("far", Forward, -5.0, "S", Reverse, 0.5),
        ]
        scores, predecessors = frontier(rows)

        kept_scores, kept_predecessors = prune(scores, predecessors, 10, 2.0)
        assert_retained(
            kept_scores,
            kept_predecessors,
            scores,
            predecessors,
            [("best", Forward), ("near", Forward), ("margin", Forward)],
        )
        Test.@test !haskey(kept_scores, ("past", Forward))
        Test.@test !haskey(kept_predecessors, ("far", Forward))

        # The margin binds inside a width that could have held every state.
        tight_scores, tight_predecessors = prune(scores, predecessors, 10, 1.0)
        assert_retained(
            tight_scores,
            tight_predecessors,
            scores,
            predecessors,
            [("best", Forward), ("near", Forward)],
        )
        Test.@test length(tight_scores) < 10
    end

    Test.@testset "equal scores break ties by vertex string then strand" begin
        # Sort keys `(string(vertex), Int(strand))`, ascending after score:
        # ("A", Forward)=("A", 1), ("M", Reverse)=("M", 0),
        # ("M", Forward)=("M", 1), ("Z", Reverse)=("Z", 0).
        # Reverse is 0 and Forward is 1, so Reverse wins a same-vertex tie.
        tied_rows = [
            ("Z", Reverse, 0.0, "S", Forward, 0.1),
            ("M", Forward, 0.0, "S", Reverse, 0.2),
            ("M", Reverse, 0.0, "S", Forward, 0.3),
            ("A", Forward, 0.0, "S", Reverse, 0.4),
        ]
        expected_top2 = [("A", Forward), ("M", Reverse)]

        for rows in (tied_rows, reverse(tied_rows))
            scores, predecessors = frontier(rows)
            kept_scores, kept_predecessors = prune(scores, predecessors, 2, Inf)
            assert_retained(
                kept_scores,
                kept_predecessors,
                scores,
                predecessors,
                expected_top2,
            )
            Test.@test kept_predecessors[("A", Forward)][2] == 0.4
            Test.@test kept_predecessors[("M", Reverse)][2] == 0.3
        end

        strand_scores, strand_preds = frontier([
            ("K", Forward, -0.5, "S", Forward, 0.11),
            ("K", Reverse, -0.5, "S", Reverse, 0.22),
        ])
        kept_strand, kept_strand_pred = prune(strand_scores, strand_preds, 1, Inf)
        assert_retained(
            kept_strand,
            kept_strand_pred,
            strand_scores,
            strand_preds,
            [("K", Reverse)],
        )

        vertex_scores, vertex_preds = frontier([
            ("B", Forward, -0.5, "S", Forward, 0.11),
            ("A", Forward, -0.5, "S", Reverse, 0.22),
        ])
        kept_vertex, kept_vertex_pred = prune(vertex_scores, vertex_preds, 1, Inf)
        assert_retained(
            kept_vertex,
            kept_vertex_pred,
            vertex_scores,
            vertex_preds,
            [("A", Forward)],
        )

        # A higher score outranks a better sort key. The remaining slot then
        # follows the sort key inside the tied group.
        mixed_rows = [
            ("A", Forward, -1.0, "S", Forward, 0.11),
            ("Z", Reverse, 0.0, "S", Reverse, 0.22),
            ("M", Reverse, -1.0, "S", Forward, 0.33),
        ]
        mixed_scores, mixed_preds = frontier(mixed_rows)
        kept_mixed, kept_mixed_pred = prune(mixed_scores, mixed_preds, 2, Inf)
        assert_retained(
            kept_mixed,
            kept_mixed_pred,
            mixed_scores,
            mixed_preds,
            [("Z", Reverse), ("A", Forward)],
        )

        # Threshold 0 keeps only the best-score group, then the width picks
        # among that tie by sort key. "M" is below the margin.
        margin_rows = [
            ("Z", Reverse, 0.0, "S", Forward, 0.1),
            ("B", Forward, 0.0, "S", Reverse, 0.2),
            ("A", Forward, 0.0, "S", Forward, 0.3),
            ("M", Reverse, -1.0, "S", Reverse, 0.4),
        ]
        margin_scores, margin_preds = frontier(margin_rows)
        kept_margin, kept_margin_pred = prune(margin_scores, margin_preds, 2, 0.0)
        assert_retained(
            kept_margin,
            kept_margin_pred,
            margin_scores,
            margin_preds,
            [("A", Forward), ("B", Forward)],
        )
        Test.@test !haskey(kept_margin, ("Z", Reverse))
        Test.@test !haskey(kept_margin_pred, ("M", Reverse))
    end

    # Callers pass a non-empty frontier and a positive width. These two cases
    # pin the helper's current behaviour when that guarantee is not met.
    Test.@testset "empty frontier throws and a non-positive width returns empty dicts" begin
        empty_scores, empty_predecessors = frontier([])
        local threw = false
        local message = ""
        try
            prune(empty_scores, empty_predecessors, 2, Inf)
        catch err
            threw = true
            Test.@test err isa BoundsError
            message = sprint(showerror, err)
        end
        Test.@test threw
        Test.@test occursin("attempt to access 0-element", message)
        Test.@test occursin("at index [1]", message)

        scores, predecessors = frontier(ranked_rows)
        for width in (0, -1)
            kept_scores, kept_predecessors = prune(scores, predecessors, width, Inf)
            Test.@test isempty(kept_scores)
            Test.@test isempty(kept_predecessors)
        end
    end
end
