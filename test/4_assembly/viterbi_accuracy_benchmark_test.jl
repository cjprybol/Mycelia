# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_accuracy_benchmark_test.jl")'
# ```

import CSV
import DataFrames
import MetaGraphsNext
import Mycelia
import Statistics
import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "viterbi_accuracy_benchmark.jl"))

# AMBIENT-ENVIRONMENT GUARD -- TERMINAL, and deliberately not a `Test.@test`.
#
# k is bound at load time from this process's argv/ENV, so an exported
# VITERBI_ACCURACY_K (or a stray --k on the test runner's command line) would
# silently change WHICH parameterization this suite measures while every
# assertion below still passed. This suite asserts the k=11 numbers.
#
# It must ERROR rather than record a failure, because a `Test.@test` is
# advisory: execution continues past it. The `main()` refusal testset below
# asserts that `main(["--k", "9", ...])` throws -- which it does only when 9
# disagrees with the load-time k. Under an exported VITERBI_ACCURACY_K=9 the
# disagreement vanishes, no error is raised, and `main` falls through to TWO
# full B8 executions that write into benchmarking/results/ inside the repo. A
# guard whose downstream consumer depends on it being enforcing cannot be
# fail-open, so this runs before any testset and stops the file.
VITERBI_ACCURACY_K == 11 || error(
    "This suite asserts the k=11 numbers but loaded k=$(VITERBI_ACCURACY_K) " *
    "(source: $(VITERBI_ACCURACY_K_SOURCE)). Refusing to run: the assertions " *
    "below would measure a different parameterization, and the main() refusal " *
    "test would execute the benchmark instead of throwing. Unset " *
    "VITERBI_ACCURACY_K and remove any --k from the test command line.")

Test.@testset "B8 Viterbi accuracy benchmark artifacts" begin
    Test.@test VITERBI_ACCURACY_K == 11
    Test.@test VITERBI_ACCURACY_K_SOURCE == "default"

    output_dir = mktempdir(prefix = "viterbi_accuracy_b8_test_")
    artifacts = run_viterbi_accuracy_benchmark(output_dir; write_plots = false)
    summary = CSV.read(artifacts.summary_csv, DataFrames.DataFrame)

    Test.@test artifacts.rows == 9
    Test.@test DataFrames.nrow(summary) == 9
    Test.@test Set(summary.dataset_id) == Set([
        "pstvd_nc002030", "phix174_nc001422", "austen_pride_prejudice_excerpt"
    ])
    Test.@test Set(summary.target_error_rate) == Set([0.01, 0.05, 0.10])
    Test.@test all(summary.baseline_edit_distance .> 0)
    Test.@test all(summary.corrected_edit_distance .<= summary.baseline_edit_distance)
    Test.@test all(summary.edit_distance_reduction .>= 0)
    Test.@test all(summary.injected_error_recall .>= 0)
    Test.@test isfile(artifacts.index)
    Test.@test isfile(artifacts.provenance)

    # k is stamped into the emitted rows and into the run id, so a reader of the
    # .csv alone (the manuscript repo reads a hand-copied .csv, not the sidecar)
    # can tell which parameterization produced it. Asserted against the LITERAL
    # 11, never against VITERBI_ACCURACY_K — comparing the column to the constant
    # that filled it is self-referential and would pass at any k, including a k
    # inherited from an ambient environment variable.
    Test.@test "benchmark_k" in DataFrames.names(summary)
    Test.@test all(summary.benchmark_k .== 11)
    Test.@test all(endswith.(summary.benchmark_run_id, "_k11"))

    # Control A — over-correction on un-corrupted input.
    overcorrection = CSV.read(artifacts.overcorrection_csv, DataFrames.DataFrame)
    Test.@test artifacts.overcorrection_rows == 9
    Test.@test DataFrames.nrow(overcorrection) == 9
    Test.@test all(overcorrection.injected_error_count .== 0)
    Test.@test all(overcorrection.over_correction_edit_distance .>= 0)
    Test.@test all(0.0 .<= overcorrection.over_correction_rate .<= 1.0)
    Test.@test "benchmark_k" in DataFrames.names(overcorrection)
    Test.@test all(overcorrection.benchmark_k .== 11)
    # Cross-field accounting invariant: a row has zero changed observations iff it
    # has zero over-correction edit distance. Exercises the numerator logic even
    # though every observed row is the (expected) zero case — a mis-indexed or
    # inverted count would break this equivalence. (review: over-correction path)
    Test.@test all(
        (overcorrection.changed_observations .== 0) .==
        (overcorrection.over_correction_edit_distance .== 0))
    # The specificity claim the manuscript needs: the corrector does not damage
    # already-correct input (abstains). Asserted as an aggregate so a single
    # tiny-fixture edit does not flake the suite.
    Test.@test Statistics.mean(overcorrection.over_correction_rate) < 0.01

    # Control B — shuffled-weight + random-rewire nulls.
    null_control = CSV.read(artifacts.null_control_csv, DataFrames.DataFrame)
    Test.@test artifacts.null_control_rows == 9
    Test.@test DataFrames.nrow(null_control) == 9
    Test.@test all(null_control.null_seed .== 20260711)
    Test.@test "benchmark_k" in DataFrames.names(null_control)
    Test.@test all(null_control.benchmark_k .== 11)
    for col in (:real_injected_error_recall, :weight_null_injected_error_recall,
        :rewire_null_injected_error_recall)
        Test.@test all(0.0 .<= null_control[!, col] .<= 1.0)
    end
    # Robust invariant: neither null may recover injected errors BETTER than the
    # real graph (a degraded graph cannot help more than the true one).
    Test.@test Statistics.mean(null_control.real_injected_error_recall) >=
               Statistics.mean(null_control.weight_null_injected_error_recall) - 1e-9
    Test.@test Statistics.mean(null_control.real_injected_error_recall) >=
               Statistics.mean(null_control.rewire_null_injected_error_recall) - 1e-9
    # POSITIVE CONTROL: without this, the "null <= real" invariant is vacuous —
    # it also passes if the corrector itself regressed to zero recovery (real ==
    # null == 0). The real graph must actually recover on these fixtures, and it
    # must never fall into the decode-collapse fallback.
    Test.@test Statistics.mean(null_control.real_injected_error_recall) > 0.5
    Test.@test all(null_control.real_decoded)
    # Fixture-scoped regression guard on the DISCRIMINATING signal: on these
    # pinned fixtures (RefSeq NC_002030.1 / NC_001422.1 + fixed text excerpt,
    # deterministic injector, fixed seed) destroying the true adjacency abolishes
    # recovery entirely. If pointed at new fixtures, re-baseline this one line;
    # the portable "null <= real" invariant above carries over unchanged.
    Test.@test Statistics.mean(null_control.rewire_null_injected_error_recall) == 0.0
    Test.@test all(.!null_control.rewire_null_decoded)
end

Test.@testset "k selection is validated, never silently defaulted" begin
    # The parser takes `args`/`env` as parameters precisely so these rules can be
    # exercised without launching a process. Every case below is a wrong-answer
    # path, not a usability nicety: k selects the fixtures, the run id, the output
    # directory and the recorded command line, so a silent fallback publishes one
    # parameterization's numbers under another's identity.
    no_env = Dict{String, String}()

    # Value forms resolve, and report which form won.
    Test.@test _viterbi_accuracy_k_at_load(["--k", "13"], no_env) == 13
    Test.@test _viterbi_accuracy_k_at_load(["--k=13"], no_env) == 13
    Test.@test _viterbi_accuracy_resolve_k(["--k", "13"], no_env).source == "--k"
    Test.@test _viterbi_accuracy_resolve_k(["--k=13"], no_env).source == "--k="

    # A bare trailing `--k` must ERROR. The realistic trigger is an unquoted
    # `--k $K` with K unset: the empty word disappears and `--k` lands last.
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(["--k"], no_env)
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(
        ["--output-dir", "/tmp/b8", "--skip-plots", "--k"], no_env)

    # Non-integer and out-of-range values error rather than defaulting.
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(["--k", "abc"], no_env)
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(["--k", "0"], no_env)
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(["--k", "-3"], no_env)
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(["--k=abc"], no_env)
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(["--k=0"], no_env)
    # `--k` followed by the NEXT flag is a missing value, not a value of "--skip-plots".
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(["--k", "--skip-plots"], no_env)

    # No selection at all is the only path to the default.
    Test.@test _viterbi_accuracy_k_at_load(String[], no_env) == 11
    Test.@test _viterbi_accuracy_resolve_k(String[], no_env).source == "default"
    Test.@test _viterbi_accuracy_k_at_load(["--skip-plots"], no_env) == 11

    # Environment variable: honoured, validated, and outranked by the flag.
    Test.@test _viterbi_accuracy_k_at_load(String[], Dict("VITERBI_ACCURACY_K" => "9")) == 9
    Test.@test _viterbi_accuracy_resolve_k(
        String[], Dict("VITERBI_ACCURACY_K" => "9")).source == "VITERBI_ACCURACY_K"
    Test.@test _viterbi_accuracy_k_at_load(
        ["--k", "13"], Dict("VITERBI_ACCURACY_K" => "9")) == 13
    # SET BUT EMPTY is an explicit selection of "", not the default path — the
    # form an unset shell variable produces under `VITERBI_ACCURACY_K="$K"`.
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(
        String[], Dict("VITERBI_ACCURACY_K" => ""))
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(
        String[], Dict("VITERBI_ACCURACY_K" => "abc"))
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(
        String[], Dict("VITERBI_ACCURACY_K" => "0"))
    # Empty INLINE form -- what unquoted `--k=$K` produces with K unset.
    Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(
        ["--k="], Dict{String, String}())
    # A REPEATED --k must error rather than resolve. Left unguarded this is
    # first-wins, so a wrapper appending --k to a command that already carries
    # one would silently discard the caller's value and then record a
    # single-flag command line that was never typed. Both spellings count, in
    # every combination, because precedence was previously form-dependent
    # rather than positional.
    for duplicate in (
        ["--k", "9", "--k", "13"],
        ["--k=9", "--k=13"],
        ["--k=9", "--k", "13"],
        ["--k", "13", "--k=9"])
        Test.@test_throws ErrorException _viterbi_accuracy_k_at_load(
            duplicate, Dict{String, String}())
    end
end

Test.@testset "value-taking flags are validated, never silently defaulted" begin
    # --output-dir is held to the same standard as --k. Each case below was
    # measured against the previous implementation and silently succeeded:
    # a bare trailing flag fell back to the default directory, an empty value
    # resolved through abspath("") to the CWD (the repository root for the
    # documented invocation), and a following flag was consumed as the path
    # while itself being dropped.
    Test.@test_throws ErrorException _viterbi_accuracy_arg_value(
        ["--skip-plots", "--output-dir"], "--output-dir", "/tmp/fallback")
    Test.@test_throws ErrorException _viterbi_accuracy_arg_value(
        ["--output-dir", ""], "--output-dir", "/tmp/fallback")
    Test.@test_throws ErrorException _viterbi_accuracy_arg_value(
        ["--output-dir", "   "], "--output-dir", "/tmp/fallback")
    Test.@test_throws ErrorException _viterbi_accuracy_arg_value(
        ["--output-dir", "--skip-plots"], "--output-dir", "/tmp/fallback")
    # Absent flag still yields the default; a real value is returned verbatim.
    Test.@test _viterbi_accuracy_arg_value(
        ["--skip-plots"], "--output-dir", "/tmp/fallback") == "/tmp/fallback"
    Test.@test _viterbi_accuracy_arg_value(
        ["--output-dir", "/tmp/real"], "--output-dir", "/tmp/fallback") ==
               "/tmp/real"
end

Test.@testset "default output directory is k-stamped unconditionally" begin
    # Every k gets its own tree, INCLUDING the default 11. The earlier k == 11
    # special case pointed the default invocation at the committed historical
    # directory, so the most common command overwrote git-tracked artifacts with
    # no --force and no backup. Expected values are pinned literally: comparing
    # the function to a re-derivation of itself would pass under that mutation.
    Test.@test basename(_viterbi_accuracy_default_output_dir(11)) ==
               "viterbi_accuracy_b8_k11"
    Test.@test basename(_viterbi_accuracy_default_output_dir(9)) ==
               "viterbi_accuracy_b8_k9"
    Test.@test _viterbi_accuracy_default_output_dir(9) !=
               _viterbi_accuracy_default_output_dir(11)
    # The historical directory is frozen: no k resolves onto it.
    Test.@test basename(VITERBI_ACCURACY_HISTORICAL_OUTPUT_DIR) == "viterbi_accuracy_b8"
    Test.@test _viterbi_accuracy_default_output_dir(11) !=
               VITERBI_ACCURACY_HISTORICAL_OUTPUT_DIR
    Test.@test _viterbi_accuracy_default_output_dir(9) !=
               VITERBI_ACCURACY_HISTORICAL_OUTPUT_DIR
end

Test.@testset "main() refuses a --k it cannot honour" begin
    # k is bound at load time from the process's argv/ENV, so a `--k` inside a
    # caller-supplied args vector arrives too late. Silently dropping it would
    # return k=11 tables to a caller who asked for 9 while `--output-dir` from
    # the SAME vector was honoured. The check runs before any benchmark work.
    Test.@test_throws ErrorException main(["--k", "9", "--skip-plots"])
    Test.@test_throws ErrorException main(["--k=9", "--skip-plots"])
    # Same validation as at load time, so a bad value cannot slip through here.
    Test.@test_throws ErrorException main(["--k"])
    Test.@test_throws ErrorException main(["--k", "abc"])
end

Test.@testset "recorded command line reflects how k was selected" begin
    # The suite runs with no --k and no VITERBI_ACCURACY_K, so the provenance
    # command line must NOT assert a --k the operator never typed.
    Test.@test VITERBI_ACCURACY_K_SOURCE == "default"
    Test.@test !("--k" in _viterbi_accuracy_command_args())
    Test.@test !any(startswith.(_viterbi_accuracy_command_args(), "--k="))
end

Test.@testset "shuffled-weight null preserves topology + weight multiset" begin
    # _shuffle_weighted_graph_weights must permute weights only: same edge set,
    # same multiset of weights, so any recovery collapse is attributable purely
    # to the weight->edge reassignment, not to a changed graph.
    fixture = first(viterbi_accuracy_fixtures())
    config = Mycelia.ViterbiCorrectionConfig(error_rate = 0.05)
    weighted = Mycelia.build_correction_weighted_graph(fixture.graph; config = config)
    before_edges = Set(collect(MetaGraphsNext.edge_labels(weighted)))
    before_weights = sort([weighted[s, d].weight for (s, d) in before_edges])

    shuffled = _shuffle_weighted_graph_weights(weighted, 20260711)
    after_edges = Set(collect(MetaGraphsNext.edge_labels(shuffled)))
    after_weights = sort([shuffled[s, d].weight for (s, d) in after_edges])

    Test.@test before_edges == after_edges
    Test.@test before_weights ≈ after_weights
    # Aliasing guard: the shuffle must NOT mutate the caller's graph (it deepcopies
    # internally), so the null can never leak into the real/Control-A arms even if
    # a future fixture ships pre-weighted. Verify `weighted` is untouched and the
    # returned graph is a distinct object. (review: aliasing-landmine)
    Test.@test shuffled !== weighted
    Test.@test sort([weighted[s, d].weight
                     for (s, d) in
                         MetaGraphsNext.edge_labels(weighted)]) ≈ before_weights
    # build_correction_weighted_graph must itself return a fresh graph for these
    # fixtures (the precondition that makes the in-place null isolated).
    Test.@test Mycelia.build_correction_weighted_graph(fixture.graph; config = config) !==
               fixture.graph
end

Test.@testset "random-rewire null preserves vertices/edge-count/weight multiset" begin
    # _random_rewired_weighted_graph must keep the same vertex set, the same edge
    # COUNT, and the same weight multiset — only the topology (which vertices an
    # edge joins) is randomized, so a recovery collapse is attributable purely to
    # the destroyed adjacency, not to a smaller/heavier graph.
    fixture = first(viterbi_accuracy_fixtures())
    config = Mycelia.ViterbiCorrectionConfig(error_rate = 0.05)
    weighted = Mycelia.build_correction_weighted_graph(fixture.graph; config = config)
    before_vertices = Set(collect(MetaGraphsNext.labels(weighted)))
    before_edges = collect(MetaGraphsNext.edge_labels(weighted))
    before_weights = sort([weighted[s, d].weight for (s, d) in before_edges])

    rewired = _random_rewired_weighted_graph(weighted, 20260711)
    after_vertices = Set(collect(MetaGraphsNext.labels(rewired)))
    after_edges = collect(MetaGraphsNext.edge_labels(rewired))
    after_weights = sort([rewired[s, d].weight for (s, d) in after_edges])

    Test.@test before_vertices == after_vertices
    Test.@test length(after_edges) == length(before_edges)
    Test.@test before_weights ≈ after_weights
    # Topology genuinely changed (not a no-op permutation): at least some edges
    # differ from the real graph. (With hundreds of vertices this is ~certain.)
    Test.@test Set(after_edges) != Set(before_edges)
    # No self-loops: the rewire rejects src == dst, so no vertex points to itself.
    # (The vertex/edge-count/weight tests above would all still pass with a
    # self-loop, so assert it explicitly.)
    Test.@test all(s != d for (s, d) in after_edges)
end
