# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/find_linear_path_set_membership_test.jl")'
# ```
#
# Guards the O(L) Vector-membership -> O(1) Set-membership fix in
# `find_linear_path` (src/rhizomorph/algorithms/contigs.jl, td-kokn). The fix is
# a PURE PERFORMANCE change: contig extraction must return byte-identical
# contigs, in identical order, with identical boundaries. This file:
#   (1) asserts `find_contigs_next` output is unchanged vs a reference
#       implementation that replicates the OLD Vector-membership logic, on
#       branchy / cyclic / multi-contig toy graphs; and
#   (2) profiles `find_contigs_next` timing at increasing contig lengths for
#       BOTH the current (Set) and reference (Vector) implementations, and
#       asserts the current implementation's log-log slope (alpha) drops toward
#       linear (~1) while the reference stays super-linear (~2).

import Test
import Mycelia
import Graphs
import MetaGraphsNext
import Statistics

const _Rz = Mycelia.Rhizomorph

# ---------------------------------------------------------------------------
# Reference implementation replicating the OLD (pre-fix) Vector-membership
# `find_linear_path`. Uses the SAME src helpers (neighbor queries, dominance
# resolution) so the ONLY difference from the current src is Vector vs Set
# membership — isolating exactly the change under test.
# ---------------------------------------------------------------------------
function _reference_find_linear_path(graph, start_vertex, visited::Set)
    if start_vertex in visited
        return Vector{typeof(start_vertex)}()
    end

    path = [start_vertex]
    current = start_vertex

    while true
        out_neighbors = _Rz.get_outgoing_neighbors(graph, current)
        valid_neighbors = [n for n in out_neighbors if !(n in visited) && !(n in path)]

        next_vertex = _Rz._select_linear_neighbor(graph, current, valid_neighbors; direction = :out)
        next_vertex === nothing && break
        _Rz._edge_is_dominant_among_incoming(graph, current, next_vertex) || break

        push!(path, next_vertex)
        current = next_vertex
    end

    current = start_vertex
    backward_path = Vector{typeof(start_vertex)}()

    while true
        in_neighbors = _Rz.get_incoming_neighbors(graph, current)
        valid_neighbors = [n
                           for n in in_neighbors
                           if !(n in visited) && !(n in path) && !(n in backward_path)]

        prev_vertex = _Rz._select_linear_neighbor(graph, current, valid_neighbors; direction = :in)
        prev_vertex === nothing && break
        _Rz._edge_is_dominant_among_outgoing(graph, prev_vertex, current) || break

        pushfirst!(backward_path, prev_vertex)
        current = prev_vertex
    end

    return [backward_path; path]
end

# Reference `find_contigs_next` mirroring the src driver loop but calling the
# reference (Vector-membership) linear path. Returns the same ContigPath vector
# shape so results can be compared field-by-field.
function _reference_find_contigs_next(graph; min_contig_length::Int = 500)
    labels = collect(MetaGraphsNext.labels(graph))
    isempty(labels) && return Mycelia.Rhizomorph.ContigPath{Any, Any}[]

    contigs = Vector{Mycelia.Rhizomorph.ContigPath}()
    visited = Set{eltype(labels)}()

    for start_vertex in labels
        start_vertex in visited && continue
        path = _reference_find_linear_path(graph, start_vertex, visited)
        if !isempty(path)
            sequence = _Rz.generate_contig_sequence(graph, path)
            coverage_profile = _Rz.generate_coverage_profile(graph, path)
            if length(sequence) >= min_contig_length
                push!(contigs, Mycelia.Rhizomorph.ContigPath(path, sequence, coverage_profile))
            end
            for vertex in path
                push!(visited, vertex)
            end
        end
    end

    return _Rz.sort_contigs_by_length(contigs)
end

# ---------------------------------------------------------------------------
# Toy-graph builders (String-labeled graphs; direct construction).
# ---------------------------------------------------------------------------
function _new_string_graph()
    return MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
        edge_data_type = Mycelia.Rhizomorph.StringEdgeData,
    )
end

function _add_vertices!(graph, labels)
    for v in labels
        graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
    end
    return graph
end

function _add_edge!(graph, src, dst)
    graph[src, dst] = Mycelia.Rhizomorph.StringEdgeData(1)
    return graph
end

# Linear chain v0001 -> v0002 -> ... of `n` vertices (single long contig; the
# worst case for the path-membership walk).
function _linear_chain_graph(n::Int)
    graph = _new_string_graph()
    labels = [string("v", lpad(i, 6, '0')) for i in 1:n]
    _add_vertices!(graph, labels)
    for i in 1:(n - 1)
        _add_edge!(graph, labels[i], labels[i + 1])
    end
    return graph
end

# Compare two ContigPath vectors field-by-field for byte-identity.
function _contigs_identical(a, b)
    length(a) == length(b) || return false
    for (ca, cb) in zip(a, b)
        ca.vertices == cb.vertices || return false
        string(ca.sequence) == string(cb.sequence) || return false
        ca.coverage_profile == cb.coverage_profile || return false
        ca.length == cb.length || return false
    end
    return true
end

Test.@testset "find_linear_path Set-membership (td-kokn)" begin
    Test.@testset "byte-identical contigs vs Vector-membership reference" begin
        graphs = Dict{String, Any}()

        # (a) single linear chain
        graphs["linear_chain"] = _linear_chain_graph(30)

        # (b) cycle with a tail: A->B->C->A, plus D->A. Path membership must stop
        #     the forward walk at the cycle closure and the backward walk at the
        #     already-in-path vertex.
        let g = _new_string_graph()
            _add_vertices!(g, ["A", "B", "C", "D"])
            _add_edge!(g, "A", "B"); _add_edge!(g, "B", "C"); _add_edge!(g, "C", "A")
            _add_edge!(g, "D", "A")
            graphs["cycle_with_tail"] = g
        end

        # (c) bubble: A->B->D and A->C->D (reconverging branch).
        let g = _new_string_graph()
            _add_vertices!(g, ["A", "B", "C", "D"])
            _add_edge!(g, "A", "B"); _add_edge!(g, "A", "C")
            _add_edge!(g, "B", "D"); _add_edge!(g, "C", "D")
            graphs["bubble"] = g
        end

        # (d) two disjoint linear components (multi-contig; visited carries across).
        let g = _new_string_graph()
            _add_vertices!(g, ["a1", "a2", "a3", "b1", "b2", "b3", "b4"])
            _add_edge!(g, "a1", "a2"); _add_edge!(g, "a2", "a3")
            _add_edge!(g, "b1", "b2"); _add_edge!(g, "b2", "b3"); _add_edge!(g, "b3", "b4")
            graphs["two_components"] = g
        end

        # (e) figure-eight: shared middle vertex on two loops, stresses membership.
        let g = _new_string_graph()
            _add_vertices!(g, ["m", "p", "q", "r", "s"])
            _add_edge!(g, "m", "p"); _add_edge!(g, "p", "q"); _add_edge!(g, "q", "m")
            _add_edge!(g, "m", "r"); _add_edge!(g, "r", "s"); _add_edge!(g, "s", "m")
            graphs["figure_eight"] = g
        end

        for (name, g) in graphs
            current = Mycelia.Rhizomorph.find_contigs_next(g; min_contig_length = 1)
            reference = _reference_find_contigs_next(g; min_contig_length = 1)
            Test.@test _contigs_identical(current, reference)
        end
    end

    Test.@testset "timing alpha drops from super-linear toward linear" begin
        sizes = [8000, 16000, 32000]
        graphs = [_linear_chain_graph(n) for n in sizes]
        # Each chain is a single contig of length N; the seed is its head, so
        # find_linear_path walks all N vertices. Time find_linear_path DIRECTLY
        # (not find_contigs_next) to isolate the membership walk — the naive
        # per-char String concatenation in generate_contig_sequence is itself
        # O(N^2) for String graphs (irrelevant to production k-mer graphs) and
        # would otherwise mask the membership scaling under test.
        seeds = [first(collect(MetaGraphsNext.labels(g))) for g in graphs]

        # Warm up (JIT) on the smallest graph so timings measure the algorithm,
        # not first-call compilation.
        Mycelia.Rhizomorph.find_linear_path(graphs[1], seeds[1], Set{String}())
        _reference_find_linear_path(graphs[1], seeds[1], Set{String}())

        # Minimum over a few trials suppresses GC/scheduler noise so the log-log
        # slope reflects the algorithm's intrinsic scaling.
        best_time(f) = minimum(@elapsed(f()) for _ in 1:5)

        t_current = Float64[]
        t_reference = Float64[]
        for (g, seed) in zip(graphs, seeds)
            push!(t_current, best_time(() -> Mycelia.Rhizomorph.find_linear_path(g, seed, Set{String}())))
            push!(t_reference, best_time(() -> _reference_find_linear_path(g, seed, Set{String}())))
        end

        # log-log slope (alpha): time ~ N^alpha.
        alpha(xs, ys) = begin
            lx = log.(Float64.(xs)); ly = log.(ys)
            mx = Statistics.mean(lx); my = Statistics.mean(ly)
            sum((lx .- mx) .* (ly .- my)) / sum((lx .- mx) .^ 2)
        end

        a_current = alpha(sizes, t_current)
        a_reference = alpha(sizes, t_reference)

        println("\n[find_linear_path Set-membership] contig-extraction scaling")
        println("  N sizes           : ", sizes)
        println("  Set (current)  t  : ", round.(t_current; sigdigits = 3), " s  alpha=", round(a_current; digits = 2))
        println("  Vector (old)   t  : ", round.(t_reference; sigdigits = 3), " s  alpha=", round(a_reference; digits = 2))

        # The current (Set) implementation is near-linear and clearly below the
        # Vector reference's super-linear (~2) slope. Thresholds are loose enough
        # to tolerate timing noise while still proving the scaling improvement:
        # the separation between the two slopes is the load-bearing assertion.
        Test.@test a_current < 1.5
        Test.@test a_reference > 1.6
        Test.@test a_reference - a_current > 0.4
    end
end
