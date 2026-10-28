# Brute-force oracle for reverse-complement handling in k-mer graphs.
#
# Included by tests with `include(joinpath(@__DIR__, "strand_oracle.jl"))`.
# Everything except `graph_counts` uses plain String operations only, so the
# expectations are independent of Mycelia, BioSequences, and Kmers.
#
# Definitions (uppercase A/C/G/T input, k-mer size k, rc = reverse complement):
#   fwd(x)     occurrences of k-mer x on the given strands
#   fwdE(a,b)  occurrences of b immediately following a in the same sequence
#   Double-strand (directed): vertices K u rc(K); coverage(x) = fwd(x) + fwd(rc(x))
#     (palindromes: 2*fwd(x)); edges E u {(rc(b), rc(a))}; weight(a,b) =
#     fwdE(a,b) + fwdE(rc(b), rc(a)), so a self-RC edge has weight 2*fwdE(a,b).
#   Canonical (undirected): canon(x) = min(x, rc(x)); coverage(c) = fwd(c) +
#     fwd(rc(c)) (palindromes: 2*fwd(c)); edges are sorted pairs
#     (canon(a), canon(b)) with weights summed over forward edges; self-loops
#     are allowed and a self-RC edge (b == rc(a)) has weight 2*fwdE(a,b).
#
# One convention for both modes: an observation that is its own reverse
# complement (a palindromic k-mer, or a self-RC edge) counts twice, once per
# strand, for vertices and edges alike.

module StrandOracle

import Random
import Mycelia
import MetaGraphsNext

const COMPLEMENT = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')

"""
    revcomp(s::String)::String

Reverse complement of an uppercase A/C/G/T string.
"""
function revcomp(s::String)::String
    return String(reverse([COMPLEMENT[c] for c in s]))
end

"""
    canon(s::String)::String

Lexicographically smaller of `s` and its reverse complement.
"""
function canon(s::String)::String
    rc = revcomp(s)
    return rc < s ? rc : s
end

"""
    fwd_counts(seqs::Vector{String}, k::Int)

Forward-strand k-mer counts and consecutive k-mer (directed edge) counts.
Returns `(vertex_counts, edge_counts)`.
"""
function fwd_counts(seqs::Vector{String}, k::Int)::Tuple{
        Dict{String, Int}, Dict{Tuple{String, String}, Int}}
    vertices = Dict{String, Int}()
    edges = Dict{Tuple{String, String}, Int}()
    for seq in seqs
        n = length(seq)
        for i in 1:(n - k + 1)
            km = seq[i:(i + k - 1)]
            vertices[km] = get(vertices, km, 0) + 1
            if i + k <= n
                nxt = seq[(i + 1):(i + k)]
                edges[(km, nxt)] = get(edges, (km, nxt), 0) + 1
            end
        end
    end
    return vertices, edges
end

"""
    doublestrand_expected(seqs::Vector{String}, k::Int)

Expected directed double-strand graph: every forward k-mer and edge is also
recorded on the reverse-complement strand. Returns `(vertex_counts, edge_counts)`.
"""
function doublestrand_expected(seqs::Vector{String},
        k::Int)::Tuple{
        Dict{String, Int}, Dict{Tuple{String, String}, Int}}
    fv, fe = fwd_counts(seqs, k)
    vertices = Dict{String, Int}()
    edges = Dict{Tuple{String, String}, Int}()
    for (x, c) in fv
        vertices[x] = get(vertices, x, 0) + c
        rx = revcomp(x)
        vertices[rx] = get(vertices, rx, 0) + c
    end
    for ((a, b), c) in fe
        edges[(a, b)] = get(edges, (a, b), 0) + c
        rev = (revcomp(b), revcomp(a))
        edges[rev] = get(edges, rev, 0) + c
    end
    return vertices, edges
end

"""
    canonical_expected(seqs::Vector{String}, k::Int)

Expected undirected canonical graph. Edge keys are sorted pairs `(min, max)`
of canonical k-mers. Returns `(vertex_counts, edge_counts)`.

Observations that are their own reverse complement are counted on both
strands: a palindromic k-mer `x` (`rc(x) == x`) contributes `2 * fwd(x)` to its
vertex, and a self-reverse-complement edge `(a, b)` with `b == rc(a)`
contributes `2 * fwdE(a, b)` to the self-loop `(canon(a), canon(a))`.
"""
function canonical_expected(seqs::Vector{String},
        k::Int)::Tuple{
        Dict{String, Int}, Dict{Tuple{String, String}, Int}}
    fv, fe = fwd_counts(seqs, k)
    vertices = Dict{String, Int}()
    edges = Dict{Tuple{String, String}, Int}()
    for (x, c) in fv
        cx = canon(x)
        multiplicity = revcomp(x) == x ? 2 : 1
        vertices[cx] = get(vertices, cx, 0) + multiplicity * c
    end
    for ((a, b), c) in fe
        ca = canon(a)
        cb = canon(b)
        key = ca <= cb ? (ca, cb) : (cb, ca)
        multiplicity = b == revcomp(a) ? 2 : 1
        edges[key] = get(edges, key, 0) + multiplicity * c
    end
    return vertices, edges
end

"""
    graph_counts(g; undirected::Bool)

Read vertex coverage and edge weights from a Mycelia Rhizomorph graph.
Vertex keys are `string(label)`; edge keys are `(string(src), string(dst))`,
sorted into `(min, max)` when `undirected` is true.
"""
function graph_counts(g; undirected::Bool)::Tuple{
        Dict{String, Int}, Dict{Tuple{String, String}, Int}}
    vertices = Dict{String, Int}()
    for label in MetaGraphsNext.labels(g)
        vertices[string(label)] = Mycelia.Rhizomorph.count_evidence(g[label])
    end
    edges = Dict{Tuple{String, String}, Int}()
    for (src, dst) in MetaGraphsNext.edge_labels(g)
        s = string(src)
        d = string(dst)
        key = (undirected && d < s) ? (d, s) : (s, d)
        edges[key] = get(edges, key, 0) + Mycelia.Rhizomorph.count_evidence(g[src, dst])
    end
    return vertices, edges
end

"""
    random_non_rc_closed(rng::Random.AbstractRNG, n::Int, k::Int)::String

Draw random length-`n` A/C/G/T strings until the k-mer set is not closed under
reverse complement (some k-mer's reverse complement is absent).
"""
function random_non_rc_closed(rng::Random.AbstractRNG, n::Int, k::Int)::String
    alphabet = collect("ACGT")
    while true
        s = String(rand(rng, alphabet, n))
        kmers = Set(s[i:(i + k - 1)] for i in 1:(n - k + 1))
        if any(x -> !(revcomp(x) in kmers), kmers)
            return s
        end
    end
end

end # module StrandOracle
