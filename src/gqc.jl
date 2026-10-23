# gqc v0.1 — reusable graph QC packet
# Spec: docs/gqc-v0.1-spec.md (in-repo). Semver 0.1.0.

const GQC_VERSION = "0.1.0"
const GQC_N_SHUFFLES = 100
const GQC_SHUFFLE_SEED = 20261021
const GQC_MIN_SHUFFLE_MEAN = 1.0
## Default BH family is k=11 only (within-graph). Paper cross-genus FDR is a post-pass.
const GQC_DEFAULT_BH_KS = [11]
## Frozen node TSV header for both `nodes/<slug>.tsv` and `formats/nodes.tsv`.
const GQC_NODE_TSV_COLUMNS = (
    "graph_id", "node_id", "index", "length", "depth",
    "degree", "in_degree", "out_degree", "sequence",
)
## Frozen summary.csv column order (every field always emitted).
const GQC_SUMMARY_COLUMNS = (
    "graph_id", "k", "n_segments", "n_links",
    "total_bp", "total_segment_bp", "density",
    "n_components", "largest_component_fraction", "beta0", "beta1",
    "n_tips", "n_linear", "n_branch", "n_isolated", "n_sources", "n_sinks",
    "n_unitigs", "unitig_n50", "unitig_l50", "unitig_max",
    "n_bubbles", "n_tangles", "n_self_loops",
    "cov_q0", "cov_q25", "cov_q50", "cov_q75", "cov_q100",
    "ref_kmer_fraction", "ref_path_fraction", "ref_absent_kmers", "ref_breaks", "ref_walkability",
    "dinucleotide_observed_repeats", "dinucleotide_null_mean", "dinucleotide_null_std",
    "dinucleotide_shuffle_ratio", "dinucleotide_z_score", "dinucleotide_shuffle_p_upper",
    "dinucleotide_null_degenerate", "dinucleotide_bh_q_within_graph",
    "codon_observed_repeats", "codon_null_mean", "codon_null_std",
    "codon_shuffle_ratio", "codon_z_score", "codon_shuffle_p_upper",
    "codon_null_degenerate",
    "null_source", "null_applicable",
    "strand_mode", "canonical", "gqc_version",
)

# =============================================================================
# Null / multiple-testing helpers (verbatim from rhizomorph-paper
# notebooks/viral-graphs.jl @ c8285d0, Oct 22 2026)
# =============================================================================

"""
    dinucleotide_shuffle(sequence::String, rng::Random.AbstractRNG)::String

Altschul–Erickson dinucleotide shuffle: preserves dinucleotide counts and the
first/last base. `U` is treated as `T` (RNA). Non-ACGT letters (e.g. `N`) are
included in the AE alphabet so they do not abort; they are preserved as symbols
in the Eulerian trail.
"""
function dinucleotide_shuffle(sequence::String, rng::Random.AbstractRNG)::String
    chars = collect(replace(uppercase(sequence), 'U' => 'T'))
    n = length(chars)
    n <= 2 && return String(chars)
    alphabet = sort(unique(chars))
    final = chars[end]
    exits = Dict{Char, Vector{Char}}(c => Char[] for c in alphabet)
    for i in 1:(n - 1)
        push!(exits[chars[i]], chars[i + 1])
    end
    tree_bases = [c for c in alphabet if c != final && !isempty(exits[c])]
    last_exit = Dict{Char, Char}()
    max_steps = max(length(alphabet), 4)
    reaches_final(c) = begin
        steps = 0
        while c != final && steps <= max_steps
            haskey(last_exit, c) || return false
            c = last_exit[c]
            steps += 1
        end
        c == final
    end
    attempts = 0
    while true
        attempts += 1
        for c in tree_bases
            last_exit[c] = rand(rng, exits[c])
        end
        all(reaches_final, tree_bases) && break
        ## Extremely constrained alphabets: fall through with original order
        attempts > 10_000 && return String(chars)
    end
    for c in alphabet
        if c in tree_bases
            deleteat!(exits[c], findfirst(==(last_exit[c]), exits[c]))
            Random.shuffle!(rng, exits[c])
            push!(exits[c], last_exit[c])
        else
            Random.shuffle!(rng, exits[c])
        end
    end
    next_exit = Dict(c => 1 for c in alphabet)
    shuffled = Vector{Char}(undef, n)
    shuffled[1] = chars[1]
    for i in 2:n
        previous = shuffled[i - 1]
        shuffled[i] = exits[previous][next_exit[previous]]
        next_exit[previous] += 1
    end
    return String(shuffled)
end

"""
    dinucleotide_shuffle(sequence::BioSequences.LongDNA, rng::Random.AbstractRNG)

Shuffle a DNA BioSequence via the string Altschul–Erickson implementation.
"""
function dinucleotide_shuffle(sequence::BioSequences.LongDNA, rng::Random.AbstractRNG)
    return BioSequences.LongDNA{4}(dinucleotide_shuffle(string(sequence), rng))
end

"""
    benjamini_hochberg(p::Vector{Float64})::Vector{Float64}

Benjamini–Hochberg q-values (step-up, monotone). Ported verbatim from
rhizomorph-paper `viral-graphs.jl`.
"""
function benjamini_hochberg(p::Vector{Float64})::Vector{Float64}
    m = length(p)
    order = sortperm(p)
    q = similar(p)
    running = 1.0
    for rank in m:-1:1
        i = order[rank]
        running = min(running, p[i] * m / rank)
        q[i] = running
    end
    return q
end

"""
    gqc_cross_bh(summary_paths; k=11, p_key="dinucleotide_shuffle_p_upper")

Cross-genus (cross-graph) Benjamini–Hochberg on per-graph p-values at a single `k`.

# Arguments
- `summary_paths`: iterable of `summary.json` paths (one packet each)

# Returns
`Vector` of `NamedTuple`s `(graph_id, k, p, q)` sorted by ascending `p`.
Packet column `dinucleotide_bh_q_within_graph` is **not** this — callers must
run this (or equivalent) for paper-style FDR across genera.
"""
function gqc_cross_bh(summary_paths;
                      k::Int=11,
                      p_key::AbstractString="dinucleotide_shuffle_p_upper")
    rows = NamedTuple{(:graph_id, :k, :p, :q), Tuple{String, Int, Float64, Float64}}[]
    collected = Tuple{String, Float64}[]
    for path in summary_paths
        data = JSON.parsefile(string(path))
        for row in data
            Int(row["k"]) == k || continue
            p = row[p_key]
            p === nothing && continue
            push!(collected, (string(row["graph_id"]), Float64(p)))
        end
    end
    isempty(collected) && return rows
    ps = Float64[p for (_, p) in collected]
    qs = benjamini_hochberg(ps)
    order = sortperm(ps)
    for i in order
        gid, p = collected[i]
        push!(rows, (graph_id=gid, k=k, p=p, q=qs[i]))
    end
    return rows
end

"""
    codon_shuffle(sequence::String, rng::Random.AbstractRNG)::String

Codon-preserving shuffle: permute in-frame codons, leave a short overhang
unchanged. Used as a second null for repeat-excess alongside dinucleotide shuffle.
"""
function codon_shuffle(sequence::String, rng::Random.AbstractRNG)::String
    n = length(sequence)
    n < 3 && return sequence
    n_codons = div(n, 3)
    overhang = n - 3 * n_codons
    codons = [sequence[(3(i - 1) + 1):(3i)] for i in 1:n_codons]
    Random.shuffle!(rng, codons)
    return join(codons) * sequence[(n - overhang + 1):n]
end

# =============================================================================
# GFA parsing (public records API)
# =============================================================================

"""
GFA segment (S line) with optional coverage tags.
"""
struct GFASegment
    id::String
    sequence::String
    length::Int
    depth::Union{Float64, Nothing}
    read_count::Union{Int, Nothing}
    kmer_count::Union{Int, Nothing}
end

"""
GFA link (L line) with orientations (`from_forward` / `to_forward` == `+`).
"""
struct GFALink
    from::String
    from_forward::Bool
    to::String
    to_forward::Bool
    overlap::String
end

"""
Parsed GFA 1.0 graph used by gqc.
"""
struct GFAGraph
    segments::OrderedCollections.OrderedDict{String, GFASegment}
    links::Vector{GFALink}
    headers::Vector{String}
    path::String
end

"""
    read_gfa_records(gfa_file::AbstractString) -> GFAGraph

Parse a GFA 1.0 file into segments and links. Reads `LN:i`, `DP:f`, `RC:i`,
`KC:i` when present. Sequence length is taken from the sequence field (or `LN:i`
if the sequence is `*`).
"""
function read_gfa_records(gfa_file::AbstractString)::GFAGraph
    segments = OrderedCollections.OrderedDict{String, GFASegment}()
    links = GFALink[]
    headers = String[]
    for line in eachline(gfa_file)
        isempty(line) && continue
        fields = split(line, '\t')
        isempty(fields) && continue
        tag = fields[1]
        if tag == "H"
            push!(headers, line)
        elseif tag == "S" && length(fields) >= 3
            sid = String(fields[2])
            seq = String(fields[3])
            depth = nothing
            read_count = nothing
            kmer_count = nothing
            ln = nothing
            for f in fields[4:end]
                if startswith(f, "DP:f:")
                    depth = parse(Float64, f[6:end])
                elseif startswith(f, "RC:i:")
                    read_count = parse(Int, f[6:end])
                elseif startswith(f, "KC:i:")
                    kmer_count = parse(Int, f[6:end])
                elseif startswith(f, "LN:i:")
                    ln = parse(Int, f[6:end])
                end
            end
            seq_len = seq == "*" ? something(ln, 0) : length(seq)
            segments[sid] = GFASegment(sid, seq, seq_len, depth, read_count, kmer_count)
        elseif tag == "L" && length(fields) >= 6
            push!(links, GFALink(
                String(fields[2]),
                fields[3] == "+",
                String(fields[4]),
                fields[5] == "+",
                String(fields[6]),
            ))
        end
    end
    return GFAGraph(segments, links, headers, String(gfa_file))
end

# =============================================================================
# Topology / size metrics
# =============================================================================

function _gqc_undirected_adjacency(g::GFAGraph)
    ids = collect(keys(g.segments))
    index = Dict(id => i for (i, id) in enumerate(ids))
    n = length(ids)
    adj = [Int[] for _ in 1:n]
    self_loops = 0
    for link in g.links
        haskey(index, link.from) || continue
        haskey(index, link.to) || continue
        u = index[link.from]
        v = index[link.to]
        if u == v
            self_loops += 1
            continue
        end
        push!(adj[u], v)
        push!(adj[v], u)
    end
    for i in 1:n
        unique!(adj[i])
    end
    return ids, adj, self_loops
end

function _gqc_connected_components(adj::Vector{Vector{Int}})
    n = length(adj)
    seen = falses(n)
    comps = Vector{Vector{Int}}()
    for start in 1:n
        seen[start] && continue
        stack = Int[start]
        seen[start] = true
        comp = Int[]
        while !isempty(stack)
            u = pop!(stack)
            push!(comp, u)
            for v in adj[u]
                if !seen[v]
                    seen[v] = true
                    push!(stack, v)
                end
            end
        end
        push!(comps, comp)
    end
    return comps
end

"""
Undirected degree for each segment id (self-loops excluded from adjacency).
"""
function _gqc_degrees(adj::Vector{Vector{Int}})
    return [length(nbrs) for nbrs in adj]
end

function _gqc_node_roles(adj::Vector{Vector{Int}})
    ## Directed view is unavailable from undirected adj alone; classify on
    ## undirected degree: tip deg==1, linear deg==2, branch/join deg>=3, isolated 0.
    tips = Int[]
    linear = Int[]
    branch = Int[]
    isolated = Int[]
    for (i, d) in enumerate(_gqc_degrees(adj))
        if d == 0
            push!(isolated, i)
        elseif d == 1
            push!(tips, i)
        elseif d == 2
            push!(linear, i)
        else
            push!(branch, i)
        end
    end
    return (; tips, linear, branch, isolated)
end

"""
Directed in/out degree from GFA links (orientation-aware endpoints).
"""
function _gqc_directed_degrees(g::GFAGraph)
    ids = collect(keys(g.segments))
    outdeg = Dict(id => 0 for id in ids)
    indeg = Dict(id => 0 for id in ids)
    for link in g.links
        haskey(outdeg, link.from) || continue
        haskey(indeg, link.to) || continue
        outdeg[link.from] += 1
        indeg[link.to] += 1
    end
    sources = [id for id in ids if outdeg[id] >= 1 && indeg[id] == 0]
    sinks = [id for id in ids if indeg[id] >= 1 && outdeg[id] == 0]
    return indeg, outdeg, sources, sinks
end

function _gqc_simple_bubbles(adj::Vector{Vector{Int}}; max_len::Int=8)
    ## Count pairs of internally-disjoint paths of length <= max_len between
    ## the same endpoints where both endpoints have undirected degree >= 3
    ## or the start has out-branching in the undirected sense (deg>=3).
    n = length(adj)
    count = 0
    for s in 1:n
        length(adj[s]) < 2 && continue
        ## BFS limited depth; record parents for path reconstruction count
        dist = fill(typemax(Int), n)
        dist[s] = 0
        queue = Int[s]
        parents = [Int[] for _ in 1:n]
        while !isempty(queue)
            u = popfirst!(queue)
            dist[u] >= max_len && continue
            for v in adj[u]
                if dist[v] == typemax(Int)
                    dist[v] = dist[u] + 1
                    push!(parents[v], u)
                    push!(queue, v)
                elseif dist[v] == dist[u] + 1
                    push!(parents[v], u)
                end
            end
        end
        for t in 1:n
            t == s && continue
            length(parents[t]) < 2 && continue
            length(adj[t]) < 2 && continue
            count += 1
        end
    end
    ## Each bubble counted twice (once from each end) roughly — divide carefully.
    ## Prefer counting unique unordered endpoint pairs via a set during a second pass.
    pairs = Set{Tuple{Int,Int}}()
    for s in 1:n
        length(adj[s]) < 2 && continue
        dist = fill(typemax(Int), n)
        dist[s] = 0
        queue = Int[s]
        parents = [Int[] for _ in 1:n]
        while !isempty(queue)
            u = popfirst!(queue)
            dist[u] >= max_len && continue
            for v in adj[u]
                if dist[v] == typemax(Int)
                    dist[v] = dist[u] + 1
                    push!(parents[v], u)
                    push!(queue, v)
                elseif dist[v] == dist[u] + 1
                    push!(parents[v], u)
                end
            end
        end
        for t in (s + 1):n
            length(parents[t]) >= 2 && length(adj[t]) >= 2 && push!(pairs, (s, t))
        end
    end
    return length(pairs)
end

function _gqc_tangle_count(adj::Vector{Vector{Int}}; degree_threshold::Int=4)
    return count(d -> d >= degree_threshold, _gqc_degrees(adj))
end

"""
Parse a GFA L-line overlap CIGAR (e.g. `30M`, `30`) to an integer bp overlap.
"""
function _gqc_parse_overlap(overlap::AbstractString; default::Int=0)::Int
    m = match(r"^(\d+)", string(overlap))
    m === nothing && return default
    return parse(Int, m.captures[1])
end

"""
Undirected (min,max) segment-id → overlap bp lookup from GFA L lines.
"""
function _gqc_overlap_lookup(g::GFAGraph; default_overlap::Int=0)
    lookup = Dict{Tuple{String, String}, Int}()
    for link in g.links
        a, b = minmax(link.from, link.to)
        ov = _gqc_parse_overlap(link.overlap; default=default_overlap)
        lookup[(a, b)] = max(get(lookup, (a, b), 0), ov)
    end
    return lookup
end

"""
Span length of a segment path: `sum(seg_len) - sum(overlap on consecutive links)`.
"""
function _gqc_path_span_bp(g::GFAGraph,
                           ids::Vector{String},
                           path::Vector{Int},
                           overlap_lookup::Dict{Tuple{String, String}, Int};
                           default_overlap::Int=0)::Int
    isempty(path) && return 0
    total = g.segments[ids[path[1]]].length
    for i in 2:length(path)
        u = ids[path[i - 1]]
        v = ids[path[i]]
        ov = get(overlap_lookup, minmax(u, v), default_overlap)
        total += g.segments[ids[path[i]]].length - ov
    end
    return max(total, 0)
end

"""
Default overlap for missing L-line CIGARs: `L-1` on uniform-length (k-mer-node)
graphs, otherwise `0`.
"""
function _gqc_default_overlap(seg_seqs::Vector{String})::Int
    isempty(seg_seqs) && return 0
    lens = length.(seg_seqs)
    all(==(first(lens)), lens) || return 0
    L = first(lens)
    return L > 1 ? L - 1 : 0
end

"""
Compact non-branching paths (undirected deg==2 chains) into unitigs.

Returns overlap-corrected unitig lengths in bp:
`sum(segment lengths) - sum(link overlaps along the path)`.
Missing overlaps use `default_overlap` (typically `k-1` for k-mer-node graphs).
"""
function _gqc_unitig_lengths(g::GFAGraph, ids::Vector{String}, adj::Vector{Vector{Int}};
                             default_overlap::Int=0,
                             overlap_lookup::Union{Nothing, Dict{Tuple{String, String}, Int}}=nothing)
    lookup = overlap_lookup === nothing ? _gqc_overlap_lookup(g; default_overlap=default_overlap) : overlap_lookup
    n = length(ids)
    deg = _gqc_degrees(adj)
    visited = falses(n)
    lengths = Int[]
    function finish!(path::Vector{Int})
        push!(lengths, _gqc_path_span_bp(g, ids, path, lookup; default_overlap=default_overlap))
        return nothing
    end
    for start in 1:n
        visited[start] && continue
        ## start unitig at tips or branch points, or any unvisited isolated/linear
        if deg[start] == 2
            continue
        end
        path = Int[start]
        visited[start] = true
        for nbr0 in adj[start]
            visited[nbr0] && continue
            prev = start
            cur = nbr0
            while true
                push!(path, cur)
                visited[cur] = true
                deg[cur] != 2 && break
                nexts = [v for v in adj[cur] if v != prev]
                isempty(nexts) && break
                prev = cur
                cur = nexts[1]
                visited[cur] && break
            end
        end
        finish!(path)
    end
    ## remaining pure cycles / unvisited linear components
    for start in 1:n
        visited[start] && continue
        path = Int[start]
        visited[start] = true
        prev = 0
        cur = start
        while true
            nexts = [v for v in adj[cur] if v != prev && !visited[v]]
            isempty(nexts) && break
            prev = cur
            cur = nexts[1]
            push!(path, cur)
            visited[cur] = true
            deg[cur] != 2 && break
        end
        finish!(path)
    end
    return lengths
end

function _gqc_n50_l50(lengths::Vector{Int})
    isempty(lengths) && return (n50=0, l50=0, max_len=0, total=0)
    sorted = sort(lengths; rev=true)
    total = sum(sorted)
    half = total / 2
    cum = 0
    for (i, L) in enumerate(sorted)
        cum += L
        if cum >= half
            return (n50=L, l50=i, max_len=first(sorted), total=total)
        end
    end
    return (n50=last(sorted), l50=length(sorted), max_len=first(sorted), total=total)
end

function _gqc_coverage_quantiles(g::GFAGraph)
    depths = Float64[s.depth for s in values(g.segments) if s.depth !== nothing]
    isempty(depths) && return Dict{String, Float64}()
    qs = (0.0, 0.25, 0.5, 0.75, 1.0)
    d = sort(depths)
    out = Dict{String, Float64}()
    for q in qs
        idx = clamp(Int(round(q * (length(d) - 1))) + 1, 1, length(d))
        out["cov_q$(Int(round(100q)))"] = d[idx]
    end
    return out
end

# =============================================================================
# Repeat excess vs nulls (paper formulas)
# =============================================================================

"""
Canonical DNA k-mer as string (lexicographic min of kmer and RC).
"""
function _gqc_canonical_kmer(kmer::AbstractString)
    rc = string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(kmer)))
    return kmer <= rc ? String(kmer) : rc
end

"""
    gqc_repeat_counts(sequences; k, circular=false, canonical=true)

Return `(kmer_positions, distinct_kmers, repeats)` where
`repeats = kmer_positions - distinct_kmers` (paper definition).
Each sequence is processed separately then counts are summed (multi-segment).
"""
function gqc_repeat_counts(sequences::Vector{String};
                           k::Int,
                           circular::Bool=false,
                           canonical::Bool=true)
    kmer_positions = 0
    distinct = Set{String}()
    for seq in sequences
        n = length(seq)
        n < k && continue
        n_pos = circular ? n : (n - k + 1)
        n_pos <= 0 && continue
        for i in 1:n_pos
            if circular
                ## wrap-around kmer
                if i + k - 1 <= n
                    raw = seq[i:(i + k - 1)]
                else
                    raw = seq[i:n] * seq[1:(k - (n - i + 1))]
                end
            else
                raw = seq[i:(i + k - 1)]
            end
            ## Count positions only for valid ACGT windows (matches reference concordance).
            occursin(r"[^ACGTacgt]", raw) && continue
            kmer_positions += 1
            raw = uppercase(raw)
            push!(distinct, canonical ? _gqc_canonical_kmer(raw) : raw)
        end
    end
    distinct_kmers = length(distinct)
    repeats = kmer_positions - distinct_kmers
    return (; kmer_positions, distinct_kmers, repeats)
end

"""
    gqc_checked_dinucleotide_shuffles(sequences, rng; n_shuffles=100)

Shuffle each segment independently `n_shuffles` times. Preserves dinucleotide
counts and end bases when a shuffle is produced. When no distinct shuffle exists
(unique Eulerian trail, AT repeats, short motifs), falls back to the original
sequence and/or allows repeated nulls — never aborts. Returns
`(shuffles, null_degenerate)`.
"""
function gqc_checked_dinucleotide_shuffles(sequences::Vector{String},
                                           rng::Random.AbstractRNG;
                                           n_shuffles::Int=GQC_N_SHUFFLES,
                                           max_attempts_per_seq::Int=1000)
    function dinuc_counts(s::String)
        c = Dict{String, Int}()
        for i in 1:(length(s) - 1)
            d = s[i:(i + 1)]
            c[d] = get(c, d, 0) + 1
        end
        return c
    end
    ## Normalize like dinucleotide_shuffle (U→T, uppercase) for comparisons.
    function normalize_seq(seq::String)
        return replace(uppercase(seq), 'U' => 'T')
    end
    ## Homopolymers / length≤2 have no distinct Altschul–Erickson shuffle.
    function shuffleable(seq::String)
        s = normalize_seq(seq)
        length(s) > 2 && length(unique(s)) > 1
    end
    results = Vector{Vector{String}}(undef, n_shuffles)
    null_degenerate = false
    for t in 1:n_shuffles
        shuffled_segs = String[]
        for seq in sequences
            norm = normalize_seq(seq)
            s = norm
            if shuffleable(seq)
                found = false
                for _ in 1:max_attempts_per_seq
                    candidate = dinucleotide_shuffle(seq, rng)
                    dinuc_counts(candidate) == dinuc_counts(norm) || continue
                    (isempty(candidate) || candidate[1] == norm[1]) || continue
                    (isempty(candidate) || candidate[end] == norm[end]) || continue
                    if candidate != norm
                        s = candidate
                        found = true
                        break
                    end
                end
                if !found
                    ## Unique trail / constrained sequence: keep original (homopolymer path).
                    s = norm
                    null_degenerate = true
                end
            else
                s = dinucleotide_shuffle(seq, rng)
                if length(norm) >= 2
                    if !(dinuc_counts(s) == dinuc_counts(norm) &&
                         s[1] == norm[1] && s[end] == norm[end])
                        s = norm
                        null_degenerate = true
                    end
                end
            end
            push!(shuffled_segs, s)
        end
        results[t] = shuffled_segs
    end
    ## Allow repeated nulls when the shuffle space is small; mark degenerate.
    if any(shuffleable, sequences)
        joined = [join(segs, '|') for segs in results]
        if length(unique(joined)) < n_shuffles
            null_degenerate = true
        end
    else
        ## Entirely non-shuffleable (homopolymers / length≤2): still degenerate.
        null_degenerate = true
    end
    return results, null_degenerate
end

"""
    gqc_repeat_null_stats(sequences; k, n_shuffles, shuffle_seed, salt, circular, canonical, min_shuffle_mean)

Paper formulas:
- `null_repeats` per shuffle = kmer_positions - distinct_kmers
- `shuffle_ratio = observed / mean(null)` if mean(null) >= min_shuffle_mean else `missing`
- `z_score` gated by the same `min_shuffle_mean` floor
- `shuffle_p_upper = (1 + #{null >= observed}) / (n_shuffles + 1)`
- `null_degenerate`: true when dinucleotide nulls fell back / had duplicates
"""
function gqc_repeat_null_stats(sequences::Vector{String};
                               k::Int,
                               n_shuffles::Int=GQC_N_SHUFFLES,
                               shuffle_seed::Int=GQC_SHUFFLE_SEED,
                               salt::Int=0,
                               genus_index::Int=1,
                               circular::Bool=false,
                               canonical::Bool=true,
                               min_shuffle_mean::Float64=GQC_MIN_SHUFFLE_MEAN,
                               null_kind::Symbol=:dinucleotide)
    observed = gqc_repeat_counts(sequences; k=k, circular=circular, canonical=canonical)
    rng = Random.Xoshiro(shuffle_seed + 1000 * salt + genus_index)
    null_repeats = Float64[]
    null_degenerate = false
    if null_kind == :dinucleotide
        shuffles, null_degenerate = gqc_checked_dinucleotide_shuffles(sequences, rng; n_shuffles=n_shuffles)
        for segs in shuffles
            r = gqc_repeat_counts(segs; k=k, circular=circular, canonical=canonical)
            push!(null_repeats, Float64(r.repeats))
        end
    elseif null_kind == :codon
        norms = [replace(uppercase(s), 'U' => 'T') for s in sequences]
        joined_orig = join(norms, '|')
        seen = Set{String}()
        for _ in 1:n_shuffles
            segs = [codon_shuffle(s, rng) for s in norms]
            joined = join(segs, '|')
            if joined == joined_orig || joined in seen
                null_degenerate = true
            end
            push!(seen, joined)
            r = gqc_repeat_counts(segs; k=k, circular=circular, canonical=canonical)
            push!(null_repeats, Float64(r.repeats))
        end
        length(seen) < n_shuffles && (null_degenerate = true)
    else
        error("unknown null_kind=$null_kind")
    end
    μ = Statistics.mean(null_repeats)
    σ = length(null_repeats) > 1 ? Statistics.std(null_repeats) : 0.0
    ratio = μ >= min_shuffle_mean ? observed.repeats / μ : missing
    z = (μ >= min_shuffle_mean && σ > 0) ? (observed.repeats - μ) / σ : missing
    p_upper = (1 + count(>=(observed.repeats), null_repeats)) / (n_shuffles + 1)
    return (;
        observed_repeats=observed.repeats,
        observed_kmer_positions=observed.kmer_positions,
        observed_distinct_kmers=observed.distinct_kmers,
        null_mean=μ,
        null_std=σ,
        shuffle_ratio=ratio,
        z_score=z,
        shuffle_p_upper=p_upper,
        null_repeats=null_repeats,
        null_degenerate=null_degenerate,
    )
end

"""
Generate null shuffles once for `null_kind`, then score repeat stats at each k.

Bit-identical to calling `gqc_repeat_null_stats` once per k with the same seed
(each of those calls reseeds to the same RNG state before shuffling).
"""
function gqc_repeat_null_stats_for_ks(sequences::Vector{String};
                                      ks::AbstractVector{<:Integer},
                                      n_shuffles::Int=GQC_N_SHUFFLES,
                                      shuffle_seed::Int=GQC_SHUFFLE_SEED,
                                      salt::Int=0,
                                      genus_index::Int=1,
                                      circular::Bool=false,
                                      canonical::Bool=true,
                                      min_shuffle_mean::Float64=GQC_MIN_SHUFFLE_MEAN,
                                      null_kind::Symbol=:dinucleotide)
    rng = Random.Xoshiro(shuffle_seed + 1000 * salt + genus_index)
    null_degenerate = false
    shuffle_sets = Vector{Vector{String}}()
    if null_kind == :dinucleotide
        shuffle_sets, null_degenerate = gqc_checked_dinucleotide_shuffles(
            sequences, rng; n_shuffles=n_shuffles)
    elseif null_kind == :codon
        norms = [replace(uppercase(s), 'U' => 'T') for s in sequences]
        joined_orig = join(norms, '|')
        seen = Set{String}()
        for _ in 1:n_shuffles
            segs = [codon_shuffle(s, rng) for s in norms]
            joined = join(segs, '|')
            (joined == joined_orig || joined in seen) && (null_degenerate = true)
            push!(seen, joined)
            push!(shuffle_sets, segs)
        end
        length(seen) < n_shuffles && (null_degenerate = true)
    else
        error("unknown null_kind=$null_kind")
    end
    out = Dict{Int, NamedTuple}()
    for k in ks
        observed = gqc_repeat_counts(sequences; k=Int(k), circular=circular, canonical=canonical)
        null_repeats = Float64[
            Float64(gqc_repeat_counts(segs; k=Int(k), circular=circular, canonical=canonical).repeats)
            for segs in shuffle_sets
        ]
        μ = Statistics.mean(null_repeats)
        σ = length(null_repeats) > 1 ? Statistics.std(null_repeats) : 0.0
        ratio = μ >= min_shuffle_mean ? observed.repeats / μ : missing
        z = (μ >= min_shuffle_mean && σ > 0) ? (observed.repeats - μ) / σ : missing
        p_upper = (1 + count(>=(observed.repeats), null_repeats)) / (n_shuffles + 1)
        out[Int(k)] = (;
            observed_repeats=observed.repeats,
            observed_kmer_positions=observed.kmer_positions,
            observed_distinct_kmers=observed.distinct_kmers,
            null_mean=μ,
            null_std=σ,
            shuffle_ratio=ratio,
            z_score=z,
            shuffle_p_upper=p_upper,
            null_repeats=null_repeats,
            null_degenerate=null_degenerate,
        )
    end
    return out
end

# =============================================================================
# Reference concordance
# =============================================================================

const GQC_REFERENCE_EXT_RE = r"\.(fa|fas|fna|ffn|fasta)(\.gz)?$"i

"""
Load reference FASTA records from `reference_dir`.

Each FASTA record becomes its own sequence (no chimeric joins across records).
Accepts `.fa`, `.fas`, `.fna`, `.ffn`, `.fasta`, optional `.gz`. `U`→`T`.
Bytes before the first `>` are ignored (not prepended to the first record).

Returns a named tuple `(sequences, n_records, files_matched, note)`.
"""
function _gqc_load_reference_sequences(reference_dir::AbstractString)
    seqs = String[]
    files_matched = String[]
    if !isdir(reference_dir)
        return (sequences=seqs, n_records=0, files_matched=files_matched,
                note="no_reference_dir")
    end
    entries = sort(readdir(reference_dir; join=true))
    any_file = false
    files_failed = String[]
    for file in entries
        isfile(file) || continue
        any_file = true
        occursin(GQC_REFERENCE_EXT_RE, file) || continue
        try
            raw_io = open(file)
            io = endswith(lowercase(file), ".gz") ? CodecZlib.GzipDecompressorStream(raw_io) : raw_io
            try
                buf = IOBuffer()
                in_record = false
                function flush_record!()
                    s = replace(String(take!(buf)), 'U' => 'T')
                    !isempty(s) && push!(seqs, s)
                    buf = IOBuffer()
                    return nothing
                end
                for line in eachline(io)
                    if startswith(line, '>')
                        in_record && flush_record!()
                        take!(buf)  ## drop any pre-header junk
                        in_record = true
                    elseif in_record
                        print(buf, uppercase(strip(line)))
                    end
                end
                in_record && flush_record!()
            finally
                close(io)
            end
            push!(files_matched, file)
        catch e
            push!(files_failed, basename(file) * ": " * sprint(showerror, e))
            @warn "gqc reference file failed; skipping" file exception=e
        end
    end
    note = if isempty(files_matched) && !isempty(files_failed)
        "reference_files_failed: " * join(files_failed, "; ")
    elseif isempty(files_matched)
        any_file ? "reference_dir_no_matching_extensions" : "reference_dir_empty"
    elseif isempty(seqs)
        "reference_files_matched_but_zero_records"
    elseif !isempty(files_failed)
        "ok_with_skipped_files"
    else
        "ok"
    end
    @info "gqc reference load" dir=reference_dir n_records=length(seqs) files=length(files_matched) note=note
    return (sequences=seqs, n_records=length(seqs), files_matched=files_matched, note=note)
end

"""
Parse GFA `L` overlap fields like `30M` / `30` into integers (best-effort).
"""
function _gqc_link_overlaps(g::GFAGraph)::Vector{Int}
    out = Int[]
    for link in g.links
        m = match(r"^(\d+)", string(link.overlap))
        m === nothing && continue
        push!(out, parse(Int, m.captures[1]))
    end
    return out
end

"""
True when every segment has the same length `L` and the graph looks like a
k-mer-node (de Bruijn) graph:

- `L >= maximum(ks)` when `ks` is non-empty, or
- all GFA link overlaps equal `L - 1` when overlaps are present, or
- `L == graph_k` when metadata k is set

Does **not** require `metadata["k"]` or that `L` appear in the sweep `ks`.
"""
function _gqc_is_kmer_node_graph(seg_seqs::Vector{String};
                                 graph_k::Int=0,
                                 ks::AbstractVector{<:Integer}=Int[],
                                 overlaps::AbstractVector{<:Integer}=Int[])
    isempty(seg_seqs) && return false
    lens = length.(seg_seqs)
    all(==(first(lens)), lens) || return false
    L = first(lens)
    L > 0 || return false
    if !isempty(ks) && L >= maximum(Int(k) for k in ks)
        return true
    end
    if !isempty(overlaps) && all(==(L - 1), overlaps)
        return true
    end
    (graph_k > 0 && L == graph_k) && return true
    return any(L == Int(k) for k in ks)
end

"""
Resolve sequences for repeat nulls.

Returns `(source, sequences, applicable, reason)`.
"""
function _gqc_resolve_null_sequences(null_source::Symbol,
                                     ref_load,
                                     seg_seqs::Vector{String};
                                     graph_k::Int=0,
                                     ks::AbstractVector{<:Integer}=Int[],
                                     overlaps::AbstractVector{<:Integer}=Int[])
    null_source in (:auto, :reference, :segments) || error(
        "null_source must be :auto, :reference, or :segments; got $null_source"
    )
    ref_seqs = ref_load.sequences
    if null_source === :reference
        isempty(ref_seqs) && error(
            "null_source=:reference but no reference records were loaded " *
            "(note=$(ref_load.note); matched=$(length(ref_load.files_matched)) files)"
        )
        return (:reference, ref_seqs, true, "explicit_reference")
    end
    if null_source === :auto && !isempty(ref_seqs)
        return (:reference, ref_seqs, true, "auto_reference_present")
    end
    ## segments path
    is_kmer = _gqc_is_kmer_node_graph(seg_seqs; graph_k=graph_k, ks=ks, overlaps=overlaps)
    reason = if null_source === :segments
        is_kmer ? "explicit_segments_kmer_node" : "explicit_segments"
    elseif null_source === :auto
        base = isempty(ref_load.files_matched) ? "auto_fallback_no_reference_files" :
               (ref_load.n_records == 0 ? "auto_fallback_zero_records ($(ref_load.note))" :
                "auto_fallback_segments")
        is_kmer ? base * ";kmer_node_guard" : base
    else
        "segments"
    end
    return (:segments, seg_seqs, !is_kmer, reason)
end

function gqc_reference_concordance(g::GFAGraph, reference_seqs::Vector{String}; k::Int, circular::Bool=false)
    isempty(reference_seqs) && return Dict{String, Any}(
        "ref_kmer_fraction" => missing,
        "ref_path_fraction" => missing,
        "ref_absent_kmers" => missing,
        "ref_breaks" => missing,
        "ref_walkability" => missing,
    )
    ## Concordance always folds RC (shared k-mer content). Known limitation vs
    ## strand_mode=:single repeat nulls — documented in docs/gqc-v0.1-spec.md (I6).
    seg_seqs = [replace(uppercase(s.sequence), 'U' => 'T')
                for s in values(g.segments) if s.sequence != "*"]
    graph_kmers = Set{String}()
    for seq in seg_seqs
        n = length(seq)
        n < k && continue
        for i in 1:(n - k + 1)
            raw = seq[i:(i + k - 1)]
            occursin(r"[^ACGT]", raw) && continue
            push!(graph_kmers, _gqc_canonical_kmer(raw))
        end
    end
    ref_kmers = Set{String}()
    ref_positions = 0
    present = 0
    for seq in reference_seqs
        n = length(seq)
        ## Skip records shorter than k (circular wrap would BoundsError when n < k-1).
        n < k && continue
        n_pos = circular ? n : (n - k + 1)
        for i in 1:n_pos
            if circular
                if i + k - 1 <= n
                    raw = seq[i:(i + k - 1)]
                else
                    need = k - (n - i + 1)
                    need <= 0 && continue
                    need > n && continue
                    raw = seq[i:n] * seq[1:need]
                end
            else
                raw = seq[i:(i + k - 1)]
            end
            length(raw) != k && continue
            occursin(r"[^ACGT]", raw) && continue
            ck = _gqc_canonical_kmer(uppercase(raw))
            push!(ref_kmers, ck)
            ref_positions += 1
            ck in graph_kmers && (present += 1)
        end
    end
    absent = length(setdiff(ref_kmers, graph_kmers))
    frac = ref_positions > 0 ? present / ref_positions : missing
    ## path fraction / walkability: consecutive ref k-mers both in the node set
    ## (does not read GFA links; known limitation I6).
    walkable = 0
    pair_total = 0
    for seq in reference_seqs
        n = length(seq)
        n < k + 1 && continue
        for i in 1:(n - k)
            a = uppercase(seq[i:(i + k - 1)])
            b = uppercase(seq[(i + 1):(i + k)])
            (occursin(r"[^ACGT]", a) || occursin(r"[^ACGT]", b)) && continue
            ca = _gqc_canonical_kmer(a)
            cb = _gqc_canonical_kmer(b)
            pair_total += 1
            if ca in graph_kmers && cb in graph_kmers
                walkable += 1
            end
        end
    end
    path_frac = pair_total > 0 ? walkable / pair_total : missing
    breaks = pair_total > 0 ? pair_total - walkable : missing
    return Dict{String, Any}(
        "ref_kmer_fraction" => frac,
        "ref_path_fraction" => path_frac,
        "ref_absent_kmers" => absent,
        "ref_breaks" => breaks,
        "ref_walkability" => path_frac,
    )
end

# =============================================================================
# Standard-format writers
# =============================================================================

"""
Resolve `graph_id` for node TSV / formats writers.

Prefer an explicit non-empty `graph_id` (what `run_gqc` passes from metadata /
`basename(input_dir)`). Otherwise fall back to `basename(dirname(g.path))`,
matching the no-metadata `run_gqc` default when the GFA lives at
`<input_dir>/graph.gfa`.
"""
function _gqc_resolve_node_tsv_graph_id(g::GFAGraph;
                                       graph_id::AbstractString="")::String
    gid = strip(string(graph_id))
    !isempty(gid) && return gid
    p = string(g.path)
    isempty(p) && return ""
    parent = basename(dirname(abspath(p)))
    return (isempty(parent) || parent == ".") ? "" : parent
end

"""
Write the unified node TSV (`GQC_NODE_TSV_COLUMNS`) for one graph.

Used by both `nodes/<slug>.tsv` and `formats/nodes.tsv`. Schema/column naming
only — values are the same quantities previously split across the two files.
`graph_id` is resolved via `_gqc_resolve_node_tsv_graph_id`.
"""
function _gqc_write_nodes_tsv(io::IO,
                              g::GFAGraph,
                              ids::Vector{String},
                              adj::Vector{Vector{Int}};
                              graph_id::AbstractString="")
    gid = _gqc_resolve_node_tsv_graph_id(g; graph_id=graph_id)
    degrees = _gqc_degrees(adj)
    indeg, outdeg, _, _ = _gqc_directed_degrees(g)
    println(io, join(GQC_NODE_TSV_COLUMNS, '\t'))
    for (i, id) in enumerate(ids)
        s = g.segments[id]
        seq = s.sequence == "*" ? "" : s.sequence
        d = s.depth === nothing ? "" : string(s.depth)
        println(io, gid, '\t', id, '\t', i, '\t', s.length, '\t', d, '\t',
                degrees[i], '\t', indeg[id], '\t', outdeg[id], '\t', seq)
    end
    return nothing
end

function gqc_write_formats(g::GFAGraph, formats_dir::AbstractString;
                           write_node_link_json::Bool=true,
                           graph_id::AbstractString="")
    mkpath(formats_dir)
    ids, adj, _ = _gqc_undirected_adjacency(g)
    n = length(ids)
    ## nodes.tsv (unified schema; same columns as nodes/<slug>.tsv)
    nodes_path = joinpath(formats_dir, "nodes.tsv")
    open(nodes_path, "w") do io
        _gqc_write_nodes_tsv(io, g, ids, adj; graph_id=graph_id)
    end
    ## edges.tsv
    edges_path = joinpath(formats_dir, "edges.tsv")
    open(edges_path, "w") do io
        println(io, "from\tfrom_orient\tto\tto_orient\toverlap")
        for link in g.links
            println(io, link.from, '\t', link.from_forward ? "+" : "-", '\t',
                    link.to, '\t', link.to_forward ? "+" : "-", '\t', link.overlap)
        end
    end
    ## Matrix Market coordinate pattern *symmetric* (undirected simple graph).
    ## Store lower triangle (row >= col), including diagonal self-loops.
    ## Directed multi-edges remain in edges.tsv / GraphML / node_link.json.
    index = Dict(id => i for (i, id) in enumerate(ids))
    edges_unique = Set{Tuple{Int,Int}}()
    for link in g.links
        haskey(index, link.from) || continue
        haskey(index, link.to) || continue
        u = index[link.from]; v = index[link.to]
        a, b = minmax(u, v)  ## a <= b
        push!(edges_unique, (b, a))  ## lower triangle: row >= col
    end
    mtx_path = joinpath(formats_dir, "adjacency.mtx")
    open(mtx_path, "w") do io
        println(io, "%%MatrixMarket matrix coordinate pattern symmetric")
        println(io, "% undirected simple graph; lower triangle (row>=col); 1-based = nodes.tsv index")
        println(io, "% directed GFA L lines are in edges.tsv / graph.graphml / node_link.json")
        println(io, n, ' ', n, ' ', length(edges_unique))
        for (row, col) in sort!(collect(edges_unique))
            println(io, row, ' ', col)
        end
    end
    ## GraphML
    graphml_path = joinpath(formats_dir, "graph.graphml")
    open(graphml_path, "w") do io
        println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
        println(io, """<graphml xmlns="http://graphml.graphdrawing.org/xmlns">""")
        println(io, """  <key id="d_len" for="node" attr.name="length" attr.type="int"/>""")
        println(io, """  <key id="d_depth" for="node" attr.name="depth" attr.type="double"/>""")
        println(io, """  <key id="d_seq" for="node" attr.name="sequence" attr.type="string"/>""")
        println(io, """  <key id="e_ov" for="edge" attr.name="overlap" attr.type="string"/>""")
        println(io, """  <key id="e_fo" for="edge" attr.name="from_orient" attr.type="string"/>""")
        println(io, """  <key id="e_to" for="edge" attr.name="to_orient" attr.type="string"/>""")
        println(io, """  <graph id="G" edgedefault="directed">""")
        for id in ids
            s = g.segments[id]
            println(io, """    <node id="$(id)">""")
            println(io, """      <data key="d_len">$(s.length)</data>""")
            if s.depth !== nothing
                println(io, """      <data key="d_depth">$(s.depth)</data>""")
            end
            if s.sequence != "*"
                println(io, """      <data key="d_seq">$(s.sequence)</data>""")
            end
            println(io, """    </node>""")
        end
        for (ei, link) in enumerate(g.links)
            println(io, """    <edge id="e$ei" source="$(link.from)" target="$(link.to)">""")
            println(io, """      <data key="e_ov">$(link.overlap)</data>""")
            println(io, """      <data key="e_fo">$(link.from_forward ? "+" : "-")</data>""")
            println(io, """      <data key="e_to">$(link.to_forward ? "+" : "-")</data>""")
            println(io, """    </edge>""")
        end
        println(io, """  </graph>""")
        println(io, """</graphml>""")
    end
    written = [nodes_path, edges_path, mtx_path, graphml_path]
    if write_node_link_json
        json_path = joinpath(formats_dir, "node_link.json")
        nodes = [Dict("id" => id, "index" => i, "length" => g.segments[id].length,
                      "depth" => g.segments[id].depth) for (i, id) in enumerate(ids)]
        links_json = [Dict("source" => l.from, "target" => l.to,
                           "from_orient" => l.from_forward ? "+" : "-",
                           "to_orient" => l.to_forward ? "+" : "-",
                           "overlap" => l.overlap) for l in g.links]
        open(json_path, "w") do io
            JSON.print(io, Dict("nodes" => nodes, "links" => links_json), 2)
        end
        push!(written, json_path)
    end
    return written
end

# =============================================================================
# Plots
# =============================================================================

function _gqc_save_plot(p, plots_dir, metric, graph_id, k)
    mkpath(joinpath(plots_dir, metric))
    base = joinpath(plots_dir, metric, "$(graph_id)__k$(k)")
    try
        Plots.savefig(p, base * ".png")
        Plots.savefig(p, base * ".svg")
    catch e
        @warn "gqc plot save failed" metric exception=e
    end
    return base
end

function gqc_make_plots(summary_row::Dict, degree_hist, unitig_lengths, coverages,
                        null_ratio_ci, plots_dir::AbstractString, graph_id::String, k::Int;
                        adj=nothing, ids=nothing, g::Union{Nothing,GFAGraph}=nothing)
    ## Safe on HPC headless: callers should set ENV["GKSwstype"]="100" (see PR docs).
    ## Catch plot *creation* errors as well as savefig failures.
    try
        if !isempty(degree_hist)
            degs = collect(keys(degree_hist))
            counts = [degree_hist[d] for d in degs]
            ## Isolates have degree 0 — log10(0) aborts; use linear axes then.
            if any(<=(0), degs) || any(<=(0), counts)
                p = Plots.plot(degs, counts; seriestype=:scatter,
                               xlabel="degree", ylabel="count", title="Degree distribution", legend=false)
            else
                p = Plots.plot(degs, counts; seriestype=:scatter, xscale=:log10, yscale=:log10,
                               xlabel="degree", ylabel="count", title="Degree distribution", legend=false)
            end
            _gqc_save_plot(p, plots_dir, "degree", graph_id, k)
        end
        if !isempty(unitig_lengths)
            sorted = sort(unitig_lengths; rev=true)
            cum = cumsum(sorted)
            p = Plots.plot(1:length(cum), cum; xlabel="unitig rank", ylabel="cumulative bp",
                           title="Cumulative unitig length", legend=false)
            Plots.hline!(p, [sum(sorted) / 2]; label="half total", linestyle=:dash)
            _gqc_save_plot(p, plots_dir, "unitig_n50", graph_id, k)
        end
        if !isempty(coverages)
            p = Plots.histogram(coverages; xlabel="depth", ylabel="count", title="Coverage", legend=false)
            _gqc_save_plot(p, plots_dir, "coverage", graph_id, k)
        end
        bars = ["tips", "bubbles", "branch", "self_loops", "tangles"]
        vals = [summary_row[k_] for k_ in ("n_tips", "n_bubbles", "n_branch", "n_self_loops", "n_tangles")]
        p = Plots.bar(bars, vals; legend=false, title="Structural counts", xrotation=30)
        _gqc_save_plot(p, plots_dir, "structure", graph_id, k)
        if null_ratio_ci !== nothing
            ratio, lo, hi = null_ratio_ci
            if ratio isa Number && lo isa Number && hi isa Number
                p = Plots.scatter([1], [ratio]; yerror=[(ratio - lo, hi - ratio)],
                                  xlabel="", ylabel="repeat ratio", title="Repeat vs null",
                                  xticks=([], []), legend=false, xlim=(0.5, 1.5))
                Plots.hline!(p, [1.0]; linestyle=:dash)
                _gqc_save_plot(p, plots_dir, "repeat_null", graph_id, k)
            end
        end
        if g !== nothing && length(g.segments) <= 5000 && adj !== nothing
            comps = _gqc_connected_components(adj)
            sizes = sort([length(c) for c in comps]; rev=true)
            p = Plots.plot(1:length(sizes), sizes; xlabel="component rank", ylabel="nodes",
                           title="Component sizes", legend=false)
            _gqc_save_plot(p, plots_dir, "components", graph_id, k)
        end
    catch e
        @warn "gqc_make_plots failed" graph_id k exception=e
    end
    return nothing
end

# =============================================================================
# Packet driver
# =============================================================================

function _gqc_read_metadata(path::AbstractString)
    isfile(path) || return Dict{String, Any}()
    return JSON.parsefile(path)
end

"""
Canonicalize params for hashing / manifest storage.

Integer vectors (`ks`, `bh_ks`) are sorted **numerically**. Other vectors are
sorted by `string`. Hash is permutation-invariant; stored order matches analysis.
"""
function _gqc_canonicalize_params(params::Dict)
    out = Dict{String, Any}()
    for (k, v) in params
        key = string(k)
        if v isa AbstractVector
            vals = collect(v)
            if !isempty(vals) && all(x -> x isa Integer, vals)
                out[key] = sort(Int[Int(x) for x in vals])
            else
                out[key] = sort(vals; by=string)
            end
        else
            out[key] = v
        end
    end
    return out
end

function _gqc_params_hash(params::Dict)
    canonical = _gqc_canonicalize_params(params)
    payload = JSON.json(sort(collect(canonical); by=first))
    return SHA.bytes2hex(SHA.sha256(payload))
end

"""
Slug `graph_id` for filesystem paths. Refuses absolute paths and `..` traversal.
"""
function _gqc_slug_graph_id(graph_id::AbstractString)::String
    s = string(graph_id)
    isempty(strip(s)) && error("graph_id is empty")
    if startswith(s, '/') || startswith(s, '\\') || occursin(r"^[A-Za-z]:[/\\]", s) ||
       occursin("..", s) || occursin('\0', s)
        error("refusing absolute or traversal graph_id: $(repr(s))")
    end
    slug = replace(s, r"[^\w.\-]+" => "_")
    slug = strip(slug, ['_', '.'])
    (isempty(slug) || slug in (".", "..")) && error("graph_id slugs to empty/invalid: $(repr(s))")
    return slug
end

"""
Resolve strand_mode / canonical (I7).

- Keyword `strand_mode === nothing` → `metadata["strand_mode"]` → default `:single`
- Keyword `canonical === nothing` → derived from strand mode (no silent override)
- `:single` → `canonical=false` unless caller passed `canonical=true` explicitly
- `:canonical` → `canonical=true`
- `:double` → honor explicit `canonical` (default `true`)
"""
function _gqc_resolve_strand(metadata::Dict, strand_mode_kw, canonical_kw)
    resolved_sm = if strand_mode_kw !== nothing
        Symbol(strand_mode_kw)
    elseif haskey(metadata, "strand_mode")
        Symbol(string(metadata["strand_mode"]))
    else
        :single
    end
    resolved_sm in (:single, :double, :canonical) || error(
        "unknown strand_mode=$resolved_sm; expected :single, :double, or :canonical"
    )
    eff = if resolved_sm === :canonical
        true
    elseif resolved_sm === :single
        canonical_kw === nothing ? false : Bool(canonical_kw)
    else
        canonical_kw === nothing ? true : Bool(canonical_kw)
    end
    return (strand_mode=resolved_sm, canonical=eff)
end

"""
Effective k-mer RC folding from resolved `strand_mode` + optional canonical flag.
"""
function _gqc_effective_canonical(strand_mode::Symbol, canonical::Bool)::Bool
    return _gqc_resolve_strand(Dict{String, Any}(), strand_mode, canonical).canonical
end

"""
Packet directory under `output_root` using slugged `graph_id`.

`overwrite`:
- `:error` — refuse non-empty existing dir (default)
- `:skip` — return existing dir if it has `manifest.json` + `CHECKSUMS.sha256`
- `:replace` — remove and recreate
"""
function _gqc_packet_dir(output_root::AbstractString, graph_id::AbstractString;
                         overwrite::Symbol=:error)
    overwrite in (:error, :skip, :replace) || error(
        "overwrite must be :error, :skip, or :replace; got $overwrite"
    )
    slug = _gqc_slug_graph_id(graph_id)
    out_dir = joinpath(abspath(output_root), slug)
    if isdir(out_dir)
        entries = filter(name -> !startswith(name, '.'), readdir(out_dir))
        if !isempty(entries)
            complete = isfile(joinpath(out_dir, "manifest.json")) &&
                       isfile(joinpath(out_dir, "CHECKSUMS.sha256"))
            if overwrite === :skip && complete
                return (dir=out_dir, slug=slug, skipped=true)
            elseif overwrite === :replace
                rm(out_dir; recursive=true, force=true)
            else
                error(
                    "refusing to overwrite existing gqc packet directory: $out_dir " *
                    "(graph_id=$(repr(graph_id)), slug=$slug). " *
                    "Pass overwrite=:skip or overwrite=:replace to resume/replace."
                )
            end
        end
    end
    return (dir=out_dir, slug=slug, skipped=false)
end

function _gqc_git_commit()
    try
        pkg_root = dirname(dirname(@__FILE__))
        return strip(read(`git -C $(pkg_root) rev-parse HEAD`, String))
    catch
        return "unknown"
    end
end

"""
    resolve_gqc_input_dirs(input_dir) -> NamedTuple{(:mode, :dirs), Tuple{Symbol, Vector{String}}}

Classify `input_dir` for gqc:

- **`:single_graph`**: `input_dir/graph.gfa` exists. Returns that one folder.
- **`:graph_set`**: no top-level `graph.gfa`, but one or more *immediate* child
  directories each contain `graph.gfa`. Returns those children (sorted by name).

A top-level `graph.gfa` always wins (single-graph mode), even if child folders
also look like graphs. Nested graph-sets deeper than one level are ignored.

# Returns
`(mode = :single_graph | :graph_set, dirs = Vector{String})` with absolute paths.
"""
function resolve_gqc_input_dirs(input_dir::AbstractString)
    root = abspath(input_dir)
    isdir(root) || error("gqc input_dir is not a directory: $root")
    if isfile(joinpath(root, "graph.gfa"))
        return (mode = :single_graph, dirs = String[root])
    end
    children = String[]
    for name in sort(readdir(root))
        child = joinpath(root, name)
        isdir(child) || continue
        if isfile(joinpath(child, "graph.gfa"))
            push!(children, child)
        end
    end
    isempty(children) && error(
        "gqc input_dir must be a single-graph folder (contains graph.gfa) or a " *
        "graph-set folder (immediate children each contain graph.gfa); got: $root"
    )
    return (mode = :graph_set, dirs = children)
end

"""
    run_gqc(input_dir, output_root; kwargs...) -> String

Run gqc v0.1 and write QC packet(s) under `output_root/<graph_set>/`.

# Arguments
- `input_dir`: either
  - **single-graph folder** — contains `graph.gfa`, optional `metadata.json`,
    optional `reference/*.fasta`; or
  - **graph-set folder** — no top-level `graph.gfa`; each immediate child that
    contains `graph.gfa` is processed as its own graph (see
    `resolve_gqc_input_dirs`)
- `output_root`: parent directory for packet output

# Keyword arguments
- `ks`: k values for repeat-null / reference concordance sweep (default: `[metadata["k"]]`)
- `n_shuffles`: dinucleotide/codon null replicates (default: `100`)
- `shuffle_seed`: base RNG seed (default: `20261021`; per-graph RNG =
  `Xoshiro(shuffle_seed + 1000*salt + genus_index)`)
- `salt`, `genus_index`: seeding offsets (paper formula). In graph-set mode,
  graph `i` (1-based, sorted child order) uses `genus_index + i - 1`.
- `min_shuffle_mean`: gates both `shuffle_ratio` and `z_score` (default: `1.0`)
- `bh_ks`: k values included in **within-graph** BH on dinucleotide `p_upper`
  (default: `[11]`). Packet `dinucleotide_bh_q_within_graph` is **not** the
  paper's cross-genus FDR — use `Mycelia.gqc_cross_bh` (or equivalent) outside
  `run_gqc`.
- `circular`: k-mer circularity (default: `metadata["circular"]` or `false`)
- `canonical`: `nothing` (default) → derived from `strand_mode`; explicit `true`/`false`
  is honored (including under `:single`)
- `strand_mode`: `nothing` (default) → `metadata["strand_mode"]` → `:single`;
  or `:single` / `:canonical` / `:double`
- `write_node_link_json`: include `formats/node_link.json` (default: `true`)
- `make_plots`: write plots under `plots/` (default: `true`; set
  `ENV["GKSwstype"]="100"` on HPC)
- `null_kinds`: `(:dinucleotide, :codon)` by default
- `null_source`: `:auto` (default; use reference records when present),
  `:reference`, or `:segments`. Uniform-length k-mer-node graphs with a segment
  source mark `null_applicable=false`.
- `overwrite`: `:error` (default), `:skip` (resume complete packets), or `:replace`

# Returns
- Single-graph mode: absolute path to the packet directory `output_root/<slug>/`.
- Graph-set mode: absolute path to `output_root` (each child under unique slug).
"""
function run_gqc(input_dir::AbstractString,
                 output_root::AbstractString;
                 ks::AbstractVector{<:Integer}=Int[],
                 n_shuffles::Int=GQC_N_SHUFFLES,
                 shuffle_seed::Int=GQC_SHUFFLE_SEED,
                 salt::Int=0,
                 genus_index::Int=1,
                 min_shuffle_mean::Float64=GQC_MIN_SHUFFLE_MEAN,
                 bh_ks::Union{Nothing, AbstractVector{<:Integer}}=nothing,
                 circular::Union{Nothing, Bool}=nothing,
                 canonical::Union{Nothing, Bool}=nothing,
                 strand_mode::Union{Nothing, Symbol}=nothing,
                 write_node_link_json::Bool=true,
                 make_plots::Bool=true,
                 null_kinds::Tuple=(:dinucleotide, :codon),
                 null_source::Symbol=:auto,
                 overwrite::Symbol=:error,
                 gqc_version::String=GQC_VERSION)::String

    resolved = resolve_gqc_input_dirs(input_dir)
    output_root = abspath(output_root)
    mkpath(output_root)
    packets = String[]
    errors = Dict{String, Any}[]
    for (i, graph_dir) in enumerate(resolved.dirs)
        try
            packet = _run_gqc_one(graph_dir, output_root;
                ks=ks,
                n_shuffles=n_shuffles,
                shuffle_seed=shuffle_seed,
                salt=salt,
                genus_index=genus_index + i - 1,
                min_shuffle_mean=min_shuffle_mean,
                bh_ks=bh_ks,
                circular=circular,
                canonical=canonical,
                strand_mode=strand_mode,
                write_node_link_json=write_node_link_json,
                make_plots=make_plots,
                null_kinds=null_kinds,
                null_source=null_source,
                overwrite=overwrite,
                gqc_version=gqc_version,
            )
            push!(packets, packet)
        catch e
            resolved.mode === :single_graph && rethrow(e)
            rec = Dict{String, Any}(
                "graph_dir" => graph_dir,
                "error" => sprint(showerror, e),
                "exception_type" => string(typeof(e)),
            )
            push!(errors, rec)
            @error "gqc graph-set child failed; continuing" graph_dir exception=e
        end
    end
    if !isempty(errors)
        err_path = joinpath(output_root, "gqc_graph_set_errors.jsonl")
        open(err_path, "w") do io
            for rec in errors
                println(io, JSON.json(rec))
            end
        end
    end
    if resolved.mode === :single_graph
        return only(packets)
    end
    return output_root
end

function _run_gqc_one(input_dir::AbstractString,
                      output_root::AbstractString;
                      ks::AbstractVector{<:Integer}=Int[],
                      n_shuffles::Int=GQC_N_SHUFFLES,
                      shuffle_seed::Int=GQC_SHUFFLE_SEED,
                      salt::Int=0,
                      genus_index::Int=1,
                      min_shuffle_mean::Float64=GQC_MIN_SHUFFLE_MEAN,
                      bh_ks::Union{Nothing, AbstractVector{<:Integer}}=nothing,
                      circular::Union{Nothing, Bool}=nothing,
                      canonical::Union{Nothing, Bool}=nothing,
                      strand_mode::Union{Nothing, Symbol}=nothing,
                      write_node_link_json::Bool=true,
                      make_plots::Bool=true,
                      null_kinds::Tuple=(:dinucleotide, :codon),
                      null_source::Symbol=:auto,
                      overwrite::Symbol=:error,
                      gqc_version::String=GQC_VERSION)::String

    input_dir = abspath(input_dir)
    output_root = abspath(output_root)
    gfa_path = joinpath(input_dir, "graph.gfa")
    isfile(gfa_path) || error("missing graph.gfa in $input_dir")
    metadata = _gqc_read_metadata(joinpath(input_dir, "metadata.json"))
    graph_id = string(get(metadata, "graph_id", basename(input_dir)))
    graph_set = string(get(metadata, "label", graph_id))
    k_meta = Int(get(metadata, "k", 0))
    circ = circular === nothing ? Bool(get(metadata, "circular", false)) : circular
    strand_res = _gqc_resolve_strand(metadata, strand_mode, canonical)
    resolved_strand = strand_res.strand_mode
    effective_canonical = strand_res.canonical
    k_list = isempty(ks) ? (k_meta > 0 ? [k_meta] : Int[]) : sort(unique(Int(k) for k in ks))
    isempty(k_list) && error("no k provided: set metadata[\"k\"] or pass ks=")
    ## Within-graph BH family defaults to k=11 only (C5). Cross-genus FDR is external.
    bh_list = sort(unique(Int(k) for k in (bh_ks === nothing ? GQC_DEFAULT_BH_KS : bh_ks)))
    bh_set = Set(bh_list)

    g = read_gfa_records(gfa_path)
    ids, adj, self_loops = _gqc_undirected_adjacency(g)
    comps = _gqc_connected_components(adj)
    degrees = _gqc_degrees(adj)
    roles = _gqc_node_roles(adj)
    indeg, outdeg, sources, sinks = _gqc_directed_degrees(g)
    bubbles = _gqc_simple_bubbles(adj)
    tangles = _gqc_tangle_count(adj)
    seg_seqs = [replace(uppercase(s.sequence), 'U' => 'T')
                for s in values(g.segments) if s.sequence != "*" && !isempty(s.sequence)]
    default_ov = _gqc_default_overlap(seg_seqs)
    overlap_lookup = _gqc_overlap_lookup(g; default_overlap=default_ov)
    unitig_lengths = _gqc_unitig_lengths(g, ids, adj;
                                         default_overlap=default_ov,
                                         overlap_lookup=overlap_lookup)
    n50stats = _gqc_n50_l50(unitig_lengths)
    cov_q = _gqc_coverage_quantiles(g)
    coverages = Float64[s.depth for s in values(g.segments) if s.depth !== nothing]
    n_seg = length(g.segments)
    n_link = length(g.links)
    total_segment_bp = sum(s.length for s in values(g.segments); init=0)
    ## Overlap-corrected assembly span (sum of unitig spans).
    total_bp = isempty(unitig_lengths) ? total_segment_bp : sum(unitig_lengths)
    density = n_seg > 1 ? (2 * n_link) / (n_seg * (n_seg - 1)) : 0.0
    largest_frac = n_seg > 0 ? maximum(length(c) for c in comps) / n_seg : 0.0
    ## Betti numbers on undirected simple graph: β0 = components, β1 = e - v + c (cyclomatic)
    n_undirected_edges = sum(length(a) for a in adj) ÷ 2
    beta0 = length(comps)
    beta1 = n_undirected_edges - n_seg + beta0

    ref_dir = joinpath(input_dir, "reference")
    ref_load = _gqc_load_reference_sequences(ref_dir)
    ref_seqs = ref_load.sequences
    overlaps = _gqc_link_overlaps(g)
    resolved_null_source, null_seqs, null_applicable, null_reason =
        _gqc_resolve_null_sequences(null_source, ref_load, seg_seqs;
                                    graph_k=k_meta, ks=k_list, overlaps=overlaps)

    packet_info = _gqc_packet_dir(output_root, graph_id; overwrite=overwrite)
    packet_info.skipped && return packet_info.dir
    out_dir = packet_info.dir
    slug = packet_info.slug
    mkpath(out_dir)
    for sub in ("nodes", "edges", "unitigs", "nulls", "plots", "formats")
        mkpath(joinpath(out_dir, sub))
    end

    ## long tables — always use slug in paths (I4); node TSV uses unified schema
    open(joinpath(out_dir, "nodes", "$(slug).tsv"), "w") do io
        _gqc_write_nodes_tsv(io, g, ids, adj; graph_id=graph_id)
    end
    open(joinpath(out_dir, "edges", "$(slug).tsv"), "w") do io
        println(io, "graph_id\tfrom\tfrom_orient\tto\tto_orient\toverlap")
        for link in g.links
            println(io, graph_id, '\t', link.from, '\t', link.from_forward ? "+" : "-", '\t',
                    link.to, '\t', link.to_forward ? "+" : "-", '\t', link.overlap)
        end
    end
    open(joinpath(out_dir, "unitigs", "$(slug).tsv"), "w") do io
        println(io, "graph_id\tunitig_index\tlength")
        for (i, L) in enumerate(sort(unitig_lengths; rev=true))
            println(io, graph_id, '\t', i, '\t', L)
        end
    end

    format_files = gqc_write_formats(g, joinpath(out_dir, "formats");
                                     write_node_link_json=write_node_link_json,
                                     graph_id=graph_id)

    ## Generate shuffles once per null kind, then score each k (I10).
    null_by_kind = Dict{Symbol, Dict{Int, NamedTuple}}()
    if null_applicable
        for kind in null_kinds
            null_by_kind[kind] = gqc_repeat_null_stats_for_ks(null_seqs;
                ks=k_list,
                n_shuffles=n_shuffles,
                shuffle_seed=shuffle_seed,
                salt=salt,
                genus_index=genus_index,
                circular=circ,
                canonical=effective_canonical,
                min_shuffle_mean=min_shuffle_mean,
                null_kind=kind,
            )
        end
    end

    ## per-k summary rows + nulls
    summary_rows = Dict{String, Any}[]
    null_p_dinuc = Float64[]
    null_k_order = Int[]

    for k in k_list
        row = Dict{String, Any}(
            "graph_id" => graph_id,
            "k" => k,
            "n_segments" => n_seg,
            "n_links" => n_link,
            "total_bp" => total_bp,
            "total_segment_bp" => total_segment_bp,
            "density" => density,
            "n_components" => beta0,
            "largest_component_fraction" => largest_frac,
            "beta0" => beta0,
            "beta1" => beta1,
            "n_tips" => length(roles.tips),
            "n_linear" => length(roles.linear),
            "n_branch" => length(roles.branch),
            "n_isolated" => length(roles.isolated),
            "n_sources" => length(sources),
            "n_sinks" => length(sinks),
            "n_unitigs" => length(unitig_lengths),
            "unitig_n50" => n50stats.n50,
            "unitig_l50" => n50stats.l50,
            "unitig_max" => n50stats.max_len,
            "n_bubbles" => bubbles,
            "n_tangles" => tangles,
            "n_self_loops" => self_loops,
            "strand_mode" => string(resolved_strand),
            "canonical" => effective_canonical,
            "gqc_version" => gqc_version,
            ## Always emit null / BH columns (frozen schema); fill below.
            "dinucleotide_observed_repeats" => missing,
            "dinucleotide_null_mean" => missing,
            "dinucleotide_null_std" => missing,
            "dinucleotide_shuffle_ratio" => missing,
            "dinucleotide_z_score" => missing,
            "dinucleotide_shuffle_p_upper" => missing,
            "dinucleotide_null_degenerate" => false,
            "dinucleotide_bh_q_within_graph" => missing,
            "codon_observed_repeats" => missing,
            "codon_null_mean" => missing,
            "codon_null_std" => missing,
            "codon_shuffle_ratio" => missing,
            "codon_z_score" => missing,
            "codon_shuffle_p_upper" => missing,
            "codon_null_degenerate" => false,
            "null_source" => string(resolved_null_source),
            "null_applicable" => null_applicable,
        )
        for cov_key in ("cov_q0", "cov_q25", "cov_q50", "cov_q75", "cov_q100")
            row[cov_key] = get(cov_q, cov_key, missing)
        end
        ref_stats = gqc_reference_concordance(g, ref_seqs; k=k, circular=circ)
        for rk in ("ref_kmer_fraction", "ref_path_fraction", "ref_absent_kmers", "ref_breaks", "ref_walkability")
            row[rk] = get(ref_stats, rk, missing)
        end

        null_ratio_ci = nothing
        for kind in null_kinds
            prefix = string(kind)
            if !null_applicable
                ## K-mer-node graph + segment source: do not emit degenerate p.
                open(joinpath(out_dir, "nulls", "$(slug)__k$(k)__$(prefix).tsv"), "w") do io
                    println(io, "graph_id\tk\tnull_kind\treplicate\tnull_repeats")
                    println(io, "# null_applicable=false reason=$(null_reason)")
                end
                continue
            end
            stats = null_by_kind[kind][k]
            row["$(prefix)_observed_repeats"] = stats.observed_repeats
            row["$(prefix)_null_mean"] = stats.null_mean
            row["$(prefix)_null_std"] = stats.null_std
            row["$(prefix)_shuffle_ratio"] = stats.shuffle_ratio
            row["$(prefix)_z_score"] = stats.z_score
            row["$(prefix)_shuffle_p_upper"] = stats.shuffle_p_upper
            row["$(prefix)_null_degenerate"] = stats.null_degenerate
            ## null long table
            open(joinpath(out_dir, "nulls", "$(slug)__k$(k)__$(prefix).tsv"), "w") do io
                println(io, "graph_id\tk\tnull_kind\treplicate\tnull_repeats")
                for (i, r) in enumerate(stats.null_repeats)
                    println(io, graph_id, '\t', k, '\t', prefix, '\t', i, '\t', r)
                end
            end
            if kind == :dinucleotide
                if k in bh_set
                    push!(null_p_dinuc, stats.shuffle_p_upper)
                    push!(null_k_order, k)
                end
                if stats.shuffle_ratio !== missing && stats.null_std > 0
                    ## rough CI on ratio via delta on null mean ± 1.96 σ/√n
                    se = stats.null_std / sqrt(n_shuffles)
                    lo = stats.null_mean - 1.96 * se
                    hi = stats.null_mean + 1.96 * se
                    r = stats.observed_repeats
                    null_ratio_ci = (
                        stats.shuffle_ratio,
                        hi > 0 ? r / hi : missing,
                        lo > 0 ? r / lo : missing,
                    )
                end
            end
        end

        ## hooks for panlingual consumers (JSON only; not in summary.csv)
        row["graph_feature_vector"] = [
            Float64(n_seg), Float64(n_link), Float64(density), Float64(beta0), Float64(beta1),
            Float64(length(roles.tips)), Float64(bubbles), Float64(get(row, "dinucleotide_shuffle_ratio", 0) isa Number ? row["dinucleotide_shuffle_ratio"] : NaN),
        ]
        row["distance_to_null"] = get(row, "dinucleotide_z_score", missing)

        push!(summary_rows, row)

        if make_plots
            degree_hist = Dict{Int, Int}()
            for d in degrees
                degree_hist[d] = get(degree_hist, d, 0) + 1
            end
            gqc_make_plots(row, degree_hist, unitig_lengths, coverages, null_ratio_ci,
                           joinpath(out_dir, "plots"), slug, k; adj=adj, ids=ids, g=g)
        end
    end

    ## Within-graph BH on dinucleotide p_upper for bh_ks only (not cross-genus FDR).
    if !isempty(null_p_dinuc)
        qvals = benjamini_hochberg(null_p_dinuc)
        qmap = Dict(null_k_order[i] => qvals[i] for i in eachindex(null_k_order))
        for row in summary_rows
            row["dinucleotide_bh_q_within_graph"] = get(qmap, row["k"], missing)
        end
    end

    if make_plots && length(summary_rows) >= 1
        try
            ks_plot = [r["k"] for r in summary_rows]
            b0 = [r["beta0"] for r in summary_rows]
            b1 = [r["beta1"] for r in summary_rows]
            p = Plots.plot(ks_plot, b0; label="β0", xlabel="k", ylabel="Betti", title="Betti vs k")
            Plots.plot!(p, ks_plot, b1; label="β1")
            _gqc_save_plot(p, joinpath(out_dir, "plots"), "betti_vs_k", slug, first(ks_plot))
        catch e
            @warn "gqc betti_vs_k plot failed" exception=e
        end
    end

    ## summary.csv — frozen ordered columns (I2)
    csv_keys = collect(String, GQC_SUMMARY_COLUMNS)
    open(joinpath(out_dir, "summary.csv"), "w") do io
        println(io, join(csv_keys, ','))
        for row in summary_rows
            rendered = String[]
            for col in csv_keys
                v = get(row, col, missing)
                if v === missing
                    push!(rendered, "")
                elseif v isa AbstractString
                    push!(rendered, '"' * replace(v, '"' => "\"\"") * '"')
                elseif v isa Bool
                    push!(rendered, string(v))
                else
                    push!(rendered, string(v))
                end
            end
            println(io, join(rendered, ','))
        end
    end
    open(joinpath(out_dir, "summary.json"), "w") do io
        JSON.print(io, summary_rows, 2)
    end

    ## report.md
    open(joinpath(out_dir, "report.md"), "w") do io
        println(io, "# gqc v$(gqc_version) report — $(graph_id)")
        println(io)
        println(io, "- segments: **$(n_seg)**; links: **$(n_link)**; total_bp (overlap-corrected): **$(total_bp)**; total_segment_bp (raw): **$(total_segment_bp)**")
        println(io, "- components (β0): **$(beta0)**; β1: **$(beta1)**; largest fraction: **$(round(largest_frac; digits=3))**")
        println(io, "- tips: **$(length(roles.tips))**; bubbles: **$(bubbles)**; tangles: **$(tangles)**")
        println(io, "- unitigs: **$(length(unitig_lengths))**; N50: **$(n50stats.n50)**; L50: **$(n50stats.l50)** (overlap-corrected spans)")
        println(io, "- strand_mode: `$(resolved_strand)` → effective canonical=$(effective_canonical)")
        println(io, "- null_source: `$(resolved_null_source)` (requested `$(null_source)`); applicable=$(null_applicable); reason=$(null_reason)")
        println(io, "- reference records loaded: **$(ref_load.n_records)** from $(length(ref_load.files_matched)) file(s)")
        println(io, "- `dinucleotide_bh_q_within_graph` uses bh_ks=[", join(bh_list, ", "),
                "] **inside this graph only** — cross-genus BH is the caller's job (`Mycelia.gqc_cross_bh`)")
        println(io, "- k sweep: ", join(k_list, ", "))
        println(io)
        println(io, "See `summary.csv`, `formats/`, and `plots/` for details.")
    end

    params = Dict{String, Any}(
        "ks" => k_list,
        "n_shuffles" => n_shuffles,
        "shuffle_seed" => shuffle_seed,
        "salt" => salt,
        "genus_index" => genus_index,
        "min_shuffle_mean" => min_shuffle_mean,
        "bh_ks" => bh_list,
        "circular" => circ,
        "canonical" => effective_canonical,
        "strand_mode" => string(resolved_strand),
        "null_kinds" => sort(collect(string.(null_kinds))),
        "null_source" => string(null_source),
        "resolved_null_source" => string(resolved_null_source),
        "null_applicable" => null_applicable,
        "null_source_reason" => null_reason,
        "reference_n_records" => ref_load.n_records,
        "write_node_link_json" => write_node_link_json,
        "make_plots" => make_plots,
        "overwrite" => string(overwrite),
    )
    input_hashes = Dict{String, String}("graph.gfa" => sha256_file(gfa_path))
    meta_path = joinpath(input_dir, "metadata.json")
    if isfile(meta_path)
        input_hashes["metadata.json"] = sha256_file(meta_path)
    end
    ## Hash every reference file the loader matched (I1).
    for ref_file in ref_load.files_matched
        input_hashes["reference/" * basename(ref_file)] = sha256_file(ref_file)
    end

    ## Hash all artifacts except CHECKSUMS / manifest (stable identity).
    artifact_entries = Tuple{String, String}[]
    for (root, dirs, files) in walkdir(out_dir)
        sort!(dirs); sort!(files)
        for f in files
            (f == "CHECKSUMS.sha256" || f == "manifest.json") && continue
            p = joinpath(root, f)
            push!(artifact_entries, (relpath(p, out_dir), sha256_file(p)))
        end
    end
    sort!(artifact_entries; by=first)
    output_hashes = OrderedCollections.OrderedDict{String, String}(artifact_entries)

    manifest_path = joinpath(out_dir, "manifest.json")
    manifest = Dict{String, Any}(
        "gqc_version" => gqc_version,
        "commit" => _gqc_git_commit(),
        "params" => _gqc_canonicalize_params(params),
        "params_hash" => _gqc_params_hash(params),
        "input_sha256" => OrderedCollections.OrderedDict(sort(collect(input_hashes); by=first)),
        "output_sha256" => output_hashes,
        "formats" => sort([relpath(f, out_dir) for f in format_files]),
        "graph_id" => graph_id,
        "graph_id_slug" => slug,
        "graph_set" => graph_set,
        "packet_dir" => basename(out_dir),
        "metadata" => metadata,
        "null_source_reason" => null_reason,
        "reference_load" => Dict{String, Any}(
            "n_records" => ref_load.n_records,
            "files_matched" => [basename(f) for f in ref_load.files_matched],
            "note" => ref_load.note,
        ),
        "bh_note" => "dinucleotide_bh_q_within_graph is BH over bh_ks inside one graph only; use Mycelia.gqc_cross_bh for cross-genus FDR",
        "length_note" => "total_bp and unitig_* use overlap-corrected spans; total_segment_bp is raw sum of segment lengths",
        "concordance_note" => "ref_walkability is shared k-mer content (always RC-folded); not GFA-link walkability (I6)",
    )
    open(manifest_path, "w") do io
        JSON.print(io, manifest, 2)
    end

    ## CHECKSUMS: sorted by path, includes manifest.json, excludes itself.
    checksum_entries = copy(artifact_entries)
    push!(checksum_entries, ("manifest.json", sha256_file(manifest_path)))
    sort!(checksum_entries; by=first)
    open(joinpath(out_dir, "CHECKSUMS.sha256"), "w") do io
        for (rel, h) in checksum_entries
            println(io, h, "  ", rel)
        end
    end

    return out_dir
end
