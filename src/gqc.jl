# gqc v0.1 — reusable graph QC packet
# Spec: docs/gqc-v0.1-spec.md (project store). Semver 0.1.0.
# Caveat: prefer SingleStrand GFA builds until write_gfa_next / PR #471 lands.

const GQC_VERSION = "0.1.0"
const GQC_N_SHUFFLES = 100
const GQC_SHUFFLE_SEED = 20261021
const GQC_MIN_SHUFFLE_MEAN = 1.0

# =============================================================================
# Null / multiple-testing helpers (verbatim from rhizomorph-paper
# notebooks/viral-graphs.jl @ c8285d0, Oct 22 2026)
# =============================================================================

"""
    dinucleotide_shuffle(sequence::String, rng::Random.AbstractRNG)::String

Altschul–Erickson dinucleotide shuffle: preserves dinucleotide counts and the
first/last base. Ported verbatim from rhizomorph-paper `viral-graphs.jl`.
"""
function dinucleotide_shuffle(sequence::String, rng::Random.AbstractRNG)::String
    chars = collect(sequence)
    n = length(chars)
    n <= 2 && return sequence
    final = chars[end]
    exits = Dict(c => Char[] for c in "ACGT")
    for i in 1:(n - 1)
        push!(exits[chars[i]], chars[i + 1])
    end
    tree_bases = [c for c in "ACGT" if c != final && !isempty(exits[c])]
    last_exit = Dict{Char, Char}()
    reaches_final(c) = begin
        steps = 0
        while c != final && steps <= 4
            c = last_exit[c]
            steps += 1
        end
        c == final
    end
    while true
        for c in tree_bases
            last_exit[c] = rand(rng, exits[c])
        end
        all(reaches_final, tree_bases) && break
    end
    for c in "ACGT"
        if c in tree_bases
            deleteat!(exits[c], findfirst(==(last_exit[c]), exits[c]))
            Random.shuffle!(rng, exits[c])
            push!(exits[c], last_exit[c])
        else
            Random.shuffle!(rng, exits[c])
        end
    end
    next_exit = Dict(c => 1 for c in "ACGT")
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
Compact non-branching paths (undirected deg==2 chains) into unitigs.
Returns vector of unitig lengths in bp (sum of segment lengths along the path).
"""
function _gqc_unitig_lengths(g::GFAGraph, ids::Vector{String}, adj::Vector{Vector{Int}})
    n = length(ids)
    deg = _gqc_degrees(adj)
    visited = falses(n)
    lengths = Int[]
    for start in 1:n
        visited[start] && continue
        ## start unitig at tips or branch points, or any unvisited isolated/linear
        if deg[start] == 2
            ## only start at an endpoint of a linear run if neighbor not both unvisited later
            continue
        end
        path = Int[start]
        visited[start] = true
        ## extend through each neighbor direction that is linear
        for nbr0 in adj[start]
            if visited[nbr0]
                continue
            end
            prev = start
            cur = nbr0
            while true
                push!(path, cur)
                visited[cur] = true
                if deg[cur] != 2
                    break
                end
                nexts = [v for v in adj[cur] if v != prev]
                isempty(nexts) && break
                prev = cur
                cur = nexts[1]
                visited[cur] && break
            end
        end
        push!(lengths, sum(g.segments[ids[i]].length for i in path))
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
        push!(lengths, sum(g.segments[ids[i]].length for i in path))
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
        kmer_positions += n_pos
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
            occursin(r"[^ACGTacgt]", raw) && continue
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

Shuffle each segment independently `n_shuffles` times. Errors unless every copy
keeps dinucleotide counts and end bases, and all copies differ from each other
and from the original (paper `checked_shuffles` contract).
"""
function gqc_checked_dinucleotide_shuffles(sequences::Vector{String},
                                           rng::Random.AbstractRNG;
                                           n_shuffles::Int=GQC_N_SHUFFLES)
    function dinuc_counts(s::String)
        c = Dict{String, Int}()
        for i in 1:(length(s) - 1)
            d = s[i:(i + 1)]
            c[d] = get(c, d, 0) + 1
        end
        return c
    end
    results = Vector{Vector{String}}(undef, n_shuffles)
    for t in 1:n_shuffles
        shuffled_segs = String[]
        for seq in sequences
            s = dinucleotide_shuffle(seq, rng)
            if length(seq) >= 2
                dinuc_counts(s) == dinuc_counts(seq) || error("dinucleotide counts not preserved")
                s[1] == seq[1] || error("first base not preserved")
                s[end] == seq[end] || error("last base not preserved")
            end
            s == seq && length(seq) > 2 && error("shuffle identical to original")
            push!(shuffled_segs, s)
        end
        results[t] = shuffled_segs
    end
    ## all shuffle copies distinct as genome concatenations
    joined = [join(segs, '|') for segs in results]
    length(unique(joined)) == n_shuffles || error("duplicate shuffles detected")
    return results
end

"""
    gqc_repeat_null_stats(sequences; k, n_shuffles, shuffle_seed, salt, circular, canonical, min_shuffle_mean)

Paper formulas:
- `null_repeats` per shuffle = kmer_positions - distinct_kmers
- `shuffle_ratio = observed / mean(null)` if mean(null) >= min_shuffle_mean else `missing`
- `shuffle_p_upper = (1 + #{null >= observed}) / (n_shuffles + 1)`
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
    if null_kind == :dinucleotide
        shuffles = gqc_checked_dinucleotide_shuffles(sequences, rng; n_shuffles=n_shuffles)
        for segs in shuffles
            r = gqc_repeat_counts(segs; k=k, circular=circular, canonical=canonical)
            push!(null_repeats, Float64(r.repeats))
        end
    elseif null_kind == :codon
        for _ in 1:n_shuffles
            segs = [codon_shuffle(s, rng) for s in sequences]
            r = gqc_repeat_counts(segs; k=k, circular=circular, canonical=canonical)
            push!(null_repeats, Float64(r.repeats))
        end
    else
        error("unknown null_kind=$null_kind")
    end
    μ = Statistics.mean(null_repeats)
    σ = length(null_repeats) > 1 ? Statistics.std(null_repeats) : 0.0
    ratio = μ >= min_shuffle_mean ? observed.repeats / μ : missing
    z = σ > 0 ? (observed.repeats - μ) / σ : missing
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
    )
end

# =============================================================================
# Reference concordance
# =============================================================================

function _gqc_load_reference_sequences(reference_dir::AbstractString)
    seqs = String[]
    isdir(reference_dir) || return seqs
    for file in readdir(reference_dir; join=true)
        occursin(r"\.(fa|fna|fasta)(\.gz)?$"i, file) || continue
        open(file) do io
            ## handle plain or gzipped lightly via Transcoding not required; assume plain for v0.1
            buf = IOBuffer()
            for line in eachline(io)
                startswith(line, '>') && continue
                print(buf, uppercase(strip(line)))
            end
            s = String(take!(buf))
            !isempty(s) && push!(seqs, s)
        end
    end
    return seqs
end

function gqc_reference_concordance(g::GFAGraph, reference_seqs::Vector{String}; k::Int, circular::Bool=false)
    isempty(reference_seqs) && return Dict{String, Any}(
        "ref_kmer_fraction" => missing,
        "ref_path_fraction" => missing,
        "ref_absent_kmers" => missing,
        "ref_breaks" => missing,
        "ref_walkability" => missing,
    )
    seg_seqs = [uppercase(s.sequence) for s in values(g.segments) if s.sequence != "*"]
    graph_kmers = Set{String}()
    for seq in seg_seqs
        n = length(seq)
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
        n_pos = circular ? n : (n - k + 1)
        for i in 1:max(n_pos, 0)
            if circular
                raw = i + k - 1 <= n ? seq[i:(i + k - 1)] : seq[i:n] * seq[1:(k - (n - i + 1))]
            else
                raw = seq[i:(i + k - 1)]
            end
            occursin(r"[^ACGT]", raw) && continue
            ck = _gqc_canonical_kmer(uppercase(raw))
            push!(ref_kmers, ck)
            ref_positions += 1
            ck in graph_kmers && (present += 1)
        end
    end
    absent = length(setdiff(ref_kmers, graph_kmers))
    frac = ref_positions > 0 ? present / ref_positions : missing
    ## path fraction / walkability: fraction of consecutive ref kmer pairs that are both in graph
    walkable = 0
    pair_total = 0
    for seq in reference_seqs
        n = length(seq)
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

function gqc_write_formats(g::GFAGraph, formats_dir::AbstractString; write_node_link_json::Bool=true)
    mkpath(formats_dir)
    ids, adj, _ = _gqc_undirected_adjacency(g)
    n = length(ids)
    ## nodes.tsv
    nodes_path = joinpath(formats_dir, "nodes.tsv")
    open(nodes_path, "w") do io
        println(io, "node_id\tindex\tlength\tdepth\tsequence")
        for (i, id) in enumerate(ids)
            s = g.segments[id]
            seq = s.sequence == "*" ? "" : s.sequence
            d = s.depth === nothing ? "" : string(s.depth)
            println(io, id, '\t', i, '\t', s.length, '\t', d, '\t', seq)
        end
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
    ## Matrix Market coordinate pattern (undirected, 1-based, row/col order = nodes.tsv)
    index = Dict(id => i for (i, id) in enumerate(ids))
    edges_unique = Set{Tuple{Int,Int}}()
    for link in g.links
        haskey(index, link.from) || continue
        haskey(index, link.to) || continue
        u = index[link.from]; v = index[link.to]
        a, b = minmax(u, v)
        push!(edges_unique, (a, b))
    end
    mtx_path = joinpath(formats_dir, "adjacency.mtx")
    open(mtx_path, "w") do io
        println(io, "%%MatrixMarket matrix coordinate pattern general")
        println(io, "% row/col order matches nodes.tsv index column (1-based)")
        println(io, n, ' ', n, ' ', length(edges_unique))
        for (a, b) in sort!(collect(edges_unique))
            println(io, a, ' ', b)
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
    ## degree log-log histogram
    if !isempty(degree_hist)
        degs = collect(keys(degree_hist))
        counts = [degree_hist[d] for d in degs]
        p = Plots.plot(degs, counts; seriestype=:scatter, xscale=:log10, yscale=:log10,
                       xlabel="degree", ylabel="count", title="Degree distribution", legend=false)
        _gqc_save_plot(p, plots_dir, "degree", graph_id, k)
    end
    ## unitig cumulative / N50
    if !isempty(unitig_lengths)
        sorted = sort(unitig_lengths; rev=true)
        cum = cumsum(sorted)
        p = Plots.plot(1:length(cum), cum; xlabel="unitig rank", ylabel="cumulative bp",
                       title="Cumulative unitig length", legend=false)
        n50 = summary_row["unitig_n50"]
        Plots.hline!(p, [sum(sorted) / 2]; label="half total", linestyle=:dash)
        _gqc_save_plot(p, plots_dir, "unitig_n50", graph_id, k)
    end
    ## coverage histogram
    if !isempty(coverages)
        p = Plots.histogram(coverages; xlabel="depth", ylabel="count", title="Coverage", legend=false)
        _gqc_save_plot(p, plots_dir, "coverage", graph_id, k)
    end
    ## tip/bubble/branch bars
    bars = ["tips", "bubbles", "branch", "self_loops", "tangles"]
    vals = [summary_row[k_] for k_ in ("n_tips", "n_bubbles", "n_branch", "n_self_loops", "n_tangles")]
    p = Plots.bar(bars, vals; legend=false, title="Structural counts", xrotation=30)
    _gqc_save_plot(p, plots_dir, "structure", graph_id, k)
    ## repeat ratio vs null CI (if present)
    if null_ratio_ci !== nothing
        ratio, lo, hi = null_ratio_ci
        p = Plots.scatter([1], [ratio]; yerror=[(ratio - lo, hi - ratio)],
                          xlabel="", ylabel="repeat ratio", title="Repeat vs null",
                          xticks=([], []), legend=false, xlim=(0.5, 1.5))
        Plots.hline!(p, [1.0]; linestyle=:dash)
        _gqc_save_plot(p, plots_dir, "repeat_null", graph_id, k)
    end
    ## compacted small multiple for tiny graphs
    if g !== nothing && length(g.segments) <= 5000 && adj !== nothing
        ## simple component-size rank plot doubles as compact overview
        comps = _gqc_connected_components(adj)
        sizes = sort([length(c) for c in comps]; rev=true)
        p = Plots.plot(1:length(sizes), sizes; xlabel="component rank", ylabel="nodes",
                       title="Component sizes", legend=false)
        _gqc_save_plot(p, plots_dir, "components", graph_id, k)
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

function _gqc_params_hash(params::Dict)
    canonical = JSON.json(sort(collect(params); by=first))
    return SHA.bytes2hex(SHA.sha256(canonical))
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
    run_gqc(input_dir, output_root; kwargs...) -> String

Run gqc v0.1 on one graph folder and write a QC packet under
`output_root/<graph_set>/`.

# Arguments
- `input_dir`: folder containing `graph.gfa`, `metadata.json`, optional `reference/*.fasta`
- `output_root`: parent directory for packet output

# Keyword arguments
- `ks`: k values for repeat-null / reference concordance sweep (default: `[metadata["k"]]`)
- `n_shuffles`: dinucleotide/codon null replicates (default: `100`)
- `shuffle_seed`: base RNG seed (default: `20261021`; per-graph RNG =
  `Xoshiro(shuffle_seed + 1000*salt + genus_index)`)
- `salt`, `genus_index`: seeding offsets (paper formula)
- `min_shuffle_mean`: minimum null mean for ratio (default: `1.0`)
- `bh_ks`: k values included in BH correction (default: all `ks`)
- `circular`, `canonical`: k-mer counting flags
- `strand_mode`: `:single` (default; recommended until PR #471) or `:double`
- `write_node_link_json`: include `formats/node_link.json` (default: `true`)
- `make_plots`: write plots under `plots/` (default: `true`)
- `null_kinds`: `(:dinucleotide, :codon)` by default

# Returns
Absolute path to the packet directory `output_root/<graph_set>/`.
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
                 canonical::Bool=true,
                 strand_mode::Symbol=:single,
                 write_node_link_json::Bool=true,
                 make_plots::Bool=true,
                 null_kinds::Tuple=(:dinucleotide, :codon),
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
    k_list = isempty(ks) ? (k_meta > 0 ? [k_meta] : Int[]) : collect(Int, ks)
    isempty(k_list) && error("no k provided: set metadata[\"k\"] or pass ks=")

    g = read_gfa_records(gfa_path)
    ids, adj, self_loops = _gqc_undirected_adjacency(g)
    comps = _gqc_connected_components(adj)
    degrees = _gqc_degrees(adj)
    roles = _gqc_node_roles(adj)
    indeg, outdeg, sources, sinks = _gqc_directed_degrees(g)
    bubbles = _gqc_simple_bubbles(adj)
    tangles = _gqc_tangle_count(adj)
    unitig_lengths = _gqc_unitig_lengths(g, ids, adj)
    n50stats = _gqc_n50_l50(unitig_lengths)
    cov_q = _gqc_coverage_quantiles(g)
    coverages = Float64[s.depth for s in values(g.segments) if s.depth !== nothing]
    n_seg = length(g.segments)
    n_link = length(g.links)
    total_bp = sum(s.length for s in values(g.segments); init=0)
    density = n_seg > 1 ? (2 * n_link) / (n_seg * (n_seg - 1)) : 0.0
    largest_frac = n_seg > 0 ? maximum(length(c) for c in comps) / n_seg : 0.0
    ## Betti numbers on undirected simple graph: β0 = components, β1 = e - v + c (cyclomatic)
    n_undirected_edges = sum(length(a) for a in adj) ÷ 2
    beta0 = length(comps)
    beta1 = n_undirected_edges - n_seg + beta0

    ref_dir = joinpath(input_dir, "reference")
    ref_seqs = _gqc_load_reference_sequences(ref_dir)
    seg_seqs = [uppercase(s.sequence) for s in values(g.segments) if s.sequence != "*" && !isempty(s.sequence)]

    out_dir = joinpath(output_root, graph_set)
    mkpath(out_dir)
    for sub in ("nodes", "edges", "unitigs", "nulls", "plots", "formats")
        mkpath(joinpath(out_dir, sub))
    end

    ## long tables
    open(joinpath(out_dir, "nodes", "$(graph_id).tsv"), "w") do io
        println(io, "graph_id\tnode_id\tlength\tdepth\tdegree\tin_degree\tout_degree")
        for (i, id) in enumerate(ids)
            s = g.segments[id]
            d = s.depth === nothing ? "" : string(s.depth)
            println(io, graph_id, '\t', id, '\t', s.length, '\t', d, '\t',
                    degrees[i], '\t', indeg[id], '\t', outdeg[id])
        end
    end
    open(joinpath(out_dir, "edges", "$(graph_id).tsv"), "w") do io
        println(io, "graph_id\tfrom\tfrom_orient\tto\tto_orient\toverlap")
        for link in g.links
            println(io, graph_id, '\t', link.from, '\t', link.from_forward ? "+" : "-", '\t',
                    link.to, '\t', link.to_forward ? "+" : "-", '\t', link.overlap)
        end
    end
    open(joinpath(out_dir, "unitigs", "$(graph_id).tsv"), "w") do io
        println(io, "graph_id\tunitig_index\tlength")
        for (i, L) in enumerate(sort(unitig_lengths; rev=true))
            println(io, graph_id, '\t', i, '\t', L)
        end
    end

    format_files = gqc_write_formats(g, joinpath(out_dir, "formats"); write_node_link_json=write_node_link_json)

    ## per-k summary rows + nulls
    summary_rows = Dict{String, Any}[]
    null_p_dinuc = Float64[]
    null_k_order = Int[]
    bh_set = bh_ks === nothing ? Set(k_list) : Set(Int(k) for k in bh_ks)

    for k in k_list
        row = Dict{String, Any}(
            "graph_id" => graph_id,
            "k" => k,
            "n_segments" => n_seg,
            "n_links" => n_link,
            "total_bp" => total_bp,
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
            "strand_mode" => string(strand_mode),
            "gqc_version" => gqc_version,
        )
        merge!(row, Dict{String, Any}(cov_q))
        ref_stats = gqc_reference_concordance(g, ref_seqs; k=k, circular=circ)
        merge!(row, ref_stats)

        null_ratio_ci = nothing
        for kind in null_kinds
            stats = gqc_repeat_null_stats(seg_seqs;
                                          k=k,
                                          n_shuffles=n_shuffles,
                                          shuffle_seed=shuffle_seed,
                                          salt=salt,
                                          genus_index=genus_index,
                                          circular=circ,
                                          canonical=canonical,
                                          min_shuffle_mean=min_shuffle_mean,
                                          null_kind=kind)
            prefix = string(kind)
            row["$(prefix)_observed_repeats"] = stats.observed_repeats
            row["$(prefix)_null_mean"] = stats.null_mean
            row["$(prefix)_null_std"] = stats.null_std
            row["$(prefix)_shuffle_ratio"] = stats.shuffle_ratio
            row["$(prefix)_z_score"] = stats.z_score
            row["$(prefix)_shuffle_p_upper"] = stats.shuffle_p_upper
            ## null long table
            open(joinpath(out_dir, "nulls", "$(graph_id)__k$(k)__$(prefix).tsv"), "w") do io
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

        ## hooks for panlingual consumers
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
                           joinpath(out_dir, "plots"), graph_id, k; adj=adj, ids=ids, g=g)
            ## beta0/beta1 vs k will be written after sweep as a combined plot
        end
    end

    ## BH across selected k for dinucleotide p_upper
    if !isempty(null_p_dinuc)
        qvals = benjamini_hochberg(null_p_dinuc)
        qmap = Dict(null_k_order[i] => qvals[i] for i in eachindex(null_k_order))
        for row in summary_rows
            row["dinucleotide_bh_q"] = get(qmap, row["k"], missing)
        end
    end

    if make_plots && length(summary_rows) >= 1
        ks_plot = [r["k"] for r in summary_rows]
        b0 = [r["beta0"] for r in summary_rows]
        b1 = [r["beta1"] for r in summary_rows]
        p = Plots.plot(ks_plot, b0; label="β0", xlabel="k", ylabel="Betti", title="Betti vs k")
        Plots.plot!(p, ks_plot, b1; label="β1")
        _gqc_save_plot(p, joinpath(out_dir, "plots"), "betti_vs_k", graph_id, first(ks_plot))
    end

    ## summary.csv / summary.json
    all_keys = OrderedCollections.OrderedSet{String}()
    for row in summary_rows
        for k in keys(row)
            push!(all_keys, k)
        end
    end
    ## drop bulky vector fields from CSV
    csv_keys = [k for k in all_keys if k != "graph_feature_vector" && k != "distance_to_null"]
    open(joinpath(out_dir, "summary.csv"), "w") do io
        println(io, join(csv_keys, ','))
        for row in summary_rows
            vals = [repr(get(row, k, "")) for k in csv_keys]
            ## simpler CSV: string without Julia repr quotes for numbers
            rendered = String[]
            for k in csv_keys
                v = get(row, k, "")
                if v === missing
                    push!(rendered, "")
                elseif v isa AbstractString
                    push!(rendered, '"' * replace(v, '"' => "\"\"") * '"')
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
        println(io, "- segments: **$(n_seg)**; links: **$(n_link)**; bp: **$(total_bp)**")
        println(io, "- components (β0): **$(beta0)**; β1: **$(beta1)**; largest fraction: **$(round(largest_frac; digits=3))**")
        println(io, "- tips: **$(length(roles.tips))**; bubbles: **$(bubbles)**; tangles: **$(tangles)**")
        println(io, "- unitigs: **$(length(unitig_lengths))**; N50: **$(n50stats.n50)**; L50: **$(n50stats.l50)**")
        println(io, "- strand_mode: `$(strand_mode)` (single-strand recommended until PR #471)")
        println(io, "- k sweep: ", join(k_list, ", "))
        println(io)
        println(io, "See `summary.csv`, `formats/`, and `plots/` for details.")
    end

    params = Dict(
        "ks" => k_list,
        "n_shuffles" => n_shuffles,
        "shuffle_seed" => shuffle_seed,
        "salt" => salt,
        "genus_index" => genus_index,
        "min_shuffle_mean" => min_shuffle_mean,
        "bh_ks" => collect(bh_set),
        "circular" => circ,
        "canonical" => canonical,
        "strand_mode" => string(strand_mode),
        "null_kinds" => collect(string.(null_kinds)),
        "write_node_link_json" => write_node_link_json,
        "make_plots" => make_plots,
    )
    input_hashes = Dict("graph.gfa" => sha256_file(gfa_path))
    meta_path = joinpath(input_dir, "metadata.json")
    if isfile(meta_path)
        input_hashes["metadata.json"] = sha256_file(meta_path)
    end

    ## checksums for outputs
    checksum_lines = String[]
    function _add_checksum!(path)
        isfile(path) || return
        rel = relpath(path, out_dir)
        h = sha256_file(path)
        push!(checksum_lines, "$(h)  $(rel)")
        return h
    end
    output_hashes = Dict{String, String}()
    for (root, _, files) in walkdir(out_dir)
        for f in files
            f == "CHECKSUMS.sha256" && continue
            f == "manifest.json" && continue
            p = joinpath(root, f)
            output_hashes[relpath(p, out_dir)] = _add_checksum!(p)
        end
    end
    open(joinpath(out_dir, "CHECKSUMS.sha256"), "w") do io
        for line in sort(checksum_lines)
            println(io, line)
        end
    end

    manifest = Dict(
        "gqc_version" => gqc_version,
        "commit" => _gqc_git_commit(),
        "params" => params,
        "params_hash" => _gqc_params_hash(params),
        "input_sha256" => input_hashes,
        "output_sha256" => output_hashes,
        "formats" => [relpath(f, out_dir) for f in format_files],
        "graph_id" => graph_id,
        "graph_set" => graph_set,
        "metadata" => metadata,
    )
    open(joinpath(out_dir, "manifest.json"), "w") do io
        JSON.print(io, manifest, 2)
    end
    ## refresh checksums to include manifest
    open(joinpath(out_dir, "CHECKSUMS.sha256"), "a") do io
        println(io, sha256_file(joinpath(out_dir, "manifest.json")), "  manifest.json")
    end

    return out_dir
end
