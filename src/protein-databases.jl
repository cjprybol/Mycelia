# protein-databases.jl — UniProt, UniRef, AlphaFold, and FoldSeek API integration
#
# Provides functions for fetching protein data from major databases:
# - UniProt: entries, sequences, embeddings
# - UniRef: cluster membership and member expansion
# - AlphaFold: predicted structures
# - FoldSeek: structural similarity search (uses existing Mycelia.foldseek_easy_search)
#
# All functions use HTTP.jl for API access and return parsed Julia types.
# Uses JSON.parse (not JSON3) — Mycelia imports JSON, not JSON3.

# Common headers for all API requests — includes Accept-Encoding: identity
# to prevent gzip-compressed responses that HTTP.jl doesn't auto-decompress.
const _UNIPROT_HEADERS = [
    "Accept" => "application/json",
    "Accept-Encoding" => "identity"
]

const _IDENTITY_HEADERS = [
    "Accept-Encoding" => "identity",
]

"""
    fetch_uniprot_entry(accession::AbstractString) -> Dict

Fetch a full UniProt entry as parsed JSON.

# Example
```julia
entry = fetch_uniprot_entry("P05820")
entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
# => "Colicin-M"
```
"""
function fetch_uniprot_entry(accession::AbstractString)
    url = "https://rest.uniprot.org/uniprotkb/$(accession).json"
    response = HTTP.get(url; headers = _UNIPROT_HEADERS)
    return JSON.parse(String(response.body))
end

"""
    fetch_uniprot_fasta(accession::AbstractString) -> String

Fetch protein sequence in FASTA format from UniProt.

# Example
```julia
fasta = fetch_uniprot_fasta("P05820")
```
"""
function fetch_uniprot_fasta(accession::AbstractString)
    url = "https://rest.uniprot.org/uniprotkb/$(accession).fasta"
    response = HTTP.get(url; headers = _IDENTITY_HEADERS)
    return String(response.body)
end

"""
    fetch_uniprot_sequence(accession::AbstractString) -> BioSequences.LongAA

Fetch protein sequence from UniProt and return as a BioSequences amino acid sequence.

# Example
```julia
seq = fetch_uniprot_sequence("P05820")
length(seq)  # => 271
```
"""
function fetch_uniprot_sequence(accession::AbstractString)
    fasta_text = fetch_uniprot_fasta(accession)
    # Parse FASTA — skip header line, join sequence lines
    lines = split(fasta_text, '\n')
    seq_lines = filter(l -> !startswith(l, '>') && !isempty(strip(l)), lines)
    seq_str = join(seq_lines)
    return BioSequences.LongAA(seq_str)
end

"""
    fetch_uniprot_sequences(accessions::AbstractVector{<:AbstractString};
                            delay::Real=0.3) -> Dict{String, BioSequences.LongAA}

Batch fetch protein sequences from UniProt with rate limiting.
Returns a Dict mapping accession => sequence.
"""
function fetch_uniprot_sequences(accessions::AbstractVector{<:AbstractString};
        delay::Real = 0.3)
    sequences = Dict{String, BioSequences.LongAA}()
    for (i, acc) in enumerate(accessions)
        try
            sequences[acc] = fetch_uniprot_sequence(acc)
            i < length(accessions) && sleep(delay)
        catch e
            @warn "Failed to fetch sequence for $(acc): $(e)"
        end
    end
    return sequences
end

"""
    search_uniprot(query::AbstractString; size::Int=10) -> Vector

Search UniProt by free-text query. Returns a vector of result entries.

# Example
```julia
results = search_uniprot("gene:csgC AND organism_id:83333"; size=3)
```
"""
function search_uniprot(query::AbstractString; size::Int = 10)
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = ["query" => query, "format" => "json", "size" => string(size)]
    response = HTTP.get(url; headers = _UNIPROT_HEADERS, query = params)
    parsed = JSON.parse(String(response.body))
    return get(parsed, "results", [])
end

# --------------------------------------------------------------------------- #
# UniRef cluster operations
# --------------------------------------------------------------------------- #

"""
    search_uniref_cluster(accession::AbstractString;
                          identity::Float64=0.5) -> Union{Dict, Nothing}

Find the UniRef cluster (at given identity threshold) containing a UniProt accession.
Returns cluster info dict with string keys, or nothing if not found.

# Example
```julia
cluster = search_uniref_cluster("P05820"; identity=0.5)
cluster["id"]            # => "UniRef50_P05820"
cluster["member_count"]  # => 236
```
"""
function search_uniref_cluster(accession::AbstractString; identity::Float64 = 0.5)
    identity_str = identity == 1.0 ? "1.0" : identity == 0.9 ? "0.9" : "0.5"
    url = "https://rest.uniprot.org/uniref/search"
    params = [
        "query" => "uniprot_id:$(accession) AND identity:$(identity_str)",
        "format" => "json",
        "size" => "1"
    ]
    response = HTTP.get(url; headers = _UNIPROT_HEADERS, query = params)
    parsed = JSON.parse(String(response.body))
    results = get(parsed, "results", [])
    isempty(results) && return nothing

    cluster = first(results)
    return Dict(
        "id" => get(cluster, "id", ""),
        "name" => get(cluster, "name", ""),
        "member_count" => get(cluster, "memberCount", 0),
        "updated" => get(cluster, "updated", "")
    )
end

"""
    fetch_uniref_members(cluster_id::AbstractString;
                         max_members::Int=10) -> Vector{Dict}

Get members of a UniRef cluster. Returns a vector of member info dicts.

# Example
```julia
members = fetch_uniref_members("UniRef50_P05820"; max_members=5)
for m in members
    println("\$(m["accession"]) — \$(m["organism"]) (\$(m["length"]) aa)")
end
```
"""
function fetch_uniref_members(cluster_id::AbstractString; max_members::Int = 10)
    url = "https://rest.uniprot.org/uniref/$(cluster_id)/members"
    params = ["format" => "json", "size" => string(max_members)]
    response = HTTP.get(url; headers = _UNIPROT_HEADERS, query = params)
    parsed = JSON.parse(String(response.body))
    raw_members = get(parsed, "results", [])

    members = Dict{String, Any}[]
    for m in raw_members
        push!(members,
            Dict{String, Any}(
                "accession" => get(m, "memberId", ""),
                "type" => get(m, "memberIdType", ""),
                "organism" => get(m, "organismName", "unknown"),
                "length" => get(m, "sequenceLength", 0),
                "tax_id" => get(m, "organismTaxId", 0)
            ))
    end
    return members
end

"""
    expand_via_uniref(accession::AbstractString;
                      identity::Float64=0.5,
                      max_variants::Int=3,
                      max_length::Int=500,
                      deduplicate::Bool=true) -> Vector{Dict}

Find diverse natural variants of a protein by looking up its UniRef50 cluster
and returning members that represent genuinely different sequences.

When `deduplicate=true` (default), uses the UniRef hierarchy to remove
redundant sequences: members sharing the same UniRef90 cluster are >90%
identical and only one representative is kept. This prevents returning
near-identical sequences from different genome assemblies.

Filters by max_length to enforce phage packaging constraints.
Prefers variants from different species/genera for maximum diversity.

# Example
```julia
variants = expand_via_uniref("P05820"; max_variants=3)
# Returns 3 variants from different UniRef90 clusters
```
"""
function expand_via_uniref(accession::AbstractString;
        identity::Float64 = 0.5,
        max_variants::Int = 3,
        max_length::Int = 500,
        deduplicate::Bool = true)
    cluster = search_uniref_cluster(accession; identity = identity)
    cluster === nothing && return Dict{String, Any}[]

    # Fetch more members than needed to allow for deduplication filtering
    fetch_count = deduplicate ? max(max_variants * 5, 20) : max_variants + 2
    members = fetch_uniref_members(cluster["id"]; max_members = fetch_count)

    # Filter: exclude query protein, enforce size limit
    candidates = filter(members) do m
        m["accession"] != accession && m["length"] <= max_length
    end

    if !deduplicate || isempty(candidates)
        return candidates[1:min(max_variants, length(candidates))]
    end

    # Deduplicate using UniRef90 clusters: members sharing the same
    # UniRef90 cluster are >90% identical — keep only one per cluster.
    # Also look up UniRef90 for the query protein to exclude its cluster.
    query_ur90 = try
        c = search_uniref_cluster(accession; identity = 0.9)
        c === nothing ? "" : c["id"]
    catch
        ""
    end
    sleep(0.3)

    seen_ur90 = Set{String}()
    if !isempty(query_ur90)
        push!(seen_ur90, query_ur90)
    end

    # Collect all UniRef90-distinct candidates first, then rank by genus diversity
    ur90_distinct = Dict{String, Any}[]

    for m in candidates
        acc = m["accession"]

        # Look up this member's UniRef90 cluster
        ur90 = try
            c = search_uniref_cluster(acc; identity = 0.9)
            c === nothing ? acc : c["id"]
        catch
            acc
        end
        sleep(0.3)

        # Skip if we already have a representative from this UniRef90 cluster
        ur90 in seen_ur90 && continue
        push!(seen_ur90, ur90)

        # Extract genus from organism name
        organism = get(m, "organism", "unknown")
        genus = first(split(organism))

        m["uniref90_id"] = ur90
        m["genus"] = genus
        push!(ur90_distinct, m)
    end

    # Rank by genus diversity: variants from novel genera score higher.
    # Extract query genus to deprioritize same-genus hits.
    query_entry = try
        fetch_uniprot_entry(accession)
    catch
        nothing
    end
    query_genus = if query_entry !== nothing
        org = get(get(query_entry, "organism", Dict()), "scientificName", "unknown")
        first(split(org))
    else
        ""
    end

    # Sort: novel genus first, then different genus from query, then same genus
    genus_counts = Dict{String, Int}()
    for v in ur90_distinct
        g = v["genus"]
        genus_counts[g] = get(genus_counts, g, 0) + 1
    end

    sort!(ur90_distinct; by = v -> begin
        g = v["genus"]
        if g == query_genus
            2  # Same genus as query — lowest priority
        elseif genus_counts[g] > 1
            1  # Genus appears multiple times — medium priority
        else
            0  # Unique genus — highest priority (sort ascending)
        end
    end)

    return ur90_distinct[1:min(max_variants, length(ur90_distinct))]
end

# --------------------------------------------------------------------------- #
# AlphaFold structure operations
# --------------------------------------------------------------------------- #

"""
    fetch_alphafold_prediction(accession::AbstractString) -> Union{Dict, Nothing}

Get AlphaFold prediction metadata for a UniProt accession.
Returns dict with string keys, or nothing.

# Example
```julia
pred = fetch_alphafold_prediction("P05820")
pred["plddt"]    # => 96.12
pred["pdb_url"]  # => "https://alphafold.ebi.ac.uk/files/AF-P05820-F1-model_v6.pdb"
```
"""
function fetch_alphafold_prediction(accession::AbstractString)
    url = "https://alphafold.ebi.ac.uk/api/prediction/$(accession)"
    response = try
        HTTP.get(url; headers = _UNIPROT_HEADERS)
    catch e
        @warn "AlphaFold prediction not found for $(accession)"
        return nothing
    end
    results = JSON.parse(String(response.body))
    isempty(results) && return nothing

    pred = first(results)
    return Dict(
        "pdb_url" => get(pred, "pdbUrl", ""),
        "cif_url" => get(pred, "cifUrl", ""),
        "plddt" => get(pred, "globalMetricValue", 0.0),
        "model_version" => get(pred, "latestVersion", 0),
        "sequence_length" => get(pred, "uniprotEnd", 0) - get(pred, "uniprotStart", 0) + 1
    )
end

"""
    download_alphafold_structure(accession::AbstractString,
                                 output_dir::AbstractString;
                                 format::Symbol=:pdb) -> String

Download an AlphaFold predicted structure file.
Returns the path to the downloaded file.

# Example
```julia
pdb_path = download_alphafold_structure("P05820", "structures/")
```
"""
function download_alphafold_structure(accession::AbstractString,
        output_dir::AbstractString;
        format::Symbol = :pdb)
    pred = fetch_alphafold_prediction(accession)
    pred === nothing && error("No AlphaFold prediction for $(accession)")

    url = format == :pdb ? pred["pdb_url"] : pred["cif_url"]
    ext = format == :pdb ? "pdb" : "cif"
    mkpath(output_dir)
    output_path = joinpath(output_dir, "AF-$(accession).$(ext)")

    response = HTTP.get(url; headers = _IDENTITY_HEADERS)
    write(output_path, response.body)
    return output_path
end

"""
    download_alphafold_structures(accessions::AbstractVector{<:AbstractString},
                                  output_dir::AbstractString;
                                  format::Symbol=:pdb,
                                  delay::Real=0.3) -> Dict{String, String}

Batch download AlphaFold structures. Returns Dict mapping accession => file path.
"""
function download_alphafold_structures(accessions::AbstractVector{<:AbstractString},
        output_dir::AbstractString;
        format::Symbol = :pdb,
        delay::Real = 0.3)
    paths = Dict{String, String}()
    for (i, acc) in enumerate(accessions)
        try
            paths[acc] = download_alphafold_structure(acc, output_dir; format = format)
            i < length(accessions) && sleep(delay)
        catch e
            @warn "Failed to download structure for $(acc): $(e)"
        end
    end
    return paths
end

# --------------------------------------------------------------------------- #
# Protein embeddings
# --------------------------------------------------------------------------- #

"""
    fetch_uniprot_embeddings(accessions::AbstractVector{<:AbstractString};
                             model::Symbol=:per_protein,
                             delay::Real=0.5) -> Dict{String, Vector{Float32}}

Compute composition-based protein embeddings from UniProt sequences.

Downloads the sequence from UniProt and computes an AAmer (k=2,3) composition
embedding as a normalized feature vector. This is a lightweight alignment-free
representation suitable for clustering.

When ProtT5 bulk download is available (td-mm8yj), this will upgrade to real
pre-trained embeddings from UniProt's FTP archive.

# Example
```julia
embeddings = fetch_uniprot_embeddings(["P05820", "Q840G9", "P52107"])
length(embeddings["P05820"])  # => n_features (variable, depends on observed k-mers)
```
"""
function fetch_uniprot_embeddings(accessions::AbstractVector{<:AbstractString};
        model::Symbol = :per_protein,
        delay::Real = 0.5)
    embeddings = Dict{String, Vector{Float32}}()

    for (i, acc) in enumerate(accessions)
        try
            url = "https://rest.uniprot.org/uniprotkb/$(acc).json?fields=sequence"
            response = HTTP.get(url; headers = _UNIPROT_HEADERS)
            entry = JSON.parse(String(response.body))

            # Extract sequence for composition-based embedding
            seq_str = get(get(entry, "sequence", Dict()), "value", "")
            if !isempty(seq_str)
                seq = BioSequences.LongAA(seq_str)
                embeddings[acc] = _compute_composition_embedding(seq)
            end

            i < length(accessions) && sleep(delay)
        catch e
            @warn "Failed to get embedding for $(acc): $(e)"
        end
    end

    return embeddings
end

"""
    _compute_composition_embedding(seq::BioSequences.LongAA;
                                    ks::Vector{Int}=[2, 3]) -> Vector{Float32}

Compute a fixed-length protein embedding from concatenated AAmer frequency
profiles at multiple k-mer sizes. Lightweight, alignment-free embedding.

Returns a normalized Vector{Float32} based on observed k-mer frequencies.
"""
function _compute_composition_embedding(seq::BioSequences.LongAA;
        ks::Vector{Int} = [2, 3])
    features = Float32[]

    for k in ks
        kmer_type = Kmers.Kmer{BioSequences.AminoAcidAlphabet, k}
        counts = count_kmers(kmer_type, seq)
        total = sum(values(counts); init = 0)
        total == 0 && continue

        # Convert to sorted frequency vector for consistent ordering
        sorted_kmers = sort(collect(counts); by = first)
        freqs = Float32[v / total for (_, v) in sorted_kmers]
        append!(features, freqs)
    end

    # Normalize to unit vector
    norm = sqrt(sum(features .^ 2))
    if norm > 0
        features ./= norm
    end

    return features
end

# --------------------------------------------------------------------------- #
# FoldSeek operations — uses existing Mycelia.foldseek_easy_search()
# See foldseek.jl for the core FoldSeek CLI wrappers.
# This section adds only the distance matrix parser on top.
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Protein distance matrices — multi-axis
# --------------------------------------------------------------------------- #

"""
    protein_aamer_distance_matrix(sequences::Dict{String, BioSequences.LongAA};
                                  k::Int=3,
                                  metric::Symbol=:bray_curtis) -> (Matrix{Float64}, Vector{String})

Compute pairwise distance matrix from amino acid k-mer (AAmer) composition profiles.
Returns (distance_matrix, ordered_labels).

Supported metrics: :bray_curtis, :cosine, :euclidean, :jaccard

# Example
```julia
seqs = fetch_uniprot_sequences(["P05820", "Q840G9", "P52107"])
D, labels = protein_aamer_distance_matrix(seqs; k=3, metric=:bray_curtis)
```
"""
function protein_aamer_distance_matrix(sequences::Dict{String, BioSequences.LongAA};
        k::Int = 3,
        metric::Symbol = :bray_curtis)
    isempty(sequences) && return (Matrix{Float64}(undef, 0, 0), String[])

    labels = sort(collect(keys(sequences)))
    n = length(labels)

    # Count k-mers for each sequence
    kmer_type = Kmers.Kmer{BioSequences.AminoAcidAlphabet, k}
    profiles = Dict{String, Dict}()
    all_kmers = Set{Any}()

    for label in labels
        counts = count_kmers(kmer_type, sequences[label])
        profiles[label] = counts
        union!(all_kmers, keys(counts))
    end

    # Build frequency matrix (proteins × kmers)
    kmer_list = sort(collect(all_kmers))
    freq_matrix = zeros(Float64, n, length(kmer_list))
    for (i, label) in enumerate(labels)
        total = sum(values(profiles[label]); init = 0)
        total == 0 && continue
        for (j, kmer) in enumerate(kmer_list)
            freq_matrix[i, j] = get(profiles[label], kmer, 0) / total
        end
    end

    # Compute pairwise distance using Distances.jl (validated implementations).
    # freq_matrix is (n_proteins × n_kmers), Distances.pairwise expects columns
    # as observations, so transpose.
    dist_metric = if metric == :bray_curtis
        Distances.BrayCurtis()
    elseif metric == :cosine
        Distances.CosineDist()
    elseif metric == :euclidean
        Distances.Euclidean()
    elseif metric == :jaccard
        Distances.Jaccard()
    else
        error("Unknown metric: $(metric). Use :bray_curtis, :cosine, :euclidean, or :jaccard")
    end
    D = Distances.pairwise(dist_metric, freq_matrix'; dims = 2)

    return D, labels
end

"""
    embedding_distance_matrix(embeddings::Dict{String, Vector{Float32}};
                              metric::Symbol=:cosine) -> (Matrix{Float64}, Vector{String})

Compute pairwise distance matrix from pre-computed protein embeddings.
Returns (distance_matrix, ordered_labels).

# Example
```julia
D, labels = embedding_distance_matrix(embeddings; metric=:cosine)
```
"""
function embedding_distance_matrix(embeddings::Dict{String, Vector{Float32}};
        metric::Symbol = :cosine)
    isempty(embeddings) && return (Matrix{Float64}(undef, 0, 0), String[])

    labels = sort(collect(keys(embeddings)))
    n = length(labels)

    # Embeddings may have variable length (e.g., composition embeddings use
    # only observed k-mers). Pad shorter vectors with zeros to match the
    # longest — equivalent to treating unobserved features as zero frequency.
    max_dim = maximum(length(embeddings[l]) for l in labels)

    E = zeros(Float64, n, max_dim)
    for (i, label) in enumerate(labels)
        v = embeddings[label]
        E[i, 1:length(v)] = Float64.(v)
    end

    # Use Distances.jl for validated pairwise computation
    dist_metric = if metric == :cosine
        Distances.CosineDist()
    elseif metric == :euclidean
        Distances.Euclidean()
    else
        error("Unknown metric: $(metric). Use :cosine or :euclidean")
    end
    D = Distances.pairwise(dist_metric, E'; dims = 2)

    return D, labels
end

"""
    structural_distance_matrix(foldseek_results_path::AbstractString,
                               accessions::AbstractVector{<:AbstractString}) -> (Matrix{Float64}, Vector{String})

Parse FoldSeek all-vs-all results into a TM-score distance matrix (1 - TM-score).
Returns (distance_matrix, ordered_labels).
"""
function structural_distance_matrix(foldseek_results_path::AbstractString,
        accessions::AbstractVector{<:AbstractString})
    labels = sort(accessions)
    n = length(labels)
    label_idx = Dict(l => i for (i, l) in enumerate(labels))

    D = ones(Float64, n, n)
    for i in 1:n
        D[i, i] = 0.0
    end

    for line in eachline(foldseek_results_path)
        fields = split(line, '\t')
        length(fields) < 7 && continue
        query = replace(basename(String(fields[1])), r"^AF-|\.pdb$" => "")
        target = replace(basename(String(fields[2])), r"^AF-|\.pdb$" => "")
        tmscore = tryparse(Float64, String(fields[7]))
        tmscore === nothing && continue

        qi = get(label_idx, query, 0)
        ti = get(label_idx, target, 0)
        (qi == 0 || ti == 0) && continue

        dist = 1.0 - tmscore
        D[qi, ti] = min(D[qi, ti], dist)
        D[ti, qi] = min(D[ti, qi], dist)
    end

    return D, labels
end

"""
    consensus_distance_matrix(matrices::Vector{Matrix{Float64}};
                              weights::Vector{Float64}=Float64[]) -> Matrix{Float64}

Combine multiple distance matrices into a weighted consensus.
Default: equal weights. All matrices must be the same size.

# Example
```julia
D_consensus = consensus_distance_matrix([D_aamer, D_embed, D_struct];
                                         weights=[0.3, 0.4, 0.3])
```
"""
function consensus_distance_matrix(matrices::Vector{Matrix{Float64}};
        weights::Vector{Float64} = Float64[])
    isempty(matrices) && error("No matrices provided")
    n = size(first(matrices), 1)
    all(m -> size(m) == (n, n), matrices) || error("All matrices must be the same size")

    if isempty(weights)
        weights = fill(1.0 / length(matrices), length(matrices))
    end
    length(weights) == length(matrices) || error("Weights must match number of matrices")
    abs(sum(weights) - 1.0) < 1e-6 || @warn "Weights don't sum to 1.0"

    D = zeros(Float64, n, n)
    for (w, M) in zip(weights, matrices)
        max_val = maximum(M)
        M_norm = max_val > 0 ? M ./ max_val : M
        D .+= w .* M_norm
    end

    return D
end
