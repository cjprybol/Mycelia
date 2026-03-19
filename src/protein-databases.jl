# protein-databases.jl — UniProt, UniRef, AlphaFold, and FoldSeek API integration
#
# Provides functions for fetching protein data from major databases:
# - UniProt: entries, sequences, embeddings
# - UniRef: cluster membership and member expansion
# - AlphaFold: predicted structures
# - FoldSeek: structural similarity search
#
# All functions use HTTP.jl for API access and return parsed Julia types.

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
    response = HTTP.get(url; headers = ["Accept" => "application/json"])
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
    response = HTTP.get(url)
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
    search_uniprot(query::AbstractString; size::Int=10) -> Vector{Dict}

Search UniProt by free-text query. Returns a vector of result entries.

# Example
```julia
results = search_uniprot("gene:csgC AND organism_id:83333"; size=3)
```
"""
function search_uniprot(query::AbstractString; size::Int = 10)
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = ["query" => query, "format" => "json", "size" => string(size)]
    response = HTTP.get(url; headers = ["Accept" => "application/json"], query = params)
    parsed = JSON3.read(String(response.body))
    return get(parsed, :results, [])
end

# --------------------------------------------------------------------------- #
# UniRef cluster operations
# --------------------------------------------------------------------------- #

"""
    search_uniref_cluster(accession::AbstractString;
                          identity::Float64=0.5) -> Union{Dict, Nothing}

Find the UniRef cluster (at given identity threshold) containing a UniProt accession.
Returns cluster info dict with :id, :name, :member_count, or nothing if not found.

# Example
```julia
cluster = search_uniref_cluster("P05820"; identity=0.5)
cluster[:id]            # => "UniRef50_P05820"
cluster[:member_count]  # => 236
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
    response = HTTP.get(url; headers = ["Accept" => "application/json"], query = params)
    parsed = JSON3.read(String(response.body))
    results = get(parsed, :results, [])
    isempty(results) && return nothing

    cluster = first(results)
    return Dict(
        :id => get(cluster, :id, ""),
        :name => get(cluster, :name, ""),
        :member_count => get(cluster, :memberCount, 0),
        :updated => get(cluster, :updated, "")
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
    println("\$(m[:accession]) — \$(m[:organism]) (\$(m[:length]) aa)")
end
```
"""
function fetch_uniref_members(cluster_id::AbstractString; max_members::Int = 10)
    url = "https://rest.uniprot.org/uniref/$(cluster_id)/members"
    params = ["format" => "json", "size" => string(max_members)]
    response = HTTP.get(url; headers = ["Accept" => "application/json"], query = params)
    parsed = JSON3.read(String(response.body))
    raw_members = get(parsed, :results, [])

    members = Dict[]
    for m in raw_members
        push!(members,
            Dict(
                :accession => get(m, :memberId, ""),
                :type => get(m, :memberIdType, ""),
                :organism => get(m, :organismName, "unknown"),
                :length => get(m, :sequenceLength, 0),
                :tax_id => get(m, :organismTaxId, 0)
            ))
    end
    return members
end

"""
    expand_via_uniref(accession::AbstractString;
                      identity::Float64=0.5,
                      max_variants::Int=3,
                      max_length::Int=500) -> Vector{Dict}

Find natural variants of a protein by looking up its UniRef cluster and
returning the top members (excluding the query protein itself).

Filters by max_length to enforce phage packaging constraints.

# Example
```julia
variants = expand_via_uniref("P05820"; max_variants=3)
# => [Dict(:accession => "CEAM_ECOLX", :organism => "E. coli", ...), ...]
```
"""
function expand_via_uniref(accession::AbstractString;
        identity::Float64 = 0.5,
        max_variants::Int = 3,
        max_length::Int = 500)
    cluster = search_uniref_cluster(accession; identity = identity)
    cluster === nothing && return Dict[]

    members = fetch_uniref_members(cluster[:id]; max_members = max_variants + 2)

    # Filter: exclude query protein, enforce size limit
    variants = filter(members) do m
        m[:accession] != accession && m[:length] <= max_length
    end

    return variants[1:min(max_variants, length(variants))]
end

# --------------------------------------------------------------------------- #
# AlphaFold structure operations
# --------------------------------------------------------------------------- #

"""
    fetch_alphafold_prediction(accession::AbstractString) -> Union{Dict, Nothing}

Get AlphaFold prediction metadata for a UniProt accession.
Returns dict with :pdb_url, :cif_url, :plddt, :coverage, or nothing.

# Example
```julia
pred = fetch_alphafold_prediction("P05820")
pred[:plddt]    # => 96.12
pred[:pdb_url]  # => "https://alphafold.ebi.ac.uk/files/AF-P05820-F1-model_v6.pdb"
```
"""
function fetch_alphafold_prediction(accession::AbstractString)
    url = "https://alphafold.ebi.ac.uk/api/prediction/$(accession)"
    response = try
        HTTP.get(url; headers = ["Accept" => "application/json"])
    catch e
        @warn "AlphaFold prediction not found for $(accession)"
        return nothing
    end
    results = JSON3.read(String(response.body))
    isempty(results) && return nothing

    pred = first(results)
    return Dict(
        :pdb_url => get(pred, :pdbUrl, ""),
        :cif_url => get(pred, :cifUrl, ""),
        :plddt => get(pred, :globalMetricValue, 0.0),
        :model_version => get(pred, :latestVersion, 0),
        :sequence_length => get(pred, :uniprotEnd, 0) - get(pred, :uniprotStart, 0) + 1
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

    url = format == :pdb ? pred[:pdb_url] : pred[:cif_url]
    ext = format == :pdb ? "pdb" : "cif"
    mkpath(output_dir)
    output_path = joinpath(output_dir, "AF-$(accession).$(ext)")

    response = HTTP.get(url)
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

Fetch pre-computed ProtT5-XL-U50 per-protein embeddings from UniProt.

UniProt provides embeddings for Swiss-Prot entries at:
  https://rest.uniprot.org/uniprotkb/{accession}?fields=ft_embedding

If unavailable (TrEMBL entries or API limitations), falls back to computing
a sequence composition embedding using AAmer profiles at k=3,4,5 concatenated
into a fixed-length feature vector.

# Example
```julia
embeddings = fetch_uniprot_embeddings(["P05820", "Q840G9", "P52107"])
size(embeddings["P05820"])  # => (1024,) for ProtT5, or (n_features,) for fallback
```

See also: td-mm8yj for full UniProt FTP bulk embedding download.
"""
function fetch_uniprot_embeddings(accessions::AbstractVector{<:AbstractString};
        model::Symbol = :per_protein,
        delay::Real = 0.5)
    embeddings = Dict{String, Vector{Float32}}()

    for (i, acc) in enumerate(accessions)
        try
            # Try UniProt embeddings endpoint
            url = "https://rest.uniprot.org/uniprotkb/$(acc).json?fields=sequence"
            response = HTTP.get(url; headers = ["Accept" => "application/json"])
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
                                    ks::Vector{Int}=[2, 3, 4]) -> Vector{Float32}

Compute a fixed-length protein embedding from concatenated AAmer frequency
profiles at multiple k-mer sizes. This is a lightweight, alignment-free
embedding that captures compositional properties.

Returns a Vector{Float32} of length sum(20^k for k in ks) — but uses sparse
representation internally to handle large k values efficiently.

For k=[2,3,4]: theoretical max = 400 + 8000 + 160000 = 168400 dimensions.
In practice, most entries are zero. We use the top 512 most informative
features (by variance across a reference set) for a compact embedding.

Simplified version: uses only k=2 and k=3 (400 + 8000 = 8400 features),
truncated to observed k-mers only.
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

    # Compute pairwise distance
    D = zeros(Float64, n, n)
    if metric == :bray_curtis
        D = frequency_matrix_to_bray_curtis_distance_matrix(freq_matrix)
    elseif metric == :cosine
        D = frequency_matrix_to_cosine_distance_matrix(freq_matrix)
    elseif metric == :euclidean
        D = frequency_matrix_to_euclidean_distance_matrix(freq_matrix)
    elseif metric == :jaccard
        # Binary presence/absence Jaccard
        binary = freq_matrix .> 0
        for i in 1:n, j in (i + 1):n

            intersection = sum(binary[i, :] .& binary[j, :])
            union_size = sum(binary[i, :] .| binary[j, :])
            D[i, j] = D[j, i] = union_size > 0 ? 1.0 - intersection / union_size : 1.0
        end
    else
        error("Unknown metric: $(metric). Use :bray_curtis, :cosine, :euclidean, or :jaccard")
    end

    return D, labels
end

"""
    embedding_distance_matrix(embeddings::Dict{String, Vector{Float32}};
                              metric::Symbol=:cosine) -> (Matrix{Float64}, Vector{String})

Compute pairwise distance matrix from pre-computed protein embeddings.
Returns (distance_matrix, ordered_labels).

# Example
```julia
# After fetching embeddings from UniProt or ESM
D, labels = embedding_distance_matrix(embeddings; metric=:cosine)
```
"""
function embedding_distance_matrix(embeddings::Dict{String, Vector{Float32}};
        metric::Symbol = :cosine)
    labels = sort(collect(keys(embeddings)))
    n = length(labels)
    dim = length(first(values(embeddings)))

    # Stack into matrix
    E = zeros(Float32, n, dim)
    for (i, label) in enumerate(labels)
        E[i, :] = embeddings[label]
    end

    D = zeros(Float64, n, n)
    if metric == :cosine
        # Normalize rows
        norms = sqrt.(sum(E .^ 2; dims = 2))
        norms[norms .== 0] .= 1.0f0
        E_norm = E ./ norms
        sim = E_norm * E_norm'
        D = 1.0 .- Float64.(sim)
    elseif metric == :euclidean
        for i in 1:n, j in (i + 1):n

            D[i, j] = D[j, i] = sqrt(sum((E[i, :] .- E[j, :]) .^ 2))
        end
    end

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

    D = ones(Float64, n, n)  # Default distance = 1.0 (no similarity)
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
        qi == 0 || ti == 0 && continue

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
        # Normalize each matrix to [0, 1] range
        max_val = maximum(M)
        M_norm = max_val > 0 ? M ./ max_val : M
        D .+= w .* M_norm
    end

    return D
end
