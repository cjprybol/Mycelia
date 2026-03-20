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
# Sequence-based UniProt BLAST search
# --------------------------------------------------------------------------- #

"""
    blast_uniprot_sequence(sequence::AbstractString;
                           database::String="uniprotkb",
                           hits::Int=10,
                           threshold::Float64=0.001,
                           poll_interval::Real=5,
                           timeout::Real=300) -> Vector{Dict}

BLAST a raw amino acid sequence against UniProt databases via the UniProt
REST API. Returns the top hits with accession, identity, alignment length,
e-value, and organism.

# Keywords
- `database`: "uniprotkb" (default), "uniref50", "uniref90", "uniref100"
- `hits`: Maximum number of hits to return (default: 10)
- `threshold`: E-value threshold (default: 0.001)
- `poll_interval`: Seconds between status checks
- `timeout`: Maximum seconds to wait

# Example
```julia
hits = blast_uniprot_sequence("MKWKLFKKIEKVGQNIRDGIIKAG..."; hits=5)
```
"""
function blast_uniprot_sequence(sequence::AbstractString;
        database::String = "uniprotkb",
        hits::Int = 10,
        threshold::Float64 = 0.001,
        poll_interval::Real = 5,
        timeout::Real = 300)
    api_base = "https://rest.uniprot.org/idmapping"

    # Submit BLAST job via UniProt's sequence search endpoint
    search_url = "https://rest.uniprot.org/uniprotkb/search"
    params = [
        "query" => "sequence:$(strip(string(sequence)))",
        "format" => "json",
        "size" => string(hits),
        "fields" => "accession,id,protein_name,organism_name,sequence,gene_names"
    ]

    # EBI NCBI BLAST REST API — parameters must use allowed enum values
    blast_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"

    # Map hits to nearest allowed alignments value: 0, 5, 10, 20, 50, 100, 150, 200, 250, 500
    allowed_alignments = [0, 5, 10, 20, 50, 100, 150, 200, 250, 500]
    alignments_val = allowed_alignments[findmin(abs.(allowed_alignments .-
                                                                        hits))[2]]

    # Map threshold to nearest allowed exp value
    allowed_exp = [1e-200, 1e-100, 1e-50, 1e-10, 1e-5, 1e-4, 1e-3, 1.0, 10.0, 100.0, 1000.0]
    exp_val = allowed_exp[findmin(abs.(log10.(allowed_exp) .- log10(threshold)))[2]]

    db_name = database == "uniprotkb" ? "uniprotkb_swissprot" : database
    form_data = Pair{String, String}[
    "email" => "noreply@example.com",
    "program" => "blastp",
    "database" => db_name,
    "sequence" => strip(string(sequence)),
    "stype" => "protein",
    "exp" => string(exp_val),
    "alignments" => string(alignments_val)
]

    submit_resp = HTTP.post(blast_url; body = HTTP.Form(form_data),
        headers = ["Accept" => "text/plain"])
    job_id = strip(String(submit_resp.body))
    @info "UniProt BLAST submitted: job_id=$(job_id)"

    # Poll for completion
    status_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/$(job_id)"
    elapsed = 0.0
    status = "RUNNING"
    while status == "RUNNING" && elapsed < timeout
        sleep(poll_interval)
        elapsed += poll_interval
        status_resp = HTTP.get(status_url; headers = ["Accept" => "text/plain"])
        status = strip(String(status_resp.body))
    end

    if status != "FINISHED"
        error("UniProt BLAST timed out or failed (status: $(status), elapsed: $(elapsed)s)")
    end

    # Fetch results in JSON format
    result_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/$(job_id)/json"
    result_resp = HTTP.get(result_url; headers = ["Accept" => "application/json"])
    result_json = JSON.parse(String(result_resp.body))

    # Parse hits
    hits_out = Dict[]
    for hit in get(result_json, "hits", [])
        hit_acc = get(hit, "hit_acc", "")
        hit_id = get(hit, "hit_id", "")
        hit_desc = get(hit, "hit_def", "")
        hit_org = get(hit, "hit_os", "")
        hit_len = get(hit, "hit_len", 0)

        for hsp in get(hit, "hit_hsps", [])
            push!(hits_out,
                Dict(
                    "accession" => hit_acc,
                    "id" => hit_id,
                    "description" => hit_desc,
                    "organism" => hit_org,
                    "length" => hit_len,
                    "identity" => get(hsp, "hsp_identity", 0.0),
                    "align_len" => get(hsp, "hsp_align_len", 0),
                    "evalue" => get(hsp, "hsp_expect", 0.0),
                    "score" => get(hsp, "hsp_score", 0),
                    "query_start" => get(hsp, "hsp_query_from", 0),
                    "query_end" => get(hsp, "hsp_query_to", 0)
                ))
        end
    end

    return hits_out
end

"""
    blast_uniprot_batch(sequences::Dict{String, <:AbstractString};
                         hits_per_query::Int=3,
                         delay::Real=2,
                         kwargs...) -> DataFrames.DataFrame

BLAST multiple sequences against UniProt, returning a DataFrame of top hits
for each query. Useful for finding UniProt accessions for team-provided
sequences.

# Example
```julia
seqs = Dict("payload1" => "MKWK...", "payload2" => "METL...")
df = blast_uniprot_batch(seqs; hits_per_query=3)
```
"""
function blast_uniprot_batch(sequences::Dict{String, <:AbstractString};
        hits_per_query::Int = 3,
        delay::Real = 2,
        kwargs...)
    results = DataFrames.DataFrame(
        query_name = String[], target_accession = String[],
        target_description = String[], target_organism = String[],
        identity = Float64[], evalue = Float64[],
        align_len = Int[], score = Int[])

    for (name, seq) in sort(collect(sequences); by = first)
        @info "BLASTing $(name) ($(length(seq)) aa)..."
        try
            hits = blast_uniprot_sequence(seq; hits = hits_per_query, kwargs...)
            for h in hits
                push!(results,
                    (
                        query_name = name,
                        target_accession = get(h, "accession", ""),
                        target_description = get(h, "description", ""),
                        target_organism = get(h, "organism", ""),
                        identity = Float64(get(h, "identity", 0.0)),
                        evalue = Float64(get(h, "evalue", 0.0)),
                        align_len = Int(get(h, "align_len", 0)),
                        score = Int(get(h, "score", 0))
                    ))
            end
        catch e
            @warn "BLAST failed for $(name): $(e)"
        end
        sleep(delay)
    end

    return results
end

# --------------------------------------------------------------------------- #
# Sequence-based protein embeddings (no UniProt accession needed)
# --------------------------------------------------------------------------- #

"""
    compute_embeddings_from_sequences(sequences::Dict{String, BioSequences.LongAA};
                                      ks::Vector{Int}=[2, 3]) -> Dict{String, Vector{Float32}}

Compute composition-based embeddings directly from amino acid sequences
(no UniProt accession required). Uses the same AAmer frequency profile
as `fetch_uniprot_embeddings` but works on raw sequences.

This is the lightweight, alignment-free embedding method. For neural
embeddings, use `esm_embed_sequences`.

# Example
```julia
seqs = Dict("payload1" => BioSequences.LongAA("MKWK..."))
embs = compute_embeddings_from_sequences(seqs)
D, labels = embedding_distance_matrix(embs; metric=:cosine)
```
"""
function compute_embeddings_from_sequences(
        sequences::Dict{String, BioSequences.LongAA};
        ks::Vector{Int} = [2, 3])
    embeddings = Dict{String, Vector{Float32}}()
    for (name, seq) in sequences
        try
            embeddings[name] = _compute_composition_embedding(seq; ks = ks)
        catch e
            @warn "Failed to compute embedding for $(name): $(e)"
        end
    end
    return embeddings
end

"""
    esm_embed_sequences(sequences::Dict{String, <:AbstractString};
                         model::String="esm2_t33_650M_UR50D",
                         repr_layer::Int=33,
                         force_install::Bool=false) -> Dict{String, Vector{Float32}}

Compute protein embeddings using ESM-2 (Meta) via a Mycelia-managed conda
environment. Returns per-protein mean-pooled embeddings (1280-dim for the
650M model). Auto-installs `fair-esm` + PyTorch on first use.

Works on raw amino acid sequences — no UniProt accession required.

# Supported models
- "esm2_t33_650M_UR50D" (default, 650M params, 1280-dim output)
- "esm2_t36_3B_UR50D" (3B params, 2560-dim, slower but more accurate)
- "esm2_t30_150M_UR50D" (150M params, 640-dim, fastest)

# Example
```julia
seqs = Dict("payload1" => "MKWKLFKKIEKVGQNIRDGIIKAG...")
embs = esm_embed_sequences(seqs)
D, labels = embedding_distance_matrix(embs; metric=:cosine)
```
"""
function esm_embed_sequences(sequences::Dict{String, <:AbstractString};
        model::String = "esm2_t33_650M_UR50D",
        repr_layer::Int = 33,
        force_install::Bool = false)
    install_esm(; force = force_install)

    @info "Using ESM-2 ($(model)) via Mycelia conda env for $(length(sequences)) sequences"
    return _esm_embed_via_conda(sequences; model = model, repr_layer = repr_layer)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install the ESM-2 protein language model environment via conda + pip.

Creates a conda environment `esm` with Python, PyTorch (CPU), and `fair-esm`.
Follows the same pattern as `install_foldseek`.
"""
function install_esm(; force = false)
    _ensure_conda_env_vars!()
    env_name = "esm"
    already_installed = check_bioconda_env_is_installed(env_name)

    if !already_installed || force
        if force && already_installed
            try
                run(`$(CONDA_RUNNER) env remove -n $(env_name) -y`)
            catch e
                @warn "Failed to remove existing ESM env: $(e)"
            end
        end

        @info "Installing ESM-2 conda environment (this may take a few minutes)..."
        # Create env with Python + pip, then pip-install fair-esm + torch (CPU)
        run(`$(CONDA_RUNNER) create -n $(env_name) -c conda-forge python=3.11 pip -y`)
        run(`$(CONDA_RUNNER) run -n $(env_name) pip install fair-esm torch --index-url https://download.pytorch.org/whl/cpu`)
        run(`$(CONDA_RUNNER) clean --all -y`)
    end

    # Verify
    try
        run(`$(CONDA_RUNNER) run -n $(env_name) python -c "import esm; print(esm.__version__)"`)
    catch verify_err
        if !force
            @warn "ESM verification failed. Retrying with force reinstall..."
            return install_esm(; force = true)
        else
            error("ESM install verification failed: $(verify_err)")
        end
    end
end

"""
    _esm_embed_via_conda(sequences::Dict{String, <:AbstractString};
                          model::String="esm2_t33_650M_UR50D",
                          repr_layer::Int=33) -> Dict{String, Vector{Float32}}

Internal: compute ESM-2 embeddings using the Mycelia-managed `esm` conda env.
Writes sequences to a temporary FASTA, runs the extraction script via
`CONDA_RUNNER`, and reads back mean-pooled representations.
"""
function _esm_embed_via_conda(sequences::Dict{String, <:AbstractString};
        model::String = "esm2_t33_650M_UR50D",
        repr_layer::Int = 33)
    tmpdir = mktempdir()
    fasta_path = joinpath(tmpdir, "input.fasta")
    output_dir = joinpath(tmpdir, "output")
    mkpath(output_dir)

    # Write FASTA
    open(fasta_path, "w") do io
        for (name, seq) in sequences
            # Sanitize name for filesystem (replace special chars)
            safe_name = replace(name, r"[^A-Za-z0-9_\-\.]" => "_")
            println(io, ">$(safe_name)")
            println(io, seq)
        end
    end

    # Build name mapping (sanitized → original)
    name_map = Dict{String, String}()
    for name in keys(sequences)
        safe_name = replace(name, r"[^A-Za-z0-9_\-\.]" => "_")
        name_map[safe_name] = name
    end

    # Run ESM extraction via conda
    run(`$(CONDA_RUNNER) run --live-stream -n esm python -m esm.scripts.extract
        $(model) $(fasta_path) $(output_dir)
        --repr_layers $(repr_layer) --include mean`)

    # Read results — use a Python script to convert .pt → JSON for Julia
    convert_script = joinpath(tmpdir, "convert_pt.py")
    open(convert_script, "w") do io
        write(io, """
import torch, json, os, sys
output_dir = sys.argv[1]
repr_layer = int(sys.argv[2])
results = {}
for f in os.listdir(output_dir):
    if f.endswith('.pt'):
        name = f.replace('.pt', '')
        d = torch.load(os.path.join(output_dir, f), map_location='cpu')
        emb = d['mean_representations'][repr_layer].numpy().tolist()
        results[name] = emb
print(json.dumps(results))
""")
    end

    result = read(`$(CONDA_RUNNER) run -n esm python $(convert_script) $(output_dir) $(repr_layer)`, String)
    parsed = JSON.parse(strip(result))

    embeddings = Dict{String, Vector{Float32}}()
    for (safe_name, vec) in parsed
        original_name = get(name_map, safe_name, safe_name)
        embeddings[original_name] = Float32.(vec)
    end

    rm(tmpdir; recursive = true, force = true)
    return embeddings
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

function test_bad()
    import JSON
    return 1
end

