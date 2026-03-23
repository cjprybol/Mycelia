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

For real pre-trained ProtT5 embeddings (1024-dim, same as UniProt's H5 files),
use `prot5_embed_sequences` instead — it combines H5 lookup with local compute
and produces embeddings in the same vector space as UniProt's pre-computed data.

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

# --------------------------------------------------------------------------- #
# ProtT5 pre-computed embeddings (UniProt FTP H5 files)
# --------------------------------------------------------------------------- #

"""
    download_uniprot_embeddings(;
        organism::String="sprot",
        output_dir::String="\$(homedir())/workspace/UniProt-embeddings",
        force::Bool=false) -> String

Download pre-computed ProtT5-XL-U50 per-protein embeddings (1024-dim) from the
UniProt FTP archive. Returns path to the local H5 file.

# Organisms
- `"sprot"` — all Swiss-Prot (~570K proteins, 1.3 GB)
- `"UP000000625_83333"` — E. coli K-12 proteome
- `"UP000005640_9606"` — Human proteome

# Example
```julia
h5_path = download_uniprot_embeddings(organism="sprot")
```
"""
function download_uniprot_embeddings(;
        organism::String = "sprot",
        output_dir::String = "$(homedir())/workspace/UniProt-embeddings",
        force::Bool = false)
    subpath = organism == "sprot" ? "uniprot_sprot" : organism
    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/$(subpath)/per-protein.h5"

    dest_dir = joinpath(output_dir, subpath)
    mkpath(dest_dir)
    output_path = joinpath(dest_dir, "per-protein.h5")

    if isfile(output_path) && filesize(output_path) > 0 && !force
        @info "ProtT5 embeddings already downloaded: $(output_path) ($(round(filesize(output_path) / 1e9; digits=2)) GB)"
        return output_path
    end

    @info "Downloading ProtT5 embeddings from UniProt FTP..."
    @info "  URL: $(url)"
    @info "  This may take several minutes (Swiss-Prot is ~1.3 GB)"

    Downloads.download(url, output_path)
    @info "Downloaded: $(output_path) ($(round(filesize(output_path) / 1e9; digits=2)) GB)"
    return output_path
end

"""
    lookup_uniprot_embeddings(
        accessions::AbstractVector{<:AbstractString};
        h5_path::String="...",
        auto_download::Bool=true) -> Dict{String, Vector{Float32}}

Look up pre-computed ProtT5-XL-U50 embeddings (1024-dim) for UniProt accessions
from a downloaded H5 file. Accessions not found in the file are skipped with a
warning.

# Example
```julia
embs = lookup_uniprot_embeddings(["P05820", "Q840G9", "P00639"])
D, labels = embedding_distance_matrix(embs; metric=:cosine)
```
"""
function lookup_uniprot_embeddings(
        accessions::AbstractVector{<:AbstractString};
        h5_path::String = joinpath(homedir(), "workspace", "UniProt-embeddings",
            "uniprot_sprot", "per-protein.h5"),
        auto_download::Bool = true)
    if !isfile(h5_path)
        if auto_download
            @info "H5 file not found, downloading..."
            download_uniprot_embeddings()
        else
            error("H5 file not found: $(h5_path)")
        end
    end

    embeddings = Dict{String, Vector{Float32}}()
    found = 0
    missed = 0

    HDF5.h5open(h5_path, "r") do fid
        available_keys = Set(keys(fid))
        for acc in accessions
            # Try bare accession, then with version suffix
            key = if acc in available_keys
                acc
            elseif "$(acc).1" in available_keys
                "$(acc).1"
            else
                nothing
            end

            if key !== nothing
                embeddings[acc] = Float32.(HDF5.read(fid[key]))
                found += 1
            else
                missed += 1
            end
        end
    end

    @info "ProtT5 H5 lookup: $(found)/$(length(accessions)) found, $(missed) missed"
    return embeddings
end

# --------------------------------------------------------------------------- #
# Composition-based embeddings (lightweight, no external tools)
# --------------------------------------------------------------------------- #

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

    # EBI API requires exact enum values — sourced from parameterdetails endpoint
    allowed_alignments = [0, 5, 10, 20, 50, 100, 150, 200, 250, 500, 750, 1000]
    alignments_val = allowed_alignments[findmin(abs.(allowed_alignments .-
                                                                        hits))[2]]

    # Exp values must be exact strings (e.g., "1e-3" not "0.001")
    allowed_exp_strings = [
        "1e-200", "1e-100", "1e-50", "1e-10", "1e-5",
        "1e-4", "1e-3", "1e-2", "1e-1", "1.0", "10", "100", "1000"]
    allowed_exp_vals = [
        1e-200, 1e-100, 1e-50, 1e-10, 1e-5,
        1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]
    exp_idx = findmin(abs.(log10.(allowed_exp_vals) .- log10(threshold)))[2]
    exp_str = allowed_exp_strings[exp_idx]

    # Valid databases include: uniprotkb, uniprotkb_swissprot, uniprotkb_bacteria,
    # uniref50, uniref90, uniref100, pdb, afdb — see EBI parameterdetails/database
    db_name = database
    form_data = Pair{String, String}[
    "email" => "noreply@example.com",
    "program" => "blastp",
    "database" => db_name,
    "sequence" => strip(string(sequence)),
    "stype" => "protein",
    "exp" => exp_str,
    "alignments" => string(alignments_val)
]

    submit_resp = HTTP.post(blast_url; body = HTTP.Form(form_data),
        headers = ["Accept" => "text/plain"])
    job_id = strip(String(submit_resp.body))
    @info "UniProt BLAST submitted: job_id=$(job_id)"

    # Poll for completion — EBI statuses: QUEUED → RUNNING → FINISHED (or ERROR/NOT_FOUND)
    status_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/$(job_id)"
    elapsed = 0.0
    status = "QUEUED"
    while status in ("QUEUED", "RUNNING") && elapsed < timeout
        sleep(poll_interval)
        elapsed += poll_interval
        status_resp = HTTP.get(status_url; headers = ["Accept" => "text/plain"])
        status = strip(String(status_resp.body))
        if elapsed > 30 && elapsed % 30 < poll_interval
            @info "  EBI BLAST $(job_id): $(status) ($(round(elapsed; digits=0))s)"
        end
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
# UniRef FASTA download + local BLAST database construction
# --------------------------------------------------------------------------- #

"""
    download_uniref_fasta(cluster_level::Int;
                           output_dir::String="\$(homedir())/workspace/UniRef",
                           force::Bool=false) -> String

Download UniRef FASTA from UniProt FTP. Returns path to the gzipped FASTA file.

# Arguments
- `cluster_level::Int`: 50, 90, or 100

# Example
```julia
fasta_path = download_uniref_fasta(90)
# => "~/workspace/UniRef/uniref90.fasta.gz"
```
"""
function download_uniref_fasta(cluster_level::Int;
        output_dir::String = "$(homedir())/workspace/UniRef",
        force::Bool = false)
    @assert cluster_level in (50, 90, 100) "cluster_level must be 50, 90, or 100"

    mkpath(output_dir)
    filename = "uniref$(cluster_level).fasta.gz"
    output_path = joinpath(output_dir, filename)

    if isfile(output_path) && filesize(output_path) > 0 && !force
        @info "UniRef$(cluster_level) FASTA already exists: $(output_path) ($(round(filesize(output_path) / 1e9; digits=1)) GB)"
        return output_path
    end

    url = "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref$(cluster_level)/$(filename)"
    @info "Downloading UniRef$(cluster_level) FASTA from $(url)..."
    @info "This is a large file (UniRef50 ~12GB, UniRef90 ~38GB, UniRef100 ~99GB)"

    Downloads.download(url, output_path)
    @info "Downloaded: $(output_path) ($(round(filesize(output_path) / 1e9; digits=1)) GB)"
    return output_path
end

"""
    download_and_build_uniref_blastdb(cluster_level::Int;
                                       uniref_dir::String="\$(homedir())/workspace/UniRef",
                                       blastdb_dir::String="\$(homedir())/workspace/blastdb",
                                       force::Bool=false) -> String

Convenience: download UniRef FASTA + build BLAST DB in one step.
Uses `download_uniref_fasta` + `ensure_blast_db` (from alignments-and-mapping.jl).

# Example
```julia
db_path = download_and_build_uniref_blastdb(90)
run_blastp_search(query_fasta="query.fa", reference_fasta=db_path)
```
"""
function download_and_build_uniref_blastdb(cluster_level::Int;
        uniref_dir::String = "$(homedir())/workspace/UniRef",
        blastdb_dir::String = "$(homedir())/workspace/blastdb",
        force::Bool = false)
    fasta_path = download_uniref_fasta(cluster_level; output_dir = uniref_dir, force = force)
    db_name = "uniref$(cluster_level)"
    ensure_blast_db(;
        fasta = fasta_path,
        dbtype = "prot",
        output_dir = blastdb_dir,
        db_name = db_name,
        title = "UniRef$(cluster_level)",
        parse_seqids = true,
        force = force)
end

# --------------------------------------------------------------------------- #
# NCBI Remote BLAST
#
# For local UniRef searching, use the existing MMSeqs2 functions:
#   run_mmseqs_easy_search(query_fasta="query.fa",
#                          target_database="~/workspace/mmseqs/UniRef50")
#   read_mmseqs_easy_search(results_path)
#
# Or for BLAST-based local search:
#   ensure_blast_db(fasta="uniref50.fasta.gz", dbtype="prot", ...)
#   run_blastp_search(query_fasta="query.fa", reference_fasta=db_path)
# --------------------------------------------------------------------------- #

"""
    blast_ncbi_sequence(sequence::AbstractString;
                         program::String="blastp",
                         database::String="nr",
                         hits::Int=10,
                         evalue::Float64=1e-3,
                         poll_interval::Real=15,
                         timeout::Real=600) -> Vector{Dict}

BLAST a raw sequence against NCBI databases via the NCBI BLAST REST API.
No local BLAST database required. Returns top hits with accession, identity,
e-value, and description.

NCBI rate limits: max 1 request per 10 seconds, max 1 poll per minute per RID.
High-volume submissions (50+) should be scheduled for nights/weekends.

# Databases
- `"nr"` — non-redundant protein sequences (default, comprehensive)
- `"swissprot"` — SwissProt curated subset
- `"refseq_protein"` — NCBI RefSeq proteins
- `"pdb"` — Protein Data Bank sequences

# Example
```julia
hits = blast_ncbi_sequence("MKWKLFKKIEKVGQNIRDGIIKAG..."; database="swissprot")
```
"""
function blast_ncbi_sequence(sequence::AbstractString;
        program::String = "blastp",
        database::String = "nr",
        hits::Int = 10,
        evalue::Float64 = 1e-3,
        poll_interval::Real = 15,
        timeout::Real = 600)
    api_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

    # Submit job
    submit_params = [
        "CMD" => "Put",
        "PROGRAM" => program,
        "DATABASE" => database,
        "QUERY" => strip(string(sequence)),
        "EXPECT" => string(evalue),
        "HITLIST_SIZE" => string(hits),
        "FORMAT_TYPE" => "JSON2",
        "TOOL" => "Mycelia.jl",
        "EMAIL" => "noreply@example.com"
    ]

    submit_resp = HTTP.post(api_base; body = HTTP.Form(submit_params))
    body = String(submit_resp.body)

    # Extract RID from response
    rid_match = match(r"RID = ([A-Z0-9\-]+)", body)
    rid_match === nothing && error("Failed to extract RID from NCBI BLAST response")
    rid = rid_match[1]

    # Extract estimated wait time
    rtoe_match = match(r"RTOE = (\d+)", body)
    wait_time = rtoe_match !== nothing ? parse(Int, rtoe_match[1]) : 30
    @info "NCBI BLAST submitted: RID=$(rid), estimated wait=$(wait_time)s"

    # Wait the estimated time before first poll
    sleep(min(wait_time, poll_interval))

    # Poll for completion (max 1 poll per minute per NCBI guidelines)
    elapsed = Float64(min(wait_time, poll_interval))
    status = "WAITING"
    while status in ("WAITING", "UNKNOWN") && elapsed < timeout
        sleep(max(poll_interval, 60.0))
        elapsed += max(poll_interval, 60.0)

        check_url = "$(api_base)?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$(rid)"
        check_resp = HTTP.get(check_url)
        check_body = String(check_resp.body)

        if occursin("Status=READY", check_body)
            status = "READY"
        elseif occursin("Status=FAILED", check_body)
            error("NCBI BLAST failed for RID=$(rid)")
        elseif occursin("Status=UNKNOWN", check_body)
            error("NCBI BLAST RID expired or invalid: $(rid)")
        end
    end

    if status != "READY"
        error("NCBI BLAST timed out (RID=$(rid), elapsed=$(elapsed)s)")
    end

    # Fetch results in JSON format
    result_url = "$(api_base)?CMD=Get&FORMAT_TYPE=JSON2&RID=$(rid)"
    result_resp = HTTP.get(result_url)
    result_json = JSON.parse(String(result_resp.body))

    # Parse hits from NCBI JSON2 format
    hits_out = Dict[]
    search_results = get(result_json, "BlastOutput2", [])
    for report in search_results
        search = get(get(report, "report", Dict()), "results", Dict())
        for hit in get(search, "search", Dict()) |> d -> get(d, "hits", [])
            hit_desc = get(hit, "description", [])
            desc_entry = isempty(hit_desc) ? Dict() : hit_desc[1]
            accession = get(desc_entry, "accession", "")
            title = get(desc_entry, "title", "")
            taxid = get(desc_entry, "taxid", 0)

            for hsp in get(hit, "hsps", [])
                push!(hits_out,
                    Dict(
                        "accession" => accession,
                        "description" => title,
                        "taxid" => taxid,
                        "identity" => get(hsp, "identity", 0),
                        "align_len" => get(hsp, "align_len", 0),
                        "evalue" => get(hsp, "evalue", 0.0),
                        "bit_score" => get(hsp, "bit_score", 0.0),
                        "query_from" => get(hsp, "query_from", 0),
                        "query_to" => get(hsp, "query_to", 0),
                        "hit_from" => get(hsp, "hit_from", 0),
                        "hit_to" => get(hsp, "hit_to", 0),
                        "gaps" => get(hsp, "gaps", 0)
                    ))
            end
        end
    end

    return hits_out
end

"""
    blast_ncbi_batch(sequences::Dict{String, <:AbstractString};
                      hits_per_query::Int=5,
                      delay::Real=15,
                      kwargs...) -> DataFrames.DataFrame

BLAST multiple sequences against NCBI, returning a DataFrame of top hits.
Respects NCBI rate limits (min 10s between submissions).

# Example
```julia
seqs = Dict("payload1" => "MKWK...", "payload2" => "METL...")
df = blast_ncbi_batch(seqs; database="swissprot", hits_per_query=5)
```
"""
function blast_ncbi_batch(sequences::Dict{String, <:AbstractString};
        hits_per_query::Int = 5,
        delay::Real = 15,
        kwargs...)
    results = DataFrames.DataFrame(
        query_name = String[], target_accession = String[],
        target_description = String[], taxid = Int[],
        identity = Int[], evalue = Float64[],
        bit_score = Float64[], align_len = Int[])

    for (name, seq) in sort(collect(sequences); by = first)
        @info "NCBI BLASTing $(name) ($(length(seq)) aa)..."
        try
            hits = blast_ncbi_sequence(seq; hits = hits_per_query, kwargs...)
            for h in hits
                push!(results,
                    (
                        query_name = name,
                        target_accession = get(h, "accession", ""),
                        target_description = get(h, "description", ""),
                        taxid = Int(get(h, "taxid", 0)),
                        identity = Int(get(h, "identity", 0)),
                        evalue = Float64(get(h, "evalue", 0.0)),
                        bit_score = Float64(get(h, "bit_score", 0.0)),
                        align_len = Int(get(h, "align_len", 0))
                    ))
            end
        catch e
            @warn "NCBI BLAST failed for $(name): $(e)"
        end
        sleep(max(delay, 10))  # NCBI requires min 10s between requests
    end

    return results
end

# --------------------------------------------------------------------------- #
# UniProt ID Mapping (cross-database accession conversion)
# --------------------------------------------------------------------------- #

"""
    map_ids_to_uniprot(ids::AbstractVector{<:AbstractString};
                        from::String="RefSeq_Protein",
                        to::String="UniProtKB",
                        poll_interval::Real=5,
                        timeout::Real=120) -> Dict{String, Vector{String}}

Map accessions from one database to UniProt (or between any supported databases)
via the UniProt ID Mapping REST API.

Returns a Dict mapping each input ID to a vector of matched UniProt accessions.
IDs with no match are omitted.

# Common `from` databases
- `"RefSeq_Protein"` — NCBI RefSeq protein accessions (e.g., WP_123456789)
- `"EMBL-GenBank-DDBJ"` — INSDC nucleotide accessions
- `"EMBL-GenBank-DDBJ_CDS"` — INSDC CDS protein accessions
- `"GI_number"` — NCBI GI numbers (legacy)
- `"PDB"` — PDB structure IDs
- `"UniProtKB_AC-ID"` — UniProt accessions (for mapping to other DBs)

# Common `to` databases
- `"UniProtKB"` — UniProt Knowledgebase (default)
- `"UniProtKB-Swiss-Prot"` — Reviewed entries only
- `"UniRef90"`, `"UniRef50"` — UniRef clusters

# Example
```julia
# Map NCBI RefSeq protein accessions to UniProt
mapping = map_ids_to_uniprot(["WP_012345678", "WP_098765432"])
# => Dict("WP_012345678" => ["P12345"], "WP_098765432" => ["Q67890"])
```
"""
function map_ids_to_uniprot(ids::AbstractVector{<:AbstractString};
        from::String = "RefSeq_Protein",
        to::String = "UniProtKB",
        poll_interval::Real = 5,
        timeout::Real = 120)
    api_base = "https://rest.uniprot.org/idmapping"

    # Submit mapping job
    form_data = Pair{String, String}[
    "from" => from,
    "to" => to,
    "ids" => join(
        ids, ",")
]

    submit_resp = HTTP.post("$(api_base)/run";
        body = HTTP.Form(form_data),
        headers = _UNIPROT_HEADERS)
    submit_json = JSON.parse(String(submit_resp.body))
    job_id = get(submit_json, "jobId", "")
    isempty(job_id) && error("UniProt ID mapping failed to return jobId: $(submit_json)")
    @info "UniProt ID mapping submitted: jobId=$(job_id)"

    # Poll for completion
    elapsed = 0.0
    while elapsed < timeout
        sleep(poll_interval)
        elapsed += poll_interval

        status_resp = HTTP.get("$(api_base)/status/$(job_id)";
            headers = _UNIPROT_HEADERS,
            redirect = true)
        status_body = String(status_resp.body)

        # Job complete when we get redirected to results or status shows complete
        if occursin("results", status_body) || status_resp.status == 303
            break
        end

        status_json = try
            JSON.parse(status_body)
        catch
            continue
        end

        job_status = get(status_json, "jobStatus", "")
        if job_status == "FINISHED"
            break
        elseif job_status == "ERROR"
            error("UniProt ID mapping failed: $(status_json)")
        end
    end

    # Fetch results
    result_resp = HTTP.get("$(api_base)/uniprotkb/results/$(job_id)?format=json&size=500";
        headers = _UNIPROT_HEADERS)
    result_json = JSON.parse(String(result_resp.body))

    # Parse into mapping dict
    mapping = Dict{String, Vector{String}}()
    for result in get(result_json, "results", [])
        from_id = get(result, "from", "")
        to_entry = get(result, "to", Dict())
        to_acc = if isa(to_entry, Dict)
            get(to_entry, "primaryAccession", "")
        else
            string(to_entry)
        end

        isempty(from_id) && continue
        isempty(to_acc) && continue

        if haskey(mapping, from_id)
            push!(mapping[from_id], to_acc)
        else
            mapping[from_id] = [to_acc]
        end
    end

    @info "ID mapping: $(length(mapping)) / $(length(ids)) IDs mapped"
    return mapping
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
# ProtT5 embeddings (same model as UniProt pre-computed H5 files)
# --------------------------------------------------------------------------- #

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install the ProtT5 protein language model environment via conda + pip.

Creates a conda environment `prot5` with Python, PyTorch (CPU), and
HuggingFace `transformers` + `sentencepiece`. Uses the same ProtT5-XL-U50
model as UniProt's pre-computed embeddings, so locally computed vectors
are directly comparable to H5 lookups.
"""
function install_prot5(; force = false)
    _ensure_conda_env_vars!()
    env_name = "prot5"
    already_installed = check_bioconda_env_is_installed(env_name)

    if !already_installed || force
        if force && already_installed
            try
                run(`$(CONDA_RUNNER) env remove -n $(env_name) -y`)
            catch e
                @warn "Failed to remove existing prot5 env: $(e)"
            end
        end

        @info "Installing ProtT5 conda environment (this may take a few minutes)..."
        run(`$(CONDA_RUNNER) create -n $(env_name) -c conda-forge python=3.11 pip -y`)
        run(`$(CONDA_RUNNER) run -n $(env_name) pip install transformers torch sentencepiece --index-url https://download.pytorch.org/whl/cpu`)
        run(`$(CONDA_RUNNER) clean --all -y`)
    end

    # Verify
    try
        run(`$(CONDA_RUNNER) run -n $(env_name) python -c "from transformers import T5EncoderModel; print('OK')"`)
    catch verify_err
        if !force
            @warn "ProtT5 verification failed. Retrying with force reinstall..."
            return install_prot5(; force = true)
        else
            error("ProtT5 install verification failed: $(verify_err)")
        end
    end
end

"""
    compute_prot5_embeddings(
        sequences::Dict{String, <:AbstractString};
        force_install::Bool=false) -> Dict{String, Vector{Float32}}

Compute ProtT5-XL-U50 embeddings (1024-dim) from raw amino acid sequences
using HuggingFace `transformers` via a Mycelia-managed conda environment.

Uses the same model (`Rostlab/prot_t5_xl_uniref50`) as UniProt's pre-computed
embeddings, so vectors are directly comparable to `lookup_uniprot_embeddings`.

First run downloads model weights (~3 GB to ~/.cache/huggingface/).
Processes sequences one at a time for memory safety on constrained machines.

# Example
```julia
seqs = Dict("DspB" => "MNCCVKGNS...", "ColicinM" => "METLTVHAPS...")
embs = compute_prot5_embeddings(seqs)
# => Dict with 1024-dim Vector{Float32} per sequence
```
"""
function compute_prot5_embeddings(
        sequences::Dict{String, <:AbstractString};
        force_install::Bool = false)
    install_prot5(; force = force_install)

    tmpdir = mktempdir()
    fasta_path = joinpath(tmpdir, "input.fasta")
    script_path = joinpath(tmpdir, "embed_prot5.py")

    # Write FASTA with sanitized names
    name_map = Dict{String, String}()
    open(fasta_path, "w") do io
        for (name, seq) in sequences
            safe_name = replace(name, r"[^A-Za-z0-9_\-\.]" => "_")
            name_map[safe_name] = name
            println(io, ">$(safe_name)")
            println(io, seq)
        end
    end

    # Write Python embedding script
    open(script_path, "w") do io
        write(io, """
import json, sys, re
from transformers import T5Tokenizer, T5EncoderModel
import torch

tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
model.to(torch.device("cpu"))

# Read FASTA
fasta_path = sys.argv[1]
sequences = {}
current_name = None
current_seq = []
with open(fasta_path) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_name:
                sequences[current_name] = "".join(current_seq)
            current_name = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_name:
        sequences[current_name] = "".join(current_seq)

# Compute embeddings one at a time (memory-safe)
results = {}
for name, seq in sequences.items():
    # ProtT5 expects space-separated amino acids
    seq_spaced = " ".join(list(re.sub(r"[UZOB]", "X", seq)))
    ids = tokenizer(seq_spaced, return_tensors="pt", add_special_tokens=True, padding=True)
    with torch.no_grad():
        output = model(input_ids=ids["input_ids"], attention_mask=ids["attention_mask"])
    # Mean pool over sequence length (exclude padding)
    mask = ids["attention_mask"].unsqueeze(-1)
    emb = (output.last_hidden_state * mask).sum(1) / mask.sum(1)
    results[name] = emb[0].tolist()
    sys.stderr.write(f"  Embedded {name}: {len(seq)} aa -> {len(results[name])}-dim\\n")

print(json.dumps(results))
""")
    end

    @info "Computing ProtT5 embeddings for $(length(sequences)) sequences..."
    result = read(`$(CONDA_RUNNER) run --live-stream -n prot5 python $(script_path) $(fasta_path)`, String)
    parsed = JSON.parse(strip(result))

    embeddings = Dict{String, Vector{Float32}}()
    for (safe_name, vec) in parsed
        original_name = get(name_map, safe_name, safe_name)
        embeddings[original_name] = Float32.(vec)
    end

    rm(tmpdir; recursive = true, force = true)
    @info "ProtT5 embeddings computed: $(length(embeddings)) sequences, $(isempty(embeddings) ? 0 : length(first(values(embeddings))))-dim"
    return embeddings
end

"""
    prot5_embed_sequences(
        sequences::Dict{String, <:AbstractString};
        accession_map::Dict{String, String}=Dict{String,String}(),
        h5_path::String="...",
        auto_download::Bool=true,
        force_install::Bool=false) -> Dict{String, Vector{Float32}}

Unified ProtT5 embedding function. All embeddings are 1024-dim ProtT5-XL-U50
vectors in the same space, regardless of source:

1. **Phase A:** Look up pre-computed embeddings from UniProt H5 file for
   sequences with known UniProt accessions (via `accession_map`)
2. **Phase B:** Compute ProtT5 locally for remaining sequences using the
   same model (`Rostlab/prot_t5_xl_uniref50`)

# Arguments
- `sequences`: name -> AA sequence string (all candidates)
- `accession_map`: name -> UniProt accession (for H5 lookup). Entries without
  a mapping are computed locally.

# Example
```julia
seqs = Dict("DspB" => "MNCCVKGNS...", "novel_payload" => "MKWKLFKK...")
acc_map = Dict("DspB" => "Q840G9")  # DspB has a known UniProt accession
embs = prot5_embed_sequences(seqs; accession_map=acc_map)
# DspB → from H5 (instant), novel_payload → computed locally
```
"""
function prot5_embed_sequences(
        sequences::Dict{String, <:AbstractString};
        accession_map::Dict{String, String} = Dict{String, String}(),
        h5_path::String = joinpath(homedir(), "workspace", "UniProt-embeddings",
            "uniprot_sprot", "per-protein.h5"),
        auto_download::Bool = true,
        force_install::Bool = false)
    embeddings = Dict{String, Vector{Float32}}()

    # Phase A: Pre-computed H5 lookup
    if !isempty(accession_map)
        accessions = unique(collect(values(accession_map)))
        @info "Phase A: Looking up $(length(accessions)) accessions in ProtT5 H5..."
        h5_embs = lookup_uniprot_embeddings(accessions;
            h5_path = h5_path, auto_download = auto_download)

        # Map back to original names
        for (name, acc) in accession_map
            if haskey(h5_embs, acc)
                embeddings[name] = h5_embs[acc]
            end
        end
        @info "  Phase A result: $(length(embeddings)) embeddings from H5"
    end

    # Phase B: Compute locally for unmatched sequences
    unmatched = Dict{String, String}()
    for (name, seq) in sequences
        if !haskey(embeddings, name)
            unmatched[name] = seq
        end
    end

    if !isempty(unmatched)
        @info "Phase B: Computing ProtT5 locally for $(length(unmatched)) unmatched sequences..."
        local_embs = compute_prot5_embeddings(unmatched; force_install = force_install)
        merge!(embeddings, local_embs)
    end

    n_from_h5 = length(embeddings) - length(unmatched)
    @info "ProtT5 summary: $(length(embeddings)) total ($(n_from_h5) from H5, $(length(unmatched)) computed locally)"
    return embeddings
end

# --------------------------------------------------------------------------- #
# FoldSeek operations — uses existing Mycelia.foldseek_easy_search()
# See foldseek.jl for the core FoldSeek CLI wrappers.
# This section adds the web API fallback and distance matrix parser.
# --------------------------------------------------------------------------- #

"""
    foldseek_web_search(pdb_path::AbstractString;
                         database::String="afdb50",
                         mode::String="3diaa",
                         poll_interval::Real=5,
                         timeout::Real=300) -> Vector{Dict}

Search a single PDB structure against FoldSeek web server databases.
Returns a vector of hit dicts with target, identity, TM-score, etc.

Uses the FoldSeek search API at https://search.foldseek.com/api/ticket

# Example
```julia
hits = foldseek_web_search("structures/AF-P05820.pdb")
```
"""
function foldseek_web_search(pdb_path::AbstractString;
        database::String = "afdb50",
        mode::String = "3diaa",
        poll_interval::Real = 5,
        timeout::Real = 300)
    pdb_content = read(pdb_path, String)

    # Submit search job
    form = HTTP.Form([
        "q" => HTTP.Multipart(basename(pdb_path), IOBuffer(pdb_content)),
        "mode" => mode,
        "database[]" => database
    ])
    submit_resp = HTTP.post("https://search.foldseek.com/api/ticket"; body = form)
    ticket = JSON.parse(String(submit_resp.body))
    ticket_id = ticket["id"]
    @info "FoldSeek web search submitted: ticket=$(ticket_id)"

    # Poll for completion
    elapsed = 0.0
    while elapsed < timeout
        sleep(poll_interval)
        elapsed += poll_interval
        status_resp = HTTP.get("https://search.foldseek.com/api/ticket/$(ticket_id)")
        status = JSON.parse(String(status_resp.body))
        if get(status, "status", "") == "COMPLETE"
            # Fetch results
            results = get(status, "result", [])
            if !isempty(results) && haskey(first(results), "alignments")
                return first(results)["alignments"]
            end
            return []
        elseif get(status, "status", "") == "ERROR"
            error("FoldSeek web search failed: $(get(status, "error", "unknown"))")
        end
    end
    error("FoldSeek web search timed out after $(timeout)s")
end

"""
    foldseek_web_allvsall(structure_dir::AbstractString,
                           output_path::AbstractString;
                           database::String="afdb50",
                           delay::Real=3) -> String

Run FoldSeek web API for all PDB files in a directory and produce a TSV
results file compatible with `structural_distance_matrix`.

This is the web API fallback when local FoldSeek CLI is not available.
Rate-limited to avoid overwhelming the public server.

# Example
```julia
foldseek_web_allvsall("structures/payloads", "figures/foldseek-web.tsv")
D, labels = structural_distance_matrix("figures/foldseek-web.tsv", accessions)
```
"""
function foldseek_web_allvsall(structure_dir::AbstractString,
        output_path::AbstractString;
        database::String = "afdb50",
        delay::Real = 3)
    pdb_files = filter(f -> endswith(f, ".pdb"), readdir(structure_dir; join = true))
    isempty(pdb_files) && error("No PDB files found in $(structure_dir)")
    @info "FoldSeek web all-vs-all: $(length(pdb_files)) structures"

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        for (i, query_pdb) in enumerate(pdb_files)
            query_name = replace(basename(query_pdb), r"\.pdb$" => "")
            @info "  Searching $(query_name) ($(i)/$(length(pdb_files)))..."
            try
                hits = foldseek_web_search(query_pdb; database = database)
                for hit in hits
                    for aln in get(hit, "alignments", [hit])
                        target = get(aln, "target", "")
                        fident = get(aln, "seqId", 0.0)
                        alnlen = get(aln, "alnLength", 0)
                        evalue = get(aln, "eval", 1.0)
                        bits = get(aln, "score", 0.0)
                        tmscore = get(aln, "prob", 0.0)  # TM-score in prob field
                        println(io, join(
                            [query_name, target, fident, alnlen,
                                evalue, bits, tmscore], "\t"))
                    end
                end
            catch e
                @warn "FoldSeek web search failed for $(query_name): $(e)"
            end
            i < length(pdb_files) && sleep(delay)
        end
    end
    @info "FoldSeek web results written to $(output_path)"
    return output_path
end

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
