"""
$(DocStringExtensions.TYPEDSIGNATURES)

Cluster protein or nucleotide sequences using MMseqs2 easy-cluster workflow.

# Arguments
- `fasta::String`: Path to input FASTA file containing sequences to cluster
- `output::String`: Base path for output files (default: input path + ".mmseqs_easy_cluster")
- `tmp::String`: Path to temporary directory (default: auto-generated temp dir)

# Returns
- `String`: Path to the output cluster TSV file containing cluster assignments

# Details
Uses MMseqs2 with minimum sequence identity threshold of 50% (-min-seq-id 0.5) and 
minimum coverage threshold of 80% (-c 0.8). The output TSV file format contains 
tab-separated cluster representative and member sequences.
"""
# --cov-mode: coverage mode (0: coverage of query and target, 1: coverage of target, 2: coverage of query)
function mmseqs_easy_cluster(;fasta, output=fasta*".mmseqs_easy_cluster", tmp=mktempdir())
    outfile = "$(output)_cluster.tsv"
    if !isfile(outfile)
        Mycelia.add_bioconda_env("mmseqs2")
        # at least 50% equivalent
        # --min-seq-id 0.5 -c 0.8
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-cluster $(fasta) $(output) $(tmp)`)
    end
    rm(tmp, recursive=true)
    return "$(output)_cluster.tsv"
end