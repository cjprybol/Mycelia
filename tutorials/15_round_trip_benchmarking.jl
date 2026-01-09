# # Round-Trip Benchmark Tutorial (Viroids)
#
# This tutorial walks through a small, reproducible round-trip benchmark:
# download a viroid reference database, build taxonomy, simulate reads, and
# screen with sketch tools. Heavy steps are gated behind `MYCELIA_RUN_EXTERNAL=true`.

# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/15_round_trip_benchmarking.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import DataFrames
import CSV
import StableRNGs
import StatsBase

run_external = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
if !run_external
    println("Set MYCELIA_RUN_EXTERNAL=true to download references and run external tools.")
end

base_outdir = "results/round_trip_tutorial"
reference_dir = joinpath(base_outdir, "reference")
mkpath(reference_dir)

reference_db = "ref_viroids_rep_genomes"
taxonomy_db = "taxdb"
reference_fasta = joinpath(reference_dir, "$(reference_db).fna.gz")

# ## Step 1: Download reference databases
#
# We use Mycelia's BLAST DB downloader for both the viroid database and taxdb.

reference_table = DataFrames.DataFrame()
if run_external
    Mycelia.download_blast_db(db=taxonomy_db)
    reference_table = Mycelia.prepare_blast_reference_table(
        blastdb=reference_db,
        blastdbs_dir=Mycelia.DEFAULT_BLASTDB_PATH,
        download_if_missing=true,
        reference_fasta=reference_fasta,
        taxonomy_map_out=joinpath(reference_dir, "$(reference_db).seqid2taxid.txt.gz"),
        table_out=joinpath(reference_dir, "viroids_reference_table.csv"),
        force=false
    )
end

# ## Step 3: Define one small scenario

depth_target = 10
n_organisms = 1
balance = :equal
readset = :illumina_pe150
replicate = 1
rng = StableRNGs.StableRNG(1234)

if run_external && !isempty(reference_table)
    available_ids = unique(String.(reference_table.sequence_id))
    selected_ids = StatsBase.sample(rng, available_ids, n_organisms; replace=false)
    weights = Mycelia.sample_abundance_weights(n_organisms=n_organisms, balance=balance, rng=rng)

    scenario_dir = joinpath(base_outdir, "depth$(depth_target)_div$(n_organisms)_$(balance)")
    sim = Mycelia.simulate_metagenome_community(
        reference_fasta=reference_fasta,
        reference_table=reference_table,
        n_organisms=n_organisms,
        depth_target=depth_target,
        abundance_profile=:custom,
        readset=readset,
        outdir=scenario_dir,
        selected_ids=selected_ids,
        weights=weights,
        rng=rng,
        replicate=replicate,
        run_simulation=true,
        emit_truth_reads=true
    )
    CSV.write(joinpath(scenario_dir, "truth_table.csv"), sim.truth_table)
end

# ## Step 4: Screen reads with sketch tools
#
# Use Mash, sourmash, and Sylph to identify supported references.

if run_external
    println("Run sketch screening in the benchmark script:")
    println("julia --project=. benchmarking/15_round_trip_benchmark.jl")
end
