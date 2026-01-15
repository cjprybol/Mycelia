# # Tutorial 17: Viroid Sketch Round Trip (Mash + sourmash + Sylph)
#
# This tutorial demonstrates a full, reproducible round-trip workflow:
# download a BLAST database, export it to FASTA, index with sketch tools,
# simulate reads from a known subset, and assert expected matches.
# Heavy steps are gated behind `MYCELIA_RUN_EXTERNAL=true`.

# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/17_viroid_sketch_round_trip.jl", "tutorials/notebooks", execute=false)'
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
import Test

run_external = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
if !run_external
    println("Set MYCELIA_RUN_EXTERNAL=true to download references and run external tools.")
end

base_outdir = "results/viroid_sketch_round_trip"
reference_dir = joinpath(base_outdir, "reference")
mkpath(reference_dir)

reference_db = "ref_viroids_rep_genomes"
reference_fasta = joinpath(reference_dir, "$(reference_db).fna.gz")
reference_table_csv = joinpath(reference_dir, "viroids_reference_table.csv")
reference_genomes_dir = joinpath(reference_dir, "genomes")

# ## Step 1: Download the viroid BLAST database

reference_table = DataFrames.DataFrame()
if run_external
    db_path = Mycelia.download_blast_db(db=reference_db, dbdir=Mycelia.DEFAULT_BLASTDB_PATH)
    Test.@test !isempty(db_path)

    # ## Step 2: Export to FASTA + build reference table
    reference_table = Mycelia.prepare_blast_reference_table(
        blastdb=db_path,
        blastdbs_dir=Mycelia.DEFAULT_BLASTDB_PATH,
        download_if_missing=false,
        reference_fasta=reference_fasta,
        taxonomy_map_out=joinpath(reference_dir, "$(reference_db).seqid2taxid.txt.gz"),
        table_out=reference_table_csv,
        force=false
    )
    Test.@test isfile(reference_fasta)
    Test.@test !isempty(reference_table)
end

# ## Step 3: Build reference indices for mash + sourmash

reference_sourmash_sig = ""
reference_mash_db = ""
reference_id_map = Dict{String,String}()
reference_genomes = String[]

if run_external && !isempty(reference_table)
    reference_id_map = Mycelia.split_fasta_by_record(
        fasta_in=reference_fasta,
        outdir=reference_genomes_dir,
        gzip=false,
        force=false
    )
    reference_genomes = Mycelia.find_fasta_files(reference_genomes_dir)
    Test.@test !isempty(reference_genomes)

    reference_sourmash_dir = joinpath(reference_dir, "sourmash")
    reference_mash_dir = joinpath(reference_dir, "mash")
    mkpath.(String[reference_sourmash_dir, reference_mash_dir])

    sourmash_ref = Mycelia.run_sourmash_sketch(
        input_files=[reference_fasta],
        outdir=reference_sourmash_dir,
        k_sizes=[31],
        scaled=1000,
        singleton=true
    )
    Test.@test !isempty(sourmash_ref.signatures)
    reference_sourmash_sig = sourmash_ref.signatures[1]

    mash_result = Mycelia.run_mash_sketch(
        input_files=reference_genomes,
        outdir=reference_mash_dir,
        k=21,
        s=10000,
        r=false
    )
    Test.@test !isempty(mash_result.sketches)
    reference_mash_db = Mycelia.run_mash_paste(
        out_file=joinpath(reference_mash_dir, "$(reference_db).msh"),
        in_files=mash_result.sketches
    )
end

# ## Step 4: Simulate reads from a known subset

depth_target = 20
n_organisms = 1
readset = :illumina_pe150
replicate = 1
rng = StableRNGs.StableRNG(1234)

sim = nothing
selected_ids = String[]
scenario_dir = joinpath(base_outdir, "depth$(depth_target)_div$(n_organisms)")

if run_external && !isempty(reference_table)
    available_ids = unique(String.(reference_table.sequence_id))
    selected_ids = StatsBase.sample(rng, available_ids, n_organisms; replace=false)
    weights = fill(1.0, n_organisms)

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

# ## Step 5: Screen reads with mash, sourmash, and Sylph

if run_external && !isnothing(sim)
    reads_forward = sim.reads.forward
    reads_reverse = sim.reads.reverse
    Test.@test !isnothing(reads_forward)

    merged_reads = reads_forward
    if !isnothing(reads_reverse)
        merged_reads = joinpath(scenario_dir, "reads_merged.fq.gz")
        Mycelia.concatenate_fastq_files(
            fastq_files=[reads_forward, reads_reverse],
            output_fastq=merged_reads,
            gzip=true,
            force=true
        )
    end

    sourmash_outdir = joinpath(scenario_dir, "sourmash")
    mash_outdir = joinpath(scenario_dir, "mash")
    sylph_outdir = joinpath(scenario_dir, "sylph")
    mkpath.(String[sourmash_outdir, mash_outdir, sylph_outdir])

    query_sig = Mycelia.run_sourmash_sketch(
        input_files=[merged_reads],
        outdir=sourmash_outdir,
        k_sizes=[31],
        scaled=1000
    ).signatures[1]
    sourmash_gather = Mycelia.run_sourmash_gather(
        query_sig=query_sig,
        database_sig=reference_sourmash_sig,
        outdir=sourmash_outdir,
        k_size=31
    )
    sourmash_df = CSV.read(sourmash_gather.results_csv, DataFrames.DataFrame)

    mash_screen = Mycelia.run_mash_screen(
        reference=reference_mash_db,
        query=merged_reads,
        outdir=mash_outdir,
        winner_takes_all=true
    )
    mash_df = Mycelia.parse_mash_screen_output(mash_screen.results_tsv)

    sylph_result = if isnothing(reads_reverse)
        Mycelia.run_sylph_profile(
            reference_fastas=reference_genomes,
            sample_reads=[reads_forward],
            outdir=sylph_outdir
        )
    else
        Mycelia.run_sylph_profile(
            reference_fastas=reference_genomes,
            first_pairs=[reads_forward],
            second_pairs=[reads_reverse],
            outdir=sylph_outdir
        )
    end

    # ## Step 6: Assert expected matches
    #
    # We expect the selected reference(s) to appear in the top hits.

    truth_ids = String.(selected_ids)

    Test.@test !isempty(mash_df)
    mash_preds = [Mycelia.normalize_reference_label(x, reference_id_map) for x in mash_df.reference]
    Test.@test truth_ids[1] in mash_preds

    Test.@test !isempty(sourmash_df)
    name_col = if "name" in DataFrames.names(sourmash_df)
        "name"
    elseif "match_name" in DataFrames.names(sourmash_df)
        "match_name"
    else
        ""
    end
    Test.@test !isempty(name_col)
    sourmash_preds = [Mycelia.normalize_reference_label(x, reference_id_map) for x in sourmash_df[!, name_col]]
    Test.@test truth_ids[1] in sourmash_preds

    Test.@test !isempty(sylph_result.table)
    sylph_col = first(filter(col -> col in DataFrames.names(sylph_result.table), ["genome", "reference", "ref", "name"]), "")
    Test.@test !isempty(sylph_col)
    sylph_preds = [Mycelia.normalize_reference_label(x, reference_id_map) for x in sylph_result.table[!, sylph_col]]
    Test.@test truth_ids[1] in sylph_preds
end
