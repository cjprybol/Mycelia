# # Round-Trip Benchmark: Viroid Metagenomics
#
# End-to-end benchmark covering reference preparation, read simulation,
# sketch-based screening, and coverage validation across multiple scenarios.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import DataFrames
import CSV
import StableRNGs
import StatsBase
import Dates
import JSON

run_external = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
if !run_external
    println("MYCELIA_RUN_EXTERNAL=false; external downloads/tools will be skipped.")
end

base_outdir = get(ENV, "MYCELIA_ROUND_TRIP_OUTDIR", "results/round_trip_benchmark")
reference_dir = joinpath(base_outdir, "reference")
mkpath(reference_dir)

reference_db = "ref_viroids_rep_genomes"
reference_fasta = joinpath(reference_dir, "$(reference_db).fna.gz")
reference_table_csv = joinpath(reference_dir, "viroids_reference_table.csv")
reference_table_arrow = joinpath(reference_dir, "viroids_reference_table.arrow")
reference_genomes_dir = joinpath(reference_dir, "genomes")
reference_id_map_json = joinpath(reference_dir, "reference_id_map.json")

reference_table = if !isfile(reference_table_csv) || !isfile(reference_fasta)
    run_external || error("Reference data missing. Set MYCELIA_RUN_EXTERNAL=true to download and build.")
    Mycelia.prepare_blast_reference_table(
        blastdb=reference_db,
        blastdbs_dir=Mycelia.DEFAULT_BLASTDB_PATH,
        download_if_missing=true,
        reference_fasta=reference_fasta,
        taxonomy_map_out=joinpath(reference_dir, "$(reference_db).seqid2taxid.txt.gz"),
        table_out=reference_table_csv,
        arrow_out=reference_table_arrow,
        force=false
    )
else
    CSV.read(reference_table_csv, DataFrames.DataFrame)
end

reference_id_map = if isfile(reference_id_map_json)
    Dict{String,String}((k, String(v)) for (k, v) in JSON.parse(String(read(reference_id_map_json))))
else
    id_map = Mycelia.split_fasta_by_record(fasta_in=reference_fasta, outdir=reference_genomes_dir, gzip=false, force=false)
    open(reference_id_map_json, "w") do io
        JSON.print(io, id_map; 4)
    end
    id_map
end

reference_genomes = Mycelia.find_fasta_files(reference_genomes_dir)
@assert !isempty(reference_genomes) "No reference genomes found in $(reference_genomes_dir)"

depths = [10, 100, 1000]
diversities = Dict("low" => 1, "medium" => 10, "high" => 100)
balances = [:equal, :random, :log_normal]
readsets = [:illumina_pe150, :illumina_se250, :nanopore, :pacbio_hifi]
n_replicates = parse(Int, get(ENV, "MYCELIA_ROUND_TRIP_REPLICATES", "1"))
base_seed = parse(Int, get(ENV, "MYCELIA_ROUND_TRIP_SEED", "1337"))

available_ids = unique(String.(reference_table.sequence_id))

reference_sketch_dir = joinpath(reference_dir, "sketches")
reference_sourmash_dir = joinpath(reference_sketch_dir, "sourmash")
reference_mash_dir = joinpath(reference_sketch_dir, "mash")

reference_sourmash_sig = joinpath(reference_sourmash_dir, "$(reference_db).sig")
if run_external
    sourmash_ref = Mycelia.run_sourmash_sketch(input_files=[reference_fasta], outdir=reference_sourmash_dir, k_sizes=[21, 31], scaled=1000, singleton=true)
    if !isempty(sourmash_ref.signatures)
        reference_sourmash_sig = sourmash_ref.signatures[1]
    end
end

reference_mash_db = joinpath(reference_mash_dir, "$(reference_db).msh")
if run_external
    mash_result = Mycelia.run_mash_sketch(input_files=reference_genomes, outdir=reference_mash_dir, k=21, s=10000, r=false)
    Mycelia.run_mash_paste(out_file=reference_mash_db, in_files=mash_result.sketches)
end

scenario_index = 0
for depth_target in depths
    for (div_label, n_organisms) in diversities
        for balance in balances
            for rep in 1:n_replicates
                scenario_index += 1
                rng = StableRNGs.StableRNG(base_seed + scenario_index)
                selected_ids = StatsBase.sample(rng, available_ids, n_organisms; replace=false)
                weights = Mycelia.sample_abundance_weights(n_organisms=n_organisms, balance=balance, rng=rng)

                scenario_id = "depth$(depth_target)_div$(div_label)_bal$(balance)_rep$(rep)"
                scenario_dir = joinpath(base_outdir, scenario_id)
                mkpath(scenario_dir)
                metadata_path = joinpath(scenario_dir, "scenario_metadata.json")
                if !isfile(metadata_path)
                    open(metadata_path, "w") do io
                        JSON.print(io, Dict(
                            "depth_target" => depth_target,
                            "diversity" => n_organisms,
                            "balance" => string(balance),
                            "replicate" => rep,
                            "selected_ids" => selected_ids,
                            "weights" => weights
                        ); 4)
                    end
                end

                truth_table_path = joinpath(scenario_dir, "truth_table.csv")
                for readset in readsets
                    readset_dir = joinpath(scenario_dir, string(readset))
                    sim = Mycelia.simulate_metagenome_community(
                        reference_fasta=reference_fasta,
                        reference_table=reference_table,
                        n_organisms=n_organisms,
                        depth_target=depth_target,
                        abundance_profile=:custom,
                        readset=readset,
                        outdir=readset_dir,
                        selected_ids=selected_ids,
                        weights=weights,
                        rng=rng,
                        replicate=rep,
                        run_simulation=run_external,
                        emit_truth_reads=true,
                        quiet=true,
                        cleanup=true
                    )
                    if !isfile(truth_table_path)
                        CSV.write(truth_table_path, sim.truth_table)
                    end
                    if !run_external
                        continue
                    end

                    reads_forward = sim.reads.forward
                    reads_reverse = sim.reads.reverse
                    @assert !isnothing(reads_forward) "Missing reads for $(scenario_id) $(readset)"

                    merged_reads = reads_forward
                    if !isnothing(reads_reverse)
                        merged_reads = joinpath(readset_dir, "reads_merged.fq.gz")
                        Mycelia.concatenate_fastq_files(fastq_files=[reads_forward, reads_reverse], output_fastq=merged_reads, gzip=true, force=true)
                    end

                    sourmash_outdir = joinpath(readset_dir, "sourmash")
                    mash_outdir = joinpath(readset_dir, "mash")
                    sylph_outdir = joinpath(readset_dir, "sylph")
                    mkpath.(String[sourmash_outdir, mash_outdir, sylph_outdir])

                    query_sig = Mycelia.run_sourmash_sketch(input_files=[merged_reads], outdir=sourmash_outdir, k_sizes=[31], scaled=1000).signatures[1]
                    sourmash_gather = Mycelia.run_sourmash_gather(query_sig=query_sig, database_sig=reference_sourmash_sig, outdir=sourmash_outdir, k_size=31)
                    sourmash_df = CSV.read(sourmash_gather.results_csv, DataFrames.DataFrame)

                    mash_screen = Mycelia.run_mash_screen(reference=reference_mash_db, query=[merged_reads], outdir=mash_outdir, winner_takes_all=true)
                    mash_df = Mycelia.parse_mash_screen_output(mash_screen.results_tsv)

                    sylph_result = if isnothing(reads_reverse)
                        Mycelia.run_sylph_profile(reference_fastas=reference_genomes, sample_reads=[reads_forward], outdir=sylph_outdir)
                    else
                        Mycelia.run_sylph_profile(reference_fastas=reference_genomes, first_pairs=[reads_forward], second_pairs=[reads_reverse], outdir=sylph_outdir)
                    end

                    truth_ids = String.(selected_ids)
                    predicted_accessions = Dict{String,Vector{String}}()
                    if !isempty(mash_df)
                        predicted_accessions["mash"] = [Mycelia.normalize_reference_label(x, reference_id_map) for x in mash_df.reference]
                    else
                        predicted_accessions["mash"] = String[]
                    end
                    if !isempty(sourmash_df)
                        name_col = if "name" in DataFrames.names(sourmash_df)
                            "name"
                        elseif "match_name" in DataFrames.names(sourmash_df)
                            "match_name"
                        else
                            ""
                        end
                        predicted_accessions["sourmash"] = isempty(name_col) ? String[] : [Mycelia.normalize_reference_label(x, reference_id_map) for x in sourmash_df[!, name_col]]
                    else
                        predicted_accessions["sourmash"] = String[]
                    end
                    if !isempty(sylph_result.table)
                        sylph_col = first(filter(col -> col in DataFrames.names(sylph_result.table), ["genome", "reference", "ref", "name"]), "")
                        predicted_accessions["sylph"] = isempty(sylph_col) ? String[] : [Mycelia.normalize_reference_label(x, reference_id_map) for x in sylph_result.table[!, sylph_col]]
                    else
                        predicted_accessions["sylph"] = String[]
                    end

                    metrics_rows = DataFrames.DataFrame(tool=String[], precision=Float64[], recall=Float64[], f1=Float64[], tp=Int[], fp=Int[], fn=Int[])
                    taxonomy_rows = DataFrames.DataFrame(tool=String[], rank=String[], precision=Float64[], recall=Float64[], f1=Float64[], tp=Int[], fp=Int[], fn=Int[])
                    for (tool, preds) in predicted_accessions
                        scores = Mycelia.presence_precision_recall_f1(truth_ids, preds)
                        DataFrames.push!(metrics_rows, (tool=tool, precision=scores.precision, recall=scores.recall, f1=scores.f1, tp=scores.tp, fp=scores.fp, fn=scores.fn))
                        tax_metrics = Mycelia.evaluate_taxonomy_presence_metrics(reference_table, truth_ids, preds)
                        tax_metrics[!, :tool] .= tool
                        taxonomy_rows = DataFrames.vcat(taxonomy_rows, tax_metrics)
                    end
                    CSV.write(joinpath(readset_dir, "metrics_accession.csv"), metrics_rows)
                    CSV.write(joinpath(readset_dir, "metrics_taxonomy.csv"), taxonomy_rows)

                    subset_fasta = joinpath(readset_dir, "reference_subset.fna.gz")
                    Mycelia.subset_fasta_by_ids(fasta_in=reference_fasta, ids=truth_ids, fasta_out=subset_fasta, force=true)

                    map_tmpdir = joinpath(readset_dir, "mapping_tmp")
                    mapping_type = readset == :nanopore ? "map-ont" : (readset == :pacbio_hifi ? "map-pb" : "sr")
                    if isnothing(reads_reverse)
                        map_res = Mycelia.minimap_merge_map_and_split(
                            reference_fasta=subset_fasta,
                            mapping_type=mapping_type,
                            single_end_fastqs=[reads_forward],
                            fastq_mode=:joint,
                            run_splitting=false,
                            write_read_map=false,
                            tmpdir=map_tmpdir
                        )
                    else
                        map_res = Mycelia.minimap_merge_map_and_split(
                            reference_fasta=subset_fasta,
                            mapping_type="sr",
                            paired_end_fastqs=[(reads_forward, reads_reverse)],
                            fastq_mode=:joint,
                            run_splitting=false,
                            write_read_map=false,
                            tmpdir=map_tmpdir
                        )
                    end
                    bam_path = map_res.merged_bam
                    Mycelia.run_mosdepth(bam_path; no_per_base=true)
                    coverm_df = Mycelia.run_coverm_contig(bam_files=[bam_path], reference_fasta=subset_fasta, outdir=joinpath(readset_dir, "coverm"), methods=["mean", "covered_fraction"])
                    CSV.write(joinpath(readset_dir, "coverm_contig.tsv"), coverm_df)

                    if !isnothing(sim.truth_reads.forward) && (readset == :illumina_pe150 || readset == :illumina_se250)
                        truth_forward = sim.truth_reads.forward
                        truth_reverse = sim.truth_reads.reverse
                        truth_tmpdir = joinpath(readset_dir, "truth_mapping_tmp")
                        if isnothing(truth_reverse)
                            truth_map = Mycelia.minimap_merge_map_and_split(
                                reference_fasta=subset_fasta,
                                mapping_type="sr",
                                single_end_fastqs=[truth_forward],
                                fastq_mode=:joint,
                                run_splitting=false,
                                write_read_map=false,
                                tmpdir=truth_tmpdir
                            )
                        else
                            truth_map = Mycelia.minimap_merge_map_and_split(
                                reference_fasta=subset_fasta,
                                mapping_type="sr",
                                paired_end_fastqs=[(truth_forward, truth_reverse)],
                                fastq_mode=:joint,
                                run_splitting=false,
                                write_read_map=false,
                                tmpdir=truth_tmpdir
                            )
                        end
                        truth_bam = truth_map.merged_bam
                        truth_coverm = Mycelia.run_coverm_contig(bam_files=[truth_bam], reference_fasta=subset_fasta, outdir=joinpath(readset_dir, "coverm_truth"), methods=["mean", "covered_fraction"])
                        CSV.write(joinpath(readset_dir, "coverm_truth.tsv"), truth_coverm)
                    end
                end
            end
        end
    end
end

println("Round-trip benchmark complete. End time: $(Dates.now())")
