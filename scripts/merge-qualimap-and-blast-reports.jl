    # # 2432 hits across 900 contigs
    @info "reading in blast report"
    blastn_results = Mycelia.parse_blast_report(blastn_report)
    # @show DataFrames.nrow(blastn_results)
    # @show length(unique(blastn_results[!, "query id"]))

    # viral_contigs = Set(unique(blastn_results[!, "query id"]));
    # # take everything that is a sensitive blastn hit against the viruses db, and then blast those against NT to verify there aren't better non-viral hits
    # viral_fasta = replace(assembled_fasta, ".fna" => ".viral.fna")
    # open(viral_fasta, "w") do io
    #     fastx_io = FASTX.FASTA.Writer(io)
    #     for record in FASTX.FASTA.Reader(open(assembled_fasta))
    #         if (FASTX.identifier(record) in viral_contigs)
    #             write(fastx_io, record)
    #         end
    #     end
    #     close(fastx_io)
    # end
    # @assert length(collect(Mycelia.open_fastx(viral_fasta))) == length(viral_contigs)

    # nt_blastdb_dir = "$(homedir())/workspace/blastdb"
    # nt_blast_db = "nt"
    # nt_blastdb_path = joinpath(nt_blastdb_dir, nt_blast_db)
    # if isdir(blastdb_dir) && !isempty(readdir(blastdb_dir))
    #     @info "blast db detected, using existing"
    # else
    #     Mycelia.download_blast_db(db=nt_blast_db, outdir=nt_blastdb_dir, source="ncbi")
    # end

    # @time blastn_nt_report = Mycelia.run_blastn(out_dir = SRR_path, fasta = viral_fasta, blast_db = nt_blastdb_path, task = "blastn")
    # errored out...

    # this actually had more viral hits than blastn against representative viruses, but only by a bit more
    # nt_megablast_report = joinpath(data_dir, "SRA/SRR6399459/blastn/final.contigs.fastg.gfa.fna.blastn.nt.megablast.txt")

    # intersect(hits, parse.(Int, blastn_results[!, "query id"]))

    # # 18741252 hits (too many!!) across 1142 contigs (not enough?!)
    # # only 17 contigs came back with viral classifications!!
    # nt_megablast_results = Mycelia.parse_blast_report(nt_megablast_report)
    # # @show DataFrames.nrow(nt_megablast_results)
    # # @show length(unique(nt_megablast_results[!, "query id"]))

    # create filtered list of contigs not identified by first pass
    # run dc-megablast to get second quick pass at genus-level identifications

    # read fasta file to determine list of contigs
    # fastx_identifiers = [FASTX.identifier(record) for record in FASTX.FASTA.Reader(open(assembled_fasta))]
    # megablast_hits = unique(map(x -> first(split(x, '\t')), Iterators.filter(l -> !occursin(r"^#", l), eachline(megablast_report))))
    # not_yet_annotated = setdiff(fastx_identifiers, megablast_hits)
    # not_yet_annotated_fasta = replace(assembled_fasta, ".fna" => ".not_yet_annotated.fna")
    # open(not_yet_annotated_fasta, "w") do io
    #     fastx_io = FASTX.FASTA.Writer(io)
    #     for record in FASTX.FASTA.Reader(open(assembled_fasta))
    #         if (FASTX.identifier(record) in not_yet_annotated)
    #             write(fastx_io, record)
    #         end
    #     end
    #     close(fastx_io)
    # end     
    # megablast_results = Mycelia.parse_blast_report(megablast_report)

    # Mycelia.run_blastn(out_dir = SRR_path, fasta = assembled_fasta, blast_db = db_path, task = "dc-megablast")

    # create filered list of contigs not identified by first or second pass

    # Mycelia.run_blastn(out_dir = SRR_path, fasta = assembled_fasta, blast_db = db_path, task = "blastn")

    # @time blast_report = Mycelia.run_blast(out_dir = SRR_path, fasta = assembled_fasta, blast_db = db_path, blast_command = "blastn")

    detected_tax_id_file = "$(assembled_fasta).detected_tax_ids.txt"
    open(detected_tax_id_file, "w") do io
        for taxid in unique(filter(!ismissing, blastn_results[!, "subject tax id"]))
            println(io, taxid)
        end
    end
    taxid_to_lineage_map = Dict(parse(Int, split_line[1]) => split_line[2] for split_line in split.(readlines(`taxonkit lineage --data-dir $(taxdump_dir) $(detected_tax_id_file)`), '\t'))

    qualimap_report_txt = joinpath(SRR_path, "megahit", "qualimap", "genome_results.txt")
    qualimap_contig_coverage_table = Mycelia.parse_qualimap_contig_coverage(qualimap_report_txt)

    blastn_results[!, "lineage"] = map(x -> get(taxid_to_lineage_map, x, ""), blastn_results[!, "subject tax id"])
    blastn_results[!, "% of subject length"] = round.(blastn_results[!, "query length"] ./ blastn_results[!, "subject length"] * 100, digits=3)
    contig_info_table = DataFrames.leftjoin(qualimap_contig_coverage_table, blastn_results, on="Contig" => "query id")

    # re-order columns based on utility
    reordered_columns = [
        "Contig",
        "Length",
        "Mapped bases",
        "Mean coverage",
        "Standard Deviation",
        "% Mapped bases",
        "subject id",
        "subject acc.",
        "subject title",
        "subject tax id",
        "lineage",
        "% identity",
        "% of subject length",
        "evalue",
        "bit score",
        "query length",
        "subject length",
        "alignment length",
        "q. start",
        "q. end",
        "s. start",
        "s. end",
        "identical",
        "mismatches"
    ]
    contig_info_table = contig_info_table[!, reordered_columns]
    sort!(contig_info_table, ["% Mapped bases", "bit score"], rev=true)
    # consider gzipping for large files!
    uCSV.write(contig_info_tsv, coalesce.(contig_info_table, ""), delim='\t')