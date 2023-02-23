# function normalize_fasta(fasta_file, outdir)
#     mkpath("$(outdir)/normalized_fasta")
#     normalized_fasta_file = "$(outdir)/normalized_fasta/$(basename(fasta_file))"
#     if !isfile(normalized_fasta_file)
#         fasta_in = FASTX.FASTA.Reader(open(fasta_file))
#         fasta_out = FASTX.FASTA.Writer(open(normalized_fasta_file, "w"))
#         for (i, record) in enumerate(fasta_in)
#             updated_record = FASTX.FASTA.Record("$(i)", FASTX.FASTA.sequence(record))
#             write(fasta_out, updated_record)
#         end
#         close(fasta_in)
#         close(fasta_out)
#     end
#     return normalized_fasta_file
# end

# function run_prokka(ID, OUT_DIR, normalized_fasta_file)
#     prokka_dir="$(OUT_DIR)/prokka"
#     if !isdir(prokka_dir)
#         mkdir(prokka_dir)
#     end
#     prokka_cmd = `prokka --force --cpus 1 --outdir $(prokka_dir) --prefix $(ID) $(normalized_fasta_file)`
#     run(pipeline(prokka_cmd, stdout="$(prokka_dir)/prokka.out"))
#     return prokka_dir
# end

# function run_mlst(ID, OUT_DIR, normalized_fasta_file)
#     mlst_dir="$(OUT_DIR)/mlst"
#     if !isdir(mlst_dir)
#         mkdir(mlst_dir)
#     end
#     p = pipeline(
#             `mlst $(normalized_fasta_file)`,
#             stdout="$(mlst_dir)/$(ID).mlst.out")
#     run(p)
#     return mlst_dir
# end

# function run_phispy(ID, OUT_DIR, prokka_dir)
#     # 1 	prophage_coordinates.tsv
#     # 2 	GenBank format output
#     # 4 	prophage and bacterial sequences
#     # 8 	prophage_information.tsv
#     # 16 	prophage.tsv
#     # 32 	GFF3 format
#     # 64 	prophage.tbl
#     # 128 	test data used in the random forest
#     # 255   for all of them

#     phispy_dir="$(OUT_DIR)/phispy"
#     if !isdir(phispy_dir)
#         mkdir(phispy_dir)
#     end
#     if isempty(readdir(phispy_dir))
#         phisphy_cmd = `PhiSpy.py $(prokka_dir)/$(ID).gbk --output_dir $(phispy_dir) --file_prefix $(ID)-phispy --output_choice 255`
#         try
#             run(pipeline(phisphy_cmd, stdout="$(phispy_dir)/phisphy.out"))
#         catch
#             if isfile("$(phispy_dir)/$(ID)-phispy_prophage.gff3")
#                 @warn "phispy errored out after prophage gff3 was written"
#             else
#                 @error "phispy prophage gff3 not written"
#             end
#         end
#     end
#     return phispy_dir
# end

# function run_trnascan(ID, out_dir, normalized_fasta_file)
#     trnascan_dir = "$(out_dir)/trnascan"
#     # trnascan doesn't like to overwrite existing things
#     if !isdir(trnascan_dir)
#         mkdir(trnascan_dir)
#     end
#     if isempty(readdir(trnascan_dir))

#         #     -B for using Bacterial
#         trnascan_cmd = 
#         `tRNAscan-SE 
#             -B 
#             --output $(trnascan_dir)/$(ID).trnascan.out 
#             --bed $(trnascan_dir)/$(ID).trnascan.bed 
#             --fasta $(trnascan_dir)/$(ID).trnascan.fasta 
#             --struct $(trnascan_dir)/$(ID).trnascan.struct
#             --stats $(trnascan_dir)/$(ID).trnascan.stats 
#             --log $(trnascan_dir)/$(ID).trnascan.log
#             $(normalized_fasta_file)`
#         run(pipeline(trnascan_cmd, stdout="$(trnascan_dir)/trnascan.out", stderr="$(trnascan_dir)/trnascan.out"))
#     end
#     return trnascan_dir
# end

# function run_counterselection_spacer_detection(strain, out_dir, normalized_fasta_file)
#     counter_selection_dir = "$(out_dir)/counter-selection"
#     if !isdir(counter_selection_dir)
#         mkdir(counter_selection_dir)
#     end

#     if isempty(readdir(counter_selection_dir))
#         regex = BioSequences.biore"TTT[CG][ACGT]{25}"dna
#         k = 29
#         KMER_TYPE = BioSequences.BigDNAMer{k}

#         spacer_table = DataFrames.DataFrame(
#             ID = [],
#             contig = [],
#             PAM_and_spacer = [],
#             spacer = [],
#             strand = [],
#             start = [],
#             stop = [],
# #             free_energy = [],
# #             visualization_url = []
#         )

#         ProgressMeter.@showprogress for record in collect(FASTX.FASTA.Reader(open(normalized_fasta_file)))
#             for (i, kmer, reverse_complement_kmer) in BioSequences.each(KMER_TYPE, FASTX.FASTA.sequence(record))
#                 strand = missing
#                 if occursin(regex, kmer)
#                     strand = "+"
#                 elseif occursin(regex, reverse_complement_kmer)
#                     strand = "-"
#                     kmer = reverse_complement_kmer
#                 end
#                 if !ismissing(strand)
#                     spacer = BioSequences.DNAMer(kmer[i] for i in 5:length(kmer)) 
# #                     RNAfold_output = read(pipeline(`echo "$(string(spacer))"`, `RNAfold --noLP`), String)

# #                     rna_sequence, structure, free_energy = match(r"([ACGU]{25})\n([.()]{25})\s\(\s*(.*?)\)", RNAfold_output).captures
# #                     url = "http://nibiru.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence=$(rna_sequence)&structure=$(structure)"

#                     kmer_range = i:i+k-1

#                     t = DataFrames.DataFrame(
#                         ID = ID,
#                         contig = FASTX.FASTA.identifier(record),
#                         PAM_and_spacer = kmer,
#                         spacer = spacer,
#                         strand = strand,
#                         start = i,
#                         stop = i+k-1,
# #                         free_energy = free_energy,
# #                         visualization_url = url
#                     )
#                     spacer_table = vcat(spacer_table, t)
#                 end
#             end
#         end


#         uCSV.write(
#             "$(counter_selection_dir)/$(ID)-cpf1-spacers.tsv",
#             delim='\t',
#             data = collect(DataFrames.eachcol(spacer_table)),
#             header = DataFrames.names(spacer_table)
#         )
#         if isfile("$(out_dir)/rna.ps")
#             rm("$(out_dir)/rna.ps")
#         end
#     end
#     return counter_selection_dir
# end



# function run_amrfinderplus(ID, out_dir, protein_fasta)
#     amrfinderplus_dir = "$(out_dir)/amrfinderplus"
#     if !isdir(amrfinderplus_dir)
#         mkdir(amrfinderplus_dir)
#     end

#     if isempty(readdir(amrfinderplus_dir))
# #         run(`amrfinder -u`)
        
#         # because the pipeline is set up with some hacky CONDA path manipulation, 
#         # explictly setting amrfinder directory path to the location in the docker host
#         amrfinder_db_path = get(ENV, "AMRFINDER_DB", "none")

#         if amrfinder_db_path != "none"
#             cmd = 
#             `amrfinder
#             -p $(protein_fasta)
#             --plus
#             --output $(amrfinderplus_dir)/$(ID).amrfinderplus.tsv
#             -d $(amrfinder_db_path)
#             `
#         else
#             cmd = 
#             `amrfinder
#             -p $(protein_fasta)
#             --plus
#             --output $(amrfinderplus_dir)/$(ID).amrfinderplus.tsv
#             `
#         end
            

#         p = pipeline(cmd, 
#                 stdout="$(amrfinderplus_dir)/$(ID).amrfinderplus.out",
#                 stderr="$(amrfinderplus_dir)/$(ID).amrfinderplus.err")
#         run(p)
#     end
#     return amrfinderplus_dir
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run diamond search, returns path to diamond results.

```jldoctest
julia> 1 + 1
2
```
"""
function run_diamond(;
        identifier,
        out_dir,
        protein_fasta,
        diamond_db,
        force=false,
        outfile="$(identifier).prodigal.faa.diamond.txt"
    )
    diamond_dir = mkpath("$(out_dir)/diamond")

    # http://www.diamondsearch.org/index.php?pages/command_line_options/
    # --block-size/-b #Block size in billions of sequence letters to be processed at a time.  
    #     This is the main pa-rameter for controlling the program’s memory usage.  
    #     Bigger numbers will increase the useof memory and temporary disk space, but also improve performance.  
    #     The program can beexpected to use roughly six times this number of memory (in GB). So for the default value of-b2.0, 
    #     the memory usage will be about 12 GB
    system_memory_in_gigabytes = Int(Sys.total_memory()) / 1e9
    # reference says 6 but let's round upwards towards 8
    gb_per_block = 8
    block_size = system_memory_in_gigabytes / gb_per_block
    
    outfile = "$(diamond_dir)/$(outfile)"
    
    if force || !isfile(outfile)
        cmd = 
        `diamond blastp
        --threads $(Sys.CPU_THREADS)
        --block-size $(block_size)
        --db $(diamond_db)
        --query $(protein_fasta)
        --out $(outfile)
        --evalue 0.001
        --iterate
        --outfmt 6 qseqid qtitle qlen sseqid sallseqid stitle salltitles slen qstart qend sstart send evalue bitscore length pident nident mismatch staxids
        `

        # --un                     file for unaligned queries
        # --al                     file or aligned queries
        # --unfmt                  format of unaligned query file (fasta/fastq)
        # --alfmt                  format of aligned query file (fasta/fastq)
        # --unal                   report unaligned queries (0=no, 1=yes)

#         Value 6 may be followed by a space-separated list of these keywords:

#         qseqid means Query Seq - id
#         qtitle means Query title
#         qlen means Query sequence length
#         sseqid means Subject Seq - id
#         sallseqid means All subject Seq - id(s), separated by a ';'
#         stitle means Subject Title
#         salltitles means All Subject Title(s), separated by a '<>'
#         slen means Subject sequence length
#         qstart means Start of alignment in query
#         qend means End of alignment in query
#         sstart means Start of alignment in subject
#         send means End of alignment in subject
#         evalue means Expect value
#         bitscore means Bit score
#         length means Alignment length
#         pident means Percentage of identical matches
#         nident means Number of identical matches
#         mismatch means Number of mismatches
#         staxids means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)
        
        @time run(pipeline(cmd))
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

```jldoctest
julia> 1 + 1
2
```
"""
function run_mmseqs_easy_search(;out_dir, query_fasta, target_database, outfile, force=false)
    out_dir = mkpath(joinpath(out_dir, "mmseqs_easy_search"))
    outfile = joinpath(out_dir, outfile * ".mmseqs_easy_search." * basename(target_database) * ".txt")
    
    format_output = "query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits"
    
    if basename(target_database) in ["UniRef100", "UniRef90", "UniRef50", "UniProtKB", "TrEMBL", "Swiss-Prot", "NR", "GTDB", "SILVA", "Kalamari"]
        format_output *= ",taxid"
    end
    
    # note: cut exhaustive-search since it was taking far too long
    # --exhaustive-search
    # killed after 11 hours w/ 16 cores on UniRef100
    # running in base mode with UniRef100 @ 16 cores = 2h12m
    # could consider the iterative sensitivity search?
    # iterative was a bit faster and found matches for all of the same proteins
    #  # Increasing sensitivity search (from 2 to 7 in 3 steps)
    # mmseqs easy-search examples/QUERY.fasta examples/DB.fasta result.m8 tmp
    # aim for sensitivities 1, 3, 5, 7
    # --start-sens 1 -s 7 --sens-steps 3
    # 1 4 7
    # --start-sens 1 -s 7 --sens-steps 2
    if force || (!force && !isfile(outfile))
        cmd = 
        `mmseqs
            easy-search
            $(query_fasta)
            $(target_database)
            $(outfile)
            $(joinpath(out_dir, "tmp"))
            --format-mode 4
            --format-output $(format_output)
            --start-sens 1 -s 7 --sens-steps 3
        `
        @time run(pipeline(cmd))
    else
        @info "target outfile $(outfile) already exists, remove it or set force=true to re-generate"
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

```jldoctest
julia> 1 + 1
2
```
"""
function run_blast(;out_dir, fasta, blast_db, blast_command, force=false)
    blast_dir = mkpath(joinpath(out_dir, blast_command))
    outfile = "$(blast_dir)/$(basename(fasta)).$(blast_command).$(basename(blast_db)).txt"

    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
    if force || (!force && !isfile(outfile))
        cmd = 
        `
        $(blast_command)
        -num_threads $(Sys.CPU_THREADS)
        -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
        -query $(fasta)
        -db $(blast_db)
        -out $(outfile)
        -max_target_seqs 100
        -evalue 0.001
        `
#         p = pipeline(cmd, 
#                 stdout="$(blastn_dir)/$(ID).blastn.out",
#                 stderr="$(blastn_dir)/$(ID).blastn.err")
        @time run(pipeline(cmd))
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

```jldoctest
julia> 1 + 1
2
```
"""
function run_prodigal(;out_dir, fasta_file)
    prodigal_dir = mkpath("$(out_dir)/prodigal")

    # $ prodigal
    # -------------------------------------
    # PRODIGAL v2.6.3 [February, 2016]         
    # Univ of Tenn / Oak Ridge National Lab
    # Doug Hyatt, Loren Hauser, et al.     
    # -------------------------------------

    # Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
    #                  [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
    #                  [-p mode] [-q] [-s start_file] [-t training_file] [-v]

    #          -a:  Write protein translations to the selected file.
    #          -c:  Closed ends.  Do not allow genes to run off edges.
    #          -d:  Write nucleotide sequences of genes to the selected file.
    #          -f:  Select output format (gbk, gff, or sco).  Default is gbk.
    #          -g:  Specify a translation table to use (default 11).
    #          -h:  Print help menu and exit.
    #          -i:  Specify FASTA/Genbank input file (default reads from stdin).
    #          -m:  Treat runs of N as masked sequence; don't build genes across them.
    #          -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
    #          -o:  Specify output file (default writes to stdout).
    #          -p:  Select procedure (single or meta).  Default is single.
    #          -q:  Run quietly (suppress normal stderr output).
    #          -s:  Write all potential genes (with scores) to the selected file.
    #          -t:  Write a training file (if none exists); otherwise, read and use
    #               the specified training file.
    #          -v:  Print version number and exit.
    
    if isempty(readdir(prodigal_dir))
        cmd = 
        `prodigal
        -o $(prodigal_dir)/$(basename(fasta_file)).prodigal.gff
        -f gff
        -m
        -p meta
        -i $(fasta_file)
        -a $(prodigal_dir)/$(basename(fasta_file)).prodigal.faa
        -d $(prodigal_dir)/$(basename(fasta_file)).prodigal.fna
        -s $(prodigal_dir)/$(basename(fasta_file)).prodigal.all_potential_gene_scores.txt
        `

        p = pipeline(cmd, 
                stdout="$(prodigal_dir)/$(basename(fasta_file)).prodigal.out",
                stderr="$(prodigal_dir)/$(basename(fasta_file)).prodigal.err")
        run(p)
    end
    return prodigal_dir
end