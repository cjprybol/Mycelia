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

function make_diamond_db(fasta_file, db_file=fasta_file)
    @time run(`diamond makedb --in $(fasta_file) -d $(db_file)`)
end

# in order to change this to be a standard blast where we don't need all pairwise hits
# just drop the parameters id, min-score, max-target-seqs
function pairwise_diamond(joint_fasta_file)
    if !isfile("$(joint_fasta_file).dmnd")
        make_diamond_db(joint_fasta_file)
    end
    n_records = count_records(joint_fasta_file)
    # max_target_seqs = Int(ceil(sqrt(n_records)))
    # @show "here!"
    sensitivity = "--iterate"
    # --block-size/-b
    # https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options
    # set block size to total memory / 8
    available_gigabytes = floor(Sys.free_memory() / 1e9)
    block_size = floor(available_gigabytes / 8)
    
    @time run(`diamond blastp $(sensitivity) --block-size $(block_size) --id 0 --min-score 0 --max-target-seqs $(n_records) --unal 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore -d $(joint_fasta_file).dmnd -q $(joint_fasta_file) -o $(joint_fasta_file).dmnd.tsv`)
    # # pairwise output is all of the alignments, super helpful!
    # # @time run(`diamond blastp $(sensitivity) --id 0 --min-score 0 --max-target-seqs $(N_RECORDS) --unal 1 --outfmt 0  -d $(joint_fasta_outfile).dmnd -q $(joint_fasta_outfile) -o $(joint_fasta_outfile).diamond.pairwise.txt`)
end


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
function run_mmseqs_easy_taxonomy(;out_dir, query_fasta, target_database, outfile, force=false)
    out_dir = mkpath(joinpath(out_dir, "mmseqs_easy_taxonomy"))
    outfile = joinpath(out_dir, outfile * ".mmseqs_easy_taxonomy." * basename(target_database) * ".txt")
    # note I tried adjusting all of the following, and none of them improved overall runtime
    # in any meaningful way
    # -s FLOAT                         Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]
    # https://github.com/soedinglab/MMseqs2/issues/577#issuecomment-1191584081
    # apparently orf-filter 1 speeds up by 50%!
    # --orf-filter INT                 Prefilter query ORFs with non-selective search
    #                               Only used during nucleotide-vs-protein classification
    #                               NOTE: Consider disabling when classifying short reads [0]
    # --lca-mode INT                   LCA Mode 1: single search LCA , 2/3: approximate 2bLCA, 4: top hit [3]
    # --lca-search BOOL                Efficient search for LCA candidates [0]
    # ^ this looks like it actually runs @ 1 with s=1.0 & --orf-filter=1
    
    # 112 days to process 600 samples at this rate....
    # 278 minutes or 4.5 hours for a single sample classification!!
    # 16688.050696 seconds (1.43 M allocations: 80.966 MiB, 0.02% gc time, 0.00% compilation time)
    # this is for default parameters
    # lowering sensitivity and taking LCA
    # 16590.725343 seconds (1.06 M allocations: 53.487 MiB, 0.00% compilation time)
    # took just as long!
    # difference was only 10 minutes
    # 15903.218456 seconds (969.92 k allocations: 48.624 MiB, 0.01% gc time)
    # use default parameters
    
    if force || (!force && !isfile(outfile))
        cmd = 
        `mmseqs
         easy-taxonomy
         $(query_fasta)
         $(target_database)
         $(outfile)
         $(joinpath(out_dir, "tmp"))
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
function run_blastn(;out_dir, fasta, blast_db, task="megablast", force=false, remote=false, wait=true)
    blast_dir = mkpath(joinpath(out_dir, "blastn"))
    outfile = "$(blast_dir)/$(basename(fasta)).blastn.$(basename(blast_db)).$(task).txt"
    # if remote
        # outfile = replace(outfile, ".txt" => ".remote.txt")
    # end
    
    need_to_run = !isfile(outfile) || (filesize(outfile) == 0)
    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
    
    # I want to speed this up more but don't know how
    # 
    # num_alignments
    # https://www.ncbi.nlm.nih.gov/books/NBK569845/
    # Windowmasker masks the over-represented sequence data and it can also mask the low complexity sequence data using the built-in dust algorithm (through the -dust option). To mask low-complexity sequences only, we will need to use dustmasker.
    # http://ftp.ncbi.nlm.nih.gov/pub/agarwala/dustmasker/README.dustmasker
    # http://ftp.ncbi.nlm.nih.gov/pub/agarwala/windowmasker/README.windowmasker
    # ./windowmasker -ustat ustat.15 -in chr1.fa -out chr1.wm -dust true
    
    # windowmasker -in hs_chr -infmt blastdb -mk_counts -parse_seqids -out hs_chr_mask.counts
    # windowmasker -in hs_chr -infmt blastdb -ustat hs_chr_mask.count -outfmt maskinfo_asn1_bin -parse_seqids -out hs_chr_mask.asnb
    
    # dustmasker -in hs_chr -infmt blastdb -parse_seqids -outfmt maskinfo_asn1_bin -out hs_chr_dust.asnb
    
    # makeblastdb -in hs_chr –input_type blastdb -dbtype nucl -parse_seqids -mask_data hs_chr_mask.asnb -out hs_chr -title "Human Chromosome, Ref B37.1"
    # blastdbcmd -db hs_chr -info
    
    # https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/winmasker/README

    # windowmasker -mk_counts -checkdup true -infmt blastdb -in nt -out nt.wm.1
    # windowmasker -ustat nt.wm.1 -dust true -in nt -infmt blastdb -out nt.wm.2
    # makeblastdb -in nt –input_type blastdb -dbtype nucl -parse_seqids -mask_data nt.wm.2 -out nt_masked_deduped -title "nt masked and deduped"
    
    # windowmasker -mk_counts -infmt blastdb -in nt -out nt.wm.no_dup_check.1
    # windowmasker -ustat nt.wm.no_dup_check.1 -dust true -in nt -infmt blastdb -out nt.wm.no_dup_check.2
    # makeblastdb -in nt –input_type blastdb -dbtype nucl -parse_seqids -mask_data nt.wm.no_dup_check.2 -out nt_masked -title "nt masked"

    # windowmasker -convert -in input_file_name -out output_file_name [-sformat output_format] [-smem available_memory]
    
    if force || need_to_run
        # if remote
        #     cmd = 
        #         `
        #         blastn
        #         -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
        #         -query $(fasta)
        #         -db $(basename(blast_db))
        #         -out $(outfile)
        #         -max_target_seqs 10
        #         -evalue 0.001
        #         -task $(task)
        #         -soft_masking true
        #         -subject_besthit
        #         -dust
        #         -remote
        #         `
        # else
        # https://www.ncbi.nlm.nih.gov/books/NBK571452/
        # cap @ 8 and also use -mt_mode = 1 based on empirical evidence from
        # above blog post
        cmd = 
        `
        blastn
        -num_threads $(min(Sys.CPU_THREADS, 8))
        -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
        -query $(fasta)
        -db $(blast_db)
        -out $(outfile)
        -max_target_seqs 10
        -subject_besthit
        -task $(task)
        -evalue 0.001
        `
        # end
#         p = pipeline(cmd, 
#                 stdout="$(blastn_dir)/$(ID).blastn.out",
#                 stderr="$(blastn_dir)/$(ID).blastn.err")
        @info "running cmd $(cmd)"
        @time run(pipeline(cmd), wait=wait)
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
function run_blast(;out_dir, fasta, blast_db, blast_command, force=false, remote=false, wait=true)
    blast_dir = mkpath(joinpath(out_dir, blast_command))
    outfile = "$(blast_dir)/$(basename(fasta)).$(blast_command).$(basename(blast_db)).txt"
    if remote
        outfile = replace(outfile, ".txt" => ".remote.txt")
    end
    
    need_to_run = !isfile(outfile) || (filesize(outfile) == 0)
    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
    if force || need_to_run
        if remote
            cmd = 
                `
                $(blast_command)
                -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
                -query $(fasta)
                -db $(basename(blast_db))
                -out $(outfile)
                -max_target_seqs 10
                -evalue 0.001
                -remote
                `
        else
            cmd = 
            `
            $(blast_command)
            -num_threads $(Sys.CPU_THREADS)
            -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
            -query $(fasta)
            -db $(blast_db)
            -out $(outfile)
            -max_target_seqs 10
            -evalue 0.001
            `
        end
#         p = pipeline(cmd, 
#                 stdout="$(blastn_dir)/$(ID).blastn.out",
#                 stderr="$(blastn_dir)/$(ID).blastn.err")
        @info "running cmd $(cmd)"
        @time run(pipeline(cmd), wait=wait)
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
function run_prodigal(;fasta_file, out_dir="", use_conda=false)
    if isempty(out_dir)
        prodigal_dir = mkpath("$(fasta_file)_prodigal")
    end

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
        if !use_conda
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
        else
            cmd = 
            `conda run --no-capture-output -n prodigal prodigal
            -o $(prodigal_dir)/$(basename(fasta_file)).prodigal.gff
            -f gff
            -m
            -p meta
            -i $(fasta_file)
            -a $(prodigal_dir)/$(basename(fasta_file)).prodigal.faa
            -d $(prodigal_dir)/$(basename(fasta_file)).prodigal.fna
            -s $(prodigal_dir)/$(basename(fasta_file)).prodigal.all_potential_gene_scores.txt
            `
        end

        p = pipeline(cmd, 
                stdout="$(prodigal_dir)/$(basename(fasta_file)).prodigal.out",
                stderr="$(prodigal_dir)/$(basename(fasta_file)).prodigal.err")
        run(p)
    end
    return prodigal_dir
end

# ```
# https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/

# GFF3 has 9 required fields, though not all are utilized (either blank or a default value of ‘.’).

#     Sequence ID
#     Source
#         Describes the algorithm or the procedure that generated this feature. Typically Genescane or Genebank, respectively.
#     Feature Type
#         Describes what the feature is (mRNA, domain, exon, etc.).
#         These terms are constrained to the [Sequence Ontology terms](http://www.sequenceontology.org/).
#     Feature Start
#     Feature End
#     Score
#         Typically E-values for sequence similarity and P-values for predictions.
#     Strand
#     Phase
#         Indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon.
#     Atributes
#         A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent . You can see the full list [here](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

# ```

"""
    create_chromosome_genedata_table(chromosome)

Take a chromosome from GenomicAnnotations.jl in GFF (and possibly genbank)
and return the formed dataframe.
"""
function create_chromosome_genedata_table(chromosome)
    # genedata is already provided as a dataframe with all of the Attributes as columns
    table = copy(chromosome.genedata)
    
    # the rest of these need to be created
    # I'm inserting them directly into their correct locations
    DataFrames.insertcols!(table, 1, "sequence-id" => fill(chromosome.name, DataFrames.nrow(table)))
    
    DataFrames.insertcols!(table, 3, "feature" => GenomicAnnotations.feature.(chromosome.genes))

    loci = getproperty.(GenomicAnnotations.locus.(chromosome.genes), :position)
    DataFrames.insertcols!(table, 4, "start" => first.(loci))
    DataFrames.insertcols!(table, 5, "stop" => last.(loci))
    
    DataFrames.insertcols!(table, 7, "strand" => .!GenomicAnnotations.iscomplement.(chromosome.genes))        
    
    return table
end

# function create_chromosome_genedata_table(chromosome)
#     table = chromosome.genedata
#     table[!, "chromosome"] .= chromosome.name
#     table[!, "feature"] = GenomicAnnotations.feature.(chromosome.genes)
#     table[!, "strand"] = .!GenomicAnnotations.iscomplement.(chromosome.genes)
#     table[!, "locus"] = getproperty.(GenomicAnnotations.locus.(chromosome.genes), :position)
#     table[!, "start"] = first.(chromosome.genedata[!, "locus"])
#     table[!, "stop"] = last.(chromosome.genedata[!, "locus"])
#     return table
# end

"""
    gff_to_table(gff)

Convert GenomicAnnotations.jl style GFF collection of chromosome vectors
into a wide-form GFF table (Attributes column expanded such that every attribute has it's own column)
"""
function gff_to_table(gff)
    # this is a no-op if already collected
    gff = collect(gff)
    table = create_chromosome_genedata_table(first(gff))
    for chromosome in gff[2:end]
        this_table = create_chromosome_genedata_table(chromosome)
        table = vcat(table, this_table)
    end
    return table
end

# TODO: switch to using GenomicAnnotations if GFF3 package isn't updated
# PR -> https://github.com/BioJulia/GFF3.jl/pull/12
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

```jldoctest
julia> 1 + 1
2
```
"""
function get_gff(;db=""::String, accession=""::String, ftp=""::String)
    if !isempty(db) && !isempty(accession)
        # API will block if we request more than 3 times per second, so set a 1/2 second sleep to set max of 2 requests per second when looping
        sleep(0.5)
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=gff3&id=$(accession)"
        return GenomicAnnotations.GFF.Reader(IOBuffer(HTTP.get(url).body))
        # return GFF3.Reader(IOBuffer(HTTP.get(url).body))
    elseif !isempty(ftp)
        # return GFF3.Reader(CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body)))
        return GenomicAnnotations.GFF.Reader(CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body)))
    else
        @error "invalid call"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

```jldoctest
julia> 1 + 1
2
```
"""
function get_genbank(;db=""::String, accession=""::String, ftp=""::String)
    if !isempty(db) && !isempty(accession)
        # API will block if we request more than 3 times per second, so set a 1/2 second sleep to set max of 2 requests per second when looping
        sleep(0.5)
        # url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=genbank&id=$(accession)&rettype=text"
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=genbank&id=$(accession)&retmode=text"
        # readgbk can't read from an io buffer, so need to download to a temp file
        # outfile = tempname()
        # open(outfile, "w") do io
        #     write(io, HTTP.get(url).body)
        # end
        # genbank_data = GenomicAnnotations.readgbk(outfile)
        # rm(outfile)
        # return genbank_data
        return GenomicAnnotations.GenBank.Reader(IOBuffer(HTTP.get(url).body))
    elseif !isempty(ftp)
        return GenomicAnnotations.GenBank.Reader(CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body)))
    else
        @error "invalid call"
    end
end