import Pkg
Pkg.activate(".")
pkgs = [
    "ArgParse"
]
Pkg.add(pkgs)
for pkg in pkgs
    eval(Meta.parse("import $pkg"))
end

function parse_arguments()
    s = ArgParse.ArgParseSettings()

    ArgParse.@add_arg_table s begin
        "--fasta"
            help = "fasta file to blast"
            arg_type = String
        "--blast_db"
            help = "path to blast_db"
            arg_type = String
            required = true
        "--out_directory"
            help = "output directory"
            arg_type = String
            default = ""
        "--task"
            help = "type of blastn, options are blastn, megablast-dc, or megablast. Default = megablast"
            arg_type = String
            default = "megablast"
        "--threads"
            help = "number of threads to run"
            arg_type = Int
            default = 1
    end

    return ArgParse.parse_args(s)
end


function run_blastn(;out_dir, fasta, blast_db, task="megablast", force=false, remote=false, wait=true, threads=min(Sys.CPU_THREADS, 8))
    blast_dir = mkpath(joinpath(out_dir, "blastn"))
    outfile = "$(blast_dir)/$(basename(fasta)).blastn.$(basename(blast_db)).$(task).txt"
    if remote
        outfile = replace(outfile, ".txt" => ".remote.txt")
    end
    
    need_to_run = !isfile(outfile) || (filesize(outfile) == 0)
    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
        
    if force || need_to_run
        # if remote
        #         -remote
        # else
        # https://www.ncbi.nlm.nih.gov/books/NBK571452/
        # cap @ 8 and also use -mt_mode = 1 based on empirical evidence from
        # above blog post
        cmd = 
        `
        /scg/apps/legacy/software/blast/2.2.31+/bin/blastn
        -num_threads $(threads)
        -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
        -query $(fasta)
        -db /scg/apps/legacy/software/blast/db/$(blast_db)
        -out $(outfile)
        -max_target_seqs 10
        -task $(task)
        -evalue 0.001
        `
        @info "running cmd $(cmd)"
        # p = pipeline(cmd, 
        #         stdout=outfile * ".out",
        #         stderr=outfile * ".err")
        # p = pipeline(`module load blast`, cmd)
        p = pipeline(cmd)
        run(p, wait=wait)
    end
    return outfile
end

args = parse_arguments()
run_blastn(
    out_dir=args["out_directory"],
    fasta=args["fasta"],
    blast_db=args["blast_db"],
    task=args["task"],
    threads=args["threads"])