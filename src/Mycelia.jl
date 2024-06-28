module Mycelia

__precompile__(false)

# import ArgParse
import BioAlignments
import BioSequences
import BioSymbols
# import CairoMakie
import Clustering
import CodecZlib
import Conda
import CSV
import DataFrames
import DataStructures
import Dates
import DelimitedFiles
import Distances
import Distributions
import DocStringExtensions
import FASTX
import FileIO
import GenomicAnnotations
import GFF3
import GraphRecipes
import Graphs
import GLM
import HTTP
import JSON
import Kmers
import LsqFit
import MetaGraphs
import Mmap
import OrderedCollections
import Plots
import Primes
import ProgressMeter
import Random
import Statistics
import StatsBase
import StatsPlots
import XAM
import uCSV
import Downloads
import SparseArrays

import Pkg

# preserve definitions between code jldoctest code blocks
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Preserving-Definitions-Between-Blocks
# use this to build up a story as we go, where outputs of earlier defined functions feed into
# tests for further downstream functions

# if things fall out of date but look correct, update them automatically
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Fixing-Outdated-Doctests

const METADATA = joinpath(dirname(dirname(pathof(Mycelia))), "docs", "metadata")
const DNA_ALPHABET = BioSymbols.ACGT
const RNA_ALPHABET = BioSymbols.ACGU

# Mycelia.NERSC_MEM * .95
# const NERSC_MEM=512
# const NERSC_MEM=480
const NERSC_MEM=460
# const NERSC_CPU=240
const NERSC_CPU=240
# const AA_ALPHABET = filter(
#     x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x) || BioSymbols.isterm(x)),
#     BioSymbols.alphabet(BioSymbols.AminoAcid))
const AA_ALPHABET = filter(
    x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x)),
    BioSymbols.alphabet(BioSymbols.AminoAcid))

function find_mamba()
    try
        CONDA_RUNNER = strip(read(`which mamba`, String))
        return CONDA_RUNNER
    catch
        CONDA_RUNNER = joinpath(Conda.BINDIR, "mamba")
        return CONDA_RUNNER
    end
end

# can add support for conda too if needed
const CONDA_RUNNER = find_mamba()
const FASTQ_REGEX = r"\.(fq\.gz|fastq\.gz|fastq|fq)$"
const FASTA_REGEX = r"\.(fa\.gz|fasta\.gz|fna\.gz|fasta|fa|fna)$"
const VCF_REGEX = r"\.(vcf|vcf\.gz)$"

function add_bioconda_env(pkg; force=false)
    try
        current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
        if !(pkg in current_environments) || force
            @info "installing conda environment $(pkg)"
            run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
            run(`$(CONDA_RUNNER) clean --all -y`)
        else
            # @info "conda environment $(pkg) already present; set force=true to update/re-install"
        end
    catch
        add_bioconda_envs()
        current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
        if !(pkg in current_environments) || force
            @info "installing conda environment $(pkg)"
            run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
            run(`$(CONDA_RUNNER) clean --all -y`)
        else
            # @info "conda environment $(pkg) already present; set force=true to update/re-install"
        end
    end
end

function add_bioconda_envs(;all=false, force=false)
    if !isfile(CONDA_RUNNER) && (basename(CONDA_RUNNER) == "mamba")
        Conda.add("mamba")
    end
    if !isfile(joinpath(Conda.BINDIR, "pigz"))
        run(`$(CONDA_RUNNER) install pigz -y`)
    end
    current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
    # https://github.com/JuliaPy/Conda.jl/issues/185#issuecomment-1145149905
    if all
        for pkg in [
            "art",
            # "bioconvert",
            "badread",
            "bcftools",
            "bedtools",
            "blast",
            "clair3-illumina",
            "clair3",    
            # "bwa",
            # "bwa-mem2",
            # "deepvariant",
            "emboss",
            "filtlong",
            # "freebayes",
            "flye",
            "gatk4",
            # "gffread",
            "htslib",
            "megahit",
            "medaka",
            "minimap2",
            "mmseqs2",
            "nanocaller",
            "nanovar",
            # "nanoq",
            # "nanosim",
            # "nanosim-h",
            "ncbi-datasets-cli",
            "pggb",
            "picard",
            # "polypolish",
            "prodigal",
            "raven-assembler",
            "rtg-tools",
            "samtools",
            "sniffles",
            "sourmash",
            "spades",
            "tabix",
            "transtermhp",
            "trim-galore",
            "vcftools",
            "vg"
            ]
            if !(pkg in current_environments) || force
                @info "installing conda environment $(pkg)"
                add_bioconda_env(pkg)
            else
                @info "conda environment $(pkg) already present; set force=true to update/re-install"
            end
        end
    end
    run(`$(CONDA_RUNNER) clean --all -y`)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a command to SLURM using sbatch
"""
function scg_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=pwd(),
        partition::String,
        account::String,
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 8,
        cmd::String
    )
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --partition=$(partition)
    --account=$(account)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a command to SLURM using sbatch

https://docs.nersc.gov/jobs/policy/
https://docs.nersc.gov/systems/perlmutter/architecture/#cpu-nodes

default is to use shared qos

use
- regular
- preempt (reduced credit usage but not guaranteed to finish)
- premium (priorty runs limited to 5x throughput)

max request is 512Gb memory and 128 cores per node

https://docs.nersc.gov/systems/perlmutter/running-jobs/#tips-and-tricks
"""
function nersc_sbatch_shared(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=pwd(),
        qos::String="shared",
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 2,
        cmd::String,
        constraint::String="cpu"
    )
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --qos=$(qos)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --constraint=cpu
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a command to SLURM using sbatch

https://docs.nersc.gov/jobs/policy/
https://docs.nersc.gov/systems/perlmutter/architecture/#cpu-nodes

default is to use shared qos

use
- regular
- preempt (reduced credit usage but not guaranteed to finish)
- premium (priorty runs limited to 5x throughput)

https://docs.nersc.gov/systems/perlmutter/running-jobs/#tips-and-tricks
"""
function nersc_sbatch_regular(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=pwd(),
        qos::String="regular",
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=Mycelia.NERSC_CPU,
        mem_gb::Int=Int(floor(Mycelia.NERSC_MEM * .9)),
        cmd::String,
        constraint::String="cpu"
    )
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --qos=$(qos)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --constraint=cpu
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a command to SLURM using sbatch

https://docs.nersc.gov/jobs/policy/
https://docs.nersc.gov/systems/perlmutter/architecture/#cpu-nodes

default is to use shared qos

use
- regular
- preempt (reduced credit usage but not guaranteed to finish)
- premium (priorty runs limited to 5x throughput)

https://docs.nersc.gov/systems/perlmutter/running-jobs/#tips-and-tricks
"""
function nersc_sbatch_premium(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=pwd(),
        qos::String="premium",
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=Mycelia.NERSC_CPU,
        mem_gb::Int=Int(floor(Mycelia.NERSC_MEM * .9)),
        cmd::String,
        constraint::String="cpu"
    )
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --qos=$(qos)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --constraint=cpu
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

function fasta_genome_size(fasta_file)
    return reduce(sum, map(record -> length(FASTX.sequence(record)), Mycelia.open_fastx(fasta_file)))
end

function gfa_to_fasta(;gfa, fasta=gfa * ".fna")
    Mycelia.add_bioconda_env("gfatools")
    p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n gfatools gfatools gfa2fa $(gfa)`, fasta)
    run(p)
    return fasta
    # collect(GraphicalFragmentAssembly.Reader(open(primary_contig_gfa)))
    # open(fasta, "w") do io
    #     fastx_io = FASTX.FASTA.Writer(io)
    #     gfa_graph = Mycelia.parse_gfa(gfa)
    #     for v in Graphs.vertices(gfa_graph)
    #         record = FASTX.FASTA.Record(gfa_graph.vprops[v][:identifier], gfa_graph.vprops[v][:sequence])
    #         write(fastx_io, record)
    #     end
    #     close(fastx_io)
    # end
end

function determine_fasta_coverage(bam)
    genome_coverage_file = bam * ".coverage.txt"
    if !isfile(genome_coverage_file)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bedtools bedtools genomecov -d -ibam $(bam)`, genome_coverage_file))
    end
    return genome_coverage_file
end

#outdir="$(homedir())/software/bandage"
# I don't think that this is very portable - assumes sudo and linux
# can make a bandage_jll to fix this longer term
function download_bandage(outdir="/usr/local/bin")
    bandage_executable = joinpath(outdir, "Bandage")
    if !isfile(bandage_executable)
        run(`wget https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip`)
        run(`unzip Bandage_Ubuntu_static_v0_8_1.zip`)
        isfile("sample_LastGraph") && rm("sample_LastGraph")
        isfile("Bandage_Ubuntu_static_v0_8_1.zip") && rm("Bandage_Ubuntu_static_v0_8_1.zip")
        run(`sudo mv Bandage $(outdir)`)
    end
    return bandage_executable
end

function annotate_fasta(;
        fasta,
        identifier = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
        basedir = pwd(),        
        mmseqsdb = "$(homedir())/workspace/mmseqs/UniRef50",
        threads=Sys.CPU_THREADS
    )
    # @show basedir
    outdir = joinpath(basedir, identifier)
    @assert outdir != fasta
    if !isdir(outdir)
        @show isdir(outdir)
        mkpath(outdir)
        f = joinpath(outdir, basename(fasta))
        # make this an rclone copy for portability
        cp(fasta, f, force=true)
        Mycelia.run_prodigal(fasta_file=f)
        nucleic_acid_fasta = f * ".prodigal.fna"
        amino_acid_fasta = f * ".prodigal.faa"
        gff_file = f * ".prodigal.gff"

        mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
        mmseqs_gff_file = Mycelia.write_gff(gff = Mycelia.update_gff_with_mmseqs(gff_file, mmseqs_outfile), outfile = mmseqs_outfile * ".gff")
        transterm_gff_file = Mycelia.transterm_output_to_gff(Mycelia.run_transterm(fasta=f))

        joint_gff = Mycelia.write_gff(
            gff=sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file)), ["#seqid", "start", "end"]),
            outfile=f * ".gff")
        annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f, gff=joint_gff, genbank = joint_gff * ".genbank")

        transterm_gff_file_raw_fasta = Mycelia.transterm_output_to_gff(Mycelia.run_transterm(fasta=f))
        joint_gff_raw_fasta = Mycelia.write_gff(
            gff=sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file_raw_fasta)), ["#seqid", "start", "end"]),
            outfile=f * ".transterm_raw.gff")
        annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f, gff=joint_gff_raw_fasta, genbank = joint_gff_raw_fasta * ".genbank")
    else
        @info "$(outdir) already present, skipping..."
    end
    return outdir
end

function rclone_list_directories(path)
    directories = [join(split(line)[5:end], " ") for line in eachline(open(`rclone lsd $(path)`))]
    directories = joinpath.(path, directories)
    return directories
end

function mmseqs_pairwise_search(;fasta, output=fasta*".mmseqs_easy_search_pairwise")
    Mycelia.add_bioconda_env("mmseqs2")
    mkpath(output)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-search
        $(fasta)
        $(fasta)
        $(output)/$(basename(fasta)).mmseqs_pairwise_search.txt $(tempdir())
        --format-mode 4
        --format-output query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits
        --start-sens 1 -s 7 --sens-steps 7 --sort-results 1 --remove-tmp-files 1 --search-type 3`)
    return output
end

function mmseqs_easy_linclust(;fasta, output=fasta*".mmseqs_easy_linclust", tmp=tempdir())
    Mycelia.add_bioconda_env("mmseqs2")   
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createdb $(fasta) $(fasta)_DB`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createindex --search-type 3 $(fasta)_DB $(tempdir())`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-linclust $(fasta)_DB $(output) $(tmp)`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createtsv $(fasta)_DB $(fasta)_DB $(output) $(output).tsv`)
    return "$(output).tsv"
end

function mmseqs_easy_cluster(;fasta, output=fasta*".mmseqs_easy_cluster", tmp=tempdir())
    Mycelia.add_bioconda_env("mmseqs2")
    run(`Mycelia.CONDA_RUNNER run --live-stream -n mmseqs2 mmseqs easy-cluster $(fasta) $(output) $(tmp) --min-seq-id 0.5 -c 0.8 --cov-mode 1`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createtsv $(fasta) $(fasta) $(output) $(output).tsv`)
    return "$(output).tsv"
end

function filter_long_reads(;
        in_fastq,
        out_fastq = replace(in_fastq, r"\.(fq\.gz|fastq\.gz|fastq|fq)$" => ".filtlong.fq.gz"),
        min_mean_q = 20,
        keep_percent = 95
    )
    Mycelia.add_bioconda_env("filtlong")
    p1 = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n filtlong filtlong --min_mean_q $(min_mean_q) --keep_percent $(keep_percent) $(in_fastq)`,
        `pigz`
    )
    p2 = pipeline(p1, out_fastq)
    return p2
end

function export_blast_db(;path_to_db, fasta = path_to_db * ".fna.gz")
    Mycelia.add_bioconda_env("blast")
    if !isfile(fasta)
        # -long_seqids adds GI identifiers - these are cross-referenceable through other means so I'm dropping
        @time run(pipeline(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd  -entry all -outfmt '%f' -db $(path_to_db)`, `pigz`), fasta))
    else
        @info "$(fasta) already present"
    end
end

function fastx_stats(fastq)
    Mycelia.add_bioconda_env("seqkit")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit stats $(fastq)`)
end

function subsample_reads_seqkit(;in_fastq::String, out_fastq::String="", n_reads::Union{Missing, Int}=missing, proportion_reads::Union{Missing, Float64}=missing)
    Mycelia.add_bioconda_env("seqkit")
    if ismissing(n_reads) && ismissing(proportion_reads)
        error("please specify the number or proportion of reads")
    elseif !ismissing(n_reads)
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit sample --two-pass --number $(n_reads) $(in_fastq)`, `gzip`)
        if isempty(out_fastq)
            out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqkit.N$(n_reads).fq.gz")
        end
    elseif !ismissing(proportion_reads)
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit sample --proportion $(proportion_reads) $(in_fastq)`, `gzip`)
        if isempty(out_fastq)
            out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqkit.P$(proportion_reads).fq.gz")
        end
    end
    @assert !isempty(out_fastq)
    if !isfile(out_fastq)
        run(pipeline(p, out_fastq))
    else
        @info "$(out_fastq) already present"
    end
    return out_fastq
end

function parse_rtg_eval_output(f)
    # import CodecZlib
    flines = readlines(CodecZlib.GzipDecompressorStream(open(f)))
    header_line = last(filter(fline -> occursin(r"^#", fline), flines))
    header = lstrip.(split(header_line, "\t"), '#')
    data_lines = filter(fline -> !occursin(r"^#", fline), flines)
    if isempty(data_lines)
        data = [Float64[] for i in 1:length(header)]
    else
        data, h = uCSV.read(IOBuffer(join(data_lines, '\n')), delim='\t')
    end
    # data = [[parse(Float64, x)] for x in split(last(flines), '\t')]
    # @show data, header
    DataFrames.DataFrame(data, header)
end

# function subsample_reads_seqtk(;in_fastq::String, out_fastq="", n_reads::Union{Missing, Int}=missing, proportion_reads::Union{Missing, Float64}=missing)
#     Mycelia.add_bioconda_env("seqtk")
#     if ismissing(n_reads) && ismissing(proportion_reads)
#         error("please specify the number or proportion of reads")
#     elseif !ismissing(n_reads)
#         p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqtk seqtk sample -2 $(n_reads) $(in_fastq)`, `gzip`)
#         if isempty(out_fastq)
#             out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqtk.N$(n_reads).fq.gz")
#         end
#     elseif !ismissing(proportion_reads)
#         p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqtk seqtk sample -2 $(proportion_reads) $(in_fastq)`, `gzip`)
#         if isempty(out_fastq)
#             out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqtk.P$(proportion_reads).fq.gz")
#         end
#     end
#     @assert !isempty(out_fastq)
#     run(pipeline(p, out_fastq))
#     return out_fastq
# end

# subsample_reads_seqtk(in_fastq = fastq, n_reads=10)


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Will write out reads as SAM and also write out an error free SAM. Choose the reads from the version you want
"""
# # ? art short read
function simulate_short_reads()
    # $(Mycelia.MAMBA) run --live-stream -n art \
    # art_illumina \
    # --samout \
    # --errfree \
    # --paired \
    # --seqSys HS25 \
    # --len 150 \
    # --mflen 500 \
    # --sdev 10 \
    # --in $(fasta_file) \
    # --out $(fasta_file).art.$(coverage)x. \
    # --rcount
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

quantity is either fold coverage, or total bases sequenced - NOT TOTAL READS

To go by total reads, do # reads * 15,000 = quantity
"""
function simulate_pacbio_reads(;fasta, quantity, outfile=replace(fasta, Mycelia.FASTA_REGEX => ".badread.$(quantity).fq.gz"))
    if !isfile(outfile)
        Mycelia.add_bioconda_env("badread")
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --error_model pacbio2021 --qscore_model pacbio2021 --identity 30,3 --reference $(fasta) --quantity $(quantity)`, `gzip`)
        run(pipeline(p, outfile))
    else
        @info "$(outfile) already exists, skipping..."
    end
    return outfile
end
    
function simulate_nanopore_reads()
# badread simulate --reference ref.fasta --quantity 50x | gzip > reads.fastq.gz
end

function simulate_nearly_perfect_long_reads()
    # badread simulate --reference ref.fasta --quantity 50x --error_model random \
    # --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
    # --identity 30,3 --length 40000,20000 --start_adapter_seq "" --end_adapter_seq "" \
    # | gzip > reads.fastq.gz
end

# conda install -c bioconda kmer-jellyfish
# count, bc, info, stats, histo, dump, merge, query, cite, mem, jf
# cap at 4 threads, 8Gb per thread by default - this should be plenty fast enough for base usage, but open it up for higher performance!
function jellyfish_count(;fastx, k, threads=min(4, Sys.CPU_THREADS), max_mem=min(threads*8e9, (Sys.total_memory() / 2)), canonical=false, outfile = ifelse(canonical, "$(fastx).k$(k).canonical.jf", "$(fastx).k$(k).jf"))
    # @show fastx
    # @show k
    # @show threads
    # @show max_mem
    # @show canonical
    # @show outfile
    Mycelia.add_bioconda_env("kmer-jellyfish")
    mem = Int(floor(max_mem))
    
    # Usage: jellyfish mem [options] file:path+

    # Give memory usage information

    # The mem subcommand gives some information about the memory usage of
    # Jellyfish when counting mers. If one replace 'count' by 'mem' in the
    # command line, it displays the amount of memory needed. All the
    # switches of the count subcommand are supported, although only the
    # meaningful one for computing the memory usage are used.

    # If the '--size' (-s) switch is omitted and the --mem switch is passed
    # with an amount of memory in bytes, then the largest size that fit in
    # that amount of memory is returned.

    # The memory usage information only takes into account the hash to store
    # the k-mers, not various buffers (e.g. in parsing the input files). But
    # typically those will be small in comparison to the hash.

    # Options (default value in (), *required):
    #  -m, --mer-len=uint32                    *Length of mer
    #  -s, --size=uint64                        Initial hash size
    #  -c, --counter-len=Length in bits         Length bits of counting field (7)
    #  -p, --reprobes=uint32                    Maximum number of reprobes (126)
    #      --mem=uint64                         Return maximum size to fit within that memory
    #      --bc=peath                           Ignored switch
    #      --usage                              Usage
    #  -h, --help                               This message
    #      --full-help                          Detailed help
    #  -V, --version                            Version
    
    jellyfish_buffer_size = parse(Int, first(split(read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish mem --mer-len $(k) --mem $(mem)`, String))))
    # @show jellyfish_buffer_size    
    # Usage: jellyfish count [options] file:path+
    # Count k-mers in fasta or fastq files
    # Options (default value in (), *required):
    #  -m, --mer-len=uint32                    *Length of mer
    #  -s, --size=uint64                       *Initial hash size
    #  -t, --threads=uint32                     Number of threads (1)
    #      --sam=PATH                           SAM/BAM/CRAM formatted input file
    #  -F, --Files=uint32                       Number files open simultaneously (1)
    #  -g, --generator=path                     File of commands generating fast[aq]
    #  -G, --Generators=uint32                  Number of generators run simultaneously (1)
    #  -S, --shell=string                       Shell used to run generator commands ($SHELL or /bin/sh)
    #  -o, --output=string                      Output file (mer_counts.jf)
    #  -c, --counter-len=Length in bitsM         Length bits of counting field (7)
    #      --out-counter-len=Length in bytes    Length in bytes of counter field in output (4)
    #  -C, --canonical                          Count both strand, canonical representation (false)
    #      --bc=peath                           Bloom counter to filter out singleton mers
    #      --bf-size=uint64                     Use bloom filter to count high-frequency mers
    #      --bf-fp=double                       False positive rate of bloom filter (0.01)
    #      --if=path                            Count only k-mers in this files
    #  -Q, --min-qual-char=string               Any base with quality below this character is changed to N
    #      --quality-start=int32                ASCII for quality values (64)
    #      --min-quality=int32                  Minimum quality. A base with lesser quality becomes an N
    #  -p, --reprobes=uint32                    Maximum number of reprobes (126)
    #      --text                               Dump in text format (false)
    #      --disk                               Disk operation. Do not do size doubling (false)
    #  -L, --lower-count=uint64                 Don't output k-mer with count < lower-count
    #  -U, --upper-count=uint64                 Don't output k-mer with count > upper-count
    #      --timing=Timing file                 Print timing information
    #      --usage                              Usage
    #  -h, --help                               This message
    #      --full-help                          Detailed help
    #  -V, --version                            Version
    if canonical
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish count --canonical --size $(jellyfish_buffer_size) --threads $(threads) --mer-len $(k) --output $(outfile) /dev/fd/0`
    else
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish count --size $(jellyfish_buffer_size) --threads $(threads) --mer-len $(k) --output $(outfile) /dev/fd/0`
    end
    if occursin(r"\.gz$", fastx)
        open_cmd = `gzip -dc $(fastx)`
    else
        open_cmd = `cat $(fastx)`
    end
    temp_tab = outfile * ".tsv"
    tabular_counts = temp_tab * ".gz"
    
    if !isfile(tabular_counts)
        if !isfile(outfile)
            run(pipeline(open_cmd, cmd))
        end
        if !isfile(temp_tab)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish dump --column --tab --output $(temp_tab) $(outfile)`)
        end
        run(`gzip $(temp_tab)`)
    end
    
    # temp_fasta = outfile * ".fna"
    # fna_counts = temp_fasta * ".gz"
    # if !isfile(fna_counts)
    #     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish dump --output $(temp_fasta) $(outfile)`)
    #     run(`gzip $(temp_fasta)`)
    # end

    # Usage: jellyfish dump [options] db:path
    # Dump k-mer counts
    # By default, dump in a fasta format where the header is the count and
    # the sequence is the sequence of the k-mer. The column format is a 2
    # column output: k-mer count.
    # Options (default value in (), *required):
    #  -c, --column                             Column format (false)
    #  -t, --tab                                Tab separator (false)
    #  -L, --lower-count=uint64                 Don't output k-mer with count < lower-count
    #  -U, --upper-count=uint64                 Don't output k-mer with count > upper-count
    #  -o, --output=string                      Output file
    #      --usage                              Usage
    #  -h, --help                               This message
    #  -V, --version                            Version

    # Usage: jellyfish histo [options] db:path
    #  -l, --low=uint64                         Low count value of histogram (1)
    #  -h, --high=uint64                        High count value of histogram (10000)
    #  -f, --full                               Full histo. Don't skip count 0. (false)
    #  -o, --output=string                      Output file
    # histogram = outfile * ".histogram"
    # if !isfile(histogram)
    #     # 10000000000
    #     # ^ runs out of memory
    #     # 10000
    #     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish histo --low 1 --high 10000 --output $(histogram) $(outfile)`)
    # end
    # return (;outfile, fna_counts, tabular_counts, histogram)
    isfile(outfile) && rm(outfile)
    return tabular_counts
end


# 7
# 0.169580 seconds (587.93 k allocations: 39.458 MiB, 76.55% compilation time)
# 11
# 4.420040 seconds (194.49 k allocations: 11.566 MiB)
# 13
# 20.690993 seconds (521.75 k allocations: 37.253 MiB, 0.58% compilation time)
# 17
# 412.670050 seconds (529.86 k allocations: 37.007 MiB, 0.03% compilation time)
function jellyfish_counts_to_kmer_frequency_histogram(jellyfish_counts_file, outfile=replace(jellyfish_counts_file, r"\.tsv\.gz$" => ".count_histogram.tsv"))
    # sorting with LC_ALL=C is the biggest speed up here of anything I've found
    if !isfile(outfile)
        io = open(pipeline(
                `gzip -dc $(jellyfish_counts_file)`,
                `cut -f2`,
                Cmd(`sort --temporary-directory . --compress-program gzip --numeric --stable`, env=Dict("LC_ALL" => "C")),
                `uniq --count`,
                `sed 's/^ *//'`,
                `sed 's/ /\t/'`
                ))
        frequency_histogram_table = CSV.read(io, DataFrames.DataFrame, header=["number of kmers", "number of observations"], delim='\t')
        CSV.write(outfile, frequency_histogram_table, delim='\t')
    else
        @info "$(outfile) already exists"
    end
    return outfile
end


# Usage: jellyfish merge [options] input:string+

# Merge jellyfish databases

# Options (default value in (), *required):
#  -o, --output=string                      Output file (mer_counts_merged.jf)
#  -m, --min                                Compute min count instead of sum (false)
#  -M, --max                                Compute max count instead of sum (false)
#  -j, --jaccard                            Compute the jaccard and weighted jaccard similarities (false)
#  -L, --lower-count=uint64                 Don't output k-mer with count < lower-count
#  -U, --upper-count=uint64                 Don't output k-mer with count > upper-count
#      --usage                              Usage
#  -h, --help                               This message
#  -V, --version                            Version

# Usage: jellyfish query [options] file:path mers:string+

# Query a Jellyfish database

# Options (default value in (), *required):
#  -s, --sequence=path                      Output counts for all mers in sequence
#  -o, --output=path                        Output file (stdout)
#  -i, --interactive                        Interactive, queries from stdin (false)
#  -l, --load                               Force pre-loading of database file into memory (false)
#  -L, --no-load                            Disable pre-loading of database file into memory (false)
#  -U, --usage                              Usage
#  -h, --help                               This message
#  -V, --version                            Version



function jitter(x, n)
    return [x + rand() / 3 * (ifelse(rand(Bool), 1, -1)) for i in 1:n]
end

function fasta_to_reference_kmer_counts(;kmer_type, fasta)
    kmer_counts = Dict{kmer_type, Int}()
    for record in Mycelia.open_fastx(fasta)
        record_sequence = BioSequences.LongDNA{2}(FASTX.sequence(record))
        forward_counts = StatsBase.countmap(kmer for (i, kmer) in Kmers.EveryKmer{kmer_type}(record_sequence))
        reverse_counts = StatsBase.countmap(kmer for (i, kmer) in Kmers.EveryKmer{kmer_type}(BioSequences.reverse_complement(record_sequence)))
        record_counts = merge(+, forward_counts, reverse_counts)
        merge!(+, kmer_counts, record_counts)
    end
    return kmer_counts
end

function filter_long_reads(;
        in_fastq,
        out_fastq = replace(in_fastq, r"\.(fq\.gz|fastq\.gz|fastq|fq)$" => ".filtlong.fq.gz"),
        min_mean_q = 20,
        keep_percent = 95
    )
    Mycelia.add_bioconda_env("filtlong")
    p1 = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n filtlong filtlong --min_mean_q $(min_mean_q) --keep_percent $(keep_percent) $(in_fastq)`,
        `pigz`
    )
    p2 = pipeline(p1, out_fastq)
    return p2
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

My standard pacbio aligning and sorting. No filtering done in this step.

Use shell_only=true to get string command to submit to SLURM
"""
function map_pacbio_reads(;
        fastq,
        reference_fasta,
        temp_sam_outfile = fastq * "." * basename(reference_fasta) * "." * "minimap2.sam",
        outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz"),
        threads = Sys.CPU_THREADS,
        memory = Sys.total_memory(),
        shell_only = false
    )
    # 4G is the default
    # smaller, higher diversity databases do better with 5+ as the denominator - w/ <=4 they run out of memory
    index_chunk_size = "$(Int(floor(memory/5e9)))G"
    @show index_chunk_size
    @show threads
    Mycelia.add_bioconda_env("minimap2")
    # Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if shell_only
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_chunk_size) -ax map-hifi $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
        return cmd
    else
        if !isfile(outfile)
            map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_chunk_size) -ax map-hifi $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
            run(map)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`)
            @assert isfile(outfile)
        else
            @info "$(outfile) already present"
        end
    end
end

# """
# My standard pacbio aligning and sorting. No filtering done in this step.

# Use shell_only=true to get string command to submit to SLURM
# """
# function minimap_index_pacbio(;
#         reference_fasta,
#         outfile = replace(reference_fasta, Mycelia.FASTA_REGEX => ".pacbio.mmi"),
#         threads = Sys.CPU_THREADS,
#         shell_only = false
#     )
#     Mycelia.add_bioconda_env("minimap2")
#     Mycelia.add_bioconda_env("samtools")
#     if shell_only
#         cmd =
#         """
#         $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -ax map-pb $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
#         && $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile) \\
#         | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -bh -o $(outfile) \\
#         && rm $(temp_sam_outfile)
#         """
#         return cmd
#     else
#         if !isfile(outfile)
#             map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -ax map-pb $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
#             run(map)
#             p = pipeline(
#                 `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile)`,
#                 `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -bh -o $(outfile)`
#             )
#             run(p)
#             rm(temp_sam_outfile)
#         else
#             @info "$(outfile) already present"
#         end
#     end
# end

# function filter_short_reads()
# end

# function map_short_reads()
# end

function download_genome_by_accession(;accession, outdir=pwd())
    temp_fasta = joinpath(outdir, accession * ".fna")
    outfile = temp_fasta * ".gz"
    if !isfile(outfile)
        try
            # pull the entire record so that if the download fails we don't leave an empty file
            fasta_records = collect(Mycelia.get_sequence(db = "nuccore", accession = accession))
            open(temp_fasta, "w") do io
                fastx_io = FASTX.FASTA.Writer(io)
                for fasta_record in fasta_records
                    write(fastx_io, fasta_record)
                end
                close(fastx_io)
                run(`gzip $(temp_fasta)`)
                @assert isfile(outfile)
            end
        catch e
            println("An error occurred: ", e)
            
        end
    end
    return outfile
end

function download_genome_by_ftp(;ftp, outdir=pwd())
    url = Mycelia.ncbi_ftp_path_to_url(ftp_path=ftp, extension="genomic.fna.gz")
    outfile = joinpath(outdir, basename(url))
    if !isfile(outfile)
        return Downloads.download(url, outfile)
    else
        return outfile
    end
end

function normalized_current_datetime()
    return replace(Dates.format(Dates.now(), Dates.ISODateTimeFormat), r"[^\w]" => "")
end

function githash(;short=false)
    git_hash = rstrip(read(`git rev-parse HEAD`, String))
    if short
        git_hash = git_hash[1:8]
    end
    return git_hash
end

# CSV is too memory inefficient, the others too slow :(
# # using uCSV
# # k=11
# # 3.444974 seconds (24.58 M allocations: 1.374 GiB, 34.65% gc time, 16.90% compilation time)
# # k=13
# # 362.285866 seconds (357.11 M allocations: 20.550 GiB, 91.60% gc time)

# # using DelimitedFiles.readdlm
# # k=11
# # 2.386620 seconds (16.11 M allocations: 632.732 MiB, 34.16% gc time, 24.25% compilation time)
# # k=13
# # 82.888552 seconds (227.49 M allocations: 8.766 GiB, 82.01% gc time)

# # CSV
# # k=11
# # 12.328422 seconds (7.62 M allocations: 732.639 MiB, 19091.67% compilation time: <1% of which was recompilation)
# # k=13
# # 37.098948 seconds (89.38 k allocations: 2.354 GiB, 93.56% gc time)

# function parse_jellyfish_counts(tabular_counts)
#     # load in the data
#     @assert occursin(r"\.gz$", tabular_counts) "this expects gzipped jellyfish tabular counts"
#     io = CodecZlib.GzipDecompressorStream(open(tabular_counts))
#     canonical_kmer_counts_table = DataFrames.DataFrame(CSV.File(io; delim='\t', header=false))
#     DataFrames.rename!(canonical_kmer_counts_table, [:Column1 => :kmer, :Column2 => :count])
    
#     # recode the kmers from strings to fixed sized kmer types
#     unique_kmer_lengths = unique(length.(canonical_kmer_counts_table[!, "kmer"]))
#     @assert length(unique_kmer_lengths) == 1
#     k = first(unique_kmer_lengths)
#     canonical_kmer_counts_table[!, "kmer"] = Kmers.DNAKmer{k}.(canonical_kmer_counts_table[!, "kmer"])
    
#     return canonical_kmer_counts_table
# end


# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/taxonomy/taxonomy/
function ncbi_taxon_summary(taxa_id)
    Mycelia.add_bioconda_env("ncbi-datasets")
    p = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets datasets summary taxonomy taxon $(taxa_id) --as-json-lines`,
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets dataformat tsv taxonomy --template tax-summary`
        )
    return DataFrames.DataFrame(uCSV.read(open(p), delim='\t', header=1))
end

function nearest_prime(n::Int)
    if n < 2
        return 2
    end
    next_p = Primes.nextprime(n)
    prev_p = Primes.prevprime(n)
    if n - prev_p <= next_p - n
        return prev_p
    else
        return next_p
    end
end

function fibonacci_numbers_less_than(n::Int)
    if n <= 0
        return []
    elseif n == 1
        return [0]
    else
        fib = [0, 1]
        next_fib = fib[end] + fib[end-1]
        while next_fib < n
            push!(fib, next_fib)
            next_fib = fib[end] + fib[end-1]
        end
        return fib
    end
end

function ks(;min=0, max=10_000)
    # flip from all odd primes to only nearest to fibonnaci primes
    flip_point = 23
    vcat(
        filter(isodd, Primes.primes(min, flip_point)),
        filter(x -> x > flip_point, nearest_prime.(fibonacci_numbers_less_than(max)))
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function rclone_copy(source, dest)
    done = false
    sleep_timer = 60
    attempts = 0
    while !done && attempts < 3
        attempts += 1
        try
            # https://forum.rclone.org/t/google-drive-uploads-failing-http-429/34147/9
            # --tpslimit                                       Limit HTTP transactions per second to this
            # --drive-chunk-size SizeSuffix                    Upload chunk size (default 8Mi)
            # --drive-upload-cutoff SizeSuffix                 Cutoff for switching to chunked upload (default 8Mi)
            # not currently using these but they may become helpful
            # --drive-pacer-burst int                          Number of API calls to allow without sleeping (default 100)
            # --drive-pacer-min-sleep Duration                 Minimum time to sleep between API calls (default 100ms)
            cmd = `rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $(source) $(dest)`
            @info "copying $(source) to $(dest) with command: $(cmd)"
            run(cmd)
            done = true
        catch
            @info "copying incomplete, sleeping $(sleep_timer) seconds and trying again..."
            sleep(sleep_timer)
            sleep_timer *= 2
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function drop_empty_columns(table)
    is_empty_column = map(col -> all(isempty.(col)), eachcol(table))
    return table[!, .!is_empty_column]
end

# # Need to add hashdeep & logging
# function tarchive(;directory, tarchive=directory * ".tar.gz")
#     run(`tar --create --gzip --verbose --file=$(tarchive) $(directory)`)
#     return tarchive
# end

function install_hashdeep()
    if Sys.which("hashdeep") !== nothing
        println("hashdeep executable found")
    else
        println("hashdeep executable not found in PATH, installing")
        try
            run(`sudo apt install hashdeep -y`)
        catch
            run(`apt install hashdeep -y`)
        end
    end
end

# Need to add hashdeep & logging
function create_tarchive(;directory, tarchive=directory * ".tar.gz")
    directory = normpath(directory)
    tarchive = normpath(tarchive)
    output_dir = mkpath(dirname(tarchive))
    
    install_hashdeep()


    working_dir, source = splitdir(directory)
    if isempty(source)
        source = working_dir
        working_dir = pwd()
    end
    if isempty(working_dir)
        working_dir = pwd()
    end
    println("output_dir: $output_dir\n" *
            "working_dir: $working_dir\n" *
            "source: $source\n" *
            "target_file: $tarchive")
    log_file = tarchive * ".log"
    hashdeep_file = tarchive * ".hashdeep.dfxml"
    
    if !isfile(tarchive)
        run(`tar --create --gzip --verbose --file=$(tarchive) $(directory)`)
    end
    if !isfile(log_file)
        run(pipeline(`tar -tvf $(tarchive)`,log_file))
    end
    if !isfile(hashdeep_file)
        run(pipeline(`hashdeep -c md5,sha1,sha256 -b -d $(tarchive)`,hashdeep_file))
    end
    return tarchive
end

function tar_extract(;tarchive, directory=dirname(tarchive))
    run(`tar --extract --gzip --verbose --file=$(tarchive) --directory=$(directory)`)
    return directory
end

function samtools_index_fasta(;fasta)
    Mycelia.add_bioconda_env("samtools")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools faidx $(fasta)`)
end

# not a very good function yet, but good enough for the pinches I need it for
function parse_gfa(gfa)
    segments = Vector{String}()
    links = Vector{Pair{String, String}}()
    paths = Dict{String, Vector{String}}()
    for l in eachline(open(gfa))
        s = split(l, '\t')
        if first(s) == "S"
            # segment
            @show "here!"
            push!(segments, string(s[2]))
        elseif first(s) == "L"
            # link
            push!(links, string(s[2]) => string(s[4]))
        elseif first(s) == "P"
            # path
            paths[string(s[2])] = string.(split(replace(s[3], r"[+-]" => ""), ','))
        elseif first(s) == "A" # hifiasm https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#output-file-formats
            continue
        else
            println(l)
            error("unexpected line encountered while parsing GFA")
        end
    end

    g = Graphs.SimpleGraph(length(segments))

    for link in links
        (u, v) = link
        ui = findfirst(segments .== u)
        vi = findfirst(segments .== v)
        Graphs.add_edge!(g, ui => vi)
    end
    return g
end
    
function find_resampling_stretches(;record_kmer_solidity, solid_branching_kmer_indices)
    indices = findall(.!record_kmer_solidity)  # Find the indices of false values
    if isempty(indices)
        return UnitRange{Int64}[]
    end
    
    diffs = diff(indices)  # Calculate the differences between consecutive indices
    # @show diffs
    range_starts = [indices[1]]  # Start with the first false index
    range_ends = Int[]
    
    for (i, d) in enumerate(diffs)
        if d > 1
            push!(range_ends, indices[i])
            push!(range_starts, indices[i+1])
        end
    end
    
    push!(range_ends, indices[end])  # Add the last false index as a range end
    
    low_quality_runs = [(start, stop) for (start, stop) in zip(range_starts, range_ends)]
    
    resampling_stretches = UnitRange{Int64}[]
    
    for low_quality_run in low_quality_runs
        unders = filter(solid_branching_kmer -> solid_branching_kmer < first(low_quality_run), solid_branching_kmer_indices)
        overs = filter(solid_branching_kmer -> solid_branching_kmer > last(low_quality_run), solid_branching_kmer_indices)
        if isempty(overs) || isempty(unders)
            continue
        else
            nearest_under = maximum(unders)
            nearest_over = minimum(overs)
            push!(resampling_stretches, nearest_under:nearest_over)
        end
    end
    if !allunique(resampling_stretches)
        resampling_stretches = unique!(resampling_stretches)
    end
    return resampling_stretches
end

function fastq_record(;identifier, sequence, quality_scores)
    # Fastx wont parse anything higher than 93
    quality_scores = min.(quality_scores, 93)
    record_string = join(["@" * identifier, sequence, "+", join([Char(x+33) for x in quality_scores])], "\n")
    return FASTX.parse(FASTX.FASTQRecord, record_string)
end

function process_fastq_record(;record, kmer_graph, yen_k_shortest_paths_and_weights, yen_k=3)
    ordered_kmers = MetaGraphs.get_prop(kmer_graph, :ordered_kmers)
    likely_valid_kmers = Set(ordered_kmers[MetaGraphs.get_prop(kmer_graph, :likely_valid_kmer_indices)])
    kmer_to_index_map = MetaGraphs.get_prop(kmer_graph, :kmer_indices)
    branching_nodes_set = MetaGraphs.get_prop(kmer_graph, :branching_nodes)
    assembly_k = MetaGraphs.get_prop(kmer_graph, :assembly_k)
    transition_likelihoods = MetaGraphs.get_prop(kmer_graph, :transition_likelihoods)
    kmer_mean_quality = MetaGraphs.get_prop(kmer_graph, :kmer_mean_quality)
    kmer_total_quality = MetaGraphs.get_prop(kmer_graph, :kmer_total_quality)
    
    new_record_identifier = FASTX.identifier(record) * ".k$(assembly_k)"
    record_sequence = BioSequences.LongDNA{4}(FASTX.sequence(record))

    kmer_type = Kmers.DNAKmer{assembly_k}
    record_kmers = last.(collect(Kmers.EveryKmer{kmer_type}(record_sequence)))
    record_quality_scores = collect(FASTX.quality_scores(record))
    record_kmer_quality_scores = [record_quality_scores[i:i+assembly_k-1] for i in 1:length(record_quality_scores)-assembly_k+1]
    
    record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
    record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
    record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    
    # trim beginning of fastq
    initial_solid_kmer = findfirst(record_kmer_solidity)
    if isnothing(initial_solid_kmer)
        return record
    elseif initial_solid_kmer > 1
        record_kmers = record_kmers[initial_solid_kmer:end]
        record_kmer_quality_scores = record_kmer_quality_scores[initial_solid_kmer:end]
        record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
        record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
        record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    end
    initial_solid_kmer = 1
    
    # trim end of fastq
    last_solid_kmer = findlast(record_kmer_solidity)
    if last_solid_kmer != length(record_kmer_solidity)
        record_kmers = record_kmers[1:last_solid_kmer]
        record_kmer_quality_scores = record_kmer_quality_scores[1:last_solid_kmer]
        record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
        record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
        record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    end
    
    # identify low quality runs and the solid branchpoints we will use for resampling
    solid_branching_kmer_indices = findall(record_solid_branching_kmers)
    resampling_stretches = find_resampling_stretches(;record_kmer_solidity, solid_branching_kmer_indices)

    # nothing to do
    if isempty(resampling_stretches)
        return record
    end
    trusted_range = 1:max(first(first(resampling_stretches))-1, 1)
    
    new_record_kmers = record_kmers[trusted_range]
    new_record_kmer_qualities = record_kmer_quality_scores[trusted_range]
    
    
    for (i, resampling_stretch) in enumerate(resampling_stretches)
        starting_solid_kmer = record_kmers[first(resampling_stretch)]
        ending_solid_kmer = record_kmers[last(resampling_stretch)]
        
        current_quality_scores = record_quality_scores[resampling_stretch]
        u = kmer_to_index_map[starting_solid_kmer]
        v = kmer_to_index_map[ending_solid_kmer]
        if !haskey(yen_k_shortest_paths_and_weights, u => v)
            yen_k_result = Graphs.yen_k_shortest_paths(kmer_graph, u, v, Graphs.weights(kmer_graph), yen_k)
            yen_k_shortest_paths_and_weights[u => v] = Vector{Pair{Vector{Int}, Float64}}()
            for path in yen_k_result.paths
                path_weight = Statistics.mean([kmer_total_quality[ordered_kmers[node]] for node in path])
                path_transition_likelihoods = 1.0
                for (a, b) in zip(path[1:end-1], path[2:end])
                    path_transition_likelihoods *= transition_likelihoods[a, b]
                end
                joint_weight = path_weight * path_transition_likelihoods
                push!(yen_k_shortest_paths_and_weights[u => v], path => joint_weight)
            end
        end
        yen_k_path_weights = yen_k_shortest_paths_and_weights[u => v]      
        if length(yen_k_path_weights) > 1
            current_distance = length(resampling_stretch)
            initial_weights = last.(yen_k_path_weights)
            path_lengths = length.(first.(yen_k_path_weights))
            deltas = map(l -> abs(l-current_distance), path_lengths)
            adjusted_weights = initial_weights .* map(d -> exp(-d * log(2)), deltas)
            # make it more severe?
            # adjusted_weights = adjusted_weights.^2
            
            # and a bonus for usually being correct
            
            selected_path_index = StatsBase.sample(StatsBase.weights(adjusted_weights))
            selected_path, selected_path_weights = yen_k_path_weights[selected_path_index]
            selected_path_kmers = [ordered_kmers[kmer_index] for kmer_index in selected_path]
            
            if last(new_record_kmers) == first(selected_path_kmers)
                selected_path_kmers = selected_path_kmers[2:end]
            end
            append!(new_record_kmers, selected_path_kmers)
            selected_kmer_qualities = [Int8.(min.(typemax(Int8), floor.(kmer_mean_quality[kmer]))) for kmer in selected_path_kmers]
            append!(new_record_kmer_qualities, selected_kmer_qualities)
        else
            selected_path_kmers = record_kmers[resampling_stretch]
            if last(new_record_kmers) == first(selected_path_kmers)
                selected_path_kmers = selected_path_kmers[2:end]
            end
            append!(new_record_kmers, selected_path_kmers)
            selected_kmer_qualities = [Int8.(min.(typemax(Int8), floor.(kmer_mean_quality[kmer]))) for kmer in selected_path_kmers]
            append!(new_record_kmer_qualities, selected_kmer_qualities)
        end
        if i < length(resampling_stretches) # append high quality gap
            next_solid_start = last(resampling_stretch)+1
            next_resampling_stretch = resampling_stretches[i+1]
            next_solid_stop = first(next_resampling_stretch)-1
            if !isempty(next_solid_start:next_solid_stop)
                selected_path_kmers = record_kmers[next_solid_start:next_solid_stop]
                append!(new_record_kmers, selected_path_kmers)
                selected_kmer_qualities = record_kmer_quality_scores[next_solid_start:next_solid_stop]
                append!(new_record_kmer_qualities, selected_kmer_qualities)
            end
        else # append remainder of sequence
            @assert i == length(resampling_stretches)
            next_solid_start = last(resampling_stretch)+1
            if next_solid_start < length(record_kmers)
                selected_path_kmers = record_kmers[next_solid_start:end]
                append!(new_record_kmers, selected_path_kmers)
                selected_kmer_qualities = record_kmer_quality_scores[next_solid_start:end]
                append!(new_record_kmer_qualities, selected_kmer_qualities)
            end
        end
    end
    
    for (a, b) in zip(new_record_kmers[1:end-1], new_record_kmers[2:end])
        @assert a != b
    end
    new_record_sequence = Mycelia.kmer_path_to_sequence(new_record_kmers)
    new_record_quality_scores = new_record_kmer_qualities[1]
    for new_record_kmer_quality in new_record_kmer_qualities[2:end]
        push!(new_record_quality_scores, last(new_record_kmer_quality))
    end
    new_record = fastq_record(identifier=new_record_identifier, sequence=new_record_sequence, quality_scores=new_record_quality_scores)
    return new_record
end

function polish_fastq(;fastq, k=1)
    kmer_graph = build_directed_kmer_graph(fastq=fastq, k=k)
    assembly_k = MetaGraphs.get_prop(kmer_graph, :assembly_k)
    @info "polishing with k = $(assembly_k)"
    revised_records = []
    yen_k_shortest_paths_and_weights = Dict{Pair{Int, Int}, Vector{Pair{Vector{Int}, Float64}}}()
    ProgressMeter.@showprogress for record in collect(Mycelia.open_fastx(fastq))
        revised_record = process_fastq_record(;record, kmer_graph, yen_k_shortest_paths_and_weights)
        push!(revised_records, revised_record)
    end
    
    fastq_out = replace(fastq, Mycelia.FASTQ_REGEX => ".k$(assembly_k).fq")
    open(fastq_out, "w") do io
        fastx_io = FASTX.FASTQ.Writer(io)
        for record in revised_records
            write(fastx_io, record)
        end
        close(fastx_io)
    end
    run(`gzip --force $(fastq_out)`)
    return (fastq = fastq_out * ".gz", k=assembly_k)
end

function build_directed_kmer_graph(;fastq, k=1, plot=false)
    if k == 1
        assembly_k = Mycelia.assess_dnamer_saturation([fastq])
    else
        @assert isodd(k)
        @assert Primes.isprime(k)
        assembly_k = k
    end
    kmer_type = Kmers.DNAKmer{assembly_k}

    # initializing the graph with kmer counts
    kmer_counts = Mycelia.count_kmers(kmer_type, fastq)
    ordered_kmers = collect(keys(kmer_counts))
    total_states = length(ordered_kmers)
    graph = MetaGraphs.MetaDiGraph(total_states)
    MetaGraphs.set_prop!(graph, :assembly_k, assembly_k)
    MetaGraphs.set_prop!(graph, :kmer_counts, kmer_counts)
    MetaGraphs.set_prop!(graph, :total_states, total_states)
    MetaGraphs.set_prop!(graph, :ordered_kmers, ordered_kmers)
    kmer_indices = sort(Dict(kmer => i for (i, kmer) in enumerate(keys(kmer_counts))))
    MetaGraphs.set_prop!(graph, :kmer_indices, kmer_indices)
    canonical_kmer_counts = Mycelia.count_canonical_kmers(kmer_type, fastq)
    MetaGraphs.set_prop!(graph, :canonical_kmer_counts, canonical_kmer_counts)
    canonical_kmer_indices = sort(Dict(kmer => i for (i, kmer) in enumerate(keys(canonical_kmer_counts))))
    MetaGraphs.set_prop!(graph, :canonical_kmer_indices, canonical_kmer_indices)
    
    
    # kmer quality and likelihoods
    
    records = collect(Mycelia.open_fastx(fastq))
    read_quality_scores = [collect(FASTX.quality_scores(record)) for record in records]
    all_kmer_quality_support = Dict{kmer_type, Vector{Float64}}()
    for record in records
        record_quality_scores = collect(FASTX.quality_scores(record))
        record_quality_score_slices = [record_quality_scores[i:i+assembly_k-1] for i in 1:length(record_quality_scores)-assembly_k+1]
        sequence = BioSequences.LongDNA{2}(FASTX.sequence(record))
        for ((i, kmer), kmer_base_qualities) in zip(Kmers.EveryKmer{kmer_type}(sequence), record_quality_score_slices)
            if haskey(all_kmer_quality_support, kmer)
                all_kmer_quality_support[kmer] = all_kmer_quality_support[kmer] .+ kmer_base_qualities
            else
                all_kmer_quality_support[kmer] = kmer_base_qualities
            end
        end
    end
    
    # strand normalization shares observational quality across strands - only relevant for non-stranded DNA genome assembly
    strand_normalized_quality_support = Dict{kmer_type, Vector{Float64}}()
    for (kmer, support) in all_kmer_quality_support
        strand_normalized_quality_support[kmer] = support
        if haskey(all_kmer_quality_support, BioSequences.reverse_complement(kmer))
            strand_normalized_quality_support[kmer] .+= all_kmer_quality_support[BioSequences.reverse_complement(kmer)]
        end
    end
    strand_normalized_quality_support
    kmer_mean_quality = sort(Dict(kmer => strand_normalized_quality_support[kmer] ./ canonical_kmer_counts[BioSequences.canonical(kmer)] for kmer in ordered_kmers))
    MetaGraphs.set_prop!(graph, :kmer_mean_quality, kmer_mean_quality)
    kmer_total_quality = sort(Dict(kmer => sum(quality_values) for (kmer, quality_values) in strand_normalized_quality_support))
    MetaGraphs.set_prop!(graph, :kmer_total_quality, kmer_total_quality)
    state_likelihoods = sort(Dict(kmer => total_quality / sum(values(kmer_total_quality)) for (kmer, total_quality) in kmer_total_quality))
    MetaGraphs.set_prop!(graph, :state_likelihoods, state_likelihoods)


    # all transition likelihood calculation
    transition_likelihoods = SparseArrays.spzeros(total_states, total_states)
    for record in records
        sequence = BioSequences.LongDNA{4}(FASTX.sequence(record))
        sources = Kmers.EveryKmer{kmer_type}(sequence[1:end-1])
        destinations = Kmers.EveryKmer{kmer_type}(sequence[2:end])
        for ((source_i, source), (destination_i, destination)) in zip(sources, destinations)
            source_index = kmer_indices[source]
            destination_index = kmer_indices[destination]
            transition_likelihoods[source_index, destination_index] += 1
        end
    end
    for source in 1:total_states
        outgoing_transition_counts = transition_likelihoods[source, :]
        if sum(outgoing_transition_counts) > 0
            transition_likelihoods[source, :] .= transition_likelihoods[source, :] ./ sum(transition_likelihoods[source, :]) 
        end
    end
    row_indices, column_indices, cell_values = SparseArrays.findnz(transition_likelihoods)
    for (row, col, value) in zip(row_indices, column_indices, cell_values)
        Graphs.add_edge!(graph, row, col)
        MetaGraphs.set_prop!(graph, row, col, :transition_likelihood, value)
    end
    MetaGraphs.set_prop!(graph, :transition_likelihoods, transition_likelihoods)

    # helpful for downstream processing
    unbranching_nodes = Set(Int[])
    for node in Graphs.vertices(graph)
        if (Graphs.indegree(graph, node) <= 1) && (Graphs.outdegree(graph, node) <= 1)
            push!(unbranching_nodes, node)
        end
    end
    branching_nodes = Set(setdiff(Graphs.vertices(graph), unbranching_nodes))
    MetaGraphs.set_prop!(graph, :unbranching_nodes, unbranching_nodes)
    MetaGraphs.set_prop!(graph, :branching_nodes, branching_nodes)
    
    
    # total_strand_normalized_quality_support = sum.(collect(values(strand_normalized_quality_support)))
    mean_total_support = Statistics.mean(collect(values(kmer_total_quality)))
    sorted_kmer_total_quality_values = collect(values(kmer_total_quality))
    mean_quality_value = Statistics.mean(sorted_kmer_total_quality_values)
    threshold = mean_quality_value

    xs = [
        [i for (i, y) in enumerate(sorted_kmer_total_quality_values) if y > threshold],
        [i for (i, y) in enumerate(sorted_kmer_total_quality_values) if y <= threshold]
        ]
    
    likely_valid_kmer_indices = xs[1]
    MetaGraphs.set_prop!(graph, :likely_valid_kmer_indices, likely_valid_kmer_indices)
    likely_sequencing_artifact_indices = xs[2]
    MetaGraphs.set_prop!(graph, :likely_sequencing_artifact_indices, likely_sequencing_artifact_indices)
    # likely_sequencing_artifact_kmers = Set(ordered_kmers[likely_sequencing_artifact_indices])
    # likely_valid_kmers = Set(ordered_kmers[likely_valid_kmer_indices])
    # kmer_to_index_map = Dict(kmer => i for (i, kmer) in enumerate(ordered_kmers))
    
    
    if plot
        ys = [
            [y for y in sorted_kmer_total_quality_values if y > threshold],
            [y for y in sorted_kmer_total_quality_values if y <= threshold]
        ]

        p = StatsPlots.scatter(
            xs,
            ys,
            title = "kmer qualities",
            ylabel = "canonical kmer cumulative QUAL value",
            label = ["above" "below"],
            legend = :outertopright,
            # size = (900, 500),
            margins=10StatsPlots.Plots.PlotMeasures.mm,
            xticks = false
        )
        p = StatsPlots.hline!(p, [mean_quality_value], label="mean")
        display(p)
    end
    return graph
end

# selected after trialing previous and next ks and finding those to be too unstable
function iterative_polishing(fastq, max_k = 89, plot=false)
    # initial polishing
    polishing_results = [polish_fastq(fastq=fastq)]
    while (!ismissing(last(polishing_results).k)) && (last(polishing_results).k < max_k)
        next_k = first(filter(k -> k > last(polishing_results).k, Mycelia.ks()))
        # @show next_k
        push!(polishing_results, polish_fastq(fastq=last(polishing_results).fastq, k=next_k))
    end
    return polishing_results
end

function export_blast_db_taxonomy_table(;path_to_db, outfile = path_to_db * ".seqid2taxid.txt.gz")
    Mycelia.add_bioconda_env("blast")
    if !isfile(outfile)
        # -long_seqids adds GI identifiers - these are cross-referenceable through other means so I'm dropping
        @time run(pipeline(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd  -entry all -outfmt "%a %T" -db $(path_to_db)`, `gzip`), outfile))
    else
        @info "$(outfile) already present"
    end
    return outfile
end

function load_blast_db_taxonomy_table(compressed_blast_db_taxonomy_table_file)
    return CSV.read(CodecZlib.GzipDecompressorStream(open(compressed_blast_db_taxonomy_table_file)), delim=' ', header=["sequence_id", "taxid"], DataFrames.DataFrame)
    # data, header = uCSV.read(CodecZlib.GzipDecompressorStream(open(compressed_blast_db_taxonomy_table_file)), delim=' ')
    # header = ["sequence_id", "taxid"]
    # DataFrames.DataFrame(data, header)
end

function parse_xam_to_mapped_records_table(xam)
    record_table = DataFrames.DataFrame(
        template = String[],
        flag = UInt16[],
        reference = String[],
        position = UnitRange{Int}[],
        mappingquality = UInt8[],
        tlen = Int[],
        alignlength = Int[],
        ismapped = Bool[],
        isprimary = Bool[],
        alignment_score = Int[],
        mismatches = Int[]
    )
    if occursin(r"\.bam$", xam)
        MODULE = XAM.BAM
        io = open(xam)
    elseif occursin(r"\.sam$", xam)
        MODULE = XAM.SAM
        io = open(xam)
    elseif occursin(r"\.sam.gz$", xam)
        MODULE = XAM.SAM
        io = CodecZlib.GzipDecompressorStream(open(xam))
    else
        error("unrecognized file extension in file: $xam")
    end
    # filter out header lines
    reader = MODULE.Reader(IOBuffer(join(Iterators.filter(line -> !startswith(line, '@'), eachline(io)), '\n')))
    # reader = MODULE.Reader(io)
    for record in reader
        if XAM.SAM.ismapped(record)
            row = (
                template = XAM.SAM.tempname(record),
                flag = XAM.flag(record),
                reference = XAM.SAM.refname(record),
                position = XAM.SAM.position(record):XAM.SAM.rightposition(record),
                mappingquality = XAM.SAM.mappingquality(record),
                tlen = XAM.SAM.templength(record),
                alignlength = XAM.SAM.alignlength(record),
                ismapped = XAM.SAM.ismapped(record),
                isprimary = XAM.SAM.isprimary(record),
                alignment_score = record["AS"],
                mismatches = record["NM"]
                )
            push!(record_table, row, promote=true)
        end
    end
    close(io)
    return record_table
end

function assess_assembly_quality(;assembly, observations, ks::Vector{Int}=filter(x -> 11 <= x <= 53, Mycelia.ks()))
    results = DataFrames.DataFrame()
    ProgressMeter.@showprogress for k in ks
        observed_canonical_kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, observations)
        assembled_canonical_kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, assembly)
        cosine_distance = kmer_counts_to_cosine_similarity(observed_canonical_kmer_counts, assembled_canonical_kmer_counts)
        js_divergence = kmer_counts_to_js_divergence(observed_canonical_kmer_counts, assembled_canonical_kmer_counts)
        qv = kmer_counts_to_merqury_qv(raw_data_counts=observed_canonical_kmer_counts, assembly_counts=assembled_canonical_kmer_counts)
        push!(results, (;k, cosine_distance, js_divergence, qv))
    end
    return results
end

# function assess_assembly_quality(;assembled_sequence::BioSequences.LongDNA{2}, fastq::String, k::Int)
#     assess_assembly_quality(assembled_sequence=assembled_sequence, fastq=fastq, ks=[k])
# end

function kmer_counts_to_cosine_similarity(kmer_counts_1, kmer_counts_2)
    sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
    a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    # Distances.cosine_dist(a, b) == Distances.cosine_dist(b, a) == Distances.cosine_dist(a ./ sum(a), b ./ sum(b)) == Distances.cosine_dist(b ./ sum(b), a ./ sum(a))
    return Distances.cosine_dist(a, b)
end

function kmer_counts_to_js_divergence(kmer_counts_1, kmer_counts_2)
    sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
    a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    a_norm = a ./ sum(a)
    b_norm = b ./ sum(b)
    # Distances.js_divergence(a ./ sum(a), b ./ sum(b)) == Distances.js_divergence(b ./ sum(b), a ./ sum(a))
    # Distances.js_divergence(a, b) != Distances.js_divergence(a ./ sum(a), b ./ sum(b))
    return Distances.js_divergence(a_norm, b_norm)
end

function jaccard_similarity(set1, set2)
    return length(intersect(set1, set2)) / length(union(set1, set2))
end

function jaccard_distance(set1, set2)
    return 1.0 - jaccard_similarity(set1, set2)
end

# function kmer_counts_to_jaccard(kmer_counts_1::AbstractDict{Kmers.DNAKmer{k}, Int64}, kmer_counts_2::AbstractDict{Kmers.DNAKmer{k}, Int64}) where k
#     # sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
#     # a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
#     # b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
#     # a_indices = findall(a .> 0)
#     # b_indices = findall(b .> 0)
#     # return Distances.jaccard(a_indices, b_indices)
#     return jaccard(collect(keys(kmer_counts_1)), collect(keys(kmer_counts_2)))
# end

function kmer_counts_to_merqury_qv(;raw_data_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}, assembly_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}) where {k,N}
    # Ktotal = # of kmers found in assembly
    Ktotal = length(keys(assembly_counts))
    # Kshared = # of shared kmers between assembly and readset
    Kshared = length(intersect(keys(raw_data_counts), keys(assembly_counts)))
    # probabilitiy_base_in_assembly_correct
    P = (Kshared/Ktotal)^(1/k)
    # # Error rate
    E = 1-P
    QV = -10log10(E)
    # return (;P, E, QV)
    return QV
end

function system_mem_to_minimap_index_size(;system_mem_gb, denominator=6)
    # smaller, higher diversity databases do better with >=5 as the denominator - w/ <=4 they run out of memory
    # denominator = 5 # produced OOM for NT on NERSC
    # denominator = 10 was only 56% efficient for NT on NERSC
    value = Int(floor(system_mem_gb/denominator))
    # 4G is the default
    # this value should be larger for larger memory machines, and smaller for smaller ones
    # it seems related to the total size of the sequences stored in memory, rather than the total size of the in-memory database
    return "$(value)G"
end
# system_mem_to_minimap_index_size(Mycelia.NERSC_MEM)

"""
Run this on the machine you intend to use to map the reads to confirm the index will fit
"""
function minimap_index(;fasta, mem_gb, mapping_type, threads, as_string=false, denominator=6)
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
    if as_string
        cmd = "$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -d $(index_file) $(fasta)"
    else
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -d $(index_file) $(fasta)`
    end
    outfile = index_file
    return (;cmd, outfile)
end

"""
aligning and compressing. No sorting or filtering.

Use shell_only=true to get string command to submit to SLURM
"""
function minimap_map(;
        fasta,
        mem_gb,
        mapping_type,
        threads,
        fastq,
        as_string=false,
        denominator=6
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    temp_sam_outfile = fastq * "." * basename(fasta) * "." * "minimap2.sam"
    outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
        cmd = [map, compress]
    end
    return (;cmd, outfile)
end


"""
aligning and compressing. No sorting or filtering.

Use shell_only=true to get string command to submit to SLURM
"""
function minimap_map_with_index(;
        fasta,
        mem_gb,
        mapping_type,
        threads,
        fastq,
        as_string=false,
        denominator=6
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
    @show index_file
    @assert isfile(index_file)
    temp_sam_outfile = fastq * "." * basename(index_file) * "." * "minimap2.sam"
    # outfile = temp_sam_outfile
    outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
        cmd = [map, compress]
    end
    return (;cmd, outfile)
end

function bam_to_fastq(;bam, fastq=bam * ".fq.gz")
    Mycelia.add_bioconda_env("samtools")
    bam_to_fastq_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools fastq $(bam)`
    gzip_cmd = `gzip`
    p = pipeline(bam_to_fastq_cmd, gzip_cmd)
    if !isfile(fastq)
        @time run(pipeline(p, fastq))
    else
        @info "$(fastq) already exists"
    end
    return fastq
end

function normalize_countmap(countmap)
    sum_total = sum(values(countmap))
    return Dict(k => v/sum_total for (k, v) in countmap)
end

function extract_pacbiosample_information(xml)
    xml_dict = XMLDict.parse_xml(read(xml, String))
    wellsample = xml_dict["ExperimentContainer"]["Runs"]["Run"]["Outputs"]["SubreadSets"]["SubreadSet"]["DataSetMetadata"]["Collections"]["CollectionMetadata"]["WellSample"]

    # Initialize empty arrays to store the data
    biosample_names = []
    barcode_names = []
    
    if haskey(wellsample, "BioSamples")
        # display(wellsample)
        for bs in wellsample["BioSamples"]["BioSample"]
            push!(biosample_names, bs[:Name])
            push!(barcode_names, bs["DNABarcodes"]["DNABarcode"][:Name])
        end
    end

    # Create the DataFrame
    df = DataFrames.DataFrame(BioSampleName=biosample_names, BarcodeName=barcode_names)
    return df
end


# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't recusively import this file
all_other_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_other_julia_files
    include(f)
end

end # module
