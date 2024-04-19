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
import DataFrames
import DataStructures
import Dates
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

function add_bioconda_env(pkg; force=false)
    try
        current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
        if !(pkg in current_environments) || force
            @info "installing conda environment $(pkg)"
            run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
            run(`$(CONDA_RUNNER) clean --all -y`)
        else
            @info "conda environment $(pkg) already present; set force=true to update/re-install"
        end
    catch
        add_bioconda_envs()
        current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
        if !(pkg in current_environments) || force
            @info "installing conda environment $(pkg)"
            run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
            run(`$(CONDA_RUNNER) clean --all -y`)
        else
            @info "conda environment $(pkg) already present; set force=true to update/re-install"
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
    run(submission)
    sleep(1)
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
"""
function nersc_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=pwd(),
        qos::String="shared",
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 4,
        cmd::String,
        constrain::String="cpu"
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
    run(submission)
    sleep(1)
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
        identifier = replace(basename(fasta), r"\.f(na|asta|a)" => ""),
        basedir = identifier,
        mmseqsdb = "$(homedir())/workspace/mmseqs/UniRef50"
    )
    mkpath(basedir)
    f = joinpath(basedir, basename(fasta))
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
    return basedir
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
    run(pipeline(p, out_fastq))
    return out_fastq
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
    jellyfish_buffer_size = parse(Int, first(split(read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish mem --mer-len $(k) --mem $(mem)`, String))))
    # @show jellyfish_buffer_size
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
    if !isfile(outfile)
        run(pipeline(open_cmd, cmd))
    end
    
    temp_fasta = outfile * ".fna"
    fna_counts = temp_fasta * ".gz"
    if !isfile(fna_counts)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish dump --output $(temp_fasta) $(outfile)`)
        run(`gzip $(temp_fasta)`)
    end

    temp_tab = outfile * ".tab"
    tabular_counts = temp_tab * ".gz"
    if !isfile(tabular_counts)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish dump --column --tab --output $(temp_tab) $(outfile)`)
        run(`gzip $(temp_tab)`)
    end
    histogram = outfile * ".histogram"
    if !isfile(histogram)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish histo --output $(histogram) $(outfile)`)
    end
    return (;outfile, fna_counts, tabular_counts, histogram)
end

# conda install -c bioconda kmer-jellyfish

# count, bc, info, stats, histo, dump, merge, query, cite, mem, jf

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

# Usage: jellyfish histo [options] db:path
#  -l, --low=uint64                         Low count value of histogram (1)
#  -h, --high=uint64                        High count value of histogram (10000)
#  -f, --full                               Full histo. Don't skip count 0. (false)
#  -o, --output=string                      Output file

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
        # outfile = replace(temp_sam_outfile, ".sam" => ".sorted.bam"),
        outfile = replace(temp_sam_outfile, ".sam" => ".sorted.sam.gz"),
        threads = Sys.CPU_THREADS,
        # 4G is the default
        # for 512Gb RAM this will ask for 102G of index
        index_chunk_size="$(Int(floor(Sys.total_memory()/5 / 1e9)))G",
        shell_only = false
    )
    @show index_chunk_size
    @show threads
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    if shell_only
        # cmd =
        # """
        # $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -ax map-hifi $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        # && $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile) \\
        # | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -bh -o $(outfile) \\
        # && rm $(temp_sam_outfile)
        # """
        # return cmd
    else
        if !isfile(outfile)
            map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_chunk_size) -ax map-hifi $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
            run(map)
            # note - mapping to NCBI NT has too many sequences and header is invalid BAM spec
            # switch to writing out gzip compressed sam instead
            # p = pipeline(
            #     `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile)`,
            #     `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -bh -o $(outfile)`
            # )
            p = pipeline(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile)`,
                `gzip`
            )
            # run(pipeline(p))
            run(pipeline(p, outfile))
            rm(temp_sam_outfile)
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

# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/taxonomy/taxonomy/
function ncbi_taxon_summary(taxa_id)
    Mycelia.add_bioconda_env("ncbi-datasets")
    p = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets datasets summary taxonomy taxon $(taxa_id) --as-json-lines`,
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets dataformat tsv taxonomy --template tax-summary`
        )
    return DataFrames.DataFrame(uCSV.read(open(p), delim='\t', header=1))
end

# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't recusively import this file
all_other_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_other_julia_files
    include(f)
end

end # module