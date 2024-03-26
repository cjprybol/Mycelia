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

const MAMBA = joinpath(Conda.BINDIR, "mamba")

function add_bioconda_env(pkg; force=false)
    try
        current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(Conda.conda) env list`))))))
        if !(pkg in current_environments) || force
            @info "installing conda environment $(pkg)"
            run(`$(MAMBA) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
            run(`$(MAMBA) clean --all -y`)
        else
            @info "conda environment $(pkg) already present; set force=true to update/re-install"
        end
    catch
        add_bioconda_envs()
        current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(Conda.conda) env list`))))))
        if !(pkg in current_environments) || force
            @info "installing conda environment $(pkg)"
            run(`$(MAMBA) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
            run(`$(MAMBA) clean --all -y`)
        else
            @info "conda environment $(pkg) already present; set force=true to update/re-install"
        end
    end
end

function add_bioconda_envs(;all=false, force=false)
    if !isfile(MAMBA)
        Conda.add("mamba")
    end
    if !isfile(joinpath(Conda.BINDIR, "pigz"))
        run(`$(MAMBA) install pigz -y`)
    end
    current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(Conda.conda) env list`))))))
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
    run(`$(MAMBA) clean --all -y`)
end

"""
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
Submit a command to SLURM using sbatch
"""
function nersc_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=pwd(),
        qos::String,
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 2,
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
    p = pipeline(`$(Mycelia.MAMBA) run --live-stream -n gfatools gfatools gfa2fa $(gfa)`, fasta)
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
        run(pipeline(`$(Mycelia.MAMBA) run --live-stream -n bedtools bedtools genomecov -d -ibam $(bam)`, genome_coverage_file))
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
    run(`$(Mycelia.MAMBA) run --live-stream -n mmseqs2 mmseqs easy-search
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
    run(`$(Mycelia.MAMBA) run --live-stream -n mmseqs2 mmseqs createdb $(fasta) $(fasta)_DB`)
    run(`$(Mycelia.MAMBA) run --live-stream -n mmseqs2 mmseqs createindex --search-type 3 $(fasta)_DB $(tempdir())`)
    run(`$(Mycelia.MAMBA) run --live-stream -n mmseqs2 mmseqs easy-linclust $(fasta)_DB $(output) $(tmp)`)
    run(`$(Mycelia.MAMBA) run --live-stream -n mmseqs2 mmseqs createtsv $(fasta)_DB $(fasta)_DB $(output) $(output).tsv`)
    return "$(output).tsv"
end

function mmseqs_easy_cluster(;fasta, output=fasta*".mmseqs_easy_cluster", tmp=tempdir())
    Mycelia.add_bioconda_env("mmseqs2")
    run(`Mycelia.MAMBA run --live-stream -n mmseqs2 mmseqs easy-cluster $(fasta) $(output) $(tmp) --min-seq-id 0.5 -c 0.8 --cov-mode 1`)
    run(`$(Mycelia.MAMBA) run --live-stream -n mmseqs2 mmseqs createtsv $(fasta) $(fasta) $(output) $(output).tsv`)
    return "$(output).tsv"
end

# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't recusively import this file
all_other_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_other_julia_files
    include(f)
end

end # module