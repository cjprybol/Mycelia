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
        "--forward_reads"
            help = "forward reads"
            required = true
            arg_type = String
            nargs = '+'
        "--reverse_reads"
            help = "reverse reads"
            arg_type = String
        "--kraken_db"
            help = "path to kraken_db"
            arg_type = String
            required = true
        "--out_directory"
            help = "output directory"
            arg_type = String
            default = ""
        "--threads"
            help = "number of threads to run"
            arg_type = Int
            default = 1
    end

    return ArgParse.parse_args(s)
end

function shared_prefix(str1::AbstractString, str2::AbstractString)
    n = min(length(str1), length(str2))
    prefix = ""
    for i in 1:n
        if str1[i] == str2[i]
            prefix *= str1[i]
        else
            break
        end
    end
    return prefix
end

args = parse_arguments()
forward_reads = args["forward_reads"]
reverse_reads = args["reverse_reads"]
out_directory = args["out_directory"]
kraken_db = args["kraken_db"]
threads = args["threads"]

shared_basename = shared_prefix(basename(first(forward_reads)), basename(reverse_reads))
shared_basename = replace(replace(shared_basename, r"\W+$" => ""), r"_+$" => "")
# @show shared_basename

if isempty(out_directory)
    # find common prefix and then add kraken to it
    all_reads = [forward_reads..., reverse_reads]
    unique_read_directories = unique(dirname.(all_reads))
    @assert length(unique_read_directories) == 1 "please specify out_directory with --out_directory flag"
    out_directory = joinpath(first(unique_read_directories), "$(shared_basename)_kraken")
end
mkpath(out_directory)
@assert isdir(out_directory)

output = joinpath(out_directory, "$(shared_basename).$(basename(kraken_db)).kraken-output.tsv")
report = joinpath(out_directory, "$(shared_basename).$(basename(kraken_db)).kraken-report.tsv")
krona_file = report * ".krona"
krona_html = krona_file * ".html"

conda_environments = filter(!isempty, readlines(`conda env list`))
matching_environments = filter(x -> x == "kraken2", first.(split.(conda_environments)))
if length(matching_environments) == 0
    run(`conda create -n kraken2 -c bioconda kraken2`)
end

if !isfile("$(output).gz") || !isfile(report)
    if !isfile(report)
        if ismissing(reverse_reads) || isempty(reverse_reads)
            forward_reads = join(forward_reads, " ")
            cmd = `conda run --live-stream --no-capture-output -n kraken2 kraken2
                    --report-zero-counts
                    --use-names
                    --threads $(threads)
                    --db $(kraken_db)
                    --output $(output)
                    --report $(report)
                    --gzip-compressed
                    $(forward_reads)
                    `
        else
            @assert length(forward_reads) == 1
            forward_reads = first(forward_reads)
            cmd = `conda run --live-stream --no-capture-output -n kraken2 kraken2
                    --report-zero-counts
                    --use-names
                    --threads $(threads)
                    --db $(kraken_db)
                    --output $(output)
                    --report $(report)
                    --gzip-compressed
                    --paired $(forward_reads) $(reverse_reads)
                    `
        end
        @time run(cmd)
    end
    run(`gzip --best $(output)`)
end

if !isfile(krona_file)
    this_script_path = Base.source_path()
    script_path = joinpath(dirname(dirname(this_script_path)), "bin/kreport2krona.py")
    run(`python $(script_path) -r $(report) -o $(krona_file)`)
end

matching_environments = filter(x -> x == "krona", first.(split.(conda_environments)))
if length(matching_environments) == 0
    run(`conda create -n krona -c bioconda krona`)
    # run(`$HOME/home/cjprybol/.conda/envs/krona/bin/ktUpdateTaxonomy.sh`)
    run(`conda run --live-stream --no-capture-output -n krona ktUpdateTaxonomy.sh`)
    
    # Krona installed.  You still need to manually update the taxonomy
    # databases before Krona can generate taxonomic reports.  The update
    # script is ktUpdateTaxonomy.sh.  The default location for storing
    # taxonomic databases is /home/cjprybol/.conda/envs/krona/opt/krona/taxonomy

    # If you would like the taxonomic data stored elsewhere, simply replace
    # this directory with a symlink.  For example:

    # rm -rf /home/cjprybol/.conda/envs/krona/opt/krona/taxonomy
    # mkdir /path/on/big/disk/taxonomy
    # ln -s /path/on/big/disk/taxonomy /home/cjprybol/.conda/envs/krona/opt/krona/taxonomy
    # ktUpdateTaxonomy.sh
    
end

if !isfile(krona_html)
    run(`conda run --live-stream --no-capture-output -n krona ktImportText $(krona_file) -o $(krona_html)`)
end