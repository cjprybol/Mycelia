# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/metadata/download_comebin_data.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/metadata/download_comebin_data.jl", "test/metadata", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Mycelia

const COMEBIN_TEST_FILE_ID = "1xWpN2z8JTaAzWW4TcOl0Lr4Y_x--Fs5s"
const COMEBIN_SAMPLE_FRACTION = 0.01
const COMEBIN_SAMPLE_SEED = 42

function _ensure_gdown(utils_env::String)
    if !isfile(Mycelia.CONDA_RUNNER)
        error("Conda runner not found at $(Mycelia.CONDA_RUNNER)")
    end
    if !Mycelia._check_conda_env_exists(utils_env)
        run(`$(Mycelia.CONDA_RUNNER) create -y -n $(utils_env) python=3.10`)
    end
    try
        run(`$(Mycelia.CONDA_RUNNER) run -n $(utils_env) python -c "import gdown"`)
    catch
        run(`$(Mycelia.CONDA_RUNNER) run -n $(utils_env) pip install gdown`)
    end
    return utils_env
end

function _ensure_samtools()
    Mycelia.add_bioconda_env("samtools")
    return "samtools"
end

function _find_files_by_extension(root::String, extensions::Vector{String})
    matches = String[]
    for (dirpath, _, filenames) in walkdir(root)
        for name in filenames
            lower = lowercase(name)
            if any(ext -> endswith(lower, ext), extensions)
                push!(matches, joinpath(dirpath, name))
            end
        end
    end
    return matches
end

function _choose_contigs_file(candidates::Vector{String})
    isempty(candidates) && error("No contigs FASTA found in extracted dataset")
    preferred = filter(path -> occursin(r"(contig|scaffold)", lowercase(basename(path))), candidates)
    return isempty(preferred) ? first(candidates) : first(preferred)
end

function _select_primary_bam(bam_files::Vector{String})
    isempty(bam_files) && error("No BAM files available for selection")
    sizes = map(filesize, bam_files)
    _, idx = findmax(sizes)
    return bam_files[idx]
end

function _downsample_bam(bam::String, output_bam::String; fraction::Float64=COMEBIN_SAMPLE_FRACTION, seed::Int=COMEBIN_SAMPLE_SEED, threads::Int=Mycelia.get_default_threads())
    if isfile(output_bam) && filesize(output_bam) > 0 && isfile(output_bam * ".bai")
        return output_bam
    end
    _ensure_samtools()
    mkpath(dirname(output_bam))
    fraction_str = replace(string(fraction), "0." => "")
    sample_arg = string(seed, ".", fraction_str)
    @info "Downsampling COMEBin BAM for tests" bam=bam output_bam=output_bam sample=sample_arg

    view_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -b -s $(sample_arg) $(bam)`
    sort_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -o $(output_bam) -`
    run(pipeline(view_cmd, sort_cmd))
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools index $(output_bam)`)
    return output_bam
end

"""
    download_and_prep_comebin_data(output_dir::String; force_download=false)

Download the COMEBin test dataset from Google Drive and extract it.
"""
function download_and_prep_comebin_data(output_dir::String; force_download::Bool=false)
    mkpath(output_dir)
    utils_env = _ensure_gdown("mycelia_utils")

    zip_path = joinpath(output_dir, "comebin_test_data.zip")
    if force_download || !isfile(zip_path)
        @info "Downloading COMEBin test data..."
        run(`$(Mycelia.CONDA_RUNNER) run -n $(utils_env) gdown --id $(COMEBIN_TEST_FILE_ID) -O $(zip_path)`)
    end

    if !isfile(zip_path)
        error("COMEBin test data download failed at $(zip_path)")
    end

    if isempty(_find_files_by_extension(output_dir, [".bam"])) ||
            isempty(_find_files_by_extension(output_dir, [".fasta", ".fa", ".fna"]))
        @info "Extracting COMEBin test data..."
        run(`unzip -o $(zip_path) -d $(output_dir)`)
    end

    bam_files = _find_files_by_extension(output_dir, [".bam"])
    isempty(bam_files) && error("No BAM files found after extraction")
    primary_bam = _select_primary_bam(bam_files)
    small_bam_dir = joinpath(dirname(primary_bam), "small_bams")
    downsampled_bam = joinpath(small_bam_dir, "comebin_test_data_small.bam")
    downsampled_bam = _downsample_bam(primary_bam, downsampled_bam)

    fasta_candidates = _find_files_by_extension(output_dir, [".fasta", ".fa", ".fna"])
    contigs = _choose_contigs_file(fasta_candidates)

    bam_path = small_bam_dir
    return (
        contigs = contigs,
        bam_path = bam_path,
        bam_files = [downsampled_bam],
        root = output_dir
    )
end
