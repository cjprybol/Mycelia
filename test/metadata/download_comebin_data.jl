import Mycelia

const COMEBIN_TEST_FILE_ID = "1xWpN2z8JTaAzWW4TcOl0Lr4Y_x--Fs5s"

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

    if isempty(filter(x -> endswith(lowercase(x), ".bam"), _find_files_by_extension(output_dir, [".bam"])))
        @info "Extracting COMEBin test data..."
        run(`unzip -o $(zip_path) -d $(output_dir)`)
    end

    bam_files = _find_files_by_extension(output_dir, [".bam"])
    isempty(bam_files) && error("No BAM files found after extraction")

    fasta_candidates = _find_files_by_extension(output_dir, [".fasta", ".fa", ".fna"])
    contigs = _choose_contigs_file(fasta_candidates)

    bam_dirs = unique(dirname.(bam_files))
    length(bam_dirs) == 1 || error("Expected BAM files in a single directory, found: $(bam_dirs)")
    bam_path = bam_dirs[1]
    return (
        contigs = contigs,
        bam_path = bam_path,
        bam_files = bam_files,
        root = output_dir
    )
end
