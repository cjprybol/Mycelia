# """
# Downloads and unpacks the desired .tar.gz prebuilt kraken index

# Go to https://benlangmead.github.io/aws-indexes/k2 and identify the appropriate ".tar.gz" url download
# """
# function download_kraken_index(;url, directory="$(homedir())/workspace/kraken")
#     @assert occursin(r"\.tar\.gz$", url)
#     filename = last(split(url, '/'))
#     output_path = joinpath(directory, filename)
#     # @show url
#     # @show filename
#     # @show output_path
#     mkpath(directory)
#     extracted_directory = replace(basename(output_path), ".tar.gz" => "")
#     if !isdir(extracted_directory)
#         mkpath(extracted_directory)
#     end
#     if isempty(readdir(extracted_directory))
#         # download the file only if needed
#         if !isfile(output_path)
#             download(url, output_path)
#             # run(`wget -P $(directory) $(url)`)
#         end
#         run(`tar -xvzf $(output_path) -C $(extracted_directory)`)
#     end
#     return extracted_directory
# end