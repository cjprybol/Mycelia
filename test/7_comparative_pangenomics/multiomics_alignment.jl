# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=test -e 'include("test/7_comparative_pangenomics/multiomics_alignment.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=test -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/multiomics_alignment.jl", "test/7_comparative_pangenomics", execute=false)'
# ````

import Pkg
if isinteractive()
    Pkg.activate("..")
end
using Revise
using Test
import Mycelia
import Random
import FASTX
import Random
import UUIDs
import CodecZlib
import CSV
import DataFrames
import Mycelia
import FASTX
import BioSequences
import Arrow

# Test merging and mapping long-read pacbio data with re-identification
@testset "Test merging and mapping long-read pacbio data with re-identification" begin
    ## simulate a fasta record of 10kb and save it to disk
    original_fasta_record = Mycelia.random_fasta_record(seed=42, L=10_000)
    @test Mycelia.seq2sha256(FASTX.sequence(original_fasta_record)) == "96f36383f772afb5f41db96c42cdaed12b8a6bc151744c108b60a33df7fd56d5"
    fasta_file = "test_fasta.fna.gz"
    @test Mycelia.write_fasta(outfile=fasta_file, records = [original_fasta_record]) == fasta_file
    @test isfile(fasta_file)
    @test filesize(fasta_file) > 0
    
    ## build a minimap2 hifi index
    minimap_index_result = Mycelia.minimap_index(
        fasta = fasta_file,
        mapping_type = "map-hifi"
    )
    if !isfile(minimap_index_result.outfile)
        @time run(minimap_index_result.cmd)
    end
    @test isfile(minimap_index_result.outfile)
    @test filesize(minimap_index_result.outfile) > 0
    
    ## simulate nearly perfect long reads at 10x coverage, 3 times
    fastq_list = [
        Mycelia.simulate_pacbio_reads(fasta = fasta_file, quantity = q) for q in ["5x", "10x", "20x"]
    ]
    @test fastq_list == [
        "test_fasta.badread.pacbio2021.5x.fq.gz",
        "test_fasta.badread.pacbio2021.10x.fq.gz",
        "test_fasta.badread.pacbio2021.20x.fq.gz"
    ]
    @test all(filesize.(fastq_list) .> 0)
    
    ## merge and map asserting all reads align and we can re-identify original files and records
    merged_mapping_results = Mycelia.merge_and_map_single_end_samples(
        fasta_reference=fasta_file,
        fastq_list=fastq_list, 
        minimap_index=minimap_index_result.outfile, 
        mapping_type="map-hifi",
    )
    @test Set(unique(merged_mapping_results.results_table.input_file)) == Set(fastq_list)
    percent_mapped = 
        count(merged_mapping_results.results_table.ismapped) / length(merged_mapping_results.results_table.ismapped)
    @test percent_mapped >= 0.9
    ## @test all(merged_mapping_results.results_table.ismapped)

    ## Cleanup
    @test isfile(fasta_file)
    rm(fasta_file)
    @test isfile(minimap_index_result.outfile)
    rm(minimap_index_result.outfile)

    @show merged_mapping_results.tsv_file
    @test isfile(merged_mapping_results.tsv_file)
    rm(merged_mapping_results.tsv_file)

    foreach(f -> isfile(f) && rm(f), fastq_list)
    @test isfile(merged_mapping_results.joint_fastq_file)
    rm(merged_mapping_results.joint_fastq_file)
    @test isfile(merged_mapping_results.fastq_id_mapping_table)
    rm(merged_mapping_results.fastq_id_mapping_table)
    @test isfile(merged_mapping_results.bam_file)
    rm(merged_mapping_results.bam_file)
    @test isfile(merged_mapping_results.jld2_file)
    rm(merged_mapping_results.jld2_file)
end

## Multi-omics alignment & mapping tests
@testset "Multi-omics alignment & mapping" begin
end