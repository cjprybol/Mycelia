import Pkg
Pkg.activate("..")
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
    fasta_file == "test_fasta.fna.gz"
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
    @test all(merged_mapping_results.results_table.ismapped)
end
    
## Multi-omics alignment & mapping tests
@testset "Multi-omics alignment & mapping" begin
end
