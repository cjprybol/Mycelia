# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/8_tool_integration/minimap_merge_map_split.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/minimap_merge_map_split.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import FASTX
import CodecZlib

Test.@testset "Minimap merge/map/split helpers" begin
    Test.@testset "prefix_fastq_reads adds unique tag and mapping" begin
        temp_fastq = tempname() * ".fq"
        write(temp_fastq, "@r1\nACGT\n+\nIIII\n@r1\nTGCA\n+\nHHHH\n")
        out_fastq = tempname() * ".fq.gz"
        map_tsv = tempname() * ".tsv.gz"

        res = Mycelia.prefix_fastq_reads(
            temp_fastq;
            sample_tag="sampleA",
            out_fastq=out_fastq,
            mapping_out=map_tsv,
            mapping_format=:tsv,
            id_delimiter="::"
        )

        # identifiers carry the tag and delimiter
        reader = Mycelia.open_fastx(res.out_fastq)
        ids = String[]
        for rec in reader
            push!(ids, String(FASTX.identifier(rec)))
        end
        close(reader)
        Test.@test all(startswith(id, "sampleA::") for id in ids)

        io = CodecZlib.GzipDecompressorStream(open(map_tsv))
        tsv_lines = readlines(io)
        close(io)
        Test.@test length(tsv_lines) == 3 # header + 2 reads
    end

    Test.@testset "Command construction without execution" begin
        ref = tempname() * ".fa"
        write(ref, ">ref\nACGT\n")
        fq1 = tempname() * ".fq"
        fq2 = tempname() * ".fq"
        write(fq1, "@a\nAAAA\n+\nIIII\n")
        write(fq2, "@b\nCCCC\n+\nIIII\n")

        pair1 = tempname() * ".fq"
        pair2 = tempname() * ".fq"
        write(pair1, "@read/1\nGGGG\n+\nIIII\n")
        write(pair2, "@read/2\nGGGG\n+\nIIII\n")

        res = Mycelia.minimap_merge_map_and_split(
            reference_fasta=ref,
            mapping_type="sr",
            single_end_fastqs=[fq1, fq2],
            paired_end_fastqs=[(pair1, pair2)],
            run_mapping=false,
            run_splitting=false,
            write_read_map=false,
            minimap_extra_args=["-N", "5"],
            as_string=true
        )

        Test.@test occursin("minimap2", res.minimap_cmd)
        Test.@test occursin("-x sr", res.minimap_cmd)
        Test.@test length(res.prefixed_fastqs) == 4
        Test.@test length(res.sample_outputs) == 3
        Test.@test occursin(basename(ref), res.merged_bam)
    end

    Test.@testset "Paired identifiers stay linked after prefixing" begin
        r1 = tempname() * ".fq"
        r2 = tempname() * ".fq"
        write(r1, "@read1/1\nAAAA\n+\nIIII\n")
        write(r2, "@read1/2\nAAAA\n+\nIIII\n")
        out1 = tempname() * ".fq.gz"
        out2 = tempname() * ".fq.gz"

        Mycelia.prefix_fastq_reads(r1; sample_tag="sampleX", out_fastq=out1, id_delimiter="::")
        Mycelia.prefix_fastq_reads(r2; sample_tag="sampleX", out_fastq=out2, id_delimiter="::")

        rec1_reader = Mycelia.open_fastx(out1)
        rec2_reader = Mycelia.open_fastx(out2)
        rec1 = first(rec1_reader)
        rec2 = first(rec2_reader)

        base1 = replace(String(FASTX.identifier(rec1)), "/1" => "")
        base2 = replace(String(FASTX.identifier(rec2)), "/2" => "")
        Test.@test base1 == base2
        close(rec1_reader)
        close(rec2_reader)
    end
end
