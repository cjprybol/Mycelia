# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/6_annotation/annotation.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/6_annotation/annotation.jl", "test/6_annotation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Annotation & feature extraction tests

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import FASTX
import BioSequences
import CSV
import DataFrames

# T4 phage
# result = Mycelia.download_genome_by_accession(accession="NC_000866")

Test.@testset "Annotation & feature extraction" begin

    Test.@testset "generate_transterm_coordinates_from_fasta" begin
        mktempdir() do dir
            fasta_file = joinpath(dir, "sample.fna")
            record = FASTX.FASTA.Record("seq1", BioSequences.LongDNA{4}("ATGCGT"))
            Mycelia.write_fasta(outfile=fasta_file, records=[record])
            coords_file = Mycelia.generate_transterm_coordinates_from_fasta(fasta_file)
            Test.@test isfile(coords_file)
            rows = [split(strip(line)) for line in readlines(coords_file) if !isempty(strip(line))]
            Test.@test length(rows) == 2
            rows_by_gene = Dict(row[1] => row for row in rows)
            Test.@test haskey(rows_by_gene, "seq1_start")
            Test.@test haskey(rows_by_gene, "seq1_stop")
            Test.@test rows_by_gene["seq1_start"][2:4] == ["1", "2", "seq1"]
            Test.@test rows_by_gene["seq1_stop"][2:4] == [string(length("ATGCGT") - 1), string(length("ATGCGT")), "seq1"]
        end
    end

    Test.@testset "generate_transterm_coordinates_from_gff" begin
        mktempdir() do dir
            gff_file = joinpath(dir, "sample.gff")
            open(gff_file, "w") do io
                println(io, "seq1\tsource\tgene\t1\t5\t0\t+\t0\tID=gene1")
            end
            coords_file = Mycelia.generate_transterm_coordinates_from_gff(gff_file)
            Test.@test isfile(coords_file)
            rows = [split(strip(line)) for line in readlines(coords_file) if !isempty(strip(line))]
            Test.@test length(rows) == 1
            row = rows[1]
            Test.@test row[1] == "gene1"
            Test.@test row[2:4] == ["0", "5", "seq1"]
        end
    end

    Test.@testset "read_gff and split attributes" begin
        mktempdir() do dir
            gff_file = joinpath(dir, "sample.gff")
            open(gff_file, "w") do io
                println(io, "seq1\tsource\tgene\t1\t5\t0\t+\t0\tID=gene1;Note=demo")
            end
            gff_df = Mycelia.read_gff(gff_file)
            Test.@test size(gff_df, 1) == 1
            Mycelia.split_gff_attributes_into_columns(gff_df)
            Test.@test gff_df.ID[1] == "gene1"
            Test.@test gff_df.Note[1] == "demo"
        end
    end

    Test.@testset "update_gff_with_mmseqs" begin
        mktempdir() do dir
            gff_file = joinpath(dir, "sample.gff")
            open(gff_file, "w") do io
                println(io, "seq1\tsrc\tgene\t1\t5\t0\t+\t0\tID=gene1")
            end
            mmseqs_file = joinpath(dir, "sample.tsv")
            df = DataFrames.DataFrame(query=["gene1"], qheader=["gene1 # ID=gene1"], qlen=[5],
                target=["t"], tlen=[5], theader=["desc protein"], nident=[5], fident=[1.0], pident=[100.0],
                alnlen=[5], mismatch=[0], gapopen=[0], qstart=[1], qend=[5], tstart=[1], tend=[5],
                evalue=[1e-5], bits=[50.0], taxid=[1], taxname=["taxon"])
            CSV.write(mmseqs_file, df; delim='\t')
            updated = Mycelia.update_gff_with_mmseqs(gff_file, mmseqs_file)
            expected = "label=\"desc__protein\";product=\"desc__protein\";ID=gene1"
            Test.@test updated[1, "attributes"] == expected
        end
    end

end
