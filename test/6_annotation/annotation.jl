# Annotation & feature extraction tests

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import FASTX
import BioSequences
import DataFrames
import CSV
import uCSV

Test.@testset "Annotation & feature extraction" begin

    Test.@testset "generate_transterm_coordinates_from_fasta" begin
        mktempdir() do dir
            fasta_file = joinpath(dir, "sample.fna")
            record = FASTX.FASTA.Record("seq1", BioSequences.DNASequence("ATGCGT"))
            Mycelia.write_fasta(outfile=fasta_file, records=[record])
            coords_file = Mycelia.generate_transterm_coordinates_from_fasta(fasta_file)
            Test.@test isfile(coords_file)
            coords_df = uCSV.read(coords_file, DataFrames.DataFrame; delim="  ")
            Test.@test size(coords_df, 1) == 2
            Test.@test coords_df.gene_id == ["seq1_start", "seq1_stop"]
            Test.@test coords_df.start == [1, length("ATGCGT")-1]
            Test.@test coords_df.stop == [2, length("ATGCGT")]
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
            df = uCSV.read(coords_file, DataFrames.DataFrame; delim="  ")
            Test.@test size(df, 1) == 1
            Test.@test df.gene_id[1] == "gene1"
            Test.@test df.start[1] == 0
            Test.@test df.stop[1] == 5
            Test.@test df."#seqid"[1] == "seq1"
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
