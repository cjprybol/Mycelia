# ORF caller wrapper tests
#
# julia --project=. --color=yes -e 'include("test/6_annotation/orf_callers.jl")'

import Test
import Mycelia

Test.@testset "ORF caller wrappers" begin
    mktempdir() do dir
        fasta_file = joinpath(dir, "contigs.fasta")
        open(fasta_file, "w") do io
            write(io, ">seq1\nATGAAATAG\n")
        end

        out_dir = joinpath(dir, "out")
        mkpath(out_dir)

        prodigal_gff = joinpath(out_dir, "$(basename(fasta_file)).prodigal.gff")
        open(prodigal_gff, "w") do io
            write(io, "")
        end
        prodigal_result = Mycelia.run_prodigal(fasta_file=fasta_file, out_dir=out_dir, translation_table=11)
        Test.@test prodigal_result.gff == prodigal_gff
        Test.@test prodigal_result.faa == joinpath(out_dir, "$(basename(fasta_file)).prodigal.faa")

        pyrodigal_gff = joinpath(out_dir, "$(basename(fasta_file)).pyrodigal.gff")
        pyrodigal_faa = joinpath(out_dir, "$(basename(fasta_file)).pyrodigal.faa")
        open(pyrodigal_gff, "w") do io
            write(io, "")
        end
        open(pyrodigal_faa, "w") do io
            write(io, "")
        end
        pyrodigal_result = Mycelia.run_pyrodigal(fasta_file=fasta_file, out_dir=out_dir, translation_table=11)
        Test.@test pyrodigal_result.gff == pyrodigal_gff
        Test.@test pyrodigal_result.faa == pyrodigal_faa

        prodigal_gv_gff = joinpath(out_dir, "$(basename(fasta_file)).prodigal-gv.gff")
        open(prodigal_gv_gff, "w") do io
            write(io, "")
        end
        prodigal_gv_result = Mycelia.run_prodigal_gv(fasta_file=fasta_file, out_dir=out_dir, translation_table=11)
        Test.@test prodigal_gv_result.gff == prodigal_gv_gff
        Test.@test prodigal_gv_result.faa == joinpath(out_dir, "$(basename(fasta_file)).prodigal-gv.faa")

        augustus_gff = joinpath(out_dir, "$(basename(fasta_file)).augustus.gff")
        open(augustus_gff, "w") do io
            write(io, "")
        end
        augustus_result = Mycelia.run_augustus(fasta_file=fasta_file, out_dir=out_dir, species="human")
        Test.@test augustus_result.gff == augustus_gff
        Test.@test augustus_result.out == joinpath(out_dir, "$(basename(fasta_file)).augustus.out")

        metaeuk_prefix = joinpath(out_dir, "$(basename(fasta_file)).metaeuk")
        metaeuk_gff = "$(metaeuk_prefix)_preds.gff"
        open(metaeuk_gff, "w") do io
            write(io, "")
        end
        metaeuk_result = Mycelia.run_metaeuk(fasta_file=fasta_file, db_file="dummy.db", out_dir=out_dir)
        Test.@test metaeuk_result.gff == metaeuk_gff
        Test.@test metaeuk_result.out_prefix == metaeuk_prefix
    end
end
