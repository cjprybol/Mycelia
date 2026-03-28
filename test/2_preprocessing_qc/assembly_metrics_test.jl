import Test
import FASTX
import BioSequences
import Mycelia

Test.@testset "assembly_metrics" begin
    Test.@testset "GC content type-safety regression (td-8ba9u)" begin
        # Regression guard: assembly_metrics must use
        # FASTX.sequence(BioSequences.LongDNA{4}, record) so BioSequences
        # symbol comparisons work. If it reverts to FASTX.sequence(record)
        # (returns String), Char comparisons against DNA_G/DNA_C are always
        # false and gc_content is always 0.0 regardless of sequence content.

        tmpdir = mktempdir()
        fasta_path = joinpath(tmpdir, "contigs.fasta")

        # Write two contigs with known GC content:
        #   contig_1: GCGCGCGC  — 100% GC (8 bases, 8 GC)
        #   contig_2: ATATATATAT — 0% GC  (10 bases, 0 GC)
        # Combined: 8 GC out of 18 bases = 44.44...%
        open(fasta_path, "w") do io
            println(io, ">contig_1 length=8")
            println(io, "GCGCGCGC")
            println(io, ">contig_2 length=10")
            println(io, "ATATATATAT")
        end

        metrics = Mycelia.assembly_metrics(fasta_path)

        Test.@test metrics !== nothing
        Test.@test metrics.n_contigs == 2
        Test.@test metrics.total_length == 18

        # Key regression assertion: GC must reflect actual content.
        # If the String-vs-BioSequences bug recurs, this will be 0.0.
        expected_gc = 8 / 18
        Test.@test metrics.gc_content ≈ expected_gc atol = 1e-6
        Test.@test metrics.gc_content > 0.0   # never silently zero
    end

    Test.@testset "all-GC contig" begin
        tmpdir = mktempdir()
        fasta_path = joinpath(tmpdir, "all_gc.fasta")
        open(fasta_path, "w") do io
            println(io, ">contig_1 length=8")
            println(io, "GCGCGCGC")
        end
        metrics = Mycelia.assembly_metrics(fasta_path)
        Test.@test metrics !== nothing
        Test.@test metrics.gc_content ≈ 1.0 atol = 1e-6
        Test.@test metrics.largest_contig == 8
        Test.@test metrics.n50 == 8
    end

    Test.@testset "nothing returned for missing file" begin
        Test.@test Mycelia.assembly_metrics("/nonexistent/path.fasta") === nothing
        Test.@test Mycelia.assembly_metrics(nothing) === nothing
    end
end
