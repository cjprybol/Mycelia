import Test
import Mycelia

Test.@testset "Annotation Parsing" begin
    Test.@testset "PGAP taxonomy report parsing" begin
        mktempdir() do dir
            report_path = joinpath(dir, "ani-tax-report.txt")
            open(report_path, "w") do io
                println(io, "ANI report for assembly: assembly_1")
                println(io, "Submitted organism: Escherichia coli (taxid = 562, rank = species, lineage = Bacteria;Proteobacteria)")
                println(io, "Predicted organism: Escherichia coli (taxid = 562, rank = species, lineage = Bacteria;Proteobacteria)")
                println(io, "Best match: Escherichia coli (taxid = 562, rank = species, lineage = Bacteria;Proteobacteria)")
                println(io, "Submitted organism has type: Yes")
                println(io, "Status: CONFIRMED")
                println(io, "Confidence: HIGH")
                println(io, "ANI (Coverages)")
                println(io, "--------------------------------")
                println(io, "95.5 (98.0 97.0) 0 0 1 OK Escherichia coli (GCF_000005845.2, Escherichia coli)")
                println(io, "Table legend:")
            end

            report = Mycelia.parse_pgap_taxonomy_report(report_path)
            Test.@test report.assembly_name == "assembly_1"
            Test.@test report.submitted_taxid == 562
            Test.@test report.status == "CONFIRMED"
            Test.@test report.confidence == "HIGH"
            Test.@test length(report.matches) == 1

            df = Mycelia.load_pgap_taxonomy_reports([report_path])
            Test.@test Mycelia.DataFrames.nrow(df) == 1
            Test.@test df[1, "top_match_organism"] == "Escherichia coli"
        end
    end

    Test.@testset "TransTerm parsing" begin
        mktempdir() do dir
            output_path = joinpath(dir, "transterm.out")
            open(output_path, "w") do io
                println(io, "SEQUENCE seq1")
                println(io, "TERM 1  10 - 20  +  F  99  -12.7 -4.0 |bidir")
            end

            parsed = Mycelia.parse_transterm_output(output_path)
            Test.@test Mycelia.DataFrames.nrow(parsed) == 1
            Test.@test parsed[1, "chromosome"] == "seq1"
            Test.@test parsed[1, "term_id"] == "TERM 1"

            gff_path = Mycelia.transterm_output_to_gff(output_path)
            Test.@test isfile(gff_path)
            gff_contents = read(gff_path, String)
            Test.@test occursin("#seqid", gff_contents)
            Test.@test occursin("terminator", gff_contents)
        end
    end

    Test.@testset "VirSorter score parsing" begin
        mktempdir() do dir
            tsv_path = joinpath(dir, "virsorter.tsv")
            open(tsv_path, "w") do io
                println(io, "seqname\tscore")
                println(io, "contig_1\t0.85")
            end

            df = Mycelia.parse_virsorter_score_tsv(tsv_path)
            Test.@test Mycelia.DataFrames.nrow(df) == 1
            Test.@test "seqname" in Mycelia.DataFrames.names(df)
        end
    end

    Test.@testset "tRNA attribute formatting" begin
        df = Mycelia.DataFrames.DataFrame(
            attributes = [
                "ID=foo-tRNA13-fMetCAT",
                "ID=bar-tRNA01-UndetNNN",
                "ID=baz-tRNA02-AlaGGC",
            ]
        )
        updated = Mycelia.add_trna_attributes(df)
        Test.@test occursin("tRNA-fMet", updated[1, "attributes"])
        Test.@test occursin("tRNA-Undet", updated[2, "attributes"])
        Test.@test occursin("tRNA-Ala", updated[3, "attributes"])
    end
end
