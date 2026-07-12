import Test
import Mycelia

Test.@testset "Stage-2 transition scoring and variant ranking" begin
    index = Mycelia.Rhizomorph.TransitionLikelihoodIndex(
        2,
        Dict(
            ("AA", "AC") => log2(0.9),
            ("AC", "CC") => log2(0.9),
            ("CC", "CG") => log2(0.9),
            ("CG", "GG") => log2(0.9),
            ("AA", "AT") => log2(0.5),
            ("AT", "TC") => log2(0.5),
            ("TC", "CG") => log2(0.5),
            ("AA", "AG") => log2(0.2),
            ("AG", "GC") => log2(0.2),
            ("GC", "CG") => log2(0.2)
        )
    )

    high = Mycelia.Rhizomorph._score_variant_sequence("AACCGG", index)
    medium = Mycelia.Rhizomorph._score_variant_sequence("AATCGG", index)
    unsupported = Mycelia.Rhizomorph._score_variant_sequence("TTTTTT", index)
    Test.@test high.scored_fraction == 1.0
    Test.@test medium.scored_fraction == 1.0
    Test.@test high.mean_log2_probability > medium.mean_log2_probability
    Test.@test unsupported.scored_fraction == 0.0
    Test.@test unsupported.mean_log2_probability == -Inf

    mktempdir() do dir
        gfa = joinpath(dir, "strain_unitigs.gfa")
        info = joinpath(dir, "phased_unitig_info_table.csv")
        write(gfa,
            "H\tVN:Z:1.0\n" *
            "S\tprimary\tAACCGG\n" *
            "S\tsecondary\tAATCGG\n" *
            "S\ttertiary\tAAGCGG\n" *
            "S\trejected\tTTTTTT\n")
        write(info,
            "unitig,coverage\n" *
            "primary,10\n" *
            "secondary,30\n" *
            "tertiary,20\n" *
            "rejected,100\n")

        variants = Mycelia.Rhizomorph._rank_stage2_variants(
            gfa, info, index; max_variants = 3, min_scored_fraction = 0.9)
        Test.@test [variant.id for variant in variants] ==
                   ["primary", "secondary", "tertiary"]
        Test.@test [variant.role for variant in variants] ==
                   [:primary, :secondary, :tertiary]
        Test.@test all(variant -> variant.scored_transition_fraction == 1.0, variants)

        outputs = Mycelia.Rhizomorph._write_ranked_variants(variants, dir)
        Test.@test isfile(outputs.primary)
        Test.@test isfile(outputs.fasta)
        Test.@test isfile(outputs.tsv)
        Test.@test occursin("role=secondary", read(outputs.fasta, String))
        Test.@test count(==('>'), read(outputs.primary, String)) == 1
        Test.@test countlines(outputs.tsv) == 4

        error_message = try
            Mycelia.Rhizomorph._rank_stage2_variants(
                gfa, info, index; max_variants = 4, min_scored_fraction = 0.9)
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin("only 3 graph-supported variants", error_message)
    end
end

Test.@testset "Stage-2 AssemblyConfig validation" begin
    config = Mycelia.Rhizomorph.AssemblyConfig(;
        k = 31,
        corrector = :iterative,
        strategy = :scalable,
        sequencing_tech = :nanopore,
        layout = :olc,
        olc_tool = :metaflye,
        output_dir = "stage2-output")
    Test.@test config.layout == :olc
    Test.@test config.olc_tool == :metaflye
    Test.@test Mycelia.Rhizomorph._resolve_olc_tool(config) == :metaflye

    error_message = try
        Mycelia.Rhizomorph.AssemblyConfig(;
            k = 31,
            corrector = :none,
            layout = :olc,
            sequencing_tech = :nanopore)
        ""
    catch error
        sprint(showerror, error)
    end
    Test.@test occursin("layout=:olc requires corrector=:iterative", error_message)
end

Test.@testset "Strainy fail-loud argument validation" begin
    both_refs_message = try
        Mycelia.run_strainy(;
            gfa_ref = "graph.gfa",
            fasta_ref = "reference.fasta",
            fastq = "reads.fastq")
        ""
    catch error
        sprint(showerror, error)
    end
    Test.@test occursin("exactly one of gfa_ref or fasta_ref", both_refs_message)

    mode_message = try
        Mycelia.run_strainy(;
            gfa_ref = "graph.gfa",
            fastq = "reads.fastq",
            read_mode = :illumina)
        ""
    catch error
        sprint(showerror, error)
    end
    Test.@test occursin("read_mode must be :nano or :hifi", mode_message)
end
