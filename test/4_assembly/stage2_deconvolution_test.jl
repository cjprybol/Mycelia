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
    partially_supported = Mycelia.Rhizomorph._score_variant_sequence("AACCGT", index)
    ambiguous = Mycelia.Rhizomorph._score_variant_sequence("AACNGG", index)
    Test.@test keytype(index.log2_probability) == UInt64
    Test.@test high.scored_fraction == 1.0
    Test.@test medium.scored_fraction == 1.0
    Test.@test high.mean_log2_probability > medium.mean_log2_probability
    Test.@test unsupported.scored_fraction == 0.0
    Test.@test unsupported.mean_log2_probability == -Inf
    Test.@test partially_supported.scored_fraction == 0.75
    Test.@test partially_supported.mean_log2_probability == -Inf
    Test.@test ambiguous.scored_fraction < 1.0
    Test.@test ambiguous.mean_log2_probability == -Inf

    graph = Mycelia.Rhizomorph.build_kmer_graph(
        [Mycelia.FASTX.FASTA.Record("candidate", "AACCGG")],
        2;
        mode = :doublestrand,
    )
    graph_index = Mycelia.Rhizomorph._transition_likelihood_index(graph, 2)
    graph_score = Mycelia.Rhizomorph._score_variant_sequence("AACCGG", graph_index)
    Test.@test graph_score.scored_fraction == 1.0
    Test.@test isfinite(graph_score.mean_log2_probability)

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
            "primary,30\n" *
            "secondary,20\n" *
            "tertiary,10\n" *
            "rejected,100\n")

        variants = Mycelia.Rhizomorph._rank_stage2_variants(
            gfa, info, index; max_variants = 3, min_scored_fraction = 0.9)
        Test.@test [variant.id for variant in variants] ==
                   ["primary", "secondary", "tertiary"]
        Test.@test [variant.likelihood_rank for variant in variants] == [1, 2, 3]
        Test.@test [variant.abundance_rank for variant in variants] == [1, 2, 3]
        Test.@test [variant.role for variant in variants] ==
                   [:primary, :secondary, :tertiary]
        Test.@test all(variant -> variant.scored_transition_fraction == 1.0, variants)

        outputs = Mycelia.Rhizomorph._write_ranked_variants(variants, dir, 3)
        Test.@test isfile(outputs.primary)
        Test.@test isfile(outputs.likelihood_fasta)
        Test.@test isfile(outputs.abundance_fasta)
        Test.@test isfile(outputs.tsv)
        Test.@test occursin(
            "likelihood_rank=2", read(outputs.likelihood_fasta, String))
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

Test.@testset "Stage-2 GFA coverage and primary agreement" begin
    index = Mycelia.Rhizomorph.TransitionLikelihoodIndex(
        2, Dict(("AA", "AC") => 0.0, ("AC", "CC") => 0.0))
    mktempdir() do dir
        gfa = joinpath(dir, "strain_contigs.gfa")
        info = joinpath(dir, "phased_unitig_info_table.csv")
        write(gfa,
            "S\tlikelihood_primary\tAACC\tdp:i:10\n" *
            "S\tabundance_primary\tAACC\tdp:i:20\n")
        write(info, "unitig,coverage\n")
        records = Mycelia.Rhizomorph._read_gfa_records(gfa)
        Test.@test records["likelihood_primary"].coverage == 10.0
        Test.@test records["abundance_primary"].coverage == 20.0

        variants = [
            Mycelia.Rhizomorph.RankedVariant(
                1, 2, :primary, "likelihood_primary", "AACC", 0.0, 1.0, 10.0),
            Mycelia.Rhizomorph.RankedVariant(
                2, 1, :secondary, "abundance_primary", "AACC", -1.0, 1.0, 20.0)
        ]
        outputs = Mycelia.Rhizomorph._write_ranked_variants(variants, dir, 3)
        Test.@test outputs.primary_status == :rank_disagreement
        Test.@test outputs.likelihood_primary_id == "likelihood_primary"
        Test.@test outputs.abundance_primary_id == "abundance_primary"
        Test.@test occursin(
            "status=rank_disagreement", read(outputs.primary, String))
        Test.@test occursin("primary_status", read(outputs.tsv, String))
    end
end

Test.@testset "Stage-2 compact transition scoring allocations" begin
    cycle = repeat("ACGT", 16)
    index = Mycelia.Rhizomorph.TransitionLikelihoodIndex(
        31,
        Dict(
            (String(cycle[index:(index + 30)]),
                String(cycle[(index + 1):(index + 31)])) => -0.1
            for index in 1:4
        ),
    )
    sequence = repeat("ACGT", 10_000)
    Test.@test isconcretetype(keytype(index.vertex_ids))
    Test.@test keytype(index.log2_probability) == UInt64
    Mycelia.Rhizomorph._score_variant_sequence(sequence, index)
    allocations = @allocated Mycelia.Rhizomorph._score_variant_sequence(
        sequence, index)
    Test.@test allocations <= 1_000_000
end

Test.@testset "Stage-2 probabilistic abundance bridge" begin
    variants = [
        Mycelia.Rhizomorph.RankedVariant(
            1, 1, :primary, "primary", "AACC", 0.0, 1.0, 30.0),
        Mycelia.Rhizomorph.RankedVariant(
            2, 2, :secondary, "secondary", "AATC", -1.0, 1.0, 20.0),
        Mycelia.Rhizomorph.RankedVariant(
            3, 3, :tertiary, "tertiary", "AAGC", -2.0, 1.0, 10.0)
    ]
    alignments = Mycelia.Rhizomorph.ReadCandidateLikelihood[
        Mycelia.Rhizomorph.ReadCandidateLikelihood("r1", "primary", -0.01),
        Mycelia.Rhizomorph.ReadCandidateLikelihood("r1", "secondary", -4.0),
        Mycelia.Rhizomorph.ReadCandidateLikelihood("r2", "primary", -0.02),
        Mycelia.Rhizomorph.ReadCandidateLikelihood("r2", "secondary", -3.0),
        Mycelia.Rhizomorph.ReadCandidateLikelihood("r3", "secondary", -0.01),
        Mycelia.Rhizomorph.ReadCandidateLikelihood("r3", "primary", -4.0),
        Mycelia.Rhizomorph.ReadCandidateLikelihood("r4", "tertiary", -0.01),
        Mycelia.Rhizomorph.ReadCandidateLikelihood("r4", "secondary", -3.0)
    ]
    noise = [
        Mycelia.Rhizomorph.ReadNoiseLikelihood("r$(index)", -8.0)
        for index in 1:4
    ]
    result = Mycelia.Rhizomorph.infer_variant_abundances(
        variants, alignments, noise)
    Test.@test result.inference.converged
    Test.@test result.primary_agreement
    Test.@test result.primary_status == :resolved
    Test.@test length(result.candidate_set_sha256) == 64
    Test.@test length(result.inference_input_sha256) == 64
    Test.@test [ranking.candidate_id for ranking in result.rankings] ==
               ["primary", "secondary", "tertiary"]
    Test.@test isapprox(sum(ranking.abundance for ranking in result.rankings) +
                        result.inference.noise_abundance, 1.0; atol = 1.0e-12)
    Test.@test all(ranking -> ranking.has_alignment_support, result.rankings)

    mktempdir() do directory
        output = Mycelia.Rhizomorph.write_probabilistic_variant_ranking(
            joinpath(directory, "posterior.tsv"), result)
        Test.@test isfile(output)
        table_text = read(output, String)
        Test.@test occursin("estimated_abundance", table_text)
        Test.@test occursin("noise_abundance", table_text)
        Test.@test occursin("final_log_likelihood", table_text)
        Test.@test occursin("trace_sha256", table_text)
        Test.@test occursin("inference_input_sha256", table_text)
        Test.@test countlines(output) == 4
    end

    tied = Mycelia.Rhizomorph.infer_variant_abundances(
        variants,
        [
            Mycelia.Rhizomorph.ReadCandidateLikelihood("tie", "primary", 0.0),
            Mycelia.Rhizomorph.ReadCandidateLikelihood("tie", "secondary", 0.0)
        ],
        [Mycelia.Rhizomorph.ReadNoiseLikelihood("tie", -10.0)],
    )
    Test.@test !tied.primary_agreement
    Test.@test tied.primary_status == :tied

    collision_left = Mycelia.Rhizomorph._stage2_inference_input_sha256(
        Mycelia.Rhizomorph.ReadCandidateLikelihood[],
        [Mycelia.Rhizomorph.ReadNoiseLikelihood("a\0-1.0\nnoise\0b", -2.0)],
    )
    collision_right = Mycelia.Rhizomorph._stage2_inference_input_sha256(
        Mycelia.Rhizomorph.ReadCandidateLikelihood[],
        [
            Mycelia.Rhizomorph.ReadNoiseLikelihood("a", -1.0),
            Mycelia.Rhizomorph.ReadNoiseLikelihood("b", -2.0),
        ],
    )
    Test.@test collision_left != collision_right
end

Test.@testset "Stage-2 AssemblyConfig validation" begin
    config = Mycelia.Rhizomorph.AssemblyConfig(;
        k = 31,
        corrector = :iterative,
        strategy = :scalable,
        sequencing_tech = :nanopore,
        output_dir = "stage2-output")
    Test.@test config.corrector == :iterative
    Test.@test config.strategy == :scalable
    Test.@test config.sequencing_tech == :nanopore
    Test.@test config.output_dir == "stage2-output"

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
