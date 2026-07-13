import Test
import Mycelia

Test.@testset "Stage-2 transition scoring and segment-candidate ranking" begin
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

    high = Mycelia.Rhizomorph._score_segment_candidate_sequence("AACCGG", index)
    medium = Mycelia.Rhizomorph._score_segment_candidate_sequence("AATCGG", index)
    unsupported = Mycelia.Rhizomorph._score_segment_candidate_sequence(
        "TTTTTT", index)
    partially_supported = Mycelia.Rhizomorph._score_segment_candidate_sequence(
        "AACCGT", index)
    ambiguous = Mycelia.Rhizomorph._score_segment_candidate_sequence("AACNGG", index)
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

    Test.@test_throws ArgumentError Mycelia.Rhizomorph.TransitionLikelihoodIndex(
        2, Dict(("AA", "AC") => 0.1))
    Test.@test_throws ArgumentError Mycelia.Rhizomorph.TransitionLikelihoodIndex(
        2, Dict(("AA", "AC") => NaN))
    Test.@test_throws ArgumentError Mycelia.Rhizomorph.TransitionLikelihoodIndex(
        2, Dict(("aa", "ac") => -1.0, ("AA", "AC") => -2.0))

    graph = Mycelia.Rhizomorph.build_kmer_graph(
        [Mycelia.FASTX.FASTA.Record("candidate", "AACCGG")],
        2;
        mode = :doublestrand,
    )
    graph_index = Mycelia.Rhizomorph._transition_likelihood_index(graph, 2)
    graph_score = Mycelia.Rhizomorph._score_segment_candidate_sequence(
        "AACCGG", graph_index)
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

        candidates = Mycelia.Rhizomorph._rank_stage2_segment_candidates(
            gfa, info, index; min_scored_fraction = 0.9)
        Test.@test [candidate.id for candidate in candidates] ==
                   ["primary", "secondary", "tertiary"]
        Test.@test [candidate.likelihood_rank for candidate in candidates] == [1, 2, 3]
        Test.@test [candidate.abundance_rank for candidate in candidates] == [1, 2, 3]
        Test.@test [candidate.role for candidate in candidates] ==
                   [:primary, :secondary, :tertiary]
        Test.@test all(candidate -> candidate.scored_transition_fraction == 1.0,
            candidates)

        outputs = Mycelia.Rhizomorph._write_ranked_variants(
            candidates, dir, nothing)
        Test.@test isfile(outputs.primary)
        Test.@test isfile(outputs.likelihood_fasta)
        Test.@test isfile(outputs.abundance_fasta)
        Test.@test isfile(outputs.tsv)
        Test.@test occursin(
            "likelihood_rank=2", read(outputs.likelihood_fasta, String))
        Test.@test occursin(
            "experimental_segment_candidate",
            read(outputs.likelihood_fasta, String),
        )
        Test.@test count(==('>'), read(outputs.primary, String)) == 1
        Test.@test countlines(outputs.tsv) == 4

        diagnostic_cap = Mycelia.Rhizomorph._rank_stage2_segment_candidates(
            gfa, info, index; max_variants = 2, min_scored_fraction = 0.9)
        Test.@test length(diagnostic_cap) == 2
        oversized_cap = Mycelia.Rhizomorph._rank_stage2_segment_candidates(
            gfa, info, index; max_variants = 4, min_scored_fraction = 0.9)
        Test.@test length(oversized_cap) == 3
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

        candidates = [
            Mycelia.Rhizomorph.RankedSegmentCandidate(
                1, 2, :primary, "likelihood_primary", "AACC", 0.0, 1.0, 10.0),
            Mycelia.Rhizomorph.RankedSegmentCandidate(
                2, 1, :secondary, "abundance_primary", "AACC", -1.0, 1.0, 20.0)
        ]
        write(joinpath(dir, "primary_segment_candidate.fasta"), ">stale\nAACC\n")
        outputs = Mycelia.Rhizomorph._write_ranked_variants(candidates, dir, 3)
        Test.@test outputs.primary_status == :rank_disagreement
        Test.@test outputs.primary === nothing
        Test.@test outputs.likelihood_primary_id == "likelihood_primary"
        Test.@test outputs.abundance_primary_id == "abundance_primary"
        Test.@test !isfile(joinpath(dir, "primary_segment_candidate.fasta"))
        Test.@test occursin("primary_status", read(outputs.tsv, String))
    end
end

Test.@testset "Stage-2 GFA and PAF validation" begin
    mktempdir() do directory
        gfa = joinpath(directory, "segments.gfa")
        invalid_cases = [
            ("S\tonly_id\n", "expected at least 3"),
            ("S\tdup\tAACC\nS\tdup\tAACC\n", "duplicate GFA segment"),
            ("S\tdup\t*\nS\tdup\tAACC\n", "duplicate GFA segment"),
            ("S\tbad\tAACC\tdp:f:-1\n", "finite and nonnegative"),
            ("S\tbad\tAACC\tdp:f:NaN\n", "finite and nonnegative"),
            ("S\tbad\tAACZ\n", "invalid DNA"),
        ]
        for (contents, expected) in invalid_cases
            write(gfa, contents)
            message = try
                Mycelia.Rhizomorph._read_gfa_records(gfa)
                ""
            catch exception
                sprint(showerror, exception)
            end
            Test.@test occursin(expected, message)
        end

        records = Dict(
            "candidate" =>
                Mycelia.Rhizomorph.GfaSegmentRecord("AACC", nothing),
            "candidate_b" =>
                Mycelia.Rhizomorph.GfaSegmentRecord("AATC", nothing),
        )
        valid = "read1\t4\t0\t4\t+\tcandidate\t4\t0\t4\t4\t4\t60\n"
        hits = Mycelia.Rhizomorph._read_stage2_paf_best_hits(
            IOBuffer(valid), records)
        Test.@test hits["read1"].target == "candidate"
        Test.@test Mycelia.Rhizomorph._stage2_support_coverages(
            records, hits)["candidate"] == 1.0

        tied = valid *
               "read1\t4\t0\t4\t+\tcandidate_b\t4\t0\t4\t4\t4\t60\n"
        tied_hits = Mycelia.Rhizomorph._read_stage2_paf_best_hits(
            IOBuffer(tied), records)
        Test.@test tied_hits["read1"].target == "candidate_b"

        invalid_paf_cases = [
            ("read1\t4\t0\n", "expected at least 12"),
            ("read1\t4\t0\t4\t+\tunknown\t4\t0\t4\t4\t4\t60\n",
                "unknown target"),
            ("read1\t4\t0\t4\t+\tcandidate\t4\t0\t4\tNaN\t4\t60\n",
                "matches is not an integer"),
            ("read1\t4\t0\t4\t+\tcandidate\t4\t0\t4\t-1\t4\t60\n",
                "support must be"),
            ("read1\t4\t0\t0\t+\tcandidate\t4\t0\t4\t4\t4\t60\n",
                "invalid query coordinates"),
            ("read1\t4\t0\t4\t+\tcandidate\t4\t0\t0\t4\t4\t60\n",
                "invalid target coordinates"),
            ("read1\t4\t0\t4\t+\tcandidate\t4\t0\t4\t0\t4\t60\n",
                "support must be positive"),
        ]
        for (contents, expected) in invalid_paf_cases
            message = try
                Mycelia.Rhizomorph._read_stage2_paf_best_hits(
                    IOBuffer(contents), records)
                ""
            catch exception
                sprint(showerror, exception)
            end
            Test.@test occursin(expected, message)
        end

        no_hits = Mycelia.Rhizomorph._read_stage2_paf_best_hits(
            IOBuffer(""), records)
        zero_message = try
            Mycelia.Rhizomorph._stage2_support_coverages(records, no_hits)
            ""
        catch exception
            sprint(showerror, exception)
        end
        Test.@test occursin("zero total mapped support", zero_message)
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
    Mycelia.Rhizomorph._score_segment_candidate_sequence(sequence, index)
    allocations = @allocated Mycelia.Rhizomorph._score_segment_candidate_sequence(
        sequence, index)
    Test.@test allocations <= 1_000_000
end

Test.@testset "Stage-2 probabilistic segment-candidate abundance bridge" begin
    candidates = [
        Mycelia.Rhizomorph.RankedSegmentCandidate(
            1, 1, :primary, "primary", "AACC", 0.0, 1.0, 30.0),
        Mycelia.Rhizomorph.RankedSegmentCandidate(
            2, 2, :secondary, "secondary", "AATC", -1.0, 1.0, 20.0),
        Mycelia.Rhizomorph.RankedSegmentCandidate(
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
    result = Mycelia.Rhizomorph.infer_segment_candidate_abundances(
        candidates, alignments, noise)
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
        output = Mycelia.Rhizomorph.write_probabilistic_segment_candidate_ranking(
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

    tied = Mycelia.Rhizomorph.infer_segment_candidate_abundances(
        candidates,
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

    threads_message = try
        Mycelia.run_strainy(;
            gfa_ref = "graph.gfa",
            fastq = "reads.fastq",
            threads = 0)
        ""
    catch error
        sprint(showerror, error)
    end
    Test.@test occursin("threads must be positive", threads_message)

    transform_message = mktempdir() do directory
        try
            Mycelia.run_strainy(;
                gfa_ref = "graph.gfa",
                fastq = "reads.fastq",
                outdir = directory,
                stage = :transform)
            ""
        catch error
            sprint(showerror, error)
        end
    end
    Test.@test occursin("requires a prior phase consensus cache", transform_message)

    phase = Mycelia._strainy_output_paths("phase-output", :phase)
    transform = Mycelia._strainy_output_paths("transform-output", :transform)
    e2e = Mycelia._strainy_output_paths("e2e-output", :e2e)
    Test.@test propertynames(phase) == propertynames(transform) == propertynames(e2e)
    Test.@test phase.alignment_phased_bam isa String
    Test.@test phase.strain_contigs_gfa === nothing
    Test.@test phase.strain_assemblies === nothing
    Test.@test transform.alignment_phased_bam === nothing
    Test.@test transform.consensus_cache isa String
    Test.@test transform.strain_contigs_gfa isa String
    Test.@test transform.experimental_segment_candidates_fasta ==
               transform.strain_assemblies
    Test.@test endswith(
        transform.experimental_segment_candidates_fasta,
        "experimental_segment_candidates.fasta",
    )
    Test.@test e2e.alignment_phased_bam isa String
    Test.@test e2e.strain_contigs_gfa isa String
    Test.@test length(Mycelia._strainy_required_outputs(phase)) == 3
    Test.@test length(Mycelia._strainy_required_outputs(transform)) == 6
    Test.@test transform.consensus_cache in
               Mycelia._strainy_required_outputs(transform)
    Test.@test length(Mycelia._strainy_required_outputs(e2e)) == 8
    Test.@test Mycelia._legacy_strainy_stage("phase") == :e2e
    Test.@test Mycelia._legacy_strainy_stage("transform") == :transform

    mktempdir() do directory
        gfa = joinpath(directory, "strain_contigs.gfa")
        fasta = joinpath(directory, "experimental_segment_candidates.fasta")
        write(gfa, "H\tVN:Z:1.0\nS\tzeta\tAACC\nS\talpha\tAATC\n")
        output = Mycelia._write_strainy_segment_candidates_fasta(gfa, fasta)
        Test.@test output == fasta
        text = read(fasta, String)
        Test.@test occursin("experimental_segment_candidate", text)
        Test.@test first(findfirst("id=zeta", text)) <
                   first(findfirst("id=alpha", text))
        Test.@test count(==('>'), text) == 2

        write(gfa, "S\tdup\tAACC\nS\tdup\tAATC\n")
        Test.@test_throws ErrorException begin
            Mycelia._write_strainy_segment_candidates_fasta(gfa, fasta)
        end

        unsafe_outdir = joinpath(directory, "unsafe ' \$(touch injected)")
        unsafe_outputs = Mycelia._strainy_output_paths(unsafe_outdir, :e2e)
        script = Mycelia._strainy_executor_script("echo safe", unsafe_outputs)
        Test.@test occursin(
            "mkdir -p $(Base.shell_escape(unsafe_outdir))", script)
        Test.@test occursin("seen[\$2]++", script)
        Test.@test occursin("test -s " * Base.shell_escape(
            unsafe_outputs.experimental_segment_candidates_fasta), script)
    end
end
