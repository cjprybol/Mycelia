import Test
import Mycelia
import StableRNGs

const COMMUNITY_ILLUMINA_CALLS = Ref{Vector{NamedTuple}}(NamedTuple[])
const COMMUNITY_LONG_READ_CALLS = Ref{Vector{NamedTuple}}(NamedTuple[])
const COMMUNITY_CONCAT_CALLS = Ref{Vector{NamedTuple}}(NamedTuple[])

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end

    function simulate_illumina_reads(; fasta::String,
            coverage::Union{Nothing, Number} = nothing,
            read_count::Union{Nothing, Int} = nothing,
            outbase::String = "",
            read_length::Int = 150,
            mflen::Int = 500,
            sdev::Int = 10,
            seqSys::String = "HS25",
            amplicon::Bool = false,
            errfree::Bool = false,
            paired::Bool = true,
            rndSeed::Int = current_unix_datetime(),
            quiet::Bool = true
    )
        push!(
            Main.COMMUNITY_ILLUMINA_CALLS[],
            (
                fasta = fasta,
                coverage = coverage,
                read_count = read_count,
                outbase = outbase,
                read_length = read_length,
                mflen = mflen,
                sdev = sdev,
                seqSys = seqSys,
                amplicon = amplicon,
                errfree = errfree,
                paired = paired,
                rndSeed = rndSeed,
                quiet = quiet
            )
        )
        if read_count === nothing && coverage === nothing
            error("Either `coverage` or `read_count` must be provided.")
        end
        if isempty(outbase)
            if read_count !== nothing
                outbase = "$(fasta).rcount_$(read_count).art"
            else
                outbase = "$(fasta).fcov_$(coverage)x.art"
            end
        end
        mkpath(dirname(outbase))
        open(joinpath(dirname(outbase), "wrapper_calls.tsv"), "w") do io
            println(io, join([
                seqSys,
                string(read_length),
                paired ? "1" : "0",
                coverage === nothing ? "" : string(coverage),
                read_count === nothing ? "" : string(read_count)
            ], '\t'))
        end
        forward = paired ? outbase * "1.fq.gz" : outbase * ".fq.gz"
        write(forward, errfree ? "truth-forward" : "forward")
        reverse = nothing
        if paired
            reverse = outbase * "2.fq.gz"
            write(reverse, errfree ? "truth-reverse" : "reverse")
        end
        sam = outbase * ".sam.gz"
        write(sam, "sam")
        error_free_sam = nothing
        if errfree
            error_free_sam = outbase * "_errFree.sam.gz"
            write(error_free_sam, "error-free-sam")
        end
        return (
            forward_reads = forward,
            reverse_reads = reverse,
            sam = sam,
            error_free_sam = error_free_sam
        )
    end

    function simulate_nanopore_reads(; fasta,
            quantity,
            outfile::String = "",
            quiet = true,
            seed = current_unix_datetime())
        if isempty(outfile)
            outfile = replace(fasta, Mycelia.FASTA_REGEX => ".badread.nanopore_r10.$(quantity).fq.gz")
        end
        push!(
            Main.COMMUNITY_LONG_READ_CALLS[],
            (platform = :nanopore, fasta = fasta, quantity = quantity, outfile = outfile, quiet = quiet, seed = seed)
        )
        mkpath(dirname(outfile))
        write(outfile, "nanopore")
        return outfile
    end

    function simulate_pacbio_reads(; fasta,
            quantity,
            outfile::String = "",
            quiet = true,
            seed = current_unix_datetime())
        if isempty(outfile)
            outfile = replace(fasta, Mycelia.FASTA_REGEX => ".badread.pacbio_hifi.$(quantity).fq.gz")
        end
        push!(
            Main.COMMUNITY_LONG_READ_CALLS[],
            (platform = :pacbio_hifi, fasta = fasta, quantity = quantity, outfile = outfile, quiet = quiet, seed = seed)
        )
        mkpath(dirname(outfile))
        write(outfile, "pacbio")
        return outfile
    end

    function concatenate_fastq_files(; fastq_files, output_fastq, gzip = true, force = false)
        push!(
            Main.COMMUNITY_CONCAT_CALLS[],
            (fastq_files = collect(fastq_files), output_fastq = output_fastq, gzip = gzip, force = force)
        )
        mkpath(dirname(output_fastq))
        write(output_fastq, join(fastq_files, "\n"))
        return output_fastq
    end
end

function write_reference_fasta(path::AbstractString)
    write(path, ">seq1\nACGTACGT\n>seq2\nTTTTGG\n>seq3\nCCCCAAAAGG\n")
    return path
end

function reset_stubbed_simulation_state!()
    empty!(COMMUNITY_ILLUMINA_CALLS[])
    empty!(COMMUNITY_LONG_READ_CALLS[])
    empty!(COMMUNITY_CONCAT_CALLS[])
    return nothing
end

Test.@testset "Simulate Metagenome Community" begin
    Test.@testset "custom selection with retained intermediates" begin
        mktempdir() do dir
            reference_fasta = write_reference_fasta(joinpath(dir, "reference.fna"))
            reference_table = Mycelia.DataFrames.DataFrame(
                sequence_id = ["seq1", "seq2", "seq3"],
                length = [8, 6, 10],
                label = ["alpha", "beta", "gamma"]
            )
            outdir = joinpath(dir, "community_custom")

            result = Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = reference_table,
                n_organisms = 2,
                depth_target = 3.0,
                abundance_profile = :custom,
                readset = :illumina_pe150,
                outdir = outdir,
                selected_ids = ["seq1", "seq3"],
                weights = [0.25, 0.75],
                rng = StableRNGs.StableRNG(21),
                replicate = 4,
                run_simulation = false,
                cleanup = false
            )

            truth_table = result.truth_table
            expected_relative_bases = [8 * 1.5, 10 * 4.5]
            expected_relative_bases ./= sum(expected_relative_bases)

            Test.@test truth_table.sequence_id == ["seq1", "seq3"]
            Test.@test truth_table.selection_order == [1, 2]
            Test.@test truth_table.replicate == [4, 4]
            Test.@test all(isapprox.(truth_table.weight, [0.25, 0.75]; atol = 1e-12))
            Test.@test all(isapprox.(truth_table.coverage, [1.5, 4.5]; atol = 1e-12))
            Test.@test all(isapprox.(
                truth_table.expected_relative_bases,
                expected_relative_bases;
                atol = 1e-12
            ))

            Test.@test result.reads.forward === nothing
            Test.@test result.reads.reverse === nothing
            Test.@test result.truth_reads.forward === nothing
            Test.@test result.truth_reads.reverse === nothing
            Test.@test Mycelia.DataFrames.nrow(result.per_genome_reads) == 0
            Test.@test sort(collect(keys(result.per_genome_fastas))) == ["seq1", "seq3"]
            Test.@test occursin(">seq1", read(result.per_genome_fastas["seq1"], String))
            Test.@test occursin("CCCCAAAAGG", read(result.per_genome_fastas["seq3"], String))
            Test.@test isdir(joinpath(outdir, "tmp"))
        end
    end

    Test.@testset "accession inputs with automatic sampling and cleanup" begin
        mktempdir() do dir
            reference_fasta = write_reference_fasta(joinpath(dir, "reference.fna"))
            reference_table = Mycelia.DataFrames.DataFrame(
                accession = Union{String, Missing}["seq1", missing, "seq3"],
                length = [8, 99, 10],
                label = ["alpha", "missing", "gamma"]
            )
            outdir = joinpath(dir, "community_accession")

            result = Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = reference_table,
                n_organisms = 2,
                depth_target = 5.0,
                abundance_profile = :equal,
                readset = :nanopore,
                outdir = outdir,
                rng = StableRNGs.StableRNG(7),
                replicate = 2,
                run_simulation = false,
                cleanup = true
            )

            Test.@test Set(String.(result.truth_table.sequence_id)) == Set(["seq1", "seq3"])
            Test.@test all(isapprox.(result.truth_table.weight, [0.5, 0.5]; atol = 1e-12))
            Test.@test all(isapprox.(result.truth_table.coverage, [5.0, 5.0]; atol = 1e-12))
            Test.@test !isdir(joinpath(outdir, "tmp"))
        end
    end

    Test.@testset "validation failures" begin
        mktempdir() do dir
            reference_fasta = write_reference_fasta(joinpath(dir, "reference.fna"))
            valid_table = Mycelia.DataFrames.DataFrame(
                sequence_id = ["seq1", "seq2"],
                length = [8, 6]
            )

            Test.@test_throws AssertionError Mycelia.simulate_metagenome_community(
                reference_fasta = joinpath(dir, "missing.fna"),
                reference_table = valid_table,
                n_organisms = 1,
                depth_target = 1.0,
                abundance_profile = :equal,
                readset = :illumina_pe150,
                outdir = joinpath(dir, "missing_reference")
            )

            Test.@test_throws ErrorException Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = Mycelia.DataFrames.DataFrame(name = ["seq1"], length = [8]),
                n_organisms = 1,
                depth_target = 1.0,
                abundance_profile = :equal,
                readset = :illumina_pe150,
                outdir = joinpath(dir, "missing_id_column")
            )

            Test.@test_throws AssertionError Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = valid_table,
                n_organisms = 3,
                depth_target = 1.0,
                abundance_profile = :equal,
                readset = :illumina_pe150,
                outdir = joinpath(dir, "too_many_organisms")
            )

            Test.@test_throws ErrorException Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = valid_table,
                n_organisms = 1,
                depth_target = 1.0,
                abundance_profile = :custom,
                readset = :illumina_pe150,
                outdir = joinpath(dir, "missing_custom_weights")
            )

            Test.@test_throws AssertionError Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = valid_table,
                n_organisms = 1,
                depth_target = 1.0,
                abundance_profile = :equal,
                readset = :illumina_pe150,
                outdir = joinpath(dir, "bad_selected_ids"),
                selected_ids = ["seq1", "seq2"]
            )
        end
    end

    Test.@testset "missing fasta ids are rejected" begin
        mktempdir() do dir
            reference_fasta = joinpath(dir, "partial_reference.fna")
            write(reference_fasta, ">seq1\nACGTACGT\n")
            reference_table = Mycelia.DataFrames.DataFrame(
                sequence_id = ["seq1", "seq3"],
                length = [8, 10]
            )

            Test.@test_throws AssertionError Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = reference_table,
                n_organisms = 1,
                depth_target = 2.0,
                abundance_profile = :custom,
                readset = :illumina_se250,
                outdir = joinpath(dir, "missing_sequence"),
                selected_ids = ["seq3"],
                weights = [1.0],
                run_simulation = false
            )
        end
    end

    Test.@testset "stubbed read simulation branches" begin
        reset_stubbed_simulation_state!()
        mktempdir() do dir
            reference_fasta = write_reference_fasta(joinpath(dir, "reference.fna"))
            reference_table = Mycelia.DataFrames.DataFrame(
                sequence_id = ["seq1", "seq2"],
                length = [8, 6]
            )
            outdir = joinpath(dir, "illumina_paired")

            result = Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = reference_table,
                n_organisms = 1,
                depth_target = 4.0,
                abundance_profile = :equal,
                readset = :illumina_pe150,
                outdir = outdir,
                selected_ids = ["seq1"],
                rng = StableRNGs.StableRNG(4),
                replicate = 3,
                run_simulation = true,
                emit_truth_reads = true,
                cleanup = true
            )

            Test.@test result.reads.forward == joinpath(outdir, "reads", "community.illumina_pe150.r3.R1.fq.gz")
            Test.@test result.reads.reverse == joinpath(outdir, "reads", "community.illumina_pe150.r3.R2.fq.gz")
            Test.@test result.truth_reads.forward == joinpath(outdir, "truth_reads", "community.illumina_pe150.r3.truth.R1.fq.gz")
            Test.@test result.truth_reads.reverse == joinpath(outdir, "truth_reads", "community.illumina_pe150.r3.truth.R2.fq.gz")
            Test.@test Mycelia.DataFrames.nrow(result.per_genome_reads) == 1
            Test.@test result.per_genome_reads[1, :forward_reads] == joinpath(outdir, "tmp", "seq1.illumina_pe150.r31.fq.gz")
            Test.@test result.per_genome_reads[1, :reverse_reads] == joinpath(outdir, "tmp", "seq1.illumina_pe150.r32.fq.gz")
            Test.@test result.per_genome_reads[1, :truth_forward_reads] == joinpath(outdir, "tmp", "seq1.illumina_pe150.r3.errfree1.fq.gz")
            Test.@test result.per_genome_reads[1, :truth_reverse_reads] == joinpath(outdir, "tmp", "seq1.illumina_pe150.r3.errfree2.fq.gz")
            Test.@test result.per_genome_reads[1, :truth_sam] == joinpath(outdir, "tmp", "seq1.illumina_pe150.r3.errfree_errFree.sam.gz")
            Test.@test !isdir(joinpath(outdir, "tmp"))
            Test.@test length(COMMUNITY_ILLUMINA_CALLS[]) == 2
            Test.@test COMMUNITY_ILLUMINA_CALLS[][1].paired
            Test.@test !COMMUNITY_ILLUMINA_CALLS[][1].errfree
            Test.@test COMMUNITY_ILLUMINA_CALLS[][1].seqSys == "HS25"
            Test.@test COMMUNITY_ILLUMINA_CALLS[][1].read_length == 150
            Test.@test COMMUNITY_ILLUMINA_CALLS[][2].errfree
            Test.@test length(COMMUNITY_CONCAT_CALLS[]) == 4
        end

        reset_stubbed_simulation_state!()
        mktempdir() do dir
            reference_fasta = write_reference_fasta(joinpath(dir, "reference.fna"))
            reference_table = Mycelia.DataFrames.DataFrame(
                accession = ["seq1", "seq2"],
                length = [8, 6]
            )
            outdir = joinpath(dir, "illumina_single")

            result = Mycelia.simulate_metagenome_community(
                reference_fasta = reference_fasta,
                reference_table = reference_table,
                n_organisms = 1,
                depth_target = 5.0,
                abundance_profile = :custom,
                readset = :illumina_se250,
                outdir = outdir,
                selected_ids = ["seq2"],
                weights = [2.0],
                rng = StableRNGs.StableRNG(8),
                replicate = 5,
                run_simulation = true,
                emit_truth_reads = false,
                cleanup = false
            )

            Test.@test result.reads.forward == joinpath(outdir, "reads", "community.illumina_se250.r5.R1.fq.gz")
            Test.@test result.reads.reverse === nothing
            Test.@test result.truth_reads.forward === nothing
            Test.@test result.truth_reads.reverse === nothing
            Test.@test Mycelia.DataFrames.nrow(result.per_genome_reads) == 1
            Test.@test result.per_genome_reads[1, :forward_reads] == joinpath(outdir, "tmp", "seq2.illumina_se250.r5.fq.gz")
            Test.@test ismissing(result.per_genome_reads[1, :reverse_reads])
            Test.@test ismissing(result.per_genome_reads[1, :truth_forward_reads])
            Test.@test ismissing(result.per_genome_reads[1, :truth_reverse_reads])
            Test.@test ismissing(result.per_genome_reads[1, :truth_sam])
            Test.@test isdir(joinpath(outdir, "tmp"))
            Test.@test length(COMMUNITY_ILLUMINA_CALLS[]) == 1
            Test.@test !COMMUNITY_ILLUMINA_CALLS[][1].paired
            Test.@test COMMUNITY_ILLUMINA_CALLS[][1].seqSys == "MSv3"
            Test.@test COMMUNITY_ILLUMINA_CALLS[][1].read_length == 250
            Test.@test length(COMMUNITY_CONCAT_CALLS[]) == 1
            Test.@test result.truth_table[1, :weight] == 1.0
        end

        for (readset, platform, suffix) in [
                (:nanopore, :nanopore, "nanopore"),
                (:pacbio_hifi, :pacbio_hifi, "pacbio_hifi")
            ]
            reset_stubbed_simulation_state!()
            mktempdir() do dir
                reference_fasta = write_reference_fasta(joinpath(dir, "reference.fna"))
                reference_table = Mycelia.DataFrames.DataFrame(
                    sequence_id = ["seq3"],
                    length = [10]
                )
                outdir = joinpath(dir, string(readset))

                result = Mycelia.simulate_metagenome_community(
                    reference_fasta = reference_fasta,
                    reference_table = reference_table,
                    n_organisms = 1,
                    depth_target = 6.0,
                    abundance_profile = :equal,
                    readset = readset,
                    outdir = outdir,
                    selected_ids = ["seq3"],
                    rng = StableRNGs.StableRNG(12),
                    replicate = 6,
                    run_simulation = true,
                    cleanup = true
                )

                Test.@test result.reads.forward == joinpath(outdir, "reads", "community.$(readset).r6.R1.fq.gz")
                Test.@test result.reads.reverse === nothing
                Test.@test result.truth_reads.forward === nothing
                Test.@test result.truth_reads.reverse === nothing
                Test.@test Mycelia.DataFrames.nrow(result.per_genome_reads) == 1
                Test.@test result.per_genome_reads[1, :forward_reads] == joinpath(outdir, "tmp", "seq3.$(suffix).r6.fq.gz")
                Test.@test result.per_genome_reads[1, :readset] == string(readset)
                Test.@test length(COMMUNITY_LONG_READ_CALLS[]) == 1
                Test.@test COMMUNITY_LONG_READ_CALLS[][1].platform == platform
                Test.@test COMMUNITY_LONG_READ_CALLS[][1].quantity == "6.0x"
                Test.@test length(COMMUNITY_CONCAT_CALLS[]) == 1
                Test.@test !isdir(joinpath(outdir, "tmp"))
            end
        end
    end
end
