# td-06er: default-CI contract tests for multi-read-set hybrid assembly.
#
# These tests use injected correction and assembler runners, so they exercise
# pairing, routing, provenance, and cleanup without installing external tools.
# The gated companion smoke test owns real third-party execution.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/multi_input_hybrid_assembly_test.jl")'

import CodecZlib
import FASTX
import Mycelia
import Test

if !isdefined(Main, :test_throws_message)
    include(joinpath(dirname(@__DIR__), "test_helpers.jl"))
end

function multi_input_fastq_record(
        identifier::AbstractString,
        sequence::AbstractString = "ACGTACGT",
)::FASTX.FASTQ.Record
    return FASTX.FASTQ.Record(
        String(identifier),
        String(sequence),
        repeat("I", length(sequence)),
    )
end

function multi_input_write_fastq(
        path::AbstractString,
        records::AbstractVector{FASTX.FASTQ.Record},
)::String
    mkpath(dirname(path))
    open(FASTX.FASTQ.Writer, path) do writer
        for record in records
            write(writer, record)
        end
    end
    return String(path)
end

function multi_input_gzip_file(
        source::AbstractString,
        destination::AbstractString,
)::String
    open(destination, "w") do io
        gzip_stream = CodecZlib.GzipCompressorStream(io)
        try
            write(gzip_stream, read(source))
        finally
            close(gzip_stream)
        end
    end
    return String(destination)
end

function multi_input_collect_fastq(reads::Any)::Vector{FASTX.FASTQ.Record}
    values = collect(reads)
    if all(value -> value isa AbstractString, values)
        records = FASTX.FASTQ.Record[]
        for path in values
            reader = Mycelia.open_fastx(String(path))
            try
                append!(records, FASTX.FASTQ.Record[record for record in reader])
            finally
                close(reader)
            end
        end
        return records
    end
    return FASTX.FASTQ.Record[record for record in values]
end

function multi_input_fake_correction_runner(
        calls::Vector{NamedTuple};
        transform::Function = (index, records) -> records,
)::Function
    function correction_runner(
            reads::Any,
            config::Mycelia.Rhizomorph.AssemblyConfig,
    )::NamedTuple
        records = multi_input_collect_fastq(reads)
        call_index = length(calls) + 1
        corrected_records = FASTX.FASTQ.Record[
            record for record in transform(call_index, records)
        ]
        ephemeral = config.output_dir === nothing
        corrected_fastq = if ephemeral
            tempname() * ".fastq"
        else
            joinpath(config.output_dir, "corrected.fastq")
        end
        multi_input_write_fastq(corrected_fastq, corrected_records)
        push!(calls, (
            technology = config.sequencing_tech,
            output_dir = config.output_dir,
            corrected_fastq = corrected_fastq,
            k = config.k,
            strategy = config.strategy,
        ))
        return (;
            corrected_fastq,
            corrected_reads = corrected_records,
            ephemeral,
        )
    end
    return correction_runner
end

function multi_input_write_fasta(
        path::AbstractString;
        records::Vector{Pair{String, String}} = ["contig_1" => "ACGTACGT"],
        gzip::Bool = false,
)::String
    mkpath(dirname(path))
    if gzip
        open(path, "w") do io
            gzip_stream = CodecZlib.GzipCompressorStream(io)
            try
                for (identifier, sequence) in records
                    write(gzip_stream, ">$(identifier)\n$(sequence)\n")
                end
            finally
                close(gzip_stream)
            end
        end
    else
        open(path, "w") do io
            for (identifier, sequence) in records
                write(io, ">$(identifier)\n$(sequence)\n")
            end
        end
    end
    return String(path)
end

function multi_input_fake_assembler_result(
        outdir::AbstractString;
        gzip::Bool = false,
        include_graph::Bool = false,
        use_contigs_key::Bool = false,
        include_polishing_artifacts::Bool = false,
)::NamedTuple
    extension = gzip ? ".fasta.gz" : ".fasta"
    assembly = multi_input_write_fasta(
        joinpath(outdir, "assembly$(extension)");
        gzip,
    )
    if include_graph
        graph = joinpath(outdir, "assembly.gfa")
        write(graph, "H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\n")
        if include_polishing_artifacts
            autocycler_assembly = multi_input_write_fasta(
                joinpath(outdir, "autocycler_assembly.fasta"),
            )
            polypolish_assembly = multi_input_write_fasta(
                joinpath(outdir, "polypolish_assembly.fasta"),
            )
            pypolca_report = joinpath(outdir, "pypolca.report")
            write(pypolca_report, "corrected_bases\t0\n")
            return (;
                assembly,
                graph,
                autocycler_assembly,
                polypolish_assembly,
                pypolca_report,
            )
        end
        return use_contigs_key ? (; contigs = assembly, graph) : (; assembly, graph)
    end
    return use_contigs_key ? (; contigs = assembly) : (; assembly)
end

const MULTI_INPUT_R1 = [
    multi_input_fastq_record("pair_a/1"),
    multi_input_fastq_record("pair_b/1", "TGCATGCA"),
]
const MULTI_INPUT_R2 = [
    multi_input_fastq_record("pair_a/2"),
    multi_input_fastq_record("pair_b/2", "TGCATGCA"),
]
const MULTI_INPUT_LONG = [
    multi_input_fastq_record("long_a", "ACGTACGTACGT"),
]

Test.@testset "multi-input hybrid assembly contracts (td-06er)" begin
    Test.@testset "typed config validation and distinct contracts" begin
        unicycler = Mycelia.Rhizomorph.UnicyclerHybridConfig(
            short_read_tech = :ultima,
            long_read_tech = :pacbio,
            correction_options = (; k = 17, strategy = :scalable),
            assembler_options = (; kmers = "21,33"),
            threads = 3,
        )
        Test.@test unicycler isa Mycelia.Rhizomorph.AbstractPairedShortLongAssemblyConfig
        Test.@test unicycler.short_read_tech == :ultima
        Test.@test unicycler.long_read_tech == :pacbio
        Test.@test Mycelia.Rhizomorph._long_read_correction_technology(unicycler) ==
                   :pacbio_clr
        Test.@test unicycler.threads == 3

        autocycler = Mycelia.Rhizomorph.AutocyclerPolishConfig(
            long_read_tech = :nanopore,
            autocycler_read_type = :ont_r9,
            threads = 4,
            jobs = 2,
            polypolish_careful = false,
            keep_intermediates = true,
            output_dir = tempname(),
        )
        Test.@test autocycler.autocycler_read_type == :ont_r9
        Test.@test !autocycler.polypolish_careful
        Test.@test autocycler.keep_intermediates
        Test.@test Mycelia.Rhizomorph.AutocyclerPolishConfig(
            long_read_tech = :pacbio,
        ).autocycler_read_type == :pacbio_clr
        Test.@test Mycelia.Rhizomorph.AutocyclerPolishConfig(
            long_read_tech = :pacbio_clr,
        ).autocycler_read_type == :pacbio_clr
        Test.@test Mycelia.Rhizomorph.AutocyclerPolishConfig(
            long_read_tech = :pacbio_hifi,
        ).autocycler_read_type == :pacbio_hifi

        test_throws_message(ArgumentError, "short_read_tech") do
            Mycelia.Rhizomorph.UnicyclerHybridConfig(short_read_tech = :nanopore)
        end
        test_throws_message(ArgumentError, "long_read_tech") do
            Mycelia.Rhizomorph.UnicyclerHybridConfig(long_read_tech = :illumina)
        end
        test_throws_message(ArgumentError, "incompatible") do
            Mycelia.Rhizomorph.AutocyclerPolishConfig(
                long_read_tech = :pacbio,
                autocycler_read_type = :ont_r10,
            )
        end
        test_throws_message(ArgumentError, "autocycler_read_type") do
            Mycelia.Rhizomorph.AutocyclerPolishConfig(
                autocycler_read_type = :unknown,
            )
        end
        test_throws_message(ArgumentError, "threads must be positive") do
            Mycelia.Rhizomorph.UnicyclerHybridConfig(threads = 0)
        end
        test_throws_message(ArgumentError, "jobs must be positive") do
            Mycelia.Rhizomorph.AutocyclerPolishConfig(jobs = 0)
        end
        test_throws_message(ArgumentError, "incompatible") do
            Mycelia.Rhizomorph.AutocyclerPolishConfig(
                long_read_tech = :pacbio_hifi,
                autocycler_read_type = :pacbio_clr,
            )
        end
        test_throws_message(ArgumentError, "incompatible") do
            Mycelia.Rhizomorph.AutocyclerPolishConfig(
                long_read_tech = :pacbio,
                autocycler_read_type = :pacbio_hifi,
            )
        end
        test_throws_message(ArgumentError, "requires a persistent output_dir") do
            Mycelia.Rhizomorph.AutocyclerPolishConfig(keep_intermediates = true)
        end
        test_throws_message(ArgumentError, "managed by the multi-input route") do
            Mycelia.Rhizomorph.UnicyclerHybridConfig(
                correction_options = (; sequencing_tech = :illumina),
            )
        end
        test_throws_message(ArgumentError, "Unicycler option :short_1") do
            Mycelia.Rhizomorph.UnicyclerHybridConfig(
                assembler_options = (; short_1 = "bypass.fastq"),
            )
        end
        test_throws_message(ArgumentError, "Unicycler option :executor") do
            Mycelia.Rhizomorph.UnicyclerHybridConfig(
                assembler_options = (; executor = :async),
            )
        end
        test_throws_message(ArgumentError, "output_dir must be a non-empty path") do
            Mycelia.Rhizomorph.UnicyclerHybridConfig(output_dir = "")
        end
        test_throws_message(ArgumentError, "output_dir must be a non-empty path") do
            Mycelia.Rhizomorph.UnicyclerHybridConfig(output_dir = "   ")
        end
    end

    Test.@testset "workflow root validates before FASTQ scans" begin
        mktempdir() do temp_dir
            stale_root = joinpath(temp_dir, "stale-root")
            mkpath(stale_root)
            stale_marker = joinpath(stale_root, "owned.txt")
            write(stale_marker, "keep\n")
            blocked_parent = joinpath(temp_dir, "blocked-parent")
            write(blocked_parent, "not a directory\n")
            blocked_root = joinpath(blocked_parent, "hybrid")
            file_root = joinpath(temp_dir, "file-root")
            write(file_root, "not a directory\n")

            correction_calls = Ref(0)
            assembler_calls = Ref(0)
            function root_validation_correction_runner(
                    reads::Any,
                    config::Mycelia.Rhizomorph.AssemblyConfig,
            )::NamedTuple
                correction_calls[] += 1
                error("correction must not run")
            end
            function root_validation_assembler_runner(
                    inputs::Any,
                    outdir::AbstractString,
            )::NamedTuple
                assembler_calls[] += 1
                error("assembler must not run")
            end
            missing_inputs = (
                (joinpath(temp_dir, "missing-r1.fastq"),),
                (joinpath(temp_dir, "missing-r2.fastq"),),
                (joinpath(temp_dir, "missing-long.fastq"),),
            )

            for (message, output_dir) in (
                    ("prevent stale hybrid assembly reuse", stale_root),
                    ("non-directory ancestor", blocked_root),
                    ("exists but is not a directory", file_root),
            )
                config = Mycelia.Rhizomorph.UnicyclerHybridConfig(; output_dir)
                test_throws_message(ArgumentError, message) do
                    Mycelia.Rhizomorph._assemble_paired_short_long(
                        (missing_inputs[1], missing_inputs[2]),
                        missing_inputs[3],
                        config,
                        :unicycler;
                        correction_runner = root_validation_correction_runner,
                        assembler_runner = root_validation_assembler_runner,
                    )
                end
            end
            Test.@test correction_calls[] == 0
            Test.@test assembler_calls[] == 0
            Test.@test isfile(stale_marker)
            Test.@test read(stale_marker, String) == "keep\n"
            Test.@test !ispath(blocked_root)
            Test.@test isfile(file_root)

            for output_dir in (
                    joinpath(temp_dir, "absent-invalid-root"),
                    mktempdir(temp_dir),
            )
                config = Mycelia.Rhizomorph.UnicyclerHybridConfig(; output_dir)
                test_throws_message(ArgumentError, "different counts") do
                    Mycelia.Rhizomorph._assemble_paired_short_long(
                        (MULTI_INPUT_R1, MULTI_INPUT_R2[1:1]),
                        MULTI_INPUT_LONG,
                        config,
                        :unicycler;
                        correction_runner = root_validation_correction_runner,
                        assembler_runner = root_validation_assembler_runner,
                    )
                end
                if isdir(output_dir)
                    Test.@test isempty(readdir(output_dir))
                else
                    Test.@test !ispath(output_dir)
                end
            end
            Test.@test correction_calls[] == 0
            Test.@test assembler_calls[] == 0
        end
    end

    Test.@testset "raw and corrected mate validation" begin
        Test.@test Mycelia.Rhizomorph._validate_paired_reads(
            MULTI_INPUT_R1,
            MULTI_INPUT_R2,
            "input",
        ) == 2

        test_throws_message(ArgumentError, "different counts") do
            Mycelia.Rhizomorph._validate_paired_reads(
                MULTI_INPUT_R1,
                MULTI_INPUT_R2[1:1],
                "input",
            )
        end
        test_throws_message(ArgumentError, "out of sync at record 1") do
            Mycelia.Rhizomorph._validate_paired_reads(
                MULTI_INPUT_R1,
                reverse(MULTI_INPUT_R2),
                "input",
            )
        end
        reversed_r1 = [multi_input_fastq_record("pair_1/2")]
        reversed_r2 = [multi_input_fastq_record("pair_1/1")]
        test_throws_message(ArgumentError, "invalid explicit mate roles") do
            Mycelia.Rhizomorph._validate_paired_reads(
                reversed_r1,
                reversed_r2,
                "input",
            )
        end
        test_throws_message(ArgumentError, "must be non-empty") do
            Mycelia.Rhizomorph._validate_paired_reads(
                FASTX.FASTQ.Record[],
                FASTX.FASTQ.Record[],
                "input",
            )
        end
        test_throws_message(ArgumentError, "sources must be distinct") do
            Mycelia.Rhizomorph._validate_paired_reads(
                MULTI_INPUT_R1,
                MULTI_INPUT_R1,
                "input",
            )
        end

        raw_dir = mktempdir()
        r1_path = multi_input_write_fastq(joinpath(raw_dir, "r1.fastq"), MULTI_INPUT_R1)
        r2_path = multi_input_write_fastq(joinpath(raw_dir, "r2.fastq"), MULTI_INPUT_R2)
        Test.@test Mycelia.Rhizomorph._validate_paired_reads(
            [r1_path],
            [r2_path],
            "input",
        ) == 2
        test_throws_message(ArgumentError, "sources must be distinct") do
            Mycelia.Rhizomorph._validate_paired_reads(
                [r1_path],
                [r1_path],
                "input",
            )
        end

        split_r1_paths = [
            multi_input_write_fastq(
                joinpath(raw_dir, "split_r1_$(index).fastq"),
                MULTI_INPUT_R1[index:index],
            ) for index in eachindex(MULTI_INPUT_R1)
        ]
        split_r2_paths = [
            multi_input_write_fastq(
                joinpath(raw_dir, "split_r2_$(index).fastq"),
                MULTI_INPUT_R2[index:index],
            ) for index in eachindex(MULTI_INPUT_R2)
        ]
        Test.@test Mycelia.Rhizomorph._validate_paired_reads(
            split_r1_paths,
            split_r2_paths,
            "input",
        ) == 2
        test_throws_message(ArgumentError, "sources must be distinct") do
            Mycelia.Rhizomorph._validate_paired_reads(
                split_r1_paths,
                [split_r2_paths[1], split_r1_paths[2]],
                "input",
            )
        end
        test_throws_message(ArgumentError, "R1=2, R2=1") do
            Mycelia.Rhizomorph._validate_paired_reads(
                split_r1_paths,
                split_r2_paths[1:1],
                "input",
            )
        end

        correction_calls = NamedTuple[]
        correction_runner = multi_input_fake_correction_runner(
            correction_calls;
            transform = function (index::Int, records::Vector{FASTX.FASTQ.Record})
                if index == 2
                    replacement = copy(records)
                    replacement[1] = multi_input_fastq_record("wrong_pair/2")
                    return replacement
                end
                return records
            end,
        )
        config = Mycelia.Rhizomorph.UnicyclerHybridConfig()
        assembler_called = Ref(false)
        assembler_runner = function (inputs::Any, outdir::AbstractString)
            assembler_called[] = true
            return multi_input_fake_assembler_result(outdir)
        end
        test_throws_message(ArgumentError, "corrected paired short reads are out of sync") do
            Mycelia.Rhizomorph._assemble_paired_short_long(
                (MULTI_INPUT_R1, MULTI_INPUT_R2),
                MULTI_INPUT_LONG,
                config,
                :unicycler;
                correction_runner,
                assembler_runner,
            )
        end
        Test.@test !assembler_called[]
        Test.@test length(correction_calls) == 3
        Test.@test all(call -> !isfile(call.corrected_fastq), correction_calls)

        for (message, transform) in (
                (
                    "short_r1 correction changed read count",
                    (records::Vector{FASTX.FASTQ.Record}) -> records[1:1],
                ),
                (
                    "short_r1 correction changed read order or identifier",
                    (records::Vector{FASTX.FASTQ.Record}) -> reverse(records),
                ),
        )
            identical_change_calls = NamedTuple[]
            identical_change_runner = multi_input_fake_correction_runner(
                identical_change_calls;
                transform = (index, records) -> index in (1, 2) ?
                                                transform(records) : records,
            )
            identical_change_assembler_called = Ref(false)
            test_throws_message(ArgumentError, message) do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (MULTI_INPUT_R1, MULTI_INPUT_R2),
                    MULTI_INPUT_LONG,
                    config,
                    :unicycler;
                    correction_runner = identical_change_runner,
                    assembler_runner = (inputs, outdir) -> begin
                        identical_change_assembler_called[] = true
                        multi_input_fake_assembler_result(outdir)
                    end,
                )
            end
            Test.@test !identical_change_assembler_called[]
            Test.@test length(identical_change_calls) == 3
            Test.@test all(
                call -> !isfile(call.corrected_fastq),
                identical_change_calls,
            )
        end

        two_long_reads = [
            multi_input_fastq_record("long_a", "ACGTACGTACGT"),
            multi_input_fastq_record("long_b", "TGCATGCATGCA"),
        ]
        reordered_long_calls = NamedTuple[]
        reordered_long_runner = multi_input_fake_correction_runner(
            reordered_long_calls;
            transform = (index, records) -> index == 3 ? reverse(records) : records,
        )
        test_throws_message(
            ArgumentError,
            "long_reads correction changed read order or identifier",
        ) do
            Mycelia.Rhizomorph._assemble_paired_short_long(
                (MULTI_INPUT_R1, MULTI_INPUT_R2),
                two_long_reads,
                config,
                :unicycler;
                correction_runner = reordered_long_runner,
                assembler_runner = (inputs, outdir) ->
                    multi_input_fake_assembler_result(outdir),
            )
        end
        Test.@test all(
            call -> !isfile(call.corrected_fastq),
            reordered_long_calls,
        )
    end

    Test.@testset "gzip inputs reach independent correction" begin
        mktempdir() do temp_dir
            r1 = multi_input_write_fastq(
                joinpath(temp_dir, "R1.fastq"),
                MULTI_INPUT_R1,
            )
            r2 = multi_input_write_fastq(
                joinpath(temp_dir, "R2.fastq"),
                MULTI_INPUT_R2,
            )
            long_reads = multi_input_write_fastq(
                joinpath(temp_dir, "long.fastq"),
                MULTI_INPUT_LONG,
            )
            r1_gzip = multi_input_gzip_file(r1, "$(r1).gz")
            r2_gzip = multi_input_gzip_file(r2, "$(r2).gz")
            long_gzip = multi_input_gzip_file(long_reads, "$(long_reads).gz")
            correction_calls = NamedTuple[]
            result = Mycelia.Rhizomorph._assemble_paired_short_long(
                (r1_gzip, r2_gzip),
                long_gzip,
                Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                :unicycler;
                correction_runner = multi_input_fake_correction_runner(
                    correction_calls,
                ),
                assembler_runner = (inputs, outdir) ->
                    multi_input_fake_assembler_result(outdir),
            )
            Test.@test [call.technology for call in correction_calls] ==
                       [:illumina, :illumina, :nanopore]
            Test.@test result.assembly_stats["input_read_counts"] == Dict(
                "short_r1" => 2,
                "short_r2" => 2,
                "long_reads" => 1,
            )
            Test.@test all(call -> !isfile(call.corrected_fastq), correction_calls)
        end
    end

    Test.@testset "lightweight Stage-1 cleanup ownership" begin
        Test.@test fieldnames(Mycelia.Rhizomorph._Stage1CleanupToken) ==
                   (:corrected_fastq, :ephemeral)
        cleanup_tokens = Mycelia.Rhizomorph._Stage1CleanupToken[]
        correction_calls = NamedTuple[]
        correction_runner = multi_input_fake_correction_runner(correction_calls)
        corrected = Mycelia.Rhizomorph._correct_read_set!(
            cleanup_tokens,
            MULTI_INPUT_LONG,
            :nanopore,
            :long_reads,
            (; k = 13, strategy = :scalable),
            mktempdir(),
            false,
            correction_runner,
        )
        Test.@test corrected.count == 1
        Test.@test length(cleanup_tokens) == 1
        Test.@test cleanup_tokens[1].corrected_fastq == corrected.path
        Mycelia.Rhizomorph._cleanup_multi_input_stages!(cleanup_tokens)
        Test.@test !isfile(corrected.path)

        counted_fastq = tempname() * ".fastq"
        multi_input_write_fastq(counted_fastq, MULTI_INPUT_LONG)
        function counted_runner(
                reads::Any,
                config::Mycelia.Rhizomorph.AssemblyConfig,
        )::NamedTuple
            return (;
                corrected_fastq = counted_fastq,
                corrected_read_count = 1,
                ephemeral = true,
                correction_graph = fill(1, 10),
            )
        end
        empty!(cleanup_tokens)
        corrected_from_count = Mycelia.Rhizomorph._correct_read_set!(
            cleanup_tokens,
            MULTI_INPUT_LONG,
            :nanopore,
            :long_reads,
            (; k = 13, strategy = :scalable),
            mktempdir(),
            false,
            counted_runner,
        )
        Test.@test corrected_from_count.count == 1
        Test.@test cleanup_tokens == [
            Mycelia.Rhizomorph._Stage1CleanupToken(counted_fastq, true),
        ]
        Mycelia.Rhizomorph._cleanup_multi_input_stages!(cleanup_tokens)
        Test.@test !isfile(counted_fastq)

        valid_fastq = multi_input_write_fastq(
            tempname() * ".fastq",
            MULTI_INPUT_LONG,
        )
        empty_fastq = tempname() * ".fastq"
        touch(empty_fastq)
        fasta_output = tempname() * ".fasta"
        write(fasta_output, ">not_fastq\nACGT\n")
        malformed_stages = [
            ("missing corrected_fastq", (; ephemeral = true)),
            (
                "no non-empty corrected FASTQ",
                (;
                    corrected_fastq = tempname() * ".fastq",
                    corrected_read_count = 1,
                    ephemeral = true,
                ),
            ),
            (
                "no non-empty corrected FASTQ",
                (;
                    corrected_fastq = empty_fastq,
                    corrected_read_count = 1,
                    ephemeral = true,
                ),
            ),
            (
                "missing corrected_read_count",
                (; corrected_fastq = valid_fastq, ephemeral = false),
            ),
            (
                "produced 0 corrected reads",
                (;
                    corrected_fastq = valid_fastq,
                    corrected_read_count = 0,
                    ephemeral = false,
                ),
            ),
            (
                "reported 7 corrected reads",
                (;
                    corrected_fastq = valid_fastq,
                    corrected_read_count = 7,
                    ephemeral = false,
                ),
            ),
            (
                "correction output is not FASTQ",
                (;
                    corrected_fastq = fasta_output,
                    corrected_read_count = 1,
                    ephemeral = false,
                ),
            ),
        ]
        for (message, stage) in malformed_stages
            tokens = Mycelia.Rhizomorph._Stage1CleanupToken[]
            bad_runner = (reads, config) -> stage
            test_throws_message(ErrorException, message) do
                Mycelia.Rhizomorph._correct_read_set!(
                    tokens,
                    MULTI_INPUT_LONG,
                    :nanopore,
                    :long_reads,
                    (; k = 13, strategy = :scalable),
                    mktempdir(),
                    false,
                    bad_runner,
                )
            end
            Mycelia.Rhizomorph._cleanup_multi_input_stages!(tokens)
        end
        rm(valid_fastq; force = true)
        rm(empty_fastq; force = true)
        rm(fasta_output; force = true)
    end

    Test.@testset "independent correction profiles and ephemeral cleanup" begin
        correction_calls = NamedTuple[]
        correction_runner = multi_input_fake_correction_runner(correction_calls)
        captured_inputs = Ref{Any}()
        captured_outdir = Ref("")
        assembler_runner = function (
                inputs::Mycelia.Rhizomorph._CorrectedPairedShortLong,
                outdir::AbstractString,
        )
            captured_inputs[] = inputs
            captured_outdir[] = String(outdir)
            return multi_input_fake_assembler_result(outdir; gzip = true)
        end
        config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
            short_read_tech = :illumina,
            long_read_tech = :nanopore,
            correction_options = (; k = 17, strategy = :scalable),
        )
        result = Mycelia.Rhizomorph._assemble_paired_short_long(
            (MULTI_INPUT_R1, MULTI_INPUT_R2),
            MULTI_INPUT_LONG,
            config,
            :unicycler;
            correction_runner,
            assembler_runner,
        )

        Test.@test [call.technology for call in correction_calls] ==
                   [:illumina, :illumina, :nanopore]
        Test.@test all(call -> call.k == 17, correction_calls)
        Test.@test all(call -> call.strategy == :scalable, correction_calls)
        Test.@test captured_inputs[].short_r1.path == correction_calls[1].corrected_fastq
        Test.@test captured_inputs[].short_r2.path == correction_calls[2].corrected_fastq
        Test.@test captured_inputs[].long_reads.path == correction_calls[3].corrected_fastq
        Test.@test result.contig_names == ["contig_1"]
        Test.@test result.contigs == ["ACGTACGT"]
        Test.@test result.assembly_stats["workflow"] == "unicycler"
        Test.@test result.assembly_stats["method"] == "HybridAssembly"
        Test.@test result.assembly_stats["assembler"] == "unicycler"
        Test.@test result.assembly_stats["input_contract"] == "paired_short_long"
        Test.@test result.assembly_stats["input_read_counts"] == Dict(
            "short_r1" => 2,
            "short_r2" => 2,
            "long_reads" => 1,
        )
        Test.@test result.assembly_stats["corrected_fastqs"] === nothing
        Test.@test result.assembly_stats["tool_artifacts"] === nothing
        Test.@test result.assembly_stats["input_technologies"] == Dict(
            "short_r1" => "illumina",
            "short_r2" => "illumina",
            "long_reads" => "nanopore",
        )
        Test.@test result.assembly_stats["workflow_settings"] == Dict(
            "workflow" => "unicycler",
            "assembler" => "unicycler",
            "threads" => config.threads,
            "correction_options" => Dict(
                "k" => 17,
                "strategy" => "scalable",
            ),
            "assembler_options" => Dict{String, Any}(),
        )
        Test.@test result.assembly_stats["correction_options"] == Dict(
            "k" => 17,
            "strategy" => "scalable",
        )
        Test.@test !result.gfa_compatible
        Test.@test all(call -> !isfile(call.corrected_fastq), correction_calls)
        Test.@test !isdir(dirname(captured_outdir[]))
    end

    Test.@testset "sibling adapters map exact corrected paths" begin
        inputs = Mycelia.Rhizomorph._CorrectedPairedShortLong(
            Mycelia.Rhizomorph._CorrectedReadSet("/corrected/r1.fastq", 2, :illumina),
            Mycelia.Rhizomorph._CorrectedReadSet("/corrected/r2.fastq", 2, :illumina),
            Mycelia.Rhizomorph._CorrectedReadSet("/corrected/long.fastq", 1, :nanopore),
        )

        unicycler_kwargs = Ref{Any}()
        unicycler_runner = function (; kwargs...)
            unicycler_kwargs[] = (; kwargs...)
            return (; assembly = "assembly.fasta", graph = "assembly.gfa")
        end
        unicycler_config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
            assembler_options = (; kmers = "21,33"),
            threads = 6,
        )
        Mycelia.Rhizomorph._run_multi_input_assembler(
            Val(:unicycler),
            inputs,
            "/workflow/unicycler",
            unicycler_config;
            runner = unicycler_runner,
        )
        Test.@test unicycler_kwargs[] == (
            kmers = "21,33",
            short_1 = "/corrected/r1.fastq",
            short_2 = "/corrected/r2.fastq",
            long_reads = "/corrected/long.fastq",
            outdir = "/workflow/unicycler",
            threads = 6,
        )

        autocycler_kwargs = Ref{Any}()
        autocycler_runner = function (; kwargs...)
            autocycler_kwargs[] = (; kwargs...)
            return (; assembly = "polished.fasta", graph = "consensus.gfa")
        end
        autocycler_config = Mycelia.Rhizomorph.AutocyclerPolishConfig(
            autocycler_read_type = :ont_r9,
            threads = 8,
            jobs = 2,
            polypolish_careful = false,
        )
        Mycelia.Rhizomorph._run_multi_input_assembler(
            Val(:autocycler_polished),
            inputs,
            "/workflow/autocycler",
            autocycler_config;
            runner = autocycler_runner,
        )
        Test.@test autocycler_kwargs[] == (
            long_reads = "/corrected/long.fastq",
            short_reads_1 = "/corrected/r1.fastq",
            short_reads_2 = "/corrected/r2.fastq",
            out_dir = "/workflow/autocycler",
            threads = 8,
            jobs = 2,
            read_type = "ont_r9",
            polypolish_careful = false,
            keep_intermediates = false,
        )

    end

    Test.@testset "Autocycler-polished workflow provenance" begin
        correction_calls = NamedTuple[]
        correction_runner = multi_input_fake_correction_runner(correction_calls)
        output_dir = mktempdir()
        requested_threads = Mycelia.AUTOCYCLER_MAX_ASSEMBLY_THREADS + 7
        assembler_runner = function (inputs::Any, outdir::AbstractString)
            Test.@test inputs.long_reads.technology == :pacbio_hifi
            result = multi_input_fake_assembler_result(
                outdir;
                include_graph = true,
                include_polishing_artifacts = true,
            )
            return merge(
                result,
                (;
                    toolchain = Dict(
                        "autocycler_script_revision" => "immutable-test-revision",
                        "packages" => Dict("autocycler" => "0.5.2"),
                    ),
                ),
            )
        end
        config = Mycelia.Rhizomorph.AutocyclerPolishConfig(;
            long_read_tech = :pacbio_hifi,
            autocycler_read_type = :pacbio_hifi,
            output_dir,
            threads = requested_threads,
            jobs = 2,
        )
        result = Mycelia.Rhizomorph._assemble_paired_short_long(
            (MULTI_INPUT_R1, MULTI_INPUT_R2),
            MULTI_INPUT_LONG,
            config,
            :autocycler_polished;
            correction_runner,
            assembler_runner,
        )
        Test.@test [call.technology for call in correction_calls] ==
                   [:illumina, :illumina, :pacbio_hifi]
        Test.@test result.assembly_stats["workflow"] == "autocycler_polished"
        Test.@test result.assembly_stats["method"] == "HybridAssembly"
        Test.@test result.assembly_stats["assembler"] == "autocycler"
        Test.@test result.assembly_stats["input_contract"] == "paired_short_long"
        Test.@test result.assembly_stats["polishers"] ==
                   ["polypolish-careful", "pypolca-careful"]
        Test.@test result.assembly_stats["input_technologies"] == Dict(
            "short_r1" => "illumina",
            "short_r2" => "illumina",
            "long_reads" => "pacbio_hifi",
        )
        Test.@test result.assembly_stats["workflow_settings"] == Dict(
            "workflow" => "autocycler_polished",
            "assembler" => "autocycler",
            "threads" => requested_threads,
            "autocycler_assembly_threads" =>
                Mycelia.AUTOCYCLER_MAX_ASSEMBLY_THREADS,
            "polishing_threads" => requested_threads,
            "jobs" => 2,
            "autocycler_read_type" => "pacbio_hifi",
            "long_read_correction_tech" => "pacbio_hifi",
            "polypolish_careful" => true,
            "pypolca_careful" => true,
            "keep_intermediates" => false,
            "correction_options" => Dict(
                "k" => 13,
                "strategy" => "scalable",
            ),
        )
        tool_artifacts = result.assembly_stats["tool_artifacts"]
        Test.@test Set(keys(tool_artifacts)) == Set([
            "final_assembly",
            "raw_graph",
            "autocycler_assembly",
            "polypolish_assembly",
            "pypolca_report",
        ])
        Test.@test all(isfile, values(tool_artifacts))
        Test.@test result.assembly_stats["raw_graph"] == tool_artifacts["raw_graph"]
        Test.@test result.assembly_stats["toolchain"] == Dict(
            "autocycler_script_revision" => "immutable-test-revision",
            "packages" => Dict{String, Any}("autocycler" => "0.5.2"),
        )
        Test.@test isempty(result.assembly_stats["retained_intermediates"])

        noncareful_calls = NamedTuple[]
        noncareful_result = Mycelia.Rhizomorph._assemble_paired_short_long(
            (MULTI_INPUT_R1, MULTI_INPUT_R2),
            MULTI_INPUT_LONG,
            Mycelia.Rhizomorph.AutocyclerPolishConfig(polypolish_careful = false),
            :autocycler_polished;
            correction_runner = multi_input_fake_correction_runner(noncareful_calls),
            assembler_runner = (inputs, outdir) ->
                multi_input_fake_assembler_result(outdir),
        )
        Test.@test noncareful_result.assembly_stats["polishers"] ==
                   ["polypolish", "pypolca-careful"]
        Test.@test noncareful_result.assembly_stats["tool_artifacts"] === nothing
    end

    Test.@testset "cleanup on assembler failure" begin
        correction_calls = NamedTuple[]
        correction_runner = multi_input_fake_correction_runner(correction_calls)
        workflow_root = Ref("")
        assembler_runner = function (inputs::Any, outdir::AbstractString)
            workflow_root[] = dirname(String(outdir))
            error("synthetic assembler failure")
        end
        test_throws_message(ErrorException, "synthetic assembler failure") do
            Mycelia.Rhizomorph._assemble_paired_short_long(
                (MULTI_INPUT_R1, MULTI_INPUT_R2),
                MULTI_INPUT_LONG,
                Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                :unicycler;
                correction_runner,
                assembler_runner,
            )
        end
        Test.@test length(correction_calls) == 3
        Test.@test all(call -> !isfile(call.corrected_fastq), correction_calls)
        Test.@test !isdir(workflow_root[])

        interrupted_correction_calls = NamedTuple[]
        successful_correction = multi_input_fake_correction_runner(
            interrupted_correction_calls,
        )
        correction_runner = function (
                reads::Any,
                config::Mycelia.Rhizomorph.AssemblyConfig,
        )
            if length(interrupted_correction_calls) == 1
                error("synthetic second-correction failure")
            end
            return successful_correction(reads, config)
        end
        test_throws_message(ErrorException, "second-correction failure") do
            Mycelia.Rhizomorph._assemble_paired_short_long(
                (MULTI_INPUT_R1, MULTI_INPUT_R2),
                MULTI_INPUT_LONG,
                Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                :unicycler;
                correction_runner,
                assembler_runner,
            )
        end
        Test.@test length(interrupted_correction_calls) == 1
        Test.@test !isfile(interrupted_correction_calls[1].corrected_fastq)
    end

    Test.@testset "caller-owned artifacts persist" begin
        output_dir = mktempdir()
        correction_calls = NamedTuple[]
        correction_runner = multi_input_fake_correction_runner(correction_calls)
        assembly_path = Ref("")
        graph_path = Ref("")
        assembler_runner = function (inputs::Any, outdir::AbstractString)
            result = multi_input_fake_assembler_result(
                outdir;
                include_graph = true,
            )
            assembly_path[] = result.assembly
            graph_path[] = result.graph
            return result
        end
        config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
            output_dir = output_dir,
        )
        result = Mycelia.Rhizomorph._assemble_paired_short_long(
            (MULTI_INPUT_R1, MULTI_INPUT_R2),
            MULTI_INPUT_LONG,
            config,
            :unicycler;
            correction_runner,
            assembler_runner,
        )

        expected_corrected = Dict(
            "short_r1" => joinpath(output_dir, "corrected", "short_r1", "corrected.fastq"),
            "short_r2" => joinpath(output_dir, "corrected", "short_r2", "corrected.fastq"),
            "long_reads" => joinpath(output_dir, "corrected", "long_reads", "corrected.fastq"),
        )
        Test.@test result.assembly_stats["corrected_fastqs"] == expected_corrected
        Test.@test all(isfile, values(expected_corrected))
        Test.@test isfile(assembly_path[])
        Test.@test isfile(graph_path[])
        Test.@test result.assembly_stats["raw_graph"] == graph_path[]
        Test.@test result.assembly_stats["output_dir"] == output_dir
        Test.@test result.assembly_stats["tool_artifacts"] == Dict(
            "final_assembly" => assembly_path[],
            "raw_graph" => graph_path[],
        )

        stale_dir = mktempdir()
        write(joinpath(stale_dir, "stale.txt"), "stale")
        stale_config = Mycelia.Rhizomorph.UnicyclerHybridConfig(output_dir = stale_dir)
        test_throws_message(ArgumentError, "prevent stale hybrid assembly reuse") do
            Mycelia.Rhizomorph._assemble_paired_short_long(
                (MULTI_INPUT_R1, MULTI_INPUT_R2),
                MULTI_INPUT_LONG,
                stale_config,
                :unicycler;
                correction_runner,
                assembler_runner,
            )
        end
        Test.@test isfile(joinpath(stale_dir, "stale.txt"))
    end

    Test.@testset "gzip FASTA wrapping and fail-loud artifacts" begin
        artifact_dir = mktempdir()
        gzip_fasta = multi_input_write_fasta(
            joinpath(artifact_dir, "contigs.fasta.gz");
            records = ["alpha" => "ACGT", "beta" => "TTAA"],
            gzip = true,
        )
        records = Mycelia.Rhizomorph._read_external_fasta(gzip_fasta)
        Test.@test [String(FASTX.identifier(record)) for record in records] ==
                   ["alpha", "beta"]

        wrapped = Mycelia.Rhizomorph._wrap_multi_input_assembly(
            (; contigs = gzip_fasta),
            :unicycler;
            input_counts = Dict(
                "short_r1" => 1,
                "short_r2" => 1,
                "long_reads" => 1,
            ),
            corrected_counts = Dict(
                "short_r1" => 1,
                "short_r2" => 1,
                "long_reads" => 1,
            ),
            corrected_paths = nothing,
            short_read_tech = :illumina,
            long_read_tech = :nanopore,
            correction_options = (; k = 13),
            output_dir = nothing,
            polishers = String[],
            workflow_settings = Dict{String, Any}(
                "workflow" => "unicycler",
                "assembler" => "unicycler",
            ),
            input_technologies = Dict(
                "short_r1" => "illumina",
                "short_r2" => "illumina",
                "long_reads" => "nanopore",
            ),
        )
        Test.@test wrapped.contig_names == ["alpha", "beta"]
        Test.@test wrapped.contigs == ["ACGT", "TTAA"]

        missing = joinpath(artifact_dir, "missing.fasta")
        test_throws_message(ErrorException, "no non-empty contigs FASTA") do
            Mycelia.Rhizomorph._read_external_fasta(missing)
        end
        zero_byte = joinpath(artifact_dir, "zero.fasta")
        touch(zero_byte)
        test_throws_message(ErrorException, "no non-empty contigs FASTA") do
            Mycelia.Rhizomorph._read_external_fasta(zero_byte)
        end

        test_throws_message(ErrorException, "neither :assembly nor :contigs") do
            Mycelia.Rhizomorph._primary_assembly_path((; outdir = artifact_dir))
        end
        test_throws_message(ErrorException, "produced no graph") do
            Mycelia.Rhizomorph._wrap_multi_input_assembly(
                (; assembly = gzip_fasta, graph = joinpath(artifact_dir, "missing.gfa")),
                :autocycler_polished;
                input_counts = Dict{String, Int}(),
                corrected_counts = Dict{String, Int}(),
                corrected_paths = nothing,
                short_read_tech = :illumina,
                long_read_tech = :nanopore,
                correction_options = (;),
                output_dir = nothing,
                polishers = ["polypolish", "pypolca-careful"],
                workflow_settings = Dict{String, Any}(
                    "workflow" => "autocycler_polished",
                    "assembler" => "autocycler",
                ),
                input_technologies = Dict(
                    "short_r1" => "illumina",
                    "short_r2" => "illumina",
                    "long_reads" => "nanopore",
                ),
            )
        end
    end

    Test.@testset "public dispatch fails before external tools" begin
        no_role_pairs = [multi_input_fastq_record("same_pair")]
        test_throws_message(ArgumentError, "sources must be distinct") do
            Mycelia.Rhizomorph.assemble_hybrid(
                (no_role_pairs, no_role_pairs),
                MULTI_INPUT_LONG;
                config = Mycelia.Rhizomorph.UnicyclerHybridConfig(),
            )
        end
        test_throws_message(ArgumentError, "Long-read input must be distinct") do
            Mycelia.Rhizomorph.assemble_hybrid(
                (MULTI_INPUT_R1, MULTI_INPUT_R2),
                MULTI_INPUT_R1;
                config = Mycelia.Rhizomorph.UnicyclerHybridConfig(),
            )
        end
        test_throws_message(ArgumentError, "different counts") do
            Mycelia.Rhizomorph.assemble_hybrid(
                (MULTI_INPUT_R1, MULTI_INPUT_R2[1:1]),
                MULTI_INPUT_LONG;
                config = Mycelia.Rhizomorph.UnicyclerHybridConfig(),
            )
        end
        Test.@test !isdefined(Mycelia.Rhizomorph, :assemble_dual_long)
        mktempdir() do temp_dir
            paired_path = multi_input_write_fastq(
                joinpath(temp_dir, "same-physical.fastq"),
                no_role_pairs,
            )
            paired_alias = joinpath(temp_dir, "same-physical-alias.fastq")
            symlink(paired_path, paired_alias)
            test_throws_message(ArgumentError, "sources must be distinct") do
                Mycelia.Rhizomorph.assemble_hybrid(
                    ([paired_path], [paired_alias]),
                    MULTI_INPUT_LONG;
                    config = Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                )
            end
            distinct_r2 = multi_input_write_fastq(
                joinpath(temp_dir, "distinct-r2.fastq"),
                [multi_input_fastq_record("same_pair")],
            )
            test_throws_message(ArgumentError, "Long-read input must be distinct") do
                Mycelia.Rhizomorph.assemble_hybrid(
                    (paired_path, distinct_r2),
                    paired_alias;
                    config = Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                )
            end

            long_path = multi_input_write_fastq(
                joinpath(temp_dir, "long.fastq"),
                MULTI_INPUT_LONG,
            )
            paired_hardlink = joinpath(temp_dir, "same-physical-hardlink.fastq")
            Base.hardlink(paired_path, paired_hardlink)
            correction_calls = Ref(0)
            assembler_calls = Ref(0)
            function alias_guard_correction_runner(
                    reads::Any,
                    config::Mycelia.Rhizomorph.AssemblyConfig,
            )::NamedTuple
                correction_calls[] += 1
                error("correction must not run")
            end
            function alias_guard_assembler_runner(
                    inputs::Any,
                    outdir::AbstractString,
            )::NamedTuple
                assembler_calls[] += 1
                error("assembler must not run")
            end
            alias_cases = (
                (
                    message = "sources must be distinct",
                    short_r1 = (paired_path,),
                    short_r2 = (path for path in Any[paired_path]),
                    long_reads = (long_path,),
                    output_name = "exact-generator-alias",
                ),
                (
                    message = "sources must be distinct",
                    short_r1 = (paired_path,),
                    short_r2 = (path for path in (paired_alias,)),
                    long_reads = (long_path,),
                    output_name = "symlink-generator-alias",
                ),
                (
                    message = "Long-read input must be distinct",
                    short_r1 = (paired_path,),
                    short_r2 = (distinct_r2,),
                    long_reads = (path for path in (paired_hardlink,)),
                    output_name = "hardlink-generator-alias",
                ),
            )
            for alias_case in alias_cases
                output_dir = joinpath(temp_dir, alias_case.output_name)
                config = Mycelia.Rhizomorph.UnicyclerHybridConfig(; output_dir)
                test_throws_message(ArgumentError, alias_case.message) do
                    Mycelia.Rhizomorph._assemble_paired_short_long(
                        (alias_case.short_r1, alias_case.short_r2),
                        alias_case.long_reads,
                        config,
                        :unicycler;
                        correction_runner = alias_guard_correction_runner,
                        assembler_runner = alias_guard_assembler_runner,
                    )
                end
                Test.@test !ispath(output_dir)
            end
            Test.@test correction_calls[] == 0
            Test.@test assembler_calls[] == 0

            excluded_outdir = joinpath(temp_dir, "excluded-mixed-metamdbg")
            test_throws_message(ArgumentError, "exactly one input technology") do
                Mycelia.run_metamdbg(
                    hifi_reads = "hifi.fastq",
                    ont_reads = "ont.fastq",
                    outdir = excluded_outdir,
                )
            end
            Test.@test !ispath(excluded_outdir)

            missing_outdir = joinpath(temp_dir, "missing-input-metamdbg")
            test_throws_message(ArgumentError, "exactly one input technology") do
                Mycelia.run_metamdbg(outdir = missing_outdir)
            end
            Test.@test !ispath(missing_outdir)

            for (label, keyword_arguments) in (
                    ("hifi_reads", (; hifi_reads = String[])),
                    ("ont_reads", (; ont_reads = "")),
                    ("ont_reads", (; ont_reads = [" "])),
            )
                empty_outdir = joinpath(temp_dir, "empty-$(label)-$(gensym())")
                test_throws_message(ArgumentError, "at least one non-empty path") do
                    Mycelia.run_metamdbg(;
                        keyword_arguments...,
                        outdir = empty_outdir,
                    )
                end
                Test.@test !ispath(empty_outdir)
            end
        end
    end
end
