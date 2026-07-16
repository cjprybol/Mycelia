# Default-CI contract tests for multi-read-set hybrid assembly.
#
# These tests use injected correction and assembler runners, so they exercise
# pairing, routing, provenance, and cleanup without installing external tools.
# The gated companion smoke test owns real third-party execution.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/multi_input_hybrid_assembly_test.jl")'

import CodecZlib
import FASTX
import Logging
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

function multi_input_reported_correction_runner(
        calls::Vector{NamedTuple},
)::Function
    base_runner = multi_input_fake_correction_runner(calls)
    function correction_runner(
            reads::Any,
            config::Mycelia.Rhizomorph.AssemblyConfig,
    )::NamedTuple
        stage = base_runner(reads, config)
        technology = config.sequencing_tech
        indel_params = if technology in (:nanopore, :pacbio_clr)
            Mycelia.IndelDecodeParams(
                0.1,
                0.4,
                0.4,
                0.2,
                0.2,
                3,
                4,
                64,
            )
        else
            nothing
        end
        max_k = config.k === nothing ? 13 : max(config.k, 13)
        substitution_error_rate = technology == :pacbio_hifi ? 0.001 : nothing
        return merge(
            stage,
            (;
                knobs = (
                    strategy = config.strategy,
                    skip_solid = technology != :nanopore,
                ),
                max_k,
                indel_params,
                substitution_error_rate,
            ),
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

function multi_input_fake_unicycler_toolchain()::Dict{String, Any}
    packages = NamedTuple[
        (
            name = "python",
            version = "3.11.9",
            build = "h955ad1f_0_cpython",
            channel = "conda-forge",
        ),
        (
            name = "spades",
            version = "4.2.0",
            build = "h5ca1c30_1",
            channel = "bioconda",
        ),
        (
            name = "unicycler",
            version = "0.5.1",
            build = "pyhdfd78af_0",
            channel = "bioconda",
        ),
    ]
    return Mycelia.Rhizomorph._unicycler_toolchain_provenance(
        inventory_reader = () -> packages,
    )
end

function multi_input_fake_autocycler_inventory(;
        build_suffix::AbstractString = "0",
)::Vector{NamedTuple}
    return NamedTuple[
        (
            name = specification.name,
            version = specification.constraint == :present ?
                      "1.0" : specification.version,
            build = "test_$(build_suffix)",
            channel = "bioconda",
        ) for specification in Mycelia.AUTOCYCLER_REQUIRED_PACKAGE_SPECS
    ]
end

function multi_input_fake_autocycler_toolchain(;
        build_suffix::AbstractString = "0",
)::Dict{String, Any}
    return Mycelia._autocycler_toolchain_metadata(
        multi_input_fake_autocycler_inventory(; build_suffix),
    )
end

function multi_input_fake_assembler_result(
        outdir::AbstractString;
        gzip::Bool = false,
        include_graph::Bool = true,
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
                toolchain = multi_input_fake_autocycler_toolchain(),
            )
        end
        toolchain = multi_input_fake_unicycler_toolchain()
        return use_contigs_key ?
               (; contigs = assembly, graph, toolchain) :
               (; assembly, graph, toolchain)
    end
    toolchain = multi_input_fake_unicycler_toolchain()
    return use_contigs_key ? (; contigs = assembly, toolchain) :
           (; assembly, toolchain)
end

function multi_input_without_field(
        result::NamedTuple,
        field::Symbol,
)::NamedTuple
    retained_fields = Tuple(filter(candidate -> candidate != field, keys(result)))
    return NamedTuple{retained_fields}(
        Tuple(getproperty(result, candidate) for candidate in retained_fields),
    )
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

Test.@testset "multi-input hybrid assembly contracts" begin
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

    Test.@testset "external smoke prerequisites fail loud in default CI" begin
        disabled = Mycelia.Rhizomorph._multi_input_hybrid_smoke_prerequisites(
            Dict{String, String}(),
        )
        Test.@test disabled == (;
            run_external = false,
            run_autocycler = false,
        )
        test_throws_message(
            ArgumentError,
            "also requires MYCELIA_RUN_EXTERNAL=true",
        ) do
            Mycelia.Rhizomorph._multi_input_hybrid_smoke_prerequisites(Dict(
                "MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true",
            ))
        end
        test_throws_message(ArgumentError, "missing:") do
            Mycelia.Rhizomorph._multi_input_hybrid_smoke_prerequisites(Dict(
                "MYCELIA_RUN_EXTERNAL" => "true",
            ))
        end

        mktempdir() do temp_dir
            r1 = multi_input_write_fastq(
                joinpath(temp_dir, "smoke-r1.fastq"),
                [multi_input_fastq_record("smoke_pair/1")],
            )
            r2 = multi_input_write_fastq(
                joinpath(temp_dir, "smoke-r2.fastq"),
                [multi_input_fastq_record("smoke_pair/2")],
            )
            long_reads = multi_input_write_fastq(
                joinpath(temp_dir, "smoke-long.fastq"),
                MULTI_INPUT_LONG,
            )
            environment = Dict(
                "MYCELIA_RUN_EXTERNAL" => "true",
                "MYCELIA_HYBRID_SHORT_R1" => r1,
                "MYCELIA_HYBRID_SHORT_R2" => r2,
                "MYCELIA_HYBRID_LONG_READS" => long_reads,
            )
            missing_file_environment = copy(environment)
            missing_file_environment["MYCELIA_HYBRID_LONG_READS"] =
                joinpath(temp_dir, "missing.fastq")
            test_throws_message(ArgumentError, "FASTQ not found") do
                Mycelia.Rhizomorph._multi_input_hybrid_smoke_prerequisites(
                    missing_file_environment,
                )
            end
            for (name, value) in (
                    "MYCELIA_HYBRID_SHORT_TECH" => "nanopore",
                    "MYCELIA_HYBRID_LONG_TECH" => "illumina",
            )
                invalid_environment = copy(environment)
                invalid_environment[name] = value
                test_throws_message(ArgumentError, "must be one of") do
                    Mycelia.Rhizomorph._multi_input_hybrid_smoke_prerequisites(
                        invalid_environment,
                    )
                end
            end

            autocycler_environment = merge(
                environment,
                Dict("MYCELIA_RUN_AUTOCYCLER_POLISHED" => "true"),
            )
            test_throws_message(ArgumentError, "READ_TYPE is required") do
                Mycelia.Rhizomorph._multi_input_hybrid_smoke_prerequisites(
                    autocycler_environment,
                )
            end
            for (name, value, message) in (
                    (
                        "MYCELIA_AUTOCYCLER_READ_TYPE",
                        "invalid",
                        "must be one of",
                    ),
                    (
                        "MYCELIA_AUTOCYCLER_TEST_JOBS",
                        "many",
                        "must be an integer",
                    ),
                    (
                        "MYCELIA_AUTOCYCLER_TEST_JOBS",
                        "5",
                        "must be between",
                    ),
            )
                invalid_environment = merge(
                    autocycler_environment,
                    Dict(
                        "MYCELIA_AUTOCYCLER_READ_TYPE" => "ont_r10",
                        name => value,
                    ),
                )
                test_throws_message(ArgumentError, message) do
                    Mycelia.Rhizomorph._multi_input_hybrid_smoke_prerequisites(
                        invalid_environment,
                    )
                end
            end
            valid_environment = merge(
                autocycler_environment,
                Dict(
                    "MYCELIA_AUTOCYCLER_READ_TYPE" => "ont_r10",
                    "MYCELIA_AUTOCYCLER_TEST_JOBS" => "2",
                ),
            )
            prerequisites =
                Mycelia.Rhizomorph._multi_input_hybrid_smoke_prerequisites(
                    valid_environment,
                )
            Test.@test prerequisites.run_external
            Test.@test prerequisites.run_autocycler
            Test.@test prerequisites.input_paths == (;
                short_r1 = abspath(r1),
                short_r2 = abspath(r2),
                long_reads = abspath(long_reads),
            )
            Test.@test prerequisites.short_read_tech == :illumina
            Test.@test prerequisites.long_read_tech == :nanopore
            Test.@test prerequisites.autocycler_read_type == :ont_r10
            Test.@test prerequisites.autocycler_jobs == 2
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

            lazy_consumptions = Ref(0)
            function lazy_missing_source(path::AbstractString)::Base.Generator
                return (
                    begin
                        lazy_consumptions[] += 1
                        value
                    end for value in (String(path),)
                )
            end
            stale_config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
                output_dir = stale_root,
            )
            test_throws_message(
                ArgumentError,
                "prevent stale hybrid assembly reuse",
            ) do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (
                        lazy_missing_source(missing_inputs[1][1]),
                        lazy_missing_source(missing_inputs[2][1]),
                    ),
                    lazy_missing_source(missing_inputs[3][1]),
                    stale_config,
                    :unicycler;
                    correction_runner = root_validation_correction_runner,
                    assembler_runner = root_validation_assembler_runner,
                )
            end
            Test.@test lazy_consumptions[] == 0

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

    Test.@testset "persistent root reservation is alias-safe" begin
        mktempdir() do temp_dir
            physical_parent = joinpath(temp_dir, "physical")
            mkpath(physical_parent)
            alias_parent = joinpath(temp_dir, "alias")
            symlink(physical_parent, alias_parent)
            physical_output = joinpath(physical_parent, "hybrid-output")
            alias_output = joinpath(alias_parent, "hybrid-output")
            physical_lock = Mycelia.Rhizomorph._multi_input_workflow_lock_path(
                physical_output,
            )
            alias_lock = Mycelia.Rhizomorph._multi_input_workflow_lock_path(
                alias_output,
            )
            Test.@test physical_lock == alias_lock
            Test.@test dirname(physical_lock) == physical_parent
            Test.@test !startswith(physical_lock, "$(physical_output)/")

            failing_lock = Mycelia.Rhizomorph._multi_input_workflow_lock_path(
                joinpath(physical_parent, "failing-workflow"),
            )
            test_throws_message(ErrorException, "synthetic locked failure") do
                Mycelia.Rhizomorph._with_multi_input_workflow_lock(failing_lock) do
                    error("synthetic locked failure")
                end
            end
            Test.@test !ispath(failing_lock)

            first_calls = NamedTuple[]
            first_base_runner = multi_input_fake_correction_runner(first_calls)
            correction_entered = Channel{Nothing}(1)
            release_correction = Channel{Nothing}(1)
            first_call_blocked = Ref(false)
            function blocking_correction_runner(
                    reads::Any,
                    config::Mycelia.Rhizomorph.AssemblyConfig,
            )::NamedTuple
                if !first_call_blocked[]
                    first_call_blocked[] = true
                    put!(correction_entered, nothing)
                    take!(release_correction)
                end
                return first_base_runner(reads, config)
            end
            function first_assembler_runner(
                    inputs::Any,
                    outdir::AbstractString,
            )::NamedTuple
                Test.@test isfile(physical_lock)
                Test.@test dirname(String(outdir)) == physical_output
                Test.@test String(outdir) != physical_lock
                return multi_input_fake_assembler_result(outdir)
            end
            first_config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
                output_dir = physical_output,
            )
            first_task = @async Mycelia.Rhizomorph._assemble_paired_short_long(
                (MULTI_INPUT_R1, MULTI_INPUT_R2),
                MULTI_INPUT_LONG,
                first_config,
                :unicycler;
                correction_runner = blocking_correction_runner,
                assembler_runner = first_assembler_runner,
            )
            wait_status = Base.timedwait(
                () -> isready(correction_entered) || istaskdone(first_task),
                30.0,
            )
            Test.@test wait_status == :ok
            if !isready(correction_entered)
                fetch(first_task)
            end
            take!(correction_entered)

            second_correction_calls = Ref(0)
            second_assembler_calls = Ref(0)
            function second_correction_runner(
                    reads::Any,
                    config::Mycelia.Rhizomorph.AssemblyConfig,
            )::NamedTuple
                second_correction_calls[] += 1
                error("second correction must not run")
            end
            function second_assembler_runner(
                    inputs::Any,
                    outdir::AbstractString,
            )::NamedTuple
                second_assembler_calls[] += 1
                error("second assembler must not run")
            end
            second_config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
                output_dir = alias_output,
            )
            try
                test_throws_message(ArgumentError, "already reserved") do
                    Mycelia.Rhizomorph._assemble_paired_short_long(
                        (MULTI_INPUT_R1, MULTI_INPUT_R2),
                        MULTI_INPUT_LONG,
                        second_config,
                        :unicycler;
                        correction_runner = second_correction_runner,
                        assembler_runner = second_assembler_runner,
                    )
                end
            finally
                put!(release_correction, nothing)
            end
            first_result = fetch(first_task)
            Test.@test length(first_calls) == 3
            Test.@test second_correction_calls[] == 0
            Test.@test second_assembler_calls[] == 0
            Test.@test first_result.assembly_stats["output_dir"] == physical_output
            Test.@test isempty(first_result.assembly_stats["retained_cleanup_files"])
            Test.@test isempty(first_result.assembly_stats["retained_cleanup_roots"])
            Test.@test !ispath(physical_lock)

            lifecycle_output = joinpath(physical_parent, "lifecycle-output")
            lifecycle_lock =
                Mycelia.Rhizomorph._multi_input_workflow_lock_path(
                    lifecycle_output,
                )
            lock_observations = Bool[]
            function lock_observing_source(values::Any)::Base.Generator
                return (
                    begin
                        push!(lock_observations, isfile(lifecycle_lock))
                        value
                    end for value in values
                )
            end
            lifecycle_calls = NamedTuple[]
            lifecycle_result =
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (
                        lock_observing_source(MULTI_INPUT_R1),
                        lock_observing_source(MULTI_INPUT_R2),
                    ),
                    lock_observing_source(MULTI_INPUT_LONG),
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(
                        output_dir = lifecycle_output,
                    ),
                    :unicycler;
                    correction_runner = multi_input_fake_correction_runner(
                        lifecycle_calls,
                    ),
                    assembler_runner = (inputs, outdir) -> begin
                        Test.@test isfile(lifecycle_lock)
                        multi_input_fake_assembler_result(outdir)
                    end,
                )
            Test.@test length(lifecycle_calls) == 3
            Test.@test !isempty(lock_observations)
            Test.@test all(lock_observations)
            Test.@test lifecycle_result.assembly_stats["output_dir"] ==
                       lifecycle_output
            Test.@test !ispath(lifecycle_lock)

            retarget_parent = joinpath(temp_dir, "retarget-physical")
            mkpath(retarget_parent)
            captured_output = joinpath(
                physical_parent,
                "retargeted-hybrid-output",
            )
            retargeted_alias_output = joinpath(
                alias_parent,
                "retargeted-hybrid-output",
            )
            captured_lock =
                Mycelia.Rhizomorph._multi_input_workflow_lock_path(
                    captured_output,
                )
            function retargeting_workflow_lock(
                    action::Function,
                    lock_path::AbstractString,
            )::Any
                Test.@test String(lock_path) == captured_lock
                return Mycelia.Rhizomorph._with_multi_input_workflow_lock(
                    lock_path,
                ) do
                    rm(alias_parent)
                    symlink(retarget_parent, alias_parent)
                    return action()
                end
            end
            retarget_correction_calls = NamedTuple[]
            retarget_result =
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (MULTI_INPUT_R1, MULTI_INPUT_R2),
                    MULTI_INPUT_LONG,
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(
                        output_dir = retargeted_alias_output,
                    ),
                    :unicycler;
                    correction_runner = multi_input_fake_correction_runner(
                        retarget_correction_calls,
                    ),
                    assembler_runner = (inputs, outdir) -> begin
                        Test.@test String(outdir) == joinpath(
                            captured_output,
                            "assembler_unicycler",
                        )
                        multi_input_fake_assembler_result(outdir)
                    end,
                    workflow_lock_runner = retargeting_workflow_lock,
                )
            Test.@test length(retarget_correction_calls) == 3
            Test.@test retarget_result.assembly_stats["output_dir"] ==
                       captured_output
            Test.@test isdir(captured_output)
            Test.@test !ispath(joinpath(
                retarget_parent,
                "retargeted-hybrid-output",
            ))
            Test.@test !ispath(captured_lock)

            hierarchy_parent = joinpath(temp_dir, "hierarchy-parent")
            hierarchy_child = joinpath(hierarchy_parent, "nested")
            parent_lock =
                Mycelia.Rhizomorph._multi_input_workflow_lock_path(
                    hierarchy_parent,
                )
            child_lock =
                Mycelia.Rhizomorph._multi_input_workflow_lock_path(
                    hierarchy_child,
                )
            hierarchy_correction_calls = Ref(0)
            hierarchy_assembler_calls = Ref(0)
            function hierarchy_correction_runner(
                    reads::Any,
                    config::Mycelia.Rhizomorph.AssemblyConfig,
            )::NamedTuple
                hierarchy_correction_calls[] += 1
                error("hierarchical correction must not run")
            end
            function hierarchy_assembler_runner(
                    inputs::Any,
                    outdir::AbstractString,
            )::NamedTuple
                hierarchy_assembler_calls[] += 1
                error("hierarchical assembler must not run")
            end
            Mycelia.Rhizomorph._with_multi_input_workflow_lock(parent_lock) do
                test_throws_message(
                    ArgumentError,
                    "active ancestor workflow reservation",
                ) do
                    Mycelia.Rhizomorph._assemble_paired_short_long(
                        (MULTI_INPUT_R1, MULTI_INPUT_R2),
                        MULTI_INPUT_LONG,
                        Mycelia.Rhizomorph.UnicyclerHybridConfig(
                            output_dir = hierarchy_child,
                        ),
                        :unicycler;
                        correction_runner = hierarchy_correction_runner,
                        assembler_runner = hierarchy_assembler_runner,
                    )
                end
            end
            Mycelia.Rhizomorph._with_multi_input_workflow_lock(child_lock) do
                test_throws_message(
                    ArgumentError,
                    "active descendant workflow reservation",
                ) do
                    Mycelia.Rhizomorph._assemble_paired_short_long(
                        (MULTI_INPUT_R1, MULTI_INPUT_R2),
                        MULTI_INPUT_LONG,
                        Mycelia.Rhizomorph.UnicyclerHybridConfig(
                            output_dir = hierarchy_parent,
                        ),
                        :unicycler;
                        correction_runner = hierarchy_correction_runner,
                        assembler_runner = hierarchy_assembler_runner,
                    )
                end
            end
            Test.@test hierarchy_correction_calls[] == 0
            Test.@test hierarchy_assembler_calls[] == 0
            Test.@test !ispath(parent_lock)
            Test.@test !ispath(child_lock)

            sibling_a = joinpath(temp_dir, "hierarchy-sibling-a")
            sibling_b = joinpath(temp_dir, "hierarchy-sibling-b")
            sibling_a_lock =
                Mycelia.Rhizomorph._multi_input_workflow_lock_path(sibling_a)
            sibling_b_lock =
                Mycelia.Rhizomorph._multi_input_workflow_lock_path(sibling_b)
            Mycelia.Rhizomorph._with_multi_input_workflow_lock(
                sibling_a_lock,
            ) do
                Mycelia.Rhizomorph._with_multi_input_workflow_lock(
                    sibling_b_lock,
                ) do
                    sibling_reservation =
                        Mycelia.Rhizomorph._require_exclusive_multi_input_workflow_domain(
                            sibling_b,
                            sibling_b_lock,
                        )
                    Test.@test sibling_reservation === nothing
                end
            end
        end
    end

    Test.@testset "raw and corrected mate validation" begin
        Test.@test Mycelia.Rhizomorph._validate_paired_reads(
            MULTI_INPUT_R1,
            MULTI_INPUT_R2,
            "input",
        ) == 2

        distinct_equal_r1 = [multi_input_fastq_record("equal_pair")]
        distinct_equal_r2 = [multi_input_fastq_record("equal_pair")]
        Test.@test distinct_equal_r1[1] == distinct_equal_r2[1]
        Test.@test distinct_equal_r1[1] !== distinct_equal_r2[1]
        Test.@test Mycelia.Rhizomorph._validate_distinct_in_memory_read_records(
            distinct_equal_r1,
            distinct_equal_r2,
            MULTI_INPUT_LONG,
        ) === nothing
        distinct_equal_calls = NamedTuple[]
        function distinct_equal_assembler_runner(
                inputs::Any,
                outdir::AbstractString,
        )::NamedTuple
            return multi_input_fake_assembler_result(outdir)
        end
        distinct_equal_result =
            Mycelia.Rhizomorph._assemble_paired_short_long(
                (distinct_equal_r1, distinct_equal_r2),
                MULTI_INPUT_LONG,
                Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                :unicycler;
                correction_runner = multi_input_fake_correction_runner(
                    distinct_equal_calls,
                ),
                assembler_runner = distinct_equal_assembler_runner,
            )
        Test.@test length(distinct_equal_calls) == 3
        expected_distinct_equal_counts = Dict(
            "short_r1" => 1,
            "short_r2" => 1,
            "long_reads" => length(MULTI_INPUT_LONG),
        )
        Test.@test distinct_equal_result.assembly_stats["input_read_counts"] ==
                   expected_distinct_equal_counts

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
        casava_r1 = [
            multi_input_fastq_record("casava_pair 1:N:0:ACGT"),
        ]
        casava_r2 = [
            multi_input_fastq_record("casava_pair 2:Y:0:ACGT"),
        ]
        Test.@test Mycelia.Rhizomorph._validate_paired_reads(
            casava_r1,
            casava_r2,
            "input",
        ) == 1
        noted_r1 = [
            multi_input_fastq_record(
                "noted_pair/1 note 2:N:0:ACGT",
            ),
        ]
        noted_r2 = [
            multi_input_fastq_record(
                "noted_pair/2 note 1:N:0:ACGT",
            ),
        ]
        Test.@test Mycelia.Rhizomorph._validate_paired_reads(
            noted_r1,
            noted_r2,
            "input",
        ) == 1
        malformed_casava_r1 = [
            multi_input_fastq_record(
                "anchored_pair/1 2:N:0:ACGT:trailing",
            ),
        ]
        malformed_casava_r2 = [
            multi_input_fastq_record(
                "anchored_pair/2 1:N:0:ACGT:trailing",
            ),
        ]
        Test.@test Mycelia.Rhizomorph._validate_paired_reads(
            malformed_casava_r1,
            malformed_casava_r2,
            "input",
        ) == 1

        invalid_casava_pairs = (
            (
                message = "invalid explicit mate roles",
                r1 = [multi_input_fastq_record("casava_pair 2:N:0:ACGT")],
                r2 = [multi_input_fastq_record("casava_pair 1:N:0:ACGT")],
            ),
            (
                message = "conflicting explicit mate roles",
                r1 = [
                    multi_input_fastq_record(
                        "conflict_pair/1 2:N:0:ACGT",
                    ),
                ],
                r2 = [
                    multi_input_fastq_record(
                        "conflict_pair/2 2:N:0:ACGT",
                    ),
                ],
            ),
        )
        workflow_configs = (
            (
                workflow = :unicycler,
                config = Mycelia.Rhizomorph.UnicyclerHybridConfig(),
            ),
            (
                workflow = :autocycler_polished,
                config = Mycelia.Rhizomorph.AutocyclerPolishConfig(),
            ),
        )
        mktempdir() do casava_dir
            valid_r1_path = multi_input_write_fastq(
                joinpath(casava_dir, "valid-r1.fastq"),
                casava_r1,
            )
            valid_r2_path = multi_input_write_fastq(
                joinpath(casava_dir, "valid-r2.fastq"),
                casava_r2,
            )
            Test.@test Mycelia.Rhizomorph._validate_paired_reads(
                [valid_r1_path],
                [valid_r2_path],
                "input",
            ) == 1

            for (workflow_index, workflow_config) in enumerate(workflow_configs)
                for (pair_index, invalid_pair) in enumerate(invalid_casava_pairs)
                    r1_path = multi_input_write_fastq(
                        joinpath(
                            casava_dir,
                            "invalid-$(workflow_index)-$(pair_index)-r1.fastq",
                        ),
                        invalid_pair.r1,
                    )
                    r2_path = multi_input_write_fastq(
                        joinpath(
                            casava_dir,
                            "invalid-$(workflow_index)-$(pair_index)-r2.fastq",
                        ),
                        invalid_pair.r2,
                    )
                    correction_calls = Ref(0)
                    assembler_calls = Ref(0)
                    correction_runner = function (
                            reads::Any,
                            correction_config::Mycelia.Rhizomorph.AssemblyConfig,
                    )
                        correction_calls[] += 1
                        error("correction must not run")
                    end
                    assembler_runner = function (
                            inputs::Any,
                            outdir::AbstractString,
                    )
                        assembler_calls[] += 1
                        error("assembler must not run")
                    end
                    test_throws_message(ArgumentError, invalid_pair.message) do
                        Mycelia.Rhizomorph._assemble_paired_short_long(
                            (r1_path, r2_path),
                            MULTI_INPUT_LONG,
                            workflow_config.config,
                            workflow_config.workflow;
                            correction_runner,
                            assembler_runner,
                        )
                    end
                    Test.@test correction_calls[] == 0
                    Test.@test assembler_calls[] == 0
                end
            end
        end
        for workflow_config in workflow_configs
            reconstructed_calls = NamedTuple[]
            correction_runner = multi_input_fake_correction_runner(
                reconstructed_calls;
                transform = function (
                        index::Int,
                        records::Vector{FASTX.FASTQ.Record},
                )
                    index == 2 || return records
                    return FASTX.FASTQ.Record[
                        multi_input_fastq_record(
                            String(FASTX.identifier(record)),
                            FASTX.sequence(String, record),
                        ) for record in records
                    ]
                end,
            )
            assembler_calls = Ref(0)
            assembler_runner = function (
                    inputs::Mycelia.Rhizomorph._CorrectedPairedShortLong,
                    outdir::AbstractString,
            )
                assembler_calls[] += 1
                corrected_r1_records = multi_input_collect_fastq([
                    inputs.short_r1.path,
                ])
                corrected_r2_records = multi_input_collect_fastq([
                    inputs.short_r2.path,
                ])
                Test.@test String(FASTX.description(only(corrected_r1_records))) ==
                           "casava_pair 1:N:0:ACGT"
                Test.@test String(FASTX.description(only(corrected_r2_records))) ==
                           "casava_pair"
                Test.@test String(FASTX.identifier(only(corrected_r1_records))) ==
                           "casava_pair"
                Test.@test String(FASTX.identifier(only(corrected_r2_records))) ==
                           "casava_pair"
                return multi_input_fake_assembler_result(
                    outdir;
                    include_polishing_artifacts =
                        workflow_config.workflow == :autocycler_polished,
                )
            end
            result = Mycelia.Rhizomorph._assemble_paired_short_long(
                (casava_r1, casava_r2),
                MULTI_INPUT_LONG,
                workflow_config.config,
                workflow_config.workflow;
                correction_runner,
                assembler_runner,
            )
            Test.@test result.assembly_stats["workflow"] ==
                       String(workflow_config.workflow)
            Test.@test length(reconstructed_calls) == 3
            Test.@test assembler_calls[] == 1
            Test.@test all(
                call -> !isfile(call.corrected_fastq),
                reconstructed_calls,
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

    Test.@testset "path-backed inputs are strictly FASTQ" begin
        mktempdir() do temp_dir
            valid_r1 = multi_input_write_fastq(
                joinpath(temp_dir, "valid-r1.fastq"),
                [multi_input_fastq_record("strict_pair/1")],
            )
            valid_r2 = multi_input_write_fastq(
                joinpath(temp_dir, "valid-r2.fastq"),
                [multi_input_fastq_record("strict_pair/2")],
            )
            valid_long = multi_input_write_fastq(
                joinpath(temp_dir, "valid-long.fastq"),
                MULTI_INPUT_LONG,
            )
            fasta_input = multi_input_write_fasta(
                joinpath(temp_dir, "reads.fasta");
                records = ["not_fastq" => "ACGT"],
            )
            mixed_input = joinpath(temp_dir, "mixed.fastq")
            write(mixed_input, read(valid_long))
            open(mixed_input, "a") do stream
                write(stream, ">mixed_fasta\nACGT\n")
            end

            invalid_sources = (
                (label = "short_r1", path = fasta_input),
                (label = "short_r2", path = mixed_input),
                (label = "long_reads", path = fasta_input),
            )
            correction_calls = Ref(0)
            assembler_calls = Ref(0)
            for (case_index, invalid_source) in enumerate(invalid_sources)
                short_r1 = invalid_source.label == "short_r1" ?
                           invalid_source.path : valid_r1
                short_r2 = invalid_source.label == "short_r2" ?
                           invalid_source.path : valid_r2
                selected_long = invalid_source.label == "long_reads" ?
                                invalid_source.path : valid_long
                output_dir = joinpath(temp_dir, "invalid-$(case_index)")
                test_throws_message(ArgumentError, "FASTQ") do
                    Mycelia.Rhizomorph._assemble_paired_short_long(
                        (short_r1, short_r2),
                        selected_long,
                        Mycelia.Rhizomorph.UnicyclerHybridConfig(; output_dir),
                        :unicycler;
                        correction_runner = (reads, config) -> begin
                            correction_calls[] += 1
                            error("correction must not run")
                        end,
                        assembler_runner = (inputs, outdir) -> begin
                            assembler_calls[] += 1
                            error("assembler must not run")
                        end,
                    )
                end
                Test.@test !ispath(output_dir)
            end
            Test.@test correction_calls[] == 0
            Test.@test assembler_calls[] == 0

            fasta_record = FASTX.FASTA.Record("not_fastq", "ACGT")
            test_throws_message(ArgumentError, "is not FASTQ") do
                Mycelia.Rhizomorph._multi_input_source_content_identity(
                    [fasta_record],
                    "long_reads",
                )
            end
        end
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
            source_contract = result.assembly_stats[
                "read_content_provenance"
            ]["source_inputs"]
            expected_paths = Dict(
                "short_r1" => r1_gzip,
                "short_r2" => r2_gzip,
                "long_reads" => long_gzip,
            )
            for (label, expected_path) in expected_paths
                identity = source_contract[label]
                Test.@test identity["kind"] == "path_set"
                Test.@test identity["source_count"] == 1
                source = only(identity["sources"])
                Test.@test source["path"] == abspath(expected_path)
                Test.@test source["canonical_path"] == realpath(expected_path)
                Test.@test source["size_bytes"] == filesize(expected_path)
                Test.@test occursin(r"^[0-9a-f]{64}$", source["sha256"])
                Test.@test occursin(r"^[0-9a-f]{64}$", identity["sha256"])
            end
            Test.@test all(call -> !isfile(call.corrected_fastq), correction_calls)
        end
    end

    Test.@testset "lightweight Stage-1 cleanup ownership" begin
        Test.@test fieldnames(Mycelia.Rhizomorph._Stage1CleanupToken) ==
                   (:corrected_fastq, :ephemeral, :identity)
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

        counted_fastq = Ref("")
        function counted_runner(
                reads::Any,
                config::Mycelia.Rhizomorph.AssemblyConfig,
        )::NamedTuple
            counted_fastq[] = multi_input_write_fastq(
                joinpath(something(config.output_dir), "counted.fastq"),
                MULTI_INPUT_LONG,
            )
            return (;
                corrected_fastq = counted_fastq[],
                corrected_read_count = 1,
                # The route, not this untrusted bit, owns cleanup authority.
                ephemeral = false,
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
            Mycelia.Rhizomorph._Stage1CleanupToken(counted_fastq[], true),
        ]
        Mycelia.Rhizomorph._cleanup_multi_input_stages!(cleanup_tokens)
        Test.@test !isfile(counted_fastq[])

        interrupted_fastq = multi_input_write_fastq(
            tempname() * ".fastq",
            MULTI_INPUT_LONG,
        )
        interrupted_tokens = [
            Mycelia.Rhizomorph._Stage1CleanupToken(interrupted_fastq, true),
        ]
        stage_interrupt = InterruptException()
        caught_stage_interrupt = try
            Mycelia.Rhizomorph._cleanup_multi_input_stages!(
                interrupted_tokens;
                remover = path -> throw(stage_interrupt),
            )
            nothing
        catch caught
            caught
        end
        Test.@test caught_stage_interrupt === stage_interrupt
        Test.@test isfile(interrupted_fastq)
        rm(interrupted_fastq; force = true)

        mktempdir() do temp_dir
            workflow_root = joinpath(temp_dir, "workflow")
            stage_dir = joinpath(workflow_root, "corrected", "long_reads")
            external_dir = joinpath(temp_dir, "external")
            corrected_fastq = multi_input_write_fastq(
                joinpath(stage_dir, "corrected.fastq"),
                MULTI_INPUT_LONG,
            )
            external_fastq = multi_input_write_fastq(
                joinpath(external_dir, "corrected.fastq"),
                MULTI_INPUT_LONG,
            )
            retargeted_tokens = [
                Mycelia.Rhizomorph._Stage1CleanupToken(
                    corrected_fastq,
                    true,
                ),
            ]
            rm(stage_dir; recursive = true)
            symlink(external_dir, stage_dir)
            retained = Test.@test_logs (
                :warn,
                r"refusing corrected FASTQ cleanup after path identity changed",
            ) min_level=Logging.Warn begin
                Mycelia.Rhizomorph._cleanup_multi_input_stages!(
                    retargeted_tokens,
                )
            end
            Test.@test retained == [corrected_fastq]
            Test.@test isfile(external_fastq)
            Test.@test read(external_fastq) == read(
                joinpath(stage_dir, "corrected.fastq"),
            )
        end

        mktempdir() do temp_dir
            corrected_fastq = multi_input_write_fastq(
                joinpath(temp_dir, "corrected.fastq"),
                MULTI_INPUT_LONG,
            )
            replacement_fastq = multi_input_write_fastq(
                joinpath(temp_dir, "replacement.fastq"),
                [multi_input_fastq_record("replacement")],
            )
            replacement_content = read(replacement_fastq)
            replaced_tokens = [
                Mycelia.Rhizomorph._Stage1CleanupToken(
                    corrected_fastq,
                    true,
                ),
            ]
            rm(corrected_fastq)
            mv(replacement_fastq, corrected_fastq)
            retained = Test.@test_logs (
                :warn,
                r"refusing corrected FASTQ cleanup after path identity changed",
            ) min_level=Logging.Warn begin
                Mycelia.Rhizomorph._cleanup_multi_input_stages!(
                    replaced_tokens,
                )
            end
            Test.@test retained == [corrected_fastq]
            Test.@test read(corrected_fastq) == replacement_content
        end

        interrupted_root = mktempdir()
        root_interrupt = InterruptException()
        caught_root_interrupt = try
            Mycelia.Rhizomorph._cleanup_multi_input_root!(
                interrupted_root;
                remover = path -> throw(root_interrupt),
            )
            nothing
        catch caught
            caught
        end
        Test.@test caught_root_interrupt === root_interrupt
        Test.@test isdir(interrupted_root)
        rm(interrupted_root; recursive = true, force = true)

        mktempdir() do temp_dir
            workflow_root = joinpath(temp_dir, "ephemeral-root")
            external_root = joinpath(temp_dir, "external-root")
            mkpath(workflow_root)
            mkpath(external_root)
            external_marker = joinpath(external_root, "preserve.txt")
            write(external_marker, "preserve\n")
            rm(workflow_root; recursive = true)
            symlink(external_root, workflow_root)
            retained = Test.@test_logs (
                :warn,
                r"refusing ephemeral root cleanup after path identity changed",
            ) min_level=Logging.Warn begin
                Mycelia.Rhizomorph._cleanup_multi_input_root!(workflow_root)
            end
            Test.@test retained
            Test.@test read(external_marker, String) == "preserve\n"
        end

        mktempdir() do temp_dir
            workflow_root = joinpath(temp_dir, "ephemeral-root")
            replacement_root = joinpath(temp_dir, "replacement-root")
            mkpath(workflow_root)
            expected_identity =
                Mycelia.Rhizomorph._workflow_path_identity(
                    workflow_root,
                    "test workflow root",
                )
            mkpath(replacement_root)
            replacement_marker = joinpath(replacement_root, "preserve.txt")
            write(replacement_marker, "replacement\n")
            rm(workflow_root; recursive = true)
            mv(replacement_root, workflow_root)
            retained = Test.@test_logs (
                :warn,
                r"refusing ephemeral root cleanup after path identity changed",
            ) min_level=Logging.Warn begin
                Mycelia.Rhizomorph._cleanup_multi_input_root!(
                    workflow_root;
                    expected_identity,
                )
            end
            Test.@test retained
            Test.@test read(
                joinpath(workflow_root, "preserve.txt"),
                String,
            ) == "replacement\n"
        end

        malformed_stages = [
            (
                "missing corrected_fastq",
                (reads, config) -> (; ephemeral = true),
            ),
            (
                "regular, non-symlink FASTQ",
                (reads, config) -> (;
                    corrected_fastq = joinpath(
                        something(config.output_dir),
                        "missing.fastq",
                    ),
                    corrected_read_count = 1,
                    ephemeral = true,
                ),
            ),
            (
                "no non-empty corrected FASTQ",
                function (reads, config)
                    empty_fastq = joinpath(
                        something(config.output_dir),
                        "empty.fastq",
                    )
                    touch(empty_fastq)
                    return (;
                        corrected_fastq = empty_fastq,
                        corrected_read_count = 1,
                        ephemeral = true,
                    )
                end,
            ),
            (
                "missing corrected_read_count",
                function (reads, config)
                    valid_fastq = multi_input_write_fastq(
                        joinpath(
                            something(config.output_dir),
                            "missing-count.fastq",
                        ),
                        MULTI_INPUT_LONG,
                    )
                    return (; corrected_fastq = valid_fastq, ephemeral = false)
                end,
            ),
            (
                "produced 0 corrected reads",
                function (reads, config)
                    valid_fastq = multi_input_write_fastq(
                        joinpath(something(config.output_dir), "zero.fastq"),
                        MULTI_INPUT_LONG,
                    )
                    return (;
                        corrected_fastq = valid_fastq,
                        corrected_read_count = 0,
                        ephemeral = false,
                    )
                end,
            ),
            (
                "reported 7 corrected reads",
                function (reads, config)
                    valid_fastq = multi_input_write_fastq(
                        joinpath(
                            something(config.output_dir),
                            "wrong-count.fastq",
                        ),
                        MULTI_INPUT_LONG,
                    )
                    return (;
                        corrected_fastq = valid_fastq,
                        corrected_read_count = 7,
                        ephemeral = false,
                    )
                end,
            ),
            (
                "correction output is not FASTQ",
                function (reads, config)
                    fasta_output = joinpath(
                        something(config.output_dir),
                        "not-fastq.fasta",
                    )
                    write(fasta_output, ">not_fastq\nACGT\n")
                    return (;
                        corrected_fastq = fasta_output,
                        corrected_read_count = 1,
                        ephemeral = false,
                    )
                end,
            ),
        ]
        for (message, bad_runner) in malformed_stages
            tokens = Mycelia.Rhizomorph._Stage1CleanupToken[]
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
    end

    Test.@testset "correction output authority is route-owned" begin
        mktempdir() do temp_dir
            source_r1 = multi_input_write_fastq(
                joinpath(temp_dir, "source-r1.fastq"),
                [multi_input_fastq_record("owned_pair/1")],
            )
            source_r2 = multi_input_write_fastq(
                joinpath(temp_dir, "source-r2.fastq"),
                [multi_input_fastq_record("owned_pair/2")],
            )
            source_long = multi_input_write_fastq(
                joinpath(temp_dir, "source-long.fastq"),
                MULTI_INPUT_LONG,
            )
            sources = (source_r1, source_r2, source_long)
            external_fastq = multi_input_write_fastq(
                joinpath(temp_dir, "external.fastq"),
                [multi_input_fastq_record("owned_pair/1")],
            )

            external_runner = (reads, config) -> (;
                corrected_fastq = external_fastq,
                corrected_read_count = 1,
                ephemeral = true,
            )
            test_throws_message(ErrorException, "outside its reserved stage") do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    ((source_r1,), (source_r2,)),
                    (source_long,),
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                    :unicycler;
                    correction_runner = external_runner,
                    assembler_runner = (inputs, outdir) ->
                        error("assembler must not run"),
                )
            end
            Test.@test isfile(external_fastq)
            Test.@test all(isfile, sources)

            function symlink_runner(reads::Any, config::Any)::NamedTuple
                corrected_fastq = joinpath(
                    something(config.output_dir),
                    "symlink.fastq",
                )
                symlink(external_fastq, corrected_fastq)
                return (;
                    corrected_fastq,
                    corrected_read_count = 1,
                    ephemeral = true,
                )
            end
            test_throws_message(ErrorException, "regular, non-symlink FASTQ") do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    ((source_r1,), (source_r2,)),
                    (source_long,),
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                    :unicycler;
                    correction_runner = symlink_runner,
                    assembler_runner = (inputs, outdir) ->
                        error("assembler must not run"),
                )
            end
            Test.@test isfile(external_fastq)
            Test.@test all(isfile, sources)

            stage_escape_root = joinpath(temp_dir, "stage-escape-root")
            escaped_stage = joinpath(
                stage_escape_root,
                "corrected",
                "long_reads",
            )
            external_stage = joinpath(temp_dir, "external-stage")
            mkpath(dirname(escaped_stage))
            mkpath(external_stage)
            symlink(external_stage, escaped_stage)
            stage_runner_called = Ref(false)
            test_throws_message(
                ArgumentError,
                "contains a symbolic-link path component",
            ) do
                Mycelia.Rhizomorph._correct_read_set!(
                    Mycelia.Rhizomorph._Stage1CleanupToken[],
                    MULTI_INPUT_LONG,
                    :nanopore,
                    :long_reads,
                    (; k = 13, strategy = :scalable),
                    stage_escape_root,
                    false,
                    (reads, config) -> begin
                        stage_runner_called[] = true
                        error("escaped correction runner must not execute")
                    end,
                )
            end
            Test.@test !stage_runner_called[]
            Test.@test isempty(readdir(external_stage))

            postcheck_root = joinpath(temp_dir, "stage-postcheck-root")
            postcheck_external = joinpath(temp_dir, "stage-postcheck-external")
            mkpath(postcheck_root)
            mkpath(postcheck_external)
            postcheck_external_fastq = joinpath(
                postcheck_external,
                "corrected.fastq",
            )
            function retargeting_stage_runner(
                    reads::Any,
                    config::Any,
            )::NamedTuple
                stage_dir = something(config.output_dir)
                rm(stage_dir; recursive = true)
                symlink(postcheck_external, stage_dir)
                multi_input_write_fastq(
                    joinpath(stage_dir, "corrected.fastq"),
                    MULTI_INPUT_LONG,
                )
                return (;
                    corrected_fastq = joinpath(stage_dir, "corrected.fastq"),
                    corrected_read_count = 1,
                    ephemeral = true,
                )
            end
            test_throws_message(
                ArgumentError,
                "contains a symbolic-link path component",
            ) do
                Mycelia.Rhizomorph._correct_read_set!(
                    Mycelia.Rhizomorph._Stage1CleanupToken[],
                    MULTI_INPUT_LONG,
                    :nanopore,
                    :long_reads,
                    (; k = 13, strategy = :scalable),
                    postcheck_root,
                    false,
                    retargeting_stage_runner,
                )
            end
            Test.@test isfile(postcheck_external_fastq)

            replacement_root = joinpath(temp_dir, "stage-replacement-root")
            mkpath(replacement_root)
            function replacing_stage_runner(
                    reads::Any,
                    config::Any,
            )::NamedTuple
                stage_dir = something(config.output_dir)
                replacement_stage = joinpath(
                    dirname(stage_dir),
                    "replacement-stage",
                )
                mkpath(replacement_stage)
                corrected_fastq = multi_input_write_fastq(
                    joinpath(replacement_stage, "corrected.fastq"),
                    MULTI_INPUT_LONG,
                )
                rm(stage_dir; recursive = true)
                mv(replacement_stage, stage_dir)
                return (;
                    corrected_fastq = joinpath(
                        stage_dir,
                        basename(corrected_fastq),
                    ),
                    corrected_read_count = 1,
                    ephemeral = true,
                )
            end
            test_throws_message(
                ErrorException,
                "correction stage changed filesystem identity",
            ) do
                Mycelia.Rhizomorph._correct_read_set!(
                    Mycelia.Rhizomorph._Stage1CleanupToken[],
                    MULTI_INPUT_LONG,
                    :nanopore,
                    :long_reads,
                    (; k = 13, strategy = :scalable),
                    replacement_root,
                    false,
                    replacing_stage_runner,
                )
            end

            function source_alias_runner(reads::Any, config::Any)::NamedTuple
                source_path = only(reads)
                corrected_fastq = joinpath(
                    something(config.output_dir),
                    "source-alias.fastq",
                )
                Base.hardlink(source_path, corrected_fastq)
                return (;
                    corrected_fastq,
                    corrected_read_count = 1,
                    ephemeral = false,
                )
            end
            test_throws_message(ErrorException, "aliases an input source") do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    ((source_r1,), (source_r2,)),
                    (source_long,),
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                    :unicycler;
                    correction_runner = source_alias_runner,
                    assembler_runner = (inputs, outdir) ->
                        error("assembler must not run"),
                )
            end
            Test.@test all(isfile, sources)

            first_corrected = Ref("")
            correction_index = Ref(0)
            function corrected_alias_runner(
                    reads::Any,
                    config::Any,
            )::NamedTuple
                correction_index[] += 1
                corrected_fastq = joinpath(
                    something(config.output_dir),
                    "corrected.fastq",
                )
                if correction_index[] == 1
                    first_corrected[] = multi_input_write_fastq(
                        corrected_fastq,
                        multi_input_collect_fastq(reads),
                    )
                else
                    Base.hardlink(first_corrected[], corrected_fastq)
                end
                return (;
                    corrected_fastq,
                    corrected_read_count = 1,
                    ephemeral = false,
                )
            end
            test_throws_message(ErrorException, "aliases another corrected read set") do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    ((source_r1,), (source_r2,)),
                    (source_long,),
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                    :unicycler;
                    correction_runner = corrected_alias_runner,
                    assembler_runner = (inputs, outdir) ->
                        error("assembler must not run"),
                )
            end
            Test.@test correction_index[] == 2
            Test.@test all(isfile, sources)

            persistent_root = joinpath(temp_dir, "persistent-lie")
            persistent_calls = NamedTuple[]
            base_runner = multi_input_fake_correction_runner(persistent_calls)
            lying_runner = (reads, config) -> merge(
                base_runner(reads, config),
                (; ephemeral = true),
            )
            persistent_result = Mycelia.Rhizomorph._assemble_paired_short_long(
                ((source_r1,), (source_r2,)),
                (source_long,),
                Mycelia.Rhizomorph.UnicyclerHybridConfig(
                    output_dir = persistent_root,
                ),
                :unicycler;
                correction_runner = lying_runner,
                assembler_runner = (inputs, outdir) ->
                    multi_input_fake_assembler_result(outdir),
            )
            Test.@test all(
                isfile,
                values(persistent_result.assembly_stats["corrected_fastqs"]),
            )
        end
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
        content_provenance =
            result.assembly_stats["read_content_provenance"]
        source_contract = content_provenance["source_inputs"]
        Test.@test source_contract["schema"] ==
                   "mycelia-paired-short-long-input-content-v1"
        for label in ("short_r1", "short_r2", "long_reads")
            identity = source_contract[label]
            Test.@test identity["kind"] == "in_memory_fastq"
            Test.@test identity["record_count"] ==
                       result.assembly_stats["input_read_counts"][label]
            Test.@test occursin(r"^[0-9a-f]{64}$", identity["sha256"])
            Test.@test occursin(
                r"^[0-9a-f]{64}$",
                identity["identifier_sha256"],
            )
        end
        corrected_contract = content_provenance["corrected_fastqs"]
        Test.@test corrected_contract["schema"] ==
                   "mycelia-corrected-fastq-content-v1"
        for label in ("short_r1", "short_r2", "long_reads")
            identity = corrected_contract[label]
            Test.@test identity["kind"] == "path"
            Test.@test isabspath(identity["path"])
            Test.@test identity["size_bytes"] > 0
            Test.@test occursin(r"^[0-9a-f]{64}$", identity["sha256"])
        end
        Test.@test result.assembly_stats["toolchain"] ==
                   multi_input_fake_unicycler_toolchain()
        correction_provenance = result.assembly_stats["correction_provenance"]
        Test.@test Set(keys(correction_provenance)) ==
                   Set(["short_r1", "short_r2", "long_reads"])
        expected_technologies = Dict(
            "short_r1" => "illumina",
            "short_r2" => "illumina",
            "long_reads" => "nanopore",
        )
        for (label, technology) in expected_technologies
            provenance = correction_provenance[label]
            Test.@test provenance["status"] == "unavailable"
            Test.@test provenance["technology"] == technology
            Test.@test all(
                available -> !available,
                values(provenance["availability"]),
            )
            for field in (
                    "knobs",
                    "max_k",
                    "indel_params",
                    "substitution_error_rate",
            )
                Test.@test provenance[field] === nothing
            end
        end
        Test.@test isempty(result.assembly_stats["retained_cleanup_files"])
        Test.@test isempty(result.assembly_stats["retained_cleanup_roots"])
        Test.@test !result.gfa_compatible
        Test.@test all(call -> !isfile(call.corrected_fastq), correction_calls)
        Test.@test !isdir(dirname(captured_outdir[]))
    end

    Test.@testset "source content mutation fails before assembler" begin
        mktempdir() do temp_dir
            path_r1 = multi_input_write_fastq(
                joinpath(temp_dir, "R1.fastq"),
                MULTI_INPUT_R1,
            )
            path_r2 = multi_input_write_fastq(
                joinpath(temp_dir, "R2.fastq"),
                MULTI_INPUT_R2,
            )
            path_long = multi_input_write_fastq(
                joinpath(temp_dir, "long.fastq"),
                MULTI_INPUT_LONG,
            )
            path_calls = NamedTuple[]
            path_base_runner = multi_input_fake_correction_runner(path_calls)
            function path_mutating_correction_runner(
                    reads::Any,
                    config::Mycelia.Rhizomorph.AssemblyConfig,
            )::NamedTuple
                stage = path_base_runner(reads, config)
                if length(path_calls) == 1
                    multi_input_write_fastq(
                        path_r1,
                        FASTX.FASTQ.Record[
                            multi_input_fastq_record("pair_a/1", "TTTTTTTT"),
                            MULTI_INPUT_R1[2],
                        ],
                    )
                end
                return stage
            end
            path_assembler_calls = Ref(0)
            function path_assembler_runner(
                    inputs::Any,
                    outdir::AbstractString,
            )::NamedTuple
                path_assembler_calls[] += 1
                return multi_input_fake_assembler_result(outdir)
            end
            test_throws_message(
                ErrorException,
                "input short_r1 content changed during independent correction",
            ) do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (path_r1, path_r2),
                    path_long,
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                    :unicycler;
                    correction_runner = path_mutating_correction_runner,
                    assembler_runner = path_assembler_runner,
                )
            end
            Test.@test length(path_calls) == 3
            Test.@test path_assembler_calls[] == 0
            Test.@test all(call -> !isfile(call.corrected_fastq), path_calls)
        end

        in_memory_r1 = copy(MULTI_INPUT_R1)
        in_memory_r2 = copy(MULTI_INPUT_R2)
        in_memory_long = copy(MULTI_INPUT_LONG)
        in_memory_calls = NamedTuple[]
        in_memory_base_runner =
            multi_input_fake_correction_runner(in_memory_calls)
        function in_memory_mutating_correction_runner(
                reads::Any,
                config::Mycelia.Rhizomorph.AssemblyConfig,
        )::NamedTuple
            stage = in_memory_base_runner(reads, config)
            if length(in_memory_calls) == 1
                reads[1] = multi_input_fastq_record("pair_a/1", "GGGGGGGG")
            end
            return stage
        end
        in_memory_assembler_calls = Ref(0)
        function in_memory_assembler_runner(
                inputs::Any,
                outdir::AbstractString,
        )::NamedTuple
            in_memory_assembler_calls[] += 1
            return multi_input_fake_assembler_result(outdir)
        end
        test_throws_message(
            ErrorException,
            "input short_r1 content changed during independent correction",
        ) do
            Mycelia.Rhizomorph._assemble_paired_short_long(
                (in_memory_r1, in_memory_r2),
                in_memory_long,
                Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                :unicycler;
                correction_runner = in_memory_mutating_correction_runner,
                assembler_runner = in_memory_assembler_runner,
            )
        end
        Test.@test length(in_memory_calls) == 3
        Test.@test in_memory_assembler_calls[] == 0
        Test.@test all(
            call -> !isfile(call.corrected_fastq),
            in_memory_calls,
        )

        corrected_mutation_calls = NamedTuple[]
        corrected_assembler_calls = Ref(0)
        function corrected_mutating_assembler_runner(
                inputs::Mycelia.Rhizomorph._CorrectedPairedShortLong,
                outdir::AbstractString,
        )::NamedTuple
            corrected_assembler_calls[] += 1
            corrected_r1_records = multi_input_collect_fastq([
                inputs.short_r1.path,
            ])
            corrected_r1_records[1] = multi_input_fastq_record(
                String(FASTX.identifier(corrected_r1_records[1])),
                "CCCCCCCC",
            )
            multi_input_write_fastq(
                inputs.short_r1.path,
                corrected_r1_records,
            )
            return multi_input_fake_assembler_result(outdir)
        end
        test_throws_message(
            ErrorException,
            "corrected short_r1 FASTQ content changed during combined-input assembly",
        ) do
            Mycelia.Rhizomorph._assemble_paired_short_long(
                (MULTI_INPUT_R1, MULTI_INPUT_R2),
                MULTI_INPUT_LONG,
                Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                :unicycler;
                correction_runner = multi_input_fake_correction_runner(
                    corrected_mutation_calls,
                ),
                assembler_runner = corrected_mutating_assembler_runner,
            )
        end
        Test.@test length(corrected_mutation_calls) == 3
        Test.@test corrected_assembler_calls[] == 1
        Test.@test all(
            call -> !isfile(call.corrected_fastq),
            corrected_mutation_calls,
        )
    end

    Test.@testset "effective per-read-set correction provenance" begin
        correction_calls = NamedTuple[]
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
            correction_runner = multi_input_reported_correction_runner(
                correction_calls,
            ),
            assembler_runner = (inputs, outdir) ->
                multi_input_fake_assembler_result(outdir),
        )
        provenance = result.assembly_stats["correction_provenance"]
        Test.@test all(
            read_set -> read_set["status"] == "reported",
            values(provenance),
        )
        Test.@test all(
            read_set -> all(values(read_set["availability"])),
            values(provenance),
        )
        Test.@test provenance["short_r1"]["knobs"] == Dict(
            "strategy" => "scalable",
            "skip_solid" => true,
        )
        Test.@test provenance["short_r2"]["max_k"] == 17
        Test.@test provenance["short_r1"]["indel_params"] === nothing
        Test.@test provenance["long_reads"]["technology"] == "nanopore"
        Test.@test provenance["long_reads"]["knobs"]["skip_solid"] == false
        Test.@test provenance["long_reads"]["substitution_error_rate"] === nothing
        Test.@test provenance["long_reads"]["indel_params"] == Dict(
            "base_error_rate" => 0.1,
            "insertion_fraction" => 0.4,
            "deletion_fraction" => 0.4,
            "insertion_extend_probability" => 0.2,
            "deletion_extend_probability" => 0.2,
            "deletion_max_run" => 3,
            "max_insertion_run" => 4,
            "band_width" => 64,
        )
        serialized_provenance = Mycelia.JSON.json(provenance)
        parsed_provenance = Mycelia.JSON.parse(serialized_provenance)
        Test.@test parsed_provenance["long_reads"]["indel_params"]["band_width"] == 64
        Test.@test all(call -> !isfile(call.corrected_fastq), correction_calls)
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
            return (;
                assembly = "assembly.fasta",
                graph = "assembly.gfa",
                toolchain = multi_input_fake_unicycler_toolchain(),
            )
        end
        unicycler_config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
            assembler_options = (; kmers = "21,33"),
            threads = 6,
        )
        unicycler_result = Mycelia.Rhizomorph._run_multi_input_assembler(
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
        Test.@test unicycler_result.toolchain ==
                   multi_input_fake_unicycler_toolchain()

        test_throws_message(ErrorException, "did not report realized toolchain") do
            Mycelia.Rhizomorph._run_multi_input_assembler(
                Val(:unicycler),
                inputs,
                "/workflow/unicycler-unreported",
                unicycler_config;
                runner = (; kwargs...) -> (;
                    assembly = "assembly.fasta",
                    graph = "assembly.gfa",
                ),
            )
        end

        package_records = Any[
            Dict(
                "name" => "unicycler",
                "version" => "0.5.1",
                "build_string" => "pyhdfd78af_0",
                "channel" => "bioconda",
            ),
            Dict(
                "name" => "spades",
                "version" => "4.2.0",
                "build_string" => "h5ca1c30_1",
                "channel" => "bioconda",
            ),
        ]
        normalized_packages =
            Mycelia.Rhizomorph._normalize_unicycler_package_inventory(
                package_records,
            )
        Test.@test [package.name for package in normalized_packages] ==
                   ["spades", "unicycler"]
        original_inventory_digest =
            Mycelia.Rhizomorph._unicycler_package_inventory_sha256(
                normalized_packages,
            )
        reversed_inventory_digest =
            Mycelia.Rhizomorph._unicycler_package_inventory_sha256(
                reverse(package_records),
            )
        Test.@test reversed_inventory_digest == original_inventory_digest
        Test.@test Mycelia.Rhizomorph._unicycler_toolchain_provenance(
            inventory_reader = () -> reverse(normalized_packages),
        ) == Mycelia.Rhizomorph._unicycler_toolchain_provenance(
            inventory_reader = () -> normalized_packages,
        )
        changed_build_records = deepcopy(package_records)
        changed_build_records[1]["build_string"] = "pyhdfd78af_1"
        changed_inventory_digest =
            Mycelia.Rhizomorph._unicycler_package_inventory_sha256(
                Mycelia.Rhizomorph._normalize_unicycler_package_inventory(
                    changed_build_records,
                ),
            )
        Test.@test original_inventory_digest != changed_inventory_digest
        Test.@test occursin(r"^[0-9a-f]{64}$", original_inventory_digest)

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
                    toolchain = multi_input_fake_autocycler_toolchain(),
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
        Test.@test result.assembly_stats["toolchain"] ==
                   multi_input_fake_autocycler_toolchain()
        Test.@test isempty(result.assembly_stats["retained_intermediates"])

        noncareful_calls = NamedTuple[]
        noncareful_result = Mycelia.Rhizomorph._assemble_paired_short_long(
            (MULTI_INPUT_R1, MULTI_INPUT_R2),
            MULTI_INPUT_LONG,
            Mycelia.Rhizomorph.AutocyclerPolishConfig(polypolish_careful = false),
            :autocycler_polished;
            correction_runner = multi_input_fake_correction_runner(noncareful_calls),
            assembler_runner = (inputs, outdir) ->
                multi_input_fake_assembler_result(
                    outdir;
                    include_polishing_artifacts = true,
                ),
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

    Test.@testset "retained cleanup artifacts are recorded" begin
        correction_calls = NamedTuple[]
        correction_runner = multi_input_fake_correction_runner(correction_calls)
        workflow_root = Ref("")
        reported_intermediate = Ref("")
        function assembler_runner(
                inputs::Any,
                outdir::AbstractString,
        )::NamedTuple
            workflow_root[] = dirname(String(outdir))
            assembler_result = multi_input_fake_assembler_result(outdir)
            reported_intermediate[] = joinpath(outdir, "retained.sam")
            write(reported_intermediate[], "retained\n")
            return merge(
                assembler_result,
                (; intermediates = [reported_intermediate[]]),
            )
        end
        function selective_corrected_remover(path::AbstractString)::Nothing
            if path == correction_calls[1].corrected_fastq
                error("synthetic corrected cleanup failure")
            end
            rm(path; force = true)
            return nothing
        end
        function failing_root_remover(path::AbstractString)::Nothing
            error("synthetic root cleanup failure: $(path)")
        end
        result = Mycelia.Rhizomorph._assemble_paired_short_long(
            (MULTI_INPUT_R1, MULTI_INPUT_R2),
            MULTI_INPUT_LONG,
            Mycelia.Rhizomorph.UnicyclerHybridConfig(),
            :unicycler;
            correction_runner,
            assembler_runner,
            corrected_fastq_remover = selective_corrected_remover,
            workflow_root_remover = failing_root_remover,
        )
        retained_file = correction_calls[1].corrected_fastq
        Test.@test result.assembly_stats["retained_cleanup_files"] ==
                   [retained_file]
        Test.@test result.assembly_stats["retained_cleanup_roots"] ==
                   [workflow_root[]]
        Test.@test result.assembly_stats["retained_intermediates"] ==
                   [reported_intermediate[]]
        Test.@test isfile(retained_file)
        Test.@test isfile(reported_intermediate[])
        Test.@test all(
            call -> !isfile(call.corrected_fastq),
            correction_calls[2:3],
        )
        Test.@test isdir(workflow_root[])
        rm(retained_file; force = true)
        rm(workflow_root[]; recursive = true, force = true)

        cleaned_correction_calls = NamedTuple[]
        cleaned_intermediate = Ref("")
        function cleaned_assembler_runner(
                inputs::Any,
                outdir::AbstractString,
        )::NamedTuple
            assembler_result = multi_input_fake_assembler_result(outdir)
            cleaned_intermediate[] = joinpath(outdir, "cleaned.sam")
            write(cleaned_intermediate[], "cleaned\n")
            return merge(
                assembler_result,
                (; intermediates = [cleaned_intermediate[]]),
            )
        end
        cleaned_result = Mycelia.Rhizomorph._assemble_paired_short_long(
            (MULTI_INPUT_R1, MULTI_INPUT_R2),
            MULTI_INPUT_LONG,
            Mycelia.Rhizomorph.UnicyclerHybridConfig(),
            :unicycler;
            correction_runner = multi_input_fake_correction_runner(
                cleaned_correction_calls,
            ),
            assembler_runner = cleaned_assembler_runner,
        )
        Test.@test length(cleaned_correction_calls) == 3
        Test.@test isempty(
            cleaned_result.assembly_stats["retained_intermediates"],
        )
        Test.@test !ispath(cleaned_intermediate[])
        Test.@test cleaned_result.assembly_stats["input_technologies"] == Dict(
            "short_r1" => "illumina",
            "short_r2" => "illumina",
            "long_reads" => "nanopore",
        )
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

    Test.@testset "assembler artifacts stay inside their reserved outdir" begin
        mktempdir() do temp_dir
            external_fasta = multi_input_write_fasta(
                joinpath(temp_dir, "external.fasta"),
            )
            external_graph = joinpath(temp_dir, "external.gfa")
            write(
                external_graph,
                "H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\n",
            )
            external_report = joinpath(temp_dir, "external.report")
            write(external_report, "corrected_bases\t0\n")
            external_paths = Dict(
                :assembly => external_fasta,
                :graph => external_graph,
                :autocycler_assembly => external_fasta,
                :polypolish_assembly => external_fasta,
                :pypolca_report => external_report,
            )

            for (case_index, field) in enumerate(keys(external_paths))
                outdir = joinpath(temp_dir, "external-field-$(case_index)")
                result = multi_input_fake_assembler_result(
                    outdir;
                    include_polishing_artifacts = true,
                )
                escaped_result = merge(
                    result,
                    NamedTuple{(field,)}((external_paths[field],)),
                )
                test_throws_message(
                    ErrorException,
                    "outside its exact reserved output directory",
                ) do
                    Mycelia.Rhizomorph._validate_multi_input_assembler_result(
                        escaped_result,
                        :autocycler_polished,
                        outdir,
                    )
                end
            end

            intermediate_outdir = joinpath(temp_dir, "external-intermediate")
            intermediate_result = merge(
                multi_input_fake_assembler_result(
                    intermediate_outdir;
                    include_polishing_artifacts = true,
                ),
                (; intermediates = [external_fasta]),
            )
            test_throws_message(
                ErrorException,
                "intermediate 1 is outside its exact reserved output directory",
            ) do
                Mycelia.Rhizomorph._validate_multi_input_assembler_result(
                    intermediate_result,
                    :autocycler_polished,
                    intermediate_outdir,
                )
            end

            symlink_outdir = joinpath(temp_dir, "symlink-artifact")
            symlink_result = multi_input_fake_assembler_result(symlink_outdir)
            symlinked_assembly = joinpath(symlink_outdir, "escaped.fasta")
            symlink(external_fasta, symlinked_assembly)
            test_throws_message(ErrorException, "regular, non-symlink file") do
                Mycelia.Rhizomorph._validate_multi_input_assembler_result(
                    merge(symlink_result, (; assembly = symlinked_assembly)),
                    :unicycler,
                    symlink_outdir,
                )
            end

            nested_symlink_outdir = joinpath(temp_dir, "nested-symlink")
            nested_symlink_result = multi_input_fake_assembler_result(
                nested_symlink_outdir,
            )
            nested_real_dir = joinpath(nested_symlink_outdir, "real")
            nested_alias_dir = joinpath(nested_symlink_outdir, "alias")
            mkpath(nested_real_dir)
            nested_assembly = multi_input_write_fasta(
                joinpath(nested_real_dir, "assembly.fasta"),
            )
            symlink(nested_real_dir, nested_alias_dir)
            test_throws_message(
                ErrorException,
                "resolves through a symbolic-link path component",
            ) do
                Mycelia.Rhizomorph._validate_multi_input_assembler_result(
                    merge(
                        nested_symlink_result,
                        (;
                            assembly = joinpath(
                                nested_alias_dir,
                                basename(nested_assembly),
                            ),
                        ),
                    ),
                    :unicycler,
                    nested_symlink_outdir,
                )
            end

            correction_calls = NamedTuple[]
            test_throws_message(
                ErrorException,
                "assembly is outside its exact reserved output directory",
            ) do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (MULTI_INPUT_R1, MULTI_INPUT_R2),
                    MULTI_INPUT_LONG,
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(),
                    :unicycler;
                    correction_runner = multi_input_fake_correction_runner(
                        correction_calls,
                    ),
                    assembler_runner = (inputs, outdir) -> merge(
                        multi_input_fake_assembler_result(outdir),
                        (; assembly = external_fasta),
                    ),
                )
            end
            workflow_output = joinpath(temp_dir, "retargeted-assembler-child")
            external_assembler_dir = joinpath(temp_dir, "external-assembler")
            mkpath(external_assembler_dir)
            retargeting_correction_calls = NamedTuple[]
            base_correction_runner = multi_input_fake_correction_runner(
                retargeting_correction_calls,
            )
            function assembler_retargeting_correction_runner(
                    reads::Any,
                    config::Mycelia.Rhizomorph.AssemblyConfig,
            )::NamedTuple
                stage = base_correction_runner(reads, config)
                if length(retargeting_correction_calls) == 3
                    symlink(
                        external_assembler_dir,
                        joinpath(workflow_output, "assembler_unicycler"),
                    )
                end
                return stage
            end
            assembler_called = Ref(false)
            test_throws_message(
                ArgumentError,
                "contains a symbolic-link path component",
            ) do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (MULTI_INPUT_R1, MULTI_INPUT_R2),
                    MULTI_INPUT_LONG,
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(
                        output_dir = workflow_output,
                    ),
                    :unicycler;
                    correction_runner =
                        assembler_retargeting_correction_runner,
                    assembler_runner = (inputs, outdir) -> begin
                        assembler_called[] = true
                        return multi_input_fake_assembler_result(outdir)
                    end,
                )
            end
            Test.@test !assembler_called[]
            Test.@test isempty(readdir(external_assembler_dir))

            postcheck_workflow = joinpath(
                temp_dir,
                "assembler-postcheck-workflow",
            )
            postcheck_assembler_dir = joinpath(
                temp_dir,
                "assembler-postcheck-external",
            )
            mkpath(postcheck_assembler_dir)
            postcheck_calls = NamedTuple[]
            test_throws_message(
                ArgumentError,
                "contains a symbolic-link path component",
            ) do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (MULTI_INPUT_R1, MULTI_INPUT_R2),
                    MULTI_INPUT_LONG,
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(
                        output_dir = postcheck_workflow,
                    ),
                    :unicycler;
                    correction_runner = multi_input_fake_correction_runner(
                        postcheck_calls,
                    ),
                    assembler_runner = (inputs, outdir) -> begin
                        symlink(postcheck_assembler_dir, outdir)
                        return multi_input_fake_assembler_result(outdir)
                    end,
                )
            end
            Test.@test isfile(joinpath(postcheck_assembler_dir, "assembly.fasta"))
            Test.@test isfile(joinpath(postcheck_assembler_dir, "assembly.gfa"))

            replaced_workflow = joinpath(temp_dir, "replaced-workflow")
            replacement_marker = joinpath(
                replaced_workflow,
                "replacement-marker.txt",
            )
            replacement_calls = NamedTuple[]
            test_throws_message(
                ErrorException,
                "multi-input workflow root changed filesystem identity",
            ) do
                Mycelia.Rhizomorph._assemble_paired_short_long(
                    (MULTI_INPUT_R1, MULTI_INPUT_R2),
                    MULTI_INPUT_LONG,
                    Mycelia.Rhizomorph.UnicyclerHybridConfig(
                        output_dir = replaced_workflow,
                    ),
                    :unicycler;
                    correction_runner = multi_input_fake_correction_runner(
                        replacement_calls,
                    ),
                    assembler_runner = (inputs, outdir) -> begin
                        staged_replacement = joinpath(
                            temp_dir,
                            "staged-replacement-workflow",
                        )
                        mkpath(staged_replacement)
                        write(
                            joinpath(
                                staged_replacement,
                                basename(replacement_marker),
                            ),
                            "preserve replacement\n",
                        )
                        rm(replaced_workflow; recursive = true)
                        mv(staged_replacement, replaced_workflow)
                        return multi_input_fake_assembler_result(outdir)
                    end,
                )
            end
            Test.@test read(replacement_marker, String) ==
                       "preserve replacement\n"
            Test.@test isfile(external_fasta)
            Test.@test isfile(external_graph)
            Test.@test isfile(external_report)
            Test.@test all(call -> !isfile(call.corrected_fastq), correction_calls)
        end
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
        gzip_graph = joinpath(artifact_dir, "assembly.gfa")
        write(
            gzip_graph,
            "H\tVN:Z:1.0\n" *
            "S\talpha\tACGT\n" *
            "S\tbeta\tTTAA\n" *
            "L\talpha\t+\tbeta\t-\t1M\n" *
            "P\tprimary\talpha+,beta-\t1M\n",
        )

        wrapped = Mycelia.Rhizomorph._wrap_multi_input_assembly(
            (;
                contigs = gzip_fasta,
                graph = gzip_graph,
                toolchain = multi_input_fake_unicycler_toolchain(),
            ),
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
            correction_provenance = Dict{String, Any}(),
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
            retained_cleanup_files = String[],
            retained_cleanup_roots = String[],
        )
        Test.@test wrapped.contig_names == ["alpha", "beta"]
        Test.@test wrapped.contigs == ["ACGT", "TTAA"]

        comma_identifier_graph = joinpath(
            artifact_dir,
            "comma-identifier.gfa",
        )
        write(
            comma_identifier_graph,
            "H\tVN:Z:1.0\n" *
            "S\tsegment,name\tACGT\n" *
            "P\tprimary\tsegment,name+\t*\n",
        )
        Test.@test Mycelia.Rhizomorph._require_valid_multi_input_gfa(
            comma_identifier_graph,
            "comma-name graph",
        ) == comma_identifier_graph

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
                correction_provenance = Dict{String, Any}(),
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
                retained_cleanup_files = String[],
                retained_cleanup_roots = String[],
            )
        end
    end

    Test.@testset "semantic assembler artifacts fail loud" begin
        mktempdir() do temp_dir
            function wrap_semantic_result(
                    result::NamedTuple,
                    workflow::Symbol,
            )::Mycelia.Rhizomorph.AssemblyResult
                assembler = workflow == :autocycler_polished ?
                            "autocycler" : "unicycler"
                polishers = workflow == :autocycler_polished ?
                            ["polypolish-careful", "pypolca-careful"] :
                            String[]
                return Mycelia.Rhizomorph._wrap_multi_input_assembly(
                    result,
                    workflow;
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
                    correction_provenance = Dict{String, Any}(),
                    output_dir = nothing,
                    polishers,
                    workflow_settings = Dict{String, Any}(
                        "workflow" => String(workflow),
                        "assembler" => assembler,
                    ),
                    input_technologies = Dict(
                        "short_r1" => "illumina",
                        "short_r2" => "illumina",
                        "long_reads" => "nanopore",
                    ),
                    retained_cleanup_files = String[],
                    retained_cleanup_roots = String[],
                )
            end

            function semantic_result_shape(
                    workflow::Symbol,
                    outdir::AbstractString;
                    assembly_records::Vector{Pair{String, String}} =
                        ["contig_1" => "ACGT"],
                    graph_contents::AbstractString =
                        "H\tVN:Z:1.0\nS\tcontig_1\tACGT\n",
                    autocycler_records::Vector{Pair{String, String}} =
                        ["contig_1" => "ACGT"],
                    polypolish_records::Vector{Pair{String, String}} =
                        ["contig_1" => "ACGT"],
            )::NamedTuple
                assembly = multi_input_write_fasta(
                    joinpath(outdir, "assembly.fasta");
                    records = assembly_records,
                )
                graph = joinpath(outdir, "assembly.gfa")
                write(graph, graph_contents)
                if workflow == :unicycler
                    return (;
                        assembly,
                        graph,
                        toolchain = multi_input_fake_unicycler_toolchain(),
                    )
                end
                autocycler_assembly = multi_input_write_fasta(
                    joinpath(outdir, "autocycler_assembly.fasta");
                    records = autocycler_records,
                )
                polypolish_assembly = multi_input_write_fasta(
                    joinpath(outdir, "polypolish_assembly.fasta");
                    records = polypolish_records,
                )
                pypolca_report = joinpath(outdir, "pypolca.report")
                write(pypolca_report, "corrected_bases\t0\n")
                return (;
                    assembly,
                    graph,
                    autocycler_assembly,
                    polypolish_assembly,
                    pypolca_report,
                    toolchain = multi_input_fake_autocycler_toolchain(),
                )
            end

            semantic_cases = (
                (
                    message = "empty FASTA identifier",
                    assembly_records = ["" => "ACGT"],
                    graph_contents = "H\tVN:Z:1.0\nS\tcontig_1\tACGT\n",
                ),
                (
                    message = "duplicate FASTA identifier",
                    assembly_records = [
                        "duplicate" => "ACGT",
                        "duplicate" => "TGCA",
                    ],
                    graph_contents = "H\tVN:Z:1.0\nS\tduplicate\tACGT\n",
                ),
                (
                    message = "empty FASTA sequence",
                    assembly_records = ["empty" => ""],
                    graph_contents = "H\tVN:Z:1.0\nS\tempty\tACGT\n",
                ),
                (
                    message = "invalid DNA at FASTA record",
                    assembly_records = ["invalid" => "ACGTZ"],
                    graph_contents = "H\tVN:Z:1.0\nS\tinvalid\tACGT\n",
                ),
                (
                    message = "invalid DNA at FASTA record",
                    assembly_records = ["gapped" => "ACGT-ACGT"],
                    graph_contents = "H\tVN:Z:1.0\nS\tgapped\tACGT\n",
                ),
                (
                    message = "unknown GFA record type",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents =
                        "not GFA\nH\tVN:Z:1.0\nS\tcontig_1\tACGT\n",
                ),
                (
                    message = "malformed GFA segment",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents = "H\tVN:Z:1.0\nS\tbroken\n",
                ),
                (
                    message = "malformed GFA link",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents =
                        "S\tcontig_1\tACGT\nS\tcontig_2\tTGCA\n" *
                        "L\tcontig_1\t+\tcontig_2\t-\n",
                ),
                (
                    message = "dangling GFA link segment reference",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents =
                        "S\tcontig_1\tACGT\n" *
                        "L\tcontig_1\t+\tmissing\t-\t1M\n",
                ),
                (
                    message = "dangling GFA path segment reference",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents =
                        "S\tcontig_1\tACGT\n" *
                        "P\tprimary\tcontig_1+,missing-\t1M\n",
                ),
                (
                    message = "no sequence-bearing GFA segments",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents = "H\tVN:Z:1.0\n",
                ),
                (
                    message = "duplicate GFA segment/path name",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents =
                        "S\tduplicate\tACGT\nS\tduplicate\tTGCA\n",
                ),
                (
                    message = "invalid DNA for GFA segment",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents = "S\tinvalid\tACGTZ\n",
                ),
                (
                    message = "invalid DNA for GFA segment",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents = "S\tgapped\tACGT-ACGT\n",
                ),
                (
                    message = "duplicate GFA segment/path name",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents =
                        "S\tshared_name\tACGT\n" *
                        "P\tshared_name\tshared_name+\t*\n",
                ),
                (
                    message = "invalid GFA segment optional tag",
                    assembly_records = ["contig_1" => "ACGT"],
                    graph_contents = "S\tcontig_1\tACGT\tbad-tag\n",
                ),
            )
            for workflow in (:unicycler, :autocycler_polished)
                for (case_index, semantic_case) in enumerate(semantic_cases)
                    outdir = joinpath(
                        temp_dir,
                        "$(workflow)-semantic-case-$(case_index)",
                    )
                    result = semantic_result_shape(
                        workflow,
                        outdir;
                        assembly_records = semantic_case.assembly_records,
                        graph_contents = semantic_case.graph_contents,
                    )
                    test_throws_message(ErrorException, semantic_case.message) do
                        wrap_semantic_result(result, workflow)
                    end
                end
            end

            unicycler_without_graph = semantic_result_shape(
                :unicycler,
                joinpath(temp_dir, "unicycler-without-graph"),
            )
            test_throws_message(ErrorException, "field :graph") do
                wrap_semantic_result(
                    multi_input_without_field(
                        unicycler_without_graph,
                        :graph,
                    ),
                    :unicycler,
                )
            end

            for field in (
                    :graph,
                    :autocycler_assembly,
                    :polypolish_assembly,
                    :pypolca_report,
            )
                complete_result = semantic_result_shape(
                    :autocycler_polished,
                    joinpath(temp_dir, "autocycler-without-$(field)"),
                )
                test_throws_message(ErrorException, "field :$(field)") do
                    wrap_semantic_result(
                        multi_input_without_field(complete_result, field),
                        :autocycler_polished,
                    )
                end
            end

            complete_autocycler_result = semantic_result_shape(
                :autocycler_polished,
                joinpath(temp_dir, "autocycler-toolchain-contract"),
            )
            test_throws_message(
                ErrorException,
                "did not report realized toolchain provenance",
            ) do
                wrap_semantic_result(
                    multi_input_without_field(
                        complete_autocycler_result,
                        :toolchain,
                    ),
                    :autocycler_polished,
                )
            end
            malformed_toolchain = copy(complete_autocycler_result.toolchain)
            malformed_toolchain["toolchain_schema"] = "unsupported-schema"
            test_throws_message(
                ErrorException,
                "incompatible toolchain_schema",
            ) do
                wrap_semantic_result(
                    merge(
                        complete_autocycler_result,
                        (; toolchain = malformed_toolchain),
                    ),
                    :autocycler_polished,
                )
            end
            build_drift_toolchain = deepcopy(
                complete_autocycler_result.toolchain,
            )
            build_drift_toolchain["packages"][1]["build"] = "drifted-build"
            test_throws_message(
                ErrorException,
                "inventory digest does not match",
            ) do
                wrap_semantic_result(
                    merge(
                        complete_autocycler_result,
                        (; toolchain = build_drift_toolchain),
                    ),
                    :autocycler_polished,
                )
            end
            alternate_build_toolchain = multi_input_fake_autocycler_toolchain(
                build_suffix = "alternate",
            )
            alternate_build_result = wrap_semantic_result(
                merge(
                    complete_autocycler_result,
                    (; toolchain = alternate_build_toolchain),
                ),
                :autocycler_polished,
            )
            Test.@test alternate_build_result.assembly_stats["toolchain"] ==
                       alternate_build_toolchain
            Test.@test alternate_build_toolchain[
                "package_inventory_sha256"
            ] != complete_autocycler_result.toolchain[
                "package_inventory_sha256"
            ]

            empty_report_result = semantic_result_shape(
                :autocycler_polished,
                joinpath(temp_dir, "empty-pypolca-report"),
            )
            write(empty_report_result.pypolca_report, "")
            test_throws_message(ErrorException, "no non-empty Pypolca report") do
                wrap_semantic_result(
                    empty_report_result,
                    :autocycler_polished,
                )
            end

            invalid_autocycler = semantic_result_shape(
                :autocycler_polished,
                joinpath(temp_dir, "invalid-autocycler-intermediate");
                autocycler_records = ["empty" => ""],
            )
            test_throws_message(ErrorException, "empty FASTA sequence") do
                wrap_semantic_result(invalid_autocycler, :autocycler_polished)
            end
            invalid_polypolish = semantic_result_shape(
                :autocycler_polished,
                joinpath(temp_dir, "invalid-polypolish-intermediate");
                polypolish_records = ["invalid" => "ACGTZ"],
            )
            test_throws_message(ErrorException, "invalid DNA at FASTA record") do
                wrap_semantic_result(invalid_polypolish, :autocycler_polished)
            end

            matched_records = [
                "contig_1" => "ACGT",
                "contig_2" => "TGCA",
            ]
            renamed_autocycler = semantic_result_shape(
                :autocycler_polished,
                joinpath(temp_dir, "renamed-autocycler-contig");
                assembly_records = matched_records,
                graph_contents =
                    "S\tcontig_1\tACGT\nS\tcontig_2\tTGCA\n",
                autocycler_records = [
                    "renamed_contig" => "ACGT",
                    "contig_2" => "TGCA",
                ],
                polypolish_records = matched_records,
            )
            test_throws_message(
                ErrorException,
                "Autocycler consensus FASTA contig identifier set does not " *
                "match final assembly",
            ) do
                wrap_semantic_result(
                    renamed_autocycler,
                    :autocycler_polished,
                )
            end

            dropped_polypolish = semantic_result_shape(
                :autocycler_polished,
                joinpath(temp_dir, "dropped-polypolish-contig");
                assembly_records = matched_records,
                graph_contents =
                    "S\tcontig_1\tACGT\nS\tcontig_2\tTGCA\n",
                autocycler_records = matched_records,
                polypolish_records = ["contig_1" => "ACGT"],
            )
            test_throws_message(
                ErrorException,
                "Polypolish assembly FASTA contig identifier set does not " *
                "match final assembly",
            ) do
                wrap_semantic_result(
                    dropped_polypolish,
                    :autocycler_polished,
                )
            end
        end
    end

    Test.@testset "public hybrid aliases are documented" begin
        unicycler_binding = Base.Docs.Binding(
            Mycelia.Rhizomorph,
            :assemble_unicycler_hybrid,
        )
        autocycler_binding = Base.Docs.Binding(
            Mycelia.Rhizomorph,
            :assemble_autocycler_polished,
        )
        Test.@test Base.Docs.doc(unicycler_binding) !== nothing
        Test.@test Base.Docs.doc(autocycler_binding) !== nothing
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
        shared_pair_record = multi_input_fastq_record("shared_pair")
        shared_pair_r1 = [shared_pair_record]
        shared_pair_r2 = copy(shared_pair_r1)
        test_throws_message(
            ArgumentError,
            "share an in-memory FASTQ record object",
        ) do
            Mycelia.Rhizomorph.assemble_hybrid(
                (shared_pair_r1, shared_pair_r2),
                MULTI_INPUT_LONG;
                config = Mycelia.Rhizomorph.UnicyclerHybridConfig(),
            )
        end
        shared_long_r1 = [multi_input_fastq_record("shared_long/1")]
        shared_long_r2 = [multi_input_fastq_record("shared_long/2")]
        test_throws_message(
            ArgumentError,
            "share an in-memory FASTQ record object",
        ) do
            Mycelia.Rhizomorph.assemble_hybrid(
                (shared_long_r1, shared_long_r2),
                copy(shared_long_r1);
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
            distinct_r2_alias = joinpath(temp_dir, "distinct-r2-alias.fastq")
            symlink(distinct_r2, distinct_r2_alias)
            long_hardlink = joinpath(temp_dir, "long-hardlink.fastq")
            Base.hardlink(long_path, long_hardlink)
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
                (
                    message = "input short_r1 sources must be distinct",
                    short_r1 = (paired_path, paired_path),
                    short_r2 = (distinct_r2,),
                    long_reads = (long_path,),
                    output_name = "duplicate-within-r1",
                ),
                (
                    message = "input short_r2 sources must be distinct",
                    short_r1 = (paired_path,),
                    short_r2 = (distinct_r2, distinct_r2_alias),
                    long_reads = (long_path,),
                    output_name = "symlink-alias-within-r2",
                ),
                (
                    message = "input long_reads sources must be distinct",
                    short_r1 = (paired_path,),
                    short_r2 = (distinct_r2,),
                    long_reads = (long_path, long_hardlink),
                    output_name = "hardlink-alias-within-long",
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
