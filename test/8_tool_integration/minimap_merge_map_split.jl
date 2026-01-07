# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/minimap_merge_map_split.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/minimap_merge_map_split.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import FASTX
import UUIDs
import XAM

Test.@testset "Minimap merge/map/split helpers" begin
    Test.@testset "prefix_fastq_reads adds unique tag and mapping" begin
        temp_fastq = tempname() * ".fq"
        write(temp_fastq, "@r1\nACGT\n+\nIIII\n@r1\nTGCA\n+\nHHHH\n")
        out_fastq = tempname() * ".fq"
        map_tsv = tempname() * ".tsv"

        res = Mycelia.prefix_fastq_reads(
            temp_fastq;
            sample_tag="sampleA",
            out_fastq=out_fastq,
            mapping_out=map_tsv,
            mapping_format=:tsv,
            id_delimiter="::"
        )

        # identifiers carry the tag and delimiter
        reader = Mycelia.open_fastx(res.out_fastq)
        ids = String[]
        for rec in reader
            push!(ids, String(FASTX.identifier(rec)))
        end
        close(reader)
        Test.@test all(startswith(id, "sampleA::") for id in ids)

        tsv_lines = readlines(map_tsv)
        Test.@test length(tsv_lines) == 3 # header + 2 reads
    end

    Test.@testset "uuid_fastq_reads rewrites ids and mapping" begin
        temp_fastq = tempname() * ".fq"
        write(temp_fastq, "@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nHHHH\n")
        out_fastq = tempname() * ".fq"
        map_tsv = tempname() * ".tsv"

        res = Mycelia.uuid_fastq_reads(
            temp_fastq;
            out_fastq=out_fastq,
            mapping_out=map_tsv,
            mapping_format=:tsv
        )

        reader = Mycelia.open_fastx(res.out_fastq)
        ids = String[]
        for rec in reader
            push!(ids, String(FASTX.identifier(rec)))
        end
        close(reader)
        Test.@test all(id -> (try UUIDs.UUID(id); true catch; false end), ids)

        tsv_lines = readlines(map_tsv)
        Test.@test tsv_lines[1] == "source_fastq\toriginal_read_id\tuuid"
        Test.@test length(tsv_lines) == 3 # header + 2 reads
    end

    Test.@testset "Command construction without execution" begin
        ref = tempname() * ".fa"
        write(ref, ">ref\nACGT\n")
        fq1 = tempname() * ".fq"
        fq2 = tempname() * ".fq"
        write(fq1, "@a\nAAAA\n+\nIIII\n")
        write(fq2, "@b\nCCCC\n+\nIIII\n")

        pair1 = tempname() * ".fq"
        pair2 = tempname() * ".fq"
        write(pair1, "@read/1\nGGGG\n+\nIIII\n")
        write(pair2, "@read/2\nGGGG\n+\nIIII\n")

        res = Mycelia.minimap_merge_map_and_split(
            reference_fasta=ref,
            mapping_type="sr",
            single_end_fastqs=[fq1, fq2],
            paired_end_fastqs=[(pair1, pair2)],
            run_mapping=false,
            run_splitting=false,
            write_read_map=false,
            gzip_prefixed_fastqs=false,
            minimap_extra_args=["-N", "5"],
            as_string=true
        )

        Test.@test occursin("minimap2", res.minimap_cmd)
        Test.@test occursin("-x sr", res.minimap_cmd)
        Test.@test length(res.prefixed_fastqs) == 4
        Test.@test length(res.sample_outputs) == 3
        Test.@test occursin(basename(ref), res.merged_bam)
    end

    Test.@testset "Accepts minimap_index without reference_fasta" begin
        idx = tempname() * ".mmi"
        write(idx, "dummy index placeholder")
        fq = tempname() * ".fq"
        write(fq, "@a\nAAAA\n+\nIIII\n")

        res = Mycelia.minimap_merge_map_and_split(
            reference_fasta=nothing,
            minimap_index=idx,
            mapping_type="sr",
            single_end_fastqs=[fq],
            run_mapping=false,
            run_splitting=false,
            gzip_prefixed_fastqs=false,
            as_string=true
        )

        Test.@test occursin(basename(idx), res.merged_bam)
        Test.@test occursin(idx, res.minimap_cmd)
    end

    Test.@testset "Paired identifiers stay linked after prefixing" begin
        r1 = tempname() * ".fq"
        r2 = tempname() * ".fq"
        write(r1, "@read1/1\nAAAA\n+\nIIII\n")
        write(r2, "@read1/2\nAAAA\n+\nIIII\n")
        out1 = tempname() * ".fq.gz"
        out2 = tempname() * ".fq.gz"

        Mycelia.prefix_fastq_reads(r1; sample_tag="sampleX", out_fastq=out1, id_delimiter="::")
        Mycelia.prefix_fastq_reads(r2; sample_tag="sampleX", out_fastq=out2, id_delimiter="::")

        rec1_reader = Mycelia.open_fastx(out1)
        rec2_reader = Mycelia.open_fastx(out2)
        rec1 = first(rec1_reader)
        rec2 = first(rec2_reader)

        base1 = replace(String(FASTX.identifier(rec1)), "/1" => "")
        base2 = replace(String(FASTX.identifier(rec2)), "/2" => "")
        Test.@test base1 == base2
        close(rec1_reader)
        close(rec2_reader)
    end

    Test.@testset "Mapping generates sorted BAMs with index-aware names" begin
        run_external = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true" ||
            lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"
        if !run_external
            @info "Skipping minimap2 integration test; set MYCELIA_RUN_EXTERNAL=true to run"
            Test.@test_skip "Requires minimap2 + samtools (external tools)"
        else
            ref_seq = repeat("ACGT", 50)
            ref = tempname() * ".fna"
            write(ref, ">ref\n", ref_seq, "\n")

            reads = FASTX.FASTQ.Record[]
            push!(reads, FASTX.FASTQ.Record("read1", ref_seq[1:80], repeat("I", 80)))
            push!(reads, FASTX.FASTQ.Record("read2", ref_seq[21:100], repeat("I", 80)))

            fq = tempname() * ".bam.fastplong.fq.gz"
            Mycelia.write_fastq(records=reads, filename=fq)

            res = Mycelia.minimap_merge_map_and_split(
                reference_fasta=ref,
                mapping_type="sr",
                single_end_fastqs=[fq],
                build_index=true,
                run_mapping=true,
                run_splitting=true,
                gzip_prefixed_fastqs=false,
                write_read_map=false
            )

            Test.@test length(res.sample_outputs) == 1
            Test.@test res.index_file !== nothing
            sample = first(res.sample_outputs)
            expected_fastq_label = Mycelia.build_output_label([fq])
            expected_bam = joinpath(dirname(fq), "$(expected_fastq_label).$(basename(res.index_file)).sorted.bam")
            Test.@test sample.output_bam == expected_bam
            Test.@test isfile(sample.output_bam)
            Test.@test Mycelia.is_bam_coordinate_sorted(sample.output_bam)
        end
    end

    Test.@testset "UUID mapping splits are recoverable (per-sample and joint)" begin
        run_external = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true" ||
            lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"
        if !run_external
            @info "Skipping minimap2 integration test; set MYCELIA_RUN_EXTERNAL=true to run"
            Test.@test_skip "Requires minimap2 + samtools (external tools)"
        else
            ref_seq = repeat("ACGT", 60)
            ref = tempname() * ".fna"
            write(ref, ">ref\n", ref_seq, "\n")

            reads_a = FASTX.FASTQ.Record[]
            push!(reads_a, FASTX.FASTQ.Record("readA1", ref_seq[1:80], repeat("I", 80)))
            push!(reads_a, FASTX.FASTQ.Record("readA2", ref_seq[41:120], repeat("I", 80)))
            fq_a = tempname() * ".fastq"
            Mycelia.write_fastq(records=reads_a, filename=fq_a)

            reads_b = FASTX.FASTQ.Record[]
            push!(reads_b, FASTX.FASTQ.Record("readB1", ref_seq[81:160], repeat("I", 80)))
            fq_b = tempname() * ".fastq"
            Mycelia.write_fastq(records=reads_b, filename=fq_b)

            index_res = Mycelia.minimap_index(
                fasta=ref,
                mapping_type="sr",
                mem_gb=1,
                threads=1,
                as_string=false
            )
            run(index_res.cmd)
            idx = index_res.outfile

            map_res = Mycelia.minimap_map_with_index(
                fasta=nothing,
                mapping_type="sr",
                fastq=fq_a,
                index_file=idx,
                sorted=true
            )
            run(map_res.cmd)
            Test.@test isfile(map_res.outfile)
            Test.@test Mycelia.is_bam_coordinate_sorted(map_res.outfile)

            per_res = Mycelia.minimap_merge_map_and_split(
                reference_fasta=nothing,
                minimap_index=idx,
                mapping_type="sr",
                single_end_fastqs=[fq_a, fq_b],
                read_id_strategy=:uuid,
                fastq_mode=:per_sample,
                read_map_format=:tsv,
                gzip_read_map_tsv=false,
                gzip_prefixed_fastqs=false,
                run_mapping=true,
                run_splitting=true
            )

            for sample in per_res.sample_outputs
                Test.@test isfile(sample.output_bam)
                Test.@test Mycelia.is_bam_coordinate_sorted(sample.output_bam)
                uuid_set = Set{String}()
                for map_path in sample.mapping_files
                    for (i, line) in enumerate(eachline(map_path))
                        i == 1 && continue
                        isempty(line) && continue
                        tab_idx = findlast('\t', line)
                        uuid = tab_idx === nothing ? line : line[tab_idx + 1:end]
                        push!(uuid_set, uuid)
                    end
                end
                reader = Mycelia.open_xam(sample.output_bam; parser=:samtools)
                for record in reader
                    Test.@test XAM.SAM.tempname(record) in uuid_set
                end
                close(reader)
            end

            joint_res = Mycelia.minimap_merge_map_and_split(
                reference_fasta=nothing,
                minimap_index=idx,
                mapping_type="sr",
                single_end_fastqs=[fq_a, fq_b],
                read_id_strategy=:uuid,
                fastq_mode=:joint,
                read_map_format=:tsv,
                gzip_read_map_tsv=false,
                gzip_prefixed_fastqs=false,
                run_mapping=true,
                run_splitting=true,
                keep_prefixed_fastqs=true,
                force=true
            )

            try
                Test.@test joint_res.joint_read_map !== nothing
                Test.@test isfile(joint_res.joint_read_map)
                uuid_to_source = Dict{String,String}()
                for sample in joint_res.sample_outputs
                    for map_path in sample.mapping_files
                        for (i, line) in enumerate(eachline(map_path))
                            i == 1 && continue
                            isempty(line) && continue
                            fields = split(line, '\t')
                            uuid = strip(fields[end])
                            source = strip(fields[1])
                            uuid_to_source[uuid] = source
                        end
                    end
                end
                for sample in joint_res.sample_outputs
                    Test.@test isfile(sample.output_bam)
                    reader = Mycelia.open_xam(sample.output_bam; parser=:samtools)
                    for record in reader
                        uuid = strip(XAM.SAM.tempname(record))
                        Test.@test uuid_to_source[uuid] in sample.source_fastqs
                    end
                    close(reader)
                end
            finally
                isdir(joint_res.tmpdir) && rm(joint_res.tmpdir; recursive=true, force=true)
            end

            fasta_res = Mycelia.minimap_merge_map_and_split(
                reference_fasta=ref,
                minimap_index="",
                mapping_type="sr",
                single_end_fastqs=[fq_a],
                read_id_strategy=:uuid,
                fastq_mode=:per_sample,
                read_map_format=:tsv,
                gzip_read_map_tsv=false,
                gzip_prefixed_fastqs=false,
                run_mapping=true,
                run_splitting=true,
                build_index=false
            )
            Test.@test isfile(first(fasta_res.sample_outputs).output_bam)
        end
    end

    Test.@testset "NCBI genome mapping consistency (short/long, index/no index)" begin
        run_all = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true"
        run_external = run_all || lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"
        if !run_external
            @info "Skipping NCBI mapping tests; set MYCELIA_RUN_EXTERNAL=true to run"
            Test.@test_skip "Requires NCBI download + minimap2 + samtools"
        else
            genome_info = Mycelia.get_test_genome_fasta(use_ncbi=true)
            try
                if genome_info.source != :ncbi
                    @info "NCBI genome download failed; skipping mapping consistency tests"
                    Test.@test_skip "NCBI download unavailable"
                else
                    function read_reference_sequence(path::String)
                        reader = Mycelia.open_fastx(path)
                        record = first(reader)
                        close(reader)
                        return FASTX.sequence(String, record)
                    end

                    function build_records(seq::String, positions::Vector{Int}, read_len::Int, prefix::String)
                        records = FASTX.FASTQ.Record[]
                        for (i, pos) in enumerate(positions)
                            subseq = seq[pos:pos + read_len - 1]
                            push!(records, FASTX.FASTQ.Record("$(prefix)_$(i)", subseq, repeat("I", read_len)))
                        end
                        return records
                    end

                    function ensure_index(ref::String, mapping_type::String)
                        idx_res = Mycelia.minimap_index(
                            fasta=ref,
                            mapping_type=mapping_type,
                            mem_gb=1,
                            threads=1,
                            as_string=false
                        )
                        if !isfile(idx_res.outfile)
                            run(idx_res.cmd)
                        end
                        return idx_res.outfile
                    end

                    function collect_primary_mappings(bam::String; id_map::Dict{String,String}=Dict{String,String}())
                        reader = Mycelia.open_xam(bam; parser=:samtools)
                        df = Mycelia.xam_to_dataframe(reader)
                        close(reader)
                        df = df[df.ismapped .& df.isprimary, :]
                        tuples = Tuple{String,String,Int,Int}[]
                        for row in eachrow(df)
                            read_key = String(row.template)
                            read_id = get(id_map, read_key, read_key)
                            push!(tuples, (read_id, row.reference, first(row.position), row.alignlength))
                        end
                        return sort(tuples)
                    end

                    function mapping_from_joint_tsv(path::String)
                        mapping_by_source = Dict{String,Dict{String,String}}()
                        for (i, line) in enumerate(eachline(path))
                            i == 1 && continue
                            isempty(line) && continue
                            fields = split(line, '\t')
                            if length(fields) < 3
                                continue
                            end
                            source = String(strip(fields[1]))
                            original = String(strip(fields[2]))
                            uuid = String(strip(fields[3]))
                            if !haskey(mapping_by_source, source)
                                mapping_by_source[source] = Dict{String,String}()
                            end
                            mapping_by_source[source][uuid] = original
                        end
                        return mapping_by_source
                    end

                    function flatten_read_map(mapping_by_source::Dict{String,Dict{String,String}})
                        flat = Dict{String,String}()
                        for sample_map in values(mapping_by_source)
                            merge!(flat, sample_map)
                        end
                        return flat
                    end

                    function strip_read_ids(tuples::Vector{Tuple{String,String,Int,Int}})
                        stripped = Tuple{String,Int,Int}[]
                        for (read_id, reference, pos, alignlen) in tuples
                            push!(stripped, (reference, pos, alignlen))
                        end
                        return sort(stripped)
                    end

                    function run_consistency_case(;ref, seq, mapping_type, read_len, use_index)
                        max_start = length(seq) - read_len + 1
                        positions_a = [1, min(1 + read_len, max_start)]
                        positions_b = [min(1 + 3 * read_len, max_start), min(1 + 5 * read_len, max_start)]
                        reads_a = build_records(seq, positions_a, read_len, "sampleA")
                        reads_b = build_records(seq, positions_b, read_len, "sampleB")

                        case_dir = mktempdir()
                        fq_a = joinpath(case_dir, "sampleA.$(mapping_type).fastq")
                        fq_b = joinpath(case_dir, "sampleB.$(mapping_type).fastq")
                        Mycelia.write_fastq(records=reads_a, filename=fq_a)
                        Mycelia.write_fastq(records=reads_b, filename=fq_b)

                        idx = use_index ? ensure_index(ref, mapping_type) : ""

                        function map_single(fq::String)
                            if use_index
                                res = Mycelia.minimap_map_with_index(
                                    fasta=nothing,
                                    mapping_type=mapping_type,
                                    fastq=fq,
                                    index_file=idx,
                                    sorted=true,
                                    threads=1,
                                    mem_gb=1
                                )
                            else
                                res = Mycelia.minimap_map(
                                    fasta=ref,
                                    fastq=fq,
                                    mapping_type=mapping_type,
                                    sorted=true,
                                    output_format="bam",
                                    threads=1,
                                    mem_gb=1
                                )
                            end
                            run(res.cmd)
                            return res.outfile
                        end

                        expected_a = collect_primary_mappings(map_single(fq_a))
                        expected_b = collect_primary_mappings(map_single(fq_b))

                        single_joint = Mycelia.minimap_merge_map_and_split(
                            reference_fasta=use_index ? nothing : ref,
                            minimap_index=use_index ? idx : "",
                            mapping_type=mapping_type,
                            single_end_fastqs=[fq_a],
                            read_id_strategy=:uuid,
                            fastq_mode=:joint,
                            read_map_format=:tsv,
                            gzip_read_map_tsv=false,
                            gzip_prefixed_fastqs=false,
                            run_mapping=true,
                            run_splitting=true,
                            threads=1,
                            mem_gb=1
                        )
                        Test.@test single_joint.joint_read_map !== nothing
                        mapping_single = mapping_from_joint_tsv(single_joint.joint_read_map)
                        flat_single = flatten_read_map(mapping_single)
                        single_sample = first(single_joint.sample_outputs)
                        single_map = collect_primary_mappings(
                            single_sample.output_bam;
                            id_map=flat_single
                        )
                        Test.@test strip_read_ids(single_map) == strip_read_ids(expected_a)

                        dual_joint = Mycelia.minimap_merge_map_and_split(
                            reference_fasta=use_index ? nothing : ref,
                            minimap_index=use_index ? idx : "",
                            mapping_type=mapping_type,
                            single_end_fastqs=[fq_a, fq_b],
                            read_id_strategy=:uuid,
                            fastq_mode=:joint,
                            read_map_format=:tsv,
                            gzip_read_map_tsv=false,
                            gzip_prefixed_fastqs=false,
                            run_mapping=true,
                            run_splitting=true,
                            threads=1,
                            mem_gb=1
                        )
                        Test.@test dual_joint.joint_read_map !== nothing
                        mapping_dual = mapping_from_joint_tsv(dual_joint.joint_read_map)
                        flat_dual = flatten_read_map(mapping_dual)
                        for sample in dual_joint.sample_outputs
                            source = first(sample.source_fastqs)
                            actual = collect_primary_mappings(sample.output_bam; id_map=flat_dual)
                            expected = source == fq_a ? expected_a : expected_b
                            Test.@test strip_read_ids(actual) == strip_read_ids(expected)
                        end
                    end

                    ref = genome_info.fasta
                    seq = read_reference_sequence(ref)
                    for (read_len, mapping_type, label) in [
                        (100, "sr", "short"),
                        (400, "map-hifi", "long")
                    ]
                        for use_index in (false, true)
                            Test.@testset "$(label) reads, $(use_index ? "prebuilt index" : "reference")" begin
                                run_consistency_case(
                                    ref=ref,
                                    seq=seq,
                                    mapping_type=mapping_type,
                                    read_len=read_len,
                                    use_index=use_index
                                )
                            end
                        end
                    end
                end
            finally
                genome_info.cleanup()
            end
        end
    end
end
