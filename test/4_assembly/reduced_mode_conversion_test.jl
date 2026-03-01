import Test
import FASTX
import BioSequences
import Graphs
import Kmers
import MetaGraphsNext
import Mycelia

Test.@testset "Reduced Memory Profile Mode Conversions" begin
    # Setup: small DNA sequences with known RC properties
    # ATGCATGC has interesting RC structure:
    # k=3 forward: ATG, TGC, GCA, CAT, ATG, TGC -> unique: ATG, TGC, GCA, CAT
    # RC of ATGCATGC = GCATGCAT
    # k=3 RC:       GCA, CAT, ATG, TGC, GCA, CAT -> same set
    records_fasta = [FASTX.FASTA.Record("seq1", "ATGCATGC")]

    # FASTQ records needed for quality profiles
    seq = "ATGCATGC"
    qual_str = String([Char(q + 33) for q in [30, 35, 32, 28, 40, 35, 30, 33]])
    records_fastq = [FASTX.FASTQ.Record("seq1", seq, qual_str)]

    non_quality_profiles = [:lightweight, :ultralight]
    quality_profiles = [:ultralight_quality, :lightweight_quality]

    Test.@testset "Doublestrand conversion" begin
        for profile in non_quality_profiles
            Test.@testset "$(profile) doublestrand" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                # Doublestrand should have >= vertices (forward + RC)
                Test.@test Graphs.nv(ds.graph) >= Graphs.nv(ss.graph)

                # Both should be directed graphs
                Test.@test ds.graph isa Graphs.DiGraph
                Test.@test ss.graph isa Graphs.DiGraph

                # Every forward kmer should be in doublestrand graph
                for label in MetaGraphsNext.labels(ss)
                    Test.@test haskey(ds, label)
                end

                # RC of every forward kmer should also be in doublestrand graph
                for label in MetaGraphsNext.labels(ss)
                    rc_label = BioSequences.reverse_complement(label)
                    Test.@test haskey(ds, rc_label)
                end

                # Doublestrand should have >= edges
                Test.@test Graphs.ne(ds.graph) >= Graphs.ne(ss.graph)
            end
        end

        for profile in quality_profiles
            Test.@testset "$(profile) doublestrand" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                Test.@test Graphs.nv(ds.graph) >= Graphs.nv(ss.graph)
                Test.@test ds.graph isa Graphs.DiGraph

                # All forward kmers present
                for label in MetaGraphsNext.labels(ss)
                    Test.@test haskey(ds, label)
                end

                # All RC kmers present
                for label in MetaGraphsNext.labels(ss)
                    rc_label = BioSequences.reverse_complement(label)
                    Test.@test haskey(ds, rc_label)
                end
            end
        end
    end

    Test.@testset "Canonical conversion" begin
        for profile in non_quality_profiles
            Test.@testset "$(profile) canonical" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :canonical
                )

                # Canonical should have <= vertices (merging RC pairs)
                Test.@test Graphs.nv(cn.graph) <= Graphs.nv(ss.graph)

                # Canonical should be undirected
                Test.@test cn.graph isa Graphs.SimpleGraph

                # Every canonical label should be the canonical form of itself
                for label in MetaGraphsNext.labels(cn)
                    Test.@test BioSequences.canonical(label) == label
                end
            end
        end

        for profile in quality_profiles
            Test.@testset "$(profile) canonical" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :canonical
                )

                Test.@test Graphs.nv(cn.graph) <= Graphs.nv(ss.graph)
                Test.@test cn.graph isa Graphs.SimpleGraph

                for label in MetaGraphsNext.labels(cn)
                    Test.@test BioSequences.canonical(label) == label
                end
            end
        end
    end

    Test.@testset "Count preservation through doublestrand conversion" begin
        simple_records = [FASTX.FASTA.Record("seq1", "ATGATG")]

        for profile in non_quality_profiles
            Test.@testset "$(profile) count preservation" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    simple_records, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    simple_records, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                ss_total = sum(
                    ss[l].total_count for l in MetaGraphsNext.labels(ss)
                )
                ds_total = sum(
                    ds[l].total_count for l in MetaGraphsNext.labels(ds)
                )

                # DS should have >= total count (forward + RC copies)
                Test.@test ds_total >= ss_total
            end
        end
    end

    Test.@testset "Count preservation through canonical conversion" begin
        simple_records = [FASTX.FASTA.Record("seq1", "ATGATG")]

        for profile in non_quality_profiles
            Test.@testset "$(profile) canonical count preservation" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    simple_records, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    simple_records, 3;
                    memory_profile = profile,
                    mode = :canonical
                )

                ss_total = sum(
                    ss[l].total_count for l in MetaGraphsNext.labels(ss)
                )
                cn_total = sum(
                    cn[l].total_count for l in MetaGraphsNext.labels(cn)
                )

                # Canonical merges RC pairs, so total count should be preserved
                # (each pair sums their counts)
                Test.@test cn_total >= ss_total
            end
        end
    end

    Test.@testset "Dataset counts preserved through conversion" begin
        records = [FASTX.FASTA.Record("seq1", "ATGATG")]

        for profile in non_quality_profiles
            Test.@testset "$(profile) dataset counts" begin
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :doublestrand,
                    dataset_id = "test_ds"
                )

                # All vertices should have dataset_counts for "test_ds"
                for label in MetaGraphsNext.labels(ds)
                    vdata = ds[label]
                    Test.@test haskey(vdata.dataset_counts, "test_ds")
                    Test.@test vdata.dataset_counts["test_ds"] > 0
                end
            end
        end
    end

    Test.@testset "Observation IDs preserved in lightweight conversions" begin
        records = [FASTX.FASTA.Record("seq1", "ATGATG")]

        for profile in [:lightweight, :lightweight_quality]
            Test.@testset "$(profile) observation IDs" begin
                input_records = profile == :lightweight_quality ?
                                [FASTX.FASTQ.Record("seq1", "ATGATG",
                    String([Char(q + 33) for q in [30, 35, 32, 28, 40, 35]]))] :
                                records

                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    input_records, 3;
                    memory_profile = profile,
                    mode = :doublestrand,
                    dataset_id = "test_ds"
                )

                # All vertices should have observation IDs
                for label in MetaGraphsNext.labels(ds)
                    vdata = ds[label]
                    if hasfield(typeof(vdata), :dataset_observations)
                        Test.@test haskey(vdata.dataset_observations, "test_ds")
                        Test.@test !isempty(vdata.dataset_observations["test_ds"])
                    end
                end
            end
        end
    end

    Test.@testset "Quality vectors present in quality profile conversions" begin
        for profile in quality_profiles
            Test.@testset "$(profile) quality in doublestrand" begin
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :doublestrand,
                    dataset_id = "test_ds"
                )

                # All vertices should have quality vectors
                for label in MetaGraphsNext.labels(ds)
                    vdata = ds[label]
                    Test.@test hasfield(typeof(vdata), :joint_quality)
                    Test.@test !isempty(vdata.joint_quality)
                end
            end

            Test.@testset "$(profile) quality in canonical" begin
                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :canonical,
                    dataset_id = "test_ds"
                )

                for label in MetaGraphsNext.labels(cn)
                    vdata = cn[label]
                    Test.@test hasfield(typeof(vdata), :joint_quality)
                    Test.@test !isempty(vdata.joint_quality)
                end
            end
        end
    end

    Test.@testset "Vertex data types preserved through conversion" begin
        type_map = Dict(
            :lightweight => Mycelia.Rhizomorph.LightweightKmerVertexData,
            :ultralight => Mycelia.Rhizomorph.UltralightKmerVertexData
        )
        quality_type_map = Dict(
            :ultralight_quality => Mycelia.Rhizomorph.UltralightQualityKmerVertexData,
            :lightweight_quality => Mycelia.Rhizomorph.LightweightQualityKmerVertexData
        )

        for (profile, expected_type) in type_map
            Test.@testset "$(profile) type preserved" begin
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )
                for label in MetaGraphsNext.labels(ds)
                    Test.@test ds[label] isa expected_type
                end

                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :canonical
                )
                for label in MetaGraphsNext.labels(cn)
                    Test.@test cn[label] isa expected_type
                end
            end
        end

        for (profile, expected_type) in quality_type_map
            Test.@testset "$(profile) type preserved" begin
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )
                for label in MetaGraphsNext.labels(ds)
                    Test.@test ds[label] isa expected_type
                end

                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :canonical
                )
                for label in MetaGraphsNext.labels(cn)
                    Test.@test cn[label] isa expected_type
                end
            end
        end
    end

    Test.@testset "Singlestrand mode unchanged by this work" begin
        # Verify singlestrand still works for all reduced profiles
        for profile in non_quality_profiles
            Test.@testset "$(profile) singlestrand" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                Test.@test Graphs.nv(ss.graph) > 0
                Test.@test ss.graph isa Graphs.DiGraph
            end
        end

        for profile in quality_profiles
            Test.@testset "$(profile) singlestrand" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                Test.@test Graphs.nv(ss.graph) > 0
                Test.@test ss.graph isa Graphs.DiGraph
            end
        end
    end

    Test.@testset "Comparison with full mode topology" begin
        # Full mode and reduced mode doublestrand should have same vertex count
        full_ds = Mycelia.Rhizomorph.build_kmer_graph(
            records_fasta, 3;
            memory_profile = :full,
            mode = :doublestrand
        )
        for profile in non_quality_profiles
            Test.@testset "$(profile) vs full doublestrand topology" begin
                reduced_ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )
                Test.@test Graphs.nv(reduced_ds.graph) == Graphs.nv(full_ds.graph)
                Test.@test Graphs.ne(reduced_ds.graph) == Graphs.ne(full_ds.graph)
            end
        end

        full_cn = Mycelia.Rhizomorph.build_kmer_graph(
            records_fasta, 3;
            memory_profile = :full,
            mode = :canonical
        )
        for profile in non_quality_profiles
            Test.@testset "$(profile) vs full canonical topology" begin
                reduced_cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fasta, 3;
                    memory_profile = profile,
                    mode = :canonical
                )
                Test.@test Graphs.nv(reduced_cn.graph) == Graphs.nv(full_cn.graph)
                Test.@test Graphs.ne(reduced_cn.graph) == Graphs.ne(full_cn.graph)
            end
        end
    end

    # Exact count verification for non-quality profiles with RC overlap.
    # When forward kmers overlap with their RC counterparts (e.g., ATG↔CAT),
    # doublestrand conversion must produce symmetric counts and correct totals.
    Test.@testset "Doublestrand count symmetry with RC overlap" begin
        # ATGCAT k=3: ATG↔CAT and TGC↔GCA are RC pairs that all appear forward
        records_overlap = [FASTX.FASTA.Record("seq1", "ATGCAT")]

        for profile in non_quality_profiles
            Test.@testset "$(profile) count symmetry" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_overlap, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_overlap, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                # Every kmer/RC pair should have equal total_count
                for label in MetaGraphsNext.labels(ds)
                    rc_label = BioSequences.reverse_complement(label)
                    if label != rc_label && haskey(ds, rc_label)
                        Test.@test ds[label].total_count == ds[rc_label].total_count
                    end
                end

                # Doublestrand total count = 2x singlestrand (no palindromes at k=3)
                ss_total = sum(
                    ss[l].total_count for l in MetaGraphsNext.labels(ss)
                )
                ds_total = sum(
                    ds[l].total_count for l in MetaGraphsNext.labels(ds)
                )
                Test.@test ds_total == 2 * ss_total
            end
        end
    end

    Test.@testset "Canonical count preservation with RC overlap" begin
        records_overlap = [FASTX.FASTA.Record("seq1", "ATGCAT")]

        for profile in non_quality_profiles
            Test.@testset "$(profile) canonical count total" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_overlap, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records_overlap, 3;
                    memory_profile = profile,
                    mode = :canonical
                )

                ss_total = sum(
                    ss[l].total_count for l in MetaGraphsNext.labels(ss)
                )
                cn_total = sum(
                    cn[l].total_count for l in MetaGraphsNext.labels(cn)
                )

                # Canonical total count == singlestrand total
                # (RC pair counts merge, not double)
                Test.@test cn_total == ss_total
            end
        end
    end

    Test.@testset "Doublestrand dataset_counts symmetry with RC overlap" begin
        records_overlap = [FASTX.FASTA.Record("seq1", "ATGCAT")]

        for profile in non_quality_profiles
            Test.@testset "$(profile) dataset count symmetry" begin
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_overlap, 3;
                    memory_profile = profile,
                    mode = :doublestrand,
                    dataset_id = "test_ds"
                )

                for label in MetaGraphsNext.labels(ds)
                    rc_label = BioSequences.reverse_complement(label)
                    if label != rc_label && haskey(ds, rc_label)
                        fwd_ds_count = ds[label].dataset_counts["test_ds"]
                        rc_ds_count = ds[rc_label].dataset_counts["test_ds"]
                        Test.@test fwd_ds_count == rc_ds_count
                    end
                end
            end
        end
    end

    Test.@testset "Doublestrand edge symmetry" begin
        # Non-overlapping case: ATGC k=3 → edges: ATG→TGC
        # RC edge: GCA→CAT (new, doesn't overlap forward)
        # So ne(ds) == 2 * ne(ss)
        records_no_overlap = [FASTX.FASTA.Record("seq1", "ATGC")]

        for profile in non_quality_profiles
            Test.@testset "$(profile) no-overlap edges 2x" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_no_overlap, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_no_overlap, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                Test.@test Graphs.ne(ds.graph) == 2 * Graphs.ne(ss.graph)
            end
        end

        # Overlapping case: ATGCAT k=3 → edges: ATG→TGC, TGC→GCA, GCA→CAT
        # RC edges: GCA→CAT, TGC→GCA, ATG→TGC — all map to existing forward edges
        # So ne(ds) == ne(ss) (all RC edges merge into forward edges)
        records_overlap = [FASTX.FASTA.Record("seq1", "ATGCAT")]

        for profile in non_quality_profiles
            Test.@testset "$(profile) self-RC-overlap edges merge" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_overlap, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_overlap, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                # Same edge count (RC edges merge into existing forward edges)
                Test.@test Graphs.ne(ds.graph) == Graphs.ne(ss.graph)

                # Every forward edge has a corresponding RC edge present
                for (src, dst) in MetaGraphsNext.edge_labels(ss)
                    rc_src = BioSequences.reverse_complement(dst)
                    rc_dst = BioSequences.reverse_complement(src)
                    Test.@test haskey(ds, rc_src, rc_dst)
                end

                # Merged edges have 2x the original count
                for (src, dst) in MetaGraphsNext.edge_labels(ss)
                    Test.@test ds[src, dst].total_count == 2 * ss[src, dst].total_count
                end
            end
        end

        # Quality profiles: same invariants with FASTQ input
        qual = UInt8[30, 35, 32, 28, 40, 35]
        qual_str = String([Char(q + 33) for q in qual])
        records_fastq_overlap = [FASTX.FASTQ.Record("seq1", "ATGCAT", qual_str)]

        for profile in quality_profiles
            Test.@testset "$(profile) self-RC-overlap edges merge" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq_overlap, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_fastq_overlap, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                Test.@test Graphs.ne(ds.graph) == Graphs.ne(ss.graph)

                for (src, dst) in MetaGraphsNext.edge_labels(ss)
                    rc_src = BioSequences.reverse_complement(dst)
                    rc_dst = BioSequences.reverse_complement(src)
                    Test.@test haskey(ds, rc_src, rc_dst)
                end

                # Merged edges: 2x count
                for (src, dst) in MetaGraphsNext.edge_labels(ss)
                    Test.@test ds[src, dst].total_count == 2 * ss[src, dst].total_count
                end
            end
        end
    end

    Test.@testset "is_self_rc_closed detects self-complementary edge sets" begin
        # ATGCAT k=3: all RC edges map to existing forward edges
        records_closed = [FASTX.FASTA.Record("seq1", "ATGCAT")]
        # ATGC k=3: RC edges don't exist in forward graph
        records_open = [FASTX.FASTA.Record("seq1", "ATGC")]

        for profile in [non_quality_profiles; quality_profiles]
            input = profile in quality_profiles ?
                    [FASTX.FASTQ.Record("seq1", "ATGCAT",
                String([Char(q + 33) for q in UInt8[30, 35, 32, 28, 40, 35]]))] :
                    records_closed
            input_open = profile in quality_profiles ?
                         [FASTX.FASTQ.Record("seq1", "ATGC",
                String([Char(q + 33) for q in UInt8[30, 35, 32, 28]]))] :
                         records_open

            Test.@testset "$(profile)" begin
                # Self-complementary sequence → closed under RC
                ss_closed = Mycelia.Rhizomorph.build_kmer_graph(
                    input, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                Test.@test Mycelia.Rhizomorph.is_self_rc_closed(ss_closed) == true

                # Non-self-complementary → not closed
                ss_open = Mycelia.Rhizomorph.build_kmer_graph(
                    input_open, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                Test.@test Mycelia.Rhizomorph.is_self_rc_closed(ss_open) == false

                # Doublestrand graphs are always closed under RC
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    input_open, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )
                Test.@test Mycelia.Rhizomorph.is_self_rc_closed(ds) == true

                # Canonical (undirected) returns false
                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    input, 3;
                    memory_profile = profile,
                    mode = :canonical
                )
                Test.@test Mycelia.Rhizomorph.is_self_rc_closed(cn) == false
            end
        end
    end

    Test.@testset "build_kmer_graph_from_files with reduced profiles" begin
        mktempdir() do dir
            # Write two FASTA files
            path1 = joinpath(dir, "sample_A.fasta")
            path2 = joinpath(dir, "sample_B.fasta")
            open(path1, "w") do io
                println(io, ">seq1")
                println(io, "ATGCATGC")
            end
            open(path2, "w") do io
                println(io, ">seq2")
                println(io, "GCATGCAT")
            end

            for profile in non_quality_profiles
                Test.@testset "$(profile) from_files doublestrand" begin
                    ds = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                        [path1, path2], 3;
                        mode = :doublestrand,
                        memory_profile = profile
                    )
                    Test.@test Graphs.nv(ds.graph) > 0
                    Test.@test ds.graph isa Graphs.DiGraph

                    # Should have data from both datasets
                    all_ds_ids = Set{String}()
                    for label in MetaGraphsNext.labels(ds)
                        for ds_id in keys(ds[label].dataset_counts)
                            push!(all_ds_ids, ds_id)
                        end
                    end
                    Test.@test length(all_ds_ids) >= 2
                end

                Test.@testset "$(profile) from_files canonical" begin
                    cn = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                        [path1, path2], 3;
                        mode = :canonical,
                        memory_profile = profile
                    )
                    Test.@test Graphs.nv(cn.graph) > 0
                    Test.@test cn.graph isa Graphs.SimpleGraph

                    for label in MetaGraphsNext.labels(cn)
                        Test.@test BioSequences.canonical(label) == label
                    end
                end

                Test.@testset "$(profile) from_files topology matches full" begin
                    full_ds = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                        [path1, path2], 3;
                        mode = :doublestrand,
                        memory_profile = :full
                    )
                    reduced_ds = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                        [path1, path2], 3;
                        mode = :doublestrand,
                        memory_profile = profile
                    )
                    Test.@test Graphs.nv(reduced_ds.graph) == Graphs.nv(full_ds.graph)
                    Test.@test Graphs.ne(reduced_ds.graph) == Graphs.ne(full_ds.graph)
                end
            end
        end
    end
end
