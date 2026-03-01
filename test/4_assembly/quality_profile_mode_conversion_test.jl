import Test
import FASTX
import BioSequences
import Graphs
import MetaGraphsNext
import Mycelia

Test.@testset "Quality Profile Mode Conversion Tests (FASTQ-based)" begin
    quality_profiles = [:ultralight_quality, :lightweight_quality]

    # --- Test 1: Quality reversal invariant for doublestrand ---
    # For every non-palindrome kmer in a doublestrand graph,
    # quality(kmer) == reverse(quality(RC(kmer)))
    Test.@testset "Doublestrand quality reversal invariant" begin
        # Simple case: no forward kmer equals any other's RC
        seq_simple = "ATGC"
        qual_simple = UInt8[30, 35, 32, 28]
        qual_str_simple = String([Char(q + 33) for q in qual_simple])
        records_simple = [FASTX.FASTQ.Record("seq1", seq_simple, qual_str_simple)]

        # Complex case: forward kmers overlap with RC kmers (ATG<->CAT, TGC<->GCA)
        seq_overlap = "ATGCAT"
        qual_overlap = UInt8[30, 35, 32, 28, 40, 35]
        qual_str_overlap = String([Char(q + 33) for q in qual_overlap])
        records_overlap = [FASTX.FASTQ.Record("seq1", seq_overlap, qual_str_overlap)]

        for profile in quality_profiles
            Test.@testset "$profile simple (no RC overlap)" begin
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_simple, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                for label in MetaGraphsNext.labels(ds)
                    rc_label = BioSequences.reverse_complement(label)
                    if label != rc_label && haskey(ds, rc_label)
                        fwd_q = collect(ds[label].joint_quality)
                        rc_q = collect(ds[rc_label].joint_quality)
                        Test.@test fwd_q == reverse(rc_q)
                    end
                end
            end

            Test.@testset "$profile complex (with RC overlap)" begin
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records_overlap, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                for label in MetaGraphsNext.labels(ds)
                    rc_label = BioSequences.reverse_complement(label)
                    if label != rc_label && haskey(ds, rc_label)
                        fwd_q = collect(ds[label].joint_quality)
                        rc_q = collect(ds[rc_label].joint_quality)
                        Test.@test fwd_q == reverse(rc_q)
                    end
                end
            end
        end
    end

    # --- Test 2: Quality vector length == k ---
    Test.@testset "Quality vector length equals k" begin
        seq = "ATGCATGC"
        qual = UInt8[30, 35, 32, 28, 40, 35, 30, 33]
        qual_str = String([Char(q + 33) for q in qual])
        records = [FASTX.FASTQ.Record("seq1", seq, qual_str)]
        k = 3

        for profile in quality_profiles
            for mode in [:singlestrand, :doublestrand, :canonical]
                Test.@testset "$profile $mode" begin
                    g = Mycelia.Rhizomorph.build_kmer_graph(
                        records, k;
                        memory_profile = profile,
                        mode = mode
                    )

                    for label in MetaGraphsNext.labels(g)
                        Test.@test length(g[label].joint_quality) == k
                    end
                end
            end
        end
    end

    # --- Test 3: Total quality sum preservation ---
    # Reversal preserves element sums; merging adds element-wise.
    # So total sum across all vertices is invariant through canonical conversion,
    # and doubles (for non-palindromic sets) through doublestrand conversion.
    Test.@testset "Total quality sum preservation" begin
        # Use two test sequences: one without RC overlap, one with
        test_cases = [
            ("no overlap", "ATGC", UInt8[30, 35, 32, 28]),
            ("with RC overlap", "ATGCAT", UInt8[30, 35, 32, 28, 40, 35])
        ]

        for (case_name, seq, qual) in test_cases
            qual_str = String([Char(q + 33) for q in qual])
            records = [FASTX.FASTQ.Record("seq1", seq, qual_str)]

            for profile in quality_profiles
                Test.@testset "$profile $case_name" begin
                    ss = Mycelia.Rhizomorph.build_kmer_graph(
                        records, 3;
                        memory_profile = profile,
                        mode = :singlestrand
                    )

                    ss_total = sum(
                        sum(Int.(ss[l].joint_quality))
                    for l in MetaGraphsNext.labels(ss)
                    )

                    # Canonical: total quality sum preserved
                    # (reversal doesn't change sum, merge adds element-wise)
                    cn = Mycelia.Rhizomorph.build_kmer_graph(
                        records, 3;
                        memory_profile = profile,
                        mode = :canonical
                    )
                    cn_total = sum(
                        sum(Int.(cn[l].joint_quality))
                    for l in MetaGraphsNext.labels(cn)
                    )
                    Test.@test cn_total == ss_total

                    # Doublestrand: total = 2x singlestrand (no palindromes at k=3)
                    ds = Mycelia.Rhizomorph.build_kmer_graph(
                        records, 3;
                        memory_profile = profile,
                        mode = :doublestrand
                    )
                    ds_total = sum(
                        sum(Int.(ds[l].joint_quality))
                    for l in MetaGraphsNext.labels(ds)
                    )
                    Test.@test ds_total == 2 * ss_total
                end
            end
        end
    end

    # --- Test 4: Forward quality preserved in singlestrand→doublestrand ---
    # When an RC kmer does NOT exist as a forward kmer, the forward quality
    # is unchanged and the RC quality is exactly the reverse.
    Test.@testset "Forward quality unchanged, RC quality exactly reversed" begin
        # ATGC k=3: ATG and TGC are forward; CAT and GCA are new RC-only
        seq = "ATGC"
        qual = UInt8[30, 35, 32, 28]
        qual_str = String([Char(q + 33) for q in qual])
        records = [FASTX.FASTQ.Record("seq1", seq, qual_str)]

        for profile in quality_profiles
            Test.@testset "$profile" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                for fwd_label in MetaGraphsNext.labels(ss)
                    rc_label = BioSequences.reverse_complement(fwd_label)

                    # Forward quality should be unchanged in doublestrand
                    fwd_ss_q = collect(ss[fwd_label].joint_quality)
                    fwd_ds_q = collect(ds[fwd_label].joint_quality)
                    Test.@test fwd_ds_q == fwd_ss_q

                    # RC quality should be exactly reverse of forward singlestrand
                    rc_ds_q = collect(ds[rc_label].joint_quality)
                    Test.@test rc_ds_q == reverse(fwd_ss_q)
                end
            end
        end
    end

    # --- Test 5: Edge quality swap in doublestrand ---
    # For an edge src→dst, the RC edge RC(dst)→RC(src) should have:
    #   RC_edge.from_quality == reverse(fwd_edge.to_quality)
    #   RC_edge.to_quality   == reverse(fwd_edge.from_quality)
    Test.@testset "Edge quality swap for RC edges" begin
        seq = "ATGC"
        qual = UInt8[30, 35, 32, 28]
        qual_str = String([Char(q + 33) for q in qual])
        records = [FASTX.FASTQ.Record("seq1", seq, qual_str)]

        for profile in quality_profiles
            Test.@testset "$profile" begin
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                # For each forward edge, check its RC edge
                for (src, dst) in MetaGraphsNext.edge_labels(ss)
                    rc_src = BioSequences.reverse_complement(dst)
                    rc_dst = BioSequences.reverse_complement(src)

                    fwd_edge = ds[src, dst]
                    Test.@test haskey(ds, rc_src, rc_dst)
                    rc_edge = ds[rc_src, rc_dst]

                    # Verify swap + reverse
                    Test.@test collect(rc_edge.from_joint_quality) ==
                               reverse(collect(fwd_edge.to_joint_quality))
                    Test.@test collect(rc_edge.to_joint_quality) ==
                               reverse(collect(fwd_edge.from_joint_quality))
                end
            end
        end
    end

    # --- Test 6: Dataset-level quality consistency ---
    # With a single dataset, dataset_joint_quality should match joint_quality
    Test.@testset "Dataset quality matches joint quality (single dataset)" begin
        seq = "ATGC"
        qual = UInt8[30, 35, 32, 28]
        qual_str = String([Char(q + 33) for q in qual])
        records = [FASTX.FASTQ.Record("seq1", seq, qual_str)]

        for profile in quality_profiles
            for mode in [:doublestrand, :canonical]
                Test.@testset "$profile $mode" begin
                    g = Mycelia.Rhizomorph.build_kmer_graph(
                        records, 3;
                        memory_profile = profile,
                        mode = mode,
                        dataset_id = "test_ds"
                    )

                    for label in MetaGraphsNext.labels(g)
                        vdata = g[label]
                        Test.@test haskey(vdata.dataset_joint_quality, "test_ds")
                        ds_qual = collect(vdata.dataset_joint_quality["test_ds"])
                        joint_qual = collect(vdata.joint_quality)
                        Test.@test ds_qual == joint_qual
                    end
                end
            end
        end
    end

    # --- Test 7: Lightweight quality preserves observations AND quality ---
    Test.@testset "Lightweight quality has observations and quality" begin
        seq = "ATGC"
        qual = UInt8[30, 35, 32, 28]
        qual_str = String([Char(q + 33) for q in qual])
        records = [FASTX.FASTQ.Record("seq1", seq, qual_str)]

        for mode in [:doublestrand, :canonical]
            Test.@testset "lightweight_quality $mode" begin
                g = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = :lightweight_quality,
                    mode = mode,
                    dataset_id = "test_ds"
                )

                for label in MetaGraphsNext.labels(g)
                    vdata = g[label]
                    # Has observations (lightweight feature)
                    Test.@test hasfield(typeof(vdata), :dataset_observations)
                    Test.@test haskey(vdata.dataset_observations, "test_ds")
                    Test.@test !isempty(vdata.dataset_observations["test_ds"])
                    # Has quality (quality feature)
                    Test.@test hasfield(typeof(vdata), :joint_quality)
                    Test.@test !isempty(vdata.joint_quality)
                    Test.@test length(vdata.joint_quality) == 3
                end
            end
        end
    end

    # --- Test 8: Canonical labels are canonical ---
    Test.@testset "All canonical graph labels are canonical form" begin
        seq = "ATGCATGC"
        qual = UInt8[30, 35, 32, 28, 40, 35, 30, 33]
        qual_str = String([Char(q + 33) for q in qual])
        records = [FASTX.FASTQ.Record("seq1", seq, qual_str)]

        for profile in quality_profiles
            Test.@testset "$profile" begin
                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :canonical
                )

                for label in MetaGraphsNext.labels(cn)
                    Test.@test BioSequences.canonical(label) == label
                end

                # Canonical graph should be undirected
                Test.@test cn.graph isa Graphs.SimpleGraph
            end
        end
    end

    # --- Test 9: Quality topology matches full mode ---
    Test.@testset "Quality profiles match full mode topology" begin
        seq = "ATGCATGC"
        qual = UInt8[30, 35, 32, 28, 40, 35, 30, 33]
        qual_str = String([Char(q + 33) for q in qual])
        records = [FASTX.FASTQ.Record("seq1", seq, qual_str)]

        full_ds = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3;
            memory_profile = :full,
            mode = :doublestrand
        )
        full_cn = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3;
            memory_profile = :full,
            mode = :canonical
        )

        for profile in quality_profiles
            Test.@testset "$profile doublestrand topology" begin
                g = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )
                Test.@test Graphs.nv(g.graph) == Graphs.nv(full_ds.graph)
                Test.@test Graphs.ne(g.graph) == Graphs.ne(full_ds.graph)
            end

            Test.@testset "$profile canonical topology" begin
                g = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :canonical
                )
                Test.@test Graphs.nv(g.graph) == Graphs.nv(full_cn.graph)
                Test.@test Graphs.ne(g.graph) == Graphs.ne(full_cn.graph)
            end
        end
    end

    # --- Test 10: Multi-record quality merging through conversion ---
    Test.@testset "Multi-record quality merging" begin
        seq1 = "ATGC"
        qual1 = UInt8[30, 35, 32, 28]
        seq2 = "ATGC"
        qual2 = UInt8[25, 30, 28, 32]
        records = [
            FASTX.FASTQ.Record("seq1", seq1, String([Char(q + 33) for q in qual1])),
            FASTX.FASTQ.Record("seq2", seq2, String([Char(q + 33) for q in qual2]))
        ]

        for profile in quality_profiles
            Test.@testset "$profile" begin
                # Singlestrand: ATG seen twice, quality merged
                ss = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :singlestrand
                )

                ss_total = sum(
                    sum(Int.(ss[l].joint_quality))
                for l in MetaGraphsNext.labels(ss)
                )

                # Doublestrand symmetry still holds with merged qualities
                ds = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :doublestrand
                )

                for label in MetaGraphsNext.labels(ds)
                    rc_label = BioSequences.reverse_complement(label)
                    if label != rc_label && haskey(ds, rc_label)
                        fwd_q = collect(ds[label].joint_quality)
                        rc_q = collect(ds[rc_label].joint_quality)
                        Test.@test fwd_q == reverse(rc_q)
                    end
                end

                # Canonical sum preservation
                cn = Mycelia.Rhizomorph.build_kmer_graph(
                    records, 3;
                    memory_profile = profile,
                    mode = :canonical
                )
                cn_total = sum(
                    sum(Int.(cn[l].joint_quality))
                for l in MetaGraphsNext.labels(cn)
                )
                Test.@test cn_total == ss_total
            end
        end
    end
end
