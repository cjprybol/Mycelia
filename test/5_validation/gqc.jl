import Test
import Mycelia
import Random
import JSON
import CodecZlib

Test.@testset "gqc v0.1 helpers and packet" begin

    Test.@testset "dinucleotide_shuffle preserves dinucs and ends" begin
        seq = "ATGCGATCGATCGTTACCGGATCCATGCGATTACGGCATGCGATCGATCGTTACCGGATCCAA"
        rng = Random.Xoshiro(20261021)
        shuffled = Mycelia.dinucleotide_shuffle(seq, rng)
        Test.@test length(shuffled) == length(seq)
        Test.@test shuffled[1] == seq[1]
        Test.@test shuffled[end] == seq[end]
        function dinucs(s)
            c = Dict{String, Int}()
            for i in 1:(length(s) - 1)
                d = s[i:(i + 1)]
                c[d] = get(c, d, 0) + 1
            end
            return c
        end
        Test.@test dinucs(shuffled) == dinucs(seq)
        copies = [Mycelia.dinucleotide_shuffle(seq, Random.Xoshiro(20261021 + i)) for i in 1:20]
        Test.@test any(c != seq for c in copies)
    end

    Test.@testset "C1 degenerate dinucleotide nulls do not abort" begin
        for seq in ("ACGT", "AT"^20, "AAAAAA")
            stats = Mycelia.gqc_repeat_null_stats([seq];
                k=3,
                n_shuffles=100,
                shuffle_seed=20261021,
                null_kind=:dinucleotide,
            )
            Test.@test length(stats.null_repeats) == 100
            Test.@test stats.null_degenerate === true
            Test.@test 0.0 < stats.shuffle_p_upper <= 1.0
        end
        ## One bad segment must not kill a diverse graph
        mixed = [
            "ATGCGATCGATCGTTACCGGATCCATGCGATTACGGCATGCGATCGATCGTTACCGGATCCAA",
            "ACGT",
        ]
        stats = Mycelia.gqc_repeat_null_stats(mixed;
            k=5, n_shuffles=100, shuffle_seed=1, null_kind=:dinucleotide)
        Test.@test length(stats.null_repeats) == 100
        Test.@test stats.null_degenerate === true
    end

    Test.@testset "C2 non-ACGT bases and valid-window kmer_positions" begin
        ## N / U must not KeyError
        Test.@test Mycelia.dinucleotide_shuffle("ATGNCGAT", Random.Xoshiro(1)) isa String
        Test.@test Mycelia.dinucleotide_shuffle("ACGUACGU", Random.Xoshiro(1)) isa String
        ## U treated as T → ACGTACGT dinucs
        u_shuf = Mycelia.dinucleotide_shuffle("ACGUACGU", Random.Xoshiro(2))
        Test.@test !occursin('U', u_shuf)
        ## positions only over valid ACGT windows (matches concordance)
        counts = Mycelia.gqc_repeat_counts(["ACGTNACGT"]; k=3, circular=false, canonical=false)
        ## windows: ACG,CGT,GTN,TNA,NAC,ACG,CGT — valid: ACG,CGT,ACG,CGT → positions=4, distinct=2, repeats=2
        Test.@test counts.kmer_positions == 4
        Test.@test counts.distinct_kmers == 2
        Test.@test counts.repeats == 2
        ## packet with N fixture must complete
        fixture = joinpath(@__DIR__, "..", "fixtures", "gqc_with_n")
        packet = Mycelia.run_gqc(fixture, mktempdir();
            ks=[5], n_shuffles=10, strand_mode=:single, make_plots=false,
            write_node_link_json=false, null_kinds=(:dinucleotide,))
        Test.@test isfile(joinpath(packet, "manifest.json"))
        Test.@test basename(packet) == "synth_with_n"
    end

    Test.@testset "benjamini_hochberg known values" begin
        p = [0.01, 0.04, 0.03, 0.20]
        q = Mycelia.benjamini_hochberg(p)
        Test.@test length(q) == 4
        Test.@test q[1] ≈ 0.04 atol=1e-12
        Test.@test q[4] ≈ 0.20 atol=1e-12
        Test.@test q[2] ≈ (0.04 * 4 / 3) atol=1e-12
        Test.@test q[3] ≈ (0.04 * 4 / 3) atol=1e-12
        order = sortperm(p)
        Test.@test issorted(q[order])
        Test.@test Mycelia.benjamini_hochberg(Float64[]) == Float64[]
        Test.@test Mycelia.benjamini_hochberg([0.05]) ≈ [0.05]
    end

    Test.@testset "I4 z-score gated by min_shuffle_mean" begin
        seqs = ["ATGCGATCGATCGTTACCGGATCCATGCGATTACGGCATGCGATCGATCGTTACCGGATCCAA"]
        stats = Mycelia.gqc_repeat_null_stats(seqs;
            k=11, n_shuffles=8, shuffle_seed=20261021, min_shuffle_mean=1.0,
            null_kind=:dinucleotide)
        if stats.null_mean < 1.0
            Test.@test stats.shuffle_ratio === missing
            Test.@test stats.z_score === missing
        end
    end

    Test.@testset "repeat null formulas" begin
        seqs = ["ATGCGATCGATCGTTACCGGATCCATGCGATTACGGCATGCGATCGATCGTTACCGGATCCAA"]
        stats = Mycelia.gqc_repeat_null_stats(seqs;
            k=5,
            n_shuffles=20,
            shuffle_seed=20261021,
            salt=0,
            genus_index=1,
            circular=false,
            canonical=true,
            min_shuffle_mean=1.0,
            null_kind=:dinucleotide,
        )
        Test.@test stats.observed_repeats == stats.observed_kmer_positions - stats.observed_distinct_kmers
        Test.@test length(stats.null_repeats) == 20
        Test.@test 0.0 < stats.shuffle_p_upper <= 1.0
        expected_p = (1 + count(>=(stats.observed_repeats), stats.null_repeats)) / (20 + 1)
        Test.@test stats.shuffle_p_upper ≈ expected_p
    end

    Test.@testset "I3 strand_mode ties effective canonical" begin
        Test.@test Mycelia._gqc_effective_canonical(:single, true) === false
        Test.@test Mycelia._gqc_effective_canonical(:canonical, false) === true
        Test.@test Mycelia._gqc_effective_canonical(:double, true) === true
        Test.@test Mycelia._gqc_effective_canonical(:double, false) === false
        ## :single vs :canonical changes repeat counts on RC-paired kmers
        seqs = ["ACGTAA"]
        ss = Mycelia.gqc_repeat_counts(seqs; k=3, canonical=Mycelia._gqc_effective_canonical(:single, true))
        can = Mycelia.gqc_repeat_counts(seqs; k=3, canonical=Mycelia._gqc_effective_canonical(:canonical, false))
        Test.@test ss.repeats != can.repeats || ss.distinct_kmers != can.distinct_kmers
    end

    Test.@testset "I1 params_hash invariant to ks/bh_ks order" begin
        a = Dict{String, Any}("ks" => [5, 11, 21], "bh_ks" => [11, 5], "n_shuffles" => 10)
        b = Dict{String, Any}("ks" => [21, 5, 11], "bh_ks" => [5, 11], "n_shuffles" => 10)
        Test.@test Mycelia._gqc_params_hash(a) == Mycelia._gqc_params_hash(b)
    end

    Test.@testset "resolve_gqc_input_dirs single vs graph-set" begin
        fixture = joinpath(@__DIR__, "..", "fixtures", "gqc_synthetic")
        single = Mycelia.resolve_gqc_input_dirs(fixture)
        Test.@test single.mode === :single_graph
        Test.@test length(single.dirs) == 1
        Test.@test single.dirs[1] == abspath(fixture)

        set_root = mktempdir()
        for (name, label) in (("graph_a", "set_a"), ("graph_b", "set_b"))
            dest = joinpath(set_root, name)
            mkpath(dest)
            cp(joinpath(fixture, "graph.gfa"), joinpath(dest, "graph.gfa"))
            meta = JSON.parsefile(joinpath(fixture, "metadata.json"))
            meta["graph_id"] = name
            meta["label"] = label
            open(joinpath(dest, "metadata.json"), "w") do io
                JSON.print(io, meta, 2)
            end
        end
        mkpath(joinpath(set_root, "notes"))
        write(joinpath(set_root, "notes", "readme.txt"), "skip me")

        resolved = Mycelia.resolve_gqc_input_dirs(set_root)
        Test.@test resolved.mode === :graph_set
        Test.@test length(resolved.dirs) == 2
        Test.@test basename.(resolved.dirs) == ["graph_a", "graph_b"]

        hybrid = mktempdir()
        cp(joinpath(fixture, "graph.gfa"), joinpath(hybrid, "graph.gfa"))
        cp(joinpath(fixture, "metadata.json"), joinpath(hybrid, "metadata.json"))
        child = joinpath(hybrid, "nested")
        mkpath(child)
        cp(joinpath(fixture, "graph.gfa"), joinpath(child, "graph.gfa"))
        hybrid_res = Mycelia.resolve_gqc_input_dirs(hybrid)
        Test.@test hybrid_res.mode === :single_graph
        Test.@test length(hybrid_res.dirs) == 1

        empty_dir = mktempdir()
        err = try
            Mycelia.resolve_gqc_input_dirs(empty_dir)
            nothing
        catch e
            e
        end
        Test.@test err isa ErrorException
        Test.@test occursin("single-graph folder", sprint(showerror, err))
    end

    Test.@testset "synthetic single-strand GFA packet" begin
        ENV["GKSwstype"] = "100"
        fixture = joinpath(@__DIR__, "..", "fixtures", "gqc_synthetic")
        out_root = mktempdir()
        packet = Mycelia.run_gqc(fixture, out_root;
            ks=[5],
            n_shuffles=10,
            shuffle_seed=20261021,
            strand_mode=:single,
            make_plots=true,
            write_node_link_json=true,
            null_kinds=(:dinucleotide,),
        )
        ## Packet path uses graph_id, not label
        Test.@test basename(packet) == "synth_ss_k5"
        Test.@test isdir(packet)
        Test.@test isfile(joinpath(packet, "manifest.json"))
        Test.@test isfile(joinpath(packet, "summary.csv"))
        Test.@test isfile(joinpath(packet, "summary.json"))
        Test.@test isfile(joinpath(packet, "report.md"))
        Test.@test isfile(joinpath(packet, "CHECKSUMS.sha256"))
        Test.@test isfile(joinpath(packet, "formats", "graph.graphml"))
        Test.@test isfile(joinpath(packet, "formats", "nodes.tsv"))
        Test.@test isfile(joinpath(packet, "formats", "edges.tsv"))
        Test.@test isfile(joinpath(packet, "formats", "adjacency.mtx"))
        Test.@test isfile(joinpath(packet, "formats", "node_link.json"))
        ## I7 PNG existence
        Test.@test isfile(joinpath(packet, "plots", "degree", "synth_ss_k5__k5.png")) ||
                   isfile(joinpath(packet, "plots", "structure", "synth_ss_k5__k5.png"))
        ## I2 frozen summary header
        header = first(readlines(joinpath(packet, "summary.csv")))
        Test.@test header == join(Mycelia.GQC_SUMMARY_COLUMNS, ',')
        ## I5 CHECKSUMS sorted by path
        checksum_paths = [split(l, r"\s+"; limit=2)[2] for l in readlines(joinpath(packet, "CHECKSUMS.sha256")) if !isempty(l)]
        Test.@test checksum_paths == sort(checksum_paths)
        ## MTX banner
        mtx = readlines(joinpath(packet, "formats", "adjacency.mtx"))
        Test.@test startswith(mtx[1], "%%MatrixMarket matrix coordinate pattern symmetric")
        nodes = readlines(joinpath(packet, "formats", "nodes.tsv"))
        n_nodes = length(nodes) - 1
        mtx_data = filter(l -> !startswith(l, '%'), mtx)
        dims = split(first(mtx_data))
        Test.@test parse(Int, dims[1]) == n_nodes
        Test.@test parse(Int, dims[2]) == n_nodes
        ## lower triangle: row >= col
        for line in mtx_data[2:end]
            r, c = parse.(Int, split(line))
            Test.@test r >= c
        end
        manifest = JSON.parsefile(joinpath(packet, "manifest.json"))
        Test.@test manifest["gqc_version"] == "0.1.0"
        Test.@test haskey(manifest, "input_sha256")
        Test.@test haskey(manifest, "params_hash")
        Test.@test manifest["params"]["bh_ks"] == [11]
        Test.@test manifest["params"]["canonical"] == false  ## :single
        Test.@test any(endswith(f, "graph.graphml") for f in manifest["formats"])
        g = Mycelia.read_gfa_records(joinpath(fixture, "graph.gfa"))
        Test.@test length(g.segments) == 5
        Test.@test length(g.links) == 4
        Test.@test g.segments["1"].depth == 10.0
    end

    Test.@testset "C3 MTX golden with reverse-only link" begin
        fixture = joinpath(@__DIR__, "..", "fixtures", "gqc_reverse_link")
        packet = Mycelia.run_gqc(fixture, mktempdir();
            ks=[5], n_shuffles=5, strand_mode=:single, make_plots=false,
            write_node_link_json=false, null_kinds=(:dinucleotide,))
        mtx = readlines(joinpath(packet, "formats", "adjacency.mtx"))
        Test.@test startswith(mtx[1], "%%MatrixMarket matrix coordinate pattern symmetric")
        ## nodes order 1..5; reverse-only L 5+ 2+ → undirected edge (5,2) → lower (5,2)
        data = filter(l -> !startswith(l, '%'), mtx)
        pairs = Set(Tuple(parse.(Int, split(l))) for l in data[2:end])
        Test.@test (5, 2) in pairs
        Test.@test !((2, 5) in pairs)  ## only lower triangle stored
        Test.@test (2, 1) in pairs
        Test.@test (3, 2) in pairs
        Test.@test (4, 3) in pairs
    end

    Test.@testset "C5 bh_ks defaults to [11] within-graph" begin
        fixture = joinpath(@__DIR__, "..", "fixtures", "gqc_synthetic")
        packet = Mycelia.run_gqc(fixture, mktempdir();
            ks=[5, 11, 21],
            n_shuffles=10,
            strand_mode=:single,
            make_plots=false,
            write_node_link_json=false,
            null_kinds=(:dinucleotide,),
        )
        rows = JSON.parsefile(joinpath(packet, "summary.json"))
        by_k = Dict(Int(r["k"]) => r for r in rows)
        ## JSON encodes missing as null → nothing
        Test.@test by_k[5]["dinucleotide_bh_q"] === nothing
        Test.@test by_k[21]["dinucleotide_bh_q"] === nothing
        Test.@test by_k[11]["dinucleotide_bh_q"] isa Real
        Test.@test by_k[11]["dinucleotide_bh_q"] ≈ by_k[11]["dinucleotide_shuffle_p_upper"]
        man = JSON.parsefile(joinpath(packet, "manifest.json"))
        Test.@test man["params"]["bh_ks"] == [11]
        Test.@test occursin("within-graph", man["bh_note"])
    end

    Test.@testset "graph-set folder runs all child graphs" begin
        fixture = joinpath(@__DIR__, "..", "fixtures", "gqc_synthetic")
        set_root = mktempdir()
        for (name, label) in (("graph_a", "same_label"), ("graph_b", "same_label"))
            dest = joinpath(set_root, name)
            mkpath(dest)
            cp(joinpath(fixture, "graph.gfa"), joinpath(dest, "graph.gfa"))
            meta = JSON.parsefile(joinpath(fixture, "metadata.json"))
            meta["graph_id"] = name
            meta["label"] = label
            open(joinpath(dest, "metadata.json"), "w") do io
                JSON.print(io, meta, 2)
            end
        end
        out_root = mktempdir()
        returned = Mycelia.run_gqc(set_root, out_root;
            ks=[5],
            n_shuffles=5,
            shuffle_seed=20261021,
            genus_index=10,
            strand_mode=:single,
            make_plots=false,
            write_node_link_json=false,
            null_kinds=(:dinucleotide,),
        )
        Test.@test returned == abspath(out_root)
        ## Unique by graph_id even when labels collide
        packet_a = joinpath(out_root, "graph_a")
        packet_b = joinpath(out_root, "graph_b")
        Test.@test isfile(joinpath(packet_a, "manifest.json"))
        Test.@test isfile(joinpath(packet_b, "manifest.json"))
        man_a = JSON.parsefile(joinpath(packet_a, "manifest.json"))
        man_b = JSON.parsefile(joinpath(packet_b, "manifest.json"))
        Test.@test man_a["params"]["genus_index"] == 10
        Test.@test man_b["params"]["genus_index"] == 11
        Test.@test man_a["graph_id"] == "graph_a"
        Test.@test man_b["graph_id"] == "graph_b"
        Test.@test man_a["graph_set"] == "same_label"
        Test.@test man_b["graph_set"] == "same_label"
    end

    Test.@testset "C4 refuse overwrite on duplicate graph_id" begin
        fixture = joinpath(@__DIR__, "..", "fixtures", "gqc_synthetic")
        out_root = mktempdir()
        Mycelia.run_gqc(fixture, out_root;
            ks=[5], n_shuffles=5, make_plots=false, write_node_link_json=false,
            null_kinds=(:dinucleotide,), strand_mode=:single)
        err = try
            Mycelia.run_gqc(fixture, out_root;
                ks=[5], n_shuffles=5, make_plots=false, write_node_link_json=false,
                null_kinds=(:dinucleotide,), strand_mode=:single)
            nothing
        catch e
            e
        end
        Test.@test err isa ErrorException
        Test.@test occursin("refusing to overwrite", sprint(showerror, err))
    end

    Test.@testset "I6 gzipped reference decode (per-record, no chimeric join)" begin
        fixture = joinpath(@__DIR__, "..", "fixtures", "gqc_synthetic")
        tmp = mktempdir()
        cp(joinpath(fixture, "graph.gfa"), joinpath(tmp, "graph.gfa"))
        cp(joinpath(fixture, "metadata.json"), joinpath(tmp, "metadata.json"))
        refdir = mkpath(joinpath(tmp, "reference"))
        open(CodecZlib.GzipCompressorStream, joinpath(refdir, "ref.fasta.gz"), "w") do io
            println(io, ">r1")
            println(io, "ATGCGATCGATCGTTACCGGATCCATGCGATTACGGCATGCGATCGATCGTTACCGGATCCAA")
            println(io, ">r2")
            println(io, "GCTTAGCTAGGCTTAACCTAGGCTTAGCTAA")
        end
        seqs = Mycelia._gqc_load_reference_sequences(refdir)
        Test.@test length(seqs) == 2
        Test.@test startswith(seqs[1], "ATGCGATCGA")
        Test.@test startswith(seqs[2], "GCTTAGCTAG")
        ## Must not concatenate across records
        Test.@test length(seqs[1]) == 63
        Test.@test length(seqs[2]) == 31
    end

    Test.@testset "null_source=:auto uses reference on k-mer-node graphs" begin
        k = 5
        ref = "ATGCGATCGATCGTTACCGGATCCATGCGATTACGGCATGCGATCGA"
        n_shuffles = 20
        seed = 20261021
        ## Canonical k-mer-node GFA: one segment per distinct canonical k-mer (len=k)
        kmers = String[ref[i:(i + k - 1)] for i in 1:(length(ref) - k + 1)]
        ckmers = String[Mycelia._gqc_canonical_kmer(km) for km in kmers]
        uniq = unique(ckmers)
        id_of = Dict(km => string(i) for (i, km) in enumerate(uniq))
        root = mktempdir()
        open(joinpath(root, "graph.gfa"), "w") do io
            println(io, "H\tVN:Z:1.0")
            for (i, km) in enumerate(uniq)
                println(io, "S\t", i, '\t', km, "\tLN:i:", k)
            end
            for i in 1:(length(ckmers) - 1)
                a = id_of[ckmers[i]]; b = id_of[ckmers[i + 1]]
                a == b && continue
                println(io, "L\t", a, "\t+\t", b, "\t+\t", k - 1, "M")
            end
        end
        open(joinpath(root, "metadata.json"), "w") do io
            JSON.print(io, Dict(
                "graph_id" => "kmer_node_k5",
                "label" => "kmer_node",
                "k" => k,
                "strand_mode" => "canonical",
                "circular" => false,
            ), 2)
        end
        mkpath(joinpath(root, "reference"))
        open(joinpath(root, "reference", "genome.fasta"), "w") do io
            println(io, ">genome")
            println(io, ref)
        end

        expected = Mycelia.gqc_repeat_null_stats([ref];
            k=k, n_shuffles=n_shuffles, shuffle_seed=seed,
            circular=false, canonical=true, null_kind=:dinucleotide)
        ## Segment-only nulls: each node is one distinct k-mer → observed repeats 0
        seg_only = Mycelia.gqc_repeat_null_stats(uniq;
            k=k, n_shuffles=n_shuffles, shuffle_seed=seed,
            circular=false, canonical=true, null_kind=:dinucleotide)
        Test.@test seg_only.observed_repeats == 0
        Test.@test expected.observed_repeats > 0

        packet = Mycelia.run_gqc(root, mktempdir();
            ks=[k],
            n_shuffles=n_shuffles,
            shuffle_seed=seed,
            strand_mode=:canonical,
            null_source=:auto,
            make_plots=false,
            write_node_link_json=false,
            null_kinds=(:dinucleotide,),
        )
        rows = JSON.parsefile(joinpath(packet, "summary.json"))
        row = only(rows)
        Test.@test row["null_source"] == "reference"
        Test.@test row["null_applicable"] == true
        Test.@test row["dinucleotide_observed_repeats"] == expected.observed_repeats
        Test.@test row["dinucleotide_null_mean"] ≈ expected.null_mean atol=1e-12
        Test.@test row["dinucleotide_shuffle_p_upper"] ≈ expected.shuffle_p_upper atol=1e-12
        man = JSON.parsefile(joinpath(packet, "manifest.json"))
        Test.@test man["params"]["null_source"] == "auto"
        Test.@test man["params"]["resolved_null_source"] == "reference"
        Test.@test man["params"]["null_applicable"] == true

        ## Explicit segments on k-mer-node graph → not_applicable (no degenerate p)
        packet2 = Mycelia.run_gqc(root, mktempdir();
            ks=[k], n_shuffles=n_shuffles, shuffle_seed=seed,
            strand_mode=:canonical, null_source=:segments,
            make_plots=false, write_node_link_json=false,
            null_kinds=(:dinucleotide,),
        )
        row2 = only(JSON.parsefile(joinpath(packet2, "summary.json")))
        Test.@test row2["null_source"] == "segments"
        Test.@test row2["null_applicable"] == false
        Test.@test row2["dinucleotide_shuffle_p_upper"] === nothing
        Test.@test row2["dinucleotide_observed_repeats"] === nothing
    end
end
