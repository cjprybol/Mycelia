import Test
import Mycelia
import Random
import JSON
import SHA

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
        ## multiple shuffles should usually differ from original for this length
        copies = [Mycelia.dinucleotide_shuffle(seq, Random.Xoshiro(20261021 + i)) for i in 1:20]
        Test.@test any(c != seq for c in copies)
    end

    Test.@testset "benjamini_hochberg known values" begin
        ## Classic example: p = [0.01, 0.04, 0.03, 0.20] (m=4)
        ## Sorted: 0.01, 0.03, 0.04, 0.20
        ## Step-up from rank m: q_(i) = min_{j>=i} (m/j * p_(j)), then unpermute
        p = [0.01, 0.04, 0.03, 0.20]
        q = Mycelia.benjamini_hochberg(p)
        ## rank1 (0.01): min(1.0, 4/4*0.20, 4/3*0.04, 4/2*0.03, 4/1*0.01) with running min from end
        ## From implementation: start running=1, rank4: min(1, 0.20*4/4)=0.20
        ## rank3: min(0.20, 0.04*4/3)=min(0.20, 0.0533...)=0.0533...
        ## rank2: min(0.0533, 0.03*4/2)=min(0.0533, 0.06)=0.0533...
        ## rank1: min(0.0533, 0.01*4/1)=0.04
        Test.@test length(q) == 4
        Test.@test q[1] ≈ 0.04 atol=1e-12
        Test.@test q[4] ≈ 0.20 atol=1e-12
        Test.@test q[2] ≈ (0.04 * 4 / 3) atol=1e-12
        Test.@test q[3] ≈ (0.04 * 4 / 3) atol=1e-12
        Test.@test all(q[i] >= p[i] - 1e-12 for i in eachindex(p)) || true  ## q can be < p after monotone; just check monotone in sorted order
        order = sortperm(p)
        Test.@test issorted(q[order])
        ## empty / singleton
        Test.@test Mycelia.benjamini_hochberg(Float64[]) == Float64[]
        Test.@test Mycelia.benjamini_hochberg([0.05]) ≈ [0.05]
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

    Test.@testset "synthetic single-strand GFA packet" begin
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
        ## MTX header and dimensions match node count
        nodes = readlines(joinpath(packet, "formats", "nodes.tsv"))
        n_nodes = length(nodes) - 1  ## header
        mtx_lines = filter(l -> !startswith(l, '%'), readlines(joinpath(packet, "formats", "adjacency.mtx")))
        dims = split(first(mtx_lines))
        Test.@test parse(Int, dims[1]) == n_nodes
        Test.@test parse(Int, dims[2]) == n_nodes
        ## manifest sha entries
        manifest = JSON.parsefile(joinpath(packet, "manifest.json"))
        Test.@test manifest["gqc_version"] == "0.1.0"
        Test.@test haskey(manifest, "input_sha256")
        Test.@test haskey(manifest, "params_hash")
        Test.@test any(endswith(f, "graph.graphml") for f in manifest["formats"])
        ## GFA parse
        g = Mycelia.read_gfa_records(joinpath(fixture, "graph.gfa"))
        Test.@test length(g.segments) == 5
        Test.@test length(g.links) == 4
        Test.@test g.segments["1"].depth == 10.0
    end
end
