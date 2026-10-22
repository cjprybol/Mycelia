# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/gfa_canonical_orientation_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/gfa_canonical_orientation_test.jl", "test/4_assembly", execute=false)'
# ```

import Test
import Mycelia
import MetaGraphsNext
import BioSequences
import FASTX
import Kmers
import Graphs

Test.@testset "write_gfa_canonical link orientation" begin
    function parse_gfa(path)
        segment_lines = String[]
        segments = Dict{String, String}()
        links = Tuple{String, Char, String, Char, Int}[]
        for line in eachline(path)
            fields = split(line, '\t')
            if fields[1] == "S"
                push!(segment_lines, line)
                segments[fields[2]] = fields[3]
            elseif fields[1] == "L"
                overlap = parse(Int, rstrip(fields[6], 'M'))
                push!(links, (
                    fields[2], only(fields[3]), fields[4], only(fields[5]), overlap))
            end
        end
        return segment_lines, segments, links
    end

    function oriented(sequence::AbstractString, orientation::Char)
        orientation == '+' && return String(sequence)
        return string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(sequence)))
    end

    function link_is_spelled(segments, link)
        (a, a_or, b, b_or, overlap) = link
        a_seq = oriented(segments[a], a_or)
        b_seq = oriented(segments[b], b_or)
        return a_seq[(end - overlap + 1):end] == b_seq[1:overlap]
    end

    flip(c::Char) = c == '+' ? '-' : '+'

    # Return the exception thrown by `f`, or `nothing`, so tests can check the
    # message as well as the type.
    function thrown(f)
        try
            f()
        catch e
            return e
        end
        return nothing
    end

    # GFA-equivalence-aware membership: (a, oa, b, ob) == (b, !ob, a, !oa).
    function has_link(links, a, a_or, b, b_or)
        return any(links) do (x, x_or, y, y_or, _)
            (x, x_or, y, y_or) == (a, a_or, b, b_or) ||
                (x, x_or, y, y_or) == (b, flip(b_or), a, flip(a_or))
        end
    end

    function canonical_graph(sequence, k)
        record = FASTX.FASTA.Record("read_1", sequence)
        return Mycelia.Rhizomorph.build_kmer_graph_canonical([record], k)
    end

    function segment_id_of(graph, kmer)
        return string(findfirst(==(kmer), collect(MetaGraphsNext.labels(graph))))
    end

    Test.@testset "reverse-complement junctions get + - / - + links" begin
        # k = 5. Canonical forms of the successive read k-mers:
        #   ACCCA(+) CCCAG(+) CCAGG(+) ACCTG(-) AACCT(-) CAACC(-) GTTGA(+)
        # so CCAGG -> CAGGT is a `+ -` junction and CAACC -> GTTGA is `- +`.
        k = 5
        sequence = "ACCCAGGTTGA"
        graph = canonical_graph(sequence, k)

        mktempdir() do dir
            gfa = joinpath(dir, "canonical.gfa")
            Test.@test Mycelia.Rhizomorph.write_gfa_canonical(graph, gfa) == gfa
            _, segments, links = parse_gfa(gfa)

            Test.@test !isempty(links)
            for link in links
                Test.@test link[5] == k - 1
                Test.@test link_is_spelled(segments, link)
            end

            # Every adjacency in the read is represented with its observed orientation.
            read_kmers = [Kmers.DNAKmer{k}(sequence[i:(i + k - 1)])
                          for i in 1:(length(sequence) - k + 1)]
            for (x, y) in zip(read_kmers[1:(end - 1)], read_kmers[2:end])
                cx = BioSequences.canonical(x)
                cy = BioSequences.canonical(y)
                x_or = x == cx ? '+' : '-'
                y_or = y == cy ? '+' : '-'
                Test.@test has_link(links, segment_id_of(graph, cx), x_or,
                    segment_id_of(graph, cy), y_or)
            end

            ccagg = segment_id_of(graph, Kmers.DNAKmer{k}("CCAGG"))
            acctg = segment_id_of(graph, Kmers.DNAKmer{k}("ACCTG"))
            caacc = segment_id_of(graph, Kmers.DNAKmer{k}("CAACC"))
            gttga = segment_id_of(graph, Kmers.DNAKmer{k}("GTTGA"))
            Test.@test has_link(links, ccagg, '+', acctg, '-')
            Test.@test has_link(links, caacc, '-', gttga, '+')
            Test.@test !has_link(links, ccagg, '+', acctg, '+')
            Test.@test any(l -> l[2] != l[4], links)

            # No GFA-equivalent duplicates.
            link_keys = map(links) do (a, a_or, b, b_or, _)
                min((a, a_or, b, b_or), (b, flip(b_or), a, flip(a_or)))
            end
            Test.@test length(unique(link_keys)) == length(link_keys)

            # Pin the exact output, including which member of each GFA-equivalent
            # pair is written (more `+` first, then lower segment index).
            Test.@test links == [
                ("1", '+', "2", '+', 4), ("2", '+', "3", '+', 4),
                ("3", '+', "4", '-', 4), ("5", '+', "4", '+', 4),
                ("6", '+', "5", '+', 4), ("6", '-', "7", '+', 4)]

            # S lines match write_gfa_next exactly, so segment IDs join.
            legacy = joinpath(dir, "legacy.gfa")
            Mycelia.Rhizomorph.write_gfa_next(graph, legacy)
            Test.@test parse_gfa(legacy)[1] == parse_gfa(gfa)[1]
        end
    end

    Test.@testset "self-loops" begin
        k = 5
        mktempdir() do dir
            # ACGCG -> CGCGT == rc(ACGCG): self-loop written `+ -`.
            graph = canonical_graph("ACGCGT", k)
            gfa = joinpath(dir, "self_rc.gfa")
            Mycelia.Rhizomorph.write_gfa_canonical(graph, gfa)
            _, segments, links = parse_gfa(gfa)
            id = segment_id_of(graph, Kmers.DNAKmer{k}("ACGCG"))
            Test.@test has_link(links, id, '+', id, '-')
            Test.@test all(link -> link_is_spelled(segments, link), links)

            # Homopolymer: AAAAA -> AAAAA self-loop written once as `+ +`.
            graph = canonical_graph("AAAAAAA", k)
            gfa = joinpath(dir, "self_fwd.gfa")
            Mycelia.Rhizomorph.write_gfa_canonical(graph, gfa)
            _, segments, links = parse_gfa(gfa)
            Test.@test links == [("1", '+', "1", '+', k - 1)]
        end
    end

    Test.@testset "forward-only cases are unaffected" begin
        k = 5
        mktempdir() do dir
            # All k-mers of AAAACCCC are already canonical, so the canonical
            # writer must agree with write_gfa_next (all `+ +`).
            graph = canonical_graph("AAAACCCC", k)
            Test.@test all(label -> label == BioSequences.canonical(label),
                MetaGraphsNext.labels(graph))
            canonical_gfa = joinpath(dir, "canonical.gfa")
            legacy_gfa = joinpath(dir, "legacy.gfa")
            Mycelia.Rhizomorph.write_gfa_canonical(graph, canonical_gfa)
            Mycelia.Rhizomorph.write_gfa_next(graph, legacy_gfa)
            _, _, canonical_links = parse_gfa(canonical_gfa)
            _, _, legacy_links = parse_gfa(legacy_gfa)
            Test.@test Set(canonical_links) == Set(legacy_links)
            Test.@test all(l -> l[2] == '+' && l[4] == '+', canonical_links)

            # Single-strand (directed) graphs keep using write_gfa_next and
            # are rejected by the canonical writer.
            record = FASTX.FASTA.Record("read_1", "ACCCAGGTTGA")
            ss_graph = Mycelia.Rhizomorph.build_kmer_graph([record], k;
                dataset_id = "dataset_01", mode = :singlestrand)
            ss_gfa = joinpath(dir, "singlestrand.gfa")
            Mycelia.Rhizomorph.write_gfa_next(ss_graph, ss_gfa)
            _, ss_segments, ss_links = parse_gfa(ss_gfa)
            Test.@test length(ss_links) == length("ACCCAGGTTGA") - k
            Test.@test all(l -> l[2] == '+' && l[4] == '+', ss_links)
            Test.@test all(link -> link_is_spelled(ss_segments, link), ss_links)
            err = thrown(() -> Mycelia.Rhizomorph.write_gfa_canonical(
                ss_graph, joinpath(dir, "rejected.gfa")))
            Test.@test err isa ArgumentError
            Test.@test occursin("expects an undirected canonical graph",
                sprint(showerror, err))
        end
    end

    Test.@testset "argument validation" begin
        graph = canonical_graph("ACCCAGGTTGA", 5)
        mktempdir() do dir
            err = thrown(() -> Mycelia.Rhizomorph.write_gfa_canonical(
                graph, joinpath(dir, "x.gfa"); segment_naming = :bogus))
            Test.@test err isa ArgumentError
            Test.@test occursin("segment_naming must be", sprint(showerror, err))
            gfa = joinpath(dir, "named.gfa")
            Mycelia.Rhizomorph.write_gfa_canonical(graph, gfa; segment_naming = :sequence)
            _, segments, links = parse_gfa(gfa)
            Test.@test all(((id, seq),) -> id == seq, segments)
            Test.@test all(link -> link_is_spelled(segments, link), links)
        end
    end

    Test.@testset "palindromic k-mer at even k" begin
        # k = 4: ACATG has k-mers ACAT and CATG; CATG is its own reverse
        # complement, so `ACAT + CATG +` and `ACAT + CATG -` are both spelled
        # and both written for the single canonical edge.
        k = 4
        graph = canonical_graph("ACATG", k)
        mktempdir() do dir
            gfa = joinpath(dir, "palindrome.gfa")
            Mycelia.Rhizomorph.write_gfa_canonical(graph, gfa)
            _, segments, links = parse_gfa(gfa)
            acat = segment_id_of(graph, Kmers.DNAKmer{k}("ACAT"))
            catg = segment_id_of(graph, Kmers.DNAKmer{k}("CATG"))
            Test.@test length(links) == 2
            Test.@test has_link(links, acat, '+', catg, '+')
            Test.@test has_link(links, acat, '+', catg, '-')
            Test.@test all(link -> link_is_spelled(segments, link), links)
        end
    end

    Test.@testset "input contract on hand-built graphs" begin
        mktempdir() do dir
            # An edge with no (k-1) overlap cannot come from the canonical
            # builders, but the public signature accepts it; it must raise
            # rather than be dropped silently.
            g = MetaGraphsNext.MetaGraph(Graphs.Graph();
                label_type = Kmers.DNAKmer{5}, vertex_data_type = Nothing,
                edge_data_type = Nothing)
            g[Kmers.DNAKmer{5}("AAAAA")] = nothing
            g[Kmers.DNAKmer{5}("CCGCC")] = nothing
            g[Kmers.DNAKmer{5}("AAAAA"), Kmers.DNAKmer{5}("CCGCC")] = nothing
            err = thrown(() -> Mycelia.Rhizomorph.write_gfa_canonical(
                g, joinpath(dir, "no_overlap.gfa")))
            Test.@test err isa ErrorException
            Test.@test occursin("no valid 4-base overlap", sprint(showerror, err))

            h = MetaGraphsNext.MetaGraph(Graphs.Graph();
                label_type = String, vertex_data_type = Nothing,
                edge_data_type = Nothing)
            h["ACGTA"] = nothing
            err = thrown(() -> Mycelia.Rhizomorph.write_gfa_canonical(
                h, joinpath(dir, "string_labels.gfa")))
            Test.@test err isa ArgumentError
            Test.@test occursin("requires nucleotide k-mer labels",
                sprint(showerror, err))
        end
    end
end
