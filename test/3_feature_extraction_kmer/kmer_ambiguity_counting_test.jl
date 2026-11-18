import Test
import BioSequences
import FASTX
import Kmers
import Mycelia

if !isdefined(Main, :test_throws_message)
    include(joinpath(dirname(@__DIR__), "test_helpers.jl"))
end

function _assert_argument_error(f, fragments...)
    try
        f()
        Test.@test false
    catch err
        Test.@test err isa ArgumentError
        message = sprint(showerror, err)
        for fragment in fragments
            Test.@test occursin(fragment, message)
        end
    end
    return nothing
end

function _dna_count(sequence::AbstractString; kwargs...)
    return Mycelia.count_kmers(
        Kmers.DNAKmer{2},
        BioSequences.LongDNA{4}(sequence);
        kwargs...,
    )
end

Test.@testset "k-mer ambiguity modes" begin
    Test.@testset "default matches UnambiguousDNAMers skip" begin
        clean = BioSequences.LongDNA{4}("ACGTAC")
        default_counts = Mycelia.count_kmers(Kmers.DNAKmer{3}, clean)
        skip_counts = Mycelia.count_kmers(Kmers.DNAKmer{3}, clean; ambiguity = :skip)
        expand_counts = Mycelia.count_kmers(Kmers.DNAKmer{3}, clean; ambiguity = :expand)
        error_counts = Mycelia.count_kmers(Kmers.DNAKmer{3}, clean; ambiguity = :error)
        Test.@test default_counts == skip_counts
        Test.@test expand_counts == default_counts
        Test.@test error_counts == default_counts
        Test.@test sum(values(default_counts)) == length(clean) - 3 + 1
        Test.@test default_counts[Kmers.DNAKmer{3}("ACG")] == 1
        Test.@test default_counts[Kmers.DNAKmer{3}("CGT")] == 1
        Test.@test default_counts[Kmers.DNAKmer{3}("GTA")] == 1
        Test.@test default_counts[Kmers.DNAKmer{3}("TAC")] == 1

        with_n = BioSequences.LongDNA{4}("ACNGT")
        skipped = Mycelia.count_kmers(Kmers.DNAKmer{2}, with_n)
        explicit_skip = Mycelia.count_kmers(Kmers.DNAKmer{2}, with_n; ambiguity = :skip)
        Test.@test skipped == explicit_skip
        Test.@test skipped[Kmers.DNAKmer{2}("AC")] == 1
        Test.@test skipped[Kmers.DNAKmer{2}("GT")] == 1
        Test.@test length(skipped) == 2
    end

    Test.@testset "expand enumerates N and other IUPAC codes" begin
        expanded = _dna_count("ACNGT"; ambiguity = :expand)
        # Windows: AC, CN, NG, GT. NG resolves to AG/CG/GG/TG, so CG is counted twice.
        Test.@test expanded[Kmers.DNAKmer{2}("AC")] == 1
        Test.@test expanded[Kmers.DNAKmer{2}("CA")] == 1
        Test.@test expanded[Kmers.DNAKmer{2}("CC")] == 1
        Test.@test expanded[Kmers.DNAKmer{2}("CG")] == 2
        Test.@test expanded[Kmers.DNAKmer{2}("CT")] == 1
        Test.@test expanded[Kmers.DNAKmer{2}("AG")] == 1
        Test.@test expanded[Kmers.DNAKmer{2}("GG")] == 1
        Test.@test expanded[Kmers.DNAKmer{2}("TG")] == 1
        Test.@test expanded[Kmers.DNAKmer{2}("GT")] == 1
        Test.@test length(expanded) == 9
        Test.@test sum(values(expanded)) == 1 + 4 + 4 + 1

        purine = _dna_count("AR"; ambiguity = :expand)
        Test.@test purine[Kmers.DNAKmer{2}("AA")] == 1
        Test.@test purine[Kmers.DNAKmer{2}("AG")] == 1
        Test.@test length(purine) == 2
        Test.@test isempty(_dna_count("AR"; ambiguity = :skip))

        single_n = Mycelia.count_kmers(
            Kmers.DNAKmer{1},
            BioSequences.LongDNA{4}("NN");
            ambiguity = :expand,
        )
        for base in ("A", "C", "G", "T")
            Test.@test single_n[Kmers.DNAKmer{1}(base)] == 2
        end
        Test.@test length(single_n) == 4
    end

    Test.@testset "error rejects N and gaps, expand rejects gaps" begin
        _assert_argument_error(
            () -> _dna_count("ACN"; ambiguity = :error),
            "ambiguous base",
            "position 3",
        )
        _assert_argument_error(
            () -> _dna_count("AC-G"; ambiguity = :error),
            "uncertain base",
            "position 3",
        )
        _assert_argument_error(
            () -> _dna_count("AC-G"; ambiguity = :expand),
            "cannot expand",
            "position 3",
        )
        skipped_gap = _dna_count("AC-G")
        Test.@test skipped_gap[Kmers.DNAKmer{2}("AC")] == 1
        Test.@test length(skipped_gap) == 1

        short_n = Mycelia.count_kmers(
            Kmers.DNAKmer{3},
            BioSequences.LongDNA{4}("N");
            ambiguity = :error,
        )
        Test.@test isempty(short_n)
        Test.@test isempty(Mycelia.count_kmers(
            Kmers.DNAKmer{3},
            BioSequences.LongDNA{4}("N");
            ambiguity = :expand,
        ))

        _assert_argument_error(
            () -> _dna_count("ACGT"; ambiguity = :keep),
            "ambiguity must be :skip, :expand, or :error",
        )
    end

    Test.@testset "RNA N uses the same three modes" begin
        rna = BioSequences.LongRNA{4}("ACGN")
        default_counts = Mycelia.count_kmers(Kmers.RNAKmer{2}, rna)
        skip_counts = Mycelia.count_kmers(Kmers.RNAKmer{2}, rna; ambiguity = :skip)
        Test.@test default_counts == skip_counts
        Test.@test skip_counts[Kmers.RNAKmer{2}("AC")] == 1
        Test.@test skip_counts[Kmers.RNAKmer{2}("CG")] == 1
        Test.@test length(skip_counts) == 2

        expanded = Mycelia.count_kmers(Kmers.RNAKmer{2}, rna; ambiguity = :expand)
        Test.@test expanded[Kmers.RNAKmer{2}("AC")] == 1
        Test.@test expanded[Kmers.RNAKmer{2}("CG")] == 1
        Test.@test expanded[Kmers.RNAKmer{2}("GA")] == 1
        Test.@test expanded[Kmers.RNAKmer{2}("GC")] == 1
        Test.@test expanded[Kmers.RNAKmer{2}("GG")] == 1
        Test.@test expanded[Kmers.RNAKmer{2}("GU")] == 1
        Test.@test sum(values(expanded)) == 6

        clean = BioSequences.LongRNA{4}("ACGU")
        clean_counts = Mycelia.count_kmers(Kmers.RNAKmer{2}, clean)
        Test.@test Mycelia.count_kmers(Kmers.RNAKmer{2}, clean; ambiguity = :expand) ==
                   clean_counts
        Test.@test Mycelia.count_kmers(Kmers.RNAKmer{2}, clean; ambiguity = :error) ==
                   clean_counts

        _assert_argument_error(
            () -> Mycelia.count_kmers(Kmers.RNAKmer{2}, rna; ambiguity = :error),
            "ambiguous base",
            "position 4",
        )
    end

    Test.@testset "records, files, and canonical counts forward the mode" begin
        record = FASTX.FASTA.Record("with_n", "ACNGT")
        from_record = Mycelia.count_kmers(Kmers.DNAKmer{2}, record; ambiguity = :expand)
        from_sequence = _dna_count("ACNGT"; ambiguity = :expand)
        Test.@test from_record == from_sequence

        records = [
            FASTX.FASTA.Record("a", "ACGT"),
            FASTX.FASTA.Record("b", "ACNG"),
        ]
        merged = Mycelia.count_kmers(Kmers.DNAKmer{2}, records; ambiguity = :skip)
        Test.@test merged[Kmers.DNAKmer{2}("AC")] == 2
        Test.@test merged[Kmers.DNAKmer{2}("CG")] == 1
        Test.@test merged[Kmers.DNAKmer{2}("GT")] == 1
        Test.@test length(merged) == 3
        _assert_argument_error(
            () -> Mycelia.count_kmers(Kmers.DNAKmer{2}, records; ambiguity = :error),
            "ambiguous base",
        )

        mktempdir() do dir
            fasta_path = joinpath(dir, "with_n.fasta")
            open(fasta_path, "w") do io
                FASTX.write(io, FASTX.FASTA.Record("with_n", "ACNG"))
            end
            from_file = Mycelia.count_kmers(
                Kmers.DNAKmer{2},
                fasta_path;
                ambiguity = :skip,
            )
            Test.@test from_file == _dna_count("ACNG"; ambiguity = :skip)
            open(fasta_path, "r") do io
                from_reader = Mycelia.count_kmers(
                    Kmers.DNAKmer{2},
                    FASTX.FASTA.Reader(io);
                    ambiguity = :expand,
                )
                Test.@test from_reader == _dna_count("ACNG"; ambiguity = :expand)
            end
        end

        raw = Mycelia.count_kmers(
            Kmers.DNAKmer{2},
            BioSequences.LongDNA{4}("ACN");
            ambiguity = :expand,
        )
        canonical = Mycelia.count_canonical_kmers(
            Kmers.DNAKmer{2},
            BioSequences.LongDNA{4}("ACN");
            ambiguity = :expand,
        )
        Test.@test canonical == Mycelia.canonicalize_kmer_counts(raw)
        Test.@test sum(values(canonical)) == sum(values(raw))
        skipped_canonical = Mycelia.count_canonical_kmers(
            Kmers.DNAKmer{2},
            BioSequences.LongDNA{4}("ACNGT"),
        )
        Test.@test skipped_canonical == Mycelia.canonicalize_kmer_counts(_dna_count("ACNGT"))
    end

    Test.@testset "amino acid default is unchanged and other modes are rejected" begin
        aa_seq = BioSequences.LongAA("ACDEFG")
        default_counts = Mycelia.count_kmers(Kmers.AAKmer{2}, aa_seq)
        skip_counts = Mycelia.count_kmers(Kmers.AAKmer{2}, aa_seq; ambiguity = :skip)
        Test.@test default_counts == skip_counts
        Test.@test sum(values(default_counts)) == length(aa_seq) - 2 + 1
        _assert_argument_error(
            () -> Mycelia.count_kmers(Kmers.AAKmer{2}, aa_seq; ambiguity = :expand),
            "amino acid",
            "only ambiguity=:skip",
        )
        _assert_argument_error(
            () -> Mycelia.count_kmers(Kmers.AAKmer{2}, aa_seq; ambiguity = :error),
            "amino acid",
        )
    end
end
