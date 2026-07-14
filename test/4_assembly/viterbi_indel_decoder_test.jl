# Indel-aware pair-HMM extension to the graph-Viterbi read corrector.
#
# Milestone A: the config carries indel/gap fields that DEFAULT to the
#   substitution-collapse values (indel_moves=false, all gap masses 0), so the
#   existing decoder is byte-identical and the whole existing suite passes with
#   no edits.
# Milestone B: an inserted-unit and a deleted-unit observation each correct back
#   to the reference under an indel-enabled config (beam_width=typemax), and the
#   substitution fixture stays byte-identical with indel_moves=false.
# Milestone C: bounding knobs (adaptive band, D_max, insertion-run cap) and a
#   small nanopore-style smoke that provably routes through the indel kernel.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/viterbi_indel_decoder_test.jl")'

import BioSequences
import FASTX
import Logging
import MetaGraphsNext
import Mycelia
import Test

function indel_decoded_label_strings(
        result::Mycelia.ViterbiCorrectionResult
)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

# A read-quality observation wrapper for the quality-aware oracle case (mirrors the
# struct used by viterbi_variable_length_graph_correction_test.jl).
struct IndelQualityObservation{S}
    sequence::S
    quality_scores::Vector{UInt8}
end

function indel_quality_record(
        identifier::AbstractString,
        sequence::AbstractString,
        phred::AbstractVector{<:Integer}
)::FASTX.FASTQ.Record
    quality = String([Char(score + 33) for score in phred])
    return FASTX.FASTQ.Record(identifier, sequence, quality)
end

# Bare-k-mer observation vector for a DNA string, threaded through the SAME
# decomposition the k-mer-graph corrector consumes (detect_alphabet ->
# typed-sequence -> record k-mer iterator). Used by the C1 base-level indel test.
function indel_string_kmer_units(seq_string::AbstractString, k::Int)
    record = FASTX.FASTA.Record("obs", seq_string)
    alphabet = Mycelia.detect_alphabet(seq_string)
    sequence_type = Mycelia.alphabet_to_biosequence_type(alphabet)
    sequence = Mycelia.extract_typed_sequence(record, sequence_type)
    return collect(Mycelia._record_kmer_iterator(sequence_type, k, sequence))
end

function indel_no_path_correction_kernel(
        ::MetaGraphsNext.MetaGraph,
        ::AbstractVector;
        config::Mycelia.ViterbiCorrectionConfig,
        weighted_graph::Any = nothing,
)::Mycelia.ViterbiCorrectionResult
    diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_indel_pair_hmm,
        :reason => :no_surviving_path,
    )
    path = Mycelia.Rhizomorph.ViterbiDecodingResult(
        nothing, -Inf, diagnostics)
    return Mycelia.ViterbiCorrectionResult(
        Any[nothing],
        Mycelia.Rhizomorph.ViterbiDecodingResult[path],
        Dict{Symbol, Any}(),
    )
end

function indel_empty_path_correction_kernel(
        ::MetaGraphsNext.MetaGraph,
        ::AbstractVector;
        config::Mycelia.ViterbiCorrectionConfig,
        weighted_graph::Any = nothing,
)::Mycelia.ViterbiCorrectionResult
    diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_indel_pair_hmm,
    )
    path = Mycelia.Rhizomorph.ViterbiDecodingResult(
        nothing, -Inf, diagnostics)
    return Mycelia.ViterbiCorrectionResult(
        Any[Any[]],
        Mycelia.Rhizomorph.ViterbiDecodingResult[path],
        Dict{Symbol, Any}(),
    )
end

function indel_partial_complete_correction_kernel(
        graph::MetaGraphsNext.MetaGraph,
        observations::AbstractVector;
        config::Mycelia.ViterbiCorrectionConfig,
        weighted_graph::Any = nothing,
)::Mycelia.ViterbiCorrectionResult
    @assert length(only(observations)) == 5
    diagnostics = Dict{Symbol, Any}(
        :algorithm => :viterbi_indel_pair_hmm,
        :truncated => false,
        :completed_columns => 1,
        :decoded_read_index => 1,
        :move_trace => Symbol[:M],
        :read_index_trace => Int[1],
        :move_counts => Dict{Symbol, Int}(:M => 1, :I => 0, :D => 0),
    )
    path = Mycelia.Rhizomorph.ViterbiDecodingResult(
        nothing, 0.0, diagnostics)
    partial_label = first(MetaGraphsNext.labels(graph))
    return Mycelia.ViterbiCorrectionResult(
        Any[Any[partial_label]],
        Mycelia.Rhizomorph.ViterbiDecodingResult[path],
        Dict{Symbol, Any}(),
    )
end

Test.@testset "Indel-aware pair-HMM Viterbi correction" begin
    Test.@testset "Milestone A: config indel fields default to substitution-collapse" begin
        cfg = Mycelia.ViterbiCorrectionConfig()
        # The gap machinery is OFF by default, so the decoder is byte-identical.
        Test.@test cfg.indel_moves == false
        Test.@test cfg.insertion_fraction == 0.0
        Test.@test cfg.deletion_fraction == 0.0
        Test.@test cfg.insertion_extend_probability == 0.0
        Test.@test cfg.deletion_extend_probability == 0.0
        Test.@test cfg.deletion_max_run == 0
        Test.@test cfg.max_insertion_run == 0
        Test.@test cfg.band_width === nothing

        # Fields are settable.
        cfg2 = Mycelia.ViterbiCorrectionConfig(
            indel_moves = true,
            insertion_fraction = 0.3,
            deletion_fraction = 0.3,
            insertion_extend_probability = 0.1,
            deletion_extend_probability = 0.1,
            deletion_max_run = 3,
            max_insertion_run = 3,
            band_width = 12
        )
        Test.@test cfg2.indel_moves == true
        Test.@test cfg2.insertion_fraction == 0.3
        Test.@test cfg2.deletion_fraction == 0.3
        Test.@test cfg2.deletion_max_run == 3
        Test.@test cfg2.max_insertion_run == 3
        Test.@test cfg2.band_width == 12
    end

    Test.@testset "length-changing acceptance compares mean k-mer support" begin
        comparable = Mycelia._comparable_likelihood_improvement
        # Ten original terms at -1 each. A deletion with nine terms and an
        # insertion with eleven terms have identical mean support, not a raw-sum
        # advantage/penalty caused solely by length.
        Test.@test comparable(-10.0, -9.0, 22, 21, 13;
            length_normalized = true) == 0.0
        Test.@test comparable(-10.0, -11.0, 22, 23, 13;
            length_normalized = true) == 0.0
        Test.@test comparable(-10.0, -8.1, 22, 21, 13;
            length_normalized = true) ≈ 0.1
        # The substitution oracle retains the historical raw-sum comparison.
        Test.@test comparable(-10.0, -9.0, 22, 21, 13;
            length_normalized = false) == 1.0
    end

    # Shared OLC fixture (mirrors viterbi_variable_length_graph_correction_test.jl):
    # a 5-vertex overlap-3 chain ATGCG -> GCGTA -> GTACC -> ACCGT -> CGTAA. The true
    # read under test covers only the FIRST THREE vertices; the graph deliberately
    # CONTINUES past that end (ACCGT, CGTAA) so a frame-shifted substitution-only
    # decode over-walks onto a wrong downstream vertex instead of coincidentally
    # stopping at the true terminal — making the insertion case discriminating.
    # The indel-enabled profile: err 0.10, indels split 30/30 of the error mass,
    # short affine gaps, exact (unbounded) beam so the assertion tests the true
    # argmax rather than the beam approximation.
    function _indel_olc_graph()
        records = [
            FASTX.FASTA.Record("read_1", BioSequences.dna"ATGCG"),
            FASTX.FASTA.Record("read_2", BioSequences.dna"GCGTA"),
            FASTX.FASTA.Record("read_3", BioSequences.dna"GTACC"),
            FASTX.FASTA.Record("read_4", BioSequences.dna"ACCGT"),
            FASTX.FASTA.Record("read_5", BioSequences.dna"CGTAA")
        ]
        return Mycelia.Rhizomorph.build_fasta_graph_olc(records; min_overlap = 3)
    end

    function _indel_config(;
            indel_moves = true,
            deletion_max_run = 3,
            max_insertion_run = 3
    )
        return Mycelia.ViterbiCorrectionConfig(
            error_rate = 0.10,
            indel_moves = indel_moves,
            insertion_fraction = 0.30,
            deletion_fraction = 0.30,
            insertion_extend_probability = 0.10,
            deletion_extend_probability = 0.10,
            deletion_max_run = deletion_max_run,
            max_insertion_run = max_insertion_run,
            beam_width = typemax(Int)
        )
    end

    Test.@testset "Milestone B: deletion-shifted read corrects back to reference" begin
        graph = _indel_olc_graph()
        # The middle read (GCGTA) is missing entirely — a unit-level deletion. The
        # terminal (GTACC) matches a graph vertex exactly, anchoring the endpoint
        # two vertices past the start. The substitution-only lock-step decoder
        # cannot reach that target in a single read step (it under-walks and never
        # reaches GTACC); the pair-HMM must take a D-move (advance the graph edge
        # ATGCG->GCGTA without consuming a read unit) and then match GTACC, landing
        # the 3-vertex reference walk.
        observed = [BioSequences.dna"ATGCG", BioSequences.dna"GTACC"]
        result = Mycelia.correct_observations(graph, [observed]; config = _indel_config())
        Test.@test indel_decoded_label_strings(result) == ["ATGCG", "GCGTA", "GTACC"]
    end

    Test.@testset "Milestone B: insertion-shifted read corrects back to reference" begin
        graph = _indel_olc_graph()
        # An extra spurious read unit (AAAAA) is inserted after the first — a
        # unit-level insertion. The last unit (GTACA) is a corrupted terminal (not a
        # graph vertex) so the decode is free-endpoint. The lock-step decoder treats
        # the junk as a real step and frame-shifts one vertex too far, over-walking
        # onto the spurious downstream vertex ACCGT ([ATGCG, GCGTA, GTACC, ACCGT]);
        # the pair-HMM must take an I-move (consume AAAAA while staying on ATGCG) and
        # recover the true 3-vertex walk, correcting GTACA back to GTACC.
        observed = [
            BioSequences.dna"ATGCG",
            BioSequences.dna"AAAAA",
            BioSequences.dna"GCGTA",
            BioSequences.dna"GTACA"
        ]
        result = Mycelia.correct_observations(graph, [observed]; config = _indel_config())
        Test.@test indel_decoded_label_strings(result) == ["ATGCG", "GCGTA", "GTACC"]
    end

    Test.@testset "Milestone B: substitution fixture stays byte-identical (indel_moves=false)" begin
        graph = _indel_olc_graph()
        # read_2 carries a substitution (GCGTT for GCGTA). Under indel_moves=false
        # the decoder is the untouched substitution path and must reproduce the
        # reference walk exactly, as the existing variable-length suite asserts.
        observed = [
            BioSequences.dna"ATGCG",
            BioSequences.dna"GCGTT",
            BioSequences.dna"GTACC"
        ]
        result = Mycelia.correct_observations(
            graph, [observed]; config = _indel_config(indel_moves = false))
        Test.@test indel_decoded_label_strings(result) == ["ATGCG", "GCGTA", "GTACC"]

        # And with indel_moves=true the same equal-length substitution fixture still
        # corrects to the reference (gap moves never beat the pure-M substitution
        # path when the read length matches the walk length).
        result_indel = Mycelia.correct_observations(
            graph, [observed]; config = _indel_config(indel_moves = true))
        Test.@test indel_decoded_label_strings(result_indel) ==
                   ["ATGCG", "GCGTA", "GTACC"]
    end

    Test.@testset "Milestone C: the decode is routed through the indel kernel" begin
        graph = _indel_olc_graph()
        del_result = Mycelia.correct_observations(
            graph, [[BioSequences.dna"ATGCG", BioSequences.dna"GTACC"]];
            config = _indel_config())
        diag = only(del_result.paths).diagnostics
        # A mis-wire to the substitution-only path would NOT stamp this algorithm.
        Test.@test diag[:algorithm] == :viterbi_indel_pair_hmm
        Test.@test diag[:indel_moves] == true
        # The deletion was bridged by >=1 D-move.
        Test.@test diag[:move_counts][:D] >= 1
        Test.@test count(==(:D), diag[:move_trace]) >= 1
        Test.@test length(diag[:move_trace]) == length(diag[:read_index_trace])

        ins_result = Mycelia.correct_observations(
            graph,
            [[BioSequences.dna"ATGCG", BioSequences.dna"AAAAA",
                BioSequences.dna"GCGTA", BioSequences.dna"GTACA"]];
            config = _indel_config())
        ins_diag = only(ins_result.paths).diagnostics
        Test.@test ins_diag[:move_counts][:I] >= 1
        Test.@test count(==(:I), ins_diag[:move_trace]) >= 1
        Test.@test length(ins_diag[:move_trace]) == length(ins_diag[:read_index_trace])
    end

    Test.@testset "Milestone C: bounding knobs actually bite" begin
        graph = _indel_olc_graph()
        deletion_observed = [BioSequences.dna"ATGCG", BioSequences.dna"GTACC"]
        insertion_observed = [
            BioSequences.dna"ATGCG", BioSequences.dna"AAAAA",
            BioSequences.dna"GCGTA", BioSequences.dna"GTACA"
        ]

        function _cfg(; kwargs...)
            base = (
                error_rate = 0.10,
                indel_moves = true,
                insertion_fraction = 0.30,
                deletion_fraction = 0.30,
                insertion_extend_probability = 0.10,
                deletion_extend_probability = 0.10,
                deletion_max_run = 3,
                max_insertion_run = 3,
                beam_width = typemax(Int)
            )
            return Mycelia.ViterbiCorrectionConfig(; base..., kwargs...)
        end

        # deletion_max_run = 0 forbids the D-move, so the anchored target GTACC
        # becomes unreachable and the decode returns no corrected path.
        no_del = Mycelia.correct_observations(
            graph, [deletion_observed]; config = _cfg(deletion_max_run = 0))
        Test.@test only(no_del.corrected_observations) === nothing
        # deletion_max_run >= 1 recovers it.
        yes_del = Mycelia.correct_observations(
            graph, [deletion_observed]; config = _cfg(deletion_max_run = 2))
        Test.@test indel_decoded_label_strings(yes_del) == ["ATGCG", "GCGTA", "GTACC"]

        # band_width = 0 forbids any net gap (|#deletions - #insertions| > 0), so the
        # deletion (net gap +1) is pruned and the target is unreachable.
        banded_out = Mycelia.correct_observations(
            graph, [deletion_observed]; config = _cfg(band_width = 0))
        Test.@test only(banded_out.corrected_observations) === nothing
        # A band wide enough to admit the single deletion recovers the reference.
        banded_in = Mycelia.correct_observations(
            graph, [deletion_observed]; config = _cfg(band_width = 4))
        Test.@test indel_decoded_label_strings(banded_in) == ["ATGCG", "GCGTA", "GTACC"]

        # max_insertion_run = 0 forbids the I-move, so the insertion cannot be
        # absorbed and the free-endpoint decode frame-shifts one vertex too far onto
        # the spurious downstream vertex ACCGT. Assert the EXACT over-walk (not just
        # "not the reference") so a compensating error that happened to differ from
        # the reference could not pass — mirroring the two-directional deletion knobs.
        no_ins = Mycelia.correct_observations(
            graph, [insertion_observed]; config = _cfg(max_insertion_run = 0))
        Test.@test indel_decoded_label_strings(no_ins) ==
                   ["ATGCG", "GCGTA", "GTACC", "ACCGT"]
        # max_insertion_run >= 1 admits the single I-move and recovers the reference.
        yes_ins = Mycelia.correct_observations(
            graph, [insertion_observed]; config = _cfg(max_insertion_run = 1))
        Test.@test indel_decoded_label_strings(yes_ins) == ["ATGCG", "GCGTA", "GTACC"]
    end

    Test.@testset "pair-HMM preserves distinct insertion-run states" begin
        graph = Mycelia.Rhizomorph.build_fasta_graph_olc(
            [
                FASTX.FASTA.Record("state_a", BioSequences.dna"ATGCG"),
                FASTX.FASTA.Record("state_b", BioSequences.dna"GCGTA"),
            ];
            min_overlap = 3,
        )
        observations = [
            BioSequences.dna"ATGCG",
            BioSequences.dna"GCGTA",
            BioSequences.dna"GCGTA",
            BioSequences.dna"AAAAA",
            BioSequences.dna"GCGTA",
        ]

        function state_alias_emission(
                observed_unit::Any,
                _node::Any,
                _alphabet::Symbol;
                quality::Any = nothing,
                error_rate::Float64 = 0.10,
        )::Float64
            return string(observed_unit) == "AAAAA" ? -Inf : 0.0
        end

        function state_alias_insertion_emission(
                _observed_unit::Any,
                _alphabet::Symbol,
        )::Float64
            return 0.0
        end

        config = Mycelia.ViterbiCorrectionConfig(
            error_rate = 0.10,
            emission_logp = state_alias_emission,
            strand_mode = :singlestrand,
            indel_moves = true,
            insertion_fraction = 0.30,
            deletion_fraction = 0.0,
            insertion_extend_probability = 0.10,
            deletion_extend_probability = 0.0,
            deletion_max_run = 0,
            max_insertion_run = 2,
            band_width = nothing,
            beam_width = typemax(Int),
            insertion_emission_logp = state_alias_insertion_emission,
        )

        result = Mycelia.correct_observations(
            graph,
            [observations];
            config = config,
        )
        path_result = only(result.paths)
        diagnostics = path_result.diagnostics

        # The only complete trace is M,I,M,I,I. At column 4, the I state reached
        # by opening a fresh run has run_length=1, while a higher-scoring competing
        # I state at the same vertex has run_length=2. They must remain distinct:
        # the run-2 state cannot consume column 5, but the run-1 state can.
        Test.@test !diagnostics[:truncated]
        Test.@test diagnostics[:decoded_read_index] == length(observations)
        Test.@test diagnostics[:completed_columns] == length(observations)
        Test.@test diagnostics[:move_trace] == [:M, :I, :M, :I, :I]
        Test.@test diagnostics[:read_index_trace] == collect(eachindex(observations))
        Test.@test diagnostics[:move_counts] ==
                   Dict{Symbol, Int}(:M => 2, :I => 3, :D => 0)
        Test.@test [
            string(step.vertex_label) for step in path_result.path.steps
        ] == ["ATGCG", "GCGTA"]
    end

    Test.@testset "Milestone C: nanopore-style read decodes through the kernel" begin
        # A real nanopore error process (indels + homopolymer boost) on a small
        # reference, threaded as k-mer observations through the clean-reference
        # k-mer graph. The point is that the indel kernel RUNS (decode-diagnostic
        # proves it, catching a silent mis-wire to the substitution-only path) and
        # completes quickly; correction quality at scale is a separate HPC concern.
        reference = Mycelia.random_fasta_record(moltype = :DNA, seed = 7, L = 60)
        k = 7
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([reference], k)

        refseq = FASTX.sequence(BioSequences.LongDNA{4}, reference)
        # observe returns a (sequence, quality) TUPLE — take [1] for the sequence.
        obs_seq, obs_quals = Mycelia.observe(refseq; error_rate = 0.05, tech = :nanopore)
        Test.@test obs_seq isa BioSequences.LongSequence

        qstr = String([Char(q + 33) for q in obs_quals])
        read = FASTX.FASTQ.Record("obs", string(obs_seq), qstr)
        sequence_string = FASTX.sequence(String, read)
        alphabet = Mycelia.detect_alphabet(sequence_string)
        sequence_type = Mycelia.alphabet_to_biosequence_type(alphabet)
        sequence = Mycelia.extract_typed_sequence(read, sequence_type)
        kmers = collect(Mycelia._record_kmer_iterator(sequence_type, k, sequence))
        Test.@test !isempty(kmers)
        quality_scores = collect(FASTX.quality_scores(read))
        n_qual = length(quality_scores)
        observations = [Mycelia.QualityObservation(
                            kmer,
                            UInt8.(quality_scores[clamp(i, 1, n_qual):clamp(i + k - 1, 1, n_qual)]))
                        for (i, kmer) in enumerate(kmers)]

        config = Mycelia.ViterbiCorrectionConfig(
            error_rate = 0.05,
            strand_mode = :singlestrand,
            indel_moves = true,
            insertion_fraction = 0.30,
            deletion_fraction = 0.30,
            insertion_extend_probability = 0.10,
            deletion_extend_probability = 0.10,
            deletion_max_run = 3,
            max_insertion_run = 3,
            band_width = 12,
            beam_width = 256
        )

        elapsed = @elapsed result = Mycelia.correct_observations(
            graph, [observations]; config = config)
        Test.@test elapsed < 60
        diag = only(result.paths).diagnostics
        Test.@test diag[:algorithm] == :viterbi_indel_pair_hmm
        Test.@test diag[:indel_moves] == true
        # A path was decoded and it consumed match moves (i.e. the kernel actually
        # walked the graph, not a degenerate empty result).
        Test.@test only(result.paths).path !== nothing
        Test.@test diag[:move_counts][:M] >= 1
    end

    Test.@testset "C1: real single-BASE indel corrects back to the reference path" begin
        # The PIVOTAL correctness test. The Milestone-B fixtures are UNIT-level
        # (whole k-mer inserted/deleted) and the nanopore smoke asserts only routing.
        # This asserts the ACTUAL nanopore failure mode: a single BASE indel in the
        # read corrupts the k k-mer windows spanning it (an indel + a local
        # substitution burst), and the decode must recover the clean reference's
        # k-mer walk EXACTLY. Start and end sit on clean k-mers (indel placed mid-
        # read) so both endpoints anchor; the indel kernel bridges the frame shift
        # with a single gap move and tolerates the corrupted-window mismatches.
        k = 7
        reference = Mycelia.random_fasta_record(moltype = :DNA, seed = 1, L = 40)
        reference_string = FASTX.sequence(String, reference)
        reference_units = indel_string_kmer_units(reference_string, k)
        # Non-repetitive reference => the graph is a simple path, so the true k-mer
        # walk is unambiguous and equality against it is a well-posed assertion.
        Test.@test length(unique(string.(reference_units))) == length(reference_units)

        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([reference], k)
        expected = [string(unit) for unit in reference_units]

        indel_position = 20
        deleted_string = reference_string[1:(indel_position - 1)] *
                         reference_string[(indel_position + 1):end]
        inserted_string = reference_string[1:indel_position] * "A" *
                          reference_string[(indel_position + 1):end]
        Test.@test length(deleted_string) == length(reference_string) - 1
        Test.@test length(inserted_string) == length(reference_string) + 1

        deleted_units = indel_string_kmer_units(deleted_string, k)
        inserted_units = indel_string_kmer_units(inserted_string, k)

        config = Mycelia.ViterbiCorrectionConfig(
            error_rate = 0.10,
            strand_mode = :singlestrand,
            indel_moves = true,
            insertion_fraction = 0.30,
            deletion_fraction = 0.30,
            insertion_extend_probability = 0.10,
            deletion_extend_probability = 0.10,
            deletion_max_run = 3,
            max_insertion_run = 3,
            beam_width = typemax(Int),
            band_width = nothing
        )

        deleted_result = Mycelia.correct_observations(
            graph, [deleted_units]; config = config)
        Test.@test indel_decoded_label_strings(deleted_result) == expected

        inserted_result = Mycelia.correct_observations(
            graph, [inserted_units]; config = config)
        Test.@test indel_decoded_label_strings(inserted_result) == expected
    end

    Test.@testset "over-correction guard: clean read takes ZERO gap moves" begin
        # A clean, error-free, equal-length read through an indel-ENABLED config
        # (with full gap capacity) must decode with NO gap moves. The Milestone-B
        # clean-read assertion checks the labels but not the gap counts, so a
        # compensating insertion+deletion pair could sneak past it; asserting zero I
        # AND zero D closes that hole.
        graph = _indel_olc_graph()
        observed = [
            BioSequences.dna"ATGCG",
            BioSequences.dna"GCGTA",
            BioSequences.dna"GTACC"
        ]
        result = Mycelia.correct_observations(graph, [observed]; config = _indel_config())
        move_counts = only(result.paths).diagnostics[:move_counts]
        Test.@test move_counts[:I] == 0
        Test.@test move_counts[:D] == 0
        Test.@test indel_decoded_label_strings(result) == ["ATGCG", "GCGTA", "GTACC"]
    end

    Test.@testset "oracle battery: substitution inputs preserved byte-for-byte" begin
        # Strengthen oracle preservation from a single fixture to a battery of
        # DIFFERENT substitution inputs, asserting each decodes IDENTICALLY under
        # (a) the DEFAULT substitution decoder, (b) the indel config with moves OFF,
        # and (c) the indel KERNEL run with all gap mass ZERO (the mathematical
        # collapse). Equality is asserted against the substitution decoder's ACTUAL
        # output, not a hand-written answer, so the collapse cannot silently diverge.
        function _assert_oracle_equivalent(graph, observation)
            substitution = Mycelia.correct_observations(
                graph, [observation];
                config = Mycelia.ViterbiCorrectionConfig(error_rate = 0.10))
            substitution_labels = indel_decoded_label_strings(substitution)

            moves_off = Mycelia.correct_observations(
                graph, [observation]; config = _indel_config(indel_moves = false))
            Test.@test indel_decoded_label_strings(moves_off) == substitution_labels

            # indel_moves=true with every gap mass 0 routes through the indel kernel
            # but collapses to the substitution path. The zero-capacity config trips
            # the new silent-no-op @warn (FIX code-1), so absorb it — the collapse is
            # a deliberate, byte-identity-proving configuration, not a misconfig.
            collapsed = Logging.with_logger(Logging.NullLogger()) do
                Mycelia.correct_observations(
                    graph, [observation];
                    config = Mycelia.ViterbiCorrectionConfig(
                        error_rate = 0.10,
                        indel_moves = true,
                        insertion_fraction = 0.0,
                        deletion_fraction = 0.0,
                        deletion_max_run = 0,
                        max_insertion_run = 0,
                        beam_width = typemax(Int)))
            end
            Test.@test indel_decoded_label_strings(collapsed) == substitution_labels
            return substitution_labels
        end

        graph = _indel_olc_graph()
        # Multi-vertex substitution (middle unit corrupted; endpoint anchored).
        _assert_oracle_equivalent(graph,
            [BioSequences.dna"ATGCG", BioSequences.dna"GCGTT", BioSequences.dna"GTACC"])
        # Free-endpoint case (last unit GTACG is not a graph vertex => endpoint free).
        _assert_oracle_equivalent(graph,
            [BioSequences.dna"ATGCG", BioSequences.dna"GCGTA", BioSequences.dna"GTACG"])

        # Quality-aware case: a FASTQ OLC graph with a read-quality observation, so
        # the emission model is :quality_aware rather than :alphabet_parameterized.
        quality_records = [
            indel_quality_record("read_1", "ATGCG", [40, 40, 40, 40, 40]),
            indel_quality_record("read_2", "GCGTA", [40, 40, 40, 40, 40]),
            indel_quality_record("read_3", "GTACC", [40, 40, 40, 40, 40])
        ]
        quality_graph = Mycelia.Rhizomorph.build_fastq_graph_olc(
            quality_records; min_overlap = 3)
        _assert_oracle_equivalent(quality_graph,
            [
                BioSequences.dna"ATGCG",
                IndelQualityObservation(BioSequences.dna"GCGTT", UInt8[2, 2, 2, 2, 2]),
                BioSequences.dna"GTACC"
            ])
    end

    Test.@testset "D_max: a deletion run exceeding the cap fails gracefully" begin
        # A read that skips THREE vertices ([ATGCG, CGTAA] on ATGCG->GCGTA->GTACC->
        # ACCGT->CGTAA) needs a deletion RUN of 2 within a read column to bridge to
        # the anchored terminal CGTAA. deletion_max_run = 1 caps the run below that,
        # so the target is unreachable and the decode returns `nothing` GRACEFULLY
        # (no crash); deletion_max_run >= 2 recovers the full 5-vertex walk.
        graph = _indel_olc_graph()
        observed = [BioSequences.dna"ATGCG", BioSequences.dna"CGTAA"]
        capped = Mycelia.correct_observations(
            graph, [observed];
            config = _indel_config(deletion_max_run = 1))
        Test.@test only(capped.corrected_observations) === nothing
        recovered = Mycelia.correct_observations(
            graph, [observed];
            config = _indel_config(deletion_max_run = 2))
        Test.@test indel_decoded_label_strings(recovered) ==
                   ["ATGCG", "GCGTA", "GTACC", "ACCGT", "CGTAA"]
    end

    Test.@testset "truncation is a first-class diagnostic" begin
        # A free-endpoint read whose frontier DIES mid-decode returns the best-so-far
        # prefix; `diagnostics[:truncated]` must make that explicit (FIX code-2). The
        # read starts on the terminal vertex CGTAA (no outgoing edges) then presents a
        # junk unit TTTTT; with insertion disabled there is no M, I, or D move, so the
        # frontier dies at read unit 1 and the decode is truncated but still returns a
        # path (the CGTAA prefix).
        graph = _indel_olc_graph()
        truncated_observed = [BioSequences.dna"CGTAA", BioSequences.dna"TTTTT"]
        truncated = Mycelia.correct_observations(
            graph, [truncated_observed];
            config = _indel_config(max_insertion_run = 0))
        truncated_diag = only(truncated.paths).diagnostics
        Test.@test only(truncated.corrected_observations) !== nothing
        Test.@test truncated_diag[:truncated] == true
        Test.@test truncated_diag[:decoded_read_index] < length(truncated_observed)

        # A normal full decode consumes every read unit => not truncated.
        full_observed = [
            BioSequences.dna"ATGCG",
            BioSequences.dna"AAAAA",
            BioSequences.dna"GCGTA",
            BioSequences.dna"GTACA"
        ]
        full = Mycelia.correct_observations(
            graph, [full_observed]; config = _indel_config())
        full_diag = only(full.paths).diagnostics
        Test.@test full_diag[:truncated] == false
        Test.@test full_diag[:decoded_read_index] == length(full_observed)
    end

    Test.@testset "truncated read corrections fail closed" begin
        reference = indel_quality_record("ref", "CGTAA", fill(40, 5))
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            [reference], 5; mode = :doublestrand)
        observed = indel_quality_record("obs", "CGTAAT", fill(40, 6))
        params = Mycelia.IndelDecodeParams(
            0.10, 0.30, 0.30, 0.10, 0.10, 1, 0, 16)
        diagnostics = Mycelia.CorrectorDiagnostics()
        soft_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        result = Mycelia.try_viterbi_path_improvement(
            observed, graph, 5; graph_mode = :doublestrand,
            beam_width = typemax(Int), soft_weights = soft_weights,
            diagnostics = diagnostics,
            indel_params = params)
        Test.@test result === nothing
        Test.@test isempty(soft_weights.weights)
        Test.@test diagnostics.structural_errors[] == 0
        Test.@test diagnostics.indel_attempts[] == 1
        Test.@test diagnostics.truncated_decodes[] == 1
        Test.@test diagnostics.trace_contract_errors[] == 0
        Test.@test diagnostics.indel_decodes[] == 0
        Test.@test diagnostics.truncated_decodes[] ==
                   diagnostics.indel_attempts[] - diagnostics.indel_decodes[]
    end

    Test.@testset "no-path and malformed results have terminal outcomes" begin
        reference = indel_quality_record("ref", "CGTAA", fill(40, 5))
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            [reference], 5; mode = :doublestrand)
        observed = indel_quality_record("obs", "CGTAA", fill(40, 5))
        params = Mycelia.IndelDecodeParams(
            0.10, 0.30, 0.30, 0.10, 0.10, 1, 1, 16)

        no_path_diagnostics = Mycelia.CorrectorDiagnostics()
        no_path = Mycelia.try_viterbi_path_improvement(
            observed,
            graph,
            5;
            graph_mode = :doublestrand,
            diagnostics = no_path_diagnostics,
            indel_params = params,
            correction_kernel = indel_no_path_correction_kernel,
        )
        Test.@test no_path === nothing
        Test.@test no_path_diagnostics.indel_attempts[] == 1
        Test.@test no_path_diagnostics.indel_decodes[] == 0
        Test.@test no_path_diagnostics.truncated_decodes[] == 1
        Test.@test no_path_diagnostics.trace_contract_errors[] == 0

        empty_path_diagnostics = Mycelia.CorrectorDiagnostics()
        empty_path = Mycelia.try_viterbi_path_improvement(
            observed,
            graph,
            5;
            graph_mode = :doublestrand,
            diagnostics = empty_path_diagnostics,
            indel_params = params,
            correction_kernel = indel_empty_path_correction_kernel,
        )
        Test.@test empty_path === nothing
        Test.@test empty_path_diagnostics.indel_attempts[] == 1
        Test.@test empty_path_diagnostics.indel_decodes[] == 0
        Test.@test empty_path_diagnostics.truncated_decodes[] == 1
        Test.@test empty_path_diagnostics.trace_contract_errors[] == 1

        # Reviewer regression: a kernel can return a nonempty one-column path
        # while falsely declaring `truncated=false` for a five-column request.
        # Reject it before path reconstruction, quality mapping, or soft-EM.
        partial_reference = indel_quality_record(
            "partial_ref", "ACGTACGTA", fill(40, 9))
        partial_graph = Mycelia.Rhizomorph.build_qualmer_graph(
            [partial_reference], 5; mode = :doublestrand)
        partial_diagnostics = Mycelia.CorrectorDiagnostics()
        partial_soft_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        partial_result = Mycelia.try_viterbi_path_improvement(
            partial_reference,
            partial_graph,
            5;
            graph_mode = :doublestrand,
            diagnostics = partial_diagnostics,
            soft_weights = partial_soft_weights,
            indel_params = params,
            correction_kernel = indel_partial_complete_correction_kernel,
        )
        Test.@test partial_result === nothing
        Test.@test isempty(partial_soft_weights.weights)
        Test.@test partial_diagnostics.indel_attempts[] == 1
        Test.@test partial_diagnostics.indel_decodes[] == 0
        Test.@test partial_diagnostics.truncated_decodes[] == 1
        Test.@test partial_diagnostics.trace_contract_errors[] == 1
        Test.@test partial_diagnostics.indel_engaged[] == 0
        Test.@test Mycelia._quality_from_indel_trace(
            repeat("I", 9), Symbol[:M], Int[1], 5, 5) === nothing

        empty_graph = deepcopy(graph)
        while !isempty(collect(MetaGraphsNext.labels(empty_graph)))
            MetaGraphsNext.rem_vertex!(empty_graph, 1)
        end
        exception_diagnostics = Mycelia.CorrectorDiagnostics()
        exception_result = Mycelia.try_viterbi_path_improvement(
            observed,
            empty_graph,
            5;
            graph_mode = :doublestrand,
            diagnostics = exception_diagnostics,
            indel_params = params,
        )
        Test.@test exception_result === nothing
        Test.@test exception_diagnostics.indel_attempts[] == 1
        Test.@test exception_diagnostics.indel_decodes[] == 0
        Test.@test exception_diagnostics.truncated_decodes[] == 1
        Test.@test exception_diagnostics.trace_contract_errors[] == 1
        Test.@test exception_diagnostics.structural_errors[] == 1
    end
end
