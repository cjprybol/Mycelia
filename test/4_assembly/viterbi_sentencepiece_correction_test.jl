# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_sentencepiece_correction_test.jl")'
# ```
#
# Optional real-SentencePiece leg (requires the conda `sentencepiece_env`):
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/4_assembly/viterbi_sentencepiece_correction_test.jl")'
# ```
#
# Wires the SentencePiece -> token-graph -> Viterbi-corrector pipeline that the
# cross-modality audit flagged as the one genuinely-untested universal claim.
#
# The always-on leg uses hand-built ASCII token sequences that mimic
# SentencePiece subword pieces so the token-graph -> corrector seam is exercised
# in default CI without the external `sentencepiece_env`. The gated leg trains a
# real tiny model, encodes sentences, and runs the same pipeline end-to-end.
#
# Finding 1 (recorded during authoring): `build_token_graph` produces
# `StringVertexData` vertices (in the corrector's accepted vertex-data union) and
# `StringEdgeData` edges (in the variable-length-edge union). The token graph is
# therefore classified as TEXT / variable-length and threads through
# `correct_observations` unchanged — NO core (viterbi-next.jl) or
# string-graphs.jl edge-data fix was needed.
#
# Finding 2 (RESOLVED, td-oj72): the TEXT emission model in viterbi-next.jl used
# to index token strings by BYTE, so a token containing a multibyte UTF-8
# character threw `StringIndexError`. Real SentencePiece emits the U+2581 ("▁",
# 3 bytes) word-boundary marker on word-initial pieces, tripping this on raw
# token graphs. The emission model (and `_hamming_like_distance`) now index by
# CHARACTER, so unstripped multibyte pieces correct end-to-end — see the
# "multibyte UTF-8 tokens correct" always-on testset below. The gated leg still
# strips U+2581 for a clean word-model vocabulary, not because it must.

import Mycelia
import Test

# SentencePiece prints this 3-byte marker before word-initial pieces. The
# correction path now handles it natively (see Finding 2 in the header); the
# gated leg still strips it only for a cleaner word-model vocabulary.
const _SP_WORD_BOUNDARY = "▁"

function _sp_decoded_labels(result::Mycelia.ViterbiCorrectionResult)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

Test.@testset "SentencePiece token-graph Viterbi correction" begin
    Test.@testset "token-graph -> corrector wiring (simulated pieces)" begin
        # ASCII subword pieces (SentencePiece-style, boundary marker stripped).
        # Two sentences share the interior token so the graph is a genuine
        # (non-trivial) token adjacency graph, not a single linear chain.
        token_sequences = [
            ["the", "quick", "brown", "fox", "jumps"],
            ["the", "lazy", "brown", "dog", "sleeps"],
        ]
        graph = Mycelia.Rhizomorph.build_token_graph(
            token_sequences; dataset_id = "sp_tokens"
        )

        # The wiring facts the corrector relies on to classify the token graph.
        Test.@test Mycelia._correction_vertex_data_type(graph) ==
                   Mycelia.Rhizomorph.StringVertexData
        Test.@test Mycelia._correction_edge_data_type(graph) ==
                   Mycelia.Rhizomorph.StringEdgeData

        # Corrupt exactly ONE interior token with a same-length substitution
        # ("brown" -> "brawn"). Start and end tokens stay valid graph vertices.
        observed = ["the", "quick", "brawn", "fox", "jumps"]
        result = Mycelia.correct_observations(graph, [observed])

        Test.@test result.diagnostics[:alphabet] == :TEXT
        Test.@test result.diagnostics[:strand_mode] == :singlestrand
        Test.@test result.diagnostics[:emission_model] == :alphabet_parameterized
        Test.@test result.diagnostics[:transition_model] == :normalized_overlap_length
        # The corrupted token is restored to "brown".
        Test.@test _sp_decoded_labels(result) ==
                   ["the", "quick", "brown", "fox", "jumps"]
    end

    Test.@testset "multibyte UTF-8 tokens correct (U+2581 boundary marker)" begin
        # Regression for the byte-indexed TEXT emission model (former Finding 2):
        # tokens carrying the raw SentencePiece U+2581 ("▁", 3 bytes) word-
        # boundary marker previously threw `StringIndexError` because the emission
        # model indexed token strings by byte. The corrector now indexes by
        # character, so unstripped multibyte pieces correct end-to-end.
        token_sequences = [
            ["▁the", "▁quick", "▁brown", "▁fox", "▁jumps"],
            ["▁the", "▁lazy", "▁brown", "▁dog", "▁sleeps"],
        ]
        graph = Mycelia.Rhizomorph.build_token_graph(
            token_sequences; dataset_id = "sp_unicode_tokens"
        )

        # Corrupt exactly ONE interior token with a same-length (same character
        # count) substitution that is off-graph: "▁brown" -> "▁brawn". The marker
        # itself stays intact, exercising the multibyte code path in emission
        # scoring for every position.
        observed = ["▁the", "▁quick", "▁brawn", "▁fox", "▁jumps"]
        result = Mycelia.correct_observations(graph, [observed])

        Test.@test result.diagnostics[:alphabet] == :TEXT
        Test.@test _sp_decoded_labels(result) ==
                   ["▁the", "▁quick", "▁brown", "▁fox", "▁jumps"]

        # A second alphabet: accented Latin + emoji, to confirm the fix is not
        # specific to U+2581. Same-length off-graph corruption of an interior
        # token restores to the supported graph token.
        emoji_sequences = [
            ["café", "☕", "☀️", "π"],
            ["café", "☕", "🌙", "π"],
        ]
        emoji_graph = Mycelia.Rhizomorph.build_token_graph(
            emoji_sequences; dataset_id = "sp_emoji_tokens"
        )
        emoji_observed = ["café", "☕", "☀️", "ρ"]
        emoji_result = Mycelia.correct_observations(emoji_graph, [emoji_observed])
        Test.@test emoji_result.diagnostics[:alphabet] == :TEXT
        Test.@test _sp_decoded_labels(emoji_result) == ["café", "☕", "☀️", "π"]
    end

    Test.@testset "TEXT token graph rejects reverse-complement strand modes" begin
        # Tokens are reverse-complement naive: a doublestrand/canonical request
        # must throw (locks the RC-naive contract for the token-graph path too).
        toks = [["alpha", "beta", "gamma"]]
        graph = Mycelia.Rhizomorph.build_token_graph(toks; dataset_id = "sp_rc_guard")
        Test.@test_throws ArgumentError Mycelia.correct_observations(
            graph,
            [["alpha", "beta", "gamma"]];
            config = Mycelia.ViterbiCorrectionConfig(strand_mode = :doublestrand)
        )
    end

    # ------------------------------------------------------------------
    # Gated real-SentencePiece leg: trains a tiny model, encodes sentences,
    # builds the token graph from the encoded pieces, corrupts one token, and
    # asserts the pipeline runs and (for a linear token path) restores it.
    # ------------------------------------------------------------------
    run_all = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true"
    run_external = run_all ||
                   lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"

    if run_external && Mycelia._check_sentencepiece_installed()
        Test.@testset "real SentencePiece encode -> token-graph -> correction" begin
            # A small corpus whose word-level tokens repeat across sentences so
            # a word-model vocabulary is learnable at tiny scale.
            corpus_sentences = [
                "the quick brown fox jumps over the lazy dog",
                "the lazy brown dog sleeps under the quick fox",
                "a quick brown fox and a lazy brown dog",
            ]
            corpus_file = tempname() * ".txt"
            model_prefix = tempname()
            try
                write(corpus_file, join(corpus_sentences, "\n") * "\n")
                trained = Mycelia.train_sentencepiece_model(
                    input_file = corpus_file,
                    model_prefix = model_prefix,
                    vocab_size = 32,
                    model_type = :word
                )
                Test.@test isfile(trained.model_file)

                encoded = Mycelia.encode_sentencepiece(
                    model_file = trained.model_file,
                    input = corpus_sentences,
                    output_format = :pieces
                )
                Test.@test encoded isa Vector{Vector{String}}
                # Strip the U+2581 word-boundary marker for a clean word-model
                # vocabulary. The corrector handles multibyte pieces natively
                # now (Finding 2 resolved); stripping is a vocabulary choice.
                token_sequences = [
                    [replace(piece, _SP_WORD_BOUNDARY => "") for piece in pieces]
                    for pieces in encoded
                ]
                token_sequences = [
                    filter(!isempty, seq) for seq in token_sequences
                ]

                graph = Mycelia.Rhizomorph.build_token_graph(
                    token_sequences; dataset_id = "sp_real_tokens"
                )
                Test.@test Mycelia._correction_vertex_data_type(graph) ==
                           Mycelia.Rhizomorph.StringVertexData

                # Pick a sentence whose token path has no repeated token (a
                # linear chain) so the corrupted-token restoration is
                # deterministic.
                linear_idx = findfirst(
                    seq -> length(seq) >= 3 && length(unique(seq)) == length(seq),
                    token_sequences
                )
                if linear_idx === nothing
                    # Still assert the pipeline runs end-to-end on the whole set.
                    result = Mycelia.correct_observations(graph, token_sequences)
                    Test.@test result.diagnostics[:alphabet] == :TEXT
                else
                    path_tokens = token_sequences[linear_idx]
                    mid = cld(length(path_tokens), 2)
                    original = path_tokens[mid]
                    # Same-length substitution: flip the final character to a
                    # different ASCII letter, guaranteeing an off-graph token.
                    last_char = original[end]
                    replacement = last_char == 'z' ? 'y' : 'z'
                    corrupted_token = original[1:(end - 1)] * string(replacement)
                    corrupted = copy(path_tokens)
                    corrupted[mid] = corrupted_token

                    result = Mycelia.correct_observations(graph, [corrupted])
                    Test.@test result.diagnostics[:alphabet] == :TEXT
                    Test.@test _sp_decoded_labels(result) == path_tokens
                end
            finally
                isfile(corpus_file) && rm(corpus_file)
                isfile(model_prefix * ".model") && rm(model_prefix * ".model")
                isfile(model_prefix * ".vocab") && rm(model_prefix * ".vocab")
            end
        end
    else
        Test.@testset "real SentencePiece leg (skipped)" begin
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true with sentencepiece_env installed"
        end
    end
end
