# # Tutorial 16: UN Corpus N-gram vs SentencePiece Token Graphs
#
# This tutorial compares two text-graph constructions over the UN Parallel Corpus:
# 1) Character n-gram graphs per language.
# 2) SentencePiece token graphs per language (token adjacency).
#
# The goals are to:
# - Build n-gram graphs across several `n` values for each language.
# - Summarize graph statistics/topology.
# - Train SentencePiece models per language and build token graphs with multiple vocab sizes.
# - Compare topology metrics between n-gram and token graphs.
#
# Generate a notebook from this Literate script:
# ```bash
# julia --project=docs -e 'using Literate; Literate.notebook("tutorials/16_un_corpus_ngram_vs_token_graphs.jl", "tutorials"; execute=false)'
# ```
#
# NOTE: UN Corpus downloads and SentencePiece training are resource-intensive and opt-in.
# Set `MYCELIA_RUN_EXTERNAL=true` to run those sections.

# ## Setup

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Dates
import Mycelia
import Test

# External tool gates (UN corpus download + SentencePiece).
const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

# ## Configuration

const DATA_ROOT = joinpath(homedir(), "workspace", "un_corpus")
const OUTDIR = joinpath(DATA_ROOT, "un_corpus_graphs")
const UN_SUBSETS = ["all"]
const NGRAM_SIZES = [3, 4, 5]
const TOKEN_VOCAB_SIZES = [2000, 8000]
const MAX_LINES_PER_LANGUAGE = nothing  # Set to an Int for quick runs

println("Workspace: ", OUTDIR)

println("""
Cleanup (optional):
  rm -rf $(OUTDIR)
""")

# ## Step 1: Download the UN Parallel Corpus (Opt-in)
#
# Uses the built-in downloader from `Mycelia.download_un_parallel_corpus`.

corpus_dir = joinpath(OUTDIR, "un_corpus")
if RUN_EXTERNAL
    mkpath(OUTDIR)
    corpus_dir = Mycelia.download_un_parallel_corpus(
        outdir=OUTDIR,
        subsets=UN_SUBSETS,
        formats=["txt"]
    )
    println("UN corpus downloaded to: ", corpus_dir)
else
    println("Skipping UN corpus download. Set MYCELIA_RUN_EXTERNAL=true to run.")
    if !isdir(corpus_dir)
        error("UN corpus not found at $(corpus_dir). Download it or set MYCELIA_RUN_EXTERNAL=true.")
    end
end

# ## Step 2: Identify per-language corpora
#
# Map each language to its raw text files in the extracted UN corpus.

function collect_un_language_files(corpus_root::String; subsets::Vector{String})
    language_files = Dict{String, Vector{String}}()
    use_all = "all" in subsets
    for (root, _dirs, files) in walkdir(corpus_root)
        for filename in files
            filepath = joinpath(root, filename)
            if endswith(filename, ".xml") || endswith(filename, ".gz") || endswith(filename, ".tar")
                continue
            end
            if !occursin("UNv1.0", filename)
                continue
            end
            if !occursin("-", filename) && !occursin("6way", filename)
                continue
            end
            parts = split(filename, '.')
            if length(parts) < 2
                continue
            end
            lang = parts[end]
            if length(lang) > 3
                continue
            end
            if !all(isletter, lang)
                continue
            end
            if !use_all
                matches_subset = any(subset -> occursin(subset, filename), subsets)
                if !matches_subset
                    continue
                end
            end
            push!(get!(Vector{String}, language_files, lang), filepath)
        end
    end
    return language_files
end

language_files = collect_un_language_files(corpus_dir; subsets=UN_SUBSETS)
if isempty(language_files)
    error("No language files found. Verify UN corpus extraction in $(corpus_dir).")
end

println("Languages discovered: ", join(sort(collect(keys(language_files))), ", "))

# ## Step 3: Build n-gram graphs per language
#
# Each language is processed independently. Each line is treated as one observation.

ngram_graphs = Dict{Tuple{String, Int}, Any}()
ngram_stats = Dict{Tuple{String, Int}, Dict{Symbol, Int}}()

for (lang, files) in language_files
    strings = String[]
    for filepath in files
        append!(strings, filter(!isempty, readlines(filepath)))
        if !isnothing(MAX_LINES_PER_LANGUAGE) && length(strings) >= MAX_LINES_PER_LANGUAGE
            strings = strings[1:MAX_LINES_PER_LANGUAGE]
            break
        end
    end
    Test.@test !isempty(strings)
    for n in NGRAM_SIZES
        dataset_id = "un_$(lang)_n$(n)"
        graph = Mycelia.Rhizomorph.build_ngram_graph(strings, n; dataset_id=dataset_id)
        ngram_graphs[(lang, n)] = graph
        ngram_stats[(lang, n)] = Mycelia.Rhizomorph.get_graph_statistics(graph)
        println("N-gram graph: $(lang) n=$(n) vertices=$(ngram_stats[(lang, n)][:vertex_count]) edges=$(ngram_stats[(lang, n)][:edge_count])")
    end
end

# ## Step 4: Summarize n-gram topology
#
# Print core metrics for quick inspection.

for ((lang, n), stats) in sort(collect(ngram_stats); by=x -> x[1])
    println("N-gram stats $(lang) n=$(n): ", stats)
end

# ## Step 5: Train SentencePiece models per language (Opt-in)
#
# SentencePiece is used to tokenize each language independently.
# This yields token sequences for graph construction.

sentencepiece_models = Dict{Tuple{String, Int}, String}()

if RUN_EXTERNAL
    for (lang, files) in language_files
        strings = String[]
        for filepath in files
            append!(strings, filter(!isempty, readlines(filepath)))
            if !isnothing(MAX_LINES_PER_LANGUAGE) && length(strings) >= MAX_LINES_PER_LANGUAGE
                strings = strings[1:MAX_LINES_PER_LANGUAGE]
                break
            end
        end
        input_file = joinpath(OUTDIR, "un_$(lang)_sentences.txt")
        open(input_file, "w") do io
            for line in strings
                println(io, line)
            end
        end
        for vocab_size in TOKEN_VOCAB_SIZES
            model_prefix = joinpath(OUTDIR, "spm_$(lang)_v$(vocab_size)")
            result = Mycelia.train_sentencepiece_model(
                input_file=input_file,
                model_prefix=model_prefix,
                vocab_size=vocab_size,
                model_type=:bpe,
                character_coverage=1.0
            )
            sentencepiece_models[(lang, vocab_size)] = result.model_file
        end
    end
else
    println("Skipping SentencePiece training. Set MYCELIA_RUN_EXTERNAL=true to run.")
end

# ## Step 6: Build token adjacency graphs
#
# Each token is a vertex; edges connect adjacent tokens in the sequence.

token_graphs = Dict{Tuple{String, Int}, Any}()
token_stats = Dict{Tuple{String, Int}, Dict{Symbol, Int}}()

if RUN_EXTERNAL
    for (lang, files) in language_files
        strings = String[]
        for filepath in files
            append!(strings, filter(!isempty, readlines(filepath)))
            if !isnothing(MAX_LINES_PER_LANGUAGE) && length(strings) >= MAX_LINES_PER_LANGUAGE
                strings = strings[1:MAX_LINES_PER_LANGUAGE]
                break
            end
        end
        Test.@test !isempty(strings)
        for vocab_size in TOKEN_VOCAB_SIZES
            model_file = sentencepiece_models[(lang, vocab_size)]
            token_sequences = Mycelia.encode_sentencepiece(
                model_file=model_file,
                input=strings,
                output_format=:pieces
            )
            graph = Mycelia.Rhizomorph.build_token_graph(
                token_sequences;
                dataset_id="un_$(lang)_tok_v$(vocab_size)"
            )
            token_graphs[(lang, vocab_size)] = graph
            token_stats[(lang, vocab_size)] = Mycelia.Rhizomorph.get_graph_statistics(graph)
            println("Token graph: $(lang) vocab=$(vocab_size) vertices=$(token_stats[(lang, vocab_size)][:vertex_count]) edges=$(token_stats[(lang, vocab_size)][:edge_count])")
        end
    end
else
    println("Skipping token graphs. SentencePiece models are required.")
end

# ## Step 7: Compare n-gram vs token graph topology
#
# Join n-gram and token graph statistics by language and parameter.

for (lang, _files) in language_files
    for n in NGRAM_SIZES
        if haskey(ngram_stats, (lang, n))
            stats = ngram_stats[(lang, n)]
            println("Compare $(lang) n=$(n): vertices=$(stats[:vertex_count]) edges=$(stats[:edge_count]) branches=$(stats[:branch_count]) joins=$(stats[:join_count])")
        end
    end
    for vocab_size in TOKEN_VOCAB_SIZES
        if haskey(token_stats, (lang, vocab_size))
            stats = token_stats[(lang, vocab_size)]
            println("Compare $(lang) vocab=$(vocab_size): vertices=$(stats[:vertex_count]) edges=$(stats[:edge_count]) branches=$(stats[:branch_count]) joins=$(stats[:join_count])")
        end
    end
end

# ## Appendix (Planned): Annotation Graphs
#
# Similar to token graphs, but vertices represent ordered protein annotations.
# Example path: UniRef_annotation_1 -> UniRef_annotation_2 -> UniRef_annotation_3.
# The annotation graph builder should support String/AA/DNA/RNA sequences and
# strand handling for DNA/RNA (singlestrand, doublestrand, canonical).
