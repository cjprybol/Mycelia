# SentencePiece integration for biological sequence tokenization
# Provides subword encoding for DNA, RNA, amino acid sequences, and text

# ============================================================================
# Environment Setup Functions
# ============================================================================

"""
    _setup_sentencepiece_environment(force_reinstall::Bool=false)

Set up the SentencePiece conda environment if it doesn't exist or if force_reinstall is true.

Creates a conda environment named 'sentencepiece_env' with Python and installs
the sentencepiece package via pip.

# Arguments
- `force_reinstall::Bool=false`: If true, removes and recreates the environment

# Details
This function is called automatically by wrapper functions when SentencePiece
is not installed. The environment uses Python 3.9 and pip to install sentencepiece.
"""
function _setup_sentencepiece_environment(force_reinstall::Bool=false)
    env_name = "sentencepiece_env"

    if force_reinstall && _check_conda_env_exists(env_name)
        println("Removing existing SentencePiece environment...")
        try
            Base.run(`$(Mycelia.CONDA_RUNNER) env remove -n $env_name -y`)
        catch e
            @warn "Failed to remove existing environment: $e"
        end
    end

    if !_check_conda_env_exists(env_name) || force_reinstall
        println("Creating SentencePiece conda environment...")
        try
            Base.run(`$(Mycelia.CONDA_RUNNER) create -y -n $env_name python=3.9`)
            println("Installing SentencePiece...")
            Base.run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $env_name pip install sentencepiece`)
            println("SentencePiece environment setup completed")
        catch e
            throw(ErrorException("Failed to setup SentencePiece environment: $e"))
        end
    else
        println("SentencePiece environment already exists")
    end
end

"""
    _check_sentencepiece_installed() -> Bool

Check if SentencePiece is installed and accessible.

Returns `true` if sentencepiece can be imported in the conda environment.
"""
function _check_sentencepiece_installed()
    env_name = "sentencepiece_env"
    if !_check_conda_env_exists(env_name)
        return false
    end
    try
        result = Base.read(`$(Mycelia.CONDA_RUNNER) run -n $env_name python -c "import sentencepiece; print('ok')"`, String)
        return occursin("ok", result)
    catch
        return false
    end
end

"""
    _ensure_sentencepiece_installed()

Ensure SentencePiece is installed, installing it if necessary.

This function is called automatically by all public SentencePiece wrapper functions.
"""
function _ensure_sentencepiece_installed()
    if !_check_sentencepiece_installed()
        _setup_sentencepiece_environment()
    end
end

# ============================================================================
# Sequence Conversion Helpers
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a BioSequence to a string for SentencePiece processing.

This is the interface point where BioSequences are converted to strings,
as required for the external SentencePiece tool.

# Arguments
- `sequence::BioSequences.LongDNA`: Input DNA sequence

# Returns
- `String`: Uppercase string representation of the sequence
"""
function sentencepiece_biosequence_to_string(sequence::BioSequences.LongDNA)::String
    return uppercase(string(sequence))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert an RNA BioSequence to a string for SentencePiece processing.
"""
function sentencepiece_biosequence_to_string(sequence::BioSequences.LongRNA)::String
    return uppercase(string(sequence))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert an amino acid BioSequence to a string for SentencePiece processing.
"""
function sentencepiece_biosequence_to_string(sequence::BioSequences.LongAA)::String
    return uppercase(string(sequence))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Detect the type of input for SentencePiece processing.

# Arguments
- `input`: A BioSequence or string

# Returns
- `:dna` for DNA sequences
- `:rna` for RNA sequences
- `:aa` for amino acid sequences
- `:ascii` for ASCII text strings
- `:unicode` for Unicode text strings
"""
function detect_sentencepiece_input_type(input::BioSequences.LongDNA)
    return :dna
end

function detect_sentencepiece_input_type(input::BioSequences.LongRNA)
    return :rna
end

function detect_sentencepiece_input_type(input::BioSequences.LongAA)
    return :aa
end

function detect_sentencepiece_input_type(input::AbstractString)
    if isascii(input)
        return :ascii
    else
        return :unicode
    end
end

# ============================================================================
# Model Training Functions
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Train a SentencePiece model on input text/sequences.

# Arguments
- `input_file::String`: Path to input file (one sequence/sentence per line)
- `model_prefix::String`: Output model file prefix (creates .model and .vocab files)
- `vocab_size::Int=8000`: Target vocabulary size
- `model_type::Symbol=:bpe`: Model type - `:bpe`, `:unigram`, `:char`, or `:word`
- `character_coverage::Float64=1.0`: Character coverage (1.0 recommended for biological sequences)
- `input_sentence_size::Int=0`: Maximum number of sentences to use (0 = all)
- `shuffle_input_sentence::Bool=true`: Shuffle input sentences
- `user_defined_symbols::Vector{String}=String[]`: User-defined symbols to include

# Returns
- `NamedTuple` with paths to the generated `.model` and `.vocab` files

# Model Types
- `:bpe` - Byte Pair Encoding (recommended for biological sequences)
- `:unigram` - Unigram language model
- `:char` - Character-level model
- `:word` - Word-level model (requires pre-tokenized input)

# Examples
```julia
result = train_sentencepiece_model(
    input_file="sequences.txt",
    model_prefix="dna_model",
    vocab_size=16000,
    model_type=:bpe,
    character_coverage=1.0
)
println("Model: ", result.model_file)
println("Vocab: ", result.vocab_file)
```
"""
function train_sentencepiece_model(;
    input_file::String,
    model_prefix::String,
    vocab_size::Int=8000,
    model_type::Symbol=:bpe,
    character_coverage::Float64=1.0,
    input_sentence_size::Int=0,
    shuffle_input_sentence::Bool=true,
    user_defined_symbols::Vector{String}=String[]
)
    # Input validation
    @assert isfile(input_file) "Input file does not exist: $input_file"
    @assert vocab_size > 0 "vocab_size must be positive"
    @assert model_type in [:bpe, :unigram, :char, :word] "model_type must be :bpe, :unigram, :char, or :word"
    @assert 0.0 < character_coverage <= 1.0 "character_coverage must be in (0, 1]"

    _ensure_sentencepiece_installed()

    # Map model type to SentencePiece string
    model_type_map = Dict(
        :bpe => "bpe",
        :unigram => "unigram",
        :char => "char",
        :word => "word"
    )

    # Build user symbols argument
    user_symbols_arg = isempty(user_defined_symbols) ? "" : "user_defined_symbols='$(join(user_defined_symbols, ","))',"

    # Escape paths for Python
    input_file_escaped = replace(input_file, "\\" => "\\\\")
    model_prefix_escaped = replace(model_prefix, "\\" => "\\\\")

    python_script = """
import sentencepiece as spm

spm.SentencePieceTrainer.train(
    input='$(input_file_escaped)',
    model_prefix='$(model_prefix_escaped)',
    vocab_size=$(vocab_size),
    model_type='$(model_type_map[model_type])',
    character_coverage=$(character_coverage),
    input_sentence_size=$(input_sentence_size),
    shuffle_input_sentence=$(shuffle_input_sentence ? "True" : "False"),
    $(user_symbols_arg)
)
print('TRAINING_COMPLETED')
"""

    # Write script to temp file and execute
    script_file = tempname() * ".py"
    try
        write(script_file, python_script)
        result = Base.read(`$(Mycelia.CONDA_RUNNER) run -n sentencepiece_env python $(script_file)`, String)

        if !occursin("TRAINING_COMPLETED", result)
            error("SentencePiece training may have failed. Output: $result")
        end
    finally
        isfile(script_file) && rm(script_file)
    end

    model_file = model_prefix * ".model"
    vocab_file = model_prefix * ".vocab"

    @assert isfile(model_file) "Model file was not created: $model_file"
    @assert isfile(vocab_file) "Vocab file was not created: $vocab_file"

    return (model_file=model_file, vocab_file=vocab_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Train a SentencePiece model directly from a vector of sequences.

# Arguments
- `sequences`: Vector of BioSequences or strings to train on
- `model_prefix::String`: Output model file prefix
- `vocab_size::Int=8000`: Target vocabulary size
- `model_type::Symbol=:bpe`: Model type
- `character_coverage::Float64=1.0`: Character coverage
- `kwargs...`: Additional arguments passed to train_sentencepiece_model

# Examples
```julia
dna_seqs = [BioSequences.LongDNA{4}("ACGTACGT"), BioSequences.LongDNA{4}("GCTAGCTA")]
result = train_sentencepiece_model_from_sequences(
    sequences=dna_seqs,
    model_prefix="dna_model",
    vocab_size=100
)
```
"""
function train_sentencepiece_model_from_sequences(;
    sequences::Vector{<:Union{BioSequences.BioSequence, AbstractString}},
    model_prefix::String,
    vocab_size::Int=8000,
    model_type::Symbol=:bpe,
    character_coverage::Float64=1.0,
    kwargs...
)
    @assert !isempty(sequences) "sequences cannot be empty"

    # Convert sequences to strings and write to temp file
    temp_input = tempname() * ".txt"
    try
        open(temp_input, "w") do io
            for seq in sequences
                if seq isa BioSequences.BioSequence
                    println(io, sentencepiece_biosequence_to_string(seq))
                else
                    println(io, string(seq))
                end
            end
        end

        return train_sentencepiece_model(;
            input_file=temp_input,
            model_prefix=model_prefix,
            vocab_size=vocab_size,
            model_type=model_type,
            character_coverage=character_coverage,
            kwargs...
        )
    finally
        isfile(temp_input) && rm(temp_input)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Train a SentencePiece model from sequences in a FASTA file.

# Arguments
- `fasta_file::String`: Path to FASTA file
- `model_prefix::String`: Output model file prefix
- `vocab_size::Int=8000`: Target vocabulary size
- `model_type::Symbol=:bpe`: Model type
- `character_coverage::Float64=1.0`: Character coverage
- `kwargs...`: Additional arguments passed to train_sentencepiece_model
"""
function train_sentencepiece_model_from_fasta(;
    fasta_file::String,
    model_prefix::String,
    vocab_size::Int=8000,
    model_type::Symbol=:bpe,
    character_coverage::Float64=1.0,
    kwargs...
)
    @assert isfile(fasta_file) "FASTA file does not exist: $fasta_file"
    @assert occursin(FASTA_REGEX, fasta_file) "File does not appear to be a FASTA file: $fasta_file"

    # Convert FASTA to temp text file (one sequence per line)
    temp_input = tempname() * ".txt"
    try
        open(temp_input, "w") do io
            for record in open_fastx(fasta_file)
                sequence_str = FASTX.sequence(String, record)
                println(io, uppercase(sequence_str))
            end
        end

        return train_sentencepiece_model(;
            input_file=temp_input,
            model_prefix=model_prefix,
            vocab_size=vocab_size,
            model_type=model_type,
            character_coverage=character_coverage,
            kwargs...
        )
    finally
        isfile(temp_input) && rm(temp_input)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Train a SentencePiece model from sequences in a FASTQ file.

# Arguments
- `fastq_file::String`: Path to FASTQ file
- `model_prefix::String`: Output model file prefix
- `vocab_size::Int=8000`: Target vocabulary size
- `model_type::Symbol=:bpe`: Model type
- `character_coverage::Float64=1.0`: Character coverage
- `kwargs...`: Additional arguments passed to train_sentencepiece_model
"""
function train_sentencepiece_model_from_fastq(;
    fastq_file::String,
    model_prefix::String,
    vocab_size::Int=8000,
    model_type::Symbol=:bpe,
    character_coverage::Float64=1.0,
    kwargs...
)
    @assert isfile(fastq_file) "FASTQ file does not exist: $fastq_file"
    @assert occursin(FASTQ_REGEX, fastq_file) "File does not appear to be a FASTQ file: $fastq_file"

    # Convert FASTQ to temp text file (one sequence per line)
    temp_input = tempname() * ".txt"
    try
        open(temp_input, "w") do io
            for record in open_fastx(fastq_file)
                sequence_str = FASTX.sequence(String, record)
                println(io, uppercase(sequence_str))
            end
        end

        return train_sentencepiece_model(;
            input_file=temp_input,
            model_prefix=model_prefix,
            vocab_size=vocab_size,
            model_type=model_type,
            character_coverage=character_coverage,
            kwargs...
        )
    finally
        isfile(temp_input) && rm(temp_input)
    end
end

# ============================================================================
# Encoding Functions
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Encode text or sequences using a trained SentencePiece model.

# Arguments
- `model_file::String`: Path to the .model file
- `input`: Single sequence/string or vector of sequences/strings to encode
- `output_format::Symbol=:pieces`: Output format - `:pieces` for subword strings, `:ids` for integer IDs
- `enable_sampling::Bool=false`: Enable subword regularization (sampling)
- `alpha::Float64=0.1`: Sampling smoothing parameter (used when enable_sampling=true)
- `nbest_size::Int=-1`: Number of sampling candidates (-1 = all, used when enable_sampling=true)

# Returns
- For single input: Vector of pieces (strings) or IDs (integers)
- For vector input: Vector of vectors (pieces or IDs for each input)

# Examples
```julia
# Single sequence
pieces = encode_sentencepiece(
    model_file="dna_model.model",
    input=BioSequences.LongDNA{4}("ACGTACGT"),
    output_format=:pieces
)

# Multiple sequences
batch_pieces = encode_sentencepiece(
    model_file="dna_model.model",
    input=["ACGT", "GCTA", "TATA"],
    output_format=:pieces
)

# With subword regularization for neural network training
sampled = encode_sentencepiece(
    model_file="dna_model.model",
    input="ACGTACGT",
    output_format=:ids,
    enable_sampling=true,
    alpha=0.1,
    nbest_size=-1
)
```
"""
function encode_sentencepiece(;
    model_file::String,
    input::Union{AbstractString, BioSequences.BioSequence, Vector{<:Union{AbstractString, BioSequences.BioSequence}}},
    output_format::Symbol=:pieces,
    enable_sampling::Bool=false,
    alpha::Float64=0.1,
    nbest_size::Int=-1
)
    @assert isfile(model_file) "Model file does not exist: $model_file"
    @assert output_format in [:pieces, :ids] "output_format must be :pieces or :ids"

    _ensure_sentencepiece_installed()

    model_file_escaped = replace(model_file, "\\" => "\\\\")

    # Handle batch vs single input
    if input isa Vector
        @assert !isempty(input) "input cannot be empty"
        return _encode_sentencepiece_batch(model_file_escaped, input, output_format, enable_sampling, alpha, nbest_size)
    else
        return _encode_sentencepiece_single(model_file_escaped, input, output_format, enable_sampling, alpha, nbest_size)
    end
end

# Internal helper for single input encoding
function _encode_sentencepiece_single(model_file_escaped::String, input, output_format::Symbol,
                                       enable_sampling::Bool, alpha::Float64, nbest_size::Int)
    # Convert BioSequence to string if needed
    input_str = input isa BioSequences.BioSequence ? sentencepiece_biosequence_to_string(input) : string(input)

    # Escape special characters for Python string
    input_escaped = replace(input_str, "\\" => "\\\\")
    input_escaped = replace(input_escaped, "'" => "\\'")
    input_escaped = replace(input_escaped, "\n" => "\\n")

    if enable_sampling
        encode_call = output_format == :pieces ?
            "sp.encode('$(input_escaped)', out_type=str, enable_sampling=True, alpha=$(alpha), nbest_size=$(nbest_size))" :
            "sp.encode('$(input_escaped)', out_type=int, enable_sampling=True, alpha=$(alpha), nbest_size=$(nbest_size))"
    else
        encode_call = output_format == :pieces ?
            "sp.encode_as_pieces('$(input_escaped)')" :
            "sp.encode_as_ids('$(input_escaped)')"
    end

    python_script = """
import sentencepiece as spm
import json

sp = spm.SentencePieceProcessor(model_file='$(model_file_escaped)')
result = $(encode_call)
print(json.dumps(result))
"""

    script_file = tempname() * ".py"
    try
        write(script_file, python_script)
        result_str = Base.read(`$(Mycelia.CONDA_RUNNER) run -n sentencepiece_env python $(script_file)`, String)
        result = JSON.parse(strip(result_str))

        if output_format == :pieces
            return Vector{String}(result)
        else
            return Vector{Int}(result)
        end
    finally
        isfile(script_file) && rm(script_file)
    end
end

# Internal helper for batch encoding
function _encode_sentencepiece_batch(model_file_escaped::String, input::Vector, output_format::Symbol,
                                      enable_sampling::Bool, alpha::Float64, nbest_size::Int)
    # Convert all inputs to strings
    input_strings = [s isa BioSequences.BioSequence ? sentencepiece_biosequence_to_string(s) : string(s) for s in input]

    # Write inputs to temp file (JSON array)
    temp_input = tempname() * ".json"

    try
        write(temp_input, JSON.json(input_strings))

        if enable_sampling
            encode_line = output_format == :pieces ?
                "results.append(sp.encode(text, out_type=str, enable_sampling=True, alpha=$(alpha), nbest_size=$(nbest_size)))" :
                "results.append(sp.encode(text, out_type=int, enable_sampling=True, alpha=$(alpha), nbest_size=$(nbest_size)))"
        else
            encode_line = output_format == :pieces ?
                "results.append(sp.encode_as_pieces(text))" :
                "results.append(sp.encode_as_ids(text))"
        end

        python_script = """
import sentencepiece as spm
import json

sp = spm.SentencePieceProcessor(model_file='$(model_file_escaped)')

with open('$(replace(temp_input, "\\" => "\\\\"))', 'r') as f:
    inputs = json.load(f)

results = []
for text in inputs:
    $(encode_line)
print(json.dumps(results))
"""

        script_file = tempname() * ".py"
        try
            write(script_file, python_script)
            result_str = Base.read(`$(Mycelia.CONDA_RUNNER) run -n sentencepiece_env python $(script_file)`, String)
            results = JSON.parse(strip(result_str))

            if output_format == :pieces
                return [Vector{String}(r) for r in results]
            else
                return [Vector{Int}(r) for r in results]
            end
        finally
            isfile(script_file) && rm(script_file)
        end
    finally
        isfile(temp_input) && rm(temp_input)
    end
end

# ============================================================================
# Decoding Functions
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Decode SentencePiece tokens back to text.

# Arguments
- `model_file::String`: Path to the .model file
- `input`: Single vector of pieces/IDs or vector of vectors

# Returns
- For single vector input: Decoded string
- For vector of vectors: Vector of decoded strings

# Examples
```julia
# Decode pieces
text = decode_sentencepiece(
    model_file="dna_model.model",
    input=["AC", "GT", "AC", "GT"]
)

# Decode IDs
text = decode_sentencepiece(
    model_file="dna_model.model",
    input=[1, 2, 3, 4]
)

# Batch decode
texts = decode_sentencepiece(
    model_file="dna_model.model",
    input=[["AC", "GT"], ["GC", "TA"]]
)
```
"""
function decode_sentencepiece(;
    model_file::String,
    input::Union{Vector{String}, Vector{Int}, Vector{Vector{String}}, Vector{Vector{Int}}}
)
    @assert isfile(model_file) "Model file does not exist: $model_file"
    @assert !isempty(input) "input cannot be empty"

    _ensure_sentencepiece_installed()

    model_file_escaped = replace(model_file, "\\" => "\\\\")

    # Dispatch based on input type
    if input isa Vector{String}
        return _decode_sentencepiece_pieces(model_file_escaped, input)
    elseif input isa Vector{Int}
        return _decode_sentencepiece_ids(model_file_escaped, input)
    elseif input isa Vector{Vector{String}}
        return [_decode_sentencepiece_pieces(model_file_escaped, v) for v in input]
    elseif input isa Vector{Vector{Int}}
        return [_decode_sentencepiece_ids(model_file_escaped, v) for v in input]
    else
        error("Unsupported input type for decode_sentencepiece")
    end
end

# Internal helper for decoding pieces
function _decode_sentencepiece_pieces(model_file_escaped::String, input::Vector{String})
    input_json = JSON.json(input)
    input_json_escaped = replace(input_json, "\\" => "\\\\")
    input_json_escaped = replace(input_json_escaped, "'" => "\\'")

    python_script = """
import sentencepiece as spm
import json

sp = spm.SentencePieceProcessor(model_file='$(model_file_escaped)')
pieces = json.loads('$(input_json_escaped)')
result = sp.decode_pieces(pieces)
print(result)
"""

    script_file = tempname() * ".py"
    try
        write(script_file, python_script)
        result = Base.read(`$(Mycelia.CONDA_RUNNER) run -n sentencepiece_env python $(script_file)`, String)
        return strip(result)
    finally
        isfile(script_file) && rm(script_file)
    end
end

# Internal helper for decoding IDs
function _decode_sentencepiece_ids(model_file_escaped::String, input::Vector{Int})
    input_json = JSON.json(input)

    python_script = """
import sentencepiece as spm
import json

sp = spm.SentencePieceProcessor(model_file='$(model_file_escaped)')
ids = json.loads('$(input_json)')
result = sp.decode_ids(ids)
print(result)
"""

    script_file = tempname() * ".py"
    try
        write(script_file, python_script)
        result = Base.read(`$(Mycelia.CONDA_RUNNER) run -n sentencepiece_env python $(script_file)`, String)
        return strip(result)
    finally
        isfile(script_file) && rm(script_file)
    end
end

# ============================================================================
# Model Loading & Utility Functions
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Check if a file is a valid SentencePiece model.

# Arguments
- `model_file::String`: Path to the potential .model file

# Returns
- `true` if the model can be loaded successfully, `false` otherwise
"""
function is_valid_sentencepiece_model(model_file::String)::Bool
    if !isfile(model_file)
        return false
    end

    _ensure_sentencepiece_installed()

    model_file_escaped = replace(model_file, "\\" => "\\\\")

    python_script = """
import sentencepiece as spm

try:
    sp = spm.SentencePieceProcessor(model_file='$(model_file_escaped)')
    print('VALID')
except Exception as e:
    print('INVALID')
"""

    script_file = tempname() * ".py"
    try
        write(script_file, python_script)
        result = Base.read(`$(Mycelia.CONDA_RUNNER) run -n sentencepiece_env python $(script_file)`, String)
        return occursin("VALID", result) && !occursin("INVALID", result)
    catch
        return false
    finally
        isfile(script_file) && rm(script_file)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load a pre-trained SentencePiece model and return its information.

# Arguments
- `model_file::String`: Path to the .model file

# Returns
- `NamedTuple` with `:model_file` and `:vocab_size`

# Examples
```julia
model_info = load_sentencepiece_model("pretrained.model")
println("Vocab size: ", model_info.vocab_size)
```
"""
function load_sentencepiece_model(model_file::String)
    @assert isfile(model_file) "Model file does not exist: $model_file"

    if !is_valid_sentencepiece_model(model_file)
        error("Invalid SentencePiece model file: $model_file")
    end

    vocab_size = sentencepiece_vocab_size(model_file)

    return (model_file=model_file, vocab_size=vocab_size)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the vocabulary size of a trained SentencePiece model.

# Arguments
- `model_file::String`: Path to the .model file

# Returns
- Integer vocabulary size
"""
function sentencepiece_vocab_size(model_file::String)::Int
    @assert isfile(model_file) "Model file does not exist: $model_file"

    _ensure_sentencepiece_installed()

    model_file_escaped = replace(model_file, "\\" => "\\\\")

    python_script = """
import sentencepiece as spm

sp = spm.SentencePieceProcessor(model_file='$(model_file_escaped)')
print(sp.get_piece_size())
"""

    script_file = tempname() * ".py"
    try
        write(script_file, python_script)
        result = Base.read(`$(Mycelia.CONDA_RUNNER) run -n sentencepiece_env python $(script_file)`, String)
        return parse(Int, strip(result))
    finally
        isfile(script_file) && rm(script_file)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the vocabulary from a trained SentencePiece model as a DataFrame.

# Arguments
- `model_file::String`: Path to the .model file

# Returns
- `DataFrames.DataFrame` with columns: `id`, `piece`, `score`
"""
function get_sentencepiece_vocab(model_file::String)
    @assert isfile(model_file) "Model file does not exist: $model_file"

    _ensure_sentencepiece_installed()

    model_file_escaped = replace(model_file, "\\" => "\\\\")

    python_script = """
import sentencepiece as spm
import json

sp = spm.SentencePieceProcessor(model_file='$(model_file_escaped)')

vocab = []
for i in range(sp.get_piece_size()):
    vocab.append({
        'id': i,
        'piece': sp.id_to_piece(i),
        'score': sp.get_score(i)
    })
print(json.dumps(vocab))
"""

    script_file = tempname() * ".py"
    try
        write(script_file, python_script)
        result_str = Base.read(`$(Mycelia.CONDA_RUNNER) run -n sentencepiece_env python $(script_file)`, String)
        vocab_data = JSON.parse(strip(result_str))

        return DataFrames.DataFrame(
            id = [v["id"] for v in vocab_data],
            piece = [v["piece"] for v in vocab_data],
            score = [v["score"] for v in vocab_data]
        )
    finally
        isfile(script_file) && rm(script_file)
    end
end
