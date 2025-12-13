# Advanced Usage Examples

This page provides advanced usage patterns and optimization techniques for Mycelia.

## Sequence Tokenization with SentencePiece

SentencePiece provides subword tokenization for biological sequences, enabling machine learning and NLP applications on DNA, RNA, amino acid sequences, and text.

### Training a Model on DNA Sequences

```julia
import Mycelia
import BioSequences

# Create training sequences
dna_sequences = [
    BioSequences.LongDNA{4}("ACGTACGTACGTACGTACGT"),
    BioSequences.LongDNA{4}("GCTAGCTAGCTAGCTAGCTA"),
    BioSequences.LongDNA{4}("ATATATATATATATATATATAT"),
    BioSequences.LongDNA{4}("GCGCGCGCGCGCGCGCGCGC"),
]

# Train a BPE model with vocabulary size 100
result = Mycelia.train_sentencepiece_model_from_sequences(
    sequences=dna_sequences,
    model_prefix="dna_tokenizer",
    vocab_size=100,
    model_type=:bpe,           # or :unigram, :char, :word
    character_coverage=1.0     # 1.0 recommended for biological sequences
)

println("Model saved to: ", result.model_file)
println("Vocab saved to: ", result.vocab_file)
```

### Training from FASTA/FASTQ Files

```julia
# Train from FASTA file
result = Mycelia.train_sentencepiece_model_from_fasta(
    fasta_file="sequences.fasta",
    model_prefix="genome_tokenizer",
    vocab_size=16000,
    model_type=:bpe
)

# Train from FASTQ file
result = Mycelia.train_sentencepiece_model_from_fastq(
    fastq_file="reads.fastq",
    model_prefix="reads_tokenizer",
    vocab_size=8000
)
```

### Encoding Sequences

```julia
# Encode a single sequence to subword pieces
pieces = Mycelia.encode_sentencepiece(
    model_file="dna_tokenizer.model",
    input="ACGTACGTACGT",
    output_format=:pieces
)
# Returns: ["AC", "GT", "AC", "GT", "AC", "GT"]

# Encode to integer IDs (for neural networks)
ids = Mycelia.encode_sentencepiece(
    model_file="dna_tokenizer.model",
    input="ACGTACGTACGT",
    output_format=:ids
)
# Returns: [5, 12, 5, 12, 5, 12]

# Encode BioSequences directly
dna = BioSequences.LongDNA{4}("ACGTACGT")
pieces = Mycelia.encode_sentencepiece(
    model_file="dna_tokenizer.model",
    input=dna,
    output_format=:pieces
)

# Batch encoding
sequences = ["ACGT", "GCTA", "TATA"]
batch_pieces = Mycelia.encode_sentencepiece(
    model_file="dna_tokenizer.model",
    input=sequences,
    output_format=:pieces
)
```

### Subword Regularization for Neural Network Training

Enable sampling for data augmentation during training:

```julia
# Each call may produce different tokenizations
for i in 1:5
    sampled = Mycelia.encode_sentencepiece(
        model_file="dna_tokenizer.model",
        input="ACGTACGT",
        output_format=:pieces,
        enable_sampling=true,
        alpha=0.1,          # Smoothing parameter
        nbest_size=-1       # Sample from all candidates
    )
    println("Sample $i: ", sampled)
end
```

### Decoding Back to Sequences

```julia
# Decode pieces back to text
pieces = ["AC", "GT", "AC", "GT"]
text = Mycelia.decode_sentencepiece(
    model_file="dna_tokenizer.model",
    input=pieces
)
# Returns: "ACGTACGT"

# Decode IDs
ids = [5, 12, 5, 12]
text = Mycelia.decode_sentencepiece(
    model_file="dna_tokenizer.model",
    input=ids
)

# Batch decode
batch_pieces = [["AC", "GT"], ["GC", "TA"]]
texts = Mycelia.decode_sentencepiece(
    model_file="dna_tokenizer.model",
    input=batch_pieces
)
```

### Working with Pre-trained Models

```julia
# Load and validate a model
model_info = Mycelia.load_sentencepiece_model("pretrained.model")
println("Vocabulary size: ", model_info.vocab_size)

# Check if a file is a valid model
if Mycelia.is_valid_sentencepiece_model("unknown.model")
    println("Valid SentencePiece model")
end

# Get vocabulary as DataFrame
vocab_df = Mycelia.get_sentencepiece_vocab("dna_tokenizer.model")
# DataFrame with columns: id, piece, score
```

### Amino Acid and RNA Sequences

```julia
# RNA sequences
rna_seqs = [BioSequences.LongRNA{4}("ACGUACGUACGU") for _ in 1:50]
result = Mycelia.train_sentencepiece_model_from_sequences(
    sequences=rna_seqs,
    model_prefix="rna_tokenizer",
    vocab_size=50
)

# Amino acid sequences
aa_seqs = [BioSequences.LongAA("ARNDCEQGHILKMFPSTWYV") for _ in 1:50]
result = Mycelia.train_sentencepiece_model_from_sequences(
    sequences=aa_seqs,
    model_prefix="protein_tokenizer",
    vocab_size=200
)
```

### Unicode Text (Non-Biological)

```julia
# Train on general text
texts = [
    "The quick brown fox jumps over the lazy dog.",
    "SentencePiece is an unsupervised text tokenizer.",
    "It supports BPE and unigram language models."
]
result = Mycelia.train_sentencepiece_model_from_sequences(
    sequences=texts,
    model_prefix="text_tokenizer",
    vocab_size=500,
    model_type=:unigram
)
```

## Performance Optimization

Tips and techniques for optimizing bioinformatics workflows.

### Memory Management

- Efficient data structures
- Streaming processing
- Memory estimation tools

### Parallel Processing

- Multi-threading strategies
- Distributed computing
- Batch processing

## Advanced Workflows

Complex workflow patterns combining multiple tools.

### Custom Pipelines

- Building custom analysis pipelines
- Workflow composition
- Error handling and recovery

## Related Documentation

- See tutorials for workflow examples
- See [API Reference](../../api-reference.md) for detailed function documentation
