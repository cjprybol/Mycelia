# ==============================================================================
# Input→workflow assembly router (Stage 0) — the capstone of the à la carte
# k-mer classification ecosystem.
#
# Premise (project author): no single k-mer classification strategy is best for
# every sequencing input. The right arm depends on the input library's
# characteristics — coverage depth, error rate, and the shape of the k-mer
# spectrum. This module CHARACTERIZES a read set and RECOMMENDS a classifier arm.
#
# This is a TRANSPARENT v0 HEURISTIC, deliberately rule-based and fully
# documented, so its decisions are auditable. It is explicitly intended to be
# REPLACED by a learned successor. The evidence base for that successor is the
# comparison matrix produced by running every `KmerClassifier` arm against a
# known-truth error set (the per-k-mer training table described in
# `kmer-classification.jl`): each arm's ROC/PR performance across a grid of
# simulated inputs becomes the labeled table on which a model can learn the
# input-features → best-arm mapping that these hand-written rules approximate.
#
# Nothing here touches the graph's decision rules — it reads the SAME two
# evidence signals the arms consume (coverage from k-mer multiplicity, quality
# from Phred), one abstraction level up: features of the whole library rather
# than of a single k-mer.
# ==============================================================================

"""
    extract_input_features(records; k) -> NamedTuple

Characterize a FASTQ read set with the library-level features a router uses to
pick a k-mer classification strategy. Computed from the records and their
(canonical) k-mer multiplicities.

# Arguments
- `records::Vector{<:FASTX.FASTQ.Record}`: the input read set.
- `k::Integer`: k-mer length used to build the multiplicity spectrum.

# Returns a `NamedTuple` with
- `estimated_coverage_depth::Float64`: mean canonical-k-mer multiplicity
  (total k-mer observations / number of distinct k-mers). This is a
  reference-free proxy for sequencing depth: for a genome of `G` distinct
  k-mers, total observations ≈ `num_reads · (mean_read_length − k + 1)`, so the
  mean multiplicity ≈ `num_reads · mean_read_length / G` — the usual
  depth = reads · read-length / genome-size relation with `G` standing in for
  genome size. `0.0` if no k-mers were observed.
- `mean_read_quality::Float64`: mean Phred score across every base of every
  read (`0.0` if there are no bases). A proxy for the library's error rate.
- `kmer_spectrum_valley::Union{Int, Missing}`: multiplicity at the first local
  minimum of the k-mer multiplicity histogram — the valley between the
  low-multiplicity error peak and the genomic peak in a well-covered bimodal
  spectrum, and a natural fixed cutoff when present. `missing` when the spectrum
  has no clear valley (typical for low coverage or a unimodal spectrum).
- `unique_to_total_kmer_ratio::Float64`: distinct k-mers / total k-mer
  observations. Near `1.0` ⇒ little redundancy (shallow / error-dominated);
  small ⇒ each k-mer seen many times (deep coverage). It is the reciprocal of
  `estimated_coverage_depth`. `0.0` if no k-mers were observed.
- `num_reads::Int`: number of records.
- `mean_read_length::Float64`: mean sequence length across records (`0.0` if
  empty).

The k-mer multiplicity spectrum is built over canonical DNA k-mers via
`Kmers.UnambiguousDNAMers{k}` (ambiguous positions skipped), matching the
`:canonical` qualmer graph the classifier arms are scored on.
"""
function extract_input_features(records::Vector{<:FASTX.FASTQ.Record}; k::Integer)
    num_reads = length(records)

    # --- read-level features (num_reads, mean_read_length, mean_read_quality) ---
    total_read_length = 0
    total_quality = 0.0
    total_bases = 0
    for record in records
        len = length(FASTX.sequence(record))
        total_read_length += len
        for q in FASTX.quality_scores(record)
            total_quality += Float64(q)
            total_bases += 1
        end
    end
    mean_read_length = num_reads > 0 ? total_read_length / num_reads : 0.0
    mean_read_quality = total_bases > 0 ? total_quality / total_bases : 0.0

    # --- k-mer multiplicity spectrum (canonical DNA k-mers) ---
    kmer_counts = Dict{Kmers.DNAKmer{k}, Int}()
    for record in records
        sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)
        for (kmer, _) in Kmers.UnambiguousDNAMers{k}(sequence)
            canonical_kmer = BioSequences.canonical(kmer)
            kmer_counts[canonical_kmer] = get(kmer_counts, canonical_kmer, 0) + 1
        end
    end
    distinct_kmers = length(kmer_counts)
    total_kmer_observations = sum(values(kmer_counts); init = 0)

    estimated_coverage_depth = distinct_kmers > 0 ?
        total_kmer_observations / distinct_kmers : 0.0
    unique_to_total_kmer_ratio = total_kmer_observations > 0 ?
        distinct_kmers / total_kmer_observations : 0.0
    kmer_spectrum_valley = _spectrum_valley(kmer_counts)

    return (
        estimated_coverage_depth = estimated_coverage_depth,
        mean_read_quality = mean_read_quality,
        kmer_spectrum_valley = kmer_spectrum_valley,
        unique_to_total_kmer_ratio = unique_to_total_kmer_ratio,
        num_reads = num_reads,
        mean_read_length = mean_read_length,
    )
end

"""
    _spectrum_valley(kmer_counts) -> Union{Int, Missing}

First local minimum of the k-mer multiplicity histogram. Scans multiplicities
`m = 2, 3, …` and returns the first `m` where the number of distinct k-mers at
`m` is `≤` its value at `m − 1` and strictly `<` its value at `m + 1` (a
descending-then-ascending valley between the error peak near `m = 1` and the
genomic peak). Returns `missing` when the histogram is too short or has no such
valley.
"""
function _spectrum_valley(kmer_counts::Dict{<:Any, Int})::Union{Int, Missing}
    isempty(kmer_counts) && return missing
    max_multiplicity = maximum(values(kmer_counts))
    max_multiplicity < 3 && return missing
    # histogram[m] = number of distinct k-mers observed exactly m times
    histogram = zeros(Int, max_multiplicity)
    for count in values(kmer_counts)
        histogram[count] += 1
    end
    for m in 2:(max_multiplicity - 1)
        if histogram[m] <= histogram[m - 1] && histogram[m] < histogram[m + 1]
            return m
        end
    end
    return missing
end

"""
    select_classifier(features) -> NamedTuple

Recommend a `KmerClassifier` arm for a read set given its
[`extract_input_features`](@ref) output. Returns
`(classifier::KmerClassifier, rationale::String)`.

HEURISTIC v0 — transparent, hand-written rules meant to be replaced by a model
trained on the comparison matrix (the per-k-mer training table in
`kmer-classification.jl`, where every arm is scored by ROC/PR against a
known-truth error set across a grid of simulated libraries). The learned
successor will map the same input features to the empirically best arm; these
rules encode the current best-guess of that mapping.

Decision rules (first match wins):
1. **Very low coverage** (`estimated_coverage_depth < 5`): a fixed count cutoff
   would discard real low-coverage k-mers, so recommend the quality-aware fused
   [`BayesianMixtureClassifier`], where strong quality can rescue a
   low-coverage real k-mer. `genomic_mean_coverage` is set to the estimated
   depth (floored at 2.0 so the Poisson component is well-defined).
2. **Elevated error rate** (`mean_read_quality < 20`, i.e. > 1% expected error):
   favor fusion that lets quality DISCOUNT coverage —
   [`EffectiveCoverageClassifier`] — so error-inflated coverage is shrunk before
   the likelihood ratio.
3. **Ample coverage with a clear spectrum valley** (`kmer_spectrum_valley`
   present): the spectrum is bimodal and the valley IS the natural cutoff, so
   the coverage-only [`FixedCoverageThreshold`] at `min_count = valley` is both
   sufficient and the most interpretable.
4. **Otherwise** (moderate coverage, decent quality, no clean valley): default
   to the principled fused [`BayesianMixtureClassifier`] with
   `genomic_mean_coverage` = estimated depth.
"""
function select_classifier(features::NamedTuple)
    depth = features.estimated_coverage_depth
    quality = features.mean_read_quality
    valley = features.kmer_spectrum_valley

    # Rule 1 — very low coverage: rescue real low-coverage k-mers with quality.
    if depth < 5.0
        genomic_mean_coverage = max(depth, 2.0)
        classifier = BayesianMixtureClassifier(;
            genomic_mean_coverage = genomic_mean_coverage,
        )
        rationale = "very low coverage (depth≈$(round(depth; digits = 2))); " *
            "a fixed count cutoff would discard real low-coverage k-mers, so " *
            "fuse quality with coverage to rescue them"
        return (classifier = classifier, rationale = rationale)
    end

    # Rule 2 — elevated error rate: let quality discount coverage.
    if quality < 20.0
        genomic_mean_coverage = max(depth, 2.0)
        classifier = EffectiveCoverageClassifier(;
            genomic_mean_coverage = genomic_mean_coverage,
        )
        rationale = "elevated error rate (mean Phred≈$(round(quality; digits = 1)) " *
            "⇒ >1% expected error); discount coverage by quality before the " *
            "likelihood ratio"
        return (classifier = classifier, rationale = rationale)
    end

    # Rule 3 — ample coverage with a clear bimodal valley: the valley is the cutoff.
    if valley !== missing
        classifier = FixedCoverageThreshold(valley)
        rationale = "ample coverage (depth≈$(round(depth; digits = 2))) with a " *
            "clear k-mer-spectrum valley at multiplicity=$(valley); a fixed " *
            "coverage cutoff at the valley cleanly separates error from genomic"
        return (classifier = classifier, rationale = rationale)
    end

    # Rule 4 — default: principled fused arm.
    genomic_mean_coverage = max(depth, 2.0)
    classifier = BayesianMixtureClassifier(;
        genomic_mean_coverage = genomic_mean_coverage,
    )
    rationale = "moderate coverage (depth≈$(round(depth; digits = 2))), decent " *
        "quality, no clean spectrum valley; principled coverage+quality fusion " *
        "as the default"
    return (classifier = classifier, rationale = rationale)
end
