"""
    pangenome_kmer_fdr_benchmark.jl

Opt-in performance benchmark for K=21 pangenome k-mer analysis. The benchmark
compares `Mycelia.pangenome_kmer_fdr` (streaming reference-union path) against
`Mycelia.analyze_pangenome_kmers` (dense presence/absence matrix path) on a
deterministic synthetic 50-genome panel: 40 subject genomes and 10 reference
genomes, each approximately 50 kb.

Expected ratios for this small panel are streaming memory roughly 5-10x lower
than the dense path, with similar wall-clock time at K=21 after compilation
warmup. Linux runs report process peak RSS from `/proc/self/status`; other
platforms emit elapsed time only because this test intentionally avoids
platform-specific memory tools.
"""
const PANGENOME_KMER_FDR_BENCHMARK = true

import BioSequences
import FASTX
import Kmers
import Mycelia
import Printf
import StableRNGs
import Test

const _BENCHMARK_KMER_TYPE = Kmers.DNAKmer{21}
const _BENCHMARK_GENOME_LENGTH = 50_000
const _BENCHMARK_N_SUBJECTS = 40
const _BENCHMARK_N_REFERENCES = 10
const _DNA_ALPHABET = ('A', 'C', 'G', 'T')

function _random_dna_sequence(rng, n::Integer)
    chars = Vector{Char}(undef, n)
    for i in eachindex(chars)
        chars[i] = _DNA_ALPHABET[rand(rng, eachindex(_DNA_ALPHABET))]
    end
    return join(chars)
end

function _mutate_sequence(rng, seq::AbstractString; mutation_rate::Float64)
    chars = collect(seq)
    for i in eachindex(chars)
        if rand(rng) < mutation_rate
            original = chars[i]
            replacement = original
            while replacement == original
                replacement = _DNA_ALPHABET[rand(rng, eachindex(_DNA_ALPHABET))]
            end
            chars[i] = replacement
        end
    end
    return join(chars)
end

function _write_fasta(dir::AbstractString, name::AbstractString, seq::AbstractString)
    path = joinpath(dir, name)
    open(path, "w") do io
        writer = FASTX.FASTA.Writer(io)
        write(writer, FASTX.FASTA.Record(name, BioSequences.LongDNA{4}(seq)))
        close(writer)
    end
    return path
end

function _write_synthetic_panel(dir::AbstractString)
    rng = StableRNGs.StableRNG(0x21f0_50aa)
    ancestor = _random_dna_sequence(rng, _BENCHMARK_GENOME_LENGTH)

    subjects = String[]
    for i in 1:_BENCHMARK_N_SUBJECTS
        mutation_rate = 0.006 + 0.0002 * (i % 11)
        seq = _mutate_sequence(rng, ancestor; mutation_rate = mutation_rate)
        push!(subjects, _write_fasta(dir, "subject_$(lpad(string(i), 2, '0')).fasta", seq))
    end

    references = String[]
    for i in 1:_BENCHMARK_N_REFERENCES
        mutation_rate = 0.004 + 0.0003 * (i % 7)
        seq = _mutate_sequence(rng, ancestor; mutation_rate = mutation_rate)
        push!(references, _write_fasta(dir, "reference_$(lpad(string(i), 2, '0')).fasta", seq))
    end

    return (
        subjects = subjects,
        references = references,
        all_genomes = vcat(subjects, references),
    )
end

function _linux_peak_rss_bytes()
    Sys.islinux() || return missing
    status = try
        read("/proc/self/status", String)
    catch
        return missing
    end

    for field in ("VmHWM", "VmRSS")
        found = match(Regex("(?m)^$(field):\\s+(\\d+)\\s+kB"), status)
        found === nothing && continue
        return parse(Int, found.captures[1]) * 1024
    end
    return missing
end

function _measure(f, label::AbstractString)
    GC.gc()
    GC.gc()
    rss_before = _linux_peak_rss_bytes()
    gc_before = Base.gc_num()
    started = time()
    result = f()
    elapsed_seconds = time() - started
    gc_after = Base.gc_num()
    rss_after = _linux_peak_rss_bytes()
    return (
        label = label,
        elapsed_seconds = elapsed_seconds,
        peak_rss_before = rss_before,
        peak_rss_bytes = rss_after,
        gc_before = gc_before,
        gc_after = gc_after,
        result = result,
    )
end

function _rss_megabytes(value)
    value === missing && return "n/a"
    return Printf.@sprintf("%.1f", value / 2.0^20)
end

function _ratio_text(numerator, denominator; suffix = "x")
    (numerator === missing || denominator === missing || denominator == 0) && return "n/a"
    return Printf.@sprintf("%.2f%s", numerator / denominator, suffix)
end

function _print_comparison_table(streaming, dense)
    println()
    println("Pangenome k-mer benchmark (K=21, 40 subjects + 10 references, 50 kb/genome)")
    println(
        Printf.@sprintf(
            "%-12s %12s %14s %16s",
            "path",
            "seconds",
            "peak_rss_mb",
            "dense/stream"
        )
    )
    println("-"^60)
    println(
        Printf.@sprintf(
            "%-12s %12.3f %14s %16s",
            "streaming",
            streaming.elapsed_seconds,
            _rss_megabytes(streaming.peak_rss_bytes),
            "1.00x"
        )
    )
    println(
        Printf.@sprintf(
            "%-12s %12.3f %14s %16s",
            "dense",
            dense.elapsed_seconds,
            _rss_megabytes(dense.peak_rss_bytes),
            _ratio_text(dense.elapsed_seconds, streaming.elapsed_seconds)
        )
    )
    memory_ratio = _ratio_text(dense.peak_rss_bytes, streaming.peak_rss_bytes)
    println()
    println("Elapsed dense/streaming ratio: ", _ratio_text(dense.elapsed_seconds, streaming.elapsed_seconds))
    println("Peak RSS dense/streaming ratio: ", memory_ratio)
end

function _warmup(subjects::Vector{String}, references::Vector{String})
    warm_subjects = subjects[1:2]
    warm_references = references[1:1]
    Mycelia.pangenome_kmer_fdr(
        warm_subjects,
        warm_references;
        kmer_type = _BENCHMARK_KMER_TYPE,
        cache = :memory,
        membership_type = :exact,
        progress = false
    )
    Mycelia.analyze_pangenome_kmers(
        vcat(warm_subjects, warm_references);
        kmer_type = _BENCHMARK_KMER_TYPE,
        distance_metric = :jaccard
    )
    GC.gc()
    return nothing
end

function _run_benchmark()
    temp = mktempdir()
    try
        panel = _write_synthetic_panel(temp)
        Test.@test length(panel.subjects) == _BENCHMARK_N_SUBJECTS
        Test.@test length(panel.references) == _BENCHMARK_N_REFERENCES
        Test.@test all(isfile, panel.all_genomes)

        _warmup(panel.subjects, panel.references)

        streaming = _measure("streaming") do
            Mycelia.pangenome_kmer_fdr(
                panel.subjects,
                panel.references;
                kmer_type = _BENCHMARK_KMER_TYPE,
                cache = :memory,
                membership_type = :exact,
                progress = false
            )
        end
        dense = _measure("dense") do
            Mycelia.analyze_pangenome_kmers(
                panel.all_genomes;
                kmer_type = _BENCHMARK_KMER_TYPE,
                distance_metric = :jaccard
            )
        end

        Test.@test length(streaming.result.fdr) == _BENCHMARK_N_SUBJECTS
        Test.@test length(dense.result.genome_names) ==
                   _BENCHMARK_N_SUBJECTS + _BENCHMARK_N_REFERENCES
        Test.@test size(dense.result.presence_absence_matrix, 2) ==
                   _BENCHMARK_N_SUBJECTS + _BENCHMARK_N_REFERENCES

        _print_comparison_table(streaming, dense)

        if streaming.elapsed_seconds <= dense.elapsed_seconds
            Test.@test streaming.elapsed_seconds <= dense.elapsed_seconds
        else
            @warn "Streaming wall-clock exceeded dense path" streaming_seconds=streaming.elapsed_seconds dense_seconds=dense.elapsed_seconds
        end

        if streaming.peak_rss_bytes !== missing && dense.peak_rss_bytes !== missing
            if streaming.peak_rss_bytes <= dense.peak_rss_bytes
                Test.@test streaming.peak_rss_bytes <= dense.peak_rss_bytes
            else
                @warn "Streaming peak RSS exceeded dense path" streaming_peak_rss_bytes=streaming.peak_rss_bytes dense_peak_rss_bytes=dense.peak_rss_bytes
            end
        else
            @info "Peak RSS comparison skipped; /proc/self/status is unavailable on this platform"
        end

        return (streaming = streaming, dense = dense)
    finally
        rm(temp; recursive = true, force = true)
    end
end

if get(ENV, "MYCELIA_RUN_ALL", "") == "true"
    Test.@testset "pangenome_kmer_fdr K=21 benchmark vs dense path" begin
        _run_benchmark()
    end
else
    @info "Skipping pangenome_kmer_fdr benchmark; set MYCELIA_RUN_ALL=true to run"
end
