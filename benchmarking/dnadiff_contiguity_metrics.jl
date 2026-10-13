# dnadiff-derived contiguity metrics: misassembly proxy + NGA50.
# ==============================================================
#
# Pure, dependency-free helpers (like quast_report_parsing.jl and
# rhizomorph_correction_accuracy_metrics.jl): no Mycelia, no network, no
# external tools — just arithmetic and text parsing over MUMmer dnadiff output.
# Kept dependency-free so the unit test exercises the logic in milliseconds
# without loading Mycelia or running dnadiff.
#
# WHY dnadiff and not QUAST: QUAST is the canonical source of NGA50 and
# misassembly counts, but its bioconda build does not solve on osx-arm64, so the
# correction-validation harness scores with MUMmer dnadiff (which does). These
# helpers recover the two contiguity signals the SOTA-parity comparison needs —
# a misassembly PROXY and NGA50 — from dnadiff's .report and .1coords outputs.
#
# The misassembly PROXY (relocations + translocations + inversions) is a
# deliberately-named proxy, NOT a QUAST-identical misassembly count: dnadiff's
# feature estimates are computed differently from QUAST's contig-break analysis.
# It is a matched, self-consistent signal for the same-reference arm-vs-arm
# comparison the harness makes; it is not cross-comparable to published QUAST
# misassembly numbers.

"""
    parse_dnadiff_feature_estimates(report) -> NamedTuple

Parse the REF column of the `[Feature Estimates]` block of a MUMmer dnadiff
`.report` into `(breakpoints, relocations, translocations, inversions)`.

Columns in the report are REF then QRY; we take REF (the reference-anchored
count). Throws if the `[Feature Estimates]` header or any of the four fields is
absent — a report missing the block is a malformed/failed run and must fail
closed rather than be read as zero misassemblies.
"""
function parse_dnadiff_feature_estimates(report::AbstractString)::NamedTuple
    occursin("[Feature Estimates]", report) ||
        throw(ArgumentError("dnadiff report has no [Feature Estimates] block"))
    function ref_count(field::AbstractString)::Int
        m = match(Regex("^" * field * raw"\s+(\d+)\s+(\d+)", "m"), report)
        m === nothing &&
            throw(ArgumentError("dnadiff [Feature Estimates] missing field: $(field)"))
        return parse(Int, m.captures[1])
    end
    return (
        breakpoints = ref_count("Breakpoints"),
        relocations = ref_count("Relocations"),
        translocations = ref_count("Translocations"),
        inversions = ref_count("Inversions")
    )
end

"""
    misassembly_proxy(feature_estimates) -> Int

The misassembly proxy = relocations + translocations + inversions. Breakpoints
are deliberately excluded: they count alignment endpoints, not misassembly
events (a single clean contig split across the reference origin yields
breakpoints with zero structural misassembly).
"""
function misassembly_proxy(fe)::Int
    return fe.relocations + fe.translocations + fe.inversions
end

"""
    nga50(block_lengths, reference_length) -> Int

NGA50: the length `L` such that reference-aligned blocks of length ≥ `L` cover
at least half of `reference_length`. Unlike N50, the denominator is the
reference genome length (not the assembly length), so a fragmented or partial
assembly is penalized — blocks are the reference-aligned segments from dnadiff's
`.1coords`, and their total can be less than half the reference (a broken
assembly), in which case NGA50 is 0 (QUAST prints "-").

The half-reference threshold is inclusive. `reference_length` must be positive
and every block length non-negative.
"""
function nga50(block_lengths, reference_length::Integer)::Int
    reference_length > 0 ||
        throw(ArgumentError("reference_length must be positive, got $(reference_length)"))
    any(<(0), block_lengths) &&
        throw(ArgumentError("block lengths must be non-negative"))
    isempty(block_lengths) && return 0
    threshold = reference_length / 2
    cumulative = 0
    for len in sort(collect(block_lengths); rev = true)
        cumulative += len
        cumulative >= threshold && return Int(len)
    end
    return 0
end

"""
    aligned_ref_block_lengths(coords_text) -> Vector{Int}

Extract reference-aligned block lengths from MUMmer dnadiff `.1coords` text
(show-coords -rclHT: tab-separated, columns `S1 E1 S2 E2 LEN1 LEN2 %IDY ...`).
Each block length is derived from the reference start/end (columns 1-2) as
`abs(E1 - S1) + 1`, which is robust to reverse-strand reference coordinates and
to any drift in the trailing columns. Lines whose first two fields are not
integers (headers, blanks) are skipped.
"""
function aligned_ref_block_lengths(coords_text::AbstractString)::Vector{Int}
    lengths = Int[]
    for line in eachline(IOBuffer(coords_text))
        fields = split(strip(line))
        length(fields) >= 2 || continue
        s1 = tryparse(Int, fields[1])
        e1 = tryparse(Int, fields[2])
        (s1 === nothing || e1 === nothing) && continue
        push!(lengths, abs(e1 - s1) + 1)
    end
    return lengths
end
