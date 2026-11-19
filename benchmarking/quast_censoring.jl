# Censoring-aware QUAST parsing and NGA50 status classification, shared by
# ont_k_sweep.jl and the D5 grid driver (td-yddpr.5). Moved verbatim from
# ont_k_sweep.jl so both harnesses classify censoring identically; this file
# parses no ARGS and has no side effects, so any harness can include() it.
#
# Requires nothing beyond Base.

# === QUAST report parsing ===
#
# Returns `missing` — never 0.0 — for any metric QUAST could not compute. This
# is the whole point of the harness: the pilot's coercion of "-" to 0.0 is what
# made a censored floor indistinguishable from a measured zero.

function parse_quast_metrics(report_tsv)
    isfile(report_tsv) && return _parse_quast_metrics(report_tsv)
    return empty_metrics()
end

function _parse_quast_metrics(report_tsv)
    values = Dict{String, Union{Missing, Float64}}()
    open(report_tsv, "r") do io
        for line in eachline(io)
            fields = split(line, '\t')
            length(fields) >= 2 || continue
            label = strip(fields[1])
            parsed = tryparse(Float64, strip(fields[2]))
            values[label] = parsed === nothing ? missing : parsed
        end
    end
    get_metric(name) = get(values, name, missing)
    return (
        quast_contigs = get_metric("# contigs"),
        total_length = get_metric("Total length"),
        N50 = get_metric("N50"),
        largest_contig = get_metric("Largest contig"),
        NGA50 = get_metric("NGA50"),
        NA50 = get_metric("NA50"),
        largest_alignment = get_metric("Largest alignment"),
        genome_fraction = get_metric("Genome fraction (%)"),
        duplication_ratio = get_metric("Duplication ratio"),
        misassemblies = get_metric("# misassemblies"),
        unaligned_contigs = get_metric("# unaligned contigs"),
        unaligned_length = get_metric("Unaligned length")
    )
end

function empty_metrics()
    return (
        quast_contigs = missing, total_length = missing, N50 = missing,
        largest_contig = missing, NGA50 = missing, NA50 = missing,
        largest_alignment = missing, genome_fraction = missing,
        duplication_ratio = missing, misassemblies = missing,
        unaligned_contigs = missing, unaligned_length = missing
    )
end

"""
    nga50_status_for(metrics, n_contigs, contigs_ge_min, quast_ran) -> String

Why NGA50 is absent, when it is. Every one of these is a CENSORED FLOOR, not a
measured zero, and they are kept distinct because they have different causes and
different implications:

  * `measured`            — QUAST computed NGA50; the number is real.
  * `no_contigs`          — the assembler emitted nothing at all.
  * `no_contigs_ge_min`   — contigs exist but NONE reaches `--min-contig`
                            (500 bp), so QUAST declines to score the assembly.
                            This is a RESULT, not a tool failure: it is the
                            signature of the low-k regime, where the graph is so
                            tangled that unitigs never grow past a few tens of
                            bases. The `asm_*` columns quantify it.
  * `censored_no_alignment`  — contigs long enough to score exist, but QUAST
                            reported NO genome fraction at all, meaning nothing
                            aligned to the reference. `# misassemblies` does NOT
                            qualify this — QUAST reports 0 misassemblies on
                            unaligned contigs, so a 0 there means non-alignment,
                            not correctness.
  * `censored_partial_alignment` — contigs ALIGNED (genome fraction is a real,
                            nonzero measurement) but NGA50 is still absent.
                            NGA50 is an NG-family statistic: it exists only once
                            aligned blocks of a given length or longer TOTAL at
                            least half the REFERENCE length. Below that the
                            statistic has no value and QUAST prints "-". This is
                            a definitional boundary, NOT an assembly that failed
                            to align, and conflating the two understates how
                            much of the genome was recovered.

                            NOTE the threshold is on the SUM OF ALIGNED BLOCK
                            LENGTHS, which counts duplicated and overlapping
                            alignments — it is NOT on genome fraction, which
                            measures UNIQUE reference coverage. The two diverge
                            whenever duplication ratio > 1, so a cell can carry
                            a defined NGA50 at well under 50% genome fraction.
                            An earlier version of this label was named
                            `censored_partial_alignment` and its docstring asserted the
                            boundary was on genome fraction; the alignment
                            threshold diagnostic in this same directory refutes
                            that directly — Lambda/ONT/k=31/30x/seed123 rescored
                            at 90% identity has genome fraction 23.762% and a
                            DEFINED NGA50 of 508. Six such rows exist. The label
                            was renamed because it asserted a condition this
                            function never tests: the check below is simply
                            "NGA50 absent, genome fraction present".
  * `quast_failed`        — QUAST genuinely errored on an assembly it should
                            have been able to score. This is the ONLY value
                            here that indicates an infrastructure problem, which
                            is why `no_contigs_ge_min` was split out of it.

Because NGA50 is a step function of genome fraction at the 50% boundary, it is
an UNSTABLE endpoint wherever cells sit near that line: two replicate seeds of
the same condition can land on opposite sides and report "measured" versus
"absent" for what is nearly the same assembly. Genome fraction is continuous
across that boundary and does not have the pathology, which is why every
summary in this harness reports it beside NGA50.

`contigs_ge_min` is computed in Julia from the contigs themselves, so this
classification never depends on QUAST having succeeded.
"""
function nga50_status_for(metrics, n_contigs::Int, contigs_ge_min::Int,
        quast_ran::Bool)
    n_contigs == 0 && return "no_contigs"
    contigs_ge_min == 0 && return "no_contigs_ge_min"
    quast_ran || return "quast_failed"
    ismissing(metrics.NGA50) || return "measured"
    # NGA50 absent: distinguish "nothing aligned" from "aligned, but under the
    # 50% genome-fraction floor below which NGA50 has no definition".
    ismissing(metrics.genome_fraction) && return "censored_no_alignment"
    return metrics.genome_fraction > 0 ? "censored_partial_alignment" :
           "censored_no_alignment"
end

"""
    contig_stats(contigs, min_contig) -> (; n, n_ge_min, max_length, total_bp)

Assembler-side contig statistics, independent of QUAST.

These are the numbers that survive when QUAST cannot score a cell, and they are
what makes the low-k regime legible: `n` large with `max_length` small is a
tangled graph, whereas `n` large with `max_length` large is a graph that
assembled long but erroneous paths. Those are opposite failure modes and the
sweep has to be able to tell them apart.
"""
function contig_stats(contigs, min_contig::Int)
    lengths = [length(c) for c in contigs]
    return (
        n = length(lengths),
        n_ge_min = count(>=(min_contig), lengths),
        max_length = isempty(lengths) ? 0 : maximum(lengths),
        total_bp = isempty(lengths) ? 0 : sum(lengths)
    )
end
