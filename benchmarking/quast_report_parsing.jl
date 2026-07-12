# QUAST report.tsv parsing + per-arm metric attribution.
# =======================================================
#
# Pure, dependency-free helpers (like rhizomorph_scale_guard.jl): no Mycelia,
# no network, no external tools — just text parsing over a QUAST report.tsv.
# Kept dependency-free so the unit test can exercise the parse/attribution
# logic in milliseconds without loading Mycelia or running QUAST.
#
# `parse_quast_metric` was extracted from mode_comparison.jl's inline parser.
# mode_comparison.jl still carries its own copy and does NOT yet include this
# file — consolidating the two is a clean follow-up. For now this file is the
# shared definition used by the correction-validation sweep and its unit test.

"""
    parse_quast_metric(report_tsv, metric) -> Union{Missing, Float64}

Return the numeric value of `metric` from a QUAST `report.tsv`, or `missing` if
the file/metric is absent or the value is non-numeric (QUAST prints "-" when a
metric could not be computed, e.g. NGA50 when nothing aligns — the empirical
signature we expect for a BROKEN assembly).
"""
function parse_quast_metric(report_tsv::String, metric::String)::Union{Missing, Float64}
    isfile(report_tsv) || return missing
    result::Union{Missing, Float64} = missing
    open(report_tsv, "r") do io
        for line in eachline(io)
            fields = split(line, '\t')
            if length(fields) >= 2 && strip(fields[1]) == metric
                parsed = tryparse(Float64, strip(fields[2]))
                result = parsed === nothing ? missing : parsed
                break
            end
        end
    end
    return result
end

"""
    empty_quast_metrics(metric_source="internal") -> NamedTuple

Single source of truth for the "no QUAST metrics" NamedTuple shape: the four
`quast_*` fields are all `missing` and `metric_source` records WHY QUAST is
absent (e.g. "internal", "internal:quast-disabled", "internal:quast-failed",
"internal:quast-skipped-empty", "internal:arm-failed"). Both this file's
absent-report branch and the sweep's inline fallback defaults call it, so the
fallback shape is defined exactly once.
"""
function empty_quast_metrics(metric_source::String = "internal")
    return (
        metric_source = metric_source,
        quast_genome_fraction = missing,
        quast_nga50 = missing,
        quast_num_misassemblies = missing,
        quast_duplication_ratio = missing
    )
end

"""
    quast_metrics_for_report(report_tsv) -> NamedTuple

Extract the alignment-validated assembly metrics for ONE assembly arm from its
QUAST `report.tsv` and label their provenance. Returns a NamedTuple with:

  * `metric_source`            — "quast" when the report exists (QUAST scored
                                 this arm), else "internal" (QUAST unavailable;
                                 the caller should fall back to its own
                                 size-ratio metrics).
  * `quast_genome_fraction`    — QUAST "Genome fraction (%)" (alignment-based)
  * `quast_nga50`              — QUAST "NGA50"
  * `quast_num_misassemblies`  — QUAST "# misassemblies"
  * `quast_duplication_ratio`  — QUAST "Duplication ratio"

Individual metrics may still be `missing` even when `metric_source == "quast"`
(QUAST prints "-" for a metric it could not compute, e.g. NGA50 when nothing
aligns); that is the honest state and the internal fallback column remains
populated regardless.
"""
function quast_metrics_for_report(report_tsv::String)
    if !isfile(report_tsv)
        return empty_quast_metrics("internal")
    end
    return (
        metric_source = "quast",
        quast_genome_fraction = parse_quast_metric(report_tsv, "Genome fraction (%)"),
        quast_nga50 = parse_quast_metric(report_tsv, "NGA50"),
        quast_num_misassemblies = parse_quast_metric(report_tsv, "# misassemblies"),
        quast_duplication_ratio = parse_quast_metric(report_tsv, "Duplication ratio")
    )
end
