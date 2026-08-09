# Corrector/Viterbi-arm counterpart to the constant-quality probe (bead td-4e19d.2)
#
# QUESTION
#   `qualmer_quality_channel_probe.jl` shows the DEFAULT greedy assembly path
#   (corrector = :none) emits byte-identical contigs under wildly different per-base
#   quality vectors. This script establishes the SCOPE of that result: is the
#   Stage-1 corrector / Viterbi decoder also quality-blind, or does it genuinely
#   consume Phred?
#
#   Answering that needs three DIFFERENT measurements, because they can disagree — and
#   on this codebase they do. Reporting only the last one would be wrong in both
#   directions, so all three are run and written out.
#
#   STAGE A — unit level. Call the emission model and the corrector's read-likelihood
#   directly with different quality vectors. Answers "do these functions read Phred at
#   all", with no gating or graph structure in the way.
#
#   STAGE B — decision level. Run the read decoder (`find_optimal_sequence_path`) over a
#   whole simulated read set under each quality condition and count how many reads come
#   back DIFFERENT. Answers "does Phred change what the decoder actually decides".
#
#   STAGE C — end-to-end. Run `assemble_genome(corrector = :iterative)` under each
#   condition and compare final assemblies, capturing the corrector diagnostics
#   (`skip_fraction_per_pass`, `decode_fraction_per_pass`, `cheap_corrections_total`)
#   alongside the verdict. Those diagnostics are REQUIRED, not decoration: an end-to-end
#   invariance means something completely different when `decode_fraction` is 0 (the
#   quality-consuming decoder never ran) than when it is nonzero (it ran and quality
#   changed nothing). Without them the invariance is uninterpretable.
#
# CONDITIONS
#   Identical to the main probe (`oracle` / `constant40` / `constant02` / `antioracle`),
#   built from the same simulator via `include`, so the two scripts compare like with
#   like. Everything is stratified by chemistry and never pooled.
#
# USAGE
#   julia --project=. benchmarking/qualmer_corrector_quality_sensitivity.jl
#   julia --project=. benchmarking/qualmer_corrector_quality_sensitivity.jl --smoke
#   julia --project=. benchmarking/qualmer_corrector_quality_sensitivity.jl --output-dir DIR

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import BioSequences
import FASTX
import SHA
import StableRNGs
import Dates

# Reuse the probe's simulator, quality-condition builders and digest. `include` rather
# than copy-paste: a divergence between the two read models would silently invalidate the
# comparison this script exists to make. PROBE_INCLUDE_ONLY suppresses the probe's own
# `main()` so including it only defines helpers.
PROBE_INCLUDE_ONLY = true
include(joinpath(@__DIR__, "qualmer_quality_channel_probe.jl"))

const CORRECTOR_CONDITIONS = ["oracle", "constant40", "constant02", "antioracle"]

const CORRECTOR_OUTPUT_DIR = something(arg_value("--output-dir"),
    joinpath(@__DIR__, "results", "qualmer_quality_channel_probe"))

# One coverage per chemistry — the corrector is expensive and the question here is
# directional (does quality reach a decision at all), not a dose-response curve.
const CORRECTOR_COVERAGE = 10
const CORRECTOR_SEED = 42

# === Stage A: does the quality machinery read Phred at all? ===

"""
Exercise the emission model and the corrector read-likelihood with high vs low quality
over IDENTICAL sequence. Any difference proves the function consumes Phred; identity
proves it does not. No graph traversal or gating is involved, so this isolates the
machinery from whether it is ever reached.
"""
function stage_a_unit_probe()
    rows = NamedTuple[]

    observed = "ACGTACGTAC"
    for (label, node) in (("match", "ACGTACGTAC"), ("mismatch", "ACGTACGTAA"))
        high = Mycelia.default_viterbi_emission_logp(
            observed, node, :dna; quality = fill(40.0, length(observed)))
        low = Mycelia.default_viterbi_emission_logp(
            observed, node, :dna; quality = fill(2.0, length(observed)))
        push!(rows,
            (stage = "A_unit", subject = "default_viterbi_emission_logp[$(label)]",
                high_quality_value = round(high; digits = 6),
                low_quality_value = round(low; digits = 6),
                consumes_quality = high != low))
    end

    # Read-likelihood over a real qualmer graph. The graph and the probe read are held
    # fixed; only the probe read's quality string changes.
    rng = StableRNGs.StableRNG(7)
    genome = BioSequences.LongDNA{4}(rand(rng, collect(DNA_BASES), 400))
    support = [FASTX.FASTQ.Record("support_$(i)", genome[1:200], fill(UInt8(30), 200))
               for i in 1:60]
    graph = Mycelia.Rhizomorph.build_qualmer_graph(support, K; mode = :canonical)

    high_likelihood = Mycelia.calculate_read_likelihood(
        FASTX.FASTQ.Record("probe", genome[1:200], fill(UInt8(40), 200)),
        graph, K; graph_mode = :canonical)
    low_likelihood = Mycelia.calculate_read_likelihood(
        FASTX.FASTQ.Record("probe", genome[1:200], fill(UInt8(2), 200)),
        graph, K; graph_mode = :canonical)
    push!(rows,
        (stage = "A_unit", subject = "calculate_read_likelihood",
            high_quality_value = round(high_likelihood; digits = 6),
            low_quality_value = round(low_likelihood; digits = 6),
            consumes_quality = high_likelihood != low_likelihood))

    return rows
end

# === Stage B: does Phred change what the decoder DECIDES? ===

"""
Decode every read under each quality condition and count how many corrected reads differ
from the `oracle` condition. This is the measurement that separates "the machinery reads
quality" (Stage A) from "quality changes an outcome".
"""
function stage_b_decision_probe(reads, profile)
    corrected = Dict{String, Vector{String}}()
    changed_counts = Dict{String, Int}()
    failure_counts = Dict{String, Int}()
    empty_counts = Dict{String, Int}()

    for condition in CORRECTOR_CONDITIONS
        records = fastq_records(reads, condition, profile)
        graph = Mycelia.Rhizomorph.build_qualmer_graph(records, K; mode = :canonical)
        outputs = String[]
        changed = 0
        failures = 0
        empty_results = 0
        for record in records
            # COUNT failures; do not let them masquerade as "the decoder chose not to
            # change this read". Passing a failed read through unchanged makes
            # `changed` stay flat, so a decoder that failed on EVERY read produces
            # `reads_differing_from_oracle = 0` in every condition — which is
            # bit-for-bit the signature of the finding this stage exists to establish
            # (that the decoder is quality-blind). Without a failure count in the
            # emitted row, a reader of the TSV cannot tell the two apart.
            improved = try
                decoded = Mycelia.find_optimal_sequence_path(
                    record, graph, K; graph_mode = :canonical)
                # `first` is OUTSIDE the decode call on purpose: a BoundsError here
                # means an EMPTY result, which is not a decoder failure and must not be
                # laundered into one.
                isempty(decoded) ? (empty_results += 1; record) : first(decoded)
            catch error
                @warn "decoder failed on read; passing through and counting" chemistry=profile.name condition exception=error
                failures += 1
                record
            end
            output = String(FASTX.sequence(improved))
            output == String(FASTX.sequence(record)) || (changed += 1)
            push!(outputs, output)
        end
        corrected[condition] = outputs
        changed_counts[condition] = changed
        failure_counts[condition] = failures
        empty_counts[condition] = empty_results
    end

    baseline = corrected["oracle"]
    rows = NamedTuple[]
    for condition in CORRECTOR_CONDITIONS
        differing = count(i -> corrected[condition][i] != baseline[i], eachindex(baseline))
        push!(rows,
            (stage = "B_decision", chemistry = profile.name,
                coverage = CORRECTOR_COVERAGE, seed = CORRECTOR_SEED,
                condition = condition, n_reads = length(baseline),
                reads_changed_by_decoder = changed_counts[condition],
                reads_differing_from_oracle = differing,
                # Without these two columns, "0 reads differed" is ambiguous between
                # "the decoder ran and quality changed nothing" and "the decoder threw
                # on every read". Only the first is a finding.
                n_decode_failures = failure_counts[condition],
                n_empty_decodes = empty_counts[condition],
                stage_b_trustworthy = failure_counts[condition] == 0))
    end
    return rows
end

# === Stage C: end-to-end, WITH the diagnostics needed to interpret it ===

function stage_c_end_to_end(reads, profile)
    rows = NamedTuple[]
    digests = String[]

    for condition in CORRECTOR_CONDITIONS
        records = fastq_records(reads, condition, profile)
        outcome = try
            elapsed = @elapsed result = Mycelia.Rhizomorph.assemble_genome(
                records; k = K, corrector = :iterative, verbose = false)
            (result = result, wall_seconds = elapsed)
        catch error
            @error "corrector assembly failed" chemistry=profile.name condition exception=(
                error, catch_backtrace())
            nothing
        end

        if outcome === nothing
            # Recorded, never silently skipped: a missing condition would otherwise read
            # as agreement in the verdict below.
            push!(digests, "ERROR")
            push!(rows,
                (stage = "C_end_to_end", chemistry = profile.name,
                    coverage = CORRECTOR_COVERAGE, seed = CORRECTOR_SEED,
                    condition = condition, n_reads = length(reads), digest = "ERROR",
                    n_contigs = -1, largest_contig = -1,
                    mean_joint_quality = "NA", skip_fraction_per_pass = "NA",
                    decode_fraction_per_pass = "NA", cheap_corrections_total = "NA",
                    wall_seconds = -1.0))
            continue
        end

        contigs = String.(outcome.result.contigs)
        lengths = length.(contigs)
        stats = outcome.result.assembly_stats
        digest = assembly_digest(contigs)
        push!(digests, digest)
        push!(rows,
            (stage = "C_end_to_end", chemistry = profile.name,
                coverage = CORRECTOR_COVERAGE, seed = CORRECTOR_SEED,
                condition = condition, n_reads = length(reads), digest = digest,
                n_contigs = length(contigs),
                largest_contig = isempty(lengths) ? 0 : maximum(lengths),
                # mean_joint_quality is the CONTROL on the whole stage: it must DIFFER
                # between conditions. If it did not, the conditions never reached the
                # graph and an "identical assemblies" verdict would be meaningless.
                mean_joint_quality = string(get(stats, "mean_joint_quality", "NA")),
                skip_fraction_per_pass = string(get(stats, "skip_fraction_per_pass", "NA")),
                decode_fraction_per_pass = string(
                    get(stats, "decode_fraction_per_pass", "NA")),
                cheap_corrections_total = string(
                    get(stats, "cheap_corrections_total", "NA")),
                wall_seconds = round(outcome.wall_seconds; digits = 3)))
    end

    return rows, digests
end

# === Verdicts ===

# Sentinel string written by every column of a failed Stage-C condition.
const CORRECTOR_ERROR = "ERROR"
const CORRECTOR_MISSING = "NA"

"""
Build one chemistry's verdict row from the Stage-B and Stage-C rows.

EXTRACTED FROM `corrector_main` ON PURPOSE. Everything upstream of this function needs
a reference genome, a read simulator and four `assemble_genome` runs, so the failure
semantics below were unreachable from any test — the only way to observe what a broken
run reports was to break a real run. `test/4_assembly/qualmer_corrector_verdict_test.jl`
now feeds it all-error rows directly.

THE RULE THIS ENCODES: a POSITIVE finding may rest on a single surviving condition, but
a NEGATIVE finding is a claim about every condition and therefore requires every
condition. When the evidence for a negative is missing, the field is `"NA"` — not
`false`.

`false` is not the safe default it looks like. `end_to_end_assemblies_identical = false`
asserts that the assemblies DIFFER, which on an all-failed run is a fabricated result in
the opposite direction; and `decoder_ever_ran = false` is the exact string the
INTERPRETATION GUIDE keys on to reclassify the whole study as a gating result. Both are
verdicts, so a failure that produces either is still a failure moving a verdict — the
defect PR #453 review found here, one layer down from the sentinel arithmetic it named.

`verdict_is_interpretable` is the single field to read first: `false` means at least one
stage lost data and no boolean below it should be quoted.
"""
function corrector_verdict(chemistry, n_reads, b_rows, c_rows, digests)
    ok_rows = filter(r -> r.digest != CORRECTOR_ERROR, c_rows)
    ok_digests = filter(!=(CORRECTOR_ERROR), digests)
    n_failed = length(digests) - length(ok_digests)
    stage_c_complete = n_failed == 0 && !isempty(ok_rows)

    # Stage B carries its own failure count per condition (`n_decode_failures`), and it
    # was previously dropped on the floor here: a decoder that threw on every read
    # returns `reads_differing_from_oracle = 0` in every condition, which is bit-for-bit
    # the signature of the quality-blindness finding this script exists to establish.
    stage_b_failures = sum(r.n_decode_failures for r in b_rows; init = 0)
    stage_b_complete = stage_b_failures == 0 && !isempty(b_rows)

    decoder_sensitive = any(
        r.condition != "oracle" && r.reads_differing_from_oracle > 0 for r in b_rows)
    decoder_decisions_depend_on_quality = decoder_sensitive ? true :
                                          stage_b_complete ? false : CORRECTOR_MISSING

    # Two DISTINCT conditions must survive for a difference to be observable at all, so
    # "no difference" from a single surviving condition is not a negative result.
    distinct_quality = length(unique(r.mean_joint_quality for r in ok_rows))
    quality_reached_graph = distinct_quality > 1 ? true :
                            stage_c_complete ? false : CORRECTOR_MISSING

    # Parse the per-pass numbers instead of pattern-matching a stringified vector. The
    # previous regex was `^`-anchored, so it inspected only pass 1 — which is 0.0 by
    # construction in every row — making the flag STRUCTURALLY incapable of returning
    # true, and it failed OPEN (true) on a non-vector value. It reported "the decoder
    # never ran" for illumina, which actually decodes on passes 2-4. This flag gates the
    # study's interpretation, so it gets parsed, not matched.
    decode_fractions(text) = text == CORRECTOR_MISSING ? Float64[] :
                             filter(!isnothing,
        [tryparse(Float64, strip(token))
         for token in split(strip(text, ['[', ']']), ',')
         if !isempty(strip(token))])
    saw_decoder_run = any(
        any(>(0.0), decode_fractions(r.decode_fraction_per_pass)) for r in ok_rows)
    decoder_ever_ran = saw_decoder_run ? true :
                       stage_c_complete ? false : CORRECTOR_MISSING
    max_decode_fraction = maximum(
        (maximum(decode_fractions(r.decode_fraction_per_pass); init = 0.0)
        for r in ok_rows); init = 0.0)

    # Identity across conditions is a statement about ALL of them, so ANY failure makes
    # it unanswerable. It cannot be `false` either: with every condition failing,
    # `unique(["ERROR", ...])` is one element and the run reported perfect agreement.
    assemblies_identical = stage_c_complete ? length(unique(ok_digests)) == 1 :
                           CORRECTOR_MISSING

    return (chemistry = chemistry, coverage = CORRECTOR_COVERAGE,
        seed = CORRECTOR_SEED, n_reads = n_reads,
        verdict_is_interpretable = stage_b_complete && stage_c_complete,
        n_conditions_failed = n_failed,
        n_stage_b_decode_failures = stage_b_failures,
        decoder_decisions_depend_on_quality = decoder_decisions_depend_on_quality,
        max_reads_differing_from_oracle = maximum(
            (r.reads_differing_from_oracle for r in b_rows); init = 0),
        end_to_end_assemblies_identical = assemblies_identical,
        end_to_end_quality_reached_graph = quality_reached_graph,
        end_to_end_decoder_ever_ran = decoder_ever_ran,
        end_to_end_max_decode_fraction = round(max_decode_fraction; digits = 4),
        end_to_end_skip_fraction = isempty(ok_rows) ? CORRECTOR_MISSING :
                                   first(ok_rows).skip_fraction_per_pass)
end

# === Main ===

function corrector_main()
    mkpath(CORRECTOR_OUTPUT_DIR)
    reference = load_reference(REFERENCE_FASTA)
    profiles = SMOKE ? CHEMISTRIES[1:1] : CHEMISTRIES

    unit_rows = stage_a_unit_probe()
    write_tsv(joinpath(CORRECTOR_OUTPUT_DIR, "qualmer_corrector_stageA_unit.tsv"), unit_rows)

    decision_rows = NamedTuple[]
    end_to_end_rows = NamedTuple[]
    verdicts = NamedTuple[]

    for profile in profiles
        reads = simulate_read_set(reference, profile, CORRECTOR_COVERAGE, CORRECTOR_SEED)
        @info "corrector cell" chemistry=profile.name n_reads=length(reads)

        b_rows = stage_b_decision_probe(reads, profile)
        append!(decision_rows, b_rows)

        c_rows, digests = stage_c_end_to_end(reads, profile)
        append!(end_to_end_rows, c_rows)

        push!(verdicts,
            corrector_verdict(profile.name, length(reads), b_rows, c_rows, digests))
    end

    write_tsv(joinpath(CORRECTOR_OUTPUT_DIR, "qualmer_corrector_stageB_decisions.tsv"),
        decision_rows)
    write_tsv(joinpath(CORRECTOR_OUTPUT_DIR, "qualmer_corrector_stageC_end_to_end.tsv"),
        end_to_end_rows)
    write_tsv(joinpath(CORRECTOR_OUTPUT_DIR, "qualmer_corrector_quality_verdicts.tsv"),
        verdicts)

    report_corrector(unit_rows, decision_rows, verdicts)
    return verdicts
end

function report_corrector(unit_rows, decision_rows, verdicts)
    println("\n=== STAGE A — does the quality machinery read Phred at all? ===")
    for row in unit_rows
        println("  $(rpad(row.subject, 46)) highQ=$(rpad(row.high_quality_value, 14)) " *
                "lowQ=$(rpad(row.low_quality_value, 14)) consumes_quality=$(row.consumes_quality)")
    end

    println("\n=== STAGE B — does Phred change what the decoder DECIDES? ===")
    for chem in unique(r.chemistry for r in decision_rows)
        println("  [$(chem)]")
        for row in filter(r -> r.chemistry == chem, decision_rows)
            println("    $(rpad(row.condition, 12)) " *
                    "reads_changed_by_decoder=$(rpad(row.reads_changed_by_decoder, 5))" *
                    "reads_differing_from_oracle=$(row.reads_differing_from_oracle)" *
                    "/$(row.n_reads)")
        end
    end

    println("\n=== STAGE C — end-to-end assemble_genome(corrector = :iterative) ===")
    for v in verdicts
        println("  [$(v.chemistry)] interpretable=$(v.verdict_is_interpretable) " *
                "failed_conditions=$(v.n_conditions_failed) " *
                "stageB_decode_failures=$(v.n_stage_b_decode_failures)")
        println("      assemblies_identical=$(v.end_to_end_assemblies_identical) " *
                "quality_reached_graph=$(v.end_to_end_quality_reached_graph) " *
                "decoder_ever_ran=$(v.end_to_end_decoder_ever_ran) " *
                "skip_fraction=$(v.end_to_end_skip_fraction)")
    end

    println("\nINTERPRETATION GUIDE")
    println("  Read verdict_is_interpretable FIRST. When it is false, a stage lost " *
            "data and the")
    println("  booleans on that row are reported as NA rather than guessed — quote " *
            "none of them.")
    println("  Stage C identity is only evidence of quality-blindness when " *
            "decoder_ever_ran=true.")
    println("  If decoder_ever_ran=false the reads were gated out before the " *
            "quality-consuming decoder,")
    println("  so the end-to-end invariance is a GATING result and says nothing about " *
            "quality-blindness.")
    return nothing
end

# Guard the driver so a test can `include()` this file for its pure verdict logic
# without a reference download, a read simulation, or four assemblies per chemistry.
# Same convention as track_a_baseline_benchmark.jl, and what makes
# test/4_assembly/qualmer_corrector_verdict_test.jl possible.
if abspath(PROGRAM_FILE) == @__FILE__
    corrector_main()
end
