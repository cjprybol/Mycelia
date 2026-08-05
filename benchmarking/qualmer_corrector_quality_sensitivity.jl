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

    for condition in CORRECTOR_CONDITIONS
        records = fastq_records(reads, condition, profile)
        graph = Mycelia.Rhizomorph.build_qualmer_graph(records, K; mode = :canonical)
        outputs = String[]
        changed = 0
        for record in records
            improved = try
                first(Mycelia.find_optimal_sequence_path(
                    record, graph, K; graph_mode = :canonical))
            catch error
                # A read the decoder cannot handle is passed through unchanged rather
                # than dropped: dropping it would shorten the vector and misalign the
                # per-read comparison against the oracle condition.
                @warn "decoder failed on read; passing through" chemistry=profile.name condition exception=error
                record
            end
            output = String(FASTX.sequence(improved))
            output == String(FASTX.sequence(record)) || (changed += 1)
            push!(outputs, output)
        end
        corrected[condition] = outputs
        changed_counts[condition] = changed
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
                reads_differing_from_oracle = differing))
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

        decoder_sensitive = any(
            r.condition != "oracle" && r.reads_differing_from_oracle > 0 for r in b_rows)
        quality_reached_graph = length(unique(r.mean_joint_quality for r in c_rows)) > 1
        decoder_ever_ran = any(
            !occursin(r"^\[?0\.0[, \]]", r.decode_fraction_per_pass) &&
            r.decode_fraction_per_pass != "NA"
        for r in c_rows)

        push!(verdicts,
            (chemistry = profile.name, coverage = CORRECTOR_COVERAGE,
                seed = CORRECTOR_SEED, n_reads = length(reads),
                decoder_decisions_depend_on_quality = decoder_sensitive,
                max_reads_differing_from_oracle = maximum(
                    r.reads_differing_from_oracle for r in b_rows),
                end_to_end_assemblies_identical = length(unique(digests)) == 1,
                end_to_end_quality_reached_graph = quality_reached_graph,
                end_to_end_decoder_ever_ran = decoder_ever_ran,
                end_to_end_skip_fraction = first(
                    r.skip_fraction_per_pass for r in c_rows)))
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
        println("  [$(v.chemistry)] assemblies_identical=$(v.end_to_end_assemblies_identical) " *
                "quality_reached_graph=$(v.end_to_end_quality_reached_graph) " *
                "decoder_ever_ran=$(v.end_to_end_decoder_ever_ran) " *
                "skip_fraction=$(v.end_to_end_skip_fraction)")
    end

    println("\nINTERPRETATION GUIDE")
    println("  Stage C identity is only evidence of quality-blindness when " *
            "decoder_ever_ran=true.")
    println("  If decoder_ever_ran=false the reads were gated out before the " *
            "quality-consuming decoder,")
    println("  so the end-to-end invariance is a GATING result and says nothing about " *
            "quality-blindness.")
    return nothing
end

corrector_main()
