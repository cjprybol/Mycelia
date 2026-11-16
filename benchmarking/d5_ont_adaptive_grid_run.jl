# D5 ONT adaptive-grid RUNNER (td-yddpr.5). CLI over d5_ont_adaptive_grid_driver.jl.
#
# This is the committed invocation that regenerates the per-cell table:
#
#   julia --project=. benchmarking/d5_ont_adaptive_grid_run.jl --list
#   julia --project=. benchmarking/d5_ont_adaptive_grid_run.jl --phase A \
#       --read-set-index <0..95> --output-dir <root>
#   julia --project=. benchmarking/d5_ont_adaptive_grid_run.jl --phase merge \
#       --output-dir <root>
#   # Only after a human has set tau and tau_hat (lock step 2):
#   julia --project=. benchmarking/d5_ont_adaptive_grid_run.jl --phase B \
#       --tau <x> --tau-hat <y> [--read-set-index i] --output-dir <root>
#
# `--toy` restricts every phase to the pre-declared toy subset (3 read sets,
# 7 cells); point it at its own --output-dir. On Perlmutter use
# run_d5_grid_nersc.sbatch, one array task per read set.
#
# Output tree:
#   <root>/run_manifest.json                  Mycelia SHA shared by every phase
#   <root>/refs/                              reference FASTA per organism
#   <root>/readsets/<read_set_id>/readset.json reads, measured identity, provenance
#   <root>/readsets/<read_set_id>/survival.json 23-rung S / S_hat on RAW reads
#   <root>/cells/<cell_id>/cell_result.json   one checkpoint per cell (atomic)
#   <root>/cells/<cell_id>/corrected_survival.json  diagnostic, :iterative only
#   <root>/d5_cells.tsv, d5_survival.tsv, illumina_none_optimum.tsv,
#   <root>/replicate_evaluability.tsv         written by --phase merge
#
# Censored metrics are written as NA (TSV) / null (JSON), never 0.0.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import BioSequences
import CSV
import DataFrames
import Dates
import FASTX
import JSON
import Mycelia
import Statistics

include(joinpath(@__DIR__, "d5_ont_adaptive_grid_driver.jl"))
import .D5OntAdaptiveGridDriver

const QUAST_THREADS = 2

# === Provenance ===

function git_provenance()
    repo = dirname(@__DIR__)
    sha = try
        strip(read(`git -C $(repo) rev-parse HEAD`, String))
    catch e
        @warn "could not read git SHA" exception = e
        "unknown"
    end
    dirty = try
        !isempty(strip(read(`git -C $(repo) status --porcelain --untracked-files=no`,
            String)))
    catch
        missing
    end
    return (; sha, dirty)
end

function tool_version(env::AbstractString, cmd::Cmd)
    try
        out = read(`$(Mycelia.CONDA_RUNNER) run -n $(env) $(cmd)`, String)
        return strip(first(split(strip(out), '\n')))
    catch e
        @warn "could not read tool version" env exception = e
        return "unknown"
    end
end

function tool_versions()
    return Dict(
        "badread" => tool_version("badread", `badread --version`),
        "minimap2" => tool_version("minimap2", `minimap2 --version`),
        "quast" => tool_version("quast", `quast.py --version`),
        "julia" => string(VERSION)
    )
end

function run_context_base()
    git = git_provenance()
    return (; mycelia_git_sha = git.sha, git_dirty = git.dirty,
        host = gethostname(),
        slurm_job_id = something(get(ENV, "SLURM_ARRAY_JOB_ID", nothing),
            get(ENV, "SLURM_JOB_ID", nothing), missing))
end

"""The registered simulation pins, as stored in (and compared against) the manifest."""
function manifest_pins()
    return Dict{String, Any}(
        "badread_identity_params" => Dict(string(k) => v
        for (k, v) in D5OntAdaptiveGridDriver.D5_BADREAD_IDENTITY_PARAMS),
        "badread_error_model" => D5OntAdaptiveGridDriver.D5_BADREAD_ERROR_MODEL,
        "badread_qscore_model" => D5OntAdaptiveGridDriver.D5_BADREAD_QSCORE_MODEL,
        "badread_length_dist" => D5OntAdaptiveGridDriver.D5_BADREAD_LENGTH_DIST,
        "illumina_params" => Dict(String(k) => v
        for (k, v) in pairs(D5OntAdaptiveGridDriver.D5_ILLUMINA_PARAMS)),
        "strategy" => String(D5OntAdaptiveGridDriver.D5_STRATEGY))
end

"""
One grid, one code version. The first writer records the Mycelia SHA and the
simulation pins; every later task (Phase A shards, merge, Phase B) must match.
An unreadable SHA, a dirty tree, or an unreadable `git status` cannot prove that,
so all refuse too. This keeps a phased run a single grid execution
(PLAN-2026-03-28 step 2). `allow_mismatch` (--allow-sha-mismatch) downgrades the
CODE-VERSION refusals to warnings. It does NOT cover changed simulation pins or a
toy/full mix: those change what the grid measures, and no flag waives them.
"""
function check_run_manifest!(root::AbstractString, git; allow_mismatch::Bool, toy::Bool)
    refuse(msg) = allow_mismatch ? (@warn msg) : error(msg)
    git.sha == "unknown" &&
        refuse("cannot read the Mycelia git SHA, so this task cannot be tied to the " *
               "grid's code version. Fix git access or pass --allow-sha-mismatch.")
    git.dirty === false ||
        refuse("the Mycelia checkout has uncommitted changes to tracked files (or git " *
               "status could not be read; dirty=$(git.dirty)); the SHA would not " *
               "describe the code that ran. Commit or pass --allow-sha-mismatch.")
    path = joinpath(root, "run_manifest.json")
    if !isfile(path)
        D5OntAdaptiveGridDriver.save_json_atomically(path,
            merge(manifest_pins(),
                Dict{String, Any}("mycelia_git_sha" => git.sha, "toy" => toy,
                    "created" => string(Dates.now()))))
    end
    manifest = JSON.parsefile(path)
    manifest["toy"] == toy ||
        error("output dir $(root) belongs to a $(manifest["toy"] ? "toy" : "full") run; " *
              "refusing to mix toy and full cells in one tree")
    recorded = manifest["mycelia_git_sha"]
    recorded == git.sha ||
        refuse("run_manifest.json records Mycelia $(recorded) but this checkout is " *
               "$(git.sha). A grid is one code version; pass --allow-sha-mismatch only " *
               "if the difference cannot touch assembly output, and say so in the results.")
    for (key, value) in manifest_pins()
        recorded_pin = get(manifest, key, nothing)
        recorded_pin == value ||
            error("run_manifest.json pins $(key) = $(repr(recorded_pin)) but this code " *
                  "has $(repr(value)); the simulation parameters changed under the grid. " *
                  "Use a fresh output root.")
    end
    return manifest
end

"""
EXACT reference length: a tolerance would accept a partially written download,
which would then feed simulation, survival and QUAST at once. `expected_bp` is
the RefSeq record length.
"""
function check_reference_length(genome_size::Integer, info, path::AbstractString)
    genome_size == info.expected_bp ||
        error("reference $(info.accession) at $(path) is $(genome_size) bp; the RefSeq " *
              "record is $(info.expected_bp) bp. Delete the file and re-stage it.")
    return nothing
end

"""A reused read set must still hold the read count it was recorded with."""
function check_read_count(records, meta, rs_dir::AbstractString)
    length(records) == meta["n_reads"] ||
        error("read set $(meta["read_set_id"]) has $(length(records)) reads on disk but " *
              "readset.json recorded $(meta["n_reads"]); delete $(rs_dir) to regenerate it")
    return records
end

# === Reference ===

function prepare_reference(organism::Symbol, refs_dir::AbstractString)
    info = D5OntAdaptiveGridDriver.D5_ORGANISM_REFERENCES[organism]
    mkpath(refs_dir)
    path = Mycelia.download_genome_by_accession(accession = info.accession,
        outdir = refs_dir, compressed = false)
    (isfile(path) && filesize(path) > 0) ||
        error("reference download failed for $(organism) ($(info.accession)): $(path)")
    seqs = read_sequences(path)
    genome_size = sum(length, seqs)
    check_reference_length(genome_size, info, path)
    return (; path, seqs, genome_size, accession = info.accession)
end

function read_records(path::AbstractString)
    reader = Mycelia.open_fastx(path)
    try
        return collect(reader)
    finally
        close(reader)
    end
end

function read_sequences(path)
    [FASTX.sequence(BioSequences.LongDNA{4}, r) for r in read_records(path)]
end

# === Read sets ===

function survival_json(rows)
    return [Dict(
                "k" => r.k, "S" => r.S, "n_ref_positions" => r.n_ref_positions,
                "n_ref_positions_observed" => r.n_ref_positions_observed,
                "n_distinct_read_kmers" => r.n_distinct_read_kmers,
                "S_hat" => r.S_hat, "s_hat_v" => r.s_hat_v, "s_hat_m" => r.s_hat_m,
                "s_hat_status" => r.s_hat_status,
                "histogram" => Dict(string(m) => n for (m, n) in r.histogram))
            for r in rows]
end

"""Realized ONT identity by minimap2 alignment (ont_read_identity.jl's scorer)."""
function measure_identity(reads_path::AbstractString, ref_path::AbstractString,
        work_dir::AbstractString)
    Mycelia.add_bioconda_env("minimap2")
    sam = joinpath(work_dir, "reads.sam")
    run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2
                  -ax map-ont --secondary=no $(ref_path) $(reads_path)`,
        stdout = sam, stderr = devnull))
    alignment = D5OntAdaptiveGridDriver.ReadIdentity.identities_from_sam(sam)
    rm(sam; force = true)
    blast = [r.blast_identity for r in alignment.rows if r.mapped]
    return Dict(
        "identity_measured_mean" => isempty(blast) ? nothing : Statistics.mean(blast),
        "identity_measured_median" => isempty(blast) ? nothing : Statistics.median(blast),
        "identity_n_primary_scored" => alignment.n_primary,
        "identity_n_unmapped" => alignment.n_unmapped,
        "identity_n_skipped" => alignment.n_skipped,
        "identity_method" => "minimap2 -ax map-ont --secondary=no; BLAST identity")
end

"""
Simulate (or reuse) one read set and compute its survival vector once.
`readset.json` is written LAST, so its presence means the read set is complete;
anything else in the directory is discarded and regenerated.
"""
function prepare_read_set(rs, ref, root::AbstractString)
    rs_id = D5OntAdaptiveGridDriver.read_set_id(rs)
    rs_dir = joinpath(root, "readsets", rs_id)
    marker = joinpath(rs_dir, "readset.json")
    if isfile(marker)
        meta = JSON.parsefile(marker)
        records = check_read_count(reduce(vcat, read_records.(meta["reads_paths"])),
            meta, rs_dir)
        return (; meta, records)
    end
    isdir(rs_dir) && rm(rs_dir; recursive = true)
    work = joinpath(rs_dir, "work")
    mkpath(work)

    meta = Dict{String, Any}("read_set_id" => rs_id,
        "technology" => String(rs.technology), "identity_target" => rs.identity,
        "coverage" => rs.coverage, "organism" => String(rs.organism),
        "seed" => rs.seed, "accession" => ref.accession,
        "genome_size" => ref.genome_size)
    if rs.technology === :ont
        # One kwargs source for the simulator AND the recorded argv.
        badread = D5OntAdaptiveGridDriver.badread_kwargs(rs, ref.path)
        argv = Mycelia._badread_nanopore_args(; badread...)
        tmp = joinpath(work, "reads.fq.gz")
        Mycelia.simulate_nanopore_reads(; badread..., outfile = tmp, quiet = true)
        final = joinpath(rs_dir, "reads.fq.gz")
        mv(tmp, final; force = true)
        meta["reads_paths"] = [final]
        meta["identity_params"] = D5OntAdaptiveGridDriver.D5_BADREAD_IDENTITY_PARAMS[rs.identity]
        meta["simulator_argv"] = join(argv, " ")
        meta["error_model"] = D5OntAdaptiveGridDriver.D5_BADREAD_ERROR_MODEL
        meta["qscore_model"] = D5OntAdaptiveGridDriver.D5_BADREAD_QSCORE_MODEL
        merge!(meta, measure_identity(final, ref.path, work))
    else
        params = D5OntAdaptiveGridDriver.D5_ILLUMINA_PARAMS
        sim = Mycelia.simulate_illumina_reads(; fasta = ref.path,
            coverage = rs.coverage, outbase = joinpath(work, "illumina"),
            rndSeed = rs.seed, quiet = true, params...)
        finals = String[]
        for (i, p) in enumerate((sim.forward_reads, sim.reverse_reads))
            p === nothing && continue
            dest = joinpath(rs_dir, "reads_R$(i).fq.gz")
            mv(p, dest; force = true)
            push!(finals, dest)
        end
        meta["reads_paths"] = finals
        meta["identity_params"] = nothing
        meta["simulator_argv"] = "art_illumina coverage=$(rs.coverage) " *
                                 "rndSeed=$(rs.seed) $(params)"
        meta["error_model"] = "art:$(params.seqSys)"
        meta["qscore_model"] = "art:$(params.seqSys)"
        merge!(meta,
            Dict("identity_measured_mean" => nothing,
                "identity_measured_median" => nothing,
                "identity_method" => "not measured (Illumina control)"))
    end
    meta["tool_versions"] = tool_versions()
    records = reduce(vcat, read_records.(meta["reads_paths"]))
    meta["n_reads"] = length(records)
    meta["read_bases"] = sum(r -> length(FASTX.sequence(r)), records; init = 0)
    # DIAGNOSTIC ONLY (DEC-2026-08-09 dispositions): recorded, never used for k.
    meta["error_rate_estimate_diagnostic"] = Mycelia.Rhizomorph.estimate_residual_error(records)

    seqs = [FASTX.sequence(BioSequences.LongDNA{4}, r) for r in records]
    t0 = time()
    rows = D5OntAdaptiveGridDriver.survival_ladder(seqs, ref.seqs)
    meta["survival_seconds"] = round(time() - t0; digits = 1)
    meta["survival_reads"] = "raw"
    D5OntAdaptiveGridDriver.save_json_atomically(joinpath(rs_dir, "survival.json"),
        Dict("read_set_id" => rs_id, "reads" => "raw", "strand" => "canonical",
            "rungs" => survival_json(rows)))
    rm(work; recursive = true, force = true)
    D5OntAdaptiveGridDriver.save_json_atomically(marker, meta)
    return (; meta, records)
end

function cell_context(meta, base)
    return (; identity_params = something(meta["identity_params"], missing),
        simulator_argv = meta["simulator_argv"], error_model = meta["error_model"],
        qscore_model = meta["qscore_model"], n_reads = meta["n_reads"],
        read_bases = meta["read_bases"],
        error_rate_estimate_diagnostic = meta["error_rate_estimate_diagnostic"],
        base.mycelia_git_sha, base.git_dirty, base.host, base.slurm_job_id)
end

# === Cells ===

"""Assemble at `k`, score with QUAST, and classify censoring. Returns a row."""
function assemble_cell(cell, k::Int, records, ref, cell_dir, context; phase::Symbol,
        row_kwargs...)
    D5OntAdaptiveGridDriver.assembler_accepts_k(k) ||
        return D5OntAdaptiveGridDriver.infeasible_row(cell, context; phase,
            k_selected = k, row_kwargs...)
    mkpath(cell_dir)
    rs = D5OntAdaptiveGridDriver.read_set_of(cell)
    iterative = cell.corrector === :iterative
    asm_dir = joinpath(cell_dir, "assembler")
    kwargs = (; k, corrector = cell.corrector,
        strategy = D5OntAdaptiveGridDriver.D5_STRATEGY,
        sequencing_tech = D5OntAdaptiveGridDriver.sequencing_tech_for(rs), verbose = false)
    if iterative
        mkpath(asm_dir)
        kwargs = merge(kwargs, (; reassembly_k_policy = :fixed, output_dir = asm_dir))
    end
    GC.gc()
    measured = D5OntAdaptiveGridDriver.timed_with_peak_rss(
        () -> Mycelia.Rhizomorph.assemble_genome(records; kwargs...))
    result = measured.value
    reassembly_k = missing
    if iterative
        reassembly_k = result.assembly_stats["reassembly_k"]
        reassembly_k == k ||
            error("reassembly_k=$(reassembly_k) != requested k=$(k) under " *
                  "reassembly_k_policy=:fixed; the comparator rung is not fixed")
    end

    contigs_path = joinpath(cell_dir, "contigs.fasta")
    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)
        end
    end
    min_contig = D5OntAdaptiveGridDriver.quast_min_contig(ref.genome_size)
    stats = D5OntAdaptiveGridDriver.contig_stats(result.contigs, min_contig)
    quast_dir = joinpath(cell_dir, "quast")
    report = joinpath(quast_dir, "report.tsv")
    quast_error = missing
    if stats.n_ge_min > 0 && filesize(contigs_path) > 0
        try
            Mycelia.run_quast([contigs_path]; outdir = quast_dir, reference = ref.path,
                min_contig, threads = QUAST_THREADS)
        catch e
            @warn "QUAST failed" cell_dir exception = (e, catch_backtrace())
            quast_error = first(sprint(showerror, e), 2000)
        end
    end
    quast_ran = isfile(report)
    metrics = quast_ran ? D5OntAdaptiveGridDriver.parse_quast_metrics(report) :
              D5OntAdaptiveGridDriver.empty_metrics()
    nga50_status = D5OntAdaptiveGridDriver.nga50_status_for(metrics, stats.n,
        stats.n_ge_min, quast_ran)
    status = nga50_status == "quast_failed" ? D5OntAdaptiveGridDriver.STATUS_QUAST_FAILED :
             D5OntAdaptiveGridDriver.STATUS_OK

    # The diagnostic must never cost the cell its assembly: a failure here is
    # recorded in its own column, and the cell keeps its QUAST result.
    corrected_status = if iterative
        try
            corrected_survival_diagnostic(asm_dir, ref, cell_dir)
        catch e
            @warn "corrected-read survival diagnostic failed" cell_dir exception = (e,
                catch_backtrace())
            D5OntAdaptiveGridDriver.save_json_atomically(
                joinpath(cell_dir, "corrected_survival.json"),
                Dict("status" => "failed", "role" => "diagnostic",
                    "error" => first(sprint(showerror, e), 2000)))
            "failed"
        end
    else
        missing
    end

    return D5OntAdaptiveGridDriver.cell_row(cell, context; phase, status,
        k_selected = k, reassembly_k,
        reassembly_k_policy = iterative ? "fixed" : missing,
        n_contigs = stats.n, contigs_ge_min = stats.n_ge_min,
        quast_min_contig = min_contig, metrics, nga50_status,
        corrected_survival_status = corrected_status, error_message = quast_error,
        wall_seconds = round(measured.wall_seconds; digits = 3),
        peak_rss_bytes = measured.peak_rss_bytes, peak_rss_method = measured.method,
        row_kwargs...)
end

"""
DIAGNOSTIC: survival on this cell's own corrected reads (user decision
2026-11-16). Never feeds k*. Returns the status recorded in the cell row
(`"measured"` or `"unavailable"`); an absent corrected FASTQ breaks the
`output_dir` contract, so it is warned for EVERY cell, not once.
"""
function corrected_survival_diagnostic(asm_dir, ref, cell_dir)
    corrected = joinpath(asm_dir, "corrected.fastq")
    out = joinpath(cell_dir, "corrected_survival.json")
    if !isfile(corrected)
        @warn "corrected FASTQ not found; corrected-read survival unavailable" corrected
        D5OntAdaptiveGridDriver.save_json_atomically(out,
            Dict("status" => "unavailable", "role" => "diagnostic"))
        return "unavailable"
    end
    seqs = read_sequences(corrected)
    rows = D5OntAdaptiveGridDriver.survival_ladder(seqs, ref.seqs)
    D5OntAdaptiveGridDriver.save_json_atomically(out,
        Dict("status" => "measured", "role" => "diagnostic", "reads" => "corrected",
            "strand" => "canonical", "rungs" => survival_json(rows)))
    rm(corrected; force = true)
    return "measured"
end

function checkpoint_path(root, cell)
    return joinpath(root, "cells", D5OntAdaptiveGridDriver.cell_id(cell), "cell_result.json")
end

"""
Skip a completed cell; retry error / quast_failed rows; record failures as rows.

`stale(prior)` returns `nothing` when a cached row is valid for this run, or a
message saying why it is not (Phase B: computed under another tau). A stale row
is an error unless `replan`, which recomputes it. `error_kwargs` are carried onto
the error row a failure produces (Phase B: tau, k*, floor), so the rows that most
need debugging keep their plan. Returns the row's status.
"""
function run_or_resume(f, root, cell, context; phase::Symbol,
        stale = prior -> nothing, replan::Bool = false, error_kwargs = (;))
    ckpt = checkpoint_path(root, cell)
    if isfile(ckpt)
        prior = try
            D5OntAdaptiveGridDriver.load_row(ckpt)
        catch e
            @warn "unreadable checkpoint; recomputing" ckpt exception = e
            nothing
        end
        if prior !== nothing &&
           !(prior["status"] in D5OntAdaptiveGridDriver.RETRYABLE_STATUSES)
            why = stale(prior)
            if why === nothing
                println("  cached  $(D5OntAdaptiveGridDriver.cell_id(cell)): $(prior["status"])")
                return prior["status"]
            end
            replan ||
                error("$(why) Pass --replan to recompute it, or use a fresh output root.")
            println("  replan  $(D5OntAdaptiveGridDriver.cell_id(cell)): $(why)")
        end
    end
    cell_dir = dirname(ckpt)
    isdir(cell_dir) && rm(cell_dir; recursive = true)
    row = try
        f(cell_dir)
    catch e
        @warn "cell failed" cell=D5OntAdaptiveGridDriver.cell_id(cell) exception=(e,
            catch_backtrace())
        D5OntAdaptiveGridDriver.cell_row(cell, context; phase,
            status = D5OntAdaptiveGridDriver.STATUS_ERROR,
            error_message = first(sprint(showerror, e), 2000), error_kwargs...)
    end
    D5OntAdaptiveGridDriver.save_json_atomically(ckpt, D5OntAdaptiveGridDriver.row_to_json(row))
    println("  $(row.status)  $(row.cell_id)  nga50_status=$(row.nga50_status) " *
            "NGA50=$(row.NGA50) GF=$(row.genome_fraction)")
    return row.status
end

"""Statuses that mean a resubmit is needed; a phase with any of these exits non-zero."""
function needs_retry(statuses)
    count(s -> s in D5OntAdaptiveGridDriver.RETRYABLE_STATUSES, statuses)
end

git_state(base) = (; sha = base.mycelia_git_sha, dirty = base.git_dirty)

# === Phases ===

function selected_read_sets(toy::Bool)
    toy ? D5OntAdaptiveGridDriver.toy_read_sets() : D5OntAdaptiveGridDriver.read_sets()
end

function phase_cells(rs, phase::Symbol, toy::Bool)
    phase === :A && toy && return D5OntAdaptiveGridDriver.toy_cells(rs)
    cells = D5OntAdaptiveGridDriver.cells_for_read_set(rs; phase)
    if toy
        # Toy Phase B covers the adaptive cells of the toy read sets only.
        cells = filter(
            c -> c.corrector in Set(t.corrector
            for t in D5OntAdaptiveGridDriver.toy_cells(rs)),
            cells)
    end
    return cells
end

function read_sets_for(cli)
    sets = selected_read_sets(cli.toy)
    cli.read_set_index === nothing && return sets
    (0 <= cli.read_set_index < length(sets)) ||
        error("--read-set-index $(cli.read_set_index) out of range 0:$(length(sets) - 1)")
    return [sets[cli.read_set_index + 1]]
end

function run_phase_a(cli, root)
    cli.read_set_index === nothing &&
        error("--phase A requires --read-set-index (one read set per array task)")
    base = run_context_base()
    check_run_manifest!(root, git_state(base);
        allow_mismatch = cli.allow_sha_mismatch, toy = cli.toy)
    statuses = String[]
    for rs in read_sets_for(cli)
        println("=== read set $(D5OntAdaptiveGridDriver.read_set_id(rs)) ===")
        ref = prepare_reference(rs.organism, joinpath(root, "refs"))
        prepared = prepare_read_set(rs, ref, root)
        context = cell_context(prepared.meta, base)
        for cell in phase_cells(rs, :A, cli.toy)
            push!(statuses,
                run_or_resume(root, cell, context; phase = :A) do cell_dir
                    assemble_cell(cell, Int(cell.k_policy), prepared.records, ref,
                        cell_dir, context; phase = :A)
                end)
        end
    end
    return statuses
end

"""
Phase B row metadata. tau and tau_hat are recorded on every row (a Phase B row
is a function of both), together with the collision floor and the S_hat reason.
"""
function phase_b_kwargs(plan, cli)
    hat = plan.k_star_hat === D5OntAdaptiveGridDriver.NO_FEASIBLE_K ?
          String(plan.k_star_hat) : plan.k_star_hat
    return (; k_star_reason = plan.reason, k_star_hat = hat,
        k_star_hat_reason = plan.reason_hat, tau = cli.tau, tau_hat = cli.tau_hat,
        k_floor = plan.floor_k)
end

function phase_b_stale(cli)
    return prior -> D5OntAdaptiveGridDriver.phase_b_checkpoint_matches(prior, cli.tau,
        cli.tau_hat) ? nothing :
                    "Phase B checkpoint $(prior["cell_id"]) was computed with " *
                    "tau=$(get(prior, "tau", missing)), tau_hat=$(get(prior, "tau_hat", missing)); " *
                    "this run has tau=$(cli.tau), tau_hat=$(cli.tau_hat)."
end

function run_phase_b(cli, root)
    base = run_context_base()
    check_run_manifest!(root, git_state(base);
        allow_mismatch = cli.allow_sha_mismatch, toy = cli.toy)
    statuses = String[]
    for rs in read_sets_for(cli)
        rs_id = D5OntAdaptiveGridDriver.read_set_id(rs)
        marker = joinpath(root, "readsets", rs_id, "readset.json")
        isfile(marker) || error("Phase A has not completed read set $(rs_id)")
        meta = JSON.parsefile(marker)
        survival = JSON.parsefile(joinpath(root, "readsets", rs_id, "survival.json"))
        context = cell_context(meta, base)
        ref = prepare_reference(rs.organism, joinpath(root, "refs"))
        records = nothing
        for cell in phase_cells(rs, :B, cli.toy)
            plan = D5OntAdaptiveGridDriver.phase_b_plan(cell, survival["rungs"], cli.tau,
                cli.tau_hat; genome_size = meta["genome_size"])
            row_kwargs = phase_b_kwargs(plan, cli)
            action = phase_b_action(plan, rs, cli.toy)
            k_sel = plan.k_star === D5OntAdaptiveGridDriver.NO_FEASIBLE_K ?
                    String(plan.k_star) : plan.k_star
            # Load and check the reads OUTSIDE the per-cell closure: a corrupt read
            # set must fail the task once (exit 1, "fix the read set"), not become
            # one retryable error row per cell that a resubmit cannot fix.
            if action === :assemble && records === nothing
                records = check_read_count(
                    reduce(vcat, read_records.(meta["reads_paths"])), meta,
                    dirname(marker))
            end
            push!(statuses,
                run_or_resume(root, cell, context; phase = :B,
                    stale = phase_b_stale(cli), replan = cli.replan,
                    error_kwargs = merge(row_kwargs, (; k_selected = k_sel))) do cell_dir
                    if action === :no_feasible_k
                        D5OntAdaptiveGridDriver.cell_row(cell, context; phase = :B,
                            status = D5OntAdaptiveGridDriver.STATUS_NO_FEASIBLE_K,
                            k_selected = String(D5OntAdaptiveGridDriver.NO_FEASIBLE_K),
                            row_kwargs...)
                    elseif action === :infeasible
                        D5OntAdaptiveGridDriver.infeasible_row(cell, context; phase = :B,
                            k_selected = plan.k_star, row_kwargs...)
                    elseif action === :reuse
                        reuse_row(root, cell, plan, context, row_kwargs)
                    else
                        assemble_cell(cell, plan.k_star, records, ref, cell_dir, context;
                            phase = :B, row_kwargs...)
                    end
                end)
        end
    end
    return statuses
end

"""
What Phase B actually does for a planned cell. The toy runs only some comparator
rungs, so a reuse whose source lies OUTSIDE the toy's Phase A cells is assembled
at the same k* instead (the identical cell, recomputed; not a fallback k). A
source INSIDE the toy, or any source in the full grid, keeps `:reuse`, so a
missing Phase A checkpoint there fails loudly rather than being papered over.
"""
function phase_b_action(plan, rs, toy::Bool)
    plan.action === :reuse && toy || return plan.action
    toy_ids = Set(D5OntAdaptiveGridDriver.cell_id.(D5OntAdaptiveGridDriver.toy_cells(rs)))
    return plan.reuse_cell_id in toy_ids ? :reuse : :assemble
end

"""Phase B row copied from the Phase A comparator assembly at k*, recorded as reuse."""
function reuse_row(root, cell, plan, context, row_kwargs)
    source = joinpath(root, "cells", plan.reuse_cell_id, "cell_result.json")
    isfile(source) ||
        error("Phase B reuse needs Phase A cell $(plan.reuse_cell_id), which has no checkpoint")
    src = D5OntAdaptiveGridDriver.load_row(source)
    src["status"] == D5OntAdaptiveGridDriver.STATUS_OK ||
        error("Phase A cell $(plan.reuse_cell_id) has status $(src["status"]); not reusable")
    metric_names = D5OntAdaptiveGridDriver.METRIC_KEYS
    metrics = NamedTuple{metric_names}(Tuple(src[String(k)] for k in metric_names))
    return D5OntAdaptiveGridDriver.cell_row(cell, context; phase = :B,
        status = D5OntAdaptiveGridDriver.STATUS_OK, k_selected = plan.k_star,
        reassembly_k = src["reassembly_k"], reassembly_k_policy = src["reassembly_k_policy"],
        reused_from_cell_id = plan.reuse_cell_id, n_contigs = src["n_contigs"],
        contigs_ge_min = src["contigs_ge_min"], quast_min_contig = src["quast_min_contig"],
        metrics, nga50_status = src["nga50_status"],
        corrected_survival_status = get(src, "corrected_survival_status", missing),
        wall_seconds = src["wall_seconds"], peak_rss_bytes = src["peak_rss_bytes"],
        peak_rss_method = src["peak_rss_method"], row_kwargs...)
end
# === Merge ===

function publish_tsv(path, df)
    mkpath(dirname(path))
    tmp = tempname(dirname(path))
    try
        CSV.write(tmp, df; delim = '\t', missingstring = "NA")
        mv(tmp, path; force = true)
    catch
        isfile(tmp) && rm(tmp; force = true)
        rethrow()
    end
    return path
end

function frame(dicts, columns)
    return DataFrames.DataFrame([Symbol(c) => [get(d, c, missing) for d in dicts]
                                 for c in columns])
end

nonmissing(xs) = collect(skipmissing(xs))
median_or_missing(xs) = (v = nonmissing(xs); isempty(v) ? missing : Statistics.median(v))
function cv_or_missing(xs)
    v = Float64.(nonmissing(xs))
    length(v) >= 2 || return missing
    m = Statistics.mean(v)
    m == 0 && return missing
    return Statistics.std(v; corrected = true) / m
end

"""
Load every checkpoint and add a row for every registered cell that lacks one.
A cell with no checkpoint gets an `absent` row; a checkpoint that cannot be
parsed gets an `unreadable` row (it exists, so it is not absent, and it is
counted once). Phase B is expected only once any Phase B checkpoint exists.
Rows are Dicts keyed by column name.
"""
function collect_merge_rows(root::AbstractString, toy::Bool)
    cells_dir = joinpath(root, "cells")
    rows = Dict{String, Any}[]
    unreadable = String[]
    for entry in sort(isdir(cells_dir) ? readdir(cells_dir) : String[])
        ckpt = joinpath(cells_dir, entry, "cell_result.json")
        isfile(ckpt) || continue
        try
            push!(rows, D5OntAdaptiveGridDriver.load_row(ckpt))
        catch e
            push!(unreadable, entry)
            @warn "unreadable checkpoint" entry exception = e
        end
    end
    # Phase B is expected once ANY Phase B checkpoint exists, readable or not:
    # a corrupt-only Phase B must not make Phase B look like it never started.
    has_b = any(r -> isequal(get(r, "phase", missing), "B"), rows) ||
            any(id -> occursin("__adaptive__", id), unreadable)
    loaded = Set(r["cell_id"] for r in rows)
    unreadable_set = Set(unreadable)
    as_row(nt) = Dict{String, Any}(String(k) => v for (k, v) in pairs(nt))
    n_absent = 0
    for c in expected_cells(toy, has_b)
        id = D5OntAdaptiveGridDriver.cell_id(c)
        id in loaded && continue
        phase = D5OntAdaptiveGridDriver.cell_phase(c)
        if id in unreadable_set
            push!(rows,
                as_row(D5OntAdaptiveGridDriver.cell_row(c,
                    D5OntAdaptiveGridDriver.EMPTY_CONTEXT; phase,
                    status = D5OntAdaptiveGridDriver.STATUS_UNREADABLE)))
        else
            n_absent += 1
            push!(rows, as_row(D5OntAdaptiveGridDriver.absent_row(c; phase)))
        end
    end
    return (; rows, n_absent, unreadable, has_b)
end

"""
Everything that makes a merged grid incomplete or internally inconsistent, as
COUNTS (no policy). Two different kinds, kept apart because they are waived
differently (see `merge_verdict`):

  * completeness: absent cells, unreadable checkpoints (every unreadable entry,
    registered or not, passed in as `n_unreadable_entries`), retryable cells.
  * consistency: more than one (tau, tau_hat) pair among Phase B rows; more
    than one Mycelia SHA among rows; rows written from a tree that was not
    clean (`git_dirty` not `false`).

Corrected-survival diagnostic failures are reported but are not problems
(the diagnostic never feeds k*).
"""
function merge_problems(rows; n_unreadable_entries::Int = 0)
    status(r) = get(r, "status", missing)
    count_status(s) = count(r -> isequal(status(r), s), rows)
    measured_rows = [r
                     for r in rows
                     if !(status(r) in (D5OntAdaptiveGridDriver.STATUS_ABSENT,
        D5OntAdaptiveGridDriver.STATUS_UNREADABLE))]
    n_retry = needs_retry(skipmissing(status.(rows)))
    tau_pairs = unique((r["tau"], r["tau_hat"])
    for r in measured_rows
    if isequal(get(r, "phase", missing), "B") && !ismissing(get(r, "tau", missing)))
    shas = unique(get(r, "mycelia_git_sha", missing) for r in measured_rows)
    n_dirty_rows = count(r -> !isequal(get(r, "git_dirty", missing), false), measured_rows)
    n_absent = count_status(D5OntAdaptiveGridDriver.STATUS_ABSENT)
    n_unreadable = max(n_unreadable_entries,
        count_status(D5OntAdaptiveGridDriver.STATUS_UNREADABLE))
    corrected(s) = count(r -> isequal(get(r, "corrected_survival_status", missing), s),
        rows)
    return (; n_absent, n_unreadable, n_retry,
        n_incomplete = n_absent + n_unreadable + n_retry,
        tau_pairs = [[t, th] for (t, th) in tau_pairs],
        n_tau_conflicts = max(0, length(tau_pairs) - 1),
        shas = collect(shas), n_sha_conflicts = max(0, length(shas) - 1),
        n_dirty_rows,
        n_corrected_survival_failed = corrected("failed"),
        n_corrected_survival_unavailable = corrected("unavailable"))
end

"""
What blocks a merge, given the problem counts and the operator's flags. Each
flag waives exactly one class, so accepting one never silently accepts another:

  * incomplete cells (absent / unreadable / retryable): `--allow-incomplete`
    (a partial view of a grid still running is legitimate);
  * rows written at more than one Mycelia SHA or from an unclean tree:
    `--allow-mixed-code` (distinct from --allow-sha-mismatch, which only lets
    THIS checkout differ from the manifest);
  * more than one (tau, tau_hat) pair: never waivable. A Phase B table is a
    function of tau (same rule as track_a_merge_hosts.jl's consistency check).

Returns `(; blocking, waived, blocking_reasons, waived_reasons)`.
"""
function merge_verdict(p; allow_incomplete::Bool, allow_mixed_code::Bool)
    blocking = 0
    waived = 0
    blocking_reasons = String[]
    waived_reasons = String[]
    function record(n, reason, waive)
        n > 0 || return
        if waive
            waived += n
            push!(waived_reasons, reason)
        else
            blocking += n
            push!(blocking_reasons, reason)
        end
    end
    record(p.n_incomplete,
        "$(p.n_incomplete) incomplete cell(s) ($(p.n_absent) absent, " *
        "$(p.n_unreadable) unreadable, $(p.n_retry) retryable)", allow_incomplete)
    record(p.n_sha_conflicts + p.n_dirty_rows,
        "mixed code versions ($(length(p.shas)) SHA(s), $(p.n_dirty_rows) row(s) " *
        "from an unclean tree)", allow_mixed_code)
    record(p.n_tau_conflicts,
        "Phase B rows from $(length(p.tau_pairs)) different (tau, tau_hat) pairs " *
        "$(p.tau_pairs); never waivable", false)
    return (; blocking, waived, blocking_reasons, waived_reasons)
end

function run_merge(cli, root)
    base = run_context_base()
    check_run_manifest!(root, git_state(base);
        allow_mismatch = cli.allow_sha_mismatch, toy = cli.toy)
    collected = collect_merge_rows(root, cli.toy)
    cells = frame(collected.rows, D5OntAdaptiveGridDriver.cell_row_columns())
    publish_tsv(joinpath(root, "d5_cells.tsv"), cells)
    survival_df = survival_table(root)
    publish_tsv(joinpath(root, "d5_survival.tsv"), survival_df)
    publish_tsv(joinpath(root, "illumina_none_optimum.tsv"), illumina_optimum(cells))
    publish_tsv(joinpath(root, "replicate_evaluability.tsv"), replicate_evaluability(cells))

    p = merge_problems(collected.rows; n_unreadable_entries = length(collected.unreadable))
    v = merge_verdict(p; cli.allow_incomplete, cli.allow_mixed_code)
    D5OntAdaptiveGridDriver.save_json_atomically(joinpath(root, "merge_provenance.json"),
        Dict{String, Any}("generated" => string(Dates.now()),
            "mycelia_git_sha" => base.mycelia_git_sha, "git_dirty" => base.git_dirty,
            "host" => base.host, "tool_versions" => tool_versions(),
            "overrides" => Dict("allow_sha_mismatch" => cli.allow_sha_mismatch,
                "allow_mixed_code" => cli.allow_mixed_code,
                "allow_incomplete" => cli.allow_incomplete),
            "n_rows" => DataFrames.nrow(cells),
            "n_expected" => length(expected_cells(cli.toy, collected.has_b)),
            "n_absent" => p.n_absent, "n_unreadable" => p.n_unreadable,
            "unreadable_checkpoints" => collected.unreadable, "n_retryable" => p.n_retry,
            "phase_b_tau_pairs" => p.tau_pairs, "row_mycelia_shas" => p.shas,
            "n_dirty_rows" => p.n_dirty_rows,
            "n_corrected_survival_failed" => p.n_corrected_survival_failed,
            "n_corrected_survival_unavailable" => p.n_corrected_survival_unavailable,
            "blocking_problems" => v.blocking, "waived_problems" => v.waived,
            "blocking_reasons" => v.blocking_reasons, "waived_reasons" => v.waived_reasons,
            "n_survival_rows" => DataFrames.nrow(survival_df),
            "invocation" => "julia --project=. benchmarking/d5_ont_adaptive_grid_run.jl " *
                            join(ARGS, " ")))
    n_ok = count(isequal("ok"), cells.status)
    println("merged $(DataFrames.nrow(cells)) rows ($(n_ok) ok) into $(root); " *
            "corrected-survival failed $(p.n_corrected_survival_failed), unavailable " *
            "$(p.n_corrected_survival_unavailable)")
    v.waived > 0 &&
        @warn "merge WAIVED $(v.waived) problem(s) by flag; the tables are not a " *
              "complete, single-version grid" waived = v.waived_reasons
    if v.blocking > 0
        @error "merge found $(v.blocking) blocking problem(s); the tables are written " *
               "but must not be used as the grid result" blocking = v.blocking_reasons
        return v.blocking
    end
    return 0
end
"""Every registered cell of this run (toy or full), Phase A plus optionally Phase B."""
function expected_cells(toy::Bool, include_phase_b::Bool)
    cells = NamedTuple[]
    for rs in selected_read_sets(toy)
        append!(cells, phase_cells(rs, :A, toy))
        include_phase_b && append!(cells, phase_cells(rs, :B, toy))
    end
    return cells
end

function survival_table(root)
    survival_rows = Dict{String, Any}[]
    readsets_dir = joinpath(root, "readsets")
    for entry in sort(isdir(readsets_dir) ? readdir(readsets_dir) : String[])
        marker = joinpath(readsets_dir, entry, "readset.json")
        isfile(marker) || continue
        meta = JSON.parsefile(marker)
        survival = JSON.parsefile(joinpath(readsets_dir, entry, "survival.json"))
        for r in survival["rungs"]
            push!(survival_rows,
                Dict{String, Any}(
                    "read_set_id" => entry, "reads" => "raw",
                    "technology" => meta["technology"],
                    "identity_target" => meta["identity_target"],
                    "identity_measured_mean" => meta["identity_measured_mean"],
                    "coverage" => meta["coverage"], "organism" => meta["organism"],
                    "seed" => meta["seed"], "k" => r["k"], "S" => r["S"],
                    "n_ref_positions" => r["n_ref_positions"],
                    "n_ref_positions_observed" => r["n_ref_positions_observed"],
                    "n_distinct_read_kmers" => r["n_distinct_read_kmers"],
                    "S_hat" => r["S_hat"], "s_hat_v" => r["s_hat_v"],
                    "s_hat_m" => r["s_hat_m"], "s_hat_status" => r["s_hat_status"],
                    "error_rate_estimate_diagnostic" => meta["error_rate_estimate_diagnostic"]))
        end
    end
    columns = ["read_set_id", "reads", "technology", "identity_target",
        "identity_measured_mean", "coverage", "organism", "seed", "k", "S",
        "n_ref_positions", "n_ref_positions_observed", "n_distinct_read_kmers", "S_hat",
        "s_hat_v", "s_hat_m", "s_hat_status", "error_rate_estimate_diagnostic"]
    return frame(
        [Dict(k => (v === nothing ? missing : v) for (k, v) in d)
         for d in survival_rows], columns)
end

"""
Illumina :none stratum, Phase A: median NGA50 per rung, so a human can read off
each organism's measured optimum when setting tau. Computes NO tau.

Every registered cell counts, whatever its status, so `n_measured` is out of
the registered seeds (x coverages for the pooled scope). Medians are over
MEASURED NGA50 only, and `is_optimum` is eligible only for rungs where every
registered replicate was measured: a rung that won on one surviving seed is
not an optimum.
"""
function illumina_optimum(cells;
        n_seeds::Int = length(
            D5OntAdaptiveGridDriver.D5OntAdaptiveGridAxes.D5_SEEDS),
        n_coverages::Int = length(D5OntAdaptiveGridDriver.D5_ILLUMINA_COVERAGES))
    calib = String(D5OntAdaptiveGridDriver.D5_ILLUMINA_CALIBRATION_CORRECTOR)
    df = DataFrames.nrow(cells) == 0 ? cells :
         cells[
    isequal.(
        cells.technology, "illumina") .& isequal.(cells.corrector, calib) .& isequal.(cells.phase, "A"), :]
    out = NamedTuple[]
    for scope in (:per_coverage, :all_coverages)
        DataFrames.nrow(df) == 0 && break
        group_cols = scope === :per_coverage ? [:organism, :coverage, :k_requested] :
                     [:organism, :k_requested]
        n_expected = scope === :per_coverage ? n_seeds : n_seeds * n_coverages
        for g in DataFrames.groupby(df, group_cols)
            measured = g[isequal.(g.nga50_status, "measured"), :]
            n_measured = DataFrames.nrow(measured)
            push!(out,
                (; scope = String(scope), organism = g.organism[1],
                    coverage = scope === :per_coverage ? string(g.coverage[1]) : "all",
                    k = g.k_requested[1], n_expected, n_rows = DataFrames.nrow(g),
                    n_measured, all_measured = n_measured == n_expected,
                    median_nga50 = median_or_missing(measured.NGA50),
                    median_genome_fraction = median_or_missing(
                        g[isequal.(g.status, "ok"), :].genome_fraction)))
        end
    end
    table = DataFrames.DataFrame(out)
    DataFrames.nrow(table) == 0 && return table
    table.is_optimum = falses(DataFrames.nrow(table))
    for g in DataFrames.groupby(table, [:scope, :organism, :coverage]; sort = false)
        vals = [r.all_measured ? coalesce(r.median_nga50, -Inf) : -Inf for r in eachrow(g)]
        best = maximum(vals)
        best == -Inf && continue
        g.is_optimum .= vals .== best
    end
    return table
end

"""
Replicate evaluability and CV per (technology, identity, coverage, corrector,
k-policy LEVEL, stratum, organism, phase), over seeds. Every registered cell
counts whatever its status (absent, error, infeasible, no_feasible_k), and the
denominator is the registered seed count, so a group with a failed or missing
seed can never read as fully measured. Adaptive cells group by policy level,
not by the k* each seed happened to select. The admissibility rule's constants
are not set here.
"""
function replicate_evaluability(cells;
        n_expected::Int = length(
            D5OntAdaptiveGridDriver.D5OntAdaptiveGridAxes.D5_SEEDS))
    out = NamedTuple[]
    DataFrames.nrow(cells) == 0 && return DataFrames.DataFrame(out)
    group_cols = [:technology, :identity_target, :coverage, :corrector, :policy_level,
        :stratum, :organism, :phase]
    status_count(g, s) = count(isequal(s), g.status)
    for g in DataFrames.groupby(cells, group_cols)
        n_measured = count(isequal("measured"), g.nga50_status)
        ok = g[isequal.(g.status, "ok"), :]
        push!(out,
            (; (k => g[1, k] for k in group_cols)..., n_expected,
                n_rows = DataFrames.nrow(g), n_ok = DataFrames.nrow(ok), n_measured,
                n_error = status_count(g, D5OntAdaptiveGridDriver.STATUS_ERROR),
                n_quast_failed = status_count(g, D5OntAdaptiveGridDriver.STATUS_QUAST_FAILED),
                n_absent = status_count(g, D5OntAdaptiveGridDriver.STATUS_ABSENT),
                n_infeasible = status_count(g, D5OntAdaptiveGridDriver.STATUS_INFEASIBLE_K),
                n_no_feasible_k = status_count(g, D5OntAdaptiveGridDriver.STATUS_NO_FEASIBLE_K),
                all_replicates_measured = n_measured == n_expected,
                cv_nga50 = cv_or_missing(g[isequal.(g.nga50_status, "measured"), :].NGA50),
                cv_genome_fraction = cv_or_missing(ok.genome_fraction)))
    end
    return DataFrames.DataFrame(out)
end
function print_counts(toy::Bool)
    a = D5OntAdaptiveGridDriver.phase_a_cells()
    b = D5OntAdaptiveGridDriver.phase_b_cells()
    ont(cs) = count(c -> c.technology === :ont, cs)
    println("registered ONT cells (n_d5_cells): ",
        D5OntAdaptiveGridDriver.D5OntAdaptiveGridAxes.n_d5_cells())
    println("Phase A: $(length(a)) = $(ont(a)) ONT + $(length(a) - ont(a)) Illumina")
    println("Phase B: $(length(b)) = $(ont(b)) ONT + $(length(b) - ont(b)) Illumina")
    infeasible = count(c -> !D5OntAdaptiveGridDriver.assembler_accepts_k(c.k_policy), a)
    println("Phase A cells above the assembler's k limit (recorded, not assembled): ",
        infeasible)
    sets = selected_read_sets(toy)
    println("read sets$(toy ? " (toy)" : ""): $(length(sets)) -> --array=0-$(length(sets) - 1)")
    if toy
        for rs in sets
            println("  $(D5OntAdaptiveGridDriver.read_set_id(rs)): ",
                join(
                    D5OntAdaptiveGridDriver.cell_id.(D5OntAdaptiveGridDriver.toy_cells(rs)),
                    ", "))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    cli = D5OntAdaptiveGridDriver.parse_cli(ARGS)
    println("=== D5 grid runner: phase $(cli.phase)$(cli.toy ? " (toy)" : "")  $(Dates.now()) ===")
    if cli.phase === :list
        print_counts(cli.toy)
    else
        cli.output_dir === nothing &&
            error("--output-dir is required for --phase $(cli.phase)")
        mkpath(abspath(cli.output_dir))
        # realpath: cell output dirs (corrector output_dir) live under this root, and
        # the corrector refuses any working root that resolves through a symlink.
        root = realpath(abspath(cli.output_dir))
        # The corrector's temp roots follow TMPDIR, not --output-dir, and it refuses
        # a root through a symlink. Pin TMPDIR to a real path under the grid root for
        # EVERY invocation (sbatch, login node, Lovelace), not only the sbatch.
        tmp = joinpath(root, "tmp")
        mkpath(tmp)
        ENV["TMPDIR"] = tmp
        # Cell failures are caught and checkpointed as rows, so the process would
        # otherwise end normally. Exit non-zero whenever a resubmit is needed, so
        # SLURM shows FAILED and the operator sees it.
        n_retry = if cli.phase === :A
            needs_retry(run_phase_a(cli, root))
        elseif cli.phase === :B
            needs_retry(run_phase_b(cli, root))
        else
            run_merge(cli, root)
        end
        if n_retry > 0
            println("=== $(n_retry) problem(s): ",
                cli.phase === :merge ? "see merge_provenance.json" :
                "cells need a resubmit of this index",
                "; exiting 3  $(Dates.now()) ===")
            exit(3)
        end
    end
    println("=== done $(Dates.now()) ===")
end
