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

"""
One grid, one Mycelia SHA. The first writer records the SHA; every later task
(Phase A shards, merge, Phase B) must match it unless explicitly overridden.
This keeps a phased run a single grid execution (PLAN-2026-03-28 step 2).
"""
function check_run_manifest!(root::AbstractString, sha::AbstractString;
        allow_mismatch::Bool, toy::Bool)
    path = joinpath(root, "run_manifest.json")
    if !isfile(path)
        D5OntAdaptiveGridDriver.save_json_atomically(path,
            Dict("mycelia_git_sha" => sha, "toy" => toy,
                "created" => string(Dates.now()),
                "badread_identity_params" => D5OntAdaptiveGridDriver.D5_BADREAD_IDENTITY_PARAMS,
                "badread_error_model" => D5OntAdaptiveGridDriver.D5_BADREAD_ERROR_MODEL,
                "badread_qscore_model" => D5OntAdaptiveGridDriver.D5_BADREAD_QSCORE_MODEL,
                "illumina_params" => Dict(String(k) => v
                for (k, v) in pairs(D5OntAdaptiveGridDriver.D5_ILLUMINA_PARAMS)),
                "strategy" => String(D5OntAdaptiveGridDriver.D5_STRATEGY)))
    end
    manifest = JSON.parsefile(path)
    manifest["toy"] == toy ||
        error("output dir $(root) belongs to a $(manifest["toy"] ? "toy" : "full") run; " *
              "refusing to mix toy and full cells in one tree")
    recorded = manifest["mycelia_git_sha"]
    if recorded != sha
        msg = "run_manifest.json records Mycelia $(recorded) but this checkout is $(sha). " *
              "A grid is one code version; pass --allow-sha-mismatch only if the " *
              "difference cannot touch assembly output, and say so in the results."
        allow_mismatch ? (@warn msg) : error(msg)
    end
    return manifest
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
    abs(genome_size - info.expected_bp) <= 0.2 * info.expected_bp ||
        error("reference $(info.accession) is $(genome_size) bp; expected ~$(info.expected_bp)")
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
        records = reduce(vcat, read_records.(meta["reads_paths"]))
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
        argv = D5OntAdaptiveGridDriver.badread_argv(rs, ref.path)
        tmp = joinpath(work, "reads.fq.gz")
        Mycelia.simulate_nanopore_reads(; fasta = ref.path,
            quantity = "$(rs.coverage)x",
            error_model = D5OntAdaptiveGridDriver.D5_BADREAD_ERROR_MODEL,
            qscore_model = D5OntAdaptiveGridDriver.D5_BADREAD_QSCORE_MODEL,
            identity = D5OntAdaptiveGridDriver.D5_BADREAD_IDENTITY_PARAMS[rs.identity],
            length_dist = D5OntAdaptiveGridDriver.D5_BADREAD_LENGTH_DIST,
            seed = rs.seed, outfile = tmp, quiet = true)
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
    if stats.n_ge_min > 0 && filesize(contigs_path) > 0
        try
            Mycelia.run_quast([contigs_path]; outdir = quast_dir, reference = ref.path,
                min_contig, threads = QUAST_THREADS)
        catch e
            @warn "QUAST failed" cell_dir exception = (e, catch_backtrace())
        end
    end
    quast_ran = isfile(report)
    metrics = quast_ran ? D5OntAdaptiveGridDriver.parse_quast_metrics(report) :
              D5OntAdaptiveGridDriver.empty_metrics()
    nga50_status = D5OntAdaptiveGridDriver.nga50_status_for(metrics, stats.n,
        stats.n_ge_min, quast_ran)
    status = nga50_status == "quast_failed" ? D5OntAdaptiveGridDriver.STATUS_QUAST_FAILED :
             D5OntAdaptiveGridDriver.STATUS_OK

    if iterative
        corrected_survival_diagnostic(asm_dir, ref, cell_dir)
    end

    return D5OntAdaptiveGridDriver.cell_row(cell, context; phase, status,
        k_selected = k, reassembly_k,
        reassembly_k_policy = iterative ? "fixed" : missing,
        n_contigs = stats.n, contigs_ge_min = stats.n_ge_min,
        quast_min_contig = min_contig, metrics, nga50_status,
        wall_seconds = round(measured.wall_seconds; digits = 3),
        peak_rss_bytes = measured.peak_rss_bytes, peak_rss_method = measured.method,
        row_kwargs...)
end

"""
DIAGNOSTIC: survival on this cell's own corrected reads (user decision
2026-11-16). Never feeds k*. Absent corrected FASTQ is recorded, not fatal.
"""
function corrected_survival_diagnostic(asm_dir, ref, cell_dir)
    corrected = joinpath(asm_dir, "corrected.fastq")
    out = joinpath(cell_dir, "corrected_survival.json")
    if !isfile(corrected)
        @warn "corrected FASTQ not found; corrected-read survival unavailable" corrected maxlog = 1
        D5OntAdaptiveGridDriver.save_json_atomically(out,
            Dict("status" => "unavailable", "role" => "diagnostic"))
        return nothing
    end
    seqs = read_sequences(corrected)
    rows = D5OntAdaptiveGridDriver.survival_ladder(seqs, ref.seqs)
    D5OntAdaptiveGridDriver.save_json_atomically(out,
        Dict("status" => "measured", "role" => "diagnostic", "reads" => "corrected",
            "strand" => "canonical", "rungs" => survival_json(rows)))
    rm(corrected; force = true)
    return nothing
end

function checkpoint_path(root, cell)
    return joinpath(root, "cells", D5OntAdaptiveGridDriver.cell_id(cell), "cell_result.json")
end

"""Skip a completed cell; retry error / quast_failed rows; record failures as rows."""
function run_or_resume(f, root, cell, context; phase::Symbol)
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
            println("  cached  $(D5OntAdaptiveGridDriver.cell_id(cell)): $(prior["status"])")
            return prior
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
            error_message = first(sprint(showerror, e), 2000))
    end
    D5OntAdaptiveGridDriver.save_json_atomically(ckpt, D5OntAdaptiveGridDriver.row_to_json(row))
    println("  $(row.status)  $(row.cell_id)  nga50_status=$(row.nga50_status) " *
            "NGA50=$(row.NGA50) GF=$(row.genome_fraction)")
    return row
end

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
    check_run_manifest!(root, base.mycelia_git_sha;
        allow_mismatch = cli.allow_sha_mismatch, toy = cli.toy)
    for rs in read_sets_for(cli)
        println("=== read set $(D5OntAdaptiveGridDriver.read_set_id(rs)) ===")
        ref = prepare_reference(rs.organism, joinpath(root, "refs"))
        prepared = prepare_read_set(rs, ref, root)
        context = cell_context(prepared.meta, base)
        for cell in phase_cells(rs, :A, cli.toy)
            run_or_resume(root, cell, context; phase = :A) do cell_dir
                assemble_cell(cell, Int(cell.k_policy), prepared.records, ref, cell_dir,
                    context; phase = :A)
            end
        end
    end
end

function run_phase_b(cli, root)
    base = run_context_base()
    check_run_manifest!(root, base.mycelia_git_sha;
        allow_mismatch = cli.allow_sha_mismatch, toy = cli.toy)
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
            plan_kwargs = (; k_star_reason = plan.reason,
                k_star_hat = plan.k_star_hat === D5OntAdaptiveGridDriver.NO_FEASIBLE_K ?
                             String(plan.k_star_hat) : plan.k_star_hat)
            run_or_resume(root, cell, context; phase = :B) do cell_dir
                if plan.action === :no_feasible_k
                    D5OntAdaptiveGridDriver.cell_row(cell, context; phase = :B,
                        status = D5OntAdaptiveGridDriver.STATUS_NO_FEASIBLE_K,
                        k_selected = String(D5OntAdaptiveGridDriver.NO_FEASIBLE_K),
                        plan_kwargs...)
                elseif plan.action === :infeasible
                    D5OntAdaptiveGridDriver.infeasible_row(cell, context; phase = :B,
                        k_selected = plan.k_star, plan_kwargs...)
                elseif plan.action === :reuse
                    reuse_row(root, cell, plan, context, plan_kwargs)
                else
                    records === nothing &&
                        (records = reduce(vcat, read_records.(meta["reads_paths"])))
                    assemble_cell(cell, plan.k_star, records, ref, cell_dir, context;
                        phase = :B, plan_kwargs...)
                end
            end
        end
    end
end

"""Phase B row copied from the Phase A comparator assembly at k*, recorded as reuse."""
function reuse_row(root, cell, plan, context, plan_kwargs)
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
        metrics, nga50_status = src["nga50_status"], wall_seconds = src["wall_seconds"],
        peak_rss_bytes = src["peak_rss_bytes"], peak_rss_method = src["peak_rss_method"],
        plan_kwargs...)
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

function run_merge(cli, root)
    base = run_context_base()
    check_run_manifest!(root, base.mycelia_git_sha;
        allow_mismatch = cli.allow_sha_mismatch, toy = cli.toy)
    cells_dir = joinpath(root, "cells")
    rows = Dict{String, Any}[]
    skipped = String[]
    for entry in sort(isdir(cells_dir) ? readdir(cells_dir) : String[])
        ckpt = joinpath(cells_dir, entry, "cell_result.json")
        isfile(ckpt) || continue
        try
            push!(rows, D5OntAdaptiveGridDriver.load_row(ckpt))
        catch e
            push!(skipped, entry)
            @warn "skipping unreadable checkpoint" entry exception = e
        end
    end
    isempty(skipped) ||
        @warn "merge omits $(length(skipped)) unreadable checkpoint(s)" skipped
    cells = frame(rows, D5OntAdaptiveGridDriver.cell_row_columns())
    publish_tsv(joinpath(root, "d5_cells.tsv"), cells)

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
                    "read_set_id" => entry, "technology" => meta["technology"],
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
    survival_columns = ["read_set_id", "technology", "identity_target",
        "identity_measured_mean", "coverage", "organism", "seed", "k", "S",
        "n_ref_positions", "n_ref_positions_observed", "n_distinct_read_kmers", "S_hat",
        "s_hat_v", "s_hat_m", "s_hat_status", "error_rate_estimate_diagnostic"]
    survival_df = frame(
        [Dict(k => (v === nothing ? missing : v) for (k, v) in d)
         for d in survival_rows],
        survival_columns)
    publish_tsv(joinpath(root, "d5_survival.tsv"), survival_df)

    write_illumina_optimum(root, cells)
    write_replicate_evaluability(root, cells)
    write_merge_provenance(root, cells, survival_df, base)
    n_ok = isempty(rows) ? 0 : count(isequal("ok"), cells.status)
    println("merged $(length(rows)) cells ($(n_ok) ok), " *
            "$(DataFrames.nrow(survival_df)) survival rows into $(root)")
end

"""
Illumina :none stratum, Phase A: median NGA50 per rung, so a human can read off
each organism's measured optimum when setting tau. This step computes NO tau.
Medians are over MEASURED NGA50 only; `n_measured` shows how many contributed.
"""
function write_illumina_optimum(root, cells)
    calib = String(D5OntAdaptiveGridDriver.D5_ILLUMINA_CALIBRATION_CORRECTOR)
    df = isempty(cells) ? cells :
         cells[(isequal.(cells.technology, "illumina")) .& (isequal.(cells.corrector, calib)) .& (isequal.(cells.phase, "A")) .& (isequal.(cells.status, "ok")), :]
    out = NamedTuple[]
    for scope in (:per_coverage, :all_coverages)
        group_cols = scope === :per_coverage ? [:organism, :coverage, :k_requested] :
               [:organism, :k_requested]
        isempty(df) && break
        for g in DataFrames.groupby(df, group_cols)
            measured = g[isequal.(g.nga50_status, "measured"), :]
            push!(out,
                (; scope = String(scope), organism = g.organism[1],
                    coverage = scope === :per_coverage ? string(g.coverage[1]) : "all",
                    k = g.k_requested[1], n = DataFrames.nrow(g),
                    n_measured = DataFrames.nrow(measured),
                    median_nga50 = median_or_missing(measured.NGA50),
                    median_genome_fraction = median_or_missing(g.genome_fraction)))
        end
    end
    table = DataFrames.DataFrame(out)
    if !isempty(table)
        table.is_optimum = falses(DataFrames.nrow(table))
        for g in DataFrames.groupby(table, [:scope, :organism, :coverage]; sort = false)
            vals = coalesce.(g.median_nga50, -Inf)
            maximum(vals) == -Inf && continue
            g.is_optimum .= vals .== maximum(vals)
        end
    end
    publish_tsv(joinpath(root, "illumina_none_optimum.tsv"), table)
end

"""
Replicate evaluability and CV per (technology, identity, coverage, corrector,
k-policy level, organism), over seeds. `all_replicates_measured` is the NGA50
evaluability the admissibility rule consumes; its constants are not set here.
"""
function write_replicate_evaluability(root, cells)
    out = NamedTuple[]
    if !isempty(cells)
        df = cells[isequal.(cells.status, "ok"), :]
        group_cols = [:technology, :identity_target, :coverage, :corrector, :k_policy,
            :k_selected, :stratum, :organism, :phase]
        for g in DataFrames.groupby(df, group_cols)
            n_measured = count(isequal("measured"), g.nga50_status)
            push!(out,
                (; (k => g[1, k] for k in group_cols)..., n = DataFrames.nrow(g), n_measured,
                    all_replicates_measured = n_measured == DataFrames.nrow(g),
                    cv_nga50 = cv_or_missing(g.NGA50),
                    cv_genome_fraction = cv_or_missing(g.genome_fraction)))
        end
    end
    publish_tsv(joinpath(root, "replicate_evaluability.tsv"), DataFrames.DataFrame(out))
end

function write_merge_provenance(root, cells, survival_df, base)
    D5OntAdaptiveGridDriver.save_json_atomically(joinpath(root, "merge_provenance.json"),
        Dict("generated" => string(Dates.now()), "mycelia_git_sha" => base.mycelia_git_sha,
            "git_dirty" => base.git_dirty, "host" => base.host,
            "tool_versions" => tool_versions(), "n_cells" => DataFrames.nrow(cells),
            "n_survival_rows" => DataFrames.nrow(survival_df),
            "invocation" => "julia --project=. benchmarking/d5_ont_adaptive_grid_run.jl " *
                            join(ARGS, " ")))
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
        root = abspath(cli.output_dir)
        mkpath(root)
        cli.phase === :A && run_phase_a(cli, root)
        cli.phase === :B && run_phase_b(cli, root)
        cli.phase === :merge && run_merge(cli, root)
    end
    println("=== done $(Dates.now()) ===")
end
