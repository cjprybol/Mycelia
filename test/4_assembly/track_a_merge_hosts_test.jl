# Unit test for the Track A cross-host checkpoint merge (bead td-bblmi).
#
# The merge's whole value is that it refuses to produce a plausible-looking
# 288-cell matrix when the inputs do not justify one. So the tests are mostly
# about the refusals:
#
#   * the same cell id with DIFFERENT content across hosts is a collision, is
#     reported with a field-level diff, is EXCLUDED from the merged matrix, and
#     fails the exit status
#   * the same cell id with IDENTICAL content is benign and reported separately
#   * missing cells are enumerated exactly, not just counted
#   * per-cell host provenance is recorded
#   * re-running is idempotent, and a previously merged cell that disagrees with
#     its source is a collision rather than an overwrite
#   * the expected-matrix constants mirrored from the harness have not drifted
#
# Dependency-free apart from JSON / SHA / DataFrames / CSV (all direct Mycelia
# deps that load without importing Mycelia). No network, no QUAST, no assembly.

import Test
import JSON
import DataFrames
import CSV

include(joinpath(@__DIR__, "..", "..", "benchmarking", "track_a_merge_hosts.jl"))

function _tam_throws_with(f, needles::Vector{String})
    try
        f()
    catch e
        msg = sprint(showerror, e)
        for needle in needles
            Test.@test occursin(needle, msg)
        end
        return msg
    end
    Test.@test false
    return ""
end

# A checkpoint shaped exactly like track_a_baseline_benchmark.jl's `cell_row`.
function _tam_cell(; organism = "Lambda", technology = "illumina", coverage = 30,
        seed = 42, arm = "qualmer", accession = "NC_001416", k = 31,
        n_reads = 1000, n_contigs = 5, NGA50 = 1316.0, misassemblies = 0.0,
        genome_fraction = 98.5, duplication_ratio = 1.0, largest_contig = 48000,
        wall_seconds = 12.5, peak_rss_bytes = 1_000_000, status = "ok")
    return Dict{String, Any}(
        "organism" => organism, "accession" => accession, "technology" => technology,
        "coverage" => coverage, "seed" => seed, "decoder_arm" => arm, "k" => k,
        "n_reads" => n_reads, "n_contigs" => n_contigs, "NGA50" => NGA50,
        "misassemblies" => misassemblies, "genome_fraction" => genome_fraction,
        "duplication_ratio" => duplication_ratio, "largest_contig" => largest_contig,
        "wall_seconds" => wall_seconds, "peak_rss_bytes" => peak_rss_bytes,
        "status" => status)
end

function _tam_write_cell(root, cell_id, cell)
    dir = joinpath(root, "cells", cell_id)
    mkpath(dir)
    open(joinpath(dir, "cell_result.json"), "w") do io
        JSON.print(io, cell, 2)
    end
    return joinpath(dir, "cell_result.json")
end

Test.@testset "Track A cross-host merge (td-bblmi)" begin
    Test.@testset "expected matrix mirrors the harness (drift guard)" begin
        # These constants are DUPLICATED from track_a_baseline_benchmark.jl, which
        # cannot be included (it runs the benchmark at include time). If the
        # harness changes its matrix or its cell-id format, "missing" would be
        # computed against the wrong set — so assert the mirror against the
        # harness SOURCE rather than trusting the duplication.
        harness = read(
            joinpath(@__DIR__, "..", "..", "benchmarking",
                "track_a_baseline_benchmark.jl"), String)
        # The harness now builds ids via `cell_id_for`, which appends `__k$(k)` only
        # for a non-default k so the pre---k trees stay resumable. Assert BOTH branches
        # of that expression, so this guard still fires if either the base format or the
        # k rule drifts away from track_a_cell_id above.
        Test.@test occursin(
            "\"\$(org)__\$(tech)__\$(cov)x__seed\$(seed)__\$(arm)\"", harness)
        Test.@test occursin(
            "\"\$(org)__\$(tech)__\$(cov)x__seed\$(seed)__\$(arm)__k\$(k)\"", harness)
        Test.@test occursin("cell_id = cell_id_for(org, tech, cov, seed, arm)", harness)

        # The ROW SCHEMA is mirrored too, and it silently fell two columns behind for
        # four review rounds: peak_rss_method and rss_baseline_bytes were added to the
        # harness so a reader can tell a sampled per-cell peak from a high-water delta,
        # and the merge kept emitting the 17-column schema — so the cross-host table,
        # the ONE table where the two hosts' methods actually mix, carried the value
        # and dropped the provenance that makes it interpretable.
        #
        # Nothing caught it because the only assertions on this constant compare the
        # merge's OUTPUT against the merge's OWN constant (below, and again in the
        # merged-table testsets). That is self-consistent by construction: it is a
        # tautology, not a drift guard. This compares against the harness instead.
        schema_ok, harness_row_keys = track_a_row_keys_match_harness()
        Test.@test schema_ok
        Test.@test harness_row_keys == TRACK_A_ROW_KEYS
        # Fails closed: an unreadable or unparseable harness must not read as agreement.
        Test.@test !isempty(harness_row_keys)
        # Absent provenance values must land as the HARNESS sentinels, not `missing`.
        # Mirroring the column NAMES was not enough: `get(cell, k, missing)` filled
        # them with `missing` for every pre-schema checkpoint (all 432 on disk), so
        # the merged table carried the names and no usable values, and the operation
        # the harness docstring prescribes — filter on peak_rss_method before
        # aggregating — threw ArgumentError on Missing. _tam_cell omits both keys, so
        # every other merge test runs this exact path and asserts nothing about it.
        # Compare against the HARNESS SOURCE, not against the literals written into
        # the constant in the same commit. The first version of these two assertions
        # did the latter, which cannot fail: drifting the harness sentinel to
        # "not-recorded"/-2 left this suite fully green while the two files disagreed
        # about exactly the columns the mirror exists to keep in step.
        defaults_ok, harness_defaults = track_a_absent_defaults_match_harness()
        Test.@test defaults_ok
        Test.@test harness_defaults == TRACK_A_ABSENT_DEFAULTS
        # Unparseable values are a MISMATCH, not a pass: the first version of this
        # guard silently DROPPED any pair it could not pattern-match, so a harness
        # gaining a third optional column this file does not mirror still compared
        # equal — the original defect, restated. 11/11 drifted harnesses read green.
        Test.@test harness_defaults == TRACK_A_ABSENT_DEFAULTS

        # The two columns whose absence was the defect.
        Test.@test "peak_rss_method" in TRACK_A_ROW_KEYS
        Test.@test "rss_baseline_bytes" in TRACK_A_ROW_KEYS

        # And the mirror holds for a non-default k in both directions.
        Test.@test track_a_cell_id("T4", "pacbio", 50, 123, "kmer"; k = 19) ==
                   "T4__pacbio__50x__seed123__kmer__k19"
        Test.@test track_a_cell_id("T4", "pacbio", 50, 123, "kmer"; k = 31) ==
                   track_a_cell_id("T4", "pacbio", 50, 123, "kmer")
        for org in TRACK_A_ORGANISMS
            Test.@test occursin("(\"$org\", ", harness)
        end
        Test.@test occursin(
            "const TECHNOLOGIES = [" *
            join(("\"$t\"" for t in TRACK_A_TECHNOLOGIES), ", ") * "]", harness)
        Test.@test occursin(
            "const COVERAGES = [" *
            join(string.(TRACK_A_COVERAGES), ", ") * "]", harness)
        Test.@test occursin(
            "const SEEDS = [" * join(string.(TRACK_A_SEEDS), ", ") * "]", harness)
        Test.@test occursin(
            "const DECODER_ARMS = [" *
            join(("\"$a\"" for a in TRACK_A_ARMS), ", ") * "]", harness)

        ids = expected_cell_ids()
        Test.@test length(ids) == 288           # the documented matrix size
        Test.@test length(unique(ids)) == 288
        Test.@test "Lambda__illumina__30x__seed42__qualmer" in ids
        Test.@test "SARS-CoV-2__ont__100x__seed456__kmer" in ids
        Test.@test track_a_cell_id("T4", "pacbio", 50, 123, "kmer") ==
                   "T4__pacbio__50x__seed123__kmer"
        # The real shard boundary the merge exists to reconcile.
        Test.@test length(expected_cell_ids(organisms = ("phi29", "SARS-CoV-2"))) == 144
    end

    Test.@testset "canonicalization: representation drift is not a collision" begin
        a = _tam_cell(n_contigs = 5, NGA50 = 1316.0)
        b = _tam_cell(n_contigs = 5.0, NGA50 = 1316)   # int/float swap only
        Test.@test cell_digest(a) == cell_digest(b)
        Test.@test isempty(differing_fields(a, b))
        # A genuine value difference always is a collision.
        c = _tam_cell(NGA50 = 1402.0)
        Test.@test cell_digest(a) != cell_digest(c)
        Test.@test differing_fields(a, c) == ["NGA50: 1316 vs 1402"]  # a's value first
        # Key-set differences are surfaced too.
        d = copy(a)
        delete!(d, "misassemblies")
        Test.@test any(occursin("misassemblies", f) for f in differing_fields(a, d))
        Test.@test any(occursin("<absent>", f) for f in differing_fields(a, d))
    end

    Test.@testset "quast_evidence is sound, not a guess" begin
        # QUAST scored it: a real Largest contig.
        Test.@test quast_evidence(_tam_cell(n_contigs = 5, largest_contig = 48_000)) ==
                   "quast:scored"
        # QUAST failed: contigs exist but empty_metrics() zeroed everything. A
        # scored run cannot report Largest contig = 0 for a non-empty assembly.
        Test.@test quast_evidence(_tam_cell(n_contigs = 5, largest_contig = 0,
            NGA50 = 0.0, genome_fraction = 0.0)) == "unknown:quast-unscored"
        # Nothing to score.
        Test.@test quast_evidence(_tam_cell(n_contigs = 0, largest_contig = 0,
            status = "empty_assembly")) == "n/a:empty-assembly"
        Test.@test quast_evidence(_tam_cell(status = "error", n_contigs = 0)) ==
                   "n/a:cell-error"
        # The signal is largest_contig, NOT "all metrics zero": a genuinely
        # unalignable assembly has genome_fraction 0 with a real largest_contig and
        # must still count as scored.
        Test.@test quast_evidence(_tam_cell(n_contigs = 3, largest_contig = 11_622,
            NGA50 = 0.0, genome_fraction = 0.0, misassemblies = 0.0)) == "quast:scored"
    end

    Test.@testset "disjoint shards merge cleanly with host provenance" begin
        mktempdir() do dir
            lovelace = joinpath(dir, "lovelace")
            lrc = joinpath(dir, "lrc")
            ids = expected_cell_ids(organisms = ("Lambda", "phi29"),
                technologies = ("illumina",), coverages = (30,), seeds = (42,))
            lambda_ids = filter(id -> startswith(id, "Lambda"), ids)
            phi29_ids = filter(id -> startswith(id, "phi29"), ids)
            for id in lambda_ids
                _tam_write_cell(lovelace, id, _tam_cell(organism = "Lambda"))
            end
            for id in phi29_ids
                _tam_write_cell(lrc,
                    id,
                    _tam_cell(organism = "phi29",
                        accession = "NC_011048", largest_contig = 19_000))
            end

            result = merge_hosts(["lovelace" => lovelace, "lrc" => lrc];
                expected_ids = ids)
            Test.@test length(result.merged) == length(ids)
            Test.@test isempty(result.collisions)
            Test.@test isempty(result.duplicates)
            Test.@test isempty(result.missing_ids)
            Test.@test isempty(result.unexpected_ids)
            # Provenance: every cell attributed to the host that produced it.
            for id in lambda_ids
                Test.@test result.origin[id] == "lovelace"
            end
            for id in phi29_ids
                Test.@test result.origin[id] == "lrc"
            end
            Test.@test all(haskey(result.digests, id) for id in ids)

            code, problems = merge_exit_status(result)
            Test.@test code == 0
            Test.@test isempty(problems)
        end
    end

    Test.@testset "CONTENT COLLISION: detected, diffed, excluded, and fatal" begin
        mktempdir() do dir
            a = joinpath(dir, "hostA")
            b = joinpath(dir, "hostB")
            cid = "Lambda__illumina__30x__seed42__qualmer"
            other = "Lambda__illumina__30x__seed42__kmer"
            _tam_write_cell(a, cid, _tam_cell(NGA50 = 1316.0, n_contigs = 5))
            # Same cell id, different result: the shards were supposed to be
            # disjoint, so this is a real anomaly (stale checkpoint / two SHAs).
            _tam_write_cell(b, cid, _tam_cell(NGA50 = 1402.0, n_contigs = 7))
            _tam_write_cell(a, other, _tam_cell())

            ids = [cid, other]
            result = merge_hosts(["hostA" => a, "hostB" => b]; expected_ids = ids)

            Test.@test length(result.collisions) == 1
            col = result.collisions[1]
            Test.@test col.cell_id == cid
            Test.@test sort(col.hosts) == ["hostA", "hostB"]
            Test.@test length(col.pairs) == 1
            diffs = col.pairs[1].diffs
            Test.@test any(occursin("NGA50", d) for d in diffs)
            Test.@test any(occursin("n_contigs", d) for d in diffs)
            # NOT last-write-wins: there is no defensible winner, so the cell does
            # not enter the analyzed matrix at all...
            Test.@test !haskey(result.merged, cid)
            # ...which means it also shows up as missing. Reported twice on purpose.
            Test.@test result.missing_ids == [cid]
            Test.@test haskey(result.merged, other)

            code, problems = merge_exit_status(result)
            Test.@test code == 1
            Test.@test Set(p.kind for p in problems) == Set([:collisions, :incomplete])
            Test.@test all(p -> p.fatal, problems)
            # Waivers are visible, not silent.
            code2,
            problems2 = merge_exit_status(result;
                allow_collisions = true, allow_incomplete = true)
            Test.@test code2 == 0
            Test.@test length(problems2) == 2
            Test.@test all(p -> !p.fatal, problems2)
            # Waiving only one is not enough.
            Test.@test first(merge_exit_status(result; allow_collisions = true)) == 1
        end
    end

    Test.@testset "identical duplicate across hosts is benign and reported" begin
        mktempdir() do dir
            a = joinpath(dir, "hostA")
            b = joinpath(dir, "hostB")
            cid = "T4__ont__100x__seed456__kmer"
            cell = _tam_cell(organism = "T4", technology = "ont", coverage = 100,
                seed = 456, arm = "kmer")
            _tam_write_cell(a, cid, cell)
            _tam_write_cell(b, cid, cell)

            result = merge_hosts(["hostA" => a, "hostB" => b]; expected_ids = [cid])
            Test.@test isempty(result.collisions)
            Test.@test length(result.duplicates) == 1
            Test.@test result.duplicates[1].cell_id == cid
            Test.@test sort(result.duplicates[1].hosts) == ["hostA", "hostB"]
            Test.@test haskey(result.merged, cid)
            Test.@test isempty(result.missing_ids)
            Test.@test first(merge_exit_status(result)) == 0
        end
    end

    Test.@testset "incompleteness enumerates exactly which cells are missing" begin
        mktempdir() do dir
            host = joinpath(dir, "lovelace")
            ids = expected_cell_ids(organisms = ("Lambda",),
                technologies = ("illumina",), coverages = (10, 30), seeds = (42,))
            Test.@test length(ids) == 4
            present = ids[1:2]
            for id in present
                _tam_write_cell(host, id, _tam_cell())
            end
            result = merge_hosts(["lovelace" => host]; expected_ids = ids)
            Test.@test sort(result.missing_ids) == sort(ids[3:4])
            Test.@test result.expected_n == 4
            code, problems = merge_exit_status(result)
            Test.@test code == 1
            Test.@test occursin("2 of 4 cells missing", problems[1].message)
        end
    end

    Test.@testset "in-progress and unparseable checkpoints are reported, not dropped" begin
        mktempdir() do dir
            host = joinpath(dir, "lrc")
            good = "phi29__illumina__30x__seed42__kmer"
            running = "phi29__illumina__30x__seed42__qualmer"
            broken = "phi29__illumina__50x__seed42__kmer"
            _tam_write_cell(host, good, _tam_cell(organism = "phi29"))
            mkpath(joinpath(host, "cells", running))            # no cell_result.json yet
            mkpath(joinpath(host, "cells", broken))
            write(joinpath(host, "cells", broken, "cell_result.json"), "{not json")

            result = merge_hosts(["lrc" => host]; expected_ids = [good, running, broken])
            Test.@test length(result.merged) == 1
            unreadable = result.per_host[1].unreadable
            Test.@test length(unreadable) == 2
            byid = Dict(unreadable)
            Test.@test occursin("no cell_result.json", byid[running])
            Test.@test occursin("unparseable", byid[broken])
            # Both are genuinely absent from the matrix.
            Test.@test sort(result.missing_ids) == sort([broken, running])
        end
    end

    Test.@testset "cell ids outside the expected matrix are flagged" begin
        mktempdir() do dir
            host = joinpath(dir, "h")
            cid = "Lambda__illumina__30x__seed42__qualmer"
            bogus = "Lambda__illumina__30x__seed99__qualmer"  # seed 99 not pre-registered
            _tam_write_cell(host, cid, _tam_cell())
            _tam_write_cell(host, bogus, _tam_cell(seed = 99))
            result = merge_hosts(["h" => host]; expected_ids = [cid])
            Test.@test result.unexpected_ids == [bogus]
        end
    end

    Test.@testset "cells_root accepts a results root or a cells/ dir" begin
        mktempdir() do dir
            root = joinpath(dir, "run")
            _tam_write_cell(root, "Lambda__illumina__30x__seed42__qualmer", _tam_cell())
            Test.@test cells_root(root) == joinpath(root, "cells")
            Test.@test cells_root(joinpath(root, "cells")) == joinpath(root, "cells")
            # Both spellings discover the same cell.
            Test.@test length(first(discover_cells(cells_root(root)))) == 1
            Test.@test length(
                first(discover_cells(cells_root(joinpath(root, "cells"))))) == 1
        end
    end

    Test.@testset "output tables carry the harness schema plus provenance" begin
        mktempdir() do dir
            a = joinpath(dir, "lovelace")
            b = joinpath(dir, "lrc")
            # Lovelace: QUAST failing (largest_contig 0 despite contigs).
            _tam_write_cell(a, "T4__illumina__30x__seed42__qualmer",
                _tam_cell(organism = "T4", n_contigs = 12, largest_contig = 0,
                    NGA50 = 0.0, genome_fraction = 0.0))
            # Lawrencium: QUAST working.
            _tam_write_cell(b, "phi29__illumina__30x__seed42__qualmer",
                _tam_cell(organism = "phi29", largest_contig = 19_000))
            ids = ["T4__illumina__30x__seed42__qualmer",
                "phi29__illumina__30x__seed42__qualmer"]
            result = merge_hosts(["lovelace" => a, "lrc" => b]; expected_ids = ids)

            results_df, prov_df = merged_tables(result)
            # Drop-in compatible with the harness's own aggregate.
            Test.@test DataFrames.names(results_df) == TRACK_A_ROW_KEYS
            Test.@test DataFrames.nrow(results_df) == 2
            # _tam_cell omits both provenance keys, i.e. the shape of all 432
            # checkpoints on disk. They must merge to the HARNESS sentinels, not
            # `missing`: with `missing` the merged table carried the column names and
            # no usable values, and the operation the harness docstring prescribes —
            # "always filter on peak_rss_method before aggregating" — threw
            # ArgumentError on Missing. Column names alone were not the fix.
            Test.@test all(results_df.peak_rss_method .== "unknown")
            Test.@test all(results_df.rss_baseline_bytes .== -1)
            # The prescribed filter must actually run on the merged table.
            Test.@test DataFrames.nrow(
                results_df[results_df.peak_rss_method .== "unknown", :]) == 2
            Test.@test "host" in DataFrames.names(prov_df)
            Test.@test "quast_evidence" in DataFrames.names(prov_df)
            Test.@test "source_digest" in DataFrames.names(prov_df)
            Test.@test Set(prov_df.host) == Set(["lovelace", "lrc"])

            # This is the item-1 <-> item-3 bridge: the merged matrix mixes metric
            # definitions BY HOST, and the guard must see it.
            breakdown = evidence_by_host(result)
            Test.@test breakdown["lovelace"]["unknown:quast-unscored"] == 1
            Test.@test breakdown["lrc"]["quast:scored"] == 1
            _tam_throws_with(
                () -> assert_single_metric_definition(prov_df;
                    columns = (:quast_evidence,), context = "merged Track A matrix"),
                ["quast_evidence", "quast:scored", "unknown:quast-unscored"])
        end
    end

    Test.@testset "report names collisions, missing cells, and the host breakdown" begin
        mktempdir() do dir
            a = joinpath(dir, "hostA")
            b = joinpath(dir, "hostB")
            cid = "Lambda__illumina__30x__seed42__qualmer"
            missing_id = "Lambda__illumina__30x__seed42__kmer"
            _tam_write_cell(a, cid, _tam_cell(NGA50 = 1316.0))
            _tam_write_cell(b, cid, _tam_cell(NGA50 = 1402.0))
            result = merge_hosts(["hostA" => a, "hostB" => b];
                expected_ids = [cid, missing_id])
            out = joinpath(dir, "merged")
            mkpath(out)
            path = write_merge_report(joinpath(out, "merge_report.md"), result;
                output_dir = out)
            text = read(path, String)
            Test.@test occursin("CONTENT COLLISIONS: 1", text)
            Test.@test occursin("Matrix complete and collision-free: NO", text)
            Test.@test occursin(cid, text)
            Test.@test occursin("NGA50: 1316 vs 1402", text)
            Test.@test occursin(missing_id, text)
            Test.@test occursin("quast_evidence x host", text)
        end
    end

    Test.@testset "materializing cells is idempotent; disagreement is a collision" begin
        mktempdir() do dir
            host = joinpath(dir, "lovelace")
            cid = "Lambda__illumina__30x__seed42__qualmer"
            _tam_write_cell(host, cid, _tam_cell())
            out = joinpath(dir, "merged")
            result = merge_hosts(["lovelace" => host]; expected_ids = [cid])

            copied, skipped, conflicts = copy_merged_cells(result, out)
            Test.@test copied == [cid]
            Test.@test isempty(skipped)
            Test.@test isempty(conflicts)
            Test.@test isfile(joinpath(out, "cells", cid, "cell_result.json"))

            # Re-run: nothing is rewritten, nothing is reported wrong. Both source
            # runs are resumable, so re-merging after new cells appear must be safe.
            copied2, skipped2, conflicts2 = copy_merged_cells(result, out)
            Test.@test isempty(copied2)
            Test.@test skipped2 == [cid]
            Test.@test isempty(conflicts2)

            # A previously merged cell that DISAGREES with its source is the same
            # class of anomaly as a cross-host collision — never a silent overwrite.
            _tam_write_cell(host, cid, _tam_cell(NGA50 = 9999.0))
            result2 = merge_hosts(["lovelace" => host]; expected_ids = [cid])
            copied3, skipped3, conflicts3 = copy_merged_cells(result2, out)
            Test.@test isempty(copied3)
            Test.@test isempty(skipped3)
            Test.@test length(conflicts3) == 1
            Test.@test conflicts3[1].cell_id == cid
            Test.@test "<merged-output>" in conflicts3[1].hosts
            # The destination on disk is untouched.
            kept = JSON.parsefile(joinpath(out, "cells", cid, "cell_result.json"))
            Test.@test kept["NGA50"] == 1316.0
        end
    end

    Test.@testset "--host spec parsing" begin
        Test.@test _parse_host_spec("lrc=/scratch/track_a") == ("lrc" => "/scratch/track_a")
        # A path containing '=' still parses: only the FIRST '=' separates.
        Test.@test _parse_host_spec("h=/a=b") == ("h" => "/a=b")
        _tam_throws_with(() -> _parse_host_spec("no-equals-sign"),
            ["--host expects LABEL=PATH"])
        _tam_throws_with(() -> _parse_host_spec("=/path"), ["label is empty"])
        _tam_throws_with(() -> _parse_host_spec("label="), ["path is empty"])
    end

    Test.@testset "C10: a nonexistent source is FATAL, not silently empty" begin
        mktempdir() do dir
            cid = "Lambda__illumina__30x__seed42__qualmer"
            result = merge_hosts(["ghost" => joinpath(dir, "does-not-exist")];
                expected_ids = [cid])
            Test.@test isempty(result.merged)
            Test.@test result.per_host[1].discovered == 0
            Test.@test result.per_host[1].exists == false
            Test.@test result.missing_ids == [cid]
            # A mistyped path or an unmounted filesystem must not read as "this
            # host contributed nothing", which silently shrinks the matrix.
            Test.@test result.unreachable_sources == ["ghost"]
            code, problems = merge_exit_status(result)
            Test.@test code == 1
            Test.@test :unreachable_source in [p.kind for p in problems]
            # NOT waivable: --allow-incomplete says the matrix is partial, not that
            # a requested source may vanish.
            code2,
            problems2 = merge_exit_status(result;
                allow_incomplete = true, allow_collisions = true)
            Test.@test code2 == 1
            Test.@test any(p -> p.kind == :unreachable_source && p.fatal, problems2)
            # Empty tables still have the right schema, so a downstream reader
            # fails on emptiness rather than on a missing column.
            results_df, prov_df = merged_tables(result)
            Test.@test DataFrames.names(results_df) == TRACK_A_ROW_KEYS
            Test.@test "quast_evidence" in DataFrames.names(prov_df)
            Test.@test DataFrames.nrow(results_df) == 0
        end
    end

    Test.@testset "C7: an unreadable field is malformed, never benign" begin
        # `_as_number -> 0.0` used to route every read failure into a benign class:
        # an unreadable n_contigs became "n/a:empty-assembly" and an unreadable
        # largest_contig became "unknown:quast-unscored", so schema drift and
        # truncated checkpoints passed at exit 0.
        truncated = _tam_cell()
        delete!(truncated, "n_contigs")
        Test.@test quast_evidence(truncated) == "malformed:unreadable(n_contigs)"
        Test.@test is_malformed_evidence(quast_evidence(truncated))

        drifted = _tam_cell()
        drifted["largest_contig"] = "n/a"          # schema change to a string
        Test.@test quast_evidence(drifted) == "malformed:unreadable(largest_contig)"

        nostatus = _tam_cell()
        delete!(nostatus, "status")
        Test.@test quast_evidence(nostatus) == "malformed:missing-field(status)"

        # A genuine empty assembly is still benign — the classes stay distinct.
        Test.@test !is_malformed_evidence(
            quast_evidence(_tam_cell(n_contigs = 0, status = "empty_assembly")))

        # And a malformed cell fails the merge rather than passing quietly.
        mktempdir() do dir
            host = joinpath(dir, "h")
            cid = "Lambda__illumina__30x__seed42__qualmer"
            _tam_write_cell(host, cid, truncated)
            result = merge_hosts(["h" => host]; expected_ids = [cid])
            Test.@test result.malformed_cells == [cid]
            code, problems = merge_exit_status(result)
            Test.@test code == 1
            Test.@test :malformed_cells in [p.kind for p in problems]
        end
    end

    Test.@testset "C10: tables are written only behind the gate, atomically" begin
        mktempdir() do dir
            out = joinpath(dir, "merged")
            mkpath(out)
            good = joinpath(out, "track_a_results.tsv")
            write(good, "organism\taccession\n" * repeat("Lambda\tNC_001416\n", 6))
            before = read(good, String)

            host = joinpath(dir, "partial")
            ids = expected_cell_ids(organisms = ("Lambda",),
                technologies = ("illumina",), coverages = (10, 30), seeds = (42,))
            _tam_write_cell(host, ids[1], _tam_cell())
            result = merge_hosts(["partial" => host]; expected_ids = ids)
            code, _ = merge_exit_status(result)
            Test.@test code == 1
            Test.@test read(good, String) == before

            # ATOMICITY, asserted directly rather than implied: a write that throws
            # part-way must leave the destination untouched and no temp behind.
            Test.@test_throws ErrorException _atomic_write(good) do io
                write(io, "partial junk")
                error("simulated interruption")
            end
            Test.@test read(good, String) == before
            Test.@test isempty(filter(f -> startswith(f, "."), readdir(out)))

            # And a successful atomic write does replace the contents.
            _atomic_write(good) do io
                write(io, "replaced\n")
            end
            Test.@test read(good, String) == "replaced\n"
        end
    end

    Test.@testset "CR: a path that exists but is not a cells tree is unreachable" begin
        mktempdir() do dir
            # A typo landing on some real-but-wrong directory used to report
            # exists=true with zero cells and escape the guard.
            notatree = joinpath(dir, "some-other-dir")
            mkpath(joinpath(notatree, "unrelated"))
            cid = "Lambda__illumina__30x__seed42__qualmer"
            result = merge_hosts(["typo" => notatree]; expected_ids = [cid])
            Test.@test result.unreachable_sources == ["typo"]
            code, problems = merge_exit_status(result)
            Test.@test code == 1
            Test.@test any(p -> occursin("not a checkpoint tree", p.message), problems)
        end
    end

    Test.@testset "CR: malformed cells are counted in the report's evidence table" begin
        mktempdir() do dir
            host = joinpath(dir, "h")
            broken = _tam_cell()
            delete!(broken, "largest_contig")
            cid = "Lambda__illumina__30x__seed42__qualmer"
            _tam_write_cell(host, cid, broken)
            result = merge_hosts(["h" => host]; expected_ids = [cid])
            path = write_merge_report(joinpath(dir, "r.md"), result; output_dir = dir)
            text = read(path, String)
            # The class carries a suffix, so a bare-name lookup rendered 0 forever.
            row = first(filter(l -> startswith(l, "| h |"), split(text, "\n")))
            Test.@test occursin("1", row)
            Test.@test !occursin("| h | 0 | 0 | 0 | 0 | 0 |", row)
        end
    end

    Test.@testset "CR: nested values canonicalize stably" begin
        # Falling through to `string(v)` rendered a nested object in Dict iteration
        # order, so two hosts with EQUAL content could digest differently.
        a = _tam_cell();
        a["extra"] = Dict("z" => 1, "a" => 2, "m" => [1, 2, 3])
        b = _tam_cell();
        b["extra"] = Dict("a" => 2, "m" => [1, 2, 3], "z" => 1)
        Test.@test cell_digest(a) == cell_digest(b)
        Test.@test isempty(differing_fields(a, b))
        c = _tam_cell();
        c["extra"] = Dict("z" => 9, "a" => 2, "m" => [1, 2, 3])
        Test.@test cell_digest(a) != cell_digest(c)
    end

    Test.@testset "CR: a non-numeric shard flag is usage error, not a stacktrace" begin
        saved = copy(ARGS)
        try
            empty!(ARGS);
            append!(ARGS, ["--coverages", "10,notanumber"])
            Test.@test _parse_int_flag("--coverages", TRACK_A_COVERAGES) === nothing
            empty!(ARGS);
            append!(ARGS, ["--coverages", "10,30"])
            Test.@test _parse_int_flag("--coverages", TRACK_A_COVERAGES) == [10, 30]
            empty!(ARGS)
            Test.@test _parse_int_flag("--seeds", TRACK_A_SEEDS) == TRACK_A_SEEDS
        finally
            empty!(ARGS);
            append!(ARGS, saved)
        end
    end
end
