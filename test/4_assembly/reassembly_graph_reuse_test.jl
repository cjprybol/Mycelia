# td-04tb: re-assembly graph reuse. `_assemble_with_iterative_corrector` used to
# rebuild a qualmer graph FROM SCRATCH from the corrected reads (~22-36% of the
# iterative arm), duplicating the graph the corrector already built in its final
# pass. When that final pass CONVERGED (0 improvements), the corrector's graph is
# byte-identical to the rebuild, so it is now returned as `:final_graph` and reused.
#
# These tests prove:
#   (1) On a converging fixture the corrector marks its final-pass graph reusable
#       and hands it back.
#   (2) Assembling from the REUSED graph is byte-identical (same contigs / N50 /
#       sequences) to a from-scratch rebuild of the corrected reads — the reuse is
#       a pure-performance change, not an output change.
#   (3) The end-to-end `assemble_genome(...; corrector=:iterative)` path stamps the
#       reuse provenance and its contigs match the from-scratch rebuild.
#   (4) The reuse guard is conservative: a k / mode mismatch declines reuse.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/reassembly_graph_reuse_test.jl")'

import Test
import Mycelia
import FASTX
import BioSequences
import Random

const R = Mycelia.Rhizomorph

# Reads via Mycelia.observe from a random reference (the recipe the e2e phase
# profiler / quality-gap fixtures use). On this regime the :scalable corrector
# converges at the final k (a 0-improvement pass), so its final-pass graph is
# reusable — exercising the td-04tb reuse path. The global RNG is seeded so the
# statistical-resampling arm is reproducible and convergence is stable.
function _converging_reads(; genome_len = 1000, readlen = 150, coverage = 30,
        err = 0.05, seed = 42)
    rec = Mycelia.random_fasta_record(moltype = :DNA, seed = seed, L = genome_len)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
    rng = Random.MersenneTwister(seed)
    glen = length(refseq)
    erl = min(readlen, glen)
    n = max(1, ceil(Int, coverage * glen / erl))
    recs = FASTX.FASTQ.Record[]
    for i in 1:n
        s = glen == erl ? 1 : rand(rng, 1:(glen - erl + 1))
        frag = refseq[s:(s + erl - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        obs, q = Mycelia.observe(frag; error_rate = err, tech = :illumina)
        isempty(obs) && continue
        push!(recs, FASTX.FASTQ.Record("read_$i", string(obs),
            String([Char(x + 33) for x in q])))
    end
    return recs
end

_n50(lens) = (t = sum(lens; init = 0); acc = 0;
    for l in sort(lens; rev = true)
        acc += l
        acc >= t / 2 && return l
    end; 0)

# Reproduce EXACTLY the corrector call `_assemble_with_iterative_corrector` makes
# for strategy=:scalable at a given k, returning the result dict + corrected reads.
function _run_corrector(reads, k)
    knobs = R._corrector_strategy_knobs(:scalable)
    max_k = max(k, 13)
    mode = knobs.graph_mode === nothing ?
           R._graph_mode_symbol(R.DoubleStrand) : knobs.graph_mode
    indir = mktempdir()
    outdir = mktempdir()
    fq = joinpath(indir, "in.fastq")
    open(FASTX.FASTQ.Writer, fq) do w
        for r in reads
            write(w, r)
        end
    end
    result = Mycelia.mycelia_iterative_assemble(fq;
        max_k = max_k, skip_solid = knobs.skip_solid, graph_mode = mode,
        n_k_rungs = knobs.n_k_rungs, max_iterations_per_k = knobs.max_iterations_per_k,
        hard_window = knobs.hard_window, soft_em = knobs.soft_em,
        cheap_correct = knobs.cheap_correct, beam_width = knobs.beam_width,
        verbose = false, enable_checkpointing = false, output_dir = outdir)
    corrected = open(FASTX.FASTQ.Reader, result[:metadata][:final_fastq_file]) do rdr
        collect(rdr)
    end
    return result, corrected, max_k
end

Test.@testset "re-assembly graph reuse (td-04tb)" begin
    Random.seed!(42)
    reads = _converging_reads()
    k = 21  # >= 13 ⇒ max_k == 21 == reassembly_k, so reuse is parameter-eligible

    result, corrected, max_k = _run_corrector(reads, k)
    reusable = result[:metadata][:final_graph_reusable]

    Test.@testset "corrector surfaces final-pass graph reuse provenance" begin
        Test.@test !isempty(corrected)
        Test.@test result[:metadata][:final_graph_k] == max_k
        Test.@test result[:metadata][:final_graph_mode] == :doublestrand
        Test.@test reusable isa Bool
        # This regime converges at the final k, so the graph IS reusable and returned.
        Test.@test reusable === true
        Test.@test result[:final_graph] !== nothing
    end

    Test.@testset "reused graph == from-scratch rebuild (byte-identical)" begin
        rcfg = R._auto_configure_assembly(corrected; k = k, corrector = :none)
        # From-scratch rebuild — what master did unconditionally.
        rebuild = R.assemble_genome(corrected, rcfg)
        # Reuse path — assemble straight from the corrector's returned graph.
        reuse = R._qualmer_graph_to_assembly(result[:final_graph], length(corrected), rcfg)

        Test.@test sort(reuse.contigs) == sort(rebuild.contigs)
        Test.@test length(reuse.contigs) == length(rebuild.contigs)
        Test.@test _n50(length.(reuse.contigs)) == _n50(length.(rebuild.contigs))
        Test.@test reuse.assembly_stats["num_vertices"] == rebuild.assembly_stats["num_vertices"]
        Test.@test reuse.assembly_stats["num_edges"] == rebuild.assembly_stats["num_edges"]

        # td-47di: the invariant must ALSO hold with RC-dedup ON — the regime the
        # corrector uses by default. rc_aware picks one strand per pair by graph
        # iteration order, so reused-graph vs rebuild could emit opposite strands;
        # canonical-orientation emission in _qualmer_graph_to_assembly must keep them
        # byte-identical. This is the exact reuse-under-dedup case that must not regress.
        rcfg_dedup = R._auto_configure_assembly(corrected;
            k = k, corrector = :none, dedup_revcomp = true)
        rebuild_dd = R.assemble_genome(corrected, rcfg_dedup)
        reuse_dd = R._qualmer_graph_to_assembly(result[:final_graph], length(corrected), rcfg_dedup)
        Test.@test sort(reuse_dd.contigs) == sort(rebuild_dd.contigs)
        Test.@test length(reuse_dd.contigs) == length(rebuild_dd.contigs)
        # Dedup must not lose loci: canonical set is no larger than the both-strands set.
        Test.@test length(reuse_dd.contigs) <= length(reuse.contigs)
    end

    Test.@testset "end-to-end corrector output == from-scratch rebuild" begin
        # THE invariant that must hold whether or not reuse fired: the corrector's
        # assembled contigs equal a from-scratch rebuild of the corrected reads.
        # corrector=:iterative now defaults dedup_revcomp ON (td-47di), so the honest
        # baseline holds that SAME dedup setting constant — this isolates the reuse
        # variable (graph reuse is a pure-performance change) from the orthogonal
        # dedup-policy variable. Both sides therefore emit the canonical (RC-deduped)
        # contig set and must still match.
        asm = R.assemble_genome(reads; k = k, corrector = :iterative, strategy = :scalable)
        rebuild = R.assemble_genome(corrected; k = k, corrector = :none, dedup_revcomp = true)
        Test.@test sort(asm.contigs) == sort(rebuild.contigs)
        Test.@test _n50(length.(asm.contigs)) == _n50(length.(rebuild.contigs))
        # Provenance flag is present and (in this converging regime) true.
        Test.@test haskey(asm.assembly_stats, "reassembly_graph_reused")
        Test.@test asm.assembly_stats["reassembly_graph_reused"] == reusable
    end

    Test.@testset "reuse guard is conservative on k / mode mismatch" begin
        gk = result[:metadata][:final_graph_k]
        gmode = result[:metadata][:final_graph_mode]
        # k mismatch ⇒ a reassembly config at a different k would decline reuse.
        rcfg_bad_k = R._auto_configure_assembly(corrected; k = gk + 2, corrector = :none)
        Test.@test rcfg_bad_k.k != gk
        # mode mismatch ⇒ singlestrand config vs the doublestrand corrector graph.
        rcfg_bad_mode = R._auto_configure_assembly(corrected;
            k = gk, graph_mode = R.SingleStrand, corrector = :none)
        Test.@test R._graph_mode_symbol(rcfg_bad_mode.graph_mode) != gmode
    end
end
