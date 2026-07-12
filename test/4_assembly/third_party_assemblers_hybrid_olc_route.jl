# td-yymj: hybrid-OLC route (a) END-TO-END (external short-read assemblers).
#
# Exercises assemble_genome(reads; corrector=:iterative, layout=:olc) — Stage-1
# correction of simulated Illumina reads handed to an external short-read layout
# assembler (megahit / metaspades) — and checks the wrapped AssemblyResult. Needs
# the bioconda assembler tools + ART (simulate_illumina_reads), so it is gated:
# the "third_party_assemblers" filename skips it in default CI, and the body also
# guards on MYCELIA_RUN_EXTERNAL so a direct include is a no-op without the tools.
#
# Run directly (from the Mycelia base directory):
#   MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/4_assembly/third_party_assemblers_hybrid_olc_route.jl")'

import Test
import Mycelia
import StableRNGs
import BioSequences
import FASTX

const R = Mycelia.Rhizomorph

Test.@testset "hybrid-OLC route (a) end-to-end (td-yymj)" begin
    run_external = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true" ||
                   lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"

    if !run_external
        @info "hybrid-OLC external run disabled; set MYCELIA_RUN_EXTERNAL=true to enable."
    else
        threads = clamp(
            something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2),
            1, 4)
        mktempdir() do dir
            # Small synthetic reference + single-end Illumina reads.
            rng = StableRNGs.StableRNG(2026)
            genome = BioSequences.randdnaseq(rng, 5_000)
            ref_fasta = joinpath(dir, "ref.fasta")
            Mycelia.write_fasta(outfile = ref_fasta,
                records = [FASTX.FASTA.Record("ref", genome)])
            sim = Mycelia.simulate_illumina_reads(
                fasta = ref_fasta, coverage = 25, read_length = 150,
                paired = false, seqSys = "HS25", rndSeed = 2026, quiet = true)
            reads = collect(Mycelia.open_fastx(sim.forward_reads))
            Test.@test !isempty(reads)

            # The corrector emits a SINGLE (single-end) corrected read set, so the
            # end-to-end assertion uses megahit, which handles single-end input
            # robustly. :metaspades is also wired (route (a) short-read arm) but
            # metaSPAdes targets PAIRED metagenomic data and is not representative
            # with single-end corrected reads, so it is exercised at the config /
            # adapter level (hybrid_olc_config_test.jl) rather than end-to-end here.
            Test.@testset "layout=:olc olc_tool=:megahit" begin
                config = R.AssemblyConfig(; k = 13, corrector = :iterative,
                    strategy = :scalable, sequencing_tech = :illumina,
                    layout = :olc, olc_tool = :megahit,
                    olc_options = (; k_list = "21", threads = threads))
                asm = R.assemble_genome(reads, config)
                Test.@test asm isa R.AssemblyResult
                # External contig-only assembly: no graph, not GFA-compatible.
                Test.@test asm.gfa_compatible == false
                Test.@test asm.graph === nothing
                # Produced at least one contig, with stamped provenance.
                Test.@test !isempty(asm.contigs)
                Test.@test length(asm.contigs) == length(asm.contig_names)
                Test.@test asm.assembly_stats["method"] == "HybridOLC"
                Test.@test asm.assembly_stats["layout"] == "olc"
                Test.@test asm.assembly_stats["olc_tool"] == "megahit"
                Test.@test asm.assembly_stats["corrected_read_count"] > 0
                # Ephemeral contract (output_dir unset): the corrected FASTQ the
                # route persisted (stamped into stats) must be CLEANED after the run.
                Test.@test !isfile(asm.assembly_stats["corrected_fastq"])
            end

            Test.@testset "layout=:olc output_dir persists artifacts" begin
                # Caller-owned output_dir: the corrected FASTQ + the assembler output
                # dir must SURVIVE the run for the Stage-2 handoff (not cleaned).
                outdir = mktempdir()
                config = R.AssemblyConfig(; k = 13, corrector = :iterative,
                    strategy = :scalable, sequencing_tech = :illumina,
                    layout = :olc, olc_tool = :megahit, output_dir = outdir,
                    olc_options = (; k_list = "21", threads = threads))
                asm = R.assemble_genome(reads, config)
                Test.@test asm isa R.AssemblyResult
                Test.@test !isempty(asm.contigs)
                # Corrected FASTQ persisted at the fixed path, left in place.
                Test.@test isfile(joinpath(outdir, "corrected.fastq"))
                Test.@test asm.assembly_stats["corrected_fastq"] ==
                           joinpath(outdir, "corrected.fastq")
                # Assembler output dir survives too.
                Test.@test isdir(joinpath(outdir, "olc_megahit"))
            end
        end
    end
end
