# Gated real-genome smoke: aggregate qualmer storage must preserve corrector
# QUALITY vs :full at real scale (td-n8ax / td-oagi). The 140bp unit toy is BELOW
# the corrector's regime and cannot answer this (it shows a spurious gap); this
# runs on phix174 where the corrector operates as designed.
#
# Gated behind MYCELIA_RUN_EXTERNAL (needs ART for read simulation, conda-
# bootstrapped, + a network fetch of the reference) so the default unit suite stays
# fast/offline. Standing regression guard: a future qualmer-storage change that
# silently degraded correction quality would keep testsets A-H green but fail here.
#
# Run: MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/4_assembly/qualmer_profile_quality_parity_test.jl")'
# Empirical baseline (2026-07-16): phix @30x k19 -> :full and :ultralight_quality
# both GF 0.987; lambda @30x -> both 0.9981. Parity holds; aggregate is faster.

import Test, Mycelia, FASTX, BioSequences

const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") in ["1", "true", "yes"]

if !RUN_EXTERNAL
    @info "Skipping qualmer profile quality parity smoke (set MYCELIA_RUN_EXTERNAL=true to run)"
else
    Test.@testset "Aggregate qualmer storage preserves corrector quality (phix)" begin
        k = 19
        cov = 10  # low coverage keeps the smoke fast; still in-regime for phix
        seed = 42

        tmp = mktempdir(prefix = "qual_parity_")
        refdir = joinpath(tmp, "ref")
        mkpath(refdir)
        ref = Mycelia.download_genome_by_accession(
            accession = "NC_001422.1", outdir = refdir, compressed = false)
        refstr = uppercase(String(FASTX.sequence(first(collect(Mycelia.open_fastx(ref))))))
        Test.@test length(refstr) > 5000  # phix ~5386 bp

        sim = Mycelia.simulate_illumina_reads(fasta = ref, coverage = cov, rndSeed = seed,
            read_length = 150, paired = true, quiet = true, outbase = joinpath(tmp, "r"))
        reads = FASTX.FASTQ.Record[]
        for p in filter(!isnothing, [sim.forward_reads, sim.reverse_reads])
            for r in Mycelia.open_fastx(p)
                push!(reads, r)
            end
        end
        Test.@test !isempty(reads)

        rc(s) = String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(String(s))))
        function ckmers(seqs)
            s = Set{String}()
            for x in String.(seqs), i in 1:(length(x) - k + 1)

                sub = x[i:(i + k - 1)]
                occursin(r"^[ACGT]+$", sub) || continue
                push!(s, min(sub, rc(sub)))
            end
            return s
        end
        refk = ckmers([refstr])
        gf(cs) = isempty(cs) ? 0.0 : length(intersect(refk, ckmers(cs))) / length(refk)

        function run_gf(profile)
            cfg = Mycelia.Rhizomorph.AssemblyConfig(; k = k, corrector = :iterative,
                qualmer_memory_profile = profile, qualmer_prefilter_min_count = 2,
                verbose = false)
            gf(Mycelia.Rhizomorph.assemble_genome(reads, cfg).contigs)
        end

        gf_full = run_gf(:full)
        gf_ultra = run_gf(:ultralight_quality)
        @info "phix genome fraction" gf_full gf_ultra
        # Both must recover the genome, and the aggregate profile must not materially
        # underperform :full (the whole point: memory-bounded AND quality-preserving).
        Test.@test gf_full >= 0.9
        Test.@test gf_ultra >= gf_full - 0.01
    end
end

println("✓ qualmer profile quality parity smoke completed")
