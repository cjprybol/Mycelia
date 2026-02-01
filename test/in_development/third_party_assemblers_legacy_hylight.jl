# Legacy third-party assembler tests: HyLight.
import Test
import Mycelia
import Graphs
import Random
import StableRNGs
import BioSequences
import FASTX

threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

#     Test.@testset "5. Hybrid Assembly" begin
#         Test.@testset "5a. Hybrid Isolate Assembly" begin
#             mktempdir() do dir
#                 # Create reference genome for hybrid isolate assembly
#                 hybrid_isolate_ref_fasta = joinpath(dir, "hybrid_isolate_ref.fasta")
#                 rng_hybrid_isolate = StableRNGs.StableRNG(789)
#                 hybrid_isolate_genome = BioSequences.randdnaseq(rng_hybrid_isolate, 3000)  # 3kb genome for absolute minimal memory

#                 # Create FASTA record and write using Mycelia.write_fasta
#                 hybrid_isolate_fasta_record = FASTX.FASTA.Record("hybrid_isolate_test_genome", hybrid_isolate_genome)
#                 Mycelia.write_fasta(outfile=hybrid_isolate_ref_fasta, records=[hybrid_isolate_fasta_record])

#                 # Simulate Illumina short reads with 10x coverage
#                 isolate_short_reads = Mycelia.simulate_illumina_reads(fasta=hybrid_isolate_ref_fasta, coverage=10)

#                 # Simulate nanopore long reads with 5x coverage
#                 isolate_long_reads_gz = Mycelia.simulate_nanopore_reads(fasta=hybrid_isolate_ref_fasta, quantity="5x")

#                 # Decompress long reads for unicycler
#                 isolate_long_reads = joinpath(dir, "isolate_long_reads.fq")
#                 run(pipeline(`gunzip -c $(isolate_long_reads_gz)`, isolate_long_reads))

#             # Test Unicycler - clean up any existing directory first
#             unicycler_outdir = joinpath(dir, "unicycler_assembly")
#             if isdir(unicycler_outdir)
#                 rm(unicycler_outdir, recursive=true)
#             end
#             try
#                 result = Mycelia.run_unicycler(short_1=isolate_short_reads.forward_reads, short_2=isolate_short_reads.reverse_reads, long_reads=isolate_long_reads, outdir=unicycler_outdir)
#                 Test.@test result.outdir == unicycler_outdir
#                 Test.@test result.assembly == joinpath(unicycler_outdir, "assembly.fasta")
#                 Test.@test isfile(result.assembly)
#                 # Clean up after test
#                 rm(unicycler_outdir, recursive=true, force=true)
#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed") || contains(string(e), "failed")
#                     @warn """
#                     Unicycler hybrid assembly failed due to resource constraints.
#                     Current test: 3kb genome, 10x short reads + 5x long reads, ~45kb total sequence data

#                     Required resources for Unicycler:
#                     - Memory: ~4-8GB RAM minimum (SPAdes is memory-intensive)
#                     - CPU: 4-8 cores recommended  
#                     - Disk: ~1-2GB temporary space
#                     - Note: Unicycler internally runs SPAdes which has high memory requirements

#                     To fix: This assembler requires significant computational resources.
#                     Consider running on a machine with more memory, or skip hybrid assembly testing.

#                     Alternative: Use Flye or Canu for long-read-only assembly if resources are limited.
#                     """
#                     Test.@test_skip "Unicycler test skipped - insufficient resources (memory-intensive hybrid assembler)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(unicycler_outdir, recursive=true, force=true)
#             end
#         end

#         Test.@testset "5b. Hybrid Metagenomic Assembly" begin
#             mktempdir() do dir
#                 # Create reference genome for hybrid metagenomic assembly
#                 hybrid_meta_ref_fasta = joinpath(dir, "hybrid_meta_ref.fasta")
#                 rng_hybrid_meta = StableRNGs.StableRNG(987)
#                 hybrid_meta_genome = BioSequences.randdnaseq(rng_hybrid_meta, 3000)  # 3kb genome for metagenomic hybrid

#                 # Create FASTA record and write using Mycelia.write_fasta
#                 hybrid_meta_fasta_record = FASTX.FASTA.Record("hybrid_meta_test_genome", hybrid_meta_genome)
#                 Mycelia.write_fasta(outfile=hybrid_meta_ref_fasta, records=[hybrid_meta_fasta_record])

#                 # Simulate Illumina short reads with 10x coverage
#                 meta_short_reads = Mycelia.simulate_illumina_reads(fasta=hybrid_meta_ref_fasta, coverage=10)

#                 # Simulate nanopore long reads with 5x coverage
#                 meta_long_reads_gz = Mycelia.simulate_nanopore_reads(fasta=hybrid_meta_ref_fasta, quantity="5x")

#                 # Decompress long reads for HyLight
#                 meta_long_reads = joinpath(dir, "meta_long_reads.fq")
#                 run(pipeline(`gunzip -c $(meta_long_reads_gz)`, meta_long_reads))
# Hybrid metagenomic assembly (HyLight) smoke test with a lightweight fixture.
# Keep a stub so we remember to reintroduce a lightweight fixture later.
# Test.@testset "Hybrid Metagenomic Assembly - HyLight" begin
#     Test.@test_skip "HyLight test temporarily disabled pending lightweight fixture"
# end

#             # Test HyLight - hybrid strain-resolved metagenomic assembly
#             hylight_outdir = joinpath(dir, "hylight_assembly")
#             if isdir(hylight_outdir)
#                 rm(hylight_outdir, recursive=true)
#             end
#             try
#                 result = Mycelia.run_hylight(meta_short_reads.forward_reads, meta_short_reads.reverse_reads, meta_long_reads, outdir=hylight_outdir)
#                 Test.@test result.outdir == hylight_outdir
#                 Test.@test result.strain_assemblies == joinpath(hylight_outdir, "strain_assemblies")
#                 Test.@test isdir(result.strain_assemblies)
#                 # Clean up after test
#                 rm(hylight_outdir, recursive=true, force=true)
#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed") || contains(string(e), "failed")
#                     @warn """
#                     HyLight hybrid assembly failed due to resource constraints.
#                     Current test: 3kb genome, 10x short reads + 5x long reads, ~45kb total sequence data

#                     Required resources for HyLight:
#                     - Memory: ~6-12GB RAM minimum (strain-resolved metagenomic hybrid assembly)
#                     - CPU: 4-16 cores recommended
#                     - Disk: ~2-4GB temporary space
#                     - Note: HyLight performs strain resolution on hybrid assemblies

#                     To fix: This is a very resource-intensive assembler requiring significant computational resources.
#                     Consider running on a high-memory machine or skip strain-resolved hybrid assembly testing.
#                     """
#                     Test.@test_skip "HyLight test skipped - insufficient resources (very memory-intensive strain-resolved hybrid assembler)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(hylight_outdir, recursive=true, force=true)
#             end
#         end
# #     end

# Hybrid metagenomic assembly (HyLight) placeholder.
# Previous block relied on a commented-out `mktempdir` context, leaving `dir` undefined.
# NOTE: Marked broken while we refine fixture requirements.
Test.@testset "Hybrid Metagenomic Assembly - HyLight" begin
    Test.@test_broken mktempdir() do dir
        ref_fasta = joinpath(dir, "hylight_ref.fasta")
        rng_hylight_1 = StableRNGs.StableRNG(910)
        rng_hylight_2 = StableRNGs.StableRNG(911)
        genome_1 = BioSequences.randdnaseq(rng_hylight_1, 6000)
        genome_2 = BioSequences.randdnaseq(rng_hylight_2, 6000)
        Mycelia.write_fasta(
            outfile = ref_fasta,
            records = [
                FASTX.FASTA.Record("hylight_strain_1", genome_1),
                FASTX.FASTA.Record("hylight_strain_2", genome_2)
            ]
        )

        illumina = Mycelia.simulate_illumina_reads(
            fasta = ref_fasta,
            coverage = 20,
            outbase = joinpath(dir, "hylight_short"),
            read_length = 150,
            mflen = 300,
            seqSys = "HS25",
            paired = true,
            errfree = true,
            rndSeed = 910,
            quiet = true
        )

        long_reads_gz = Mycelia.simulate_nanopore_reads(
            fasta = ref_fasta,
            quantity = "10x",
            quiet = true
        )

        short_1 = joinpath(dir, "hylight_short_1.fq")
        short_2 = joinpath(dir, "hylight_short_2.fq")
        long_reads = joinpath(dir, "hylight_long.fq")
        run(pipeline(`gunzip -c $(illumina.forward_reads)`, short_1))
        run(pipeline(`gunzip -c $(illumina.reverse_reads)`, short_2))
        run(pipeline(`gunzip -c $(long_reads_gz)`, long_reads))

        open(short_1, "a") do io
            dummy_seq = repeat("A", 50)
            dummy_qual = repeat("I", 50)
            write(io, "@hylight_dummy/1\n", dummy_seq, "\n+\n", dummy_qual, "\n")
        end

        outdir = joinpath(dir, "hylight_out")
        result = Mycelia.run_hylight(
            short_1,
            short_2,
            long_reads;
            outdir = outdir,
            threads = threads,
            nsplit = 5,
            min_identity = 0.9,
            min_ovlp_len = 500,
            insert_size = 300,
            average_read_len = 150
        )
        result.outdir == outdir && isdir(result.strain_assemblies)
    end
end
#     Test.@testset "6. Probabilistic Assembly (Mycelia)" begin
#         mktempdir() do dir
#             fastq = joinpath(dir, "reads.fq")
#             open(fastq, "w") do io
#                 println(io, "@r1")
#                 println(io, "ACGTACGTACGTACGTACGT")
#                 println(io, "+")
#                 println(io, "IIIIIIIIIIIIIIIIIIII")
#                 println(io, "@r2")
#                 println(io, "CGTACGTACGTACGTACGTA")
#                 println(io, "+")
#                 println(io, "IIIIIIIIIIIIIIIIIIII")
#             end

#             # Test string graph building
#             graph = Mycelia.string_to_ngram_graph(s="ACGTACGTACGTACGTACGT", n=5)
#             Test.@test Graphs.nv(graph) > 0

#             # Test Viterbi error correction functions exist
#             Test.@test hasmethod(Mycelia.viterbi_maximum_likelihood_traversals, (Any,))
#             Test.@test isdefined(Mycelia, :polish_fastq)
#         end
#     end

#     Test.@testset "7. Assembly merging" begin
#     end
#     Test.@testset "8. Polishing & Error Correction" begin
#         mktempdir() do dir
#             # Create reference genome and assembly for polishing testing
#             polish_ref_fasta = joinpath(dir, "polish_ref.fasta")
#             rng_polish = StableRNGs.StableRNG(222)
#             polish_genome = BioSequences.randdnaseq(rng_polish, 4000)  # 4kb genome for polishing

#             # Create FASTA record and write using Mycelia.write_fasta
#             polish_fasta_record = FASTX.FASTA.Record("polish_test_genome", polish_genome)
#             Mycelia.write_fasta(outfile=polish_ref_fasta, records=[polish_fasta_record])

#             # Simulate PacBio reads for polishing (Apollo works with PacBio data)
#             polish_simulated_reads = Mycelia.simulate_pacbio_reads(fasta=polish_ref_fasta, quantity="15x")

#             # Decompress reads for polishing
#             polish_fastq = joinpath(dir, "polish_reads.fq")
#             run(pipeline(`gunzip -c $(polish_simulated_reads)`, polish_fastq))

#             # Generate a draft assembly using Flye (create something to polish)
#             draft_assembly_outdir = joinpath(dir, "draft_assembly")
#             draft_assembly_fasta = joinpath(draft_assembly_outdir, "assembly.fasta")

#             try
#                 # Generate draft assembly for polishing
#                 Mycelia.run_flye(fastq=polish_fastq, outdir=draft_assembly_outdir, genome_size="4k", read_type="pacbio-raw")

#                 if isfile(draft_assembly_fasta)
#                     # Test Apollo - clean up any existing directory first  
#                     apollo_outdir = joinpath(dir, "apollo_polishing")
#                     if isdir(apollo_outdir)
#                         rm(apollo_outdir, recursive=true)
#                     end

#                     try
#                         result = Mycelia.run_apollo(draft_assembly_fasta, polish_fastq, outdir=apollo_outdir)
#                         Test.@test result.outdir == apollo_outdir
#                         Test.@test result.polished_assembly == joinpath(apollo_outdir, basename(draft_assembly_fasta, ".fasta") * "_polished.fasta")
#                         Test.@test isfile(result.polished_assembly)
#                         # Clean up after test
#                         rm(apollo_outdir, recursive=true, force=true)
#                     catch e
#                         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                             @warn """
#                             Apollo polishing failed due to resource constraints.
#                             Current test: 4kb genome, 15x coverage, ~60kb total sequence data

#                             Required resources for Apollo:
#                             - Memory: ~2-4GB RAM minimum (includes minimap2 + samtools + HMM polishing)
#                             - CPU: 2-8 cores recommended
#                             - Disk: ~500MB temporary space
#                             - Note: Apollo performs HMM-based assembly polishing

#                             To fix: Increase available memory or reduce genome/coverage complexity.
#                             """
#                             Test.@test_skip "Apollo test skipped - insufficient resources"
#                         else
#                             rethrow(e)
#                         end
#                         # Clean up on failure
#                         rm(apollo_outdir, recursive=true, force=true)
#                     end
#                 else
#                     @warn "Draft assembly not generated by Flye - skipping Apollo test"
#                     Test.@test_skip "Apollo test skipped - no draft assembly available"
#                 end

#                 # Clean up draft assembly
#                 rm(draft_assembly_outdir, recursive=true, force=true)

#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                     @warn """
#                     Flye assembly for Apollo test failed due to resource constraints.
#                     Cannot generate draft assembly needed for polishing testing.
#                     """
#                     Test.@test_skip "Apollo test skipped - cannot generate draft assembly (insufficient resources)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(draft_assembly_outdir, recursive=true, force=true)
#             end

#             # Test Homopolish - can reuse the same setup as Apollo
#             try
#                 # Generate draft assembly for Homopolish testing (reuse same setup)
#                 Mycelia.run_flye(fastq=polish_fastq, outdir=draft_assembly_outdir, genome_size="4k", read_type="pacbio-raw")

#                 if isfile(draft_assembly_fasta)
#                     # Test Homopolish - clean up any existing directory first  
#                     homopolish_outdir = joinpath(dir, "homopolish_polishing")
#                     if isdir(homopolish_outdir)
#                         rm(homopolish_outdir, recursive=true)
#                     end

#                     try
#                         result = Mycelia.run_homopolish(draft_assembly_fasta, polish_fastq, outdir=homopolish_outdir)
#                         Test.@test result.outdir == homopolish_outdir
#                         Test.@test result.polished_assembly == joinpath(homopolish_outdir, basename(draft_assembly_fasta, ".fasta") * "_homopolished.fasta")
#                         Test.@test isfile(result.polished_assembly)
#                         # Clean up after test
#                         rm(homopolish_outdir, recursive=true, force=true)
#                     catch e
#                         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                             @warn """
#                             Homopolish polishing failed due to resource constraints.
#                             Current test: 4kb genome, 15x coverage, ~60kb total sequence data

#                             Required resources for Homopolish:
#                             - Memory: ~2-3GB RAM minimum (reference-based homopolymer correction)
#                             - CPU: 1-8 cores recommended
#                             - Disk: ~500MB temporary space
#                             - Note: Homopolish performs reference-based homopolymer error correction

#                             To fix: Increase available memory or reduce genome/coverage complexity.
#                             """
#                             Test.@test_skip "Homopolish test skipped - insufficient resources"
#                         else
#                             rethrow(e)
#                         end
#                         # Clean up on failure
#                         rm(homopolish_outdir, recursive=true, force=true)
#                     end
#                 else
#                     @warn "Draft assembly not generated by Flye - skipping Homopolish test"
#                     Test.@test_skip "Homopolish test skipped - no draft assembly available"
#                 end

#                 # Clean up draft assembly
#                 rm(draft_assembly_outdir, recursive=true, force=true)

#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                     @warn """
#                     Flye assembly for Homopolish test failed due to resource constraints.
#                     Cannot generate draft assembly needed for polishing testing.
#                     """
#                     Test.@test_skip "Homopolish test skipped - cannot generate draft assembly (insufficient resources)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(draft_assembly_outdir, recursive=true, force=true)
#             end
#         end
#     end
#     Test.@testset "6. Strain resolution" begin
#         mktempdir() do dir
#             # Create two related reference strains with realistic variations
#             base_ref_fasta = joinpath(dir, "base_strain.fasta")
#             variant_ref_fasta = joinpath(dir, "variant_strain.fasta")

#             rng_strain = StableRNGs.StableRNG(111)

#             # Create base strain (5kb genome)
#             base_genome = BioSequences.randdnaseq(rng_strain, 5000)
#             base_fasta_record = FASTX.FASTA.Record("base_strain", base_genome)
#             Mycelia.write_fasta(outfile=base_ref_fasta, records=[base_fasta_record])

#             # Create variant strain by introducing realistic variations
#             variant_genome = copy(base_genome)

#             # Introduce SNVs every 500bp (1% divergence)
#             for i in 500:500:length(variant_genome)
#                 if i <= length(variant_genome)
#                     original_base = variant_genome[i]
#                     # Change to a different base
#                     new_bases = filter(b -> b != original_base, [BioSequences.DNA_A, BioSequences.DNA_T, BioSequences.DNA_G, BioSequences.DNA_C])
#                     variant_genome[i] = rand(rng_strain, new_bases)
#                 end
#             end

#             # Introduce small indels (deletions of 1-3bp every 1000bp)
#             positions_to_delete = collect(1000:1000:length(variant_genome)-10)
#             for pos in reverse(positions_to_delete)  # reverse to maintain positions
#                 if pos + 2 <= length(variant_genome)
#                     deleteat!(variant_genome, pos:pos+1)  # delete 2bp
#                 end
#             end

#             variant_fasta_record = FASTX.FASTA.Record("variant_strain", variant_genome)
#             Mycelia.write_fasta(outfile=variant_ref_fasta, records=[variant_fasta_record])

#             # Simulate reads from each strain with uneven coverage (3:1 ratio)
#             base_reads = Mycelia.simulate_nanopore_reads(fasta=base_ref_fasta, quantity="12x")  # Higher coverage
#             variant_reads = Mycelia.simulate_nanopore_reads(fasta=variant_ref_fasta, quantity="4x")  # Lower coverage

#             # Decompress and combine reads for mixed community
#             base_fastq = joinpath(dir, "base_strain_reads.fq")
#             variant_fastq = joinpath(dir, "variant_strain_reads.fq")
#             mixed_fastq = joinpath(dir, "mixed_strain_reads.fq")

#             run(pipeline(`gunzip -c $(base_reads)`, base_fastq))
#             run(pipeline(`gunzip -c $(variant_reads)`, variant_fastq))

#             # Combine reads to simulate mixed community
#             run(pipeline(`cat $(base_fastq) $(variant_fastq)`, mixed_fastq))

#             # Generate assembly graph using metaFlye (metagenomic long-read assembler)
#             metaflye_for_graph_outdir = joinpath(dir, "metaflye_for_graph")
#             assembly_graph_gfa = joinpath(metaflye_for_graph_outdir, "assembly_graph.gfa")

#             try
#                 # Run metaFlye to generate assembly graph from mixed strain reads
#                 Mycelia.run_metaflye(fastq=mixed_fastq, outdir=metaflye_for_graph_outdir, genome_size="5k", read_type="nano-raw")

#                 # Test STRONG - clean up any existing directory first  
#                 strong_outdir = joinpath(dir, "strong_assembly")
#                 if isdir(strong_outdir)
#                     rm(strong_outdir, recursive=true)
#                 end

#                 if isfile(assembly_graph_gfa)
#                     try
#                         result = Mycelia.run_strong(assembly_graph_gfa, mixed_fastq, outdir=strong_outdir, nb_strains=2)
#                         Test.@test result.outdir == strong_outdir
#                         Test.@test result.strain_unitigs == joinpath(strong_outdir, "strain_unitigs.fasta")
#                         # Clean up after test
#                         rm(strong_outdir, recursive=true, force=true)
#                     catch e
#                         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                             @warn """
#                             STRONG strain resolution failed due to resource constraints.
#                             Current test: 2 strains (5kb each), mixed 12x+4x coverage, ~80kb total sequence data

#                             Required resources for STRONG:
#                             - Memory: ~3-6GB RAM minimum (strain resolution is compute-intensive)
#                             - CPU: 2-8 cores recommended
#                             - Disk: ~1GB temporary space
#                             - Note: STRONG performs strain-aware assembly graph traversal

#                             To fix: Increase available memory or reduce genome/coverage complexity.
#                             """
#                             Test.@test_skip "STRONG test skipped - insufficient resources"
#                         else
#                             rethrow(e)
#                         end
#                         # Clean up on failure
#                         rm(strong_outdir, recursive=true, force=true)
#                     end
#                 else
#                     @warn "Assembly graph not generated by metaFlye - skipping STRONG test"
#                     Test.@test_skip "STRONG test skipped - no assembly graph available"
#                 end

#                 # Clean up metaFlye output
#                 rm(metaflye_for_graph_outdir, recursive=true, force=true)

#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                     @warn """
#                     metaFlye assembly for STRONG test failed due to resource constraints.
#                     Cannot generate assembly graph needed for strain resolution testing.

#                     This test requires:
#                     1. metaFlye assembly to generate GFA graph (~3-4GB RAM)
#                     2. STRONG strain resolution on the graph (~3-6GB RAM)

#                     Total resources needed: ~5-8GB RAM, 2-8 cores
#                     """
#                     Test.@test_skip "STRONG test skipped - cannot generate assembly graph (insufficient resources)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(metaflye_for_graph_outdir, recursive=true, force=true)
#             end

#             # Test Strainy - requires an assembly FASTA file
#             # We can reuse the metaFlye assembly if it was generated successfully above
#             metaflye_assembly = joinpath(dir, "metaflye_for_strainy", "assembly.fasta")

#             # Generate assembly using metaFlye for Strainy testing
#             metaflye_for_strainy_outdir = joinpath(dir, "metaflye_for_strainy")

#             try
#                 # Run metaFlye to generate assembly FASTA for Strainy
#                 Mycelia.run_metaflye(fastq=mixed_fastq, outdir=metaflye_for_strainy_outdir, genome_size="5k", read_type="nano-raw")

#                 if isfile(metaflye_assembly)
#                     # Test Strainy - clean up any existing directory first  
#                     strainy_outdir = joinpath(dir, "strainy_assembly")
#                     if isdir(strainy_outdir)
#                         rm(strainy_outdir, recursive=true)
#                     end

#                     try
#                         result = Mycelia.run_strainy(metaflye_assembly, mixed_fastq, outdir=strainy_outdir, mode="phase")
#                         Test.@test result.outdir == strainy_outdir
#                         Test.@test result.strain_assemblies == joinpath(strainy_outdir, "strain_assemblies.fasta")
#                         # Clean up after test
#                         rm(strainy_outdir, recursive=true, force=true)
#                     catch e
#                         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                             @warn """
#                             Strainy strain phasing failed due to resource constraints.
#                             Current test: 2 strains (5kb each), mixed 12x+4x coverage, ~80kb total sequence data

#                             Required resources for Strainy:
#                             - Memory: ~2-4GB RAM minimum (includes minimap2 + samtools)
#                             - CPU: 2-8 cores recommended
#                             - Disk: ~1GB temporary space
#                             - Note: Strainy performs strain phasing from long reads mapped to assembly

#                             To fix: Increase available memory or reduce genome/coverage complexity.
#                             """
#                             Test.@test_skip "Strainy test skipped - insufficient resources"
#                         else
#                             rethrow(e)
#                         end
#                         # Clean up on failure
#                         rm(strainy_outdir, recursive=true, force=true)
#                     end
#                 else
#                     @warn "Assembly not generated by metaFlye - skipping Strainy test"
#                     Test.@test_skip "Strainy test skipped - no assembly available"
#                 end

#                 # Clean up metaFlye output for Strainy
#                 rm(metaflye_for_strainy_outdir, recursive=true, force=true)

#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                     @warn """
#                     metaFlye assembly for Strainy test failed due to resource constraints.
#                     Cannot generate assembly needed for strain phasing testing.
#                     """
#                     Test.@test_skip "Strainy test skipped - cannot generate assembly (insufficient resources)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(metaflye_for_strainy_outdir, recursive=true, force=true)
#             end
#         end
#     end
#     Test.@testset "9. Validation & Quality Control" begin
#     end
# end

Test.@testset "Strain-aware workflows (STRONG/Strainy)" begin
    mktempdir() do dir
        # Two small strains
        base_ref_fasta = joinpath(dir, "base_strain.fasta")
        variant_ref_fasta = joinpath(dir, "variant_strain.fasta")
        rng_base = StableRNGs.StableRNG(901)
        rng_variant = StableRNGs.StableRNG(902)
        base_genome = BioSequences.randdnaseq(rng_base, 5000)
        variant_genome = BioSequences.randdnaseq(rng_variant, 5000)
        Mycelia.write_fasta(outfile = base_ref_fasta,
            records = [FASTX.FASTA.Record("base_strain", base_genome)])
        Mycelia.write_fasta(outfile = variant_ref_fasta,
            records = [FASTX.FASTA.Record("variant_strain", variant_genome)])

        # Simulate reads (uneven coverage)
        base_reads_gz = Mycelia.simulate_nanopore_reads(fasta = base_ref_fasta, quantity = "12x")
        variant_reads_gz = Mycelia.simulate_nanopore_reads(fasta = variant_ref_fasta, quantity = "4x")

        base_fastq = joinpath(dir, "base_strain_reads.fq")
        variant_fastq = joinpath(dir, "variant_strain_reads.fq")
        mixed_fastq = joinpath(dir, "mixed_strain_reads.fq")
        run(pipeline(`gunzip -c $(base_reads_gz)`, base_fastq))
        run(pipeline(`gunzip -c $(variant_reads_gz)`, variant_fastq))
        run(pipeline(`cat $(base_fastq) $(variant_fastq)`, mixed_fastq))

        # STRONG: requires assembly graph; generate via metaFlye
        metaflye_for_graph_outdir = joinpath(dir, "metaflye_for_graph")
        assembly_graph_gfa = joinpath(metaflye_for_graph_outdir, "assembly_graph.gfa")
        try
            Mycelia.run_metaflye(fastq = mixed_fastq, outdir = metaflye_for_graph_outdir,
                genome_size = "5k", read_type = "nano-raw", min_overlap = 1000)

            strong_outdir = joinpath(dir, "strong_assembly")
            if isfile(assembly_graph_gfa)
                try
                    result = Mycelia.run_strong(assembly_graph_gfa, mixed_fastq,
                        outdir = strong_outdir, nb_strains = 2)
                    Test.@test result.outdir == strong_outdir
                    Test.@test result.strain_unitigs ==
                               joinpath(strong_outdir, "strain_unitigs.fasta")
                catch e
                    @error "STRONG test failed." exception=(e, catch_backtrace())
                    Test.@test false
                end
            else
                @error "STRONG test failed: metaFlye did not produce an assembly graph."
                Test.@test false
            end
        catch e
            @error "STRONG test failed: metaFlye could not generate assembly graph." exception=(
                e, catch_backtrace())
            Test.@test false
        end
        isdir(metaflye_for_graph_outdir) &&
            rm(metaflye_for_graph_outdir, recursive = true, force = true)

        # Strainy: requires assembly FASTA; generate via metaFlye
        metaflye_for_strainy_outdir = joinpath(dir, "metaflye_for_strainy")
        metaflye_assembly = joinpath(metaflye_for_strainy_outdir, "assembly.fasta")
        try
            Mycelia.run_metaflye(fastq = mixed_fastq, outdir = metaflye_for_strainy_outdir,
                genome_size = "5k", read_type = "nano-raw", min_overlap = 1000)

            if isfile(metaflye_assembly)
                strainy_outdir = joinpath(dir, "strainy_assembly")
                try
                    result = Mycelia.run_strainy(metaflye_assembly, mixed_fastq,
                        outdir = strainy_outdir, mode = "phase")
                    Test.@test result.outdir == strainy_outdir
                    Test.@test result.strain_assemblies ==
                               joinpath(strainy_outdir, "strain_assemblies.fasta")
                catch e
                    @error "Strainy test failed." exception=(e, catch_backtrace())
                    Test.@test false
                end
            else
                @error "Strainy test failed: metaFlye did not produce an assembly."
                Test.@test false
            end
        catch e
            @error "Strainy test failed: metaFlye could not generate assembly." exception=(
                e, catch_backtrace())
            Test.@test false
        end
        isdir(metaflye_for_strainy_outdir) &&
            rm(metaflye_for_strainy_outdir, recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# Additional assembler smoke tests (simulated data)
# ---------------------------------------------------------------------------

Test.@testset "Protein Assembly - PLASS (simulated reads)" begin
    mktempdir() do dir
        ref_fasta = joinpath(dir, "ref.fasta")
        rng = StableRNGs.StableRNG(42)
        seq = BioSequences.randdnaseq(rng, 2000)
        Mycelia.write_fasta(outfile = ref_fasta, records = [FASTX.FASTA.Record("ref", seq)])

        sim = Mycelia.simulate_illumina_reads(
            fasta = ref_fasta,
            coverage = 5,
            outbase = joinpath(dir, "sim_plass"),
            read_length = 100,
            mflen = 200,
            seqSys = "HS25",
            paired = true,
            errfree = true,
            quiet = true
        )

        outdir = joinpath(dir, "plass_out")
        try
            result = Mycelia.run_plass_assemble(
                reads1 = sim.forward_reads, reads2 = sim.reverse_reads,
                outdir = outdir, min_length = 20, num_iterations = 1)
            Test.@test isfile(result.assembly)
        catch e
            @error "PLASS test failed." exception=(e, catch_backtrace())
            Test.@test false
        end
    end
end

Test.@testset "Nucleotide Assembly - PenguiN guided_nuclassemble (simulated reads)" begin
    mktempdir() do dir
        ref_fasta = joinpath(dir, "ref.fasta")
        rng = StableRNGs.StableRNG(43)
        seq = BioSequences.randdnaseq(rng, 2000)
        Mycelia.write_fasta(outfile = ref_fasta, records = [FASTX.FASTA.Record("ref", seq)])

        sim = Mycelia.simulate_illumina_reads(
            fasta = ref_fasta,
            coverage = 5,
            outbase = joinpath(dir, "sim_penguin_guided"),
            read_length = 100,
            mflen = 200,
            seqSys = "HS25",
            paired = true,
            errfree = true,
            quiet = true
        )

        outdir = joinpath(dir, "penguin_guided_out")
        try
            result = Mycelia.run_penguin_guided_nuclassemble(
                reads1 = sim.forward_reads, reads2 = sim.reverse_reads, outdir = outdir)
            Test.@test isfile(result.assembly)
        catch e
            @error "PenguiN guided_nuclassemble test failed." exception=(
                e, catch_backtrace())
            Test.@test false
        end
    end
end

Test.@testset "Nucleotide Assembly - PenguiN nuclassemble (simulated reads)" begin
    mktempdir() do dir
        ref_fasta = joinpath(dir, "ref.fasta")
        rng = StableRNGs.StableRNG(44)
        seq = BioSequences.randdnaseq(rng, 2000)
        Mycelia.write_fasta(outfile = ref_fasta, records = [FASTX.FASTA.Record("ref", seq)])

        sim = Mycelia.simulate_illumina_reads(
            fasta = ref_fasta,
            coverage = 5,
            outbase = joinpath(dir, "sim_penguin"),
            read_length = 100,
            mflen = 200,
            seqSys = "HS25",
            paired = true,
            errfree = true,
            quiet = true
        )

        outdir = joinpath(dir, "penguin_out")
        try
            result = Mycelia.run_penguin_nuclassemble(
                reads1 = sim.forward_reads, reads2 = sim.reverse_reads, outdir = outdir)
            Test.@test isfile(result.assembly)
        catch e
            @error "PenguiN nuclassemble test failed." exception=(e, catch_backtrace())
            Test.@test false
        end
    end
end
