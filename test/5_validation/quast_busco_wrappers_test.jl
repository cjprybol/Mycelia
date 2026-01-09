# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/5_validation/quast_busco_wrappers_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/5_validation/quast_busco_wrappers_test.jl", "test/5_validation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia

"""
QUAST and BUSCO wrapper tests.
Execution is opt-in for extended runs:
- Enable with `MYCELIA_RUN_EXTERNAL=true` (uses auto-lineage by default for BUSCO).
Wrappers auto-install required Bioconda envs.
"""

Test.@testset "QUAST wrapper" begin
    # Input validation (no external tool execution)
    Test.@test_throws ErrorException Mycelia.run_quast(["/nonexistent.fasta"])
    Test.@test_throws ErrorException Mycelia.run_quast([tempname() * ".fasta"], reference="/nonexistent_ref.fasta")

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        mktempdir() do _
            # Simulated genome smoke (portable)
            sim_info = Mycelia.get_test_genome_fasta(use_ncbi=false)
            asm_sim = sim_info.fasta
            try
                sim_outdir = Mycelia.run_quast([asm_sim]; min_contig=1, threads=1)
                expected_sim = joinpath(dirname(asm_sim), replace(basename(asm_sim), Mycelia.FASTA_REGEX => "") * "_quast")
                Test.@test sim_outdir == expected_sim
                Test.@test isfile(joinpath(sim_outdir, "report.txt"))
                Test.@test isfile(joinpath(sim_outdir, "report.tsv"))
            finally
                sim_info.cleanup()
            end

            # Real small genome (phiX174 by default) via NCBI helper; falls back to simulated if download fails
            real_info = Mycelia.get_test_genome_fasta(use_ncbi=true)
            asm_real = real_info.fasta
            try
                real_outdir = Mycelia.run_quast([asm_real]; min_contig=1, threads=1)
                expected_real = joinpath(dirname(asm_real), replace(basename(asm_real), Mycelia.FASTA_REGEX => "") * "_quast")
                Test.@test real_outdir == expected_real
                Test.@test isfile(joinpath(real_outdir, "report.txt"))
                Test.@test isfile(joinpath(real_outdir, "report.tsv"))
            finally
                real_info.cleanup()
            end
        end
    else
        Test.@test_skip "QUAST run skipped (set MYCELIA_RUN_EXTERNAL=true to enable)"
    end
end

Test.@testset "BUSCO wrapper" begin
    # Input validation (no external tool execution)
    Test.@test_throws ErrorException Mycelia.run_busco(["/nonexistent.fasta"])

    sample_list_output = """
    Available lineage datasets:
      - bacteria_odb10 (Bacteria)
      - archaea_odb10 (Archaea)
      - bacteria_odb10 (duplicate to test unique)
    """
    Test.@test Mycelia._parse_busco_dataset_list(sample_list_output) == ["bacteria_odb10", "archaea_odb10"]

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        mktempdir() do _
            genome_info = Mycelia.get_test_genome_fasta(use_ncbi=true, accession="GCF_000005845.2")
            asm = genome_info.fasta
            try
                if genome_info.source != :ncbi
                    Test.@test_skip "BUSCO run skipped (NCBI genome download failed)"
                else
                    result_dir = Mycelia.run_busco([asm]; force=true, lineage="bacteria_odb10", auto_lineage=false)
                    expected_outdir = joinpath(dirname(asm), replace(basename(asm), Mycelia.FASTA_REGEX => "") * "_busco")
                    Test.@test result_dir == expected_outdir
                    # BUSCO writes summary inside outdir/assembly_name/
                    assembly_name = replace(basename(asm), Mycelia.FASTA_REGEX => "")
                    summary_dir = joinpath(expected_outdir, assembly_name)
                    has_summary = any(isfile, [
                        joinpath(summary_dir, "short_summary.txt"),
                        joinpath(summary_dir, "short_summary.specific.auto_lineage.txt"),
                        joinpath(summary_dir, "short_summary.specific.auto_lineage_prok.txt"),
                        joinpath(summary_dir, "short_summary.specific.auto_lineage_euk.txt"),
                        joinpath(summary_dir, "short_summary.specific.bacteria_odb10.txt")
                    ])
                    if !has_summary && isdir(summary_dir)
                        for (_, _, files) in walkdir(summary_dir)
                            if any(name -> occursin(r"^short_summary.*\.txt$", name), files)
                                has_summary = true
                                break
                            end
                        end
                    end
                    Test.@test has_summary
                end
            finally
                genome_info.cleanup()
            end

            list_info = Mycelia.list_busco_datasets()
            Test.@test !isempty(strip(list_info.raw))
            Test.@test isa(list_info.datasets, Vector{String})
        end
    else
        Test.@test_skip "BUSCO run skipped (set MYCELIA_RUN_EXTERNAL=true to enable)"
    end
end
