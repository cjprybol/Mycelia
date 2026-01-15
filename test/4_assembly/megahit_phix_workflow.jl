# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/megahit_phix_workflow.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/megahit_phix_workflow.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia

const PHIX_ACCESSION = "NC_001422"

"""
Download PhiX174 reference to `outdir` using existing acquisition helper.
"""
function ensure_phix(outdir)
    mkpath(outdir)
    return Mycelia.download_genome_by_accession(accession=PHIX_ACCESSION, outdir=outdir, compressed=false)
end

Test.@testset "MEGAHIT + Bandage + Qualimap on PhiX" begin
    mktempdir() do dir
        phix_fasta = ensure_phix(dir)

        sim = Mycelia.simulate_illumina_reads(
            fasta=phix_fasta,
            coverage=20,
            read_length=100,
            paired=false,
            quiet=true,
            rndSeed=174
        )

        asm = Mycelia.run_megahit(fastq1=sim.forward_reads, outdir=joinpath(dir, "megahit"))
        Test.@test isfile(asm.contigs)
        Test.@test isfile(asm.fastg)
        Test.@test isfile(asm.gfa)

        map_result = Mycelia.minimap_map(
            fasta=asm.contigs,
            fastq=sim.forward_reads,
            mapping_type="sr",
            output_format="bam",
            sorted=true,
            quiet=true
        )
        run(map_result.cmd)
        Test.@test isfile(map_result.outfile)

        qc = Mycelia.run_qualimap_bamqc(bam=map_result.outfile, outdir=joinpath(dir, "qualimap"), threads=1)
        Test.@test isfile(qc.report_pdf)
        Test.@test isfile(qc.report_txt)

        coverage_df = Mycelia.parse_qualimap_contig_coverage(qc.report_txt)
        Test.@test !isempty(coverage_df)
    end
end
