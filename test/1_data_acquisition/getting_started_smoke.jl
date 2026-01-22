import Mycelia
import Test

Test.@testset "Getting Started Smoke Test" begin
    tmp_dir = mktempdir()
    fasta_path = joinpath(tmp_dir, "test_genome.fasta")

    record = Mycelia.random_fasta_record(moltype=:DNA, seed=1, L=1000)
    Mycelia.write_fasta(outfile=fasta_path, records=[record])
    Test.@test isfile(fasta_path)

    if lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"
        reads_path = Mycelia.simulate_pacbio_reads(
            fasta=fasta_path,
            quantity="1x",
            quiet=true
        )
        Test.@test isfile(reads_path)
    end
end
