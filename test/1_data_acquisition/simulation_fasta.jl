# ```bash
# julia --project=. --color=yes -e 'include("test/1_data_acquisition/simulation_fasta.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run from the Mycelia base directory:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/1_data_acquisition/simulation_fasta.jl", "test/1_data_acquisition", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import FASTX
import Random
import CodecZlib

const SEED = 42

Test.@testset "FASTA simulation and acquisition" begin
    Test.@testset "get_base_extension" begin
        Test.@test Mycelia.get_base_extension("foo.fasta.gz") == ".fasta"
        for ext in ["fna", "faa", "frn"]
            Test.@test Mycelia.get_base_extension("foo.$ext.gz") == "." * ext
        end

        Test.@test Mycelia.get_base_extension("foo.fasta") == ".fasta"
        for ext in ["fna", "faa", "frn"]
            Test.@test Mycelia.get_base_extension("foo.$ext") == "." * ext
        end
    end
    
    Test.@testset "random_fasta_record" begin
        for molecule in [:DNA, :RNA, :AA]
            a = Mycelia.random_fasta_record(moltype=molecule, seed=42, L=10)
            b = Mycelia.random_fasta_record(moltype=molecule, seed=42, L=10)
            Test.@test typeof(a) == typeof(b) == FASTX.FASTA.Record
            Test.@test length(FASTX.sequence(a)) == 10
            Test.@test FASTX.sequence(a) == FASTX.sequence(b)
            Test.@test FASTX.identifier(a) != FASTX.identifier(b)
            if molecule == :DNA
                Test.@test FASTX.sequence(a) == "CCGCCGCTCA"
            elseif molecule == :RNA
                Test.@test FASTX.sequence(a) == "CCGCCGCUCA"
            elseif molecule == :AA
                Test.@test FASTX.sequence(a) == "VATAGWWITI"
            end
        end
    end
    Test.@testset "dna record" begin
        dna_record = Mycelia.random_fasta_record(moltype=:DNA, seed=SEED, L = 10)
        Test.@test FASTX.sequence(dna_record) == "CCGCCGCTCA"
    end
    Test.@testset "rna record" begin
        rna_record = Mycelia.random_fasta_record(moltype=:RNA, seed=SEED, L = 10)
        Test.@test FASTX.sequence(rna_record) == "CCGCCGCUCA"
    end
    Test.@testset "aa record" begin
        aa_record = Mycelia.random_fasta_record(moltype=:AA, seed=SEED, L = 10)
        Test.@test FASTX.sequence(aa_record) == "VATAGWWITI"
    end
    Test.@testset "deterministic sequence" begin
        rec1 = Mycelia.random_fasta_record(seed=SEED, L=15)
        rec2 = Mycelia.random_fasta_record(seed=SEED, L=15)
        Test.@test FASTX.sequence(rec1) == FASTX.sequence(rec2)
        Test.@test FASTX.identifier(rec1) != FASTX.identifier(rec2)
    end
    Test.@testset "write_fasta" begin
        record = Mycelia.random_fasta_record(seed=SEED, L=10)
        fasta_file = "temp_write_fasta.fna.gz"
        try
            result = Mycelia.write_fasta(outfile=fasta_file, records=[record])
            Test.@test result == fasta_file
            Test.@test isfile(fasta_file)
            FASTX.FASTA.Reader(CodecZlib.GzipDecompressorStream(open(fasta_file))) do reader
                rec = first(reader)
                Test.@test FASTX.sequence(rec) == FASTX.sequence(record)
                Test.@test FASTX.identifier(rec) == FASTX.identifier(record)
            end
        finally
            isfile(fasta_file) && rm(fasta_file)
        end
    end
    Test.@testset "virus phiX174" begin
        phiX174_accession_id = "NC_001422"
        phiX174_assembly_id = "GCF_000007125.1"
        genome_result = Mycelia.download_genome_by_accession(accession=phiX174_accession_id)
        Test.@test basename(genome_result) == phiX174_accession_id * ".fna.gz"
        Test.@test Mycelia.get_base_extension(genome_result) == ".fna"
        rm(genome_result)
        phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="gff3,rna,cds,protein,genome,seq-report")
        # Test.@test basename(phiX174_assembly_dataset.genome) == phiX174_assembly_id * "_ViralProj14015_genomic.fna"
        Test.@test basename(phiX174_assembly_dataset.genome) == phiX174_assembly_id * "_ASM712v1_genomic.fna"
        Test.@test Mycelia.get_base_extension(phiX174_assembly_dataset.genome) == ".fna"
        Test.@test Mycelia.get_base_extension(phiX174_assembly_dataset.protein) == ".faa"
        rm(phiX174_assembly_id, recursive=true)
    end
end
