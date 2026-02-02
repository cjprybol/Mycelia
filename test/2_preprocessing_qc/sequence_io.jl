# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/sequence_io.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/sequence_io.jl", "test/2_preprocessing_qc", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import FASTX
import BioSequences

Test.@testset "sequence IO" begin
    Test.@testset "detect alphabet type" begin
        Test.@test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :DNA))) == :DNA
        Test.@test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :RNA))) == :RNA
        Test.@test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :AA))) == :AA
    end
    Test.@testset "autoconvert sequences" begin
        Test.@test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :DNA)))) <: BioSequences.LongDNA{4}
        Test.@test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :RNA)))) <: BioSequences.LongRNA{4}
        Test.@test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :AA)))) <: BioSequences.LongAA
    end
    Test.@testset "detect sequence extension" begin
        Test.@test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :DNA))) == ".fna"
        Test.@test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :RNA))) == ".frn"
        Test.@test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(
            L = 100, moltype = :AA))) == ".faa"
    end
end
