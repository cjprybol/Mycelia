# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/amino_acid_fastq_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/amino_acid_fastq_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import BioSequences
import FASTX

# Simple test to verify amino acid issue is fixed (no * characters)
Test.@testset "Amino Acid Fix Test" begin
    # Create a simple amino acid sequence
    reference_seq = BioSequences.LongAA("ACDLMVFGHY")

    # Test that the observe function doesn't produce invalid characters for FASTQ
    try
        observed_seq, quality_scores = Mycelia.observe(reference_seq, error_rate=0.1)
        observed_seq_str = string(observed_seq)

        # Check that no termination characters (*) are present
        Test.@test !occursin('*', observed_seq_str)

        # Check that the sequence only contains valid amino acids (no gaps or ambiguous characters)
        Test.@test all(c in Mycelia.UNAMBIGUOUS_AA_CHARSET for c in observed_seq_str)

        # Test that we can create a FASTQ record (this was failing before)
        quality_string = String([Char(q + 33) for q in quality_scores])
        record = FASTX.FASTQ.Record("test", observed_seq_str, quality_string)
        Test.@test record isa FASTX.FASTQ.Record

        println("Successfully created amino acid FASTQ record: ", observed_seq_str)

    catch e
        println("Error: ", e)
        Test.@test false
    end
end
