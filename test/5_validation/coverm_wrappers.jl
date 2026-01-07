# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/5_validation/coverm_wrappers.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/5_validation/coverm_wrappers.jl", "test/5_validation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import DataFrames

Test.@testset "CoverM argument validation and parsing" begin
    Test.@testset "Builder helpers" begin
        contig_args = Mycelia._build_coverm_contig_args(
            bam_files=["sample.bam"],
            reference_fasta="ref.fa",
            methods=["mean", "covered_fraction"],
            threads=4,
            min_read_percent_identity=95.0,
            min_covered_fraction=0.8,
            out_path="out.tsv",
            additional_args=["--extra-flag", "value"]
        )
        Test.@test contig_args == [
            "contig",
            "--bam-files",
            "sample.bam",
            "--reference",
            "ref.fa",
            "--methods",
            "mean",
            "covered_fraction",
            "--threads",
            "4",
            "--min-read-percent-identity",
            "95.0",
            "--min-covered-fraction",
            "0.8",
            "--output-file",
            "out.tsv",
            "--extra-flag",
            "value",
        ]

        genome_args = Mycelia._build_coverm_genome_args(
            bam_files=["a.bam", "b.bam"],
            genome_fasta_files=["g1.fa"],
            genome_directory=nothing,
            genome_extension="fa",
            methods=["relative_abundance", "mean_coverage"],
            threads=2,
            out_path="genome.tsv",
            additional_args=["--min-read-aligned-percent", "50"]
        )
        Test.@test genome_args == [
            "genome",
            "--bam-files",
            "a.bam",
            "b.bam",
            "--genome-fasta-files",
            "g1.fa",
            "--methods",
            "relative_abundance",
            "mean_coverage",
            "--threads",
            "2",
            "--output-file",
            "genome.tsv",
            "--min-read-aligned-percent",
            "50",
        ]
    end

    Test.@testset "Validation failures" begin
        Test.@test_throws AssertionError Mycelia.run_coverm_contig(bam_files=String[])

        bam_path = tempname() * ".bam"
        open(bam_path, "w") do io
            write(io, "bam")
        end

        Test.@test_throws ErrorException Mycelia.run_coverm_genome(bam_files=[bam_path])
        Test.@test_throws ErrorException Mycelia.run_coverm_genome(
            bam_files=[bam_path],
            genome_fasta_files=["a.fa"],
            genome_directory="genomes"
        )
    end

    Test.@testset "Parsing cached outputs without invoking CoverM" begin
        bam_path = tempname() * ".bam"
        open(bam_path, "w") do io
            write(io, "bam")
        end
        fasta_path = tempname() * ".fa"
        open(fasta_path, "w") do io
            write(io, ">c1\nACGT\n")
        end

        contig_output = tempname() * ".tsv"
        open(contig_output, "w") do io
            write(io, "contig\tmean\tcovered_fraction\nc1\t1.0\t0.5\n")
        end
        df_contig = Mycelia.run_coverm_contig(
            bam_files=[bam_path],
            reference_fasta=fasta_path,
            output_tsv=contig_output,
            quiet=true
        )
        Test.@test df_contig isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df_contig) == 1

        genome_output = tempname() * ".tsv"
        open(genome_output, "w") do io
            write(io, "genome\trelative_abundance\tmean_coverage\nbin1\t0.8\t5.0\n")
        end
        df_genome = Mycelia.run_coverm_genome(
            bam_files=[bam_path],
            genome_fasta_files=[fasta_path],
            output_tsv=genome_output,
            quiet=true
        )
        Test.@test df_genome isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df_genome) == 1
    end
end
