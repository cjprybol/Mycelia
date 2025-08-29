# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/preprocessing.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/preprocessing.jl", "test/2_preprocessing_qc", execute=false)'
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
import CSV

const phiX174_accession_id = "NC_001422.1"

Test.@testset "Preprocessing" begin
    Test.@testset "FASTX stats" begin
        srr_identifier = "SRR31812976"
        outdir = mkpath("fastx-stats-test")
        fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
        Test.@test fasterq_dump_result.unpaired_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier).fastq.gz"
        table = Mycelia.fastx_stats(fasterq_dump_result.unpaired_reads)
        io = IOBuffer()
        CSV.write(io, table)
        Test.@test String(take!(io)) == 
        """
        file,format,type,num_seqs,sum_len,min_len,avg_len,max_len,Q1,Q2,Q3,sum_gap,N50,N50_num,Q20(%),Q30(%),AvgQual,GC(%),sum_n,N90
        fastx-stats-test/SRR31812976/SRR31812976.fastq.gz,FASTQ,DNA,58982,34921275,80,592.1,2315,341.0,544.0,811.0,0,743,608,50.26,9.75,11.4,39.62,0,329
        """
        rm(outdir, recursive=true)
    end
    Test.@testset "fastx2normalized_table" begin
        genome_result = Mycelia.download_genome_by_accession(accession=phiX174_accession_id)
        fastx_table = Mycelia.fastx2normalized_table(genome_result)
        io = IOBuffer()
        CSV.write(io, fastx_table, delim='\t')
        Test.@test String(take!(io)) == 
        """
        fastx_path\tfastx_sha256\trecord_identifier\trecord_description\trecord_sha256\trecord_quality\trecord_alphabet\trecord_type\tmean_record_quality\tmedian_record_quality\trecord_length\trecord_sequence
        NC_001422.1.fna.gz\tae3c70dc4d180aab0aa8d6ab81c0048998575c05e402639ffe6b99e82815f953\tNC_001422.1\tNC_001422.1 Escherichia phage phiX174, complete genome\t97038c7e1edea2297667d7f0426ba942b322c74cb30e072ec66ba47f9c0448d0\t\tACGT\tDNA\t\t5386\tGAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGATGCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTTTCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTTCGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTACGGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTACGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAGTGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAAGGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGCCGCTTTTCTTGGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACAACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTTTCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAACCTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTTGATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACAGACCTATAAACATTCTGTGCCGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGACTAAAGAGATTCAGTACCTTAACGCTAAAGGTGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTGTATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGAAGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCGCCACCATGATTATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGATTGAGTGTGAGGTTATAACGCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGCTTAGGAGTTTAATCATGTTTCAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGTTCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTATATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTGAATGGTCGCCATGATGGTGGTTATTATACCGTCAAGGACTGTGTGACTATTGACGTCCTTCCCCGTACGCCGGGCAATAACGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAATCAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTGCTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAGGCGGTCAAAAAGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCTGGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACCCTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGATAAGTTGCTTGATTTGGTTGGACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTATCTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGGTTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTTCACGCCAGAATACGAAAGACCAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTATGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCAAACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGACTTAGTTCATCAGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTTCTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAGCTGTTGCCGATACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCGTCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTCCAAATCTTGGAGGCTTTTTTATGGTTCGTTCTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTATTGAGGCTTGTGGCATTTCTACTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGCATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATGTTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGAATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTATAGACCACCGCCCCGAAGGGGACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCCCTCTTAAGGATATTCGCGATGAGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATTGCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAGGCTCATGCTGATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTTATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGGTCGCAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGCCGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTCCTACAGGTAGCGTTGACCCTAATTTTGGTCGTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAGCCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATATGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACTTCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGCAGAACGCAAAAAGAGAGATGAGATTGAGGCTGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACCTGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA
        """
        rm(genome_result)
    end
    Test.@testset "Read Quality Control" begin
        Test.@test true # https://github.com/OpenGene/fastp
        Test.@test true # https://github.com/FelixKrueger/TrimGalore
        Test.@test true # https://github.com/rrwick/Filtlong
        Test.@test true # https://github.com/OpenGene/fastplong
        Test.@test true # https://github.com/wdecoster/chopper
        Test.@test true  # placeholder
    end
    Test.@testset "Read Statistics" begin
        Test.@test true  # placeholder
    end
    Test.@testset "fastx2normalized_table - compressed and raw input" begin
        # Prepare test records
        fasta_file = "test-normalized.fasta"
        fastq_file = "test-normalized.fastq"
        gz_fasta = fasta_file * ".gz"
        gz_fastq = fastq_file * ".gz"
        records = [
            Mycelia.random_fasta_record(moltype=:DNA, L=10, seed=1),
            Mycelia.random_fasta_record(moltype=:DNA, L=10, seed=2)
        ]
        # Write FASTA
        open(fasta_file, "w") do io
            for rec in records
                FASTX.write(io, rec)
            end
        end
        # Write FASTQ
        open(fastq_file, "w") do io
            for rec in records
                qual = "I"^length(FASTX.sequence(rec))
                fqrec = FASTX.FASTQ.Record(FASTX.identifier(rec), FASTX.sequence(rec), qual)
                FASTX.write(io, fqrec)
            end
        end
        # Compress files
        run(pipeline(`gzip -c $fasta_file`, stdout=gz_fasta))
        run(pipeline(`gzip -c $fastq_file`, stdout=gz_fastq))

        # Test raw and gzipped FASTA
        for file in (fasta_file, gz_fasta)
            table = Mycelia.fastx2normalized_table(file)
            Test.@test nrow(table) == 2
            # Case-insensitive: sequences should match after uppercasing
            seqs = [uppercase(row.record_sequence) for row in eachrow(table)]
            Test.@test length(unique(seqs)) == 2
        end

        # Test raw and gzipped FASTQ
        for file in (fastq_file, gz_fastq)
            table = Mycelia.fastx2normalized_table(file)
            Test.@test nrow(table) == 2
            # Case-sensitive: sequences should match exactly
            seqs = [row.record_sequence for row in eachrow(table)]
            Test.@test length(unique(seqs)) == 2
        end

        # Clean up
        rm(fasta_file, force=true)
        rm(fastq_file, force=true)
        rm(gz_fasta, force=true)
        rm(gz_fastq, force=true)
    end

    Test.@testset "fastx2normalized_table - identifier normalization" begin
        # Two records with different identifiers but same sequence (case-insensitive)
        fasta_file = "test-idnorm.fasta"
        open(fasta_file, "w") do io
            FASTX.write(io, FASTX.FASTA.Record("id1", "acgtACGT"))
            FASTX.write(io, FASTX.FASTA.Record("id2", "ACGTacgt"))
        end
        table = Mycelia.fastx2normalized_table(fasta_file)
        seqs = [uppercase(row.record_sequence) for row in eachrow(table)]
        Test.@test length(unique(seqs)) == 1
        rm(fasta_file, force=true)
    end

    Test.@testset "fastx_stats and fastx2normalized_table - fixtures" begin
        mktempdir() do dir
            fasta_path = joinpath(dir, "small.fasta")
            fastq_path = joinpath(dir, "small.fastq")
            records = [
                FASTX.FASTA.Record("seq1", "ACGT"),
                FASTX.FASTA.Record("seq2", "GGGG"),
            ]
            open(fasta_path, "w") do io
                for rec in records
                    FASTX.write(io, rec)
                end
            end
            open(fastq_path, "w") do io
                for rec in records
                    qual = repeat("I", length(FASTX.sequence(rec)))
                    fqrec = FASTX.FASTQ.Record(FASTX.identifier(rec), FASTX.sequence(rec), qual)
                    FASTX.write(io, fqrec)
                end
            end

            expected_cols = [
                "file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len",
                "max_len", "Q1", "Q2", "Q3", "sum_gap", "N50", "N50_num", "Q20(%)",
                "Q30(%)", "AvgQual", "GC(%)", "sum_n", "N90",
            ]

            stats_fasta = Mycelia.fastx_stats(fasta_path)
            Test.@test names(stats_fasta) == expected_cols
            Test.@test stats_fasta.num_seqs[1] == 2
            Test.@test stats_fasta.sum_len[1] == 8
            Test.@test isapprox(stats_fasta[1, "GC(%)"], 75.0; atol=0.01)

            stats_fastq = Mycelia.fastx_stats(fastq_path)
            Test.@test stats_fastq.num_seqs[1] == 2
            Test.@test isapprox(stats_fastq.AvgQual[1], 40.0; atol=0.01)

            norm_cols = [
                "fastx_path", "fastx_sha256", "record_identifier", "record_description",
                "record_sha256", "record_quality", "record_alphabet", "record_type",
                "mean_record_quality", "median_record_quality", "record_length", "record_sequence",
            ]

            nfasta = Mycelia.fastx2normalized_table(fasta_path)
            nfastq = Mycelia.fastx2normalized_table(fastq_path)
            Test.@test names(nfasta) == norm_cols
            Test.@test names(nfastq) == norm_cols
            Test.@test nrow(nfasta) == 2
            gc_val = sum(count(c -> c in ['G','C'], row.record_sequence) for row in eachrow(nfasta)) / sum(nfasta.record_length) * 100
            Test.@test isapprox(gc_val, 75.0; atol=0.01)
            Test.@test isapprox(mean(skipmissing(nfastq.mean_record_quality)), 40.0; atol=0.01)
        end
    end
end
