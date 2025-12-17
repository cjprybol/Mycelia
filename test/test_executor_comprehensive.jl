import Mycelia
import Test
import Dates

# Mock Conda environment setup to avoid errors
# Mock Conda environment setup to avoid errors
Mycelia.add_bioconda_env(pkg) = nothing

# Mock samtools indexing to avoid runtime errors
function Mycelia.samtools_index_fasta(;fasta::String)
    touch(fasta * ".fai")
    return nothing
end

Test.@testset "Comprehensive Executor Dry Run" begin
    
    # Setup dummy input files
    tmpdir = mktempdir()
    dummy_fq1 = joinpath(tmpdir, "reads_1.fq")
    dummy_fq2 = joinpath(tmpdir, "reads_2.fq")
    dummy_fa = joinpath(tmpdir, "ref.fasta")
    dummy_sig = joinpath(tmpdir, "query.sig")
    dummy_db = joinpath(tmpdir, "db_dir")
    mkpath(dummy_db)
    
    touch(dummy_fq1)
    touch(dummy_fq2)
    # Write some dummy content to avoid "empty file" assertions
    write(dummy_fq1, ">seq1\nACGT\n")
    write(dummy_fq2, ">seq1\nACGT\n")
    write(dummy_fa, ">ref\nACGTACGT\n")
    write(dummy_sig, "signature")
    
    # Create dummy taxonomy files for Metabuli
    taxonomy_dir = joinpath(tmpdir, "taxonomy")
    mkpath(taxonomy_dir)
    touch(joinpath(taxonomy_dir, "names.dmp"))
    touch(joinpath(taxonomy_dir, "nodes.dmp"))
    
    executor = Mycelia.Execution.CollectExecutor()
    
    Test.@testset "Assembly Tools" begin
        # Reset jobs
        empty!(executor.jobs)
        
        # MEGAHIT
        Mycelia.run_megahit(fastq1=dummy_fq1, fastq2=dummy_fq2, outdir=joinpath(tmpdir, "megahit"), executor=executor)
        
        # SPAdes & MetaSPAdes
        Mycelia.run_spades(fastq1=dummy_fq1, fastq2=dummy_fq2, outdir=joinpath(tmpdir, "spades"), executor=executor)
        Mycelia.run_metaspades(fastq1=dummy_fq1, fastq2=dummy_fq2, outdir=joinpath(tmpdir, "metaspades"), executor=executor)
        
        # Unicycler
        Mycelia.run_unicycler(short_1=dummy_fq1, short_2=dummy_fq2, long_reads=dummy_fq1, outdir=joinpath(tmpdir, "unicycler"), executor=executor)
        
        # Flye & MetaFlye
        Mycelia.run_flye(fastq=dummy_fq1, outdir=joinpath(tmpdir, "flye"), executor=executor)
        Mycelia.run_metaflye(fastq=dummy_fq1, outdir=joinpath(tmpdir, "metaflye"), executor=executor)
        
        # SKESA
        Mycelia.run_skesa(fastq1=dummy_fq1, fastq2=dummy_fq2, outdir=joinpath(tmpdir, "skesa"), executor=executor)
        
        # Canu (requires rewriting wrapper to return content properly or handle dry run)
        # Note: canu wrapper might need attention regarding return values for verification
        Mycelia.run_canu(fastq=dummy_fq1, outdir=joinpath(tmpdir, "canu"), genome_size="5m", executor=executor)
        
        # hifiasm
        Mycelia.run_hifiasm(fastq=dummy_fq1, outdir=joinpath(tmpdir, "hifiasm"), executor=executor)
        Mycelia.run_hifiasm_meta(fastq=dummy_fq1, outdir=joinpath(tmpdir, "hifiasm_meta"), executor=executor)
        
        # Plass
        Mycelia.run_plass_assemble(reads1=dummy_fq1, outdir=joinpath(tmpdir, "plass"), executor=executor)
        
        # Penguin
        Mycelia.run_penguin_guided_nuclassemble(reads1=dummy_fq1, outdir=joinpath(tmpdir, "penguin_guided"), executor=executor)
        Mycelia.run_penguin_nuclassemble(reads1=dummy_fq1, outdir=joinpath(tmpdir, "penguin_denovo"), executor=executor)
        
        # Strainy (requires BAM and FASTA)
        Mycelia.run_strainy(dummy_fa, dummy_fq1; outdir=joinpath(tmpdir, "strainy"), executor=executor)
        
        # MetaMDBG (requires reads)
        Mycelia.run_metamdbg(hifi_reads=dummy_fq1, outdir=joinpath(tmpdir, "metamdbg"), executor=executor)
        
        # Verify all tools added jobs
        # 15 tools called above
        Test.@test length(executor.jobs) == 15
        
        # Check specific job content
        megahit_job = executor.jobs[1]
        Test.@test occursin("megahit", string(megahit_job.cmd))
        Test.@test megahit_job.mem == "64G" # Default
        
        metaspades_job = executor.jobs[3]
        Test.@test occursin("metaspades", string(metaspades_job.cmd))
        Test.@test metaspades_job.cpus == Mycelia.get_default_threads()
    end
    
    Test.@testset "Pangenome Analysis" begin
        empty!(executor.jobs)
        
        # PGGB
        # construct_pangenome_pggb(genome_files::Vector{String}, output_dir::String; ...)
        Mycelia.construct_pangenome_pggb([dummy_fa, dummy_fa], joinpath(tmpdir, "pggb"), executor=executor)
        
        # Cactus
        # construct_pangenome_cactus(genome_files, genome_names, output_dir, reference_name; ...)
        Mycelia.construct_pangenome_cactus([dummy_fa, dummy_fa], ["ref1", "ref2"], joinpath(tmpdir, "cactus"), "ref1", executor=executor)
        
        Test.@test length(executor.jobs) >= 2
    end
    
    Test.@testset "Clustering" begin
        empty!(executor.jobs)
        
        # MMseqs2 easy-cluster
        Mycelia.mmseqs_easy_cluster(fasta=dummy_fa, output=joinpath(tmpdir, "clustering", "out"), executor=executor)
        
        Test.@test length(executor.jobs) == 1
        Test.@test occursin("mmseqs easy-cluster", string(executor.jobs[1].cmd))
    end
    
    Test.@testset "Alignments" begin
        empty!(executor.jobs)
        
        # Diamond
        Mycelia.run_diamond_search(query_fasta=dummy_fa, reference_fasta=dummy_fa, output_dir=joinpath(tmpdir, "diamond"), executor=executor)
        
        # MMseqs search
        Mycelia.run_mmseqs_search(query_fasta=dummy_fa, reference_fasta=dummy_fa, output_dir=joinpath(tmpdir, "mmseqs_search_2"), executor=executor)
        
        # BLASTP
        Mycelia.run_blastp_search(query_fasta=dummy_fa, reference_fasta=dummy_fa, output_dir=joinpath(tmpdir, "blastp"), executor=executor)
        
        Test.@test length(executor.jobs) >= 3
    end
    
    Test.@testset "Classification" begin
        empty!(executor.jobs)
        
        # Sourmash
        Mycelia.run_sourmash_sketch(input_files=[dummy_fa], outdir=joinpath(tmpdir, "sourmash_sketch"), executor=executor)
        # Search needs sig file
        Mycelia.run_sourmash_search(query_sig=dummy_sig, database_sig=dummy_sig, outdir=joinpath(tmpdir, "sourmash_search"), executor=executor)
        Mycelia.run_sourmash_gather(query_sig=dummy_sig, database_sig=dummy_sig, outdir=joinpath(tmpdir, "sourmash_gather"), executor=executor)
        
        # MetaPhlAn
        Mycelia.run_metaphlan(input_file=dummy_fq1, outdir=joinpath(tmpdir, "metaphlan"), executor=executor)
        
        # Metabuli
        Mycelia.run_metabuli_classify(input_files=[dummy_fq1], database_path=tmpdir, outdir=joinpath(tmpdir, "metabuli_classify"), executor=executor)
        # Metabuli (requires taxonomy)
        Mycelia.run_metabuli_build_db(reference_fasta=dummy_fa, taxonomy_dir=taxonomy_dir, outdir=joinpath(tmpdir, "metabuli_db"), executor=executor)
        
        Test.@test length(executor.jobs) == 6
        Test.@test occursin("sourmash sketch", string(executor.jobs[1].cmd))
    end
    
    println("Verified $(length(executor.jobs)) jobs.")
end
