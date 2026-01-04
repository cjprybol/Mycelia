import Mycelia
import Test

# Mock Mycelia.add_bioconda_env to avoid errors if conda is not set up
Mycelia.add_bioconda_env(pkg) = nothing

Test.@testset "Executor Dry Runs" begin
    
    executor = Mycelia.Execution.CollectExecutor()
    tmpdir = mktempdir()
    fq1 = joinpath(tmpdir, "reads_1.fq")
    fq2 = joinpath(tmpdir, "reads_2.fq")
    write(fq1, ">r1\nACGT\n")
    write(fq2, ">r2\nTGCA\n")
    sour1 = joinpath(tmpdir, "in1.fa")
    sour2 = joinpath(tmpdir, "in2.fa")
    write(sour1, ">s1\nACGT\n")
    write(sour2, ">s2\nTGCA\n")
    
    Mycelia.Execution.with_executor(executor) do
        println("Testing run_megahit with collect executor...")
        
        # Test 1: Megahit
        res_mega = Mycelia.run_megahit(
            fastq1=fq1,
            fastq2=fq2, 
            outdir=joinpath(tmpdir, "test_megahit_out"),
            executor=executor,
            job_name="test_megahit",
            time="02:00:00",
            mem="32G"
        )
        
        Test.@test length(executor.jobs) == 1
        job = executor.jobs[1]
        println("Job 1 Name: ", job.name)
        println("Job 1 Cmd: \n", job.cmd)
        Test.@test job.name == "test_megahit"
        Test.@test job.time == "02:00:00"
        Test.@test occursin("megahit", job.cmd)
        
        # Test 2: MetaSPAdes
        println("\nTesting run_metaspades...")
        res_spades = Mycelia.run_metaspades(
            fastq1=fq1,
            fastq2=fq2,
            outdir=joinpath(tmpdir, "test_spades_out"),
            executor=executor,
            time="1-00:00:00"
        )
        
        Test.@test length(executor.jobs) == 2
        job2 = executor.jobs[2]
        println("Job 2 Name: ", job2.name)
        Test.@test job2.time == "1-00:00:00"
        
        # Test 3: Sourmash Sketch (loop)
        println("\nTesting run_sourmash_sketch...")
        res_sketch = Mycelia.run_sourmash_sketch(
            input_files=[sour1, sour2],
            outdir=joinpath(tmpdir, "test_sourmash_out"),
            executor=executor,
            job_name="sketch_test"
        )
        
        # Sourmash sketch refactored to submit ONE job with a script loop (or sequence)
        Test.@test length(executor.jobs) == 3
        job3 = executor.jobs[3]
        println("Job 3 Cmd:\n", job3.cmd)
        Test.@test occursin("sourmash sketch", job3.cmd)
        Test.@test occursin("in1.fa", job3.cmd)
        Test.@test occursin("in2.fa", job3.cmd)
    end
end
