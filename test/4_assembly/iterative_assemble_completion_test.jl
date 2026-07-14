# End-to-end completion invariant for mycelia_iterative_assemble.
#
# The function must run to completion (through more than one k-step) without
# error on a toy input. This guards the `current_reads`-scoping crash class: the
# variable was local to the inner while-loop and undefined at the post-loop
# "Final improvement rate" line, so the pipeline crashed at the FIRST k-transition
# and never completed on master. A benchmark happening to run past k=3 is what
# surfaced it; this test makes "the pipeline completes at all" a guarded invariant.

import Test
import Mycelia
import FASTX
import Random

const BASES = ['A', 'C', 'G', 'T']

function write_toy_fastq(path, rng; reflen = 1000, n_reads = 100, readlen = 80, err = 0.01)
    ref = join(rand(rng, BASES, reflen))
    open(path, "w") do io
        for i in 1:n_reads
            s = rand(rng, 1:(reflen - readlen + 1))
            seq = collect(ref[s:(s + readlen - 1)])
            for j in 1:readlen
                rand(rng) < err && (seq[j] = rand(rng, filter(!=(seq[j]), BASES)))
            end
            println(io, "@r$i"); println(io, join(seq)); println(io, "+")
            println(io, String(fill('I', readlen)))
        end
    end
    return path
end

Test.@testset "mycelia_iterative_assemble completes end-to-end" begin
    dir = mktempdir()
    fastq = write_toy_fastq(joinpath(dir, "reads.fastq"), Random.MersenneTwister(42))

    # Both skip modes must run to completion through more than one k-step
    # (length(k_progression) >= 2 proves it survived the k-transition where the
    # current_reads crash used to fire).
    for skip in (false, true)
        result = Mycelia.mycelia_iterative_assemble(fastq;
            max_k = 7, max_iterations_per_k = 1, verbose = false,
            enable_checkpointing = false, skip_solid = skip,
            output_dir = joinpath(dir, "asm_$(skip)"))
        Test.@test result isa AbstractDict
        Test.@test haskey(result, :final_assembly)
        Test.@test haskey(result, :k_progression)
        Test.@test length(result[:k_progression]) >= 2   # completed past the first k
    end
end

Test.@testset "iterative finalization supports a disk-backed result" begin
    output_dir = mktempdir()
    timestamp = "20260713_120000"
    final_fastq = joinpath(output_dir, "reads_k3_iter1_$(timestamp).fastq")
    write(
        final_fastq,
        "@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nIIII\n",
    )
    history = Dict{Int, Vector{Dict{Symbol, Any}}}(
        3 => [Dict{Symbol, Any}(
            :timestamp => timestamp,
            :improvements_made => 0,
        )],
    )

    disk_backed = Mycelia.finalize_iterative_assembly(
        output_dir,
        [3],
        history,
        0.0;
        verbose = false,
        materialize_final_assembly = false,
    )
    Test.@test disk_backed[:final_assembly] === nothing
    Test.@test disk_backed[:metadata][:final_fastq_file] == final_fastq

    historical_default = Mycelia.finalize_iterative_assembly(
        output_dir,
        [3],
        history,
        0.0;
        verbose = false,
    )
    Test.@test historical_default[:final_assembly] == ["ACGT", "TGCA"]

    rm(final_fastq)
    missing_message = try
        Mycelia.finalize_iterative_assembly(
            output_dir,
            [3],
            history,
            0.0;
            verbose = false,
            materialize_final_assembly = false,
        )
        ""
    catch error
        sprint(showerror, error)
    end
    Test.@test occursin("final FASTQ does not exist for disk-backed result", missing_message)
end
