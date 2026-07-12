# Re-assembly-k must ADAPT to the residual error of the corrected reads, not pin
# to the corrector's k ceiling. High-error long reads (nanopore) leave a shattered
# k-mer spectrum at high k, so the re-assembly k must drop to keep the contig graph
# connected; clean (Illumina) corrected reads keep a high k for specificity.
#
# The criterion is COVERAGE-AWARE CONNECTIVITY (not a residual-error → survival
# chain): median_solid_kmer_multiplicity(reads, k) measures C·(1-e)^k directly (the
# median count of k-mers that recur >= 2x, i.e. the genomic peak). select_reassembly_k
# returns the LARGEST prime k (ceiling honored as-is) whose solid median stays above
# connectivity_floor. Clean 30x reads keep the ceiling; high-error nanopore reads
# see the solid median collapse at high k and drop to a lower prime.
import Test
import BioSequences
import FASTX
import Random
import Mycelia

# --- deterministic clean + nanopore-noisy read sets from a synthetic reference ---
Random.seed!(42)
const _RK_BASES = ['A', 'C', 'G', 'T']
_rk_ref = String(rand(_RK_BASES, 3000))

function _rk_tile(seq::String, readlen::Int, cov::Int)
    n = max(1, round(Int, length(seq) * cov / readlen))
    starts = rand(1:(length(seq) - readlen + 1), n)
    return [seq[s:(s + readlen - 1)] for s in starts]
end

function _rk_noisy(clean::Vector{String}, err::Float64)
    out = String[]
    for r in clean
        obs = Mycelia.observe(BioSequences.LongDNA{4}(r); error_rate = err, tech = :nanopore)
        push!(out, String(obs isa Tuple ? obs[1] : obs))
    end
    return out
end

_rk_clean = _rk_tile(_rk_ref, 1000, 30)
_rk_n05 = _rk_noisy(_rk_clean, 0.05)
_rk_n10 = _rk_noisy(_rk_clean, 0.10)

Test.@testset "Rhizomorph Re-assembly K Selection" begin
    Test.@testset "residual-error estimate rises with injected error (k-mer spectrum)" begin
        e0 = Mycelia.Rhizomorph.estimate_residual_error(_rk_clean)
        e5 = Mycelia.Rhizomorph.estimate_residual_error(_rk_n05)
        e10 = Mycelia.Rhizomorph.estimate_residual_error(_rk_n10)
        Test.@test e0 < e5
        Test.@test e5 < e10
        Test.@test e0 < 0.02      # clean reads are ~error-free
        Test.@test e10 > 0.04     # 10% nanopore leaves substantial residual error
    end

    Test.@testset "reassembly-k drops as error rises; bounded; adapted k is prime" begin
        k0 = Mycelia.Rhizomorph.select_reassembly_k(_rk_clean, 21)
        k5 = Mycelia.Rhizomorph.select_reassembly_k(_rk_n05, 21)
        k10 = Mycelia.Rhizomorph.select_reassembly_k(_rk_n10, 21)
        for k in (k0, k5, k10)
            Test.@test 7 <= k <= 21
        end
        # Clean corrected reads honor the requested ceiling exactly (identical to
        # legacy behavior; keeps graph-reuse eligible). Non-prime ceiling is allowed
        # here because it is the caller's explicit k, not an adaptive choice.
        Test.@test k0 == 21
        # High-error reads adapt DOWN, and the adapted k is prime.
        Test.@test k5 < 21
        Test.@test k10 < 21
        Test.@test Mycelia.Primes.isprime(k5)
        Test.@test Mycelia.Primes.isprime(k10)
        Test.@test k0 >= k5      # monotonic non-increasing
        Test.@test k5 >= k10
        Test.@test k10 <= 11     # heavy nanopore residual drops k low
    end

    Test.@testset "never exceeds the requested ceiling" begin
        Test.@test Mycelia.Rhizomorph.select_reassembly_k(_rk_clean, 11) <= 11
        Test.@test Mycelia.Rhizomorph.select_reassembly_k(_rk_n10, 13) <= 13
    end

    Test.@testset "connectivity signal: solid-backbone median separates clean vs shattered" begin
        # The criterion driving select_reassembly_k is median_solid_kmer_multiplicity:
        # the median count of the SOLID (count>=2) canonical k-mers = backbone coverage
        # C·(1-e)^k. A well-covered backbone stays high; a shattered one sits low. The
        # 6.0 connectivity floor lives between the two clusters (calibration finding).
        ms(reads, k) = Mycelia.Rhizomorph.median_solid_kmer_multiplicity(reads, k)
        # Clean backbone recurs ~coverage even at high k — safely above the floor.
        Test.@test ms(_rk_clean, 21) >= 6.0
        Test.@test ms(_rk_clean, 21) > 10.0
        # Shattered (raw high-error) backbone barely survives — below the floor at the
        # ceiling, which is exactly what forces the k drop.
        Test.@test ms(_rk_n05, 21) < 6.0
        Test.@test ms(_rk_n10, 21) < 6.0
        # Clean dominates noisy by a wide margin (drives the keep-vs-drop split).
        Test.@test ms(_rk_clean, 21) > ms(_rk_n05, 21)
    end

    Test.@testset "effective_kmer_coverage (secondary diagnostic) shows error-driven decay" begin
        # effective_kmer_coverage (mean over ALL distinct canonical k-mers, singletons
        # included) is retained as an interpretable cross-check: it decays smoothly with
        # k as errors fragment long windows. NOT the primary signal — its clean-vs-
        # shattered range is too narrow for a byte-identical oracle (see docstring).
        eff(reads, k) = Mycelia.Rhizomorph.effective_kmer_coverage(reads, k)
        Test.@test eff(_rk_clean, 21) > eff(_rk_n05, 21)
        Test.@test eff(_rk_n10, 7) > eff(_rk_n10, 13)
        Test.@test eff(_rk_n10, 13) > eff(_rk_n10, 21)
        Test.@test eff(_rk_n05, 7) > eff(_rk_n05, 21)
    end

    Test.@testset "Q-values inform residual error when present (FASTQ)" begin
        # identical clean sequence, differing only in per-base quality: the estimate
        # must reflect the Q-value signal (Q10 ~ 0.1 err) above the k-mer signal (~0).
        seq = _rk_ref[1:200]
        hi_q = [FASTX.FASTQ.Record("hi$i", seq, repeat("I", length(seq))) for i in 1:20]  # Q40
        lo_q = [FASTX.FASTQ.Record("lo$i", seq, repeat("+", length(seq))) for i in 1:20]  # Q10
        e_hi = Mycelia.Rhizomorph.estimate_residual_error(hi_q)
        e_lo = Mycelia.Rhizomorph.estimate_residual_error(lo_q)
        Test.@test e_lo > e_hi
    end

    Test.@testset "genome-size-aware uniqueness floor (td-jt7r)" begin
        R = Mycelia.Rhizomorph
        # A multi-kb genome's distinct k-mer set exceeds the 7-mer space (4^7=16384),
        # so the drop-down floor must rise above 7 to keep genomic k-mers ~unique.
        f_multi = R._genome_size_floor_k(_rk_clean, 7, 21)
        Test.@test 7 < f_multi <= 21
        # A tiny genome (~150 bp) fits comfortably in the 7-mer space ⇒ floor stays 7.
        tiny = [FASTX.FASTQ.Record("t$i", _rk_ref[1:150], repeat("I", 150)) for i in 1:30]
        Test.@test R._genome_size_floor_k(tiny, 7, 21) == 7
        # Never exceeds the ceiling.
        Test.@test R._genome_size_floor_k(_rk_clean, 7, 11) <= 11
        # ORACLE: clean reads still honor the ceiling — the floor is inert on the
        # no-drop path, byte-identical with the floor on or off.
        Test.@test R.select_reassembly_k(_rk_clean, 21) == 21
        Test.@test R.select_reassembly_k(_rk_clean, 21; genome_size_floor = false) == 21
        # Shattered reads never drop BELOW the genome-size floor (no k=7 collapse).
        Test.@test R.select_reassembly_k(_rk_n10, 21) >=
                   R._genome_size_floor_k(_rk_n10, 7, 21)
    end
end
