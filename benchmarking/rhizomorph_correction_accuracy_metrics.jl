# Pure per-base correction-accuracy metrics — NO Mycelia dependency.
# ==================================================================
#
# Factored out of rhizomorph_correction_accuracy_benchmark.jl so the metric
# classification math can be unit-tested in milliseconds without loading Mycelia
# (the benchmark's simulator + corrector integration is tested separately). This
# file depends only on Julia Base.
#
# `per_base_metrics` joins corrected reads to per-read truth BY READ ID and
# classifies every base of every scored read. Reads absent from the corrected
# set (`reads_dropped`) and reads whose corrected length != truth length
# (`reads_excluded_len`, a Viterbi indel path) are counted and EXCLUDED from the
# positional metric rather than silently miscounted. `corrected_unjoined` counts
# corrected reads whose id is NOT in truth — a non-zero value alongside
# reads_scored==0 is the signature of an ID-FORMAT mismatch (the corrector
# renamed reads) rather than genuine drops, which the caller warns on.
#
# Position classification for a scored read (truth T, observed O, corrected C, |T|=|O|=|C|):
#   injected error   : O[p] != T[p]
#   edit             : C[p] != O[p]
#   true_fix         : O[p] != T[p]  AND  C[p] == T[p]     (error corrected to truth)
#   mis_fix          : O[p] != T[p]  AND  C[p] != O[p]  AND  C[p] != T[p]  (error edited to a WRONG base)
#   over_correction  : O[p] == T[p]  AND  C[p] != T[p]     (a correct base made wrong)
#
# These partition every edit exactly: total_edits == true_fixes + mis_fixes + over_corrections.
# (An injected error left unchanged, C[p]==O[p], is a miss: it is neither an edit
# nor a fix, so it depresses recall without touching precision.)
#
# Metrics (aggregated over all scored reads):
#   recall          = true_fixes / injected_errors
#   precision       = true_fixes / total_edits
#   over_correction = over_corrections / correct_positions
#   correction_rate = total_edits / total_bases            (compare to error_rate)
# Zero-denominator cells return NaN (an "undefined here — denominator was zero"
# signal, NOT a measured value); downstream readers must treat NaN as missing.

"""
    per_base_metrics(truth_by_id, observed_by_id, corrected_by_id) -> NamedTuple

Compute per-base correction metrics by joining `corrected_by_id` to
`truth_by_id` / `observed_by_id` on read id. See the file header for the exact
position classification and the `total_edits == true_fixes + mis_fixes +
over_corrections` partition invariant. Pure (Julia Base only).
"""
function per_base_metrics(truth_by_id::Dict{String, String},
        observed_by_id::Dict{String, String}, corrected_by_id::Dict{String, String})
    tp = 0                 # true fixes  (O!=T, C==T)
    mis = 0                # mis-fixes   (O!=T, edited, C!=T)
    over = 0               # over-corrections (O==T, C!=T)
    total_edits = 0        # positions where C != O
    injected = 0           # positions where O != T
    correct_positions = 0  # positions where O == T
    total_bases = 0
    reads_scored = 0
    reads_excluded_len = 0
    reads_dropped = 0
    for (rid, T) in truth_by_id
        O = observed_by_id[rid]
        if !haskey(corrected_by_id, rid)
            reads_dropped += 1
            continue
        end
        C = corrected_by_id[rid]
        # collect to Char vectors: ACGT are ASCII, but this is robust to any
        # codeunit subtlety and lets us index positionally.
        Tc = collect(T)
        Oc = collect(O)
        Cc = collect(C)
        if length(Cc) != length(Tc)
            reads_excluded_len += 1
            continue
        end
        reads_scored += 1
        L = length(Tc)
        total_bases += L
        for p in 1:L
            terr = Oc[p] != Tc[p]
            edited = Cc[p] != Oc[p]
            terr ? (injected += 1) : (correct_positions += 1)
            edited && (total_edits += 1)
            if terr
                if Cc[p] == Tc[p]
                    tp += 1               # error corrected to truth
                elseif edited
                    mis += 1              # error edited to a wrong base
                end
            elseif edited                 # !terr && C!=O  ==>  C!=T (since O==T)
                over += 1                 # a correct base made wrong
            end
        end
    end
    # Corrected reads whose id is NOT in truth: with an exact 1:1 substitution
    # simulation this is 0; a large value alongside reads_scored==0 means the
    # corrector renamed ids (join is broken), which the caller diagnoses.
    corrected_unjoined = count(rid -> !haskey(truth_by_id, rid), keys(corrected_by_id))
    n_reads = length(truth_by_id)
    scored_fraction = n_reads == 0 ? NaN : reads_scored / n_reads
    recall = injected == 0 ? NaN : tp / injected
    precision = total_edits == 0 ? NaN : tp / total_edits
    over_rate = correct_positions == 0 ? NaN : over / correct_positions
    correction_rate = total_bases == 0 ? NaN : total_edits / total_bases
    return (; tp, mis_fixes = mis, over, total_edits, injected, correct_positions,
        total_bases, recall, precision, over_rate, correction_rate,
        reads_scored, reads_excluded_len, reads_dropped, corrected_unjoined,
        n_reads, scored_fraction)
end
