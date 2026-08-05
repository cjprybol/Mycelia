# Constant-quality positive control for the qualmer decoder arm (bead td-4e19d.2)
#
# QUESTION
#   Does per-base Phred quality reach any decision in the DEFAULT (greedy) qualmer
#   assembly path, or is the quality channel inert?
#
# DESIGN
#   Hold the read SEQUENCES fixed and vary ONLY the per-base quality string, then
#   compare the assemblies byte-for-byte. Because the sequences are identical across
#   conditions, any difference in output is attributable to quality alone; conversely
#   byte-identical output across wildly different quality vectors proves quality never
#   reaches a decision.
#
#   Conditions per cell (identical reads, four quality vectors):
#     oracle      per-base Phred that is CORRECT by construction — a base that was
#                 corrupted by the simulator gets a low Q, an uncorrupted base a high Q.
#     constant40  every base Q40.
#     constant02  every base Q2.
#     antioracle  the oracle inverted — corrupted bases get HIGH Q, correct bases LOW Q.
#
#   `oracle` is strictly MORE informative than any real instrument's Phred estimate, and
#   `antioracle` is maximally adversarial. An assembler that consumes quality at all
#   cannot produce the same contigs under all four. This makes the control decisive in
#   the positive direction: identity across conditions is proof of inertness, not merely
#   absence of evidence.
#
#   SENSITIVITY CONTROL: a fifth run mutates a single base in a single READ (quality left
#   at oracle). Its digest MUST differ from the oracle digest, otherwise the comparator is
#   not measuring anything and every "identical" verdict above would be vacuous.
#
#   CROSS-ARM CHECK: the same reads are also assembled with quality stripped to FASTA
#   (the Track-A `kmer` arm), reproducing in-process the arm identity observed in the
#   committed 66-pair Track-A pilot table.
#
# WHY A SELF-CONTAINED READ SIMULATOR
#   The Track-A harness calls ART / pbsim / badread through Conda. Those are unavailable
#   on hosts without a Conda root, and they do not let the caller specify which bases were
#   corrupted — which is exactly what the oracle-quality condition requires. This script
#   therefore simulates reads in pure Julia with StableRNGs, so it is deterministic,
#   dependency-free, and reproducible by anyone with the repo checked out. Error profiles
#   are stratified by chemistry and NEVER pooled.
#
# USAGE
#   julia --project=. benchmarking/qualmer_quality_channel_probe.jl            # full grid
#   julia --project=. benchmarking/qualmer_quality_channel_probe.jl --smoke    # 1 cell
#   julia --project=. benchmarking/qualmer_quality_channel_probe.jl --output-dir DIR

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import BioSequences
import FASTX
import SHA
import StableRNGs
import Statistics
import Dates

const REFERENCE_FASTA = joinpath(
    @__DIR__, "fixtures", "viterbi_accuracy", "phix174_nc001422.fasta")

const K = 31

# Per-chemistry error/length profiles. Deliberately separate rows in every output table:
# ONT and Illumina are different measurement processes and are never pooled.
#
# `sub`/`ins`/`del` are per-base probabilities. Illumina is substitution-dominated;
# ONT and PacBio CLR carry substantial indels, which is the regime where a quality-aware
# decoder would be expected to help most.
const CHEMISTRIES = (
    (name = "illumina", read_length = 150, sub = 0.002, ins = 0.0000, del = 0.0000,
        q_correct = 37, q_error = 12),
    (name = "pacbio", read_length = 2000, sub = 0.015, ins = 0.008, del = 0.008,
        q_correct = 25, q_error = 6),
    (name = "ont", read_length = 3000, sub = 0.030, ins = 0.015, del = 0.015,
        q_correct = 18, q_error = 4)
)

# Low coverage is included on purpose: it is where per-observation quality carries the
# most decision weight, because there are too few observations for counts alone to
# separate signal from error.
const COVERAGES = [5, 10, 30]
const SEEDS = [42, 123]

const CONDITIONS = ["oracle", "constant40", "constant02", "antioracle"]

# === Argument parsing ===

arg_value(flag) =
    let i = findfirst(==(flag), ARGS)
        (i !== nothing && i < length(ARGS)) ? ARGS[i + 1] : nothing
    end

const SMOKE = "--smoke" in ARGS
const OUTPUT_DIR = something(arg_value("--output-dir"),
    joinpath(@__DIR__, "results", "qualmer_quality_channel_probe"))

chemistries = SMOKE ? CHEMISTRIES[1:1] : CHEMISTRIES
coverages = SMOKE ? [10] : COVERAGES
seeds = SMOKE ? [42] : SEEDS

# === Reference ===

function load_reference(path::String)
    isfile(path) || error("reference FASTA not found: $(path)")
    reader = FASTX.FASTA.Reader(open(path))
    record = first(reader)
    close(reader)
    return BioSequences.LongDNA{4}(FASTX.sequence(record))
end

# === Read simulation with a per-base corruption oracle ===

const DNA_BASES = (BioSequences.DNA_A, BioSequences.DNA_C,
    BioSequences.DNA_G, BioSequences.DNA_T)

"""
Simulate one read from `reference`, returning the read sequence together with a
per-emitted-base `corrupted` mask. The mask is the ground truth the oracle-quality
condition encodes; it is available only because the simulator lives here rather than
in an external binary.

Insertions emit a base with no reference support (marked corrupted); deletions skip a
reference base (emitting nothing); substitutions emit a wrong base (marked corrupted).
"""
function simulate_read(rng, reference, profile; artifact_loci = nothing)
    reflen = length(reference)
    span = min(profile.read_length, reflen)
    start = rand(rng, 1:(reflen - span + 1))
    # Reads come off both strands, matching a real library.
    forward = rand(rng, Bool)

    template = copy(reference[start:(start + span - 1)])

    # Optional SYSTEMATIC artifact (td-2tg8): a fixed reference locus mis-called the same
    # way in EVERY read covering it. Applied here — in reference coordinates, before
    # strand selection and before random error injection — because that is what makes it
    # correlated across reads rather than just more random error. Applying it at a
    # read-relative position instead (e.g. each read's midpoint) would scatter it across
    # different loci and produce independent errors, which is the opposite of the point.
    #
    # `artifact_loci === nothing` draws no random numbers and mutates nothing, so the
    # default path is byte-identical to the pre-existing simulator.
    artifact_mask = falses(span)
    if artifact_loci !== nothing
        for locus in artifact_loci
            if start <= locus <= (start + span - 1)
                offset = locus - start + 1
                original = template[offset]
                template[offset] = original == BioSequences.DNA_A ?
                                   BioSequences.DNA_C : BioSequences.DNA_A
                artifact_mask[offset] = true
            end
        end
    end

    if !forward
        template = BioSequences.reverse_complement(template)
        artifact_mask = reverse(artifact_mask)
    end

    bases = BioSequences.DNA[]
    corrupted = Bool[]

    for (index, base) in enumerate(template)
        # Deletion: emit nothing for this reference position.
        if rand(rng) < profile.del
            continue
        end
        # An artifact base is emitted as-is and flagged corrupted, so oracle quality
        # reports it as low-confidence. It deliberately does NOT consume a substitution
        # draw: the artifact is systematic, not a random error at that position.
        if artifact_mask[index]
            push!(bases, base)
            push!(corrupted, true)
        elseif rand(rng) < profile.sub
            # Substitution: emit a different base.
            alternatives = filter(!=(base), collect(DNA_BASES))
            push!(bases, rand(rng, alternatives))
            push!(corrupted, true)
        else
            push!(bases, base)
            push!(corrupted, false)
        end
        # Insertion: emit an extra, unsupported base.
        if rand(rng) < profile.ins
            push!(bases, rand(rng, DNA_BASES))
            push!(corrupted, true)
        end
    end

    return (sequence = BioSequences.LongDNA{4}(bases), corrupted = corrupted)
end

"""
Simulate a read set at `coverage` x over `reference` for one chemistry.
Returns a vector of (identifier, sequence, corrupted-mask) triples. Quality is NOT
assigned here — each condition derives its own quality vector from the same masks, so
every condition sees byte-identical read sequences.
"""
function simulate_read_set(reference, profile, coverage::Int, seed::Int;
        artifact_loci = nothing)
    rng = StableRNGs.StableRNG(seed)
    target_bases = length(reference) * coverage
    reads = NamedTuple{
        (:identifier, :sequence, :corrupted),
        Tuple{String, BioSequences.LongDNA{4}, Vector{Bool}}}[]
    emitted = 0
    i = 0
    while emitted < target_bases
        i += 1
        read = simulate_read(rng, reference, profile; artifact_loci = artifact_loci)
        isempty(read.sequence) && continue
        push!(reads,
            (identifier = "$(profile.name)_read_$(i)",
                sequence = read.sequence, corrupted = read.corrupted))
        emitted += length(read.sequence)
    end
    return reads
end

# === Quality assignment per condition ===

"""
Build the per-base Phred vector for one read under one condition.

`oracle` and `antioracle` both consume the ground-truth corruption mask; the constants
ignore it entirely. All four return a vector the same length as the read, so the only
thing that varies downstream is the numeric quality itself.
"""
function quality_for(condition::String, corrupted::Vector{Bool}, profile)
    n = length(corrupted)
    if condition == "oracle"
        return UInt8[c ? profile.q_error : profile.q_correct for c in corrupted]
    elseif condition == "antioracle"
        return UInt8[c ? profile.q_correct : profile.q_error for c in corrupted]
    elseif condition == "constant40"
        return fill(UInt8(40), n)
    elseif condition == "constant02"
        return fill(UInt8(2), n)
    else
        error("unknown condition: $(condition)")
    end
end

function fastq_records(reads, condition::String, profile)
    return [FASTX.FASTQ.Record(r.identifier, r.sequence,
                quality_for(condition, r.corrupted, profile))
            for r in reads]
end

fasta_records(reads) = [FASTX.FASTA.Record(r.identifier, r.sequence) for r in reads]

# === Assembly + digest ===

"""
Digest an assembly by its contig CONTENT, independent of emission order.

Contigs are canonicalised (each replaced by the lexicographic minimum of itself and its
reverse complement) and sorted before hashing, so a difference in digest means a genuine
difference in assembled sequence, not a reshuffle or a strand flip.
"""
function assembly_digest(contigs::Vector{String})
    canonical = map(contigs) do contig
        seq = BioSequences.LongDNA{4}(contig)
        rc = string(BioSequences.reverse_complement(seq))
        min(contig, rc)
    end
    sort!(canonical)
    return SHA.bytes2hex(SHA.sha256(join(canonical, "\n")))
end

function n50(lengths::Vector{Int})
    isempty(lengths) && return 0
    sorted = sort(lengths; rev = true)
    half = sum(sorted) / 2
    running = 0
    for len in sorted
        running += len
        running >= half && return len
    end
    return sorted[end]
end

function assemble_and_summarise(records)
    elapsed = @elapsed result = Mycelia.Rhizomorph.assemble_genome(
        records; k = K, verbose = false)
    contigs = String.(result.contigs)
    lengths = length.(contigs)
    return (
        digest = assembly_digest(contigs),
        n_contigs = length(contigs),
        total_length = sum(lengths; init = 0),
        largest_contig = isempty(lengths) ? 0 : maximum(lengths),
        n50 = n50(lengths),
        wall_seconds = elapsed
    )
end

# === Main ===

function main()
    mkpath(OUTPUT_DIR)
    reference = load_reference(REFERENCE_FASTA)
    @info "probe start" reference_bp=length(reference) k=K output_dir=OUTPUT_DIR

    rows = NamedTuple[]
    verdicts = NamedTuple[]

    for profile in chemistries, coverage in coverages, seed in seeds
        cell = "$(profile.name)_cov$(coverage)_seed$(seed)"
        reads = simulate_read_set(reference, profile, coverage, seed)
        corrupted_fraction = let total = sum(length(r.corrupted) for r in reads)
            total == 0 ? 0.0 : sum(count(r.corrupted) for r in reads) / total
        end
        @info "cell" cell n_reads=length(reads) corrupted_fraction

        summaries = Dict{String, Any}()
        for condition in CONDITIONS
            records = fastq_records(reads, condition, profile)
            summary = assemble_and_summarise(records)
            summaries[condition] = summary
            push!(rows,
                (chemistry = profile.name, coverage = coverage, seed = seed,
                    arm = "qualmer", condition = condition,
                    n_reads = length(reads),
                    corrupted_base_fraction = round(corrupted_fraction; digits = 5),
                    digest = summary.digest, n_contigs = summary.n_contigs,
                    total_length = summary.total_length,
                    largest_contig = summary.largest_contig, n50 = summary.n50,
                    wall_seconds = round(summary.wall_seconds; digits = 3)))
        end

        # Cross-arm: quality stripped to FASTA == Track-A `kmer` arm.
        kmer_summary = assemble_and_summarise(fasta_records(reads))
        summaries["kmer_arm_fasta"] = kmer_summary
        push!(rows,
            (chemistry = profile.name, coverage = coverage, seed = seed,
                arm = "kmer", condition = "fasta_no_quality",
                n_reads = length(reads),
                corrupted_base_fraction = round(corrupted_fraction; digits = 5),
                digest = kmer_summary.digest, n_contigs = kmer_summary.n_contigs,
                total_length = kmer_summary.total_length,
                largest_contig = kmer_summary.largest_contig, n50 = kmer_summary.n50,
                wall_seconds = round(kmer_summary.wall_seconds; digits = 3)))

        # SENSITIVITY CONTROL: perturb one base of one read, keep oracle quality.
        # If this does NOT change the digest the comparator is blind and no verdict
        # above can be trusted.
        mutated = deepcopy(reads)
        target = findfirst(r -> length(r.sequence) > 0, mutated)
        control_digest = if target === nothing
            "NA"
        else
            seq = copy(mutated[target].sequence)
            pos = cld(length(seq), 2)
            alternatives = filter(!=(seq[pos]), collect(DNA_BASES))
            seq[pos] = first(alternatives)
            mutated[target] = (identifier = mutated[target].identifier,
                sequence = seq, corrupted = mutated[target].corrupted)
            assemble_and_summarise(fastq_records(mutated, "oracle", profile)).digest
        end
        push!(rows,
            (chemistry = profile.name, coverage = coverage, seed = seed,
                arm = "qualmer", condition = "sensitivity_control_1base",
                n_reads = length(mutated),
                corrupted_base_fraction = round(corrupted_fraction; digits = 5),
                digest = control_digest, n_contigs = -1, total_length = -1,
                largest_contig = -1, n50 = -1, wall_seconds = -1.0))

        oracle_digest = summaries["oracle"].digest
        push!(verdicts,
            (chemistry = profile.name, coverage = coverage, seed = seed,
                n_reads = length(reads),
                corrupted_base_fraction = round(corrupted_fraction; digits = 5),
                quality_conditions_identical = all(
                    summaries[c].digest == oracle_digest for c in CONDITIONS),
                kmer_arm_identical = summaries["kmer_arm_fasta"].digest == oracle_digest,
                # THREE states, not two. `control_digest == "NA"` means the control
                # never RAN, and `"NA" != oracle_digest` is `true` — so the one failure
                # mode of the control that validates every other verdict here used to
                # report as PASS. A control that did not run is not a control that
                # fired; it invalidates the cell.
                sensitivity_control_state = control_digest == "NA" ? "not_run" :
                                            (control_digest != oracle_digest ? "fired" :
                                             "did_not_fire"),
                sensitivity_control_detected_change = control_digest != "NA" &&
                                                      control_digest != oracle_digest,
                oracle_n_contigs = summaries["oracle"].n_contigs,
                oracle_n50 = summaries["oracle"].n50))
    end

    write_tsv(joinpath(OUTPUT_DIR, "qualmer_quality_channel_probe_runs.tsv"), rows)
    write_tsv(joinpath(OUTPUT_DIR, "qualmer_quality_channel_probe_verdicts.tsv"), verdicts)

    report(verdicts)
    return verdicts
end

function write_tsv(path::String, rows::Vector{<:NamedTuple})
    isempty(rows) && return nothing
    open(path, "w") do io
        cols = keys(rows[1])
        println(io, join(String.(collect(cols)), "\t"))
        for row in rows
            println(io, join((string(getproperty(row, c)) for c in cols), "\t"))
        end
    end
    @info "wrote" path rows=length(rows)
    return nothing
end

function report(verdicts)
    println("\n=== qualmer quality-channel probe — verdicts (stratified by chemistry) ===")
    for chem in unique(v.chemistry for v in verdicts)
        cells = filter(v -> v.chemistry == chem, verdicts)
        println("\n[$(chem)]  cells=$(length(cells))")
        for v in cells
            println("  cov=$(lpad(v.coverage, 3))x seed=$(v.seed) " *
                    "err=$(rpad(v.corrupted_base_fraction, 7)) " *
                    "quality_inert=$(v.quality_conditions_identical) " *
                    "kmer_arm_identical=$(v.kmer_arm_identical) " *
                    "control_fired=$(v.sensitivity_control_detected_change) " *
                    "(contigs=$(v.oracle_n_contigs) N50=$(v.oracle_n50))")
        end
    end

    controls_ok = all(v.sensitivity_control_detected_change for v in verdicts)
    all_inert = all(v.quality_conditions_identical for v in verdicts)
    println("\nsensitivity controls all fired: $(controls_ok)")
    println("quality channel inert in every cell: $(all_inert)")
    if !controls_ok
        println("WARNING: a sensitivity control did NOT fire — inertness verdicts above " *
                "are NOT trustworthy for those cells.")
    end
    return nothing
end

# Run on direct invocation. `qualmer_corrector_quality_sensitivity.jl` includes this file
# to reuse the simulator, quality conditions and digest, and sets PROBE_INCLUDE_ONLY first
# so that including it defines the helpers WITHOUT launching a second full grid.
if !(@isdefined(PROBE_INCLUDE_ONLY) && PROBE_INCLUDE_ONLY)
    main()
end
