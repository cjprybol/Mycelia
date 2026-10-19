# D5 ONT adaptive-grid AXIS REGISTRY — not the per-cell driver.
#
# Source of truth for the registered cell count (td-yddpr.4 step 2):
#   PLAN-2026-08-07-ont-adaptive-grid.md
#   DEC-2026-08-07 D5, amended 2026-09-11 (td-h3yc)
#
# 1152 ONT grid cells = identity(3) × coverage(4) × corrector(2) × k-policy(8)
#                     × organism(2) × seeds(3)
# plus 96 validation cells = dropped-rungs(16) × organism(2) × seeds(3)
# on the pre-declared stratum (identity ~99%, coverage 100x, corrector :iterative).
# Total 1248. Illumina is a separate control and is NOT in this product.
#
# This file exists so a committed invocation can name the matrix in CODE, not in
# a comment. The Track A harness (track_a_baseline_benchmark.jl) does not yet
# expose identity / corrector / k-policy; do not sbatch this grid until a driver
# imports these constants and threads them through simulate→assemble→QUAST.
#
# Do NOT run assemblies from this file. It has no simulate/assemble path.

module D5OntAdaptiveGridAxes

export D5_IDENTITIES, D5_COVERAGES, D5_CORRECTORS, D5_K_POLICIES,
       D5_ORGANISMS, D5_SEEDS, D5_COMPARATOR_RUNGS, D5_DROPPED_RUNGS,
       D5_CANDIDATE_LADDER, n_ont_grid_cells, n_validation_cells, n_d5_cells,
       ont_grid_cell, validation_cell

const D5_IDENTITIES = (90, 95, 99)          # mean percent; Badread strings live in the driver
const D5_COVERAGES = (10, 30, 50, 100)
const D5_CORRECTORS = (:none, :iterative)
# 8-level k-policy: adaptive level-set + seven fixed comparator rungs (td-h3yc).
const D5_COMPARATOR_RUNGS = (13, 15, 17, 19, 21, 23, 31)
const D5_K_POLICIES = (:adaptive, D5_COMPARATOR_RUNGS...)
const D5_ORGANISMS = (:Lambda, :T4)
const D5_SEEDS = (42, 123, 456)             # same three as Track A until D4 revises

# Full selection ladder: primes in [13, 101] plus {15, 21}. 23 rungs.
const D5_CANDIDATE_LADDER = (
    13, 15, 17, 19, 21, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101
)
const D5_DROPPED_RUNGS = Tuple(r for r in D5_CANDIDATE_LADDER if r ∉ D5_COMPARATOR_RUNGS)

n_ont_grid_cells() = length(D5_IDENTITIES) * length(D5_COVERAGES) * length(D5_CORRECTORS) *
                     length(D5_K_POLICIES) * length(D5_ORGANISMS) * length(D5_SEEDS)

n_validation_cells() = length(D5_DROPPED_RUNGS) * length(D5_ORGANISMS) * length(D5_SEEDS)

n_d5_cells() = n_ont_grid_cells() + n_validation_cells()

"""One ONT grid cell. `k_policy` is `:adaptive` or a fixed-k integer rung."""
function ont_grid_cell(identity, coverage, corrector, k_policy, organism, seed)
    return (; identity, coverage, corrector, k_policy, organism, seed, stratum = :grid)
end

"""Validation cell: dropped rung on the pre-declared high-identity stratum."""
function validation_cell(dropped_rung, organism, seed)
    return (;
        identity = 99,
        coverage = 100,
        corrector = :iterative,
        k_policy = dropped_rung,
        organism,
        seed,
        stratum = :validation,
    )
end

end # module
