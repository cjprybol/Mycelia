# H2 K-shortest-paths exact-control benchmark summary

Scope: exact-control tier ONLY (bubble graph + repeat-choice graph),
from PLAN-2026-04-28-h2-kshortest-alternative-assemblies.md. The
realistic diploid (SK1/S288C HiFi) and repeat-rich (phage T4) tiers
named in that plan were NOT run in this script.

## Decision-rule verdicts (K=10, evidence-weighted)

- **bubble**: PASS (alt_recall_at_k=1.0, mrr=0.75)
- **repeat_choice**: PASS (alt_recall_at_k=1.0, mrr=0.75)

## Random-walk baseline (1000 walks/seed, seeds [42, 123, 456])

- bubble seed=42: 1000/1000 walks reached sink; recovers_true_alt=true; best_walk_score=0.5754
- bubble seed=123: 1000/1000 walks reached sink; recovers_true_alt=true; best_walk_score=0.5754
- bubble seed=456: 1000/1000 walks reached sink; recovers_true_alt=true; best_walk_score=0.5754
- repeat_choice seed=42: 1000/1000 walks reached sink; recovers_true_alt=true; best_walk_score=0.6931
- repeat_choice seed=123: 1000/1000 walks reached sink; recovers_true_alt=true; best_walk_score=0.6931
- repeat_choice seed=456: 1000/1000 walks reached sink; recovers_true_alt=true; best_walk_score=0.6931

## Scope caveats

- Exact-control tier only; realistic diploid and repeat-rich tiers are unrun and out of scope for this script.
- `Rhizomorph.probabilistic_walk_next` records the raw (un-normalized) edge weight in each `WalkStep`'s probability/cumulative_probability fields instead of the normalized transition probability it correctly uses internally for sampling; its `GraphPath.total_probability` is therefore not a valid probability (observed ~7.29e6 for a 7-edge bubble-graph walk instead of ~0.5625). This benchmark does not modify that library function; instead it independently recomputes each baseline walk's path probability from the realized vertex sequence using the same weight/total-outgoing-weight rule `k_shortest_paths` uses, so the baseline comparison in the decision rule is valid.
- Deduplication uses exact joined-vertex-label string matching, not SHA-256 hashing plus 99.5%-identity near-duplicate clustering as in the full plan's general deduplication section. This is sufficient and exact here because both graphs are constructed so that true and decoy arms are vertex-disjoint by design.
- Truth-matching uses exact vertex-subset matching against known arm-interior vertices, not sequence alignment (appropriate only for these hand-constructed toy graphs, per the plan's exact-control truth-matching plan).
- `decision_rule_pass` is only meaningful (and only computed as non-trivially-false) on the K=10, evidence-weighted row per graph; see `is_primary_decision_row`.
