# Demoted-read reroute, BASE vs FIX: pre-registered analysis

**Verdict: VOID (G0 failed)**

The outcome statistic below is reported descriptively only; it is not a verdict.

## Inputs

- BASE label `td2rxh@ca9eb7366361+manifest85b75e3eb5c8` (16 shards, 192 rows)
- FIX label `td2rxh@dae40feaf697+manifest85b75e3eb5c8` (16 shards, 192 rows)
- Expected heads BASE `ca9eb7366361`, FIX `dae40feaf697`, shared manifest `85b75e3eb5c8`
- Label problems: none
- Unpaired (cell, seed): BASE 0, FIX 0

## G0 controls

- G0a, illumina identical in all 6 cells: PASS
- G0b, nanopore identical in the 2000 bp cells: FAIL (61 mismatching fields)
- G0c, oracle gate rows (except the historical-sha row) pass at both commits and illumina sha equal: PASS
- Columns excluded from identity: label, wall_seconds, wall_seconds_per_read

| run | rows_passed | gates_pass | failed_gates | illumina_sha | equals_historical | code_sha |
| --- | --- | --- | --- | --- | --- | --- |
| BASE | 19/19 | true |  | d36e3b6a10685346aa7b0238b48b4ab7fcefbed88f82cad7d959b0a831cdd311 | true | ca9eb736636101e6a4b3e6adf87687574608466b |
| FIX | 19/19 | true |  | d36e3b6a10685346aa7b0238b48b4ab7fcefbed88f82cad7d959b0a831cdd311 | true | dae40feaf6973c087ce956d08a25a7f7a46d5bfd |

All 61 mismatching fields (a = BASE, b = FIX):

| gate | cell | seed | arm | column | a | b |
| --- | --- | --- | --- | --- | --- | --- |
| G0b | 6000/2000/20x/e0.01 | 2 | nanopore | n_contigs | 439 | 437 |
| G0b | 6000/2000/20x/e0.01 | 2 | nanopore | total_contig_bases | 32578 | 32519 |
| G0b | 6000/2000/20x/e0.01 | 7 | nanopore | n_contigs | 461 | 462 |
| G0b | 6000/2000/20x/e0.01 | 7 | nanopore | total_contig_bases | 33503 | 33541 |
| G0b | 6000/2000/20x/e0.01 | 10 | nanopore | total_contig_bases | 34871 | 34873 |
| G0b | 6000/2000/20x/e0.01 | 11 | nanopore | n_contigs | 529 | 532 |
| G0b | 6000/2000/20x/e0.01 | 11 | nanopore | total_contig_bases | 38089 | 38179 |
| G0b | 6000/2000/20x/e0.01 | 14 | nanopore | total_contig_bases | 33493 | 33494 |
| G0b | 6000/2000/20x/e0.02 | 1 | nanopore | n_contigs | 1022 | 1025 |
| G0b | 6000/2000/20x/e0.02 | 1 | nanopore | total_contig_bases | 53369 | 53439 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | substitution_length_divergences | 93 | 90 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | window_anchor_rejections | 3 | 6 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | indel_engaged | 3 | 4 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | reassembly_k | 19 | 17 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | n_contigs | 1101 | 1077 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | largest_contig | 4021 | 3869 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | largest_ratio | 0.6702 | 0.6448 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | total_contig_bases | 46671 | 41500 |
| G0b | 6000/2000/20x/e0.02 | 2 | nanopore | best_contig_reference_coverage | 0.670167 | 0.644833 |
| G0b | 6000/2000/20x/e0.02 | 3 | nanopore | substitution_length_divergences | 98 | 96 |
| G0b | 6000/2000/20x/e0.02 | 3 | nanopore | window_anchor_rejections | 6 | 11 |
| G0b | 6000/2000/20x/e0.02 | 3 | nanopore | indel_engaged | 8 | 12 |
| G0b | 6000/2000/20x/e0.02 | 3 | nanopore | indel_requested | 1206 | 1207 |
| G0b | 6000/2000/20x/e0.02 | 3 | nanopore | n_contigs | 980 | 982 |
| G0b | 6000/2000/20x/e0.02 | 3 | nanopore | total_contig_bases | 67816 | 67743 |
| G0b | 6000/2000/20x/e0.02 | 8 | nanopore | substitution_length_divergences | 78 | 77 |
| G0b | 6000/2000/20x/e0.02 | 8 | nanopore | n_contigs | 1097 | 1108 |
| G0b | 6000/2000/20x/e0.02 | 8 | nanopore | total_contig_bases | 56071 | 56454 |
| G0b | 6000/2000/20x/e0.02 | 9 | nanopore | n_contigs | 1105 | 1108 |
| G0b | 6000/2000/20x/e0.02 | 9 | nanopore | total_contig_bases | 57418 | 57502 |
| G0b | 6000/2000/20x/e0.02 | 11 | nanopore | n_contigs | 1014 | 1012 |
| G0b | 6000/2000/20x/e0.02 | 11 | nanopore | total_contig_bases | 70359 | 70285 |
| G0b | 6000/2000/20x/e0.02 | 12 | nanopore | total_contig_bases | 45460 | 45459 |
| G0b | 6000/2000/20x/e0.02 | 13 | nanopore | n_contigs | 1124 | 1122 |
| G0b | 6000/2000/20x/e0.02 | 13 | nanopore | total_contig_bases | 43249 | 43215 |
| G0b | 6000/2000/20x/e0.05 | 1 | nanopore | substitution_length_divergences | 589 | 582 |
| G0b | 6000/2000/20x/e0.05 | 1 | nanopore | indel_engaged | 21 | 14 |
| G0b | 6000/2000/20x/e0.05 | 1 | nanopore | indel_requested | 1347 | 1344 |
| G0b | 6000/2000/20x/e0.05 | 1 | nanopore | n_contigs | 3166 | 3158 |
| G0b | 6000/2000/20x/e0.05 | 1 | nanopore | total_contig_bases | 74586 | 74250 |
| G0b | 6000/2000/20x/e0.05 | 3 | nanopore | indel_engaged | 15 | 13 |
| G0b | 6000/2000/20x/e0.05 | 3 | nanopore | indel_requested | 1322 | 1321 |
| G0b | 6000/2000/20x/e0.05 | 3 | nanopore | n_contigs | 3282 | 3280 |
| G0b | 6000/2000/20x/e0.05 | 3 | nanopore | total_contig_bases | 74913 | 74882 |
| G0b | 6000/2000/20x/e0.05 | 4 | nanopore | n_contigs | 3317 | 3318 |
| G0b | 6000/2000/20x/e0.05 | 4 | nanopore | total_contig_bases | 76234 | 76236 |
| G0b | 6000/2000/20x/e0.05 | 8 | nanopore | substitution_length_divergences | 575 | 579 |
| G0b | 6000/2000/20x/e0.05 | 8 | nanopore | indel_requested | 1331 | 1335 |
| G0b | 6000/2000/20x/e0.05 | 8 | nanopore | n_contigs | 3294 | 3305 |
| G0b | 6000/2000/20x/e0.05 | 8 | nanopore | total_contig_bases | 75474 | 75630 |
| G0b | 6000/2000/20x/e0.05 | 9 | nanopore | substitution_length_divergences | 496 | 498 |
| G0b | 6000/2000/20x/e0.05 | 9 | nanopore | n_contigs | 3253 | 3262 |
| G0b | 6000/2000/20x/e0.05 | 9 | nanopore | total_contig_bases | 74892 | 75081 |
| G0b | 6000/2000/20x/e0.05 | 12 | nanopore | substitution_length_divergences | 358 | 357 |
| G0b | 6000/2000/20x/e0.05 | 12 | nanopore | indel_engaged | 12 | 13 |
| G0b | 6000/2000/20x/e0.05 | 12 | nanopore | indel_requested | 978 | 979 |
| G0b | 6000/2000/20x/e0.05 | 12 | nanopore | n_contigs | 3413 | 3399 |
| G0b | 6000/2000/20x/e0.05 | 12 | nanopore | total_contig_bases | 76681 | 76452 |
| G0b | 6000/2000/20x/e0.05 | 15 | nanopore | indel_requested | 1322 | 1324 |
| G0b | 6000/2000/20x/e0.05 | 15 | nanopore | n_contigs | 3156 | 3139 |
| G0b | 6000/2000/20x/e0.05 | 15 | nanopore | total_contig_bases | 72926 | 72610 |

## G2 outcome (x = |d_BASE| − |d_FIX|, H1: x > 0)

| cell | gated | status | n_pairs | n_nonzero | n_zero | statistic | pvalue_one_sided | median_x | median_d_fix_minus_d_base | n_x_positive | median_d_base | median_d_fix | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6000/1000/20x/e0.01 | true | COMPLETE | 16 | 16 | 0 | -124.0 | 0.9998 | -0.1912 | -0.1912 | 2 | -0.4261 | -0.57 | exact one-sided sign-flip (65536 assignments) |
| 6000/1000/20x/e0.02 | false | COMPLETE | 16 | 16 | 0 | -64.0 | 0.9533 | -0.03267 | -0.09458 | 5 | -0.1085 | -0.2208 | exact one-sided sign-flip (65536 assignments) |
| 6000/1000/20x/e0.05 | false | COMPLETE | 16 | 15 | 1 | -40.0 | 0.8738 | -0.0095 | -0.007084 | 5 | 0.02033 | 0.02425 | exact one-sided sign-flip (32768 assignments) |

alpha 0.05.

### Per seed, 6000/1000/20x/e0.01

| cell | seed | d_base | d_fix | x | d_fix_minus_d_base |
| --- | --- | --- | --- | --- | --- |
| 6000/1000/20x/e0.01 | 1 | -0.4218 | -0.6357 | -0.2138 | -0.2138 |
| 6000/1000/20x/e0.01 | 2 | -0.2783 | -0.559 | -0.2807 | -0.2807 |
| 6000/1000/20x/e0.01 | 3 | -0.4508 | -0.676 | -0.2252 | -0.2252 |
| 6000/1000/20x/e0.01 | 4 | -0.3847 | -0.3078 | 0.07683 | 0.07683 |
| 6000/1000/20x/e0.01 | 5 | -0.2742 | -0.6758 | -0.4017 | -0.4017 |
| 6000/1000/20x/e0.01 | 6 | -0.657 | -0.7567 | -0.09967 | -0.09967 |
| 6000/1000/20x/e0.01 | 7 | -0.5698 | -0.5357 | 0.03417 | 0.03417 |
| 6000/1000/20x/e0.01 | 8 | -0.4942 | -0.7722 | -0.278 | -0.278 |
| 6000/1000/20x/e0.01 | 9 | -0.4983 | -0.5207 | -0.02233 | -0.02233 |
| 6000/1000/20x/e0.01 | 10 | -0.051 | -0.5047 | -0.4537 | -0.4537 |
| 6000/1000/20x/e0.01 | 11 | -0.4337 | -0.6023 | -0.1687 | -0.1687 |
| 6000/1000/20x/e0.01 | 12 | -0.4303 | -0.5193 | -0.089 | -0.089 |
| 6000/1000/20x/e0.01 | 13 | -0.389 | -0.7068 | -0.3178 | -0.3178 |
| 6000/1000/20x/e0.01 | 14 | -0.3865 | -0.462 | -0.0755 | -0.0755 |
| 6000/1000/20x/e0.01 | 15 | -0.2338 | -0.5138 | -0.28 | -0.28 |
| 6000/1000/20x/e0.01 | 16 | -0.4317 | -0.581 | -0.1493 | -0.1493 |

## Descriptive counters (not gates)

demoted_whole_read_routes: NOT CAPTURED. The td2rxh harness does not export it (amendment 5).

wall_seconds_per_read is confounded with node placement (BASE and FIX tasks ran on different nodes) and is not attributable to the code change.

| cell | arm | metric | median_base | median_fix | median_fix_minus_base | n_pairs |
| --- | --- | --- | --- | --- | --- | --- |
| 6000/1000/20x/e0.01 | illumina | substitution_length_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.01 | illumina | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.01 | illumina | window_anchor_rejections | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.01 | illumina | wall_seconds_per_read | 0.2532 | 0.2196 | -0.00585 | 16 |
| 6000/1000/20x/e0.01 | nanopore | substitution_length_divergences | 83.5 | 107.5 | 23.5 | 16 |
| 6000/1000/20x/e0.01 | nanopore | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.01 | nanopore | window_anchor_rejections | 10.5 | 10.0 | -1.5 | 16 |
| 6000/1000/20x/e0.01 | nanopore | wall_seconds_per_read | 5.578 | 3.862 | -0.1017 | 16 |
| 6000/1000/20x/e0.02 | illumina | substitution_length_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.02 | illumina | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.02 | illumina | window_anchor_rejections | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.02 | illumina | wall_seconds_per_read | 0.1964 | 0.19 | 0.00205 | 16 |
| 6000/1000/20x/e0.02 | nanopore | substitution_length_divergences | 177.0 | 189.0 | 15.5 | 16 |
| 6000/1000/20x/e0.02 | nanopore | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.02 | nanopore | window_anchor_rejections | 1.0 | 1.0 | 0.0 | 16 |
| 6000/1000/20x/e0.02 | nanopore | wall_seconds_per_read | 5.206 | 3.974 | -0.1183 | 16 |
| 6000/1000/20x/e0.05 | illumina | substitution_length_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.05 | illumina | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.05 | illumina | window_anchor_rejections | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.05 | illumina | wall_seconds_per_read | 0.2006 | 0.1837 | -0.00635 | 16 |
| 6000/1000/20x/e0.05 | nanopore | substitution_length_divergences | 595.5 | 359.5 | -241.5 | 16 |
| 6000/1000/20x/e0.05 | nanopore | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/1000/20x/e0.05 | nanopore | window_anchor_rejections | 1.5 | 1.0 | 0.0 | 16 |
| 6000/1000/20x/e0.05 | nanopore | wall_seconds_per_read | 5.479 | 4.606 | -0.1623 | 16 |
| 6000/2000/20x/e0.01 | illumina | substitution_length_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.01 | illumina | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.01 | illumina | window_anchor_rejections | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.01 | illumina | wall_seconds_per_read | 0.4683 | 0.4378 | -0.0109 | 16 |
| 6000/2000/20x/e0.01 | nanopore | substitution_length_divergences | 56.5 | 56.5 | 0.0 | 16 |
| 6000/2000/20x/e0.01 | nanopore | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.01 | nanopore | window_anchor_rejections | 23.0 | 23.0 | 0.0 | 16 |
| 6000/2000/20x/e0.01 | nanopore | wall_seconds_per_read | 9.995 | 7.085 | -0.1182 | 16 |
| 6000/2000/20x/e0.02 | illumina | substitution_length_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.02 | illumina | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.02 | illumina | window_anchor_rejections | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.02 | illumina | wall_seconds_per_read | 0.3998 | 0.3855 | 0.00205 | 16 |
| 6000/2000/20x/e0.02 | nanopore | substitution_length_divergences | 107.0 | 107.0 | 0.0 | 16 |
| 6000/2000/20x/e0.02 | nanopore | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.02 | nanopore | window_anchor_rejections | 5.5 | 6.5 | 0.0 | 16 |
| 6000/2000/20x/e0.02 | nanopore | wall_seconds_per_read | 9.58 | 7.02 | -0.1628 | 16 |
| 6000/2000/20x/e0.05 | illumina | substitution_length_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.05 | illumina | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.05 | illumina | window_anchor_rejections | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.05 | illumina | wall_seconds_per_read | 0.4128 | 0.4024 | 0.00785 | 16 |
| 6000/2000/20x/e0.05 | nanopore | substitution_length_divergences | 564.5 | 564.5 | 0.0 | 16 |
| 6000/2000/20x/e0.05 | nanopore | window_divergences | 0.0 | 0.0 | 0.0 | 16 |
| 6000/2000/20x/e0.05 | nanopore | window_anchor_rejections | 4.0 | 4.0 | 0.0 | 16 |
| 6000/2000/20x/e0.05 | nanopore | wall_seconds_per_read | 11.33 | 8.965 | -0.5166 | 16 |

## Context: earlier run vs BASE (drift; enters no gate)

| cell | arm | n_seeds_identical | n_seeds_compared | median_coverage_base_minus_context |
| --- | --- | --- | --- | --- |
| 6000/1000/20x/e0.01 | illumina | 16 | 16 | 0.0 |
| 6000/1000/20x/e0.01 | nanopore | 10 | 16 | 0.0 |
| 6000/1000/20x/e0.02 | illumina | 16 | 16 | 0.0 |
| 6000/1000/20x/e0.02 | nanopore | 7 | 16 | 0.0 |
| 6000/1000/20x/e0.05 | illumina | 16 | 16 | 0.0 |
| 6000/1000/20x/e0.05 | nanopore | 11 | 16 | 0.0 |
| 6000/2000/20x/e0.01 | illumina | 16 | 16 | 0.0 |
| 6000/2000/20x/e0.01 | nanopore | 13 | 16 | 0.0 |
| 6000/2000/20x/e0.02 | illumina | 16 | 16 | 0.0 |
| 6000/2000/20x/e0.02 | nanopore | 13 | 16 | 0.0 |
| 6000/2000/20x/e0.05 | illumina | 16 | 16 | 0.0 |
| 6000/2000/20x/e0.05 | nanopore | 13 | 16 | 0.0 |

